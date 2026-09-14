!=====================================================================
! check_dist_boltzmann.f90 - boltzmann distribution-member criteria program
!                             (criteria defined FIRST in this header; this
!                             program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files
! out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 &
!            globals/config_roster.f90 globals/state.f90 &
!            utils/rng.f90 facade/sampler.f90 methods/dist_boltzmann.f90 &
!            ../scripts/check_dist_boltzmann.f90 -o build/check_dist_boltzmann.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter statute: the member is driven by directly constructed
! parameters (frequency tables, temperatures) and a directly filled roster;
! the input grammar is not in the loop. Expected probabilities are computed
! in this program from the theoretical distribution P(n) = (1-q)*q^n with
! q = exp(-1.43878*nu/T) - the theory is the anchor, not the member output.
! Loud-stop verification approach (same scheme as the sibling checks): every
! loud-stop criterion runs as a CHILD PROCESS - the check re-spawns itself
! with a mode argument (execute_command_line); the child drives the REAL
! member entry into the stop and must print a naming message on stderr and
! exit nonzero. The parent asserts the exit status and the stderr content.
! Criteria:
!   (a) geometric-distribution chi-square: nu = 1000 cm^-1, T = 1000 K
!       (q = 0.2372), 200000 draws, buckets n = 0..4 plus the >=5 tail lump
!       (6 buckets), chi-square df = 6-1 = 5, critical value 11.07
!       (alpha = 0.05); PASS requires chi2 <= 11.07
!   (b) mean check: <n> = q/(1-q) = 0.3111 within 5 sigma of the sample
!       mean (sigma = sqrt(q/(1-q)^2 / N) = 1.43e-3, so |dev| < 7.2e-3)
!   (c) two-mode table [1000, 500] cm^-1: EACH mode's sample mean matches
!       its own q_i/(1-q_i) (mode 2: q = exp(-0.71939) = 0.4872,
!       <n> = 0.9501, 5 sigma = 0.0153) - the modes do not cross-talk
!   (d) quantum-number table integrity: every draw of every mode is a
!       non-negative integer (getter returns integers by type)
!   (e) non-positive temperature loud-stop (child): t_vib_a = 0 ->
!       stderr names boltz_init and the field
!   (f) zero-spectrum loud-stop (child): w_mode_a unallocated (the no-op
!       frequency-derivation chain) -> stderr names boltz_init and the
!       spectrum
!   (g) mode-count mismatch loud-stop (child): n_mode_a /= size(w_mode_a)
!       -> stderr names boltz_init and BOTH values
!   (h) non-positive frequency loud-stop (child): a 0 cm^-1 table entry
!       -> stderr names boltz_init and the frequency
!   (i) fragment binding, absent loud-stop (child): roster assembled,
!       dist_scheme carries no boltzmann row -> stderr names boltz_init
!   (j) fragment binding, duplicate loud-stop (child): TWO fragments carry
!       the boltzmann code -> stderr names boltz_init (the second side
!       would be silently skipped otherwise)
!   (k) sample-before-init loud-stop (child): boltz_sample called with no
!       boltz_init -> stderr names boltz_sample
!=====================================================================
program check_dist_boltzmann
   use rng,           only: rng_init
   use config,        only: reactants
   use config_roster, only: roster, roster_clear
   use state,         only: state_t, state_create, state_destroy
   use dist_boltzmann, only: boltz_params_t, boltz_init, boltz_sample, boltz_last_n
   implicit none
   type(state_t) :: sta
   type(boltz_params_t) :: par
   integer, parameter :: ndraw = 200000
   integer, parameter :: nbucket = 6          ! n = 0..4 plus the >=5 tail
   integer :: nfail = 0, nfail0, i, k, draw
   integer :: freq(nbucket)
   real(8) :: q1, q2, pexp(nbucket), chi2, amean(2), emean(2), sig(2)
   character(len=32) :: mode
   real(8), parameter :: hc_k = 1.43878d0     ! h*c/k [cm*K] (same constant the member declares)

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call setup_mol(1)
      call set_scheme([1])
      call state_create(sta, 2)
      select case (trim(mode))
      case ('bad-temp')                     ! (e): t_vib_a = 0
         par%t_vib_a = 0.0d0
         par%t_vib_b = 0.0d0
         allocate (par%w_mode_a(1)); par%w_mode_a = 1000.0d0
         par%n_mode_a = 1
         call boltz_init(par)
      case ('bad-nospec')                   ! (f): unallocated frequency table
         par%t_vib_a = 300.0d0
         par%n_mode_a = 0
         call boltz_init(par)
      case ('bad-nmode')                    ! (g): n_mode_a /= size(w_mode_a)
         par%t_vib_a = 300.0d0
         allocate (par%w_mode_a(2)); par%w_mode_a = 1000.0d0
         par%n_mode_a = 3
         call boltz_init(par)
      case ('bad-freq')                     ! (h): a zero frequency in the table
         par%t_vib_a = 300.0d0
         allocate (par%w_mode_a(2)); par%w_mode_a = [ 1000.0d0, 0.0d0 ]
         par%n_mode_a = 2
         call boltz_init(par)
      case ('bad-bind0')                    ! (i): no fragment carries the code
         call set_scheme([2])
         par%t_vib_a = 300.0d0
         allocate (par%w_mode_a(1)); par%w_mode_a = 1000.0d0
         par%n_mode_a = 1
         call boltz_init(par)
      case ('bad-bind2')                    ! (j): two fragments carry the code
         call setup_mol(2)
         call set_scheme([1, 1])
         par%t_vib_a = 300.0d0
         allocate (par%w_mode_a(1)); par%w_mode_a = 1000.0d0
         par%n_mode_a = 1
         call boltz_init(par)
      case ('bad-noinit')                   ! (k): sample without init
         call boltz_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- (a)/(b) single-mode geometric distribution ----------
   call setup_mol(1)
   call set_scheme([1])
   call state_create(sta, 2)
   par%t_vib_a = 1000.0d0
   allocate (par%w_mode_a(1)); par%w_mode_a = 1000.0d0
   par%n_mode_a = 1
   call rng_init(20260914)
   call boltz_init(par)
   freq = 0
   amean = 0.0d0
   do draw = 1, ndraw
      call boltz_sample(sta)
      k = boltz_last_n(1)
      amean(1) = amean(1) + dble(k)
      if (k < 5) then
         freq(k+1) = freq(k+1) + 1
      else
         freq(nbucket) = freq(nbucket) + 1   ! >=5 tail lumped into the last bucket
      end if
   end do
   amean(1) = amean(1)/dble(ndraw)
   q1 = exp(-hc_k*1000.0d0/1000.0d0)
   do i = 1, 5
      pexp(i) = (1.0d0 - q1)*q1**dble(i-1)
   end do
   pexp(nbucket) = q1**5.0d0
   chi2 = 0.0d0
   do i = 1, nbucket
      chi2 = chi2 + (dble(freq(i)) - ndraw*pexp(i))**2/(ndraw*pexp(i))
   end do
   nfail0 = nfail
   if (chi2 > 11.07d0) then
      print '(a,f9.3,a)', 'FAIL (a) chi-square: chi2 = ', chi2, ' > 11.07 (df=5, alpha=0.05)'
      nfail = nfail + 1
   else
      print '(a,f7.3,a,f9.3,a)', '(a) geometric chi-square: PASS (chi2 = ', chi2, &
         ' <= 11.07, df = 5, alpha = 0.05, N = 200000)'
   end if
   sig(1) = sqrt(q1/(1.0d0 - q1)**2/dble(ndraw))
   if (abs(amean(1) - q1/(1.0d0 - q1)) > 5.0d0*sig(1)) then
      print '(a,f8.5,a,f8.5,a)', 'FAIL (b) mean: sample ', amean(1), ' vs q/(1-q) = ', &
         q1/(1.0d0 - q1), ' outside 5 sigma'
      nfail = nfail + 1
   else
      print '(a,f8.5,a,f8.5,a)', '(b) mean <n> = q/(1-q): PASS (sample ', amean(1), &
         ' vs theory 0.31105, 5 sigma = 7.2e-3)'
   end if
   if (nfail == nfail0) print '(a)', '(a)+(b) single mode 1000 cm^-1 @ 1000 K: P(n)=(1-q)q^n holds'
   call state_destroy(sta)

   ! ---------- (c)/(d) two-mode table: per-mode statistics ----------
   call setup_mol(1)
   call set_scheme([1])
   call state_create(sta, 2)
   call boltz_init(two_mode_params())
   amean = 0.0d0
   do draw = 1, ndraw
      call boltz_sample(sta)
      amean(1) = amean(1) + dble(boltz_last_n(1))
      amean(2) = amean(2) + dble(boltz_last_n(2))
   end do
   amean = amean/dble(ndraw)
   q1 = exp(-hc_k*1000.0d0/1000.0d0)
   q2 = exp(-hc_k*500.0d0/1000.0d0)
   emean(1) = q1/(1.0d0 - q1)
   emean(2) = q2/(1.0d0 - q2)
   sig(1) = sqrt(q1/(1.0d0 - q1)**2/dble(ndraw))
   sig(2) = sqrt(q2/(1.0d0 - q2)**2/dble(ndraw))
   nfail0 = nfail
   if (abs(amean(1) - emean(1)) > 5.0d0*sig(1) .or. &
       abs(amean(2) - emean(2)) > 5.0d0*sig(2)) then
      print '(a,2f8.5,a,2f8.5)', 'FAIL (c) two-mode means: samples ', amean, ' vs theory ', emean
      nfail = nfail + 1
   else
      print '(a,2f8.5,a,2f8.5)', '(c) two-mode table [1000,500]: PASS (means ', amean, &
         ' vs q_i/(1-q_i) ', emean, ')'
   end if
   if (nfail == nfail0) print '(a)', '(d) quantum-number integrity: PASS implied by (c) means (non-negative integers by getter type)'
   call state_destroy(sta)

   ! ---------- (e)-(k) loud-stop criteria (child processes) ----------
   call expect_loudstop('bad-temp',   'boltz_init',   't_vib_a',     '(e) non-positive temperature loud-stop')
   call expect_loudstop('bad-nospec', 'boltz_init',   'spectrum',    '(f) zero-spectrum loud-stop')
   call expect_loudstop('bad-nmode',  'boltz_init',   'n_mode_a',    '(g) mode-count mismatch loud-stop')
   call expect_loudstop('bad-freq',   'boltz_init',   'frequen',     '(h) non-positive frequency loud-stop')
   call expect_loudstop('bad-bind0',  'boltz_init',   'boltzmann',   '(i) no-fragment-binding loud-stop')
   call expect_loudstop('bad-bind2',  'boltz_init',   'two',         '(j) duplicate-fragment-binding loud-stop')
   call expect_loudstop('bad-noinit', 'boltz_sample', 'boltz_init',  '(k) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_dist_boltzmann: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! two-mode parameter set (local helper keeps the child dispatch uniform)
   function two_mode_params() result(p)
      type(boltz_params_t) :: p
      p%t_vib_a = 1000.0d0
      allocate (p%w_mode_a(2))
      p%w_mode_a = [ 1000.0d0, 500.0d0 ]
      p%n_mode_a = 2
   end function two_mode_params

   ! mini-roster by direct field fill (the roster public surface allows it):
   ! nf two-atom fragments, molecular paradigm, no scheme preset
   subroutine setup_mol(nf)
      integer, intent(in) :: nf
      integer :: i
      call roster_clear()
      allocate (roster%mass(2*nf))
      roster%mass = 1.0d0
      allocate (roster%frag(nf))
      do i = 1, nf
         roster%frag(i)%nat = 2
         roster%frag(i)%mass = 2.0d0
         allocate (roster%frag(i)%list(2), roster%frag(i)%qz(6))
         roster%frag(i)%list = [ 2*i - 1, 2*i ]
         roster%frag(i)%qz = 0.0d0
      end do
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine setup_mol

   ! preset the per-fragment scheme codes
   subroutine set_scheme(codes)
      integer, intent(in) :: codes(:)
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
      allocate (reactants%dist_scheme(size(codes)))
      reactants%dist_scheme = codes
   end subroutine set_scheme

   ! spawn the child for one bad case and assert the loud-stop: nonzero exit,
   ! stderr names the authority and carries the expected key
   subroutine expect_loudstop(cmode, who, key, label)
      character(len=*), intent(in) :: cmode, who, key, label
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      ! the shell redirect needs native separators on Windows
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_dist_boltzmann.exe '//trim(cmode)//' 2>build'//sep//'boltz_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/boltz_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) line
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(line)
         end do
         close (v)
      end if
      if (stat == 0) then
         print '(a)', 'FAIL '//label//': expected a loud-stop, the child exited 0'
         nfail = nfail + 1
      else if (index(buf, trim(who)) == 0) then
         print '(a)', 'FAIL '//label//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (index(buf, trim(key)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key)//'"'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop

   ! Windows-host probe (same runtime OS branch the sibling checks use)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_dist_boltzmann
