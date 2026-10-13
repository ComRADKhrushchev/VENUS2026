!=====================================================================
! check_samp_boltzmann.f90 - boltzmann distribution-member criteria program
!                             (criteria defined FIRST in this header; this
!                             program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files
! out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 &
!            globals/config_atoms.f90 globals/state.f90 utils/rng.f90 utils/geometry.f90 &
!            utils/linalg.f90 utils/hessian.f90 utils/geom_init.f90 utils/spectrum.f90 &
!            interface/elec_interface.f90 interface/force_interface.f90 dynamics/propagator.f90 &
!            interface/bath.f90 interface/sampler.f90 methods/samp_boltzmann.f90 &
!            ../scripts/check_samp_boltzmann.f90 -o build/check_samp_boltzmann.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter rule: the member is driven by directly constructed
! parameters (frequency tables, temperatures) and a directly filled list_atoms;
! the input grammar is not in the loop. Expected probabilities are computed
! in this program from the theoretical distribution P(n) = (1-q)*q^n with
! q = exp(-1.43878*nu/T) - the theory is the anchor, not the member output.
! Named-abort verification approach (same scheme as the sibling checks): every
! named-abort criterion runs as a CHILD PROCESS - the check re-spawns itself
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
!   (e) non-positive temperature named abort (child): t_vib = 0 ->
!       stderr names boltz_init and the field
!   (f) zero-spectrum named abort (child): w_mode_c unallocated (the no-op
!       frequency-derivation chain) -> stderr names boltz_init and the
!       spectrum
!   (g) inner-table unallocated named abort (child): the table set carries
!       an unallocated component table -> stderr names boltz_init and the
!       unallocated carrier table
!   (h) non-positive frequency named abort (child): a 0 cm^-1 table entry
!       -> stderr names boltz_init and the frequency
!   (i) fragment binding, absent named abort (child): list_atoms assembled,
!       dist_scheme carries no boltzmann row -> stderr names boltz_init
!   (j) table-count mismatch named abort (child): TWO fragments carry the
!       boltzmann code but the parameter set carries ONE table -> stderr
!       names boltz_init and the carrier count (a same-word double
!       selection is legal only with one table per carrier)
!   (k) sample-before-init named abort (child): boltz_sample called with no
!       boltz_init -> stderr names boltz_draw (the two-phase guard site)
!   (l) dual-carrier equal treatment (positive): two fragments, two ragged
!       tables ([1000] and [1000, 500] cm^-1) at ONE temperature -> each
!       carrier's per-mode sample mean matches its own q_i/(1-q_i) (the
!       carriers do not cross-talk; one parameter set serves both)
!=====================================================================
program check_samp_boltzmann
   use rng,           only: rng_init
   use config,        only: reactants
   use config_atoms, only: list_atoms, list_atoms_clear
   use state,         only: state_t, state_create, state_destroy
   use samp_boltzmann, only: boltz_params_t, boltz_init, boltz_draw, boltz_realize, &
                           boltz_sample, boltz_last_n
   use spectrum,      only: w_to_wvn
   use consts,        only: hbar_act
   implicit none
   type(state_t) :: sta
   type(boltz_params_t) :: par
   integer, parameter :: ndraw = 200000
   integer, parameter :: nbucket = 6          ! n = 0..4 plus the >=5 tail
   integer :: nfail = 0, nfail0, i, k, draw
   integer :: freq(nbucket)
   real(8) :: q1, q2, pexp(nbucket), chi2, amean(2), emean(2), sig(2)
   real(8) :: qm(1), pm(1), w_use, e_shell, e_theory, p_tot(3)
   character(len=32) :: mode
   real(8), parameter :: hc_k = 1.43878d0     ! h*c/k [cm*K] (same constant the member declares)

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   allocate (par%t_vib(1))   ! single-carrier check face
   par%t_vib(1) = 0.0d0
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call setup_mol(1)
      call set_scheme([1])
      call state_create(sta, 2)
      select case (trim(mode))
      case ('bad-temp')                     ! (e): t_vib = 0
         par%t_vib(1) = 0.0d0
         call fill_tbl(par, [ 1000.0d0 ])
         call boltz_init(par)
      case ('bad-nospec')                   ! (f): unallocated table set
         par%t_vib(1) = 300.0d0
         call boltz_init(par)
      case ('bad-tbl')                      ! (g): an unallocated inner table
         par%t_vib(1) = 300.0d0
         allocate (par%tbl(2))
         call fill_tbl(par, [ 1000.0d0 ])   ! (fills row 1; row 2 stays unallocated)
         call boltz_init(par)
      case ('bad-freq')                     ! (h): a zero frequency in the table
         par%t_vib(1) = 300.0d0
         call fill_tbl(par, [ 1000.0d0, 0.0d0 ])
         call boltz_init(par)
      case ('bad-bind0')                    ! (i): no fragment carries the code
         call set_scheme([2])
         par%t_vib(1) = 300.0d0
         call fill_tbl(par, [ 1000.0d0 ])
         call boltz_init(par)
      case ('bad-bind2')                    ! (j): two carriers, one table
         call setup_mol(2)
         call set_scheme([1, 1])
         par%t_vib = [ 300.0d0, 300.0d0 ]
         call fill_tbl(par, [ 1000.0d0 ])
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
   par%t_vib(1) = 1000.0d0
   call fill_tbl(par, [ 1000.0d0 ])
   call rng_init(20260914)
   call boltz_init(par)
   freq = 0
   amean = 0.0d0
   do draw = 1, ndraw
      call boltz_sample(sta)
      k = boltz_last_n(1, 1)
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

   ! ---------- (i) the quasi-classical realization: on-shell per-mode
   ! energies (n + 1/2)*hbar*w, zero net momentum, displaced geometry ----------
   nfail0 = nfail
   call rng_init(20260915)
   call boltz_draw()
   call boltz_realize(sta)
   qm(1) = dot_product(sta%q(1:3) - sta%q(4:6), [ 1.0d0, 0.0d0, 0.0d0 ])/sqrt(2.0d0)
   pm(1) = dot_product(sta%p(1:3) - sta%p(4:6), [ 1.0d0, 0.0d0, 0.0d0 ])/sqrt(2.0d0)
   w_use = 1000.0d0/w_to_wvn(1.0d0)
   e_shell = 0.5d0*(pm(1)**2 + w_use**2*qm(1)**2)
   e_theory = (dble(boltz_last_n(1, 1)) + 0.5d0)*hbar_act*w_use
   p_tot(1:3) = sta%p(1:3) + sta%p(4:6)
   if (abs(e_shell - e_theory) > 1.0d-9*e_theory .or. norm2(p_tot(1:3)) > 1.0d-9 .or. &
       all(sta%q == 0.0d0)) then
      print '(a)', 'FAIL (i) realization: mode off shell, net momentum nonzero or q undisplaced'
      nfail = nfail + 1
   else
      print '(a,es10.3,a)', '(i) realization: PASS (mode exactly on its (n+1/2)*hbar*w '// &
         'shell, E = ', e_theory, '; zero net momentum; geometry displaced)'
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
      amean(1) = amean(1) + dble(boltz_last_n(1, 1))
      amean(2) = amean(2) + dble(boltz_last_n(1, 2))
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
   if (nfail == nfail0) print '(a)', '(d) quantum-number integrity: PASS implied by (c) means (on-negative ints by getter type)'
   call state_destroy(sta)

   ! ---------- (l) dual-carrier equal treatment: two fragments, two ragged
   ! tables, one temperature - each carrier's modes follow their own q ----
   call setup_mol(2)
   call set_scheme([1, 1])
   call state_create(sta, 4)
   call boltz_init(dual_params())
   amean = 0.0d0
   do draw = 1, ndraw
      call boltz_sample(sta)
      amean(1) = amean(1) + dble(boltz_last_n(1, 1))   ! carrier 1, mode 1000
      amean(2) = amean(2) + dble(boltz_last_n(2, 2))   ! carrier 2, mode 500
   end do
   amean = amean/dble(ndraw)
   q1 = exp(-hc_k*1000.0d0/1000.0d0)
   q2 = exp(-hc_k*500.0d0/1000.0d0)
   nfail0 = nfail
   if (abs(amean(1) - q1/(1.0d0 - q1)) > 5.0d0*sqrt(q1/(1.0d0 - q1)**2/dble(ndraw)) .or. &
       abs(amean(2) - q2/(1.0d0 - q2)) > 5.0d0*sqrt(q2/(1.0d0 - q2)**2/dble(ndraw))) then
      print '(a,2f8.5)', 'FAIL (l) dual-carrier means: ', amean
      nfail = nfail + 1
   else
      print '(a,2f8.5,a)', '(l) dual-carrier equal treatment: PASS (means ', amean, &
         ' - each carrier follows its own table at the one temperature)'
   end if
   call state_destroy(sta)

   ! ---------- (e)-(k) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-temp',   'boltz_init',   't_vib',      '(e) non-positive temperature loud-stop')
   call expect_loudstop('bad-nospec', 'boltz_init',   'spectrum',    '(f) zero-spectrum loud-stop')
   call expect_loudstop('bad-tbl',    'boltz_init',   'unallocated', '(g) inner-table unallocated loud-stop')
   call expect_loudstop('bad-freq',   'boltz_init',   'frequen',     '(h) non-positive frequency loud-stop')
   call expect_loudstop('bad-bind0',  'boltz_init',   'boltzmann',   '(i) no-fragment-binding loud-stop')
   call expect_loudstop('bad-bind2',  'boltz_init',   'carriers',    '(j) table-count/carrier-count mismatch loud-stop')
   call expect_loudstop('bad-noinit', 'boltz_draw', 'boltz_init',  '(k) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_samp_boltzmann: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! full-table builder row 1: k-th stretch column (x/y/z per mode index)
   ! for the equal-mass 2-atom carrier; w chosen so w_to_wvn(w) = nu
   subroutine fill_tbl(par, nus)
      type(boltz_params_t), intent(inout) :: par
      real(8), intent(in) :: nus(:)
      integer :: k, m
      if (.not. allocated(par%tbl)) allocate (par%tbl(1))
      m = size(nus)
      allocate (par%tbl(1)%w(m), par%tbl(1)%c(6, m))
      par%tbl(1)%c = 0.0d0
      do k = 1, m
         par%tbl(1)%w(k) = nus(k)/w_to_wvn(1.0d0)
         par%tbl(1)%c(k, k) = 1.0d0/sqrt(2.0d0)
         par%tbl(1)%c(k+3, k) = -1.0d0/sqrt(2.0d0)
      end do
   end subroutine fill_tbl

   ! two-mode parameter set (local helper keeps the child dispatch uniform)
   function two_mode_params() result(p)
      type(boltz_params_t) :: p
      allocate (p%t_vib(1))   ! single-carrier check face
      p%t_vib(1) = 0.0d0
      p%t_vib(1) = 1000.0d0
      allocate (p%tbl(1))
      call fill_tbl(p, [ 1000.0d0, 500.0d0 ])
   end function two_mode_params

   ! (l) two-carrier ragged tables: carrier 1 [1000], carrier 2 [1000, 500]
   function dual_params() result(p)
      type(boltz_params_t) :: p
      allocate (p%t_vib(2))   ! dual-carrier check face
      p%t_vib = 1000.0d0
      allocate (p%tbl(2))
      call fill_tbl(p, [ 1000.0d0 ])                    ! carrier 1: one mode
      call fill_tbl2(p, [ 1000.0d0, 500.0d0 ])          ! carrier 2: two modes
   end function dual_params

   ! row-2 variant of the table builder (the ragged dual-carrier helper)
   subroutine fill_tbl2(par, nus)
      type(boltz_params_t), intent(inout) :: par
      real(8), intent(in) :: nus(:)
      integer :: k, m
      m = size(nus)
      allocate (par%tbl(2)%w(m), par%tbl(2)%c(6, m))
      par%tbl(2)%c = 0.0d0
      do k = 1, m
         par%tbl(2)%w(k) = nus(k)/w_to_wvn(1.0d0)
         par%tbl(2)%c(k, k) = 1.0d0/sqrt(2.0d0)
         par%tbl(2)%c(k+3, k) = -1.0d0/sqrt(2.0d0)
      end do
   end subroutine fill_tbl2

   ! mini-list_atoms by direct field fill (the atom list public surface allows it):
   ! nf two-atom fragments, molecular paradigm, no scheme preset
   subroutine setup_mol(nf)
      integer, intent(in) :: nf
      integer :: i
      call list_atoms_clear()
      allocate (list_atoms%mass(2*nf))
      list_atoms%mass = 1.0d0
      allocate (list_atoms%frag(nf))
      do i = 1, nf
         list_atoms%frag(i)%nat = 2
         list_atoms%frag(i)%mass = 2.0d0
         allocate (list_atoms%frag(i)%list(2), list_atoms%frag(i)%qz(6))
         list_atoms%frag(i)%list = [ 2*i - 1, 2*i ]
         list_atoms%frag(i)%qz = 0.0d0
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

   ! spawn the child for one bad case and assert the named abort: nonzero exit,
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
      cmd = 'build'//sep//'check_samp_boltzmann.exe '//trim(cmode)//' 2>build'//sep//'boltz_err.txt'
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

end program check_samp_boltzmann
