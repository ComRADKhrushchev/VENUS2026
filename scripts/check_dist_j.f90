!=====================================================================
! check_dist_j.f90 - j distribution-member criteria program (criteria
!                    defined FIRST in this header; this program executes
!                    them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 &
!            globals/config_roster.f90 globals/state.f90 &
!            utils/rng.f90 utils/linalg.f90 facade/sampler.f90 methods/dist_j.f90 &
!            ../scripts/check_dist_j.f90 -o build/check_dist_j.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter statute: the member is driven by directly constructed
! parameters and a directly filled roster; the input grammar is not in the
! loop. Expected probabilities are computed in this program from the
! theoretical distribution P(J) prop. (2J+1)*exp(-E_J/kT) with
! E_J = hbar^2*J(J+1)/(2*I), I the moment of inertia the CHECK computes
! independently from the roster geometry - the theory is the anchor, not
! the member output.
! Loud-stop verification approach (same scheme as the sibling checks):
! every loud-stop criterion runs as a CHILD PROCESS - the check re-spawns
! itself with a mode argument (execute_command_line); the child drives the
! REAL member entry into the stop and must print a naming message on
! stderr and exit nonzero. The parent asserts the exit status and the
! stderr content.
! Criteria:
!   (a) Boltzmann-J chi-square: linear carrier I = 0.5 amu*Ang^2 (two
!       unit masses at +/-0.5 Ang on x), T = 300 K -> c = hbar^2/(2*I*kb*T)
!       = 0.1617, peak J* = sqrt(1/(2c)) - 1/2 = 1.26; 200000 draws,
!       buckets J = 0..6 plus the >=7 tail lump (8 buckets), chi-square
!       df = 8-1 = 7, critical value 14.07 (alpha = 0.05); PASS requires
!       chi2 <= 14.07
!   (b) mean J: theory mean = sum(J*P~(J)) over the truncated support,
!       sample mean within 5 sigma (sigma = sqrt(Var/N), Var computed from
!       the same theoretical table)
!   (c) moment-of-inertia correctness: the member's exponent constant
!       getter j_c_rot() equals hbar^2/(2*I*kb*T) with I computed here
!       from the roster geometry, relative tolerance 1e-12 (a wrong I in
!       the member would shift the whole distribution)
!   (d) quantum-number integrity: every draw is an integer >= 0 and
!       bounded by the support bound j_j_max() (uniform draw + rejection
!       cannot leave [0, j_max])
!   (e) non-positive temperature loud-stop (child): t_rot = 0 -> stderr
!       names j_init and the field
!   (f) single-atom carrier loud-stop (child): nat = 1 -> stderr names
!       j_init and the atom count (no rotational degree of freedom)
!   (g) nonlinear carrier loud-stop (child): a bent three-atom carrier ->
!       stderr names j_init and the linearity requirement
!   (h) degenerate geometry loud-stop (child): two coincident carrier
!       atoms (I = 0) -> stderr names j_init and the moment of inertia
!   (i) no-fragment-binding loud-stop (child): dist_scheme carries no j
!       row -> stderr names j_init
!   (j) duplicate-fragment-binding loud-stop (child): two fragments carry
!       the j code -> stderr names j_init
!   (k) sample-before-init loud-stop (child): j_sample with no j_init ->
!       stderr names j_sample
!   (l) support-scan cap loud-stop (child): t_rot = 1e15 K makes the
!       exponent constant c ~ 5e-17, the Boltzmann-J ladder effectively
!       flat over the whole 1e7 cap (peak J* ~ 1e8 lies beyond it) ->
!       stderr names j_init, the cap and t_rot
!=====================================================================
program check_dist_j
   use consts,        only: hbar_code, kb_code
   use rng,           only: rng_init
   use config,        only: reactants
   use config_roster, only: roster, roster_clear
   use state,         only: state_t, state_create, state_destroy
   use dist_j,        only: j_params_t, j_init, j_sample, j_last_j, j_c_rot, j_j_max
   implicit none
   type(state_t) :: sta
   type(j_params_t) :: par
   integer, parameter :: ndraw = 200000
   integer, parameter :: nbucket = 8          ! J = 0..6 plus the >=7 tail
   integer :: nfail = 0, nfail0, i, k, draw, jmax
   integer :: freq(nbucket)
   real(8) :: c_theory, pnorm, pexp(nbucket), chi2
   real(8) :: amean, emean, varr, sig, jmax_seen
   real(8) :: com(3), r(3), i_perp
   character(len=32) :: mode

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 2)
      par%t_rot = 300.0d0
      select case (trim(mode))
      case ('bad-temp')                     ! (e): t_rot = 0
         call setup_mol(1, linear_geometry())
         call set_scheme([4])
         par%t_rot = 0.0d0
         call j_init(par)
      case ('bad-scan')                     ! (l): ladder flat to the 1e7 cap
         call setup_mol(1, linear_geometry())
         call set_scheme([4])
         par%t_rot = 1.0d15
         call j_init(par)
      case ('bad-atom')                     ! (f): single-atom carrier
         call setup_mol(1, atom_geometry())
         call set_scheme([4])
         call j_init(par)
      case ('bad-bent')                     ! (g): bent three-atom carrier
         call setup_mol(1, bent_geometry())
         call set_scheme([4])
         call j_init(par)
      case ('bad-degen')                    ! (h): coincident atoms, I = 0
         call setup_mol(1, degenerate_geometry())
         call set_scheme([4])
         call j_init(par)
      case ('bad-bind0')                    ! (i): no fragment carries the code
         call setup_mol(1, linear_geometry())
         call set_scheme([1])
         call j_init(par)
      case ('bad-bind2')                    ! (j): two fragments carry the code
         call setup_mol(2, linear_geometry())
         call set_scheme([4, 4])
         call j_init(par)
      case ('bad-noinit')                   ! (k): sample without init
         call setup_mol(1, linear_geometry())
         call set_scheme([4])
         call j_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- (a)/(b)/(c)/(d) single linear carrier ----------
   call setup_mol(1, linear_geometry())
   call set_scheme([4])
   call state_create(sta, 2)
   par%t_rot = 300.0d0
   call rng_init(7319)
   call j_init(par)

   ! theory: I from the roster geometry computed here (COM-frame mass sum,
   ! the linear-rotor perpendicular moment = half the inertia-tensor trace)
   com = 0.0d0
   do i = 1, 2
      com = com + roster%mass(i)*roster%frag(1)%qz(3*i-2:3*i)
   end do
   com = com/sum(roster%mass(1:2))
   i_perp = 0.0d0
   do i = 1, 2
      r = roster%frag(1)%qz(3*i-2:3*i) - com
      i_perp = i_perp + roster%mass(i)*dot_product(r, r)
   end do
   ! (trace of the inertia tensor = 2*sum m|r|^2, and the linear rotor's
   ! two equal perpendicular moments make trace = 2*I_perp - so I_perp is
   ! the mass sum itself, no further factor)
   c_theory = hbar_code**2/(2.0d0*i_perp*kb_code*300.0d0)

   nfail0 = nfail
   if (abs(j_c_rot() - c_theory) > 1.0d-12*c_theory) then
      print '(a,es12.4,a,es12.4)', 'FAIL (c) exponent constant: member ', j_c_rot(), &
         ' vs theory ', c_theory
      nfail = nfail + 1
   else
      print '(a,es12.4,a,f6.3,a)', '(c) exponent constant: PASS (member c = ', j_c_rot(), &
         ' = hbar^2/(2*I*kb*T), I = ', i_perp, ' amu*Ang^2)'
   end if

   freq = 0
   amean = 0.0d0
   jmax_seen = 0.0d0
   do draw = 1, ndraw
      call j_sample(sta)
      k = j_last_j()
      amean = amean + dble(k)
      jmax_seen = max(jmax_seen, dble(k))
      if (k < 7) then
         freq(k+1) = freq(k+1) + 1
      else
         freq(nbucket) = freq(nbucket) + 1     ! >=7 tail lumped into the last bucket
      end if
   end do
   amean = amean/dble(ndraw)

   pnorm = 0.0d0
   do i = 0, 1000                              ! truncated-support normalization
      pnorm = pnorm + p_raw(i)
   end do
   do i = 1, 7
      pexp(i) = p_raw(i-1)/pnorm
   end do
   pexp(nbucket) = 0.0d0
   do i = 7, 1000
      pexp(nbucket) = pexp(nbucket) + p_raw(i)/pnorm
   end do
   chi2 = 0.0d0
   do i = 1, nbucket
      chi2 = chi2 + (dble(freq(i)) - ndraw*pexp(i))**2/(ndraw*pexp(i))
   end do
   if (chi2 > 14.07d0) then
      print '(a,f9.3,a)', 'FAIL (a) chi-square: chi2 = ', chi2, ' > 14.07 (df=7, alpha=0.05)'
      nfail = nfail + 1
   else
      print '(a,f7.3,a)', '(a) Boltzmann-J chi-square: PASS (chi2 = ', chi2, &
         ' <= 14.07, df = 7, alpha = 0.05, N = 200000)'
   end if

   emean = 0.0d0
   varr = 0.0d0
   do i = 0, 1000
      emean = emean + dble(i)*p_raw(i)/pnorm
   end do
   do i = 0, 1000
      varr = varr + (dble(i) - emean)**2*p_raw(i)/pnorm
   end do
   sig = sqrt(varr/dble(ndraw))
   if (abs(amean - emean) > 5.0d0*sig) then
      print '(a,f8.4,a,f8.4,a)', 'FAIL (b) mean J: sample ', amean, ' vs theory ', emean, &
         ' outside 5 sigma'
      nfail = nfail + 1
   else
      print '(a,f8.4,a,f8.4,a,f7.4,a)', '(b) mean J: PASS (sample ', amean, ' vs theory ', &
         emean, ', 5 sigma = ', 5.0d0*sig, ')'
   end if

   jmax = j_j_max()
   if (jmax < 1 .or. jmax_seen > dble(jmax) .or. jmax > 1000000) then
      print '(a,i0)', 'FAIL (d) quantum-number integrity: draws outside [0, j_max] or runaway j_max: ', jmax
      nfail = nfail + 1
   else
      print '(a,i0,a,f6.1,a)', '(d) quantum-number integrity: PASS (all draws integers in [0, ', jmax, &
         '], largest observed ', jmax_seen, ')'
   end if
   call state_destroy(sta)

   ! ---------- (e)-(k) loud-stop criteria (child processes) ----------
   call expect_loudstop('bad-temp',   'j_init',   't_rot',      '(e) non-positive temperature loud-stop')
   call expect_loudstop('bad-atom',   'j_init',   'atom count', '(f) single-atom carrier loud-stop')
   call expect_loudstop('bad-bent',   'j_init',   'linear',     '(g) nonlinear carrier loud-stop')
   call expect_loudstop('bad-degen',  'j_init',   'inertia',    '(h) degenerate-geometry loud-stop')
   call expect_loudstop('bad-bind0',  'j_init',   'j scheme',   '(i) no-fragment-binding loud-stop')
   call expect_loudstop('bad-bind2',  'j_init',   'two',        '(j) duplicate-fragment-binding loud-stop')
   call expect_loudstop('bad-noinit', 'j_sample', 'j_init',     '(k) sample-before-init loud-stop')
   call expect_loudstop('bad-scan',  'j_init',   'cap',        '(l) support-scan cap loud-stop')

   print '(a,i0,a)', '=== check_dist_j: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! unnormalized Boltzmann-J weight (2J+1)*exp(-c*J*(J+1))
   real(8) function p_raw(j)
      integer, intent(in) :: j
      p_raw = dble(2*j + 1)*exp(-c_theory*dble(j*(j + 1)))
   end function p_raw

   ! linear two-atom geometry: unit masses at +/-0.5 Ang on x (I_perp = 0.5)
   function linear_geometry() result(g)
      real(8) :: g(6)
      g = 0.0d0
      g(1) = -0.5d0
      g(4) = 0.5d0
   end function linear_geometry

   ! single-atom "fragment"
   function atom_geometry() result(g)
      real(8) :: g(3)
      g = 0.0d0
   end function atom_geometry

   ! bent three-atom geometry (not collinear)
   function bent_geometry() result(g)
      real(8) :: g(9)
      g = 0.0d0
      g(1) = -0.6d0
      g(4) = 0.6d0
      g(8) = 0.6d0
   end function bent_geometry

   ! two coincident atoms (I = 0)
   function degenerate_geometry() result(g)
      real(8) :: g(6)
      g = 0.0d0
      g(4) = 0.0d0
   end function degenerate_geometry

   ! mini-roster by direct field fill: nf fragments whose geometry comes
   ! from the argument (3*nat values packed); unit masses, molecular paradigm
   subroutine setup_mol(nf, g)
      integer, intent(in) :: nf
      real(8), intent(in) :: g(:)
      integer :: nat, i, k
      call roster_clear()
      nat = size(g)/3
      allocate (roster%mass(nat*nf))
      roster%mass = 1.0d0
      allocate (roster%frag(nf))
      do i = 1, nf
         roster%frag(i)%nat = nat
         roster%frag(i)%mass = dble(nat)
         allocate (roster%frag(i)%list(nat), roster%frag(i)%qz(size(g)))
         roster%frag(i)%list = [ ((i - 1)*nat + k, k = 1, nat) ]
         roster%frag(i)%qz = g
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
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_dist_j.exe '//trim(cmode)//' 2>build'//sep//'dist_j_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/dist_j_err.txt', status='old', action='read', iostat=ios)
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

end program check_dist_j
