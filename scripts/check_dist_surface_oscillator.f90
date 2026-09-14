!=====================================================================
! check_dist_surface_oscillator.f90 - surface_oscillator distribution-
!                                     member criteria program (criteria
!                                     defined FIRST in this header; this
!                                     program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 &
!            globals/config_roster.f90 globals/state.f90 &
!            utils/rng.f90 facade/sampler.f90 methods/dist_surface_oscillator.f90 &
!            ../scripts/check_dist_surface_oscillator.f90 -o build/check_dist_surface_oscillator.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter statute: the member is driven by directly
! constructed parameters and a directly filled roster (a two-atom
! atom-anchor carrier: masses 1 amu and 5000 amu, bond 2.3 Ang along a
! skewed direction - a heavy anchor standing in for a surface atom); the
! input grammar is not in the loop. Expected distributions are computed
! in this program from theory, not from the member.
! Loud-stop verification approach (same scheme as the sibling checks):
! every loud-stop criterion runs as a CHILD PROCESS - the check re-spawns
! itself with a mode argument (execute_command_line); the child drives
! the REAL member entry into the stop and must print a naming message on
! stderr and exit nonzero. The parent asserts the exit status and the
! stderr content.
! Criteria (carrier: mu = 1*5000/5001 = 0.99980 amu; k_osc = 42.0
! kcal/mol/Ang^2 -> omega = sqrt(k*e_conv/mu) = 1.3258 rad/(10 fs);
! anchors in internal energy units):
!   (a) quantum-mode energy exact: n_level = 3 -> E = 3.5*hbar*omega;
!       per draw the back-projected energy
!       E_bp = 0.5*k_int*delta^2 + mu*v^2/2 (delta = bond-length change,
!       v = relative velocity along the bond) matches 3.5*hbar*omega to
!       1e-10 relative, 5000 draws; the getter surface_oscillator_last_e
!       matches the same value
!   (b) thermal-mode exponential energy: t_osc = 180 K, E = -kb*T*ln(u);
!       equal-CDF-bin chi-square, 10 bins over the energy quantiles
!       F(E) = 1-exp(-E/(kb*T)), N = 50000, df = 9, critical 16.92
!       (alpha = 0.05); mean = kb*T within 5 sigma (sigma = kb*T/sqrt(N))
!   (c) phase-uniformity identities (quantum mode): mean displacement =
!       0 within 5 sigma (sigma = A/sqrt(2N), N = 20000); mean kinetic
!       energy = E/2 within 5 sigma (sigma = E/sqrt(8N) - the sin^2(phi)
!       variance)
!   (d) momentum structure: p_atom = -p_anchor exactly (1e-14 relative);
!       both momenta along the bond direction (|p x u_bond| = 0) and the
!       COM stays fixed (atom-anchor COM equals the equilibrium COM to
!       1e-12 Ang)
!   (e) bond placement: the written bond vector is the equilibrium
!       direction exactly; the bond length is r_eq + delta with
!       delta = A*cos(phi) in [-A, A]; the displacement is purely along
!       the bond (perpendicular deviation 0 to 1e-12 Ang)
!   (f) bad energy-mode selector loud-stop (child): n_osc = 2 -> stderr
!       names surface_oscillator_init and n_osc
!   (g) non-positive temperature loud-stop (child): mode 0 with t_osc =
!       0 -> stderr names surface_oscillator_init and t_osc
!   (h) negative quantum number loud-stop (child): mode 1 with n_level =
!       -1 -> stderr names surface_oscillator_init and n_level
!   (i) non-positive force constant loud-stop (child): k_osc = 0 ->
!       stderr names surface_oscillator_init and k_osc
!   (j) non-diatomic carrier loud-stop (child): a three-atom carrier ->
!       stderr names surface_oscillator_init and the atom count
!   (k) degenerate bond loud-stop (child): coincident atoms -> stderr
!       names surface_oscillator_init and the bond
!   (l) no-fragment-binding loud-stop (child) and (m) duplicate-carrier
!       loud-stop (child): stderr names surface_oscillator_init
!   (n) sample-before-init loud-stop (child): stderr names
!       surface_oscillator_sample
!=====================================================================
program check_dist_surface_oscillator
   use consts,        only: hbar_code, kb_code, e_conv, two_pi
   use rng,           only: rng_init
   use config,        only: reactants
   use config_roster, only: roster, roster_clear
   use state,         only: state_t, state_create, state_destroy
   use dist_surface_oscillator, only: surf_osc_params_t, surface_oscillator_init, &
                                      surface_oscillator_sample, surface_oscillator_last_e
   implicit none
   type(state_t) :: sta
   type(surf_osc_params_t) :: par
   integer, parameter :: ndraw_exact = 5000
   integer, parameter :: ndraw_stat = 50000
   integer, parameter :: ndraw_phase = 20000
   integer, parameter :: nbin = 10
   integer :: nfail = 0, nfail0, i, draw, ib
   integer :: freq(nbin)
   real(8) :: m_atom, m_anch, mu, k_int, omega, e_quant
   real(8) :: r_eq, u_bond(3), q_eq(6), com_eq(3)
   real(8) :: bond(3), blen, delta, vrel, vproj
   real(8) :: p_atom(3), p_anch(3), com_new(3)
   real(8) :: e_bp, amean, emean, sig, chi2, ebin, cdf
   real(8) :: dmean, kmean, ksig, dvar
   real(8) :: a_quant
   character(len=32) :: mode
   logical :: ok_all

   m_atom = 1.0d0
   m_anch = 5000.0d0
   mu = m_atom*m_anch/(m_atom + m_anch)
   k_int = 42.0d0*e_conv                       ! internal force constant
   omega = sqrt(k_int/mu)
   e_quant = 3.5d0*hbar_code*omega              ! n_level = 3 zero-point-inclusive level
   a_quant = sqrt(2.0d0*e_quant/k_int)          ! classical amplitude at that energy

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 3)
      par%n_osc = 1
      par%n_level = 3
      par%k_osc = 42.0d0
      select case (trim(mode))
      case ('bad-nosc')                    ! (f): n_osc = 2
         call setup_pair()
         call set_scheme([5])
         par%n_osc = 2
         call surface_oscillator_init(par)
      case ('bad-tosc')                    ! (g): mode 0, t_osc = 0
         call setup_pair()
         call set_scheme([5])
         par%n_osc = 0
         par%t_osc = 0.0d0
         call surface_oscillator_init(par)
      case ('bad-nlev')                    ! (h): mode 1, n_level = -1
         call setup_pair()
         call set_scheme([5])
         par%n_level = -1
         call surface_oscillator_init(par)
      case ('bad-k')                       ! (i): k_osc = 0
         call setup_pair()
         call set_scheme([5])
         par%k_osc = 0.0d0
         call surface_oscillator_init(par)
      case ('bad-nat')                     ! (j): three-atom carrier
         call setup_trio()
         call set_scheme([5])
         call surface_oscillator_init(par)
      case ('bad-degen')                   ! (k): coincident atoms
         call setup_degen()
         call set_scheme([5])
         call surface_oscillator_init(par)
      case ('bad-bind0')                   ! (l): no carrier
         call setup_pair()
         call set_scheme([1])
         call surface_oscillator_init(par)
      case ('bad-bind2')                   ! (m): two carriers
         call setup_two()
         call set_scheme([5, 5])
         call surface_oscillator_init(par)
      case ('bad-noinit')                  ! (n): sample without init
         call setup_pair()
         call set_scheme([5])
         call surface_oscillator_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! equilibrium data of the standard carrier
   q_eq = pair_geom()
   r_eq = norm2(q_eq(4:6) - q_eq(1:3))
   u_bond = (q_eq(4:6) - q_eq(1:3))/r_eq
   com_eq = (m_atom*q_eq(1:3) + m_anch*q_eq(4:6))/(m_atom + m_anch)

   ! ---------- (a) quantum-mode energy exactness ----------
   call setup_pair()
   call set_scheme([5])
   call state_create(sta, 2)
   par%n_osc = 1
   par%n_level = 3
   par%k_osc = 42.0d0
   call rng_init(7319)
   call surface_oscillator_init(par)
   ok_all = .true.
   do draw = 1, ndraw_exact
      call surface_oscillator_sample(sta)
      bond = sta%q(4:6) - sta%q(1:3)
      blen = norm2(bond)
      delta = blen - r_eq
      vrel = dot_product(sta%p(4:6)/m_anch - sta%p(1:3)/m_atom, u_bond)  ! d/dt bond length
      e_bp = 0.5d0*k_int*delta*delta + 0.5d0*mu*vrel*vrel
      if (abs(e_bp - e_quant) > 1.0d-10*e_quant) ok_all = .false.
      if (abs(surface_oscillator_last_e() - e_quant) > 1.0d-12*e_quant) ok_all = .false.
   end do
   nfail0 = nfail
   if (ok_all) then
      print '(a,es12.4,a)', '(a) quantum energy exact: PASS (E_bp = 3.5*hbar*omega = ', e_quant, &
         ' internal, 5000 draws, back-projection 1e-10)'
   else
      print '(a)', 'FAIL (a) quantum energy exact: back-projected energy or getter off'
      nfail = nfail + 1
   end if

   ! ---------- (c)/(d)/(e) structure and phase identities (same init) ----------
   dmean = 0.0d0
   kmean = 0.0d0
   ok_all = .true.
   do draw = 1, ndraw_phase
      call surface_oscillator_sample(sta)
      bond = sta%q(4:6) - sta%q(1:3)
      blen = norm2(bond)
      delta = blen - r_eq
      vrel = dot_product(sta%p(4:6)/m_anch - sta%p(1:3)/m_atom, u_bond)
      dmean = dmean + delta
      kmean = kmean + 0.5d0*mu*vrel*vrel
      ! (d) momentum structure
      p_atom = sta%p(1:3)
      p_anch = sta%p(4:6)
      if (norm2(p_atom + p_anch) > 1.0d-14*norm2(p_atom)) ok_all = .false.
      if (norm2(cross(p_atom, u_bond)) > 1.0d-14*norm2(p_atom)) ok_all = .false.
      com_new = (m_atom*sta%q(1:3) + m_anch*sta%q(4:6))/(m_atom + m_anch)
      if (norm2(com_new - com_eq) > 1.0d-12) ok_all = .false.
      ! (e) bond geometry: direction exact, displacement along the bond
      if (norm2(cross(bond/blen, u_bond)) > 1.0d-12) ok_all = .false.
      if (abs(delta) > a_quant*(1.0d0 + 1.0d-12)) ok_all = .false.
   end do
   dmean = dmean/dble(ndraw_phase)
   kmean = kmean/dble(ndraw_phase)
   dvar = 0.5d0*a_quant*a_quant/dble(ndraw_phase)         ! Var(mean of A cos phi)
   if (abs(dmean) > 5.0d0*sqrt(dvar)) ok_all = .false.
   ksig = e_quant/sqrt(8.0d0*dble(ndraw_phase))           ! Var(E sin^2 phi) = E^2/8
   if (abs(kmean - 0.5d0*e_quant) > 5.0d0*ksig) ok_all = .false.
   nfail0 = nfail
   if (ok_all) then
      print '(a,f9.5,a,f8.4,a)', '(c)-(e) structure/phase: PASS (mean delta = ', dmean, &
         ', <KE> = ', kmean, ' vs E/2; p = +/-mu v u_bond, COM fixed, bond direction exact)'
   else
      print '(a,f9.5,a,f8.4,a)', 'FAIL (c)-(e): phase or structure identity broke (mean delta ', dmean, &
         ', <KE> ', kmean, ')'
      nfail = nfail + 1
   end if
   call state_destroy(sta)

   ! ---------- (b) thermal-mode exponential energy ----------
   call setup_pair()
   call set_scheme([5])
   call state_create(sta, 2)
   par%n_osc = 0
   par%t_osc = 180.0d0
   call surface_oscillator_init(par)
   call rng_init(7321)                       ! fresh seed: rejection-free draw, but the
   ! stream window differs from (a) - per-criterion reseeding statute
   freq = 0
   amean = 0.0d0
   do draw = 1, ndraw_stat
      call surface_oscillator_sample(sta)
      e_bp = surface_oscillator_last_e()
      amean = amean + e_bp
      cdf = 1.0d0 - exp(-e_bp/(kb_code*180.0d0))
      ib = min(nbin, max(1, int(cdf*dble(nbin)) + 1))
      freq(ib) = freq(ib) + 1
   end do
   amean = amean/dble(ndraw_stat)
   chi2 = 0.0d0
   do ib = 1, nbin
      chi2 = chi2 + (dble(freq(ib)) - ndraw_stat/dble(nbin))**2/(ndraw_stat/dble(nbin))
   end do
   emean = kb_code*180.0d0
   sig = emean/sqrt(dble(ndraw_stat))
   nfail0 = nfail
   if (chi2 > 16.92d0) then
      print '(a,f9.3,a)', 'FAIL (b) thermal chi-square: chi2 = ', chi2, ' > 16.92 (df=9, alpha=0.05)'
      nfail = nfail + 1
   else if (abs(amean - emean) > 5.0d0*sig) then
      print '(a,es12.4,a,es12.4)', 'FAIL (b) thermal mean: sample ', amean, ' vs kb*T = ', emean
      nfail = nfail + 1
   else
      print '(a,f7.3,a,es11.4,a)', '(b) thermal energy: PASS (chi2 = ', chi2, ' <= 16.92, df = 9; mean ', &
         amean, ' = kb*T, N = 50000)'
   end if
   call state_destroy(sta)

   ! ---------- (f)-(n) loud-stop criteria (child processes) ----------
   call expect_loudstop('bad-nosc',   'surface_oscillator_init',   'n_osc',   '(f) bad energy-mode selector loud-stop')
   call expect_loudstop('bad-tosc',   'surface_oscillator_init',   't_osc',   '(g) non-positive temperature loud-stop')
   call expect_loudstop('bad-nlev',   'surface_oscillator_init',   'n_level', '(h) negative quantum number loud-stop')
   call expect_loudstop('bad-k',      'surface_oscillator_init',   'k_osc',   '(i) non-positive force constant loud-stop')
   call expect_loudstop('bad-nat',    'surface_oscillator_init',   'atom count', '(j) non-diatomic carrier loud-stop')
   call expect_loudstop('bad-degen',  'surface_oscillator_init',   'bond',    '(k) degenerate bond loud-stop')
   call expect_loudstop('bad-bind0',  'surface_oscillator_init',   'surface_oscillator scheme', '(l) no-carrier loud-stop')
   call expect_loudstop('bad-bind2',  'surface_oscillator_init',   'two',     '(m) duplicate-carrier loud-stop')
   call expect_loudstop('bad-noinit', 'surface_oscillator_sample', 'surface_oscillator_init', '(n) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_dist_surface_oscillator: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   pure function cross(a, b) result(c)
      real(8), intent(in) :: a(3), b(3)
      real(8) :: c(3)
      c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
   end function cross

   ! atom-anchor pair along a skewed direction, bond 2.3 Ang
   function pair_geom() result(g)
      real(8) :: g(6), u(3)
      u = [ 0.4d0, -0.5d0, sqrt(1.0d0 - 0.16d0 - 0.25d0) ]
      g(1:3) = -2.3d0*u*(m_anch/(m_atom + m_anch))
      g(4:6) = 2.3d0*u*(m_atom/(m_atom + m_anch))
   end function pair_geom

   subroutine setup_pair()
      call fill_roster(2, pair_geom(), [ m_atom, m_anch ])
   end subroutine setup_pair

   subroutine setup_trio()
      real(8) :: g(9)
      g = 0.0d0
      g(1) = -0.6d0
      g(4) = 0.6d0
      g(8) = 0.5d0
      call fill_roster(3, g, [ 1.0d0, 1.0d0, 1.0d0 ])
   end subroutine setup_trio

   subroutine setup_degen()
      real(8) :: g(6)
      g = 0.35d0
      call fill_roster(2, g, [ m_atom, m_anch ])
   end subroutine setup_degen

   subroutine setup_two()
      real(8) :: g(6)
      integer :: k
      g = pair_geom()
      call roster_clear()
      allocate (roster%mass(4))
      roster%mass = [ m_atom, m_anch, m_atom, m_anch ]
      allocate (roster%frag(2))
      do k = 1, 2
         roster%frag(k)%nat = 2
         roster%frag(k)%mass = m_atom + m_anch
         allocate (roster%frag(k)%list(2), roster%frag(k)%qz(6))
         roster%frag(k)%list = [ 2*k - 1, 2*k ]
         roster%frag(k)%qz = g
      end do
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine setup_two

   ! mini-roster by direct field fill: one fragment
   subroutine fill_roster(nat, g, mss)
      integer, intent(in) :: nat
      real(8), intent(in) :: g(:), mss(:)
      integer :: k
      call roster_clear()
      allocate (roster%mass(nat))
      roster%mass = mss
      allocate (roster%frag(1))
      roster%frag(1)%nat = nat
      roster%frag(1)%mass = sum(mss)
      allocate (roster%frag(1)%list(nat), roster%frag(1)%qz(3*nat))
      roster%frag(1)%list = [ (k, k = 1, nat) ]
      roster%frag(1)%qz = g
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine fill_roster

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
      cmd = 'build'//sep//'check_dist_surface_oscillator.exe '//trim(cmode)// &
            ' 2>build'//sep//'surf_osc_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/surf_osc_err.txt', status='old', action='read', iostat=ios)
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

end program check_dist_surface_oscillator
