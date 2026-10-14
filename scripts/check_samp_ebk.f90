!=====================================================================
! check_samp_ebk.f90 - ebk distribution-member criteria program (criteria
!                       defined FIRST in this header; this program executes
!                       them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files
! out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 &
!            globals/config_atoms.f90 globals/state.f90 &
!            utils/rng.f90 utils/geometry.f90 utils/linalg.f90 utils/hessian.f90 utils/spectrum.f90 &
!            interface/elec_interface.f90 interface/force_interface.f90 interface/spectrum_interface.f90 &
!            dynamics/propagator.f90 utils/geom_init.f90 interface/bath.f90 &
!            interface/sampler.f90 methods/samp_ebk.f90 &
!            ../scripts/check_samp_ebk.f90 -o build/check_samp_ebk.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter rule: the member is driven by directly constructed
! parameters (spectroscopic constants, a directly filled two-atom list_atoms);
! the input grammar is not in the loop. The theory anchors:
!   - a Morse oscillator built from (w_e, w_ex_e, mu) has its WKB action
!     exactly quantized: the Bohr-Sommerfeld levels are EXACTLY
!     E_n = w_e*(n+1/2) - w_ex_e*(n+1/2)^2 [cm^-1] for J = 0;
!   - every accepted draw must satisfy the classical energy identity
!     pr^2/(2*mu) + V_eff(r) = e_level exactly (algebraic anchor).
! System: HCl-like synthetic constants w_e = 2990.9 cm^-1, w_ex_e = 52.8 cm^-1,
!   b_rot = 10.59 cm^-1, masses 1.008/35.45 amu, r_e = 1.275 Angstrom (list_atoms
!   geometry), state (n = 3, J = 2).
! Named-abort verification approach (same scheme as the sibling checks): every
! named-abort criterion runs as a CHILD PROCESS - the check re-spawns itself
! with a mode argument; the child drives the REAL member entry into the stop
! and must print a naming message on stderr and exit nonzero.
! Criteria:
!   (a) action quantization, J = 0 Morse exactness: with j_rot = 0 the
!       converged level energy equals the analytic Dunham form
!       wvn2ev*ev_dyn*(w_e*x - w_ex_e*x^2), x = n+1/2, within 5e-4 relative
!       (the expected floor is the turning-point-inset action sliver,
!       order 1e-4 - the fixed point shifts e to quantize the INSET
!       interval)
!   (b) action quantization, rotating state: for (n=3, J=2) the action
!       integral (1/(pi*hbar)) * integral of pr dr over the member's
!       sampled interval, recomputed HERE with an independent 100000-panel
!       composite Simpson rule from the member's published Morse
!       parameters, equals n+1/2 = 3.5 within 1e-3
!   (c) per-draw classical energy identity, reconstructed FROM THE STATE:
!       r from the written bond length, pr from p.(bond axis): |pr^2/(2*mu)
!       + V_eff(r) - e_level| <= 1e-9*e_level for every draw (20000 draws)
!   (d) written-state structure: fragment center of mass at the origin,
!       total momentum zero, |L| = sqrt(J(J+1))*hbar within 1e-9, L
!       perpendicular to the bond axis within 1e-12
!   (e) r-distribution shape: the internuclear distance is distributed with
!       density proportional to 1/|pr| on the sampled interval; 20 bins,
!       200000 draws, expected fractions from an independent 20000-panel
!       trapezoid CDF; chi-square df = 20-1 = 19, critical value 30.14
!       (alpha = 0.05)
!   (f) zero anharmonic constant named abort (child): w_ex_e = 0 -> stderr
!       names ebk_init and the field
!   (g) negative quantum number named abort (child): n_vib = -1 -> stderr
!       names ebk_init
!   (h) non-diatomic fragment named abort (child): a 3-atom carrier -> stderr
!       names ebk_init and the atom count
!   (i) fragment binding, absent named abort (child): no list_atoms fragment
!       carries the ebk scheme -> stderr names ebk_init
!   (j) level above the dissociation asymptote named abort (child):
!       n_vib = 20, j_rot = 30 (rotational term pushes the level above
!       D_e) -> stderr names ebk_init and the asymptote
!   (k) sample-before-init named abort (child): ebk_sample with no ebk_init
!       -> stderr names ebk_draw (the fused wrapper delegates; the two-phase
!       split renamed the abort sites - draw fires first)
!   (l) dual-carrier equal treatment (positive): two identical diatomic
!       carriers both carry the word (scheme [2,2]) -> one parameter set
!       binds both (the carrier-2 Morse set and level equal the carrier-1
!       set, 1e-12) and one sample realizes EACH carrier: both fragments'
!       q/p come off the sentinel and each carrier's per-draw energy
!       identity pr^2/(2*mu) + V_eff(r) = e_level holds to 1e-9 relative
!       (r/pr from ebk_last(k), carrier 2 explicit)
!=====================================================================
program check_samp_ebk
   use consts,        only: pi, hbar_act, wvn2ev, ev_dyn
   use rng,           only: rng_init
   use config,        only: reactants
   use config_atoms, only: list_atoms, list_atoms_clear
   use state,         only: state_t, state_create, state_destroy
   use samp_ebk,      only: ebk_car_t, ebk_params_t, ebk_init, ebk_sample, ebk_e_level, &
                            ebk_r_turn, ebk_morse, ebk_last
   implicit none
   type(state_t) :: sta
   type(ebk_params_t) :: par
   integer, parameter :: ndraw = 200000, ncheck = 20000, nbin = 20
   integer :: nfail = 0, nfail0, draw, ib, i
   integer :: hist(nbin)
   real(8) :: a_m, d_e, mu, r_e, al_m, e_lvl, r1, r2, last_r, last_pr
   real(8) :: cv, x, e_dunham, chi2, dum, rr2, s_w
   real(8) :: cdf(0:nbin)
   real(8) :: r_sta, pr_sta, v_eff, e_id, l_vec(3), q1(3), q2(3), p1(3), p2(3)
   real(8) :: ax(3), com(3), ptot(3)
   real(8) :: a2, d2, mu2, re2, al2
   logical :: ok_id, ok_struct
   character(len=32) :: mode
   real(8), parameter :: tol_a1 = 5.0d-4, tol_a2 = 1.0d-3
   real(8), parameter :: sentinel = -777.0d0

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   allocate (par%car(1))   ! single-carrier check face (component defaults stand)
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 2)
      select case (trim(mode))
      case ('bad-exe')                    ! (f): w_ex_e = 0
         call setup_diatom()
         call set_scheme([2])
         par = hcl_params()
         par%car(1)%w_ex_e = 0.0d0
         call ebk_init(par)
      case ('bad-n')                      ! (g): n_vib = -1
         call setup_diatom()
         call set_scheme([2])
         par = hcl_params()
         par%car(1)%n_vib = -1
         call ebk_init(par)
      case ('bad-trimer')                 ! (h): 3-atom carrier
         call setup_trimer()
         call set_scheme([2])
         par = hcl_params()
         call ebk_init(par)
      case ('bad-bind')                   ! (i): no carrier
         call setup_diatom()
         call set_scheme([1])
         par = hcl_params()
         call ebk_init(par)
      case ('bad-high')                   ! (j): level above the asymptote
         call setup_diatom()
         call set_scheme([2])
         par = hcl_params()
         par%car(1)%n_vib = 20
         par%car(1)%j_rot = 30
         call ebk_init(par)
      case ('bad-noinit')                 ! (k): sample without init
         call setup_diatom()
         call set_scheme([2])
         call ebk_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   cv = wvn2ev*ev_dyn              ! internal energy per cm^-1 [amu*Ang^2/fs^2]
   call setup_diatom()
   call set_scheme([2])
   call state_create(sta, 2)
   call rng_init(4219)

   ! ---------- (a) J = 0 Morse-WKB exactness ----------
   par = hcl_params()
   par%car(1)%j_rot = 0
   call ebk_init(par)
   e_lvl = ebk_e_level()
   x = dble(par%car(1)%n_vib) + 0.5d0
   e_dunham = cv*(par%car(1)%w_e*x - par%car(1)%w_ex_e*x*x)
   nfail0 = nfail
   if (abs(e_lvl - e_dunham) > tol_a1*abs(e_dunham)) then
      print '(a,es12.4,a,es12.4,a,es9.2)', 'FAIL (a) J=0 exactness: level ', e_lvl, &
         ' vs Dunham ', e_dunham, ' outside 5e-4 relative'
      nfail = nfail + 1
   else
      print '(a,es12.4,a,es12.4,a,es9.2)', '(a) J=0 Morse-WKB exactness: PASS (level ', e_lvl, &
         ' vs analytic ', e_dunham, ', rel dev = inset-floor scale)'
   end if
   if (nfail == nfail0) then
      dum = abs(e_lvl - e_dunham)/abs(e_dunham)
      print '(a,es9.2,a)', '(a) relative deviation ', dum, ' (tolerance 5e-4; the inset sliver floor)'
   end if

   ! ---------- main state (n=3, J=2) ----------
   par = hcl_params()
   call ebk_init(par)
   e_lvl = ebk_e_level()
   call ebk_morse(1, a_m, d_e, mu, r_e, al_m)
   call ebk_r_turn(1, r1, r2)
   x = dble(par%car(1)%n_vib) + 0.5d0

   ! ---------- (b) action quantization of the rotating state ----------
   ! independent 100000-panel composite Simpson rule over the member's sampled
   ! interval (an independent quadrature - not the member's own rule)
   dum = 0.0d0
   do i = 0, 200000
      rr2 = r1 + (r2 - r1)*dble(i)/dble(200000)
      s_w = 2.0d0
      if (i == 0 .or. i == 200000) s_w = 1.0d0
      if (mod(i, 2) == 1) s_w = 4.0d0
      dum = dum + s_w*sqrt(max(0.0d0, 2.0d0*mu*(e_lvl - v_eff_at(rr2))))
   end do
   dum = dum*(r2 - r1)/dble(200000)/3.0d0/(pi*hbar_act)
   nfail0 = nfail
   if (abs(dum - x) > tol_a2) then
      print '(a,f10.6,a,f6.2,a)', 'FAIL (b) rotating action: ', dum, ' vs n+1/2 = ', x, &
         ' outside 1e-3'
      nfail = nfail + 1
   else
      print '(a,f10.6,a,f6.2,a)', '(b) rotating-state action: PASS (J/pi*hbar = ', dum, &
         ' vs n+1/2 = ', x, ', tol 1e-3)'
   end if

   ! ---------- (c)/(d) per-draw identity + written-state structure ----------
   ok_id = .true.
   ok_struct = .true.
   hist = 0
   do draw = 1, ndraw
      call ebk_sample(sta)
      call ebk_last(1, last_r, last_pr)
      if (draw <= ncheck) then
         q1 = sta%q(1:3); q2 = sta%q(4:6)
         p1 = sta%p(1:3); p2 = sta%p(4:6)
         ax = (q2 - q1)/norm2(q2 - q1)          ! bond axis from the written state
         r_sta = norm2(q2 - q1)
         pr_sta = dot_product(p2, ax)
         v_eff = v_eff_at(r_sta)
         e_id = pr_sta**2/(2.0d0*mu) + v_eff
         com = (1.008d0*q1 + 35.45d0*q2)/(1.008d0 + 35.45d0)
         ptot = p1 + p2
         l_vec = cross(q1, p1) + cross(q2, p2)
         if (abs(e_id - e_lvl) > 1.0d-9*e_lvl .or. abs(pr_sta + dot_product(p1, ax)) > 1.0d-9) &
            ok_id = .false.
         if (norm2(com) > 1.0d-12 .or. norm2(ptot) > 1.0d-12 .or. &
             abs(norm2(l_vec) - al_m) > 1.0d-9*al_m .or. &
             abs(dot_product(l_vec, ax)) > 1.0d-12) ok_struct = .false.
      end if
      ib = min(nbin, max(1, int(nbin*(last_r - r1)/(r2 - r1)) + 1))
      hist(ib) = hist(ib) + 1
   end do
   nfail0 = nfail
   if (.not. ok_id) then
      print '(a)', 'FAIL (c) energy identity: pr^2/2mu + V_eff(r) /= e_level in some draw'
      nfail = nfail + 1
   else
      print '(a)', '(c) per-draw energy identity: PASS (20000 draws, 1e-9 relative, reconstructed from state)'
   end if
   if (.not. ok_struct) then
      print '(a)', 'FAIL (d) structure: COM/total momentum/|L|/perpendicularity violated'
      nfail = nfail + 1
   else
      print '(a)', '(d) written-state structure: PASS (COM 0, total p 0, |L| = sqrt(J(J+1))hbar, L perp bond)'
   end if

   ! ---------- (e) r-distribution shape (density prop. to 1/|pr|) ----------
   cdf(0) = 0.0d0
   do ib = 1, nbin
      cdf(ib) = trapz_cdf(r1 + (r2 - r1)*dble(ib)/dble(nbin))
   end do
   cdf = cdf/cdf(nbin)             ! normalize by the full integral (cancellation-free bin shares)
   chi2 = 0.0d0
   do ib = 1, nbin
      dum = ndraw*(cdf(ib) - cdf(ib-1))
      chi2 = chi2 + (dble(hist(ib)) - dum)**2/dum
   end do
   nfail0 = nfail
   if (chi2 > 30.14d0) then
      print '(a,f9.3,a)', 'FAIL (e) r-shape chi-square: chi2 = ', chi2, ' > 30.14 (df=19, alpha=0.05)'
      nfail = nfail + 1
   else
      print '(a,f7.3,a)', '(e) r-distribution 1/|pr| shape: PASS (chi2 = ', chi2, &
         ' <= 30.14, df = 19, alpha = 0.05, N = 200000)'
   end if
   call state_destroy(sta)

   ! ---------- (l) dual-carrier equal treatment: one Morse parameter set,
   ! two carriers - each carrier realizes its own draw ----------
   call setup_two_diatoms()
   call set_scheme([2, 2])
   call state_create(sta, 4)
   call rng_init(8677)
   par = hcl_params()
   block                                  ! duplicate the row onto the second
      type(ebk_car_t) :: tmp              ! carrier (identical rows - the equal-
      tmp = par%car(1)                    ! treatment face of the per-carrier form)
      deallocate (par%car)
      allocate (par%car(2))
      par%car(1) = tmp
      par%car(2) = tmp
   end block
   call ebk_init(par)
   e_lvl = ebk_e_level(1)
   call ebk_morse(1, a_m, d_e, mu, r_e, al_m)
   call ebk_morse(2, a2, d2, mu2, re2, al2)
   ok_id = .true.
   sta%q = sentinel
   sta%p = sentinel
   do draw = 1, ncheck
      call ebk_sample(sta)
      if (draw == 1 .and. (any(sta%q == sentinel) .or. any(sta%p == sentinel))) ok_id = .false.
      call ebk_last(1, last_r, last_pr)             ! carrier 1 identity
      e_id = last_pr**2/(2.0d0*mu) + v_eff_at(last_r)
      if (abs(e_id - e_lvl) > 1.0d-9*e_lvl) ok_id = .false.
      call ebk_last(2, last_r, last_pr)             ! carrier 2 identity
      e_id = last_pr**2/(2.0d0*mu) + v_eff_at(last_r)
      if (abs(e_id - e_lvl) > 1.0d-9*e_lvl) ok_id = .false.
   end do
   nfail0 = nfail
   if (.not. ok_id .or. abs(ebk_e_level(2) - e_lvl) > 1.0d-12*e_lvl .or. &
       abs(a2 - a_m) > 1.0d-12*a_m .or. abs(d2 - d_e) > 1.0d-12*d_e .or. &
       abs(mu2 - mu) > 1.0d-12*mu .or. abs(re2 - r_e) > 1.0d-12*r_e .or. &
       abs(al2 - al_m) > 1.0d-12*max(al_m, 1.0d-300)) then
      print '(a)', 'FAIL (l) dual-carrier equal treatment: a carrier identity broke or the sets differ'
      nfail = nfail + 1
   else
      print '(a)', '(l) dual-carrier equal treatment: PASS (both fragments written; each carrier'// &
         ' realizes pr^2/2mu + V_eff = e_level; one Morse set, equal levels)'
   end if
   call state_destroy(sta)

   ! ---------- (f)-(k) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-exe',     'ebk_init',   'w_ex_e',       '(f) zero-anharmonic loud-stop')
   call expect_loudstop('bad-n',       'ebk_init',   'n_vib',        '(g) negative-quantum-number loud-stop')
   call expect_loudstop('bad-trimer',  'ebk_init',   'diatomic',     '(h) non-diatomic-carrier loud-stop')
   call expect_loudstop('bad-bind',    'ebk_init',   'ebk',          '(i) no-fragment-binding loud-stop')
   call expect_loudstop('bad-high',    'ebk_init',   'asymptote',    '(j) above-asymptote loud-stop')
   call expect_loudstop('bad-noinit',  'ebk_draw', 'ebk_init',     '(k) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_samp_ebk: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! V_eff(r) rebuilt HERE from the member's published Morse set (independent
   ! of the member's private copy): Morse from the well bottom + centrifugal
   real(8) function v_eff_at(r)
      real(8), intent(in) :: r
      v_eff_at = d_e*(1.0d0 - exp(-a_m*(r - r_e)))**2 + al_m**2/(2.0d0*mu*r*r)
   end function v_eff_at

   ! trapezoid CDF of the density 1/sqrt(2*mu*(e_level - V_eff)) up to r_edge
   ! (20000 panels over the whole interval; the node density of the whole-grid
   ! integral serves every edge by accumulation - recomputed per call for
   ! simplicity, the grid is cheap)
   real(8) function trapz_cdf(r_edge)
      real(8), intent(in) :: r_edge
      integer, parameter :: npan = 20000
      real(8) :: r_lo, r_hi, w, f_lo, f_hi, rr, acc
      integer :: k
      acc = 0.0d0
      w = (r2 - r1)/dble(npan)
      r_lo = r1
      f_lo = 1.0d0/sqrt(max(1.0d-300, 2.0d0*mu*(e_lvl - v_eff_at(r_lo))))
      do k = 1, npan
         rr = min(r_edge, r1 + w*dble(k))
         f_hi = 1.0d0/sqrt(max(1.0d-300, 2.0d0*mu*(e_lvl - v_eff_at(rr))))
         acc = acc + 0.5d0*(f_lo + f_hi)*(rr - r_lo)
         if (rr >= r_edge) exit
         r_lo = rr
         f_lo = f_hi
      end do
      trapz_cdf = acc
   end function trapz_cdf

   function cross(a, b) result(c)
      real(8), intent(in) :: a(3), b(3)
      real(8) :: c(3)
      c = (/ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) /)
   end function cross

   ! HCl-like synthetic parameter set
   function hcl_params() result(p)
      type(ebk_params_t) :: p
      allocate (p%car(1))   ! single-carrier check face (component defaults stand)
      p%car(1)%n_vib   = 3
      p%car(1)%j_rot   = 2
      p%car(1)%w_e     = 2990.9d0     ! [cm^-1]
      p%car(1)%w_ex_e  = 52.8d0       ! [cm^-1]
      p%car(1)%b_rot   = 10.59d0      ! [cm^-1]
      p%car(1)%ai_rot  = 1.5918d0     ! [amu Ang^2]
   end function hcl_params

   ! two-atom list_atoms by direct field fill: H(1.008) at origin, Cl(35.45) on +x
   subroutine setup_diatom()
      call list_atoms_clear()
      allocate (list_atoms%mass(2))
      list_atoms%mass = [ 1.008d0, 35.45d0 ]
      allocate (list_atoms%frag(1))
      list_atoms%frag(1)%nat = 2
      list_atoms%frag(1)%mass = 36.458d0
      allocate (list_atoms%frag(1)%list(2), list_atoms%frag(1)%qz(6))
      list_atoms%frag(1)%list = [ 1, 2 ]
      list_atoms%frag(1)%qz = [ 0.0d0, 0.0d0, 0.0d0, 1.275d0, 0.0d0, 0.0d0 ]
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine setup_diatom

   ! two-fragment list_atoms (the dual-carrier form): two identical HCl-like
   ! diatomics, fragment 2 on atoms 3/4
   subroutine setup_two_diatoms()
      integer :: k
      call list_atoms_clear()
      allocate (list_atoms%mass(4))
      list_atoms%mass = [ 1.008d0, 35.45d0, 1.008d0, 35.45d0 ]
      allocate (list_atoms%frag(2))
      do k = 1, 2
         list_atoms%frag(k)%nat = 2
         list_atoms%frag(k)%mass = 36.458d0
         allocate (list_atoms%frag(k)%list(2), list_atoms%frag(k)%qz(6))
         list_atoms%frag(k)%list = [ 2*k - 1, 2*k ]
         list_atoms%frag(k)%qz = [ 0.0d0, 0.0d0, 0.0d0, 1.275d0, 0.0d0, 0.0d0 ]
      end do
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine setup_two_diatoms

   ! three-atom list_atoms (the non-diatomic carrier for (h))
   subroutine setup_trimer()
      call list_atoms_clear()
      allocate (list_atoms%mass(3))
      list_atoms%mass = 1.0d0
      allocate (list_atoms%frag(1))
      list_atoms%frag(1)%nat = 3
      list_atoms%frag(1)%mass = 3.0d0
      allocate (list_atoms%frag(1)%list(3), list_atoms%frag(1)%qz(9))
      list_atoms%frag(1)%list = [ 1, 2, 3 ]
      list_atoms%frag(1)%qz = 0.0d0
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine setup_trimer

   ! preset the per-fragment scheme codes
   subroutine set_scheme(codes)
      integer, intent(in) :: codes(:)
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
      allocate (reactants%dist_scheme(size(codes)))
      reactants%dist_scheme = codes
   end subroutine set_scheme

   ! spawn the child for one bad case and assert the named abort
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
      cmd = 'build'//sep//'check_samp_ebk.exe '//trim(cmode)//' 2>build'//sep//'ebk_err.txt'
      if (windows_host()) then
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      else
         ! Linux: the exe suffix is absent - build the matching child name
         block
            character(len=512) :: cmd2
            cmd2 = 'build/check_samp_ebk '//trim(cmode)//' 2>build/ebk_err.txt'
            call execute_command_line(trim(cmd2), wait=.true., exitstat=stat)
         end block
      end if
      open (newunit=v, file='build/ebk_err.txt', status='old', action='read', iostat=ios)
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

   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_samp_ebk
