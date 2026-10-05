!=====================================================================
! check_inc_pair.f90 - the pair (two-reactant) incident-member criteria
!                       program (criteria defined FIRST in this header;
!                       this program executes them)
! Compile (from src/, linked against the tree's module objects):
!   gfortran -Jbuild/gfortran -c ../scripts/check_inc_pair.f90 -o build/ck_inc_pair.o
!   gfortran build/ck_inc_pair.o <tree objects> -o build/check_inc_pair.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter rule: the member is driven by directly constructed
! beam_params_t and a directly filled list_atoms; e_rel/b_max/r_sep go
! through config%reactants exactly as the member contract prescribes. All
! quantities are RECONSTRUCTED FROM THE WRITTEN STATE - no member getters.
! System: gas-phase pair, projectile mass 1.008 amu (i_proj = 1 row), target
!   mass 35.45 amu (i_targ = 2 row), reduced mass mu = 0.9801 amu.
! Named-abort verification: every abort criterion runs as a CHILD PROCESS
! (execute_command_line re-spawn); the child drives the REAL member entry
! into the stop. Stream note: statistical criteria RE-SEED a fixed seed
! (7319 everywhere except (f), which draws from 7321 - see the legacy
! header note about the shuffle-table generator's window wobble).
! Criteria:
!   (a) b^2 uniform on [0, b_max^2] (n_b = 1, fixed theta/chi): 200000
!       draws, 20 bins, chi-square df = 19, critical 30.14 (alpha = 0.05)
!   (b) cos(theta) linear (uniform) on [cos(theta_max), 1]: 20 bins,
!       df = 19, critical 30.14, from the projectile momentum direction
!   (c) azimuth chi uniform on [0, 2*pi): 20 bins, df = 19, critical 30.14
!   (d) fixed collision energy (mode 1): E = |p1|^2/(2*mu) equals
!       e_conv*e_rel to 1e-9 relative (20000 draws)
!   (e) thermal collision energy (mode 0, t_trans = 300 K): mean/variance
!       within 5 sigma of the gamma(2) moments (200000 draws)
!   (f) density-weighted collision speed (mode 2, seed 7321): equal-CDF
!       bins against an independent 20000-panel trapezoid CDF, df = 19,
!       critical 30.14
!   (g) fixed-mode exactness battery: b = b_max, cos(theta), chi, E all
!       exact (1e-12) on 20000 draws
!   (h) placement + momentum structure (fixed modes): target COM at the
!       drop point, projectile COM at b(cos,si) + r_sep*sin(theta)*(cos,si)
!       with height r_sep exactly, p1 = mu*v*u_inc, p2 = -p1 (1e-12)
!   (j) unit-cell drop (n_aim = 1, box 3x3, skew 60 deg): drop points fill
!       the parallelogram, both generator margials uniform (pooled chi2,
!       df = 38, critical 53.38)
!   (t) role mirror: swapped role rows exchange the absolute placement,
!       same site and relative momentum (1e-12)
!   (u) stale-carry wipe: a second assembly on an unwiped state reproduces
!       the fresh placement and COM momenta exactly (the assignment, never
!       accumulation, statute)
!   (k)-(r) named aborts (children): illegal energy mode / mode 0 without
!       temperature / mode 2 without width / illegal polar mode /
!       degenerate sampling range / missing unit cell / three-fragment
!       geometry / draw-before-init (the beam_laws guard site)
!   (w) surface-paradigm rejection (child): N_AIM = 2 (the fixed-site drop)
!       aborts in the pair member - it is a surface-paradigm option
!=====================================================================
program check_inc_pair
   use consts,       only: pi, two_pi, e_conv, r_kcal
   use rng,          only: rng_init
   use config,       only: reactants
   use config_atoms, only: list_atoms, list_atoms_clear
   use state,        only: state_t, state_create, state_destroy
   use beam_laws,    only: beam_params_t
   use inc_pair,     only: inc_pair_init, inc_pair_sample
   implicit none
   type(state_t) :: sta
   type(beam_params_t) :: par
   integer, parameter :: ndraw = 200000, ncheck = 20000, nbin = 20
   integer :: nfail = 0, draw, i, k
   integer :: hist(nbin), edge_bin(nbin + 1)
   real(8) :: mu, m1, m2, v, vhat(3), costh, sinth, chi_v, b_rec, e_rec
   real(8) :: lat(2), chi2, dum, sd
   real(8) :: amean, avar, rx0, ry0, s1, s2
   real(8) :: v0, vwid, cdfv(0:20000), vv, zedge(nbin + 1)
   real(8) :: p1v(3), p2v(3), lat_exp(2), h_exp(3)
   real(8) :: q_ref(3), p_ref(3)
   logical :: ok_exact, ok_place
   character(len=32) :: mode
   real(8), parameter :: thta_fix = 1.1d0, chi_fix = 0.9d0

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call setup_gas(1.008d0, 35.45d0)
      call set_react(5.0d0, 3.0d0, 8.0d0)
      par = base_params()
      select case (trim(mode))
      case ('bad-erel3')                  ! (k): n_e_rel = 3
         par%n_e_rel = 3
         call inc_pair_init(par)
      case ('bad-ttrans')                 ! (l): mode 0 without a temperature
         par%n_e_rel = 0
         par%t_trans = 0.0d0
         call inc_pair_init(par)
      case ('bad-vwidth')                 ! (m): mode 2 without a width
         par%n_e_rel = 2
         par%v_width = 0.0d0
         call inc_pair_init(par)
      case ('bad-nthta')                  ! (n): n_thta = 2
         par%n_thta = 2
         call inc_pair_init(par)
      case ('bad-thtamax')                ! (o): degenerate sampling range
         par%n_thta = 1
         par%thta_max = 0.0d0
         call inc_pair_init(par)
      case ('bad-cell')                   ! (p): n_aim = 1 with no cell buffered
         par%n_aim = 1
         call inc_pair_init(par)
      case ('bad-geom')                   ! (q): three gas-phase fragments
         call setup_trimer_gas()
         call inc_pair_init(par)
      case ('bad-noinit')                 ! (r): sample without init
         call state_create(sta, 2)
         call fresh_state(sta)
         call inc_pair_sample(sta)
      case ('bad-aim2')                   ! (w): the surface fixed-site mode
         par%n_aim = 2
         par%aim_x = 1.0d0
         par%aim_y = 1.0d0
         par%aim_x_given = .true.
         par%aim_y_given = .true.
         call inc_pair_init(par)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   m1 = 1.008d0
   m2 = 35.45d0
   mu = m1*m2/(m1 + m2)
   call setup_gas(m1, m2)
   call set_react(5.0d0, 3.0d0, 8.0d0)
   call state_create(sta, 2)
   call rng_init(7319)

   ! ---------- (a) + (d): b^2 uniform with a fixed-mode backdrop ----------
   par = base_params()
   par%n_thta = 0
   par%thta_max = 0.7853981633974483d0   ! pi/4 fixed beam
   par%n_chi = 0
   par%chi = 0.7d0
   par%n_b = 1
   call inc_pair_init(par)
   hist = 0
   ok_exact = .true.
   do draw = 1, ndraw
      call fresh_state(sta)
      call inc_pair_sample(sta)
      if (draw <= ncheck) then
         e_rec = dot_product(sta%p(1:3), sta%p(1:3))/(2.0d0*mu)/e_conv
         if (abs(e_rec - 5.0d0) > 1.0d-9*5.0d0) ok_exact = .false.
      end if
      lat = sta%q(1:2) - sta%q(4:5) - 8.0d0*sin(0.7853981633974483d0)* &
            [ cos(0.7d0), sin(0.7d0) ]
      b_rec = sqrt(dot_product(lat, lat))
      hist(min(nbin, max(1, int(nbin*b_rec*b_rec/9.0d0) + 1))) = &
         hist(min(nbin, max(1, int(nbin*b_rec*b_rec/9.0d0) + 1))) + 1
   end do
   chi2 = chisq_eq(hist, ndraw)
   if (chi2 > 30.14d0) then
      print '(a,f9.3,a)', 'FAIL (a) b^2 uniform: chi2 = ', chi2, ' > 30.14 (df=19, alpha=0.05)'
      nfail = nfail + 1
   else
      print '(a,f7.3,a)', '(a) b^2 uniform on [0,b_max^2]: PASS (chi2 = ', chi2, &
         ' <= 30.14, df = 19, alpha = 0.05, N = 200000)'
   end if
   if (.not. ok_exact) then
      print '(a)', 'FAIL (d) fixed energy: reconstructed E deviates from e_conv*e_rel'
      nfail = nfail + 1
   else
      print '(a)', '(d) fixed collision energy: PASS (20000 draws, 1e-9 relative, from |p1|^2/2mu)'
   end if

   ! ---------- (b) cos(theta) linear ----------
   call inc_pair_init(cos_mode_params())
   call rng_init(7319)
   hist = 0
   do draw = 1, ndraw
      call fresh_state(sta)
      call inc_pair_sample(sta)
      v = sqrt(dot_product(sta%p(1:3), sta%p(1:3)))
      costh = -sta%p(3)/v
      hist(min(nbin, max(1, int(nbin*costh) + 1))) = hist(min(nbin, max(1, int(nbin*costh) + 1))) + 1
   end do
   chi2 = chisq_eq(hist, ndraw)
   if (chi2 > 30.14d0) then
      print '(a,f9.3,a)', 'FAIL (b) cos(theta) linear: chi2 = ', chi2, ' > 30.14 (df=19, alpha=0.05)'
      nfail = nfail + 1
   else
      print '(a,f7.3,a)', '(b) cos(theta) linear on [0,1]: PASS (chi2 = ', chi2, &
         ' <= 30.14, df = 19, alpha = 0.05, N = 200000)'
   end if

   ! ---------- (c) azimuth uniform ----------
   call inc_pair_init(azi_mode_params())
   call rng_init(7319)
   hist = 0
   do draw = 1, ndraw
      call fresh_state(sta)
      call inc_pair_sample(sta)
      v = sqrt(dot_product(sta%p(1:3), sta%p(1:3)))
      chi_v = atan2(-sta%p(2)/v, -sta%p(1)/v)
      if (chi_v < 0.0d0) chi_v = chi_v + two_pi
      hist(min(nbin, max(1, int(nbin*chi_v/two_pi) + 1))) = &
         hist(min(nbin, max(1, int(nbin*chi_v/two_pi) + 1))) + 1
   end do
   chi2 = chisq_eq(hist, ndraw)
   if (chi2 > 30.14d0) then
      print '(a,f9.3,a)', 'FAIL (c) azimuth uniform: chi2 = ', chi2, ' > 30.14 (df=19, alpha=0.05)'
      nfail = nfail + 1
   else
      print '(a,f7.3,a)', '(c) azimuth uniform on [0,2pi): PASS (chi2 = ', chi2, &
         ' <= 30.14, df = 19, alpha = 0.05, N = 200000)'
   end if

   ! ---------- (e) thermal energy moments (mode 0) ----------
   par = base_params()
   par%n_e_rel = 0
   par%t_trans = 300.0d0
   call inc_pair_init(par)
   call rng_init(7319)
   amean = 0.0d0
   avar = 0.0d0
   do draw = 1, ndraw
      call fresh_state(sta)
      call inc_pair_sample(sta)
      e_rec = dot_product(sta%p(1:3), sta%p(1:3))/(2.0d0*mu)/e_conv
      amean = amean + e_rec
      avar = avar + e_rec*e_rec
   end do
   amean = amean/ndraw
   avar = (avar - ndraw*amean*amean)/(ndraw - 1)
   sd = sqrt(2.0d0/dble(ndraw))*r_kcal*300.0d0
   dum = 2.0d0*(r_kcal*300.0d0)**2*sqrt(5.0d0/dble(ndraw))
   if (abs(amean - 2.0d0*r_kcal*300.0d0) > 5.0d0*sd) then
      print '(a,f9.5,a,f9.5)', 'FAIL (e) thermal mean: ', amean, ' vs ', 2.0d0*r_kcal*300.0d0
      nfail = nfail + 1
   else if (abs(avar - 2.0d0*(r_kcal*300.0d0)**2) > 5.0d0*dum) then
      print '(a,f9.5,a,f9.5)', 'FAIL (e) thermal variance: ', avar, ' vs ', 2.0d0*(r_kcal*300.0d0)**2
      nfail = nfail + 1
   else
      print '(a,f9.5,a,f8.5)', '(e) thermal E = r_kcal*gamma(2)*T: PASS (mean ', amean, &
         ' var ', avar, '; 5-sigma bands held)'
   end if

   ! ---------- (f) DENS speed shape (mode 2), equal-CDF bins ----------
   par = base_params()
   par%n_e_rel = 2
   par%v_width = 0.3d0
   call set_react(10.0d0, 3.0d0, 8.0d0)
   call inc_pair_init(par)
   call rng_init(7321)               ! stream note: variable-draw rejection loop
   v0 = sqrt(2.0d0*e_conv*10.0d0/mu)
   vwid = 0.3d0
   cdfv(0) = 0.0d0
   do i = 1, 20000
      vv = 8.0d0*dble(i)/dble(20000)
      cdfv(i) = cdfv(i-1) + 0.5d0*(dens(vv) + dens(8.0d0*dble(i-1)/dble(20000)))*8.0d0/20000.0d0
   end do
   cdfv = cdfv/cdfv(20000)
   do k = 0, nbin
      zedge(k+1) = dble(k)/dble(nbin)
   end do
   edge_bin = 1
   i = 0
   do k = 1, nbin + 1
      do while (i < 20000 .and. cdfv(i+1) < zedge(k))
         i = i + 1
      end do
      edge_bin(k) = i
   end do
   hist = 0
   do draw = 1, ndraw
      call fresh_state(sta)
      call inc_pair_sample(sta)
      v = sqrt(dot_product(sta%p(1:3), sta%p(1:3)))/mu
      do k = 1, nbin
         if (8.0d0*dble(edge_bin(k))/dble(20000) <= v .and. &
             v < 8.0d0*dble(edge_bin(k+1))/dble(20000) + merge(1.0d0, 0.0d0, k == nbin)) then
            hist(k) = hist(k) + 1
            exit
         end if
      end do
   end do
   chi2 = chisq_eq(hist, ndraw)
   if (chi2 > 30.14d0) then
      print '(a,f9.3,a)', 'FAIL (f) DENS shape: chi2 = ', chi2, ' > 30.14 (df=19, alpha=0.05)'
      nfail = nfail + 1
   else
      print '(a,f7.3,a)', '(f) DENS v^3*exp(-((v-v0)/a)^2) shape: PASS (chi2 = ', chi2, &
         ' <= 30.14, df = 19, alpha = 0.05, N = 200000, equal-CDF bins)'
   end if

   ! ---------- (g) fixed-mode exactness battery ----------
   call set_react(5.0d0, 3.0d0, 8.0d0)
   par = base_params()
   par%n_thta = 0
   par%thta_max = thta_fix
   par%n_chi = 0
   par%chi = chi_fix
   par%n_b = 0
   call inc_pair_init(par)
   ok_exact = .true.
   do draw = 1, ncheck
      call fresh_state(sta)
      call inc_pair_sample(sta)
      v = sqrt(dot_product(sta%p(1:3), sta%p(1:3)))
      costh = -sta%p(3)/v
      chi_v = atan2(-sta%p(2)/v, -sta%p(1)/v)
      if (chi_v < 0.0d0) chi_v = chi_v + two_pi
      sinth = sqrt(1.0d0 - costh*costh)
      lat = sta%q(1:2) - sta%q(4:5) - 8.0d0*sinth*[ cos(chi_fix), sin(chi_fix) ]
      b_rec = sqrt(dot_product(lat, lat))
      e_rec = v*v/(2.0d0*mu)/e_conv
      if (abs(b_rec - 3.0d0) > 1.0d-12 .or. abs(costh - cos(thta_fix)) > 1.0d-12 .or. &
          abs(chi_v - chi_fix) > 1.0d-12 .or. abs(e_rec - 5.0d0) > 1.0d-9*5.0d0) ok_exact = .false.
   end do
   if (.not. ok_exact) then
      print '(a)', 'FAIL (g) fixed-mode battery: b/cos(theta)/chi/E not exact'
      nfail = nfail + 1
   else
      print '(a)', '(g) fixed modes: PASS (b = b_max, cos(theta), chi, E all exact on 20000 draws)'
   end if

   ! ---------- (h) placement + momentum structure (gas pair) ----------
   call inc_pair_init(par)                     ! same fixed-mode configuration
   ok_place = .true.
   do draw = 1, ncheck
      call fresh_state(sta)
      call inc_pair_sample(sta)
      p1v = sta%p(1:3)
      p2v = sta%p(4:6)
      v = sqrt(dot_product(p1v, p1v))/mu    ! relative speed (|p1| = mu*v)
      vhat = p1v/(mu*v)
      costh = -vhat(3)
      sinth = sqrt(1.0d0 - costh*costh)
      lat_exp = 3.0d0*[ cos(chi_fix), sin(chi_fix) ] + 8.0d0*sinth*[ cos(chi_fix), sin(chi_fix) ]
      h_exp = (/ lat_exp(1), lat_exp(2), 8.0d0 /)
      if (sqrt(dot_product(sta%q(4:6), sta%q(4:6))) > 1.0d-12 .or. &
          norm2(sta%q(1:3) - h_exp) > 1.0d-12 .or. &
          norm2(p1v + p2v) > 1.0d-12 .or. &
          abs(v - sqrt(2.0d0*e_conv*5.0d0/mu)) > 1.0d-9*v) ok_place = .false.
   end do
   if (.not. ok_place) then
      print '(a)', 'FAIL (h) placement/momentum structure: gas geometry wrong'
      nfail = nfail + 1
   else
      print '(a)', '(h) placement + momentum structure: PASS (target COM at drop, projectile '// &
         'COM formula, p1 = mu*v*u, p2 = -p1)'
   end if

   ! ---------- (j) unit-cell drop (n_aim = 1) ----------
   call setup_gas(m1, m2)
   call state_destroy(sta)
   call state_create(sta, 2)
   list_atoms%box_lx = 3.0d0
   list_atoms%box_ly = 3.0d0
   list_atoms%skew = pi/3.0d0
   par = base_params()
   par%n_aim = 1
   call inc_pair_init(par)
   call rng_init(7319)
   hist = 0
   ok_exact = .true.
   do draw = 1, ndraw
      call fresh_state(sta)
      call inc_pair_sample(sta)
      rx0 = sta%q(4)
      ry0 = sta%q(5)
      s2 = ry0/(3.0d0*sin(pi/3.0d0))
      s1 = (rx0 - s2*3.0d0*cos(pi/3.0d0))/3.0d0
      if (s1 < -1.0d-12 .or. s1 > 1.0d0 + 1.0d-12 .or. s2 < -1.0d-12 .or. s2 > 1.0d0 + 1.0d-12) &
         ok_exact = .false.
      hist(min(nbin, max(1, int(nbin*s1) + 1))) = hist(min(nbin, max(1, int(nbin*s1) + 1))) + 1
      hist(min(nbin, max(1, int(nbin*s2) + 1))) = hist(min(nbin, max(1, int(nbin*s2) + 1))) + 1
   end do
   chi2 = chisq_eq(hist, 2*ndraw)
   if (chi2 > 53.38d0 .or. .not. ok_exact) then
      print '(a,f9.3,a)', 'FAIL (j) unit-cell drop: chi2 = ', chi2, ' or points outside the cell'
      nfail = nfail + 1
   else
      print '(a,f7.3,a)', '(j) unit-cell drop (skew parallelogram): PASS (chi2 = ', chi2, &
         ' <= 53.38, df = 38, alpha = 0.05, both generator marginals pooled)'
   end if

   ! ---------- (t) role mirror ----------
   call setup_gas(m1, m2)
   call state_destroy(sta)
   call state_create(sta, 2)
   par = base_params()
   par%n_thta = 0
   par%thta_max = thta_fix
   par%n_chi = 0
   par%chi = chi_fix
   par%n_b = 0
   call inc_pair_init(par)
   call fresh_state(sta)
   call inc_pair_sample(sta)
   q_ref = sta%q(1:3)
   p_ref = sta%p(1:3)
   list_atoms%i_proj = 2
   list_atoms%i_targ = 1
   call inc_pair_init(par)
   call fresh_state(sta)
   call inc_pair_sample(sta)
   if (norm2(sta%q(4:6) - q_ref) > 1.0d-12 .or. &
       norm2(sta%p(4:6) - p_ref) > 1.0d-12 .or. &
       norm2(sta%q(1:3)) > 1.0d-12 .or. &
       norm2(sta%p(1:3) + p_ref) > 1.0d-12) then
      print '(a)', 'FAIL (t) role mirror: swapped roles did not exchange the placement'
      nfail = nfail + 1
   else
      print '(a)', '(t) role mirror: PASS (swapped role rows: same site and relative '// &
         'momentum, the two fragments exchange absolute placement)'
   end if

   ! ---------- (u) stale-carry wipe (the assignment statute) ----------
   call setup_gas(m1, m2)
   call set_react(5.0d0, 3.0d0, 8.0d0)
   call state_destroy(sta)
   call state_create(sta, 2)
   par = base_params()
   par%n_thta = 0
   par%thta_max = thta_fix
   par%n_chi = 0
   par%chi = chi_fix
   par%n_b = 0
   call inc_pair_init(par)
   call fresh_state(sta)
   call inc_pair_sample(sta)
   q_ref = sta%q(1:3)
   p_ref = sta%p(1:3)
   sta%p = 3.7d0
   sta%q = sta%q + 5.1d0
   call inc_pair_sample(sta)
   if (norm2(sta%q(1:3) - q_ref) > 1.0d-12 .or. &
       norm2(sta%p(1:3) - p_ref) > 1.0d-12 .or. &
       norm2(sta%p(4:6) + p_ref) > 1.0d-12) then
      print '(a)', 'FAIL (u) stale-carry wipe: a repeated assembly accumulated onto the carry'
      nfail = nfail + 1
   else
      print '(a)', '(u) stale-carry wipe: PASS (a second assembly on the unwiped state'// &
         ' reproduces the fresh-state placement and COM momenta exactly)'
   end if
   call state_destroy(sta)

   ! ---------- (k)-(r)+(w) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-erel3',  'inc_pair_init', 'n_e_rel',   '(k) illegal-energy-mode loud-stop')
   call expect_loudstop('bad-ttrans', 'inc_pair_init', 't_trans',   '(l) mode-0-no-temperature loud-stop')
   call expect_loudstop('bad-vwidth', 'inc_pair_init', 'v_width',   '(m) mode-2-no-width loud-stop')
   call expect_loudstop('bad-nthta',  'inc_pair_init', 'n_thta',    '(n) illegal-polar-mode loud-stop')
   call expect_loudstop('bad-thtamax','inc_pair_init', 'thta_max',  '(o) degenerate-range loud-stop')
   call expect_loudstop('bad-cell',   'inc_pair_init', 'unit cell', '(p) missing-unit-cell loud-stop')
   call expect_loudstop('bad-geom',   'inc_pair_init', 'fragments', '(q) collision-geometry loud-stop')
   call expect_loudstop('bad-noinit', 'beam_draw',   'beam_init', '(r) draw-before-init loud-stop')
   call expect_loudstop('bad-aim2',   'inc_pair_init', 'surface-paradigm', '(w) fixed-site-mode rejection loud-stop')

   print '(a,i0,a)', '=== check_inc_pair: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   subroutine fresh_state(sta)
      type(state_t), intent(inout) :: sta
      sta%q = 0.0d0
      sta%p = 0.0d0
   end subroutine fresh_state

   real(8) function dens(vv)
      real(8), intent(in) :: vv
      dens = vv**3*exp(-((vv - v0)/vwid)**2)
   end function dens

   real(8) function chisq_eq(counts, n)
      integer, intent(in) :: counts(:), n
      real(8) :: ev
      integer :: i
      ev = dble(n)/dble(size(counts))
      chisq_eq = 0.0d0
      do i = 1, size(counts)
         chisq_eq = chisq_eq + (dble(counts(i)) - ev)**2/ev
      end do
   end function chisq_eq

   function base_params() result(p)
      type(beam_params_t) :: p
      p%n_e_rel = 1
      p%t_trans = 0.0d0
      p%v_width = 0.0d0
      p%n_thta = 0
      p%thta_max = 0.7853981633974483d0
      p%n_chi = 0
      p%chi = 0.0d0
      p%n_b = 0
      p%n_aim = 0
   end function base_params

   function cos_mode_params() result(p)
      type(beam_params_t) :: p
      p = base_params()
      p%n_thta = 1
      p%thta_max = pi/2.0d0
      p%n_chi = 0
      p%chi = 0.3d0
      p%n_b = 0
   end function cos_mode_params

   function azi_mode_params() result(p)
      type(beam_params_t) :: p
      p = base_params()
      p%n_thta = 0
      p%thta_max = 0.7853981633974483d0
      p%n_chi = 1
      p%n_b = 0
   end function azi_mode_params

   subroutine setup_gas(w1, w2)
      real(8), intent(in) :: w1, w2
      call list_atoms_clear()
      list_atoms%i_proj = 1
      list_atoms%i_targ = 2
      allocate (list_atoms%mass(2))
      list_atoms%mass = [ w1, w2 ]
      allocate (list_atoms%frag(2))
      list_atoms%frag(1)%nat = 1
      list_atoms%frag(1)%mass = w1
      allocate (list_atoms%frag(1)%list(1), list_atoms%frag(1)%qz(3))
      list_atoms%frag(1)%list = 1
      list_atoms%frag(1)%qz = 0.0d0
      list_atoms%frag(2)%nat = 1
      list_atoms%frag(2)%mass = w2
      allocate (list_atoms%frag(2)%list(1), list_atoms%frag(2)%qz(3))
      list_atoms%frag(2)%list = 2
      list_atoms%frag(2)%qz = 0.0d0
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
      allocate (reactants%dist_scheme(2))
      reactants%dist_scheme = [ 1, 2 ]
   end subroutine setup_gas

   subroutine setup_trimer_gas()
      integer :: i
      call list_atoms_clear()
      list_atoms%i_proj = 1
      list_atoms%i_targ = 2
      allocate (list_atoms%mass(3))
      list_atoms%mass = 1.0d0
      allocate (list_atoms%frag(3))
      do i = 1, 3
         list_atoms%frag(i)%nat = 1
         list_atoms%frag(i)%mass = 1.0d0
         allocate (list_atoms%frag(i)%list(1), list_atoms%frag(i)%qz(3))
         list_atoms%frag(i)%list = i
         list_atoms%frag(i)%qz = 0.0d0
      end do
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
      allocate (reactants%dist_scheme(3))
      reactants%dist_scheme = 1
   end subroutine setup_trimer_gas

   subroutine set_react(e_rel_v, b_max_v, r_sep_v)
      real(8), intent(in) :: e_rel_v, b_max_v, r_sep_v
      reactants%e_rel = e_rel_v
      reactants%b_max = b_max_v
      reactants%r_sep = r_sep_v
   end subroutine set_react

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
      cmd = 'build'//sep//'check_inc_pair.exe '//trim(cmode)//' 2>build'//sep//'inc_pair_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/inc_pair_err.txt', status='old', action='read', iostat=ios)
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

end program check_inc_pair
