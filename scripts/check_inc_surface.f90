!=====================================================================
! check_inc_surface.f90 - the surface-paradigm incident-member criteria
!                          program (criteria defined FIRST in this header;
!                          this program executes them)
! Compile (from src/, linked against the tree's module objects):
!   gfortran -Jbuild/gfortran -c ../scripts/check_inc_surface.f90 -o build/ck_inc_surface.o
!   gfortran build/ck_inc_surface.o <tree objects> -o build/check_inc_surface.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter rule: the member is driven by directly constructed
! beam_params_t and a directly filled list_atoms (surface_model = 1, cell
! fields, i_proj resolved); e_rel/b_max/r_sep through config%reactants.
! All quantities are RECONSTRUCTED FROM THE WRITTEN STATE.
! Scenarios: (S1) single 2-atom projectile fragment (mass 2 amu, i_proj = 1,
!   i_targ = 0) - the case_surface_window form; (S2) projectile fragment +
!   substrate carrier fragment (the 2026-09-28 relaxed statute).
! Criteria:
!   (a) placement + momentum structure (fixed modes, S1): projectile COM
!       height r_sep exactly, lateral offset (b + r_sep*sin(theta))*
!       (cos(chi),sin(chi)), total momentum m1*v*u_inc with the fixed
!       direction, no momentum split (the surface absorbs the recoil)
!   (v) fixed drop site (n_aim = 2, AIM_X/AIM_Y, on-site B_MAX = 0 with
!       N_B = 0): COM lateral = aim + r_sep*sin(theta)*(cos(chi),sin(chi))
!       exactly, b identically zero, on 20000 draws (the fixed-incidence
!       option of the surface paradigm)
!   (w) substrate-carrier relaxation (S2): a two-fragment surface roster -
!       the projectile fragment is placed per the formula, the substrate
!       fragment's q/p stay BITWISE at their pre-sampling fill (the beam
!       assembly never touches a substrate carrier)
!   (j) unit-cell drop (n_aim = 1, box 3x3 skew 60 deg, S1): drop points
!       fill the parallelogram, both generator marginals uniform (pooled
!       chi2, df = 38, critical 53.38)
!   (k) A_LAT commensurability loud-stop (child): a_lat = 2.0 against the
!       3x3 cell -> stderr names inc_surface_init and the commensurability
!   (l) missing A_LAT loud-stop (child): a_lat unset -> stderr names
!       inc_surface_init and the lattice constant
!   (m) unresolved-projectile loud-stop (child): i_proj = 0 -> stderr
!       names inc_surface_init and the projectile role
!   (n) on-site-with-disk loud-stop (child): B_MAX = 0 with N_B = 1 ->
!       stderr names inc_surface_init and N_B = 0
!   (p) missing unit cell loud-stop (child): n_aim = 1 with no cell ->
!       stderr names inc_surface_init and the unit cell
!   (q) fixed-site-without-coordinates loud-stop (child): n_aim = 2 with
!       no AIM_X/AIM_Y -> stderr names inc_surface_init and AIM_X
!   (r) draw-before-init loud-stop (child): the beam_laws guard site
!=====================================================================
program check_inc_surface
   use consts,       only: pi, two_pi, e_conv
   use rng,          only: rng_init
   use config,       only: reactants
   use config_atoms, only: list_atoms, list_atoms_clear
   use state,        only: state_t, state_create, state_destroy
   use beam_laws,    only: beam_params_t
   use inc_surface,  only: inc_surface_init, inc_surface_sample
   implicit none
   type(state_t) :: sta
   type(beam_params_t) :: par
   integer, parameter :: ndraw = 400000, ncheck = 20000, nbin = 20
   integer :: nfail = 0, draw
   integer :: hist(nbin)
   real(8) :: v, costh, sinth, chi_v, chi2, rx0, ry0, s1, s2
   real(8) :: proj_com(3), p1v(3), lat_exp(2)
   logical :: ok_place
   character(len=32) :: mode
   real(8), parameter :: thta_fix = 1.1d0, chi_fix = 0.9d0
   real(8), parameter :: aim_x_fix = 1.3d0, aim_y_fix = 0.7d0

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call setup_surface1()
      call set_react(5.0d0, 3.0d0, 8.0d0)
      par = surf_params()
      select case (trim(mode))
      case ('bad-alat')                   ! (k): incommensurate lattice constant
         list_atoms%a_lat = 2.0d0
         call inc_surface_init(par)
      case ('bad-noalat')                 ! (l): no lattice constant buffered
         list_atoms%a_lat = 0.0d0
         call inc_surface_init(par)
      case ('bad-noproj')                 ! (m): unresolved projectile role
         list_atoms%i_proj = 0
         call inc_surface_init(par)
      case ('bad-b0disk')                 ! (n): on-site beam with the disk mode
         call set_react(5.0d0, 0.0d0, 8.0d0)
         par%n_b = 1
         call inc_surface_init(par)
      case ('bad-cell')                   ! (p): n_aim = 1 with no cell buffered
         list_atoms%box_lx = 0.0d0
         list_atoms%box_ly = 0.0d0
         par%n_aim = 1
         call inc_surface_init(par)
      case ('bad-aim-nocoord')            ! (q): fixed-site mode without coordinates
         par%n_aim = 2
         call inc_surface_init(par)
      case ('bad-noinit')                 ! (r): sample without init
         call state_create(sta, 2)
         sta%q = 0.0d0
         sta%p = 0.0d0
         call inc_surface_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- (a) placement + momentum structure (single-fragment form) ----------
   call setup_surface1()
   call set_react(5.0d0, 3.0d0, 8.0d0)
   call state_create(sta, 2)
   par = surf_params()
   par%n_thta = 0
   par%thta_max = thta_fix
   par%n_chi = 0
   par%chi = chi_fix
   par%n_b = 0
   call inc_surface_init(par)
   ok_place = .true.
   do draw = 1, ncheck
      call fresh_state(sta)
      call inc_surface_sample(sta)
      proj_com = (sta%q(1:3) + sta%q(4:6))/2.0d0
      p1v = sta%p(1:3) + sta%p(4:6)
      v = sqrt(dot_product(p1v, p1v))
      costh = -p1v(3)/v
      chi_v = atan2(-p1v(2)/v, -p1v(1)/v)
      if (chi_v < 0.0d0) chi_v = chi_v + two_pi
      sinth = sqrt(1.0d0 - costh*costh)
      if (abs(proj_com(3) - 8.0d0) > 1.0d-12 .or. &
          norm2(proj_com(1:2) - (3.0d0 + 8.0d0*sinth)*[ cos(chi_fix), sin(chi_fix) ]) > 1.0d-12 .or. &
          abs(v - 2.0d0*sqrt(2.0d0*e_conv*5.0d0/2.0d0)) > 1.0d-9*v .or. &
          abs(costh - cos(thta_fix)) > 1.0d-12 .or. &
          abs(chi_v - chi_fix) > 1.0d-12) ok_place = .false.
   end do
   if (.not. ok_place) then
      print '(a)', 'FAIL (a) placement/momentum structure: surface geometry wrong'
      nfail = nfail + 1
   else
      print '(a)', '(a) placement + momentum structure: PASS (COM height r_sep, lateral '// &
         'formula, total p = m1*v*u with the fixed direction - the surface absorbs the recoil)'
   end if
   call state_destroy(sta)

   ! ---------- (v) fixed drop site + on-site beam ----------
   call set_react(5.0d0, 0.0d0, 8.0d0)
   call state_create(sta, 2)
   par = surf_params()
   par%n_thta = 0
   par%thta_max = thta_fix
   par%n_chi = 0
   par%chi = chi_fix
   par%n_b = 0
   par%n_aim = 2
   par%aim_x = aim_x_fix
   par%aim_y = aim_y_fix
   par%aim_x_given = .true.
   par%aim_y_given = .true.
   call inc_surface_init(par)
   ok_place = .true.
   do draw = 1, ncheck
      call fresh_state(sta)
      call inc_surface_sample(sta)
      proj_com = (sta%q(1:3) + sta%q(4:6))/2.0d0
      p1v = sta%p(1:3) + sta%p(4:6)
      v = sqrt(dot_product(p1v, p1v))
      costh = -p1v(3)/v
      sinth = sqrt(1.0d0 - costh*costh)
      lat_exp = [ aim_x_fix, aim_y_fix ] + 8.0d0*sinth*[ cos(chi_fix), sin(chi_fix) ]
      if (abs(proj_com(3) - 8.0d0) > 1.0d-12 .or. &
          norm2(proj_com(1:2) - lat_exp) > 1.0d-12 .or. &
          abs(v - 2.0d0*sqrt(2.0d0*e_conv*5.0d0/2.0d0)) > 1.0d-9*v) ok_place = .false.
   end do
   if (.not. ok_place) then
      print '(a)', 'FAIL (v) fixed drop site: the placement does not follow AIM_X/AIM_Y'
      nfail = nfail + 1
   else
      print '(a)', '(v) fixed drop site: PASS (n_aim = 2 + AIM_X/AIM_Y with the on-site '// &
         'B_MAX = 0 beam: COM = aim + r_sep*sin(theta)*(cos,si) exactly, no disk offset)'
   end if
   call state_destroy(sta)

   ! ---------- (w) substrate-carrier relaxation (two-fragment surface roster) ----------
   call setup_surface2()
   call set_react(5.0d0, 3.0d0, 8.0d0)
   call state_create(sta, 3)
   par = surf_params()
   par%n_thta = 0
   par%thta_max = thta_fix
   par%n_chi = 0
   par%chi = chi_fix
   par%n_b = 0
   call inc_surface_init(par)
   ok_place = .true.
   do draw = 1, ncheck
      ! the substrate atoms carry a NONZERO sentinel fill: the beam assembly
      ! must leave them bitwise untouched
      sta%q = 0.0d0
      sta%p = 0.0d0
      sta%q(4:9) = [ 0.11d0, 0.22d0, 0.33d0, 0.44d0, 0.55d0, 0.66d0 ]
      sta%p(4:9) = [ 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0 ]
      call inc_surface_sample(sta)
      if (any(sta%q(4:9) /= [ 0.11d0, 0.22d0, 0.33d0, 0.44d0, 0.55d0, 0.66d0 ]) .or. &
          any(sta%p(4:9) /= [ 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0 ])) ok_place = .false.
      ! projectile placement: monoatomic fragment 1 at the site with full beam velocity
      v = sqrt(dot_product(sta%p(1:3), sta%p(1:3)))/12.0d0
      costh = -sta%p(3)/(12.0d0*v)
      sinth = sqrt(1.0d0 - costh*costh)
      lat_exp = 3.0d0*[ cos(chi_fix), sin(chi_fix) ] + 8.0d0*sinth*[ cos(chi_fix), sin(chi_fix) ]
      if (norm2(sta%q(1:2) - lat_exp) > 1.0d-12 .or. abs(sta%q(3) - 8.0d0) > 1.0d-12 .or. &
          abs(v - sqrt(2.0d0*e_conv*5.0d0/12.0d0)) > 1.0d-9*v) ok_place = .false.
   end do
   if (.not. ok_place) then
      print '(a)', 'FAIL (w) substrate relaxation: the substrate carriers moved or the projectile is misplaced'
      nfail = nfail + 1
   else
      print '(a)', '(w) substrate-carrier relaxation: PASS (the two-fragment surface roster: '// &
         'projectile placed with the full beam velocity, substrate q/p bitwise untouched)'
   end if
   call state_destroy(sta)

   ! ---------- (j) unit-cell drop (n_aim = 1) ----------
   call setup_surface1()
   call set_react(5.0d0, 3.0d0, 8.0d0)
   call state_create(sta, 2)
   list_atoms%box_lx = 3.0d0
   list_atoms%box_ly = 3.0d0
   list_atoms%skew = pi/3.0d0
   par = surf_params()
   par%n_aim = 1
   call inc_surface_init(par)
   call rng_init(7319)
   hist = 0
   do draw = 1, ndraw
      call fresh_state(sta)
      call inc_surface_sample(sta)
      proj_com = (sta%q(1:3) + sta%q(4:6))/2.0d0
      ! strip the FIXED beam offsets (n_thta = 0 at pi/4, n_chi = 0, n_b = 0
      ! at b_max) to recover the bare drop point
      rx0 = proj_com(1) - (3.0d0 + 8.0d0*sin(0.7853981633974483d0))*cos(0.0d0)
      ry0 = proj_com(2) - (3.0d0 + 8.0d0*sin(0.7853981633974483d0))*sin(0.0d0)
      s2 = ry0/(3.0d0*sin(pi/3.0d0))
      s1 = (rx0 - s2*3.0d0*cos(pi/3.0d0))/3.0d0
      if (s1 < -1.0d-12 .or. s1 > 1.0d0 + 1.0d-12 .or. s2 < -1.0d-12 .or. s2 > 1.0d0 + 1.0d-12) &
         ok_place = .false.
      hist(min(nbin, max(1, int(nbin*s1) + 1))) = hist(min(nbin, max(1, int(nbin*s1) + 1))) + 1
      hist(min(nbin, max(1, int(nbin*s2) + 1))) = hist(min(nbin, max(1, int(nbin*s2) + 1))) + 1
   end do
   chi2 = chisq_eq(hist, 2*ndraw)
   if (chi2 > 53.38d0) then
      print '(a,f9.3,a)', 'FAIL (j) unit-cell drop: chi2 = ', chi2, ' > 53.38'
      nfail = nfail + 1
   else
      print '(a,f7.3,a)', '(j) unit-cell drop: PASS (chi2 = ', chi2, &
         ' <= 53.38, df = 38, alpha = 0.05, both generator marginals pooled)'
   end if
   call state_destroy(sta)

   ! ---------- (k)-(r) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-alat',        'inc_surface_init', 'commensurate', '(k) A_LAT-commensurability loud-stop')
   call expect_loudstop('bad-noalat',      'inc_surface_init', 'lattice constant', '(l) missing-A_LAT loud-stop')
   call expect_loudstop('bad-noproj',      'inc_surface_init', 'projectile role', '(m) unresolved-projectile loud-stop')
   call expect_loudstop('bad-b0disk',      'inc_surface_init', 'N_B = 0', '(n) on-site-with-disk loud-stop')
   call expect_loudstop('bad-cell',        'inc_surface_init', 'unit cell', '(p) missing-unit-cell loud-stop')
   call expect_loudstop('bad-aim-nocoord', 'inc_surface_init', 'AIM_X', '(q) fixed-site-no-coordinates loud-stop')
   call expect_loudstop('bad-noinit',      'beam_draw', 'beam_init', '(r) draw-before-init loud-stop')

   print '(a,i0,a)', '=== check_inc_surface: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   subroutine fresh_state(sta)
      type(state_t), intent(inout) :: sta
      sta%q = 0.0d0
      sta%p = 0.0d0
   end subroutine fresh_state

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

   function surf_params() result(p)
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
      p%aim_x = 0.0d0
      p%aim_y = 0.0d0
      p%aim_x_given = .false.
      p%aim_y_given = .false.
   end function surf_params

   ! S1: single 2-atom projectile fragment (the case_surface_window form)
   subroutine setup_surface1()
      call list_atoms_clear()
      list_atoms%i_proj = 1
      list_atoms%i_targ = 0
      allocate (list_atoms%mass(2))
      list_atoms%mass = 1.0d0
      allocate (list_atoms%frag(1))
      list_atoms%frag(1)%nat = 2
      list_atoms%frag(1)%mass = 2.0d0
      allocate (list_atoms%frag(1)%list(2), list_atoms%frag(1)%qz(6))
      list_atoms%frag(1)%list = [ 1, 2 ]
      list_atoms%frag(1)%qz = 0.0d0
      list_atoms%box_lx = 3.0d0
      list_atoms%box_ly = 3.0d0
      list_atoms%skew = pi/2.0d0
      list_atoms%a_lat = 3.0d0
      reactants%surface_model = 1
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
      allocate (reactants%dist_scheme(1))
      reactants%dist_scheme = 1
   end subroutine setup_surface1

   ! S2: monoatomic projectile fragment (mass 12) + two-atom substrate carrier
   ! (the 2026-09-28 relaxed statute: every non-projectile fragment is substrate)
   subroutine setup_surface2()
      call list_atoms_clear()
      list_atoms%i_proj = 1
      list_atoms%i_targ = 0
      allocate (list_atoms%mass(3))
      list_atoms%mass = [ 12.0d0, 196.97d0, 196.97d0 ]
      allocate (list_atoms%frag(2))
      list_atoms%frag(1)%nat = 1
      list_atoms%frag(1)%mass = 12.0d0
      allocate (list_atoms%frag(1)%list(1), list_atoms%frag(1)%qz(3))
      list_atoms%frag(1)%list = 1
      list_atoms%frag(1)%qz = 0.0d0
      list_atoms%frag(2)%nat = 2
      list_atoms%frag(2)%mass = 2.0d0*196.97d0
      allocate (list_atoms%frag(2)%list(2), list_atoms%frag(2)%qz(6))
      list_atoms%frag(2)%list = [ 2, 3 ]
      list_atoms%frag(2)%qz = 0.0d0
      list_atoms%box_lx = 3.0d0
      list_atoms%box_ly = 3.0d0
      list_atoms%skew = pi/2.0d0
      list_atoms%a_lat = 3.0d0
      reactants%surface_model = 1
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
      allocate (reactants%dist_scheme(2))
      reactants%dist_scheme = [ 1, 5 ]   ! boltzmann on the projectile, surface_oscillator on the substrate
   end subroutine setup_surface2

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
      cmd = 'build'//sep//'check_inc_surface.exe '//trim(cmode)//' 2>build'//sep//'inc_surf_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/inc_surf_err.txt', status='old', action='read', iostat=ios)
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

end program check_inc_surface
