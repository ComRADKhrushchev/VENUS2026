!=====================================================================
! check_dist_rotation.f90 - rotation distribution-member criteria program
!                            (criteria defined FIRST in this header; this
!                            program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/config.f90 utils/specio.f90 &
!            facade/sysdef_facade.f90 globals/config_roster.f90 globals/state.f90 &
!            utils/rng.f90 utils/linalg.f90 facade/sampler.f90 methods/dist_rotation.f90 &
!            ../scripts/check_dist_rotation.f90 -o build/check_dist_rotation.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter statute: the member is driven by directly
! constructed parameters and a directly filled roster. The anchors below
! are exact algebraic identities of the rigid-rotor construction
! (L = sum r x p, I at the WRITTEN geometry, COM frame), evaluated in
! this program independently of the member's internal path.
! Loud-stop verification approach (same scheme as the sibling checks):
! every loud-stop criterion runs as a CHILD PROCESS - the check re-spawns
! itself with a mode argument (execute_command_line); the child drives
! the REAL member entry into the stop and must print a naming message on
! stderr and exit nonzero. The parent asserts the exit status and the
! stderr content.
! Criteria:
!   (a) angular-momentum magnitude exact: for J = 0, 1, 2, 3, 4, 5 over a
!       LINEAR carrier (masses 1/2 amu at +/-0.6 Ang), 500 draws each:
!       |sum r x p| = sqrt(J(J+1))*hbar_code to 1e-10 relative (J = 0:
!       identically zero momenta)
!   (b) zero total momentum exact: |sum p| <= 1e-9 hbar_code-scale for
!       every draw of (a) (the construction leaves no residual COM drift)
!   (c) linear-rotor energy identity: E = sum |p|^2/(2m) equals
!       |L|^2/(2*I_perp) with I_perp computed here at the WRITTEN
!       geometry, to 1e-10 relative (E in internal units; a wrong
!       inertia or omega construction breaks the identity)
!   (d) L perpendicular to the bond: |L.e_bond| <= 1e-9*|L| for the
!       linear carrier (rotation about the molecular axis is no degree
!       of freedom of the linear rotor)
!   (e) nonlinear carrier (3 atoms, masses 1/2/3, non-symmetric
!       geometry), J = 2, 500 draws: |L| exact (1e-10) AND the rigid-body
!       energy identity E = L^T I^-1 L / 2 with the inertia tensor
!       computed here at the written geometry (1e-10 relative; the
!       principal-axis decomposition of omega is what this pins)
!   (f) rigid-geometry preservation: every interatomic distance of the
!       written fragment matches the equilibrium geometry to 1e-10
!       Angstrom, over all draws of (e)
!   (g) orientation isotropy of the bond axis (linear carrier, J = 1,
!       20000 draws): the squared norm of the mean axis unit vector is
!       tested chi-square(3)-calibrated - 3N*|mean axis|^2 <= 16.27 (the
!       99.9% upper critical of chi^2_3; per-component sigma is
!       sqrt(1/(3N)) = 4.08e-3, and a 3-vector norm needs the joint, not
!       the per-component, distribution), and the mean of the squared z
!       component within 5 sigma of 1/3 (sigma_z2 = sqrt((1/5 - 1/9)/N) =
!       sqrt(4/(45N)) = 6.67e-3)
!   (h) negative J loud-stop (child): j_rot = -1 -> stderr names
!       rotation_init and j_rot
!   (i) single-atom carrier loud-stop (child): nat = 1 -> stderr names
!       rotation_init and the atom count
!   (j) degenerate geometry loud-stop (child): coincident atoms ->
!       stderr names rotation_init and the moment of inertia
!   (k) no-fragment-binding loud-stop (child): dist_scheme carries no
!       rotation row -> stderr names rotation_init
!   (l) duplicate-fragment-binding loud-stop (child): two carriers ->
!       stderr names rotation_init
!   (m) sample-before-init loud-stop (child): rotation_sample with no
!       rotation_init -> stderr names rotation_sample
!=====================================================================
program check_dist_rotation
   use consts,        only: hbar_code, e_conv
   use rng,           only: rng_init
   use config,        only: reactants
   use config_roster, only: roster, roster_clear
   use state,         only: state_t, state_create, state_destroy
   use dist_rotation, only: rotation_params_t, rotation_init, rotation_sample
   use linalg,        only: eig_sym
   implicit none
   type(state_t) :: sta
   type(rotation_params_t) :: par
   integer, parameter :: ndraw_id = 500
   integer, parameter :: ndraw_iso = 20000
   integer :: nfail = 0, nfail0, i, k, draw, jval
   integer :: nat, at(3), list(3)
   real(8) :: l_vec(3), l_mag, e_kin, e_id, i_perp, rel
   real(8) :: r(3), tens(3,3), ev(3), vecs(3,3), tin(3,3)
   real(8) :: com(3), wt, axis(3), axlen, asum(3), az2
   real(8) :: d_eq(3,3), d_new, dev_max
   real(8) :: chi_iso
   real(8) :: qz_lin(6), qz_non(9), mass_non(3)
   character(len=32) :: mode
   logical :: ok_all

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 3)
      par%j_rot = 2
      select case (trim(mode))
      case ('bad-jneg')                    ! (h): j_rot = -1
         call setup_lin()
         call set_scheme([3])
         par%j_rot = -1
         call rotation_init(par)
      case ('bad-atom')                    ! (i): single-atom carrier
         call setup_atom()
         call set_scheme([3])
         call rotation_init(par)
      case ('bad-degen')                   ! (j): coincident atoms, I = 0
         call setup_degen()
         call set_scheme([3])
         call rotation_init(par)
      case ('bad-bind0')                   ! (k): no rotation row
         call setup_lin()
         call set_scheme([1])
         call rotation_init(par)
      case ('bad-bind2')                   ! (l): two carriers
         call setup_two()
         call set_scheme([3, 3])
         call rotation_init(par)
      case ('bad-noinit')                  ! (m): sample without init
         call setup_lin()
         call set_scheme([3])
         call rotation_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- (a)-(d) linear carrier identities ----------
   qz_lin = lin_geom()
   do jval = 0, 5                           ! J = 0..5 (header (a) set)
      call setup_lin()
      call set_scheme([3])
      call state_create(sta, 2)
      par%j_rot = jval
      call rng_init(104729 + jval)
      call rotation_init(par)
      l_mag = sqrt(dble(jval*(jval + 1)))*hbar_code
      ok_all = .true.
      do draw = 1, ndraw_id
         call rotation_sample(sta)
         call frag_lxpr(sta, 2, l_vec, e_kin)
         if (abs(norm2(l_vec) - l_mag) > 1.0d-10*max(l_mag, 1.0d-12)) ok_all = .false.  ! (a)
         if (jval == 0 .and. any(sta%p(1:6) /= 0.0d0)) ok_all = .false.                ! (a) J=0
         if (norm2(sta%p(1:3) + sta%p(4:6)) > 1.0d-9*hbar_code) ok_all = .false.  ! (b)
         ! (c) linear-rotor energy identity at the written geometry
         i_perp = 0.0d0
         com = 0.0d0
         wt = sum(roster%mass(1:2))
         do i = 1, 2
            com = com + roster%mass(i)*sta%q(3*i-2:3*i)
         end do
         com = com/wt
         do i = 1, 2
            r = sta%q(3*i-2:3*i) - com
            i_perp = i_perp + roster%mass(i)*dot_product(r, r)
         end do
         e_id = l_mag*l_mag/(2.0d0*i_perp)
         if (abs(e_kin - e_id) > 1.0d-10*e_id) ok_all = .false.
         ! (d) L perpendicular to the bond axis
         axis = sta%q(4:6) - sta%q(1:3)
         axlen = norm2(axis)
         if (abs(dot_product(l_vec, axis/axlen)) > 1.0d-9*max(l_mag, 1.0d-12)) ok_all = .false.
      end do
      nfail0 = nfail
      if (ok_all) then
         print '(a,i0,a)', '(a)-(d) J = ', jval, ': PASS (|L|, zero P, E=|L|^2/(2I), L perp bond; 500 draws)'
      else
         print '(a,i0,a)', 'FAIL (a)-(d) J = ', jval, ': an identity broke (|L| / total P / energy / perpendicularity)'
         nfail = nfail + 1
      end if
      call state_destroy(sta)
   end do

   ! ---------- (e)/(f) nonlinear carrier: |L| + rigid-body energy identity ----------
   qz_non = non_geom()
   mass_non = [ 1.0d0, 2.0d0, 3.0d0 ]
   call setup_non()
   call set_scheme([3])
   call state_create(sta, 3)
   par%j_rot = 2
   call rng_init(15485863)
   call rotation_init(par)
   l_mag = sqrt(6.0d0)*hbar_code
   ! equilibrium pair distances (rigid-geometry reference)
   do i = 1, 3
      do k = i + 1, 3
         d_eq(i,k) = norm2(qz_non(3*k-2:3*k) - qz_non(3*i-2:3*i))
      end do
   end do
   ok_all = .true.
   dev_max = 0.0d0
   do draw = 1, ndraw_id
      call rotation_sample(sta)
      call frag_lxpr(sta, 3, l_vec, e_kin)
      if (abs(norm2(l_vec) - l_mag) > 1.0d-10*l_mag) ok_all = .false.
      ! inertia tensor at the written geometry, I^-1 L, energy identity
      nat = 3
      at = [ 1, 2, 3 ]
      com = 0.0d0
      wt = sum(mass_non)
      do i = 1, nat
         com = com + mass_non(i)*sta%q(3*at(i)-2:3*at(i))
      end do
      com = com/wt
      tens = 0.0d0
      do i = 1, nat
         r = sta%q(3*at(i)-2:3*at(i)) - com
         tens(1,1) = tens(1,1) + mass_non(i)*(r(2)**2 + r(3)**2)
         tens(2,2) = tens(2,2) + mass_non(i)*(r(1)**2 + r(3)**2)
         tens(3,3) = tens(3,3) + mass_non(i)*(r(1)**2 + r(2)**2)
         tens(1,2) = tens(1,2) - mass_non(i)*r(1)*r(2)
         tens(1,3) = tens(1,3) - mass_non(i)*r(1)*r(3)
         tens(2,3) = tens(2,3) - mass_non(i)*r(2)*r(3)
      end do
      tens(2,1) = tens(1,2)
      tens(3,1) = tens(1,3)
      tens(3,2) = tens(2,3)
      tin = tens
      call eig_sym(tin, ev)                 ! tin now eigenvectors by column
      vecs = tin
      e_id = 0.0d0
      do k = 1, 3
         e_id = e_id + dot_product(l_vec, vecs(:,k))**2/ev(k)
      end do
      e_id = e_id/2.0d0
      if (abs(e_kin - e_id) > 1.0d-10*e_id) ok_all = .false.
      ! (f) rigid geometry
      do i = 1, 3
         do k = i + 1, 3
            d_new = norm2(sta%q(3*at(k)-2:3*at(k)) - sta%q(3*at(i)-2:3*at(i)))
            dev_max = max(dev_max, abs(d_new - d_eq(i,k)))
         end do
      end do
   end do
   if (dev_max > 1.0d-10) ok_all = .false.
   nfail0 = nfail
   if (ok_all) then
      print '(a,f8.2,a)', '(e)+(f) nonlinear carrier J = 2: PASS (|L| exact, E = L^T I^-1 L/2, geometry rigid; max dev = ', &
         dev_max, ' Ang)'
   else
      print '(a,f8.2,a)', 'FAIL (e)/(f) nonlinear carrier: identity or rigidity broke; max dev = ', dev_max
      nfail = nfail + 1
   end if
   call state_destroy(sta)

   ! ---------- (g) orientation isotropy of the bond axis ----------
   call setup_lin()
   call set_scheme([3])
   call state_create(sta, 2)
   par%j_rot = 1
   call rng_init(32452843)
   call rotation_init(par)
   asum = 0.0d0
   az2 = 0.0d0
   do draw = 1, ndraw_iso
      call rotation_sample(sta)
      axis = sta%q(4:6) - sta%q(1:3)
      axis = axis/norm2(axis)
      asum = asum + axis
      az2 = az2 + axis(3)**2
   end do
   asum = asum/dble(ndraw_iso)
   az2 = az2/dble(ndraw_iso)
   nfail0 = nfail
   chi_iso = 3.0d0*dble(ndraw_iso)*dot_product(asum, asum)
   if (chi_iso > 16.27d0 .or. &
       abs(az2 - 1.0d0/3.0d0) > 5.0d0*sqrt(4.0d0/(45.0d0*dble(ndraw_iso)))) then
      print '(a,3f8.4,a,f8.4)', 'FAIL (g) bond-axis isotropy: mean ', asum, ' mean z^2 ', az2
      nfail = nfail + 1
   else
      print '(a,f6.2,a,3f8.4,a,f8.4)', '(g) bond-axis isotropy: PASS (chi^2_3 = ', chi_iso, &
         ' <= 16.27; mean axis ', asum, ', mean z^2 = ', az2, ' vs 1/3; N = 20000)'
   end if
   call state_destroy(sta)

   ! ---------- (h)-(m) loud-stop criteria (child processes) ----------
   call expect_loudstop('bad-jneg',    'rotation_init',   'j_rot',     '(h) negative J loud-stop')
   call expect_loudstop('bad-atom',    'rotation_init',   'atom',      '(i) single-atom carrier loud-stop')
   call expect_loudstop('bad-degen',   'rotation_init',   'inertia',   '(j) degenerate-geometry loud-stop')
   call expect_loudstop('bad-bind0',   'rotation_init',   'rotation scheme', '(k) no-fragment-binding loud-stop')
   call expect_loudstop('bad-bind2',   'rotation_init',   'two',       '(l) duplicate-fragment-binding loud-stop')
   call expect_loudstop('bad-noinit',  'rotation_sample', 'rotation_init', '(m) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_dist_rotation: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! fragment angular momentum and kinetic energy (COM-frame L = sum r x p
   ! with r about the fragment COM; E internal = sum |p|^2/(2m))
   subroutine frag_lxpr(st, n_at, lv, ek)
      type(state_t), intent(in) :: st
      integer, intent(in) :: n_at
      real(8), intent(out) :: lv(3), ek
      real(8) :: cm(3), rr(3), wtt
      integer :: m
      cm = 0.0d0
      wtt = sum(roster%mass(1:n_at))
      do m = 1, n_at
         cm = cm + roster%mass(m)*st%q(3*m-2:3*m)
      end do
      cm = cm/wtt
      lv = 0.0d0
      ek = 0.0d0
      do m = 1, n_at
         rr = st%q(3*m-2:3*m) - cm
         lv = lv + cross(rr, st%p(3*m-2:3*m))
         ek = ek + dot_product(st%p(3*m-2:3*m), st%p(3*m-2:3*m))/(2.0d0*roster%mass(m))
      end do
   end subroutine frag_lxpr

   pure function cross(a, b) result(c)
      real(8), intent(in) :: a(3), b(3)
      real(8) :: c(3)
      c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
   end function cross

   ! linear two-atom geometry: masses 1/2 amu at +/-0.6 Ang on x
   function lin_geom() result(g)
      real(8) :: g(6)
      g = 0.0d0
      g(1) = -0.6d0
      g(4) = 0.6d0
   end function lin_geom

   ! bent three-atom geometry (non-principal, non-symmetric)
   function non_geom() result(g)
      real(8) :: g(9)
      g = 0.0d0
      g(1) = -0.7d0
      g(4) = 0.5d0
      g(8) = 0.45d0
      g(9) = 0.2d0
   end function non_geom

   subroutine setup_lin()
      call fill_roster(2, lin_geom(), [ 1.0d0, 2.0d0 ])
   end subroutine setup_lin

   subroutine setup_non()
      call fill_roster(3, non_geom(), [ 1.0d0, 2.0d0, 3.0d0 ])
   end subroutine setup_non

   subroutine setup_atom()
      real(8) :: g(3)
      g = 0.0d0
      call fill_roster(1, g, [ 1.0d0 ])
   end subroutine setup_atom

   subroutine setup_degen()
      real(8) :: g(6)
      g = 0.0d0
      call fill_roster(2, g, [ 1.0d0, 2.0d0 ])
   end subroutine setup_degen

   subroutine setup_two()
      call roster_clear()
      allocate (roster%mass(4))
      roster%mass = [ 1.0d0, 2.0d0, 1.0d0, 2.0d0 ]
      allocate (roster%frag(2))
      roster%frag(1)%nat = 2
      roster%frag(2)%nat = 2
      roster%frag(1)%mass = 3.0d0
      roster%frag(2)%mass = 3.0d0
      allocate (roster%frag(1)%list(2), roster%frag(1)%qz(6))
      roster%frag(1)%list = [ 1, 2 ]
      roster%frag(1)%qz = lin_geom()
      allocate (roster%frag(2)%list(2), roster%frag(2)%qz(6))
      roster%frag(2)%list = [ 3, 4 ]
      roster%frag(2)%qz = lin_geom()
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine setup_two

   ! mini-roster by direct field fill: one fragment of nat atoms
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
      cmd = 'build'//sep//'check_dist_rotation.exe '//trim(cmode)//' 2>build'//sep//'rotation_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/rotation_err.txt', status='old', action='read', iostat=ios)
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

end program check_dist_rotation
