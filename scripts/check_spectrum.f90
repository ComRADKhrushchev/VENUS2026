!=====================================================================
! check_spectrum.f90 - the member-seam spectrum criterion program
!                       (criteria defined FIRST in this header; this
!                       program executes them)
! Criteria (2026-10-05 spectrum-source legislation: SPECTRUM_SOURCE
!   selects ONE producer - the internal derivation (COMPUTE, default)
!   or the MANUAL spectrum-table file; the container export channel is
!   withdrawn, docs/plans/2026-10-05-spectrum-source.md):
!   (a) the unit-gradient calibration: a container delivering
!       g = 1 eV/A measures the interface fold through the probe entry
!       (the eV conversion literal never lands on the check side - the
!       measured fold feeds every expected value below)
!   (b) a synthetic Morse diatomic (D = 5.0 eV, beta = 1.8 1/A,
!       re = 0.9 A): the COMPUTE pull reproduces omega_analytic =
!       sqrt(2*D*beta^2*fold/mu) to 2e-7 relative (mu = m_H/2; the
!       residual is the Hessian central-difference truncation)
!   (c) the BKMP2 H2 pair limit at the buffer probe geometry
!       (spectator at atom-scale isolation): the pull reproduces the
!       frequency from the INDEPENDENT energy second difference
!       V'' = 35.933611286509404 eV/A^2 (recorded against the
!       archived surface; spectator-position independent) to 2e-7
!       relative; nu = 4403.2 cm^-1, against the literature H2 w_e =
!       4401.21 cm^-1 at 0.045% (the Schwenke harmonic curvature -
!       documented, not asserted); the H fragment pulls an
!       allocated-EMPTY table (the monoatomic rule)
!   (d) a synthetic tetrahedral spring molecule (4 H atoms, 6 equal
!       springs k = 30 eV/A^2): 6 modes with the exact eigenvalue
!       pattern {1,1,2,2,2,4}*k/m (2e-7 relative) - the nonlinear
!       3N-6 split, the degeneracy structure, ascending order
!   (e) a non-decaying inter-fragment coupling: the COMPUTE pull
!       named-aborts (child) on the decayed-tail phrase, pointing at
!       SPECTRUM_SOURCE=MANUAL
!   (f) PARADIGM NEUTRALITY (the old surface blanket ban is
!       withdrawn): surface_model = 1 with the same decaying Morse
!       pair DERIVES through the unchanged COMPUTE path - the pull
!       returns the leg-(b) frequency (a surface label carries no
!       vote; only the PES's own decay behavior does)
!   (g) the MANUAL producer, end to end: a spectrum-table file staged
!       with stretch_mode's own column and constants (w_e distinct
!       from the WVN row - the two channels must not cross) round
!       trips exactly: the mode table reproduces the written WVN, the
!       column reproduces stretch_mode bit for bit, the DIAT record
!       reproduces (w_e, w_e*x_e, B_e); a duplicate composition in the
!       file named-aborts (child); a carrier with no matching record
!       named-aborts (child); MANUAL without SPECTRUM_FILE
!       named-aborts (child); the EBK diatomic pull under COMPUTE
!       named-aborts with the no-derivation phrase (child)
!   (h) a well-pinned "molecule" (no zero modes) named-aborts on the
!       ambiguous split (child); stretch_mode is unit-norm with the
!       antisymmetric mass-weighted structure; the wavenumber
!       conversion round-trips (parent, pure)
! Compile (from src/; modules in use-order first, -Jbuild):
!   gfortran -c -Jbuild -o build/bkmp2_pot.o ../systems/inactive/h3_bkmp2/bkmp2_pot.f90
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 utils/rng.f90 &
!            globals/state.f90 globals/config_atoms.f90 &
!            interface/elec_interface.f90 methods/elec_adiabatic.f90 &
!            utils/linalg.f90 utils/hessian.f90 utils/spectrum.f90 &
!            interface/force_interface.f90 interface/spectrum_interface.f90 &
!            build/bkmp2_pot.o ../scripts/check_spectrum.f90 -o build/check_spectrum.exe
!   ./build/check_spectrum.exe
!=====================================================================
module check_containers
   ! the synthetic containers of the criterion legs (adapters over the
   ! force-slot contract; eV/A delivery, exactly like a real container)
   use bkmp2_pot, only: bkmp2_vg
   implicit none
   private
   public :: unit_pes, morse_pes, bkmp2_adapter, tetra_pes, escal_pes, well_pes
   real(8), parameter :: mh = 1.008d0                            ! model H mass [amu]
   real(8), parameter :: d_mo = 5.0d0, b_mo = 1.8d0, re_mo = 0.9d0 ! the Morse set [eV/1-A/A]
   real(8), parameter :: k_tet = 30.0d0                          ! spring [eV/A^2]
   real(8), parameter :: k_well = 5.0d0                          ! well pinning [eV/A^2]
   public :: mh, d_mo, b_mo, re_mo, k_tet
contains
   ! unit_pes - g(1) = 1 eV/A, everything else zero (the fold calibration)
   subroutine unit_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      g = 0.0d0
      g(1) = 1.0d0
      v = 0.0d0
   end subroutine unit_pes

   ! morse_pes - one H2 Morse pair (equilibrium at the buffered bond)
   subroutine morse_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      real(8) :: r, e1, dvdr, u(3)
      r = norm2(q(4:6) - q(1:3))
      e1 = exp(-b_mo*(r - re_mo))
      v = d_mo*(e1*e1 - 2.0d0*e1)
      dvdr = 2.0d0*d_mo*b_mo*(e1 - e1*e1)
      u = (q(4:6) - q(1:3))/r
      g(1:3) = -dvdr*u
      g(4:6) = dvdr*u
   end subroutine morse_pes

   ! bkmp2_adapter - the archived BKMP2 surface behind the slot contract
   subroutine bkmp2_adapter(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      call bkmp2_vg(3, q, v, g)
   end subroutine bkmp2_adapter

   ! tetra_pes - a 4-H tetrahedron with 6 equal edge springs
   subroutine tetra_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      integer, parameter :: pr(2, 6) = reshape([1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, 6])
      real(8) :: u(3), dr, a0, r
      integer :: p, i, j
      a0 = 2.0d0*sqrt(2.0d0)                    ! the unit-cube tetra edge [A]
      g = 0.0d0
      v = 0.0d0
      do p = 1, 6
         i = pr(1, p)
         j = pr(2, p)
         u = q(3*j-2:3*j) - q(3*i-2:3*i)
         r = norm2(u)
         u = u/r
         dr = r - a0
         v = v + 0.5d0*k_tet*dr*dr
         g(3*i-2:3*i) = g(3*i-2:3*i) - k_tet*dr*u
         g(3*j-2:3*j) = g(3*j-2:3*j) + k_tet*dr*u
      end do
   end subroutine tetra_pes

   ! escal_pes - two H2 Morse pairs + a NON-decaying inter-fragment
   ! coupling of the two bond coordinates (constant with separation)
   subroutine escal_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      real(8), parameter :: kc = 5.0d0          ! the undecayed coupling [eV/A^2]
      real(8) :: dc
      g = 0.0d0
      v = 0.0d0
      call add_pair(q, g, v, 1)                 ! atoms 1,2
      call add_pair(q, g, v, 2)                 ! atoms 3,4
      dc = (q(10) - q(7)) - (q(4) - q(1))       ! bondB_x - bondA_x
      v = v + 0.5d0*kc*dc*dc
      g(1) = g(1) + kc*dc
      g(4) = g(4) - kc*dc
      g(7) = g(7) - kc*dc
      g(10) = g(10) + kc*dc
   end subroutine escal_pes

   ! well_pes - per-atom Cartesian wells (NO zero modes - the ambiguous-
   ! split abort leg)
   subroutine well_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      integer :: i
      g = 0.0d0
      v = 0.0d0
      do i = 1, size(q)/3
         v = v + 0.5d0*k_well*dot_product(q(3*i-2:3*i), q(3*i-2:3*i))
         g(3*i-2:3*i) = k_well*q(3*i-2:3*i)
      end do
   end subroutine well_pes

   ! add_pair(q, g, v, i) - one Morse pair over atoms (2i-1, 2i)
   subroutine add_pair(q, g, v, i)
      real(8), intent(in)    :: q(:)
      real(8), intent(inout) :: g(:), v
      integer, intent(in)    :: i
      real(8) :: r, e1, dv, u(3)
      u = q(6*i-2:6*i) - q(6*i-5:6*i-3)   ! atom 2i minus atom 2i-1
      r = norm2(u)
      u = u/r
      e1 = exp(-b_mo*(r - re_mo))
      dv = 2.0d0*d_mo*b_mo*(e1 - e1*e1)
      v = v + d_mo*(e1*e1 - 2.0d0*e1)
      g(6*i-5:6*i-3) = g(6*i-5:6*i-3) - dv*u
      g(6*i-2:6*i) = g(6*i-2:6*i) + dv*u
   end subroutine add_pair
end module check_containers

program check_spectrum
   use config,         only: reactants
   use config_atoms,  only: list_atoms, list_atoms_clear
   use force_interface,   only: container_bind_pes, container_probe_force
   use spectrum,       only: stretch_mode, w_to_wvn, wvn_to_w
   use spectrum_interface, only: spectrum_frag_w, spectrum_frag_wvn, spectrum_frag_diatomic
   use check_containers
   implicit none
   integer :: nfail = 0, npass = 0
   character(len=256) :: mode
   real(8) :: fold, wexp, tol
   real(8), allocatable :: wt(:), ct(:,:), nut(:)
   integer :: k, ios
   character(len=512) :: line
   logical :: hit

   ! ---------- child dispatch: one scenario per spawn ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('fold')
         call fill_list_atoms([2])
         call container_bind_pes(unit_pes)
         block
            real(8) :: q(6), f(6)
            q = 0.0d0
            call container_probe_force(q, f)
            print '(a,es24.16)', 'FOLD ', -f(1)
         end block
      case ('morse')
         call fill_list_atoms([2])
         call container_bind_pes(morse_pes)
         call spectrum_frag_w(1, wt, ct)
         print '(a,i0)', 'N ', size(wt)
         do k = 1, size(wt)
            print '(a,i0,a,es24.16)', 'W ', k, ' ', wt(k)
         end do
      case ('bkmp2')
         call fill_list_atoms([2, 1], 1.4014719756078895d0*0.529177210903d0)
         call container_bind_pes(bkmp2_adapter)
         call spectrum_frag_wvn(1, nut)
         print '(a,i0)', 'NNU ', size(nut)
         do k = 1, size(nut)
            print '(a,i0,a,es24.16)', 'NU ', k, ' ', nut(k)
         end do
         call spectrum_frag_w(2, wt, ct)
         print '(a,i0)', 'NATOM ', size(wt)
      case ('tetra')
         call fill_list_atoms([4])
         call container_bind_pes(tetra_pes)
         call spectrum_frag_w(1, wt, ct)
         print '(a,i0)', 'N ', size(wt)
         do k = 1, size(wt)
            print '(a,i0,a,es24.16)', 'W ', k, ' ', wt(k)
         end do
      case ('escal')
         call fill_list_atoms([2, 2])
         call container_bind_pes(escal_pes)
         call spectrum_frag_w(1, wt, ct)            ! must abort (the gate)
      case ('surf')
         ! (f) paradigm neutrality: the surface label with a decaying PES derives
         reactants%surface_model = 1
         call fill_list_atoms([2])
         call container_bind_pes(morse_pes)
         call spectrum_frag_w(1, wt, ct)            ! must SUCCEED now
         print '(a,i0)', 'N ', size(wt)
         do k = 1, size(wt)
            print '(a,i0,a,es24.16)', 'W ', k, ' ', wt(k)
         end do
      case ('manual')
         ! (g) the MANUAL round trip: file staged with stretch_mode's own column,
         ! WVN and DIAT constants deliberately DISTINCT (channel separation)
         call fill_list_atoms([2])
         call write_manual_tbl('build/spec_manual.tbl')
         reactants%spectrum_src = 2
         reactants%spectrum_file = 'build/spec_manual.tbl'
         call spectrum_frag_w(1, wt, ct)
         print '(a,es24.16)', 'NU1 ', w_to_wvn(wt(1))
         block
            real(8) :: c6(6), dmax
            call stretch_mode(list_atoms%frag(1)%qz, list_atoms%mass(list_atoms%frag(1)%list), c6)
            dmax = maxval(abs(ct(:, 1) - c6))
            print '(a,es24.16)', 'CMAX ', dmax
         end block
         block
            real(8) :: we, wx, be
            call spectrum_frag_diatomic(1, we, wx, be)
            print '(a,es24.16)', 'DIATWE ', we
            print '(a,es24.16)', 'DIATWX ', wx
            print '(a,es24.16)', 'DIATBE ', be
         end block
      case ('dup')
         call fill_list_atoms([2])
         call write_dup_tbl('build/spec_dup.tbl')
         reactants%spectrum_src = 2
         reactants%spectrum_file = 'build/spec_dup.tbl'
         call spectrum_frag_w(1, wt, ct)            ! must abort (duplicate)
      case ('miss')
         call fill_list_atoms([2])
         list_atoms%symb = 'D'                      ! the file stages an H,H record
         call write_manual_tbl('build/spec_manual.tbl')
         reactants%spectrum_src = 2
         reactants%spectrum_file = 'build/spec_manual.tbl'
         call spectrum_frag_w(1, wt, ct)            ! must abort (no record)
      case ('manempty')
         call fill_list_atoms([2])
         call container_bind_pes(morse_pes)         ! irrelevant: MANUAL never probes
         reactants%spectrum_src = 2
         reactants%spectrum_file = ''
         call spectrum_frag_w(1, wt, ct)            ! must abort (empty file key)
      case ('diatc')
         call fill_list_atoms([2])
         reactants%spectrum_src = 1                 ! COMPUTE
         block
            real(8) :: we, wx, be
            call spectrum_frag_diatomic(1, we, wx, be)   ! must abort (no derivation)
         end block
      case ('ambig')
         call fill_list_atoms([2])
         call container_bind_pes(well_pes)
         call spectrum_frag_w(1, wt, ct)            ! must abort
      case default
         print '(a)', 'unknown mode'
         call exit(2)
      end select
      stop
   end if

   ! ---------- (a) the fold calibration ----------
   fold = child_real('fold', 'FOLD')
   if (fold > 0.9d0 .and. fold < 1.1d0) then       ! ~0.96485 sanity
      npass = npass + 1
      print '(a,es12.6)', '(a) fold calibration        PASS  fold = ', fold
   else
      nfail = nfail + 1
      print '(a,es12.6)', '(a) fold calibration        FAIL  fold = ', fold
   end if

   ! ---------- (b) the synthetic Morse diatomic ----------
   block
      real(8) :: got(8)
      integer :: n
      call child_ws('morse', got, n)
      wexp = sqrt(2.0d0*d_mo*b_mo*b_mo*fold/(mh/2.0d0))
      if (n >= 1) then
         tol = abs(got(1) - wexp)/wexp
      else
         tol = huge(1.0d0)
      end if
      if (n == 1 .and. tol <= 2.0d-7) then
         npass = npass + 1
         print '(a,es9.2)', '(b) Morse diatomic         PASS  rel dev = ', tol
      else
         nfail = nfail + 1
         print '(a,es9.2,i4)', '(b) Morse diatomic         FAIL  rel dev, n = ', tol, n
      end if
   end block

   ! ---------- (c) the BKMP2 H2 pair limit + the empty atom table ----------
   block
      real(8), parameter :: vpp = 35.933611286509404d0   ! eV/A^2 (recorded)
      real(8) :: got(8)
      integer :: n, nat
      call child_run('bkmp2')
      n = 0
      nat = -1
      block
         integer :: iw
         open (10, file='build/spectrum_child.txt', status='old', action='read')
         do
            read (10, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:3) == 'NU ') then
               read (line(4:), *) iw, got(min(n + 1, 8))
               n = n + 1
            else if (line(1:6) == 'NATOM ') then
               read (line(7:), *) nat
            end if
         end do
         close (10)
      end block
      wexp = w_to_wvn(sqrt(vpp*fold/(mh/2.0d0)))
      if (n >= 1) then
         tol = abs(got(1) - wexp)/wexp
      else
         tol = huge(1.0d0)
      end if
      if (n == 1 .and. tol <= 2.0d-7 .and. nat == 0) then
         npass = npass + 1
         print '(a,es9.2,a,f9.1)', '(c) BKMP2 H2 + empty atom  PASS  rel dev = ', tol, &
            '  nu = ', got(1)
      else
         nfail = nfail + 1
         print '(a,es9.2,2i4)', '(c) BKMP2 H2 + empty atom  FAIL  rel dev, n, nat = ', tol, n, nat
      end if
   end block

   ! ---------- (d) the tetrahedral spring molecule ----------
   block
      real(8) :: got(8), expw(6)
      integer :: n
      real(8), parameter :: lam(6) = [1.0d0, 1.0d0, 2.0d0, 2.0d0, 2.0d0, 4.0d0]
      call child_ws('tetra', got, n)
      expw = sqrt(lam*k_tet*fold/mh)
      if (n == 6) then
         tol = maxval(abs(got(1:6) - expw))/expw(6)
      else
         tol = huge(1.0d0)
      end if
      if (n == 6 .and. tol <= 2.0d-7) then
         npass = npass + 1
         print '(a,es9.2)', '(d) tetrahedral molecule    PASS  rel dev = ', tol
      else
         nfail = nfail + 1
         print '(a,es9.2,i4)', '(d) tetrahedral molecule    FAIL  rel dev, n = ', tol, n
      end if
   end block

   ! ---------- (e)-(g): the gate abort, paradigm neutrality, MANUAL ----------
   call abort_leg('(e) non-decayed coupling     ', 'escal', 'decayed', hit)
   if (hit) then
      npass = npass + 1
   else
      nfail = nfail + 1
   end if

   ! (f) the surface label derives through the unchanged COMPUTE path
   block
      real(8) :: got(8)
      integer :: n
      call child_ws('surf', got, n)
      wexp = sqrt(2.0d0*d_mo*b_mo*b_mo*fold/(mh/2.0d0))
      if (n >= 1) then
         tol = abs(got(1) - wexp)/wexp
      else
         tol = huge(1.0d0)
      end if
      if (n == 1 .and. tol <= 2.0d-7) then
         npass = npass + 1
         print '(a,es9.2)', '(f) surface derives        PASS  rel dev = ', tol
      else
         nfail = nfail + 1
         print '(a,es9.2,i4)', '(f) surface derives        FAIL  rel dev, n = ', tol, n
      end if
   end block

   ! (g) the MANUAL round trip (constants deliberately channel-distinct)
   block
      real(8), parameter :: wvn_row = 2000.0d0                  ! the WVN row [cm^-1]
      real(8), parameter :: dwe = 4401.21d0, dwx = 121.33d0, dbe = 60.853d0
      real(8) :: nu1, cmax, we, wx, be
      nu1 = child_real('manual', 'NU1')
      cmax = child_real2('manual', 'CMAX')
      we = child_real2('manual', 'DIATWE')
      wx = child_real2('manual', 'DIATWX')
      be = child_real2('manual', 'DIATBE')
      if (abs(nu1 - wvn_row) <= 1.0d-9 .and. cmax <= 0.0d0 .and. &
          we == dwe .and. wx == dwx .and. be == dbe) then
         npass = npass + 1
         print '(a)', '(g) MANUAL round trip      PASS'
      else
         nfail = nfail + 1
         print '(a,es12.4,es9.2)', '(g) MANUAL round trip      FAIL  nu1, cmax = ', nu1, cmax
      end if
   end block
   call abort_leg('(g) duplicate record       ', 'dup', 'one composition', hit)
   if (hit) then
      npass = npass + 1
   else
      nfail = nfail + 1
   end if
   call abort_leg('(g) carrier without record ', 'miss', 'no MANUAL spectrum record', hit)
   if (hit) then
      npass = npass + 1
   else
      nfail = nfail + 1
   end if
   call abort_leg('(g) MANUAL without file    ', 'manempty', 'MANUAL-mode requirement', hit)
   if (hit) then
      npass = npass + 1
   else
      nfail = nfail + 1
   end if
   call abort_leg('(g) EBK under COMPUTE      ', 'diatc', 'no generic derivation', hit)
   if (hit) then
      npass = npass + 1
   else
      nfail = nfail + 1
   end if

   call abort_leg('(h) ambiguous zero split    ', 'ambig', 'ambiguous', hit)
   if (hit) then
      npass = npass + 1
   else
      nfail = nfail + 1
   end if
   ! stretch_mode unit norm + structure; wavenumber round trip (pure)
   block
      real(8) :: c6(6), qz6(6), m2(2)
      qz6 = [0.1d0, -0.2d0, 0.3d0, 0.9d0, 0.1d0, 0.4d0]
      m2 = [1.008d0, 1.008d0]
      call stretch_mode(qz6, m2, c6)
      tol = abs(dot_product(c6, c6) - 1.0d0) + abs(c6(1) + c6(4)) + abs(c6(2) + c6(5)) &
            + abs(c6(3) + c6(6))
      if (tol <= 1.0d-12 .and. abs(wvn_to_w(w_to_wvn(7.3d0)) - 7.3d0) <= 1.0d-12) then
         npass = npass + 1
         print '(a)', '(h) stretch mode + wvn rt   PASS'
      else
         nfail = nfail + 1
         print '(a,es9.2)', '(h) stretch mode + wvn rt   FAIL  dev = ', tol
      end if
   end block

   print '(/,a,i0,a,i0,a)', 'check_spectrum: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) call exit(1)

contains
   ! fill_list_atoms(nats) - H atoms in fragments of the given sizes; the
   ! buffered geometries: 2 atoms = the Morse bond along x, 4 atoms = the
   ! unit-cube tetrahedron vertices, 1 atom = the origin
   subroutine fill_list_atoms(nats, bond)
      integer, intent(in) :: nats(:)
      real(8), intent(in), optional :: bond
      real(8) :: bl
      integer :: i, j, off, nt
      real(8) :: qz3(3, 4)
      call list_atoms_clear()
      nt = sum(nats)
      allocate (list_atoms%mass(nt), list_atoms%symb(nt))
      list_atoms%mass = mh
      list_atoms%symb = 'H'
      allocate (list_atoms%frag(size(nats)))
      bl = re_mo
      if (present(bond)) bl = bond
      qz3 = 0.0d0
      qz3(:, 1) = [1.0d0, 1.0d0, 1.0d0]
      qz3(:, 2) = [1.0d0, -1.0d0, -1.0d0]
      qz3(:, 3) = [-1.0d0, 1.0d0, -1.0d0]
      qz3(:, 4) = [-1.0d0, -1.0d0, 1.0d0]
      off = 0
      do i = 1, size(nats)
         associate (fr => list_atoms%frag(i))
            fr%nat = nats(i)
            allocate (fr%list(fr%nat), fr%qz(3*fr%nat))
            do j = 1, fr%nat
               fr%list(j) = off + j
            end do
            if (fr%nat == 2) then
               fr%qz = [0.0d0, 0.0d0, 0.0d0, bl, 0.0d0, 0.0d0]
            else if (fr%nat == 1) then
               fr%qz = 0.0d0
            else
               fr%qz = reshape(qz3(:, 1:fr%nat), [3*fr%nat])
            end if
            fr%mass = sum(list_atoms%mass(fr%list))
            off = off + fr%nat
         end associate
      end do
   end subroutine fill_list_atoms

   ! write_manual_tbl(path) - the (g) fixture: one H,H record whose WVN row
   ! (2000 cm^-1) is deliberately DISTINCT from its DIAT constants (the
   ! Herzberg H2 set) - a crossed channel cannot pass the round trip
   subroutine write_manual_tbl(path)
      character(len=*), intent(in) :: path
      real(8) :: c6(6)
      integer :: u
      call stretch_mode(list_atoms%frag(1)%qz, list_atoms%mass(list_atoms%frag(1)%list), c6)
      open (newunit=u, file=trim(path), status='replace', action='write')
      write (u, '(a)') '# check_spectrum manual fixture'
      write (u, '(a)') 'FRAG H H'
      write (u, '(a)') 'NMODE 1'
      write (u, '(a,1x,es24.16)') 'WVN', 2000.0d0
      write (u, '(a)') 'COL 1'
      write (u, '(6es25.16)') c6
      write (u, '(a,1x,3(1x,es24.16))') 'DIAT', 4401.21d0, 121.33d0, 60.853d0
      close (u)
   end subroutine write_manual_tbl

   ! write_dup_tbl(path) - the duplicate-composition abort fixture
   subroutine write_dup_tbl(path)
      character(len=*), intent(in) :: path
      real(8) :: c6(6)
      integer :: u
      call stretch_mode(list_atoms%frag(1)%qz, list_atoms%mass(list_atoms%frag(1)%list), c6)
      open (newunit=u, file=trim(path), status='replace', action='write')
      write (u, '(a)') 'FRAG H H'
      write (u, '(a)') 'NMODE 1'
      write (u, '(a,1x,es24.16)') 'WVN', 2000.0d0
      write (u, '(a)') 'COL 1'
      write (u, '(6es25.16)') c6
      write (u, '(a)') 'FRAG H H'
      write (u, '(a)') 'NMODE 1'
      write (u, '(a,1x,es24.16)') 'WVN', 3000.0d0
      write (u, '(a)') 'COL 1'
      write (u, '(6es25.16)') c6
      close (u)
   end subroutine write_dup_tbl

   ! child_run(mode) - respawn self with a mode; stdout/stderr to scratch
   subroutine child_run(mode)
      character(len=*), intent(in) :: mode
      character(len=256) :: self
      integer :: st
      call get_command_argument(0, self)
      call execute_command_line('"'//trim(self)//'" '//trim(mode)// &
                                ' > build/spectrum_child.txt 2> build/spectrum_child_err.txt', &
                                wait=.true., exitstat=st)
   end subroutine child_run

   ! child_real(mode, tag) - one tagged real from the child stdout
   function child_real(mode, tag) result(val)
      character(len=*), intent(in) :: mode, tag
      real(8) :: val
      integer :: ios, tl
      call child_run(mode)
      val = -huge(1.0d0)
      tl = len_trim(tag) + 1
      open (10, file='build/spectrum_child.txt', status='old', action='read')
      do
         read (10, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > tl) then
            if (line(1:tl) == trim(tag)//' ') then
               read (line(tl + 1:), *) val
               exit
            end if
         end if
      end do
      close (10)
   end function child_real

   ! child_real2(mode, tag) - like child_real but keeps the LAST match
   ! (multi-tag children: run once, read many)
   real(8) function child_real2_cached(tag) result(val)
      character(len=*), intent(in) :: tag
      integer :: ios, tl
      val = -huge(1.0d0)
      tl = len_trim(tag) + 1
      open (10, file='build/spectrum_child.txt', status='old', action='read')
      do
         read (10, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > tl) then
            if (line(1:tl) == trim(tag)//' ') then
               read (line(tl + 1:), *) val
            end if
         end if
      end do
      close (10)
   end function child_real2_cached

   real(8) function child_real2(mode, tag) result(val)
      character(len=*), intent(in) :: mode, tag
      integer :: udummy
      call child_run(mode)
      val = child_real2_cached(tag)
   end function child_real2

   ! child_ws(mode, got, n) - the W table from the child stdout (each W
   ! line carries its index and value; n cross-checked with the N line)
   subroutine child_ws(mode, got, n)
      character(len=*), intent(in) :: mode
      real(8), intent(out) :: got(8)
      integer, intent(out) :: n
      integer :: ios, m, iw
      real(8) :: val
      call child_run(mode)
      n = 0
      m = 0
      got = -huge(1.0d0)
      open (10, file='build/spectrum_child.txt', status='old', action='read')
      do
         read (10, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:2) == 'W ') then
            read (line(3:), *) iw, val
            if (iw >= 1 .and. iw <= 8) got(iw) = val
            n = n + 1
         else if (line(1:2) == 'N ') then
            read (line(3:), *) m
         end if
      end do
      close (10)
      if (m > 0 .and. m == n) n = m
   end subroutine child_ws

   ! abort_leg(label, mode, phrase, hit) - the child must exit nonzero
   ! with the phrase on stderr
   subroutine abort_leg(label, mode, phrase, hit)
      character(len=*), intent(in) :: label, mode, phrase
      logical, intent(out) :: hit
      integer :: st, ios
      hit = .false.
      call child_run(mode)
      inquire (file='build/spectrum_child_err.txt', exist=hit)
      if (.not. hit) return
      hit = .false.
      open (11, file='build/spectrum_child_err.txt', status='old', action='read')
      do
         read (11, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, trim(phrase)) > 0) hit = .true.
      end do
      close (11)
      if (hit) then
         ! nonzero exit implied by the stderr fatal line; still verify
         call child_exit(mode, st)
         if (st == 0) hit = .false.
      end if
      if (.not. hit) then
         print '(a,a,a)', trim(label), 'FAIL  (child exit/phrase missing: ', trim(mode), ')'
      end if
   end subroutine abort_leg

   ! child_exit(mode, st) - the child's exit status
   subroutine child_exit(mode, st)
      character(len=*), intent(in) :: mode
      integer, intent(out) :: st
      character(len=256) :: self
      call get_command_argument(0, self)
      call execute_command_line('"'//trim(self)//'" '//trim(mode)// &
                                ' > nul 2> nul', wait=.true., exitstat=st)
   end subroutine child_exit
end program check_spectrum
