!=====================================================================
! probe_pes.f90 - the universal standalone PES scan engine: sweep a
!                  container's potential face along 1 or 2 scan axes,
!                  write the V grid + projected analytic slopes, verify
!                  them against central differences, and hand the grid
!                  to the plotter (plot_pes.py) for the image
! Design:
!   The engine is container-agnostic: every container binds through a
!   pes_hook translation unit (the hook contract below) exposing the
!   container's POTENTIAL FACE - the same member code the production
!   chain reaches through the box aggregate (the box adds only the
!   atom-count check; the check programs hook the same faces - see
!   check_fh2_leps / check_h3_bkmp2 / check_h2ag111_pes). One hook is
!   linked per executable; the engine never changes.
!
!   Hook contract (a translation unit module named pes_hook):
!     pes_desc            - character string constant (container name;
!                           lands in the dat header and the image)
!     pes_init()          - no-argument setup: stage the container's
!                           parameter defaults / load its weight files
!                           (empty body when the face needs nothing)
!     pes_entry(q,mass,g,v) - the container_pes_i delivery contract:
!                           q [A] flattened 3N in, mass [amu] in,
!                           g = dV/dq [eV/A] out, V [eV] out
!   Hook template (fill the three marked spots):
!     module pes_hook
!        use <container pot module>, only: <face>, [<buffer calls>]
!        implicit none
!        private
!        public :: pes_init, pes_entry, pes_desc
!        character(len=*), parameter :: pes_desc = '<container + face>'
!     contains
!        subroutine pes_init()
!           ! container-specific buffer (parameter defaults, weight
!           ! files); may be empty
!        end subroutine pes_init
!        subroutine pes_entry(q, mass, g, v)
!           real(8), intent(in)  :: q(:), mass(:)
!           real(8), intent(out) :: g(:), v
!           call <face>(size(mass), q, v, g)   ! or the face's own shape
!        end subroutine pes_entry
!     end module pes_hook
!
! Control file (pes_scan.txt in cwd, or the exe's first argument; the
!   house KEYWORD = VALUE grammar - one keyword per line, '#' starts a
!   comment, blank lines ignored, unknown or duplicated keywords stop
!   with a fatal message):
!     TITLE   = <free text>                     (default 'PES scan')
!     NAT     = <n>                             (required; atoms)
!     MASS    = <n values>                      (required; amu - rides
!                                               the delivery contract)
!     Q0      = <3n values>                     (required; base
!                                               geometry [A], x y z
!                                               per atom)
!     AXIS1   = BOND i j      scan the |r_j - r_i| distance: atom j moves
!                             to distance s along the pair's Q0 direction
!                             (fixed at the base geometry - the scan
!                             stretches the bond as Q0 orients it), every
!                             other atom (i included) stays
!             = CART n i1 .. in ux uy uz
!                             place the listed atoms' mean coordinate
!                             along the (normalized) direction u at the
!                             ABSOLUTE value s (a rigid translation of
!                             the group from its base geometry - e.g. a
!                             molecule's height z above a surface)
!     RANGE1  = <min> <max> <npts>              (required; npts >= 2)
!     AXIS2 / RANGE2 - optional; a second axis makes the 2D surface
!     H_FD    = <h>      central-difference step [A] (default 1e-4)
!     FD_CHECK = T|F     gradient-vs-CD verification (default T)
!     OUT     = <stem>   output stem (default 'pes_scan' ->
!                        pes_scan.dat / pes_scan.png)
!     PLOT    = T|F      spawn the plotter at the end (default T)
!     PLOT_SCRIPT = <path>  plotter invocation (default 'plot_pes.py';
!                        resolved from cwd)
!   Axis application order: axis2 first, then axis1 (s1 is the fast
!   index). Shared-atom chains (A-B-C): the axis that MOVES the shared
!   middle atom must sit in the OUTER slot - e.g. an F-H-H cut over
!   (rHH, rFH) uses AXIS1 = BOND 2 3 and AXIS2 = BOND 1 2 (the inner
!   axis then moves only the end atom and never disturbs the outer
!   distance; the reverse order shears the grid). The projected
!   analytic slopes g.(dq/ds) carried in the dat are exact for any two
!   axes: axis2's pattern is pushed through axis1's transform Jacobian
!   (the chain_slope step - a BOND rewrite maps a perturbation of its
!   pair to the radial projection; a CART translate is position-free),
!   and the CD check audits the whole construction against V itself.
!
! Outputs (cwd):
!   <stem>.dat  - metadata header (# KEY VALUE lines; parsed by the
!                 plotter) then rows: 1D 's1 V dV/ds1'; 2D 's2 s1 V
!                 dV/ds1 dV/ds2' (s2 outer, s1 fast)
!   stdout      - scan summary: V range and argmin, the CD verification
!                 max deviation, file names
!   <stem>.png  - via the spawned plotter (python plot_pes.py <stem>)
!
! Compile (from scripts/; one hook linked in - recipes for the shipped
!   hooks; every recipe also needs plot_pes.py beside the run):
!   reference (the active container's well face; K buffered at the box
!   default 0.04 eV/A^2 - the box's run-directory parameter file is
!   NOT read at probe level; stage another constant in the hook):
!     gfortran -O2 -o probe_pes.exe ../systems/reference/reference_pot.f90 hook_reference.f90 probe_pes.f90
!   F+H2 modified LEPS (literature default set):
!     gfortran -O2 -o probe_pes.exe ../systems/inactive/fh2_leps/fh2_pot.f90 hook_fh2_leps.f90 probe_pes.f90
!   H3 BKMP2:
!     gfortran -O2 -o probe_pes.exe ../systems/inactive/h3_bkmp2/bkmp2_pot.f90 hook_h3_bkmp2.f90 probe_pes.f90
!   H2/Ag(111) NN (pes_init reads weights-h2ag.txt + biases-h2ag.txt
!   from the RUN directory - copy them from the container folder):
!     gfortran -O2 -o probe_pes.exe ../systems/inactive/h2ag111PES/h2ag111_pot.f90 hook_h2ag111.f90 probe_pes.f90
!   The mean-field containers (cau_*): bind the ground sheet face
!   (cau_<sys>_ground after cau_<sys>_load - the composed entry equals
!   it before arming; file set per the check_cau_* recipe).
!
! Run: <exe> [control file]   (default control file pes_scan.txt)
!=====================================================================
program probe_pes
   use pes_hook, only: pes_init, pes_entry, pes_desc
   implicit none

   ! ---- the scan axis ----
   type :: axis_t
      integer :: kind = 0             ! 1 = BOND, 2 = CART
      integer, allocatable :: at(:)   ! BOND: (i, j) - j moves; CART: moved-atom list
      real(8) :: u(3) = 0.0d0         ! CART direction (normalized at parse)
      real(8) :: s0 = 0.0d0           ! CART: the group's base mean coordinate along u
      real(8) :: dn(3) = 0.0d0        ! BOND: the pair's radial direction, fixed by Q0
      real(8) :: smin = 0.0d0, smax = 0.0d0
      integer :: n = 0
      character(len=96) :: label = ''
   end type axis_t

   ! ---- control state ----
   character(len=512) :: title, out_stem, plot_script, ctl_file
   character(len=512) :: line, key, val
   integer :: nat, ios, u, p, k
   real(8), allocatable :: mass(:), q0(:)
   real(8) :: h_fd
   logical :: fd_check, do_plot
   logical :: got_title, got_mass, got_q0, got_nat, got_out, got_hfd, &
              got_fd, got_plot, got_pscript
   logical :: got_ax1, got_r1, got_ax2, got_r2
   type(axis_t) :: ax1, ax2
   logical :: two_d

   ! ---- scan state ----
   real(8), allocatable :: s1v(:), s2v(:), vv(:, :), d1(:, :), d2(:, :)
   real(8), allocatable :: q(:), g(:), dq1(:), dq2(:)
   real(8) :: v, vmin, vmax, vmins1, vmins2, vmaxs1, vmaxs2
   real(8) :: fdmax, fdval, fd_s1, fd_s2
   integer :: fdax, nfd, i1, i2, st1, st2

   ! parser workspace
   character(len=64) :: tok(80)
   integer :: ntok

   ! 1. defaults
   title = 'PES scan'
   out_stem = 'pes_scan'
   plot_script = 'plot_pes.py'
   ctl_file = 'pes_scan.txt'
   h_fd = 1.0d-4
   fd_check = .true.
   do_plot = .true.
   got_title = .false.; got_mass = .false.; got_q0 = .false.; got_nat = .false.
   got_out = .false.; got_hfd = .false.; got_fd = .false.
   got_plot = .false.; got_pscript = .false.
   got_ax1 = .false.; got_r1 = .false.; got_ax2 = .false.; got_r2 = .false.
   nat = 0

   ! 2. the optional command-line control-file override
   call get_command_argument(1, ctl_file)
   if (len_trim(ctl_file) == 0) ctl_file = 'pes_scan.txt'

   ! 3. read + dispatch the control file (house KEYWORD = VALUE grammar)
   open (newunit=u, file=trim(ctl_file), status='old', action='read', iostat=ios)
   if (ios /= 0) call bail('cannot open the control file '//trim(ctl_file))
   do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      p = index(line, '#')
      if (p > 0) line = line(:p - 1)
      if (len_trim(line) == 0) cycle
      p = index(line, '=')
      if (p == 0) call bail('a control line without "=" : '//trim(line))
      key = upper(adjustl(line(:p - 1)))
      val = adjustl(line(p + 1:))
      if (len_trim(key) == 0) call bail('an empty keyword : '//trim(line))
      call split(val, tok, ntok)
      if (ntok == 0) call bail('an empty value for '//trim(key))

      select case (trim(key))
      case ('TITLE')
         if (got_title) call bail('duplicated TITLE')
         got_title = .true.
         title = val
      case ('NAT')
         if (got_nat) call bail('duplicated NAT')
         got_nat = .true.
         if (ntok /= 1) call bail('NAT needs exactly one integer')
         read (tok(1), *, iostat=ios) nat
         if (ios /= 0 .or. nat < 1) call bail('NAT must be a positive integer')
      case ('MASS')
         if (got_mass) call bail('duplicated MASS')
         got_mass = .true.
         allocate (mass(ntok))
         do k = 1, ntok
            read (tok(k), *, iostat=ios) mass(k)
            if (ios /= 0 .or. mass(k) <= 0) call bail('MASS needs positive reals')
         end do
      case ('Q0')
         if (got_q0) call bail('duplicated Q0')
         got_q0 = .true.
         allocate (q0(ntok))
         do k = 1, ntok
            read (tok(k), *, iostat=ios) q0(k)
            if (ios /= 0) call bail('Q0 needs reals')
         end do
      case ('AXIS1')
         if (got_ax1) call bail('duplicated AXIS1')
         got_ax1 = .true.
         call parse_axis(val, tok, ntok, ax1, 1)
      case ('AXIS2')
         if (got_ax2) call bail('duplicated AXIS2')
         got_ax2 = .true.
         call parse_axis(val, tok, ntok, ax2, 2)
      case ('RANGE1')
         if (got_r1) call bail('duplicated RANGE1')
         got_r1 = .true.
         call parse_range(tok, ntok, ax1)
      case ('RANGE2')
         if (got_r2) call bail('duplicated RANGE2')
         got_r2 = .true.
         call parse_range(tok, ntok, ax2)
      case ('H_FD')
         if (got_hfd) call bail('duplicated H_FD')
         got_hfd = .true.
         read (tok(1), *, iostat=ios) h_fd
         if (ios /= 0 .or. h_fd <= 0) call bail('H_FD must be a positive real')
      case ('FD_CHECK')
         if (got_fd) call bail('duplicated FD_CHECK')
         got_fd = .true.
         fd_check = parse_bool(tok(1), 'FD_CHECK')
      case ('PLOT')
         if (got_plot) call bail('duplicated PLOT')
         got_plot = .true.
         do_plot = parse_bool(tok(1), 'PLOT')
      case ('OUT')
         if (got_out) call bail('duplicated OUT')
         got_out = .true.
         out_stem = val
      case ('PLOT_SCRIPT')
         if (got_pscript) call bail('duplicated PLOT_SCRIPT')
         got_pscript = .true.
         plot_script = val
      case default
         call bail('unknown keyword '//trim(key)//' in '//trim(ctl_file))
      end select
   end do
   close (u)

   ! 4. consistency gate
   two_d = got_ax2 .or. got_r2
   if (two_d .and. .not. (got_ax2 .and. got_r2)) &
      call bail('a second axis needs both AXIS2 and RANGE2')
   if (.not. (got_nat .and. got_mass .and. got_q0 .and. got_ax1 .and. got_r1)) &
      call bail('required keys: NAT, MASS, Q0, AXIS1, RANGE1')
   if (size(mass) /= nat) call bail('MASS needs exactly NAT values')
   if (size(q0) /= 3*nat) call bail('Q0 needs exactly 3*NAT values (x y z per atom)')
   call check_axis(ax1, nat, 'AXIS1')
   if (two_d) call check_axis(ax2, nat, 'AXIS2')
   call axis_base(ax1, q0)
   if (two_d) call axis_base(ax2, q0)

   ! 5. the grid axes
   allocate (s1v(ax1%n))
   do i1 = 1, ax1%n
      s1v(i1) = ax1%smin + (ax1%smax - ax1%smin)*real(i1 - 1, 8)/real(ax1%n - 1, 8)
   end do
   if (two_d) then
      allocate (s2v(ax2%n))
      do i2 = 1, ax2%n
         s2v(i2) = ax2%smin + (ax2%smax - ax2%smin)*real(i2 - 1, 8)/real(ax2%n - 1, 8)
      end do
   else
      allocate (s2v(1))
      s2v(1) = 0.0d0
   end if
   allocate (vv(ax1%n, size(s2v)), d1(ax1%n, size(s2v)), d2(ax1%n, size(s2v)))
   allocate (q(3*nat), g(3*nat), dq1(3*nat), dq2(3*nat))

   ! 6. container init + the scan (axis2 outer/slow, axis1 inner/fast)
   call pes_init()
   write (*, '(a)') 'probe_pes: container = '//trim(pes_desc)
   write (*, '(a)') 'probe_pes: title     = '//trim(title)
   write (*, '(a,i0,a,i0,a)') 'probe_pes: grid     = ', ax1%n, ' x ', size(s2v), ' points'
   do i2 = 1, size(s2v)
      do i1 = 1, ax1%n
         call eval_point(s1v(i1), s2v(i2), v, dq1, dq2)
         vv(i1, i2) = v
         d1(i1, i2) = dot_product(g, dq1)
         if (two_d) d2(i1, i2) = dot_product(g, dq2)
      end do
   end do

   ! 7. V range + argmin / argmax
   vmin = vv(1, 1); vmax = vv(1, 1)
   vmins1 = s1v(1); vmins2 = s2v(1); vmaxs1 = s1v(1); vmaxs2 = s2v(1)
   do i2 = 1, size(s2v)
      do i1 = 1, ax1%n
         if (vv(i1, i2) < vmin) then
            vmin = vv(i1, i2); vmins1 = s1v(i1); vmins2 = s2v(i2)
         end if
         if (vv(i1, i2) > vmax) then
            vmax = vv(i1, i2); vmaxs1 = s1v(i1); vmaxs2 = s2v(i2)
         end if
      end do
   end do

   ! 8. the CD verification: analytic projected slope vs central
   !    difference of V along each axis (interior points; 2D subsampled
   !    to at most ~9x9 probe points per axis)
   fdmax = -1.0d0; fdax = 0; fd_s1 = 0; fd_s2 = 0; nfd = 0
   if (fd_check) then
      st1 = max(1, (ax1%n - 1)/8); st2 = max(1, (size(s2v) - 1)/8)
      do i2 = 1, size(s2v), st2
         if (two_d .and. (i2 == 1 .or. i2 == size(s2v))) cycle
         do i1 = 1, ax1%n
            if (i1 == 1 .or. i1 == ax1%n) cycle
            if (two_d .and. mod(i1 - 1, st1) /= 0) cycle
            ! axis 1
            call eval_point(s1v(i1) + h_fd, s2v(i2), v, dq1, dq2)
            fdval = v
            call eval_point(s1v(i1) - h_fd, s2v(i2), v, dq1, dq2)
            fdval = (fdval - v)/(2*h_fd)
            fdval = abs(fdval - d1(i1, i2))
            nfd = nfd + 1
            if (fdval > fdmax) then
               fdmax = fdval; fdax = 1; fd_s1 = s1v(i1); fd_s2 = s2v(i2)
            end if
            ! axis 2
            if (two_d) then
               call eval_point(s1v(i1), s2v(i2) + h_fd, v, dq1, dq2)
               fdval = v
               call eval_point(s1v(i1), s2v(i2) - h_fd, v, dq1, dq2)
               fdval = (fdval - v)/(2*h_fd)
               fdval = abs(fdval - d2(i1, i2))
               nfd = nfd + 1
               if (fdval > fdmax) then
                  fdmax = fdval; fdax = 2; fd_s1 = s1v(i1); fd_s2 = s2v(i2)
               end if
            end if
         end do
      end do
   end if

   ! 9. write the dat (metadata header parsed by plot_pes.py)
   open (newunit=u, file=trim(out_stem)//'.dat', status='replace', action='write')
   write (u, '(a)') '# PROBE_PES 1'
   write (u, '(a)') '# CONTAINER '//trim(pes_desc)
   write (u, '(a)') '# TITLE '//trim(title)
   write (u, '(a,i0)') '# NAXES ', merge(2, 1, two_d)
   write (u, '(a,i0)') '# NAT ', nat
   write (u, '(a)') '# LABEL1 '//trim(ax1%label)
   if (two_d) write (u, '(a)') '# LABEL2 '//trim(ax2%label)
   write (u, '(a,i0)') '# N1 ', ax1%n
   if (two_d) write (u, '(a,i0)') '# N2 ', ax2%n
   write (u, '(a,2(1x,es16.8))') '# RANGE1 ', ax1%smin, ax1%smax
   if (two_d) write (u, '(a,2(1x,es16.8))') '# RANGE2 ', ax2%smin, ax2%smax
   write (u, '(a,es16.8)') '# FD_H ', h_fd
   if (fd_check) then
      if (fdmax >= 0.0d0) then
         write (u, '(a,i0,a,es16.8,a,i0,a,2(1x,es16.8))') &
            '# FD_MAX ', nfd, ' points; max |g.proj - CD| =', fdmax, &
            ' eV/A at axis ', fdax, ' s =', fd_s1, fd_s2
      else
         write (u, '(a)') '# FD_MAX none (no interior probe points)'
      end if
   else
      write (u, '(a)') '# FD_MAX none (FD_CHECK = F)'
   end if
   write (u, '(a)') '# UNITS s [A]; V [eV]; dV/ds [eV/A]'
   if (two_d) then
      write (u, '(a)') '# COLUMNS s2 s1 V dV/ds1 dV/ds2 (s2 outer, s1 fast)'
      do i2 = 1, size(s2v)
         do i1 = 1, ax1%n
            write (u, '(3(1x,es16.8),2(1x,es15.7))') s2v(i2), s1v(i1), vv(i1, i2), &
               d1(i1, i2), d2(i1, i2)
         end do
      end do
   else
      write (u, '(a)') '# COLUMNS s1 V dV/ds1'
      do i1 = 1, ax1%n
         write (u, '(2(1x,es16.8),1x,es15.7)') s1v(i1), vv(i1, 1), d1(i1, 1)
      end do
   end if
   close (u)

   ! 10. stdout summary
   write (*, '(a)') 'probe_pes: scan of '//trim(ax1%label)// &
      merge(' x '//trim(ax2%label), '', two_d)
   write (*, '(a,es16.8,a,es16.8)') 'probe_pes: V range  [eV]  min ', vmin, '  max ', vmax
   write (*, '(a,es16.8,a,es16.8)') 'probe_pes: argmin at s1 ', vmins1, '  s2 ', vmins2
   write (*, '(a,es16.8,a,es16.8)') 'probe_pes: argmax at s1 ', vmaxs1, '  s2 ', vmaxs2
   if (fd_check) then
      write (*, '(a,i0,a,es12.4,a,i0,a)') 'probe_pes: CD check  ', nfd, &
         ' probes, max |g.proj - CD| = ', fdmax, ' eV/A (axis ', fdax, ')'
   else
      write (*, '(a)') 'probe_pes: CD check  skipped (FD_CHECK = F)'
   end if
   write (*, '(a)') 'probe_pes: wrote '//trim(out_stem)//'.dat'

   ! 11. hand the grid to the plotter (auxiliary - a failure is a
   !     notice, never a scan failure)
   if (do_plot) then
      call run_plotter(trim(plot_script), trim(out_stem))
   else
      write (*, '(a)') 'probe_pes: plotting skipped (PLOT = F)'
   end if

! ---- internal: geometry assembly + one container evaluation ----
contains

   ! eval_point(s1, s2, v_out, dq1, dq2) - build q(s1,s2) from q0 (axis2
   !   first, axis1 last), evaluate the container face, and hand back
   !   the slope patterns of both axes (axis2's chained through axis1)
   subroutine eval_point(s1, s2, v_out, dq_1, dq_2)
      real(8), intent(in) :: s1, s2
      real(8), intent(out) :: v_out
      real(8), intent(out) :: dq_1(:), dq_2(:)
      q = q0
      dq_2 = 0.0d0
      if (two_d) then
         call apply_axis(ax2, s2, q)
         call slope_pattern(ax2, dq_2)
      end if
      call apply_axis(ax1, s1, q)
      if (two_d) call chain_slope(ax1, dq_2)
      call slope_pattern(ax1, dq_1)
      call pes_entry(q, mass, g, v_out)
   end subroutine eval_point

   ! chain_slope(a, dq) - push a slope pattern through the Jacobian of
   !   a LATER-applied axis transform (axis2's pattern meets axis1 after
   !   axis1 rewrote the geometry): a BOND rewrite r_j = r_i + dn*s turns
   !   a pre-image perturbation (dpi, dpj) of its pair into
   !   dpj' = dpi + (I - dn dn^T)(dpj - dpi) (dn fixed by Q0); a CART
   !   translate is position-free (identity). Other atoms are untouched
   !   by either transform
   subroutine chain_slope(a, dq)
      type(axis_t), intent(in) :: a
      real(8), intent(inout) :: dq(:)
      real(8) :: pi3(3), rel(3)
      integer :: ia, ib
      if (a%kind /= 1) return
      ia = 3*a%at(1) - 2; ib = 3*a%at(2) - 2
      pi3 = dq(ia:ia + 2)
      rel = dq(ib:ib + 2) - pi3
      rel = rel - a%dn*dot_product(a%dn, rel)
      dq(ib:ib + 2) = pi3 + rel
   end subroutine chain_slope

   ! apply_axis(a, s, q) - rewrite q by the axis transform at scan value s
   subroutine apply_axis(a, s, q_io)
      type(axis_t), intent(in) :: a
      real(8), intent(in) :: s
      real(8), intent(inout) :: q_io(:)
      integer :: ia, ib
      if (a%kind == 1) then
         ia = 3*a%at(1) - 2; ib = 3*a%at(2) - 2
         q_io(ib:ib + 2) = q_io(ia:ia + 2) + a%dn*s
      else
         do k = 1, size(a%at)
            ib = 3*a%at(k) - 2
            q_io(ib:ib + 2) = q_io(ib:ib + 2) + a%u*(s - a%s0)
         end do
      end if
   end subroutine apply_axis

   ! slope_pattern(a, dq) - dq/ds of the axis transform (BOND: the Q0
   !   radial direction on the moving atom; CART: u on every listed atom)
   subroutine slope_pattern(a, dq)
      type(axis_t), intent(in) :: a
      real(8), intent(out) :: dq(:)
      integer :: ib
      dq = 0.0d0
      if (a%kind == 1) then
         ib = 3*a%at(2) - 2
         dq(ib:ib + 2) = a%dn
      else
         do k = 1, size(a%at)
            ib = 3*a%at(k) - 2
            dq(ib:ib + 2) = a%u
         end do
      end if
   end subroutine slope_pattern

   ! ---- control-file parsing helpers ----

   ! axis_base(a, q) - the anchors from the BASE geometry: a CART axis
   !   fixes the listed group's mean coordinate along u (s is absolute
   !   from here); a BOND axis fixes the pair's radial direction dn (the
   !   scan stretches the bond along its Q0 orientation - dn never
   !   re-derives from the swept geometry, so no axis can flip another
   !   axis's direction mid-scan)
   subroutine axis_base(a, q_in)
      type(axis_t), intent(inout) :: a
      real(8), intent(in) :: q_in(:)
      real(8) :: d(3), rn
      integer :: ib
      select case (a%kind)
      case (1)
         d = q_in(3*a%at(2) - 2:3*a%at(2)) - q_in(3*a%at(1) - 2:3*a%at(1))
         rn = sqrt(dot_product(d, d))
         if (rn <= 0.0d0) call bail('a BOND axis has a zero base distance in Q0')
         a%dn = d/rn
      case (2)
         a%s0 = 0.0d0
         do k = 1, size(a%at)
            ib = 3*a%at(k) - 2
            a%s0 = a%s0 + dot_product(a%u, q_in(ib:ib + 2))
         end do
         a%s0 = a%s0/size(a%at)
      end select
   end subroutine axis_base

   subroutine parse_axis(val_s, tok_a, nt, a, idx)
      character(len=*), intent(in) :: val_s
      character(len=64), intent(in) :: tok_a(:)
      integer, intent(in) :: nt, idx
      type(axis_t), intent(out) :: a
      character(len=8) :: word
      real(8) :: un
      integer :: m, kk, ios2
      if (nt < 1) call bail('an empty AXIS')
      word = upper(adjustl(tok_a(1)))
      select case (trim(word))
      case ('BOND')
         if (nt /= 3) call bail('AXIS BOND needs: BOND i j')
         a%kind = 1
         allocate (a%at(2))
         read (tok_a(2), *, iostat=ios2) a%at(1)
         if (ios2 /= 0) call bail('AXIS BOND needs two integer atom indices')
         read (tok_a(3), *, iostat=ios2) a%at(2)
         if (ios2 /= 0) call bail('AXIS BOND needs two integer atom indices')
         write (a%label, '(a,i0,a,i0,a)') 'r(a', a%at(2), '-a', a%at(1), ') [A]'
      case ('CART')
         if (nt < 6) call bail('AXIS CART needs: CART n i1..in ux uy uz')
         read (tok_a(2), *, iostat=ios2) m
         if (ios2 /= 0 .or. m < 1 .or. nt /= 5 + m) &
            call bail('AXIS CART token count mismatch (CART n i1..in ux uy uz)')
         a%kind = 2
         allocate (a%at(m))
         do kk = 1, m
            read (tok_a(2 + kk), *, iostat=ios2) a%at(kk)
            if (ios2 /= 0) call bail('AXIS CART needs integer atom indices')
         end do
         read (tok_a(3 + m), *, iostat=ios2) a%u(1)
         read (tok_a(4 + m), *, iostat=ios2) a%u(2)
         read (tok_a(5 + m), *, iostat=ios2) a%u(3)
         if (ios2 /= 0) call bail('AXIS CART needs a direction vector ux uy uz')
         un = sqrt(dot_product(a%u, a%u))
         if (un <= 0.0d0) call bail('AXIS CART direction is the zero vector')
         a%u = a%u/un
         write (a%label, '(a,i0,a)') 'group coord (', m, ' atoms) [A]'
      case default
         call bail('AXIS must start with BOND or CART, got '//trim(val_s))
      end select
      if (a%kind == 1 .and. a%at(1) == a%at(2)) call bail('a BOND axis needs two distinct atoms')
      if (idx == 1 .and. a%kind == 0) call bail('AXIS1 parse failure')
   end subroutine parse_axis

   subroutine parse_range(tok_a, nt, a)
      character(len=64), intent(in) :: tok_a(:)
      integer, intent(in) :: nt
      type(axis_t), intent(inout) :: a
      integer :: ios2
      if (nt /= 3) call bail('a RANGE needs: min max npts')
      read (tok_a(1), *, iostat=ios2) a%smin
      read (tok_a(2), *, iostat=ios2) a%smax
      if (ios2 /= 0 .or. a%smax <= a%smin) call bail('a RANGE needs min < max')
      read (tok_a(3), *, iostat=ios2) a%n
      if (ios2 /= 0 .or. a%n < 2) call bail('a RANGE needs npts >= 2')
   end subroutine parse_range

   subroutine check_axis(a, n_at, name)
      type(axis_t), intent(in) :: a
      integer, intent(in) :: n_at
      character(len=*), intent(in) :: name
      if (a%kind == 0) call bail(name//' was not parsed')
      if (a%n < 2) call bail(name//' has no RANGE')
      do k = 1, size(a%at)
         if (a%at(k) < 1 .or. a%at(k) > n_at) &
            call bail(name//' references an atom outside 1..NAT')
      end do
   end subroutine check_axis

   logical function parse_bool(w, name)
      character(len=*), intent(in) :: w, name
      select case (trim(upper(w)))
      case ('T', '.TRUE.', 'TRUE')
         parse_bool = .true.
      case ('F', '.FALSE.', 'FALSE')
         parse_bool = .false.
      case default
         call bail(name//' needs T or F')
      end select
   end function parse_bool

   ! split(v, tok, n) - whitespace/comma tokenizer
   subroutine split(v, tok_o, n)
      character(len=*), intent(in) :: v
      character(len=64), intent(out) :: tok_o(:)
      integer, intent(out) :: n
      character(len=len(v)) :: w
      integer :: p, c
      n = 0
      w = adjustl(v)
      do while (len_trim(w) > 0)
         p = scan(w, ' ,')
         if (p == 0) p = len_trim(w) + 1
         n = n + 1
         if (n > size(tok_o)) call bail('a control value carries too many tokens')
         tok_o(n) = adjustl(w(:p - 1))
         w = adjustl(w(p:))
         ! drop the separator characters just consumed
         do while (len_trim(w) > 0 .and. verify(w(1:1), ' ,') == 0)
            w = adjustl(w(2:))
         end do
         c = len_trim(w)
         if (c == 0) exit
      end do
   end subroutine split

   function upper(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, ic
      do i = 1, len(s)
         ic = iachar(s(i:i))
         if (ic >= 97 .and. ic <= 122) then
            t(i:i) = achar(ic - 32)
         else
            t(i:i) = s(i:i)
         end if
      end do
   end function upper

   ! run_plotter - spawn 'python <script> <stem>'; a notice on failure
   subroutine run_plotter(script, stem)
      character(len=*), intent(in) :: script, stem
      character(len=1024) :: cmd
      integer :: cs, es
      cmd = 'python "'//trim(script)//'" "'//trim(stem)//'"'
      write (*, '(a)') 'probe_pes: plotting -> '//trim(cmd)
      call execute_command_line(trim(cmd), wait=.true., cmdstat=cs, exitstat=es)
      if (cs /= 0) then
         write (*, '(a)') 'probe_pes: notice: the plotter could not start (python '// &
            'not found or script missing) - the dat stands alone; plot '// &
            'later with:  python plot_pes.py '//trim(stem)
      else if (es /= 0) then
         write (*, '(a,i0)') 'probe_pes: notice: the plotter exited nonzero (', es, &
            ') - the dat stands alone'
      else
         write (*, '(a)') 'probe_pes: image written ('//trim(stem)//'.png)'
      end if
   end subroutine run_plotter

   ! bail(who, msg) - the named fatal channel
   subroutine bail(msg)
      character(len=*), intent(in) :: msg
      write (0, '(a)') 'probe_pes: fatal: '//trim(msg)
      write (0, '(a)') 'probe_pes: fatal: the scan did not run (control file: '// &
         trim(ctl_file)//')'
      stop 1
   end subroutine bail

end program probe_pes
