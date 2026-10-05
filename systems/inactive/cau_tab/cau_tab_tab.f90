!=====================================================================
! cau_tab_tab.f90 - the DFT-point table storage: loader with
!   self-checks + per-site 1D Akima interpolation (values + slopes)
! Design:
!   The table file (format: TABLE_FORMAT.md in this folder) holds five
!   channel sheets (E0 / h / U / J / VCP) x three lateral sites (top /
!   bridge / hollow), each a z-ascending point list in A / eV. The
!   channels are read DIRECTLY (2026-10-01 ruling): the former
!   four-diabat subtraction (J = E_1D-E_3P etc.) stacked each sheet's
!   interpolation error and C1 kinks into the channels - oscillatory;
!   the GS DFT data is reference-only, the mechanical base is the E0
!   sheet. Loading
!   builds the Akima (1970) point slopes per (sheet, site); tab_eval
!   evaluates one sheet-site pair (value + dz slope) by bisection +
!   cubic Hermite. Range policy: below the site's first z -> fatal
!   (never extrapolate physics the hand-edited table does not carry);
!   above the last z -> hold the end value with zero slope (the
!   asymptote hold, consistent with the 6 A freeze-line ruling of
!   docs/plans/2026-09-28-coupling-kernel-fix.md).
!   The table file IS the experiment knob (2026-09-30 ruling): a run
!   points TABLE_FILE at a hand-edited copy; the load echo gives every
!   run's stdout its provenance.
!=====================================================================
module cau_tab_tab
   implicit none
   private
   public :: tab_load, tab_eval, tab_loaded, tab_lattice, tab_info_out
   public :: sh_E0, sh_h, sh_U, sh_J, sh_VCP
   public :: si_top, si_bri, si_hol

   integer, parameter :: sh_E0 = 1, sh_h = 2, sh_U = 3
   integer, parameter :: sh_J = 4, sh_VCP = 5
   integer, parameter :: si_top = 1, si_bri = 2, si_hol = 3
   integer, parameter :: n_sheet = 5, n_site = 3, max_pts = 1024

   character(len=8), parameter :: sh_name(n_sheet) = ['E0 ', 'h  ', 'U  ', 'J  ', 'VCP']
   character(len=8), parameter :: si_name(n_site) = ['top   ', 'bridge', 'hollow']

   real(8) :: lattice = 0.0d0                   ! lateral period [A] (the LATTICE row)
   real(8) :: site_xy(n_site, 2) = 0.0d0        ! site coordinates [A] (cross-sheet checked)
   integer :: t_n(n_sheet, n_site) = 0          ! point counts
   real(8), allocatable :: t_z(:, :, :)         ! (sheet, site, point) z [A]
   real(8), allocatable :: t_e(:, :, :)         ! energies [eV]
   real(8), allocatable :: t_m(:, :, :)         ! Akima point slopes [eV/A]
   logical :: loaded = .false.

contains

   !------------------------------------------------------------------
   ! tab_load(file) - parse + self-check + build the Akima slopes
   !                 (re-loadable: the criterion program loads twice)
   !------------------------------------------------------------------
   subroutine tab_load(file, quiet)
      character(len=*), intent(in) :: file
      logical, intent(in), optional :: quiet
      logical :: echo_on
      character(len=512) :: line
      character(len=16) :: kw, name
      real(8) :: wz(n_sheet, n_site, max_pts), we(n_sheet, n_site, max_pts)
      real(8) :: xy(n_site, 2)
      logical :: seen_sh(n_sheet), seen_si(n_sheet, n_site), xy_known(n_site)
      integer :: wn(n_sheet, n_site)
      integer :: u, ios, cur_sh, cur_si, ish, isi, k, mx
      real(8) :: aval, xv, yv, zv, ev
      logical :: fexist

      echo_on = .true.
      if (present(quiet)) echo_on = .not. quiet
      inquire (file=trim(file), exist=fexist)
      if (.not. fexist) then
         call stop_tab('tab_load', 'cannot find the table file "'//trim(file)//'"')
      end if
      wz = 0.0d0
      we = 0.0d0
      xy = 0.0d0
      wn = 0
      seen_sh = .false.
      seen_si = .false.
      xy_known = .false.
      lattice = 0.0d0
      cur_sh = 0
      cur_si = 0
      open (newunit=u, file=trim(file), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call stop_tab('tab_load', 'cannot open the table file "'//trim(file)//'"')
      end if
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         read (line, *, iostat=ios) kw
         if (ios /= 0) cycle                          ! unreachable blank guard
         if (trim(kw) == 'LATTICE') then
            if (lattice > 0.0d0) then
               call stop_tab('tab_load', 'duplicate LATTICE row')
            end if
            read (line, *, iostat=ios) kw, aval
            if (ios /= 0 .or. aval <= 0.0d0) then
               call stop_tab('tab_load', 'LATTICE needs one positive period [A]')
            end if
            lattice = aval
         else if (trim(kw) == 'SHEET') then
            read (line, *, iostat=ios) kw, name
            if (ios /= 0) call stop_tab('tab_load', 'SHEET needs a name (E0/h/U/J/VCP)')
            cur_sh = sheet_index(name)
            if (cur_sh == 0) then
               call stop_tab('tab_load', 'unknown sheet name "'//trim(name)//'"')
            end if
            if (seen_sh(cur_sh)) then
               call stop_tab('tab_load', 'duplicate SHEET '//trim(sh_name(cur_sh)))
            end if
            seen_sh(cur_sh) = .true.
            cur_si = 0
         else if (trim(kw) == 'SITE') then
            read (line, *, iostat=ios) kw, name, xv, yv
            if (ios /= 0) call stop_tab('tab_load', 'SITE needs a label and x y [A]')
            if (cur_sh == 0) then
               call stop_tab('tab_load', 'SITE row outside any SHEET block')
            end if
            cur_si = site_index(name)
            if (cur_si == 0) then
               call stop_tab('tab_load', 'unknown site label "'//trim(name)//'" (top/bridge/hollow)')
            end if
            if (seen_si(cur_sh, cur_si)) then
               call stop_tab('tab_load', 'duplicate SITE '//trim(si_name(cur_si))//' in sheet '// &
                             trim(sh_name(cur_sh)))
            end if
            seen_si(cur_sh, cur_si) = .true.
            if (xy_known(cur_si)) then
               if (abs(xy(cur_si, 1) - xv) > 1.0d-3 .or. abs(xy(cur_si, 2) - yv) > 1.0d-3) then
                  write (0, '(a,2f10.4,a,2f10.4)') 'tab_load: first-seen ', xy(cur_si, :), &
                     ' vs this row ', xv, yv
                  call stop_tab('tab_load', 'site '//trim(si_name(cur_si))// &
                                ' has inconsistent coordinates across sheets')
               end if
            else
               xy(cur_si, :) = [xv, yv]
               xy_known(cur_si) = .true.
            end if
         else
            ! 1. a data row: z [A]  E [eV]
            read (line, *, iostat=ios) zv, ev
            if (ios /= 0) then
               call stop_tab('tab_load', 'unrecognized row "'//trim(line)//'"')
            end if
            if (cur_sh == 0 .or. cur_si == 0) then
               call stop_tab('tab_load', 'data row outside a SITE block')
            end if
            if (wn(cur_sh, cur_si) >= max_pts) then
               call stop_tab('tab_load', 'more than max_pts points in one site block')
            end if
            if (wn(cur_sh, cur_si) > 0) then
               if (zv <= wz(cur_sh, cur_si, wn(cur_sh, cur_si))) then
                  call stop_tab('tab_load', 'z not strictly ascending in sheet '// &
                                trim(sh_name(cur_sh))//' site '//trim(si_name(cur_si)))
               end if
               if (zv - wz(cur_sh, cur_si, wn(cur_sh, cur_si)) < 1.0d-4) then
                  write (0, '(a,2f12.6)') 'tab_load: the segment is thinner than 1e-4 A: ', &
                     wz(cur_sh, cur_si, wn(cur_sh, cur_si)), zv
                  call stop_tab('tab_load', 'sheet '//trim(sh_name(cur_sh))//' site '// &
                                trim(si_name(cur_si))//' carries a near-degenerate z pair'// &
                                ' (two calculation passes interleaved? merge them - a'// &
                                ' hairline segment blows the Akima slopes)')
               end if
            end if
            wn(cur_sh, cur_si) = wn(cur_sh, cur_si) + 1
            wz(cur_sh, cur_si, wn(cur_sh, cur_si)) = zv
            we(cur_sh, cur_si, wn(cur_sh, cur_si)) = ev
         end if
      end do
      close (u)

      ! 2. completeness + geometry self-checks
      if (lattice <= 0.0d0) then
         call stop_tab('tab_load', 'the LATTICE row is missing')
      end if
      do ish = 1, n_sheet
         if (.not. seen_sh(ish)) then
            call stop_tab('tab_load', 'the SHEET '//trim(sh_name(ish))//' block is missing')
         end if
         do isi = 1, n_site
            if (.not. seen_si(ish, isi)) then
               call stop_tab('tab_load', 'the SITE '//trim(si_name(isi))//' block of sheet '// &
                             trim(sh_name(ish))//' is missing')
            end if
            if (wn(ish, isi) < 4) then
               call stop_tab('tab_load', 'sheet '//trim(sh_name(ish))//' site '// &
                             trim(si_name(isi))//' carries fewer than 4 points')
            end if
         end do
      end do

      ! 3. move into the module state + build the Akima slopes
      if (allocated(t_z)) deallocate (t_z, t_e, t_m)
      mx = maxval(wn)
      allocate (t_z(n_sheet, n_site, mx), t_e(n_sheet, n_site, mx), t_m(n_sheet, n_site, mx))
      t_n = wn
      site_xy = xy
      do ish = 1, n_sheet
         do isi = 1, n_site
            t_z(ish, isi, 1:wn(ish, isi)) = wz(ish, isi, 1:wn(ish, isi))
            t_e(ish, isi, 1:wn(ish, isi)) = we(ish, isi, 1:wn(ish, isi))
            call akima_build(wn(ish, isi), t_z(ish, isi, :), t_e(ish, isi, :), t_m(ish, isi, :))
         end do
      end do
      loaded = .true.

      ! 4. the load echo (every run's stdout carries its table provenance;
      !    the sparse-grid maker reloads per greedy round and silences it)
      if (echo_on) then
         write (*, '(a,a,a,f8.4,a)') 'cau_tab_tab: loaded "', trim(file), '"  lattice ', lattice, ' A'
         do ish = 1, n_sheet
            do isi = 1, n_site
               write (*, '(a,a,a,a,a,i0,a,f7.3,a,f7.3,a)') '  sheet ', sh_name(ish), '  site ', &
                  si_name(isi), ': ', t_n(ish, isi), ' pts  z [', t_z(ish, isi, 1), ', ', &
                  t_z(ish, isi, t_n(ish, isi)), '] A'
            end do
         end do
      end if
   end subroutine tab_load

   !------------------------------------------------------------------
   ! tab_eval(isheet, isite, z, e, de) - one sheet-site Akima evaluation
   !                 (value + dz slope; bisection on the ascending grid)
   !------------------------------------------------------------------
   subroutine tab_eval(isheet, isite, z, e, de)
      integer, intent(in) :: isheet, isite
      real(8), intent(in) :: z
      real(8), intent(out) :: e, de
      integer :: n, lo, hi, mid, i
      real(8) :: h, t, z0, e0, m0, z1, e1, m1
      if (.not. loaded) then
         call stop_tab('tab_eval', 'the table is not loaded (seam order error)')
      end if
      n = t_n(isheet, isite)
      if (z < t_z(isheet, isite, 1)) then
         write (0, '(a,f10.4,a,f10.4)') 'tab_eval: requested z = ', z, &
            ' below the table minimum ', t_z(isheet, isite, 1)
         call stop_tab('tab_eval', 'sheet '//trim(sh_name(isheet))//' site '// &
                       trim(si_name(isite))//' has no data below its first point - '// &
                       'extend the table (the container never extrapolates)')
      end if
      if (z >= t_z(isheet, isite, n)) then
         ! the asymptote hold: end value, zero slope
         e = t_e(isheet, isite, n)
         de = 0.0d0
         return
      end if
      lo = 1
      hi = n - 1
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (z >= t_z(isheet, isite, mid)) then
            lo = mid
         else
            hi = mid
         end if
      end do
      i = lo
      z0 = t_z(isheet, isite, i)
      z1 = t_z(isheet, isite, i + 1)
      e0 = t_e(isheet, isite, i)
      e1 = t_e(isheet, isite, i + 1)
      m0 = t_m(isheet, isite, i)
      m1 = t_m(isheet, isite, i + 1)
      h = z1 - z0
      t = (z - z0)/h
      ! cubic Hermite on the unit segment + its dz derivative
      e = (2.0d0*t**3 - 3.0d0*t**2 + 1.0d0)*e0 + (t**3 - 2.0d0*t**2 + t)*h*m0 &
        + (-2.0d0*t**3 + 3.0d0*t**2)*e1 + (t**3 - t**2)*h*m1
      de = ((6.0d0*t**2 - 6.0d0*t)*e0 + (3.0d0*t**2 - 4.0d0*t + 1.0d0)*h*m0 &
           + (-6.0d0*t**2 + 6.0d0*t)*e1 + (3.0d0*t**2 - 2.0d0*t)*h*m1)/h
   end subroutine tab_eval

   logical function tab_loaded()
      tab_loaded = loaded
   end function tab_loaded

   real(8) function tab_lattice()
      tab_lattice = lattice
   end function tab_lattice

   !------------------------------------------------------------------
   ! tab_info_out(lat, npts) - the criterion introspection (loader legs)
   !------------------------------------------------------------------
   subroutine tab_info_out(lat, npts)
      real(8), intent(out) :: lat
      integer, intent(out) :: npts(n_sheet, n_site)
      lat = lattice
      npts = t_n
   end subroutine tab_info_out

   !------------------------------------------------------------------
   ! akima_build(n, z, e, m) - the Akima (1970) local slopes: secant
   !                 slopes d(i) between points i/i+1, edge-extended by
   !                 the original two-point linear extrapolation, then
   !                 m(i) = (|d(i+1)-d(i)|*d(i-1) + |d(i-1)-d(i-2)|*d(i))
   !                 / (sum); flat tie -> the secant average
   !------------------------------------------------------------------
   subroutine akima_build(n, z, e, m)
      integer, intent(in) :: n
      real(8), intent(in) :: z(n), e(n)
      real(8), intent(out) :: m(n)
      real(8) :: d(-1:n + 1)
      real(8) :: w1, w2
      integer :: i
      do i = 1, n - 1
         d(i) = (e(i + 1) - e(i))/(z(i + 1) - z(i))
      end do
      d(0) = 2.0d0*d(1) - d(2)
      d(-1) = 2.0d0*d(0) - d(1)
      d(n) = 2.0d0*d(n - 1) - d(n - 2)
      d(n + 1) = 2.0d0*d(n) - d(n - 1)
      do i = 1, n
         w1 = abs(d(i + 1) - d(i))
         w2 = abs(d(i - 1) - d(i - 2))
         if (w1 + w2 > 0.0d0) then
            m(i) = (w1*d(i - 1) + w2*d(i))/(w1 + w2)
         else
            m(i) = 0.5d0*(d(i - 1) + d(i))
         end if
      end do
   end subroutine akima_build

   integer function sheet_index(name)
      character(len=*), intent(in) :: name
      integer :: i
      sheet_index = 0
      do i = 1, n_sheet
         if (trim(name) == trim(sh_name(i))) sheet_index = i
      end do
   end function sheet_index

   integer function site_index(name)
      character(len=*), intent(in) :: name
      integer :: i
      site_index = 0
      do i = 1, n_site
         if (trim(name) == trim(si_name(i))) site_index = i
      end do
   end function site_index

   ! stop_tab(who, msg) - the named-abort channel (stderr + STOP 1)
   subroutine stop_tab(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_tab_tab error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container table layer is not usable)'
      stop 1
   end subroutine stop_tab
end module cau_tab_tab
