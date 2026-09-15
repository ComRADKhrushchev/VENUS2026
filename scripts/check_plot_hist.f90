!=====================================================================
! check_plot_hist.f90 - criteria program for the five histogram/impact-map
!                        kinds of plot_hist (binning, the normalization
!                        parameterization, the five-class dispatch, the
!                        lattice channel; criteria defined FIRST in this
!                        header; this program executes them)
! Compile (from src/; -Jbuild keeps the .mod files out of the source tree):
!   gfortran -Jbuild output/plot_hist.f90 ../scripts/check_plot_hist.f90 -o build/check_plot_hist.exe
! Run (from src/): build\check_plot_hist.exe
!
! Constructed material (the anchor derivations):
!   The reference lattice throughout: a_lat = 2.0 A, skew = pi/3 (the
!   oblique cell ex = (2,0), ey = (1,sqrt(3))), q_ref = (0.5, 0.25).
!   Drop points are built as q_ref + (f1 + n1)*ex + (f2 + n2)*ey with
!   interior fractional parts f (away from the fold boundaries by >=
!   0.05) and integer lattice shifts n (the folding is exercised); the
!   expected folded position is f1*ex + f2*ey.
!   Site anchor points (fractional, in the reference cell):
!     top  (0.1 ,0.1 ) - nearest the cell corner c1;
!     brg  (0.5 ,0.05) - nearest the c1-c2 bridge midpoint;
!     hcp  (1/3 ,1/3 ) - exactly the c1-c3-c2 face centroid;
!     fcc  (2/3 ,2/3 ) - exactly the c2-c4-c3 face centroid
!   (both centroids sit at interior fractional positions for every skew -
!   f = (1/3,1/3) and (2/3,2/3) are skew-independent). Margins to the
!   nearest other candidate exceed 0.15 A in every case.
!   Histogram values: bin centers and counts are exact-binary and <= 6
!   decimals, so the (2F12.6) print round-trips them exactly; the
!   Cartesian folded coordinates carry sqrt(3) and are compared at the
!   print-quantization scale (see Tolerances).
!
! Tolerances:
!   - histogram centers/counts: abs 1e-9 (exact-binary anchors; a wrong
!     bin edge or index deviates by O(1));
!   - norm-line TOTAL and P: abs 2e-6 (F14.8/F12.6 print quantization is
!     5e-7/5e-7; a wrong denominator - trajectory count vs group count -
!     deviates at the 1e-2 level);
!   - impact-map Cartesian columns: abs 2e-6 (sqrt(3) values quantize at
!     5e-7 on print; a wrong fold or lattice shifts a column by >= half
!     a lattice constant ~ 1 A).
!
! Criteria:
!   (a) uniform-density histogram (kind 1): two samples at each of the
!       10 bin centers of [0,10) plus four out-of-range rows (below the
!       range, exactly at v_hi, above the range) -> every bin count 2,
!       centers on v_lo + (i-1/2)*d, the out-of-range rows not counted
!       (the v == v_hi edge excluded - the [v_lo,v_hi) contract)
!   (b) the normalization parameterization: a ramp data set (bin i
!       carries i samples; 55 in-range of 57 rows) plotted with
!       norm_p off and on -> the bin/center lines IDENTICAL, the
!       normalization adds exactly one trailing line whose TOTAL is the
!       in-range count (the histogram area) and whose P is that count
!       over the trajectory count (the fraction); the peak bin is
!       untouched by the switch
!   (c) one minimal anchor per kind: kind 1 = criterion (a); kind 2 =
!       angle histogram ([0,180), 4 bins, one sample per bin, one
!       out-of-range angle, norm on: TOTAL 4, P 4/5); kind 3 =
!       impact-position map (three folded rows, e_rel and z_min passed
!       through verbatim); kind 4 = site-grouped histograms (one row per
!       site group at the anchor fractional points plus a second top row
!       with an out-of-range energy -> per-group files, the top group's
!       P = 1/2 anchors the per-group denominator, the other groups'
!       P = 1); kind 5 = bounce-grouped histograms (merge_at = 3: classes
!       0,1,2,ge3; an n_bnc = -1 row skipped; n_bnc = 5 folded into the
!       ge class)
!   (d) the lattice channel: (d1) the same absolute drop point folded
!       with q_ref shifted by +0.1 A in x -> the map column shifts by
!       exactly that (the fractional remainder moves); (d2) the same top
!       point assigned with q_ref shifted by half a lattice vector ->
!       the assignment flips to the bridge group (the empty top file and
!       the populated bridge file are both checked)
! Registered notes (part of this check's contract):
!   - kind 2 honors the normalization line (the normalization is a
!     parameter of the unified dispatch, not a property of one kind);
!   - the site grouping divides by the GROUP's trajectory count (a
!     whole-table denominator would give 1/5 where 1/2 is anchored).
!=====================================================================
program check_plot_hist
   use hist, only: hist_spec_t, hist_data_t, lat_t, plot_hist, &
                   k_e_hist, k_t_hist, k_pos, k_site, k_bnc
   implicit none
   integer :: nfail = 0
   integer :: nfail0                     ! per-criterion failure watermark
   real(8), parameter :: pi = 3.141592653589793d0

   call crit_a()      ! (a) uniform-density binning (kind 1)
   call crit_b()      ! (b) the normalization parameterization
   call crit_c()      ! (c) one minimal anchor per kind (kinds 2..5)
   call crit_d()      ! (d) the lattice channel

   if (nfail == 0) then
      print '(a)', 'check_plot_hist: ALL PASS'
   else
      print '(a,i0,a)', 'check_plot_hist: ', nfail, ' FAILURE(S)'
      stop 1
   end if
contains
   !------------------------------------------------------------------
   ! crit_a() - (a): the uniform-density histogram (kind 1)
   subroutine crit_a()
      type(hist_data_t) :: d
      type(hist_spec_t) :: sp
      real(8) :: want(10)
      integer :: i, k
      nfail0 = nfail
      call rows_alloc(d, 24)
      k = 1
      do i = 1, 10                        ! two samples per bin, at the centers
         d%e_rel(k) = dble(i) - 0.5d0
         k = k + 1
         d%e_rel(k) = dble(i) - 0.5d0
         k = k + 1
      end do
      d%e_rel(21:24) = [-0.5d0, 10.0d0, 12.75d0, -3.25d0]   ! below / at v_hi / above
      sp = hist_spec_t(k_e_hist, 'build/ph_uniform.dat', 0.0d0, 10.0d0, 10, 0, .false.)
      call plot_hist(sp, d, lat_t())
      want = 2.0d0
      call expect_hist('(a) uniform kind 1', 'build/ph_uniform.dat', 10, 0.0d0, 1.0d0, &
                       want, .false., 0.0d0, 0.0d0)
      if (nfail == nfail0) then
         print '(a)', '(a): PASS (ten bins at count 2; out-of-range rows and the v_hi'// &
            ' edge not counted; centers v_lo + (i-1/2)*d)'
      end if
   end subroutine crit_a

   !------------------------------------------------------------------
   ! crit_b() - (b): the normalization parameterization (area/fraction
   !             semantics, counts untouched by the switch)
   subroutine crit_b()
      type(hist_data_t) :: d
      type(hist_spec_t) :: sp
      real(8) :: want(10)
      integer :: i, j, k
      nfail0 = nfail
      call rows_alloc(d, 57)
      k = 1
      do i = 1, 10                        ! the ramp: bin i carries i samples
         do j = 1, i
            d%e_rel(k) = dble(i) - 0.5d0
            k = k + 1
         end do
      end do
      d%e_rel(56:57) = [11.5d0, -0.25d0]  ! two out-of-range rows
      sp = hist_spec_t(k_e_hist, 'build/ph_ramp_off.dat', 0.0d0, 10.0d0, 10, 0, .false.)
      call plot_hist(sp, d, lat_t())
      sp%file_out = 'build/ph_ramp_on.dat'
      sp%norm_p = .true.
      call plot_hist(sp, d, lat_t())
      want = [(dble(i), i=1, 10)]
      call expect_hist('(b) ramp, norm off', 'build/ph_ramp_off.dat', 10, 0.0d0, 1.0d0, &
                       want, .false., 0.0d0, 0.0d0)
      call expect_hist('(b) ramp, norm on', 'build/ph_ramp_on.dat', 10, 0.0d0, 1.0d0, &
                       want, .true., 55.0d0, 55.0d0/57.0d0)
      if (nfail == nfail0) then
         print '(a)', '(b): PASS (bin/center lines identical in both states; the'// &
            ' normalization adds one trailing line: area TOTAL = 55 in-range rows,'// &
            ' fraction P = 55/57; the peak bin untouched)'
      end if
   end subroutine crit_b

   !------------------------------------------------------------------
   ! crit_c() - (c): one minimal anchor per kind (2..5; kind 1 = crit_a)
   subroutine crit_c()
      type(hist_data_t) :: d
      type(hist_spec_t) :: sp
      type(lat_t) :: lat
      real(8) :: want(10), rows(4, 3)
      real(8), parameter :: s3 = 1.7320508075688772d0
      integer :: i
      logical :: ok
      nfail0 = nfail
      lat = lat_t(a_lat=2.0d0, skew=pi/3.0d0, q_ref=[0.5d0, 0.25d0])
      ! ---- kind 2: the angle histogram ----
      call rows_alloc(d, 5)
      d%thta = [22.5d0, 67.5d0, 112.5d0, 157.5d0, 200.0d0]
      sp = hist_spec_t(k_t_hist, 'build/ph_theta.dat', 0.0d0, 180.0d0, 4, 0, .true.)
      call plot_hist(sp, d, lat_t())
      call expect_hist('(c) kind 2 angle histogram', 'build/ph_theta.dat', 4, 0.0d0, 45.0d0, &
                       [1.0d0, 1.0d0, 1.0d0, 1.0d0], .true., 4.0d0, 0.8d0)
      ! ---- kind 3: the impact-position map ----
      call rows_alloc(d, 3)
      call set_point(d, 1, lat, [0.3d0, 0.2d0], [3, -2])
      call set_point(d, 2, lat, [0.6d0, 0.6d0], [0, 0])
      call set_point(d, 3, lat, [0.05d0, 0.5d0], [-4, 2])
      d%e_rel = [1.5d0, 9.25d0, 0.0625d0]
      d%z_min = [-0.25d0, 0.125d0, 3.0d0]
      sp = hist_spec_t(k_pos, 'build/ph_pos.dat', 0.0d0, 0.0d0, 0, 0, .false.)
      call plot_hist(sp, d, lat)
      call read_pos('build/ph_pos.dat', 3, rows, ok)
      if (.not. ok) then
         call fail('(c) kind 3 impact map: unreadable or mis-shaped file')
      else
         ! the folded Cartesian anchors: f1*ex + f2*ey with ex=(2,0), ey=(1,sqrt(3))
         call chk6('(c) kind 3 row 1 folded x', rows(1, 1), 0.3d0*2.0d0 + 0.2d0)
         call chk6('(c) kind 3 row 1 folded y', rows(2, 1), 0.2d0*s3)
         call chk6('(c) kind 3 row 2 folded x', rows(1, 2), 0.6d0*2.0d0 + 0.6d0)
         call chk6('(c) kind 3 row 2 folded y', rows(2, 2), 0.6d0*s3)
         call chk6('(c) kind 3 row 3 folded x', rows(1, 3), 0.05d0*2.0d0 + 0.5d0)
         call chk6('(c) kind 3 row 3 folded y', rows(2, 3), 0.5d0*s3)
         call chk6('(c) kind 3 row 1 e_rel passthrough', rows(3, 1), 1.5d0)
         call chk6('(c) kind 3 row 2 e_rel passthrough', rows(3, 2), 9.25d0)
         call chk6('(c) kind 3 row 3 e_rel passthrough', rows(3, 3), 0.0625d0)
         call chk6('(c) kind 3 row 1 z_min passthrough', rows(4, 1), -0.25d0)
         call chk6('(c) kind 3 row 2 z_min passthrough', rows(4, 2), 0.125d0)
         call chk6('(c) kind 3 row 3 z_min passthrough', rows(4, 3), 3.0d0)
      end if
      ! ---- kind 4: the site-grouped histograms ----
      call rows_alloc(d, 5)
      call set_point(d, 1, lat, [0.1d0, 0.1d0], [2, 1])        ! top
      call set_point(d, 2, lat, [0.05d0, 0.05d0], [0, 0])      ! top (out-of-range energy)
      call set_point(d, 3, lat, [0.5d0, 0.05d0], [-1, 3])      ! brg
      call set_point(d, 4, lat, [1.0d0/3.0d0, 1.0d0/3.0d0], [5, -2])   ! hcp
      call set_point(d, 5, lat, [2.0d0/3.0d0, 2.0d0/3.0d0], [0, 0])    ! fcc
      d%e_rel = [0.5d0, 99.0d0, 3.5d0, 7.25d0, 9.75d0]
      sp = hist_spec_t(k_site, 'build/ph_site', 0.0d0, 10.0d0, 10, 0, .true.)
      call plot_hist(sp, d, lat)
      want = 0.0d0
      want(1) = 1.0d0
      call expect_hist('(c) kind 4 group top', 'build/ph_site_top.dat', 10, 0.0d0, 1.0d0, &
                       want, .true., 1.0d0, 0.5d0)
      want = 0.0d0
      want(4) = 1.0d0
      call expect_hist('(c) kind 4 group brg', 'build/ph_site_brg.dat', 10, 0.0d0, 1.0d0, &
                       want, .true., 1.0d0, 1.0d0)
      want = 0.0d0
      want(8) = 1.0d0
      call expect_hist('(c) kind 4 group hcp', 'build/ph_site_hcp.dat', 10, 0.0d0, 1.0d0, &
                       want, .true., 1.0d0, 1.0d0)
      want = 0.0d0
      want(10) = 1.0d0
      call expect_hist('(c) kind 4 group fcc', 'build/ph_site_fcc.dat', 10, 0.0d0, 1.0d0, &
                       want, .true., 1.0d0, 1.0d0)
      ! ---- kind 5: the bounce-grouped histograms ----
      call rows_alloc(d, 6)
      d%n_bnc = [0, 0, 1, 3, -1, 5]
      d%e_rel = [0.25d0, 0.5d0, 1.5d0, 3.25d0, 9.9d0, 7.5d0]
      sp = hist_spec_t(k_bnc, 'build/ph_bnc', 0.0d0, 10.0d0, 10, 3, .false.)
      call plot_hist(sp, d, lat_t())
      want = 0.0d0
      want(1) = 2.0d0
      call expect_hist('(c) kind 5 bounce class 0', 'build/ph_bnc_bounce_0.dat', 10, 0.0d0, &
                       1.0d0, want, .false., 0.0d0, 0.0d0)
      want = 0.0d0
      want(2) = 1.0d0
      call expect_hist('(c) kind 5 bounce class 1', 'build/ph_bnc_bounce_1.dat', 10, 0.0d0, &
                       1.0d0, want, .false., 0.0d0, 0.0d0)
      want = 0.0d0
      call expect_hist('(c) kind 5 bounce class 2 (empty)', 'build/ph_bnc_bounce_2.dat', 10, &
                       0.0d0, 1.0d0, want, .false., 0.0d0, 0.0d0)
      want = 0.0d0
      want(4) = 1.0d0
      want(8) = 1.0d0
      call expect_hist('(c) kind 5 bounce class ge3 (merged)', 'build/ph_bnc_bounce_ge3.dat', &
                       10, 0.0d0, 1.0d0, want, .false., 0.0d0, 0.0d0)
      if (nfail == nfail0) then
         print '(a)', '(c): PASS (angle histogram with its normalization line; the'// &
            ' impact map folded with e_rel/z_min verbatim; the four site groups with'// &
            ' the per-group denominator; the bounce classes with the -1 skip and the'// &
            ' ge-merge)'
      end if
   end subroutine crit_c

   !------------------------------------------------------------------
   ! crit_d() - (d): the lattice channel
   subroutine crit_d()
      type(hist_data_t) :: d
      type(hist_spec_t) :: sp
      type(lat_t) :: lat1, lat2
      real(8) :: rows(4, 1), want(10)
      real(8) :: pq(2)
      real(8), parameter :: s3 = 1.7320508075688772d0
      logical :: ok
      nfail0 = nfail
      lat1 = lat_t(a_lat=2.0d0, skew=pi/3.0d0, q_ref=[0.5d0, 0.25d0])
      ! ---- (d1): same drop point, q_ref shifted by +0.1 A in x ----
      lat2 = lat_t(a_lat=2.0d0, skew=pi/3.0d0, q_ref=[0.6d0, 0.25d0])
      pq = lat_point(lat1, [0.3d0, 0.2d0], [3, -2])
      call rows_alloc(d, 1)
      d%x0(1) = pq(1)
      d%y0(1) = pq(2)
      d%e_rel(1) = 1.5d0
      d%z_min(1) = -0.25d0
      sp = hist_spec_t(k_pos, 'build/ph_d1_ref.dat', 0.0d0, 0.0d0, 0, 0, .false.)
      call plot_hist(sp, d, lat1)
      sp%file_out = 'build/ph_d1_shift.dat'
      call plot_hist(sp, d, lat2)
      call read_pos('build/ph_d1_ref.dat', 1, rows, ok)
      if (.not. ok) then
         call fail('(d1) the reference-latitude map is unreadable')
      else
         call chk6('(d1) folded x under the reference cell', rows(1, 1), 0.3d0*2.0d0 + 0.2d0)
         call chk6('(d1) folded y under the reference cell', rows(2, 1), 0.2d0*s3)
      end if
      call read_pos('build/ph_d1_shift.dat', 1, rows, ok)
      if (.not. ok) then
         call fail('(d1) the shifted-lattice map is unreadable')
      else
         ! the shift moves the fractional remainder: f1 0.3 -> 0.25, so x -0.1
         call chk6('(d1) folded x under the shifted q_ref', rows(1, 1), 0.25d0*2.0d0 + 0.2d0)
         call chk6('(d1) folded y unchanged by the x shift', rows(2, 1), 0.2d0*s3)
      end if
      ! ---- (d2): same top point, q_ref shifted by half a lattice vector ----
      lat2 = lat_t(a_lat=2.0d0, skew=pi/3.0d0, q_ref=[1.5d0, 0.25d0])
      pq = lat_point(lat1, [0.1d0, 0.1d0], [0, 0])
      call rows_alloc(d, 1)
      d%x0(1) = pq(1)
      d%y0(1) = pq(2)
      d%e_rel(1) = 0.5d0
      sp = hist_spec_t(k_site, 'build/ph_d2', 0.0d0, 10.0d0, 10, 0, .true.)
      call plot_hist(sp, d, lat2)
      want = 0.0d0
      want(1) = 1.0d0
      call expect_hist('(d2) the flipped group (bridge)', 'build/ph_d2_brg.dat', 10, 0.0d0, &
                       1.0d0, want, .true., 1.0d0, 1.0d0)
      want = 0.0d0
      call expect_hist('(d2) the vacated group (top)', 'build/ph_d2_top.dat', 10, 0.0d0, &
                       1.0d0, want, .true., 0.0d0, 0.0d0)
      if (nfail == nfail0) then
         print '(a)', '(d): PASS (the q_ref shift moves the folded remainder; the'// &
            ' half-lattice shift flips the site assignment top -> bridge)'
      end if
   end subroutine crit_d

   !------------------------------------------------------------------
   ! construction helpers
   !------------------------------------------------------------------

   ! rows_alloc - a full row table of n trajectories, zeroed (any previous
   !              table dropped first - allocating an already-allocated
   !              array is a runtime error)
   subroutine rows_alloc(d, n)
      type(hist_data_t), intent(inout) :: d
      integer, intent(in) :: n
      if (allocated(d%e_rel)) deallocate (d%e_rel)
      if (allocated(d%thta)) deallocate (d%thta)
      if (allocated(d%x0)) deallocate (d%x0)
      if (allocated(d%y0)) deallocate (d%y0)
      if (allocated(d%z_min)) deallocate (d%z_min)
      if (allocated(d%v_fin)) deallocate (d%v_fin)
      if (allocated(d%n_bnc)) deallocate (d%n_bnc)
      if (allocated(d%site)) deallocate (d%site)
      allocate (d%e_rel(n), d%thta(n), d%x0(n), d%y0(n), d%z_min(n), d%v_fin(n), &
                d%n_bnc(n), d%site(n))
      d%n_traj = n
      d%e_rel = 0.0d0
      d%thta = 0.0d0
      d%x0 = 0.0d0
      d%y0 = 0.0d0
      d%z_min = 0.0d0
      d%v_fin = 0.0d0
      d%n_bnc = 0
      d%site = ''
   end subroutine rows_alloc

   ! lat_point - q_ref + (f + n) riding the lattice vectors (a drop-point
   !             builder: f = the interior fractional part, n = the shift)
   function lat_point(lat, f, nn) result(q)
      type(lat_t), intent(in) :: lat
      real(8), intent(in) :: f(2)
      integer, intent(in) :: nn(2)
      real(8) :: q(2)
      q = lat%q_ref + (f(1) + dble(nn(1)))*[lat%a_lat, 0.0d0] &
                    + (f(2) + dble(nn(2)))*[lat%a_lat*cos(lat%skew), &
                                             lat%a_lat*sin(lat%skew)]
   end function lat_point

   ! set_point - stage one row's drop point from fractional coordinates
   subroutine set_point(d, i, lat, f, nn)
      type(hist_data_t), intent(inout) :: d
      integer, intent(in) :: i
      type(lat_t), intent(in) :: lat
      real(8), intent(in) :: f(2)
      integer, intent(in) :: nn(2)
      real(8) :: q(2)
      q = lat_point(lat, f, nn)
      d%x0(i) = q(1)
      d%y0(i) = q(2)
   end subroutine set_point

   !------------------------------------------------------------------
   ! file readers / assertion helpers (watermark discipline)
   !------------------------------------------------------------------

   ! read_hist - one histogram file: nb (center,count) lines, optionally
   !             one TOTAL line, nothing more (ok = the shape holds)
   subroutine read_hist(fname, nb, centers, counts, norm, sum_c, p_val, ok)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: nb
      real(8), intent(out) :: centers(:), counts(:)
      logical, intent(out) :: norm
      real(8), intent(out) :: sum_c, p_val
      logical, intent(out) :: ok
      character(len=132) :: line
      integer :: u, ios, i
      norm = .false.
      sum_c = 0.0d0
      p_val = 0.0d0
      ok = .false.
      open (newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do i = 1, nb
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            close (u)
            return
         end if
         read (line, '(2f12.6)') centers(i), counts(i)
      end do
      read (u, '(a)', iostat=ios) line
      if (ios == 0) then
         if (index(line, 'TOTAL & P') == 0) then
            close (u)
            return
         end if
         norm = .true.
         read (line, '(9x,f14.8,f12.6)') sum_c, p_val
         read (u, '(a)', iostat=ios) line
         if (ios == 0) then
            close (u)
            return
         end if
      end if
      close (u)
      ok = .true.
   end subroutine read_hist

   ! expect_hist - a whole histogram file against its wants
   subroutine expect_hist(label, fname, nb, v_lo, dlt, want_counts, want_norm, &
                          want_sum, want_p)
      character(len=*), intent(in) :: label, fname
      integer, intent(in) :: nb
      real(8), intent(in) :: v_lo, dlt, want_counts(:)
      logical, intent(in) :: want_norm
      real(8), intent(in) :: want_sum, want_p
      real(8) :: centers(1000), counts(1000), s_c, p_v
      logical :: norm, ok
      integer :: i
      call read_hist(fname, nb, centers, counts, norm, s_c, p_v, ok)
      if (.not. ok) then
         call fail(trim(label)//': unreadable or mis-shaped file '//trim(fname))
         return
      end if
      do i = 1, nb
         call chk9(trim(label)//': bin center '//i2s(i), centers(i), &
                   v_lo + (dble(i) - 0.5d0)*dlt)
         call chk9(trim(label)//': bin count '//i2s(i), counts(i), want_counts(i))
      end do
      if (norm .neqv. want_norm) then
         call fail(trim(label)//': normalization line present = '//l2s(norm)//', want '// &
                   l2s(want_norm))
      end if
      if (want_norm) then
         call chk6(trim(label)//': TOTAL (the area)', s_c, want_sum)
         call chk6(trim(label)//': P (the fraction)', p_v, want_p)
      end if
   end subroutine expect_hist

   ! read_pos - one impact-map file: n four-column rows, nothing more
   subroutine read_pos(fname, n, rows, ok)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: n
      real(8), intent(out) :: rows(:, :)
      logical, intent(out) :: ok
      character(len=132) :: line
      integer :: u, ios, i
      ok = .false.
      open (newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do i = 1, n
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            close (u)
            return
         end if
         read (line, '(4f12.6)') rows(1, i), rows(2, i), rows(3, i), rows(4, i)
      end do
      read (u, '(a)', iostat=ios) line
      if (ios == 0) then
         close (u)
         return
      end if
      close (u)
      ok = .true.
   end subroutine read_pos

   ! chk9 - exact-binary anchor check (histogram fields; see Tolerances)
   subroutine chk9(label, got, want)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got, want
      if (abs(got - want) > 1.0d-9) then
         call fail(trim(label)//': got '//d2s(got)//', want '//d2s(want))
      end if
   end subroutine chk9

   ! chk6 - print-quantization-scale check (sqrt(3) columns / norm fields)
   subroutine chk6(label, got, want)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got, want
      if (abs(got - want) > 2.0d-6) then
         call fail(trim(label)//': got '//d2s(got)//', want '//d2s(want))
      end if
   end subroutine chk6

   function l2s(b)
      logical, intent(in) :: b
      character(len=5) :: l2s
      if (b) then
         l2s = 'true '
      else
         l2s = 'false'
      end if
   end function l2s

   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      print '(a)', 'FAIL '//trim(msg)
      nfail = nfail + 1
   end subroutine fail

   function i2s(i)
      integer, intent(in) :: i
      character(len=12) :: i2s
      write (i2s, '(i0)') i
   end function i2s
   function d2s(x)
      real(8), intent(in) :: x
      character(len=24) :: d2s
      write (d2s, '(es12.4)') x
   end function d2s
end program check_plot_hist
