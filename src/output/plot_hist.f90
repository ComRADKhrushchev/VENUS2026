!=====================================================================
! analysis/plot_hist.f90 - L4 analysis: five histogram/impact-map kinds (normalization
!                           parameterized; lattice constants
!                           all passed as arguments - no system constants)
! (module name hist: Fortran forbids a module sharing its contained routine's name; the
!  routine name plot_hist is the one fixed by the plan
!  §1 batch-4 table - the module yields the name to the routine)
! Reference: src_VENUS/HIST.f90:26-68 (plot_Erel_hist: bin centers lo+(i-1/2)d :48-53,
!       [lo,hi) binning int((v-lo)/d)+1 :56-60, TOTAL & P normalization line :66),
!       :70-93 (plot_Erel_impact_pos: cell-folded drop points + (x,y,e_rel,z_min) four
!       columns written :87-92),
!       :95-131 (plot_theta_hist - same-form angle version), :134-230
!       (plot_site_histograms: site-grouped double histograms - top/brg/hcp/fcc four groups,
!       per-group TOTAL & P :219-221, n_bin<=1000 guard :167-172), :232-311
!       (plot_Erel_by_bounce: bounce-count grouping + merge_at merge threshold :265-274,
!       >=merge_at folded into the last bin :283-296, n_bnc<0 skipped :284), :313-373 (surfcell:
!       Wigner-Seitz-type nearest-site assignment - minloc over the four cell corners / five
!       bridge-edge midpoints / two face-center candidate points
!       :340-372), :375-404 (Dist_cell: unit-cell PBC folding - fractional lattice-vector
!       coordinates rounded
!       and folded back :393-403);
!       src_VENUS/VENUS.f90:286-299 (the five call sites + HISTO directory layout - to the
!       statistics-phase driver)
!       (see previous line)
! Excision warning (hardwired system knowledge - the new version passes everything as
! arguments; K2 plan sourcing note):
!   (1) auro and PI/3 lattice constants hardwired (HIST.f90:19, :89, :188 - Au(111) system)
!       -> passed via lat_t (same source as the roster unit-cell fields - a_lat/skew);
!   (2) site names 'top'/'brg'/'hcp'/'fcc' and triangular-lattice corner/bridge/center
!       geometry hardwired (:165,
!       :340-345) -> site candidate sets as data (finalized in C8; the skeleton keeps the four
!       names for cross-checking)
!=====================================================================
module hist
   implicit none
   private
   public :: hist_spec_t, hist_data_t, lat_t, plot_hist
   public :: k_e_hist, k_t_hist, k_pos, k_site, k_bnc

   ! Five spec keys (values of spec%kind - the parameterized unification of the old five
   ! plot_* functions)
   integer, parameter :: k_e_hist = 1  ! kind 1: energy histogram (old plot_Erel_hist)
   integer, parameter :: k_t_hist = 2  ! kind 2: angle histogram (old plot_theta_hist)
   integer, parameter :: k_pos    = 3  ! kind 3: impact-position map (old plot_Erel_impact_pos)
   integer, parameter :: k_site   = 4  ! kind 4: site-grouped histograms (old plot_site_histograms)
   integer, parameter :: k_bnc    = 5  ! kind 5: bounce-count-grouped histograms (old plot_Erel_by_bounce)

   ! Histogram specification (range / bin count / normalization - all input data)
   type :: hist_spec_t
      integer :: kind = 0                  ! five-kind key (k_*)
      character(len=64) :: file_out = ''   ! output file name (old 'HISTO/...' - the directory
                                           ! layout belongs to the driver; grouped kinds derive
                                           ! names as
      ! file_out prefix + group name)
      real(8) :: v_lo = 0.0d0              ! range lower bound [units of the quantity] (old Erel_min/theta_min)
      real(8) :: v_hi = 0.0d0              ! range upper bound (old Erel_max/theta_max)
      integer :: n_bin = 0                 ! bin count [count] (<=1000 guard, old :167-172/:265-268)
      integer :: merge_at = 0              ! kind-5 merge threshold [count] (bounce counts >= this fold into the last bin)
      logical :: norm_p = .false.          ! final TOTAL & P normalization line (old :66/:220)
   end type hist_spec_t

   ! Trajectory row table (aggregation of recorder final-state rows - the old AErelS/AInitX
   ! etc. global arrays turned into members)
   type :: hist_data_t
      integer :: n_traj = 0                    ! trajectory row count [row] (old NT)
      real(8), allocatable :: e_rel(:)         ! final relative energy column [kcal/mol] (old AErelS)
      real(8), allocatable :: thta(:)          ! scattering-angle column [deg] (old AThtaS = ANG(19))
      real(8), allocatable :: x0(:)            ! initial drop-point x column [Å] (old AInitX = RX0)
      real(8), allocatable :: y0(:)            ! initial drop-point y column [Å] (old AInitY = RY0)
      real(8), allocatable :: z_min(:)         ! minimum-height column [Å] (old AZminS)
      real(8), allocatable :: v_fin(:)         ! final potential-energy column [kcal/mol] (old AVfinS)
      integer, allocatable :: n_bnc(:)         ! bounce-count column [count] (old NBncS; -1 = invalid placeholder)
      character(len=3), allocatable :: site(:) ! site column [-] (old site - record_all :19)
   end type hist_data_t

   ! Lattice constants (all passed as arguments - no system constants; same source as the
   ! roster unit-cell fields)
   type :: lat_t
      real(8) :: a_lat = 0.0d0     ! lattice constant [Å] (old hardwired auro - excision warning (1))
      real(8) :: skew = 0.0d0      ! unit-cell skew angle [rad] (old hardwired PI/3 - excision warning (1))
      real(8) :: q_ref(2) = 0.0d0  ! reference site coordinates [Å] (old Qref - unit-cell
                                   ! origin;
                                   ! the old call site hardwired Qr2 = 2.55476, VENUS.f90:290)
   end type lat_t
contains
   !------------------------------------------------------------------
   ! plot_hist(spec, data, lattice) - dispatch of the five histogram/impact-map kinds
   ! (normalization parameterized)
   ! Reference: HIST.f90:26-311 (the five functions), :313-404 (the surfcell/Dist_cell tools);
   !       until the site candidate sets are data-driven, kind 4 / the site column runs with
   !       the four names for cross-checking - finalized in C8
   subroutine plot_hist(spec, data, lattice)
      type(hist_spec_t), intent(in) :: spec      ! histogram specification (kind key / range / bin count / normalization)
      type(hist_data_t), intent(in) :: data      ! trajectory row table (recorder aggregation)
      type(lat_t), intent(in) :: lattice         ! lattice constants (arguments - no system constants)
      integer :: nb, mrg                         ! guarded bin count / merge threshold
      if (spec%kind < k_e_hist .or. spec%kind > k_bnc) then
         write (0, '(a,i0)') 'plot_hist: unknown histogram kind ', spec%kind
         write (0, '(a)') 'plot_hist: loud-stop (the kind key is outside the five-kind set)'
         stop 1
      end if
      if (spec%kind /= k_pos .and. spec%n_bin < 1) then
         write (0, '(a,i0)') 'plot_hist: nonpositive bin count ', spec%n_bin
         write (0, '(a)') 'plot_hist: loud-stop (the bin geometry is not computable)'
         stop 1
      end if
      nb = spec%n_bin
      if (nb > 1000) nb = 1000                   ! the fixed-count array bound
      mrg = spec%merge_at
      if (mrg < 1) mrg = 1
      select case (spec%kind)
      case (k_e_hist)
         call column_hist(spec%file_out, data%e_rel, data%n_traj, spec, nb)
      case (k_t_hist)
         call column_hist(spec%file_out, data%thta, data%n_traj, spec, nb)
      case (k_pos)
         call pos_map(spec%file_out, data, lattice)
      case (k_site)
         call site_hists(spec, data, lattice, nb)
      case (k_bnc)
         call bounce_hists(spec, data, nb, mrg)
      end select
   end subroutine plot_hist

   !------------------------------------------------------------------
   ! column_hist(fname, vals, n_traj, spec, nb) - one plain (center, count)
   !              histogram file, optionally closed by the TOTAL & P line
   !              (the two plain kinds share everything but their column)
   subroutine column_hist(fname, vals, n_traj, spec, nb)
      character(len=*), intent(in) :: fname      ! output file name (caller-staged)
      real(8), intent(in) :: vals(:)             ! the binned column
      integer, intent(in) :: n_traj              ! trajectory row count (the P denominator)
      type(hist_spec_t), intent(in) :: spec      ! histogram specification
      integer, intent(in) :: nb                  ! guarded bin count
      real(8) :: dlt, cnt(1000)
      integer :: u, i, bi
      dlt = (spec%v_hi - spec%v_lo)/dble(nb)
      cnt = 0.0d0
      do i = 1, n_traj
         if (vals(i) >= spec%v_lo .and. vals(i) < spec%v_hi) then
            bi = int((vals(i) - spec%v_lo)/dlt) + 1
            if (bi >= 1 .and. bi <= nb) cnt(bi) = cnt(bi) + 1.0d0
         end if
      end do
      open (newunit=u, file=fname, status='replace', action='write')
      do i = 1, nb
         write (u, '(2f12.6)') spec%v_lo + (dble(i) - 0.5d0)*dlt, cnt(i)
      end do
      if (spec%norm_p) then
         write (u, '(a9,f14.8,f12.6)') 'TOTAL & P:', sum(cnt(1:nb)), &
                                       sum(cnt(1:nb))/dble(max(n_traj, 1))
      end if
      close (u)
   end subroutine column_hist

   !------------------------------------------------------------------
   ! pos_map(fname, data, lattice) - the impact-position map: every row's
   !              drop point folded into the reference cell, four columns
   !              (x, y, e_rel, z_min)
   subroutine pos_map(fname, data, lattice)
      character(len=*), intent(in) :: fname      ! output file name (caller-staged)
      type(hist_data_t), intent(in) :: data      ! trajectory row table
      type(lat_t), intent(in) :: lattice         ! lattice constants
      real(8) :: f(2)
      integer :: u, i
      open (newunit=u, file=fname, status='replace', action='write')
      do i = 1, data%n_traj
         call cell_fold([data%x0(i), data%y0(i)], lattice, f)
         write (u, '(4f12.6)') f(1), f(2), data%e_rel(i), data%z_min(i)
      end do
      close (u)
   end subroutine pos_map

   !------------------------------------------------------------------
   ! site_hists(spec, data, lattice, nb) - site-grouped histograms: every
   !              row folded into the reference cell, assigned its nearest
   !              site candidate (cell corners = top / bridge-edge midpoints
   !              = brg / the two face centroids = hcp, fcc - candidates
   !              constructed from a_lat and skew, the group labels the
   !              four conventional names), then one derived-name histogram
   !              file per group with the GROUP's row count as the P
   !              denominator
   subroutine site_hists(spec, data, lattice, nb)
      type(hist_spec_t), intent(in) :: spec      ! histogram specification (file_out = the name prefix)
      type(hist_data_t), intent(in) :: data      ! trajectory row table
      type(lat_t), intent(in) :: lattice         ! lattice constants
      integer, intent(in) :: nb                  ! guarded bin count
      character(len=3), parameter :: gname(4) = ['top', 'brg', 'hcp', 'fcc']
      real(8) :: cell(2, 4), edge(2, 5), cand(2, 11)
      real(8) :: dist(11), f(2), dlt, cnt(4, 1000)
      integer :: count_all(4), mi, g, u, i, bi
      dlt = (spec%v_hi - spec%v_lo)/dble(nb)
      ! the eleven site candidates of the unit cell (relative to the cell origin)
      cell(:, 1) = [0.0d0, 0.0d0]
      cell(:, 2) = [lattice%a_lat, 0.0d0]
      cell(:, 3) = [lattice%a_lat*cos(lattice%skew), lattice%a_lat*sin(lattice%skew)]
      cell(:, 4) = [lattice%a_lat*(1.0d0 + cos(lattice%skew)), lattice%a_lat*sin(lattice%skew)]
      do i = 1, 4
         edge(:, i) = (cell(:, i) + cell(:, mod(i, 4) + 1))/2.0d0
      end do
      edge(:, 5) = (cell(:, 1) + cell(:, 3))/2.0d0
      cand(:, 1:4) = cell                     ! top
      cand(:, 5:9) = edge                     ! brg
      cand(:, 10) = (cell(:, 1) + cell(:, 3) + cell(:, 2))/3.0d0   ! hcp
      cand(:, 11) = (cell(:, 2) + cell(:, 4) + cell(:, 3))/3.0d0   ! fcc
      cnt = 0.0d0
      count_all = 0
      do i = 1, data%n_traj
         call cell_fold([data%x0(i), data%y0(i)], lattice, f)
         do mi = 1, 11
            dist(mi) = norm2(f - cand(:, mi))
         end do
         mi = minloc(dist, 1)
         if (mi <= 4) then
            g = 1
         else if (mi <= 9) then
            g = 2
         else if (mi == 10) then
            g = 3
         else
            g = 4
         end if
         count_all(g) = count_all(g) + 1
         if (data%e_rel(i) >= spec%v_lo .and. data%e_rel(i) < spec%v_hi) then
            bi = int((data%e_rel(i) - spec%v_lo)/dlt) + 1
            if (bi >= 1 .and. bi <= nb) cnt(g, bi) = cnt(g, bi) + 1.0d0
         end if
      end do
      do g = 1, 4
         open (newunit=u, file=trim(spec%file_out)//'_'//gname(g)//'.dat', status='replace', &
               action='write')
         do i = 1, nb
            write (u, '(2f12.6)') spec%v_lo + (dble(i) - 0.5d0)*dlt, cnt(g, i)
         end do
         if (spec%norm_p) then
            write (u, '(a9,f14.8,f12.6)') 'TOTAL & P:', sum(cnt(g, 1:nb)), &
                                          sum(cnt(g, 1:nb))/dble(max(count_all(g), 1))
         end if
         close (u)
      end do
   end subroutine site_hists

   !------------------------------------------------------------------
   ! bounce_hists(spec, data, nb, mrg) - bounce-count-grouped histograms:
   !              invalid rows skipped, bounce counts at/above the merge
   !              threshold folded into the ge class, one derived-name
   !              histogram file per class
   subroutine bounce_hists(spec, data, nb, mrg)
      type(hist_spec_t), intent(in) :: spec      ! histogram specification (file_out = the name prefix)
      type(hist_data_t), intent(in) :: data      ! trajectory row table
      integer, intent(in) :: nb, mrg             ! guarded bin count / merge threshold
      real(8), allocatable :: cnt(:, :)
      real(8) :: dlt
      character(len=256) :: fname
      integer :: u, i, bi, ci
      dlt = (spec%v_hi - spec%v_lo)/dble(nb)
      allocate (cnt(mrg + 1, nb))
      cnt = 0.0d0
      do i = 1, data%n_traj
         if (data%n_bnc(i) < 0) cycle
         if (data%n_bnc(i) >= mrg) then
            ci = mrg + 1
         else
            ci = data%n_bnc(i) + 1
         end if
         if (data%e_rel(i) >= spec%v_lo .and. data%e_rel(i) < spec%v_hi) then
            bi = int((data%e_rel(i) - spec%v_lo)/dlt) + 1
            if (bi >= 1 .and. bi <= nb) cnt(ci, bi) = cnt(ci, bi) + 1.0d0
         end if
      end do
      do ci = 1, mrg + 1
         if (ci < mrg + 1) then
            write (fname, '(a,i0,a)') trim(spec%file_out)//'_bounce_', ci - 1, '.dat'
         else
            write (fname, '(a,i0,a)') trim(spec%file_out)//'_bounce_ge', mrg, '.dat'
         end if
         open (newunit=u, file=trim(fname), status='replace', action='write')
         do i = 1, nb
            write (u, '(2f12.6)') spec%v_lo + (dble(i) - 0.5d0)*dlt, cnt(ci, i)
         end do
         close (u)
      end do
      deallocate (cnt)
   end subroutine bounce_hists

   !------------------------------------------------------------------
   ! cell_fold(q, lattice, r) - fold one drop point into the reference cell:
   !              fractional lattice coordinates truncated to the enclosing
   !              cell, returning the Cartesian offset of q within it
   !              (relative to the cell origin at q_ref; the folded remainder
   !              lies in [0,1] - an exactly integer non-positive coordinate
   !              folds to the upper boundary 1.0 under the -1 adjustment)
   subroutine cell_fold(q, lattice, r)
      real(8), intent(in) :: q(2)                ! the drop point [A]
      type(lat_t), intent(in) :: lattice         ! lattice constants
      real(8), intent(out) :: r(2)               ! the folded offset within the cell [A]
      real(8) :: a1, a2, q_op(2)
      integer :: n1, n2
      a2 = (q(2) - lattice%q_ref(2))/(lattice%a_lat*sin(lattice%skew))
      a1 = (q(1) - lattice%q_ref(1) - a2*lattice%a_lat*cos(lattice%skew))/lattice%a_lat
      if (a1 > 0.0d0) then
         n1 = int(a1)
      else
         n1 = int(a1) - 1
      end if
      if (a2 > 0.0d0) then
         n2 = int(a2)
      else
         n2 = int(a2) - 1
      end if
      q_op = n1*[lattice%a_lat, 0.0d0] &
           + n2*[lattice%a_lat*cos(lattice%skew), lattice%a_lat*sin(lattice%skew)] &
           + lattice%q_ref
      r = q - q_op
   end subroutine cell_fold
end module hist
