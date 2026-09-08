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
      ! Pseudocode:
      ! 1. Bin geometry: d = (v_hi - v_lo)/n_bin; bin centers c_i = v_lo + (i-0.5)*d
      !    (old :48-53/:175-180); guards: n_bin>1000 truncated, merge_at<1 raised to 1
      !    (old :167-172/:265-273)
      ! 2. Dispatch per spec%kind:
      ! 2a. k_e_hist / k_t_hist (old :50-67/:115-130): take data%e_rel or data%thta,
      !     bin into [v_lo,v_hi) with idx = int((v-v_lo)/d)+1 (out-of-range not counted,
      !     old :56-60); write (c_i, count) per bin; with norm_p a final
      !     'TOTAL & P' line (sum(count), sum(count)/n_traj - old :66)
      ! 2b. k_pos (old :87-92): per trajectory fold the drop point into the unit cell (as the
      !     first step of 2c), then write
      !     the (x, y, e_rel, z_min) four columns (old write(98,'(4F12.6)') :90)
      ! 2c. k_site (old :185-229): per trajectory first fold into the reference cell -
      !     lattice vectors
      !     e_x=(a_lat,0), e_y=(a_lat*cos(skew), a_lat*sin(skew)), fractional coordinates
      !     a1/a2 rounded and folded back (old Dist_cell :393-403); then nearest-site
      !     assignment - minloc over
      !     the cell-corner / bridge-edge-midpoint / face-center candidate distances, mapped to
      !     group names (old surfcell
      !     :340-372; candidate sets as data in C8); per site group, bin and write a
      !     derived-name file;
      !     with norm_p a TOTAL & P line per group (divided by that group's trajectory count -
      !     old :219-221)
      ! 2d. k_bnc (old :282-309): n_bnc<0 skipped (old :284); n_bnc >= merge_at
      !     folded into the last bin (old :285-289); write a histogram file per bounce class
      !     (class names
      !     _bounce_{0..merge_at-1} and _bounce_ge{merge_at} - old :298-309)
      ! 3. Output unit open(file_out, status='replace') ... close (old :62-67 sites);
      !    the HISTO directory layout belongs to the statistics-phase driver (old VENUS.f90:293)
   end subroutine plot_hist
end module hist
