!=====================================================================
! plot_hist.f90 - histogram kinds over trajectory final-state rows
! Design:
!   One dispatcher (plot_hist) unifies the kinds: the two plain column
!   histograms and the bounce-count-grouped histograms. Normalization is
!   parameterized; no system constants enter. The a_lat-based site
!   statistics (impact-position map, site-grouped histograms) retired with
!   the A_LAT key - 2026-10-13 ruling: such surface-site analysis belongs
!   to post-processing. Module named hist: Fortran forbids a module sharing
!   its name with a contained subroutine (the dispatcher is plot_hist).
!=====================================================================
module hist
   implicit none
   private
   public :: hist_spec_t, hist_data_t, plot_hist
   public :: k_e_hist, k_t_hist, k_bnc

   integer, parameter :: k_e_hist = 1  ! kind 1: energy histogram
   integer, parameter :: k_t_hist = 2  ! kind 2: angle histogram
   integer, parameter :: k_bnc    = 3  ! kind 3: bounce-count-grouped histograms

   ! Histogram specification (range / bin count / normalization - all input data)
   type :: hist_spec_t
      integer :: kind = 0                  ! kind key (k_*)
      character(len=64) :: file_out = ''   ! output file name (grouped kinds derive names as
                                           ! this prefix + group name)
      real(8) :: v_lo = 0.0d0              ! range lower bound [units of the binned quantity]
      real(8) :: v_hi = 0.0d0              ! range upper bound [units of the binned quantity]
      integer :: n_bin = 0                 ! bin count [count] (capped at 1000)
      integer :: merge_at = 0              ! bounce-kind merge threshold [count] (bounce counts >= this fold into the last bin)
      logical :: norm_p = .false.          ! write the final TOTAL & P normalization line
   end type hist_spec_t

   ! Trajectory row table (aggregation of recorder final-state rows)
   type :: hist_data_t
      integer :: n_traj = 0                    ! trajectory row count [row]
      real(8), allocatable :: e_rel(:)         ! final relative energy column [eV]
      real(8), allocatable :: thta(:)          ! scattering-angle column [deg]
      real(8), allocatable :: v_fin(:)         ! final potential-energy column [eV]
      integer, allocatable :: n_bnc(:)         ! bounce-count column [count] (-1 = invalid placeholder)
   end type hist_data_t
contains
   !------------------------------------------------------------------
   ! plot_hist(spec, data) - dispatch of the histogram kinds (normalization
   ! parameterized)
   !------------------------------------------------------------------
   subroutine plot_hist(spec, data)
      type(hist_spec_t), intent(in) :: spec      ! histogram specification (kind key / range / bin count / normalization)
      type(hist_data_t), intent(in) :: data      ! trajectory row table (recorder aggregation)
      integer :: nb, mrg                         ! guarded bin count / merge threshold
      ! 1. guards: known kind; positive bin count and nonempty range
      if (spec%kind < k_e_hist .or. spec%kind > k_bnc) then
         write (0, '(a,i0)') 'plot_hist: unknown histogram kind ', spec%kind
         write (0, '(a)') 'plot_hist: fatal (the kind key is outside the kind set)'
         stop 1
      end if
      if (spec%n_bin < 1) then
         write (0, '(a,i0)') 'plot_hist: nonpositive bin count ', spec%n_bin
         write (0, '(a)') 'plot_hist: fatal (the bin geometry is not computable)'
         stop 1
      end if
      if (spec%v_hi <= spec%v_lo) then
         write (0, '(2(a,es12.4))') 'plot_hist: empty histogram range v_hi <= v_lo (v_lo = ', &
            spec%v_lo, ', v_hi = ', spec%v_hi, ')'
         write (0, '(a)') 'plot_hist: fatal (the bin width would be nonpositive -'// &
            ' the range must satisfy v_hi > v_lo)'
         stop 1
      end if
      ! 2. clamp the bin count to the fixed-count array bound; floor the merge threshold
      nb = spec%n_bin
      if (nb > 1000) nb = 1000                   ! the fixed-count array bound
      mrg = spec%merge_at
      if (mrg < 1) mrg = 1
      ! 3. dispatch to the kind routine
      select case (spec%kind)
      case (k_e_hist)
         call column_hist(spec%file_out, data%e_rel, data%n_traj, spec, nb)
      case (k_t_hist)
         call column_hist(spec%file_out, data%thta, data%n_traj, spec, nb)
      case (k_bnc)
         call bounce_hists(spec, data, nb, mrg)
      end select
   end subroutine plot_hist

   !------------------------------------------------------------------
   ! column_hist(fname, vals, n_traj, spec, nb) - one plain (center, count)
   ! histogram file, optionally closed by the TOTAL & P line (the two plain
   ! kinds share everything but their column)
   !------------------------------------------------------------------
   subroutine column_hist(fname, vals, n_traj, spec, nb)
      character(len=*), intent(in) :: fname      ! output file name (caller-buffered)
      real(8), intent(in) :: vals(:)             ! the binned column
      integer, intent(in) :: n_traj              ! trajectory row count (the P denominator)
      type(hist_spec_t), intent(in) :: spec      ! histogram specification
      integer, intent(in) :: nb                  ! guarded bin count
      real(8) :: dlt, cnt(1000)
      integer :: u, i, bi
      dlt = (spec%v_hi - spec%v_lo)/dble(nb)
      cnt = 0.0d0
      ! 1. bin the rows ([lo,hi) half-open)
      do i = 1, n_traj
         if (vals(i) >= spec%v_lo .and. vals(i) < spec%v_hi) then
            bi = int((vals(i) - spec%v_lo)/dlt) + 1
            if (bi >= 1 .and. bi <= nb) cnt(bi) = cnt(bi) + 1.0d0
         end if
      end do
      open (newunit=u, file=fname, status='replace', action='write')
      ! 2. write bin centers and counts
      do i = 1, nb
         write (u, '(2f12.6)') spec%v_lo + (dble(i) - 0.5d0)*dlt, cnt(i)
      end do
      ! 3. optional TOTAL & P normalization line
      if (spec%norm_p) then
         write (u, '(a9,f14.8,f12.6)') 'TOTAL & P:', sum(cnt(1:nb)), &
                                       sum(cnt(1:nb))/dble(max(n_traj, 1))
      end if
      close (u)
   end subroutine column_hist

   !------------------------------------------------------------------
   ! bounce_hists(spec, data, nb, mrg) - bounce-count-grouped histograms:
   ! invalid rows skipped, bounce counts at/above the merge threshold folded
   ! into the ge class, one derived-name histogram file per class
   !------------------------------------------------------------------
   subroutine bounce_hists(spec, data, nb, mrg)
      type(hist_spec_t), intent(in) :: spec      ! histogram specification (file_out = the name prefix)
      type(hist_data_t), intent(in) :: data      ! trajectory row table
      integer, intent(in) :: nb, mrg             ! guarded bin count / merge threshold
      real(8), allocatable :: cnt(:, :)
      real(8) :: dlt
      character(len=256) :: fname
      integer :: u, i, bi, ci
      dlt = (spec%v_hi - spec%v_lo)/dble(nb)
      ! 1. mrg exact bounce classes + one merged >= mrg class
      allocate (cnt(mrg + 1, nb))
      cnt = 0.0d0
      ! 2. per-row class selection + [lo,hi) binning (invalid rows skipped)
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
      ! 3. one derived-name file per class
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
end module hist
