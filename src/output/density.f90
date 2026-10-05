!=====================================================================
! density.f90 - quantum density of states and formatted eigenvalue/vector output
! Design:
!   dens_states: Beyer-Swinehart direct counting of the cumulative state count
!   N(e); density taken by central difference at e +/- 50 bins. eig_out: writes
!   eigenvalues with their vectors in 6-column blocks (vectors temporarily
!   transposed to rows, restored after writing). Energies in unit bins (del = 1).
!=====================================================================
module density
   implicit none
   private
   public :: dens_states, eig_out
contains
   !------------------------------------------------------------------
   ! dens_states(w_mode, e_probe, rho) - density of states dN/dE at e_probe by
   ! central difference of the cumulative count at e +/- 50 bins
   !------------------------------------------------------------------
   subroutine dens_states(w_mode, e_probe, rho)
      real(8), intent(in) :: w_mode(:)  ! normal-mode frequency table [energy bins] (effective
                                        ! frequencies in units of the bin width del)
      real(8), intent(in) :: e_probe    ! probe energy [energy bins]
      real(8), intent(out) :: rho       ! density of states dN/dE [1/bin]
      real(8) :: s_hi, s_lo             ! the cumulative counts at the two probes
      call sum_count(w_mode, e_probe + 50.0d0, s_hi)
      call sum_count(w_mode, e_probe - 50.0d0, s_lo)
      rho = (s_hi - s_lo)/100.0d0
   end subroutine dens_states

   !------------------------------------------------------------------
   ! sum_count(w_mode, e, n_states) - cumulative quantum state count N(e) by
   ! Beyer-Swinehart direct counting
   !------------------------------------------------------------------
   subroutine sum_count(w_mode, e, n_states)
      real(8), intent(in) :: w_mode(:)     ! normal-mode frequency table [energy bins]
      real(8), intent(in) :: e             ! probe energy [energy bins]
      real(8), intent(out) :: n_states     ! cumulative count N(e) [count]
      real(8), parameter :: del = 1.0d0    ! bin width [bin] (the in-bin contract)
      real(8), allocatable :: t(:)         ! counting work array (exact integers below 2^53)
      real(8) :: frac                      ! the interpolation remainder within one bin
      integer :: ne(size(w_mode))          ! per-mode effective integer frequencies [bin]
      integer :: maxt, nd, i, k
      ! 1. negative probe: no states below the zero-point level
      if (e < 0.0d0) then
         n_states = 0.0d0
         return
      end if
      ! 2. counting grid (probe rounded up + 20 margin bins) + per-mode effective
      !    integer frequencies
      maxt = int(e/del + 0.499999d0) + 20  ! grid up, plus 20 bins of margin
      do i = 1, size(w_mode)
         ne(i) = int(w_mode(i)/del + 0.499999d0)
      end do
      allocate (t(maxt))
      t(1) = 1.0d0                         ! the zero-point state
      t(2:maxt) = 0.0d0
      ! 3. counting core: one forward pass per mode, t(ne_i + k) += t(k)
      do i = 1, size(w_mode)
         do k = 1, maxt - ne(i)
            t(ne(i) + k) = t(ne(i) + k) + t(k)
         end do
      end do
      ! 4. prefix sum: exact levels -> N(E)
      do i = 2, maxt
         t(i) = t(i) + t(i - 1)
      end do
      ! 5. linear interpolation of N(E) within the bin
      nd = int(e/del)
      frac = e/del - dble(nd)
      n_states = t(nd + 1) + (t(nd + 2) - t(nd + 1))*frac
      deallocate (t)
   end subroutine sum_count

   !------------------------------------------------------------------
   ! eig_out(a, eig, unit) - formatted eigenvalue/vector output in 6-column
   ! blocks (the matrix's column convention restored after writing)
   !------------------------------------------------------------------
   subroutine eig_out(a, eig, unit)
      real(8), intent(inout) :: a(:,:)  ! eigenvector matrix (entry: vectors by column, the
                                         ! linalg eig_sym convention; temporarily transposed
                                         ! to rows before writing, restored in place after)
      real(8), intent(in) :: eig(:)     ! eigenvalues (ascending, aligned with the columns of a)
      integer, intent(in) :: unit       ! output unit number [-]
      integer :: n, i, l, k, k6
      ! 1. guards: the matrix must carry at least an n x n block
      n = size(eig)
      if (size(a, 1) < n .or. size(a, 2) < n) then
         write (0, '(a,i0,a,i0,a,i0)') 'eig_out: the eigenvector matrix is ', size(a, 1), 'x', &
            size(a, 2), ' but the eigenvalue table carries ', n, ' values'
         write (0, '(a)') 'eig_out: fatal (the table and the values disagree)'
         stop 1
      end if
      call swap_triangles(a, n)         ! column vectors temporarily row-stored
      ! 2. one 6-column block per pass: eigenvalues, then one row per component
      k = 0
      do while (k < n)
         l = k + 1
         k = min(k + 6, n)              ! one 6-column block l..k
         write (unit, '(//2x,6(9x,i3))') (i, i=l, k)
         write (unit, '(5x,1p6e12.4)') (eig(i), i=l, k)
         write (unit, '(/)')
         do i = 1, n
            write (unit, '(i3,1x,6f12.6)') i, (a(i, k6), k6=l, k)
         end do
      end do
      call swap_triangles(a, n)         ! the caller's column convention restored
   end subroutine eig_out

   !------------------------------------------------------------------
   ! swap_triangles(a, n) - in-place transpose of the leading n x n block
   !                         (upper/lower triangle swap; the diagonal rides along)
   subroutine swap_triangles(a, n)
      real(8), intent(inout) :: a(:,:)
      integer, intent(in) :: n
      real(8) :: tmp
      integer :: i, j
      do i = 1, n
         do j = 1, i
            tmp = a(i, j)
            a(i, j) = a(j, i)
            a(j, i) = tmp
         end do
      end do
   end subroutine swap_triangles
end module density
