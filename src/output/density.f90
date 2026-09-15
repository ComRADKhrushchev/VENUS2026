!=====================================================================
! analysis/density.f90 - L4 analysis: quantum density of states (Beyer-Swinehart direct
!                         counting) and
!                         formatted eigenvalue/vector output
! Reference: src_VENUS/DENQ.f:10-18 (DENQ: cumulative state counts at the two probe energies
!       E +/- 50, density =
!       (SUM1-SUM2)/100 central difference), :25-96 (SUMQ direct counting: energy bin width
!       DEL = 1.0
!       :40, +20 bins of margin :50-51, per-mode effective integer frequencies
!       NEFIN = int(W/DEL+0.499999)
!       :57-59, T-array addition sequence T(j+ne_i)+=T(j) :63-77, cumulation T(i)+=T(i-1)
!       :81-83, E/DEL linear interpolation :86-89, negative-probe early return :44-47, T
!       allocation/release :53/:93);
!       src_VENUS/EIGOUT.f:1-39 (eigenvalue output: upper/lower triangle swap to take row
!       vectors :13-19,
!       6-column blocks writing eigenvalues EIG and per-row vector components :21-31,
!       post-write transpose restore :33-38)
!=====================================================================
module density
   implicit none
   private
   public :: dens_states, eig_out
contains
   !------------------------------------------------------------------
   ! dens_states(w_mode, e_probe, rho) - density of states (Beyer-Swinehart direct counting +
   !                                      central difference)
   ! (see previous line)
   ! Reference: DENQ.f:10-18 (the differencing shell); DENQ.f:25-96 (SUMQ direct counting -
   ! inlined as step 1)
   subroutine dens_states(w_mode, e_probe, rho)
      real(8), intent(in) :: w_mode(:)  ! normal-mode frequency table [energy bins] (old WK -
                                        ! effective frequencies in units of the bin width DEL;
                                        ! cm^-1 -> bin conversion done via
      ! consts in the assembly phase - finalized in phase C)
      real(8), intent(in) :: e_probe    ! probe energy [energy bins] (old EN)
      real(8), intent(out) :: rho       ! density of states dN/dE [1/bin] (old RHO)
      real(8) :: s_hi, s_lo             ! the cumulative counts at the two probes
      call sum_count(w_mode, e_probe + 50.0d0, s_hi)
      call sum_count(w_mode, e_probe - 50.0d0, s_lo)
      rho = (s_hi - s_lo)/100.0d0
   end subroutine dens_states

   !------------------------------------------------------------------
   ! sum_count(w_mode, e, n_states) - cumulative quantum state count N(e) by
   !                                  Beyer-Swinehart direct counting (the
   !                                  differencing shell's step 1, inlined here)
   subroutine sum_count(w_mode, e, n_states)
      real(8), intent(in) :: w_mode(:)     ! normal-mode frequency table [energy bins]
      real(8), intent(in) :: e             ! probe energy [energy bins]
      real(8), intent(out) :: n_states     ! cumulative count N(e) [count]
      real(8), parameter :: del = 1.0d0    ! bin width [bin] (the in-bin contract)
      real(8), allocatable :: t(:)         ! counting work array (exact integers below 2^53)
      real(8) :: frac                      ! the interpolation remainder within one bin
      integer :: ne(size(w_mode))          ! per-mode effective integer frequencies [bin]
      integer :: maxt, nd, i, k
      if (e < 0.0d0) then
         n_states = 0.0d0                  ! no states below the zero-point level
         return
      end if
      maxt = int(e/del + 0.499999d0) + 20  ! grid up, plus 20 bins of margin
      do i = 1, size(w_mode)
         ne(i) = int(w_mode(i)/del + 0.499999d0)
      end do
      allocate (t(maxt))
      t(1) = 1.0d0                         ! the zero-point state
      t(2:maxt) = 0.0d0
      do i = 1, size(w_mode)               ! one forward pass per mode: the
         do k = 1, maxt - ne(i)            ! counting core t(ne_i + k) += t(k)
            t(ne(i) + k) = t(ne(i) + k) + t(k)
         end do
      end do
      do i = 2, maxt                       ! prefix sum: exact levels -> N(E)
         t(i) = t(i) + t(i - 1)
      end do
      nd = int(e/del)
      frac = e/del - dble(nd)
      n_states = t(nd + 1) + (t(nd + 2) - t(nd + 1))*frac
      deallocate (t)
   end subroutine sum_count

   !------------------------------------------------------------------
   ! eig_out(a, eig, unit) - formatted eigenvalue/vector output (6-column blocks; restored after writing)
   ! Reference: EIGOUT.f:13-39
   subroutine eig_out(a, eig, unit)
      real(8), intent(inout) :: a(:,:)  ! eigenvector matrix (entry: vectors by column - the
                                         ! linalg eig_sym convention; temporarily transposed to
                                         ! rows before writing, restored in place after)
      real(8), intent(in) :: eig(:)     ! eigenvalues (ascending, aligned with the columns of a;
                                         ! the old global EIG now passed as a dummy argument)
      integer, intent(in) :: unit       ! output unit number [-] (old dummy argument IP)
      integer :: n, i, l, k, k6
      n = size(eig)
      if (size(a, 1) < n .or. size(a, 2) < n) then
         write (0, '(a,i0,a,i0,a,i0)') 'eig_out: the eigenvector matrix is ', size(a, 1), 'x', &
            size(a, 2), ' but the eigenvalue table carries ', n, ' values'
         write (0, '(a)') 'eig_out: loud-stop (the table and the values disagree)'
         stop 1
      end if
      call swap_triangles(a, n)         ! column vectors temporarily row-stored
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
