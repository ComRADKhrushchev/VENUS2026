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
      ! Pseudocode:
      ! 1. Cumulative state count N(E) (Beyer-Swinehart direct counting, old SUMQ :25-96):
      ! 1a. n_grid = int(e_probe/DEL + 0.499999) + 20 (DEL = 1.0 bin; +20 bins
      !     of margin, old :40-51; e_probe<0 -> N=0 early return, old :44-47)
      ! 1b. Per-mode effective integer frequency ne_i = int(w_mode(i)/DEL + 0.499999) (old :57-59)
      ! 1c. Work array t(1..n_grid): t(1)=1, rest zeroed (old :63-66)
      ! 1d. Per-mode addition sequence: t(j+ne_i) += t(j), j = 1..n_grid-ne_i
      !     (single forward pass per mode - the Beyer-Swinehart core, old :70-77)
      ! 1e. Cumulation: t(i) += t(i-1) (per-bin prefix sum - turning it into N(E), old :81-83)
      ! 1f. N(E) = linear interpolation of t at e_probe/DEL (old :86-89);
      !     release t (old :93)
      ! 2. rho = (N(e_probe+50) - N(e_probe-50))/100 (central difference over +/- 50 bins
      !     divided by the 100-bin window width, old DENQ.f:10-18; step 1 runs once per probe)
   end subroutine dens_states

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
      ! Pseudocode:
      ! 1. Swap the upper/lower triangles of a (column-stored vectors temporarily row-stored -
      !    old :13-19)
      ! 2. Loop over 6-column blocks (old :21-31): block-header column-number line; eigenvalue
      !    line (eig(block columns));
      !    per row i write the row number + vector components a(i, block columns)
      ! 3. Swap the triangles of a again to restore (old :33-38)
   end subroutine eig_out
end module density
