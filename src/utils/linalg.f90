!=====================================================================
! linalg.f90 - L5 instrument: symmetric-matrix eigendecomposition (eigenvalues ascending +
!              eigenvectors by column)
! Reference: src_VENUS/EIGN.f:1-248 (EIGN - Givens-Householder tridiagonalization + shifted QL)
!       consumers: inertia-tensor diagonalization venus_input.f90:553/:814; force-constant
!       matrix FGMTRX.f:35/:60
!=====================================================================
module linalg
   implicit none
   private
   public :: eig_sym
contains
   !------------------------------------------------------------------
   ! eig_sym(a, d) - symmetric-matrix eigenvalues/vectors (a read by lower triangle and destroyed;
   ! eigenvalues ascending)
   ! Reference: src_VENUS/EIGN.f:27-99 (tridiagonalization); :103-186 (adjoint identity + shifted QL
   ! iteration);
   !       :192-229 (transposed sort + vector reordering); :231-247 (back-substitution
   !       accumulating transforms)
   subroutine eig_sym(a, d)
      real(8), intent(inout) :: a(:,:)  ! symmetric matrix n x n (lower triangle on input,
                                        ! destroyed;
      ! on exit, eigenvectors by column, ascending to
      ! match d)
      real(8), intent(out)   :: d(:)    ! eigenvalues (n of them, algebraically ascending)
      ! Pseudocode:
      ! 1. n = size(d); n=0 returns immediately (old :20)
      ! 2. Givens-Householder tridiagonalization: build the W vector and scale factor D
      ! column by column,
      !    the P/Q vector method forms PAP, giving the diagonal gamma, off-diagonal beta, and
      !    the transform-factor array
      !    (uses and destroys only the lower triangle, old :27-99)
      ! 3. Adjoint identity matrix; shifted QL iteration (the shift is the root of the corner
      ! 2x2 closest to the lower diagonal element,
      !    old :109-186); the convergence threshold (old dummy RHO) takes an internal default,
      !    to be promoted to a dummy argument in phase C
      !    (see previous line)
      ! 4. Sort the transposed eigenvalues ascending, swapping eigenvector columns in step (old :192-229)
      ! 5. Back-substitute the eigenvectors accumulating the Householder transforms (old :231-247)
      ! 6. d <- the sorted eigenvalues; vectors stored by column in a
   end subroutine eig_sym
end module linalg
