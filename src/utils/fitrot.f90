!=====================================================================
! fitrot.f90 - mass-weighted quaternion best-fit rotation + rotational
!              Eckart residual (pure array form, no system knowledge)
! Design:
!   quat_fit builds the 4x4 quaternion quadric form (legacy QTRFIT's C
!   matrix, Cyclic Jacobi retired) and takes its LARGEST eigenvector via
!   linalg's eig_sym (ascending -> column 4); the rotation matrix follows
!   the legacy Q2MAT form. Inputs must be COM-subtracted by the caller
!   (translation is not fitted). eckart_rot reports the rotational Eckart
!   condition residual; rot3 applies u to one 3-vector in place.
!=====================================================================
module fitrot
   use linalg, only: eig_sym
   implicit none
   private
   public :: quat_fit, eckart_rot, rot3
contains
   !------------------------------------------------------------------
   ! quat_fit(q_ref, q, mass, u, rms) - minimize sum m_i |u q_i - q_ref_i|^2
   !------------------------------------------------------------------
   subroutine quat_fit(q_ref, q, mass, u, rms)
      real(8), intent(in)  :: q_ref(:,:), q(:,:), mass(:)
      real(8), intent(out) :: u(3,3), rms
      integer :: i
      real(8) :: s(3,3), c4(4,4), d(4), q4(4), wsum, r
      ! 1. correlation S = sum m q qref^T
      s = 0.0d0
      do i = 1, size(mass)
         s = s + mass(i)*outer(q(:,i), q_ref(:,i))
      end do
      ! 2. quadric form, lower triangle only (eig_sym's read contract)
      c4 = 0.0d0
      c4(1,1) =  s(1,1)+s(2,2)+s(3,3)
      c4(2,1) =  s(3,2)-s(2,3)
      c4(2,2) =  s(1,1)-s(2,2)-s(3,3)
      c4(3,1) =  s(1,3)-s(3,1)
      c4(3,2) =  s(1,2)+s(2,1)
      c4(3,3) = -s(1,1)+s(2,2)-s(3,3)
      c4(4,1) =  s(2,1)-s(1,2)
      c4(4,2) =  s(1,3)+s(3,1)
      c4(4,3) =  s(2,3)+s(3,2)
      c4(4,4) = -s(1,1)-s(2,2)+s(3,3)
      ! 3. largest eigenvector -> quaternion -> rotation matrix
      call eig_sym(c4, d)
      q4 = c4(:,4)
      u = q2mat(q4)
      ! 4. residual (mass-weighted RMS over the total mass)
      wsum = sum(mass)
      r = 0.0d0
      do i = 1, size(mass)
         r = r + mass(i)*sum((matmul(u, q(:,i)) - q_ref(:,i))**2)
      end do
      rms = sqrt(r/wsum)
   end subroutine quat_fit

   !------------------------------------------------------------------
   ! eckart_rot(q_ref, q, mass) - rotational Eckart condition residual
   !                               [amu A^2]: |sum m qref x (q - qref)|
   !------------------------------------------------------------------
   function eckart_rot(q_ref, q, mass) result(res)
      real(8), intent(in) :: q_ref(:,:), q(:,:), mass(:)
      real(8) :: res
      integer :: i
      real(8) :: l(3)
      l = 0.0d0
      do i = 1, size(mass)
         l = l + mass(i)*cross3(q_ref(:,i), q(:,i) - q_ref(:,i))
      end do
      res = sqrt(sum(l**2))
   end function eckart_rot

   !------------------------------------------------------------------
   ! rot3(u, v) - rotate one 3-vector in place (v = u*v)
   !------------------------------------------------------------------
   subroutine rot3(u, v)
      real(8), intent(in)    :: u(3,3)
      real(8), intent(inout) :: v(3)
      real(8) :: w(3)
      w = matmul(u, v)
      v = w
   end subroutine rot3

   !------------------------------------------------------------------
   ! q2mat(q4) - rotation matrix of a normalized quaternion (legacy Q2MAT,
   !             0-based rewritten 1-based; u*q column convention)
   !------------------------------------------------------------------
   function q2mat(q4) result(u)
      real(8), intent(in) :: q4(4)
      real(8) :: u(3,3)
      u(1,1) = q4(1)**2 + q4(2)**2 - q4(3)**2 - q4(4)**2
      u(2,1) = 2.0d0*(q4(2)*q4(3) - q4(1)*q4(4))
      u(3,1) = 2.0d0*(q4(2)*q4(4) + q4(1)*q4(3))
      u(1,2) = 2.0d0*(q4(3)*q4(2) + q4(1)*q4(4))
      u(2,2) = q4(1)**2 - q4(2)**2 + q4(3)**2 - q4(4)**2
      u(3,2) = 2.0d0*(q4(3)*q4(4) - q4(1)*q4(2))
      u(1,3) = 2.0d0*(q4(4)*q4(2) - q4(1)*q4(3))
      u(2,3) = 2.0d0*(q4(4)*q4(3) + q4(1)*q4(2))
      u(3,3) = q4(1)**2 - q4(2)**2 - q4(3)**2 + q4(4)**2
   end function q2mat

   !------------------------------------------------------------------
   ! outer(a, b) - dyad a*b^T
   !------------------------------------------------------------------
   function outer(a, b) result(c)
      real(8), intent(in) :: a(3), b(3)
      real(8) :: c(3,3)
      integer :: j, k
      do j = 1, 3
         do k = 1, 3
            c(k,j) = a(k)*b(j)
         end do
      end do
   end function outer

   !------------------------------------------------------------------
   ! cross3(a, b) - vector product a x b
   !------------------------------------------------------------------
   function cross3(a, b) result(c)
      real(8), intent(in) :: a(3), b(3)
      real(8) :: c(3)
      c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
   end function cross3

end module fitrot
