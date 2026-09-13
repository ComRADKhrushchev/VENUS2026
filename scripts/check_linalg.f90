!=====================================================================
! check_linalg.f90 - C1.2 acceptance criteria program (criteria defined FIRST in
!                    docs/plans/C1.md)
! Compile: gfortran -Jbuild src/utils/linalg.f90 scripts/check_linalg.f90 -o build/check_linalg.exe
! Criteria:
!   (a) analytic diagonal 3x3 - known eigenpairs, machine-precision residual
!   (b) analytic non-trivial 2x2 - the closed-form eigenvalues/vectors
!   (c) constructed symmetric 10x10 - V^T V = I, A v = lambda v residuals < 1e-12,
!       eigenvalues algebraically ascending
!   (d) degenerate case (I3, triple eigenvalue) - orthonormality and A v = v still hold
!=====================================================================
program check_linalg
   use linalg, only: eig_sym
   implicit none
   real(8) :: a(10, 10), d(10), r(10, 10), res
   real(8) :: lam2(2), cth, sth, aa, bb, cc, disc
   integer :: i, j, k, nfail = 0

   ! ---------- (a) diagonal 3x3 ----------
   a = 0.0d0
   a(1,1) = 1.0d0; a(2,2) = 2.0d0; a(3,3) = 3.0d0
   call eig_sym(a(1:3,1:3), d(1:3))
   if (abs(d(1)-1.0d0) > 1.0d-14 .or. abs(d(2)-2.0d0) > 1.0d-14 .or. &
       abs(d(3)-3.0d0) > 1.0d-14) then
      print '(a,3f12.8)', 'FAIL (a) eigenvalues: ', d(1:3)
      nfail = nfail + 1
   else
      print '(a)', '(a) diagonal 3x3: PASS (eigenvalues 1,2,3 at machine precision)'
   end if

   ! ---------- (b) analytic 2x2 ----------
   aa = 2.0d0; bb = 5.0d0; cc = 1.5d0
   a = 0.0d0
   a(1,1) = aa; a(2,2) = bb; a(2,1) = cc      ! lower triangle (the contract's read side)
   call eig_sym(a(1:2,1:2), d(1:2))
   disc = sqrt(((aa-bb)/2.0d0)**2 + cc*cc)
   lam2(1) = (aa+bb)/2.0d0 - disc
   lam2(2) = (aa+bb)/2.0d0 + disc
   if (abs(d(1)-lam2(1)) > 1.0d-12 .or. abs(d(2)-lam2(2)) > 1.0d-12) then
      print '(a,2f14.10,a,2f14.10)', 'FAIL (b): got ', d(1:2), ' expect ', lam2
      nfail = nfail + 1
   else
      ! eigenvector of the lower branch: (c, lam-a) normalized, stored as column 1
      cth = cc/sqrt(cc*cc + (lam2(1)-aa)**2)
      sth = (lam2(1)-aa)/sqrt(cc*cc + (lam2(1)-aa)**2)
      if (abs(abs(a(1,1)*cth + a(2,1)*sth) - 1.0d0) > 1.0d-10) then
         print '(a)', 'FAIL (b): eigenvector direction'
         nfail = nfail + 1
      else
         print '(a,2f14.10)', '(b) analytic 2x2: PASS (lambda = ', d(1:2), ')'
      end if
   end if

   ! ---------- (c) constructed symmetric 10x10 ----------
   do j = 1, 10
      do i = j, 10
         a(i,j) = sin(real(i*11 + j*7, 8))*0.5d0 + cos(real(i+j, 8))*0.3d0
      end do
   end do
   call eig_sym(a, d)
   ! ascending?
   do i = 1, 9
      if (d(i) > d(i+1) + 1.0d-13) then
         print '(a,i0)', 'FAIL (c): not ascending at ', i
         nfail = nfail + 1
         exit
      end if
   end do
   ! V^T V = I (columns orthonormal)
   r = matmul(transpose(a), a)
   res = 0.0d0
   do j = 1, 10
      do i = 1, 10
         res = max(res, abs(r(i,j) - merge(1.0d0, 0.0d0, i == j)))
      end do
   end do
   if (res > 1.0d-12) then
      print '(a,es10.2)', 'FAIL (c): V^T V residual = ', res
      nfail = nfail + 1
   end if
   ! A v = lambda v (rebuild the original matrix - a was destroyed; mirror the lower
   ! triangle to the upper - dot_product needs the full symmetric matrix)
   do j = 1, 10
      do i = j, 10
         r(i,j) = sin(real(i*11 + j*7, 8))*0.5d0 + cos(real(i+j, 8))*0.3d0
      end do
   end do
   do j = 1, 10
      do i = 1, j - 1
         r(i,j) = r(j,i)
      end do
   end do
   res = 0.0d0
   do k = 1, 10
      do i = 1, 10
         res = max(res, abs(dot_product(r(i,:), a(:,k)) - d(k)*a(i,k)))
      end do
   end do
   if (res > 1.0d-12) then
      print '(a,es10.2)', 'FAIL (c): A v - lambda v residual = ', res
      nfail = nfail + 1
   end if
   if (nfail == 0) print '(a,es10.2,a,es10.2)', '(c) 10x10: PASS (VtV res ', &
      maxval(abs(r)), ' omitted; Av-lv and VtV residuals < 1e-12)'

   ! ---------- (d) degenerate I3 ----------
   a = 0.0d0
   a(1,1) = 1.0d0; a(2,2) = 1.0d0; a(3,3) = 1.0d0
   call eig_sym(a(1:3,1:3), d(1:3))
   res = 0.0d0
   do j = 1, 3
      do i = 1, 3
         res = max(res, abs(dot_product(a(1:3,i), a(1:3,j)) - merge(1.0d0, 0.0d0, i == j)))
      end do
   end do
   if (res > 1.0d-12 .or. abs(d(1)-1.0d0) > 1.0d-13 .or. abs(d(3)-1.0d0) > 1.0d-13) then
      print '(a,es10.2)', 'FAIL (d): degenerate case residual = ', res
      nfail = nfail + 1
   else
      print '(a)', '(d) degenerate I3: PASS (orthonormal basis, lambda = 1 x3)'
   end if

   print '(a,i0)', '=== check_linalg: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0
end program check_linalg
