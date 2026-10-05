!=====================================================================
! check_rng.f90 - random-stream criteria program (criteria defined FIRST in
!                 this header; this program executes them)
! Compile: gfortran scripts/check_rng.f90 src/utils/rng.f90 -o build/check_rng.exe
! Criteria:
!   (a) determinism  - same seed, two init+draw rounds, 1e6 draws bit-identical
!   (b) uniformity   - |mean-0.5| < 1e-3; 20-bucket chi-square within df=19 crit 30.14
!   (c) gaussianity  - |mean| < 5/sqrt(N); |var-1| < 1e-2; 10 equiprobable quantile
!                      buckets chi-square < df=9 crit 16.92
!   (d) gamma(2,1)   - |mean-2| and |var-2| within 5-sigma sample error
!   (e) cache health - gaussian stats un-degraded when interleaved with rng_u draws
!=====================================================================
program check_rng
   use rng, only: rng_init, rng_u, rng_gauss, rng_gamma
   implicit none
   integer, parameter :: n = 1000000, nb = 20, nd = 10
   real(8) :: x(n), s, mean, varr, chi2, expct
   integer :: i, k, bucket(nb)
   ! standard-normal decile edges Phi^{-1}(k/10), k=1..9 (4-digit accuracy suffices)
   real(8), parameter :: zedge(9) = (/ -1.2816d0, -0.8416d0, -0.5244d0, -0.2533d0, &
                                        0.0d0, 0.2533d0, 0.5244d0, 0.8416d0, 1.2816d0 /)
   integer :: nfail = 0

   ! ---------- (a) determinism ----------
   call rng_init(12345)
   do i = 1, n
      x(i) = rng_u()
   end do
   call rng_init(12345)
   do i = 1, n
      s = rng_u()
      if (s /= x(i)) then
         print '(a,i0)', 'FAIL (a) determinism: first mismatch at draw ', i
         nfail = nfail + 1
         exit
      end if
   end do
   if (nfail == 0) print '(a)', '(a) determinism: PASS (1e6 draws bit-identical, seed 12345)'

   ! ---------- (b) uniformity (reuse the x array from (a)) ----------
   mean = sum(x)/real(n, 8)
   varr = sum(x*x)/real(n, 8) - mean*mean
   bucket = 0
   do i = 1, n
      k = int(x(i)*real(nb, 8)) + 1
      if (k > nb) k = nb
      bucket(k) = bucket(k) + 1
   end do
   expct = real(n, 8)/real(nb, 8)
   chi2 = 0.0d0
   do k = 1, nb
      chi2 = chi2 + (real(bucket(k), 8) - expct)**2/expct
   end do
   print '(a,f12.8,a,f10.8,a,f8.3)', '(b) uniformity: mean = ', mean, '  var = ', varr, &
      '  chi2(20) = ', chi2
   if (abs(mean - 0.5d0) >= 1.0d-3) then
      print '(a)', 'FAIL (b): |mean-0.5| >= 1e-3'
      nfail = nfail + 1
   end if
   if (chi2 > 30.14d0) then
      print '(a)', 'FAIL (b): chi2 > 30.14 (df=19, alpha=0.05)'
      nfail = nfail + 1
   end if

   ! ---------- (c) gaussianity ----------
   call rng_init(777)
   do i = 1, n
      x(i) = rng_gauss()
   end do
   mean = sum(x)/real(n, 8)
   varr = sum(x*x)/real(n, 8) - mean*mean
   bucket = 0
   do i = 1, n
      k = 1
      do while (k <= 9 .and. x(i) > zedge(k))
         k = k + 1
      end do
      bucket(k) = bucket(k) + 1
   end do
   expct = real(n, 8)/real(nd, 8)
   chi2 = 0.0d0
   do k = 1, nd
      chi2 = chi2 + (real(bucket(k), 8) - expct)**2/expct
   end do
   print '(a,f12.8,a,f10.8,a,f8.3)', '(c) gaussianity: mean = ', mean, '  var = ', varr, &
      '  chi2(10 deciles) = ', chi2
   if (abs(mean) >= 5.0d0/sqrt(real(n, 8))) then
      print '(a)', 'FAIL (c): |mean| >= 5/sqrt(N)'
      nfail = nfail + 1
   end if
   if (abs(varr - 1.0d0) >= 1.0d-2) then
      print '(a)', 'FAIL (c): |var-1| >= 1e-2'
      nfail = nfail + 1
   end if
   if (chi2 > 16.92d0) then
      print '(a)', 'FAIL (c): chi2 > 16.92 (df=9, alpha=0.05)'
      nfail = nfail + 1
   end if

   ! ---------- (d) gamma(2,1) moments ----------
   call rng_init(31415)
   do i = 1, n
      x(i) = rng_gamma(2)
   end do
   mean = sum(x)/real(n, 8)
   varr = sum(x*x)/real(n, 8) - mean*mean
   print '(a,f12.8,a,f10.8)', '(d) gamma(2,1): mean = ', mean, '  var = ', varr
   if (abs(mean - 2.0d0) >= 5.0d0*sqrt(2.0d0/real(n, 8))) then
      print '(a)', 'FAIL (d): |mean-2| >= 5-sigma'
      nfail = nfail + 1
   end if
   if (abs(varr - 2.0d0) >= 5.0d0*2.0d0*sqrt(2.0d0/real(n, 8))) then
      print '(a)', 'FAIL (d): |var-2| >= 5-sigma'
      nfail = nfail + 1
   end if

   ! ---------- (e) cache health under interleaving ----------
   call rng_init(999)
   do i = 1, n/2
      s = rng_u()
      x(i) = rng_gauss()
   end do
   mean = sum(x(1:n/2))/real(n/2, 8)
   varr = sum(x(1:n/2)**2)/real(n/2, 8) - mean*mean
   print '(a,f12.8,a,f10.8)', '(e) interleaved gauss: mean = ', mean, '  var = ', varr
   if (abs(mean) >= 5.0d0/sqrt(real(n/2, 8)) .or. abs(varr - 1.0d0) >= 2.0d-2) then
      print '(a)', 'FAIL (e): gaussian stats degraded under interleaving'
      nfail = nfail + 1
   end if

   print '(a,i0)', '=== check_rng: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0
end program check_rng
