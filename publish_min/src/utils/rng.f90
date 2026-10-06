!=====================================================================
! rng.f90 - uniform random stream (Knuth multiplicative congruential +
!           shuffle table) and distribution transforms (normal / gamma)
! Design:
!   rng_init establishes the stream from a decimal seed (base-256 conversion,
!   warm-up, shuffle-table fill); rng_u draws the next uniform in [0,1);
!   rng_gauss (Box-Muller polar form, cached pair) and rng_gamma draw from
!   the same single stream. Same seed -> same uniform sequence (the
!   reproducibility anchor).
!=====================================================================
module rng
   implicit none
   private
   public :: rng_init, rng_u, rng_gauss, rng_gamma

   ! Stream state (internalized in the module)
   integer, parameter :: dig_base = 256     ! digit base 2^8
   integer :: seed_dig(9) = 0               ! stream seed as 8 base-256 digits + the 9th
                                            ! carry-absorber digit (the doubling carry may
                                            ! reach digit 9 - absorbed, never read back)
   real(8) :: tbl_shuf(100) = 0.0d0         ! shuffle table of 100 entries
   ! Box-Muller cached-pair state:
   integer :: i_cache = 0                   ! cache flag (0 = none, 1 = cached)
   real(8) :: g_cache = 0.0d0               ! cached second normal deviate
contains
   !------------------------------------------------------------------
   ! rng_step() - one multiplicative-congruential step: seed_dig <- A*seed_dig (mod 2^64),
   !              returned as [0,1) (private here)
   !------------------------------------------------------------------
   real(8) function rng_step()
      integer :: ia(8), id(16), ip, i, j, k
      real(8), parameter :: bi = 3.90625d-3 ! 1/256
      ia = (/ 45, 127, 149, 76, 45, 244, 81, 88 /)   ! base-256 multiplier digits (A = Knuth's
                                                    ! 6364136223846793005)
      id(1:8) = 0                                     ! additive constant C = 0
      id(9:16) = 0
      ! form A*seed_dig + C: digit-wise convolution with carry
      do j = 1, 8
         do i = 1, 9 - j
            k = j + i - 1
            ip = ia(j)*seed_dig(i)
            do
               ip = ip + id(k)
               id(k) = mod(ip, dig_base)
               ip = ip/dig_base
               if (ip == 0 .or. k == 8) exit
               k = k + 1
            end do
         end do
      end do
      seed_dig(1:8) = id(1:8)
      ! fold the 8 digits as base 256 into [0,1)
      rng_step = dble(seed_dig(1))
      do i = 2, 8
         rng_step = dble(seed_dig(i)) + rng_step*bi
      end do
      rng_step = rng_step*bi
   end function rng_step

   !------------------------------------------------------------------
   ! rng_init(i_seed) - establish a new random stream from a decimal seed (base-256
   !                    conversion + warm-up + filling the shuffle table)
   !------------------------------------------------------------------
   subroutine rng_init(i_seed)
      integer, intent(in) :: i_seed   ! decimal seed (0 <= i_seed <= 2^31-1)
      integer :: is, i
      ! 1. decimal -> 8 base-256 digits
      is = i_seed
      do i = 1, 8
         seed_dig(i) = mod(is, dig_base)
         is = is/dig_base
      end do
      seed_dig(9) = 0
      ! 2. double digit-wise with carry (an even value)
      do i = 1, 8
         seed_dig(i) = seed_dig(i) + seed_dig(i)
      end do
      do i = 1, 8
         do while (seed_dig(i) >= dig_base)
            seed_dig(i) = seed_dig(i) - dig_base
            seed_dig(i + 1) = seed_dig(i + 1) + 1
         end do
      end do
      ! 3. +1 (an odd stream seed)
      seed_dig(1) = seed_dig(1) + 1
      do i = 1, 8
         do while (seed_dig(i) >= dig_base)
            seed_dig(i) = seed_dig(i) - dig_base
            seed_dig(i + 1) = seed_dig(i + 1) + 1
         end do
      end do
      ! 4. warm-up: fill the shuffle table TWICE in a row (same loop body, different results)
      do i = 1, 100
         tbl_shuf(i) = rng_step()
      end do
      do i = 1, 100
         tbl_shuf(i) = rng_step()
      end do
   end subroutine rng_init

   !------------------------------------------------------------------
   ! rng_u() - draw the next uniform deviate in [0,1) (shuffle-table draw +
   !           multiplicative congruential step to refill)
   !------------------------------------------------------------------
   real(8) function rng_u()
      integer :: j
      j = int(99.0d0*tbl_shuf(100)) + 1     ! shuffle pointer
      rng_u = tbl_shuf(100)                 ! the draw
      tbl_shuf(100) = tbl_shuf(j)           ! replace the drawn slot with table entry j
      tbl_shuf(j) = rng_step()              ! refill through the congruential step
   end function rng_u

   !------------------------------------------------------------------
   ! rng_gauss() - standard normal deviate N(0,1) (Box-Muller polar form, cached pair;
   !               uniform source = the single stream, one rng_u per component)
   !------------------------------------------------------------------
   real(8) function rng_gauss()
      real(8) :: fac, rsq, v1, v2
      if (i_cache == 0) then
         do
            v1 = 2.0d0*rng_u() - 1.0d0
            v2 = 2.0d0*rng_u() - 1.0d0
            rsq = v1*v1 + v2*v2
            if (rsq < 1.0d0 .and. rsq > 0.0d0) exit
         end do
         fac = sqrt(-2.0d0*log(rsq)/rsq)
         g_cache = v1*fac
         rng_gauss = v2*fac
         i_cache = 1
      else
         rng_gauss = g_cache
         i_cache = 0
      end if
   end function rng_gauss

   !------------------------------------------------------------------
   ! rng_gamma(order) - integer-order gamma deviate (order=1 degenerates to the
   !                    unit exponential; uniform source = the stream)
   !------------------------------------------------------------------
   real(8) function rng_gamma(order)
      integer, intent(in) :: order    ! gamma distribution order (positive integer; 1 = exponential)
      integer :: j
      real(8) :: dum
      dum = 1.0d0
      do j = 1, order
         dum = dum*rng_u()
      end do
      rng_gamma = -log(dum)
   end function rng_gamma
end module rng
