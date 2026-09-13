!=====================================================================
! rng.f90 - L5 instrument: uniform random stream (Knuth multiplicative congruential + shuffle) and
!           distribution transforms (normal / gamma)
! Reference: src_VENUS/RANDST.f:47-105 (stream initialization); src_VENUS/RAND0.f:19-23
! (shuffle-table draws);
!       src_VENUS/RAND1.f:59-106 (multiplicative congruential stepping, 8 base-256 digits);
!       src_TEST/interface_TEST.f90:112-136 (gasdev normal pairs);
!       src_VENUS/GAMA.f:7-13 (gamma/exponential)
! Port notes (C1, per docs/plans/C1.md D1/D2): the RANDST/RAND0/RAND1 arithmetic is carried
!       over VERBATIM (same seed -> same uniform sequence as the old library, the V-phase
!       same-ISEED anchor); the two RANDOM_NUMBER consumers are switched onto the single
!       stream (gasdev: two rng_u draws per pair instead of one random_number(u(2)); GAMA:
!       rng_u instead of the intrinsic) - their sequences therefore differ from the old
!       library's, the P-C1 ruling carried into the new tree; the old init-time banner
!       print (RANDST :910) and the gasdev NaN debug print are NOT carried over.
!=====================================================================
module rng
   implicit none
   private
   public :: rng_init, rng_u, rng_gauss, rng_gamma

   ! Stream state (internalized in the module - K1 plan §2: the old RANLST/ISEED3/IBFCTR become
   ! module-private)
   integer, parameter :: dig_base = 256     ! digit base 2^8 (old IBFCTR)
   integer :: seed_dig(9) = 0               ! stream seed as 8 base-256 digits + the 9th
                                            ! carry-absorber digit (old ISEED3(1:9); RANDST's
                                            ! doubling carry may reach digit 9 - absorbed,
                                            ! never read back, the old behavior verbatim)
   real(8) :: tbl_shuf(100) = 0.0d0         ! shuffle table of 100 entries (old RANLST)
   ! Box-Muller cached-pair state (the iset/gset of old gasdev):
   integer :: i_cache = 0                   ! cache flag (0 = none, 1 = cached; old iset)
   real(8) :: g_cache = 0.0d0               ! cached second normal deviate (old gset)
contains
   !------------------------------------------------------------------
   ! rng_step() - one multiplicative-congruential step: seed_dig <- A*seed_dig (mod 2^64),
   !              returned as [0,1) (the old RAND1, private here)
   ! Reference: src_VENUS/RAND1.f:59-106 (algorithm + constants verbatim; A = Knuth's
   !       6364136223846793005 in base 256, C = 0)
   real(8) function rng_step()
      integer :: ia(8), id(16), ip, i, j, k
      real(8), parameter :: bi = 3.90625d-3 ! 1/256 (old BI)
      ia = (/ 45, 127, 149, 76, 45, 244, 81, 88 /)   ! old DATA IA (base-256 multiplier)
      id(1:8) = 0                                     ! C = 0 (old DATA IC)
      id(9:16) = 0
      ! form IA*ISEED+IC: digit-wise convolution with carry (old :66-96)
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
      ! fold the 8 digits as base 256 into [0,1) (old :100-106)
      rng_step = dble(seed_dig(1))
      do i = 2, 8
         rng_step = dble(seed_dig(i)) + rng_step*bi
      end do
      rng_step = rng_step*bi
   end function rng_step

   !------------------------------------------------------------------
   ! rng_init(i_seed) - establish a new random stream from a decimal seed (base-256 conversion +
   !                    warm-up + filling the shuffle table)
   ! Reference: src_VENUS/RANDST.f:47-105 (verbatim arithmetic; the banner print :910 not carried)
   subroutine rng_init(i_seed)
      integer, intent(in) :: i_seed   ! decimal seed (0 <= i_seed <= 2^31-1; old scalar ISEED,
                                      ! the control%i_seed conversion mapping finalized in C4
      integer :: is, i
      ! 1. decimal -> 8 base-256 digits (old :52-56)
      is = i_seed
      do i = 1, 8
         seed_dig(i) = mod(is, dig_base)
         is = is/dig_base
      end do
      seed_dig(9) = 0
      ! 2. double digit-wise with carry (an even value; old :65-76)
      do i = 1, 8
         seed_dig(i) = seed_dig(i) + seed_dig(i)
      end do
      do i = 1, 8
         do while (seed_dig(i) >= dig_base)
            seed_dig(i) = seed_dig(i) - dig_base
            seed_dig(i + 1) = seed_dig(i + 1) + 1
         end do
      end do
      ! 3. +1 (an odd stream seed; old :81-90)
      seed_dig(1) = seed_dig(1) + 1
      do i = 1, 8
         do while (seed_dig(i) >= dig_base)
            seed_dig(i) = seed_dig(i) - dig_base
            seed_dig(i + 1) = seed_dig(i + 1) + 1
         end do
      end do
      ! 4. warm-up: fill the shuffle table TWICE in a row (same loop body, different results -
      !    the old two-loop note; old :92-105)
      do i = 1, 100
         tbl_shuf(i) = rng_step()
      end do
      do i = 1, 100
         tbl_shuf(i) = rng_step()
      end do
   end subroutine rng_init

   !------------------------------------------------------------------
   ! rng_u() - draw the next uniform deviate in [0,1) (shuffle-table draw + multiplicative congruential
   !           step to refill)
   ! Reference: src_VENUS/RAND0.f:19-23; src_VENUS/RAND1.f:66-106
   real(8) function rng_u()
      integer :: j
      j = int(99.0d0*tbl_shuf(100)) + 1     ! shuffle pointer (old RAND0.f:19)
      rng_u = tbl_shuf(100)                 ! the draw (old :20)
      tbl_shuf(100) = tbl_shuf(j)           ! (old :21)
      tbl_shuf(j) = rng_step()              ! refill through the congruential step (old :22)
   end function rng_u

   !------------------------------------------------------------------
   ! rng_gauss() - standard normal deviate N(0,1) (Box-Muller polar form + cached pair)
   ! Reference: src_TEST/interface_TEST.f90:112-136 (gasdev; uniform source switched to the
   !       single stream per D2 - one rng_u per component; the NaN debug print not carried)
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
   ! rng_gamma(order) - integer-order gamma deviate (order=1 degenerates to the unit exponential)
   ! Reference: src_VENUS/GAMA.f:7-13 (uniform source switched to the stream per D2)
   real(8) function rng_gamma(order)
      integer, intent(in) :: order    ! gamma distribution order (positive integer; 1 = exponential, old dummy IA)
      integer :: j
      real(8) :: dum
      dum = 1.0d0
      do j = 1, order
         dum = dum*rng_u()
      end do
      rng_gamma = -log(dum)
   end function rng_gamma
end module rng
