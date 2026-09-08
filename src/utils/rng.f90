!=====================================================================
! rng.f90 - L5 instrument: uniform random stream (Knuth multiplicative congruential + shuffle) and
!           distribution transforms (normal / gamma)
! Reference: src_VENUS/RANDST.f:47-105 (stream initialization); src_VENUS/RAND0.f:19-23
! (shuffle-table draws);
!       src_VENUS/RAND1.f:59-106 (multiplicative congruential stepping, 8 base-256 digits);
!       src_TEST/interface_TEST.f90:112-136 (gasdev normal pairs);
!       src_VENUS/GAMA.f:7-13 (gamma/exponential)
!=====================================================================
module rng
   implicit none
   private
   public :: rng_init, rng_u, rng_gauss, rng_gamma

   ! Stream state (internalized in the module - K1 plan §2: the old RANLST/ISEED3/IBFCTR become
   ! module-private)
   integer, parameter :: dig_base = 256     ! digit base 2^8 (old IBFCTR)
   integer :: seed_dig(8) = 0               ! stream seed as 8 base-256 digits (old ISEED3)
   real(8) :: tbl_shuf(100) = 0.0d0         ! shuffle table of 100 entries (old RANLST)
   ! Box-Muller cached-pair state (the iset/gset of old gasdev):
   integer :: i_cache = 0                   ! cache flag (0 = none, 1 = cached; old iset)
   real(8) :: g_cache = 0.0d0               ! cached second normal deviate (old gset)
contains
   !------------------------------------------------------------------
   ! rng_init(i_seed) - establish a new random stream from a decimal seed (base-256 conversion +
   ! warm-up + filling the shuffle table)
   ! Reference: src_VENUS/RANDST.f:47-105
   subroutine rng_init(i_seed)
      integer, intent(in) :: i_seed   ! decimal seed (0 <= i_seed <= 2^31-1; old scalar ISEED,
                                      ! the control%i_seed conversion mapping finalized in C4
      ! Pseudocode:
      ! 1. Convert i_seed digit by digit (divide by dig_base) into 8 digits stored in seed_dig (old :52-56)
      ! 2. Double seed_dig and normalize digit-wise (subtract and carry up when >= dig_base), ensuring an
      !    even value (old :65-76)
      ! 3. Increment seed_dig(1) and normalize, ensuring an odd stream seed (old :81-90)
      ! 4. Draw 100 numbers twice in a row to fill the shuffle table tbl_shuf (generator warm-up; both
      !    rounds share the loop body, old :92-102)
      !    but different results, old :100-105)
   end subroutine rng_init

   !------------------------------------------------------------------
   ! rng_u() - draw the next uniform deviate in [0,1) (shuffle-table draw + multiplicative congruential
   ! step to refill)
   ! Reference: src_VENUS/RAND0.f:19-23; src_VENUS/RAND1.f:66-106
   real(8) function rng_u()
      ! Result: a uniform random number in [0,1) [dimensionless]
      ! Pseudocode:
      ! 1. j = int(99*tbl_shuf(100)) + 1 (shuffle pointer, old RAND0.f:19)
      ! 2. Result = tbl_shuf(100); tbl_shuf(100) = tbl_shuf(j) (old :20-21)
      ! 3. tbl_shuf(j) = the multiplicative-congruential step value (old :22 calls RAND1):
      ! 3a. Seed as 8 digits and multiplier as 8 digits (the base-256 decomposition of Knuth's multiplier
      ! 6364136223846793005, old RAND1.f:59)
      !     base-256 form (/45,127,149,76,45,244,81,88/), old RAND1.f:59)
      !     digit-wise convolution multiply-add, normalizing and carrying at dig_base
      ! 3b. Write the low 8 digits back into seed_dig (old :96-98)
      ! 3c. Fold the 8 digits as base 256 into a float, then multiply by 1/dig_base to get [0,1) (old :102-106)
   end function rng_u

   !------------------------------------------------------------------
   ! rng_gauss() - standard normal deviate N(0,1) (Box-Muller polar form + cached pair)
   ! Reference: src_TEST/interface_TEST.f90:112-136 (gasdev)
   real(8) function rng_gauss()
      ! Result: a standard normal deviate [dimensionless]
      ! Pseudocode:
      ! 1. Cache flag set: result = g_cache, clear the flag and return (old :131-134)
      ! 2. Otherwise rejection sampling: draw pairs v = 2*rng_u()-1 until r2 = v1^2+v2^2 falls in (0,1)
      !    (uniform points inside the unit circle, old :118-125)
      ! 3. fac = sqrt(-2*ln(r2)/r2) (old :127)
      ! 4. Set g_cache = v1*fac with the flag; result = v2*fac (a cached pair, old :128-130)
   end function rng_gauss

   !------------------------------------------------------------------
   ! rng_gamma(order) - integer-order gamma deviate (order=1 degenerates to the unit exponential)
   ! Reference: src_VENUS/GAMA.f:7-13
   real(8) function rng_gamma(order)
      integer, intent(in) :: order    ! gamma distribution order (positive integer; 1 = exponential, old dummy IA)
      ! Result: a Gamma(order,1) deviate [dimensionless]
      ! Pseudocode:
      ! 1. dum = product of order rng_u() draws (old :8-11; the old version mixed in the intrinsic
      !    the new version draws uniformly through the rng_u stream, keeping the stream reproducible)
      ! 2. Result = -log(dum) (old :12)
   end function rng_gamma
end module rng
