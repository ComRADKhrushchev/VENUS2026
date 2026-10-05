!=====================================================================
! inc_single.f90 - the trivial single-reactant incident-channel member
! Design:
!   A single-fragment gas-phase run has no collision partner: the
!   per-fragment members write the whole initial state and the incident
!   channel owns NOTHING. This member formalizes that degenerate case so
!   the sampler dispatches the (paradigm-resolved) incident member
!   uniformly - draw consumes no RNG, realize touches no state (both are
!   bitwise no-ops; the historical skip-guard produced the identical
!   stream). The E_REL / B_MAX / R_SEP keys remain required program keys
!   (placeholders in this paradigm - they carry no prescription here).
!=====================================================================
module inc_single
   use state, only: state_t
   implicit none
   private
   public :: inc_single_draw, inc_single_realize, inc_single_sample
contains
   !------------------------------------------------------------------
   ! inc_single_draw() - nothing to draw (no collision prescription)
   !------------------------------------------------------------------
   subroutine inc_single_draw()
   end subroutine inc_single_draw

   !------------------------------------------------------------------
   ! inc_single_realize(sta) - nothing to assemble (the per-fragment
   !                 members own the whole initial state)
   !------------------------------------------------------------------
   subroutine inc_single_realize(sta)
      type(state_t), intent(inout) :: sta  ! untouched by contract (owned by nobody here)
   end subroutine inc_single_realize

   !------------------------------------------------------------------
   ! inc_single_sample(sta) - fused convenience (the check-program entry)
   !------------------------------------------------------------------
   subroutine inc_single_sample(sta)
      type(state_t), intent(inout) :: sta
      call inc_single_draw()
      call inc_single_realize(sta)
   end subroutine inc_single_sample
end module inc_single
