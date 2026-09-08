!=====================================================================
! elec_adiabatic.f90 - L3 class-B degenerate member: adiabatic electron law (ds/dt == 0, force rule =
!                      identity -
!                      under a single-surface representation the composite collapses to a pure
!                      PES; the interface is well-defined at the degenerate point)
! Reference: docs/design_spec_completeness.md §3.2-2 (two independent axes / well-defined
!       degenerate point /
!       adiabatic: F = -grad V_alpha - consuming a single-surface representation collapses the
!       composite to a pure PES);
!       docs/final_completion_plan.md §5 N1 (class-B member degenerate point - ~20 lines of
!       explicit existence)
! Registration: facade_reg_elec(name, proc) is executed by the assembly phase (driver); this
! degenerate member's two
! routines are adapted to the class-B interface slot in C2 (E9 extended the class-B signature
! with dt, E10 further with ws -
! elec_proc_i(s, f, dt, ws); the adapter ignores dt and ws: the adiabatic law has no hop
! decision and
! does not consume multi-surface constructs - the degenerate point is the identity along the
! whole channel)
!=====================================================================
module elec_adiabatic
   use state, only: elec_t
   implicit none
   private
   public :: ad_dsdt, ad_force_rule
contains
   !------------------------------------------------------------------
   ! ad_dsdt(s) - electronic-state propagation rule: adiabatic degenerate point ds/dt == 0 (the
   ! electronic state is frozen on the current surface)
   ! Reference: docs/design_spec_completeness.md §3.2-2 (adiabatic is a degenerate point)
   subroutine ad_dsdt(s)
      type(elec_t), intent(inout) :: s   ! electronic state (unchanged under the adiabatic law -
      ! single surface, no coupling, no hops)
      ! Pseudocode:
      ! 1. Identity degenerate point: s keeps the current surface i_surf unchanged (ds/dt == 0 -
      ! the propagator advances s not at all)
   end subroutine ad_dsdt

   !------------------------------------------------------------------
   ! ad_force_rule(f) - electron-law force modification: adiabatic degenerate point identity (the summed
   ! single-surface forces are the final nuclear force)
   ! Reference: docs/design_spec_completeness.md §3.2-2 (adiabatic: F = -grad V_alpha)
   subroutine ad_force_rule(f)
      real(8), intent(inout) :: f(:)   ! nuclear forces [internal energy/Å] (identity in and out)
      ! Pseudocode:
      ! 1. Identity degenerate point: f passes through unchanged (no surface-weight averaging,
      ! no surface-switch modification)
   end subroutine ad_force_rule
end module elec_adiabatic
