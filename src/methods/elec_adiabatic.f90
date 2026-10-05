!=====================================================================
! elec_adiabatic.f90 - the adiabatic degenerate electronic-method package
! Design:
!   The all-no-op built-in: prop = no-op (no amplitudes exist under
!   adiabatic), match = no-op (occ stays [1] constant), no sync member
!   (the null pointer is the legal 'no deliverable cache' form); force
!   composition reads the single-surface PES directly, so this package
!   contributes nothing to any entry.
!   Members conform to the interface's abstract interfaces, checked at
!   the method_reg call site in reg_adiabatic.f90.
! Units: dt in 10 fs internal units.
!=====================================================================
module elec_adiabatic
   use state,       only: state_t
   use elec_interface, only: elec_prop_i, mqc_match_i   ! contract documentation only (the
   ! members below conform; conformance is checked at the method_reg call site)
   implicit none
   private
   public :: ad_prop, ad_match
contains
   !------------------------------------------------------------------
   ! ad_prop(sta, dt) - adiabatic propagation = no-op (no amplitudes
   !                    exist; the electronic state is frozen on the
   !                    current surface)
   subroutine ad_prop(sta, dt)
      type(state_t), intent(inout) :: sta  ! physical state (unchanged - no ontic
                                           ! container to advance)
      real(8), intent(in) :: dt            ! current step [10 fs] (ignored - no hop
                                           ! probability exists under adiabatic)
      ! 1. No-op: single surface, no coupling, no hops - the empty body IS the law
   end subroutine ad_prop

   !------------------------------------------------------------------
   ! ad_match(sta) - adiabatic matching = no-op (occ stays [1] constant,
   !                  no rng use - a hop has no meaning on one surface)
   subroutine ad_match(sta)
      type(state_t), intent(inout) :: sta  ! physical state (unchanged - occ already [1],
                                           ! nothing to sample; n_hop stays 0)
      ! 1. No-op: the state passes through untouched
   end subroutine ad_match
end module elec_adiabatic
