!=====================================================================
! elec_adiabatic.f90 - L3 method package (methods/ library seed, E14 spec §5): the adiabatic
!                      degenerate package - prop = no-op (no amplitudes exist under adiabatic),
!                      match = no-op (occ stays [1] constant), no sync member (neither ontic
!                      container is declared, spec §4) - the all-no-op degenerate point of the
!                      method table's built-in set
! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §5 (the built-in method
!       library: 'adiabatic' (all no-ops - the degenerate package) is the K2 seed, the other
!       three built-ins land in phase C; two exposure levels - level 1 applies as
!       config%electronic%method = 'adiabatic' binds the package with zero method code;
!       level 2 does not apply, the package exposes no reusable components); §3 (entries 2-3
!       semantics: propagation forms include no-op (LZSH, adiabatic); the matching entry
!       writes occ/p - here neither); §4 (adiabatic allocates neither a nor rho; occ = [1]
!       constant under adiabatic);
!       docs/design_spec_completeness.md §3.2-2 (two independent axes / well-defined
!       degenerate point / adiabatic: F = -grad V_alpha - consuming a single-surface
!       representation collapses the composite to a pure PES);
!       docs/final_completion_plan.md §5 N1 (class-B member degenerate point - ~20 lines
!       of explicit existence)
! Registration: method_reg('adiabatic', '', ad_prop, ad_match) - ONE call in reg_adiabatic.f90
!       (the thin reg file, spec §5 atomic registration unit); no optional sync - the null
!       pointer is the legal "no deliverable cache" form. The members conform to
!       elec_facade's abstract interfaces elec_prop_i / mqc_match_i - they take state_t
!       only (E14: ws is conveyed through the facade-internal channel, never a
!       public-signature dummy; conformance is machine-checked at the method_reg call
!       site - the E14 member-signature freeze lands in C2)
! Supersession note (E14): the old ad_dsdt/ad_force_rule adapter pair (the E9/E10-era
!       elec_proc_i(s, f, dt, ws) adaptation) is withdrawn - E14 splits the electronic
!       phases out of force composition (the four-entry contract), so the old adapter's
!       force side (identity force rule) is vacuous: force_eval composes the
!       single-surface PES kernel directly and the adiabatic package contributes nothing
!       to any entry. E9-note carried: the degenerate-equivalence verification criterion
!       (C2) now compares the four-entry orchestration (prop_step -> elec_prop ->
!       mqc_statistic -> rec_frame) against a pure-PES run - the two no-ops must leave
!       the trajectory bit-identical to the classical limit
!=====================================================================
module elec_adiabatic
   use state,       only: state_t
   use elec_facade, only: elec_prop_i, mqc_match_i   ! contract documentation only (the
   ! members below conform to these abstract interfaces; conformance is checked at the
   ! method_reg call site in reg_adiabatic.f90 - the E14 member-signature freeze lands
   ! in C2 together with the facade-internal ws conveyance wrapper)
   implicit none
   private
   public :: ad_prop, ad_match
contains
   !------------------------------------------------------------------
   ! ad_prop(sta, dt) - layer-2 member (elec_prop_i): adiabatic propagation = no-op - no
   !                    amplitudes exist under adiabatic (the method declares neither ontic
   !                    container, spec §4; ds/dt == 0 - the electronic state is frozen on
   !                    the current surface)
   ! Reference: docs/design_spec_completeness.md §3.2-2 (adiabatic degenerate point);
   !       spec §3 entry 2 (propagation forms: no-op (LZSH, adiabatic))
   subroutine ad_prop(sta, dt)
      type(state_t), intent(inout) :: sta  ! physical state (unchanged - no ontic container
                                           ! to advance, nothing read under this law)
      real(8), intent(in) :: dt            ! current step [10 fs] (ignored - no hop
                                           ! probability g_ij(dt) exists under adiabatic;
                                           ! E13 wording carried: the degenerate point is
                                           ! the identity along the whole channel)
      ! 1. Return (no-op; ignores ws and dt - the member reads nothing from the container:
      !    single surface, no coupling, no hops; the real dispatch target of elec_prop -
      !    the empty body IS the degenerate law)
   end subroutine ad_prop

   !------------------------------------------------------------------
   ! ad_match(sta) - layer-3 member (mqc_match_i): adiabatic matching = no-op - occ stays
   !                  [1] constant (the degenerate occupation, state_create seed; a hop has
   !                  no meaning on one surface), no rng use (the sole evolution-phase rng
   !                  consumption point consumes nothing here - the §4.1 ruling holds
   !                  trivially)
   ! Reference: spec §3 entry 3 + §4 (occ = [1] constant under adiabatic);
   !       docs/design_spec_completeness.md §3.2-2 (adiabatic degenerate point)
   subroutine ad_match(sta)
      type(state_t), intent(inout) :: sta  ! physical state (unchanged - occ already [1],
                                           ! nothing to sample; no momentum projection;
                                           ! n_hop untouched at 0)
      ! 1. Return (no-op - single surface, nothing to sample; the real dispatch
      !    target of mqc_statistic - occ stays [1] in the statistics holder,
      !    n_hop stays 0, the state passes through untouched)
   end subroutine ad_match
end module elec_adiabatic
