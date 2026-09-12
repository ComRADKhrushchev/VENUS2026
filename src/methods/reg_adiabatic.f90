!=====================================================================
! reg_adiabatic.f90 - thin library-level registration module (methods/ seed, E14 spec §5):
!                     ONE method_reg call fills the whole 'adiabatic' package row - the reg
!                     file is the atomic registration unit (drag in a new file = drag in a
!                     new method; zero edits to existing files)
! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §5 (registration data
!       flow, library level: each method file is registered through its own thin reg_*.f90
!       in methods/; prebuild aggregates every reg_*.f90 under methods/ and systems/ into
!       registry_gen.f90 - K3 scope; the K2 skeleton wires reg_adiabatic_all from the
!       driver pseudocode transitionally); §5 built-in set ('adiabatic' - all no-ops, the
!       degenerate package; the K2 seed, reg_adiabatic.f90 is the sole reg file of the
!       library seed - the dist_* reg files land at C5 with their sample routines, per the
!       K2 E3 slimming ruling)
!=====================================================================
module reg_adiabatic
   use elec_facade,   only: method_reg
   use elec_adiabatic, only: ad_prop, ad_match   ! the package members - passing them as
   ! method_reg's procedure dummies checks their conformance to the abstract interfaces
   ! elec_prop_i / mqc_match_i (the E14 member-signature freeze lands in C2)
   implicit none
   private
   public :: reg_adiabatic_all
contains
   !------------------------------------------------------------------
   ! reg_adiabatic_all() - register the 'adiabatic' method package (one call fills the
   !                       whole row - the row is the atomic registration unit, spec §5)
   ! Reference: spec §5 (level-1 registration: 'adiabatic' in the built-in set; one
   !       registration vouches for the whole package - the pairing legality is this
   !       author's knowledge, the facade checks only name-table hygiene, E13 carried)
   subroutine reg_adiabatic_all()
      ! Pseudocode:
      ! 1. call method_reg('adiabatic', '', ad_prop, ad_match) (no optional sync - adiabatic
      !    allocates neither ontic container, nothing to synchronize; the null pointer is
      !    the legal "no deliverable cache" form, spec §4. Executed once in the assembly
      !    phase - transitional wiring from the driver pseudocode until the K3 prebuild
      !    aggregates this file into registry_gen, then the driver carries zero
      !    registration lines, E11)
   end subroutine reg_adiabatic_all
end module reg_adiabatic
