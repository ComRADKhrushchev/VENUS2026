!=====================================================================
! reg_adiabatic.f90 - thin library-level registration module (methods/ seed):
!                     ONE method_reg call fills the whole 'adiabatic' package
!                     row - the reg file is the atomic registration unit (drag
!                     in a new file = drag in a new method; zero edits to
!                     existing files), aggregated into the generated registry
!                     module by the prebuild wildcard
! Contract:
!   reg_adiabatic_all() - register the 'adiabatic' method package (one call
!     fills the whole row - the row is the atomic registration unit): the
!     all-no-op degenerate package of the method table's built-in set (prop =
!     no-op - no amplitudes exist under adiabatic; match = no-op - occ stays
!     [1] constant; no sync member - the package declares neither ontic
!     container, and the null pointer is the legal 'no deliverable cache'
!     form). Passing ad_prop / ad_match as method_reg's procedure dummies
!     checks their conformance to the abstract interfaces elec_prop_i /
!     mqc_match_i at compile time (the member-signature freeze lands with the
!     interface phase). Executed once in the assembly phase - the driver binds
!     nothing by hand; this row is also the idle-chain method (the assembly
!     binds the degenerate package even when no container was dragged in, so
!     the body phases have something to run against).
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_adiabatic
   use elec_facade,   only: method_reg
   use elec_adiabatic, only: ad_prop, ad_match   ! the package members - passing them as
   ! method_reg's procedure dummies checks their conformance to the abstract interfaces
   ! elec_prop_i / mqc_match_i (the member-signature freeze lands with the interface
   ! phase)
   implicit none
   private
   public :: reg_adiabatic_all
contains
   !------------------------------------------------------------------
   ! reg_adiabatic_all() - register the 'adiabatic' method package (one call fills the
   !                       whole row - the row is the atomic registration unit)
   subroutine reg_adiabatic_all()
      call method_reg('adiabatic', '', ad_prop, ad_match)
   end subroutine reg_adiabatic_all
end module reg_adiabatic
