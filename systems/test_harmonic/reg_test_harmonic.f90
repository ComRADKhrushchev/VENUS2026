!=====================================================================
! reg_test_harmonic.f90 - the verification container's drag-in file
!                         (VERIFICATION MATERIAL: binds the pure-math
!                         harmonic well of harmonic_box into the container
!                         slots - zero system knowledge on either side of
!                         the binding; the folder's presence in the build
!                         is the whole installation act)
! Contract:
!   reg_test_harmonic_all() - the class-A binding, unconditional (one
!     container per run; the drag-in does not consult anything): the
!     isotropic well in scalar mode - the degenerate single-surface
!     delivery set, no declaration argument
!   reg_test_harmonic_init() - the assembly seam: load the parameter file
!     harmonic_pars.txt from the run directory, then bind the OPTIONAL
!     slots the file provides - Q_TERM -> the termination judgment,
!     E_SPLIT -> the classification judgment. Absent keys leave the slot
!     unbound (the safety-net-only, no-classification run). Guard: no
!     staged input (system_dir empty) -> skip, nothing loaded, nothing
!     bound
! Units: none (registration carries code references and the parameter
!   file path only).
!=====================================================================
module reg_test_harmonic
   use force_facade, only: container_bind_pes, container_bind_term
   use final_state,  only: container_bind_classify
   use config,       only: reactants
   use harmonic_box, only: harmonic_load, harmonic_pes, harmonic_term, &
                           harmonic_classify, harmonic_term_on, harmonic_class_on
   implicit none
   private
   public :: reg_test_harmonic_all, reg_test_harmonic_init
contains
   !------------------------------------------------------------------
   ! reg_test_harmonic_all() - bind the container aggregate class-A entry
   subroutine reg_test_harmonic_all()
      call container_bind_pes(harmonic_pes)
   end subroutine reg_test_harmonic_all

   !------------------------------------------------------------------
   ! reg_test_harmonic_init() - the assembly seam: parameters first, then
   !                       the optional-slot bindings the parameters provide
   subroutine reg_test_harmonic_init()
      if (len_trim(reactants%system_dir) == 0) return      ! no input staged - skip
      call harmonic_load()
      if (harmonic_term_on()) call container_bind_term(harmonic_term)
      if (harmonic_class_on()) call container_bind_classify(harmonic_classify)
   end subroutine reg_test_harmonic_init
end module reg_test_harmonic
