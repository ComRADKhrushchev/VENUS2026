!=====================================================================
! reg_reference.f90 - the reference container's drag-in file (REFERENCE
!                     MATERIAL: binds the pure-math harmonic well of
!                     reference_box / reference_pot into the container
!                     slots - zero system knowledge on either side of the
!                     binding; the folder's presence in the build is the
!                     whole installation act)
! Contract:
!   reg_reference_all() - the class-A binding, unconditional (one
!     container per run; the drag-in does not consult anything): the
!     isotropic well in scalar mode, delivery declared 'scalar'. The
!     member init runs BEFORE the bind: the well's default constant is
!     staged into the potential module, so the bound PES is computable
!     the moment the slot is filled (the degenerate posture without a
!     parameter file runs on the default constant)
!   reg_reference_init() - the assembly seam: load the parameter file
!     reference_pars.txt from the run directory, then bind the OPTIONAL
!     slots the file provides - Q_TERM -> the termination judgment,
!     E_SPLIT -> the classification judgment. Absent keys leave the slot
!     unbound (the safety-net-only, no-classification run). Guard: no
!     staged input (system_dir empty) -> skip, nothing loaded, nothing
!     bound
! Units: none (registration carries code references and the parameter
!   file path only).
!=====================================================================
module reg_reference
   use force_facade, only: container_bind_pes, container_bind_term
   use final_state,  only: container_bind_classify
   use config,       only: reactants
   use reference_box, only: reference_load, reference_pes, reference_term, &
                            reference_classify, reference_term_on, reference_class_on, &
                            reference_k_default
   use reference_pot, only: harm_set_constant
   implicit none
   private
   public :: reg_reference_all, reg_reference_init
contains
   !------------------------------------------------------------------
   ! reg_reference_all() - member init (stage the default constant into the
   !                       potential module), then bind the container
   !                       aggregate class-A entry with its delivery declared
   subroutine reg_reference_all()
      call harm_set_constant(reference_k_default())
      call container_bind_pes(reference_pes, 'scalar')
   end subroutine reg_reference_all

   !------------------------------------------------------------------
   ! reg_reference_init() - the assembly seam: parameters first, then
   !                       the optional-slot bindings the parameters provide
   subroutine reg_reference_init()
      if (len_trim(reactants%system_dir) == 0) return      ! no input staged - skip
      call reference_load()
      if (reference_term_on()) call container_bind_term(reference_term)
      if (reference_class_on()) call container_bind_classify(reference_classify)
   end subroutine reg_reference_init
end module reg_reference
