!=====================================================================
! reg_reference.f90 - the reference container's drag-in file (REFERENCE
!   MATERIAL: hands the pure-math harmonic well of reference_box /
!   reference_pot into the schedule's wiring slots - zero system
!   knowledge on either side of each binding; the folder's presence in
!   the build is the whole installation act)
! Design:
!   reg_reference_slots() is the single handover: it fills the
!   schedule slots with this folder's wiring slots, grouped below as
!   slot 1 input / slot 2 regular / slot 3 initialization / slot 4
!   output. The schedule (interface/container_sched) owns the run
!   order and the initialization-assembly guard; every conditional
!   (key gates) lives inside the slot bodies. Absent keys leave the
!   optional judgment slots unbound - the safety-net-only,
!   no-classification run.
! Units: none (registration carries code references only)
!=====================================================================
module reg_reference
   ! --- the schedule (the single handover target) ---
   use container_sched, only: container_sched_bind
   ! --- the wiring call sites ---
   use force_interface, only: container_bind_pes, container_bind_term
   use final_state,     only: container_bind_classify
   use input,           only: input_declare_keys, pull
   ! --- the container's own material ---
   use reference_box, only: reference_params_t, reference_keys, reference_load, &
                            reference_pes, reference_term, reference_classify, &
                            reference_term_on, reference_class_on, reference_k_default
   use reference_pot, only: harm_set_constant
   implicit none
   private
   public :: reg_reference_slots
contains
   !------------------------------------------------------------------
   ! reg_reference_slots() - the handover: fill the schedule's slots
   !------------------------------------------------------------------
   subroutine reg_reference_slots()
      call container_sched_bind(keys     =declare_keys, &
                                load     =load_params, &
                                force    =bind_force, &
                                term     =bind_term, &
                                classify =bind_classify)
   end subroutine reg_reference_slots

   !===== input wiring =======================================

   ! declare_keys() - declare the container key vocabulary (pure data;
   !               runs before the parse)
   subroutine declare_keys()
      call input_declare_keys(reference_keys())
   end subroutine declare_keys

   ! load_params() - pull the buffered container keys into the parameter
   !               receiver, then load them into the module state
   subroutine load_params()
      type(reference_params_t) :: p  ! parameter receiver (field defaults)
      call pull('K_HARM',  p%k_harm)
      call pull('Q_TERM',  p%q_term,  p%q_term_found)
      call pull('E_SPLIT', p%e_split, p%e_split_found)
      call reference_load(p)
   end subroutine load_params

   !===== regular wiring =====================================

   ! bind_force() - stage the well's default constant into the
   !                potential module (the bound PES is computable the
   !                moment the slot is filled), then fill the force
   !                slot with the container's aggregate entry
   subroutine bind_force()
      call harm_set_constant(reference_k_default())
      call container_bind_pes(reference_pes)
   end subroutine bind_force

   !===== initialization wiring ==============================

   ! bind_term() - bind the termination judgment the Q_TERM key provides
   subroutine bind_term()
      if (reference_term_on()) call container_bind_term(reference_term)
   end subroutine bind_term

   !===== output wiring ======================================

   ! bind_classify() - bind the classification judgment the E_SPLIT key provides
   subroutine bind_classify()
      if (reference_class_on()) call container_bind_classify(reference_classify)
   end subroutine bind_classify
end module reg_reference
