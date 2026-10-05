!=====================================================================
! reg_h3_bkmp2.f90 - the H+H2 BKMP2 container's drag-in file (hands
!   the bkmp2_pot / bkmp2_box material into the schedule's wiring
!   slots - zero system knowledge on either side of each binding; the
!   folder's presence in the build is the whole installation act)
! Design:
!   reg_h3_bkmp2_slots() is the single handover (slot 1 input / slot 2
!   regular / slot 3 initialization / slot 4 output; the schedule in
!   interface/container_sched owns the run order and the guard). The
!   archived fit constants are staged in the potential module itself,
!   so the force slot needs no preparation. Absent keys leave the
!   optional judgment slots unbound. The BKMP2 PES has a proper
!   asymptote, so the internal derivation (SPECTRUM_SOURCE=COMPUTE,
!   the default) serves the harmonic mode table; the exact literature
!   constants (and the EBK diatomic set, which has no derivation
!   either way) belong in the MANUAL spectrum-table file (the
!   container export slot is withdrawn, 2026-10-05 spectrum-source
!   legislation).
! Units: none (registration carries code references only)
!=====================================================================
module reg_h3_bkmp2
   ! --- the schedule (the single handover target) ---
   use container_sched, only: container_sched_bind
   ! --- the wiring call sites ---
   use force_interface, only: container_bind_pes, container_bind_term
   use final_state,     only: container_bind_classify
   use input,           only: input_declare_keys, pull
   ! --- the container's own material ---
   use bkmp2_box, only: bkmp2_params_t, bkmp2_keys, bkmp2_load, bkmp2_pes, &
                       bkmp2_term, bkmp2_classify, bkmp2_term_on, bkmp2_class_on
   implicit none
   private
   public :: reg_h3_bkmp2_slots
contains
   !------------------------------------------------------------------
   ! reg_h3_bkmp2_slots() - the handover: fill the schedule's slots
   !------------------------------------------------------------------
   subroutine reg_h3_bkmp2_slots()
      call container_sched_bind(keys     =declare_keys, &
                                load     =load_params, &
                                force    =bind_force, &
                                term     =bind_term, &
                                classify =bind_classify)
   end subroutine reg_h3_bkmp2_slots

   !===== input wiring =======================================

   ! declare_keys() - declare the container key vocabulary (pure data;
   !               runs before the parse)
   subroutine declare_keys()
      call input_declare_keys(bkmp2_keys())
   end subroutine declare_keys

   ! load_params() - pull the buffered container keys into the parameter
   !               receiver, then load them into the module state
   subroutine load_params()
      type(bkmp2_params_t) :: p     ! parameter receiver (field defaults)
      call pull('R_TERM', p%r_term, p%r_term_found)
      call pull('R_H2',   p%r_h2,   p%r_h2_found)
      call bkmp2_load(p)
   end subroutine load_params

   !===== regular wiring =====================================

   ! bind_force() - fill the force slot with the container's aggregate
   !                entry (the archived fit constants stand by default;
   !                the slot needs no preparation)
   subroutine bind_force()
      call container_bind_pes(bkmp2_pes)
   end subroutine bind_force

   !===== initialization wiring ==============================

   ! bind_term() - bind the termination judgment the R_TERM key provides
   subroutine bind_term()
      if (bkmp2_term_on()) call container_bind_term(bkmp2_term)
   end subroutine bind_term

   !===== output wiring ======================================

   ! bind_classify() - bind the channel-classification judgment the R_H2 key provides
   subroutine bind_classify()
      if (bkmp2_class_on()) call container_bind_classify(bkmp2_classify)
   end subroutine bind_classify
end module reg_h3_bkmp2
