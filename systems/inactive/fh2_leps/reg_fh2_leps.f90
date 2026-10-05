!=====================================================================
! reg_fh2_leps.f90 - the F+H2 modified-LEPS container's drag-in file
!   (hands the fh2_pot / fh2_box material into the schedule's wiring
!   slots - zero system knowledge on either side of each binding; the
!   folder's presence in the build is the whole installation act)
! Design:
!   reg_fh2_leps_slots() is the single handover (slot 1 input / slot 2
!   regular / slot 3 initialization / slot 4 output; the schedule in
!   interface/container_sched owns the run order and the guard). The
!   modified-LEPS surface's legislated constants are staged in the
!   potential module itself, so the force slot needs no preparation.
!   Absent keys leave the optional judgment slots unbound.
! Units: none (registration carries code references only)
!=====================================================================
module reg_fh2_leps
   ! --- the schedule (the single handover target) ---
   use container_sched, only: container_sched_bind
   ! --- the wiring call sites ---
   use force_interface, only: container_bind_pes, container_bind_term
   use final_state,     only: container_bind_classify
   use input,           only: input_declare_keys, pull
   ! --- the container's own material ---
   use fh2_box, only: fh2_params_t, fh2_keys, fh2_load, fh2_pes, fh2_term, &
                      fh2_classify, fh2_term_on, fh2_class_on
   implicit none
   private
   public :: reg_fh2_leps_slots
contains
   !------------------------------------------------------------------
   ! reg_fh2_leps_slots() - the handover: fill the schedule's slots
   !------------------------------------------------------------------
   subroutine reg_fh2_leps_slots()
      call container_sched_bind(keys     =declare_keys, &
                                load     =load_params, &
                                force    =bind_force, &
                                term     =bind_term, &
                                classify =bind_classify)
   end subroutine reg_fh2_leps_slots

   !===== input wiring =======================================

   ! declare_keys() - declare the container key vocabulary (pure data;
   !               runs before the parse)
   subroutine declare_keys()
      call input_declare_keys(fh2_keys())
   end subroutine declare_keys

   ! load_params() - pull the buffered container keys into the parameter
   !               receiver, then load them into the module state
   subroutine load_params()
      type(fh2_params_t) :: p     ! parameter receiver (field defaults)
      call pull('R_TERM',  p%r_term, p%r_term_found)
      call pull('R_HF',    p%r_hf,   p%r_hf_found)
      call pull('D_FH',    p%surf(1), p%surf_found(1))
      call pull('R0_FH',   p%surf(2), p%surf_found(2))
      call pull('BETA_FH', p%surf(3), p%surf_found(3))
      call pull('D_HH',    p%surf(4), p%surf_found(4))
      call pull('R0_HH',   p%surf(5), p%surf_found(5))
      call pull('BETA_HH', p%surf(6), p%surf_found(6))
      call pull('DL0_HH',  p%surf(7), p%surf_found(7))
      call pull('DL0_FH',  p%surf(8), p%surf_found(8))
      call pull('AMOD',    p%surf(9), p%surf_found(9))
      call fh2_load(p)
   end subroutine load_params

   !===== regular wiring =====================================

   ! bind_force() - fill the force slot with the container's aggregate
   !                entry (the legislated constants stand by default;
   !                the slot needs no preparation)
   subroutine bind_force()
      call container_bind_pes(fh2_pes)
   end subroutine bind_force

   !===== initialization wiring ==============================

   ! bind_term() - bind the termination judgment the R_TERM key provides
   subroutine bind_term()
      if (fh2_term_on()) call container_bind_term(fh2_term)
   end subroutine bind_term

   !===== output wiring ======================================

   ! bind_classify() - bind the channel-classification judgment the R_HF key provides
   subroutine bind_classify()
      if (fh2_class_on()) call container_bind_classify(fh2_classify)
   end subroutine bind_classify
end module reg_fh2_leps
