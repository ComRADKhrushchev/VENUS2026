!=====================================================================
! reg_h2ag111.f90 - the H2/Ag(111) container's drag-in file (hands the
!   h2ag111_pot / h2ag111_box material into the schedule's wiring
!   slots - zero system knowledge on either side of each binding; the
!   folder's presence in the build is the whole installation act)
! Design:
!   reg_h2ag111_slots() is the single handover (slot 1 input / slot 2
!   regular / slot 3 initialization / slot 4 output; the schedule in
!   interface/container_sched owns the run order and the guard). The
!   network is NOT computable at registration: its one-time
!   initialization is run-material dependent and happens inside
!   h2ag111_load at the seam; force_eval's first call comes after the
!   seam, so the bound PES is computable by then. The NN PES carries
!   no asymptotic tail, so the isolated-H2 spectrum cannot be probed
!   (COMPUTE named-aborts at the decay gate or samples a non-literature
!   curvature) - the literature gas-phase set belongs in the MANUAL
!   spectrum-table file (SPECTRUM_SOURCE=MANUAL + SPECTRUM_FILE; the
!   container export slot is withdrawn, 2026-10-05 spectrum-source
!   legislation).
! Units: none (registration carries code references only)
!=====================================================================
module reg_h2ag111
   ! --- the schedule (the single handover target) ---
   use container_sched, only: container_sched_bind
   ! --- the wiring call sites ---
   use force_interface, only: container_bind_pes, container_bind_term
   use final_state,     only: container_bind_classify
   use input,           only: input_declare_keys, pull
   ! --- the container's own material ---
   use h2ag111_box, only: h2ag111_params_t, h2ag111_keys, h2ag111_load, h2ag111_pes, &
                         h2ag111_term, h2ag111_classify, h2ag111_term_on, h2ag111_class_on
   implicit none
   private
   public :: reg_h2ag111_slots
contains
   !------------------------------------------------------------------
   ! reg_h2ag111_slots() - the handover: fill the schedule's slots
   !------------------------------------------------------------------
   subroutine reg_h2ag111_slots()
      call container_sched_bind(keys     =declare_keys, &
                                load     =load_params, &
                                force    =bind_force, &
                                term     =bind_term, &
                                classify =bind_classify)
   end subroutine reg_h2ag111_slots

   !===== input wiring =======================================

   ! declare_keys() - declare the container key vocabulary (pure data;
   !               runs before the parse)
   subroutine declare_keys()
      call input_declare_keys(h2ag111_keys())
   end subroutine declare_keys

   ! load_params() - pull the buffered container keys into the parameter
   !               receiver, then load them into the module state (the
   !               load runs the network's one-time initialization)
   subroutine load_params()
      type(h2ag111_params_t) :: p     ! parameter receiver (field defaults)
      call pull('Z_TERM', p%z_term, p%z_term_found)
      call pull('R_H2',   p%r_h2,   p%r_h2_found)
      call h2ag111_load(p)
   end subroutine load_params

   !===== regular wiring =====================================

   ! bind_force() - fill the force slot with the container's aggregate
   !                entry (computable once the seam's load has run;
   !                force_eval's first call comes after the seam)
   subroutine bind_force()
      call container_bind_pes(h2ag111_pes)
   end subroutine bind_force

   !===== initialization wiring ==============================

   ! bind_term() - bind the termination judgment the Z_TERM key provides
   subroutine bind_term()
      if (h2ag111_term_on()) call container_bind_term(h2ag111_term)
   end subroutine bind_term

   !===== output wiring ======================================

   ! bind_classify() - bind the channel-classification judgment the R_H2 key provides
   subroutine bind_classify()
      if (h2ag111_class_on()) call container_bind_classify(h2ag111_classify)
   end subroutine bind_classify
end module reg_h2ag111
