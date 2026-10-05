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
!   optional judgment slots unbound.
! Units: none (registration carries code references only)
!=====================================================================
module reg_h3_bkmp2
   ! --- the schedule (the single handover target) ---
   use container_sched, only: container_sched_bind
   ! --- the wiring call sites ---
   use force_interface, only: container_bind_pes, container_bind_term
   use final_state,     only: container_bind_classify
   use input,           only: input_declare_keys, pull
   use config_atoms,    only: list_atoms
   use spectrum,        only: wvn_to_w, stretch_mode
   use spectrum_export, only: spectrum_export_put, spectrum_export_put_diatomic
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
                                export   =export_spectrum, &
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

   ! export_spectrum() - registration-time spectrum export: export the gas-phase
   !                w_e of the H2 target fragment (Herzberg constants of the
   !                record) plus the diatomic constants for the EBK member.
   !                The BKMP2 PES has a proper asymptote (unlike surface NN),
   !                so the generic derivation could also work for the harmonic
   !                w, but the constants below are exact and the EBK channel
   !                is export-only either way
   subroutine export_spectrum()
      real(8), parameter :: h2_we_wvn = 4401.21d0   ! gas-phase H2 w_e [cm^-1]
      real(8), parameter :: h2_wx = 121.33d0, h2_be = 60.853d0
      real(8), parameter :: d2_we_wvn = 3115.5d0    ! gas-phase D2 w_e [cm^-1]
      real(8), parameter :: d2_wx = 61.82d0,  d2_be = 30.429d0
      real(8) :: w(1), c(6, 1)
      integer :: i
      if (.not. allocated(list_atoms%frag)) return
      do i = 1, size(list_atoms%frag)
         if (list_atoms%frag(i)%nat == 2) then
            if (is_sym(i, 1, 'h') .and. is_sym(i, 2, 'h')) then
               w(1) = wvn_to_w(h2_we_wvn)
               call stretch_mode(list_atoms%frag(i)%qz, list_atoms%mass(list_atoms%frag(i)%list), c(:, 1))
               call spectrum_export_put([character(len=4)::'H', 'H'], w, c)
               call spectrum_export_put_diatomic([character(len=4)::'H', 'H'], &
                                                 h2_we_wvn, h2_wx, h2_be)
               return
            end if
            if (is_sym(i, 1, 'd') .and. is_sym(i, 2, 'd')) then
               w(1) = wvn_to_w(d2_we_wvn)
               call stretch_mode(list_atoms%frag(i)%qz, list_atoms%mass(list_atoms%frag(i)%list), c(:, 1))
               call spectrum_export_put([character(len=4)::'D', 'D'], w, c)
               call spectrum_export_put_diatomic([character(len=4)::'D', 'D'], &
                                                 d2_we_wvn, d2_wx, d2_be)
               return
            end if
         end if
      end do
   end subroutine export_spectrum

   ! is_sym(fr, j, s) - list_atoms symbol of fragment fr's j-th atom folds to s
   pure logical function is_sym(fr, j, s)
      character(len=*), intent(in) :: s
      integer, intent(in) :: fr, j
      character(len=4) :: sym
      integer :: k, cc
      is_sym = .false.
      sym = list_atoms%symb(list_atoms%frag(fr)%list(j))
      do k = 1, len_trim(sym)
         cc = iachar(sym(k:k))
         if (cc >= iachar('A') .and. cc <= iachar('Z')) sym(k:k) = char(cc + 32)
      end do
      is_sym = (trim(sym) == s)
   end function is_sym
end module reg_h3_bkmp2
