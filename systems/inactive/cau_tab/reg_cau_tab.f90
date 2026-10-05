!=====================================================================
! reg_cau_tab.f90 - the table-PES container's drag-in file (hands the
!   composed force entry and the folder's 'tdhf_hub' method package
!   into the schedule's wiring slots)
! Design:
!   reg_cau_tab_slots() is the single handover (slot 1 input / slot 2
!   regular / slot 3 initialization / slot 4 output; the schedule in
!   interface/container_sched owns the run order and the guard). The
!   method row registers from this folder (ELEC_METHOD = tdhf_hub
!   selects it; adiabatic keeps the pure ground-table reference face);
!   the TDHF arming and its recorder columns run only when the run
!   selected the package - the gate lives inside the slot bodies.
!   Same wiring shape as reg_cau_peri (2026-09-30); the only physics
!   difference is the PES source (TABLE_FILE vs the analytic/NN pair).
! Units: none (registration carries code references only)
!=====================================================================
module reg_cau_tab
   ! --- the schedule (the single handover target) ---
   use container_sched, only: container_sched_bind
   ! --- the wiring call sites ---
   use force_interface, only: container_bind_pes, container_bind_term
   use elec_interface,  only: method_reg
   use config,          only: electronic
   use input,           only: input_declare_keys, pull
   use recorder,        only: rec_reg_col
   ! --- the container's own material ---
   use cau_tab_box,  only: cau_tab_params_t, cau_tab_keys, cau_tab_load, &
                           cau_tab_term, cau_tab_term_on, &
                           cau_tab_pes_r, cau_tab_prop_r, &
                           cau_pars_start, cau_pars_norb, cau_pars_bw, cau_pars_telec, &
                           cau_pars_msbstp, cau_pars_nsbstp, cau_pars_pmode, cau_pars_pamp
   use cau_tab_tdhf, only: cau_tab_match, cau_tab_tdhf_setup
   use cau_tab_tdhf, only: cau_tab_sync, cau_tab_col_nup, cau_tab_col_ndn, &
                           cau_tab_col_mag, cau_tab_col_eelec, cau_tab_col_emf
   implicit none
   private
   public :: reg_cau_tab_slots
contains
   !------------------------------------------------------------------
   ! reg_cau_tab_slots() - the handover: fill the schedule's slots
   !------------------------------------------------------------------
   subroutine reg_cau_tab_slots()
      call container_sched_bind(keys    =declare_keys, &
                                load    =load_params, &
                                force   =bind_force, &
                                method  =register_method, &
                                term    =bind_term, &
                                arm     =arm_method, &
                                columns =register_columns)
   end subroutine reg_cau_tab_slots

   !===== input wiring =======================================

   ! declare_keys() - declare the container key vocabulary (pure data;
   !               runs before the parse)
   subroutine declare_keys()
      call input_declare_keys(cau_tab_keys())
   end subroutine declare_keys

   ! load_params() - pull the buffered container keys into the parameter
   !               receiver, then load them into the module state
   subroutine load_params()
      type(cau_tab_params_t) :: p   ! parameter receiver (field defaults)
      call pull('START_U0',     p%start)
      call pull('N_ORB',        p%norb)
      call pull('BAND_WIDTH',   p%bw)
      call pull('TELECTRO',     p%telec)
      call pull('MSBSTP',       p%msbstp)
      call pull('NSBSTP',       p%nsbstp)
      call pull('PERTURB_MODE', p%pmode)
      call pull('PERTURB_AMP',  p%pamp)
      call pull('K_SURF',       p%ksurf)
      call pull('K_GHOST',      p%kghost)
      call pull('Q0_GHOST',     p%q0ghost)
      call pull('TABLE_FILE',   p%table_file)
      call pull('Z_TERM',       p%z_term, p%z_term_found)
      call cau_tab_load(p)
   end subroutine load_params

   !===== regular wiring =====================================

   ! bind_force() - fill the force slot with the container's composed entry
   !                (the roster-order boundary wrapper: [C,S,Y] -> engine [S,Y,C])
   subroutine bind_force()
      call container_bind_pes(cau_tab_pes_r)
   end subroutine bind_force

   ! register_method() - register the folder's electronic method package
   !                 (one method row; ELEC_METHOD selects it at the seam)
   ! the propagation member rides the same roster-order wrapper
   ! (match/sync read only the rho state)
   subroutine register_method()
      call method_reg('tdhf_hub', cau_tab_prop_r, cau_tab_match, cau_tab_sync)
   end subroutine register_method

   !===== initialization wiring ==============================

   ! bind_term() - bind the termination judgment the Z_TERM key provides
   subroutine bind_term()
      if (cau_tab_term_on()) call container_bind_term(cau_tab_term)
   end subroutine bind_term

   ! arm_method() - arm the TDHF package when the run selected it
   subroutine arm_method()
      if (trim(electronic%method) == 'tdhf_hub') then
         call cau_tab_tdhf_setup(cau_pars_norb, cau_pars_bw, cau_pars_telec, &
                                 cau_pars_start, cau_pars_msbstp, cau_pars_nsbstp, &
                                 cau_pars_pmode, cau_pars_pamp)
      end if
   end subroutine arm_method

   !===== output wiring ======================================

   ! register_columns() - the container observables: obs_list-only recorder
   !                columns (the run opts in via OBS_LIST = n_imp_up,
   !                n_imp_dn, mag_mom, ...); registration order = column
   !                order. Selection-gated like the arm slot
   subroutine register_columns()
      if (trim(electronic%method) == 'tdhf_hub') then
         call rec_reg_col('n_imp_up', cau_tab_col_nup)
         call rec_reg_col('n_imp_dn', cau_tab_col_ndn)
         call rec_reg_col('mag_mom', cau_tab_col_mag)
         call rec_reg_col('e_elec_ev', cau_tab_col_eelec)
         call rec_reg_col('e_mf_ev', cau_tab_col_emf)
      end if
   end subroutine register_columns
end module reg_cau_tab
