!=====================================================================
! cau_tab_wires.f90 - the Tier 1 reserved-name wrappers (the folder's
!   auto-weave contract: prebuild weaves these names into registry_gen
!   when the folder sits reg-less at systems/ top level; each delegates
!   to the public API the check programs and probe hooks keep using).
!   Replaces reg_cau_tab.f90 (the reg form retired 2026-10-14,
!   docs/plans/2026-10-14-reg-retirement-observables.md).
!   A dedicated wires file (not in cau_tab_box) keeps the box free of
!   use input/recorder/config - the five criterion recipes compile the
!   box standalone without those modules.
! Units: none
!=====================================================================
module cau_tab_wires
   use state,          only: state_t, elec_t
   use config,         only: electronic
   use input,          only: input_declare_keys, pull
   use recorder,       only: rec_reg_col
   use cau_tab_box,    only: cau_tab_params_t, cau_tab_keys, cau_tab_load, &
                             cau_tab_term, cau_tab_term_on, &
                             cau_tab_pes_r, cau_tab_prop_r, &
                             cau_pars_start, cau_pars_norb, cau_pars_bw, cau_pars_telec, &
                             cau_pars_msbstp, cau_pars_nsbstp, cau_pars_pmode, cau_pars_pamp
   use cau_tab_tdhf,   only: cau_tab_match, cau_tab_tdhf_setup, cau_tab_sync, &
                             cau_tab_col_nup, cau_tab_col_ndn, cau_tab_col_mag, &
                             cau_tab_col_eelec, cau_tab_col_emf, &
                             cau_tab_col_eimp, cau_tab_col_ebath, cau_tab_col_ehyb
   implicit none
   private
   public :: pes, keys, load, term
   public :: prop, match, sync, method_name, init, observables
contains
   !------------------------------------------------------------------
   ! the pure-PES four
   !------------------------------------------------------------------
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      call cau_tab_pes_r(q, mass, g, v)   ! the roster-order boundary wrapper
   end subroutine pes

   function keys() result(kk)
      character(len=16) :: kk(13)
      kk = cau_tab_keys()
   end function keys

   subroutine load()
      type(cau_tab_params_t) :: p
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
   end subroutine load

   logical function term(sta)
      type(state_t), intent(in) :: sta
      if (.not. cau_tab_term_on()) then   ! the reg-era bind gate, carried
         term = .false.                   ! by the wrapper (the body does
         return                           ! not self-gate; an unbound run
      end if                              ! must answer .false.)
      term = cau_tab_term(sta)
   end function term

   !------------------------------------------------------------------
   ! the method package (tdhf_hub)
   !------------------------------------------------------------------
   subroutine prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
      call cau_tab_prop_r(sta, dt)
   end subroutine prop

   subroutine match(sta)
      type(state_t), intent(inout) :: sta
      call cau_tab_match(sta)
   end subroutine match

   subroutine sync(s)
      type(elec_t), intent(inout) :: s
      call cau_tab_sync(s)
   end subroutine sync

   function method_name() result(s)
      character(len=32) :: s
      s = 'tdhf_hub'
   end function method_name

   subroutine init()
      if (trim(electronic%method) == 'tdhf_hub') then
         call cau_tab_tdhf_setup(cau_pars_norb, cau_pars_bw, cau_pars_telec, &
                                 cau_pars_start, cau_pars_msbstp, cau_pars_nsbstp, &
                                 cau_pars_pmode, cau_pars_pamp)
      end if
   end subroutine init

   !------------------------------------------------------------------
   ! the record columns (the 11th reserved name, 2026-10-14)
   !------------------------------------------------------------------
   subroutine observables()
      if (trim(electronic%method) == 'tdhf_hub') then
         call rec_reg_col('n_imp_up', cau_tab_col_nup)
         call rec_reg_col('n_imp_dn', cau_tab_col_ndn)
         call rec_reg_col('mag_mom', cau_tab_col_mag)
         call rec_reg_col('e_elec_ev', cau_tab_col_eelec)
         call rec_reg_col('e_mf_ev', cau_tab_col_emf)
         call rec_reg_col('e_imp_ev', cau_tab_col_eimp)
         call rec_reg_col('e_bath_ev', cau_tab_col_ebath)
         call rec_reg_col('e_hyb_ev', cau_tab_col_ehyb)
      end if
   end subroutine observables
end module cau_tab_wires
