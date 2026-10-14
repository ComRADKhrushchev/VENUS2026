!=====================================================================
! cau_emt2_wires.f90 - the Tier 1 reserved-name wrappers (the folder's
!   auto-weave contract: prebuild weaves these names into registry_gen
!   when the folder sits reg-less at systems/ top level; each delegates
!   to the public API the check programs and probe hooks keep using).
!   Replaces reg_cau_emt2.f90 (the reg form retired 2026-10-14,
!   docs/plans/2026-10-14-reg-retirement-observables.md).
!   A dedicated wires file (not in cau_emt2_box) keeps the box free of
!   use input/recorder/config - the criterion recipes compile the
!   box standalone without those modules.
! Units: none
!=====================================================================
module cau_emt2_wires
   use state,          only: state_t, elec_t
   use config,         only: electronic
   use input,          only: input_declare_keys, pull
   use recorder,       only: rec_reg_col
   use cau_emt2_box,    only: cau_emt2_params_t, cau_emt2_keys, cau_emt2_load, &
                             cau_emt2_term, cau_emt2_term_on, &
                             cau_pars_start, cau_pars_norb, cau_pars_bw, cau_pars_telec, &
                             cau_pars_msbstp, cau_pars_pmode, cau_pars_pamp
   use cau_emt2_tdhf,   only: cau_emt2_pes, cau_emt2_prop, cau_emt2_match, cau_emt2_tdhf_setup, cau_emt2_sync, &
                             cau_emt2_col_nup, cau_emt2_col_ndn, cau_emt2_col_mag, &
                             cau_emt2_col_eelec, cau_emt2_col_emf
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
      call cau_emt2_pes(q, mass, g, v)
   end subroutine pes

   function keys() result(kk)
      character(len=16) :: kk(12)
      kk = cau_emt2_keys()
   end function keys

   subroutine load()
      type(cau_emt2_params_t) :: p
      call pull('START_U0',     p%start)
      call pull('N_ORB',        p%norb)
      call pull('BAND_WIDTH',   p%bw)
      call pull('TELECTRO',     p%telec)
      call pull('MSBSTP',       p%msbstp)
      call pull('PERTURB_MODE', p%pmode)
      call pull('PERTURB_AMP',  p%pamp)
      call pull('EMT_FILE',     p%emt_file)
      call pull('BOXLX',        p%boxl(1))
      call pull('BOXLY',        p%boxl(2))
      call pull('SKEW',         p%boxl(3))
      call pull('Z_TERM',       p%z_term, p%z_term_found)
      call cau_emt2_load(p)
   end subroutine load

   logical function term(sta)
      type(state_t), intent(in) :: sta
      if (.not. cau_emt2_term_on()) then   ! the reg-era bind gate, carried
         term = .false.                   ! by the wrapper (the body does
         return                           ! not self-gate; an unbound run
      end if                              ! must answer .false.)
      term = cau_emt2_term(sta)
   end function term

   !------------------------------------------------------------------
   ! the method package (tdhf_hub)
   !------------------------------------------------------------------
   subroutine prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
      call cau_emt2_prop(sta, dt)
   end subroutine prop

   subroutine match(sta)
      type(state_t), intent(inout) :: sta
      call cau_emt2_match(sta)
   end subroutine match

   subroutine sync(s)
      type(elec_t), intent(inout) :: s
      call cau_emt2_sync(s)
   end subroutine sync

   function method_name() result(s)
      character(len=32) :: s
      s = 'tdhf_hub'
   end function method_name

   subroutine init()
      if (trim(electronic%method) == 'tdhf_hub') then
         call cau_emt2_tdhf_setup(cau_pars_norb, cau_pars_bw, cau_pars_telec, &
                                  cau_pars_start, cau_pars_msbstp, &
                                  cau_pars_pmode, cau_pars_pamp)
      end if
   end subroutine init

   !------------------------------------------------------------------
   ! the record columns (the 11th reserved name, 2026-10-14)
   !------------------------------------------------------------------
   subroutine observables()
      if (trim(electronic%method) == 'tdhf_hub') then
         call rec_reg_col('n_imp_up', cau_emt2_col_nup)
         call rec_reg_col('n_imp_dn', cau_emt2_col_ndn)
         call rec_reg_col('mag_mom', cau_emt2_col_mag)
         call rec_reg_col('e_elec_ev', cau_emt2_col_eelec)
         call rec_reg_col('e_mf_ev', cau_emt2_col_emf)
      end if
   end subroutine observables
end module cau_emt2_wires
