!=====================================================================
! cau_rst_box.f90 - the container's parameter channel + the
!   trajectory-termination judgment (problem knowledge, container-side)
! Design:
!   Declared container keys of input_qct.txt (the reg file declares the
!   vocabulary and pulls the buffered rows into cau_rst_params_t at the
!   assembly seam). Engine files: SLAB_FILE / WEIGHTS_FILE / BVK_FILE
!   name the Slab.xyz + nn_weights_rst_hubbard.txt + Analytic_Potential.txt
!   copies the run directory carries (the container folder stores
!   reference copies; an absent BVK_FILE key falls to the engine
!   defaults). The TDHF/bath keys and the optional Z_TERM termination
!   height mirror cau_emt2_box; PERTURB_MODE defaults 0 = off.
!   cau_rst_term: the adsorbate height above the highest slab atom
!   exceeds Z_TERM -> the trajectory ends (scatter judgment) and the
!   rho lifecycle reseeds for the next trajectory.
!=====================================================================
module cau_rst_box
   use state, only: state_t
   implicit none
   private
   public :: cau_rst_params_t, cau_rst_keys, cau_rst_load, cau_rst_term, cau_rst_term_on
   public :: cau_pars_start, cau_pars_norb, cau_pars_bw, cau_pars_telec, &
             cau_pars_msbstp, cau_pars_pmode, cau_pars_pamp

   ! the parameter receiver: field defaults ARE the container defaults (an
   ! absent input key leaves its field untouched at the pull; the found
   ! flag carries the optional-slot provision semantics)
   type :: cau_rst_params_t
      character(len=2)  :: start = 'GS'                            ! START_U0
      integer           :: norb = 20                               ! N_ORB
      real(8)           :: bw = 7.0d0                              ! BAND_WIDTH [eV]
      real(8)           :: telec = 300.0d0                         ! TELECTRO [K]
      integer           :: msbstp = 4                              ! MSBSTP
      integer           :: pmode = 0                               ! PERTURB_MODE
      real(8)           :: pamp = 0.05d0                           ! PERTURB_AMP
      character(len=64) :: slab_file = 'Slab.xyz'                  ! SLAB_FILE
      character(len=64) :: weights_file = 'nn_weights_rst_hubbard.txt'  ! WEIGHTS_FILE
      character(len=64) :: bvk_file = 'Analytic_Potential.txt'     ! BVK_FILE
      real(8)           :: z_term = 0.0d0                          ! Z_TERM [A]
      logical           :: z_term_found = .false.
   end type cau_rst_params_t

   ! -- buffered seam values (read by the reg _init, passed to setup calls) --
   character(len=2)  :: cau_pars_start = 'GS'
   integer           :: cau_pars_norb = 20
   real(8)           :: cau_pars_bw = 7.0d0
   real(8)           :: cau_pars_telec = 300.0d0
   integer           :: cau_pars_msbstp = 4
   integer           :: cau_pars_pmode = 0
   real(8)           :: cau_pars_pamp = 0.05d0

   ! -- termination state --
   real(8)   :: z_term = 0.0d0
   logical   :: term_on = .false.
   logical   :: loaded = .false.

contains

   !------------------------------------------------------------------
   ! cau_rst_keys() - the container's declared-key vocabulary (the reg
   !                 file hands it to input_declare_keys before the parse)
   !------------------------------------------------------------------
   function cau_rst_keys() result(keys)
      character(len=16) :: keys(11)
      keys = [character(len=16) :: 'START_U0', 'N_ORB', 'BAND_WIDTH', 'TELECTRO', &
              'MSBSTP', 'PERTURB_MODE', 'PERTURB_AMP', 'SLAB_FILE', 'WEIGHTS_FILE', &
              'BVK_FILE', 'Z_TERM']
   end function cau_rst_keys

   !------------------------------------------------------------------
   ! cau_rst_load(p) - take the pulled parameter set into the module
   !                 state + run the engine init (slab + weights + BVK
   !                 + tensors)
   !------------------------------------------------------------------
   subroutine cau_rst_load(p)
      use cau_rst_pot, only: cau_rst_load_engine => cau_rst_load
      type(cau_rst_params_t), intent(in) :: p
      if (loaded) then
         call stop_box('cau_rst_load', 'the container was already loaded (double seam call)')
      end if
      if (.not. (p%start == 'GS' .or. p%start == '1D')) then
         call stop_box('cau_rst_load', 'START_U0 must be GS or 1D (got "'//trim(p%start)//'")')
      end if
      ! 1. buffered seam values + termination state
      cau_pars_start = p%start
      cau_pars_norb = p%norb
      cau_pars_bw = p%bw
      cau_pars_telec = p%telec
      cau_pars_msbstp = p%msbstp
      cau_pars_pmode = p%pmode
      cau_pars_pamp = p%pamp
      z_term = p%z_term
      term_on = p%z_term_found
      ! 2. engine one-time initialization (slab + weights + BVK + tensors)
      call cau_rst_load_engine(p%slab_file, p%weights_file, p%bvk_file)
      loaded = .true.
   end subroutine cau_rst_load

   logical function cau_rst_term_on()
      cau_rst_term_on = term_on
   end function cau_rst_term_on

   !------------------------------------------------------------------
   ! cau_rst_term(sta) - scatter judgment (height above the slab)
   !------------------------------------------------------------------
   logical function cau_rst_term(sta)
      use cau_rst_tdhf, only: cau_rst_tdhf_end_traj
      type(state_t), intent(in) :: sta
      integer :: i, natoms
      real(8) :: z_top
      natoms = size(sta%q)/3
      z_top = sta%q(6)
      do i = 2, natoms
         z_top = max(z_top, sta%q(3*i))
      end do
      cau_rst_term = sta%q(3) - z_top > z_term
      if (cau_rst_term) call cau_rst_tdhf_end_traj()
   end function cau_rst_term

   subroutine stop_box(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_rst_box error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container parameter channel is not usable)'
      stop 1
   end subroutine stop_box
end module cau_rst_box
