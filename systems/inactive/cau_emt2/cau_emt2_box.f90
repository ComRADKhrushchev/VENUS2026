!=====================================================================
! cau_emt2_box.f90 - the container's parameter channel + the
!   trajectory-termination judgment (problem knowledge, container-side)
! Design:
!   Declared container keys of input_qct.txt (the reg file declares the
!   vocabulary and pulls the buffered rows into cau_emt2_params_t at the
!   assembly seam): the TDHF/bath model keys (START_U0 / N_ORB /
!   BAND_WIDTH / TELECTRO / MSBSTP / PERTURB_MODE / PERTURB_AMP), the
!   periodic cell (BOXLX / BOXLY / SKEW), the EMT parameter file name
!   (EMT_FILE), and the optional Z_TERM termination height. Absent keys
!   take the production defaults of the C@Au_TDHF/2D run (N_ORB = 20,
!   BAND_WIDTH = 7 eV, TELECTRO = 300 K, MSBSTP = 4, hexagonal 17.7 A
!   cell); PERTURB_MODE defaults 0 = off (the old production default 1
!   is a non-energy-conserving stochastic law - enable it explicitly
!   per run). Z_TERM absent = the MAX_STEPS-only run (termination slot
!   unbound).
!   cau_emt2_term: the adsorbate height above the highest slab atom
!   exceeds Z_TERM -> the trajectory ends (scatter judgment) and the
!   rho lifecycle reseeds for the next trajectory.
!=====================================================================
module cau_emt2_box
   use state, only: state_t
   implicit none
   private
   public :: cau_emt2_params_t, cau_emt2_keys, cau_emt2_load, cau_emt2_term, cau_emt2_term_on
   public :: cau_pars_start, cau_pars_norb, cau_pars_bw, cau_pars_telec, &
             cau_pars_msbstp, cau_pars_pmode, cau_pars_pamp, &
             cau_pars_emt_file, cau_pars_boxl

   ! the parameter receiver: field defaults ARE the container defaults (an
   ! absent input key leaves its field untouched at the pull; the found
   ! flag carries the optional-slot provision semantics)
   type :: cau_emt2_params_t
      character(len=2)  :: start = 'GS'                       ! START_U0
      integer           :: norb = 20                          ! N_ORB
      real(8)           :: bw = 7.0d0                         ! BAND_WIDTH [eV]
      real(8)           :: telec = 300.0d0                    ! TELECTRO [K]
      integer           :: msbstp = 4                         ! MSBSTP
      integer           :: pmode = 0                          ! PERTURB_MODE
      real(8)           :: pamp = 0.05d0                      ! PERTURB_AMP
      character(len=64) :: emt_file = 'emt_2spec_params.txt'  ! EMT_FILE
      real(8)           :: boxl(3) = [17.7d0, 17.7d0, 60.0d0] ! BOXLX/BOXLY/SKEW [A, A, deg]
      real(8)           :: z_term = 0.0d0                     ! Z_TERM [A]
      logical           :: z_term_found = .false.
   end type cau_emt2_params_t

   ! -- buffered seam values (read by the reg _init, passed to setup calls) --
   character(len=2)  :: cau_pars_start = 'GS'
   integer           :: cau_pars_norb = 20
   real(8)           :: cau_pars_bw = 7.0d0
   real(8)           :: cau_pars_telec = 300.0d0
   integer           :: cau_pars_msbstp = 4
   integer           :: cau_pars_pmode = 0
   real(8)           :: cau_pars_pamp = 0.05d0
   character(len=64) :: cau_pars_emt_file = 'emt_2spec_params.txt'
   real(8)           :: cau_pars_boxl(3) = [17.7d0, 17.7d0, 60.0d0]  ! lx, ly, skew [A, A, deg]

   ! -- termination state --
   real(8)   :: z_term = 0.0d0
   logical   :: term_on = .false.
   logical   :: loaded = .false.

contains

   !------------------------------------------------------------------
   ! cau_emt2_keys() - the container's declared-key vocabulary (the reg
   !                 file hands it to input_declare_keys before the parse)
   !------------------------------------------------------------------
   function cau_emt2_keys() result(keys)
      character(len=16) :: keys(12)
      keys = [character(len=16) :: 'START_U0', 'N_ORB', 'BAND_WIDTH', 'TELECTRO', &
              'MSBSTP', 'PERTURB_MODE', 'PERTURB_AMP', 'EMT_FILE', 'BOXLX', &
              'BOXLY', 'SKEW', 'Z_TERM']
   end function cau_emt2_keys

   !------------------------------------------------------------------
   ! cau_emt2_load(p) - take the pulled parameter set into the module
   !                 state + load the EMT parameter file + hand the cell
   !                 to the physics layer
   !------------------------------------------------------------------
   subroutine cau_emt2_load(p)
      use cau_emt2_pot, only: cau_emt2_load_params, cau_emt2_set_box
      type(cau_emt2_params_t), intent(in) :: p
      if (loaded) then
         call stop_box('cau_emt2_load', 'the container was already loaded (double seam call)')
      end if
      if (.not. (p%start == 'GS' .or. p%start == '1D')) then
         call stop_box('cau_emt2_load', 'START_U0 must be GS or 1D (got "'//trim(p%start)//'")')
      end if
      ! 1. buffered seam values + termination state
      cau_pars_start = p%start
      cau_pars_norb = p%norb
      cau_pars_bw = p%bw
      cau_pars_telec = p%telec
      cau_pars_msbstp = p%msbstp
      cau_pars_pmode = p%pmode
      cau_pars_pamp = p%pamp
      cau_pars_emt_file = p%emt_file
      cau_pars_boxl = p%boxl
      z_term = p%z_term
      term_on = p%z_term_found
      ! 2. physics layer: EMT parameters + the periodic cell
      call cau_emt2_load_params(cau_pars_emt_file)
      call cau_emt2_set_box(cau_pars_boxl(1), cau_pars_boxl(2), cau_pars_boxl(3))
      loaded = .true.
   end subroutine cau_emt2_load

   logical function cau_emt2_term_on()
      cau_emt2_term_on = term_on
   end function cau_emt2_term_on

   !------------------------------------------------------------------
   ! cau_emt2_term(sta) - scatter judgment: adsorbate height above the
   !                 highest slab atom > Z_TERM; a true judgment also
   !                 reseeds the rho lifecycle for the next trajectory
   !------------------------------------------------------------------
   logical function cau_emt2_term(sta)
      use cau_emt2_tdhf, only: cau_emt2_tdhf_end_traj
      type(state_t), intent(in) :: sta
      integer :: i, natoms
      real(8) :: z_top
      natoms = size(sta%q)/3
      z_top = sta%q(6)
      do i = 2, natoms
         z_top = max(z_top, sta%q(3*i))
      end do
      cau_emt2_term = sta%q(3) - z_top > z_term
      if (cau_emt2_term) call cau_emt2_tdhf_end_traj()
   end function cau_emt2_term

   ! stop_box(who, msg) - the named-abort channel (stderr + STOP 1)
   subroutine stop_box(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_emt2_box error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container parameter channel is not usable)'
      stop 1
   end subroutine stop_box
end module cau_emt2_box
