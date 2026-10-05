!=====================================================================
! cau_peri_box.f90 - the container's parameter channel + the
!   trajectory-termination judgment (problem knowledge, container-side)
! Design:
!   Declared container keys of input_qct.txt (the reg file declares the
!   vocabulary and pulls the buffered rows into cau_peri_params_t at the
!   assembly seam): the TDHF/bath model keys (START_U0 / N_ORB /
!   BAND_WIDTH / TELECTRO / MSBSTP / PERTURB_MODE / PERTURB_AMP), the
!   S-Y oscillator constants (K_SURF / K_GHOST, 3-vectors [eV/A^2]), the
!   ghost equilibrium (Q0_GHOST [A]), the NN path (USE_NN 0/1 + NN_FILE,
!   run-directory file when 1; default 0 = analytic READAP defaults) and
!   the optional Z_TERM termination height (measures the incident atom
!   height above the anchor; absent = the MAX_STEPS-only run, slot
!   unbound). Absent keys take the buffered defaults of the analytic
!   production run; PERTURB_MODE defaults 0 = off (the old production
!   default 1 is a non-energy-conserving stochastic law - enable it
!   explicitly per run).
!   ATOM ORDER: the roster follows the family name-sort convention
!   (C.xyz before Slab.poscar) -> [atom 1 = incident C, atom 2 = anchor S,
!   atom 3 = ghost Y]. The ENGINE's internal layout is [S, Y, C] (the
!   periodic kernels index the host at q(1:3)), so the boundary wrappers
!   cau_peri_pes_r / cau_peri_prop_r permute roster -> engine (and the
!   gradient back), and the termination judgment reads roster indices
!   (the incident C above the anchor: q(3) - q(6)).
!   SURFACE PARADIGM (2026-09-28): the Slab rides a POSCAR (hexagonal
!   cell a = 2.95 A = the engine period) so the run derives
!   surface_model = 1 - PROJECTILE = C (no TARGET), A_LAT = 2.95, the
!   substrate fragment sampled by its DIST_SCHEME word; the incident
!   member is incident_surface (drop site: N_AIM 0/1/2, the fixed site
!   AIM_X/AIM_Y [A] in the POSCAR frame).
!=====================================================================
module cau_peri_box
   use state, only: state_t
   implicit none
   private
   public :: cau_peri_params_t, cau_peri_keys, cau_peri_load, cau_peri_term, cau_peri_term_on
   public :: cau_peri_pes_r, cau_peri_prop_r
   public :: cau_pars_start, cau_pars_norb, cau_pars_bw, cau_pars_telec, &
             cau_pars_msbstp, cau_pars_pmode, cau_pars_pamp

   ! the parameter receiver: field defaults ARE the container defaults (an
   ! absent input key leaves its field untouched at the pull; the found
   ! flag carries the optional-slot provision semantics)
   type :: cau_peri_params_t
      character(len=2)  :: start = 'GS'                       ! START_U0
      integer           :: norb = 20                          ! N_ORB
      real(8)           :: bw = 7.0d0                         ! BAND_WIDTH [eV]
      real(8)           :: telec = 300.0d0                    ! TELECTRO [K]
      integer           :: msbstp = 4                         ! MSBSTP
      integer           :: pmode = 0                          ! PERTURB_MODE
      real(8)           :: pamp = 0.05d0                      ! PERTURB_AMP
      real(8)           :: ksurf(3) = [1.0d0, 1.0d0, 1.0d0]   ! K_SURF [eV/A^2]
      real(8)           :: kghost(3) = [1.0d0, 1.0d0, 1.0d0]  ! K_GHOST [eV/A^2]
      real(8)           :: q0ghost(3) = [0.0d0, 0.0d0, 0.0d0] ! Q0_GHOST [A]
      integer           :: use_nn = 0                         ! USE_NN
      character(len=64) :: nn_file = 'nn_weights.txt'         ! NN_FILE
      real(8)           :: z_term = 0.0d0                     ! Z_TERM [A]
      logical           :: z_term_found = .false.
   end type cau_peri_params_t

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
   real(8)   :: ksurf_d(3) = [1.0d0, 1.0d0, 1.0d0]
   real(8)   :: kghost_d(3) = [1.0d0, 1.0d0, 1.0d0]
   real(8)   :: q0g_d(3) = [0.0d0, 0.0d0, 0.0d0]
   logical   :: term_on = .false.
   logical   :: loaded = .false.

contains

   !------------------------------------------------------------------
   ! cau_peri_keys() - the container's declared-key vocabulary (the reg
   !                 file hands it to input_declare_keys before the parse)
   !------------------------------------------------------------------
   function cau_peri_keys() result(keys)
      character(len=16) :: keys(13)
      keys = [character(len=16) :: 'START_U0', 'N_ORB', 'BAND_WIDTH', 'TELECTRO', &
              'MSBSTP', 'PERTURB_MODE', 'PERTURB_AMP', 'K_SURF', 'K_GHOST', &
              'Q0_GHOST', 'USE_NN', 'NN_FILE', 'Z_TERM']
   end function cau_peri_keys

   !------------------------------------------------------------------
   ! cau_peri_load(p) - take the pulled parameter set into the module
   !                 state + run the engine init (oscillator constants,
   !                 then the one-time engine setup, NN when selected)
   !------------------------------------------------------------------
   subroutine cau_peri_load(p)
      use cau_peri_pot, only: cau_peri_load_engine => cau_peri_load
      use surface_potentials, only: k_surf, k_ghost, Q0_ghost
      type(cau_peri_params_t), intent(in) :: p
      if (loaded) then
         call stop_box('cau_peri_load', 'the container was already loaded (double seam call)')
      end if
      if (.not. (p%start == 'GS' .or. p%start == '1D')) then
         call stop_box('cau_peri_load', 'START_U0 must be GS or 1D (got "'//trim(p%start)//'")')
      end if
      ! 1. buffered seam values + termination state
      cau_pars_start = p%start
      cau_pars_norb = p%norb
      cau_pars_bw = p%bw
      cau_pars_telec = p%telec
      cau_pars_msbstp = p%msbstp
      cau_pars_pmode = p%pmode
      cau_pars_pamp = p%pamp
      ksurf_d = p%ksurf
      kghost_d = p%kghost
      q0g_d = p%q0ghost
      z_term = p%z_term
      term_on = p%z_term_found
      ! 2. oscillator constants + engine one-time initialization
      k_surf = ksurf_d
      k_ghost = kghost_d
      Q0_ghost = q0g_d
      if (p%use_nn == 1) then
         call cau_peri_load_engine(p%nn_file)
      else
         call cau_peri_load_engine()
      end if
      loaded = .true.
   end subroutine cau_peri_load

   logical function cau_peri_term_on()
      cau_peri_term_on = term_on
   end function cau_peri_term_on

   !------------------------------------------------------------------
   ! cau_peri_term(sta) - scatter judgment (roster indices: the incident
   !                 C, atom 1, above the anchor S, atom 2)
   !------------------------------------------------------------------
   logical function cau_peri_term(sta)
      use cau_peri_tdhf, only: cau_peri_tdhf_end_traj
      type(state_t), intent(in) :: sta
      cau_peri_term = sta%q(3) - sta%q(6) > z_term
      if (cau_peri_term) call cau_peri_tdhf_end_traj()
   end function cau_peri_term

   !------------------------------------------------------------------
   ! cau_peri_pes_r(q, mass, g, v) - the roster-order force entry: permute
   !                 [C, S, Y] -> the engine layout [S, Y, C], call the
   !                 composed entry, permute the gradient back. The bound
   !                 PES contract is fully assigned in roster order
   !------------------------------------------------------------------
   subroutine cau_peri_pes_r(q, mass, g, v)
      use cau_peri_tdhf, only: cau_peri_pes
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      real(8) :: qe(9), ge(9), me(3)
      if (size(q) /= 9 .or. size(mass) /= 3) then
         call stop_box('cau_peri_pes_r', 'the container needs exactly 3 atoms (9 dofs)'// &
                       ' - the system folder must define C then the two-atom slab'// &
                       ' (anchor S, ghost Y) in this order')
      end if
      qe(1:3) = q(4:6)      ! engine atom 1 = roster atom 2 (the anchor S)
      qe(4:6) = q(7:9)      ! engine atom 2 = roster atom 3 (the ghost Y)
      qe(7:9) = q(1:3)      ! engine atom 3 = roster atom 1 (the incident C)
      me(1) = mass(2)
      me(2) = mass(3)
      me(3) = mass(1)
      call cau_peri_pes(qe, me, ge, v)
      g(1:3) = ge(7:9)      ! scatter back to roster order
      g(4:6) = ge(1:3)
      g(7:9) = ge(4:6)
   end subroutine cau_peri_pes_r

   !------------------------------------------------------------------
   ! cau_peri_prop_r(sta, dt) - the roster-order propagation entry: the
   !                 same permutation on a local geometry copy (the rho
   !                 lifecycle lives in the engine's module state; sta
   !                 itself is handed through untouched)
   !------------------------------------------------------------------
   subroutine cau_peri_prop_r(sta, dt)
      use cau_peri_tdhf, only: cau_peri_prop
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
      type(state_t) :: loc
      if (size(sta%q) /= 9) then
         call stop_box('cau_peri_prop_r', 'the container needs exactly 3 atoms (9 dofs)'// &
                       ' - the system folder must define C then the two-atom slab'// &
                       ' (anchor S, ghost Y) in this order')
      end if
      loc = sta
      loc%q(1:3) = sta%q(4:6)
      loc%q(4:6) = sta%q(7:9)
      loc%q(7:9) = sta%q(1:3)
      call cau_peri_prop(loc, dt)
   end subroutine cau_peri_prop_r

   subroutine stop_box(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_peri_box error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container parameter channel is not usable)'
      stop 1
   end subroutine stop_box
end module cau_peri_box
