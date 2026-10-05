!=====================================================================
! beam_laws.f90 - shared beam-law component for the incident-channel
!   members (orientation, azimuth, aiming, energy draw layers + the
!   COM placement/velocity assignment helpers)
! Design:
!   A level-2 component (no reg file - it registers nothing): the three
!   incident members (inc_surface / inc_pair) each host one beam_t
!   instance, so the stochastic law lives in ONE copy (registry-family
!   members never depend on each other; a shared component is the
!   sanctioned shape). beam_init validates + derives (mu_eff is the
!   CALLER's effective mass - projectile mass for the surface member,
!   reduced mass for the pair); beam_draw consumes the RNG into the
!   instance's record (draw order: orientation, aiming, drop point,
!   speed - the historical order); the realize-side helpers assign COM
!   velocities (an assignment, never an accumulation - see beam_set_com_vel).
!   Drop-point modes: 0 = origin, 1 = uniform in the unit cell (needs the
!   buffered cell), 2 = the fixed site (aim_x, aim_y) [Angstrom, the
!   POSCAR frame] - the surface member's fixed-incidence option.
! Units: angles [rad]; t_trans [K]; e_rel [kcal/mol] via e_conv; speeds [Ang/(10 fs)]; p [amu*Ang/(10 fs)].
!=====================================================================
module beam_laws
   use consts,       only: pi, two_pi, e_conv, r_kcal
   use rng,          only: rng_u, rng_gamma
   use state,        only: state_t
   use config_atoms, only: list_atoms
   implicit none
   private
   public :: beam_params_t, beam_t, beam_init, beam_draw
   public :: beam_com_shift, beam_set_com_vel

   integer, parameter :: max_reject = 100000  ! rejection-sampling draw cap [-]

   ! Beam shape parameters (the collision energy / impact parameter / separation
   ! are system-level prescriptions in config%reactants - not duplicated here)
   type :: beam_params_t
      integer :: n_e_rel  = 1      ! collision-energy prescription mode [-] (0 = thermal
                                   ! population, 1 = fixed value, 2 = density-weighted)
      real(8) :: t_trans  = 0.0d0  ! translational temperature [K] (mode 0: E_rel =
                                   ! r_kcal*gamma(2)*t_trans)
      real(8) :: v_width  = 0.0d0  ! velocity width of the mode-2 density [Angstrom/(10 fs)]
      integer :: n_thta   = 0      ! polar-angle mode [-] (0 = fixed beam theta = thta_max,
                                   ! 1 = P(theta) prop. to sin(theta) sampling)
      real(8) :: thta_max = 0.0d0  ! polar-angle bound [rad] (fixed value / sampling bound)
      integer :: n_chi    = 0      ! azimuth mode [-] (0 = fixed chi, 1 = uniform [0, 2*pi))
      real(8) :: chi      = 0.0d0  ! azimuth [rad] (fixed value)
      integer :: n_b      = 0      ! aiming mode [-] (0 = fixed b = b_max, 1 = b^2 uniform
                                   ! on [0, b_max^2] - the collision disk)
      integer :: n_aim    = 0      ! drop-point mode [-] (0 = origin, 1 = uniform in the
                                   ! unit cell, 2 = the fixed site aim_x/aim_y)
      real(8) :: aim_x    = 0.0d0  ! fixed drop-point x [Angstrom] (mode 2)
      real(8) :: aim_y    = 0.0d0  ! fixed drop-point y [Angstrom] (mode 2)
      logical :: aim_x_given = .false. ! AIM_X presence flag (mode-2 validation)
      logical :: aim_y_given = .false. ! AIM_Y presence flag (mode-2 validation)
   end type

   ! One beam instance: derived law state + the drawn prescription record
   type :: beam_t
      logical :: inited = .false.   ! init completed [flag]
      integer :: mode_e = 1         ! collision-energy mode copy [-]
      integer :: mode_th = 0        ! polar-angle mode copy [-]
      integer :: mode_ch = 0        ! azimuth mode copy [-]
      integer :: mode_b = 0         ! aiming mode copy [-]
      integer :: mode_aim = 0       ! drop-point mode copy [-]
      real(8) :: cos_tmax = 1.0d0   ! cos(thta_max) (fixed polar value / sampling bound) [-]
      real(8) :: chi_fix = 0.0d0    ! fixed azimuth [rad]
      real(8) :: t_trans_c = 0.0d0  ! translational temperature copy [K]
      real(8) :: e_rel_i = 0.0d0    ! collision energy [internal] (modes 1 and 2)
      real(8) :: b_max_c = 0.0d0    ! collision-disk radius [Angstrom]
      real(8) :: r_sep_c = 0.0d0    ! drop height / beam back-off scale [Angstrom]
      real(8) :: mu_eff = 0.0d0     ! effective mass of the relative motion [amu]
      real(8) :: v0_c = 0.0d0       ! mode-2 center speed [Angstrom/(10 fs)]
      real(8) :: alpha_c = 0.0d0    ! mode-2 velocity width [Angstrom/(10 fs)]
      real(8) :: v_star = 0.0d0     ! mode-2 density peak [Angstrom/(10 fs)]
      real(8) :: v_hi = 0.0d0       ! mode-2 proposal upper bound [Angstrom/(10 fs)]
      real(8) :: f_star = 0.0d0     ! mode-2 envelope height DENS(v*) [-]
      real(8) :: a1_vec(2) = 0.0d0  ! unit-cell edge a1 [Angstrom]
      real(8) :: a2_vec(2) = 0.0d0  ! unit-cell edge a2 [Angstrom]
      real(8) :: aim_x_c = 0.0d0    ! fixed drop-point x copy [Angstrom] (mode 2)
      real(8) :: aim_y_c = 0.0d0    ! fixed drop-point y copy [Angstrom] (mode 2)
      ! the drawn prescription record (the draw-phase product realize consumes)
      real(8) :: costh_s = 0.0d0    ! drawn cos(theta) [-]
      real(8) :: sinth_s = 0.0d0    ! drawn sin(theta) [-]
      real(8) :: chi_s = 0.0d0      ! drawn azimuth chi [rad]
      real(8) :: bval_s = 0.0d0     ! drawn aiming radius [Angstrom]
      real(8) :: phi_s = 0.0d0      ! drawn disk azimuth [rad]
      real(8) :: rx0_s = 0.0d0      ! drawn drop point x [Angstrom]
      real(8) :: ry0_s = 0.0d0      ! drawn drop point y [Angstrom]
      real(8) :: vv_s = 0.0d0       ! drawn collision speed [Angstrom/(10 fs)]
   end type
contains
   !------------------------------------------------------------------
   ! beam_init(b, p, mu_eff, e_rel, b_max, r_sep, caller [, onsite_ok]) -
   !                 validate the shape parameters, adopt the system-level
   !                 prescriptions, precompute the mode-2 envelope and the
   !                 cell vectors
   !------------------------------------------------------------------
   subroutine beam_init(b, p, mu_eff, e_rel, b_max, r_sep, caller, onsite_ok)
      type(beam_t), intent(inout) :: b
      type(beam_params_t), intent(in) :: p
      real(8), intent(in) :: mu_eff          ! effective mass [amu] (caller's paradigm)
      real(8), intent(in) :: e_rel           ! collision energy [kcal/mol]
      real(8), intent(in) :: b_max           ! collision-disk radius [Angstrom]
      real(8), intent(in) :: r_sep           ! initial separation [Angstrom]
      character(len=*), intent(in) :: caller ! abort prefix (the member's init name)
      logical, intent(in), optional :: onsite_ok ! B_MAX = 0 legality (surface: on-site beam)
      real(8) :: dum
      logical :: b_zero_ok

      b_zero_ok = .false.
      if (present(onsite_ok)) b_zero_ok = onsite_ok

      ! 1. parameter validation (fail-loud, naming the field)
      if (p%n_e_rel < 0 .or. p%n_e_rel > 2) then
         call stop_beami(caller, 'collision-energy mode n_e_rel must be 0 (thermal), '// &
            '1 (fixed) or 2 (density-weighted), got', p%n_e_rel)
      end if
      if (p%n_e_rel == 0 .and. p%t_trans <= 0.0d0) then
         dum = p%t_trans
         call stop_beamz(caller, 'thermal mode 0 requires a translational temperature '// &
            't_trans [K] > 0, got', dum)
      end if
      if (p%n_e_rel == 2 .and. p%v_width <= 0.0d0) then
         dum = p%v_width
         call stop_beamz(caller, 'density-weighted mode 2 requires a velocity width '// &
            'v_width > 0, got', dum)
      end if
      if (p%n_thta < 0 .or. p%n_thta > 1) then
         call stop_beami(caller, 'polar-angle mode n_thta must be 0 (fixed) or 1 (sampled), got', p%n_thta)
      end if
      if (p%n_thta == 1) then
         if (p%thta_max <= 0.0d0 .or. p%thta_max > 0.5d0*pi) then
            dum = p%thta_max
            call stop_beamz(caller, 'sampled polar bound thta_max [rad] must lie in '// &
               '(0, pi/2] (a degenerate or sub-surface range), got', dum)
         end if
      else
         if (p%thta_max < 0.0d0 .or. p%thta_max > 0.5d0*pi) then
            dum = p%thta_max
            call stop_beamz(caller, 'fixed polar angle thta_max [rad] must lie in '// &
               '[0, pi/2], got', dum)
         end if
      end if
      if (p%n_chi < 0 .or. p%n_chi > 1) then
         call stop_beami(caller, 'azimuth mode n_chi must be 0 (fixed) or 1 (sampled), got', p%n_chi)
      end if
      if (p%n_b < 0 .or. p%n_b > 1) then
         call stop_beami(caller, 'aiming mode n_b must be 0 (fixed b_max) or 1 (collision disk), got', p%n_b)
      end if
      if (p%n_aim < 0 .or. p%n_aim > 2) then
         call stop_beami(caller, 'drop-point mode n_aim must be 0 (origin), 1 (unit cell) '// &
            'or 2 (fixed site), got', p%n_aim)
      end if
      if (p%n_aim == 1) then
         if (list_atoms%box_lx <= 0.0d0 .or. list_atoms%box_ly <= 0.0d0 .or. sin(list_atoms%skew) <= 0.0d0) then
            call stop_beamx(caller, 'drop-point mode 1 requires a unit cell in the atom list '// &
               '(box_lx, box_ly, skew - none buffered)')
         end if
      end if
      if (p%n_aim == 2) then
         if (.not. (p%aim_x_given .and. p%aim_y_given)) then
            call stop_beamx(caller, 'drop-point mode 2 requires both AIM_X and AIM_Y '// &
               '(the fixed drop site [Angstrom, POSCAR frame])')
         end if
      end if

      ! 2. system-level incident prescriptions (validated here, adopted once)
      if ((p%n_e_rel == 1 .or. p%n_e_rel == 2) .and. e_rel <= 0.0d0) then
         dum = e_rel
         call stop_beamz(caller, 'collision-energy modes 1 and 2 require a positive collision '// &
            'energy e_rel [kcal/mol] in the reactants provision, got', dum)
      end if
      if (b_max < 0.0d0 .or. (b_max == 0.0d0 .and. .not. b_zero_ok)) then
         dum = b_max
         call stop_beamz(caller, 'the collision-disk radius b_max [Angstrom] must be positive, got', dum)
      end if
      if (b_max == 0.0d0 .and. p%n_b /= 0) then
         call stop_beamx(caller, 'an on-site beam (B_MAX = 0) requires the fixed-aiming '// &
            'mode N_B = 0 (the collision-disk draw would be degenerate)')
      end if
      if (r_sep <= 0.0d0) then
         dum = r_sep
         call stop_beamz(caller, 'the initial separation r_sep [Angstrom] must be positive, got', dum)
      end if

      ! 3. precompute (mode copies, trigonometry, mode-2 envelope, cell vectors)
      b%mode_e = p%n_e_rel
      b%mode_th = p%n_thta
      b%mode_ch = p%n_chi
      b%mode_b = p%n_b
      b%mode_aim = p%n_aim
      b%cos_tmax = cos(p%thta_max)
      b%chi_fix = p%chi
      b%t_trans_c = p%t_trans
      b%e_rel_i = e_conv*e_rel
      b%b_max_c = b_max
      b%r_sep_c = r_sep
      b%mu_eff = mu_eff
      b%aim_x_c = p%aim_x
      b%aim_y_c = p%aim_y
      if (b%mode_e == 2) then
         b%alpha_c = p%v_width
         b%v0_c = sqrt(2.0d0*b%e_rel_i/mu_eff)
         b%v_star = 0.5d0*b%v0_c + sqrt(0.25d0*b%v0_c*b%v0_c + 1.5d0*b%alpha_c*b%alpha_c)
         b%v_hi = 5.0d0*b%v_star
         b%f_star = dens(b, b%v_star)
      end if
      b%a1_vec = (/ list_atoms%box_lx, 0.0d0 /)
      b%a2_vec = list_atoms%box_ly*(/ cos(list_atoms%skew), sin(list_atoms%skew) /)
      b%inited = .true.
   end subroutine beam_init

   !------------------------------------------------------------------
   ! beam_draw(b) - draw the whole collision prescription (the historical
   !                 order): beam orientation, aiming radius + disk azimuth,
   !                 drop point, collision speed; no state access
   !------------------------------------------------------------------
   subroutine beam_draw(b)
      type(beam_t), intent(inout) :: b
      real(8) :: s1, s2
      integer :: n_try
      if (.not. b%inited) then
         call stop_beamx('beam_draw', 'draw called before beam_init - the instance carries '// &
            'no parameter set')
      end if

      ! 1. beam orientation: cos(theta) (fixed or P(theta) prop. sin(theta)),
      !    azimuth chi (fixed or uniform)
      if (b%mode_th == 1) then
         b%costh_s = 1.0d0 - rng_u()*(1.0d0 - b%cos_tmax)
      else
         b%costh_s = b%cos_tmax
      end if
      b%sinth_s = sqrt(max(0.0d0, 1.0d0 - b%costh_s*b%costh_s))
      if (b%mode_ch == 1) then
         b%chi_s = two_pi*rng_u()
      else
         b%chi_s = b%chi_fix
      end if

      ! 2. aiming radius (fixed b_max or b^2 uniform) and disk azimuth phi
      !    (uniform, or chi in the fixed-b mode)
      if (b%mode_b == 1) then
         b%bval_s = b%b_max_c*sqrt(rng_u())
         b%phi_s = two_pi*rng_u()
      else
         b%bval_s = b%b_max_c
         b%phi_s = b%chi_s
      end if

      ! 3. drop point (origin / unit cell / fixed site)
      b%rx0_s = 0.0d0
      b%ry0_s = 0.0d0
      select case (b%mode_aim)
      case (1)
         s1 = rng_u()
         s2 = rng_u()
         b%rx0_s = s1*b%a1_vec(1) + s2*b%a2_vec(1)
         b%ry0_s = s1*b%a1_vec(2) + s2*b%a2_vec(2)
      case (2)
         b%rx0_s = b%aim_x_c
         b%ry0_s = b%aim_y_c
      end select

      ! 4. collision speed per energy mode
      select case (b%mode_e)
      case (0)
         b%vv_s = sqrt(2.0d0*e_conv*r_kcal*b%t_trans_c*rng_gamma(2)/b%mu_eff)
      case (1)
         b%vv_s = sqrt(2.0d0*b%e_rel_i/b%mu_eff)
      case (2)
         do n_try = 1, max_reject
            b%vv_s = b%v_hi*rng_u()
            if (rng_u()*b%f_star <= dens(b, b%vv_s)) exit
            if (n_try == max_reject) then
               call stop_beamx('beam_draw', 'rejection sampling of the collision speed '// &
                  'exceeded its draw cap (the density ratio is inconsistent with the envelope)')
            end if
         end do
      end select
   end subroutine beam_draw

   !------------------------------------------------------------------
   ! beam_com_shift(sta, list, com_new) - translate the fragment's atoms so
   !                 their mass center lands on com_new (internal geometry
   !                 untouched)
   !------------------------------------------------------------------
   subroutine beam_com_shift(sta, list, com_new)
      type(state_t), intent(inout) :: sta
      integer, intent(in) :: list(:)
      real(8), intent(in) :: com_new(3)
      real(8) :: com(3), shift(3)
      integer :: ia
      com = 0.0d0
      do ia = 1, size(list)
         com = com + list_atoms%mass(list(ia))*sta%q(3*list(ia) - 2:3*list(ia))
      end do
      com = com/sum(list_atoms%mass(list))
      shift = com_new - com
      do ia = 1, size(list)
         sta%q(3*list(ia) - 2:3*list(ia)) = sta%q(3*list(ia) - 2:3*list(ia)) + shift
      end do
   end subroutine beam_com_shift

   !------------------------------------------------------------------
   ! beam_set_com_vel(sta, list, vel) - place the fragment's COM velocity at
   !                 vel. An assignment, not an accumulation - the
   !                 per-fragment members write zero-net-momentum interiors
   !                 (a plain add suffices for them and is bitwise
   !                 preserved), but a monoatomic projectile has no interior
   !                 writer: its carry from the previous trajectory must be
   !                 wiped or the beam momentum accumulates once per
   !                 trajectory
   !------------------------------------------------------------------
   subroutine beam_set_com_vel(sta, list, vel)
      type(state_t), intent(inout) :: sta
      integer, intent(in) :: list(:)
      real(8), intent(in) :: vel(3)
      real(8) :: wt, v_com(3)
      integer :: ia, j
      wt = sum(list_atoms%mass(list))
      v_com = 0.0d0
      do ia = 1, size(list)
         v_com = v_com + sta%p(3*list(ia) - 2:3*list(ia))
      end do
      v_com = v_com/wt
      do ia = 1, size(list)
         do j = 1, 3
            sta%p(3*list(ia) - 3 + j) = sta%p(3*list(ia) - 3 + j) &
               + list_atoms%mass(list(ia))*(vel(j) - v_com(j))
         end do
      end do
   end subroutine beam_set_com_vel

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! dens(b, v) - mode-2 speed density v^3*exp(-((v-v0)/alpha)^2) (unnormalized;
   !              the rejection envelope f_star = dens(v_star) is its global maximum)
   pure real(8) function dens(b, v)
      type(beam_t), intent(in) :: b
      real(8), intent(in) :: v
      dens = v*v*v*exp(-((v - b%v0_c)/b%alpha_c)**2)
   end function dens

   ! named-abort helpers (name the caller; STOP 1)
   subroutine stop_beamx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_beamx

   subroutine stop_beamz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_beamz

   subroutine stop_beami(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_beami
end module beam_laws
