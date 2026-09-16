!=====================================================================
! dist_incident.f90 - distribution member: incident-channel collision
!   prescription (energy, impact parameter, beam orientation, drop
!   point), sampled after the per-fragment members write the interiors
! Design:
!   One placement formula for both paradigms (surface superset):
!   projectile COM = target site + b offset + r_sep back-off; interiors
!   ride along via COM shifts, momenta are additive. e_rel / b_max /
!   r_sep are system-level prescriptions from config%reactants - member
!   parameters carry only distribution shapes.
!   Units: angles [rad]; t_trans [K]; e_rel [kcal/mol] via e_conv; speeds [Ang/(10 fs)]; p [amu*Ang/(10 fs)].
!=====================================================================
module dist_incident
   use consts,        only: pi, two_pi, e_conv, r_kcal
   use rng,           only: rng_u, rng_gamma
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   implicit none
   private
   public :: inc_params_t, inc_init, inc_sample

   integer, parameter :: max_reject = 100000  ! rejection-sampling draw cap [-]

   ! Member-owned parameters (the collision energy / impact parameter / separation
   ! are system-level prescriptions in config%reactants - not duplicated here)
   type :: inc_params_t
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
                                   ! unit cell, 2 = uniform in the cluster box; the box
                                   ! half-width is roster%box_pad, the cell edges/skew are
                                   ! roster fields)
   end type

   ! Derived member state (module-private; rebuilt by every inc_init)
   logical :: inited = .false.     ! init completed [flag]
   logical :: surf_case = .false.  ! surface paradigm (surface_model /= 0) [flag]
   integer :: mode_e = 1           ! collision-energy mode copy [-]
   integer :: mode_th = 0          ! polar-angle mode copy [-]
   integer :: mode_ch = 0          ! azimuth mode copy [-]
   integer :: mode_b = 0           ! aiming mode copy [-]
   integer :: mode_aim = 0         ! drop-point mode copy [-]
   real(8) :: cos_tmax = 1.0d0     ! cos(thta_max) (fixed polar value / sampling bound) [-]
   real(8) :: chi_fix = 0.0d0      ! fixed azimuth [rad]
   real(8) :: t_trans_c = 0.0d0    ! translational temperature copy [K]
   real(8) :: e_rel_i = 0.0d0      ! collision energy [internal] (modes 1 and 2)
   real(8) :: b_max_c = 0.0d0      ! collision-disk radius [Angstrom]
   real(8) :: r_sep_c = 0.0d0      ! drop height / beam back-off scale [Angstrom]
   real(8) :: mu_eff = 0.0d0       ! effective mass of the relative motion [amu]
   real(8) :: w1_frac = 0.0d0      ! projectile COM speed fraction m2/(m1+m2) [-]
   real(8) :: w2_frac = 0.0d0      ! target COM speed fraction m1/(m1+m2) [-]
   real(8) :: v0_c = 0.0d0         ! mode-2 center speed [Angstrom/(10 fs)]
   real(8) :: alpha_c = 0.0d0      ! mode-2 velocity width [Angstrom/(10 fs)]
   real(8) :: v_star = 0.0d0       ! mode-2 density peak [Angstrom/(10 fs)]
   real(8) :: v_hi = 0.0d0         ! mode-2 proposal upper bound [Angstrom/(10 fs)]
   real(8) :: f_star = 0.0d0       ! mode-2 envelope height DENS(v*) [-]
   real(8) :: a1_vec(2) = 0.0d0    ! unit-cell edge a1 [Angstrom]
   real(8) :: a2_vec(2) = 0.0d0    ! unit-cell edge a2 [Angstrom]
   real(8) :: pad_c = 0.0d0        ! cluster-box half-width [Angstrom]
   integer, allocatable :: at_proj(:)  ! projectile fragment atom list [global index]
   integer, allocatable :: at_targ(:)  ! target fragment atom list [global index] (gas pair)
   integer :: nat_hi = 0           ! highest atom index touched [global index]
contains
   !------------------------------------------------------------------
   ! inc_init(p) - validate the parameters and the collision geometry,
   !               precompute the mode-2 envelope and the cell vectors
   subroutine inc_init(p)
      type(inc_params_t), intent(in) :: p ! member parameters (mode switches / shape parameters)
      real(8) :: m1, m2, mt, dum
      integer :: nf

      ! 1. parameter validation (fail-loud, naming the field)
      if (p%n_e_rel < 0 .or. p%n_e_rel > 2) then
         call stop_inci('inc_init', 'collision-energy mode n_e_rel must be 0 (thermal), '// &
            '1 (fixed) or 2 (density-weighted), got', p%n_e_rel)
      end if
      if (p%n_e_rel == 0 .and. p%t_trans <= 0.0d0) then
         dum = p%t_trans
         call stop_incz('inc_init', 'thermal mode 0 requires a translational temperature '// &
            't_trans [K] > 0, got', dum)
      end if
      if (p%n_e_rel == 2 .and. p%v_width <= 0.0d0) then
         dum = p%v_width
         call stop_incz('inc_init', 'density-weighted mode 2 requires a velocity width '// &
            'v_width > 0, got', dum)
      end if
      if (p%n_thta < 0 .or. p%n_thta > 1) then
         call stop_inci('inc_init', 'polar-angle mode n_thta must be 0 (fixed) or 1 (sampled), got', p%n_thta)
      end if
      if (p%n_thta == 1) then
         if (p%thta_max <= 0.0d0 .or. p%thta_max > 0.5d0*pi) then
            dum = p%thta_max
            call stop_incz('inc_init', 'sampled polar bound thta_max [rad] must lie in '// &
               '(0, pi/2] (a degenerate or sub-surface range), got', dum)
         end if
      else
         if (p%thta_max < 0.0d0 .or. p%thta_max > 0.5d0*pi) then
            dum = p%thta_max
            call stop_incz('inc_init', 'fixed polar angle thta_max [rad] must lie in '// &
               '[0, pi/2], got', dum)
         end if
      end if
      if (p%n_chi < 0 .or. p%n_chi > 1) then
         call stop_inci('inc_init', 'azimuth mode n_chi must be 0 (fixed) or 1 (sampled), got', p%n_chi)
      end if
      if (p%n_b < 0 .or. p%n_b > 1) then
         call stop_inci('inc_init', 'aiming mode n_b must be 0 (fixed b_max) or 1 (collision disk), got', p%n_b)
      end if
      if (p%n_aim < 0 .or. p%n_aim > 2) then
         call stop_inci('inc_init', 'drop-point mode n_aim must be 0 (origin), 1 (unit cell) '// &
            'or 2 (cluster box), got', p%n_aim)
      end if
      if (p%n_aim == 1) then
         if (roster%box_lx <= 0.0d0 .or. roster%box_ly <= 0.0d0 .or. sin(roster%skew) <= 0.0d0) then
            call stop_incx('inc_init', 'drop-point mode 1 requires a unit cell in the roster '// &
               '(box_lx, box_ly, skew - none staged)')
         end if
      end if
      if (p%n_aim == 2 .and. roster%box_pad <= 0.0d0) then
         dum = roster%box_pad
         call stop_incz('inc_init', 'drop-point mode 2 requires a positive cluster-box '// &
            'half-width box_pad [Angstrom], got', dum)
      end if

      ! 2. collision geometry (roster fragments + paradigm)
      if (.not. allocated(roster%frag)) then
         call stop_incx('inc_init', 'no roster fragment table - the roster must be assembled '// &
            'before member init (assembly-order error)')
      end if
      nf = size(roster%frag)
      surf_case = (reactants%surface_model /= 0)
      if (surf_case) then
         if (nf /= 1) then
            call stop_inci('inc_init', 'a surface-paradigm incident channel carries exactly one '// &
               'projectile fragment (the surface is not a roster fragment); fragment count =', nf)
         end if
         at_proj = roster%frag(1)%list
         mu_eff = roster%frag(1)%mass
         w1_frac = 1.0d0
         w2_frac = 0.0d0
      else
         if (nf /= 2) then
            call stop_inci('inc_init', 'a gas-phase incident channel carries exactly two fragments '// &
               '(fragment 1 projectile, fragment 2 target); fragment count =', nf)
         end if
         at_proj = roster%frag(1)%list
         at_targ = roster%frag(2)%list
         m1 = roster%frag(1)%mass
         m2 = roster%frag(2)%mass
         mt = m1 + m2
         mu_eff = m1*m2/mt
         w1_frac = m2/mt
         w2_frac = m1/mt
      end if
      nat_hi = maxval(at_proj)
      if (allocated(at_targ)) nat_hi = max(nat_hi, maxval(at_targ))

      ! 3. system-level incident prescriptions (config%reactants, read once here)
      if ((p%n_e_rel == 1 .or. p%n_e_rel == 2) .and. reactants%e_rel <= 0.0d0) then
         dum = reactants%e_rel
         call stop_incz('inc_init', 'collision-energy modes 1 and 2 require a positive collision '// &
            'energy e_rel [kcal/mol] in the reactants provision, got', dum)
      end if
      if (reactants%b_max <= 0.0d0) then
         dum = reactants%b_max
         call stop_incz('inc_init', 'the collision-disk radius b_max [Angstrom] must be positive, got', dum)
      end if
      if (reactants%r_sep <= 0.0d0) then
         dum = reactants%r_sep
         call stop_incz('inc_init', 'the initial separation r_sep [Angstrom] must be positive, got', dum)
      end if
      e_rel_i = e_conv*reactants%e_rel
      b_max_c = reactants%b_max
      r_sep_c = reactants%r_sep

      ! 4. precompute (mode copies, trigonometry, mode-2 envelope, cell vectors)
      mode_e = p%n_e_rel
      mode_th = p%n_thta
      mode_ch = p%n_chi
      mode_b = p%n_b
      mode_aim = p%n_aim
      cos_tmax = cos(p%thta_max)
      chi_fix = p%chi
      t_trans_c = p%t_trans
      if (mode_e == 2) then
         alpha_c = p%v_width
         v0_c = sqrt(2.0d0*e_rel_i/mu_eff)
         v_star = 0.5d0*v0_c + sqrt(0.25d0*v0_c*v0_c + 1.5d0*alpha_c*alpha_c)
         v_hi = 5.0d0*v_star
         f_star = dens(v_star)
      end if
      a1_vec = (/ roster%box_lx, 0.0d0 /)
      a2_vec = roster%box_ly*(/ cos(roster%skew), sin(roster%skew) /)
      pad_c = roster%box_pad
      inited = .true.
   end subroutine inc_init

   !------------------------------------------------------------------
   ! inc_sample(sta) - draw the collision prescription and assemble the
   !                   placement + relative motion into the state
   subroutine inc_sample(sta)
      type(state_t), intent(inout) :: sta  ! initial state (fragment COMs are placed,
                                           ! translational momenta are added; the fragment
                                           ! interiors written by the per-fragment members
                                           ! ride along untouched)
      real(8) :: costh, sinth, chi_v, bval, phi, rx0, ry0, vv
      real(8) :: u_inc(3), p_site(3), vel(3)
      real(8) :: s1, s2
      integer :: n_try

      if (.not. inited) then
         call stop_incx('inc_sample', 'sample called before inc_init - the member carries '// &
            'no parameter set')
      end if
      if (size(sta%q) < 3*nat_hi .or. size(sta%p) < 3*nat_hi) then
         call stop_incx('inc_sample', 'state too small for the collision fragments - '// &
            'state_create must follow the roster')
      end if

      ! 1. beam orientation: cos(theta) (fixed or P(theta) prop. sin(theta)),
      !    azimuth chi (fixed or uniform)
      if (mode_th == 1) then
         costh = 1.0d0 - rng_u()*(1.0d0 - cos_tmax)
      else
         costh = cos_tmax
      end if
      sinth = sqrt(max(0.0d0, 1.0d0 - costh*costh))
      if (mode_ch == 1) then
         chi_v = two_pi*rng_u()
      else
         chi_v = chi_fix
      end if

      ! 2. aiming radius (fixed b_max or b^2 uniform) and disk azimuth phi
      !    (uniform, or chi in the fixed-b mode)
      if (mode_b == 1) then
         bval = b_max_c*sqrt(rng_u())
         phi = two_pi*rng_u()
      else
         bval = b_max_c
         phi = chi_v
      end if

      ! 3. drop point (origin / unit cell / cluster box)
      rx0 = 0.0d0
      ry0 = 0.0d0
      select case (mode_aim)
      case (1)
         s1 = rng_u()
         s2 = rng_u()
         rx0 = s1*a1_vec(1) + s2*a2_vec(1)
         ry0 = s1*a1_vec(2) + s2*a2_vec(2)
      case (2)
         rx0 = pad_c*(2.0d0*rng_u() - 1.0d0)
         ry0 = pad_c*(2.0d0*rng_u() - 1.0d0)
      end select

      ! 4. collision speed per energy mode
      select case (mode_e)
      case (0)
         vv = sqrt(2.0d0*e_conv*r_kcal*t_trans_c*rng_gamma(2)/mu_eff)
      case (1)
         vv = sqrt(2.0d0*e_rel_i/mu_eff)
      case (2)
         do n_try = 1, max_reject
            vv = v_hi*rng_u()
            if (rng_u()*f_star <= dens(vv)) exit
            if (n_try == max_reject) then
               call stop_incx('inc_sample', 'rejection sampling of the collision speed '// &
                  'exceeded its draw cap (the density ratio is inconsistent with the envelope)')
            end if
         end do
      end select

      ! 5. assembly: placement by COM shift, momenta additive. Incident
      !    direction points at the target; drop height = r_sep exactly;
      !    lateral beam back-off = r_sep*sin(theta).
      u_inc = (/ -sinth*cos(chi_v), -sinth*sin(chi_v), -costh /)
      p_site = (/ rx0 + bval*cos(phi) + r_sep_c*sinth*cos(chi_v), &
                  ry0 + bval*sin(phi) + r_sep_c*sinth*sin(chi_v), r_sep_c /)
      if (surf_case) then
         call com_shift(sta, at_proj, p_site)
         vel = vv*u_inc
         call add_mom(sta, at_proj, vel)
      else
         call com_shift(sta, at_targ, (/ rx0, ry0, 0.0d0 /))
         call com_shift(sta, at_proj, p_site)
         vel = w1_frac*vv*u_inc
         call add_mom(sta, at_proj, vel)
         vel = -w2_frac*vv*u_inc
         call add_mom(sta, at_targ, vel)
      end if
   end subroutine inc_sample

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! dens(v) - mode-2 speed density v^3*exp(-((v-v0)/alpha)^2) (unnormalized;
   !           the rejection envelope f_star = dens(v_star) is its global maximum)
   pure real(8) function dens(v)
      real(8), intent(in) :: v
      dens = v*v*v*exp(-((v - v0_c)/alpha_c)**2)
   end function dens

   ! com_shift(sta, list, com_new) - translate the fragment's atoms so their mass
   !                                center lands on com_new (internal geometry untouched)
   subroutine com_shift(sta, list, com_new)
      type(state_t), intent(inout) :: sta
      integer, intent(in) :: list(:)
      real(8), intent(in) :: com_new(3)
      real(8) :: com(3), shift(3)
      integer :: ia
      com = 0.0d0
      do ia = 1, size(list)
         com = com + roster%mass(list(ia))*sta%q(3*list(ia) - 2:3*list(ia))
      end do
      com = com/sum(roster%mass(list))
      shift = com_new - com
      do ia = 1, size(list)
         sta%q(3*list(ia) - 2:3*list(ia)) = sta%q(3*list(ia) - 2:3*list(ia)) + shift
      end do
   end subroutine com_shift

   ! add_mom(sta, list, vel) - add the fragment COM velocity vel to every atom's
   !                           momentum (translational momenta are additive)
   subroutine add_mom(sta, list, vel)
      type(state_t), intent(inout) :: sta
      integer, intent(in) :: list(:)
      real(8), intent(in) :: vel(3)
      integer :: ia, j
      do ia = 1, size(list)
         do j = 1, 3
            sta%p(3*list(ia) - 3 + j) = sta%p(3*list(ia) - 3 + j) + roster%mass(list(ia))*vel(j)
         end do
      end do
   end subroutine add_mom

   ! loud-stop helpers (name the caller; STOP 1)
   subroutine stop_incx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_incx

   subroutine stop_incz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_incz

   subroutine stop_inci(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_inci

end module dist_incident
