!=====================================================================
! inc_pair.f90 - the two-reactant (gas-pair) incident-channel member
!   (energy, impact parameter, beam orientation, drop point), sampled
!   after the per-fragment members write the interiors
! Design:
!   Projectile COM = target site + b offset + r_sep back-off; interiors
!   ride along via COM shifts, fragment COM velocities are ASSIGNED (an
!   add would pile the beam momentum onto a monoatomic projectile's stale
!   carry once per trajectory) with the reduced-mass split
!   w1 = m2/(m1+m2) on the projectile and w2 = m1/(m1+m2) on the target.
!   Exactly two fragments, both named by the PROJECTILE / TARGET role
!   keys. The fixed-site drop mode (N_AIM = 2) is a surface-paradigm
!   option - here the drop point is the origin (0) or the collision-disk
!   anchor per N_B (1).
! Units: angles [rad]; t_trans [K]; e_rel [kcal/mol] via e_conv; speeds [Ang/(10 fs)]; p [amu*Ang/(10 fs)].
!=====================================================================
module inc_pair
   use state,        only: state_t
   use config,       only: reactants
   use config_atoms, only: list_atoms
   use beam_laws,    only: beam_params_t, beam_t, beam_init, beam_draw, &
                           beam_com_shift, beam_set_com_vel
   implicit none
   private
   public :: inc_pair_init, inc_pair_draw, inc_pair_realize, inc_pair_sample

   ! derived member state (module-private; rebuilt by every inc_pair_init)
   type(beam_t) :: beam                   ! the shared beam-law instance
   real(8) :: w1_frac = 0.0d0             ! projectile COM speed fraction m2/(m1+m2) [-]
   real(8) :: w2_frac = 0.0d0             ! target COM speed fraction m1/(m1+m2) [-]
   integer, allocatable :: at_proj(:)     ! projectile fragment atom list [global index]
   integer, allocatable :: at_targ(:)     ! target fragment atom list [global index]
   integer :: nat_hi = 0                  ! highest atom index touched [global index]
contains
   !------------------------------------------------------------------
   ! inc_pair_init(p) - validate the collision geometry (two role-named
   !                 fragments), bind both, init the beam instance with
   !                 the reduced mass
   !------------------------------------------------------------------
   subroutine inc_pair_init(p)
      type(beam_params_t), intent(in) :: p  ! beam shape parameters
      real(8) :: m1, m2, mt, mu
      integer :: nf
      character(len=*), parameter :: here = 'inc_pair_init'

      ! 1. collision geometry (list_atoms fragments + both roles)
      if (.not. allocated(list_atoms%frag)) then
         call stop_pairx(here, 'no list_atoms fragment table - the atom list must be assembled '// &
            'before member init (assembly-order error)')
      end if
      nf = size(list_atoms%frag)
      if (reactants%surface_model /= 0) then
         call stop_pairx(here, 'the pair incident channel requires the molecular paradigm '// &
            '(no cell-carrying system file) - the surface member owns the surface paradigm')
      end if
      if (nf /= 2) then
         call stop_pairi(here, 'a gas-phase incident channel carries exactly two fragments '// &
            '(the PROJECTILE and TARGET role keys name them); fragment count =', nf)
      end if
      if (list_atoms%i_proj < 1 .or. list_atoms%i_proj > nf .or. list_atoms%i_targ < 1 .or. &
          list_atoms%i_targ > nf .or. list_atoms%i_proj == list_atoms%i_targ) then
         call stop_pairx(here, 'no resolved projectile/target pair (list_atoms%i_proj / '// &
            'i_targ unset, out of range, or equal) - the list_atoms_load role reconcile must '// &
            'resolve both role keys before member init (assembly-order error)')
      end if
      if (p%n_aim == 2) then
         call stop_pairx(here, 'drop-point mode 2 (the fixed site AIM_X/AIM_Y) is a '// &
            'surface-paradigm option - the pair member drops at the origin (0) or per N_B (1)')
      end if
      at_proj = list_atoms%frag(list_atoms%i_proj)%list
      at_targ = list_atoms%frag(list_atoms%i_targ)%list
      m1 = list_atoms%frag(list_atoms%i_proj)%mass
      m2 = list_atoms%frag(list_atoms%i_targ)%mass
      mt = m1 + m2
      mu = m1*m2/mt
      w1_frac = m2/mt
      w2_frac = m1/mt
      nat_hi = max(maxval(at_proj), maxval(at_targ))

      ! 2. beam laws (the pair keeps the positive-b_max statute)
      call beam_init(beam, p, mu, reactants%e_rel, reactants%b_max, reactants%r_sep, here)
   end subroutine inc_pair_init

   !------------------------------------------------------------------
   ! inc_pair_draw() - draw the whole collision prescription (the shared
   !                 historical order); no state access
   !------------------------------------------------------------------
   subroutine inc_pair_draw()
      call beam_draw(beam)
   end subroutine inc_pair_draw

   !------------------------------------------------------------------
   ! inc_pair_realize(sta) - assembly from the drawn prescription: both
   !                 COMs placed (target at the drop anchor, projectile at
   !                 the backed-off site), velocities assigned with the
   !                 reduced-mass split. No RNG draws
   !------------------------------------------------------------------
   subroutine inc_pair_realize(sta)
      type(state_t), intent(inout) :: sta
      real(8) :: u_inc(3), p_site(3), vel(3)
      if (.not. beam%inited) then
         call stop_pairx('inc_pair_realize', 'realize called before inc_pair_init - '// &
            'the member carries no parameter set')
      end if
      if (size(sta%q) < 3*nat_hi .or. size(sta%p) < 3*nat_hi) then
         call stop_pairx('inc_pair_realize', 'state too small for the collision fragments - '// &
            'state_create must follow the atom list')
      end if
      u_inc = (/ -beam%sinth_s*cos(beam%chi_s), -beam%sinth_s*sin(beam%chi_s), -beam%costh_s /)
      p_site = (/ beam%rx0_s + beam%bval_s*cos(beam%phi_s) + beam%r_sep_c*beam%sinth_s*cos(beam%chi_s), &
                  beam%ry0_s + beam%bval_s*sin(beam%phi_s) + beam%r_sep_c*beam%sinth_s*sin(beam%chi_s), &
                  beam%r_sep_c /)
      call beam_com_shift(sta, at_targ, (/ beam%rx0_s, beam%ry0_s, 0.0d0 /))
      call beam_com_shift(sta, at_proj, p_site)
      vel = w1_frac*beam%vv_s*u_inc
      call beam_set_com_vel(sta, at_proj, vel)
      vel = -w2_frac*beam%vv_s*u_inc
      call beam_set_com_vel(sta, at_targ, vel)
   end subroutine inc_pair_realize

   !------------------------------------------------------------------
   ! inc_pair_sample(sta) - fused convenience: draw then realize (the
   !                 check-program entry)
   !------------------------------------------------------------------
   subroutine inc_pair_sample(sta)
      type(state_t), intent(inout) :: sta
      call inc_pair_draw()
      call inc_pair_realize(sta)
   end subroutine inc_pair_sample

   ! named-abort helpers (name the caller; STOP 1)
   subroutine stop_pairx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_pairx

   subroutine stop_pairi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_pairi
end module inc_pair
