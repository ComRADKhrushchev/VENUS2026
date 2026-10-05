!=====================================================================
! inc_surface.f90 - the surface-paradigm incident-channel member
!   (energy, impact parameter, beam orientation, drop point), sampled
!   after the per-fragment members write the interiors
! Design:
!   Projectile COM = drop site + b offset + r_sep back-off; the surface
!   absorbs recoil (the beam momentum is not split - the projectile COM
!   velocity is assigned in full, w = 1). Exactly one PROJECTILE fragment
!   (i_proj names it); every OTHER fragment is a substrate carrier (its
!   own dist_scheme member samples it; the incident channel never touches
!   it) - the surface's crystal frame is the POSCAR cell, a system-level
!   fact, so a substrate fragment is never re-placed or rotated by the
!   beam assembly. The drop site: origin (N_AIM = 0), uniform in the unit
!   cell (1), or the fixed site AIM_X/AIM_Y [Angstrom, POSCAR frame] (2);
!   B_MAX must NOT be written (ruled 2026-10-01 - lateral placement is
!   N_AIM's alone; the unwritten default 0 = the on-site beam). The A_LAT key
!   must be commensurate with the buffered cell (integer ratio, checked).
! Units: angles [rad]; t_trans [K]; e_rel [kcal/mol] via e_conv; speeds [Ang/(10 fs)]; p [amu*Ang/(10 fs)].
!=====================================================================
module inc_surface
   use state,        only: state_t
   use config,       only: reactants
   use config_atoms, only: list_atoms
   use input,        only: buffer_n_rows, buffer_key, buffer_line
   use beam_laws,    only: beam_params_t, beam_t, beam_init, beam_draw, &
                           beam_com_shift, beam_set_com_vel
   implicit none
   private
   public :: inc_surface_init, inc_surface_draw, inc_surface_realize, inc_surface_sample

   ! derived member state (module-private; rebuilt by every inc_surface_init)
   type(beam_t) :: beam                   ! the shared beam-law instance
   integer, allocatable :: at_proj(:)     ! projectile fragment atom list [global index]
   integer :: nat_hi = 0                  ! highest atom index touched [global index]
contains
   !------------------------------------------------------------------
   ! inc_surface_init(p) - validate the collision geometry (one projectile
   !                 fragment + substrate carriers), bind the projectile,
   !                 init the beam instance (on-site beam legal), check the
   !                 A_LAT/cell commensurability
   !------------------------------------------------------------------
   subroutine inc_surface_init(p)
      type(beam_params_t), intent(in) :: p  ! beam shape parameters
      real(8) :: mu, rlx, rly
      integer :: nf, i
      character(len=*), parameter :: here = 'inc_surface_init'

      ! 1. collision geometry (list_atoms fragments + the projectile role)
      if (.not. allocated(list_atoms%frag)) then
         call stop_surfx(here, 'no list_atoms fragment table - the atom list must be assembled '// &
            'before member init (assembly-order error)')
      end if
      nf = size(list_atoms%frag)
      if (reactants%surface_model == 0) then
         call stop_surfi(here, 'the surface incident channel requires the surface paradigm '// &
            '(a cell-carrying system file); surface_model =', reactants%surface_model)
      end if
      if (list_atoms%i_proj < 1 .or. list_atoms%i_proj > nf) then
         call stop_surfx(here, 'no resolved projectile role (list_atoms%i_proj unset or '// &
            'out of range) - the list_atoms_load role reconcile must resolve the '// &
            'PROJECTILE key before member init (assembly-order error)')
      end if
      at_proj = list_atoms%frag(list_atoms%i_proj)%list
      mu = list_atoms%frag(list_atoms%i_proj)%mass
      nat_hi = maxval(at_proj)

      ! 2. B_MAX interception (raised 2026-09-28, ruled 2026-10-01): a
      !    periodic surface's lateral placement is owned EXCLUSIVELY by the
      !    N_AIM drop-point modes (origin / uniform in the unit cell / the
      !    fixed site); the collision-disk radius is a gas-phase pair
      !    concept and must not appear in a surface input at all
      do i = 1, buffer_n_rows()
         if (trim(buffer_key(i)) == 'B_MAX' .and. buffer_line(i) > 0) then
            call stop_surfx(here, 'B_MAX must not be written in a surface-paradigm input - '// &
               'the lateral placement belongs to the N_AIM drop-point modes (origin /'// &
               'unit cell / fixed site); the collision-disk radius is a gas-phase'// &
               'pair concept (raised 2026-09-28, ruled 2026-10-01)')
         end if
      end do

      ! 3. beam laws (the unwritten default b_max = 0 = the on-site beam)
      call beam_init(beam, p, mu, reactants%e_rel, reactants%b_max, reactants%r_sep, &
                     here, onsite_ok=.true.)

      ! 3. A_LAT / cell commensurability (the lattice constant must divide
      !    both in-plane edges of the buffered cell - the histogram lattice
      !    and the drop-point cell must describe the same surface)
      if (list_atoms%a_lat <= 0.0d0) then
         call stop_surfx(here, 'no surface lattice constant (list_atoms%a_lat unset) - '// &
            'the surface paradigm requires the A_LAT key')
      end if
      rlx = list_atoms%box_lx/list_atoms%a_lat
      rly = list_atoms%box_ly/list_atoms%a_lat
      if (abs(rlx - nint(rlx)) > 1.0d-6*max(1.0d0, abs(rlx)) .or. &
          abs(rly - nint(rly)) > 1.0d-6*max(1.0d0, abs(rly))) then
         call stop_surfx(here, 'the buffered cell is not commensurate with A_LAT '// &
            '(box_lx/A_LAT, box_ly/A_LAT must be integers)')
      end if
   end subroutine inc_surface_init

   !------------------------------------------------------------------
   ! inc_surface_draw() - draw the whole collision prescription (the shared
   !                 historical order); no state access
   !------------------------------------------------------------------
   subroutine inc_surface_draw()
      call beam_draw(beam)
   end subroutine inc_surface_draw

   !------------------------------------------------------------------
   ! inc_surface_realize(sta) - assembly from the drawn prescription: the
   !                 projectile COM is placed at the site, its COM velocity
   !                 assigned in full (the surface absorbs the recoil);
   !                 substrate fragments ride untouched. No RNG draws
   !------------------------------------------------------------------
   subroutine inc_surface_realize(sta)
      type(state_t), intent(inout) :: sta
      real(8) :: u_inc(3), p_site(3)
      if (.not. beam%inited) then
         call stop_surfx('inc_surface_realize', 'realize called before inc_surface_init - '// &
            'the member carries no parameter set')
      end if
      if (size(sta%q) < 3*nat_hi .or. size(sta%p) < 3*nat_hi) then
         call stop_surfx('inc_surface_realize', 'state too small for the collision fragments - '// &
            'state_create must follow the atom list')
      end if
      u_inc = (/ -beam%sinth_s*cos(beam%chi_s), -beam%sinth_s*sin(beam%chi_s), -beam%costh_s /)
      p_site = (/ beam%rx0_s + beam%bval_s*cos(beam%phi_s) + beam%r_sep_c*beam%sinth_s*cos(beam%chi_s), &
                  beam%ry0_s + beam%bval_s*sin(beam%phi_s) + beam%r_sep_c*beam%sinth_s*sin(beam%chi_s), &
                  beam%r_sep_c /)
      call beam_com_shift(sta, at_proj, p_site)
      call beam_set_com_vel(sta, at_proj, beam%vv_s*u_inc)
   end subroutine inc_surface_realize

   !------------------------------------------------------------------
   ! inc_surface_sample(sta) - fused convenience: draw then realize (the
   !                 check-program entry)
   !------------------------------------------------------------------
   subroutine inc_surface_sample(sta)
      type(state_t), intent(inout) :: sta
      call inc_surface_draw()
      call inc_surface_realize(sta)
   end subroutine inc_surface_sample

   ! named-abort helpers (name the caller; STOP 1)
   subroutine stop_surfx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_surfx

   subroutine stop_surfi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_surfi
end module inc_surface
