!   The member carries its own registration (member_register/ member_init - the
!   fixed-name entries of the 2026-10-14 weave batch; the retired
!   reg_reg_inc_surface.f90's bodies live here now, paradigm guard included).
!=====================================================================
! inc_surface.f90 - the surface-paradigm incident-channel member
!   (energy, impact parameter, beam orientation, drop point), sampled
!   after the per-fragment members write the interiors
! Design:
!   Projectile COM = drop site + b offset + r_sep back-off; the surface
!   absorbs recoil (the beam momentum is not split - the projectile COM
!   velocity is assigned in full, w = 1). Exactly one PROJECTILE fragment
!   (i_proj names it); TARGET (optional, the 2026-10-16 ruling) names the
!   explicit substrate fragment - the channel RESETS its atoms to the file
!   lattice positions (q0) before placing the projectile, closing the
!   placement responsibility chain the surface paradigm had lost; every
!   OTHER fragment is a substrate carrier (its own dist_scheme member
!   samples it; the incident channel never touches it) - the surface's
!   crystal frame is the POSCAR cell, a system-level fact, so a
!   non-TARGET substrate fragment is never re-placed or rotated by the
!   beam assembly. The drop site: origin (N_AIM = 0), uniform in the unit
!   cell (1), or the fixed site AIM_X/AIM_Y [Angstrom, POSCAR frame] (2);
!   B_MAX must NOT be written (ruled 2026-10-01 - lateral placement is
!   N_AIM's alone; the unwritten default 0 = the on-site beam). The A_LAT key
!   must be commensurate with the buffered cell (integer ratio, checked).
! Units: angles [rad]; t_trans [K]; e_rel [eV]; speeds [Ang/fs]; p [amu*Ang/fs].
!=====================================================================
module inc_surface
   use state,        only: state_t
   use config,       only: reactants
   use config_atoms, only: list_atoms
   use input,        only: buffer_key, buffer_line, buffer_n_rows, buffer_val
   use sampler, only: samp_reg
   use seam_cast, only: seam_int, seam_real
   use beam_laws,    only: beam_params_t, beam_t, beam_init, beam_draw, &
                           beam_com_shift, beam_set_com_vel
   use diag, only: diag_fatal, diag_cat_sampler
   implicit none
   private
   public :: inc_surface_init, inc_surface_draw, inc_surface_realize, inc_surface_sample, &
             member_register, member_init

   ! derived member state (module-private; rebuilt by every inc_surface_init)
   type(beam_t) :: beam                   ! the shared beam-law instance
   integer, allocatable :: at_proj(:)     ! projectile fragment atom list [global index]
   integer :: nat_hi = 0                  ! highest atom index touched [global index]
contains
   !------------------------------------------------------------------
   ! inc_surface_init(p) - validate the collision geometry (one projectile
   !                 fragment + substrate carriers), bind the projectile,
   !                 init the beam instance (on-site beam legal)
   !------------------------------------------------------------------
   subroutine inc_surface_init(p)
      type(beam_params_t), intent(in) :: p  ! beam shape parameters
      real(8) :: mu
      integer :: nf, i
      character(len=*), parameter :: here = 'inc_surface_init'

      ! 1. collision geometry (list_atoms fragments + the projectile role)
      if (.not. allocated(list_atoms%frag)) then
         call stop_surfx(here, 'no list_atoms fragment table - the atom list must be assembled '// &
            'before member init (assembly-order error)')
      end if
      nf = size(list_atoms%frag)
      if (reactants%surface_model == 0) then
         call stop_surfi(here, 'the surface incident channel requires the surface class '// &
            '(SCATTERING = SURFACE); surface_model =', reactants%surface_model)
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
            call stop_surfx(here, 'B_MAX must not be written in a surface-scattering input - '// &
               'the lateral placement belongs to the N_AIM drop-point modes (origin /'// &
               'unit cell / fixed site); the collision-disk radius is a gas-phase'// &
               'pair concept (raised 2026-09-28, ruled 2026-10-01)')
         end if
      end do

      ! 3. beam laws (the unwritten default b_max = 0 = the on-site beam;
      !    the drop-site frame is the CELL-derived box - the site statistics
      !    lattice itself retired with the A_LAT key, 2026-10-13 ruling)
      call beam_init(beam, p, mu, reactants%e_rel, reactants%b_max, reactants%r_sep, &
                     here, onsite_ok=.true.)
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
      integer :: fa, g_atom, i3
      if (.not. beam%inited) then
         call stop_surfx('inc_surface_realize', 'realize called before inc_surface_init - '// &
            'the member carries no parameter set')
      end if
      if (size(sta%q) < 3*nat_hi .or. size(sta%p) < 3*nat_hi) then
         call stop_surfx('inc_surface_realize', 'state too small for the collision fragments - '// &
            'state_create must follow the atom list')
      end if
      ! 1. TARGET placement (2026-10-16 ruling): the explicit substrate fragment's
      !    atoms reset to the file lattice positions (q0) - the placement
      !    responsibility chain returns to the surface paradigm. Ghost rows land on
      !    their declared positions; mobile rows of an un-bathed target start from
      !    the lattice - the bath provision and this reset together define the
      !    substrate's per-trajectory initial state
      if (list_atoms%i_targ >= 1 .and. list_atoms%i_targ <= size(list_atoms%frag)) then
         do fa = 1, list_atoms%frag(list_atoms%i_targ)%nat
            g_atom = list_atoms%frag(list_atoms%i_targ)%list(fa)
            i3 = 3*(g_atom - 1)
            sta%q(i3+1:i3+3) = list_atoms%q0(i3+1:i3+3)
         end do
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

   ! named-abort helpers (name the caller; abort via diag_fatal)
   subroutine stop_surfx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      call diag_fatal(diag_cat_sampler, 'inc_surface.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_surfx

   subroutine stop_surfi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      call diag_fatal(diag_cat_sampler, 'inc_surface.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_surfi
   !------------------------------------------------------------------
   ! member_register - the fixed-name weave entry (prebuild binds it into
   !                   registry_gen_all; the retired reg file's _all body,
   !                   paradigm guard included - exactly one incident member
   !                   is armed per process)
   !------------------------------------------------------------------
   subroutine member_register()
      if (.not. armed_here()) return
      call samp_reg('incident_surface', inc_surface_draw, inc_surface_realize)
   end subroutine member_register

   !------------------------------------------------------------------
   ! member_init - the fixed-name init seam (registry_gen_init calls it;
   !               the retired reg file's _init body moved verbatim)
   !------------------------------------------------------------------
   subroutine member_init()
      type(beam_params_t) :: p        ! parameter receiver (field defaults)
      character(len=512) :: sval      ! one buffered raw value
      integer :: i
      if (.not. armed_here()) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('N_E_REL');   call seam_int(buffer_key(i), sval, buffer_line(i), p%n_e_rel)
         case ('T_TRANS');   call seam_real(buffer_key(i), sval, buffer_line(i), p%t_trans)
         case ('V_WIDTH');   call seam_real(buffer_key(i), sval, buffer_line(i), p%v_width)
         case ('N_THTA');    call seam_int(buffer_key(i), sval, buffer_line(i), p%n_thta)
         case ('THTA_MAX');  call seam_real(buffer_key(i), sval, buffer_line(i), p%thta_max)
         case ('N_CHI');     call seam_int(buffer_key(i), sval, buffer_line(i), p%n_chi)
         case ('CHI');       call seam_real(buffer_key(i), sval, buffer_line(i), p%chi)
         case ('N_B');       call seam_int(buffer_key(i), sval, buffer_line(i), p%n_b)
         case ('N_AIM');     call seam_int(buffer_key(i), sval, buffer_line(i), p%n_aim)
         case ('AIM_X');     call seam_real(buffer_key(i), sval, buffer_line(i), p%aim_x)
                            p%aim_x_given = .true.
         case ('AIM_Y');     call seam_real(buffer_key(i), sval, buffer_line(i), p%aim_y)
                            p%aim_y_given = .true.
         end select
      end do
      call inc_surface_init(p)
   end subroutine member_init

   !------------------------------------------------------------------
   ! (moved verbatim from the retired reg file - the seam's own helper)
   !------------------------------------------------------------------
logical function armed_here()
      armed_here = allocated(list_atoms%frag) .and. reactants%surface_model /= 0
   end function armed_here

end module inc_surface
