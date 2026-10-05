!=====================================================================
! reg_inc_surface.f90 - register the 'incident_surface' member (guarded:
!   the surface paradigm only)
! Design:
!   One samp_reg call fills the row when the run derived the surface
!   paradigm (a cell-carrying system file) - the paradigm-guarded
!   registration keeps exactly ONE incident member armed per process
!   (the members= witness and the archives stay stable). _init is the
!   assembly seam: buffered beam keys (+ the fixed-site AIM_X/AIM_Y
!   pair, presence-flagged) -> params receiver -> the member's own init.
!=====================================================================
module reg_inc_surface
   use sampler,       only: samp_reg
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use input,         only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,     only: seam_real, seam_int
   use beam_laws,     only: beam_params_t
   use inc_surface,   only: inc_surface_draw, inc_surface_realize, inc_surface_init
   implicit none
   private
   public :: reg_inc_surface_all
   public :: reg_inc_surface_init
contains
   !------------------------------------------------------------------
   ! reg_inc_surface_all() - fill the member-table row (surface paradigm only)
   !------------------------------------------------------------------
   subroutine reg_inc_surface_all()
      if (.not. armed_here()) return
      call samp_reg('incident_surface', inc_surface_draw, inc_surface_realize)
   end subroutine reg_inc_surface_all

   !------------------------------------------------------------------
   ! reg_inc_surface_init() - assembly seam: buffered keys -> receiver ->
   !                          member init (surface paradigm only)
   !------------------------------------------------------------------
   subroutine reg_inc_surface_init()
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
   end subroutine reg_inc_surface_init

   ! armed_here() - the paradigm guard (surface only; the atom list must
   !                be assembled - the mirror of the dispatch rule)
   logical function armed_here()
      armed_here = allocated(list_atoms%frag) .and. reactants%surface_model /= 0
   end function armed_here
end module reg_inc_surface
