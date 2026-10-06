!=====================================================================
! reg_samp_glo_target.f90 - register the 'glo_target' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: buffered input keys -> params receiver
!   (absent keys keep field defaults) -> the member's own init (the
!   driver, input and member know nothing of each other); an
!   unselected member skips.
!=====================================================================
module reg_samp_glo_target
   use sampler,         only: samp_reg, samp_code
   use config,          only: reactants
   use input,           only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,       only: seam_real
   use samp_glo_target, only: glo_target_draw, glo_target_realize, glo_target_params_t, glo_target_init
   implicit none
   private
   public :: reg_samp_glo_target_all
   public :: reg_samp_glo_target_init
contains
   !------------------------------------------------------------------
   ! reg_samp_glo_target_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_glo_target_all()
      call samp_reg('glo_target', glo_target_draw, glo_target_realize)
   end subroutine reg_samp_glo_target_all

   !------------------------------------------------------------------
   ! reg_samp_glo_target_init() - assembly seam: buffered keys -> receiver ->
   !                              member init
   subroutine reg_samp_glo_target_init()
      type(glo_target_params_t) :: p     ! member-parameter receiver (field defaults)
      character(len=512) :: sval         ! one buffered raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('glo_target')) == 0) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('T_GLO');     call seam_real(buffer_key(i), sval, buffer_line(i), p%t_glo)
         case ('GAMMA_GLO'); call seam_real(buffer_key(i), sval, buffer_line(i), p%gamma_glo)
         case ('W_GHOST');   call seam_real(buffer_key(i), sval, buffer_line(i), p%w_ghost)
         end select
      end do
      call glo_target_init(p)
   end subroutine reg_samp_glo_target_init
end module reg_samp_glo_target
