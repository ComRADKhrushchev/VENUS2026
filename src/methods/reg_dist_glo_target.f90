!=====================================================================
! reg_dist_glo_target.f90 - register the 'glo_target' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: staged input keys -> params receiver
!   (absent keys keep field defaults) -> the member's own init (the
!   driver, input and member know nothing of each other); an
!   unselected member skips.
!=====================================================================
module reg_dist_glo_target
   use sampler,         only: samp_reg, samp_code
   use config,          only: reactants
   use input,           only: staging_n_rows, staging_key, staging_val, staging_line
   use seam_cast,       only: seam_real
   use dist_glo_target, only: glo_target_sample, glo_target_params_t, glo_target_init
   implicit none
   private
   public :: reg_dist_glo_target_all
   public :: reg_dist_glo_target_init
contains
   !------------------------------------------------------------------
   ! reg_dist_glo_target_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_glo_target_all()
      call samp_reg('glo_target', glo_target_sample)
   end subroutine reg_dist_glo_target_all

   !------------------------------------------------------------------
   ! reg_dist_glo_target_init() - assembly seam: staged keys -> receiver ->
   !                              member init
   subroutine reg_dist_glo_target_init()
      type(glo_target_params_t) :: p     ! member-parameter receiver (field defaults)
      character(len=512) :: sval         ! one staged raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('glo_target')) == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('T_GLO');     call seam_real(staging_key(i), sval, staging_line(i), p%t_glo)
         case ('GAMMA_GLO'); call seam_real(staging_key(i), sval, staging_line(i), p%gamma_glo)
         case ('W_GHOST');   call seam_real(staging_key(i), sval, staging_line(i), p%w_ghost)
         end select
      end do
      call glo_target_init(p)
   end subroutine reg_dist_glo_target_init
end module reg_dist_glo_target
