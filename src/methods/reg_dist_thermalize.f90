!=====================================================================
! reg_dist_thermalize.f90 - register the 'thermalize' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: staged input keys -> params receiver
!   (absent keys keep field defaults) -> the member's own init; an
!   unselected member skips.
!=====================================================================
module reg_dist_thermalize
   use sampler,         only: samp_reg, samp_code
   use config,          only: reactants
   use input,           only: staging_n_rows, staging_key, staging_val, staging_line
   use seam_cast,       only: seam_real, seam_int
   use dist_thermalize, only: thermalize_sample, thermalize_params_t, thermalize_init
   implicit none
   private
   public :: reg_dist_thermalize_all
   public :: reg_dist_thermalize_init
contains
   !------------------------------------------------------------------
   ! reg_dist_thermalize_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_thermalize_all()
      call samp_reg('thermalize', thermalize_sample)
   end subroutine reg_dist_thermalize_all

   !------------------------------------------------------------------
   ! reg_dist_thermalize_init() - assembly seam: staged keys -> receiver ->
   !                              member init
   subroutine reg_dist_thermalize_init()
      type(thermalize_params_t) :: p   ! member-parameter receiver (field defaults)
      character(len=512) :: sval       ! one staged raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('thermalize')) == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('T_EQ');  call seam_real(staging_key(i), sval, staging_line(i), p%t_eq)
         case ('N_EQ');  call seam_int(staging_key(i), sval, staging_line(i), p%n_eq)
         case ('DT_EQ'); call seam_real(staging_key(i), sval, staging_line(i), p%dt_eq)
         end select
      end do
      call thermalize_init(p)
   end subroutine reg_dist_thermalize_init
end module reg_dist_thermalize
