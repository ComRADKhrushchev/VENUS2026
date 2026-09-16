!=====================================================================
! reg_dist_rotation.f90 - register the 'rotation' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: staged input keys -> params receiver
!   (absent keys keep field defaults) -> the member's own init; an
!   unselected member skips.
!=====================================================================
module reg_dist_rotation
   use sampler,       only: samp_reg, samp_code
   use config,        only: reactants
   use input,         only: staging_n_rows, staging_key, staging_val, staging_line
   use seam_cast,     only: seam_int
   use dist_rotation, only: rotation_sample, rotation_params_t, rotation_init
   implicit none
   private
   public :: reg_dist_rotation_all
   public :: reg_dist_rotation_init
contains
   !------------------------------------------------------------------
   ! reg_dist_rotation_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_rotation_all()
      call samp_reg('rotation', rotation_sample)
   end subroutine reg_dist_rotation_all

   !------------------------------------------------------------------
   ! reg_dist_rotation_init() - assembly seam: staged keys -> receiver ->
   !                            member init
   subroutine reg_dist_rotation_init()
      type(rotation_params_t) :: p     ! member-parameter receiver (field defaults)
      character(len=512) :: sval       ! one staged raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('rotation')) == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('J_ROT'); call seam_int(staging_key(i), sval, staging_line(i), p%j_rot)
         end select
      end do
      call rotation_init(p)
   end subroutine reg_dist_rotation_init
end module reg_dist_rotation
