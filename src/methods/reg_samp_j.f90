!=====================================================================
! reg_samp_j.f90 - register the 'j' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: buffered input keys -> params receiver
!   (absent keys keep field defaults) -> the member's own init; an
!   unselected member skips.
!=====================================================================
module reg_samp_j
   use sampler, only: samp_reg, samp_code
   use config,  only: reactants
   use input,   only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast, only: seam_real
   use samp_j,  only: j_draw, j_realize, j_params_t, j_init
   implicit none
   private
   public :: reg_samp_j_all
   public :: reg_samp_j_init
contains
   !------------------------------------------------------------------
   ! reg_samp_j_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_j_all()
      call samp_reg('j', j_draw, j_realize)
   end subroutine reg_samp_j_all

   !------------------------------------------------------------------
   ! reg_samp_j_init() - assembly seam: buffered keys -> receiver -> member init
   subroutine reg_samp_j_init()
      type(j_params_t) :: p            ! member-parameter receiver (field defaults)
      character(len=512) :: sval       ! one buffered raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('j')) == 0) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('T_ROT'); call seam_real(buffer_key(i), sval, buffer_line(i), p%t_rot)
         end select
      end do
      call j_init(p)
   end subroutine reg_samp_j_init
end module reg_samp_j
