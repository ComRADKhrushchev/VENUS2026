!=====================================================================
! reg_samp_rotation.f90 - register the 'rotation' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the carrier's bracket parameters
!   (J_ROT - one row per carrier, the 2026-10-13 bracket ruling; a bare
!   word keeps the field default) -> the member's own init; an
!   unselected member skips.
!=====================================================================
module reg_samp_rotation
   use sampler,       only: samp_reg, samp_code
   use config,        only: reactants
   use input,         only: par_int
   use samp_rotation, only: rotation_draw, rotation_realize, rotation_params_t, rotation_init
   implicit none
   private
   public :: reg_samp_rotation_all
   public :: reg_samp_rotation_init
contains
   !------------------------------------------------------------------
   ! reg_samp_rotation_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_rotation_all()
      call samp_reg('rotation', rotation_draw, rotation_realize)
   end subroutine reg_samp_rotation_all

   !------------------------------------------------------------------
   ! reg_samp_rotation_init() - assembly seam: per-carrier bracket keys ->
   !                            receiver -> member init
   !------------------------------------------------------------------
   subroutine reg_samp_rotation_init()
      type(rotation_params_t) :: p     ! member-parameter receiver (field defaults)
      integer :: i, k, code
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('rotation')
      if (count(reactants%dist_scheme == code) == 0) return
      allocate (p%j_rot(count(reactants%dist_scheme == code)))
      p%j_rot = 0                        ! absent bracket name keeps the default
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         call par_int('ROTATION', i, 'J_ROT', p%j_rot(k))
      end do
      call rotation_init(p)
   end subroutine reg_samp_rotation_init
end module reg_samp_rotation
