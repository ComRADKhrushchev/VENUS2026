!=====================================================================
! reg_samp_j.f90 - register the 'j' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the carrier's bracket parameters
!   (T_ROT - one row per carrier, the 2026-10-13 bracket ruling; a bare
!   word keeps the field default and the member's own validation
!   speaks) -> the member's own init; an unselected member skips.
!=====================================================================
module reg_samp_j
   use sampler, only: samp_reg, samp_code
   use config,  only: reactants
   use input,   only: par_real
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
   ! reg_samp_j_init() - assembly seam: per-carrier bracket keys -> receiver
   !                     -> member init
   !------------------------------------------------------------------
   subroutine reg_samp_j_init()
      type(j_params_t) :: p            ! member-parameter receiver (field defaults)
      integer :: i, k, code
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('j')
      if (count(reactants%dist_scheme == code) == 0) return
      allocate (p%t_rot(count(reactants%dist_scheme == code)))
      p%t_rot = 0.0d0                   ! absent bracket name keeps the default;
                                        ! j_init's own validation rejects it
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         call par_real('J', i, 'T_ROT', p%t_rot(k))
      end do
      call j_init(p)
   end subroutine reg_samp_j_init
end module reg_samp_j
