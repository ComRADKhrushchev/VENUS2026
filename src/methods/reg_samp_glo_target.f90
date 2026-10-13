!=====================================================================
! reg_samp_glo_target.f90 - register the 'glo_target' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the carrier's bracket parameters
!   (T_GLO/GAMMA_GLO/W_GHOST - one row per carrier, the 2026-10-13
!   bracket ruling; a bare word keeps the field defaults and the
!   member's own validation speaks) -> the member's own init (the
!   driver, input and member know nothing of each other); an
!   unselected member skips.
!=====================================================================
module reg_samp_glo_target
   use sampler,         only: samp_reg, samp_code
   use config,          only: reactants
   use input,           only: par_real
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
   ! reg_samp_glo_target_init() - assembly seam: per-carrier bracket keys ->
   !                              receiver -> member init
   !------------------------------------------------------------------
   subroutine reg_samp_glo_target_init()
      type(glo_target_params_t) :: p     ! member-parameter receiver (field defaults)
      integer :: i, k, code
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('glo_target')
      if (count(reactants%dist_scheme == code) == 0) return
      allocate (p%t_glo(count(reactants%dist_scheme == code)), &
                p%gamma_glo(count(reactants%dist_scheme == code)), &
                p%w_ghost(count(reactants%dist_scheme == code)))
      p%t_glo     = 0.0d0               ! absent bracket names keep the defaults;
      p%gamma_glo = 0.0d0               ! glo_target_init's own validation rejects
      p%w_ghost   = 0.0d0               ! them
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         call par_real('GLO_TARGET', i, 'T_GLO',     p%t_glo(k))
         call par_real('GLO_TARGET', i, 'GAMMA_GLO', p%gamma_glo(k))
         call par_real('GLO_TARGET', i, 'W_GHOST',   p%w_ghost(k))
      end do
      call glo_target_init(p)
   end subroutine reg_samp_glo_target_init
end module reg_samp_glo_target
