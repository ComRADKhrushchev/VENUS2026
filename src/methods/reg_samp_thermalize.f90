!=====================================================================
! reg_samp_thermalize.f90 - register the 'thermalize' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the carrier's bracket parameters
!   (T_EQ/N_EQ/DT_EQ - one row per carrier, the 2026-10-13 bracket
!   ruling; a bare word keeps the field defaults and the member's own
!   validation speaks) -> the member's own init; an unselected member
!   skips.
!=====================================================================
module reg_samp_thermalize
   use sampler,         only: samp_reg, samp_code
   use config,          only: reactants
   use input,           only: par_real, par_int
   use samp_thermalize, only: thermalize_draw, thermalize_realize, thermalize_params_t, thermalize_init
   implicit none
   private
   public :: reg_samp_thermalize_all
   public :: reg_samp_thermalize_init
contains
   !------------------------------------------------------------------
   ! reg_samp_thermalize_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_thermalize_all()
      call samp_reg('thermalize', thermalize_draw, thermalize_realize)
   end subroutine reg_samp_thermalize_all

   !------------------------------------------------------------------
   ! reg_samp_thermalize_init() - assembly seam: per-carrier bracket keys ->
   !                              receiver -> member init (the member owns the
   !                              one-shared-run constraint on N_EQ/DT_EQ)
   !------------------------------------------------------------------
   subroutine reg_samp_thermalize_init()
      type(thermalize_params_t) :: p   ! member-parameter receiver (field defaults)
      integer :: i, k, code
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('thermalize')
      if (count(reactants%dist_scheme == code) == 0) return
      allocate (p%t_eq(count(reactants%dist_scheme == code)), &
                p%n_eq(count(reactants%dist_scheme == code)), &
                p%dt_eq(count(reactants%dist_scheme == code)))
      p%t_eq  = 0.0d0                    ! absent bracket names keep the defaults;
      p%n_eq  = 0                         ! thermalize_init's own validation rejects
      p%dt_eq = 0.0d0                     ! them (and enforces the shared-run rule)
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         call par_real('THERMALIZE', i, 'T_EQ',  p%t_eq(k))
         call par_int('THERMALIZE', i, 'N_EQ',  p%n_eq(k))
         call par_real('THERMALIZE', i, 'DT_EQ', p%dt_eq(k))
      end do
      call thermalize_init(p)
   end subroutine reg_samp_thermalize_init
end module reg_samp_thermalize
