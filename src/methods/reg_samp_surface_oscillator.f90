!=====================================================================
! reg_samp_surface_oscillator.f90 - register the 'surface_oscillator'
!                                    distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the carrier's bracket parameters
!   (N_OSC/T_OSC/N_LEVEL/K_OSC - one row per carrier, the 2026-10-13
!   bracket ruling; absent names keep the member's field defaults) ->
!   the member's own init; an unselected member skips.
!=====================================================================
module reg_samp_surface_oscillator
   use sampler,                 only: samp_reg, samp_code
   use config,                  only: reactants
   use input,                   only: par_real, par_int
   use samp_surface_oscillator, only: surface_oscillator_draw, surface_oscillator_realize, &
                                      surf_osc_params_t, surface_oscillator_init
   implicit none
   private
   public :: reg_samp_surface_oscillator_all
   public :: reg_samp_surface_oscillator_init
contains
   !------------------------------------------------------------------
   ! reg_samp_surface_oscillator_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_surface_oscillator_all()
      call samp_reg('surface_oscillator', surface_oscillator_draw, surface_oscillator_realize)
   end subroutine reg_samp_surface_oscillator_all

   !------------------------------------------------------------------
   ! reg_samp_surface_oscillator_init() - assembly seam: per-carrier
   !                                       bracket keys (defaults first -
   !                                       absent names keep them) -> member init
   !------------------------------------------------------------------
   subroutine reg_samp_surface_oscillator_init()
      type(surf_osc_params_t) :: p     ! member-parameter receiver (field defaults)
      integer :: i, k, code, n_carry
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('surface_oscillator')
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) return
      allocate (p%n_osc(n_carry), p%t_osc(n_carry), p%n_level(n_carry), p%k_osc(n_carry))
      p%n_osc   = 0                       ! the member's field defaults (a bare
      p%t_osc   = 300.0d0                 ! word or an absent bracket name keeps
      p%n_level = 0                       ! them)
      p%k_osc   = 1.30092582d0
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         call par_int('SURFACE_OSCILLATOR', i, 'N_OSC',   p%n_osc(k))
         call par_real('SURFACE_OSCILLATOR', i, 'T_OSC',   p%t_osc(k))
         call par_int('SURFACE_OSCILLATOR', i, 'N_LEVEL', p%n_level(k))
         call par_real('SURFACE_OSCILLATOR', i, 'K_OSC',   p%k_osc(k))
      end do
      call surface_oscillator_init(p)
   end subroutine reg_samp_surface_oscillator_init
end module reg_samp_surface_oscillator
