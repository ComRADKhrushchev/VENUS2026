!=====================================================================
! reg_samp_barrier_excitation.f90 - register the 'barrier_excitation' member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the carrier's bracket parameters
!   (E_STAB/N_E_BAR/E_BAR/T_BAR - one row per carrier, the 2026-10-13
!   bracket ruling; absent names keep the member's field defaults) ->
!   the member's own init; an unselected member skips.
!=====================================================================
module reg_samp_barrier_excitation
   use sampler,                only: samp_reg, samp_code
   use config,                 only: reactants
   use input,                  only: par_real, par_int
   use samp_barrier_excitation, only: barrier_excitation_draw, barrier_excitation_realize, &
                                      barrier_excitation_params_t, barrier_excitation_init
   implicit none
   private
   public :: reg_samp_barrier_excitation_all
   public :: reg_samp_barrier_excitation_init
contains
   !------------------------------------------------------------------
   ! reg_samp_barrier_excitation_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_barrier_excitation_all()
      call samp_reg('barrier_excitation', barrier_excitation_draw, barrier_excitation_realize)
   end subroutine reg_samp_barrier_excitation_all

   !------------------------------------------------------------------
   ! reg_samp_barrier_excitation_init() - assembly seam: per-carrier bracket
   !                                       keys (defaults first - absent names
   !                                       keep them) -> member init
   !------------------------------------------------------------------
   subroutine reg_samp_barrier_excitation_init()
      type(barrier_excitation_params_t) :: p   ! member-parameter receiver (field defaults)
      integer :: i, k, code
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('barrier_excitation')
      if (count(reactants%dist_scheme == code) == 0) return
      allocate (p%e_stab(count(reactants%dist_scheme == code)), &
                p%n_e_bar(count(reactants%dist_scheme == code)), &
                p%e_bar(count(reactants%dist_scheme == code)), &
                p%t_bar(count(reactants%dist_scheme == code)))
      p%e_stab = 0.34691355d0           ! the member's field defaults (the
      p%n_e_bar = 0                     ! pre-2026-10-06 kcal-era values
      p%e_bar = 0.21682097d0            ! re-expressed; absent names keep them)
      p%t_bar = 900.0d0
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         call par_real('BARRIER_EXCITATION', i, 'E_STAB',  p%e_stab(k))
         call par_int('BARRIER_EXCITATION', i, 'N_E_BAR',  p%n_e_bar(k))
         call par_real('BARRIER_EXCITATION', i, 'E_BAR',   p%e_bar(k))
         call par_real('BARRIER_EXCITATION', i, 'T_BAR',   p%t_bar(k))
      end do
      call barrier_excitation_init(p)
   end subroutine reg_samp_barrier_excitation_init
end module reg_samp_barrier_excitation
