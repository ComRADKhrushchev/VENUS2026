!=====================================================================
! reg_samp_surface_oscillator.f90 - register the 'surface_oscillator'
!                                    distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: buffered input keys -> params receiver
!   (absent keys keep field defaults) -> the member's own init; an
!   unselected member skips.
!=====================================================================
module reg_samp_surface_oscillator
   use sampler,                 only: samp_reg, samp_code
   use config,                  only: reactants
   use input,                   only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,               only: seam_real, seam_int
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
   ! reg_samp_surface_oscillator_init() - assembly seam: buffered keys ->
   !                                       receiver -> member init
   subroutine reg_samp_surface_oscillator_init()
      type(surf_osc_params_t) :: p     ! member-parameter receiver (field defaults)
      character(len=512) :: sval       ! one buffered raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('surface_oscillator')) == 0) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('N_OSC');   call seam_int(buffer_key(i), sval, buffer_line(i), p%n_osc)
         case ('T_OSC');   call seam_real(buffer_key(i), sval, buffer_line(i), p%t_osc)
         case ('N_LEVEL'); call seam_int(buffer_key(i), sval, buffer_line(i), p%n_level)
         case ('K_OSC');   call seam_real(buffer_key(i), sval, buffer_line(i), p%k_osc)
         end select
      end do
      call surface_oscillator_init(p)
   end subroutine reg_samp_surface_oscillator_init
end module reg_samp_surface_oscillator
