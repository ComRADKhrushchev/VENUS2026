!=====================================================================
! reg_samp_barrier_excitation.f90 - register the 'barrier_excitation' member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: buffered input keys -> params receiver
!   (absent keys keep defaults) -> the member's own init; an unselected
!   member skips. The w_mode/c_mode tables are assembly data, not input
!   keys - a selected member aborts with a named error at its own unallocated check.
!=====================================================================
module reg_samp_barrier_excitation
   use sampler,                only: samp_reg, samp_code
   use config,                 only: reactants
   use input,                  only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,              only: seam_real, seam_int
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
   ! reg_samp_barrier_excitation_init() - assembly seam: buffered keys -> receiver ->
   !                                       member init
   subroutine reg_samp_barrier_excitation_init()
      type(barrier_excitation_params_t) :: p   ! member-parameter receiver (field defaults)
      character(len=512) :: sval               ! one buffered raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('barrier_excitation')) == 0) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('E_STAB');  call seam_real(buffer_key(i), sval, buffer_line(i), p%e_stab)
         case ('N_E_BAR'); call seam_int(buffer_key(i), sval, buffer_line(i), p%n_e_bar)
         case ('E_BAR');   call seam_real(buffer_key(i), sval, buffer_line(i), p%e_bar)
         case ('T_BAR');   call seam_real(buffer_key(i), sval, buffer_line(i), p%t_bar)
         end select
      end do
      call barrier_excitation_init(p)
   end subroutine reg_samp_barrier_excitation_init
end module reg_samp_barrier_excitation
