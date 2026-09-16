!=====================================================================
! reg_dist_boltzmann.f90 - register the 'boltzmann' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: staged input keys -> params receiver
!   (absent keys keep defaults) -> the member's own init; an unselected
!   member skips. The w_mode_a table is assembly data, not an input
!   key - a selected member loud-stops at its own unallocated check.
!=====================================================================
module reg_dist_boltzmann
   use sampler,        only: samp_reg, samp_code
   use config,         only: reactants
   use input,          only: staging_n_rows, staging_key, staging_val, staging_line
   use seam_cast,      only: seam_real
   use dist_boltzmann, only: boltz_sample, boltz_params_t, boltz_init
   implicit none
   private
   public :: reg_dist_boltzmann_all
   public :: reg_dist_boltzmann_init
contains
   !------------------------------------------------------------------
   ! reg_dist_boltzmann_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_boltzmann_all()
      call samp_reg('boltzmann', boltz_sample)
   end subroutine reg_dist_boltzmann_all

   !------------------------------------------------------------------
   ! reg_dist_boltzmann_init() - assembly seam: staged keys -> receiver ->
   !                              member init
   subroutine reg_dist_boltzmann_init()
      type(boltz_params_t) :: p          ! member-parameter receiver (field defaults)
      character(len=512) :: sval         ! one staged raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('boltzmann')) == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('T_VIB_A'); call seam_real(staging_key(i), sval, staging_line(i), p%t_vib_a)
         case ('T_VIB_B'); call seam_real(staging_key(i), sval, staging_line(i), p%t_vib_b)
         end select
      end do
      call boltz_init(p)
   end subroutine reg_dist_boltzmann_init
end module reg_dist_boltzmann
