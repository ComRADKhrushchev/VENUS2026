!=====================================================================
! reg_dist_thermostat.f90 - register the 'thermostat' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: staged input keys -> params receiver
!   (absent keys keep field defaults) -> the member's own init; an
!   unselected member skips.
!=====================================================================
module reg_dist_thermostat
   use sampler,        only: samp_reg, samp_code
   use config,         only: reactants
   use input,          only: staging_n_rows, staging_key, staging_val, staging_line
   use seam_cast,      only: seam_real
   use dist_thermostat, only: thermostat_sample, thermostat_params_t, thermostat_init
   implicit none
   private
   public :: reg_dist_thermostat_all
   public :: reg_dist_thermostat_init
contains
   !------------------------------------------------------------------
   ! reg_dist_thermostat_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_thermostat_all()
      call samp_reg('thermostat', thermostat_sample)
   end subroutine reg_dist_thermostat_all

   !------------------------------------------------------------------
   ! reg_dist_thermostat_init() - assembly seam: staged keys -> receiver ->
   !                              member init
   subroutine reg_dist_thermostat_init()
      type(thermostat_params_t) :: p   ! member-parameter receiver (field defaults)
      character(len=512) :: sval       ! one staged raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('thermostat')) == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('T_BATH');  call seam_real(staging_key(i), sval, staging_line(i), p%t_bath)
         case ('NU_COLL'); call seam_real(staging_key(i), sval, staging_line(i), p%nu_coll)
         end select
      end do
      call thermostat_init(p)
   end subroutine reg_dist_thermostat_init
end module reg_dist_thermostat
