!=====================================================================
! reg_dist_incident.f90 - register the 'incident' distribution member
! Design:
!   One samp_reg call fills the row; incident runs once after the
!   per-fragment loop (dispatch rule, not a dist_scheme row); the reg
!   file is the atomic registration unit (new member = method file +
!   one thin reg file; zero edits to existing files).
!   _init is the assembly seam: staged input keys -> params receiver
!   (absent keys keep defaults) -> the member's own init; the guard
!   mirrors the dispatch rule - a single-fragment run with no surface
!   model never dispatches incident, so the member skips.
!=====================================================================
module reg_dist_incident
   use sampler,        only: samp_reg
   use config,         only: reactants
   use config_roster,  only: roster
   use input,          only: staging_n_rows, staging_key, staging_val, staging_line
   use seam_cast,      only: seam_real, seam_int
   use dist_incident,  only: inc_sample, inc_params_t, inc_init
   implicit none
   private
   public :: reg_dist_incident_all
   public :: reg_dist_incident_init
contains
   !------------------------------------------------------------------
   ! reg_dist_incident_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_incident_all()
      call samp_reg('incident', inc_sample)
   end subroutine reg_dist_incident_all

   !------------------------------------------------------------------
   ! reg_dist_incident_init() - assembly seam: staged keys -> receiver ->
   !                             member init
   subroutine reg_dist_incident_init()
      type(inc_params_t) :: p           ! member-parameter receiver (field defaults)
      character(len=512) :: sval        ! one staged raw value
      integer :: i
      integer :: nf                     ! roster fragment count (the dispatch-rule mirror)
      ! 1. dispatch-rule guard: single fragment + no surface model -> member skips
      if (.not. allocated(roster%frag)) return
      nf = size(roster%frag)
      if (nf <= 1 .and. reactants%surface_model == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('N_E_REL');   call seam_int(staging_key(i), sval, staging_line(i), p%n_e_rel)
         case ('T_TRANS');   call seam_real(staging_key(i), sval, staging_line(i), p%t_trans)
         case ('V_WIDTH');   call seam_real(staging_key(i), sval, staging_line(i), p%v_width)
         case ('N_THTA');    call seam_int(staging_key(i), sval, staging_line(i), p%n_thta)
         case ('THTA_MAX');  call seam_real(staging_key(i), sval, staging_line(i), p%thta_max)
         case ('N_CHI');     call seam_int(staging_key(i), sval, staging_line(i), p%n_chi)
         case ('CHI');       call seam_real(staging_key(i), sval, staging_line(i), p%chi)
         case ('N_B');       call seam_int(staging_key(i), sval, staging_line(i), p%n_b)
         case ('N_AIM');     call seam_int(staging_key(i), sval, staging_line(i), p%n_aim)
         end select
      end do
      call inc_init(p)
   end subroutine reg_dist_incident_init
end module reg_dist_incident
