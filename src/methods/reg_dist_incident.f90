!=====================================================================
! reg_dist_incident.f90 - registration module (methods/): ONE samp_reg
!                          call registers the 'incident' distribution
!                          member row - the reg file is the atomic
!                          registration unit (a new member = its method
!                          file + this one thin file; zero edits to
!                          existing files)
! Contract:
!   reg_dist_incident_all() - register the incident member's sampling
!     routine under its scheme word 'incident' (the word the per-
!     fragment dispatch never resolves - the incident member runs once
!     after the per-fragment loop, selected by the dispatch rule, not by
!     a dist_scheme row). Passing inc_sample as samp_reg's procedure
!     dummy checks its conformance to the dist_proc_i interface at
!     compile time. Executed once in the assembly phase, before any
!     samp_run.
!   reg_dist_incident_init() - the ASSEMBLY SEAM entry: load this
!     member's staged input parameters (the input framework's public
!     staging rows) into the member's params receiver - absent keys keep
!     the field defaults - then run the member's own init (validation).
!     Guard mirrors the dispatch rule (not the scheme table): a run with
!     a single fragment and no surface model never dispatches incident,
!     so the member skips (its init would otherwise demand collision
!     parameters on every such run). Three-party zero knowledge: the
!     driver knows no member, the input knows no member, the member
!     never touches the driver.
! Units: none (registration carries names and code references only).
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
   ! reg_dist_incident_init() - the assembly seam: staged member keys ->
   !                             params receiver -> the member's own init
   subroutine reg_dist_incident_init()
      type(inc_params_t) :: p           ! member-parameter receiver (field defaults)
      character(len=512) :: sval        ! one staged raw value
      integer :: i
      integer :: nf                     ! roster fragment count (the dispatch-rule mirror)
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
