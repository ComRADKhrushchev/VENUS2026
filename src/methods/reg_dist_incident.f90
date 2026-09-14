!=====================================================================
! reg_dist_incident.f90 - registration module (methods/): ONE samp_reg
!                         call registers the 'incident' distribution
!                         member row - the reg file is the atomic
!                         registration unit (a new member = its method
!                         file + this one thin file; zero edits to
!                         existing files)
! Contract:
!   reg_dist_incident_all() - register the incident member's sampling
!     routine under the name 'incident'. The incident channel is the
!     one member NOT selected by a per-fragment DIST_SCHEME row (the
!     word 'incident' is withdrawn from that grammar - an illegal value
!     there): the driver passes this member's name as samp_run's scheme
!     argument, and this registration row is what that lookup resolves.
!     Passing inc_sample as samp_reg's procedure dummy checks its
!     conformance to the dist_proc_i interface at compile time.
!     Executed once in the assembly phase, before any samp_run.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_incident
   use sampler,       only: samp_reg
   use dist_incident, only: inc_sample
   implicit none
   private
   public :: reg_dist_incident_all
contains
   !------------------------------------------------------------------
   ! reg_dist_incident_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_incident_all()
      call samp_reg('incident', inc_sample)
   end subroutine reg_dist_incident_all
end module reg_dist_incident
