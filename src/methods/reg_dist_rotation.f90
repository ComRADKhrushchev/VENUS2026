!=====================================================================
! reg_dist_rotation.f90 - registration module (methods/): ONE samp_reg
!                         call registers the 'rotation' distribution
!                         member row - the reg file is the atomic
!                         registration unit (a new member = its method
!                         file + this one thin file; zero edits to
!                         existing files)
! Contract:
!   reg_dist_rotation_all() - register the rotation member's sampling
!     routine under its scheme word 'rotation' (the word the input
!     grammar accepts for DIST_SCHEME and the per-fragment dispatch
!     resolves). Passing rotation_sample as samp_reg's procedure dummy
!     checks its conformance to the dist_proc_i interface at compile
!     time. Executed once in the assembly phase, before any samp_run.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_rotation
   use sampler,       only: samp_reg
   use dist_rotation, only: rotation_sample
   implicit none
   private
   public :: reg_dist_rotation_all
contains
   !------------------------------------------------------------------
   ! reg_dist_rotation_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_rotation_all()
      call samp_reg('rotation', rotation_sample)
   end subroutine reg_dist_rotation_all
end module reg_dist_rotation
