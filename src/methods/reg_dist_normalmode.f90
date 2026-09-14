!=====================================================================
! reg_dist_normalmode.f90 - registration module (methods/): ONE samp_reg
!                         call registers the 'normalmode' distribution
!                         member row - the reg file is the atomic
!                         registration unit (a new member = its method
!                         file + this one thin file; zero edits to
!                         existing files)
! Contract:
!   reg_dist_normalmode_all() - register the normal-mode member's
!     sampling routine under its scheme word 'normalmode' (the word the
!     input grammar accepts for DIST_SCHEME and the per-fragment
!     dispatch resolves). Passing normalmode_sample as samp_reg's
!     procedure dummy checks its conformance to the dist_proc_i
!     interface at compile time. Executed once in the assembly phase,
!     before any samp_run.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_normalmode
   use sampler,          only: samp_reg
   use dist_normalmode,  only: normalmode_sample
   implicit none
   private
   public :: reg_dist_normalmode_all
contains
   !------------------------------------------------------------------
   ! reg_dist_normalmode_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_normalmode_all()
      call samp_reg('normalmode', normalmode_sample)
   end subroutine reg_dist_normalmode_all
end module reg_dist_normalmode
