!=====================================================================
! reg_dist_thermalize.f90 - registration module (methods/): ONE samp_reg
!                         call registers the 'thermalize' distribution
!                         member row - the reg file is the atomic
!                         registration unit (a new member = its method
!                         file + this one thin file; zero edits to
!                         existing files)
! Contract:
!   reg_dist_thermalize_all() - register the MD-equilibration member's
!     sampling routine under its scheme word 'thermalize' (the word the
!     input grammar accepts for DIST_SCHEME and the per-fragment
!     dispatch resolves). Passing thermalize_sample as samp_reg's
!     procedure dummy checks its conformance to the dist_proc_i
!     interface at compile time. Executed once in the assembly phase,
!     before any samp_run.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_thermalize
   use sampler,        only: samp_reg
   use dist_thermalize, only: thermalize_sample
   implicit none
   private
   public :: reg_dist_thermalize_all
contains
   !------------------------------------------------------------------
   ! reg_dist_thermalize_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_thermalize_all()
      call samp_reg('thermalize', thermalize_sample)
   end subroutine reg_dist_thermalize_all
end module reg_dist_thermalize
