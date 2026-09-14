!=====================================================================
! reg_dist_boltzmann.f90 - registration module (methods/): ONE samp_reg
!                          call registers the 'boltzmann' distribution
!                          member row - the reg file is the atomic
!                          registration unit (a new member = its method
!                          file + this one thin file; zero edits to
!                          existing files)
! Contract:
!   reg_dist_boltzmann_all() - register the boltzmann member's sampling
!     routine under its scheme word 'boltzmann' (the word the input
!     grammar accepts for DIST_SCHEME and the per-fragment dispatch
!     resolves). Passing boltz_sample as samp_reg's procedure dummy
!     checks its conformance to the dist_proc_i interface at compile
!     time. Executed once in the assembly phase, before any samp_run.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_boltzmann
   use sampler,        only: samp_reg
   use dist_boltzmann, only: boltz_sample
   implicit none
   private
   public :: reg_dist_boltzmann_all
contains
   !------------------------------------------------------------------
   ! reg_dist_boltzmann_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_boltzmann_all()
      call samp_reg('boltzmann', boltz_sample)
   end subroutine reg_dist_boltzmann_all
end module reg_dist_boltzmann
