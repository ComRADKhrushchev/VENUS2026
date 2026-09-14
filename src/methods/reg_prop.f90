!=====================================================================
! reg_prop.f90 - thin closed-family registration module (methods/): ONE
!                prop_reg call registers the verlet row of the integrator
!                closed family - the reg file is the atomic registration
!                unit (a new member = its method file + this one thin
!                file; zero edits to existing files), aggregated into the
!                generated registry module by the prebuild wildcard
! Contract:
!   reg_prop_all() - register the verlet member under integrator id 1
!     (the key of config%propagator%id): velocity Verlet at order=1 /
!     third-order Beeman at order=2, one row - order is the per-row
!     runtime switch of config%propagator%order. Passing verlet_step as
!     prop_reg's procedure dummy checks its conformance to the
!     integ_proc_i interface at compile time. The symple/radau rows land
!     with their own members' concretization tasks - one row per task,
!     each registered exactly once at assembly.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_prop
   use propagator, only: prop_reg, verlet_step
   implicit none
   private
   public :: reg_prop_all
contains
   !------------------------------------------------------------------
   ! reg_prop_all() - fill the verlet row of the closed-family registry (one call, one row)
   subroutine reg_prop_all()
      call prop_reg(1, verlet_step)
   end subroutine reg_prop_all
end module reg_prop
