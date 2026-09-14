!=====================================================================
! reg_prop.f90 - thin closed-family registration module (methods/): TWO
!                prop_reg calls register the verlet and symple rows of the
!                integrator closed family - the reg file is the atomic
!                registration unit (a new member = its method file + this
!                one thin file; zero edits to existing files), aggregated
!                into the generated registry module by the prebuild
!                wildcard
! Contract:
!   reg_prop_all() - register the verlet member under integrator id 1
!     (velocity Verlet at order=1 / third-order Beeman at order=2) and
!     the symple member under id 2 (the 4/6/8th-order drift-kick
!     composition) - order is the per-row runtime switch of
!     config%propagator%order. Passing each *_step as prop_reg's
!     procedure dummy checks its conformance to the integ_proc_i
!     interface at compile time. The radau row lands with its own
!     member's concretization task - one row per member, each registered
!     exactly once at assembly.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_prop
   use propagator, only: prop_reg, verlet_step, symple_step
   implicit none
   private
   public :: reg_prop_all
contains
   !------------------------------------------------------------------
   ! reg_prop_all() - fill the verlet and symple rows of the closed-family registry
   subroutine reg_prop_all()
      call prop_reg(1, verlet_step)
      call prop_reg(2, symple_step)
   end subroutine reg_prop_all
end module reg_prop
