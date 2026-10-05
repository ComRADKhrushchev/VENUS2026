!=====================================================================
! reg_prop.f90 - register the verlet, symple and radau integrator rows
! Design:
!   Three prop_reg calls fill the closed integrator family: verlet
!   under id 1, symple under id 2, radau under id 3; order is the
!   per-row runtime switch of config%propagator%order (for radau,
!   ss = 10^-order). One row per member, registered exactly once at
!   assembly; the reg file is the atomic registration unit (a new
!   member = its method file + one thin reg file; zero edits to
!   existing files).
!=====================================================================
module reg_prop
   use propagator, only: prop_reg, verlet_step, symple_step, radau_step
   implicit none
   private
   public :: reg_prop_all
contains
   !------------------------------------------------------------------
   ! reg_prop_all() - fill the verlet (1), symple (2) and radau (3) rows
   subroutine reg_prop_all()
      call prop_reg(1, verlet_step)
      call prop_reg(2, symple_step)
      call prop_reg(3, radau_step)
   end subroutine reg_prop_all
end module reg_prop
