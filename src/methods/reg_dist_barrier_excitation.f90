!=====================================================================
! reg_dist_barrier_excitation.f90 - registration module (methods/): ONE
!                         samp_reg call registers the
!                         'barrier_excitation' distribution member row -
!                         the reg file is the atomic registration unit
!                         (a new member = its method file + this one thin
!                         file; zero edits to existing files)
! Contract:
!   reg_dist_barrier_excitation_all() - register the barrier-excitation
!     member's sampling routine under its scheme word
!     'barrier_excitation' (the word the input grammar accepts for
!     DIST_SCHEME and the per-fragment dispatch resolves). Passing
!     barrier_excitation_sample as samp_reg's procedure dummy checks its
!     conformance to the dist_proc_i interface at compile time. Executed
!     once in the assembly phase, before any samp_run.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_barrier_excitation
   use sampler,                   only: samp_reg
   use dist_barrier_excitation,   only: barrier_excitation_sample
   implicit none
   private
   public :: reg_dist_barrier_excitation_all
contains
   !------------------------------------------------------------------
   ! reg_dist_barrier_excitation_all() - fill the member-table row (one
   ! call, one row)
   subroutine reg_dist_barrier_excitation_all()
      call samp_reg('barrier_excitation', barrier_excitation_sample)
   end subroutine reg_dist_barrier_excitation_all
end module reg_dist_barrier_excitation
