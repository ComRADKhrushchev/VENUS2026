!=====================================================================
! pot.f90 - the auto_sample fixture's PES group: origin-spring well
!   (Tier 1 pure-PES container example; pes/keys/term magic names)
! Design:
!   V = 0.5*k*sum(q_i^2), g = k*q - every atom springs to the origin,
!   geometry-agnostic (any nat), closed-form checkable. k stages from
!   K_SPRING (default 5.0 eV/A^2) via auto_sample_set_params (io.f90's
!   load hands it over); r_term/r_found stage from R_TERM (absent key =
!   the MAX_STEPS-only run, term stays false).
! Units: q[A], g[eV/A], v[eV]
!=====================================================================
module auto_sample_pot
   use state, only: state_t
   implicit none
   private
   public :: pes, keys, term
   public :: auto_sample_set_params   ! fixture wiring: io.f90 load's handover
   real(8) :: k_spring = 5.0d0        ! spring constant [eV/A^2] (K_SPRING)
   real(8) :: r_term = 0.0d0          ! termination radius [A] (R_TERM)
   logical :: r_term_found = .false.  ! R_TERM presence flag [-]
contains
   !------------------------------------------------------------------
   ! keys() - the container key vocabulary (pure data, pre-parse)
   !------------------------------------------------------------------
   function keys() result(kk)
      character(len=16) :: kk(2)
      kk = [character(len=16) :: 'K_SPRING', 'R_TERM']
   end function keys

   !------------------------------------------------------------------
   ! pes(q, mass, g, v) - the force entry (the container force delivery, flat dof)
   !------------------------------------------------------------------
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)      ! coordinates [A] (3*nat, flattened)
      real(8), intent(in)  :: mass(:)   ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)      ! dV/dq [eV/A]
      real(8), intent(out) :: v         ! V [eV]
      v = 0.5d0*k_spring*dot_product(q, q)
      g = k_spring*q
   end subroutine pes

   !------------------------------------------------------------------
   ! term(sta) - termination: any coordinate leaves the R_TERM ball
   !             (key absent -> always false = MAX_STEPS-only)
   !------------------------------------------------------------------
   logical function term(sta)
      type(state_t), intent(in) :: sta
      term = .false.
      if (r_term_found) term = any(abs(sta%q) > r_term)
   end function term

   !------------------------------------------------------------------
   ! auto_sample_set_params(k, r, found) - the load handover (io.f90)
   !------------------------------------------------------------------
   subroutine auto_sample_set_params(k, r, found)
      real(8), intent(in) :: k, r
      logical, intent(in) :: found
      k_spring = k
      r_term = r
      r_term_found = found
   end subroutine auto_sample_set_params
end module auto_sample_pot
