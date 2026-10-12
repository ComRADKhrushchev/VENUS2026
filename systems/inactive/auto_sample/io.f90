!=====================================================================
! io.f90 - the auto_sample fixture's input/output group: load/classify
!   magic names (the input/final_state dependency layer)
! Design:
!   load pulls K_SPRING/R_TERM buffered keys and hands them to pot's
!   module state via auto_sample_set_params (absent key keeps the
!   default; found carries the R_TERM optional-slot semantics).
!   classify is the wiring-exercise stub: constant 'bound' (fin unread)
!   - the drill asserts it end-to-end through the real binary.
! Units: none (judgments carry code references only)
!=====================================================================
module auto_sample_io
   use input,       only: pull
   use final_state, only: fin_t
   use auto_sample_pot, only: auto_sample_set_params
   implicit none
   private
   public :: load, classify
contains
   !------------------------------------------------------------------
   ! load() - pull the buffered container keys, stage pot's state
   !------------------------------------------------------------------
   subroutine load()
      real(8) :: k, r
      logical :: found
      k = 5.0d0                          ! field defaults ARE the container defaults
      r = 0.0d0
      found = .false.
      call pull('K_SPRING', k)
      call pull('R_TERM', r, found)
      call auto_sample_set_params(k, r, found)
   end subroutine load

   !------------------------------------------------------------------
   ! classify(fin) - wiring-exercise stub: the bound channel (fin unread)
   !------------------------------------------------------------------
   function classify(fin) result(name)
      type(fin_t), intent(in) :: fin     ! final-state quantities (unused)
      character(len=32) :: name
      name = 'bound'
   end function classify
end module auto_sample_io
