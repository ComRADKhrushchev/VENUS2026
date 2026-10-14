!=====================================================================
! io.f90 - the auto_sample fixture's input/output group: load/classify/
!   observables reserved names (the input/final_state/recorder layer)
! Design:
!   load pulls K_SPRING/R_TERM buffered keys and hands them to pot's
!   module state via auto_sample_set_params (absent key keeps the
!   default; found carries the R_TERM optional-slot semantics).
!   classify is the wiring-exercise stub: constant 'bound' (fin unread)
!   - the drill asserts it end-to-end through the real binary.
!   observables is the 11th reserved name's exercise (2026-10-14):
!   one stub recorder column, constant 1.0.
! Units: none (judgments carry code references only)
!=====================================================================
module auto_sample_io
   use state,      only: state_t
   use input,      only: pull
   use final_state, only: fin_t
   use recorder,   only: rec_reg_col
   use auto_sample_pot, only: auto_sample_set_params
   implicit none
   private
   public :: load, classify, observables
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

   !------------------------------------------------------------------
   ! observables() - the 11th reserved name's exercise: one stub column
   !------------------------------------------------------------------
   subroutine observables()
      call rec_reg_col('stub_count', stub_count)
   end subroutine observables

   real(8) function stub_count(sta)
      type(state_t), intent(in) :: sta
      stub_count = 1.0d0
   end function stub_count
end module auto_sample_io
