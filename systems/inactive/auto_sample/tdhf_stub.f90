!=====================================================================
! tdhf_stub.f90 - the auto_sample fixture's method group: a no-op
!   package named auto_tdhf (prop/match/sync/init/method_name)
! Design:
!   The no-op propagation and matching forms are the legal degenerate
!   members (elec_prop_i lists no-op; match without a hop writes
!   nothing - occ stays the seeded [1]). init writes the witness file
!   auto_init_witness.txt into the run directory - the drill's
!   selection-gate evidence (present iff the run ran the package init).
! Units: none
!=====================================================================
module auto_sample_met
   use state, only: state_t, elec_t
   implicit none
   private
   public :: prop, match, sync, init, method_name
   logical :: inited = .false.            ! witness latch [-]
contains
   !------------------------------------------------------------------
   ! method_name() - the registration name (lowercase, author data)
   !------------------------------------------------------------------
   function method_name() result(s)
      character(len=32) :: s
      s = 'auto_tdhf'
   end function method_name

   !------------------------------------------------------------------
   ! prop(sta, dt) - no-op propagation (legal degenerate form)
   !------------------------------------------------------------------
   subroutine prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine prop

   !------------------------------------------------------------------
   ! match(sta) - no-op matching (no hop, occ stays the seeded [1])
   !------------------------------------------------------------------
   subroutine match(sta)
      type(state_t), intent(inout) :: sta
   end subroutine match

   !------------------------------------------------------------------
   ! sync(s) - no-op container sync (no deliverable cache)
   !------------------------------------------------------------------
   subroutine sync(s)
      type(elec_t), intent(inout) :: s
   end subroutine sync

   !------------------------------------------------------------------
   ! init() - the selected-package witness (writes the drill's file once)
   !------------------------------------------------------------------
   subroutine init()
      integer :: u                        ! witness unit (newunit, the repo convention)
      if (inited) return
      inited = .true.
      open (newunit=u, file='auto_init_witness.txt', status='replace', action='write')
      write (u, '(a)') 'auto_tdhf init ran'
      close (u)
   end subroutine init
end module auto_sample_met
