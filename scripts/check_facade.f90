!=====================================================================
! check_facade.f90 - facade-mechanism criteria program, first batch:
!                    the frozen ws holder association (criteria defined
!                    FIRST in this header; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/config.f90 globals/state.f90 facade/elec_ws.f90 &
!            facade/elec_facade.f90 ../scripts/check_facade.f90 -o build/check_facade.exe
! All artifacts are written under build/ - the source tree stays clean.
! Posture: the subject is the holder mechanism itself, so no force chain
!   is assembled - a degenerate no-op method row registers + binds the
!   name, which seeds the statistics holder and makes the occ_get channel
!   live; the members are never dispatched by this batch. No loud-stop
!   criteria here, hence no child-process machinery.
! Criteria:
!   (a) handle identity: two ws_handle() calls yield pointers to the SAME
!       target (pointer association comparison), and a field written
!       through one handle is visible through the other handle AND through
!       the ws_get snapshot channel (the alias is onto the holder itself,
!       not onto a copy)
!   (b) copy-out semantics: ws_get returns a snapshot - scribbling the
!       snapshot (a scalar component and an array component) leaves the
!       holder intact, and a fresh ws_get reproduces the holder's
!       pre-scribble content
!   (c) reset semantics: with the holder dirtied (repr set, scalar set,
!       every per-mode array allocated and filled), ws_reset returns it to
!       the DEFINED empty state - a subsequent ws_get shows repr = r_empty,
!       every per-mode array unallocated, the scalar zeroed (the one-step
!       survival contract: representation constructs live one integration
!       step, then die - not merely become unreadable)
!   (d) snapshot-style vs handle-style coexistence: the existing consumer
!       pattern of the scalar force branch (occ_get request -> ws_get
!       snapshot -> compose into the snapshot) runs unchanged and leaves
!       the HOLDER untouched (zero interference), while the handle-style
!       compose (the fill-site posture that writes the holder itself)
!       reads back through the snapshot channel without retroactively
!       perturbing an earlier snapshot. The regression half of this
!       criterion is external: check_propagator (the force_eval scalar
!       branch end-to-end on real forces) and check_driver (the full
!       four-phase loop) re-run green on the rewritten holder.
!=====================================================================
program check_facade
   use state,       only: state_t
   use elec_ws,     only: ws_t, r_empty, r_scalar, r_diab
   use elec_facade, only: ws_handle, ws_get, ws_reset, occ_get, &
                          method_reg, method_bind
   implicit none
   type(ws_t) :: w, w2                     ! snapshot instances (ws_get dummies)
   type(ws_t), pointer :: p1, p2           ! handle aliases
   integer, allocatable :: occ_req(:)      ! the delivery request (occ_get)
   integer :: nfail = 0

   ! degenerate no-op method row (the sibling checks' posture): registers +
   ! binds the name and seeds the statistics holder - occ_get is live for
   ! criterion (d); the members are never dispatched by this batch
   call method_reg('adiabatic', '', noop_prop, noop_match)
   call method_bind('adiabatic')

   ! ---------- (a) handle identity + write-through observability ----------
   call ws_reset()
   p1 => ws_handle()
   p2 => ws_handle()
   allocate (p1%g_add(3))
   p1%g_add = [1.0d0, 2.0d0, 3.0d0]
   p1%v_add = 3.25d0
   call ws_get(w)
   if (.not. (associated(p1) .and. associated(p2) .and. associated(p1, p2)) .or. &
       p2%v_add /= 3.25d0 .or. any(p2%g_add /= [1.0d0, 2.0d0, 3.0d0]) .or. &
       w%v_add /= 3.25d0 .or. any(w%g_add /= [1.0d0, 2.0d0, 3.0d0])) then
      print '(a)', 'FAIL (a) handle identity: the two handles do not alias one '// &
                   'target, or a handle write is invisible to the other handle / ws_get'
      nfail = nfail + 1
   else
      print '(a)', '(a) handle identity: PASS (two ws_handle() calls alias one '// &
                   'target; the write through one is seen by the other and by ws_get)'
   end if

   ! ---------- (b) copy-out semantics: the snapshot is not the holder ----------
   w%v_add = -999.0d0                      ! scribble on the snapshot taken above
   w%g_add = -999.0d0                      ! (a scalar and an array component)
   call ws_get(w2)                         ! a fresh snapshot of the holder
   if (w2%v_add /= 3.25d0 .or. any(w2%g_add /= [1.0d0, 2.0d0, 3.0d0])) then
      print '(a)', 'FAIL (b) copy-out: writing the snapshot perturbed the holder'
      nfail = nfail + 1
   else
      print '(a)', '(b) copy-out: PASS (scribbling the snapshot leaves the holder '// &
                   'intact; a fresh ws_get reproduces the held content)'
   end if

   ! ---------- (c) reset semantics: full clear to the defined empty state ----------
   call ws_reset()                         ! clean slate so every allocation below is fresh
   p1 => ws_handle()
   p1%repr = r_diab
   p1%v_add = -777.0d0
   allocate (p1%i_want(2))
   p1%i_want = [1, 2]
   allocate (p1%g_add(3))
   p1%g_add = 1.0d0
   allocate (p1%h_ij(2, 2))
   p1%h_ij = 1.0d0
   allocate (p1%dh_ij(3, 2, 2))
   p1%dh_ij = 1.0d0
   allocate (p1%v_surf(2))
   p1%v_surf = 1.0d0
   allocate (p1%g_surf(3, 2))
   p1%g_surf = 1.0d0
   allocate (p1%d_ij(3, 2))
   p1%d_ij = 1.0d0
   call ws_reset()
   call ws_get(w)
   if (w%repr /= r_empty .or. w%v_add /= 0.0d0 .or. allocated(w%i_want) .or. &
       allocated(w%g_add) .or. allocated(w%h_ij) .or. allocated(w%dh_ij) .or. &
       allocated(w%v_surf) .or. allocated(w%g_surf) .or. allocated(w%d_ij)) then
      print '(a)', 'FAIL (c) reset: the holder did not return to the defined empty '// &
                   'state (repr, scalar, or a per-mode array survived)'
      nfail = nfail + 1
   else
      print '(a)', '(c) reset: PASS (repr = r_empty, all seven per-mode arrays '// &
                   'gone, scalar zeroed - the constructs die at step end)'
   end if

   ! ---------- (d) snapshot-style vs handle-style coexistence ----------
   call ws_reset()
   call occ_get(occ_req)
   ! the snapshot fill pattern of the scalar force branch, verbatim: request
   ! through occ_get, snapshot through ws_get, compose into the snapshot
   call ws_get(w)
   if (allocated(w%i_want)) deallocate (w%i_want)
   allocate (w%i_want(size(occ_req)))
   w%i_want = occ_req
   allocate (w%g_add(3))
   w%g_add = 0.5d0
   w%v_add = -1.5d0
   w%repr = r_scalar
   ! zero interference: the snapshot fill must leave the holder in its empty
   ! state (the pre-freeze consumer behavior, unchanged by the rewrite)
   p1 => ws_handle()
   if (size(occ_req) /= 1 .or. occ_req(1) /= 1 .or. &
       p1%repr /= r_empty .or. p1%v_add /= 0.0d0 .or. allocated(p1%i_want) .or. &
       allocated(p1%g_add)) then
      print '(a)', 'FAIL (d) coexistence: the occ channel broke, or the snapshot '// &
                   'fill trespassed on the holder'
      nfail = nfail + 1
   else
      print '(a)', '(d) coexistence part 1: PASS (occ request seeds [1]; the '// &
                   'snapshot fill leaves the holder empty - zero interference)'
   end if
   ! the handle-style compose (the fill-site posture onto the frozen holder):
   ! writes the holder itself, then reads back through the snapshot channel
   p1%repr = r_scalar
   allocate (p1%i_want(size(occ_req)))
   p1%i_want = occ_req
   allocate (p1%g_add(3))
   p1%g_add = 0.25d0
   p1%v_add = -2.5d0
   call ws_get(w2)
   if (p1%repr /= r_scalar .or. p1%v_add /= -2.5d0 .or. any(p1%g_add /= 0.25d0) .or. &
       w2%repr /= r_scalar .or. w2%v_add /= -2.5d0 .or. any(w2%g_add /= 0.25d0) .or. &
       any(w2%i_want /= 1) .or. &
       w%v_add /= -1.5d0 .or. any(w%g_add /= 0.5d0)) then
      print '(a)', 'FAIL (d) coexistence: the handle read-back broke, or the earlier '// &
                   'snapshot was retroactively perturbed'
      nfail = nfail + 1
   else
      print '(a)', '(d) coexistence part 2: PASS (the handle write reads back '// &
                   'through ws_get; the earlier snapshot keeps its own content)'
   end if
   call ws_reset()                         ! leave the holder in the empty state

   print '(a,i0,a)', '=== check_facade: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! degenerate no-op method members (the all-noop package this check registers
   ! so the method-name lookup resolves and the statistics holder is seeded)
   subroutine noop_prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine noop_prop

   subroutine noop_match(sta)
      type(state_t), intent(inout) :: sta
   end subroutine noop_match

end program check_facade
