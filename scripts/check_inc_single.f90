!=====================================================================
! check_inc_single.f90 - the trivial single-reactant incident-member
!                          criteria program (criteria defined FIRST in
!                          this header; this program executes them)
! Compile (from src/, linked against the tree's module objects):
!   gfortran -Jbuild/gfortran -c ../scripts/check_inc_single.f90 -o build/ck_inc_single.o
!   gfortran build/ck_inc_single.o <tree objects> -o build/check_inc_single.exe
! All artifacts are written under build/ - the source tree stays clean.
! The member contract: a single-fragment gas-phase run has no collision
! partner - draw consumes NO RNG, realize touches NO state (both bitwise
! no-ops; the historical skip-guard produced the identical stream).
! Criteria:
!   (a) state immutability: a sentinel-filled state passes through
!       inc_single_sample BITWISE unchanged (q, p, f)
!   (b) stream immutability: the draw consumes no RNG - the next rng_u
!       after a draw equals the control stream's next value bitwise
!=====================================================================
program check_inc_single
   use rng,          only: rng_init, rng_u
   use state,        only: state_t, state_create, state_destroy
   use inc_single,   only: inc_single_draw, inc_single_realize, inc_single_sample
   implicit none
   type(state_t) :: sta
   real(8) :: u_after, u_ctrl
   real(8), parameter :: sentinel = -7.31d0
   integer :: nfail = 0

   ! ---------- (a) state immutability ----------
   call state_create(sta, 3)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call inc_single_sample(sta)
   if (any(sta%q /= sentinel) .or. any(sta%p /= sentinel) .or. any(sta%f /= sentinel)) then
      print '(a)', 'FAIL (a) state immutability: the member wrote a slot it does not own'
      nfail = nfail + 1
   else
      print '(a)', '(a) state immutability: PASS (sentinel q/p/f bitwise unchanged)'
   end if
   call state_destroy(sta)

   ! ---------- (b) stream immutability ----------
   call rng_init(4271)
   call inc_single_draw()
   u_after = rng_u()
   call rng_init(4271)
   u_ctrl = rng_u()
   if (u_after /= u_ctrl) then
      print '(a)', 'FAIL (b) stream immutability: the draw consumed RNG'
      nfail = nfail + 1
   else
      print '(a)', '(b) stream immutability: PASS (the next rng_u after a draw equals '// &
         'the control stream bitwise - no RNG consumed)'
   end if

   print '(a,i0,a)', '=== check_inc_single: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0
end program check_inc_single
