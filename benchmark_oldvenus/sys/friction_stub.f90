!------------------------------------------------------------------
! friction_stub.f90 - loud-stop replacements for the GLO ghost-atom
! friction entries (FRICTION / FRICFORCE). The clean benchmark tree
! serves gas-phase 3-atom calibrations: NGLO must be 0 and these
! entries are never reached; if a run ever strays into the branch
! it stops with a named message instead of computing silently.
!------------------------------------------------------------------
subroutine FRICTION(NFC)
   implicit none
   integer :: NFC
   write (0, '(a)') 'FRICTION: the GLO/ghost-atom branch is not part of the clean benchmark tree'
   write (0, '(a)') 'FRICTION: fatal (run with NGLO = 0)'
   stop 1
end subroutine FRICTION

subroutine FRICFORCE(NFC, TELEC)
   implicit double precision (a-h, o-z)
   integer :: NFC
   write (0, '(a)') 'FRICFORCE: the GLO/ghost-atom branch is not part of the clean benchmark tree'
   write (0, '(a)') 'FRICFORCE: fatal (run with NGLO = 0)'
   stop 1
end subroutine FRICFORCE

!------------------------------------------------------------------
! readff / swap - loud no-ops for the molecular-mechanics bookkeeping
! the variant added for CH4+O (AmberFF.f90). The clean tree carries
! no force field: the PES seam owns the whole potential, so the
! parameter read is unnecessary and SWAP's atom exchange (gated on
! an MM parameter that is never set here) stays inert.
!------------------------------------------------------------------
subroutine readff
   implicit none
end subroutine readff

subroutine swap
   implicit none
end subroutine swap

!------------------------------------------------------------------
! TSH_Algorithm - no-op stub for the surface-hopping hook the
! variant's VERLET calls each step (the original gates itself on
! the TSH flags; the clean tree runs single-surface NSURF=0 and
! the real algorithm is retired with src_NN)
!------------------------------------------------------------------
subroutine TSH_Algorithm
   implicit none
end subroutine TSH_Algorithm
