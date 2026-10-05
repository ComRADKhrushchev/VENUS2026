!------------------------------------------------------------------
! PES_h2ag.f90 - the user-PES seam: H2/Ag(111) Jiang&Guo PIP-NN surface.
! Engine: h2ag_engine.f - the distributed h2ag111.f, entries renamed
! H2AG_POT/H2AG_GRAD/H2AG_INIT, statement text untouched (bit-audited
! against the archived original by check_h2ag111).
! Seam contract (stock VERLET header): POT0(N,V) reads Q [A] from
! COMMON/QPDOT and returns V [kcal/mol*C1]; DPESHON(N) fills PDOT with
! the physical force F = -dV/dq [internal/A]. The engine speaks A/eV and
! the central-difference gradient, so the only work here is the unit
! fold  eV -> kcal/mol (23.0605) -> internal (C1) = 0.96485132 per eV
! and the force sign flip. NATOMS must be 2 (the Ag(111) slab is rigid
! inside the PES; the 4 ghost cell-corner atoms live beyond NATOMS and
! are never touched by the force engine).
!------------------------------------------------------------------
subroutine POTPRE(nPES)
   implicit none
   integer :: nPES
   if (nPES /= 1) then
      write (0, '(a,i0)') 'POTPRE: the h2ag seam carries exactly one PES, got nPES = ', nPES
      stop 1
   end if
   call H2AG_INIT()             ! one-time: reads weights/biases-h2ag.txt from cwd
end subroutine POTPRE

subroutine POT0(NATOMS, V)
   implicit double precision (a-h, o-z)
   parameter (NDA = 160, NDA3 = NDA*3)
   double precision :: kcal_ev
   dimension q6(6)
   common/QPDOT/Q(NDA3), PDOT(NDA3), FCOEF(NDA3, NDA3)
   common/CONSTN/C1, C2, C3, C4, C5, C6, C7, PI, HALFPI, TWOPI
   if (NATOMS /= 2) then
      write (0, '(a,i0)') 'POT0: the h2ag seam serves exactly 2 atoms, got NATOMS = ', NATOMS
      stop 1
   end if
   kcal_ev = 23.0605d0          ! 1 eV [kcal/mol]
   q6(1:6) = Q(1:6)
   call H2AG_POT(2, q6, V)
   V = V*kcal_ev*C1
end subroutine POT0

subroutine DPESHON(NATOMS)
   implicit double precision (a-h, o-z)
   parameter (NDA = 160, NDA3 = NDA*3)
   double precision :: kcal_ev
   dimension q6(6), g6(6)
   integer :: k
   common/QPDOT/Q(NDA3), PDOT(NDA3), FCOEF(NDA3, NDA3)
   common/CONSTN/C1, C2, C3, C4, C5, C6, C7, PI, HALFPI, TWOPI
   if (NATOMS /= 2) then
      write (0, '(a,i0)') 'DPESHON: the h2ag seam serves exactly 2 atoms, got NATOMS = ', NATOMS
      stop 1
   end if
   kcal_ev = 23.0605d0
   q6(1:6) = Q(1:6)
   call H2AG_GRAD(2, q6, g6)
   do k = 1, 6
      PDOT(k) = -g6(k)*kcal_ev*C1   ! the machinery advances q with the force
   end do
end subroutine DPESHON
