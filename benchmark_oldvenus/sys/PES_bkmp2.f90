!------------------------------------------------------------------
! PES_bkmp2.f90 - the user-PES seam of the clean benchmark tree.
! The three entries the VENUS machinery calls:
!   POTPRE(nPES)  one-time init at startup (VENUS.f); single surface.
!   POT0(N,V)     energy:  reads Q [A] from COMMON/QPDOT, returns V
!                 in internal units (kcal/mol * C1).
!   DPESHON(N)    forces:  fills PDOT = -dV/dq [internal/A].
! Sign contract (stock VERLET.f header): F(T) = -dV/dq = PDOT -
! the integrators advance q WITH PDOT, so PDOT carries the physical
! force, NOT the gradient.
! Occupant: the official BKMP2 H3 surface (surface950621; the
! verbatim engine sits beside this file as bkmp2_engine.f; distances
! in bohr, energies in hartree at the engine boundary).
!------------------------------------------------------------------
subroutine POTPRE(nPES)
   implicit none
   integer :: nPES
   if (nPES /= 1) then
      write (0, '(a,i0)') 'POTPRE: the clean tree carries exactly one PES, got nPES = ', nPES
      stop 1
   end if
end subroutine POTPRE

subroutine POT0(NATOMS, V)
   implicit double precision (a-h, o-z)
   parameter (NDA = 160, NDA3 = NDA*3)
   double precision :: bohr_a, ha_ev, kcal_ev
   dimension r(3), dv(3)
   common/QPDOT/Q(NDA3), PDOT(NDA3), FCOEF(NDA3, NDA3)
   common/CONSTN/C1, C2, C3, C4, C5, C6, C7, PI, HALFPI, TWOPI
   bohr_a = 0.529177210903d0        ! 1 bohr [A]
   ha_ev = 27.211386245988d0        ! 1 hartree [eV]
   kcal_ev = 23.0605d0              ! 1 eV [kcal/mol]
   V = 0.0d0
   if (NATOMS /= 3) then
      write (0, '(a,i0)') 'POT0: the BKMP2 seam serves exactly 3 atoms, got NATOMS = ', NATOMS
      stop 1
   end if
   call h3_pairs(r)
   call bkmp2(r, V, dv, 1)
   V = V*ha_ev*kcal_ev*C1
end subroutine POT0

subroutine DPESHON(NATOMS)
   implicit double precision (a-h, o-z)
   parameter (NDA = 160, NDA3 = NDA*3)
   double precision :: bohr_a, ha_ev, kcal_ev, conv
   dimension r(3), dv(3), g(9)
   double precision :: u12, u23, u13
   integer :: k
   common/QPDOT/Q(NDA3), PDOT(NDA3), FCOEF(NDA3, NDA3)
   common/CONSTN/C1, C2, C3, C4, C5, C6, C7, PI, HALFPI, TWOPI
   bohr_a = 0.529177210903d0
   ha_ev = 27.211386245988d0
   kcal_ev = 23.0605d0
   conv = ha_ev*kcal_ev*C1            ! hartree/A -> internal/A (the A-over-bohr
                                       ! chain-rule factors already carry the
                                       ! bohr->A conversion: the unit vectors are
                                       ! built from A differences over bohr lengths)
   if (NATOMS /= 3) then
      write (0, '(a,i0)') 'DPESHON: the BKMP2 seam serves exactly 3 atoms, got NATOMS = ', NATOMS
      stop 1
   end if
   call h3_pairs(r)
   call bkmp2(r, V, dv, 1)
   ! cartesian chain rule: scale dv to [hartree/A] first (engine gives
   ! hartree/bohr), then contract with the A-unit pair vectors
   do k = 1, 3
      u12 = (Q(k) - Q(3 + k))/(r(1)*bohr_a)
      u23 = (Q(3 + k) - Q(6 + k))/(r(2)*bohr_a)
      u13 = (Q(k) - Q(6 + k))/(r(3)*bohr_a)
      g(k) = (dv(1)*u12 + dv(3)*u13)/bohr_a
      g(3 + k) = (-dv(1)*u12 + dv(2)*u23)/bohr_a
      g(6 + k) = (-dv(3)*u13 - dv(2)*u23)/bohr_a
   end do
   do k = 1, 9
      PDOT(k) = -g(k)*conv           ! the force, not the gradient
   end do
end subroutine DPESHON

! h3_pairs(r) - the three pair distances [bohr] of atoms 1,2,3
subroutine h3_pairs(r)
   implicit double precision (a-h, o-z)
   parameter (NDA = 160, NDA3 = NDA*3)
   double precision :: bohr_a
   dimension r(3)
   common/QPDOT/Q(NDA3), PDOT(NDA3), FCOEF(NDA3, NDA3)
   bohr_a = 0.529177210903d0
   r(1) = sqrt((Q(4) - Q(1))**2 + (Q(5) - Q(2))**2 + (Q(6) - Q(3))**2)/bohr_a
   r(2) = sqrt((Q(7) - Q(4))**2 + (Q(8) - Q(5))**2 + (Q(9) - Q(6))**2)/bohr_a
   r(3) = sqrt((Q(7) - Q(1))**2 + (Q(8) - Q(2))**2 + (Q(9) - Q(3))**2)/bohr_a
   if (minval(r) < 0.2d0) then
      write (0, '(a,3f10.4)') 'h3_pairs: pair distance below the BKMP2 domain (0.2 bohr): ', r
      stop 1
   end if
end subroutine h3_pairs
