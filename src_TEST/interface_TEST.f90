!******************************************************************************
! interface_TEST — PES interface for the TEST potential suite.
!
! Provides the same three-symbol contract every PES interface supplies to the
! VENUS core (POTPRE / POT0 / DPESHON); conversion chain:
!   POT0    : V[eV] * 23.0605            -> kcal/mol
!   DPESHON : PDOT += -grad[eV/bohr] * 23.0605 * C1
! No IESH/TDHF symbols: classical (ELEC_METHOD=ADIABATIC) use only.
! PES selection is a runtime keyword (TEST_PES), so one build serves all
! three test potentials — no build-variant proliferation.
!******************************************************************************
subroutine POTPRE
   use test_potentials
   implicit none
   call init_test_potentials()
end subroutine POTPRE

subroutine POT0(NDUM, Vpot, Qarr)
   use venus_data
   use test_potentials
   implicit real*8 (a-h, o-z)
   real*8, intent(in) :: Qarr(*)
   real*8 :: v_ev
   ! evaluate on live Q (Q, not the dummy Qarr)
   call test_pot_v(NATOMS, Q(1:3*NATOMS), v_ev)
   Vpot = v_ev * 23.0605d0
   ! Raw engine E_int (before the VZERO reactant shift) in kcal/mol, so the
   ! fort.1001 'E0(eV)' diagnostic column matches the 2D interface behaviour.
   E0 = Vpot
end subroutine POT0

subroutine DPESHON(NDUM, Qarr)
   use venus_params, only: C1, NDA3
   use venus_data
   use test_potentials
   implicit real*8 (a-h, o-z)
   real*8, intent(in) :: Qarr(*)
   real*8 :: v_ev, g(3*NDA3)
   integer :: i
   call test_pot_vg(NATOMS, Q(1:3*NATOMS), v_ev, g)
   do i = 1, 3*NATOMS
      PDOT(i) = PDOT(i) - g(i)*23.0605d0*C1
   end do
end subroutine DPESHON


! FIXROTDATM: align fragment A's atom1-atom2 axis with the z direction
! (surface-oscillator / surface-bond builds). Ported from the 2D interface.
! Degenerate case (coincident atoms, e.g. a GLO spring pair at equilibrium)
! is a no-op - there is no axis to align.
subroutine FIXROTDATM(I)
   use venus_data, only: Q, W
   integer :: I
   real(8) :: VECZ(3), BDVC(3), VECX(3), BDVCXY(3)
   real(8) :: COSPHI, PHI, COSTHTA, nxy, n3
   VECX = (/1.0d0, 0.0d0, 0.0d0/)
   VECZ = (/0.0d0, 0.0d0, 1.0d0/)
   BDVC = Q(4:6) - Q(1:3)
   n3 = sqrt(dot_product(BDVC, BDVC))
   if (n3 < 1.0d-8) return          ! coincident pair: nothing to align
   BDVCXY(1:2) = BDVC(1:2)
   nxy = sqrt(dot_product(BDVCXY, BDVCXY))
   if (nxy > 1.0d-8) then
      COSPHI = dot_product(BDVCXY, VECX) / nxy
      PHI = acos(max(-1.0d0, min(1.0d0, COSPHI)))
      ! rotate about z by -+PHI to bring the pair into the xz plane
      call ALIGN_ZSTEP(PHI)
   end if
   BDVC = Q(4:6) - Q(1:3)
   COSTHTA = dot_product(BDVC, VECZ) / sqrt(dot_product(BDVC, BDVC))
   ! tilt onto the z axis (rotation about y by the complement angle)
   call ALIGN_YSTEP(COSTHTA)
   write(6,*) '#### ENFORCED ROTATION ACTIVATED ####'
end subroutine FIXROTDATM

! Rotate atoms 1-2 about the z axis by angle PHI (rigid, about origin).
subroutine ALIGN_ZSTEP(PHI)
   use venus_data, only: Q, P
   real(8) :: PHI, CP, SP, x, y, px, py
   integer :: J
   CP = cos(PHI); SP = sin(PHI)
   do J = 1, 2
      x = Q(3*J-2); y = Q(3*J-1)
      Q(3*J-2) =  x*CP + y*SP
      Q(3*J-1) = -x*SP + y*CP
      px = P(3*J-2); py = P(3*J-1)
      P(3*J-2) =  px*CP + py*SP
      P(3*J-1) = -px*SP + py*CP
   end do
end subroutine ALIGN_ZSTEP

! Rotate atoms 1-2 about the y axis so that the pair axis cos matches
! COSTHTA with the z axis (sign-aware tilt).
subroutine ALIGN_YSTEP(COSTHTA)
   use venus_data, only: Q, P
   real(8) :: COSTHTA, SNT, CST, x, z, px, pz
   integer :: J
   CST = COSTHTA
   SNT = sqrt(max(0.0d0, 1.0d0 - CST*CST))
   do J = 1, 2
      x = Q(3*J-2); z = Q(3*J)
      Q(3*J-2) = x*CST - z*SNT
      Q(3*J)   = x*SNT + z*CST
      px = P(3*J-2); pz = P(3*J)
      P(3*J-2) = px*CST - pz*SNT
      P(3*J)   = px*SNT + pz*CST
   end do
end subroutine ALIGN_YSTEP

! GASDEV: normal deviate (Box-Muller). Called by THERMO/GLO and the VERLET
! thermostat path; src_VENUS/GASDEV.f provides the equivalent grandom().
function gasdev() RESULT(gauss_random)
real*8 :: gauss_random
INTEGER, SAVE :: iset = 0
real*8, SAVE :: gset
real*8 :: fac, rsq, v1, v2, u(2)

if (iset == 0) then
    do
        call random_number(u)
        v1 = 2.0*u(1) - 1.0
        v2 = 2.0*u(2) - 1.0
        rsq = v1*v1 + v2*v2
        if (rsq < 1.0 .and. rsq > 0.0) exit
    end do

    fac = sqrt(-2.0*log(rsq)/rsq)
    gset = v1*fac
    gauss_random = v2*fac
    iset = 1
else
    gauss_random = gset
    iset = 0
endif
if(isnan(gauss_random)) write(6,*) 'GASDEV_RETURNED_NAN iset_was=',1-iset,' fac=',fac,' rsq=',rsq,' v1=',v1,' v2=',v2
END FUNCTION gasdev
