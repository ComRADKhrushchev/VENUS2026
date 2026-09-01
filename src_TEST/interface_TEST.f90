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


! FIXROTDATM: surface-bond alignment (surface builds only). TEST systems
! have no surface bond; reaching this means NZDOWN=1 was set with a TEST PES,
! which is not a supported combination — stop loudly.
subroutine FIXROTDATM(I)
   integer :: I
   write(6,*) 'ERROR: FIXROTDATM stub reached in TEST build ', &
              '(NZDOWN=1 is surface-specific; not valid with TEST_PES)'
   stop
end subroutine FIXROTDATM

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
