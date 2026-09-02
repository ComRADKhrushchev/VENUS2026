subroutine INITEBK(N, J, RMIN, RMAX, DH, RMASS, ENJ, PTEST, AL)
  use venus_params, only: C7
  use venus_data, only: W, LBOND
  implicit none
  integer,  intent(in)    :: N, J
  real(8),  intent(in)    :: DH
  real(8),  intent(inout) :: RMIN, RMAX, ENJ
  real(8),  intent(out)   :: RMASS, PTEST, AL
  real(8)  :: BN, HNU, AM, DUM, AN, AJ
  integer  :: ICOUNT
!
!         INITIALIZE PARAMETERS FOR AN OSCILLATOR WITH GIVEN
!         QUANTUM NUMBERS N AND J BY SEMICLASSICAL EBK QUANTIZATION
!
  BN = dble(N)
  HNU = ENJ/(BN+0.5D0)
  AM = sqrt(dble(J*(J+1)))
  AL = AM*C7
  RMASS = W(LBOND(1))*W(LBOND(2))/(W(LBOND(1))+W(LBOND(2)))
  DUM = 1.0D0
!
!         SOLVE FOR ENJ BY FIXED POINT APPROACH
!
  ICOUNT = 0
  do while (abs(DUM) > 1.0D-6)
     call FINLNJ(ENJ, AM, RMIN, RMAX, DH, AN, AJ)
     DUM = BN-AN
     ENJ = ENJ+DUM*HNU
     ICOUNT = ICOUNT+1

!  icount changed from 200 to 1000 8/1/2010
     if (ICOUNT > 1000) stop
  end do
!
  RMIN = RMIN+0.001D0
  RMAX = RMAX-0.001D0
  PTEST = sqrt(0.0001D0*2.0D0*RMASS*ENJ)
end subroutine INITEBK
