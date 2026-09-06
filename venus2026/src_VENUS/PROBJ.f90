subroutine PROBJ(T, AI, ISEED, JFINAL)
  use venus_data, only: WGT
  implicit none
  real(8),  intent(in)    :: T, AI
  integer,  intent(inout) :: ISEED
  integer,  intent(out)   :: JFINAL
  real(8)  :: B, AJPEAK, FAC, PJMPQ, DUM, AJ, PJQ, PF, PCOMP
  integer  :: JMAX, J
  real(8), external :: RAND0
!     COMMON/WNS/WGT — now in venus_data

  B = 48.5085/(2*AI*T)
  call JMAXCALC(T, AI, JMAX, AJPEAK)
!-->    FAC give the maximum weight(odd for orthoH2,even for orthoD2), zlj
  FAC = max(WGT(0), WGT(1))
  PJMPQ = FAC*(2*AJPEAK+1)*exp(-AJPEAK*(AJPEAK+1)*B)
  do
     DUM = RAND0(ISEED)
!-->    modified by Bin, 2/29/2020
     AJ = DUM*dble(JMAX+1)-0.5
     J = nint(AJ)
     if (J < 0) J = 0
     if (J > JMAX) J = JMAX
!-->    sample the integer from -0.5 to JMAX+0.5, zlj
     PJQ = WGT(J)*(2*J+1)*exp(-J*(J+1)*B)
     PF = PJQ/PJMPQ
     PCOMP = RAND0(ISEED)
     if (PF >= PCOMP) exit
  end do
  JFINAL = J
!-->    end
end subroutine PROBJ
