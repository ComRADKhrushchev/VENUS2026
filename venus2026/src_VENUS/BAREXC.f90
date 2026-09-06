subroutine BAREXC(WT, C, AM, EROT, N, NM)
  use venus_params
  use venus_data
  implicit none
  real(8),  intent(in)    :: WT, C(NDA3yf,NDA3yf), EROT
  real(8),  intent(inout) :: AM(4)
  integer,  intent(in)    :: N, NM
  real(8)  :: QCM(3), VCM(3), RAND, SDUM, PBAR, DUM, DT, DUM1, EINT
  integer  :: I, K, JJ, J, II, ISEED
  real(8), external :: RAND0
!
!         BARRIER EXCITATION
!         SELECT MOMENTA FOR THE REACTION COORDINATE
!
5 format(/4x,'REACTION COORDINATE ENERGY =',1pe18.9,' KCAL/MOL'/)
15 format(15x,'INTERNAL ENERGY =',1pe18.9,' KCAL/MOL'/)
!
!         CALCULATE VELOCITY FOR REACTION COORDINATE
!         THE REACTION COORDINATE IS HELD FIXED
!
!                 FIXED ENERGY FOR REACTION COORDINATE
!
  RAND = RAND0(ISEED)
  if (NBAR /= 2) then
     SDUM = C1*EBAR
     PBAR = sqrt(2*SDUM)
!
!                 FIXED TEMPERATURE FOR REACTION COORDINATE
!
  else if (NBAR == 2) then
     DUM = log(1.0D0-RAND)
     PBAR = sqrt(-2.0D0*C5*TBAR*DUM)
  end if
!
  write(6,5) 0.5D0*PBAR*PBAR/C1

  if (IJDIR /= 0) then
!
!        DIRECT TRAJECTORY TOWARDS PRODUCTS IF IJDIR = 1
!        AND TOWARDS REACTANTS IF IJDIR = -1
!
     DT = 0.1D0
!
!        DT IS AN ARBITRARY, SMALL STEP IN THE "POSITIVE" IRC
!        DIRECTION
!
     do I = 1, N
        do K = 1, 3
           JJ = 3*I+1-K
           J = 3*LL(I)+1-K
           QQ(J) = Q(J)+C(JJ,NM)*DT
        end do
     end do
     DUM1 = 0.0D0
     DUM = 0.0D0
     I = IDIR
     J = JDIR
     II = 3*I
     JJ = 3*J
     do K = 1, 3
        DUM1 = DUM1+(QQ(II)-QQ(JJ))**2
        DUM = DUM+(Q(II)-Q(JJ))**2
        II = II-1
        JJ = JJ-1
     end do
     DUM = sqrt(DUM1)-sqrt(DUM)
     if (IJDIR < 0.0D0 .and. DUM > 0.0D0) PBAR = -PBAR
     if (IJDIR > 0.0D0 .and. DUM < 0.0D0) PBAR = -PBAR
  else
     if (RAND >= 0.5D0) PBAR = -PBAR
  end if
!
!         ADD IRC MOMENTUM AND TRANSFORM FROM NORMAL MODE TO CARTESIAN VELOCITY
!
  do I = 1, N
     do K = 1, 3
        JJ = 3*I+1-K
        J = 3*LL(I)+1-K
        P(J) = P(J)+C(JJ,NM)*PBAR*W(LL(I))
     end do
  end do
!
!         CALCULATE CENTER OF MASS COORDINATES QQ AND MOMENTA PP
!
  call CENMAS(WT, QCM, VCM, N)
!
!         MOVE PP ARRAY TO P ARRAY AND QQ ARRAY TO Q ARRAY
!
  do I = 1, N
     J = 3*LL(I)+1
     do K = 1, 3
        Q(J-K) = QQ(J-K)
        P(J-K) = PP(J-K)
     end do
  end do
  call DVDQ_1
  call ENERGY_1
!
  call ROTN(AM, EROT, N)
!
  EINT = H-EZERO
  write(6,15) EINT
end subroutine BAREXC
