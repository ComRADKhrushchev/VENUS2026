subroutine ORTHAN(AM, WT, ENMT, HSCALE, QMAX, QMIN, PMAX, PSCALE, EROT, N)
  use venus_params
  use venus_data
  implicit none
  real(8),  intent(inout) :: AM(4), EROT
  real(8),  intent(in)    :: WT, ENMT, HSCALE, PSCALE
  real(8),  intent(inout) :: QMAX(NDA3), QMIN(NDA3), PMAX(NDA)
  integer,  intent(in)    :: N
  real(8)  :: RV(3*NDA3), QCM(3), VCM(3)
  real(8)  :: DUM1, DUM2, DUM3, SDUM, SUMM, SUM, RAND, XS, PRO
  integer  :: I, J, K, NDIM, JJ, ISEED
  real(8), external :: RAND0
!
!         INITIALIZE COORDINATES AND MOMENTA FROM ORTHANT SAMPLING
!
27 format(15x,'INTERNAL ENERGY=',1pe18.9,' KCAL/MOL')
!
  call DVDQ_1
  call ENERGY_1
  EZERO = H
  NSCALE = 0
!
!         STORE THE ANGULAR MOMENTUM VECTOR FROM SELECT
!
  DUM1 = AM(1)
  DUM2 = AM(2)
  DUM3 = AM(3)
!
!         CALCULATE INITIAL CONDITIONS USING ORTHANT SAMPLING.
!         FIRST TIME THROUGH SELECT NSFLAG=0, AND THE QMAX, QMIN, PMAX,
!         AND PMIN ARRAYS MUST BE CALCULATED.
!
  if (NSFLAG /= 1) then
     do I = 1, N
        J = 3*LL(I)+1
        do K = 1, 3
           P(J-K) = 0.0D0
           Q(J-K) = QZ(J-K)
        end do
     end do
     do I = 1, N
        J = 3*LL(I)-3
        do K = 1, 3
           do
              Q(J+K) = Q(J+K)+0.1D0
              call DVDQ_1
              call TEST
              if (NTEST == 0) then
                 call ENERGY_1
                 H = H-EZERO
                 if (H < ENMT) cycle
              end if
              exit
           end do
           QMAX(J+K) = Q(J+K)
           Q(J+K) = QZ(J+K)
           do
              Q(J+K) = Q(J+K)-0.1D0
              call DVDQ_1
              call TEST
              if (NTEST == 0) then
                 call ENERGY_1
                 H = H-EZERO
                 if (H < ENMT) cycle
              end if
              exit
           end do
           QMIN(J+K) = Q(J+K)
           Q(J+K) = QZ(J+K)
        end do
     end do
     SDUM = ENMT*C1*2.0D0
     do I = 1, N
        PMAX(LL(I)) = sqrt(SDUM*W(LL(I)))*PSCALE
     end do
  end if
!
!         CALCULATE 3N DIMENSIONAL RANDOM UNIT VECTOR
!
  NDIM = 6*N-2
  SUMM = 1.0D0
  SUM = 1.0D0
  J = 1
  do
     RAND = RAND0(ISEED)
     XS = RAND*SUM
     SDUM = SUMM-XS*XS
     PRO = dble(NDIM-1)/2.0D0
     PRO = (SDUM/SUMM)**PRO
     RAND = RAND0(ISEED)
     if (PRO < RAND) cycle
     RV(J) = XS
     SUMM = SDUM
     SUM = sqrt(SUMM)
     NDIM = NDIM-1
     J = J+1
     if (NDIM > 0) cycle
     exit
  end do
  RAND = RAND0(ISEED)
  XS = SUM*sin(HALFPI*RAND)
  RV(J) = XS
  J = J+1
  RV(J) = sqrt(SUMM-XS*XS)
!
!         SELECT MOMENTA
!
  J = 1
  do I = 1, N
     K = 3*LL(I)
     P(K-2) = RV(J)*PMAX(LL(I))
     J = J+1
     P(K-1) = RV(J)*PMAX(LL(I))
     J = J+1
     P(K) = RV(J)*PMAX(LL(I))
     J = J+1
  end do
  do I = 1, N
     J = 3*LL(I)+1
     do K = 1, 3
        RAND = RAND0(ISEED)
        if (RAND < 0.50D0) P(J-K) = -P(J-K)
     end do
  end do
!
!         SELECT COORDINATES
!
  JJ = 3*N+1
  do I = 1, N
     J = 3*LL(I)-3
     do K = 1, 3
        RAND = RAND0(ISEED)
        if (RAND >= 0.50D0) then
           Q(J+K) = (QMAX(J+K)-QZ(J+K))*RV(JJ)+QZ(J+K)
        else
           Q(J+K) = (QMIN(J+K)-QZ(J+K))*RV(JJ)+QZ(J+K)
        end if
        JJ = JJ+1
     end do
  end do
!
!         SUBTRACT OFF CENTER OF MASS VELOCITY
!
  do
     call CENMAS(WT, QCM, VCM, N)
     do I = 1, N
        J = 3*LL(I)+1
        do K = 1, 3
           P(J-K) = PP(J-K)
        end do
     end do
!
!         ADD ANGULAR VELOCITY VECTOR TO THE MOLECULE TO FIT THE
!         TOTAL ROTATIONAL ENERGY
!
!         ADD ANGULAR MOMENTUM VECTOR FROM SELECT TO THE MOLECULE.
!         CALCULATE THE REQUIRED ANGULAR VELOCITY AND ADD IT TO THE
!         MOLECULE.
!
     call ROTN(AM, EROT, N)
     AM(1) = DUM1-AM(1)
     AM(2) = DUM2-AM(2)
     AM(3) = DUM3-AM(3)
     NAM = 1
     call ROTN(AM, EROT, N)
     NAM = 0
     WX = -WX
     WY = -WY
     WZ = -WZ
     call ANGVEL(N)
!
!         SCALE COORDINATES AND MOMENTA TO FIT THE TOTAL ENERGY
!
     call DVDQ_1
     call ENERGY_1
     H = H-EZERO
     write(6,27) H
     SDUM = abs(HSCALE-H)/HSCALE
     if (SDUM < 0.001D0) exit
     NSCALE = NSCALE+1
     if (NSCALE > 50) stop
     SDUM = sqrt(HSCALE/H)
     do I = 1, N
        J = 3*LL(I)+1
        do K = 1, 3
           P(J-K) = P(J-K)*SDUM
           Q(J-K) = (Q(J-K)-QZ(J-K))*SDUM+QZ(J-K)
        end do
     end do
  end do
!
  H = H+EZERO
!
!         RANDOMLY ROTATE THE MOLECULE ABOUT ITS CENTER OF MASS BY
!         EULER'S ANGLES.  CENTER OF MASS COORDINATES QQ AND MOMENTA PP
!         ARE PASSED FROM SUBROUTINES CENMAS AND ANGVEL THROUGH COMMON
!         BLOCK WASTE.
!
  if (N /= NATOMS) call ROTATE(N)
!
  call ROTN(AM, EROT, N)
end subroutine ORTHAN
