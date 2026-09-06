!THIS SUBROUTINE IS ADDED BY BIN FOR THE GENERALIZED LANGEVIN MODEL
!2017/2/5
      SUBROUTINE GLOINIT
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FACTOR1=413.4D0                                         !-->TIME^-1, CONVERT FRICTION COEFFICIENT TO VENUS'S UNIT (10FS)^-1
      FACTOR2=413.4D0                                         !-->CONVERT FREQUENCY TO VENUS'S UNIT (10FS)^-1
      FCG=FCG*FACTOR1
      DO I=1,3
         WS1(I)=WS1(I)*FACTOR2
         WS2(I)=WS1(I)**2
      ENDDO
      DO I=1,3
         WG1(I)=WG1(I)*FACTOR2
         WG2(I)=WG1(I)**2
      ENDDO
      DO I=1,3
         WGS1(I)=WS1(I)
         WGS2(I)=WS1(I)*WG1(I)
         WEFF(I)=DSQRT(2D0*WS2(I)-WGS1(I)**2)
      ENDDO
      END

      SUBROUTINE GLOSELECT(TVIBB)
      use venus_params
      use venus_data, only: TVIBB_H => TVIBB
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BOLTZ=8.6173D-5*23.0605D0*C1
      RBOHR=0.5291772D0
      FACTOR=DSQRT(TVIBB*BOLTZ*W(NATOMA(1)+1))
      FACTOR2=DSQRT(627.509D0*C1*1822.89D0)
      DO I=1,3
        G1=GASDEV()
        G2=GASDEV()
        G3=GASDEV()
        G4=GASDEV()
         J=3*(NATOMA(1))+I
         K=3*(NATOMA(1)+1)+I
         Q(J)=FACTOR*G1/WEFF(I)/W(NATOMA(1)+1)
         P(J)=FACTOR*G2
         Q(K)=FACTOR*G3/WG1(I)/W(NATOMA(1)+1)+WGS1(I)*Q(J)/WG1(I)
         P(K)=FACTOR*G4
      ENDDO
c      FCG=0D0   !REMOVE DAMPING FORCE AND WHITE NOISE FOR CHECKING THE CONSERVATION, OTHERWISE, COMMENT THIS LINE
      COEFA=(1-FCG*ATIME*0.5D0)/(1+FCG*ATIME*0.5D0)
      COEFB=1D0/(1+FCG*ATIME*0.5D0)
      GSW=DSQRT(2D0*BOLTZ*TVIBB*FCG*W(NATOMA(1)+1)/ATIME)
      END

      SUBROUTINE GLOEQU(TVIBB)
      use venus_params
      use venus_data, only: TVIBB_H => TVIBB
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FT(NDA3)
      SAVE FT

      N=NATOMA(1)

      NC=0

      NSE=1000 ! FOR EQULIBRATION

!.....FOR GLO MODEL, ADDITIONAL FROCES
      CALL DVDQ_1
      DO K=1,3*N
         PDOT(K)=0.0D0
      ENDDO

      DO II=1,3
        J=3*(NATOMA(1))+II
        K=3*(NATOMA(1)+1)+II
        I1=NATOMA(1)+1
        I2=NATOMA(1)+2
        GN(K)=GASDEV()*GSW
        PDOT(J)=PDOT(J)-2D0*WS2(II)*W(I1)*Q(J)+WGS2(II)*W(I1)*Q(K)
        PDOT(K)=PDOT(K)-2D0*WG2(II)*W(I2)*Q(K)+WGS2(II)*W(I2)*Q(J)
        PDOT(K)=PDOT(K)-FCG*P(K)+GN(K)
      ENDDO

      TB=0.D0
      DO I=1,NSE
         NC=NC+1

          DO K=NI-5,NI
            KK=(K+2)/3
            FT(K)=PDOT(K)
            Q(K)=Q(K)+P(K)*ATIME/W(KK)+PDOT(K)/(2.D0*W(KK))*ATIME**2
          ENDDO

!...UPDATE VELOCITY (MOMENTUM)
         CALL DVDQ_1
         DO K=1,3*N
            PDOT(K)=0.0D0
         ENDDO

!.....FOR GLO MODEL
         DO II=1,3
           J=3*(NATOMA(1))+II
           K=3*(NATOMA(1)+1)+II
           I1=NATOMA(1)+1
           I2=NATOMA(1)+2
           GN(K)=GASDEV()*GSW
           PDOT(J)=PDOT(J)-2D0*WS2(II)*W(I1)*Q(J)+WGS2(II)*W(I1)*Q(K)
           PDOT(K)=PDOT(K)-2D0*WG2(II)*W(I2)*Q(K)+WGS2(II)*W(I2)*Q(J)
           PDOT(K)=PDOT(K)-FCG*P(K)+GN(K)
         ENDDO
C
C   P=MV SO P(T+DT)=P(T)+(F(T+DT)+F(T))*DT/2
C
          DO K=1,NI
            P(K)=P(K)+(FT(K)+PDOT(K))*ATIME*0.5D0
          ENDDO

         DO II=N+1,NATOMS
           J3=II*3
           J2=J3-1
           J1=J3-2
           TB=TB+(P(J1)**2+P(J2)**2+P(J3)**2)/W(II)
         ENDDO
         TEMP=TB/(3.0D0*DBLE(NATOMS-N)*0.00198717D0*C1)/NC
         IF (NC.EQ.NC/1000*1000) THEN
            WRITE(30,'(A,2F12.2)')'SELECTED AND CURRENT TEMPERATURE='
     &      ,TVIBB,TEMP
         ENDIF
      ENDDO
      WRITE(30,*)
C
      WRITE(6,*)'EQUALIBRATION FOR GLO MODEL IS NOW OVER'
C
      NC=0
C
      END SUBROUTINE GLOEQU
!.....END BIN, 2017/2/5

      SUBROUTINE  DISTRIBUTE(X, N, RANGE, M, BUCKET)
      IMPLICIT NONE
      REAL*8, DIMENSION(1:N) :: X     ! INPUT SCORE
      INTEGER               :: N     ! # OF SCORES
      REAL*8, DIMENSION(1:M) :: RANGE ! RANGE ARRAY
      INTEGER                :: M     ! # OF RANGES
      INTEGER                            :: I, J  
      REAL*8, DIMENSION(1:M+1) :: BUCKET! COUNTING BUCKET

      DO I = 1, M+1                     ! CLEAR BUCKETS
         BUCKET(I) = 0
      END DO

      DO I = 1, N                       ! FOR EACH INPUT SCORE
         DO J = 1, M-1                  ! DETERMINE THE BUCKET
            IF (X(I).GE.RANGE(J).AND.X(I).LT.RANGE(J+1)) THEN
               BUCKET(J) = BUCKET(J) + 1
            END IF               
         END DO                         ! DON'T FORGET THE LAST BUCKET
         IF (X(I).GE.RANGE(M))  BUCKET(M+1) = BUCKET(M+1)+1
      END DO

      END SUBROUTINE  DISTRIBUTE
