      SUBROUTINE ENMODE(EN,NMA)
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         CALCULATE NORMAL MODE ENERGIES
C
      DIMENSION QNCR(NDA),DQNCR(NDA),EN(NDA),A(NDA3yf,NDA3yf)
      N=NATOMA(1)
      I3N=3*N
C
C         SET UP APPROPRIATE TRANSFORMATION MATRIX A
C
      DO I=1,NMA
         DO J=1,N
            K=0
   90       KK=3*J-K
            A(KK,I)=CA(KK,I)*W(J)
            K=K+1
            IF (K.LT.3) GOTO 90
         ENDDO
      ENDDO
C
      DO I=1,N
         J=0
   30    II=3*I-J
         P(II)=P(II)/W(I)
         J=J+1
         IF (J.LT.3) GOTO 30
      ENDDO
C
C         TRANSFORM FROM CARTESIAN COORDINATES AND VELOCITIES
C         INTO NORMAL MODE COORDINATES AND VELOCITIES
C
      DO I=1,NMA
         QNCR(I)=0.0D0
         DQNCR(I)=0.0D0
         DO J=1,I3N
            QNCR(I)=QNCR(I)+A(J,I)*(Q(J)-QZA(1,J))
            DQNCR(I)=DQNCR(I)+A(J,I)*P(J)
         ENDDO
C
C         CALCULATE FOR THE NORMAL MODE ENERGIES
C
         EN(I)=(DQNCR(I)*DQNCR(I)+WWA(I)*WWA(I)*QNCR(I)*QNCR(I))/2.0D0
         EN(I)=EN(I)/C1
      ENDDO
C
C         TRANSFORM ALL P'S BACK AND STORE THEM
C
      DO I=1,N
         J=0
   60    II=3*I-J
         P(II)=P(II)*W(I)
         J=J+1
         IF (J.LT.3) GOTO 60
      ENDDO
      RETURN
      END
