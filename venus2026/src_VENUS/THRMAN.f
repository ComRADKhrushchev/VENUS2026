      SUBROUTINE THRMAN(WW,ANQ,T,NM)
      use venus_params
      use venus_data, only: T_H => T
      use venus_data
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C         CALCULATE VIBRATIONAL QUANTUM NUMBERS FROM A THERMAL
C         (BOLTZMAN) DISTRIBUTION
C
      DIMENSION WW(NDA3),ANQ(NDA3)
C
      DUM1=C7/C5/T
      DO I=1,NM
         DUM=DUM1*WW(I)
         N=GAMA(1,ISEED)/DUM 
         ANQ(I)=DBLE(N)
      ENDDO
      RETURN
      END
