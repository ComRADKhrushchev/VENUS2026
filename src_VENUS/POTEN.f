      SUBROUTINE POTEN(THETAA,PHI,QCM,N)
      use venus_params
      use venus_data, only: PHI_H => PHI, L => LL
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         ROTATE MOLECULE BY THETAA AND PHI ON ITS CENTER OF MASS AND
C         CALCULATE POTENTIAL ENERGY
C
C     COMMON/QPDOT/ — now in venus_data
C     COMMON/WASTE/ — now in venus_data
C     COMMON/COORS/ — now in venus_data
C     COMMON/PRLIST/ — now in venus_data
C     COMMON/CONSTN/ — constants now in venus_params
      DIMENSION QCM(3)
C
C          CALCULATE ANGLES
C
      THETH=THETAA-HALFPI
      STHET=SIN(THETH)
      CTHET=COS(THETH)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
C
C          ROTATE FRAGMENT ABOUT ITS CENTER OF MASS.  FIRST ROTATE BY
C          THETAA ABOUT THE X-AXIS AND THEN BY PHI ABOUT THE Z-AXIS.
C          THETAA = 0-PI, AND PHI = 0-TWOPI.
C          THE REACTION PATH IS ALONG THE Y-AXIS.
C
      DO I=1,N
         J=3*L(I)
         Q(J)=STHET*QQ(J-1)+CTHET*QQ(J)+QCM(3)
         Q(J-1)=-SPHI*QQ(J-2)+CPHI*CTHET*QQ(J-1)-CPHI*STHET*QQ(J)
     *          +QCM(2)
         Q(J-2)=CPHI*QQ(J-2)+SPHI*CTHET*QQ(J-1)-SPHI*STHET*QQ(J)+QCM(1)
      ENDDO
      CALL DVDQ_1
      CALL ENERGY_1
      RETURN
      END
