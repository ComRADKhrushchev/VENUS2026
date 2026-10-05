      SUBROUTINE POTENZ(II)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      
C         
C         SETS THE COORDINATES FOR REACTANTS OR PRODUCTS TO THEIR
C         EQUILIBRIUM VALUES, DISPLACES A AND B, AND CALCULATES THE
C         POTENTIAL ENERGY
C
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/HFIT/PSCALA,PSCALB,VZERO
!-->..Added by Bin 04/24/2014
      COMMON/HDIAG/HTMIN,HTMAX,ZASYM,DELM
!-->..End
C
C         INITIALIZE Q AND P ARRAYS
C

!-->..Adapted by Bin 03/22/2016
c      J=3*NATOMA(II)
c      DO I=1,J
c         Q(I)=QZA(II,I)
c         P(I)=0.0D0
c         write(*,*)Q(I)
c      ENDDO
c      K=3*NATOMB(II)       
c      DO I=1,K
c         Q(J+I)=QZB(II,I)
c         P(J+I)=0.0D0
c         write(*,*)Q(J+I)
c      ENDDO
      DO I=1,NATOMA(II)
         J1=3*LA(II,I)
         J2=J1-1
         J3=J1-2
         K1=3*I
         K2=K1-1
         K3=K1-2
         Q(J1)=QZA(II,K1)
         Q(J2)=QZA(II,K2)
         Q(J3)=QZA(II,K3)
         P(J1)=0.0D0
         P(J2)=0.0D0
         P(J3)=0.0D0
      ENDDO
      DO I=1,NATOMB(II)
         J1=3*LB(II,I)
         J2=J1-1
         J3=J1-2
         K1=3*I
         K2=K1-1
         K3=K1-2
         Q(J1)=QZB(II,K1)
         Q(J2)=QZB(II,K2)
         Q(J3)=QZB(II,K3)
         P(J1)=0.0D0
         P(J2)=0.0D0
         P(J3)=0.0D0
      ENDDO
!-->..End

C
C         SEPARATE A AND B BY 1000 ANGSTROMS
C

      IF(NATOMB(II).GT.0)THEN 
C        DO I=1,NATOMB(II)
C          J3=3*LB(II,I)
C          Q(J3)=Q(J3)+1000.0D0
C        ENDDO
!-->    modified by Bin, 7/30/2016
         DO I=1,NATOMA(II)
          J3=3*LA(II,I)
          J2=J3-1
          J1=J2-1
          Q(J3)=Q(J3)+ZASYM
         ENDDO
!-->    end Bin's changes
      ENDIF
C
C         CALCULATE THE A + B POTENTIAL ENERGY, WITH A + B SEPARATED
C         BY 1000 ANGSTROMS.
C       
C      CALL DVDQ        ! commented by bin, 2016/10/02
      CALL ENERGY
C
C         REMOVE THE 1000 ANGSTROM SEPARATION BETWEEN A AND B
C
C  Kyoyeon 11/25/09
      IF(NATOMB(II).GT.0)THEN 
C        DO I=1,NATOMB(II)
C          J3=3*LB(II,I)
C          Q(J3)=Q(J3)-1000.0D0
C        ENDDO
!-->    modified by Bin, 7/30/2016
         DO I=1,NATOMA(II)
          J3=3*LA(II,I)
          J2=J3-1
          J1=J2-1
          Q(J3)=Q(J3)-ZASYM
         ENDDO
!-->    end Bin's changes
      ENDIF
C
      RETURN
      END
