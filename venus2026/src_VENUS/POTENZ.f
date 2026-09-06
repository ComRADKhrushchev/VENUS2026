      SUBROUTINE POTENZ(II)
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      ZASYM = 10d0
C
C         SETS THE COORDINATES FOR REACTANTS OR PRODUCTS TO THEIR
C         EQUILIBRIUM VALUES, DISPLACES A AND B, AND CALCULATES THE
C         POTENTIAL ENERGY
C
!!    ADDED BY BIN, 2018/07/09
!!    COMMON NSID,NDELH(NDP),MPATH now in venus_data
!!-->.END


C
C         INITIALIZE Q AND P ARRAYS
C

!-->..Adapted by Bin 09/07/2018
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

      IF(NATOMB(II).GT.0.AND.NDELH(II).EQ.0)THEN 
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
      CALL ENERGY_1
      !ADDED BY Meng
      write(666,'(4I6)') II,NATOMA(II),NATOMB(II),NDELH(II)
      write(666,'(3F15.8)') Q(1),Q(2),Q(3)
      write(666,'(3F15.8)') Q(4),Q(5),Q(6)
      write(666,'(3F15.8)') U0,U1,V01
      write(666,'(F15.8)') VZERO
      write(666,'(F15.8)') V
      write(666,*)
      !ENDED
C
C         REMOVE THE 1000 ANGSTROM SEPARATION BETWEEN A AND B
C
C  Kyoyeon 11/25/09
      IF(NATOMB(II).GT.0.AND.NDELH(II).EQ.0)THEN 
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
