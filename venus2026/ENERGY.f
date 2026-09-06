      SUBROUTINE ENERGY_1
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         CALCULATE POTENTIAL, KINETIC AND TOTAL ENERGY OF THE
C         MOLECULAR SYSTEM
C     All data formerly in COMMON blocks — now in venus_data / venus_params
C
      T=0.0D0
      V=0.0D0
        
C
C       NOW: CALCULATE USER CUSTOMED POTENTIAL ENERGY
C     ADDED BY BIN 2016/10/1
      IF (NSURF.EQ.1) THEN
c         call pbc(natoms)
         CALL POT0(NATOMS,VV)
!#################################
! Added by Zexing Qu 2023.7.10
!         CALL POT_SH(NATOMS,VV,SH_ENG)
!#################################
      ELSE
         IF (NGLO.EQ.0) THEN
            CALL POT0(NATOMS,VV)

!#################################
! Added by Zexing Qu 2023.7.10
!         CALL POT_SH(NATOMS,VV,SH_ENG)
!#################################

         ELSE
            CALL POT0(NATOMS,VV)

!#################################
! Added by Zexing Qu 2023.7.10
!         CALL POT_SH(NATOMS,VV,SH_ENG)
!#################################

            DO I=1,3
            I1=NATOMA(1)+1
            I2=NATOMA(1)+2
            J=3*(NATOMA(1))+I
            K=3*(NATOMA(1)+1)+I
            VV=VV+WS2(I)*W(I1)*Q(J)**2
     &   +WG2(I)*W(I2)*Q(K)**2-WGS2(I)*W(I1)*Q(J)*Q(K)
            ENDDO
         ENDIF
      endif
C       BY BIN 12/18/2013
C
C         ADD VZERO TO THE POTENTIAL ENERGY
C
      V=VV+VZERO
C
C         CALCULATE KINETIC ENERGY
C
      J=1
      DO I=1,NATOMS
         T=T+(P(J)**2+P(J+1)**2+P(J+2)**2)/2.0/W(I)
         J=J+3
      ENDDO
C
C         CONVERT ENERGY TO KCAL/MOLE
C
      T=T/C1
  !   V=V/C1  !Now POT0 directly gives kcal/mol  2024/7/19
      H=T+V
      RETURN
      END
