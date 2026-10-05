      SUBROUTINE ENERGY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CALCULATE POTENTIAL, KINETIC AND TOTAL ENERGY OF THE
C         MOLECULAR SYSTEM
C
C      PARAMETER(ND1=2000)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/GLOP/WS1(3),WG1(3),WS2(3),WG2(3),WGS1(3),WGS2(3),WEFF(3)
     &,GSW,FCG,COEFA,COEFB,GN(NDA3)
      COMMON/HFIT/PSCALA,PSCALB,VZERO
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/FORCES/NATOMS,I3N,NFC,NGLO
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/VRSCAL/THERMOTEMP,NSEL,NSCALE,NEQUAL,NRGD
!############################
! Addedn by Zexing Qu 2023.7.10
      COMMON/SH/SH_ENG
!############################
C
	!  real*8::V_MM,ZMat(3,NDA) !#### Added by CZZ, energy contribution of molecular force field
      T=0.0D0
      V=0.0D0
        
C
C       NOW: CALCULATE USER CUSTOMED POTENTIAL ENERGY
C     ADDED BY BIN 2016/10/1
      IF (NSURF.EQ.2) THEN
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
      V=V/C1
	  !#### Added by CZZ, Calculate Z matrix
	!  CALL Zmatrix(ZMat) 
	!  CALL V_MFF(V_MM,ZMat)
	!  V = V + V_MM
	  !####2023/10/30
      H=T+V
      RETURN
      END
