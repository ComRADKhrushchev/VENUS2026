      SUBROUTINE INITQP(WW,A,C,AM,WT,EINT,EROTS,AI,EROT,PHASE,N,NM,
     *NPHASE,IPHASE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         INITIALIZE COORDINATES AND MOMENTA FROM NORMAL MODE
C         PARAMETERS (FREQUENCY, AMPLITUDE,...)
C
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/FINALB/EROTA,EROTB,EA(3),EB(3),AMA(4),AMB(4),AN,AJ,BN,BJ,
     *OAM(4),EREL,ERELSQ,ETCM,BF,SDA,SDB,DELH(NDP),ANG(NDG),NFINAL
      COMMON/LMODEB/ENON,EDELTA,RWANT,PWANT,NEXM,NLEV,JFLAG
      COMMON/FORCES/NATOMS,I3N,NFC,NGLO
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/SADDLE/EBAR,TBAR,EZERO,NBAR,IDIR,JDIR,IJDIR 
!-->..Added by Bin 04/26/2014
      COMMON/HDIAG/HTMIN,HTMAX,ZASYM,DELM
      COMMON/ALIGN/NTHTA
!-->..End
!--> added by Bin, 2014/6/18
      COMMON/NROTAEQ2/JROTA,KROTA,JROTB,KROTB
!--> end
      DIMENSION WW(NDA3),A(NDA3),C(NDA3yf,NDA3yf),QCM(3),VCM(3),AM(4),
     *AI(3),COOR(NDA3),DCOOR(NDA3),PHASE(5),IPHASE(5)
    5 FORMAT('  A-B INTERACTION ENERGY WHEN ENTERING INITQP =',
     *1PE18.9,' KCAL/MOL'/)
   15 FORMAT(15X,'INTERNAL ENERGY =',1PE18.9,' KCAL/MOL')
   25 FORMAT(/,10X,15HCHOSEN:  EROT =,F7.3,9H KCAL/MOL,/,
     *10X,10HJX,JY,JZ =,3D13.5,6H H-BAR,/)
   35 FORMAT(9X,'VIBRATIONAL ANGULAR MOMENTUM =',1PE18.9,' H-BAR')
   45 FORMAT(/,15X,'ONLY KINETIC ENERGY FOR FIRST ',I3,' MODES')
C
      ESEL=EINT
C
C         CALCULATE THE TOTAL ENERGY, WHICH IS THE REFERENCE ENERGY
C         (EZERO) WITH RESPECT TO ADDING EINT.
C
C      CALL DVDQ        ! by bin, 2016/10/02
      CALL ENERGY
      EZERO=H
      WRITE(6,5)EZERO
C
C         SET COUNTER FOR NUMBER OF SCALING ATTEMPTS
C
      NSCALE=0
C
C         SET IFLAG WHICH IS USED FOR LOCAL MODE EXCITATION BETWEEN
C         ATOMS NONI AND NONJ
C
      IFLAG=0
      IF (NACTA.EQ.4) CALL LMODE(1,ENU,EDELTU,ENL,EDELTL)
C
C         STORE THE ANGULAR MOMENTUM VECTOR FROM SELECT
C

      DUM1=AM(1)
      DUM2=AM(2)
      DUM3=AM(3)
C
C         CALCULATE NORMAL MODE COORDINATES AND VELOCITIES
C
C
C!-->    modified by Bin, 09/14/2015
1000  CONTINUE
      NSCALE=0
C!-->    end
C         SET COOR(I) AND DCOOR(I) FOR MODES WHICH ONLY RECEIVE
C         KINETIC ENERGY SPECIFIED BY NUMP
C
      NUMP=0
  48  IF (NUMP .GE. 1) THEN
         DO I = 1,NUMP
            COOR(I)=0.D0
            DCOOR(I)=-WW(I)*A(I)
         ENDDO
      ENDIF
C
C         FOR OTHER MODES
C
C         TO SET THE VIBRATIONAL ANGULAR MOMENTUM FOR DEGENERATE BENDS,
C         THE PHASE OF THE SECOND MODE OF THE DEGENERATE PAIR IS SHIFTED
C         FROM THAT OF THE FIRST MODE BY THE ANGLE "PHASE".
C
      DO I=1+NUMP,NM
          IF (NPHASE.GT.0) THEN
           DO J=1, NPHASE
            IF (IPHASE(J).EQ.I) THEN
            DUM = DUM + PHASE(J)
            GOTO 49
            ENDIF
           ENDDO
          ENDIF
         RAND=RAND0(ISEED)
         DUM=TWOPI*RAND
 49   CONTINUE
         COOR(I)=A(I)*COS(DUM)
         DCOOR(I)=-WW(I)*A(I)*SIN(DUM)
      ENDDO
c
C
C         TRANSFORM FROM NORMAL MODE TO CARTESIAN COORDINATES AND VELOCIT
C
      DO II=1,N
         DO K=1,3
            JJ=3*II+1-K
            J=3*L(II)+1-K
            Q(J)=0.0D0
            P(J)=0.0D0
            DO I=1,NM
               Q(J)=Q(J)+C(JJ,I)*COOR(I)
               P(J)=P(J)+C(JJ,I)*DCOOR(I)
            ENDDO
            P(J)=P(J)*W(L(II))
            Q(J)=Q(J)+QZ(J)
         ENDDO
      ENDDO
C
C         CALCULATE CENTER OF MASS COORDINATES QQ AND MOMENTA PP
C
   50 CALL CENMAS(WT,QCM,VCM,N)

C
C         MOVE PP ARRAY TO P ARRAY AND QQ ARRAY TO Q ARRAY
C
      DO I=1,N
         J=3*L(I)+1
         DO K=1,3
            Q(J-K)=QQ(J-K)
            P(J-K)=PP(J-K)
         ENDDO
      ENDDO

C
C         ADD ANGULAR MOMENTUM VECTOR FROM SELECT TO THE MOLECULE.
C         CALCULATE THE REQUIRED ANGULAR VELOCITY AND ADD IT TO THE
C         MOLECULE.
C
C         IF NPHASE > 0 THERE IS VIBRATIONAL ANGULAR MOMENTUM ABOUT THE
C         X-AXIS OF A LINEAR MOLECULE.  THIS IS ADDED BY SETTING THE
C         QUANTUM NUMBERS AND PHASE FOR THE DEGENERATE BENDS.  THIS
C         ANGULAR MOMENTUM IS NOT SPURIOUS AND SHOULD NOT BE SUBTRACTED.
C
      CALL ROTN(AM,EROT,N)
      AM(1)=DUM1-AM(1)
      AM(2)=DUM2-AM(2)
      AM(3)=DUM3-AM(3)
      NAM=1
      CALL ROTN(AM,EROT,N)
      NAM=0
      WX=-WX
      WY=-WY
      WZ=-WZ
      IF (NPHASE.GT.0) THEN
      WX = 0.0D0
      DUM = AM(1)/C7
      WRITE(6,35)DUM
      ENDIF

      CALL ANGVEL(N)

C
C         SCALE COORDINATES AND MOMENTA TO FIT THE TOTAL ENERGY
C         THE INITIAL CONDITION IS ACCEPTED IF THE CALCULATED AND 
C         DESIRED ENERGY AGREE TO WITHIN 0.1 PER-CENT.
C
C         EINT IS THE SUM OF THE SELECTED INTERNAL VIBRATIONAL AND
C         ROTATIONAL ENERGIES.  DO NOT SCALE IF THE SELECTED VIBRATIONAL
C         ENERGY IS ZERO.
C
      IF (IFLAG.NE.1) THEN
         DDD=EINT-EROTS
         IF (DDD.NE.0.0D0) THEN
!-->    modified by Bin, 4/29/2014
!-->    shift z coorindates up in order to calculate the potential energy
            IF (NSURF.GT.0) then
            J=N
            DO I=1,J
               K3=3*I
               Q(K3)=Q(K3)+ZASYM
            ENDDO
            ENDIF
!-->    end Bin's changes
C            CALL DVDQ  ! by bin, 2016/10/02
            CALL ENERGY
!-->    modified by Bin, 4/29/2014
!-->    shift z coorindates back after the potential energy
            IF (NSURF.GT.0) then
            J=N
            DO I=1,J
               K3=3*I
               Q(K3)=Q(K3)-ZASYM
            ENDDO
            ENDIF
!-->    end Bin's changes
 
            ESEL=H-EZERO
            SDUM=ABS(EINT-ESEL)/EINT
C
C         TEST TO SEE IF CALCULATED ENERGY IS OUT OF RANGE FOR ACCURATE
C         SCALING. IF SO, ONLY ADD KINETIC ENERGY TO AN ADDITIONAL MODE
C
C.........Modified by Bin 12/18/2013
C.........The following four lines are commented by Bin in order to match the Venus96's results of Ar-H2O collision
C.........Their influence on the entire trajectory calculations is not clear 
C            IF (SDUM.GE.0.1D0 .AND. NUMP.LT.NM) THEN
C               NUMP=NUMP+1
C               GOTO 48
C            ENDIF
C.........End

            WRITE(6,15)ESEL
            IF (SDUM.GE.0.001D0) THEN
               NSCALE=NSCALE+1
C               IF (NSCALE.GT.50) STOP
C.........Modified by Bin 09/14/2015
               IF (NSCALE.GT.50) THEN
                  WRITE(6,*)'IMPROPER RANDOM NUMBER, USE THE NEXT ONE'
                  GOTO 1000 
               ENDIF
C.........end
               SDUM=SQRT(EINT/ESEL)
               DO I=1,N
                  J=3*L(I)+1
                  DO K=1,3
                     P(J-K)=P(J-K)*SDUM
                     Q(J-K)=(Q(J-K)-QZ(J-K))*SDUM+QZ(J-K)
                  ENDDO
               ENDDO
               GOTO 50
            ENDIF
            IF (NUMP.GT.0) WRITE(6,45)NUMP
         ENDIF
         IF (NACTA.NE.4) GOTO 120
      ENDIF
C
C         CHOOSE CONDITIONS FOR LOCAL MODE(NACT=4)
C
      CALL LMEXCT
      IFLAG=1
      IF (JFLAG.EQ.0) GOTO 50
C
C         CALCULATE THE TOTAL ENERGY
C

!-->    modified by Bin, 4/29/2014
!-->    shift z coorindates up in order to calculate the potential energy
      IF (NSURF.GT.0) then
      J=N
      DO I=1,J
         K3=3*I
         Q(K3)=Q(K3)+ZASYM
      ENDDO
      ENDIF
!-->    end Bin's changes

C      CALL DVDQ        ! by bin, 2016/10/02
      CALL ENERGY

!-->    modified by Bin, 4/29/2014
!-->    shift z coorindates back after the potential energy
      IF (NSURF.GT.0) then
      J=N
      DO I=1,J
         K3=3*I
         Q(K3)=Q(K3)-ZASYM
      ENDDO
      ENDIF
!-->    end Bin's changes

      ESEL=H-EZERO
C
C         CALCULATE THE ROTATIONAL ENERGY
C
  120 CALL ROTN(AM,EROT,N)
      DUM1=AM(1)/C7
      DUM2=AM(2)/C7
      DUM3=AM(3)/C7
C
      WRITE(6,25)EROT,DUM1,DUM2,DUM3
      WRITE(6,15)ESEL
      WRITE(6,*)
C
      EINT=ESEL

!-->..Modified by Bin 4/24/2014
c      IF (N.EQ.NATOMS) RETURN
c      IF ((NSURF.NE.0).AND.(N.EQ.NATOMB(1))) RETURN

      IF (N.EQ.NATOMS.AND.NSURF.NE.3) RETURN
!-->..end

C
C         RANDOMLY ROTATE THE MOLECULE ABOUT ITS CENTER OF MASS
C         BY EULER'S ANGLES.
C         CENTER OF MASS COORDINATES QQ AND MOMENTA PP ARE PASSED FROM
C         SUBROUTINES CENMAS AND ANGVEL THROUGH COMMON BLOCK WASTE.
C
!-->..Modified by Bin 4/24/2014
!-->..Adapt this to control a specific orientation of initial state
!-->
        IF (NTHTA.GE.0) THEN
           THLAS=90d0
           IF (NTHTA.GT.JROTA) THEN
             WRITE(6,*)'NTHTA CAN NOT EXCEED JROTA'
             STOP
           ENDIF
           CALL ROTATEJKM(NTHTA,JROTA,AM,N,THLAS)
        ELSE
           CALL ROTATE(N)
        ENDIF
C
C         RECALCULATE THE ROTATIONAL ENERGY AND ANGULAR MOMENTUM
C
      CALL ROTN(AM,EROT,N)

C
      RETURN
      END

