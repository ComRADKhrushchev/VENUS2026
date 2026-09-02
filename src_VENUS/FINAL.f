      SUBROUTINE FINAL
      use venus_params
      use venus_data, only: L => LBOND
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         SUBROUTINE FINAL CALCULATES THE PRODUCT ENERGIES
C         AND SCATTERING ANGLES
C
!!    ADDED BY BIN, 2018/07/09
!!    COMMON NSID,NDELH(NDP),MPATH
!!-->.END

      DIMENSION QCMA(3),VCMA(3),QCMB(3),VCMB(3),QR(3),VR(3),QTEM(NDA3)
C
C         IF NFINAL=1, ONLY CALCULATE A AND B EROT

C.....modified by Bin, 7/10/2018
      IF (NSURF.GT.0.AND.NDELH(NPATH).NE.0) GOTO 999
C.....end
C
C         SET PROPERTIES FOR ATOM A
C
      IF (NATOMA(NPATH).LE.1) THEN
         EA(1)=0.0D0
         EA(2)=0.0D0
         AMA(1)=0.0D0
         AMA(2)=0.0D0
         AMA(3)=0.0D0
         AMA(4)=0.0D0
         EROTA=0.0D0
         GOTO 50
      ENDIF
C
C         CALCULATE CENTER OF MASS FOR FRAGMENT B
C
         WT=WTB(NPATH)
         N=NATOMB(NPATH)
         DO I=1,N
            LL(I)=LB(NPATH,I)
         ENDDO
C.....modified by Bin, 9/30/2016
      IF (NSURF.EQ.0) THEN
         CALL CENMAS(WT,QCMB,VCMB,N)
C
C            SET Q FOR FRAGMENT B TO QZB + QCMB.  THE ORIGIN FOR QZB IS 0.
C            PUT COORDINATES FOR B IN TEMPORARY STORAGE QTEM.
C
         DO I=1,N
            K3=3*I
            K2=K3-1
            K1=K2-1
            J3=3*LL(I)
            J2=J3-1
            J1=J2-1
            QTEM(J3)=Q(J3)
            QTEM(J2)=Q(J2)
            QTEM(J1)=Q(J1)
            Q(J3)=QZB(NPATH,K3)+QCMB(3)+ZASYM !+ 500.0D0
            Q(J2)=QZB(NPATH,K2)+QCMB(2)+ZASYM !+ 500.0D0
            Q(J1)=QZB(NPATH,K1)+QCMB(1)+ZASYM !+ 500.0D0
         ENDDO
C.....modified by Bin, 2/4/2017
       ELSE 
         DO I=1,N
            K3=3*I
            K2=K3-1
            K1=K2-1
            J3=3*LL(I)
            J2=J3-1
            J1=J2-1
            QTEM(J3)=Q(J3)
            QTEM(J2)=Q(J2)
            QTEM(J1)=Q(J1)
            Q(J3)=QZB(NPATH,K3) !+ 500.0D0
            Q(J2)=QZB(NPATH,K2) !+ 500.0D0
            Q(J1)=QZB(NPATH,K1) !+ 500.0D0
         ENDDO
         IF (NGLO.NE.0) THEN
            DO I=NATOMA(1)+1,NATOMS
               J3=3*I
               J2=J3-1
               J1=J2-1
               QTEM(J3)=Q(J3)
               QTEM(J2)=Q(J2)
               QTEM(J1)=Q(J1)
               Q(J3)=0d0 
               Q(J2)=0d0 
               Q(J1)=0d0 
            ENDDO
          ENDIF
      ENDIF
C.....end

C
C         CALCULATE PROPERTIES FOR FRAGMENT A
C
C             ANGULAR MOMENTUM AND ROTATIONAL ENERGY
C
      WT=WTA(NPATH)
      N=NATOMA(NPATH)
      DO I=1,N
         LL(I)=LA(NPATH,I)
      ENDDO
      CALL CENMAS(WT,QCMA,VCMA,N)
      CALL ROTN(AMA,EROTA,N)
      DO I=1,4
         AMA(I)=AMA(I)/C7
      ENDDO
C
C             INTERNAL KINETIC ENERGY
C
      IF (NFINAL.NE.1) THEN
         EA(1)=0.0D0
         DO I=1,N
            J=LL(I)
            J3=3*J
            J2=J3-1
            J1=J2-1
            EA(1)=EA(1)+(PP(J1)**2+PP(J2)**2+PP(J3)**2)/W(J)
         ENDDO
         EA(1)=EA(1)/C1/2.0D0
C
C             INTERNAL POTENTIAL ENERGY AND TOTAL ENERGY
C
         CALL DVDQ_1
         CALL ENERGY_1

         EA(2)=V-DELH(NPATH)
         EA(3)=EA(1)+EA(2)
         write(777,'(3F15.8)') Q(1),Q(2),Q(3)
         write(777,'(3F15.8)') Q(4),Q(5),Q(6)
         write(777,'(3F15.8)') U0,U1,V01
         write(777,'(I6,5F15.8)') NPATH,DELH(NPATH),V,EA(1),EA(2),EROT
         write(777,*)

C
C             CALCULATE n AND j IF A IS A DIATOM
C
         IF ((N.EQ.2).AND.(NTEST.EQ.2)) THEN
            DUM=EA(3)*C1
            DUM1=AMA(4)
            CALL FINLNJ(DUM,DUM1,RMIN,RDUM,DELH(NPATH),AN,AJ)
         ENDIF
C
C.....Normal mode analysis added by Bin, 3/23/2015
C
C
C         SUBTRACT OFF THE ANGULAR VELOCITY
C
            DO I=1,N
               J=LL(I)
               J3=3*J
               J2=J3-1
               J1=J2-1
               PP(J1)=PP(J1)-(QQ(J3)*WY-QQ(J2)*WZ)*W(J)
               PP(J2)=PP(J2)-(QQ(J1)*WZ-QQ(J3)*WX)*W(J)
               PP(J3)=PP(J3)-(QQ(J2)*WX-QQ(J1)*WY)*W(J)
            ENDDO
C
C         SAVE THE FINAL COORDINATES AND MOMENTA
C
            DO I=1,N
               J=LL(I)
               J3=J*3
               J2=J3-1
               J1=J2-1
               QFINAL(J1)=QQ(J1)
               QFINAL(J2)=QQ(J2)
               QFINAL(J3)=QQ(J3)
               PFINAL(J1)=PP(J1)
               PFINAL(J2)=PP(J2)
               PFINAL(J3)=PP(J3)
            ENDDO
C.....End
      ENDIF
!C.....commented by Bin, 2016/10/10
!C.....modified by Bin, 3/19/2015
!      IF (NSURF.EQ.3) THEN
!      DO I=1,N
!         J3=3*LA(NPATH,I)
!         J2=J3-1
!         J1=J3-2
!         Q(J3)=QTEM(J3)
!         Q(J2)=QTEM(J2)
!         Q(J1)=QTEM(J1)
!      ENDDO
!      ENDIF
!C.....end
!
C
C         RESET Q ARRAY TO CURRENT B COORDINATES
C
C.....modified by Bin, 2/4/2017
      IF (NGLO.NE.0) THEN
         DO I=NATOMA(1)+1,NATOMS
            J3=3*I
            J2=J3-1
            J1=J3-2
            Q(J3)=QTEM(J3)
            Q(J2)=QTEM(J2)
            Q(J1)=QTEM(J1)
         ENDDO
      ENDIF

      N=NATOMB(NPATH)
      DO I=1,N
         J3=3*LB(NPATH,I)
         J2=J3-1
         J1=J2-1
         Q(J3)=QTEM(J3)
         Q(J2)=QTEM(J2)
         Q(J1)=QTEM(J1)
      ENDDO
C.....end
      IF (NSURF.GT.0) GOTO 60
C
C         CALCULATE CENTER OF MASS FOR FRAGMENT A
C
   50 CONTINUE
      WT=WTA(NPATH)
      N=NATOMA(NPATH)
      DO I=1,N
         LL(I)=LA(NPATH,I)
      ENDDO
      CALL CENMAS(WT,QCMA,VCMA,N)
C
C         SET PROPERTIES FOR ATOM B
C
      IF (NATOMB(NPATH).LE.1) THEN
         EB(1)=0.0D0
         EB(2)=0.0D0
         AMB(1)=0.0D0
         AMB(2)=0.0D0
         AMB(3)=0.0D0
         AMB(4)=0.0D0
         EROTB=0.0D0
         IF (NFINAL.EQ.1) GOTO 999
         GOTO 60
      ENDIF
C
C         SET Q FOR FRAGMENT A TO QZA + QCMA.  THE ORIGIN FOR QZA IS 0.
C         PUT COORDINATES FOR A IN TEMPORARY STORAGE QTEM.
C
      DO I=1,N
         K3=3*I
         K2=K3-1
         K1=K2-1
         J3=3*LL(I)
         J2=J3-1
         J1=J2-1
         QTEM(J3)=Q(J3)
         QTEM(J2)=Q(J2)
         QTEM(J1)=Q(J1)
         Q(J3)=QZA(NPATH,K3)+QCMA(3)+ZASYM !+ 500.0D0
         Q(J2)=QZA(NPATH,K2)+QCMA(2)+ZASYM !+ 500.0D0
         Q(J1)=QZA(NPATH,K1)+QCMA(1)+ZASYM !+ 500.0D0
      ENDDO
C
C         CALCULATE PROPERTIES FOR FRAGMENT B
C
C             ANGULAR MOMENTUM AND ROTATIONAL ENERGY
C
      WT=WTB(NPATH)
      N=NATOMB(NPATH)
      DO I=1,N
         LL(I)=LB(NPATH,I)
      ENDDO
      CALL CENMAS(WT,QCMB,VCMB,N)
      CALL ROTN(AMB,EROTB,N)
      DO I=1,4
         AMB(I)=AMB(I)/C7
      ENDDO
C
C             INTERNAL KINETIC ENERGY
C
      IF (NFINAL.NE.1) THEN
         EB(1)=0.0D0
         DO I=1,N
            J=LL(I)
            J3=3*J
            J2=J3-1
            J1=J2-1
            EB(1)=EB(1)+(PP(J1)**2+PP(J2)**2+PP(J3)**2)/W(J)
         ENDDO
         EB(1)=EB(1)/C1/2.0D0
C
C             INTERNAL POTENTIAL ENERGY AND TOTAL ENERGY
C
         CALL DVDQ_1
         CALL ENERGY_1
         EB(2)=V-DELH(NPATH)
         EB(3)=EB(1)+EB(2)
C
C             CALCULATE FINAL n AND j IF B IS A DIATOM
C
         IF ((N.EQ.2).AND.(NTEST.EQ.2)) THEN
            DUM=EB(3)*C1
            DUM1=AMB(4)
            CALL FINLNJ(DUM,DUM1,RMIN,RDUM,DELH(NPATH),BN,BJ)
         ENDIF
C
C.....Normal mode analysis added by Bin, 3/23/2015
C
C
C         SUBTRACT OFF THE ANGULAR VELOCITY
C
            DO I=1,N
               J=LL(I)
               J3=3*J
               J2=J3-1
               J1=J2-1
               PP(J1)=PP(J1)-(QQ(J3)*WY-QQ(J2)*WZ)*W(J)
               PP(J2)=PP(J2)-(QQ(J1)*WZ-QQ(J3)*WX)*W(J)
               PP(J3)=PP(J3)-(QQ(J2)*WX-QQ(J1)*WY)*W(J)
            ENDDO
C
C         SAVE THE FINAL COORDINATES AND MOMENTA
C
            DO I=1,N
               J=LL(I)
               J3=J*3
               J2=J3-1
               J1=J2-1
               QFINAL(J1)=QQ(J1)
               QFINAL(J2)=QQ(J2)
               QFINAL(J3)=QQ(J3)
               PFINAL(J1)=PP(J1)
               PFINAL(J2)=PP(J2)
               PFINAL(J3)=PP(J3)
            ENDDO
C.....End
      ENDIF
C
C         RESET Q ARRAY TO CURRENT A COORDINATES
C
      N=NATOMA(NPATH)
      DO I=1,N
         J3=3*LA(NPATH,I)
         J2=J3-1
         J1=J2-1
         Q(J3)=QTEM(J3)
         Q(J2)=QTEM(J2)
         Q(J1)=QTEM(J1)
      ENDDO
      IF (NFINAL.EQ.1) GOTO 999
C
C         CALCULATE RELATIVE A+B PROPERTIES
C
C
   60 CONTINUE
C
C            CALCULATE CENTER OF MASS TRANSLATIONAL ENERGY FOR THE 
C            COMPLETE A + B SYSTEM
C
      DUMX=0.0D0
      DUMY=0.0D0
      DUMZ=0.0D0
      DUMW=0.0D0
      DO I=1,NATOMS
         J3=3*I
         J2=J3-1
         J1=J2-1
         DUMX=DUMX + P(J1)
         DUMY=DUMY + P(J2)
         DUMZ=DUMZ + P(J3)
         DUMW=DUMW + W(I)
      ENDDO
      DUMX=DUMX**2/DUMW/2.0D0/C1
      DUMY=DUMY**2/DUMW/2.0D0/C1
      DUMZ=DUMZ**2/DUMW/2.0D0/C1
      ETCM=DUMX+DUMY+DUMZ
C

C.....modified by Bin, 3/19/2015
      IF (NSURF.GT.0) THEN
         RDMASS=WTA(NPATH)
         QCMB(1)=QCMA(1)-QCMA(3)*VCMA(1)/VCMA(3)
         QCMB(2)=QCMA(2)-QCMA(3)*VCMA(2)/VCMA(3)
         QCMB(3)=0.0D0
         VCMB=0.0D0
      ELSE
         RDMASS=WTA(NPATH)*WTB(NPATH)/(WTA(NPATH)+WTB(NPATH))
      ENDIF
C.....end

      RCM=0.0D0
      DO I=1,3
         QR(I)=QCMA(I)-QCMB(I)
         VR(I)=VCMA(I)-VCMB(I)
         RCM=RCM+QR(I)**2
      ENDDO
      RCM=SQRT(RCM)
      VREL=0.0D0
      VRELSQ=0.0D0
      DO I=1,3
         VREL=VREL+VR(I)*QR(I)
         VRELSQ=VRELSQ+VR(I)**2
      ENDDO
      VREL=VREL/RCM
      EREL=RDMASS*VREL**2/2.0D0/C1
      ERELSQ=RDMASS*VRELSQ/2.0D0/C1
      OAM(1)=(QR(2)*VR(3)-QR(3)*VR(2))*RDMASS
      OAM(2)=(QR(3)*VR(1)-QR(1)*VR(3))*RDMASS
      OAM(3)=(QR(1)*VR(2)-QR(2)*VR(1))*RDMASS
      OAM(4)=SQRT(OAM(1)**2+OAM(2)**2+OAM(3)**2)
      BF=OAM(4)/RDMASS/SQRT(VRELSQ)
      DO I=1,4
         OAM(I)=OAM(I)/C7
      ENDDO

C.....modified by Bin, 3/19/2015
      DO I=1,NDG
         ANG(I)=0.0D0
      ENDDO
      IF (NSURF.GT.0) GOTO 70
c.....end
C
C         GAS PHASE STUDY: 
C
C         CALCULATE SCATTERING ANGLES
C
C         VI AND VF ARE THE INITIAL AND FINAL RELATIVE VELOCITY
C         LI AND LF ARE THE INITIAL AND FINAL ORBITAL ANG. MOM.
C         JAI AND JAF ARE THE INITIAL AND FINAL A ROTATIONAL ANG. MOM.
C         JBI AND JBF ARE THE INITIAL AND FINAL B ROTATIONAL ANG. MOM.
C
C             LF,JAF:  ANG(1)
C             LF,JBF:  ANG(2)
C             JAF,JBF: ANG(3)
C
      IF (NATOMA(NPATH).NE.1) THEN
         DUM=(OAM(1)*AMA(1)+OAM(2)*AMA(2)+OAM(3)*AMA(3))
     *       /OAM(4)/AMA(4)
         IF (DUM.GT.1.00D0) DUM=1.00D0
         IF (DUM.LT.-1.00D0) DUM=-1.00D0
         ANG(1)=ACOS(DUM)/C4
      ENDIF
      IF (NATOMB(NPATH).NE.1) THEN
         DUM=(OAM(1)*AMB(1)+OAM(2)*AMB(2)+OAM(3)*AMB(3))
     *         /OAM(4)/AMB(4)
         ANG(2)=ACOS(DUM)/C4
         IF (NATOMA(NPATH).NE.1) THEN
            DUM=(AMA(1)*AMB(1)+AMA(2)*AMB(2)
     *          +AMA(3)*AMB(3))/AMA(4)/AMB(4)
            ANG(3)=ACOS(DUM)/C4
         ENDIF
      ENDIF
C
C             VI,VF:   ANG(4)
C             LI,LF:   ANG(5)
C             LI,JAF:  ANG(6)
C             LI,JBF:  ANG(7)
C
      IF (NATOMB(1).NE.0) THEN
         VRELSQ=SQRT(VRELSQ)
         DUM=(VI(1)*VR(1)+VI(2)*VR(2)+VI(3)*VR(3))/VI(4)/VRELSQ
         VRELSQ=VRELSQ**2
         ANG(4)=ACOS(DUM)/C4
         IF (OAMI(4).GE.1.00D-05) THEN
            DUM=(OAMI(1)*OAM(1)+OAMI(2)*OAM(2)+OAMI(3)*OAM(3))
     *          /OAMI(4)/OAM(4)
            ANG(5)=ACOS(DUM)/C4
            IF (NATOMA(NPATH).NE.1) THEN
               DUM=(OAMI(1)*AMA(1)+OAMI(2)*AMA(2)+OAMI(3)*AMA(3))
     *             /OAMI(4)/AMA(4)
               ANG(6)=ACOS(DUM)/C4
            ENDIF
            IF (NATOMB(NPATH).NE.1) THEN
               DUM=(OAMI(1)*AMB(1)+OAMI(2)*AMB(2)+OAMI(3)*AMB(3))
     *              /OAMI(4)/AMB(4)
               ANG(7)=ACOS(DUM)/C4
            ENDIF
         ENDIF
C
C             JAI,LI:  ANG(8)
C             JAI,LF:  ANG(9)
C             JAI,JAF: ANG(10)
C             JAI,JBF: ANG(11)
C             JAI,JBI: ANG(12)
C
         IF (AMAI(4).GE.1.00D-05) THEN
            IF (OAMI(4).GE.1.00D-05) THEN
               DUM=(AMAI(1)*OAMI(1)+AMAI(2)*OAMI(2)+
     *              AMAI(3)*OAMI(3))/AMAI(4)/OAMI(4)
               ANG(8)=ACOS(DUM)/C4
            ENDIF
            DUM=(AMAI(1)*OAM(1)+AMAI(2)*OAM(2)+AMAI(3)*OAM(3))
     *           /AMAI(4)/OAM(4)
            ANG(9)=ACOS(DUM)/C4
            IF (NATOMA(NPATH).NE.1) THEN
               DUM=(AMAI(1)*AMA(1)+AMAI(2)*AMA(2)+AMAI(3)*AMA(3))
     *              /AMAI(4)/AMA(4)
               ANG(10)=ACOS(DUM)/C4
            ENDIF
            IF (NATOMB(NPATH).NE.1) THEN
               DUM=(AMAI(1)*AMB(1)+AMAI(2)*AMB(2)+AMAI(3)*AMB(3))
     *              /AMAI(4)/AMB(4)
               ANG(11)=ACOS(DUM)/C4
            ENDIF
            IF (AMBI(4).LT.1.00D-05) GOTO 70
            DUM=(AMAI(1)*AMBI(1)+AMAI(2)*AMBI(2)
     *          +AMAI(3)*AMBI(3))/AMAI(4)/AMBI(4)
            ANG(12)=ACOS(DUM)/C4
         ENDIF
C
C             JBI,LI:  ANG(13)
C             JBI,LF:  ANG(14)
C             JBI,JAF: ANG(15)
C             JBI,JBF: ANG(16)
C
         IF (AMBI(4).GE.1.00D-05) THEN
            IF (OAMI(4).GE.1.00D-05) THEN
               DUM=(AMBI(1)*OAMI(1)+AMBI(2)*OAMI(2)
     *             +AMBI(3)*OAMI(3))/AMBI(4)/OAMI(4)
               ANG(13)=ACOS(DUM)/C4
            ENDIF
            DUM=(AMBI(1)*OAM(1)+AMBI(2)*OAM(2)
     *           +AMBI(3)*OAM(3))/AMBI(4)/OAM(4)
            ANG(14)=ACOS(DUM)/C4
            IF (NATOMA(NPATH).NE.1) THEN
               DUM=(AMBI(1)*AMA(1)+AMBI(2)*AMA(2)
     *             +AMBI(3)*AMA(3))/AMBI(4)/AMA(4)
               ANG(15)=ACOS(DUM)/C4
            ENDIF
            IF (NATOMB(NPATH).NE.1) THEN
               DUM=(AMBI(1)*AMB(1)+AMBI(2)*AMB(2)
     *             +AMBI(3)*AMB(3))/AMBI(4)/AMB(4)
               ANG(16)=ACOS(DUM)/C4
            ENDIF
         ENDIF

C
C              SURFACE STUDY: 
C              ANG(17): SCATTERING ANLGE WITH RESPECT TO SURFACE NORMAL
C              ANG(18): AZIMUTHAL ANGLE
C
         IF (NSURF.GE.2) THEN
C
C           DUMX,DUMY,DUMZ ARE THE COORDINATE OF A UNIT VECTOR ORTHOGONAL
C           TO THE SURFACE PLANE THAT HAS BEEN DEFINED BY THE USER
C
            DUMX1=Q(3*NN2-2)-Q(3*NN1-2)
            DUMY1=Q(3*NN2-1)-Q(3*NN1-1)
            DUMZ1=Q(3*NN2)-Q(3*NN1)
            DUMX2=Q(3*NN3-2)-Q(3*NN1-2)
            DUMY2=Q(3*NN3-1)-Q(3*NN1-1)
            DUMZ2=Q(3*NN3)-Q(3*NN1)
            DUMX=DUMY1*DUMZ2-DUMZ1*DUMY2
            DUMY=DUMX2*DUMZ1-DUMX1*DUMZ2
            DUMZ=DUMX1*DUMY2-DUMY1*DUMX2
            PNORM=SQRT(DUMX*DUMX+DUMY*DUMY+DUMZ*DUMZ)
            DUMX=DUMX/PNORM
            DUMY=DUMY/PNORM
            DUMZ=DUMZ/PNORM
C
C           PROJECTION OF THE INITIAL AND FINAL VELOCITY OF THE FRAGMENT ON THE
C           SURFACE PLANE NORM (DUMX,DUMY,DUMZ) 
C
            DUNOI=DUMX*VI(1)+DUMY*VI(2)+DUMZ*VI(3)
            DUNOF=DUMX*VCMA(1)+DUMY*VCMA(2)+DUMZ*VCMA(3)
C
C           SCATTERING ANGLE WITH RESPECT TO THE SURFACE NORMAL
C
            VCMASQ=SQRT(VCMA(1)**2+VCMA(2)**2+VCMA(3)**2)
            ANG(17)=ACOS(DUNOF/VCMASQ)/C4
C
C           PROJECTION OF THE INITIAL AND FINAL VELOCITY OF THE FRAGMENT ON
C           THE SURFACE PLANE (i.e. REMOVE SURFACE NORMAL PART)
C
            DUPAIX=VI(1)-DUNOI*DUMX
            DUPAIY=VI(2)-DUNOI*DUMY
            DUPAIZ=VI(3)-DUNOI*DUMZ
            DUPAI=SQRT(DUPAIX**2+DUPAIY**2+DUPAIZ**2)
C
            DUPAFX=VCMA(1)-DUNOF*DUMX
            DUPAFY=VCMA(2)-DUNOF*DUMY
            DUPAFZ=VCMA(3)-DUNOF*DUMZ
            DUPAF=SQRT(DUPAFX**2+DUPAFY**2+DUPAFZ**2)
C
C           AZIMUTHAL ANGLE CHANGE RELATIVE TO INITIAL
C
            DUMAZY=DUPAIX*DUPAFX+DUPAIY*DUPAFY+DUPAIZ*DUPAFZ
            ANG(18)=ACOS(DUMAZY/(DUPAI*DUPAF))/C4
C
         ENDIF
      ENDIF

C.....modified by Bin, 3/19/2015
   70 CONTINUE
      IF (NSURF.GT.0) THEN
C
C              SURFACE STUDY:
C              ANG(17): SCATTERING ANLGE WITH RESPECT TO SURFACE NORMAL
C              ANG(18): AZIMUTHAL ANGLE
C
C
C           INCIDENT ANGLES ANG(17):THTA, ANG(18):CHI
C
            VISQ=SQRT(VI(1)**2+VI(2)**2+VI(3)**2)
            VISQ1=SQRT(VI(1)**2+VI(2)**2)
            ANG(17)=ACOS(-VI(3)/VISQ)/C4

            IF (VISQ1.LE.1D-6) THEN
               ANG(18)=0D0
            ELSE
               IF (-VI(2).GE.0D0) THEN
                  ANG(18)=ACOS(-VI(1)/VISQ1)/C4
               ELSE
                  ANG(18)=(2*PI-ACOS(-VI(1)/VISQ1))/C4
               ENDIF
            ENDIF
C
C           SCATTERING ANGLES ANG(19):THTAf, ANG(20):CHIf
C
            IF (NZDOWN.EQ.1) THEN
!              Surface oscillator model: use fragment B velocity
               VCMX = P(7)/W(3)
               VCMY = P(8)/W(3)
               VCMZ = P(9)/W(3)
            ELSE
               VCMX = VCMA(1)
               VCMY = VCMA(2)
               VCMZ = VCMA(3)
            END IF
            VCMASQ=SQRT(VCMX**2+VCMY**2+VCMZ**2)
            VCMASQ1=SQRT(VCMX**2+VCMY**2)
            ANG(19)=ACOS(VCMZ/VCMASQ)/C4
            IF (VCMASQ1.LE.1D-6) THEN
               ANG(20)=0D0
            ELSE
               IF (VCMY.GE.0D0) THEN
                  ANG(20)=ACOS(VCMX/VCMASQ1)/C4
               ELSE
                  ANG(20)=(2*PI-ACOS(VCMX/VCMASQ1))/C4
               ENDIF
            ENDIF
      ENDIF
C.....end

   25 NFINAL=1
C


  999 RETURN
      END
