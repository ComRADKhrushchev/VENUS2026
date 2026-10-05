      SUBROUTINE SELECT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      PARAMETER(NSTEP=5000)
C
C         SELECT INITIAL CONDITIONS FOR COORDINATES AND MOMENTA
C
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/TRANSB/TRANS,NREL
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/FORCES/NATOMS,I3N,NFC,NGLO
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/INTEGR/ATIME,NI,NID
      COMMON/HFIT/PSCALA,PSCALB,VZERO
      COMMON/PRFLAG/NFQP,NCOOR,NFR,NUMR,NFB,NUMB,NFA,NUMA,NFTAU,NUMTAU,
     *NFTET,NUMTET,NFDH,NUMDH,NFHT,NUMHT
      COMMON/FINALB/EROTA,EROTB,EA(3),EB(3),AMA(4),AMB(4),AN,AJ,BN,BJ,
     *OAM(4),EREL,ERELSQ,ETCM,BF,SDA,SDB,DELH(NDP),ANG(NDG),NFINAL
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/CHEMAC/WWA(NDA3),CA(NDA3YF,NDA3YF),AI(3),ENMTA,
     *AMPA(NDA3),WWB(NDA3),CB(NDA3YF,NDA3YF),BI(3),ENMTB,
     *AMPB(NDA3),SEREL,S,BMAX,TROTA,TROTB,ANQA(NDA3),ANQB(NDA3),
     *TVIBA,TVIBB,NROTA,NROTB,NOB
!--> ADDED BY BIN, 2014/6/18
      COMMON/NROTAEQ2/JROTA,KROTA,JROTB,KROTB
!--> END
      COMMON/DIATB/NNA,JA,NNB,JB
      COMMON/VECTB/VI(4),OAMI(4),AMAI(4),AMBI(4),ETAI,ERAI,ETBI,ERBI
      COMMON/ARRAYS/A(NDA3YF,NDA3YF),DA(NDA3),B(NDA3YF,NDA3YF),DB(NDA3)
      COMMON/EIGVL/EIG(NDA3YF)
      COMMON/SADDLE/EBAR,TBAR,EZERO,NBAR,IDIR,JDIR,IJDIR 
      COMMON/VANGB/PHASEA(5),PHASEB(5),IPHASA(5),IPHASB(5),NPHASA,NPHASB
      COMMON/VRSCAL/THERMOTEMP,NSEL,NSCALE,NEQUAL,NRGD
!...added by Bin 2016/10/21
      COMMON/QPSCAL/QTEMP(NDA3),PTEMP(NDA3)
!...end
      COMMON/THERMOBATH/NTHERMB,NRSCL,NTHMID(NDA)
      COMMON/TMPNJ/TRVA,AIA,TRVB,AIB
      COMMON/WNJ/WD1,WD2
      COMMON/SFEQUIL/INTEGRATOR,LLL,NITER
!-->..ADDED BY BIN 04/24/2014
      COMMON/HDIAG/HTMIN,HTMAX,ZASYM,DELM
!-->..END


C
C     ADDED FOR MICROCANONICAL SAMPLING FOR A CONICAL INTERSECTION
C     01/06/2011 KYOYEON
C
      COMMON/PAR1CI/SXCI,SYCI,GGCI,HHCI
      COMMON/PAR2CI/PXMAX,PXMIN,PYMAX,PYMIN,XXMAX,XXMIN,YYMAX,
     *YYMIN,PXLEN,PYLEN,XXLEN,YYLEN,CONEXX,CONEYY,CONEPX,CONEPY,N3NM8
C
      DIMENSION QMAXA(NDA3),QMINA(NDA3),PMAXA(NDA),QMAXB(NDA3),
     *QMINB(NDA3),PMAXB(NDA)
      DIMENSION ENLOW(NDA3),ENMOD(NDA3),SUM(0:NSTEP)
      DIMENSION QCM(3),VCM(3)
      SAVE NMBAR, NMA, NMB
      SAVE WWASTORE,WWBSTORE

C
   45 FORMAT(//5X,'SELECT:NORMAL MODE QUANTUM NUMBERS')
   46 FORMAT(5X,10F10.2)
  100 FORMAT(5X,'DIATOM A FREQUENCY =',F7.1,' CM-1, AND ENERGY =',
     *F7.2,' KCAL/MOL'/)
  106 FORMAT(5X,'SELECT:   JXA,JYA,JZA=',1P3D13.5,' H-BAR'/)
  117 FORMAT(5X,'REACTANT A')
  124 FORMAT(/5X,'SELECT:    EROTA =',F7.3,' KCAL/MOL'/
     *15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
  128 FORMAT(5X,'DIATOM B FREQUENCY =',F7.1,' CM-1, AND ENERGY =',
     *F7.2,' KCAL/MOL'/)
  136 FORMAT(5X,'SELECT:   JXB,JYB,JZB=',1P3D13.5,' H-BAR'/)
  158 FORMAT(//5X,'REACTANT B')
  162 FORMAT(/15X,'IMPACT PARAMETER=',F7.3,' A')
  163 FORMAT(/5X,'SELECT:   EROTB =',F7.3,' KCAL/MOL'/
     *15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
  196 FORMAT(/5X,'RELATIVE TRANSLATIONAL ENERGY SELECTED: ',F7.2,
     *' KCAL/MOL'/)
  198 FORMAT(/5X,'CHOSEN:   LX,LY,LZ =',1P3D13.5,' H-BAR')
  200 FORMAT(/5X,'CHOSEN:   EROT =',F7.3,' KCAL/MOL'/
     *15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
  206 FORMAT(/5X,'CHOSEN:   EROTA =',F7.3,' KCAL/MOL'/
     *15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
  207 FORMAT(/5X,'CHOSEN:   EROTB =',F7.3,' KCAL/MOL'/
     *15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
  208 FORMAT(/5X,'CHOSEN:   EVIBA =',F7.3,' KCAL/MOL'/)
  209 FORMAT(/5X,'CHOSEN:   EVIBB =',F7.3,' KCAL/MOL'/)
 
 201  FORMAT('THE ZERO POINT ENERGY CANNOT BE LARGER THAN',
     &           ' THE AVAILABLE ENERGY')
 202  FORMAT('ZERO POINT ENERGY ',F7.2)
 203  FORMAT('AVAILABLE ENERGY ',F7.2)
 205  FORMAT('INCREASE THE VALUE OF PARAMETER NSTEP TO AT LEAST ',I8)
C
C         NSELT=-1 PROGRAM DOES NORMAL MODE ANALYSIS
C         NSELT=0  Q'S AND P'S ARE READ IN
C         NSELT=1  PROGRAM FINDS MINIMUM ENERGY GEOMETRY
C                  (Q'S AND P'S ARE READ IN)
C         NSELT=2  CHOOSE INITIAL CONDITIONS FOR ONE OR TWO MOLECULES
C         NSELT=3  CHOOSE INITIAL CONDITIONS FROM POTENTIAL BARRIER
C               NBAR=1  MICROCANONICAL SAMPLING
C                       WITH FIXED REACTION COORDINATE ENERGY
C               NBAR=2  THERMAL SAMPLING
C               NBAR=3  MICROCANONICAL QUASICLASSICAL SAMPLING
C                       INCLUDING REACTION COORDINATE ENERGY
C                   ENMTA IS THE TOTAL AVAILABLE ENERGY
C                         + ZERO POINT ENERGY AT TS.
C
C         NACT=1  ACTIVATE WITH ORTHANT SAMPLING
C         NACT=2  ACTIVATE WITH MICROCANONICAL NORMAL MODE SAMPLING
C         NACT=3  ACTIVATE WITH NORMAL MODE SAMPLING
C         NACT=4  ACTIVATE WITH LOCAL MODE EXCITATION
C         NACT=5  ACTIVATE WITH BOLTZMANN VIBRATIONAL DISTRIBUTION
C         NACT=6   SAME AS NBAR=3 (CHANGES NACT TO 6 WHEN NBAR=3)
C         NACT=7   MOLECULAR DYNAMICS SAMPLING BY RESCALING VELOCITIES
C   RESCALE THE VELOCITIES OF THE SYSTEM ACCORDING TO GIVEN TEMP.
C   REFERENCE:  "MOLECULAR DYNAMICS SIMULATION" BY JIM HAILE  P.458
C   ONLY FRAGMENT B CAN HAVE THIS WHEN IT IS SURFACE.
C
      IF (NSELT.NE.2.AND.NSELT.NE.3) THEN
         READ(5,*)(Q(I),I=1,I3N)
         IF(NSELT.EQ.0) THEN
           READ(5,*)(P(I),I=1,I3N)
         ELSE
           DO I=1,I3N
             P(I)=0.D0
           ENDDO
         ENDIF
         CALL DVDQ    
         CALL ENERGY
         RETURN
      ENDIF
C
C         IN THE MAIN PROGRAM THE COORDINATES Q ARE SET EQUAL TO QZ
C         AND THE MOMENTA P ARE SET EQUAL TO ZERO (FOR NSELT=2).
C
C         INITIAL CONDITIONS FOR REACTANT A
C
      IF (NATOMA(1).LE.1) THEN
C
C             IF A IS AN ATOM ZERO ITS Q, P ELEMENTS
C
         DO K=1,3
            J=3*LA(1,1)-3+K
            Q(K)=0.0D0
            P(K)=0.0D0
            QQ(K)=0.0D0
            PP(K)=0.0D0
         ENDDO
         GOTO 126
      ENDIF

!-->    MODIFIED BY BIN, 4/24/2014
      IF (NSURF.GT.0) THEN      !-->IN THIS CASE, ONLY REACTANT A IS MOVING, SURFACE IS FIXED AS REACTANT B 
         J=NATOMA(1)
         DO I=1,J
            K3=3*I
            Q(K3)=QZA(1,K3)+ZASYM
         ENDDO
      ELSE 
!-->    END
C
C             DISPLACE REACTANT B BY 1000.0 ANGSTROMS
C
      WRITE(6,117)
      N=NATOMB(1)
      IF (N.NE.0) THEN
         DO I=1,N
            J3=3*LB(1,I)
            K3=3*I
!-->    MODIFIED BY BIN, 4/24/2014
C            Q(J3)=QZB(1,K3)+1000.0D0
            Q(J3)=QZB(1,K3)+ZASYM
!-->    END BIN'S CHANGES
         ENDDO
      ENDIF

!-->    MODIFIED BY BIN, 1/16/2014
      ENDIF
!-->    END

C
C             ENERGY REFERENCE FOR SEPARATED REACTANTS
C
C      CALL DVDQ        by bin, 2016/10/02
      CALL ENERGY
!-->    MODIFIED BY BIN, 1/16/2014
C      WRITE(*,*)
C      WRITE(*,*)'POTENTIAL CHECK=',V,'KCAL/MOL'
C      WRITE(*,*)
!-->    END
      DH=V

C
C             SELECT INITIAL Q'S AND P'S FOR REACTANT A
C
      IF (NATOMA(1).GT.2.OR.NLINA.EQ.0) GOTO 120
      L(1)=LA(1,1)
      L(2)=LA(1,2)
C
C          DIATOM A IS TREATED SEMICLASSICALLY
C
      IF (NSFLAG.NE.1) THEN
         IF(NTZ.EQ.1)THEN
          N=NATOMA(1)
!-->    MODIFIED BY BIN, 3/31/2015
          WRITE(26,*)'NORMAL MODES FOR FRAGMENT A IN PATH ',1
          WRITE(26,*)
          CALL NMODE(N,0)
!-->    END
          DUM=EIG(6)
         ENDIF

         IF(TRVA.GE.0.0D0) THEN

            IF (NACTA.EQ.0) THEN
!!!   adapted by Bin, 2016/10/10, for maxwell-boltzmann sampling
            DESKET=SQRT(0.00198717D0*TRVA*C1) 
            WT=WTA(1)
            DO I = 1,NATOMA(1)
               J3=3*LA(1,I)
               J2=J3-1
               J1=J2-1
               J=LA(1,I)
               P(J1)=GASDEV(ISEED)*DESKET*SQRT(W(J))
               P(J2)=GASDEV(ISEED)*DESKET*SQRT(W(J))
               P(J3)=GASDEV(ISEED)*DESKET*SQRT(W(J))
            ENDDO

            CALL CENMAS(WT,QCM,VCM,NATOMA(1))
C
            DO I=1,NATOMA(1)
               J=3*LA(1,I)
               DO K=0,2
                  Q(J-K)=QQ(J-K)      !Coordinates shifted to center of mass frame
                  PP(J-K)=P(J-K)      !This is important because we actually do not want to remove the center of mass velocity
               ENDDO
            ENDDO
    
            CALL ROTATE(NATOMA(1))

!!!   remove velocities along z axis which are replaced by the
!translational energy given later 
            DO I=1,NATOMA(1)
               J=3*LA(1,I)
               P(J)=0D0
               PP(J)=0D0
            ENDDO

            GOTO 126

            ENDIF
!!!   end Bin

           WRITE(6,*)'N AND J CHOSEN BASED ON TEMPERATURE'
           WRITE(6,*)'TRVA = ',TRVA
C          CALCULATE N
           IF(NTZ.EQ.1)THEN
            WWA(1)=EIG(6)*C6
            WWASTORE=WWA(1)
           ELSE
            WWA(1)=WWASTORE
            DUM=WWA(1)/C6  
           ENDIF
           NMBAR=1
!-->    MODIFIED BY BIN, 5/14/2014
           CALL THRMAN(WWA,ANNA,TRVA,NMBAR)
           NNA=NINT(ANNA)
!-->    END
           WRITE(6,*)'NNA = ',NNA
C          CALCULATE J
           WD1=W(L(1))
           WD2=W(L(2))
           CALL PROBJ(TRVA,AIA,ISEED,JA)
           WRITE(6,*)'JA = ',JA
         ELSE
           WRITE(6,*)'N AND J USED AS INPUT'
         ENDIF

         ENJA=(DBLE(NNA)+0.5D0)*DUM*0.0028591*C1
         CALL INITEBK(NNA,JA,RMINA,RMAXA,DH,RMASSA,ENJA,PTESTA,ALA)
         SDUM=ENJA/C1
         WRITE(6,100)DUM,SDUM
      ENDIF
C
C             SELECT INITIAL RELATIVE COORDINATE AND MOMENTUM.
C
      DUM=ALA**2/2.0D0/RMASSA
  102 RAND=RAND0(ISEED)
      R=RMINA+(RMAXA-RMINA)*RAND
      Q(3*L(1))=-0.5D0*R
      Q(3*L(1)-1)=0.0D0
      Q(3*L(1)-2)=0.0D0
      Q(3*L(2))= 0.5D0*R
      Q(3*L(2)-1)=0.0D0
      Q(3*L(2)-2)=0.0D0
!-->    MODIFIED BY BIN, 4/24/2014
!-->    SHIFT Z COORINDATES UP IN ORDER TO CALCULATE THE POTENTIAL ENERGY
      IF (NSURF.GT.0) THEN
         J=NATOMA(1)
         DO I=1,J
            K3=3*I
            Q(K3)=Q(K3)+ZASYM
         ENDDO
      ENDIF
!-->    END BIN'S CHANGES
!      CALL DVDQ        commented by bin
      CALL ENERGY
!-->    MODIFIED BY BIN, 4/24/2014
!-->    SHIFT Z COORINDATES BACK AFTER THE POTENTIAL ENERGY
      IF (NSURF.GT.0) THEN
      J=NATOMA(1)
      DO I=1,J
         K3=3*I
         Q(K3)=Q(K3)-ZASYM
      ENDDO
      ENDIF
!-->    END BIN'S CHANGES
      VDUM=(V-DH)*C1
      SUMM=ENJA-DUM/R**2-VDUM
      IF (SUMM.LE.0.0D0) THEN
         SUMM=0.0D0
         PR=0.0D0
      ELSE
         PR=SQRT(2.0D0*RMASSA*SUMM)
         SDUM=PTESTA/PR
         RAND=RAND0(ISEED)
         IF (SDUM.LT.RAND) GOTO 102
      ENDIF
      RAND=RAND0(ISEED)
      IF (RAND.LT.0.5D0) PR=-PR
C
C             CHOOSE INITIAL CARTESIAN COORDINATES AND MOMENTA, AND
C             ANGULAR MOMENTUM.  DIATOM LIES ALONG THE X-AXIS.
C             THEN RANDOMLY ROTATE THE CARTESIAN COORDINATES AND
C             MOMENTA IN THE CENTER OF MASS FRAME.
C
      CALL HOMOQP(R,PR,ALA,AMA,RMASSA,AI)
      AMAI(1)=AMA(1)/C7
      AMAI(2)=AMA(2)/C7
      AMAI(3)=AMA(3)/C7
      AMAI(4)=AMA(4)/C7
      CALL ROTN(AMA,EROTA,2)
      WRITE(6,106)AMAI(1),AMAI(2),AMAI(3)
      GOTO 126
C
  120 CONTINUE
C
C             FRAGMENT A IS A POLYATOMIC
C
       IF (NACTA.EQ.0) THEN
!!!   adapted by Bin, 2016/10/10, for maxwell-boltzmann sampling
          DESKET=SQRT(0.00198717D0*TVIBA*C1) 
          WT=WTA(1)
          DO I = 1,NATOMA(1)
             L(I)=LA(1,I)
             J3=3*LA(1,I)
             J2=J3-1
             J1=J2-1
             J=LA(1,I)
             P(J1)=GASDEV(ISEED)*DESKET*SQRT(W(J))
             P(J2)=GASDEV(ISEED)*DESKET*SQRT(W(J))
             P(J3)=GASDEV(ISEED)*DESKET*SQRT(W(J))
          ENDDO

            CALL CENMAS(WT,QCM,VCM,NATOMA(1))
C
            DO I=1,NATOMA(1)
               J=3*LA(1,I)
               DO K=0,2
                  Q(J-K)=QQ(J-K)      !Coordinates shifted to center of mass frame
                  PP(J-K)=P(J-K)      !This is important because we actually do not want to remove the center of mass velocity
               ENDDO
            ENDDO
    
          CALL ROTATE(NATOMA(1))

!!!   remove velocities along z axis which are replaced by the
!translational energy given later 
            DO I=1,NATOMA(1)
               J=3*LA(1,I)
               P(J)=0D0
               PP(J)=0D0
            ENDDO
 
          GOTO 126

          ENDIF
!!!   end Bin


      IF (NACTA.EQ.1) GOTO 123
C
C             CALCULATE NORMAL MODE EIGENVALUES AND EIGENVECTORS
C
      IF(NSFLAG.EQ.1)THEN
        IF(NACTA.EQ.5.OR.NACTA.EQ.8.OR.NACTA.EQ.9)THEN
          GOTO 113
        ELSE
          GOTO 123
        ENDIF
      ENDIF
C
C     IF (NSFLAG.EQ.1.AND.NACTA.NE.5) GOTO 123
C     IF (NSFLAG.EQ.1.AND.NACTA.EQ.5.OR.
C                         NACTA.EQ.8.OR.NACTA.EQ.9) GOTO 113
C
      N=NATOMA(1)
      K=3*N
      M=6-NLINA
      NMA=K-M
C
C             TRANSITION STATE ONLY HAS 3N-7 NORMAL MODES FOR A 
C             NONLINEAR MOLECULE, AND 3N-6 MODES FOR A LINEAR ONE
C
      IF (NSELT.EQ.3) THEN
         NMBAR=NMA-1
         IBARR=1
      ELSE
         NMBAR=NMA
         IBARR=0
      ENDIF
C
!-->    MODIFIED BY BIN, 3/31/2015
      WRITE(26,*)'NORMAL MODES FOR FRAGMENT A IN PATH ',1
      WRITE(26,*)
      CALL NMODE(N,0)
!-->    END
      DO I=1,NMBAR
         WWA(I)=EIG(I+M+IBARR)*C6
         DO J=1,K
            CA(J,I)=A(I+M+IBARR,J)
         ENDDO
      ENDDO
      IF (NSELT.EQ.3) THEN
         WWA(NMA)=EIG(1)*C6
         DO J=1,K
            CA(J,NMA)=A(1,J)
         ENDDO
      ENDIF
C
C     THIS PART IS FOR MICROCANONICAL SAMPLING.
C     CHECK THAT THE TOTAL VIBRATIONAL ENERGY DOES NOT EXCEED THE
C     MAXIMUM TOTAL ENERGY THAT IS USED IN BEYER-SWINEHARDT COUNTING
C
      IF(NACTA.EQ.6) THEN
         NBAR=1                   ! FIXED ENERGY IN REACTION COORDINATE
         DUM=ENMTA*CAL2CM
         ZPE=0.0D0
         DO I=1,NMBAR
            ENMOD(I)=WWA(I)/C6
            ZPE=ZPE+WWA(I)/C6
            ENLOW(I)=ENMOD(I)*0.5D0*0.0028591*C1
         ENDDO
         ZPE=ZPE*0.5D0
         DUM=DUM-ZPE
         STEP=ENMOD(1)*0.1D0      ! ENMOD(1) IS THE LOWEST FREQUENCY
         IF(DUM/STEP.GE.DBLE(NSTEP)) THEN
            WRITE(6,205)IDNINT(DUM/STEP)+1
            STOP
         ENDIF
         IF(ZPE/CAL2CM.GT.ENMTA) THEN
            WRITE(6,201)
            WRITE(6,202)ZPE*0.0028591
            WRITE(6,203)ENMTA
            STOP
         ENDIF
         DUMNAC6=DUM
      ENDIF
C
      IF ((NACTA.EQ.2).OR.(NACTA.EQ.6)) GOTO 123
C
C             CHOOSE NORMAL MODE QUANTUM NUMBERS FROM A THERMAL
C             DISTRIBUTION IF NACT=5
C
  113 CONTINUE
      IF (NACTA.EQ.5) THEN
         CALL THRMAN(WWA,ANQA,TVIBA,NMBAR)
         WRITE(6,45)
         WRITE(6,46)(ANQA(I),I=1,NMBAR)
         WRITE(9,46)(ANQA(I),I=1,NMBAR)
      ENDIF      
C
C             SELECT CLASSICAL 2D CONE ENERGY AND 
C             QUANTUM ENERGY FOR REMAINING 3N-8 NORMAL MODES
C
      IF (NACTA.EQ.9) CALL MICROCI(ECONE,EQNM)
C
C             SELECT QUANTUM NUMBERS FOR 
C             A QUANTUM MICROCANONICAL ENSEMBLE OF NORMAL MODES
C
      IF (NACTA.EQ.8.OR.NACTA.EQ.9) THEN
         IF (NACTA.EQ.8) THEN
            EQNM=ENMTA
            CALL QMMICRO(WWA,ANQA,EQNM,NMA)
            WRITE(6,45)
            WRITE(6,46)(ANQA(I),I=1,NMA)
            WRITE(9,46)(ANQA(I),I=1,NMA)
         ELSEIF (NACTA.EQ.9) THEN
            CALL QMMICRO(WWA,ANQA,EQNM,NM3N8)
            WRITE(6,45)
            WRITE(6,46)(ANQA(I),I=1,NM3N8)
            WRITE(9,46)(ANQA(I),I=1,NM3N8)
         ENDIF
      ENDIF

C
C             CALCULATE NORMAL MODE ENERGIES AND AMPLITUDES 
C             FOR NACTA=3,4,5,8,9.
C             ENMTA CAN NOT CHANGE FOR NACTA=8 OR 9
C
      ENMDUM=0.0D0
      DO I=1,NMBAR
         DUM=(ANQA(I)+0.5D0)*WWA(I)/(C6*CAL2CM)
         ENMDUM=ENMDUM+DUM
         DUM=DUM*C1
         AMPA(I)=SQRT(2.0D0*DUM)/WWA(I)
      ENDDO
      IF (NACTA.EQ.3.OR.NACTA.EQ.4.OR.NACTA.EQ.5) ENMTA=ENMDUM
      IF (NACTA.EQ.8.OR.NACTA.EQ.9) ENM89=ENMDUM
C
C             CHOOSE THE ANGULAR MOMENTUM VECTOR
C
  123 CONTINUE

      CALL ROTEN(AMA,AI,TROTA,EROTA,NROTA,NLINA,JROTA,KROTA)

C
C             SAVE INITIAL ROTATIONAL ANGULAR MOMENTUM
C
      AMAI(1)=AMA(1)/C7
      AMAI(2)=AMA(2)/C7
      AMAI(3)=AMA(3)/C7
      AMAI(4)=SQRT(AMAI(1)**2+AMAI(2)**2+AMAI(3)**2)
      WRITE(6,124)EROTA,AMAI(1),AMAI(2),AMAI(3)
C
C             CALCULATE THE TOTAL ENERGY
C
      ETAI=EROTA+ENMTA
      IF (NACTA.EQ.8.OR.NACTA.EQ.9) ETAI=EROTA+ENM89
C
C             CHOOSE THE INITIAL COORDINATES AND MOMENTA
C
      N=NATOMA(1)
      DO I=1,N
         L(I)=LA(1,I)
         J3=3*L(I)
         J2=J3-1
         J1=J2-1
         K3=3*I
         K2=K3-1
         K1=K2-1
         QZ(J1)=QZA(1,K1)
         QZ(J2)=QZA(1,K2)
         QZ(J3)=QZA(1,K3)
      ENDDO
      DUM1=WTA(1)
C
C           SELECT CARTESIAN COORDINATES AND MOMENTA 
C           USING ORTHANT SAMPLING
C
      IF (NACTA.EQ.1) THEN
         CALL ORTHAN(AMA,DUM1,ENMTA,ETAI,QMAXA,QMINA,PMAXA,
     *               PSCALA,ERAI,N)
         GOTO 126
      ENDIF
C
C           CALCULATE NORMAL MODE ENERGIES AND AMPLITUDES FOR A
C           CLASSICAL MICROCANONICAL ENSEMBLE OF NORMAL MODES 
C
      IF (NACTA.EQ.2) THEN
         DUM=ENMTA*C1
         NN=NMBAR-1
         DO I=1,NN
            RAND=RAND0(ISEED)
            SDUM=1.0D0/DBLE(NMBAR-I)
            SDUM=DUM*(1.0D0-RAND**SDUM)
            DUM=DUM-SDUM
            AMPA(I)=SQRT(2.0D0*SDUM)/WWA(I)
         ENDDO
         AMPA(NMBAR)=SQRT(2.0D0*DUM)/WWA(NMBAR)
      ENDIF
C
      IF (NACTA.EQ.6) THEN
         DUM=DUMNAC6     ! TOTAL ENERGY ABOVE ZPE IN 1/CM
         DO II=1,NMBAR   ! LOWEST TO HIGHEST FREQUENCY FOR EFFICIENCY
            IMAXST=IDINT(DUM/ENMOD(II)) ! MAXIMUM QUANTUM NO.
            IF(II.NE.NMBAR) THEN
               STEP=ENMOD(II+1)*0.1D0
C
C     THIS IS THE BEYER-SWINEHARDT ALGORITHM
C
               IEND=IDNINT(DUM/STEP)
               DO J=0,IEND
                  SUM(J)=1.0D0
               ENDDO
               DO J=II+1,NMBAR
                  ISTART=IDNINT(ENMOD(J)/STEP)
                  DO K=ISTART,IEND
                     SUM(K)=SUM(K)+SUM(K-ISTART)
                  ENDDO
               ENDDO
C     
C     
               PROB=0.0D0
               DO J=IMAXST,0,-1
                  ENDUM=DUM-DBLE(J)*ENMOD(II) ! ENERGY IN REST OF MOLECULE
                  IDDUM=IDNINT(ENDUM/STEP)
                  PROB=PROB+SUM(IDDUM)
               ENDDO
            ELSE
               STEP=ENMOD(II)
               PROB=0.0D0
               DO J=IMAXST,0,-1
                  ENDUM=DUM-DBLE(J)*ENMOD(II)
                  IDDUM=IDNINT(ENDUM/STEP)
                  SUM(IDDUM)=1.0D0
                  PROB=PROB+SUM(IDDUM)
               ENDDO
            ENDIF
            RAND=RAND0(ISEED)*PROB
            IF(RAND.GT.PROB)RAND=PROB
            PROB=0.0D0
            DO J=IMAXST,0,-1
               ENDUM=DUM-DBLE(J)*ENMOD(II)
               IDDUM=IDNINT(ENDUM/STEP)
               PROB=SUM(IDDUM)+PROB
               IF(RAND.LE.PROB) THEN
                  SDUM=ENMOD(II)*DBLE(J)
                  GOTO 1115
               ENDIF
            ENDDO
 1115       CONTINUE
            DUM=DUM-SDUM
            SDUM=SDUM*0.0028591*C1
            AMPA(II)=SQRT(2.0D0*(SDUM+ENLOW(II)))/WWA(II)
         ENDDO
C             GIVE THE REMAINDER OF THE ENERGY TO THE REACTION COORDINATE
         EBAR=DUM*0.0028591
         ETAI=ETAI-EBAR
      ENDIF
C
      CALL INITQP(WWA,AMPA,CA,AMA,DUM1,ETAI,EROTA,AI,ERAI,PHASEA,N,
     *NMBAR,NPHASA,IPHASA)
C
C             INITIAL CONDITION FOR REACTION COORDINATE
C
      IF (NSELT.EQ.3) CALL BAREXC(DUM1,CA,AMA,ERAI,N,NMA)
C
C             INITIAL CONDITIONS FOR REACTANT B
C
  126 CONTINUE

!-->    MODIFIED BY BIN, 1/16/2014
      IF (NSURF.EQ.3) GOTO 160 
!-->    END

      IF (NATOMB(1).EQ.0) THEN
         IF (NATOMA(1).NE.2) GOTO 999
         CALL DVDQ
         CALL ENERGY
         GOTO 999
      ENDIF
      IF (NATOMB(1).GT.1) GOTO 129
C
C             IF FRAGMENT B IS AN ATOM ZERO ITS Q ELEMENTS
C
      J3=3*LB(1,1)
      J2=J3-1
      J1=J2-1
      Q(J1)=0.0D0
      Q(J2)=0.0D0
      Q(J3)=0.0D0
      GOTO 160
C
C             SELECT INITIAL Q'S AND P'S FOR REACTANT B.
C
C             SAVE Q'S AND P'S OF A IN TEMPORARY STORAGE. (IN FACT,
C             THIS IS UNNECESSORY AS WE CAN STORE THEM FROM QQ,PP)
C             SET P ARRAY TO ZERO AND EQUATE THE Q AND QZ ARRAYS FOR A
C
  129 CONTINUE
      WRITE(6,158)
      N=NATOMA(1)
      DO I=1,N
         J3=3*LA(1,I)
         J2=J3-1
         J1=J2-1
         K3=3*I
         K2=K3-1
         K1=K2-1
         Q(J1)=QZA(1,K1)
         Q(J2)=QZA(1,K2)
!-->    MODIFIED BY BIN, 4/24/2014
C         Q(J3)=QZA(1,K3)+1000.0D0
         Q(J3)=QZA(1,K3)+ZASYM
!-->    END BIN'S CHANGES
         P(J1)=0.0D0
         P(J2)=0.0D0
         P(J3)=0.0D0
      ENDDO
C
C             RESET Q ARRAY TO QZ ARRAY FOR B
C
      N=NATOMB(1)
      DO I=1,N
         J3=3*LB(1,I)
         K3=3*I
         Q(J3)=QZB(1,K3)
      ENDDO
C
C    ADDED BY BIN, 2016/9/30
      IF (NATOMB(1).GT.2.AND.NACTB.EQ.7) GOTO 111
C    END

      IF (NATOMB(1).GT.2.OR.NLINB.EQ.0) GOTO 150
      L(1)=LB(1,1)
      L(2)=LB(1,2)
C
C             DIATOM B IS TREATED SEMICLASSICALLY
C
      IF (NSFLAG.NE.1) THEN
         IF(NTZ.EQ.1)THEN
          I=NATOMA(1)
          N=NATOMB(1)
!       MODIFIED BY BIN, 3/31/2015
          WRITE(26,*)'NORMAL MODES FOR FRAGMENT B IN PATH ',1
          WRITE(26,*)
          CALL NMODE(N,I)
!       END
          DUM=EIG(6)
         ENDIF

         IF(TRVB.GE.0.0D0)THEN
           WRITE(6,*)'N AND J CHOSEN BASED ON TEMPERATURE'
           WRITE(6,*)'TRVB = ',TRVB
C          CALCULATE N
           IF(NTZ.EQ.1)THEN
            WWB(1)=EIG(6)*C6
            WWBSTORE=WWB(1)
           ELSE
            WWB(1)=WWBSTORE
            DUM=WWB(1)/C6
           ENDIF
           NMBAR=1
!-->    MODIFIED BY BIN, 5/14/2014
           CALL THRMAN(WWB,ANNB,TRVB,NMBAR)
           NNB=NINT(ANNB)
!-->    END
           WRITE(6,*)'NNB = ',NNB
C          CALCULATE J
           WD1=W(L(1))
           WD2=W(L(2))
           CALL PROBJ(TRVB,AIB,ISEED,JB)
           WRITE(6,*)'JB = ',JB
         ELSE
           WRITE(6,*)'N AND J USED AS INPUT'
         ENDIF

         ENJB=(DBLE(NNB)+0.5D0)*DUM*0.0028591*C1
		 DH = DELH(1) !CZZ
         CALL INITEBK(NNB,JB,RMINB,RMAXB,DH,RMASSB,ENJB,PTESTB,ALB)
         SDUM=ENJB/C1
         WRITE(6,128)DUM,SDUM
      ENDIF
C
C             SELECT INITIAL RELATIVE COORDINATE AND MOMENTUM
C
      DUM=ALB**2/2.0D0/RMASSB
  132 RAND=RAND0(ISEED)
      R=RMINB+(RMAXB-RMINB)*RAND
      Q(3*L(1))=-0.5D0*R
      Q(3*L(1)-1)=0.0D0
      Q(3*L(1)-2)=0.0D0
      Q(3*L(2))= 0.5D0*R
      Q(3*L(2)-1)=0.0D0
      Q(3*L(2)-2)=0.0D0
      CALL DVDQ
      CALL ENERGY
      VDUM=(V-DH)*C1
      SUMM=ENJB-DUM/R**2-VDUM
      IF (SUMM.LE.0.0D0) THEN
         SUMM=0.0D0
         PR=0.0D0
      ELSE
         PR=SQRT(2.0D0*RMASSB*SUMM)
         SDUM=PTESTB/PR
         RAND=RAND0(ISEED)
         IF (SDUM.LT.RAND) GOTO 132
      ENDIF
      RAND=RAND0(ISEED)
      IF (RAND.LT.0.5D0) PR=-PR
C
C             CHOOSE INITIAL CARTESIAN COORDINATES AND MOMENTA, AND
C             ANGULAR MOMENTUM.  DIATOM LIES ALONG THE X-AXIS.
C             THEN RANDOMLY ROTATE THE CARTESIAN COORDINATES AND
C             MOMENTA IN THE CENTER OF MASS FRAME.
C
  142 CONTINUE
      CALL HOMOQP(R,PR,ALB,AMB,RMASSB,BI)
      AMBI(1)=AMB(1)/C7
      AMBI(2)=AMB(2)/C7
      AMBI(3)=AMB(3)/C7
      AMBI(4)=AMB(4)/C7
      CALL ROTN(AMB,EROTB,2)
      WRITE(6,136)AMBI(1),AMBI(2),AMBI(3)
      GOTO 155
C
  150 CONTINUE
C
C             FRAGMENT B IS A POLYATOMIC
C
      IF (NACTB.EQ.1) GOTO 153
      IF (NACTB.EQ.7) GOTO 111
C
C             CALCULATE NORMAL MODE EIGENVALUES AND EIGENVECTORS
C

      IF (NSFLAG.EQ.1) THEN
         IF(NACTB.NE.5.AND.NACTB.NE.8) GOTO 153
      END IF

      IF (NSFLAG.NE.1.OR.NACTB.NE.5) THEN
         N=NATOMB(1)
         K=3*N
         M=6-NLINB
         NMB=K-M
         I=NATOMA(1)
!       MODIFIED BY BIN, 3/31/2015
         WRITE(26,*)'NORMAL MODES FOR FRAGMENT B IN PATH ',1
         WRITE(26,*)
         CALL NMODE(N,I)
!       END
         DO I=1,NMB
            WWB(I)=EIG(I+M)*C6
            DO J=1,K
               CB(J,I)=A(I+M,J)
            ENDDO
         ENDDO
      ENDIF
C
      IF (NACTB.EQ.2) GOTO 153
C
C             CHOOSE NORMAL MODE QUANTUM NUMBERS FROM A THERMAL
C             DISTRIBUTION IF NACT=5
C
      IF (NACTB.EQ.5) THEN
         CALL THRMAN(WWB,ANQB,TVIBB,NMB)
         WRITE(6,45)
         WRITE(6,46)(ANQB(I),I=1,NMB)
         WRITE(9,46)(ANQB(I),I=1,NMB)
      ENDIF

C
C             SELECT QUANTUM NUMBERS FOR 
C             A QUANTUM MICROCANONICAL ENSEMBLE OF NORMAL MODES
C
      IF (NACTB.EQ.8) THEN
         EQNM=ENMTB
         CALL QMMICRO(WWB,ANQB,EQNM,NMB)
         WRITE(6,45)
         WRITE(6,46)(ANQB(I),I=1,NMB)
         WRITE(9,46)(ANQB(I),I=1,NMB)
      ENDIF
C
C             CALCULATE NORMAL MODE ENERGIES AND AMPLITIDES 
C             FOR NACTB=3,4,5,8.
C             ENMTB CAN NOT CHANGE FOR NACTB=8 
C
      ENMDUM=0.0D0
      DO I=1,NMB
         DUM=(ANQB(I)+0.5D0)*WWB(I)/(C6*CAL2CM)
         ENMDUM=ENMDUM+DUM
         DUM=DUM*C1
         AMPB(I)=SQRT(2.0D0*DUM)/WWB(I)
      ENDDO
      IF (NACTB.EQ.3.OR.NACTB.EQ.4.OR.NACTB.EQ.5) ENMTB=ENMDUM
      IF (NACTB.EQ.8) ENM89=ENMDUM
C
C             CHOOSE THE ANGULAR MOMUNTUM VECTOR
C
  153 CONTINUE
      CALL ROTEN(AMB,BI,TROTB,EROTB,NROTB,NLINB,JROTB,KROTB)
C
C             SAVE INITIAL ROTATIONAL ANGULAR MOMENTUM
C
      AMBI(1)=AMB(1)/C7
      AMBI(2)=AMB(2)/C7
      AMBI(3)=AMB(3)/C7
      AMBI(4)=SQRT(AMBI(1)**2+AMBI(2)**2+AMBI(3)**2)
      WRITE(6,163)EROTB,AMBI(1),AMBI(2),AMBI(3)
C
C             CALCULATE THE TOTAL ENERGY
C
      ETBI=EROTB+ENMTB
      IF (NACTB.EQ.8) ETBI=EROTB+ENM89
C
C             CHOOSE THE INITIAL COORDINATES AND MOMENTA
C
  111 CONTINUE
      N=NATOMB(1)
      DO I=1,N
         L(I)=LB(1,I)
         J3=3*L(I)
         J2=J3-1
         J1=J2-1
         K3=3*I
         K2=K3-1
         K1=K2-1
         QZ(J1)=QZB(1,K1)
         QZ(J2)=QZB(1,K2)
         QZ(J3)=QZB(1,K3)
      ENDDO
      IF(NACTB.EQ.7) GOTO 155
      DUM1=WTB(1)
C
C       SELECT CARTESIAN COORDINATES AND MOMENTA USING
C       ORTHANT SAMPLING
C
      IF (NACTB.EQ.1) THEN
         CALL ORTHAN(AMB,DUM1,ENMTB,ETBI,QMAXB,QMINB,PMAXB,
     *               PSCALB,ERBI,N)
         GOTO 155
      ENDIF
C
C             CALCULATE NORMAL MODE ENERGIES AND AMPLITUDES FOR A
C             CLASSICAL MICROCANONICAL ENSEMBLE OF NORMAL MODES 
C
      IF (NACTB.EQ.2) THEN
         DUM=ENMTB*C1
         NN=NMB-1
         DO I=1,NN
            RAND=RAND0(ISEED)
            SDUM=1.0D0/DBLE(NMB-I)
            SDUM=DUM*(1.0D0-RAND**SDUM)
            DUM=DUM-SDUM
            AMPB(I)=SQRT(2.0D0*SDUM)/WWB(I)
         ENDDO
         AMPB(NMB)=SQRT(2.0D0*DUM)/WWB(NMB)
      ENDIF
      CALL INITQP(WWB,AMPB,CB,AMB,DUM1,ETBI,EROTB,BI,ERBI,PHASEB,N,
     *NMB,NPHASB,IPHASB)
C
C         RESTORE Q'S AND P'S FROM TEMPORARY STORAGE
C         (IN FACT WE CAN STORE FROM QQ, PP)
C
  155 CONTINUE
C--------------------------------------------------------------------
C   INITIALIZE A PHASE SPACE POINT BY RANDOMLY ASSIGNING VELOCITIES
C   AT EQUALIBRIUM COORDINATES

!...rewrite by Bin, 2016/10/11
      IF (NACTB.EQ.7) THEN
        TELEC=THERMOTEMP
C
C       STORE THE COORDINATES AND MOMENTA OF REACTANT A FIRST
C
        IF (THERMOTEMP.LE.0.01D0) THEN      
           N=NATOMB(1)-NRGD
           DO I = 1,N
              J3=3*LB(1,I)
              J2=J3-1
              J1=J2-1
              J=LB(1,I)
              P(J1)=0D0
              P(J2)=0D0
              P(J3)=0D0
           ENDDO
           GOTO 159
         ENDIF

         N=NATOMA(1)

         DO I=1,N
            J3=3*LA(1,I)
            J2=J3-1
            J1=J2-1
            QTEMP(J3)=Q(J3)
            QTEMP(J2)=Q(J2)
            QTEMP(J1)=Q(J1)
            P(J3) = 0.0D0
            P(J2) = 0.0D0
            P(J1) = 0.0D0
         END DO

         NCOORORG=NCOOR
         NCOOR=0
         N=NATOMB(1)-NRGD
C
C  TO RUN EQUILIBRATION, USE A TIME STEP OF 0.5 FS IS USUALLY SUFFICIENT
C
        TIMEORG=TIME
        ATIMEORG=ATIME
        TIME=0.05D0
        IF (INTEGRATOR.EQ.0) THEN
          ATIME=TIME/1440.0D0
        ELSE
          ATIME=TIME
        ENDIF
        IF (INTEGRATOR.EQ.1) THEN
           WRITE(*,*)'INTEGRATOR=1 IS NOT SUPPORTED IN NACTB=7'
           WRITE(*,*)'SEE SELECT.f FOR DETAILS'
           STOP
        ENDIF

        NC=0
C
C   RANDOMLY ASSIGNING VELOCITIES FIRSTLY THEN WILL BE RESCALED
C   TO THE TEMPERATURE WANTED. TWO STEPS
C      1. GENERATE RANDOM NUMBERS: GENERATE 3N RANDOM NUMBERS FROM A NORMAL 
C         DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE. 
C      2. GENERATE MOMENTA: MULTIPLY THESE RANDOM NUMBERS BY SQRT(MKT) TO OBTAIN 
C         THE X, Y, AND Z COMPONENTS OF MOMENTA FOR THE N PARTICLES
C   DESKET = SQRT(KB*T), KB - BOLTZMANN CONST. T - TEMPERATURE
C
        IF (NTZ.EQ.1) THEN      
           DESKET=SQRT(0.00198717D0*THERMOTEMP*C1) 
           DO I = 1,N
              J3=3*LB(1,I)
              J2=J3-1
              J1=J2-1
              J=LB(1,I)
              P(J1)=GASDEV(ISEED)*DESKET*SQRT(W(J))
              P(J2)=GASDEV(ISEED)*DESKET*SQRT(W(J))
              P(J3)=GASDEV(ISEED)*DESKET*SQRT(W(J))
           ENDDO
C
C    MAKE SURE THAT THE TOTAL LINEAR MOMENTUM EQUALS ZERO
C
           SUMX = 0.00
           SUMY = 0.00
           SUMZ = 0.00
           DO I=1,N
             J3 = 3 * LB(1,I)
             J2 = J3 - 1
             J1 = J2 - 1
             SUMX = SUMX + P( J1 )
             SUMY = SUMY + P( J2 )
             SUMZ = SUMZ + P( J3 )
           END DO
           SUMX = SUMX / REAL( N )
           SUMY = SUMY / REAL( N )
           SUMZ = SUMZ / REAL( N )
           DO I=1,N
             J3 = 3 * LB(1,I)
             J2 = J3 - 1
             J1 = J2 - 1
             P(J1) = P(J1) - SUMX
             P(J2) = P(J2) - SUMY
             P(J3) = P(J3) - SUMZ
           END DO
C
C   RESCALE THE VELOCITIES OF THE SYSTEM ACCORDING TO GIVEN TEMP.
C   REFERENCE:  "MOLECULAR DYNAMICS SIMULATION" BY JIM HAILE  P.458
C   K = 1.38066 * 10(-23) J/K = 0.00198624 KCAL/MOL K
C
           IF (NSCALE .GE. 0) THEN
              CALL THERMO(0)
           ENDIF

           NSE=NSCALE+NEQUAL

        ELSE

           DO I=1,N
             J3=3*LB(1,I)
             J2=J3-1
             J1=J2-1
             Q(J3)=QTEMP(J3)
             Q(J2)=QTEMP(J2)
             Q(J1)=QTEMP(J1)
             P(J3)=PTEMP(J3)
             P(J2)=PTEMP(J2)
             P(J1)=PTEMP(J1)
           END DO

           CALL DVDQ
           RAND=RAND0(ISEED)
           NSE=INT(RAND*20)+1

        ENDIF
          
        IF (NSE .GT. 0) THEN
           IF (INTEGRATOR.EQ.0) THEN
C
C         PERFORM THE FIRST SIX INTEGRATION CYCLES REQUIRED FOR ADAMSM
C
              CALL PARTI
              DO NC=1,6
                 CALL RUNGEK
                 DO IA=1,NATOMA(1)
                    J3=3*LA(1,IA)
                    J2=J3-1
                    J1=J2-1
                    Q(J3)=QTEMP(J3)
                    Q(J2)=QTEMP(J2)
                    Q(J1)=QTEMP(J1)
                    P(J3) = 0.0D0
                    P(J2) = 0.0D0
                    P(J1) = 0.0D0
                 END DO
                 IF (NSEL.EQ.1) THEN 
                    TB=0.D0
                    DO II=1,N
                      J=LB(1,II)
                      J3=3*J
                      J2=J3-1
                      J1=J2-1
                      TB=TB+(P(J1)**2+P(J2)**2+P(J3)**2)/W(J)
                    ENDDO
                    TEMPINIT=TB/(3.0D0*DBLE(N)*0.00198717D0*C1)
                    WRITE(30,*)'SYSTEM TEMPERATURE=',TEMPINIT
                 ENDIF
              ENDDO
              NSE=NSCALE+NEQUAL-6
           ELSEIF (INTEGRATOR.EQ.2) THEN
              CALL DVDQ
           ENDIF
C
C        RUN THE TRAJECTORY WITH NO A-B TRANSLATIONAL ENERGY TO
C        GET TO EQUILIBRIUM
C
           DO I=1,NSE
              NC=NC+1
              IF (INTEGRATOR.EQ.0) THEN
                 CALL ADAMSM 
              ELSEIF (INTEGRATOR.EQ.2) THEN
                 CALL SYMPLE(LLL)
              ELSEIF (INTEGRATOR.EQ.3) THEN
                 CALL VERLET(LLL,TELEC)
              ENDIF
              DO IA=1,NATOMA(1)
                 J3=3*LA(1,IA)
                 J2=J3-1
                 J1=J2-1
                 Q(J3)=QTEMP(J3)
                 Q(J2)=QTEMP(J2)
                 Q(J1)=QTEMP(J1)
                 P(J3) = 0.0D0
                 P(J2) = 0.0D0
                 P(J1) = 0.0D0
              END DO
              CALL ENERGY
              IF (NSEL.EQ.1) THEN 
                 TB=0.D0
                 DO II=1,N
                   J=LB(1,II)
                   J3=3*J
                   J2=J3-1
                   J1=J2-1
                   TB=TB+(P(J1)**2+P(J2)**2+P(J3)**2)/W(J)
                 ENDDO
                 TEMPINIT=TB/(3.0D0*DBLE(N)*0.00198717D0*C1)
                 WRITE(30,*)'SYSTEM TEMPERATURE=',TEMPINIT
              ENDIF
              IF (NC.LE.NSCALE.AND.MOD(NC,300).EQ.0) THEN
                 CALL THERMO(NC)
              ENDIF
           ENDDO
        ENDIF
C
C     SAVE THE COORDINATES AND MOMENTA FOR THE SURFACE, WHICH
C     WILL BE USED AS INITIAL VALUES FOR THE NEXT TRAJECTORY.
C
        DO I=1,N
           J3=3*LB(1,I)
           J2=J3-1
           J1=J2-1
           PTEMP(J3)=P(J3)
           PTEMP(J2)=P(J2)
           PTEMP(J1)=P(J1)
           QTEMP(J3)=Q(J3)
           QTEMP(J2)=Q(J2)
           QTEMP(J1)=Q(J1)
        END DO
C
C     SET THE COORDINATES AND MOMENTA FOR REACTANT A BACK TO
C     THE INITIAL CHOICE RECALL THAT THE MOMENTA WERE EQUATED
C     TO ZERO AND STORED IN PP.  THE COORDINATES HAVE THE Z
C     COMPONENT OF FRAGRMENT B AT 1000 A ABOVE THE SURFACE.
C     THIS WILL BE CORRECTED LATER IN THE CODE.
C
        N=NATOMA(1)
        DO I=1,N
           J3=3*LA(1,I)
           J2=J3-1
           J1=J2-1
           P(J3)=PP(J3)
           P(J2)=PP(J2)
           P(J1)=PP(J1)
           Q(J3)=QTEMP(J3)
           Q(J2)=QTEMP(J2)
           Q(J1)=QTEMP(J1)
        END DO

C-------end Bin, 2016/10/11

        WRITE(6,*)'EQUALIBRATION FOR NACTB EQ 7 IS NOW OVER'
C
C       RECOVER THE ORGINAL VALUES BEFORE SAMPLING
C
        TIME=TIMEORG
        ATIME=ATIMEORG
C
        NFINAL=0
        NCOOR=NCOORORG
        NC=0
159     CONTINUE
C
      ENDIF
C---------------------------------------------------------------
C
      N=NATOMA(1)
      DO I=1,N
         DO K=1,3
            J=3*LA(1,I)-3+K
            Q(J)=QQ(J)
            P(J)=PP(J)
         ENDDO 
      ENDDO
  160 CONTINUE
C
C         CHOOSE MOMENTA AND RELATIVE POSITIONS FOR REACTANTS.
C         USE THESE TO FIND REACTANTS' INITIAL P AND Q
C
C             GAS/SURFACE COLLISION
C
      IF (NSURF.NE.0) THEN
         CALL SURF(NSURF)
!    ADDED BY BIN, 2016/8/7
         IF (NSURF.EQ.3.AND.NGLO.NE.0) THEN
            CALL GLOSELECT(TVIBB)
c            CALL GLOEQU(TVIBB)  !ADDED BY BIN, 2017/2/5
         ENDIF
         CALL DVDQ
         CALL ENERGY
!    END BIN 
      ELSE
C
C             GAS PHASE COLLISION
C             SELECT IMPACT PARAMETER.  FIX POSITIONS OF A AND
C             B FOR CHOSEN IMPACT PARAMETER(B) AND SEPARATION(S).
C
         SB=BMAX
         IF (NOB.NE.1) THEN
            RAND=RAND0(ISEED)
            SB=BMAX*SQRT(RAND)
         ENDIF
         WRITE(6,162)SB
         DUM1=SQRT(S*S-SB*SB)
         N=NATOMB(1)
         DO 164 I=1,N
            J=3*LB(1,I)
            Q(J)=Q(J)+DUM1
            Q(J-1)=Q(J-1)+SB
  164    CONTINUE
C
C             IF NREL = 0 FOR GAS PHASE COLLISION
C             CHOOSE RELATIVE ENERGY FROM BOLTZMANN DISTRIBUTION
C
         IF (NREL.EQ.0) THEN
            DUM=GAMA(2,ISEED)
            SEREL=0.00198717D0*DUM*TRANS
C            WRITE(6,196)SEREL
            SEREL=SEREL*C1
         ENDIF
C
C             ADD RELATIVE TRANSLATIONAL ENERGY FOR GAS PHASE
C             COLLISION
C
         WT=WTA(1)+WTB(1)
         SDUM=WTA(1)*WTB(1)/WT
         DUM=SQRT(2.0D0*SEREL/SDUM)
         VELA=DUM*WTB(1)/WT
         VELB=VELA-DUM
         N=NATOMA(1)
         DO I=1,N
            J=3*LA(1,I)
            P(J)=P(J)+VELA*W(LA(1,I))
         ENDDO
         N=NATOMB(1)
         DO I=1,N
            J=3*LB(1,I)
            P(J)=P(J)+VELB*W(LB(1,I))
         ENDDO
         CALL DVDQ 
         CALL ENERGY
C
C             SAVE INITIAL RELATIVE VELOCITY AND ORBITAL ANGULAR
C             MOMENTUM FOR GAS PHASE COLLISION
C
         VI(1)=0.0D0
         VI(2)=0.0D0
         VI(3)=DUM
         VI(4)=DUM
         OAMI(1)=-SB*DUM*SDUM/C7
         OAMI(2)=0.0D0
         OAMI(3)=0.0D0
         OAMI(4)=ABS(OAMI(1))
         WRITE(6,198)(OAMI(I),I=1,3)
      ENDIF
C
  999 CONTINUE

!-->..MODIFIED BY BIN, 4/24/2014
C
C     NORMAL MODE ANALYSIS FOR EACH PRODUCT CHANNEL
C
      IF (NSURF.EQ.0) THEN
         IF (NSFLAG.NE.1) THEN
            IF(NTZ.EQ.1)THEN
C
C     SAVE INITIAL COORDINATES IN TEMPORARY STORAGE QTEMP.
C
                  N=NATOMA(1)+NATOMB(1)
                  DO J=1,3*N
                     QTEMP(J)=Q(J)
                  ENDDO
C
C     READ PRODUCT EQUILIBRIUM COORDINATES AND SEPARATE PRODUCTS
C
                  DO I=2,NPATHS+1
                     N=NATOMA(I)
                     DO K=1,3*N
                        Q(K)=QZA(I,K)
                     ENDDO
                     M=NATOMB(I)
                     DO K=1,3*M
                        Q(3*N+K)=QZB(I,K)
                     ENDDO
                     DO K=1,N
                        K3=3*K
                        Q(K3)=QZA(I,K3)+ZASYM
                     ENDDO
C
C     NORMAL MODE ANALYSIS FOR BOTH PRODUCTS
C
                     IF (N.GE.2) THEN
                     WRITE(26,*)'NORMAL MODES FOR FRAGMENT A IN PATH ',I
                     WRITE(26,*)
                     CALL NMODE(N,0)
                     ENDIF
                     IF (M.GE.2) THEN
                     WRITE(26,*)'NORMAL MODES FOR FRAGMENT B IN PATH ',I
                     WRITE(26,*)
                     CALL NMODE(M,N)
                     ENDIF
                  ENDDO
C
C     RESTORE INITIAL COORDINATES 
C
                  N=NATOMA(1)+NATOMB(1)
                  DO J=1,3*N
                     Q(J)=QTEMP(J)
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
!-->..END
C
C             SAVE CHOSEN INITIAL ROTATIONAL ANGULAR MOMENTUM
C
      IF (NATOMA(1).GT.1) THEN 
         AMAI(1)=AMA(1)/C7
         AMAI(2)=AMA(2)/C7
         AMAI(3)=AMA(3)/C7
         AMAI(4)=SQRT(AMAI(1)**2+AMAI(2)**2+AMAI(3)**2)
!-->..MODIFIED BY BIN, 2/4/2017
         IF (NATOMA(1).EQ.2) THEN
            WRITE(6,206)EROTA,AMAI(1),AMAI(2),AMAI(3)
            WRITE(6,208)ENJA/C1
         ELSE
            WRITE(6,206)ERAI,AMAI(1),AMAI(2),AMAI(3)
            WRITE(6,208)ETAI-ERAI
         ENDIF
!.....END
      ENDIF

      IF (NATOMB(1).GT.1) THEN
         AMBI(1)=AMB(1)/C7
         AMBI(2)=AMB(2)/C7
         AMBI(3)=AMB(3)/C7
         AMBI(4)=SQRT(AMBI(1)**2+AMBI(2)**2+AMBI(3)**2)
!-->..MODIFIED BY BIN, 2/4/2017
         IF (NATOMA(1).EQ.2) THEN
            WRITE(6,207)EROTB,AMBI(1),AMBI(2),AMBI(3)
            WRITE(6,209)ENJB/C1
         ELSE
           WRITE(6,207)ERBI,AMBI(1),AMBI(2),AMBI(3)
           WRITE(6,209)ETBI-ERBI
         ENDIF
!.....END
      ENDIF
C
!-->..MODIFIED BY BIN, 9/11/2014
      WRITE(6,196)SEREL/C1
!.....END
      NSFLAG=1
      RETURN
      END
