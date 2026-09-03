program VENUS
  use input_parser, only: get_real
  use venus_params
  use venus_data
  use venus_input, only: read_venus_input, NMA, ENU, EDELTU, ENL, EDELTL, &
                           ISEED, NIP, NCROT, MPLOT, ENSAV, ENSQ
  implicit double precision (A-H,O-Z)
!***********************************************************************
!                                                                      *
!                              VENUS10                                 *
!                                                                      *
!              A GENERAL CHEMICAL DYNAMICS COMPUTER PROGRAM            *
!                                                                      *
!                                 BY                                   *
!                                                                      *
!             WILLIAM L. HASE, KIM BOLTON, NAVDEEP CHAWLA,             *
!        PASCAL DE SAINTE CLAIRE, RONALD J. DUCHOVIC, XICHE HU,        *
!       ANDREW KOMORNICKI,GUOSHENG LI, KIERAN F. LIM, DA-HONG LU,      *
!     SAMY MEROUEH, KYOYEON PARK, GILLES H. PESLHERBE, KIHYUNG SONG,   *
!          LIPENG SUN, KANDADAI N. SWAMY, SCOTT R. VANDE LINDE,        *
!       ANTONIO VARANDAS, HAOBIN WANG, RALPH J. WOLF, MINGYING XUE,    *
!                            AND TIANYING YAN                          *
!                                                                      *
!                               JUNE, 2010                             *
!                                                                      *
!***********************************************************************
! NOTE: THIS IS A MODIFIED VERSION OF VENUS, AIMING TO REMOVE THOSE
!      OLD SUBROUTINES THAT ARE NOT USEFULL FOR OUR PURPOSE AND ADD SOME
!      SOME NEW FEATURES THAT ARE NOT AVAILABLE IN ORIGINAL VENUS.
!
!      MORE READABLE, WHICH IS USFUL FOR FUTURE EXPANSION.
! BY BIN JIANG, 7/30/2016
!
! --- REFACTORED: All GOTO statements replaced with structured loops  ---
!
! --- All COMMON blocks moved to venus_data module ---
!
  dimension ERAVA(500),ERAVB(500)
  dimension ETIM(2000),ESAV(2000),ESQ(2000),NEVIB(2000), &
     NEVIBU(2000),NEVIBL(2000),NEVIBN(2000)
  dimension EBM(NDA),ENM(NDA),EBSAV(2,2000),EBSQ(2,2000), &
     MNLM(50)
  dimension QDUM(NDA3),TABDUM(42*NDA),GDUM(NDP),NFDUM(NDP)
  dimension QSAVE(NDA3),PSAVE(NDA3)

  character*30 ADATE
  save TUSER1,TUSER2,TTOTAL1,TTOTAL2
  save XL

  logical :: trajectory_done

  integer :: I,J,K,NCZ,NXPLOT,ISAV,N,M,IPATH
  integer :: Nbin, idplt
  integer :: NDUM, NEFF, NSEARCH, NSELECT
  integer :: IBAR_LOCAL, KRE
  double precision :: START_TIME, FINISH_TIME, TUSER1, TUSER2, TTOTAL1, TTOTAL2
  double precision :: FACTOR1, FACTOR2, BOLTZ, RBOHR, FACTOR, FACTOR3
  double precision :: TI, ENN, EDELTN, UTMP
  double precision :: Ebin, Qr1, Qr2
  double precision :: P_z0

  891 FORMAT(/'   THE OPTION NCHKP=-1 IS NOT AVAILABLE FOR ', &
     'TRAJECTORY CALCULATIONS'/)
  515 FORMAT(//' USER AND REAL TIME FOR VENUS:', &
     2F8.1,' SECONDS'//)
  520 FORMAT(//' USER AND REAL TIME FOR SELECT:', &
     2F8.1,' SECONDS'//)
!
!        INITIALIZE ARRAYS AND PARAMETERS.
!
  CALL CPU_TIME(START_TIME)

  DO I=1,2000
     NEVIB(I)=0
     NEVIBU(I)=0
     NEVIBL(I)=0
     NEVIBN(I)=0
     ESAV(I)=0.0D0
     ESQ(I)=0.0D0
  ENDDO

!
!        CONSTANTS — now in venus_params (C1..C8, PI, HALFPI, TWOPI)
!        RKUTTA coefficients computed below
!
  RAA1=1.0D0/2.0D0
  RA1=1.0D0-SQRT(2.0D0)/2.0D0
  RA2=2.0D0*RA1
  RA3=2.0D0-3.0D0*SQRT(2.0D0)/2.0D0
  RB1=2.0D0-RA1
  RB2=2.0D0*RB1
  RB3=4.0D0-RA3
  RC1=1.0D0/6.0D0
  RC2=1.0D0/3.0D0
!
  GAO=0.0D0
  NSAD=0
  NCBA=0
  NCAB=0
  IBAR=0
  NPTS=2
  HINC=0.001D0
  TRVA = -1.0
  TRVB = -1.0

  NONREACT=0
  NREACT=0
  NLONG=0
  ZASYM=10.0D0
  DELM=1.0D4
  NCHKP=0
!
!        READ INPUT AND SETUP INITIAL CONFIGURATION
!
      call read_venus_input()

! Output for product normal mode analysis
  WRITE(26,*)' XXXXXXXXXXXXXXXXXXXXXXXX', &
     ' PRODUCT NORMAL MODE ANALYSIS', &
     ' XXXXXXXXXXXXXXXXXXXXXXXX'
  WRITE(26,*)
  WRITE(26,*)'NSURF=',NSURF,'NPATHS=',MPATH
  WRITE(26,*)

  IF (NSURF.GT.0) THEN
     WRITE(26,*)'FOR GAS-SURFACE REACTIONS, ONLY THESE PATHS WITH AN', &
        ' GAS PHASE SPECIES ARE ANALYZED'
  ELSE
     WRITE(26,*)'FOR GASEOUS REACTIONS, ALL PATHS WILL BE ANALYZED'
  ENDIF

  WRITE(26,*)
  WRITE(26,'("   MASSES OF ATOMS:",I5," ATOMS"/)' ) NATOMS
  WRITE(26,'(6F11.5)' ) (W(I),I=1,NATOMS)
  WRITE(26,*)

  DO I=1, NPATHS+1
     IF (NDELH(I).EQ.0) THEN
        WRITE(26,*)'EQUILIBRIUM COORDINATES FOR PATH ',I
        WRITE(26,*)
        WRITE(26,*)'NATOMA(I)=',NATOMA(I),'NATOMB(I)=',NATOMB(I)
        WRITE(26,*)
        WRITE(26,'(5X,"EQUILIBRIUM COORDINATES FOR A:")')
        K=3*NATOMA(I)
        WRITE(26,'(7X,F10.6,4X,F10.6,4X,F10.6)' )(QZA(I,J),J=1,K)
        IF (NSURF.LE.1) THEN
           WRITE(26,'(5X,"EQUILIBRIUM COORDINATES FOR B:")')
           K=3*NATOMB(I)
           WRITE(26,'(7X,F10.6,4X,F10.6,4X,F10.6)' )(QZB(I,J),J=1,K)
        ENDIF
        WRITE(26,*)
     ENDIF
  ENDDO
!
! --- Integration setup ---
!
  IF(NSELT.GE.0)THEN
     ATIME=TIME
     IF (INTEGRATOR.EQ.1) THEN
        XL=ATIME
        IF(LLL.LT.0) THEN
           TF=XL
        ELSE
           TF=ATIME*DBLE(NS)
        ENDIF
     ENDIF
  ENDIF
  ! Original 432 CONTINUE
!
!         SET FLAGS AND PARAMETERS
!
  NAM=0
  NI=I3N-3*NRGD
  NID=2*NI
  NTZ=0

  CALL CPUSEC(TUSER1,TTOTAL1)
  SECADD=0.0D0
!
!         CHECKPOINT HANDLING
!
  IF (NCHKP == 0) THEN
     ! Original 223: initialize random numbers
     IF (NSELT.EQ.2.OR.NSELT.EQ.3) CALL RANDST(ISEED)
  ELSE
     ! NCHKP != 0: checkpoint restart logic
     IF (NCHKP.EQ.-2) THEN
        ! Read random number array (currently commented out)
        ! Original GOTO 451: skip to trajectory increment
        ! (handled by the trajectory loop below)
     ELSE
        KRE=1
        CALL DVDQ_1
        CALL ENERGY_1

        IF (NSELT.EQ.0) THEN
           NDUM=MAX0(1,NTZ)
           IF (NC.EQ.NS) NDUM=NDUM-1
           DO J=1,NDUM
              READ(5,*)(QDUM(I),I=1,I3N)
              READ(5,*)(QDUM(I),I=1,I3N)
           ENDDO
           IF (NCHKP.EQ.-1) THEN
              NC=0
              NTZ=1
              NX=NIP
           ENDIF
        ENDIF
        CALL FLUSH(6)

        IF (NC.EQ.NS.AND.NCHKP.EQ.1) THEN
           NTZ=NTZ-1
           ! Original GOTO 451 (handled by loop below)
        ELSE IF (NSELT.EQ.0.AND.NCHKP.EQ.-1) THEN
           ! Original GOTO 425: skip to integration
        ELSE IF (NSELT.EQ.0) THEN
           ! Original GOTO 400: start integration
        ELSE
           IF (NCHKP.EQ.-1) THEN
              WRITE(6,891)
              STOP
           ENDIF
           ! Original GOTO 400
        ENDIF
        ! After checkpoint handling, fall through to trajectory loop
     ENDIF
  ENDIF

  ! Initialize random sequence if not done
  IF (NCHKP == 0) THEN
     IF (NSELT.EQ.2.OR.NSELT.EQ.3) CALL RANDST(ISEED)
  ENDIF

!========================================================================
!  MAIN TRAJECTORY LOOP  (replaces: 451 NTZ=NTZ+1 ... GOTO 451)
!========================================================================


  TRAJECTORY_LOOP: DO
     NTZ = NTZ + 1

     ! Original 451: check for completion
     IF (NTZ > NT) THEN
        CALL CPUSEC(TUSER2,TTOTAL2)
        WRITE(6,515) (SECADD+TUSER2-TUSER1),(SECADD+TTOTAL2-TTOTAL1)
        TUSER1=TUSER2
        TTOTAL1=TTOTAL2

        CALL CPU_TIME(FINISH_TIME)
        WRITE(6,*)'+++++++++++++++++++++++++++++++++++++++++'
        WRITE(6,*)'FINAL! ','TOTAL TRAJECTORIES:',NT
        WRITE(6,*)'NON-REACTIVE TRAJECTORY NUMBERS FOR PATH 1 ',NONREACT
        NEFF=NONREACT
        DO IPATH=1,NPATHS
           WRITE(6,'(A,I6,I10)')' REACTIVE TRAJECTORY NUMBERS FOR PATH ', &
              IPATH+1,NREACT(IPATH)
           NEFF=NEFF+NREACT(IPATH)
        ENDDO
        WRITE(6,*)'TRAJECTORIES WITH TOO LONG TIME ',NLONG
        WRITE(6,*)'TOTAL EFFECTIVE TRAJECTORIES ',NEFF
        WRITE(6,*)
        WRITE(6,*)'+++++++++++++++++++++++++++++++++++++++++'
        WRITE(6,*)'TOTAL TIME COST ',FINISH_TIME-START_TIME,' SECONDS'
        WRITE(6,*)'TOTAL TRAJECTORY NUMBER ', NT
        WRITE(6,*)'MEAN TIME COST PER ONE TRAJECTORY:', &
           (FINISH_TIME-START_TIME)/NT, 'SECONDS'

        WRITE(999,*)'+++++++++++++++++++++++++++++++++++++++++'
        WRITE(999,*)'FINAL! ','TOTAL TRAJECTORIES:',NT
        WRITE(999,*)'NON-REACTIVE TRAJECTORY NUMBERS FOR PATH 1 ',NONREACT
        NEFF=NONREACT
        DO IPATH=1,NPATHS
           WRITE(999,'(A,I6,I10)')' REACTIVE TRAJECTORY NUMBERS FOR PATH ', &
              IPATH+1,NREACT(IPATH)
           NEFF=NEFF+NREACT(IPATH)
        ENDDO
        WRITE(999,*)'TRAJECTORIES WITH TOO LONG TIME ',NLONG
        WRITE(999,*)'TOTAL EFFECTIVE TRAJECTORIES ',NEFF
        WRITE(999,*)
        WRITE(999,*)'+++++++++++++++++++++++++++++++++++++++++'
        WRITE(999,*)'TOTAL TIME COST ',FINISH_TIME-START_TIME,' SECONDS'
        WRITE(999,*)'TOTAL TRAJECTORY NUMBER ', NT
        WRITE(999,*)'MEAN TIME COST PER ONE TRAJECTORY:', &
           (FINISH_TIME-START_TIME)/NT, 'SECONDS'

        ! Histogram plots
        Ebin = 40d0
        Nbin = 40
        Qr1 = 2.94998
        Qr2 = 2.55476
        idplt = 0
        if(idplt == 1) then
           call SYSTEM('mkdir HISTO')
           call plot_Erel_hist(0d0, Ebin, Nbin)
           call plot_theta_hist(0d0, 90d0, 30)
           call plot_Erel_impact_pos()
           call plot_site_histograms(0d0, Ebin, Nbin,0d0,90d0, 30,Qr1,Qr2)
           call plot_Erel_by_bounce(0d0, Ebin, Nbin,5)
        end if
        STOP
     ENDIF

     ! Save random seeds
     DO I=1,8
        ISEED0(I)=ISEED3(I)
     ENDDO
     VRELO=0.0D0
     INTST=0
     NAST=2
     NFINAL=0
     if (TRVA >= 0.0d0) NSFLAG = 0
     if (TRVB >= 0.0d0) NSFLAG = 0
     KRE=1
     VMAX=-1.0D20

     ! Initialize Q and P arrays
     J=3*NATOMA(1)
     DO I=1,J
        Q(I)=QZA(1,I)
        P(I)=0.0D0
     ENDDO
     IF (NSURF.EQ.2.AND.NGLO.NE.0) THEN
        J=3*(NATOMA(1)+2)
        K=3*NATOMB(1)
        DO I=1,K
           Q(J+I)=QZB(1,I)
           P(J+I)=0.0D0
        ENDDO
     ELSE
        K=3*NATOMB(1)
        DO I=1,K
           Q(J+I)=QZB(1,I)
           P(J+I)=0.0D0
        ENDDO
     ENDIF

     ! Enforce surface alignment before SELECT modifies Q.
     ! Only meaningful for a polyatomic fragment A (aligns its surface bond);
     ! a monatomic A is spherically symmetric - skip the alignment stub.
     if (NZDOWN == 1 .and. NATOMA(1) > 1) call FIXROTDATM(1)

     ! If user specified a fixed incident angle (NTHTA >= 0), ensure
     ! THTA is correctly set in radians before SELECT overwrites it.
     if (NSURF >= 1 .and. NTHTA >= 0) then
        THTA = get_real('THTA', 0.0d0) * C4
     end if

     ! Select initial conditions
     CALL SELECT
     NSELECT = 1
     NSEARCH=0
     MINENE=1D5

     ! Compute time for SELECT
     CALL CPUSEC(TUSER2,TTOTAL2)
     WRITE(6,520) (SECADD+TUSER2-TUSER1),(SECADD+TTOTAL2-TTOTAL1)
     TUSER1=TUSER2
     TTOTAL1=TTOTAL2
     CALL FLUSH(6)

     ! Mode populations initialization
     TI=0.0D0
     IF (MPLOT == 1) THEN
        ISAV=1
        ETIM(ISAV)=TI
        IF (NMA.NE.0) THEN
           IF (MNTR.NE.0) THEN
              NEVIBN(ISAV)=NEVIBN(ISAV)+1
              I=MNTR
              WWA(I)=WWA(I)/C6*CM2CAL
              ENN=(ANQA(I)+0.5D0)*WWA(I)
              EDELTN=0.5D0*WWA(I)
              WWA(I)=WWA(I)*C6/CM2CAL
           ENDIF
           CALL ENMODE(ENM,NMA)
           DO I=1,NMA
              ENSAV(I,ISAV)=ENSAV(I,ISAV)+ENM(I)
              ENSQ(I,ISAV)=ENSQ(I,ISAV)+ENM(I)*ENM(I)
           ENDDO
        ENDIF
        IF (NACTA.EQ.4) THEN
           CALL EBOND(EBCH,EKCH,RCH,NEXM)
           NEVIB(ISAV)=NEVIB(ISAV)+1
           ESAV(ISAV)=ESAV(ISAV)+EBCH
           ESQ(ISAV)=ESQ(ISAV)+EBCH**2
           IF (NLM.NE.0) THEN
              DO I=1,NLM
                 J=MNLM(I)
                 CALL EBOND(EBM(I),EK,RCH,J)
                 EBSAV(I,ISAV)=EBSAV(I,ISAV)+EBM(I)
                 EBSQ(I,ISAV)=EBSQ(I,ISAV)+EBM(I)*EBM(I)
              ENDDO
           ENDIF
        ENDIF
     ENDIF

     ! Original 425 CONTINUE
     ! Electronic structure initialization
     z_pos_min = S
     BOUNCE = 0

     CALL ENERGY_1
     CALL DVDQ_1

     CALL GWRITE

     NC=0
     NX=NIP
     NXPLOT=NPLOT

! Save initial coordinates and momenta for restoration upon termination
     QSAVE = Q
     PSAVE = P

!=====================================================================
!  INTEGRATION LOOP  (replaces: 402/400 ... GOTO 400)
!=====================================================================
     trajectory_done = .false.
     INTEGRATION_LOOP: DO WHILE (.NOT. trajectory_done)
        ! Original 402: integration step
        NCZ=NC

        IF (INTEGRATOR.EQ.2) THEN
           CALL DVDQ_1
        ENDIF

        ! Original 400: NC=NC+1
        NC=NC+1

        IF (INTEGRATOR.EQ.1) THEN
           IF(LLL.GE.0) THEN
              WRITE(6,*)'DOING TIME VARIABLE INTEGRATION'
              CALL RADAU(TF,XL,LLL,NIP,NITER)
              IF (NTZ.GT.NT) THEN
                 CALL CPUSEC(TUSER2,TTOTAL2)
                 WRITE(6,515) (SECADD+TUSER2-TUSER1), &
                    (SECADD+TTOTAL2-TTOTAL1)
                 TUSER1=TUSER2
                 TTOTAL1=TTOTAL2
                 STOP
              ENDIF
              ! Original GOTO 451: next trajectory
              trajectory_done = .true.
              CYCLE INTEGRATION_LOOP
           ELSE
              CALL RADAU(TF,XL,LLL,NIP,NITER)
           ENDIF
        ELSEIF (INTEGRATOR.EQ.2) THEN
           CALL SYMPLE(LLL,NSELECT)
           IF (NCOOR.EQ.1.AND.MOD(NC,NIP).EQ.0) CALL GWRITE
           if(H.lt.HTMIN) HTMIN=H
           if(H.gt.HTMAX) HTMAX=H
        ELSEIF (INTEGRATOR.EQ.3) THEN
           P_z0 = P(3)
           CALL VERLET(LLL,0.0d0,NSELECT)

           if(p(3)*P_z0 < 0 .AND. P(3)<0) BOUNCE = BOUNCE+1
           if(Q(3)<z_pos_min) z_pos_min = Q(3)
           IF (NCOOR.EQ.1.AND.MOD(NC,NIP).EQ.0) CALL GWRITE
           if(H.lt.HTMIN) HTMIN=H
           if(H.gt.HTMAX) HTMAX=H
           if(HTMAX-HTMIN>8d0) then
              ! Energy not conserved, but proceed
           end if
        ENDIF

        ! Mode population tracking (original <= 426)
        IF (MPLOT == 1 .AND. NC >= NXPLOT) THEN
           NXPLOT=NXPLOT+NPLOT
           TI=DBLE(NC)*TIME
           IF (ISAV < 2000) THEN
              ISAV=ISAV+1
              ETIM(ISAV)=TI
              IF (NMA.NE.0) THEN
                 CALL ENMODE(ENM,NMA)
                 IF (MNTR.NE.0) THEN
                    I=MNTR
                    IF (ABS(ENN-ENM(I)).LT.EDELTN) NEVIBN(ISAV)=NEVIBN(ISAV)+1
                 ENDIF
                 DO I=1,NMA
                    ENSAV(I,ISAV)=ENSAV(I,ISAV)+ENM(I)
                    ENSQ(I,ISAV)=ENSQ(I,ISAV)+ENM(I)*ENM(I)
                 ENDDO
              ENDIF
              IF (NACTA.EQ.4) THEN
                 CALL EBOND(EBCH,EKCH,RCH,NEXM)
                 IF (ABS(ENON-EBCH).LT.EDELTA) NEVIB(ISAV)=NEVIB(ISAV)+1
                 IF (ABS(ENU-EBCH).LT.EDELTU) NEVIBU(ISAV)=NEVIBU(ISAV)+1
                 IF (ABS(ENL-EBCH).LT.EDELTL) NEVIBL(ISAV)=NEVIBL(ISAV)+1
                 ESAV(ISAV)=ESAV(ISAV)+EBCH
                 ESQ(ISAV)=ESQ(ISAV)+EBCH*EBCH
                 IF (NLM.NE.0) THEN
                    DO I=1,NLM
                       J=MNLM(I)
                       CALL EBOND(EBM(I),EK,RCH,J)
                       EBSAV(I,ISAV)=EBSAV(I,ISAV)+EBM(I)
                       EBSQ(I,ISAV)=EBSQ(I,ISAV)+EBM(I)*EBM(I)
                    ENDDO
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
        ! Original 426 CONTINUE

        ! Check for NFINAL == 1 (final analysis triggered)
        IF (NFINAL == 1) THEN
           ! Original 414: final analysis
           CALL FINAL

           ! Rotational energy averaging
           ERAVA(KRE)=EROTA
           ERAVB(KRE)=EROTB
           KRE=KRE+1
           IF (NSURF.GT.0.AND.NPATH.NE.1) THEN
              SDA=EROTA
              SDB=EROTB
           ELSEIF (KRE <= NCROT) THEN
              ! Need more rotational samples — return to integration
              CYCLE INTEGRATION_LOOP
           ELSE
              EROTA=0.0D0
              EROTB=0.0D0
              DO I=1,NCROT
                 EROTA=EROTA+ERAVA(I)
                 EROTB=EROTB+ERAVB(I)
              ENDDO
              EROTA=EROTA/DBLE(NCROT)
              EROTB=EROTB/DBLE(NCROT)
              SDA=0.0D0
              SDB=0.0D0
              DO I=1,NCROT
                 SDA=SDA+(EROTA-ERAVA(I))**2
                 SDB=SDB+(EROTB-ERAVB(I))**2
              ENDDO
              IF (NCROT.GT.1) THEN
                 SDA=SQRT(SDA/(NCROT-1))
                 SDB=SQRT(SDB/(NCROT-1))
              ELSE
                 SDA=0.0D0
                 SDB=0.0D0
              ENDIF
           ENDIF
           CALL GFINAL
           AVfinS(NTZ)= V
           NBncS(NTZ) = BOUNCE
           AZminS(NTZ) = z_pos_min
           AInitX(NTZ) = RX0
           AInitY(NTZ) = RY0
           call POT0(1,E0)
           call record_all

           ! Reaction counting
           WRITE(6,*)
           WRITE(6,'(I7,E12.5,A,I7)')NTZ,HTMAX-HTMIN,' KCAL/MOL ',NPATH
           WRITE(6,*)
           CALL FLUSH(6)
           WRITE(999,'(I7,E12.5,A,I7)')NTZ,HTMAX-HTMIN,' KCAL/MOL ',NPATH
           CALL FLUSH(999)
           IF(NFC.EQ.0.AND.NGLO.EQ.0) THEN
              IF((HTMAX-HTMIN).LE.DELM) THEN
                 IF(NPATH.EQ.1) NONREACT=NONREACT+1
                 DO IPATH=1,NPATHS
                    IF(NPATH.EQ.IPATH+1) NREACT(IPATH)=NREACT(IPATH)+1
                 ENDDO
              ENDIF
           ELSE
              IF(NPATH.EQ.1) NONREACT=NONREACT+1
              DO IPATH=1,NPATHS
                 IF(NPATH.EQ.IPATH+1) NREACT(IPATH)=NREACT(IPATH)+1
              ENDDO
           ENDIF
           IF (NCOOR.EQ.1) CLOSE(1000+NTZ)

           trajectory_done = .true.
           CYCLE INTEGRATION_LOOP
        ENDIF

        ! Check for reaction events via TEST
        CALL TEST

        IF (NTEST == 2) THEN
           ! Original 410: reaction event detected
           IF (NPATH == 1 .AND. NTEST == 2) THEN
              ! Stable point search (disabled in original via goto 414)
              ! Original: goto 414 bypasses the search below
              ! (stable point search code omitted — dead code in original)
           END IF
           if(NPATH==1) then
              !Q=QSAVE
              !P=PSAVE
           end if

           CALL GWRITE

           ! Original 414: Final analysis
           CALL FINAL

           ERAVA(KRE)=EROTA
           ERAVB(KRE)=EROTB
           KRE=KRE+1
           IF (NSURF.GT.0.AND.NPATH.NE.1) THEN
              SDA=EROTA
              SDB=EROTB
           ELSEIF (KRE <= NCROT) THEN
              ! Need more rotational samples — return to integration
       !       CYCLE INTEGRATION_LOOP  2026.5.12
           ELSE
              EROTA=0.0D0
              EROTB=0.0D0
              DO I=1,NCROT
                 EROTA=EROTA+ERAVA(I)
                 EROTB=EROTB+ERAVB(I)
              ENDDO
              EROTA=EROTA/DBLE(NCROT)
              EROTB=EROTB/DBLE(NCROT)
              SDA=0.0D0
              SDB=0.0D0
              DO I=1,NCROT
                 SDA=SDA+(EROTA-ERAVA(I))**2
                 SDB=SDB+(EROTB-ERAVB(I))**2
              ENDDO
              IF (NCROT.GT.1) THEN
                 SDA=SQRT(SDA/(NCROT-1))
                 SDB=SQRT(SDB/(NCROT-1))
              ELSE
                 SDA=0.0D0
                 SDB=0.0D0
              ENDIF
           ENDIF
           CALL GFINAL
           AVfinS(NTZ)= V
           NBncS(NTZ) = BOUNCE
           AZminS(NTZ) = z_pos_min
           AInitX(NTZ) = RX0
           AInitY(NTZ) = RY0
           call POT0(1,E0)
           call record_all

           WRITE(6,*)
           WRITE(6,'(I7,E12.5,A,I7)')NTZ,HTMAX-HTMIN,' KCAL/MOL ',NPATH
           WRITE(6,*)
           CALL FLUSH(6)
           WRITE(999,'(I7,E12.5,A,I7)')NTZ,HTMAX-HTMIN,' KCAL/MOL ',NPATH
           CALL FLUSH(999)
           IF(NFC.EQ.0.AND.NGLO.EQ.0) THEN
              IF((HTMAX-HTMIN).LE.DELM) THEN
                 IF(NPATH.EQ.1) NONREACT=NONREACT+1
                 DO IPATH=1,NPATHS
                    IF(NPATH.EQ.IPATH+1) NREACT(IPATH)=NREACT(IPATH)+1
                 ENDDO
              ENDIF
           ELSE
              IF(NPATH.EQ.1) NONREACT=NONREACT+1
              DO IPATH=1,NPATHS
                 IF(NPATH.EQ.IPATH+1) NREACT(IPATH)=NREACT(IPATH)+1
              ENDDO
           ENDIF
           IF (NCOOR.EQ.1) CLOSE(1000+NTZ)

           trajectory_done = .true.
           CYCLE INTEGRATION_LOOP
        ENDIF

        ! Update state flags
        IF (NTEST == 0 .AND. NAST /= 0) NAST=0
        IF (NTEST == 1 .AND. NAST == 2) NAST=1

        ! Check for max integration cycles (original 450)
        IF (NC >= NS) THEN
           ! Long trajectory
           AERELS(NTZ)= -1d0
           ATHTAS(NTZ)= -1d0
           AVfinS(NTZ)= -1d0
           NBncS(NTZ) = BOUNCE
           AZminS(NTZ) = z_pos_min
           AInitX(NTZ) = RX0
           AInitY(NTZ) = RY0
           call POT0(1,E0)
           call record_all

           IF(NFC.EQ.0.AND.NGLO.EQ.0) THEN
              IF((HTMAX-HTMIN).LE.DELM) NLONG=NLONG+1
           ELSE
              NLONG=NLONG+1
           ENDIF

           ! Original 450: write long trajectory info
           CALL ENERGY_1
           WRITE(999,'(I7,E12.5,A,I7)')NTZ,HTMAX-HTMIN,' -KCAL/MOL- ',NPATH
           CALL FLUSH(999)
           WRITE(6,*)
           WRITE(6,*)'TRAJECTORIES WITH LONG LIFE TIME'
           WRITE(6,'(I7,E12.5,A,I7)')NTZ,HTMAX-HTMIN,' KCAL/MOL ',NPATH
           CALL FLUSH(6)
           IF (NCOOR.EQ.1) CLOSE(1000+NTZ)

           trajectory_done = .true.
           CYCLE INTEGRATION_LOOP
        ENDIF

        ! Continue integration (original GOTO 400)
     END DO INTEGRATION_LOOP

     ! End of one trajectory; loop back for next (original GOTO 451)
  END DO TRAJECTORY_LOOP

end program VENUS
