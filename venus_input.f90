module venus_input
  use input_parser
  use venus_params
  use venus_data
  use param_mapping
  use presets
  implicit none
  private

  public :: read_venus_input
  public :: NMA, ENU, EDELTU, ENL, EDELTL
  public :: ISEED, NIP, NCROT, MPLOT, ENSAV, ENSQ

  ! Persistent variables needed by main trajectory loop
  integer :: NMA = 0
  real(8) :: ENU   = 0.0d0
  real(8) :: EDELTU = 0.0d0
  real(8) :: ENL   = 0.0d0
  real(8) :: EDELTL = 0.0d0
  integer :: ISEED = 0, NIP = 10, NCROT = 1, MPLOT = 0
  real(8) :: ENSAV(30,2000), ENSQ(30,2000)

  ! Private scratch variables
  integer :: I, J, K, II, JJ, JL, IL, ICB, M, N
  real(8) :: DUM, DUM1, QZDUM, WT, CBMAX, DETCB, DETCB1, DETCB2, DETCB3
  real(8) :: RRA2, RRB2
  integer :: J1, J2, J3, J4, J5, J6
  real(8) :: QCM(3), VCM(3), AM(4)


contains

  subroutine read_venus_input()
    ! Local variables (only needed within this subroutine)
    character(len=64) :: kwbuf
    character(len=512) :: errmsg
    integer :: ios
    real(8) :: RSEED
      ! Input-only variables (not shared with main program)
      integer :: NCVIB, NERTEST
      real(8) :: EROT, RHO
      ! F9 fix (ported from isolated copy 7eb55ae): local 3x3 staging for
      ! EIGN — module CA/CB are (NDA3yf,NDA3yf) work arrays (also used by
      ! SELECT/ENMODE at full width), but EIGN declares its dummies as
      ! A(NN,NN); passing the big arrays with NN=3 mismatches column-major
      ! layouts — EIGN writes eigenvectors at linear offsets 1-9 (its 3x3
      ! view) while this routine reads CB(1:3,2)/CB(1:3,3) at offsets
      ! 91-93/181-183 (the 90-wide view), never-written static storage = 0.
      ! Downstream: QZA recomputed as CB·QQ collapses to zero -> main
      ! program starts every trajectory at the origin -> NaN Hessian.
      real(8) :: CA3(3,3), CB3(3,3)
      ! Intrinsics that need explicit declaration under implicit none
      intrinsic :: dble, mod
  ! FORMAT statements
  802 FORMAT(/1X,A30/)
  809 FORMAT(3X,8F9.2)
  810 FORMAT(6F11.5)
  811 FORMAT('   NT=',I5,'  NS=',I10,'  NIP=',I10,'  NCROT=',I6)
  812 FORMAT('   ISEED=',I12,'  TIME=',F9.5)
  813 FORMAT('   ACTIVATE WITH ORTHANT SAMPLING'/)
  814 FORMAT(' NUMBER OF ATOMS=',I5/ &
     ' NUMBER OF HARMONIC STRETCHES=',I5/ &
     ' NUMBER OF MORSE STRETCHES=',I5/ &
     ' NUMBER OF HARMONIC BENDS=',I5/ &
     ' NUMBER OF ALPHA BENDS=',I5/ &
     ' NUMBER OF LENNARD-JONES INTERACTIONS=',I5/ &
     ' NUMBER OF TORSIONS=',I5/ &
     ' NUMBER OF REPULSIONS=',I5/ &
     ' NUMBER OF GHOST PAIRS=',I5/ &
     ' NUMBER OF TETRAHEDRAL CENTERS=',I5/ &
     ' NUMBER OF R-R COUPLINGS=',I5/ &
     ' NUMBER OF R-THETA COUPLINGS=',I5/ &
     ' NUMBER OF THETA-THETA COUPLINGS=',I5/ &
     ' NUMBER OF DIHEDRAL ANGLES=',I5/ &
     ' NUMBER OF AXILROD-TELLER THREE-BODY INTERACTIONS=',I5/ &
     ' CHOICE OF SN2 POTENTIALS=',I5/ &
     ' NUMBER OF RYDBERG POTENTIAL FUNCTIONS=',I5/ &
     ' NUMBER OF HARTREE-FOCK DISPERSION FUNCTIONS=',I5/ &
     ' NUMBER OF GENERALIZED LEPS(A) FUNCTIONS=',I5/ &
     ' NUMBER OF GENERALIZED LEPS(B) FUNCTIONS=',I5/ &
     ' CHOICE OF DMBE POTENTIAL ENERGY SURFACE=',I5/ &
     ' NUMBER OF RELAXATION TERMS FOR H+DIAMOND INTERACTIONS=',I5/ &
     ' NUMBER OF H---H NON-BONDED INTERACTION TERMS=',I5/ &
     ' NUMBER OF ATOMS TREATED BY M.O. CALCULATIONS=',I5/ &
     ' CHOICE OF CRCO6 POTENTIAL=', I5//)
  821 FORMAT(6X,20I4)
  826 FORMAT('   ACTIVATE WITH A BOLTZMANN DISTRIBUTION OF VIBRATIONAL', &
     ' ENERGIES')
  827 FORMAT(28X,'NSELT = ',I2,4X,'NSURF = ',I2/)
  828 FORMAT('   MASSES OF ATOMS:',I5,' ATOMS'/)
  835 FORMAT('   ACTIVATE WITH CLASSICAL MICROCANONICAL NORMAL', &
     ' MODE SAMPLING'/)
  836 FORMAT(7X,F10.6,4X,F10.6,4X,F10.6)
  837 FORMAT('   VENUS CALCULATES AND SETS VZERO AND DELH VALUES SO THAT', &
     ' THE POTENTIAL ENERGY OF THE REACTANT(S) IS ZERO')
  838 FORMAT('   VZERO=',F14.3,' KCAL/MOLE')
  839 FORMAT('   PARAMETERS FOR PATH ',I1,':'/5X,'RMAX=',F6.2,'  RBAR=', &
     F6.2,'  NATOMA=',I3,'  NATOMB=',I3,'  DELH=',F8.2)
  840 FORMAT(5X,'INDICES FOR ATOMS OF FRAGMENT A:')
  841 FORMAT(5X,'INDICES FOR ATOMS OF FRAGMENT B:')
  842 FORMAT(5X,'DISTANCE BETWEEN THESE ATOMS DEFINES R.C. :',2I4,2X, &
     'AND',2X,2I4)
  843 FORMAT('   NUMBER OF ADDITIONAL REACTION PATHS=',I2/)
  844 FORMAT('   PARAMETERS FOR REACTANT A')
  845 FORMAT(5X,'NORMAL MODE QUANTUM NUMBERS:')
  846 FORMAT('     REACTANT IS A DIATOMIC (TREATED SEMICLASSICALLY)')
  847 FORMAT(5X,'TRV:',F10.3,5X,'TROT:',F10.3,5X,'AI:',F10.3,5X,'N:',I3, &
     8X,'J:',I3)
  848 FORMAT('RMAX=',F8.2,4X,'RBAR=',F8.2,4X,'DELH=',F8.3)
  851 FORMAT('   ACTIVATE WITH MICROCANONICAL QUASICLASSICAL NORMAL ', &
     'MODE SAMPLING '/'      EBAR=',1PE13.5)
  852 FORMAT('   REACTANT B CANNOT BE INITIALIZED AT A BARRIER')
  853 FORMAT(5X,'L-ATOM=',I3,4X,'K-ATOM=',I3)
  855 FORMAT(/,15X,'VIBRATIONAL ANGULAR MOMENTUM',/,3X, &
     'SECOND MODE OF DEGENERATE PAIR',5X,'PHASE DIFFERENCE')
  856 FORMAT(I16,20X,F6.1)
  894 FORMAT(/5X,'COORDINATES AND MOMENTA ARE READ IN FROM ', &
     'CHECKPOINT FILE'/5X,'READING FROM UNIT 50'/)
  895 FORMAT(/'   CONTINUING RANDOM NUMBER SEQUENCE - READING FROM', &
     ' UNIT 50'/)
  896 FORMAT(5X,'SAVE RANDOM NUMBER SEQUENCE IN UNIT 50', &
     /,5X,12HNEXT SEED IS,1X,8I4,/)
  897 FORMAT(5X,47HCALCULATIONS ARE RESTARTED FROM CHECKPOINT FILE, &
     /5X,'READING FROM UNIT 50'//)
  900 FORMAT(5X,'MOMENTS OF INERTIA IX, IY AND IZ(AMU-A**2):',3F9.3)
  901 FORMAT('   PARAMETERS FOR REACTANT B')
  902 FORMAT(5X,'RELATIVE ENERGY(KCAL):',F7.2,5X, &
     'INITIAL SEPARATION(A):',F6.2)
  903 FORMAT(5X,'NOB=',I2,5X,'BMAX(A)=',F5.1)
  904 FORMAT(5X,'NROT=',I2,5X,'TROT=',F9.2)
  906 FORMAT(5X,'EQUILIBRIUM COORDINATES FOR A:')
  907 FORMAT(5X,'EQUILIBRIUM COORDINATES FOR B:')
  909 FORMAT(6X,25I4)
  910 FORMAT(10X,'REACTION OCCURRED FOR PATH',I3)
  914 FORMAT('   ACTIVATE WITH NORMAL MODE SAMPLING'/)
  916 FORMAT('   ACTIVATE WITH QUANTUM MICROCANONICAL NORMAL MODE', &
     ' SAMPLING'/)
  917 FORMAT('   ACTIVATE CI WITH QUANTUM MICROCANONICAL SAMPLING AND', &
     ' CONE TREATED CLASSICALLY'/)
  918 FORMAT(A)
  919 FORMAT(1X,A)
  920 FORMAT(5X,'HSCALE=',F8.3,'  PSCALE=',F5.2)
  921 FORMAT(5X,'HSCALE=',F8.3)
  922 FORMAT(5X,'ENERGY AT CI=',F8.3)
  925 FORMAT('   HINC=',F9.6,5X,'NPTS=',I2//)
  955 FORMAT('   ACTIVATE WITH LOCAL MODE EXCITATION'/)
  956 FORMAT('   PARAMETERS FOR LOCAL MODE EXCITATION'/6X,I3, &
     ' THE MORSE STRETCH TO BE EXCITED INTO N =',I2,' LEVEL'/6X, &
     'LOCAL MODE ENERGY =',F8.3,6X,' EDELTA(BOXING)=',F6.3/)
  957 FORMAT(10X,'BOXING (N+1) AND (N-1) LEVELS'/ &
     10X,'E(N+1)=',F10.3,'  EDEL=',F8.3/ &
     10X,'E(N-1)=',F10.3,'  EDEL=',F8.3)
  958 FORMAT(5X,'MPLOT=',I2,'  NPLOT=',I3,'  NLM=',I3,'  MNTR=',I3)
  960 FORMAT(5X,'IJDIR=',I2,' IDIR=',I3,' JDIR=',I3)
  961 FORMAT('   THE INDICES FOR THE MORSE OSCILLATORS ', &
     'TO BE MONITORED ',4I3)
  962 FORMAT(5X,'VIBRATIONAL TEMPERATURE=',F7.1)
  963 FORMAT(/5X,'DIATOM CANNOT HAVE BARRIER EXCITATION!'/)
  964 FORMAT(12X,'INITIATE TRAJECTORIES AT A POTENTIAL BARRIER'/)
  965 FORMAT(17X,'REACTION COORDINATE HAS FIXED ENERGY  ',F7.3, &
     '  KCAL/MOL'//)
  966 FORMAT(17X,'REACTION COORDINATE HAS FIXED TEMPERATURE ',F7.1, &
     '  K'//)
  967 FORMAT(/8X,'ORTHANT SAMPLING CANNOT HAVE BARRIER EXCITATION!!!'/)
  968 FORMAT(14X,'INITIAL CONDITIONS ARE CHOSEN RANDOMLY'/)
  969 FORMAT(23X,'MINIMUM ENERGY SEARCH'//)
  970 FORMAT(19X,'INITIAL CONDITIONS ARE READ IN'//)
  971 FORMAT(23X,'NORMAL MODE ANALYSIS'//)
  972 FORMAT(22X,'REACTION PATH FOLLOWING'/)
  973 FORMAT(22X,'TRAJECTORY CALCULATIONS')
  974 FORMAT(/'   ACTIVATE WITH A BOLTZMAN DISTRIBUTION OF ', &
     'TRANSLATIONAL ENERGIES'/5X, &
     'TRANSLATIONAL TEMPERATURE = ',F10.1,' KELVIN', &
     /5X,'INITIAL SEPARATION (A) = : ', F6.2/)
  975 FORMAT(17X,'QUASICLASSICAL MICROCANONICAL ENSEMBLE AT BARRIER'//)
  982 FORMAT(/'   SURFACE PARAMETERS:'/5X, &
     'ATOMS DEFINING THE REFERENCE PLANE IN THE SURFACE:',3I5)
  984 FORMAT(5X,'THTA=',F7.3, &
     '  NCHI=',I2,'  CHI=',F7.3,' (DEG)')
  986 FORMAT(5X,'DISTANCES DEFINING THE ORIGIN:', &
     ' RX0=',F8.3,' RY0=',F8.3,' RZ0=',F8.3,' (A)'/ &
     5X,'NTHET=',I2,'  THET=',F7.3,'   NPHI1=',I2,'  PHI1=', &
     F7.3,'   NPHI2=',I2,'  PHI2=',F7.3,' (DEG)')
 1000 FORMAT(/'   VZERO SET TO ',6F15.8,' KCAL/MOL SO THE REACTANT(S)', &
     'POTENTIAL ENERGY EQUALS ZERO'//)

      RHO = 1.0D-8
    ! Load input file and apply MODEL preset
    call load_input('input_qct.txt')

    if (has_keyword('MODEL')) then
       call apply_preset(trim(get_str('MODEL','')), ios, errmsg)
       if (ios /= 0) then
          write(6,*) trim(errmsg)
          stop
       end if
    end if

    ! Safety check: reject old sequential format
    if (.not. is_keyword_format()) then
       write(6,*) 'ERROR: Old sequential input format is no longer supported.'
       write(6,*) 'Please convert your input to keyword=value format.'
       stop
    end if

    ! ********************************************************************
    ! === NEW KEYWORD-BASED INPUT (VASP INCAR style) ===
    ! ********************************************************************
    TITLE1 = get_str('TITLE1', ' ')
    TITLE2 = get_str('TITLE2', ' ')
    WRITE(6,*) TITLE1
    WRITE(6,*) TITLE2

    if (has_keyword('ELEC_METHOD')) then
       call map_elec_method(trim(get_str('ELEC_METHOD','ADIABATIC')), CALTYP, ios, errmsg)
       if (ios /= 0) then
          write(6,*) trim(errmsg)
          stop
       end if
    else
       CALTYP = get_int('CALTYP', -1)
    end if
    NEDC   = get_int('NEDC', 0)

    NATOMS   = get_int('NATOMS', 0)
    NFC      = get_int('NFC', 0)
    NGLO     = get_int('NGLO', 0)

    call allocate_venus_data(NATOMS)

    GWRITE_LEVEL = get_int('GWRITE_LEVEL', 2)

    NMA=0
    NATOMA(1)=0
    NATOMB(1)=0
    NT=1
    NPATHS=0
    NABJ(1)=1
    NABK(1)=1
    RMAX(1)=100.D0
    RBAR(1)=100.D0
    NERTEST=0
    NSID=0
    NDELH=0
    MPATH=0
    WRITE(6,814) NATOMS, NFC, NGLO

    ! Masses: read ATOM_MASSES; if only 2 values given, broadcast W(2)
    call get_real_arr('ATOM_MASSES', W, NATOMS)
    if (NATOMS > 2) then
       if (count_keyword_values('ATOM_MASSES') <= 2) then
          do I = 3, NATOMS
             W(I) = W(2)
          end do
       end if
    end if
    if (NGLO /= 0) then
       call get_real_arr('WS1', WS1, 3)
       call get_real_arr('WG1', WG1, 3)
       FCG = get_real('FCG', 0.0d0)
       ! GLO layout has NATOMB=0, so the NACTB==5 branch never runs:
       ! read the oscillator sampling temperature here instead.
       TVIBB = get_real('TVIB_B', 0.0d0)
       WRITE(6,*)'GENERALIZED LANGEVIN OSCILLATOR MODEL'
       WRITE(6,*)'FREQUENCY AND FRICTION PARAMETERS IN AU'
       WRITE(6,'(7F10.5)') WS1(1:3), WG1(1:3), FCG
       WRITE(6,*)
       CALL GLOINIT
    else
       FCG = get_real('FCG', 0.0d0)
    end if
    I3N = 3*NATOMS

    NVZERO = 1
    VZERO  = 0.0D0
    WRITE(6,837)

    if (has_keyword('TASK')) then
       call map_task(trim(get_str('TASK','TRAJECTORY')), NSELT, ios, errmsg)
       if (ios /= 0) then
          write(6,*) trim(errmsg)
          stop
       end if
    else
       NSELT = get_int('NSELT', 0)
       ! NSELT 合法域 {-1,0,2,3,4}；静止点(-3)/反应路径(-2)/MIN-ENERGY(1) 已移除；
       ! -1（NORMAL-MODE）为死代码入口，读入即强制终止
       if (NSELT /= -1 .and. NSELT /= 0 .and. NSELT /= 2 &
           .and. NSELT /= 3 .and. NSELT /= 4) then
          write(6,*) 'ERROR: NSELT must be -1, 0, 2, 3, or 4 (got ', NSELT, &
                     '); stationary-point/reaction-path/min-energy removed'
          stop
       end if
    end if
    if (has_keyword('SURFACE_MODEL')) then
       call map_surface_model(trim(get_str('SURFACE_MODEL','NONE')), NSURF, ios, errmsg)
       if (ios /= 0) then
          write(6,*) trim(errmsg)
          stop
       end if
    else
       NSURF = get_int('NSURF', 0)
       NSURF = map_old_nsurf(NSURF)
    end if
    NTHTA  = get_int('NTHTA', -1)
    NRNDXY = get_int('NRNDXY', 1)
    RND_BOX = get_real('RND_BOX', 2.0d0)
    if (NRNDXY == 0) then
       WRITE(6,*)'MANUALLY GIVE XY COORDINATES'
       RX0 = get_real('RX0', 0.0d0)
       RY0 = get_real('RY0', 0.0d0)
    end if
    IF (NSURF.EQ.0) NTHTA=-1
    IF (NSURF.GT.0.AND.NTHTA.GE.0) THEN
       WRITE(6,*)'INITIAL ORIENTATION OF REACTANT ARE SPECIFIED'
    ELSE
       WRITE(6,*)'INITIAL ORIENTATION OF REACTANT ARE RANDOMLY SAMPLED'
    ENDIF
    IF (NSURF.GT.0 .AND. NSELT.NE.2 .AND. NSELT.NE.4) THEN
       WRITE(6,*)'OPTION NOT AVAILABLE, CHANGE NSURF OR NSELT'
       STOP
    ENDIF

    WRITE(6,*)
    WRITE(6,828) NATOMS
    WRITE(6,810) (W(I),I=1,NATOMS)
    WRITE(6,*)
    WRITE(6,827) NSELT, NSURF

    if (has_keyword('INTEGRATOR')) then
       call map_integrator(trim(get_str('INTEGRATOR','VERLET')), INTEGRATOR, LLL, ios, errmsg)
       if (ios /= 0) then
          write(6,*) trim(errmsg)
          stop
       end if
    else
       INTEGRATOR = get_int('INTEGRATOR', 3)
       LLL        = get_int('LLL', 0)
       if (INTEGRATOR == 0) then
          write(6,*) 'ERROR: INTEGRATOR=0 (ADAMS-MOULTON) has been removed; &
               &use VERLET/BEEMAN/RK4/RADAU/SYMPLECTIC'
          stop
       end if
    end if
    NITER      = get_int('NITER', 0)
    NZDOWN     = get_int('NZDOWN', 0)
    IF (NGLO.NE.0) THEN
       INTEGRATOR=3
       NITER=0
       WRITE(6,'(/A)')' FOR THE GLO MODEL, ONLY THE VERLET/BEEMAN', &
          ' ALGORITHM IS AVAILABLE'
    ENDIF
    IF (NFC.NE.0) THEN
       INTEGRATOR=3
       NITER=0
       WRITE(6,'(/A)')' FOR THE FRICTION MODEL, ONLY THE VERLET/BEEMAN', &
          ' ALGORITHM IS AVAILABLE'
    ENDIF
    WRITE(6,*)'INTEGRATOR=',INTEGRATOR,'  LLL=',LLL

    IF (INTEGRATOR.EQ.1) THEN
       WRITE(6,'(/A)')' GAUSS-RADAU INTEGRATION'
       IF(LLL.EQ.0) THEN
          WRITE (6,'(A)')'   LLL CAN NOT BE ZERO'
          STOP
       ELSEIF(LLL.LT.0) THEN
          WRITE(6,'(A)')'   FIXED TIME STEP'
       ELSEIF(LLL.GT.0) THEN
          WRITE(6,'(A)')'   VARIABLE TIME STEP'
          WRITE(6,'(A,I3)')'   CONVERGING DIGITS ARE ',LLL
       ENDIF
    ELSEIF (INTEGRATOR.EQ.2) THEN
       IF (LLL.EQ.4.OR.LLL.EQ.6.OR.LLL.EQ.8)THEN
          WRITE(6,'(/I3,A)')LLL,'TH ORDER SYMPLECTIC INTEGRATION'
       ELSE
          WRITE(6,'(A)')' *** WRONG ORDER FOR SYMPLECTIC INTEGRATION'
          WRITE(6,'(A)')'     LLL MUST BE 4, 6, OR 8***'
          STOP
       ENDIF
    ELSEIF (INTEGRATOR.EQ.3) THEN
       IF(LLL.EQ.1) THEN
          WRITE(6,'(/A)')' VELOCITY VERLET INTEGRATION'
       ELSE
          WRITE(6,'(/A)')" BEEMAN'S THIRD ORDER INTEGRATION"
       ENDIF
    ENDIF

    NT    = get_int('NT', 1)
    NS    = get_int('NS', 100000)
    NIP   = get_int('NIP', 10)
    NCROT = get_int('NCROT', 1)
    NCVIB = get_int('NCVIB', 0)
    IF (NCROT.LE.1) NCROT=1
    TIME  = get_real('DT', 0.01d0)
    SYBTI = 10.0D0*TIME*DBLE(NIP)

    HINC = get_real('HINC', 0.001D0)
    NPTS = get_int('NPTS', 2)
    IF (NSELT.GT.0) WRITE(6,925) HINC, NPTS

    IF (NSELT.EQ.2.OR.NSELT.EQ.3) THEN
       if (has_keyword('INIT_SAMPLING_A')) then
          call map_init_sampling(trim(get_str('INIT_SAMPLING_A','MB')), NACTA, ios, errmsg)
          if (ios /= 0) then
             write(6,*) trim(errmsg)
             stop
          end if
       else
          NACTA = get_int('NACTA', 0)
          ! 整数回退同封闭：4=LOCAL-MODE（F20）/9=CI-QM-MICRO（F23）已移除，
          ! 读入即显式终止（照 NSELT 域校验先例）
          if (NACTA == 4 .or. NACTA == 9) then
             write(6,*) 'ERROR: NACTA=', NACTA, &
                '（LOCAL-MODE/CI-QM-MICRO 分支已移除（F20/F23），不接受该关键字）；'// &
                'use INIT_SAMPLING_A keyword instead'
             stop
          end if
       end if
       if (has_keyword('INIT_SAMPLING_B')) then
          call map_init_sampling(trim(get_str('INIT_SAMPLING_B','MB')), NACTB, ios, errmsg)
          if (ios /= 0) then
             write(6,*) trim(errmsg)
             stop
          end if
       else
          NACTB = get_int('NACTB', 0)
          if (NACTB == 4 .or. NACTB == 9) then
             write(6,*) 'ERROR: NACTB=', NACTB, &
                '（LOCAL-MODE/CI-QM-MICRO 分支已移除（F20/F23），不接受该关键字）；'// &
                'use INIT_SAMPLING_B keyword instead'
             stop
          end if
       end if
       ISEED = get_int('ISEED', 0)
    END IF

    IF (ISEED.LT.0) THEN
       ISEED=-ISEED
    ELSE
       CALL RANDOM_NUMBER(RSEED)
       ISEED=RSEED*ISEED
    ENDIF
    IF (MOD(ISEED,2).EQ.0) ISEED=ISEED+1

    IF (NSELT.EQ.3.AND.NACTA.EQ.1) THEN
       WRITE(6,967)
       STOP
    ENDIF

    IF (NSELT.EQ.3) THEN
       NBAR = get_int('NBAR', 0)
       EBAR = get_real('EBAR', 0.0d0)
       TBAR = get_real('TBAR', 0.0d0)
       WRITE(6,964)
       IF (NBAR.EQ.0) WRITE(6,975)
       IF (NBAR.EQ.1) WRITE(6,965) EBAR
       IF (NBAR.EQ.2) WRITE(6,966) TBAR
       IJDIR = get_int('IJDIR', 0)
       IDIR  = get_int('IDIR', 0)
       JDIR  = get_int('JDIR', 0)
       WRITE(6,960) IJDIR,IDIR,JDIR
    ENDIF
    IF (NSELT.EQ.2) THEN
       WRITE(6,968)
       IF (NSURF.EQ.0) WRITE(6,'(/)')
       IF (NSURF.EQ.1) THEN
          WRITE(6,'(19X,"GAS/SURFACE COLLISION --- RELAXED SURFACE")')
       ELSEIF (NSURF.EQ.2) THEN
          WRITE(6,'(19X,"GAS/SURFACE COLLISION --- RIGID SURFACE")')
       ENDIF
    ENDIF
    IF (NSELT.EQ.1) WRITE(6,969)
    IF (NSELT.EQ.0) WRITE(6,970)

    WRITE(6,811) NT,NS,NIP,NCROT
    WRITE(6,812) ISEED,TIME

    ! --- REACTANT A ---
    NATOMA(1) = get_int('NATOMA', 0)
    NLINA     = get_int('NLINA', 1)

    CALL POTPRE
    IF (NATOMA(1).LE.2) NLINA=1

    IF (NSELT.EQ.3)THEN
       IF(NATOMA(1).EQ.2) THEN
          WRITE(6,963)
          STOP
       ELSEIF(NATOMA(1).LT.2) THEN
          WRITE(6,*)"ERROR INPUT FOR BARRIER EXCITATION!!!"
          STOP
       ENDIF
    ENDIF

    IF (NACTA.EQ.1) WRITE(6,813)
    IF (NACTA.EQ.2) WRITE(6,835)
    IF (NACTA.EQ.3) WRITE(6,914)
    IF (NACTA.EQ.4) WRITE(6,955)
    IF (NACTA.EQ.5) WRITE(6,826)
    IF (NACTA.EQ.6) WRITE(6,*) ' FIXED ENERGY SAMPLING INCLUDING RC'
    IF (NACTA.EQ.7) THEN
       WRITE(6,*)' MD SAMPLING FOR REACTANT A IS NOT SUPPORTED'
       STOP
    ENDIF
    IF (NACTA.EQ.8) WRITE(6,916)
    IF (NACTA.EQ.9) WRITE(6,917)

    K=NATOMA(1)
    WTA(1)=0.0D0
    DO J=1,K
       LA(1,J)=J
       LL(J)=LA(1,J)
       WTA(1)=WTA(1)+W(J)
    ENDDO
    K=3*K
    call get_real_arr('QZA_EQ', QZA(1,1:K), K)
    WRITE(6,906)
    WRITE(6,836)(QZA(1,J),J=1,K)


    ! Branch: single atom, diatomic, or polyatomic?
    IF (NATOMA(1) > 1) THEN
       ! --- Coordinate transform to center of mass / principal axes ---
       DUM=0.0D0
       DO J=1,K/3
          DUM1=QZA(1,3*J)-QZA(1,3)
          DUM=DUM+DUM1*DUM1
       ENDDO
       IF (DUM.LE.1.0D-10) THEN
          DO J=1,K/3
             QZDUM=QZA(1,3*J)
             QZA(1,3*J)=QZA(1,3*J-1)
             QZA(1,3*J-1)=QZA(1,3*J-2)
             QZA(1,3*J-2)=QZDUM
          ENDDO
       ENDIF

       DO J=1,K
          Q(J)=QZA(1,J)
       ENDDO
       WT=WTA(1)
       N=NATOMA(1)

       CALL CENMAS(WT,QCM,VCM,N)

       CALL ROTN(AM,EROT,N)
       CA3(1,1)=AIXX
       CA3(2,1)=-AIXY
       CA3(2,2)=AIYY
       CA3(3,1)=-AIXZ
       CA3(3,2)=-AIYZ
       CA3(3,3)=AIZZ
       CALL EIGN(CA3,CB3,3,RHO)
       AI(1)=EIG(1)
       AI(2)=EIG(2)
       AI(3)=EIG(3)

       ! Transpose matrix CB
       DO II=1,3
          DO JJ=1,II
             DUM=CB3(II,JJ)
             CB3(II,JJ)=CB3(JJ,II)
             CB3(JJ,II)=DUM
          END DO
       END DO

       ! Set phase of the eigenvector
       DO I=1,3
          CBMAX=0.0D0
          DO J=1,3
             IF (ABS(CB3(J,I)).GE.ABS(CBMAX)) CBMAX=CB3(J,I)
          ENDDO
          IF (CBMAX.LT.0.0D0) THEN
             DO J=1,3
                CB3(J,I)=-CB3(J,I)
             ENDDO
          ENDIF
       ENDDO

       ! Check for right-handed cartesian system
       DETCB1=CB3(1,1)*(CB3(2,2)*CB3(3,3)-CB3(2,3)*CB3(3,2))
       DETCB2=CB3(1,2)*(CB3(2,3)*CB3(3,1)-CB3(2,1)*CB3(3,3))
       DETCB3=CB3(1,3)*(CB3(2,1)*CB3(3,2)-CB3(2,2)*CB3(3,1))
       DETCB=DETCB1+DETCB2+DETCB3
       IF (DETCB.LT.0.0) THEN
          DO J=1,3
             CB3(J,1)=-CB3(J,1)
          ENDDO
       END IF

       if (NZDOWN /= 1) then
          DO J=1,K/3
             DO JL=1,3
                DUM=0.0D0
                DO IL=1,3
                   DUM=DUM+CB3(IL,JL)*QQ(3*J+IL-3)
                ENDDO
                QZA(1,3*J+JL-3)=DUM
             ENDDO
          ENDDO
       end if


       ! --- Diatomic or polyatomic parameters ---
       IF (NATOMA(1) == 2) THEN
          ! Diatomic: EBK quantization
          TRVA  = get_real('TRV_A', 0.0d0)
          TROTA = get_real('TROT_A', 0.0d0)
          AIA   = get_real('AI_A', 0.0d0)
          NNA   = get_int('NN_A', 0)
          JA    = get_int('J_A', 0)
          RRA2=(Q(1)-Q(4))**2+(Q(2)-Q(5))**2+(Q(3)-Q(6))**2
          AIA=W(1)*W(2)*RRA2/WTA(1)
          WRITE(6,846)
          WRITE(6,847)TRVA,TROTA,AIA,NNA,JA
          if (NGLO == 0 .and. FCG > 0.0d0) then
             FCG = FCG * 413.4d0
             GSW_ghost = dsqrt(4.0d0*0.00198717d0*C1*TRVA*FCG*W(2)/TIME)
             WRITE(6,*)'GHOST LANGEVIN OSCILLATOR (GLO-SURFACE)'
             WRITE(6,'(A,F12.5,A,F12.5)')' FCG =',FCG,' (code)  GSW =',GSW_ghost
          end if
       ELSE
          ! Polyatomic (NATOMA(1) > 2)
          IF (NACTA.EQ.1.OR.NACTA.EQ.2.OR.NACTA.EQ.6.OR.NACTA.EQ.8.OR.NACTA.EQ.9) THEN
             ENMTA  = get_real('ENMT_A', 0.0d0)
             PSCALA = get_real('PSCALE_A', 0.0d0)
             IF (NACTA.EQ.1) WRITE(6,920)ENMTA,PSCALA
             IF (NACTA.EQ.2.OR.NACTA.EQ.6.OR.NACTA.EQ.8) WRITE(6,921)ENMTA
             IF (NACTA.EQ.9) WRITE(6,922)ENMTA
          ELSEIF (NACTA.EQ.3 .OR. NACTA.EQ.4) THEN
             J=K-6+NLINA
             IF (NSELT.NE.3) THEN
                call get_real_arr('ANQ_A', ANQA, J)
                WRITE(6,845)
                WRITE(6,809)(ANQA(I),I=1,J)
             ELSE
                call get_real_arr('ANQ_A', ANQA, J-1)
                WRITE(6,845)
                WRITE(6,809)(ANQA(I),I=1,J-1)
             ENDIF
             NPHASA = get_int('NPHAS_A', 0)
             IF (NPHASA.GT.0) THEN
                call get_int_arr('IPHAS_A', IPHASA, NPHASA)
                call get_real_arr('PHAS_A', PHASEA, NPHASA)
                WRITE(6,855)
                DO I = 1, NPHASA
                   WRITE(6,856)IPHASA(I),PHASEA(I)
                   PHASEA(I) = PHASEA(I)*C4
                ENDDO
             ENDIF
          ELSEIF(NACTA.EQ.5.OR.NACTA.EQ.0) THEN
             TVIBA = get_real('TVIB_A', 0.0d0)
             WRITE(6,962)TVIBA
          ENDIF

          ! NMA init (polyatomic only)
          NMA=K-6+NLINA
          DO I=1,NMA
             DO J=1,2000
                ENSAV(I,J)=0.0D0
                ENSQ(I,J)=0.0D0
             ENDDO
          ENDDO

          ! NROTA reading
          NROTA = get_int('NROT_A', 1)
          TROTA = get_real('TROT_A', 0.0d0)
          JROTA = get_int('JROT_A', 0)
          KROTA = get_int('KROT_A', 0)
          WRITE(6,904)NROTA,TROTA
          WRITE(6,900)(AI(I),I=1,3)
          IF (NROTA.EQ.2) THEN
             WRITE(6,*)"NROTA=2, SAMPLE ROTATIONAL ENERGY AS A SYSMETRIC TOP"
             WRITE(6,*)"JROTA AND KROTA =",JROTA,KROTA
          ELSE
             NTHTA=-1
             NBFZ=1
          ENDIF

          ! Local mode excitation
          IF (NACTA.EQ.4) THEN
             WRITE(6,*)
             NEXM = get_int('NEXM', 0)
             NLEV = get_int('NLEV', 0)
             CALL LMODE(0,ENU,EDELTU,ENL,EDELTL)
             WRITE(6,956)NEXM,NLEV,ENON,EDELTA
             WRITE(6,957)ENU,EDELTU,ENL,EDELTL
          ENDIF
       ENDIF
    ENDIF

    MPLOT = 0

  ! =====================================================================
  ! --- REACTANT B ---
  ! =====================================================================
    if (has_keyword('NATOMB')) NATOMB(1) = get_int('NATOMB', 0)
    NLINB     = get_int('NLINB', 1)
    IF (NATOMB(1).LE.2) NLINB=1

    IF (NSURF.EQ.2) THEN
       IF (NGLO.EQ.0) THEN
          IF (NATOMA(1).NE.NATOMS .AND. NATOMB(1).EQ.0) STOP 'NATOMA(1).NE.NATOMS'
          IF (NATOMB(1).NE.4 .AND. NATOMB(1).LT.10) STOP 'NATOMB(1).NE.4'
       ELSE
          IF(NATOMA(1)+2.NE.NATOMS) THEN
             WRITE(*,*)'NATOMS SHOULD EQUAL TO NATOMA(1)+2 IN GLO MODEL'
             STOP
          ENDIF
       ENDIF
    ENDIF

    IF(NSELT.EQ.3.AND.NATOMB(1).NE.0)THEN
       WRITE(6,*)"ERROR INPUT FOR BARRIER EXCITATION!!!"
       STOP
    ENDIF

    IF (NATOMB(1) > 0) THEN
       WRITE(6,*)
       WRITE(6,901)
       IF (NACTB.EQ.1) WRITE(6,813)
       IF (NACTB.EQ.2) WRITE(6,835)
       IF (NACTB.EQ.3) WRITE(6,914)
       IF (NACTB.EQ.4) WRITE(6,955)
       IF (NACTB.EQ.5) WRITE(6,826)
       IF (NACTB.EQ.6.AND.NSELT.NE.3) THEN
          WRITE(6,*)' MICROCANONICAL QUASICLASSICAL SAMPLING FOR B?'
       ENDIF
       IF (NACTB.EQ.7)THEN
          IF(NSURF.EQ.1)THEN
             WRITE(6,*)' MD SAMPLING FOR SURFACE'
          ELSE
             WRITE(6,*)' MD SAMPLING FOR MOLECULE IS NOT SUPPORTED'
             STOP
          ENDIF
       ENDIF
       IF (NACTB.EQ.8) WRITE(6,916)
       IF (NACTB.EQ.9) THEN
          WRITE(6,*)'   ACTIVATION AT A CI NOT ALLOWED FOR REACTANT B'
          STOP
       ENDIF

       K=NATOMB(1)
       WTB(1)=0.0D0
       DO J=1,K
          M=J+NATOMA(1)
          IF (NSURF.EQ.2.AND.NGLO.NE.0) THEN
             M=J+NATOMA(1)+2
          ENDIF
          LB(1,J)=M
          LL(J)=LB(1,J)
          IF (NSURF.EQ.2) W(M)=1D30
          WTB(1)=WTB(1)+W(M)
       ENDDO
       K=3*K
       WRITE(6,907)

       IF (NSURF.EQ.2) THEN
          ! F24/D2 fix: the rigid-surface slab coordinates come from QZB_EQ
          ! in the absolute Cartesian frame (no COM/principal-axes transform
          ! — the slab defines the lab frame). Previously QZB was only read
          ! for NSURF=0, leaving every B atom at the origin. Must be read
          ! BEFORE the POTENZ(II) call below so that VZERO/DELH reference
          ! the actual slab geometry.
          call get_real_arr('QZB_EQ', QZB(1,1:K), K)
          WRITE(6,836)(QZB(1,J),J=1,K)
          II=3*NATOMA(1)
          DO J=1,K
             Q(J+II)=QZB(1,J)
          ENDDO
       ENDIF

       ! VZERO determination
       IF(NVZERO.EQ.1)THEN
          II=1
          NDELH(II)=0
          CALL POTENZ(II)
          VZERO=-V
          WRITE(6,1000)VZERO
          WRITE(6,*)
       ENDIF

       IF (NSURF.EQ.0) THEN
          ! Transform QZB to CM/principal axes frame
          call get_real_arr('QZB_EQ', QZB(1,1:K), K)
          DUM=0.0D0
          DO J=1,K/3
             DUM1=QZB(1,3*J)-QZB(1,3)
             DUM=DUM+DUM1*DUM1
          ENDDO
          IF (DUM.LE.1.0D-10) THEN
             DO J=1,K/3
                QZDUM=QZB(1,3*J)
                QZB(1,3*J)=QZB(1,3*J-1)
                QZB(1,3*J-1)=QZB(1,3*J-2)
                QZB(1,3*J-2)=QZDUM
             ENDDO
          ENDIF

          I=3*NATOMA(1)
          DO J=1,K
             Q(J+I)=QZB(1,J)
          ENDDO
          WT=WTB(1)
          N=NATOMB(1)
          CALL CENMAS(WT,QCM,VCM,N)
          CALL ROTN(AM,EROT,N)
          CA3(1,1)=AIXX
          CA3(2,1)=-AIXY
          CA3(2,2)=AIYY
          CA3(3,1)=-AIXZ
          CA3(3,2)=-AIYZ
          CA3(3,3)=AIZZ
          CALL EIGN(CA3,CB3,3,RHO)
          BI(1)=EIG(1)
          BI(2)=EIG(2)
          BI(3)=EIG(3)

          DO II=1,3
             DO JJ=1,II
                DUM=CB3(II,JJ)
                CB3(II,JJ)=CB3(JJ,II)
                CB3(JJ,II)=DUM
             END DO
          END DO

          DO ICB=1,3
             CBMAX=0.0D0
             DO J=1,3
                IF (ABS(CB3(J,ICB)).GE.ABS(CBMAX)) CBMAX=CB3(J,ICB)
             ENDDO
             IF (CBMAX.LT.0.0D0) THEN
                DO J=1,3
                   CB3(J,ICB)=-CB3(J,ICB)
                ENDDO
             ENDIF
          ENDDO

          DETCB1=CB3(1,1)*(CB3(2,2)*CB3(3,3)-CB3(2,3)*CB3(3,2))
          DETCB2=CB3(1,2)*(CB3(2,3)*CB3(3,1)-CB3(2,1)*CB3(3,3))
          DETCB3=CB3(1,3)*(CB3(2,1)*CB3(3,2)-CB3(2,2)*CB3(3,1))
          DETCB=DETCB1+DETCB2+DETCB3
          IF (DETCB.LT.0.0) THEN
             DO J=1,3
                CB3(J,1)=-CB3(J,1)
             ENDDO
          ENDIF

          DO J=1,K/3
             DO JL=1,3
                DUM=0.0D0
                DO IL=1,3
                   DUM=DUM+CB3(IL,JL)*QQ(I+3*J+IL-3)
                ENDDO
                QZB(1,3*J+JL-3)=DUM
             ENDDO
          ENDDO
       ENDIF

       ! Diatomic or polyatomic B parameters
       IF (NATOMB(1) > 1) THEN
          IF (NATOMB(1) <= 2) THEN
             ! Diatomic B
             IF (NATOMB(1).NE.1) THEN
                TRVB  = get_real('TRV_B', 0.0d0)
                TROTB = get_real('TROT_B', 0.0d0)
                AIB   = get_real('AI_B', 0.0d0)
                NNB   = get_int('NN_B', 0)
                JB    = get_int('J_B', 0)
                J1=LB(1,1)*3-2
                J2=LB(1,1)*3-1
                J3=LB(1,1)*3
                J4=LB(1,2)*3-2
                J5=LB(1,2)*3-1
                J6=LB(1,2)*3
                RRB2=(Q(J1)-Q(J4))**2+(Q(J2)-Q(J5))**2 &
                   +(Q(J3)-Q(J6))**2
                AIB=W(LB(1,1))*W(LB(1,2))*RRB2/WTB(1)
                WRITE(6,846)
                WRITE(6,847)TRVB,TROTB,AIB,NNB,JB
                WRITE(6,*)
             ENDIF
          ELSE
             ! Polyatomic B
             IF(NACTB.EQ.7)THEN
                THERMOTEMP = get_real('THERMOTEMP', 300.0d0)
                NSCALE     = get_int('NSCALE', 0)
                NEQUAL     = get_int('NEQUAL', 0)
                NSEL       = get_int('NSEL_B', 0)
                NRGD       = get_int('NRGD_B', 0)
                WRITE(6,*)' DOING MD SAMPLING'
                WRITE(6,*)'NSEL =',NSEL,', NSCALE = ',NSCALE, &
                   ', NEQUAL =', NEQUAL, ', THERMOTEMP =',THERMOTEMP, &
                   ', NRGD =',NRGD
                NTHERMB=0
                NRSCL=0
             ELSE
                NSEL=0
                NSCALE=0
                NEQUAL=0
                NRGD=0
             ENDIF

             IF (NACTB.EQ.1.OR.NACTB.EQ.2.OR.NACTB.EQ.8) THEN
                ENMTB  = get_real('ENMT_B', 0.0d0)
                PSCALB = get_real('PSCALE_B', 0.0d0)
                IF(NACTB.EQ.1)WRITE(6,920)ENMTB,PSCALB
                IF(NACTB.EQ.2.OR.NACTB.EQ.8)WRITE(6,921)ENMTB
             ELSEIF (NACTB.EQ.3.OR.NACTB.EQ.4) THEN
                J=K-6+NLINB
                call get_real_arr('ANQ_B', ANQB, J)
                WRITE(6,845)
                WRITE(6,809)(ANQB(I),I=1,J)
                NPHASB = get_int('NPHAS_B', 0)
                IF (NPHASB.GT.0) THEN
                   call get_int_arr('IPHAS_B', IPHASB, NPHASB)
                   call get_real_arr('PHAS_B', PHASEB, NPHASB)
                   WRITE(6,855)
                   DO I = 1, NPHASB
                      WRITE(6,856)IPHASB(I),PHASEB(I)
                      PHASEB(I)=PHASEB(I)*C4
                   ENDDO
                ENDIF
             ELSEIF(NACTB.EQ.5)THEN
                TVIBB = get_real('TVIB_B', 0.0d0)
                WRITE(6,962)TVIBB
             ENDIF

             IF(NSURF.EQ.0)THEN
                NROTB = get_int('NROT_B', 1)
                TROTB = get_real('TROT_B', 0.0d0)
                JROTB = get_int('JROT_B', 0)
                KROTB = get_int('KROT_B', 0)
                WRITE(6,904)NROTB,TROTB
                WRITE(6,900)(BI(I),I=1,3)
                IF (NROTB.EQ.2) THEN
                   WRITE(6,*)"NROTB=2, SAMPLE ROTATIONAL ENERGY AS A", &
                      " SYSMETRIC TOP"
                   WRITE(6,*)"JROTB AND KROTB =",JROTB,KROTB
                ENDIF
             ELSE
                NROTB=1
                TROTB=0.D0
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    RMAX(1)  = get_real('RMAX', 8.0d0)
    RBAR(1)  = get_real('RBAR', 5.0d0)
    NDELH(1) = get_int('NDELH', 0)

    IF (NDELH(1).EQ.0) THEN
       MPATH=1
       IF (NVZERO.EQ.1) THEN
          II=1
          CALL POTENZ(II)
          DELH(1)=V
       ENDIF
    ENDIF

    WRITE(6,848)RMAX(1),RBAR(1),DELH(1)
    WRITE(6,*)

    NREL = get_int('NREL', 1)
    IF (NREL.EQ.1) THEN
       SEREL = get_real_element('EREL', 1, 0.0d0)
       S     = get_real_element('EREL', 2, 6.0d0)
       WRITE(6,902)SEREL,S
       SEREL=SEREL*C1
    ELSEIF (NREL.EQ.0) THEN
       SEREL = get_real_element('EREL', 1, 0.0d0)
       S     = get_real_element('EREL', 2, 6.0d0)
       TRANS=SEREL
       WRITE(6,974)TRANS,S
    ELSEIF (NREL.EQ.2) THEN
       SEREL  = get_real_element('EREL', 1, 0.0d0)
       S      = get_real_element('EREL', 2, 6.0d0)
       ALPHAA = get_real_element('EREL', 3, 0.0d0) / 1.0D4
       TRANS=SEREL
    ENDIF

    NOB  = get_int('NOB', 1)
    BMAX = get_real('BMAX', 0.0d0)
    WRITE(6,903)NOB,BMAX

    IF (NSURF.GT.0) THEN
       THTA = get_real('THTA', 0.0d0)
       NCHI = get_int('NCHI', 0)
       CHI  = get_real('CHI', 0.0d0)
       WRITE(6,984)THTA,NCHI,CHI
       THTA=THTA*C4
       CHI=CHI*C4
       IF (NSURF.EQ.1) THEN
          WRITE(6,*)'PERIODIC BOUNDARY CONDITIONS'
          SKEW  = get_real('SKEW',  90.0d0)
         BOXLX = get_real('BOXLX', 32.5d0)
         BOXLY = get_real('BOXLY', 30.7d0)
         SKEW  = SKEW * C4
         WRITE(6,*)'THE SIZE AND SHAPE OF THE CELL',SKEW,BOXLX,BOXLY
          WRITE(6,*)
       ENDIF
    ENDIF

    WRITE(6,*)
    NPATHS = get_int('NPATHS', 0)
    IF(NPATHS.LT.0)THEN
       NSID=1
       NPATHS=-NPATHS
       WRITE(6,*)' SURFACE DOES NOT REACT AND FRAGMENTATION OF A'
       WRITE(6,*)' ONLY WILL BE CONSIDERED'
    ELSE
       NSID=0
    ENDIF
    WRITE(6,843)NPATHS
    IF (NPATHS.NE.0) THEN
       M=NPATHS+1
       DO I=2,M
          write(kwbuf,'(A,I0,A)') 'PATH_', I, '_RMAX'
          RMAX(I) = get_real(trim(kwbuf), 8.0d0)
          write(kwbuf,'(A,I0,A)') 'PATH_', I, '_RBAR'
          RBAR(I) = get_real(trim(kwbuf), 5.0d0)
          write(kwbuf,'(A,I0,A)') 'PATH_', I, '_NATOMA'
          NATOMA(I) = get_int(trim(kwbuf), 0)
          write(kwbuf,'(A,I0,A)') 'PATH_', I, '_NATOMB'
          NATOMB(I) = get_int(trim(kwbuf), 0)
          write(kwbuf,'(A,I0,A)') 'PATH_', I, '_NDELH'
          NDELH(I) = get_int(trim(kwbuf), 0)

          K=NATOMA(I)
          write(kwbuf,'(A,I0,A)') 'PATH_', I, '_LA'
          call get_int_arr(trim(kwbuf), LA(I,1:K), K)
          WTA(I)=0.0D0
          DO J=1,K
             WTA(I)=WTA(I)+W(LA(I,J))
          ENDDO
          WRITE(6,840)
          WRITE(6,909)(LA(I,J),J=1,K)
          K=3*NATOMA(I)
          write(kwbuf,'(A,I0,A)') 'PATH_', I, '_QZA'
          call get_real_arr(trim(kwbuf), QZA(I,1:K), K)
          WRITE(6,906)
          WRITE(6,836)(QZA(I,J),J=1,K)
          K=NATOMB(I)
          IF (K.NE.0) THEN
             write(kwbuf,'(A,I0,A)') 'PATH_', I, '_LB'
             call get_int_arr(trim(kwbuf), LB(I,1:K), K)
             WTB(I)=0.0D0
             DO J=1,K
                WTB(I)=WTB(I)+W(LB(I,J))
             ENDDO
             WRITE(6,841)
             WRITE(6,909)(LB(I,J),J=1,K)
             K=3*NATOMB(I)
             write(kwbuf,'(A,I0,A)') 'PATH_', I, '_QZB'
             call get_real_arr(trim(kwbuf), QZB(I,1:K), K)
             WRITE(6,907)
             WRITE(6,836)(QZB(I,J),J=1,K)
             WRITE(6,*)
          ENDIF
          IF (NVZERO.EQ.1) THEN
             II=I
             CALL POTENZ(II)
             DELH(I) = V
          ENDIF
          IF (NDELH(I).EQ.0) MPATH=MPATH+1
          WRITE(6,839)I,RMAX(I),RBAR(I),NATOMA(I),NATOMB(I),DELH(I)
       ENDDO
    ENDIF

    NFQP   = get_int('NFQP', 0)
    NCOOR  = get_int('NCOOR', 0)
    CALL PRINFO(NFQP, NCOOR)

  end subroutine read_venus_input

end module venus_input
