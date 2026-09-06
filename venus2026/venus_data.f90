!***********************************************************************
!  VENUS_DATA — Global trajectory state (replaces ALL COMMON blocks)
!
!  All arrays formerly sized by compile-time NDA/NDA3 are now
!  ALLOCATABLE and sized at runtime by the actual number of atoms.
!  Fixed-size small arrays and scalars remain as-is.
!***********************************************************************
module venus_data
  use venus_params
  implicit none
  public

  ! =====================================================================
  ! Runtime dimensions
  ! =====================================================================
  integer :: natoms = 0
  integer :: i3n   = 0

  ! =====================================================================
  ! Electronic method selector.
  ! CALTYP = -1 (adiabatic classical); non-adiabatic values (>=0) are
  ! rejected by map_elec_method (loud-stop), so CALTYP stays -1 at runtime.
  ! NEDC retained only for input-keyword compatibility.
  ! =====================================================================
  integer :: CALTYP = -1
  integer :: NEDC   = 0

  ! =====================================================================
  ! PRLIST   — printing / trajectory list
  ! =====================================================================
  real(8) :: T, V, H, TIME
  integer :: NTZ, NT
  integer :: ISEED0(8)
  integer :: NC, NX

  ! =====================================================================
  ! INTEGR   — integration timing
  ! =====================================================================
  real(8) :: ATIME
  integer :: NI, NID

  ! =====================================================================
  ! TABLEB   — lookup table
  ! =====================================================================
  real(8), allocatable :: TABLE(:)

  ! =====================================================================
  ! PRFLAG   — print flags
  ! =====================================================================
  integer :: NFQP, NCOOR, NFR, NUMR, NFB, NUMB, NFA, NUMA
  integer :: NFTAU, NUMTAU, NFTET, NUMTET, NFDH, NUMDH, NFHT, NUMHT
  integer :: GWRITE_LEVEL = 2                     ! verbosity: 0=min .. 3=max

  ! =====================================================================
  ! PARRAY   — print arrays
  ! =====================================================================
  integer :: KR(300), JR(300), KB(300), MB(300), IB(300), IA(300)
  integer :: ITAU(300), ITET(300), IDH(300), IHT(300)

  ! =====================================================================
  ! SELTB    — selection / surface
  ! =====================================================================
  real(8), allocatable :: QZ(:)
  integer :: NSELT, NSFLAG, NACTA, NACTB, NLINA, NLINB, NSURF, NRNDXY

  ! =====================================================================
  ! TRANSB   — translation
  ! =====================================================================
  real(8) :: TRANS
  integer :: NREL

  ! =====================================================================
  ! QPDOT    — coordinates, forces (trajectory core)
  ! =====================================================================
  real(8), allocatable, target :: Q(:),Q_old(:)
  real(8), allocatable :: PDOT(:)
  real(8), allocatable :: FCOEF(:,:)

  ! =====================================================================
  ! PQDOT    — momenta, velocities, masses
  ! =====================================================================
  real(8), allocatable :: P(:)
  real(8), allocatable :: QDOT(:)
  real(8), allocatable :: W(:)

  ! =====================================================================
  ! HFIT     — scaling
  ! =====================================================================
  real(8) :: PSCALA, PSCALB, VZERO
  integer :: NVZERO

  ! =====================================================================
  ! ENPRNT   — energy printing (TDHF/NAMD)
  ! =====================================================================
  real(8) :: E0, E_elec, htil, E_elec_init, elec_diag_count
  real(8) :: H_prev = 0.0d0, T_prev = 0.0d0      ! GWRITE prev-step tracking
  real(8) :: E0_initial = 0.0d0                   ! E0 at first GWRITE call
  real(8) :: E_elec_initial = 0.0d0               ! E_elec at first GWRITE call
  real(8) :: E_mf_initial = 0.0d0                 ! E_mf at first GWRITE call

  ! =====================================================================
  ! ENDBG    — debug energies/forces (NAMD/TDHF)
  ! =====================================================================
  real(8) :: DB_ENGY, DB_LNGY
  real(8), allocatable :: DB_FRCE(:)

  ! =====================================================================
  ! PSN2     — potential surface
  ! =====================================================================
  real(8) :: PESN2, GA, RA, RB

  ! =====================================================================
  ! FORCES   — force field dimensions
  ! =====================================================================
  integer :: NFC, NGLO

  ! =====================================================================
  ! CUBEB    — cubic spline / bonding
  ! =====================================================================
  real(8) :: S3(4), DS3(4), CBIC(15,6), ANG1(20,6,4), GN4(20)

  ! =====================================================================
  ! TESTIN   — test input
  ! =====================================================================
  real(8) :: VRELO
  integer :: INTST

  ! =====================================================================
  ! COORS    — internal coordinates
  ! =====================================================================
  real(8), allocatable :: RBOND(:)
  real(8) :: THETA(ND03), ALPHA(ND04), CTAU(ND06)
  real(8) :: GR(ND08,5), TT(ND09,6), DANG(ND13I)

  ! =====================================================================
  ! FRAGB    — fragments
  ! =====================================================================
  real(8), allocatable :: WTA(:), WTB(:)
  integer, allocatable :: LA(:,:), LB(:,:)
  real(8), allocatable :: QZA(:,:), QZB(:,:)
  integer, allocatable :: NATOMA(:), NATOMB(:)

  ! =====================================================================
  ! TESTB    — reaction path / barrier
  ! =====================================================================
  real(8), allocatable :: RMAX(:), RBAR(:)
  integer :: NTEST, NPATHS
  integer, allocatable :: NABJ(:), NABK(:), NABL(:), NABM(:)
  integer :: NPATH, NAST

  ! =====================================================================
  ! TESTSN2  — saddle point
  ! =====================================================================
  real(8) :: GAO
  integer :: NSAD, NCBA, NCAB, IBAR

  ! =====================================================================
  ! FINALB   — final analysis
  ! =====================================================================
  real(8) :: EROTA, EROTB, EA(3), EB(3)
  real(8) :: AMA(4), AMB(4)
  real(8) :: AN, AJ, BN, BJ
  real(8) :: OAM(4)
  real(8) :: EREL, ERELSQ, ETCM, BF, SDA, SDB
  real(8), allocatable :: DELH(:)
  real(8) :: ANG(NDG)
  integer :: NFINAL

  ! =====================================================================
  ! CHEMAC   — chemical activation
  ! =====================================================================
  real(8), allocatable :: WWA(:), CA(:,:)
  real(8) :: AI(3), ENMTA
  real(8), allocatable :: AMPA(:)
  real(8), allocatable :: WWB(:), CB(:,:)
  real(8) :: BI(3), ENMTB
  real(8), allocatable :: AMPB(:)
  real(8) :: SEREL, S, BMAX, TROTA, TROTB
  real(8), allocatable :: ANQA(:), ANQB(:)
  real(8) :: TVIBA, TVIBB
  integer :: NROTA, NROTB, NOB

  ! =====================================================================
  ! NROTAEQ2 — rotational quantum numbers
  ! =====================================================================
  integer :: JROTA, KROTA, JROTB, KROTB

  ! =====================================================================
  ! ALIGN    — alignment
  ! =====================================================================
  integer :: NTHTA, NBFZ

  ! =====================================================================
  ! GLOP     — generalized Langevin oscillator
  ! =====================================================================
  real(8) :: WS1(3), WG1(3), WS2(3), WG2(3), WGS1(3), WGS2(3), WEFF(3)
  real(8) :: GSW, FCG, COEFA, COEFB
  real(8), allocatable :: GN(:)

      ! =====================================================================
      ! GHOST ANCHOR — Langevin ghost (atom 2 friction + noise)
      ! =====================================================================
      real(8) :: GSW_ghost           ! noise strength for ghost anchor
      real(8) :: F_lang_save(3)      ! Langevin force from first half-step

  ! =====================================================================
  ! WASTE    — workspace copies of Q, P
  ! =====================================================================
  real(8), allocatable :: QQ(:), PP(:)
  real(8) :: WX, WY, WZ
  integer, allocatable :: LBOND(:), LL(:)
  integer :: NAM

  ! =====================================================================
  ! DIATB    — diatomic parameters
  ! =====================================================================
  integer :: NNA, JA, NNB, JB

  ! =====================================================================
  ! RSTART   — restart / step control
  ! =====================================================================
  real(8) :: HINC
  integer :: NPTS

  ! =====================================================================
  ! RKUTTA   — Runge-Kutta coefficients
  ! =====================================================================
  real(8) :: RAA1, RA1, RA2, RA3, RB1, RB2, RB3, RC1, RC2

  ! =====================================================================
  ! LMODEB   — local mode excitation
  ! =====================================================================
  real(8) :: ENON, EDELTA, RWANT, PWANT
  integer :: NEXM, NLEV, JFLAG

  ! =====================================================================
  ! RANCOM   — random number state
  ! =====================================================================
  real(8) :: RANLST(100)
  integer :: ISEED3(8), IBFCTR

  ! =====================================================================
  ! GPATHB   — general path
  ! =====================================================================
  real(8), allocatable :: WM(:)
  real(8), allocatable :: TEMP(:), GTEMP(:)
  real(8) :: AI1D(5), AAI(2), BBI(2), SYMM(5)
  real(8) :: SYMA, SYMB
  integer, allocatable :: NFLAG(:)
  integer :: N1DR, N2DR

  ! =====================================================================
  ! SYBB     — system block
  ! =====================================================================
  character(len=80) :: TITLE1, TITLE2
  real(8) :: SYBTI

  ! =====================================================================
  ! SADDLE   — saddle point for barrier excitation
  ! =====================================================================
  real(8) :: EBAR, TBAR, EZERO
  integer :: NBAR, IDIR, JDIR, IJDIR

  ! =====================================================================
  ! INERT    — inertia tensor
  ! =====================================================================
  real(8) :: UXX, UXY, UXZ, UYY, UYZ, UZZ
  real(8) :: AIXX, AIXY, AIXZ, AIYY, AIYZ, AIZZ

  ! =====================================================================
  ! EIGVL    — eigenvalues
  ! =====================================================================
  real(8), allocatable :: EIG(:)

  ! =====================================================================
  ! VECTB    — vector block
  ! =====================================================================
  real(8) :: VI(4), OAMI(4), AMAI(4), AMBI(4)
  real(8) :: ETAI, ERAI, ETBI, ERBI

  ! =====================================================================
  ! GAUSS    — Gaussian job info
  ! =====================================================================
  integer :: ICHRG, IMULTP
  integer :: IAN(NDAyf)
  real(8) :: VGAUSS

  ! =====================================================================
  ! SURFB    — surface / orientation
  ! =====================================================================
  real(8) :: THTA, PHI, CHI, RX0, RY0, RZ0, THET, PHI1, PHI2
  integer :: NN1, NN2, NN3, NN4, NCHI, NTHET, NPHI1, NPHI2

  ! =====================================================================
  ! VMAXB    — velocity maximum
  ! =====================================================================
  real(8), allocatable :: QVMAX(:), PVMAX(:)
  real(8) :: VMAX
  integer :: NCVMAX

  ! =====================================================================
  ! VANGB    — vibrational angle / phase
  ! =====================================================================
  real(8) :: PHASEA(5), PHASEB(5)
  integer :: IPHASA(5), IPHASB(5)
  integer :: NPHASA, NPHASB

  ! =====================================================================
  ! GDD      — direct dynamics
  ! =====================================================================
  real(8) :: DIATM(MXDIATOM,4), GSTOP(MXPTH,6)
  integer, allocatable :: IATMP(:,:)
  integer :: IDIATM(MXDIATOM,2), IATSP(MXPTH,3)
  integer, allocatable :: ISOTOPE(:)
  integer :: NPTH, NDIATOM, NMEM, NDISK, NS, NCHKP, NHESS, NCPU
  integer :: NISOTOPE
  character(len=30) :: METHOD

  ! =====================================================================
  ! VRSCAL   — velocity rescaling / thermostat
  ! =====================================================================
  real(8) :: THERMOTEMP
  integer :: NSEL, NSCALE, NEQUAL, NRGD

  ! =====================================================================
  ! QPSCAL   — coordinate/momentum scaling
  ! =====================================================================
  real(8), allocatable :: QTEMP(:), PTEMP(:)

  ! =====================================================================
  ! THERMOBATH — thermostat bath
  ! =====================================================================
  integer :: NTHERMB, NRSCL
  integer, allocatable :: NTHMID(:)

  ! =====================================================================
  ! TMPNJ    — temperature / diatomic
  ! =====================================================================
  real(8) :: TRVA, AIA, TRVB, AIB

  ! =====================================================================
  ! SFEQUIL  — surface equilibrium
  ! =====================================================================
  integer :: INTEGRATOR, LLL, NITER

  ! =====================================================================
  ! CELL     — periodic cell
  ! =====================================================================
  real(8) :: SKEW, BOXLX, BOXLY, auro, RND_BOX

  ! =====================================================================
  ! HDIAG    — Hamiltonian diagonalization
  ! =====================================================================
  real(8) :: HTMIN, HTMAX, ZASYM, DELM

  ! =====================================================================
  ! STATIS   — statistics
  ! =====================================================================
  integer :: NONREACT, NREACT(100), NLONG

  ! =====================================================================
  ! DENSWIT  — density switching
  ! =====================================================================
  real(8) :: ALPHAA

  ! =====================================================================
  ! unnamed  — NSID, NDELH, MPATH
  ! =====================================================================
  integer :: NSID, MPATH
  integer, allocatable :: NDELH(:)

  ! =====================================================================
  ! HIST_ATOM — history / analysis
  ! =====================================================================
  real(8), allocatable :: AErelS(:), AInitX(:), AInitY(:)
  integer, allocatable :: NBncS(:)
  real(8), allocatable :: AZminS(:), AVfinS(:), AThtaS(:)

  ! =====================================================================
  ! HIST_ATOM_NAD — non-adiabatic history
  ! =====================================================================
  real(8), allocatable :: AEphS(:), AEelS(:), AEu0S(:)
  integer, allocatable :: NCrstS(:)
  real(8) :: E_surf_initial, E_C_initial

  ! =====================================================================
  ! MORSEB   — Morse bond parameters (shared: EBOND, WEBOND, FINLNJ, etc.)
  ! =====================================================================
  real(8) :: RMZ(ND02), B(ND02), D(ND02)
  integer :: N2J(ND02), N2K(ND02)
  real(8) :: CM1(ND02), CM2(ND02), CM3(ND02), CM4(ND02), CM5(ND02)

  ! =====================================================================
  ! FORCE01  — force arrays for non-adiabatic dynamics (NAMDf3.f90)
  ! =====================================================================
  real(8), allocatable :: PDOT0(:), PDOT1(:), PDOT01(:)

  ! =====================================================================
  ! FINALQP  — final coordinates and momenta for normal mode analysis
  ! =====================================================================
  real(8), allocatable :: QFINAL(:), PFINAL(:)

  ! =====================================================================
  ! GAUSS2   — Gaussian Hessian (DVDQ.f90, FMTRX.f)
  ! =====================================================================
  real(8) :: GAUHES(NDIHE)

  ! =====================================================================
  ! WNJ/WNS  — diatomic parameters for JMAXCALC, PROBJ, SELECT
  ! =====================================================================
  real(8) :: WD1, WD2
  real(8) :: WGT(0:500)

  ! =====================================================================
  ! STRETB   — stretch bond parameters (FINLNJ.f)
  ! =====================================================================
  real(8) :: RSZ(ND01), FS(ND01)
  integer :: N1J(ND01), N1K(ND01)

  ! =====================================================================

  ! =====================================================================
  real(8) :: NZDOWN

  ! =====================================================================
  ! FREQFLAG — frequency flag (FMTRX.f)
  ! =====================================================================
  integer :: IARB

  ! =====================================================================
  ! diabatic_energy — shared with NAMDf3.f90, FINAL.f, POTENZ.f
  ! =====================================================================
  real(8) :: U0, U1, V01, Eg, U0_ini, V01_ini

  ! =====================================================================
  ! FRQC — frequency workspace, shared between GPATH.f and MPATH.f
  ! =====================================================================
  real(8), allocatable :: WC(:), WCO(:)
  real(8) :: EZP, EZPO, EZM
  integer :: NRQ, NVQ

  ! =====================================================================
  ! FRQR — internal rotation workspace for MPATH.f
  ! =====================================================================
  real(8), allocatable :: WR1(:), WR2(:)
  real(8) :: TXM, TYM, TZM, TXR1, TYR1, TZR1, TXR2, TYR2, TZR2
  real(8) :: WTM, WTR1, WTR2
  integer :: NCHN, NWM, NWR1, NWR2, NIRM, NIRR1, NIRR2

  ! =====================================================================
  ! FR2 / ARRAYS — local workspace for MPATH.f (normal mode analysis)
  !   Named A_MPATH/B_MPATH to avoid conflict with B(ND02) Morse parameter
  ! =====================================================================
  real(8), allocatable :: DG(:,:), DIM(:)
  real(8), allocatable :: A_MPATH(:,:), DA_MPATH(:), B_MPATH(:,:), DB_MPATH(:)

  ! =====================================================================
  ! Additional variables used locally in VENUS.f90 (not in COMMON)
  ! but needed by subroutines via modules
  ! =====================================================================
  real(8) :: C8_venus        ! local copy, set from venus_params%C8

  ! NDAyf-related — keep for reference
  integer :: NDAYF_val = NDAyf

contains

  ! -------------------------------------------------------------------
  ! Allocate all arrays based on actual number of atoms
  ! Must be called AFTER reading natoms from input.
  ! -------------------------------------------------------------------
  subroutine allocate_venus_data(n)
    integer, intent(in) :: n
    integer :: alloc_err

    natoms = n
    i3n    = 3 * n

    ! Zero-fill scalars that need init
    T = 0.0d0; V = 0.0d0; H = 0.0d0; TIME = 0.0d0
    NTZ = 0; NT = 0; NC = 0; NX = 0
    ISEED0 = 0
    ATIME = 0.0d0; NI = 0; NID = 0

    NFQP = 0; NCOOR = 0; NFR = 0; NUMR = 0; NFB = 0; NUMB = 0
    NFA = 0; NUMA = 0; NFTAU = 0; NUMTAU = 0
    NFTET = 0; NUMTET = 0; NFDH = 0; NUMDH = 0; NFHT = 0; NUMHT = 0
    KR = 0; JR = 0; KB = 0; MB = 0; IB = 0; IA = 0
    ITAU = 0; ITET = 0; IDH = 0; IHT = 0

    NSELT = 0; NSFLAG = 0; NACTA = 0; NACTB = 0
    NLINA = 0; NLINB = 0; NSURF = 0; NRNDXY = 1
    TRANS = 0.0d0; NREL = 0

    PSCALA = 0.0d0; PSCALB = 0.0d0; VZERO = 0.0d0
    NVZERO = 1
    E0 = 0.0d0; E_elec = 0.0d0; htil = 0.0d0
    E_elec_init = 0.0d0; elec_diag_count = 0.0d0
    GWRITE_LEVEL = 2; H_prev = 0.0d0; T_prev = 0.0d0
    E0_initial = 0.0d0; E_elec_initial = 0.0d0; E_mf_initial = 0.0d0
    PESN2 = 0.0d0; GA = 0.0d0; RA = 0.0d0; RB = 0.0d0
    continue  ! NFC, NGLO already set from input

    S3 = 0.0d0; DS3 = 0.0d0; CBIC = 0.0d0; ANG1 = 0.0d0; GN4 = 0.0d0
    VRELO = 0.0d0; INTST = 0

    THETA = 0.0d0; ALPHA = 0.0d0; CTAU = 0.0d0
    GR = 0.0d0; TT = 0.0d0; DANG = 0.0d0

    NTEST = 0; NPATHS = 0; NPATH = 0; NAST = 0
    GAO = 0.0d0; NSAD = 0; NCBA = 0; NCAB = 0; IBAR = 0

    EROTA = 0.0d0; EROTB = 0.0d0; EA = 0.0d0; EB = 0.0d0
    AMA = 0.0d0; AMB = 0.0d0
    AN = 0.0d0; AJ = 0.0d0; BN = 0.0d0; BJ = 0.0d0
    OAM = 0.0d0
    EREL = 0.0d0; ERELSQ = 0.0d0; ETCM = 0.0d0
    BF = 0.0d0; SDA = 0.0d0; SDB = 0.0d0
    ANG = 0.0d0; NFINAL = 0

    AI = 0.0d0; ENMTA = 0.0d0; BI = 0.0d0; ENMTB = 0.0d0
    SEREL = 0.0d0; S = 0.0d0; BMAX = 0.0d0; TROTA = 0.0d0; TROTB = 0.0d0
    TVIBA = 0.0d0; TVIBB = 0.0d0
    NROTA = 0; NROTB = 0; NOB = 0
    JROTA = 0; KROTA = 0; JROTB = 0; KROTB = 0
    NTHTA = 0; NBFZ = 0

    WS1 = 0.0d0; WG1 = 0.0d0; WS2 = 0.0d0; WG2 = 0.0d0
    WGS1 = 0.0d0; WGS2 = 0.0d0; WEFF = 0.0d0
    GSW = 0.0d0; FCG = 0.0d0; COEFA = 0.0d0; COEFB = 0.0d0

    WX = 0.0d0; WY = 0.0d0; WZ = 0.0d0; NAM = 0
    NNA = 0; JA = 0; NNB = 0; JB = 0

    HINC = 0.0d0; NPTS = 0
    RAA1 = 0.0d0; RA1 = 0.0d0; RA2 = 0.0d0; RA3 = 0.0d0
    RB1 = 0.0d0; RB2 = 0.0d0; RB3 = 0.0d0; RC1 = 0.0d0; RC2 = 0.0d0

    ENON = 0.0d0; EDELTA = 0.0d0; RWANT = 0.0d0; PWANT = 0.0d0
    NEXM = 0; NLEV = 0; JFLAG = 0

    RANLST = 0.0d0; ISEED3 = 0; IBFCTR = 0

    AI1D = 0.0d0; AAI = 0.0d0; BBI = 0.0d0; SYMM = 0.0d0
    SYMA = 0.0d0; SYMB = 0.0d0; N1DR = 0; N2DR = 0

    SYBTI = 0.0d0  ! TITLE1, TITLE2 already set from input

    EBAR = 0.0d0; TBAR = 0.0d0; EZERO = 0.0d0
    NBAR = 0; IDIR = 0; JDIR = 0; IJDIR = 0

    UXX = 0.0d0; UXY = 0.0d0; UXZ = 0.0d0
    UYY = 0.0d0; UYZ = 0.0d0; UZZ = 0.0d0
    AIXX = 0.0d0; AIXY = 0.0d0; AIXZ = 0.0d0
    AIYY = 0.0d0; AIYZ = 0.0d0; AIZZ = 0.0d0

    VI = 0.0d0; OAMI = 0.0d0; AMAI = 0.0d0; AMBI = 0.0d0
    ETAI = 0.0d0; ERAI = 0.0d0; ETBI = 0.0d0; ERBI = 0.0d0

    ICHRG = 0; IMULTP = 0; IAN = 0; VGAUSS = 0.0d0

    THTA = 0.0d0; PHI = 0.0d0; CHI = 0.0d0
    RX0 = 0.0d0; RY0 = 0.0d0; RZ0 = 0.0d0
    THET = 0.0d0; PHI1 = 0.0d0; PHI2 = 0.0d0
    NN1 = 0; NN2 = 0; NN3 = 0; NN4 = 0
    NCHI = 0; NTHET = 0; NPHI1 = 0; NPHI2 = 0

    VMAX = 0.0d0; NCVMAX = 0

    PHASEA = 0.0d0; PHASEB = 0.0d0
    IPHASA = 0; IPHASB = 0; NPHASA = 0; NPHASB = 0

    DIATM = 0.0d0; GSTOP = 0.0d0
    IDIATM = 0; IATSP = 0
    NPTH = 0; NDIATOM = 0; NMEM = 0; NDISK = 0
    NS = 0; NCHKP = 0; NHESS = 0; NCPU = 0; NISOTOPE = 0
    METHOD = ' '

    THERMOTEMP = 0.0d0; NSEL = 0; NSCALE = 0; NEQUAL = 0; NRGD = 0
    NTHERMB = 0; NRSCL = 0

    TRVA = 0.0d0; AIA = 0.0d0; TRVB = 0.0d0; AIB = 0.0d0
    INTEGRATOR = 0; LLL = 0; NITER = 0

    SKEW = 0.0d0; BOXLX = 0.0d0; BOXLY = 0.0d0; auro = 2.94998d0
    RND_BOX = 2.0d0
    HTMIN = 0.0d0; HTMAX = 0.0d0
    ! ZASYM, DELM already set before call (ZASYM=10.0D0, DELM=1.0D4)

    NONREACT = 0; NREACT = 0; NLONG = 0
    ALPHAA = 0.0d0

    NSID = 0; MPATH = 0

    U0 = 0.0d0; U1 = 0.0d0; V01 = 0.0d0
    Eg = 0.0d0; U0_ini = 0.0d0; V01_ini = 0.0d0

    RMZ = 0.0d0; B = 0.0d0; D = 0.0d0
    N2J = 0; N2K = 0
    CM1 = 0.0d0; CM2 = 0.0d0; CM3 = 0.0d0; CM4 = 0.0d0; CM5 = 0.0d0

    ! --- Allocatable arrays ---
    allocate(Q(i3n),Q_old(i3n), PDOT(i3n), FCOEF(i3n,i3n), &
             P(i3n), QDOT(i3n), W(n), &
             QQ(i3n), PP(i3n), LBOND(n), LL(n), &
             QZ(i3n), GN(i3n), &
             QVMAX(i3n), PVMAX(i3n), &
             QTEMP(i3n), PTEMP(i3n), NTHMID(n), &
             TABLE(42*n), &
             RBOND(n*(n+1)/2), &
             WTA(NDP), WTB(NDP), &
             LA(NDP,n), LB(NDP,n), &
             QZA(NDP,i3n), QZB(NDP,i3n), &
             NATOMA(NDP), NATOMB(NDP), &
             RMAX(NDP), RBAR(NDP), &
             NABJ(NDP), NABK(NDP), NABL(NDP), NABM(NDP), &
             DELH(NDP), NDELH(NDP), &
             WWA(i3n), CA(NDA3yf,NDA3yf), AMPA(i3n), &
             WWB(i3n), CB(NDA3yf,NDA3yf), AMPB(i3n), &
             ANQA(i3n), ANQB(i3n), &
             EIG(NDA3yf), &
             WM(i3n), TEMP(NDP), GTEMP(NDP), NFLAG(NDP), &
             IATMP(mxpth,n), ISOTOPE(n), &
             AErelS(99999), NBncS(99999), AInitX(99999), AInitY(99999), &
             AZminS(99999), AVfinS(99999), AThtaS(99999), &
             AEphS(99999), AEelS(99999), AEu0S(99999), NCrstS(99999), &
             PDOT0(i3n), PDOT1(i3n), PDOT01(i3n), &
             QFINAL(i3n), PFINAL(i3n), &
             WC(i3n), WCO(i3n), &
             WR1(i3n), WR2(i3n), &
             DG(i3n,2), DIM(i3n), &
             A_MPATH(i3n,i3n), DA_MPATH(i3n), B_MPATH(i3n,i3n), DB_MPATH(i3n), &
             DB_FRCE(i3n), &
             stat=alloc_err)
    if (alloc_err /= 0) then
       write(0,*) 'VENUS_DATA: ERROR: allocation failed for natoms=', n
       stop
    end if

    ! Initialize all allocatable arrays to zero
    Q = 0.0d0; PDOT = 0.0d0; FCOEF = 0.0d0
    P = 0.0d0; QDOT = 0.0d0; W = 0.0d0
    QQ = 0.0d0; PP = 0.0d0; LBOND = 0; LL = 0
    QZ = 0.0d0; GN = 0.0d0
    QVMAX = 0.0d0; PVMAX = 0.0d0
    QTEMP = 0.0d0; PTEMP = 0.0d0; NTHMID = 0
    TABLE = 0.0d0
    RBOND = 0.0d0
    WTA = 0.0d0; WTB = 0.0d0
    LA = 0; LB = 0
    QZA = 0.0d0; QZB = 0.0d0
    NATOMA = 0; NATOMB = 0
    RMAX = 0.0d0; RBAR = 0.0d0
    NABJ = 0; NABK = 0; NABL = 0; NABM = 0
    DELH = 0.0d0; NDELH = 0
    WWA = 0.0d0; CA = 0.0d0; AMPA = 0.0d0
    WWB = 0.0d0; CB = 0.0d0; AMPB = 0.0d0
    ANQA = 0.0d0; ANQB = 0.0d0
    EIG = 0.0d0
    WM = 0.0d0; TEMP = 0.0d0; GTEMP = 0.0d0; NFLAG = 0
    IATMP = 0; ISOTOPE = 0
    AErelS = 0.0d0; NBncS = 0
    AInitX = 0.0d0; AInitY = 0.0d0
    AZminS = 0.0d0; AVfinS = 0.0d0; AThtaS = 0.0d0
    AEphS = 0.0d0; AEelS = 0.0d0; AEu0S = 0.0d0; NCrstS = 0
    E_surf_initial = 0.0d0; E_C_initial = 0.0d0
    PDOT0 = 0.0d0; PDOT1 = 0.0d0; PDOT01 = 0.0d0
    QFINAL = 0.0d0; PFINAL = 0.0d0
    WC = 0.0d0; WCO = 0.0d0
    EZP = 0.0d0; EZPO = 0.0d0; EZM = 0.0d0
    NRQ = 0; NVQ = 0
    WR1 = 0.0d0; WR2 = 0.0d0
    TXM = 0.0d0; TYM = 0.0d0; TZM = 0.0d0
    TXR1 = 0.0d0; TYR1 = 0.0d0; TZR1 = 0.0d0
    TXR2 = 0.0d0; TYR2 = 0.0d0; TZR2 = 0.0d0
    WTM = 0.0d0; WTR1 = 0.0d0; WTR2 = 0.0d0
    NCHN = 0; NWM = 0; NWR1 = 0; NWR2 = 0
    NIRM = 0; NIRR1 = 0; NIRR2 = 0
    DG = 0.0d0; DIM = 0.0d0
    A_MPATH = 0.0d0; DA_MPATH = 0.0d0
    B_MPATH = 0.0d0; DB_MPATH = 0.0d0
    GAUHES = 0.0d0
    DB_FRCE = 0.0d0; DB_ENGY = 0.0d0; DB_LNGY = 0.0d0
    WD1 = 0.0d0; WD2 = 0.0d0; WGT = 0.0d0
    RSZ = 0.0d0; FS = 0.0d0; N1J = 0; N1K = 0
    IARB = 0

  end subroutine allocate_venus_data

end module venus_data
