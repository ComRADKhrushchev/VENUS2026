!***********************************************************************
!  VENUS_PARAMS — Compile-time constants (replaces SIZES + CONSTN)
!
!  This module consolidates all PARAMETER dimensions and physical
!  constants that were formerly in the SIZES include file and the
!  CONSTN COMMON block.
!***********************************************************************
module venus_params
  implicit none
  public

  ! =====================================================================
  ! VENUS dimensions (formerly SIZES)
  !
  ! NDA     — maximum number of atoms (compile-time limit for local arrays)
  ! NDA3    — total degrees of freedom (NDA*3)
  ! NDAyf   — max atoms in fragment A
  ! NDA3yf  — 3*NDAyf
  ! NDP     — number of reaction paths
  ! NDG     — number of angles calculated in FINAL
  ! =====================================================================
  integer, parameter :: maxk    = 50000
  integer, parameter :: NDA     = 1000
  integer, parameter :: NDA3    = NDA * 3
  integer, parameter :: NDP     = 10
  integer, parameter :: NDG     = 20
  integer, parameter :: NDAyf   = 30
  integer, parameter :: NDA3yf  = NDAyf * 3

  ! Potential function term dimensions
  integer, parameter :: ND01   = 6000
  integer, parameter :: ND02   = 12000
  integer, parameter :: ND03   = 10000
  integer, parameter :: ND04   = 20
  integer, parameter :: ND05   = 200000
  integer, parameter :: ND06   = 20
  integer, parameter :: ND07   = 2000000
  integer, parameter :: ND08   = 20
  integer, parameter :: ND09   = 20
  integer, parameter :: ND10   = 500
  integer, parameter :: ND11   = 1200
  integer, parameter :: ND12   = 1200
  integer, parameter :: ND13I  = 4000
  integer, parameter :: ND13J  = 12
  integer, parameter :: ND21   = 2
  integer, parameter :: ND22   = 100

  ! Gauss-Legendre quadrature
  integer, parameter :: NCF    = 300

  ! Hessian storage
  integer, parameter :: NDIHE  = NDA3yf * (NDA3yf + 1) / 2 + 2 * NDA3yf

  ! EBK quantization
  integer, parameter :: ne     = 300
  integer, parameter :: mxlvl  = 600

  ! Direct dynamics
  integer, parameter :: mxpth    = 10
  integer, parameter :: mxdiatom = 5

  ! Energy conversion
  real(8), parameter :: cal2cm = 349.755d0
  real(8), parameter :: cm2cal = 1.0d0 / cal2cm

  ! =====================================================================
  ! Physical / conversion constants (formerly COMMON/CONSTN/)
  ! =====================================================================
  real(8), parameter :: C1     = 0.04184D0        ! kcal/mol → a.u.
  real(8), parameter :: C2     = 6.022045D0
  real(8), parameter :: C3     = 6.022045D0
  real(8), parameter :: C4     = 0.01745329D0     ! degrees → radians
  real(8), parameter :: C5     = 0.083144D-3
  real(8), parameter :: C6     = 1.8836518D-3
  real(8), parameter :: C7     = 0.063508D0
  real(8), parameter :: C8     = 1.9872198404D-3
  real(8), parameter :: PI     = 3.14159265358979323846D0
  real(8), parameter :: HALFPI = PI / 2.0D0
  real(8), parameter :: TWOPI  = 2.0D0 * PI

  ! =====================================================================
  ! Neural Network potential parameters
  !
  ! Per-type architecture is read from NN_weights_all.txt at runtime.
  ! Each potential type declares its own NH=<hidden_count> in the file.
  ! Weights and metadata are allocated by load_all_nn_weights.
  ! =====================================================================
  integer, parameter :: NN_POT_TYPES   = 7
  integer, parameter :: NN_MAX_HIDDEN  = 20   ! sanity bound for NH

  ! Potential type indices for NN weight array
  integer, parameter :: NN_E0   = 1   ! Ground state (Rydberg)
  integer, parameter :: NN_1D   = 2   ! Singlet D state (Morse)
  integer, parameter :: NN_3P   = 3   ! Triplet P state (Morse)
  integer, parameter :: NN_mC   = 4   ! Negative ion state (Morse)
  integer, parameter :: NN_pC   = 5   ! Positive ion state (Morse)
  integer, parameter :: NN_COUP = 6   ! Coupling V_k
  integer, parameter :: NN_Eh   = 7   ! Energy correction for orbital energies

end module venus_params
