!=====================================================================
! consts.f90 - shared conversion and mathematical constants (fixed by decision)
! Design:
!   Minimal shared set; structural upper limits travel with the consumers.
!   e_conv = 0.04184 is the by-decision fixed kcal/mol -> internal-energy
!   conversion (fixed rounding - analytic anchors must use it, not re-derived
!   CODATA values); e2wvn/wvn2e convert kcal/mol <-> cm^-1 (wvn = wavenumber),
!   wvn2e the exact reciprocal.
!=====================================================================
module consts
  implicit none
  private
  public :: e_conv, r_kcal, kb_code, hbar_code, pi, half_pi, two_pi, dtor, e2wvn, wvn2e

  ! Item by item (the pi family is declared first so dtor can be written as pi/180):
  real(8), parameter :: e_conv    = 0.04184d0                ! kcal/mol to internal energy units
  real(8), parameter :: r_kcal    = 1.9872198404d-3          ! gas constant [kcal/mol/K]
  real(8), parameter :: kb_code   = 0.083144d-3              ! Boltzmann constant in internal units
  real(8), parameter :: hbar_code = 0.063508d0               ! hbar in internal units
  real(8), parameter :: pi        = 3.14159265358979323846d0 ! pi
  real(8), parameter :: half_pi   = pi/2.0d0                 ! pi family: pi/2
  real(8), parameter :: two_pi    = 2.0d0*pi                 ! pi family: 2*pi
  real(8), parameter :: dtor      = pi/180.0d0               ! pi/180, degrees to radians
  real(8), parameter :: e2wvn     = 349.755d0                ! kcal/mol <-> cm^-1: forward conversion
  real(8), parameter :: wvn2e     = 1.0d0/e2wvn              ! kcal/mol <-> cm^-1: reciprocal of e2wvn
end module consts
