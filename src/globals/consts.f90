! consts.f90 - unit contract: the minimal shared set of conversion and mathematical constants
!             (structural upper limits travel with the consumers)
module consts
  implicit none
  private
  public :: e_conv, r_kcal, kb_code, hbar_code, pi, half_pi, two_pi, dtor, e2wvn, wvn2e

  ! Item by item per the §1.1 table (old names in parentheses at line end; the pi family
  ! is declared first so dtor can be written as pi/180):
  real(8), parameter :: e_conv    = 0.04184d0                ! kcal/mol to internal energy units (old name C1)
  real(8), parameter :: r_kcal    = 1.9872198404d-3          ! gas constant [kcal/mol/K] (old name C8)
  real(8), parameter :: kb_code   = 0.083144d-3              ! Boltzmann constant in internal units (old name C5)
  real(8), parameter :: hbar_code = 0.063508d0               ! hbar in internal units (old name C7)
  real(8), parameter :: pi        = 3.14159265358979323846d0 ! pi (old name PI)
  real(8), parameter :: half_pi   = pi/2.0d0                 ! pi family: pi/2 (old name HALFPI)
  real(8), parameter :: two_pi    = 2.0d0*pi                 ! pi family: 2*pi (old name TWOPI)
  real(8), parameter :: dtor      = pi/180.0d0               ! pi/180, degrees to radians (old name C4)
  real(8), parameter :: e2wvn     = 349.755d0                ! kcal/mol <-> cm^-1: forward conversion (old name cal2cm)
  real(8), parameter :: wvn2e     = 1.0d0/e2wvn              ! kcal/mol <-> cm^-1: reciprocal of e2wvn (old name cm2cal)
end module consts
