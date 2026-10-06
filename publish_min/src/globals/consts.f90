!=====================================================================
! consts.f90 - shared conversion and mathematical constants
! Unit system (ruling 2026-10-06; replaces the kcal/mol-anchored set):
!   energy [eV], length [A], mass [amu], time [fs]. Derived units:
!   momenta [amu*A/fs], forces [amu*A/fs^2], angular frequency
!   [rad/fs], action [amu*A^2/fs]. Values frozen from CODATA 2018;
!   each derivation is stated inline and reproduced in
!   logs/constants_freeze.log (awk one-liner recomputation).
! Bridge pair: eV-family energies and dynamical quantities are linked
!   by ev_dyn/dyn_ev alone (single kinematic fold, LAMMPS-metal style):
!   E[eV] = (p*p/2m)*dyn_ev   and   f[amu*A/fs^2] = g[eV/A]*ev_dyn.
! Face rules: hbar_evfs/kb_ev serve energy contexts (results in eV);
!   hbar_act/kb_act serve action contexts (sigma = sqrt(m*kb_act*T),
!   |L| = sqrt(J(J+1))*hbar_act, WKB action integrals); any eV ->
!   dynamical conversion is always a factor ev_dyn; any dynamical ->
!   eV conversion is always a factor dyn_ev.
!=====================================================================
module consts
  implicit none
  private
  public :: dyn_ev, ev_dyn, hbar_evfs, hbar_act, kb_ev, kb_act, pi, half_pi, two_pi, dtor, ev2wvn, wvn2ev

  ! Item by item (dependencies order: pi before dtor, ev_dyn before *_act):
  real(8), parameter :: dyn_ev    = 103.6426965d0            ! 1 amu*A^2/fs^2 in eV = 1.66053906660e-17 J / 1.602176634e-19 J
                                                              ! (the dynamical -> eV output factor)
  real(8), parameter :: ev_dyn    = 1.0d0/dyn_ev             ! 1 eV in amu*A^2/fs^2 (and 1 eV/A in amu*A/fs^2;
                                                              ! the eV -> dynamical input factor)
  real(8), parameter :: hbar_evfs = 0.6582119569d0           ! hbar [eV*fs] = 6.582119569e-16 eV*s / 1e-15 s
  real(8), parameter :: hbar_act  = hbar_evfs*ev_dyn         ! hbar [amu*A^2/fs] (action context)
  real(8), parameter :: kb_ev     = 8.617333262d-5           ! Boltzmann [eV/K] = 1.380649e-23 J/K / eV
  real(8), parameter :: kb_act    = kb_ev*ev_dyn             ! Boltzmann [amu*A^2/fs^2/K] (Maxwell width context)
  real(8), parameter :: pi        = 3.14159265358979323846d0 ! pi
  real(8), parameter :: half_pi   = pi/2.0d0                 ! pi family: pi/2
  real(8), parameter :: two_pi    = 2.0d0*pi                 ! pi family: 2*pi
  real(8), parameter :: dtor      = pi/180.0d0               ! pi/180, degrees to radians
  real(8), parameter :: ev2wvn    = 8065.543937d0            ! eV <-> cm^-1: forward conversion (1 eV = e/(h*c))
  real(8), parameter :: wvn2ev    = 1.0d0/ev2wvn             ! eV <-> cm^-1: reciprocal of ev2wvn
end module consts
