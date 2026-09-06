module harmonic_sampling
  ! Classical harmonic oscillator Boltzmann energy sampling
  ! for the surface oscillator (Au-Au stretch mode).
  !
  ! Samples ONE vibrational energy from P(E) ~ exp(-E/kT), then projects
  ! the amplitude onto the 3D bond direction. Both atoms are displaced
  ! with mass-weighted amplitudes.
  use venus_params, only: C1
  use venus_data, only: Q, P, W, NATOMS
  implicit none
  private

  public :: sample_surface_oscillator

  real(8), parameter :: kB_kcal = 0.00198717d0   ! kcal/mol/K
  real(8), parameter :: eV2kcal = 23.0605d0

contains

  ! ---------------------------------------------------------------------------
  ! Sample the surface harmonic oscillator from a classical Boltzmann
  ! distribution at temperature t_vib (Kelvin).
  !
  ! On entry: Q, P hold current coordinates and momenta.
  ! On exit:  Q and P for BOTH idx_atom and idx_anchor are modified.
  !
  ! k_surf(3): surface oscillator force constants (eV/Å²) for x, y, z.
  !            The effective k along the bond direction is:
  !              k_bond = k_surf(1)*dr(1)² + k_surf(2)*dr(2)² + k_surf(3)*dr(3)²
  ! ---------------------------------------------------------------------------
  subroutine sample_surface_oscillator(idx_atom, idx_anchor, k_surf, t_vib, iseed)
    integer, intent(in)    :: idx_atom     ! surface Au atom index (1-based)
    integer, intent(in)    :: idx_anchor   ! anchor/ghost atom index (1-based)
    real(8), intent(in)    :: k_surf(3)    ! force constants (eV/Å²)
    real(8), intent(in)    :: t_vib        ! vibrational temperature (K)
    integer, intent(inout) :: iseed        ! random seed

    real(8) :: r0(3), dr(3), k_bond, omega, mu, energy, amplitude, phase
    real(8) :: urand1, urand2
    integer :: kk
    real(8), external :: RAND0

    ! ---- bond direction from current equilibrium positions ----
    r0(1) = Q(3*idx_atom-2) - Q(3*idx_anchor-2)
    r0(2) = Q(3*idx_atom-1) - Q(3*idx_anchor-1)
    r0(3) = Q(3*idx_atom)   - Q(3*idx_anchor)
    dr(1:3) = r0(1:3) / dsqrt(r0(1)**2 + r0(2)**2 + r0(3)**2)

    ! ---- effective force constant along bond direction (eV/Å²) ----
    k_bond = k_surf(1)*dr(1)**2 + k_surf(2)*dr(2)**2 + k_surf(3)*dr(3)**2
    if (k_bond <= 0.0d0) k_bond = 1.0d0

    ! ---- reduced mass (amu) ----
    mu = W(idx_atom) * W(idx_anchor) / (W(idx_atom) + W(idx_anchor))

    ! ---- harmonic frequency in code units (1/(10 fs)) ----
    ! k_code = k_bond * eV2kcal * C1  (eV/Å² → code_energy/Å²)
    omega = dsqrt(k_bond * eV2kcal * C1 / mu)

    ! ---- classical Boltzmann energy sampling ----
    urand1 = RAND0(iseed)
    if (urand1 < 1.0d-15) urand1 = 1.0d-15
    energy = -kB_kcal * t_vib * dlog(urand1)          ! kcal/mol
    energy = energy * C1                               ! code units

    ! ---- amplitude: E = 0.5 * k_code * A² ----
    amplitude = dsqrt(2.0d0 * energy / (k_bond * eV2kcal * C1))  ! Å

    ! ---- random phase in [0, 2π) ----
    urand2 = RAND0(iseed)
    phase = 2.0d0 * 3.141592653589793d0 * urand2

    ! ---- mass-weighted displacement and velocity for BOTH atoms ----
    ! For a two-body harmonic oscillator:
    !   Δx₁ = +(μ/m₁) * A * cos(φ) * dr
    !   Δx₂ = -(μ/m₂) * A * cos(φ) * dr
    !   P₁ = +μ * v_rel   (v_rel = -ω*A*sin(φ))
    !   P₂ = -μ * v_rel
    do kk = 1, 3
       Q(3*idx_atom-3+kk)   = Q(3*idx_atom-3+kk)   &
            + (mu/W(idx_atom)) * amplitude * dcos(phase) * dr(kk)
       Q(3*idx_anchor-3+kk) = Q(3*idx_anchor-3+kk) &
            - (mu/W(idx_anchor)) * amplitude * dcos(phase) * dr(kk)

       P(3*idx_atom-3+kk)   = P(3*idx_atom-3+kk)   &
            + mu * (-omega * amplitude * dsin(phase)) * dr(kk)
       P(3*idx_anchor-3+kk) = P(3*idx_anchor-3+kk) &
            - mu * (-omega * amplitude * dsin(phase)) * dr(kk)
    end do
  end subroutine sample_surface_oscillator

end module harmonic_sampling
