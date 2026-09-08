!=====================================================================
! dist_boltzmann.f90 - L3 example distribution member (B' registry family): Boltzmann
!                       vibration distribution -
!                       a member prototype that samples each normal-mode quantum number via a
!                       geometric distribution set by the vibrational temperature (example 1/3)
! Reference: src_VENUS/THRMAN.f:1-19 (thermal population -> vibrational quantum numbers:
!       constant step DUM = h*nu/k_BT
!       = 1.43878*nu/T at :12-14, geometric-distribution truncation :15-16),
!       src_VENUS/GAMA.f:1-14
!       (integer-order Gamma deviates - first order is the exponential deviate -ln U, :7-12),
!       src_VENUS/SELECT.f90:308-336 (temperature/frequency sources: the TRV_A temperature and
!       the WWA harmonic frequencies
!       = EIG(6)*C6 - Hessian eigenvalues obtained via EIGN at the equilibrium geometry, :321),
!       :649 (the polyatomic
!       normal-mode branch TVIBA), :863/:988 (the two B-side TRVB/TVIBB call sites for
!       comparison),
!       cases/boltzmann_vib/README.md:16-20 (P(n)=(1-q)*q^n, q=exp(-1.43878*nu/T)
!       - theoretical cross-check); rotational-J temperature sampling (PROBJ) belongs to the
!       rotation member - the remaining 9 -> C5
! Registration: this member enters the sampler distribution member table (one samp_reg line -
! executed in the assembly phase; C5 wiring)
!=====================================================================
module dist_boltzmann
   implicit none
   private
   public :: boltz_params_t, boltz_init

   ! Member-owned parameters (E2: member = code + its own data; fields listed per what the
   ! old THRMAN actually reads,
   ! input parameters enumerated, units noted in comments, blanks allowed - finalized in
   ! phase C)
   type :: boltz_params_t
      real(8) :: t_vib_a = 0.0d0          ! fragment-A vibrational temperature [K] (old TRVA/TVIBA -
                                          ! a negative value prescribes fixed quantum numbers,
                                          ! routing to the ebk member, SELECT.f90:308/:337-339 fork)
      real(8) :: t_vib_b = 0.0d0          ! fragment-B vibrational temperature [K] (old TRVB/TVIBB -
                                          ! same form on the B side, SELECT.f90:863/:988)
      real(8), allocatable :: w_mode_a(:) ! fragment-A normal-mode frequencies [cm^-1] (old WWA -
                                          ! assembled from Hessian eigenvalues, SELECT.f90:321;
                                          ! spectrum computation belongs to the assembly phase,
      ! the member only stores the table)
      integer :: n_mode_a = 0             ! fragment-A sampled-mode count [count] (old NMBAR -
      ! 1 for diatomics, 3N-6/3N-5 for polyatomics)
      ! Left blank: B-side frequency tables (w_mode_b/n_mode_b), per-mode geometric constants
      ! and energy-cap caches,
      ! multi-fragment extension (K1 risk note max_frag) - finalized when the C5 members are
      ! concretized
   end type
contains
   !------------------------------------------------------------------
   ! boltz_init(p) - initialize the member (validate the temperature/frequency tables,
   ! precompute the geometric-distribution constants)
   ! Reference: src_VENUS/THRMAN.f:12-17 (constant step and geometric truncation - see the
   ! module header)
   subroutine boltz_init(p)
      type(boltz_params_t), intent(in) :: p ! member parameters (vibrational temperatures /
      ! frequency tables)
      ! Pseudocode:
      ! 1. Validate parameters (illegal -> loud-stop): t_vib_a/t_vib_b > 0, w_mode_a
      !    allocated with all frequencies positive, n_mode_a matching the table length
      ! 2. Precompute per-mode geometric-distribution constants: dum_i = h*nu_i/k_BT =
      !    1.43878*nu_i/T
      !    (the old :12-14 C7/C5/T constant combination - the new implementation regroups the
      !    consts constants, exact
      !    conversion chain finalized in C5) and q_i = exp(-dum_i), mean <n_i> = q_i/(1-q_i)
   end subroutine boltz_init
   ! -- sample direction (finalized in C5) -- per mode n_i = trunc(GAMA(1)/dum_i) (first-order
   !   Gamma, i.e. the rng exponential deviate rng_gamma; hence P(n)=(1-q)*q^n - cross-check
   !   with cases), write the
   !   fragment's vibrational initial-state components
   ! sample routines: finalized when the C5 members are concretized (not declared at this stage)
   ! Registration-line note: samp_reg registers this member (name 'boltzmann') in the sampler
   ! distribution member table
   ! - C5 wiring
end module dist_boltzmann
