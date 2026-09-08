!=====================================================================
! dist_ebk.f90 - L3 example distribution member (B' registry family): EBK semiclassical state
!                 distribution - fixed
!                 (n,J) quantum numbers, a diatomic initial-state member (the B' second
!                 implementation prototype; example 2/3)
! Reference: src_VENUS/INITEBK.f90:1-38 (energy-level initial value and frequency step :15-16,
!       rotational-energy moment
!       AL = sqrt(J(J+1))*C7 :17-18, reduced mass :19, fixed-point iteration :24-33,
!       turning-point
!       inset :35-36, momentum scale PTEST :37), src_VENUS/FINLNJ.f:1-104 (the action
!       integral's RKR inversion - effective potential VEFF = V-DH+centrifugal term :48,
!       turning-point scan :56-81,
!       Gauss-Legendre quadrature :85-95), src_VENUS/SELECT.f90:337-351 (the fixed-quantum-
!       number branch
!       TRV_A<0 and the energy-level initial value ENJA=(n+1/2)*nu*cm2cal*C1 :347),
!       src_VENUS/HOMOQP.f90 (r/pr -> diatomic coordinate/momentum assembly - called at
!       SELECT.f90:373),
!       cases/ebk_fixed_nj/README.md:18-25 (energy-level formula hc[w_e(n+1/2) - w_exe(n+1/2)^2]
!       + hc*B_e*J(J+1) against the Morse spectroscopic constants),
!       cases/ebk_fixed_nj/input_qct.txt:36-39
!       (the fixed-quantum-number / moment-of-inertia input slots)
! Registration: this member enters the sampler distribution member table (one samp_reg line -
! executed in the assembly phase; C5 wiring)
!=====================================================================
module dist_ebk
   implicit none
   private
   public :: ebk_params_t, ebk_init

   ! Member-owned parameters (E2; listed per what the old INITEBK and its call chain actually
   ! read, blanks allowed;
   ! —— final shape in phase C. Reduced mass is not stored — it is computed
   ! from the roster mass table (old W(LBOND),
   ! INITEBK.f90:19); the dissociation energy is not stored - the old SELECT.f90:121 DH=V is
   ! evaluated by the facade,
   ! energy belongs to the facade)
   type :: ebk_params_t
      integer :: n_vib  = 0        ! fixed vibrational quantum number n [-] (old NNA - the
                                   ! TRV_A<0 prescription, input_qct.txt:37)
      integer :: j_rot  = 0        ! fixed rotational quantum number J [-] (old JA - same, :38)
      real(8) :: w_e    = 0.0d0    ! spectroscopic constant w_e [cm^-1] (Morse harmonic term -
                                   ! the linear (n+1/2) term; cross-check value 2326.4)
      real(8) :: w_ex_e = 0.0d0    ! anharmonic constant w_exe [cm^-1] (the (n+1/2)^2 term;
                                   ! cross-check value 35.35 - derived from the potential
      ! parameters, finalized in C5)
      real(8) :: b_rot  = 0.0d0    ! rotational constant B_e [cm^-1] (rigid-rotor term J(J+1))
      real(8) :: ai_rot = 0.0d0    ! moment of inertia [amu·Å²] (old AIA - for assembly/orientation;
      ! cross-check value 0.98925, input_qct.txt:39)
      ! Left blank: turning-point initial values r_min/r_max [Å], fixed-point iteration
      ! tolerance (old 1e-6) and cap
      ! (old 1000 iterations, INITEBK.f90:31-32), quadrature node count (old 50) - finalized
      ! in C5
   end type
contains
   !------------------------------------------------------------------
   ! ebk_init(p) - initialize the member (validate quantum numbers / spectroscopic constants,
   ! precompute the energy level and rotational moment)
   ! Reference: src_VENUS/INITEBK.f90:15-19 (HNU step, AL = sqrt(J(J+1))*C7, reduced mass)
   subroutine ebk_init(p)
      type(ebk_params_t), intent(in) :: p ! member parameters (quantum numbers / spectroscopic
      ! constants)
      ! Pseudocode:
      ! 1. Validate parameters (illegal -> loud-stop): n_vib >= 0, j_rot >= 0, w_e > 0,
      !    w_ex_e >= 0, b_rot >= 0, ai_rot > 0
      ! 2. Precompute the state energy-level initial value [kcal/mol]:
      !    e_nj = (w_e*(n+1/2) - w_ex_e*(n+1/2)^2 + b_rot*J(J+1))*wvn2e
      !    (cm^-1 -> kcal/mol via consts; the old initial value (n+1/2)*nu*cm2cal*C1,
      !    SELECT.f90:347 - energy-level formula cross-check cases/ebk_fixed_nj/README.md:20)
      ! 3. Precompute the rotational moment al = sqrt(J(J+1))*hbar_code [hbar units] (old
      !    AL = AM*C7,
      !    INITEBK.f90:17-18) and the frequency step hnu = e_nj/(n+1/2) (old HNU :16)
   end subroutine ebk_init
   ! -- sample direction (finalized in C5) -- fixed-point iteration e_nj <- e_nj + (n-an)*hnu
   !   until |n-an| <= 1e-6 (an is a FINLNJ-type action integral - WKB density p(r) ~ 1/|p_r|,
   !   rejection-sample the internuclear distance r within [r_min,r_max], radial momentum
   !   pr = sqrt(2*mu*summ) with the sign
   !   flipped with probability 1/2, then a HOMOQP-type assembly writes the initial-state
   !   components)
   ! sample routines: finalized when the C5 members are concretized (not declared at this stage)
   ! Registration-line note: samp_reg registers this member (name 'ebk') in the sampler
   ! distribution member table - C5 wiring
end module dist_ebk
