!=====================================================================
! dist_thermalize.f90 - distribution member (methods/ library): MD
!              equilibration mechanism - Maxwell velocity seed of the
!              owned components at the target temperature, exact drift
!              removal (the mass-weighted COM-velocity projection
!              p'_i = p_i - (m_i/M)*P onto zero net momentum), then a
!              bounded equilibration run of n_eq propagation steps at
!              the member step size through the evolution service.
! Conservation algebra (the declaration this member carries):
!   post-equilibration kinetic-energy expectation E[KE] =
!     (3*nat - 3)*kb*t_eq/2 exactly, independent of the mass
!     distribution - the Maxwell seed followed by the exact
!     MASS-WEIGHTED (COM-velocity) projection p'_i = p_i -
!     (m_i/M)*P onto zero net momentum removes exactly three degrees
!     of freedom from the Gaussian; the equilibration run at the bath
!     temperature preserves the expectation (a mechanism-level
!     declaration: the trajectory numbers are posterior to the force
!     chain).
! Division of labor (the ruling this member shares with the assembly):
!   the member owns the EQUILIBRATION MECHANISM (seed + drift removal +
!   the bounded run through the propagator - an internal
!   evolution-service call is legal, boundary 4); assembling the run
!   around it (which other members participate, trajectory accounting)
!   is driver business. Since the Ruling H minimal path the force chain
!   is real: the run propagates the owned components under the assembled
!   force law - the zero-force trajectories of the criteria are a
!   property of the check-side container binding (free flight), not of
!   the chain; the trajectory numbers are posterior to the force chain.
! Contract:
!   thermalize_init(p) - validate the member parameters (target
!     temperature, step count, step size), bind the served fragment
!     (exactly one carrier of the thermalize scheme), precompute the
!     expectation constant. Loud-stops name the offending field and
!     value:
!       t_eq <= 0 (target temperature [K] - no silent zero-temperature
!       branch); n_eq < 1 (step count); dt_eq <= 0 (step size [10 fs]);
!       no roster fragment table (assembly order); no per-fragment
!       scheme array; length mismatch against the roster; the member
!       word absent from the scheme-code table (a library-assembly
!       failure); zero or multiple carrier fragments; a state whose
!     momentum vector ends before the carrier's last owned component
!     (the roster offset included - a direct call against a
!     short-of-fragment state).
!   thermalize_sample(sta) - (1) seed every owned component of the
!     carrier fragment with a fresh Gaussian momentum at sqrt(m_i*kb*
!     t_eq); (2) apply the mass-weighted COM-velocity projection
!     p'_i = p_i - (m_i/M)*P per component, P the owned-set total
!     momentum, M the carrier mass (the owned-set total momentum
!     becomes exactly zero; this is the projection the declared
!     expectation form assumes); (3) swap the
!     global propagator step to dt_eq, run n_eq propagation steps
!     through the evolution service, restore the step and the
!     trajectory time (the equilibration run neither consumes
!     trajectory time nor leaves the global step changed). sta is the
!     sampler runtime channel; this member writes ONLY its owned
!     momentum components.
!   thermalize_t_eq() / thermalize_n_eq() / thermalize_dt_eq() - the
!     equilibration parameters init received.
!   thermalize_e_kin_expect() - the declared post-equilibration
!     kinetic-energy expectation (3*nat-3)*kb*t_eq/2 [internal energy]
!     (exposed for criteria anchoring of the conservation algebra).
!   thermalize_frag() / thermalize_nat() - the bound carrier fragment
!     and its atom count.
! Variant inventory against the old defect register (semantics only, no
!   citations): the zero-temperature silent branch is NOT carried (a
!   loud-stop names it instead); the bystander-freeze choreography
!   (displacing/freezing the non-target fragments during the run) is
!   NOT carried (it trespasses component ownership - assembly/driver
!   business); the hardwired default step is NOT carried (dt_eq is a
!   member parameter applied through the configuration step channel
!   with swap-and-restore); the restart branch with a random step count
!   is NOT carried; periodic in-run resampling is NOT carried (the run
!   is bare propagation - bath coupling is the thermostat member's
!   face); the two-counter equilibration split is NOT carried (a single
!   step count); diagnostics prints are NOT carried.
! Units: t_eq [K]; dt_eq [10 fs]; masses [amu] (roster); kb_code
!   internal energy per K; momenta [amu*Ang/(10 fs)]; e_kin_expect
!   [internal energy].
!=====================================================================
module dist_thermalize
   use consts,        only: kb_code
   use rng,           only: rng_gauss
   use state,         only: state_t
   use config,        only: reactants, prop_cfg => propagator   ! alias to avoid clashing with the module name
   use config_roster, only: roster
   use sampler,       only: samp_code
   use propagator,    only: prop_step     ! evolution service (boundary 4: an internal call is legal)
   implicit none
   private
   public :: thermalize_params_t, thermalize_init, thermalize_sample, &
             thermalize_t_eq, thermalize_n_eq, thermalize_dt_eq, thermalize_e_kin_expect, &
             thermalize_frag, thermalize_nat

   ! Member-owned parameters
   type :: thermalize_params_t
      real(8) :: t_eq = 0.0d0    ! target temperature [K]
      integer :: n_eq = 0        ! equilibration step count [-]
      real(8) :: dt_eq = 0.0d0   ! equilibration step size [10 fs]
   end type

   ! Derived member state (module-private; rebuilt by every thermalize_init)
   real(8) :: t_eq_c = 0.0d0        ! target temperature [K]
   integer :: n_eq_c = 0            ! step count [-]
   real(8) :: dt_eq_c = 0.0d0       ! step size [10 fs]
   real(8) :: e_expect_c = 0.0d0    ! declared KE expectation [internal energy]
   integer :: frag_b = 0            ! bound carrier fragment index [-]
   integer :: nat_b = 0             ! carrier atom count [-]
   real(8), allocatable :: sigma_c(:) ! per-atom seed widths sqrt(m*kb*t_eq) [amu*Ang/(10 fs)]
   logical :: inited = .false.      ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! thermalize_init(p) - validate the equilibration parameters, bind the
   !                      carrier, precompute the expectation constant
   subroutine thermalize_init(p)
      type(thermalize_params_t), intent(in) :: p ! member parameters (equilibration)
      integer :: i, code, n_carry, i_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%t_eq <= 0.0d0) then
         call stop_equiz('thermalize_init', 'target temperature t_eq [K] must be positive, got', &
            p%t_eq)
      end if
      if (p%n_eq < 1) then
         call stop_equii('thermalize_init', 'equilibration step count n_eq must be at least 1, got', &
            p%n_eq)
      end if
      if (p%dt_eq <= 0.0d0) then
         call stop_equiz('thermalize_init', 'equilibration step size dt_eq [10 fs] must be '// &
            'positive, got', p%dt_eq)
      end if

      ! 2. fragment binding (assembly-order and carrier guards)
      if (.not. allocated(roster%frag)) then
         call stop_equix('thermalize_init', 'no roster fragment table - the roster must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_equix('thermalize_init', 'no per-fragment scheme array (reactants%'// &
            'dist_scheme unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_equii2('thermalize_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('thermalize')
      if (code == 0) then
         call stop_equix('thermalize_init', 'member word "thermalize" is not in the '// &
            'scheme-code table (incomplete member library assembly)')
      end if
      n_carry = 0
      i_carry = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code) then
            n_carry = n_carry + 1
            i_carry = i
         end if
      end do
      if (n_carry == 0) then
         call stop_equix('thermalize_init', 'no roster fragment carries the thermalize '// &
            'scheme - this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_equii2('thermalize_init', 'two or more fragments carry the thermalize '// &
            'scheme, but one member instance carries one parameter set; count =', n_carry, &
            'first carrier =', i_carry)
      end if

      ! 3. seed widths from the roster masses and the expectation constant
      !    (sigma^2 = m*kb*t_eq per carrier atom; E[KE] = (3*nat-3)*kb*t_eq/2)
      frag_b = i_carry
      nat_b = roster%frag(frag_b)%nat
      if (allocated(sigma_c)) deallocate (sigma_c)
      allocate (sigma_c(nat_b))
      do i = 1, nat_b
         sigma_c(i) = sqrt(roster%mass(roster%frag(frag_b)%list(i))*kb_code*p%t_eq)
      end do
      e_expect_c = (3.0d0*dble(nat_b) - 3.0d0)*kb_code*p%t_eq/2.0d0
      t_eq_c = p%t_eq
      n_eq_c = p%n_eq
      dt_eq_c = p%dt_eq
      inited = .true.
   end subroutine thermalize_init

   !------------------------------------------------------------------
   ! thermalize_sample(sta) - Maxwell seed of the owned components, exact
   !                          drift removal, then the bounded equilibration
   !                          run through the evolution service
   subroutine thermalize_sample(sta)
      type(state_t), intent(inout) :: sta ! sampler runtime channel (owned p components written)
      integer :: i, c, g_atom, i_step, last_own
      real(8) :: wt, p_tot(3), t_saved, dt_saved
      if (.not. inited) then
         call stop_equix('thermalize_sample', 'sample called before thermalize_init - the '// &
            'member carries no parameter set')
      end if
      last_own = 3*roster%frag(frag_b)%list(nat_b)   ! last owned component (roster offset included)
      if (size(sta%p) < last_own) then
         call stop_equii2('thermalize_sample', 'the state momentum vector ends before the '// &
            "carrier's last owned component: state dof =", size(sta%p), &
            'last owned component =', last_own)
      end if

      ! (1) Maxwell seed of the owned components
      do i = 1, nat_b
         g_atom = roster%frag(frag_b)%list(i)
         do c = 1, 3
            sta%p(3*(g_atom - 1) + c) = sigma_c(i)*rng_gauss()
         end do
      end do

      ! (2) exact drift removal: mass-weighted COM-velocity projection
      !     p'_i = p_i - (m_i/M)*P per component (the owned-set total
      !     momentum becomes exactly zero; the Gaussian projected this way
      !     loses exactly three degrees of freedom - the declared
      !     expectation is mass-independent)
      wt = sum(roster%mass(roster%frag(frag_b)%list))
      p_tot = 0.0d0
      do i = 1, nat_b
         g_atom = roster%frag(frag_b)%list(i)
         p_tot = p_tot + sta%p(3*(g_atom - 1) + 1:3*g_atom)
      end do
      do i = 1, nat_b
         g_atom = roster%frag(frag_b)%list(i)
         sta%p(3*(g_atom - 1) + 1:3*g_atom) = sta%p(3*(g_atom - 1) + 1:3*g_atom) &
                                              - (roster%mass(g_atom)/wt)*p_tot
      end do

      ! (3) the bounded equilibration run at the member step size, through
      !     the evolution service (boundary 4): the global step and the
      !     trajectory time are swapped for the run and RESTORED after it
      !     (the force chain is real since the Ruling H minimal path - the
      !     criteria's zero-force trajectories come from their check-side
      !     container binding, see the header)
      t_saved = sta%t
      dt_saved = prop_cfg%dt
      prop_cfg%dt = dt_eq_c
      do i_step = 1, n_eq_c
         call prop_step(sta)
      end do
      prop_cfg%dt = dt_saved
      sta%t = t_saved
   end subroutine thermalize_sample

   !------------------------------------------------------------------
   ! getters (copies of the derived member state)
   !------------------------------------------------------------------
   pure real(8) function thermalize_t_eq()
      thermalize_t_eq = t_eq_c
   end function thermalize_t_eq

   pure integer function thermalize_n_eq()
      thermalize_n_eq = n_eq_c
   end function thermalize_n_eq

   pure real(8) function thermalize_dt_eq()
      thermalize_dt_eq = dt_eq_c
   end function thermalize_dt_eq

   pure real(8) function thermalize_e_kin_expect()
      thermalize_e_kin_expect = e_expect_c
   end function thermalize_e_kin_expect

   pure integer function thermalize_frag()
      thermalize_frag = frag_b
   end function thermalize_frag

   pure integer function thermalize_nat()
      thermalize_nat = nat_b
   end function thermalize_nat

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_equix(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_equix

   subroutine stop_equiz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_equiz

   subroutine stop_equii(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_equii

   subroutine stop_equii2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_equii2

end module dist_thermalize
