!=====================================================================
! dist_thermostat.f90 - distribution member (methods/ library): bath
!              thermostat - the velocity-reset mechanism of an Andersen
!              bath: one reset draws FRESH Gaussian momenta for the owned
!              components at the bath temperature,
!              p(i,c) = sigma(i)*g, sigma(i)^2 = m_i*kb*T_bath
!              (the Maxwell-shape declaration this member carries).
! Division of labor (the ruling this member shares with the dynamics
!   side): the member owns the RESET MECHANISM and the bath parameters;
!   the collision schedule (when a reset fires during a run) and the
!   Langevin integration are force-chain/propagator business. The
!   collision frequency is loaded bath data for that side, exposed
!   through thermostat_nu - loaded, not dropped.
! Contract:
!   thermostat_init(p) - validate the member parameters (bath
!     temperature, collision frequency), bind the served fragment
!     (exactly one carrier of the thermostat scheme), precompute the
!     per-atom Maxwell widths from the roster masses. Loud-stops name
!     the offending field and value:
!       t_bath <= 0 (bath temperature [K]); nu_coll <= 0 (collision
!       frequency [1/(10 fs)]); no roster fragment table (assembly
!       order); no per-fragment scheme array; length mismatch against
!       the roster; the member word absent from the scheme-code table
!       (a library-assembly failure); zero or multiple carrier
!       fragments (one member instance carries one parameter set).
!   thermostat_sample guards: sample before init; a state whose
!     momentum vector ends before the carrier's last owned component
!     (the roster offset included - a direct call against a
!     short-of-fragment state). One reset: every owned component of the
!     carrier fragment receives a fresh Gaussian momentum at the
!     declared width. sta is the sampler runtime channel; this member
!     writes ONLY its owned momentum components (q/f untouched, and the
!     atoms outside the carrier untouched).
!   thermostat_t_bath() / thermostat_nu() - the bath parameters init
!     received (collision frequency = carried bath data).
!   thermostat_frag() / thermostat_nat() - the bound carrier fragment
!     and its atom count.
!   thermostat_sigma(i) - the declared Maxwell width sqrt(m_i*kb*T) of
!     carrier-atom i [amu*Ang/(10 fs)] (exposed for criteria anchoring
!     of the conservation algebra).
! Variant inventory against the old defect register (semantics only, no
!   citations): the brute rescale-to-exact-temperature face is NOT
!   carried (a deterministic rescale biases the declared Maxwell
!   marginals - the reset is the exact Gaussian draw); the periodic
!   in-run resampling schedule is NOT carried (dynamics-side business,
!   the member is the reset mechanism); net-momentum zeroing is NOT
!   carried here (drift removal is the equilibration face's business -
!   doing it inside a reset would bias the declared marginals); a
!   rigid-atom exclusion list is NOT carried (the roster is the sole
!   composition channel, the carrier is whole); diagnostics prints are
!   NOT carried.
! Units: t_bath [K]; nu_coll [1/(10 fs)]; masses [amu] (roster); kb_code
!   internal energy per K; sigma [amu*Ang/(10 fs)] = internal momentum.
!=====================================================================
module dist_thermostat
   use consts,        only: kb_code
   use rng,           only: rng_gauss
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use sampler,       only: samp_code
   implicit none
   private
   public :: thermostat_params_t, thermostat_init, thermostat_sample, &
             thermostat_t_bath, thermostat_nu, thermostat_frag, thermostat_nat, &
             thermostat_sigma

   ! Member-owned parameters
   type :: thermostat_params_t
      real(8) :: t_bath = 0.0d0   ! bath temperature [K]
      real(8) :: nu_coll = 0.0d0  ! collision frequency [1/(10 fs)] (bath data for the dynamics side)
   end type

   ! Derived member state (module-private; rebuilt by every thermostat_init)
   real(8) :: t_bath_c = 0.0d0       ! bath temperature [K]
   real(8) :: nu_coll_c = 0.0d0      ! collision frequency [1/(10 fs)]
   integer :: frag_b = 0             ! bound carrier fragment index [-]
   integer :: nat_b = 0              ! carrier atom count [-]
   real(8), allocatable :: sigma_c(:) ! per-atom Maxwell widths sqrt(m*kb*T) [amu*Ang/(10 fs)]
   logical :: inited = .false.       ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! thermostat_init(p) - validate the bath parameters, bind the carrier
   !                      fragment, precompute the Maxwell widths
   subroutine thermostat_init(p)
      type(thermostat_params_t), intent(in) :: p ! member parameters (bath)
      integer :: i, code, n_carry, i_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%t_bath <= 0.0d0) then
         call stop_thermoz('thermostat_init', 'bath temperature t_bath [K] must be positive, got', &
            p%t_bath)
      end if
      if (p%nu_coll <= 0.0d0) then
         call stop_thermoz('thermostat_init', 'collision frequency nu_coll [1/(10 fs)] must be '// &
            'positive, got', p%nu_coll)
      end if

      ! 2. fragment binding (assembly-order and carrier guards)
      if (.not. allocated(roster%frag)) then
         call stop_thermox('thermostat_init', 'no roster fragment table - the roster must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_thermox('thermostat_init', 'no per-fragment scheme array (reactants%'// &
            'dist_scheme unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_thermoi2('thermostat_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('thermostat')
      if (code == 0) then
         call stop_thermox('thermostat_init', 'member word "thermostat" is not in the '// &
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
         call stop_thermox('thermostat_init', 'no roster fragment carries the thermostat '// &
            'scheme - this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_thermoi2('thermostat_init', 'two or more fragments carry the thermostat '// &
            'scheme, but one member instance carries one parameter set; count =', n_carry, &
            'first carrier =', i_carry)
      end if

      ! 3. Maxwell widths from the roster masses (the conservation algebra:
      !    sigma^2 = m*kb*T_bath, per carrier atom)
      frag_b = i_carry
      nat_b = roster%frag(frag_b)%nat
      if (allocated(sigma_c)) deallocate (sigma_c)
      allocate (sigma_c(nat_b))
      do i = 1, nat_b
         sigma_c(i) = sqrt(roster%mass(roster%frag(frag_b)%list(i))*kb_code*p%t_bath)
      end do
      t_bath_c = p%t_bath
      nu_coll_c = p%nu_coll
      inited = .true.
   end subroutine thermostat_init

   !------------------------------------------------------------------
   ! thermostat_sample(sta) - one velocity reset: fresh Gaussian momenta
   !                          for every owned component of the carrier
   subroutine thermostat_sample(sta)
      type(state_t), intent(inout) :: sta ! sampler runtime channel (owned p components written)
      integer :: i, c, g_atom, last_own
      if (.not. inited) then
         call stop_thermox('thermostat_sample', 'sample called before thermostat_init - the '// &
            'member carries no parameter set')
      end if
      last_own = 3*roster%frag(frag_b)%list(nat_b)   ! last owned component (roster offset included)
      if (size(sta%p) < last_own) then
         call stop_thermoi2('thermostat_sample', 'the state momentum vector ends before the '// &
            "carrier's last owned component: state dof =", size(sta%p), &
            'last owned component =', last_own)
      end if
      do i = 1, nat_b
         g_atom = roster%frag(frag_b)%list(i)
         do c = 1, 3
            sta%p(3*(g_atom - 1) + c) = sigma_c(i)*rng_gauss()
         end do
      end do
   end subroutine thermostat_sample

   !------------------------------------------------------------------
   ! getters (copies of the derived member state)
   !------------------------------------------------------------------
   pure real(8) function thermostat_t_bath()
      thermostat_t_bath = t_bath_c
   end function thermostat_t_bath

   pure real(8) function thermostat_nu()
      thermostat_nu = nu_coll_c
   end function thermostat_nu

   pure integer function thermostat_frag()
      thermostat_frag = frag_b
   end function thermostat_frag

   pure integer function thermostat_nat()
      thermostat_nat = nat_b
   end function thermostat_nat

   pure real(8) function thermostat_sigma(i)
      integer, intent(in) :: i ! carrier-atom index [-]
      thermostat_sigma = sigma_c(i)
   end function thermostat_sigma

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_thermox(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_thermox

   subroutine stop_thermoz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_thermoz

   subroutine stop_thermoi2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_thermoi2

end module dist_thermostat
