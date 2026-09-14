!=====================================================================
! dist_boltzmann.f90 - distribution member (methods/ library): Boltzmann
!                      vibrational distribution - each normal mode's
!                      vibrational quantum number is drawn from the thermal
!                      geometric distribution P(n) = (1-q)*q^n with
!                      q = exp(-(h*c/k)*nu/T)
! Contract:
!   boltz_init(p) - validate the member parameters, precompute the per-mode
!     geometric constants and bind the served fragment. Loud-stops name the
!     offending field and value:
!       t_vib_a <= 0 (temperature [K]); w_mode_a unallocated (the zero
!       spectrum of a no-op frequency-derivation chain - legal loud-stop,
!       not a silent zero table); a non-positive entry in w_mode_a
!       [cm^-1]; n_mode_a /= size(w_mode_a).
!     Binding: the served fragment is the roster fragment whose
!     dist_scheme row selects THIS member (the member word looked up in the
!     sampler scheme-code table). Exactly one such fragment is required:
!     zero -> loud-stop; two or more -> loud-stop (one member instance
!     carries one parameter set, so a second carrier would be silently
!     skipped - the multi-fragment face is the retained t_vib_b field and
!     the un-finalized second-side frequency table).
!   boltz_sample(sta) - one draw per mode: n_i = floor(g/dum_i) with g a
!     unit-exponential deviate and dum_i = (h*c/k)*nu_i/T, which gives
!     P(n) = (1-q)*q^n exactly. The draw lands in the member's quantum-
!     number table (read back through boltz_last_n). sta is the sampler
!     runtime channel; this member writes NO q/p components - its runtime
!     product is the quantum-number table, because Cartesian normal-mode
!     assembly additionally needs the mode eigenvectors of the
!     frequency-derivation chain (mass-weighted Hessian eigen-problem at
!     the fragment equilibrium geometry, an assembly-time mechanism that
!     hands this member its frequency table; through a no-op chain the
!     spectrum is absent and boltz_init loud-stops - the degradation
!     posture is fail-loud, never a silent zero table).
!   boltz_last_n(m) - copy of the last drawn quantum-number table entry.
! Frequency derivation (assembly-time mechanism, not member code): the
!   w_mode_a table arrives as member data from the Hessian eigenvalue chain
!   (force procedure -> mass-weighted Hessian -> eigenvalues -> sqrt scaling
!   -> [cm^-1]); the member only stores and consumes it.
! Units: nu [cm^-1], T [K]; h*c/k = 1.43878 cm*K (member constant hc_k);
!   quantum numbers are dimensionless counts.
!=====================================================================
module dist_boltzmann
   use rng,           only: rng_gamma
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use sampler,       only: samp_code
   implicit none
   private
   public :: boltz_params_t, boltz_init, boltz_sample, boltz_last_n

   real(8), parameter :: hc_k = 1.43878d0   ! h*c/k [cm*K] - the geometric step constant

   ! Member-owned parameters (units in field comments; the second-side
   ! temperature t_vib_b is retained for the multi-fragment face and is
   ! NOT consumed while no second-side frequency table is finalized)
   type :: boltz_params_t
      real(8) :: t_vib_a = 0.0d0          ! fragment-A vibrational temperature [K]
      real(8) :: t_vib_b = 0.0d0          ! fragment-B vibrational temperature [K] (retained, unused)
      real(8), allocatable :: w_mode_a(:) ! fragment-A normal-mode frequencies [cm^-1]
      integer :: n_mode_a = 0             ! fragment-A sampled-mode count [count]
   end type

   ! Derived member state (module-private; rebuilt by every boltz_init)
   real(8), allocatable :: dum_mode(:)  ! per-mode geometric step (h*c/k)*nu/T [-]
   integer, allocatable :: last_n(:)    ! last drawn quantum-number table [-]
   integer :: frag_served = 0           ! bound roster fragment index [index]
   logical :: inited = .false.          ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! boltz_init(p) - validate parameters, precompute geometric constants,
   !                 bind the served fragment
   subroutine boltz_init(p)
      type(boltz_params_t), intent(in) :: p ! member parameters (vibrational temperature /
                                            ! frequency table)
      integer :: i, code, n_carry, i_carry
      real(8) :: val

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%t_vib_a <= 0.0d0) then
         val = p%t_vib_a
         call stop_boltz('boltz_init', 'vibrational temperature t_vib_a must be positive, got', val)
      end if
      if (.not. allocated(p%w_mode_a)) then
         call stop_boltx('boltz_init', 'no frequency spectrum (w_mode_a unallocated) - the '// &
            'frequency-derivation chain produced no table; a zero spectrum is not a sampleable state')
      end if
      if (size(p%w_mode_a) /= p%n_mode_a) then
         call stop_bolti2('boltz_init', 'mode count mismatch: n_mode_a =', p%n_mode_a, &
            'but size(w_mode_a) =', size(p%w_mode_a))
      end if
      do i = 1, size(p%w_mode_a)
         if (p%w_mode_a(i) <= 0.0d0) then
            val = p%w_mode_a(i)
            call stop_bolti('boltz_init', 'mode', i, 'carries a non-positive frequency [cm^-1]:', val)
         end if
      end do

      ! 2. fragment binding (the roster fragment whose dist_scheme row selects
      !    this member; assembly-order and count guards first)
      if (.not. allocated(roster%frag)) then
         call stop_boltx('boltz_init', 'no roster fragment table - the roster must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_boltx('boltz_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_bolti2('boltz_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('boltzmann')
      if (code == 0) then
         call stop_boltx('boltz_init', 'member word "boltzmann" is not in the scheme-code table '// &
            '(incomplete member library assembly)')
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
         call stop_boltx('boltz_init', 'no roster fragment carries the boltzmann scheme - '// &
            'this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_bolti2('boltz_init', 'two or more fragments carry the boltzmann scheme, but '// &
            'one member instance carries one parameter set (second-side face not wired); count =', &
            n_carry, 'first carrier =', i_carry)
      end if
      frag_served = i_carry

      ! 3. per-mode geometric constants: dum_i = (h*c/k)*nu_i/T, q_i = exp(-dum_i)
      if (allocated(dum_mode)) deallocate (dum_mode)
      if (allocated(last_n)) deallocate (last_n)
      allocate (dum_mode(p%n_mode_a), last_n(p%n_mode_a))
      do i = 1, p%n_mode_a
         dum_mode(i) = hc_k*p%w_mode_a(i)/p%t_vib_a
      end do
      last_n = 0
      inited = .true.
   end subroutine boltz_init

   !------------------------------------------------------------------
   ! boltz_sample(sta) - draw one quantum-number table (all modes of the
   !                     bound fragment); sta passes through untouched
   subroutine boltz_sample(sta)
      type(state_t), intent(inout) :: sta  ! sampler runtime channel (no q/p components owned)
      integer :: m
      if (.not. inited) then
         call stop_boltx('boltz_sample', 'sample called before boltz_init - the member carries '// &
            'no parameter set')
      end if
      do m = 1, size(dum_mode)
         last_n(m) = int(rng_gamma(1)/dum_mode(m))   ! floor(g/dum): P(n)=(1-q)*q^n exactly
      end do
   end subroutine boltz_sample

   !------------------------------------------------------------------
   ! boltz_last_n(m) - copy of the last drawn quantum number of mode m
   pure integer function boltz_last_n(m)
      integer, intent(in) :: m
      boltz_last_n = 0
      if (inited .and. m >= 1 .and. m <= size(last_n)) boltz_last_n = last_n(m)
   end function boltz_last_n

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_boltx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_boltx

   subroutine stop_boltz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_boltz

   subroutine stop_bolti(where, msg, ival, tail, val)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival
      real(8), intent(in) :: val
      write (0, '(a,a,i0,a,es12.4)') trim(where)//': ', trim(msg), ival, ' '//trim(tail), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_bolti

   subroutine stop_bolti2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_bolti2

end module dist_boltzmann
