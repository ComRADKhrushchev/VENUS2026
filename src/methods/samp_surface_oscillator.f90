!=====================================================================
! samp_surface_oscillator.f90 - distribution member: analytic harmonic
!   stretch sampling of a two-atom atom-anchor carrier - one energy
!   draw (thermal or quantum level) on the E-shell with a uniform
!   phase, COM fixed, zero total momentum
! Design:
!   delta = A*cos(phi) along the equilibrium bond, v = -omega*A*sin(phi);
!   positions partition by mass fractions so the COM does not move.
!   One scalar force constant; the member overwrites (never adds to)
!   its carrier's q/p components.
!   Units: k_osc [kcal/mol/Ang^2] via e_conv; E [internal]; omega [rad/(10 fs)]; q [Ang]; p [amu*Ang/(10 fs)].
!=====================================================================
module samp_surface_oscillator
   use consts,        only: hbar_code, kb_code, e_conv, two_pi
   use rng,           only: rng_u
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code
   implicit none
   private
   public :: surf_osc_params_t, surface_oscillator_init, surface_oscillator_draw, surface_oscillator_realize, &
             surface_oscillator_sample, surface_oscillator_last_e

   ! Member-owned parameters
   type :: surf_osc_params_t
      integer :: n_osc = 0         ! energy-law selector: 0 thermal, 1 quantum level [-]
      real(8) :: t_osc = 300.0d0   ! oscillator temperature [K] (mode 0)
      integer :: n_level = 0       ! vibrational quantum number n [-] (mode 1)
      real(8) :: k_osc = 30.0d0    ! harmonic force constant [kcal/mol/Ang^2]
   end type

   ! Derived member state (module-private; rebuilt by every surface_oscillator_init) -
   ! one row per carrier; the law constants are shared (equal treatment)
   type :: carry_t
      integer :: frag = 0                   ! carrier list_atoms fragment row [index]
      integer, allocatable :: list(:)       ! carrier atom index list [global atom index]
      real(8) :: m_a = 0.0d0, m_b = 0.0d0   ! carrier masses, atom then anchor [amu]
      real(8) :: mu_c = 0.0d0               ! reduced mass [amu]
      real(8) :: qz_a(3) = 0.0d0, qz_b(3) = 0.0d0  ! equilibrium positions [Angstrom]
      real(8) :: u_bond(3) = 0.0d0          ! equilibrium bond unit direction (a -> b) [-]
      real(8) :: fr_a = 0.0d0, fr_b = 0.0d0 ! COM partition fractions m_b/m_tot, m_a/m_tot [-]
      real(8) :: om_c = 0.0d0               ! harmonic frequency [rad/(10 fs)]
      real(8) :: e_quant = 0.0d0            ! fixed quantum energy (mode 1) [internal]
      real(8) :: last_e = 0.0d0             ! energy of the last realization [internal]
      real(8) :: del_s = 0.0d0, vel_s = 0.0d0 ! the drawn stretch/velocity pair (post-draw record) [Ang], [Ang/(10 fs)]
   end type carry_t
   type(carry_t), allocatable :: carry(:)   ! carrier rows (list_atoms order)
   real(8) :: k_int = 0.0d0                 ! force constant [internal/Ang^2]
   real(8) :: t_c = 0.0d0                   ! validated oscillator temperature (mode 0) [K]
   integer :: n_osc_c = 0                   ! validated energy-law selector [-]
   logical :: inited = .false.              ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! surface_oscillator_init(p) - validate the parameters, bind the
   !                              atom-anchor carrier, precompute the
   !                              harmonic constants
   subroutine surface_oscillator_init(p)
      type(surf_osc_params_t), intent(in) :: p  ! member parameters
      real(8) :: r0(3), r_len, m_tot
      integer :: i, kc, code, n_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%n_osc /= 0 .and. p%n_osc /= 1) then
         call stop_soi('surface_oscillator_init', 'energy-law selector n_osc must be 0 (thermal) '// &
            'or 1 (quantum level), got', p%n_osc)
      end if
      if (p%n_osc == 0 .and. p%t_osc <= 0.0d0) then
         call stop_soiz('surface_oscillator_init', 'oscillator temperature t_osc [K] must be '// &
            'positive in the thermal mode, got', p%t_osc)
      end if
      if (p%n_osc == 1 .and. p%n_level < 0) then
         call stop_soi('surface_oscillator_init', 'vibrational quantum number n_level must be '// &
            '>= 0, got', p%n_level)
      end if
      if (p%k_osc <= 0.0d0) then
         call stop_soiz('surface_oscillator_init', 'force constant k_osc [kcal/mol/Ang^2] must '// &
            'be positive (a silent substitute is not drawn), got', p%k_osc)
      end if

      ! 2. fragment binding (assembly-order and carrier guards first)
      if (.not. allocated(list_atoms%frag)) then
         call stop_sox('surface_oscillator_init', 'no list_atoms fragment table - the atom list must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_sox('surface_oscillator_init', 'no per-fragment scheme array '// &
            '(reactants%dist_scheme unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_soi2('surface_oscillator_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('surface_oscillator')
      if (code == 0) then
         call stop_sox('surface_oscillator_init', 'member word "surface_oscillator" is not in the '// &
            'scheme-code table (incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_sox('surface_oscillator_init', 'no list_atoms fragment carries the '// &
            'surface_oscillator scheme - this member would never be dispatched')
      end if
      k_int = p%k_osc*e_conv
      n_osc_c = p%n_osc
      if (allocated(carry)) deallocate (carry)
      allocate (carry(n_carry))
      kc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         kc = kc + 1
         associate (fr => list_atoms%frag(i), cr => carry(kc))
            cr%frag = i
            if (fr%nat /= 2) then
               call stop_soi('surface_oscillator_init', 'the surface_oscillator carrier must be '// &
                  'exactly a two-atom atom-anchor pair (one bond); atom count =', fr%nat)
            end if
            allocate (cr%list(2))
            cr%list = fr%list
            cr%m_a = list_atoms%mass(cr%list(1))
            cr%m_b = list_atoms%mass(cr%list(2))
            m_tot = cr%m_a + cr%m_b
            cr%mu_c = cr%m_a*cr%m_b/m_tot
            cr%qz_a = fr%qz(1:3)
            cr%qz_b = fr%qz(4:6)
            r0 = cr%qz_b - cr%qz_a
            r_len = norm2(r0)
            if (r_len <= 1.0d-10) then
               call stop_soiz('surface_oscillator_init', 'degenerate carrier bond - equilibrium '// &
                  'bond length [Ang] =', r_len)
            end if
            cr%u_bond = r0/r_len
            cr%fr_a = cr%m_b/m_tot                ! position partition: r_a shifts by -fr_a*delta
            cr%fr_b = cr%m_a/m_tot                ! r_b shifts by +fr_b*delta
            cr%om_c = sqrt(k_int/cr%mu_c)
            cr%e_quant = (dble(p%n_level) + 0.5d0)*hbar_code*cr%om_c
         end associate
      end do
      t_c = p%t_osc
      inited = .true.
   end subroutine surface_oscillator_init

   !------------------------------------------------------------------
   ! surface_oscillator_sample(sta) - one realization per carrier (list_atoms
   !                                  order) on the energy shell
   !------------------------------------------------------------------
   ! surface_oscillator_draw() - per carrier: one energy draw (thermal
   !                              exponential or the fixed quantum level) and
   !                              one uniform phase; the derived stretch/
   !                              velocity pair lands in the record; no state
   !                              access
   !------------------------------------------------------------------
   subroutine surface_oscillator_draw()
      real(8) :: u, e_draw, amp, ph
      integer :: kc
      if (.not. inited) then
         call stop_sox('surface_oscillator_draw', 'draw called before surface_oscillator_init '// &
            '- the member carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         if (n_osc_c == 0) then
            u = rng_u()
            if (u <= 0.0d0) u = 1.0d-300          ! degenerate-draw clamp (idle in practice)
            e_draw = -kb_code*t_c*log(u)
         else
            e_draw = cr%e_quant
         end if
         amp = sqrt(2.0d0*e_draw/k_int)
         ph = two_pi*rng_u()
         cr%del_s = amp*cos(ph)
         cr%vel_s = -cr%om_c*amp*sin(ph)
         cr%last_e = e_draw
         end associate
      end do
   end subroutine surface_oscillator_draw

   !------------------------------------------------------------------
   ! surface_oscillator_realize(sta) - write the drawn stretch/velocity about
   !                                   the fixed COM (no RNG draws)
   !------------------------------------------------------------------
   subroutine surface_oscillator_realize(sta)
      type(state_t), intent(inout) :: sta  ! initial state (every carrier's q/p
                                           ! components are overwritten, list_atoms order)
      integer :: kc, i3
      if (.not. inited) then
         call stop_sox('surface_oscillator_realize', 'realize called before surface_oscillator_init '// &
            '- the member carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         if (size(sta%q) < 3*maxval(cr%list) .or. size(sta%p) < 3*maxval(cr%list)) then
            call stop_sox('surface_oscillator_realize', 'state too small for the carrier fragment '// &
               'atoms - state_create must follow the atom list')
         end if
         i3 = 3*(cr%list(1) - 1)
         sta%q(i3+1:i3+3) = cr%qz_a - cr%fr_a*cr%del_s*cr%u_bond
         sta%p(i3+1:i3+3) = -cr%mu_c*cr%vel_s*cr%u_bond
         i3 = 3*(cr%list(2) - 1)
         sta%q(i3+1:i3+3) = cr%qz_b + cr%fr_b*cr%del_s*cr%u_bond
         sta%p(i3+1:i3+3) = cr%mu_c*cr%vel_s*cr%u_bond
         end associate
      end do
   end subroutine surface_oscillator_realize

   !------------------------------------------------------------------
   ! surface_oscillator_sample(sta) - fused convenience: draw then realize
   !                                  (the check-program entry)
   !------------------------------------------------------------------
   subroutine surface_oscillator_sample(sta)
      type(state_t), intent(inout) :: sta
      call surface_oscillator_draw()
      call surface_oscillator_realize(sta)
   end subroutine surface_oscillator_sample

   !------------------------------------------------------------------
   ! surface_oscillator_last_e(k) - energy of carrier k's most recent
   !                                realization (k optional, default 1)
   real(8) function surface_oscillator_last_e(k)
      integer, intent(in), optional :: k
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      surface_oscillator_last_e = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) surface_oscillator_last_e = carry(kc)%last_e
   end function surface_oscillator_last_e

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; STOP 1)
   subroutine stop_sox(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_sox

   subroutine stop_soiz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_soiz

   subroutine stop_soi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_soi

   subroutine stop_soi2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_soi2

end module samp_surface_oscillator
