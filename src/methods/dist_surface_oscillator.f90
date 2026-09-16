!=====================================================================
! dist_surface_oscillator.f90 - distribution member: analytic harmonic
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
module dist_surface_oscillator
   use consts,        only: hbar_code, kb_code, e_conv, two_pi
   use rng,           only: rng_u
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use sampler,       only: samp_code
   implicit none
   private
   public :: surf_osc_params_t, surface_oscillator_init, surface_oscillator_sample, &
             surface_oscillator_last_e

   ! Member-owned parameters
   type :: surf_osc_params_t
      integer :: n_osc = 0         ! energy-law selector: 0 thermal, 1 quantum level [-]
      real(8) :: t_osc = 300.0d0   ! oscillator temperature [K] (mode 0)
      integer :: n_level = 0       ! vibrational quantum number n [-] (mode 1)
      real(8) :: k_osc = 30.0d0    ! harmonic force constant [kcal/mol/Ang^2]
   end type

   ! Derived member state (module-private; rebuilt by every surface_oscillator_init)
   integer :: frag_served = 0               ! bound roster fragment index [index]
   integer, allocatable :: list_c(:)        ! carrier atom index list [global atom index]
   real(8) :: m_a = 0.0d0, m_b = 0.0d0      ! carrier masses, atom then anchor [amu]
   real(8) :: mu_c = 0.0d0                  ! reduced mass [amu]
   real(8) :: qz_a(3) = 0.0d0, qz_b(3) = 0.0d0  ! equilibrium positions [Angstrom]
   real(8) :: u_bond(3) = 0.0d0             ! equilibrium bond unit direction (a -> b) [-]
   real(8) :: fr_a = 0.0d0, fr_b = 0.0d0    ! COM partition fractions m_b/m_tot, m_a/m_tot [-]
   real(8) :: k_int = 0.0d0                 ! force constant [internal/Ang^2]
   real(8) :: om_c = 0.0d0                  ! harmonic frequency [rad/(10 fs)]
   real(8) :: t_c = 0.0d0                   ! validated oscillator temperature (mode 0) [K]
   real(8) :: e_quant = 0.0d0               ! fixed quantum energy (mode 1) [internal]
   integer :: n_osc_c = 0                   ! validated energy-law selector [-]
   real(8) :: last_e = 0.0d0                ! energy of the last realization [internal]
   logical :: inited = .false.              ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! surface_oscillator_init(p) - validate the parameters, bind the
   !                              atom-anchor carrier, precompute the
   !                              harmonic constants
   subroutine surface_oscillator_init(p)
      type(surf_osc_params_t), intent(in) :: p  ! member parameters
      real(8) :: r0(3), r_len, m_tot
      integer :: i, code, n_carry, i_carry

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
      if (.not. allocated(roster%frag)) then
         call stop_sox('surface_oscillator_init', 'no roster fragment table - the roster must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_sox('surface_oscillator_init', 'no per-fragment scheme array '// &
            '(reactants%dist_scheme unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_soi2('surface_oscillator_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('surface_oscillator')
      if (code == 0) then
         call stop_sox('surface_oscillator_init', 'member word "surface_oscillator" is not in the '// &
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
         call stop_sox('surface_oscillator_init', 'no roster fragment carries the '// &
            'surface_oscillator scheme - this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_soi2('surface_oscillator_init', 'two or more fragments carry the '// &
            'surface_oscillator scheme, but one member instance carries one parameter set; '// &
            'count =', n_carry, 'first carrier =', i_carry)
      end if
      frag_served = i_carry
      associate (fr => roster%frag(i_carry))
         if (fr%nat /= 2) then
            call stop_soi('surface_oscillator_init', 'the surface_oscillator carrier must be '// &
               'exactly a two-atom atom-anchor pair (one bond); atom count =', fr%nat)
         end if
         if (allocated(list_c)) deallocate (list_c)
         allocate (list_c(2))
         list_c = fr%list
         m_a = roster%mass(list_c(1))
         m_b = roster%mass(list_c(2))
         m_tot = m_a + m_b
         mu_c = m_a*m_b/m_tot
         qz_a = fr%qz(1:3)
         qz_b = fr%qz(4:6)
         r0 = qz_b - qz_a
         r_len = norm2(r0)
         if (r_len <= 1.0d-10) then
            call stop_soiz('surface_oscillator_init', 'degenerate carrier bond - equilibrium '// &
            'bond length [Ang] =', r_len)
         end if
         u_bond = r0/r_len
         fr_a = m_b/m_tot                      ! position partition: r_a shifts by -fr_a*delta
         fr_b = m_a/m_tot                      ! r_b shifts by +fr_b*delta
         k_int = p%k_osc*e_conv
         om_c = sqrt(k_int/mu_c)
         n_osc_c = p%n_osc
         t_c = p%t_osc
         e_quant = (dble(p%n_level) + 0.5d0)*hbar_code*om_c
      end associate
      inited = .true.
   end subroutine surface_oscillator_init

   !------------------------------------------------------------------
   ! surface_oscillator_sample(sta) - one realization on the energy shell
   subroutine surface_oscillator_sample(sta)
      type(state_t), intent(inout) :: sta  ! initial state (the carrier's q/p
                                           ! components are overwritten)
      real(8) :: u, e_draw, amp, ph, del, vel
      integer :: i3

      if (.not. inited) then
         call stop_sox('surface_oscillator_sample', 'sample called before surface_oscillator_init '// &
            '- the member carries no parameter set')
      end if
      if (size(sta%q) < 3*maxval(list_c) .or. size(sta%p) < 3*maxval(list_c)) then
         call stop_sox('surface_oscillator_sample', 'state too small for the carrier fragment '// &
            'atoms - state_create must follow the roster')
      end if

      ! 1. energy: thermal exponential draw or the fixed quantum level
      if (n_osc_c == 0) then
         u = rng_u()
         if (u <= 0.0d0) u = 1.0d-300          ! degenerate-draw clamp (idle in practice)
         e_draw = -kb_code*t_c*log(u)
      else
         e_draw = e_quant
      end if

      ! 2./3. amplitude and uniform phase
      amp = sqrt(2.0d0*e_draw/k_int)
      ph = two_pi*rng_u()
      del = amp*cos(ph)
      vel = -om_c*amp*sin(ph)

      ! 4. positions about the fixed COM, momenta along the bond
      i3 = 3*(list_c(1) - 1)
      sta%q(i3+1:i3+3) = qz_a - fr_a*del*u_bond
      sta%p(i3+1:i3+3) = -mu_c*vel*u_bond
      i3 = 3*(list_c(2) - 1)
      sta%q(i3+1:i3+3) = qz_b + fr_b*del*u_bond
      sta%p(i3+1:i3+3) = mu_c*vel*u_bond
      last_e = e_draw
   end subroutine surface_oscillator_sample

   !------------------------------------------------------------------
   ! surface_oscillator_last_e() - energy of the most recent realization
   real(8) function surface_oscillator_last_e()
      surface_oscillator_last_e = last_e
   end function surface_oscillator_last_e

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_sox(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_sox

   subroutine stop_soiz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_soiz

   subroutine stop_soi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_soi

   subroutine stop_soi2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_soi2

end module dist_surface_oscillator
