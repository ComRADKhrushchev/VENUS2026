!=====================================================================
! dist_ebk.f90 - distribution member (methods/ library): EBK semiclassical
!                state distribution - a fixed-(n,J) diatomic initial state
!                from Morse spectroscopic constants
! Contract:
!   ebk_init(p) - validate the member parameters, derive the Morse curve
!     from the spectroscopic constants, locate the (n,J) energy level by
!     WKB action quantization, and bind the served fragment. Loud-stops
!     name the offending field/value: n_vib < 0; j_rot < 0; w_e <= 0;
!     w_ex_e <= 0; b_rot < 0; ai_rot <= 0 (retained for the orientation
!     members, validated on arrival); a carrier fragment that is not a
!     two-atom roster row; a degenerate equilibrium geometry (zero bond
!     length); a level that admits no classical region (at/above the
!     dissociation asymptote or below the effective-potential minimum);
!     a fixed-point iteration that fails to converge within its cap.
!     Binding: the served fragment is the roster fragment whose
!     dist_scheme row selects THIS member (exactly one carrier required).
!   ebk_sample(sta) - one draw: internuclear distance r with WKB
!     probability density prop. to 1/|pr| on the inset turning-point
!     interval (exact rejection: envelope constant = the smaller endpoint
!     momentum, positive on the inset interval), radial momentum pr with a
!     random sign, then the diatomic assembly writes the carrier
!     fragment's q/p components: COM at the origin, zero total momentum,
!     bond along the equilibrium axis direction, |L| = sqrt(J(J+1))*hbar
!     perpendicular to the bond with a random azimuth. Translation (where
!     the fragment sits, how fast it flies) belongs to the incident
!     channel - this member owns the internal degrees of freedom only.
!   Getters (copies): ebk_e_level (converged level energy [internal]),
!     ebk_r_turn (the sampled interval [Angstrom]), ebk_morse (the derived
!     curve parameters a [1/Angstrom], D_e [internal], reduced mass
!     [amu], equilibrium bond length [Angstrom], rotational moment
!     [internal angular momentum unit]), ebk_last (the last draw's r and
!     signed pr).
! Physics of the derivation (assembly-time mechanism note): the Morse
!   curve is DERIVED from the spectroscopic constants of the parameter
!   type - D_e = w_e^2/(4*w_ex_e) [cm^-1], and the range parameter from
!   the requirement that the WKB spectrum reproduce the constants:
!   a = sqrt(2*mu*e_conv*wvn2e*w_ex_e)/hbar [1/Angstrom]. With this
!   construction the Bohr-Sommerfeld action of the pure Morse is exactly
!   quantized (the J = 0 levels are EXACTLY w_e*(n+1/2) - w_ex_e*(n+1/2)^2
!   - the anchor the member's criteria program rides). The rotational
!   term b_rot*J*(J+1) enters the level seed; the level is then re-located
!   by fixed-point quantization of the full effective potential (Morse
!   plus centrifugal al^2/(2*mu*r^2)) - e <- e + (x - action/pi*hbar)*hnu,
!   hnu = e_seed/x, tolerance 1e-6, cap 1000 iterations. Frequencies may
!   equivalently arrive from the Hessian eigenvalue chain at assembly
!   time; through a no-op chain nothing arrives and the parameter
!   validation is the loud-stop (fail-loud degradation posture).
! Units: w_e/w_ex_e/b_rot [cm^-1], n_vib/j_rot [-], ai_rot [amu*Ang^2];
!   energies internal (= e_conv*wvn2e*[cm^-1]); coordinates [Angstrom];
!   momenta [amu*Ang/(10 fs)]; hbar = hbar_code [internal angular
!   momentum unit]. Turning-point inset 0.001 [Angstrom] (both the
!   quantization integral and the sampling interval use it - one
!   consistent convention).
!=====================================================================
module dist_ebk
   use consts,        only: pi, two_pi, hbar_code, e_conv, wvn2e
   use rng,           only: rng_u
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use sampler,       only: samp_code
   use geometry,      only: gl_nodes
   implicit none
   private
   public :: ebk_params_t, ebk_init, ebk_sample, ebk_e_level, ebk_r_turn, &
             ebk_morse, ebk_last

   real(8), parameter :: eps_turn = 1.0d-3   ! turning-point inset [Angstrom]
   integer, parameter :: n_quad = 50         ! Gauss-Legendre nodes for the action [-]
   integer, parameter :: max_iter = 1000     ! fixed-point iteration cap [-]
   integer, parameter :: max_reject = 100000 ! rejection-sampling draw cap [-]
   real(8), parameter :: tol_quant = 1.0d-6  ! quantization tolerance [quantum]

   ! Member-owned parameters (reduced mass, bond length and orientation are
   ! roster-derived; the dissociation energy is DERIVED, not stored)
   type :: ebk_params_t
      integer :: n_vib  = 0        ! fixed vibrational quantum number n [-]
      integer :: j_rot  = 0        ! fixed rotational quantum number J [-]
      real(8) :: w_e    = 0.0d0    ! spectroscopic constant w_e [cm^-1]
      real(8) :: w_ex_e = 0.0d0    ! anharmonic constant w_ex_e [cm^-1]
      real(8) :: b_rot  = 0.0d0    ! rotational constant B_e [cm^-1]
      real(8) :: ai_rot = 0.0d0    ! moment of inertia [amu*Ang^2] (orientation-member input)
   end type

   ! Derived member state (module-private; rebuilt by every ebk_init)
   real(8) :: a_morse = 0.0d0     ! Morse range parameter [1/Angstrom]
   real(8) :: d_e = 0.0d0         ! Morse depth from the well bottom [internal]
   real(8) :: mu_red = 0.0d0      ! reduced mass [amu]
   real(8) :: r_eq = 0.0d0        ! equilibrium bond length [Angstrom]
   real(8) :: al_mom = 0.0d0      ! rotational moment |L| [hbar_code units]
   real(8) :: e_lvl = 0.0d0       ! converged level energy [internal]
   real(8) :: r_lo = 0.0d0        ! sampled interval lower bound (inset) [Angstrom]
   real(8) :: r_hi = 0.0d0        ! sampled interval upper bound (inset) [Angstrom]
   real(8) :: p_env = 0.0d0       ! rejection envelope = smaller endpoint momentum [-]
   real(8) :: ax_bond(3) = 0.0d0  ! equilibrium bond axis unit vector [-]
   real(8) :: e_perp1(3) = 0.0d0  ! bond-perpendicular triad, first axis [-]
   real(8) :: e_perp2(3) = 0.0d0  ! bond-perpendicular triad, second axis [-]
   integer :: at1 = 0, at2 = 0    ! carrier fragment atom indices [global index]
   real(8) :: last_r = 0.0d0      ! last drawn internuclear distance [Angstrom]
   real(8) :: last_pr = 0.0d0     ! last drawn signed radial momentum [-]
   logical :: inited = .false.    ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! ebk_init(p) - validate, derive the Morse curve, quantize the level,
   !               bind the carrier fragment
   subroutine ebk_init(p)
      type(ebk_params_t), intent(in) :: p ! member parameters (quantum numbers /
                                          ! spectroscopic constants)
      real(8) :: node(n_quad), wts(n_quad)
      real(8) :: x, e_seed, hnu, an, dum, dumv
      real(8) :: m1, m2, dv(3), big(3), small
      integer :: i, code, n_carry, i_carry, it, k_small

      ! 1. parameter validation (fail-loud, naming the field)
      if (p%n_vib < 0) call stop_ebki('ebk_init', 'vibrational quantum number n_vib must be >= 0, got', p%n_vib)
      if (p%j_rot < 0) call stop_ebki('ebk_init', 'rotational quantum number j_rot must be >= 0, got', p%j_rot)
      if (p%w_e <= 0.0d0) then
         dumv = p%w_e
         call stop_ebkz('ebk_init', 'spectroscopic constant w_e [cm^-1] must be positive, got', dumv)
      end if
      if (p%w_ex_e <= 0.0d0) then
         dumv = p%w_ex_e
         call stop_ebkz('ebk_init', 'anharmonic constant w_ex_e [cm^-1] must be positive, got', dumv)
      end if
      if (p%b_rot < 0.0d0) then
         dumv = p%b_rot
         call stop_ebkz('ebk_init', 'rotational constant b_rot [cm^-1] must be >= 0, got', dumv)
      end if
      if (p%ai_rot <= 0.0d0) then
         dumv = p%ai_rot
         call stop_ebkz('ebk_init', 'moment of inertia ai_rot [amu*Ang^2] must be positive, got', dumv)
      end if

      ! 2. fragment binding (assembly-order and carrier guards, then the
      !    roster-derived diatomic data)
      if (.not. allocated(roster%frag)) then
         call stop_ebkx('ebk_init', 'no roster fragment table - the roster must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_ebkx('ebk_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_ebki2('ebk_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('ebk')
      if (code == 0) then
         call stop_ebkx('ebk_init', 'member word "ebk" is not in the scheme-code table '// &
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
         call stop_ebkx('ebk_init', 'no roster fragment carries the ebk scheme - '// &
            'this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_ebki2('ebk_init', 'two or more fragments carry the ebk scheme, but one '// &
            'member instance carries one parameter set; count =', n_carry, 'first carrier =', i_carry)
      end if
      associate (fr => roster%frag(i_carry))
         if (fr%nat /= 2) then
            call stop_ebki('ebk_init', 'the ebk carrier must be a diatomic fragment; fragment '// &
               'atom count =', fr%nat)
         end if
         at1 = fr%list(1)
         at2 = fr%list(2)
         m1 = roster%mass(at1)
         m2 = roster%mass(at2)
         mu_red = m1*m2/(m1 + m2)
         dv = fr%qz(4:6) - fr%qz(1:3)
         r_eq = norm2(dv)
         if (r_eq <= 0.0d0) then
            dumv = r_eq
            call stop_ebkz('ebk_init', 'degenerate carrier geometry: equilibrium bond length ', dumv)
         end if
         ax_bond = dv/r_eq
      end associate

      ! bond-perpendicular triad: perturb the smallest axis component, orthonormalize
      big = abs(ax_bond)
      k_small = minloc(big, dim=1)
      dv = 0.0d0
      dv(k_small) = 1.0d0
      e_perp1 = dv - dot_product(dv, ax_bond)*ax_bond
      e_perp1 = e_perp1/norm2(e_perp1)
      e_perp2 = cross_p(ax_bond, e_perp1)

      ! 3. Morse derivation from the spectroscopic constants (see the module
      !    header for the construction and its WKB-exactness anchor)
      d_e = e_conv*wvn2e*(p%w_e*p%w_e/(4.0d0*p%w_ex_e))
      a_morse = sqrt(2.0d0*mu_red*(e_conv*wvn2e*p%w_ex_e))/hbar_code
      al_mom = sqrt(dble(p%j_rot*(p%j_rot + 1)))*hbar_code

      ! 4. level seed and WKB quantization by fixed point over the effective
      !    potential (Morse from the well bottom + centrifugal)
      x = dble(p%n_vib) + 0.5d0
      e_seed = e_conv*wvn2e*(p%w_e*x - p%w_ex_e*x*x + p%b_rot*dble(p%j_rot*(p%j_rot + 1)))
      hnu = e_seed/x
      e_lvl = e_seed
      do it = 1, max_iter
         call quantize_step(e_lvl, node, wts, an)
         dum = x - an
         if (abs(dum) <= tol_quant) exit
         e_lvl = e_lvl + dum*hnu
         if (it == max_iter) then
            call stop_ebkx('ebk_init', 'WKB quantization failed to converge within the '// &
               'iteration cap (level energy drifting without settling)')
         end if
      end do

      ! 5. final turning points, inset interval, envelope
      call turn_bounds(e_lvl, r_lo, r_hi)
      r_lo = r_lo + eps_turn
      r_hi = r_hi - eps_turn
      if (r_hi <= r_lo) then
         call stop_ebkx('ebk_init', 'the turning-point inset eats the whole classical region '// &
            '(classical interval narrower than the inset)')
      end if
      small = p_abs(r_lo)
      dum = p_abs(r_hi)
      if (dum < small) small = dum
      p_env = small
      if (p_env <= 0.0d0) then
         call stop_ebkx('ebk_init', 'non-positive rejection envelope on the inset interval '// &
            '(turning-point location failure)')
      end if
      last_r = 0.0d0
      last_pr = 0.0d0
      inited = .true.
   end subroutine ebk_init

   !------------------------------------------------------------------
   ! ebk_sample(sta) - draw (r, sign(pr)) and write the carrier fragment's
   !                   internal q/p components
   subroutine ebk_sample(sta)
      type(state_t), intent(inout) :: sta  ! initial state (the carrier fragment's q/p
                                           ! components are overwritten)
      real(8) :: r, pr_mag, pr, delta, om(3), m1, m2, mt
      real(8) :: q1v(3), q2v(3), w1(3), w2(3)
      integer :: j, n_try

      if (.not. inited) then
         call stop_ebkx('ebk_sample', 'sample called before ebk_init - the member carries '// &
            'no parameter set')
      end if
      if (size(sta%q) < 3*max(at1, at2) .or. size(sta%p) < 3*max(at1, at2)) then
         call stop_ebkx('ebk_sample', 'state too small for the carrier fragment atoms - '// &
            'state_create must follow the roster')
      end if

      ! rejection: r uniform on the inset interval, accept prob p_env/|pr|
      ! (exact: p_env is the smaller endpoint momentum, the |pr| minimum there)
      do n_try = 1, max_reject
         r = r_lo + (r_hi - r_lo)*rng_u()
         pr_mag = p_abs(r)
         if (rng_u()*pr_mag <= p_env) exit
         if (n_try == max_reject) then
            call stop_ebkx('ebk_sample', 'rejection sampling exceeded its draw cap (the '// &
               'density ratio is inconsistent with the envelope)')
         end if
      end do
      if (rng_u() < 0.5d0) then
         pr = -pr_mag
      else
         pr = pr_mag
      end if
      last_r = r
      last_pr = pr

      ! diatomic assembly along the equilibrium bond axis: COM at the origin,
      ! zero total momentum, radial momentum pr, moment L perpendicular
      m1 = roster%mass(at1)
      m2 = roster%mass(at2)
      mt = m1 + m2
      delta = two_pi*rng_u()
      om = al_mom/(mu_red*r*r)*(cos(delta)*e_perp1 + sin(delta)*e_perp2)
      q1v = -r*(m2/mt)*ax_bond
      q2v = r*(m1/mt)*ax_bond
      w1 = m1*cross_p(om, q1v)
      w2 = m2*cross_p(om, q2v)
      do j = 1, 3
         sta%q(3*at1 - 3 + j) = q1v(j)
         sta%q(3*at2 - 3 + j) = q2v(j)
         sta%p(3*at1 - 3 + j) = -pr*ax_bond(j) + w1(j)
         sta%p(3*at2 - 3 + j) = pr*ax_bond(j) + w2(j)
      end do
   end subroutine ebk_sample

   !------------------------------------------------------------------
   ! getters (copies of the derived member state)
   !------------------------------------------------------------------
   pure real(8) function ebk_e_level()
      ebk_e_level = e_lvl
   end function ebk_e_level

   subroutine ebk_r_turn(r1, r2)
      real(8), intent(out) :: r1, r2
      r1 = r_lo
      r2 = r_hi
   end subroutine ebk_r_turn

   subroutine ebk_morse(a_out, d_out, mu_out, re_out, al_out)
      real(8), intent(out) :: a_out, d_out, mu_out, re_out, al_out
      a_out = a_morse
      d_out = d_e
      mu_out = mu_red
      re_out = r_eq
      al_out = al_mom
   end subroutine ebk_morse

   subroutine ebk_last(r_out, pr_out)
      real(8), intent(out) :: r_out, pr_out
      r_out = last_r
      pr_out = last_pr
   end subroutine ebk_last

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! v_eff(r) - effective potential: Morse from the well bottom + centrifugal
   pure real(8) function v_eff(r)
      real(8), intent(in) :: r
      v_eff = d_e*(1.0d0 - exp(-a_morse*(r - r_eq)))**2 + al_mom**2/(2.0d0*mu_red*r*r)
   end function v_eff

   ! p_abs(r) - radial momentum magnitude sqrt(2*mu*(e - V_eff)) (0 when
   !             outside the classical region - defensive, never sampled there)
   pure real(8) function p_abs(r)
      real(8), intent(in) :: r
      p_abs = sqrt(max(0.0d0, 2.0d0*mu_red*(e_lvl - v_eff(r))))
   end function p_abs

   ! quantize_step(e, node, wts, an) - one action evaluation: turning points
   !   at energy e, inset, 50-node Gauss-Legendre, an = action/(pi*hbar)
   subroutine quantize_step(e, node, wts, an)
      real(8), intent(in) :: e
      real(8), intent(out) :: node(n_quad), wts(n_quad), an
      real(8) :: r1, r2, asum
      integer :: i
      call turn_bounds(e, r1, r2)
      r1 = r1 + eps_turn
      r2 = r2 - eps_turn
      call gl_nodes(r1, r2, node, wts)
      asum = 0.0d0
      do i = 1, n_quad
         asum = asum + wts(i)*sqrt(max(0.0d0, 2.0d0*mu_red*(e - v_eff(node(i)))))
      end do
      an = asum/(pi*hbar_code)
   end subroutine quantize_step

   ! turn_bounds(e, r1, r2) - classical interval at energy e: coarse scan for
   !   the effective-potential minimum, walk out to brackets, bisect both sides
   subroutine turn_bounds(e, r1, r2)
      real(8), intent(in) :: e
      real(8), intent(out) :: r1, r2
      integer, parameter :: n_scan = 1000
      real(8) :: r_well, v_well, r_try, lo, hi
      integer :: i
      if (e >= d_e) then
         call stop_ebkz('ebk_init', 'level energy at or above the dissociation asymptote - '// &
            'no outer turning point exists; level [internal] =', e)
      end if
      r_well = r_eq
      v_well = v_eff(r_eq)
      do i = 1, n_scan
         r_try = 0.1d0*r_eq + (5.0d0*r_eq - 0.1d0*r_eq)*dble(i - 1)/dble(n_scan - 1)
         if (v_eff(r_try) < v_well) then
            v_well = v_eff(r_try)
            r_well = r_try
         end if
      end do
      if (v_well >= e) then
         call stop_ebkz('ebk_init', 'level energy below the effective-potential minimum - '// &
            'no classical region exists; level [internal] =', e)
      end if
      lo = r_well
      do i = 1, 60
         lo = lo/1.5d0
         if (v_eff(lo) > e) exit
         if (i == 60) call stop_ebkx('ebk_init', 'no inner turning point (inner wall not rising)')
      end do
      r1 = bisect(lo, r_well, e)
      hi = r_well
      do i = 1, 100
         hi = 1.5d0*hi + 0.1d0
         if (v_eff(hi) > e) exit
         if (i == 100) call stop_ebkx('ebk_init', 'no outer turning point below the asymptote')
      end do
      r2 = bisect(r_well, hi, e)
   end subroutine turn_bounds

   ! bisect - root of V_eff(r) = e on [a, b] (V_eff(a) > e > V_eff(b) or vice
   !          versa; 200 iterations drive the bracket under 1e-12 relative)
   real(8) function bisect(a, b, e)
      real(8), intent(in) :: a, b, e
      real(8) :: lo, hi, mid
      integer :: i
      lo = a
      hi = b
      do i = 1, 200
         mid = 0.5d0*(lo + hi)
         if (v_eff(mid) > e) then
            if (v_eff(lo) > e) then
               lo = mid
            else
               hi = mid
            end if
         else
            if (v_eff(lo) > e) then
               hi = mid
            else
               lo = mid
            end if
         end if
         if (hi - lo < 1.0d-12*max(1.0d0, hi)) exit
      end do
      bisect = 0.5d0*(lo + hi)
   end function bisect

   ! cross_p - vector cross product (pure, allocatable-free for the hot path)
   pure function cross_p(a, b) result(c)
      real(8), intent(in) :: a(3), b(3)
      real(8) :: c(3)
      c = (/ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) /)
   end function cross_p

   ! loud-stop helpers (name the caller; STOP 1)
   subroutine stop_ebkx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_ebkx

   subroutine stop_ebkz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_ebkz

   subroutine stop_ebki(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_ebki

   subroutine stop_ebki2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_ebki2

end module dist_ebk
