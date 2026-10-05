!=====================================================================
! samp_ebk.f90 - distribution member: EBK semiclassical fixed-(n,J)
!   diatomic initial state from Morse spectroscopic constants
! Design:
!   One parameter set serves every carrier equally; the Morse curve,
!   the quantized level, the sampled interval and the bond triad are
!   per-carrier derived state (each carrier is its own diatomic). The
!   Morse curve is derived so the pure-Morse WKB spectrum is exactly
!   quantized; the level is re-located by fixed-point quantization of
!   the full effective potential. Sample draws r with density prop. to
!   1/|pr| and writes the carrier's internal q/p only - translation is
!   the incident channel's business.
!   Units: w_e/w_ex_e/b_rot [cm^-1]; energies internal; q [Angstrom];
!   p [amu*Ang/(10 fs)]; hbar = hbar_code.
!=====================================================================
module samp_ebk
   use consts,        only: pi, two_pi, hbar_code, e_conv, wvn2e
   use rng,           only: rng_u
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code
   use geometry,      only: gl_nodes
   implicit none
   private
   public :: ebk_params_t, ebk_init, ebk_draw, ebk_realize, ebk_sample, ebk_e_level, ebk_r_turn, &
             ebk_morse, ebk_last

   real(8), parameter :: eps_turn = 1.0d-3   ! turning-point inset [Angstrom]
   integer, parameter :: n_quad = 50         ! Gauss-Legendre nodes for the action [-]
   integer, parameter :: max_iter = 1000     ! fixed-point iteration cap [-]
   integer, parameter :: max_reject = 100000 ! rejection-sampling draw cap [-]
   real(8), parameter :: tol_quant = 1.0d-6  ! quantization tolerance [quantum]

   ! Member-owned parameters (reduced mass, bond length and orientation are
   ! list_atoms-derived; the dissociation energy is DERIVED, not stored)
   type :: ebk_params_t
      integer :: n_vib  = 0        ! fixed vibrational quantum number n [-]
      integer :: j_rot  = 0        ! fixed rotational quantum number J [-]
      real(8) :: w_e    = 0.0d0    ! spectroscopic constant w_e [cm^-1]
      real(8) :: w_ex_e = 0.0d0    ! anharmonic constant w_ex_e [cm^-1]
      real(8) :: b_rot  = 0.0d0    ! rotational constant B_e [cm^-1]
      real(8) :: ai_rot = 0.0d0    ! moment of inertia [amu*Ang^2] (orientation-member input)
   end type

   ! Derived member state (module-private; rebuilt by every ebk_init) - one
   ! row per carrier: each carrier is its own diatomic with its own Morse
   ! curve, level and triad
   type :: carry_t
      integer :: at1 = 0, at2 = 0    ! carrier fragment atom indices [global index]
      real(8) :: mu_red = 0.0d0      ! reduced mass [amu]
      real(8) :: r_eq = 0.0d0        ! equilibrium bond length [Angstrom]
      real(8) :: ax_bond(3) = 0.0d0  ! equilibrium bond axis unit vector [-]
      real(8) :: e_perp1(3) = 0.0d0  ! bond-perpendicular triad, first axis [-]
      real(8) :: e_perp2(3) = 0.0d0  ! bond-perpendicular triad, second axis [-]
      real(8) :: a_morse = 0.0d0     ! Morse range parameter [1/Angstrom]
      real(8) :: d_e = 0.0d0         ! Morse depth from the well bottom [internal]
      real(8) :: al_mom = 0.0d0      ! rotational moment |L| [hbar_code units]
      real(8) :: e_lvl = 0.0d0       ! converged level energy [internal]
      real(8) :: r_lo = 0.0d0        ! sampled interval lower bound (inset) [Angstrom]
      real(8) :: r_hi = 0.0d0        ! sampled interval upper bound (inset) [Angstrom]
      real(8) :: p_env = 0.0d0       ! rejection envelope = smaller endpoint momentum [-]
      real(8) :: last_r = 0.0d0      ! last drawn internuclear distance [Angstrom]
      real(8) :: last_pr = 0.0d0     ! last drawn signed radial momentum [-]
      real(8) :: last_delta = 0.0d0    ! last drawn L-direction angle in the bond plane [rad]
   end type carry_t
   type(carry_t), allocatable :: carry(:)  ! carrier rows (list_atoms order)
   logical :: inited = .false.    ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! ebk_init(p) - validate, derive the Morse curve, quantize the level and
   !               bind every carrier fragment (one shared parameter set)
   subroutine ebk_init(p)
      type(ebk_params_t), intent(in) :: p ! member parameters (quantum numbers /
                                          ! spectroscopic constants)
      type(carry_t) :: cr
      real(8) :: node(n_quad), wts(n_quad)
      real(8) :: x, e_seed, hnu, an, dum, dumv
      real(8) :: m1, m2, dv(3), big(3), small
      integer :: i, kc, code, n_carry, it, k_small

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

      ! 2. carrier binding (assembly-order guards, then the atom list-derived
      !    diatomic data of every carrier; list_atoms order)
      if (.not. allocated(list_atoms%frag)) then
         call stop_ebkx('ebk_init', 'no list_atoms fragment table - the atom list must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_ebkx('ebk_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its carriers')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_ebki2('ebk_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('ebk')
      if (code == 0) then
         call stop_ebkx('ebk_init', 'member word "ebk" is not in the scheme-code table '// &
            '(incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_ebkx('ebk_init', 'no list_atoms fragment carries the ebk scheme - '// &
            'this member would never be dispatched')
      end if
      if (allocated(carry)) deallocate (carry)
      allocate (carry(n_carry))
      kc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         kc = kc + 1
         cr = carry_t()
         associate (fr => list_atoms%frag(i))
            if (fr%nat /= 2) then
               call stop_ebki('ebk_init', 'the ebk carrier must be a diatomic fragment; fragment '// &
                  'atom count =', fr%nat)
            end if
            cr%at1 = fr%list(1)
            cr%at2 = fr%list(2)
            m1 = list_atoms%mass(cr%at1)
            m2 = list_atoms%mass(cr%at2)
            cr%mu_red = m1*m2/(m1 + m2)
            dv = fr%qz(4:6) - fr%qz(1:3)
            cr%r_eq = norm2(dv)
            if (cr%r_eq <= 0.0d0) then
               dumv = cr%r_eq
               call stop_ebkz('ebk_init', 'degenerate carrier geometry: equilibrium bond length ', dumv)
            end if
            cr%ax_bond = dv/cr%r_eq
         end associate

         ! bond-perpendicular triad: perturb the smallest axis component, orthonormalize
         big = abs(cr%ax_bond)
         k_small = minloc(big, dim=1)
         dv = 0.0d0
         dv(k_small) = 1.0d0
         cr%e_perp1 = dv - dot_product(dv, cr%ax_bond)*cr%ax_bond
         cr%e_perp1 = cr%e_perp1/norm2(cr%e_perp1)
         cr%e_perp2 = cross_p(cr%ax_bond, cr%e_perp1)

         ! 3. Morse derivation from the spectroscopic constants (see the module
         !    header for the construction and its WKB-exactness anchor)
         cr%d_e = e_conv*wvn2e*(p%w_e*p%w_e/(4.0d0*p%w_ex_e))
         cr%a_morse = sqrt(2.0d0*cr%mu_red*(e_conv*wvn2e*p%w_ex_e))/hbar_code
         cr%al_mom = sqrt(dble(p%j_rot*(p%j_rot + 1)))*hbar_code

         ! 4. level seed and WKB quantization by fixed point over the effective
         !    potential (Morse from the well bottom + centrifugal)
         x = dble(p%n_vib) + 0.5d0
         e_seed = e_conv*wvn2e*(p%w_e*x - p%w_ex_e*x*x + p%b_rot*dble(p%j_rot*(p%j_rot + 1)))
         hnu = e_seed/x
         cr%e_lvl = e_seed
         do it = 1, max_iter
            call quantize_step(cr, cr%e_lvl, node, wts, an)
            dum = x - an
            if (abs(dum) <= tol_quant) exit
            cr%e_lvl = cr%e_lvl + dum*hnu
            if (it == max_iter) then
               call stop_ebkx('ebk_init', 'WKB quantization failed to converge within the '// &
                  'iteration cap (level energy drifting without settling)')
            end if
         end do

         ! 5. final turning points, inset interval, envelope
         call turn_bounds(cr, cr%e_lvl, cr%r_lo, cr%r_hi)
         cr%r_lo = cr%r_lo + eps_turn
         cr%r_hi = cr%r_hi - eps_turn
         if (cr%r_hi <= cr%r_lo) then
            call stop_ebkx('ebk_init', 'the turning-point inset eats the whole classical region '// &
               '(classical interval narrower than the inset)')
         end if
         small = p_abs(cr, cr%r_lo)
         dum = p_abs(cr, cr%r_hi)
         if (dum < small) small = dum
         cr%p_env = small
         if (cr%p_env <= 0.0d0) then
            call stop_ebkx('ebk_init', 'non-positive rejection envelope on the inset interval '// &
               '(turning-point location failure)')
         end if
         carry(kc) = cr
      end do
      inited = .true.
   end subroutine ebk_init

   !------------------------------------------------------------------
   ! ebk_draw() - per carrier: the rejection draw of r on the inset interval,
   !              the signed radial momentum (fresh uniform), then the
   !              L-direction angle in the bond plane - the historical order;
   !              no state access
   !------------------------------------------------------------------
   subroutine ebk_draw()
      real(8) :: r, pr_mag, pr
      integer :: kc, n_try
      if (.not. inited) then
         call stop_ebkx('ebk_draw', 'draw called before ebk_init - the member carries '// &
            'no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         do n_try = 1, max_reject
            r = cr%r_lo + (cr%r_hi - cr%r_lo)*rng_u()
            pr_mag = p_abs(cr, r)
            if (rng_u()*pr_mag <= cr%p_env) exit
            if (n_try == max_reject) then
               call stop_ebkx('ebk_draw', 'rejection sampling exceeded its draw cap (the '// &
                  'density ratio is inconsistent with the envelope)')
            end if
         end do
         if (rng_u() < 0.5d0) then
            pr = -pr_mag
         else
            pr = pr_mag
         end if
         cr%last_r = r
         cr%last_pr = pr
         cr%last_delta = two_pi*rng_u()
         end associate
      end do
   end subroutine ebk_draw

   !------------------------------------------------------------------
   ! ebk_realize(sta) - per carrier (list_atoms order): diatomic assembly
   !                    along the equilibrium bond axis from the drawn record
   !                    (COM at the origin, zero total momentum, radial
   !                    momentum pr, moment L perpendicular) - no RNG draws
   !------------------------------------------------------------------
   subroutine ebk_realize(sta)
      type(state_t), intent(inout) :: sta  ! initial state (every carrier fragment's
                                           ! q/p components are overwritten, list_atoms order)
      real(8) :: r, pr, delta, om(3), m1, m2, mt
      real(8) :: q1v(3), q2v(3), w1(3), w2(3)
      integer :: j, kc, at1, at2
      if (.not. inited) then
         call stop_ebkx('ebk_realize', 'realize called before ebk_init - the member carries '// &
            'no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         at1 = cr%at1
         at2 = cr%at2
         if (size(sta%q) < 3*max(at1, at2) .or. size(sta%p) < 3*max(at1, at2)) then
            call stop_ebkx('ebk_realize', 'state too small for the carrier fragment atoms - '// &
               'state_create must follow the atom list')
         end if
         r = cr%last_r
         pr = cr%last_pr
         delta = cr%last_delta
         m1 = list_atoms%mass(at1)
         m2 = list_atoms%mass(at2)
         mt = m1 + m2
         om = cr%al_mom/(cr%mu_red*r*r)*(cos(delta)*cr%e_perp1 + sin(delta)*cr%e_perp2)
         q1v = -r*(m2/mt)*cr%ax_bond
         q2v = r*(m1/mt)*cr%ax_bond
         w1 = m1*cross_p(om, q1v)
         w2 = m2*cross_p(om, q2v)
         do j = 1, 3
            sta%q(3*at1 - 3 + j) = q1v(j)
            sta%q(3*at2 - 3 + j) = q2v(j)
            sta%p(3*at1 - 3 + j) = -pr*cr%ax_bond(j) + w1(j)
            sta%p(3*at2 - 3 + j) = pr*cr%ax_bond(j) + w2(j)
         end do
         end associate
      end do
   end subroutine ebk_realize

   !------------------------------------------------------------------
   ! ebk_sample(sta) - fused convenience: draw then realize (the
   !                   check-program entry)
   !------------------------------------------------------------------
   subroutine ebk_sample(sta)
      type(state_t), intent(inout) :: sta
      call ebk_draw()
      call ebk_realize(sta)
   end subroutine ebk_sample

   !------------------------------------------------------------------
   ! getters (copies of the per-carrier derived state; k = carrier row)
   !------------------------------------------------------------------
   pure real(8) function ebk_e_level(k)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      ebk_e_level = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) ebk_e_level = carry(kc)%e_lvl
   end function ebk_e_level

   subroutine ebk_r_turn(k, r1, r2)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      real(8), intent(out) :: r1, r2
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      r1 = 0.0d0
      r2 = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) then
         r1 = carry(kc)%r_lo
         r2 = carry(kc)%r_hi
      end if
   end subroutine ebk_r_turn

   subroutine ebk_morse(k, a_out, d_out, mu_out, re_out, al_out)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      real(8), intent(out) :: a_out, d_out, mu_out, re_out, al_out
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      a_out = 0.0d0
      d_out = 0.0d0
      mu_out = 0.0d0
      re_out = 0.0d0
      al_out = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) then
         a_out = carry(kc)%a_morse
         d_out = carry(kc)%d_e
         mu_out = carry(kc)%mu_red
         re_out = carry(kc)%r_eq
         al_out = carry(kc)%al_mom
      end if
   end subroutine ebk_morse

   subroutine ebk_last(k, r_out, pr_out)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      real(8), intent(out) :: r_out, pr_out
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      r_out = 0.0d0
      pr_out = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) then
         r_out = carry(kc)%last_r
         pr_out = carry(kc)%last_pr
      end if
   end subroutine ebk_last

   !------------------------------------------------------------------
   ! private helpers (all carrier-row-aware: the curve constants live in the
   ! carrier row under evaluation)
   !------------------------------------------------------------------

   ! v_eff(cr, r) - effective potential: Morse from the well bottom + centrifugal
   pure real(8) function v_eff(cr, r)
      type(carry_t), intent(in) :: cr
      real(8), intent(in) :: r
      v_eff = cr%d_e*(1.0d0 - exp(-cr%a_morse*(r - cr%r_eq)))**2 &
              + cr%al_mom**2/(2.0d0*cr%mu_red*r*r)
   end function v_eff

   ! p_abs(cr, r) - radial momentum magnitude sqrt(2*mu*(e - V_eff)) (0 when
   !                 outside the classical region - defensive, never sampled there)
   pure real(8) function p_abs(cr, r)
      type(carry_t), intent(in) :: cr
      real(8), intent(in) :: r
      p_abs = sqrt(max(0.0d0, 2.0d0*cr%mu_red*(cr%e_lvl - v_eff(cr, r))))
   end function p_abs

   ! quantize_step(cr, e, node, wts, an) - one action evaluation: turning points
   !   at energy e, inset, 50-node Gauss-Legendre, an = action/(pi*hbar)
   subroutine quantize_step(cr, e, node, wts, an)
      type(carry_t), intent(in) :: cr
      real(8), intent(in) :: e
      real(8), intent(out) :: node(n_quad), wts(n_quad), an
      real(8) :: r1, r2, asum
      integer :: i
      call turn_bounds(cr, e, r1, r2)
      r1 = r1 + eps_turn
      r2 = r2 - eps_turn
      call gl_nodes(r1, r2, node, wts)
      asum = 0.0d0
      do i = 1, n_quad
         asum = asum + wts(i)*sqrt(max(0.0d0, 2.0d0*cr%mu_red*(e - v_eff(cr, node(i)))))
      end do
      an = asum/(pi*hbar_code)
   end subroutine quantize_step

   ! turn_bounds(cr, e, r1, r2) - classical interval at energy e: coarse scan for
   !   the effective-potential minimum, walk out to brackets, bisect both sides
   subroutine turn_bounds(cr, e, r1, r2)
      type(carry_t), intent(in) :: cr
      real(8), intent(in) :: e
      real(8), intent(out) :: r1, r2
      integer, parameter :: n_scan = 1000
      real(8) :: r_well, v_well, r_try, lo, hi
      integer :: i
      if (e >= cr%d_e) then
         call stop_ebkz('ebk_init', 'level energy at or above the dissociation asymptote - '// &
            'no outer turning point exists; level [internal] =', e)
      end if
      r_well = cr%r_eq
      v_well = v_eff(cr, r_well)
      do i = 1, n_scan
         r_try = 0.1d0*cr%r_eq + (5.0d0*cr%r_eq - 0.1d0*cr%r_eq)*dble(i - 1)/dble(n_scan - 1)
         if (v_eff(cr, r_try) < v_well) then
            v_well = v_eff(cr, r_try)
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
         if (v_eff(cr, lo) > e) exit
         if (i == 60) call stop_ebkx('ebk_init', 'no inner turning point (inner wall not rising)')
      end do
      r1 = bisect(cr, lo, r_well, e)
      hi = r_well
      do i = 1, 100
         hi = 1.5d0*hi + 0.1d0
         if (v_eff(cr, hi) > e) exit
         if (i == 100) call stop_ebkx('ebk_init', 'no outer turning point below the asymptote')
      end do
      r2 = bisect(cr, r_well, hi, e)
   end subroutine turn_bounds

   ! bisect - root of V_eff(r) = e on [a, b] (V_eff(a) > e > V_eff(b) or vice
   !          versa; 200 iterations drive the bracket under 1e-12 relative)
   real(8) function bisect(cr, a, b, e)
      type(carry_t), intent(in) :: cr
      real(8), intent(in) :: a, b, e
      real(8) :: lo, hi, mid
      integer :: i
      lo = a
      hi = b
      do i = 1, 200
         mid = 0.5d0*(lo + hi)
         if (v_eff(cr, mid) > e) then
            if (v_eff(cr, lo) > e) then
               lo = mid
            else
               hi = mid
            end if
         else
            if (v_eff(cr, lo) > e) then
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

   ! named-abort helpers (name the caller; STOP 1)
   subroutine stop_ebkx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_ebkx

   subroutine stop_ebkz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_ebkz

   subroutine stop_ebki(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_ebki

   subroutine stop_ebki2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_ebki2

end module samp_ebk
