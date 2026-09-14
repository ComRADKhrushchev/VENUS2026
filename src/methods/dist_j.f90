!=====================================================================
! dist_j.f90 - distribution member (methods/ library): rotational-state
!              sampling - the draw of the rotational quantum number J
!              itself from the thermal Boltzmann-J distribution over a
!              LINEAR rotor's ladder
!              P(J) prop. (2J+1)*exp(-E_J/kT), E_J = hbar^2*J(J+1)/(2*I)
!              with I the perpendicular moment of inertia of the carrier's
!              equilibrium geometry (roster-derived, never a duplicated
!              constant).
! Division of labor (the ruling this member and dist_rotation share): the
!   J-draw and the orientation/angular-velocity realization at a GIVEN J
!   are separate members - THIS member draws J and writes NO q/p; its
!   runtime product is the integer J (read back through j_last_j, the
!   same table-getter posture as the boltzmann quantum-number member).
!   dist_rotation realizes the fragment's orientation and angular
!   momentum at the J it is given. Wiring J from this member's table into
!   rotation's parameter is assembly-phase business, not member code.
! Contract:
!   j_init(p) - validate the member parameter, derive the carrier's
!     perpendicular moment of inertia from the roster equilibrium
!     geometry, precompute the exponent constant, the support bound and
!     the rejection envelope, and bind the served fragment. Loud-stops
!     name the offending field and value:
!       t_rot <= 0 (temperature [K]); no roster fragment table (assembly
!       order); no per-fragment scheme array; length mismatch against the
!       roster; the member word absent from the scheme-code table (a
!       library-assembly failure); zero or multiple carrier fragments
!       (one member instance carries one parameter set); a carrier with
!       fewer than two atoms (a single atom has no rotational degree of
!       freedom); a carrier geometry whose inertia tensor vanishes
!       (coincident atoms); a nonlinear carrier (this member samples the
!       LINEAR rotor ladder - the symmetric-top face with its K structure
!       is not carried: it needs a symmetry-axis designation the roster
!       does not hold); a rotational temperature so high the support scan
!       reaches its cap without the tail cut (the ladder is effectively
!       flat to the cap - the support bound must never silently clamp).
!   j_sample(sta) - one draw by rejection: J uniform on [0, j_max]
!     (integer), accepted when u*p_env <= P(J); the envelope is the
!     distribution maximum (exact - P is unimodal with the peak located
!     by the scan). The draw lands in the member's quantum-number slot
!     (j_last_j). sta is the sampler runtime channel; this member writes
!     NO q/p components.
!   j_last_j() - copy of the last drawn rotational quantum number.
!   j_c_rot() - the exponent constant hbar^2/(2*I*kb*T) [-] (exposed for
!     criteria anchoring of the roster-derived moment of inertia).
!   j_j_max() - the support bound of the rejection table [-].
! Variant inventory against the old defect register (semantics only, no
!   citations): the thermal symmetric-top face and the equipartition face
!   of the old rotational selector are NOT carried - both duplicate the
!   J-draw business this member owns under the division ruling, the first
!   needing symmetry-axis designation the roster lacks; the nuclear-spin
!   statistical weights are NOT carried (system-specific spin data - a
!   knowledge-ownership exclusion, not a member face).
! Units: t_rot [K]; I [amu*Ang^2] (roster); hbar_code, kb_code internal;
!   the exponent constant and J are dimensionless. E_J = hbar^2*J(J+1)/
!   (2*I) in internal energy units matches kb_code*T by construction.
!=====================================================================
module dist_j
   use consts,        only: hbar_code, kb_code
   use rng,           only: rng_u
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use sampler,       only: samp_code
   use linalg,        only: eig_sym
   implicit none
   private
   public :: j_params_t, j_init, j_sample, j_last_j, j_c_rot, j_j_max

   integer, parameter :: max_reject = 1000000   ! rejection-sampling draw cap [-]
   integer, parameter :: max_support = 10000000 ! support-scan cap [-]
   real(8), parameter :: tail_cut = 1.0d-10     ! support tail cutoff vs the peak [-]
   real(8), parameter :: lin_tol = 1.0d-8       ! linearity eigenvalue ratio tolerance [-]

   ! Member-owned parameters (the moment of inertia is roster-derived)
   type :: j_params_t
      real(8) :: t_rot = 0.0d0     ! rotational temperature [K]
   end type

   ! Derived member state (module-private; rebuilt by every j_init)
   real(8) :: c_rot = 0.0d0        ! exponent constant hbar^2/(2*I*kb*T) [-]
   real(8) :: p_env = 0.0d0        ! rejection envelope = max P(J) [-]
   integer :: j_max = 0            ! support bound [quantum number]
   integer :: last_j = 0           ! last drawn rotational quantum number [-]
   logical :: inited = .false.     ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! j_init(p) - validate the temperature, derive the linear-rotor moment
   !             of inertia from the carrier geometry, precompute the
   !             distribution table, bind the carrier
   subroutine j_init(p)
      type(j_params_t), intent(in) :: p ! member parameters (rotational temperature)
      real(8) :: com(3), r(3), dv(3)
      real(8) :: tens(3,3), ev(3)
      real(8) :: i_perp, wt, pj, p_max
      integer :: i, j, k, code, n_carry, i_carry, j_peak
      logical :: cut

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%t_rot <= 0.0d0) then
         call stop_jz('j_init', 'rotational temperature t_rot [K] must be positive, got', p%t_rot)
      end if

      ! 2. fragment binding (assembly-order and carrier guards first)
      if (.not. allocated(roster%frag)) then
         call stop_jx('j_init', 'no roster fragment table - the roster must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_jx('j_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_ji2('j_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('j')
      if (code == 0) then
         call stop_jx('j_init', 'member word "j" is not in the scheme-code table '// &
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
         call stop_jx('j_init', 'no roster fragment carries the j scheme - this member '// &
            'would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_ji2('j_init', 'two or more fragments carry the j scheme, but one member '// &
            'instance carries one parameter set; count =', n_carry, 'first carrier =', i_carry)
      end if
      associate (fr => roster%frag(i_carry))
         if (fr%nat < 2) then
            call stop_ji('j_init', 'the j carrier must carry at least two atoms (a single '// &
               'atom has no rotational degree of freedom); atom count =', fr%nat)
         end if

         ! 3. linear-rotor moment of inertia from the equilibrium geometry:
         !    COM frame, inertia tensor, eigenvalue test for linearity
         wt = sum(roster%mass(fr%list))
         com = 0.0d0
         do i = 1, fr%nat
            com = com + roster%mass(fr%list(i))*fr%qz(3*i-2:3*i)
         end do
         com = com/wt
         tens = 0.0d0
         do i = 1, fr%nat
            r = fr%qz(3*i-2:3*i) - com
            do j = 1, 3
               do k = 1, 3
                  tens(j,k) = tens(j,k) + roster%mass(fr%list(i))* &
                     ((dot_product(r,r))*merge(1,0,j == k) - r(j)*r(k))
               end do
            end do
         end do
         i_perp = tens(1,1) + tens(2,2) + tens(3,3)   ! trace
         if (i_perp <= 0.0d0) then
            call stop_jz('j_init', 'degenerate carrier geometry: moment of inertia [amu*Ang^2] =', &
               i_perp)
         end if
         if (fr%nat == 2) then
            i_perp = i_perp/2.0d0        ! two points are always collinear: trace = 2*I_perp
         else
            call eig_sym(tens, ev)
            if (ev(1) > lin_tol*ev(3)) then
               call stop_jz('j_init', 'the j carrier must be a linear fragment (this member '// &
                  'samples the linear-rotor ladder; smallest-to-largest inertia eigenvalue '// &
                  'ratio =', ev(1)/ev(3))
            end if
            i_perp = 0.5d0*(ev(2) + ev(3))
         end if
         if (i_perp <= 0.0d0) then
            call stop_jz('j_init', 'degenerate carrier geometry: moment of inertia [amu*Ang^2] =', &
               i_perp)
         end if
      end associate

      ! 4. exponent constant, support scan, envelope
      c_rot = hbar_code*hbar_code/(2.0d0*i_perp*kb_code*p%t_rot)
      p_max = 0.0d0
      j_peak = 0
      j_max = max_support
      cut = .false.
      do j = 0, max_support
         pj = dble(2*j + 1)*exp(-c_rot*dble(j*(j + 1)))
         if (pj > p_max) then
            p_max = pj
            j_peak = j
         end if
         if (j > j_peak .and. pj < tail_cut*p_max) then
            j_max = j
            cut = .true.
            exit
         end if
      end do
      if (.not. cut) then
         call stop_jz('j_init', 'support scan reached the cap 10000000 without the tail cut - '// &
            'the Boltzmann-J ladder is effectively flat to the cap; rotational temperature '// &
            't_rot [K] =', p%t_rot)
      end if
      p_env = p_max
      last_j = 0
      inited = .true.
   end subroutine j_init

   !------------------------------------------------------------------
   ! j_sample(sta) - one rejection draw of J; sta passes through untouched
   subroutine j_sample(sta)
      type(state_t), intent(inout) :: sta  ! sampler runtime channel (no q/p components owned)
      integer :: jtry, n_try
      real(8) :: pj
      if (.not. inited) then
         call stop_jx('j_sample', 'sample called before j_init - the member carries no '// &
            'parameter set')
      end if
      do n_try = 1, max_reject
         jtry = int(rng_u()*dble(j_max + 1))
         if (jtry > j_max) jtry = j_max          ! defensive clamp (rng_u < 1 keeps this idle)
         pj = dble(2*jtry + 1)*exp(-c_rot*dble(jtry*(jtry + 1)))
         if (rng_u()*p_env <= pj) exit
         if (n_try == max_reject) then
            call stop_jx('j_sample', 'rejection sampling exceeded its draw cap (the J table '// &
               'is inconsistent with the envelope)')
         end if
      end do
      last_j = jtry
   end subroutine j_sample

   !------------------------------------------------------------------
   ! getters (copies of the derived member state)
   !------------------------------------------------------------------
   pure integer function j_last_j()
      j_last_j = last_j
   end function j_last_j

   pure real(8) function j_c_rot()
      j_c_rot = c_rot
   end function j_c_rot

   pure integer function j_j_max()
      j_j_max = j_max
   end function j_j_max

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_jx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_jx

   subroutine stop_jz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_jz

   subroutine stop_ji(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_ji

   subroutine stop_ji2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_ji2

end module dist_j
