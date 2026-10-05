!=====================================================================
! samp_j.f90 - distribution member: draw of the rotational quantum
!   number J from the thermal linear-rotor ladder
!   P(J) prop. (2J+1)*exp(-E_J/kT), E_J = hbar^2*J(J+1)/(2*I)
! Design:
!   One temperature serves every carrier equally; the moment of inertia
!   is derived per carrier from its list_atoms equilibrium geometry (never
!   duplicated), so the ladder constants, support and envelope are
!   per-carrier state. Each carrier must be linear. This member draws J
!   only (read back via j_last_j) and writes no q/p; dist_rotation
!   realizes the orientation at a given J, and wiring J between them is
!   assembly-phase business.
!   Units: t_rot [K]; I [amu*Ang^2]; hbar_code/kb_code internal.
!=====================================================================
module samp_j
   use consts,        only: hbar_code, kb_code
   use rng,           only: rng_u
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code
   use linalg,        only: eig_sym
   implicit none
   private
   public :: j_params_t, j_init, j_draw, j_realize, j_sample, j_last_j, j_c_rot, j_j_max

   integer, parameter :: max_reject = 1000000   ! rejection-sampling draw cap [-]
   integer, parameter :: max_support = 10000000 ! support-scan cap [-]
   real(8), parameter :: tail_cut = 1.0d-10     ! support tail cutoff vs the peak [-]
   real(8), parameter :: lin_tol = 1.0d-8       ! linearity eigenvalue ratio tolerance [-]

   ! Member-owned parameters (the moments of inertia are list_atoms-derived)
   type :: j_params_t
      real(8) :: t_rot = 0.0d0     ! rotational temperature [K] (every carrier)
   end type

   ! Derived member state (module-private; rebuilt by every j_init) - one row
   ! per carrier: the ladder is set by the carrier's own moment of inertia
   type :: carry_t
      integer :: frag = 0          ! carrier list_atoms fragment row [-]
      real(8) :: c_rot = 0.0d0     ! exponent constant hbar^2/(2*I*kb*T) [-]
      real(8) :: p_env = 0.0d0     ! rejection envelope = max P(J) [-]
      integer :: j_max = 0         ! support bound [quantum number]
      integer :: last_j = 0        ! last drawn rotational quantum number [-]
   end type carry_t
   type(carry_t), allocatable :: carry(:)  ! carrier rows (list_atoms order)
   logical :: inited = .false.     ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! j_init(p) - validate the temperature, derive the linear-rotor moment
   !             of inertia of every carrier, precompute the per-carrier
   !             distribution tables, bind the carriers
   subroutine j_init(p)
      type(j_params_t), intent(in) :: p ! member parameters (rotational temperature)
      real(8) :: com(3), r(3)
      real(8) :: tens(3,3), ev(3)
      real(8) :: i_perp, wt, pj, p_max
      integer :: i, j, k, m, kc, code, n_carry, j_peak
      logical :: cut

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%t_rot <= 0.0d0) then
         call stop_jz('j_init', 'rotational temperature t_rot [K] must be positive, got', p%t_rot)
      end if

      ! 2. carrier binding (assembly-order guards; the atom list fragments whose
      !    dist_scheme rows select this member, list_atoms order)
      if (.not. allocated(list_atoms%frag)) then
         call stop_jx('j_init', 'no list_atoms fragment table - the atom list must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_jx('j_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its carriers')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_ji2('j_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('j')
      if (code == 0) then
         call stop_jx('j_init', 'member word "j" is not in the scheme-code table '// &
            '(incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_jx('j_init', 'no list_atoms fragment carries the j scheme - this member '// &
            'would never be dispatched')
      end if
      if (allocated(carry)) deallocate (carry)
      allocate (carry(n_carry))
      kc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         kc = kc + 1
         carry(kc)%frag = i
         associate (fr => list_atoms%frag(i), cr => carry(kc))
            if (fr%nat < 2) then
               call stop_ji('j_init', 'the j carrier must carry at least two atoms (a single '// &
                  'atom has no rotational degree of freedom); atom count =', fr%nat)
            end if

            ! 3. linear-rotor moment of inertia from the equilibrium geometry:
            !    COM frame, inertia tensor, eigenvalue test for linearity
            wt = sum(list_atoms%mass(fr%list))
            com = 0.0d0
            do j = 1, fr%nat
               com = com + list_atoms%mass(fr%list(j))*fr%qz(3*j-2:3*j)
            end do
            com = com/wt
            tens = 0.0d0
            do j = 1, fr%nat
               r = fr%qz(3*j-2:3*j) - com
               do k = 1, 3
                  do m = 1, 3
                     tens(k,m) = tens(k,m) + list_atoms%mass(fr%list(j))* &
                        ((dot_product(r,r))*merge(1,0,k == m) - r(k)*r(m))
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

            ! 4. exponent constant, support scan, envelope (per carrier)
            cr%c_rot = hbar_code*hbar_code/(2.0d0*i_perp*kb_code*p%t_rot)
            p_max = 0.0d0
            j_peak = 0
            cr%j_max = max_support
            cut = .false.
            do j = 0, max_support
               pj = dble(2*j + 1)*exp(-cr%c_rot*dble(j*(j + 1)))
               if (pj > p_max) then
                  p_max = pj
                  j_peak = j
               end if
               if (j > j_peak .and. pj < tail_cut*p_max) then
                  cr%j_max = j
                  cut = .true.
                  exit
               end if
            end do
            if (.not. cut) then
               call stop_jz('j_init', 'support scan reached the cap 10000000 without the tail cut - '// &
                  'the Boltzmann-J ladder is effectively flat to the cap; rotational temperature '// &
                  't_rot [K] =', p%t_rot)
            end if
            cr%p_env = p_max
            cr%last_j = 0
         end associate
      end do
      inited = .true.
   end subroutine j_init

   !------------------------------------------------------------------
   ! j_sample(sta) - one rejection draw of J per carrier (list_atoms order);
   !                  sta passes through untouched
   subroutine j_draw()
      integer :: kc, jtry, n_try
      real(8) :: pj
      if (.not. inited) then
         call stop_jx('j_draw', 'draw called before j_init - the member carries no '// &
            'parameter set')
      end if
      do kc = 1, size(carry)
         do n_try = 1, max_reject
            jtry = int(rng_u()*dble(carry(kc)%j_max + 1))
            if (jtry > carry(kc)%j_max) jtry = carry(kc)%j_max   ! defensive clamp (rng_u < 1 keeps this idle)
            pj = dble(2*jtry + 1)*exp(-carry(kc)%c_rot*dble(jtry*(jtry + 1)))
            if (rng_u()*carry(kc)%p_env <= pj) exit
            if (n_try == max_reject) then
               call stop_jx('j_draw', 'rejection sampling exceeded its draw cap (the J table '// &
                  'is inconsistent with the envelope)')
            end if
         end do
         carry(kc)%last_j = jtry
      end do
   end subroutine j_draw

   !------------------------------------------------------------------
   ! j_realize(sta) - the registered gap: the orientation/momenta realization
   !                   of the drawn J (the rotation member's business) has no
   !                   wiring yet - the member's product stays last_j
   !------------------------------------------------------------------
   subroutine j_realize(sta)
      type(state_t), intent(inout) :: sta  ! sampler runtime channel (nothing written)
      if (.not. inited) then
         call stop_jx('j_realize', 'realize called before j_init - the member carries no '// &
            'parameter set')
      end if
   end subroutine j_realize

   !------------------------------------------------------------------
   ! j_sample(sta) - fused convenience: draw then realize (the check-program
   !                  entry; production dispatch calls the two halves)
   !------------------------------------------------------------------
   subroutine j_sample(sta)
      type(state_t), intent(inout) :: sta
      call j_draw()
      call j_realize(sta)
   end subroutine j_sample

   !------------------------------------------------------------------
   ! getters (copies of the per-carrier derived state; k = carrier row)
   !------------------------------------------------------------------
   pure integer function j_last_j(k)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      j_last_j = 0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) j_last_j = carry(kc)%last_j
   end function j_last_j

   pure real(8) function j_c_rot(k)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      j_c_rot = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) j_c_rot = carry(kc)%c_rot
   end function j_c_rot

   pure integer function j_j_max(k)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      j_j_max = 0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) j_j_max = carry(kc)%j_max
   end function j_j_max

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; STOP 1)
   subroutine stop_jx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_jx

   subroutine stop_jz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_jz

   subroutine stop_ji(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_ji

   subroutine stop_ji2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_ji2

end module samp_j
