!=====================================================================
! samp_thermalize.f90 - distribution member: MD equilibration mechanism
!   - Maxwell seed of the owned components, exact drift removal, then a
!   bounded run of n_eq propagation steps through the evolution service
! Design:
!   The COM-velocity projection p'_i = p_i - (m_i/M)*P removes exactly
!   three degrees of freedom per carrier, so E[KE] = sum over carriers of
!   (3*nat-3)*kb*t_eq/2 exactly. One parameter set serves every carrier
!   equally: each carrier is seeded and drift-removed in list_atoms order,
!   then ONE shared equilibration run propagates the whole state (no
!   bystander freezing - co-selected carriers ride the same bath run; no
!   in-run bath coupling, that is the thermostat's face). The member
!   writes its carriers' owned momenta, swaps the global step to dt_eq
!   for the run and restores step + time after.
!   Units: t_eq [K]; dt_eq [10 fs]; p [amu*Ang/(10 fs)]; E [internal].
!=====================================================================
module samp_thermalize
   use consts,        only: kb_code
   use rng,           only: rng_gauss
   use state,         only: state_t
   use config,        only: reactants, prop_cfg => propagator   ! alias to avoid clashing with the module name
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code
   use propagator,    only: prop_step     ! evolution service (an internal call is legal)
   implicit none
   private
   public :: thermalize_params_t, thermalize_init, thermalize_draw, thermalize_realize, thermalize_sample, &
             thermalize_t_eq, thermalize_n_eq, thermalize_dt_eq, thermalize_e_kin_expect, &
             thermalize_frag, thermalize_nat

   ! Member-owned parameters
   type :: thermalize_params_t
      real(8) :: t_eq = 0.0d0    ! target temperature [K] (every carrier)
      integer :: n_eq = 0        ! equilibration step count [-] (one shared run)
      real(8) :: dt_eq = 0.0d0   ! equilibration step size [10 fs]
   end type

   ! Derived member state (module-private; rebuilt by every thermalize_init) -
   ! one row per carrier plus the shared scalar copies
   type :: carry_t
      integer :: frag = 0             ! carrier list_atoms fragment row [-]
      integer :: nat = 0              ! carrier atom count [-]
      real(8), allocatable :: sigma(:) ! per-atom seed widths sqrt(m*kb*t_eq) [amu*Ang/(10 fs)]
      real(8), allocatable :: p_seed(:) ! the drawn drift-removed seed momenta (3*nat record) [amu*Ang/(10 fs)]
   end type carry_t
   type(carry_t), allocatable :: carry(:)  ! carrier rows (list_atoms order)
   real(8) :: t_eq_c = 0.0d0        ! target temperature [K]
   integer :: n_eq_c = 0            ! step count [-]
   real(8) :: dt_eq_c = 0.0d0       ! step size [10 fs]
   real(8) :: e_expect_c = 0.0d0    ! declared KE expectation, all carriers [internal energy]
   logical :: inited = .false.      ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! thermalize_init(p) - validate the equilibration parameters, bind every
   !                      carrier, precompute the expectation constant
   subroutine thermalize_init(p)
      type(thermalize_params_t), intent(in) :: p ! member parameters (equilibration)
      integer :: i, k, code, n_carry

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

      ! 2. carrier binding (assembly-order guards; the atom list fragments whose
      !    dist_scheme rows select this member, list_atoms order - one parameter
      !    set serves them all equally)
      if (.not. allocated(list_atoms%frag)) then
         call stop_equix('thermalize_init', 'no list_atoms fragment table - the atom list must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_equix('thermalize_init', 'no per-fragment scheme array (reactants%'// &
            'dist_scheme unallocated) - member init cannot bind its carriers')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_equii2('thermalize_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('thermalize')
      if (code == 0) then
         call stop_equix('thermalize_init', 'member word "thermalize" is not in the '// &
            'scheme-code table (incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_equix('thermalize_init', 'no list_atoms fragment carries the thermalize '// &
            'scheme - this member would never be dispatched')
      end if
      if (allocated(carry)) deallocate (carry)
      allocate (carry(n_carry))
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code) then
            k = k + 1
            carry(k)%frag = i
         end if
      end do

      ! 3. seed widths from the atom list masses and the expectation constant
      !    (sigma^2 = m*kb*t_eq per carrier atom; E[KE] = sum over carriers of
      !    (3*nat-3)*kb*t_eq/2)
      e_expect_c = 0.0d0
      do k = 1, n_carry
         carry(k)%nat = list_atoms%frag(carry(k)%frag)%nat
         allocate (carry(k)%sigma(carry(k)%nat), carry(k)%p_seed(3*carry(k)%nat))
         do i = 1, carry(k)%nat
            carry(k)%sigma(i) = sqrt(list_atoms%mass(list_atoms%frag(carry(k)%frag)%list(i))*kb_code*p%t_eq)
         end do
         e_expect_c = e_expect_c + (3.0d0*dble(carry(k)%nat) - 3.0d0)*kb_code*p%t_eq/2.0d0
      end do
      t_eq_c = p%t_eq
      n_eq_c = p%n_eq
      dt_eq_c = p%dt_eq
      inited = .true.
   end subroutine thermalize_init

   !------------------------------------------------------------------
   ! thermalize_sample(sta) - Maxwell seed + exact drift removal of every
   !                          carrier (list_atoms order), then the ONE bounded
   !                          equilibration run through the evolution service
   !------------------------------------------------------------------
   ! thermalize_draw() - per carrier (list_atoms order, the historical
   !                     consumption order: carrier, atom, component): the
   !                     Maxwell seed into the record, then the exact drift
   !                     removal (the mass-weighted COM-velocity projection
   !                     p'_i = p_i - (m_i/M)*P - the owned-set total
   !                     momentum becomes exactly zero; the Gaussian projected
   !                     this way loses exactly three degrees of freedom - the
   !                     declared expectation is mass-independent); no state
   !                     access
   !------------------------------------------------------------------
   subroutine thermalize_draw()
      integer :: i, c, k
      real(8) :: wt, p_tot(3)
      if (.not. inited) then
         call stop_equix('thermalize_draw', 'draw called before thermalize_init - the member '// &
            'carries no parameter set')
      end if
      do k = 1, size(carry)
         do i = 1, carry(k)%nat
            do c = 1, 3
               carry(k)%p_seed(3*(i - 1) + c) = carry(k)%sigma(i)*rng_gauss()
            end do
         end do
         wt = sum(list_atoms%mass(list_atoms%frag(carry(k)%frag)%list))
         p_tot = 0.0d0
         do i = 1, carry(k)%nat
            p_tot = p_tot + carry(k)%p_seed(3*i - 2:3*i)
         end do
         do i = 1, carry(k)%nat
            carry(k)%p_seed(3*i - 2:3*i) = carry(k)%p_seed(3*i - 2:3*i) &
                                           - (list_atoms%mass(list_atoms%frag(carry(k)%frag)%list(i))/wt)*p_tot
         end do
      end do
   end subroutine thermalize_draw

   !------------------------------------------------------------------
   ! thermalize_realize(sta) - write the drawn seed momenta into the owned
   !                           components, then the ONE bounded equilibration
   !                           run at the member step size through the
   !                           evolution service (every carrier rides it): the
   !                           global step and the trajectory time are swapped
   !                           for the run and RESTORED after it (the run
   !                           propagates under the real assembled force law;
   !                           any zero-force criteria trajectories come from
   !                           their check-side container binding, not from
   !                           this member) - no RNG draws
   !------------------------------------------------------------------
   subroutine thermalize_realize(sta)
      type(state_t), intent(inout) :: sta ! sampler runtime channel (owned p components written)
      integer :: i, k, g_atom, i_step, last_own
      real(8) :: t_saved, dt_saved
      if (.not. inited) then
         call stop_equix('thermalize_realize', 'realize called before thermalize_init - the '// &
            'member carries no parameter set')
      end if
      last_own = 0
      do k = 1, size(carry)
         last_own = max(last_own, 3*list_atoms%frag(carry(k)%frag)%list(carry(k)%nat))
      end do
      if (size(sta%p) < last_own) then
         call stop_equii2('thermalize_realize', 'the state momentum vector ends before a '// &
            "carrier's last owned component: state dof =", size(sta%p), &
            'last owned component =', last_own)
      end if
      do k = 1, size(carry)
         do i = 1, carry(k)%nat
            g_atom = list_atoms%frag(carry(k)%frag)%list(i)
            sta%p(3*(g_atom - 1) + 1:3*g_atom) = carry(k)%p_seed(3*i - 2:3*i)
         end do
      end do
      t_saved = sta%t
      dt_saved = prop_cfg%dt
      prop_cfg%dt = dt_eq_c
      do i_step = 1, n_eq_c
         call prop_step(sta)
      end do
      prop_cfg%dt = dt_saved
      sta%t = t_saved
   end subroutine thermalize_realize

   !------------------------------------------------------------------
   ! thermalize_sample(sta) - fused convenience: draw then realize (the
   !                          check-program entry)
   !------------------------------------------------------------------
   subroutine thermalize_sample(sta)
      type(state_t), intent(inout) :: sta
      call thermalize_draw()
      call thermalize_realize(sta)
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
      thermalize_e_kin_expect = e_expect_c      ! the sum over all carriers
   end function thermalize_e_kin_expect

   pure integer function thermalize_frag(k)
      integer, intent(in), optional :: k     ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      thermalize_frag = 0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) thermalize_frag = carry(kc)%frag
   end function thermalize_frag

   pure integer function thermalize_nat(k)
      integer, intent(in), optional :: k     ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      thermalize_nat = 0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) thermalize_nat = carry(kc)%nat
   end function thermalize_nat

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; STOP 1)
   subroutine stop_equix(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_equix

   subroutine stop_equiz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_equiz

   subroutine stop_equii(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_equii

   subroutine stop_equii2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_equii2

end module samp_thermalize
