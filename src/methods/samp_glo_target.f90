!=====================================================================
! samp_glo_target.f90 - distribution member: ghost Langevin-oscillator
!   target selection, bath-coupling loading (one ghost per carrier
!   atom), and the stationary ghost draw
! Design:
!   Loads strengths by the fluctuation-dissipation algebra
!   (sigma_noise^2 = 2*gamma*m*kb*t/dt; sigma_p^2 = m*kb*t;
!   sigma_q^2 = kb*t/(m*w^2)) and draws the exact stationary ghost
!   state (member product via getters; Langevin is dynamics-side).
!   Units: t_glo [K]; gamma_glo, w_ghost, dt as declared; sigma_p
!   [amu*Ang/(10 fs)]; sigma_q [Ang]; sigma_noise [amu*Ang/(10 fs)^2].
!=====================================================================
module samp_glo_target
   use consts,        only: kb_code
   use rng,           only: rng_gauss
   use state,         only: state_t
   use config,        only: reactants, prop_cfg => propagator   ! alias to avoid clashing with the module name
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code
   implicit none
   private
   public :: glo_target_params_t, glo_target_init, glo_target_draw, glo_target_realize, glo_target_sample, &
             glo_target_t_glo, glo_target_gamma, glo_target_w, glo_target_dt_loaded, &
             glo_target_frag, glo_target_nat, glo_target_ghost_size, &
             glo_target_last_q, glo_target_last_p, &
             glo_target_sigma_noise, glo_target_sigma_p, glo_target_sigma_q

   ! Member-owned parameters
   type :: glo_target_params_t
      real(8) :: t_glo = 0.0d0     ! bath temperature [K]
      real(8) :: gamma_glo = 0.0d0 ! damping coefficient [1/(10 fs)]
      real(8) :: w_ghost = 0.0d0   ! ghost oscillator frequency [1/(10 fs)]
   end type

   ! Derived member state (module-private; rebuilt by every glo_target_init) -
   ! one row per carrier; the bath scalars are shared (equal treatment)
   type :: carry_t
      integer :: frag = 0             ! carrier list_atoms fragment row [-]
      integer :: nat = 0              ! carrier atom count [-]
      real(8), allocatable :: sig_noise(:) ! per-atom noise amplitudes [amu*Ang/(10 fs)^2]
      real(8), allocatable :: sig_p(:)     ! per-atom stationary momentum widths [amu*Ang/(10 fs)]
      real(8), allocatable :: sig_q(:)     ! per-atom stationary coordinate widths [Ang]
      real(8), allocatable :: ghost_q(:)   ! last drawn ghost coordinates [Ang] (3 per target atom)
      real(8), allocatable :: ghost_p(:)   ! last drawn ghost momenta [amu*Ang/(10 fs)]
   end type carry_t
   type(carry_t), allocatable :: carry(:)  ! carrier rows (list_atoms order)
   real(8) :: t_glo_c = 0.0d0        ! bath temperature [K]
   real(8) :: gamma_c = 0.0d0        ! damping coefficient [1/(10 fs)]
   real(8) :: w_c = 0.0d0            ! ghost frequency [1/(10 fs)]
   real(8) :: dt_c = 0.0d0           ! the propagator step read at init [10 fs]
   logical :: inited = .false.       ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! glo_target_init(p) - validate the bath parameters and the propagator
   !                      step, bind the carrier, load the coupling strengths
   subroutine glo_target_init(p)
      type(glo_target_params_t), intent(in) :: p ! member parameters (bath coupling)
      integer :: i, j, kc, code, n_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%t_glo <= 0.0d0) then
         call stop_gloz('glo_target_init', 'bath temperature t_glo [K] must be positive, got', p%t_glo)
      end if
      if (p%gamma_glo <= 0.0d0) then
         call stop_gloz('glo_target_init', 'damping coefficient gamma_glo [1/(10 fs)] must be '// &
            'positive, got', p%gamma_glo)
      end if
      if (p%w_ghost <= 0.0d0) then
         call stop_gloz('glo_target_init', 'ghost oscillator frequency w_ghost [1/(10 fs)] must '// &
            'be positive, got', p%w_ghost)
      end if
      if (prop_cfg%dt <= 0.0d0) then
         call stop_gloz('glo_target_init', 'the propagator step dt [10 fs] the bath strengths are '// &
            'calibrated against must be positive, got', prop_cfg%dt)
      end if

      ! 2. fragment binding (assembly-order and carrier guards)
      if (.not. allocated(list_atoms%frag)) then
         call stop_glox('glo_target_init', 'no list_atoms fragment table - the atom list must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_glox('glo_target_init', 'no per-fragment scheme array (reactants%'// &
            'dist_scheme unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_gloi2('glo_target_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('glo_target')
      if (code == 0) then
         call stop_glox('glo_target_init', 'member word "glo_target" is not in the '// &
            'scheme-code table (incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_glox('glo_target_init', 'no list_atoms fragment carries the glo_target '// &
            'scheme - this member would never be dispatched')
      end if

      ! 3. carrier loop (list_atoms order): coupling-strength loading (the
      !    conservation algebra: sigma_noise^2 = 2*gamma*m*kb*t/dt;
      !    sigma_p^2 = m*kb*t; sigma_q^2 = kb*t/(m*w^2), per carrier atom)
      !    and the zeroed ghost state
      if (allocated(carry)) deallocate (carry)
      allocate (carry(n_carry))
      kc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         kc = kc + 1
         carry(kc)%frag = i
         carry(kc)%nat = list_atoms%frag(i)%nat
         allocate (carry(kc)%sig_noise(carry(kc)%nat), carry(kc)%sig_p(carry(kc)%nat), &
                   carry(kc)%sig_q(carry(kc)%nat))
         do j = 1, carry(kc)%nat
            associate (m_i => list_atoms%mass(list_atoms%frag(carry(kc)%frag)%list(j)))
               carry(kc)%sig_noise(j) = sqrt(2.0d0*p%gamma_glo*m_i*kb_code*p%t_glo/prop_cfg%dt)
               carry(kc)%sig_p(j) = sqrt(m_i*kb_code*p%t_glo)
               carry(kc)%sig_q(j) = sqrt(kb_code*p%t_glo/(m_i*p%w_ghost**2))
            end associate
         end do
         allocate (carry(kc)%ghost_q(3*carry(kc)%nat), carry(kc)%ghost_p(3*carry(kc)%nat))
         carry(kc)%ghost_q = 0.0d0
         carry(kc)%ghost_p = 0.0d0
      end do
      t_glo_c = p%t_glo
      gamma_c = p%gamma_glo
      w_c = p%w_ghost
      dt_c = prop_cfg%dt
      inited = .true.
   end subroutine glo_target_init

   !------------------------------------------------------------------
   ! glo_target_sample(sta) - draw the ghost stationary state per carrier
   !                          (list_atoms order; one fresh Gaussian per ghost
   !                          component); sta untouched
   subroutine glo_target_draw()
      integer :: i, c, kc
      if (.not. inited) then
         call stop_glox('glo_target_draw', 'draw called before glo_target_init - the '// &
            'member carries no parameter set')
      end if
      do kc = 1, size(carry)
         do i = 1, carry(kc)%nat
            do c = 1, 3
               carry(kc)%ghost_q(3*(i - 1) + c) = carry(kc)%sig_q(i)*rng_gauss()
               carry(kc)%ghost_p(3*(i - 1) + c) = carry(kc)%sig_p(i)*rng_gauss()
            end do
         end do
      end do
   end subroutine glo_target_draw

   !------------------------------------------------------------------
   ! glo_target_realize(sta) - the ghost state is a member product (getters),
   !                            not a state write; nothing to realize
   !------------------------------------------------------------------
   subroutine glo_target_realize(sta)
      type(state_t), intent(inout) :: sta ! sampler runtime channel (nothing written)
      if (.not. inited) then
         call stop_glox('glo_target_realize', 'realize called before glo_target_init - the '// &
            'member carries no parameter set')
      end if
   end subroutine glo_target_realize

   !------------------------------------------------------------------
   ! glo_target_sample(sta) - fused convenience: draw then realize (the
   !                           check-program entry)
   !------------------------------------------------------------------
   subroutine glo_target_sample(sta)
      type(state_t), intent(inout) :: sta
      call glo_target_draw()
      call glo_target_realize(sta)
   end subroutine glo_target_sample

   !------------------------------------------------------------------
   ! getters (copies of the derived member state; the carrier index is
   ! optional and defaults to the first carrier - the single-carrier face)
   !------------------------------------------------------------------
   pure real(8) function glo_target_t_glo()
      glo_target_t_glo = t_glo_c
   end function glo_target_t_glo

   pure real(8) function glo_target_gamma()
      glo_target_gamma = gamma_c
   end function glo_target_gamma

   pure real(8) function glo_target_w()
      glo_target_w = w_c
   end function glo_target_w

   pure real(8) function glo_target_dt_loaded()
      glo_target_dt_loaded = dt_c
   end function glo_target_dt_loaded

   pure integer function glo_target_frag(k)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      glo_target_frag = 0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) glo_target_frag = carry(kc)%frag
   end function glo_target_frag

   pure integer function glo_target_nat(k)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      glo_target_nat = 0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) glo_target_nat = carry(kc)%nat
   end function glo_target_nat

   pure integer function glo_target_ghost_size(k)
      integer, intent(in), optional :: k    ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(k)) kc = k
      glo_target_ghost_size = 0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) glo_target_ghost_size = 3*carry(kc)%nat
   end function glo_target_ghost_size

   pure real(8) function glo_target_last_q(k, kc_in)
      integer, intent(in) :: k               ! ghost component index (1..3*nat) [-]
      integer, intent(in), optional :: kc_in ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(kc_in)) kc = kc_in
      glo_target_last_q = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) glo_target_last_q = carry(kc)%ghost_q(k)
   end function glo_target_last_q

   pure real(8) function glo_target_last_p(k, kc_in)
      integer, intent(in) :: k               ! ghost component index (1..3*nat) [-]
      integer, intent(in), optional :: kc_in ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(kc_in)) kc = kc_in
      glo_target_last_p = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) glo_target_last_p = carry(kc)%ghost_p(k)
   end function glo_target_last_p

   pure real(8) function glo_target_sigma_noise(i, kc_in)
      integer, intent(in) :: i                ! carrier-atom index [-]
      integer, intent(in), optional :: kc_in  ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(kc_in)) kc = kc_in
      glo_target_sigma_noise = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) glo_target_sigma_noise = carry(kc)%sig_noise(i)
   end function glo_target_sigma_noise

   pure real(8) function glo_target_sigma_p(i, kc_in)
      integer, intent(in) :: i                ! carrier-atom index [-]
      integer, intent(in), optional :: kc_in  ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(kc_in)) kc = kc_in
      glo_target_sigma_p = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) glo_target_sigma_p = carry(kc)%sig_p(i)
   end function glo_target_sigma_p

   pure real(8) function glo_target_sigma_q(i, kc_in)
      integer, intent(in) :: i                ! carrier-atom index [-]
      integer, intent(in), optional :: kc_in  ! carrier row (default 1)
      integer :: kc
      kc = 1
      if (present(kc_in)) kc = kc_in
      glo_target_sigma_q = 0.0d0
      if (inited .and. kc >= 1 .and. kc <= size(carry)) glo_target_sigma_q = carry(kc)%sig_q(i)
   end function glo_target_sigma_q

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; STOP 1)
   subroutine stop_glox(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_glox

   subroutine stop_gloz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_gloz

   subroutine stop_gloi2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_gloi2

end module samp_glo_target
