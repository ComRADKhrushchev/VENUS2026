!=====================================================================
! dist_glo_target.f90 - distribution member: ghost Langevin-oscillator
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
module dist_glo_target
   use consts,        only: kb_code
   use rng,           only: rng_gauss
   use state,         only: state_t
   use config,        only: reactants, prop_cfg => propagator   ! alias to avoid clashing with the module name
   use config_roster, only: roster
   use sampler,       only: samp_code
   implicit none
   private
   public :: glo_target_params_t, glo_target_init, glo_target_sample, &
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

   ! Derived member state (module-private; rebuilt by every glo_target_init)
   real(8) :: t_glo_c = 0.0d0        ! bath temperature [K]
   real(8) :: gamma_c = 0.0d0        ! damping coefficient [1/(10 fs)]
   real(8) :: w_c = 0.0d0            ! ghost frequency [1/(10 fs)]
   real(8) :: dt_c = 0.0d0           ! the propagator step read at init [10 fs]
   integer :: frag_b = 0             ! bound carrier fragment index [-]
   integer :: nat_b = 0              ! carrier atom count [-]
   real(8), allocatable :: sig_noise_c(:) ! per-atom noise amplitudes [amu*Ang/(10 fs)^2]
   real(8), allocatable :: sig_p_c(:)     ! per-atom stationary momentum widths [amu*Ang/(10 fs)]
   real(8), allocatable :: sig_q_c(:)     ! per-atom stationary coordinate widths [Ang]
   real(8), allocatable :: ghost_q(:)     ! last drawn ghost coordinates [Ang] (3 per target atom)
   real(8), allocatable :: ghost_p(:)     ! last drawn ghost momenta [amu*Ang/(10 fs)]
   logical :: inited = .false.       ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! glo_target_init(p) - validate the bath parameters and the propagator
   !                      step, bind the carrier, load the coupling strengths
   subroutine glo_target_init(p)
      type(glo_target_params_t), intent(in) :: p ! member parameters (bath coupling)
      integer :: i, code, n_carry, i_carry

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
      if (.not. allocated(roster%frag)) then
         call stop_glox('glo_target_init', 'no roster fragment table - the roster must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_glox('glo_target_init', 'no per-fragment scheme array (reactants%'// &
            'dist_scheme unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_gloi2('glo_target_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('glo_target')
      if (code == 0) then
         call stop_glox('glo_target_init', 'member word "glo_target" is not in the '// &
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
         call stop_glox('glo_target_init', 'no roster fragment carries the glo_target '// &
            'scheme - this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_gloi2('glo_target_init', 'two or more fragments carry the glo_target '// &
            'scheme, but one member instance carries one parameter set; count =', n_carry, &
            'first carrier =', i_carry)
      end if

      ! 3. target selection + coupling-strength loading (the conservation
      !    algebra: sigma_noise^2 = 2*gamma*m*kb*t/dt; sigma_p^2 = m*kb*t;
      !    sigma_q^2 = kb*t/(m*w^2), per carrier atom)
      frag_b = i_carry
      nat_b = roster%frag(frag_b)%nat
      if (allocated(sig_noise_c)) deallocate (sig_noise_c)
      if (allocated(sig_p_c)) deallocate (sig_p_c)
      if (allocated(sig_q_c)) deallocate (sig_q_c)
      allocate (sig_noise_c(nat_b), sig_p_c(nat_b), sig_q_c(nat_b))
      do i = 1, nat_b
         associate (m_i => roster%mass(roster%frag(frag_b)%list(i)))
            sig_noise_c(i) = sqrt(2.0d0*p%gamma_glo*m_i*kb_code*p%t_glo/prop_cfg%dt)
            sig_p_c(i) = sqrt(m_i*kb_code*p%t_glo)
            sig_q_c(i) = sqrt(kb_code*p%t_glo/(m_i*p%w_ghost**2))
         end associate
      end do
      t_glo_c = p%t_glo
      gamma_c = p%gamma_glo
      w_c = p%w_ghost
      dt_c = prop_cfg%dt
      if (allocated(ghost_q)) deallocate (ghost_q)
      if (allocated(ghost_p)) deallocate (ghost_p)
      allocate (ghost_q(3*nat_b), ghost_p(3*nat_b))
      ghost_q = 0.0d0
      ghost_p = 0.0d0
      inited = .true.
   end subroutine glo_target_init

   !------------------------------------------------------------------
   ! glo_target_sample(sta) - draw the ghost stationary state (one fresh
   !                          Gaussian per ghost component); sta untouched
   subroutine glo_target_sample(sta)
      type(state_t), intent(inout) :: sta ! sampler runtime channel (no q/p/f components owned)
      integer :: i, c
      if (.not. inited) then
         call stop_glox('glo_target_sample', 'sample called before glo_target_init - the '// &
            'member carries no parameter set')
      end if
      do i = 1, nat_b
         do c = 1, 3
            ghost_q(3*(i - 1) + c) = sig_q_c(i)*rng_gauss()
            ghost_p(3*(i - 1) + c) = sig_p_c(i)*rng_gauss()
         end do
      end do
   end subroutine glo_target_sample

   !------------------------------------------------------------------
   ! getters (copies of the derived member state)
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

   pure integer function glo_target_frag()
      glo_target_frag = frag_b
   end function glo_target_frag

   pure integer function glo_target_nat()
      glo_target_nat = nat_b
   end function glo_target_nat

   pure integer function glo_target_ghost_size()
      glo_target_ghost_size = 3*nat_b
   end function glo_target_ghost_size

   pure real(8) function glo_target_last_q(k)
      integer, intent(in) :: k ! ghost component index (1..3*nat) [-]
      glo_target_last_q = ghost_q(k)
   end function glo_target_last_q

   pure real(8) function glo_target_last_p(k)
      integer, intent(in) :: k ! ghost component index (1..3*nat) [-]
      glo_target_last_p = ghost_p(k)
   end function glo_target_last_p

   pure real(8) function glo_target_sigma_noise(i)
      integer, intent(in) :: i ! carrier-atom index [-]
      glo_target_sigma_noise = sig_noise_c(i)
   end function glo_target_sigma_noise

   pure real(8) function glo_target_sigma_p(i)
      integer, intent(in) :: i ! carrier-atom index [-]
      glo_target_sigma_p = sig_p_c(i)
   end function glo_target_sigma_p

   pure real(8) function glo_target_sigma_q(i)
      integer, intent(in) :: i ! carrier-atom index [-]
      glo_target_sigma_q = sig_q_c(i)
   end function glo_target_sigma_q

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_glox(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_glox

   subroutine stop_gloz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_gloz

   subroutine stop_gloi2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_gloi2

end module dist_glo_target
