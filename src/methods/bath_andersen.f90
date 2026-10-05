!=====================================================================
! bath_andersen.f90 - bath species: Andersen velocity-reset - one sweep
!   draws, per atom, one collision trial (probability nu*dt) and on hit
!   replaces the atom's full velocity with a fresh Maxwell draw at
!   p(i,c) = sigma(i)*g, sigma(i)^2 = m_i*kb*T_bath
! Design:
!   The species covers the WHOLE system (one temperature for every atom -
!   the slab atoms of the surface paradigm are reachable precisely because
!   the scope is not the per-fragment scheme channel). It owns the exact
!   per-atom collision primitive only (no deterministic rescale, no
!   net-momentum zeroing - the loop's one whole-system drift removal is
!   bath.f90's business). The sweep signature is the bath-sweep contract
!   (sta + dt); the loop, provisioning and drift removal live in bath.f90.
!   Units: t_bath [K]; nu_coll [1/(10 fs)]; sigma [amu*Ang/(10 fs)].
!=====================================================================
module bath_andersen
   use consts,        only: kb_code
   use rng,           only: rng_u, rng_gauss
   use state,         only: state_t
   use config_atoms, only: list_atoms
   implicit none
   private
   public :: andersen_params_t, andersen_init, andersen_sweep, &
             andersen_t_bath, andersen_nu, andersen_sigma, andersen_nat

   ! Species-owned parameters
   type :: andersen_params_t
      real(8) :: t_bath = 0.0d0   ! bath temperature [K] (every atom)
      real(8) :: nu_coll = 0.0d0  ! collision frequency [1/(10 fs)]
   end type andersen_params_t

   ! Derived species state (module-private; rebuilt by every andersen_init)
   real(8), allocatable :: sig(:)  ! per-atom Maxwell widths sqrt(m*kb*T) [amu*Ang/(10 fs)]
   integer :: nat_c = 0            ! atom count [-]
   real(8) :: t_bath_c = 0.0d0     ! bath temperature [K]
   real(8) :: nu_coll_c = 0.0d0    ! collision frequency [1/(10 fs)]
   logical :: inited = .false.     ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! andersen_init(p) - validate the bath parameters, derive the
   !                      whole-system Maxwell widths from the atom list
   !------------------------------------------------------------------
   subroutine andersen_init(p)
      type(andersen_params_t), intent(in) :: p ! species parameters (bath)
      integer :: i

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%t_bath <= 0.0d0) then
         call stop_andz('andersen_init', 'bath temperature t_bath [K] must be positive, got', &
            p%t_bath)
      end if
      if (p%nu_coll <= 0.0d0) then
         call stop_andz('andersen_init', 'collision frequency nu_coll [1/(10 fs)] must be '// &
            'positive, got', p%nu_coll)
      end if

      ! 2. whole-system widths from the atom list masses (sigma^2 = m*kb*T)
      if (.not. allocated(list_atoms%mass)) then
         call stop_andx('andersen_init', 'no atom mass table - the atom list must be '// &
            'assembled before species init (assembly-order error)')
      end if
      nat_c = size(list_atoms%mass)
      if (allocated(sig)) deallocate (sig)
      allocate (sig(nat_c))
      do i = 1, nat_c
         sig(i) = sqrt(list_atoms%mass(i)*kb_code*p%t_bath)
      end do
      t_bath_c = p%t_bath
      nu_coll_c = p%nu_coll
      inited = .true.
   end subroutine andersen_init

   !------------------------------------------------------------------
   ! andersen_sweep(sta, dt) - one Andersen collision sweep over every
   !                             atom: one uniform trial per atom, on hit a
   !                             fresh Maxwell velocity for all three
   !                             components (the trial order is atom order)
   !------------------------------------------------------------------
   subroutine andersen_sweep(sta, dt, mask)
      type(state_t), intent(inout) :: sta ! physical state (collided atoms' momenta replaced)
      real(8), intent(in) :: dt           ! the sweep's step size [10 fs]
      logical, intent(in) :: mask(:)      ! per-atom collision scope [flag]
      integer :: i, i3
      if (.not. inited) then
         call stop_andx('andersen_sweep', 'sweep called before andersen_init - the '// &
            'species carries no parameter set')
      end if
      if (size(sta%p) < 3*nat_c) then
         call stop_andx('andersen_sweep', 'the state momentum vector ends before the last '// &
            'atom - state_create must follow the atom list')
      end if
      if (size(mask) < nat_c) then
         call stop_andx('andersen_sweep', 'the collision scope mask ends before the last '// &
            'atom - the bath stage must build it over the whole atom list')
      end if
      do i = 1, nat_c
         if (.not. mask(i)) cycle
         if (rng_u() < nu_coll_c*dt) then
            i3 = 3*(i - 1)
            sta%p(i3+1) = sig(i)*rng_gauss()
            sta%p(i3+2) = sig(i)*rng_gauss()
            sta%p(i3+3) = sig(i)*rng_gauss()
         end if
      end do
   end subroutine andersen_sweep

   !------------------------------------------------------------------
   ! getters (copies of the derived species state)
   !------------------------------------------------------------------
   pure real(8) function andersen_t_bath()
      andersen_t_bath = t_bath_c
   end function andersen_t_bath

   pure real(8) function andersen_nu()
      andersen_nu = nu_coll_c
   end function andersen_nu

   pure integer function andersen_nat()
      andersen_nat = nat_c
   end function andersen_nat

   pure real(8) function andersen_sigma(i)
      integer, intent(in) :: i             ! global atom index [-]
      andersen_sigma = 0.0d0
      if (inited .and. i >= 1 .and. i <= nat_c) andersen_sigma = sig(i)
   end function andersen_sigma

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; STOP 1)
   subroutine stop_andx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the bath loop is not usable)'
      stop 1
   end subroutine stop_andx

   subroutine stop_andz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the bath loop is not usable)'
      stop 1
   end subroutine stop_andz

end module bath_andersen
