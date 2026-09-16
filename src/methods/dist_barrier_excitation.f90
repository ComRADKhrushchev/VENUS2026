!=====================================================================
! dist_barrier_excitation.f90 - distribution member: saddle sampling
!   (e_stab split over given stable modes + reaction-coordinate
!   momentum kick, fresh-uniform sign, zero rc displacement)
! Design:
!   Saddle data arrives as assembly data; init validates it loudly and
!   sample writes the carrier's q/p about its COM origin. The barrier-
!   momentum sign is a fresh uniform independent of the magnitude.
!   Units: e_stab/e_bar [kcal/mol], t_bar [K], E [internal], w_mode
!   [rad/(10 fs)], q [Angstrom], p [amu*Ang/(10 fs)], m [amu].
!=====================================================================
module dist_barrier_excitation
   use consts,        only: kb_code, e_conv, two_pi
   use rng,           only: rng_u, rng_gauss, rng_gamma
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use sampler,       only: samp_code
   implicit none
   private
   public :: barrier_excitation_params_t, barrier_excitation_init, barrier_excitation_sample

   real(8), parameter :: orth_tol = 1.0d-8    ! C^T*C - I / |c_rc| - 1 tolerance [-]
   real(8), parameter :: trans_tol = 1.0d-8   ! translation-orthogonality tolerance [-]

   ! Member-owned parameters (the saddle data arrives as assembly data)
   type :: barrier_excitation_params_t
      real(8) :: e_stab = 8.0d0               ! stabilization energy [kcal/mol] (>= 0)
      integer :: n_e_bar = 0                  ! barrier energy law: 0 fixed e_bar, 1 thermal t_bar
      real(8) :: e_bar = 5.0d0                ! barrier energy [kcal/mol] (mode 0)
      real(8) :: t_bar = 900.0d0              ! barrier temperature [K] (mode 1)
      real(8), allocatable :: w_mode(:)       ! stable-mode frequencies [rad/(10 fs)]
      real(8), allocatable :: c_mode(:,:)     ! (3*nat x n_mode) orthonormal mass-weighted columns
      real(8), allocatable :: c_rc(:)         ! unit mass-weighted reaction coordinate (3*nat)
   end type

   ! Derived member state (module-private; rebuilt by every init)
   integer :: frag_served = 0                 ! bound roster fragment index [index]
   integer :: nat_c = 0                       ! carrier atom count [count]
   integer :: n_mode_c = 0                    ! stable-mode count [count]
   integer, allocatable :: list_c(:)          ! carrier atom index list [global atom index]
   real(8), allocatable :: qz_com(:)          ! COM-shifted seed coordinates [Angstrom]
   real(8), allocatable :: sqm(:)             ! sqrt(mass) per component [sqrt(amu)]
   real(8), allocatable :: w_c(:), c_c(:,:), rc_c(:)
   real(8) :: e_stab_int = 0.0d0              ! stabilization energy [internal]
   integer :: n_e_bar_c = 0                   ! validated barrier energy law [-]
   real(8) :: pbar_fix = 0.0d0                ! mode-0 barrier momentum magnitude [mass-weighted]
   real(8) :: t_bar_c = 0.0d0                 ! validated barrier temperature [K]
   logical :: inited = .false.                ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! barrier_excitation_init(p) - validate the parameters and saddle
   !                               data, bind the carrier, precompute
   subroutine barrier_excitation_init(p)
      type(barrier_excitation_params_t), intent(in) :: p  ! member parameters
      real(8) :: com(3), wt, dev, tv(3), tscale
      integer :: i, k, code, n_carry, i_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%e_stab < 0.0d0) then
         call stop_bez('barrier_excitation_init', 'stabilization energy e_stab [kcal/mol] must '// &
            'be >= 0 (zero = pure barrier kick), got', p%e_stab)
      end if
      if (p%n_e_bar /= 0 .and. p%n_e_bar /= 1) then
         call stop_bei('barrier_excitation_init', 'barrier energy-law selector n_e_bar must be 0 '// &
            '(fixed e_bar) or 1 (thermal t_bar), got', p%n_e_bar)
      end if
      if (p%n_e_bar == 0 .and. p%e_bar < 0.0d0) then
         call stop_bez('barrier_excitation_init', 'barrier energy e_bar [kcal/mol] must be >= 0 '// &
            'in the fixed mode, got', p%e_bar)
      end if
      if (p%n_e_bar == 1 .and. p%t_bar <= 0.0d0) then
         call stop_bez('barrier_excitation_init', 'barrier temperature t_bar [K] must be positive '// &
            'in the thermal mode, got', p%t_bar)
      end if

      ! 2. fragment binding (assembly-order and carrier guards first)
      if (.not. allocated(roster%frag)) then
         call stop_bex('barrier_excitation_init', 'no roster fragment table - the roster must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_bex('barrier_excitation_init', 'no per-fragment scheme array '// &
            '(reactants%dist_scheme unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_bei2('barrier_excitation_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('barrier_excitation')
      if (code == 0) then
         call stop_bex('barrier_excitation_init', 'member word "barrier_excitation" is not in the '// &
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
         call stop_bex('barrier_excitation_init', 'no roster fragment carries the '// &
            'barrier_excitation scheme - this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_bei2('barrier_excitation_init', 'two or more fragments carry the '// &
            'barrier_excitation scheme, but one member instance carries one parameter set; '// &
            'count =', n_carry, 'first carrier =', i_carry)
      end if
      frag_served = i_carry
      nat_c = roster%frag(i_carry)%nat

      ! 3. saddle-data validation (shape, contents, orthonormality, rc)
      if (.not. allocated(p%w_mode) .or. .not. allocated(p%c_mode) .or. .not. allocated(p%c_rc)) then
         call stop_bex('barrier_excitation_init', 'saddle data unallocated - the stable-mode table '// &
            'and the reaction coordinate are assembly data the member requires')
      end if
      if (size(p%c_mode, 2) /= size(p%w_mode)) then
         call stop_bei2('barrier_excitation_init', 'mode table mismatch: eigenvector column '// &
            'count =', size(p%c_mode, 2), 'but frequency count =', size(p%w_mode))
      end if
      if (size(p%c_mode, 1) /= 3*nat_c) then
         call stop_bei2('barrier_excitation_init', 'mode table shape mismatch: eigenvector row '// &
            'count =', size(p%c_mode, 1), 'but the carrier carries 3*nat =', 3*nat_c)
      end if
      n_mode_c = size(p%w_mode)
      if (n_mode_c < 0) then
         call stop_bei('barrier_excitation_init', 'negative stable-mode count =', n_mode_c)
      end if
      if (n_mode_c > 3*nat_c - 4) then
         call stop_bei2('barrier_excitation_init', 'stable-mode count =', n_mode_c, &
            'leaves no room for the three translations and the reaction coordinate (max =', &
            3*nat_c - 4)
      end if
      if (p%e_stab > 0.0d0 .and. n_mode_c < 1) then
         call stop_bex('barrier_excitation_init', 'stabilization energy e_stab > 0 with an empty '// &
            'stable-mode table - no stable support to absorb it')
      end if
      if (size(p%c_rc) /= 3*nat_c) then
         call stop_bei2('barrier_excitation_init', 'reaction-coordinate length =', size(p%c_rc), &
            'but the carrier carries 3*nat =', 3*nat_c)
      end if
      do k = 1, n_mode_c
         if (p%w_mode(k) <= 0.0d0) then
            call stop_bez('barrier_excitation_init', 'stable-mode frequency w_mode [rad/(10 fs)] '// &
               'must be positive; entry =', p%w_mode(k))
         end if
      end do
      tscale = sqrt(sum(roster%mass(roster%frag(i_carry)%list)))
      do k = 1, n_mode_c
         do i = 1, n_mode_c
            dev = dot_product(p%c_mode(:, k), p%c_mode(:, i)) - merge(1.0d0, 0.0d0, k == i)
            if (abs(dev) > orth_tol) then
               call stop_bez('barrier_excitation_init', 'stable-mode table not orthonormal '// &
                  '(C^T*C - I entry, tolerance 1e-8) =', dev)
            end if
         end do
         tv = 0.0d0
         do i = 1, nat_c
            tv = tv + sqrt(roster%mass(roster%frag(i_carry)%list(i)))*p%c_mode(3*i-2:3*i, k)
         end do
         if (norm2(tv) > trans_tol*tscale) then
            call stop_bez('barrier_excitation_init', 'stable-mode column not orthogonal to the '// &
               'mass-weighted translations (residual norm, tolerance 1e-8*sqrt(sum m)) =', norm2(tv))
         end if
      end do
      if (abs(norm2(p%c_rc) - 1.0d0) > orth_tol) then
         call stop_bez('barrier_excitation_init', 'reaction coordinate not unit (|c_rc| - 1, '// &
            'tolerance 1e-8) =', norm2(p%c_rc) - 1.0d0)
      end if
      tv = 0.0d0
      do i = 1, nat_c
         tv = tv + sqrt(roster%mass(roster%frag(i_carry)%list(i)))*p%c_rc(3*i-2:3*i)
      end do
      if (norm2(tv) > trans_tol*tscale) then
         call stop_bez('barrier_excitation_init', 'reaction coordinate not orthogonal to the '// &
            'mass-weighted translations (residual norm, tolerance 1e-8*sqrt(sum m)) =', norm2(tv))
      end if
      do k = 1, n_mode_c
         dev = dot_product(p%c_mode(:, k), p%c_rc)
         if (abs(dev) > orth_tol) then
            call stop_bez('barrier_excitation_init', 'reaction coordinate overlaps the stable-mode '// &
            'span (|C^T*c_rc| entry, tolerance 1e-8) =', dev)
         end if
      end do

      ! 4. precompute (COM-shifted seed, sqrt(mass) vector, constants)
      associate (fr => roster%frag(i_carry))
         if (allocated(list_c)) deallocate (list_c)
         if (allocated(qz_com)) deallocate (qz_com)
         if (allocated(sqm)) deallocate (sqm)
         if (allocated(w_c)) deallocate (w_c)
         if (allocated(c_c)) deallocate (c_c)
         if (allocated(rc_c)) deallocate (rc_c)
         allocate (list_c(nat_c), qz_com(3*nat_c), sqm(3*nat_c))
         allocate (w_c(max(n_mode_c, 1)), c_c(3*nat_c, max(n_mode_c, 1)), rc_c(3*nat_c))
         list_c = fr%list
         wt = sum(roster%mass(list_c))
         com = 0.0d0
         do i = 1, nat_c
            com = com + roster%mass(list_c(i))*fr%qz(3*i-2:3*i)
         end do
         com = com/wt
         do i = 1, nat_c
            qz_com(3*i-2:3*i) = fr%qz(3*i-2:3*i) - com
            sqm(3*i-2:3*i) = sqrt(roster%mass(list_c(i)))
         end do
      end associate
      if (n_mode_c >= 1) then
         w_c(1:n_mode_c) = p%w_mode
         c_c(:, 1:n_mode_c) = p%c_mode
      end if
      rc_c = p%c_rc
      e_stab_int = p%e_stab*e_conv
      n_e_bar_c = p%n_e_bar
      pbar_fix = sqrt(2.0d0*p%e_bar*e_conv)
      t_bar_c = p%t_bar
      inited = .true.
   end subroutine barrier_excitation_init

   !------------------------------------------------------------------
   ! barrier_excitation_sample(sta) - one realization (stable modes +
   !                              barrier momentum, fresh-uniform sign)
   subroutine barrier_excitation_sample(sta)
      type(state_t), intent(inout) :: sta  ! initial state (the carrier's q/p
                                           ! components are overwritten)
      real(8) :: g_vec(n_mode_c), un, e_i, amp, ph, pbar
      real(8) :: qmode(n_mode_c), pmode(n_mode_c), x_vec(3*nat_c), y_vec(3*nat_c)
      integer :: k, i, i3

      if (.not. inited) then
         call stop_bex('barrier_excitation_sample', 'sample called before barrier_excitation_init '// &
            '- the member carries no parameter set')
      end if
      if (size(sta%q) < 3*maxval(list_c) .or. size(sta%p) < 3*maxval(list_c)) then
         call stop_bex('barrier_excitation_sample', 'state too small for the carrier fragment '// &
            'atoms - state_create must follow the roster')
      end if

      ! 1. stable-mode microcanonical split (as the normalmode member)
      if (n_mode_c >= 1) then
         do k = 1, n_mode_c
            g_vec(k) = rng_gauss()
         end do
         un = norm2(g_vec)
         if (un < 1.0d-300) then
            g_vec = 0.0d0                  ! degenerate-draw guard (idle in practice)
            g_vec(1) = 1.0d0
            un = 1.0d0
         end if
         do k = 1, n_mode_c
            e_i = e_stab_int*(g_vec(k)/un)**2
            amp = sqrt(2.0d0*e_i)/w_c(k)
            ph = two_pi*rng_u()
            qmode(k) = amp*cos(ph)
            pmode(k) = -w_c(k)*amp*sin(ph)
         end do
      end if

      ! 2. barrier momentum along the reaction coordinate: fixed or
      !    one-dimensional Maxwell; the SIGN from a fresh uniform
      !    (independent of the magnitude deviate)
      if (n_e_bar_c == 0) then
         pbar = pbar_fix
      else
         pbar = sqrt(2.0d0*kb_code*t_bar_c*rng_gamma(1))
      end if
      if (rng_u() < 0.5d0) pbar = -pbar

      ! 3. back-transformation about the COM origin (zero displacement
      !    along the reaction coordinate)
      x_vec = 0.0d0
      y_vec = pbar*rc_c
      if (n_mode_c >= 1) then
         x_vec = matmul(c_c(:, 1:n_mode_c), qmode)
         y_vec = y_vec + matmul(c_c(:, 1:n_mode_c), pmode)
      end if
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         sta%q(i3+1:i3+3) = qz_com(3*i-2:3*i) + x_vec(3*i-2:3*i)/sqm(3*i-2:3*i)
         sta%p(i3+1:i3+3) = sqm(3*i-2:3*i)*y_vec(3*i-2:3*i)
      end do
   end subroutine barrier_excitation_sample

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_bex(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_bex

   subroutine stop_bez(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_bez

   subroutine stop_bei(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_bei

   subroutine stop_bei2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_bei2

end module dist_barrier_excitation
