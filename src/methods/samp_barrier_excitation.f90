!=====================================================================
! samp_barrier_excitation.f90 - distribution member: saddle sampling
!   (e_stab split over given stable modes + reaction-coordinate
!   momentum kick, fresh-uniform sign, zero rc displacement)
! Design:
!   One parameter set serves every carrier equally; saddle data arrives
!   as assembly data - one table per carrier - and init validates it
!   loudly per carrier. Sample writes each carrier's q/p about its COM
!   origin (list_atoms order). The barrier-momentum sign is a fresh uniform
!   independent of the magnitude.
!   Units: e_stab/e_bar [kcal/mol], t_bar [K], E [internal], w_mode
!   [rad/(10 fs)], q [Angstrom], p [amu*Ang/(10 fs)], m [amu].
!=====================================================================
module samp_barrier_excitation
   use consts,        only: kb_code, e_conv, two_pi
   use rng,           only: rng_u, rng_gauss, rng_gamma
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code
   implicit none
   private
   public :: saddle_tbl_t, barrier_excitation_params_t, barrier_excitation_init, &
             barrier_excitation_draw, barrier_excitation_realize, barrier_excitation_sample

   real(8), parameter :: orth_tol = 1.0d-8    ! C^T*C - I / |c_rc| - 1 tolerance [-]
   real(8), parameter :: trans_tol = 1.0d-8   ! translation-orthogonality tolerance [-]

   ! One carrier's saddle data (assembly data; ragged per carrier)
   type :: saddle_tbl_t
      real(8), allocatable :: w(:)       ! stable-mode frequencies [rad/(10 fs)]
      real(8), allocatable :: c(:,:)     ! (3*nat x n_mode) orthonormal mass-weighted columns
      real(8), allocatable :: rc(:)      ! unit mass-weighted reaction coordinate (3*nat)
   end type saddle_tbl_t

   ! Member-owned parameters (the saddle tables arrive as assembly data)
   type :: barrier_excitation_params_t
      real(8) :: e_stab = 8.0d0               ! stabilization energy [kcal/mol] (>= 0, every carrier)
      integer :: n_e_bar = 0                  ! barrier energy law: 0 fixed e_bar, 1 thermal t_bar
      real(8) :: e_bar = 5.0d0                ! barrier energy [kcal/mol] (mode 0)
      real(8) :: t_bar = 900.0d0              ! barrier temperature [K] (mode 1)
      type(saddle_tbl_t), allocatable :: tbl(:)  ! per-carrier saddle tables (list_atoms order)
   end type

   ! Derived member state (module-private; rebuilt by every init) - one row
   ! per carrier; the energy-law constants are shared (equal treatment)
   type :: carry_t
      integer :: frag = 0                     ! carrier list_atoms fragment row [index]
      integer :: nat = 0                      ! carrier atom count [count]
      integer :: n_mode = 0                   ! carrier stable-mode count [count]
      integer, allocatable :: list(:)         ! carrier atom index list [global atom index]
      real(8), allocatable :: qz_com(:)       ! COM-shifted seed coordinates [Angstrom]
      real(8), allocatable :: sqm(:)          ! sqrt(mass) per component [sqrt(amu)]
      real(8), allocatable :: w_c(:), c_c(:,:), rc_c(:)
      real(8), allocatable :: x_s(:), y_s(:)  ! the drawn Cartesian pair incl. the rc kick (post-draw record)
   end type carry_t
   type(carry_t), allocatable :: carry(:)     ! carrier rows (list_atoms order)
   real(8) :: e_stab_int = 0.0d0              ! stabilization energy [internal]
   integer :: n_e_bar_c = 0                   ! validated barrier energy law [-]
   real(8) :: pbar_fix = 0.0d0                ! mode-0 barrier momentum magnitude [mass-weighted]
   real(8) :: t_bar_c = 0.0d0                 ! validated barrier temperature [K]
   logical :: inited = .false.                ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! barrier_excitation_init(p) - validate the parameters and every
   !                               carrier's saddle data, bind the carriers
   subroutine barrier_excitation_init(p)
      type(barrier_excitation_params_t), intent(in) :: p  ! member parameters
      real(8) :: com(3), wt, dev, tv(3), tscale
      integer :: i, k, d, kc, code, n_carry

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

      ! 2. carrier binding (assembly-order guards; list_atoms order)
      if (.not. allocated(list_atoms%frag)) then
         call stop_bex('barrier_excitation_init', 'no list_atoms fragment table - the atom list must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_bex('barrier_excitation_init', 'no per-fragment scheme array '// &
            '(reactants%dist_scheme unallocated) - member init cannot bind its carriers')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_bei2('barrier_excitation_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('barrier_excitation')
      if (code == 0) then
         call stop_bex('barrier_excitation_init', 'member word "barrier_excitation" is not in the '// &
            'scheme-code table (incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_bex('barrier_excitation_init', 'no list_atoms fragment carries the '// &
            'barrier_excitation scheme - this member would never be dispatched')
      end if
      if (.not. allocated(p%tbl) .or. size(p%tbl) /= n_carry) then
         call stop_bei2('barrier_excitation_init', 'saddle-table set mismatch: the parameter '// &
            'carries', merge(size(p%tbl), 0, allocated(p%tbl)), &
            'tables but the selection carries carriers', n_carry)
      end if
      if (allocated(carry)) deallocate (carry)
      allocate (carry(n_carry))
      kc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         kc = kc + 1
         associate (fr => list_atoms%frag(i), cr => carry(kc))
            cr%frag = i
            cr%nat = fr%nat

            ! 3. saddle-data validation for THIS carrier (shape, contents,
            !    orthonormality, reaction coordinate)
            if (.not. allocated(p%tbl(kc)%w) .or. .not. allocated(p%tbl(kc)%c) .or. &
                .not. allocated(p%tbl(kc)%rc)) then
               call stop_bex('barrier_excitation_init', 'saddle data unallocated - the stable-mode '// &
                  'table and the reaction coordinate are assembly data the member requires')
            end if
            if (size(p%tbl(kc)%c, 2) /= size(p%tbl(kc)%w)) then
               call stop_bei2('barrier_excitation_init', 'mode table mismatch: eigenvector column '// &
                  'count =', size(p%tbl(kc)%c, 2), 'but frequency count =', size(p%tbl(kc)%w))
            end if
            if (size(p%tbl(kc)%c, 1) /= 3*cr%nat) then
               call stop_bei2('barrier_excitation_init', 'mode table shape mismatch: eigenvector row '// &
                  'count =', size(p%tbl(kc)%c, 1), 'but the carrier carries 3*nat =', 3*cr%nat)
            end if
            cr%n_mode = size(p%tbl(kc)%w)
            if (cr%n_mode < 0) then
               call stop_bei('barrier_excitation_init', 'negative stable-mode count =', cr%n_mode)
            end if
            if (cr%n_mode > 3*cr%nat - 4) then
               call stop_bei2('barrier_excitation_init', 'stable-mode count =', cr%n_mode, &
                  'leaves no room for the three translations and the reaction coordinate (max =', &
                  3*cr%nat - 4)
            end if
            if (p%e_stab > 0.0d0 .and. cr%n_mode < 1) then
               call stop_bex('barrier_excitation_init', 'stabilization energy e_stab > 0 with an empty '// &
                  'stable-mode table - no stable support to absorb it')
            end if
            if (size(p%tbl(kc)%rc) /= 3*cr%nat) then
               call stop_bei2('barrier_excitation_init', 'reaction-coordinate length =', size(p%tbl(kc)%rc), &
                  'but the carrier carries 3*nat =', 3*cr%nat)
            end if
            do k = 1, cr%n_mode
               if (p%tbl(kc)%w(k) <= 0.0d0) then
                  call stop_bez('barrier_excitation_init', 'stable-mode frequency w_mode [rad/(10 fs)] '// &
                     'must be positive; entry =', p%tbl(kc)%w(k))
               end if
            end do
            tscale = sqrt(sum(list_atoms%mass(fr%list)))
            do k = 1, cr%n_mode
               do d = 1, cr%n_mode
                  dev = dot_product(p%tbl(kc)%c(:, k), p%tbl(kc)%c(:, d)) &
                        - merge(1.0d0, 0.0d0, k == d)
                  if (abs(dev) > orth_tol) then
                     call stop_bez('barrier_excitation_init', 'stable-mode table not orthonormal '// &
                        '(C^T*C - I entry, tolerance 1e-8) =', dev)
                  end if
               end do
               tv = 0.0d0
               do d = 1, cr%nat
                  tv = tv + sqrt(list_atoms%mass(fr%list(d)))*p%tbl(kc)%c(3*d-2:3*d, k)
               end do
               if (norm2(tv) > trans_tol*tscale) then
                  call stop_bez('barrier_excitation_init', 'stable-mode column not orthogonal to the '// &
                     'mass-weighted translations (residual norm, tolerance 1e-8*sqrt(sum m)) =', norm2(tv))
               end if
            end do
            if (abs(norm2(p%tbl(kc)%rc) - 1.0d0) > orth_tol) then
               call stop_bez('barrier_excitation_init', 'reaction coordinate not unit (|c_rc| - 1, '// &
                  'tolerance 1e-8) =', norm2(p%tbl(kc)%rc) - 1.0d0)
            end if
            tv = 0.0d0
            do d = 1, cr%nat
               tv = tv + sqrt(list_atoms%mass(fr%list(d)))*p%tbl(kc)%rc(3*d-2:3*d)
            end do
            if (norm2(tv) > trans_tol*tscale) then
               call stop_bez('barrier_excitation_init', 'reaction coordinate not orthogonal to the '// &
                  'mass-weighted translations (residual norm, tolerance 1e-8*sqrt(sum m)) =', norm2(tv))
            end if
            do k = 1, cr%n_mode
               dev = dot_product(p%tbl(kc)%c(:, k), p%tbl(kc)%rc)
               if (abs(dev) > orth_tol) then
                  call stop_bez('barrier_excitation_init', 'reaction coordinate overlaps the stable-mode '// &
                     'span (|C^T*c_rc| entry, tolerance 1e-8) =', dev)
               end if
            end do

            ! 4. precompute (COM-shifted seed, sqrt(mass) vector, table copies)
            allocate (cr%list(cr%nat), cr%qz_com(3*cr%nat), cr%sqm(3*cr%nat), &
                      cr%x_s(3*cr%nat), cr%y_s(3*cr%nat))
            allocate (cr%w_c(max(cr%n_mode, 1)), cr%c_c(3*cr%nat, max(cr%n_mode, 1)), &
                      cr%rc_c(3*cr%nat))
            cr%list = fr%list
            wt = sum(list_atoms%mass(cr%list))
            com = 0.0d0
            do d = 1, cr%nat
               com = com + list_atoms%mass(cr%list(d))*fr%qz(3*d-2:3*d)
            end do
            com = com/wt
            do d = 1, cr%nat
               cr%qz_com(3*d-2:3*d) = fr%qz(3*d-2:3*d) - com
               cr%sqm(3*d-2:3*d) = sqrt(list_atoms%mass(cr%list(d)))
            end do
            if (cr%n_mode >= 1) then
               cr%w_c(1:cr%n_mode) = p%tbl(kc)%w
               cr%c_c(:, 1:cr%n_mode) = p%tbl(kc)%c
            end if
            cr%rc_c = p%tbl(kc)%rc
         end associate
      end do
      e_stab_int = p%e_stab*e_conv
      n_e_bar_c = p%n_e_bar
      pbar_fix = sqrt(2.0d0*p%e_bar*e_conv)
      t_bar_c = p%t_bar
      inited = .true.
   end subroutine barrier_excitation_init

   !------------------------------------------------------------------
   ! barrier_excitation_sample(sta) - one realization per carrier (list_atoms
   !                               order): stable modes + barrier momentum,
   !                               fresh-uniform sign
   !------------------------------------------------------------------
   ! barrier_excitation_draw() - per carrier: the stable-mode sphere point
   !                              (Gaussian block), the per-mode uniform
   !                              phases, then the barrier momentum (fixed or
   !                              one-dimensional Maxwell) with its fresh
   !                              uniform sign - the historical order - plus
   !                              the deterministic back-transformation into
   !                              the record; no state access
   !------------------------------------------------------------------
   subroutine barrier_excitation_draw()
      real(8), allocatable :: g_vec(:), qmode(:), pmode(:)
      real(8) :: un, e_i, amp, ph, pbar
      integer :: k, kc
      if (.not. inited) then
         call stop_bex('barrier_excitation_draw', 'draw called before barrier_excitation_init '// &
            '- the member carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         allocate (g_vec(max(cr%n_mode, 1)), qmode(max(cr%n_mode, 1)), pmode(max(cr%n_mode, 1)))
         if (cr%n_mode >= 1) then
            do k = 1, cr%n_mode
               g_vec(k) = rng_gauss()
            end do
            un = norm2(g_vec(1:cr%n_mode))
            if (un < 1.0d-300) then
               g_vec = 0.0d0                  ! degenerate-draw guard (idle in practice)
               g_vec(1) = 1.0d0
               un = 1.0d0
            end if
            do k = 1, cr%n_mode
               e_i = e_stab_int*(g_vec(k)/un)**2
               amp = sqrt(2.0d0*e_i)/cr%w_c(k)
               ph = two_pi*rng_u()
               qmode(k) = amp*cos(ph)
               pmode(k) = -cr%w_c(k)*amp*sin(ph)
            end do
         end if
         if (n_e_bar_c == 0) then
            pbar = pbar_fix
         else
            pbar = sqrt(2.0d0*kb_code*t_bar_c*rng_gamma(1))
         end if
         if (rng_u() < 0.5d0) pbar = -pbar
         cr%x_s = 0.0d0
         cr%y_s = pbar*cr%rc_c
         if (cr%n_mode >= 1) then
            cr%x_s = matmul(cr%c_c(:, 1:cr%n_mode), qmode(1:cr%n_mode))
            cr%y_s = cr%y_s + matmul(cr%c_c(:, 1:cr%n_mode), pmode(1:cr%n_mode))
         end if
         deallocate (g_vec, qmode, pmode)
         end associate
      end do
   end subroutine barrier_excitation_draw

   !------------------------------------------------------------------
   ! barrier_excitation_realize(sta) - write the drawn pair incl. the rc kick
   !                                   about the COM origin (no RNG draws)
   !------------------------------------------------------------------
   subroutine barrier_excitation_realize(sta)
      type(state_t), intent(inout) :: sta  ! initial state (every carrier's q/p
                                           ! components are overwritten, list_atoms order)
      integer :: i, i3, kc
      if (.not. inited) then
         call stop_bex('barrier_excitation_realize', 'realize called before barrier_excitation_init '// &
            '- the member carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         if (size(sta%q) < 3*maxval(cr%list) .or. size(sta%p) < 3*maxval(cr%list)) then
            call stop_bex('barrier_excitation_realize', 'state too small for the carrier fragment '// &
               'atoms - state_create must follow the atom list')
         end if
         do i = 1, cr%nat
            i3 = 3*(cr%list(i) - 1)
            sta%q(i3+1:i3+3) = cr%qz_com(3*i-2:3*i) + cr%x_s(3*i-2:3*i)/cr%sqm(3*i-2:3*i)
            sta%p(i3+1:i3+3) = cr%sqm(3*i-2:3*i)*cr%y_s(3*i-2:3*i)
         end do
         end associate
      end do
   end subroutine barrier_excitation_realize

   !------------------------------------------------------------------
   ! barrier_excitation_sample(sta) - fused convenience: draw then realize
   !                                  (the check-program entry)
   !------------------------------------------------------------------
   subroutine barrier_excitation_sample(sta)
      type(state_t), intent(inout) :: sta
      call barrier_excitation_draw()
      call barrier_excitation_realize(sta)
   end subroutine barrier_excitation_sample

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; STOP 1)
   subroutine stop_bex(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_bex

   subroutine stop_bez(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_bez

   subroutine stop_bei(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_bei

   subroutine stop_bei2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_bei2

end module samp_barrier_excitation
