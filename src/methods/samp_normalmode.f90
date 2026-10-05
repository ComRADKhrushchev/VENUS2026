!=====================================================================
! samp_normalmode.f90 - distribution member: microcanonical normal-mode
!   sampling - fixed e_vib split over a given orthonormal mass-weighted
!   mode table (uniform sphere point, uniform phase), COM-origin output
! Design:
!   One e_vib serves every carrier equally; the mode tables arrive as
!   assembly data - one per carrier (each fragment is its own spectrum);
!   init validates shape and orthonormality loudly per carrier. Sample
!   places the full energy exactly once per carrier - no angular-momentum
!   stripping or orientation randomization (the rotation member's channel).
!   Units: e_vib [kcal/mol] via e_conv; E [internal]; w_mode
!   [rad/(10 fs)]; q [Angstrom]; p [amu*Ang/(10 fs)]; m [amu].
!=====================================================================
module samp_normalmode
   use consts,        only: e_conv, two_pi
   use rng,           only: rng_u, rng_gauss
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code
   implicit none
   private
   public :: mode_tbl_t, normalmode_params_t, normalmode_init, normalmode_draw, normalmode_realize, normalmode_sample

   real(8), parameter :: orth_tol = 1.0d-8    ! C^T*C - I tolerance [-]
   real(8), parameter :: trans_tol = 1.0d-8   ! translation-orthogonality tolerance [-]

   ! One carrier's mode table (assembly data; fragments carry different
   ! spectra, so the tables are ragged per carrier)
   type :: mode_tbl_t
      real(8), allocatable :: w(:)       ! mode frequencies [rad/(10 fs)]
      real(8), allocatable :: c(:,:)     ! (3*nat x n_mode) orthonormal mass-weighted columns
   end type mode_tbl_t

   ! Member-owned parameters (the mode tables arrive as assembly data)
   type :: normalmode_params_t
      real(8) :: e_vib = 10.0d0               ! total vibrational energy [kcal/mol] (every carrier)
      type(mode_tbl_t), allocatable :: tbl(:) ! per-carrier mode tables (list_atoms order)
   end type

   ! Derived member state (module-private; rebuilt by every normalmode_init) -
   ! one row per carrier; e_vib is shared (equal treatment)
   type :: carry_t
      integer :: frag = 0                     ! carrier list_atoms fragment row [index]
      integer :: nat = 0                      ! carrier atom count [count]
      integer :: n_mode = 0                   ! carrier mode count [count]
      integer, allocatable :: list(:)         ! carrier atom index list [global atom index]
      real(8), allocatable :: qz_com(:)       ! COM-shifted seed coordinates [Angstrom]
      real(8), allocatable :: sqm(:)          ! sqrt(mass) replicated per component [sqrt(amu)]
      real(8), allocatable :: w_c(:)          ! validated frequencies [rad/(10 fs)]
      real(8), allocatable :: c_c(:,:)        ! validated eigenvector table
      real(8), allocatable :: x_s(:), y_s(:)  ! the drawn Cartesian pair (post-draw record) [sqrt(amu)*Ang], [amu*Ang/(10 fs)]
   end type carry_t
   type(carry_t), allocatable :: carry(:)     ! carrier rows (list_atoms order)
   real(8) :: e_tot = 0.0d0                   ! total vibrational energy [internal]
   logical :: inited = .false.                ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! normalmode_init(p) - validate the parameters and every carrier's mode
   !                      table, bind the carriers, precompute the seeds
   subroutine normalmode_init(p)
      type(normalmode_params_t), intent(in) :: p  ! member parameters
      real(8) :: com(3), wt, dev, tv(3), tscale
      integer :: i, k, d, kc, code, n_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%e_vib <= 0.0d0) then
         call stop_nmz('normalmode_init', 'total vibrational energy e_vib [kcal/mol] must be '// &
            'positive, got', p%e_vib)
      end if

      ! 2. carrier binding (assembly-order guards; the atom list fragments whose
      !    dist_scheme rows select this member, list_atoms order)
      if (.not. allocated(list_atoms%frag)) then
         call stop_nmx('normalmode_init', 'no list_atoms fragment table - the atom list must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_nmx('normalmode_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its carriers')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_nmi2('normalmode_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('normalmode')
      if (code == 0) then
         call stop_nmx('normalmode_init', 'member word "normalmode" is not in the scheme-code '// &
            'table (incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_nmx('normalmode_init', 'no list_atoms fragment carries the normalmode scheme - '// &
            'this member would never be dispatched')
      end if
      if (.not. allocated(p%tbl) .or. size(p%tbl) /= n_carry) then
         call stop_nmi2('normalmode_init', 'mode-table set mismatch: the parameter carries', &
            merge(size(p%tbl), 0, allocated(p%tbl)), 'tables but the selection carries carriers', n_carry)
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

            ! 3. mode-table validation for THIS carrier (shape, contents,
            !    orthonormality, translation orthogonality)
            if (.not. allocated(p%tbl(kc)%w) .or. .not. allocated(p%tbl(kc)%c)) then
               call stop_nmx('normalmode_init', 'mode table unallocated - the frequencies and the '// &
                  'mass-weighted eigenvector columns are assembly data the member requires')
            end if
            if (size(p%tbl(kc)%c, 2) /= size(p%tbl(kc)%w)) then
               call stop_nmi2('normalmode_init', 'mode table mismatch: eigenvector column count =', &
                  size(p%tbl(kc)%c, 2), 'but frequency count =', size(p%tbl(kc)%w))
            end if
            if (size(p%tbl(kc)%c, 1) /= 3*cr%nat) then
               call stop_nmi2('normalmode_init', 'mode table shape mismatch: eigenvector row count =', &
                  size(p%tbl(kc)%c, 1), 'but the carrier carries 3*nat =', 3*cr%nat)
            end if
            cr%n_mode = size(p%tbl(kc)%w)
            if (cr%n_mode < 1) then
               call stop_nmi('normalmode_init', 'empty mode table - no vibrational support to '// &
                  'distribute over; mode count =', cr%n_mode)
            end if
            if (cr%n_mode > 3*cr%nat - 3) then
               call stop_nmi2('normalmode_init', 'mode count =', cr%n_mode, &
                  'leaves no room for the three translations (max 3*nat - 3 =', 3*cr%nat - 3)
            end if
            do k = 1, cr%n_mode
               if (p%tbl(kc)%w(k) <= 0.0d0) then
                  call stop_nmz('normalmode_init', 'mode frequency w_mode [rad/(10 fs)] must be '// &
                     'positive; entry =', p%tbl(kc)%w(k))
               end if
            end do
            do k = 1, cr%n_mode
               do d = 1, cr%n_mode
                  dev = dot_product(p%tbl(kc)%c(:, k), p%tbl(kc)%c(:, d)) &
                        - merge(1.0d0, 0.0d0, k == d)
                  if (abs(dev) > orth_tol) then
                     call stop_nmz('normalmode_init', 'mode table not orthonormal (C^T*C - I entry, '// &
                        'tolerance 1e-8) =', dev)
                  end if
               end do
            end do
            tscale = sqrt(sum(list_atoms%mass(fr%list)))
            do k = 1, cr%n_mode
               tv = 0.0d0
               do d = 1, cr%nat
                  tv = tv + sqrt(list_atoms%mass(fr%list(d)))*p%tbl(kc)%c(3*d-2:3*d, k)
               end do
               if (norm2(tv) > trans_tol*tscale) then
                  call stop_nmz('normalmode_init', 'mode column not orthogonal to the mass-weighted '// &
                     'translations (residual norm, tolerance 1e-8*sqrt(sum m)) =', norm2(tv))
               end if
            end do

            ! 4. precompute (COM-shifted seed, sqrt(mass) vector, table copies)
            allocate (cr%list(cr%nat), cr%qz_com(3*cr%nat), cr%sqm(3*cr%nat), &
                      cr%w_c(cr%n_mode), cr%c_c(3*cr%nat, cr%n_mode), cr%x_s(3*cr%nat), cr%y_s(3*cr%nat))
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
            cr%w_c = p%tbl(kc)%w
            cr%c_c = p%tbl(kc)%c
         end associate
      end do
      e_tot = p%e_vib*e_conv
      inited = .true.
   end subroutine normalmode_init

   !------------------------------------------------------------------
   ! normalmode_sample(sta) - one microcanonical realization per carrier
   !                          (list_atoms order; the same e_vib on every carrier)
   !------------------------------------------------------------------
   ! normalmode_draw() - per carrier: the uniform sphere point (one Gaussian
   !                     block), then the per-mode uniform phases (the
   !                     historical order), plus the deterministic
   !                     back-transformation into the record; no state access
   !------------------------------------------------------------------
   subroutine normalmode_draw()
      real(8), allocatable :: g_vec(:), qmode(:), pmode(:)
      real(8) :: un, e_i, amp, ph
      integer :: k, kc
      if (.not. inited) then
         call stop_nmx('normalmode_draw', 'draw called before normalmode_init - the member '// &
            'carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         allocate (g_vec(cr%n_mode), qmode(cr%n_mode), pmode(cr%n_mode))
         do k = 1, cr%n_mode
            g_vec(k) = rng_gauss()
         end do
         un = norm2(g_vec)
         if (un < 1.0d-300) then
            g_vec = 0.0d0                     ! degenerate-draw guard (idle in practice)
            g_vec(1) = 1.0d0
            un = 1.0d0
         end if
         do k = 1, cr%n_mode
            e_i = e_tot*(g_vec(k)/un)**2
            amp = sqrt(2.0d0*e_i)/cr%w_c(k)
            ph = two_pi*rng_u()
            qmode(k) = amp*cos(ph)
            pmode(k) = -cr%w_c(k)*amp*sin(ph)
         end do
         cr%x_s = matmul(cr%c_c, qmode)
         cr%y_s = matmul(cr%c_c, pmode)
         deallocate (g_vec, qmode, pmode)
         end associate
      end do
   end subroutine normalmode_draw

   !------------------------------------------------------------------
   ! normalmode_realize(sta) - write the drawn Cartesian pair about the COM
   !                           origin (no RNG draws)
   !------------------------------------------------------------------
   subroutine normalmode_realize(sta)
      type(state_t), intent(inout) :: sta  ! initial state (every carrier's q/p
                                           ! components are overwritten, list_atoms order)
      integer :: i, i3, kc
      if (.not. inited) then
         call stop_nmx('normalmode_realize', 'realize called before normalmode_init - the '// &
            'member carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         if (size(sta%q) < 3*maxval(cr%list) .or. size(sta%p) < 3*maxval(cr%list)) then
            call stop_nmx('normalmode_realize', 'state too small for the carrier fragment atoms - '// &
               'state_create must follow the atom list')
         end if
         do i = 1, cr%nat
            i3 = 3*(cr%list(i) - 1)
            sta%q(i3+1:i3+3) = cr%qz_com(3*i-2:3*i) + cr%x_s(3*i-2:3*i)/cr%sqm(3*i-2:3*i)
            sta%p(i3+1:i3+3) = cr%sqm(3*i-2:3*i)*cr%y_s(3*i-2:3*i)
         end do
         end associate
      end do
   end subroutine normalmode_realize

   !------------------------------------------------------------------
   ! normalmode_sample(sta) - fused convenience: draw then realize (the
   !                          check-program entry)
   !------------------------------------------------------------------
   subroutine normalmode_sample(sta)
      type(state_t), intent(inout) :: sta
      call normalmode_draw()
      call normalmode_realize(sta)
   end subroutine normalmode_sample

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; STOP 1)
   subroutine stop_nmx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_nmx

   subroutine stop_nmz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_nmz

   subroutine stop_nmi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_nmi

   subroutine stop_nmi2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_nmi2

end module samp_normalmode
