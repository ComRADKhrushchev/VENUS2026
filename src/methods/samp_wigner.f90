!=====================================================================
! samp_wigner.f90 - distribution member: Wigner normal-mode
!   distribution - independent Gaussian q and p per mode with the
!   virial-matched widths sigma_q = sqrt(E)/w, sigma_p = sqrt(E),
!   E = e_scale*(hbar/2)*w (the ground state at e_scale = 1)
! Design:
!   One e_scale serves every carrier equally; the mode tables arrive as
!   assembly data - one per carrier (each fragment is its own spectrum);
!   init validates shape and orthonormality loudly per carrier. A Wigner
!   draw is an energy DISTRIBUTION, not a fixed shell: the per-mode
!   energy E*(g1^2+g2^2)/2 is exponential with mean E, so the virial
!   pair <T> = <V> = E/2 holds in the mean (e_scale < 1 serves the
!   sub-zero-point sampling used to curb zero-point leak). Sample writes
!   q/p about the COM origin; the mode span excludes translations, so
!   the net momentum is zero - no angular-momentum stripping or
!   orientation randomization (the rotation member's channel).
!   Units: e_scale [-]; E [internal]; w_mode [rad/(10 fs)]; q [Angstrom];
!   p [amu*Ang/(10 fs)]; m [amu]; hbar = hbar_code [internal*(10 fs)].
!=====================================================================
module samp_wigner
   use consts,        only: hbar_code
   use rng,           only: rng_gauss
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code
   implicit none
   private
   public :: wigner_tbl_t, wigner_params_t, wigner_init, wigner_draw, wigner_realize, wigner_sample

   real(8), parameter :: orth_tol = 1.0d-8    ! C^T*C - I tolerance [-]
   real(8), parameter :: trans_tol = 1.0d-8   ! translation-orthogonality tolerance [-]

   ! One carrier's mode table (assembly data; fragments carry different
   ! spectra, so the tables are ragged per carrier)
   type :: wigner_tbl_t
      real(8), allocatable :: w(:)       ! mode frequencies [rad/(10 fs)]
      real(8), allocatable :: c(:,:)     ! (3*nat x n_mode) orthonormal mass-weighted columns
   end type wigner_tbl_t

   ! Member-owned parameters (the mode tables arrive as assembly data)
   type :: wigner_params_t
      real(8) :: e_scale = 1.0d0              ! per-mode energy scale, fraction of the zero-point energy [-] (every carrier)
      type(wigner_tbl_t), allocatable :: tbl(:) ! per-carrier mode tables (list_atoms order)
   end type wigner_params_t

   ! Derived member state (module-private; rebuilt by every wigner_init) -
   ! one row per carrier; the Gaussian widths fold the shared e_scale
   type :: carry_t
      integer :: frag = 0                     ! carrier list_atoms fragment row [index]
      integer :: nat = 0                      ! carrier atom count [count]
      integer :: n_mode = 0                   ! carrier mode count [count]
      integer, allocatable :: list(:)         ! carrier atom index list [global atom index]
      real(8), allocatable :: qz_com(:)       ! COM-shifted seed coordinates [Angstrom]
      real(8), allocatable :: sqm(:)          ! sqrt(mass) replicated per component [sqrt(amu)]
      real(8), allocatable :: c_c(:,:)        ! validated eigenvector table
      real(8), allocatable :: sig_q(:)        ! per-mode coordinate width [sqrt(amu)*Angstrom]
      real(8), allocatable :: sig_p(:)        ! per-mode momentum width [amu*Ang/(10 fs)]
      real(8), allocatable :: x_s(:), y_s(:)  ! the drawn Cartesian pair (post-draw record) [sqrt(amu)*Ang], [amu*Ang/(10 fs)]
   end type carry_t
   type(carry_t), allocatable :: carry(:)     ! carrier rows (list_atoms order)
   logical :: inited = .false.                ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! wigner_init(p) - validate the parameters and every carrier's mode
   !                  table, bind the carriers, precompute the widths
   subroutine wigner_init(p)
      type(wigner_params_t), intent(in) :: p  ! member parameters
      real(8) :: com(3), wt, dev, tv(3), tscale, e_mode
      integer :: i, k, d, kc, code, n_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%e_scale <= 0.0d0) then
         call stop_wgz('wigner_init', 'energy scale e_scale (fraction of the zero-point energy) must be '// &
            'positive, got', p%e_scale)
      end if

      ! 2. carrier binding (assembly-order guards; the atom list fragments whose
      !    dist_scheme rows select this member, list_atoms order)
      if (.not. allocated(list_atoms%frag)) then
         call stop_wgx('wigner_init', 'no list_atoms fragment table - the atom list must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_wgx('wigner_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its carriers')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_wgi2('wigner_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('wigner')
      if (code == 0) then
         call stop_wgx('wigner_init', 'member word "wigner" is not in the scheme-code table '// &
            '(incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_wgx('wigner_init', 'no list_atoms fragment carries the wigner scheme - '// &
            'this member would never be dispatched')
      end if
      if (.not. allocated(p%tbl) .or. size(p%tbl) /= n_carry) then
         call stop_wgi2('wigner_init', 'mode-table set mismatch: the parameter carries', &
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
               call stop_wgx('wigner_init', 'mode table unallocated - the frequencies and the '// &
                  'mass-weighted eigenvector columns are assembly data the member requires')
            end if
            if (size(p%tbl(kc)%c, 2) /= size(p%tbl(kc)%w)) then
               call stop_wgi2('wigner_init', 'mode table mismatch: eigenvector column count =', &
                  size(p%tbl(kc)%c, 2), 'but frequency count =', size(p%tbl(kc)%w))
            end if
            if (size(p%tbl(kc)%c, 1) /= 3*cr%nat) then
               call stop_wgi2('wigner_init', 'mode table shape mismatch: eigenvector row count =', &
                  size(p%tbl(kc)%c, 1), 'but the carrier carries 3*nat =', 3*cr%nat)
            end if
            cr%n_mode = size(p%tbl(kc)%w)
            if (cr%n_mode < 1) then
               call stop_wgi('wigner_init', 'empty mode table - no vibrational support to '// &
                  'distribute over; mode count =', cr%n_mode)
            end if
            if (cr%n_mode > 3*cr%nat - 3) then
               call stop_wgi2('wigner_init', 'mode count =', cr%n_mode, &
                  'leaves no room for the three translations (max 3*nat - 3 =', 3*cr%nat - 3)
            end if
            do k = 1, cr%n_mode
               if (p%tbl(kc)%w(k) <= 0.0d0) then
                  call stop_wgz('wigner_init', 'mode frequency w_mode [rad/(10 fs)] must be '// &
                     'positive; entry =', p%tbl(kc)%w(k))
               end if
            end do
            do k = 1, cr%n_mode
               do d = 1, cr%n_mode
                  dev = dot_product(p%tbl(kc)%c(:, k), p%tbl(kc)%c(:, d)) &
                        - merge(1.0d0, 0.0d0, k == d)
                  if (abs(dev) > orth_tol) then
                     call stop_wgz('wigner_init', 'mode table not orthonormal (C^T*C - I entry, '// &
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
                  call stop_wgz('wigner_init', 'mode column not orthogonal to the mass-weighted '// &
                     'translations (residual norm, tolerance 1e-8*sqrt(sum m)) =', norm2(tv))
               end if
            end do

            ! 4. precompute (COM-shifted seed, sqrt(mass) vector, widths)
            allocate (cr%list(cr%nat), cr%qz_com(3*cr%nat), cr%sqm(3*cr%nat), &
                      cr%c_c(3*cr%nat, cr%n_mode), cr%sig_q(cr%n_mode), cr%sig_p(cr%n_mode), &
                      cr%x_s(3*cr%nat), cr%y_s(3*cr%nat))
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
            cr%c_c = p%tbl(kc)%c
            do k = 1, cr%n_mode
               e_mode = p%e_scale*hbar_code*p%tbl(kc)%w(k)/2.0d0   ! E = e_scale*(hbar/2)*w
               cr%sig_p(k) = sqrt(e_mode)                          ! sigma_p = sqrt(E)
               cr%sig_q(k) = cr%sig_p(k)/p%tbl(kc)%w(k)            ! sigma_q = sqrt(E)/w
            end do
         end associate
      end do
      inited = .true.
   end subroutine wigner_init

   !------------------------------------------------------------------
   ! wigner_sample(sta) - one Wigner realization per carrier (list_atoms
   !                      order; independent Gaussian q and p per mode)
   !------------------------------------------------------------------
   ! wigner_draw() - per carrier: independent Gaussian mode q/p draws (the
   !                 historical interleaved order: one q then one p per mode)
   !                 plus the deterministic back-transformation into the
   !                 record; no state access
   !------------------------------------------------------------------
   subroutine wigner_draw()
      real(8), allocatable :: qmode(:), pmode(:)
      integer :: k, kc
      if (.not. inited) then
         call stop_wgx('wigner_draw', 'draw called before wigner_init - the member '// &
            'carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         allocate (qmode(cr%n_mode), pmode(cr%n_mode))
         do k = 1, cr%n_mode
            qmode(k) = cr%sig_q(k)*rng_gauss()
            pmode(k) = cr%sig_p(k)*rng_gauss()
         end do
         cr%x_s = matmul(cr%c_c, qmode)
         cr%y_s = matmul(cr%c_c, pmode)
         deallocate (qmode, pmode)
         end associate
      end do
   end subroutine wigner_draw

   !------------------------------------------------------------------
   ! wigner_realize(sta) - write the drawn Cartesian pair about the COM
   !                       origin (no RNG draws)
   !------------------------------------------------------------------
   subroutine wigner_realize(sta)
      type(state_t), intent(inout) :: sta  ! initial state (every carrier's q/p
                                           ! components are overwritten, list_atoms order)
      integer :: i, i3, kc
      if (.not. inited) then
         call stop_wgx('wigner_realize', 'realize called before wigner_init - the member '// &
            'carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         if (size(sta%q) < 3*maxval(cr%list) .or. size(sta%p) < 3*maxval(cr%list)) then
            call stop_wgx('wigner_realize', 'state too small for the carrier fragment atoms - '// &
               'state_create must follow the atom list')
         end if
         do i = 1, cr%nat
            i3 = 3*(cr%list(i) - 1)
            sta%q(i3+1:i3+3) = cr%qz_com(3*i-2:3*i) + cr%x_s(3*i-2:3*i)/cr%sqm(3*i-2:3*i)
            sta%p(i3+1:i3+3) = cr%sqm(3*i-2:3*i)*cr%y_s(3*i-2:3*i)
         end do
         end associate
      end do
   end subroutine wigner_realize

   !------------------------------------------------------------------
   ! wigner_sample(sta) - fused convenience: draw then realize (the
   !                       check-program entry)
   !------------------------------------------------------------------
   subroutine wigner_sample(sta)
      type(state_t), intent(inout) :: sta
      call wigner_draw()
      call wigner_realize(sta)
   end subroutine wigner_sample

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; STOP 1)
   subroutine stop_wgx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_wgx

   subroutine stop_wgz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_wgz

   subroutine stop_wgi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_wgi

   subroutine stop_wgi2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_wgi2

end module samp_wigner
