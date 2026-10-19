!   The member carries its own registration (member_register /
!   member_init - the fixed-name entries of the 2026-10-14 weave batch;
!   the retired reg_samp_*.f90's _all/_init bodies live here now).
!=====================================================================
! samp_boltzmann.f90 - distribution member: Boltzmann vibrational
!   distribution - per-mode quantum numbers drawn from the thermal
!   geometric law P(n) = (1-q)*q^n, q = exp(-(h*c/k)*nu/T), then the
!   quasi-classical realization: per mode E_k = (n_k + 1/2)*hbar*w_k on
!   the energy shell with a uniform phase, back-transformed about the COM
! Design:
!   One temperature row per word-carrying fragment (the 2026-10-13 bracket
!   ruling realized per carrier - two BOLTZMANN fragments may carry different
!   T_VIB; a monoatomic word-carrier's row is validated but carries no
!   table); the FULL mode tables (frequencies in rad/fs plus orthonormal
!   mass-weighted columns) arrive as assembly data - one per polyatomic
!   carrier - and init validates shape and orthonormality loudly per
!   carrier (the normalmode member's contract). The draw law keeps the
!   legislated constant hc_k verbatim over wavenumbers derived from the
!   table (w_to_wvn); the draw consumes [one gamma per mode] then [one
!   uniform phase per mode]. The realize writes q/p from the drawn record
!   (no RNG): mode amplitudes sqrt(2E)/w at the drawn phases - each mode
!   exactly on its (n+1/2)*hbar*w shell (QCT standard: the zero point is
!   included), COM-origin output, zero net momentum (the mode span
!   excludes translations).
!   Units: nu [cm^-1], w [rad/fs], T [K], E [eV]; hbar =
!   hbar_evfs; hc/k = 1.43878 cm*K.
!=====================================================================
module samp_boltzmann
   use consts,        only: hbar_evfs, ev_dyn
   use rng,           only: rng_gamma, rng_u
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code, samp_reg
   use spectrum,      only: w_to_wvn
   use input, only: par_real
   use spectrum_interface, only: spectrum_frag_w
   use diag,        only: diag_fatal, diag_cat_sampler
   implicit none
   private
   public :: mode_tbl_t, boltz_params_t, boltz_init, boltz_draw, boltz_realize, &
             boltz_sample, boltz_last_n, member_register, member_init

   real(8), parameter :: hc_k = 1.43878d0   ! h*c/k [cm*K] - the geometric step constant
   real(8), parameter :: orth_tol = 1.0d-8  ! C^T*C - I tolerance [-]
   real(8), parameter :: trans_tol = 1.0d-8 ! translation-orthogonality tolerance [-]
   real(8), parameter :: two_pi = 6.283185307179586476925286766559d0

   ! One carrier's mode table (assembly data; fragments carry different
   ! mode counts, so the tables are ragged per carrier)
   type :: mode_tbl_t
      real(8), allocatable :: w(:)       ! mode frequencies [rad/fs]
      real(8), allocatable :: c(:,:)     ! (3*nat x n_mode) orthonormal mass-weighted columns
   end type mode_tbl_t

   ! Member-owned parameters (one temperature row per word-carrying fragment -
   ! the bracket ruling's per-carrier form; 0 = the name was absent and the
   ! member's own validation rejects it)
   type :: boltz_params_t
      real(8), allocatable :: t_vib(:)    ! vibrational temperatures [K] (word-carrier rows)
      type(mode_tbl_t), allocatable :: tbl(:) ! per-polyatomic-carrier mode tables (list_atoms order)
   end type boltz_params_t

   ! Derived member state (module-private; rebuilt by every boltz_init) -
   ! one row per carrier: the carrier fragment, its ladder constants, its
   ! geometry seeds, and the last drawn record (quantum numbers + phases)
   type :: carry_t
      integer :: frag = 0                   ! carrier list_atoms fragment row [-]
      integer :: nat = 0                    ! carrier atom count [count]
      integer, allocatable :: list(:)       ! carrier atom index list [global atom index]
      real(8), allocatable :: qz_com(:)     ! COM-shifted seed coordinates [Angstrom]
      real(8), allocatable :: sqm(:)        ! sqrt(mass) per component [sqrt(amu)]
      real(8), allocatable :: w_c(:)        ! validated frequencies [rad/fs]
      real(8), allocatable :: c_c(:,:)      ! validated eigenvector table
      real(8), allocatable :: dum(:)        ! per-mode geometric step (h*c/k)*nu/T [-]
      integer, allocatable :: n_last(:)     ! last drawn quantum-number table [-]
      real(8), allocatable :: ph_last(:)    ! last drawn phases [rad]
   end type carry_t
   type(carry_t), allocatable :: carry(:)   ! carrier rows (list_atoms order)
   logical :: inited = .false.              ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! boltz_init(p) - validate parameters, validate every carrier's mode
   !                 table, bind the carriers, precompute the seeds and
   !                 the per-mode geometric constants
   !------------------------------------------------------------------
   subroutine boltz_init(p)
      type(boltz_params_t), intent(in) :: p ! member parameters (one temperature
                                            ! row per word carrier / one mode
                                            ! table per polyatomic carrier)
      integer :: i, k, d, kc, wc, code, n_carry, n_word
      real(8) :: com(3), wt, dev, tv(3), tscale

      ! 1. parameter validation (fail-loud, naming field + value): one positive
      !    temperature per word carrier, one nonempty well-shaped table per
      !    selected polyatomic carrier
      if (.not. allocated(p%t_vib)) then
         call stop_boltx('boltz_init', 'no temperature rows (t_vib unallocated) - the reg '// &
            'must stage one row per word-carrying fragment')
      end if
      do wc = 1, size(p%t_vib)
         if (p%t_vib(wc) <= 0.0d0) then
            call stop_boltz('boltz_init', 'vibrational temperature t_vib [K] must be positive '// &
               '(the carrier wrote no T_VIB bracket name); row =', p%t_vib(wc))
         end if
      end do
      if (.not. allocated(p%tbl)) then
         call stop_boltx('boltz_init', 'no mode tables (tbl unallocated) - the frequency-'// &
            'derivation chain produced no tables; a zero spectrum is not a sampleable state')
      end if
      do k = 1, size(p%tbl)
         if (.not. allocated(p%tbl(k)%w) .or. .not. allocated(p%tbl(k)%c)) then
            call stop_boltx('boltz_init', 'mode table unallocated - the frequencies and the '// &
               'mass-weighted eigenvector columns are assembly data the member requires')
         end if
         if (size(p%tbl(k)%c, 2) /= size(p%tbl(k)%w)) then
            call stop_bolti2('boltz_init', 'mode table mismatch: eigenvector column count =', &
               size(p%tbl(k)%c, 2), 'but frequency count =', size(p%tbl(k)%w))
         end if
         do i = 1, size(p%tbl(k)%w)
            if (p%tbl(k)%w(i) <= 0.0d0) then
               call stop_boltz('boltz_init', 'mode frequency w [rad/fs] must be '// &
                  'positive; entry =', p%tbl(k)%w(i))
            end if
         end do
      end do

      ! 2. carrier binding (the list_atoms fragments whose dist_scheme rows
      !    select this member, list_atoms order; the temperature rows index
      !    every word carrier, the tables only the polyatomic ones)
      if (.not. allocated(list_atoms%frag)) then
         call stop_boltx('boltz_init', 'no list_atoms fragment table - the atom list must be '// &
            'assembled before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_boltx('boltz_init', 'no per-fragment scheme array (reactants%'// &
            'dist_scheme unallocated) - member init cannot bind its carriers')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_bolti2('boltz_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('boltzmann')
      if (code == 0) then
         call stop_boltx('boltz_init', 'member word "boltzmann" is not in the scheme-code '// &
            'table (incomplete member library assembly)')
      end if
      n_word = count(reactants%dist_scheme == code)
      if (n_word == 0) then
         call stop_boltx('boltz_init', 'no list_atoms fragment carries the boltzmann scheme - '// &
            'this member would never be dispatched')
      end if
      if (size(p%t_vib) /= n_word) then
         call stop_bolti2('boltz_init', 'the parameter set carries', size(p%t_vib), &
            'temperature rows but the selection carries word carriers', n_word)
      end if
      ! count only polyatomic carriers (a monoatomic fragment has no internal
      ! DOFs; the sampler skips it and the reg builds no table for it)
      n_carry = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code .and. list_atoms%frag(i)%nat > 1) then
            n_carry = n_carry + 1
         end if
      end do
      if (size(p%tbl) /= n_carry) then
         call stop_bolti2('boltz_init', 'the parameter set carries', size(p%tbl), &
            'mode tables but the selection carries carriers', n_carry)
      end if
      if (allocated(carry)) deallocate (carry)
      allocate (carry(n_carry))
      kc = 0
      wc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         wc = wc + 1                          ! the temperature row of THIS word carrier
         if (list_atoms%frag(i)%nat < 2) cycle   ! monoatomic: no table, skip
         kc = kc + 1
         associate (fr => list_atoms%frag(i), cr => carry(kc))
            cr%frag = i
            cr%nat = fr%nat

            ! 3. mode-table validation for THIS carrier (shape, contents,
            !    orthonormality, translation orthogonality)
            if (size(p%tbl(kc)%c, 1) /= 3*cr%nat) then
               call stop_bolti2('boltz_init', 'mode table shape mismatch: eigenvector row '// &
                  'count =', size(p%tbl(kc)%c, 1), 'but the carrier carries 3*nat =', 3*cr%nat)
            end if
            if (size(p%tbl(kc)%w) < 1) then
               call stop_bolti('boltz_init', 'empty mode table - no vibrational support to '// &
                  'distribute over; mode count =', size(p%tbl(kc)%w))
            end if
            do k = 1, size(p%tbl(kc)%w)
               do d = 1, size(p%tbl(kc)%w)
                  dev = dot_product(p%tbl(kc)%c(:, k), p%tbl(kc)%c(:, d)) &
                        - merge(1.0d0, 0.0d0, k == d)
                  if (abs(dev) > orth_tol) then
                     call stop_boltz('boltz_init', 'mode table not orthonormal (C^T*C - I '// &
                        'entry, tolerance 1e-8) =', dev)
                  end if
               end do
            end do
            tscale = sqrt(sum(list_atoms%mass(fr%list)))
            do k = 1, size(p%tbl(kc)%w)
               tv = 0.0d0
               do d = 1, cr%nat
                  tv = tv + sqrt(list_atoms%mass(fr%list(d)))*p%tbl(kc)%c(3*d-2:3*d, k)
               end do
               if (norm2(tv) > trans_tol*tscale) then
                  call stop_boltz('boltz_init', 'mode column not orthogonal to the '// &
                     'mass-weighted translations (residual norm, tol 1e-8*sqrt(sum m)) =', &
                     norm2(tv))
               end if
            end do

            ! 4. precompute (COM-shifted seed, sqrt(mass) vector, table
            !    copies, ladder constants)
            allocate (cr%list(cr%nat), cr%qz_com(3*cr%nat), cr%sqm(3*cr%nat), &
                      cr%w_c(size(p%tbl(kc)%w)), cr%c_c(3*cr%nat, size(p%tbl(kc)%w)), &
                      cr%dum(size(p%tbl(kc)%w)), cr%n_last(size(p%tbl(kc)%w)), &
                      cr%ph_last(size(p%tbl(kc)%w)))
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
            do k = 1, size(cr%w_c)
               cr%dum(k) = hc_k*w_to_wvn(cr%w_c(k))/p%t_vib(wc)
            end do
            cr%n_last = 0
            cr%ph_last = 0.0d0
         end associate
      end do
      inited = .true.
   end subroutine boltz_init

   !------------------------------------------------------------------
   ! boltz_draw() - one draw per carrier (list_atoms order): the quantum
   !                numbers first (one gamma per mode - the historical
   !                consumption order), then the phases (one uniform per
   !                mode); no state access
   !------------------------------------------------------------------
   subroutine boltz_draw()
      integer :: k, kc
      if (.not. inited) then
         call stop_boltx('boltz_draw', 'draw called before boltz_init - the member carries '// &
            'no parameter set')
      end if
      do kc = 1, size(carry)
         do k = 1, size(carry(kc)%dum)
            carry(kc)%n_last(k) = int(rng_gamma(1)/carry(kc)%dum(k))   ! floor(g/dum): P(n)=(1-q)*q^n exactly
         end do
      end do
      do kc = 1, size(carry)
         do k = 1, size(carry(kc)%dum)
            carry(kc)%ph_last(k) = two_pi*rng_u()
         end do
      end do
   end subroutine boltz_draw

   !------------------------------------------------------------------
   ! boltz_realize(sta) - the quasi-classical realization per carrier:
   !                       every mode exactly on its (n + 1/2)*hbar*w shell
   !                       at the drawn phase, back-transformed about the
   !                       COM origin (no RNG draws)
   !------------------------------------------------------------------
   subroutine boltz_realize(sta)
      type(state_t), intent(inout) :: sta  ! initial state (every carrier's q/p
                                           ! components are overwritten, list_atoms order)
      real(8), allocatable :: qmode(:), pmode(:), x_vec(:), y_vec(:)
      real(8) :: e_i, amp
      integer :: k, i, i3, kc
      if (.not. inited) then
         call stop_boltx('boltz_realize', 'realize called before boltz_init - the member '// &
            'carries no parameter set')
      end if
      do kc = 1, size(carry)
         associate (cr => carry(kc))
         if (size(sta%q) < 3*maxval(cr%list) .or. size(sta%p) < 3*maxval(cr%list)) then
            call stop_boltx('boltz_realize', 'state too small for the carrier fragment '// &
               'atoms - state_create must follow the atom list')
         end if
         allocate (qmode(size(cr%w_c)), pmode(size(cr%w_c)), x_vec(3*cr%nat), y_vec(3*cr%nat))
         do k = 1, size(cr%w_c)
            e_i = (dble(cr%n_last(k)) + 0.5d0)*hbar_evfs*cr%w_c(k)
            amp = sqrt(2.0d0*e_i*ev_dyn)/cr%w_c(k)
            qmode(k) = amp*cos(cr%ph_last(k))
            pmode(k) = -cr%w_c(k)*amp*sin(cr%ph_last(k))
         end do
         x_vec = matmul(cr%c_c, qmode)
         y_vec = matmul(cr%c_c, pmode)
         do i = 1, cr%nat
            i3 = 3*(cr%list(i) - 1)
            sta%q(i3+1:i3+3) = cr%qz_com(3*i-2:3*i) + x_vec(3*i-2:3*i)/cr%sqm(3*i-2:3*i)
            sta%p(i3+1:i3+3) = cr%sqm(3*i-2:3*i)*y_vec(3*i-2:3*i)
         end do
         deallocate (qmode, pmode, x_vec, y_vec)
         end associate
      end do
   end subroutine boltz_realize

   !------------------------------------------------------------------
   ! boltz_sample(sta) - fused convenience: draw then realize (the
   !                       check-program entry; production dispatch calls the
   !                       two halves separately)
   !------------------------------------------------------------------
   subroutine boltz_sample(sta)
      type(state_t), intent(inout) :: sta
      call boltz_draw()
      call boltz_realize(sta)
   end subroutine boltz_sample

   !------------------------------------------------------------------
   ! boltz_last_n(k, m) - copy of the last drawn quantum number of carrier
   !                      k's mode m (0 guards: not inited / out of range)
   !------------------------------------------------------------------
   pure integer function boltz_last_n(k, m)
      integer, intent(in) :: k, m
      boltz_last_n = 0
      if (inited .and. k >= 1 .and. k <= size(carry)) then
         if (m >= 1 .and. m <= size(carry(k)%n_last)) boltz_last_n = carry(k)%n_last(m)
      end if
   end function boltz_last_n

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; abort via diag_fatal)
   subroutine stop_boltx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      call diag_fatal(diag_cat_sampler, 'samp_boltzmann.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_boltx

   subroutine stop_boltz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      call diag_fatal(diag_cat_sampler, 'samp_boltzmann.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_boltz

   subroutine stop_bolti(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      call diag_fatal(diag_cat_sampler, 'samp_boltzmann.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_bolti

   subroutine stop_bolti2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      call diag_fatal(diag_cat_sampler, 'samp_boltzmann.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_bolti2

   !------------------------------------------------------------------
   ! member_register - the fixed-name weave entry (prebuild binds it into
   !                   registry_gen_all; the retired reg file's _all body)
   !------------------------------------------------------------------
   subroutine member_register()
      call samp_reg('boltzmann', boltz_draw, boltz_realize)
   end subroutine member_register

   !------------------------------------------------------------------
   ! member_init - the fixed-name init seam (registry_gen_init calls it;
   !               the retired reg file's _init body moved verbatim)
   !------------------------------------------------------------------
   subroutine member_init()
      type(boltz_params_t) :: p          ! member-parameter receiver (field defaults)
      integer :: i, k, wc, code, n_carry, n_word
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('boltzmann')
      n_word = count(reactants%dist_scheme == code)
      if (n_word == 0) return
      ! build tables only for polyatomic carriers; a monoatomic carrier rides
      ! the incident channel alone (BOLTZMANN on it = "no internal sampling")
      n_carry = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code .and. list_atoms%frag(i)%nat > 1) then
            n_carry = n_carry + 1
         end if
      end do
      allocate (p%t_vib(n_word))
      p%t_vib = 0.0d0                     ! absent bracket name keeps the default;
                                           ! boltz_init's own validation rejects it
      wc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         wc = wc + 1
         call par_real('BOLTZMANN', i, 'T_VIB', p%t_vib(wc))
      end do
      allocate (p%tbl(n_carry))
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code .and. list_atoms%frag(i)%nat > 1) then
            k = k + 1
            call spectrum_frag_w(i, p%tbl(k)%w, p%tbl(k)%c)
         end if
      end do
      call boltz_init(p)
   end subroutine member_init

end module samp_boltzmann
