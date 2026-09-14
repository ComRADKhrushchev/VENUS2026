!=====================================================================
! dist_normalmode.f90 - distribution member (methods/ library): the
!                     microcanonical normal-mode sampling of a fragment
!                     - a fixed total vibrational energy e_vib is split
!                     over a GIVEN orthonormal mass-weighted mode table
!                     (frequencies + eigenvector columns) by a uniform
!                     point on the unit sphere in mode space
!                     (E_i = e_vib*u_i^2, u = g/|g| with g gaussian -
!                     the sphere-uniform generator), each mode then
!                     placed at a uniform phase on its own ellipse:
!                     Q_i = A_i*cos(phi_i), P_i = -w_i*A_i*sin(phi_i)
!                     with A_i = sqrt(2*E_i)/w_i. Coordinates and
!                     momenta are written about the carrier's
!                     center-of-mass origin (the family posture):
!                     q = qz_com + (C*Q)/sqrt(m), p = sqrt(m)*(C*P).
! Contract:
!   normalmode_init(p) - validate the member parameters and bind the
!     served fragment. Loud-stops name the offending field and value:
!       e_vib <= 0 (total vibrational energy); the binding battery (no
!       roster fragment table / no scheme array / length mismatch / the
!       member word absent from the scheme-code table / zero or multiple
!       carriers); an unallocated mode table (either half); a frequency
!       table and an eigenvector table of different column counts; an
!       eigenvector table whose row count is not 3*(carrier atom count);
!       an empty mode table (no vibrational support to distribute over);
!       a mode count that leaves no room for the three translations
!       (n_mode > 3*nat - 3); a non-positive frequency; a table that is
!       not orthonormal to 1e-8; a table column not orthogonal to the
!       mass-weighted translations to 1e-8 (such a column would give the
!       fragment a net momentum).
!     Precompute: the COM-shifted seed coordinates, the replicated
!       sqrt(mass) vector, and e_tot = e_vib*e_conv.
!   normalmode_sample(sta) - one realization: gaussian vector g (one per
!     mode), u = g/|g| (degenerate-draw guard idles in practice), the
!     energy split, per-mode amplitudes and uniform phases, then the
!     back-transformation onto the carrier's q/p components.
!     Identities exact by construction: sum_i 0.5*(P_i^2 + w_i^2*Q_i^2)
!     = e_tot; total momentum = sqrt(m)*(C*P) summed over atoms = 0
!     (translation-orthogonality); the COM of the written q is the
!     origin; the written state lies in span(C).
!     The frequencies and the mode table ARRIVE as assembly data - the
!     derivation of a physical table from a Hessian is assembly-phase
!     mechanism (the hessian instrument exists for it); this member
!     consumes any table meeting the contract above (its criteria
!     program drives a synthetic table).
! Rulings carried in the member's design (variant inventory against the
!   old defect register, semantics only, no citations): the LOCAL-MODE
!   variant is NOT implemented (its parameters had no source - the old
!   path produced NaN; deleted per the pre-deletion rule); the collinear
!   FIXED-ENERGY trigger variants are NOT implemented (the collinear
!   pseudo-mode overflowed its step bound and stopped the run; deleted
!   per the pre-deletion rule) - the legal microcanonical body is
!   retained; the vibrational angular-momentum stripping of the old
!   initializer is NOT carried (it distorted the microcanonical energy
!   split - the sphere split here carries the full energy; angular
!   momentum is the rotation member's channel, not subtracted here); the
!   orientation randomization after assembly is NOT carried (the mode
!   table is expressed in the roster equilibrium frame; orientation is
!   the rotation member's business); the energy-scale resampling loop of
!   the old initializer is NOT carried (the microcanonical total is the
!   parameter and is placed exactly once); the orthant-restriction face
!   (mode signs confined to one orthant to direct the flow) is NOT
!   carried (it biases the microcanonical measure and its trigger
!   semantics were unclear - dubious-variant deletion).
! Units: e_vib [kcal/mol] (input, converted with e_conv); E [internal]
!   (= kcal/mol * e_conv); w_mode [rad/(10 fs)] (so 0.5*hbar*w is
!   internal); c_mode columns orthonormal in mass-weighted coordinates
!   (components [1/sqrt(amu)]-scaled cartesian directions); q
!   [Angstrom]; p [amu*Ang/(10 fs)]; masses [amu].
!=====================================================================
module dist_normalmode
   use consts,        only: e_conv, two_pi
   use rng,           only: rng_u, rng_gauss
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use sampler,       only: samp_code
   implicit none
   private
   public :: normalmode_params_t, normalmode_init, normalmode_sample

   real(8), parameter :: orth_tol = 1.0d-8    ! C^T*C - I tolerance [-]
   real(8), parameter :: trans_tol = 1.0d-8   ! translation-orthogonality tolerance [-]

   ! Member-owned parameters (the mode table arrives as assembly data)
   type :: normalmode_params_t
      real(8) :: e_vib = 10.0d0               ! total vibrational energy [kcal/mol]
      real(8), allocatable :: w_mode(:)       ! mode frequencies [rad/(10 fs)]
      real(8), allocatable :: c_mode(:,:)     ! (3*nat x n_mode) orthonormal mass-weighted columns
   end type

   ! Derived member state (module-private; rebuilt by every normalmode_init)
   integer :: frag_served = 0                 ! bound roster fragment index [index]
   integer :: nat_c = 0                       ! carrier atom count [count]
   integer :: n_mode_c = 0                    ! mode count [count]
   integer, allocatable :: list_c(:)          ! carrier atom index list [global atom index]
   real(8), allocatable :: qz_com(:)          ! COM-shifted seed coordinates [Angstrom]
   real(8), allocatable :: sqm(:)             ! sqrt(mass) replicated per component [sqrt(amu)]
   real(8), allocatable :: w_c(:)             ! validated frequencies [rad/(10 fs)]
   real(8), allocatable :: c_c(:,:)           ! validated eigenvector table
   real(8) :: e_tot = 0.0d0                   ! total vibrational energy [internal]
   logical :: inited = .false.                ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! normalmode_init(p) - validate the parameters and mode table, bind
   !                      the carrier, precompute the seed and constants
   subroutine normalmode_init(p)
      type(normalmode_params_t), intent(in) :: p  ! member parameters
      real(8) :: com(3), wt, dev, tv(3), tscale
      integer :: i, k, d, code, n_carry, i_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%e_vib <= 0.0d0) then
         call stop_nmz('normalmode_init', 'total vibrational energy e_vib [kcal/mol] must be '// &
            'positive, got', p%e_vib)
      end if

      ! 2. fragment binding (assembly-order and carrier guards first)
      if (.not. allocated(roster%frag)) then
         call stop_nmx('normalmode_init', 'no roster fragment table - the roster must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_nmx('normalmode_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_nmi2('normalmode_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('normalmode')
      if (code == 0) then
         call stop_nmx('normalmode_init', 'member word "normalmode" is not in the scheme-code '// &
            'table (incomplete member library assembly)')
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
         call stop_nmx('normalmode_init', 'no roster fragment carries the normalmode scheme - '// &
            'this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_nmi2('normalmode_init', 'two or more fragments carry the normalmode scheme, '// &
            'but one member instance carries one parameter set; count =', n_carry, &
            'first carrier =', i_carry)
      end if
      frag_served = i_carry
      nat_c = roster%frag(i_carry)%nat

      ! 3. mode-table validation (shape, contents, orthonormality)
      if (.not. allocated(p%w_mode) .or. .not. allocated(p%c_mode)) then
         call stop_nmx('normalmode_init', 'mode table unallocated - the frequencies and the '// &
            'mass-weighted eigenvector columns are assembly data the member requires')
      end if
      if (size(p%c_mode, 2) /= size(p%w_mode)) then
         call stop_nmi2('normalmode_init', 'mode table mismatch: eigenvector column count =', &
            size(p%c_mode, 2), 'but frequency count =', size(p%w_mode))
      end if
      if (size(p%c_mode, 1) /= 3*nat_c) then
         call stop_nmi2('normalmode_init', 'mode table shape mismatch: eigenvector row count =', &
            size(p%c_mode, 1), 'but the carrier carries 3*nat =', 3*nat_c)
      end if
      n_mode_c = size(p%w_mode)
      if (n_mode_c < 1) then
         call stop_nmi('normalmode_init', 'empty mode table - no vibrational support to '// &
            'distribute over; mode count =', n_mode_c)
      end if
      if (n_mode_c > 3*nat_c - 3) then
         call stop_nmi2('normalmode_init', 'mode count =', n_mode_c, &
            'leaves no room for the three translations (max 3*nat - 3 =', 3*nat_c - 3)
      end if
      do k = 1, n_mode_c
         if (p%w_mode(k) <= 0.0d0) then
            call stop_nmz('normalmode_init', 'mode frequency w_mode [rad/(10 fs)] must be '// &
            'positive; entry =', p%w_mode(k))
         end if
      end do
      do k = 1, n_mode_c
         do i = 1, n_mode_c
            dev = dot_product(p%c_mode(:, k), p%c_mode(:, i)) - merge(1.0d0, 0.0d0, k == i)
            if (abs(dev) > orth_tol) then
               call stop_nmz('normalmode_init', 'mode table not orthonormal (C^T*C - I entry, '// &
                  'tolerance 1e-8) =', dev)
            end if
         end do
      end do
      tscale = sqrt(sum(roster%mass(roster%frag(i_carry)%list)))
      do k = 1, n_mode_c
         tv = 0.0d0
         do i = 1, nat_c
            tv = tv + sqrt(roster%mass(roster%frag(i_carry)%list(i)))*p%c_mode(3*i-2:3*i, k)
         end do
         if (norm2(tv) > trans_tol*tscale) then
            call stop_nmz('normalmode_init', 'mode column not orthogonal to the mass-weighted '// &
               'translations (residual norm, tolerance 1e-8*sqrt(sum m)) =', norm2(tv))
         end if
      end do

      ! 4. precompute (COM-shifted seed, sqrt(mass) vector, constants)
      associate (fr => roster%frag(i_carry))
         if (allocated(list_c)) deallocate (list_c)
         if (allocated(qz_com)) deallocate (qz_com)
         if (allocated(sqm)) deallocate (sqm)
         if (allocated(w_c)) deallocate (w_c)
         if (allocated(c_c)) deallocate (c_c)
         allocate (list_c(nat_c), qz_com(3*nat_c), sqm(3*nat_c), w_c(n_mode_c), c_c(3*nat_c, n_mode_c))
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
      w_c = p%w_mode
      c_c = p%c_mode
      e_tot = p%e_vib*e_conv
      inited = .true.
   end subroutine normalmode_init

   !------------------------------------------------------------------
   ! normalmode_sample(sta) - one microcanonical realization
   subroutine normalmode_sample(sta)
      type(state_t), intent(inout) :: sta  ! initial state (the carrier's q/p
                                           ! components are overwritten)
      real(8) :: g_vec(n_mode_c), un, e_i, amp, ph
      real(8) :: qmode(n_mode_c), pmode(n_mode_c), x_vec(3*nat_c), y_vec(3*nat_c)
      integer :: k, i, i3

      if (.not. inited) then
         call stop_nmx('normalmode_sample', 'sample called before normalmode_init - the member '// &
            'carries no parameter set')
      end if
      if (size(sta%q) < 3*maxval(list_c) .or. size(sta%p) < 3*maxval(list_c)) then
         call stop_nmx('normalmode_sample', 'state too small for the carrier fragment atoms - '// &
            'state_create must follow the roster')
      end if

      ! 1. uniform point on the unit sphere in mode space (gaussian/radial)
      do k = 1, n_mode_c
         g_vec(k) = rng_gauss()
      end do
      un = norm2(g_vec)
      if (un < 1.0d-300) then
         g_vec = 0.0d0                     ! degenerate-draw guard (idle in practice)
         g_vec(1) = 1.0d0
         un = 1.0d0
      end if

      ! 2. energy split + uniform phase per mode
      do k = 1, n_mode_c
         e_i = e_tot*(g_vec(k)/un)**2
         amp = sqrt(2.0d0*e_i)/w_c(k)
         ph = two_pi*rng_u()
         qmode(k) = amp*cos(ph)
         pmode(k) = -w_c(k)*amp*sin(ph)
      end do

      ! 3. back-transformation about the COM origin
      x_vec = matmul(c_c, qmode)
      y_vec = matmul(c_c, pmode)
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         sta%q(i3+1:i3+3) = qz_com(3*i-2:3*i) + x_vec(3*i-2:3*i)/sqm(3*i-2:3*i)
         sta%p(i3+1:i3+3) = sqm(3*i-2:3*i)*y_vec(3*i-2:3*i)
      end do
   end subroutine normalmode_sample

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_nmx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_nmx

   subroutine stop_nmz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_nmz

   subroutine stop_nmi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_nmi

   subroutine stop_nmi2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_nmi2

end module dist_normalmode
