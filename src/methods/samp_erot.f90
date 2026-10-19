!   The member carries its own registration (member_register /
!   member_init - the fixed-name entries of the 2026-10-14 weave batch).
!=====================================================================
! samp_erot.f90 - distribution member: NONLINEAR-fragment thermal
!   rotation - three independent principal-axis classical thermal
!   rotational energies E_k = 0.5*kb*T_ROT*Exp(1) (ANT-manual A1 law:
!   one exponential draw per inertia principal axis), a random
!   orientation, and the rigid-body realization p_i = m_i*(omega x r_i)
!   with omega_k = sign_k*sqrt(2*E_k/I_k) along the principal axes
! Design:
!   Nonlinear carriers only (the linear rotor's thermal channel is the
!   J word + ROTATION realization - the J ladder is quantized there;
!   init rejects a linear carrier by the inertia eigenvalue ratio).
!   The word self-orients (the ROTATION member's design), so the
!   sampler's RANDOM_ORIENT eligibility mask excludes it (no double
!   rotation). Expected total rotational energy 3/2*kb*T_ROT. Draw
!   order per carrier: orientation matrix (three Euler uniforms), then
!   the three energies (one uniform each), then the three axis signs
!   (one uniform each) - nine uniforms per carrier; realize draws no
!   RNG and reconstructs everything at the WRITTEN geometry (the
!   orientation-invariance of the inertia eigenvalues is what makes
!   the drawn per-axis energies geometric-portable).
!   Units: t_rot [K]; E_k [amu*Ang^2/fs^2] (kb_act face); q [Ang];
!   p [amu*Ang/fs]; omega [rad/fs]; inertia [amu*Ang^2].
!=====================================================================
module samp_erot
   use consts,        only: kb_act
   use rng,           only: rng_u
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use sampler,       only: samp_code, samp_reg
   use geometry,      only: rot_rand_rmat, rot_apply
   use linalg,        only: eig_sym
   use input, only: par_real
   use diag,        only: diag_fatal, diag_cat_sampler
   implicit none
   private
   public :: erot_params_t, erot_init, erot_draw, erot_realize, &
             erot_sample, member_register, member_init

   real(8), parameter :: lin_tol = 1.0d-8   ! linearity eigenvalue ratio tolerance [-]

   ! Member-owned parameters, one row per carrier (the 2026-10-13 bracket
   ! ruling realized per fragment)
   type :: erot_params_t
      real(8), allocatable :: t_rot(:)  ! rotational temperature [K] (must be > 0;
                                        ! absent bracket name keeps 0 - the init
                                        ! guard rejects it loudly)
   end type erot_params_t

   ! Derived member state (module-private; rebuilt by every erot_init) -
   ! one row per carrier
   type :: carry_t
      integer :: nat = 0                     ! carrier atom count [count]
      integer, allocatable :: list(:)        ! carrier atom index list [global atom index]
      real(8), allocatable :: qz_com(:)      ! COM-shifted seed coordinates [Angstrom]
      real(8) :: t_rot = 0.0d0               ! the carrier's temperature [K]
      real(8) :: rmat_s(3,3) = 0.0d0         ! the drawn orientation matrix (post-draw record) [-]
      real(8) :: e_s(3) = 0.0d0              ! the drawn per-axis energies [amu*Ang^2/fs^2]
      real(8) :: sgn_s(3) = 0.0d0            ! the drawn per-axis omega signs [-]
   end type carry_t
   type(carry_t), allocatable :: carry(:)    ! carrier rows (list_atoms order)
   logical :: inited = .false.               ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! erot_init(p) - validate the parameters, bind the carrier, seed the
   !                COM-shifted geometry and reject linear carriers
   subroutine erot_init(p)
      type(erot_params_t), intent(in) :: p   ! member parameters (the temperature,
                                              ! one row per carrier)
      real(8) :: com(3), r(3), wt, trace_i, tens(3,3), ev(3)
      integer :: i, j, k, code, n_carry

      ! 1. parameter validation, per carrier row (fail-loud, naming the field)
      if (.not. allocated(p%t_rot)) then
         call stop_erx('erot_init', 'no parameter rows (t_rot unallocated) - the reg '// &
            'must stage one row per erot carrier')
      end if
      do k = 1, size(p%t_rot)
         if (p%t_rot(k) <= 0.0d0) then
            call stop_erz('erot_init', 'rotational temperature t_rot [K] must be positive, got', &
               p%t_rot(k))
         end if
      end do

      ! 2. fragment binding (assembly-order and carrier guards first)
      if (.not. allocated(list_atoms%frag)) then
         call stop_erx('erot_init', 'no list_atoms fragment table - the atom list must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_erx('erot_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(list_atoms%frag)) then
         call stop_eri2('erot_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the atom list fragment count', size(list_atoms%frag))
      end if
      code = samp_code('erot')
      if (code == 0) then
         call stop_erx('erot_init', 'member word "erot" is not in the scheme-code table '// &
            '(incomplete member library assembly)')
      end if
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) then
         call stop_erx('erot_init', 'no list_atoms fragment carries the erot scheme - '// &
            'this member would never be dispatched')
      end if
      if (size(p%t_rot) /= n_carry) then
         call stop_eri2('erot_init', 'parameter row count', size(p%t_rot), &
            'does not match the erot carrier count', n_carry)
      end if
      if (allocated(carry)) deallocate (carry)
      allocate (carry(n_carry))
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         associate (fr => list_atoms%frag(i), cr => carry(k))
            if (fr%nat < 2) then
               call stop_eri('erot_init', 'the erot carrier must carry at least two atoms (a '// &
                  'single atom has no rotational degree of freedom); atom count =', fr%nat)
            end if
            cr%nat = fr%nat
            cr%t_rot = p%t_rot(k)
            allocate (cr%list(cr%nat), cr%qz_com(3*cr%nat))
            cr%list = fr%list
            ! COM shift of the seed geometry + degeneracy guard, then the
            ! linearity rejection (the eigenvalue ratio at the seed geometry)
            wt = sum(list_atoms%mass(cr%list))
            com = 0.0d0
            do j = 1, cr%nat
               com = com + list_atoms%mass(cr%list(j))*fr%qz(3*j-2:3*j)
            end do
            com = com/wt
            tens = 0.0d0
            trace_i = 0.0d0
            do j = 1, cr%nat
               r = fr%qz(3*j-2:3*j) - com
               cr%qz_com(3*j-2:3*j) = r
               trace_i = trace_i + list_atoms%mass(cr%list(j))*dot_product(r, r)
               tens(1,1) = tens(1,1) + list_atoms%mass(cr%list(j))*(r(2)**2 + r(3)**2)
               tens(2,2) = tens(2,2) + list_atoms%mass(cr%list(j))*(r(1)**2 + r(3)**2)
               tens(3,3) = tens(3,3) + list_atoms%mass(cr%list(j))*(r(1)**2 + r(2)**2)
               tens(1,2) = tens(1,2) - list_atoms%mass(cr%list(j))*r(1)*r(2)
               tens(1,3) = tens(1,3) - list_atoms%mass(cr%list(j))*r(1)*r(3)
               tens(2,3) = tens(2,3) - list_atoms%mass(cr%list(j))*r(2)*r(3)
            end do
            if (trace_i <= 0.0d0) then
               call stop_erz('erot_init', 'degenerate carrier geometry: moment of inertia '// &
                  '[amu*Ang^2] (trace/2) =', trace_i/2.0d0)
            end if
            tens(2,1) = tens(1,2)
            tens(3,1) = tens(1,3)
            tens(3,2) = tens(2,3)
            call eig_sym(tens, ev)            ! eigenvalues ascending (vectors unused here)
            if (ev(1) <= lin_tol*ev(3)) then
               call stop_erz('erot_init', 'the erot carrier is LINEAR (smallest inertia '// &
                  'eigenvalue ratio below tolerance - the linear-rotor thermal channel is '// &
                  'the J word + ROTATION realization); ratio =', ev(1)/ev(3))
            end if
         end associate
      end do
      inited = .true.
   end subroutine erot_init

   !------------------------------------------------------------------
   ! erot_draw() - per carrier (the legislated stream order): the
   !               orientation matrix first (three Euler uniforms), then
   !               the three principal-axis energies (one uniform each),
   !               then the three axis signs (one uniform each); no
   !               state access
   !------------------------------------------------------------------
   subroutine erot_draw()
      integer :: kc, k
      if (.not. inited) then
         call stop_erx('erot_draw', 'draw called before erot_init - the member '// &
            'carries no parameter set')
      end if
      do kc = 1, size(carry)
         call rot_rand_rmat(carry(kc)%rmat_s)
         do k = 1, 3
            ! E_k = 0.5*kb*T*Exp(1): the -log(1-u) form is finite for u in [0,1)
            carry(kc)%e_s(k) = -0.5d0*kb_act*carry(kc)%t_rot*log(1.0d0 - rng_u())
            if (rng_u() < 0.5d0) then
               carry(kc)%sgn_s(k) = -1.0d0
            else
               carry(kc)%sgn_s(k) = 1.0d0
            end if
         end do
      end do
   end subroutine erot_draw

   !------------------------------------------------------------------
   ! erot_realize(sta) - per carrier (list_atoms order): seed, apply the
   !                     drawn orientation, rebuild the principal axes at
   !                     the written geometry, write the rigid-body
   !                     momenta from the drawn energies - no RNG draws
   !------------------------------------------------------------------
   subroutine erot_realize(sta)
      type(state_t), intent(inout) :: sta  ! initial state (every carrier fragment's
                                           ! q/p components are overwritten, list_atoms order)
      real(8) :: tens(3,3), ev(3), vecs(3,3), om(3), r(3), com(3), wt
      integer :: i, k, i3, kc

      if (.not. inited) then
         call stop_erx('erot_realize', 'realize called before erot_init - the member '// &
            'carries no parameter set')
      end if

      do kc = 1, size(carry)
      associate (nat_c => carry(kc)%nat, list_c => carry(kc)%list, qz_com => carry(kc)%qz_com)
      if (size(sta%q) < 3*maxval(list_c) .or. size(sta%p) < 3*maxval(list_c)) then
         call stop_erx('erot_realize', 'state too small for the carrier fragment atoms - '// &
            'state_create must follow the atom list')
      end if

      ! 1. seed the carrier (COM-shifted geometry, zero momenta)
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         sta%q(i3+1:i3+3) = qz_com(3*i-2:3*i)
         sta%p(i3+1:i3+3) = 0.0d0
      end do

      ! 2. apply the drawn orientation about the fragment COM
      call rot_apply(sta, list_atoms%mass, list_c, carry(kc)%rmat_s)

      ! 3. inertia tensor at the written orientation, principal axes
      com = 0.0d0
      wt = sum(list_atoms%mass(list_c))
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         com = com + list_atoms%mass(list_c(i))*sta%q(i3+1:i3+3)
      end do
      com = com/wt
      tens = 0.0d0
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         r = sta%q(i3+1:i3+3) - com
         tens(1,1) = tens(1,1) + list_atoms%mass(list_c(i))*(r(2)**2 + r(3)**2)
         tens(2,2) = tens(2,2) + list_atoms%mass(list_c(i))*(r(1)**2 + r(3)**2)
         tens(3,3) = tens(3,3) + list_atoms%mass(list_c(i))*(r(1)**2 + r(2)**2)
         tens(1,2) = tens(1,2) - list_atoms%mass(list_c(i))*r(1)*r(2)
         tens(1,3) = tens(1,3) - list_atoms%mass(list_c(i))*r(1)*r(3)
         tens(2,3) = tens(2,3) - list_atoms%mass(list_c(i))*r(2)*r(3)
      end do
      tens(2,1) = tens(1,2)
      tens(3,1) = tens(1,3)
      tens(3,2) = tens(2,3)
      vecs = tens
      call eig_sym(vecs, ev)                ! eigenvalues ascending, vectors BY COLUMN

      ! 4. omega along the principal axes from the DRAWN energies (the
      !    eigenvalues are orientation-invariant, so the seed-geometry
      !    linearity rejection carries; the defensive guard keeps a
      !    degenerate diagonalization honest)
      om = 0.0d0
      do k = 1, 3
         if (ev(k) <= 0.0d0) then
            call stop_erz('erot_realize', 'nonpositive inertia eigenvalue at the written '// &
               'geometry [amu*Ang^2] =', ev(k))
         end if
         om = om + carry(kc)%sgn_s(k)*sqrt(2.0d0*carry(kc)%e_s(k)/ev(k))*vecs(:,k)
      end do

      ! 5. rigid-body momenta (total momentum zero: sum m*(omega x r) with
      !    sum m*r = 0 in the COM frame)
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         r = sta%q(i3+1:i3+3) - com
         sta%p(i3+1) = list_atoms%mass(list_c(i))*(om(2)*r(3) - om(3)*r(2))
         sta%p(i3+2) = list_atoms%mass(list_c(i))*(om(3)*r(1) - om(1)*r(3))
         sta%p(i3+3) = list_atoms%mass(list_c(i))*(om(1)*r(2) - om(2)*r(1))
      end do
      end associate
      end do
   end subroutine erot_realize

   !------------------------------------------------------------------
   ! erot_sample(sta) - fused convenience: draw then realize (the
   !                     check-program entry)
   !------------------------------------------------------------------
   subroutine erot_sample(sta)
      type(state_t), intent(inout) :: sta
      call erot_draw()
      call erot_realize(sta)
   end subroutine erot_sample

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; abort via diag_fatal)
   subroutine stop_erx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      call diag_fatal(diag_cat_sampler, 'samp_erot.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_erx

   subroutine stop_erz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      call diag_fatal(diag_cat_sampler, 'samp_erot.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_erz

   subroutine stop_eri(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      call diag_fatal(diag_cat_sampler, 'samp_erot.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_eri

   subroutine stop_eri2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      call diag_fatal(diag_cat_sampler, 'samp_erot.f90', trim(where), &
         'the initial state is not sampled')
   end subroutine stop_eri2

   !------------------------------------------------------------------
   ! member_register - the fixed-name weave entry (prebuild binds it into
   !                   registry_gen_all)
   !------------------------------------------------------------------
   subroutine member_register()
      call samp_reg('erot', erot_draw, erot_realize)
   end subroutine member_register

   !------------------------------------------------------------------
   ! member_init - the fixed-name init seam (registry_gen_init calls it);
   !               T_ROT rides the EROT bracket, absent name keeps 0
   !               (the init guard rejects it loudly)
   !------------------------------------------------------------------
   subroutine member_init()
      type(erot_params_t) :: p         ! member-parameter receiver (field defaults)
      integer :: i, k, code
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('erot')
      if (count(reactants%dist_scheme == code) == 0) return
      allocate (p%t_rot(count(reactants%dist_scheme == code)))
      p%t_rot = 0.0d0                    ! absent bracket name keeps the default
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         call par_real('EROT', i, 'T_ROT', p%t_rot(k))
      end do
      call erot_init(p)
   end subroutine member_init

end module samp_erot
