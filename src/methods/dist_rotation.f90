!=====================================================================
! dist_rotation.f90 - distribution member: fragment rotation realization
!   at a GIVEN J - random orientation about the COM plus
!   |L| = sqrt(J(J+1))*hbar isotropic, written as rigid-body momenta
!   p_i = m_i*(omega x r_i)
! Design:
!   This member takes J as its parameter; dist_j draws J itself, and
!   wiring a drawn J between them is assembly-phase business. A linear
!   carrier (smallest inertia eigenvalue ~0) gets L perpendicular to
!   the molecular axis; translation is the incident channel's business.
!   Units: q [Ang]; p [amu*Ang/(10 fs)]; |L|/hbar_code internal; omega [rad/(10 fs)]; inertia [amu*Ang^2].
!=====================================================================
module dist_rotation
   use consts,        only: hbar_code
   use rng,           only: rng_u, rng_gauss
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use sampler,       only: samp_code
   use geometry,      only: rot_euler
   use linalg,        only: eig_sym
   implicit none
   private
   public :: rotation_params_t, rotation_init, rotation_sample

   real(8), parameter :: lin_tol = 1.0d-8   ! linearity eigenvalue ratio tolerance [-]

   ! Member-owned parameters
   type :: rotation_params_t
      integer :: j_rot = 0        ! rotational quantum number J [-] (the GIVEN J - the
                                  ! draw of J is dist_j's business)
   end type

   ! Derived member state (module-private; rebuilt by every rotation_init)
   integer :: frag_served = 0               ! bound roster fragment index [index]
   integer :: nat_c = 0                     ! carrier atom count [count]
   integer, allocatable :: list_c(:)        ! carrier atom index list [global atom index]
   real(8), allocatable :: qz_com(:)        ! COM-shifted seed coordinates [Angstrom]
   real(8) :: al_mom = 0.0d0                ! angular-momentum magnitude |L| [hbar_code unit]
   logical :: inited = .false.              ! init completed [flag]
contains
   !------------------------------------------------------------------
   ! rotation_init(p) - validate the parameter, bind the carrier, seed the
   !                   COM-shifted geometry and the angular-momentum magnitude
   subroutine rotation_init(p)
      type(rotation_params_t), intent(in) :: p ! member parameters (the given J)
      real(8) :: com(3), r(3), wt, trace_i
      integer :: i, j, code, n_carry, i_carry

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%j_rot < 0) then
         call stop_roti('rotation_init', 'rotational quantum number j_rot must be >= 0, got', p%j_rot)
      end if

      ! 2. fragment binding (assembly-order and carrier guards first)
      if (.not. allocated(roster%frag)) then
         call stop_rotx('rotation_init', 'no roster fragment table - the roster must be assembled '// &
            'before member init (assembly-order error)')
      end if
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_rotx('rotation_init', 'no per-fragment scheme array (reactants%dist_scheme '// &
            'unallocated) - member init cannot bind its fragment')
      end if
      if (size(reactants%dist_scheme) /= size(roster%frag)) then
         call stop_roti2('rotation_init', 'dist_scheme length', size(reactants%dist_scheme), &
            'does not match the roster fragment count', size(roster%frag))
      end if
      code = samp_code('rotation')
      if (code == 0) then
         call stop_rotx('rotation_init', 'member word "rotation" is not in the scheme-code table '// &
            '(incomplete member library assembly)')
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
         call stop_rotx('rotation_init', 'no roster fragment carries the rotation scheme - '// &
            'this member would never be dispatched')
      end if
      if (n_carry > 1) then
         call stop_roti2('rotation_init', 'two or more fragments carry the rotation scheme, but one '// &
            'member instance carries one parameter set; count =', n_carry, 'first carrier =', i_carry)
      end if
      frag_served = i_carry
      associate (fr => roster%frag(i_carry))
         if (fr%nat < 2) then
            call stop_roti('rotation_init', 'the rotation carrier must carry at least two atoms (a '// &
               'single atom has no rotational degree of freedom); atom count =', fr%nat)
         end if
         nat_c = fr%nat
         if (allocated(list_c)) deallocate (list_c)
         if (allocated(qz_com)) deallocate (qz_com)
         allocate (list_c(nat_c), qz_com(3*nat_c))
         list_c = fr%list
         ! COM shift of the seed geometry + degeneracy guard (the inertia
         ! trace of a coincident-atom geometry vanishes)
         wt = sum(roster%mass(list_c))
         com = 0.0d0
         do i = 1, nat_c
            com = com + roster%mass(list_c(i))*fr%qz(3*i-2:3*i)
         end do
         com = com/wt
         trace_i = 0.0d0
         do i = 1, nat_c
            r = fr%qz(3*i-2:3*i) - com
            qz_com(3*i-2:3*i) = r
            trace_i = trace_i + roster%mass(list_c(i))*dot_product(r, r)
         end do
         if (trace_i <= 0.0d0) then
            call stop_rotz('rotation_init', 'degenerate carrier geometry: moment of inertia '// &
               '[amu*Ang^2] (trace/2) =', trace_i/2.0d0)
         end if
      end associate
      al_mom = sqrt(dble(p%j_rot*(p%j_rot + 1)))*hbar_code
      inited = .true.
   end subroutine rotation_init

   !------------------------------------------------------------------
   ! rotation_sample(sta) - one realization: seed, orient, build omega, write p
   subroutine rotation_sample(sta)
      type(state_t), intent(inout) :: sta  ! initial state (the carrier fragment's q/p
                                           ! components are overwritten)
      real(8) :: tens(3,3), ev(3), vecs(3,3), u(3), l_vec(3), om(3), r(3), com(3), wt, un
      integer :: i, j, k, i3

      if (.not. inited) then
         call stop_rotx('rotation_sample', 'sample called before rotation_init - the member '// &
            'carries no parameter set')
      end if
      if (size(sta%q) < 3*maxval(list_c) .or. size(sta%p) < 3*maxval(list_c)) then
         call stop_rotx('rotation_sample', 'state too small for the carrier fragment atoms - '// &
            'state_create must follow the roster')
      end if

      ! 1. seed the carrier (COM-shifted geometry, zero momenta)
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         sta%q(i3+1:i3+3) = qz_com(3*i-2:3*i)
         sta%p(i3+1:i3+3) = 0.0d0
      end do

      ! 2. random orientation about the fragment COM (geometry instrument)
      call rot_euler(sta, roster%mass, list_c)

      ! 3. inertia tensor at the written orientation, principal axes
      com = 0.0d0
      wt = sum(roster%mass(list_c))
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         com = com + roster%mass(list_c(i))*sta%q(i3+1:i3+3)
      end do
      com = com/wt
      tens = 0.0d0
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         r = sta%q(i3+1:i3+3) - com
         tens(1,1) = tens(1,1) + roster%mass(list_c(i))*(r(2)**2 + r(3)**2)
         tens(2,2) = tens(2,2) + roster%mass(list_c(i))*(r(1)**2 + r(3)**2)
         tens(3,3) = tens(3,3) + roster%mass(list_c(i))*(r(1)**2 + r(2)**2)
         tens(1,2) = tens(1,2) - roster%mass(list_c(i))*r(1)*r(2)
         tens(1,3) = tens(1,3) - roster%mass(list_c(i))*r(1)*r(3)
         tens(2,3) = tens(2,3) - roster%mass(list_c(i))*r(2)*r(3)
      end do
      tens(2,1) = tens(1,2)
      tens(3,1) = tens(1,3)
      tens(3,2) = tens(2,3)
      vecs = tens
      call eig_sym(vecs, ev)                ! eigenvalues ascending, vectors BY COLUMN

      ! 4. L: |L| along an isotropic direction; the linear carrier projects
      !    the direction perpendicular to the molecular axis (eigenvector 1)
      u = [ rng_gauss(), rng_gauss(), rng_gauss() ]
      un = norm2(u)
      if (un < 1.0d-300) u = [ 1.0d0, 0.0d0, 0.0d0 ]   ! degenerate draw guard (idle in practice)
      u = u/un
      om = 0.0d0
      if (ev(1) <= lin_tol*ev(3)) then
         ! linear carrier: strip the axis component, skip the axis eigenvalue
         u = u - dot_product(u, vecs(:,1))*vecs(:,1)
         un = norm2(u)
         if (un < 1.0d-12) then
            u = vecs(:,2)    ! the draw landed on the axis: use a perpendicular axis (already
         else                ! unit - the stale tiny un must NOT renormalize it)
            u = u/un
         end if
         l_vec = al_mom*u
         do k = 2, 3
            om = om + (dot_product(l_vec, vecs(:,k))/ev(k))*vecs(:,k)
         end do
      else
         l_vec = al_mom*u
         do k = 1, 3
            om = om + (dot_product(l_vec, vecs(:,k))/ev(k))*vecs(:,k)
         end do
      end if

      ! 5. rigid-body momenta (total momentum zero: sum m*(omega x r) with
      !    sum m*r = 0 in the COM frame)
      do i = 1, nat_c
         i3 = 3*(list_c(i) - 1)
         r = sta%q(i3+1:i3+3) - com
         sta%p(i3+1) = roster%mass(list_c(i))*(om(2)*r(3) - om(3)*r(2))
         sta%p(i3+2) = roster%mass(list_c(i))*(om(3)*r(1) - om(1)*r(3))
         sta%p(i3+3) = roster%mass(list_c(i))*(om(1)*r(2) - om(2)*r(1))
      end do
   end subroutine rotation_sample

   !------------------------------------------------------------------
   ! private loud-stop helpers (name the caller; value variants; STOP 1)
   subroutine stop_rotx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_rotx

   subroutine stop_rotz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_rotz

   subroutine stop_roti(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_roti

   subroutine stop_roti2(where, msg, ival1, tail, ival2)
      character(len=*), intent(in) :: where, msg, tail
      integer, intent(in) :: ival1, ival2
      write (0, '(a,a,i0,a,i0)') trim(where)//': ', trim(msg), ival1, ' '//trim(tail), ival2
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_roti2

end module dist_rotation
