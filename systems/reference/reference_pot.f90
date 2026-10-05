!=====================================================================
! reference_pot.f90 - the reference potentials: pure-math test material
!                     (an anisotropic harmonic well, a Morse diatomic
!                     pair, a triatomic LEPS with its two-surface
!                     diabatic form) - no program knowledge, no system
!                     knowledge, zero dependencies
! Material status: every routine is system-independent pure mathematics.
!   Positions in angstrom, energies in eV, gradients g = dV/dq in
!   eV/angstrom (the member unit contract of the delivery channel: the
!   eV-to-internal conversion belongs to the force interface alone). The well acts on
!   every atom; the pair/triatomic forms read atoms 1-2 / 1-2-3 and
!   leave the remaining atoms at zero gradient.
! Parameters: module data carrying LITERATURE REFERENCE CONSTANTS as
!   defaults - the H2 spectroscopy set (De = 4.746 eV, re = 1.401 A,
!   a = 1.028 1/A) and the H3 LEPS Sato parameter (delta = 0.164).
!   These are reference values of the analytic forms, not system
!   assets: a run overrides them through the container's parameter
!   channel (the box file of reference_box), never by editing this
!   file. The harmonic defaults are the unit well k = (1, 1, 1) eV/A^2
!   centered at the origin.
! Forms (the analytic contracts every consumer may rely on):
!   harm_vg   V = sum_i 0.5*[kx*(x_i-x0)^2 + ky*(y_i-y0)^2 + kz*(z_i-z0)^2]
!             on every atom i; g = k*(q-x0) component-wise
!   morse_vg  V = De*(em-1)^2 - De, em = exp(-a(r-re)) on the atom 1-2
!             distance r => V(re) = -De, V(inf) = 0 (zero of energy at
!             dissociation; well depth De); dV/dr = 2*De*(em-1)*(-a*em)
!   leps_vg   the CGM formulation on atoms 1-2-3, pairs p = 1:(1,2),
!             2:(2,3), 3:(1,3):
!               Q_p(r) = De/2*[(3+d)e^-2a(r-re) - (2+2d)e^-a(r-re)]   (triplet)
!               A_p(r) = De/2*[(3-d)e^-2a(r-re) + (2-2d)e^-a(r-re)] + De*d (singlet)
!               V = (Q1+Q2+Q3)/2 - sqrt((A1-Q1)^2 + (A2-Q2)^2 + (A3-Q3)^2)/2
!             (the lower London sheet)
!   leps_diab the SAME physics as the two-surface diabatic matrix whose
!             eigenvalues are the two London sheets:
!               H = [ s/2   -rt/2 ]   s = Q1+Q2+Q3
!                   [ -rt/2   s/2 ]   rt = sqrt(sum_p (A_p-Q_p)^2)
!             lambda_min(H) = s/2 - rt/2 = the scalar leps_vg surface;
!             the upper sheet s/2 + rt/2 rides along. Derivatives are
!             analytic: dH_11 = dH_22 = (1/2) sum_p dQ_p dr_p; dH_12 =
!             -(1/2) sum_p (A_p-Q_p)(dA_p-dQ_p) dr_p / rt, chained to
!             Cartesian dofs through the pair units. Delivery is the
!             FULL MATRIX (both triangles of H and of every dh face) -
!             the workspace contract reads both triangles and assumes
!             no symmetry on read.
!   pair_dist / pair_unit - the pair-distance and pair-unit helpers the
!             LEPS chain rule rides on (public: check material reuses
!             them).
!=====================================================================
module reference_pot
   implicit none
   private
   public :: harm_vg, morse_vg, leps_vg, leps_diab
   public :: pair_dist, pair_unit
   public :: harm_set_constant

   ! HARMONIC: independent anisotropic well on every atom
   real(8), save :: h_k(3) = [ 1.0d0, 1.0d0, 1.0d0 ]  ! force constants [eV/A^2] (unit well)
   real(8), save :: h_x0(3) = 0.0d0                    ! well centers [A]

   ! MORSE: diatomic pair potential on atoms 1-2
   real(8), save :: m_de = 4.746d0      ! eV        (H2 literature value)
   real(8), save :: m_re = 1.401d0      ! angstrom  (H2)
   real(8), save :: m_a  = 1.028d0      ! 1/angstrom (H2)

   ! LEPS: triatomic on atoms 1-2-3, classic H3 parameters
   real(8), save :: l_de    = 4.746d0   ! eV        (H3 literature values)
   real(8), save :: l_re    = 1.401d0   ! angstrom
   real(8), save :: l_a     = 1.028d0   ! 1/angstrom
   real(8), save :: l_delta = 0.164d0   ! Sato parameter
contains
   !------------------------------------------------------------------
   ! harm_set_constant(k) - apply one isotropic force constant to every
   !                 axis of the well (the container's parameter channel;
   !                 the centers stay at the origin). The value's validity
   !                 is the caller's gate - the box file parser rejects
   !                 non-positive constants before this runs
   subroutine harm_set_constant(k)
      real(8), intent(in) :: k          ! isotropic force constant [eV/A^2]
      h_k = k
   end subroutine harm_set_constant

   !------------------------------------------------------------------
   ! harm_vg(natom, q, v, g) - the anisotropic well on every atom:
   !                 energy [eV] + gradient [eV/A] (g = dV/dq)
   subroutine harm_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      integer :: i, k
      v = 0.0d0
      g = 0.0d0
      do i = 1, natom
         do k = 1, 3
            v = v + 0.5d0*h_k(k)*(q(3*(i-1)+k) - h_x0(k))**2
            g(3*(i-1)+k) = h_k(k)*(q(3*(i-1)+k) - h_x0(k))
         end do
      end do
   end subroutine harm_vg

   !------------------------------------------------------------------
   ! morse_vg(natom, q, v, g) - the Morse pair on atoms 1-2: energy [eV] +
   !                 gradient [eV/A]; atoms beyond the pair carry zero
   subroutine morse_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      real(8) :: d(3), r, em, dvdr
      g = 0.0d0
      d = q(4:6) - q(1:3)
      r = sqrt(dot_product(d, d))
      em = exp(-m_a*(r - m_re))
      v = m_de*(em - 1.0d0)**2 - m_de          ! V(re) = -De, V(inf) = 0
      ! V = De*(em-1)^2 - De  =>  dV/dr = 2*De*(em-1)*(-a*em)
      dvdr = -2.0d0*m_de*(em - 1.0d0)*m_a*em
      g(1:3) = -dvdr*d/r                        ! atom 1
      g(4:6) =  dvdr*d/r                        ! atom 2
   end subroutine morse_vg

   !------------------------------------------------------------------
   ! leps_vg(natom, q, v, g) - the lower London sheet on atoms 1-2-3:
   !                 energy [eV] + gradient [eV/A]
   subroutine leps_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      real(8) :: r(3), u(3,3)        ! pair distances; u(k,p) = d r_p / d(cart k of atom B-A scheme below)
      real(8) :: qq(3), aa(3), dqq(3), daa(3)
      real(8) :: s2, rt, dvdq(3), dvda(3), dvdr(3)
      integer :: p, k
      integer :: pa(3), pb(3)
      data pa /1, 2, 1/  ! pair atoms (A<B): (1,2),(2,3),(1,3)
      data pb /2, 3, 3/

      do p = 1, 3
         call pair_dist(q(3*pa(p)-2:3*pa(p)), q(3*pb(p)-2:3*pb(p)), r(p))
         call leps_branch(r(p), qq(p), aa(p), dqq(p), daa(p))
      end do

      s2 = sum((aa - qq)**2)
      rt = sqrt(s2)
      v = 0.5d0*sum(qq) - 0.5d0*rt
      ! d V/d Q_p = 1/2 + (A_p-Q_p)/(2 rt);  d V/d A_p = -(A_p-Q_p)/(2 rt)
      do p = 1, 3
         dvdq(p) = 0.5d0 + (aa(p) - qq(p))/(2.0d0*rt)
         dvda(p) = -(aa(p) - qq(p))/(2.0d0*rt)
         dvdr(p) = dvdq(p)*dqq(p) + dvda(p)*daa(p)
      end do

      g = 0.0d0
      do p = 1, 3
         ! r_p = |x_B - x_A|: dr/dx_A = -u, dr/dx_B = +u
         call pair_unit(q(3*pa(p)-2:3*pa(p)), q(3*pb(p)-2:3*pb(p)), u(:,p))
         do k = 1, 3
            g(3*(pa(p)-1)+k) = g(3*(pa(p)-1)+k) - dvdr(p)*u(k,p)
            g(3*(pb(p)-1)+k) = g(3*(pb(p)-1)+k) + dvdr(p)*u(k,p)
         end do
      end do
   end subroutine leps_vg

   !------------------------------------------------------------------
   ! leps_diab(natom, q, h_ij, dh_ij) - the two-surface diabatic form of
   !                 the same LEPS: the FULL symmetric H [eV] and its
   !                 Cartesian derivative faces dh_ij(3*natom, 2, 2)
   !                 [eV/A], both triangles delivered (the workspace
   !                 contract reads both and assumes no symmetry)
   subroutine leps_diab(natom, q, h_ij, dh_ij)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: h_ij(2, 2)
      real(8), intent(out) :: dh_ij(3*natom, 2, 2)
      real(8) :: r(3), u(3, 3)       ! pair units u(k,p) - columns contiguous at the call
      real(8) :: qq(3), aa(3), dqq(3), daa(3)
      real(8) :: cg(3), rt
      integer :: p, k
      integer :: pa(3), pb(3)
      data pa /1, 2, 1/  ! pair atoms (A<B): (1,2),(2,3),(1,3)
      data pb /2, 3, 3/

      h_ij = 0.0d0
      dh_ij = 0.0d0
      do p = 1, 3
         call pair_dist(q(3*pa(p)-2:3*pa(p)), q(3*pb(p)-2:3*pb(p)), r(p))
         call leps_branch(r(p), qq(p), aa(p), dqq(p), daa(p))
      end do

      cg = aa - qq                      ! the singlet-triplet gaps A_p - Q_p
      rt = sqrt(dot_product(cg, cg))
      h_ij(1, 1) = 0.5d0*sum(qq)        ! H = [ s/2  -rt/2 ]
      h_ij(2, 2) = h_ij(1, 1)           !     [ -rt/2  s/2 ]
      h_ij(1, 2) = -0.5d0*rt
      h_ij(2, 1) = h_ij(1, 2)

      do p = 1, 3
         ! r_p = |x_B - x_A|: dr/dx_A = -u, dr/dx_B = +u
         call pair_unit(q(3*pa(p)-2:3*pa(p)), q(3*pb(p)-2:3*pb(p)), u(:,p))
         do k = 1, 3
            ! dH_11/dr_p = dH_22/dr_p = (1/2) dQ_p/dr_p
            dh_ij(3*(pa(p)-1)+k, 1, 1) = dh_ij(3*(pa(p)-1)+k, 1, 1) - 0.5d0*dqq(p)*u(k, p)
            dh_ij(3*(pb(p)-1)+k, 1, 1) = dh_ij(3*(pb(p)-1)+k, 1, 1) + 0.5d0*dqq(p)*u(k, p)
            ! dH_12/dr_p = -(1/2) (A_p-Q_p)(dA_p/dr_p - dQ_p/dr_p) / rt
            dh_ij(3*(pa(p)-1)+k, 1, 2) = dh_ij(3*(pa(p)-1)+k, 1, 2) &
                                         + 0.5d0*cg(p)*(daa(p) - dqq(p))*u(k, p)/rt
            dh_ij(3*(pb(p)-1)+k, 1, 2) = dh_ij(3*(pb(p)-1)+k, 1, 2) &
                                         - 0.5d0*cg(p)*(daa(p) - dqq(p))*u(k, p)/rt
         end do
      end do
      dh_ij(:, 2, 1) = dh_ij(:, 1, 2)   ! full-matrix delivery: the second
      dh_ij(:, 2, 2) = dh_ij(:, 1, 1)   ! triangle of every face, complete
   end subroutine leps_diab

   !------------------------------------------------------------------
   ! pair_dist(a, b, r) - the distance of two atomic positions [A]
   subroutine pair_dist(a, b, r)
      real(8), intent(in)  :: a(3), b(3)
      real(8), intent(out) :: r
      real(8) :: d(3)
      d = b - a
      r = sqrt(dot_product(d, d))
   end subroutine pair_dist

   !------------------------------------------------------------------
   ! pair_unit(a, b, u) - the unit vector along b - a [-]
   subroutine pair_unit(a, b, u)
      real(8), intent(in)  :: a(3), b(3)
      real(8), intent(out) :: u(3)
      real(8) :: d(3)
      d = b - a
      u = d/sqrt(dot_product(d, d))
   end subroutine pair_unit

   !------------------------------------------------------------------
   ! leps_branch(r, qe, ae, dqe, dae) - the Q (triplet) and A (singlet)
   !                 branch energies [eV] and their r-derivatives [eV/A]
   subroutine leps_branch(r, qe, ae, dqe, dae)
      real(8), intent(in)  :: r
      real(8), intent(out) :: qe, ae, dqe, dae
      real(8) :: em, dem
      em  = exp(-l_a*(r - l_re))
      dem = -l_a*em
      qe = 0.5d0*l_de*((3.0d0 + l_delta)*em*em - (2.0d0 + 2.0d0*l_delta)*em)
      ae = 0.5d0*l_de*((3.0d0 - l_delta)*em*em + (2.0d0 - 2.0d0*l_delta)*em) &
           + l_delta*l_de
      dqe = 0.5d0*l_de*(2.0d0*(3.0d0 + l_delta)*em*dem &
                        - (2.0d0 + 2.0d0*l_delta)*dem)
      dae = 0.5d0*l_de*(2.0d0*(3.0d0 - l_delta)*em*dem &
                        + (2.0d0 - 2.0d0*l_delta)*dem)
   end subroutine leps_branch
end module reference_pot
