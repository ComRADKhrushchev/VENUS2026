!=====================================================================
! geometry.f90 - L5 instrument: geometric transforms (orientation rotations / center-of-mass
!                separation / angular velocity and angular momentum / alignment steps / quadrature nodes)
! Reference: src_VENUS/ROTATE.f90:19-54 (random Euler orientation); src_VENUS/ROTATEJKM.f:21-108 ((J,M) orientation);
!       src_VENUS/ArbitraryAxisRotation.f:18-75 (arbitrary-axis rotation); src_VENUS/CENMAS.f90:16-50 (center-of-mass separation);
!       src_VENUS/ANGVEL.f:9-20 (angular-velocity subtraction); src_VENUS/ROTN.f90:18-104 (angular momentum / rotational energy);
!       src_TEST/interface_TEST.f90:76-108 (rigid-body alignment steps Z/Y); src_VENUS/GLPAR.f:23-69 (Gauss-Legendre)
!=====================================================================
module geometry
   use consts, only: pi, two_pi
   use state,  only: state_t
   implicit none
   private
   public :: rot_euler, rot_jkm, rot_axis, cenmas, angvel, amom, align_step, gl_nodes

   ! Center-of-mass frame workspace (derived quantities do not enter state - K1 decision D2; the old
   ! QQ/PP common block internalized as module-private;
   ! sizes grow with the fragment atom count, lifetimes finalized in phase C):
   real(8), allocatable :: q_cm_frm(:)  ! center-of-mass-frame coordinates [Å] (old QQ; built by cenmas,
   real(8), allocatable :: p_cm_frm(:)  ! center-of-mass-frame momenta [amu·Å/(10 fs)] (old PP)
contains
   !------------------------------------------------------------------
   ! rot_euler(st, mass, list) - rotate a fragment about its center of mass by random Euler
   ! angles (coordinates and momenta transformed together)
   ! Reference: src_VENUS/ROTATE.f90:19-54
   subroutine rot_euler(st, mass, list)
      type(state_t), intent(inout) :: st  ! physical state (the in-list q/p components are
      ! rotated in place)
      real(8), intent(in) :: mass(:)      ! per-atom masses [amu] (roster%mass, indexed by global atom)
      integer, intent(in) :: list(:)      ! fragment atom index list [global atom index] (old LL)
      ! Pseudocode:
      ! 1. Separate the center of mass of the fragment's atoms (reuse the cenmas steps to get
      ! the COM-frame workspace)
      ! 2. Draw u three times (rng_u): phi = two_pi*u, costht = 2*u-1, chi = two_pi*u;
      !    thta = acos(costht) (old :19-25)
      ! 3. Expand the nine elements of the 3x3 rotation matrix from (phi,thta,chi) (old :26-39)
      ! 4. Right-multiply the COM-frame coordinates and momenta of every atom by that matrix,
      ! write back to st%q/st%p and sync the workspace (old :40-54)
   end subroutine rot_euler

   !------------------------------------------------------------------
   ! rot_jkm(st, mass, list, j_vec, j_rot, m_proj) - orient a fragment by (J,M) rotations
   ! (rotate J onto the z axis + projection angle + random azimuth about J + random rotation
   ! about the normal)
   ! Reference: src_VENUS/ROTATEJKM.f:21-108; src_VENUS/ROTATEY.f:13-29;
   !       src_VENUS/ROTATEZ.f:13-27; src_VENUS/ArbitraryAxisRotation.f:63-72
   subroutine rot_jkm(st, mass, list, j_vec, j_rot, m_proj)
      type(state_t), intent(inout) :: st  ! physical state (the in-list q/p components are
      ! rotated in place)
      real(8), intent(in) :: mass(:)      ! per-atom masses [amu]
      integer, intent(in) :: list(:)      ! fragment atom index list [global atom index]
      real(8), intent(in) :: j_vec(3)     ! fragment angular momentum vector [amu·Å²/(10 fs)] (old AM(1:3), computed by amom)
      integer, intent(in) :: j_rot        ! rotational quantum number J [dimensionless] (old JROT)
      integer, intent(in) :: m_proj       ! projection quantum number M [dimensionless] (the old first dummy carried the M value)
      ! Pseudocode:
      ! 1. Separate the center of mass of the fragment's atoms (same precondition as rot_euler)
      ! 2. Rotate about z to bring j_vec into the xz plane (old :21-41)
      ! 3. Rotate within the xz plane to bring j_vec onto the z axis (old :43-62)
      ! 4. Rotate about y so that j makes the angle
      !    acos(m_proj/sqrt(j_rot*(j_rot+1))) with the surface normal (z) (old :64-70, the
      !    about-y step as in ROTATEY)
      ! 5. Random azimuth about the current j axis (through the origin, direction
      !    (sin(thta),0,cos(thta))):
      !    rotate by phi = two_pi*u, calling rot_axis on every atom's coordinates and momenta
      !    (old :72-103; the old version passed angles in degrees, the new one unifies on radians)
      ! 6. Uniform random rotation about z so the final distribution is cylindrically
      !    symmetric about the normal (old :105-108, the about-z step as in ROTATEZ)
   end subroutine rot_jkm

   !------------------------------------------------------------------
   ! rot_axis(pt_a, pt_b, pt_in, ang, pt_out) - rotate a point about an arbitrary axis (pt_a->pt_b),
   ! Murray closed form
   ! Reference: src_VENUS/ArbitraryAxisRotation.f:18-75 (closed-form expansion :63-72; axis
   ! normalization :32-35)
   subroutine rot_axis(pt_a, pt_b, pt_in, ang, pt_out)
      real(8), intent(in)  :: pt_a(3)    ! axis start coordinates [Å] (old p1)
      real(8), intent(in)  :: pt_b(3)    ! axis end coordinates [Å] (old p2; axis direction = pt_b - pt_a)
      real(8), intent(in)  :: pt_in(3)   ! coordinates before rotation [Å] (old x1)
      real(8), intent(in)  :: ang        ! rotation angle [rad] (counterclockwise; old th in degrees, the
      real(8), intent(out) :: pt_out(3)  ! coordinates after rotation [Å] (old x2)
      ! Pseudocode:
      ! 1. Axis unit vector (u,v,w) = (pt_b - pt_a)/|pt_b - pt_a| (old :32-35)
      ! 2. Expand pt_out by the three Murray closed-form components (cos(ang)/sin(ang) combinations,
      !    point-to-axis projection term + rotation term, old :63-72)
   end subroutine rot_axis

   !------------------------------------------------------------------
   ! cenmas(st, mass, list, q_cm, v_cm) - center-of-mass position/velocity computation + separation of
   ! center-of-mass-frame coordinates and momenta
   ! Reference: src_VENUS/CENMAS.f90:16-50
   subroutine cenmas(st, mass, list, q_cm, v_cm)
      type(state_t), intent(in) :: st    ! physical state (read-only; in-list atoms enter the sums)
      real(8), intent(in) :: mass(:)     ! per-atom masses [amu]
      integer, intent(in) :: list(:)     ! fragment atom index list [global atom index]
      real(8), intent(out) :: q_cm(3)    ! center-of-mass coordinates [Å]
      real(8), intent(out) :: v_cm(3)    ! center-of-mass velocity [Å/(10 fs)] (the old name VCM is
      ! Pseudocode:
      ! 1. v_cm = sum p(list)/sum m, q_cm = sum (m.q)(list)/sum m (old :21-37)
      ! 2. Write the separated quantities into the workspace: q_cm_frm = q - q_cm, p_cm_frm = p - m*v_cm
      !    (old :39-50; the workspace is reused by rot_euler/rot_jkm/angvel/amom - the old QQ/PP)
   end subroutine cenmas

   !------------------------------------------------------------------
   ! angvel(st, mass, list, omega) - subtract the rigid-body rotation contribution m*(omega x r) from
   ! the momenta
   ! Reference: src_VENUS/ANGVEL.f:9-20
   subroutine angvel(st, mass, list, omega)
      type(state_t), intent(inout) :: st  ! physical state (the in-list st%p components are
      ! subtracted in place)
      real(8), intent(in) :: mass(:)      ! per-atom masses [amu]
      integer, intent(in) :: list(:)      ! fragment atom index list [global atom index]
      real(8), intent(in) :: omega(3)     ! angular velocity [rad/(10 fs)] (old WX/WY/WZ, returned by amom)
      ! Pseudocode:
      ! 1. Precondition: cenmas has built the COM-frame workspace (old QQ)
      ! 2. Per atom st%p -= mass*(omega x r_cm) (three-component expansion, old :14-16)
      ! 3. Sync the COM-frame momentum workspace to the subtracted p (p_cm_frm = st%p, old
      ! :17-19)
   end subroutine angvel

   !------------------------------------------------------------------
   ! amom(st, mass, list, j_vec, e_rot, omega) - angular momentum / rotational energy (including inertia-
   ! tensor inversion)
   ! Reference: src_VENUS/ROTN.f90:18-104
   subroutine amom(st, mass, list, j_vec, e_rot, omega)
      type(state_t), intent(in) :: st          ! physical state (read-only; r, p taken via the
      ! COM-frame workspace)
      real(8), intent(in) :: mass(:)           ! per-atom masses [amu]
      integer, intent(in) :: list(:)           ! fragment atom index list [global atom index]
      real(8), intent(out) :: j_vec(3)         ! angular momentum vector [amu·Å²/(10 fs)] (old AM(1:3))
      real(8), intent(out) :: e_rot            ! rotational energy [kcal/mol] (old EROT = internal units /
      real(8), intent(out), optional :: omega(3) ! optional return of angular velocity [rad/(10 fs)] (old
      ! Pseudocode:
      ! 1. Precondition: cenmas has built the center-of-mass workspace
      ! 2. j_vec = sum r x p (three center-of-mass components, old :25-33); one-atom fragments get zero
      ! 3. Diatomic fragment: e_rot = |J|^2/(2*I), I = sum m*r^2 (linear special case, old :35-47)
      ! 4. Polyatomic: six inertia-tensor components (old :52-69) -> determinant -> inverse by the adjugate
      ! (old :70-81)
      ! 5. omega = I^-1 * j_vec (old :83-87); e_rot = (omega . j_vec)/2/e_conv (old :104)
      ! 6. Singular determinant (|det|<0.01) falls back to the linear special-case formula (old :88-102)
   end subroutine amom

   !------------------------------------------------------------------
   ! align_step(st, i_axis, ang, list) - one rigid-body alignment step (about z or about y, through the
   ! origin; coordinates and momenta transformed together)
   ! Reference: src_TEST/interface_TEST.f90:76-90 (Z step); :94-108 (Y step)
   subroutine align_step(st, i_axis, ang, list)
      type(state_t), intent(inout) :: st  ! physical state (the in-list q/p components are
      ! rotated in place)
      integer, intent(in) :: i_axis       ! rotation axis (3 = about z, old ALIGN_ZSTEP; 2 = about y, old
      real(8), intent(in) :: ang          ! about z: rotation angle phi [rad]; about y: target direction
      integer, intent(in) :: list(:)      ! alignment atom index list (old: fixed to the fragment's first
      ! Pseudocode:
      ! 1. About-z step: c=cos(ang), s=sin(ang), (x,y) -> (x*c+y*s, -x*s+y*c)
      !    (old :81-89, applied in parallel to coordinates and momenta)
      ! 2. About-y step: c=ang, s=sqrt(1-c^2), (x,z) -> (x*c-z*s, x*s+z*c)
      !    (old :98-107, applied in parallel to coordinates and momenta)
   end subroutine align_step

   !------------------------------------------------------------------
   ! gl_nodes(x_lo, x_hi, node, wts) - Gauss-Legendre quadrature nodes and weights (symmetric halving +
   ! Newton refinement)
   ! Reference: src_VENUS/GLPAR.f:23-50 (node/weight mapping); :55-69 (Legendre polynomial recursion PL)
   subroutine gl_nodes(x_lo, x_hi, node, wts)
      real(8), intent(in)  :: x_lo       ! integration lower bound [integrand-variable units]
      real(8), intent(in)  :: x_hi       ! integration upper bound [integrand-variable units]
      real(8), intent(out) :: node(:)    ! quadrature nodes [integrand-variable units] (n = size(node))
      real(8), intent(out) :: wts(:)     ! quadrature weights [integrand-variable units] (aligned with node)
      ! Pseudocode:
      ! 1. n = size(node); m = (n+1)/2 (symmetric halving, old :23)
      ! 2. Root initial guesses x0 = cos(pi*(i-0.25)/(n+0.5)) (old :28)
      ! 3. Newton-Raphson refinement x0 <- x0 - P_n(x0)/dP_n(x0), P_n by three-term recursion,
      !    dP_n = n*(x0*P_n - P_{n-1})/(x0^2-1), until |delta| < 1e-13 (old :32-35 and PL
      !    :60-68)
      ! 4. Map to [x_lo,x_hi]: midpoint +/- half-width*x0; weight = half-width/((1-x0^2)*dP_n^2);
      !    the symmetric half mirror-fills the other half (old :43-50)
   end subroutine gl_nodes
end module geometry
