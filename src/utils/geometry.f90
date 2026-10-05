!=====================================================================
! geometry.f90 - geometric transforms: orientation rotations, center-of-mass
!                separation, angular momentum/velocity, alignment steps, quadrature
! Design:
!   Stateless routines over (state, mass, list) triples; all angles in RADIANS.
!   The COM-frame workspace q_cm_frm/p_cm_frm is module-private and rebuilt by
!   cenmas on every call (derived quantities do not enter state) - cenmas is the
!   precondition of rot_apply/rot_euler/rot_jkm/angvel/amom. rot_euler is the
!   fused draw+apply convenience; rot_rand_rmat draws the matrix (the only RNG
!   site) and rot_apply writes the state (no RNG - the realize half).
!=====================================================================
module geometry
   use consts, only: pi, two_pi, e_conv
   use state,  only: state_t
   use rng,    only: rng_u
   implicit none
   private
   public :: rot_euler, rot_rand_rmat, rot_apply, rot_jkm, rot_axis, cenmas, angvel, amom, align_step, gl_nodes

   ! Center-of-mass frame workspace (derived quantities do not enter state; rebuilt
   ! by cenmas on every call):
   real(8), allocatable :: q_cm_frm(:)  ! center-of-mass-frame coordinates [Å]
   real(8), allocatable :: p_cm_frm(:)  ! center-of-mass-frame momenta [amu·Å/(10 fs)]
contains
   !------------------------------------------------------------------
   ! cenmas(st, mass, list, q_cm, v_cm) - center-of-mass position/velocity computation + separation of
   ! center-of-mass-frame coordinates and momenta
   subroutine cenmas(st, mass, list, q_cm, v_cm)
      type(state_t), intent(in) :: st    ! physical state (read-only; in-list atoms enter the sums)
      real(8), intent(in) :: mass(:)     ! per-atom masses [amu]
      integer, intent(in) :: list(:)     ! fragment atom index list [global atom index]
      real(8), intent(out) :: q_cm(3)    ! center-of-mass coordinates [Å]
      real(8), intent(out) :: v_cm(3)    ! center-of-mass velocity [Å/(10 fs)]
      integer :: i, j, j3, j2, j1
      real(8) :: wt
      ! 1. center-of-mass position and velocity sums over the atom list
      wt = sum(mass(list))
      do i = 1, 3
         v_cm(i) = 0.0d0
         q_cm(i) = 0.0d0
      end do
      do i = 1, size(list)
         j = list(i)
         j3 = 3*j; j2 = j3 - 1; j1 = j2 - 1
         v_cm(1) = v_cm(1) + st%p(j1)
         v_cm(2) = v_cm(2) + st%p(j2)
         v_cm(3) = v_cm(3) + st%p(j3)
         q_cm(1) = q_cm(1) + mass(j)*st%q(j1)
         q_cm(2) = q_cm(2) + mass(j)*st%q(j2)
         q_cm(3) = q_cm(3) + mass(j)*st%q(j3)
      end do
      do i = 1, 3
         v_cm(i) = v_cm(i)/wt
         q_cm(i) = q_cm(i)/wt
      end do
      ! 2. rebuild the COM-frame workspace (q/p minus the COM position/velocity)
      if (allocated(q_cm_frm)) deallocate (q_cm_frm)
      if (allocated(p_cm_frm)) deallocate (p_cm_frm)
      allocate (q_cm_frm(size(st%q)), p_cm_frm(size(st%q)))
      do i = 1, size(list)
         j = list(i)
         j3 = 3*j; j2 = j3 - 1; j1 = j2 - 1
         p_cm_frm(j1) = st%p(j1) - mass(j)*v_cm(1)
         p_cm_frm(j2) = st%p(j2) - mass(j)*v_cm(2)
         p_cm_frm(j3) = st%p(j3) - mass(j)*v_cm(3)
         q_cm_frm(j1) = st%q(j1) - q_cm(1)
         q_cm_frm(j2) = st%q(j2) - q_cm(2)
         q_cm_frm(j3) = st%q(j3) - q_cm(3)
      end do
   end subroutine cenmas

   !------------------------------------------------------------------
   ! rot_rand_rmat(rmat) - draw one isotropic uniform rotation matrix
   ! (three Euler-angle uniforms; the draw half of rot_euler)
   subroutine rot_rand_rmat(rmat)
      real(8), intent(out) :: rmat(3,3)  ! the drawn rotation matrix [-]
      real(8) :: rand, rphi, csthta, rchi, rthta, snthta, snphi, csphi, snchi, cschi
      real(8) :: rxx, rxy, rxz, ryx, ryy, ryz, rzx, rzy, rzz
      rand = rng_u()
      rphi = two_pi*rand
      rand = rng_u()
      csthta = 2.0d0*rand - 1.0d0
      rand = rng_u()
      rchi = two_pi*rand
      rthta = acos(csthta)
      snthta = sin(rthta)
      snphi = sin(rphi)
      csphi = cos(rphi)
      snchi = sin(rchi)
      cschi = cos(rchi)
      ! compute into scalars first (the historical inline form - keeps compiler
      ! contraction patterns stable across the draw/apply split), then store
      rxx = csthta*csphi*cschi - snphi*snchi
      rxy = -csthta*csphi*snchi - snphi*cschi
      rxz = snthta*csphi
      ryx = csthta*snphi*cschi + csphi*snchi
      ryy = -csthta*snphi*snchi + csphi*cschi
      ryz = snthta*snphi
      rzx = -snthta*cschi
      rzy = snthta*snchi
      rzz = csthta
      rmat(1,1) = rxx
      rmat(1,2) = rxy
      rmat(1,3) = rxz
      rmat(2,1) = ryx
      rmat(2,2) = ryy
      rmat(2,3) = ryz
      rmat(3,1) = rzx
      rmat(3,2) = rzy
      rmat(3,3) = rzz
   end subroutine rot_rand_rmat

   !------------------------------------------------------------------
   ! rot_apply(st, mass, list, rmat) - rotate a fragment about its center of
   ! mass by rmat (coordinates and momenta transformed together, written back
   ! COM-re-centered with zero net fragment momentum; the apply half of
   ! rot_euler - no RNG draws)
   subroutine rot_apply(st, mass, list, rmat)
      type(state_t), intent(inout) :: st  ! physical state (the in-list q/p components rotated in place)
      real(8), intent(in) :: mass(:)      ! per-atom masses [amu]
      integer, intent(in) :: list(:)      ! fragment atom index list
      real(8), intent(in) :: rmat(3,3)    ! the rotation to apply [-]
      real(8) :: rxx, rxy, rxz, ryx, ryy, ryz, rzx, rzy, rzz
      real(8) :: qcm(3), vcm(3)
      integer :: i, j
      ! copy to scalars first (the historical inline form - keeps compiler
      ! contraction patterns stable across the draw/apply split)
      rxx = rmat(1,1); rxy = rmat(1,2); rxz = rmat(1,3)
      ryx = rmat(2,1); ryy = rmat(2,2); ryz = rmat(2,3)
      rzx = rmat(3,1); rzy = rmat(3,2); rzz = rmat(3,3)
      call cenmas(st, mass, list, qcm, vcm)   ! the COM-frame workspace precondition
      do i = 1, size(list)
         j = 3*list(i)
         st%q(j-2) = q_cm_frm(j-2)*rxx + q_cm_frm(j-1)*rxy + q_cm_frm(j)*rxz
         st%q(j-1) = q_cm_frm(j-2)*ryx + q_cm_frm(j-1)*ryy + q_cm_frm(j)*ryz
         st%q(j)   = q_cm_frm(j-2)*rzx + q_cm_frm(j-1)*rzy + q_cm_frm(j)*rzz
         st%p(j-2) = p_cm_frm(j-2)*rxx + p_cm_frm(j-1)*rxy + p_cm_frm(j)*rxz
         st%p(j-1) = p_cm_frm(j-2)*ryx + p_cm_frm(j-1)*ryy + p_cm_frm(j)*ryz
         st%p(j)   = p_cm_frm(j-2)*rzx + p_cm_frm(j-1)*rzy + p_cm_frm(j)*rzz
         q_cm_frm(j-2:j) = st%q(j-2:j)
         p_cm_frm(j-2:j) = st%p(j-2:j)
      end do
   end subroutine rot_apply

   !------------------------------------------------------------------
   ! rot_euler(st, mass, list) - rotate a fragment about its center of mass by random Euler
   ! angles (coordinates and momenta transformed together) - the fused
   ! draw+apply convenience atop rot_rand_rmat/rot_apply
   subroutine rot_euler(st, mass, list)
      type(state_t), intent(inout) :: st  ! physical state (the in-list q/p components rotated in place)
      real(8), intent(in) :: mass(:)      ! per-atom masses [amu]
      integer, intent(in) :: list(:)      ! fragment atom index list
      real(8) :: rmat(3,3)
      call rot_rand_rmat(rmat)
      call rot_apply(st, mass, list, rmat)
   end subroutine rot_euler

   !------------------------------------------------------------------
   ! rot_jkm(st, mass, list, j_vec, j_rot, m_proj) - orient a fragment by (J,M) rotations
   ! (rotate J onto the z axis + projection angle + random azimuth about J + random rotation
   ! about the normal)
   subroutine rot_jkm(st, mass, list, j_vec, j_rot, m_proj)
      type(state_t), intent(inout) :: st  ! physical state (the in-list q/p components rotated in place)
      real(8), intent(in) :: mass(:)      ! per-atom masses [amu]
      integer, intent(in) :: list(:)      ! fragment atom index list
      real(8), intent(in) :: j_vec(3)     ! fragment angular momentum vector [amu·Å²/(10 fs)]
      integer, intent(in) :: j_rot        ! rotational quantum number J
      integer, intent(in) :: m_proj       ! projection quantum number M
      real(8) :: xt2, yt2, zt2, xt3, yt3, zt3, xt4, zt4, rxx, rxy, ryx, ryy, rxz, rzx, rzz
      real(8) :: csthta, snthta, phi, qcm(3), vcm(3), p1(3), p2(3), x1(3), x2(3)
      integer :: i, j, k
      call cenmas(st, mass, list, qcm, vcm)   ! the COM-frame workspace precondition
      ! 1. rotate J about z into the xz plane
      xt2 = j_vec(1)
      yt2 = j_vec(2)
      zt2 = j_vec(3)
      rxx = xt2/sqrt(xt2**2 + yt2**2)
      rxy = yt2/sqrt(xt2**2 + yt2**2)
      ryx = -rxy
      ryy = rxx
      do i = 1, size(list)
         j = 3*list(i)
         st%q(j-2) = q_cm_frm(j-2)*rxx + q_cm_frm(j-1)*rxy
         st%q(j-1) = q_cm_frm(j-2)*ryx + q_cm_frm(j-1)*ryy
         st%p(j-2) = p_cm_frm(j-2)*rxx + p_cm_frm(j-1)*rxy
         st%p(j-1) = p_cm_frm(j-2)*ryx + p_cm_frm(j-1)*ryy
         q_cm_frm(j-2:j-1) = st%q(j-2:j-1)
         p_cm_frm(j-2:j-1) = st%p(j-2:j-1)
      end do
      xt3 = xt2*rxx + yt2*rxy
      yt3 = xt2*ryx + yt2*ryy
      zt3 = zt2
      ! 2. rotate J within the xz plane onto the z axis
      rxx = zt2/sqrt(xt2**2 + yt2**2 + zt2**2)
      rxz = -sqrt(xt2**2 + yt2**2)/sqrt(xt2**2 + yt2**2 + zt2**2)
      rzx = -rxz
      rzz = rxx
      do i = 1, size(list)
         j = 3*list(i)
         st%q(j-2) = q_cm_frm(j-2)*rxx + q_cm_frm(j)*rxz
         st%q(j)   = q_cm_frm(j-2)*rzx + q_cm_frm(j)*rzz
         st%p(j-2) = p_cm_frm(j-2)*rxx + p_cm_frm(j)*rxz
         st%p(j)   = p_cm_frm(j-2)*rzx + p_cm_frm(j)*rzz
         q_cm_frm(j-2) = st%q(j-2)
         q_cm_frm(j)   = st%q(j)
         p_cm_frm(j-2) = st%p(j-2)
         p_cm_frm(j)   = st%p(j)
      end do
      xt4 = xt3*rxx + zt3*rxz
      zt4 = xt3*rzx + zt3*rzz
      ! 3. J now on z; rotate about y so J makes acos(m/sqrt(J(J+1))) with the normal
      csthta = dble(m_proj)/sqrt(dble(j_rot*j_rot + j_rot))
      snthta = sqrt(1.0d0 - csthta**2)
      do i = 1, size(list)
         j = 3*list(i)
         st%q(j-2) = q_cm_frm(j-2)*csthta + q_cm_frm(j)*snthta
         st%q(j)   = -q_cm_frm(j-2)*snthta + q_cm_frm(j)*csthta
         st%p(j-2) = p_cm_frm(j-2)*csthta + p_cm_frm(j)*snthta
         st%p(j)   = -p_cm_frm(j-2)*snthta + p_cm_frm(j)*csthta
         q_cm_frm(j-2) = st%q(j-2)
         q_cm_frm(j)   = st%q(j)
         p_cm_frm(j-2) = st%p(j-2)
         p_cm_frm(j)   = st%p(j)
      end do
      ! 4. random azimuth about the current J axis (direction (snthta,0,csthta)); radians
      p1 = (/ 0.0d0, 0.0d0, 0.0d0 /)
      p2 = (/ snthta, 0.0d0, csthta /)
      phi = two_pi*rng_u()
      do i = 1, size(list)
         j = 3*list(i)
         do k = 1, 3
            x1(k) = q_cm_frm(j-3+k)
         end do
         call rot_axis(p1, p2, x1, phi, x2)
         do k = 1, 3
            st%q(j-3+k) = x2(k)
         end do
         do k = 1, 3
            x1(k) = p_cm_frm(j-3+k)
         end do
         call rot_axis(p1, p2, x1, phi, x2)
         do k = 1, 3
            st%p(j-3+k) = x2(k)
         end do
         q_cm_frm(j-2:j) = st%q(j-2:j)
         p_cm_frm(j-2:j) = st%p(j-2:j)
      end do
      ! 5. random rotation about the surface normal (cylindrical symmetry)
      phi = two_pi*rng_u()
      do i = 1, size(list)
         j = 3*list(i)
         do k = 1, 3
            x1(k) = q_cm_frm(j-3+k)
         end do
         call rot_axis(p1, (/ 0.0d0, 0.0d0, 1.0d0 /), x1, phi, x2)
         do k = 1, 3
            st%q(j-3+k) = x2(k)
         end do
         do k = 1, 3
            x1(k) = p_cm_frm(j-3+k)
         end do
         call rot_axis(p1, (/ 0.0d0, 0.0d0, 1.0d0 /), x1, phi, x2)
         do k = 1, 3
            st%p(j-3+k) = x2(k)
         end do
         q_cm_frm(j-2:j) = st%q(j-2:j)
         p_cm_frm(j-2:j) = st%p(j-2:j)
      end do
   end subroutine rot_jkm

   !------------------------------------------------------------------
   ! rot_axis(pt_a, pt_b, pt_in, ang, pt_out) - rotate a point about an arbitrary axis (pt_a->pt_b),
   ! Murray closed form
   subroutine rot_axis(pt_a, pt_b, pt_in, ang, pt_out)
      real(8), intent(in)  :: pt_a(3)    ! axis start coordinates [Å]
      real(8), intent(in)  :: pt_b(3)    ! axis end coordinates [Å] (axis direction = pt_b - pt_a)
      real(8), intent(in)  :: pt_in(3)   ! coordinates before rotation [Å]
      real(8), intent(in)  :: ang        ! rotation angle [rad] (counterclockwise)
      real(8), intent(out) :: pt_out(3)  ! coordinates after rotation [Å]
      real(8) :: a, b, c, u, v, w, dist, x, y, z
      a = pt_a(1)
      b = pt_a(2)
      c = pt_a(3)
      dist = sqrt((pt_b(1)-pt_a(1))**2 + (pt_b(2)-pt_a(2))**2 + (pt_b(3)-pt_a(3))**2)
      u = (pt_b(1)-pt_a(1))/dist
      v = (pt_b(2)-pt_a(2))/dist
      w = (pt_b(3)-pt_a(3))/dist
      x = pt_in(1)
      y = pt_in(2)
      z = pt_in(3)
      pt_out(1) = (a*(v*v + w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1.0d0 - cos(ang)) &
                  + x*cos(ang) + (-c*v + b*w - w*y + v*z)*sin(ang)
      pt_out(2) = (b*(u*u + w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1.0d0 - cos(ang)) &
                  + y*cos(ang) + (c*u - a*w + w*x - u*z)*sin(ang)
      pt_out(3) = (c*(u*u + v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1.0d0 - cos(ang)) &
                  + z*cos(ang) + (-b*u + a*v - v*x + u*y)*sin(ang)
   end subroutine rot_axis

   !------------------------------------------------------------------
   ! angvel(st, mass, list, omega) - subtract the rigid-body rotation contribution m*(omega x r) from
   ! the momenta
   subroutine angvel(st, mass, list, omega)
      type(state_t), intent(inout) :: st  ! physical state (the in-list st%p components subtracted in place)
      real(8), intent(in) :: mass(:)      ! per-atom masses [amu]
      integer, intent(in) :: list(:)      ! fragment atom index list
      real(8), intent(in) :: omega(3)     ! angular velocity [rad/(10 fs)] (from amom)
      integer :: i, j, j3, j2, j1
      ! precondition: cenmas has built the COM-frame workspace
      do i = 1, size(list)
         j = list(i)
         j3 = 3*j; j2 = j3 - 1; j1 = j2 - 1
         st%p(j1) = st%p(j1) - (q_cm_frm(j3)*omega(2) - q_cm_frm(j2)*omega(3))*mass(j)
         st%p(j2) = st%p(j2) - (q_cm_frm(j1)*omega(3) - q_cm_frm(j3)*omega(1))*mass(j)
         st%p(j3) = st%p(j3) - (q_cm_frm(j2)*omega(1) - q_cm_frm(j1)*omega(2))*mass(j)
         p_cm_frm(j1) = st%p(j1)
         p_cm_frm(j2) = st%p(j2)
         p_cm_frm(j3) = st%p(j3)
      end do
   end subroutine angvel

   !------------------------------------------------------------------
   ! amom(st, mass, list, j_vec, e_rot, omega) - angular momentum / rotational energy (including inertia-
   ! tensor inversion)
   subroutine amom(st, mass, list, j_vec, e_rot, omega)
      type(state_t), intent(in) :: st          ! physical state (read-only; r, p via the COM workspace)
      real(8), intent(in) :: mass(:)           ! per-atom masses [amu]
      integer, intent(in) :: list(:)           ! fragment atom index list
      real(8), intent(out) :: j_vec(3)         ! angular momentum vector [amu·Å²/(10 fs)]
      real(8), intent(out) :: e_rot            ! rotational energy [kcal/mol]
      real(8), intent(out), optional :: omega(3) ! angular velocity [rad/(10 fs)]
      real(8) :: aixx, aiyy, aizz, aixy, aixz, aiyz, det, sr
      real(8) :: uxx, uxy, uxz, uyy, uyz, uzz, wx, wy, wz
      integer :: i, j, j3, j2, j1, k
      j_vec = 0.0d0
      e_rot = 0.0d0
      if (present(omega)) omega = 0.0d0
      if (size(list) == 1) return
      ! 1. angular-momentum vector from the COM-frame workspace
      do i = 1, size(list)
         j = list(i)
         j3 = 3*j; j2 = j3 - 1; j1 = j2 - 1
         j_vec(1) = j_vec(1) + (q_cm_frm(j2)*p_cm_frm(j3) - q_cm_frm(j3)*p_cm_frm(j2))
         j_vec(2) = j_vec(2) + (q_cm_frm(j3)*p_cm_frm(j1) - q_cm_frm(j1)*p_cm_frm(j3))
         j_vec(3) = j_vec(3) + (q_cm_frm(j1)*p_cm_frm(j2) - q_cm_frm(j2)*p_cm_frm(j1))
      end do
      if (size(list) == 2) then
         ! the linear special case: E = |J|^2/(2I), I = sum m r^2
         aixx = 0.0d0
         do i = 1, size(list)
            j = 3*list(i) + 1
            sr = 0.0d0
            do k = 1, 3
               sr = sr + q_cm_frm(j-k)**2
            end do
            aixx = aixx + sr*mass(list(i))
         end do
         e_rot = (j_vec(1)**2 + j_vec(2)**2 + j_vec(3)**2)/aixx/2.0d0/e_conv
         return
      end if
      ! 2. the inertia tensor
      aixx = 0.0d0; aiyy = 0.0d0; aizz = 0.0d0
      aixy = 0.0d0; aixz = 0.0d0; aiyz = 0.0d0
      do i = 1, size(list)
         j = list(i)
         j3 = 3*j; j2 = j3 - 1; j1 = j2 - 1
         aixx = aixx + mass(j)*(q_cm_frm(j2)**2 + q_cm_frm(j3)**2)
         aiyy = aiyy + mass(j)*(q_cm_frm(j1)**2 + q_cm_frm(j3)**2)
         aizz = aizz + mass(j)*(q_cm_frm(j1)**2 + q_cm_frm(j2)**2)
         aixy = aixy + mass(j)*q_cm_frm(j1)*q_cm_frm(j2)
         aixz = aixz + mass(j)*q_cm_frm(j1)*q_cm_frm(j3)
         aiyz = aiyz + mass(j)*q_cm_frm(j2)*q_cm_frm(j3)
      end do
      det = aixx*(aiyy*aizz - aiyz*aiyz) - aixy*(aixy*aizz + aiyz*aixz) - &
            aixz*(aixy*aiyz + aiyy*aixz)
      if (abs(det) >= 0.01d0) then
         ! inverse by the adjugate + omega = I^-1 J
         uxx = (aiyy*aizz - aiyz*aiyz)/det
         uxy = (aixy*aizz + aixz*aiyz)/det
         uxz = (aixy*aiyz + aixz*aiyy)/det
         uyy = (aixx*aizz - aixz*aixz)/det
         uyz = (aixx*aiyz + aixz*aixy)/det
         uzz = (aixx*aiyy - aixy*aixy)/det
         wx = uxx*j_vec(1) + uxy*j_vec(2) + uxz*j_vec(3)
         wy = uxy*j_vec(1) + uyy*j_vec(2) + uyz*j_vec(3)
         wz = uxz*j_vec(1) + uyz*j_vec(2) + uzz*j_vec(3)
      else
         ! near-singular (linear geometry): the linear special-case formula
         aixx = 0.0d0
         do i = 1, size(list)
            j = 3*list(i) + 1
            sr = 0.0d0
            do k = 1, 3
               sr = sr + q_cm_frm(j-k)**2
            end do
            aixx = aixx + sr*mass(list(i))
         end do
         e_rot = (j_vec(1)**2 + j_vec(2)**2 + j_vec(3)**2)/aixx/2.0d0/e_conv
         return
      end if
      if (present(omega)) omega = (/ wx, wy, wz /)
      e_rot = (wx*j_vec(1) + wy*j_vec(2) + wz*j_vec(3))/2.0d0/e_conv
   end subroutine amom

   !------------------------------------------------------------------
   ! align_step(st, i_axis, ang, list) - one rigid-body alignment step (about z or about y, through the
   ! origin; coordinates and momenta transformed together)
   subroutine align_step(st, i_axis, ang, list)
      type(state_t), intent(inout) :: st  ! physical state (the in-list q/p components rotated in place)
      integer, intent(in) :: i_axis       ! rotation axis (3 = about z; 2 = about y, ang = target cos)
      real(8), intent(in) :: ang          ! about z: rotation angle [rad]; about y: target direction cosine
      integer, intent(in) :: list(:)      ! alignment atom index list
      real(8) :: cp, sp, x, y, z, px, py, pz
      integer :: i, j
      if (i_axis == 3) then
         cp = cos(ang)
         sp = sin(ang)
         do i = 1, size(list)
            j = 3*list(i)
            x = st%q(j-2); y = st%q(j-1)
            st%q(j-2) = x*cp + y*sp
            st%q(j-1) = -x*sp + y*cp
            px = st%p(j-2); py = st%p(j-1)
            st%p(j-2) = px*cp + py*sp
            st%p(j-1) = -px*sp + py*cp
         end do
      else if (i_axis == 2) then
         cp = ang
         sp = sqrt(max(0.0d0, 1.0d0 - cp*cp))
         do i = 1, size(list)
            j = 3*list(i)
            x = st%q(j-2); z = st%q(j)
            st%q(j-2) = x*cp - z*sp
            st%q(j) = x*sp + z*cp
            px = st%p(j-2); pz = st%p(j)
            st%p(j-2) = px*cp - pz*sp
            st%p(j) = px*sp + pz*cp
         end do
      else
         write (0, '(a,i0)') 'align_step: unknown axis code ', i_axis
         stop 1
      end if
   end subroutine align_step

   !------------------------------------------------------------------
   ! gl_nodes(x_lo, x_hi, node, wts) - Gauss-Legendre quadrature nodes and weights (symmetric halving +
   ! Newton refinement)
   subroutine gl_nodes(x_lo, x_hi, node, wts)
      real(8), intent(in)  :: x_lo       ! integration lower bound
      real(8), intent(in)  :: x_hi       ! integration upper bound
      real(8), intent(out) :: node(:)    ! quadrature nodes (n = size(node))
      real(8), intent(out) :: wts(:)     ! quadrature weights (aligned with node)
      integer :: n, m, i
      real(8) :: xx, dpl, dum, xm, xl
      n = size(node)
      m = (n + 1)/2
      xm = 0.5d0*(x_hi + x_lo)
      xl = 0.5d0*(x_hi - x_lo)
      do i = 1, m
         dum = 1.0d0
         xx = cos(pi*(i - 0.25d0)/(n + 0.5d0))
         do while (abs(dum) > 1.0d-13)
            dum = pl(xx, n, dpl)
            xx = xx - dum/dpl
         end do
         node(i) = xm - xl*xx
         node(n+1-i) = xm + xl*xx
         wts(i) = xl/((1.0d0 - xx*xx)*dpl*dpl)/0.5d0
         wts(n+1-i) = wts(i)
      end do
   contains
      ! Legendre polynomial of order n + its derivative
      real(8) function pl(x, n, dpl)
         real(8), intent(in) :: x
         integer, intent(in) :: n
         real(8), intent(out) :: dpl
         real(8) :: p1, p2
         integer :: k
         p2 = x
         pl = 1.5d0*x*x - 0.5d0
         do k = 2, n - 1
            p1 = p2
            p2 = pl
            pl = (dble(2*k + 1)*p2*x - dble(k)*p1)/dble(k + 1)
         end do
         dpl = dble(n)*(x*pl - p2)/(x*x - 1.0d0)
      end function pl
   end subroutine gl_nodes
end module geometry
