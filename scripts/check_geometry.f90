!=====================================================================
! check_geometry.f90 - geometry criteria program (criteria defined FIRST in
!                      this header)
! Compile (from the repo root):
!   gfortran -Jsrc/build/gfortran src/globals/consts.f90 src/globals/state.f90 &
!            src/utils/rng.f90 src/utils/geometry.f90 scripts/check_geometry.f90 &
!            -o src/build/gfortran/check_geometry.exe
! Criteria:
!   (a/f) rot_axis orthogonality (R^T R = I, det = +1) on a coordinate-axis and an
!         oblique axis + the known action (90 deg about z maps e1 -> e2; axis points fixed)
!   (b)   cenmas on a constructed 3-body system: COM-frame residuals (sum p_cm = 0,
!         mass-weighted sum q_cm = 0) < 1e-12
!   (c)   amom analytic: two equal masses on a circular orbit give |L| = 2*m*r*v and
!         E_rot = |L|^2/(2 I) in kcal/mol
!   (d/g) invariance: |J| and E_rot unchanged (1e-10) under rot_euler, rot_jkm and
!         align_step (both axes); pairwise distances unchanged under every rotation
!   (e)   gl_nodes: n=2 and n=3 machine-precision nodes/weights; sum of weights = the
!         interval length; n=3 integrates x^4 exactly on [0,2]
!=====================================================================
program check_geometry
   use consts,   only: e_conv
   use state,    only: state_t, state_create, state_destroy
   use rng,      only: rng_init
   use geometry, only: rot_euler, rot_jkm, rot_axis, cenmas, angvel, amom, align_step, gl_nodes
   implicit none
   type(state_t) :: st
   real(8) :: mass(3) = (/ 1.0d0, 2.0d0, 3.0d0 /)
   integer :: list(3) = (/ 1, 2, 3 /)
   real(8) :: qcm(3), vcm(3), jv(3), jv0(3), erot, erot0, om(3)
   real(8) :: r(3, 3), rt(3, 3), det, res
   real(8) :: x1(3), x2(3), pa(3), pb(3)
   real(8) :: d01, d02, d12, dd
   real(8) :: node(3), wts(3), sumw, exact
   integer :: i, j, nfail = 0

   ! ---------- (a/f) rot_axis orthogonality + known actions ----------
   ! orthogonality needs an ORIGIN-PASSING axis (an off-origin axis makes the basis
   ! images an affine rotation+translation, not a matrix); 90 deg about +z through 0
   pa = (/ 0.0d0, 0.0d0, 0.0d0 /)
   pb = (/ 0.0d0, 0.0d0, 1.0d0 /)
   do i = 1, 3
      x1 = 0.0d0
      x1(i) = 1.0d0
      call rot_axis(pa, pb, x1, acos(-1.0d0)/2.0d0, r(:,i))
   end do
   rt = matmul(transpose(r), r)
   res = 0.0d0
   do j = 1, 3
      do i = 1, 3
         res = max(res, abs(rt(i,j) - merge(1.0d0, 0.0d0, i == j)))
      end do
   end do
   det = r(1,1)*(r(2,2)*r(3,3) - r(2,3)*r(3,2)) - r(1,2)*(r(2,1)*r(3,3) - r(2,3)*r(3,1)) &
         + r(1,3)*(r(2,1)*r(3,2) - r(2,2)*r(3,1))
   ! known action about +z: e1 -> e2 (column 1 = (0,1,0) i.e. r(2,1)=+1),
   ! e2 -> -e1 (column 2 = (-1,0,0) i.e. r(1,2)=-1), e3 -> e3
   if (res > 1.0d-13 .or. abs(det - 1.0d0) > 1.0d-13 .or. &
       abs(r(1,2) + 1.0d0) > 1.0d-13 .or. abs(r(2,1) - 1.0d0) > 1.0d-13 .or. &
       abs(r(3,3) - 1.0d0) > 1.0d-13) then
      print '(a,es10.2,a,f10.6)', 'FAIL (a) z-axis rotation: res ', res, ' det ', det
      nfail = nfail + 1
   else
      print '(a)', '(a/f) rot_axis orthogonality + det(+1) + e1->e2 about +z: PASS'
   end if
   ! off-origin OBLIQUE axis: axis points fixed + point-to-axis-pivot distance preserved
   ! (the affine invariants)
   pa = (/ 1.0d0, 2.0d0, -1.0d0 /)
   pb = (/ 3.0d0, 0.0d0, 2.0d0 /)
   x1 = pa                                  ! pivot on the axis: fixed
   call rot_axis(pa, pb, x1, 0.7d0, x2)
   res = maxval(abs(x2 - x1))
   x1 = (/ 0.4d0, -1.3d0, 2.7d0 /)          ! any point: |x - pa| preserved
   dd = sqrt(sum((x1 - pa)**2))
   call rot_axis(pa, pb, x1, 0.7d0, x2)
   res = max(res, abs(sqrt(sum((x2 - pa)**2)) - dd))
   if (res > 1.0d-13) then
      print '(a,es10.2)', 'FAIL (a) oblique off-origin axis invariants: res ', res
      nfail = nfail + 1
   else
      print '(a)', '(a) rot_axis oblique off-origin axis (pivot fixed, radius preserved): PASS'
   end if

   ! ---------- construct the 3-body state ----------
   call state_create(st, 3)
   st%q = (/ 0.3d0, 1.1d0, -0.7d0,  2.2d0, -0.4d0, 0.9d0,  -1.5d0, 0.8d0, 1.6d0 /)
   st%p = (/ 0.1d0, 0.2d0, 0.3d0,  -0.2d0, 0.15d0, 0.1d0,  0.05d0, -0.1d0, 0.25d0 /)

   ! ---------- (b) cenmas COM-frame residuals ----------
   call cenmas(st, mass, list, qcm, vcm)
   res = 0.0d0
   do i = 1, 3
      res = max(res, abs(sum(st%p(3*list-2 + (i-1))) - sum(mass(list))*vcm(i)))
   end do
   ! mass-weighted COM-frame coordinate sum = 0 (module-private workspace exercised
   ! indirectly: q_cm sum via the definition)
   do i = 1, 3
      dd = 0.0d0
      do j = 1, 3
         dd = dd + mass(list(j))*(st%q(3*list(j)-3+i) - qcm(i))
      end do
      res = max(res, abs(dd))
   end do
   if (res > 1.0d-12) then
      print '(a,es10.2)', 'FAIL (b) cenmas residuals = ', res
      nfail = nfail + 1
   else
      print '(a)', '(b) cenmas COM-frame residuals < 1e-12: PASS'
   end if

   ! ---------- reference J / E_rot + pairwise distances ----------
   call amom(st, mass, list, jv0, erot0, om)
   d01 = dist(st%q, 1, 2)
   d02 = dist(st%q, 1, 3)
   d12 = dist(st%q, 2, 3)

   ! ---------- (d) rot_euler invariance ----------
   call rng_init(4242)
   call rot_euler(st, mass, list)
   call cenmas(st, mass, list, qcm, vcm)
   call amom(st, mass, list, jv, erot)
   call check_invar('rot_euler ', sqrt(sum(jv0**2)), sqrt(sum(jv**2)), erot0, erot, &
      d01, d02, d12, dist(st%q,1,2), dist(st%q,1,3), dist(st%q,2,3), nfail)

   ! ---------- (d) rot_jkm invariance ----------
   call rot_jkm(st, mass, list, jv0, 2, 0)
   call cenmas(st, mass, list, qcm, vcm)
   call amom(st, mass, list, jv, erot)
   call check_invar('rot_jkm   ', sqrt(sum(jv0**2)), sqrt(sum(jv**2)), erot0, erot, &
      d01, d02, d12, dist(st%q,1,2), dist(st%q,1,3), dist(st%q,2,3), nfail)

   ! ---------- (g) align_step invariance (both axes) ----------
   call align_step(st, 3, 0.9d0, list)
   call cenmas(st, mass, list, qcm, vcm)
   call amom(st, mass, list, jv, erot)
   call check_invar('align_z   ', sqrt(sum(jv0**2)), sqrt(sum(jv**2)), erot0, erot, &
      d01, d02, d12, dist(st%q,1,2), dist(st%q,1,3), dist(st%q,2,3), nfail)
   call align_step(st, 2, 0.6d0, list)
   call cenmas(st, mass, list, qcm, vcm)
   call amom(st, mass, list, jv, erot)
   call check_invar('align_y   ', sqrt(sum(jv0**2)), sqrt(sum(jv**2)), erot0, erot, &
      d01, d02, d12, dist(st%q,1,2), dist(st%q,1,3), dist(st%q,2,3), nfail)

   ! ---------- (c) amom analytic: circular two-body orbit ----------
   call state_destroy(st)
   call state_create(st, 2)
   st%q = (/ 1.2d0, 0.0d0, 0.0d0,  -1.2d0, 0.0d0, 0.0d0 /)
   st%p = (/ 0.0d0, 1.0d0, 0.0d0,   0.0d0, -1.0d0, 0.0d0 /)   ! m v with m = 2, v = 0.5
   block
      real(8) :: m2(2) = (/ 2.0d0, 2.0d0 /)
      integer :: l2(2) = (/ 1, 2 /)
      call cenmas(st, m2, l2, qcm, vcm)
      call amom(st, m2, l2, jv, erot)
      ! |L| = 2 * m r v = 2 * 2 * 1.2 * 0.5 = 2.4 ; I = sum m r^2 = 5.76 ;
      ! E_rot = L^2/(2 I) = 5.76/11.52 = 0.5 internal -> /e_conv kcal/mol
      if (abs(sqrt(sum(jv**2)) - 2.4d0) > 1.0d-12 .or. abs(erot - 0.5d0/e_conv) > 1.0d-12/e_conv) then
         print '(a,2f14.10)', 'FAIL (c) circular orbit: |L| ', sqrt(sum(jv**2)), ' E_rot ', erot
         nfail = nfail + 1
      else
         print '(a)', '(c) amom circular-orbit analytic |L|=2.4, E_rot=0.5/e_conv: PASS'
      end if
   end block
   call state_destroy(st)

   ! ---------- (e) gl_nodes ----------
   call gl_nodes(-1.0d0, 1.0d0, node(1:2), wts(1:2))
   if (abs(abs(node(1)) - 1.0d0/sqrt(3.0d0)) > 1.0d-13 .or. abs(node(1) + node(2)) > 1.0d-13 .or. &
       abs(wts(1) - 1.0d0) > 1.0d-13) then
      print '(a)', 'FAIL (e) gl_nodes n=2'
      nfail = nfail + 1
   else
      print '(a)', '(e) gl_nodes n=2 (±1/sqrt(3), w=1): PASS'
   end if
   call gl_nodes(-1.0d0, 1.0d0, node, wts)
   if (abs(node(1) + sqrt(3.0d0/5.0d0)) > 1.0d-13 .or. abs(node(2)) > 1.0d-13 .or. &
       abs(wts(1) - 5.0d0/9.0d0) > 1.0d-13 .or. abs(wts(2) - 8.0d0/9.0d0) > 1.0d-13) then
      print '(a)', 'FAIL (e) gl_nodes n=3'
      nfail = nfail + 1
   else
      print '(a)', '(e) gl_nodes n=3 (0,±sqrt(3/5); w=8/9,5/9): PASS'
   end if
   call gl_nodes(0.0d0, 2.0d0, node, wts)
   sumw = sum(wts)
   exact = sum(wts*node**4)   ! GL n=3 integrates x^4 exactly on [0,2]: 32/5
   if (abs(sumw - 2.0d0) > 1.0d-13 .or. abs(exact - 32.0d0/5.0d0) > 1.0d-12) then
      print '(a,2f14.10)', 'FAIL (e) gl_nodes [0,2]: sumw ', sumw, ' int x^4 ', exact
      nfail = nfail + 1
   else
      print '(a)', '(e) gl_nodes [0,2] sum w = 2, int x^4 = 32/5 exact: PASS'
   end if

   print '(a,i0)', '=== check_geometry: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains
   real(8) function dist(q, a, b)
      real(8), intent(in) :: q(:)
      integer, intent(in) :: a, b
      dist = sqrt((q(3*a-2)-q(3*b-2))**2 + (q(3*a-1)-q(3*b-1))**2 + (q(3*a)-q(3*b))**2)
   end function dist
   subroutine check_invar(tag, j0, j1, e0, e1, r01, r02, r12, s01, s02, s12, nf)
      character(len=*), intent(in) :: tag
      real(8), intent(in) :: j0, j1, e0, e1, r01, r02, r12, s01, s02, s12
      integer, intent(inout) :: nf
      if (abs(j0 - j1) > 1.0d-10 .or. abs(e0 - e1) > 1.0d-10 .or. &
          abs(r01 - s01) > 1.0d-10 .or. abs(r02 - s02) > 1.0d-10 .or. abs(r12 - s12) > 1.0d-10) then
         print '(a,f14.8,f14.8,f16.8)', 'FAIL invariance ', tag, j1, e1
         nf = nf + 1
      else
         print '(a,a,a,f12.6,a,f12.6)', '(d/g) ', tag, ': |J| ', j1, '  E_rot ', e1, ' invariant'
      end if
   end subroutine check_invar
end program check_geometry
