!=====================================================================
! check_hessian.f90 - hessian instrument criteria program (criteria defined
!                     FIRST in this header; this program executes them)
! Compile (from src/; -Jbuild keeps the .mod files out of the source tree):
!   gfortran -Jbuild utils/hessian.f90 ../scripts/check_hessian.f90 -o build/check_hessian.exe
! All artifacts are written under build/ - the source tree stays clean.
! Loud-stop verification (same scheme as the sibling criteria programs):
!   every loud-stop criterion runs as a CHILD PROCESS - the check re-spawns
!   itself with a mode argument (execute_command_line); the child drives the
!   REAL hessian_cd into the stop and must print a naming message on stderr
!   and exit nonzero. A child that returns without stopping exits 0 and
!   fails the criterion.
! The test potentials are pure mathematics - no system knowledge enters.
! Criteria:
!   (a) diagonal accuracy, explicit step: V = 1/2 sum k_i q_i^2 (n = 6,
!       distinct k_i) at a generic nonzero point, caller step h = 5e-4 -
!       every diagonal element within 1e-8 RELATIVE of k_i, every
!       off-diagonal element below 1e-9 absolute
!   (b) off-diagonal coupling, DEFAULT step: V = 1/2 k(x^2+y^2) + c xy
!       + 1/2 k3 z^2 (n = 3) at a generic point - H = [[k,c,0],[c,k,0],
!       [0,0,k3]] within 1e-8 relative (nonzero elements) / below 1e-9
!       absolute (zero elements), and the returned matrix is EXACTLY
!       symmetric (bitwise H = H^T)
!   (c) quartic convergence order: V = x^4 + y^4, H = diag(12x^2, 12y^2)
!       at two generic points; halving the step (1e-1 -> 5e-2 -> 2.5e-2)
!       quarters the error - every consecutive error ratio within 4 +/-
!       0.05 (second-order central differences; the analytic error is
!       exactly 4 h^2, so a first-order scheme would ratio 2, a zeroth-
!       order one 1)
!   (d) protection loud-stops (children):
!     (d1) force returns NaN in component 1 exactly at the coordinate-2
!          probe -> nonzero exit, stderr names hessian_cd, coordinate 2,
!          component 1 and the +h probe
!     (d2) force returns +Inf in component 2 exactly at the coordinate-1
!          MINUS probe -> nonzero exit, stderr names the -h probe and
!          coordinate 1
!     (d3) hess shape 3 x 3 against q length 4 -> nonzero exit, stderr
!          names hessian_cd and both sizes
!     (d4) step h = 0 -> nonzero exit, stderr names the step requirement
!=====================================================================
module hessian_case
   implicit none
   private
   public :: force_harm_diag, force_coupled, force_quartic, force_poison_nan, &
             force_poison_inf, poison_setup, k_diag, k_cpl, c_cpl, k3_cpl
   real(8), parameter :: k_diag(6) = (/ 0.5d0, 1.25d0, 2.0d0, 0.75d0, 1.5d0, 3.0d0 /)
   real(8), parameter :: k_cpl = 1.4d0, c_cpl = 0.6d0, k3_cpl = 2.2d0
   ! poison plumbing: the poison force procs return finite forces everywhere
   ! except at one chosen probe, recognized by displacement from the base point
   real(8) :: poison_base(2) = 0.0d0
   real(8) :: poison_h = 1.0d-3
   real(8) :: val_zero = 0.0d0, val_one = 1.0d0
contains
   ! V = 1/2 sum k_i q_i^2  ->  F_i = -k_i q_i (any n <= 6)
   subroutine force_harm_diag(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      f = -k_diag(1:size(q))*q
   end subroutine force_harm_diag

   ! V = 1/2 k(x^2+y^2) + c xy + 1/2 k3 z^2  ->  F = (-(k x + c y), -(k y + c x), -k3 z)
   subroutine force_coupled(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      f(1) = -(k_cpl*q(1) + c_cpl*q(2))
      f(2) = -(k_cpl*q(2) + c_cpl*q(1))
      f(3) = -k3_cpl*q(3)
   end subroutine force_coupled

   ! V = sum q_i^4  ->  F_i = -4 q_i^3 (any n)
   subroutine force_quartic(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      f = -4.0d0*q**3
   end subroutine force_quartic

   ! finite harmonic forces EXCEPT: NaN in component 1 when coordinate 2 is
   ! displaced from the base point (either sign of the coordinate-2 probe)
   subroutine force_poison_nan(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      f = -k_cpl*q
      if (abs(q(2) - poison_base(2)) > 0.5d0*poison_h) f(1) = val_zero/val_zero
   end subroutine force_poison_nan

   ! finite harmonic forces EXCEPT: +Inf in component 2 when coordinate 1 is
   ! displaced BELOW the base point (the -h probe of coordinate 1)
   subroutine force_poison_inf(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      f = -k_cpl*q
      if (q(1) < poison_base(1)) f(2) = val_one/val_zero
   end subroutine force_poison_inf

   ! arm the poison recognition (call right before hessian_cd with the same
   ! base point and the same step hessian_cd will use)
   subroutine poison_setup(base, h_step)
      real(8), intent(in) :: base(2), h_step
      poison_base = base
      poison_h = h_step
   end subroutine poison_setup
end module hessian_case

program check_hessian
   use hessian, only: hessian_cd
   use hessian_case
   implicit none
   real(8), allocatable :: hess(:,:), hess_h1(:,:), hess_h2(:,:), hess_h3(:,:)
   real(8) :: q0(6), q4(4), qc(3), qp(2)
   real(8) :: h_ref(3,3), pts(2,2)
   real(8) :: err_seq(3), ex_err(3), ratio, worst_dev, worst_rel, worst_abs, worst_err, asym
   real(8), parameter :: h_seq(3) = (/ 1.0d-1, 5.0d-2, 2.5d-2 /)
   integer :: nfail = 0, nfail0, i, j, p
   character(len=256) :: mode

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('bad-nan')     ! (d1): NaN in component 1 at the coordinate-2 probe
         qp = (/ 0.4d0, -0.9d0 /)
         call poison_setup(qp, 1.0d-3)
         allocate (hess(2,2))
         call hessian_cd(force_poison_nan, qp, hess, 1.0d-3)
      case ('bad-inf')     ! (d2): +Inf in component 2 at the coordinate-1 minus probe
         qp = (/ 0.4d0, -0.9d0 /)
         call poison_setup(qp, 1.0d-3)
         allocate (hess(2,2))
         call hessian_cd(force_poison_inf, qp, hess, 1.0d-3)
      case ('bad-shape')   ! (d3): 3 x 3 output against a length-4 q
         q4 = (/ 0.1d0, -0.7d0, 0.9d0, 0.3d0 /)
         allocate (hess(3,3))
         call hessian_cd(force_quartic, q4, hess)
      case ('bad-step')    ! (d4): step h = 0
         qp = (/ 0.4d0, -0.9d0 /)
         allocate (hess(2,2))
         call hessian_cd(force_quartic, qp, hess, 0.0d0)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (hessian_cd returned - bug)'
      stop 0
   end if

   ! ---------- (a) diagonal accuracy, explicit step ----------
   nfail0 = nfail
   q0 = (/ 0.1d0, -0.7d0, 0.9d0, 0.3d0, -1.2d0, 0.55d0 /)
   allocate (hess(6,6))
   call hessian_cd(force_harm_diag, q0, hess, 5.0d-4)
   worst_rel = 0.0d0
   worst_abs = 0.0d0
   do i = 1, 6
      do j = 1, 6
         if (i == j) then
            worst_rel = max(worst_rel, abs(hess(i,i) - k_diag(i))/k_diag(i))
         else
            worst_abs = max(worst_abs, abs(hess(i,j)))
         end if
      end do
   end do
   print '(a,es10.3,a,es10.3)', '(a) harmonic diagonal: worst rel = ', worst_rel, &
      '  off-diag abs = ', worst_abs
   if (worst_rel >= 1.0d-8 .or. worst_abs >= 1.0d-9) then
      print '(a)', 'FAIL (a): diagonal accuracy out of tolerance'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(a) harmonic diagonal: PASS (6 distinct k_i within 1e-8 rel at h = 5e-4, off-diagonals < 1e-9)'

   ! ---------- (b) off-diagonal coupling, DEFAULT step ----------
   nfail0 = nfail
   deallocate (hess)
   allocate (hess(3,3))
   qc = (/ 0.4d0, -0.9d0, 1.1d0 /)
   call hessian_cd(force_coupled, qc, hess)
   h_ref = reshape( (/ k_cpl, c_cpl, 0.0d0, c_cpl, k_cpl, 0.0d0, 0.0d0, 0.0d0, k3_cpl /), (/ 3, 3 /) )
   worst_rel = 0.0d0
   worst_abs = 0.0d0
   do j = 1, 3
      do i = 1, 3
         if (h_ref(i,j) /= 0.0d0) then
            worst_rel = max(worst_rel, abs(hess(i,j) - h_ref(i,j))/abs(h_ref(i,j)))
         else
            worst_abs = max(worst_abs, abs(hess(i,j)))
         end if
      end do
   end do
   asym = maxval(abs(hess - transpose(hess)))
   print '(a,es10.3,a,es10.3,a,es10.3)', '(b) coupling: worst rel = ', worst_rel, &
      '  zero-elem abs = ', worst_abs, '  asymmetry = ', asym
   if (worst_rel >= 1.0d-8 .or. worst_abs >= 1.0d-9 .or. asym /= 0.0d0) then
      print '(a)', 'FAIL (b): coupling accuracy or symmetry out of tolerance'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(b) coupling: PASS (H(1,2) = c within 1e-8 rel, zero block < 1e-9, bitwise symmetric, default step)'

   ! ---------- (c) quartic convergence order ----------
   nfail0 = nfail
   pts(1,:) = (/ 0.3d0, 1.7d0 /)
   pts(2,:) = (/ 2.9d0, -1.1d0 /)
   allocate (hess_h1(2,2), hess_h2(2,2), hess_h3(2,2))
   worst_dev = 0.0d0
   worst_err = 0.0d0
   ex_err = 0.0d0
   do p = 1, 2
      qp = pts(p,:)
      call hessian_cd(force_quartic, qp, hess_h1, h_seq(1))
      call hessian_cd(force_quartic, qp, hess_h2, h_seq(2))
      call hessian_cd(force_quartic, qp, hess_h3, h_seq(3))
      do i = 1, 2
         err_seq(1) = abs(hess_h1(i,i) - 12.0d0*qp(i)**2)
         err_seq(2) = abs(hess_h2(i,i) - 12.0d0*qp(i)**2)
         err_seq(3) = abs(hess_h3(i,i) - 12.0d0*qp(i)**2)
         if (p == 1 .and. i == 1) ex_err = err_seq
         ratio = err_seq(1)/err_seq(2)
         worst_dev = max(worst_dev, abs(ratio - 4.0d0))
         ratio = err_seq(2)/err_seq(3)
         worst_dev = max(worst_dev, abs(ratio - 4.0d0))
         worst_err = max(worst_err, err_seq(3))
      end do
   end do
   print '(a,3es10.3)', '(c) quartic order: sample error decay (point 1, component 1) ', ex_err
   print '(a,es10.3,a,es10.3)', '(c) quartic order: worst |ratio-4| = ', worst_dev, &
      '  smallest-h error = ', worst_err
   if (worst_dev > 0.05d0) then
      print '(a)', 'FAIL (c): halving the step does not quarter the error'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(c) quartic order: PASS (h halved 1e-1 -> 5e-2 -> 2.5e-2, error ratio 4 at 2 points x 2 components)'

   ! ---------- (d) protection loud-stops (children) ----------
   call expect_loudstop('bad-nan',   'hessian_cd', 'coordinate 2', 'component 1', '+h probe', &
                        '(d1) NaN-force loud-stop')
   call expect_loudstop('bad-inf',   'hessian_cd', 'coordinate 1', 'component 2', '-h probe', &
                        '(d2) Inf-force loud-stop')
   call expect_loudstop('bad-shape', 'hessian_cd', 'shape 3 x 3',  'length 4',    '', &
                        '(d3) shape-mismatch loud-stop')
   call expect_loudstop('bad-step',  'hessian_cd', 'step',         'positive',    '', &
                        '(d4) non-positive-step loud-stop')

   print '(a,i0,a)', '=== check_hessian: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! spawn the child for one bad case and assert the loud-stop: nonzero exit,
   ! stderr names the expected authority (hessian_cd) and carries the expected
   ! key(s); key2/key3 = '' skip that key
   subroutine expect_loudstop(cmode, who, key1, key2, key3, label)
      character(len=*), intent(in) :: cmode, who, key1, key2, key3, label
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      ! the shell redirect needs native separators on Windows (forward slashes
      ! make cmd misparse the redirect target)
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_hessian.exe '//trim(cmode)//' 2>build'//sep//'hessian_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/hessian_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) line
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(line)
         end do
         close (v)
      end if
      if (stat == 0) then
         print '(a)', 'FAIL '//label//': expected a loud-stop, the child exited 0'
         nfail = nfail + 1
      else if (index(buf, who) == 0) then
         print '(a)', 'FAIL '//label//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (index(buf, trim(key1)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key1)//'"'
         nfail = nfail + 1
      else if (len_trim(key2) > 0 .and. index(buf, trim(key2)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key2)//'"'
         nfail = nfail + 1
      else if (len_trim(key3) > 0 .and. index(buf, trim(key3)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key3)//'"'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop

   ! Windows-host probe (same runtime OS branch the sibling checks use)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_hessian
