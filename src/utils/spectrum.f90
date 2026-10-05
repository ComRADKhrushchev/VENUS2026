!=====================================================================
! spectrum.f90 - fragment normal-mode spectrum from a Hessian block:
!                mass-weighting, eigensolve, zero-mode separation
! Design:
!   Pure mathematics over the internal unit family (amu / A / 10 fs):
!   with the force folded to internal-per-A (the force interface's unique
!   fold), the mass-weighted Hessian eigenvalue IS omega^2 in
!   [rad/(10 fs)] - no conversion literal lives here. The wavenumber
!   face goes through the shared consts pair e2wvn*hbar_code/e_conv
!   (numerically 2*pi*c*1e-14, the SI chain). Zero-mode count derives
!   from geometry linearity (principal moments at the COM: one
!   vanishing moment = linear, 3N-5 zeros; else 3N-6); the split is
!   validated by magnitude (a zero that is not small, or a vibrational
!   eigenvalue that is not positive - the probe point is not a minimum
!   - is a named abort, never a silent wrong table).
!=====================================================================
module spectrum
   use consts, only: e2wvn, hbar_code, e_conv
   use linalg, only: eig_sym
   implicit none
   private
   public :: spectrum_modes, stretch_mode, w_to_wvn, wvn_to_w

   ! split validation: the largest "zero" must sit this far below the
   ! smallest vibrational eigenvalue (relative). The buffered geometry is
   ! hand-written (4-6 significant digits), so a rounding-level placement
   ! error leaves a string-tension elevation of order 1e-5 on the zero
   ! modes - tolerated; a genuinely off-minimum buffer (a 0.16 A error
   ! on H2 sits at 0.1) still fails by two orders
   real(8), parameter :: zero_rel = 1.0d-3

   ! linear-fragment test: smallest/middle principal moment ratio
   real(8), parameter :: lin_tol = 1.0d-4
contains
   !------------------------------------------------------------------
   ! spectrum_modes(hess, mass, qz, w, c) - one fragment's vibrational
   !                 spectrum. hess: 3nat x 3nat Cartesian Hessian block
   !                 [internal/A^2] (destroyed); mass(nat) [amu];
   !                 qz(3nat) the fragment geometry [A] (linearity only).
   !                 Returns w(n_mode) [rad/(10 fs)] ascending and the
   !                 matching mass-weighted orthonormal columns
   !                 c(3nat, n_mode). A monoatomic fragment returns an
   !                 allocated-empty table (no evaluation)
   !------------------------------------------------------------------
   subroutine spectrum_modes(hess, mass, qz, w, c)
      real(8), intent(inout) :: hess(:,:)
      real(8), intent(in)    :: mass(:), qz(:)
      real(8), allocatable, intent(out) :: w(:), c(:,:)
      integer :: nat, n, n_zero, i
      real(8) :: d(3*size(mass)), sqmv(3*size(mass))
      nat = size(mass)
      n = 3*nat
      if (nat <= 1) then               ! monoatomic: no internal dof
         allocate (w(0), c(n, 0))
         return
      end if
      ! 1. mass-weight: D = H / sqrt(m_i m_j); eigenvectors come out in
      !    mass-weighted coordinates (the member table convention)
      do i = 1, nat
         sqmv(3*i-2:3*i) = sqrt(mass(i))
      end do
      hess = hess/spread(sqmv, 1, n)
      hess = hess/spread(sqmv, 2, n)
      call eig_sym(hess, d)              ! ascending eigenvalues, columns
      ! 2. zero-mode count from linearity (two atoms are always linear)
      if (nat == 2) then
         n_zero = 5
      else if (is_linear(mass, qz)) then
         n_zero = 5
      else
         n_zero = 6
      end if
      ! 3. split validation + vibrational set
      if (abs(d(n_zero)) > zero_rel*d(n_zero + 1)) then
         call stop_spec('the zero-mode/vibrational split is ambiguous: the largest zero '// &
                        'eigenvalue is not small against the smallest vibrational one - '// &
                        'the probe point is not a fragment minimum')
      end if
      if (d(n_zero + 1) <= 0.0d0) then
         call stop_spec('a vibrational eigenvalue is non-positive - the probe geometry is '// &
                        'not a minimum of the fragment potential')
      end if
      allocate (w(n - n_zero), c(n, n - n_zero))
      do i = 1, n - n_zero
         w(i) = sqrt(d(n_zero + i))
         c(:, i) = hess(:, n_zero + i)
      end do
   end subroutine spectrum_modes

   !------------------------------------------------------------------
   ! is_linear(mass, qz) - one vanishing principal moment of inertia at
   !                 the COM (eigenvalues of the inertia tensor)
   !------------------------------------------------------------------
   logical function is_linear(mass, qz)
      real(8), intent(in) :: mass(:), qz(:)
      real(8) :: inertia(3, 3), mom(3), rc(3), r(3)
      real(8) :: wt
      integer :: i, nat
      nat = size(mass)
      wt = sum(mass)
      rc = 0.0d0
      do i = 1, nat
         rc = rc + mass(i)*qz(3*i-2:3*i)
      end do
      rc = rc/wt
      inertia = 0.0d0
      do i = 1, nat
         r = qz(3*i-2:3*i) - rc
         inertia(1, 1) = inertia(1, 1) + mass(i)*(r(2)**2 + r(3)**2)
         inertia(2, 2) = inertia(2, 2) + mass(i)*(r(1)**2 + r(3)**2)
         inertia(3, 3) = inertia(3, 3) + mass(i)*(r(1)**2 + r(2)**2)
         inertia(1, 2) = inertia(1, 2) - mass(i)*r(1)*r(2)
         inertia(2, 3) = inertia(2, 3) - mass(i)*r(2)*r(3)
         inertia(1, 3) = inertia(1, 3) - mass(i)*r(1)*r(3)
      end do
      inertia(2, 1) = inertia(1, 2)
      inertia(3, 2) = inertia(2, 3)
      inertia(3, 1) = inertia(1, 3)
      call eig_sym(inertia, mom)
      is_linear = (mom(1) <= lin_tol*mom(2))
   end function is_linear

   !------------------------------------------------------------------
   ! stretch_mode(qz, mass, c) - the diatomic stretch column aligned with
   !                 the buffered bond: unit mass-weighted eigenvector of
   !                 the relative coordinate
   !------------------------------------------------------------------
   subroutine stretch_mode(qz, mass, c)
      real(8), intent(in)  :: qz(6)      ! the two atom positions [A]
      real(8), intent(in)  :: mass(2)    ! [amu]
      real(8), intent(out) :: c(6)       ! unit mass-weighted column [-]
      real(8) :: u(3), mt, mu
      mt = mass(1) + mass(2)
      mu = mass(1)*mass(2)/mt
      u = qz(4:6) - qz(1:3)
      u = u/norm2(u)
      c(1:3) = -u*sqrt(mass(1))*mass(2)/(mt*sqrt(mu))
      c(4:6) = u*sqrt(mass(2))*mass(1)/(mt*sqrt(mu))
   end subroutine stretch_mode

   ! w_to_wvn(w) / wvn_to_w(nu) - internal angular frequency [rad/(10 fs)]
   !                 <-> wavenumber [cm^-1] through the shared consts
   pure function w_to_wvn(w) result(nu)
      real(8), intent(in) :: w
      real(8) :: nu
      nu = e2wvn*hbar_code*w/e_conv
   end function w_to_wvn

   pure function wvn_to_w(nu) result(w)
      real(8), intent(in) :: nu
      real(8) :: w
      w = nu*e_conv/(e2wvn*hbar_code)
   end function wvn_to_w

   ! stop_spec(msg) - the named-abort channel
   subroutine stop_spec(msg)
      character(len=*), intent(in) :: msg
      write (0, '(a)') 'spectrum_modes: '//trim(msg)
      write (0, '(a)') 'spectrum_modes: fatal (the spectrum table is unusable)'
      stop 1
   end subroutine stop_spec
end module spectrum
