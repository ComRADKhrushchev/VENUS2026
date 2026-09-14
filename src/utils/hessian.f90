!=====================================================================
! hessian.f90 - L5 instrument: symmetric Hessian of a potential assembled
!               from a force procedure by second-order central differences
! Responsibility: given a position vector q and a procedure that returns
!   forces at displaced positions, build the symmetric second-derivative
!   matrix of the underlying potential. Pure form: no system knowledge, no
!   facade, no state - the module owns no units of its own; q, f, h and
!   hess all carry the caller's units (h is an increment in the units of
!   q; hess carries [force unit]/[length unit]).
! Contract:
!   force      - procedure conforming to the force_proc interface returning
!                f = -dV/dq at q (force is the NEGATIVE gradient). The
!                output satisfies H(i,j) = d2V/dq_i dq_j = -dF_i/dq_j,
!                assembled column-wise as the central difference
!                    H_raw(:,j) = -( F(q + h e_j) - F(q - h e_j) ) / (2h)
!                (2n force calls: coordinates ascending, the plus probe
!                then the minus probe per coordinate; no call at the bare
!                base point).
!   step       - ONE uniform scalar step for every coordinate (not
!                per-coordinate). Default 1e-4 in the units of q - the
!                double-precision balance for O(1) coordinates (truncation
!                O(h^2) ~ 1e-8 against subtractive roundoff O(eps/h) ~
!                1e-12); the optional h argument overrides it for other
!                coordinate scales. h must be positive and finite.
!   output     - hess is n x n with n = size(q) (3N x 3N for N atoms in
!                the downstream Cartesian use), fully defined on exit,
!                and SYMMETRIZED: H = (H_raw + H_raw^T)/2 - the raw
!                central difference of an exact force is symmetric only
!                up to rounding; the average is exactly symmetric by
!                construction.
! Error posture (the sibling L5 instrument precedent, geometry): an
!   invalid condition writes a naming message to stderr (unit 0) and
!   stops with exit code 1 - hess not n x n (both sizes named), step h
!   non-positive or non-finite, or a NaN/Inf force component at a probe
!   (component, probe sign and coordinate named).
!=====================================================================
module hessian
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: force_proc, hessian_cd

   ! the force procedure contract: f = -dV/dq at q (force is the NEGATIVE
   ! gradient); both arrays carry the caller's units, both of length n
   abstract interface
      subroutine force_proc(q, f)
         real(8), intent(in)  :: q(:)   ! position vector (n = size(q))
         real(8), intent(out) :: f(:)   ! forces at q
      end subroutine force_proc
   end interface

   ! default uniform step: an increment in the units of q (see header)
   real(8), parameter :: h_default = 1.0d-4
contains
   !------------------------------------------------------------------
   ! hessian_cd(force, q, hess, h) - central-difference Hessian over a force
   ! procedure, symmetrized (see the module header for the full contract)
   subroutine hessian_cd(force, q, hess, h)
      procedure(force_proc)         :: force    ! force provider, f = -dV/dq
      real(8), intent(in)           :: q(:)     ! base position vector
      real(8), intent(out)          :: hess(:,:)! symmetric Hessian, n x n
      real(8), intent(in), optional :: h        ! uniform step override (units of q)
      real(8) :: h_step
      real(8) :: q_probe(size(q)), f_plus(size(q)), f_minus(size(q))
      integer :: n, j
      n = size(q)
      if (size(hess,1) /= n .or. size(hess,2) /= n) then
         write (0, '(a,i0,a,i0,a,i0)') 'hessian_cd: hess shape ', size(hess,1), ' x ', &
            size(hess,2), ' does not match q length ', n
         stop 1
      end if
      h_step = h_default
      if (present(h)) then
         if (.not. ieee_is_finite(h) .or. h <= 0.0d0) then
            write (0, '(a,es13.5)') 'hessian_cd: step h must be positive and finite, got ', h
            stop 1
         end if
         h_step = h
      end if
      do j = 1, n
         q_probe = q
         q_probe(j) = q_probe(j) + h_step
         call force(q_probe, f_plus)
         call require_finite(f_plus, j, plus=.true.)
         q_probe = q
         q_probe(j) = q_probe(j) - h_step
         call force(q_probe, f_minus)
         call require_finite(f_minus, j, plus=.false.)
         ! H(i,j) = -dF_i/dq_j (force is the negative gradient)
         hess(:,j) = -(f_plus - f_minus)/(2.0d0*h_step)
      end do
      ! symmetrize in place: H <- (H + H^T)/2 (see the module header)
      hess = 0.5d0*(hess + transpose(hess))
   end subroutine hessian_cd

   !------------------------------------------------------------------
   ! require_finite(f, j, plus) - probe-result gate: every component finite,
   ! else loud-stop naming the component, the probe sign and the coordinate
   subroutine require_finite(f, j, plus)
      real(8), intent(in) :: f(:)     ! force vector returned by one probe
      integer, intent(in) :: j        ! displaced coordinate index
      logical, intent(in) :: plus     ! .true. = the +h probe, .false. = the -h probe
      character(len=2) :: sgn
      integer :: i
      sgn = '-h'
      if (plus) sgn = '+h'
      do i = 1, size(f)
         if (.not. ieee_is_finite(f(i))) then
            write (0, '(a,i0,a,a,a,i0)') 'hessian_cd: force component ', i, &
               ' is non-finite at the ', sgn, ' probe of coordinate ', j
            stop 1
         end if
      end do
   end subroutine require_finite
end module hessian
