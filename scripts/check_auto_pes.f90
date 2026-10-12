!=====================================================================
! check_auto_pes.f90 - the auto_sample fixture's criterion program
!   (criteria defined FIRST in this header; this program executes them)
! Criteria:
!   (a) pes vs closed form over a battery (nat=1..4, k = the staged
!       default 5.0 eV/A^2): v = 0.5*k*sum(q^2), g = k*q,
!       max |d| <= 1e-12 (exact arithmetic)
!   (b) pes vs central-difference gradient (h = 1e-5 A, nat=3 fixed
!       battery): max |g - g_fd| <= 1e-8 eV/A (truncation ~ k*h^2/6)
!   (c) keys() is exactly ['K_SPRING', 'R_TERM'] (count + strings)
!   (d) full assignment: g/v prefilled with huge() come back finite
!   (e) method_name() == 'auto_tdhf' and is lowercase
! Compile (cwd = scripts/):
!   mkdir -p build && gfortran -c -Jbuild -o build/state.o ../src/globals/state.f90
!   gfortran -c -Jbuild -Ibuild -o build/auto_sample_pot.o ../systems/inactive/auto_sample/pot.f90
!   gfortran -c -Jbuild -Ibuild -o build/auto_sample_met.o ../systems/inactive/auto_sample/tdhf_stub.f90
!   gfortran -c -Jbuild -Ibuild -o build/check_auto_pes.o check_auto_pes.f90
!   gfortran -o build/check_auto_pes.exe build/state.o build/auto_sample_pot.o build/auto_sample_met.o build/check_auto_pes.o
!=====================================================================
program check_auto_pes
   use auto_sample_pot, only: pes, keys
   use auto_sample_met,  only: method_name
   implicit none
   integer :: nfail
   nfail = 0
   call crit_a; call crit_b; call crit_c; call crit_d; call crit_e
   if (nfail == 0) then
      print '(a)', 'check_auto_pes: ALL PASS'
   else
      print '(a,i0,a)', 'check_auto_pes: ', nfail, ' FAILURES'
      stop 1
   end if
contains
   subroutine crit_a
      integer, parameter :: nats(4) = [1, 2, 3, 4]
      real(8), parameter :: k_ref = 5.0d0   ! the staged default (criterion constant)
      real(8), allocatable :: q(:), m(:), g(:)
      real(8) :: v, dmax
      integer :: it, i
      dmax = 0.0d0
      do it = 1, size(nats)
         allocate (q(3*nats(it)), m(nats(it)), g(3*nats(it)))
         q = [(0.37d0*real(i, 8) + 0.11d0*real(it, 8), i = 1, 3*nats(it))]
         m = real(nats(it), 8)
         call pes(q, m, g, v)
         dmax = max(dmax, abs(v - 0.5d0*k_ref*dot_product(q, q)), &
                    maxval(abs(g - k_ref*q)))
         deallocate (q, m, g)
      end do
      if (dmax <= 1.0d-12) then
         print '(a,es12.3)', 'PASS: (a) closed form, max|d| = ', dmax
      else
         nfail = nfail + 1; print '(a,es12.3)', 'FAIL: (a) closed form, max|d| = ', dmax
      end if
   end subroutine crit_a
   subroutine crit_b
      real(8), parameter :: h = 1.0d-5
      real(8) :: q(9), m(3), g(9), gp(9), gd(9), gfd(9), v, v_hi, v_lo, dmax
      integer :: i
      q = [0.51d0, -1.32d0, 0.74d0, 2.05d0, -0.18d0, 1.44d0, -0.92d0, 0.33d0, 1.11d0]
      m = 1.0d0
      call pes(q, m, g, v)
      do i = 1, 9
         q(i) = q(i) + h;  call pes(q, m, gp, v);  v_hi = v
         q(i) = q(i) - 2*h; call pes(q, m, gd, v);  v_lo = v
         q(i) = q(i) + h
         gfd(i) = (v_hi - v_lo)/(2*h)           ! g = dV/dq (the positive gradient;
                                                 !  the interface folds F = -g)
      end do
      dmax = maxval(abs(g - gfd))
      if (dmax <= 1.0d-8) then
         print '(a,es12.3)', 'PASS: (b) central difference, max|d| = ', dmax
      else
         nfail = nfail + 1; print '(a,es12.3)', 'FAIL: (b) central difference, max|d| = ', dmax
      end if
   end subroutine crit_b
   subroutine crit_c
      character(len=16) :: kk(2)
      kk = keys()
      if (trim(kk(1)) == 'K_SPRING' .and. trim(kk(2)) == 'R_TERM') then
         print '(a)', 'PASS: (c) keys() vocabulary exact'
      else
         nfail = nfail + 1; print '(a)', 'FAIL: (c) keys() vocabulary'
      end if
   end subroutine crit_c
   subroutine crit_d
      real(8) :: q(3), m(1), g(3), v
      q = [0.1d0, 0.2d0, 0.3d0]; m = 1.0d0
      g = huge(1.0d0); v = huge(1.0d0)         ! garbage prefill
      call pes(q, m, g, v)
      if (all(abs(g) < huge(1.0d0)) .and. abs(v) < huge(1.0d0)) then
         print '(a)', 'PASS: (d) full assignment (garbage overwritten)'
      else
         nfail = nfail + 1; print '(a)', 'FAIL: (d) full assignment'
      end if
   end subroutine crit_d
   subroutine crit_e
      character(len=32) :: s, lows
      s = method_name()
      lows = s
      call to_lower(lows)
      if (trim(s) == 'auto_tdhf' .and. s == lows) then
         print '(a)', 'PASS: (e) method_name() == ''auto_tdhf'' (lowercase)'
      else
         nfail = nfail + 1; print '(a)', 'FAIL: (e) method_name()'
      end if
   end subroutine crit_e
   subroutine to_lower(str)
      character(len=*), intent(inout) :: str
      integer :: i, ic
      do i = 1, len_trim(str)
         ic = iachar(str(i:i))
         if (ic >= iachar('A') .and. ic <= iachar('Z')) str(i:i) = achar(ic + 32)
      end do
   end subroutine to_lower
end program check_auto_pes
