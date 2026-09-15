!=====================================================================
! check_density.f90 - criteria program for the quantum state density
!                      dens_states (Beyer-Swinehart direct counting inside a
!                      central-difference shell) and the formatted eigenvalue
!                      output eig_out (criteria defined FIRST in this header;
!                      this program executes them)
! Compile (from src/; -Jbuild keeps the .mod files out of the source tree):
!   gfortran -Jbuild output/density.f90 ../scripts/check_density.f90 -o build/check_density.exe
! Run (from src/): build\check_density.exe
!
! Anchor derivations:
!   The counting core quantizes every frequency to an effective integer
!   ne = int(w + 0.499999) and counts on unit bins; at INTEGER probe
!   energies the counted array holds exact small integers (double adds of
!   integer values stay exact far below 2^53), so the shell's two probe
!   evaluations (E +/- 50 bins) and the 100-bin quotient are exact
!   arithmetic. The independent anchors are brute-force quantum-number
!   enumerations (n_exact: recursive product over the mode quantum
!   numbers, counting the combinations whose energy does not exceed the
!   probe) - a different algorithm than the recursion under test, never a
!   re-run of it. Window edge cases ride the same enumeration: E - 50 < 0
!   contributes no states (the early empty return), E < -50 empties both
!   probes and leaves rho bitwise zero.
!   The eigenvalue contract: values and vector components are chosen
!   exactly representable in binary with at most five significant decimal
!   digits, so the 6-column block print round-trips them bitwise; the
!   block structure is checked line by line. Blank-line layout (the
!   contract of the block formats): TWO blank lines open a block (the
!   leading double slash), and TWO blank lines follow the eigenvalue line
!   (the lone-slash record: the slash terminates an empty record and the
!   statement end terminates the record opened after it - a conforming
!   processor writes both).
!
! Tolerances:
!   - state-density scalars: rel 1e-10 with floor 1.0 (integer-valued
!     quantities; a wrong recursion, ordering, or interpolation deviates
!     at the O(1) level).
!   - eigenvalue/component round-trip: bitwise (exact-binary anchors);
!     the matrix restoration bitwise.
!
! Criteria:
!   (a) single- and few-mode harmonic counting: modes [1], [3], [5] and
!       the triple [2,3,5] at several integer probes, each against the
!       enumeration anchor (interior windows and windows clipped at zero
!       below); the unit frequency gives rho = 1 exactly; a probe below
!       -50 leaves rho bitwise zero (both probes empty)
!   (b) multi-mode combination density: the coprime pair [7,11] (no
!       common period below the window scale) at four probes including
!       one below the first combined excitation; two pairs of
!       non-integer frequencies quantized by the counting core to
!       effective integers - the anchors enumerate the EFFECTIVE
!       integers (registered reading: a genuinely irrational ratio has
!       no exact anchor against the counting because the core quantizes
!       first; the pairs exercise the quantization branch and are
!       anchored on their effective image)
!   (c) the eig_out output contract: an 8x8 table (two 6-column blocks,
!       6+2) with known eigenvalues and vectors -> the written file
!       parsed back: block structure (blank-pair / column header /
!       eigenvalue line / blank-pair / rows), bitwise value round-trip,
!       bitwise matrix restoration
! Registered notes (part of this check's contract):
!   - pointwise anchors sit at INTEGER probe energies (the shell's
!     linear interpolation between integer bins is the counting's own
!     interpolation - anchoring it against itself would be circular);
!   - the window half-width 50 bins and the bin width 1.0 are the frozen
!     constants of the differencing shell's contract (in-bin units).
!=====================================================================
program check_density
   use density, only: dens_states, eig_out
   implicit none
   integer :: nfail = 0
   integer :: nfail0                     ! per-criterion failure watermark

   call crit_a()      ! (a) single/few-mode counts vs enumeration
   call crit_b()      ! (b) multi-mode combination density
   call crit_c()      ! (c) the eig_out output contract

   if (nfail == 0) then
      print '(a)', 'check_density: ALL PASS'
   else
      print '(a,i0,a)', 'check_density: ', nfail, ' FAILURE(S)'
      stop 1
   end if
contains
   !------------------------------------------------------------------
   ! crit_a() - (a): single- and few-mode counting against enumeration
   subroutine crit_a()
      real(8) :: rho
      nfail0 = nfail
      call rho_chk('(a) w=[1] E=120 (interior window)', [1.0d0], 120.0d0)
      call rho_chk('(a) w=[1] E=40 (clipped below zero)', [1.0d0], 40.0d0)
      call rho_chk('(a) w=[3] E=97 (interior window)', [3.0d0], 97.0d0)
      call rho_chk('(a) w=[3] E=40 (clipped below zero)', [3.0d0], 40.0d0)
      call rho_chk('(a) w=[5] E=200 (interior window)', [5.0d0], 200.0d0)
      call rho_chk('(a) w=[2,3,5] E=100 (few-mode)', [2.0d0, 3.0d0, 5.0d0], 100.0d0)
      call rho_chk('(a) w=[2,3,5] E=37 (few-mode, low)', [2.0d0, 3.0d0, 5.0d0], 37.0d0)
      call dens_states([3.0d0], -60.0d0, rho)
      call chk_z('(a) E=-60 (both probes empty)', rho)
      if (nfail == nfail0) then
         print '(a)', '(a): PASS (single/few-mode counts on the enumeration anchor at'// &
            ' integer probes; the zero window below -50 is bitwise zero)'
      end if
   end subroutine crit_a

   !------------------------------------------------------------------
   ! crit_b() - (b): multi-mode combination density
   subroutine crit_b()
      nfail0 = nfail
      call rho_chk('(b) w=[7,11] E=55', [7.0d0, 11.0d0], 55.0d0)
      call rho_chk('(b) w=[7,11] E=132', [7.0d0, 11.0d0], 132.0d0)
      call rho_chk('(b) w=[7,11] E=210', [7.0d0, 11.0d0], 210.0d0)
      call rho_chk('(b) w=[7,11] E=8 (below the first combined level)', [7.0d0, 11.0d0], 8.0d0)
      call rho_chk('(b) w=[1,sqrt(2)] E=83 (quantized pair)', &
                   [1.0d0, 1.4142135623730951d0], 83.0d0)
      call rho_chk('(b) w=[2.4,3.6] E=90 (quantized pair)', [2.4d0, 3.6d0], 90.0d0)
      if (nfail == nfail0) then
         print '(a)', '(b): PASS (the coprime pair across four probes; non-integer'// &
            ' pairs anchored on their effective integers - the quantization branch)'
      end if
   end subroutine crit_b

   !------------------------------------------------------------------
   ! crit_c() - (c): the eig_out output contract (8x8 table, blocks 6+2)
   subroutine crit_c()
      real(8) :: a(8, 8), ac(8, 8), eig(8)
      real(8) :: ev(6), comp(6)
      real(8), parameter :: eig_vals(8) = [ 0.5d0, -1.0d0, 2.0d0, -0.25d0, &
                                           4.0d0, 10.0d0, -7.5d0, 0.125d0 ]
      integer :: cols(6), inum, ios
      character(len=132) :: lines(64)
      integer :: nlines, il, l, k, m, ir, i, u
      nfail0 = nfail
      eig = eig_vals
      do i = 1, 8                          ! vectors by column (the linalg convention)
         do m = 1, 8
            a(m, i) = 0.25d0*mod(i*m, 5)   ! exact-binary, at most 6 decimals
         end do
      end do
      ac = a
      open (newunit=u, file='build/dens_eig.txt', status='replace', action='write')
      call eig_out(a, eig, u)
      close (u)
      do i = 1, 8
         do m = 1, 8
            if (a(m, i) /= ac(m, i)) then
               call fail('(c) the matrix was not restored bitwise at ('//i2s(m)//','// &
                         i2s(i)//')')
            end if
         end do
      end do
      if (nfail /= nfail0) return
      call load_lines('build/dens_eig.txt', lines, nlines)
      il = 1
      l = 1
      do while (l <= 8)
         k = min(l + 5, 8)
         call expect_blank(lines, nlines, il, '(c) the two blank lines opening a block')
         call expect_blank(lines, nlines, il, '(c) the two blank lines opening a block')
         if (il > nlines) then
            call fail('(c) the file ended inside a block header')
            return
         end if
         read (lines(il), '(2x,6(9x,i3))', iostat=ios) cols
         if (ios /= 0) then
            call fail('(c) the column header does not parse: "'//trim(lines(il))//'"')
            return
         end if
         do m = 1, k - l + 1
            if (cols(m) /= l + m - 1) then
               call fail('(c) column header '//i2s(m)//': got '//i2s(cols(m))//', want '// &
                         i2s(l + m - 1))
            end if
         end do
         il = il + 1
         if (il > nlines) then
            call fail('(c) the file ended before the eigenvalue line')
            return
         end if
         read (lines(il), '(5x,6e12.4)') ev
         do m = 1, k - l + 1
            if (ev(m) /= eig(l + m - 1)) then
               call fail('(c) eigenvalue '//i2s(l + m - 1)//' round-trip: got '//d2s(ev(m))// &
                        ', want '//d2s(eig(l + m - 1)))
            end if
         end do
         il = il + 1
         call expect_blank(lines, nlines, il, '(c) the blank pair after the eigenvalues')
         call expect_blank(lines, nlines, il, '(c) the blank pair after the eigenvalues')
         do ir = 1, 8
            if (il > nlines) then
               call fail('(c) the file ended inside the row lines')
               return
            end if
            read (lines(il), '(i3,1x,6f12.6)') inum, comp
            if (inum /= ir) then
               call fail('(c) row number: got '//i2s(inum)//', want '//i2s(ir))
            end if
            do m = 1, k - l + 1
               if (comp(m) /= ac(l + m - 1, ir)) then
                  call fail('(c) component (col '//i2s(l + m - 1)//', row '//i2s(ir)// &
                            ') round-trip: got '//d2s(comp(m))//', want '// &
                            d2s(ac(l + m - 1, ir)))
               end if
            end do
            il = il + 1
         end do
         l = k + 1
      end do
      if (il /= nlines + 1) then
         call fail('(c) trailing lines after the last block: the file has '//i2s(nlines)// &
                   ' lines, the parse consumed '//i2s(il - 1))
      end if
      if (nfail == nfail0) then
         print '(a)', '(c): PASS (two 6-column blocks: blank-pair/header/eigenvalue/'// &
            'blank/rows structure; values round-trip bitwise; the matrix restored bitwise)'
      end if
   end subroutine crit_c

   !------------------------------------------------------------------
   ! rho_chk - one dens_states call against the enumeration anchor
   subroutine rho_chk(label, w, e)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: w(:), e
      real(8) :: rho, want
      want = (dble(n_exact(eff_int(w), e + 50.0d0)) - &
              dble(n_exact(eff_int(w), e - 50.0d0)))/100.0d0
      call dens_states(w, e, rho)
      call chk_r(label, rho, want)
   end subroutine rho_chk

   ! eff_int - the counting core's frequency quantization (anchor side)
   pure function eff_int(w) result(ne)
      real(8), intent(in) :: w(:)
      integer :: ne(size(w))
      ne = int(w + 0.499999d0)
   end function eff_int

   ! n_exact - exact state count by quantum-number enumeration
   recursive function n_exact(ne, e) result(n)
      integer, intent(in) :: ne(:)          ! effective integer frequencies [bin]
      real(8), intent(in) :: e              ! probe energy [bin]
      integer :: n, j
      n = 0
      if (e < 0.0d0) return
      if (size(ne) == 0) then
         n = 1
         return
      end if
      do j = 0, int(e)/ne(1)
         n = n + n_exact(ne(2:), e - dble(j*ne(1)))
      end do
   end function n_exact

   !------------------------------------------------------------------
   ! assertion helpers (watermark discipline: chk_* record failures only)
   !------------------------------------------------------------------

   ! chk_r - relative scalar check, tol 1e-10 with floor 1.0 (header derivation)
   subroutine chk_r(label, got, want)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got, want
      if (abs(got - want) > 1.0d-10*max(abs(want), 1.0d0)) then
         call fail(trim(label)//': got '//d2s(got)//', want '//d2s(want))
      end if
   end subroutine chk_r

   ! chk_z - bitwise-zero scalar check (empty-window contract)
   subroutine chk_z(label, got)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got
      if (got /= 0.0d0) then
         call fail(trim(label)//': must be bitwise zero, got '//d2s(got))
      end if
   end subroutine chk_z

   ! expect_blank - one blank-structure line at the cursor
   subroutine expect_blank(lines, nlines, il, label)
      character(len=*), intent(in) :: lines(:), label
      integer, intent(in) :: nlines
      integer, intent(inout) :: il
      if (il > nlines) then
         call fail(trim(label)//': the file ended where a blank line was expected')
         return
      end if
      if (len_trim(lines(il)) /= 0) then
         call fail(trim(label)//': non-blank line "'//trim(lines(il))//'"')
      end if
      il = il + 1
   end subroutine expect_blank

   ! load_lines - one text file into a line table (missing file = zero lines)
   subroutine load_lines(fname, lines, n)
      character(len=*), intent(in) :: fname
      character(len=*), intent(out) :: lines(:)
      integer, intent(out) :: n
      integer :: u, ios
      n = 0
      open (newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do while (n < size(lines))
         read (u, '(a)', iostat=ios) lines(n + 1)
         if (ios /= 0) exit
         n = n + 1
      end do
      close (u)
   end subroutine load_lines

   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      print '(a)', 'FAIL '//trim(msg)
      nfail = nfail + 1
   end subroutine fail

   function i2s(i)
      integer, intent(in) :: i
      character(len=12) :: i2s
      write (i2s, '(i0)') i
   end function i2s
   function d2s(x)
      real(8), intent(in) :: x
      character(len=24) :: d2s
      write (d2s, '(es12.4)') x
   end function d2s
end program check_density
