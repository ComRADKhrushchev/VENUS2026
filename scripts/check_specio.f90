!=====================================================================
! check_specio.f90 - C4.1 acceptance criteria program (criteria defined FIRST in
!                    docs/plans/C4.md; this program executes them)
! Compile (from src/; module first, -Jbuild keeps the .mod out of the source tree):
!   gfortran -Jbuild utils/specio.f90 ../scripts/check_specio.f90 -o build/check_specio.exe
! All sample files are written under build/ by this program - the source tree stays clean.
! Criteria (C4.md C4.1):
!   (a) read_xyz clean sample: per-field names/coordinates check (exact values),
!       errmsg blank on the clean parse
!   (b) read_xyz trailing columns: extra fields after "Symbol x y z" are read past,
!       never interpreted - clean parse, coordinates intact
!   (c) elem_mass spot checks: H = 1.008, C = 12.011, Au = 196.97, He = 4.003
!       (IUPAC 2013 standard atomic weights at the quoted rounding, tol 5e-3);
!       case/blank normalization ('AU', ' c' resolve to the 'Au'/'C' rows)
!   (d) read_poscar clean Direct sample: title / scaled cell / counts-expanded names
!       (block order preserved) / fractional-to-Cartesian conversion, per field;
!       errmsg blank
!   (e) read_poscar Cartesian mode: the universal factor scales the cell AND the
!       coordinates
!   (f) read_poscar negative scale (target volume): |det(cell)| = |s| = 27 and the
!       cell rows are the raw vectors rescaled by (|s|/|det(raw vectors)|)^(1/3)
!   (g) read_poscar Selective-dynamics line skipped, trailing T/F flags read past,
!       the coordinate block stops before the velocity block (never read)
!   (h) bad samples -> ok = .false. + non-blank errmsg:
!       xyz    - non-numeric count / zero count / short coordinate row / non-numeric
!                coordinate field / symbol longer than 4 characters / file ends early
!       POSCAR - non-numeric scaling factor / zero scaling factor / symbols-counts
!                field mismatch / zero atoms / unknown mode tag / short coordinate row
!   (i) read_poscar skewed (non-diagonal) cell with Direct coordinates: the row-stored
!       lattice vectors a1=(1,0,0), a2=(0.5,2,0), a3=(0.2,0.3,1.5) combined as
!       r = 0.25*a1 + 0.50*a2 + 0.75*a3 = (0.65, 1.225, 1.125) - distinguishes
!       matmul(fr, cell) (correct linear combination) from matmul(cell, fr)
!       (row-contraction, coincides only for diagonal cells; fix-round-1 regression)
!=====================================================================
program check_specio
   use specio, only: read_poscar, read_xyz, elem_mass
   implicit none
   character(len=4), allocatable :: names(:)
   real(8), allocatable :: xyz(:,:)
   real(8) :: cell(3,3), res, detc, sexp
   character(len=128) :: title, emsg
   logical :: ok
   integer :: u, nfail = 0, nfail0
   real(8), parameter :: wxyz_a(3,3) = reshape( [ 0.10d0,0.20d0,0.30d0, &
      1.10d0,2.20d0,3.30d0, -1.00d0,0.50d0,2.00d0 ], [ 3,3 ])
   real(8), parameter :: wxyz_b(3,2) = reshape( [ 0.50d0,0.0d0,0.0d0, &
      -0.50d0,0.0d0,0.0d0 ], [ 3,2 ])
   real(8), parameter :: wcell_d(3,3) = reshape( [ 2.5d0,0.0d0,0.0d0, &
      0.0d0,3.0d0,0.0d0, 0.0d0,0.0d0,2.25d0 ], [ 3,3 ])
   real(8), parameter :: wxyz_d(3,3) = reshape( [ 0.0d0,0.0d0,0.0d0, &
      1.25d0,1.50d0,1.125d0, 0.625d0,2.25d0,0.225d0 ], [ 3,3 ])
   real(8), parameter :: wcell_e(3,3) = reshape( [ 2.0d0,0.0d0,0.0d0, &
      0.0d0,2.0d0,0.0d0, 0.0d0,0.0d0,2.0d0 ], [ 3,3 ])
   real(8), parameter :: wxyz_e(3,2) = reshape( [ 1.0d0,1.0d0,1.0d0, &
      2.0d0,4.0d0,6.0d0 ], [ 3,2 ])
   real(8), parameter :: wxyz_g(3,1) = reshape( [ 1.0d0,1.0d0,1.0d0 ], [ 3,1 ])
   real(8), parameter :: wxyz_i(3,1) = reshape( [ 0.65d0,1.225d0,1.125d0 ], [ 3,1 ])

   ! ---------- sample files (all written under build/) ----------
   call write_sample('build/specio_xyz_ok.xyz', [ character(len=64) :: &
      '3', 'c4.1 xyz clean sample', &
      'C 0.10 0.20 0.30', 'O 1.10 2.20 3.30', 'Au -1.00 0.50 2.00' ])
   call write_sample('build/specio_xyz_extra.xyz', [ character(len=64) :: &
      '2', 'c4.1 xyz trailing-columns sample (extra fields read past)', &
      'C 0.50 0.00 0.00 0.10 0.20 0.30', 'O -0.50 0.00 0.00 -0.10 -0.20 -0.30' ])
   call write_sample('build/specio_poscar_ok.poscar', [ character(len=64) :: &
      'c4.1 direct-mode sample', '2.5', &
      '1.0 0.0 0.0', '0.0 1.2 0.0', '0.0 0.0 0.9', &
      'Si O', '2 1', 'Direct', &
      '0.00 0.00 0.00', '0.50 0.50 0.50', '0.25 0.75 0.10' ])
   call write_sample('build/specio_poscar_cart.poscar', [ character(len=64) :: &
      'c4.1 cartesian-mode sample', '2.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'H', '2', 'Cartesian', &
      '0.5 0.5 0.5', '1.0 2.0 3.0' ])
   call write_sample('build/specio_poscar_vol.poscar', [ character(len=64) :: &
      'c4.1 negative-scale (target-volume) sample', '-27.0', &
      '1.0 0.0 0.0', '0.0 2.0 0.0', '0.0 0.0 1.5', &
      'H', '1', 'Direct', &
      '0.10 0.20 0.30' ])
   call write_sample('build/specio_poscar_sel.poscar', [ character(len=64) :: &
      'c4.1 selective-dynamics sample', '1.0', &
      '2.0 0.0 0.0', '0.0 2.0 0.0', '0.0 0.0 2.0', &
      'Au', '1', 'Selective dynamics', 'Direct', &
      '0.50 0.50 0.50 T T F', &
      '0.01 0.01 0.01' ])
   ! bad samples (each must fail with ok = .false. and a non-blank errmsg)
   call write_sample('build/specio_xyz_bcount.xyz', [ character(len=64) :: &
      'three', 'comment', 'C 0.0 0.0 0.0' ])
   call write_sample('build/specio_xyz_bzero.xyz', [ character(len=64) :: &
      '0', 'comment' ])
   call write_sample('build/specio_xyz_bshort.xyz', [ character(len=64) :: &
      '2', 'comment', 'C 0.0 0.0 0.0', 'O 1.0 2.0' ])
   call write_sample('build/specio_xyz_bnum.xyz', [ character(len=64) :: &
      '1', 'comment', 'O 1.0 two 3.0' ])
   call write_sample('build/specio_xyz_bsymb.xyz', [ character(len=64) :: &
      '1', 'comment', 'Dummy 0.0 0.0 0.0' ])
   call write_sample('build/specio_xyz_beof.xyz', [ character(len=64) :: &
      '3', 'comment', 'C 0.0 0.0 0.0' ])
   call write_sample('build/specio_poscar_skew.poscar', [ character(len=64) :: &
      'c4.1 skewed-cell Direct sample', '1.0', &
      '1.0 0.0 0.0', '0.5 2.0 0.0', '0.2 0.3 1.5', &
      'H', '1', 'Direct', &
      '0.25 0.50 0.75' ])
   call write_sample('build/specio_poscar_bscale0.poscar', [ character(len=64) :: &
      'c4.1 bad: zero scaling factor', '0.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'H', '1', 'Direct', '0.0 0.0 0.0' ])
   call write_sample('build/specio_poscar_bscale.poscar', [ character(len=64) :: &
      'c4.1 bad: non-numeric scaling factor', '2.oops', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'H', '1', 'Direct', '0.0 0.0 0.0' ])
   call write_sample('build/specio_poscar_bmism.poscar', [ character(len=64) :: &
      'c4.1 bad: symbols-counts field mismatch', '1.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'Si O', '2', 'Direct', '0.0 0.0 0.0', '0.5 0.5 0.5' ])
   call write_sample('build/specio_poscar_bzero.poscar', [ character(len=64) :: &
      'c4.1 bad: zero atoms', '1.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'Si O', '0 0', 'Direct' ])
   call write_sample('build/specio_poscar_bmode.poscar', [ character(len=64) :: &
      'c4.1 bad: unknown mode tag', '1.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'H', '1', 'Fractional', '0.0 0.0 0.0' ])
   call write_sample('build/specio_poscar_bcoord.poscar', [ character(len=64) :: &
      'c4.1 bad: short coordinate row', '1.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'Si', '2', 'Direct', '0.0 0.0 0.0', '0.25 0.75' ])

   ! ---------- (a) read_xyz clean sample ----------
   nfail0 = nfail
   open (newunit=u, file='build/specio_xyz_ok.xyz', status='old', action='read')
   call read_xyz(u, names, xyz, ok, errmsg=emsg)
   close (u)
   res = -1.0d0
   if (allocated(names) .and. allocated(xyz)) then
      if (size(names) == 3 .and. size(xyz,2) == 3 .and. all(shape(xyz) == [3,3])) then
         res = maxval(abs(xyz - wxyz_a))
      end if
   end if
   if (.not. ok .or. len_trim(emsg) /= 0 .or. res /= 0.0d0 .or. &
       .not. allocated(names) .or. size(names) /= 3 .or. &
       names(1) /= 'C' .or. names(2) /= 'O' .or. names(3) /= 'Au') then
      print '(a,es10.2)', 'FAIL (a) read_xyz clean sample: max coord deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', '(a) read_xyz clean sample: PASS (names C/O/Au, coords exact, errmsg blank)'

   ! ---------- (b) read_xyz trailing columns ----------
   nfail0 = nfail
   open (newunit=u, file='build/specio_xyz_extra.xyz', status='old', action='read')
   call read_xyz(u, names, xyz, ok, errmsg=emsg)
   close (u)
   res = -1.0d0
   if (allocated(xyz)) then
      if (size(xyz,2) == 2) res = maxval(abs(xyz - wxyz_b))
   end if
   if (.not. ok .or. len_trim(emsg) /= 0 .or. res /= 0.0d0 .or. &
       .not. allocated(names) .or. size(names) /= 2 .or. names(1) /= 'C' .or. names(2) /= 'O') then
      print '(a,es10.2)', 'FAIL (b) read_xyz trailing columns: max coord deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', '(b) read_xyz trailing columns: PASS (extra fields read past, coords intact)'

   ! ---------- (c) elem_mass spot checks ----------
   nfail0 = nfail
   res = max(abs(elem_mass('H') - 1.008d0), abs(elem_mass('C') - 12.011d0), &
             abs(elem_mass('Au') - 196.97d0), abs(elem_mass('He') - 4.003d0))
   print '(a,es10.2)', '(c) elem_mass spot check: max deviation from IUPAC quotes = ', res
   if (res > 5.0d-3) then
      print '(a)', 'FAIL (c): elem_mass deviates from the IUPAC 2013 values'
      nfail = nfail + 1
   end if
   if (elem_mass('AU') /= elem_mass('Au') .or. elem_mass(' c') /= elem_mass('C')) then
      print '(a)', 'FAIL (c): case/blank normalization broken'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', '(c) elem_mass: PASS (H/C/Au/He IUPAC values, case normalization)'

   ! ---------- (d) read_poscar clean Direct sample ----------
   nfail0 = nfail
   open (newunit=u, file='build/specio_poscar_ok.poscar', status='old', action='read')
   call read_poscar(u, title, cell, names, xyz, ok, errmsg=emsg)
   close (u)
   res = -1.0d0
   if (allocated(xyz)) then
      if (size(xyz,2) == 3) res = maxval(abs(xyz - wxyz_d))
   end if
   if (.not. ok .or. len_trim(emsg) /= 0 .or. trim(title) /= 'c4.1 direct-mode sample' .or. &
       maxval(abs(cell - wcell_d)) > 0.0d0 .or. res /= 0.0d0 .or. &
       .not. allocated(names) .or. size(names) /= 3 .or. &
       names(1) /= 'Si' .or. names(2) /= 'Si' .or. names(3) /= 'O') then
      print '(a,es10.2)', 'FAIL (d) read_poscar Direct sample: max coord deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(d) read_poscar Direct sample: PASS (title/cell/names Si-Si-O/fractional->Cartesian, errmsg blank)'

   ! ---------- (e) read_poscar Cartesian mode ----------
   nfail0 = nfail
   open (newunit=u, file='build/specio_poscar_cart.poscar', status='old', action='read')
   call read_poscar(u, title, cell, names, xyz, ok, errmsg=emsg)
   close (u)
   res = -1.0d0
   if (allocated(xyz)) then
      if (size(xyz,2) == 2) res = maxval(abs(xyz - wxyz_e))
   end if
   if (.not. ok .or. len_trim(emsg) /= 0 .or. maxval(abs(cell - wcell_e)) > 0.0d0 .or. &
       res /= 0.0d0 .or. .not. allocated(names) .or. size(names) /= 2) then
      print '(a,es10.2)', 'FAIL (e) read_poscar Cartesian mode: max coord deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(e) read_poscar Cartesian mode: PASS (universal factor scales cell and coordinates)'

   ! ---------- (f) read_poscar negative scale = target volume ----------
   nfail0 = nfail
   open (newunit=u, file='build/specio_poscar_vol.poscar', status='old', action='read')
   call read_poscar(u, title, cell, names, xyz, ok, errmsg=emsg)
   close (u)
   sexp = (27.0d0/3.0d0)**(1.0d0/3.0d0)          ! expected rescale = (|s|/|det(raw)|)^(1/3)
   detc = cell(1,1)*(cell(2,2)*cell(3,3) - cell(2,3)*cell(3,2)) &
        - cell(1,2)*(cell(2,1)*cell(3,3) - cell(2,3)*cell(3,1)) &
        + cell(1,3)*(cell(2,1)*cell(3,2) - cell(2,2)*cell(3,1))
   res = -1.0d0
   if (allocated(xyz)) then
      if (size(xyz,2) == 1) then
         res = max(abs(xyz(1,1) - 0.10d0*sexp), abs(xyz(2,1) - 0.40d0*sexp), &
                   abs(xyz(3,1) - 0.45d0*sexp))
      end if
   end if
   print '(a,f14.9)', '(f) negative scale: det(cell) = ', detc
   if (.not. ok .or. len_trim(emsg) /= 0 .or. abs(detc - 27.0d0) > 1.0d-9 .or. &
       abs(cell(1,1) - sexp) > 1.0d-12 .or. abs(cell(2,2) - 2.0d0*sexp) > 1.0d-12 .or. &
       abs(cell(3,3) - 1.5d0*sexp) > 1.0d-12 .or. res > 1.0d-12) then
      print '(a,es10.2)', 'FAIL (f) negative scale (target volume): max coord deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(f) negative scale: PASS (|det(cell)| = |s| = 27, isotropic rescale, Direct coords)'

   ! ---------- (g) read_poscar Selective dynamics + velocity block ----------
   nfail0 = nfail
   open (newunit=u, file='build/specio_poscar_sel.poscar', status='old', action='read')
   call read_poscar(u, title, cell, names, xyz, ok, errmsg=emsg)
   close (u)
   res = -1.0d0
   if (allocated(xyz)) then
      if (size(xyz,2) == 1) res = maxval(abs(xyz - wxyz_g))
   end if
   if (.not. ok .or. len_trim(emsg) /= 0 .or. res /= 0.0d0 .or. &
       .not. allocated(names) .or. size(names) /= 1 .or. names(1) /= 'Au' .or. &
       abs(cell(1,1) - 2.0d0) > 0.0d0 .or. abs(cell(2,2) - 2.0d0) > 0.0d0) then
      print '(a,es10.2)', 'FAIL (g) Selective dynamics: max coord deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(g) Selective dynamics: PASS (line skipped, T/F flags read past, velocity block untouched)'

   ! ---------- (h) bad samples -> ok/errmsg channel ----------
   call expect_bad_xyz('build/specio_xyz_bcount.xyz', '(h1) xyz non-numeric count')
   call expect_bad_xyz('build/specio_xyz_bzero.xyz',  '(h2) xyz zero atom count')
   call expect_bad_xyz('build/specio_xyz_bshort.xyz', '(h3) xyz short coordinate row')
   call expect_bad_xyz('build/specio_xyz_bnum.xyz',   '(h4) xyz non-numeric coordinate field')
   call expect_bad_xyz('build/specio_xyz_bsymb.xyz',  '(h5) xyz symbol longer than 4 characters')
   call expect_bad_xyz('build/specio_xyz_beof.xyz',   '(h6) xyz file ends early')
   call expect_bad_poscar('build/specio_poscar_bscale.poscar',  '(h7) POSCAR non-numeric scaling factor')
   call expect_bad_poscar('build/specio_poscar_bscale0.poscar', '(h8) POSCAR zero scaling factor')
   call expect_bad_poscar('build/specio_poscar_bmism.poscar',   '(h9) POSCAR symbols-counts field mismatch')
   call expect_bad_poscar('build/specio_poscar_bzero.poscar',   '(h10) POSCAR zero atoms')
   call expect_bad_poscar('build/specio_poscar_bmode.poscar',   '(h11) POSCAR unknown mode tag')
   call expect_bad_poscar('build/specio_poscar_bcoord.poscar',  '(h12) POSCAR short coordinate row')

   ! ---------- (i) skewed cell, Direct coordinates (fix-round-1 regression) ----------
   ! r = 0.25*a1 + 0.50*a2 + 0.75*a3 with a1=(1,0,0), a2=(0.5,2,0), a3=(0.2,0.3,1.5):
   ! x = 0.25 + 0.25 + 0.15 = 0.65;  y = 1.00 + 0.225 = 1.225;  z = 1.125
   ! (matmul(cell, fr) would instead give the row contractions (0.25, 1.125, 1.325))
   nfail0 = nfail
   open (newunit=u, file='build/specio_poscar_skew.poscar', status='old', action='read')
   call read_poscar(u, title, cell, names, xyz, ok, errmsg=emsg)
   close (u)
   res = -1.0d0
   if (allocated(xyz)) then
      if (size(xyz,2) == 1) res = maxval(abs(xyz - wxyz_i))
   end if
   if (.not. ok .or. len_trim(emsg) /= 0 .or. res /= 0.0d0 .or. &
       .not. allocated(names) .or. size(names) /= 1 .or. names(1) /= 'H' .or. &
       abs(cell(2,1) - 0.5d0) > 0.0d0 .or. abs(cell(3,1) - 0.2d0) > 0.0d0 .or. &
       abs(cell(3,2) - 0.3d0) > 0.0d0) then
      print '(a,es10.2)', 'FAIL (i) skewed-cell Direct: max coord deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(i) skewed-cell Direct: PASS (r = 0.65/1.225/1.125, off-diagonal cell rows verified)'

   print '(a,i0,a)', '=== check_specio: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! write one sample file, one array element per record (samples live under build/)
   subroutine write_sample(fname, lines)
      character(len=*), intent(in) :: fname
      character(len=*), intent(in) :: lines(:)
      integer :: v, i
      open (newunit=v, file=fname, status='replace', action='write')
      do i = 1, size(lines)
         write (v, '(a)') trim(lines(i))
      end do
      close (v)
   end subroutine write_sample

   ! a bad .xyz sample must come back ok = .false. with a non-blank errmsg
   subroutine expect_bad_xyz(fname, label)
      character(len=*), intent(in) :: fname, label
      character(len=4), allocatable :: nm(:)
      real(8), allocatable :: xc(:,:)
      character(len=128) :: em
      logical :: o
      integer :: v
      open (newunit=v, file=fname, status='old', action='read')
      call read_xyz(v, nm, xc, o, errmsg=em)
      close (v)
      if (o .or. len_trim(em) == 0) then
         print '(a)', 'FAIL '//label//': expected ok=.false. + errmsg, got neither'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(em)//')'
      end if
   end subroutine expect_bad_xyz

   ! a bad POSCAR sample must come back ok = .false. with a non-blank errmsg
   subroutine expect_bad_poscar(fname, label)
      character(len=*), intent(in) :: fname, label
      character(len=128) :: tt, em
      character(len=4), allocatable :: nm(:)
      real(8), allocatable :: xc(:,:)
      real(8) :: cc(3,3)
      logical :: o
      integer :: v
      open (newunit=v, file=fname, status='old', action='read')
      call read_poscar(v, tt, cc, nm, xc, o, errmsg=em)
      close (v)
      if (o .or. len_trim(em) == 0) then
         print '(a)', 'FAIL '//label//': expected ok=.false. + errmsg, got neither'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(em)//')'
      end if
   end subroutine expect_bad_poscar

end program check_specio
