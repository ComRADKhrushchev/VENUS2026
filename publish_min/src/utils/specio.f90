!=====================================================================
! specio.f90 - system-definition file readers (POSCAR / .xyz) + the element
!              mass table
! Design:
!   Parsers, mass table and derivation rule are generic code with zero system
!   knowledge - the file CONTENTS are container material. This module knows
!   FORMATS and ELEMENTS, not systems: nothing branches on system identity;
!   a new format extends the reader member set, never the existing
!   signatures. Readers are pure unit consumers (the caller opens/positions/
!   closes); malformed input is data - the soft ok/errmsg channel, never a
!   named abort in a reader.
!=====================================================================
module specio
   implicit none
   private
   public :: read_poscar, read_xyz, elem_mass

   ! Element mass table (module-private paired arrays; pure physical-constant data):
   ! symbol + standard atomic weight per row, H..Og. Values are the IUPAC standard
   ! atomic weights, current values incl. CIAAW revisions (not the frozen 2013 list,
   ! e.g. Yb/Hg) [amu]; elements without a stable isotope carry the mass number of
   ! the recognized isotope. DATA-initialized constant table, never mutated
   integer, parameter :: n_elem = 118              ! periodic-table row count [row]
   integer, parameter :: max_sp = 64               ! POSCAR species-line capacity [species]
   character(len=2) :: tbl_symb(n_elem)            ! element symbols, normalized case [-]
   real(8)          :: tbl_mass(n_elem)            ! standard atomic weights [amu]
   data tbl_symb / 'H','He','Li','Be','B','C','N','O','F','Ne', &
                   'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca', &
                   'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                   'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr', &
                   'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                   'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                   'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                   'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg', &
                   'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                   'Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
                   'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
                   'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' /
   data tbl_mass / 1.008d0, 4.002602d0, 6.94d0, 9.0121831d0, 10.81d0, &
                   12.011d0, 14.007d0, 15.999d0, 18.998403163d0, 20.1797d0, &
                   22.98976928d0, 24.305d0, 26.9815385d0, 28.085d0, 30.973761998d0, &
                   32.06d0, 35.45d0, 39.948d0, 39.0983d0, 40.078d0, &
                   44.955908d0, 47.867d0, 50.9415d0, 51.9961d0, 54.938044d0, &
                   55.845d0, 58.933194d0, 58.6934d0, 63.546d0, 65.38d0, &
                   69.723d0, 72.630d0, 74.921595d0, 78.971d0, 79.904d0, &
                   83.798d0, 85.4678d0, 87.62d0, 88.90584d0, 91.224d0, &
                   92.90637d0, 95.95d0, 98.0d0, 101.07d0, 102.90550d0, &
                   106.42d0, 107.8682d0, 112.414d0, 114.818d0, 118.710d0, &
                   121.760d0, 127.60d0, 126.90447d0, 131.293d0, 132.90545196d0, &
                   137.327d0, 138.90547d0, 140.116d0, 140.90766d0, 144.242d0, &
                   145.0d0, 150.36d0, 151.964d0, 157.25d0, 158.92535d0, &
                   162.500d0, 164.93033d0, 167.259d0, 168.93422d0, 173.045d0, &
                   174.9668d0, 178.49d0, 180.94788d0, 183.84d0, 186.207d0, &
                   190.23d0, 192.217d0, 195.084d0, 196.966569d0, 200.592d0, &
                   204.38d0, 207.2d0, 208.98040d0, 209.0d0, 210.0d0, &
                   222.0d0, 223.0d0, 226.0d0, 227.0d0, 232.0377d0, &
                   231.03588d0, 238.02891d0, 237.0d0, 244.0d0, 243.0d0, &
                   247.0d0, 247.0d0, 251.0d0, 252.0d0, 257.0d0, &
                   258.0d0, 259.0d0, 266.0d0, 267.0d0, 268.0d0, &
                   269.0d0, 270.0d0, 269.0d0, 278.0d0, 281.0d0, &
                   282.0d0, 285.0d0, 286.0d0, 289.0d0, 290.0d0, &
                   293.0d0, 294.0d0, 294.0d0 /
contains
   !------------------------------------------------------------------
   ! read_poscar(u, title, cell, names, xyz, ok [, errmsg]) - parse one POSCAR: comment
   !              line, universal scaling factor, three lattice vectors, element-symbols
   !              line, per-species counts line, Direct/Cartesian tag, coordinate rows;
   !              a velocity block following the coordinates is read past, never
   !              interpreted (step 7)
   !------------------------------------------------------------------
   subroutine read_poscar(u, title, cell, names, xyz, ok, errmsg)
      integer, intent(in) :: u    ! unit number (the caller opens/positions and closes the
                                   ! file - this reader is a pure unit consumer, no file
                                   ! names)
      character(len=*), intent(out) :: title  ! line-1 comment (the reactant title; free
                                   ! text by definition, carried out for diagnostics and
                                   ! list_atoms naming - never interpreted)
      real(8), intent(out) :: cell(3,3)        ! lattice vectors by row [Å] (universal
                                   ! scaling applied)
      character(len=4), allocatable, intent(out) :: names(:)  ! per-atom element symbols
                                   ! [-] (species blocks expanded by the counts line;
                                   ! len=4 matches list_atoms%symb)
      real(8), allocatable, intent(out) :: xyz(:,:)  ! per-atom Cartesian coordinates [Å]
                                   ! (3 x natoms, one column per atom)
      logical, intent(out) :: ok                ! clean-parse flag (malformed input is data,
                                   ! not a program fault - soft error channel, no
                                   ! named abort in a reader)
      character(len=*), intent(out), optional :: errmsg  ! failure description (names the
                                   ! offending line; meaningful only when ok = .false.)
      character(len=512) :: line
      character(len=16) :: toks(max_sp)
      character(len=4) :: symb(max_sp)
      character(len=96) :: emsg_row   ! wide enough for the longest missing-row message
      integer :: cnt(max_sp)
      integer :: ios, i, k, nat, nsp, ntok, iat
      real(8) :: s, s_eff, vec(3,3), fr(3), det0
      logical :: direct

      ok = .false.
      if (present(errmsg)) errmsg = ''
      title = ''
      cell = 0.0d0

      ! 1. comment line -> title (free text, never interpreted)
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) then
         call fail(errmsg, 'POSCAR: empty file (comment line missing)')
         return
      end if
      title = trim(line)

      ! 2. universal scaling factor s (a negative value means the target cell VOLUME per
      !    the format definition)
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) then
         call fail(errmsg, 'POSCAR: scaling-factor line missing')
         return
      end if
      read (line, *, iostat=ios) s
      if (ios /= 0) then
         call fail(errmsg, 'POSCAR scaling-factor line: non-numeric: '//trim(line))
         return
      end if
      if (s == 0.0d0) then
         call fail(errmsg, 'POSCAR scaling-factor line: zero scaling factor (all-zero cell)')
         return
      end if

      ! 3. three lattice-vector rows -> vec; cell = effective scale * vectors (row-wise)
      do i = 1, 3
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            write (emsg_row, '(a,i0)') 'POSCAR: lattice-vector row ', i + 2, ' missing'
            call fail(errmsg, trim(emsg_row))
            return
         end if
         read (line, *, iostat=ios) vec(i,:)
         if (ios /= 0) then
            write (emsg_row, '(a,i0,a)') 'POSCAR lattice-vector row ', i + 2, ': short or non-numeric: '
            call fail(errmsg, trim(emsg_row)//' '//trim(line))
            return
         end if
      end do
      if (s < 0.0d0) then
         det0 = det3(vec)
         if (det0 == 0.0d0) then
            call fail(errmsg, 'POSCAR: singular lattice vectors (negative scaling factor = target volume)')
            return
         end if
         s_eff = (abs(s)/abs(det0))**(1.0d0/3.0d0)
      else
         s_eff = s
      end if
      do i = 1, 3
         cell(i,:) = s_eff * vec(i,:)
      end do

      ! 4. element-symbols line, then per-species counts line; expand names(:) block-wise
      !    by the counts (concatenation order preserved)
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) then
         call fail(errmsg, 'POSCAR: element-symbols line missing')
         return
      end if
      call split_tokens(line, toks, ntok)
      if (ntok < 1 .or. ntok > max_sp) then
         call fail(errmsg, 'POSCAR element-symbols line: unreadable or too many species: '//trim(line))
         return
      end if
      nsp = ntok
      do k = 1, nsp
         if (len_trim(toks(k)) > 4) then
            call fail(errmsg, 'POSCAR element-symbols line: symbol longer than 4 characters: '//trim(toks(k)))
            return
         end if
         symb(k) = toks(k)
      end do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) then
         call fail(errmsg, 'POSCAR: counts line missing')
         return
      end if
      call split_tokens(line, toks, ntok)
      if (ntok /= nsp) then
         call fail(errmsg, 'POSCAR counts line: field count does not match the element-symbols line: '//trim(line))
         return
      end if
      do k = 1, nsp
         read (toks(k), *, iostat=ios) cnt(k)
         if (ios /= 0) then
            call fail(errmsg, 'POSCAR counts line: non-integer count: '//trim(toks(k)))
            return
         end if
         if (cnt(k) < 0) then
            call fail(errmsg, 'POSCAR counts line: negative count: '//trim(toks(k)))
            return
         end if
      end do
      nat = sum(cnt(1:nsp))
      if (nat <= 0) then
         call fail(errmsg, 'POSCAR: zero atoms (the counts line sums to zero)')
         return
      end if
      allocate (names(nat), xyz(3,nat))
      iat = 0
      do k = 1, nsp
         do i = 1, cnt(k)
            iat = iat + 1
            names(iat) = symb(k)
         end do
      end do

      ! 5. mode tag: first letter D/d = Direct (fractional), C/c = Cartesian; an optional
      !    "Selective dynamics" line may precede it (per-coordinate T/F flags - read past
      !    and NOT interpreted: freezing dof is system knowledge)
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) then
         call fail(errmsg, 'POSCAR: mode-tag line missing')
         return
      end if
      line = adjustl(line)
      if (line(1:1) == 'S' .or. line(1:1) == 's') then
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            call fail(errmsg, 'POSCAR: mode-tag line missing after the Selective-dynamics line')
            return
         end if
         line = adjustl(line)
      end if
      select case (line(1:1))
      case ('D', 'd')
         direct = .true.
      case ('C', 'c')
         direct = .false.
      case default
         call fail(errmsg, 'POSCAR: unknown mode tag (need Direct or Cartesian): '//trim(line))
         return
      end select

      ! 6. one coordinate row per atom: Direct -> xyz(:,i) = fr_i * cell (the
      !    fractional components LINEARLY COMBINE the row-stored lattice vectors:
      !    r = fr(1)*a1 + fr(2)*a2 + fr(3)*a3 = matmul(fr, cell) - matmul(cell, fr)
      !    would contract row 1 with fr instead, coinciding only for diagonal cells);
      !    Cartesian -> xyz(:,i) = s_eff * cartesian_i (the universal factor scales
      !    Cartesian coordinates too, per the format definition); trailing
      !    per-coordinate T/F flags read past; the coordinate block stops here -
      !    a velocity block, if present, is handled by step 7
      do iat = 1, nat
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            write (emsg_row, '(a,i0,a)') 'POSCAR: coordinate row ', iat, ' missing (file ends inside the coordinate block)'
            call fail(errmsg, trim(emsg_row))
            return
         end if
         read (line, *, iostat=ios) fr
         if (ios /= 0) then
            write (emsg_row, '(a,i0,a)') 'POSCAR coordinate row ', iat, ': short or non-numeric: '
            call fail(errmsg, trim(emsg_row)//' '//trim(line))
            return
         end if
         if (direct) then
            xyz(:,iat) = matmul(fr, cell)
         else
            xyz(:,iat) = s_eff * fr
         end if
      end do

      ! 7. velocity block (first variant): the format allows one velocity row per
      !    atom after the coordinate block - probe one read: end-of-file right
      !    after the coordinates means no block (a clean file ends there); a
      !    present block is READ PAST without interpreting it (velocities are
      !    initial-condition material owned elsewhere - this reader returns
      !    coordinates and cell only); a block cut short (fewer than nat rows) is
      !    a soft error naming the missing velocity row
      do iat = 1, nat
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            if (iat == 1) exit      ! no velocity block - the file ends with the coordinates
            write (emsg_row, '(a,i0,a)') 'POSCAR: velocity row ', iat, &
               ' missing (file ends inside the velocity block)'
            call fail(errmsg, trim(emsg_row))
            return
         end if
      end do

      ! 8. clean parse
      ok = .true.
   end subroutine read_poscar

   !------------------------------------------------------------------
   ! read_xyz(u, names, xyz, ok [, errmsg]) - parse one .xyz frame: atom-count line,
   !              comment line, "Symbol x y z" rows (no cell - by format definition)
   !------------------------------------------------------------------
   subroutine read_xyz(u, names, xyz, ok, errmsg)
      integer, intent(in) :: u    ! unit number (the caller opens/positions and closes the
                                   ! file - pure unit consumer, no file names)
      character(len=4), allocatable, intent(out) :: names(:)  ! per-atom element symbols
                                   ! [-] (len=4 matches list_atoms%symb)
      real(8), allocatable, intent(out) :: xyz(:,:)  ! per-atom Cartesian coordinates [Å]
                                   ! (3 x natoms, one column per atom; Å per the common
                                   ! convention)
      logical, intent(out) :: ok                ! clean-parse flag (soft error channel,
                                   ! same rule as read_poscar)
      character(len=*), intent(out), optional :: errmsg  ! failure description (names the
                                   ! offending row; meaningful only when ok = .false.)
      character(len=512) :: line
      character(len=64) :: emsg_row
      character(len=16) :: sfield
      integer :: ios, n, i
      real(8) :: p(3)

      ok = .false.
      if (present(errmsg)) errmsg = ''

      ! 1. count line: leading integer -> n; short read or n <= 0 -> soft error
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) then
         call fail(errmsg, 'xyz: empty file (atom-count line missing)')
         return
      end if
      read (line, *, iostat=ios) n
      if (ios /= 0) then
         call fail(errmsg, 'xyz atom-count line: not an integer: '//trim(line))
         return
      end if
      if (n <= 0) then
         call fail(errmsg, 'xyz atom-count line: non-positive count')
         return
      end if
      allocate (names(n), xyz(3,n))

      ! 2. comment line (discarded - free text; the format carries no cell, hence
      !    cell-less systems only)
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) then
         call fail(errmsg, 'xyz: comment line missing')
         return
      end if

      ! 3. n "Symbol x y z" rows -> names(i), xyz(:,i); trailing extended columns read
      !    past, never interpreted
      do i = 1, n
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            write (emsg_row, '(a,i0,a)') 'xyz: coordinate row ', i, ' missing (file ends early)'
            call fail(errmsg, trim(emsg_row))
            return
         end if
         read (line, *, iostat=ios) sfield, p(1), p(2), p(3)
         if (ios /= 0) then
            write (emsg_row, '(a,i0,a)') 'xyz coordinate row ', i, ': short or non-numeric: '
            call fail(errmsg, trim(emsg_row)//' '//trim(line))
            return
         end if
         if (len_trim(sfield) > 4) then
            write (emsg_row, '(a,i0,a)') 'xyz coordinate row ', i, ': symbol longer than 4 characters: '
            call fail(errmsg, trim(emsg_row)//' '//trim(sfield))
            return
         end if
         names(i) = sfield
         xyz(:,i) = p
      end do

      ! 4. clean parse
      ok = .true.
   end subroutine read_xyz

   !------------------------------------------------------------------
   ! elem_mass(symbol) - element-symbol -> atomic mass lookup [amu] (generic
   !                     physical-constant table)
   !------------------------------------------------------------------
   function elem_mass(symbol) result(m)
      character(len=*), intent(in) :: symbol  ! element symbol (any case, e.g. 'AU'/'Au')
      real(8) :: m                            ! standard atomic mass [amu]
      character(len=2) :: norm
      integer :: i
      ! 1. normalize the case (first letter upper, rest lower - the stored table is
      !    normalized; surrounding blanks stripped)
      norm = trim(adjustl(symbol))
      norm(1:1) = upcase(norm(1:1))
      norm(2:2) = dncase(norm(2:2))
      ! 2. the deuterium alias (an isotope, not an element - kept out of the
      !    118-row table to preserve the Z indexing; 2.014 amu is a constant
      !    of the record, never re-derived)
      if (norm == 'D') then
         m = 2.014d0
         return
      end if
      ! 3. linear scan of tbl_symb for the match; hit -> m = tbl_mass(row)
      do i = 1, n_elem
         if (tbl_symb(i) == norm) then
            m = tbl_mass(i)
            return
         end if
      end do
      ! 4. miss -> named abort naming the symbol (the pure-function channel has no
      !    soft-error return and an unusable symbol admits no default fallback - unlike
      !    the readers' ok/errmsg channel, this is a hard stop by design)
      write (0, '(a)') 'elem_mass: unknown element symbol: '//trim(symbol)
      stop 1
   end function elem_mass

   !------------------------------------------------------------------
   ! private helpers (pure string/number mechanics - no format knowledge of their own)
   !------------------------------------------------------------------

   ! fail(errmsg, msg) - record the failure description when the caller passed the
   ! optional channel (the ok = .false. flag itself is set by the caller)
   subroutine fail(errmsg, msg)
      character(len=*), intent(out), optional :: errmsg
      character(len=*), intent(in) :: msg
      if (present(errmsg)) errmsg = msg
   end subroutine fail

   ! split_tokens(str, toks, ntok) - blank/tab-separated tokenizer; ntok = size(toks)+1
   ! signals overflow (a token longer than len(toks) is truncated - callers that need
   ! length fidelity check len_trim on their own)
   subroutine split_tokens(str, toks, ntok)
      character(len=*), intent(in) :: str
      character(len=*), intent(out) :: toks(:)
      integer, intent(out) :: ntok
      integer :: i, n, len_tk
      n = len_trim(str)
      ntok = 0
      i = 1
      do while (i <= n .and. ntok <= size(toks))
         do while (i <= n .and. (str(i:i) == ' ' .or. iachar(str(i:i)) == 9))
            i = i + 1
         end do
         if (i > n) exit
         len_tk = 0
         do while (i <= n .and. str(i:i) /= ' ' .and. iachar(str(i:i)) /= 9)
            i = i + 1
            len_tk = len_tk + 1
         end do
         ntok = ntok + 1
         if (ntok <= size(toks)) toks(ntok) = str(i-len_tk:i-1)
      end do
   end subroutine split_tokens

   ! det3(a) - 3x3 determinant (cofactor expansion; used for the negative-scaling-factor
   ! volume rule)
   pure function det3(a) result(d)
      real(8), intent(in) :: a(3,3)
      real(8) :: d
      d = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
        - a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))
   end function det3

   ! upcase(c)/dncase(c) - single-character ASCII case folding for the symbol lookup
   pure function upcase(c) result(uc)
      character(len=1), intent(in) :: c
      character(len=1) :: uc
      uc = c
      if (c >= 'a' .and. c <= 'z') uc = char(iachar(c) - 32)
   end function upcase

   pure function dncase(c) result(dc)
      character(len=1), intent(in) :: c
      character(len=1) :: dc
      dc = c
      if (c >= 'A' .and. c <= 'Z') dc = char(iachar(c) + 32)
   end function dncase

end module specio
