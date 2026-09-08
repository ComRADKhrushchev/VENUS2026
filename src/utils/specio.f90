!=====================================================================
! specio.f90 - L5 instrument: system-definition file readers (POSCAR / .xyz) + the element
!              mass table (E14 spec §6: parsers / mass table / derivation rule are generic
!              code with zero system knowledge; the file CONTENTS are blackbox material)
! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §6 (formats: POSCAR
!       (with cell) and .xyz (no cell) initially; a new format = a new parser member,
!       grammar versioning; data flow: scan -> generic parser (utils) -> element symbols
!       against the generic mass table -> fill roster);
!       POSCAR public format definition (comment line / universal scaling factor / three
!       lattice vectors / element-symbols line / per-species counts line / Direct|
!       Cartesian tag / coordinate rows; the scaling factor applies to the cell AND, in
!       Cartesian mode, to the coordinates);
!       .xyz public format definition (atom-count line / comment line / "Symbol x y z"
!       rows; no cell by definition);
!       element mass table: IUPAC standard atomic weights (public physical-constant data;
!       no legacy citation - the old library holds no element-mass table, masses were
!       input data there, venus2026/venus_data.f90 carries no symbol/mass array)
! Ownership note (E12/E13 audit line): this module knows FORMATS and ELEMENTS, not systems
! or methods - nothing here branches on system identity; a new format extends the reader
! member set (grammar versioning), never the existing signatures
!=====================================================================
module specio
   implicit none
   private
   public :: read_poscar, read_xyz, elem_mass

   ! Element mass table (module-private paired arrays; pure physical-constant data, zero
   ! system knowledge): symbol + standard atomic weight per row, H..Og. Values are IUPAC
   ! standard atomic weights [amu], populated as data statements in the C-phase
   ! concretization (the K2 statute is the lookup + loud-stop contract below, not the
   ! numeric table)
   integer, parameter :: n_elem = 118              ! periodic-table row count [row]
   character(len=2) :: tbl_symb(n_elem) = ''       ! element symbols, normalized case [-]
   real(8)          :: tbl_mass(n_elem) = 0.0d0    ! standard atomic weights [amu]
contains
   !------------------------------------------------------------------
   ! read_poscar(u, title, cell, names, xyz, ok [, errmsg]) - parse one POSCAR: comment
   !              line, universal scaling factor, three lattice vectors, element-symbols
   !              line, per-species counts line, Direct/Cartesian tag, coordinate rows
   ! Reference: POSCAR public format definition (VASP wiki); E14 spec §6 (formats row)
   subroutine read_poscar(u, title, cell, names, xyz, ok, errmsg)
      integer, intent(in) :: u    ! unit number (the caller opens/positions and closes the
                                   ! file - this reader is a pure unit consumer, no file
                                   ! names)
      character(len=*), intent(out) :: title  ! line-1 comment (the reactant title; free
                                   ! text by definition, carried out for diagnostics and
                                   ! roster naming - never interpreted)
      real(8), intent(out) :: cell(3,3)        ! lattice vectors by row [Å] (universal
                                   ! scaling applied)
      character(len=4), allocatable, intent(out) :: names(:)  ! per-atom element symbols
                                   ! [-] (species blocks expanded by the counts line;
                                   ! len=4 matches roster%symb)
      real(8), allocatable, intent(out) :: xyz(:,:)  ! per-atom Cartesian coordinates [Å]
                                   ! (3 x natoms, one column per atom)
      logical, intent(out) :: ok                ! clean-parse flag (malformed input is data,
                                   ! not a program fault - soft error channel, no
                                   ! loud-stop in a reader)
      character(len=*), intent(out), optional :: errmsg  ! failure description (names the
                                   ! offending line; meaningful only when ok = .false.)
      ! Pseudocode:
      ! 1. Read the comment line -> title
      ! 2. Read the universal scaling factor s (a negative value means the target cell
      !    VOLUME per the format definition: s <- -(|s|/det(lattice vectors))^(1/3))
      ! 3. Read three lattice-vector rows -> cell = s * vectors (row-wise scaling)
      ! 4. Read the element-symbols line, then the per-species counts line; expand
      !    names(:) block-wise by the counts (concatenation order preserved); allocate
      !    xyz(3, sum(counts)); the counts-less single-line variant some writers emit is a
      !    grammar-versioning matter, not this reader's business
      ! 5. Read the mode tag: first letter D/d = Direct (fractional), C/c = Cartesian; an
      !    optional "Selective dynamics" line may precede it (per-coordinate T/F flags -
      !    read past and NOT interpreted: freezing dof is system knowledge)
      ! 6. Read one coordinate row per atom: Direct -> xyz(:,i) = fractional_i * cell
      !    (the fractional components contract the scaled lattice rows); Cartesian ->
      !    xyz(:,i) = s * cartesian_i (the universal factor scales Cartesian coordinates
      !    too, per the format definition); trailing per-coordinate T/F flags read past
      ! 7. Any short read / non-numeric field / unknown mode tag / zero atoms ->
      !    ok=.false. + errmsg naming the offending line; a clean parse -> ok=.true.
      !    (errmsg left blank)
   end subroutine read_poscar

   !------------------------------------------------------------------
   ! read_xyz(u, names, xyz, ok [, errmsg]) - parse one .xyz frame: atom-count line,
   !              comment line, "Symbol x y z" rows (no cell - by format definition)
   ! Reference: .xyz public format definition; E14 spec §6 (formats row)
   subroutine read_xyz(u, names, xyz, ok, errmsg)
      integer, intent(in) :: u    ! unit number (the caller opens/positions and closes the
                                   ! file - pure unit consumer, no file names)
      character(len=4), allocatable, intent(out) :: names(:)  ! per-atom element symbols
                                   ! [-] (len=4 matches roster%symb)
      real(8), allocatable, intent(out) :: xyz(:,:)  ! per-atom Cartesian coordinates [Å]
                                   ! (3 x natoms, one column per atom; Å per the common
                                   ! convention)
      logical, intent(out) :: ok                ! clean-parse flag (soft error channel,
                                   ! same statute as read_poscar)
      character(len=*), intent(out), optional :: errmsg  ! failure description (names the
                                   ! offending row; meaningful only when ok = .false.)
      ! Pseudocode:
      ! 1. Read the count line: leading integer -> n; allocate names(n), xyz(3,n); a short
      !    read or n <= 0 -> ok=.false. + errmsg
      ! 2. Read the comment line (discarded - free text by definition; the format carries
      !    no cell, hence cell-less systems only, spec §6)
      ! 3. Read n "Symbol x y z" rows -> names(i), xyz(:,i); trailing extended columns read
      !    past, never interpreted
      ! 4. Any short read / non-numeric field / symbol longer than len=4 -> ok=.false. +
      !    errmsg naming the offending row; a clean parse -> ok=.true. (errmsg left blank)
   end subroutine read_xyz

   !------------------------------------------------------------------
   ! elem_mass(symbol) - element-symbol -> atomic mass lookup [amu] (generic
   !              physical-constant table)
   ! Reference: IUPAC standard atomic weights (public physical-constant data; no legacy
   !       citation - the old library holds no element-mass table, masses were input data)
   function elem_mass(symbol) result(m)
      character(len=*), intent(in) :: symbol  ! element symbol (any case, e.g. 'AU'/'Au')
      real(8) :: m                            ! standard atomic mass [amu]
      ! Pseudocode:
      ! 1. Normalize the case (first letter upper, rest lower - the stored table is
      !    normalized; surrounding blanks stripped)
      ! 2. Linear scan of tbl_symb for the match; hit -> m = tbl_mass(row)
      ! 3. Miss -> loud-stop naming the symbol (the pure-function channel has no
      !    soft-error return and an unusable symbol admits no default fallback - unlike
      !    the readers' ok/errmsg channel, this is a hard stop by design)
   end function elem_mass
end module specio
