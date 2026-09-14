!=====================================================================
! sysdef_facade.f90 - L2 system-definition facade (E14 spec §6): one assembly-phase entry,
!                    sysdef_load(scan_dir) — folder discovery, parser dispatch, mass
!                    lookup, staging fill, paradigm derivation
! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §6 (granularity: one
!       file per reactant, file base name = reactant name; discovery: assembly-phase scan
!       of the container folder for known extensions, no manifest layer; formats: POSCAR
!       (with cell) / .xyz (no cell) initially, a new format = a new parser member,
!       grammar versioning; derivation rule: cell present -> surface, absent -> molecular
!       - a generic formal rule lawfully in generic code; data flow: scan -> generic
!       parser (utils) -> element symbols against the generic mass table -> staging fill
!       -> derive paradigm; ownership: file contents = container material, parsers /
!       mass table / derivation rule = generic code), §3 (assembly-phase entries:
!       sysdef_load is consumed by the driver assembly phase - no evolution-loop
!       presence), §9 (affected-files rows: config surface_model derived-not-input;
!       config_roster roster_load data source = sysdef facade products);
!       src/utils/specio.f90 (the readers + mass table dispatched into);
!       src/globals/config_roster.f90 (roster_t - the transfer-assembly target, NOT
!       used here - see the S2 ruling below);
!       venus2026/src_VENUS/SELECT.f90:1300-1334 setup_surface_collision / :1337-1404
!       setup_gas_collision (citation anchor: the two routes the derived paradigm selects
!       - the old paradigm fork was code, the new one is data)
!
! -- Data-material contract, not a fifth entry (spec §3/§6) --
!   The evolution loop's four entries (force_eval / elec_prop / mqc_statistic / predicate)
!   run per step; sysdef_load runs ONCE in the assembly phase and produces DATA (staged
!   system products + the derived paradigm code), not callable physics. The geometric
!   definition of the reactants is system knowledge - container material in the E11/E12
!   sense: drag in a folder and the system is fully defined (no manifest to edit, no
!   registration lines).
!
! -- Ownership split (spec §6 ownership row; E12/E13 audit line) --
!   Generic here (zero system knowledge): the extension scan, per-extension dispatch,
!   symbol->mass lookup, the staging mechanics, the cell-presence paradigm rule. The
!   system author owns only the file CONTENTS (geometries, cells, elements). Nothing in
!   this module branches on system identity; the paradigm code it writes is a formal fact
!   of the files (cell present/absent), not a physics judgment.
!
! -- S2 staging-transfer ruling (plan C4.md S2; E14/T11 open point closed) --
!   sysdef_load's parse products are stored in MODULE-PRIVATE staging below; this module
!   NEVER uses config_roster (the inter-facade zero-dependency statute reads "no reverse
!   reach" - the transfer is one-way: config_roster pulls, ruling A). Staging usage for
!   roster_load (C4.3):
!     guard      if (.not. staging_loaded()) -> loud-stop (assembly-order error)
!     fragments  n = staging_frag_count(); per row staging_frag_name(i) /
!                staging_frag_nat(i) (rows in case-insensitive reactant-name order)
!     per-atom   staging_atom_symbols() / staging_atom_masses() / staging_atom_geometry()
!                (files concatenated in the same order; geometry 3 x n_atoms columns [A])
!     unit cell  staging_cell_present() / staging_unit_cell() (at most ONE cell-carrying
!                file per system - the unit cell is a single system-level fact)
!   Every getter returns a fresh copy - staging itself stays intact for further readers.
!   reactants%n_frag (file count, D1) and reactants%surface_model (D2) are written here
!   directly into config - representation-layer data, not roster.
!
! -- Discovery rule (D1; generic formal logic) --
!   Scan the container folder for FILES (subdirectories are not part of the drag-in
!   contract): extension .xyz (case-folded) -> cell-less xyz format; extension .poscar
!   or a base name starting POSCAR* (case-folded, e.g. a bare "POSCAR" or "POSCAR.au")
!   -> cell-carrying POSCAR format; a recognized extension wins over the prefix rule.
!   File base name (extension stripped) = reactant name; n_frag = file count (NOT an
!   input key - E14 "one reactant one file"). An empty folder, an unknown-extension
!   file, an extension-only (empty-reactant-name) file, a duplicate reactant name
!   (fold-case - the comparison matches the sort key), a reader soft failure, or more
!   than one cell-carrying file is an ASSEMBLY ERROR -> loud-stop (the ok/errmsg channel belongs
!   to the specio readers; this facade converts a bad return into a hard stop).
!
! -- Field notes (grammar lands C4; the field statutes are Task 11's) --
!   scan_dir is a character dummy of this entry; the config-side path field (system_dir -
!   the container folder path as config data) reaches this dummy at the driver wiring.
!   The paradigm field reactants%surface_model is an integer code (config.f90: 0 = gas
!   phase, 1/2 = the two surface paradigms): the molecular/surface fork is derivable
!   from file shape and is written here (derived, not input - D2 lands "any file carries
!   a cell -> 1; none -> 0"); the 1-vs-2 surface-paradigm selection is NOT derivable
!   from cell presence alone and arrives with the C5 surface-member window.
!=====================================================================
module sysdef_facade
   use config, only: reactants
   use specio, only: read_poscar, read_xyz, elem_mass
   implicit none
   private
   public :: sysdef_load
   public :: staging_loaded, staging_frag_count, staging_frag_name, staging_frag_nat, &
             staging_atom_symbols, staging_atom_masses, staging_atom_geometry, &
             staging_cell_present, staging_unit_cell

   integer, parameter :: name_len  = 64    ! reactant/file-name capacity [character]
   integer, parameter :: path_len  = 512   ! path/error-message capacity [character]
   integer, parameter :: entry_cap = 256   ! container-folder file capacity [file]

   ! staging (module-private, S2): the sysdef_load products roster_load transfers from.
   ! Filled fresh by every sysdef_load call (reset first - no cross-load leakage);
   ! read out only through the public getters above.
   logical                              :: stg_loaded = .false.
   integer                              :: stg_nfrag  = 0
   character(len=name_len), allocatable :: stg_name(:)   ! reactant names, name-sorted [-]
   integer,                 allocatable :: stg_nat(:)    ! per-fragment atom counts [count]
   character(len=4),        allocatable :: stg_symb(:)   ! per-atom symbols, concatenated [-]
   real(8),                 allocatable :: stg_mass(:)   ! per-atom masses [amu]
   real(8),                 allocatable :: stg_xyz(:,:)  ! per-atom geometry [A] (3 x nat)
   logical                              :: stg_has_cell = .false.
   real(8)                              :: stg_cell(3,3) = 0.0d0  ! unit cell [A]

   ! per-file parse material (private; one row per container-folder file)
   integer, parameter :: fmt_xyz = 1, fmt_poscar = 2
   type :: material_t
      character(len=name_len)     :: file = ''          ! folder entry as listed
      character(len=name_len)     :: name = ''          ! reactant name (extension stripped)
      integer                     :: fmt  = 0           ! fmt_xyz / fmt_poscar
      integer                     :: nat  = 0           ! atom count [count]
      character(len=4), allocatable :: symb(:)          ! per-atom symbols [-]
      real(8),            allocatable :: xyz(:,:)       ! per-atom geometry [A] (3 x nat)
      real(8)                     :: cell(3,3) = 0.0d0  ! lattice vectors [A] (POSCAR only)
   end type
contains
   !------------------------------------------------------------------
   ! sysdef_load(scan_dir) - assembly-phase materialization of the system definition (the
   !              spec §6 data flow end to end; consumed once by the driver assembly phase)
   ! Reference: spec §6 (data flow row); §3 (assembly-phase entries); §9 (config_roster
   !       data source = sysdef facade products); implemented steps mirror collect_system
   subroutine sysdef_load(scan_dir)
      character(len=*), intent(in) :: scan_dir  ! system-definition folder (the container
                                                ! folder; the config-side path field
                                                ! system_dir reaches this dummy at the
                                                ! driver wiring)
      logical :: ok
      character(len=path_len) :: errmsg
      call collect_system(scan_dir, ok, errmsg)
      if (.not. ok) then
         write (0, '(a)') 'sysdef_load: assembly error: '//trim(errmsg)
         write (0, '(a)') 'sysdef_load: loud-stop (malformed container material is unusable)'
         stop 1
      end if
   end subroutine sysdef_load

   !------------------------------------------------------------------
   ! staging accessors (public; the C4.3 roster_load consumption interface - see the
   !                    S2 header block; all getters return fresh copies)
   !------------------------------------------------------------------
   logical pure function staging_loaded()
      staging_loaded = stg_loaded
   end function staging_loaded

   pure function staging_frag_count() result(n)
      integer :: n
      n = 0
      if (stg_loaded) n = stg_nfrag
   end function staging_frag_count

   function staging_frag_name(i) result(nm)
      integer, intent(in) :: i
      character(len=name_len) :: nm
      nm = ''
      if (stg_loaded .and. i >= 1 .and. i <= stg_nfrag) nm = stg_name(i)
   end function staging_frag_name

   function staging_frag_nat(i) result(n)
      integer, intent(in) :: i
      integer :: n
      n = 0
      if (stg_loaded .and. i >= 1 .and. i <= stg_nfrag) n = stg_nat(i)
   end function staging_frag_nat

   function staging_atom_symbols() result(s)
      character(len=4), allocatable :: s(:)
      allocate (s(0))
      if (stg_loaded .and. allocated(stg_symb)) s = stg_symb
   end function staging_atom_symbols

   function staging_atom_masses() result(m)
      real(8), allocatable :: m(:)
      allocate (m(0))
      if (stg_loaded .and. allocated(stg_mass)) m = stg_mass
   end function staging_atom_masses

   function staging_atom_geometry() result(g)
      real(8), allocatable :: g(:,:)
      allocate (g(3,0))
      if (stg_loaded .and. allocated(stg_xyz)) g = stg_xyz
   end function staging_atom_geometry

   logical pure function staging_cell_present()
      staging_cell_present = stg_loaded .and. stg_has_cell
   end function staging_cell_present

   function staging_unit_cell() result(c)
      real(8) :: c(3,3)
      c = 0.0d0
      if (stg_loaded .and. stg_has_cell) c = stg_cell
   end function staging_unit_cell

   !------------------------------------------------------------------
   ! collect_system(scan_dir, ok, errmsg) - the whole assembly flow on the soft
   !              ok/errmsg channel (detection layer; sysdef_load converts a bad
   !              return into the loud-stop): reset staging -> folder existence ->
   !              listing -> classify -> name sort -> duplicate check -> per-file
   !              parse (specio readers, unit-number contract) -> single-cell check
   !              -> concatenate + elem_mass -> paradigm write into config
   ! Reference: spec §6 data flow; D1/D2 (plan C4.md)
   subroutine collect_system(scan_dir, ok, errmsg)
      character(len=*), intent(in) :: scan_dir
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=name_len) :: entries(entry_cap)
      character(len=path_len) :: path, title, cell_file
      integer :: n_entry, i, j, iat, nat_tot, u, ios
      logical :: ex, ok_sub, any_cell
      character(len=128) :: emsg_sub
      type(material_t) :: mat(entry_cap), swap

      ok = .false.
      errmsg = ''
      call reset_staging()

      ! 1. folder existence (drag-in contract: the container folder must exist)
      inquire (file=trim(scan_dir), exist=ex)
      if (.not. ex) then
         errmsg = 'system-definition folder not found: '//trim(scan_dir)
         return
      end if

      ! 2. listing (D1 scan; files only, no manifest layer)
      call list_folder(scan_dir, entries, n_entry, ok_sub, emsg_sub)
      if (.not. ok_sub) then
         errmsg = emsg_sub
         return
      end if
      if (n_entry == 0) then
         errmsg = 'no system-definition files in the container folder: '//trim(scan_dir)
         return
      end if

      ! 3. classify every entry (extension / POSCAR* prefix -> format + reactant name)
      do i = 1, n_entry
         call classify_entry(entries(i), mat(i)%fmt, mat(i)%name, emsg_sub)
         if (len_trim(emsg_sub) > 0) then
            errmsg = emsg_sub
            return
         end if
         mat(i)%file = entries(i)
      end do

      ! 4. name order (case-insensitive insertion sort - a deterministic fragment
      !    order independent of the file system's listing order)
      do i = 2, n_entry
         swap = mat(i)
         j = i - 1
         do while (j >= 1 .and. fold_case(trim(mat(j)%name)) > fold_case(trim(swap%name)))
            mat(j+1) = mat(j)
            j = j - 1
         end do
         mat(j+1) = swap
      end do

      ! 5. duplicate reactant names (one reactant one file - the same name in two
      !    formats would define the reactant twice); the comparison is fold-case LIKE
      !    THE SORT KEY - Co.xyz + co.xyz is one reactant spelled twice, not two
      do i = 2, n_entry
         if (fold_case(trim(mat(i)%name)) == fold_case(trim(mat(i-1)%name))) then
            errmsg = 'duplicate reactant name: '//trim(mat(i)%name)// &
                     ' ('//trim(mat(i-1)%file)//' and '//trim(mat(i)%file)//')'
            return
         end if
      end do

      ! 6. per-file parse through the specio readers (unit-number contract: open here,
      !    hand the unit over, close; a reader soft failure is an assembly error)
      any_cell = .false.
      do i = 1, n_entry
         path = trim(scan_dir)//'/'//trim(mat(i)%file)
         open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
         if (ios /= 0) then
            errmsg = 'cannot open system-definition file: '//trim(path)
            return
         end if
         if (mat(i)%fmt == fmt_xyz) then
            call read_xyz(u, mat(i)%symb, mat(i)%xyz, ok_sub, errmsg=emsg_sub)
         else
            call read_poscar(u, title, mat(i)%cell, mat(i)%symb, mat(i)%xyz, ok_sub, errmsg=emsg_sub)
         end if
         close (u)
         if (.not. ok_sub) then
            errmsg = trim(mat(i)%file)//': '//trim(emsg_sub)
            return
         end if
         mat(i)%nat = size(mat(i)%symb)
         ! 6b. single-cell statute: the roster carries ONE unit cell - a second
         !     cell-carrying file is a half-periodic definition no paradigm covers
         if (mat(i)%fmt == fmt_poscar) then
            if (any_cell) then
               errmsg = 'multiple cell-carrying files ('//trim(cell_file)//' and '// &
                        trim(mat(i)%file)//'): the unit cell is a single system-level fact'
               return
            end if
            stg_cell = mat(i)%cell
            cell_file = mat(i)%file
            any_cell = .true.
         end if
      end do

      ! 7. staging fill: concatenate per-atom rows in name order; masses through the
      !    generic mass table (an unknown symbol loud-stops inside elem_mass naming it)
      nat_tot = 0
      do i = 1, n_entry
         nat_tot = nat_tot + mat(i)%nat
      end do
      allocate (stg_name(n_entry), stg_nat(n_entry))
      allocate (stg_symb(nat_tot), stg_mass(nat_tot), stg_xyz(3,nat_tot))
      stg_nfrag = n_entry
      iat = 0
      do i = 1, n_entry
         stg_name(i) = mat(i)%name
         stg_nat(i) = mat(i)%nat
         do j = 1, mat(i)%nat
            iat = iat + 1
            stg_symb(iat) = mat(i)%symb(j)
            stg_mass(iat) = elem_mass(mat(i)%symb(j))
            stg_xyz(:,iat) = mat(i)%xyz(:,j)
         end do
      end do
      stg_has_cell = any_cell

      ! 8. paradigm derivation (D2) + fragment count (D1): any file carries a cell ->
      !    surface class (=1; the 1-vs-2 split is the C5 window); none -> molecular (=0)
      reactants%n_frag = n_entry
      reactants%surface_model = merge(1, 0, any_cell)

      stg_loaded = .true.
      ok = .true.
   end subroutine collect_system

   !------------------------------------------------------------------
   ! private helpers (pure folder/string mechanics - no system knowledge)
   !------------------------------------------------------------------

   ! reset_staging() - return staging to the pre-load state (every sysdef_load starts
   ! from a clean slate - repeated loads in one process never leak)
   subroutine reset_staging()
      stg_loaded = .false.
      stg_nfrag = 0
      stg_has_cell = .false.
      stg_cell = 0.0d0
      if (allocated(stg_name)) deallocate (stg_name)
      if (allocated(stg_nat)) deallocate (stg_nat)
      if (allocated(stg_symb)) deallocate (stg_symb)
      if (allocated(stg_mass)) deallocate (stg_mass)
      if (allocated(stg_xyz)) deallocate (stg_xyz)
   end subroutine reset_staging

   ! list_folder(dirp, entries, n_entry, ok, errmsg) - directory listing through the
   ! shell (no portable Fortran intrinsic exists): 'dir /b /a-d' on Windows, 'ls -1'
   ! elsewhere, output to a scratch file in the system temp directory. Files only;
   ! the listing failure of an EXISTING empty folder is indistinguishable from a
   ! listing of zero files - both surface as n_entry = 0 and the caller reports the
   ! empty-container error (folder existence was gated beforehand)
   subroutine list_folder(dirp, entries, n_entry, ok, errmsg)
      character(len=*), intent(in) :: dirp
      character(len=*), intent(out) :: entries(:)
      integer, intent(out) :: n_entry
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=path_len) :: cmd, tmpfile, dirnorm, line, tmpdir
      integer :: u, ios, i

      ok = .false.
      errmsg = ''
      n_entry = 0
      entries = ''

      ! normalize: trim, strip trailing separators; Windows shell wants backslashes
      dirnorm = trim(adjustl(dirp))
      do while (len_trim(dirnorm) > 1)
         i = len_trim(dirnorm)
         if (dirnorm(i:i) == '/' .or. dirnorm(i:i) == '\') then
            dirnorm(i:i) = ' '
         else
            exit
         end if
      end do
      if (windows_host()) dirnorm = to_backslash(dirnorm)

      ! scratch listing file in the system temp directory (never in the source tree
      ! or the scanned folder)
      tmpdir = ''
      call get_environment_variable('TEMP', tmpdir)
      if (len_trim(tmpdir) == 0) call get_environment_variable('TMPDIR', tmpdir)
      if (len_trim(tmpdir) == 0) tmpdir = '.'
      if (windows_host()) then
         tmpfile = trim(to_backslash(tmpdir))//'\venus2026_sysdef_scan.txt'
      else
         tmpfile = trim(tmpdir)//'/venus2026_sysdef_scan.txt'
      end if

      if (windows_host()) then
         cmd = 'dir /b /a-d "'//trim(dirnorm)//'" 2>nul > "'//trim(tmpfile)//'"'
      else
         cmd = 'ls -1 "'//trim(dirnorm)//'" > "'//trim(tmpfile)//'"'
      end if
      call execute_command_line(trim(cmd), wait=.true.)

      open (newunit=u, file=trim(tmpfile), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         errmsg = 'cannot list the system-definition folder: '//trim(dirp)
         return
      end if
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         n_entry = n_entry + 1
         if (n_entry > size(entries)) then
            errmsg = 'system-definition folder holds more files than the scan capacity: '//trim(dirp)
            close (u)
            return
         end if
         entries(n_entry) = trim(line)
      end do
      close (u)

      ! tidy the scratch listing (best effort)
      if (windows_host()) then
         call execute_command_line('del "'//trim(tmpfile)//'"', wait=.true.)
      else
         call execute_command_line('rm -f "'//trim(tmpfile)//'"', wait=.true.)
      end if
      ok = .true.
   end subroutine list_folder

   ! classify_entry(entry, fmt, rname, err) - one folder entry -> format code +
   !              reactant name; err non-blank = unclassifiable (assembly error)
   ! Rule (D1): recognized extension first (.xyz / .poscar, case-folded), then the
   ! bare-name POSCAR* prefix (POSCAR, POSCAR.au, poscar_au111, ...); base name with
   ! the extension stripped = reactant name
   subroutine classify_entry(entry, fmt, rname, err)
      character(len=*), intent(in) :: entry
      integer, intent(out) :: fmt
      character(len=*), intent(out) :: rname, err
      character(len=name_len) :: base, ext, head
      character(len=len(entry)) :: fname
      integer :: idot
      err = ''
      rname = ''
      fmt = 0
      fname = trim(entry)
      idot = index(trim(fname), '.', back=.true.)
      if (idot > 0) then
         base = fname(:idot-1)
         ext = fold_case(fname(idot+1:))
      else
         base = fname
         ext = ''
      end if
      select case (trim(ext))
      case ('xyz')
         fmt = fmt_xyz
      case ('poscar')
         fmt = fmt_poscar
      case default
         if (len_trim(base) >= 6) then
            head = fold_case(base)
            if (head(1:6) == 'poscar') then
               fmt = fmt_poscar
            end if
         end if
         if (fmt == 0) then
            if (len_trim(base) == 0) then
               err = 'empty reactant name (extension-only file name): '//trim(entry)
            else
               err = 'unknown system-definition file type (need .xyz / .poscar / POSCAR*): '//trim(entry)
            end if
            return
         end if
      end select
      ! the extension-only name has no reactant identity on ANY path - recognized
      ! extensions included (a hidden file named '.xyz' reaches here with fmt set but
      ! an empty base; loud-stop rather than staging an anonymous fragment)
      if (len_trim(base) == 0) then
         err = 'empty reactant name (extension-only file name): '//trim(entry)
         return
      end if
      rname = base
   end subroutine classify_entry

   ! fold_case(s) - ASCII lowercase of a string (name-order key, extension folding)
   pure function fold_case(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i
      t = s
      do i = 1, len(t)
         if (t(i:i) >= 'A' .and. t(i:i) <= 'Z') t(i:i) = char(iachar(t(i:i)) + 32)
      end do
   end function fold_case

   ! to_backslash(p) - forward slashes to Windows backslashes (shell path form)
   pure function to_backslash(p) result(q)
      character(len=*), intent(in) :: p
      character(len=len(p)) :: q
      integer :: i
      q = p
      do i = 1, len_trim(q)
         if (q(i:i) == '/') q(i:i) = '\'
      end do
   end function to_backslash

   ! windows_host() - runtime OS probe (shell-command fork: dir vs ls)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end module sysdef_facade
