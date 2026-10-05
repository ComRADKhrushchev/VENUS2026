!=====================================================================
! sysdef_interface.f90 - system-definition interface: the assembly-phase entry
!   sysdef_load(scan_dir) - discovery, parse, mass lookup, buffer, paradigm
! Design:
!   One file per reactant in the container folder (base name = reactant
!   name; .xyz = no cell, .poscar or POSCAR* = cell-carrying); drag in a
!   folder and the system is fully defined. Generic code owns scan /
!   dispatch / mass lookup / buffer / the cell-presence paradigm rule
!   (any cell -> surface, none -> molecular); the author owns only the
!   file contents. Products sit in module-private buffer behind
!   fresh-copy getters; config_atoms is never used here (list_atoms_load pulls).
!=====================================================================
module sysdef_interface
   use config, only: reactants
   use specio, only: read_poscar, read_xyz, elem_mass
   implicit none
   private
   public :: sysdef_load
   public :: buffer_loaded, buffer_frag_count, buffer_frag_name, buffer_frag_nat, &
             buffer_atom_symbols, buffer_atom_masses, buffer_atom_geometry, &
             buffer_cell_present, buffer_unit_cell

   integer, parameter :: name_len  = 64    ! reactant/file-name capacity [character]
   integer, parameter :: path_len  = 512   ! path/error-message capacity [character]
   integer, parameter :: entry_cap = 256   ! container-folder file capacity [file]

   ! buffer (module-private): the sysdef_load products list_atoms_load transfers from.
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
   ! sysdef_load(scan_dir) - assembly-phase materialization of the system
   !              definition (consumed once by the driver assembly phase)
   !------------------------------------------------------------------
   subroutine sysdef_load(scan_dir)
      character(len=*), intent(in) :: scan_dir  ! system-definition folder (the config-side
                                                ! path field system_dir reaches this
                                                ! dummy at the driver wiring)
      logical :: ok
      character(len=path_len) :: errmsg
      call collect_system(scan_dir, ok, errmsg)
      if (.not. ok) then
         write (0, '(a)') 'sysdef_load: assembly error: '//trim(errmsg)
         write (0, '(a)') 'sysdef_load: fatal (malformed container material is unusable)'
         stop 1
      end if
   end subroutine sysdef_load

   !------------------------------------------------------------------
   ! buffer accessors (public; the list_atoms_load consumption interface;
   !                    all getters return fresh copies)
   !------------------------------------------------------------------
   logical pure function buffer_loaded()
      buffer_loaded = stg_loaded
   end function buffer_loaded

   pure function buffer_frag_count() result(n)
      integer :: n
      n = 0
      if (stg_loaded) n = stg_nfrag
   end function buffer_frag_count

   function buffer_frag_name(i) result(nm)
      integer, intent(in) :: i
      character(len=name_len) :: nm
      nm = ''
      if (stg_loaded .and. i >= 1 .and. i <= stg_nfrag) nm = stg_name(i)
   end function buffer_frag_name

   function buffer_frag_nat(i) result(n)
      integer, intent(in) :: i
      integer :: n
      n = 0
      if (stg_loaded .and. i >= 1 .and. i <= stg_nfrag) n = stg_nat(i)
   end function buffer_frag_nat

   function buffer_atom_symbols() result(s)
      character(len=4), allocatable :: s(:)
      allocate (s(0))
      if (stg_loaded .and. allocated(stg_symb)) s = stg_symb
   end function buffer_atom_symbols

   function buffer_atom_masses() result(m)
      real(8), allocatable :: m(:)
      allocate (m(0))
      if (stg_loaded .and. allocated(stg_mass)) m = stg_mass
   end function buffer_atom_masses

   function buffer_atom_geometry() result(g)
      real(8), allocatable :: g(:,:)
      allocate (g(3,0))
      if (stg_loaded .and. allocated(stg_xyz)) g = stg_xyz
   end function buffer_atom_geometry

   logical pure function buffer_cell_present()
      buffer_cell_present = stg_loaded .and. stg_has_cell
   end function buffer_cell_present

   function buffer_unit_cell() result(c)
      real(8) :: c(3,3)
      c = 0.0d0
      if (stg_loaded .and. stg_has_cell) c = stg_cell
   end function buffer_unit_cell

   !------------------------------------------------------------------
   ! collect_system(scan_dir, ok, errmsg) - the whole assembly flow on the
   !              soft ok/errmsg channel: reset buffer -> folder existence
   !              -> listing -> classify -> name sort -> duplicate check ->
   !              per-file parse -> single-cell check -> concatenate +
   !              elem_mass -> paradigm write into config (sysdef_load
   !              converts a bad return into the named abort)
   !------------------------------------------------------------------
   subroutine collect_system(scan_dir, ok, errmsg)
      character(len=*), intent(in) :: scan_dir
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=name_len) :: entries(entry_cap)
      character(len=path_len) :: path, title, cell_file
      integer :: n_entry, i, j, iat, nat_tot, u, ios
      logical :: ok_sub, any_cell
      character(len=128) :: emsg_sub
      type(material_t) :: mat(entry_cap), swap

      ok = .false.
      errmsg = ''
      call reset_buffer()

      ! 1. folder existence (drag-in contract: the container folder must exist);
      !    probed through the shell - the Intel runtime answers inquire(file=)
      !    with .false. for directories, and gfortran does not accept the F2018
      !    inquire(directory=), so no Fortran-level probe is portable here
      if (.not. folder_exists(scan_dir)) then
         errmsg = 'system-definition folder not found: '//trim(scan_dir)
         return
      end if

      ! 2. listing (files only, no manifest layer)
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
         ! 6b. single-cell rule: the atom list carries ONE unit cell - a second
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

      ! 7. buffer fill: concatenate per-atom rows in name order; masses through the
      !    generic mass table (an unknown symbol aborts with a named error inside elem_mass naming it)
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

      ! 8. paradigm derivation + fragment count: any file carries a cell ->
      !    surface class (=1; the 1-vs-2 split arrives with the surface-member
      !    window); none -> molecular (=0)
      reactants%n_frag = n_entry
      reactants%surface_model = merge(1, 0, any_cell)

      stg_loaded = .true.
      ok = .true.
   end subroutine collect_system

   !------------------------------------------------------------------
   ! private helpers (pure folder/string mechanics - no system knowledge)
   !------------------------------------------------------------------

   ! reset_buffer() - return buffer to the pre-load state (every sysdef_load starts
   ! from a clean slate - repeated loads in one process never leak)
   subroutine reset_buffer()
      stg_loaded = .false.
      stg_nfrag = 0
      stg_has_cell = .false.
      stg_cell = 0.0d0
      if (allocated(stg_name)) deallocate (stg_name)
      if (allocated(stg_nat)) deallocate (stg_nat)
      if (allocated(stg_symb)) deallocate (stg_symb)
      if (allocated(stg_mass)) deallocate (stg_mass)
      if (allocated(stg_xyz)) deallocate (stg_xyz)
   end subroutine reset_buffer

   ! folder_exists(dirp) - directory-existence probe through the shell (the
   ! one portable authority: see the collect_system step-1 note). 'dir /b /ad'
   ! succeeds exactly for an existing directory - a missing path and a regular
   ! file both fail - and 'test -d' answers the same question elsewhere; the
   ! trailing-separator strip and backslash form mirror list_folder
   logical function folder_exists(dirp)
      character(len=*), intent(in) :: dirp
      character(len=path_len) :: dirnorm
      integer :: i, stat
      dirnorm = trim(adjustl(dirp))
      do while (len_trim(dirnorm) > 1)
         i = len_trim(dirnorm)
         if (dirnorm(i:i) == '/' .or. dirnorm(i:i) == '\') then
            dirnorm(i:i) = ' '
         else
            exit
         end if
      end do
      if (windows_host()) then
         dirnorm = to_backslash(dirnorm)
         call execute_command_line('dir /b /ad "'//trim(dirnorm)//'" > nul 2>nul', wait=.true., exitstat=stat)
      else
         call execute_command_line('test -d "'//trim(dirnorm)//'"', wait=.true., exitstat=stat)
      end if
      folder_exists = (stat == 0)
   end function folder_exists

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
         tmpfile = trim(tmpdir)//'/venus2026_sysdef_scan_'//trim(pid_str())//'.txt'
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
   ! Rule: recognized extension first (.xyz / .poscar, case-folded), then the
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
      ! an empty base; named abort rather than buffer an anonymous fragment)
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


   !------------------------------------------------------------------
   ! pid_str() - the process id as a string (the folder-listing scratch
   !                 file must be per-process: a fixed name raced across
   !                 concurrent venus.e runs in one ensemble job)
   !------------------------------------------------------------------
   function pid_str() result(s)
      use iso_c_binding, only: c_int
      character(len=16) :: s
      interface
         function c_getpid() bind(C, name="getpid")
            import :: c_int
            integer(c_int) :: c_getpid
         end function c_getpid
      end interface
      write (s, '(i0)') int(c_getpid())
   end function pid_str

end module sysdef_interface
