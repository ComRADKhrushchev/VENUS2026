!=====================================================================
! check_sysdef.f90 - C4.2 acceptance criteria program (criteria defined FIRST in
!                    docs/plans/C4.md; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files out
! of the source tree):
!   gfortran -Jbuild globals/config.f90 utils/specio.f90 facade/sysdef_facade.f90 &
!            ../scripts/check_sysdef.f90 -o build/check_sysdef.exe
! All sample folders/files are written under build/ - the source tree stays clean.
! Loud-stop verification approach: the loud-stop criteria (d)-(i) run as CHILD
! PROCESSES - the check re-spawns itself with a mode argument
! (execute_command_line); the child calls the REAL sysdef_load on a prepared bad
! folder and must print a naming message on stderr and exit nonzero. The parent
! asserts the exit status and the stderr content (the print+STOP behavior itself
! is verified, not a softened internal flag); a child that returns without
! stopping exits 0 and fails the criterion.
! Criteria (C4.md C4.2):
!   (a) molecular folder (two .xyz, created OUT of name order): staging in place -
!       staging_loaded() = .false. beforehand; frag count 2; reactant names = the
!       stripped base names in case-insensitive name order (au_cluster before co);
!       per-fragment nat = 2/2; concatenated symbols Au-Au-C-O; masses = elem_mass
!       values exactly; geometry exact; reactants%surface_model = 0 and
!       reactants%n_frag = 2 (D1: n_frag = file count, not input); no cell staged
!   (b) surface folder (poscar_slab.poscar + probe.xyz): reactants%surface_model = 1
!       derived by the D2 cell-presence rule; cell present and exact; both files'
!       atom data concatenated in name order (slab first, probe second)
!   (c) reactant-name rule + case-insensitive extension: a folder with na.XYZ loads
!       clean; reactant name 'na' (base name; extension case folded)
!   (d) empty folder -> loud-stop (child): nonzero exit, stderr names sysdef_load
!       and reports no system-definition files
!   (e) unknown-extension file (notes.txt next to a good .xyz) -> loud-stop (child)
!       naming the offending file
!   (f) malformed .xyz file -> loud-stop (child) echoing the reader errmsg with the
!       file name
!   (g) reactant name given in both formats (co.xyz + co.poscar) -> loud-stop
!       (child) naming the duplicate
!   (h) two cell-carrying files (poscar_a.poscar + bare POSCAR_b) -> loud-stop
!       (child) - the unit cell is a single system-level fact (also exercises the
!       bare-name POSCAR* prefix rule)
!   (i) unknown element symbol (<= 4 characters, not in the table) -> loud-stop
!       (child) via elem_mass naming the symbol
!   (j) re-load resets staging: after the surface load of (b), a molecular re-load
!       leaves exactly the second system (frag count back to 2, surface_model 0,
!       no cell) - no cross-load leakage
!=====================================================================
program check_sysdef
   use sysdef_facade
   use config, only: reactants
   use specio, only: elem_mass
   implicit none
   character(len=4), allocatable :: symb(:)
   real(8), allocatable :: mass(:), geom(:,:)
   real(8) :: cell(3,3), res
   integer :: nfail = 0, nfail0
   character(len=256) :: mode
   ! expected staging geometry of the molecular sample (columns: Au, Au, C, O) [A]
   real(8), parameter :: wq_mol(3,4) = reshape( [ 0.0d0,0.0d0,0.0d0, &
      3.0d0,0.0d0,0.0d0, 0.0d0,0.0d0,0.0d0, 1.20d0,0.0d0,0.0d0 ], [ 3,4 ])
   ! expected staging geometry of the surface sample (columns: Si, Si, H) [A]
   real(8), parameter :: wq_surf(3,3) = reshape( [ 0.0d0,0.0d0,0.0d0, &
      2.5d0,2.5d0,2.5d0, 0.0d0,0.0d0,6.0d0 ], [ 3,3 ])
   ! expected staged unit cell of the surface sample (5.0 * identity) [A]
   real(8), parameter :: wcell_surf(3,3) = reshape( [ 5.0d0,0.0d0,0.0d0, &
      0.0d0,5.0d0,0.0d0, 0.0d0,0.0d0,5.0d0 ], [ 3,3 ])

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('bad-empty');  call sysdef_load('build/sysdef_case_empty')
      case ('bad-ext');    call sysdef_load('build/sysdef_case_ext')
      case ('bad-parse');  call sysdef_load('build/sysdef_case_parse')
      case ('bad-dup');    call sysdef_load('build/sysdef_case_dup')
      case ('bad-2cell');  call sysdef_load('build/sysdef_case_2cell')
      case ('bad-symbol'); call sysdef_load('build/sysdef_case_symb')
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (sysdef_load returned - bug)'
      stop 0
   end if

   ! ---------- sample folders (all created under build/) ----------
   ! molecular: created in NON-name order (co first) to prove the name sort
   call make_dir('build/sysdef_case_mol')
   call write_sample('build/sysdef_case_mol/co.xyz', [ character(len=64) :: &
      '2', 'check_sysdef co sample (created first, sorts second)', &
      'C 0.00 0.00 0.00', 'O 1.20 0.00 0.00' ])
   call write_sample('build/sysdef_case_mol/au_cluster.xyz', [ character(len=64) :: &
      '2', 'check_sysdef au-cluster sample', &
      'Au 0.00 0.00 0.00', 'Au 3.00 0.00 0.00' ])
   ! surface: POSCAR slab + xyz probe (names: poscar_slab < probe)
   call make_dir('build/sysdef_case_surf')
   call write_sample('build/sysdef_case_surf/probe.xyz', [ character(len=64) :: &
      '1', 'check_sysdef probe sample', 'H 0.00 0.00 6.00' ])
   call write_sample('build/sysdef_case_surf/poscar_slab.poscar', [ character(len=64) :: &
      'check_sysdef slab sample', '5.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'Si', '2', 'Direct', '0.00 0.00 0.00', '0.50 0.50 0.50' ])
   ! uppercase-extension folder (extension case rule)
   call make_dir('build/sysdef_case_case')
   call write_sample('build/sysdef_case_case/na.XYZ', [ character(len=64) :: &
      '1', 'check_sysdef uppercase-extension sample', 'Na 0.00 0.00 0.00' ])
   ! bad samples (each folder feeds one child loud-stop case)
   call make_dir('build/sysdef_case_empty')
   call make_dir('build/sysdef_case_ext')
   call write_sample('build/sysdef_case_ext/co.xyz', [ character(len=64) :: &
      '2', 'check_sysdef good file with a bad neighbor', &
      'C 0.00 0.00 0.00', 'O 1.20 0.00 0.00' ])
   call write_sample('build/sysdef_case_ext/notes.txt', [ character(len=64) :: &
      'not a system-definition file' ])
   call make_dir('build/sysdef_case_parse')
   call write_sample('build/sysdef_case_parse/broken.xyz', [ character(len=64) :: &
      '2', 'check_sysdef malformed sample', 'C 0.0 0.0 0.0', 'O 1.0 2.0' ])
   call make_dir('build/sysdef_case_dup')
   call write_sample('build/sysdef_case_dup/co.xyz', [ character(len=64) :: &
      '2', 'check_sysdef duplicate-name sample (xyz side)', &
      'C 0.00 0.00 0.00', 'O 1.20 0.00 0.00' ])
   call write_sample('build/sysdef_case_dup/co.poscar', [ character(len=64) :: &
      'check_sysdef duplicate-name sample (poscar side)', '1.0', &
      '3.0 0.0 0.0', '0.0 3.0 0.0', '0.0 0.0 3.0', &
      'C O', '1 1', 'Direct', '0.00 0.00 0.00', '0.10 0.10 0.10' ])
   call make_dir('build/sysdef_case_2cell')
   call write_sample('build/sysdef_case_2cell/poscar_a.poscar', [ character(len=64) :: &
      'check_sysdef first cell-carrying file', '2.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'Si', '1', 'Direct', '0.00 0.00 0.00' ])
   call write_sample('build/sysdef_case_2cell/POSCAR_b', [ character(len=64) :: &
      'check_sysdef second cell-carrying file (bare POSCAR* name)', '3.0', &
      '1.0 0.0 0.0', '0.0 1.0 0.0', '0.0 0.0 1.0', &
      'Si', '1', 'Direct', '0.00 0.00 0.00' ])
   call make_dir('build/sysdef_case_symb')
   call write_sample('build/sysdef_case_symb/odd.xyz', [ character(len=64) :: &
      '1', 'check_sysdef unknown-symbol sample', 'Zz 0.0 0.0 0.0' ])

   ! ---------- (a) molecular folder ----------
   nfail0 = nfail
   if (staging_loaded()) then
      print '(a)', 'FAIL (a): staging_loaded() = .true. before any load'
      nfail = nfail + 1
   end if
   call sysdef_load('build/sysdef_case_mol')
   symb = staging_atom_symbols()
   mass = staging_atom_masses()
   geom = staging_atom_geometry()
   res = -1.0d0
   if (allocated(geom) .and. all(shape(geom) == [3,4])) res = maxval(abs(geom - wq_mol))
   if (.not. staging_loaded() .or. staging_frag_count() /= 2 .or. &
       trim(staging_frag_name(1)) /= 'au_cluster' .or. trim(staging_frag_name(2)) /= 'co' .or. &
       staging_frag_nat(1) /= 2 .or. staging_frag_nat(2) /= 2 .or. &
       .not. allocated(symb) .or. size(symb) /= 4 .or. &
       symb(1) /= 'Au' .or. symb(2) /= 'Au' .or. symb(3) /= 'C' .or. symb(4) /= 'O' .or. &
       .not. allocated(mass) .or. size(mass) /= 4 .or. &
       mass(1) /= elem_mass('Au') .or. mass(2) /= elem_mass('Au') .or. &
       mass(3) /= elem_mass('C') .or. mass(4) /= elem_mass('O') .or. &
       res /= 0.0d0 .or. reactants%surface_model /= 0 .or. reactants%n_frag /= 2 .or. &
       staging_cell_present()) then
      print '(a,es10.2)', 'FAIL (a) molecular folder: staging mismatch, max geom deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(a) molecular folder: PASS (2 frags name-sorted au_cluster+co, symb/mass/geom exact, surface_model=0, no cell)'

   ! ---------- (b) surface folder ----------
   nfail0 = nfail
   call sysdef_load('build/sysdef_case_surf')
   symb = staging_atom_symbols()
   mass = staging_atom_masses()
   geom = staging_atom_geometry()
   cell = staging_unit_cell()
   res = -1.0d0
   if (allocated(geom) .and. all(shape(geom) == [3,3])) res = maxval(abs(geom - wq_surf))
   if (.not. staging_loaded() .or. staging_frag_count() /= 2 .or. &
       trim(staging_frag_name(1)) /= 'poscar_slab' .or. trim(staging_frag_name(2)) /= 'probe' .or. &
       staging_frag_nat(1) /= 2 .or. staging_frag_nat(2) /= 1 .or. &
       .not. allocated(symb) .or. size(symb) /= 3 .or. &
       symb(1) /= 'Si' .or. symb(2) /= 'Si' .or. symb(3) /= 'H' .or. &
       .not. allocated(mass) .or. size(mass) /= 3 .or. &
       mass(1) /= elem_mass('Si') .or. mass(3) /= elem_mass('H') .or. &
       res /= 0.0d0 .or. reactants%surface_model /= 1 .or. reactants%n_frag /= 2 .or. &
       .not. staging_cell_present() .or. maxval(abs(cell - wcell_surf)) > 0.0d0) then
      print '(a,es10.2)', 'FAIL (b) surface folder: staging mismatch, max geom deviation ', res
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(b) surface folder: PASS (surface_model=1 derived, cell 5*I exact, Si-Si-H concatenated in name order)'

   ! ---------- (c) reactant-name rule + case-insensitive extension ----------
   nfail0 = nfail
   call sysdef_load('build/sysdef_case_case')
   symb = staging_atom_symbols()
   mass = staging_atom_masses()
   if (.not. staging_loaded() .or. staging_frag_count() /= 1 .or. &
       trim(staging_frag_name(1)) /= 'na' .or. staging_frag_nat(1) /= 1 .or. &
       .not. allocated(symb) .or. size(symb) /= 1 .or. symb(1) /= 'Na' .or. &
       .not. allocated(mass) .or. mass(1) /= elem_mass('Na') .or. &
       reactants%surface_model /= 0 .or. staging_cell_present()) then
      print '(a)', 'FAIL (c) name rule / uppercase extension: staging mismatch'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(c) name rule: PASS (na.XYZ -> reactant name na, extension case folded, molecular paradigm)'

   ! ---------- (d)-(i) loud-stop criteria (child processes) ----------
   ! third argument = the naming authority expected on stderr (the unknown-symbol stop
   ! is elem_mass's own hard stop by the specio C4.1 statute, not the sysdef banner)
   call expect_loudstop('bad-empty',  'no system-definition files', 'sysdef_load', &
                        '(d) empty folder loud-stop')
   call expect_loudstop('bad-ext',    'notes.txt', 'sysdef_load', &
                        '(e) unknown-extension loud-stop')
   call expect_loudstop('bad-parse',  'broken.xyz', 'sysdef_load', &
                        '(f) malformed-file loud-stop')
   call expect_loudstop('bad-dup',    'duplicate reactant name', 'sysdef_load', &
                        '(g) duplicate-name loud-stop')
   call expect_loudstop('bad-2cell',  'multiple cell-carrying files', 'sysdef_load', &
                        '(h) two-cell loud-stop')
   call expect_loudstop('bad-symbol', 'Zz', 'elem_mass', &
                        '(i) unknown-symbol loud-stop')

   ! ---------- (j) re-load resets staging ----------
   nfail0 = nfail
   call sysdef_load('build/sysdef_case_mol')
   if (staging_frag_count() /= 2 .or. &
       trim(staging_frag_name(1)) /= 'au_cluster' .or. trim(staging_frag_name(2)) /= 'co' .or. &
       reactants%surface_model /= 0 .or. staging_cell_present()) then
      print '(a)', 'FAIL (j) re-load: staging leaked across loads'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(j) re-load: PASS (staging fully reset - no cross-load leakage)'

   print '(a,i0,a)', '=== check_sysdef: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! make one sample folder (idempotent - a second run tolerates its own leftovers)
   subroutine make_dir(path)
      character(len=*), intent(in) :: path
      if (windows_host()) then
         call execute_command_line('mkdir "'//trim(path)//'"', wait=.true.)
      else
         call execute_command_line('mkdir -p "'//trim(path)//'"', wait=.true.)
      end if
   end subroutine make_dir

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

   ! spawn the child for one bad case and assert the loud-stop: nonzero exit,
   ! stderr names the expected authority (sysdef_load / elem_mass) and carries
   ! the expected key
   subroutine expect_loudstop(cmode, key, who, label)
      character(len=*), intent(in) :: cmode, key, who, label
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      ! the shell redirect needs native separators on Windows (forward slashes make
      ! cmd misparse the redirect target)
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_sysdef.exe '//trim(cmode)//' 2>build'//sep//'sysdef_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/sysdef_err.txt', status='old', action='read', iostat=ios)
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
      else if (index(buf, key) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key)//'"'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop

   ! Windows-host probe (same runtime OS branch the module uses for dir/ls)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_sysdef
