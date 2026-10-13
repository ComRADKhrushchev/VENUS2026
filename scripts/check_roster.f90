!=====================================================================
! check_roster.f90 - list_atoms criteria program (criteria defined FIRST in
!                    this header; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files out
! of the source tree; control/input enter because the window keys travel the real
! input route - read_input stages them, list_atoms_load pulls them):
!   gfortran -Jbuild globals/control.f90 globals/config.f90 utils/specio.f90 &
!            interface/sysdef_interface.f90 input.f90 globals/config_atoms.f90 &
!            ../scripts/check_roster.f90 -o build/check_roster.exe
! All sample folders/files are written under build/ - the source tree stays clean.
! Named-abort verification approach (same scheme as check_sysdef): the named abort
! criteria (d)/(e)/(f) run as CHILD PROCESSES - the check re-spawns itself with a
! mode argument (execute_command_line); the child runs the REAL
! sysdef_load -> list_atoms_load chain on a prepared bad case and must print a naming
! message on stderr and exit nonzero. The parent asserts the exit status and the
! stderr content; a child that returns without stopping exits 0 and fails the
! criterion.
! Negative-sample note (D4 nat/list/mass rules + the fill-time bounds guard;
! implementer's choice, stated):
!   the buffer accessors are READ-ONLY and consistent by construction (nat runs
!   sum to the atom count; masses come from elem_mass, all positive) - no folder
!   and no public setter can produce a violating list_atoms, so those three D4 rules
!   AND list_atoms_load's pre-access fill guard (offset+nat <= n_atoms, review fix
!   round 1) are NOT triggerable through the public buffer surface; forcing a
!   backdoor would enlarge that surface, so they are tested as POSITIVE invariants
!   on both assembled end-to-end atom lists (c) instead. The named abort branches
!   remain in list_atoms_load as fail-safe defensive guards for later fill routes
!   (the phase-C surface rows carry their own content); since the fix-round-1
!   reorder the guard fires BEFORE any geometry access, so an oversized buffered
!   nat can no longer surface as an out-of-bounds read ahead of the named stop.
! Criteria:
!   (a) end-to-end molecular chain (acceptance 2): folder (au_cluster.xyz + co.xyz,
!       created OUT of name order) -> read_input (the molecular input stages the
!       role keys PROJECTILE = co / TARGET = au_cluster) -> sysdef_load ->
!       list_atoms_load -> every list_atoms
!       field exact: symb Au-Au-C-O; mass = elem_mass values exactly; frag(:) in
!       the E15 arrayed shape (2 homogeneous rows, nat 2/2, list = contiguous runs
!       1-2 / 3-4, mass = row sums, qz = packed geometry, name = the file stem);
!       q0 = the 3N packed
!       vector; no cell -> skew/box_lx/box_ly 0; size(frag)
!       = reactants%n_frag; surface_model 0; the role rows resolve through the
!       names (i_proj = 2 (co), i_targ = 1 (au_cluster) - the name sort puts
!       au_cluster first)
!   (b) end-to-end surface chain: skewed-cell POSCAR (a1=(4,0,0), a2=(1,3,0),
!       Cartesian mode - see the sample note below) + probe.xyz; the input route
!       supplies the class and cell (read_input stages SCATTERING = SURFACE and
!       CELL = the same skewed vectors - input-owned, validated not derived)
!       and the surface role key (PROJECTILE = probe) ->
!       box_lx = 4, box_ly = sqrt(10),
!       skew = acos(1/sqrt(10)) exactly (law of cosines on the CELL rows);
!       slab rows before probe (name order: a_slab < probe)
!   (c) D4 invariants positive-path: on BOTH assembled atom lists every fragment
!       nat <= atom count, every list index within 1..n_atoms, all masses positive
!       (the molecular re-load here re-reads the MOLECULAR input first - the input
!       buffer persists module-side and must not carry the surface window rows
!       into a molecular load)
!   (d) guard named abort (child): list_atoms_load with no prior sysdef_load -> nonzero
!       exit, stderr names list_atoms_load and the assembly-order error
!   (e) missing window key named abort (child): surface class staged, no
!       CELL buffered -> nonzero exit, stderr names list_atoms_load and the missing
!       required key CELL (a system parameter with no generic default)
!   (f) degenerate-cell named abort (child): a CELL with a zero in-plane
!       vector -> nonzero exit, stderr names list_atoms_load and the degenerate cell
!       (old SURF.f:74-79 rule)
!   (g) list_atoms_clear + re-load: clear deallocates every member and resets every
!       scalar to its default (skew/box_lx/box_ly 0, role rows
!       i_proj/i_targ 0); a
!       molecular re-load (its own input re-read first - same buffer-persistence
!       rule as (c)) after the surface load leaves exactly the molecular
!       list_atoms - no cross-load leakage (the box scalars included: the surface
!       load left box_lx 4, the molecular re-load must return it to 0)
!=====================================================================
program check_roster
   use config_atoms, only: list_atoms, list_atoms_load, list_atoms_clear
   use sysdef_interface, only: sysdef_load
   use config, only: reactants
   use specio, only: elem_mass
   use input, only: read_input
   implicit none
   integer :: nfail = 0, nfail0
   character(len=256) :: mode
   real(8) :: dev
   ! expected molecular q0 (Au, Au, C, O packed; [A]) and surface q0 (Si, Si, H)
   real(8), parameter :: wq0_mol(12) = [ 0.0d0,0.0d0,0.0d0, 3.0d0,0.0d0,0.0d0, &
      0.0d0,0.0d0,0.0d0, 1.20d0,0.0d0,0.0d0 ]
   real(8), parameter :: wq0_surf(9) = [ 2.0d0,0.0d0,0.0d0, 0.5d0,1.5d0,0.0d0, &
      0.0d0,0.0d0,6.0d0 ]
   ! expected surface cell fields: a1=(4,0,0), a2=(1,3,0) -> |a1|=4, |a2|=sqrt(10),
   ! cos(skew) = a1.a2/(|a1||a2|) = 4/(4 sqrt 10) = 1/sqrt(10)
   real(8), parameter :: w_lx = 4.0d0, w_ly = sqrt(10.0d0)
   real(8), parameter :: w_skew = acos(1.0d0/sqrt(10.0d0))

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('c-guard');  call list_atoms_load()
      case ('c-cell');   call read_input('build/roster_case_surf_nocell_input.txt')
                        call sysdef_load('build/roster_case_surf'); call list_atoms_load()
      case ('c-degen');  call read_input('build/roster_case_degen_input.txt')
                        call sysdef_load('build/roster_case_degen'); call list_atoms_load()
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (list_atoms_load returned - bug)'
      stop 0
   end if

   ! ---------- sample folders (all created under build/) ----------
   ! molecular: created in NON-name order (co first) to prove the chain honors the
   ! sysdef name sort (au_cluster before co)
   call make_dir('build/roster_case_mol')
   call write_sample('build/roster_case_mol/co.xyz', [ character(len=64) :: &
      '2', 'check_roster co sample (created first, sorts second)', &
      'C 0.00 0.00 0.00', 'O 1.20 0.00 0.00' ])
   call write_sample('build/roster_case_mol/au_cluster.xyz', [ character(len=64) :: &
      '2', 'check_roster au-cluster sample', &
      'Au 0.00 0.00 0.00', 'Au 3.00 0.00 0.00' ])
   ! surface: skewed-cell slab + xyz probe (names: a_slab < probe). Cartesian mode:
   ! a latent specio defect transposes Direct coordinates on non-symmetric cells
   ! (matmul(cell,fr) contracts the columns; the vectors are stored by row) - reported
   ! to the controller, not fixed here (this check does not touch other modules); Cartesian
   ! mode sidesteps it while still buffer the skewed cell for the law-of-cosines test
   call make_dir('build/roster_case_surf')
   call write_sample('build/roster_case_surf/probe.xyz', [ character(len=64) :: &
      '1', 'check_roster probe sample', 'H 0.00 0.00 6.00' ])
   call write_sample('build/roster_case_surf/a_slab.poscar', [ character(len=64) :: &
      'check_roster skewed slab sample', '1.0', &
      '4.0 0.0 0.0', '1.0 3.0 0.0', '0.0 0.0 5.0', &
      'Si', '2', 'Cartesian', '2.00 0.00 0.00', '0.50 1.50 0.00' ])
   ! degenerate cell: zero-length first in-plane lattice vector
   call make_dir('build/roster_case_degen')
   call write_sample('build/roster_case_degen/flat.poscar', [ character(len=64) :: &
      'check_roster degenerate-cell sample', '1.0', &
      '0.0 0.0 0.0', '0.0 2.0 0.0', '0.0 0.0 2.0', &
      'Si', '1', 'Direct', '0.00 0.00 0.00' ])
   ! (b) surface-window input: the six required keys + the class/cell keys (the
   ! box the surface chain must land in the atom list - input-owned, validated
   ! not derived) + the surface role key
   call write_sample('build/roster_case_surf_input.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/roster_case_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'SCATTERING = SURFACE', &
      'CELL = 4.0 0.0 0.0, 1.0 3.0 0.0, 0.0 0.0 5.0', &
      'PROJECTILE = probe' ])
   ! (e) surface class with the CELL key removed (the missing-window child input)
   call write_sample('build/roster_case_surf_nocell_input.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/roster_case_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'SCATTERING = SURFACE' ])
   ! (f) degenerate CELL (zero-length first in-plane vector; the degenerate child input)
   call write_sample('build/roster_case_degen_input.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/roster_case_degen', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'SCATTERING = SURFACE', &
      'CELL = 0.0 0.0 0.0, 0.0 2.0 0.0, 0.0 0.0 2.0' ])
   ! the molecular-window input: the six required keys + the role keys (the
   ! molecular paradigm admits no window keys; a two-fragment list_atoms resolves
   ! its roles through the names)
   call write_sample('build/roster_case_mol_input.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/roster_case_mol', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'PROJECTILE = co', 'TARGET = au_cluster' ])

   ! ---------- (a) end-to-end molecular chain ----------
   nfail0 = nfail
   call read_input('build/roster_case_mol_input.txt')  ! the input route stages the role keys
   call sysdef_load('build/roster_case_mol')
   call list_atoms_load()
   dev = -1.0d0
   if (allocated(list_atoms%q0)) dev = maxval(abs(list_atoms%q0 - wq0_mol))
   if (.not. molecular_list_atoms_ok() .or. dev /= 0.0d0) then
      print '(a,es10.2)', 'FAIL (a) molecular chain: list_atoms mismatch, max q0 deviation ', dev
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(a) molecular chain: PASS (symb/mass exact, 2 frag rows contiguous 1-2/3-4, q0 packed, no cell)'

   ! ---------- (b) end-to-end surface chain ----------
   nfail0 = nfail
   call read_input('build/roster_case_surf_input.txt')   ! the real input route stages the window keys
   call sysdef_load('build/roster_case_surf')
   call list_atoms_load()
   dev = -1.0d0
   if (allocated(list_atoms%q0)) dev = maxval(abs(list_atoms%q0 - wq0_surf))
   if (.not. surface_list_atoms_ok() .or. dev /= 0.0d0) then
      print '(a,es10.2)', 'FAIL (b) surface chain: list_atoms mismatch, max q0 deviation ', dev
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(b) surface chain: PASS (box_lx=4, box_ly=sqrt(10), skew=acos(1/sqrt(10)) exact from the CELL rows)'

   ! ---------- (c) D4 invariants positive-path (both paradigms) ----------
   nfail0 = nfail
   if (.not. d4_invariants_hold()) then
      print '(a)', 'FAIL (c): D4 invariants violated on the assembled surface list_atoms'
      nfail = nfail + 1
   end if
   call read_input('build/roster_case_mol_input.txt')   ! buffer re-read: the surface window rows must not survive
   call sysdef_load('build/roster_case_mol')
   call list_atoms_load()
   if (.not. d4_invariants_hold()) then
      print '(a)', 'FAIL (c): D4 invariants violated on the assembled molecular list_atoms'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(c) D4 invariants: PASS (nat<=n_atoms, list in range, masses positive - both atom lists)'

   ! ---------- (d)-(f) named-abort criteria (child processes) ----------
   call expect_loudstop('c-guard', 'sysdef_load must run first', 'list_atoms_load', &
                        '(d) no-buffer guard loud-stop')
   call expect_loudstop('c-cell',  'CELL', 'list_atoms_load', &
                        '(e) missing-CELL loud-stop (the required window key)')
   call expect_loudstop('c-degen', 'degenerate unit cell', 'list_atoms_load', &
                        '(f) degenerate-cell loud-stop')

   ! ---------- (g) list_atoms_clear + re-load ----------
   nfail0 = nfail
   call list_atoms_clear()
   if (allocated(list_atoms%mass) .or. allocated(list_atoms%symb) .or. &
       allocated(list_atoms%frag) .or. allocated(list_atoms%q0) .or. &
       list_atoms%skew /= 0.0d0 .or. list_atoms%box_lx /= 0.0d0 .or. list_atoms%box_ly /= 0.0d0 .or. &
       list_atoms%i_proj /= 0 .or. list_atoms%i_targ /= 0) then
      print '(a)', 'FAIL (g) list_atoms_clear: members not deallocated / scalars not reset'
      nfail = nfail + 1
   end if
   call read_input('build/roster_case_mol_input.txt')   ! buffer re-read (same rule as (c))
   call sysdef_load('build/roster_case_mol')
   call list_atoms_load()
   dev = -1.0d0
   if (allocated(list_atoms%q0)) dev = maxval(abs(list_atoms%q0 - wq0_mol))
   if (.not. molecular_list_atoms_ok() .or. dev /= 0.0d0 .or. reactants%surface_model /= 0) then
      print '(a,es10.2)', 'FAIL (g) re-load: cross-load leakage after the surface system', dev
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(g) clear + re-load: PASS (full reset, molecular re-load leaves exactly the molecular list_atoms)'

   print '(a,i0,a)', '=== check_roster: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! field-by-field assertion of the assembled molecular list_atoms (sample of (a))
   logical function molecular_list_atoms_ok()
      molecular_list_atoms_ok = .false.
      if (.not. allocated(list_atoms%mass) .or. size(list_atoms%mass) /= 4) return
      if (.not. allocated(list_atoms%symb) .or. size(list_atoms%symb) /= 4) return
      if (.not. allocated(list_atoms%q0) .or. size(list_atoms%q0) /= 12) return
      if (.not. allocated(list_atoms%frag) .or. size(list_atoms%frag) /= 2) return
      if (reactants%n_frag /= 2 .or. reactants%surface_model /= 0) return
      if (list_atoms%symb(1) /= 'Au' .or. list_atoms%symb(2) /= 'Au' .or. &
          list_atoms%symb(3) /= 'C' .or. list_atoms%symb(4) /= 'O') return
      if (list_atoms%mass(1) /= elem_mass('Au') .or. list_atoms%mass(2) /= elem_mass('Au') .or. &
          list_atoms%mass(3) /= elem_mass('C') .or. list_atoms%mass(4) /= elem_mass('O')) return
      if (list_atoms%frag(1)%nat /= 2 .or. list_atoms%frag(2)%nat /= 2) return
      if (.not. equal_ints(list_atoms%frag(1)%list, [1,2]) .or. &
          .not. equal_ints(list_atoms%frag(2)%list, [3,4])) return
      if (list_atoms%frag(1)%mass /= elem_mass('Au') + elem_mass('Au') .or. &
          list_atoms%frag(2)%mass /= elem_mass('C') + elem_mass('O')) return
      if (trim(list_atoms%frag(1)%name) /= 'au_cluster' .or. &
          trim(list_atoms%frag(2)%name) /= 'co') return
      if (list_atoms%i_proj /= 2 .or. list_atoms%i_targ /= 1) return   ! PROJECTILE=co / TARGET=au_cluster
      if (.not. equal_reals(list_atoms%frag(1)%qz, [0.0d0,0.0d0,0.0d0,3.0d0,0.0d0,0.0d0])) return
      if (.not. equal_reals(list_atoms%frag(2)%qz, [0.0d0,0.0d0,0.0d0,1.20d0,0.0d0,0.0d0])) return
      if (list_atoms%skew /= 0.0d0 .or. list_atoms%box_lx /= 0.0d0 .or. list_atoms%box_ly /= 0.0d0) return
      molecular_list_atoms_ok = .true.
   end function molecular_list_atoms_ok

   ! field-by-field assertion of the assembled surface list_atoms (sample of (b))
   logical function surface_list_atoms_ok()
      surface_list_atoms_ok = .false.
      if (.not. allocated(list_atoms%mass) .or. size(list_atoms%mass) /= 3) return
      if (.not. allocated(list_atoms%symb) .or. size(list_atoms%symb) /= 3) return
      if (.not. allocated(list_atoms%q0) .or. size(list_atoms%q0) /= 9) return
      if (.not. allocated(list_atoms%frag) .or. size(list_atoms%frag) /= 2) return
      if (reactants%n_frag /= 2 .or. reactants%surface_model /= 1) return
      if (list_atoms%symb(1) /= 'Si' .or. list_atoms%symb(2) /= 'Si' .or. list_atoms%symb(3) /= 'H') return
      if (list_atoms%mass(1) /= elem_mass('Si') .or. list_atoms%mass(2) /= elem_mass('Si') .or. &
          list_atoms%mass(3) /= elem_mass('H')) return
      if (list_atoms%frag(1)%nat /= 2 .or. list_atoms%frag(2)%nat /= 1) return
      if (.not. equal_ints(list_atoms%frag(1)%list, [1,2]) .or. &
          .not. equal_ints(list_atoms%frag(2)%list, [3])) return
      if (list_atoms%frag(1)%mass /= 2.0d0*elem_mass('Si') .or. &
          list_atoms%frag(2)%mass /= elem_mass('H')) return
      if (.not. equal_reals(list_atoms%frag(1)%qz, [2.0d0,0.0d0,0.0d0,0.5d0,1.5d0,0.0d0])) return
      if (.not. equal_reals(list_atoms%frag(2)%qz, [0.0d0,0.0d0,6.0d0])) return
      if (trim(list_atoms%frag(1)%name) /= 'a_slab' .or. &
          trim(list_atoms%frag(2)%name) /= 'probe') return
      if (list_atoms%i_proj /= 2 .or. list_atoms%i_targ /= 0) return   ! PROJECTILE=probe; surface carries no target row
      if (list_atoms%box_lx /= w_lx .or. list_atoms%box_ly /= w_ly .or. abs(list_atoms%skew - w_skew) > 0.0d0) return
      surface_list_atoms_ok = .true.
   end function surface_list_atoms_ok

   ! D4 positive invariants on the assembled list_atoms (nat bounds / list range /
   ! mass positivity - see the negative-sample note in the header)
   logical function d4_invariants_hold()
      integer :: n, i, iat
      d4_invariants_hold = .false.
      if (.not. allocated(list_atoms%mass) .or. .not. allocated(list_atoms%frag)) return
      n = size(list_atoms%mass)
      if (n <= 0 .or. size(list_atoms%frag) <= 0) return
      do i = 1, size(list_atoms%frag)
         if (list_atoms%frag(i)%nat > n) return
         do iat = 1, list_atoms%frag(i)%nat
            if (list_atoms%frag(i)%list(iat) < 1 .or. list_atoms%frag(i)%list(iat) > n) return
         end do
      end do
      do iat = 1, n
         if (list_atoms%mass(iat) <= 0.0d0) return
      end do
      d4_invariants_hold = .true.
   end function d4_invariants_hold

   ! exact integer-array comparison helper
   logical function equal_ints(a, b)
      integer, intent(in) :: a(:), b(:)
      equal_ints = (size(a) == size(b))
      if (equal_ints) equal_ints = all(a == b)
   end function equal_ints

   ! exact real-array comparison helper
   logical function equal_reals(a, b)
      real(8), intent(in) :: a(:), b(:)
      equal_reals = (size(a) == size(b))
      if (equal_reals) equal_reals = all(a == b)
   end function equal_reals

   ! make one sample folder (idempotent - a second run tolerates its own leftovers)
   subroutine make_dir(path)
      character(len=*), intent(in) :: path
      if (windows_host()) then
         call execute_command_line('mkdir "'//trim(path)//'" 2>nul', wait=.true.)
      else
         call execute_command_line('mkdir -p "'//trim(path)//'" 2>/dev/null', wait=.true.)
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

   ! spawn the child for one bad case and assert the named abort: nonzero exit,
   ! stderr names the authority (list_atoms_load) and carries the expected key
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
      cmd = 'build'//sep//'check_roster.exe '//trim(cmode)//' 2>build'//sep//'roster_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/roster_err.txt', status='old', action='read', iostat=ios)
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

   ! Windows-host probe (same runtime OS branch the modules use for dir/ls)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_roster
