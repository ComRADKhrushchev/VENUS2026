!=====================================================================
! check_roster.f90 - C4.3 acceptance criteria program (criteria defined FIRST in
!                    docs/plans/C4.md; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files out
! of the source tree):
!   gfortran -Jbuild globals/config.f90 utils/specio.f90 facade/sysdef_facade.f90 &
!            globals/config_roster.f90 ../scripts/check_roster.f90 -o build/check_roster.exe
! All sample folders/files are written under build/ - the source tree stays clean.
! Loud-stop verification approach (same scheme as check_sysdef): the loud-stop
! criteria (d)/(e)/(f) run as CHILD PROCESSES - the check re-spawns itself with a
! mode argument (execute_command_line); the child runs the REAL
! sysdef_load -> roster_load chain on a prepared bad case and must print a naming
! message on stderr and exit nonzero. The parent asserts the exit status and the
! stderr content; a child that returns without stopping exits 0 and fails the
! criterion.
! Negative-sample note (D4 nat/list/mass rules; implementer's choice, stated):
!   the staging accessors are READ-ONLY and consistent by construction (nat runs
!   sum to the atom count; masses come from elem_mass, all positive) - no folder
!   and no public setter can produce a violating roster, so those three D4 rules
!   are tested as POSITIVE invariants on both assembled end-to-end rosters (c);
!   their loud-stop branches remain in roster_load as defensive guards for later
!   fill routes (the phase-C surface rows carry their own content).
! Criteria (C4.md C4.3 + D4):
!   (a) end-to-end molecular chain (acceptance 2): folder (au_cluster.xyz + co.xyz,
!       created OUT of name order) -> sysdef_load -> roster_load -> every roster
!       field exact: symb Au-Au-C-O; mass = elem_mass values exactly; frag(:) in
!       the E15 arrayed shape (2 homogeneous rows, nat 2/2, list = contiguous runs
!       1-2 / 3-4, mass = row sums, qz = packed geometry); q0 = the 3N packed
!       vector; no cell -> skew/box_lx/box_ly 0; box_pad 2.0; a_lat 0; size(frag)
!       = reactants%n_frag; surface_model 0
!   (b) end-to-end surface chain: skewed-cell POSCAR (a1=(4,0,0), a2=(1,3,0),
!       Cartesian mode - see the sample note below) + probe.xyz, roster%a_lat = 3
!       preset BEFORE roster_load (the input route stands in - a_lat is input-owned,
!       validated not derived) -> box_lx = 4, box_ly = sqrt(10),
!       skew = acos(1/sqrt(10)) exactly (law of cosines on the staged cell vectors);
!       a_lat carried through untouched (3.0); slab rows before probe (name order:
!       a_slab < probe)
!   (c) D4 invariants positive-path: on BOTH assembled rosters every fragment
!       nat <= atom count, every list index within 1..n_atoms, all masses positive
!   (d) guard loud-stop (child): roster_load with no prior sysdef_load -> nonzero
!       exit, stderr names roster_load and the assembly-order error
!   (e) D4 a_lat loud-stop (child): surface folder loaded clean, a_lat left at its
!       default 0 -> nonzero exit, stderr names roster_load and the missing
!       surface lattice constant (E12 ruling text: no generic default)
!   (f) degenerate-cell loud-stop (child): POSCAR with a zero in-plane lattice
!       vector -> nonzero exit, stderr names roster_load and the degenerate cell
!       (old SURF.f:74-79 statute)
!   (g) roster_clear + re-load: clear deallocates every member and resets every
!       scalar to its default (skew/box_lx/box_ly/a_lat 0, box_pad 2.0); a
!       molecular re-load after the surface load leaves exactly the molecular
!       roster - no cross-load leakage
!=====================================================================
program check_roster
   use config_roster, only: roster, roster_load, roster_clear
   use sysdef_facade, only: sysdef_load
   use config, only: reactants
   use specio, only: elem_mass
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

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('c-guard');  call roster_load()
      case ('c-alat');   call sysdef_load('build/roster_case_surf'); call roster_load()
      case ('c-degen');  call sysdef_load('build/roster_case_degen'); call roster_load()
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (roster_load returned - bug)'
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
   ! a latent C4.1 specio defect transposes Direct coordinates on non-symmetric cells
   ! (matmul(cell,fr) contracts the columns; the vectors are stored by row) - reported
   ! to the controller, not fixed here (C4.3 does not touch other modules); Cartesian
   ! mode sidesteps it while still staging the skewed cell for the law-of-cosines test
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

   ! ---------- (a) end-to-end molecular chain ----------
   nfail0 = nfail
   call sysdef_load('build/roster_case_mol')
   call roster_load()
   dev = -1.0d0
   if (allocated(roster%q0)) dev = maxval(abs(roster%q0 - wq0_mol))
   if (.not. molecular_roster_ok() .or. dev /= 0.0d0) then
      print '(a,es10.2)', 'FAIL (a) molecular chain: roster mismatch, max q0 deviation ', dev
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(a) molecular chain: PASS (symb/mass exact, 2 frag rows contiguous 1-2/3-4, q0 packed, no cell)'

   ! ---------- (b) end-to-end surface chain ----------
   nfail0 = nfail
   call sysdef_load('build/roster_case_surf')
   roster%a_lat = 3.0d0            ! the input route stands in (a_lat input-owned)
   call roster_load()
   dev = -1.0d0
   if (allocated(roster%q0)) dev = maxval(abs(roster%q0 - wq0_surf))
   if (.not. surface_roster_ok() .or. dev /= 0.0d0) then
      print '(a,es10.2)', 'FAIL (b) surface chain: roster mismatch, max q0 deviation ', dev
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(b) surface chain: PASS (box_lx=4, box_ly=sqrt(10), skew=acos(1/sqrt(10)) exact, a_lat carried)'

   ! ---------- (c) D4 invariants positive-path (both paradigms) ----------
   nfail0 = nfail
   if (.not. d4_invariants_hold()) then
      print '(a)', 'FAIL (c): D4 invariants violated on the assembled surface roster'
      nfail = nfail + 1
   end if
   call sysdef_load('build/roster_case_mol')
   call roster_load()
   if (.not. d4_invariants_hold()) then
      print '(a)', 'FAIL (c): D4 invariants violated on the assembled molecular roster'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(c) D4 invariants: PASS (nat<=n_atoms, list in range, masses positive - both rosters)'

   ! ---------- (d)-(f) loud-stop criteria (child processes) ----------
   call expect_loudstop('c-guard', 'sysdef_load must run first', 'roster_load', &
                        '(d) no-staging guard loud-stop')
   call expect_loudstop('c-alat',  'lattice constant', 'roster_load', &
                        '(e) missing-a_lat loud-stop (D4)')
   call expect_loudstop('c-degen', 'degenerate unit cell', 'roster_load', &
                        '(f) degenerate-cell loud-stop')

   ! ---------- (g) roster_clear + re-load ----------
   nfail0 = nfail
   call roster_clear()
   if (allocated(roster%mass) .or. allocated(roster%symb) .or. &
       allocated(roster%frag) .or. allocated(roster%q0) .or. &
       roster%skew /= 0.0d0 .or. roster%box_lx /= 0.0d0 .or. roster%box_ly /= 0.0d0 .or. &
       roster%box_pad /= 2.0d0 .or. roster%a_lat /= 0.0d0) then
      print '(a)', 'FAIL (g) roster_clear: members not deallocated / scalars not reset'
      nfail = nfail + 1
   end if
   call sysdef_load('build/roster_case_mol')
   call roster_load()
   dev = -1.0d0
   if (allocated(roster%q0)) dev = maxval(abs(roster%q0 - wq0_mol))
   if (.not. molecular_roster_ok() .or. dev /= 0.0d0 .or. reactants%surface_model /= 0) then
      print '(a,es10.2)', 'FAIL (g) re-load: cross-load leakage after the surface system', dev
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(g) clear + re-load: PASS (full reset, molecular re-load leaves exactly the molecular roster)'

   print '(a,i0,a)', '=== check_roster: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! field-by-field assertion of the assembled molecular roster (sample of (a))
   logical function molecular_roster_ok()
      molecular_roster_ok = .false.
      if (.not. allocated(roster%mass) .or. size(roster%mass) /= 4) return
      if (.not. allocated(roster%symb) .or. size(roster%symb) /= 4) return
      if (.not. allocated(roster%q0) .or. size(roster%q0) /= 12) return
      if (.not. allocated(roster%frag) .or. size(roster%frag) /= 2) return
      if (reactants%n_frag /= 2 .or. reactants%surface_model /= 0) return
      if (roster%symb(1) /= 'Au' .or. roster%symb(2) /= 'Au' .or. &
          roster%symb(3) /= 'C' .or. roster%symb(4) /= 'O') return
      if (roster%mass(1) /= elem_mass('Au') .or. roster%mass(2) /= elem_mass('Au') .or. &
          roster%mass(3) /= elem_mass('C') .or. roster%mass(4) /= elem_mass('O')) return
      if (roster%frag(1)%nat /= 2 .or. roster%frag(2)%nat /= 2) return
      if (.not. equal_ints(roster%frag(1)%list, [1,2]) .or. &
          .not. equal_ints(roster%frag(2)%list, [3,4])) return
      if (roster%frag(1)%mass /= elem_mass('Au') + elem_mass('Au') .or. &
          roster%frag(2)%mass /= elem_mass('C') + elem_mass('O')) return
      if (.not. equal_reals(roster%frag(1)%qz, [0.0d0,0.0d0,0.0d0,3.0d0,0.0d0,0.0d0])) return
      if (.not. equal_reals(roster%frag(2)%qz, [0.0d0,0.0d0,0.0d0,1.20d0,0.0d0,0.0d0])) return
      if (roster%skew /= 0.0d0 .or. roster%box_lx /= 0.0d0 .or. roster%box_ly /= 0.0d0 .or. &
          roster%box_pad /= 2.0d0 .or. roster%a_lat /= 0.0d0) return
      molecular_roster_ok = .true.
   end function molecular_roster_ok

   ! field-by-field assertion of the assembled surface roster (sample of (b))
   logical function surface_roster_ok()
      surface_roster_ok = .false.
      if (.not. allocated(roster%mass) .or. size(roster%mass) /= 3) return
      if (.not. allocated(roster%symb) .or. size(roster%symb) /= 3) return
      if (.not. allocated(roster%q0) .or. size(roster%q0) /= 9) return
      if (.not. allocated(roster%frag) .or. size(roster%frag) /= 2) return
      if (reactants%n_frag /= 2 .or. reactants%surface_model /= 1) return
      if (roster%symb(1) /= 'Si' .or. roster%symb(2) /= 'Si' .or. roster%symb(3) /= 'H') return
      if (roster%mass(1) /= elem_mass('Si') .or. roster%mass(2) /= elem_mass('Si') .or. &
          roster%mass(3) /= elem_mass('H')) return
      if (roster%frag(1)%nat /= 2 .or. roster%frag(2)%nat /= 1) return
      if (.not. equal_ints(roster%frag(1)%list, [1,2]) .or. &
          .not. equal_ints(roster%frag(2)%list, [3])) return
      if (roster%frag(1)%mass /= 2.0d0*elem_mass('Si') .or. &
          roster%frag(2)%mass /= elem_mass('H')) return
      if (.not. equal_reals(roster%frag(1)%qz, [2.0d0,0.0d0,0.0d0,0.5d0,1.5d0,0.0d0])) return
      if (.not. equal_reals(roster%frag(2)%qz, [0.0d0,0.0d0,6.0d0])) return
      if (roster%box_lx /= w_lx .or. roster%box_ly /= w_ly .or. abs(roster%skew - w_skew) > 0.0d0) return
      if (roster%a_lat /= 3.0d0) return      ! input-owned value carried through
      surface_roster_ok = .true.
   end function surface_roster_ok

   ! D4 positive invariants on the assembled roster (nat bounds / list range /
   ! mass positivity - see the negative-sample note in the header)
   logical function d4_invariants_hold()
      integer :: n, i, iat
      d4_invariants_hold = .false.
      if (.not. allocated(roster%mass) .or. .not. allocated(roster%frag)) return
      n = size(roster%mass)
      if (n <= 0 .or. size(roster%frag) <= 0) return
      do i = 1, size(roster%frag)
         if (roster%frag(i)%nat > n) return
         do iat = 1, roster%frag(i)%nat
            if (roster%frag(i)%list(iat) < 1 .or. roster%frag(i)%list(iat) > n) return
         end do
      end do
      do iat = 1, n
         if (roster%mass(iat) <= 0.0d0) return
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

   ! spawn the child for one bad case and assert the loud-stop: nonzero exit,
   ! stderr names the authority (roster_load) and carries the expected key
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
