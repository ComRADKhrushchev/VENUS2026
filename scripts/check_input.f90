!=====================================================================
! check_input.f90 - input-grammar criteria program (criteria defined FIRST in
!                    this header; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files out
! of the source tree; specio/sysdef_interface/config_atoms enter for the (o) window
! criteria, which drive the real assembly route input -> sysdef -> list_atoms):
!   gfortran -Jbuild globals/control.f90 globals/config.f90 utils/specio.f90 &
!            interface/sysdef_interface.f90 input.f90 globals/config_atoms.f90 &
!            ../scripts/check_input.f90 -o build/check_input.exe
! All sample input files are written under build/ - the source tree stays clean.
! Named-abort verification approach (same scheme as check_sysdef): the named abort
! criteria (e)/(f)/(g)/(h)/(k)/(l)/(o2)-(o4)/(p) run as CHILD PROCESSES - the
! check re-spawns itself with a mode argument (execute_command_line); the child
! calls the REAL read_input on a prepared bad file (the surface/molecular window
! children continue into sysdef_load + list_atoms_load - the paradigm condition lives
! at assembly time) and must print a naming message on stderr and
! exit nonzero. The parent asserts the exit status and the stderr content (the
! print+STOP behavior itself is verified, not a softened internal flag); a child
! that returns without stopping exits 0 and fails the criterion.
! Criteria:
!   (a) full key sweep - one good file carrying EVERY key of the first-version
!       table (both E21 blocks); every value lands in the exact slot:
!       control%(task/n_traj/title/i_seed/max_steps), propagator%(id/order/dt),
!       electronic%(method), observables%(rec_level/obs_list),
!       reactants%(system_dir/e_rel/b_max/r_sep/dist_scheme(:))
!       (DIST_SCHEME as a two-word list -> code array); ISEED stored as word 1 of
!       control%i_seed only (scalar rule - the 8-word conversion is rng-side)
!   (b) syntax mechanics (old input_parser semantics): blank lines and '#','!'
!       comment lines skipped; lowercase keys accepted (uppercased); trailing-
!       comma continuation merges the DIST_SCHEME list across lines - a messy
!       reordered rewrite of file (a) maps to the IDENTICAL slots
!   (c) per-key default injection (E22 - no preset track): a minimal file with
!       only the six required keys -> every Ruling-D default injected into its
!       slot (N_TRAJ=1, ISEED=1, MAX_STEPS=100000, TITLE='', TASK=TRAJECTORY,
!       INTEGRATOR=VERLET, ORDER=1, DT_MIN=0.0/DT_MAX=0.0 (unclamped - the
!       adaptive bounds' legal <=0 default, the reason both keys are kind_real
!       not kind_posreal), ELEC_METHOD='adiabatic',
!       GWRITE_LEVEL=2 -> rec_level (decision-fixed default -
!       old-library compatible, sec_qp on), OBSERVABLES=''; SCATTERING = GAS;
!       RANDOM_ORIENT off, BATH none); buffer holds 25 rows (6 user + 19 injected)
!   (d) defaults never override (set_default semantics): minimal + explicit
!       N_TRAJ/INTEGRATOR/ORDER/ISEED -> the user values win, untouched defaults
!       still injected
!   (e) missing required key (Ruling D - system-bound keys have no generic
!       default): each of SYSTEM_DIR/E_REL/R_SEP/DT/DIST_SCHEME removed in
!       turn -> named abort (child) naming the key. B_MAX left the required
!       set (2026-09-29): omission injects 0.0 - the surface paradigm's
!       on-site beam; the pair member enforces a positive radius at assembly
!   (f) unknown keys - the standalone grammar has NO fallback track (E21/E22):
!       MODEL (preset key withdrawn), old key names CALTYP/NACTA, withdrawn
!       event keys TERMINATE/OUTCOME (E18/E19), old member key TRV_A (no dispatch
!       row carries it - the semantic member keys replace it) -> named abort
!       (child) naming key + line number
!   (g) bad values -> named abort (child) naming key + line: bad integer N_TRAJ,
!       illegal enums TASK/INTEGRATOR/DIST_SCHEME, GWRITE_LEVEL outside {0,1,2,5},
!       bad logical RANDOM_ORIENT, bad real DT, bad real DT_MIN/DT_MAX (the over-width
!       reading for REAL keys: a scalar real carries no char value slot, so the
!       width column of criterion (p) has no member - the formal-domain stop a
!       wide garbage value hits is the bad-real parse itself)
!   (h) parse errors -> named abort (child): a data line without '=' (line number
!       reported; no sequential-grammar fallback), missing input file
!   (i) E21 two-block partition visible in the parse buffer: after the (a) load
!       every buffered row carries its key-table block (18 block-one rows, 5
!       block-two rows; line numbers of the source preserved)
   !   (j) member-key admission (the central dispatch table's rows): a file with
   !       the 13 channel-level member keys (N_E_REL/T_TRANS/V_WIDTH/N_THTA/
   !       THTA_MAX/N_CHI/CHI/N_B/N_AIM incident_surface/pair; AIM_X/AIM_Y
   !       incident_surface; T_BATH/NU_COLL bath_andersen) parses WITHOUT a stop
   !       and every key stages in block two with its value VERBATIM and the
   !       source line number kept. The 21 per-fragment member keys retired
   !       with the bracket ruling (s3) no longer stage
!   (k) Ruling F (the INCIDENT word withdrawal): DIST_SCHEME = incident ->
!       illegal-value named abort (child) naming DIST_SCHEME (the incident
!       channel is selected by the samp_run scheme argument, never by a
!       per-fragment scheme row)
!   (l) empty member-parameter value named abort (child): 'THTA_MAX =' with no
!       value -> named abort naming THTA_MAX + line
!   (m) the five new per-fragment scheme words parse to their codes:
!       DIST_SCHEME = rotation, j, surface_oscillator, normalmode,
!       barrier_excitation -> reactants%dist_scheme = [3,4,5,6,7] (the same
!       word/code pairing the sampler scheme-code table carries)
!   (n) the two dynamics-family scheme words parse to their codes:
!       DIST_SCHEME = glo_target, thermalize ->
!       reactants%dist_scheme = [8,9] (the word/code pairing shared with
!       the sampler scheme-code table, now ten words - thermostat left the
!       registry for the BATH provision 2026-09-22)
!   (o) scattering-class window keys (2026-10-13 ruling: the SCATTERING master
!       key + the CELL unit-cell key; the paradigm is an INPUT decision, the
!       system files carry geometry only): the window condition executes at
!       ASSEMBLY time (the sampler count-check precedent - the input layer
!       holds no system knowledge), so these criteria drive the real chain
!       read_input -> sysdef_load -> list_atoms_load:
!       (o1) surface input (SCATTERING = SURFACE + CELL + the role key) -> the
!            box lands in list_atoms (box_lx/box_ly/skew from the CELL rows),
!            surface_model = 1, the SCATTERING row tagged block two
!       (o2) SCATTERING = SURFACE with CELL removed -> named abort (child)
!            naming list_atoms_load and the missing required key CELL
!       (o3) SCATTERING = GAS (the default) with CELL present -> named abort
!            (child) naming list_atoms_load and the class conflict
!       (o4) malformed CELL (two rows instead of three) -> named abort
!            (child) at read_input naming key + line
!   (p) slot-width column (every char key registers its value-slot width = the
!       target field's width): the four families TITLE wider than 80,
!       ELEC_METHOD wider than 32, SYSTEM_DIR wider than 128, OBSERVABLES wider
!       than 256 -> named abort (child) naming key + line + the actual width
!       (silent truncation closed)
!   (q) collision-role keys (paradigm-conditional, assembly-side - the same
!       real-chain route as (o)): PROJECTILE/TARGET name list_atoms fragments by
!       their system-definition file stems (frag%name, case-insensitive fold):
!       (q1) gas pair + both keys -> list_atoms%i_proj/i_targ resolve to the named
!            rows, frag%name carries the stems (TARGET value lowercase - the
!            fold is case-insensitive)
!       (q2)/(q3) gas pair with PROJECTILE / TARGET removed -> named abort
!            (child) naming list_atoms_load, the missing key, and 'missing'
!       (q4) gas pair + unknown name -> named abort (child) naming list_atoms_load,
!            the key, and 'no list_atoms fragment named'
!       (q5) both keys naming the SAME fragment -> named abort (child) with
!            'distinct'
!       (q6) single-fragment molecular + PROJECTILE -> named abort (child)
!            'no collision roles'
!   (r) declared container keys (the E28 channel: declare -> admit -> pull
!       -> audit):
!       (r1) round-trip, in-process: a declared vocabulary (CK_REAL/CK_INT/
!            CK_CHR/CK_R3/CK_STRAY) admits the keys of a good file; typed
!            pulls land the values (real / int / char / 3-vector), the found
!            flag carries presence, a never-buffered key keeps its caller
!            default with found = F, and input_audit_declared passes clean
!       (r2) buffered-but-unconsumed: a declared key buffered but never pulled
!            -> input_audit_declared names the key (child)
!       (r3) undeclared pull: pulling a name outside the declared vocabulary
!            -> pull aborts 'not a declared container key' (child)
!       (r4) double pull: one declared key pulled twice -> pull aborts (child)
!       (r5) empty container value: 'CK_REAL =' with nothing -> read_input
!            aborts naming the key (child)
!       (q7) surface + TARGET -> named abort (child) 'gas-pair' (the surface is
!            not an atom list fragment)
!       (q8) surface missing PROJECTILE -> named abort (child) naming
!            list_atoms_load and 'missing'
   ! Sweep-along note: the map_slots key-kind select's case-default guard
   !       (the stale-ok leak fix, fourth site of the family) is UNREACHABLE
   !       from any input file - every key in the table carries kind 1..8, so
   !       no file can drive a non-listed kind into that branch; the guard
   !       completes the fixed family and is recorded here rather than
   !       exercised by a child.
   !   (s) per-fragment bracket parameters (2026-10-13 ruling: scheme-word
   !       parameters ride ONLY in brackets on the DIST_SCHEME item - no
   !       top-level member key, no fallback; a bracketed carrier is
   !       self-contained, a bare word rides the member's own field defaults):
   !       (s1) bracket form parses and stages per fragment: DIST_SCHEME =
   !            BOLTZMANN[T_VIB=350.0], EBK[N_VIB=1, N_ROT=2], ROTATION ->
   !            codes [1,2,3] and the pairs retrievable through the par
   !            accessors (typed casts, found flags, lowercase folded);
   !            continuation across lines still works (trailing comma)
   !       (s2) a bare word carries an empty pair set: pulls on it return
   !            found = F and keep the caller default (no fallback anywhere)
   !       (s3) the 21 retired per-fragment member keys are UNKNOWN KEYS as
   !            top-level lines (one child per scheme word: T_VIB boltzmann,
   !            T_ROT j, J_ROT rotation, N_OSC surface_oscillator, E_VIB
   !            normalmode, E_STAB barrier_excitation, T_GLO glo_target,
   !            T_EQ thermalize, E_SCALE wigner, N_VIB ebk) - the no-fallback
   !            statute is enforced by the key tables themselves
   !       (s4) malformed brackets abort naming DIST_SCHEME + the reason:
   !            unbalanced "[", empty "[]", a pair without "=", a duplicate
   !            parameter name, trailing text after "]", a nested "["
   !       (s5) bracket audit: a buffered pair the owning member never pulled
   !            -> input_audit_dist_par names word + fragment + parameter
   !            (child; mirrors the container-key audit)
   !       (s6) an empty parameter value "EBK[N_VIB=]" aborts at read_input
   !            naming the parameter (child)
   !       (s7) the value cast is member-side: par_int on "EBK[N_VIB=abc]"
   !            aborts naming the word, the fragment and the parameter (child)
!=====================================================================
program check_input
   use input,   only: read_input, buffer_n_rows, buffer_key, buffer_block, &
                      buffer_line, buffer_val, input_declare_keys, &
                      input_audit_declared, pull, par_int, par_real, &
                      input_audit_dist_par
   use control, only: task, n_traj, title, i_seed, max_steps
   use config,  only: propagator, electronic, reactants, observables
   use sysdef_interface, only: sysdef_load
   use config_atoms, only: list_atoms, list_atoms_load
   implicit none
   integer :: nfail = 0, nfail0, i, k, nrows
   character(len=32) :: mode
   character(len=32) :: kname
   logical :: ok_arr, ok_part
   character(len=16) :: blk1(13), blk2(12)
   character(len=16) :: mkey(13), mval(13)
   real(8) :: ck_r, ck_v3(3), ck_stray = 9.9d0   ! (r) pull targets (the stray
   integer :: ck_i                                ! default proves the absent-key
   character(len=32) :: ck_c                      ! posture)
   logical :: ck_found, ck_stray_found
   real(8) :: bp_r                                ! (s1) bracket-pair pull targets
   integer :: bp_i1, bp_i2, bp_j
   logical :: bp_found

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   ! bad-file naming convention: build/input_case_<mode>.txt (bad-nofile excepted);
   ! the w-surf-/w-mol- prefixed modes drive the full window chain (the paradigm
   ! condition lives at assembly time, so the child continues into list_atoms_load)
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      if (trim(mode) == 'bad-nofile') then
         call read_input('build/input_no_such_file.txt')
      else if (index(trim(mode), 'w-surf-') == 1) then
         call read_input('build/input_case_'//trim(mode)//'.txt')
         call sysdef_load('build/window_surf')
         call list_atoms_load()
      else if (index(trim(mode), 'w-mol-') == 1) then
         call read_input('build/input_case_'//trim(mode)//'.txt')
         call sysdef_load('build/window_mol')
         call list_atoms_load()
      else if (index(trim(mode), 'w-pr-') == 1) then
         call read_input('build/input_case_'//trim(mode)//'.txt')
         call sysdef_load('build/window_pair')
         call list_atoms_load()
      else if (index(trim(mode), 'ck-') == 1) then
         call child_ck(trim(mode))
      else if (trim(mode) == 'bp-audit') then
         call read_input('build/input_case_bp-audit.txt')
         call input_audit_dist_par()          ! a pair no member pull consumed
      else if (trim(mode) == 'bp-cast') then
         call read_input('build/input_case_bp-cast.txt')
         call par_int('EBK', 1, 'N_VIB', i)   ! the member-side value cast
      else
         call read_input('build/input_case_'//trim(mode)//'.txt')
      end if
      print '(a)', 'child: NO LOUD-STOP (the child chain returned - bug)'
      stop 0
   end if

   ! ---------- good sample files (all created under build/) ----------
   ! (a) full sweep: every key of the table, block one then block two; the
   ! DIST_SCHEME list continues onto the next line (trailing comma)
   call write_sample('build/input_case_full.txt', [ character(len=64) :: &
      'TASK = TRAJECTORY', &
      'N_TRAJ = 7', &
      'TITLE = full sweep sample', &
      'ISEED = 123456', &
      'MAX_STEPS = 5000', &
      'INTEGRATOR = symplectic', &
      'DT = 0.025', &
      'ORDER = 4', &
      'DT_MIN = 0.05', &
      'DT_MAX = 0.5', &
      'ELEC_METHOD = Adiabatic', &
      'GWRITE_LEVEL = 5', &
      'OBSERVABLES = t_kin,e_tot', &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 12.5', &
      'B_MAX = 3.0', &
      'R_SEP = 9.75', &
      'DIST_SCHEME = ebk,', &
      '   boltzmann', &
      'RANDOM_ORIENT = T', &
      'BATH = ANDERSEN', &
      'N_BATH = 5', &
      'DT_BATH = 0.025', &
      'SCATTERING = GAS' ])
   ! (b) messy rewrite of the same content: block two FIRST, lowercase keys,
   ! comments and blank lines interleaved, continuation kept
   call write_sample('build/input_case_messy.txt', [ character(len=64) :: &
      '# check_input messy-syntax sample (same content as file (a))', &
      'system_dir=build/sysdef_case_mol', &
      '', &
      'e_rel=12.5', &
      'b_max=3.0', &
      'r_sep=9.75', &
      'dist_scheme=ebk,', &
      '   boltzmann', &
      'random_orient=t', &
      'bath=andersen', &
      'n_bath=5', &
      'dt_bath=0.025', &
      'task=trajectory', &
      '! interleaved comment', &
      'n_traj=7', &
      'title=full sweep sample', &
      'iseed=123456', &
      'max_steps=5000', &
      'integrator=SYMPLECTIC', &
      'dt=0.025', &
      'order=4', &
      'dt_min=0.05', &
      'dt_max=0.5', &
      'elec_method=Adiabatic', &
      'gwrite_level=5', &
      'observables=t_kin,e_tot' ])
   ! (c) minimal file: the six required keys ONLY (Ruling D required set)
   call write_sample('build/input_case_min.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 1.0', &
      'B_MAX = 2.0', &
      'R_SEP = 8.0', &
      'DT = 0.01', &
      'DIST_SCHEME = boltzmann' ])
   ! (d) minimal + explicit user values that must WIN over the injected defaults
   call write_sample('build/input_case_over.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 1.0', &
      'B_MAX = 2.0', &
      'R_SEP = 8.0', &
      'DT = 0.02', &
      'DIST_SCHEME = ebk', &
      'N_TRAJ = 7', &
      'INTEGRATOR = radau', &
      'ORDER = 6', &
      'DT_MIN = 0.002', &
      'DT_MAX = 0.2', &
      'ISEED = 42' ])
   ! (j) the 13 channel-level member keys: minimal six + every member key
   ! (the 21 per-fragment member keys retired with the bracket ruling ride
   !  DIST_SCHEME brackets now - the (s3) children prove their rejection)
   call write_sample('build/input_case_mkeys.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 5.0', &
      'B_MAX = 3.0', &
      'R_SEP = 8.0', &
      'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'N_E_REL = 2', &
      'T_TRANS = 0.0', &
      'V_WIDTH = 0.3', &
      'N_THTA = 0', &
      'THTA_MAX = 1.1', &
      'N_CHI = 0', &
      'CHI = 0.9', &
      'N_B = 0', &
      'N_AIM = 0', &
      'AIM_X = 1.5', &
      'AIM_Y = 2.5', &
      'T_BATH = 475.0', &
      'NU_COLL = 0.137' ])
   ! (m) the five new per-fragment scheme words in one list (order kept)
   call write_sample('build/input_case_newdist.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 5.0', &
      'B_MAX = 3.0', &
      'R_SEP = 8.0', &
      'DT = 0.01', &
      'DIST_SCHEME = rotation, j, surface_oscillator,', &
      '   normalmode, barrier_excitation' ])
   ! (n) the three bath-member scheme words in one list (order kept)
   call write_sample('build/input_case_bathdist.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 5.0', &
      'B_MAX = 3.0', &
      'R_SEP = 8.0', &
      'DT = 0.01', &
      'DIST_SCHEME = glo_target, thermalize' ])
   ! (o) window-chain folders (dedicated names - the sysdef_case_* folders belong
   ! to check_sysdef's samples): a molecular window (single .xyz) and a surface
   ! window (a slab POSCAR) - since the 2026-10-13 ruling the folders are pure
   ! geometry (the class comes from the SCATTERING key); plus a gas-PAIR window
   ! for the (q) role criteria (two .xyz fragments)
   call make_dir('build/window_mol')
   call write_sample('build/window_mol/probe.xyz', [ character(len=64) :: &
      '1', 'window molecular sample', 'Ar 0.00 0.00 0.00' ])
   call make_dir('build/window_surf')
   call write_sample('build/window_surf/slab.poscar', [ character(len=64) :: &
      'window surface sample', '1.0', &
      '2.9 0.0 0.0', '0.0 2.9 0.0', '0.0 0.0 18.6', &
      'Au', '1', 'Direct', '0.00 0.00 0.10' ])
   call make_dir('build/window_pair')
   call write_sample('build/window_pair/frag_a.xyz', [ character(len=64) :: &
      '2', 'pair fragment A', 'H 0.00 0.00 0.00', 'H 0.74 0.00 0.00' ])
   call write_sample('build/window_pair/frag_b.xyz', [ character(len=64) :: &
      '2', 'pair fragment B', 'H 0.00 0.00 0.00', 'H 0.74 0.00 0.00' ])
   ! (o1) surface window: the six required keys + the class/cell keys + the role key
   call write_sample('build/input_case_surfwin.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'SCATTERING = SURFACE', &
      'CELL = 3.25 0 0, 0 3.25 0, 0 0 18.6', &
      'PROJECTILE = slab' ])
   ! (o2) surface class with the CELL key removed
   call write_sample('build/input_case_w-surf-nocell.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'SCATTERING = SURFACE' ])
   ! (o3) gas class (the SCATTERING default) with the CELL key present
   call write_sample('build/input_case_w-mol-cell.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_mol', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'CELL = 3.25 0 0, 0 3.25 0, 0 0 18.6' ])
   ! (q) role-key chain materials: gas-pair variants over window_pair, plus the
   ! surface / single-fragment-molecular negatives (w-surf-noproj reuses the
   ! full window set; window_reconcile passes, the role reconcile fires)
   call write_sample('build/input_case_w-pr-good.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_pair', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann, ebk', &
      'PROJECTILE = frag_a', 'TARGET = Frag_B' ])
   call write_sample('build/input_case_w-pr-noproj.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_pair', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann, ebk', &
      'TARGET = frag_b' ])
   call write_sample('build/input_case_w-pr-notarg.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_pair', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann, ebk', &
      'PROJECTILE = frag_a' ])
   call write_sample('build/input_case_w-pr-unknown.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_pair', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann, ebk', &
      'PROJECTILE = frag_x', 'TARGET = frag_b' ])
   call write_sample('build/input_case_w-pr-same.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_pair', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann, ebk', &
      'PROJECTILE = frag_a', 'TARGET = frag_a' ])
   call write_sample('build/input_case_w-mol-proj.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_mol', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', 'PROJECTILE = probe' ])
   call write_sample('build/input_case_w-surf-targ.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'SCATTERING = SURFACE', &
      'CELL = 3.25 0 0, 0 3.25 0, 0 0 18.6', &
      'PROJECTILE = slab', 'TARGET = slab' ])
   call write_sample('build/input_case_w-surf-noproj.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'SCATTERING = SURFACE', &
      'CELL = 3.25 0 0, 0 3.25 0, 0 0 18.6' ])
   ! (o4) malformed CELL (two rows instead of three; value-domain abort at read_input)
   call write_case('v-cell', 'CELL = 3.25 0 0, 0 3.25 0')
   ! (p) slot-width families: the offending key at line 1 with a value w
   ! characters wide, its minimal-six row replaced (a duplicate required key
   ! would named abort first and mask the width stop)
   call write_width_case('w-title', 'TITLE',       81)
   call write_width_case('w-meth', 'ELEC_METHOD',  33)
   call write_width_case('w-sysd', 'SYSTEM_DIR',  129)
   call write_width_case('w-obs',  'OBSERVABLES', 257)

   ! ---------- bad sample files (one per child named abort case) ----------
   ! (e) one required key removed each time
   call write_sample('build/input_case_miss-system_dir.txt', [ character(len=64) :: &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann' ])
   call write_sample('build/input_case_miss-e_rel.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', 'B_MAX = 2.0', 'R_SEP = 8.0', &
      'DT = 0.01', 'DIST_SCHEME = boltzmann' ])
   call write_sample('build/input_case_miss-b_max.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', 'E_REL = 1.0', 'R_SEP = 8.0', &
      'DT = 0.01', 'DIST_SCHEME = boltzmann' ])
   call write_sample('build/input_case_miss-r_sep.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', 'E_REL = 1.0', 'B_MAX = 2.0', &
      'DT = 0.01', 'DIST_SCHEME = boltzmann' ])
   call write_sample('build/input_case_miss-dt.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', 'E_REL = 1.0', 'B_MAX = 2.0', &
      'R_SEP = 8.0', 'DIST_SCHEME = boltzmann' ])
   call write_sample('build/input_case_miss-dist_scheme.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', 'E_REL = 1.0', 'B_MAX = 2.0', &
      'R_SEP = 8.0', 'DT = 0.01' ])
   ! (f)/(g)/(h) offending line FIRST (line 1), then the six required lines
   call write_case('u-model',     'MODEL = rigid_surface')
   call write_case('u-caltyp',    'CALTYP = -1')
   call write_case('u-nacta',     'NACTA = 2,1')
   call write_case('u-terminate', 'TERMINATE = 1000')
   call write_case('u-outcome',   'OUTCOME = table.txt')
   call write_case('u-trv_a',     'TRV_A = 300.0')
   call write_case('v-ntraj',     'N_TRAJ = abc')
   call write_case('v-task',      'TASK = BOGUS')
   call write_case('v-integ',     'INTEGRATOR = QUANTUM')
   call write_case('v-level',     'GWRITE_LEVEL = 9')
   call write_case('v-rorient',   'RANDOM_ORIENT = maybe')
   call write_case('p-noeq',      'N_TRAJ 7')
   ! DT / DIST_SCHEME bad values: the offending line REPLACES the minimal copy (those
   ! keys are part of the minimal six - a duplicate would named abort first)
   call write_sample('build/input_case_v-dist.txt', [ character(len=64) :: &
      'DIST_SCHEME = boltzmann, weird', 'SYSTEM_DIR = build/sysdef_case_mol', 'E_REL = 1.0', &
      'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01' ])
   call write_sample('build/input_case_v-dt.txt', [ character(len=64) :: &
      'DT = fast', 'SYSTEM_DIR = build/sysdef_case_mol', 'E_REL = 1.0', 'B_MAX = 2.0', &
      'R_SEP = 8.0', 'DIST_SCHEME = boltzmann' ])
   ! DT_MIN/DT_MAX bad reals (the real-key over-width stop - header note of (g))
   call write_case('v-dtmin', 'DT_MIN = narrow')
   call write_case('v-dtmax', 'DT_MAX = wide')
   ! (k) Ruling F: DIST_SCHEME = incident is an ILLEGAL per-fragment scheme word
   ! (offending line replaces the minimal copy - DIST_SCHEME is one of the six)
   call write_sample('build/input_case_v-incid.txt', [ character(len=64) :: &
      'DIST_SCHEME = incident', 'SYSTEM_DIR = build/sysdef_case_mol', 'E_REL = 1.0', &
      'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01' ])
   ! (l) a member key with an EMPTY value (line 1, then the minimal six)
   call write_case('m-empty', 'THTA_MAX =')
   ! (s3) the 21 retired per-fragment member keys are unknown as top-level lines
   ! (one representative per scheme word; the other 11 retired names share the
   ! same removed dispatch rows - N_ROT/K_OSC/T_OSC/N_LEVEL/N_E_BAR/E_BAR/T_BAR/
   ! GAMMA_GLO/W_GHOST/N_EQ/DT_EQ)
   call write_case('u-tvib',  'T_VIB = 350.0')
   call write_case('u-trot',  'T_ROT = 300.0')
   call write_case('u-jrot',  'J_ROT = 2')
   call write_case('u-nosc',  'N_OSC = 1')
   call write_case('u-evib',  'E_VIB = 12.5')
   call write_case('u-estab', 'E_STAB = 8.0')
   call write_case('u-tglo',  'T_GLO = 550.0')
   call write_case('u-teq',   'T_EQ = 325.0')
   call write_case('u-esc',   'E_SCALE = 1.0')
   call write_case('u-nvib',  'N_VIB = 1')
   ! (s4) malformed brackets / (s6) an empty parameter value: the DIST_SCHEME
   ! line REPLACES the minimal copy (a duplicate would named abort first)
   call write_dist_case('bp-unbal', 'DIST_SCHEME = EBK[N_VIB=1')
   call write_dist_case('bp-empty', 'DIST_SCHEME = EBK[]')
   call write_dist_case('bp-noeq',  'DIST_SCHEME = EBK[3]')
   call write_dist_case('bp-dup',   'DIST_SCHEME = EBK[N_VIB=1, N_VIB=2]')
   call write_dist_case('bp-trail', 'DIST_SCHEME = EBK[N_VIB=1] x')
   call write_dist_case('bp-nest',  'DIST_SCHEME = EBK[N_VIB=[1]]')
   call write_dist_case('bp-noval', 'DIST_SCHEME = EBK[N_VIB=]')
   ! (s5) a bracket parameter the owning member does not know (ebk carries
   ! N_VIB/N_ROT - T_VIB is boltzmann's) -> the audit names word+fragment+name
   call write_dist_case('bp-audit', 'DIST_SCHEME = EBK[T_VIB=300.0]')
   ! (s7) the value cast is member-side (par_int on a non-integer value)
   call write_dist_case('bp-cast', 'DIST_SCHEME = EBK[N_VIB=abc]')
   ! (s1)/(s2) the good bracket file: a bracketed boltzmann, a bracketed ebk
   ! (lowercase, two parameters, trailing comma continuing the list) and a
   ! BARE rotation word - three fragments, three codes
   call write_sample('build/input_case_brgood.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 1.0', &
      'B_MAX = 2.0', &
      'R_SEP = 8.0', &
      'DT = 0.01', &
      'DIST_SCHEME = boltzmann[T_VIB=350.0], ebk[n_vib=1, n_rot=2],', &
      '   rotation' ])

   ! ---------- (a) full key sweep ----------
   call read_input('build/input_case_full.txt')
   call assert_full('(a) full key sweep')

   ! ---------- (i) E21 two-block partition in the parse buffer ----------
   ! (run on the (a) buffer: 23 user rows + the 2 spectrum-key defaults injected)
   nfail0 = nfail
   blk1 = [ character(len=16) :: 'TASK', 'N_TRAJ', 'TITLE', 'ISEED', 'MAX_STEPS', &
            'INTEGRATOR', 'DT', 'ORDER', 'DT_MIN', 'DT_MAX', 'ELEC_METHOD', &
            'GWRITE_LEVEL', 'OBSERVABLES' ]
   blk2 = [ character(len=16) :: 'SYSTEM_DIR', 'E_REL', 'B_MAX', 'R_SEP', 'DIST_SCHEME', &
            'RANDOM_ORIENT', 'BATH', 'N_BATH', 'DT_BATH', 'SCATTERING', &
            'SPECTRUM_SOURCE', 'SPECTRUM_FILE' ]
   nrows = buffer_n_rows()
   ok_part = (nrows == 25)
   do i = 1, nrows
      kname = buffer_key(i)
      if (any(blk1 == trim(kname))) then
         if (buffer_block(i) /= 1) ok_part = .false.
      else if (any(blk2 == trim(kname))) then
         if (buffer_block(i) /= 2) ok_part = .false.
      else
         ok_part = .false.
      end if
   end do
   ! line numbers preserved into buffer (SYSTEM_DIR sits on source line 14)
   do i = 1, nrows
      if (trim(buffer_key(i)) == 'SYSTEM_DIR') exit
   end do
   if (i > nrows) ok_part = .false.
   if (ok_part .and. buffer_line(i) == 14) then
      print '(a)', '(i) E21 two-block buffer: PASS (13 block-one + 12 block-two rows tagged, line numbers kept)'
   else
      print '(a)', 'FAIL (i) E21 two-block buffer: partition or line numbers wrong'
      nfail = nfail + 1
   end if

   ! ---------- (b) syntax mechanics ----------
   call read_input('build/input_case_messy.txt')
   call assert_full('(b) syntax mechanics')

   ! ---------- (c) per-key default injection ----------
   nfail0 = nfail
   call read_input('build/input_case_min.txt')
   ok_arr = allocated(reactants%dist_scheme)
   if (ok_arr) ok_arr = size(reactants%dist_scheme) == 1 .and. &
                        reactants%dist_scheme(1) == 1
   if (task /= 1 .or. n_traj /= 1 .or. len_trim(title(1)) /= 0 .or. &
       i_seed(1) /= 1 .or. any(i_seed(2:8) /= 0) .or. max_steps /= 100000 .or. &
       propagator%id /= 1 .or. propagator%order /= 1 .or. propagator%dt /= 0.01d0 .or. &
       propagator%dt_min /= 0.0d0 .or. propagator%dt_max /= 0.0d0 .or. &
       trim(electronic%method) /= 'adiabatic' .or. observables%rec_level /= 2 .or. &
       len_trim(observables%obs_list) /= 0 .or. .not. ok_arr .or. &
       reactants%random_orient .or. reactants%bath /= 1 .or. reactants%n_bath /= 0 .or. &
       reactants%dt_bath /= 0.0d0 .or. &
       buffer_n_rows() /= 25) then
      print '(a)', 'FAIL (c) default injection: a default value is wrong (or buffer row count off 25)'
      nfail = nfail + 1
   else
      print '(a)', '(c) default injection: PASS (19 Ruling-D defaults injected, none preset-wide, buffer 6+19 rows)'
   end if

   ! ---------- (d) defaults never override user-set keys ----------
   nfail0 = nfail
   call read_input('build/input_case_over.txt')
   ok_arr = allocated(reactants%dist_scheme)
   if (ok_arr) ok_arr = size(reactants%dist_scheme) == 1 .and. &
                        reactants%dist_scheme(1) == 2
   if (n_traj /= 7 .or. propagator%id /= 3 .or. propagator%order /= 6 .or. &
       i_seed(1) /= 42 .or. propagator%dt /= 0.02d0 .or. &
       propagator%dt_min /= 0.002d0 .or. propagator%dt_max /= 0.2d0 .or. .not. ok_arr .or. &
       task /= 1 .or. max_steps /= 100000 .or. &
       trim(electronic%method) /= 'adiabatic' .or. observables%rec_level /= 2) then
      print '(a)', 'FAIL (d) no-override: user values lost or defaults clobbered them'
      nfail = nfail + 1
   else
      print '(a)', '(d) no-override: PASS (N_TRAJ=7/INTEGRATOR=radau/ORDER=6/DT_MIN=0.002/'// &
         'DT_MAX=0.2/ISEED=42 win; defaults hold elsewhere)'
   end if

   ! ---------- (e) missing required keys (children) ----------
   call expect_loudstop('miss-system_dir',   'SYSTEM_DIR',   '', '(e1) missing SYSTEM_DIR loud-stop')
   call expect_loudstop('miss-e_rel',        'E_REL',        '', '(e2) missing E_REL loud-stop')
   ! (e3) B_MAX omission contract: optional since 2026-09-29 - the injected
   ! default 0.0 is the surface on-site beam (pair enforces positivity later)
   call read_input('build/input_case_miss-b_max.txt')
   if (reactants%b_max == 0.0d0 .and. reactants%e_rel == 1.0d0) then
      print '(a)', '(e3) B_MAX omission: PASS (default 0.0 injected - the surface on-site beam)'
   else
      print '(a)', 'FAIL (e3) B_MAX omission: default 0.0 not injected'
      nfail = nfail + 1
   end if
   call expect_loudstop('miss-r_sep',        'R_SEP',        '', '(e4) missing R_SEP loud-stop')
   call expect_loudstop('miss-dt',           'DT',           '', '(e5) missing DT loud-stop')
   call expect_loudstop('miss-dist_scheme',  'DIST_SCHEME',  '', '(e6) missing DIST_SCHEME loud-stop')

   ! ---------- (f) unknown keys - no fallback track (children) ----------
   call expect_loudstop('u-model',     'MODEL',     'line 1', '(f1) MODEL intercepted as unknown key (E22)')
   call expect_loudstop('u-caltyp',    'CALTYP',    'line 1', '(f2) old key CALTYP intercepted (E21)')
   call expect_loudstop('u-nacta',     'NACTA',     'line 1', '(f3) old integer key NACTA intercepted (E21)')
   call expect_loudstop('u-terminate', 'TERMINATE', 'line 1', '(f4) withdrawn TERMINATE= stays out (E18)')
   call expect_loudstop('u-outcome',   'OUTCOME',   'line 1', '(f5) withdrawn OUTCOME= stays out (E19)')
   call expect_loudstop('u-trv_a',     'TRV_A',     'line 1', '(f6) old member key TRV_A intercepted (no dispatch row - semantic names)')

   ! ---------- (g) bad values (children) ----------
   call expect_loudstop('v-ntraj', 'N_TRAJ',       'line 1', '(g1) bad integer N_TRAJ')
   call expect_loudstop('v-task',  'TASK',         'line 1', '(g2) illegal enum TASK')
   call expect_loudstop('v-integ', 'INTEGRATOR',   'line 1', '(g3) illegal enum INTEGRATOR')
   call expect_loudstop('v-dist',  'DIST_SCHEME',  'line 1', '(g4) illegal enum item in DIST_SCHEME list')
   call expect_loudstop('v-level', 'GWRITE_LEVEL', 'line 1', '(g5) GWRITE_LEVEL outside {0,1,2,5}')
   call expect_loudstop('v-rorient', 'RANDOM_ORIENT', 'line 1', '(g6) bad logical RANDOM_ORIENT')
   call expect_loudstop('v-dt',    'DT',           'line 1', '(g7) bad real DT')
   call expect_loudstop('v-dtmin', 'DT_MIN',       'line 1', '(g8) bad real DT_MIN')
   call expect_loudstop('v-dtmax', 'DT_MAX',       'line 1', '(g9) bad real DT_MAX')

   ! ---------- (h) parse errors (children) ----------
   call expect_loudstop('p-noeq',     'line 1',  '', '(h1) no-= line loud-stop (line number reported)')
   call expect_loudstop('bad-nofile', 'not found', '', '(h2) missing input file loud-stop')

   ! ---------- (j) member-key admission (dispatch rows stage, verbatim) ----------
   nfail0 = nfail
   call read_input('build/input_case_mkeys.txt')
   mkey = [ character(len=16) :: 'N_E_REL', 'T_TRANS', 'V_WIDTH', &
            'N_THTA', 'THTA_MAX', 'N_CHI', 'CHI', 'N_B', 'N_AIM', 'AIM_X', 'AIM_Y', &
            'T_BATH', 'NU_COLL' ]
   mval = [ character(len=16) :: '2', '0.0', '0.3', &
            '0', '1.1', '0', '0.9', '0', '0', '1.5', '2.5', '475.0', '0.137' ]
   ok_part = (buffer_n_rows() == 38)      ! 19 user rows (six registry + the 13
                                            ! channel-level member keys) + the 19
                                            ! injected per-key defaults
   do i = 1, size(mkey)
      k = staged_row(trim(mkey(i)))
      if (k == 0 .or. buffer_block(k) /= 2 .or. trim(buffer_val(k)) /= trim(mval(i))) then
         ok_part = .false.
         exit
      end if
   end do
   k = staged_row('N_E_REL')              ! line numbers preserved (line 7 of the file)
   if (k == 0 .or. buffer_line(k) /= 7) ok_part = .false.
   if (ok_part) then
      print '(a)', '(j) member-key admission: PASS (13 channel-level member keys buffered '// &
         'block two, values verbatim)'
   else
      print '(a)', 'FAIL (j) member-key admission: a member key missing, wrong block, or value altered'
      nfail = nfail + 1
   end if

   ! ---------- (k)/(l) member-key aborts with a named error (children) ----------
   call expect_loudstop('v-incid',  'DIST_SCHEME', 'line 1', &
                        '(k) DIST_SCHEME=INCIDENT illegal value (Ruling F word withdrawal)')
   call expect_loudstop('m-empty',  'THTA_MAX',    'line 1', &
                        '(l) empty member-parameter value loud-stop')

   ! ---------- (m) the five new scheme words parse to codes 3..7 ----------
   nfail0 = nfail
   call read_input('build/input_case_newdist.txt')
   ok_arr = allocated(reactants%dist_scheme)
   if (ok_arr) ok_arr = size(reactants%dist_scheme) == 5 .and. &
                        all(reactants%dist_scheme == [ 3, 4, 5, 6, 7 ])
   if (ok_arr) then
      print '(a)', '(m) new scheme words: PASS (rotation/j/surface_oscillator/'// &
         'normalmode/barrier_excitation -> codes 3,4,5,6,7)'
   else
      print '(a)', 'FAIL (m) new scheme words: DIST_SCHEME list did not map to codes 3..7'
      nfail = nfail + 1
   end if

   ! ---------- (n) the dynamics-family scheme words parse to codes 8..9 ----------
   nfail0 = nfail
   call read_input('build/input_case_bathdist.txt')
   ok_arr = allocated(reactants%dist_scheme)
   if (ok_arr) ok_arr = size(reactants%dist_scheme) == 2 .and. &
                        all(reactants%dist_scheme == [ 8, 9 ])
   if (ok_arr) then
      print '(a)', '(n) dynamics-family scheme words: PASS (glo_target/'// &
         'thermalize -> codes 8,9; thermostat left for the BATH provision)'
   else
      print '(a)', 'FAIL (n) dynamics-family scheme words: DIST_SCHEME list did not map to codes 8..9'
      nfail = nfail + 1
   end if

   ! ---------- (o1) surface window chain: input -> sysdef -> list_atoms ----------
   nfail0 = nfail
   call read_input('build/input_case_surfwin.txt')
   call sysdef_load('build/window_surf')
   call list_atoms_load()
   k = staged_row('SCATTERING')
   if (abs(list_atoms%box_lx - 3.25d0) > 1.0d-12 .or. &
       abs(list_atoms%box_ly - 3.25d0) > 1.0d-12 .or. &
       abs(list_atoms%skew - acos(0.0d0)) > 1.0d-12 .or. &
       reactants%surface_model /= 1 .or. k == 0 .or. buffer_block(k) /= 2 .or. &
       list_atoms%i_proj /= 1 .or. trim(list_atoms%frag(1)%name) /= 'slab') then
      print '(a)', 'FAIL (o1) surface window chain: the box, role row, or fragment name did not land'
      nfail = nfail + 1
   else
      print '(a)', '(o1) surface window chain: PASS (CELL=3.25 orthogonal box -> '// &
         'box_lx=box_ly=3.25, skew=pi/2; PROJECTILE=slab -> i_proj=1, frag%name=slab;'// &
         ' surface_model=1, SCATTERING row tagged block two)'
   end if

   ! ---------- (o2)-(o4) window aborts with a named error (children) ----------
   call expect_loudstop_at('list_atoms_load', 'w-surf-nocell', 'CELL', 'missing', &
        '(o2) surface class missing CELL loud-stop (required under SCATTERING = SURFACE)')
   call expect_loudstop_at('list_atoms_load', 'w-mol-cell', 'CELL', 'SCATTERING = GAS', &
        '(o3) gas class + CELL loud-stop (the window key is surface-only)')
   call expect_loudstop('v-cell', 'CELL', 'line 1', '(o4) malformed CELL loud-stop (value domain)')

   ! ---------- (q) collision-role keys (one positive chain + children) ----------
   nfail0 = nfail
   call read_input('build/input_case_w-pr-good.txt')
   call sysdef_load('build/window_pair')
   call list_atoms_load()
   if (list_atoms%i_proj /= 1 .or. list_atoms%i_targ /= 2 .or. &
       trim(list_atoms%frag(1)%name) /= 'frag_a' .or. &
       trim(list_atoms%frag(2)%name) /= 'frag_b') then
      print '(a)', 'FAIL (q1) role chain: the role rows or fragment names did not resolve'
      nfail = nfail + 1
   else
      print '(a)', '(q1) role chain: PASS (PROJECTILE=frag_a / TARGET=Frag_B folded -> '// &
         'i_proj=1, i_targ=2; frag%name carries the file stems)'
   end if
   call expect_loudstop_at('list_atoms_load', 'w-pr-noproj', 'PROJECTILE', 'missing', &
        '(q2) gas pair missing PROJECTILE loud-stop (no positional default)')
   call expect_loudstop_at('list_atoms_load', 'w-pr-notarg', 'TARGET', 'missing', &
        '(q3) gas pair missing TARGET loud-stop')
   call expect_loudstop_at('list_atoms_load', 'w-pr-unknown', 'PROJECTILE', 'no list_atoms fragment named', &
        '(q4) unknown fragment name loud-stop (available names listed)')
   call expect_loudstop_at('list_atoms_load', 'w-pr-same', 'TARGET', 'distinct', &
        '(q5) PROJECTILE and TARGET naming the same fragment loud-stop')
   call expect_loudstop_at('list_atoms_load', 'w-mol-proj', 'collision roles', '', &
        '(q6) single-fragment molecular + PROJECTILE loud-stop (no collision roles)')
   call expect_loudstop_at('list_atoms_load', 'w-surf-targ', 'TARGET', 'gas-pair', &
        '(q7) surface + TARGET loud-stop (the surface is not an atom list fragment)')
   call expect_loudstop_at('list_atoms_load', 'w-surf-noproj', 'PROJECTILE', 'missing', &
        '(q8) surface missing PROJECTILE loud-stop')

   ! ---------- (p) slot-width families (children) ----------
   call expect_loudstop('w-title', 'TITLE',       'width 81',  '(p1) TITLE over its 80-wide slot loud-stop')
   call expect_loudstop('w-meth', 'ELEC_METHOD',  'width 33',  '(p2) ELEC_METHOD over its 32-wide slot loud-stop')
   call expect_loudstop('w-sysd', 'SYSTEM_DIR',   'width 129', '(p3) SYSTEM_DIR over its 128-wide slot loud-stop')
   call expect_loudstop('w-obs',  'OBSERVABLES',  'width 257', '(p4) OBSERVABLES over its 256-wide slot loud-stop')

   ! ---------- (r) declared container keys: declare -> admit -> pull -> audit ----------
   call input_declare_keys([character(len=16) :: 'CK_REAL', 'CK_INT', 'CK_CHR', &
                            'CK_R3', 'CK_STRAY'])
   call write_sample('build/input_case_ckgood.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 1.0', &
      'B_MAX = 2.0', &
      'R_SEP = 8.0', &
      'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'CK_REAL = 1.5', &
      'CK_INT = 7', &
      'CK_CHR = weights.txt', &
      'CK_R3 = 0.1, 0.2, 0.3' ])
   call read_input('build/input_case_ckgood.txt')
   ck_r = -1.0d0
   ck_i = -1
   ck_c = ''
   ck_v3 = -1.0d0
   call pull('CK_REAL', ck_r, ck_found)
   call pull('CK_INT', ck_i)
   call pull('CK_CHR', ck_c)
   call pull('CK_R3', ck_v3)
   call pull('CK_STRAY', ck_stray, ck_stray_found)   ! never buffered: default + found = F
   call input_audit_declared()                       ! every buffered declared key consumed
   if (ck_found .and. ck_r == 1.5d0 .and. ck_i == 7 .and. &
       trim(ck_c) == 'weights.txt' .and. &
       all(abs(ck_v3 - [0.1d0, 0.2d0, 0.3d0]) < 1.0d-12) .and. &
       .not. ck_stray_found .and. ck_stray == 9.9d0) then
      print '(a)', '(r1) declared-key channel: PASS (admission -> typed pulls -> audit clean)'
   else
      print '(a)', 'FAIL (r1) declared-key round-trip'
      nfail = nfail + 1
   end if
   call write_case_ck('ck-audit', 'CK_STRAY = 5.0')
   call expect_loudstop_at('input_audit_declared', 'ck-audit', 'CK_STRAY', '', &
      '(r2) buffered-but-unconsumed key audit abort')
   call write_case_ck('ck-undeclared', '')
   call expect_loudstop_at('pull', 'ck-undeclared', 'not a declared container key', '', &
      '(r3) undeclared pull abort')
   call write_case_ck('ck-dbl', 'CK_REAL = 1.5')
   call expect_loudstop_at('pull', 'ck-dbl', 'twice', '', '(r4) double-pull abort')
   call write_case_ck('ck-empty', 'CK_REAL =')
   call expect_loudstop('ck-empty', 'CK_REAL', 'empty container-parameter', &
        '(r5) empty container-value abort')

   ! ---------- (s3) retired per-fragment member keys are unknown (children) ----------
   call expect_loudstop('u-tvib',  'T_VIB',    'unknown key', &
        '(s3a) retired T_VIB (boltzmann) intercepted as unknown key')
   call expect_loudstop('u-trot',  'T_ROT',    'unknown key', &
        '(s3b) retired T_ROT (j) intercepted as unknown key')
   call expect_loudstop('u-jrot',  'J_ROT',    'unknown key', &
        '(s3c) retired J_ROT (rotation) intercepted as unknown key')
   call expect_loudstop('u-nosc',  'N_OSC',    'unknown key', &
        '(s3d) retired N_OSC (surface_oscillator) intercepted as unknown key')
   call expect_loudstop('u-evib',  'E_VIB',    'unknown key', &
        '(s3e) retired E_VIB (normalmode) intercepted as unknown key')
   call expect_loudstop('u-estab', 'E_STAB',   'unknown key', &
        '(s3f) retired E_STAB (barrier_excitation) intercepted as unknown key')
   call expect_loudstop('u-tglo',  'T_GLO',    'unknown key', &
        '(s3g) retired T_GLO (glo_target) intercepted as unknown key')
   call expect_loudstop('u-teq',   'T_EQ',     'unknown key', &
        '(s3h) retired T_EQ (thermalize) intercepted as unknown key')
   call expect_loudstop('u-esc',   'E_SCALE',  'unknown key', &
        '(s3i) retired E_SCALE (wigner) intercepted as unknown key')
   call expect_loudstop('u-nvib',  'N_VIB',    'unknown key', &
        '(s3j) retired N_VIB (ebk) intercepted as unknown key')

   ! ---------- (s4)/(s6) malformed brackets (children) ----------
   call expect_loudstop('bp-unbal', 'DIST_SCHEME', 'unbalanced', &
        '(s4a) unbalanced bracket abort')
   call expect_loudstop('bp-empty', 'DIST_SCHEME', 'bare', &
        '(s4b) empty bracket abort (a word without parameters is written bare)')
   call expect_loudstop('bp-noeq',  'DIST_SCHEME', 'NAME=VALUE', &
        '(s4c) a bracket item without = aborts (parameters are NAME=VALUE)')
   call expect_loudstop('bp-dup',   'DIST_SCHEME', 'duplicate', &
        '(s4d) a duplicate parameter name aborts (one value per carrier)')
   call expect_loudstop('bp-trail', 'DIST_SCHEME', 'trailing', &
        '(s4e) trailing text after the closing bracket aborts')
   call expect_loudstop('bp-nest',  'DIST_SCHEME', 'nested', &
        '(s4f) a nested bracket aborts')
   call expect_loudstop('bp-noval', 'DIST_SCHEME', 'empty parameter', &
        '(s6) an empty parameter value aborts at read_input')

   ! ---------- (s1)/(s2) bracket form stages per fragment ----------
   nfail0 = nfail
   call read_input('build/input_case_brgood.txt')
   ok_arr = allocated(reactants%dist_scheme)
   if (ok_arr) ok_arr = size(reactants%dist_scheme) == 3 .and. &
                        all(reactants%dist_scheme == [ 1, 2, 3 ])
   bp_r = -1.0d0
   bp_i1 = -1
   bp_i2 = -1
   call par_real('BOLTZMANN', 1, 'T_VIB', bp_r)          ! typed cast, uppercase fold
   call par_int('EBK', 2, 'n_vib', bp_i1)                ! lowercase parameter name folds
   call par_int('EBK', 2, 'N_ROT', bp_i2)
   bp_j = -1
   call par_int('ROTATION', 3, 'J_ROT', bp_j, bp_found)  ! bare word: no fallback
   call input_audit_dist_par()                           ! every staged pair consumed
   if (ok_arr .and. bp_r == 350.0d0 .and. bp_i1 == 1 .and. bp_i2 == 2 .and. &
       .not. bp_found .and. bp_j == -1) then
      print '(a)', '(s1)+(s2) bracket staging: PASS (codes [1,2,3]; BOLTZMANN T_VIB=350,'// &
         ' EBK N_VIB=1/N_ROT=2 cast per carrier; bare ROTATION keeps the caller default)'
   else
      print '(a)', 'FAIL (s1)+(s2) bracket staging: codes, casts, or the bare-word default wrong'
      nfail = nfail + 1
   end if

   ! ---------- (s5)/(s7) bracket seam aborts (children) ----------
   call expect_loudstop_at('input_audit_dist_par', 'bp-audit', 'EBK', 'T_VIB', &
        '(s5) unconsumed bracket parameter audit abort (names word + fragment)')
   call expect_loudstop_at('par', 'bp-cast', 'N_VIB', 'bad integer', &
        '(s7) the member-side value cast aborts naming word + parameter')

   print '(a,i0,a)', '=== check_input: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! staged_row(key) - buffer row index of a key (0 = not buffered; a scan over the
   ! public accessors, mirroring what the assembly wiring does to find member keys)
   function staged_row(key) result(idx)
      character(len=*), intent(in) :: key
      integer :: idx, j
      idx = 0
      do j = 1, buffer_n_rows()
         if (trim(buffer_key(j)) == trim(key)) then
            idx = j
            return
         end if
      end do
   end function staged_row

   ! the slot assertions shared by (a) and (b): identical content, different syntax
   subroutine assert_full(label)
      character(len=*), intent(in) :: label
      nfail0 = nfail
      ok_arr = allocated(reactants%dist_scheme)
      if (ok_arr) ok_arr = size(reactants%dist_scheme) == 2 .and. &
                           all(reactants%dist_scheme == (/ 2, 1 /))
      if (task /= 1 .or. n_traj /= 7 .or. &
          trim(title(1)) /= 'full sweep sample' .or. len_trim(title(2)) /= 0 .or. &
          i_seed(1) /= 123456 .or. any(i_seed(2:8) /= 0) .or. max_steps /= 5000 .or. &
          propagator%id /= 2 .or. propagator%order /= 4 .or. propagator%dt /= 0.025d0 .or. &
          propagator%dt_min /= 0.05d0 .or. propagator%dt_max /= 0.5d0 .or. &
          trim(electronic%method) /= 'adiabatic' .or. observables%rec_level /= 5 .or. &
          trim(observables%obs_list) /= 't_kin,e_tot' .or. &
          trim(reactants%system_dir) /= 'build/sysdef_case_mol' .or. &
          reactants%e_rel /= 12.5d0 .or. reactants%b_max /= 3.0d0 .or. &
          reactants%r_sep /= 9.75d0 .or. .not. reactants%random_orient .or. &
          reactants%bath /= 2 .or. reactants%n_bath /= 5 .or. reactants%dt_bath /= 0.025d0 .or. &
          .not. ok_arr) then
         print '(a)', 'FAIL '//trim(label)//': slot mismatch'
         nfail = nfail + 1
      else
         print '(a)', trim(label)//': PASS (27 keys -> exact slots; lowercase words, ISEED word-1-only)'
      end if
   end subroutine assert_full

   ! write one bad case file: offending line first (line 1), minimal six after
   ! (a LOCAL array, not an array-constructor actual: a constructor mixing a len=*
   ! dummy with literals trips a heap-corruption at the temp free on this gfortran
   ! 15.2/Windows target)
   subroutine write_case(cmode, first)
      character(len=*), intent(in) :: cmode, first
      character(len=64) :: lines(7), fname
      lines(1) = first
      lines(2) = 'SYSTEM_DIR = build/sysdef_case_mol'
      lines(3) = 'E_REL = 1.0'
      lines(4) = 'B_MAX = 2.0'
      lines(5) = 'R_SEP = 8.0'
      lines(6) = 'DT = 0.01'
      lines(7) = 'DIST_SCHEME = boltzmann'
      fname = 'build/input_case_'//trim(cmode)//'.txt'
      call write_sample(fname, lines)
   end subroutine write_case

   ! write one bracket-family bad case: the DIST_SCHEME line itself is the
   ! offender and REPLACES the minimal copy (a second DIST_SCHEME row would
   ! named abort as a duplicate before the bracket stop under test)
   subroutine write_dist_case(cmode, distline)
      character(len=*), intent(in) :: cmode, distline
      character(len=64) :: lines(6), fname
      lines(1) = 'SYSTEM_DIR = build/sysdef_case_mol'
      lines(2) = 'E_REL = 1.0'
      lines(3) = 'B_MAX = 2.0'
      lines(4) = 'R_SEP = 8.0'
      lines(5) = 'DT = 0.01'
      lines(6) = distline
      fname = 'build/input_case_'//trim(cmode)//'.txt'
      call write_sample(fname, lines)
   end subroutine write_dist_case

   ! child_ck(mode) - the (r) children: declare the test vocabulary, read the
   !                 case input, then drive the one bad channel of the mode
   subroutine child_ck(mode)
      character(len=*), intent(in) :: mode
      real(8) :: v
      call input_declare_keys([character(len=16) :: 'CK_REAL', 'CK_INT', 'CK_CHR', &
                               'CK_R3', 'CK_STRAY'])
      call read_input('build/input_case_'//trim(mode)//'.txt')
      select case (trim(mode))
      case ('ck-audit')
         call input_audit_declared()          ! CK_STRAY buffered, never pulled
      case ('ck-undeclared')
         call pull('CK_NOPE', v)              ! outside the declared vocabulary
      case ('ck-dbl')
         call pull('CK_REAL', v)
         call pull('CK_REAL', v)              ! one key, one consumer
      case ('ck-empty')
         continue                             ! the abort fired inside read_input
      end select
      print '(a)', 'child: NO LOUD-STOP (the ck chain returned - bug)'
      stop 0
   end subroutine child_ck

   ! write_case_ck(cmode, extra) - the six required keys + one optional extra
   !                 line (an empty extra leaves a blank line the parser skips)
   subroutine write_case_ck(cmode, extra)
      character(len=*), intent(in) :: cmode, extra
      character(len=64) :: lines(7)
      lines(1) = 'SYSTEM_DIR = build/sysdef_case_mol'
      lines(2) = 'E_REL = 1.0'
      lines(3) = 'B_MAX = 2.0'
      lines(4) = 'R_SEP = 8.0'
      lines(5) = 'DT = 0.01'
      lines(6) = 'DIST_SCHEME = boltzmann'
      lines(7) = extra
      call write_sample('build/input_case_'//trim(cmode)//'.txt', lines)
   end subroutine write_case_ck

   ! write one sample file, one array element per record (samples live under build/)
   subroutine write_sample(fname, lines)
      character(len=*), intent(in) :: fname
      character(len=*), intent(in) :: lines(:)
      integer :: v, j
      open (newunit=v, file=fname, status='replace', action='write')
      do j = 1, size(lines)
         write (v, '(a)') trim(lines(j))
      end do
      close (v)
   end subroutine write_sample

   ! spawn the child for one bad case and assert the named abort: nonzero exit,
   ! stderr names the authority (read_input) and carries the expected key /
   ! extra needle (an empty needle is skipped)
   subroutine expect_loudstop(cmode, key, needle, label)
      character(len=*), intent(in) :: cmode, key, needle, label
      call expect_loudstop_at('read_input', cmode, key, needle, label)
   end subroutine expect_loudstop

   ! the authority-parameterized form: the window criteria stop inside
   ! list_atoms_load (the paradigm condition is assembly-side), so the asserted
   ! stderr authority is the assembly face, not read_input
   subroutine expect_loudstop_at(who, cmode, key, needle, label)
      character(len=*), intent(in) :: who, cmode, key, needle, label
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
      cmd = 'build'//sep//'check_input.exe '//trim(cmode)//' 2>build'//sep//'input_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/input_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) line
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(line)
         end do
         close (v)
      end if
      if (stat == 0) then
         print '(a)', 'FAIL '//trim(label)//': expected a loud-stop, the child exited 0'
         nfail = nfail + 1
      else if (index(buf, trim(who)) == 0) then
         print '(a)', 'FAIL '//trim(label)//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (len_trim(key) > 0 .and. index(buf, trim(key)) == 0) then
         print '(a)', 'FAIL '//trim(label)//': stderr missing the key "'//trim(key)//'"'
         nfail = nfail + 1
      else if (len_trim(needle) > 0 .and. index(buf, trim(needle)) == 0) then
         print '(a)', 'FAIL '//trim(label)//': stderr missing "'//trim(needle)//'"'
         nfail = nfail + 1
      else
         print '(a)', trim(label)//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop_at

   ! write one width-family bad case: the offending key at line 1 with a value w
   ! characters wide, then the minimal six with the offending key's own row
   ! replaced (a duplicate required key would named abort first and mask the stop
   ! under test)
   subroutine write_width_case(cmode, key, w)
      character(len=*), intent(in) :: cmode, key
      integer, intent(in) :: w
      character(len=256) :: six(6), keep(6), fname
      character(len=300) :: offend
      integer :: n, m, v
      offend = trim(key)//' = '//repeat('x', w)
      six = [ character(len=256) :: 'SYSTEM_DIR = build/sysdef_case_mol', 'E_REL = 1.0', &
             'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', 'DIST_SCHEME = boltzmann' ]
      m = 0
      do n = 1, 6
         if (index(six(n), trim(key)//' =') == 1) cycle    ! its row is replaced
         m = m + 1
         keep(m) = six(n)
      end do
      fname = 'build/input_case_'//trim(cmode)//'.txt'
      open (newunit=v, file=fname, status='replace', action='write')
      write (v, '(a)') trim(offend)
      do n = 1, m
         write (v, '(a)') trim(keep(n))
      end do
      close (v)
   end subroutine write_width_case

   ! make one sample folder (idempotent - a second run tolerates its own leftovers)
   subroutine make_dir(path)
      character(len=*), intent(in) :: path
      if (windows_host()) then
         call execute_command_line('mkdir "'//trim(path)//'" 2>nul', wait=.true.)
      else
         call execute_command_line('mkdir -p "'//trim(path)//'" 2>/dev/null', wait=.true.)
      end if
   end subroutine make_dir

   ! Windows-host probe (same runtime OS branch check_sysdef uses for redirects)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_input
