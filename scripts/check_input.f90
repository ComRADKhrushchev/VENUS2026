!=====================================================================
! check_input.f90 - C4.4 acceptance criteria program (criteria defined FIRST in
!                    docs/plans/C4.md; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files out
! of the source tree; specio/sysdef_facade/config_roster enter for the (o) window
! criteria, which drive the real assembly route input -> sysdef -> roster):
!   gfortran -Jbuild globals/control.f90 globals/config.f90 utils/specio.f90 &
!            facade/sysdef_facade.f90 input.f90 globals/config_roster.f90 &
!            ../scripts/check_input.f90 -o build/check_input.exe
! All sample input files are written under build/ - the source tree stays clean.
! Loud-stop verification approach (same scheme as check_sysdef): the loud-stop
! criteria (e)/(f)/(g)/(h)/(k)/(l)/(o2)-(o7)/(p) run as CHILD PROCESSES - the
! check re-spawns itself with a mode argument (execute_command_line); the child
! calls the REAL read_input on a prepared bad file (the surface/molecular window
! children continue into sysdef_load + roster_load - the paradigm condition lives
! at assembly time) and must print a naming message on stderr and
! exit nonzero. The parent asserts the exit status and the stderr content (the
! print+STOP behavior itself is verified, not a softened internal flag); a child
! that returns without stopping exits 0 and fails the criterion.
! Criteria (C4.md C4.4 step + task-4 brief Ruling D):
!   (a) full key sweep - one good file carrying EVERY key of the first-version
!       table (both E21 blocks); every value lands in the exact slot:
!       control%(task/n_traj/title/i_seed/max_steps), propagator%(id/order/dt),
!       electronic%(method/n_surf/pes_main), observables%(rec_level/obs_list/
!       e_kin/e_pot/e_tot), reactants%(system_dir/e_rel/b_max/r_sep/dist_scheme(:))
!       (DIST_SCHEME as a two-word list -> code array); ISEED stored as word 1 of
!       control%i_seed only (scalar statute - the 8-word conversion is rng-side)
!   (b) syntax mechanics (old input_parser semantics): blank lines and '#','!'
!       comment lines skipped; lowercase keys accepted (uppercased); trailing-
!       comma continuation merges the DIST_SCHEME list across lines - a messy
!       reordered rewrite of file (a) maps to the IDENTICAL slots
!   (c) per-key default injection (E22 - no preset track): a minimal file with
!       only the six required keys -> every Ruling-D default injected into its
!       slot (N_TRAJ=1, ISEED=1, MAX_STEPS=100000, TITLE='', TASK=TRAJECTORY,
!       INTEGRATOR=VERLET, ORDER=1, DT_MIN=0.0/DT_MAX=0.0 (unclamped - the
!       adaptive bounds' legal <=0 default, the reason both keys are kind_real
!       not kind_posreal), ELEC_METHOD='adiabatic', N_SURF=0,
!       PES_MAIN='1', GWRITE_LEVEL=2 -> rec_level (E26 legislated default -
!       old-library compatible, sec_qp on), OBSERVABLES='' and
!       E_KIN/E_POT/E_TOT off); staging holds 23 rows (6 user + 17 injected)
!   (d) defaults never override (set_default semantics): minimal + explicit
!       N_TRAJ/INTEGRATOR/ORDER/ISEED -> the user values win, untouched defaults
!       still injected
!   (e) missing required key (Ruling D - system-bound keys have no generic
!       default): each of SYSTEM_DIR/E_REL/B_MAX/R_SEP/DT/DIST_SCHEME removed in
!       turn -> loud-stop (child) naming the key
!   (f) unknown keys - the standalone grammar has NO fallback track (E21/E22):
!       MODEL (preset key withdrawn), old key names CALTYP/NACTA, withdrawn
!       event keys TERMINATE/OUTCOME (E18/E19), old member key TRV_A (no dispatch
!       row carries it - the semantic member keys replace it) -> loud-stop
!       (child) naming key + line number
!   (g) bad values -> loud-stop (child) naming key + line: bad integer N_TRAJ,
!       illegal enums TASK/INTEGRATOR/DIST_SCHEME, GWRITE_LEVEL outside {0,1,2,5},
!       bad logical E_KIN, bad real DT, bad real DT_MIN/DT_MAX (the over-width
!       reading for REAL keys: a scalar real carries no char value slot, so the
!       width column of criterion (p) has no member - the formal-domain stop a
!       wide garbage value hits is the bad-real parse itself)
!   (h) parse errors -> loud-stop (child): a data line without '=' (line number
!       reported; no sequential-grammar fallback), missing input file
!   (i) E21 two-block partition visible in the parse staging: after the (a) load
!       every staged row carries its key-table block (18 block-one rows, 5
!       block-two rows; line numbers of the source preserved)
!   (j) member-key admission (the central dispatch table's rows): a file with
!       all 30 member keys (T_VIB_A/T_VIB_B boltzmann; N_E_REL/T_TRANS/V_WIDTH/
!       N_THTA/THTA_MAX/N_CHI/CHI/N_B/N_AIM incident; T_ROT j; J_ROT rotation;
!       N_OSC/T_OSC/N_LEVEL/K_OSC surface_oscillator; E_VIB normalmode;
!       E_STAB/N_E_BAR/E_BAR/T_BAR barrier_excitation; T_BATH/NU_COLL
!       thermostat; T_GLO/GAMMA_GLO/W_GHOST glo_target; T_EQ/N_EQ/DT_EQ
!       thermalize) parses WITHOUT a stop - no longer unknown keys - and every
!       key stages in block two with its value VERBATIM (the raw trimmed
!       string: '350.0', '2', ...) and the source line number kept
!   (k) Ruling F (the INCIDENT word withdrawal): DIST_SCHEME = incident ->
!       illegal-value loud-stop (child) naming DIST_SCHEME (the incident
!       channel is selected by the samp_run scheme argument, never by a
!       per-fragment scheme row)
!   (l) empty member-parameter value loud-stop (child): 'THTA_MAX =' with no
!       value -> loud-stop naming THTA_MAX + line
!   (m) the five new per-fragment scheme words parse to their codes:
!       DIST_SCHEME = rotation, j, surface_oscillator, normalmode,
!       barrier_excitation -> reactants%dist_scheme = [3,4,5,6,7] (the same
!       word/code pairing the sampler scheme-code table carries)
!   (n) the three bath-member scheme words parse to their codes:
!       DIST_SCHEME = thermostat, glo_target, thermalize ->
!       reactants%dist_scheme = [8,9,10] (the word/code pairing shared with
!       the sampler scheme-code table, now ten words)
!   (o) surface window keys (paradigm-conditional): A_LAT [A] (positive real,
!       the surface lattice constant) and BOXPAD [A] (positive, the periodic-box
!       fill half-width) stage in block two; the paradigm condition executes at
!       ASSEMBLY time (the sampler count-check precedent - the input layer holds
!       no system knowledge), so these criteria drive the real chain
!       read_input -> sysdef_load -> roster_load:
!       (o1) surface folder + both keys -> the values land in roster%a_lat /
!            roster%box_pad (the registered slot path: config carries no relay
!            field, the values stay in the input staging and roster_load pulls
!            them through the staging accessors - the S2-shaped transfer),
!            surface_model = 1, the A_LAT row tagged block two
!       (o2) surface folder + A_LAT removed -> loud-stop (child) naming
!            roster_load and the missing required key A_LAT
!       (o3) surface folder + BOXPAD removed -> loud-stop (child) naming
!            roster_load and BOXPAD
!       (o4)/(o5) molecular folder + A_LAT / BOXPAD present -> loud-stop
!            (child) naming roster_load and the paradigm conflict (the
!            molecular paradigm has no surface window)
!       (o6)/(o7) non-positive A_LAT / BOXPAD -> loud-stop (child) at
!            read_input naming key + line (the positive-real value domain)
!   (p) slot-width column (every char key registers its value-slot width = the
!       target field's width): the four families TITLE wider than 80,
!       ELEC_METHOD wider than 32, SYSTEM_DIR wider than 128, OBSERVABLES wider
!       than 256 -> loud-stop (child) naming key + line + the actual width
!       (silent truncation closed; PES_MAIN shares the 32 family)
! Sweep-along note: the map_slots key-kind select's case-default guard
!       (the stale-ok leak fix, fourth site of the family) is UNREACHABLE
!       from any input file - every key in the table carries kind 1..7, so
!       no file can drive a non-listed kind into that branch; the guard
!       completes the fixed family and is recorded here rather than
!       exercised by a child.
!=====================================================================
program check_input
   use input,   only: read_input, staging_n_rows, staging_key, staging_block, &
                      staging_line, staging_val
   use control, only: task, n_traj, title, i_seed, max_steps
   use config,  only: propagator, electronic, reactants, observables
   use sysdef_facade, only: sysdef_load
   use config_roster, only: roster, roster_load
   implicit none
   integer :: nfail = 0, nfail0, i, k, nrows
   character(len=32) :: mode
   character(len=32) :: kname
   logical :: ok_arr, ok_part
   character(len=16) :: blk1(18), blk2(5)
   character(len=16) :: mkey(30), mval(30)

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   ! bad-file naming convention: build/input_case_<mode>.txt (bad-nofile excepted);
   ! the w-surf-/w-mol- prefixed modes drive the full window chain (the paradigm
   ! condition lives at assembly time, so the child continues into roster_load)
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      if (trim(mode) == 'bad-nofile') then
         call read_input('build/input_no_such_file.txt')
      else if (index(trim(mode), 'w-surf-') == 1) then
         call read_input('build/input_case_'//trim(mode)//'.txt')
         call sysdef_load('build/window_surf')
         call roster_load()
      else if (index(trim(mode), 'w-mol-') == 1) then
         call read_input('build/input_case_'//trim(mode)//'.txt')
         call sysdef_load('build/window_mol')
         call roster_load()
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
      'N_SURF = 0', &
      'PES_MAIN = 1', &
      'GWRITE_LEVEL = 5', &
      'OBSERVABLES = t_kin,e_tot', &
      'E_KIN = T', &
      'E_POT = false', &
      'E_TOT = 1', &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 12.5', &
      'B_MAX = 3.0', &
      'R_SEP = 9.75', &
      'DIST_SCHEME = ebk,', &
      '   boltzmann' ])
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
      'n_surf=0', &
      'pes_main=1', &
      'gwrite_level=5', &
      'observables=t_kin,e_tot', &
      'e_kin=T', &
      'e_pot=false', &
      'e_tot=1' ])
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
   ! (j) the 11 dispatch-table member keys: minimal six + every member key
   call write_sample('build/input_case_mkeys.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/sysdef_case_mol', &
      'E_REL = 5.0', &
      'B_MAX = 3.0', &
      'R_SEP = 8.0', &
      'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'T_VIB_A = 350.0', &
      'T_VIB_B = 300.0', &
      'N_E_REL = 2', &
      'T_TRANS = 0.0', &
      'V_WIDTH = 0.3', &
      'N_THTA = 0', &
      'THTA_MAX = 1.1', &
      'N_CHI = 0', &
      'CHI = 0.9', &
      'N_B = 0', &
      'N_AIM = 0', &
      'T_ROT = 300.0', &
      'J_ROT = 2', &
      'N_OSC = 1', &
      'T_OSC = 180.0', &
      'N_LEVEL = 3', &
      'K_OSC = 42.0', &
      'E_VIB = 12.5', &
      'E_STAB = 8.0', &
      'N_E_BAR = 0', &
      'E_BAR = 5.0', &
      'T_BAR = 900.0', &
      'T_BATH = 475.0', &
      'NU_COLL = 0.137', &
      'T_GLO = 550.0', &
      'GAMMA_GLO = 0.09', &
      'W_GHOST = 0.0021', &
      'T_EQ = 325.0', &
      'N_EQ = 250', &
      'DT_EQ = 0.0413' ])
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
      'DIST_SCHEME = thermostat, glo_target, thermalize' ])
   ! (o) window-chain folders (dedicated names - the sysdef_case_* folders belong
   ! to check_sysdef's samples): a molecular window (single .xyz, no cell) and a
   ! surface window (a slab POSCAR carrying the cell) - the paradigm is DERIVED
   ! from the cell's presence, so the two folders are the two paradigms
   call make_dir('build/window_mol')
   call write_sample('build/window_mol/probe.xyz', [ character(len=64) :: &
      '1', 'window molecular sample', 'Ar 0.00 0.00 0.00' ])
   call make_dir('build/window_surf')
   call write_sample('build/window_surf/slab.poscar', [ character(len=64) :: &
      'window surface sample', '1.0', &
      '2.9 0.0 0.0', '0.0 2.9 0.0', '0.0 0.0 18.6', &
      'Au', '1', 'Direct', '0.00 0.00 0.10' ])
   ! (o1) surface window: the six required keys + both window keys
   call write_sample('build/input_case_surfwin.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', &
      'A_LAT = 3.25', 'BOXPAD = 5.5' ])
   ! (o2)/(o3) surface window with one window key removed each time
   call write_sample('build/input_case_w-surf-noalat.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', 'BOXPAD = 5.5' ])
   call write_sample('build/input_case_w-surf-nopad.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_surf', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', 'A_LAT = 3.25' ])
   ! (o4)/(o5) molecular window with one window key present each time
   call write_sample('build/input_case_w-mol-alat.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_mol', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', 'A_LAT = 3.25' ])
   call write_sample('build/input_case_w-mol-bpad.txt', [ character(len=64) :: &
      'SYSTEM_DIR = build/window_mol', &
      'E_REL = 1.0', 'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.01', &
      'DIST_SCHEME = boltzmann', 'BOXPAD = 5.5' ])
   ! (o6)/(o7) non-positive window values (value-domain loud-stop at read_input)
   call write_case('v-alat', 'A_LAT = -3.0')
   call write_case('v-bpad', 'BOXPAD = 0.0')
   ! (p) slot-width families: the offending key at line 1 with a value w
   ! characters wide, its minimal-six row replaced (a duplicate required key
   ! would loud-stop first and mask the width stop)
   call write_width_case('w-title', 'TITLE',       81)
   call write_width_case('w-meth', 'ELEC_METHOD',  33)
   call write_width_case('w-sysd', 'SYSTEM_DIR',  129)
   call write_width_case('w-obs',  'OBSERVABLES', 257)

   ! ---------- bad sample files (one per child loud-stop case) ----------
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
   call write_case('v-ekin',      'E_KIN = maybe')
   call write_case('p-noeq',      'N_TRAJ 7')
   ! DT / DIST_SCHEME bad values: the offending line REPLACES the minimal copy (those
   ! keys are part of the minimal six - a duplicate would loud-stop first)
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

   ! ---------- (a) full key sweep ----------
   call read_input('build/input_case_full.txt')
   call assert_full('(a) full key sweep')

   ! ---------- (i) E21 two-block partition in the parse staging ----------
   ! (run on the (a) staging: 23 user rows, no injected defaults among them)
   nfail0 = nfail
   blk1 = [ character(len=16) :: 'TASK', 'N_TRAJ', 'TITLE', 'ISEED', 'MAX_STEPS', &
            'INTEGRATOR', 'DT', 'ORDER', 'DT_MIN', 'DT_MAX', 'ELEC_METHOD', 'N_SURF', &
            'PES_MAIN', 'GWRITE_LEVEL', 'OBSERVABLES', 'E_KIN', 'E_POT', 'E_TOT' ]
   blk2 = [ character(len=16) :: 'SYSTEM_DIR', 'E_REL', 'B_MAX', 'R_SEP', 'DIST_SCHEME' ]
   nrows = staging_n_rows()
   ok_part = (nrows == 23)
   do i = 1, nrows
      kname = staging_key(i)
      if (any(blk1 == trim(kname))) then
         if (staging_block(i) /= 1) ok_part = .false.
      else if (any(blk2 == trim(kname))) then
         if (staging_block(i) /= 2) ok_part = .false.
      else
         ok_part = .false.
      end if
   end do
   ! line numbers preserved into staging (SYSTEM_DIR sits on source line 19)
   do i = 1, nrows
      if (trim(staging_key(i)) == 'SYSTEM_DIR') exit
   end do
   if (i > nrows) ok_part = .false.
   if (ok_part .and. staging_line(i) == 19) then
      print '(a)', '(i) E21 two-block staging: PASS (18 block-one + 5 block-two rows tagged, line numbers kept)'
   else
      print '(a)', 'FAIL (i) E21 two-block staging: partition or line numbers wrong'
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
       trim(electronic%method) /= 'adiabatic' .or. electronic%n_surf /= 0 .or. &
       trim(electronic%pes_main) /= '1' .or. observables%rec_level /= 2 .or. &
       len_trim(observables%obs_list) /= 0 .or. observables%e_kin .or. &
       observables%e_pot .or. observables%e_tot .or. .not. ok_arr .or. &
       staging_n_rows() /= 23) then
      print '(a)', 'FAIL (c) default injection: a default value is wrong (or staging row count off 23)'
      nfail = nfail + 1
   else
      print '(a)', '(c) default injection: PASS (17 Ruling-D defaults injected, none preset-wide, staging 6+17 rows)'
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
   call expect_loudstop('miss-b_max',        'B_MAX',        '', '(e3) missing B_MAX loud-stop')
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
   call expect_loudstop('v-ekin',  'E_KIN',        'line 1', '(g6) bad logical E_KIN')
   call expect_loudstop('v-dt',    'DT',           'line 1', '(g7) bad real DT')
   call expect_loudstop('v-dtmin', 'DT_MIN',       'line 1', '(g8) bad real DT_MIN')
   call expect_loudstop('v-dtmax', 'DT_MAX',       'line 1', '(g9) bad real DT_MAX')

   ! ---------- (h) parse errors (children) ----------
   call expect_loudstop('p-noeq',     'line 1',  '', '(h1) no-= line loud-stop (line number reported)')
   call expect_loudstop('bad-nofile', 'not found', '', '(h2) missing input file loud-stop')

   ! ---------- (j) member-key admission (dispatch rows stage, verbatim) ----------
   nfail0 = nfail
   call read_input('build/input_case_mkeys.txt')
   mkey = [ character(len=16) :: 'T_VIB_A', 'T_VIB_B', 'N_E_REL', 'T_TRANS', 'V_WIDTH', &
            'N_THTA', 'THTA_MAX', 'N_CHI', 'CHI', 'N_B', 'N_AIM', 'T_ROT', 'J_ROT', &
            'N_OSC', 'T_OSC', 'N_LEVEL', 'K_OSC', 'E_VIB', 'E_STAB', 'N_E_BAR', &
            'E_BAR', 'T_BAR', 'T_BATH', 'NU_COLL', 'T_GLO', 'GAMMA_GLO', 'W_GHOST', &
            'T_EQ', 'N_EQ', 'DT_EQ' ]
   mval = [ character(len=16) :: '350.0', '300.0', '2', '0.0', '0.3', &
            '0', '1.1', '0', '0.9', '0', '0', '300.0', '2', '1', '180.0', &
            '3', '42.0', '12.5', '8.0', '0', '5.0', '900.0', '475.0', '0.137', &
            '550.0', '0.09', '0.0021', '325.0', '250', '0.0413' ]
   ok_part = (staging_n_rows() == 53)      ! 36 user rows (six registry + 30 member
                                            ! keys) + the 17 injected per-key defaults
   do i = 1, size(mkey)
      k = staged_row(trim(mkey(i)))
      if (k == 0 .or. staging_block(k) /= 2 .or. trim(staging_val(k)) /= trim(mval(i))) then
         ok_part = .false.
         exit
      end if
   end do
   k = staged_row('T_VIB_A')               ! line numbers preserved (line 7 of the file)
   if (k == 0 .or. staging_line(k) /= 7) ok_part = .false.
   if (ok_part) then
      print '(a)', '(j) member-key admission: PASS (30 dispatch keys staged block two, values verbatim)'
   else
      print '(a)', 'FAIL (j) member-key admission: a member key missing, wrong block, or value altered'
      nfail = nfail + 1
   end if

   ! ---------- (k)/(l) member-key loud-stops (children) ----------
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

   ! ---------- (n) the three bath-member scheme words parse to codes 8..10 ----------
   nfail0 = nfail
   call read_input('build/input_case_bathdist.txt')
   ok_arr = allocated(reactants%dist_scheme)
   if (ok_arr) ok_arr = size(reactants%dist_scheme) == 3 .and. &
                        all(reactants%dist_scheme == [ 8, 9, 10 ])
   if (ok_arr) then
      print '(a)', '(n) bath-member scheme words: PASS (thermostat/glo_target/'// &
         'thermalize -> codes 8,9,10)'
   else
      print '(a)', 'FAIL (n) bath-member scheme words: DIST_SCHEME list did not map to codes 8..10'
      nfail = nfail + 1
   end if

   ! ---------- (o1) surface window chain: input -> sysdef -> roster ----------
   nfail0 = nfail
   call read_input('build/input_case_surfwin.txt')
   call sysdef_load('build/window_surf')
   call roster_load()
   k = staged_row('A_LAT')
   if (roster%a_lat /= 3.25d0 .or. roster%box_pad /= 5.5d0 .or. &
       reactants%surface_model /= 1 .or. k == 0 .or. staging_block(k) /= 2) then
      print '(a)', 'FAIL (o1) surface window chain: the window scalars did not land in the roster'
      nfail = nfail + 1
   else
      print '(a)', '(o1) surface window chain: PASS (A_LAT=3.25/BOXPAD=5.5 pulled by '// &
         'roster_load from the staging; surface_model=1, rows tagged block two)'
   end if

   ! ---------- (o2)-(o7) window loud-stops (children) ----------
   call expect_loudstop_at('roster_load', 'w-surf-noalat', 'A_LAT', 'missing', &
        '(o2) surface window missing A_LAT loud-stop (required under the surface paradigm)')
   call expect_loudstop_at('roster_load', 'w-surf-nopad', 'BOXPAD', 'missing', &
        '(o3) surface window missing BOXPAD loud-stop')
   call expect_loudstop_at('roster_load', 'w-mol-alat', 'A_LAT', 'molecular paradigm', &
        '(o4) molecular paradigm + A_LAT loud-stop (the window keys are surface-only)')
   call expect_loudstop_at('roster_load', 'w-mol-bpad', 'BOXPAD', 'molecular paradigm', &
        '(o5) molecular paradigm + BOXPAD loud-stop')
   call expect_loudstop('v-alat', 'A_LAT', 'line 1', '(o6) non-positive A_LAT loud-stop (value domain)')
   call expect_loudstop('v-bpad', 'BOXPAD', 'line 1', '(o7) non-positive BOXPAD loud-stop (value domain)')

   ! ---------- (p) slot-width families (children) ----------
   call expect_loudstop('w-title', 'TITLE',       'width 81',  '(p1) TITLE over its 80-wide slot loud-stop')
   call expect_loudstop('w-meth', 'ELEC_METHOD',  'width 33',  '(p2) ELEC_METHOD over its 32-wide slot loud-stop')
   call expect_loudstop('w-sysd', 'SYSTEM_DIR',   'width 129', '(p3) SYSTEM_DIR over its 128-wide slot loud-stop')
   call expect_loudstop('w-obs',  'OBSERVABLES',  'width 257', '(p4) OBSERVABLES over its 256-wide slot loud-stop')

   print '(a,i0,a)', '=== check_input: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! staged_row(key) - staging row index of a key (0 = not staged; a scan over the
   ! public accessors, mirroring what the assembly wiring does to find member keys)
   function staged_row(key) result(idx)
      character(len=*), intent(in) :: key
      integer :: idx, j
      idx = 0
      do j = 1, staging_n_rows()
         if (trim(staging_key(j)) == trim(key)) then
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
          trim(electronic%method) /= 'adiabatic' .or. electronic%n_surf /= 0 .or. &
          trim(electronic%pes_main) /= '1' .or. observables%rec_level /= 5 .or. &
          trim(observables%obs_list) /= 't_kin,e_tot' .or. .not. observables%e_kin .or. &
          observables%e_pot .or. .not. observables%e_tot .or. &
          trim(reactants%system_dir) /= 'build/sysdef_case_mol' .or. &
          reactants%e_rel /= 12.5d0 .or. reactants%b_max /= 3.0d0 .or. &
          reactants%r_sep /= 9.75d0 .or. .not. ok_arr) then
         print '(a)', 'FAIL '//trim(label)//': slot mismatch'
         nfail = nfail + 1
      else
         print '(a)', trim(label)//': PASS (23 keys -> exact slots; lowercase words, ISEED word-1-only)'
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

   ! spawn the child for one bad case and assert the loud-stop: nonzero exit,
   ! stderr names the authority (read_input) and carries the expected key /
   ! extra needle (an empty needle is skipped)
   subroutine expect_loudstop(cmode, key, needle, label)
      character(len=*), intent(in) :: cmode, key, needle, label
      call expect_loudstop_at('read_input', cmode, key, needle, label)
   end subroutine expect_loudstop

   ! the authority-parameterized form: the window criteria stop inside
   ! roster_load (the paradigm condition is assembly-side), so the asserted
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
   ! replaced (a duplicate required key would loud-stop first and mask the stop
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
