!=====================================================================
! input.f90 - L0 input: KEYWORD=VALUE grammar (line-by-line parsing -> per-key default
!             injection -> enumeration registry mapping -> E21 two-block slot mapping ->
!             validation; realized in C4.4 - the STANDALONE grammar: no old-integer-key
!             fallback (E21), no MODEL/preset track (E22))
! Reference: venus2026/input_parser.f90:67-152 (load_input: blank and '#','!' comment
!       skipping :107-108, '=' key/value splitting :125-130, trailing ',' continuation
!       absorption :134-145, key uppercasing :170); :489-502 (set_default - a default is
!       injected only when the user has not set the key; the ONLY survivor of the preset
!       framework, E22); venus2026/param_mapping.f90:34-160 (the mapping-table family -
!       the table idea survives as the enum registries below: task / integrator /
!       dist-scheme word -> internal code, illegal values loud-stop with the name
!       printed; the old map_surface_model role withdrawn - surface_model is DERIVED by
!       sysdef_load, E14 spec §6); venus2026/presets.f90:33-100 (apply_preset - WITHDRAWN
!       by E22: the MODEL key and the five-preset framework are not ported; a stray MODEL
!       is an unknown-key loud-stop); venus2026/venus_input.f90:187-269 (prototype of the
!       read-in main flow), :607-681 (the scattered TRV_A/NACTA-branch member parameters -
!       superseded by the central dispatch-table structure below);
! Statutes realized (C4 plan S3):
!   E21 (two blocks, no fallback): key_tab is organized block one = run control (how
!       THIS run goes) / block two = system definition (WHAT system and initial
!       conditions); every staging row carries its block tag; the old sequential /
!       integer-key grammar has NO fallback track - a data line without '=' is a parse
!       error and old key names are unknown keys.
!   E22 (no MODEL/presets): defaults are PER-KEY (set_default semantics - never
!       overriding user-set keys); the keys formerly backstopped by presets are required
!       or per-key-defaulted per the controller Ruling D assignment (required = no
!       system-independent default exists: SYSTEM_DIR/E_REL/B_MAX/R_SEP/DT/DIST_SCHEME).
!   E18/E19: TERMINATE= / OUTCOME= do not exist (termination and classification are
!       container-side judgments) - both arrive as unknown keys.
!   E14 generalization: container-defined material (composition / masses / geometry /
!       unit cell / surface paradigm - sysdef material) is NOT part of the grammar.
!   E12 reverse stance: the input layer holds ZERO system knowledge - the required/
!       default split follows system-independence only; the sole value-domain check is
!       the formal GWRITE_LEVEL recording-level set {0,1,2,5}.
!   i_seed: rng_init eats a SCALAR decimal seed (src/utils/rng.f90) - the 8-word
!       base-256 conversion lives rng-side, so the scalar is stored directly into
!       control%i_seed(1) (words 2-8 stay zero; the field's (8) shape is the rng-side
!       digit form).
!   Member-parameter dispatch: the central table below is framework-only THIS round
!       (empty - rows arrive with the members in C2(interface)/C5); an empty table
!       means such keys (TRV_A ...) are intercepted as UNKNOWN keys at parse time -
!       no silent swallowing.
!=====================================================================
module input
   use control, only: task, n_traj, title, i_seed, max_steps
   use config,  only: propagator, electronic, reactants, observables
   implicit none
   private
   public :: read_input
   public :: staging_n_rows, staging_key, staging_block, staging_line

   integer, parameter :: key_len  = 32     ! keyword capacity [character]
   integer, parameter :: val_len  = 512    ! value capacity [character]
   integer, parameter :: line_len = 1024   ! physical line buffer [character]
   integer, parameter :: msg_len  = 512    ! error-message capacity [character]

   ! E21 blocks
   integer, parameter :: blk_flow = 1      ! block one: run control
   integer, parameter :: blk_system = 2    ! block two: system definition

   ! value kinds
   integer, parameter :: kind_int = 1      ! integer scalar
   integer, parameter :: kind_real = 2     ! real scalar
   integer, parameter :: kind_char = 3     ! character (stored trimmed, raw)
   integer, parameter :: kind_log = 4      ! logical (T/F/TRUE/FALSE/.TRUE./.FALSE./1/0)
   integer, parameter :: kind_level = 5    ! integer restricted to {0,1,2,5}
   integer, parameter :: kind_enum = 6     ! registry word -> internal code (scalar)
   integer, parameter :: kind_elist = 7    ! comma list of registry words -> code array

   ! enumeration families (which registry a kind_enum/kind_elist key maps through)
   integer, parameter :: enum_none = 0, enum_task = 1, enum_integ = 2, enum_dist = 3

   ! Parse staging (KEYWORD=VALUE row table - module-private; every read_input starts
   ! from a clean slate; read out through the public staging accessors only)
   type :: kv_row_t
      character(len=key_len) :: key = ''   ! keyword (uppercased)
      character(len=val_len) :: val = ''   ! value section (continuations already merged)
      integer :: ln  = 0                   ! source line number (0 = default-injected row)
      integer :: blk = 0                   ! E21 block tag (from the key table)
   end type kv_row_t
   type(kv_row_t), allocatable :: kv_tbl(:)  ! row table (module-private)
   integer :: n_kv = 0                       ! number of rows [row]

   ! Key registry: the first-version key set, block one first (the whole grammar in
   ! one table - a new key is one row, the framework does not change). def = the
   ! Ruling-D per-key default (injected when the user has not set the key);
   ! required = no system-independent default exists (the row's def stays unused)
   type :: key_spec_t
      character(len=key_len) :: name = ''
      integer :: blk = 0
      integer :: kind = 0
      integer :: enum_id = enum_none
      logical :: required = .false.
      character(len=val_len) :: def = ''
   end type key_spec_t

   type :: enum_row_t
      character(len=key_len) :: word = ''
      integer :: code = 0
   end type enum_row_t

   type(key_spec_t), parameter :: key_tab(21) = [ &
      key_spec_t('TASK',         blk_flow,   kind_enum,  enum_task,  .false., 'TRAJECTORY'), &
      key_spec_t('N_TRAJ',       blk_flow,   kind_int,   enum_none,  .false., '1'), &
      key_spec_t('TITLE',        blk_flow,   kind_char,  enum_none,  .false., ''), &
      key_spec_t('ISEED',        blk_flow,   kind_int,   enum_none,  .false., '1'), &
      key_spec_t('MAX_STEPS',    blk_flow,   kind_int,   enum_none,  .false., '100000'), &
      key_spec_t('INTEGRATOR',   blk_flow,   kind_enum,  enum_integ, .false., 'VERLET'), &
      key_spec_t('DT',           blk_flow,   kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('ORDER',        blk_flow,   kind_int,   enum_none,  .false., '1'), &
      key_spec_t('ELEC_METHOD',  blk_flow,   kind_char,  enum_none,  .false., 'ADIABATIC'), &
      key_spec_t('N_SURF',       blk_flow,   kind_int,   enum_none,  .false., '0'), &
      key_spec_t('PES_MAIN',     blk_flow,   kind_char,  enum_none,  .false., '1'), &
      key_spec_t('GWRITE_LEVEL', blk_flow,   kind_level, enum_none,  .false., '1'), &
      key_spec_t('OBSERVABLES',  blk_flow,   kind_char,  enum_none,  .false., ''), &
      key_spec_t('E_KIN',        blk_flow,   kind_log,   enum_none,  .false., 'F'), &
      key_spec_t('E_POT',        blk_flow,   kind_log,   enum_none,  .false., 'F'), &
      key_spec_t('E_TOT',        blk_flow,   kind_log,   enum_none,  .false., 'F'), &
      key_spec_t('SYSTEM_DIR',   blk_system, kind_char,  enum_none,  .true.,  ''), &
      key_spec_t('E_REL',        blk_system, kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('B_MAX',        blk_system, kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('R_SEP',        blk_system, kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('DIST_SCHEME',  blk_system, kind_elist, enum_dist,  .true.,  '') ]

   ! Enumeration registries (param_mapping precedent, table-driven). Codes:
   !   task  - 1 = TRAJECTORY (the one task this program runs; other task words arrive
   !            with their own drivers if ever needed - an unlisted word loud-stops)
   !   integ - the propagator closed-family codes (prop_reg ids, dynamics/propagator.f90)
   !   dist  - the dist_* member families (boltzmann/ebk/incident, methods/); the
   !            integer-code -> member-name mapping is wired at C5
   type(enum_row_t), parameter :: task_tbl(1) = [ &
      enum_row_t('TRAJECTORY', 1) ]
   type(enum_row_t), parameter :: integ_tbl(3) = [ &
      enum_row_t('VERLET', 1), enum_row_t('SYMPLECTIC', 2), enum_row_t('RADAU', 3) ]
   type(enum_row_t), parameter :: dist_tbl(3) = [ &
      enum_row_t('BOLTZMANN', 1), enum_row_t('EBK', 2), enum_row_t('INCIDENT', 3) ]

   ! GWRITE_LEVEL recording-level set (the four-level semantics, N3.1 - the formal
   ! domain of observables%rec_level; recorder opens sec_qp/sec_frc at >=2/>=5)
   integer, parameter :: level_set(4) = [ 0, 1, 2, 5 ]

   ! Member-parameter central dispatch table (one row per member parameter key -> the
   ! owning member's init; row-fill rounds C2(interface)/C5). EMPTY this round: such
   ! keys never pass the unknown-key gate, so nothing is silently swallowed
   type :: dispatch_row_t
      character(len=key_len) :: key = ''    ! member parameter key (e.g. 'TRV_A')
      character(len=32) :: member = ''      ! owning member name (init target)
   end type dispatch_row_t
   type(dispatch_row_t), parameter :: dispatch_tbl(0) = [ dispatch_row_t :: ]
contains
   !------------------------------------------------------------------
   ! read_input(file) - read and map the input (parse -> required check -> per-key
   !                   default injection -> enumeration mapping -> two-block slot
   !                   mapping -> member-parameter dispatch -> validation)
   ! Reference: venus_input.f90:187-269 (prototype of the read-in main flow)
   subroutine read_input(file)
      character(len=*), intent(in) :: file  ! input file name (the old venus_input.f90:189
                                            ! hardwired name, now a dummy argument)
      logical :: ok
      character(len=msg_len) :: errmsg
      call collect_input(file, ok, errmsg)
      if (.not. ok) then
         write (0, '(a)') 'read_input: input error: '//trim(errmsg)
         write (0, '(a)') 'read_input: loud-stop (unusable input - the standalone grammar has no fallback track, E21/E22)'
         stop 1
      end if
   end subroutine read_input

   !------------------------------------------------------------------
   ! staging accessors (public, read-only; the E21 block partition and the source
   !                   line numbers are inspectable for the check - all getters
   !                   return copies, staging itself stays intact)
   !------------------------------------------------------------------
   pure integer function staging_n_rows()
      staging_n_rows = n_kv
   end function staging_n_rows

   function staging_key(i) result(k)
      integer, intent(in) :: i
      character(len=key_len) :: k
      k = ''
      if (i >= 1 .and. i <= n_kv) k = kv_tbl(i)%key
   end function staging_key

   pure integer function staging_block(i)
      integer, intent(in) :: i
      staging_block = 0
      if (i >= 1 .and. i <= n_kv) staging_block = kv_tbl(i)%blk
   end function staging_block

   pure integer function staging_line(i)
      integer, intent(in) :: i
      staging_line = 0
      if (i >= 1 .and. i <= n_kv) staging_line = kv_tbl(i)%ln
   end function staging_line

   !------------------------------------------------------------------
   ! collect_input(file, ok, errmsg) - the whole flow on the soft ok/errmsg channel
   !              (detection layer; read_input converts a bad return into the loud-stop):
   !              reset staging -> parse -> required check -> default injection ->
   !              slot mapping (block one first, block two after - the E21 order)
   subroutine collect_input(file, ok, errmsg)
      character(len=*), intent(in) :: file
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      ok = .false.
      errmsg = ''
      call reset_staging()
      call parse_file(file, ok, errmsg)
      if (.not. ok) return
      call check_required(ok, errmsg)
      if (.not. ok) return
      call inject_defaults()
      call map_slots(ok, errmsg)
      if (.not. ok) return
      ok = .true.
   end subroutine collect_input

   !------------------------------------------------------------------
   ! parse_file(file, ok, errmsg) - line-by-line parse into the staging (old
   !              input_parser semantics): skip blank and '#','!' whole-line comments;
   !              split on '='; a trailing ',' in the value absorbs continuation lines;
   !              uppercase keys; line numbers feed the error channel. A data line
   !              without '=' is a PARSE ERROR (E21: no sequential-grammar fallback);
   !              unknown keys (old names, MODEL, member keys ahead of their members)
   !              and duplicate keys loud-stop with key + line number
   ! Reference: input_parser.f90:97-148 (the loop this modernizes)
   subroutine parse_file(file, ok, errmsg)
      character(len=*), intent(in) :: file
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=line_len) :: line, kw, vals
      integer :: u, ios, eqp, ln, kln, ri, dupe, blk
      logical :: ex
      ok = .false.
      errmsg = ''
      inquire (file=trim(file), exist=ex)
      if (.not. ex) then
         errmsg = 'input file not found: '//trim(file)
         return
      end if
      open (newunit=u, file=trim(file), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         errmsg = 'cannot open input file: '//trim(file)
         return
      end if
      ln = 0
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         ln = ln + 1
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle
         eqp = index(line, '=')
         if (eqp == 0) then
            errmsg = 'line '//trim(i2s(ln))//': not a KEYWORD=VALUE line: '//trim(line)// &
                     ' (the standalone grammar has no sequential fallback - E21)'
            close (u)
            return
         end if
         kw = to_upper(adjustl(line(:eqp-1)))
         if (len_trim(kw) == 0) then
            errmsg = 'line '//trim(i2s(ln))//': empty keyword before ='
            close (u)
            return
         end if
         kln = ln                              ! the key's OWN line (continuations do not move it)
         vals = adjustl(line(eqp+1:))
         ! trailing ',' absorbs continuation lines (old :134-145)
         do while (len_trim(vals) > 0 .and. vals(len_trim(vals):len_trim(vals)) == ',')
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) then
               errmsg = 'line '//trim(i2s(ln))//': continuation runs past the end of the file (key '// &
                        trim(kw)//')'
               close (u)
               return
            end if
            ln = ln + 1
            vals = trim(vals)//' '//adjustl(line)
         end do
         if (len_trim(vals) > val_len) then
            errmsg = 'key '//trim(kw)//' (line '//trim(i2s(kln))//'): value longer than '// &
                     trim(i2s(val_len))//' characters'
            close (u)
            return
         end if
         ! registry membership: key table or member-parameter dispatch table
         ri = key_index(trim(kw))
         if (ri == 0 .and. dispatch_index(trim(kw)) == 0) then
            errmsg = 'unknown key '//trim(kw)//' (line '//trim(i2s(kln))//'): not in the key table '// &
                     '(old key names, withdrawn keys and MODEL are intercepted as unknown - '// &
                     'E21/E22; member-parameter keys arrive with their members, C2/C5)'
            close (u)
            return
         end if
         dupe = staged_index(trim(kw))
         if (dupe > 0) then
            errmsg = 'duplicate key '//trim(kw)//' (lines '//trim(i2s(kv_tbl(dupe)%ln))//' and '// &
                     trim(i2s(kln))//') - one occurrence per key'
            close (u)
            return
         end if
         if (ri > 0) then
            blk = key_tab(ri)%blk
         else
            blk = blk_system   ! member-parameter keys sit in block two (the E21 table's last row)
         end if
         ! capacity guard (the known tripwire of C4.4's deferred list): kv_tbl is sized
         ! to the key table - the moment dispatch rows go real (C2/C5), member keys
         ! reach here too and must NOT write past the table (same defect class as the
         ! phase's other OOB bugs; loud-stop, never a silent overflow)
         if (n_kv >= size(kv_tbl)) then
            errmsg = 'staging capacity exhausted at key '//trim(kw)//' (line '//trim(i2s(kln))// &
                     '): distinct staged keys exceed size(kv_tbl) (sized to the key table - '// &
                     'member-parameter keys need their dispatch rows counted, C2/C5)'
            close (u)
            return
         end if
         n_kv = n_kv + 1
         kv_tbl(n_kv)%key = trim(kw)
         kv_tbl(n_kv)%val = trim(vals)
         kv_tbl(n_kv)%ln = kln
         kv_tbl(n_kv)%blk = blk
      end do
      close (u)
      ok = .true.
   end subroutine parse_file

   !------------------------------------------------------------------
   ! check_required(ok, errmsg) - every required key staged with a non-empty value
   !              (Ruling D: a required key is one with NO system-independent default;
   !              an empty value counts as missing)
   subroutine check_required(ok, errmsg)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      integer :: r, k
      ok = .false.
      errmsg = ''
      do r = 1, size(key_tab)
         if (.not. key_tab(r)%required) cycle
         k = staged_index(trim(key_tab(r)%name))
         if (k == 0 .or. len_trim(kv_tbl(k)%val) == 0) then
            errmsg = 'missing required key '//trim(key_tab(r)%name)//' - no system-independent '// &
                     'default exists for it (per-key-default statute E22, Ruling D)'
            return
         end if
      end do
      ok = .true.
   end subroutine check_required

   !------------------------------------------------------------------
   ! inject_defaults() - set_default semantics (input_parser.f90:489-502): each
   !              non-required key absent from the user's rows gets its per-key default
   !              staged (ln = 0 marks a default-injected row); user-set keys are
   !              NEVER overridden
   subroutine inject_defaults()
      integer :: r
      do r = 1, size(key_tab)
         if (key_tab(r)%required) cycle
         if (staged_index(trim(key_tab(r)%name)) > 0) cycle
         ! capacity guard (defensive, the phase's fail-safe stance): unreachable while
         ! the user table holds fewer rows than the registry (a full user table has no
         ! absent keys to inject) - guarded anyway so no future fill route writes past
         ! the table (same defect class as the parse_file tripwire above)
         if (n_kv >= size(kv_tbl)) then
            write (0, '(a)') 'inject_defaults: staging capacity exhausted at key '// &
                             trim(key_tab(r)%name)//': rows would exceed size(kv_tbl)'
            write (0, '(a)') 'inject_defaults: loud-stop (defensive capacity guard - '// &
                             'unreachable through the current fill routes)'
            stop 1
         end if
         n_kv = n_kv + 1
         kv_tbl(n_kv)%key = key_tab(r)%name
         kv_tbl(n_kv)%val = key_tab(r)%def
         kv_tbl(n_kv)%ln = 0
         kv_tbl(n_kv)%blk = key_tab(r)%blk
      end do
   end subroutine inject_defaults

   !------------------------------------------------------------------
   ! map_slots(ok, errmsg) - the E21 two-block slot mapping, in key-table order
   !              (block one first, block two after): for every key the staged value is
   !              parsed per its kind (bad values loud-stop with key + line) and written
   !              into its control/config slot; the member-parameter dispatch loop
   !              follows (empty table this round)
   subroutine map_slots(ok, errmsg)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      type(key_spec_t) :: sp
      character(len=val_len) :: sval
      integer :: r, k, ln, ival, d
      integer, allocatable :: ivec(:)
      real(8) :: rval
      logical :: lval
      ok = .false.
      errmsg = ''
      do r = 1, size(key_tab)
         sp = key_tab(r)
         k = staged_index(trim(sp%name))
         ! guaranteed staged: required keys were checked, per-key defaults injected (the
         ! guard is fail-safe against a future key-table row flipping the flags wrongly)
         if (k == 0) then
            errmsg = 'internal: key '//trim(sp%name)//' neither user-set nor defaulted'
            return
         end if
         sval = kv_tbl(k)%val
         ln = kv_tbl(k)%ln
         select case (sp%kind)
         case (kind_int)
            call parse_int(sp, sval, ln, ival, ok, errmsg)
         case (kind_real)
            call parse_real(sp, sval, ln, rval, ok, errmsg)
         case (kind_char)
            ok = .true.                     ! raw trimmed string, no parsing
         case (kind_log)
            call parse_log(sp, sval, ln, lval, ok, errmsg)
         case (kind_level)
            call parse_level(sp, sval, ln, ival, ok, errmsg)
         case (kind_enum)
            call parse_enum(sp, sval, ln, ival, ok, errmsg)
         case (kind_elist)
            call parse_elist(sp, sval, ln, ivec, ok, errmsg)
         case default
            errmsg = 'internal: unhandled key kind for '//trim(sp%name)
            return
         end select
         if (.not. ok) return
         select case (trim(sp%name))
         case ('TASK')                                  ! block one: run control
            task = ival
         case ('N_TRAJ')
            n_traj = ival
         case ('TITLE')
            title(1) = sval                              ! one TITLE key -> line 1 (old TITLE1/2 withdrawn)
            title(2) = ''
         case ('ISEED')
            i_seed = 0                                   ! scalar statute: word 1 only, the
            i_seed(1) = ival                             ! 8-word conversion is rng-side
         case ('MAX_STEPS')
            max_steps = ival
         case ('INTEGRATOR')
            propagator%id = ival
         case ('DT')
            propagator%dt = rval
         case ('ORDER')
            propagator%order = ival
         case ('ELEC_METHOD')
            electronic%method = to_lower(sval)           ! package-name string (E13/E14:
                                                         ! lookup at the facade, not here)
         case ('N_SURF')
            electronic%n_surf = ival
         case ('PES_MAIN')
            electronic%pes_main = sval                   ! name string, stored as written
         case ('GWRITE_LEVEL')
            observables%rec_level = ival
         case ('OBSERVABLES')
            observables%obs_list = sval                  ! raw comma list (C4.5 rec_sel_lst input)
         case ('E_KIN')
            observables%e_kin = lval
         case ('E_POT')
            observables%e_pot = lval
         case ('E_TOT')
            observables%e_tot = lval
         case ('SYSTEM_DIR')                             ! block two: system definition
            reactants%system_dir = sval
         case ('E_REL')
            reactants%e_rel = rval
         case ('B_MAX')
            reactants%b_max = rval
         case ('R_SEP')
            reactants%r_sep = rval
         case ('DIST_SCHEME')
            if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
            allocate (reactants%dist_scheme(size(ivec)))
            reactants%dist_scheme = ivec                 ! per-fragment code array (E15; the
                                                         ! size-vs-n_frag reconcile is C5 wiring)
         case default
            errmsg = 'internal: unhandled key slot for '//trim(sp%name)
            return
         end select
      end do
      ! member-parameter dispatch (central table -> member init; EMPTY this round - see
      ! the dispatch_tbl declaration; the loop site is the fixed framework)
      do d = 1, size(dispatch_tbl)
         ! row-fill round (C2(interface)/C5): stage the key's value into dispatch_tbl(d)
      end do
      ok = .true.
   end subroutine map_slots

   !------------------------------------------------------------------
   ! value parsers (soft ok/errmsg channel; every error names the key and line)
   !------------------------------------------------------------------
   subroutine parse_int(sp, sval, ln, ival, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      integer, intent(out) :: ival
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      integer :: ios
      if (index(sval, ',') > 0) then
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': list value not allowed here ("'// &
                  trim(sval)//'")'
         return
      end if
      read (sval, *, iostat=ios) ival
      if (ios /= 0) then
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': bad integer value "'//trim(sval)//'"'
         return
      end if
      ok = .true.
   end subroutine parse_int

   subroutine parse_real(sp, sval, ln, rval, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      real(8), intent(out) :: rval
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      integer :: ios
      if (index(sval, ',') > 0) then
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': list value not allowed here ("'// &
                  trim(sval)//'")'
         return
      end if
      read (sval, *, iostat=ios) rval
      if (ios /= 0) then
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': bad real value "'//trim(sval)//'"'
         return
      end if
      ok = .true.
   end subroutine parse_real

   subroutine parse_log(sp, sval, ln, lval, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      logical, intent(out) :: lval
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      select case (to_upper(trim(sval)))
      case ('T', 'TRUE', '.TRUE.', '1')
         lval = .true.
      case ('F', 'FALSE', '.FALSE.', '0')
         lval = .false.
      case default
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': bad logical value "'//trim(sval)// &
                  '" (T/F/TRUE/FALSE/.TRUE./.FALSE./1/0)'
         return
      end select
      ok = .true.
   end subroutine parse_log

   subroutine parse_level(sp, sval, ln, ival, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      integer, intent(out) :: ival
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      call parse_int(sp, sval, ln, ival, ok, errmsg)
      if (.not. ok) return
      if (any(ival == level_set)) then
         ok = .true.
      else
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': value '//trim(i2s(ival))// &
                  ' outside the recording-level set 0/1/2/5 (GWRITE_LEVEL semantics)'
      end if
   end subroutine parse_level

   subroutine parse_enum(sp, sval, ln, ival, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      integer, intent(out) :: ival
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      if (index(sval, ',') > 0) then
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': list value not allowed here ("'// &
                  trim(sval)//'")'
         return
      end if
      ival = enum_lookup(sp%enum_id, to_upper(trim(sval)))
      if (ival == 0) then
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': illegal value "'//trim(sval)// &
                  '" (legal: '//trim(legal_words(sp%enum_id))//')'
         return
      end if
      ok = .true.
   end subroutine parse_enum

   subroutine parse_elist(sp, sval, ln, ivec, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      integer, allocatable, intent(out) :: ivec(:)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=len_trim(sval)) :: rest
      character(len=key_len) :: tok
      integer :: n, i, ic, code
      ok = .false.
      errmsg = ''
      n = 1 + count([(sval(i:i) == ',', i = 1, len_trim(sval))])
      allocate (ivec(n))
      n = 0
      rest = adjustl(sval)
      do while (len_trim(rest) > 0)
         ic = index(trim(rest), ',')
         if (ic == 0) then
            tok = adjustl(trim(rest))
            rest = ''
         else
            tok = adjustl(rest(:ic-1))
            rest = adjustl(rest(ic+1:))
         end if
         if (len_trim(tok) == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': empty item in the value list'
            return
         end if
         code = enum_lookup(sp%enum_id, to_upper(trim(tok)))
         if (code == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': illegal value "'//trim(tok)// &
                     '" (legal: '//trim(legal_words(sp%enum_id))//')'
            return
         end if
         n = n + 1
         ivec(n) = code
      end do
      if (n == 0) then
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': empty value list'
         return
      end if
      ivec = ivec(1:n)
      ok = .true.
   end subroutine parse_elist

   !------------------------------------------------------------------
   ! private helpers (pure string/table mechanics - no system knowledge)
   !------------------------------------------------------------------
   subroutine reset_staging()
      if (allocated(kv_tbl)) deallocate (kv_tbl)
      allocate (kv_tbl(size(key_tab)))   ! one row per registry key at most (unknown
      n_kv = 0                           ! keys rejected, duplicates rejected); BOTH
                                         ! fill sites (parse_file / inject_defaults)
                                         ! guard their increment with a capacity
                                         ! loud-stop - capacity is code, not memory
   end subroutine reset_staging

   ! key_index(key) - registry row of a key (0 = not a key-table key)
   pure integer function key_index(key)
      character(len=*), intent(in) :: key
      integer :: i
      key_index = 0
      do i = 1, size(key_tab)
         if (trim(key_tab(i)%name) == trim(key)) then
            key_index = i
            exit
         end if
      end do
   end function key_index

   ! dispatch_index(key) - dispatch-table row of a member-parameter key (0 = none; the
   ! table is empty this round - the site where member rows will be looked up)
   pure integer function dispatch_index(key)
      character(len=*), intent(in) :: key
      integer :: i
      dispatch_index = 0
      do i = 1, size(dispatch_tbl)
         if (trim(dispatch_tbl(i)%key) == trim(key)) then
            dispatch_index = i
            exit
         end if
      end do
   end function dispatch_index

   ! staged_index(key) - staging row of a key (0 = not staged)
   pure integer function staged_index(key)
      character(len=*), intent(in) :: key
      integer :: i
      staged_index = 0
      do i = 1, n_kv
         if (trim(kv_tbl(i)%key) == trim(key)) then
            staged_index = i
            exit
         end if
      end do
   end function staged_index

   ! enum_lookup(enum_id, word) - registry code of an (uppercased) word (0 = not found)
   pure integer function enum_lookup(enum_id, word)
      integer, intent(in) :: enum_id
      character(len=*), intent(in) :: word
      integer :: i
      enum_lookup = 0
      select case (enum_id)
      case (enum_task)
         do i = 1, size(task_tbl)
            if (trim(task_tbl(i)%word) == trim(word)) then
               enum_lookup = task_tbl(i)%code
               exit
            end if
         end do
      case (enum_integ)
         do i = 1, size(integ_tbl)
            if (trim(integ_tbl(i)%word) == trim(word)) then
               enum_lookup = integ_tbl(i)%code
               exit
            end if
         end do
      case (enum_dist)
         do i = 1, size(dist_tbl)
            if (trim(dist_tbl(i)%word) == trim(word)) then
               enum_lookup = dist_tbl(i)%code
               exit
            end if
         end do
      end select
   end function enum_lookup

   ! legal_words(enum_id) - the legal-value list for an error message
   function legal_words(enum_id) result(s)
      integer, intent(in) :: enum_id
      character(len=128) :: s
      integer :: i
      s = ''
      select case (enum_id)
      case (enum_task)
         do i = 1, size(task_tbl)
            if (i > 1) s = trim(s)//', '
            s = trim(s)//trim(task_tbl(i)%word)
         end do
      case (enum_integ)
         do i = 1, size(integ_tbl)
            if (i > 1) s = trim(s)//', '
            s = trim(s)//trim(integ_tbl(i)%word)
         end do
      case (enum_dist)
         do i = 1, size(dist_tbl)
            if (i > 1) s = trim(s)//', '
            s = trim(s)//trim(dist_tbl(i)%word)
         end do
      end select
   end function legal_words

   ! line_tag(ln) - the ' (line n)' error suffix ('' for a default-injected row)
   function line_tag(ln) result(tag)
      integer, intent(in) :: ln
      character(len=24) :: tag
      if (ln > 0) then
         tag = ' (line '//trim(i2s(ln))//')'
      else
         tag = ' (default)'
      end if
   end function line_tag

   ! i2s(n) - decimal string of an integer
   function i2s(n) result(s)
      integer, intent(in) :: n
      character(len=16) :: s
      write (s, '(i0)') n
   end function i2s

   ! to_upper(s) - ASCII uppercase of a trimmed string
   pure function to_upper(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len_trim(s)) :: t
      integer :: i, ch
      t = trim(adjustl(s))
      do i = 1, len_trim(t)
         ch = iachar(t(i:i))
         if (ch >= iachar('a') .and. ch <= iachar('z')) t(i:i) = achar(ch - 32)
      end do
   end function to_upper

   ! to_lower(s) - ASCII lowercase of a trimmed string (the config method-name form)
   pure function to_lower(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len_trim(s)) :: t
      integer :: i, ch
      t = trim(adjustl(s))
      do i = 1, len_trim(t)
         ch = iachar(t(i:i))
         if (ch >= iachar('A') .and. ch <= iachar('Z')) t(i:i) = achar(ch + 32)
      end do
   end function to_lower

end module input
