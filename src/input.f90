!=====================================================================
! input.f90 - KEYWORD=VALUE input grammar (parse -> required check -> per-key default
!             injection -> enumeration mapping -> two-block slot mapping -> validation)
! Design:
!   Standalone grammar, no fallback track: a data line without '=' is a parse error;
!   old key names and MODEL are unknown keys. Two blocks - run control / system
!   definition; member-parameter keys stage raw (central dispatch table) for the
!   member inits at assembly; declared container keys (the dragged-in container's
!   parameter vocabulary, input_declare_keys BEFORE the parse) stage the same way
!   and are consumed at the assembly seam by pull (absent key keeps the target's
!   own default; input_audit_declared fatals a buffered-but-unconsumed declared key);
!   this layer holds zero system knowledge (formal checks
!   only: kind domains, slot widths). Defaults are per-key, never overriding
!   user-set keys; required keys are those with no system-independent default
!   (SYSTEM_DIR, E_REL, B_MAX, R_SEP, DT, DIST_SCHEME). GWRITE_LEVEL defaults to 2;
!   RANDOM_ORIENT defaults to F (the archived five-case stream).
!=====================================================================
module input
   use control, only: task, n_traj, title, i_seed, max_steps
   use config,  only: propagator, electronic, reactants, observables
   implicit none
   private
   public :: read_input
   public :: buffer_n_rows, buffer_key, buffer_block, buffer_line, buffer_val
   public :: input_declare_keys, input_audit_declared, pull

   integer, parameter :: key_len  = 32     ! keyword capacity [character]
   integer, parameter :: val_len  = 512    ! value capacity [character]
   integer, parameter :: line_len = 1024   ! physical line buffer [character]
   integer, parameter :: msg_len  = 512    ! error-message capacity [character]

   ! grammar blocks
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
   integer, parameter :: kind_posreal = 8  ! positive real scalar (parses as real AND > 0 -
                                           ! the surface lattice length)

   ! enumeration families (which registry a kind_enum/kind_elist key maps through)
   integer, parameter :: enum_none = 0, enum_task = 1, enum_integ = 2, enum_dist = 3, enum_bath = 4

   ! Parse buffer (KEYWORD=VALUE row table - module-private, clean slate per
   ! read_input, read out through the public buffer accessors only)
   type :: kv_row_t
      character(len=key_len) :: key = ''   ! keyword (uppercased)
      character(len=val_len) :: val = ''   ! value section (continuations already merged)
      integer :: ln  = 0                   ! source line number (0 = default-injected row)
      integer :: blk = 0                   ! block tag (from the key table)
      logical :: used = .false.            ! a declared-key row consumed by a pull
   end type kv_row_t
   type(kv_row_t), allocatable :: kv_tbl(:)  ! row table (module-private)
   integer :: n_kv = 0                       ! number of rows [row]

   ! Declared container keys (the dragged-in container's parameter vocabulary):
   ! its reg _all hands the key list over BEFORE read_input runs, the parse admits
   ! those names, and the container's reg _init consumes each buffered row by pull.
   ! Process-lifetime state (not reset by read_input - declaration precedes parse).
   integer, parameter :: max_declared = 128          ! declared-key capacity [key]
   character(len=key_len) :: declared_tbl(max_declared) = ''
   integer :: n_declared = 0                          ! declared-key count [key]

   ! pull(key, target [, found]) - the declared-key consumption entry at the
   ! assembly seam: one call, one key, type-dispatched on the target
   interface pull
      module procedure pull_real, pull_int, pull_char, pull_r3
   end interface pull

   ! Key registry: the whole grammar in one table, block one first (a new key is one
   ! row, the framework does not change). def = per-key default (injected when the
   ! user has not set the key); required = no system-independent default exists (the
   ! row's def stays unused); slotw = the value-slot width - an over-wide value
   ! aborts with a named error at slot mapping instead of silently truncating (0 = no fixed slot);
   ! cond = paradigm-conditional key (the surface window): never default-injected,
   ! never unconditionally required - the paradigm condition runs at assembly time
   ! (the value stays buffered for the assembly-time pull)
   type :: key_spec_t
      character(len=key_len) :: name = ''
      integer :: blk = 0
      integer :: kind = 0
      integer :: enum_id = enum_none
      logical :: required = .false.
      character(len=val_len) :: def = ''
      integer :: slotw = 0
      logical :: cond = .false.
   end type key_spec_t

   type :: enum_row_t
      character(len=key_len) :: word = ''
      integer :: code = 0
   end type enum_row_t

   ! DT_MIN/DT_MAX (the adaptive member's sequence-length bounds): kind_real, NOT
   ! kind_posreal - a clamp value <= 0 is the legal unclamped default (0 = off);
   ! consumed by the radau member only (fixed-step members ignore them)
   type(key_spec_t), parameter :: key_tab(30) = [ &
      key_spec_t('TASK',         blk_flow,   kind_enum,  enum_task,  .false., 'TRAJECTORY'), &
      key_spec_t('N_TRAJ',       blk_flow,   kind_int,   enum_none,  .false., '1'), &
      key_spec_t('TITLE',        blk_flow,   kind_char,  enum_none,  .false., '', slotw=80), &
      key_spec_t('ISEED',        blk_flow,   kind_int,   enum_none,  .false., '1'), &
      key_spec_t('MAX_STEPS',    blk_flow,   kind_int,   enum_none,  .false., '100000'), &
      key_spec_t('INTEGRATOR',   blk_flow,   kind_enum,  enum_integ, .false., 'VERLET'), &
      key_spec_t('DT',           blk_flow,   kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('ORDER',        blk_flow,   kind_int,   enum_none,  .false., '1'), &
      key_spec_t('DT_MIN',       blk_flow,   kind_real,  enum_none,  .false., '0.0'), &
      key_spec_t('DT_MAX',       blk_flow,   kind_real,  enum_none,  .false., '0.0'), &
      key_spec_t('ELEC_METHOD',  blk_flow,   kind_char,  enum_none,  .false., 'ADIABATIC', slotw=32), &
      key_spec_t('N_SURF',       blk_flow,   kind_int,   enum_none,  .false., '0'), &
      key_spec_t('PES_MAIN',     blk_flow,   kind_char,  enum_none,  .false., '1', slotw=32), &
      key_spec_t('GWRITE_LEVEL', blk_flow,   kind_level, enum_none,  .false., '2'), &
      key_spec_t('OBSERVABLES',  blk_flow,   kind_char,  enum_none,  .false., '', slotw=256), &
      key_spec_t('E_KIN',        blk_flow,   kind_log,   enum_none,  .false., 'F'), &
      key_spec_t('E_POT',        blk_flow,   kind_log,   enum_none,  .false., 'F'), &
      key_spec_t('E_TOT',        blk_flow,   kind_log,   enum_none,  .false., 'F'), &
      key_spec_t('SYSTEM_DIR',   blk_system, kind_char,  enum_none,  .true.,  '', slotw=128), &
      key_spec_t('E_REL',        blk_system, kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('B_MAX',        blk_system, kind_real,  enum_none,  .false., '0.0'), & ! optional (PAIR channel only): the surface member loud-stops a WRITTEN B_MAX - lateral placement is N_AIM's alone (ruled 2026-10-01); the pair member enforces a positive radius at assembly (beam_init)
      key_spec_t('R_SEP',        blk_system, kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('DIST_SCHEME',  blk_system, kind_elist, enum_dist,  .true.,  ''), &
      key_spec_t('RANDOM_ORIENT', blk_system, kind_log,  enum_none,  .false., 'F'), &
      key_spec_t('BATH',          blk_system, kind_enum, enum_bath,  .false., 'NONE'), &
      key_spec_t('N_BATH',        blk_system, kind_int,  enum_none,  .false., '0'), &
      key_spec_t('DT_BATH',       blk_system, kind_real, enum_none,  .false., '0.0'), &
      key_spec_t('PROJECTILE',   blk_system, kind_char,  enum_none,  .false., '', slotw=64, cond=.true.), &
      key_spec_t('TARGET',       blk_system, kind_char,  enum_none,  .false., '', slotw=64, cond=.true.), &
      key_spec_t('A_LAT',        blk_system, kind_posreal, enum_none, .false., '', cond=.true.) ]

   ! Enumeration registries (table-driven word -> code). task: 1 = TRAJECTORY (the one
   ! task this program runs - an unlisted word aborts with a named error); integ: the propagator
   ! closed-family codes (prop_reg ids); dist: the PER-FRAGMENT scheme words the
   ! sampler's scheme-code lookup resolves. 'INCIDENT' is withdrawn from the dist
   ! registry - the incident channel is selected by the samp_run scheme argument, so
   ! as a per-fragment word it is an ILLEGAL value
   type(enum_row_t), parameter :: task_tbl(1) = [ &
      enum_row_t('TRAJECTORY', 1) ]
   type(enum_row_t), parameter :: integ_tbl(3) = [ &
      enum_row_t('VERLET', 1), enum_row_t('SYMPLECTIC', 2), enum_row_t('RADAU', 3) ]
   type(enum_row_t), parameter :: dist_tbl(10) = [ &
      enum_row_t('BOLTZMANN', 1), enum_row_t('EBK', 2), &
      enum_row_t('ROTATION', 3), enum_row_t('J', 4), &
      enum_row_t('SURFACE_OSCILLATOR', 5), enum_row_t('NORMALMODE', 6), &
      enum_row_t('BARRIER_EXCITATION', 7), &
      enum_row_t('GLO_TARGET', 8), &
      enum_row_t('THERMALIZE', 9), enum_row_t('WIGNER', 10) ]

   ! bath provision words (THERMOSTAT left the dist registry 2026-09-22: the
   ! Andersen species is a bath, not a per-fragment distribution word - the
   ! BATH selector arms the pre-evolution bath loop inside samp_run)
   type(enum_row_t), parameter :: bath_tbl(2) = [ &
      enum_row_t('NONE', 1), enum_row_t('ANDERSEN', 2) ]


   ! GWRITE_LEVEL recording-level set - the formal domain of observables%rec_level
   ! (recorder opens sec_qp/sec_frc at >=2/>=5)
   integer, parameter :: level_set(4) = [ 0, 1, 2, 5 ]

   ! Member-parameter central dispatch table (one row per member parameter key -> the
   ! owning member's init). Key names are the semantic uppercase of the member's field
   ! names (T_VIB_A for boltzmann's t_vib_a, ...); a member key outside this table is
   ! an unknown key at parse time (no silent swallowing)
   type :: dispatch_row_t
      character(len=key_len) :: key = ''    ! member parameter key (e.g. 'T_VIB_A')
      character(len=32) :: member = ''      ! owning member name (init target)
   end type dispatch_row_t
   type(dispatch_row_t), parameter :: dispatch_tbl(35) = [ &
      dispatch_row_t('T_VIB_A',  'boltzmann'), &
      dispatch_row_t('T_VIB_B',  'boltzmann'), &
      dispatch_row_t('N_E_REL',  'incident_surface/pair'), &
      dispatch_row_t('T_TRANS',  'incident_surface/pair'), &
      dispatch_row_t('V_WIDTH',  'incident_surface/pair'), &
      dispatch_row_t('N_THTA',   'incident_surface/pair'), &
      dispatch_row_t('THTA_MAX', 'incident_surface/pair'), &
      dispatch_row_t('N_CHI',    'incident_surface/pair'), &
      dispatch_row_t('CHI',      'incident_surface/pair'), &
      dispatch_row_t('N_B',      'incident_surface/pair'), &
      dispatch_row_t('N_AIM',    'incident_surface/pair'), &
      dispatch_row_t('AIM_X',    'incident_surface'), &
      dispatch_row_t('AIM_Y',    'incident_surface'), &
      dispatch_row_t('T_ROT',    'j'), &
      dispatch_row_t('J_ROT',    'rotation'), &
      dispatch_row_t('N_OSC',    'surface_oscillator'), &
      dispatch_row_t('T_OSC',    'surface_oscillator'), &
      dispatch_row_t('N_LEVEL',  'surface_oscillator'), &
      dispatch_row_t('K_OSC',    'surface_oscillator'), &
      dispatch_row_t('E_VIB',    'normalmode'), &
      dispatch_row_t('E_STAB',   'barrier_excitation'), &
      dispatch_row_t('N_E_BAR',  'barrier_excitation'), &
      dispatch_row_t('E_BAR',    'barrier_excitation'), &
      dispatch_row_t('T_BAR',    'barrier_excitation'), &
      dispatch_row_t('T_BATH',   'bath_andersen'), &
      dispatch_row_t('NU_COLL',  'bath_andersen'), &
      dispatch_row_t('T_GLO',    'glo_target'), &
      dispatch_row_t('GAMMA_GLO','glo_target'), &
      dispatch_row_t('W_GHOST',  'glo_target'), &
      dispatch_row_t('T_EQ',     'thermalize'), &
      dispatch_row_t('N_EQ',     'thermalize'), &
      dispatch_row_t('DT_EQ',    'thermalize'), &
      dispatch_row_t('E_SCALE',  'wigner'), &
      dispatch_row_t('N_VIB',    'ebk'), &
      dispatch_row_t('N_ROT',    'ebk') ]
contains
   !------------------------------------------------------------------
   ! read_input(file) - read and map the input; a bad return from the detection
   !                   layer converts into the named abort
   !------------------------------------------------------------------
   subroutine read_input(file)
      character(len=*), intent(in) :: file  ! input file name (a dummy argument)
      logical :: ok
      character(len=msg_len) :: errmsg
      call collect_input(file, ok, errmsg)
      if (.not. ok) then
         write (0, '(a)') 'read_input: input error: '//trim(errmsg)
         write (0, '(a)') 'read_input: fatal (unusable input - the standalone grammar has no fallback track, E21/E22)'
         stop 1
      end if
   end subroutine read_input

   !------------------------------------------------------------------
   ! buffer accessors (public, read-only; all getters return copies, buffer
   !                   itself stays intact)
   !------------------------------------------------------------------
   pure integer function buffer_n_rows()
      buffer_n_rows = n_kv
   end function buffer_n_rows

   function buffer_key(i) result(k)
      integer, intent(in) :: i
      character(len=key_len) :: k
      k = ''
      if (i >= 1 .and. i <= n_kv) k = kv_tbl(i)%key
   end function buffer_key

   pure integer function buffer_block(i)
      integer, intent(in) :: i
      buffer_block = 0
      if (i >= 1 .and. i <= n_kv) buffer_block = kv_tbl(i)%blk
   end function buffer_block

   pure integer function buffer_line(i)
      integer, intent(in) :: i
      buffer_line = 0
      if (i >= 1 .and. i <= n_kv) buffer_line = kv_tbl(i)%ln
   end function buffer_line

   function buffer_val(i) result(v)
      integer, intent(in) :: i
      character(len=val_len) :: v
      v = ''
      if (i >= 1 .and. i <= n_kv) v = kv_tbl(i)%val
   end function buffer_val

   !------------------------------------------------------------------
   ! declared-key channel (container parameter keys): declare before the
   ! parse (the container reg _all), consume at the assembly seam (pull),
   ! audit after the seams (input_audit_declared)
   !------------------------------------------------------------------
   subroutine input_declare_keys(keys)
      character(len=*), intent(in) :: keys(:)  ! the container's key vocabulary
      integer :: i
      do i = 1, size(keys)
         if (len_trim(keys(i)) == 0) then
            write (0, '(a)') 'input_declare_keys: an empty key name in the container list'
            write (0, '(a)') 'input_declare_keys: fatal (a container key must be a name)'
            stop 1
         end if
         if (len_trim(keys(i)) > key_len) then
            write (0, '(a)') 'input_declare_keys: key '//trim(keys(i))//' exceeds the '// &
                             'keyword width '//trim(i2s(key_len))
            write (0, '(a)') 'input_declare_keys: fatal (a container key must fit the grammar)'
            stop 1
         end if
         if (key_index(trim(keys(i))) > 0 .or. dispatch_index(trim(keys(i))) > 0) then
            write (0, '(a)') 'input_declare_keys: key '//trim(keys(i))//' collides with a '// &
                             'program key or a member-parameter key'
            write (0, '(a)') 'input_declare_keys: fatal (a container key must be a new name)'
            stop 1
         end if
         if (declared_index(trim(keys(i))) > 0) then
            write (0, '(a)') 'input_declare_keys: key '//trim(keys(i))//' declared twice'
            write (0, '(a)') 'input_declare_keys: fatal (one declaration per container key)'
            stop 1
         end if
         if (n_declared >= max_declared) then
            write (0, '(a)') 'input_declare_keys: declared-key capacity '//trim(i2s(max_declared))// &
                             ' exhausted at key '//trim(keys(i))
            write (0, '(a)') 'input_declare_keys: fatal (defensive capacity guard)'
            stop 1
         end if
         n_declared = n_declared + 1
         declared_tbl(n_declared) = to_upper(trim(keys(i)))   ! the parse uppercases the
      end do                                                 ! user's rows; compare in kind
   end subroutine input_declare_keys

   ! pull_row(key) - the shared consumption locate: buffer row of a declared key
   !              (0 = never buffered, the caller default stands); an undeclared pull
   !              or a second pull of one key is a container-authoring named abort
   function pull_row(key) result(k)
      character(len=*), intent(in) :: key
      integer :: k
      if (declared_index(trim(key)) == 0) then
         write (0, '(a)') 'pull: key '//trim(key)//' is not a declared container key '// &
                          '(input_declare_keys runs in the container reg _all)'
         write (0, '(a)') 'pull: fatal (an undeclared pull is a container-authoring error)'
         stop 1
      end if
      k = staged_index(trim(key))
      if (k == 0) return                       ! never buffered: the caller default stands
      if (kv_tbl(k)%used) then
         write (0, '(a)') 'pull: key '//trim(key)//' pulled twice (one key, one consumer)'
         write (0, '(a)') 'pull: fatal (a double pull is a container-authoring error)'
         stop 1
      end if
      kv_tbl(k)%used = .true.
   end function pull_row

   ! pull_fatal(key, errmsg) - a failed cast of a buffered value (the parse-side
   !                         parsers already named key and line in errmsg)
   subroutine pull_fatal(errmsg)
      character(len=*), intent(in) :: errmsg
      write (0, '(a)') 'pull: input error: '//trim(errmsg)
      write (0, '(a)') 'pull: fatal (a buffered container parameter is not consumable)'
      stop 1
   end subroutine pull_fatal

   !------------------------------------------------------------------
   ! pull(key, target [, found]) - one declared key, one call: the buffered value
   !              lands in the target (absent key leaves the caller's default;
   !              found reports presence for the optional-slot semantics)
   !------------------------------------------------------------------
   subroutine pull_real(key, target, found)
      character(len=*), intent(in) :: key
      real(8), intent(inout) :: target
      logical, intent(out), optional :: found
      type(key_spec_t) :: sp
      character(len=val_len) :: sval
      character(len=msg_len) :: emsg
      integer :: k, ln
      logical :: ok
      k = pull_row(key)
      if (present(found)) found = k > 0
      if (k == 0) return
      sval = kv_tbl(k)%val
      ln = kv_tbl(k)%ln
      sp%name = trim(key)                     ! the parsers read only the key name
      call parse_real(sp, sval, ln, target, ok, emsg)
      if (.not. ok) call pull_fatal(emsg)
   end subroutine pull_real

   subroutine pull_int(key, target, found)
      character(len=*), intent(in) :: key
      integer, intent(inout) :: target
      logical, intent(out), optional :: found
      type(key_spec_t) :: sp
      character(len=val_len) :: sval
      character(len=msg_len) :: emsg
      integer :: k, ln
      logical :: ok
      k = pull_row(key)
      if (present(found)) found = k > 0
      if (k == 0) return
      sval = kv_tbl(k)%val
      ln = kv_tbl(k)%ln
      sp%name = trim(key)
      call parse_int(sp, sval, ln, target, ok, emsg)
      if (.not. ok) call pull_fatal(emsg)
   end subroutine pull_int

   subroutine pull_char(key, target, found)
      character(len=*), intent(in) :: key
      character(len=*), intent(inout) :: target
      logical, intent(out), optional :: found
      integer :: k
      k = pull_row(key)
      if (present(found)) found = k > 0
      if (k == 0) return
      if (len_trim(kv_tbl(k)%val) > len(target)) then
         write (0, '(a)') 'pull: key '//trim(key)//': value width '// &
                          trim(i2s(len_trim(kv_tbl(k)%val)))//' exceeds its slot width '// &
                          trim(i2s(len(target)))//' (silent truncation closed - fatal)'
         stop 1
      end if
      target = trim(kv_tbl(k)%val)
   end subroutine pull_char

   subroutine pull_r3(key, target, found)
      character(len=*), intent(in) :: key
      real(8), intent(inout) :: target(3)
      logical, intent(out), optional :: found
      integer :: k, ln, ios, i, j
      k = pull_row(key)
      if (present(found)) found = k > 0
      if (k == 0) return
      ln = kv_tbl(k)%ln
      ! exactly three comma-separated items (the strict form: no fourth item,
      ! no empty item - a list-directed read of three reals fails the rest)
      i = count([(kv_tbl(k)%val(j:j) == ',', j = 1, len_trim(kv_tbl(k)%val))])
      if (i /= 2) then
         write (0, '(a)') 'pull: key '//trim(key)//trim(line_tag(ln))// &
                          ': bad 3-vector value "'//trim(kv_tbl(k)%val)//'" (exactly three '// &
                          'comma-separated numbers expected)'
         write (0, '(a)') 'pull: fatal (a buffered container parameter is not consumable)'
         stop 1
      end if
      read (kv_tbl(k)%val, *, iostat=ios) target(1), target(2), target(3)
      if (ios /= 0) then
         write (0, '(a)') 'pull: key '//trim(key)//trim(line_tag(ln))// &
                          ': bad 3-vector value "'//trim(kv_tbl(k)%val)//'"'
         write (0, '(a)') 'pull: fatal (a buffered container parameter is not consumable)'
         stop 1
      end if
   end subroutine pull_r3

   !------------------------------------------------------------------
   ! input_audit_declared() - after the assembly seams: every buffered declared-key
   !              row was pulled exactly once; a row nobody consumed is a named
   !              fatal (a pull typo in the container seam, or the seam skipped -
   !              e.g. container keys buffered while SYSTEM_DIR was not)
   !------------------------------------------------------------------
   subroutine input_audit_declared()
      integer :: k
      logical :: bad
      bad = .false.
      do k = 1, n_kv
         if (declared_index(trim(kv_tbl(k)%key)) == 0) cycle
         if (kv_tbl(k)%used) cycle
         if (.not. bad) then
            bad = .true.
            write (0, '(a)') 'input_audit_declared: buffered container key(s) never pulled:'
         end if
         write (0, '(a)') '  '//trim(kv_tbl(k)%key)//trim(line_tag(kv_tbl(k)%ln))// &
                          ' - the dragged-in container did not consume it (a pull typo, '// &
                          'or the seam skipped: SYSTEM_DIR unstaged?)'
      end do
      if (bad) then
         write (0, '(a)') 'input_audit_declared: fatal (a buffered container key outside '// &
                          'the container seam is a configuration or authoring error)'
         stop 1
      end if
   end subroutine input_audit_declared

   !------------------------------------------------------------------
   ! collect_input(file, ok, errmsg) - the whole flow on the soft ok/errmsg channel
   !              (detection layer; read_input converts a bad return into the
   !              named abort): reset buffer -> parse -> required check -> default
   !              injection -> slot mapping (block one first, block two after)
   subroutine collect_input(file, ok, errmsg)
      character(len=*), intent(in) :: file
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      ok = .false.
      errmsg = ''
      call reset_buffer()
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
   ! parse_file(file, ok, errmsg) - line-by-line parse into the buffer: skip blank
   !              and '#','!' whole-line comments; split on '='; a trailing ',' in
   !              the value absorbs continuation lines; uppercase keys; line numbers
   !              feed the error channel. A data line without '=' is a PARSE ERROR
   !              (no sequential-grammar fallback); unknown keys (old names, MODEL,
   !              member keys outside the dispatch table) and duplicate keys
   !              abort with a named error (key + line) number
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
         ! 1. split on '='; a data line without '=' is a parse error
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
         ! 2. a trailing ',' absorbs continuation lines
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
         ! 3. registry membership: key table, member-parameter dispatch table, or
         !    the declared container-key set
         ri = key_index(trim(kw))
         if (ri == 0 .and. dispatch_index(trim(kw)) == 0 .and. declared_index(trim(kw)) == 0) then
            errmsg = 'unknown key '//trim(kw)//' (line '//trim(i2s(kln))//'): not in the key table '// &
                     '(old key names, withdrawn keys and MODEL are intercepted as unknown - '// &
                     'E21/E22; member keys outside the dispatch table and container keys '// &
                     'outside the declared set likewise)'
            close (u)
            return
         end if
         ! 4. one occurrence per key
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
            blk = blk_system   ! member-parameter keys sit in block two
         end if
         ! 5. capacity guard: kv_tbl is sized to the key table PLUS the dispatch rows
         ! and the declared container keys (all of them stage here) - named abort,
         ! never a silent overflow
         if (n_kv >= size(kv_tbl)) then
            errmsg = 'buffer capacity exhausted at key '//trim(kw)//' (line '//trim(i2s(kln))// &
                     '): distinct buffered keys exceed size(kv_tbl) (sized to the key table '// &
                     'plus the dispatch rows and the declared container keys)'
            close (u)
            return
         end if
         ! 6. stage the row
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
   ! check_required(ok, errmsg) - every required key buffered with a non-empty value
   !              (an empty value counts as missing)
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
   ! inject_defaults() - each non-required key absent from the user's rows gets its
   !              per-key default buffered (ln = 0 marks a default-injected row);
   !              user-set keys are NEVER overridden
   subroutine inject_defaults()
      integer :: r
      do r = 1, size(key_tab)
         if (key_tab(r)%required) cycle
         if (key_tab(r)%cond) cycle     ! conditional keys never inject (a buffered row
                                       ! is always a user row)
         if (staged_index(trim(key_tab(r)%name)) > 0) cycle
         ! defensive capacity guard (unreachable through the current fill routes):
         ! named abort rather than write past the table
         if (n_kv >= size(kv_tbl)) then
            write (0, '(a)') 'inject_defaults: buffer capacity exhausted at key '// &
                             trim(key_tab(r)%name)//': rows would exceed size(kv_tbl)'
            write (0, '(a)') 'inject_defaults: fatal (defensive capacity guard - '// &
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
   ! map_slots(ok, errmsg) - two-block slot mapping, in key-table order (block one
   !              first, block two after): for every key the buffered value is parsed
   !              per its kind (bad values abort with a named error (key + line)) and written
   !              into its control/config slot; the member-parameter dispatch loop
   !              follows (buffered member keys must be non-empty)
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
      ! 1. slot mapping loop: parse each buffered value per its kind and write its slot
      do r = 1, size(key_tab)
         sp = key_tab(r)
         k = staged_index(trim(sp%name))
         ! guaranteed buffered: required keys were checked, defaults injected (the guard
         ! is fail-safe against a future key-table row flipping the flags wrongly)
         if (k == 0) then
            if (sp%cond) cycle    ! conditional key absent: nothing to map - the
                                  ! paradigm condition is assembly-side
            ok = .false.
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
         case (kind_posreal)
            call parse_posreal(sp, sval, ln, rval, ok, errmsg)
         case default
            ok = .false.                          ! stale-ok leak guard
            errmsg = 'internal: unhandled key kind for '//trim(sp%name)
            return
         end select
         if (.not. ok) return
         ! slot-width guard: an over-wide value aborts with a named error (key + line + actual
         ! width) instead of silently truncating at the slot write below
         if (sp%slotw > 0 .and. len_trim(sval) > sp%slotw) then
            ok = .false.
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': value width '// &
                     trim(i2s(len_trim(sval)))//' exceeds its slot width '//trim(i2s(sp%slotw))// &
                     ' (silent truncation closed - fatal)'
            return
         end if
         select case (trim(sp%name))
         case ('TASK')                                  ! block one: run control
            task = ival
         case ('N_TRAJ')
            n_traj = ival
         case ('TITLE')
            title(1) = sval                              ! one TITLE key -> line 1
            title(2) = ''
         case ('ISEED')
            i_seed = 0                                   ! scalar rule: word 1 only, the
            i_seed(1) = ival                             ! 8-word conversion is rng-side
         case ('MAX_STEPS')
            max_steps = ival
         case ('INTEGRATOR')
            propagator%id = ival
         case ('DT')
            propagator%dt = rval
         case ('ORDER')
            propagator%order = ival
         case ('DT_MIN')
            propagator%dt_min = rval
         case ('DT_MAX')
            propagator%dt_max = rval
         case ('ELEC_METHOD')
            electronic%method = to_lower(sval)           ! package-name string (lookup at
                                                         ! the interface, not here)
         case ('N_SURF')
            electronic%n_surf = ival
         case ('PES_MAIN')
            electronic%pes_main = sval                   ! name string, stored as written
         case ('GWRITE_LEVEL')
            observables%rec_level = ival
         case ('OBSERVABLES')
            observables%obs_list = sval                  ! raw comma list (rec_sel_lst input)
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
            reactants%dist_scheme = ivec                 ! per-fragment code array (the
                                                         ! size-vs-n_frag reconcile is assembly wiring)
         case ('RANDOM_ORIENT')
            reactants%random_orient = lval                ! universal orientation stage switch
                                                         ! (the sampler stage reads it; default F
                                                         ! keeps the archived stream)
         case ('BATH')
            reactants%bath = ival                         ! bath provision selector (arms the
                                                         ! pre-evolution bath loop; NONE default)
         case ('N_BATH')
            reactants%n_bath = ival                       ! bath-loop step count
         case ('DT_BATH')
            reactants%dt_bath = rval                      ! bath-loop step size [10 fs]

         case ('A_LAT', 'PROJECTILE', 'TARGET') ! surface-window / role keys: no slot
            continue                                     ! write here - the values stay buffered and
                                                         ! list_atoms_load pulls them at assembly time
                                                         ! (paradigm condition included)
         case default
            ok = .false.
            errmsg = 'internal: unhandled key slot for '//trim(sp%name)
            return
         end select
      end do
      ! 2. member-parameter dispatch (central table): a buffered member key must carry a
      ! non-empty value (key + line named on failure - the value semantics themselves
      ! are member-side knowledge; the assembly wiring hands the buffered values to the
      ! member inits)
      do d = 1, size(dispatch_tbl)
         k = staged_index(trim(dispatch_tbl(d)%key))
         if (k > 0 .and. len_trim(kv_tbl(k)%val) == 0) then
            ok = .false.
            errmsg = 'key '//trim(dispatch_tbl(d)%key)//trim(line_tag(kv_tbl(k)%ln))// &
                     ': empty member-parameter value (the owning member "'// &
                     trim(dispatch_tbl(d)%member)//'" would receive nothing)'
            return
         end if
      end do
      ! 3. declared container keys: same non-empty rule (the value semantics are
      !    container-side knowledge; the seam pull hands the buffered value over)
      do d = 1, n_declared
         k = staged_index(trim(declared_tbl(d)))
         if (k > 0 .and. len_trim(kv_tbl(k)%val) == 0) then
            ok = .false.
            errmsg = 'key '//trim(declared_tbl(d))//trim(line_tag(kv_tbl(k)%ln))// &
                     ': empty container-parameter value (the container pull '// &
                     'would receive nothing)'
            return
         end if
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

   subroutine parse_posreal(sp, sval, ln, rval, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      real(8), intent(out) :: rval
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      call parse_real(sp, sval, ln, rval, ok, errmsg)
      if (.not. ok) return
      if (rval <= 0.0d0) then
         ok = .false.
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': value must be positive ("'// &
                  trim(sval)//'" - the surface-window keys carry lengths)'
         return
      end if
      ok = .true.
   end subroutine parse_posreal

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
   subroutine reset_buffer()
      if (allocated(kv_tbl)) deallocate (kv_tbl)
      allocate (kv_tbl(size(key_tab) + size(dispatch_tbl) + n_declared))   ! one row per
      n_kv = 0                           ! registry key PLUS one per dispatch row (member
                                         ! keys stage too) PLUS one per declared container
                                         ! key; both fill sites guard their increment
                                         ! with a capacity named abort
   end subroutine reset_buffer

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

   ! dispatch_index(key) - dispatch-table row of a member-parameter key (0 = none)
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

   ! declared_index(key) - declared-set row of a container parameter key (0 = none;
   !                       stored uppercased, compared in kind)
   pure integer function declared_index(key)
      character(len=*), intent(in) :: key
      integer :: i
      declared_index = 0
      do i = 1, n_declared
         if (trim(declared_tbl(i)) == trim(key)) then
            declared_index = i
            exit
         end if
      end do
   end function declared_index

   ! staged_index(key) - buffer row of a key (0 = not buffered)
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
      case (enum_bath)
         do i = 1, size(bath_tbl)
            if (trim(bath_tbl(i)%word) == trim(word)) then
               enum_lookup = bath_tbl(i)%code
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
      case (enum_bath)
         do i = 1, size(bath_tbl)
            if (i > 1) s = trim(s)//', '
            s = trim(s)//trim(bath_tbl(i)%word)
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
