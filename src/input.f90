!=====================================================================
! input.f90 - KEYWORD=VALUE input grammar (parse -> required check -> per-key default
!             injection -> enumeration mapping -> two-block slot mapping -> validation)
! Design:
!   Standalone grammar, no fallback track: a data line without '=' is a parse error;
!   old key names and MODEL are unknown keys. Two blocks - run control / system
!   definition; CHANNEL-LEVEL member-parameter keys (the incident beam laws and
!   the bath species) stage raw through the central dispatch table for their member
!   inits at assembly, while a scheme word's own parameters ride ONLY its
!   DIST_SCHEME bracket (WORD[NAME=VALUE, ...] - the 2026-10-13 bracket ruling:
!   no top-level member key for them, no fallback, a bare word keeps the member's
!   field defaults); declared container keys (the dragged-in container's parameter
!   vocabulary, input_declare_keys BEFORE the parse) stage the same way and are
!   consumed at the assembly seam by pull (absent key keeps the target's own
!   default; input_audit_declared fatals a buffered-but-unconsumed declared key,
!   input_audit_dist_par the bracket analogue); this layer holds zero system
!   knowledge (formal checks only: kind domains, slot widths, bracket FORM - which
!   parameter names a word carries is member-side). Defaults are per-key, never
!   overriding user-set keys; required keys are those with no system-independent
!   default (SYSTEM_DIR, E_REL, B_MAX, R_SEP, DT, DIST_SCHEME). GWRITE_LEVEL
!   defaults to 2; RANDOM_ORIENT defaults to F (the archived five-case stream).
!   SPECTRUM_SOURCE defaults to COMPUTE (the internal derivation; MANUAL stages
!   the spectrum-table file - 2026-10-05 spectrum-source legislation).
!=====================================================================
module input
   use control, only: task, n_traj, title, i_seed, max_steps
   use config,  only: propagator, electronic, reactants, observables
   implicit none
   private
   public :: read_input
   public :: buffer_n_rows, buffer_key, buffer_block, buffer_line, buffer_val
   public :: input_declare_keys, input_audit_declared, pull
   public :: par_int, par_real, input_audit_dist_par

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
   integer, parameter :: kind_cell = 8     ! 3x3 unit cell: three comma-separated rows of
                                           ! three reals (the surface CELL key)

   ! enumeration families (which registry a kind_enum/kind_elist key maps through)
   integer, parameter :: enum_none = 0, enum_task = 1, enum_integ = 2, enum_dist = 3, &
                          enum_bath = 4, enum_spectrum = 5, enum_scattering = 6

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
   type(key_spec_t), parameter :: key_tab(28) = [ &
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
      key_spec_t('GWRITE_LEVEL', blk_flow,   kind_level, enum_none,  .false., '2'), &
      key_spec_t('OBSERVABLES',  blk_flow,   kind_char,  enum_none,  .false., '', slotw=256), &
      key_spec_t('SYSTEM_DIR',   blk_system, kind_char,  enum_none,  .true.,  '', slotw=128), &
      key_spec_t('E_REL',        blk_system, kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('B_MAX',        blk_system, kind_real,  enum_none,  .false., '0.0'), & ! optional (PAIR channel only): the surface member loud-stops a WRITTEN B_MAX - lateral placement is N_AIM's alone (ruled 2026-10-01); the pair member enforces a positive radius at assembly (beam_init)
      key_spec_t('R_SEP',        blk_system, kind_real,  enum_none,  .true.,  ''), &
      key_spec_t('DIST_SCHEME',  blk_system, kind_elist, enum_dist,  .true.,  ''), &
      key_spec_t('RANDOM_ORIENT', blk_system, kind_log,  enum_none,  .false., 'F'), &
      key_spec_t('BATH',          blk_system, kind_enum, enum_bath,  .false., 'NONE'), &
      key_spec_t('N_BATH',        blk_system, kind_int,  enum_none,  .false., '0'), &
      key_spec_t('DT_BATH',       blk_system, kind_real, enum_none,  .false., '0.0'), &
      key_spec_t('SCATTERING',    blk_system, kind_enum, enum_scattering, .false., 'GAS'), &
      key_spec_t('CELL',          blk_system, kind_cell, enum_none,  .false., '', cond=.true.), &
      key_spec_t('SPECTRUM_SOURCE', blk_system, kind_enum, enum_spectrum, .false., 'COMPUTE'), &
      key_spec_t('SPECTRUM_FILE',  blk_system, kind_char, enum_none,  .false., '', slotw=128), &
      key_spec_t('PROJECTILE',   blk_system, kind_char,  enum_none,  .false., '', slotw=64, cond=.true.), &
      key_spec_t('TARGET',       blk_system, kind_char,  enum_none,  .false., '', slotw=64, cond=.true.) ]

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

   ! spectrum producer words (2026-10-05 spectrum-source legislation: the container
   ! export slot is withdrawn - the producer is the internal derivation (default) or
   ! the manual spectrum-table file; the consumer seam is interface/spectrum_interface)
   type(enum_row_t), parameter :: spectrum_tbl(2) = [ &
      enum_row_t('COMPUTE', 1), enum_row_t('MANUAL', 2) ]

   ! scattering-class words (2026-10-13 ruling: the paradigm is a main-program input
   ! decision, no longer derived from the system files; GAS keeps the archived stream)
   type(enum_row_t), parameter :: scattering_tbl(2) = [ &
      enum_row_t('GAS', 1), enum_row_t('SURFACE', 2) ]


   ! GWRITE_LEVEL recording-level set - the formal domain of observables%rec_level
   ! (recorder opens sec_qp/sec_frc at >=2/>=5)
   integer, parameter :: level_set(4) = [ 0, 1, 2, 5 ]

   ! Member-parameter central dispatch table (one row per CHANNEL-LEVEL member
   ! parameter key -> the owning member's init). The incident beam laws and the
   ! bath species keys stage here (they belong to the one-per-process channel,
   ! not to a fragment); a scheme word's own parameters ride its DIST_SCHEME
   ! bracket instead (the 2026-10-13 ruling retired their 21 dispatch rows - a
   ! top-level line of one of those names is an unknown key). Key names are the
   ! semantic uppercase of the member's field names; a member key outside this
   ! table is an unknown key at parse time (no silent swallowing)
   type :: dispatch_row_t
      character(len=key_len) :: key = ''    ! member parameter key (e.g. 'N_B')
      character(len=32) :: member = ''      ! owning member name (init target)
   end type dispatch_row_t
   type(dispatch_row_t), parameter :: dispatch_tbl(13) = [ &
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
      dispatch_row_t('T_BATH',   'bath_andersen'), &
      dispatch_row_t('NU_COLL',  'bath_andersen') ]

   ! Per-fragment bracket parameters (the 2026-10-13 ruling's delivery form):
   ! each DIST_SCHEME item may carry WORD[NAME=VALUE, ...]. The parse splits the
   ! FORM (names/values staged uppercased, duplicate names rejected) and owns no
   ! member semantics - which parameter names a word carries and their domains
   ! are the owning member's init business; the member consumes its carriers'
   ! pairs at the assembly seam through par_int/par_real, and
   ! input_audit_dist_par fatals a pair nobody consumed (the container-key
   ! audit's bracket analogue)
   type :: par_pair_t
      character(len=key_len) :: name = ''   ! parameter name (uppercased)
      character(len=val_len) :: val = ''    ! raw value text (uppercased)
      logical :: used = .false.             ! a member pull consumed the pair
   end type par_pair_t
   type :: frag_par_t
      character(len=key_len) :: word = ''        ! the item's scheme word (uppercased)
      type(par_pair_t), allocatable :: pair(:)   ! its bracket pairs (unallocated = bare word)
      integer :: n_pair = 0                      ! staged pair count [-]
   end type frag_par_t
   type(frag_par_t), allocatable :: dist_pars(:) ! one row per DIST_SCHEME item
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
   ! par_int / par_real(word, ifrag, name, target [, found]) - one carrier's
   !              bracket-parameter consumption at the member-seam: the typed
   !              value lands in the target; an absent name keeps the caller's
   !              default (found = F - no fallback anywhere); the owning
   !              member's init validates domains. One pull marks the pair
   !              consumed (input_audit_dist_par after the seams)
   !------------------------------------------------------------------
   subroutine par_int(word, ifrag, name, target, found)
      character(len=*), intent(in) :: word, name
      integer, intent(in) :: ifrag
      integer, intent(inout) :: target
      logical, intent(out), optional :: found
      integer :: j, ios
      call par_row(word, ifrag, name, j)
      if (present(found)) found = j > 0
      if (j == 0) return
      read (dist_pars(ifrag)%pair(j)%val, *, iostat=ios) target
      if (ios /= 0) call par_fatal(word, ifrag, name, 'bad integer value "'// &
         trim(dist_pars(ifrag)%pair(j)%val)//'"')
   end subroutine par_int

   subroutine par_real(word, ifrag, name, target, found)
      character(len=*), intent(in) :: word, name
      integer, intent(in) :: ifrag
      real(8), intent(inout) :: target
      logical, intent(out), optional :: found
      integer :: j, ios
      call par_row(word, ifrag, name, j)
      if (present(found)) found = j > 0
      if (j == 0) return
      read (dist_pars(ifrag)%pair(j)%val, *, iostat=ios) target
      if (ios /= 0) call par_fatal(word, ifrag, name, 'bad real value "'// &
         trim(dist_pars(ifrag)%pair(j)%val)//'"')
   end subroutine par_real

   ! par_row(word, ifrag, name, j) - the shared consumption locate: pair row of
   !              a fragment's parameter (0 = not staged, the caller default
   !              stands); ownership and double-pull are member-authoring
   !              named aborts (one name, one consumer)
   subroutine par_row(word, ifrag, name, j)
      character(len=*), intent(in) :: word, name
      integer, intent(in) :: ifrag
      integer, intent(out) :: j
      integer :: k
      if (.not. allocated(dist_pars)) then
         write (0, '(a)') 'par: no DIST_SCHEME pair table staged (read_input must run first)'
         write (0, '(a)') 'par: fatal (the bracket seam was consulted outside the parse)'
         stop 1
      end if
      if (ifrag < 1 .or. ifrag > size(dist_pars)) then
         write (0, '(a)') 'par: fragment index '//trim(i2s(ifrag))//' lies outside the '// &
                          'DIST_SCHEME list of '//trim(i2s(size(dist_pars)))//' item(s)'
         write (0, '(a)') 'par: fatal (a member pulling a fragment it does not own)'
         stop 1
      end if
      if (trim(dist_pars(ifrag)%word) /= to_upper(trim(word))) then
         write (0, '(a)') 'par: fragment '//trim(i2s(ifrag))//' carries scheme word '// &
                          trim(dist_pars(ifrag)%word)//', not '//trim(to_upper(word))
         write (0, '(a)') 'par: fatal (a member pulling a fragment it does not own)'
         stop 1
      end if
      do j = 1, dist_pars(ifrag)%n_pair
         if (trim(dist_pars(ifrag)%pair(j)%name) == to_upper(trim(name))) then
            if (dist_pars(ifrag)%pair(j)%used) then
               write (0, '(a)') 'par: parameter '//trim(name)//' of fragment '// &
                                trim(i2s(ifrag))//' pulled twice (one name, one consumer)'
               write (0, '(a)') 'par: fatal (a double pull is a member-authoring error)'
               stop 1
            end if
            dist_pars(ifrag)%pair(j)%used = .true.
            return
         end if
      end do
      j = 0
   end subroutine par_row

   ! par_fatal(word, ifrag, name, why) - a failed cast of a staged pair
   subroutine par_fatal(word, ifrag, name, why)
      character(len=*), intent(in) :: word, name, why
      integer, intent(in) :: ifrag
      write (0, '(a)') 'par: DIST_SCHEME fragment '//trim(i2s(ifrag))//' ('//trim(to_upper(word))// &
                       '): parameter '//trim(name)//': '//trim(why)
      write (0, '(a)') 'par: fatal (a bracket parameter is not consumable)'
      stop 1
   end subroutine par_fatal

   !------------------------------------------------------------------
   ! input_audit_dist_par() - after the member seams: every staged bracket
   !              pair was pulled exactly once; a pair nobody consumed is a
   !              named fatal (a parameter name the owning member does not
   !              know, or the member's seam skipped)
   !------------------------------------------------------------------
   subroutine input_audit_dist_par()
      integer :: i, j
      logical :: bad
      if (.not. allocated(dist_pars)) return   ! no DIST_SCHEME staged: clean
      bad = .false.
      do i = 1, size(dist_pars)
         do j = 1, dist_pars(i)%n_pair
            if (dist_pars(i)%pair(j)%used) cycle
            if (.not. bad) then
               bad = .true.
               write (0, '(a)') 'input_audit_dist_par: unconsumed bracket parameter(s):'
            end if
            write (0, '(a)') '  '//trim(dist_pars(i)%word)//' fragment '//trim(i2s(i))//': '// &
                             trim(dist_pars(i)%pair(j)%name)//' - no member pull consumed it '// &
                             '(a parameter the owning member does not know?)'
         end do
      end do
      if (bad) then
         write (0, '(a)') 'input_audit_dist_par: fatal (a bracket parameter outside its '// &
                          'member seam is a configuration or authoring error)'
         stop 1
      end if
   end subroutine input_audit_dist_par

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
      real(8) :: rval, cvec(3,3)
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
         case (kind_cell)
            call parse_cell(sp, sval, ln, cvec, ok, errmsg)
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
         case ('GWRITE_LEVEL')
            observables%rec_level = ival
         case ('OBSERVABLES')
            observables%obs_list = sval                  ! raw comma list (rec_sel_lst input)
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
            reactants%dt_bath = rval                      ! bath-loop step size [fs]

         case ('SPECTRUM_SOURCE')                     ! spectrum producer selector (consumed
            reactants%spectrum_src = ival             ! at the member seam, spectrum_interface)
         case ('SPECTRUM_FILE')
            reactants%spectrum_file = sval            ! spectrum-table path (MANUAL's records;
                                                     ! a staged path under COMPUTE never reaches
                                                     ! a pull - the cross rule below stops it)

         case ('SCATTERING')                       ! paradigm master key (2026-10-13 ruling):
            reactants%surface_model = ival - 1      ! GAS -> 0, SURFACE -> 1 (the input owns
                                                     ! the decision; the system files carry
                                                     ! geometry only)
         case ('CELL')                             ! unit cell by rows [A] (surface-only; the
            reactants%cell = cvec                  ! window reconcile requires it under
            reactants%cell_present = .true.        ! SURFACE and rejects it under GAS; the
                                                     ! box derivation consumes it at list_atoms)

         case ('PROJECTILE', 'TARGET') ! role keys: no slot
            continue                                 ! write here - the values stay buffered and
                                                     ! list_atoms_load pulls them at assembly time
                                                     ! (paradigm condition included)
         case default
            ok = .false.
            errmsg = 'internal: unhandled key slot for '//trim(sp%name)
            return
         end select
      end do
      ! 1b. formal cross-key rule (formal, no system knowledge): a staged
      !     SPECTRUM_FILE under COMPUTE would sit unused (no pull happens) - a
      !     silent unused key is not this grammar's way
      if (len_trim(reactants%spectrum_file) > 0 .and. reactants%spectrum_src /= 2) then
         ok = .false.
         errmsg = 'SPECTRUM_FILE is staged while SPECTRUM_SOURCE is COMPUTE - a '// &
                  'spectrum-table file requires SPECTRUM_SOURCE=MANUAL'
         return
      end if
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

   subroutine parse_cell(sp, sval, ln, cvec, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      real(8), intent(out) :: cvec(3,3)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=len_trim(sval)) :: rest, grp
      integer :: i, ic, ios
      ok = .false.
      ! exactly three comma-separated groups (the three lattice-vector rows)
      if (count([(sval(i:i) == ',', i = 1, len_trim(sval))]) /= 2) then
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': bad CELL value "'//trim(sval)// &
                  '" (exactly three comma-separated vector rows expected)'
         return
      end if
      rest = adjustl(sval)
      do i = 1, 3
         ic = index(trim(rest), ',')
         if (ic == 0) then
            grp = adjustl(trim(rest))
            rest = ''
         else
            grp = adjustl(rest(:ic-1))
            rest = adjustl(rest(ic+1:))
         end if
         read (grp, *, iostat=ios) cvec(i,1), cvec(i,2), cvec(i,3)
         if (ios /= 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': bad CELL row '//trim(i2s(i))// &
                     ' "'//trim(grp)//'" (three numbers per row expected)'
            return
         end if
      end do
      ok = .true.
   end subroutine parse_cell

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

   !------------------------------------------------------------------
   ! parse_elist(sp, ...) - comma list of registry words -> code array. The
   !              split runs at TOP-LEVEL commas only (depth 0 outside
   !              brackets): under enum_dist each item may carry
   !              WORD[NAME=VALUE, ...] - the word maps to its code and the
   !              bracket pairs stage per item in dist_pars (form only); a
   !              bracket under any other enum family is an illegal value
   !              (brackets are the dist family's form)
   subroutine parse_elist(sp, sval, ln, ivec, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: sval
      integer, intent(in) :: ln
      integer, allocatable, intent(out) :: ivec(:)
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=len_trim(sval)) :: rest, item
      integer :: n, i, ic, depth, code, nw
      type(frag_par_t), allocatable :: pars(:)
      ok = .false.
      errmsg = ''
      ! 1. whole-value bracket balance BEFORE splitting (a value-wide net-open
      !    bracket is a form defect no per-item check would place correctly)
      depth = 0
      do i = 1, len_trim(sval)
         select case (sval(i:i))
         case ('['); depth = depth + 1
         case (']'); depth = depth - 1
         end select
      end do
      if (depth /= 0) then
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': unbalanced bracket in "'// &
                  trim(sval)//'" (every bracket opens and closes inside its list item)'
         return
      end if
      ! 2. split at top-level commas; one code row and one pair row per item
      n = 1 + count_top_commas(sval)
      allocate (ivec(n), pars(n))
      nw = 0
      rest = adjustl(sval)
      do while (len_trim(rest) > 0)
         ic = top_comma(trim(rest))
         if (ic == 0) then
            item = adjustl(trim(rest))
            rest = ''
         else
            item = adjustl(rest(:ic-1))
            rest = adjustl(rest(ic+1:))
         end if
         if (len_trim(item) == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': empty item in the value list'
            return
         end if
         nw = nw + 1
         call scheme_item(sp, item, ln, code, pars(nw), ok, errmsg)
         if (.not. ok) return
         ivec(nw) = code
      end do
      if (nw == 0) then
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': empty value list'
         return
      end if
      if (allocated(dist_pars)) deallocate (dist_pars)
      call move_alloc (pars, dist_pars)
      ok = .true.
   end subroutine parse_elist

   ! scheme_item(sp, item, ...) - one list item: WORD or WORD[PARAMS] -> the
   !              word's enum code + the item's staged pair row (form checks
   !              here: word present, bracket closed once, no nesting, no
   !              trailing text; the pairs themselves split in parse_par_pairs)
   subroutine scheme_item(sp, item, ln, code, prow, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: item
      integer, intent(in) :: ln
      integer, intent(out) :: code
      type(frag_par_t), intent(out) :: prow
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=len_trim(item)) :: tail
      character(len=val_len) :: raw
      character(len=key_len) :: word
      integer :: br, lt
      ok = .false.
      errmsg = ''
      raw = ''
      br = index(item, '[')
      if (br == 0) then
         word = to_upper(trim(item))
      else
         if (sp%enum_id /= enum_dist) then     ! brackets are the dist family's form
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': illegal value "'//trim(item)// &
                     '" (legal: '//trim(legal_words(sp%enum_id))//')'
            return
         end if
         word = to_upper(adjustl(trim(item(:br-1))))
         if (len_trim(word) == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': no scheme word before the '// &
                     'bracket in "'//trim(item)//'"'
            return
         end if
         tail = adjustl(trim(item(br+1:)))
         lt = len_trim(tail)
         if (lt == 0 .or. tail(lt:lt) /= ']') then
            ! the whole value balanced (step 1), so the bracket closed early:
            ! text trails the closing bracket
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': trailing text after the '// &
                     'closing bracket in "'//trim(item)//'"'
            return
         end if
         raw = tail(:lt-1)
         if (len_trim(raw) == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': empty bracket in "'//trim(item)// &
                     '" (a word without parameters is written bare)'
            return
         end if
         if (scan(raw, '[') > 0 .or. scan(raw, ']') > 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': nested bracket in "'//trim(item)// &
                     '" (bracket parameters are flat NAME=VALUE pairs)'
            return
         end if
      end if
      code = enum_lookup(sp%enum_id, word)
      if (code == 0) then
         errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': illegal value "'//trim(item)// &
                  '" (legal: '//trim(legal_words(sp%enum_id))//')'
         return
      end if
      prow%word = word
      if (len_trim(raw) > 0) then
         call parse_par_pairs(sp, item, raw, ln, prow, ok, errmsg)
         if (.not. ok) return
      end if
      ok = .true.
   end subroutine scheme_item

   ! parse_par_pairs(sp, item, raw, ...) - split the bracket content on commas
   !              into NAME=VALUE pairs (uppercased; empty names/values,
   !              non-identifier names, over-wide names and duplicate names are
   !              named form errors)
   subroutine parse_par_pairs(sp, item, raw, ln, prow, ok, errmsg)
      type(key_spec_t), intent(in) :: sp
      character(len=*), intent(in) :: item, raw
      integer, intent(in) :: ln
      type(frag_par_t), intent(inout) :: prow
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg
      character(len=len_trim(raw)) :: rest, pair
      character(len=key_len) :: name
      character(len=val_len) :: val
      integer :: np, ic, eq, j, k, m
      ok = .false.
      np = 1 + count([(raw(j:j) == ',', j = 1, len_trim(raw))])
      allocate (prow%pair(np))
      prow%n_pair = 0
      rest = adjustl(raw)
      do while (len_trim(rest) > 0)
         ic = index(trim(rest), ',')
         if (ic == 0) then
            pair = adjustl(trim(rest))
            rest = ''
         else
            pair = adjustl(rest(:ic-1))
            rest = adjustl(rest(ic+1:))
         end if
         if (len_trim(pair) == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': empty parameter in "'// &
                     trim(item)//'" (bracket parameters are NAME=VALUE)'
            return
         end if
         eq = index(pair, '=')
         if (eq == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': bracket item "'//trim(pair)// &
                     '" of "'//trim(item)//'" is not a NAME=VALUE parameter'
            return
         end if
         name = to_upper(adjustl(trim(pair(:eq-1))))
         val = to_upper(adjustl(trim(pair(eq+1:))))
         if (len_trim(name) == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': bracket item "'//trim(pair)// &
                     '" of "'//trim(item)//'" is not a NAME=VALUE parameter'
            return
         end if
         if (len_trim(name) > key_len) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': parameter name "'//trim(name)// &
                     '" exceeds the name width '//trim(i2s(key_len))//' in "'//trim(item)//'"'
            return
         end if
         if (len_trim(val) == 0) then
            errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': empty parameter value in "'// &
                     trim(item)//'" (the parameter '//trim(name)//' would receive nothing)'
            return
         end if
         do k = 1, len_trim(name)           ! identifier charset: A-Z/0-9/_ (folded)
            m = iachar(name(k:k))
            if (.not. ((m >= iachar('A') .and. m <= iachar('Z')) .or. &
                       (m >= iachar('0') .and. m <= iachar('9')) .or. &
                       m == iachar('_'))) then
               errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': bad parameter name "'// &
                        trim(name)//'" in "'//trim(item)//'" (names are letters/digits/underscore)'
               return
            end if
         end do
         do j = 1, prow%n_pair
            if (trim(prow%pair(j)%name) == trim(name)) then
               errmsg = 'key '//trim(sp%name)//trim(line_tag(ln))//': duplicate parameter '// &
                        trim(name)//' in "'//trim(item)//'" (one value per carrier)'
               return
            end if
         end do
         prow%n_pair = prow%n_pair + 1
         prow%pair(prow%n_pair)%name = trim(name)
         prow%pair(prow%n_pair)%val = trim(val)
      end do
      ok = .true.
   end subroutine parse_par_pairs

   ! count_top_commas(s) - commas at bracket depth 0 [-] (the list separators)
   pure integer function count_top_commas(s)
      character(len=*), intent(in) :: s
      integer :: i, d
      count_top_commas = 0
      d = 0
      do i = 1, len_trim(s)
         select case (s(i:i))
         case ('['); d = d + 1
         case (']'); d = d - 1
         case (','); if (d == 0) count_top_commas = count_top_commas + 1
         end select
      end do
   end function count_top_commas

   ! top_comma(s) - position of the first depth-0 comma (0 = none)
   pure integer function top_comma(s)
      character(len=*), intent(in) :: s
      integer :: i, d
      top_comma = 0
      d = 0
      do i = 1, len_trim(s)
         select case (s(i:i))
         case ('['); d = d + 1
         case (']'); d = d - 1
         case (',')
            if (d == 0) then
               top_comma = i
               return
            end if
         end select
      end do
   end function top_comma

   !------------------------------------------------------------------
   ! private helpers (pure string/table mechanics - no system knowledge)
   !------------------------------------------------------------------
   subroutine reset_buffer()
      if (allocated(kv_tbl)) deallocate (kv_tbl)
      allocate (kv_tbl(size(key_tab) + size(dispatch_tbl) + n_declared))   ! one row per
      n_kv = 0                           ! registry key PLUS one per dispatch row (channel-
                                         ! level member keys stage too) PLUS one per declared
                                         ! container key; both fill sites guard their increment
                                         ! with a capacity named abort
      if (allocated(dist_pars)) deallocate (dist_pars)   ! fresh pair table per parse
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
      case (enum_spectrum)
         do i = 1, size(spectrum_tbl)
            if (trim(spectrum_tbl(i)%word) == trim(word)) then
               enum_lookup = spectrum_tbl(i)%code
               exit
            end if
         end do
      case (enum_scattering)
         do i = 1, size(scattering_tbl)
            if (trim(scattering_tbl(i)%word) == trim(word)) then
               enum_lookup = scattering_tbl(i)%code
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
      case (enum_spectrum)
         do i = 1, size(spectrum_tbl)
            if (i > 1) s = trim(s)//', '
            s = trim(s)//trim(spectrum_tbl(i)%word)
         end do
      case (enum_scattering)
         do i = 1, size(scattering_tbl)
            if (i > 1) s = trim(s)//', '
            s = trim(s)//trim(scattering_tbl(i)%word)
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
