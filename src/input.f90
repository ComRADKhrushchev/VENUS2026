!=====================================================================
! input.f90 - L0 input: KEYWORD=VALUE framework (line-by-line parsing -> control/config slot mapping
!             and member-parameter dispatch; grammar finalized in C4 - the skeleton sets up the frame and mapping slots first)
! Reference: venus2026/input_parser.f90:67-152 (load_input: blank lines and '#'/'!' comment
!       skipping :107-108, '=' splitting key/value :125-130, trailing ',' continuation
!       absorption :134-145, key uppercasing :170, grammar detection by '=' presence in the first data line :111-120),
!       :196-239 (extract_nth, comma-separated nth value), :268-458 (get_int /
!       get_real/get_str/get_real_arr/get_int_arr/has_keyword/
!       count_keyword_values/get_*_element getters), :489-502 (set_default -
!       default keys are injected only when the user has not set them);
!       venus2026/param_mapping.f90:34-160 (enumeration mapping-table family: map_elec_method
!       :34-55, map_surface_model :60-75, map_task :80-101, map_integrator
!       :106-121+, map_init_sampling, map_old_nsurf - string enumerations -> integer codes,
!       illegal values stop with the name printed);
!       venus2026/presets.f90:33-100 (apply_preset: MODEL key -> five-preset dispatch,
!       preset body = a set_default default-key table);
!       venus2026/venus_input.f90:187-204 (loading + MODEL preset + old sequential grammar
!       rejection site :200-204), :209-231 (TITLE x2 / ELEC_METHOD new keys and CALTYP old key
!       dual-track fallback :214-222 / NATOMS / GWRITE_LEVEL), :248-256 (ATOM_MASSES mass
!       table + two-value broadcast rule), :403-413 (INIT_SAMPLING_A new key and NACTA old integer
!       dual-track fallback - prototype of member-parameter dispatch), :495-506 (NACTA branch printout -
!       dispatch-site embryo), :607-681 (per-branch reading of TRV_A and other member parameters), :996-1000 (SKEW /
!       BOXLX/BOXLY unit-cell section);
!       grammar finalized in C4 (K2 plan §1 batch-4 table note) - the new key set, the TERMINATE=/OUTCOME=
!       table grammar and the old-integer-key compatibility table are finalized in C4; this skeleton sets up
!       the frame and mapping slots without locking the key set
!=====================================================================
module input
   use control, only: task, n_traj, title, i_seed
   use config,  only: propagator, electronic, reactants, event, observables
   implicit none
   private
   public :: read_input

   ! Parse staging (KEYWORD=VALUE row table - module-private; the consumption interface for
   ! roster_load and member init is finalized in C4)
   type :: kv_row_t
      character(len=32)  :: key = ''   ! keyword (uppercased; shortened from old MAX_KEY_LEN)
      character(len=256) :: val = ''   ! value section (comma-separated; continuations already merged back)
   end type kv_row_t
   type(kv_row_t), allocatable :: kv_tbl(:)  ! row table (module-private)
   integer :: n_kv = 0                          ! number of rows [row]
contains
   !------------------------------------------------------------------
   ! read_input(file) - read and map the input (parse -> presets -> enumerations -> slot mapping ->
   !                   roster staging -> member-parameter dispatch -> validation)
   ! Reference: venus_input.f90:187-269 (prototype of the read-in main flow); input_parser.f90:67-152
   subroutine read_input(file)
      character(len=*), intent(in) :: file  ! input file name (e.g. 'input_qct.txt' -
                                            ! the old venus_input.f90:189 hardwired name,
                                            ! now passed as a dummy argument)
      ! Pseudocode:
      ! 1. Parse line by line into kv_tbl: skip blank and comment lines ('#'/'!' first character);
      !    split on '=' into key/value; a trailing ',' in the value absorbs continuation lines;
      !    uppercase keys (old input_parser.f90:97-148; line numbers feed the error channel)
      ! 2. Grammar detection and old integer-keyword fallback: no '=' in the first data line -> old
      !    sequential grammar - read back as the old integer-keyword stream (old venus_input.f90:200-204
      !    hard-rejected instead; whether the new version falls back reading or prints a conversion note
      !    is a C4 decision - the skeleton provides the "detect + fallback slot" without hardwiring a policy)
      ! 3. Preset expansion: MODEL key -> preset dispatch (five presets - gas phase / relaxed surface / rigid
      !    surface / full surface / GLO surface); the preset body injects default keys (never overriding
      !    user-set keys - set_default semantics, input_parser.f90:489-502; old presets.f90:33-100)
      ! 4. Enumeration mapping (mapping-table family idea - param_mapping.f90): string enumeration ->
      !    integer code, illegal values stop with the name printed - electron law -> electronic%method (old
      !    map_elec_method), configuration paradigm -> reactants%surface_model, task -> task,
      !    integrator + order -> propagator%id/order (old map_integrator), sampling scheme ->
      !    reactants%dist_scheme (old map_init_sampling); when no new enumeration key matches, the old
      !    integer keys are read directly for compatibility (CALTYP/NSURF/INTEGRATOR/NACTA.NACTB - old
      !    venus_input.f90:214-222/:403-413 dual-track fallback sites)
      ! 5. Control slot mapping: task, n_traj (old NT), title(2) (old TITLE1/2),
      !    i_seed converted to 8 words (old scalar ISEED -> rng seed conversion finalized in C4 - consumed
      !    by rng_init)
      ! 6. config slot mapping:
      !    propagator <- id/order/dt (old INTEGRATOR/LLL/HTMIN family);
      !    electronic <- method/n_surf/pes_main (old CALTYP/NSURF/TEST_PES);
      !    reactants <- n_frag/surface_model/dist_scheme(2)/e_rel/b_max/r_sep
      !    (old NACTA.B/NSURF/SEREL/BMAX/S);
      !    event <- term_expr/outcome_tbl (TERMINATE=/OUTCOME= new grammar keys -
      !    the predicate module-header grammar; wired in C3/C7);
      !    observables <- rec_level/obs_list/e_kin/e_pot/e_tot (old
      !    GWRITE_LEVEL -> rec_level)
      ! 7. Raw roster-section staging: ATOM_MASSES mass table + two-value broadcast rule (if only 2 values
      !    are given they broadcast to the remaining atoms - old venus_input.f90:248-256), fragment-A table
      !    and equilibrium geometry (old :508-518), unit-cell parameters SKEW/BOXLX/BOXLY (old :996-1000);
      !    the staged data is consumed by roster_load (consumption interface finalized in C4 - roster does
      !    not use this module; the dependency is one-way)
      ! 8. Member-parameter dispatch: parameter keys owned by dist_*/pes_* members -> dispatched to each
      !    member's init (member = code + its own data; old TRV_A/TROT_A/AI_A etc. were scattered over the
      !    NACTA branches in venus_input.f90:607-681 - the new version has one central dispatch table,
      !    filled in row by row with the members in C2/C5)
      ! 9. Validation: missing required keys / numeric-domain violations / illegal enumerations ->
      !    loud-stop naming the key and line number (error channel same as pred_parse's ios/errmsg)
   end subroutine read_input
end module input
