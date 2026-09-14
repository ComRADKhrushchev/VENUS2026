!=====================================================================
! recorder.f90 - L4 recording: column-style recorder of trajectory frames / final states
!                (column = name + extractor registration mechanism; column selection by
!                preset level or list - system knowledge parameterized via roster/list)
! Reference: src_VENUS/GWRITE.f90:25-31 (frame header: trajectory number / step number / time
!       + T/V/H energy line),
!       :49-60 (fort.1NNN molecule files: step header + E0(eV)/T(eV)/H(eV) lines + the
!       designated atom's
!       Q/V triple), :106-117 (all-atom Q/P flattened section - GWRITE_LEVEL>=2),
!       :122-127 (per-atom force section - GWRITE_LEVEL>=5), :133-143 (NFQP full-coordinate
!       stdout section), :145 (FLUSH checkpoint flush);
!       src_VENUS/HIST.f90:8-24 (record_all final-state record columns: Ntz/Nb/Site/RX0/RY0/
!       Zmin/Vfinal/Erel/Theta/Ephi - the source of the rec_final column group; the site column
!       via surfcell
!       :19 - lattice parameterized into the plot_hist module);
!       src_VENUS/VENUS.f90:405/:450/:459 (GWRITE call sites - the old orchestration of the
!       NIP checkpoint interval,
!       new version: the driver four-phase loop calls rec_frame directly after mqc_statistic,
!       spec §8 - hooks withdrawn), :572 (the CLOSE(1000+NTZ) unit
!       close-out - orchestration belongs to the driver);
!       src_TEST/interface_TEST.f90:18-30 (POT0 energy convention - the energy-column
!       extractor convention)
! Excision warning (four spots of system knowledge in GWRITE - do not carry over when
! referencing; K2 plan sourcing note):
!   (1) Hardwired atom labels: 'Q(C)='/'V(C)=' treat atom 1 as the system's designated atom
!       (:59-60)
!       -> the observed object is designated via roster fragments / lists;
!   (2) Mass-threshold element names: W(i)<50 -> 'C '/'Au ' labels (:110-114) -> labels
!       taken from roster%symb;
!       (see previous line)
!   (3) nearest-Au diagnostics: distance/coordinates/component scan from the designated atom
!       to the nearest Au (:68-99)
!       -> nearest-neighbor diagnostic column parameterized (target atom and candidate set via
!       the roster, finalized in C8);
!   (4) Hardwired PBC box: derived from the QZB(1,.) extrema + 6.0 Å padding (:69-87) ->
!       box parameters via the roster unit-cell fields (box_lx/box_ly/skew);
!   Also: the E0/23.0605 eV display conversion (:57) is not scattered in this module - the eV
!   convention belongs to the unit contract
!   (the single force_facade ev_kcal point, same source as T4); the old random-number echo
!   section (:30) and
!   the old temperature echo section (:33-45) belong to old task branches - not carried over
! Observable-source statute (E14, spec §4 - the electronic observable column groups):
!   population columns read the sta%s%rho DIAGONALS, coherence columns the OFF-diagonals -
!   the recorder NEVER touches the semantics of a (its row/column index semantics are
!   method knowledge, E11/E12 route; rho is the observable face of the electronic state,
!   refreshed to the method's declared deliverable state by sync_proc). The occupation / hop
!   columns read the elec_facade statistics holder occ_get / n_hop_get (E16: the classical-
!   side statistical constructs relocated out of state_t - sampled constructs have no
!   thing-in-itself status; occ(k) = state of electron k - a hop = one component change
!   of the occupation vector; n_hop is a recorder-only read)
!
! -- C4.5 implementation notes (the landed semantics; criteria: scripts/check_recorder.f90) --
!   Column set, first version (plan D3): framework-preregistered columns = the energy three
!   (t_kin/v_pot/e_tot) + the structural flattened sections sec_qp (q/p six-tuple rows,
!   per atom) / sec_frc (per-atom force rows) + the frame header (trajectory/step/time).
!   pop_i / coh_ij electronic columns are NOT registered yet: their index range arrives
!   with the class-B method declaration (a row/column count is member knowledge, E11/E12);
!   the two index-free electronic constructs land now as the occ / n_hop placeholder
!   columns (E16/E19). Member/analysis-quantity columns (e_rel, z_min, n_bnc, site ...) are
!   registered by their own members and analysis - C8 wiring.
!   N3.1 four-level preset map (DATA, not branches - the GWRITE_LEVEL semantics of
!   GWRITE.f90:11-15): every registry row carries lvl_min, the minimum preset level that
!   includes the column; rec_sel_lvl selects by a single scan over that field, and the two
!   section gates open at the lvl_qp / lvl_frc parameters. Landed map:
!     levels 0 / 1 -> the energy three (minimal frame: header + energy line; the old
!                     level-1 diagnostic block was excised system knowledge, warnings (3)(4)
!                     - no generic level-1 columns exist)
!     level  2     -> + sec_qp (old :106-117; the backward-compatible default)
!     level  5     -> + sec_frc (old :122-127)
!   The level domain is the closed set {0,1,2,5} (the C4.4 GWRITE_LEVEL input domain);
!   anything else loud-stops. occ / n_hop carry lvl_min = huge (no preset group - they are
!   obs_list-only until their group membership is legislated).
!   Coexistence semantics of the two selection entries (plan N3 wording "coexisting",
!   realized as UNION): the level entry SETS the selection - it replaces any previous
!   selection with the preset column group and sets both section gates by its level; the
!   list entry APPENDS its named columns to the current selection (registry duplicates are
!   not re-listed; lookup is case-insensitive) and opens the reserved section gates
!   'all_qp' / 'all_frc'. Driver pattern: rec_select(observables%rec_level) first, then
!   rec_select(observables%obs_list) when the list key is non-empty - the recorded set is
!   the union of the preset group and the named list, the gates the union of both entries.
!   Units (the unit contract): the energy columns report kcal/mol - t_kin =
!   sum p^2/(2*roster%mass)/e_conv (old ENERGY.f:58-66 convention), v_pot / e_tot via
!   force_eval (the kcal/mol e_tot of the force_facade N1 contract; v_pot = e_tot - t_kin).
!   v_pot / e_tot call force_eval on a private COPY of the snapshot (force_eval takes an
!   inout state; the recorded snapshot itself stays read-only).
!   E16/E19 placeholder statute: occ / n_hop are registered-but-INERT - the elec_facade
!   getters are no-op skeletons in C4 (occ_get would hand back an unallocated set,
!   n_hop_get an undefined value), so the extractors return a defined 0.0 placeholder and
!   the getter call lands when the statistics holder goes real (C-phase); their values are
!   a phase-V posterior.
!   Frame-index channel: the signature froze two dummies (sta, unit), so the trajectory /
!   step counters live as module-private state (the "closure host" of the skeleton note):
!   step_i advances per rec_frame call, traj_i advances per rec_final call (which also
!   resets step_i - one trajectory of frames then one final row, the driver loop shape).
!   rec_final writes the GENERIC final-state row (traj, t, selected column values; the
!   name header once, at traj_i = 1 - old HIST.f90:18 shape). The member columns of the
!   old record_all (e_rel/z_min/site ...) are C8 wiring through final_state/fin_t.
!   Assembly-order guards (loud-stop, formal): recording before any rec_select; recording
!   with an unassembled roster (mass/symb unallocated) or a state whose size is not an
!   EXACT fit to the roster (both directions - col_tkin walks roster%mass up to
!   size(sta%p)/3, so an oversized state would otherwise read mass out of bounds) - the
!   two checks cover every extractor / section read the frame path makes.
!   Input-key consumption note: the E_KIN / E_POT / E_TOT input keys land in
!   observables%e_kin / e_pot / e_tot (C4.4 map_slots) and no code reads them yet - the
!   energy columns here register at lvl_min >= 0 regardless of those flags; their
!   consumption destination is the C8 / run-statistics wiring (routed keys, not dead).
!=====================================================================
module recorder
   use state,         only: state_t
   use config,        only: observables
   use config_roster, only: roster
   use force_facade,  only: force_eval
   use consts,        only: e_conv
   implicit none
   private
   public :: col_val_i, rec_reg_col, rec_select, rec_frame, rec_final

   ! Column-extractor interface: extract one scalar observable from the physical state
   ! (the measurement class consuming evolution services - energy columns obtain the potential
   ! energy internally via force_eval, a legal direction; spec §2.2 measurement-class
   ! definition)
   abstract interface
      real(8) function col_val_i(sta)
         import :: state_t
         type(state_t), intent(in) :: sta   ! physical state (read-only snapshot)
      end function col_val_i
   end interface

   ! Column registry (column = name + extractor; registered at assembly time, read-only
   ! during recording)
   integer, parameter :: max_col = 64        ! registry capacity [column]
   type :: col_entry_t
      character(len=16) :: name = ''         ! column name (obs_list list key; 1-16 characters)
      integer :: lvl_min = huge(1)           ! N3.1 preset map: minimum preset level including
                                             ! this column (huge = no preset group - obs_list-only)
      procedure(col_val_i), pointer, nopass :: val => null() ! extractor
   end type col_entry_t
   type(col_entry_t) :: col_tbl(max_col)     ! column registry (module-private)
   integer :: n_col = 0                      ! registered column count [column]
   integer :: sel_idx(max_col) = 0           ! registry indices of the selected columns [-]
   integer :: n_sel = 0                      ! selected column count [column] (rec_select assembly)
   logical :: sec_qp  = .false.              ! all-atom Q/P flattened-section gate (rec_level>=2)
   logical :: sec_frc = .false.              ! per-atom force section gate (rec_level>=5)
   logical :: sel_done = .false.             ! rec_select ran (assembly-order guard for the
                                             ! write entries)
   logical :: fw_done = .false.              ! framework preregistration flag (lazy, first use)
   integer :: traj_i = 1                     ! trajectory counter [-] (frame header / final row)
   integer :: step_i = 0                     ! per-trajectory frame counter [-]

   ! N3.1 section-gate levels (data companions of the lvl_min field; rec_sel_lvl derives
   ! both gates from these constants - no branch table)
   integer, parameter :: lvl_qp  = 2         ! sec_qp opens at this preset level
   integer, parameter :: lvl_frc = 5         ! sec_frc opens at this preset level
   ! the closed preset-level domain (the C4.4 GWRITE_LEVEL input set)
   integer, parameter :: lvl_domain(4) = [ 0, 1, 2, 5 ]

   ! Two entries of the column-selection generic: preset level or list (the four-level
   ! GWRITE_LEVEL semantics of GWRITE.f90:11-15
   ! paired with the config%observables%rec_level/obs_list keys - the K1 §1.4 observable
   ! settings slot)
   interface rec_select
      module procedure rec_sel_lvl   ! by preset level (the rec_level key)
      module procedure rec_sel_lst   ! by list (the obs_list key - comma-separated column names)
   end interface rec_select
contains
   !------------------------------------------------------------------
   ! rec_reg_col(name, proc) - register an observable column (one row each of name + extractor)
   ! Reference: docs/final_completion_plan.md §5 (the column-registration mechanism - C4
   ! recorder);
   !       interface_TEST.f90:18-30 (energy-column convention - reference for the
   !       framework-preregistered columns)
   ! Implemented: the framework columns preregister themselves lazily on first use (fw_reg
   ! below - registration at assembly precedes selection because every public entry calls
   ! fw_reg first); an external row enters with lvl_min = huge (no preset group - external
   ! columns are obs_list-only; preset membership is framework legislation, N3.1). Duplicate
   ! name or full table loud-stops naming the column
   subroutine rec_reg_col(name, proc)
      character(len=*), intent(in) :: name  ! column name (obs_list list key)
      procedure(col_val_i) :: proc          ! column extractor
      call fw_reg()
      call reg_row(name, proc, huge(1))
   end subroutine rec_reg_col

   !------------------------------------------------------------------
   ! rec_select(level / list) - column-selection generic (two entries in the interface block
   ! of the module declaration part)
   !
   ! rec_level preset column-group data table (the config%observables%rec_level key - the
   ! parameterized restatement of the GWRITE_LEVEL semantics, GWRITE.f90:11-15):
   !   Level 0 / 1 = minimal: frame header (trajectory/step/time) + t_kin/v_pot/e_tot energy
   !   line (old :28/:31; no generic level-1 diagnostic columns exist - that old block was
   !   excised system knowledge, excision warnings (3)(4))
   !   Level 2 = + the all-atom Q/P flattened section sec_qp (default - old :106-117, the
   !          backward-compatible level; labels from roster%symb - excision warning (2))
   !   Level 5 = + the per-atom force section sec_frc (old :122-127)
   !
   ! Coexistence: see the C4.5 implementation notes in the module header (union semantics).

   !------------------------------------------------------------------
   ! rec_sel_lvl(level) - select columns by preset level (sets sel_idx and the two section
   ! gates; replaces any previous selection - the level entry is the base of the union)
   ! Reference: GWRITE.f90:11-15 (four-level semantics)
   subroutine rec_sel_lvl(level)
      integer, intent(in) :: level   ! preset level [-] (0/1/2/5; the default call passes
                                     ! observables%rec_level)
      integer :: j
      call fw_reg()
      if (all(level /= lvl_domain)) then
         call stop_rec('rec_sel_lvl', 'unknown preset level - the recording-level set is '// &
                       '0/1/2/5 (GWRITE_LEVEL semantics, N3.1)')
      end if
      n_sel = 0
      do j = 1, n_col
         if (col_tbl(j)%lvl_min <= level) then
            n_sel = n_sel + 1
            sel_idx(n_sel) = j
         end if
      end do
      sec_qp  = (level >= lvl_qp)
      sec_frc = (level >= lvl_frc)
      sel_done = .true.
   end subroutine rec_sel_lvl

   !------------------------------------------------------------------
   ! rec_sel_lst(list) - select columns by list (comma-separated column names looked up in
   ! the registry one by one; appends to the current selection - the list entry is the
   ! additive half of the union)
   ! Reference: config%observables%obs_list (K1 §1.4 - the list key)
   subroutine rec_sel_lst(list)
      character(len=*), intent(in) :: list  ! column-name list (comma-separated;
                                            ! the reserved section names 'all_qp'/'all_frc'
                                            ! open the corresponding flattened section)
      character(len=len(list)) :: rest
      character(len=16) :: tok
      integer :: ic, j
      call fw_reg()
      if (len_trim(list) == 0) then
         sel_done = .true.                   ! empty list key = no columns to add (legal no-op)
         return
      end if
      rest = adjustl(list)
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
            call stop_rec('rec_sel_lst', 'empty item in the column list (obs_list grammar: '// &
                          'comma-separated names, no empty items)')
         end if
         if (to_upper(trim(tok)) == 'ALL_QP') then
            sec_qp = .true.
         else if (to_upper(trim(tok)) == 'ALL_FRC') then
            sec_frc = .true.
         else
            j = col_lookup(tok)
            if (j == 0) then
               call stop_rec('rec_sel_lst', 'unknown column name "'//trim(tok)// &
                             '" (not in the recorder column registry)')
            end if
            if (.not. any(sel_idx(1:n_sel) == j)) then
               if (n_sel >= max_col) then
                  call stop_rec('rec_sel_lst', 'selection table full ('// &
                                'too many selected columns)')
               end if
               n_sel = n_sel + 1
               sel_idx(n_sel) = j
            end if
         end if
      end do
      sel_done = .true.
   end subroutine rec_sel_lst

   !------------------------------------------------------------------
   ! rec_frame(sta, unit) - write one checkpoint frame (per selected column values + gated
   ! flattened sections)
   ! Reference: GWRITE.f90:25-31 (frame header / energy line), :49-60 (molecule-file step
   ! sections),
   !       :106-127 (the two flattened sections), :145 (flush)
   ! Implemented frame layout (generic quantities only - zero system knowledge):
   !   '# frame traj=<t> step=<s> t=<time>'          frame header (old :27-28, one line)
   !   '<name>=<val> <name>=<val> ...'               selected column values, selection order
   !   '# sec_qp begin' / per-atom rows / end        gated at rec_level>=2 (old :106-117)
   !   '# sec_frc begin' / per-atom rows / end       gated at rec_level>=5 (old :122-127)
   ! atom row = '<index> <symb> <3 reals> [<3 reals>]' (labels from roster%symb - excision
   ! fix (2); no designated atom, excision fix (1)); flush(unit) at the end (old :145).
   ! Units: t = sta%t [10 fs] as held; energies [kcal/mol] per the unit contract above
   subroutine rec_frame(sta, unit)
      type(state_t), intent(in) :: sta   ! physical state (read-only snapshot - the post-mqc
      ! self-consistent state: the driver four-phase loop calls rec_frame AFTER
      ! mqc_statistic (spec §8 ordering - the held occupation set already updated for the
      ! step (E16 statistics holder) and rho refreshed to the declared deliverable state,
      ! repairing the one-step occupation lag of the old checkpoint-timed frames)
      integer, intent(in) :: unit        ! output unit number [-] (driver orchestration - unit
      ! assignment / opening / closing belongs to the driver, this routine only writes)
      integer :: j, k, nat
      call fw_reg()
      call check_ready(sta, 'rec_frame')
      nat = size(roster%mass)
      step_i = step_i + 1
      write (unit, '(a,i0,a,i0,a,es24.16)') '# frame traj=', traj_i, ' step=', step_i, ' t=', sta%t
      if (n_sel > 0) then
         do j = 1, n_sel
            if (j == 1) then
               write (unit, '(a)', advance='no') trim(trim(col_tbl(sel_idx(j))%name)//'='// &
                  num_str(col_tbl(sel_idx(j))%val(sta)))
            else
               write (unit, '(a)', advance='no') trim(' '//trim(col_tbl(sel_idx(j))%name)//'='// &
                  num_str(col_tbl(sel_idx(j))%val(sta)))
            end if
         end do
         write (unit, '(a)') ''
      end if
      if (sec_qp) then
         write (unit, '(a)') '# sec_qp begin'
         do k = 1, nat
            write (unit, '(i0,1x,a,1x,6es24.16)') k, roster%symb(k), &
               sta%q(3*k-2:3*k), sta%p(3*k-2:3*k)
         end do
         write (unit, '(a)') '# sec_qp end'
      end if
      if (sec_frc) then
         write (unit, '(a)') '# sec_frc begin'
         do k = 1, nat
            write (unit, '(i0,1x,a,1x,3es24.16)') k, roster%symb(k), sta%f(3*k-2:3*k)
         end do
         write (unit, '(a)') '# sec_frc end'
      end if
      flush (unit)
   end subroutine rec_frame

   !------------------------------------------------------------------
   ! rec_final(sta, unit) - write the trajectory final-state row (the record_all column group
   ! - the per-trajectory record consumed by the statistics phase)
   ! Reference: HIST.f90:8-24 (record_all: final-state row columns + first-line header), :18
   !       (header column-name order)
   ! Implemented (C4 first version - the GENERIC row): '# final: traj t <column names>'
   ! written once (traj_i = 1), then one row per trajectory = traj_i, sta%t, the selected
   ! column values at the final state; afterwards traj_i advances and step_i resets (the
   ! frame/final counters of the closure host). The member columns of the old record_all
   ! (e_rel/thta/z_min/v_fin/n_bnc/x0/y0/site via final_state fin_t + the plot_hist site
   ! tool, lattice parameterized) are C8 wiring - the old VENUS.f90:572 CLOSE site belongs
   ! to the driver in any case
   subroutine rec_final(sta, unit)
      type(state_t), intent(in) :: sta   ! physical state (final-state snapshot - read-only)
      integer, intent(in) :: unit        ! output unit number [-]
      integer :: j
      call fw_reg()
      call check_ready(sta, 'rec_final')
      if (traj_i == 1) then
         write (unit, '(a)', advance='no') '# final: traj t'
         do j = 1, n_sel
            write (unit, '(a)', advance='no') ' '//trim(col_tbl(sel_idx(j))%name)
         end do
         write (unit, '(a)') ''
      end if
      write (unit, '(i0,1x,es24.16)', advance='no') traj_i, sta%t
      do j = 1, n_sel
         write (unit, '(1x,es24.16)', advance='no') col_tbl(sel_idx(j))%val(sta)
      end do
      write (unit, '(a)') ''
      flush (unit)
      traj_i = traj_i + 1
      step_i = 0
   end subroutine rec_final

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! fw_reg() - lazy framework preregistration (runs once, before any registry use; the
   ! assembly-time "called first" of the skeleton note realized as first-use initialization)
   ! Framework columns (D3 first version): the energy three at preset level 0 (minimal
   ! frame) + the occ/n_hop placeholders (lvl_min = huge - obs_list-only)
   subroutine fw_reg()
      if (fw_done) return
      fw_done = .true.
      call reg_row('t_kin', col_tkin, 0)
      call reg_row('v_pot', col_vpot, 0)
      call reg_row('e_tot', col_etot, 0)
      call reg_row('occ',   col_occ,  huge(1))
      call reg_row('n_hop', col_nhop, huge(1))
   end subroutine fw_reg

   ! reg_row(name, proc, lvl) - append one registry row (shared by the public entry and the
   ! framework preregistration); duplicate name / full table / overlong name loud-stop
   subroutine reg_row(name, proc, lvl)
      character(len=*), intent(in) :: name
      procedure(col_val_i) :: proc
      integer, intent(in) :: lvl
      if (len_trim(name) == 0 .or. len_trim(name) > 16) then
         call stop_rec('rec_reg_col', 'illegal column name length ("'//trim(name)// &
                       '"; 1-16 characters)')
      end if
      if (col_lookup(name) /= 0) then
         call stop_rec('rec_reg_col', 'duplicate column name "'//trim(name)// &
                       '" (already in the recorder column registry)')
      end if
      if (n_col >= max_col) then
         call stop_rec('rec_reg_col', 'column registry full ('// &
                       'too many registered columns; max_col exceeded)')
      end if
      n_col = n_col + 1
      col_tbl(n_col)%name = name
      col_tbl(n_col)%lvl_min = lvl
      col_tbl(n_col)%val => proc
   end subroutine reg_row

   ! col_lookup(name) - case-insensitive name lookup; 0 = not found
   integer function col_lookup(name)
      character(len=*), intent(in) :: name
      integer :: j
      col_lookup = 0
      do j = 1, n_col
         if (to_upper(trim(col_tbl(j)%name)) == to_upper(trim(name))) then
            col_lookup = j
            return
         end if
      end do
   end function col_lookup

   ! check_ready(sta, who) - the two formal assembly-order guards of the write path:
   ! selection made, roster assembled and size-EXACTLY-matched to the state (equality,
   ! not a one-way lower bound: col_tkin walks roster%mass(k) for k = 1..size(sta%p)/3,
   ! so a state LARGER than the roster would read mass out of bounds - the guard must
   ! fire by name ahead of that, in both directions; exact fit stays legal)
   subroutine check_ready(sta, who)
      type(state_t), intent(in) :: sta
      character(len=*), intent(in) :: who
      integer :: nat
      if (.not. sel_done) then
         call stop_rec(who, 'called before any rec_select (assembly-order error: '// &
                       'column selection precedes recording)')
      end if
      if (.not. allocated(roster%mass) .or. .not. allocated(roster%symb)) then
         call stop_rec(who, 'roster not assembled (roster_load must run before recording)')
      end if
      nat = size(roster%mass)
      if (size(sta%q) /= 3*nat .or. size(sta%p) /= 3*nat .or. size(sta%f) /= 3*nat) then
         call stop_rec(who, 'state/roster size mismatch (state carries '//trim(i8_str(size(sta%q)))// &
                       ' q-dofs, the roster carries '//trim(i8_str(nat))//' atoms - exact fit required)')
      end if
   end subroutine check_ready

   ! i8_str(i) - integer rendered as a string (the size-mismatch loud-stop names both sides)
   pure function i8_str(i) result(s)
      integer, intent(in) :: i
      character(len=12) :: s
      write (s, '(i0)') i
      s = trim(adjustl(s))
   end function i8_str

   ! stop_rec(who, msg) - the loud-stop channel (same semantics as config_roster /
   ! sysdef_facade: name the caller, print the named message, STOP 1)
   subroutine stop_rec(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': recording error: '//trim(msg)
      write (0, '(a)') trim(who)//': loud-stop (the recording is unusable)'
      stop 1
   end subroutine stop_rec

   ! num_str(v) - a real rendered as a blank-free ES string (the column line's
   ! name=value tokens stay single-space-separated: an unadjusted ES field would leave
   ! the format's leading blank inside the token)
   pure function num_str(v) result(s)
      real(8), intent(in) :: v
      character(len=32) :: s
      write (s, '(es24.16)') v
      s = trim(adjustl(s))
   end function num_str

   ! to_upper(s) - ASCII uppercase of a trimmed string (the obs_list lookup is
   ! case-insensitive; input.f90 stores the list raw)
   pure function to_upper(s)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: to_upper
      integer :: i, ic
      do i = 1, len(s)
         ic = iachar(s(i:i))
         if (ic >= iachar('a') .and. ic <= iachar('z')) then
            to_upper(i:i) = achar(ic - 32)
         else
            to_upper(i:i) = s(i:i)
         end if
      end do
   end function to_upper

   !------------------------------------------------------------------
   ! framework column extractors (the preregistered set)
   !------------------------------------------------------------------

   ! t_kin - kinetic energy T = sum_atoms p.p/(2 m) / e_conv [kcal/mol] (old ENERGY.f:58-66
   ! convention; masses from the roster - read-only)
   real(8) function col_tkin(sta) result(v)
      type(state_t), intent(in) :: sta
      integer :: k
      v = 0.0d0
      do k = 1, size(sta%p)/3
         v = v + dot_product(sta%p(3*k-2:3*k), sta%p(3*k-2:3*k))/(2.0d0*roster%mass(k))
      end do
      v = v/e_conv
   end function col_tkin

   ! v_pot - potential energy via force_eval [kcal/mol]: V = H - T on a private copy of the
   ! snapshot (force_eval takes inout state; C4: the facade is a no-op skeleton, so the
   ! value stays at its initialized 0 until C2(interface) - a phase-V numeric posterior)
   real(8) function col_vpot(sta) result(v)
      type(state_t), intent(in) :: sta
      type(state_t) :: wrk
      real(8) :: e
      wrk = sta
      e = 0.0d0
      call force_eval(wrk, e)
      v = e - col_tkin(sta)
   end function col_vpot

   ! e_tot - total energy H = T(P) + V via force_eval [kcal/mol] (the facade N1 contract;
   ! same copy / no-op-skeleton note as col_vpot)
   real(8) function col_etot(sta) result(v)
      type(state_t), intent(in) :: sta
      type(state_t) :: wrk
      wrk = sta
      v = 0.0d0
      call force_eval(wrk, v)
   end function col_etot

   ! occ - E16/E19 placeholder column: the occupation set is held in the elec_facade
   ! statistics holder, read through occ_get - which is a NO-OP SKELETON in C4 (it would
   ! hand back an unallocated set), so this extractor is registered-but-inert: a defined
   ! placeholder until the holder goes real, its value a phase-V posterior
   real(8) function col_occ(sta) result(v)
      type(state_t), intent(in) :: sta
      v = 0.0d0
   end function col_occ

   ! n_hop - E16/E19 placeholder column: the hop counter via n_hop_get (same no-op-skeleton
   ! note as col_occ - n_hop_get would return an undefined value; inert placeholder)
   real(8) function col_nhop(sta) result(v)
      type(state_t), intent(in) :: sta
      v = 0.0d0
   end function col_nhop
end module recorder
