!=====================================================================
! recorder.f90 - L4 recording: column-style recorder of trajectory frames / final states
!                (column = name + extractor registration mechanism; column selection by
!                preset level or list - system knowledge parameterized via roster/list)
! (see previous line)
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
!=====================================================================
module recorder
   use state,         only: state_t
   use config,        only: observables
   use config_roster, only: roster
   use force_facade,  only: force_eval
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
      procedure(col_val_i), pointer, nopass :: val => null() ! extractor
   end type col_entry_t
   type(col_entry_t) :: col_tbl(max_col)     ! column registry (module-private)
   integer :: n_col = 0                      ! registered column count [column]
   integer :: sel_idx(max_col) = 0           ! registry indices of the selected columns [-]
   integer :: n_sel = 0                      ! selected column count [column] (rec_select assembly)
   logical :: sec_qp  = .false.              ! all-atom Q/P flattened-section gate (rec_level>=2)
   logical :: sec_frc = .false.              ! per-atom force section gate (rec_level>=5)

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
   subroutine rec_reg_col(name, proc)
      character(len=*), intent(in) :: name  ! column name (obs_list list key)
      procedure(col_val_i) :: proc          ! column extractor
      ! Pseudocode:
      ! 1. Duplicate name or full table -> loud-stop (naming the column)
      ! 2. Register the row {name, proc} (executed in the assembly phase - the driver / member
      ! inits each register rows)
      ! 3. Framework-preregistered columns (called first when this module is assembled):
      !    t_kin (sum p^2/2*roster%mass/e_conv,
      !    the old T convention ENERGY.f:58-66), v_pot/e_tot (via force_eval(sta, e)
      !    - the observables%e_kin/e_pot/e_tot switches each map to a column, the old
      !    GWRITE.f90:31
      !    T/V/H line); electronic observable columns (E14, spec §4 - the module-header
      !    observable-source statute): pop_i reads the diagonal sta%s%rho(i,i) (population),
      !    coh_ij the off-diagonal sta%s%rho(i,j) (coherence, reduced to a real scalar -
      !    modulus or real part, a display convention), occ_k the component of the held
      !    occupation set (occ_get, E16 statistics holder - classical-side statistical
      !    construct; a hop = one component change), n_hop the counter n_hop_get()
      !    (recorder-only read) - these columns never
      !    touch the semantics of a and form the demand side of the blackbox_init
      !    observable-vs-container formal check (spec §4: a column demanding a container no
      !    assembled package member delivers -> loud-stop); preset-group membership in the
      !    rec_level table deferred - selectable via obs_list, finalized in C4;
      !    member/analysis-quantity columns (e_rel, z_min, n_bnc etc.) are registered by
      !    their own members and analysis - C8 wiring
   end subroutine rec_reg_col

   !------------------------------------------------------------------
   ! rec_select(level / list) - column-selection generic (two entries in the interface block
   ! of the module declaration part)
   !
   ! rec_level preset column-group data table (the config%observables%rec_level key - the
   ! parameterized restatement of the GWRITE_LEVEL semantics, GWRITE.f90:11-15):
   !   Level 0 = minimal: frame header (trajectory/step/time) + t_kin/v_pot/e_tot energy line
   !   (old :28/:31)
   !   Level 1 = + structural diagnostic columns (nearest-neighbor distances / coordinates
   !   etc. - the old :63-99 diagnostic section;
   !          system knowledge parameterized into columns, excision warnings (3)(4))
   !   Level 2 = + the all-atom Q/P flattened section sec_qp (default - old :106-117, the
    !          backward-compatible level;
       !          labels from roster%symb - excision warning (2))
       !   Level 5 = + the per-atom force section sec_frc (old :122-127)

   !------------------------------------------------------------------
   ! rec_sel_lvl(level) - select columns by preset level (sets sel_idx and the two section
   ! gates)
   ! Reference: GWRITE.f90:11-15 (four-level semantics)
   subroutine rec_sel_lvl(level)
      integer, intent(in) :: level   ! preset level [-] (0/1/2/5; the default call passes
                                     ! observables%rec_level)
      ! Pseudocode:
      ! 1. Set sel_idx per the data table above (level n includes all columns of levels < n)
      !    and the section gates
      !    sec_qp/sec_frc (open at levels >=2/>=5)
      ! 2. A table column not registered -> loud-stop naming it (assembly-order error -
      !    registration precedes column selection)
   end subroutine rec_sel_lvl

   !------------------------------------------------------------------
   ! rec_sel_lst(list) - select columns by list (comma-separated column names looked up in
   ! the registry one by one)
   ! Reference: config%observables%obs_list (K1 §1.4 - the list key)
   subroutine rec_sel_lst(list)
      character(len=*), intent(in) :: list  ! column-name list (comma-separated;
                                            ! a list containing the reserved section names
                                            ! 'all_qp'/'all_frc'
      ! opens the corresponding flattened section -
      ! finalized in C4)
      ! Pseudocode:
      ! 1. Split on commas and look up each name in col_tbl (unregistered -> loud-stop
      !    naming it)
      ! 2. Matching row indices enter sel_idx in order, counting n_sel; section names open
      !    the section gates as reserved words
   end subroutine rec_sel_lst

   !------------------------------------------------------------------
   ! rec_frame(sta, unit) - write one checkpoint frame (per selected column values + gated
   ! flattened sections)
   ! Reference: GWRITE.f90:25-31 (frame header / energy line), :49-60 (molecule-file step
   ! sections),
   !       :106-127 (the two flattened sections), :145 (flush)
   subroutine rec_frame(sta, unit)
      type(state_t), intent(in) :: sta   ! physical state (read-only snapshot - the post-mqc
      ! self-consistent state: the driver four-phase loop calls rec_frame AFTER
      ! mqc_statistic (spec §8 ordering - the held occupation set already updated for the
      ! step (E16 statistics holder) and rho
      ! refreshed to the declared deliverable state, repairing the one-step occupation lag
      ! of the old checkpoint-timed frames)
      integer, intent(in) :: unit        ! output unit number [-] (driver orchestration - the
      ! old 6/8/1000+NTZ family, unit assignment finalized
      ! in C4)
      ! Pseudocode:
      ! 1. Frame header: trajectory number traj_i / step number step_i / time sta%t (old
      !    :27-28 - the index channel arrives
      !    via the closure host or an optional dummy argument; the skeleton signature keeps
      !    the planned two dummies - finalized in C4)
      ! 2. Per selected column j = 1..n_sel: val = col_tbl(sel_idx(j))%val(sta),
      !    format-write one line/section per column (energy-line convention old :31)
      ! 3. sec_qp open: per atom one q/p six-tuple line (labels roster%symb - excision
      !    warning (2);
      !    old :109-115)
      ! 4. sec_frc open: per atom one force triple (old :122-127; sta%f)
      ! 5. flush(unit) (old :145)
   end subroutine rec_frame

   !------------------------------------------------------------------
   ! rec_final(sta, unit) - write the trajectory final-state row (the record_all column group
   ! - the per-trajectory record
   ! consumed by the statistics phase)
   ! (see previous line)
   ! Reference: HIST.f90:8-24 (record_all: final-state row columns + first-line header), :18
   !       (header column order Ntz/Nb/Site/RX0/RY0/Zmin/Vfinal/Erel/Theta/Ephi)
   subroutine rec_final(sta, unit)
      type(state_t), intent(in) :: sta   ! physical state (final-state snapshot - read-only)
      integer, intent(in) :: unit        ! output unit number [-]
      ! Pseudocode:
      ! 1. When traj_i = 1 write the header line (old HIST.f90:18 - column-name order)
      ! 2. Final-state row column values: e_rel/thta/z_min/v_fin/n_bnc/x0/y0/site - final
      !    state quantities via
      !    analysis/final_state's fin_t (column-extractor registration - C8 wiring;
      !    the member-ization of the old record_all reading globals like AErelS directly);
      !    the site column via plot_hist's
      !    site-assignment tool (lattice parameterized - the parameterized form of surfcell
      !    HIST.f90:313-373,
      !    module-to-module wiring in C8)
      ! 3. flush(unit); this routine only writes lines - unit opening/closing and file-name
      !    layout belong to the driver
      !    (the old VENUS.f90:572 CLOSE site)
   end subroutine rec_final
end module recorder
