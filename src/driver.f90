!=====================================================================
! driver.f90 - program venus: three-phase orchestration (assembly -> evolution -> statistics)
! Design:
!   Assembly order is binding: registry_gen_keys -> read_input -> rng_init ->
!   sysdef_load -> list_atoms_load -> state_create -> registry_gen_all ->
!   registry_gen_init (the container seam FIRST - a load-dependent PES must
!   be computable before any member seam that probes its force, e.g. the
!   spectrum derivation; then buffered member keys -> member inits, dt
!   already final) -> input_audit_declared (every buffered container key was
!   pulled) -> container_init -> rec_select -> hand-off allocation.
!   The key-declaration aggregate runs BEFORE the parse so the dragged-in
!   container's parameter vocabulary precedes the parse admission (a
!   <module>_keys entry is pure vocabulary data - reads no config, no
!   list_atoms; _all and _init keep their established timing).
!   Evolution: per trajectory samp_run -> force_eval -> first frame, then
!   the step loop prop_step -> elec_prop -> mqc_statistic -> rec_frame ->
!   container_term, bounded by MAX_STEPS (termination and outcome naming
!   are the container's judgments). Statistics: one fin_eval pass -> outcome
!   counts + histograms + drift summary; output: trajectories.txt, hist_e_rel.dat, hist_thta.dat (cwd-relative).
!=====================================================================
program venus
   use control,        only: n_traj, i_seed, max_steps
   use config,         only: observables, reactants, electronic, &
                             prop_cfg => propagator  ! alias: the config instance 'propagator'
                             ! clashes with the module name
   use config_atoms,  only: list_atoms, list_atoms_load
   use state,          only: state_t, state_create, state_destroy
   use rng,            only: rng_init
   use propagator,     only: prop_step
   use force_interface,   only: container_init, force_eval, container_term
   use elec_interface,    only: elec_prop, mqc_statistic
   use sysdef_interface,  only: sysdef_load
   use registry_gen,   only: registry_gen_all, registry_gen_init, registry_gen_keys   ! the build-generated
                                   ! aggregate (scripts/prebuild.sh emits it from every
                                   ! reg_*.f90 under methods/ and systems/ - zero
                                   ! hand-written registration lines in the driver)
   use sampler,        only: samp_run, samp_n_dist, samp_incident_word
   use recorder,       only: rec_select, rec_frame, rec_final, rec_n_sel, rec_fin_val
   use final_state,    only: fin_t, fin_eval, container_classify, container_classify_bound
   use density,        only: dens_states, eig_out
   use hist,           only: hist_spec_t, hist_data_t, lat_t, plot_hist, &
                             k_e_hist, k_t_hist
   use input,          only: read_input, input_audit_declared
   implicit none

   ! -- orchestration counters and flags (locals) --
   integer :: traj_i        ! trajectory index [count] (1..n_traj)
   integer :: step_i        ! step index [step] (the MAX_STEPS safety-net counter - the
                            ! driver loop condition step_i < max_steps)
   logical :: term_flag     ! termination flag for this trajectory (written by the
                            ! container_term call of the step loop - the container's judgment)
   real(8) :: e0            ! total energy of the initial state [kcal/mol] (drift reference)
   real(8) :: dt_eff        ! effective step length [10 fs] (buffered once from prop_cfg%dt;
                            ! consumed by elec_prop - the hop-probability channel g_ij(dt)
                            ! needs the step's dt)
   real(8) :: t_cpu(2)      ! timing pair [s] (printed in the statistics section)
   integer :: u_out         ! recording unit number [-] (opened on 'trajectories.txt' at the
                            ! assembly tail, closed at the statistics wrap-up)
   character(len=32) :: outc_name  ! outcome name of this trajectory [-] (the container's
                                   ! classification output)
   integer, parameter :: max_outc = 32   ! outcome-name table capacity [name]
   character(len=32) :: outc_names(max_outc) ! outcome-name table [-] (first-seen insertion
                                             ! order, counts parallel in cnt_outc)
   integer :: n_outc        ! distinct outcome names counted so far [name]
   integer :: n_valid       ! valid-trajectory count [count] (a final row archived)
   integer :: k             ! scratch index [-]
   integer :: row_outc      ! outcome-table row scratch [-]
   integer :: ios           ! open-status scratch [-] (the recording-unit go-live check below)
   real(8), allocatable :: e0_traj(:)     ! per-trajectory reference energy [kcal/mol] (the
                                          ! head force_eval output - the drift reference)
   real(8), allocatable :: e_fin_traj(:)  ! per-trajectory final e_tot [kcal/mol] (the
                                          ! recorder final-state row)
   real(8), allocatable :: vpot_fin_traj(:) ! per-trajectory final v_pot [kcal/mol] (the
                                          ! recorder final-state row - the histogram v_fin column)
   real(8) :: drift_max     ! max |e_fin - e0| over trajectories [kcal/mol]
   ! -- data structures handed between phases --
   type(state_t) :: sta                      ! physical state (written exclusively in the
                                             ! evolution phase, read-only in statistics)
   type(state_t), allocatable :: sta_fin(:)  ! final-state snapshots per trajectory [count]
                                             ! (evolution -> statistics hand-off)
   type(fin_t) :: fin                        ! final-state quantities (the argument type of
                                             ! the container-classification slot)
   integer, allocatable :: cnt_outc(:)       ! per-outcome counts [count] (name-keyed over the
                                             ! names container_classify returns)
   type(hist_spec_t) :: spec_hist            ! histogram specification (kind, file name, the
                                             ! [0, column max + 1) range rule, 40 bins, normalized)
   type(hist_data_t) :: data_hist            ! trajectory row table (aggregated from fin truth +
                                             ! the recorder final-state rows; system-diagnostic
                                             ! columns stay unaggregated - no generic producer)
   type(lat_t) :: lattice                    ! lattice constants (mapped from the atom list unit
                                             ! cell; consumed by the lattice-aware histogram kinds)

   ! ===== Phase 1: assembly (the order below is binding) =====
   call registry_gen_keys()                        ! 1. key-declaration aggregate (the
                                                   !    dragged-in container's parameter
                                                   !    vocabulary - pure data, the parse
                                                   !    admits these key names)
   call read_input('input_qct.txt')                ! 2. input framework (QCT = quasi-classical
                                                   !    trajectory; KEYWORD=VALUE -> control/config
                                                   !    slots; member and container keys stage raw
                                                   !    for the seams at steps 8-9)
   call rng_init(i_seed(1))                        ! 3. random stream (word 1 of i_seed seeds it)
   call sysdef_load(reactants%system_dir)          ! 4. system definition (folder scan -> parser
                                                   !    dispatch -> mass table; a bad folder
                                                   !    aborts with a named error inside sysdef_load)
   call list_atoms_load()                              ! 5. list_atoms (fragments / masses / equilibrium
                                                   !    placement, from the sysdef products)
   call state_create(sta, size(list_atoms%mass), max(electronic%n_surf, 1))
                                                   ! 6. physical state (natoms from the atom list;
                                                   !    the max() is the N_SURF reconcile - the
                                                   !    input default 0 folds to the adiabatic
                                                   !    one-surface point)
   if (size(sta%q) /= 3*size(list_atoms%mass)) then    ! size-equality guard (drift backstop; both
                                                   ! sides named)
      write (0, '(a,i0,a,i0,a)') 'driver: size guard - the state holds ', size(sta%q), &
         ' coordinates but the atom list mass table holds ', 3*size(list_atoms%mass), &
         ' (assembly-choreography error, formal)'
      write (0, '(a)') 'driver: fatal (the state and the atom list disagree)'
      stop 1
   end if
   call registry_gen_all()                         ! 7. registration aggregate (method packages,
                                                   !    distribution members, container bindings,
                                                   !    and the integrator closed set)
   call registry_gen_init()                        ! 8. the assembly: the container seam FIRST
                                                   !    (parameter load / term / classify / arm /
                                                   !    columns - a load-dependent PES is
                                                   !    computable from here on), then buffered
                                                   !    member keys -> each member's own init (an
                                                   !    unselected member skips inside its own
                                                   !    guard); the spectrum seams in here may
                                                   !    probe the bound force (SPECTRUM_SOURCE
                                                   !    legislation, 2026-10-05); dt was
                                                   !    finalized at step 2, so a member
                                                   !    reconciling against dt reads the final
                                                   !    value - the assembly order is the contract
   call input_audit_declared()                    ! 9. container-key consumption audit (every
                                                   !    buffered declared key pulled exactly once;
                                                   !    a stranded row is a configuration or
                                                   !    authoring error, named here)
   call container_init()                          ! 10. container assembly + formal validation
                                                   !    (legal with no container dragged in - the
                                                   !    degenerate method row still binds;
                                                   !    named-abort authority = force_eval's first
                                                   !    call)
   call rec_select(observables%rec_level)         ! 11a. recording preset group (level entry)
   if (len_trim(observables%obs_list) > 0) then   ! 11b. named-column list entry (the recorded
      call rec_select(trim(observables%obs_list)) !    set is the union of both entries)
   end if
   allocate (sta_fin(n_traj), e0_traj(n_traj), e_fin_traj(n_traj), vpot_fin_traj(n_traj))
                                                   ! 12. inter-phase hand-off structures (the
                                                   !     snapshots + the drift/column pairs the
                                                   !     statistics pass consumes)
   dt_eff = prop_cfg%dt                            !     buffered once (the fixed-step source -
                                                   !     see the local's note)
   call cpu_time(t_cpu(1))                         !     timing start
   open (newunit=u_out, file='trajectories.txt', status='replace', action='write', iostat=ios)
                                                   !     the recording unit goes live on the
                                                   !     registered output file 'trajectories.txt'
                                                   !     (cwd-relative); the iostat channel keeps
                                                   !     an unwritable run directory a NAMED
                                                   !     named abort, never a bare runtime open error
   if (ios /= 0) then
      write (0, '(a,i0,a)') 'driver: the registered output file trajectories.txt cannot be'// &
         ' opened for writing in the run directory (iostat ', ios, ')'
      write (0, '(a)') 'driver: fatal (the recording unit cannot go live - no run)'
      stop 1
   end if
   write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
      'assembly complete: trajectories=', n_traj, ' atoms=', size(list_atoms%mass), &
      ' fragments=', size(list_atoms%frag), ' paradigm=', reactants%surface_model, &
      ' members=', samp_n_dist(), ' init=ok columns=', rec_n_sel()

   ! ===== Phase 2: evolution (trajectory loop x integration loop) =====
   ! Inside the integration loop the direct entries are exactly the four-phase
   ! choreography below; force_eval runs inside the propagator member per the
   ! endpoint-force contract, and the sampler is a trajectory-loop-head entry.
   n_valid = 0
   do traj_i = 1, n_traj
      sta%t = 0.0d0                     ! fresh-trajectory reset: the time origin and the
                                        ! buffer boot key (boot tests read t == 0; the
                                        ! thermalize member saves and restores t around its
                                        ! own equilibration)
      call samp_run(samp_incident_word(), sta) ! per-trajectory initial-state preparation (the
                                        ! per-fragment schemes arrive inside samp_run via the
                                        ! dist_scheme data; the incident member is
                                        ! paradigm-resolved - the trivial single-reactant
                                        ! member is a bitwise no-op)
      call elec_prop(sta, 0.0d0)        ! trajectory-head rho reseed: the prop member's
                                        ! t-backwards guard fires on the reset t (dt=0
                                        ! evolves nothing) - e0 below must reference the
                                        ! SEEDED state, not the previous trajectory's
                                        ! leftover rho (stale-e0 artifact, 2026-10-04)
      call force_eval(sta, e0)          ! initial force + the reference energy of the
                                        ! initial state (the pure force entry)
      e0_traj(traj_i) = e0
      call rec_frame(sta, u_out)        ! first frame (the state the reference energy names)
      step_i = 0
      term_flag = .false.
      do while (.not. term_flag .and. step_i < max_steps)
         step_i = step_i + 1
         call prop_step(sta)               ! 1. nuclear propagation
         call elec_prop(sta, dt_eff)       ! 2. electronic evolution
         call mqc_statistic(sta)           ! 3. MQC statistics
         call rec_frame(sta, u_out)        ! 4. frame recording
         term_flag = container_term(sta)   ! 5. termination check (container judgment)
      end do
      sta_fin(traj_i) = sta             ! loop tail: the final-state snapshot hand-off
      n_valid = n_valid + 1             ! placeholder semantics: no validity judgment
                                        ! exists yet - every archived trajectory counts
      call rec_final(sta, u_out)        ! the final-state row (the recorder advances its
                                        ! own trajectory counter here)
      e_fin_traj(traj_i) = rec_fin_val('e_tot')
                                        ! the drift datum: the recorder's final-state row;
                                        ! read here because the row stash lives only
                                        ! until the next rec_final
      vpot_fin_traj(traj_i) = rec_fin_val('v_pot')
                                        ! the same row's potential-energy column - the
                                        ! v_fin entry of the histogram row table
   end do

   ! ===== Phase 3: statistics (post-processing class) =====
   call cpu_time(t_cpu(2))
   allocate (cnt_outc(max_outc))
   cnt_outc = 0
   n_outc = 0
   ! the aggregation pass: ONE fin_eval per archived trajectory feeds both consumers -
   ! the guarded classification counting and the data_hist row table. Only the GENERIC
   ! columns land (the final relative energy and scattering angle from fin truth, the
   ! final potential energy from the recorder row); the remaining skeleton columns carry
   ! no generic producer in this tree and stay UNAGGREGATED - unallocated members, no
   ! guessed values
   allocate (data_hist%e_rel(n_traj), data_hist%thta(n_traj), data_hist%v_fin(n_traj))
   data_hist%n_traj = n_traj
   do traj_i = 1, n_traj
      call fin_eval(sta_fin(traj_i), fin)
      data_hist%e_rel(traj_i) = fin%e_rel     ! the classification's own energy datum,
                                             ! reused as the histogram column
      data_hist%thta(traj_i) = fin%ang(19)   ! the scattering-angle class of the angle
                                             ! array (bitwise 0 in a gas-phase run)
      data_hist%v_fin(traj_i) = vpot_fin_traj(traj_i)
                                             ! the recorder final row's potential column
      if (container_classify_bound()) then   ! bound guard: the no-classification run is
                                             ! legal and simply never calls the slot
         outc_name = container_classify(fin)
         row_outc = 0
         do k = 1, n_outc
            if (trim(outc_names(k)) == trim(outc_name)) row_outc = k
         end do
         if (row_outc == 0) then
            if (n_outc >= max_outc) then
               write (0, '(a,i0,a)') 'driver: outcome-name table full (', max_outc, &
                  ' distinct names - the statistics capacity is exhausted)'
               write (0, '(a)') 'driver: fatal (the statistics cannot proceed)'
               stop 1
            end if
            n_outc = n_outc + 1
            outc_names(n_outc) = outc_name
            row_outc = n_outc
         end if
         cnt_outc(row_outc) = cnt_outc(row_outc) + 1
      end if
   end do
   drift_max = 0.0d0                          ! data source: the recorder final-state rows
   do traj_i = 1, n_traj                      ! (e_fin_traj) against the head reference
      drift_max = max(drift_max, abs(e_fin_traj(traj_i) - e0_traj(traj_i)))
   end do
   write (*, '(a,i0,a,i0)') 'statistics: trajectories=', n_traj, ' valid=', n_valid
   if (n_outc == 0) then
      write (*, '(a)') 'outcome counts: none (no classification bound - skipped)'
   else
      write (*, '(a)', advance='no') 'outcome counts:'
      do k = 1, n_outc
         write (*, '(a,i0)', advance='no') ' '//trim(outc_names(k))//'=', cnt_outc(k)
      end do
      write (*, '(a)') ''
   end if
   write (*, '(a,i0,a,es12.4,a)') 'conservation drift: trajectories=', n_traj, &
      ' max|e_fin-e0|=', drift_max, ' kcal/mol'
   if (n_valid >= 1) then
      ! the two plain histogram kinds over the aggregated columns; range rule
      ! [0, column max + 1) keeps every aggregated row in range by construction
      ! (the v_hi > v_lo invariant is owned by the guard inside plot_hist). The
      ! lattice feeds only the lattice-aware kinds - the two plain kinds never
      ! read it
      lattice%a_lat = list_atoms%a_lat
      lattice%skew = list_atoms%skew
      lattice%q_ref = 0.0d0
      spec_hist = hist_spec_t(k_e_hist, 'hist_e_rel.dat', 0.0d0, &
                              maxval(data_hist%e_rel) + 1.0d0, 40, 0, .true.)
      call plot_hist(spec_hist, data_hist, lattice)
      spec_hist = hist_spec_t(k_t_hist, 'hist_thta.dat', 0.0d0, &
                              maxval(data_hist%thta) + 1.0d0, 40, 0, .true.)
      call plot_hist(spec_hist, data_hist, lattice)
      write (*, '(a,i0,a)') 'histograms: rows=', n_valid, &
         ' files=hist_e_rel.dat hist_thta.dat'
   end if
   write (*, '(a,f8.2,a)') 'timing: cpu=', t_cpu(2) - t_cpu(1), ' s'
   ! Not wired: dens_states / eig_out consume a mode-frequency set and an eigen-pair
   ! table - products of a Hessian / normal-mode service that does not exist in this
   ! tree, so the two calls stay out until one does; the use line above keeps the
   ! dependency face explicit
   close (u_out)
   call state_destroy(sta)
   do traj_i = 1, n_traj
      call state_destroy(sta_fin(traj_i))
   end do
   deallocate (sta_fin, e0_traj, e_fin_traj, vpot_fin_traj, cnt_outc)
end program venus
