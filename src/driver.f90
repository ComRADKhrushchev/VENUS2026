!=====================================================================
! driver.f90 - L0 control kernel: program venus - three-phase orchestration
!              (assembly -> evolution -> statistics); the integration loop
!              runs the four-phase step choreography (E14 spec §8); holder of
!              the physical state (state_t passed by argument)
! Phase constitution legislated here (binding for the whole program):
!   - three phases, E7 phase restriction: assembly = operational calls +
!     initializing measurements; evolution = evolution + measurement calls
!     inside the loops; statistics = post-processing calls. The assembly
!     order below is binding legislation (input -> random stream -> system
!     definition -> roster -> state -> registration -> member-parameter seam
!     -> container -> recording selection -> hand-off): every consumer of a
!     datum runs after that datum's producer.
!   - E18 termination ownership: termination is the container's own judgment
!     (container_bind_term at assembly); an unbound slot reads .false. - the
!     legal safety-net-only run. MAX_STEPS is the driver loop condition -
!     the one judgment never delegated to the container.
!   - E19 classification ownership: the outcome name is the container's own
!     judgment (container_bind_classify at assembly); the statistics phase
!     calls container_classify only when container_classify_bound() says a
!     judgment is bound - the legal no-classification run guards and skips.
!   - E20 idle-chain posture: with no container dragged in, assembly
!     proceeds unbound (container_init still binds the degenerate adiabatic
!     method row); the loud-stop authority is force_eval's first call.
!   - excision-list spirit (inherited from the rewrite sourcing): no system
!     branches in the driver (the paradigm lives in data - surface_model);
!     no duplicated per-exit outcome blocks (one per-trajectory statistics
!     pass); no checkpoint-restart branch; no dead integrator-coefficient
!     block; bounce/mode-population diagnostics belong to recorder columns;
!     rotational-resampling averaging and the energy-conservation gate
!     belong to statistics (the gate verdict to the container's
!     classification).
!   - zero hand-written registration: every method package, distribution
!     member, container binding - and the integrator closed set - arrives
!     through the BUILD-GENERATED registry_gen aggregate (scripts/prebuild.sh
!     scans methods/ and systems/ for reg_*.f90; dragging a folder in is the
!     whole installation act). The program skeleton once hand-wrote the
!     three propagator registrations as a closed-family exception; that
!     exception is superseded - the integrator registrations arrive through
!     the aggregate too (methods/reg_prop.f90).
!   - the assembly seam, three-party zero knowledge: the driver knows no
!     member (it calls the registry aggregate), the input knows no member
!     (member keys stage raw), the member never touches the driver (its reg
!     file loads its own staged keys into its own init).
! Contracts staged here:
!   - size-equality guard: right after state creation the driver loud-stops unless
!     size(sta%q) == 3*size(roster%mass) - the state and the mass table
!     disagreeing is an assembly-choreography error (both sides named in the
!     message). With state_create consuming size(roster%mass) itself the
!     guard is a drift backstop, not triggerable through the public face
!     today - kept for naming the invariant it protects.
!   - dt-ordering contract: dt is finalized by read_input BEFORE the
!     member-parameter seam runs, so a member whose init reconciles against
!     dt (e.g. an equilibration step size) reads the final value. Naturally
!     satisfied by the assembly order below; legislated here so a future
!     reordering cannot silently break it.
!   - N_SURF reconcile: the input default N_SURF=0 means 'no multi-surface
!     declared' - the driver passes max(electronic%n_surf, 1) to
!     state_create (0 folds to 1: the adiabatic degenerate point, one
!     surface, the 1x1 ontic seed, E9/T1). The input key's own default
!     stays untouched.
!   - output file registration: the recording unit writes 'trajectories.txt' in the
!     run directory (cwd-relative, the semantic file name chosen for this program;
!     every trajectory's frames, gated flattened sections and final-state rows
!     land there; the statistics summary lines go to stdout). The name is
!     registered here - criteria and tooling grep this contract.
! Units: none (orchestration; quantities carry their own units).
!=====================================================================
program venus
   use control,        only: n_traj, i_seed, max_steps
   use config,         only: observables, reactants, electronic, &
                             prop_cfg => propagator  ! alias as in propagator.f90 - the config
                             ! instance 'propagator' clashes with the module name; the dt_eff source
   use config_roster,  only: roster, roster_load
   use state,          only: state_t, state_create, state_destroy
   use rng,            only: rng_init
   use propagator,     only: prop_step
   use force_facade,   only: container_init, force_eval, container_term
   use elec_facade,    only: elec_prop, mqc_statistic
   use sysdef_facade,  only: sysdef_load
   use registry_gen,   only: registry_gen_all, registry_gen_init   ! the BUILD-GENERATED
                                   ! aggregate (scripts/prebuild.sh emits it from every
                                   ! reg_*.f90 under methods/ and systems/ - drag-in; zero
                                   ! hand-written registration lines in the driver)
   use sampler,        only: samp_run, samp_n_dist
   use recorder,       only: rec_select, rec_frame, rec_final, rec_n_sel, rec_fin_val
   use final_state,    only: fin_t, fin_eval, container_classify, container_classify_bound
   use density,        only: dens_states, eig_out
   use hist,           only: hist_spec_t, hist_data_t, lat_t, plot_hist
   use input,          only: read_input
   implicit none

   ! -- orchestration counters and flags (locals - E7) --
   integer :: traj_i        ! trajectory index [count] (1..n_traj)
   integer :: step_i        ! step index [step] (E18: the counter of the MAX_STEPS safety
                            ! net - the driver loop condition step_i < max_steps)
   logical :: term_flag     ! termination flag for this trajectory (written by the phase-4
                            ! container_term call of the integration loop; E18: the
                            ! container's own judgment)
   real(8) :: e0            ! total energy of the initial state [kcal/mol] (force_eval output -
                            ! the reference for conservation drift)
   real(8) :: dt_eff        ! effective step length [10 fs] (staged once from prop_cfg%dt - the
                            ! fixed-step source: verlet/symple take it directly, radau seeds its
                            ! adaptive sequence from it, so under radau sub-stepping the outer
                            ! step length is the channel value; consumed by elec_prop(sta, dt_eff)
                            ! - the hop-probability channel g_ij(dt) needs the step's dt)
   real(8) :: t_cpu(2)      ! timing pair [s] (printed in the statistics section)
   integer :: u_out         ! recording unit number [-] (shared by frames and final-state rows;
                            ! opened on the registered output file 'trajectories.txt' at the
                            ! assembly tail, closed at the statistics wrap-up)
   character(len=32) :: outc_name  ! outcome name of this trajectory [-] (statistics phase -
                                   ! the container's classification output, E19)
   integer, parameter :: max_outc = 32   ! outcome-name table capacity [name]
   character(len=32) :: outc_names(max_outc) ! outcome-name table [-] (name-keyed counting:
                                             ! first-seen insertion order, counts parallel in
                                             ! cnt_outc - the shape fixed with the first real
                                             ! container)
   integer :: n_outc        ! distinct outcome names counted so far [name]
   integer :: n_valid       ! valid-trajectory count [count] (a final row archived)
   integer :: k             ! scratch index [-]
   integer :: row_outc      ! outcome-table row scratch [-]
   real(8), allocatable :: e0_traj(:)     ! per-trajectory reference energy [kcal/mol] (the
                                          ! head force_eval output - the drift reference)
   real(8), allocatable :: e_fin_traj(:)  ! per-trajectory final e_tot [kcal/mol] (the
                                          ! recorder final-state row - see the statistics note)
   real(8) :: drift_max     ! max |e_fin - e0| over trajectories [kcal/mol]
   ! -- data structures handed between phases --
   type(state_t) :: sta                      ! physical state (held by the driver; written
                                             ! exclusively in the evolution phase, read-only in
                                             ! statistics)
   type(state_t), allocatable :: sta_fin(:)  ! final-state snapshots per trajectory [count]
                                             ! (evolution -> statistics hand-off)
   type(fin_t) :: fin                        ! final-state quantities (taken over per
                                             ! trajectory in statistics; the argument type of
                                             ! the container-classification slot, E19)
   integer, allocatable :: cnt_outc(:)       ! per-outcome counts [count] (name-keyed over the
                                             ! names container_classify returns, E19; the
                                             ! counting shape - name-keyed vs the container
                                             ! declaring its outcome-name list - decided with
                                             ! the first real container)
   type(hist_spec_t) :: spec_hist            ! histogram specification (assembled in statistics)
   type(hist_data_t) :: data_hist            ! trajectory row table (aggregated from recorder
                                             ! final-state rows)
   type(lat_t) :: lattice                    ! lattice constants (mapped from the roster unit
                                             ! cell - statistics phase)

   ! ===== Phase 1: assembly (real - the order below is binding legislation) =====
   call read_input('input_qct.txt')                ! 1. input framework (KEYWORD=VALUE ->
                                                   !    control/config slots; member keys stage
                                                   !    raw for the seam at step 7)
   call rng_init(i_seed(1))                        ! 2. random stream (the ISEED scalar statute:
                                                   !    word 1 seeds the stream)
   call sysdef_load(reactants%system_dir)          ! 3. system definition (folder scan -> parser
                                                   !    dispatch -> mass table; a bad folder
                                                   !    loud-stops inside sysdef_load - the
                                                   !    driver adds no judgment of its own)
   call roster_load()                              ! 4. roster (fragments / masses / equilibrium
                                                   !    placement, from the sysdef products
                                                   !    through the staging accessors)
   call state_create(sta, size(roster%mass), max(electronic%n_surf, 1))
                                                   ! 5. physical state (natoms from the roster;
                                                   !    the max() is the N_SURF reconcile - the
                                                   !    input default 0 folds to the adiabatic
                                                   !    degenerate point, 1 surface, E9/T1 1x1
                                                   !    seed)
   if (size(sta%q) /= 3*size(roster%mass)) then    ! size-equality guard (drift backstop; both sides
                                                   ! named)
      write (0, '(a,i0,a,i0,a)') 'driver: size guard - the state holds ', size(sta%q), &
         ' coordinates but the roster mass table holds ', 3*size(roster%mass), &
         ' (assembly-choreography error, formal)'
      write (0, '(a)') 'driver: loud-stop (the state and the roster disagree)'
      stop 1
   end if
   call registry_gen_all()                         ! 6. registration aggregate (method packages,
                                                   !    distribution members, container bindings,
                                                   !    and the integrator closed set - no
                                                   !    hand-written registration lines here)
   call registry_gen_init()                        ! 7. the assembly seam: staged member keys ->
                                                   !    each member's own init (an unselected
                                                   !    member skips inside its own guard).
                                                   !    dt-ordering: dt was finalized at step
                                                   !    1, so a member reconciling against dt
                                                   !    reads the final value - the assembly
                                                   !    order is the contract
   call container_init()                           ! 8. container assembly + formal validation
                                                   !    (E20: legal with no container dragged in -
                                                   !    the degenerate method row still binds;
                                                   !    loud-stop authority = force_eval's first
                                                   !    call)
   call rec_select(observables%rec_level)          ! 9a. recording preset group (level entry)
   if (len_trim(observables%obs_list) > 0) then    ! 9b. named-column list entry (union
      call rec_select(trim(observables%obs_list))  !    semantics with the level entry - the
   end if                                          !    recorded set is the union of both)
   allocate (sta_fin(n_traj), e0_traj(n_traj), e_fin_traj(n_traj))
                                                   ! 10. inter-phase hand-off structures (the
                                                   !     snapshots + the drift reference/final
                                                   !     pairs the statistics pass consumes)
   dt_eff = prop_cfg%dt                            !     staged once (the fixed-step source -
                                                   !     see the local's note)
   call cpu_time(t_cpu(1))                         !     timing start
   open (newunit=u_out, file='trajectories.txt', status='replace', action='write')
                                                   !     the recording unit goes live on the
                                                   !     registered output file (cwd-relative,
                                                   !     the semantic name stated in the module
                                                   !     header; frames, flattened sections and
                                                   !     final-state rows all land here)
   write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
      'assembly complete: trajectories=', n_traj, ' atoms=', size(roster%mass), &
      ' fragments=', size(roster%frag), ' paradigm=', reactants%surface_model, &
      ' members=', samp_n_dist(), ' init=ok columns=', rec_n_sel()

   ! ===== Phase 2: evolution (trajectory loop x integration loop) =====
   ! Phase restriction (E7): inside the integration loop the direct entries are
   ! exactly the four-phase choreography - prop_step, elec_prop, mqc_statistic,
   ! rec_frame, container_term (the body below); force_eval runs inside the
   ! propagator member per the endpoint-force contract, and the sampler is a
   ! trajectory-loop-head operational entry, never an integration-loop entry.
   n_valid = 0
   do traj_i = 1, n_traj
      sta%t = 0.0d0                     ! fresh-trajectory reset: the time origin and the
                                        ! fresh-trajectory staging key (a propagated state's
                                        ! boot tests read t == 0; one reset before the sampler
                                        ! suffices - the thermalize member saves and restores
                                        ! t around its own equilibration)
      call samp_run('incident', sta)    ! per-trajectory initial-state preparation (the
                                        ! per-fragment schemes arrive inside samp_run through
                                        ! the dist_scheme data; 'incident' is the one
                                        ! orchestration member not selected per-fragment)
      call force_eval(sta, e0)          ! initial force + the reference energy of the
                                        ! initial state (the pure force entry; no dt
                                        ! argument)
      e0_traj(traj_i) = e0
      call rec_frame(sta, u_out)        ! first frame (the state the reference energy
                                        ! names)
      step_i = 0
      term_flag = .false.
      do while (.not. term_flag .and. step_i < max_steps)
         step_i = step_i + 1
         call prop_step(sta)
         call elec_prop(sta, dt_eff)
         call mqc_statistic(sta)
         call rec_frame(sta, u_out)
         term_flag = container_term(sta)
      end do
      sta_fin(traj_i) = sta             ! loop tail: the final-state snapshot hand-off
      n_valid = n_valid + 1
      call rec_final(sta, u_out)        ! the final-state row (the recorder advances its
                                        ! own trajectory counter here)
      e_fin_traj(traj_i) = rec_fin_val('e_tot')
                                        ! the drift datum: the recorder's final-state row
                                        ! (a force_eval reading at the final geometry);
                                        ! read here because the row stash lives only until
                                        ! the next rec_final - the data-source choice is
                                        ! stated at the drift summary below
   end do

   ! ===== Phase 3: statistics (post-processing class) - first version =====
   call cpu_time(t_cpu(2))
   allocate (cnt_outc(max_outc))
   cnt_outc = 0
   n_outc = 0
   do traj_i = 1, n_traj
      call fin_eval(sta_fin(traj_i), fin)
      if (container_classify_bound()) then   ! the bound guard: the no-classification run
                                             ! is legal and simply never calls the slot
         outc_name = container_classify(fin)
         row_outc = 0
         do k = 1, n_outc
            if (trim(outc_names(k)) == trim(outc_name)) row_outc = k
         end do
         if (row_outc == 0) then
            if (n_outc >= max_outc) then
               write (0, '(a,i0,a)') 'driver: outcome-name table full (', max_outc, &
                  ' distinct names - the statistics capacity is exhausted)'
               write (0, '(a)') 'driver: loud-stop (the statistics cannot proceed)'
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
   write (*, '(a,f8.2,a)') 'timing: cpu=', t_cpu(2) - t_cpu(1), ' s'
   ! Analysis call sites, prepared but not wired (the histogram machinery is a
   ! later landing; the locals above keep the signatures honest):
   ! call plot_hist(spec_hist, data_hist, lattice)
   ! call dens_states(spec_hist)
   ! call eig_out(spec_hist)
   close (u_out)
   call state_destroy(sta)
   do traj_i = 1, n_traj
      call state_destroy(sta_fin(traj_i))
   end do
   deallocate (sta_fin, e0_traj, e_fin_traj, cnt_outc)
end program venus
