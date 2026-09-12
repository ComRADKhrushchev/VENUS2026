!=====================================================================
! driver.f90 - L0 control kernel: program venus - three-phase orchestration (assembly -> evolution -> statistics;
!               the integration loop runs the E14 four-phase step choreography, spec §8);
!               holder of the physical state (state_t passed by argument - K1 plan D1)
! Reference: venus2026/VENUS.f90:71-115 (array zero-init + read-in - prototype of the assembly phase),
!       :157-167 (integrator parameter assembly - old site of the propagator slot), :182-232 (random-stream
!       initialization RANDST :184/:231 - the rng_init site), :239-705 (TRAJECTORY_LOOP x
!       INTEGRATION_LOOP two-loop structure - prototype of the evolution-phase skeleton), :316-349 (initial-value
!       loading and the SELECT call :349 - the samp_run site), :402-413 (initial force / first-frame GWRITE /
!       QSAVE-PSAVE staging - prototype of the trajectory-loop head), :419-465 (integration-loop INTEGRATOR
!       three-branch dispatch :423-465 - unified by prop_step), :505-576 and :581-662 (the outcome blocks of the
!       NFINAL/NTEST two exits), :669-699 (long-time exit), :243-301 (completion section:
!       outcome summary and the five histogram calls :286-299 - prototype of the statistics phase);
!       docs/superpowers/specs/2026-09-08-facade-e14-design.md §8 (driver orchestration - the
!       four-phase integration-loop body, the wiring step deleted, the E7 phase-restriction
!       restatement), §3 (assembly-phase entries: sysdef_load / container_bind_pes /
!       container_init), §6 (the sysdef step - system-definition materialization), §4 (allocation
!       timing - the state_create 1x1-seed reconcile note)
! Excision list (do not carry over when referencing; K2 plan sourcing note - system branches / duplicated
!   outcome blocks / dead blocks are all excised):
!   - RKUTTA coefficient dead block :86-94 (RKUTTA is on the K1 §2 dead-block list - never entered)
!   - checkpoint-restart branch :182-232 (NCHKP - left blank in control, introduced later if needed)
!   - system branch :322-346 (NSURF==2.NGLO surface four-atom loading / FIXROTDATM alignment /
!     THTA degree-radian fork - configuration paradigm moved into data as
!     config%reactants%surface_model; details belong to sampler/roster members)
!     (kept for context; see previous line)
!   - duplicated outcome blocks :505-576 and :581-662 (FINAL/GFINAL/record_all of the two entrances are
!     nearly line-identical - merged into a single per-trajectory pass in this skeleton's statistics phase)
!   - bounce/minimum-height tracking :453-464 (BOUNCE/z_pos_min - system diagnostics, moved to recorder
!     column registration)
!   - mode-population tracking :363-395/:468-502 (MPLOT/ENMODE/EBOND - deferred to analysis, C8)
!   - rotational-energy resampling average :511-542/:598-629 (NCROT ERAVA/ERAVB duplicated -
!     fin_t left blank / statistics phase)
!   - energy-conservation tracking and gate :451-452/:461-464/:559-565 (HTMAX-HTMIN - moved to recorder
!     energy-column aggregation + an outcome-table row)
!   - duplicated STOP/timing prints :243-301 (999/6 double channel, line-by-line repeats - merged into one statistics-section pass)
!=====================================================================
program venus
   use control,        only: task, n_traj, title, i_seed, max_steps
   use config,         only: event, observables, reactants, &
                             prop_cfg => propagator  ! alias as in propagator.f90 - the config
                             ! instance 'propagator' clashes with the module name; the dt_eff source
   use config_roster,  only: roster, roster_load
   use state,          only: state_t, state_create, state_destroy
   use rng,            only: rng_init
   use propagator,     only: prop_reg, prop_step, &
                             verlet_step, symple_step, radau_step
   use force_facade,   only: container_init, force_eval, container_term
   use elec_facade,    only: elec_prop, mqc_statistic
   use sysdef_facade,  only: sysdef_load
   use reg_adiabatic,  only: reg_adiabatic_all   ! E14 T5: the thin reg module replaces the raw
                                   ! member use line (ad_dsdt/ad_force_rule withdrawn - one
                                   ! method_reg call fills the whole 'adiabatic' package)
   use sampler,        only: samp_run
   use recorder,       only: rec_select, rec_frame, rec_final
   use final_state,    only: fin_t, fin_eval
   use outcome,        only: outc_tbl_t, outc_classify
   use density,        only: dens_states, eig_out
   use hist,           only: hist_spec_t, hist_data_t, lat_t, plot_hist
   use input,          only: read_input
   implicit none

   ! -- orchestration counters and flags (local variables - E7; the old global NTZ/NC/NX bookkeeping lands here) --
   integer :: traj_i        ! trajectory index [count] (old NTZ; 1..n_traj)
   integer :: step_i        ! step index [step] (old NC; E18: the counter of the MAX_STEPS safety
                            ! net - the driver loop condition step_i < max_steps, control side)
   logical :: term_flag     ! termination flag for this trajectory (written by the phase-4
                            ! container_term call of the integration loop - E7 termination
                            ! channel; E18: the container's own judgment, spec §8/E18)
   real(8) :: e0            ! total energy of the initial state [kcal/mol] (force_eval output -
                            ! the reference for conservation drift, old E0)
   real(8) :: dt_eff        ! effective step length [10 fs] (staged once from prop_cfg%dt - the fixed-step
                            ! source: verlet/symple take it directly, radau seeds its adaptive sequence
                            ! from it, so under radau sub-stepping the outer step length is the channel
                            ! value; consumed by elec_prop(sta, dt_eff) - the MQC hop-probability channel
                            ! g_ij(dt) needs the step's dt; sibling note: the propagator header's
                            ! endpoint-force contract, the same dt source)
   real(8) :: t_cpu(2)      ! timing pair [s] (old START_TIME/FINISH_TIME - printed in the statistics section)
   integer :: u_out         ! recording unit number [-] (shared by frames and final-state rows; old 6/8/1000+NTZ
                            ! family of unit assignments, finalized in C4)
   character(len=32) :: outc_name  ! outcome name of this trajectory [-] (statistics phase - outc_classify output)
   ! -- data structures handed between phases --
   type(state_t) :: sta                      ! physical state (held by the driver - K1 D1;
                                             ! written exclusively in the evolution phase, read-only in statistics)
   type(state_t), allocatable :: sta_fin(:)  ! final-state snapshots per trajectory [count] (evolution -> statistics
                                             ! hand-off; same source as the old QSAVE/PSAVE staging)
   type(outc_tbl_t) :: outc_tbl              ! outcome-classification table (parsed from outcome_tbl - C7; the
                                             ! expression family flagged for the L4 review, E18)
   type(fin_t) :: fin                        ! final-state quantities (taken over per trajectory in statistics)
   type(hist_spec_t) :: spec_hist            ! histogram specification (assembled in statistics - C8)
   type(hist_data_t) :: data_hist            ! trajectory row table (aggregated from recorder final-state rows - C8)
   type(lat_t) :: lattice                    ! lattice constants (mapped from the roster unit cell - statistics phase)
   integer, allocatable :: cnt_outc(:)       ! per-outcome counts [count] (aligned with outc_tbl rows;
                                             ! replaces the old NONREACT/NREACT/NLONG counters)

   ! Pseudocode (three phases - E7; spec §2.5 phase restriction table: assembly = operational class +
   !    initializing measurements, evolution = evolution class + measurement class, statistics = post-processing class):
   !
   ! ===== Phase 1: assembly =====
   ! 1. read_input('input_qct.txt') - input framework (KEYWORD=VALUE -> control/config
   !    slots / member parameters + the a_lat-class residual keys; grammar finalized in C4;
   !    old VENUS.f90:115 read_venus_input site)
   ! 2. rng_init(converted i_seed) - establish the random stream (old :184/:231 RANDST(ISEED);
   !    the control%i_seed 8-word -> i_seed conversion is finalized in C4 - input step 5)
   ! 3. sysdef_load(reactants%system_dir) - system-definition materialization (E14 spec §6: folder
   !    discovery, parser dispatch, mass lookup, roster fill, paradigm derivation, sampler hand-off;
   !    an assembly-phase data-material entry, spec §3 - no evolution-loop presence). Runs FIRST
   !    ahead of roster_load by the config_roster ordering statute (sysdef fills the roster, so
   !    roster_load shrinks to allocation/algebra; a roster_load ahead of the sysdef products is an
   !    assembly-order error). The write-hand question - whether sysdef writes roster%* directly or
   !    roster_load transfers staged products - stays C4-open
   ! 4. roster_load() - assemble and validate the roster (E14-shrunken to allocation / algebra /
   !    validation: per-field fill from the sysdef products, fragment mass sums, equilibrium-geometry
   !    placement, unit-cell derivation - the old input-keyword source is superseded, spec §9;
   !    the roster-ification of the old :316-335 initial-value loading site)
   ! 5. state_create(sta, size(roster%mass), electronic%n_surf) - allocate the physical state
   !    (natoms from the roster; n_surf passes the config value - the ontic containers a/rho are
   !    seeded at the degenerate 1x1 shape, a(1,1)=rho(1,1)=(1,0), default n_surf=1 the adiabatic
   !    degenerate point; E9)
   !    (T1 reconcile, spec §4: the 1x1 seed is the degenerate-legal DEFAULT occupation -
   !    consistent with occ=[1]; "adiabatic/LZSH/MDEF allocate neither" governs the METHOD
   !    members, which never allocate larger ontic containers - the member-init step does that
   !    (spec §4 allocation timing), at container_init)
   ! 6. member registration (the registry-family assembly-time convention; once E11 is
   !    legislated this step is replaced by the build artifact registry_gen_all(): each
   !    system-container folder carries its own reg_<system>_all module and a prebuild step
   !    aggregates them - the driver's final form has zero registration lines; effective at
   !    K3, this pseudocode is the transitional form):
   !    prop_reg(1, verlet_step), prop_reg(2, symple_step),
   !    prop_reg(3, radau_step) (the closed three-member set - X4, unifying the old :423-465 three branches);
   !    reg_adiabatic_all() - the 'adiabatic' degenerate method package (E14: the thin reg file's one
   !    method_reg('adiabatic', '', ad_prop, ad_match) call fills the whole package row - '' the
   !    empty Hamiltonian preset of the degenerate method, E17(5); the member-signature freeze
   !    lands in C2); the container's class-A binding - the system reg module's ONE
   !    container_bind_pes(container_pes [, h_deliver]) call (E17: the aggregate entry + its
   !    H delivery declaration - the member sum executes inside the container; members are
   !    system material, E12: each systems/<system>/ folder carries self-contained functional
   !    form + parameters); dist_* registration is absent by design - C5 per the T6 note
   !    (no samp_reg lines in this transitional form: the batch-3 skeleton members boltzmann/ebk/
   !    incident wait for the sampler wiring window)
   ! 7. container_init() - assembly + purely formal validation spanning BOTH facades (the
   !    single assembly entry point of the whole container, E14/E17 naming - the entry names
   !    the action; slot bound, H forms matched (h_form subset-of h_deliver, E17(5)), method
   !    slot filled via method_bind, the n_surf>=2-vs-scalar-delivery and
   !    observable-vs-container contradictions loud-stop; the POTPRE-type member-init
   !    routines run container-side in its reg module - where larger ontic containers are
   !    reallocated, spec §4 allocation timing)
   ! 8. rec_select(observables%rec_level) (or rec_select(observables%obs_list)
   !    list entry) - recording column-group selection (registration precedes column selection)
   ! 9. (E18: pred_parse withdrawn - termination is the container's own judgment, bound by
   !    its reg module's container_bind_term at step 6; no program-side grammar remains)
   !    remaining assembly parsing: event%outcome_tbl -> outc_tbl (the C7 grammar
   !    deliverable; the outcome expression family awaits the L4 review - flagged as the
   !    same complexity question); spec_hist / lattice assembly (histogram specification
   !    and roster unit-cell mapping - C8)
   ! 10. allocate(sta_fin(n_traj), cnt_outc(outc_tbl%n_row), data_hist columns -
   !     assemble the inter-phase hand-off structures; dt_eff = prop_cfg%dt (staged once - the
   !     fixed-step source, see the local's note); timing start t_cpu(1)
   !
   ! ===== Phase 2: evolution (trajectory loop x integration loop) =====
   ! [Phase restriction statement (E7 - spec §2.5, restated by E14 spec §8 with an explicit
   !    call list)]: inside the integration loop (the inner loop below) the direct calls are
   !    prop_step (propagator class), the facade phase entries elec_prop / mqc_statistic (facade
   !    class; force_eval runs inside the propagator member per the endpoint-force contract, and
   !    optionally inside mqc_statistic as an intra-box force refresh), rec_frame (recorder class),
   !    and container_term (container judgment, E18) - the four-phase orchestration is explicit in the driver.
   !    No sampler / analysis / input / roster calls inside the loop; the orchestration counters
   !    traj_i/step_i are locals.
   ! Trajectory loop: do traj_i = 1, n_traj (old TRAJECTORY_LOOP :239)
   ! 1. samp_run(scheme, sta) - per-trajectory initial-state preparation (an operational-class call - spec §2.3
   !    boundary 4: actual calls in the assembly phase include distribution members, and preparation members
   !    count as assembly-like calls; located at the trajectory-loop head, outside the integration loop;
   !    old :349 SELECT site. scheme takes the member name mapped from reactants%dist_scheme(k) - C5 wiring)
   ! 2. force_eval(sta, e0) - initial force and total energy of the initial state (old :402-403
   !    DVDQ_1/ENERGY_1; E14 naming statute - the pure force entry names the action; no dt
   !    argument)
   ! 3. rec_frame(sta, u_out) - first frame (old :405 GWRITE site)
   ! 4. step_i = 0, term_flag = .false. (orchestration reset)
   ! 5. Integration loop: do while (.not. term_flag .and. step_i < max_steps) (old
   !    INTEGRATION_LOOP :419; E18: the MAX_STEPS safety net is the driver condition - the
   !    only judgment never delegated to the container) - the spec §8
   !    four-phase step body, verbatim (the spec's dt dummy arrives here as dt_eff - the local
   !    staged from prop_cfg%dt; "occupation lag" per the landed T10 precedent - the spec's own
   !    "i_surf lag" wording is its historical citation):
   !    step_i = step_i + 1
   !    call prop_step(sta)                  ! phase 1: nuclear propagation (multiple force_eval
   !                                         ! calls inside; endpoint-force contract obliges the
   !                                         ! final ws to match the returned sta%q)
   !    call elec_prop(sta, dt_eff)          ! phase 2: electronic propagation
   !    call mqc_statistic(sta)              ! phase 3: matching (hop / momentum rescaling /
   !                                         ! optional force refresh)
   !    call rec_frame(sta, u_out)           ! frame output (post-mqc self-consistent state -
   !                                         ! repairs the prior one-step occupation lag in frames)
   !    term_flag = container_term(sta)               ! phase 4: termination (E18: the
   !                                                  ! container's own judgment; unbound slot =
   !                                                  ! .false., the safety-net-only run; end do)
   ! 6. Loop tail (outside the integration loop, at the trajectory-loop end): sta_fin(traj_i) = sta - archive
   !    the final-state snapshot (same source as old :411-413 QSAVE/PSAVE); rec_final(sta, u_out) -
   !    the final-state row (measurement class; recorder final-state rows carry system-diagnostic columns - bounce/z_min etc.,
   !    the landing spot of excision-list item 5)
   !    (the old duplicated outcome blocks :505-576/:581-662 - FINAL/GFINAL/record_all/reaction counting
   !    and the long-time exit :669-699 - all deferred to the statistics phase, see the excision list)
   !
   ! ===== Phase 3: statistics (post-processing class - outcome/analysis calls) =====
   ! 1. Per trajectory: do traj_i = 1, n_traj
   ! 1a. fin_eval(sta_fin(traj_i), fin) - final-state quantities (the generalized part of old FINAL)
   ! 1b. outc_name = outc_classify(fin, outc_tbl); increment cnt_outc for the matched row
   !     (old :552-571/:645-657 reaction counting and :669-685 long-time counting - now data-driven)
   ! 1c. data_hist aggregation (fin / recorder final-state rows -> e_rel/thta/x0/y0/z_min/
   !     v_fin/n_bnc/site columns - the assembly site of old GFINAL.f:199-202 and record_all,
   !     C8 wiring)
   ! 2. analysis calls: plot_hist(spec_hist, data_hist, lattice) (the five histogram kinds
   !    - old :286-299 five plot_* call sites; HISTO directory layout lives here); when needed,
   !    dens_states/eig_out (density of states / eigenvalue output)
   ! 3. Summary printout: outcome count table / valid-trajectory count / energy-conservation drift summary /
   !    timing (old
   !    :243-284 double channel merged into one pass; the conservation-gate verdict goes to an outcome-table row)
   ! 4. state_destroy(sta), deallocate sta_fin/cnt_outc/data_hist - wrap up
end program venus
