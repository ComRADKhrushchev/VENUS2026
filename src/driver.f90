!=====================================================================
! driver.f90 - L0 control kernel: program venus - three-phase orchestration (assembly -> evolution -> statistics);
!               holder of the physical state (state_t passed by argument - K1 plan D1)
! Reference: venus2026/VENUS.f90:71-115 (array zero-init + read-in - prototype of the assembly phase),
!       :157-167 (integrator parameter assembly - old site of the propagator slot), :182-232 (random-stream
!       initialization RANDST :184/:231 - the rng_init site), :239-705 (TRAJECTORY_LOOP x
!       INTEGRATION_LOOP two-loop structure - prototype of the evolution-phase skeleton), :316-349 (initial-value
!       loading and the SELECT call :349 - the samp_run site), :402-413 (initial force / first-frame GWRITE /
!       QSAVE-PSAVE staging - prototype of the trajectory-loop head), :419-465 (integration-loop INTEGRATOR
!       three-branch dispatch :423-465 - unified by prop_step), :505-576 and :581-662 (the outcome blocks of the
!       NFINAL/NTEST two exits), :669-699 (long-time exit), :243-301 (completion section:
!       outcome summary and the five histogram calls :286-299 - prototype of the statistics phase)
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
   use control,        only: task, n_traj, title, i_seed
   use config,         only: event, observables, reactants
   use config_roster,  only: roster, roster_load
   use state,          only: state_t, state_create, state_destroy
   use rng,            only: rng_init
   use predicate,      only: expr_t, pred_parse, pred_eval
   use propagator,     only: prop_reg, prop_step, &
                             verlet_step, symple_step, radau_step
   use force_facade,   only: pes_reg, blackbox_init, force_eval
   use reg_adiabatic,  only: reg_adiabatic_all   ! E14 T5: the thin reg module replaces the raw
                                   ! member use line (ad_dsdt/ad_force_rule withdrawn - one
                                   ! method_reg call fills the whole 'adiabatic' package)
   use sampler,        only: samp_reg, samp_run
   use recorder,       only: rec_select, rec_frame, rec_final
   use final_state,    only: fin_t, fin_eval
   use outcome,        only: outc_tbl_t, outc_classify
   use density,        only: dens_states, eig_out
   use hist,           only: hist_spec_t, hist_data_t, lat_t, plot_hist
   use input,          only: read_input
   implicit none

   ! -- orchestration counters and flags (local variables - E7; the old global NTZ/NC/NX bookkeeping lands here) --
   integer :: traj_i        ! trajectory index [count] (old NTZ; 1..n_traj)
   integer :: step_i        ! step index [step] (old NC - the MAX_STEPS channel of pred_eval)
   logical :: term_flag     ! termination flag for this trajectory (held by the hooks%term closure host - E7
                            ! evolution-loop termination channel, same source as propagator.f90 note 3)
   real(8) :: e0            ! total energy of the initial state [kcal/mol] (facade_force output -
                            ! the reference for conservation drift, old E0)
   real(8) :: t_cpu(2)      ! timing pair [s] (old START_TIME/FINISH_TIME - printed in the statistics section)
   integer :: u_out         ! recording unit number [-] (shared by frames and final-state rows; old 6/8/1000+NTZ
                            ! family of unit assignments, finalized in C4)
   character(len=32) :: outc_name  ! outcome name of this trajectory [-] (statistics phase - outc_classify output)
   ! -- data structures handed between phases --
   type(state_t) :: sta                      ! physical state (held by the driver - K1 D1;
                                             ! written exclusively in the evolution phase, read-only in statistics)
   type(state_t), allocatable :: sta_fin(:)  ! final-state snapshots per trajectory [count] (evolution -> statistics
                                             ! hand-off; same source as the old QSAVE/PSAVE staging)
   type(expr_t)  :: expr_term                ! termination-predicate tree (produced by pred_parse in the assembly phase)
   type(outc_tbl_t) :: outc_tbl              ! outcome-classification table (parsed from outcome_tbl - C7)
   type(fin_t) :: fin                        ! final-state quantities (taken over per trajectory in statistics)
   type(hist_spec_t) :: spec_hist            ! histogram specification (assembled in statistics - C8)
   type(hist_data_t) :: data_hist            ! trajectory row table (aggregated from recorder final-state rows - C8)
   type(lat_t) :: lattice                    ! lattice constants (mapped from the roster unit cell - statistics phase)
   integer, allocatable :: cnt_outc(:)       ! per-outcome counts [count] (aligned with outc_tbl rows;
                                             ! replaces the old NONREACT/NREACT/NLONG counters)

   ! Pseudocode (three phases - E7; spec §2.5 phase restriction table: assembly = operational class +
   ! initializing measurements, evolution = evolution class + measurement class, statistics = post-processing class):
   !
   ! ===== Phase 1: assembly =====
   ! 1. read_input('input_qct.txt') - input framework (KEYWORD=VALUE -> control/config
   !    slots / roster staging / member parameters; grammar finalized in C4; old VENUS.f90:115 read_venus_input site)
   ! 2. rng_init(converted i_seed) - establish the random stream (old :184/:231 RANDST(ISEED);
   !    the control%i_seed 8-word -> i_seed conversion is finalized in C4 - input step 5)
   ! 3. roster_load() - assemble the roster (mass table / fragment-A table / equilibrium geometry / unit cell;
   !    the roster-ification of the old :316-335 initial-value loading site)
   ! 4. state_create(sta, size(roster%mass), electronic%n_surf) - allocate the physical state
   !    (natoms from the roster; n_surf passes the config value - electronic amplitudes are sized by the
   !    number of surfaces, default 1, the adiabatic degenerate point; E9)
   !    (see previous line)
   ! 5. member registration (one line each - the registry-family assembly-time convention; once
   !    E11 is legislated this step is replaced by the build artifact registry_gen_all(): each
   !    system-blackbox folder carries its own reg_<system>_all module and a prebuild step
   !    aggregates them - the driver's final form has zero registration lines; effective at
   !    K3, this pseudocode is the transitional form):
   !    prop_reg(1, verlet_step), prop_reg(2, symple_step),
   !    prop_reg(3, radau_step) (the closed three-member set - X4, unifying the old :423-465 three branches);
   !    reg_adiabatic_all() - the 'adiabatic' degenerate method package (E14: the thin reg file's one
   !    method_reg(ad_prop, ad_match) call fills the whole package row - the member-signature freeze
   !    lands in C2); pes_reg(main-PES name, ...) - pes_*
   !    members are system material (E12: each systems/<system>/ folder carries self-contained functional
   !    form + parameters; single-surface members register with the default r_scalar, multi-surface main
   !    PES with optional repr=r_diab/r_adiab - the E10 unified single signature, still one registration
   !    line per member, T2 unchanged); samp_reg(distribution member name, ...) - dist_* member wiring in C5
   !    (batch-3 skeleton boltzmann/ebk/incident + the remaining 9)
   ! 6. facade_init() - member assembly + capability-triple validation (illegal combinations loud-stop)
   ! 7. rec_select(observables%rec_level) (or rec_select(observables%obs_list)
   !    list entry) - recording column-group selection (registration precedes column selection)
   ! 8. pred_parse(event%term_expr, expr_term) - assemble the termination tree; with no TERMINATE key,
   !    the default expression (same-parameter version - built by the assembly layer, predicate.f90 note 5);
   !    event%outcome_tbl -> outc_tbl parsing (the C7 grammar deliverable); spec_hist /
   !    lattice assembly (histogram specification and roster unit-cell mapping - C8)
   ! 9. hooks wiring (the E7 termination channel - propagator.f90 note 3): hooks%out => output
   !    closure (host internally calls rec_frame(sta, u_out) - nopass pointer + unit number held by the host;
   !    closure-host routine shape finalized in C3); hooks%term => termination closure (internally calls
   !    pred_eval(sta, expr_term, step_i) and sets term_flag when true - the flag is held by this program)
   ! 10. allocate(sta_fin(n_traj), cnt_outc(outc_tbl%n_row), data_hist columns -
   !     assemble the inter-phase hand-off structures; timing start t_cpu(1)
   !
   ! ===== Phase 2: evolution (trajectory loop x integration loop) =====
   ! [Phase restriction statement (E7 - spec §2.5)]: inside the integration loop (the inner loop below)
   !    only the four call classes facade/propagator/recorder/predicate appear - the only direct call is
   !    prop_step (propagator class); facade_force (facade class) is evaluated inside the propagator member,
   !    rec_frame (recorder class) comes through the hooks%out checkpoint callback,
   !    pred_eval (predicate class) through the hooks%term closure. No sampler /
   !    analysis/input/roster calls inside the loop; the orchestration counters traj_i/step_i are locals.
   ! Trajectory loop: do traj_i = 1, n_traj (old TRAJECTORY_LOOP :239)
   ! 1. samp_run(scheme, sta) - per-trajectory initial-state preparation (an operational-class call - spec §2.3
   !    boundary 4: actual calls in the assembly phase include distribution members, and preparation members
   !    count as assembly-like calls; located at the trajectory-loop head, outside the integration loop;
   !    old :349 SELECT site. scheme takes the member name mapped from reactants%dist_scheme - C5 wiring)
   ! 2. facade_force(sta, e0) - initial force and total energy of the initial state (old :402-403 DVDQ_1/ENERGY_1)
   ! 3. rec_frame(sta, u_out) - first frame (old :405 GWRITE site)
   ! 4. step_i = 0, term_flag = .false. (orchestration reset)
   ! 5. Integration loop: do while (.not. term_flag) (old INTEGRATION_LOOP :419)
   ! 5a. step_i = step_i + 1 (orchestration counter - the MAX_STEPS safety-net channel)
   ! 5b. prop_step(sta, hooks) - advance one checkpoint (inside the member: facade_force called
   !     multiple times + checkpoint callbacks, first out (rec_frame) then term (pred_eval -> term_flag);
   !     the unified site of the old :423-465 three branches + :450/:459 interval output)
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
