# E14 Design: Four-Entry Blackbox Contract — facade Rebuild

- Date: 2026-09-08
- Status: proposed (awaiting owner review)
- Scope: delete and rebuild the facade layer per the blackbox contract requirements;
  promote interface_electronic/ + interface_distributions/ to a flat built-in method library.

## Ruling revision ledger (what this supersedes)

| Prior ruling | Revision |
|---|---|
| E13 "facade_force is the orchestration of the three" | Withdrawn. The blackbox contract becomes four entries; orchestration belongs to the driver. |
| E10 axiom B "outside the box only (Q,P,s[,dt])->(F,E) is visible" | Outside the box four signatures are visible: `force_eval` / `elec_prop` / `mqc_statistic` / termination primitives (composed through predicate). ws still never leaves the box. |
| E10 axiom C "ws lifetime = one evaluation, declared and dying inside facade_force" | ws lifetime extends to one integration step: filled at force evaluation, consumed by elec_prop / mqc_statistic, reset at step end. Holder: a module-private instance in elec_facade (still invisible outside the box). |
| E9 "class-B dt channel = third argument of facade_force" | dt is removed from the facade_force signature; it becomes the argument of `elec_prop(sta, dt)` — hop probability g_ij(dt) belongs to the electronic phases. |
| E7 "inside the integration loop the only direct call is prop_step" | Revised: the loop body directly calls prop_step and the facade phase entries, plus rec_frame/pred_eval — the four-phase orchestration is explicit in the driver. |
| N2 "the predicate primitive set is closed (four primitives)" | The primitive set becomes an open registered set (system termination primitives as grammar leaves); the boolean combinators remain a closed set. |
| E12 "interface_electronic/ remains the class-B home" | Role upgraded: the two interface_* directories are rebuilt as the flat built-in method library `methods/` (pre-stored generic methods, library-level registration, config-key dispatch). |
| E11 "prebuild scans systems/*/reg_*.f90" | Scan rule generalized: prebuild aggregates every reg_*.f90 under methods/ and systems/ into registry_gen.f90. The reg file is the atomic unit of registration (drag in a new file = drag in a new method or system; zero edits to existing files). |
| E12 "state.f90 n_surf/i_surf representation semantics" | i_surf is ruled a statistical variable, not an electronic one (the discrete output of mqc's stochastic sampling of the populations — an irreversible classical-side projection). i_surf and n_hop relocate from elec_t to state_t top level (classical-side statistical fields); elec_t becomes purely the quantum-state ontic layer {n_surf, a, rho}. |
| (naming) facade_-prefixed entry names | Withdrawn. Entry names state the physical role: force_eval / elec_prop / mqc_statistic (evolution), blackbox_init / blackbox_info / sysdef_load (assembly), pes_reg / method_reg / pred_reg_prim (registration, noun_reg per the prop_reg/samp_reg convention). Module names keep *_facade (architectural position). |

## 1. Background and motivation

The prior facade (force_facade.f90, E13 form) exposes a single entry `facade_force(sta, e_tot, dt)` that internally orchestrates force composition, unit conversion, and a monolithic class-B call `elec_proc(s, f, dt, ws)`. Against the three-layer decomposition of MQC propagation (design_spec §4.1: nuclear propagation / electronic propagation / state-statistics matching):

1. layers 2 and 3 are fused in one interface — the compositionality of Corollary 1 (Ehrenfest and FSSH share all electronic-propagation code; a method = force rule × matching rule) cannot be realized;
2. the matching layer needs the electronic statistics in ws (d_ij for momentum projection), but ws dies when facade_force returns — an independent matching entry has nothing to read;
3. the termination criterion lives entirely outside the blackbox contract (predicate + hooks%term);
4. the electronic state container (complex vector amp) cannot naturally host multi-electron methods (IESH electron×state matrices) or density-matrix-native methods.

E14 therefore rebuilds the facade layer as a four-entry contract with driver-side orchestration.

## 2. Directory end-state

```
src/
  facade/
    force_facade.f90    REWRITTEN  — force entry + class-A pes registry + unit contract
    elec_facade.f90     NEW        — elec_prop / mqc_statistic entries + method-package
                                     registry + ws holder
    sysdef_facade.f90   NEW        — system-definition-file contract: discovery, parsing dispatch,
                                     roster fill, paradigm derivation, hand-off to sampler
    elec_ws.f90         REVISED    — ws_t container contract (comments: one-step lifetime, holder)
    sampler.f90         KEPT       — interface-B' facade (initial-state distributions)
  methods/              REBUILT    — flat built-in method library (from interface_electronic/
                                     + interface_distributions/), file prefixes self-describe the
                                     family: elec_adiabatic.f90, elec_fssh.f90 (phase C),
                                     dist_boltzmann.f90, dist_ebk.f90, dist_incident.f90, ...
                                     shared amplitude integrators live here as library modules
                                     (amp_rk.f90 / amp_expm.f90, phase C)
  systems/<system>/     unchanged   — each folder: PES members + system-definition files
                                     (one per reactant) + reg_<system>_all.f90
  dynamics/             unchanged   — propagator.f90 (closed family, hooks withdrawn) + predicate.f90
                                     (grammar interpreter + primitive registry)
```

Termination-primitive registration lives in predicate.f90 (leaf-node types and evaluation are already there — cohesive), not in facade/.

## 3. The four-entry contract

Naming statute: module names keep the *_facade suffix (architectural position — residents of
facade/); entry names state the physical role (force_eval / elec_prop / mqc_statistic). The
module names the place, the entry names the action.

```fortran
! ---- force_facade.f90 — entry 1: force phase (force supply for nuclear propagation) ----
subroutine force_eval(sta, e_tot)
   ! (Q,P,s) -> (F,E): class-A members deliver via ws -> facade composes -> converts (the
   ! unique unit-contract point, ev_kcal) -> writes sta%f, e_tot. As a side effect the ws of the
   ! current geometry is left in place for the later electronic phases (one-step lifetime).
   ! The dt argument is removed (E9 channel withdrawn — dt now belongs to elec_prop).

! ---- elec_facade.f90 — entry 2: electronic propagation (§4.1 layer 2) ----
subroutine elec_prop(sta, dt)
   ! Reads the held ws + sta%q; writes the propagation container (sta%s%a or sta%s%rho,
   ! whichever the method declares as its ontic container). Forms: amplitude ODE
   ! (TSH/IESH/Ehrenfest), mapping flow (SQC), analytic elimination (MDEF: no-op),
   ! no-op (LZSH, adiabatic). Pure-form tail call: sync_proc if registered.

! ---- elec_facade.f90 — entry 3: state-statistics matching (§4.1 layer 3) ----
subroutine mqc_statistic(sta)
   ! Input: the electronic state (+ ws statistics). Output: the state — writes sta%i_surf
   ! (discrete reassignment — a classical-side statistical field, see §4), sta%p (momentum
   ! projection), optionally sta%f (post-hop force refresh — the member decides and may call
   ! force_eval itself; intra-box calls are legal), optionally the propagation container
   ! (decay-of-mixing amplitude resets).
   ! The sole consumption point of the random stream during the evolution phase (§4.1 ruling).
   ! Pure-form tail call: sync_proc if registered.

! ---- entry 4: termination — two-track composition ----
! predicate.f90: pred_reg_prim(name, proc) registers a system primitive
!   (logical function proc(sta)); grammar leaf PRIM('name') participates in boolean
!   composition on equal footing with COM_DISTANCE/Z_HEIGHT/BOUND/MAX_STEPS.
! The built-in grammar (four primitives + AND/OR/NOT, closed combinator set) stays;
! the leaf set opens up. MAX_STEPS remains the ever-present safety net.
```

### Assembly-phase entries (not part of the evolution loop)

```fortran
! force_facade.f90: pes_reg(name, proc [, repr])   — class-A registration (E10 single
!                    signature, unchanged; noun_reg per the prop_reg/samp_reg convention);
!                    blackbox_init() — assembly + purely formal validation, now spanning both
!                    facades (names resolve, main-PES slot filled, method slot filled, additive
!                    members r_scalar, n_surf>=2 vs r_scalar contradiction, observable-vs-
!                    container delivery contradiction);
!                    blackbox_info() — diagnostics (member lists, assembly state).
! elec_facade.f90:  method_reg(name, prop_proc, match_proc [, sync_proc]) — the method
!                    package registration (one call fills the whole package).
! predicate.f90:    pred_reg_prim(name, proc) — termination primitives.
! sysdef_facade.f90: sysdef_load() — consumed by the driver assembly phase (discovery, parsing
!                    dispatch, roster fill, paradigm derivation, sampler hand-off); no
!                    evolution-loop presence.
```

### ws holding and the endpoint-force contract

- The ws instance is a module-private variable of elec_facade (saved, not a dummy); it never
  appears in any public signature. force_facade fills it (use elec_facade to obtain the
  instance); elec_prop / mqc_statistic consume it; reset to r_empty at each step end.
- Dependency direction: force_facade → elec_facade → elec_ws; members (pes_*/elec_*/dist_*)
  → elec_ws (type visibility only).
- Propagator member obligation (legislated in propagator.f90 header): prop_step returns only
  after the endpoint force has been evaluated (ws corresponds to the returned sta%q). Under MQC
  the velocity-Verlet right-end force cannot be folded into the next step's left end — the
  amplitude propagation of phase 2 needs the endpoint H.

## 4. state_t revision: symmetric dual electronic containers + statistical fields on the classical side

```fortran
type :: elec_t
   integer  :: n_surf = 1                  ! unchanged (E12 semantics statute retained)
   complex(8), allocatable :: a(:,:)       ! UPGRADED from complex vector amp(:): coefficient
                                           ! matrix. Row/column index semantics defined by the
                                           ! class-B method member (E11/E12 route: surface number /
                                           ! electron×state / configuration cardinality / mapping
                                           ! z = x + ip encoding).
   complex(8), allocatable :: rho(:,:)     ! NEW: density matrix.
end type
! elec_t is now purely the quantum-state ontic layer (state + dimension label).

type :: state_t
   ! ... q/p/f/t unchanged ...
   type(elec_t) :: s                       ! electronic state (ontic: a / rho)
   integer :: i_surf = 1                   ! RELOCATED from elec_t: active-surface index — a
                                           ! CLASSICAL-SIDE STATISTICAL field, not an electronic
                                           ! variable: the discrete output of mqc_statistic's
                                           ! stochastic sampling of the populations (irreversible
                                           ! projection of quantum statistics onto the classical
                                           ! side). Read by force_eval (adiab-mode i_want
                                           ! request source); written by mqc_statistic.
   integer :: n_hop  = 0                   ! RELOCATED from elec_t: hop counter — same statistical
                                           ! nature (pure observable; read by recorder only).
end type
```

Conceptual gain: "mqc_statistic takes the electronic state as input and returns the state"
now reads literally true — its outputs are classical-side fields of state_t (i_surf, p),
while the electronic ontic state (a/rho) is touched only when a method lawfully resets it
(decay-of-mixing).

- The two containers are **symmetric ontic candidates**: a method declares which one it
  propagates (amplitude ODE methods work on a; density-matrix-dynamics methods work on rho;
  adiabatic/LZSH/MDEF allocate neither). state legislates only the existence of two
  allocatable complex matrices — no subordination between them.
- Representation covariance ≠ representation construct: a and rho are the electronic state
  itself and its observable face (trajectory lifetime, legitimately in state_t); H_ij/d_ij
  remain ws constructs (one-step lifetime).
- Synchronization contract (digesting the dual-container risk with a formal guarantee): the
  method-package third slot `sync_proc(s)` — "bring the containers the method declares as
  deliverables up to date". Refresh direction and form (a→rho, rho→a, or none) are entirely
  the method's knowledge; facade calls it as a pure-form tail of elec_prop / mqc_statistic
  when registered.
  The assembly-time check is purely formal: an observable column demanding a container that no
  assembled package member delivers -> loud-stop naming the conflict (derivable from the
  container contract alone).
- Allocation timing: state_create allocates the degenerate shape (single surface / 1×1);
  a method package whose ontic container is larger (e.g. IESH n_e × n_s) reallocates in its
  member-init hook (the blackbox_init step — POTPRE-type, per the existing assembly spec).
- recorder reads rho diagonals (populations) / off-diagonals (coherences) as observable
  columns — never touching the semantics of a.

## 5. The built-in method library (methods/)

Architecture rationale:

1. Ownership: MQC methods are parameter-free generic algorithms with zero system knowledge —
   the generic side is their lawful home (E12 assigned system material to systems/; methods
   were always generic).
2. Corollary 1 realized in-library: shared amplitude-integrator modules (RK / matrix
   exponential) are library modules; FSSH and Ehrenfest share the same prop member and differ
   only in match — the method count drops multiplicatively.
3. Testability: library methods pair with the Tully benchmark series (R8-5) for per-layer
   acceptance — far more testable than copies scattered across systems.
4. Blackbox authoring drops to configuration level: NO+Au IESH = the pre-stored 'iesh' package
   + the system PES. The system author writes PES members and system-definition files;
   the method is a config key.

### Degree of system-independence of a mature MQC method (FSSH worked example)

Per-component analysis of FSSH: amplitude propagation (dc/dt = -iHc along the trajectory —
a matrix ODE integration) — 100% generic; hop probability g_ij = 2·dt·Re[c_i* c_j (d_ij·v)]/|c_j|^2
— 100% generic formula; hop acceptance + momentum-projection energy budget — 100% generic;
momentum rescaling p' = p - beta·d_ij — 100% generic; frustrated-hop policy and decoherence
corrections — generic formulas with method parameters (sub-options, not system knowledge).
All system dependence is pushed to the input side and isolated by three contractual layers:
(1) the ws container — the sole channel through which the system delivers H/d_ij (the method
imports no system module, it only eats the container); (2) the representation-support
declaration — standard FSSH eats r_adiab (v_surf/g_surf/d_ij); diabatic input goes through the
facade diagonalization view (branch 4'a ruling); unsupported repr -> loud-stop (E10 contract);
(3) method parameters ride config. Tully 1990 itself is presented this way — models (systems)
and algorithms (methods) separated, the algorithm written purely in H_ij/d_ij notation; the
R8-5 benchmark series validates the same algorithm against different models, per layer.
IESH is isomorphic (independent-electron n_e × n_s matrix algebra, still container-only).
Honest boundary: numerical-stability knowledge (amplitude sub-stepping policy etc.) is
generic-side-lawful; system authors' parameter tuning rides config and does not break the
single-copy property of the algorithm.

Registration data flow (two levels):

```
Library level:  each method file is registered through its own thin reg_*.f90 in methods/
                (e.g. reg_fssh.f90 calling method_reg('fssh', prop, match [, sync])).
                Built-in set: 'adiabatic' (all no-ops — the degenerate package), 'ehrenfest',
                'fssh', 'iesh' (the latter three land in phase C; the K2 skeleton carries
                elec_adiabatic + the three batch-3 dist members as the library seed).
System level:   systems/<sys>/reg_<sys>_all registers only PES members and termination
                primitives; the method arrives as config%electronic%method = 'fssh' —
                assembly is a pure name-lookup into the method table. Zero method
                registration lines for the system author.
New method    = one new file in methods/ + one thin reg_*.f90 (drag in and it compiles —
                the Makefile wildcard picks it up; prebuild aggregates every reg_*.f90
                under methods/ and systems/ into registry_gen.f90).
System-specific method variants: the algorithm stays generic and single-copy; parameters
                arrive via config / system data at runtime (decoupling wins — E12 spirit).
```

Method-package × ws coupling retains the E10 contract: a member loud-stops on a ws%repr it
does not support (e.g. iesh rejects r_adiab), and on dt semantics it cannot live with (the
E13 open point, carried to phase C unchanged).

## 6. System-definition files

The geometric definition of the reactants is system knowledge — it belongs to the blackbox
material (the natural extension of E11/E12: drag in a folder and the system is fully defined).

| Item | Contract |
|---|---|
| Granularity | one file per reactant; file name = reactant name (A.poscar, B.xyz, ...) |
| Discovery | the assembly phase scans the blackbox folder for known extensions (POSCAR / .poscar / .xyz) — no manifest layer |
| Formats | POSCAR (with cell) and .xyz (no cell) initially; a new format = a new parser member (grammar versioning) |
| Derivation rule | cell present -> config%reactants%surface_model = surface; absent -> molecular. A generic formal rule lawfully in generic code (same class as E13's retained "n_surf >= 2 with r_scalar main PES -> loud-stop") |
| Data flow | scan -> generic parser (utils) -> element symbols against the generic mass table -> fill roster (geometry / cell / masses) -> derive paradigm -> hand equilibrium geometries to sampler (B' placement) |
| Ownership | file contents = blackbox material; parsers / mass table / derivation rule = generic code (zero system knowledge) |

sysdef_facade.f90 owns the contract: discovery, parser dispatch, roster fill, paradigm
derivation, sampler hand-off. This is not a fifth callable entry of the evolution loop — it
is a data-material contract consumed in the assembly phase.

## 7. Termination: two-track composition

- Track 1 (system material): each blackbox registers its termination primitives via
  pred_reg_prim(name, proc) — e.g. "adsorbate height above surface" written by the system
  author. System knowledge lives in the blackbox.
- Track 2 (program-side grammar): boolean composition (AND/OR/NOT, closed set) and the
  MAX_STEPS safety net stay in predicate. Composition logic is program knowledge (pure
  Boolean algebra — the §3 closure argument survives verbatim).
- Assembly: TERMINATE= expressions compose leaves of both tracks on equal footing,
  e.g. AND(PRIM('desorb_z'), MAX_STEPS(500000)).

## 8. Driver orchestration and hooks withdrawal

```fortran
do while (.not. term_flag)
   step_i = step_i + 1
   call prop_step(sta)                  ! phase 1: nuclear propagation (multiple force_eval
                                        ! calls inside; endpoint-force contract obliges the
                                        ! final ws to match the returned sta%q)
   call elec_prop(sta, dt)              ! phase 2: electronic propagation
   call mqc_statistic(sta)              ! phase 3: matching (hop / momentum rescaling /
                                        ! optional force refresh)
   call rec_frame(sta, u_out)           ! frame output (post-mqc self-consistent state —
                                        ! repairs the prior one-step i_surf lag in frames)
   term_flag = pred_eval(sta, expr_term, step_i)   ! phase 4: termination
end do
```

- hooks_t is withdrawn entirely (hook_out_i / hook_term_i / hooks_t / the second prop_step
  argument are deleted): the propagator reverts to a pure integrator; frame output and
  termination timing become explicit. This also repairs a latent timing defect: the old
  hooks%term checkpoint fired after nuclear propagation but before electronic evolution —
  a termination there skipped the step's electronic phases, leaving an inconsistent state.
- E7 phase restriction restated: inside the integration loop the direct calls are prop_step,
  the facade phase entries, rec_frame, and pred_eval; no sampler / analysis / input / roster
  calls (unchanged restriction, now with an explicit call list).

## 9. Affected files

| File | Change |
|---|---|
| src/facade/force_facade.f90 | deleted and rewritten (force_eval entry + pes registry; all class-B content excised) |
| src/facade/elec_facade.f90 | new (elec_prop / mqc_statistic entries, method-package registry, ws holder, sync tail calls) |
| src/facade/sysdef_facade.f90 | new (system-definition contract: discovery/dispatch/roster fill/derivation/hand-off) |
| src/facade/elec_ws.f90 | comments revised (one-step lifetime, holder wording) |
| src/facade/sampler.f90 | kept (comment: initial geometries arrive from sysdef hand-off) |
| src/globals/state.f90 | amp(:) -> a(:,:); rho(:,:) added; i_surf/n_hop relocated from elec_t to state_t top level (classical-side statistical fields); comment statutes (symmetric ontic containers, member-defined semantics, statistical vs ontic distinction) |
| src/dynamics/propagator.f90 | hooks withdrawn; endpoint-force contract in header; dt-channel note rewritten |
| src/dynamics/predicate.f90 | primitive registry (pred_reg_prim + k_prim leaf + PRIM grammar leaf) |
| src/driver.f90 | four-phase orchestration pseudocode; hooks wiring removed; assembly phase gains sysdef step |
| src/globals/config.f90 | method-key comment (name-lookup into the method table); surface_model comment (derived, not input) |
| src/globals/config_roster.f90 | roster_load data source: sysdef facade products (was input keywords) |
| src/input.f90 | comment updates (TERMINATE grammar gains PRIM; surface_model no longer an input key) |
| src/utils/geometry.f90 (or new parser file) | POSCAR / .xyz readers + element mass table (generic) |
| interface_electronic/, interface_distributions/ | rebuilt as flat methods/ (files moved, headers rewritten; elec_adiabatic + dist_* are the seed set) |
| Makefile / prebuild | wildcard picks up methods/; prebuild aggregates reg_*.f90 under methods/ and systems/ |
| READMEs (src/, src/interface_electronic/ -> src/methods/), docs/plans/K2.md, docs/final_completion_plan.md | E14 ruling entries; §5 N1 supersession |

## 10. Deliverable form and verification

The project is in the K2 skeleton stage: the deliverable of this redesign is the
comment-legislated skeleton of the new contract (pseudocode bodies), with physical members
(elec_fssh etc.) still landing in phase C per the existing plan. Verification = the full tree
(24+ units, grown by elec_facade/sysdef_facade and the methods/ rebuild) compiles under both
gfortran and ifx with zero errors, plus a whole-tree comment-consistency audit (the same
standard as the E13 closure).

## 11. Open points (deferred, carried from prior rulings)

- dt-rejection mechanism for fixed-step-only laws (member-declared prerequisite flag vs
  runtime first-step comparison) — phase C with elec_tsh (E13 open point, unchanged).
- U-delivery tension (ws adiab-view cache field vs re-diagonalization inside the class-B
  member at double O(n_surf^3) cost) — phase C; the consumer is now elec_prop /
  mqc_statistic rather than an in-force_eval call (E12 open point, re-anchored).
- PRIM leaf grammar details (quoting, argument passing into system primitives) — the C7
  grammar deliverable window.
- The sysdef parser set grows per format (POSCAR/.xyz first); CIF/gro etc. are grammar-
  versioning matters.
