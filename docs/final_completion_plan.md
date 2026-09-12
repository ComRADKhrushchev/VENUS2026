# VENUS2026 Program Completion — Final Plan

> **Executor's note**: This plan is organized as "phase–task–step", with steps tracked by checkboxes (`- [ ]`).
> Every task ends with "acceptance command passes + git commit"; commit messages reference the audit entry name
> (e.g. `P-B1: GWRITE observables demoted to configuration (audit: GWRITE)`). Before executing, read
> `docs/design_spec_completeness.md` (adjudication standard) and `venus2026/CLAUDE.md` (build/unit conventions).
> **Presentation protocol** (user ruling 2026-09-06): every presentation of a design/skeleton/implementation follows
> "outline first, confirm, then expand" — present only the current layer's outline at a time, and expand the next
> layer only after user confirmation; dumping all concrete content at once is forbidden.
> **Execution strategy** (user ruling, same day): skeleton-first — Phase K first generates all placeholder files
> (real definitions of global variables + input/output signatures and pseudocode for every routine); after the
> structure review confirms, concretize file by file.
> **Standalone new program** (second revision of the user ruling, same day): the goal is to deliver a standalone
> new program; the legacy code (venus2026/) serves purely as reference (provenance of algorithms and formulas)
> and is read-only throughout — no modification, no cleanup, no migration, no retirement.
> **Per-phase plan files** (user ruling, same day): before executing each small phase (e.g. K1, C2), first submit
> a separate plan file for that phase, `docs/plans/<phase>.md` (goals/deliverables/file list/steps/acceptance),
> and execute only after confirmation; this plan retains only the sequence, the mapping, and the yardsticks.

Status: **execution baseline** (this document is the sole basis of execution)
Date: 2026-09-06 (same day, remediation_plan.md and audit_brief.md were deleted — the former's purification steps
have been fully merged into this document, the latter was a completed review brief; this avoids dual-track
execution and document confusion)
Dependency chain: design_spec_completeness.md (spec) → completeness_audit.md/.json (material inventory, 99 entries)
→ optimal_structure.md (target structure and material mapping, **the baseline of this plan**)

**Goal:** deliver a **standalone new program** — the target structure is the L0-L5 layering of optimal_structure.md
§1, with variable code concentrated in interface A/B′ members and closed-family registration lines. The legacy code
(venus2026/) serves purely as reference. **Verification system rebuild** (user ruling 2026-09-06): all 10 old cases
are discarded as baselines; new verification cases and acceptance criteria will be drafted one by one alongside the
C components once the architecture (K skeleton) is in place, then assembled into a set and fully executed in Phase V.

**Architecture:** standalone new program + skeleton-first (user ruling 2026-09-06, twice revised). The new program
is its own source tree (with its own build); the legacy codebase venus2026/ is a read-only reference. Phase K
generates the full L0-L5 placeholder skeleton (real definitions of global variables + signature/pseudocode per
routine; the new tree compiles independently), then — after the review gate — Phase C concretizes file by file
(following the §3 mapping back to the legacy reference algorithms, with component-level old-vs-new cross-checks),
and finally Phase V executes the rebuilt verification case system (criterion coverage check + full runs). Before
each small phase executes, a separate plan file is produced first.

**Tech stack:** Fortran (ifx, verified with oneAPI 2025.2), make (single-entry Makefile), git.
Testing-method note: this project's "test loop" = new verification cases (minimal input + analytical/statistical
anchors + conservation-law criteria, drafted alongside the C components) + ad-hoc same-quantity cross-checks
(comparison against the old program when it can compute the same physical quantity; not a baseline), not carried
by a unit-test framework.

**Spec:** docs/design_spec_completeness.md; structure mapping details in docs/optimal_structure.md §2.

## Global constraints (implicitly obeyed by every task)

1. Internal unit system Å/amu/(10 fs), conversion constants C1=0.04184, C7=0.063508, 23.0605 (eV→kcal/mol),
   C8=1.987e-3; conversions happen only at the force_facade boundary (the point where the unit contract is honored).
2. Build with `make` (under venus2026/; on Windows run in Git Bash, linker path variable overrides see CLAUDE.md);
   flags `-r8 -double-size=64 -i8 -O2 -w`.
3. Input accepts only the `KEYWORD=VALUE` grammar (input_parser); new keywords always go through the
   param_mapping registry pattern; fallback to the old integer keywords is retained.
4. No operator-class calls during the evolution phase (phase restriction, spec 2.5); preparation members internally
   calling evolution services is legal (boundary 4).
5. Single random stream: the new program is single-stream from C1 onward (rng.f90, Knuth RAND0 family, ISEED
   state); no second stream anywhere in the codebase.
6. After editing a module file, delete stale `.mod` before rebuilding (CLAUDE.md convention).
7. Line-number references follow the 2026-09-04 audit; if line numbers drift during execution, re-locate by
   grepping symbol names.
8. Each C step closes with "minimal verification input + analytical/statistical criteria" (criteria first: define
   the criteria before writing the implementation); V/S run the new case set over the full chain; the old program
   serves only as an ad-hoc same-quantity reference, not as a gate.

---

## 0. Baseline verification (2026-09-06, differences from the audit)

Verified before execution; the following are established facts, and the plan adjusts accordingly:

| Verified item | Current state | Implication for the plan |
|---|---|---|
| `src_VENUS/GLO.f` (uncompiled leftover copy) | **no longer in the tree** | the corresponding P-A1 entry is naturally complete |
| `ROTATEX_b/ROTATEZ_b` | **no longer in the tree** | same |
| GASDEV.f/VFDATE.f/POTEN.f/EBOND.f/THERMBATH.f/WEBOND/WLBOND/WENMOD/MICROCI/VOLPSCONE/LMODE/LMEXCT | still present | P-A1/A2 delete as planned |
| gasdev lives at interface_TEST.f90:112 (function) | still present | P-C1 proceed as planned |
| GWRITE :59/95-99/110-114, TEST :88-94, VERLET :47,74,137 ghost/GLO embedding | still present, line numbers consistent with the audit | P-B/P-D proceed as planned |
| Makefile MKL three-library linking | still present | P-A5 proceed as planned |
| `venus2026/.git` | **does not exist** | Phase 0 must establish version control first |
| cases/ count | **10** (boltzmann_vib, ebk_fixed_nj, fixed_normalmode_qnums, integrator_matrix, leps_potential, mb_thermal, rigid_surface, rst_beam_scattering, rst_surface_md, twobody_collision) | **all discarded as verification baselines** (user ruling 2026-09-06: the old cases are "old-behavior preservation" use cases, with uneven grounds for judgment and imprinted with defects such as F15/dual random streams); the directory stays in the legacy tree as historical reference; the new case system will be drafted after the K skeleton. CLAUDE.md's "16" is outdated; P0 corrects it |
| docs/completeness_audit.json | present | subject of the S2 re-audit |

## 1. Target end state (definition of the acceptance state)

Structure (excerpt from optimal_structure §1; the file names are the placement targets):

```
L0 control kernel        driver.f90          the assembly→evolution→statistics call sequence, no physics
L1 thing-in-itself container state.f90       Q, P, s, F, t (allocatable; dimensions = roster content parameters)
L1 representation layer  config.f90          scenario, roster, distributions, predicates, outcome classification, observables configuration
L2 interface A facade    force_facade.f90    (Q,P,s)→(F,E); internally Σ force-contribution members × single-choice electron propagation law
L2 interface B′ framework sampler.f90        distribution→(Q0,P0,s0) generic sampler
L3 registry-family members  pes_*.f90 / dist_*.f90
L3 closed families       propagator.f90 / predicate.f90
L4 measurement           recorder.f90        trajectory serialization + predicate evaluation
L4 post-processing       analysis/           outcome classification execution, statistics, plotting (offline-capable)
L5 instruments           rng.f90, linalg.f90, geometry.f90
```

End-state criteria (T1-T7, verified item by item at closing S):

- **T1 file structure**: all files in the table above present; Makefile SRC ordered by the new dependency order;
- **T2 variable-code location**: adding a new PES/distribution member = a new file + one registration line, with
  zero touching of existing code (verified by the "add a pes_leps variant test piece" drill);
- **T3 dual-implementation matrix**: A axis HARMONIC+MORSE+LEPS+EMT-NN+RST (≥2 ✓); B′ axis
  THRMAN+INITEBK+dist_surface_oscillator (≥3); at the configuration level, molecule vs surface, two paradigms
  (each with a representative among the new verification cases);
- **T4 unit contract**: grep the whole codebase — 23.0605 and the C1 force conversion appear only at the
  force_facade boundary;
- **T5 single random stream**: no `RANDOM_NUMBER` residue anywhere (rng.f90 is the only stream in the codebase);
- **T6 three phases**: the driver contains only the assembly→evolution→statistics call sequence; inside the
  evolution loop only facade/propagator/recorder/predicate;
- **T7 completeness re-audit**: re-audit every routine of the new program against the spec's three questions,
  producing a new completeness audit (md+json); nothing remains outside the invariant and member classes. The old
  audit is kept as historical reference.

## 2. Execution sequence (user ruling 2026-09-06: standalone new program + skeleton-first + per-phase plan files)

```
Phase 0 preparation (git + new-program tree setup)
   0.1 git init (workspace root; the legacy code is read-only from here on)
   0.2 new-program tree setup (src/: Makefile skeleton + README + analysis/ placeholder)
   0.3 legacy CLAUDE.md case-count correction
   ↓
Phase K skeleton generation (all L0-L5 files as placeholders in the new tree)
   K1 real definitions of global variables (state.f90 / config.f90 / consts;
       block-by-block placement of venus_data's ~40 representation blocks)
   K2 placeholder per routine: input/output signature + pseudocode + provenance comment (legacy file:line)
   K3 the new tree compiles independently (stub main program; the compiler verifies interface consistency)
   ↓ [review gate: show only the skeleton outline (file tree + one-line routine list per file);
        release after layer-by-layer confirmation; pseudocode expanded per file under the presentation protocol]
Phase C concretization (file by file; per file, outline first then expand after confirmation; component-level old-vs-new cross-checks)
   C1 L5 instruments: rng (single random stream) / linalg / geometry   ← reference RAND0 family/gasdev/EIGN/ROTATE family
   C2 force_facade + elec_adiabatic + pes_* members       ← reference interface_TEST/emt_nn/rst/VERLET embedded forces
   C3 predicate grammar + termination-criteria assembly   ← reference TEST.f:106-146
   C4 recorder + config_roster + input grammar            ← reference GWRITE/HIST (stripping the violations listed in the audit)
   C5 sampler + dist_* member table                       ← reference SELECT/THRMAN/INITEBK/SURF/…
   C6 propagator closed family (verlet/symple/radau)      ← reference VERLET/SYMPLE/RADAU (orchestration stripped)
   C7 driver + outcome-classification data                ← reference VENUS.f90 three-phase orchestration
   C8 analysis/ offline family                            ← reference FINAL/GFINAL/EIGOUT/plot_*
   ↓
Phase V full-chain verification (the C-phase per-component verification inputs and criteria assembled into the
   formal case set; criterion coverage check + full runs + same-ISEED reproducibility spot checks)
   ↓
S closing (new-program completeness re-audit / document wrap-up / spec §4 backfill)
```

Key points:
- **Legacy code read-only**: venus2026/ is not modified, its dead code is not cleaned, nothing is retired (the
  original C0/T old-tree operations of both kinds are voided); the audit's file:line references serve only to
  locate reference code.
- **Per-phase plan files**: before executing each small phase, first submit `docs/plans/<phase>.md`
  (P0/K1/K2/K3/C1…C8/V/S) and execute after confirmation; this plan retains the sequence, the mapping, and the
  yardsticks.
- **Verification system rebuild** (user ruling 2026-09-06): all 10 old cases are discarded as baselines; each C
  step drafts its own minimal verification input + criteria (analytical/statistical anchors > conservation-law/
  distribution checks; the criterion design may borrow the anchor ideas from the old case READMEs without reusing
  their inputs and output baselines), and Phase V assembles them into the formal case set. Ad-hoc same-quantity
  cross-check: when the old program can compute the same physical quantity it may be used for comparison (not a
  baseline, not a gate).
- **Input grammar**: extend the old KEYWORD=VALUE as the base with constructor new keywords
  (TERMINATE/OBSERVABLES/ELEMENTS/…, grammar finalized at C4); when temporarily cross-checking the old program,
  old inputs may be reused.
- **Skeleton drift ban**: if a signature turns out to be unsuitable during concretization, go back to the review
  gate, present the difference, and revise after confirmation.
- The C order is the dependency order: instruments→facade→predicate/recorder→sampling→propagation→driver→
  post-processing; C2 needs C1 (members internally call rng), C7 needs C2-C6 (the driver only orchestrates).
- §5 (formerly N1-N3) = the interface specs of the three key skeleton files; §6 (formerly P-A…P-G) = the
  stripping/operation points when porting (the audit violations are exactly the "do not carry over when
  referencing" list); §7 (formerly the X track) = the placement check for the new tree.

## 3. Placement mapping table (legacy concern → new-tree placement)

(Second revision, 2026-09-06: in the standalone-new-program context, this table's semantics = where each concern
lands in the new tree; the "original placement" column keeps the historical provenance. The legacy tree
venus2026/ is read-only; all target files below refer to the new tree.)

| Step | Original placement | Final placement |
|---|---|---|
| B1 GWRITE observables as configuration | GWRITE internal table | recorder.f90 column registration (the mechanism produced by N3) |
| B2 TEST generic criteria registration | TEST.f internal | predicate.f90 grammar member (produced by N2) |
| C1 gasdev/GAMA migration | "src_VENUS standalone instrument module" | rng.f90 (L5, co-located with RAND0/RAND1/RANDST) |
| D1/D2/D4 stochastic-law members | interface A members (location undetermined) | pes_glo_spring / pes_langevin_ghost / pes_ldfa, registered into force_facade.f90 |
| D3 ghost parameter assembly | "member init hook" | facade member init protocol (POTPRE-style) |
| E1 outcome block merge + classification | VENUS.f90 single path | driver calls the analysis outcome-classification executor + △6 classification data |
| E2 observable accumulation | "measurement recorder module" | recorder.f90 final-state columns |
| E3 system branching | "assembly layer" | configuration assembler (roster, the N3 foundation) |
| G1/G3 representation write-back turned into delivery | "input file or configuration assembler" | configuration assembler + config grammar |
| G4 RADOUT callback | "callback to the control kernel" | the driver's output hook + predicate evaluation point |

Conflict adjudication (where the two documents disagree, the baseline optimal_structure prevails; after the second
revision of 2026-09-06, adjudications 1/2 dissolve automatically — the old tree is untouched and the new tree maps
directly per the baseline, so there is no delete/keep dispute):

1. **ENMODE**: remediation A2 deletes it (MPLOT dead); optimal_structure judges ◐ to be placed with
   dist_normalmode. Adjudication: **do not delete it along with A2**; its MPLOT dead call chain is deleted, and at
   X3 it is placed as dist_normalmode's assembly-time analysis tool; if at that point the whole codebase truly has
   no consumer, delete it per the audit's original ruling and record that.
2. **FRICFORCE/DIAG2**: A2's deletion (dead entry) and D4's revival (LDFA member) coexist. Adjudication: in A2,
   **delete** from VERLET.f90 (together with DIAG2/LINEQ/LUDCMP/LUBKSB/FEHLER); the physics is reborn in D4 as
   the pes_ldfa member in the facade (the material is retrievable via git history).
3. **case count**: 16 → 10 (see §0).
4. **harmonic_sampling**: the audit rules "wire into B′ or delete". Adjudication: wire it in — X2 registers it as
   dist_surface_oscillator (the B′ third implementation) and adds a row to the SELECT member table.

Numbering correspondence: what optimal_structure §5 calls R0→R7 corresponds to the original (deleted)
remediation_plan's P-A…P-G + closing; this document keeps the P-X numbering and adds N-X/X-X/S-X.

---

## 4. Phase 0 preparation (version control and new-program tree setup)

**Files:** workspace root `.git` (newly created), `docs/plans/` (per-phase plan file directory),
new program tree `src/` (naming finalized in the P0 plan file).

- [x] **0.1 git initialization**: workspace root `git init && git add -A && git commit` (docs/, the venus2026/
  reference, and README under version control; the legacy code is read-only from here on). Optionally connect a
  remote (ask for user consent before pushing).
- [x] **0.2 new-program tree setup**: `src/Makefile` (ifx, flags as per the legacy conventions, the Windows
  linker-masking section ported from the legacy Makefile's platform section) + `src/README.md` (L0-L5 overview,
  placeholder note: generated in Phase K) + `src/analysis/README.md` (L4 purpose).
- [x] **0.3 legacy doc number correction**: venus2026/CLAUDE.md "16 example cases"→10 (avoiding a misleading
  reference; this is the only allowed change to the legacy tree).
- [x] Commit `P0: version control and new-program tree setup`. Execution details in docs/plans/P0.md (v2).

## 5. Skeleton interface specs: facade / predicate / recorder+roster (formerly track B N1-N3)

(Revised 2026-09-06: N1-N3 are merged into Phase K as the interface specs of the three key skeleton files; the
check_* verification programs move to the acceptance of the corresponding C steps. The signatures and grammar
definitions in this section take direct effect for K.)

### N1 force_facade.f90 + elec_adiabatic (△1+△5)

**Files:**
- Create: `src/force_facade.f90` (new tree L2)
- Create: `src/elec_adiabatic.f90`
- Test: `scripts/check_facade.f90` (component-level cross-check, independently compiled)

**Interfaces (produced; later tasks depend on them):**
```fortran
module force_facade
  public :: facade_register_pes(name, pes_proc)      ! class-A member registration (additive)
  public :: facade_register_elec(name, elec_proc)    ! class-B member registration (single choice)
  public :: facade_init(config)                      ! assemble the member table + formal validation (E13: names resolve, slots filled, container-contract contradictions caught; physical-combination legality is the members' own affair)
  public :: facade_force(Q, P, s, F, E)              ! (Q,P,s)→(F,E); F=Σ member contributions; E=T(P)+ΣV
  public :: facade_info                             ! member list / capability levels (diagnostic output)
end module
! Class-A member signature: Q,W → F_add, E_add; stochastic-law members call rng internally (the signature acknowledges stochastic laws)
! [Unit-contract revision (2026-09-06 K2 batch-2 ruling, confirmed by the overall review)] Members deliver an
! eV force core (g_add[eV/Å], v_add[eV]); the ×23.0605×C1 conversion is centralized in the facade (ev_kcal being
! the sole occurrence point in the codebase, the accumulation form identical to DPESHON) — N1's original "code
! units" annotation is to be understood as "after facade conversion"; T4 follows from this. C2 implementers take
! this as authoritative.
! Class-B member signature (degenerate point): elec_adiabatic: ds/dt≡0, identity force rule — explicitly present, ~20 lines
```
[E13 revision 2026-09-08 — supersedes the capability-level negotiation table] The facade is a purely
formal authority: name->pointer dispatch plus checks derivable from the container contract alone. There is
no legal-combination table inside facade_init — physical-combination legality is the members' own knowledge
(the class-B member loud-stops on a ws%repr it does not support and on dt semantics it cannot live with;
fixed-step-only laws such as TSH/IESH exclude radau this way). Blackbox authors are assumed to know their
system, not this program; a facade-side table would also break E11 drag-in (a new combination would require
editing force_facade). config%electronic%method is a name string ('adiabatic'/'tsh'/'iesh'), same form as
pes_main.

[E14 revision 2026-09-08 — supersedes the single-entry orchestration; docs/superpowers/specs/2026-09-08-facade-e14-design.md
§3/§8] The blackbox contract becomes four entries — `force_eval` (state→forces), `elec_prop` (electronic
evolution, owns dt), `mqc_statistic` (electronic-state→state statistics), termination composed through
predicate's two tracks — and the orchestration is explicit in the driver (the four-phase integration-loop
step prop_step -> elec_prop -> mqc_statistic -> rec_frame -> pred_eval; hooks withdrawn). The
signature/registration wording above is historical: facade_register_pes / facade_register_elec and
facade_force(Q, P, s, F, E) are superseded by pes_reg / method_reg (method packages registered by thin
reg_*.f90 files; the K3 prebuild aggregates every reg_*.f90 into registry_gen.f90) and by the state-based
entries force_eval(sta, e_tot) / elec_prop(sta, dt) / mqc_statistic(sta), with blackbox_init /
blackbox_info as the assembly entries. All negotiation/registration contract wording in this section now
points to the spec as the ruling authority.

- [ ] **N1.1** Write force_facade.f90: the member registry (name→procedure array), the additive composition
  loop, formal validation (E13 — no capability-negotiation table; physical-combination legality is checked
  by the members themselves). The conversion expression exactly matches DPESHON in form
  (`-grad×23.0605×C1` accumulated), guaranteeing bit-level comparability.
- [ ] **N1.2** Write elec_adiabatic.f90 (degenerate member: ds/dt≡0, identity force rule).
- [ ] **N1.3** Write the five class-A member files: pes_reference (the system-independent reference potentials
  HARMONIC/MORSE/LEPS, referencing test_potentials.f90), pes_emt_nn, pes_rst (algorithms ported from the
  corresponding legacy files; externalized data and path keywords carried over), registered into the facade.
- [ ] **N1.4** Write scripts/check_facade.f90: build a temporary input, run the old program for a reference
  trajectory (≥100 frames; or hand-compute analytical gradients for the reference potentials), compute F/E at
  the same Q via both `DVDQ_1→DPESHON` (old path) and `facade_force` (new path), and diff element by element.
- [ ] **N1.5 Acceptance**: `make` passes; check_facade gives max|ΔF|=0 and max|ΔE|=0 over all frames (bit level;
  if nonzero, check the conversion ordering and fix until zero). Commit `N1: force_facade + degenerate electron
  member (bypass verification, bit-level identical)`.

### N2 predicate.f90 grammar interpreter (△2)

[E18 withdrawn 2026-09-12, K2 review gate: the grammar interpreter is deleted — termination is
the container's own judgment (container_bind_term / container_term, logical term_i(sta)); MAX_STEPS
remains as control%max_steps in the driver loop condition. Section kept for the historical
interface intent; phase C3 is withdrawn likewise — TEST.f:106-146 becomes the provenance of the
C2 container example term function. The N2 block below is the pre-E14 historical text.]

**Files:**
- Create: `src/predicate.f90` (new tree L3 closed family)
- Test: `scripts/check_predicate.f90`

**Interfaces (produced):**
```fortran
module predicate
  public :: predicate_parse(str, expr, ios, errmsg)  ! input-grammar expression → expr tree
  public :: predicate_eval(Q, P, W, roster, expr) → logical
  ! [E14 annotation 2026-09-08] These signatures were never implemented as written: the landed module
  ! (src/dynamics/predicate.f90, K2) exports pred_parse / pred_eval — pred_parse matches the first line's
  ! argument list modulo the name; pred_eval is state-based, pred_eval(sta, expr [, i_step]) (one state_t
  ! argument, no Q/P/W/roster dummies; i_step = the MAX_STEPS safety-net channel) — plus pred_reg_prim
  ! (E14 track 1). Block kept for the historical interface intent.
end module
```
Input grammar (new keyword `TERMINATE =`, registered in param_mapping):
```
expression := atom | AND(e1,e2) | OR(e1,e2) | NOT(e)
atom   := COM_DISTANCE(fragment A, fragment B, threshold Å[, recession-velocity threshold])   ← TEST.f:106-146 prototype
        | Z_HEIGHT(fragment, threshold Å[, OUTWARD])                  ← replaces the TEST.f:88-94 hardcoded criteria
        | MAX_STEPS(count)                                   ← replaces the NC magic number (always a safety net)
        | BOUND(fragment A, fragment B, lower, upper)        ← bounded window (constrained scattering)
        | PRIM('name')                                       ← [E14] registered system primitive (two-track composition,
                                                                  must match a pred_reg_prim row; quoting/argument
                                                                  conventions finalized at C7)
fragment   := roster fragment name (A/B/…) or ATOM(i)
```
[E14 revision 2026-09-08 — supersedes the closed-primitive-set wording; spec §7] The primitive set is an
open registered set: system termination primitives enter as `PRIM('name')` leaves (pred_reg_prim) composing
with the four built-ins on equal footing (e.g. AND(PRIM('desorb_z'), MAX_STEPS(500000))); the boolean
combinators (AND/OR/NOT) remain a closed set and MAX_STEPS stays the ever-present safety net. Nonterminal
divergence note: the landed grammar (src/dynamics/predicate.f90) names the leaf nonterminal `primitive`
(`expression := primitive | AND(e1,e2) | OR(e1,e2) | NOT(e)`) — this section's `atom :=` is the historical
spelling of the same language.

- [ ] **N2.1** Write the grammar data structures (atom/boolean-combination nodes) + parser + evaluator (the COM
  distance and recession-velocity algorithms taken line by line from TEST.f:106-146, guaranteeing numerical
  provenance).
- [ ] **N2.2** The default expression when no `TERMINATE` keyword is present = a same-parameter version of the
  current TEST generic criteria (existing parameters such as RCRIT mapped into the expression), guaranteeing
  unchanged behavior before B2.
- [ ] **N2.3** Write check_predicate.f90: replay over a temporarily generated old-program trajectory (or a
  constructed state sequence), asserting that the boolean at every step matches the current TEST.f output.
- [ ] **N2.4 Acceptance**: `make` passes; check_predicate consistent over the whole trajectory. Commit
  `N2: predicate grammar interpreter`.

### N3 recorder.f90 column registration + roster (△4+△3 foundation)

**Files:**
- Create: `src/recorder.f90`, `src/config_roster.f90` (new tree)
- Test: the component-level cross-check happens when the wiring is done at C4; this item's self-test is compile + smoke

**Interfaces (produced):**
```fortran
module recorder
  public :: recorder_register_col(name, col_proc)    ! column = name + extractor
  public :: recorder_select(cols_spec)                ! OBSERVABLES= list or GWRITE_LEVEL preset group
  public :: recorder_frame(Q, P, F, unit)             ! write one frame of the selected columns
end module
module config_roster
  public :: atom_spec(symbol, mass, fragment)         ! roster element structure
  public :: roster_load(ELEMENTS= list) / roster_symbol(i) → 'C '/'Au ' (original format width)
end module
```
New keywords: `ELEMENTS = C, AU(111), …` (order corresponds to the atomic index; the fallback default is `AT<i>`
labels); `OBSERVABLES = column names,…` (coexisting with the GWRITE_LEVEL preset groups, the latter referencing
the column-group table).

- [ ] **N3.1** Write the recorder column-registration mechanism (the column-group preset table as data: level
  0/1/2/3 → column lists).
- [ ] **N3.2** Write config_roster: the element symbol table (symbol→2-char label; masses still go through the
  current W assembly), and the ELEMENTS keyword hookup.
- [ ] **N3.3 Acceptance**: `make` passes; all case regressions unchanged (nothing consumes the new module yet).
  Commit `N3: recorder column registration and roster data structures`.

## 6. Material operation list: the original purification main chain (P-A…P-G)

(Second revision, 2026-09-06: the legacy tree is read-only, and the steps in this section are not executed in the
legacy tree. Their role = **the stripping/operation points when porting for reference** — the violations listed in
the audit are exactly the "do not carry into the new tree when referencing" list, and the line numbers and entry
names serve to locate the reference code. The whole P-A dead-code-removal section is void: the legacy tree is not
cleaned and the new tree has no dead code to clean (the audit A8 list = the list of functionality the new tree
need not implement). The "acceptance" lines embedded in each task are the old plan's old-tree acceptances,
superseded by the new acceptance approach of §2/§9 (component criterion verification + the V full-chain new case
system); only the step essentials remain valid.)

### P-A dead-code removal (zero physics risk)

**Files:** delete the files listed in A1/A2; edit VENUS.f90, venus_data.f90, VERLET.f90, GLO.f90, Makefile.

- [ ] **A1** Delete files: `src_VENUS/{GASDEV.f, VFDATE.f, POTEN.f, EBOND.f, THERMBATH.f, WEBOND.f, WLBOND.f,
  WENMOD.f}` (GLO.f and ROTATEX_b/ROTATEZ_b are already out of the tree; verify the Makefile SRC has no residual
  references).
- [ ] **A2** Delete dead entry points and their dependency chains: MICROCI.f, VOLPSCONE.f, LMODE.f, LMEXCT.f
  (EBOND already in A1); delete ENMODE.f's MPLOT call chain (the file itself kept pending the X3 adjudication,
  see §3 conflict 1), and clean venus_input.f90:692's MPLOT=0 hardcoded keyword as well; PBC in SURF.f;
  FRICFORCE+DIAG2+LINEQ+LUDCMP+LUBKSB+FEHLER in VERLET.f90 (the physics material is reborn in D4); GLOEQU and
  DISTRIBUTE in GLO.f90.
- [ ] **A3** VENUS.f90 dead blocks: the NCHKP branch (:182-232), the RKUTTA coefficients (:86-94), the idplt
  histogram dead block (:286-299, including the Qr1/Qr2 site constants and the SYSTEM shell call), the NCVIB dead
  keyword.
- [ ] **A4** venus_data dead blocks: MORSEB, STRETB, CUBEB, COORS, GAUSS2/GAUHES, GDD, HIST_ATOM_NAD, U0/U1/V01,
  the RKUTTA block.
- [ ] **A5** Makefile: remove the MKL three-library linking (no calls anywhere) and the SRC entries of deleted
  files.
- [ ] **Acceptance**: `make` passes; `bash scripts/regress.sh` — all 10 cases PASS (bit level). Commit
  `P-A: dead-code removal (audit A8 list)`.

### P-B measurement-layer purification (observable selection demoted to configuration)

**Files:** edit GWRITE.f90, TEST.f, GFINAL.f, FINAL.f, POTENZ.f, HIST.f90; wire in recorder/predicate/roster.

- [ ] **B1** GWRITE stripped of system knowledge: delete the atom1=C special-case lines (:59-60), the
  mass-threshold element names (:110-114, use roster_symbol), the nearest-Au diagnostics (:95-99, moved to
  analysis as an optional offline tool), the +6.0 Å PBC box (:78-79, made configurable), and the :43 inline kB
  (use the venus_params constant); the per-level observable sets of GWRITE_LEVEL go through the recorder
  column-group table.
- [ ] **B2** TEST.f: delete the live C/Au criteria (:88-94) and the two commented system-knowledge blocks
  (:54-78, :149-192); the termination criteria are assembled via predicate.f90 (without a TERMINATE keyword, the
  N2.2 default expression applies). The affected cases (rst_beam_scattering, rst_surface_md, rigid_surface and
  the other NSURF>0 ones) get a TERMINATE line added to input_qct.txt that is equivalent to the current criteria
  (thresholds at current values; the equivalent translation follows the current TEST.f text).
- [ ] **B3** GFINAL.f:80 — delete the hardcoded atom-3 EREL override and adopt FINAL's generalized COM relative
  energy; also delete the :66-69 uninitialized current_U0 print.
- [ ] **B4** Delete the fort.777/666 U0/U1/V01 prints (FINAL.f:131-135, POTENZ.f:76-82 — no writers anywhere).
- [ ] **B5** record_all/HIST: auro=2.94998 and the hexagonal lattice (PI/3, top/brg/hcp/fcc) become slab roster
  parameters (new keyword `SURFACE_LATTICE = HEX(auro)|FCC(a)`, into config).
- [ ] **B6** (addendum, audit FINAL row) delete the FINAL.f:616-620 NZDOWN hardcoded atom-3 dead branch (the
  input-layer loud-stop stub is kept).
- [ ] **Acceptance**: regress all PASS — rst-type case trajectories numerically identical except for the
  element-label source change; GWRITE element names consistent with the ELEMENTS roster (after adding the
  ELEMENTS line to the rst case inputs, the diff is clean). Commit `P-B: measurement-layer purification (audit:
  GWRITE/TEST/GFINAL/FINAL/record_all)`.

### P-C random-stream unification (reproducibility)

**Files:** create rng.f90; edit the consumers interface_TEST.f90, GAMA.f, VERLET.f90, VENUS.f90, GLO.f90,
THERMO.f, etc.

- [ ] **C1** rng.f90 (L5): RAND0/RAND1/RANDST migrated as-is; gasdev moved in from interface_TEST.f90:112 and
  switched to the main stream (internal uniform source RANDOM_NUMBER→RAND1; Box-Muller and the cache state
  kept); GAMA moved in as well (the IDUM dummy argument made effective again).
- [ ] **C2** Single-stream verification over the whole codebase: grep `RANDOM_NUMBER` — only one occurrence
  allowed, at rng.f90's internal conversion point (zero if everything goes through RAND1); all consumers (VERLET
  noise, GLOSELECT, THERMO Andersen, QMMICRO, etc.) go through rng.f90.
- [ ] **C3** Fix the VENUS.f90:184,231 RANDST double call into a single call after read-in (F15 closed).
- [ ] **Acceptance**: the noise/thermostat cases (rigid_surface, mb_thermal, rst_surface_md — re-verify the list
  against each input's NFC/NGLO/NTHERM at execution time) are bit-level identical across two runs with the same
  ISEED; the other cases are bit-level identical with the baseline. **Expected baseline replacement**: the
  baseline snapshots of those three cases are updated, and the README is annotated "P-C stream switch, old
  stream not reproducible". Commit `P-C: single random stream (audit A2; F15 closed)`.

### P-D stochastic laws into interface A (members registered into force_facade)

**Preconditions**: N1 done; P-B/P-C done (constants parameterized first, stream switch first).

**Files:** create pes_glo_spring.f90, pes_langevin_ghost.f90, pes_ldfa.f90; edit VERLET.f90, ENERGY.f,
venus_input.f90, force_facade.f90 (registration lines).

- [ ] **D1** VERLET.f90's ghost Langevin (:44-50,71-79,153-161) and GLO spring+random+friction (:137-148)
  extracted into members: pes_langevin_ghost (stochastic law), pes_glo_spring (spring+random+friction); the atom
  indices parameterized via the member ctx (no longer hardcoded atom2/NATOMA+1/+2); the members call rng
  internally. The registration lines go into facade_init (mapped from the NGLO/NFC switches to member assembly,
  grammar backward compatible).
- [ ] **D2** Delete the ENERGY.f:40-47 GLO spring potential; the energy side is summed via the facade
  (E=T(P)+ΣV, the kcal/mol reporting convention unchanged: same scale after T/C1).
- [ ] **D3** venus_input.f90:616-621 ghost parameter assembly (W(2) hardcoded, the 413.4 conversion) moved into
  the member init hook (called at facade_init).
- [ ] **D4** pes_ldfa.f90: the FRICFORCE physics (NFC=1/2, two laws) reborn as the third stochastic-law member.
- [ ] **D5** Clean VERLET's residual CALTYP branches (:34-36,60-62,123-125).
- [ ] **Acceptance**: the rigid_surface energy distribution is bit-level identical to the post-P-C baseline (same
  stream, same arithmetic order); `grep -n "NGLO\|ghost" VERLET.f90` gives zero hits — the integrator proper
  carries no system knowledge. Commit `P-D: stochastic-law members into interface A (audit:
  VERLET/ENERGY/FRICFORCE)`.

### P-E control-kernel cleanup (VENUS.f90)

**Files:** edit VENUS.f90; wire in recorder/config_roster/predicate.

- [ ] **E1** Merge the duplicated outcome-handling + statistics blocks (:506-576 and :581-662) into a single
  outcome path.
- [ ] **E2** Observable accumulation (BOUNCE/z_pos_min hardcoded atom 3 :399-400,454-458; the AE* arrays
  :544-548,631-635,671-677) moved into recorder final-state columns; the atom indices demoted to configuration.
- [ ] **E3** System branching (GLO layout :322-335, surface alignment :340, the THTA special case :344-346)
  handed over to the configuration assembler (roster/assembly time).
- [ ] **E4** Remove the empty check statements (:462-464).
- [ ] **E5** (=△6) the outcome-classification if-chain made data: `OUTCOME = name: predicate expression` rows ×
  count buckets (reusing the N2 grammar); VENUS.f90's outcome path only does a table lookup and count.
- [ ] **Acceptance**: regress all PASS; inside VENUS.f90 grep finds no lattice constants/element names/hardcoded
  atom numbers; only the assembly→evolution→statistics orchestration + a few loop variables remain. Commit
  `P-E: control-kernel purification (audit: VENUS)`.

### P-F global-state partitioning (venus_data)

- [ ] **F1** venus_data explicitly partitioned in three via comments: the thing-in-itself zone
  (Q/P/PDOT/QDOT/Q_old — annotated →state) / the representation zone (→config) / the instrument zone (→rng), with
  each block's line comment stating its X-track destination.
- [ ] **F2** The NN_POT_TYPES indices (venus_params.f90:88-98) migrate into the emt_nn member file; auro already
  moved out with B5 — verify.
- [ ] **F3** Dead-block cleanup cross-checked jointly with A4.
- [ ] **Acceptance**: `make` passes + regress all PASS. Commit `P-F: global-state partition annotations`.

### P-G input grammar and configuration-assembly wrap-up

- [ ] **G1** test_potentials.f90's members writing back representation (BOXLX/BOXLY/SKEW/NATOMB/QZB,
  :152-163,260-293) changed to delivery via the input/configuration assembler; the PES members keep only Q→V/F
  and data reading.
- [ ] **G2** JMAXCALC.f:35-48's H2/D2 nuclear-spin statistics (mass discrimination) demoted to input data
  (`NUCLEAR_SPIN_WEIGHTS =` table or data file).
- [ ] **G3** The BOXLX/BOXLY default geometry (venus_input.f90:997-998) and the NATOMB hardcoded domain
  (:703-704) made mandatory (a missing value loud-stops naming the keyword) or roster-derived; SURF.f:48-50's
  SKEW=π/4 default magic number made mandatory in sync; venus_input.f90:752's W=1D30 rigid mass delivered
  explicitly by the assembler per RIGID-model semantics.
- [ ] **G4** RADAU.f's RADOUT embedded orchestration (:438-498) changed to a callback to the control kernel
  (output hook + predicate evaluation point, received on the driver side); 'CALLING RADOUT' (:112,380,428) and
  the leftover debug prints in SELECT/SURF/ROTen/THERMO/RADAU/FINLNJ cleaned; PRINFO's ~70 lines of commented
  residue cleaned; ROTEN :101-105 debug prints cleaned; the ROTATEJKM trailing laser comment block cleaned.
- [ ] **Acceptance**: regress all PASS; the EMT-NN/RST cases numerically unchanged (member-internal legality
  preserved). Commit `P-G: input grammar and configuration assembly wrap-up (audit:
  test_potentials/JMAXCALC/RADAU)`.

## 7. Placement checklist: the former convergence track X

(Second revision, 2026-09-06: X1 merged into K1 (native state/config definitions); X2/X3/X4/X5/X6 merged into
C7/C5/C6/C1/C8; X7 merged into P0/C7 (new-tree build and EXTENDING). This section is kept as an item-by-item
new-tree integrity check ahead of V verification.)

- [ ] **X1 state/config physical split**: execute per the F1 annotations. Two commits: first extract state.f90
  (the thing-in-itself 5 groups + dimension parameters; mechanical replacement of the use sites across the
  codebase), regress; then consolidate config.f90 (the representation blocks), regress. venus_params' C1..C8/PI
  are kept (unit contract, L1 constants).
- [ ] **X2 driver finalized**: VENUS.f90 → driver.f90 (the three-phase call sequence); harmonic_sampling
  registered as dist_surface_oscillator (the B′ third implementation, activating the currently dead compiled
  unit).
- [ ] **X3 sampler member table**: SELECT's NACTA/NACTB integer dispatch becomes dist_* member-table
  registration. Roster: dist_boltzmann(THRMAN), dist_ebk(INITEBK), dist_surface_oscillator,
  dist_rotation(ROTEN), dist_j(PROBJ+JMAXCALC), dist_incident(SURF+DENS), dist_barrier_excitation(BAREXC),
  dist_thermostat(THERMO), dist_glo_target(GLOSELECT), dist_thermalize(md_equilibrate),
  dist_normalmode(INITQP/HOMOQP/ORThAN + NMODE/FMTRX/FGMTRX[±ENMODE per adjudication]), dist_qmmicro(QMMICRO).
  Incidentally: parameterize the +100 Å magic number in SELECT :1139-1143 md_equilibrate (default unchanged,
  behavior unchanged).
- [ ] **X4 propagator.f90**: the integrators unified under an init/step/output-callback protocol; the
  verlet/symple/radau three members registered (orchestration already stripped by G4).
- [ ] **X5 L5 placement**: linalg.f90(EIGN), geometry.f90(ROTATE family/ROTATENO/ROTATEJKM/ROTATEX/Y/Z/CENMAS/
  ANGVEL/ROTN/ALIGN_ZSTEP/ALIGN_YSTEP/rotaline/GLPAR+PL); rng.f90 already in (C1); CPUSEC into the instrument
  slot alongside rng (VFDATE already dropped).
- [ ] **X6 analysis/ made offline**: the FINAL/GFINAL generalized COM/EREL/OAM parts (the outcome-classification
  executor, consuming E5 data), EIGOUT, DENQ+SUMQ (DEL=1 cm⁻¹ incidentally made configurable), HIST's five
  plot_* functions (lattice constants parameterized), the nearest-Au diagnostics (moved over from B1), and
  EIGOUT's IP=77/fort.77 channel-reuse split; conserv_check and smoke_rst_bvk brought into the build (scripts
  farewell to orphan command lines).
- [ ] **X7 build and docs**: Makefile SRC reordered by the L5→L1→L2→L3→L0 dependency order; EXTENDING.md
  rewritten (extension points = pes_* members/dist_* members/predicate primitives/outcome-classification rows/
  integrator registration).
- [ ] **Per-task acceptance**: `make` + regress all PASS (bit level). After X7, run the T2 drill: adding a
  `pes_leps` variant member (e.g. a reparameterized LEPSII) touches only a new file + a registration line; diff
  confirms zero touching of existing files.

## 8. Closing S

- [ ] **S1** Full execution and wrap-up of the new verification case system + per-case results.txt generation
  (produced by the new program, including records of criteria and tolerances).
- [ ] **S2** New-program completeness re-audit: review all routines against the spec's three questions,
  producing a new completeness audit (md+json); the old audit kept as historical reference.
- [ ] **S3** Document wrap-up: new-tree CLAUDE.md/README/EXTENDING.md (extension points = pes_* members/dist_*
  members/predicate primitives/outcome-classification rows/integrator registration); design_spec_completeness.md
  §4 blanks filled (interface A/B′ formal signatures + degenerate-point well-definedness checks +
  dual-implementation verification records: A axis RST+EMT-NN+GLO stochastic member+reference potentials, B′
  axis THRMAN+INITEBK+surface_oscillator, configuration twobody vs rst, electron axis adiabatic degenerate +
  TDHF-FSSH placeholder); optimal_structure.md status changed to "implemented (new program)"; root README layout
  updated (src/ new program + venus2026/ reference).
- [ ] **S4** Notion "VENUS-TDHF generalization" (Phase P1 architecture) task reconciliation: the checkbox
  states of this plan's phases synced.

## 9. Acceptance yardsticks and regression matrix

**Verification matrix** (third revision, 2026-09-06: old cases discarded, verification system rebuilt; 2026-09-05
addendum introducing the R8 literature benchmarks — external anchors take priority over self-set criteria):

| Time point | Scope | Criterion hierarchy |
|---|---|---|
| C1-C8 (per component) | each component's minimal verification input + criteria (criteria before implementation) | analytical anchors > statistical distributions (KS etc.) > conservation laws/closed quantities |
| V | the C-step verification inputs assembled into the formal case set; criterion coverage check (every registry-family member and every closed-family member covered by at least one case); full runs + same-ISEED reproducibility spot checks | same as above + reproducibility (two runs with the same ISEED bit-level identical) |
| S | full wrap-up + results.txt | same as above |

**R8 literature benchmark series** (selected 2026-09-05; criteria: mutually distinct aspects, reproducible, clear
numerical anchors, analytic potential functions preferred):

| ID | Benchmark work | Potential surface | Verification surface | Anchor data | Hook-up |
|---|---|---|---|---|---|
| R8-1 | F+H₂ QCT, Aguado et al., Chem. Phys. 187, 227 (1994) | modified LEPS (Takayanagi-Sato 1988, analytic) | gas-phase scenario + initial-state quantization + final-state analysis | σ_R(E_T, j=0-2) excitation-function table; DCS; rate constants | isomorphic with the LEPS reference potential, just parameter substitution; can start once C2 is done |
| R8-2 | H+H₂ QCT (Bohm-Truhlar M5 chain, LSTH/SLTH) | LSTH analytic (London form) | gas-phase scenario (optional, three identical atoms) | rate-constant table (author archive at comp.chem.umn.edu) | Truhlar group website download channel verified |
| R8-3 | GLO surface oscillator, Zhou & Jiang JCP 150, 024105 (2019) | analytic oscillator + friction (all parameters analytic) | force-contribution stochastic-law member (C2 stochastic-law acceptance) | T_s-dependent sticking curves | literature library id=83; GLO parameters of the same kind |
| R8-4 | H₂O/Ni(111) QCT, Jiang & Guo PRL 114, 166101 (2015) | PIP-NN (non-analytic) | surface scenario + outcome classification (methodology-level comparison, not strict reproduction) | η efficiency ratios 1.50-1.85/1.19-1.46/1.40-1.60; top-site anomaly; θ dependence 9×@573K | literature library id=318 |
| R8-5 | Tully 1990 three model problems, JCP 93, 1062 | 1D two-state analytic diabatic | electron-axis class-B member acceptance piece (first stop of the FSSH regression) | transition/transmission probability vs momentum curves (incl. the LZ analytic limit) | in parallel with C2 elec_* member development; MLatom 2024 (arXiv 2512.19152)/CRIT 2025 as modern comparisons |
| R8-6 | H₂/Cu(111) SRP-DFT (Leiden series) | SRP functional (non-analytic) | surface scenario + statistical accuracy (long-term) | molecular beam sticking curves, chemical accuracy | awaiting the ML-PES pipeline |

R8 implementation order: R8-1 → R8-5 (parallel with C2) → R8-3 (C2 stochastic-law acceptance) → R8-4
(methodology comparison) → R8-2/R8-6 (optional/long-term). Relation of R8 to the verification matrix: R8 is the
external-anchor source for the V-phase formal case set, each entry bound to the criterion design of the
corresponding C step.

Literature acquisition TODO: Tully 1990, Aguado 1994, Takayanagi-Sato 1988 (download into the library).

**Structural yardsticks** (S phase, item by item): T1-T7 (§1); invariant three-question spot checks (spec 2.6) —
sample one routine each from driver/force_facade/recorder/predicate/sampler and ask "does it only handle Q/P/s/F,
does it request services per contract, does it contain system branches/constants internally".

## 10. Risks and rollback

| Risk | Mitigation |
|---|---|
| Reference-code defect contagion (audit violations carried into the new tree) | the §6 stripping list checked item by item; check_* component-level cross-checks catch numerical deviations |
| Behavioral drift from rewriting from scratch (same algorithm, different operation order → floating-point differences) | criteria mainly analytical/statistical/conservation-based (insensitive to operation order); where bit-level precision is needed (e.g. conversion expressions identical in form), declared explicitly in the component plan; ad-hoc same-quantity comparison with the old program aids attribution |
| Skeleton at odds with reality (a signature found unsuitable in C) | skeleton drift ban: back to the review gate with the difference, revise the skeleton after confirmation; silently changing signatures forbidden |
| Old/new confusion (venus2026/ reference vs src/ new tree) | legacy read-only; the new tree independently built from P0 onward; every new-tree file header notes the reference provenance (legacy file:line) |
| Insufficient verification-criterion coverage (new risk after discarding the old cases) | criteria-first principle: each C step defines criteria before implementation; V-phase criterion coverage check (registry-family/closed-family members each covered); a gap means an added case |
| C5 sampler has the most members and the largest surface | member-table cross-check of the NACTA/NACTB dispatch table branch by branch; integrator_matrix and other cases cover multiple branches |
| Windows/Linux dual-platform build | the new-tree Makefile continuously buildable on both platforms from K3 onward; one run on each platform at the end of V (if feasible) |
| Rollback | each small phase gets an independent plan file + an independent commit group; `git revert` restores at single-phase granularity |

## 11. Deliverables list (completion markers)

1. The full standalone new-program tree (src/, L0-L5, T1) + the new verification case system fully passing (V/S1);
2. The new-program completeness audit report (S2);
3. The full document wrap-up set (S3): spec §4 filled in, new-tree CLAUDE/README/EXTENDING, the
   optimal_structure status update, root README layout;
4. Process trail: docs/plans/ with one plan file per small phase + git commit groups, fully traceable.
