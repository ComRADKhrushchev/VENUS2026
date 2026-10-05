# AGENTS.md — VENUS2026 workspace guide

VENUS2026 is a completion-oriented Fortran rewrite of the VENUS chemical dynamics program
(adiabatic QCT). The active codebase is `src/`; `venus2026/` is the OLD minimal form kept as a
read-only reference — never edit or migrate it. Repo docs and commit messages are largely
Chinese; match the surrounding language when writing docs/comments.

## Authority documents (read before touching sensitive areas)

- `docs/final_completion_plan.md` — sole execution authority (target end state, sequence, acceptance).
- `docs/design_spec_completeness.md` — architecture adjudication basis.
- `docs/optimal_structure.md` — target structure and material mapping.
- `src/README.md` — layer table (L0–L5), layout, and the assembly call-chain statute.
- `docs/plans/` — per-phase task records (C1–C8, V, K1–K3...); history lives here, not in code comments.

## Build & run

All builds go through `src/Makefile` (single entry):

```
cd src
make              # gfortran (default) -> build/gfortran/venus.e
make FC=ifx       # ifx (project's conventional flags: -r8 -double-size=64 -i8 -O2 -w)
make clean        # removes the whole build/ tree (both compiler subtrees)
```

- Registration is a BUILD ARTIFACT: `scripts/prebuild.sh` regenerates `build/<FC>/registry_gen.f90`
  and `depend.mk` on every make. A new method = one file in `src/methods/` (or `systems/<name>/`)
  plus a thin `reg_*.f90`; the wildcard picks it up — never hand-edit the Makefile or registry.
  Naming contract (2026-09-22): a methods/ reg hosts `<mod>_all` (+ optional `<mod>_init` —
  the member-parameter load); a systems/ (container) reg hosts `<mod>_slots` only (the
  container-schedule handover — see the E-block below); the aggregates keep their three
  driver-facing names (`registry_gen_keys` / `_all` / `_init`).
- Container storage (2026-09-21 ruling): only ONE container folder sits active in
  `systems/` per build (one-container-per-process). Verified containers are stored in
  the repo-root `containers/` library (outside the build glob — nothing there compiles);
  work-in-progress containers park in `systems/inactive/` (one level deeper than the glob).
  Moving a folder back into `systems/` is the whole activation act. Their criterion
  programs' compile recipes point at `containers/`.
- Artifacts land only in `build/$(FC)/` (per-compiler trees; .mod formats are incompatible between
  compilers). Stale `src/*.mod` or `*.obj` in the source tree cause "Unexpected EOF" — never compile
  ad hoc into `src/`.
- Windows ifx quirks are already handled in the Makefile (MSVC link.exe path shadowing, `-module:`).
- Run: `cd cases/<name> && ../../src/build/<FC>/venus.e` — input is hardcoded `input_qct.txt`
  in the run directory; outputs are `fort.*` / `stdout*` (gitignored).
- Container shelf (user ruling 2026-09-21): all non-reference containers live in
  `systems/inactive/` (two levels deep — outside the build glob). `systems/` top level keeps
  only `reference/` (the production container of the five cases); to test a container, move
  its folder to `systems/` top level, `make` twice, test, then move it back. Concurrent
  containers double-bind at assembly (STOP 1) — the check_* criterion programs compile
  containers standalone and are unaffected; running venus.e on a tested container requires
  temporarily shelving reference too (or the still-parked selection mechanism).

## Architecture rules that matter for edits

- Layers L0–L5 depend only downward (see `src/README.md`). The control kernel (`driver.f90`)
  sees physics only through `interface/` (E14 four-entry contract: force_eval / elec_prop /
  mqc_statistic; termination via the container_term slot).
- **E27 (2026-09-16) — the physical-system container as composition point**: the container
  (`systems/<sys>/`) delivers the composed force evaluation (`container_pes_i(q, mass, g, v)`
  in eV, fully assigned); the former ws interchange and the h_form/h_deliver matchmaking are
  retired (historical terms — docs/plans/S.md, E27). System-specific
  electronic method packages live in the system folder and register from there
  (`method_reg(name, prop, match [, sync])`); statistics write-back via `occ_set`/`n_hop_set`.
- **E28 (2026-09-22) — container parameters ride the program input**: the run-directory
  `*_pars.txt` channel is RETIRED. A container's parameter vocabulary = DECLARED container keys
  of input_qct.txt: the box owns `*_keys()` + `*_params_t` (field defaults are the container
  defaults); the reg file's `declare_keys` body (under its `reg_<sys>_slots` handover — pure
  vocabulary data) declares them (`input_declare_keys`), and its `load_params` body pulls one
  line per key
  (`pull(key, target [, found])` — absent key keeps the default, `found` carries optional-slot
  presence) into the receiver, then `*_load(p)` takes it into module state and validates.
  `input_audit_declared` (driver, after `registry_gen_init`) fatals a buffered-but-unconsumed
  declared key. A declared key must not collide with a program/member key. Ruling record:
  docs/plans/T.md.
- **Container schedule (2026-09-22) — the reg file carries content, the schedule carries
  time**: a systems/ (container) reg file hosts exactly ONE public entry,
  `reg_<sys>_slots()` — the single handover that fills the nine wiring slots of
  `interface/container_sched` via one keyword call. The slot bodies below it are private and
  grouped in four fixed wiring groups — input (`declare_keys`, `load_params`), regular
  (`bind_force`, `register_method`, `export_spectrum`), initialization (`bind_term`,
  `arm_method`), output (`bind_classify`, `register_columns`) — a null slot is a legal skip.
  The schedule owns the lifecycle order (keys before the parse; force/method/export at
  registration; load/term/classify/arm/columns at the initialization assembly) and the
  seam guard ONCE (`system_dir` empty → no seam wiring runs); every gate (key found, method
  selected) lives inside the container's own slot bodies. Slots are single — a second fill
  aborts (one container per process). Ruling record: docs/plans/2026-09-22-container-schedule.md.
- **Reactant definitions live in the container folder (2026-09-22)**: each system-specific
  container carries its reactant configuration files in a `system/` subfolder (one file per
  reactant, .xyz / .poscar; the subfolder is scan-pure — the loud-stop scan must not see the
  code/engine files at the folder root; `SYSTEM_DIR` points here). Name sort = list_atoms order,
  which matches the aggregate entries' atom layouts (C+Slab → atom 1 = incident C then the
  slab; F+H2 → F, H, H). `system/Slab.xyz` (cau_emtnn / cau_rst) is dual-role: engine input
  (the `SLAB_FILE` key / check path constants) and the target reactant definition. The
  system-agnostic reference container carries none — the five cases keep their case-side
  `system/` folders.
- Registry-family members (methods/, systems/) never depend on each other; they may use utils/.
- **Collision roles & equal sampling (2026-09-21)**: `PROJECTILE`/`TARGET` name list_atoms
  fragments by file stem (required for gas pairs and the surface paradigm — no positional
  default; list_atoms_load's role reconcile is the single write site of `list_atoms%i_proj/i_targ`;
  single-fragment gas rejects both keys). A `DIST_SCHEME` word may select several fragments:
  the sampler dispatches each distinct member once (first list_atoms occurrence) and the member
  serves all its carriers with ONE parameter set (equal treatment; per-side parameter
  differentiation stays a registered future item).
- **Universal orientation stage (2026-09-22)**: `RANDOM_ORIENT = T` (default `F` — the archived
  five-case stream stays byte-identical) inserts, inside `samp_run` after the per-fragment members
  and before the incident channel, one `rot_euler` pass per nat≥2 fragment whose scheme word is not
  `rotation` (that member already draws its orientation at fixed J). `rot_euler` writes the fragment
  back COM-re-centered with zero net momentum — safe at that point because the incident channel
  re-places the COMs by prescription and members own zero net momentum by contract.
- **Sampler two-phase dispatch (2026-09-22, plan A)**: every distribution member splits into
  `xxx_draw()` (consumes the RNG into a module-private sample record, touches no state) and
  `xxx_realize(sta)` (writes q/p from the record, draws no RNG); registration is
  `samp_reg(name, draw, realize)`. `samp_run` runs one DRAW phase (member draws in list_atoms
  order, then the orientation-stage matrices, then the incident prescription) then one REALIZE
  phase (member realizations in the SAME order - a compatibility statute; reordering would change
  what the equilibration members' runs see - then orientation applies, then the incident
  assembly). The fused `xxx_sample` wrappers remain as the check-program entries. Because
  realizations draw no RNG, the single draw phase reproduces the historical fused stream exactly
  (verified bitwise under gfortran for the archived cases and the RANDOM_ORIENT = T path). Known
  cost: under ifx the records' module-memory passage shifts default FP contraction - trajectories
  differ from the pre-split ifx binary at ULP level (conservation statistics identical); the
  byte-identity statute holds under gfortran, the cases' archive compiler.
- **Incident split per paradigm (2026-09-28 ruling)**: the incident channel is three members
  under the `inc_` prefix - `inc_surface` (word `incident_surface`; exactly one PROJECTILE
  fragment, every other fragment is a SUBSTRATE carrier the beam assembly never re-places or
  rotates; drop site N_AIM 0/1/2 with `AIM_X`/`AIM_Y` [Angstrom, POSCAR frame] for mode 2;
  B_MAX = 0 is the legal on-site beam under N_B = 0; A_LAT must divide both in-plane cell
  edges; the projectile COM velocity is assigned in full - the surface is the momentum sink),
  `inc_pair` (word `incident_pair`; the gas-pair branch verbatim - two role-named fragments,
  reduced-mass split), and `inc_single` (word `incident_single`; the trivial no-op owner of
  the single-fragment gas case - the sampler dispatches the incident member UNIFORMLY).
  `samp_incident_word()` resolves the member from (surface_model, fragment count); the regs
  register under PARADIGM GUARDS so exactly one incident member is armed per process (the
  members= witness stays 11). Shared stochastic law lives in the level-2 component
  `methods/beam_laws.f90` (a beam_t instance per member - members never depend on each other).
  Same ruling: the per-fragment distribution family was renamed `dist_*` -> `samp_*`
  (registration words unchanged). Ruling record: docs/plans/2026-09-28-incident-split-surface-paradigm.md.
- **Bath provision (2026-09-22 ruling: thermostat left the dist family)**: THERMOSTAT is no
  longer a DIST_SCHEME word (the registry is ten words now); the Andersen species lives in
  `methods/bath_andersen.f90` and registers `bath_reg('andersen', andersen_sweep)` from
  `reg_bath_andersen.f90`. `BATH = NONE|ANDERSEN` (default NONE - bitwise inert) arms the
  pre-evolution bath loop owned by `interface/bath.f90`: inside `samp_run`, after the member
  realizations and BEFORE the orientation stage and the incident channel. An armed loop
  implies, in order: (1) implicit seeding of the OWNERLESS atoms (no roster fragment claims
  them - the slab of the surface paradigm) from `list_atoms%q0` with zero momenta; (2) the
  FREEZE STATUTE - every scheme-carrying fragment is member-owned: its sampled q/p is saved
  bitwise, the fragment is placed at a PES-validated non-interacting separation
  (`utils/geom_init.f90` geom_sep_probe: ladder scan until the fragment's perturbation of the
  environment force stays under tolerance - a never-decaying interaction aborts naming the
  suspect PES, the cheap PES sanity check), and restored bitwise after; (3) N_BATH steps of
  [prop_step at DT_BATH + one collision sweep over the ownerless atoms only]; (4) one exact
  whole-system COM-momentum projection, then the unfreeze. The loop's in-loop RNG draws are
  the one sanctioned exception to the realize-draws-no-RNG rule. Species keys T_BATH/NU_COLL
  are unchanged.
- Unit conversion happens only at the force_interface boundary; internal units Å/amu/10 fs.
- One single random stream for the whole library (`utils/rng.f90`) — never add another RNG.
- "Physical constants" (values fixed by recorded decision — e.g. the conversion pair
  23.0605 × 0.04184 = 0.96485132) are used verbatim in analytic checks, never re-derived
  from CODATA.
- Input grammar is KEYWORD=VALUE (`src/input.f90`); `GWRITE_LEVEL` bare default is 2,
  `RANDOM_ORIENT` bare default is F.

## Conventions & verification

- **Data reports must link their raw files (user ruling, re-established 2026-10-04)**: every
  number quoted in a report (ensemble statistics, drift readings, flip rates, distribution
  metrics) carries a clickable path to the raw file it was computed from — local files as
  markdown links to the absolute path, cloud-side artifacts as the full server path, archive
  references as the archived file path. A number without its source file is not reportable;
  derived tables cite both the extraction output file and the upstream run directory.
- **Session-start memory check (2026-09-22)**: interactive sessions get NO auto-injected
  cross-session context — before the first analysis or edit, query the claude-mem database
  once for this project's recent state (`mcp-search` `search` with `project="VENUS2026"` and a
  recent date window, or the mem-search skill's 3-layer workflow; `session_start_context`
  renders the worker-mode digest). Concurrent sessions on this tree are common and the memory
  database is the only channel that sees them. Pair the check with a fresh `git status` before
  editing anything whose map came from earlier turns or sessions (mid-session tree reorgs have
  happened; a stale path map is how collisions start).
- Terminology: the canonical register is `docs/glossary.md` (Chinese/English term tables, code-namespace
  registry, grep audits). Key rulings (2026-09-17): 物理系统容器 (physical-system container), 接口
  (the four interface modules), 受力 (the composed force evaluation), 演化循环 (evolution loop), 阶段
  (stage), 决定 (decision), 物理常数 (fixed-value physical constants — never re-derive from CODATA).
  Retired coinages (do not reintroduce): 停放/空转/组合根/锚/loud-stop, class-A/class-B
  (design-phase member-class codes — say force members / electronic method members instead;
  historical docs keep the original words) — write plain statements
  instead (runtime abort messages use the `fatal:` prefix). Cross-document code citations carry the
  owner prefix (C5-R1, check_facade (p), case_b_layer T2, venus_params C1). Each phase close includes
  a terminology-sync pass against the glossary.
- File headers: one-line purpose + short Design note (≤12 lines); subroutine headers 1–2 lines
  wrapped in `---`; numbered `! 1.` `! 2.` blocks for main logic. No historical E-number
  adjudication narratives, old-repo line references, or task-phase tags in source comments —
  those live in docs/plans (registered exception: criterion provenance notes inside check_*
  headers; the five E-number string literals in input.f90 / elec_interface.f90 error messages).
- `scripts/check_*.f90` are criterion/verification drivers (gfortran, self-contained). Verify
  changes with the relevant check program plus a clean dual build (`make` and `make FC=ifx`).
- Commit messages in this repo are detailed Chinese narratives describing evidence and
  verification legs; follow that style.
