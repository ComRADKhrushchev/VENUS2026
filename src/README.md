# src/ - VENUS2026 completion rewrite

Completion-oriented rewrite of the VENUS chemical dynamics program: adiabatic QCT with L0-L5 layering (layer
concept in the table below;
directory names carry no layer numbers).

**Sole execution authority: [docs/final_completion_plan.md](../docs/final_completion_plan.md)**
(target end state §1, execution sequence §2, acceptance criteria §9). Architecture adjudication
standards are in docs/design_spec_completeness.md; target structure and material mapping in docs/optimal_structure.md.
The old library venus2026/ is a read-only reference (source of algorithms and formulas); it is neither migrated nor modified.

## Layout (entry points at the root; facade/ = the physics container boundary, methods/ = the built-in
method library, systems/ = per-system container material)

```
src/
  driver.f90  input.f90                          entry points: control-kernel orchestration (assembly -> evolution -> statistics) + input grammar
  globals/     consts  state  control  config  config_roster
                                                  variables and provisions (the root everything else uses)
  facade/      force_facade  elec_facade  sysdef_facade  sampler
                                                  the container boundary (E14 four-entry contract:
                                                  force_eval / elec_prop / mqc_statistic + the
                                                  container_term slot; E13 purely formal authority;
                                                  E27 composition root - the container delivers the
                                                  composed force kernel, the ws interchange retired) + system-definition
                                                  contract (discovery / parsing dispatch / roster fill / paradigm
                                                  derivation / sampler hand-off) + sampling framework
  dynamics/    propagator                    closed family (three integrators; fixed membership -
                                                  E18: the predicate grammar interpreter withdrawn,
                                                  termination is the container-side container_term slot)
  methods/     elec_adiabatic  reg_adiabatic  dist_boltzmann  dist_ebk  dist_incident  README.md
                                                  flat built-in method library (rebuilt from the former
                                                  interface_electronic/ + interface_distributions/, E14):
                                                  elec_* method packages (one method_reg row each, filled by a
                                                  thin reg_*.f90) + dist_* initial-state distribution members;
                                                  README.md = the library charter (E27: system-specific
                                                  method packages live in their systems/ folder)
  output/      recorder  final_state  density  plot_hist
                                                  online recording (column registry) + offline post-processing
  utils/       rng  linalg  geometry  specio  hessian
                                                  instruments (single random stream, diagonalization, rotations,
                                                  POSCAR/.xyz system-definition readers, central-difference Hessian)
  Makefile  README.md  CLAUDE.md  EXTENDING.md
```

| Layer | Directory | Responsibility |
|---|---|---|
| L0 control kernel | root driver.f90 | assembly -> evolution -> statistics call sequence, no physics |
| L1 variables & provisions | globals/ (input stays at root) | physical-state container, control flow, freezing-derived provisions, roster |
| L2 facades & frameworks | facade/ | the E14 four-entry container contract (force_eval / elec_prop / mqc_statistic + the container_term termination slot) + method-package registry (method_reg; E13 purely formal authority - name->pointer dispatch, physical-combination legality left to the members; rows may come from methods/ or a systems/ folder alike - E27 composition root) + run-lifetime statistics holders (occ / n_hop with member write-back) + system-definition contract sysdef_facade (E14) + sampling framework |
| L3 closed families | dynamics/ | propagators only (fixed choices, not extension points; E18: the predicate grammar interpreter withdrawn - termination is the container-side container_term slot in force_facade) |
| L3 method library | methods/ | flat built-in library (E14, rebuilt from the two former interface_* homes): elec_* MQC method packages - parameter-free generic algorithms bound by the config method key (level 1) or consumed as reusable components via use (level 2) - and dist_* initial-state distributions; a new method = one file + one thin reg_*.f90 (drag in and it compiles). System-specific members (pes_* PES + termination primitives) are system material, self-contained under systems/ (E12) |
| L4 recording & post-processing | output/ | trajectory serialization (column registry), final state + the container classification slot (E19)/density of states/plotting (offline-capable) |
| L5 instruments | utils/ | random stream (single), diagonalization, rotations |

Layering essentials: each layer depends only downward, and registry-family members never depend on each other
(members may call utils instruments internally);
the facade is all the physics the control kernel sees; unit conversion happens only at the force_facade boundary;
one single random stream for the whole library (rng.f90).

Call graph (E14, spec §8): the assembly phase chains read_input -> sysdef_load FIRST -> roster_load ->
registration (prop_reg / reg_*_all / pes_reg; the K3 prebuild's registry_gen_all() replaces the hand-written
lines) -> container_init - the config_roster ordering statute (sysdef fills the roster, so roster_load shrinks
to allocation/algebra/validation; a roster_load ahead of the sysdef products is an assembly-order error).
The integration-loop body is the four-phase step, explicit in the driver:
prop_step -> elec_prop -> mqc_statistic -> rec_frame -> container_term
(force_eval runs inside the prop_step members per the propagator's endpoint-force contract, and optionally
inside mqc_statistic as an intra-box force refresh; the container composes its own active force - the E27
delivery contract, no interchange workspace).

## Build (K3 final form)

```bash
make              # gfortran (default) -> build/gfortran/venus.e
make FC=ifx       # ifx (the project's conventional flags) -> build/ifx/venus.e
make clean        # remove the whole build/ tree (both compiler subtrees)
```

- Single-entry Makefile: wildcard source pickup (src/*.f90, src/*/*.f90, systems/*/*.f90) +
  a prebuild step (`scripts/prebuild.sh`, runs on every make) that regenerates
  `build/<FC>/registry_gen.f90` (the aggregate registration module — registration is a BUILD
  ARTIFACT, the driver carries zero hand-written registration lines) and `depend.mk`
  (object-level use-dependencies). Dragging a new methods/ or systems/ file in is the whole
  installation act; the Makefile never changes. Known residual: a tree-structure change
  (new reg file) needs a SECOND make (depend.mk is one run stale).
- Artifacts land only in per-compiler `build/<FC>/` subtrees (.mod formats are mutually
  incompatible); a stale cwd-level .mod/obj in the source tree shadows the build dir and kills
  downstream compiles — never compile ad hoc into src/.
- The Windows ifx section (MSVC link.exe PATH shadowing, `-module:` colon form, LIB export) is
  built in; per-machine overrides: `make FC=ifx MSVC_BASE=... WINKIT_LIB=... INTEL_BASE=...`.
- Working guidance (cases, criteria programs, contracts): [CLAUDE.md](CLAUDE.md);
  extension points: [EXTENDING.md](EXTENDING.md).
