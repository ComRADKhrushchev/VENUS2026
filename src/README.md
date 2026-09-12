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
  facade/      force_facade  elec_facade  sysdef_facade  elec_ws  sampler
                                                  the container boundary (E14 four-entry contract:
                                                  force_eval / elec_prop / mqc_statistic, termination composed
                                                  through predicate; E13 purely formal authority, E10 ws container
                                                  held one integration step in elec_facade) + system-definition
                                                  contract (discovery / parsing dispatch / roster fill / paradigm
                                                  derivation / sampler hand-off) + sampling framework
  dynamics/    predicate  propagator             closed families (predicate grammar - four built-in leaves +
                                                  the open registered primitive set, three integrators; fixed membership)
  methods/     elec_adiabatic  reg_adiabatic  dist_boltzmann  dist_ebk  dist_incident  README.md
                                                  flat built-in method library (rebuilt from the former
                                                  interface_electronic/ + interface_distributions/, E14):
                                                  elec_* method packages (one method_reg row each, filled by a
                                                  thin reg_*.f90) + dist_* initial-state distribution members;
                                                  README.md = the library charter
  output/      recorder  final_state  outcome  density  plot_hist
                                                  online recording (column registry) + offline post-processing
  utils/       rng  linalg  geometry  specio     instruments (single random stream, diagonalization, rotations,
                                                  POSCAR/.xyz system-definition readers)
  Makefile  README.md
```

| Layer | Directory | Responsibility |
|---|---|---|
| L0 control kernel | root driver.f90 | assembly -> evolution -> statistics call sequence, no physics |
| L1 variables & provisions | globals/ (input stays at root) | physical-state container, control flow, freezing-derived provisions, roster |
| L2 facades & frameworks | facade/ | the E14 four-entry container contract (force_eval / elec_prop / mqc_statistic; termination two-track through predicate) + method-package registry (method_reg; E13 purely formal authority - name->pointer dispatch, physical-combination legality left to the members) + ws holder (E10 container, lifetime one integration step) + system-definition contract sysdef_facade (E14) + sampling framework |
| L3 closed families | dynamics/ | propagators only (fixed choices, not extension points; E18: the predicate grammar interpreter withdrawn - termination is the container-side container_term slot in force_facade) |
| L3 method library | methods/ | flat built-in library (E14, rebuilt from the two former interface_* homes): elec_* MQC method packages - parameter-free generic algorithms bound by the config method key (level 1) or consumed as reusable components via use (level 2) - and dist_* initial-state distributions; a new method = one file + one thin reg_*.f90 (drag in and it compiles). System-specific members (pes_* PES + termination primitives) are system material, self-contained under systems/ (E12) |
| L4 recording & post-processing | output/ | trajectory serialization (column registry), final state/outcome/density of states/plotting (offline-capable) |
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
inside mqc_statistic as an intra-box force refresh; ws is filled at force evaluation, consumed by
elec_prop / mqc_statistic, and reset at step end).

## Build (P0 skeleton)

```bash
make            # one ifx invocation compiles + links, producing venus.e (effective once SRC is filled in)
make clean      # remove venus.e and .mod/.obj/.pdb/.ilk
```

- Compiler: Intel `ifx` (validated on oneAPI 2025.2); flags `-r8 -double-size=64 -i8 -O2 -w`.
- `SRC` is an empty list for now; phase K3 (independent build of the new tree) fills in the fixed
  skeleton list in dependency order (utils -> globals -> facade -> dynamics -> methods -> output ->
  input -> driver) plus a wildcard pickup of systems/ (E11 "drag-in" ruling: each system container sits
  in its own folder with its self-contained members + a reg_<system>_all module - E12; a prebuild
  step aggregates every reg_*.f90 under systems/ and methods/ into registry_gen.f90 (E11 generalized
  by E14), so dragging a folder in is enough - zero edits to
  Makefile or driver).
- The new program has no MKL dependency; the Windows linker-shadowing handling (`-Qlocation,link`, `LIB` export,
  `MSVC_BASE`/`WINKIT_LIB`/`INTEL_BASE` override mechanism) is ported from the old library's Makefile platform block,
  see the Build section of venus2026/CLAUDE.md for usage.
