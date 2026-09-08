# src/ - VENUS2026 completion rewrite

Completion-oriented rewrite of the VENUS chemical dynamics program: adiabatic QCT with L0-L5 layering (layer
concept in the table below;
directory names carry no layer numbers).

**Sole execution authority: [docs/final_completion_plan.md](../docs/final_completion_plan.md)**
(target end state §1, execution sequence §2, acceptance criteria §9). Architecture adjudication
standards are in docs/design_spec_completeness.md; target structure and material mapping in docs/optimal_structure.md.
The old library venus2026/ is a read-only reference (source of algorithms and formulas); it is neither migrated nor modified.

## Layout (entry points at the root; interface_ prefix = extensible registry-family members)

```
src/
  driver.f90  input.f90                          entry points: control-kernel orchestration (assembly -> evolution -> statistics) + input grammar
  globals/     consts  state  control  config  config_roster
                                                  variables and provisions (the root everything else uses)
  facade/      force_facade  elec_ws  sampler     single physics entry: force facade + electronic workspace
                                                  contract (ws container, E10) + sampling framework
  dynamics/    predicate  propagator             closed families (predicate grammar, three integrators; fixed membership)
  interface_electronic/    elec_adiabatic        home of the class-B electron-law axis (parameter-free
                                                  generic algorithms; class-A pes_* members are system
                                                  material - self-contained in systems/, E12)
  interface_distributions/ dist_boltzmann  dist_ebk  dist_incident
                                                  initial-state distribution members (B' extension; 3 examples now, the rest in C5)
  output/      recorder  final_state  outcome  density  plot_hist
                                                  online recording (column registry) + offline post-processing
  utils/       rng  linalg  geometry             instruments (single random stream, diagonalization, rotations)
  Makefile  README.md
```

| Layer | Directory | Responsibility |
|---|---|---|
| L0 control kernel | root driver.f90 | assembly -> evolution -> statistics call sequence, no physics |
| L1 variables & provisions | globals/ (input stays at root) | physical-state container, control flow, freezing-derived provisions, roster |
| L2 facades & frameworks | facade/ | (Q,P,s)->(F,E) facade (E10 single class-A signature pes_proc_i(q,mass,ws) - representation demoted to a ws fill mode; E13 purely formal authority - name->pointer dispatch, physical-combination legality left to the members) + electronic workspace contract elec_ws (ws_t container); sampling framework |
| L3 closed families | dynamics/ | predicate grammar, propagators (fixed choices, not extension points) |
| L3 registry-family members | interface_*/ | two families: electronic (class-B elec_* electron propagation laws - parameter-free generic algorithms; the class-A pes_* axis registers through the same facade but its members are system material, self-contained under systems/ - E12) / distributions (B' initial-state distributions) - a new file plus a one-line registration plugs a member in |
| L4 recording & post-processing | output/ | trajectory serialization (column registry), final state/outcome/density of states/plotting (offline-capable) |
| L5 instruments | utils/ | random stream (single), diagonalization, rotations |

Layering essentials: each layer depends only downward, and registry-family members never depend on each other
(members may call utils instruments internally);
the facade is all the physics the control kernel sees; unit conversion happens only at the force_facade boundary;
one single random stream for the whole library (rng.f90).

## Build (P0 skeleton)

```bash
make            # one ifx invocation compiles + links, producing venus.e (effective once SRC is filled in)
make clean      # remove venus.e and .mod/.obj/.pdb/.ilk
```

- Compiler: Intel `ifx` (validated on oneAPI 2025.2); flags `-r8 -double-size=64 -i8 -O2 -w`.
- `SRC` is an empty list for now; phase K3 (independent build of the new tree) fills in the fixed
  skeleton list in dependency order (utils -> globals -> facade -> dynamics -> interface_* -> output ->
  input -> driver) plus a wildcard pickup of systems/ (E11 "drag-in" ruling: each system blackbox sits
  in its own folder with its self-contained members + a reg_<system>_all module - E12; a prebuild
  step generates registry_gen.f90 from them, so dragging a folder in is enough - zero edits to
  Makefile or driver).
- The new program has no MKL dependency; the Windows linker-shadowing handling (`-Qlocation,link`, `LIB` export,
  `MSVC_BASE`/`WINKIT_LIB`/`INTEL_BASE` override mechanism) is ported from the old library's Makefile platform block,
  see the Build section of venus2026/CLAUDE.md for usage.
