# src/ - VENUS2026 completion rewrite

Completion-oriented rewrite of the VENUS chemical dynamics program: adiabatic QCT with L0-L5 layering (layer
concept in the table below;
directory names carry no layer numbers).

**Sole execution authority: [docs/final_completion_plan.md](../docs/final_completion_plan.md)**
(target end state §1, execution sequence §2, acceptance criteria §9). Architecture adjudication
standards are in docs/design_spec_completeness.md; target structure and material mapping in docs/optimal_structure.md.
The old library venus2026/ is a read-only reference (source of algorithms and formulas); it is neither migrated nor modified.

## Layout (entry points at the root; interface/ = the four interface modules, methods/ = the built-in
method library, systems/ = physical-system container material)

```
src/
  driver.f90  input.f90                          entry points: control-kernel orchestration (assembly -> evolution -> statistics) + input grammar
  globals/     consts  state  control  config  config_atoms
                                                  variables and provisions (the root everything else uses)
  interface/      force_interface  elec_interface  sysdef_interface  sampler
                                                  container_sched  spectrum_interface
                                                  the four interface modules (E14 four-entry contract:
                                                  force_eval / elec_prop / mqc_statistic + the
                                                  container_term slot; E13 purely formal authority;
                                                  E27 - the physical-system container delivers the
                                                  composed force evaluation; the ws interchange is
                                                  retired, a historical term) + system-definition
                                                  contract (discovery / parsing dispatch / list_atoms fill / paradigm
                                                  derivation / sampler hand-off) + sampling framework
                                                  + the container schedule (container_sched: eight
                                                  wiring slots + the initialization-assembly guard,
                                                  2026-09-22; the export slot withdrawn 2026-10-05)
                                                  + the member-seam spectrum producer
                                                  (spectrum_interface: the SPECTRUM_SOURCE-selected
                                                  derivation or manual table file, 2026-10-05)
  dynamics/    propagator                    closed family (three integrators; fixed membership -
                                                  E18: the predicate grammar interpreter withdrawn,
                                                  termination is the container-side container_term slot)
  methods/     elec_adiabatic  reg_adiabatic  samp_boltzmann  samp_ebk  inc_surface  inc_pair  inc_single  beam_laws  README.md
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
| L1 variables & provisions | globals/ (input stays at root) | physical-state holder, control flow, freezing-derived provisions, list_atoms |
| L2 interfaces & frameworks | interface/ | the E14 four-entry interface contract (force_eval / elec_prop / mqc_statistic + the container_term termination slot) + method-package registry (method_reg; E13 purely formal authority - name->pointer dispatch, physical-combination legality left to the members; rows may come from methods/ or a systems/ folder alike - E27) + run-lifetime statistics holders (occ / n_hop with member write-back) + system-definition contract sysdef_interface (E14) + sampling framework + the container schedule container_sched (eight wiring slots + the initialization-assembly guard; the container reg's single `reg_<sys>_slots` handover fills them - 2026-09-22; the registration-time export slot withdrawn 2026-10-05) + the member-seam spectrum producer (spectrum_interface: the SPECTRUM_SOURCE-selected internal derivation or manual spectrum-table file - 2026-10-05 spectrum-source legislation) |
| L3 closed families | dynamics/ | propagators only (fixed choices, not extension points; E18: the predicate grammar interpreter withdrawn - termination is the container-side container_term slot in force_interface) |
| L3 method library | methods/ | flat built-in library (E14, rebuilt from the two former interface_* homes): elec_* MQC method packages - parameter-free generic algorithms bound by the config method key (level 1) or consumed as reusable components via use (level 2) - the samp_* initial-state distributions, the inc_* incident-channel members (paradigm-split, paradigm-guarded registration) and their shared beam_laws component; a new method = one file + one thin reg_*.f90 (drag in and it compiles). System-specific members (pes_* PES + termination primitives) are system material, self-contained under systems/ (E12) |
| L4 recording & post-processing | output/ | trajectory serialization (column registry), final state + the container classification slot (E19)/density of states/plotting (offline-capable) |
| L5 instruments | utils/ | random stream (single), diagonalization, rotations |

Layering essentials: each layer depends only downward, and registry-family members never depend on each other
(members may call utils instruments internally);
the interface is all the physics the control kernel sees; unit conversion happens only at the force_interface boundary;
one single random stream for the whole library (rng.f90).

Call graph (E14, spec §8): the assembly phase chains read_input -> sysdef_load FIRST -> list_atoms_load ->
registration (methods `reg_*_all`; the container's `reg_*_slots` handover fills the
container_sched slots - the K3 prebuild's registry_gen_keys/_all/_init aggregates drive
keys-before-parse, registration and the initialization assembly) -> container_init - the config_atoms ordering statute (sysdef fills the atom list, so list_atoms_load shrinks
to allocation/algebra/validation plus the window and role reconciles - the A_LAT/BOXPAD surface-window keys
and the PROJECTILE/TARGET collision-role keys are paradigm-conditional there, their single write site).
The evolution loop body, explicit in the driver (five call entries):
prop_step -> elec_prop -> mqc_statistic -> rec_frame -> container_term
(force_eval runs inside the prop_step members per the propagator's endpoint-force contract, and optionally
inside mqc_statistic as an intra-box force refresh; the physical-system container composes its own active
force - the E27 delivery contract, no interchange workspace). Per trajectory the sampler runs a two-phase
dispatch: a DRAW phase (each DISTINCT selected member's draw once - same word on several fragments
= one call serving all carriers, equal treatment - then the orientation-stage matrices, then the
incident prescription) and a REALIZE phase (the member realizations in the same order, then the
orientation applies, then the incident assembly). Members split into draw (RNG -> module-private
record, no state access) and realize (record -> q/p, no RNG); the fused xxx_sample wrappers stay
as the check-program entries. Between the realizations and the orientation stage sits the bath loop
(interface/bath.f90; BATH = ANDERSEN arms it, default NONE is bitwise inert): ownerless atoms
(the slab) are seeded from list_atoms%q0, scheme-carrying fragments are frozen bitwise and
placed at a PES-validated non-interacting separation (utils/geom_init.f90 - the prober doubles
as the cheap PES sanity check), N_BATH steps of propagation + one ownerless-scoped Andersen
sweep (methods/bath_andersen.f90 - thermostat left the DIST_SCHEME registry for this) run,
then one whole-system drift removal and the bitwise unfreeze. The
orientation stage (RANDOM_ORIENT = T only, default F keeps the archived stream) orients every
nat>=2 fragment whose scheme word is not 'rotation', and must
sit before it: after it the fragments carry the beam momenta and a rotation would bend the collision
prescription; rot_euler writes the fragment back COM-re-centered with zero net momentum - incident
re-places the COMs and members own zero net momentum by contract).

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
