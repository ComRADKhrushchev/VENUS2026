# VENUS2026 Optimal Program Structure and Material Mapping

Status: design draft (derived from the completeness specification, unconstrained by the existing code organization)
Date: 2026-09-05
Basis: design_spec_completeness.md (the specification), completeness_audit.md (the material inventory, 99 entries)
Method: first derive the target structure from the standard, then place the existing subroutines as material. Placement outcomes fall into four categories: placeable (as-is or with light modification), placeable after refactoring, no slot (eliminated), structure lacking material (must be newly written).

---

## 1. Target Structure (Derived from the Specification)

### 1.1 Layering

```
L0 control kernel      driver.f90          call sequence assembly→evolution→statistics, no physics
L1 thing-in-itself container  state.f90    Q, P, s, F, t (allocatable, dimension = roster content parameter)
L1 representation layer  config.f90        input-grammar product: scenario, roster, distributions, predicates, outcome classification, observable configuration
L2 Interface A facade  force_facade.f90    (Q,P,s) → (F,E); internally assembles Σ force-contribution members × single-choice electron propagation law
L2 Interface B′ framework  sampler.f90     distribution → (Q0,P0,s0) generic sampler
L3 registry-family members  pes_*.f90 / dist_*.f90   swappable implementation bodies
L3 closed families     propagator.f90 / predicate.f90    limited choices + grammar interpreter
L4 measurement         recorder.f90        trajectory serialization + predicate evaluation (all that resides in the evolution phase)
L4 post-processing     analysis/           outcome-classification execution, statistics, plotting (can be offline)
L5 instruments         rng.f90, linalg.f90, geometry.f90  random stream, diagonalization, rotations
```

Key points:
- each layer depends only downward, and L3 members do not depend on each other (registry-family mechanics);
- the L2 facade is all the physics the control kernel sees — swapping the PES, the method, or the scenario never touches L0/L4;
- L5 instruments form an independent layer, correcting inversions of hierarchy of the "gasdev lives in a PES file" kind;
- the single venus_data global block is dissolved: the thing-in-itself goes into L1 state, representations into L1 config, instrument state (the RNG stream) into L5.

### 1.2 Formal Definition of Each Layer (Corresponding Specification Clauses)

| Layer | Specification identity | Signature | Mutability |
|---|---|---|---|
| driver | control kernel | orchestration | invariant |
| state | thing-in-itself | data | invariant (form) |
| config | representation | data | invariant (grammar extends monotonically) |
| force_facade | Interface A facade | (Q,P,s)→(F,E) | invariant (framework code) |
| pes_* | Interface A class-A members | Q→(V,∇V), additive | registration-open |
| elec_* | Interface A class-B members | (H,s)→(ds/dt, force rule), single choice | registration-open (adiabatic = degenerate member) |
| dist_* | Interface B′ members | representation→distribution | registration-open |
| sampler | operation-class framework | distribution→(Q0,P0) | invariant |
| propagator | closed family | force functional→step | closed (adding a member = one registration line) |
| predicate | closed grammar | state→Boolean | invariant (interpreter) |
| recorder | measurement | thing-in-itself stream→records | invariant |
| analysis | post-processing | records→statistics | append law |

---

## 2. Material Mapping Master Table

Disposition symbols: ● = placeable (as-is / light modification); ◐ = placeable after refactoring; ✕ = no slot, eliminated; △ = structure lacking material (must be newly written; existing reusable material cited in parentheses).
Source line numbers omitted; see completeness_audit.md.

### L0 control kernel

| Material | Disposition | Notes |
|---|---|---|
| VENUS.f90 main-loop skeleton | ◐ | after stripping the admixtures listed in the E1-E4 inventory, the remaining orchestration is the driver embryo |
| (E1 outcome-block duplicated copy) | ◐ | once merged, becomes the statistics-phase entry point |

### L1 state / config

| Material | Disposition | Notes |
|---|---|---|
| venus_data Q/Q_old/P/PDOT/QDOT block | ◐ | migrates into state.f90; the other ~40 representation blocks all migrate to config |
| venus_params C1..C8/PI | ● | unit contract retained (L1 constants) |
| venus_params NN_POT_TYPES index | ◐ | migrates inside the pes_emt_nn member |
| input_parser / param_mapping / presets | ● | config's grammar layer retained as-is (the five mapping tables = closed-family registration lines) |
| venus_input reading skeleton | ◐ | the geometry/roster/ghost assembly segments are handed over to the configuration assembler (see △3) |
| △ roster configuration assembler | △ | missing. Material: the venus_input geometry segment + the test_potentials back-write segment (BOXLX/SKEW/NATOMB/QZB) + SURF surface preparation + FIXROTDATM alignment |

### L2 Interface A

| Material | Disposition | Notes |
|---|---|---|
| DVDQ_1 force dispatch | ◐ | its "zero PDOT → call DPESHON" logic is absorbed into force_facade |
| interface_TEST's POTPRE/POT0/DPESHON | ◐ | the POTPRE-style initialization hook becomes the facade protocol; the ×23.0605, ×C1 conversions stay at the boundary (the unit-contract fulfillment point) |
| test_potentials: HARMONIC/MORSE/LEPS | ● | pes_reference members (system-independent reference potentials, for dual-implementation testing) |
| emt_nn_pes_venus | ● | pes_emt_nn member (self-contained + externalized data, already the correct shape) |
| rst_pes + rst_slab_bvk | ● | pes_rst member (with PBC, its own projection logic) |
| (GLO spring force, currently embedded in VERLET) | ◐ | pes_glo_spring member |
| (ghost Langevin, currently embedded in VERLET) | ◐ | pes_langevin_ghost member (stochastic law) |
| FRICFORCE+NFC | ◐ | pes_ldfa member (revived) |
| ENERGY_1's NSURF/NGLO dispatch and GLO potential energy | ◐ | the energy side is absorbed into the facade (E = T(P) + ΣV is summed by the facade) |
| △ force_facade framework | △ | missing. Material: the dispatch shells of DVDQ + interface_TEST; to be newly written = member-table data structure, additive composition, capability-level negotiation |
| △ elec_adiabatic degenerate member | △ | missing (trivial: ds/dt≡0, identity force rule, but must exist explicitly — well-definedness at the degenerate point) |

### L2 Interface B′

| Material | Disposition | Notes |
|---|---|---|
| SELECT sampler framework | ◐ | once the NACTA/NACTB integer dispatch becomes a dist_* member table, it is the sampler |
| THRMAN | ● | dist_boltzmann (one of the dual implementations) |
| INITEBK | ● | dist_ebk (the second of the dual implementations) |
| harmonic_sampling | ◐ | dist_surface_oscillator (alive upon registration, instead of deletion) |
| ROTEN | ● | dist_rotation (the three NROT laws) |
| PROBJ + JMAXCALC | ◐ | dist_j; to be merged once the H2/D2 nuclear-spin tables are turned into data |
| SURF (incident channel) + DENS | ◐ | dist_incident (relative-motion distribution) |
| BAREXC | ● | dist_barrier_excitation |
| THERMO | ◐ | dist_thermostat (Andersen / velocity-scaling laws; assembly-phase use only) |
| ROTATE/ROTATENO/ROTATEJKM/ROTATEX/Y/Z | ● | L5 geometry (orientation implementations, called by the various dist_*) |
| HOMOQP/INITQP/ORThAN/QMMICRO | ◐ | normal-mode coordinate assembly is merged into dist_normalmode / dist_qmmicro; the NACTA=4 dead branch is dropped in the move |
| GLOSELECT | ● | dist_glo_target (GLO target preparation) |
| md_equilibrate (inside SELECT) | ● | dist_thermalize (the member internally calls evolution services, permitted by boundary 4 of 2.3) |
| NMODE/FMTRX/FGMTRX/ENMODE | ◐ | assembly-phase analysis tools, placed alongside dist_normalmode (the GAUSS dead branch is dropped) |

### L3 closed families

| Material | Disposition | Notes |
|---|---|---|
| VERLET | ◐ | after stripping the embedded GLO/Langevin = the closed propagator verlet |
| SYMPLE / RADAU | ◐ | RADAU joins after RADOUT's orchestration is stripped (switched to callbacks) |
| TEST's generic COM-criterion segment | ◐ | raw material for the and/or/com_distance primitives of the predicate grammar |
| △ predicate.f90 grammar interpreter | △ | missing. Material: TEST:106-146; to be newly written = grammar data structure (threshold comparison / Boolean combination / bounded window) + evaluator |
| △ outcome-classification configuration | △ | missing (currently a hard-coded if chain inside the VENUS.f90 outcome block). Material: the classification logic of the merged E1, turned into data |

### L4 measurement / post-processing

| Material | Disposition | Notes |
|---|---|---|
| GWRITE | ◐ | after stripping four sites of system knowledge = the recorder (serializer); GWRITE_LEVEL is retained as a configuration entry |
| record_all/HIST recording columns | ◐ | merged into the recorder's record table (column selection made configurable) |
| FINAL/GFINAL | ◐ | migrates to analysis/; GFINAL:80's EREL override is dropped; the generalized COM/EREL/OAM parts are the outcome-classification executor material |
| EIGOUT/DENQ+SUMQ | ● | analysis (density-of-states / eigen output) |
| HIST's five plot_* functions | ◐ | offline plotting in analysis (the generic surfcell/Dist_cell implementations are retained; the lattice constants become passed parameters) |
| conserv_check / smoke_rst_bvk | ● | analysis acceptance tools retained as-is |

### L5 instruments

| Material | Disposition | Notes |
|---|---|---|
| RAND0/RAND1/RANDST | ● | the rng.f90 main stream (ISEED3 state) |
| gasdev/GAMA | ◐ | migrates into rng.f90 and switches to the main stream (single random stream) |
| EIGN | ● | linalg.f90 |
| GLPAR+PL / rotaline | ● | mathematical tools |
| CPUSEC | ● | timing (VFDATE dropped — a SYSTEM shell call) |
| CENMAS/ANGVEL/ROTN | ● | geometry/kinematics |
| ALIGN_ZSTEP/ALIGN_YSTEP | ● | geometry |

### No Slot, Eliminated (Consistent with the Audit's Redundancy List)

The RK4 coefficients and dead entries, ADAMS (already loud-stop), MICROCI+VOLPSCONE, LMODE+LMEXCT+EBOND, ENMODE (MPLOT dead), POTEN, POTENZ ghost printing, THERMBATH, GLOEQU, DISTRIBUTE, PBC, LINEQ/LUDCMP/LUBKSB/FEHLER, WEBOND/WLBOND/WENMOD, ROTATEX_b/ROTATEZ_b, grandom/ran1/gasdev_0, the GLO.f stranded copy, the three dead blocks of VENUS.f90, the venus_data dead-block cluster, the MKL linkage, and the debug-print cluster (SELECT/SURF/ROTen/THERMO/RADAU/FINLNJ leftovers).

---

## 3. Gap List (Material the Structure Needs but the Library Lacks)

Ordered by necessity:

1. **force_facade framework** — the facade is the body of Interface A; the current library has only dispatch embryos (DVDQ+POTPRE). New code: member registry (name→procedure array), additive-composition loop, capability-level negotiation (H representation × propagation law × integrator compatibility triple check). Roughly 200-300 lines.
2. **predicate grammar interpreter** — termination predicates go from hard-coded ifs to data-driven evaluation. Material: TEST:106-146. New code: grammar parsing (predicate expressions in the input file) + evaluator. Roughly 150 lines.
3. **roster configuration assembler** — the single assembly point for roster / geometry / element tables. The material is scattered across venus_input/SURF/test_potentials. New code: roster data structure + element tables. Roughly 150 lines.
4. **recorder observable-configuration table** — GWRITE_LEVEL goes from "hard-coded levels" to "user-selectable columns". New code: a column-registration mechanism. Roughly 50-100 lines.
5. **elec_adiabatic degenerate member** — trivial but must exist explicitly (about 20 lines), reserving the place for TDHF-FSSH's future return as a class-B member.
6. **outcome-classification data** — after the E1 merge, turn the if chain into a table lookup (about 50 lines).

## 4. Structure Self-Check (Against the Specification)

- Phase restriction: the driver's three-phase call sequence; during the evolution phase only facade/propagator/recorder/predicate;
- Factorization: any pes_* pairs freely with any elec_* (currently only the degenerate member); no bundled routines;
- Dual implementations present: A axis HARMONIC+MORSE+LEPS+EMT-NN+RST (≥2); B′ axis THRMAN+INITEBK (+surface_oscillator, ≥3 after registration); configuration level, the two paradigms of molecular/surface (twobody vs rst_* cases);
- Unit contract: conversions only at the facade boundary (×23.0605, ×C1); inside state, always code units;
- Single random stream: rng.f90; the whole ISEED stream is reproducible.

## 5. Migration-Strategy Implications

The material mapping shows: the eliminated items coincide with the R0 list; most ◐ items correspond to steps R1-R6; the △ gaps are pure additions (no old code touched — they can be developed in parallel with the purification and verified by bypass). The recommended implementation order is unchanged (R0→R7), but △1/△2 can start early — they are the foundation of the new structure and depend on no purification step being completed.
