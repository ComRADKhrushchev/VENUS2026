# VENUS2026 Subroutine Completeness Audit Report

- Date: 2026-09-04
- Adjudication standard: `docs/design_spec_completeness.md` (the thing-in-itself/representation dichotomy, four classes of pure subroutines, control kernel, interfaces A/B′, closed layer, three-valued verdicts)
- Audit method: read-only static review. Every routine was read through across the root sources (10 files), src_VENUS (61 files), src_TEST (6 files), and 2 scripts tools; for every suspected dead entry point, all call sites in the codebase were cross-checked (grep-verified).
- Machine-readable master table: `docs/completeness_audit.json` (99 entries).
- Verdict vocabulary: **invariant class** / **interface or closed-family member** / **redundancy candidate or awaiting purification**.
- Category vocabulary: evolution (thing-in-itself→thing-in-itself), measurement (thing-in-itself→representation), operation (representation→thing-in-itself), post-processing (representation→representation), control kernel, infrastructure (representation-layer infrastructure), instrument, mathematics.

---

## 0. Comparison against known calibration points (conclusions up front)

| Calibration point | Audit conclusion |
|---|---|
| GWRITE.f90 historically hardcodes atom3=C | **Not clean**. The current form has evolved into **atom1=C** (`Q(C)=`/`V(C)=` at `GWRITE.f90:59-60`), plus two new spots of system knowledge: the mass-threshold element labeling `W(i)<50 amu → 'C '` else `'Au '` (:110-114), and the `Au_nearest`/`Au_xyz` nearest-neighbor diagnostics (:95-99). The violation went from "hardcoded atom 3" to "hardcoded atom 1 + mass-threshold element names"; the nature is unchanged. |
| RK4 dead entry, broken ADAMS → redundancy candidates | **Confirmed**. RK4: `param_mapping.f90` still accepts `RK4` (→INTEGRATOR=1, LLL=0), but `venus_input.f90:363-365` rejects LLL=0 → a dead entry that stops on read-in; the RKUTTA coefficients RAA1..RC2 computed at `VENUS.f90:86-94` have no consumers anywhere. ADAMS: INTEGRATOR=0 is explicitly loud-stopped at `venus_input.f90:339-343` — upgraded from "broken" to "explicitly rejected"; the elimination conclusion stands. |
| TEST.f termination criteria → should be closed predicate-grammar members | **Confirmed**. `TEST.f:88-94` has **live** system-specific branches (when NSURF>0, scattering judged by the hardcoded atom-1 height Q(3), the NC>50 magic number, P(3)>0 outward velocity — C/Au(111) criteria); :54-78 (F+H2O) and :149-192 (CH4/Ni(111)) are two commented-out blocks of dead system knowledge; the generic COM distance + recession-velocity criteria at :106-146 are the correct prototype of a closed grammar member. |
| No compile-time PES variants (TEST_PES runtime selection) → already conforms to interface A | **Confirmed**. The `interface_TEST.f90` header comment explicitly declares single-build multi-potential; the five members HARMONIC/MORSE/LEPS/EMT-NN/RST are registered at runtime via the `TEST_PES` keyword, with no compile-time variant proliferation. |

The conclusions on all four calibration points agree with the audit; there are no discrepancies needing explanation.

---

## 1. Master table (full routine coverage)

Verdict column abbreviations: inv=invariant class; mem=interface/closed-family member; red=redundancy candidate/awaiting purification. Violations are given as file:line.

### 1.1 Root sources (10 files)

| Routine | File | Four-class membership | Verdict | Violations / notes |
|---|---|---|---|---|
| VENUS (program) | VENUS.f90 | control kernel | red | :86-94 RKUTTA coefficients computed dead; :182-232 NCHKP dead branch (no input source); :184,231 RANDST double call (F15); :286-299 histogram dead block (idplt=0 hardcoded) + the site constants Qr1=2.94998/Qr2=2.55476 + SYSTEM shell call; :322-335 GLO layout branch; :340 surface-alignment call; :399-400,454-458 BOUNCE/z_pos_min observable accumulation hardcoded to atom 3; :462-464 empty check statements; :544-548,631-635,671-677 AE* observation recording embedded; :506-576 and :581-662 entire outcome blocks duplicated |
| venus_params | venus_params.f90 | infrastructure (unit contract/constants) | red | :88-98 NN_POT_TYPES=7 and the E0/1D/3P/ion-state indices — C+Au-family NN potential knowledge hardcoded in a constants module; C1..C8/PI are the legitimate unit contract |
| venus_data + allocate_venus_data | venus_data.f90 | global state | red | Mixed occupancy: the thing-in-itself (Q/P/PDOT/QDOT/Q_old) and ~40 representation blocks in one module with no separating comments; :585 auro=2.94998 (Au lattice constant); :382-384 HIST_ATOM_NAD with no writer; :389-391,417-418 MORSEB/STRETB dead data; dead blocks such as CUBEB/COORS/GDD/GAUSS2/U0-U1-V01/RKUTTA |
| read_venus_input | venus_input.f90 | operation (representation read-in and assembly) | red | :616-621 ghost Langevin assembly hardcodes W(2) and the 413.4 conversion; :703-704 NATOMB domain hardcoded 4/10; :997-998 BOXLX=32.5/BOXLY=30.7 default geometry; :752 W=1D30 rigid mass; :692 MPLOT=0 hardcoded (condemning the entire mode-population pathway to death); :392 NCVIB dead keyword. The read-in grammar itself is fine |
| input_parser | input_parser.f90 | infrastructure (representation grammar) | inv | Pure grammar (keywords/continuations/case/constructors); old formats loud-stop |
| param_mapping | param_mapping.f90 | infrastructure (enumeration registry) | inv | Five string→integer mapping tables = closed-family registration lines; removed items loud-stop |
| presets | presets.f90 | infrastructure (freezing presets) | inv | Five set_default presets made data — the correct direction of "moving system prescriptions from code into data" |
| sample_surface_oscillator | harmonic_sampling.f90 | operation (B′ member material) | red | No call site in the codebase (compiled but dead); should be wired into interface B′ or deleted |
| ENERGY_1 | ENERGY.f | measurement (T/V/H functionals) | red | :16-49 NSURF/NGLO dispatch; :40-47 GLO spring potential embedded (hardcoded NATOMA+1/+2) — interface A content out of bounds; POT_SH dead call in comments |
| DVDQ_1 | DVDQ.f | evolution (force-evaluation dispatch) | inv | Zeroes PDOT→calls DPESHON; the sole consumption point of interface A, treating the PES as a black box |

### 1.2 src_TEST (6 files)

| Routine | File | Four-class membership | Verdict | Violations / notes |
|---|---|---|---|---|
| POTPRE | interface_TEST.f90 | evolution (interface A registration) | inv | Runtime TEST_PES single-build multi-potential — conforms to interface A |
| POT0 | interface_TEST.f90 | evolution (interface A boundary) | inv | eV→kcal/mol (×23.0605) conversion centralized at the boundary (unit contract honored); a side effect writes E0 for diagnostics (minor) |
| DPESHON | interface_TEST.f90 | evolution (interface A boundary) | inv | −grad×23.0605×C1 accumulated into PDOT; the signature acknowledges additive contributions |
| FIXROTDATM | interface_TEST.f90 | operation (orientation alignment) | mem | :58-72 hardcodes atoms 1-2 along z; parameterized, it becomes a closed member |
| ALIGN_ZSTEP / ALIGN_YSTEP | interface_TEST.f90 | operation (kinematics) | inv | Rigid-body rotation stepping |
| gasdev | interface_TEST.f90 | instrument (normal random law) | red | :112-136 a core instrument living in a PES interface file (layer inversion, surviving only via Makefile order sensitivity); runs on the RANDOM_NUMBER stream (coexisting with the RAND0 stream; ISEED irreproducible). Live code |
| test_potentials (whole module family) | test_potentials.f90 | evolution (interface A member family) | mem | :152-158,287-293 members writing back representation (BOXLX/BOXLY/SKEW); :163,260 members writing the roster NATOMB; :210-213,274-281 members writing assembly geometry QZB. HARMONIC/MORSE/LEPS are system-independent reference members; EMT-NN/RST carrying C/Au(111) is legal inside members |
| emt_nn_pes_venus | emt_nn_pes.f90 | evolution (interface A member) | mem | Self-contained engine (does not use venus_data), data externalized under data/emt_nn/ with an overridable path keyword — the correct shape of a member carrying external data |
| rst_pes | rst_pes.f90 | evolution (interface A member auxiliary) | mem | Zero-dependency module; FCC(111) slab generation/weights IO |
| rst_slab_bvk | rst_slab_bvk.f90 | evolution (interface A member) | mem | BVK elastic slab + RST GS; rst_total_vg delivers V and the gradient; contains PBC |

### 1.3 src_VENUS (61 files)

| Routine | File | Four-class membership | Verdict | Violations / notes |
|---|---|---|---|---|
| ANGVEL | ANGVEL.f | operation (kinematics) | inv | Angular-velocity subtraction |
| rotaline | ArbitraryAxisRotation.f | mathematics | inv | Rotation about an arbitrary axis |
| BAREXC | BAREXC.f90 | operation (barrier-excitation member) | mem | NSELT=3 reaction-coordinate momentum; IRC direction via configuration |
| CENMAS | CENMAS.f90 | operation (kinematics) | inv | COM separation; QQ/PP workspaces |
| CPUSEC | CPUSEC.f | instrument | inv | Timing |
| DENQ / SUMQ | DENQ.f | post-processing (density of states) | inv | Beyer-Swinehart; DEL=1 cm⁻¹ hardcoded (could be made configurable) |
| EBOND | EBOND.f | measurement (local-mode energy) | red | Callers are only the MPLOT dead branch and LMEXCT (NACTA=4 already rejected) → dead entry; depends on MORSEB, which has no writer |
| EIGN | EIGN.f | mathematics | inv | Givens-Householder diagonalization |
| EIGOUT | EIGOUT.f | measurement (output) | inv | IP=77 and fort.77 reused (minor) |
| ENMODE | ENMODE.f | measurement (normal-mode energy) | red | The only call is inside the MPLOT==1 block; MPLOT is hardcoded 0 → dead entry |
| FGMTRX | FGMTRX.f | operation (normal-mode analysis) | inv | :26-40 Gaussian Hessian dead branch (GAUHES has no writer) |
| FINAL | FINAL.f | measurement (final-state properties) | red | :131-135 fort.777 writes U0/U1/V01 (no writers anywhere; prints undefined values); :616-620 NZDOWN branch hardcodes atom 3 (P(7:9)/W(3)). The generalized COM/EREL/OAM/ANG1-20 parts are fine |
| FINLNJ | FINLNJ.f | measurement (RKR inversion n,J) | inv | :60-63 debug-assignment residue; takes the potential via ENERGY_1 as a black box |
| FMTRX | FMTRX.f | operation (force-constant differencing) | inv | :25-60 Gaussian branch dead (GAUHES has no writer; nmo uninitialized) |
| GAMA | GAMA.f | instrument (Gamma distribution) | inv | :9 already on the RANDOM_NUMBER intrinsic stream — split from the RAND0 stream |
| grandom | GASDEV.f | instrument | red | No callers (a duplicate implementation of interface_TEST's gasdev) |
| ran1 | GASDEV.f | instrument | red | No callers |
| gasdev_0 | GASDEV.f | instrument | red | :48-71 function name does not match the return name — broken dead code |
| GLOINIT | GLO.f90 | operation (GLO unit conversion) | mem | 413.4 conversion |
| GLOSELECT | GLO.f90 | operation (GLO target-preparation member) | mem | Member on the interface B′ side |
| GLOEQU | GLO.f90 | operation | red | No callers (superseded by md_equilibrate) |
| DISTRIBUTE | GLO.f90 | post-processing | red | No callers |
| (whole file) | GLO.f | — | red | A near-duplicate of GLO.f90 not listed in the Makefile SRC, left lingering in the source tree |
| GLPAR (+PL) | GLPAR.f | mathematics | inv | Gauss-Legendre quadrature |
| GWRITE | GWRITE.f90 | measurement (trajectory serialization) | red | See the §0 calibration point: :59-60 atom1=C; :95-99 Au_nearest; :110-114 mass-threshold element names; :78-79 QZB-extremes +6.0 Å PBC box; :43 inline kB; the GWRITE_LEVEL layering is configuration (good), but each level's observable set is hardcoded — an incomplete measurement subroutine; should be purified |
| record_all | HIST.f90 | measurement (final-state recording) | red | :5-6 self-declared "C/Au only"; :19 auro+PI/3 hexagonal lattice hardcoded; :18-22 recording columns include the dead NAMD quantities AEphS |
| plot_Erel_hist and 4 others | HIST.f90 | post-processing (histograms) | red | :165,340-343 FCC(111) site classification hardcoded; the only entry is the VENUS.f90 idplt=0 dead switch → currently dead code; bin ranges/counts are parameters (good) |
| surfcell / Dist_cell | HIST.f90 | post-processing (geometry) | red | Implementation is generic (le/skew passed as arguments); the callers hardcode auro and PI/3 |
| HOMOQP | HOMOQP.f90 | operation (diatomic assembly) | inv | Orientation dispatch via NTHTA configuration |
| INITEBK | INITEBK.f90 | operation (EBK member) | mem | One of the interface B′ dual-implementation materials (Boltzmann+EBK) |
| INITQP | INITQP.f90 | operation (normal-mode assembly) | inv | :46,227-233 NACTA=4 dead branch; energy-scale resampling loop |
| JMAXCALC | JMAXCALC.f | operation (J-distribution upper bound) | red | :35-48 hardcodes H2/D2 nuclear-spin statistics by mass 1.008/2.014 — system knowledge that should be demoted to input data |
| LMEXCT | LMEXCT.f | operation (local modes) | red | F20: rejected at input; depends on MORSEB, which has no writer |
| LMODE | LMODE.f | operation (local modes) | red | Same as above; both call sites dead |
| MICROCI | MICROCI.f | operation (CI microcanonical) | red | F23: infinite loop, rejected at input |
| VOLPSCONE | VOLPSCONE.f | operation | red | Called only by MICROCI → transitively dead code |
| NMODE | NMODE.f90 | operation (analysis driver) | inv | Coordinate/gradient save-restore wrapping FMTRX |
| ORTHAN | ORTHAN.f90 | operation (orthant member) | mem | Uses the TEST predicate as a geometric bound at assembly time (initialization-time measurement; phase-legal) |
| POTEN | POTEN.f | measurement (orientation scan) | red | No callers (left over after an old task was removed) |
| POTENZ | POTENZ.f | operation (energy reference) | inv | :76-82 fort.666 writes U0/U1/V01, which have no writers |
| PRINFO | PRINFO.f | infrastructure (echo) | inv | ~70 lines of commented-out old-read residue |
| PROBJ | PROBJ.f90 | operation (J sampling) | inv | Nuclear-spin weights via WGT (filled by JMAXCALC) |
| RADAU | RADAU.f | evolution (RADAU15 closed propagator) | mem | :438-498 RADOUT embedded orchestration (calls GWRITE/TEST/FINAL/GFINAL + maintains the NAST state machine — orchestration overreach); :112,380,428 'CALLING RADOUT' debug prints. The algorithm proper is invariant |
| RADOUT | RADAU.f | measurement + orchestration | red | Should become a hook calling back to the control kernel |
| RAND0 / RAND1 / RANDST | RAND0.f / RAND1.f / RANDST.f | instrument (Knuth stream) | inv | The reproducible main stream; the F15 double call is in the caller VENUS.f90 |
| ROTATE / ROTATENO | ROTATE.f90 | operation (orientation) | inv / mem | ROTATENO is a field-induced orientation-distribution member (comment residue of NO/CO context; implementation generic) |
| ROTATEX / ROTATEY / ROTATEZ | same-named .f | operation (kinematics) | inv | |
| ROTATEX_b / ROTATEZ_b | ROTATEX.f / ROTATEZ.f | operation (kinematics) | red | No callers (2024.11 fixed-angle variants) |
| ROTATEJKM / ROTATEJM | same-named .f | operation ((J,M) orientation) | inv | ROTATEJKM has a trailing laser-polarization comment block (dead-knowledge note) |
| ROTEN | ROTEN.f90 | operation (rotational-sampling member) | mem | :101-105 debug prints + consistency-stop residue; the three laws NROT=0/1/2 |
| ROTN | ROTN.f90 | measurement (angular momentum/rotational energy) | inv | |
| SELECT | SELECT.f90 | operation (sampler framework + preparation-law dispatch) | inv | :156-160,362-400,417,1321-1326 debug-print residue; :1139-1143 MD equilibration +100 Å magic number. The NACTA/NACTB integer dispatch = closed-family registration lines; md_equilibrate calls evolution services at assembly time (allowed by 2.3 boundary 4) |
| SURF | SURF.f | operation (incident-channel member) | mem | :191 write(88) debug; :48-50 SKEW=π/4 magic number; the NSURF=1/2 dual-implementation check is in place |
| PBC | SURF.f | evolution (PBC folding) | red | The only call (ENERGY.f:17) is commented out → dead entry; PBC is actually implemented separately by rst_slab_bvk/GWRITE |
| DENS | SURF.f | operation (NREL=2 velocity distribution) | mem | Live path |
| SYMPLE | SYMPLE.f | evolution (symplectic integrator, closed propagator) | mem | :94-96 CALTYP branch residue (always true in the adiabatic version); Schlier coefficients |
| TEST | TEST.f | measurement (termination predicate) | red | See the §0 calibration point: :88-94 live C/Au criteria; :54-78/:149-192 two commented system-knowledge blocks; :106-146 generic COM criteria are the grammar-member prototype |
| THERMO | THERMO.f | operation (thermostat member) | mem | :56,115-122 debug-output residue; the two laws info=0/1 (scaling/Andersen) |
| THERMBATH | THERMBATH.f | operation (thermostat bath) | red | No callers; NTHERMB hardcoded 0, NTHMID has no writer |
| THRMAN | THRMAN.f | operation (Boltzmann quantum-number member) | mem | One of the interface B′ dual-implementation materials |
| VFDATE | VFDATE.f | instrument | red | No callers; via SYSTEM(date) shell call |
| VERLET | VERLET.f90 | evolution (VV/Beeman closed propagator) | red | :44-50,71-79 ghost Langevin hardcodes atom 2 (k=3+i); :137-148 GLO spring+random+friction embedded in the Beeman branch — interface A force contributions out of bounds into the closed layer; :153-161 same; :34-36,60-62,123-125 CALTYP residue. The propagation core proper is invariant |
| FRICFORCE | VERLET.f90 | evolution (LDFA friction+noise candidate member) | red | :130-133 call site commented out → dead entry; NFC=1/2, ready-made stochastic-law member material |
| DIAG2 | VERLET.f90 | mathematics | red | Called only by the dead entry FRICFORCE |
| LINEQ / LUDCMP / LUBKSB / FEHLER | VERLET.f90 | mathematics | red | No callers anywhere; FEHLER is an empty shell |
| WEBOND / WLBOND / WENMOD | same-named .f | post-processing (population output) | red | No callers |
| QMMICRO | QMMICRO.f | operation (quantum microcanonical member) | mem | NACTA/B=8 live path |

### 1.4 scripts (standalone tools, not in the main build)

| Routine | File | Membership | Verdict | Notes |
|---|---|---|---|---|
| conserv_check | scripts/conserv_check.f90 | post-processing | inv | Offline per-frame energy reconstruction; exposes the interface gap that the main program does not refresh T/V/H inside the evolution loop |
| smoke_rst_bvk | scripts/smoke_rst_bvk.f | post-processing | inv | Analytical/numerical gradient acceptance for the RST/BVK members — an interface A dual-implementation verification program |

**Statistics**: 99 entries — invariant class 39, interface/closed-family members 19, redundancy candidates/awaiting purification 41.

---

## 2. Architecture-level findings

### A1 Mixed occupancy of global state (venus_data)

The thing-in-itself zone (Q, Q_old, P, PDOT, QDOT, :75-84) and ~40 representation blocks (print flags, roster/geometry QZA/QZB/LA/LB, GLO parameters, thermostat parameters, GLOP ghost, the statistics arrays HIST_ATOM, etc.) live in one module with no comments distinguishing the two kinds — violating the "mixed occupancy must be reported" criterion. Additionally: `auro=2.94998` (Au nearest-neighbor distance, :585, actively used by HIST/record_all); a dozen-plus dead blocks with no writers or no consumers (MORSEB, STRETB, CUBEB, COORS, GAUSS2/GAUHES, GDD, HIST_ATOM_NAD, U0/U1/V01, the RKUTTA block).

### A2 Dual random-stream coexistence (reproducibility gap)

- Reproducible main stream: RAND0/RAND1/RANDST (Knuth 64-bit multiplicative congruential + shuffle, ISEED3 state) — SELECT/SURF/PROBJ/INITQP and other sampling run on this stream.
- Compiler intrinsic stream: `gasdev` (interface_TEST.f90:112), `GAMA` (GAMA.f:9), and `grandom` run on `RANDOM_NUMBER` — VERLET's Langevin noise, GLOSELECT, Maxwell-Boltzmann momenta, and THERMO's Andersen resampling are all on this stream.
- Consequence: ISEED cannot control the second stream, so trajectories involving noise/thermostats are irreproducible under the same ISEED; moreover **gasdev, a core instrument, physically lives in src_TEST/interface_TEST.f90 (a PES interface file)** — src_VENUS core depends in reverse on src_TEST for an instrument (layer inversion, holding up only via Makefile order sensitivity).

### A3 Interface A conformance

- **Conforms**: no compile-time PES variants (runtime TEST_PES); unit conversion centralized at the interface_TEST boundary (×23.0605, ×C1, honoring the unit contract's "conversion centralized at the interface boundary"); five members registered with a one-line case each; the EMT-NN/RST members self-contained, data externalized, path overridable.
- **Violations**:
  1. The GLO spring force and random/friction forces are embedded in `VERLET.f90:137-148`, the GLO spring potential in `ENERGY.f:40-47`, and ghost Langevin in `VERLET.f90:44-50,71-79,153-161` — per spec 3.2, "stochastic laws belong to the force-contribution family (requiring the interface A signature to acknowledge stochastic laws)", these should be interface A members, and the atom indices (atom2, NATOMA+1/+2) are hardcoded;
  2. `test_potentials.f90` members writing back representation and roster (BOXLX/BOXLY/SKEW overriding the input cell, NATOMB(1), QZB slab geometry) — members stepping outside the Q→V/F signature; geometry/roster should come via the input file or the configuration assembler;
  3. FRICFORCE (LDFA friction+noise) is ready-made stochastic-law member material, but its call site is commented out into a dead entry.

### A4 Phase-restriction check (no operations during evolution)

**Satisfied**. The evolution loop (VENUS.f90 INTEGRATION_LOOP and inside the integrators) calls only: DVDQ_1 (evolution), TEST/GWRITE/FINAL/GFINAL/ENERGY_1/record_all (measurement, including the termination predicate), the integrators. The operation class (SELECT/SURF/GLOSELECT/ROTATE* etc.) sits entirely in the assembly phase. Assembly-time `SELECT→NMODE/DVDQ/ENERGY` and `md_equilibrate→VERLET/SYMPLE/THERMO` count as "preparation members internally calling evolution services", explicitly allowed by spec 2.3 boundary 4. The only overreach is `RADOUT` (RADAU.f:438-498) duplicating control-kernel orchestration inside the propagator (calling FINAL/GFINAL, maintaining the NAST state machine, returning on its own) — an **orchestration overreach**, not a phase violation.

### A5 Compile-time coupling

- The Makefile compiles in order in a single ifx invocation (order-sensitive); gasdev's layer inversion (see A2) depends on this order.
- `src_VENUS/GLO.f` is not in SRC — an uncompiled near-duplicate file lingering in the source tree.
- The three MKL libraries (explicit Windows paths / Linux -qmkl) are linked, but no Fortran source anywhere contains any MKL/LAPACK call — a **dead build dependency** (diagonalization all goes through the in-house EIGN/DIAG2).
- The two scripts tools float outside the build system (their own compile command lines).

### A6 Incomplete control kernel (VENUS.f90)

Mixed into the control kernel: system branching (GLO layout :322-335, surface alignment :340, the THTA special case :344-346); observable selection and accumulation (BOUNCE/z_pos_min hardcoded to atom 3, AE* array recording); the outcome-handling + statistics blocks duplicated entire (:506-576 and :581-662, including the EROTA averages and reaction counting); three blocks of dead code (the NCHKP branch, the RKUTTA coefficients, the idplt histogram block). By the decomposition theorem, the control kernel should retain only the "assembly→evolution→statistics" call sequence.

### A7 Incomplete measurement layer (observable selection not demoted to configuration)

- GWRITE: atom1=C, mass-threshold element names, nearest-Au, the +6.0 Å PBC box (calibration-point details in §0).
- TEST: the live C/Au height criteria (system branching).
- record_all/HIST: the Au(111) hexagonal lattice (auro, PI/3, top/brg/hcp/fcc).
- **GFINAL.f:80**: `EREL=NORM2(P(7:9))**2/2/W(3)/C1` — overriding FINAL's computed generalized relative energy with the kinetic energy of hardcoded atom 3 (the comment admits "Erel for 1D model"); the heaviest spot of outcome-classification data polluted by system knowledge.
- FINAL.f:131-135 and POTENZ.f:76-82: fort.777/666 writing U0/U1/V01 — an NAMD leftover with no writers anywhere, printing undefined values.
- GFINAL.f:66-69: current_U0 undefined (an implicit local printed uninitialized).

### A8 Redundancy total (~26 spots)

Dead files: GLO.f. Dead entries/dead code: harmonic_sampling (no calls), VFDATE, POTEN, EBOND, ENMODE (MPLOT dead), WEBOND, WLBOND, WENMOD, MICROCI+VOLPSCONE (F23), LMODE+LMEXCT (F20), GLOEQU, DISTRIBUTE, THERMBATH, PBC, FRICFORCE+DIAG2, LINEQ+LUDCMP+LUBKSB, FEHLER, ROTATEX_b, ROTATEZ_b, grandom, ran1, gasdev_0, the RKUTTA coefficients, the NCHKP branch, the idplt histogram block, the NCVIB dead keyword, the MKL linking, the venus_data dead-block crowd. F20/F23/F14 etc. are already loud-stopped at the input layer (misuse protection), but the implementation bodies still linger in compilation.

---

## 3. Purification recommendation list (current state → target location)

1. **TEST.f:88-94 surface carbon-atom criteria** → closed predicate-grammar members: boolean combinations of `z_height(atom, threshold, outward velocity)` and `com_distance(fragment pair, threshold, recession velocity)`; system criteria assembled via the input file (:106-146 generic criteria already the prototype); the two commented system-knowledge blocks deleted.
2. **GWRITE hardcoded observables (:59-60,95-99,110-114)** → the serializer writes only the Q/P/F arrays and unit labels; element names supplied by the input roster; the C-atom/nearest-Au diagnostics demoted to post-processing configuration.
3. **GFINAL.f:80 EREL override (hardcoded atom 3)** → delete the override and adopt FINAL's generalized COM relative energy; where a subset of degrees of freedom is truly needed, let the outcome-classification configuration specify it.
4. **GLO/ghost forces in VERLET.f90:44-50,137-161 and ENERGY.f:40-47** → interface A force-contribution family members (a GLO-spring member, a Langevin-ghost member, the signature acknowledging stochastic laws and friction); the integrator/energy routines return to pure propagation/pure functionals. FRICFORCE can be revived along the way as the LDFA member.
5. **venus_input.f90:616-621 ghost parameter assembly (W(2) hardcoded)** → moved into the GLO/Langevin members' own init hooks (POTPRE-style).
6. **venus_data.f90:585 auro and the HIST hexagonal lattice (PI/3, top/brg/hcp/fcc)** → representation input parameters / surface-lattice-type configuration (slab roster data).
7. **venus_params.f90:88-98 NN potential-type indices** → constants internal to the NN PES member file.
8. **JMAXCALC.f:35-48 H2/D2 nuclear-spin statistics (mass discrimination)** → input data (a nuclear-spin weight table).
9. **venus_input.f90:997-998 BOXLX/BOXLY default geometry and :703-704 NATOMB hardcoded domain** → mandatory inputs or derived from the slab roster; validation generalized into a roster-consistency check.
10. **VENUS.f90 observable accumulation and the duplicated outcome blocks** → a measurement-recorder module + a single outcome path; the atom indices of BOUNCE/z_pos_min demoted to configuration; the dead blocks (NCHKP/RKUTTA/idplt) deleted.
11. **gasdev's location and random stream** → moved into a standalone instrument module in src_VENUS; a single random stream unified codebase-wide (either all on RAND0 or all on the intrinsic RNG), restoring full-stream ISEED reproducibility; GAMA's IDUM dummy argument made effective again.
12. **test_potentials.f90:152-163,260-293 members writing BOXLX/NATOMB/QZB** → a configuration assembler (roster/geometry delivered via the input file), the PES members keeping only the Q→V/F mapping and data-file reading.
13. **RADAU.f RADOUT embedded orchestration** → the propagator calling back to the control kernel (output hook + termination-predicate evaluation point), orchestration returned to the control kernel; the three 'CALLING RADOUT' debug prints deleted.
14. **~26 spots of dead code (the A8 list)** → delete; if harmonic_sampling is kept, register it formally as an interface B′ member (the surface-oscillator Boltzmann law), forming a third implementation with THRMAN/INITEBK.
15. **venus_data partitioning** → at minimum, comment-delimit the thing-in-itself zone/representation zone/instrument zone, with the dead blocks cleaned as deleted; the MKL linking removed (or its LAPACK used to replace EIGN — either one, eliminating the "linked but unconsumed" state).

---

## 4. Concluding remarks

This codebase has already reached the closed form the spec requires in **interface A's registration shape** (runtime TEST_PES, unit conversion centralized at the boundary, self-contained members + externalized data) and in the **input grammar** (keyword=value + constructor mapping + presets made data); interface B′'s dual-implementation materials (THRMAN=Boltzmann, INITEBK=EBK) and the target-preparation members (GLOSELECT, md_equilibrate, the SURF incident channel) are all in place. The main deficits concentrate in three places: (1) **stochastic laws not in interface A** (the GLO/ghost/LDFA forces hardcoded inside the integrator and energy routines); (2) **the measurement layer's observable selection not demoted to configuration** (four spots of system knowledge in GWRITE/TEST/HIST/GFINAL, including GFINAL.f:80's hardcoded override of the generalized EREL); (3) **mixed occupancy and dead code in the control kernel and global state** (~26 spots of redundancy, dual random streams, the duplicated outcome blocks). After purification per the §3 list, the variable code will concentrate in interface A/B′ members and closed-family registration lines, aligning with spec 2.3's variable-code distribution conclusion.
