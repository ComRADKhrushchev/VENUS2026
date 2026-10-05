# VENUS Streamlining and Generalization: Program Completeness Design Specification

Status: design-rationale document (a derivation, not an implementation document)
Date: 2026-09-03 (first revision round of the month)
Scope conditions: this document discusses only the program's abstract architecture — the categories of physical quantities, the complete classification of subroutines, and the theoretical form of the extension interfaces. Variables related to program execution (random-stream state, restart/checkpointing, process management) and unit-system conventions belong to the empirical implementation layer and are excluded from this derivation; the unit contract is attached as a one-line ancillary obligation at the end of Section 3.

---

## 1. Definition and Derivation Goal of the QCT Program

**Definition**: a QCT program is a program that, within the Born-Oppenheimer framework, samples initial states by Monte-Carlo, propagates trajectories according to the Newtonian equations of motion (supplemented, when electronic degrees of freedom are present, by electronic evolution equations), and performs statistics on the final states. This definition fixes the program's skeleton: ensemble sampling → trajectory propagation → final-state statistics. The skeleton does not change with any specific system.

Under this definition, the QCT program is a **family**: specific computational systems (gas-phase bimolecular, surface scattering, scattering with a dissipative bath, ...) are members of this family. The program's task is to parameterize the family members.

**Derivation constraints**: only a pure category analysis of physical quantities is used. Physical quantities divide into thing-in-itself and representation.

**Definition of the thing-in-itself (irreducibility criterion)**: a thing-in-itself is a variable that at no moment can be further reduced, i.e. expressed as a compositional equivalent of other variables. By this criterion, quantities such as the kinetic energy are not thing-in-itself — it can be written as a function of the momentum and is a functional of the thing-in-itself. The complete list of thing-in-itself is:

- Q: nuclear coordinates;
- P: nuclear momenta;
- s: electronic state (under the Born-Oppenheimer approximation this can be further specified as an occupation function; in the single adiabatic surface case it is the trivial point thereof);
- F: force field — reducing the computation of the force field from coordinates requires a system specification, so the force field is not reducible to a functional of other thing-in-itself variables.

The trajectory time t is held by the control kernel as an orchestration variable and is not included in the classification of physical quantities.

**Definition of representation**: all stipulations and characterization data for the system — observables are a subset of representation; the quantum numbers, temperature, collision energy, and similar characterization quantities of the initial-state stipulation likewise belong to representation; reducible quantities such as the kinetic energy also belong to representation (as functionals of the thing-in-itself). This definition of representation is derived from the completeness requirement of Section 2; see the repair note in that section.

The derivation requires two system-independence assertions — the form of the thing-in-itself is independent of the system specification, and the general form of scattering is independent of the system specification — which are established as lemmas in Section 2.3 rather than used as unproven premises.

**Derivation goal**: a completeness norm — maximize, relative to the designated extension interfaces, the amount of code that is encapsulable (unchanged as the system is extended). Concretely it decomposes into two propositions:

1. Completeness: the four classes of pure subroutine are formally sufficient to implement a full QCT program;
2. Optimality: the all-pure-subroutine architecture attains the maximum amount of encapsulable code.

Goal 1 is established by the decomposition theorem of 2.2 and strengthened in 2.3 into a conclusion about the class distribution of variable code; goal 2 is established by the exchange argument of 2.4.

---

## 2. Completeness of the Four Classes of Subroutine

### 2.1 Type system

A subroutine is defined as a single type-conversion unit: it takes some variables as input and produces some variables as output. Variables are dichotomized into thing-in-itself/representation, and subroutines are classified into four classes by the combination of input and output types:

| Class | Input → Output | Physical content |
|------|------------|---------|
| Evolution | thing-in-itself → thing-in-itself | all computation that changes the physical state |
| Measurement | thing-in-itself → representation | translating the physical state into characterization data |
| Operation | representation → thing-in-itself | assembling a stipulation into a physical state |
| Post-processing | representation → representation | statistics and secondary processing, never touching the physical state |

**Repair note**: if representation were taken strictly in the sense of "observables", the operation class would be empty (the inputs of initial-state assembly — quantum numbers, temperature, collision energy — are neither thing-in-itself nor observables). The completeness requirement forces the extension of representation to be corrected to "all stipulations and characterization data". This does not introduce a third class of variable; it is a self-consistentiation of the second class's definition.

**Enumerative completeness**: the 2×2 matrix over the dichotomized variables exhausts all types of pure subroutine; a fifth class does not exist.

### 2.2 Division of the Program's Internal Functions

**Proposition (decomposition theorem)**: any QCT program can be decomposed into three functional layers:

- (a) the **control kernel** — the call sequence of assembly, looping (evolution step → termination test → periodic recording), and statistics, containing no physical knowledge;
- (b) the **four classes of pure subroutine** — the entire computational content of the program;
- (c) **interface member implementations** — the members among the four classes of subroutine that are replaced as the system specification changes.

(c) is not a fourth kind of code: interface members themselves still belong to one of the four classes (e.g. the potential energy surface is an evolution-class service; the preparation law, after the relocation of 2.3, is representation-layer code). The distinction between (b) and (c) is not in the kind of code but in its mutability.

**Proof (structural)**: by the enumerative completeness of 2.1, any subroutine in the program must be one of the four classes. Physical knowledge can reside only in computation that changes the physical state (necessarily evolution class) and computation that translates the physical state (necessarily measurement class); the physical content of assembly (quantization rules, etc.) must be in the operation class; statistics does not touch the physical state and must be in the post-processing class. What remains after stripping these away is only the call sequence itself, i.e. the control kernel, whose physical ignorance is guaranteed by construction. ∎

**Establishment of completeness**: this decomposition simultaneously establishes goal 1 of Section 1 — the control kernel contains no physics and all computation falls into the four classes, so the four classes of pure subroutine (plus control-kernel orchestration) are formally sufficient to implement a full QCT program.

**Delimitation of the main loop**: the control kernel is not a subroutine, just as the main loop of a differential-equation solver is not the right-hand-side function. It is the orchestration layer — minimal, constant, insensitive to extension; its existence does not compromise the maximization of encapsulability.

### 2.3 Encapsulation Theorem for the Non-Evolution Classes

**Forward-reference note**: this section uses Interface A (the force-contribution family — the part of the evolution class that varies with the system specification) and Interface B′ (the initial-state distribution family — the part of the representation layer that varies with the system specification); their formal definitions and irreducibility arguments are in Section 3, and they are used here in the intuitive sense just stated.

Write the system specification R = (content specification D, formal specification F).

**Proposition**: there exist fixed implementations of the measurement and operation classes of pure subroutine such that, for any system specification R in the QCT family, the program's functional integrity is preserved (computing the same cross sections and statistics as the unclosed program), and the code of these two classes does not change with R. The code that changes with R (added or modified) is concentrated in the Interface A/B members and in new members of the closed selection families (including their registration lines).

The post-processing class is no longer proved as part of the base program: post-processing is secondary processing of measurement records; the triviality of its code follows directly from the self-closure of the representation layer (representation → representation) and, per this specification, is not separately argued.

The derivation proceeds from two system-independence lemmas.

**Lemma 1 (form-independence of the thing-in-itself)**: the form of the thing-in-itself variables Q, P, s, F — their structure as abstract types — is independent of R.

Proof: the dichotomy of form and content. The type structures of Q, P, s (coordinate array, momentum array, electronic-state object) are fixed; the degree-of-freedom roster merely determines their dimension. The dimension is a content parameter the type carries at runtime, not a change in how the type is constructed — just as an integer does not change its type by taking a different value. The type structure of F (a mapping from coordinates to forces) is likewise fixed; its concrete functional content is stipulated by R, but the form of the mapping is untouched by R. Hence the form of the four thing-in-itself items is the same for all R. ∎

**Lemma 2 (independence of the form of scattering)**: the general form of the scattering process — with the cross section at the core of the observable concept — is independent of R.

Proof: cross section = the quotient of event rate by incident flux; its definition uses only three ingredients — the initial-state distribution, the classification of trajectory outcomes, and ensemble counting — all three given by the QCT definition. R determines only the concrete membership of the event classification (what counts as scattering, what as adsorption) and does not change the form of the cross section itself. ∎

**Theorem (form invariance of the measurement class)**:

1. The input of the measurement class must be a subset of the sequence {Q,P,s,F}; by Lemma 1 the form of the sequence is fixed, independent of R; hence the purely formal operations in the serialization code contain no knowledge of R;
2. A measurement can be written as M({Q,P,s,F}), where M may contain only purely formal operations, e.g. f(Q(1),P(1),s(1),F(1)); for system-dependent measured quantities, these can always be written as operations in a post-processing step; therefore the measurement class in the program can possess form invariance, and hence closability.

**Theorem (closability of the operation class)**:

1. By definition, the input of the operation class is a quantity independent of the current system's physical information (otherwise it would be evolution); hence it must be a stipulation of the initial state (otherwise it would break self-consistent physical evolution and be physically meaningless);
2. It follows that any initial-state stipulation can be delivered to the program only in the form of a distribution, so the operation class contracts to a single fixed construction: sampling (Q₀, P₀, s₀) from a given distribution, where the distribution is stipulated by physical law and is therefore necessarily system-independent;

**Conclusion (summary)**: the entire non-evolution part is independent of R — the form invariance of the measurement code is given by Lemma 1 via the closure of serialization and predicate interpretation; the fixity of the operation class is given by the intermediary of the initial-state distribution; post-processing is trivial. All variable code lies in the evolution class (Interface A) and the representation layer (Interface B′). This conclusion is complementary to the exchange argument of 2.4: the exchange argument gives global optimality, while this conclusion gives the class distribution of the variable code.

**Boundaries (honest declaration)**:

1. Mutability is relocated, not eliminated — the irreducible lower bound, namely the interface member implementations, remains; what this conclusion rewrites is their address (evolution-class members + representation-layer members);
2. The defense in Lemma 1 that "the dimension is a content parameter" presupposes that type structures can carry a runtime dimension (satisfied by Fortran allocatable arrays);
3. The coverage of the predicate grammar and of distribution representations is an empirical premise that must be checked by enumeration against the known families of criteria and preparation laws;
4. Phase correction: the actual invocation during the assembly phase = operation-class framework + Interface B′ distribution member (representation-layer code); a preparation member may internally call evolution services — this is internal member structure and does not violate the directionality restriction that "no operation-class calls occur during the evolution phase".

### 2.4 Optimality (Exchange Argument)

The argument proceeds in two steps.

Step one (feasibility): by the decomposition theorem of 2.2, the all-pure architecture can implement a full QCT program — without feasibility, optimality cannot be discussed.

Step two (optimality): there is an irreducible lower bound on the amount of variable code — the implementations of the interface members (every new formal specification necessarily brings new implementation code). The all-pure architecture exactly attains this lower bound: the code modified upon extension is exactly the interface implementations. Any impure subroutine deposits physical knowledge outside the interfaces and, when the corresponding specification changes, is forced to modify the encapsulated region, giving a variable-code amount strictly greater than the lower bound.

**Pure-completeness corollary**: the all-pure-subroutine architecture is the architecture with the maximum encapsulable code — completeness (2.2) guarantees its feasibility, and the exchange argument guarantees its optimality.

### 2.5 Phase Restriction (Derivational Status of the Three-Phase Picture)

The three phases (assembly → evolution → statistics) are not an independent assumption but a corollary of the four-class classification plus the phase restriction:

| Phase | Classes permitted to be called |
|------|--------------|
| Assembly | operation class (plus initializing measurements) |
| Evolution | evolution class and measurement class (termination test and recording; the termination Boolean is a degenerate member of representation) |
| Statistics | post-processing class |

No operation-class call occurs during the evolution phase — runtime intervention would destroy the purity of the trajectory flow; this is the precise expression, within the QCT concept, of "a trajectory, once launched, is autonomous".

### 2.6 Invariant Criteria (Acceptance Gauge)

The completeness standard is the norm; the criteria are the gauge. For existing code, each routine is examined with three questions:

1. Does it handle only the physical quantities fixed by the QCT concept (Q, P, s, F as an abstract state)?
2. Does it request services only through interface contracts, treating the service provider as a black box?
3. Does its interior contain no branches or constants hard-wired for a specific system / potential energy surface / method?

A routine passing the three questions is of invariant class; one carrying component knowledge is an incomplete subroutine (e.g. of the GWRITE type: a measurement-class unit with embedded system knowledge), whose fate is not deletion but purification — the system knowledge is parameterized and enters through an interface. The reliability of the criteria presupposes the completeness of the interface contracts; contract changes are rare, versioned events.

---

## 3. Definition and Completeness of the Interfaces

### 3.1 Definition of Interfaces and the Freezing Operation

**Definition of an interface**: code that does not possess invariance under arbitrary changes of the computational system. (The three layers of meaning of invariance: the signature is invariant under all changes; the closed member family is invariant under replacement within the family and variable under family extension; an open interface is not even closed as a family.)

Terminology convention: "interface" is taken broadly — both closed selection families and open interfaces fall under it; the "irreducible interfaces" of Section 3.3 refer to the open part thereof.

**Definition of an input file**: a stipulation by which a system can be defined without changing code. The **quantities** (numerical parameters) and **qualities** (choices among prestored schemes) in an input file correspond to the parameters and constructors of a fixed reading grammar: system = C(q₁,…,qₙ).

**The freezing operation**: for those parts of the interface with a bounded functional scope (e.g. the nuclear propagator), several choices are written in advance and closed. Freezing does not eliminate mutability; it addresses it: switching a preset changes only the input file (zero code touched); adding a preset only adds a member body plus one registration line (zero touching of existing code). The essence of freezing is to move the system specification from code into data.

**Completeness of the interfaces** (a subset of the completeness corollary): interface design takes the minimization of code changes via input files as its completeness standard — the closed layer covers all enumerable functional families, and the actually open interfaces (those requiring new code to extend) shrink to the minimal set, namely the irreducible interfaces of Section 3.3.

### 3.2 Decomposition of the Residual Category

A category analysis of quasi-classical systems distributes the mutability over three levels (absolute completeness):

1. **Nuclear level**: the nuclear propagator. The mutability is trivial (the Verlet family is small and stable); closed.

2. **Electronic level**: the mode of nuclear-electron coupling. By the factorization theorem it decomposes into two independent axes:

   - **H-field provider** (generalization of the potential energy surface): delivers Q ↦ (H(Q), ∇H(Q) or equivalent), independent of the electronic method. The delivered representation is stratified by **capability level**: single surface / multiple adiabatic surfaces / full diabatic matrix / independent-electron tensor.
   - **Electron propagator** (mixed quantum-classical methods): the consumer of that representation (adiabatic / Ehrenfest / FSSH ...). Adiabatic is a degenerate point (when consuming the single-surface representation, the composition collapses to a pure PES); the interface must be well defined at degenerate points.

   The two axes are independent (any PES × any method is a legal pairing — a Cartesian product structure), composed through a fixed intermediary functional (adiabatic: F=−∇V_α; Ehrenfest: F=−⟨ψ|∇H|ψ⟩; FSSH: current surface −∇V_α). Bundled implementations (writing the knowledge of H and the electronic method into the same routine) are architecturally excluded by the factorization theorem; the number of implementations drops from |potentials|×|methods| to |potentials|+|methods|. The legality of capability levels and methods is type-checked when the input is read; illegal combinations are explicitly rejected with the two conflicting axes named.

3. **Scattering level**: the a priori stipulation of scattering — necessarily an interaction event belonging to several reactants; stipulating the two stipulates the scattering. Accordingly it is decomposed into:

   - **Stipulation of the reactants = pure configuration stipulation**. This is not a background remark but the residue after exhaustion: the dynamical content (forces) has already been assigned to the PES by factorization, the statistical content to the initial-state distribution; what remains on the reactants is only identity and geometry. The configuration stipulation further dichotomizes into the molecular formal specification (roster, geometry) and the initial-state statistical distribution, both determined by the input file. The degree-of-freedom roster (rigid surface, movable slab, GLO ghosts) is housed here. Empirically there are two stipulation paradigms — molecular scattering and surface scattering — which constitute the dual-implementation test of this level; roster and geometry are pure data, not code degrees of freedom;
   - **Stipulation of the event**, fully decomposed into three blocks: the incident channel (the distribution of relative motion: collision energy, impact parameter, orientation — assigned to the initial-state statistical distribution), the termination predicate (a closed grammar, see below), and the outcome classification (what counts as scattering/adsorption/dissociation — by the measurement form-invariance theorem, this is post-processing configuration, data not code).

   Consequence of the decomposition: the former "scattering scenario descriptor" is no longer a whole standing above the other slots — roster and geometry recede to pure data, the noise law goes to the force-contribution family (requiring the signature of Interface A to admit stochastic-law members, e.g. the GLO Langevin random force), the termination predicate goes through the closed grammar, and the outcome classification goes to post-processing configuration; the sole surviving open code degree of freedom is the initial-state distribution family (Interface B′, see 3.3). The electronic initial state s₀ (the initial stipulation of the occupation function) is explicitly assigned to the initial-state statistical distribution; its representation grows with the capability level and is a matter of grammar versioning.

### 3.3 The Irreducible Interfaces

After the above decomposition and freezing, the actually open interfaces are irreducibly two. Both are registry families (new member body + one registration line, zero touching of existing code), but they differ in type and cannot be merged:

**Interface A: the force-contribution family (generalization of the PES)**
The total force on any system is a sum of contributions (surface interaction, GLO spring term, constraint projection, external-field driving, random force, ...). The open interface is not "one potential slot" but an additively composable family of contributions: each member = code + external data (NN weights, fitted parameters). The type signature is Q ↦ F (evolution class). Root cause of non-closure: the function space is unbounded and members carry data; the signature must admit stochastic-law members.

**Interface B′: the initial-state distribution family**
The member family taking stipulations → initial-state distributions (representation layer): internal-state quantization laws (Boltzmann, EBK, Wigner, adiabatic switching, ...), target preparation (lattice thermalization, oscillator equilibrium), and the stipulation of the electronic initial state s₀. The type signature is representation → distribution. The root cause of non-closure is isomorphic to A's and irreducible: the quantization rules are a research-open family; there is no combinatory closure of the "grammar of distributions" kind — criterion comparison: the closure of the termination predicate rests on the compositionality of Boolean algebra, that of the nuclear propagator on the boundedness of its functional scope; the distribution family has neither. A and B′ cannot be merged: their type signatures differ, and merging would destroy the decomposition structure of 2.2 (evolution-class services vs representation-layer stipulations).

**Closed layer (for contrast)**: the nuclear propagator and the electron propagator are selection families of bounded functional scope; the termination predicate is a Boolean-combination grammar; the outcome classification is post-processing data. None belong to the open interfaces.

**Irreducibility note**: the force-contribution family and the initial-state distribution family are each not further reducible — the former by the unboundedness of the function space, the latter by the research-openness of the quantization rules. The open/closed distinction can here only be stated empirically: registry families are research-open; closed families are functionally bounded.

**Unit contract (ancillary obligation)**: the semantics of an interface signature are defined for the physical quantities themselves, independent of the units representing them; but at implementation time a single unit convention must be pinned down (inside this program Å/amu/10 fs, see CLAUDE.md), with conversions concentrated at interface boundaries. The unit system does not shape the architecture and therefore does not enter the main derivation.

---

## 4. Formal Definitions of the Irreducible Interfaces (Partially Filled In)

### 4.1 Internal Structure of Class-B Members (Three-Layer Decomposition of MQC Propagation)

**Proposition**: any MQC dynamical propagation can be split into three successive sub-phases — nuclear propagation, electronic propagation, and state-statistics matching. This is a miniature replica, inside the evolution step, of the control kernel's three phases (assembly → evolution → statistics).

- **Nuclear propagation**: advance (Q, P) given the forces. A closed-family member (Verlet/Symple/Radau), independent of all MQC methods — this layer was separated long ago in the decomposition, and the three-way division confirms its status;
- **Electronic propagation**: advance the electronic variables given Q(t) — the amplitude ODE of wavefunction methods (TDSE integration along the trajectory), the mapping-phase-space flow of mapping methods, or analytic elimination (friction/noise kernels);
- **State-statistics matching**: read the electronic statistics (populations, gaps, derivatives) and write the classical side (discrete reassignment of the active-surface index) and the nuclear side (momentum-projection adjustment) — the enforcement mechanism of quantum-classical statistical consistency. The sole consumption point of the random stream during the evolution phase lies in this layer.

**Generalized definition of state-statistics matching**: every method enforces quantum-classical statistical consistency by some mechanism of its own — FSSH by hopping rules, MDEF by the fluctuation-dissipation relation (borne by a noise-law member), SQC by zero-point constraints. Under this definition the third layer is complete for the MQC family; Ehrenfest is a degenerate point that actively abandons consistency (at the price of mean-field illusion).

**Method × three-layer degenerate-point table**:

| Method | Nuclear propagation | Electronic propagation | State-statistics matching |
|---|---|---|---|
| FSSH | closed family | amplitude ODE | fewest-switches hopping + momentum adjustment |
| Ehrenfest | closed family | the same amplitude ODE | identity (no active surface) |
| IESH | closed family | independent-electron amplitude ODE | per-electron hopping |
| LZSH | closed family | degenerate (amplitudes not integrated) | LZ probability at crossing points |
| SQC/CSDM | closed family | mapping-phase-space flow | final-state window binning |
| MDEF | closed family | analytic elimination (TDPT) | fluctuation-dissipation (noise-law member) |

Corollary 1 (decomposition of the method concept): Ehrenfest and FSSH share all electronic-propagation code, differing only in the force rule and the matching rule — a "method" decomposes into fewer reusable parts, and the method count drops multiplicatively to a combination of "force rule × matching rule".

**Three restrictions**:

1. The three layers are not a free Cartesian product but a stratification with dependency declarations: the matching rule declares its requirements on electronic propagation (FSSH needs the amplitudes and their time derivatives; LZSH only the gap and slopes in the crossing region; mapping methods only the final state). This is isomorphic to the capability-level mechanism: each layer declares its capability requirement on the H-field provider (adiabatic single surface → add gradients → add nonadiabatic coupling vectors — the momentum projection of the matching layer raises the requirement to the highest level), and the facade takes the union at assembly time for a legality check;
2. Decoherence corrections straddle the layer boundary between electronic propagation and matching: DoM-type schemes modify the amplitude equation (a declarative sub-option of electronic propagation); decay-of-mixing-type schemes reset the amplitudes after a failed hop (post-processing of matching) — decoherence models enter as sub-options and do not break the decomposition;
3. Inter-trajectory coupling methods (CT-MQC, AIMS) break the premise that "matching reads the statistics of a single trajectory by itself" and remain outside the QCT family definition (consistent with the scope-boundary clause of 3.3).

**Architectural placement**: the electronic method members (electron propagation laws) are further decomposed internally into two registration units — the amplitude propagator (cross-method shared infrastructure; RK/matrix-exponential integrators for a given representation, bounded functional scope, closable) and the matching rule (a method-specific small unit, the consumption point of the random stream). Refined diagnosis of the historical pathology: the over-reacher that forced VERLET in TDHF-FSSH is the matching layer (the momentum projection requires force evaluation in step with integration); it should be handled as a dependency declaration from the matching layer to the nuclear-propagation layer, negotiated at assembly time, not as a hard constraint.

**Benchmarking feedback**: the three Tully 1990 models (R8-5) degrade nuclear propagation to the trivial (one-dimensional translation) and isolate electronic propagation and matching for testing — model 1 tests the asymptotic correctness of matching (the LZ limit); model 2 tests the time-step coupling of the two layers. The acceptance targets of R8-5 are stated per layer.

### 4.2 Formal Signatures (Backfilled 2026-09-16 at S/T7; E27 container contract)

[E27 revision 2026-09-16 — the ws interchange and the h_form/h_deliver matchmaking are RETIRED
(superseding the relevant parts of E10/E13/E14): the physical-system container is the composition
point of system-specific code. It delivers
the composed force evaluation directly; system-specific method packages may live in their
systems/ folder and share private modules with the H definition — no lingua franca. The entries
below are the landed forms (src/facade/), verified by the S-phase audit
(docs/completeness_audit_src.md) and the V/S case set.]

[E28 revision 2026-09-22 — the container parameter channel rides the program input: the
run-directory *_pars.txt channel is RETIRED. A container's parameter vocabulary is a set of
declared container keys of input_qct.txt (the box owns the key list and a params_t whose field
defaults are the container defaults; the reg file declares the vocabulary before the parse via
registry_gen_keys and consumes the staged rows at its seam by pull(key, target [, found]);
input_audit_declared fatals a staged-but-unconsumed declared key). E21/E22's closed-grammar
posture is unchanged — admission stays parse-time, the membership set gains the declared names;
a declared key may not collide with a program or member key. Ruling record: docs/plans/T.md.]

**Interface A (the force boundary, force_facade.f90)**:

- `container_pes_i(q, mass, g, v)` — q [Å, all dof flattened], mass [amu] in;
  `g = dV/dq` of the composed ACTIVE configuration [eV/Å] and `v` [eV] out, BOTH fully assigned
  by the container (the active-surface fold, the member sum, and stochastic laws are
  container-side; rng single-stream inside; an H-shaped container folds its own eigenproblem).
- `force_eval(sta, e_tot)` — the one evolution-phase force entry: container call → the unique
  unit conversion `f = -g*ev_kcal*e_conv` (23.0605 × 0.04184 lives only here, T4) → E = T(P)+V.
- Assembly: `container_bind_pes(proc)` (single slot, one container per run),
  `container_bind_term(proc)` (optional termination judgment, read-only sta),
  `container_init()` (method name lookup + the no-container-bound posture; purely formal, E13).

**Interface B (the electronic axis, elec_facade.f90)**:

- `elec_prop_i(sta, dt)` / `mqc_match_i(sta)` / `sync_proc_i(s)` — method-package members;
  data reach is the member's own channels (its folder's private modules / level-2 components);
  the match member is the sole evolution-phase rng point and writes back through
  `occ_set(occ)` / `n_hop_set(n)` (run-lifetime holders; getters for the recorder).
- Registration: `method_reg(name, prop, match [, sync])` from methods/ (generic) or a
  systems/ folder (container-composed) alike; `method_bind(name)` a pure name lookup
  (ELEC_METHOD key; the degenerate 'adiabatic' row is the default — an adiabatic reference
  run of ANY system stays expressible).

**Interface B′ (initial-state sampling, sampler.f90)**: member pairs registered by
`samp_reg(name, draw, realize)` (the 2026-09-22 two-phase split: draw consumes the RNG into a
module-private record, realize writes q/p and draws no RNG), selected per fragment by
DIST_SCHEME codes; `samp_run(scheme, sta)` dispatches fragment members in roster order, then
the ONE paradigm-resolved incident member once - the incident channel is three members
(`inc_surface`/`inc_pair`/`inc_single`, 2026-09-28 split; `samp_incident_word()` resolves,
the regs register under paradigm guards, and the trivial `inc_single` is a bitwise no-op that
replaced the historical skip-guard).

**Well-definedness at degenerate points (verified)**: (i) the adiabatic degenerate package —
all-no-op members, occ seeded [1], zero rng consumption (check_facade (h)/(j)); (ii) the
single-fragment surface paradigm — incident skipped, degenerate histogram channels pinned
(case_surface_window A4, dual to the gas-phase thta degeneracy); (iii) no container bound —
assembling without a container dragged in succeeds and aborts with a named error only at
force_eval's first call (check_facade posture notes).

**Dual-implementation verification records (this repository's carriers; the master plan's
RST/EMT-NN/GLO-stochastic names are full-form assets outside this repo, kept as historical
targets)**: Interface A — three analytic potentials (harmonic/morse/leps: gradient
self-consistency vs central difference, plus the OLD-CHAIN BIT-LEVEL comparison
max|dF| = max|dE| = 0 per potential and geometry, re-verified green AFTER E27); Interface
B′ — eleven dist_* members covered by the five-case set and eleven check_dist_* criteria
programs (mechanism level) + on-chain numeric posteriors (b_layer χ² exact bands, chain_freq
spectra bitwise to 1e-13); a twelfth member dist_wigner (ground-state Wigner normal-mode
distribution, 2026-09-21) joined after the S-phase record — criterion-covered
(check_dist_wigner, sixteen legs green), five-case coverage pending; configuration
level — the two paradigms (gas-phase four cases,
surface POSCAR-cell case); electronic axis — the adiabatic degenerate point (the nonadiabatic
axis lands with the first container-composed package, per E27).

---

## Appendix: Relation to Existing Tasks

- Deletion verdicts: any functionality that belongs neither to the control kernel and the four classes of pure subroutine nor to expressibility as an Interface A/B member or a closed selection-family member is redundant (the RK4 dead entry, ADAMS, and compile-time PES variants are thereby eliminated);
- Input keyword refactoring: the existing keyword=value input moves toward a "constructor + parameters" grammar; quality parameters select closed-family members, quantity parameters fill in numerical values;
- Task tracking: Notion "VENUS-TDHF generalization" (phase P1 architecture).
