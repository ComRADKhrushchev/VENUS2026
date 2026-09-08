# methods/ - the flat built-in method library

The library home of the generic methods (E14 spec §5, rebuilt from the former
interface_electronic/ + interface_distributions/): pre-stored, parameter-free method
packages and reusable algorithm components that any system binds by name. Everything
here is physics-generic code — the same FSSH machinery serves any system; system
specificity enters only through the ws container, state, and config at runtime (E12
ruling: system material lives in systems/<system>/; a method that embeds system
parameters or system knowledge does not belong here — see also the
knowledge-ownership audit: system/method-specific knowledge must not enter generic
code).

## File prefixes (self-describing families)

| Prefix | Family | Registered by | Present / planned |
|---|---|---|---|
| `elec_*` | MQC method packages — {prop member, match member [, sync member]}, one `method_reg` row each | its thin `reg_<package>.f90` (assembly phase) | `elec_adiabatic` (present — the all-no-op degenerate seed); `elec_tsh` / `elec_iesh` / `elec_ehrenfest` phase C |
| `dist_*` | initial-state distribution members (the B' sampler family; parameters via config — E12 wording carried, zero embedded system data) | thin `reg_dist_*.f90` at C5 together with the sample routines (K2 E3 slimming ruling) | boltzmann / ebk / incident (present — migrated from the former interface_distributions/, E14); the remaining 9 arrive C5 |
| `amp_*` | shared amplitude integrators — library modules, not packages (Corollary 1: FSSH and Ehrenfest share the same prop code; the method count drops multiplicatively) | none (components, consumed via use — level 2) | `amp_rk` / `amp_expm` phase C |
| `reg_*` | thin registration modules — build material, not methods | (they are the registration) | `reg_adiabatic` (present — the sole reg file of the K2 seed) |

## Two exposure levels of the library ("ready-made semantics", spec §5)

**Level 1 — whole packages** (the standard path, zero code):
`config%electronic%method = 'fssh'` binds the assembled package at blackbox_init; the
mqc_statistic / elec_prop facade dispatch lands directly in the package's match/prop
members. The blackbox author writes no method code at all. Dispatch chain of
mqc_statistic: facade table lookup -> match member (reads sta%s%a + sta%p + ws; g_ij
from the library formula; rng-sampled hop; writes sta%occ/sta%p; optional force_eval
self-call) -> sync tail.

**Level 2 — reusable components** (the assembly path): the library modules expose
algorithm-level public routines that a custom package's members may use to compose
variants — `amp_rk` / `amp_expm` (amplitude ODE integrators, propagation-layer
components), `fssh_g_ij(a, d_ij, v, dt)` (the fewest-switches probability matrix,
generic formula), `hop_decide(g_ij)` (the rng sampling decision — the sole
evolution-phase rng point), `momentum_rescale(p, d_ij, dE)` (the projection solve
p' = p - beta*d_ij), `dom_reset(a, ...)` (decay-of-mixing amplitude reset). A
system-tuned variant (e.g. a custom decoherence schedule) = a thin match member that
uses these components — the physics formulas stay single-copy; the author writes only
the scheduling logic.

## Drag-in rule (a new method = zero edits to existing files)

A new method is one new file in methods/ + one thin `reg_<package>.f90` calling
`method_reg` — drag in and it compiles (the Makefile wildcard picks it up; the K2
skeleton's explicit compile list is the transitional form until the K3 build). The
prebuild step aggregates every reg_*.f90 under methods/ and systems/ into
registry_gen.f90 (E11 drag-in generalized by E14; effective at K3), so the driver's
final form carries zero registration lines. Until then the K2 skeleton wires
`reg_adiabatic_all()` from the driver pseudocode transitionally.

## Registration and pairing legality (E13 carried)

The facade is a purely formal authority: `method_reg` fills one row and checks only
name-table hygiene (duplicate name / full table). That one registration vouches for
the whole package — which prop member pairs with which match member is the
registering author's own knowledge, not the facade's. A member loud-stops at runtime
on a ws%repr it does not support (E10 contract) and on dt semantics it cannot live
with (E13 open point — fixed-step-only laws, shape landed in phase C with elec_tsh).

## The dist_* family (initial-state distributions — the B' sampler registry)

`dist_*` members are pre-stored generic initial-state distribution methods: each
draws one fragment of the trajectory's initial state per its prescription
(Boltzmann normal-mode quantum numbers, fixed (n,J) EBK levels, incident-channel
collision energy / impact parameter / orientation). E12 wording carried: the
members are physics-generic code with zero embedded system data — temperatures,
frequencies, quantum numbers, and mode switches all arrive as member params
filled from config at runtime; system material lives in systems/<system>/.

Registration timing (K2 E3 slimming ruling): `samp_reg` registers a
`dist_proc_i(sta)` sample routine, and the K2 skeletons carry only `*_init` —
the sample routines land at C5 — so each member's thin `reg_dist_*.f90` and its
`samp_reg` line land at C5 together with the sample routine (drag-in rule
unchanged). `reg_adiabatic.f90` remains the sole reg file of the K2 library
seed.

Present (migrated from the former interface_distributions/ in this same rebuild,
E14): `dist_boltzmann` ('boltzmann'), `dist_ebk` ('ebk'), `dist_incident`
('incident'). The remaining 9 members arrive at C5: surface_oscillator,
rotation, j (rotational-J temperature), barrier_excitation, thermostat,
glo_target, thermalize, normalmode, qmmicro.

Reference material (old library, read-only): the MQC propagation sites in
venus2026/src_VENUS/*.f — citation pass during C-phase implementation, no
skeleton-time citation.
