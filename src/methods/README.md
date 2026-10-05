# methods/ - the flat built-in method library

The library home of the generic methods (E14 spec §5, rebuilt from the former
interface_electronic/ + interface_distributions/): pre-stored, parameter-free method
packages and reusable algorithm components. Everything here is physics-generic code;
system specificity enters only through the method's own data channels at runtime (E12
ruling: system material lives in systems/<system>/; a method that embeds system
parameters or system knowledge does not belong here — see also the
knowledge-ownership audit: system/method-specific knowledge must not enter generic
code). E27: a SYSTEM-SPECIFIC method package (bespoke Hamiltonian + bespoke evolution
law) is physical-system container material - it lives in its systems/<system>/
folder, shares private modules with the container's H definition, and registers its
method row from a reg_*.f90 in that same folder (one registration site; the former
ws interchange between the two homes is retired - a historical term, docs/plans/S.md
E27).

## File prefixes (self-describing families)

| Prefix | Family | Registered by | Present / planned |
|---|---|---|---|
| `elec_*` | MQC method packages — {prop member, match member [, sync member]}, one `method_reg` row each | its thin `reg_<package>.f90` (assembly phase; in methods/ for generic packages, in systems/<system>/ for container-composed ones - E27) | `elec_adiabatic` (present — the all-no-op degenerate seed); `elec_tsh` / `elec_iesh` / `elec_ehrenfest` planned (R8-5) |
| `dist_*` | initial-state distribution members (the B′ sampler family; parameters via config — E12 wording carried, zero embedded system data) | thin `reg_dist_*.f90`, landed at C5 together with the sample routines (K2 E3 slimming ruling) | 12 members present (boltzmann, ebk, incident, surface_oscillator, rotation, j, barrier_excitation, thermostat, glo_target, thermalize, normalmode, wigner); qmmicro dropped (old-library defect variant, C5-R4) |
| `amp_*` | shared amplitude integrators — library modules, not packages (Corollary 1: FSSH and Ehrenfest share the same prop code; the method count drops multiplicatively) | none (components, consumed via use — level 2) | `amp_rk` / `amp_expm` planned |
| `reg_*` | thin registration modules — build material, not methods | (they are the registration) | `reg_adiabatic` (present — the sole reg file of the K2 seed) |

## Two exposure levels of the library ("ready-made semantics", spec §5)

**Level 1 — whole packages** (the standard path, zero code):
`config%electronic%method = 'fssh'` binds the assembled package at container_init; the
mqc_statistic / elec_prop facade dispatch lands directly in the package's match/prop
members. A generic package (methods/) binds any system; a container-composed package
(E27, systems/) binds its own folder's system. Dispatch chain of
mqc_statistic: facade table lookup -> match member (reads sta%s%a + sta%p + its own
data channels; g_ij from the library formula; rng-sampled hop; writes the held occ
via occ_set / the counter via n_hop_set / sta%p; optional force_eval
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
registering author's own knowledge, not the facade's. A member aborts at runtime with
a named error on dt semantics it cannot live with (E13 open point — fixed-step-only
laws; the shape lands with elec_tsh).

## The samp_* family (initial-state distributions — the B′ sampler registry)

`samp_*` members are pre-stored generic initial-state distribution methods: each
draws one fragment of the trajectory's initial state per its prescription
(Boltzmann normal-mode quantum numbers, fixed (n,J) EBK levels, surface-oscillator
stretch ...). E12 wording carried: the members are physics-generic code with zero
embedded system data — temperatures, frequencies, quantum numbers, and mode
switches all arrive as member params filled from config at runtime; system
material lives in systems/<system>/. The family was renamed `dist_*` -> `samp_*`
on 2026-09-22..28 (registration words unchanged); the incident channel left the
family for its own `inc_*` prefix (below).

Present — 10 members: `samp_boltzmann` / `samp_ebk` (migrated from the former
interface_distributions/ at the E14 rebuild) plus `samp_surface_oscillator`,
`samp_rotation`, `samp_j`, `samp_barrier_excitation`, `samp_glo_target`,
`samp_thermalize`, `samp_normalmode` (landed at C5; scheme-code table in
interface/sampler.f90; `samp_thermostat` was deleted when the thermostat became
the BATH provision, 2026-09-22), and `samp_wigner` (2026-09-21: the ground-state
Wigner normal-mode distribution — independent Gaussian q/p per mode with
virial-matched widths; the per-mode energy scale E_SCALE defaults to the
zero-point level and may be set below it for sub-ZPE sampling; criterion
program scripts/check_samp_wigner.f90). `qmmicro` is not
carried over (old-library defect variant, dropped by C5-R4).

## The inc_* incident family (2026-09-28 split) + beam_laws

The incident channel is THREE members, one per assembly paradigm, registered
under paradigm guards (exactly one armed per process): `inc_surface`
(`incident_surface` — one PROJECTILE fragment + substrate carriers; drop site
N_AIM 0/1/2, the fixed site AIM_X/AIM_Y [A, POSCAR frame]; on-site B_MAX = 0
under N_B = 0; A_LAT/cell commensurability check; full COM velocity - the
surface is the momentum sink), `inc_pair` (`incident_pair` — the gas-pair
branch verbatim, reduced-mass split), `inc_single` (`incident_single` — the
trivial no-op for single-fragment gas; the sampler dispatches uniformly).
The shared stochastic law (orientation/azimuth/aiming/energy draws, parameter
validation, COM placement/velocity assignment) lives in the level-2 component
`beam_laws.f90` (a `beam_t` instance per member — registry-family members never
depend on each other). Criterion programs: scripts/check_inc_surface.f90,
check_inc_pair.f90, check_inc_single.f90.

Reference material (old library, read-only): the MQC propagation sites in
venus2026/src_VENUS/*.f — citation pass during C-phase implementation, no
skeleton-time citation.
