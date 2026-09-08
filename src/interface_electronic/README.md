# interface_electronic/ - class-B electron-law members (the electronic propagation axis)

Home of the class-B members of spec §3.2-2 — the quantum-classical electron propagation laws
(`elec_*` prefix). These are parameter-free generic algorithms (the same FSSH/IESH machinery for
any system; system specificity enters only through the `ws` container and `state` at runtime), so
they are physics-generic code, not system blackbox material.

| Axis | Prefix | Interface | Present / planned |
|---|---|---|---|
| Electron propagators (quantum-classical methods) | elec_* | Class B: facade_reg_elec (single select). Signature elec_proc_i(s, f, dt, ws) - E9 adds the step (needed for the hop probability g_ij(dt)), E10 adds the workspace (read-only consumption: interpret the container per ws%repr; a repr the member does not support -> loud-stop). Adiabatic = degenerate point (ds/dt == 0, identity force rule, ignores dt and ws); the interface is well-defined at the degenerate point | elec_adiabatic (present); elec_tsh / elec_iesh later (fixed-step integrators only, per triad_tbl) |
| H-field providers (additive contribution family) | pes_* | **Not hosted here (E12)** — class-A members are system-specific: each systems/<system>/ folder carries its own complete members (self-contained functional form + parameters, copy-and-modify friendly). The class-A axis itself is unchanged: facade_reg_pes single entry (optional repr=r_scalar [default]/r_diab/r_adiab - E10: representation is registration data, not a registration fork), signature pes_proc_i(q, mass, ws), members deliver in eV per their own ws fill mode (scalar fills g_add/v_add; diab fills H_ij + matrix-element gradients; adiab fills V_i/force kernels + the requested derivative-coupling row via ws%i_want); the x23.0605xC1 conversion is centralized in the facade (unit contract, master plan §5 N1); multi-surface main PES is not additively summed | built from C2 on as self-contained members of the example systems (systems/<system>/) |

**The two axes are independent** (any PES x any method is a legal pairing, a Cartesian-product
structure), guaranteed by code mechanisms, orthogonal to directory grouping: the modules are
mutually independent (pes_* and elec_* share zero references), the two registries
(reg_pes / reg_elec are separate), and the facade_init capability negotiation (illegal
combinations loud-stop, naming both conflicting axes).

**E9 multi-surface upgrade -> E10 workspace revision (2026-09-07)**: E9 introduced the class-A
dual signature (diabatic/adiabatic) + the class-B dt argument + the negotiation table as data
triad_tbl (four rows: current single-surface adiabatic + multi-surface tsh/iesh pre-allowed).
E10's axiomatic correction (three axioms: interface invariance / black-box property /
thing-in-itself principle) - the class-A signature is fully unified as pes_proc_i(q, mass, ws):
representation constructs (H_ij/d_ij/V_i) are housed in the ws_t container (facade/elec_ws.f90 -
representation demoted from a signature fork to the ws%repr data tag); they do not enter state_t,
do not leave the facade box, and live for one step; the dual registration routines are withdrawn
(facade_reg_pes is a single entry with optional repr). **radau stays false on multi-surface rows**
unchanged (adaptive sequence boundaries have no well-defined "one step", the hop probability
g_ij(dt) loses its basis, so nonadiabatic only pairs with fixed-step integrators).

**E12 single-layer self-contained blackbox (2026-09-08, amends E9's "both axes in one home" for
class A)**: the program splits into system-specific code (the system blackbox) and physics-generic
code. Class-A potential members embed system parameters (functional form + numbers), so they are
system material: no generic method-skeleton directory is kept ("parameterless LEPS/EMT-NN classes
here, parameterized instances under systems/" was rejected) — each systems/<system>/ folder carries
complete, self-contained members. This directory therefore hosts class B only; the E3 note "the six
pes_* files land here in C2" is superseded. interface_distributions/ is unaffected (dist_* members
take parameters via config at runtime — zero embedded system data, legitimately generic).

**E11 drag-in build (2026-09-08)**: once K3 lands, a new system = a new folder under systems/
carrying its own reg_<system>_all module - the prebuild-generated registry_gen.f90 aggregates the
registration and the driver keeps zero registration lines ("drag in and it works": registration
demoted from hand-written code to a build artifact). Plugging in a *generic* member (this
directory) remains: a new file + one registration line in the assembly phase, zero touch to
existing code (the T2 criterion exercise point).

Reference material (old library, read-only): interface_TEST.f90 / EMT-NN / RST / VERLET.f90
inlined forces.
