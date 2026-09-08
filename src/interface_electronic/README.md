# interface_electronic/ - extensible electronic-layer members (both axes in one home)

The interchangeable members of the two independent electronic-layer axes of spec §3.2-2 live in
this directory, distinguished by file prefix:

| Axis | Prefix | Interface | Present / planned |
|---|---|---|---|
| H-field providers (additive contribution family) | pes_* | Class A: facade_reg_pes (single entry, optional repr=r_scalar [default]/r_diab/r_adiab - E10 unification: representation is registration data, not a registration fork). Signature pes_proc_i(q, mass, ws) - members deliver per their own ws fill mode: scalar fills g_add [eV/Å]/v_add [eV]; diab fills the H matrix + matrix-element gradients (eV, eV/Å); adiab fills per-surface V_i/force kernels + the current-surface derivative-coupling row d_i,. (via the ws%i_want request - i_surf stays out of the signature). The x23.0605xC1 conversion is centralized in the facade (unit contract, master plan §5 N1); multi-surface main PES is not additively summed | Six in C2: reference / emt_nn / rst / glo_spring / langevin_ghost / ldfa; multi-surface members arrive with the elec_tsh/elec_iesh consumers |
| Electron propagators (quantum-classical methods) | elec_* | Class B: facade_reg_elec (single select). Signature elec_proc_i(s, f, dt, ws) - E9 adds the step (needed for the hop probability g_ij(dt)), E10 adds the workspace (read-only consumption: interpret the container per ws%repr; a repr the member does not support -> loud-stop). Adiabatic = degenerate point (ds/dt == 0, identity force rule, ignores dt and ws); the interface is well-defined at the degenerate point | elec_adiabatic (first); IESH / FSSH will live here later |

**The two axes are independent** (any PES x any method is a legal pairing, a Cartesian-product
structure), guaranteed by code mechanisms,
orthogonal to directory grouping: the modules are mutually independent (pes_* and elec_* share
zero references), the two registries
(reg_pes / reg_elec are separate), and the facade_init capability negotiation (illegal
combinations loud-stop,
naming both conflicting axes).

**E9 multi-surface upgrade -> E10 workspace revision (2026-09-07)**: E9 introduced the class-A
dual signature (diabatic/
adiabatic) + the class-B dt argument + the negotiation table as data triad_tbl (four rows:
current single-surface adiabatic + multi-surface
tsh/iesh pre-allowed). E10's axiomatic correction (three axioms: interface invariance /
black-box property / thing-in-itself principle)
- the class-A signature is fully unified as pes_proc_i(q, mass, ws): representation constructs
(H_ij/d_ij/V_i)
are housed in the ws_t container (facade/elec_ws.f90 - representation demoted from a signature
fork to the ws%repr data
tag); they do not enter state_t, do not leave the facade box, and live for one step; the dual
registration routines are withdrawn
(facade_reg_pes is a single entry with optional repr). **radau stays false on multi-surface
rows** unchanged (adaptive
sequence boundaries have no well-defined "one step", the hop probability g_ij(dt) loses its
basis, so nonadiabatic only pairs with fixed-step integrators).

Plugging in = a new file + one registration line in the assembly phase, zero touch to existing
code (the T2 criterion exercise point).

Reference material (old library, read-only): interface_TEST.f90 / EMT-NN / RST / VERLET.f90
inlined forces.
