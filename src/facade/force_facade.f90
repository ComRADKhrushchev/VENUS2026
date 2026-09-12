!=====================================================================
! force_facade.f90 - L2 interface-A facade (E14 four-entry form, E17 container form): entry 1
!                    - force_eval, (Q,P,s)->(F,E) (the unit-contract fulfillment point -
!                    x23.0605 x e_conv appears only here) + the container single-slot
!                    binding + container_init, the single assembly entry point spanning
!                    BOTH facades
! Reference: venus2026/DVDQ.f:10-22 (PDOT zeroing + DPESHON dispatch - the old form of additive composition);
!       src_TEST/interface_TEST.f90:12-16 (POTPRE member-init hook), :18-30
!       (POT0 energy conversion x23.0605), :32-44 (DPESHON gradient accumulation
!       PDOT(i)=PDOT(i)-g(i)*23.0605*C1 - the new facade conversion is character-for-character identical,
!       N1.1 bit-comparable); ENERGY.f:54-68 (T=sum p^2/2m/C1, V=POT0 convention, H=T+V);
!       docs/final_completion_plan.md §5 N1 (registration/assembly spec - directly in force; its
!       negotiation-table clause superseded by E13, see the E13 note below the slot)
! Member unit contract: class-A members deliver force kernels/energies in eV (eV is the exchange unit -
! following the old interface_TEST v/g chain); eV -> internal-unit conversion is centralized in this
! facade (the T4 grep requirement that 23.0605 appear only in this file - hence the "F_add code units" note
! in the N1 member-signature spec reads as "after facade conversion")
! (see previous line)
! (see previous line)
! Member tables are filled in C2 (main PES + additive stochastic-law members; E12 ruling: pes members
! are system material - each systems/<system>/ folder carries self-contained functional form + parameters,
! the "six pes_* files land in interface_electronic/" plan of K2 E3 is superseded; the facade itself is
! untouched)
!
! -- E10 workspace revision (user adjudication 2026-09-07, three axioms enacted - amending E9) --
!   A Interface invariance: the abstract interfaces are invariant across systems - signatures carry only
!     the container ws_t (facade/elec_ws.f90), never shape arrays like h_ij/d_ij; the representation
!     (diabatic/adiabatic) is demoted from a signature fork to the ws%repr data tag. E9's pes_diab_i/
!     pes_adiab_i dual signatures and pes_adiab_i's i_surf argument are hereby withdrawn (representation
!     and electronic-state knowledge stay out of the class-A signature).
!     (see previous line)
!     (see previous line)
!     (see previous line)
!   B Black-box property: class A + class B form one system black box - outside the box only
!     (Q,P,s[,dt])->(F,E) is visible (zero change to the facade_force signature); ws is a facade_force
!     local variable, declared inside the facade and dying inside the facade.
!     (see previous line)
!     (see previous line)
!   C Thing-in-itself principle: representation constructs such as H_ij/d_ij/V_i have no thing-in-itself
!     status - they do not enter state_t, do not leave the box, and live for one step (reset and rebuilt
!     at every evaluation).
!     (see previous line)
! E9 retained assets: elec_t amplitudes/hop counter, the class-B dt channel (the negotiation table
! triad_tbl is withdrawn by E13 below)
!
! -- E13 formal-facade revision (user ruling 2026-09-08) --
!   The facade is a purely FORMAL authority: signature/unit/container contracts and name->pointer
!   dispatch - nothing more. Physical-combination legality (which electron law tolerates which
!   representation or integrator - e.g. the hop probability g_ij(dt) losing its basis at adaptive
!   sequence boundaries) is the members' own knowledge: the class-B member loud-stops on a ws%repr
!   it does not support (E10 contract, unchanged) and on dt semantics it cannot live with (shape
!   fixed in phase C with elec_tsh). Rationale: blackbox authors are assumed to know their system,
!   not this program - a facade-side legal-combination table is knowledge displaced from its owners,
!   and it collides with E11 drag-in (a new combination would require editing this file to add a
!   row). The blackbox machine does three things (legislated wording): compute forces from the
!   state (class A), evolve the electronic state (class B), and feed the electronic state back
!   into the state to enforce quantum statistics (class B - surface switching, momentum rescaling);
!   facade_force is the orchestration of the three.
!
! -- E14 four-entry rebuild (user rulings 2026-09-08; docs/superpowers/specs/
!    2026-09-08-facade-e14-design.md §3) --
!   The E13 closing sentence "facade_force is the orchestration of the three" is withdrawn: the
!   three things the blackbox machine does (compute forces from the state (class A), evolve the
!   electronic state, feed the electronic state back into the state to enforce quantum statistics)
!   are now three DRIVER-orchestrated entries - force_eval (this module, entry 1), elec_prop /
!   mqc_statistic (elec_facade, entries 2-3; spec §8 choreography: prop_step -> elec_prop ->
!   mqc_statistic -> rec_frame -> pred_eval). Consequences here: facade_force-the-orchestrator is
!   withdrawn (force_eval is a pure force entry - no electron-law consumption inside); the class-B
!   slot (elec_proc_i / facade_reg_elec / elec_name / elec_proc) is excised from this module -
!   that content now lives in elec_facade as method packages (method_reg / method_bind); the dt
!   argument is removed from the entry signature (the E9 channel re-anchored onto elec_prop - the
!   hop probability g_ij(dt) belongs to the electronic phases, not to the force evaluation); and ws
!   is no longer a facade_force local variable - force_eval fills the elec_facade-held instance
!   through ws_get and leaves it in place for the later electronic phases (one-integration-step
!   lifetime; E10 axioms B/C as amended by E14). Entry names drop the facade_ prefix (spec §3
!   naming statute - the module names the place, the entry names the action): facade_reg_pes ->
!   pes_reg, facade_init -> blackbox_init, facade_force -> force_eval, facade_info ->
!   blackbox_info; the module keeps its *_facade name (architectural position).
!
! -- E17 container revision (review-gate rulings 2026-09-12: the singleton argument + the
!    container-is-the-concrete-Hamiltonian axiom) --
!   (1) Terminology: "blackbox" renamed CONTAINER throughout the live code and the live
!       legislation comments (blackbox_init -> container_init). The historical ruling blocks
!       above keep their original wording - the E14 naming ledger records the older map; a
!       grep for "blackbox" therefore still hits the historical blocks (audit grep basis).
!   (2) The singleton argument: exactly ONE container exists per run - the registry/lookup
!       indirection built for multi-instance selection is overengineering on the class-A
!       side. pes_reg / pes_proc_i (the public contract interface) / blackbox_info are
!       hereby WITHDRAWN (the diagnostic entry dies with the registry it printed).
!   (3) Class-A composition ownership returns to the container (E13 completed): the sum over
!       members (main PES + additive stochastic laws) executes INSIDE the container's single
!       aggregate entry, bound at assembly through container_bind_pes(proc [, h_deliver]) -
!       a single slot, no name table (a module-private abstract interface remains as the
!       pointer's type declaration - a Fortran implementation detail, no longer a public
!       contract).
!   (4) The container's essence (legislated axiom, user wording 2026-09-12): the container IS
!       the system's concrete (generalized) Hamiltonian / dynamics generator. Its members are
!       H's evaluators (forces = -grad H, energies = H's values); the ws deliveries are H's
!       electronic part in concrete representations; the system-definition files are H's
!       parameter set (masses enter T_nuc; composition/cell fix H_el's functional form); the
!       stochastic/dissipative members (GLO/Langevin ghost/LDFA) are its generalized
!       (non-symplectic) part. The lone exception: termination primitives - problem
!       knowledge, not H knowledge, carried by the container only because equally
!       system-specific. Corollary - the ownership criterion: "would this knowledge change
!       with a different system?" yes -> container; no -> the generic side (algorithms /
!       formal legislation). Program = generic algorithm machine (integration / propagation /
!       matching / sampling / recording) + container (the concrete generalized Hamiltonian);
!       state_t is their conversation site (Q/P/s/F/t).
!   (5) H-form declaration channel (E14's deferred deliverable-declaration field activated
!       early): the method row gains h_form - the method's internal Hamiltonian PRESET (which
!       H forms it consumes: comma-list over the r_* mode names; empty = consumes no H
!       constructs, the degenerate method). The container declares its delivery set through
!       container_bind_pes's h_deliver argument; container_init matches demand
!       subset-of-supply at ASSEMBLY time (both sides named in the loud-stop) - the E10
!       runtime repr-mismatch stop moves earlier, it does not go away.
!   (6) ws construction ownership restated: the ws CONTENT is constructed entirely inside the
!       container (fill modes self-declared by the filling member); the generic side owns the
!       request channel (i_want <- occ) and the formal container contract only. The ws_t type
!       stays in the generic facade as FORMAL legislation (E13 authority - the facade
!       legislates form; this is method-spectrum knowledge, not system knowledge; a new H
!       form = a generic-program version step, the same precedent as a new predicate
!       primitive class).
!=====================================================================
module force_facade
   use consts,        only: e_conv
   use state,         only: state_t
   use config,        only: electronic
   use config_roster, only: roster
   use elec_ws,       only: ws_t, r_scalar, r_diab, r_adiab, r_empty
   use elec_facade,   only: ws_get, ws_reset, method_bind, occ_get, method_h_form
   implicit none
   private
   public :: force_eval, container_bind_pes, container_init

   ! eV to kcal/mol conversion (unit-contract literal — T4: appears only in this file)
   real(8), parameter :: ev_kcal = 23.0605d0

   ! The container aggregate entry - E17 single-slot contract (module-PRIVATE interface:
   ! the pointer's type declaration, a Fortran implementation detail - no longer a public
   ! member contract; the container's ONE class-A entry, the member sum executing inside
   ! it, E17(3)). Delivery contract: the scalar channel g_add/v_add = the container-side
   ! SUM in eV (the main PES when in scalar mode + every additive stochastic-law /
   ! constraint member; multi-surface mode: the additive sum alone - the main PES rides
   ! the construct channel); the construct channel = the main PES's multi-surface
   ! deliveries (h_ij/dh_ij or v_surf/g_surf/d_ij [the i_want rows]) per its
   ! self-declared fill mode
   abstract interface
      subroutine container_pes_i(q, mass, ws)
         import :: ws_t
         real(8), intent(in)  :: q(:)     ! coordinates [Å] (all degrees of freedom, flattened)
         real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (used by stochastic/mass-law
                                          ! members; pure-PES members may ignore it)
         type(ws_t), intent(inout) :: ws  ! workspace (inout: the facade pre-writes the request
                                          ! field i_want(:) - the occupation set (E14/E16);
                                          ! the container reads it and fills per its own mode
                                          ! - see the container contract in elec_ws.f90)
      end subroutine container_pes_i
   end interface
   ! Note: stochastic-law members call rng internally (N1 - the signature acknowledged
   ! stochastic laws; the single stream is guaranteed by rng) - this holds inside the
   ! container exactly as it held per member

   ! The bound slot (E17: a single slot replaces pes_tbl/pes_entry_t/n_pes - no registry,
   ! no name lookup; one container exists per run, the singleton argument)
   procedure(container_pes_i), pointer :: container_pes => null()
   character(len=32) :: h_deliver = ''  ! the container's declared H delivery set (comma-list
                                        ! of r_* mode names, e.g. 'scalar' / 'diab' /
                                        ! 'scalar,adiab'; the assembly-time match target of
                                        ! the method row's h_form, E17(5))

   ! E13/E17: the facade is a purely FORMAL authority - no capability-combination table
   ! lives here. The assembly-time checks kept in container_init are purely formal: the
   ! slot is bound; the method name resolves; n_surf >= 2 without a multi-surface mode in
   ! h_deliver loud-stops (derivable from the container contract alone); the method row's
   ! h_form demand is a subset of the container's h_deliver supply (E17(5)). Physical-
   ! combination legality beyond form (which integrator a law tolerates, dt semantics)
   ! remains the members' own knowledge - they loud-stop at runtime (E13 carried).

contains
   !------------------------------------------------------------------
   ! container_bind_pes(proc [, h_deliver]) - E17 single-slot assembly injection: bind the
   !                 container's aggregate class-A entry + its declared H delivery set.
   !                 No name table, no registry (the singleton argument); the caller is
   !                 the container's own reg module (reg_<system>_all), aggregated by the
   !                 K3 prebuild into registry_gen - drag-in unchanged (a new system = a
   !                 new folder with its own reg, zero edits to existing files, E11)
   ! Reference: this file's E17 header block; E11 (drag-in - the reg file is the atomic
   !       registration unit); E13 (formal authority - the injection checks slot hygiene
   !       only); T2 (the upgraded criterion: a new system = a new folder)
   subroutine container_bind_pes(proc, h_deliver)
      procedure(container_pes_i) :: proc         ! the container aggregate entry (the single
                                                 ! class-A slot)
      character(len=*), intent(in), optional :: h_deliver ! the declared delivery set
                                                          ! (comma-list of r_* names; default
                                                          ! '' = scalar only - the degenerate
                                                          ! single-surface container)
      ! Pseudocode:
      ! 1. Already bound (associated(container_pes)) -> loud-stop (container_bind_pes ran
      !    twice - an assembly-choreography error, formal)
      ! 2. Bind the slot, store h_deliver (present-check defaults to '') - executed once
      !    in the assembly phase by the container's reg module
   end subroutine container_bind_pes

   !------------------------------------------------------------------
   ! container_init() - assembly + purely formal validation spanning BOTH facades (E13
   !                 form carried through E14/E17; renamed facade_init -> blackbox_init
   !                 (E14) -> container_init (E17); the single assembly entry point of
   !                 the whole container; composition legality is the container's own
   !                 affair - it already executed inside the bound aggregate entry)
   ! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §3 (assembly-phase
   !       entries - the E17 form: slot bound, H forms matched, method slot filled,
   !       observable-vs-container delivery contradiction); src_TEST/interface_TEST.f90:12-16
   !       (POTPRE-type member-init hook - now container-internal, executed by its reg
   !       module); docs/design_spec_completeness.md §3.2-2
   subroutine container_init()
      ! Pseudocode:
      ! 1. Slot check: container_pes unbound (container_bind_pes not yet run - an empty
      !    container folder / a missing reg module) -> loud-stop naming the missing binding
      ! 2. Container-contract formal checks: multi-surface system (config%electronic%n_surf
      !    >= 2) while h_deliver names no multi-surface mode (diab/adiab) -> loud-stop
      !    naming the conflict (E13 kept: derivable from the container contract alone)
      ! 3. Method-slot assembly (E14): call method_bind(config%electronic%method) - a pure
      !    name lookup into method_tbl (the empty-table / not-found / already-bound
      !    loud-stops live in method_bind itself - no duplicate here; the table lives in
      !    elec_facade, the init choreography stays here as the single assembly entry
      !    point)
      ! 4. H-form match (E17(5)): parse the bound method row's h_form (method_h_form())
      !    and the container's h_deliver; h_form must be a subset of h_deliver -> else
      !    loud-stop naming BOTH sides (the method's Hamiltonian preset vs the container's
      !    delivery set); an empty h_form (the degenerate method) always passes
      ! 5. Container-side member init: executed by the container's reg module around its
      !    container_bind_pes call (POTPRE-type hooks are container-internal now - the
      !    facade no longer choreographs them; parameters arrive through the members' own
      !    params, D3)
      ! 6. Observable-vs-container formal check (spec §4): an observable column demanding
      !    a container no assembled package member delivers (a/rho declared undelivered by
      !    the bound method row) -> loud-stop naming the conflict (derivable from the
      !    container contract alone, E13); integrator x electron-law legality is NOT
      !    checked here (a law's own prerequisites are its author's knowledge - the member
      !    loud-stops at runtime, e.g. hop probabilities g_ij(dt) losing their basis at
      !    adaptive sequence boundaries)
   end subroutine container_init

   !------------------------------------------------------------------
   ! force_eval(sta, e_tot) - entry 1 (spec §3), (Q,P,s)->(F,E): the container aggregate
   !                            entry delivers via ws -> facade converts (the unique
   !                            unit-contract point, ev_kcal) -> writes sta%f, e_tot;
   !                            E = T(P) + V (the E10 unified path; renamed
   !                            facade_force -> force_eval by E14; dt removed from the
   !                            signature - the E9 channel re-anchored onto elec_prop)
   ! Reference: venus2026/DVDQ.f:10-22 (zeroing + dispatch); src_TEST/interface_TEST.f90:32-44
   !       (identical-form accumulation), :18-30 (energy convention); ENERGY.f:54-68 (T/H reporting
   !       convention); spec §3 entry-1 contract (as a side effect the ws of the current geometry is
   !       left in place for the later electronic phases)
   subroutine force_eval(sta, e_tot)
      type(state_t), intent(inout) :: sta  ! physical state (q/p/s read in; f written out - the composed force)
      real(8), intent(out) :: e_tot        ! total energy H = T(P)+V [kcal/mol] (N1 convention;
                                           ! multi-surface mode V = current surface / law-member convention, see 4')
      ! Workspace (E10 axioms B/C as amended by E14): NOT a local variable of this routine anymore -
      ! the elec_facade-held instance, obtained through ws_get (the fill site), composed into here,
      ! and left in place for elec_prop / mqc_statistic (one-integration-step lifetime; still
      ! invisible outside the box - never an argument of any public signature, zero footprint in
      ! the force_eval signature)
      ! Pseudocode:
      ! 1. Zero the force-kernel workspace g(:), v_sum = 0 (old DVDQ.f:12-14 PDOT zeroing site);
      !    ws_reset() (the head reset - every evaluation begins from the empty container; with
      !    elec_facade's step-end reset the two sites together make the one-step lifetime robust to
      !    skipped phases); dt no longer appears here (E9 channel withdrawn - the hop probability
      !    g_ij(dt) is consumed in elec_prop)
      ! 2. Unbound slot (container_pes not associated - container_init not yet run) ->
      !    loud-stop (an evolution-phase entry called before assembly - assembly-order
      !    error, formal)
      ! 3. ONE call into the container (E17: the member sum - main PES + additive
      !    stochastic laws - executes INSIDE it): request the occupation set -
      !    ws%i_want(:) = occ_get() (the adiab-mode delivery request, E16: generalized
      !    from the scalar i_surf; length 1 = the single-active-surface degenerate form;
      !    the occupation set stays out of the class-A signature and is conveyed through
      !    the request field, E10); via the ws_get instance:
      !    call container_pes(sta%q, roster%mass, ws) - the container fills ws per its
      !    own mode (scalar channel: the SUMMED g_add/v_add [eV]; construct channel:
      !    h_ij/dh_ij or v_surf/g_surf/d_ij [the i_want rows]; inter-member ws
      !    reset/folding is container-internal)
      ! 4. Fold the scalar channel: g(:) = ws%g_add, v_sum = ws%v_add (the container
      !    already summed - the facade folds ONCE; no member loop remains here)
      ! 5. Unit-contract fulfillment point (T4 - unique in the whole library):
      !    after pre-zeroing sta%f, element-wise sta%f(i) = sta%f(i) - g(i)*ev_kcal*e_conv
      !    (character-for-character identical to the old DPESHON "PDOT(i)=PDOT(i)-g(i)*23.0605*C1" -
      !    N1.1 bit-comparable; e_conv is the old C1; how g is composed in multi-surface mode, see 4')
      !    (see previous line)
      ! 6. NO electron-law consumption here (E14: withdrawn to elec_prop / mqc_statistic in
      !    elec_facade) - the ws of the current geometry is left in place for the later phases
      !    (the spec §3 side-effect contract: what force_eval writes into the holder is what
      !    elec_prop / mqc_statistic next read)
      ! 7. Energy: v_pot = v_sum*ev_kcal [kcal/mol] (old POT0 :26 convention);
      !    t_kin = sum p^2/(2*roster%mass)/e_conv [kcal/mol] (per-atom fold -
      !    old ENERGY.f:58-66 convention); e_tot = t_kin + v_pot (old :68 H=T+V)
      ! ————————————————————————————————————————————————
      ! 4'. Multi-surface branch (the container's declared mode = diab/adiab - a
      !     composition swap of steps 4-7):
      !    a. Force-kernel composition: r_diab -> diagonalize ws%h_ij via utils/linalg to compute the
      !       adiabatic view on the spot (eigenvalues = per-surface V, eigenvectors = representation
      !       transform U); the current-occupation adiabatic force kernel takes the transformed
      !       current rows; the dh_ij.U product needed for the hop probability g_ij(dt) is fetched by
      !       the electronic phases themselves (the facade does not pre-build it; the diagonalization
      !       itself is system-independent generic physics and lawfully stays in the facade - not
      !       system-specific code. Delivery of U: an adiab-view cache field in ws, written by the
      !       facade and read-only to the electronic phases (one-step lifetime, axiom C intact),
      !       or a re-diagonalization inside the consuming member - fixed in phase C with the
      !       elec_tsh consumer, noting the double O(n_surf^3) cost of the latter; the consumer is
      !       now elec_prop / mqc_statistic, no longer an in-force_eval call); r_adiab -> g(:) =
      !       fold of the occupation-set rows ws%g_surf(:, occ) (the current-occupation force kernel,
      !       E14 - generalized from the single current-surface row; length 1 = the degenerate
      !       single-row form; after a hop mqc_statistic updates the held occ (E16), recomputing with
      !       new rows belongs to the next evaluation); in BOTH modes the additive scalar kernel
      !       g_add (the container-side additive sum) folds on top of the construct-channel kernel -
      !       stochastic laws / constraint potentials are electronic-state-independent, they feed
      !       only the nuclear-force channel
      !       (see previous line)
      !       (see previous line)
      !       (see previous line)
      !       (see previous line)
      !    b. v_pot convention: adiab -> v_surf(occ)*ev_kcal (the current-occupation surfaces); diab ->
      !       population-weighted or current-surface (class-B law-member convention, fixed by the
      !       consumer); the additive scalar v_add folds in as usual;
      !       (see previous line)
      !       e_tot = t_kin + v_pot
      !    c. radau note: fixed-step only on multi-surface paths (E13: enforced inside the class-B
      !       member - not by a facade table) - adaptive substeps have no
      !       well-defined "one step" dt semantics; this path only ever sees fixed steps (same source as
      !       the propagator header note)
      !       (see previous line)
   end subroutine force_eval
end module force_facade
