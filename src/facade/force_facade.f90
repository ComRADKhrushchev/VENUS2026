!=====================================================================
! force_facade.f90 - L2 interface-A facade (E14 four-entry form): entry 1 - force_eval,
!                    (Q,P,s)->(F,E) (inside: the sum of class-A force members (additive);
!                    the unit-contract fulfillment point - x23.0605 x e_conv appears only
!                    here) + the class-A member registry + blackbox_init, the single
!                    assembly entry point spanning BOTH facades
! Reference: venus2026/DVDQ.f:10-22 (PDOT zeroing + DPESHON dispatch - the old form of additive composition);
!       src_TEST/interface_TEST.f90:12-16 (POTPRE member-init hook), :18-30
!       (POT0 energy conversion x23.0605), :32-44 (DPESHON gradient accumulation
!       PDOT(i)=PDOT(i)-g(i)*23.0605*C1 - the new facade conversion is character-for-character identical,
!       N1.1 bit-comparable); ENERGY.f:54-68 (T=sum p^2/2m/C1, V=POT0 convention, H=T+V);
!       docs/final_completion_plan.md §5 N1 (registration/assembly spec - directly in force; its
!       negotiation-table clause superseded by E13, see the E13 note below the registries)
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
!=====================================================================
module force_facade
   use consts,        only: e_conv
   use state,         only: state_t
   use config,        only: electronic
   use config_roster, only: roster
   use elec_ws,       only: ws_t, r_scalar, r_diab, r_adiab, r_empty
   use elec_facade,   only: ws_get, ws_reset, method_bind
   implicit none
   private
   public :: force_eval, pes_reg, blackbox_init, blackbox_info

   ! eV to kcal/mol conversion (unit-contract literal — T4: appears only in this file)
   real(8), parameter :: ev_kcal = 23.0605d0

   ! Class-A member (additive) interface, E10-unified: (q,mass) -> ws (invariant in both system and
   ! representation - single-surface scalar and multi-surface diabatic/adiabatic unified as ws fill modes;
   ! shape arrays stay out of the signature)
   ! (see previous line)
   abstract interface
      subroutine pes_proc_i(q, mass, ws)
         import :: ws_t
         real(8), intent(in)  :: q(:)     ! coordinates [Å] (all degrees of freedom, flattened)
         real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (used by stochastic/mass-law members;
                                          ! pure-PES members may ignore it)
         type(ws_t), intent(inout) :: ws  ! workspace (inout: the facade pre-writes the request field
                                          ! i_want(:) - the occupation set (E14: generalized from the
                                          ! scalar; length 1 = the degenerate form), the member reads
                                          ! it and fills per its own mode the
                                          ! delivery fields - scalar fills g_add/v_add;
                                          ! diab fills h_ij/dh_ij; adiab fills
                                          ! v_surf/g_surf/d_ij [the i_want rows];
                                          ! see the container contract in elec_ws.f90)
      end subroutine pes_proc_i
      ! Note: stochastic-law members call rng internally (N1 - the signature acknowledges stochastic laws;
      ! the single stream is guaranteed by rng)
   end interface

   ! Class-A member registry (additive table; registered in the assembly phase, read-only
   ! during evolution - E14: the class-B slot lives in elec_facade now)
   integer, parameter :: max_pes = 16     ! class-A member-table capacity [row] (main PES + additive stochastic laws; filled in C2)
   type :: pes_entry_t
      character(len=32) :: name = ''      ! member name (the config%electronic%pes_main assembly key)
      integer :: kind = r_scalar          ! member ws fill mode (r_* - registration data;
                                          ! E10: representation is data, not an interface fork)
      procedure(pes_proc_i), pointer, nopass :: proc => null() ! member evaluation routine
                                          ! (single signature - scalar/diab/adiab unified)
   end type pes_entry_t
   type(pes_entry_t) :: pes_tbl(max_pes)  ! class-A registry (module-private)
   integer :: n_pes = 0                   ! registered class-A row count [row]

   ! E13 (user ruling 2026-09-08): no capability-combination table lives here - triad_tbl and the
   ! m_* electron-law codes are withdrawn. The facade is a purely formal authority;
   ! physical-combination legality belongs to the members: the class-B member loud-stops on a
   ! ws%repr it does not support (E10 contract) and on dt semantics it cannot live with
   ! (fixed-step-only laws - e.g. hop probabilities g_ij(dt) losing their basis at adaptive
   ! sequence boundaries; shape fixed in phase C with elec_tsh). A facade-side enumeration of
   ! legal combinations would also collide with E11 drag-in: a blackbox author bringing a new
   ! combination would have to edit this file to add a row - exactly what "drag in and it works"
   ! forbids. The assembly-time checks kept in facade_init are purely formal: names resolve, the
   ! main-PES slot is filled, additive members are r_scalar (the additive channel of the container
   ! is scalar by construction), and contradictions derivable from the container contract alone
   ! (n_surf >= 2 with an r_scalar main PES) still loud-stop.

contains
   !------------------------------------------------------------------
   ! pes_reg(name, proc[, repr]) - register a class-A force member (additive: main PES + stochastic
   ! laws / constraint potentials; E10-unified single entry - representation via the repr data tag, no
   ! longer a registration fork; renamed facade_reg_pes -> pes_reg by E14 - noun_reg per the
   ! prop_reg/samp_reg convention)
   ! (see previous line)
   ! Reference: docs/final_completion_plan.md §5 N1 (registration spec); T2 (new member = new file +
   !       one registration line); E10 adjudication (single signature fully unified - E9 dual signatures withdrawn)
   subroutine pes_reg(name, proc, repr)
      character(len=*), intent(in) :: name  ! member name (assembly key - config%electronic%pes_main)
      procedure(pes_proc_i) :: proc         ! member evaluation routine (single signature - scalar/diab/
                                            ! adiab unified, delivered via the ws container)
      integer, intent(in), optional :: repr ! ws fill mode r_* (default r_scalar -
                                            ! single-surface members need zero registration-line
                                            ! changes; multi-surface main PES registers with
      ! r_diab/r_adiab)
      ! (see previous line)
      ! Pseudocode:
      ! 1. Duplicate name or full table -> loud-stop
      ! 2. Register the row {name, kind=repr (present-check defaults to r_scalar), proc} (executed in the
      !    assembly phase by the driver; multi-surface pes members arrive with the elec_tsh consumers)
   end subroutine pes_reg

   !------------------------------------------------------------------
   ! blackbox_init() - member assembly + purely formal validation spanning BOTH facades (E13 form
   !                 carried; E14: renamed facade_init -> blackbox_init - the single assembly
   !                 entry point of the whole blackbox; physical-combination legality remains
   !                 the members' own affair)
   ! Reference: docs/final_completion_plan.md §5 N1 (assembly spec; the negotiation-table clause is
   !       superseded by E13); docs/superpowers/specs/2026-09-08-facade-e14-design.md §3 (assembly-
   !       phase entries: names resolve, main-PES slot filled, method slot filled, additive members
   !       r_scalar, n_surf>=2 vs r_scalar contradiction, observable-vs-container delivery
   !       contradiction); src_TEST/interface_TEST.f90:12-16 (POTPRE-type member-init hook - the D3
   !       ghost parameter assembly site); docs/design_spec_completeness.md §3.2-2
   !       (type-checking of capability and method legality; illegal combinations rejected
   !        explicitly with both conflicting axes named)
   subroutine blackbox_init()
      ! Pseudocode:
      ! 1. Class-A assembly: freeze the main-PES member row from pes_tbl per config%electronic%pes_main
      !    + the additive member rows (additive member table - a config%electronic field left blank,
      !    filled in C2 from the old NGLO/NFC switches for grammar backward compatibility); an empty
      !    assembly table (no main-PES member) -> loud-stop.
      !    Multi-surface systems (n_surf >= 2): the main-PES row's ws fill mode must be r_diab/r_adiab
      !    (an r_scalar main PES with multiple surfaces -> loud-stop naming the conflict); additive
      !    members are restricted to r_scalar
      !    (stochastic laws / constraint potentials are independent of the electronic state - they feed
      !    only the nuclear-force channel; the multi-surface construct channel is carried by the main PES
      !    alone, see force_eval step 4')
      !    (see previous line)
      ! 2. Method-slot assembly (E14 - the old class-B slot is now a method PACKAGE in elec_facade):
      !    call method_bind(config%electronic%method) - a pure name lookup into method_tbl (formal; the
      !    empty-table / not-found / already-bound loud-stops live in method_bind itself, ratified - no
      !    duplicate here. The table lives in elec_facade; the init choreography stays here as the
      !    single assembly entry point)
      ! 3. Member-init hooks (POTPRE type): call each assembled member's init (parameters arrive through
      !    the member's own params - D3; the hook-slot wiring shape is finalized in C2)
      !    (see previous line)
      ! 4. Observable-vs-container formal check (spec §4): an observable column demanding a container
      !    no assembled package member delivers (a/rho declared undelivered by the bound method row)
      !    -> loud-stop naming the conflict (derivable from the container contract alone, E13);
      !    integrator x electron-law legality is NOT checked here (a law's own prerequisites are its
      !    author's knowledge - the member loud-stops at runtime, e.g. hop probabilities g_ij(dt)
      !    losing their basis at adaptive sequence boundaries)
   end subroutine blackbox_init

   !------------------------------------------------------------------
   ! force_eval(sta, e_tot) - entry 1 (spec §3), (Q,P,s)->(F,E): class-A members deliver via ws ->
   !                            facade composes -> converts (the unique unit-contract point, ev_kcal)
   !                            -> writes sta%f, e_tot; E = T(P) + V (the E10 unified path; renamed
   !                            facade_force -> force_eval by E14; dt removed from the signature -
   !                            the E9 channel re-anchored onto elec_prop)
   ! Reference: venus2026/DVDQ.f:10-22 (zeroing + dispatch); src_TEST/interface_TEST.f90:32-44
   !       (identical-form accumulation), :18-30 (energy convention); ENERGY.f:54-68 (T/H reporting
   !       convention); spec §3 entry-1 contract (as a side effect the ws of the current geometry is
   !       left in place for the later electronic phases)
   subroutine force_eval(sta, e_tot)
      type(state_t), intent(inout) :: sta  ! physical state (q/p/s/occ read in; f written out - the composed force)
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
      ! 2. Main-PES evaluation (the assembled main-PES member - one unified call): request the
      !    occupation set - ws%i_want(:) = occ_get() (the adiab-mode delivery request, E16:
      !    generalized from the scalar i_surf; length 1 = the single-active-surface degenerate
      !    form; the occupation set stays out of the class-A signature and is conveyed through
      !    the request field, E10); via the ws_get instance: proc(sta%q, roster%mass, ws) -
      !    the member fills ws per its own mode
      !    (scalar: g_add/v_add; diab: h_ij/dh_ij; adiab: v_surf/g_surf/
      !    d_ij [the i_want rows])
      ! 3. Fold scalar contributions into g(:)/v_sum: when the main PES is in scalar mode, fold its
      !    g_add/v_add;
      !    then for each additive member (restricted to r_scalar rows - stochastic laws / constraint
      !    potentials are independent of the electronic state, feeding only the nuclear-force channel):
      !    call proc(sta%q, roster%mass, ws) and fold in (accumulation in eV; the facade resets ws
      !    between members - deliveries do not stack, folding is done explicitly by the facade;
      !    stochastic-law members call rng internally; the multi-surface construct channel is carried
      !    by the main PES alone)
      !    (see previous line)
      !    (see previous line)
      ! 4. Unit-contract fulfillment point (T4 - unique in the whole library):
      !    after pre-zeroing sta%f, element-wise sta%f(i) = sta%f(i) - g(i)*ev_kcal*e_conv
      !    (character-for-character identical to the old DPESHON "PDOT(i)=PDOT(i)-g(i)*23.0605*C1" -
      !    N1.1 bit-comparable; e_conv is the old C1; how g is composed in multi-surface mode, see 4')
      !    (see previous line)
      ! 5. NO electron-law consumption here (E14: withdrawn to elec_prop / mqc_statistic in
      !    elec_facade) - the ws of the current geometry is left in place for the later phases
      !    (the spec §3 side-effect contract: what force_eval writes into the holder is what
      !    elec_prop / mqc_statistic next read)
      ! 6. Energy: v_pot = v_sum*ev_kcal [kcal/mol] (old POT0 :26 convention);
      !    t_kin = sum p^2/(2*roster%mass)/e_conv [kcal/mol] (per-atom fold -
      !    old ENERGY.f:58-66 convention); e_tot = t_kin + v_pot (old :68 H=T+V)
      ! ————————————————————————————————————————————————
      ! 4'. Multi-surface branch (main-PES row mode = r_diab/r_adiab - a composition swap of steps 3-6):
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
      !       new rows belongs to the next evaluation)
      !       (see previous line)
      !       (see previous line)
      !       (see previous line)
      !       (see previous line)
      !    b. v_pot convention: adiab -> v_surf(occ)*ev_kcal (the current-occupation surfaces); diab ->
      !       population-weighted or current-surface (class-B law-member convention, fixed by the
      !       consumer); scalar additive contributions fold in as usual;
      !       (see previous line)
      !       e_tot = t_kin + v_pot
      !    c. radau note: fixed-step only on multi-surface paths (E13: enforced inside the class-B
      !       member - not by a facade table) - adaptive substeps have no
      !       well-defined "one step" dt semantics; this path only ever sees fixed steps (same source as
      !       the propagator header note)
      !       (see previous line)
   end subroutine force_eval

   !------------------------------------------------------------------
   ! blackbox_info() - diagnostic output: member list and assembly state (renamed
   !                   facade_info -> blackbox_info by E14)
   ! Reference: docs/final_completion_plan.md §5 N1 (facade_info - member list / capability levels)
   subroutine blackbox_info()
      ! Pseudocode:
      ! 1. Print the class-A registered/assembled member list (name + main-PES/additive tag + ws fill
      !    mode r_*) and the assembled method-package name (the bound row of the elec_facade method
      !    table - E14: the old class-B member-name slot is gone)
      !    (see previous line)
      ! 2. Print the current capability combination: H-capability n_surf x main-PES ws mode x method
      !    x integrator (from the config sources - descriptive output only; legality is
      !    the members' own affair, E13)
      !    (see previous line)
      ! 3. Print the unit-contract constants (ev_kcal/e_conv) and the assembly state - for debugging
      !    assembly errors
   end subroutine blackbox_info
end module force_facade
