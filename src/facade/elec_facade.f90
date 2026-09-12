!=====================================================================
! elec_facade.f90 - L2 electronic-phase facade: entries 2-3 of the E14 four-entry blackbox
!                   contract (elec_prop = §4.1 layer 2 electronic propagation; mqc_statistic =
!                   §4.1 layer 3 state-statistics matching) + the method-package registry +
!                   the ws instance holder
! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §3 (the four-entry
!       contract, entries 2-3; ws holding and the endpoint-force contract), §4 (dual ontic
!       containers + the sync contract), §5 (the built-in method library: method-package
!       registration, the two-level data flow, Corollary 1);
!       docs/design_spec_completeness.md §4.1 (three-layer decomposition — layer semantics:
!       electronic propagation advances the electronic variables given Q(t) — amplitude
!       ODE / mapping flow / analytic elimination; state-statistics matching reads the
!       electronic statistics and writes the classical + nuclear sides, and is the sole
!       consumption point of the random stream during the evolution phase);
!       old-library propagation sites venus2026/src_VENUS/*.f (citation pass during C-phase
!       implementation — no skeleton-time citation)
!
! -- E14 four-entry rebuild (user rulings 2026-09-08; spec §3) --
!   The E13 sentence "facade_force is the orchestration of the three" is withdrawn: the
!   blackbox contract becomes four entries and the step choreography is explicit in the
!   driver (prop_step -> elec_prop -> mqc_statistic -> rec_frame -> pred_eval, spec §8).
!   This module owns the two electronic-phase entries; the dt argument retired from
!   facade_force lands on elec_prop (the E9 channel re-anchored — the hop probability
!   g_ij(dt) belongs to the electronic phases, not to the force evaluation).
!   (see previous line)
!
! -- ws holding (E10 axioms B/C amended by E14; spec §3) --
!   The ws instance is this module's module-private variable (saved, not a dummy): filled
!   by force_eval (through ws_get — the fill site), consumed by elec_prop / mqc_statistic,
!   reset to r_empty at step end (ws_reset). One-integration-step lifetime (extended from
!   E10's one evaluation), still invisible outside the box — ws never appears in any public
!   signature. The E14 motivation is exactly this read: the matching layer needs d_ij for
!   the momentum projection after a facade-force-local ws would already have died.
!   (see previous line)
!
! -- E13 carried: the facade is a purely FORMAL authority --
!   Name->pointer dispatch + container/sync contracts — nothing more. No capability-
!   combination table lives here: which prop member pairs with which match member is the
!   registering author's own knowledge (one method_reg call vouches for the whole
!   package); which ws%repr a member eats and which dt semantics it can live with are the
!   member's own knowledge (unsupported repr -> member loud-stop, E10 contract;
!   dt-semantics loud-stop — E13 open point, carried to phase C with elec_tsh). A
!   facade-side legality table would also collide with E11 drag-in (a new combination
!   would require editing this file to add a row — exactly what "drag in and it works"
!   forbids).
!   (see previous line)
!
! -- Method packages (spec §5, Corollary 1 realized in-library) --
!   A method = {prop member, match member [, sync member]} — the row is the atomic
!   registration unit; 'ehrenfest' and 'fssh' share a prop member and differ only in
!   match (the method count drops multiplicatively). Library-level reg_*.f90 files under
!   methods/ call method_reg; prebuild aggregates every reg_*.f90 under methods/ and
!   systems/ into registry_gen (E11 drag-in: a new method = one new file + one thin reg
!   file, zero edits to existing files). Assembly is a pure name lookup:
!   method_bind(config%electronic%method) — the system author writes zero method
!   registration lines (the method arrives as a config key).
!   (see previous line)
!
! -- E16 statistics holder (review-gate ruling 2026-09-10) --
!   The occupation vector occ(:) and the hop counter n_hop relocated here from state_t:
!   sampled constructs have no thing-in-itself status (the E10 axiom family extended —
!   occ is the discrete output of mqc_statistic's stochastic sampling, n_hop a pure
!   observable accumulator; state_t is purely Q/P/s/F/t). Run-lifetime module-private
!   holder, the ws sibling (ws = one-step representation construct, occ/n_hop = run-time
!   statistical constructs — same privacy axioms, different lifetimes). Access: members
!   read/write through the facade-internal conveyance channel (same vehicle as ws,
!   wrapper shape C2); force_eval (the adiab-mode delivery-request source) and recorder
!   (the n_hop column) read through the public read-only getter pair occ_get / n_hop_get.
!   The degenerate occupation occ=[1] is seeded at method_bind — it arrives with the
!   method declaration, not with the container (the adiabatic degenerate point is method
!   knowledge); dimension n_e set at member init, reallocated together with a.
!   (supersedes spec §4's "occ and n_hop relocate to state_t top level" clause)
!
! -- E17 container revision (review-gate rulings 2026-09-12; force_facade's E17 header
!    block is the primary text - this note carries the method-row share only) --
!   The method row gains h_form - the method's internal Hamiltonian PRESET (which H
!   forms it consumes: comma-list over the r_* mode names; empty = consumes no H
!   constructs, the degenerate method - 'adiabatic' registers ''). This is E14's
!   deferred deliverable-declaration field activated early: container_init matches the
!   bound row's h_form against the container's h_deliver (the container_bind_pes
!   declaration) at assembly time - the E10 runtime repr-mismatch stop moves earlier,
!   it does not go away. Terminology: blackbox -> container (blackbox_init ->
!   container_init); historical blocks above keep their original wording.
!=====================================================================
module elec_facade
   use state,   only: state_t, elec_t
   use config,  only: electronic
   use elec_ws, only: ws_t, r_empty
   implicit none
   private
   public :: elec_prop, mqc_statistic            ! evolution-phase entries (spec §3, entries 2-3)
   public :: method_reg, method_bind             ! package registration + assembly lookup (spec §5)
   public :: method_h_form                       ! E17: the bound row's H-form declaration (the
                                                  ! assembly-time match source for container_init)
   public :: ws_get, ws_reset                    ! ws-holder accessor pair (spec §3 ws holding)
   public :: occ_get, n_hop_get                  ! E16 statistics-holder read-only getters
   public :: elec_prop_i, mqc_match_i, sync_proc_i  ! member contract interfaces (E14 T5: made
   ! explicit for method-package authors - the members/ files use them as contract
   ! documentation; conformance is machine-checked at the method_reg call site)

   ! Method-package member interfaces (spec §5 — one registration fills the whole package;
   ! the registering author knows the pairing is legal — no combination table, E13 carried).
   ! Layer semantics per design_spec §4.1
   abstract interface
      subroutine elec_prop_i(sta, dt)      ! layer 2: propagation member (electronic propagation)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! reads ws (through the facade-internal channel —
         ! see the conveyance note below) + sta%q; writes the method's declared ontic container
         ! (sta%s%a or sta%s%rho, spec §4 symmetric candidates; adiabatic/LZSH/MDEF allocate
         ! neither). Forms: amplitude ODE (TSH/IESH/Ehrenfest), mapping flow (SQC), analytic
         ! elimination (MDEF: no-op), no-op (LZSH, adiabatic)
         real(8), intent(in) :: dt             ! current step [10 fs] (the hop-probability channel
         ! g_ij(dt) — E9/E14: dt belongs to the electronic phases, not to the force evaluation;
         ! this layer consumes no rng — the stream belongs to layer 3)
      end subroutine elec_prop_i
      subroutine mqc_match_i(sta)          ! layer 3: matching member (state-statistics matching)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! reads the electronic state + ws statistics (d_ij
         ! for the momentum projection — the one-step ws read this entry exists for); writes
         ! the held occ / sta%p (occupation reassignment through the conveyance channel,
         ! E16 statistics holder + momentum projection — classical-side statistical
         ! constructs; a hop = one component change of the occupation vector),
         ! optionally sta%f (post-hop force refresh: the member may call force_eval itself —
         ! intra-box calls are legal), optionally the ontic container (decay-of-mixing amplitude
         ! resets). Sole consumption point of the random stream during the evolution phase
         ! (§4.1 ruling — one stream, guaranteed by rng)
      end subroutine mqc_match_i
      subroutine sync_proc_i(s)            ! container-synchronization member (spec §4 sync contract)
         import :: elec_t
         type(elec_t), intent(inout) :: s  ! brings the method's declared deliverable containers
         ! up to date (rho = the observable face of a, or a rebuilt from rho, per the method's
         ! own declaration); direction/form (a→rho, rho→a, none) is method knowledge — the
         ! facade calls it as a pure-form tail and never interprets it
      end subroutine sync_proc_i
   end interface

   ! ws conveyance to members (E10 axiom A intact under E14): the member-facing interfaces
   ! above take only (sta[, dt]) — ws is never a public-signature argument. Members read
   ! the held ws through the FACADE-INTERNAL channel: elec_prop / mqc_statistic fetch the
   ! module instance (the holder declaration below) and convey it into the member through
   ! an internal wrapper signature carrying ws_t; the wrapper layer (the exact internal
   ! signature the method members implement) is finalized in C2 together with the
   ! member-signature freeze. The skeleton legislates the channel's existence and privacy,
   ! not its shape (class-A pes members, by contrast, carry ws in their signature — the two
   ! routes are the "signature or facade-internal channel" pair of the elec_ws.f90 note).
   ! (see previous line)

   ! Method-package table (single-select by name at assembly; rows registered via reg_*.f90
   ! under methods/ — spec §5 two-level registration; read-only during evolution)
   integer, parameter :: max_method = 16  ! package-table capacity [row] (K2 seed: 'adiabatic';
                                          ! 'ehrenfest'/'fssh'/'iesh' land in phase C)
   type :: method_entry_t
      character(len=32) :: name = ''      ! package name (the config%electronic%method assembly key)
      character(len=32) :: h_form = ''    ! the method's Hamiltonian preset - which H forms it
                                          ! consumes (comma-list of r_* mode names; empty = no H
                                          ! constructs, the degenerate method; E17(5) - matched
                                          ! against the container's h_deliver at assembly)
      procedure(elec_prop_i), pointer, nopass :: prop  => null()  ! layer-2 member (propagation)
      procedure(mqc_match_i), pointer, nopass :: match => null()  ! layer-3 member (matching)
      procedure(sync_proc_i), pointer, nopass :: sync  => null()  ! container sync (null = the
                                          ! method declares no deliverable cache — legal, spec §4)
   end type method_entry_t
   type(method_entry_t) :: method_tbl(max_method)  ! method-package registry (module-private)
   integer :: n_method = 0                ! registered row count [row]
   integer :: i_assembled = 0             ! assembled row index [-] (0 = not bound; bound at
                                          ! container_init from config%electronic%method — the
                                          ! single assembly entry point stays in force_facade,
                                          ! spec §3 assembly-phase entries)

   ! The ws instance — E10 axiom C / E14: module-private holder (saved, not a dummy),
   ! one-integration-step lifetime, invisible outside the box (never a dummy of any public
   ! signature; the accessor pair below is the only sanctioned handle)
   type(ws_t) :: ws

   ! The statistics holder — E16: run-lifetime module-private (the ws sibling; see the
   ! header note). occ written by mqc_statistic's member (through the conveyance channel),
   ! read by force_eval / recorder through occ_get; n_hop incremented at hop detection,
   ! read through n_hop_get. Seeded occ=[1] / n_hop=0 at method_bind (the degenerate
   ! occupation arrives with the method declaration)
   integer, allocatable :: occ(:)   ! occupation vector (E14 generalization of i_surf:
                                    ! occ(k) = state of electron k / the configuration
                                    ! index — length 1 = the degenerate single-active-
                                    ! surface form; dimension n_e set at member init)
   integer :: n_hop = 0             ! hop counter (pure observable, recorder-only read)

contains
   !------------------------------------------------------------------
   ! ws_get(w) - hand the held ws instance to force_facade (the fill site of the one-step
   !             lifetime); members never call this — they receive ws through the
   !             facade-internal conveyance channel (see the note above)
   ! Reference: spec §3 "ws holding and the endpoint-force contract" (force_facade fills it:
   !       use elec_facade to obtain the instance); elec_ws.f90 axiom C comment (holder)
   subroutine ws_get(w)
      type(ws_t), intent(inout) :: w  ! the fill handle on the held instance (force_eval writes
                                      ! the step's composition through it — the request field
                                      ! i_want from the held occ (E16: same module, direct
                                      ! read — no getter round-trip), then the per-mode fields
                                      ! per the elec_ws.f90 container contract — and leaves it
                                      ! in place for the electronic phases)
      ! Pseudocode:
      ! 1. Hand the held instance out through the dummy (the only sanctioned handle on the
      !    holder — force_eval calls this at its head, composes into the container, and
      !    returns; the instance stays module-private, the caller never owns it)
      ! 2. The exact association mechanism belongs to the C2 member-signature freeze (the
      !    same window that fixes the member-side conveyance wrapper); the K2 statute is
      !    the contract: what force_eval writes here is what elec_prop / mqc_statistic next
      !    read, and nothing outside the box ever holds the instance
   end subroutine ws_get

   !------------------------------------------------------------------
   ! ws_reset() - reset the held ws to the empty state (repr = r_empty; each step end + the
   !              force_eval head re-request — axiom C one-step lifetime, E14 form)
   ! Reference: elec_ws.f90 (r_empty tag; axiom C comment: filled at force evaluation,
   !       consumed by elec_prop / mqc_statistic, reset to r_empty at step end); spec §3
   subroutine ws_reset()
      ! Pseudocode:
      ! 1. ws%repr = r_empty (the stale per-mode arrays become unreadable — repr gates every
      !    consumer read; the next force_eval re-requests from scratch)
      ! 2. Invoked at each step end (after mqc_statistic has consumed the statistics) and
      !    re-issued as the force_eval head reset (every evaluation begins from the empty
      !    container — together the two sites make the one-step lifetime robust to skipped
      !    phases); allocation/reuse of the per-mode arrays is finalized in C2 with the
      !    facade concretization (elec_ws.f90 tail note — the skeleton legislates the repr
      !    reset only)
   end subroutine ws_reset

   !------------------------------------------------------------------
   ! occ_get(occ) - read-only hand-out of the held occupation set (E16): force_eval's
   !                adiab-mode delivery-request source + the recorder columns
   ! Reference: this module's E16 header note (statistics holder); spec §4 occupation
   !       semantics (the location clause superseded — the holder is here)
   subroutine occ_get(occ)
      integer, allocatable, intent(out) :: occ ! read-only copy of the held set (the
                                               ! caller never receives the holder itself —
                                               ! no allocation handle: the single writer
                                               ! stays mqc_statistic's member)
      ! Pseudocode:
      ! 1. Unallocated holder (pre-assembly window — method_bind not yet run) ->
      !    loud-stop; otherwise copy the held occ into the dummy (allocate + assign)
   end subroutine occ_get

   !------------------------------------------------------------------
   ! n_hop_get() - read-only hand-out of the hop counter (E16): the recorder column
   ! Reference: this module's E16 header note; spec §4 (n_hop a pure observable)
   function n_hop_get() result(n)
      integer :: n                            ! the current hop count [count]
      ! Pseudocode:
      ! 1. n = the held n_hop
   end function n_hop_get

   !------------------------------------------------------------------
   ! method_reg(name, h_form, prop, match [, sync]) - register one method package into the table
   !                 (the whole row at once — the row is the atomic registration unit, spec §5)
   ! Reference: spec §5 (registration data flow: library level = each method file registered
   !       through its own thin reg_*.f90 in methods/; built-in set 'adiabatic' (all no-ops —
   !       the degenerate package) / 'ehrenfest' / 'fssh' / 'iesh'; new method = one new file
   !       + one thin reg file — drag in and it compiles); Corollary 1 (Ehrenfest and FSSH
   !       share all electronic-propagation code — the same prop row, different match);
   !       E17(5) (the h_form Hamiltonian-preset field)
   subroutine method_reg(name, h_form, prop, match, sync)
      character(len=*), intent(in) :: name  ! package name (the config%electronic%method key,
                                            ! e.g. 'adiabatic'/'ehrenfest'/'fssh'/'iesh')
      character(len=*), intent(in) :: h_form ! the method's Hamiltonian preset - comma-list of
                                            ! the r_* mode names it consumes (empty = no H
                                            ! constructs, the degenerate method; E17(5))
      procedure(elec_prop_i) :: prop        ! layer-2 member (propagation)
      procedure(mqc_match_i) :: match       ! layer-3 member (matching)
      procedure(sync_proc_i), optional :: sync ! container sync (absent = the method declares
                                            ! no deliverable cache)
      ! Pseudocode:
      ! 1. Duplicate name in method_tbl or full table (n_method = max_method) -> loud-stop
      !    (purely formal name-table hygiene; the E13 statute: pairing legality is NOT
      !    checked here — one registration vouches for the whole package, the registering
      !    author knows their pairing is legal)
      ! 2. Register the row {name, h_form, prop, match, sync} (sync absent -> null pointer,
      !    legal — spec §4: null = no deliverable container to keep up to date); executed by
      !    the reg_*.f90 files during the assembly phase (the driver calls the prebuilt
      !    registry_gen, which aggregates every reg file under methods/ and systems/ —
      !    zero registration lines in this file, E11)
   end subroutine method_reg

   !------------------------------------------------------------------
   ! method_h_form() - read-only hand-out of the bound row's Hamiltonian preset (E17(5)):
   !                   the match source consumed by container_init's assembly-time
   !                   h_form-vs-h_deliver check (force_facade)
   ! Reference: this module's E17 header note; force_facade's E17 header block (5)
   function method_h_form() result(form)
      character(len=32) :: form              ! the bound row's h_form (comma-list / empty)
      ! Pseudocode:
      ! 1. Unbound table (i_assembled = 0 - method_bind not yet run) -> loud-stop (a
      !    match check before assembly is an assembly-order error, formal); otherwise
      !    form = method_tbl(i_assembled)%h_form
   end function method_h_form

   !------------------------------------------------------------------
   ! method_bind(name) - assembly-time name lookup into method_tbl (single-select): bind
   !                 i_assembled — called from container_init, which stays in force_facade as
   !                 the single assembly entry point spanning both facades
   ! Reference: spec §3 (assembly-phase entries: container_init spans both facades — names
   !       resolve, method slot filled); §5 (system level: the method arrives as
   !       config%electronic%method = 'fssh' — assembly is a pure name-lookup into the
   !       method table, zero method registration lines for the system author)
   subroutine method_bind(name)
      character(len=*), intent(in) :: name ! package name (config%electronic%method — a name
                                           ! string, same form as pes_main, E13)
      ! Pseudocode:
      ! 1. Empty table (n_method = 0 — no library registered, registry_gen missing) or name
      !    not found -> loud-stop naming the requested method (a purely formal check: the
      !    name is data, its physical legality is not judged here, E13)
      ! 2. Already bound (i_assembled /= 0) -> loud-stop (container_init ran twice — an
      !    assembly-choreography error, formal)
      ! 3. Bind i_assembled to the matched row index (single-select; from here the
      !    evolution entries dispatch through this row — the table is read-only during
      !    evolution)
      ! 4. Seed the statistics holder (E16): occ = [1] (the degenerate occupation arrives
      !    with the method declaration — the adiabatic degenerate point is method
      !    knowledge, not container knowledge), n_hop = 0; larger occupation dimensions
      !    (n_e) reallocated by the member-init hook together with a
   end subroutine method_bind

   !------------------------------------------------------------------
   ! elec_prop(sta, dt) - entry 2 (spec §3; §4.1 layer 2): electronic propagation — table
   !                 dispatch to the assembled package's prop member + the pure-form sync tail
   ! Reference: spec §3 (entry 2: reads the held ws + sta%q, writes the propagation
   !       container; pure-form tail call sync_proc if registered), §4 (sync contract);
   !       §3 endpoint-force contract (prop_step returns only after the endpoint force has
   !       been evaluated — the held ws corresponds to the returned sta%q; under MQC the
   !       velocity-Verlet right-end force cannot be folded into the next step's left end:
   !       the amplitude propagation of this phase needs the endpoint H);
   !       docs/design_spec_completeness.md §4.1 (layer-2 semantics)
   subroutine elec_prop(sta, dt)
      type(state_t), intent(inout) :: sta  ! physical state (sta%q read by the member; the
                                           ! method's declared ontic container sta%s%a /
                                           ! sta%s%rho written — spec §4 symmetric candidates)
      real(8), intent(in) :: dt            ! current step [10 fs] (the g_ij(dt) channel — E9
                                           ! re-anchored onto this entry by E14; the value is
                                           ! the propagator's actual step, not a config copy)
      ! Pseudocode:
      ! 1. Unassembled table (i_assembled = 0 — method_bind not yet run) -> loud-stop (an
      !    evolution-phase entry called before container_init — assembly-order error, formal)
      ! 2. Convey the held ws into the member through the facade-internal channel (fetch
      !    the module instance; the internal wrapper signature carrying ws_t — the
      !    member-facing interface takes only (sta, dt), E10 axiom A intact; wrapper shape
      !    finalized in C2 with the member-signature freeze)
      ! 3. Call the assembled row's prop member: electronic propagation given Q(t) (§4.1 —
      !    amplitude ODE along the trajectory / mapping flow / analytic elimination / no-op;
      !    the member reads ws interpreted per ws%repr — an unsupported repr loud-stops
      !    inside the member, E10/E13 contract — and writes the declared ontic container;
      !    dt semantics the member cannot live with loud-stop there too, E13 open point)
      !    (see previous line)
      ! 4. Sync tail (pure-form): if the row's sync pointer is attached, call sync(sta%s) —
      !    bring the declared deliverable containers up to date (a→rho / rho→a as the
      !    method declared; the facade never interprets the direction, spec §4)
   end subroutine elec_prop

   !------------------------------------------------------------------
   ! mqc_statistic(sta) - entry 3 (spec §3; §4.1 layer 3): state-statistics matching — table
   !                 dispatch to the assembled package's match member + the pure-form sync tail
   ! Reference: spec §3 (entry 3: input the electronic state + ws statistics, output the
   !       state — occ/p are classical-side statistical fields, spec §4; pure-form tail call
   !       sync_proc if registered), §4 (sync tail), §5 level 1 (the dispatch chain: facade
   !       table lookup -> match member — reads sta%s%a + sta%p + ws, g_ij from the library
   !       formula, rng-sampled hop, writes the held occ/sta%p (E16), optional force_eval self-call ->
   !       sync tail);
   !       docs/design_spec_completeness.md §4.1 (layer-3 semantics + the rng ruling)
   subroutine mqc_statistic(sta)
      type(state_t), intent(inout) :: sta  ! physical state (sta%p written; sta%f
                                           ! and the ontic container optionally — the
                                           ! member's declaration, see step 3; the held occ
                                           ! (E16 statistics holder, this module) feeds the
                                           ! next force_eval's adiab-mode delivery request)
      ! Pseudocode:
      ! 1. Unassembled table (i_assembled = 0) -> loud-stop (same formal check as elec_prop)
      ! 2. Convey the held ws into the member through the facade-internal channel (the
      !    one-step ws read this entry exists for: d_ij survives the force evaluation —
      !    under the old facade-local ws an independent matching entry had nothing to read)
      ! 3. Call the assembled row's match member: read the electronic statistics
      !    (populations / gaps / d_ij from ws + the ontic container), sample, and write the
      !    classical side — the held occ (E16: the statistics holder through the
      !    conveyance channel; a hop = one component change, n_hop incremented at
      !    detection), sta%p (momentum projection p' = p - beta*d_ij), optionally sta%f
      !    (post-hop force refresh: the member may call force_eval itself — intra-box,
      !    legal) and optionally the ontic container (decay-of-mixing amplitude resets).
      !    Sole consumption point of the random stream during the evolution phase (§4.1
      !    ruling — hop_decide is the rng leaf, spec §5 level 2)
      !    (see previous line)
      ! 4. Sync tail (pure-form): if the row's sync pointer is attached, call sync(sta%s) —
      !    same contract as elec_prop step 4 (both entries carry the tail: post-hop
      !    container resets and post-propagation refreshes are indistinguishable to the
      !    facade)
   end subroutine mqc_statistic
end module elec_facade
