!=====================================================================
! force_facade.f90 - L2 interface-A facade (E14 four-entry form, E17 container form): entry 1
!                    - force_eval, (Q,P,s)->(F,E) (the unit-contract fulfillment point -
!                    x23.0605 x e_conv appears only here) + the container single-slot
!                    binding + container_init, the single assembly entry point spanning
!                    BOTH facades
! Reference: docs/final_completion_plan.md §5 N1 (registration/assembly spec - directly in force; its
! negotiation-table clause superseded by E13, see the E13 note below the slot)
! Member unit contract: class-A members deliver force kernels/energies in eV (eV is the exchange unit);
! eV -> internal-unit conversion is centralized in this facade (the unit-contract grep
! requirement that the eV-to-kcal/mol literal appear only in this file - hence the "F_add code
! units" note in the member-signature spec reads as "after facade conversion")
! Member tables are filled when the interface phase lands the real containers (main PES + additive
! stochastic-law members; E12 ruling: pes members are system material - each systems/<system>/
! folder carries self-contained functional form + parameters; the facade itself is untouched)
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
!   fixed with the first real dt-sensitive members). Rationale: blackbox authors are assumed to know their system,
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
!
! -- E18 termination ownership (review-gate ruling 2026-09-12: "predicate is overcomplex -
!    let the container write its own termination condition; the output being whether the
!    trajectory terminates is enough") --
!   The predicate grammar interpreter is WITHDRAWN (pred_parse expression trees / the four
!   built-in primitives / PRIM open registration / the boolean closed combinators / the k_*
!   node codes / frag_ref addressing - the whole file deleted; a small language built for
!   one judgment is overengineering, the E17 singleton argument applied to the problem
!   side). Termination is problem knowledge, equally system-specific (E17(4)'s lone non-H
!   exception): the container registers ONE judgment function - container_bind_term(proc),
!   the private contract logical term_i(sta) - and the driver phase-4 call becomes
!   container_term(sta) (unbound slot = legal: the safety-net-only run). MAX_STEPS stays
!   the ever-present safety net as control%max_steps - program safety, generic side, the
!   only judgment never delegated. Geometry criteria (COM distance / recession velocity /
!   height) are written container-side on utils/geometry instruments (cenmas etc.).
!   The outcome expression family (output/outcome.f90) is the same complexity question -
!   flagged for its L4 review.
!=====================================================================
module force_facade
   use consts,        only: e_conv
   use state,         only: state_t
   use config,        only: electronic
   use config_roster, only: roster
   use elec_ws,       only: ws_t, r_scalar, r_diab, r_adiab, r_empty
   use elec_facade,   only: ws_handle, ws_reset, method_bind, occ_get, method_h_form
   use linalg,        only: eig_sym
   implicit none
   private
   public :: force_eval, container_bind_pes, container_init
   public :: container_bind_term, container_term  ! E18: the termination slot + its phase-4
                                                  ! evaluation entry

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

   ! The termination slot - E18: the container's own trajectory-termination judgment
   ! (problem knowledge; the predicate grammar interpreter withdrawn). Contract: read-only
   ! sta; geometry criteria written container-side with utils/geometry instruments.
   ! Unbound slot = legal (the safety-net-only run: fixed-step adiabatic cases; the
   ! ever-present net is control%max_steps, checked in the driver loop condition)
   abstract interface
      logical function term_i(sta)
         import :: state_t
         type(state_t), intent(in) :: sta  ! physical state (read-only snapshot)
      end function term_i
   end interface
   procedure(term_i), pointer :: container_term_proc => null()

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
      if (associated(container_pes)) then
         call stop_ff('container_bind_pes', 'the class-A slot is already bound ('// &
                      'container_bind_pes ran twice - an assembly-choreography error, formal)')
      end if
      container_pes => proc
      if (present(h_deliver)) then
         call store_deliver(h_deliver)
      else
         call store_deliver('')                  ! default: scalar only, the degenerate container
      end if
   end subroutine container_bind_pes

   !------------------------------------------------------------------
   ! container_bind_term(proc) - E18 single-slot assembly injection: bind the container's
   !                 trajectory-termination judgment (the caller: the container's reg
   !                 module, alongside container_bind_pes; no grammar, no table)
   ! Reference: this file's E18 header note; E17(4) (termination = problem knowledge, the
   !       lone non-H container exception)
   subroutine container_bind_term(proc)
      procedure(term_i) :: proc              ! the container termination judgment
      if (associated(container_term_proc)) then
         call stop_ff('container_bind_term', 'the termination slot is already bound ('// &
                      'container_bind_term ran twice - an assembly-choreography error, formal)')
      end if
      container_term_proc => proc            ! executed once in the assembly phase by the
                                             ! container's reg module
   end subroutine container_bind_term

   !------------------------------------------------------------------
   ! container_term(sta) - phase-4 evolution-loop entry (E18): the trajectory-termination
   !                 boolean (the driver's fourth phase; supersedes pred_eval in the
   !                 E7/E14 four-phase call list)
   ! Reference: this file's E18 header note; control%max_steps (the ever-present safety
   !       net - the driver loop's step_i < max_steps condition, the only judgment never
   !       delegated)
   logical function container_term(sta)
      type(state_t), intent(in) :: sta       ! physical state (read-only snapshot)
      ! unbound slot = the safety-net-only run (legal): .false. - the driver loop's
      ! step_i < max_steps condition remains the ever-present net; a bound slot runs
      ! the container's own judgment
      container_term = .false.
      if (.not. associated(container_term_proc)) return
      container_term = container_term_proc(sta)
   end function container_term

   !------------------------------------------------------------------
   ! container_init() - assembly + purely formal validation spanning BOTH facades (E13
   !                 form carried through E14/E17; renamed facade_init -> blackbox_init
   !                 (E14) -> container_init (E17); the single assembly entry point of
   !                 the whole container; composition legality is the container's own
   !                 affair - it already executed inside the bound aggregate entry)
   ! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §3 (assembly-phase
   !       entries - the E17 form: slot bound, H forms matched, method slot filled,
   !       observable-vs-container delivery contradiction); docs/design_spec_completeness.md
   !       §3.2-2 (member-init hooks are container-internal, executed by the container's own
   !       reg module around its binding)
   subroutine container_init()
      character(len=32) :: hf                 ! the bound method row's Hamiltonian preset
      ! 1. Method-slot assembly (E14) - runs in BOTH postures: the idle chain (no
      !    container dragged in, E20) still binds the degenerate no-op method row, so
      !    the body phases have something to run against. A pure name lookup - the
      !    empty-table / not-found / already-bound loud-stops live in method_bind
      !    itself (the table lives in elec_facade; the init choreography stays here
      !    as the single assembly entry point)
      call method_bind(electronic%method)
      ! 2. Idle-chain posture (E20): an unbound class-A slot means no container was
      !    dragged in (an empty container folder / a missing reg module) - the
      !    assembly proceeds, neither binding nor crashing; the loud-stop authority
      !    for the missing slot is force_eval's FIRST CALL (an evolution-phase entry
      !    called before any binding names the slot there). The checks below read the
      !    container's own declarations and run only when one is bound
      if (.not. associated(container_pes)) return
      ! 3. Container-contract formal check: a multi-surface system needs a multi-surface
      !    delivery mode (E13 kept: derivable from the container contract alone)
      if (electronic%n_surf >= 2) then
         if (.not. names_mode(h_deliver, 'diab') .and. .not. names_mode(h_deliver, 'adiab')) then
            call stop_ff('container_init', 'multi-surface system (n_surf >= 2) but the '// &
                         'declared delivery set names no multi-surface mode: "'//trim(h_deliver)//'"')
         end if
      end if
      ! 4. H-form match (E17(5)): the method row's demand must be a subset of the
      !    container's supply; an empty h_form (the degenerate method) always passes
      hf = method_h_form()
      if (.not. subset_of(hf, h_deliver)) then
         call stop_ff('container_init', 'the method Hamiltonian preset "'//trim(hf)// &
                      '" is not covered by the container delivery set "'//trim(h_deliver)//'"')
      end if
      ! 5. Container-side member init: executed by the container's reg module around its
      !    container_bind_pes call (member-init hooks are container-internal; parameters
      !    arrive through the members' own params, D3) - not choreographed here
      ! 6. Observable-vs-container formal check (spec §4): deferred to the interface phase
      !    with the observable declarations; integrator x electron-law legality is NOT
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
   ! Reference: spec §3 entry-1 contract (as a side effect the ws of the current geometry is
   !       left in place for the later electronic phases); the module-header unit contract
   !       (this routine is the unique conversion point of the whole library)
   subroutine force_eval(sta, e_tot)
      type(state_t), intent(inout) :: sta  ! physical state (q/p/s read in; f written out - the composed force)
      real(8), intent(out) :: e_tot        ! total energy H = T(P)+V [kcal/mol] (N1 convention;
                                           ! multi-surface mode V = current surface / law-member convention, see 4')
      ! Workspace (E10 axioms B/C as amended by E14; holder association frozen as the
      ! pointer-handle form): the elec_facade-held instance, aliased through the pointer
      ! function ws_handle - the ONE sanctioned alias (an ordinary dummy would detach
      ! the fill from the holder by copy semantics). The fill site composes the step's
      ! representation constructs INTO THE HOLDER: the request field is pre-written
      ! through the alias, the container call receives the holder through the signature
      ! dummy, and the constructs stay in place for elec_prop / mqc_statistic of the
      ! same integration step (one-step lifetime; still invisible outside the box -
      ! never an argument of any public signature, zero footprint in the force_eval
      ! signature)
      type(ws_t), pointer :: w             ! the holder alias (ws_handle)
      integer, allocatable :: occ_req(:)   ! the delivery request - the held occupation set
      integer :: k, r, i, j, n
      real(8) :: g(size(sta%q))            ! folded force kernel [eV/A]
      real(8) :: v_sum, v_pot, t_kin       ! energy channels [eV / kcal/mol]
      real(8), allocatable :: hh(:,:), hev(:)  ! the diagonalization pair: the H copy
                                              ! (destroyed by eig_sym - it returns the
                                              ! eigenvectors BY COLUMN) + the eigenvalues
      ! 1. Zero the force-kernel accumulators; ws_reset() (the head reset - every
      !    evaluation begins from the empty container; with elec_facade's step-end reset
      !    the two sites together make the one-step lifetime robust to skipped phases);
      !    dt does not appear here (the E9 channel withdrawn - the hop probability
      !    g_ij(dt) is consumed in elec_prop)
      call ws_reset()
      g = 0.0d0
      v_sum = 0.0d0
      ! 2. Unbound slot -> loud-stop (an evolution-phase entry called before assembly -
      !    assembly-order error, formal)
      if (.not. associated(container_pes)) then
         call stop_ff('force_eval', 'the class-A slot is not bound (container_bind_pes not '// &
                      'yet run - an evolution-phase entry called before assembly)')
      end if
      ! 3. ONE call into the container (E17: the member sum - main PES + additive
      !    stochastic laws - executes INSIDE it), THROUGH THE HOLDER ALIAS: the request
      !    field i_want(:) = occ_get() (the delivery request, E16: generalized from the
      !    scalar i_surf; length 1 = the single-active-surface degenerate form; the
      !    occupation set stays out of the class-A signature and is conveyed through
      !    the request field, E10) is pre-written into the holder, and the container
      !    receives the holder itself, filling it per its own declared mode
      call occ_get(occ_req)
      w => ws_handle()
      allocate (w%i_want(size(occ_req)))
      w%i_want = occ_req
      call container_pes(sta%q, roster%mass, w)   ! the container fills ws per its own mode
      ! 4. Fold the delivered channel, branching on the container's SELF-DECLARED fill
      !    mode ws%repr (representation is data, E10; the assembly-time h_form match of
      !    container_init moved the mode mismatch earlier - the runtime tag remains the
      !    authority): every not-delivered or not-yet-consumed mode is a NAMED loud-stop,
      !    never an undefined read
      select case (w%repr)
      case (r_scalar)
         ! scalar channel: g(:) = ws%g_add, v_sum = ws%v_add (the container already
         ! summed - the facade folds ONCE; no member loop remains here)
         if (.not. allocated(w%g_add)) then
            call stop_ff('force_eval', 'scalar-mode delivery carries no g_add kernel '// &
                         '(repr = r_scalar but the array was not delivered - a container-contract violation)')
         end if
         if (size(w%g_add) /= size(sta%q)) then
            call stop_ff('force_eval', 'scalar-mode g_add is not of ndof length '// &
                         '(a container-contract violation, formal shape check)')
         end if
         g = w%g_add
         v_sum = w%v_add
      case (r_diab)
         ! 4'. Multi-surface diabatic branch: diagonalize the delivered H_ij on the
         !     spot (the adiabatic view: eigenvalues = per-surface V, eigenvectors =
         !     the representation transform U) and fold the OCCUPATION ROWS of dh_ij:
         !     d(lambda_k)/dq = sum over ordered pairs (i,j) of U(i,k) U(j,k)
         !     (dH/dq)_ij - the occupied eigenvalue's gradient (a symmetric dh_ij
         !     delivery counts every off-diagonal twice, the correct 2 U_i U_j d_ij;
         !     eigenvector sign ambiguity cancels in the U U products). The
         !     diagonalization is system-independent generic physics and lawfully
         !     stays in the facade; the dh_ij.U product for the hop probability is
         !     fetched by the electronic phases themselves (the facade does not
         !     pre-build it). v_pot convention: the current-occupation surfaces (the
         !     population-weighted variant is the class-B member's own convention,
         !     fixed by the consumer). The additive scalar kernel folds ON TOP of
         !     the construct-channel kernel (electronic-state-independent members
         !     feed the nuclear-force channel only) and its energy in as usual
         if (.not. allocated(w%h_ij) .or. .not. allocated(w%dh_ij)) then
            call stop_ff('force_eval', 'diab-mode delivery carries no h_ij/dh_ij '// &
                         'constructs (repr = r_diab but the arrays were not delivered - a container-contract violation)')
         end if
         n = size(w%h_ij, 1)
         if (size(w%h_ij, 2) /= n .or. size(w%dh_ij, 1) /= size(sta%q) .or. &
             size(w%dh_ij, 2) /= n .or. size(w%dh_ij, 3) /= n) then
            call stop_ff('force_eval', 'diab-mode construct shapes disagree (h_ij not '// &
                         'square or dh_ij not ndof x n_surf^2 - formal shape checks)')
         end if
         do r = 1, size(occ_req)
            if (occ_req(r) < 1 .or. occ_req(r) > n) then
               call stop_ff('force_eval', 'the held occupation set names a surface '// &
                            'outside the delivered H (formal index check)')
            end if
         end do
         allocate (hh(n, n), hev(n))
         hh = w%h_ij          ! a COPY: eig_sym destroys its matrix argument, and the
                              ! delivered constructs must SURVIVE in the holder for
                              ! the electronic phases of this step
         call eig_sym(hh, hev)      ! hev = per-surface V ascending; hh = eigenvectors by column
         do r = 1, size(occ_req)
            k = occ_req(r)
            v_sum = v_sum + hev(k)
            do j = 1, n
               do i = 1, n
                  g = g + hh(i, k)*hh(j, k)*w%dh_ij(:, i, j)
               end do
            end do
         end do
         if (allocated(w%g_add)) then         ! the additive kernel is OPTIONAL in
            if (size(w%g_add) /= size(sta%q)) then    ! multi-surface mode (a pure-
               call stop_ff('force_eval', 'the additive g_add is not of ndof length '// &
                            '(a container-contract violation, formal shape check)')
            end if
            g = g + w%g_add                    ! construct container leaves it absent)
         end if
         v_sum = v_sum + w%v_add
         deallocate (hh, hev)
      case (r_adiab)
         ! the adiabatic construct channel (v_surf/g_surf/d_ij) is not consumed by
         ! this entry yet - the adiabatic aggregate fold is unrealized; a named stop,
         ! never an undefined read of the scalar channel
         call stop_ff('force_eval', 'the adiabatic construct channel (v_surf/g_surf/d_ij) '// &
                      'is not consumed by this entry (the adiabatic aggregate fold is unrealized)')
      case default
         call stop_ff('force_eval', 'the container left the workspace empty (repr still '// &
                      'r_empty after the call - no delivery)')
      end select
      ! 5. Unit-contract fulfillment point (unique in the whole library): after
      !    pre-zeroing sta%f, element-wise sta%f(i) = sta%f(i) - g(i)*ev_kcal*e_conv
      !    (the one eV->internal conversion; how g is composed in multi-surface mode: 4')
      sta%f = 0.0d0
      sta%f = sta%f - g*ev_kcal*e_conv
      ! 6. NO electron-law consumption here (E14: withdrawn to elec_prop / mqc_statistic
      !    in elec_facade) - the composed ws of the current geometry is what the later
      !    phases next read (the spec §3 side-effect contract; holder association: the
      !    note at the declaration of w)
      ! 7. Energy [kcal/mol]: v_pot = v_sum*ev_kcal; t_kin = sum p^2/(2*roster%mass)/e_conv
      !    (per-atom fold); e_tot = t_kin + v_pot (H = T + V)
      v_pot = v_sum*ev_kcal
      t_kin = 0.0d0
      do k = 1, size(roster%mass)
         t_kin = t_kin + dot_product(sta%p(3*k-2:3*k), sta%p(3*k-2:3*k))/(2.0d0*roster%mass(k))
      end do
      t_kin = t_kin/e_conv
      e_tot = t_kin + v_pot
   end subroutine force_eval

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! store_deliver(s) - assign the delivery-set field (the public dummy name shadows the
   ! module field inside container_bind_pes; this helper reaches it through host association)
   subroutine store_deliver(s)
      character(len=*), intent(in) :: s     ! the declared delivery set (comma-list / '')
      h_deliver = s
   end subroutine store_deliver

   ! names_mode(list, mode) - does the comma-list name this single mode?
   logical function names_mode(list, mode)
      character(len=*), intent(in) :: list, mode
      names_mode = subset_of(mode, list)
   end function names_mode

   ! subset_of(demand, supply) - is every comma-list token of demand a token of supply?
   ! (an empty demand is the degenerate case - always covered)
   logical function subset_of(demand, supply)
      character(len=*), intent(in) :: demand, supply
      character(len=32) :: dtok, stok
      integer :: dp, sp
      logical :: found
      subset_of = .true.
      if (len_trim(demand) == 0) return
      dp = 1
      do while (next_token(demand, dp, dtok))
         found = .false.
         sp = 1
         do while (next_token(supply, sp, stok))
            if (trim(stok) == trim(dtok)) then
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            subset_of = .false.
            return
         end if
      end do
   end function subset_of

   ! next_token(s, pos, tok) - comma-list tokenizer: fetch the next non-empty token at or
   ! after pos, advancing pos past it; .false. when the string is exhausted
   logical function next_token(s, pos, tok)
      character(len=*), intent(in) :: s
      integer, intent(inout) :: pos
      character(len=*), intent(out) :: tok
      integer :: c
      next_token = .false.
      do
         if (pos > len(s)) return
         c = index(s(pos:), ',')
         if (c == 0) then
            tok = adjustl(s(pos:))
            pos = len(s) + 1
         else
            tok = adjustl(s(pos:pos + c - 1))
            pos = pos + c
         end if
         if (len_trim(tok) > 0) exit               ! skip empty items
      end do
      next_token = .true.
   end function next_token

   ! stop_ff(who, msg) - the loud-stop channel (name the caller, print the named
   ! message, STOP 1)
   subroutine stop_ff(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': container error: '//trim(msg)
      write (0, '(a)') trim(who)//': loud-stop (the container is not usable)'
      stop 1
   end subroutine stop_ff
end module force_facade
