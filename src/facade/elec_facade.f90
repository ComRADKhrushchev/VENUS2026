!=====================================================================
! elec_facade.f90 - electronic-phase facade: elec_prop (electronic
!   propagation), mqc_statistic (state-statistics matching), the
!   method-package registry, and the ws / statistics holders
! Design:
!   The two evolution-phase entries dispatch through the assembled method
!   package {prop, match [, sync]} - a pure name lookup at assembly;
!   pairing legality is the registering author's knowledge. The ws
!   instance (one-step lifetime) and the run-lifetime statistics (occ,
!   n_hop) are module-private holders, never public-signature arguments;
!   members receive ws through the facade-internal conveyance channel.
!=====================================================================
module elec_facade
   use state,   only: state_t, elec_t
   use config,  only: electronic
   use elec_ws, only: ws_t, r_empty
   implicit none
   private
   public :: elec_prop, mqc_statistic            ! evolution-phase entries
   public :: method_reg, method_bind             ! package registration + assembly lookup
   public :: method_h_form                       ! the bound row's H-form declaration (the
                                                  ! assembly-time match source for container_init)
   public :: ws_get, ws_reset                    ! ws-holder accessor pair
   public :: ws_handle                           ! the pointer alias onto the held ws
                                                  ! (box-internal consumers only - see below)
   public :: occ_get, n_hop_get                  ! statistics-holder read-only getters
   public :: elec_prop_i, mqc_match_i, sync_proc_i  ! member contract interfaces (used by
   ! method-package authors as contract documentation; conformance is machine-checked at
   ! the method_reg call site)

   ! Method-package member interfaces (one registration fills the whole package; the
   ! registering author knows the pairing is legal - no combination table).
   abstract interface
      subroutine elec_prop_i(sta, dt)      ! electronic-propagation member
         import :: state_t
         type(state_t), intent(inout) :: sta  ! reads ws (facade-internal channel) + sta%q;
         ! writes the method's declared ontic container (sta%s%a or sta%s%rho;
         ! adiabatic/LZSH/MDEF allocate neither). Forms: amplitude ODE, mapping flow,
         ! analytic elimination, no-op
         real(8), intent(in) :: dt             ! current step [10 fs] (the g_ij(dt)
         ! hop-probability channel; this member consumes no rng - the stream
         ! belongs to the matching member)
      end subroutine elec_prop_i
      subroutine mqc_match_i(sta)          ! state-statistics matching member
         import :: state_t
         type(state_t), intent(inout) :: sta  ! reads the electronic state + ws statistics
         ! (d_ij for the momentum projection); writes the held occ / sta%p (a hop
         ! = one component change of the occupation vector), optionally sta%f
         ! (post-hop refresh: intra-box force_eval calls are legal), optionally
         ! the ontic container. Sole consumption point of the random stream
         ! during the evolution phase (one stream, guaranteed by rng)
      end subroutine mqc_match_i
      subroutine sync_proc_i(s)            ! container-synchronization member
         import :: elec_t
         type(elec_t), intent(inout) :: s  ! brings the method's declared deliverable
         ! containers up to date (a->rho / rho->a / none, method knowledge);
         ! the facade calls it as a pure-form tail and never interprets it
      end subroutine sync_proc_i
   end interface

   ! ws conveyance to members: the member-facing interfaces take only (sta[, dt]) -
   ! ws is never a public-signature argument. The planned facade-internal wrapper
   ! carrying ws_t lands with the first stateful members; today the dispatch is
   ! real but the conveyance is vacuous (the degenerate member reads nothing), so
   ! the entries call the members through the frozen (sta[, dt]) interfaces directly

   ! Method-package table (single-select by name at assembly; read-only during evolution)
   integer, parameter :: max_method = 16  ! package-table capacity [row]
   type :: method_entry_t
      character(len=32) :: name = ''      ! package name (the config%electronic%method assembly key)
      character(len=32) :: h_form = ''    ! the method's Hamiltonian preset (comma-list of
                                          ! r_* mode names consumed; empty = none; matched
                                          ! against the container's h_deliver at assembly)
      procedure(elec_prop_i), pointer, nopass :: prop  => null()  ! propagation member
      procedure(mqc_match_i), pointer, nopass :: match => null()  ! matching member
      procedure(sync_proc_i), pointer, nopass :: sync  => null()  ! container sync (null =
                                          ! no deliverable cache, legal)
   end type method_entry_t
   type(method_entry_t) :: method_tbl(max_method)  ! method-package registry (module-private)
   integer :: n_method = 0                ! registered row count [row]
   integer :: i_assembled = 0             ! assembled row index [-] (0 = not bound; bound at
                                          ! container_init from config%electronic%method)

   ! The ws instance - module-private holder (target + save, not a dummy),
   ! one-integration-step lifetime, invisible outside the box (never a dummy of
   ! the public signatures). The holder carries TARGET and the pointer function
   ! ws_handle below hands out the one sanctioned alias through which the fill
   ! site writes the holder itself; ws_get / ws_reset ride on the same holder
   type(ws_t), target, save :: ws

   ! The statistics holder - run-lifetime module-private (the ws sibling):
   ! occ written by mqc_statistic's member (through the conveyance channel),
   ! read by force_eval / recorder through occ_get; n_hop incremented at hop
   ! detection, read through n_hop_get. Seeded occ=[1] / n_hop=0 at method_bind
   integer, allocatable :: occ_held(:)   ! occupation vector (occ(k) = the configuration
                                        ! index; length 1 = the degenerate single-active-
                                        ! surface form; dimension n_e set at member init)
   integer :: n_hop = 0             ! hop counter (pure observable, recorder-only read)

contains
   !------------------------------------------------------------------
   ! ws_handle() - the pointer alias onto the held ws instance (the fill
   !               site composes into the holder through it; the electronic
   !               phases read the same instance inside the box)
   !------------------------------------------------------------------
   function ws_handle() result(h)
      type(ws_t), pointer :: h            ! the alias onto the module-private holder
                                           ! (never disassociated - static storage)
      ! 1. Point at the holder: the one sanctioned alias (an ordinary dummy cannot
      !    alias a module-private variable - copy semantics would detach the fill
      !    from the held instance). Consumers: the two facades of the box only -
      !    public only because Fortran privacy ends at the module boundary; the
      !    driver / propagator / output never call it
      h => ws
   end function ws_handle

   !------------------------------------------------------------------
   ! ws_get(w) - copy-out convenience: a snapshot of the held ws (writers
   !             use ws_handle; members receive ws through the conveyance
   !             channel, never either accessor)
   !------------------------------------------------------------------
   subroutine ws_get(w)
      type(ws_t), intent(inout) :: w  ! a deep intrinsic-assignment copy of the held
                                      ! instance; writing w never touches the holder
      ! 1. Copy the holder out through the dummy (composing into the snapshot
      !    never reaches the holder - the one-step conveyance to the
      !    electronic phases rides the ws_handle alias)
      w = ws
   end subroutine ws_get

   !------------------------------------------------------------------
   ! ws_reset() - reset the held ws to the defined empty state (repr =
   !              r_empty, arrays deallocated, scalar zeroed); invoked at
   !              each step end + the force_eval head re-request
   !------------------------------------------------------------------
   subroutine ws_reset()
      ! 1. Full clear: each fill allocates afresh (no cross-step reuse); repr
      !    = r_empty gates every consumer read against stale tags. Invoked at
      !    each step end (after mqc_statistic has consumed the statistics) and
      !    re-issued as the force_eval head reset - the two sites together
      !    make the one-step lifetime robust to skipped phases
      if (allocated(ws%i_want)) deallocate (ws%i_want)
      if (allocated(ws%g_add))  deallocate (ws%g_add)
      if (allocated(ws%h_ij))   deallocate (ws%h_ij)
      if (allocated(ws%dh_ij))  deallocate (ws%dh_ij)
      if (allocated(ws%v_surf)) deallocate (ws%v_surf)
      if (allocated(ws%g_surf)) deallocate (ws%g_surf)
      if (allocated(ws%d_ij))   deallocate (ws%d_ij)
      ws%v_add = 0.0d0
      ws%repr = r_empty
   end subroutine ws_reset

   !------------------------------------------------------------------
   ! occ_get(occ) - read-only hand-out of the held occupation set
   !                (force_eval's delivery-request source + recorder columns)
   !------------------------------------------------------------------
   subroutine occ_get(occ)
      integer, allocatable, intent(out) :: occ(:) ! read-only copy of the held set (the
                                                  ! single writer stays mqc_statistic's
                                                  ! member)
      ! 1. Unallocated holder (pre-assembly window - method_bind not yet run) ->
      !    loud-stop; otherwise copy the held occ into the dummy (allocate + assign)
      if (.not. allocated(occ_held)) then
         call stop_ele('occ_get', 'the occupation holder is not seeded (method_bind not '// &
                       'yet run - an evolution-phase read before assembly)')
      end if
      allocate (occ(size(occ_held)))
      occ = occ_held
   end subroutine occ_get

   !------------------------------------------------------------------
   ! n_hop_get() - read-only hand-out of the hop counter (the recorder column)
   !------------------------------------------------------------------
   function n_hop_get() result(n)
      integer :: n                            ! the current hop count [count]
      ! 1. Hand out the held counter (a plain integer, always defined: seeded 0
      !    at method_bind; a pre-assembly read hands back the initialized 0)
      n = n_hop
   end function n_hop_get

   !------------------------------------------------------------------
   ! method_reg(name, h_form, prop, match [, sync]) - register one method
   !                 package row into the table (the atomic registration unit)
   !------------------------------------------------------------------
   subroutine method_reg(name, h_form, prop, match, sync)
      character(len=*), intent(in) :: name  ! package name (the config%electronic%method key)
      character(len=*), intent(in) :: h_form ! the method's Hamiltonian preset - comma-list of
                                            ! the r_* mode names it consumes (empty = none)
      procedure(elec_prop_i) :: prop        ! propagation member
      procedure(mqc_match_i) :: match       ! matching member
      procedure(sync_proc_i), optional :: sync ! container sync (absent = no deliverable cache)
      integer :: j
      ! 1. Duplicate name in method_tbl or full table (n_method = max_method) -> loud-stop
      !    (purely formal name-table hygiene; pairing legality is NOT checked here -
      !    one registration vouches for the whole package)
      do j = 1, n_method
         if (trim(method_tbl(j)%name) == trim(name)) then
            call stop_ele('method_reg', 'duplicate method name "'//trim(name)//'" (each '// &
                          'package registers exactly once, at assembly)')
         end if
      end do
      if (n_method >= max_method) then
         call stop_ele('method_reg', 'method table full (max_method rows exceeded - the '// &
                       'built-in set and the dragged-in packages exceeded the capacity)')
      end if
      ! 2. Register the row {name, h_form, prop, match, sync} (sync absent -> null
      !    pointer, legal: no deliverable container to keep up to date); executed
      !    by the reg files aggregated into registry_gen during the assembly phase
      n_method = n_method + 1
      method_tbl(n_method)%name = name
      method_tbl(n_method)%h_form = h_form
      method_tbl(n_method)%prop => prop
      method_tbl(n_method)%match => match
      if (present(sync)) then
         method_tbl(n_method)%sync => sync
      else
         method_tbl(n_method)%sync => null()
      end if
   end subroutine method_reg

   !------------------------------------------------------------------
   ! method_h_form() - the bound row's Hamiltonian preset (the match source
   !                   for container_init's assembly-time check)
   !------------------------------------------------------------------
   function method_h_form() result(form)
      character(len=32) :: form              ! the bound row's h_form (comma-list / empty)
      ! 1. Unbound table (i_assembled = 0 - method_bind not yet run) -> loud-stop
      !    (a match check before assembly is an assembly-order error, formal);
      !    otherwise form = method_tbl(i_assembled)%h_form
      if (i_assembled == 0) then
         call stop_ele('method_h_form', 'no method row bound (method_bind not yet run - '// &
                       'a match check before assembly, formal)')
      end if
      form = method_tbl(i_assembled)%h_form
   end function method_h_form

   !------------------------------------------------------------------
   ! method_bind(name) - assembly-time name lookup into method_tbl
   !                 (single-select): bind i_assembled (called from
   !                 container_init)
   !------------------------------------------------------------------
   subroutine method_bind(name)
      character(len=*), intent(in) :: name ! package name (config%electronic%method - a
                                           ! name string)
      integer :: j, i_row
      ! 1. Empty table (n_method = 0 - no library registered, registry_gen missing) or
      !    name not found -> loud-stop naming the requested method (a purely formal
      !    check: the name is data, its physical legality is not judged here)
      if (n_method == 0) then
         call stop_ele('method_bind', 'the method table is empty (no package registered - '// &
                       'the aggregated registry module is missing or registered nothing)')
      end if
      i_row = 0
      do j = 1, n_method
         if (trim(method_tbl(j)%name) == trim(name)) i_row = j
      end do
      if (i_row == 0) then
         call stop_ele('method_bind', 'the requested method "'//trim(name)//'" names no '// &
                       'registered row (a purely formal name lookup, E13)')
      end if
      ! 2. Already bound (i_assembled /= 0) -> loud-stop (container_init ran twice - an
      !    assembly-choreography error, formal)
      if (i_assembled /= 0) then
         call stop_ele('method_bind', 'a method row is already bound (container_init ran '// &
                       'twice - an assembly-choreography error, formal)')
      end if
      ! 3. Bind i_assembled to the matched row index (single-select; from here the
      !    evolution entries dispatch through this row — the table is read-only during
      !    evolution)
      i_assembled = i_row
      ! 4. Seed the statistics holder: occ = [1] (the degenerate occupation arrives
      !    with the method declaration - the adiabatic degenerate point is method
      !    knowledge, not container knowledge), n_hop = 0; larger occupation
      !    dimensions (n_e) reallocated by the member-init hook together with a
      if (allocated(occ_held)) deallocate (occ_held)
      allocate (occ_held(1))
      occ_held(1) = 1
      n_hop = 0
   end subroutine method_bind

   !------------------------------------------------------------------
   ! elec_prop(sta, dt) - electronic propagation: dispatch to the
   !                 assembled package's prop member + the sync tail
   !------------------------------------------------------------------
   subroutine elec_prop(sta, dt)
      type(state_t), intent(inout) :: sta  ! physical state (sta%q read by the member; the
                                           ! method's declared ontic container sta%s%a /
                                           ! sta%s%rho written)
      real(8), intent(in) :: dt            ! current step [10 fs] (the g_ij(dt) channel;
                                           ! the propagator's actual step, not a config copy)
      ! 1. Assembly-order guard: an unassembled table (i_assembled = 0 - method_bind
      !    not yet run) loud-stops (an evolution-phase entry called before assembly)
      if (i_assembled == 0) then
         call stop_ele('elec_prop', 'no method row bound (method_bind not yet run - '// &
                       'an evolution-phase entry called before assembly, formal)')
      end if
      ! 2. ws conveyance: vacuous in the current posture - the frozen member
      !    interfaces take only (sta, dt) and the one registered member (the
      !    degenerate package) reads nothing, so no wrapper call happens; the
      !    holder already carries this step's constructs (filled by force_eval -
      !    prop_step returns only after the endpoint force, so the held ws
      !    corresponds to sta%q). The wrapper signature carrying ws_t lands
      !    together with the first stateful members
      ! 3. Dispatch: call the assembled row's prop member - electronic propagation
      !    given Q(t) (amplitude ODE along the trajectory / mapping flow /
      !    analytic elimination / no-op). The member reads ws interpreted per
      !    ws%repr - an unsupported repr or dt semantics the member cannot live
      !    with loud-stops INSIDE the member (the degenerate member ignores dt)
      call method_tbl(i_assembled)%prop(sta, dt)
      ! 4. Sync tail: if the row's sync pointer is attached, call sync(sta%s) -
      !    bring the declared deliverable containers up to date (a->rho /
      !    rho->a as the method declared; the facade never interprets the
      !    direction; null pointer = no deliverable cache, legal)
      if (associated(method_tbl(i_assembled)%sync)) then
         call method_tbl(i_assembled)%sync(sta%s)
      end if
   end subroutine elec_prop

   !------------------------------------------------------------------
   ! mqc_statistic(sta) - state-statistics matching: dispatch to the
   !                 assembled package's match member + the sync tail
   !------------------------------------------------------------------
   subroutine mqc_statistic(sta)
      type(state_t), intent(inout) :: sta  ! physical state (sta%p written; sta%f and
                                           ! the ontic container optionally - the
                                           ! member's declaration, see step 3; the
                                           ! held occ feeds the next force_eval's
                                           ! delivery request)
      ! 1. Assembly-order guard: same formal check as elec_prop (an evolution-phase
      !    entry called before assembly)
      if (i_assembled == 0) then
         call stop_ele('mqc_statistic', 'no method row bound (method_bind not yet run - '// &
                       'an evolution-phase entry called before assembly, formal)')
      end if
      ! 2. ws conveyance: vacuous in the current posture - same note as
      !    elec_prop step 2 (the one-step ws read this entry exists for: d_ij
      !    survives the force evaluation; the wrapper carrying ws_t lands with
      !    the first stateful members)
      ! 3. Dispatch: call the assembled row's match member - the member owns the
      !    whole statistics write-back. It reads the electronic statistics
      !    (populations / gaps / d_ij from ws + the ontic container), samples
      !    (sole consumption point of the random stream during the evolution
      !    phase) and writes the classical side: the held occ (a hop = one
      !    component change of the occupation vector, n_hop incremented at
      !    detection; run-lifetime channels, never touched by the step-end
      !    reset below), sta%p (momentum projection), optionally sta%f (the
      !    post-hop force-refresh hook: the member may call force_eval itself -
      !    intra-box, legal) and optionally the ontic container. The degenerate
      !    member writes nothing: occ stays [1], n_hop stays 0, no rng
      !    consumption
      call method_tbl(i_assembled)%match(sta)
      ! 4. Sync tail: if the row's sync pointer is attached, call sync(sta%s) -
      !    same contract as elec_prop step 4 (both entries carry the tail:
      !    post-hop container resets and post-propagation refreshes are
      !    indistinguishable to the facade)
      if (associated(method_tbl(i_assembled)%sync)) then
         call method_tbl(i_assembled)%sync(sta%s)
      end if
      ! 5. STEP-END RESET: ws_reset() after the match member has consumed the
      !    statistics - the closing call of the one-step lifetime (the site is
      !    this entry's tail, inside the box - the driver never calls ws_reset
      !    itself; force_eval's head re-request makes the lifetime robust to
      !    skipped phases)
      call ws_reset()
   end subroutine mqc_statistic

   !------------------------------------------------------------------
   ! private helper: the loud-stop channel of this facade (name the caller, print the
   ! named message, STOP 1)
   !------------------------------------------------------------------
   subroutine stop_ele(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': electronic-facade error: '//trim(msg)
      write (0, '(a)') trim(who)//': loud-stop (the method table is not usable)'
      stop 1
   end subroutine stop_ele
end module elec_facade
