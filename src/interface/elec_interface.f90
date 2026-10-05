!=====================================================================
! elec_interface.f90 - electronic-phase interface: elec_prop (electronic
!   propagation), mqc_statistic (state-statistics matching; MQC = mixed
!   quantum-classical), the
!   method-package registry, and the run-lifetime statistics holders
! Design:
!   The two evolution-phase entries dispatch through the assembled method
!   package {prop, match [, sync]} - a pure name lookup at assembly;
!   pairing legality is the registering author's knowledge. Method rows
!   register from methods/ (generic library packages) or from a systems/
!   folder (a container-composed package: the
!   H definition and a system-specific evolution law live in ONE folder
!   and share private modules - no interchange workspace).
!   Members reach system data through their own folder's modules or the
!   level-2 library components - never through this interface. The
!   statistics holders (occ, n_hop) are module-private, run-lifetime,
!   written by the match member through occ_set / n_hop_set and read by
!   force-side consumers / the recorder through the getters.
!=====================================================================
module elec_interface
   use state,   only: state_t, elec_t
   use config,  only: electronic
   implicit none
   private
   public :: elec_prop, mqc_statistic            ! evolution-phase entries
   public :: method_reg, method_bind             ! package registration + assembly lookup
   public :: occ_get, occ_set, n_hop_get, n_hop_set  ! statistics-holder access pair
                                                  ! (getters: force-side / recorder reads;
                                                  ! setters: the match member's write-back)
   public :: elec_prop_i, mqc_match_i, sync_proc_i  ! member contract interfaces (used by
   ! method-package authors as contract documentation; conformance is machine-checked at
   ! the method_reg call site)

   ! Method-package member interfaces (one registration fills the whole package; the
   ! registering author knows the pairing is legal - no combination table).
   abstract interface
      subroutine elec_prop_i(sta, dt)      ! electronic-propagation member
         import :: state_t
         type(state_t), intent(inout) :: sta  ! reads sta%q + the member's own data channels
         ! (its container folder's private modules / level-2 components); writes the
         ! method's declared ontic container (sta%s%a or sta%s%rho; adiabatic/LZSH/MDEF
         ! allocate neither). Forms: amplitude ODE, mapping flow, analytic elimination, no-op
         real(8), intent(in) :: dt             ! current step [10 fs] (the g_ij(dt)
                                              ! hop-probability channel; this member consumes no rng - the stream
                                              ! belongs to the matching member)
      end subroutine elec_prop_i
      subroutine mqc_match_i(sta)          ! state-statistics matching member
         import :: state_t
         type(state_t), intent(inout) :: sta  ! reads the electronic state + the member's
         ! own data channels; writes the held occ (occ_set - a hop = one component
         ! change of the occupation vector) and the hop counter (n_hop_set), sta%p
         ! (momentum projection), optionally sta%f (post-hop refresh: intra-box
         ! force_eval calls are legal), optionally the ontic container. Sole
         ! consumption point of the random stream during the evolution phase
      end subroutine mqc_match_i
      subroutine sync_proc_i(s)            ! container-synchronization member
         import :: elec_t
         type(elec_t), intent(inout) :: s  ! brings the method's declared deliverable
         ! containers up to date (a->rho / rho->a / none, method knowledge);
         ! the interface calls it as a pure-form tail and never interprets it
      end subroutine sync_proc_i
   end interface

   ! Method-package table (single-select by name at assembly; read-only during evolution)
   integer, parameter :: max_method = 16  ! package-table capacity [row]
   type :: method_entry_t
      character(len=32) :: name = ''      ! package name (the config%electronic%method assembly key)
      procedure(elec_prop_i), pointer, nopass :: prop  => null()  ! propagation member
      procedure(mqc_match_i), pointer, nopass :: match => null()  ! matching member
      procedure(sync_proc_i), pointer, nopass :: sync  => null()  ! container sync (null =
                                          ! no deliverable cache, legal)
   end type method_entry_t
   type(method_entry_t) :: method_tbl(max_method)  ! method-package registry (module-private)
   integer :: n_method = 0                ! registered row count [row]
   integer :: i_assembled = 0             ! assembled row index [-] (0 = not bound; bound at
                                          ! container_init from config%electronic%method)

   ! The statistics holder - run-lifetime module-private: occ written by the match
   ! member through occ_set, read by the recorder through occ_get; n_hop written
   ! through n_hop_set, read through n_hop_get. Seeded occ=[1] / n_hop=0 at method_bind
   integer, allocatable :: occ_held(:)   ! occupation vector (occ(k) = the configuration
                                        ! index; length 1 = the degenerate single-active-
                                        ! surface form; larger dimensions set at the
                                        ! container-side member init through occ_set)
   integer :: n_hop = 0             ! hop counter (pure observable, recorder-only read)

contains
   !------------------------------------------------------------------
   ! occ_get(occ) - read-only hand-out of the held occupation set
   !                (recorder columns)
   !------------------------------------------------------------------
   subroutine occ_get(occ)
      integer, allocatable, intent(out) :: occ(:) ! read-only copy of the held set (the
                                                  ! single writer stays the match member)
      ! 1. Unallocated holder (pre-assembly window - method_bind not yet run) ->
      ! named abort; otherwise copy the held occ into the dummy (allocate + assign)
      if (.not. allocated(occ_held)) then
         call stop_ele('occ_get', 'the occupation holder is not seeded (method_bind not '// &
                       'yet run - an evolution-phase read before assembly)')
      end if
      allocate (occ(size(occ_held)))
      occ = occ_held
   end subroutine occ_get

   !------------------------------------------------------------------
   ! occ_set(occ) - the match member's occupation write-back (a hop = one
   !                component change; reallocation legal - a multi-occupation
   !                method sets its dimension once at member init)
   !------------------------------------------------------------------
   subroutine occ_set(occ)
      integer, intent(in) :: occ(:)       ! the new occupation set (indices >= 1)
      integer :: k
      ! 1. Pre-assembly window (holder unseeded - method_bind not yet run) -> named abort;
      !    formal hygiene: every component names a state (>= 1)
      if (.not. allocated(occ_held)) then
         call stop_ele('occ_set', 'the occupation holder is not seeded (method_bind not '// &
                       'yet run - a statistics write before assembly)')
      end if
      do k = 1, size(occ)
         if (occ(k) < 1) then
            call stop_ele('occ_set', 'an occupation component names no state (index < 1 -'// &
                          ' a formal shape check)')
         end if
      end do
      if (allocated(occ_held)) deallocate (occ_held)
      allocate (occ_held(size(occ)))
      occ_held = occ
   end subroutine occ_set

   !------------------------------------------------------------------
   ! n_hop_get() - read-only hand-out of the hop counter (the recorder column)
   !------------------------------------------------------------------
   function n_hop_get() result(n)
      integer :: n                            ! the current hop count [count]
      ! 1. Hand out the held counter (a plain integer, always defined: seeded 0
      ! at method_bind; a pre-assembly read hands back the initialized 0)
      n = n_hop
   end function n_hop_get

   !------------------------------------------------------------------
   ! n_hop_set(n) - the match member's hop-counter write-back (the run
   !                accumulator; the member passes its own running count)
   !------------------------------------------------------------------
   subroutine n_hop_set(n)
      integer, intent(in) :: n                ! the hop count [count] (>= 0)
      if (n < 0) then
         call stop_ele('n_hop_set', 'a negative hop count (a formal check)')
      end if
      n_hop = n
   end subroutine n_hop_set

   !------------------------------------------------------------------
   ! method_reg(name, prop, match [, sync]) - append one row to the
   !                 method-package table method_tbl, storing the passed
   !                 procedures as the row's prop / match / sync procedure
   !                 pointers (the atomic registration unit; no row is
   !                 reachable until method_bind selects it)
   !------------------------------------------------------------------
   subroutine method_reg(name, prop, match, sync)
      character(len=*), intent(in) :: name  ! package name (the config%electronic%method key)
      procedure(elec_prop_i) :: prop        ! propagation member
      procedure(mqc_match_i) :: match       ! matching member
      procedure(sync_proc_i), optional :: sync ! container sync (absent = no deliverable cache)
      integer :: j
      ! 1. Duplicate name in method_tbl or full table (n_method = max_method) -> named abort
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
      ! 2. Register the row {name, prop, match, sync} (sync absent -> null
      !    pointer, legal: no deliverable container to keep up to date); executed
      !    by the reg files aggregated into registry_gen during the assembly phase
      n_method = n_method + 1
      method_tbl(n_method)%name = name
      method_tbl(n_method)%prop => prop
      method_tbl(n_method)%match => match
      if (present(sync)) then
         method_tbl(n_method)%sync => sync
      else
         method_tbl(n_method)%sync => null()
      end if
   end subroutine method_reg

   !------------------------------------------------------------------
   ! method_bind(name) - the assembly-time selection: resolve name
   !                 against the registered rows of method_tbl and record
   !                 the matched row index in the module-private integer
   !                 i_assembled; from then on elec_prop / mqc_statistic
   !                 dispatch through that row's prop / match / sync
   !                 procedure pointers. Called once from container_init;
   !                 a second call aborts. Also seeds the statistics
   !                 holders (occ = [1], n_hop = 0)
   !------------------------------------------------------------------
   subroutine method_bind(name)
      character(len=*), intent(in) :: name ! package name (config%electronic%method - a
                                           ! name string; rows come from methods/ or a
                                           ! systems/ folder alike - one registration
                                           ! site either way)
      integer :: j, i_row
      ! 1. Empty table (n_method = 0 - no library registered, registry_gen missing) or
      !    name not found -> named abort naming the requested method (a purely formal
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
      ! 2. Already bound (i_assembled /= 0) -> named abort (container_init ran twice - an
         ! assembly-choreography error, formal)
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
      !    dimensions are set by the container-side member init through occ_set
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
      !    not yet run) aborts with a named error (an evolution-phase entry called before assembly)
      if (i_assembled == 0) then
         call stop_ele('elec_prop', 'no method row bound (method_bind not yet run - '// &
                       'an evolution-phase entry called before assembly, formal)')
      end if
      ! 2. Dispatch: call the assembled row's prop member - electronic propagation
      !    given Q(t) (amplitude ODE along the trajectory / mapping flow /
      !    analytic elimination / no-op). The member reaches its Hamiltonian data
      !    through its own channels (container-folder private modules / level-2
      !    components) - dt semantics it cannot live with named abort INSIDE the
      !    member (the degenerate member ignores dt)
      call method_tbl(i_assembled)%prop(sta, dt)
      ! 3. Sync tail: if the row's sync pointer is attached, call sync(sta%s) -
      !    bring the declared deliverable containers up to date (a->rho /
      !    rho->a as the method declared; the interface never interprets the
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
                                           ! member's declaration, see step 2; the
                                           ! held occ feeds the statistics read-back)
      ! 1. Assembly-order guard: same formal check as elec_prop (an evolution-phase
      !    entry called before assembly)
      if (i_assembled == 0) then
         call stop_ele('mqc_statistic', 'no method row bound (method_bind not yet run - '// &
                       'an evolution-phase entry called before assembly, formal)')
      end if
      ! 2. Dispatch: call the assembled row's match member - the member owns the
      !    whole statistics write-back. It reads the electronic statistics from its
      !    own channels, samples (sole consumption point of the random stream during
      !    the evolution phase) and writes the classical side: the held occ through
      !    occ_set (a hop = one component change of the occupation vector) and the
      !    hop counter through n_hop_set (run-lifetime channels), sta%p (momentum
      !    projection), optionally sta%f (the post-hop force-refresh hook: the member
      !    may call force_eval itself - intra-box, legal) and optionally the ontic
      !    container. The degenerate member writes nothing: occ stays [1], n_hop
      !    stays 0, no rng consumption
      call method_tbl(i_assembled)%match(sta)
      ! 3. Sync tail: if the row's sync pointer is attached, call sync(sta%s) -
      !    same contract as elec_prop step 3 (both entries carry the tail:
      !    post-hop container resets and post-propagation refreshes are
      !    indistinguishable to the interface)
      if (associated(method_tbl(i_assembled)%sync)) then
         call method_tbl(i_assembled)%sync(sta%s)
      end if
   end subroutine mqc_statistic

   !------------------------------------------------------------------
   ! private helper: the named-abort channel of this interface (name the caller, print the
   ! named message, STOP 1)
   !------------------------------------------------------------------
   subroutine stop_ele(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': electronic-interface error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the method table is not usable)'
      stop 1
   end subroutine stop_ele
end module elec_interface
