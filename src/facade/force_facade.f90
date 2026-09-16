!=====================================================================
! force_facade.f90 - container class-A facade: force_eval (Q,P,s)->(F,E),
!   container slot binding, and the single assembly entry container_init
! Design:
!   force_eval is a pure force entry: one call into the container's bound
!   aggregate PES entry, then F and E = T(P)+V are folded from the delivered
!   workspace. Members deliver in eV (the exchange unit); the eV ->
!   internal-unit conversion (23.0605 * e_conv) appears only in this file.
!   The facade is a purely formal authority (signatures/units/containers/
!   name dispatch) - physical-combination legality is member knowledge.
!   One container per run: a single class-A slot + optional termination slot.
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
   public :: container_bind_term, container_term  ! the termination slot + its evaluation entry

   ! eV to kcal/mol conversion (the unit-contract literal; appears only in this file)
   real(8), parameter :: ev_kcal = 23.0605d0

   ! The container aggregate entry - the single class-A slot (module-private interface:
   ! the pointer's type declaration, a Fortran implementation detail). Delivery
   ! contract: the scalar channel g_add/v_add = the container-side sum in eV (the main
   ! PES in scalar mode + every additive stochastic-law / constraint member; in
   ! multi-surface mode the additive sum alone - the main PES rides the construct
   ! channel); the construct channel = the main PES's multi-surface deliveries
   ! (h_ij/dh_ij or v_surf/g_surf/d_ij) per its self-declared fill mode
   abstract interface
      subroutine container_pes_i(q, mass, ws)
         import :: ws_t
         real(8), intent(in)  :: q(:)     ! coordinates [Å] (all degrees of freedom, flattened)
         real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (stochastic/mass-law members;
                                          ! pure-PES members may ignore it)
         type(ws_t), intent(inout) :: ws  ! workspace (inout: the i_want(:) request field is
                                          ! pre-written; the container fills per its own mode)
      end subroutine container_pes_i
   end interface
   ! Note: stochastic-law members call rng internally (the single stream is guaranteed
   ! by rng) - this holds inside the container exactly as it held per member

   ! The bound slot (a single slot, no registry, no name lookup; one container per run)
   procedure(container_pes_i), pointer :: container_pes => null()
   character(len=32) :: h_deliver = ''  ! the container's declared H delivery set (comma-list
                                        ! of r_* mode names, e.g. 'scalar,adiab'; the
                                        ! assembly-time match target of the method row's
                                        ! h_form)

   ! The termination slot - the container's own trajectory-termination judgment
   ! (problem knowledge). Contract: read-only sta; geometry criteria are written
   ! container-side. Unbound slot = legal (the safety-net-only run: the ever-present
   ! net is control%max_steps, checked in the driver loop condition)
   abstract interface
      logical function term_i(sta)
         import :: state_t
         type(state_t), intent(in) :: sta  ! physical state (read-only snapshot)
      end function term_i
   end interface
   procedure(term_i), pointer :: container_term_proc => null()

   ! Purely formal authority: no capability-combination table lives here. The
   ! assembly-time checks kept in container_init are purely formal (the slot is bound;
   ! the method name resolves; a multi-surface system names a multi-surface delivery
   ! mode; the method row's h_form demand is a subset of the container's h_deliver
   ! supply). Physical-combination legality beyond form remains the members' own
   ! knowledge - they loud-stop at runtime.

contains
   !------------------------------------------------------------------
   ! container_bind_pes(proc [, h_deliver]) - assembly injection: bind the
   !                 container's aggregate class-A entry + its declared H
   !                 delivery set (the caller is the container's reg module)
   !------------------------------------------------------------------
   subroutine container_bind_pes(proc, h_deliver)
      procedure(container_pes_i) :: proc         ! the container aggregate entry (the single class-A slot)
      character(len=*), intent(in), optional :: h_deliver ! declared delivery set (comma-list of
                                                          ! r_* names; default '' = scalar only)
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
   ! container_bind_term(proc) - assembly injection: bind the container's
   !                 trajectory-termination judgment (the caller: the
   !                 container's reg module, alongside container_bind_pes)
   !------------------------------------------------------------------
   subroutine container_bind_term(proc)
      procedure(term_i) :: proc              ! the container termination judgment
      if (associated(container_term_proc)) then
         call stop_ff('container_bind_term', 'the termination slot is already bound ('// &
                      'container_bind_term ran twice - an assembly-choreography error, formal)')
      end if
      container_term_proc => proc            ! executed once at assembly by the
                                             ! container's reg module
   end subroutine container_bind_term

   !------------------------------------------------------------------
   ! container_term(sta) - evolution-loop termination entry: .false. when
   !                 unbound (the safety-net-only run), else the
   !                 container's own judgment
   !------------------------------------------------------------------
   logical function container_term(sta)
      type(state_t), intent(in) :: sta       ! physical state (read-only snapshot)
      ! unbound slot = the safety-net-only run (legal): .false.; a bound slot
      ! runs the container's own judgment
      container_term = .false.
      if (.not. associated(container_term_proc)) return
      container_term = container_term_proc(sta)
   end function container_term

   !------------------------------------------------------------------
   ! container_init() - the single assembly entry point spanning both
   !                 facades: method lookup + purely formal validation
   !------------------------------------------------------------------
   subroutine container_init()
      character(len=32) :: hf                 ! the bound method row's Hamiltonian preset
      ! 1. Method-slot assembly - runs in both postures (an idle chain with no
      !    container still binds the degenerate no-op method row); a pure name
      !    lookup - the loud-stops live in method_bind itself
      call method_bind(electronic%method)
      ! 2. Idle-chain posture: an unbound class-A slot means no container was
      !    dragged in - assembly proceeds without binding; the loud-stop for a
      !    missing slot is force_eval's first call. The checks below run only
      !    when a container is bound
      if (.not. associated(container_pes)) return
      ! 3. Container-contract formal check: a multi-surface system needs a
      !    multi-surface delivery mode
      if (electronic%n_surf >= 2) then
         if (.not. names_mode(h_deliver, 'diab') .and. .not. names_mode(h_deliver, 'adiab')) then
            call stop_ff('container_init', 'multi-surface system (n_surf >= 2) but the '// &
                         'declared delivery set names no multi-surface mode: "'//trim(h_deliver)//'"')
         end if
      end if
      ! 4. H-form match: the method row's demand must be a subset of the
      !    container's supply; an empty h_form (the degenerate method) always
      !    passes
      hf = method_h_form()
      if (.not. subset_of(hf, h_deliver)) then
         call stop_ff('container_init', 'the method Hamiltonian preset "'//trim(hf)// &
                      '" is not covered by the container delivery set "'//trim(h_deliver)//'"')
      end if
      ! 5. Container-side member init: executed by the container's reg module
      !    around its container_bind_pes call - not choreographed here
      ! 6. Observable-vs-container formal check: deferred to the interface
      !    phase; integrator x electron-law legality is not checked here (a
      !    law's own prerequisites are its author's knowledge - the member
      !    loud-stops at runtime)
   end subroutine container_init

   !------------------------------------------------------------------
   ! force_eval(sta, e_tot) - the pure force entry (Q,P,s)->(F,E): the
   !                            container delivers via ws -> facade converts
   !                            (the unique unit-contract point, ev_kcal) ->
   !                            writes sta%f, e_tot; E = T(P) + V
   !------------------------------------------------------------------
   subroutine force_eval(sta, e_tot)
      type(state_t), intent(inout) :: sta  ! physical state (q/p/s read in; f written out - the composed force)
      real(8), intent(out) :: e_tot        ! total energy H = T(P)+V [kcal/mol] (multi-surface
                                           ! mode: V = current-surface convention, see 4')
      ! Workspace: the elec_facade-held instance, aliased through the pointer function
      ! ws_handle - the one sanctioned alias (an ordinary dummy would detach the fill
      ! from the holder by copy semantics). The constructs stay in place for
      ! elec_prop / mqc_statistic of the same integration step (one-step lifetime;
      ! never an argument of any public signature)
      type(ws_t), pointer :: w             ! the holder alias (ws_handle)
      integer, allocatable :: occ_req(:)   ! the delivery request - the held occupation set
      integer :: k, r, i, j, n
      real(8) :: g(size(sta%q))            ! folded force kernel [eV/A]
      real(8) :: v_sum, v_pot, t_kin       ! energy channels [eV / kcal/mol]
      real(8), allocatable :: hh(:,:), hev(:)  ! the diagonalization pair: the H copy
                                              ! (destroyed by eig_sym - it returns the
                                              ! eigenvectors BY COLUMN) + the eigenvalues
      ! 1. Zero the accumulators; ws_reset() - every evaluation begins from the
      !    empty container (with elec_facade's step-end reset the two sites
      !    together make the one-step lifetime robust to skipped phases)
      call ws_reset()
      g = 0.0d0
      v_sum = 0.0d0
      ! 2. Unbound slot -> loud-stop (an evolution-phase entry called before
      !    assembly - an assembly-order error, formal)
      if (.not. associated(container_pes)) then
         call stop_ff('force_eval', 'the class-A slot is not bound (container_bind_pes not '// &
                      'yet run - an evolution-phase entry called before assembly)')
      end if
      ! 3. One call into the container (the member sum executes inside it),
      !    through the holder alias: the request field i_want(:) = occ_get()
      !    (the delivery request; length 1 = the single-active-surface
      !    degenerate form) is pre-written into the holder, and the container
      !    fills it per its own declared mode
      call occ_get(occ_req)
      w => ws_handle()
      allocate (w%i_want(size(occ_req)))
      w%i_want = occ_req
      call container_pes(sta%q, roster%mass, w)   ! the container fills ws per its own mode
      ! 4. Fold the delivered channel, branching on the container's
      !    self-declared fill mode ws%repr (representation is data; the runtime
      !    tag remains the authority): every not-delivered or not-consumed
      !    mode is a named loud-stop, never an undefined read
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
         ! 4'. Multi-surface diabatic branch: diagonalize the delivered H_ij on
         !     the spot (eigenvalues = per-surface V, eigenvectors = the
         !     representation transform U) and fold the OCCUPATION ROWS of dh_ij:
         !     d(lambda_k)/dq = sum over ordered pairs (i,j) of U(i,k) U(j,k)
         !     (dH/dq)_ij - the occupied eigenvalue's gradient (a symmetric dh_ij
         !     delivery counts every off-diagonal twice; eigenvector sign
         !     ambiguity cancels in the U U products). The diagonalization is
         !     system-independent generic physics and stays in the facade. v_pot
         !     convention: the current-occupation surfaces; the additive scalar
         !     kernel folds on top of the construct-channel kernel
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
                              ! delivered constructs must survive in the holder for
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
      ! 5. Unit-contract fulfillment point (unique in the whole library):
      !    sta%f = -g*ev_kcal*e_conv (the one eV->internal conversion; how g
      !    is composed in multi-surface mode: 4')
      sta%f = 0.0d0
      sta%f = sta%f - g*ev_kcal*e_conv
      ! 6. No electron-law consumption here (that lives in elec_prop /
      !    mqc_statistic in elec_facade) - the composed ws of the current
      !    geometry is what the later phases next read
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
