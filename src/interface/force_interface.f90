!=====================================================================
! force_interface.f90 - force interface: force_eval (Q,P,s)->(F,E), the two
!   container slot-binding subroutines, and the single assembly entry
!   container_init
! Design:
!   The slots are module-private procedure pointers, each filled once at
!   assembly by the container's reg module: container_pes (the force
!   slot - the single force provider of the run, one container per
!   process) and container_term_proc (the optional termination slot).
!   The evaluation entries call through these pointers. force_eval is a
!   pure force entry: one call into the bound pointer's target, which
!   delivers the composed force evaluation and potential energy in eV
!   (the container-delivery contract: the container owns the force
!   composition - the active-surface fold, the member sum, stochastic
!   laws; an H-shaped container may call the level-2 library components
!   for the generic eigen-fold). The interface converts once and folds
!   E = T(P)+V. Members deliver in eV (the exchange unit); the
!   eV -> internal-unit conversion (23.0605 * e_conv) appears only in
!   this file. The interface is a purely formal authority
!   (signatures/units/slots/name dispatch) - physical-combination
!   legality is member knowledge.
!=====================================================================
module force_interface
   use consts,        only: e_conv
   use state,         only: state_t
   use config,        only: electronic
   use config_atoms, only: list_atoms
   use elec_interface,   only: method_bind
   implicit none
   private
   public :: force_eval, container_bind_pes, container_init
   public :: container_bind_term, container_term  ! the termination slot + its evaluation entry
   public :: container_probe_force                ! assembly-time point-force probe (spectra)

   ! eV to kcal/mol conversion (the unit-contract literal; appears only in this file)
   real(8), parameter :: ev_kcal = 23.0605d0

   ! The signature of what the force slot accepts - the type declaration of
   ! the procedure pointer container_pes below (a Fortran implementation
   ! detail; the pointer itself is the channel force_eval calls through).
   ! Delivery contract: g = the COMPOSED force kernel dV/dq of the active
   ! configuration [eV/A] (the container-side sum: the main PES's active
   ! surface + every additive stochastic-law / constraint member) and
   ! v = the matching potential-energy contribution [eV]. The bound entry
   ! FULLY ASSIGNS both (no read of the incoming values); stochastic
   ! members call rng internally inside the container (the single stream
   ! is guaranteed by rng)
   abstract interface
      subroutine container_pes_i(q, mass, g, v)
         real(8), intent(in)  :: q(:)     ! coordinates [Å] (all degrees of freedom, flattened)
         real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (stochastic/mass-law members;
                                          ! pure-PES members may ignore it)
         real(8), intent(out) :: g(:)     ! the composed force kernel dV/dq [eV/Å]
         real(8), intent(out) :: v        ! the potential-energy contribution [eV]
      end subroutine container_pes_i
   end interface

   ! The force slot: a module-private procedure pointer, the only channel
   ! through which force_eval reaches a container. Filled exactly once by
   ! container_bind_pes (assembly), called at every force_eval step. No
   ! registry, no name lookup - one container per process
   procedure(container_pes_i), pointer :: container_pes => null()

   ! The termination slot: an optional module-private procedure pointer holding
   ! the container's own trajectory-termination judgment (problem knowledge).
   ! Contract: read-only sta; geometry criteria are written container-side.
   ! Left null = legal (the MAX_STEPS-only run: the ever-present net is
   ! control%max_steps, checked in the driver loop condition)
   abstract interface
      logical function term_i(sta)
         import :: state_t
         type(state_t), intent(in) :: sta  ! physical state (read-only snapshot)
      end function term_i
   end interface
   procedure(term_i), pointer :: container_term_proc => null()

   ! Purely formal authority: no capability-combination table lives here. The
   ! assembly-time checks kept in container_init are purely formal (the slot is
   ! bound; the method name resolves). Physical-combination legality beyond form
   ! remains the members' own knowledge - they named abort at runtime.

contains
   !------------------------------------------------------------------
   ! container_bind_pes(proc) - fill the force slot: store the container's
   !                 aggregate entry proc into the module-private procedure
   !                 pointer container_pes; every later force_eval call
   !                 dispatches through that pointer. Executed once, at
   !                 assembly, by the container's reg module
   !                 (systems/<system>/reg_*.f90); a second call aborts
   !                 (one container per process)
   !------------------------------------------------------------------
   subroutine container_bind_pes(proc)
      procedure(container_pes_i) :: proc         ! the container's aggregate force entry (the pointer's target)
      if (associated(container_pes)) then
         call stop_ff('container_bind_pes', 'the force slot is already bound ('// &
                      'container_bind_pes ran twice - an assembly-choreography error, formal)')
      end if
      container_pes => proc
   end subroutine container_bind_pes

   !------------------------------------------------------------------
   ! container_bind_term(proc) - fill the termination slot: store the
   !                 container's trajectory-termination judgment proc into
   !                 the module-private procedure pointer
   !                 container_term_proc; the driver loop consults it
   !                 through container_term once per step. Optional - a
   !                 run that never fills it ends at max_steps only.
   !                 Executed once, at assembly, by the container's reg
   !                 module, alongside container_bind_pes; a second call
   !                 aborts
   !------------------------------------------------------------------
   subroutine container_bind_term(proc)
      procedure(term_i) :: proc              ! the container's termination judgment (the pointer's target)
      if (associated(container_term_proc)) then
         call stop_ff('container_bind_term', 'the termination slot is already bound ('// &
                      'container_bind_term ran twice - an assembly-choreography error, formal)')
      end if
      container_term_proc => proc
   end subroutine container_bind_term

   !------------------------------------------------------------------
   ! container_term(sta) - evolution-loop termination entry: .false. when
   !                 unbound (the MAX_STEPS-only run), else the
   !                 container's own judgment
   !------------------------------------------------------------------
   logical function container_term(sta)
      type(state_t), intent(in) :: sta       ! physical state (read-only snapshot)
      ! unbound slot = the MAX_STEPS-only run (legal): .false.; a bound slot
      ! runs the container's own judgment
      container_term = .false.
      if (.not. associated(container_term_proc)) return
      container_term = container_term_proc(sta)
   end function container_term

   !------------------------------------------------------------------
   ! container_init() - the single assembly entry point spanning both
   !                 interfaces: method lookup + purely formal validation
   !------------------------------------------------------------------
   subroutine container_init()
      ! 1. Method-slot assembly - runs in both postures (an idle chain with no
      !    container still binds the degenerate no-op method row); a pure name
      !    lookup - the named aborts live in method_bind itself
      call method_bind(electronic%method)
      ! 2. Idle-chain posture: an unbound force slot means no container was
      !    dragged in - assembly proceeds without binding; the named abort for a
      !    missing slot is force_eval's first call
      if (.not. associated(container_pes)) return
      ! 3. Container-side member init: executed by the container's reg module
      !    around its container_bind_pes call - not choreographed here
      ! 4. Integrator x electron-law legality is not checked here (a law's own
      !    prerequisites are its author's knowledge - the member aborts with a named error at
      !    runtime)
   end subroutine container_init

   !------------------------------------------------------------------
   ! force_eval(sta, e_tot) - the pure force entry (Q,P,s)->(F,E): the
   !                            container delivers the composed kernel ->
   !                            interface converts (the unique unit-contract
   !                            point, ev_kcal) -> writes sta%f, e_tot;
   !                            E = T(P) + V
   !------------------------------------------------------------------
   subroutine force_eval(sta, e_tot)
      type(state_t), intent(inout) :: sta  ! physical state (q/p/s read in; f written out - the composed force)
      real(8), intent(out) :: e_tot        ! total energy H = T(P)+V [kcal/mol]
      real(8) :: g(size(sta%q))            ! the delivered composed kernel [eV/A]
      real(8) :: v_sum                     ! the delivered potential contribution [eV]
      real(8) :: v_pot, t_kin              ! energy channels [kcal/mol]
      integer :: k
      ! 1. Unbound slot -> named abort (an evolution-phase entry called before
      !    assembly - an assembly-order error, formal)
      if (.not. associated(container_pes)) then
         call stop_ff('force_eval', 'the force slot is not bound (container_bind_pes not '// &
                      'yet run - an evolution-phase entry called before assembly)')
      end if
      ! 2. One call into the container (the composed sum executes inside it -
      !    the active-surface fold and every additive member are container-side;
      !    the delivery contract fully assigns g and v
      call container_pes(sta%q, list_atoms%mass, g, v_sum)
      ! 3. Unit-contract fulfillment point (unique in the whole library):
      !    sta%f = -g*ev_kcal*e_conv (the one eV->internal conversion)
      sta%f = 0.0d0
      sta%f = sta%f - g*ev_kcal*e_conv
      ! 4. No electron-law consumption here (that lives in elec_prop /
      !    mqc_statistic in elec_interface)
      ! 5. Energy [kcal/mol]: v_pot = v_sum*ev_kcal; t_kin = sum p^2/(2*list_atoms%mass)/e_conv
      !    (per-atom fold); e_tot = t_kin + v_pot (H = T + V)
      v_pot = v_sum*ev_kcal
      t_kin = 0.0d0
      do k = 1, size(list_atoms%mass)
         t_kin = t_kin + dot_product(sta%p(3*k-2:3*k), sta%p(3*k-2:3*k))/(2.0d0*list_atoms%mass(k))
      end do
      t_kin = t_kin/e_conv
      e_tot = t_kin + v_pot
   end subroutine force_eval

   !------------------------------------------------------------------
   ! container_probe_force(q, f) - the assembly-time point-force probe of
   !                 the bound slot: one call through container_pes, the
   !                 result folded EXACTLY like force_eval (f = -g with
   !                 the unique eV->internal fold staying in this file).
   !                 Consumed by the spectrum derivation (Hessian
   !                 differencing needs a point force in internal-per-A,
   !                 not a state); list_atoms masses ride along per the
   !                 container delivery contract
   !------------------------------------------------------------------
   subroutine container_probe_force(q, f)
      real(8), intent(in)  :: q(:)        ! coordinates [A] (all dof, size = 3*list_atoms atoms)
      real(8), intent(out) :: f(:)        ! forces [internal/A] (the fold of -g [eV/A])
      real(8) :: g(size(q)), v
      if (.not. associated(container_pes)) then
         call stop_ff('container_probe_force', 'the force slot is not bound (the probe is'// &
                      ' an assembly-phase entry called before assembly)')
      end if
      if (size(q) /= 3*size(list_atoms%mass) .or. size(f) /= 3*size(list_atoms%mass)) then
         call stop_ff('container_probe_force', 'the probe vector must cover the atom list'// &
                      ' atoms exactly (3 per atom)')
      end if
      call container_pes(q, list_atoms%mass, g, v)
      f = -g*ev_kcal*e_conv
   end subroutine container_probe_force

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! stop_ff(who, msg) - the named-abort channel (name the caller, print the named
   ! message, STOP 1)
   subroutine stop_ff(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': container error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container is not usable)'
      stop 1
   end subroutine stop_ff
end module force_interface
