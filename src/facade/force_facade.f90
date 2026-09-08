!=====================================================================
! force_facade.f90 - L2 interface-A facade: (Q,P,s)->(F,E). Inside: sum of class-A force members
!                    (additive) x the single-selected class-B electron-law member; the unit-contract
!                    fulfillment point (x23.0605 x e_conv appears only here)
!                    (see previous line)
! Reference: venus2026/DVDQ.f:10-22 (PDOT zeroing + DPESHON dispatch - the old form of additive composition);
!       src_TEST/interface_TEST.f90:12-16 (POTPRE member-init hook), :18-30
!       (POT0 energy conversion x23.0605), :32-44 (DPESHON gradient accumulation
!       PDOT(i)=PDOT(i)-g(i)*23.0605*C1 - the new facade conversion is character-for-character identical,
!       N1.1 bit-comparable); ENERGY.f:54-68 (T=sum p^2/2m/C1, V=POT0 convention, H=T+V);
!       docs/final_completion_plan.md §5 N1 (registration/assembly/capability negotiation-table spec - directly in force)
! Member unit contract: class-A members deliver force kernels/energies in eV (eV is the exchange unit -
! following the old interface_TEST v/g chain); eV -> internal-unit conversion is centralized in this
! facade (the T4 grep requirement that 23.0605 appear only in this file - hence the "F_add code units" note
! in the N1 member-signature spec reads as "after facade conversion")
! (see previous line)
! (see previous line)
! Member tables are filled in C2 (the six pes_* files deferred - K2 plan E3; main PES + additive stochastic-law members)
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
! E9 retained assets: elec_t amplitudes/hop counter, the class-B dt channel, the negotiation table as
! data (triad_tbl)
! (the pes_kind column now means the ws fill mode repr; the radau-false-on-multi-surface-rows constraint
!   stays - at an adaptive sequence boundary "one step" is ill-defined, so the hop probability g_ij(dt)
!   loses its basis)
!   (see previous line)
!=====================================================================
module force_facade
   use consts,        only: e_conv
   use state,         only: state_t, elec_t
   use config,        only: electronic, propagator
   use config_roster, only: roster
   use elec_ws,       only: ws_t, r_scalar, r_diab, r_adiab
   implicit none
   private
   public :: facade_reg_pes
   public :: facade_reg_elec, facade_init, facade_force, facade_info

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
                                          ! i_want, the member reads it and fills per its own mode the
                                          ! delivery fields - scalar fills g_add/v_add;
                                          ! diab fills h_ij/dh_ij; adiab fills
                                          ! v_surf/g_surf/d_ij [the i_want row];
                                          ! see the container contract in elec_ws.f90)
      end subroutine pes_proc_i
      ! Note: stochastic-law members call rng internally (N1 - the signature acknowledges stochastic laws;
      ! the single stream is guaranteed by rng)
   end interface

   ! Class-B member (single-select) interface: electron propagation law - electronic-state evolution +
   ! nuclear-force modification (adiabatic = two degenerate points)
   abstract interface
      subroutine elec_proc_i(s, f, dt, ws)
         import :: elec_t, ws_t
         type(elec_t), intent(inout) :: s  ! electronic state (the law member evolves amp/i_surf/
                                           ! n_hop per its law; adiabatic: unchanged)
         real(8), intent(inout) :: f(:)    ! nuclear forces [internal energy/Å] (modified by the law member -
                                           ! surface selection/averaging/switching and kinetic-energy
         ! rescaling; adiabatic: identity)
         real(8), intent(in)    :: dt      ! current actual step [10 fs] (needed for the hop probability
         ! g_ij(dt) - E9; adiabatic members ignore it)
         type(ws_t), intent(in)  :: ws     ! workspace (read-only consumption: interpret the container per
                                           ! ws%repr - in-box wiring, E10; a repr the member does not
                                           ! support -> loud-stop, e.g. adiab x iesh)
      end subroutine elec_proc_i
   end interface

   ! Registries (additive class-A table + single-select class-B slot; registered in the assembly phase,
   ! read-only during evolution)
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
   character(len=32) :: elec_name = ''    ! class-B member name (module-private; single-select slot; empty = not assembled)
   procedure(elec_proc_i), pointer :: elec_proc ! class-B routine (module-private; attached by
                                                ! facade_init in the assembly phase; nopass is a component attribute only)

   ! Legal capability-combination table (E9 as data; E10: the pes_kind column now means the ws fill mode -
   ! adding a row = adding data, not code). Four axes: H-capability n_surf interval x main-PES repr mode
   ! x electron law x integrator (the integrator axis is carried by three booleans). The radau column is
   ! permanently false on multi-surface rows: "one step" is ill-defined at an adaptive sequence boundary,
   ! so the hop probability g_ij(dt) loses its basis - nonadiabatic only pairs with fixed-step integrators
   ! (verlet/symple).
   ! (see previous line)
   ! (see previous line)
   ! (see previous line)
   ! (see previous line)
   integer, parameter :: m_adiabatic = 1  ! electron-law code: adiabatic degenerate point (config%electronic%method key)
   integer, parameter :: m_tsh      = 2   ! electron-law code: TSH/FSSH (future elec_tsh registration key)
   integer, parameter :: m_iesh     = 3   ! electron-law code: IESH (future elec_iesh registration key)
   integer, parameter :: max_triad = 4    ! legal-combination table capacity [row]
   integer, parameter :: n_surf_unbounded = huge(1) ! capability upper bound on multi-surface rows (meaning unbounded)
   type :: triad_row_t
      integer :: n_surf_lo = 1            ! H-capability lower bound (n_surf >= this)
      integer :: n_surf_hi = 1            ! H-capability upper bound (single-surface rows 1..1; multi-surface rows .. unbounded)
      integer :: pes_kind  = r_scalar     ! main-PES ws fill mode (r_* - the E10 repr semantics)
      integer :: method    = m_adiabatic  ! electron-law code (m_*)
      logical :: ok_verlet = .true.       ! integrator 1 (verlet/beeman) allow bit
      logical :: ok_symple = .true.       ! integrator 2 (symple) allow bit
      logical :: ok_radau  = .true.       ! integrator 3 (radau) allow bit (false on multi-surface rows)
   end type triad_row_t
   type(triad_row_t), parameter :: triad_tbl(max_triad) = (/ &
      ! Row 1 (current): single-surface scalar mode x adiabatic x all three integrators allowed
      triad_row_t(1, 1, r_scalar, m_adiabatic, .true., .true., .true.), &
      ! Row 2: multi-surface diabatic mode x TSH x fixed-step integrators
      triad_row_t(2, n_surf_unbounded, r_diab,  m_tsh, .true., .true., .false.), &
      ! Row 3: multi-surface adiabatic mode x TSH x fixed-step integrators
      triad_row_t(2, n_surf_unbounded, r_adiab, m_tsh, .true., .true., .false.), &
      ! Row 4: multi-surface diabatic mode x IESH x fixed-step integrators (IESH naturally consumes the
      ! diabatic container - model Hamiltonians; adiabatic x IESH is not allowed)
      ! (see previous line)
      triad_row_t(2, n_surf_unbounded, r_diab,  m_iesh, .true., .true., .false.) &
      /)
   ! Note: the method of rows 2-4 (tsh/iesh) currently has no class-B member registered - the negotiation
   ! table allows first, and facade_init's class-B assembly (slot lookup by method) loud-stops as the
   ! backstop: a double gate (see the triad_tbl note)
   ! (see previous line)
contains
   !------------------------------------------------------------------
   ! facade_reg_pes(name, proc[, repr]) - register a class-A force member (additive: main PES + stochastic
   ! laws / constraint potentials; E10-unified single entry - representation via the repr data tag, no
   ! longer a registration fork)
   ! (see previous line)
   ! Reference: docs/final_completion_plan.md §5 N1 (registration spec); T2 (new member = new file +
   !       one registration line); E10 adjudication (single signature fully unified - E9 dual signatures withdrawn)
   subroutine facade_reg_pes(name, proc, repr)
      character(len=*), intent(in) :: name  ! member name (assembly key - config%electronic%pes_main)
      procedure(pes_proc_i) :: proc         ! member evaluation routine (single signature - scalar/diab/
                                            ! adiab unified, delivered via the ws container)
      integer, intent(in), optional :: repr ! ws fill mode r_* (default r_scalar -
                                            ! the six single-surface members need zero registration-line
                                            ! changes; multi-surface main PES registers with
      ! r_diab/r_adiab)
      ! (see previous line)
      ! Pseudocode:
      ! 1. Duplicate name or full table -> loud-stop
      ! 2. Register the row {name, kind=repr (present-check defaults to r_scalar), proc} (executed in the
      !    assembly phase by the driver; multi-surface pes members arrive with the elec_tsh consumers)
   end subroutine facade_reg_pes

   !------------------------------------------------------------------
   ! facade_reg_elec(name, proc) - register a class-B electron-law member (single select)
   ! Reference: docs/final_completion_plan.md §5 N1 (class-B single-select spec)
   subroutine facade_reg_elec(name, proc)
      character(len=*), intent(in) :: name  ! member name (the config%electronic%method mapping key)
      procedure(elec_proc_i) :: proc        ! member routine
      ! Pseudocode:
      ! 1. Slot already occupied -> loud-stop (class B is single-select - duplicate registration is an
      !    assembly error)
      ! 2. Fill the slot {name, proc} (the elec_adiabatic degenerate member is registered in the assembly
      !    phase)
   end subroutine facade_reg_elec

   !------------------------------------------------------------------
   ! facade_init() - member assembly + capability-combination validation (illegal combinations loud-stop)
   ! Reference: docs/final_completion_plan.md §5 N1 (negotiation-table spec); E9 (table as data -
   !       triad_tbl four rows: current single-surface adiabatic + multi-surface tsh/iesh pre-allowed);
   !       src_TEST/interface_TEST.f90:12-16 (POTPRE-type member-init hook - the D3
   !       ghost parameter assembly site); docs/design_spec_completeness.md §3.2-2
   !       (type-checking of capability and method legality; illegal combinations rejected
   !        explicitly with both conflicting axes named)
   subroutine facade_init()
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
      !    alone, see facade_force step 4')
      !    (see previous line)
      ! 2. Class-B assembly: select the law member per config%electronic%method (1 = adiabatic ->
      !    the elec_adiabatic degenerate member; 2/3 = tsh/iesh registered later); unregistered -> loud-stop
      !    (negotiation table allows first, empty registration slot blocks after - a double gate, see the
      !    triad_tbl note)
      ! 3. Member-init hooks (POTPRE type): call each assembled member's init (parameters arrive through
      !    the member's own params - D3; the hook-slot wiring shape is finalized in C2)
      !    (see previous line)
      ! 4. Capability-combination validation (E9 table lookup): match the 4-tuple
      !    (n_surf in row interval) x (main-PES ws mode = row pes_kind) x (method = row method)
      !    x (integrator id allow bits - config%propagator%id mapped to ok_verlet/ok_symple/
      !    ok_radau) against triad_tbl row by row; no matching row -> print the current values of the
      !    conflicting axes and STOP
      !    (see previous line)
      !    (spec §3.2-2 type check - both axes named; typical interception: multi-surface x radau -
      !    the hop probability g_ij(dt) loses its basis at adaptive sequence boundaries)
   end subroutine facade_init

   !------------------------------------------------------------------
   ! facade_force(sta, e_tot[, dt]) - (Q,P,s)->(F,E): members deliver via ws -> facade composes ->
   !                            converts -> electron-law consumption; E = T(P) + V (the E10 unified path)
   ! Reference: venus2026/DVDQ.f:10-22 (zeroing + dispatch); src_TEST/interface_TEST.f90:32-44
   !       (identical-form accumulation), :18-30 (energy convention); ENERGY.f:54-68 (T/H reporting convention)
   subroutine facade_force(sta, e_tot, dt)
      type(state_t), intent(inout) :: sta  ! physical state (q/p/s read in; f written out - the composed force)
      real(8), intent(out) :: e_tot        ! total energy H = T(P)+V [kcal/mol] (N1 convention;
                                           ! multi-surface mode V = current surface / law-member convention, see 4')
      real(8), intent(in), optional :: dt  ! current actual step [10 fs] (the hop probability g_ij(dt)
                                           ! channel - E9; the default takes config%propagator%dt,
                                           ! fixed-step semantics; the adiabatic law ignores it)
                                           ! (see previous line)
      ! Facade-local workspace (E10 axioms B/C: declared in the box, dies in the box, one-step lifetime -
      ! invisible outside the box (driver/propagator/output), zero footprint in the facade_force signature)
      ! (see previous line)
      type(ws_t) :: ws                     ! electronic workspace (main-PES multi-surface constructs +
      ! scalar-member contribution container; member filling and law
      ! consumption share one channel)
      ! (see previous line)
      ! Pseudocode:
      ! 1. Zero the force-kernel workspace g(:), v_sum = 0 (old DVDQ.f:12-14 PDOT zeroing site);
      !    ws%repr = r_empty (reset every step - axiom C one-step lifetime); dt_eff = the default
      !    config%propagator%dt (alias prop_cfg) or the passed-in value
      ! 2. Main-PES evaluation (the assembled main-PES member - one unified call): ws%i_want = sta%s%i_surf
      !    (the current-surface coupling-row request in adiab mode - i_surf stays out of the class-A
      !    signature and is conveyed through the request field, E10); proc(sta%q, roster%mass, ws) -
      !    the member fills ws per its own mode
      !    (see previous line)
      !    (scalar: g_add/v_add; diab: h_ij/dh_ij; adiab: v_surf/g_surf/
      !    d_ij [the i_want row])
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
      ! 5. Electron-law consumption: elec_proc(sta%s, sta%f, dt_eff, ws) (read-only consumption -
      !    interpret the container per ws%repr; a repr the member does not support -> loud-stop. The
      !    adiabatic member is the identity and ignores ws - under a single-surface representation the
      !    composite collapses to a pure PES, spec §3.2-2; a nonadiabatic law completes the amplitude
      !    propagation, the hop decision (single rng stream), and the force modification / kinetic-energy
      !    rescaling after the i_surf switch here - consumption details are fixed with the elec_tsh/
      !    elec_iesh members (phase C))
      !    (see previous line)
      !    (see previous line)
      !    (see previous line)
      !    (see previous line)
      ! 6. Energy: v_pot = v_sum*ev_kcal [kcal/mol] (old POT0 :26 convention);
      !    t_kin = sum p^2/(2*roster%mass)/e_conv [kcal/mol] (per-atom fold -
      !    old ENERGY.f:58-66 convention); e_tot = t_kin + v_pot (old :68 H=T+V)
      ! ————————————————————————————————————————————————
      ! 4'. Multi-surface branch (main-PES row mode = r_diab/r_adiab - a composition swap of steps 3-6):
      !    a. Force-kernel composition: r_diab -> diagonalize ws%h_ij via utils/linalg to compute the
      !       adiabatic view on the spot (eigenvalues = per-surface V, eigenvectors = representation
      !       transform U); the current-surface adiabatic force kernel takes the transformed current row;
      !       the dh_ij.U product needed for the hop probability g_ij(dt) is fetched by the class-B member
      !       itself (the facade does not pre-build it); r_adiab -> g = ws%g_surf(:, sta%s%i_surf) (the
      !       current-surface row - after a hop the class-B member updates i_surf, and recomputing with
      !       the new row belongs to the next evaluation)
      !       (see previous line)
      !       (see previous line)
      !       (see previous line)
      !       (see previous line)
      !    b. v_pot convention: adiab -> v_surf(i_surf)*ev_kcal; diab -> population-weighted or
      !       current-surface (class-B law-member convention, fixed by the consumer); scalar additive
      !       contributions fold in as usual;
      !       (see previous line)
      !       e_tot = t_kin + v_pot
      !    c. radau note: permanently false on multi-surface rows (triad_tbl) - adaptive substeps have no
      !       well-defined "one step" dt semantics; this path only ever sees fixed steps (same source as
      !       the propagator header note)
      !       (see previous line)
   end subroutine facade_force

   !------------------------------------------------------------------
   ! facade_info() - diagnostic output: member list and capability combination
   ! Reference: docs/final_completion_plan.md §5 N1 (facade_info - member list / capability levels)
   subroutine facade_info()
      ! Pseudocode:
      ! 1. Print the class-A registered/assembled member list (name + main-PES/additive tag + ws fill
      !    mode r_*) and the class-B member name
      !    (see previous line)
      ! 2. Print the current capability combination: H-capability n_surf x main-PES ws mode x electron
      !    law x integrator (from the three config sources) + the matched triad_tbl row number (shown
      !    empty if not assembled)
      !    (see previous line)
      ! 3. Print the unit-contract constants (ev_kcal/e_conv) and the assembly state - for debugging
      !    assembly errors
   end subroutine facade_info
end module force_facade
