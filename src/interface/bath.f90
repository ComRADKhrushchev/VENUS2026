!=====================================================================
! bath.f90 - the pre-evolution bath loop: species registration,
!   provisioning, and the loop itself (propagation + collision sweeps
!   + one exact whole-system drift removal)
! Design:
!   A bath is NOT a per-fragment distribution member: its sampling power
!   comes from repetition in time, so it runs as ONE loop inside samp_run's
!   trajectory head - after the member realizations (the interiors exist)
!   and before the orientation stage and the incident channel (the bath
!   must not see the beam momenta; a rotation after it would be neutral,
!   the incident placement after it re-places the COMs). bath_reg(name,
!   sweep) registers a species (the methods/ wildcard picks its reg file);
!   bath_setup reads the BATH provision (BATH = NONE stays bitwise inert:
!   no RNG draws, no state writes - the archived stream is untouched);
!   bath_run(sta) executes N_BATH steps of [prop_step at DT_BATH + the
!   species sweep], saves/restores the global step and trajectory time
!   around it, and ends with ONE exact whole-system COM-momentum
!   projection (the collisions inject a random net momentum; the incident
!   channel prescribes the beam translation, so the residual is stripped).
!   Units: dt [10 fs]; p [amu*Ang/(10 fs)].
!=====================================================================
module bath
   use state,      only: state_t
   use config,     only: reactants, prop_cfg => propagator
   use config_atoms, only: list_atoms
   use propagator, only: prop_step
   use geom_init,  only: geom_sep_probe
   use force_interface, only: container_probe_force
   implicit none
   private
   public :: bath_reg, bath_setup, bath_run, bath_armed

   integer, parameter :: max_bath = 4     ! species-table capacity [row]
   real(8), parameter :: sep_tol = 1.0d-6  ! fragment-separation interaction
                                           ! tolerance [internal ~ eV/Angstrom]

   ! Bath species interface: one collision sweep at step size dt over the
   ! masked atoms (the mask is the BATH_SCOPE resolution - all atoms, or the
   ! non-fragment atoms of the surface paradigm)
   abstract interface
      subroutine bath_sweep_i(sta, dt, mask)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! physical state (the masked atoms' momenta are
                                             ! reset on collision)
         real(8), intent(in) :: dt            ! the sweep's step size [10 fs] (the collision
                                             ! probability is nu*dt)
         logical, intent(in) :: mask(:)       ! per-atom collision scope [flag]
      end subroutine bath_sweep_i
   end interface

   ! Species table (the registry family - a new species = a new methods/
   ! file + one registration line)
   type :: bath_entry_t
      character(len=32) :: name = ''                       ! species name (the BATH word)
      procedure(bath_sweep_i), pointer, nopass :: sweep => null() ! one collision sweep
   end type bath_entry_t
   type(bath_entry_t) :: bath_tbl(max_bath)  ! species table (module-private)
   integer :: n_species = 0                  ! registered row count [row]

   ! Provisioned state (module-private; set by bath_setup)
   procedure(bath_sweep_i), pointer :: sweep_c => null()          ! the bound species sweep
   integer :: n_step_c = 0                    ! provisioned step count [-]
   real(8) :: dt_c = 0.0d0                    ! provisioned step size [10 fs]
   logical :: armed = .false.                 ! the BATH provision armed this bath [flag]
   logical :: setup_done = .false.            ! bath_setup ran [flag]
contains
   !------------------------------------------------------------------
   ! bath_reg(name, sweep) - register a bath species (the BATH enum word
   !                         names it; duplicate names and capacity abort)
   !------------------------------------------------------------------
   subroutine bath_reg(name, sweep)
      character(len=*), intent(in) :: name   ! species name (must equal a BATH word)
      procedure(bath_sweep_i) :: sweep        ! one collision sweep
      integer :: i
      do i = 1, n_species
         if (trim(bath_tbl(i)%name) == trim(name)) then
            call stop_bath('bath_reg', 'bath species "'//trim(name)// &
               '" is already registered (duplicate registration line)')
         end if
      end do
      if (n_species >= max_bath) then
         call stop_bath('bath_reg', 'the bath species table is full ('// &
            trim(i8_str(max_bath))//' rows) - species "'//trim(name)//'" not registered')
      end if
      n_species = n_species + 1
      bath_tbl(n_species)%name = name
      bath_tbl(n_species)%sweep => sweep
   end subroutine bath_reg

   !------------------------------------------------------------------
   ! bath_setup() - provision from the BATH/N_BATH/DT_BATH input slots;
   !                called by the species' assembly seam AFTER its own init
   !                (the species parameters are loaded first, the loop shape
   !                is validated here). BATH = NONE stays inert
   !------------------------------------------------------------------
   subroutine bath_setup()
      integer :: row
      armed = .false.
      sweep_c => null()
      n_step_c = 0
      dt_c = 0.0d0
      setup_done = .true.
      if (reactants%bath /= 2) return      ! 1 = NONE (or 0 = unset): inert
      row = find_species('andersen')
      if (row == 0) then
         call stop_bath('bath_setup', 'the BATH provision selects the Andersen species '// &
            'but no bath species registered it (incomplete bath library assembly)')
      end if
      if (reactants%n_bath < 0) then
         call stop_bathi('bath_setup', 'bath step count n_bath must be >= 0, got', reactants%n_bath)
      end if
      if (reactants%n_bath > 0 .and. reactants%dt_bath <= 0.0d0) then
         call stop_bathz('bath_setup', 'a positive bath step size dt_bath [10 fs] is required '// &
            'when n_bath > 0, got', reactants%dt_bath)
      end if
      sweep_c => bath_tbl(row)%sweep
      n_step_c = reactants%n_bath
      dt_c = reactants%dt_bath
      armed = .true.
   end subroutine bath_setup

   !------------------------------------------------------------------
   ! bath_armed() - the provision probe (true = the bath loop will run)
   !------------------------------------------------------------------
   pure logical function bath_armed()
      bath_armed = armed .and. n_step_c > 0
   end function bath_armed

   !------------------------------------------------------------------
   ! bath_run(sta) - the bath loop: N_BATH steps of [propagation at DT_BATH
   !                 + one collision sweep], then ONE exact whole-system
   !                 drift removal. Unarmed or zero-step: bitwise no-op.
   !                 The global step and the trajectory time are swapped
   !                 for the run and restored after it (the thermalize
   !                 member's precedent)
   !------------------------------------------------------------------
   subroutine bath_run(sta)
      type(state_t), intent(inout) :: sta    ! initial state (bathed in place)
      integer :: i_step, i, i3, nat, nfrag
      integer :: fi, fa
      real(8) :: t_saved, dt_saved, wt, p_tot(3), d_use
      real(8), allocatable :: q_save(:), p_save(:)
      logical, allocatable :: in_frag(:), mask(:)
      if (.not. armed .or. n_step_c <= 0) return
      if (.not. associated(sweep_c)) then
         call stop_bath('bath_run', 'no species sweep bound (bath_setup must run first)')
      end if
      nat = size(sta%q)/3

      ! 1. precondition - seed the ownerless atoms (those no roster fragment
      !    claims: the slab of the surface paradigm) from the whole-system
      !    equilibrium geometry list_atoms%q0 with zero momenta. An armed
      !    bath implies this seeding: bathing atoms without a configuration
      !    is meaningless, and the geometry file data would otherwise sit
      !    unconsumed. Fragment atoms are NOT seeded - their members own
      !    their interiors
      if (.not. allocated(list_atoms%q0)) then
         call stop_bath('bath_run', 'no whole-system equilibrium geometry (list_atoms%'// &
            'q0 unallocated) - the atom list must be assembled first (assembly-order error)')
      end if
      if (size(list_atoms%q0) /= 3*nat) then
         call stop_bath('bath_run', 'the equilibrium geometry size disagrees with the state '// &
            '- state_create must follow the atom list')
      end if
      allocate (in_frag(nat))
      in_frag = .false.
      nfrag = 0
      if (allocated(list_atoms%frag)) nfrag = size(list_atoms%frag)
      do fi = 1, nfrag
         do fa = 1, list_atoms%frag(fi)%nat
            in_frag(list_atoms%frag(fi)%list(fa)) = .true.
         end do
      end do
      do i = 1, nat
         if (in_frag(i)) cycle
         i3 = 3*(i - 1)
         sta%q(i3+1:i3+3) = list_atoms%q0(i3+1:i3+3)
         sta%p(i3+1:i3+3) = 0.0d0
      end do

      ! 2. the freeze statute: every scheme-carrying fragment is
      !    member-owned - its sampled interior is saved bitwise and restored
      !    bitwise after the loop; the bath leaves no footprint on it. For
      !    the loop the fragment is PLACED at a PES-validated non-interacting
      !    separation (geom_sep_probe: ladder scan, interaction force under
      !    sep_tol - doubles as the cheap PES sanity check, a non-decaying
      !    interaction aborts naming the suspect); the incident channel
      !    re-places it by prescription afterwards anyway
      allocate (q_save(3*nat), p_save(3*nat))
      q_save = sta%q
      p_save = sta%p
      do fi = 1, nfrag
         call geom_sep_probe(probe_wrap, list_atoms%mass, list_atoms%frag(fi)%list, &
                             sta%q, reactants%r_sep, sep_tol, d_use)
      end do

      ! 3. the loop (collision scope = the ownerless atoms, by the freeze
      !    statute; scheme-carrying fragments are never swept)
      allocate (mask(nat))
      mask = .not. in_frag
      t_saved = sta%t
      dt_saved = prop_cfg%dt
      prop_cfg%dt = dt_c
      do i_step = 1, n_step_c
         call prop_step(sta)
         call sweep_c(sta, dt_c, mask)
      end do
      prop_cfg%dt = dt_saved
      sta%t = t_saved

      ! one exact whole-system COM-momentum projection: p_i <- p_i - (m_i/M)*P
      ! (the collisions inject a random net momentum; the incident channel
      ! prescribes the beam translation afterwards, so the residual goes)
      wt = sum(list_atoms%mass)
      p_tot = 0.0d0
      do i3 = 1, size(sta%p), 3
         p_tot = p_tot + sta%p(i3:i3+2)
      end do
      do i3 = 1, size(sta%p), 3
         sta%p(i3:i3+2) = sta%p(i3:i3+2) - (list_atoms%mass((i3+2)/3)/wt)*p_tot
      end do

      ! 5. unfreeze: bitwise restore of every scheme-carrying fragment ONLY
      !    (geometry AND momenta - exactly what its member sampled); the
      !    ownerless atoms keep their bathed q/p - that snapshot is their
      !    trajectory-starting configuration
      do i = 1, nat
         if (.not. in_frag(i)) cycle
         i3 = 3*(i - 1)
         sta%q(i3+1:i3+3) = q_save(i3+1:i3+3)
         sta%p(i3+1:i3+3) = p_save(i3+1:i3+3)
      end do
      deallocate (q_save, p_save, in_frag, mask)

   end subroutine bath_run

   ! probe_wrap(q, f) - the bound container through the interface probe
   ! (the same channel spectrum_interface rides; the unit fold stays at
   ! its unique site)
   subroutine probe_wrap(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      call container_probe_force(q, f)
   end subroutine probe_wrap

   ! find_species(name) - species-table row of name (0 = absent)
   function find_species(name) result(row)
      character(len=*), intent(in) :: name
      integer :: row, i
      row = 0
      do i = 1, n_species
         if (trim(bath_tbl(i)%name) == trim(name)) then
            row = i
            return
         end if
      end do
   end function find_species

   ! i8_str(v) - integer to decimal string for named abort messages
   function i8_str(v) result(s)
      integer, intent(in) :: v
      character(len=32) :: s
      write (s, '(i0)') v
   end function i8_str

   ! stop_bath(where, msg) - the named-abort channel
   subroutine stop_bath(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the bath loop is not usable)'
      stop 1
   end subroutine stop_bath

   subroutine stop_bathz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      write (0, '(a)') trim(where)//': fatal (the bath loop is not usable)'
      stop 1
   end subroutine stop_bathz

   subroutine stop_bathi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      write (0, '(a)') trim(where)//': fatal (the bath loop is not usable)'
      stop 1
   end subroutine stop_bathi

end module bath
