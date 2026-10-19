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
!   sweep [, reset]) registers a species (the optional reset is the stateful
!   species' per-trajectory zeroing, called at the loop head - 2026-10-15
!   legislation for the thermostat chains); bath_setup resolves the species
!   BY the BATH word (any word past NONE - 2026-10-15 generalization, was
!   hardcoded 'andersen'); BATH = NONE stays bitwise inert:
!   no RNG draws, no state writes - the archived stream is untouched);
!   bath_run(sta) executes N_BATH steps of [prop_step at DT_BATH + the
!   species sweep], saves/restores the global step and trajectory time
!   around it, and ends with ONE exact whole-system COM-momentum
!   projection (the collisions inject a random net momentum; the incident
!   channel prescribes the beam translation, so the residual is stripped).
!   Units: dt [fs]; p [amu*Ang/fs].
!=====================================================================
module bath
   use state,      only: state_t
   use config,     only: reactants, prop_cfg => propagator
   use config_atoms, only: list_atoms
   use propagator, only: prop_step
   use geom_init,  only: geom_sep_probe
   use config_atoms, only: substrate_mask, list_atoms
   use force_interface, only: container_probe_force
   use input,      only: input_bath_word
   use diag,       only: diag_fatal, diag_cat_sampler
   implicit none
   private
   public :: bath_reg, bath_setup, bath_run, bath_armed

   integer, parameter :: max_bath = 4     ! species-table capacity [row]
   real(8), parameter :: sep_tol = 1.0d-8  ! fragment-separation interaction
                                           ! tolerance [amu*Ang/fs^2] (~1.04e-6 eV/Angstrom;
                                           ! same physical value as the pre-2026-10-06
                                           ! 1e-6 in old internal-per-A units)

   ! Bath species interface: one collision sweep at step size dt over the
   ! masked atoms (the mask is the BATH_SCOPE resolution - all atoms, or the
   ! non-fragment atoms of the surface paradigm). A species registers EITHER
   ! a sweep (the loop propagates: [prop_step; sweep] per bath step - the
   ! Andersen form) OR a step (the species takes the bath-step integration
   ! over ENTIRELY - the extended-variable thermostat form: the chain must
   ! interleave with the kick/drift beats, which an opaque prop_step cannot
   ! offer, so the species carries its own splitting; ANT manual p.71: NH
   ! cannot ride a plain/position-Verlet/Beeman integrator, only a Liouville
   ! VV - the takeover legislates that requirement)
   abstract interface
      subroutine bath_sweep_i(sta, dt, mask)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! physical state (the masked atoms' momenta are
                                             ! reset on collision)
         real(8), intent(in) :: dt            ! the sweep's step size [fs] (the collision
                                             ! probability is nu*dt)
         logical, intent(in) :: mask(:)       ! per-atom collision scope [flag]
      end subroutine bath_sweep_i
      subroutine bath_reset_i()
      end subroutine bath_reset_i
      subroutine bath_seed_i(sta, mask)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! physical state (the masked atoms'
                                             ! momenta are Maxwell-seeded here)
         logical, intent(in) :: mask(:)       ! per-atom seed scope [flag]
      end subroutine bath_seed_i
   end interface

   ! Species table (the registry family - a new species = a new methods/
   ! file + one registration line)
   type :: bath_entry_t
      character(len=32) :: name = ''                       ! species name (the BATH word)
      procedure(bath_sweep_i), pointer, nopass :: sweep => null() ! one collision sweep
      procedure(bath_sweep_i), pointer, nopass :: step => null()  ! the self-contained bath
                                             ! step (takes the integration over)
      procedure(bath_reset_i), pointer, nopass :: reset => null() ! optional per-trajectory
                                             ! state reset (a stateful species - e.g. a
                                             ! thermostat chain - zeroes its carried
                                             ! variables; a stateless species omits it)
      procedure(bath_seed_i), pointer, nopass :: seed => null()   ! loop-head one-shot
                                             ! seeding (the thermalize form: Maxwell
                                             ! over the masked atoms; the loop body is
                                             ! then a bare run - no sweep, no takeover)
   end type bath_entry_t
   type(bath_entry_t) :: bath_tbl(max_bath)  ! species table (module-private)
   integer :: n_species = 0                  ! registered row count [row]

   ! Provisioned state (module-private; set by bath_setup)
   procedure(bath_sweep_i), pointer :: sweep_c => null()          ! the bound species sweep
   procedure(bath_sweep_i), pointer :: step_c => null()           ! the bound species step (takeover)
   procedure(bath_reset_i), pointer :: reset_c => null()          ! the bound species reset (nullable)
   procedure(bath_seed_i), pointer :: seed_c => null()            ! the bound species seed
                                             ! (the bare-run form's loop-head seeding)
   integer :: n_step_c = 0                    ! provisioned step count [-]
   real(8) :: dt_c = 0.0d0                    ! provisioned step size [fs]
   logical :: armed = .false.                 ! the BATH provision armed this bath [flag]
   logical :: setup_done = .false.            ! bath_setup ran [flag]
contains
   !------------------------------------------------------------------
   ! bath_reg(name [, sweep] [, step] [, reset]) - register a bath
   !                         species (the BATH enum word names it; duplicate
   !                         names and capacity abort; EXACTLY ONE of
   ! sweep/step must be present - the loop-propagated collision form or
   ! the self-contained takeover form; the optional reset is the
   ! per-trajectory state reset of a stateful species - method_reg's
   ! optional-sync shape)
   !------------------------------------------------------------------
   subroutine bath_reg(name, sweep, step, seed, reset)
      character(len=*), intent(in) :: name   ! species name (must equal a BATH word)
      procedure(bath_sweep_i), optional :: sweep  ! one collision sweep
      procedure(bath_sweep_i), optional :: step   ! the self-contained bath step
      procedure(bath_seed_i), optional :: seed    ! loop-head one-shot seeding (bare run)
      procedure(bath_reset_i), optional :: reset  ! per-trajectory state reset
      integer :: i, n_forms
      n_forms = count([present(sweep), present(step), present(seed)])
      if (n_forms /= 1) then
         call stop_bath('bath_reg', 'species "'//trim(name)//'" must register exactly '// &
            'one of sweep (loop-propagated collisions), step (self-contained takeover) '// &
            'or seed (one-shot seeding + bare run)')
      end if
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
      if (present(sweep)) then
         bath_tbl(n_species)%sweep => sweep
      else if (present(step)) then
         bath_tbl(n_species)%step => step
      else
         bath_tbl(n_species)%seed => seed
      end if
      if (present(reset)) then
         bath_tbl(n_species)%reset => reset
      else
         bath_tbl(n_species)%reset => null()
      end if
   end subroutine bath_reg

   !------------------------------------------------------------------
   ! bath_setup() - provision from the BATH/N_BATH/DT_BATH input slots;
   !                called by the species' assembly seam AFTER its own init
   !                (the species parameters are loaded first, the loop shape
   !                is validated here). BATH = NONE stays inert
   !------------------------------------------------------------------
   subroutine bath_setup()
      integer :: row
      character(len=32) :: word
      armed = .false.
      sweep_c => null()
      step_c => null()
      reset_c => null()
      n_step_c = 0
      dt_c = 0.0d0
      setup_done = .true.
      if (reactants%bath <= 1) return    ! 1 = NONE (or 0 = unset): inert
      ! the species is named by the BATH word itself (any word past NONE)
      word = input_bath_word(reactants%bath)
      if (len_trim(word) == 0) then
         call stop_bath('bath_setup', 'the BATH provision code '//trim(i8_str(reactants%bath))// &
            ' is outside the bath word table (grammar/table drift)')
      end if
      row = find_species(word)
      if (row == 0) then
         call stop_bath('bath_setup', 'the BATH provision selects species "'//trim(word)// &
            '" but no bath species registered it (incomplete bath library assembly)')
      end if
      if (reactants%n_bath < 0) then
         call stop_bathi('bath_setup', 'bath step count n_bath must be >= 0, got', reactants%n_bath)
      end if
      if (reactants%n_bath > 0 .and. reactants%dt_bath <= 0.0d0) then
         call stop_bathz('bath_setup', 'a positive bath step size dt_bath [fs] is required '// &
            'when n_bath > 0, got', reactants%dt_bath)
      end if
      sweep_c => bath_tbl(row)%sweep
      step_c => bath_tbl(row)%step
      reset_c => bath_tbl(row)%reset
      seed_c => bath_tbl(row)%seed
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
      logical, allocatable :: in_frag(:), mask(:), frozen(:)
      if (.not. armed .or. n_step_c <= 0) return
      if (.not. associated(sweep_c) .and. .not. associated(step_c) .and. &
          .not. associated(seed_c)) then
         call stop_bath('bath_run', 'no species sweep/step/seed bound (bath_setup must run first)')
      end if
      nat = size(sta%q)/3

      ! 1. precondition - place the substrate geometry: ownerless atoms
      !    (those no roster fragment claims) from list_atoms%q0 with zero
      !    momenta, and the surface TARGET fragment's atoms at their q0
      !    lattice rows (q only - the member-sampled momenta stand until the
      !    species' seed/sweep rewrites them; a NONE target's rows would
      !    otherwise sit at the state-create zeros, and this run's forces
      !    would see a collapsed geometry). An armed bath implies this
      !    placement: bathing atoms without a configuration is meaningless,
      !    and the geometry file data would otherwise sit unconsumed. Other
      !    fragment atoms are NOT placed - their members own their interiors
      !    (the incident channel places the projectile afterwards)
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
         if (.not. in_frag(i)) then
            i3 = 3*(i - 1)
            sta%q(i3+1:i3+3) = list_atoms%q0(i3+1:i3+3)
            sta%p(i3+1:i3+3) = 0.0d0
         end if
      end do
      if (nfrag > 0 .and. list_atoms%i_targ >= 1 .and. list_atoms%i_targ <= nfrag) then
         do fa = 1, list_atoms%frag(list_atoms%i_targ)%nat
            i = list_atoms%frag(list_atoms%i_targ)%list(fa)
            i3 = 3*(i - 1)
            sta%q(i3+1:i3+3) = list_atoms%q0(i3+1:i3+3)
         end do
      end if
      ! zero-state atoms (no member wrote their geometry - e.g. a NONE projectile
      ! before the incident placement): reset to the file baseline. The separation
      ! prober below measures forces at the current geometry, so a still-zero
      ! projectile would evaluate the PES on a collapsed site; a sampled atom
      ! (nonzero q) keeps its member's work - its placement authority already ran
      do i = 1, nat
         i3 = 3*(i - 1)
         if (all(sta%q(i3+1:i3+3) == 0.0d0) .and. &
             .not. all(list_atoms%q0(i3+1:i3+3) == 0.0d0)) then
            sta%q(i3+1:i3+3) = list_atoms%q0(i3+1:i3+3)
         end if
      end do

      ! 2. the freeze statute (2026-10-16 narrowed): every fragment EXCEPT the
      !    projectile (i_proj - re-placed by the incident channel afterwards
      !    anyway) and the surface TARGET (i_targ - the bath's own domain; its
      !    post-bath state IS the trajectory-starting substrate state) is
      !    member-owned - sampled interior saved bitwise, restored bitwise
      !    after the loop. Every placed fragment (including the projectile)
      !    still rides geom_sep_probe out to a PES-validated non-interacting
      !    separation for the run (ladder scan under sep_tol - doubles as the
      !    cheap PES sanity check; a non-decaying interaction aborts naming
      !    the suspect)
      allocate (q_save(3*nat), p_save(3*nat))
      q_save = sta%q
      p_save = sta%p
      do fi = 1, nfrag
         ! the TARGET fragment is the bath's OWN domain - separating it makes no
         ! sense and its post-bath state IS the product (2026-10-16: the probe
         ! loop rides the freeze statute - projectile out, everything else frozen)
         if (list_atoms%i_targ >= 1 .and. fi == list_atoms%i_targ) cycle
         if (reactants%surface_model /= 0) then
            ! a periodic surface: separate ALONG z (the vacuum) - the lateral
            ! landscape never decays along x (see geom_sep_probe's dir note)
            call geom_sep_probe(probe_wrap, list_atoms%mass, list_atoms%frag(fi)%list, &
                                sta%q, reactants%r_sep, sep_tol, d_use, [0.0d0, 0.0d0, 1.0d0])
         else
            call geom_sep_probe(probe_wrap, list_atoms%mass, list_atoms%frag(fi)%list, &
                                sta%q, reactants%r_sep, sep_tol, d_use)
         end if
      end do
      ! the frozen set: fragment atoms minus the projectile and TARGET rows
      allocate (frozen(nat))
      frozen = .false.
      do fi = 1, nat
         frozen(fi) = in_frag(fi)
      end do
      if (nfrag > 0) then
         if (list_atoms%i_proj >= 1 .and. list_atoms%i_proj <= nfrag) then
            do fa = 1, list_atoms%frag(list_atoms%i_proj)%nat
               frozen(list_atoms%frag(list_atoms%i_proj)%list(fa)) = .false.
            end do
         end if
         if (list_atoms%i_targ >= 1 .and. list_atoms%i_targ <= nfrag) then
            do fa = 1, list_atoms%frag(list_atoms%i_targ)%nat
               frozen(list_atoms%frag(list_atoms%i_targ)%list(fa)) = .false.
            end do
         end if
      end if

      ! 3. the loop. A SWEEP species rides the loop's propagation:
      !    [prop_step at DT_BATH + one sweep] per step (the collision scope
      !    = the ownerless atoms, by the freeze statute; scheme-carrying
      !    fragments are never swept). A STEP species TAKES the integration
      !    over entirely (the extended-variable thermostat form - the chain
      !    needs the kick/drift beats). A STATEFUL species (a thermostat
      !    chain) resets its carried variables HERE - once per trajectory,
      !    at the loop head before the first step (a stateless species
      !    registers no reset and skips this)
      ! bath domain (2026-10-16 ruling): the substrate domain - ownerless atoms
      ! UNION the surface TARGET fragment's mobile atoms (ghost rows excluded).
      ! With no TARGET and no ownerless atoms this is the historical empty set
      ! (every fragment atom froze); with one, the bath owns the substrate's
      ! thermal state while the projectile stays frozen out at a separated site
      allocate (mask(nat))
      mask = substrate_mask()
      if (associated(reset_c)) call reset_c()
      t_saved = sta%t
      dt_saved = prop_cfg%dt
      prop_cfg%dt = dt_c
      if (associated(seed_c)) call seed_c(sta, mask)
      if (associated(step_c)) then
         do i_step = 1, n_step_c
            call step_c(sta, dt_c, mask)
         end do
      else if (associated(sweep_c)) then
         do i_step = 1, n_step_c
            call prop_step(sta)
            call sweep_c(sta, dt_c, mask)
         end do
      else
         do i_step = 1, n_step_c              ! the bare run (seeded form): free
            call prop_step(sta)               ! equilibration, no in-loop coupling
         end do
      end if
      prop_cfg%dt = dt_saved
      sta%t = t_saved

      ! one exact COM-momentum projection over the BATH DOMAIN (2026-10-16:
      ! scoped - the historical whole-system form charged the frozen fragments
      ! with the domain's residual; the projectile stays out of the bath's
      ! business, the incident channel prescribes its translation afterwards):
      ! p_i <- p_i - (m_i/M_domain)*P_domain over the masked atoms
      wt = 0.0d0
      p_tot = 0.0d0
      i = 0
      do i3 = 1, size(sta%p), 3
         if (.not. mask((i3+2)/3)) cycle
         i = i + 1
         wt = wt + list_atoms%mass((i3+2)/3)
         p_tot = p_tot + sta%p(i3:i3+2)
      end do
      ! a SINGLE-atom domain has no drift to remove - its whole momentum IS its
      ! thermal motion (projecting would freeze it dead; the multi-atom domain
      ! sheds its rigid translation, keeping the internal thermal motion)
      if (wt > 0.0d0 .and. i >= 2) then
         do i3 = 1, size(sta%p), 3
            if (.not. mask((i3+2)/3)) cycle
            sta%p(i3:i3+2) = sta%p(i3:i3+2) - (list_atoms%mass((i3+2)/3)/wt)*p_tot
         end do
      end if

      ! 5. unfreeze: bitwise restore of the FROZEN set only (geometry AND
      !    momenta - exactly what its member sampled); the bath-domain atoms
      !    (ownerless + TARGET fragment) keep their bathed q/p - that snapshot
      !    is the trajectory-starting substrate state; the projectile is
      !    re-placed by the incident channel's prescription
      do i = 1, nat
         if (.not. frozen(i)) cycle
         i3 = 3*(i - 1)
         sta%q(i3+1:i3+3) = q_save(i3+1:i3+3)
         sta%p(i3+1:i3+3) = p_save(i3+1:i3+3)
      end do
      deallocate (q_save, p_save, in_frag, mask, frozen)

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
      call diag_fatal(diag_cat_sampler, 'bath.f90', trim(where), &
         'the bath loop is not usable')
   end subroutine stop_bath

   subroutine stop_bathz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      call diag_fatal(diag_cat_sampler, 'bath.f90', trim(where), &
         'the bath loop is not usable')
   end subroutine stop_bathz

   subroutine stop_bathi(where, msg, ival)
      character(len=*), intent(in) :: where, msg
      integer, intent(in) :: ival
      write (0, '(a,a,i0)') trim(where)//': ', trim(msg), ival
      call diag_fatal(diag_cat_sampler, 'bath.f90', trim(where), &
         'the bath loop is not usable')
   end subroutine stop_bathi

end module bath
