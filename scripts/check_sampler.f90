!=====================================================================
! check_sampler.f90 - sampler framework criteria program (criteria defined
!                     FIRST in this header; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files
! out of the source tree):
!   gfortran -Jbuild globals/control.f90 globals/config.f90 globals/consts.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 &
!            globals/config_atoms.f90 globals/state.f90 utils/rng.f90 &
!            utils/geometry.f90 utils/linalg.f90 utils/hessian.f90 &
!            utils/geom_init.f90 utils/spectrum.f90 interface/elec_interface.f90 &
!            interface/force_interface.f90 dynamics/propagator.f90 interface/bath.f90 &
!            interface/sampler.f90 &
!            ../scripts/check_sampler.f90 -o build/check_sampler.exe
! All artifacts are written under build/ - the source tree stays clean.
! Named-abort verification approach (same scheme as the other criteria
! programs): every named-abort criterion runs as a CHILD PROCESS - the check
! re-spawns itself with a mode argument (execute_command_line); the child
! drives the REAL sampler entry into the stop and must print a naming message
! on stderr and exit nonzero. The parent asserts the exit status and the
! stderr content; a child that returns without stopping exits 0 and fails
! the criterion.
! The distribution members are exercised through LOCAL stand-in procedures
! with the member interface, registered under the real member names - this
! program tests the FRAMEWORK (registration, dispatch order, incident
! channel, count reconcile); the member physics carries its own checks.
! Criteria:
!   (a) per-fragment dispatch order + member ownership: two-fragment list_atoms,
!       dist_scheme = the boltzmann,ebk codes, samp_run('incident', sta) calls
!       the fragment members in list_atoms order then the incident member once;
!       each stand-in writes ONLY its own slots (q runs 11/22, p 33) and f
!       keeps its pre-fill sentinel (nobody writes f); a second samp_run call
!       repeats the sequence (the orchestration is per call, incident once
!       each)
!   (b) incident-channel dispatch: a single-fragment gas-phase system
!       (surface_model = 0) dispatches the TRIVIAL single-reactant member
!       (a bitwise no-op - log carries its draw/realize pair, p untouched);
!       a single-fragment surface-paradigm system dispatches the writing
!       member (log grows by its tag, p written)
!   (c) count-reconcile named abort (child): dist_scheme length 1 against a
!       two-fragment list_atoms -> nonzero exit, stderr names samp_run and BOTH
!       values (length 1, fragment count 2)
!   (d) missing-scheme named abort (child): list_atoms assembled, dist_scheme
!       unallocated -> stderr names samp_run and dist_scheme
!   (e) unregistered scheme name named abort (child): the scheme argument names
!       a member that is not registered -> stderr names samp_run and the name
!   (f) unregistered fragment member named abort (child): a valid scheme code
!       whose member name is not registered -> stderr names samp_run, the
!       member name and the fragment
!   (g) out-of-range scheme code named abort (child): code 12 in fragment 2 ->
!       stderr names samp_run, the code and the fragment (sweep-along: the
!       probe moved off 7 when the scheme-code table grew to seven words,
!       then off 8 at ten words, then off 11 at eleven words)
!   (h) duplicate-registration named abort (child): samp_reg twice with the
!       same name -> stderr names samp_reg and the duplicate
!   (i) full-table named abort (child): 16 rows registered, the 17th -> stderr
!       names samp_reg and the capacity
!   (j) list_atoms assembly-order named abort (child): no list_atoms fragment table ->
!       stderr names samp_run and list_atoms_load
!   (k) state-size named abort (child): a state sized for 1 atom against a
!       4-atom list_atoms -> stderr names samp_run and BOTH sizes
!   (l) per-member dispatch dedup: DIST_SCHEME carrying the SAME word on two
!       fragments dispatches that member ONCE per samp_run call (the member
!       serves both carriers inside the one call - the multi-carrier stand-in
!       writes both fragments' q runs); distinct members keep list_atoms
!       first-occurrence order
!   (m) orientation-stage neutrality (option off): random_orient = .false.
!       (the default and the archived-stream contract), scheme [3,4] - the
!       members' q writes pass through bitwise (runs 44/55 and 66/77/88
!       untouched, the unwritten atom keeps its sentinel), the incident
!       channel still writes p = 33 - the stage consumes nothing and touches
!       nothing
!   (n) orientation-stage semantics (option on): same list_atoms, random_orient
!       = .true. with a fixed seed - the 'rotation'-word fragment (frag1) is
!       SKIPPED (its q runs stay exactly 44/55); the non-rotation fragment
!       (frag2) is rotated: its unwritten-sentinel atom leaves the sentinel
!       (a rotation moved it), the pair distance |q3-q4| is preserved, the
!       fragment COM re-centers at the coordinate origin (the geometry
!       instrument's write-back contract), and the rotated q differs from
!       the unrotated write
!   (o) orientation-stage nat guard: a single-fragment single-atom gas
!       list_atoms with the option on - the monoatomic fragment dispatches
!       no member and draws no orientation matrix (q/p/f bitwise sentinel;
!       only the no-op single member's dispatch pair logs)
!   (p) paradigm word resolution: samp_incident_word returns incident_single
!       / incident_pair / incident_surface for (single-fragment gas, gas
!       pair, surface paradigm)
!   (q) substrate orientation exemption: surface paradigm, option on - only
!       the projectile fragment is orientation-eligible (the substrate
!       carrier's q run stays bitwise)
!=====================================================================
module sampler_standin
   use state, only: state_t
   implicit none
   private
   public :: stand_frag1_d, stand_frag1_r, stand_frag2_d, stand_frag2_r, &
             stand_incident_d, stand_incident_r, stand_both_d, stand_both_r, &
             stand_frag4_d, stand_frag4_r, stand_atom_d, stand_atom_r, &
             stand_single_d, stand_single_r, log_reset, log_is
   integer, parameter :: max_log = 32
   integer :: log_seq(max_log) = 0
   integer :: n_log = 0
contains
   ! call-order log: one tag per member invocation (1/2 = the fragment
   ! stand-ins, 3 = the incident stand-in)
   subroutine log_reset()
      n_log = 0
      log_seq = 0
   end subroutine log_reset

   logical function log_is(seq)
      integer, intent(in) :: seq(:)
      integer :: i
      log_is = (n_log == size(seq))
      if (log_is) then
         do i = 1, n_log
            if (log_seq(i) /= seq(i)) then
               log_is = .false.
               exit
            end if
         end do
      end if
   end function log_is

   subroutine push(tag)
      integer, intent(in) :: tag
      n_log = n_log + 1
      if (n_log <= max_log) log_seq(n_log) = tag
   end subroutine push

   ! stand-in members as draw/realize pairs: each REALIZE writes ONLY the
   ! slots it owns (a fragment member owns one q run, the incident member owns
   ! all of p); the f array is owned by nobody and must keep its sentinel
   ! fill. Both halves push the SAME tag, so the logged sequence shows the
   ! two-phase dispatch: all draws first (list_atoms order), then the
   ! realizations in the same order, then the incident draw+realize pair
   subroutine stand_frag1_d()
      call push(1)
   end subroutine stand_frag1_d

   subroutine stand_frag1_r(sta)
      type(state_t), intent(inout) :: sta
      call push(1)
      sta%q(1:3) = 11.0d0
   end subroutine stand_frag1_r

   subroutine stand_frag2_d()
      call push(2)
   end subroutine stand_frag2_d

   subroutine stand_frag2_r(sta)
      type(state_t), intent(inout) :: sta
      call push(2)
      sta%q(4:6) = 22.0d0
   end subroutine stand_frag2_r

   subroutine stand_incident_d()
      call push(3)
   end subroutine stand_incident_d

   subroutine stand_incident_r(sta)
      type(state_t), intent(inout) :: sta
      call push(3)
      sta%p = 33.0d0
   end subroutine stand_incident_r

   ! a multi-carrier-shaped stand-in: ONE realize call writes BOTH fragments'
   ! q runs (the member form a same-word double selection expects)
   subroutine stand_both_d()
      call push(4)
   end subroutine stand_both_d

   subroutine stand_both_r(sta)
      type(state_t), intent(inout) :: sta
      call push(4)
      sta%q(1:3) = 44.0d0
      sta%q(4:6) = 55.0d0
   end subroutine stand_both_r

   ! a NON-UNIFORM writer for fragment 2's first atom (the orientation stage's
   ! rigid-rotation observable: a uniform run would hide the rotation)
   subroutine stand_frag4_d()
      call push(5)
   end subroutine stand_frag4_d

   subroutine stand_frag4_r(sta)
      type(state_t), intent(inout) :: sta
      call push(5)
      sta%q(7:9) = [ 66.0d0, 77.0d0, 88.0d0 ]
   end subroutine stand_frag4_r

   ! a monoatomic-fragment writer: q AND p non-uniform (the nat guard's
   ! observable - a rotation would move both)
   subroutine stand_atom_d()
      call push(6)
   end subroutine stand_atom_d

   subroutine stand_atom_r(sta)
      type(state_t), intent(inout) :: sta
      call push(6)
      sta%q(1:3) = [ 1.0d0, 2.0d0, 3.0d0 ]
      sta%p(1:3) = [ 4.0d0, 5.0d0, 6.0d0 ]
   end subroutine stand_atom_r

   ! the trivial single-reactant stand-in (the inc_single shape): dispatched
   ! but owns NOTHING - draws no RNG, writes no state
   subroutine stand_single_d()
      call push(7)
   end subroutine stand_single_d

   subroutine stand_single_r(sta)
      type(state_t), intent(inout) :: sta
      call push(7)
   end subroutine stand_single_r
end module sampler_standin

program check_sampler
   use sampler, only: samp_reg, samp_run, samp_incident_word
   use sampler_standin
   use config_atoms, only: list_atoms, list_atoms_clear
   use config, only: reactants
   use state, only: state_t, state_create, state_destroy
   use rng, only: rng_init
   implicit none
   type(state_t) :: sta
   integer :: nfail = 0, nfail0, i
   logical :: ok_gas, ok_surf
   character(len=256) :: mode, name
   real(8), parameter :: sentinel = -777.0d0
   real(8) :: com2(3), dist_pre, dist_post

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('bad-count')    ! (c): length 1 vs 2 fragments
         call reg_base()
         call setup_mol(2)
         call set_scheme([1])
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-noscheme') ! (d): dist_scheme unallocated
         call reg_base()
         call setup_mol(2)
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-scheme')   ! (e): 'incident' not registered
         call samp_reg('boltzmann', stand_frag1_d, stand_frag1_r)
         call samp_reg('ebk', stand_frag2_d, stand_frag2_r)
         call setup_mol(2)
         call set_scheme([1, 2])
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-member')   ! (f): code 2 = 'ebk' not registered
         call samp_reg('boltzmann', stand_frag1_d, stand_frag1_r)
         call samp_reg('incident', stand_incident_d, stand_incident_r)
         call setup_mol(2)
         call set_scheme([1, 2])
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-code')     ! (g): code 11 in fragment 2 (one past the ten-word table)
         call reg_base()
         call setup_mol(2)
         call set_scheme([1, 11])
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-dup')      ! (h): same name twice
         call samp_reg('boltzmann', stand_frag1_d, stand_frag1_r)
         call samp_reg('boltzmann', stand_frag2_d, stand_frag2_r)
      case ('bad-full')     ! (i): 16 rows then the 17th
         do i = 1, 16
            write (name, '(a,i0)') 'stand', i
            call samp_reg(trim(name), stand_frag1_d, stand_frag1_r)
         end do
         call samp_reg('overflow', stand_frag2_d, stand_frag2_r)
      case ('bad-order')    ! (j): list_atoms%frag left unallocated
         call reg_base()
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-state')    ! (k): 3-dof state vs a 4-atom list_atoms
         call reg_base()
         call setup_mol(2)
         call set_scheme([1, 2])
         call state_create(sta, 1)
         call samp_run('incident', sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the sampler returned - bug)'
      stop 0
   end if

   ! ---------- stand-in registration (the framework's member table) ----------
   call reg_base()

   ! ---------- (a) per-fragment dispatch order + member ownership ----------
   nfail0 = nfail
   call setup_mol(2)
   call set_scheme([1, 2])
   call state_create(sta, 4)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call log_reset()
   call samp_run('incident', sta)
   call samp_run('incident', sta)   ! second trajectory: sequence repeats
   if (.not. log_is([1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]) .or. &
       any(sta%q(1:3) /= 11.0d0) .or. any(sta%q(4:6) /= 22.0d0) .or. &
       any(sta%p /= 33.0d0) .or. any(sta%f /= sentinel)) then
      print '(a)', 'FAIL (a) per-fragment dispatch: order/ownership mismatch'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(a) per-fragment dispatch: PASS (order boltzmann,ebk,incident per call; q runs 11/22, p 33, f sentinel)'
   call state_destroy(sta)

   ! ---------- (b) incident-channel dispatch: the trivial single-reactant
   ! member is dispatched UNIFORMLY (a bitwise no-op - the historical
   ! skip-guard produced the identical stream); the surface member writes ----------
   nfail0 = nfail
   call setup_mol(1)
   call set_scheme([1])
   call state_create(sta, 2)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call log_reset()
   call samp_run('incident_single', sta)   ! single-fragment gas: no-op dispatch
   ok_gas = log_is([1, 7, 1, 7]) .and. all(sta%q(1:3) == 11.0d0) .and. all(sta%q(4:6) == sentinel) .and. &
            all(sta%p == sentinel)
   reactants%surface_model = 1
   sta%q = sentinel
   sta%p = sentinel
   call log_reset()
   call samp_run('incident_surface', sta)  ! single-fragment surface: channel on
   ok_surf = log_is([1, 3, 1, 3]) .and. all(sta%p == 33.0d0)
   reactants%surface_model = 0
   if (.not. ok_gas .or. .not. ok_surf) then
      print '(a)', 'FAIL (b) incident dispatch: the no-op single dispatch or the surface dispatch wrong'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(b) incident dispatch: PASS (single-fragment gas dispatches the no-op single member '// &
      'bitwise untouched; the surface paradigm dispatches the writing member)'
   call state_destroy(sta)

   ! ---------- (p) paradigm word resolution: samp_incident_word returns the
   ! member the regs arm for the derived paradigm ----------
   nfail0 = nfail
   call setup_mol(1)
   ok_gas = (trim(samp_incident_word()) == 'incident_single')
   call setup_mol(2)
   ok_surf = (trim(samp_incident_word()) == 'incident_pair')
   reactants%surface_model = 1
   if (trim(samp_incident_word()) /= 'incident_surface') then
      print '(a)', 'FAIL (p) paradigm word: the surface resolution wrong'
      nfail = nfail + 1
   end if
   reactants%surface_model = 0
   if (.not. ok_gas .or. .not. ok_surf) then
      print '(a)', 'FAIL (p) paradigm word: the molecular resolutions wrong'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(p) paradigm word resolution: PASS (single-fragment gas -> incident_single, '// &
      'gas pair -> incident_pair, surface -> incident_surface)'

   ! ---------- (l) per-member dispatch dedup: a same-word double selection
   ! is ONE call (the member serves both carriers inside it) ----------
   nfail0 = nfail
   call setup_mol(2)
   call set_scheme([3, 3])              ! the rotation word on BOTH fragments
   call state_create(sta, 4)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call log_reset()
   call samp_run('incident', sta)
   call samp_run('incident', sta)       ! second trajectory: one call each again
   if (.not. log_is([4, 3, 4, 3, 4, 3, 4, 3]) .or. &
       any(sta%q(1:3) /= 44.0d0) .or. any(sta%q(4:6) /= 55.0d0) .or. &
       any(sta%p /= 33.0d0) .or. any(sta%f /= sentinel)) then
      print '(a)', 'FAIL (l) same-word dedup: call count or carrier coverage wrong'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(l) per-member dispatch dedup: PASS (same word on two fragments = ONE call serving '// &
      'both carriers; the call order and incident channel unchanged)'
   call state_destroy(sta)

   ! ---------- (m) orientation-stage neutrality: option OFF leaves the
   ! members' writes bitwise untouched (the archived-stream contract) ----------
   nfail0 = nfail
   call setup_mol(2)
   call set_scheme([3, 4])              ! frag1 carries the rotation WORD, frag2 'j'
   call state_create(sta, 4)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   reactants%random_orient = .false.
   call log_reset()
   call samp_run('incident', sta)
   if (.not. log_is([4, 5, 3, 4, 5, 3]) .or. &
       any(sta%q(1:3) /= 44.0d0) .or. any(sta%q(4:6) /= 55.0d0) .or. &
       any(sta%q(7:9) /= [ 66.0d0, 77.0d0, 88.0d0 ]) .or. &
       any(sta%q(10:12) /= sentinel) .or. &
       any(sta%p /= 33.0d0) .or. any(sta%f /= sentinel)) then
      print '(a)', 'FAIL (m) option-off neutrality: a member write moved or the log wrong'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(m) option-off neutrality: PASS (members q runs 44/55 and 66/77/88 bitwise, '// &
      'sentinel atom untouched, p = 33, f sentinel)'
   call state_destroy(sta)

   ! ---------- (n) orientation-stage semantics: rotation-word fragment skipped,
   ! non-rotation fragment rotated rigidly and re-centered ----------
   nfail0 = nfail
   dist_pre = norm2([ 66.0d0 + 777.0d0, 77.0d0 + 777.0d0, 88.0d0 + 777.0d0 ])
   call state_create(sta, 4)
   call rng_init(20260922)              ! the rotation draw is deterministic
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   reactants%random_orient = .true.
   call log_reset()
   call samp_run('incident', sta)
   com2 = 0.5d0*(sta%q(7:9) + sta%q(10:12))
   dist_post = norm2(sta%q(7:9) - sta%q(10:12))
   if (.not. log_is([4, 5, 3, 4, 5, 3]) .or. &
       any(sta%q(1:3) /= 44.0d0) .or. any(sta%q(4:6) /= 55.0d0) .or. &
       all(sta%q(7:9) == [ 66.0d0, 77.0d0, 88.0d0 ]) .or. &
       all(sta%q(10:12) == sentinel) .or. &
       norm2(com2) > 1.0d-9 .or. &
       abs(dist_post - dist_pre) > 1.0d-9*dist_pre .or. &
       any(sta%f /= sentinel)) then
      print '(a)', 'FAIL (n) option-on semantics: skip rule, rigidity or re-centering wrong'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(n) option-on semantics: PASS (rotation-word fragment skipped bitwise; the '// &
      'other rotated rigidly - pair distance preserved, COM re-centered at the origin)'
   reactants%random_orient = .false.
   call state_destroy(sta)

   ! ---------- (o) orientation-stage nat guard: a monoatomic fragment is
   ! skipped exactly (q/p bitwise the member writes) ----------
   nfail0 = nfail
   call list_atoms_clear()
   allocate (list_atoms%mass(1))
   list_atoms%mass = 1.0d0
   allocate (list_atoms%frag(1))
   list_atoms%frag(1)%nat = 1
   list_atoms%frag(1)%mass = 1.0d0
   allocate (list_atoms%frag(1)%list(1), list_atoms%frag(1)%qz(3))
   list_atoms%frag(1)%list = 1
   list_atoms%frag(1)%qz = 0.0d0
   reactants%surface_model = 0
   call set_scheme([9])                 ! the 'thermalize' word -> the monoatomic stand-in
   call state_create(sta, 1)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   reactants%random_orient = .true.
   call log_reset()
   call samp_run('incident_single', sta)  ! single-fragment gas: no-op dispatch
   if (.not. log_is([7, 7]) .or. &
       any(sta%q(1:3) /= sentinel) .or. &
       any(sta%p(1:3) /= sentinel) .or. &
       any(sta%f /= sentinel)) then
      print '(a)', 'FAIL (o) nat guard: the monoatomic fragment was touched'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(o) nat guard: PASS (the monoatomic fragment dispatches no member - '// &
      'q/p/f bitwise sentinel; only the no-op single pair logs)'
   reactants%random_orient = .false.
   call state_destroy(sta)

   ! ---------- (q) substrate orientation exemption: in the surface paradigm
   ! ONLY the projectile fragment is orientation-eligible (a substrate
   ! fragment carries the crystal frame, which must not rotate) ----------
   nfail0 = nfail
   call setup_mol(2)
   list_atoms%i_proj = 1
   call set_scheme([1, 2])              ! non-rotation words on both fragments
   call state_create(sta, 4)
   sta%q = sentinel
   sta%p = sentinel
   reactants%surface_model = 1
   reactants%random_orient = .true.
   call log_reset()
   call samp_run('incident_surface', sta)
   ! frag1 (the projectile) was rotated -> its q run left 11.0; frag2 (the
   ! substrate carrier) was NOT rotated -> its q run stays exactly 22.0
   ! frag1 (the projectile, atoms 1-2) was rotated -> its members' q writes
   ! (11 on atom 1, 22 on atom 2) both left their values; frag2 (the substrate
   ! carrier, atoms 3-4) was NOT rotated -> its atoms keep the sentinel
   if (any(sta%q(1:3) == 11.0d0) .or. any(sta%q(4:6) == 22.0d0) .or. &
       any(sta%q(7:12) /= sentinel)) then
      print '(a)', 'FAIL (q) substrate exemption: the substrate rotated or the projectile did not'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(q) substrate orientation exemption: PASS (surface paradigm: the projectile '// &
      'fragment rotates, the substrate carrier stays bitwise)'
   reactants%random_orient = .false.
   reactants%surface_model = 0
   call state_destroy(sta)

   ! ---------- (c)-(k) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-count',    'samp_run', 'length 1',         'count 2', &
                        '(c) count-reconcile loud-stop')
   call expect_loudstop('bad-noscheme', 'samp_run', 'dist_scheme',      '', &
                        '(d) missing-scheme loud-stop')
   call expect_loudstop('bad-scheme',   'samp_run', 'incident',         'not a registered', &
                        '(e) unregistered-scheme loud-stop')
   call expect_loudstop('bad-member',   'samp_run', 'ebk',              'fragment 2', &
                        '(f) unregistered-member loud-stop')
   call expect_loudstop('bad-code',     'samp_run', 'code 11',          'fragment 2', &
                        '(g) out-of-range-code loud-stop')
   call expect_loudstop('bad-dup',      'samp_reg', 'boltzmann',        'already registered', &
                        '(h) duplicate-registration loud-stop')
   call expect_loudstop('bad-full',     'samp_reg', '16 rows',          'overflow', &
                        '(i) full-table loud-stop')
   call expect_loudstop('bad-order',    'samp_run', 'list_atoms_load',      'assembly-order', &
                        '(j) list_atoms-order loud-stop')
   call expect_loudstop('bad-state',    'samp_run', '3 dof',            '12 dof', &
                        '(k) state-size loud-stop')

   print '(a,i0,a)', '=== check_sampler: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! register the stand-ins under the real member names
   subroutine reg_base()
      call samp_reg('boltzmann', stand_frag1_d, stand_frag1_r)
      call samp_reg('ebk', stand_frag2_d, stand_frag2_r)
      call samp_reg('incident', stand_incident_d, stand_incident_r)
      call samp_reg('incident_surface', stand_incident_d, stand_incident_r)
      call samp_reg('incident_single', stand_single_d, stand_single_r)
      call samp_reg('rotation', stand_both_d, stand_both_r)
      call samp_reg('j', stand_frag4_d, stand_frag4_r)
      call samp_reg('thermalize', stand_atom_d, stand_atom_r)
   end subroutine reg_base

   ! mini-list_atoms built by direct field fill (the atom list public surface allows
   ! it; the sysdef->list_atoms_load chain has its own criteria program): two
   ! atoms per fragment, molecular paradigm, no scheme preset
   subroutine setup_mol(nf)
      integer, intent(in) :: nf
      integer :: i
      call list_atoms_clear()
      allocate (list_atoms%mass(2*nf))
      list_atoms%mass = 1.0d0
      allocate (list_atoms%frag(nf))
      do i = 1, nf
         list_atoms%frag(i)%nat = 2
         list_atoms%frag(i)%mass = 2.0d0
         allocate (list_atoms%frag(i)%list(2), list_atoms%frag(i)%qz(6))
         list_atoms%frag(i)%list = [ 2*i - 1, 2*i ]
         list_atoms%frag(i)%qz = 0.0d0
      end do
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine setup_mol

   ! preset the per-fragment scheme codes
   subroutine set_scheme(codes)
      integer, intent(in) :: codes(:)
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
      allocate (reactants%dist_scheme(size(codes)))
      reactants%dist_scheme = codes
   end subroutine set_scheme

   ! spawn the child for one bad case and assert the named abort: nonzero exit,
   ! stderr names the expected authority (samp_reg / samp_run) and carries the
   ! expected key(s); key2 = '' skips the second key
   subroutine expect_loudstop(cmode, who, key1, key2, label)
      character(len=*), intent(in) :: cmode, who, key1, key2, label
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      ! the shell redirect needs native separators on Windows (forward slashes
      ! make cmd misparse the redirect target)
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_sampler.exe '//trim(cmode)//' 2>build'//sep//'sampler_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/sampler_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) line
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(line)
         end do
         close (v)
      end if
      if (stat == 0) then
         print '(a)', 'FAIL '//label//': expected a loud-stop, the child exited 0'
         nfail = nfail + 1
      else if (index(buf, who) == 0) then
         print '(a)', 'FAIL '//label//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (index(buf, trim(key1)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key1)//'"'
         nfail = nfail + 1
      else if (len_trim(key2) > 0 .and. index(buf, trim(key2)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key2)//'"'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop

   ! Windows-host probe (same runtime OS branch the sibling checks use)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_sampler
