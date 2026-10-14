!=====================================================================
! check_samp_thermalize.f90 - thermalize distribution-member criteria
!                    program (criteria defined FIRST in this header;
!                    this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree; the evolution chain enters because the
! member drives the real propagator through its equilibration run):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 &
!            globals/config_atoms.f90 globals/state.f90 &
!            utils/rng.f90 interface/elec_interface.f90 &
!            utils/linalg.f90 interface/force_interface.f90 dynamics/propagator.f90 &
!            utils/hessian.f90 utils/geom_init.f90 utils/spectrum.f90 interface/bath.f90 &
!            methods/reg_prop.f90 utils/geometry.f90 interface/sampler.f90 &
!            methods/samp_thermalize.f90 methods/seam_cast.f90 &
!            ../scripts/check_samp_thermalize.f90 -o build/check_samp_thermalize.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter rule: the member is driven by directly constructed
! parameters and a directly filled list_atoms; the input grammar is not in the
! loop. This is a MECHANISM-level check: no trajectory statistics are
! asserted (the equilibration RUN is exercised through the real verlet
! row on a check-bound ZERO-force container - free flight is the lawful
! mechanism-level posture, real potentials are posterior; the free-flight
! landing spot is asserted exactly); the kinetic-energy expectation and
! the drift identity are verified ALGEBRAICALLY from the member-public
! parameters and the exact transformation the member applies.
! The incident channel is not this check's subject: a LOCAL no-op stand-in
! is registered under the name 'incident' only because samp_run resolves
! the scheme argument before dispatch (the channel itself is skipped for
! the single-fragment gas-phase list_atoms used here).
! Named-abort verification approach (same scheme as the sibling checks):
! every named-abort criterion runs as a CHILD PROCESS - the check re-spawns
! itself with a mode argument (execute_command_line); the child drives the
! REAL member entry into the stop and must print a naming message on
! stderr and exit nonzero. The parent asserts the exit status and the
! stderr content.
! Criteria:
!   (a) registration + dispatch chain: the member registers through its
!       reg file row, samp_code('thermalize') = 10; a single-fragment
!       gas-phase list_atoms with dist_scheme = [10] and samp_run('incident',
!       sta) REACHES the member - the owned p components come off their
!       sentinels (the Maxwell seed) while q lands exactly on the
!       free-flight mimic and f exactly 0 (the zero-force container); the
!       incident channel is not dispatched
!   (b) post-equilibration kinetic-energy expectation algebra (analytic,
!       no sampling): the declared expectation obeys
!       e_kin_expect = (3*nat - 3)*kb*t_eq/2 for the bound carrier (the
!       mass-weighted COM-velocity Gaussian projection onto zero net
!       momentum removes exactly three degrees of freedom - independent
!       of the mass distribution), relative 1e-12
!   (c) drift identity (exact algebra of the transformation): after a
!       direct sample call the owned-set momentum sums to zero PER
!       COMPONENT (the mass-weighted projection p'_i = p_i - (m_i/M)*P
!       carries the total to zero exactly) - relative residual <= 1e-10
!       (floating-point roundoff of the subtraction, not a statistical
!       claim)
!   (d) projection-form replay (exact, no distribution assertion): the
!       seeded Gaussians replayed here in the member's draw order and
!       projected by the DECLARED mass-weighted formula reproduce the
!       member's output component-for-component, relative 1e-12 - an
!       unweighted mean subtraction would shift every component of the
!       unequal-mass carrier (1/4 amu), so the projection form itself
!       is pinned, not just the getter constant
!   (e) parameter-loading integrity: the t_eq/n_eq/dt_eq getters return
!       the values init received
!   (f) footprint confinement: with the propagator step preset to a
!       distinctive value and sta%t preset to a sentinel, a sample call
!       (which swaps the step to dt_eq for the internal run) RESTORES
!       both - the equilibration run neither consumes trajectory time nor
!       leaves the global step changed; on the zero-force container q
!       lands exactly on the free-flight mimic (p constant through the
!       run) and f exactly 0
!   (g) non-positive temperature named abort (child): t_eq = 0 -> stderr
!       names thermalize_init and t_eq
!   (h) non-positive step count named abort (child): n_eq = 0 -> stderr
!       names thermalize_init and n_eq
!   (i) non-positive step size named abort (child): dt_eq = 0 -> stderr
!       names thermalize_init and dt_eq
!   (j) no-fragment-binding named abort (child): no dist_scheme row carries
!       the code -> stderr names thermalize_init
!   (k) dual-carrier equal treatment (positive): two fragments both carry
!       the code (scheme [10,10]) -> one init binds BOTH carriers
!       (frag/nat getters, carrier 2 explicit) and e_kin_expect = the SUM
!       of the two (3*nat-3)*kb*t_eq/2 terms (1e-12); one sample seeds
!       every carrier (all owned p off sentinel), drift-removes each
!       carrier (per-component sums zero on both, 1e-10), then the ONE
!       shared equilibration run lands both fragments on the free-flight
!       mimic with step/time restored (retired with multi-carrier
!       support: the duplicate-binding abort - a same-word double
!       selection is the lawful form now)
!   (l) sample-before-init named abort (child): thermalize_sample with no
!       thermalize_init -> stderr names thermalize_draw (the two-phase guard site)
!   (m) list_atoms-offset short-state named abort (child): a 6-dof state
!       against a SECOND-fragment carrier (2x2-atom list_atoms, scheme
!       [1,10]: the owned components are 7..12) -> stderr names
!       thermalize_sample and the state (the guard must see the
!       fragment's list_atoms offset, not just the atom count)
!=====================================================================
module equi_standin
   use state, only: state_t
   implicit none
   private
   public :: stand_incident_d, stand_incident_r
contains
   ! the incident-channel stand-in: registered so samp_run's scheme-argument
   ! lookup resolves; NOT dispatched in these scenarios (single-fragment gas
   ! phase) - this check's subject is the thermalize member
   subroutine stand_incident_d()
   end subroutine stand_incident_d

   subroutine stand_incident_r(sta)
      type(state_t), intent(inout) :: sta
   end subroutine stand_incident_r
end module equi_standin

program check_samp_thermalize
   use consts,        only: kb_act, kb_ev
   use rng,           only: rng_init, rng_gauss
   use sampler,       only: samp_run, samp_code
   use config,        only: reactants, propagator
   use config_atoms, only: list_atoms, list_atoms_clear
   use state,         only: state_t, state_create, state_destroy
   use elec_interface,   only: method_reg
   use force_interface,  only: container_bind_pes, container_init
   use reg_prop,      only: reg_prop_all
   use samp_thermalize, only: thermalize_params_t, thermalize_init, thermalize_sample, &
                              thermalize_t_eq, thermalize_n_eq, thermalize_dt_eq, &
                              thermalize_e_kin_expect, thermalize_frag, thermalize_nat
   use samp_thermalize, only: thermalize_register => member_register
   use equi_standin, only: stand_incident_d, stand_incident_r
   implicit none
   type(state_t) :: sta
   type(thermalize_params_t) :: par
   integer :: nfail = 0, i, c, frag_b, istep
   real(8) :: e_expect, comp_sum, comp_abs, rel, rel_max
   real(8) :: exp_p(6), qe(12), wt, ptot
   real(8), parameter :: t_use = 325.0d0
   real(8), parameter :: dt_use = 0.0413d0
   real(8), parameter :: dt_live = 0.031d0
   integer, parameter :: n_use = 3
   real(8), parameter :: sentinel = -777.0d0
   character(len=32) :: mode

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   allocate (par%t_eq(1), par%n_eq(1), par%dt_eq(1))   ! single-carrier check face
   par%t_eq(1) = 0.0d0
   par%n_eq(1) = 0
   par%dt_eq(1) = 0.0d0
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 2)
      par%t_eq(1) = t_use
      par%n_eq(1) = n_use
      par%dt_eq(1) = dt_use
      propagator%dt = dt_live
      select case (trim(mode))
      case ('bad-temp')                     ! (g): t_eq = 0
         call setup_mol(1)
         call set_scheme([9])
         par%t_eq(1) = 0.0d0
         call thermalize_init(par)
      case ('bad-nsteps')                   ! (h): zero propagation steps
         call setup_mol(1)
         call set_scheme([9])
         par%n_eq(1) = 0
         call thermalize_init(par)
      case ('bad-dt')                       ! (i): zero equilibration step
         call setup_mol(1)
         call set_scheme([9])
         par%dt_eq(1) = 0.0d0
         call thermalize_init(par)
      case ('bad-bind0')                    ! (j): no fragment carries the code
         call setup_mol(1)
         call set_scheme([1])
         call thermalize_init(par)
      case ('bad-noinit')                   ! (l): sample without init
         call setup_mol(1)
         call set_scheme([9])
         call thermalize_sample(sta)
      case ('bad-offset')                   ! (m): 6-dof state vs a SECOND-fragment carrier
         call setup_mol(2)
         call set_scheme([1, 9])
         call rng_init(7737)                ! the draw half runs before the realize guard -
                                            ! the child must seed (per-criterion reseeding rule)
         call thermalize_init(par)
         call thermalize_sample(sta)        ! sta holds 2 atoms (6 dof) - owned 7..12
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- single carrier: algebra + registration/dispatch ----------
   call setup_mol(1)
   call set_scheme([9])
   par%t_eq(1) = t_use
   par%n_eq(1) = n_use
   par%dt_eq(1) = dt_use
   propagator%dt = dt_live
   call rng_init(7717)
   call thermalize_init(par)

   ! the evolution chain the member drives (real since the propagator
   ! landed): degenerate method row + check-local ZERO-force container +
   ! the verlet registration row - the equilibration run inside every
   ! sample below is genuine free-flight propagation, assembled CHECK-SIDE
   ! (the tree's default binding posture stays untouched)
   call method_reg('adiabatic', noop_prop, noop_match)
   call container_bind_pes(zero_pes)
   call container_init()
   call reg_prop_all()
   propagator%id = 1
   propagator%order = 1

   ! (e) parameter-loading integrity
   if (thermalize_t_eq() /= t_use .or. thermalize_n_eq() /= n_use .or. &
       thermalize_dt_eq() /= dt_use .or. thermalize_frag() /= 1 .or. thermalize_nat() /= 2) then
      print '(a)', 'FAIL (e) parameter integrity: an equilibration parameter or the carrier was dropped'
      nfail = nfail + 1
   else
      print '(a)', '(e) parameter integrity: PASS (t_eq/n_eq/dt_eq round-trip; carrier fragment 1, 2 atoms)'
   end if

   ! (b) post-equilibration KE expectation algebra: (3*nat-3)*kb*t_eq/2
   frag_b = thermalize_frag()
   e_expect = (3.0d0*dble(list_atoms%frag(frag_b)%nat) - 3.0d0)*kb_ev*t_use/2.0d0
   if (abs(thermalize_e_kin_expect() - e_expect) > 1.0d-12*e_expect) then
      print '(a,es12.4,a,es12.4)', 'FAIL (b) KE expectation: member ', thermalize_e_kin_expect(), &
         ' vs theory ', e_expect
      nfail = nfail + 1
   else
      print '(a,es12.4,a)', '(b) KE expectation algebra: PASS (e_kin_expect = (3*nat-3)*kb*t_eq/2 = ', &
         e_expect, ')'
   end if

   ! (f) footprint confinement: the sample swaps the propagator step to dt_eq
   ! for the internal run and must RESTORE step and trajectory time; on the
   ! zero-force container q lands exactly on the free-flight mimic (p stays
   ! constant through the run) and f lands exactly 0
   call state_create(sta, 2)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   sta%t = sentinel
   call thermalize_sample(sta)
   qe = sentinel
   do istep = 1, n_use
      do i = 1, 2
         qe(3*i - 2:3*i) = qe(3*i - 2:3*i) + sta%p(3*i - 2:3*i)*(dt_use/list_atoms%mass(i))
      end do
   end do
   if (propagator%dt /= dt_live .or. sta%t /= sentinel .or. any(abs(sta%q - qe(1:6)) > 1.0d-10) .or. &
       any(sta%f /= 0.0d0) .or. any(sta%p == sentinel)) then
      print '(a)', 'FAIL (f) footprint: step/time not restored, q off the free-flight mimic, or f nonzero'
      nfail = nfail + 1
   else
      print '(a)', '(f) footprint confinement: PASS (propagator dt and sta%t restored; q on the'// &
         ' free-flight mimic; f = 0; owned p seeded)'
   end if

   ! (c) drift identity: owned-set momentum sums to zero per component
   ! (exact algebra of the mass-weighted projection - roundoff-level residual)
   rel_max = 0.0d0
   do i = 1, 3
      comp_sum = sta%p(i) + sta%p(i + 3)
      comp_abs = abs(sta%p(i)) + abs(sta%p(i + 3))
      rel = abs(comp_sum)/max(comp_abs, 1.0d-300)
      rel_max = max(rel_max, rel)
   end do
   if (rel_max > 1.0d-10) then
      print '(a,es12.4)', 'FAIL (c) drift identity: max relative component-sum residual ', rel_max
      nfail = nfail + 1
   else
      print '(a,es12.4,a)', '(c) drift identity: PASS (per-component owned momentum sums zero to ', rel_max, ')'
   end if

   ! (d) projection-form replay: the seeded Gaussians replayed here in the
   ! member's draw order (atom-major, component within atom) and projected
   ! by the DECLARED mass-weighted formula must reproduce the member's
   ! output component-for-component (unequal masses 1/4 amu make an
   ! unweighted mean subtraction visibly different)
   call rng_init(7717)
   do i = 1, 2
      do c = 1, 3
         exp_p(3*(i - 1) + c) = sqrt(list_atoms%mass(list_atoms%frag(1)%list(i))*kb_act*t_use)*rng_gauss()
      end do
   end do
   wt = sum(list_atoms%mass(list_atoms%frag(1)%list))
   do c = 1, 3
      ptot = exp_p(c) + exp_p(3 + c)
      exp_p(c) = exp_p(c) - (list_atoms%mass(list_atoms%frag(1)%list(1))/wt)*ptot
      exp_p(3 + c) = exp_p(3 + c) - (list_atoms%mass(list_atoms%frag(1)%list(2))/wt)*ptot
   end do
   call rng_init(7717)
   call thermalize_sample(sta)
   rel = 0.0d0
   do i = 1, 6
      rel = max(rel, abs(sta%p(i) - exp_p(i))/max(abs(exp_p(i)), 1.0d-300))
   end do
   if (rel > 1.0d-12) then
      print '(a,es12.4)', 'FAIL (d) projection-form replay: max relative deviation ', rel
      nfail = nfail + 1
   else
      print '(a,es12.4,a)', '(d) projection-form replay: PASS (declared mass-weighted formula reproduces output to ', rel, ')'
   end if
   call state_destroy(sta)

   ! (a) registration + dispatch chain through samp_run (the run below is
   ! real free flight on the zero-force container: q on the mimic, f = 0)
   call thermalize_register()
   call samp_reg_incident()
   call state_create(sta, 2)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call samp_run('incident', sta)
   qe = sentinel
   do istep = 1, n_use
      do i = 1, 2
         qe(3*i - 2:3*i) = qe(3*i - 2:3*i) + sta%p(3*i - 2:3*i)*(dt_use/list_atoms%mass(i))
      end do
   end do
   if (samp_code('thermalize') /= 9 .or. any(abs(sta%q - qe(1:6)) > 1.0d-10) .or. &
       any(sta%f /= 0.0d0) .or. any(sta%p == sentinel) .or. propagator%dt /= dt_live) then
      print '(a)', 'FAIL (a) dispatch chain: samp_run did not reach the member cleanly'
      nfail = nfail + 1
   else
      print '(a)', '(a) dispatch chain: PASS (reg row -> code 9 -> samp_run dispatch; owned p seeded,'// &
         ' free-flight q, f = 0)'
   end if
   call state_destroy(sta)

   ! ---------- (k) dual-carrier equal treatment: per-carrier seed + drift
   ! removal, then the ONE shared equilibration run (free flight) ----------
   call setup_mol(2)
   call set_scheme([9, 9])
   par%t_eq = [ t_use, t_use ]          ! one temperature set, two rows
   par%n_eq = [ n_use, n_use ]
   par%dt_eq = [ dt_use, dt_use ]
   call thermalize_init(par)               ! re-bind: two carrier rows
   call state_create(sta, 4)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   sta%t = sentinel
   call rng_init(7727)                     ! fresh window: per-criterion reseeding rule
   call thermalize_sample(sta)
   qe = sentinel                            ! the (f)-style free-flight mimic, both fragments
   do istep = 1, n_use
      do i = 1, 4
         qe(3*i - 2:3*i) = qe(3*i - 2:3*i) + sta%p(3*i - 2:3*i)*(dt_use/list_atoms%mass(i))
      end do
   end do
   rel_max = 0.0d0                          ! the (c)-style drift identity, both carriers
   do i = 1, 3
      rel = abs(sta%p(i) + sta%p(i + 3))/max(abs(sta%p(i)) + abs(sta%p(i + 3)), 1.0d-300)
      rel_max = max(rel_max, rel)
      rel = abs(sta%p(6 + i) + sta%p(9 + i))/max(abs(sta%p(6 + i)) + abs(sta%p(9 + i)), 1.0d-300)
      rel_max = max(rel_max, rel)
   end do
   e_expect = 2.0d0*(3.0d0*dble(list_atoms%frag(1)%nat) - 3.0d0)*kb_ev*t_use/2.0d0
   if (thermalize_frag(1) /= 1 .or. thermalize_frag(2) /= 2 .or. thermalize_nat(2) /= 2 .or. &
       abs(thermalize_e_kin_expect() - e_expect) > 1.0d-12*e_expect .or. any(sta%p == sentinel) .or. &
       rel_max > 1.0d-10 .or. any(abs(sta%q - qe) > 1.0d-10) .or. propagator%dt /= dt_live .or. &
       sta%t /= sentinel) then
      print '(a)', 'FAIL (k) dual-carrier equal treatment: a carrier was missed or the shared run left footprints'
      nfail = nfail + 1
   else
      print '(a,es12.4,a)', '(k) dual-carrier equal treatment: PASS (both carriers bound/seeded/'// &
         'drift-removed; e_kin_expect = sum over carriers = ', e_expect, &
         '; one shared free-flight run, step/time restored)'
   end if
   call state_destroy(sta)

   ! ---------- (g)-(m) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-temp',   'thermalize_init', 't_eq',     '(g) non-positive temperature loud-stop')
   call expect_loudstop('bad-nsteps', 'thermalize_init', 'n_eq',     '(h) non-positive step-count loud-stop')
   call expect_loudstop('bad-dt',     'thermalize_init', 'dt_eq',    '(i) non-positive step-size loud-stop')
   call expect_loudstop('bad-bind0',  'thermalize_init', 'scheme',   '(j) no-fragment-binding loud-stop')
   call expect_loudstop('bad-noinit', 'thermalize_draw', 'thermalize_init', '(l) sample-before-init loud-stop')
   call expect_loudstop('bad-offset', 'thermalize_realize', 'state dof =6', '(m) list_atoms-offset short-state loud-stop')

   print '(a,i0,a)', '=== check_samp_thermalize: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! register the incident stand-in (see the header rule)
   subroutine samp_reg_incident()
      use sampler, only: samp_reg
      call samp_reg('incident', stand_incident_d, stand_incident_r)
   end subroutine samp_reg_incident

   ! the check-local zero container: empty kernel and energy (eV member
   ! units) - the free-flight posture of every equilibration run below
   subroutine zero_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A] (unused - zero force)
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A] (zero)
      real(8), intent(out) :: v        ! V [eV] (zero)
      g = 0.0d0
      v = 0.0d0
   end subroutine zero_pes

   ! degenerate no-op method members (the row container_init resolves)
   subroutine noop_prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine noop_prop

   subroutine noop_match(sta)
      type(state_t), intent(inout) :: sta
   end subroutine noop_match

   ! mini-list_atoms by direct field fill: nf two-atom fragments, unequal masses
   ! 1/4 amu (the expectation form must not care, and the drift identity is
   ! mass-blind), molecular paradigm, no scheme preset
   subroutine setup_mol(nf)
      integer, intent(in) :: nf
      integer :: i
      call list_atoms_clear()
      allocate (list_atoms%mass(2*nf))
      do i = 1, nf
         list_atoms%mass(2*i - 1) = 1.0d0
         list_atoms%mass(2*i) = 4.0d0
      end do
      allocate (list_atoms%frag(nf))
      do i = 1, nf
         list_atoms%frag(i)%nat = 2
         list_atoms%frag(i)%mass = 5.0d0
         allocate (list_atoms%frag(i)%list(2), list_atoms%frag(i)%qz(6))
         list_atoms%frag(i)%list = [ 2*i - 1, 2*i ]
         list_atoms%frag(i)%qz = 0.0d0
         list_atoms%frag(i)%qz(1) = -0.5d0
         list_atoms%frag(i)%qz(4) = 0.5d0
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
   ! stderr names the authority and carries the expected key
   subroutine expect_loudstop(cmode, who, key, label)
      character(len=*), intent(in) :: cmode, who, key, label
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_samp_thermalize.exe '//trim(cmode)// &
            ' 2>build'//sep//'equi_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/equi_err.txt', status='old', action='read', iostat=ios)
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
      else if (index(buf, trim(who)) == 0) then
         print '(a)', 'FAIL '//label//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (index(buf, trim(key)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key)//'"'
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

end program check_samp_thermalize
