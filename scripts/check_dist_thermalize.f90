!=====================================================================
! check_dist_thermalize.f90 - thermalize distribution-member criteria
!                    program (criteria defined FIRST in this header;
!                    this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree; the evolution chain enters because the
! member drives the propagator through its no-op skeleton):
!   gfortran -Jbuild globals/consts.f90 globals/config.f90 utils/specio.f90 &
!            facade/sysdef_facade.f90 globals/config_roster.f90 globals/state.f90 &
!            utils/rng.f90 facade/elec_ws.f90 facade/elec_facade.f90 &
!            facade/force_facade.f90 dynamics/propagator.f90 facade/sampler.f90 &
!            methods/dist_thermalize.f90 methods/reg_dist_thermalize.f90 &
!            ../scripts/check_dist_thermalize.f90 -o build/check_dist_thermalize.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter statute: the member is driven by directly constructed
! parameters and a directly filled roster; the input grammar is not in the
! loop. This is a MECHANISM-level check: no trajectory statistics are
! asserted (the equilibration RUN is exercised through the no-op force
! chain - the zero-force trajectory is the lawful mechanism-level
! posture, real dynamics are posterior); the kinetic-energy expectation and
! the drift identity are verified ALGEBRAICALLY from the member-public
! parameters and the exact transformation the member applies.
! The incident channel is not this check's subject: a LOCAL no-op stand-in
! is registered under the name 'incident' only because samp_run resolves
! the scheme argument before dispatch (the channel itself is skipped for
! the single-fragment gas-phase roster used here).
! Loud-stop verification approach (same scheme as the sibling checks):
! every loud-stop criterion runs as a CHILD PROCESS - the check re-spawns
! itself with a mode argument (execute_command_line); the child drives the
! REAL member entry into the stop and must print a naming message on
! stderr and exit nonzero. The parent asserts the exit status and the
! stderr content.
! Criteria:
!   (a) registration + dispatch chain: the member registers through its
!       reg file row, samp_code('thermalize') = 10; a single-fragment
!       gas-phase roster with dist_scheme = [10] and samp_run('incident',
!       sta) REACHES the member - the owned p components come off their
!       sentinels (the Maxwell seed) while q/f keep theirs (the no-op
!       force chain advances nothing - zero-force trajectory); the
!       incident channel is not dispatched
!   (b) post-equilibration kinetic-energy expectation algebra (analytic,
!       no sampling): the declared expectation obeys
!       e_kin_expect = (3*nat - 3)*kb*t_eq/2 for the bound carrier (the
!       Gaussian projection onto zero net momentum removes exactly three
!       degrees of freedom - independent of the mass distribution),
!       relative 1e-12
!   (c) drift identity (exact algebra of the transformation): after a
!       direct sample call the owned-set momentum sums to zero PER
!       COMPONENT (each component had its owned-set mean subtracted) -
!       relative residual <= 1e-10 (floating-point roundoff of the
!       subtraction, not a statistical claim)
!   (d) parameter-loading integrity: the t_eq/n_eq/dt_eq getters return
!       the values init received
!   (e) footprint confinement: with the propagator step preset to a
!       distinctive value and sta%t preset to a sentinel, a sample call
!       (which swaps the step to dt_eq for the internal run) RESTORES
!       both - the equilibration run neither consumes trajectory time nor
!       leaves the global step changed
!   (f) non-positive temperature loud-stop (child): t_eq = 0 -> stderr
!       names thermalize_init and t_eq
!   (g) non-positive step count loud-stop (child): n_eq = 0 -> stderr
!       names thermalize_init and n_eq
!   (h) non-positive step size loud-stop (child): dt_eq = 0 -> stderr
!       names thermalize_init and dt_eq
!   (i) no-fragment-binding loud-stop (child): no dist_scheme row carries
!       the code -> stderr names thermalize_init
!   (j) duplicate-binding loud-stop (child): two fragments carry the
!       code -> stderr names thermalize_init
!   (k) sample-before-init loud-stop (child): thermalize_sample with no
!       thermalize_init -> stderr names thermalize_sample
!=====================================================================
module equi_standin
   use state, only: state_t
   implicit none
   private
   public :: stand_incident
contains
   ! the incident-channel stand-in: registered so samp_run's scheme-argument
   ! lookup resolves; NOT dispatched in these scenarios (single-fragment gas
   ! phase) - this check's subject is the thermalize member
   subroutine stand_incident(sta)
      type(state_t), intent(inout) :: sta
   end subroutine stand_incident
end module equi_standin

program check_dist_thermalize
   use consts,        only: kb_code
   use rng,           only: rng_init
   use sampler,       only: samp_run, samp_code
   use config,        only: reactants, propagator
   use config_roster, only: roster, roster_clear
   use state,         only: state_t, state_create, state_destroy
   use dist_thermalize, only: thermalize_params_t, thermalize_init, thermalize_sample, &
                              thermalize_t_eq, thermalize_n_eq, thermalize_dt_eq, &
                              thermalize_e_kin_expect, thermalize_frag, thermalize_nat
   use reg_dist_thermalize, only: reg_dist_thermalize_all
   use equi_standin, only: stand_incident
   implicit none
   type(state_t) :: sta
   type(thermalize_params_t) :: par
   integer :: nfail = 0, i, frag_b
   real(8) :: e_expect, comp_sum, comp_abs, rel, rel_max
   real(8), parameter :: t_use = 325.0d0
   real(8), parameter :: dt_use = 0.0413d0
   real(8), parameter :: dt_live = 0.031d0
   integer, parameter :: n_use = 3
   real(8), parameter :: sentinel = -777.0d0
   character(len=32) :: mode

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 2)
      par%t_eq = t_use
      par%n_eq = n_use
      par%dt_eq = dt_use
      propagator%dt = dt_live
      select case (trim(mode))
      case ('bad-temp')                     ! (f): t_eq = 0
         call setup_mol(1)
         call set_scheme([10])
         par%t_eq = 0.0d0
         call thermalize_init(par)
      case ('bad-nsteps')                   ! (g): zero propagation steps
         call setup_mol(1)
         call set_scheme([10])
         par%n_eq = 0
         call thermalize_init(par)
      case ('bad-dt')                       ! (h): zero equilibration step
         call setup_mol(1)
         call set_scheme([10])
         par%dt_eq = 0.0d0
         call thermalize_init(par)
      case ('bad-bind0')                    ! (i): no fragment carries the code
         call setup_mol(1)
         call set_scheme([1])
         call thermalize_init(par)
      case ('bad-bind2')                    ! (j): two fragments carry the code
         call setup_mol(2)
         call set_scheme([10, 10])
         call thermalize_init(par)
      case ('bad-noinit')                   ! (k): sample without init
         call setup_mol(1)
         call set_scheme([10])
         call thermalize_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- single carrier: algebra + registration/dispatch ----------
   call setup_mol(1)
   call set_scheme([10])
   par%t_eq = t_use
   par%n_eq = n_use
   par%dt_eq = dt_use
   propagator%dt = dt_live
   call rng_init(7717)
   call thermalize_init(par)

   ! (d) parameter-loading integrity
   if (thermalize_t_eq() /= t_use .or. thermalize_n_eq() /= n_use .or. &
       thermalize_dt_eq() /= dt_use .or. thermalize_frag() /= 1 .or. thermalize_nat() /= 2) then
      print '(a)', 'FAIL (d) parameter integrity: an equilibration parameter or the carrier was dropped'
      nfail = nfail + 1
   else
      print '(a)', '(d) parameter integrity: PASS (t_eq/n_eq/dt_eq round-trip; carrier fragment 1, 2 atoms)'
   end if

   ! (b) post-equilibration KE expectation algebra: (3*nat-3)*kb*t_eq/2
   frag_b = thermalize_frag()
   e_expect = (3.0d0*dble(roster%frag(frag_b)%nat) - 3.0d0)*kb_code*t_use/2.0d0
   if (abs(thermalize_e_kin_expect() - e_expect) > 1.0d-12*e_expect) then
      print '(a,es12.4,a,es12.4)', 'FAIL (b) KE expectation: member ', thermalize_e_kin_expect(), &
         ' vs theory ', e_expect
      nfail = nfail + 1
   else
      print '(a,es12.4,a)', '(b) KE expectation algebra: PASS (e_kin_expect = (3*nat-3)*kb*t_eq/2 = ', &
         e_expect, ')'
   end if

   ! (e) footprint confinement: the sample swaps the propagator step to dt_eq
   ! for the internal run and must RESTORE step and trajectory time
   call state_create(sta, 2)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   sta%t = sentinel
   call thermalize_sample(sta)
   if (propagator%dt /= dt_live .or. sta%t /= sentinel .or. any(sta%q /= sentinel) .or. &
       any(sta%f /= sentinel) .or. any(sta%p == sentinel)) then
      print '(a)', 'FAIL (e) footprint: step/time not restored or q/f trespassed'
      nfail = nfail + 1
   else
      print '(a)', '(e) footprint confinement: PASS (propagator dt and sta%t restored; q/f untouched; owned p seeded)'
   end if

   ! (c) drift identity: owned-set momentum sums to zero per component
   ! (exact algebra of the mean subtraction - roundoff-level residual)
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
   call state_destroy(sta)

   ! (a) registration + dispatch chain through samp_run
   call reg_dist_thermalize_all()
   call samp_reg_incident()
   call state_create(sta, 2)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call samp_run('incident', sta)
   if (samp_code('thermalize') /= 10 .or. any(sta%q /= sentinel) .or. &
       any(sta%f /= sentinel) .or. any(sta%p == sentinel) .or. propagator%dt /= dt_live) then
      print '(a)', 'FAIL (a) dispatch chain: samp_run did not reach the member cleanly'
      nfail = nfail + 1
   else
      print '(a)', '(a) dispatch chain: PASS (reg row -> code 10 -> samp_run dispatch; owned p seeded, zero-force run)'
   end if
   call state_destroy(sta)

   ! ---------- (f)-(k) loud-stop criteria (child processes) ----------
   call expect_loudstop('bad-temp',   'thermalize_init', 't_eq',     '(f) non-positive temperature loud-stop')
   call expect_loudstop('bad-nsteps', 'thermalize_init', 'n_eq',     '(g) non-positive step-count loud-stop')
   call expect_loudstop('bad-dt',     'thermalize_init', 'dt_eq',    '(h) non-positive step-size loud-stop')
   call expect_loudstop('bad-bind0',  'thermalize_init', 'scheme',   '(i) no-fragment-binding loud-stop')
   call expect_loudstop('bad-bind2',  'thermalize_init', 'two',      '(j) duplicate-binding loud-stop')
   call expect_loudstop('bad-noinit', 'thermalize_sample', 'thermalize_init', '(k) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_dist_thermalize: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! register the incident stand-in (see the header statute)
   subroutine samp_reg_incident()
      use sampler, only: samp_reg
      call samp_reg('incident', stand_incident)
   end subroutine samp_reg_incident

   ! mini-roster by direct field fill: nf two-atom fragments, unequal masses
   ! 1/4 amu (the expectation form must not care, and the drift identity is
   ! mass-blind), molecular paradigm, no scheme preset
   subroutine setup_mol(nf)
      integer, intent(in) :: nf
      integer :: i
      call roster_clear()
      allocate (roster%mass(2*nf))
      do i = 1, nf
         roster%mass(2*i - 1) = 1.0d0
         roster%mass(2*i) = 4.0d0
      end do
      allocate (roster%frag(nf))
      do i = 1, nf
         roster%frag(i)%nat = 2
         roster%frag(i)%mass = 5.0d0
         allocate (roster%frag(i)%list(2), roster%frag(i)%qz(6))
         roster%frag(i)%list = [ 2*i - 1, 2*i ]
         roster%frag(i)%qz = 0.0d0
         roster%frag(i)%qz(1) = -0.5d0
         roster%frag(i)%qz(4) = 0.5d0
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

   ! spawn the child for one bad case and assert the loud-stop: nonzero exit,
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
      cmd = 'build'//sep//'check_dist_thermalize.exe '//trim(cmode)// &
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

end program check_dist_thermalize
