!=====================================================================
! check_dist_thermostat.f90 - thermostat distribution-member criteria
!                    program (criteria defined FIRST in this header;
!                    this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 &
!            globals/config_roster.f90 globals/state.f90 &
!            utils/rng.f90 facade/sampler.f90 methods/dist_thermostat.f90 &
!            methods/seam_cast.f90 methods/reg_dist_thermostat.f90 &
!            ../scripts/check_dist_thermostat.f90 -o build/check_dist_thermostat.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter statute: the member is driven by directly constructed
! parameters and a directly filled roster; the input grammar is not in the
! loop. This is a MECHANISM-level check (the force-chain-dependent numerics
! are posterior): no sampling statistics are asserted - the distribution
! declarations are verified ALGEBRAICALLY from the member-public parameters.
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
!       reg file row, samp_code('thermostat') = 8; a single-fragment
!       gas-phase roster with dist_scheme = [8] and samp_run('incident',
!       sta) REACHES the member - the owned p components come off their
!       sentinels while q/f keep theirs (the member writes only its owned
!       momentum components; the incident channel is not dispatched)
!   (b) Maxwell-shape algebra (analytic, no sampling): the declared
!       per-component momentum width obeys sigma(i)^2 = m_i*kb*t_bath for
!       every owned atom (unequal masses 1/4 amu, m_i from the roster, kb
!       from consts, sigma from the member getter), relative 1e-12;
!       parameter-loading integrity: the t_bath and nu getters return the
!       values init received (the collision frequency is carried bath
!       data for the dynamics side - loaded, not dropped)
!   (c) component ownership: a two-fragment roster with only fragment 1
!       carrying the code; a DIRECT thermostat_sample call leaves q and f
!       on their sentinels everywhere, leaves the unbound fragment's p on
!       its sentinel, and draws the owned fragment's p off the sentinel
!   (d) non-positive temperature loud-stop (child): t_bath = 0 -> stderr
!       names thermostat_init and t_bath
!   (e) zero collision frequency loud-stop (child): nu_coll = 0 -> stderr
!       names thermostat_init and nu_coll
!   (f) no-fragment-binding loud-stop (child): no dist_scheme row carries
!       the code -> stderr names thermostat_init
!   (g) duplicate-binding loud-stop (child): two fragments carry the
!       code -> stderr names thermostat_init
!   (h) sample-before-init loud-stop (child): thermostat_sample with no
!       thermostat_init -> stderr names thermostat_sample
!   (i) undersized-state loud-stop (child): a 3-dof state against a
!       2-atom carrier -> stderr names thermostat_sample and the state
!   (j) roster-offset short-state loud-stop (child): a 6-dof state
!       against a SECOND-fragment carrier (2x2-atom roster, scheme
!       [1,8]: the owned components are 7..12) -> stderr names
!       thermostat_sample and the state (the guard must see the
!       fragment's roster offset, not just the atom count)
!=====================================================================
module bath_standin
   use state, only: state_t
   implicit none
   private
   public :: stand_incident
contains
   ! the incident-channel stand-in: registered so samp_run's scheme-argument
   ! lookup resolves; NOT dispatched in these scenarios (single-fragment gas
   ! phase) - this check's subject is the thermostat member
   subroutine stand_incident(sta)
      type(state_t), intent(inout) :: sta
   end subroutine stand_incident
end module bath_standin

program check_dist_thermostat
   use consts,        only: kb_code
   use rng,           only: rng_init
   use sampler,       only: samp_run, samp_code
   use config,        only: reactants
   use config_roster, only: roster, roster_clear
   use state,         only: state_t, state_create, state_destroy
   use dist_thermostat, only: thermostat_params_t, thermostat_init, thermostat_sample, &
                              thermostat_t_bath, thermostat_nu, thermostat_frag, &
                              thermostat_nat, thermostat_sigma
   use reg_dist_thermostat, only: reg_dist_thermostat_all
   use bath_standin, only: stand_incident
   implicit none
   type(state_t) :: sta
   type(thermostat_params_t) :: par
   integer :: nfail = 0, i, g_atom, frag_b
   real(8) :: sig_theory, dev_max
   real(8), parameter :: t_use = 450.0d0
   real(8), parameter :: nu_use = 0.137d0
   real(8), parameter :: sentinel = -777.0d0
   character(len=32) :: mode

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 2)
      par%t_bath = t_use
      par%nu_coll = nu_use
      select case (trim(mode))
      case ('bad-temp')                     ! (d): t_bath = 0
         call setup_mol(1)
         call set_scheme([8])
         par%t_bath = 0.0d0
         call thermostat_init(par)
      case ('bad-nu')                       ! (e): collision frequency zero
         call setup_mol(1)
         call set_scheme([8])
         par%nu_coll = 0.0d0
         call thermostat_init(par)
      case ('bad-bind0')                    ! (f): no fragment carries the code
         call setup_mol(1)
         call set_scheme([1])
         call thermostat_init(par)
      case ('bad-bind2')                    ! (g): two fragments carry the code
         call setup_mol(2)
         call set_scheme([8, 8])
         call thermostat_init(par)
      case ('bad-noinit')                   ! (h): sample without init
         call setup_mol(1)
         call set_scheme([8])
         call thermostat_sample(sta)
      case ('bad-state')                    ! (i): 3-dof state vs a 2-atom carrier
         call setup_mol(1)
         call set_scheme([8])
         call thermostat_init(par)
         call state_destroy(sta)
         call state_create(sta, 1)
         call thermostat_sample(sta)
      case ('bad-offset')                   ! (j): 6-dof state vs a SECOND-fragment carrier
         call setup_mol(2)
         call set_scheme([1, 8])
         call thermostat_init(par)
         call thermostat_sample(sta)        ! sta holds 2 atoms (6 dof) - owned 7..12
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- single carrier: algebra + registration/dispatch ----------
   call setup_mol(1)
   call set_scheme([8])
   par%t_bath = t_use
   par%nu_coll = nu_use
   call rng_init(5501)
   call thermostat_init(par)

   ! (b) Maxwell-shape algebra: sigma(i)^2 = m_i*kb*t_bath, plus the
   ! parameter-loading integrity of t_bath and the collision frequency
   frag_b = thermostat_frag()
   dev_max = 0.0d0
   do i = 1, thermostat_nat()
      g_atom = roster%frag(frag_b)%list(i)
      sig_theory = sqrt(roster%mass(g_atom)*kb_code*t_use)
      dev_max = max(dev_max, abs(thermostat_sigma(i) - sig_theory)/sig_theory)
   end do
   if (dev_max > 1.0d-12 .or. thermostat_frag() /= 1 .or. thermostat_nat() /= 2 .or. &
       thermostat_t_bath() /= t_use .or. thermostat_nu() /= nu_use) then
      print '(a,es12.4)', 'FAIL (b) Maxwell algebra: max relative sigma deviation ', dev_max
      nfail = nfail + 1
   else
      print '(a,es12.4,a)', '(b) Maxwell-shape algebra: PASS (sigma^2 = m*kb*t_bath to ', dev_max, &
         '; t_bath/nu loaded intact)'
   end if

   ! (a) registration + dispatch chain through samp_run
   call reg_dist_thermostat_all()
   call samp_reg_incident()
   call state_create(sta, 2)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call samp_run('incident', sta)
   if (samp_code('thermostat') /= 8 .or. any(sta%q /= sentinel) .or. &
       any(sta%f /= sentinel) .or. any(sta%p == sentinel)) then
      print '(a)', 'FAIL (a) dispatch chain: samp_run did not reach the member cleanly'
      nfail = nfail + 1
   else
      print '(a)', '(a) dispatch chain: PASS (reg row -> code 8 -> samp_run dispatch; only owned p written)'
   end if
   call state_destroy(sta)

   ! ---------- (c) component ownership: two fragments, direct sample ----------
   call setup_mol(2)
   call set_scheme([8, 1])
   call thermostat_init(par)               ! re-bind to fragment 1 of the new roster
   call state_create(sta, 4)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call thermostat_sample(sta)
   if (any(sta%q /= sentinel) .or. any(sta%f /= sentinel) .or. &
       any(sta%p(7:12) /= sentinel) .or. any(sta%p(1:6) == sentinel)) then
      print '(a)', 'FAIL (c) ownership: trespass outside the owned p components'
      nfail = nfail + 1
   else
      print '(a)', '(c) component ownership: PASS (q/f untouched; unbound fragment p untouched; owned p drawn)'
   end if
   call state_destroy(sta)

   ! ---------- (d)-(i) loud-stop criteria (child processes) ----------
   call expect_loudstop('bad-temp',   'thermostat_init',   't_bath',          '(d) non-positive temperature loud-stop')
   call expect_loudstop('bad-nu',     'thermostat_init',   'nu_coll',         '(e) zero collision frequency loud-stop')
   call expect_loudstop('bad-bind0',  'thermostat_init',   'thermostat scheme', '(f) no-fragment-binding loud-stop')
   call expect_loudstop('bad-bind2',  'thermostat_init',   'two or more',     '(g) duplicate-binding loud-stop')
   call expect_loudstop('bad-noinit', 'thermostat_sample', 'thermostat_init', '(h) sample-before-init loud-stop')
   call expect_loudstop('bad-state',  'thermostat_sample', 'state',           '(i) undersized-state loud-stop')
   call expect_loudstop('bad-offset', 'thermostat_sample', 'state dof =6',    '(j) roster-offset short-state loud-stop')

   print '(a,i0,a)', '=== check_dist_thermostat: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! register the incident stand-in (see the header statute)
   subroutine samp_reg_incident()
      use sampler, only: samp_reg
      call samp_reg('incident', stand_incident)
   end subroutine samp_reg_incident

   ! mini-roster by direct field fill: nf two-atom fragments, unequal masses
   ! 1/4 amu (mass-dependent widths make the algebra sharp), molecular
   ! paradigm, no scheme preset
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
      cmd = 'build'//sep//'check_dist_thermostat.exe '//trim(cmode)// &
            ' 2>build'//sep//'thermo_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/thermo_err.txt', status='old', action='read', iostat=ios)
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

end program check_dist_thermostat
