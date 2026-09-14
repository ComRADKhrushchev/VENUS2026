!=====================================================================
! check_sampler.f90 - sampler framework criteria program (criteria defined
!                     FIRST in this header; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files
! out of the source tree):
!   gfortran -Jbuild globals/config.f90 utils/specio.f90 facade/sysdef_facade.f90 &
!            globals/config_roster.f90 globals/state.f90 facade/sampler.f90 &
!            ../scripts/check_sampler.f90 -o build/check_sampler.exe
! All artifacts are written under build/ - the source tree stays clean.
! Loud-stop verification approach (same scheme as the other criteria
! programs): every loud-stop criterion runs as a CHILD PROCESS - the check
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
!   (a) per-fragment dispatch order + member ownership: two-fragment roster,
!       dist_scheme = the boltzmann,ebk codes, samp_run('incident', sta) calls
!       the fragment members in roster order then the incident member once;
!       each stand-in writes ONLY its own slots (q runs 11/22, p 33) and f
!       keeps its pre-fill sentinel (nobody writes f); a second samp_run call
!       repeats the sequence (the orchestration is per call, incident once
!       each)
!   (b) incident-channel applicability: a single-fragment gas-phase system
!       (surface_model = 0) dispatches NO incident member (log carries only
!       the fragment member, p untouched); a single-fragment surface-paradigm
!       system (surface_model = 1) dispatches it (log grows by the incident
!       tag, p written)
!   (c) count-reconcile loud-stop (child): dist_scheme length 1 against a
!       two-fragment roster -> nonzero exit, stderr names samp_run and BOTH
!       values (length 1, fragment count 2)
!   (d) missing-scheme loud-stop (child): roster assembled, dist_scheme
!       unallocated -> stderr names samp_run and dist_scheme
!   (e) unregistered scheme name loud-stop (child): the scheme argument names
!       a member that is not registered -> stderr names samp_run and the name
!   (f) unregistered fragment member loud-stop (child): a valid scheme code
!       whose member name is not registered -> stderr names samp_run, the
!       member name and the fragment
!   (g) out-of-range scheme code loud-stop (child): code 8 in fragment 2 ->
!       stderr names samp_run, the code and the fragment (sweep-along: the
!       probe moved off 7 when the scheme-code table grew to seven words)
!   (h) duplicate-registration loud-stop (child): samp_reg twice with the
!       same name -> stderr names samp_reg and the duplicate
!   (i) full-table loud-stop (child): 16 rows registered, the 17th -> stderr
!       names samp_reg and the capacity
!   (j) roster assembly-order loud-stop (child): no roster fragment table ->
!       stderr names samp_run and roster_load
!   (k) state-size loud-stop (child): a state sized for 1 atom against a
!       4-atom roster -> stderr names samp_run and BOTH sizes
!=====================================================================
module sampler_standin
   use state, only: state_t
   implicit none
   private
   public :: stand_frag1, stand_frag2, stand_incident, log_reset, log_is
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

   ! stand-in members: each writes ONLY the slots it owns (a fragment member
   ! owns one q run, the incident member owns all of p); the f array is owned
   ! by nobody and must keep its sentinel fill
   subroutine stand_frag1(sta)
      type(state_t), intent(inout) :: sta
      call push(1)
      sta%q(1:3) = 11.0d0
   end subroutine stand_frag1

   subroutine stand_frag2(sta)
      type(state_t), intent(inout) :: sta
      call push(2)
      sta%q(4:6) = 22.0d0
   end subroutine stand_frag2

   subroutine stand_incident(sta)
      type(state_t), intent(inout) :: sta
      call push(3)
      sta%p = 33.0d0
   end subroutine stand_incident
end module sampler_standin

program check_sampler
   use sampler, only: samp_reg, samp_run
   use sampler_standin
   use config_roster, only: roster, roster_clear
   use config, only: reactants
   use state, only: state_t, state_create, state_destroy
   implicit none
   type(state_t) :: sta
   integer :: nfail = 0, nfail0, i
   logical :: ok_gas, ok_surf
   character(len=256) :: mode, name
   real(8), parameter :: sentinel = -777.0d0

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
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
         call samp_reg('boltzmann', stand_frag1)
         call samp_reg('ebk', stand_frag2)
         call setup_mol(2)
         call set_scheme([1, 2])
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-member')   ! (f): code 2 = 'ebk' not registered
         call samp_reg('boltzmann', stand_frag1)
         call samp_reg('incident', stand_incident)
         call setup_mol(2)
         call set_scheme([1, 2])
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-code')     ! (g): code 8 in fragment 2
         call reg_base()
         call setup_mol(2)
         call set_scheme([1, 8])
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-dup')      ! (h): same name twice
         call samp_reg('boltzmann', stand_frag1)
         call samp_reg('boltzmann', stand_frag2)
      case ('bad-full')     ! (i): 16 rows then the 17th
         do i = 1, 16
            write (name, '(a,i0)') 'stand', i
            call samp_reg(trim(name), stand_frag1)
         end do
         call samp_reg('overflow', stand_frag2)
      case ('bad-order')    ! (j): roster%frag left unallocated
         call reg_base()
         call state_create(sta, 4)
         call samp_run('incident', sta)
      case ('bad-state')    ! (k): 3-dof state vs a 4-atom roster
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
   if (.not. log_is([1, 2, 3, 1, 2, 3]) .or. &
       any(sta%q(1:3) /= 11.0d0) .or. any(sta%q(4:6) /= 22.0d0) .or. &
       any(sta%p /= 33.0d0) .or. any(sta%f /= sentinel)) then
      print '(a)', 'FAIL (a) per-fragment dispatch: order/ownership mismatch'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(a) per-fragment dispatch: PASS (order boltzmann,ebk,incident per call; q runs 11/22, p 33, f sentinel)'
   call state_destroy(sta)

   ! ---------- (b) incident-channel applicability ----------
   nfail0 = nfail
   call setup_mol(1)
   call set_scheme([1])
   call state_create(sta, 2)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call log_reset()
   call samp_run('incident', sta)   ! single-fragment gas phase: no incident channel
   ok_gas = log_is([1]) .and. all(sta%q(1:3) == 11.0d0) .and. all(sta%q(4:6) == sentinel) .and. &
            all(sta%p == sentinel)
   reactants%surface_model = 1
   sta%q = sentinel
   sta%p = sentinel
   call log_reset()
   call samp_run('incident', sta)   ! single-fragment surface paradigm: channel on
   ok_surf = log_is([1, 3]) .and. all(sta%p == 33.0d0)
   reactants%surface_model = 0
   if (.not. ok_gas .or. .not. ok_surf) then
      print '(a)', 'FAIL (b) incident applicability: gas-phase skip or surface dispatch wrong'
      nfail = nfail + 1
   end if
   if (nfail == nfail0) print '(a)', &
      '(b) incident applicability: PASS (single-fragment gas phase skips the channel, surface paradigm dispatches it)'
   call state_destroy(sta)

   ! ---------- (c)-(k) loud-stop criteria (child processes) ----------
   call expect_loudstop('bad-count',    'samp_run', 'length 1',         'count 2', &
                        '(c) count-reconcile loud-stop')
   call expect_loudstop('bad-noscheme', 'samp_run', 'dist_scheme',      '', &
                        '(d) missing-scheme loud-stop')
   call expect_loudstop('bad-scheme',   'samp_run', 'incident',         'not a registered', &
                        '(e) unregistered-scheme loud-stop')
   call expect_loudstop('bad-member',   'samp_run', 'ebk',              'fragment 2', &
                        '(f) unregistered-member loud-stop')
   call expect_loudstop('bad-code',     'samp_run', 'code 8',           'fragment 2', &
                        '(g) out-of-range-code loud-stop')
   call expect_loudstop('bad-dup',      'samp_reg', 'boltzmann',        'already registered', &
                        '(h) duplicate-registration loud-stop')
   call expect_loudstop('bad-full',     'samp_reg', '16 rows',          'overflow', &
                        '(i) full-table loud-stop')
   call expect_loudstop('bad-order',    'samp_run', 'roster_load',      'assembly-order', &
                        '(j) roster-order loud-stop')
   call expect_loudstop('bad-state',    'samp_run', '3 dof',            '12 dof', &
                        '(k) state-size loud-stop')

   print '(a,i0,a)', '=== check_sampler: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! register the three stand-ins under the real member names
   subroutine reg_base()
      call samp_reg('boltzmann', stand_frag1)
      call samp_reg('ebk', stand_frag2)
      call samp_reg('incident', stand_incident)
   end subroutine reg_base

   ! mini-roster built by direct field fill (the roster public surface allows
   ! it; the sysdef->roster_load chain has its own criteria program): two
   ! atoms per fragment, molecular paradigm, no scheme preset
   subroutine setup_mol(nf)
      integer, intent(in) :: nf
      integer :: i
      call roster_clear()
      allocate (roster%mass(2*nf))
      roster%mass = 1.0d0
      allocate (roster%frag(nf))
      do i = 1, nf
         roster%frag(i)%nat = 2
         roster%frag(i)%mass = 2.0d0
         allocate (roster%frag(i)%list(2), roster%frag(i)%qz(6))
         roster%frag(i)%list = [ 2*i - 1, 2*i ]
         roster%frag(i)%qz = 0.0d0
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
