!=====================================================================
! check_dist_glo_target.f90 - glo_target distribution-member criteria
!                    program (criteria defined FIRST in this header;
!                    this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/config.f90 utils/specio.f90 &
!            facade/sysdef_facade.f90 globals/config_roster.f90 globals/state.f90 &
!            utils/rng.f90 facade/sampler.f90 methods/dist_glo_target.f90 &
!            methods/reg_dist_glo_target.f90 &
!            ../scripts/check_dist_glo_target.f90 -o build/check_dist_glo_target.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter statute: the member is driven by directly constructed
! parameters and a directly filled roster; the input grammar is not in the
! loop. This is a MECHANISM-level check: no sampling statistics are
! asserted - the fluctuation-dissipation and stationary-width declarations
! are verified ALGEBRAICALLY from the member-public parameters (the ghost
! Langevin INTEGRATION belongs to the force-chain side and is not this
! member's product).
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
!       reg file row, samp_code('glo_target') = 9; a single-fragment
!       gas-phase roster with dist_scheme = [9] and samp_run('incident',
!       sta) REACHES the member while leaving sta COMPLETELY unchanged
!       (q/p/f all keep their sentinels - ghosts are member product, not
!       roster atoms; the roster is the sole composition channel) and the
!       ghost-state getters come back populated; a second dispatch
!       redraws the ghost state (fresh draw per sample call)
!   (b) selection rule: one ghost per carrier atom, isotropic three
!       components each - the ghost-state size is 3*nat of the bound
!       carrier, and the carrier/atom getters name the bound fragment
!   (c) fluctuation-dissipation algebra (analytic, no sampling): the
!       loaded bath noise strength obeys
!       sigma_noise(i)^2 = 2*gamma_glo*m_i*kb*t_glo/dt for every carrier
!       atom (unequal masses 1/4 amu, dt the propagator step the member
!       read at init), relative 1e-12
!   (d) stationary-width algebra (analytic): sigma_p(i)^2 = m_i*kb*t_glo
!       and sigma_q(i)^2 = kb*t_glo/(m_i*w_ghost^2) for every carrier
!       atom (the stationary state of the damped oscillator at bath
!       temperature), relative 1e-12
!   (e) parameter-loading integrity: the t_glo/gamma_glo/w_ghost getters
!       return the values init received, and the dt getter returns the
!       propagator step the member read at init (the coupling parameters
!       are loaded bath data for the dynamics side - loaded, not dropped)
!   (f) non-positive bath temperature loud-stop (child): t_glo = 0 ->
!       stderr names glo_target_init and t_glo
!   (g) zero damping loud-stop (child): gamma_glo = 0 -> stderr names
!       glo_target_init and gamma_glo
!   (h) zero ghost frequency loud-stop (child): w_ghost = 0 -> stderr
!       names glo_target_init and w_ghost
!   (i) zero propagator step loud-stop (child): propagator dt = 0 at
!       init -> stderr names glo_target_init and the step
!   (j) no-fragment-binding loud-stop (child): no dist_scheme row carries
!       the code -> stderr names glo_target_init
!   (k) duplicate-binding loud-stop (child): two fragments carry the
!       code -> stderr names glo_target_init
!   (l) sample-before-init loud-stop (child): glo_target_sample with no
!       glo_target_init -> stderr names glo_target_sample
!=====================================================================
module glo_standin
   use state, only: state_t
   implicit none
   private
   public :: stand_incident
contains
   ! the incident-channel stand-in: registered so samp_run's scheme-argument
   ! lookup resolves; NOT dispatched in these scenarios (single-fragment gas
   ! phase) - this check's subject is the glo_target member
   subroutine stand_incident(sta)
      type(state_t), intent(inout) :: sta
   end subroutine stand_incident
end module glo_standin

program check_dist_glo_target
   use consts,        only: kb_code
   use rng,           only: rng_init
   use sampler,       only: samp_run, samp_code
   use config,        only: reactants, propagator
   use config_roster, only: roster, roster_clear
   use state,         only: state_t, state_create, state_destroy
   use dist_glo_target, only: glo_target_params_t, glo_target_init, glo_target_sample, &
                              glo_target_t_glo, glo_target_gamma, glo_target_w, &
                              glo_target_dt_loaded, glo_target_frag, glo_target_nat, &
                              glo_target_ghost_size, glo_target_last_q, glo_target_last_p, &
                              glo_target_sigma_noise, glo_target_sigma_p, glo_target_sigma_q
   use reg_dist_glo_target, only: reg_dist_glo_target_all
   use glo_standin, only: stand_incident
   implicit none
   type(state_t) :: sta
   type(glo_target_params_t) :: par
   integer :: nfail = 0, i, g_atom, frag_b, nghost
   real(8) :: s_theory, dev_max, ghost_sum1, ghost_sum2
   real(8) :: draw1(64), draw2(64)
   real(8), parameter :: t_use = 550.0d0
   real(8), parameter :: g_use = 0.09d0
   real(8), parameter :: w_use = 0.0021d0
   real(8), parameter :: dt_use = 0.02d0
   real(8), parameter :: sentinel = -777.0d0
   character(len=32) :: mode

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 2)
      par%t_glo = t_use
      par%gamma_glo = g_use
      par%w_ghost = w_use
      propagator%dt = dt_use
      select case (trim(mode))
      case ('bad-temp')                     ! (f): t_glo = 0
         call setup_mol(1)
         call set_scheme([9])
         par%t_glo = 0.0d0
         call glo_target_init(par)
      case ('bad-gamma')                    ! (g): damping zero
         call setup_mol(1)
         call set_scheme([9])
         par%gamma_glo = 0.0d0
         call glo_target_init(par)
      case ('bad-w')                        ! (h): ghost frequency zero
         call setup_mol(1)
         call set_scheme([9])
         par%w_ghost = 0.0d0
         call glo_target_init(par)
      case ('bad-dt')                       ! (i): propagator step zero at init
         call setup_mol(1)
         call set_scheme([9])
         propagator%dt = 0.0d0
         call glo_target_init(par)
      case ('bad-bind0')                    ! (j): no fragment carries the code
         call setup_mol(1)
         call set_scheme([1])
         call glo_target_init(par)
      case ('bad-bind2')                    ! (k): two fragments carry the code
         call setup_mol(2)
         call set_scheme([9, 9])
         call glo_target_init(par)
      case ('bad-noinit')                   ! (l): sample without init
         call setup_mol(1)
         call set_scheme([9])
         call glo_target_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- single carrier: registration/dispatch + algebra ----------
   call setup_mol(1)
   call set_scheme([9])
   par%t_glo = t_use
   par%gamma_glo = g_use
   par%w_ghost = w_use
   propagator%dt = dt_use
   call rng_init(6613)
   call glo_target_init(par)

   ! (e) parameter-loading integrity: t_glo/gamma/w round-trip and the dt
   ! the member read at init is the propagator step it saw
   if (glo_target_t_glo() /= t_use .or. glo_target_gamma() /= g_use .or. &
       glo_target_w() /= w_use .or. glo_target_dt_loaded() /= dt_use) then
      print '(a)', 'FAIL (e) parameter integrity: a bath parameter or the loaded dt was dropped'
      nfail = nfail + 1
   else
      print '(a)', '(e) parameter integrity: PASS (t_glo/gamma_glo/w_ghost round-trip; dt loaded from propagator)'
   end if

   ! (b) selection rule: one ghost per carrier atom x 3 components
   if (glo_target_frag() /= 1 .or. glo_target_nat() /= 2 .or. &
       glo_target_ghost_size() /= 6) then
      print '(a)', 'FAIL (b) selection rule: ghost set does not match the carrier atom set'
      nfail = nfail + 1
   else
      print '(a)', '(b) selection rule: PASS (one ghost x 3 components per carrier atom: 2 atoms -> 6)'
   end if

   ! (c) fluctuation-dissipation algebra: sigma_noise^2 = 2*gamma*m*kb*T/dt
   frag_b = glo_target_frag()
   dev_max = 0.0d0
   do i = 1, glo_target_nat()
      g_atom = roster%frag(frag_b)%list(i)
      s_theory = sqrt(2.0d0*g_use*roster%mass(g_atom)*kb_code*t_use/dt_use)
      dev_max = max(dev_max, abs(glo_target_sigma_noise(i) - s_theory)/s_theory)
   end do
   if (dev_max > 1.0d-12) then
      print '(a,es12.4)', 'FAIL (c) FDT algebra: max relative noise-strength deviation ', dev_max
      nfail = nfail + 1
   else
      print '(a,es12.4,a)', '(c) FDT algebra: PASS (sigma_noise^2 = 2*gamma*m*kb*t/dt to ', dev_max, ')'
   end if

   ! (d) stationary-width algebra: sigma_p^2 = m*kb*T, sigma_q^2 = kb*T/(m*w^2)
   dev_max = 0.0d0
   do i = 1, glo_target_nat()
      g_atom = roster%frag(frag_b)%list(i)
      s_theory = sqrt(roster%mass(g_atom)*kb_code*t_use)
      dev_max = max(dev_max, abs(glo_target_sigma_p(i) - s_theory)/s_theory)
      s_theory = sqrt(kb_code*t_use/(roster%mass(g_atom)*w_use**2))
      dev_max = max(dev_max, abs(glo_target_sigma_q(i) - s_theory)/s_theory)
   end do
   if (dev_max > 1.0d-12) then
      print '(a,es12.4)', 'FAIL (d) stationary-width algebra: max relative deviation ', dev_max
      nfail = nfail + 1
   else
      print '(a,es12.4,a)', '(d) stationary-width algebra: PASS (sigma_p^2 = m*kb*t; sigma_q^2 = kb*t/(m*w^2) to ', dev_max, ')'
   end if

   ! (a) registration + dispatch chain through samp_run: sta COMPLETELY
   ! unchanged (the ghost state is member product), ghost getters populated,
   ! and a second dispatch redraws the ghost state
   call reg_dist_glo_target_all()
   call samp_reg_incident()
   call state_create(sta, 2)
   sta%q = sentinel
   sta%p = sentinel
   sta%f = sentinel
   call samp_run('incident', sta)
   nghost = glo_target_ghost_size()
   ghost_sum1 = 0.0d0
   do i = 1, nghost
      draw1(i) = glo_target_last_q(i)
      draw2(i) = glo_target_last_p(i)
      ghost_sum1 = ghost_sum1 + abs(draw1(i)) + abs(draw2(i))
   end do
   call samp_run('incident', sta)
   ghost_sum2 = 0.0d0
   do i = 1, nghost
      ghost_sum2 = ghost_sum2 + abs(glo_target_last_q(i) - draw1(i)) + &
                   abs(glo_target_last_p(i) - draw2(i))
   end do
   if (samp_code('glo_target') /= 9 .or. any(sta%q /= sentinel) .or. &
       any(sta%p /= sentinel) .or. any(sta%f /= sentinel) .or. &
       nghost /= 6 .or. ghost_sum1 <= 0.0d0 .or. ghost_sum2 <= 0.0d0) then
      print '(a)', 'FAIL (a) dispatch chain: sta trespassed, ghost state not populated, or no redraw'
      nfail = nfail + 1
   else
      print '(a)', '(a) dispatch chain: PASS (reg row -> code 9 -> samp_run dispatch; sta untouched; ghosts drawn and redrawn)'
   end if
   call state_destroy(sta)

   ! ---------- (f)-(l) loud-stop criteria (child processes) ----------
   call expect_loudstop('bad-temp',   'glo_target_init', 't_glo',     '(f) non-positive bath temperature loud-stop')
   call expect_loudstop('bad-gamma',  'glo_target_init', 'gamma_glo', '(g) zero damping loud-stop')
   call expect_loudstop('bad-w',      'glo_target_init', 'w_ghost',   '(h) zero ghost frequency loud-stop')
   call expect_loudstop('bad-dt',     'glo_target_init', 'step',      '(i) zero propagator step loud-stop')
   call expect_loudstop('bad-bind0',  'glo_target_init', 'scheme',    '(j) no-fragment-binding loud-stop')
   call expect_loudstop('bad-bind2',  'glo_target_init', 'two',       '(k) duplicate-binding loud-stop')
   call expect_loudstop('bad-noinit', 'glo_target_sample', 'glo_target_init', '(l) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_dist_glo_target: ', nfail, ' failure(s) ==='
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
      cmd = 'build'//sep//'check_dist_glo_target.exe '//trim(cmode)// &
            ' 2>build'//sep//'glo_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/glo_err.txt', status='old', action='read', iostat=ios)
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

end program check_dist_glo_target
