!=====================================================================
! check_propagator.f90 - criteria program for the propagator framework + the verlet
!                        member (criteria defined FIRST in this header; this program
!                        executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files out
! of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 globals/config_roster.f90 &
!            globals/state.f90 facade/elec_ws.f90 facade/elec_facade.f90 &
!            facade/force_facade.f90 dynamics/propagator.f90 methods/reg_prop.f90 &
!            ../scripts/check_propagator.f90 -o build/check_propagator.exe
! Posture: the force chain is assembled CHECK-SIDE - a local analytic harmonic test
!   container is bound through container_bind_pes (a degenerate no-op method row
!   registered through method_reg, then container_init), so every criterion below
!   runs on REAL forces; the default binding posture of the tree is untouched (the
!   driver binds nothing). The harmonic test potential (per-dof distinct stiffness
!   and centers over 2 atoms of distinct mass) is a local procedure of this check -
!   it never enters src/.
! Unit notes: q [Ang]; p [amu*Ang/(10 fs)]; t and dt [10 fs]; mass [amu] (roster);
!   facade force [internal energy units per Ang]. The eV->kcal/mol conversion
!   literal lives only in force_facade - this check NEVER writes it: the composed
!   conversion constant is measured once by a calibrated force_eval call (cal below)
!   and reused everywhere.
! Loud-stop verification approach (same scheme as the sibling checks): every
!   loud-stop criterion runs as a CHILD PROCESS - the check re-spawns itself with a
!   mode argument (execute_command_line, stderr redirected); the child drives the
!   real stop path and must print a naming message on stderr and exit nonzero.
! Criteria:
!   (a) framework loud-stops (children): (a1) duplicate integrator id;
!       (a2) registry full (max_integ exceeded); (a3) dispatch with an unregistered
!       id names the id; (a4) illegal verlet order (= 3) names the order
!   (b) velocity-Verlet energy-error scaling: harmonic oscillator, max |E(t)-E(0)|
!       over the run drops ~4x per dt halving (dt 0.02 -> 0.01, same T = 12.56;
!       ratio window [3.5, 4.5] - second-order integrator, dt^2 energy error)
!   (c) velocity-Verlet time reversibility: 200 steps, p -> -p, 200 steps: q returns
!       to the start and p to -p0 (absolute tolerance 1e-11 on coordinates of
!       magnitude ~0.3 and momenta ~0.2 - roundoff-level for 400 steps)
!   (d) endpoint-force contract: after prop_step returns (an order=1 step; an
!       order=2 boot step; an order=2 THIRD step - genuine Beeman differencing),
!       sta%f equals the analytic force of the test potential at the RETURNED q,
!       component-wise (rel 1e-13 - the same arithmetic, the measured conversion)
!   (e) Beeman position accuracy: third-order scaling - the q(T) error against the
!       analytic solution drops ~8x per dt halving (dt 0.08/0.04/0.02, T = 12.56,
!       N = 157/314/628; ratio window [6.0, 10.5] per halving; anchor = the closed-
!       form harmonic solution with the measured-conversion frequencies)
!   (f) Beeman first-step boot + force staging (registered semantics: no usable
!       history -> the step runs the velocity-Verlet sequence and stages its
!       left-end force as the next step's history):
!       (f1) the boot step is IDENTICAL to an order=1 step from the same state
!            (q, p, f compared exactly - both run the same shared sequence);
!       (f2) from step 2 on the force differencing is correct - a 5-step order=2
!            run at dt 0.04 beats the order=1 position error by more than 2x
!            (third- vs second-order positions; a mis-staged history degrades the
!            error toward or below the order=1 level)
!   (g) force path through the facade (bound slot): force_eval returns the analytic
!       force component-wise at two geometries (proportionality to the gradient -
!       the conversion constant cancels), and e_tot decomposes as t_kin(p) + v_pot(q)
!       with the derived conversion consistent across both geometries and momentum
!       sets (rel 1e-12; e_conv from consts, never the facade literal)
!   (h) unbound class-A slot loud-stop (child): force_eval called before any
!       container_bind_pes names the slot on stderr
!=====================================================================
program check_propagator
   use consts,        only: e_conv
   use config,        only: prop_cfg => propagator
   use config_roster, only: roster
   use state,         only: state_t, state_create, state_destroy
   use elec_ws,       only: ws_t, r_scalar
   use elec_facade,   only: method_reg
   use force_facade,  only: force_eval, container_bind_pes, container_init
   use propagator,    only: prop_step, prop_reg, verlet_step
   use reg_prop,      only: reg_prop_all
   implicit none
   integer :: nfail = 0, nfail0
   character(len=64) :: mode
   type(state_t) :: sta, sa, sb
   real(8) :: e, ea, eb, cal, tk
   real(8) :: om(6), fa(6), fb(6), qa(6), u2(6)
   real(8) :: errc, errf, r1, r2, errh(3), errb, errv, tmax
   integer :: i, n
   ! harmonic test-potential parameters (per-dof, 6 dofs / 2 atoms of mass 1 and 4 amu)
   real(8), parameter :: k_eV(6)  = [ 1.00d0, 1.10d0, 0.90d0, 1.20d0, 0.80d0, 1.05d0 ] ! stiffness [eV/A^2]
   real(8), parameter :: q_ctr(6) = [ 0.10d0,-0.20d0, 0.30d0,-0.15d0, 0.25d0,-0.05d0 ] ! centers [A]
   real(8), parameter :: u0(6)    = [ 0.32d0,-0.24d0, 0.28d0,-0.20d0, 0.30d0,-0.16d0 ] ! turning-point offsets [A]
   real(8), parameter :: p0c(6)   = [ 0.21d0,-0.14d0, 0.10d0,-0.18d0, 0.12d0, 0.26d0 ] ! moving-run momenta
   ! criterion-(g) geometries / momenta (distinct from every run above)
   real(8), parameter :: uA(6)    = [ 0.30d0,-0.21d0, 0.18d0,-0.12d0, 0.24d0,-0.09d0 ]
   real(8), parameter :: uB(6)    = [-0.15d0, 0.27d0,-0.06d0, 0.33d0,-0.18d0, 0.12d0 ]
   real(8), parameter :: pB(6)    = [ 1.1d0,-0.7d0, 0.5d0,-0.9d0, 0.6d0, 1.3d0 ]
   tmax = 12.56d0

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('c-dup')                       ! (a1): the verlet row registered twice
         call reg_prop_all()
         call prop_reg(1, verlet_step)
      case ('c-full')                      ! (a2): 1 shipped row + 7 fills = 8, then one more
         call reg_prop_all()
         do i = 1, 7
            call prop_reg(10 + i, verlet_step)
         end do
         call prop_reg(20, verlet_step)
      case ('c-unreg')                     ! (a3): dispatch with an unregistered id
         call reg_prop_all()
         call state_create(sa, 2)
         prop_cfg%id = 9
         prop_cfg%order = 1
         call prop_step(sa)
      case ('c-order')                     ! (a4): illegal verlet order
         call reg_prop_all()
         call state_create(sa, 2)
         prop_cfg%id = 1
         prop_cfg%order = 3
         call prop_step(sa)
      case ('c-unbound')                   ! (h): force_eval before any binding
         call state_create(sa, 2)
         call force_eval(sa, e)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (returned - bug)'
      stop 0
   end if

   ! ---------- parent: assemble the force chain check-side ----------
   allocate (roster%mass(2))
   roster%mass = [ 1.0d0, 4.0d0 ]
   call state_create(sta, 2)
   call state_create(sa, 2)
   call state_create(sb, 2)

   ! degenerate no-op method row + the local harmonic container + the single assembly
   ! entry (fills the method slot, seeds the statistics holder), then the shipped
   ! verlet registration row - the family every criterion below dispatches through
   call method_reg('adiabatic', '', noop_prop, noop_match)
   call container_bind_pes(harm_pes)
   call container_init()
   call reg_prop_all()
   prop_cfg%id = 1

   ! ---------- calibration: measure the composed conversion constant ----------
   ! one force_eval call at a known single-dof offset (p = 0): f = -k*u*cal, with cal
   ! the eV->internal product the facade applies - the literal lives only in
   ! force_facade, never in this check; the per-dof frequencies follow
   ! omega^2 = k*cal/m
   sa%q = q_ctr
   sa%q(1) = q_ctr(1) + 0.30d0
   sa%p = 0.0d0
   call force_eval(sa, e)
   cal = -sa%f(1)/(k_eV(1)*0.30d0)
   if (cal <= 0.0d0) then
      print '(a)', 'FAIL (cal) measured conversion constant non-positive: '//trim(d2s(cal))
      nfail = nfail + 1
   end if
   do i = 1, 6
      om(i) = sqrt(k_eV(i)*cal/roster%mass((i + 2)/3))
   end do

   ! ---------- (a)+(h) framework loud-stops (children) ----------
   call expect_loudstop('c-dup',     'duplicate',     'prop_reg',    '(a1) duplicate-id loud-stop')
   call expect_loudstop('c-full',    'full',          'prop_reg',    '(a2) registry-full loud-stop')
   call expect_loudstop('c-unreg',   'id 9',          'prop_step',   '(a3) unregistered-id loud-stop')
   call expect_loudstop('c-order',   'order',         'verlet_step', '(a4) illegal-order loud-stop')
   call expect_loudstop('c-unbound', 'slot',          'force_eval',  '(h) unbound-slot loud-stop')

   ! ---------- (b) velocity-Verlet energy-error scaling ----------
   errc = max_e_err(0.02d0, nint(tmax/0.02d0))
   errf = max_e_err(0.01d0, nint(tmax/0.01d0))
   nfail0 = nfail
   r1 = 0.0d0
   if (errf <= 0.0d0) then
      call fail('(b) the fine-step energy error vanished - scaling ratio undefined')
   else
      r1 = errc/errf
      if (r1 < 3.5d0 .or. r1 > 4.5d0) then
         call fail('(b) energy-error scaling out of window: ratio '//trim(d2s(r1))// &
                   ' (want [3.5, 4.5]; coarse '//trim(d2s(errc))//' fine '//trim(d2s(errf))//')')
      end if
   end if
   if (nfail == nfail0) print '(a)', '(b) VV energy scaling: PASS (max|E-E0| ratio '//trim(d2s(r1))// &
      ' per dt halving, dt 0.02->0.01, T = 12.56)'

   ! ---------- (c) velocity-Verlet time reversibility ----------
   prop_cfg%order = 1
   prop_cfg%dt = 0.01d0
   sta%q = q_ctr + u0
   sta%p = p0c
   sta%t = 0.0d0
   do i = 1, 200
      call prop_step(sta)
   end do
   sta%p = -sta%p
   do i = 1, 200
      call prop_step(sta)
   end do
   nfail0 = nfail
   if (maxval(abs(sta%q - (q_ctr + u0))) > 1.0d-11) then
      call fail('(c) reversibility: q did not return to the start (max dev '// &
                trim(d2s(maxval(abs(sta%q - (q_ctr + u0)))))//')')
   end if
   if (maxval(abs(sta%p + p0c)) > 1.0d-11) then
      call fail('(c) reversibility: p did not return to -p0 (max dev '// &
                trim(d2s(maxval(abs(sta%p + p0c))))//')')
   end if
   if (nfail == nfail0) print '(a)', '(c) VV reversibility: PASS (400 steps round trip, q/p back to 1e-11)'

   ! ---------- (d) endpoint-force contract ----------
   nfail0 = nfail
   prop_cfg%order = 1
   prop_cfg%dt = 0.02d0
   sta%q = q_ctr + u0
   sta%p = p0c
   sta%t = 0.0d0
   call prop_step(sta)
   call check_endpoint_force('(d) order=1 step')
   prop_cfg%order = 2
   sta%q = q_ctr + u0
   sta%p = p0c
   sta%t = 0.0d0
   call prop_step(sta)
   call check_endpoint_force('(d) order=2 boot step')
   call prop_step(sta)                       ! genuine Beeman (history staged)
   call prop_step(sta)                       ! genuine Beeman
   call check_endpoint_force('(d) order=2 third step')
   if (nfail == nfail0) print '(a)', '(d) endpoint force: PASS (sta%f = analytic force at the returned q,'// &
      ' order=1 / boot / third Beeman step, rel 1e-13)'

   ! ---------- (e) Beeman third-order position scaling ----------
   errh(1) = run_q_err(2, 0.08d0, nint(tmax/0.08d0))
   errh(2) = run_q_err(2, 0.04d0, nint(tmax/0.04d0))
   errh(3) = run_q_err(2, 0.02d0, nint(tmax/0.02d0))
   nfail0 = nfail
   r1 = 0.0d0
   r2 = 0.0d0
   if (errh(2) <= 0.0d0 .or. errh(3) <= 0.0d0) then
      call fail('(e) a Beeman position error vanished - scaling ratio undefined')
   else
      r1 = errh(1)/errh(2)
      r2 = errh(2)/errh(3)
      if (r1 < 6.0d0 .or. r1 > 10.5d0 .or. r2 < 6.0d0 .or. r2 > 10.5d0) then
         call fail('(e) Beeman scaling out of window: ratios '//trim(d2s(r1))//' / '//trim(d2s(r2))// &
                   ' (want [6.0, 10.5]; errors '//trim(d2s(errh(1)))//' / '//trim(d2s(errh(2)))//' / '// &
                   trim(d2s(errh(3)))//')')
      end if
   end if
   if (nfail == nfail0) print '(a)', '(e) Beeman position scaling: PASS (error ratios '//trim(d2s(r1))//' / '// &
      trim(d2s(r2))//' per dt halving, dt 0.08/0.04/0.02, T = 12.56, vs analytic)'

   ! ---------- (f) Beeman boot + force staging ----------
   ! (f1) the boot step is identical to an order=1 step from the same state
   prop_cfg%dt = 0.02d0
   sa%q = q_ctr + u0
   sa%p = p0c
   sa%t = 0.0d0
   prop_cfg%order = 2
   call prop_step(sa)
   sb%q = q_ctr + u0
   sb%p = p0c
   sb%t = 0.0d0
   prop_cfg%order = 1
   call prop_step(sb)
   nfail0 = nfail
   if (any(sa%q /= sb%q) .or. any(sa%p /= sb%p) .or. any(sa%f /= sb%f)) then
      call fail('(f1) the boot step differs from an order=1 step (registered semantics: VV start)')
   end if
   if (nfail == nfail0) print '(a)', '(f1) boot step = VV step: PASS (q/p/f identical from the same state)'
   ! (f2) from step 2 on the differencing is correct: the 5-step Beeman error beats VV
   errb = run_q_err(2, 0.04d0, 5)
   errv = run_q_err(1, 0.04d0, 5)
   nfail0 = nfail
   if (errv <= 0.0d0) then
      call fail('(f2) the order=1 position error vanished - comparison undefined')
   else if (errb >= errv/2.0d0) then
      call fail('(f2) Beeman 5-step error not below half the order=1 error: '//trim(d2s(errb))//' vs '// &
                trim(d2s(errv))//' (mis-staged force history?)')
   end if
   if (nfail == nfail0) print '(a)', '(f2) Beeman staging: PASS (5-step error '//trim(d2s(errb))//' < half the '// &
      'order=1 error '//trim(d2s(errv))//', dt 0.04)'

   ! ---------- (g) force path through the facade ----------
   nfail0 = nfail
   sa%q = q_ctr + uA
   sa%p = 0.0d0
   call force_eval(sa, ea)
   sb%q = q_ctr + uB
   sb%p = pB
   call force_eval(sb, eb)
   fa = -(k_eV*(sa%q - q_ctr))*cal
   fb = -(k_eV*(sb%q - q_ctr))*cal
   do i = 1, 6
      if (abs(sa%f(i) - fa(i)) > 1.0d-13*max(abs(fa(i)), 1.0d0)) then
         call fail('(g) analytic-force mismatch, geometry A dof '//trim(i2s(i))//': '// &
                   trim(d2s(sa%f(i)))//' vs '//trim(d2s(fa(i))))
      end if
      if (abs(sb%f(i) - fb(i)) > 1.0d-13*max(abs(fb(i)), 1.0d0)) then
         call fail('(g) analytic-force mismatch, geometry B dof '//trim(i2s(i))//': '// &
                   trim(d2s(sb%f(i)))//' vs '//trim(d2s(fb(i))))
      end if
   end do
   ! decomposition e_tot = t_kin(p) + v_pot(q) [kcal/mol], conversion derived, never literal
   u2 = uA*uA
   tk = 0.0d0
   do i = 1, 2
      tk = tk + dot_product(sa%p(3*i - 2:3*i), sa%p(3*i - 2:3*i))/(2.0d0*roster%mass(i))
   end do
   tk = tk/e_conv
   if (abs(ea - (tk + 0.5d0*dot_product(k_eV, u2)*cal/e_conv)) > 1.0d-12*max(abs(ea), 1.0d0)) then
      call fail('(g) e_tot decomposition mismatch, geometry A: '//trim(d2s(ea))//' vs '//trim(d2s(tk)))
   end if
   u2 = uB*uB
   tk = 0.0d0
   do i = 1, 2
      tk = tk + dot_product(sb%p(3*i - 2:3*i), sb%p(3*i - 2:3*i))/(2.0d0*roster%mass(i))
   end do
   tk = tk/e_conv
   if (abs(eb - (tk + 0.5d0*dot_product(k_eV, u2)*cal/e_conv)) > 1.0d-12*max(abs(eb), 1.0d0)) then
      call fail('(g) e_tot decomposition mismatch, geometry B: '//trim(d2s(eb))//' vs '//trim(d2s(tk)))
   end if
   if (nfail == nfail0) print '(a)', '(g) facade force path: PASS (bound force_eval = analytic gradient at two'// &
      ' geometries; e_tot = t_kin(p) + v_pot(q) via the derived conversion)'

   call state_destroy(sta)
   call state_destroy(sa)
   call state_destroy(sb)
   print '(a,i0,a)', '=== check_propagator: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! the local harmonic test container: per-dof distinct stiffness/centers, the
   ! gradient and energy in eV (the member unit contract); fills the scalar channel
   subroutine harm_pes(q, mass, ws)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused - a pure PES)
      type(ws_t), intent(inout) :: ws  ! workspace (scalar channel out)
      integer :: n
      n = size(q)
      if (allocated(ws%g_add)) then
         if (size(ws%g_add) /= n) deallocate (ws%g_add)
      end if
      if (.not. allocated(ws%g_add)) allocate (ws%g_add(n))
      ws%g_add = k_eV*(q - q_ctr)                            ! dV/dq [eV/A]
      ws%v_add = 0.5d0*dot_product(k_eV, (q - q_ctr)*(q - q_ctr))  ! V [eV]
      ws%repr = r_scalar
   end subroutine harm_pes

   ! degenerate no-op method members (the all-noop package this check registers so
   ! container_init resolves the default method name)
   subroutine noop_prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine noop_prop

   subroutine noop_match(sta)
      type(state_t), intent(inout) :: sta
   end subroutine noop_match

   ! max |E(t) - E(0)| over an order=1 run from the turning-point ICs (energy
   ! sampled after every step through the same force entry the member uses)
   real(8) function max_e_err(dt, n)
      real(8), intent(in) :: dt
      integer, intent(in) :: n
      integer :: i
      real(8) :: e0l, el
      prop_cfg%order = 1
      prop_cfg%dt = dt
      sta%q = q_ctr + u0
      sta%p = 0.0d0
      sta%t = 0.0d0
      call force_eval(sta, e0l)
      max_e_err = 0.0d0
      do i = 1, n
         call prop_step(sta)
         call force_eval(sta, el)
         max_e_err = max(max_e_err, abs(el - e0l))
      end do
   end function max_e_err

   ! max |q(T) - q_analytic(T)| for an n-step run at dt from the turning-point ICs
   ! (each dof an independent harmonic oscillator, omega from the measured conversion)
   real(8) function run_q_err(order, dt, n)
      integer, intent(in) :: order
      real(8), intent(in) :: dt
      integer, intent(in) :: n
      integer :: i
      prop_cfg%order = order
      prop_cfg%dt = dt
      sta%q = q_ctr + u0
      sta%p = 0.0d0
      sta%t = 0.0d0
      do i = 1, n
         call prop_step(sta)
      end do
      qa = q_ctr + u0*cos(om*sta%t)
      run_q_err = maxval(abs(sta%q - qa))
   end function run_q_err

   ! assert sta%f equals the analytic force at the RETURNED sta%q (endpoint contract)
   subroutine check_endpoint_force(tag)
      character(len=*), intent(in) :: tag
      integer :: i
      fa = -(k_eV*(sta%q - q_ctr))*cal
      do i = 1, 6
         if (abs(sta%f(i) - fa(i)) > 1.0d-13*max(abs(fa(i)), 1.0d0)) then
            call fail(tag//': endpoint force mismatch at dof '//trim(i2s(i))//': '// &
                      trim(d2s(sta%f(i)))//' vs analytic '//trim(d2s(fa(i))))
         end if
      end do
   end subroutine check_endpoint_force

   ! fail(msg) - record one criterion failure
   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      print '(a)', 'FAIL '//trim(msg)
      nfail = nfail + 1
   end subroutine fail

   ! spawn the child for one bad case and assert the loud-stop: nonzero exit, stderr
   ! names the authority and carries the expected key
   subroutine expect_loudstop(cmode, who, key, label)
      character(len=*), intent(in) :: cmode, who, key, label
      character(len=512) :: cmd, ln
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_propagator.exe '//trim(cmode)//' 2>build'//sep//'prop_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/prop_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) ln
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(ln)
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

   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

   ! small int/real to string helpers for failure messages
   function i2s(i)
      integer, intent(in) :: i
      character(len=12) :: i2s
      write (i2s, '(i0)') i
   end function i2s
   function d2s(x)
      real(8), intent(in) :: x
      character(len=24) :: d2s
      write (d2s, '(es12.4)') x
   end function d2s

end program check_propagator
