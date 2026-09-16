!=====================================================================
! check_propagator.f90 - criteria program for the propagator framework + the verlet
!                        member (criteria defined FIRST in this header; this program
!                        executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files out
! of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 globals/config_roster.f90 &
!            globals/state.f90 facade/elec_facade.f90 &
!            utils/linalg.f90 facade/force_facade.f90 dynamics/propagator.f90 methods/reg_prop.f90 &
!            ../scripts/check_propagator.f90 -o build/check_propagator.exe
! Posture: the force chain is assembled CHECK-SIDE - a local analytic test container is
!   bound through container_bind_pes (a degenerate no-op method row registered through
!   method_reg, then container_init), so every criterion below runs on REAL forces; the
!   default binding posture of the tree is untouched (the driver binds nothing). The
!   bound container branches on a check-side selector: harmonic (the default - per-dof
!   distinct stiffness and centers over 2 atoms of distinct mass), kepler two-body
!   attraction, a uniform field, free flight, per-dof Morse, and a per-dof quartic
!   wall (the radau adaptive-response branches); every branch is a local procedure of
!   this check - none of them ever enters src/.
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
!   ---------- symple criteria (the composition member, id = 2) ----------
!   (s-a) coefficient loading (criteria stated on the composed array a(0:2N), N = 5/8/17
!       for orders 4/6/8 - measured through the stepping, the cache being module-private):
!       (s-a1) drift normalization - on free flight q advances by dt*p/m exactly,
!           i.e. the drift coefficients a(1,3,..,2N-1) sum to 1 (rel 1e-12);
!       (s-a2) kick normalization - in a uniform field p advances by dt*f exactly,
!           i.e. the kick coefficients a(0,2,..,2N) sum to 1 (rel 1e-12);
!       (s-a3) palindromic mirror - the single-step momentum-reversal retrace returns
!           to the start (R S R S = I, R: p -> -p, holds iff a(k) = a(2N-k); an
!           asymmetric table retraces with an O(dt) defect; abs 1e-12);
!       (s-a4) bit-level values - q and p bitwise equal after EVERY of 30 steps to an
!           independent transcription of the historical Schlier tables driving the same
!           arithmetic stream (the reference's expression grouping mirrors the member's;
!           any differing coefficient bit perturbs the stream and flips trajectory bits;
!           the order-4 row of the reference carries the SIGN CORRECTION legislated in
!           the propagator header - the historical a(2) = +1/3 set is second order)
!   (s-b) per-order convergence: harmonic position error at fixed T = 25.12 against the
!       analytic solution drops ~2^order per dt halving, orders 4/6/8 each; dt sets per
!       order keep every error well above the accumulated-time roundoff floor (sta%t
!       rounding grows with step count - coarser chains, fewer steps): order 4 dt
!       0.04/0.02/0.01 (window [13, 19]), order 6 dt 0.5024/0.2512/0.1256 (window
!       [40, 80]), order 8 dt 0.5024/0.2512 (window [180, 340])
!   (s-c) symplecticity signature: 50000 steps at dt = 0.2, order 6 (T = 10^4); block
!       maxima of |E-E0| over 10 equal blocks: bounded oscillation - the second half's
!       worst block <= 1.5x the first half's worst; the non-symplectic contrast
!       (classical RK4 on the same bound force, dt and length, check-local) drifts
!       secularly, linearly in step count - its LAST block's max >= 5x its FIRST
!       block's max (a linear drift gives ~10x over ten blocks) AND its overall worst
!       >= 100x the composition's overall worst
!   (s-d) kepler two-body (V = -k/|q1-q2|, k set so GM of the relative motion = 1;
!       a = 1 A, e = 0.5, prograde): the apsidal precession - the mean per-orbit
!       advance of the Runge-Lenz eccentricity-vector angle of the RELATIVE motion,
!       sampled at successive perihelion passages (radial-velocity zero crossings -
!       to - +; per-crossing angles from the step-end state directly, deltas wrapped
!       to (-pi, pi]; the interpolated-atan2 form is NOT accurate enough - its O(dt^2)
!       interpolation bias swamps the order-6/8 signals) over 60 orbits: shrinks
!       monotonically with order 4 -> 6 -> 8 at dt = 0.05 (each at least 8x below the
!       previous) and each order's own advance scales with its dt halving (order 4:
!       dt 0.05 -> 0.025, window [8, 32] x 1; order 6: 0.05 -> 0.025, window
!       [32, 128]; order 8: 0.10 -> 0.05, window [128, 512])
!   (s-e) reversibility: 5 kepler orbits forward + back at order 6, dt = 0.02: q
!       returns to the start and p to -p0 (roundoff-level abs tol)
!   (s-f) endpoint-force contract: after 1 and 3 symple steps (each order, harmonic;
!       order 8 kepler) sta%f equals the analytic force of the branch at the RETURNED
!       q, component-wise (rel 1e-13, as (d))
!   (s-ls) illegal symple order (= 5) loud-stop (child): names the order
!   ---------- radau criteria (the RA15 adaptive member, id = 3) ----------
!   Test potentials beyond the harmonic default: branch 4 = per-dof Morse
!   (V = D(e^{-2au} - 2e^{-au}), D = 5 eV, a = 1.5/A - the anharmonic bound
!   well), branch 5 = per-dof quartic wall (V = K u^4, K = 50 eV/A^4 - the
!   STIFFNESS-GRADIENT branch: curvature 12Ku^2 grows two orders of magnitude
!   across the orbit, so the adaptive length must respond; the implementer's
!   documented second choice after Morse).
!   (r-a) coefficient loading, two layers:
!       (r-a1) node fingerprint - the first-sequence substep geometries observed
!           through the bound container (the check logs the q of every force
!           site) equal the closed form q0 + h*(tt*v) + h^2*tt^2*f1/2 evaluated
!           at an INDEPENDENT transcription of the published interior nodes
!           (uniform branch: f_j = f1 exactly, so B stays 0 and the iteration-1
!           predictors are that closed form; rel 1e-12), the seven nodes are
!           distinct and increasing, and the mapped nodes x = 2h-1 are roots of
!           the Gauss-Radau defining combination P8 + P7 (Legendre recurrence,
!           |P8+P7| < 1e-12 - the value at a non-node midpoint is ~0.27);
!       (r-a2) machinery exactness - ONE sequence on the harmonic branch
!           (dt 0.01, order 8): the 15th-order collocation leaves roundoff only
!           - q and p at the sequence end equal the closed-form map to 1e-13
!           (a wrong cc/d triangle entry anywhere breaks the corrector
!           convergence and the error jumps above the window by orders of
!           magnitude; the w/u close-out is exact by construction)
!   (r-b) Morse energy conservation, long run: 6-dof Morse from rest at u ~ 2 A
!       (bound, per-dof amplitudes), order 8, seed dt 0.02, run to T = 100
!       (~350 well periods, hundreds of adaptive sequences): max |E(t)-E(0)| <
!       1e-7 kcal/mol (the measured bound is ~1e-9 - the window carries two
!       decades of margin; E itself is ~1e2 kcal/mol, so the relative drift is
!       ~1e-11 at ss = 1e-8)
!   (r-c) accuracy-ss subordination: harmonic branch, seed dt 0.05, run to
!       T = 25.12 at order 6 (ss = 1e-6) vs order 7 (ss = 1e-7): the end
!       position error against the closed form improves by at least one decade.
!       MODEL BOUNDS (the window's quantitative anchor): under per-sequence
!       error CONTROL the accumulated error scales as ss^(8/9) - each sequence
!       contributes ~ss and the sequence count grows as ss^(-1/9) - so one ss
!       decade buys ~8x; under pure 15th-order TRUNCATION (sequence length
!       h ~ ss^(1/9), per-sequence error ~ h^15) a decade buys 10^(15/9) ~ 46,
!       50-60 with controller overshoot riding in. The measured ~29 sits
!       between the bounds; the window [4, 40] carries margin below the control
!       bound and stops short of the truncation bound - a pair whose ratio
!       approaches the upper edge is measuring truncation, not control.
!       FLOOR CLEARANCE (asserted alongside the window): both end errors must
!       exceed 100x the anchor-rounding scale 1e-14 (the measured order-8 floor
!       is ~5e-14; the retired 7-vs-8 pair's ratio ~27 was window-compatible
!       yet floor-contaminated - invisible to the window alone, caught here)
!   (r-d) restart test path (verified by COUNTING bound-container calls - one
!       call per force_eval): a mild first sequence (harmonic, dt 0.05,
!       order 8) accepts on the first trial: exactly 44 calls (1 boot head +
!       6 iterations x 7 substeps + 1 endpoint close-out); a harsh one
!       (quartic wall, dt 0.2, order 12) triggers tp/tt <= 1 and restarts:
!       44 + 43 per restart (42 retry substeps + 1 re-head; the measured case
!       restarts twice = 130 calls), and the accepted length lands BELOW the
!       seed (the discard-retry path ran, not a silent accept). SEED CHOICE:
!       0.2, not a larger one - at seed 0.5 the iteration-1 predictor
!       overshoots the quartic wall catastrophically (geometries past 1e50;
!       gfortran's corrector rescues the trial, ifx's never accepts - the
!       diverged-trial guard names it at the hv estimate, the estimate being
!       where the divergence first shows) - the restart path is proven on a
!       seed whose trials stay finite under BOTH compilers, and the divergent
!       cousin is the (r-ls) c-rdiv criterion's domain
!   (r-e) dt_min/dt_max clamping (magnitude bounds; post-clamp step = clamp):
!       upper - free flight (hv = 0, tp = Inf) with dt_max = 0.2, seed 0.05:
!           the increments ramp through the 1.4 growth limiter and from the
!           sixth sequence on EQUAL 0.2 exactly (1e-12);
!       lower - quartic wall, order 15 (ss = 1e-15), dt_min = 0.01: the run's
!           smallest increments EQUAL 0.01 (1e-12) while the UNCLAMPED twin run
!           (same ICs, dt_min = 0) dips strictly below 0.01 - the clamp is
!           active, not vacuous
!   (r-f) endpoint-force contract (as (d)/(s-f)): after the 1st and 3rd
!       sequences (harmonic and Morse branches, order 8) sta%f equals the
!       analytic branch force at the RETURNED q (rel 1e-13)
!   (r-g) dispatch: id = 3 through prop_step advances sta%t by the seed on the
!       first sequence and by a DIFFERENT (adaptive, growth-limited) length on
!       the second - the member is reachable and internally adaptive
!   (r-ls) loud-stops (children): illegal radau order (= 17) names the order;
!       dt = 0 names dt; a positive dt_min > dt_max pair names the pair; a
!       floor-pinned diverged trial names the guard - quartic wall at seed
!       dt = dt_min = 0.7, order 12: the diverged trial's probe clamps back UP
!       to the trial length and the equality path of the accept test would
!       admit the garbage with no naming (the retry cap only sees FAILED
!       comparisons on this path). The child asserts the shared 'diverged'
!       name of the TERMINATION GUARD family's two arms and must fire at ANY
!       optimization level: the estimate arm (hv non-finite / > 1e30, named
!       at the hv computation) catches the divergence early where rounding
!       lets it see the NaN; the close-out STATE arm (per-component q/v
!       sanity - no reduction, order-independent) is the backstop for the
!       measured -O0 instance where the hv max() keeps the finite side and
!       the garbage state carries NaN behind moderate maxval magnitudes
!   Cross-sequence carry: every multi-sequence criterion above (r-b/r-c/r-e/r-f
!       step-3/r-g) runs consecutive sequences through one member instance -
!       the B-carry path (2-iteration later sequences) is exercised throughout
!=====================================================================
program check_propagator
   use consts,        only: e_conv, two_pi
   use config,        only: prop_cfg => propagator
   use config_roster, only: roster
   use state,         only: state_t, state_create, state_destroy
   use elec_facade,   only: method_reg
   use force_facade,  only: force_eval, container_bind_pes, container_init
   use propagator,    only: prop_step, prop_reg, verlet_step
   use reg_prop,      only: reg_prop_all
   implicit none
   integer :: nfail = 0, nfail0
   character(len=64) :: mode
   type(state_t) :: sta, sa, sb, sr
   real(8) :: e, ea, eb, cal, tk
   real(8) :: om(6), fa(6), fb(6), qa(6), u2(6)
   real(8) :: errc, errf, r1, r2, errh(3), errb, errv, tmax
   integer :: i, n
   ! ---------- symple-section state ----------
   integer :: iord
   real(8) :: dqv(3), fu(6)
   real(8) :: serr(3), srat(2)
   real(8) :: bmax(10), rkmax(10)
   real(8) :: adv(3), advh(3), qk0(6), pk0(6)
   ! ---------- radau-section state ----------
   integer :: fcount = 0                  ! bound-container call count (= force_eval count;
                                         ! one container call per evaluation, facade step 3)
   logical :: rlog_on = .false.           ! substep-geometry logging flag (test_pes appends)
   integer :: rlog_n = 0                  ! logged force-site count
   real(8), allocatable :: rlog_q(:, :)   ! logged force-site geometries (6, site)
   integer :: jj
   real(8) :: tprev, inc1, inc2, hmin_c, hmin_u, mor_drift
   real(8) :: rerr(2), pa(6), vk(6)
   ! the independent transcription of the published RA15 interior nodes (the seven
   ! nonzero Gauss-Radau points of the 8-point left rule) - criterion (r-a1) anchors
   ! the member's loaded values through the OBSERVED geometries against THESE
   real(8), parameter :: rad_node(7) = [ 0.05626256053692215d0, 0.18024069173689236d0, &
                                        0.35262471711316964d0, 0.54715362633055538d0, &
                                        0.73421017721541053d0, 0.88532094683909577d0, &
                                        0.97752061356128750d0 ]
   ! Morse / quartic-wall branch parameters (per-dof; u = q - q_ctr)
   real(8), parameter :: D_eV = 5.0d0     ! Morse well depth [eV]
   real(8), parameter :: a_mor = 1.5d0    ! Morse range parameter [1/A]
   real(8), parameter :: K_wall = 50.0d0  ! quartic stiffness [eV/A^4]
   real(8), parameter :: u_mor(6) = [ 2.0d0, 1.7d0, 2.3d0, 1.5d0, 2.1d0, 1.8d0 ] ! Morse ICs [A]
   real(8), parameter :: u_wall(6) = [ 1.2d0, 1.0d0, 1.4d0, 0.8d0, 1.3d0, 1.1d0 ] ! wall ICs [A]
   ! bound-container branch selector + branch parameters: 1 harmonic (the default, every
   ! pre-existing criterion) / 2 kepler / 3 uniform field / otherwise free flight
   integer :: pot_sel = 1
   real(8) :: kep_k = 0.0d0               ! kepler coupling [eV*A] (V = -kep_k/|q1-q2|; set
                                          ! after calibration so GM of the relative motion
                                          ! = gm_want in internal units)
   real(8), parameter :: gm_want = 1.0d0  ! target GM of the relative motion [internal]
   real(8), parameter :: kep_a = 1.0d0    ! orbit semi-major axis [A]
   real(8), parameter :: kep_ecc = 0.5d0  ! orbit eccentricity [-]
   real(8), parameter :: f_uni(6) = [ 0.30d0,-0.20d0, 0.25d0,-0.15d0, 0.10d0,-0.22d0 ] ! uniform [eV/A]
   ! the reference composition of (s-a4): coefficient values from an INDEPENDENT
   ! transcription of the historical tables (the assignment forms mirror the historical
   ! block, mirror copies included); ref_n stages, ref_a(0:2*ref_n)
   real(8) :: ref_a(0:34)
   integer :: ref_n = 0
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
      case ('c-full')                      ! (a2): 3 shipped rows + 5 fills = 8, then one more
         call reg_prop_all()
         do i = 1, 5
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
      case ('c-sorder')                    ! (s-ls): illegal symple order (= 5)
         call reg_prop_all()
         call state_create(sa, 2)
         prop_cfg%id = 2
         prop_cfg%order = 5
         call prop_step(sa)
      case ('c-rorder')                    ! (r-ls): illegal radau order (= 17, domain 1..15)
         call reg_prop_all()
         call state_create(sa, 2)
         prop_cfg%id = 3
         prop_cfg%order = 17
         prop_cfg%dt = 0.05d0
         call prop_step(sa)
      case ('c-rdt0')                      ! (r-ls): radau with dt = 0 (no seed/direction)
         call reg_prop_all()
         call state_create(sa, 2)
         prop_cfg%id = 3
         prop_cfg%order = 8
         prop_cfg%dt = 0.0d0
         call prop_step(sa)
      case ('c-rmm')                       ! (r-ls): contradictory positive clamp pair
         call reg_prop_all()               ! (dt_min > dt_max - named at boot)
         call state_create(sa, 2)
         prop_cfg%id = 3
         prop_cfg%order = 8
         prop_cfg%dt = 0.05d0
         prop_cfg%dt_min = 0.5d0
         prop_cfg%dt_max = 0.1d0
         call prop_step(sa)
      case ('c-rdiv')                      ! (r-ls): floor-pinned diverged trial - the
         allocate (roster%mass(2))         ! quartic branch must be bound (the trial runs
         roster%mass = [ 1.0d0, 4.0d0 ]    ! real forces before the estimate fires)
         call state_create(sa, 2)
         call method_reg('adiabatic', noop_prop, noop_match)
         call container_bind_pes(test_pes)
         call container_init()
         call reg_prop_all()
         pot_sel = 5
         prop_cfg%id = 3
         prop_cfg%order = 12
         prop_cfg%dt = 0.7d0
         prop_cfg%dt_min = 0.7d0           ! the floor pins the retry probe to the trial
         sa%q = q_ctr + u_wall             ! length - without the guard the equality
         sa%p = 0.0d0                      ! path of the accept test admits the garbage
         sa%t = 0.0d0
         call prop_step(sa)
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
   call state_create(sr, 2)

   ! degenerate no-op method row + the local branch container + the single assembly
   ! entry (fills the method slot, seeds the statistics holder), then the shipped
   ! registration rows - the family every criterion below dispatches through
   call method_reg('adiabatic', noop_prop, noop_match)
   call container_bind_pes(test_pes)
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
   call expect_loudstop('c-sorder', 'order',          'symple_step', '(s-ls) illegal-symple-order loud-stop')

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

   ! ================= symple section (the composition member, id = 2) =================
   ! every sub-criterion sets its own id/order/dt and container branch; pot_sel is
   ! restored to the harmonic default at the section end
   prop_cfg%id = 2

   ! ---------- (s-a1) drift normalization: free flight ----------
   ! form: with F = 0 the kicks vanish and the composed drift is q1 = q0 +
   ! dt*(a(1)+a(3)+...+a(2N-1))*p0/m - the drift coefficients must sum to 1
   pot_sel = 0
   nfail0 = nfail
   do iord = 4, 8, 2
      sa%q = q_ctr + uA
      sa%p = pB
      sa%t = 0.0d0
      prop_cfg%order = iord
      prop_cfg%dt = 0.02d0
      call prop_step(sa)
      if (any(sa%p /= pB)) then
         call fail('(s-a1) free-flight step changed p (order '//trim(i2s(iord))//')')
      end if
      do i = 1, 2
         dqv = 0.02d0*pB(3*i - 2:3*i)/roster%mass(i)       ! dt*p0/m per atom
         if (maxval(abs(sa%q(3*i - 2:3*i) - (q_ctr(3*i - 2:3*i) + uA(3*i - 2:3*i) + dqv))) > &
             1.0d-12*max(1.0d0, maxval(abs(dqv)))) then
            call fail('(s-a1) drift sum /= 1 (order '//trim(i2s(iord))//', atom '//trim(i2s(i))//')')
         end if
      end do
   end do
   if (nfail == nfail0) print '(a)', '(s-a1) drift normalization: PASS (free flight advances q by'// &
      ' dt*p/m exactly, orders 4/6/8 - drift coefficients sum to 1)'

   ! ---------- (s-a2) kick normalization: uniform field ----------
   ! form: with F = f independent of q every kick adds dt*a(2i)*f and the composed
   ! momentum is p1 = p0 + dt*(a(0)+a(2)+...+a(2N))*f - the kick coefficients must sum
   ! to 1 (fu = the branch's uniform force in internal units, via the measured cal)
   pot_sel = 3
   fu = -(f_uni*cal)
   nfail0 = nfail
   do iord = 4, 8, 2
      sa%q = q_ctr + uB
      sa%p = 0.0d0
      sa%t = 0.0d0
      prop_cfg%order = iord
      prop_cfg%dt = 0.02d0
      call prop_step(sa)
      do i = 1, 6
         if (abs(sa%p(i) - 0.02d0*fu(i)) > 1.0d-12*max(1.0d0, abs(0.02d0*fu(i)))) then
            call fail('(s-a2) kick sum /= 1 (order '//trim(i2s(iord))//', dof '//trim(i2s(i))//')')
         end if
      end do
   end do
   if (nfail == nfail0) print '(a)', '(s-a2) kick normalization: PASS (uniform field advances p by'// &
      ' dt*f exactly, orders 4/6/8 - kick coefficients sum to 1)'

   ! ---------- (s-a3) palindromic mirror: single-step retrace ----------
   ! form: R S R S = I (R: p -> -p) holds to roundoff iff the coefficient sequence is
   ! palindromic, a(k) = a(2N-k) - an asymmetric table retraces with an O(dt) defect
   pot_sel = 1
   nfail0 = nfail
   do iord = 4, 8, 2
      sa%q = q_ctr + u0
      sa%p = p0c
      sa%t = 0.0d0
      prop_cfg%order = iord
      prop_cfg%dt = 0.02d0
      call prop_step(sa)
      sa%p = -sa%p
      call prop_step(sa)
      if (maxval(abs(sa%q - (q_ctr + u0))) > 1.0d-12 .or. &
          maxval(abs(sa%p + p0c)) > 1.0d-12) then
         call fail('(s-a3) single-step retrace failed (order '//trim(i2s(iord))//': q dev '// &
                   trim(d2s(maxval(abs(sa%q - (q_ctr + u0)))))//', p dev '// &
                   trim(d2s(maxval(abs(sa%p + p0c))))//' - asymmetric coefficients?)')
      end if
   end do
   if (nfail == nfail0) print '(a)', '(s-a3) palindromic mirror: PASS (single-step momentum-reversal'// &
      ' retrace back to 1e-12, orders 4/6/8 - a(k) = a(2N-k))'

   ! ---------- (s-a4) bit-level coefficient values: trajectory equality ----------
   ! form: the reference runs the SAME arithmetic stream (identical per-beat expression
   ! grouping) on coefficient values from an independent transcription of the historical
   ! tables - the order-4 row carrying the sign correction legislated in the propagator
   ! header (the historical +1/3 set is second order; both this reference and the member
   ! carry the corrected -1/3); q and p must be bitwise equal after every step - any
   ! differing coefficient bit perturbs a product at the ulp level and flips trajectory
   ! bits
   nfail0 = nfail
   do iord = 4, 8, 2
      call ref_load(iord)
      sa%q = q_ctr + u0
      sa%p = p0c
      sa%t = 0.0d0
      sr%q = q_ctr + u0
      sr%p = p0c
      sr%t = 0.0d0
      prop_cfg%order = iord
      prop_cfg%dt = 0.02d0
      call force_eval(sr, e)               ! seed the reference's carried force at q0
      do i = 1, 30
         call prop_step(sa)
         call ref_step(sr, 0.02d0)
         if (any(sa%q /= sr%q) .or. any(sa%p /= sr%p)) then
            call fail('(s-a4) member vs transcription trajectories differ (order '//trim(i2s(iord))// &
                      ', step '//trim(i2s(i))//' - a coefficient value is not bit-equal)')
            exit
         end if
      end do
   end do
   if (nfail == nfail0) print '(a)', '(s-a4) bit-level coefficients: PASS (q/p bitwise equal to the'// &
      ' independent transcription over 30 steps, orders 4/6/8)'

   ! ---------- (s-b) per-order convergence: harmonic position error O(dt^order) ----------
   ! anchor = the closed-form solution with the measured-conversion frequencies (as (e));
   ! T = 25.12 for every run (comparison at the ACCUMULATED sta%t, so the dt need not
   ! divide T exactly); the per-order dt sets keep the finest error above the roundoff
   ! floor (~2e-15) and the coarsest inside the asymptotic regime. PASS line prints the
   ! measured errors/ratios
   serr(1) = run_q_err(4, 0.04d0, 628, 2)
   serr(2) = run_q_err(4, 0.02d0, 1256, 2)
   serr(3) = run_q_err(4, 0.01d0, 2512, 2)
   nfail0 = nfail
   r1 = 0.0d0
   r2 = 0.0d0
   if (serr(2) <= 0.0d0 .or. serr(3) <= 0.0d0) then
      call fail('(s-b) an order-4 position error vanished - scaling ratio undefined')
   else
      r1 = serr(1)/serr(2)
      r2 = serr(2)/serr(3)
      if (r1 < 13.0d0 .or. r1 > 19.0d0 .or. r2 < 13.0d0 .or. r2 > 19.0d0) then
         call fail('(s-b) order-4 scaling out of window: ratios '//trim(d2s(r1))//' / '//trim(d2s(r2))// &
                   ' (want [13, 19]; errors '//trim(d2s(serr(1)))//' / '//trim(d2s(serr(2)))//' / '// &
                   trim(d2s(serr(3)))//')')
      end if
   end if
   if (nfail == nfail0) print '(a)', '(s-b) order-4 scaling: PASS (error ratios '//trim(d2s(r1))//' / '// &
      trim(d2s(r2))//' per dt halving, dt 0.04/0.02/0.01, T = 25.12; errors '//trim(d2s(serr(1)))// &
      ' / '//trim(d2s(serr(2)))//' / '//trim(d2s(serr(3)))//')'
   serr(1) = run_q_err(6, 0.5024d0, 50, 2)
   serr(2) = run_q_err(6, 0.2512d0, 100, 2)
   serr(3) = run_q_err(6, 0.1256d0, 200, 2)
   nfail0 = nfail
   r1 = 0.0d0
   r2 = 0.0d0
   if (serr(2) <= 0.0d0 .or. serr(3) <= 0.0d0) then
      call fail('(s-b) an order-6 position error vanished - scaling ratio undefined')
   else
      r1 = serr(1)/serr(2)
      r2 = serr(2)/serr(3)
      if (r1 < 40.0d0 .or. r1 > 80.0d0 .or. r2 < 40.0d0 .or. r2 > 80.0d0) then
         call fail('(s-b) order-6 scaling out of window: ratios '//trim(d2s(r1))//' / '//trim(d2s(r2))// &
                   ' (want [40, 80]; errors '//trim(d2s(serr(1)))//' / '//trim(d2s(serr(2)))//' / '// &
                   trim(d2s(serr(3)))//')')
      end if
   end if
   if (nfail == nfail0) print '(a)', '(s-b) order-6 scaling: PASS (error ratios '//trim(d2s(r1))//' / '// &
      trim(d2s(r2))//' per dt halving, dt 0.5024/0.2512/0.1256, T = 25.12; errors '//trim(d2s(serr(1)))// &
      ' / '//trim(d2s(serr(2)))//' / '//trim(d2s(serr(3)))//')'
   serr(1) = run_q_err(8, 0.5024d0, 50, 2)
   serr(2) = run_q_err(8, 0.2512d0, 100, 2)
   nfail0 = nfail
   r1 = 0.0d0
   if (serr(2) <= 0.0d0) then
      call fail('(s-b) the order-8 position error vanished - scaling ratio undefined')
   else
      r1 = serr(1)/serr(2)
      if (r1 < 180.0d0 .or. r1 > 340.0d0) then
         call fail('(s-b) order-8 scaling out of window: ratio '//trim(d2s(r1))// &
                   ' (want [180, 340]; errors '//trim(d2s(serr(1)))//' / '//trim(d2s(serr(2)))//')')
      end if
   end if
   if (nfail == nfail0) print '(a)', '(s-b) order-8 scaling: PASS (error ratio '//trim(d2s(r1))// &
      ' per dt halving, dt 0.5024/0.2512, T = 25.12; errors '//trim(d2s(serr(1)))//' / '//trim(d2s(serr(2)))//')'

   ! ---------- (s-c) symplecticity signature: bounded oscillation, no drift ----------
   ! form: 50000 steps at dt = 0.2, order 6 (T = 10^4); |E-E0| block maxima over 10
   ! equal blocks through the same force entry. Bounded: worst second-half block <= 1.5x
   ! worst first-half block. Contrast: check-local classical RK4 on the same ICs/force/
   ! dt/length drifts secularly - LINEARLY in step count, so the block maxima march
   ! upward block over block: LAST block max >= 5x FIRST block max (a linear drift over
   ! ten equal blocks puts the last block's max ~10x the first's) AND its overall worst
   ! >= 100x the composition's overall worst
   call sym_e_blocks(6, 0.2d0, 50000, bmax)
   call rk4_e_blocks(0.2d0, 50000, rkmax)
   nfail0 = nfail
   if (maxval(bmax(6:10)) > 1.5d0*maxval(bmax(1:5))) then
      call fail('(s-c) composition energy error grows: worst second-half block '//trim(d2s(maxval(bmax(6:10))))// &
                ' vs first-half '//trim(d2s(maxval(bmax(1:5))))//' (bounded oscillation expected)')
   end if
   if (rkmax(10) < 5.0d0*rkmax(1)) then
      call fail('(s-c) the RK4 contrast does not drift: last block max '//trim(d2s(rkmax(10)))// &
                ' vs first block max '//trim(d2s(rkmax(1)))//' (linear secular drift expected)')
   end if
   if (maxval(rkmax) < 100.0d0*maxval(bmax)) then
      call fail('(s-c) contrast too weak: RK4 overall worst '//trim(d2s(maxval(rkmax)))// &
                ' vs composition overall worst '//trim(d2s(maxval(bmax))))
   end if
   if (nfail == nfail0) print '(a)', '(s-c) symplecticity: PASS (50000 steps, dt 0.2, order 6: worst'// &
      ' |E-E0| block first half '//trim(d2s(maxval(bmax(1:5))))//' second half '//trim(d2s(maxval(bmax(6:10))))// &
      ' - bounded oscillation; RK4 last/first block max ratio '//trim(d2s(rkmax(10)/rkmax(1)))// &
      ', drifts to '//trim(d2s(maxval(rkmax)))//')'

   ! ---------- (s-d) kepler apsidal precession converges with order ----------
   ! potential: V = -kep_k/|q1-q2| (the gravitational two-body form), kep_k set so GM
   ! of the relative motion = 1 internal; orbit a = 1 A, e = 0.5, prograde. Instrument
   ! (header note): the Runge-Lenz e-vector angle of the relative motion at successive
   ! perihelion passages, mean wrapped advance per orbit over 60 orbits. adv(i): orders
   ! 4/6/8 at dt 0.05 - the convergence chain; dt-halving pairs: order 4 0.05->0.025,
   ! order 6 0.05->0.025, order 8 0.10->0.05 (its dt 0.05 advance sits ~300x above
   ! the double-precision floor - halving once more would fall through it)
   pot_sel = 2
   kep_k = gm_want*(roster%mass(1)*roster%mass(2)/sum(roster%mass))/cal
   adv(1) = kep_adv(4, 0.05d0, 60)
   adv(2) = kep_adv(6, 0.05d0, 60)
   adv(3) = kep_adv(8, 0.05d0, 60)
   advh(1) = kep_adv(4, 0.025d0, 60)
   advh(2) = kep_adv(6, 0.025d0, 60)
   advh(3) = kep_adv(8, 0.10d0, 60)
   nfail0 = nfail
   if (abs(adv(2)) >= abs(adv(1))/8.0d0 .or. abs(adv(3)) >= abs(adv(2))/8.0d0) then
      call fail('(s-d) precession does not converge with order: advances '// &
                trim(d2s(adv(1)))//' / '//trim(d2s(adv(2)))//' / '//trim(d2s(adv(3)))// &
                ' rad/orbit at dt 0.05 (each order at least 8x below the previous expected)')
   end if
   if (abs(advh(1)) <= 0.0d0 .or. abs(advh(2)) <= 0.0d0 .or. abs(advh(3)) <= 0.0d0) then
      call fail('(s-d) a dt-halved advance vanished - scaling undefined')
   else
      if (abs(adv(1))/abs(advh(1)) < 8.0d0 .or. abs(adv(1))/abs(advh(1)) > 32.0d0) then
         call fail('(s-d) order-4 precession scaling out of window: ratio '// &
                   trim(d2s(abs(adv(1))/abs(advh(1))))//' (want [8, 32] about 2^4)')
      end if
      if (abs(adv(2))/abs(advh(2)) < 32.0d0 .or. abs(adv(2))/abs(advh(2)) > 128.0d0) then
         call fail('(s-d) order-6 precession scaling out of window: ratio '// &
                   trim(d2s(abs(adv(2))/abs(advh(2))))//' (want [32, 128] about 2^6)')
      end if
      if (abs(advh(3))/abs(adv(3)) < 128.0d0 .or. abs(advh(3))/abs(adv(3)) > 512.0d0) then
         call fail('(s-d) order-8 precession scaling out of window: ratio '// &
                   trim(d2s(abs(advh(3))/abs(adv(3))))//' (want [128, 512] about 2^8; pair dt 0.10->0.05)')
      end if
   end if
   if (nfail == nfail0) print '(a)', '(s-d) kepler precession: PASS (advances per orbit '// &
      trim(d2s(adv(1)))//' / '//trim(d2s(adv(2)))//' / '//trim(d2s(adv(3)))//' rad at dt 0.05,'// &
      ' orders 4/6/8; dt-halving ratios '//trim(d2s(abs(adv(1))/abs(advh(1))))//' / '// &
      trim(d2s(abs(adv(2))/abs(advh(2))))//' / '//trim(d2s(abs(advh(3))/abs(adv(3))))//')'

   ! ---------- (s-e) reversibility: kepler round trip ----------
   ! 5 orbits forward at order 6, dt = 0.02 (1571 steps), p -> -p, 5 orbits back: q
   ! returns to the start and p to -p0 (roundoff-level tolerance)
   prop_cfg%order = 6
   prop_cfg%dt = 0.02d0
   call kep_ic(sa)
   qk0 = sa%q
   pk0 = sa%p
   do i = 1, 1571
      call prop_step(sa)
   end do
   sa%p = -sa%p
   do i = 1, 1571
      call prop_step(sa)
   end do
   nfail0 = nfail
   if (maxval(abs(sa%q - qk0)) > 1.0d-10) then
      call fail('(s-e) reversibility: q did not return to the start (max dev '// &
                trim(d2s(maxval(abs(sa%q - qk0))))//')')
   end if
   if (maxval(abs(sa%p + pk0)) > 1.0d-10) then
      call fail('(s-e) reversibility: p did not return to -p0 (max dev '// &
                trim(d2s(maxval(abs(sa%p + pk0))))//')')
   end if
   if (nfail == nfail0) print '(a)', '(s-e) symple reversibility: PASS (3142-step round trip over 5'// &
      ' kepler orbits, q/p back to 1e-10)'

   ! ---------- (s-f) endpoint-force contract ----------
   ! after 1 and 3 steps per order (harmonic) plus order-8 kepler: sta%f equals the
   ! analytic branch force at the RETURNED q (rel 1e-13, as (d))
   pot_sel = 1
   nfail0 = nfail
   do iord = 4, 8, 2
      prop_cfg%order = iord
      prop_cfg%dt = 0.02d0
      sta%q = q_ctr + u0
      sta%p = p0c
      sta%t = 0.0d0
      call prop_step(sta)
      call check_endpoint_force('(s-f) harmonic order='//trim(i2s(iord))//' step 1')
      call prop_step(sta)
      call prop_step(sta)
      call check_endpoint_force('(s-f) harmonic order='//trim(i2s(iord))//' step 3')
   end do
   pot_sel = 2
   prop_cfg%order = 8
   prop_cfg%dt = 0.02d0
   call kep_ic(sta)
   call prop_step(sta)
   call check_kep_endpoint('(s-f) kepler order=8 step 1')
   call prop_step(sta)
   call prop_step(sta)
   call check_kep_endpoint('(s-f) kepler order=8 step 3')
   if (nfail == nfail0) print '(a)', '(s-f) endpoint force: PASS (sta%f = analytic branch force at the'// &
      ' returned q, orders 4/6/8 + kepler order 8, steps 1 and 3, rel 1e-13)'

   pot_sel = 1   ! restore the harmonic default for anything that follows
   ! ================= end symple section =================

   ! ================= radau section (the RA15 adaptive member, id = 3) =================
   ! every sub-criterion sets its own order/dt/clamps and container branch; the clamps
   ! are restored to unclamped at the section end. pot_sel = 1 here (the symple restore)
   prop_cfg%id = 3
   prop_cfg%dt_min = 0.0d0
   prop_cfg%dt_max = 0.0d0

   ! ---------- (r-a1) node fingerprint: observed geometries + the P8+P7 identity ----------
   ! uniform branch: f_j = f1 exactly, so B stays identically zero and every iteration-1
   ! predictor is the closed form yy = q0 + h*(tt*v) + h^2*tt^2*(f1/2). The logged force
   ! sites: site 1 = the boot head at q0, sites 2..8 = the iteration-1 substep geometries
   ! at the seven interior nodes - they must equal the closed form at the INDEPENDENTLY
   ! transcribed published nodes (the member's loaded hh never enters the check side)
   pot_sel = 3
   prop_cfg%order = 8
   prop_cfg%dt = 0.05d0
   allocate (rlog_q(6, 128))
   rlog_n = 0
   rlog_on = .true.
   sa%q = q_ctr + uA
   sa%p = pB
   sa%t = 0.0d0
   call prop_step(sa)
   rlog_on = .false.
   do i = 1, 6
      vk(i) = pB(i)/roster%mass((i + 2)/3)    ! folded IC velocity
      fa(i) = -(f_uni(i)*cal)/roster%mass((i + 2)/3)  ! folded constant acceleration f1
   end do
   nfail0 = nfail
   if (maxval(abs(rlog_q(:, 1) - (q_ctr + uA))) > 1.0d-12) then
      call fail('(r-a1) the boot head site is not the origin geometry')
   end if
   do jj = 1, 7
      do i = 1, 6
         qa(i) = q_ctr(i) + uA(i) + rad_node(jj)*0.05d0*vk(i) &
                 + rad_node(jj)**2*(0.05d0*0.05d0)*fa(i)*0.5d0
         if (abs(rlog_q(i, jj + 1) - qa(i)) > 1.0d-12*max(1.0d0, abs(qa(i)))) then
            call fail('(r-a1) substep geometry at node '//trim(i2s(jj))//' dof '//trim(i2s(i))// &
                      ': '//trim(d2s(rlog_q(i, jj + 1)))//' vs closed form '//trim(d2s(qa(i))))
         end if
      end do
   end do
   do jj = 1, 7                              ! Gauss-Radau identity + distinctness
      if (abs(lp_sum(2.0d0*rad_node(jj) - 1.0d0)) > 1.0d-12) then
         call fail('(r-a1) mapped node '//trim(i2s(jj))//' is not a P8+P7 root: '// &
                   trim(d2s(lp_sum(2.0d0*rad_node(jj) - 1.0d0))))
      end if
      if (jj > 1) then
         if (rad_node(jj) - rad_node(jj - 1) < 1.0d-3) then
            call fail('(r-a1) nodes not distinct/increasing at '//trim(i2s(jj)))
         end if
      end if
   end do
   if (abs(lp_sum(0.0d0)) < 1.0d-2) then     ! contrast: a non-node point is NOT a root
      call fail('(r-a1) P8+P7 near zero at a non-node point - the identity is vacuous')
   end if
   if (nfail == nfail0) print '(a)', '(r-a1) node fingerprint: PASS (first-sequence substep'// &
      ' geometries = the closed form at the published nodes to 1e-12; P8+P7 = 0 at every'// &
      ' mapped node (non-node value '//trim(d2s(lp_sum(0.0d0)))//'); nodes distinct)'

   ! ---------- (r-a2) machinery exactness: one harmonic sequence ----------
   ! one sequence at dt 0.01, order 8: the 15th-order collocation leaves roundoff only -
   ! q AND p at the sequence end equal the closed-form map (a wrong cc/d triangle entry
   ! anywhere breaks the corrector convergence; the exact-moment close-out is structural)
   pot_sel = 1
   prop_cfg%order = 8
   prop_cfg%dt = 0.01d0
   sa%q = q_ctr + u0
   sa%p = 0.0d0
   sa%t = 0.0d0
   call prop_step(sa)
   qa = q_ctr + u0*cos(om*sa%t)
   do i = 1, 6
      pa(i) = -u0(i)*om(i)*sin(om(i)*sa%t)*roster%mass((i + 2)/3)  ! p = m*dq/dt
   end do
   nfail0 = nfail
   if (abs(sa%t - 0.01d0) > 1.0d-14) then
      call fail('(r-a2) the first sequence did not consume the seed (t = '//trim(d2s(sa%t))//')')
   end if
   do i = 1, 6
      if (abs(sa%q(i) - qa(i)) > 1.0d-13*max(1.0d0, abs(qa(i)))) then
         call fail('(r-a2) harmonic one-sequence q mismatch, dof '//trim(i2s(i))//': '// &
                   trim(d2s(sa%q(i)))//' vs closed form '//trim(d2s(qa(i))))
      end if
      if (abs(sa%p(i) - pa(i)) > 1.0d-13*max(1.0d0, abs(pa(i)))) then
         call fail('(r-a2) harmonic one-sequence p mismatch, dof '//trim(i2s(i))//': '// &
                   trim(d2s(sa%p(i)))//' vs closed form '//trim(d2s(pa(i))))
      end if
   end do
   if (nfail == nfail0) print '(a)', '(r-a2) machinery exactness: PASS (one sequence, dt 0.01,'// &
      ' order 8: q and p = closed-form map to 1e-13; worst q dev '// &
      trim(d2s(maxval(abs(sa%q - qa))))//')'

   ! ---------- (r-b) Morse energy conservation (long adaptive run) ----------
   ! bound 6-dof Morse from rest at u ~ 2 A, order 8 (ss = 1e-8), seed 0.02, to T = 100
   ! (~350 well periods, hundreds of sequences through the B-carry); |E(t)-E(0)| is
   ! sampled after every sequence through the same force entry (e_tot channel)
   pot_sel = 4
   prop_cfg%order = 8
   prop_cfg%dt = 0.02d0
   sa%q = q_ctr + u_mor
   sa%p = 0.0d0
   sa%t = 0.0d0
   call force_eval(sa, ea)
   mor_drift = 0.0d0
   do while (sa%t < 100.0d0)
      call prop_step(sa)
      call force_eval(sa, eb)
      mor_drift = max(mor_drift, abs(eb - ea))
   end do
   nfail0 = nfail
   if (mor_drift >= 1.0d-7) then
      call fail('(r-b) Morse energy drift out of tolerance: max|E-E0| = '//trim(d2s(mor_drift))// &
                ' kcal/mol over T = '//trim(d2s(sa%t))//' (want < 1e-7)')
   end if
   if (nfail == nfail0) print '(a)', '(r-b) Morse conservation: PASS (max|E-E0| = '// &
      trim(d2s(mor_drift))//' kcal/mol over T = '//trim(d2s(sa%t))//', order 8, hundreds of'// &
      ' adaptive sequences)'

   ! ---------- (r-c) accuracy-ss subordination (orders 6 vs 7) ----------
   ! harmonic, seed 0.05, run to T = 25.12: ss 1e-6 vs 1e-7 - the end position error
   ! against the closed form improves by about one decade (error ~ ss per unit time).
   ! Orders 6/7, not 7/8: the order-8 run's end error sits at the analytic-anchor
   ! roundoff floor (~5e-14: sta%t accumulation rounding mapped through dq/dt), which
   ! floored the 7-vs-8 ratio at ~27 - one pair lower keeps both errors honest
   pot_sel = 1
   do jj = 1, 2
      prop_cfg%order = 5 + jj                ! 6 then 7
      prop_cfg%dt = 0.05d0
      sa%q = q_ctr + u0
      sa%p = 0.0d0
      sa%t = 0.0d0
      do while (sa%t < 25.12d0)
         call prop_step(sa)
      end do
      qa = q_ctr + u0*cos(om*sa%t)
      rerr(jj) = maxval(abs(sa%q - qa))
   end do
   nfail0 = nfail
   if (rerr(2) <= 0.0d0) then
      call fail('(r-c) the order-7 end error vanished - ratio undefined')
   else if (rerr(1)/rerr(2) < 4.0d0 .or. rerr(1)/rerr(2) > 40.0d0) then
      call fail('(r-c) ss-subordination out of window: err(6)/err(7) = '// &
                trim(d2s(rerr(1)/rerr(2)))//' (want [4, 40]; errors '//trim(d2s(rerr(1)))// &
                ' / '//trim(d2s(rerr(2)))//')')
   end if
   if (rerr(1) <= 100.0d0*1.0d-14 .or. rerr(2) <= 100.0d0*1.0d-14) then
      call fail('(r-c) floor clearance: an end error sits at or below 100x the'// &
                ' anchor-rounding scale 1e-14: '//trim(d2s(rerr(1)))//' / '//trim(d2s(rerr(2)))// &
                ' (floor-contaminated - the ratio window alone cannot see it)')
   end if
   if (nfail == nfail0) print '(a)', '(r-c) ss-subordination: PASS (order 6 vs 7 errors '// &
      trim(d2s(rerr(1)))//' / '//trim(d2s(rerr(2)))//' - about one decade per ss decade, both'// &
      ' > 100x the 1e-14 anchor-rounding scale)'

   ! ---------- (r-d) restart test path (verified by counting container calls) ----------
   ! one container call per force_eval: a mild first sequence (harmonic, dt 0.05,
   ! order 8) accepts on the first trial - exactly 44 calls (1 boot head + 6x7 substeps
   ! + 1 endpoint); a harsh one (quartic wall, dt 0.2, order 12) hits tp/tt <= 1 and
   ! RESTARTS: 44 + 43 per restart (42 retry substeps + 1 re-head), and the accepted
   ! length lands below the seed (the discard-retry path ran, not a silent accept)
   pot_sel = 1
   prop_cfg%order = 8
   prop_cfg%dt = 0.05d0
   sa%q = q_ctr + u0
   sa%p = p0c
   sa%t = 0.0d0
   fcount = 0
   call prop_step(sa)
   nfail0 = nfail
   if (fcount /= 44) then
      call fail('(r-d) mild first sequence: '//trim(i2s(fcount))//' container calls (want 44 - '// &
                '1 boot head + 6x7 substeps + 1 endpoint)')
   end if
   if (nfail == nfail0) print '(a)', '(r-d) mild first sequence: PASS ('//trim(i2s(fcount))// &
      ' container calls = 1 boot head + 6x7 substeps + 1 endpoint - accepted first trial)'
   pot_sel = 5
   prop_cfg%order = 12
   prop_cfg%dt = 0.2d0
   sa%q = q_ctr + u_wall
   sa%p = 0.0d0
   sa%t = 0.0d0
   fcount = 0
   call prop_step(sa)
   nfail0 = nfail
   if (fcount == 44 .or. mod(fcount - 44, 43) /= 0) then
      call fail('(r-d) harsh first sequence: '//trim(i2s(fcount))//' container calls (want'// &
                ' 44 + 43 per restart - the count pattern does not decode)')
   else if (sa%t >= 0.2d0) then
      call fail('(r-d) harsh first sequence accepted the full seed (t = '//trim(d2s(sa%t))// &
                ' - the restart shrink did not run)')
   end if
   if (nfail == nfail0) print '(a)', '(r-d) restart path: PASS ('//trim(i2s(fcount))// &
      ' container calls = 44 + '//trim(i2s((fcount - 44)/43))//' restart(s); accepted length '// &
      trim(d2s(sa%t))//' < seed 0.2)'

   ! ---------- (r-e) dt_min/dt_max clamping (post-clamp step = clamp) ----------
   ! upper: free flight (hv = 0 -> tp = Inf) with dt_max = 0.2, seed 0.05 - the
   ! increments ramp through the 1.4 growth limiter (0.07, 0.098, ...) and from the
   ! sixth sequence on EQUAL the clamp exactly (the cap is respected throughout)
   pot_sel = 0
   prop_cfg%order = 8
   prop_cfg%dt = 0.05d0
   prop_cfg%dt_min = 0.0d0
   prop_cfg%dt_max = 0.2d0
   sa%q = q_ctr + uA
   sa%p = pB
   sa%t = 0.0d0
   tprev = 0.0d0
   nfail0 = nfail
   do i = 1, 10
      call prop_step(sa)
      if (i >= 6) then
         if (abs((sa%t - tprev) - 0.2d0) > 1.0d-12) then
            call fail('(r-e) dt_max clamp: increment '//trim(i2s(i))//' is '// &
                      trim(d2s(sa%t - tprev))//' (want the clamp 0.2)')
         end if
      end if
      tprev = sa%t
   end do
   if (nfail == nfail0) print '(a)', '(r-e) dt_max clamp: PASS (increments 6..10 all = 0.2'// &
      ' exactly after the growth-limited ramp, free flight at seed 0.05)'
   ! lower: quartic wall, order 15 (ss = 1e-15), dt_min = 0.01 - the run's smallest
   ! increments EQUAL the floor while the UNCLAMPED twin dips strictly below it
   pot_sel = 5
   prop_cfg%order = 15
   prop_cfg%dt = 0.02d0
   prop_cfg%dt_min = 0.01d0
   prop_cfg%dt_max = 0.0d0
   sa%q = q_ctr + u_wall
   sa%p = 0.0d0
   sa%t = 0.0d0
   hmin_c = huge(1.0d0)
   tprev = 0.0d0
   do i = 1, 40
      call prop_step(sa)
      hmin_c = min(hmin_c, sa%t - tprev)
      tprev = sa%t
   end do
   prop_cfg%dt_min = 0.0d0                   ! the unclamped twin (same ICs, same order)
   sa%q = q_ctr + u_wall
   sa%p = 0.0d0
   sa%t = 0.0d0
   hmin_u = huge(1.0d0)
   tprev = 0.0d0
   do i = 1, 40
      call prop_step(sa)
      hmin_u = min(hmin_u, sa%t - tprev)
      tprev = sa%t
   end do
   nfail0 = nfail
   if (abs(hmin_c - 0.01d0) > 1.0d-12) then
      call fail('(r-e) dt_min clamp: smallest increment '//trim(d2s(hmin_c))//' (want the'// &
                ' floor 0.01 - the run must touch it)')
   else if (hmin_u >= 0.01d0) then
      call fail('(r-e) dt_min clamp vacuous: the unclamped twin never dips below 0.01 '// &
                '(min '//trim(d2s(hmin_u))//')')
   end if
   if (nfail == nfail0) print '(a)', '(r-e) dt_min clamp: PASS (clamped run floor = 0.01;'// &
      ' unclamped twin dips to '//trim(d2s(hmin_u))//' - the bound is active)'
   prop_cfg%dt_min = 0.0d0
   prop_cfg%dt_max = 0.0d0

   ! ---------- (r-f) endpoint-force contract ----------
   ! after the 1st and 3rd sequences (harmonic and Morse, order 8): sta%f equals the
   ! analytic branch force at the RETURNED q (rel 1e-13, as (d)/(s-f))
   pot_sel = 1
   prop_cfg%order = 8
   prop_cfg%dt = 0.01d0
   nfail0 = nfail
   sta%q = q_ctr + u0
   sta%p = p0c
   sta%t = 0.0d0
   call prop_step(sta)
   call check_endpoint_force('(r-f) harmonic seq 1')
   call prop_step(sta)
   call prop_step(sta)
   call check_endpoint_force('(r-f) harmonic seq 3')
   pot_sel = 4
   sta%q = q_ctr + u_mor
   sta%p = 0.0d0
   sta%t = 0.0d0
   call prop_step(sta)
   call check_morse_endpoint('(r-f) Morse seq 1')
   call prop_step(sta)
   call prop_step(sta)
   call check_morse_endpoint('(r-f) Morse seq 3')
   if (nfail == nfail0) print '(a)', '(r-f) endpoint force: PASS (sta%f = analytic branch'// &
      ' force at the returned q, harmonic + Morse, sequences 1 and 3, rel 1e-13)'

   ! ---------- (r-g) dispatch: id = 3 through prop_step, adaptive increments ----------
   pot_sel = 1
   prop_cfg%order = 8
   prop_cfg%dt = 0.05d0
   sa%q = q_ctr + u0
   sa%p = 0.0d0
   sa%t = 0.0d0
   call prop_step(sa)
   inc1 = sa%t                                ! the seed consumed by sequence 1
   call prop_step(sa)
   inc2 = sa%t - inc1                         ! the adaptive length of sequence 2
   nfail0 = nfail
   if (abs(inc1 - 0.05d0) > 1.0d-14) then
      call fail('(r-g) the first sequence did not consume the seed (inc1 = '//trim(d2s(inc1))//')')
   else if (inc2 <= inc1 + 1.0d-12 .or. inc2 > 1.4d0*inc1 + 1.0d-12) then
      call fail('(r-g) the second increment is not adaptive-growth-limited: inc2 = '// &
                trim(d2s(inc2))//' vs inc1 = '//trim(d2s(inc1))//' (want (inc1, 1.4*inc1])')
   end if
   if (nfail == nfail0) print '(a)', '(r-g) id=3 dispatch: PASS (t advances; seed increment '// &
      trim(d2s(inc1))//', adaptive increment '//trim(d2s(inc2))//' - growth-limited adaptivity)'

   ! ---------- (r-ls) radau loud-stops (children) ----------
   call expect_loudstop('c-rorder', 'order', 'radau_step', '(r-ls) illegal-order loud-stop')
   call expect_loudstop('c-rdt0', 'dt = 0', 'radau_step', '(r-ls) zero-dt loud-stop')
   call expect_loudstop('c-rmm', 'dt_min exceeds dt_max', 'radau_step', &
                        '(r-ls) contradictory-clamp-pair loud-stop')
   call expect_loudstop('c-rdiv', 'diverged', 'radau_step', &
                        '(r-ls) floor-pinned-diverged-trial loud-stop')

   pot_sel = 1   ! restore the harmonic default for anything that follows
   ! ================= end radau section =================

   call state_destroy(sta)
   call state_destroy(sa)
   call state_destroy(sb)
   call state_destroy(sr)
   print '(a,i0,a)', '=== check_propagator: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! the local branch test container (bound once; the check-side selector pot_sel picks
   ! the analytic branch): harmonic per-dof distinct stiffness/centers (the default),
   ! kepler two-body attraction, uniform field, free flight. The kernel and energy in
   ! eV (the E27 delivery contract: g and v fully assigned)
   subroutine test_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused - pure PES branches)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A]
      real(8), intent(out) :: v        ! V [eV]
      integer :: n
      real(8) :: rr(3), rn, ex
      n = size(q)
      fcount = fcount + 1                     ! one container call per force_eval (facade
                                             ! step 2) - the radau restart instrument (r-d)
      if (rlog_on) then                       ! substep-geometry log (r-a1): append the
         if (rlog_n < size(rlog_q, 2)) then   ! q of every force site while armed
            rlog_n = rlog_n + 1
            rlog_q(:, rlog_n) = q
         end if
      end if
      select case (pot_sel)
      case (1)                                 ! harmonic: V = 1/2 sum k (q-qc)^2
         g = k_eV*(q - q_ctr)                            ! dV/dq [eV/A]
         v = 0.5d0*dot_product(k_eV, (q - q_ctr)*(q - q_ctr))  ! V [eV]
      case (2)                                 ! kepler: V = -kep_k/|q1-q2|
         rr = q(1:3) - q(4:6)
         rn = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)
         g(1:3) = kep_k*rr/rn**3               ! dV/dq1 = +k r/|r|^3 (attraction: F = -g*cal)
         g(4:6) = -g(1:3)                      ! dV/dq2 = -dV/dq1 (central two-body form)
         v = -kep_k/rn
      case (3)                                 ! uniform field: V = -f.q
         g = f_uni
         v = -dot_product(f_uni, q)
      case (4)                                 ! Morse: V = D(e^{-2au} - 2e^{-au}) per dof
         v = 0.0d0
         do n = 1, size(q)
            ex = exp(-a_mor*(q(n) - q_ctr(n)))
            g(n) = 2.0d0*a_mor*D_eV*(ex - ex*ex)             ! dV/du [eV/A]
            v = v + D_eV*(ex*ex - 2.0d0*ex)
         end do
      case (5)                                 ! quartic wall: V = K u^4 (curvature 12Ku^2 -
         v = 0.0d0                             ! the stiffness-gradient branch)
         do n = 1, size(q)
            g(n) = 4.0d0*K_wall*(q(n) - q_ctr(n))**3         ! dV/du [eV/A]
            v = v + K_wall*(q(n) - q_ctr(n))**4
         end do
      case default                             ! free flight: V = 0
         g = 0.0d0
         v = 0.0d0
      end select
   end subroutine test_pes

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
   ! (each dof an independent harmonic oscillator, omega from the measured conversion);
   ! id_sel picks the family member (default 1 = verlet; the symple section passes 2)
   real(8) function run_q_err(order, dt, n, id_sel)
      integer, intent(in) :: order
      real(8), intent(in) :: dt
      integer, intent(in) :: n
      integer, intent(in), optional :: id_sel
      integer :: i
      prop_cfg%id = 1
      if (present(id_sel)) prop_cfg%id = id_sel
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

   ! assert sta%f equals the analytic kepler force at the RETURNED sta%q (as above)
   subroutine check_kep_endpoint(tag)
      character(len=*), intent(in) :: tag
      integer :: i
      real(8) :: rr(3), rn3
      rr = sta%q(1:3) - sta%q(4:6)
      rn3 = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)**3
      fa(1:3) = -(kep_k*rr/rn3)*cal           ! internal force on atom 1 (attraction)
      fa(4:6) = -fa(1:3)
      do i = 1, 6
         if (abs(sta%f(i) - fa(i)) > 1.0d-13*max(abs(fa(i)), 1.0d0)) then
            call fail(tag//': kepler endpoint force mismatch at dof '//trim(i2s(i))//': '// &
                      trim(d2s(sta%f(i)))//' vs analytic '//trim(d2s(fa(i))))
         end if
      end do
   end subroutine check_kep_endpoint

   ! lp_sum(x) - the Legendre combination P8(x) + P7(x) by the Bonnet recurrence
   ! (the Gauss-Radau defining identity: the mapped interior nodes 2h-1 are its roots)
   real(8) function lp_sum(x)
      real(8), intent(in) :: x
      real(8) :: pm, pc, nx
      integer :: nn
      pm = 1.0d0                              ! P0
      pc = x                                  ! P1
      do nn = 2, 8
         nx = ((2*nn - 1)*x*pc - (nn - 1)*pm)/nn
         pm = pc
         pc = nx
      end do
      lp_sum = pc + pm                        ! P8 + P7
   end function lp_sum

   ! assert sta%f equals the analytic Morse force at the RETURNED sta%q (branch 4 of
   ! the container: V = D(e^{-2au} - 2e^{-au}), dV/du = 2aD(e^{-au} - e^{-2au}); the
   ! facade applies F = -g*cal - the measured conversion, as check_endpoint_force)
   subroutine check_morse_endpoint(tag)
      character(len=*), intent(in) :: tag
      integer :: i
      real(8) :: ex
      do i = 1, 6
         ex = exp(-a_mor*(sta%q(i) - q_ctr(i)))
         fa(i) = -(2.0d0*a_mor*D_eV*(ex - ex*ex))*cal
         if (abs(sta%f(i) - fa(i)) > 1.0d-13*max(abs(fa(i)), 1.0d0)) then
            call fail(tag//': Morse endpoint force mismatch at dof '//trim(i2s(i))//': '// &
                      trim(d2s(sta%f(i)))//' vs analytic '//trim(d2s(fa(i))))
         end if
      end do
   end subroutine check_morse_endpoint

   ! kepler ICs: relative orbit a = kep_a, e = kep_ecc, GM = gm_want; the center of mass
   ! at rest at the origin, perihelion on +x, prograde in xy (atom positions/momenta
   ! from the COM/relative split: q1 = (m2/M) r, p1 = mu*v_r with P_com = 0)
   subroutine kep_ic(sx)
      type(state_t), intent(inout) :: sx
      real(8) :: m1, m2, mt, mu, rp, vp
      m1 = roster%mass(1)
      m2 = roster%mass(2)
      mt = m1 + m2
      mu = m1*m2/mt
      rp = kep_a*(1.0d0 - kep_ecc)
      vp = sqrt(gm_want*(1.0d0 + kep_ecc)/(kep_a*(1.0d0 - kep_ecc)))
      sx%q = 0.0d0
      sx%p = 0.0d0
      sx%q(1) = (m2/mt)*rp
      sx%q(4) = -(m1/mt)*rp
      sx%p(2) = mu*vp
      sx%p(5) = -mu*vp
      sx%t = 0.0d0
   end subroutine kep_ic

   ! mean apsidal advance per orbit (the precession), via the Runge-Lenz instrument
   ! (header note of (s-d)): perihelion = the inward-to-outward zero crossing of the
   ! radial velocity r.v_rel of the RELATIVE motion; at each crossing the angle of the
   ! eccentricity vector e = (v x h)/GM - r/|r| (the orbit plane is z = 0 exactly - the
   ! planar IC and the central force keep it there bitwise) is taken from the step-end
   ! state directly - no interpolation; deltas between consecutive crossings are
   ! wrapped to (-pi, pi] (each delta is the per-orbit precession, |delta| << pi); the
   ! mean over the crossings of an n_orb-orbit run is returned. vrprev is seeded from
   ! the TRUE initial radial velocity (zero at the perihelion IC) so the first step
   ! cannot fake a crossing
   real(8) function kep_adv(order, dt, n_orb)
      integer, intent(in) :: order
      real(8), intent(in) :: dt
      integer, intent(in) :: n_orb
      integer :: i, n, ncross
      real(8) :: r(3), v(3), vr, vrprev, hz, rn, ex, ey, thA, lastA, dth, acc
      prop_cfg%id = 2
      prop_cfg%order = order
      prop_cfg%dt = dt
      call kep_ic(sta)
      ncross = 0
      acc = 0.0d0
      lastA = 0.0d0
      r = sta%q(1:3) - sta%q(4:6)
      v = sta%p(1:3)/roster%mass(1) - sta%p(4:6)/roster%mass(2)
      vrprev = r(1)*v(1) + r(2)*v(2) + r(3)*v(3)
      n = nint(n_orb*two_pi/dt)
      do i = 1, n
         call prop_step(sta)
         r = sta%q(1:3) - sta%q(4:6)
         v = sta%p(1:3)/roster%mass(1) - sta%p(4:6)/roster%mass(2)
         vr = r(1)*v(1) + r(2)*v(2) + r(3)*v(3)
         if (vrprev < 0.0d0 .and. vr >= 0.0d0) then
            hz = r(1)*v(2) - r(2)*v(1)          ! orbit-plane angular momentum [internal]
            rn = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
            ex = v(2)*hz/gm_want - r(1)/rn      ! e-vector (GM = gm_want = 1 here)
            ey = -v(1)*hz/gm_want - r(2)/rn
            thA = atan2(ey, ex)
            ncross = ncross + 1
            if (ncross >= 2) then
               dth = thA - lastA
               dth = dth - two_pi*anint(dth/two_pi)   ! wrap to (-pi, pi]
               acc = acc + dth
            end if
            lastA = thA
         end if
         vrprev = vr
      end do
      kep_adv = 0.0d0
      if (ncross >= 2) kep_adv = acc/(ncross - 1)
   end function kep_adv

   ! per-block maxima of |E - E0| over 10 equal blocks of a symple run on the harmonic
   ! branch (turning-point ICs; energy sampled after every step through the same force
   ! entry the member uses)
   subroutine sym_e_blocks(order, dt, n, bmax)
      integer, intent(in) :: order
      real(8), intent(in) :: dt
      integer, intent(in) :: n
      real(8), intent(out) :: bmax(10)
      integer :: i, b, bps
      real(8) :: e0l, el
      bps = n/10
      prop_cfg%id = 2
      prop_cfg%order = order
      prop_cfg%dt = dt
      sta%q = q_ctr + u0
      sta%p = 0.0d0
      sta%t = 0.0d0
      call force_eval(sta, e0l)
      bmax = 0.0d0
      do i = 1, n
         call prop_step(sta)
         call force_eval(sta, el)
         b = min((i - 1)/bps + 1, 10)
         bmax(b) = max(bmax(b), abs(el - e0l))
      end do
   end subroutine sym_e_blocks

   ! the non-symplectic contrast of (s-c): classical RK4 on the state (q,p) with the
   ! same bound force entry, ICs, dt and run length (check-local reference; sb is the
   ! stage-geometry scratch - its energy/momenta are never consumed); per-block maxima
   ! of |E - E0| as in sym_e_blocks
   subroutine rk4_e_blocks(dt, n, bmax)
      real(8), intent(in) :: dt
      integer, intent(in) :: n
      real(8), intent(out) :: bmax(10)
      integer :: i, b, bps, k
      real(8) :: e0l, el, h2, h6
      real(8) :: q0(6), p0(6), k1q(6), k1p(6), k2q(6), k2p(6), k3q(6), k3p(6), k4q(6), k4p(6)
      h2 = dt/2.0d0
      h6 = dt/6.0d0
      bps = n/10
      sta%q = q_ctr + u0
      sta%p = 0.0d0
      sta%t = 0.0d0
      call force_eval(sta, e0l)
      bmax = 0.0d0
      do i = 1, n
         q0 = sta%q
         p0 = sta%p
         call force_eval(sta, el)              ! k1 at q0
         k1p = sta%f
         do k = 1, 2
            k1q(3*k - 2:3*k) = p0(3*k - 2:3*k)/roster%mass(k)
         end do
         sb%q = q0 + h2*k1q
         call force_eval(sb, el)               ! k2
         k2p = sb%f
         do k = 1, 2
            k2q(3*k - 2:3*k) = (p0(3*k - 2:3*k) + h2*k1p(3*k - 2:3*k))/roster%mass(k)
         end do
         sb%q = q0 + h2*k2q
         call force_eval(sb, el)               ! k3
         k3p = sb%f
         do k = 1, 2
            k3q(3*k - 2:3*k) = (p0(3*k - 2:3*k) + h2*k2p(3*k - 2:3*k))/roster%mass(k)
         end do
         sb%q = q0 + dt*k3q
         call force_eval(sb, el)               ! k4
         k4p = sb%f
         do k = 1, 2
            k4q(3*k - 2:3*k) = (p0(3*k - 2:3*k) + dt*k3p(3*k - 2:3*k))/roster%mass(k)
         end do
         sta%q = q0 + h6*(k1q + 2.0d0*k2q + 2.0d0*k3q + k4q)
         sta%p = p0 + h6*(k1p + 2.0d0*k2p + 2.0d0*k3p + k4p)
         sta%t = sta%t + dt
         call force_eval(sta, el)
         b = min((i - 1)/bps + 1, 10)
         bmax(b) = max(bmax(b), abs(el - e0l))
      end do
   end subroutine rk4_e_blocks

   ! ref_load(order) - fill the reference coefficient table for one order: an
   ! INDEPENDENT transcription of the historical tables (same values, the historical
   ! assignment structure with explicit mirror copies; order 4 carries the sign
   ! correction legislated in the propagator header) - the member's own literals were
   ! written from the same source; agreement is asserted bitwise by (s-a4)
   subroutine ref_load(order)
      integer, intent(in) :: order
      select case (order)
      case (4)                                 ! historical assignment structure with
         ref_n = 5                             ! the SIGN-CORRECTED a(2) (propagator
         ref_a(0) = 0.5d0                      ! header note: the historical +1/3 set is
         ref_a(1) = -1.0d0/48.0d0              ! second order)
         ref_a(2) = -1.0d0/3.0d0
         ref_a(3) = 3.0d0/8.0d0
         ref_a(4) = -ref_a(2)
         ref_a(5) = 7.0d0/24.0d0
         ref_a(6) = -ref_a(2)
         ref_a(7) = ref_a(3)
         ref_a(8) = ref_a(2)
         ref_a(9) = ref_a(1)
         ref_a(10) = ref_a(0)
      case (6)
         ref_n = 8
         ref_a(0) = 0.06942944346252987735848865824703402d0
         ref_a(1) = 0.2848783771728008405274534645665783d0
         ref_a(2) = -0.1331551983159820940996130995137351d0
         ref_a(3) = 0.3278397575961294541205467836732554d0
         ref_a(4) = 0.0012903891798107897423048174644328d0
         ref_a(5) = -0.3812210427193262947562278437421127d0
         ref_a(6) = 0.4224353656736414269988196238022683d0
         ref_a(7) = 0.268502907950396000108227595502279d0
         ref_a(8) = 0.28d0
         ref_a(9) = ref_a(7)
         ref_a(10) = ref_a(6)
         ref_a(11) = ref_a(5)
         ref_a(12) = ref_a(4)
         ref_a(13) = ref_a(3)
         ref_a(14) = ref_a(2)
         ref_a(15) = ref_a(1)
         ref_a(16) = ref_a(0)
      case (8)
         ref_n = 17
         ref_a(0) = 0.04463795052359022755913999625733590d0
         ref_a(1) = 0.13593258071690959145543264213495574d0
         ref_a(2) = 0.2198844042714707225445535069606167d0
         ref_a(3) = 0.13024946780523828601621193778196846d0
         ref_a(4) = 0.10250365693975069608261241007779814d0
         ref_a(5) = 0.43234521869358547487983257884877035d0
         ref_a(6) = -0.00477482916916881658022489063962934d0
         ref_a(7) = -0.58253476904040845493112837930861213d0
         ref_a(8) = -0.03886264282111817697737420875189743d0
         ref_a(9) = 0.31548728537940479698273603797274199d0
         ref_a(10) = 0.18681583743297155471526153503972746d0
         ref_a(11) = 0.26500275499062083398346002963079872d0
         ref_a(12) = -0.02405084735747361993573587982407554d0
         ref_a(13) = -0.45040492499772251180922896712151891d0
         ref_a(14) = -0.05897433015592386914575323926766330d0
         ref_a(15) = -0.02168476171861335324934388684707580d0
         ref_a(16) = 0.07282080033590128173761892641234244d0
         ref_a(17) = 0.55121429634197067334405601381594315d0
         ref_a(18) = ref_a(16)
         ref_a(19) = ref_a(15)
         ref_a(20) = ref_a(14)
         ref_a(21) = ref_a(13)
         ref_a(22) = ref_a(12)
         ref_a(23) = ref_a(11)
         ref_a(24) = ref_a(10)
         ref_a(25) = ref_a(9)
         ref_a(26) = ref_a(8)
         ref_a(27) = ref_a(7)
         ref_a(28) = ref_a(6)
         ref_a(29) = ref_a(5)
         ref_a(30) = ref_a(4)
         ref_a(31) = ref_a(3)
         ref_a(32) = ref_a(2)
         ref_a(33) = ref_a(1)
         ref_a(34) = ref_a(0)
      case default
         print '(a)', 'ref_load: unknown order (reference transcription: 4/6/8 only)'
         stop 2
      end select
   end subroutine ref_load

   ! ref_step(sx, dt) - the reference composition step: the historical main-loop form
   ! (the force of the PREVIOUS evaluation is carried in sx%f; N fresh evaluations per
   ! step, at the drift ends). The expression grouping mirrors the member's per-beat
   ! arithmetic exactly - same operations in the same order on the same values - so
   ! equal coefficient tables yield bit-identical trajectories ((s-a4))
   subroutine ref_step(sx, dt)
      type(state_t), intent(inout) :: sx
      real(8), intent(in) :: dt
      integer :: i, k
      real(8) :: ed
      ed = 0.0d0
      do i = 0, ref_n - 1
         do k = 1, size(roster%mass)
            sx%p(3*k - 2:3*k) = sx%p(3*k - 2:3*k) + sx%f(3*k - 2:3*k)*(dt*ref_a(2*i))
         end do
         do k = 1, size(roster%mass)
            sx%q(3*k - 2:3*k) = sx%q(3*k - 2:3*k) + sx%p(3*k - 2:3*k)*(dt*ref_a(2*i + 1)/roster%mass(k))
         end do
         call force_eval(sx, ed)
      end do
      do k = 1, size(roster%mass)
         sx%p(3*k - 2:3*k) = sx%p(3*k - 2:3*k) + sx%f(3*k - 2:3*k)*(dt*ref_a(2*ref_n))
      end do
      sx%t = sx%t + dt
   end subroutine ref_step

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
