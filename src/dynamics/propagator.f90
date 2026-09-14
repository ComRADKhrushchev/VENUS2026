!=====================================================================
! propagator.f90 - L3 closed family: nuclear propagators (verlet/symple/radau three-member registration +
!                   unified dispatch stepping - the integrator body holds no system knowledge;
!                   a pure integrator: frame output and termination are the driver's affair)
! Family algorithms (semantics): verlet = velocity Verlet (order=1) and the
!   third-order Beeman predictor-corrector (order=2); symple = the 4/6/8th-order
!   alternating drift-kick composition (Schlier coefficients); radau = the RA15
!   Gauss-Radau adaptive sequence. Every force evaluation goes through force_eval
!   exclusively - stochastic/dissipative laws are facade members, never inlined
!   integrator code; no debug prints live in the stepping path.
!
! -- E14 hooks withdrawal (legislated; spec §8) --
!   hooks_t / hook_out_i / hook_term_i / the second prop_step argument are deleted: the propagator
!   reverts to a pure integrator. The integration loop directly calls prop_step and the facade
!   phase entries (elec_prop / mqc_statistic) plus rec_frame / container_term - the four-phase
!   orchestration is explicit in the driver. Frame output and termination move to the driver;
!   this also repairs the latent timing defect of the hooks form - a hooks%term checkpoint fired
!   after nuclear propagation but before the electronic phases, so a termination there skipped
!   the step's electronic evolution and left an inconsistent state.
!
! -- endpoint-force contract (spec §3 - the propagator member obligation, legislated here) --
!   prop_step returns only after the endpoint force has been evaluated: the ws held by elec_facade
!   corresponds to the returned sta%q. Under MQC the velocity-Verlet right-end force cannot be
!   folded into the next step's left end - the amplitude propagation of phase 2 (elec_prop) needs
!   the endpoint H - and between the two force sites the electronic phases run (a hop may change
!   the held occupation set, E16), so each step evaluates its own left-end force afresh.
!
! -- dt channel note (E9, re-anchored by E14; spec §3) --
!   dt is no longer a force_eval argument - force_eval(sta, e_tot) is the pure force entry (two
!   dummies; the E9 g_ij(dt) channel belongs to elec_prop's dt dummy, the electronic phase). The
!   member's step is the fixed dt from config (prop_cfg%dt): verlet/symple take it directly; radau
!   seeds its adaptive sequence from it (substep sites are interior - moot now that the force entry
!   carries no dt). The multi-surface x radau conflict stays caught by the method-package member
!   itself (E13, carried: a fixed-step-only law loud-stops on dt semantics it cannot live with -
!   adaptive sequence boundaries have no well-defined "one step" dt; no facade table; same source
!   as the force_facade multi-surface note), so the multi-surface path only ever sees fixed steps;
!   this note pins the adiabatic-path behavior to stay unchanged.
!
! -- Beeman force history (module-private staging - the force-differencing state of the order=2
!    branch, legislated by the family design) --
!   The Beeman corrector differences the forces of two consecutive step starts. The staged
!   history lives as module-private statics, carried across prop_step calls: a
!   continuous-propagation assumption, matching the driver's trajectory loop (one propagator
!   instance advances one trajectory; the staging is rebuilt on every boot). FIRST-STEP BOOT,
!   registered semantics: with no usable history - nothing staged yet, a size change of the
!   state, or a fresh trajectory (sta%t = 0) - the step runs the full velocity-Verlet sequence
!   and stages its left-end force as the next step's history; from the second step on, genuine
!   Beeman differencing. A velocity-Verlet step (order=1) CLEARS the staging: mode switches
!   never inherit a history they did not build.
!
!   CORRECTOR FORM (registered design decision): the velocity corrector is the two-step
!   Adams-Moulton quadrature p += dt(5F(t+dt) + 8F(t) - F(t-dt))/12; positions ride the
!   third-order predictor untouched. The classic Beeman corrector p += dt(2F(t+dt) +
!   F(t))/6 - algebraically the same, once the predictor is substituted, as the
!   position-difference form p = m(q-q_prev)/dt + dt(2F(t+dt)+F(t))/6 - is
!   exactly SECOND order: the substitution telescopes the position sequence to the Stoermer
!   recursion q_{n+1} = 2q_n - q_{n-1} + dt^2 F_n (measured scaling ratio ~4, not 8). The
!   family design legislates third-order positions at order=2, so the corrector carries the
!   third-order-consistent coefficients; staging, boot semantics, and the two-evaluations-
!   per-step cost are identical to the classic form.
!
! -- symple composition cache (module-private staging - the per-order coefficient table of
!    the drift-kick member, legislated by the family design) --
!   Coefficients: Christoph Schlier's 1999 composition sets (order 4: N = 5 stages;
!   order 6: N = 8; order 8: N = 17). The half table a(0..N) is literal and the second
!   half is mirror-copied, a(k) = a(2N-k) - the palindromic structure IS the exact
!   single-step time reversibility of the composition. The symplecticity normalization
!   reads directly off the array: the kick coefficients (even sites, a(0)+a(2)+...+
!   a(2N)) sum to 1 and the drift coefficients (odd sites) sum to 1 - together a unit
!   advance of p and q. The cache keys on order: a switch of config%propagator%order
!   reloads the table (the table is a pure function of order - no trajectory state is
!   held, contrast the Beeman staging). HYGIENE: a symple step builds no force history,
!   so it clears the Beeman staging on entry - mode switches never inherit a history
!   they did not build (the order=1 rule, applied symmetrically). A step makes N+1
!   force evaluations (one per loop head i = 0..N-1 plus the final beat); the last
!   evaluation sits at the endpoint q - the header endpoint-force contract holds with
!   no extra evaluation.
!
!   ORDER-4 SIGN CORRECTION (documented deviation from the historical table): the
!   historical order-4 half table carries the four middle kicks with flipped signs -
!   its a(2) = +1/3 with the mirror a(4) = a(6) = -1/3, a(8) = +1/3. That set is
!   exactly SECOND order: the modified-frequency series of the composed harmonic map
!   has theta = h(1 + h^2/6 + ...) - the order-4 condition (vanishing h^2 term) fails
!   in exact rational arithmetic, and the measured one-step error scales dt^3 on
!   quartic and Kepler potentials (a validated harness reproduces order 4 for
!   Forest-Ruth and order 2 for velocity Verlet, so the instrument is sound). With
!   a(2) = -1/3 the h^2 term vanishes IDENTICALLY (rational arithmetic) and the
!   measured one-step error scales dt^5 - genuine fourth order; the orders 6/8
!   tables are verified genuine as transcribed (their series vanish through h^6/h^8
!   and the measured scalings are dt^7/dt^9). The correction flips one sign of the
!   historical table; every other value, the palindromic structure and both
!   normalizations are untouched (the kick sum is sign-invariant: +-(1/3) cancels
!   either way). Criteria (s-b)/(s-d) of check_propagator anchor the fourth-order
!   behavior; (s-a4) anchors the loaded values bit-level against an independently
!   transcribed copy of the CORRECTED table.
!=====================================================================
module propagator
   use state,         only: state_t
   use config,        only: prop_cfg => propagator   ! alias to avoid clashing with the module name
                                                      ! (prop_cfg is the config instance 'propagator')
   use config_roster, only: roster
   use force_facade,  only: force_eval               ! the pure force entry (sta, e_tot) - no dt
                                                      ! dummy, see the header note
   implicit none
   private
   public :: prop_reg, prop_step
   public :: verlet_step, symple_step, radau_step

   ! Integrator member interface: advance one fixed step (verlet/symple) or one adaptive sequence (radau)
   abstract interface
      subroutine integ_proc_i(sta)
         import :: state_t
         type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      end subroutine integ_proc_i
   end interface

   ! Closed-family registry (each member registers one row at assembly time; capacity has slack)
   integer, parameter :: max_integ = 8        ! registry capacity [row]
   type :: integ_entry_t
      integer :: id = 0                       ! integrator code (the key of config%propagator%id)
      procedure(integ_proc_i), pointer, nopass :: proc => null() ! member stepping routine
   end type integ_entry_t
   type(integ_entry_t) :: reg_tbl(max_integ)  ! registry (module-private)
   integer :: n_reg = 0                       ! registered row count [row]

   ! Beeman force history (see the header note): the staged F(t) of the previous order=2
   ! step + its readiness flag
   real(8), allocatable :: f_prev(:)          ! F(t-dt) of the upcoming Beeman step [internal
                                              ! force units] (sta%f units)
   logical :: f_prev_ok = .false.             ! history usable [flag]

   ! symple composition cache (see the header note): the coefficient table a(0:2N) +
   ! the stage count N + the order key the table was loaded for (reload on change)
   real(8) :: sym_a(0:34) = 0.0d0             ! composition coefficients a(k) [dimensionless]
                                              ! (k = 0..2N; 2N = 10/16/34 at orders 4/6/8)
   integer :: sym_n = 0                       ! stage count N [stage]
   integer :: sym_order_cached = 0            ! order key of the loaded table [order]
contains
   !------------------------------------------------------------------
   ! prop_reg(id, proc) - register an integrator member into the closed-family registry
   ! Reference: docs/final_completion_plan.md §7 X4 (unified integrator protocol + three-member registration)
   subroutine prop_reg(id, proc)
      integer, intent(in) :: id           ! integrator code [-] (1=verlet/2=symple/3=radau -
                                          ! the key of config%propagator%id)
      procedure(integ_proc_i) :: proc     ! member stepping routine
      integer :: j
      do j = 1, n_reg
         if (reg_tbl(j)%id == id) then
            call stop_propi('prop_reg', 'duplicate integrator id', id, &
                            '(each closed-family member registers exactly once, at assembly)')
         end if
      end do
      if (n_reg >= max_integ) then
         call stop_prop('prop_reg', 'integrator registry full ('// &
                        'max_integ rows exceeded - the closed family has three members)')
      end if
      n_reg = n_reg + 1
      reg_tbl(n_reg)%id = id
      reg_tbl(n_reg)%proc => proc
   end subroutine prop_reg

   !------------------------------------------------------------------
   ! prop_step(sta) - step dispatch: call the registered member per config%propagator%id
   ! Reference: spec §8 (the driver four-phase loop calls this as phase 1; the orchestration
   !       lives in the driver - this module keeps only the stepping)
   subroutine prop_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (advanced one step/sequence by the member)
      integer :: j
      do j = 1, n_reg
         if (reg_tbl(j)%id == prop_cfg%id) then
            call reg_tbl(j)%proc(sta)
            return
         end if
      end do
      call stop_propi('prop_step', 'unregistered integrator id', prop_cfg%id, &
                      '(config%propagator%id names no registered closed-family row)')
   end subroutine prop_step

   !------------------------------------------------------------------
   ! verlet_step(sta) - velocity Verlet (order=1) / third-order Beeman (order=2) single step
   ! Both branches satisfy the endpoint-force contract literally: the LAST force
   ! evaluation of the step sits at the returned sta%q (velocity Verlet: the right-end
   ! half-kick force; Beeman: the post-predictor re-evaluation - the corrector that
   ! follows moves only p, never q). sta%t advances by dt; q/p/mass per the unit
   ! contract (masses from the roster, expanded per atom over its three components).
   subroutine verlet_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      integer :: k
      real(8) :: dt, e
      real(8) :: f_tm1(size(sta%q))         ! F(t-dt) of this Beeman step (saved across restaging)
      e = 0.0d0                             ! force_eval fills it; the integrator does not consume it
      dt = prop_cfg%dt

      select case (prop_cfg%order)
      case (1)                              ! velocity Verlet
         call beeman_forget()               ! mode-switch hygiene: VV steps carry no history
         call vv_advance(sta, dt)
         sta%t = sta%t + dt

      case (2)                              ! third-order Beeman
         if (beeman_boot_needed(sta)) then
            ! FIRST-STEP BOOT (registered semantics, header note): no usable history - run
            ! the velocity-Verlet sequence, staging its left-end force as the next step's
            ! history (the extra head evaluation is deterministic and unused by the update)
            call beeman_forget()
            call force_eval(sta, e)         ! F(t): left end, evaluated afresh (header contract)
            allocate (f_prev(size(sta%f)))
            f_prev = sta%f
            f_prev_ok = .true.
            call vv_advance(sta, dt)
            sta%t = sta%t + dt
         else
            f_tm1 = f_prev                  ! F(t-dt): the staged history, saved before restaging
            call force_eval(sta, e)         ! F(t): left end, evaluated afresh (header contract)
            f_prev = sta%f                  ! stage F(t) as the next step's history
            do k = 1, size(roster%mass)
               ! position predictor (third-order): q += (p/m)dt + dt^2(4F(t)-F(t-dt))/(6m)
               sta%q(3*k-2:3*k) = sta%q(3*k-2:3*k) + sta%p(3*k-2:3*k)*(dt/roster%mass(k)) &
                  + (4.0d0*sta%f(3*k-2:3*k) - f_tm1(3*k-2:3*k))*((dt*dt)/(6.0d0*roster%mass(k)))
            end do
            call force_eval(sta, e)         ! F(t+dt) at the predicted (final) q - the endpoint
            do k = 1, size(roster%mass)
               ! momentum corrector (two-step Adams-Moulton, header note): p += dt(5F(t+dt)
               ! + 8F(t) - F(t-dt))/12 (q untouched - positions stay on the predictor)
               sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) &
                  + (5.0d0*sta%f(3*k-2:3*k) + 8.0d0*f_prev(3*k-2:3*k) - f_tm1(3*k-2:3*k)) &
                  *(dt/12.0d0)
            end do
            sta%t = sta%t + dt
         end if

      case default
         call stop_propi('verlet_step', 'illegal propagator order for verlet', prop_cfg%order, &
                         '(1 = velocity Verlet, 2 = Beeman)')
      end select
   end subroutine verlet_step

   !------------------------------------------------------------------
   ! symple_step(sta) - 4th/6th/8th-order symplectic single step (alternating drift-kick
   !                   composition, the Schlier coefficient sets of the header note)
   ! The legislated sequence: per-order coefficient load on first use (the cache keys on
   ! order), then the alternating loop - force_eval BEFORE each kick (N evaluations),
   ! p += dt*a(2i)*f, q += dt*a(2i+1)*p/m - then the final beat (force_eval at the
   ! endpoint q, p += dt*a(2N)*f) and sta%t += dt. Endpoint-force contract: the final
   ! beat's evaluation is the LAST of the step and sits at the returned sta%q. sta%f on
   ! return is the composed endpoint force; q/p/mass per the unit contract (as verlet_step).
   subroutine symple_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      integer :: i, k
      real(8) :: dt, e
      e = 0.0d0                             ! force_eval fills it; the integrator does not consume it
      dt = prop_cfg%dt

      call beeman_forget()                  ! hygiene: this member builds no force history
                                           ! (the header-note rule, applied symmetrically)
      if (sym_order_cached /= prop_cfg%order) call symple_load(prop_cfg%order)

      do i = 0, sym_n - 1                   ! alternating loop: kick a(2i), drift a(2i+1)
         call force_eval(sta, e)            ! F at the current q, before the kick
         do k = 1, size(roster%mass)
            sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) + sta%f(3*k-2:3*k)*(dt*sym_a(2*i))
         end do
         do k = 1, size(roster%mass)
            sta%q(3*k-2:3*k) = sta%q(3*k-2:3*k) + sta%p(3*k-2:3*k)*(dt*sym_a(2*i + 1)/roster%mass(k))
         end do
      end do
      call force_eval(sta, e)               ! the final beat's force site = the endpoint q
      do k = 1, size(roster%mass)
         sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) + sta%f(3*k-2:3*k)*(dt*sym_a(2*sym_n))
      end do
      sta%t = sta%t + dt
   end subroutine symple_step

   !------------------------------------------------------------------
   ! radau_step(sta) - one RA15 Gauss-Radau sequence (adaptive step)
   ! Not yet realized - the adaptive member arrives with its own concretization task;
   ! the legislated sequence is kept as the roadmap:
   ! 1. On first call load the constants: Radau nodes hh(2:8); w/u weights; triangle
   !    coefficients c/d/ra
   ! 2. Fold velocities v = p/m; iteration count 6 for the first sequence, 2 afterwards
   ! 3. Substep loop within the sequence (iterations x 7 substeps): predict substep
   !    coordinates from the seventh-order B-coefficient expansion; force_eval at the
   !    substep coordinates (f_j = f/m fold); difference correction - divided differences
   !    of the g coefficients + triangle accumulation of the B coefficients
   ! 4. Sequence close-out: advance q/p by the sequence length tt; last-term hv estimate;
   !    restart test - if tp/tt <= 1 discard the sequence and restart
   ! 5. Adaptive step: tp = tcorr*(ss/hv)^(1/9) (ss = 10^(-order), direction follows tf);
   !    step clamping dt_min/dt_max through the deferred config fields
   ! 6. sta%t += sequence length tt
   ! 7. Endpoint-force close-out: force_eval(sta, e) at the returned sta%q - the substep
   !    force sites sit at interior node coordinates, not the sequence endpoint, so the
   !    member closes with an explicit endpoint evaluation (the header contract)
   subroutine radau_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
   end subroutine radau_step

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! vv_advance(sta, dt) - the velocity-Verlet update, shared by the order=1 branch and
   ! the Beeman boot (the boot step is thereby IDENTICAL to an order=1 step from the
   ! same state): left-end force afresh, half kick, full drift, endpoint force, half kick
   subroutine vv_advance(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
      integer :: k
      real(8) :: e
      e = 0.0d0
      call force_eval(sta, e)               ! F(t): left end, evaluated afresh
      do k = 1, size(roster%mass)
         sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) + sta%f(3*k-2:3*k)*(0.5d0*dt)
      end do
      do k = 1, size(roster%mass)
         sta%q(3*k-2:3*k) = sta%q(3*k-2:3*k) + sta%p(3*k-2:3*k)*(dt/roster%mass(k))
      end do
      call force_eval(sta, e)               ! F(t+dt): the endpoint (header contract)
      do k = 1, size(roster%mass)
         sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) + sta%f(3*k-2:3*k)*(0.5d0*dt)
      end do
   end subroutine vv_advance

   ! beeman_forget() - clear the staged Beeman history (boot hygiene / VV mode switch)
   subroutine beeman_forget()
      if (allocated(f_prev)) deallocate (f_prev)
      f_prev_ok = .false.
   end subroutine beeman_forget

   ! beeman_boot_needed(sta) - the registered boot condition: no usable history
   ! (nothing staged yet, or a size change of the state) or a fresh trajectory
   ! (sta%t = 0). Split into its own function so size(f_prev) is evaluated ONLY
   ! when the staging flag says the array is allocated - Fortran does not
   ! guarantee short-circuit .or., and a size() on an unallocated array is
   ! nonconforming (the flag and the allocation move together: forget clears
   ! both, the boot sets both)
   logical function beeman_boot_needed(sta)
      type(state_t), intent(in) :: sta
      beeman_boot_needed = .true.
      if (f_prev_ok) then
         beeman_boot_needed = (size(f_prev) /= size(sta%q) .or. sta%t == 0.0d0)
      end if
   end function beeman_boot_needed

   ! symple_load(order) - fill the composition cache for one order (4/6/8; loud-stop
   ! otherwise, named for the member). The half table a(0..N) is literal, the second
   ! half mirror-copied (a(2N-k) = a(k), k = 0..N-1) - Schlier's 1999 sets, the values
   ! verified bit-level against an independently transcribed historical copy (the
   ! header note); the table is a pure function of order, so the cache reloads on a
   ! change of the order key
   subroutine symple_load(order)
      integer, intent(in) :: order          ! composition order [-] (4, 6 or 8)
      integer :: k
      select case (order)
      case (4)                              ! N = 5 (SIGN-CORRECTED vs the historical
                                           ! table - the correction is legislated in the
                                           ! module header note; criterion (s-b) anchors it)
         sym_n = 5
         sym_a(0) = 0.5d0
         sym_a(1) = -1.0d0/48.0d0
         sym_a(2) = -1.0d0/3.0d0
         sym_a(3) = 3.0d0/8.0d0
         sym_a(4) = 1.0d0/3.0d0
         sym_a(5) = 7.0d0/24.0d0
      case (6)                              ! N = 8
         sym_n = 8
         sym_a(0) = 0.06942944346252987735848865824703402d0
         sym_a(1) = 0.2848783771728008405274534645665783d0
         sym_a(2) = -0.1331551983159820940996130995137351d0
         sym_a(3) = 0.3278397575961294541205467836732554d0
         sym_a(4) = 0.0012903891798107897423048174644328d0
         sym_a(5) = -0.3812210427193262947562278437421127d0
         sym_a(6) = 0.4224353656736414269988196238022683d0
         sym_a(7) = 0.268502907950396000108227595502279d0
         sym_a(8) = 0.28d0
      case (8)                              ! N = 17
         sym_n = 17
         sym_a(0) = 0.04463795052359022755913999625733590d0
         sym_a(1) = 0.13593258071690959145543264213495574d0
         sym_a(2) = 0.2198844042714707225445535069606167d0
         sym_a(3) = 0.13024946780523828601621193778196846d0
         sym_a(4) = 0.10250365693975069608261241007779814d0
         sym_a(5) = 0.43234521869358547487983257884877035d0
         sym_a(6) = -0.00477482916916881658022489063962934d0
         sym_a(7) = -0.58253476904040845493112837930861213d0
         sym_a(8) = -0.03886264282111817697737420875189743d0
         sym_a(9) = 0.31548728537940479698273603797274199d0
         sym_a(10) = 0.18681583743297155471526153503972746d0
         sym_a(11) = 0.26500275499062083398346002963079872d0
         sym_a(12) = -0.02405084735747361993573587982407554d0
         sym_a(13) = -0.45040492499772251180922896712151891d0
         sym_a(14) = -0.05897433015592386914575323926766330d0
         sym_a(15) = -0.02168476171861335324934388684707580d0
         sym_a(16) = 0.07282080033590128173761892641234244d0
         sym_a(17) = 0.55121429634197067334405601381594315d0
      case default
         call stop_propi('symple_step', 'illegal propagator order for symple', order, &
                         '(4, 6 or 8 - the composition orders)')
      end select
      do k = 0, sym_n - 1                   ! palindromic mirror of the second half
         sym_a(2*sym_n - k) = sym_a(k)
      end do
      sym_order_cached = order
   end subroutine symple_load

   ! loud-stop channel (name the caller; value variant; STOP 1)
   subroutine stop_prop(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': propagator error: '//trim(msg)
      write (0, '(a)') trim(who)//': loud-stop (the state is not advanced)'
      stop 1
   end subroutine stop_prop

   subroutine stop_propi(who, msg, ival, tail)
      character(len=*), intent(in) :: who, msg, tail
      integer, intent(in) :: ival
      write (0, '(a,a,1x,i0,a,a)') trim(who)//': ', trim(msg), ival, ' ', trim(tail)
      write (0, '(a)') trim(who)//': loud-stop (the state is not advanced)'
      stop 1
   end subroutine stop_propi

end module propagator
