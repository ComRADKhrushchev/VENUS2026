!=====================================================================
! propagator.f90 - closed family of nuclear propagators (verlet/symple/radau + prop_step dispatch)
! Design:
!   Pure integrators: all forces come through force_eval(sta, e_tot) - the eV ->
!   internal-unit conversion lives in force_interface; frame output and termination
!   are the driver's affair. Contract: prop_step returns only after the endpoint
!   force at the returned sta%q has been evaluated. Members hold module-private
!   buffer (Beeman force history, symple coefficient cache, radau cross-sequence
!   carry), rebuilt on a fresh trajectory / size change / step-size change and
!   cleared on member switch. Units: time [10 fs], length [A]; masses from the
!   list_atoms, expanded per atom over its three components.
!=====================================================================
module propagator
   use state,         only: state_t
   use config,        only: prop_cfg => propagator   ! alias: the config instance 'propagator'
                                                      ! clashes with the module name
   use config_atoms, only: list_atoms
   use force_interface,  only: force_eval               ! the pure force entry (sta, e_tot)
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

   ! Beeman force history: the buffered F(t) of the previous order=2 step + readiness flag
   real(8), allocatable :: f_prev(:)          ! buffered F(t-dt) of the upcoming Beeman step
                                              ! (sta%f units)
   logical :: f_prev_ok = .false.             ! history usable [flag]
   real(8) :: dt_staged = 0.0d0               ! step size the buffered f_prev is spaced at [10 fs]
                                              ! (a dt change re-boots the equal-spacing AM corrector)

   ! symple composition cache: coefficient table a(0:2N), stage count N, order key
   real(8) :: sym_a(0:34) = 0.0d0             ! composition coefficients a(k), k = 0..2N
                                              ! (2N = 10/16/34 at orders 4/6/8)
   integer :: sym_n = 0                       ! stage count N [stage]
   integer :: sym_order_cached = 0            ! order key of the loaded table [order]

   ! radau constants: loaded once - nodes, w/u weights, 21-entry triangle tables
   integer, parameter :: rad_nw(8) = [ 0, 0, 1, 3, 6, 10, 15, 21 ]  ! column-start map:
                                              ! entries nw(k)+1..nw(k+1) form triangle column k
   real(8) :: rad_hh(8) = 0.0d0               ! Gauss-Radau nodes hh(1..8), hh(1) = 0 [s]
   real(8) :: rad_w(7) = 0.0d0                ! position close-out weights 1/(n(n+1)), n = j+1
   real(8) :: rad_u(7) = 0.0d0                ! velocity close-out weights 1/n, n = j+1
   real(8) :: rad_cc(21) = 0.0d0              ! triangle: B-corrector accumulation coefficients
   real(8) :: rad_d(21) = 0.0d0               ! triangle: monomial-B -> differenced-g conversion
   real(8) :: rad_ra(21) = 0.0d0              ! triangle: difference recursion = 1/(node diffs)
   logical :: rad_const_ok = .false.          ! constants loaded [flag]

   ! radau cross-sequence carry: the B-coefficient buffer of the upcoming sequence
   ! (trajectory state, carried across prop_step calls, rebuilt on every boot)
   real(8), allocatable :: rad_bb(:, :)       ! monomial B coefficients (7, n_dof)
   real(8), allocatable :: rad_ee(:, :)       ! binomial extrapolation onto the next scale
   real(8), allocatable :: rad_bd(:, :)       ! running correction B_converged - ee_previous
   real(8), allocatable :: rad_f1(:)          ! left-end folded acceleration of the sequence
   real(8), allocatable :: rad_v(:)           ! folded velocity (the working representation)
   real(8) :: rad_tp = 0.0d0                  ! upcoming sequence length [10 fs]
   logical :: rad_seq_ok = .false.            ! a sequence has been accepted (restart test off)
   integer :: rad_ns = 0                      ! accepted-sequence counter (bd rule: 1st keeps 0)
contains
   !------------------------------------------------------------------
   ! prop_reg(id, proc) - register an integrator member into the closed-family registry
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
   ! Both branches end with the endpoint force at the returned sta%q; sta%t advances by dt.
   subroutine verlet_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      integer :: k
      real(8) :: dt, e
      real(8) :: f_tm1(size(sta%q))         ! F(t-dt) of this Beeman step (saved across the swap)
      e = 0.0d0                             ! force_eval fills it; the integrator does not consume it
      dt = prop_cfg%dt

      call rad_forget()                     ! hygiene: this member builds no radau carry

      select case (prop_cfg%order)
      case (1)                              ! 1. velocity Verlet
         call beeman_forget()               ! mode-switch hygiene: VV steps carry no history
         call vv_advance(sta, dt)
         sta%t = sta%t + dt

      case (2)                              ! 2. third-order Beeman
         if (beeman_boot_needed(sta)) then
            ! FIRST-STEP BOOT: no usable history - run the velocity-Verlet sequence,
            ! buffer its left-end force as the next step's history
            call beeman_forget()
            call force_eval(sta, e)         ! F(t): left end, evaluated afresh (endpoint contract)
            allocate (f_prev(size(sta%f)))
            f_prev = sta%f
            f_prev_ok = .true.
            dt_staged = dt                  ! the history is spaced at THIS step's dt
            call vv_advance(sta, dt)
            sta%t = sta%t + dt
         else
            f_tm1 = f_prev                  ! F(t-dt): the buffered history, saved before the swap
            call force_eval(sta, e)         ! F(t): left end, evaluated afresh (header contract)
            f_prev = sta%f                  ! stage F(t) as the next step's history
            do k = 1, size(list_atoms%mass)
               ! position predictor (third-order): q += (p/m)dt + dt^2(4F(t)-F(t-dt))/(6m)
               sta%q(3*k-2:3*k) = sta%q(3*k-2:3*k) + sta%p(3*k-2:3*k)*(dt/list_atoms%mass(k)) &
                  + (4.0d0*sta%f(3*k-2:3*k) - f_tm1(3*k-2:3*k))*((dt*dt)/(6.0d0*list_atoms%mass(k)))
            end do
            call force_eval(sta, e)         ! F(t+dt) at the predicted (final) q - the endpoint
            do k = 1, size(list_atoms%mass)
               ! momentum corrector (two-step Adams-Moulton): p += dt(5F(t+dt)+8F(t)-F(t-dt))/12
               ! (q untouched - positions stay on the predictor)
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
   !                   composition, Schlier coefficient sets)
   ! Load the per-order coefficient table on first use, run the alternating loop
   ! (force_eval before each kick, p += dt*a(2i)*f, q += dt*a(2i+1)*p/m), then the
   ! final beat at the endpoint q and sta%t += dt - the final evaluation is the last
   ! of the step and sits at the returned sta%q.
   subroutine symple_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      integer :: i, k
      real(8) :: dt, e
      e = 0.0d0                             ! force_eval fills it; the integrator does not consume it
      dt = prop_cfg%dt

      call beeman_forget()                  ! hygiene: this member builds no force history
      call rad_forget()                     ! and no radau carry
      if (sym_order_cached /= prop_cfg%order) call symple_load(prop_cfg%order)

      do i = 0, sym_n - 1                   ! 1. alternating loop: kick a(2i), drift a(2i+1)
         call force_eval(sta, e)            ! F at the current q, before the kick
         do k = 1, size(list_atoms%mass)
            sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) + sta%f(3*k-2:3*k)*(dt*sym_a(2*i))
         end do
         do k = 1, size(list_atoms%mass)
            sta%q(3*k-2:3*k) = sta%q(3*k-2:3*k) + sta%p(3*k-2:3*k)*(dt*sym_a(2*i + 1)/list_atoms%mass(k))
         end do
      end do
      call force_eval(sta, e)               ! 2. the final beat's force site = the endpoint q
      do k = 1, size(list_atoms%mass)
         sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) + sta%f(3*k-2:3*k)*(dt*sym_a(2*sym_n))
      end do
      sta%t = sta%t + dt
   end subroutine symple_step

   !------------------------------------------------------------------
   ! radau_step(sta) - one RA15 Gauss-Radau sequence (adaptive step)
   ! Boot as needed, then ONE accepted sequence per call: q/p/t advance by the
   ! sequence length tt, the close-out evaluates the endpoint force explicitly
   ! (sta%f on return is the force at the returned sta%q), and the B-carry/next
   ! length are buffered for the following call. p is written back from the folded
   ! working velocity at close-out.
   subroutine radau_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      integer :: m, j, k, l, i, jd, ni, nq, nrst
      real(8) :: e, s, gk, acc, temp, tt, t2, tval, hv, dir, ss, tpn, qr
      real(8) :: gg(7, size(sta%q))         ! differenced g coefficients of the sequence
      real(8) :: yy(size(sta%q))            ! predicted substep geometry
      real(8) :: fj(size(sta%q))            ! folded substep force (f/m)
      real(8) :: qs(size(sta%q))            ! sequence-origin geometry save/restore
      e = 0.0d0                             ! force_eval fills it; the integrator does not consume it

      ! 1. order = the accuracy exponent (ss = 10^-order); dt gives the direction and the
      !    first-sequence seed - both illegal values are named aborts
      if (prop_cfg%order < 1 .or. prop_cfg%order > 15) then
         call stop_propi('radau_step', 'illegal propagator order for radau', prop_cfg%order, &
                         '(1..15 - the accuracy exponent of ss = 10^-order)')
      end if
      if (prop_cfg%dt == 0.0d0) then
         call stop_prop('radau_step', 'dt = 0 (dt seeds the first sequence and gives '// &
                        'the integration direction - a zero seed has neither)')
      end if
      dir = 1.0d0
      if (prop_cfg%dt < 0.0d0) dir = -1.0d0
      ss = 10.0d0**(-prop_cfg%order)

      ! 2. hygiene + boot as needed; re-fold the working velocity from sta%p so an
      !    external write between sequences survives the close-out write-back
      call beeman_forget()                  ! hygiene: this member builds no Beeman history
      if (.not. rad_const_ok) call radau_load()
      if (radau_boot_needed(sta)) call radau_boot(sta)
      nq = size(sta%q)
      do k = 1, nq                          ! re-fold the working velocity from sta%p
         rad_v(k) = sta%p(k)/list_atoms%mass((k + 2)/3)
      end do
      ni = 6                                ! first sequence iterates six times,
      if (rad_seq_ok) ni = 2                ! later sequences twice (the B-carry)

      ! ---- 3. sequence trial loop (restart-capable; the accept test fires on the
      ! first sequence only - once a sequence has been accepted the carried B makes
      ! rejection pointless) ----
      nrst = 0
      do
         ! g coefficients from the carried B (triangle-d conversion)
         do k = 1, nq
            do i = 1, 7
               acc = rad_bb(i, k)
               do j = i + 1, 7
                  acc = acc + rad_d(rad_nw(j) + i)*rad_bb(j, k)
               end do
               gg(i, k) = acc
            end do
         end do
         tt = rad_tp
         t2 = tt*tt
         tval = real(abs(tt), 4)            ! the historical single-precision floor of hv
         hv = 0.0d0
         do m = 1, ni
            do j = 2, 8
               jd = j - 1
               s = rad_hh(j)
               do k = 1, nq
                  ! predictor: the seventh-order B expansion at the node s (Horner form)
                  acc = rad_w(3)*rad_bb(3, k) + s*(rad_w(4)*rad_bb(4, k) &
                        + s*(rad_w(5)*rad_bb(5, k) + s*(rad_w(6)*rad_bb(6, k) + s*rad_w(7)*rad_bb(7, k))))
                  yy(k) = sta%q(k) + s*(tt*rad_v(k) + t2*s*(rad_f1(k)*0.5d0 &
                          + s*(rad_w(1)*rad_bb(1, k) + s*(rad_w(2)*rad_bb(2, k) + s*acc))))
               end do
               ! force at the predicted substep geometry (the geometry swap rides sta%q;
               ! sta%f is the force workspace, restored before the fold consumes it)
               qs = sta%q
               sta%q = yy
               call force_eval(sta, e)
               do k = 1, nq
                  sta%q(k) = qs(k)
                  fj(k) = sta%f(k)/list_atoms%mass((k + 2)/3)
               end do
               do k = 1, nq
                  ! corrector: difference the folded substep force against the head value,
                  ! walk the ra recursion, accumulate through the cc triangle into B
                  temp = gg(jd, k)
                  gk = (fj(k) - rad_f1(k))/s
                  if (jd == 1) then
                     gg(1, k) = gk
                  else
                     acc = gk - gg(1, k)
                     do i = 2, jd - 1
                        acc = acc*rad_ra(rad_nw(jd) + i - 1) - gg(i, k)
                     end do
                     gg(jd, k) = acc*rad_ra(rad_nw(jd) + jd - 1)
                  end if
                  temp = gg(jd, k) - temp
                  rad_bb(jd, k) = rad_bb(jd, k) + temp
                  do l = 1, jd - 1
                     rad_bb(l, k) = rad_bb(l, k) + rad_cc(rad_nw(jd) + l)*temp
                  end do
               end do
            end do
            if (m == ni) then               ! last-iteration error estimate (hv)
               do k = 1, nq
                  hv = max(hv, abs(rad_bb(7, k)))
               end do
               hv = hv*rad_w(7)/tval**7
               ! diverged-trial guard: hv non-finite (NaN or Inf) or beyond the
               ! divergence bound 1e30 aborts with a named error - a diverged trial under a dt_min
               ! floor would otherwise clamp back up to the trial length and accept
               ! garbage through the equality path, unnamed (legitimate measured hv
               ! peaks near 9e8, so the bound keeps >20 decades of margin)
               if (hv /= hv .or. hv > 1.0d30) then
                  call stop_prop('radau_step', 'diverged trial sequence: the error'// &
                                 ' estimate hv is non-finite or beyond the divergence'// &
                                 ' bound 1e30 (the dt seed or the dt_min floor exceeds'// &
                                 ' what this force law supports)')
               end if
            end if
         end do

         ! ---- first-sequence accept test: adapt from hv, clamp, accept - or discard
         ! and retry at 0.8*tp (the converged B survives as the retry's predictor; a
         ! probe floored back to the trial length accepts at equality, else it would
         ! retry the same length forever) ----
         if (rad_seq_ok) exit               ! an accepted past sequence: no restart path
         tpn = (ss/hv)**(1.0d0/9.0d0)*dir
         call rad_clamp(tpn, dir)
         rad_tp = tpn
         if (rad_tp/tt >= 1.0d0) exit       ! accept (>= : equality is the floor-pinned
         rad_tp = 0.8d0*rad_tp              ! case - in the unclamped world tp == tt is
         call rad_clamp(rad_tp, dir)        ! measure-zero, so >= is the published >)
         nrst = nrst + 1                    ! termination guard: a diverging trial never
         if (nrst > 10) then                ! accepts - a named abort beats a silent spin
            call stop_prop('radau_step', & ! forever
                           'first sequence failed to accept after 10 restarts (the'// &
                           ' error estimate diverges - the dt seed is too large for'// &
                           ' this force law, or the trial overflowed)')
         end if
         ni = 6                             ! the retry runs the full iteration count
         call force_eval(sta, e)            ! head force afresh at the same origin
         do k = 1, nq                       ! (the retry's f1; the state never advanced)
            rad_f1(k) = sta%f(k)/list_atoms%mass((k + 2)/3)
         end do
      end do

      ! ---- 4. sequence close-out: advance by tt through the exact B moments (w/u),
      ! guard the closed-out state, evaluate the endpoint force and stage the next
      ! sequence's head values ----
      do k = 1, nq
         sta%q(k) = sta%q(k) + rad_v(k)*tt + t2*(rad_f1(k)*0.5d0 + rad_bb(1, k)*rad_w(1) &
                    + rad_bb(2, k)*rad_w(2) + rad_bb(3, k)*rad_w(3) + rad_bb(4, k)*rad_w(4) &
                    + rad_bb(5, k)*rad_w(5) + rad_bb(6, k)*rad_w(6) + rad_bb(7, k)*rad_w(7))
         rad_v(k) = rad_v(k) + tt*(rad_f1(k) + rad_bb(1, k)*rad_u(1) + rad_bb(2, k)*rad_u(2) &
                    + rad_bb(3, k)*rad_u(3) + rad_bb(4, k)*rad_u(4) + rad_bb(5, k)*rad_u(5) &
                    + rad_bb(6, k)*rad_u(6) + rad_bb(7, k)*rad_u(7))
      end do

      ! ---- close-out state guard (per component, no reduction): non-finite or
      ! physically-absurd values abort with a named error - the max()-reduced hv estimate cannot
      ! catch NaN components (Fortran max with a NaN argument is
      ! processor-dependent). Bounds: |q| > 1e8 A or |v| > 1e10 A/(10 fs)
      ! (legitimate runs peak |q| ~ 5, |v| ~ 50; p = m*v rides the v test) ----
      do k = 1, nq
         if (sta%q(k) /= sta%q(k) .or. rad_v(k) /= rad_v(k) &
             .or. abs(sta%q(k)) > 1.0d8 .or. abs(rad_v(k)) > 1.0d10) then
            call stop_prop('radau_step', 'diverged sequence state: a closed-out'// &
                           ' q/p component is non-finite or beyond the physical-'// &
                           'absurdity bound (q 1e8 A, v 1e10 A/(10 fs)) - the'// &
                           ' accepted sequence was garbage (the dt seed or the'// &
                           ' dt_min floor exceeds what this force law supports)')
         end if
      end do

      sta%t = sta%t + tt
      rad_seq_ok = .true.
      rad_ns = rad_ns + 1
      call force_eval(sta, e)               ! the endpoint force (endpoint contract)
      do k = 1, nq
         rad_f1(k) = sta%f(k)/list_atoms%mass((k + 2)/3)   ! next sequence's head acceleration
         sta%p(k) = rad_v(k)*list_atoms%mass((k + 2)/3)    ! p written back from the working velocity
      end do

      ! ---- 5. adaptive length of the next sequence: (ss/hv)^(1/9), growth-limited
      ! to 1.4*tt, then clamped (the limiter shapes the request, the bounds are hard) ----
      tpn = dir*(ss/hv)**(1.0d0/9.0d0)
      if (tpn/tt > 1.4d0) tpn = tt*1.4d0
      call rad_clamp(tpn, dir)
      rad_tp = tpn

      ! ---- 6. B-carry onto the next scale: binomial shift of the converged B by the
      ! length ratio, plus the running correction (bd stays zero through the first
      ! accepted sequence - no earlier extrapolation to difference against) ----
      qr = rad_tp/tt
      do k = 1, nq
         if (rad_ns /= 1) then
            do l = 1, 7
               rad_bd(l, k) = rad_bb(l, k) - rad_ee(l, k)
            end do
         end if
         rad_ee(1, k) = qr*(rad_bb(1, k) + 2.0d0*rad_bb(2, k) + 3.0d0*rad_bb(3, k) &
                          + 4.0d0*rad_bb(4, k) + 5.0d0*rad_bb(5, k) + 6.0d0*rad_bb(6, k) &
                          + 7.0d0*rad_bb(7, k))
         rad_ee(2, k) = qr*qr*(rad_bb(2, k) + 3.0d0*rad_bb(3, k) + 6.0d0*rad_bb(4, k) &
                              + 10.0d0*rad_bb(5, k) + 15.0d0*rad_bb(6, k) + 21.0d0*rad_bb(7, k))
         rad_ee(3, k) = (qr**3)*(rad_bb(3, k) + 4.0d0*rad_bb(4, k) + 10.0d0*rad_bb(5, k) &
                               + 20.0d0*rad_bb(6, k) + 35.0d0*rad_bb(7, k))
         rad_ee(4, k) = (qr**4)*(rad_bb(4, k) + 5.0d0*rad_bb(5, k) + 15.0d0*rad_bb(6, k) &
                               + 35.0d0*rad_bb(7, k))
         rad_ee(5, k) = (qr**5)*(rad_bb(5, k) + 6.0d0*rad_bb(6, k) + 21.0d0*rad_bb(7, k))
         rad_ee(6, k) = (qr**6)*(rad_bb(6, k) + 7.0d0*rad_bb(7, k))
         rad_ee(7, k) = (qr**7)*rad_bb(7, k)
         do l = 1, 7
            rad_bb(l, k) = rad_ee(l, k) + rad_bd(l, k)
         end do
      end do
   end subroutine radau_step

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! vv_advance(sta, dt) - the velocity-Verlet update, shared by the order=1 branch and
   ! the Beeman boot (the boot step is thereby identical to an order=1 step from the
   ! same state): left-end force afresh, half kick, full drift, endpoint force, half kick
   subroutine vv_advance(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
      integer :: k
      real(8) :: e
      e = 0.0d0
      call force_eval(sta, e)               ! F(t): left end, evaluated afresh
      do k = 1, size(list_atoms%mass)
         sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) + sta%f(3*k-2:3*k)*(0.5d0*dt)
      end do
      do k = 1, size(list_atoms%mass)
         sta%q(3*k-2:3*k) = sta%q(3*k-2:3*k) + sta%p(3*k-2:3*k)*(dt/list_atoms%mass(k))
      end do
      call force_eval(sta, e)               ! F(t+dt): the endpoint (header contract)
      do k = 1, size(list_atoms%mass)
         sta%p(3*k-2:3*k) = sta%p(3*k-2:3*k) + sta%f(3*k-2:3*k)*(0.5d0*dt)
      end do
   end subroutine vv_advance

   ! beeman_forget() - clear the buffered Beeman history (boot hygiene / VV mode switch)
   subroutine beeman_forget()
      if (allocated(f_prev)) deallocate (f_prev)
      f_prev_ok = .false.
      dt_staged = 0.0d0                     ! the flag, the array, and its spacing move together
   end subroutine beeman_forget

   ! beeman_boot_needed(sta) - the boot condition: no usable history (nothing buffered,
   ! a size change of the state) or a fresh trajectory (sta%t = 0). Own function so
   ! size(f_prev) is evaluated only when the buffer flag says the array is allocated
   ! (Fortran does not guarantee short-circuit .or.; a size() on an unallocated
   ! array is nonconforming)
   logical function beeman_boot_needed(sta)
      type(state_t), intent(in) :: sta
      beeman_boot_needed = .true.
      if (f_prev_ok) then
         ! the dt key: a dt swap re-boots - the equal-spacing AM quadrature would
         ! silently degrade on history buffered at another spacing
         beeman_boot_needed = (size(f_prev) /= size(sta%q) .or. sta%t == 0.0d0 &
                               .or. prop_cfg%dt /= dt_staged)
      end if
   end function beeman_boot_needed

   ! symple_load(order) - fill the composition cache for one order (4/6/8; named abort
   ! otherwise). The half table a(0..N) is literal and the second half mirror-copied
   ! (a(2N-k) = a(k)) - Schlier's 1999 sets, verified bit-level by the check program;
   ! the cache reloads on a change of the order key (a pure function of order)
   subroutine symple_load(order)
      integer, intent(in) :: order          ! composition order [-] (4, 6 or 8)
      integer :: k
      select case (order)
      case (4)                              ! N = 5 (sign-corrected vs the historical table:
                                           ! the historical a(2) is exactly second order)
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

   ! radau_load() - fill the constant cache once: the Gauss-Radau nodes (canonical
   ! published values, hh(1) = 0), the closed-form w/u weights, and the 21-entry
   ! cc/d/ra triangles from Everhart's column recurrences (column k fills entries
   ! nw(k)+1..nw(k+1) from column k-1; column 2 is seeded first; ra entries are
   ! reciprocals of node differences by construction)
   subroutine radau_load()
      integer :: k, l, la, lb, lc, ld, le, n
      rad_hh(1) = 0.0d0
      rad_hh(2) = 0.05626256053692215d0
      rad_hh(3) = 0.18024069173689236d0
      rad_hh(4) = 0.35262471711316964d0
      rad_hh(5) = 0.54715362633055538d0
      rad_hh(6) = 0.73421017721541053d0
      rad_hh(7) = 0.88532094683909577d0
      rad_hh(8) = 0.97752061356128750d0
      do n = 2, 8                           ! w(j) = 1/(n(n+1)), u(j) = 1/n (n = j+1) -
         rad_w(n - 1) = 1.0d0/(n + n*n)    ! exact moments of s^j on [0,1]
         rad_u(n - 1) = 1.0d0/n
      end do
      rad_cc(1) = -rad_hh(2)                ! column 2 seed
      rad_d(1) = rad_hh(2)
      rad_ra(1) = 1.0d0/(rad_hh(3) - rad_hh(2))
      la = 1
      lc = 1
      do k = 3, 7                           ! columns 3..7
         lb = la
         la = lc + 1
         lc = rad_nw(k + 1)
         rad_cc(la) = -rad_hh(k)*rad_cc(lb)
         rad_cc(lc) = rad_cc(la - 1) - rad_hh(k)
         rad_d(la) = rad_hh(2)*rad_d(lb)
         rad_d(lc) = -rad_cc(lc)
         rad_ra(la) = 1.0d0/(rad_hh(k + 1) - rad_hh(2))
         rad_ra(lc) = 1.0d0/(rad_hh(k + 1) - rad_hh(k))
         do l = 4, k                        ! column interior
            ld = la + l - 3
            le = lb + l - 4
            rad_cc(ld) = rad_cc(le) - rad_hh(k)*rad_cc(le + 1)
            rad_d(ld) = rad_d(le) + rad_hh(l - 1)*rad_d(le + 1)
            rad_ra(ld) = 1.0d0/(rad_hh(k + 1) - rad_hh(l - 1))
         end do
      end do
      rad_const_ok = .true.
   end subroutine radau_load

   ! radau_boot(sta) - first-call buffer: zero the B-carry, fold the working velocity
   ! from p, seed the sequence length from dt (clamped), and evaluate the head force
   ! at the current geometry. A positive dt_min/dt_max pair with dt_min > dt_max is
   ! contradictory - named here
   subroutine radau_boot(sta)
      type(state_t), intent(inout) :: sta
      integer :: k
      real(8) :: e, dir
      e = 0.0d0
      if (prop_cfg%dt_min > 0.0d0 .and. prop_cfg%dt_max > 0.0d0 .and. &
          prop_cfg%dt_min > prop_cfg%dt_max) then
         call stop_prop('radau_step', 'dt_min exceeds dt_max (both clamps positive - '// &
                        'no sequence length can satisfy the pair)')
      end if
      if (allocated(rad_bb)) deallocate (rad_bb, rad_ee, rad_bd, rad_f1, rad_v)
      allocate (rad_bb(7, size(sta%q)), rad_ee(7, size(sta%q)), rad_bd(7, size(sta%q)), &
                rad_f1(size(sta%q)), rad_v(size(sta%q)))
      rad_bb = 0.0d0
      rad_ee = 0.0d0
      rad_bd = 0.0d0
      rad_seq_ok = .false.
      rad_ns = 0
      dir = 1.0d0
      if (prop_cfg%dt < 0.0d0) dir = -1.0d0
      rad_tp = prop_cfg%dt                  ! the seed (clamped - a bound on lengths)
      call rad_clamp(rad_tp, dir)
      do k = 1, size(sta%q)                 ! fold the working velocity
         rad_v(k) = sta%p(k)/list_atoms%mass((k + 2)/3)
      end do
      call force_eval(sta, e)               ! head force at the boot geometry
      do k = 1, size(sta%q)
         rad_f1(k) = sta%f(k)/list_atoms%mass((k + 2)/3)
      end do
   end subroutine radau_boot

   ! radau_boot_needed(sta) - the boot condition: no buffer yet, a size change of the
   ! state, or a fresh trajectory (sta%t = 0). allocated() is unconditionally safe;
   ! the size probe runs only under it
   logical function radau_boot_needed(sta)
      type(state_t), intent(in) :: sta
      radau_boot_needed = .true.
      if (allocated(rad_bb)) then
         if (size(rad_bb, 2) == size(sta%q) .and. sta%t /= 0.0d0) radau_boot_needed = .false.
      end if
   end function radau_boot_needed

   ! rad_forget() - clear the radau cross-sequence buffer (member-switch hygiene):
   ! deallocating rad_bb drives radau_boot_needed back to true; the carry counters
   ! reset with it (a switch 3 -> 1/2 -> 3 at t /= 0 must boot afresh, not resume a
   ! stale carry)
   subroutine rad_forget()
      if (allocated(rad_bb)) deallocate (rad_bb, rad_ee, rad_bd, rad_f1, rad_v)
      rad_seq_ok = .false.
      rad_ns = 0
   end subroutine rad_forget

   ! rad_clamp(tp, dir) - apply the dt_min/dt_max MAGNITUDE bounds to one sequence
   ! length (either sign; <=0 bound = unclamped) - no consumed length crosses a bound
   subroutine rad_clamp(tp, dir)
      real(8), intent(inout) :: tp
      real(8), intent(in) :: dir
      if (prop_cfg%dt_max > 0.0d0 .and. abs(tp) > prop_cfg%dt_max) tp = dir*prop_cfg%dt_max
      if (prop_cfg%dt_min > 0.0d0 .and. abs(tp) < prop_cfg%dt_min) tp = dir*prop_cfg%dt_min
   end subroutine rad_clamp

   ! named-abort channel (name the caller; value variant; STOP 1)
   subroutine stop_prop(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': propagator error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the state is not advanced)'
      stop 1
   end subroutine stop_prop

   subroutine stop_propi(who, msg, ival, tail)
      character(len=*), intent(in) :: who, msg, tail
      integer, intent(in) :: ival
      write (0, '(a,a,1x,i0,a,a)') trim(who)//': ', trim(msg), ival, ' ', trim(tail)
      write (0, '(a)') trim(who)//': fatal (the state is not advanced)'
      stop 1
   end subroutine stop_propi

end module propagator
