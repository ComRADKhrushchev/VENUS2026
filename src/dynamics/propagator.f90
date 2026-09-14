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
   ! symple_step(sta) - 4th/6th/8th-order symplectic single step (alternating drift-kick)
   ! Not yet realized - the composition member arrives with its own concretization task;
   ! the legislated sequence is kept as the roadmap:
   ! 1. On first call load the composition coefficients per order (4/6/8, prop_cfg%order)
   !    into a(0:2N) (N = 5/8/17; the symmetric second half mirror-copied) - a
   !    module-private static cache with a first-call flag
   ! 2. Alternating loop i = 0..N-1: force_eval (force before the kick);
   !    p += dt*a(2i)*f (momentum kick); q += dt*a(2i+1)*p/m (position drift)
   ! 3. Final beat: force_eval; p += dt*a(2N)*f
   ! 4. sta%t += dt (the final beat's force site is the endpoint q after the last drift -
   !    the header contract holds with no extra evaluation)
   subroutine symple_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
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
