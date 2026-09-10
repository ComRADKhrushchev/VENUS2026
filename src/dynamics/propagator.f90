!=====================================================================
! propagator.f90 - L3 closed family: nuclear propagators (verlet/symple/radau three-member registration +
!                   unified dispatch stepping - the integrator body holds no system knowledge;
!                   a pure integrator: frame output and termination are the driver's affair, E14)
! Reference: src_VENUS/VERLET.f90:31-79 (velocity Verlet: force first :34-41, full position step
!       :52-56, force re-evaluation + half step :60-67); :90-170 (third-order Beeman: position predictor
!       :98-103, velocity predictor :105-109, force staging and re-evaluation :111-125, velocity
!       corrector :163-168);
!       src_VENUS/SYMPLE.f:14-86 (Schlier 4th/6th/8th-order coefficient loading), :88-100 (alternating
!       drift-kick main loop); src_VENUS/RADAU.f:90-94 (Radau nodes hh), :138-175
!       (w/u weights and triangle coefficients c/d/ra), :183-195 (6 iterations for the first sequence),
!       :197-311 (substep predictor-corrector), :313-336 (restart test and sequence advance), :390-396
!       (adaptive step tp=dir*tcorr*(ss/hv)^(1/9)); RADOUT :438-498 (inlined orchestration - not carried
!       into this module; the GWRITE/TEST call sites become the driver's four-phase loop, E14 spec §8);
!       Excision list (do not carry over when referencing): VERLET.f90:44-50/:71-79/:151-161 (ghost Langevin),
!       :130-148 (GLO/LDFA inlined forces), :179-270 (FRICFORCE/DIAG2 chain) - stochastic laws always go
!       through force_facade member aggregation (P-D1/D4); the 'CALLING RADOUT' debug prints (RADAU.f:112/
!       :380/:428) and the PRINFO comment remnants are not carried over (G4)
!
! -- E14 hooks withdrawal (user rulings 2026-09-08; spec §8) --
!   hooks_t / hook_out_i / hook_term_i / the second prop_step argument are deleted: the propagator
!   reverts to a pure integrator. The E7 clause "inside the integration loop the only direct call is
!   prop_step" is withdrawn - the loop body directly calls prop_step and the facade phase entries
!   (elec_prop / mqc_statistic) plus rec_frame / pred_eval, the four-phase orchestration being explicit
!   in the driver (prop_step -> elec_prop -> mqc_statistic -> rec_frame -> pred_eval). Frame output
!   and termination move to the driver; this also repairs the latent timing defect of the hooks form -
!   a hooks%term checkpoint fired after nuclear propagation but before the electronic phases, so a
!   termination there skipped the step's electronic evolution and left an inconsistent state.
!
! -- endpoint-force contract (spec §3 - the propagator member obligation, legislated here) --
!   prop_step returns only after the endpoint force has been evaluated: the ws held by elec_facade
!   corresponds to the returned sta%q. Under MQC the velocity-Verlet right-end force cannot be folded
!   into the next step's left end - the amplitude propagation of phase 2 (elec_prop) needs the endpoint
!   H - and between the two force sites the electronic phases run (a hop may change the held
!   occupation set, E16), so each step evaluates its own left-end force afresh.
!
! -- dt channel note (E9, re-anchored by E14; spec §3) --
!   dt is no longer a force_eval argument - force_eval(sta, e_tot) is the pure force entry (two
!   dummies; the E9 g_ij(dt) channel belongs to elec_prop's dt dummy, the electronic phase). The
!   member's step is the fixed dt from config (prop_cfg%dt): verlet/symple take it directly; radau
!   seeds its adaptive sequence from it (substep sites are interior - moot now that the force entry
!   carries no dt). The multi-surface x radau conflict stays caught by the method-package member
!   itself (E13, carried: a fixed-step-only law loud-stops on dt semantics it cannot live with -
!   adaptive sequence boundaries have no well-defined "one step" dt; no facade table; same source as
!   the force_facade 4'c note), so the multi-surface path only ever sees fixed steps; this note pins
!   the adiabatic-path behavior to stay unchanged
!=====================================================================
module propagator
   use state,         only: state_t
   use config,        only: prop_cfg => propagator   ! alias to avoid clashing with the module name
                                                      ! (prop_cfg is the config instance 'propagator')
   use config_roster, only: roster
   use force_facade,  only: force_eval               ! T4 rename carried: the pure force entry
                                                      ! (sta, e_tot) - no dt dummy, see the header note
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

   ! Closed-family registry (each of the three members registers one row at assembly time; capacity has slack)
   integer, parameter :: max_integ = 8        ! registry capacity [row]
   type :: integ_entry_t
      integer :: id = 0                       ! integrator code (the key of config%propagator%id)
      procedure(integ_proc_i), pointer, nopass :: proc => null() ! member stepping routine
   end type integ_entry_t
   type(integ_entry_t) :: reg_tbl(max_integ)  ! registry (module-private)
   integer :: n_reg = 0                       ! registered row count [row]
contains
   !------------------------------------------------------------------
   ! prop_reg(id, proc) - register an integrator member into the closed-family registry
   ! Reference: docs/final_completion_plan.md §7 X4 (unified integrator protocol + three-member registration)
   subroutine prop_reg(id, proc)
      integer, intent(in) :: id           ! integrator code [-] (1=verlet/2=symple/3=radau -
                                          ! the key of config%propagator%id)
      procedure(integ_proc_i) :: proc     ! member stepping routine
      ! Pseudocode:
      ! 1. Duplicate id or full table -> loud-stop (the closed family has exactly three members, each
      !    registered once at assembly time)
      ! 2. Register the row {id, proc} (executed in the driver's assembly phase - E7)
   end subroutine prop_reg

   !------------------------------------------------------------------
   ! prop_step(sta) - step dispatch: call the registered member per config%propagator%id
   ! Reference: src_VENUS/RADAU.f:438-498 (the old inline orchestration - now the driver's four-phase
   !       loop, E14 spec §8; this module keeps only the stepping)
   subroutine prop_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (advanced one step/sequence by the member)
      ! Pseudocode:
      ! 1. Look up reg_tbl by prop_cfg%id (unregistered id -> loud-stop naming the id)
      ! 2. Call the member proc(sta): advance sta by prop_cfg%dt (verlet/symple one fixed step, radau
      !    one adaptive sequence; the member calls force_eval as often as it needs) - the member
      !    returns only after the endpoint force has been evaluated (the header contract: the held ws
      !    corresponds to the returned sta%q)
   end subroutine prop_step

   !------------------------------------------------------------------
   ! verlet_step(sta) - velocity Verlet (order=1) / third-order Beeman (order=2) single step
   ! Reference: src_VENUS/VERLET.f90:31-79 (VV); :90-170 (Beeman)
   subroutine verlet_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      ! Local work arrays are declared in C6 (old FT2/QTP1 force and coordinate staging - force
      ! differencing is state of the Beeman branch)
      ! Pseudocode:
      ! 1. dt = prop_cfg%dt, mass = roster%mass, order = prop_cfg%order
      ! 2. order=1 (velocity Verlet, old :31-79):
      ! 2a. force_eval(sta, e) (force F(t) at step start - old :34-41 force site)
      ! 2b. p += f*dt/2 (half-step momentum, old :39-41)
      ! 2c. q += (p/m)*dt (full position step, old :52-56)
      ! 2d. force_eval(sta, e) (force F(t+dt) at step end - old :60-62 force site)
      ! 2e. p += f*dt/2 (completing the half step, old :65-67)
      ! 3. order=2 (third-order Beeman, old :90-170):
      ! 3a. position predictor: q += (p/m)*dt + dt^2*(4f(t)-f(t-dt))/(6m) (old :98-103)
      ! 3b. momentum predictor: p += dt*(3f(t)-f(t-dt))/2 (old :105-109)
      ! 3c. stage force f(t)->f(t-dt) (old :111-114); force_eval(sta, e) (force re-evaluation, old :123-125)
      ! 3d. momentum corrector: p = m*(q-q_prev)/dt + dt*(2f(t+dt)+f(t))/(6m)
      !     (old :163-168; q_prev is the q before 3a, staged in the old QTP1)
      ! 4. The old ghost Langevin (:44-50/:71-79/:151-161) and GLO/LDFA inlined forces
      !    (:130-148, FRICFORCE :179-270) are not carried over - stochastic laws go through force_facade
      !    member aggregation (P-D1/D4); the integrator body holds no system knowledge
      ! 5. sta%t += dt (return: the last force site - 2d / 3c - is the endpoint, the header contract
      !    holds with no extra evaluation)
   end subroutine verlet_step

   !------------------------------------------------------------------
   ! symple_step(sta) - 4th/6th/8th-order symplectic single step (Schlier coefficients, alternating drift-kick)
   ! Reference: src_VENUS/SYMPLE.f:14-86 (coefficient loading), :88-100 (main loop)
   subroutine symple_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      ! The coefficient cache is declared in C6 (old A(0:34) + FIRST flag - module-private statics)
      ! Pseudocode:
      ! 1. On first call load the Schlier coefficients per order (4/6/8, config%propagator%order)
      !    into a(0:2N) (N=5/8/17; symmetric second half mirror-copied, the FIRST block of old :14-86)
      ! 2. Alternating loop i=0..N-1 (old :88-97):
      ! 2a. force_eval(sta, e) (force before the kick - old :94-96 force site)
      ! 2b. p += dt*a(2i)*f (momentum kick, old :91)
      ! 2c. q += dt*a(2i+1)*p/m (position drift, old :92)
      ! 3. Final beat: force_eval(sta, e); p += dt*a(2N)*f (old :98-100)
      ! 4. sta%t += dt (return: the final beat's force site is the endpoint q after the last drift -
      !    the header contract holds with no extra evaluation)
   end subroutine symple_step

   !------------------------------------------------------------------
   ! radau_step(sta) - one Everhart RA15 Gauss-Radau sequence (adaptive step)
   ! Reference: src_VENUS/RADAU.f:90-94 (nodes hh), :138-175 (w/u weights and c/d/ra),
   !       :183-195 (6 iterations for the first sequence + force at sequence start), :197-311 (substep
   !       predictor-difference corrector),
   !       :313-336 (restart test and sequence advance), :390-396 (next-sequence step);
   !       RADOUT :438-498 (orchestration - now the driver's four-phase loop, E14 spec §8)
   subroutine radau_step(sta)
      type(state_t), intent(inout) :: sta   ! physical state (q/p/t advanced; f is the force workspace)
      ! Caches for the constants and B/G/E coefficient arrays are declared in C6 (module-private statics -
      ! carried across sequences)
      ! Pseudocode:
      ! 1. On first call load the constants: Radau nodes hh(2:8) (old :92-94), w/u weights (old :138-143),
      !    triangle coefficients c/d/ra (old :153-175)
      ! 2. Fold velocities v = p/m (old :103-107); iteration count 6 for the first sequence, 2 afterwards
      !    (old :185, :419)
      ! 3. Substep loop within the sequence (iterations x 7 substeps, old :197-311):
      ! 3a. Predict substep coordinates yy from the seventh-order B-coefficient expansion (old :220-245)
      ! 3b. force_eval evaluated at substep coordinates (old :250-260 force site; f_j = f/m fold)
      ! 3c. Difference correction: divided differences of the g coefficients + triangle accumulation of the
      !    B coefficients (old :262-311)
      ! 4. Sequence close-out: advance q/p by the sequence length tt (old :329-336); last-term hv estimate
      !    (old :313-317); restart test - if tp/tt<=1 discard the sequence and restart (old :319-327)
      ! 5. Adaptive step: tp = tcorr*(ss/hv)^(1/9) (old :390; the 1/9 exponent is old :74;
      !    ss = 10^(-order) is the accuracy parameter, direction dir follows the sign of tf - old :133-135)
      ! 6. Step clamping dt_min/dt_max - deferred fields of prop_cfg, introduced in C6
      !    (K1 plan §1.3 deferred note; until then accuracy is controlled by ss alone)
      ! 7. sta%t += sequence length tt
      ! 8. Endpoint-force close-out: force_eval(sta, e) at the returned sta%q - the substep force sites
      !    (3b) sit at interior node coordinates, not the sequence endpoint, so the member closes with
      !    an explicit endpoint evaluation (the header contract: the held ws corresponds to the
      !    returned q); the old NIP interval output (old :99-101 prtime) left with the hooks - frame
      !    output timing is the driver/recorder affair (E14 spec §8)
   end subroutine radau_step
end module propagator
