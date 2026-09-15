!=====================================================================
! check_facade.f90 - facade-mechanism criteria program, batches 1-4:
!                    the frozen ws holder association (a)-(d) + the
!                    force_eval fill chain (e)-(g) + the electronic-phase
!                    dispatch batch (h)-(l) + the reference-container
!                    batch (m)-(p) (criteria defined FIRST in this
!                    header; this program executes them)
! Compile (from src/; run make FIRST - the aggregated registry module is
! a build product; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 &
!            globals/config_roster.f90 globals/state.f90 facade/elec_ws.f90 &
!            facade/elec_facade.f90 methods/elec_adiabatic.f90 &
!            methods/reg_adiabatic.f90 utils/linalg.f90 facade/force_facade.f90 &
!            facade/sampler.f90 utils/geometry.f90 dynamics/propagator.f90 &
!            methods/dist_boltzmann.f90 methods/dist_ebk.f90 methods/dist_incident.f90 &
!            methods/dist_j.f90 methods/dist_rotation.f90 methods/dist_surface_oscillator.f90 &
!            methods/dist_normalmode.f90 methods/dist_barrier_excitation.f90 &
!            methods/dist_thermostat.f90 methods/dist_glo_target.f90 &
!            methods/dist_thermalize.f90 methods/seam_cast.f90 &
!            methods/reg_prop.f90 &
!            methods/reg_dist_boltzmann.f90 methods/reg_dist_ebk.f90 &
!            methods/reg_dist_incident.f90 methods/reg_dist_j.f90 &
!            methods/reg_dist_rotation.f90 methods/reg_dist_surface_oscillator.f90 &
!            methods/reg_dist_normalmode.f90 methods/reg_dist_barrier_excitation.f90 &
!            methods/reg_dist_thermostat.f90 methods/reg_dist_glo_target.f90 &
!            methods/reg_dist_thermalize.f90 output/final_state.f90 &
!            ../systems/reference/reference_pot.f90 &
!            ../systems/reference/reference_box.f90 &
!            ../systems/reference/reg_reference.f90 &
!            build/gfortran/registry_gen.f90 ../scripts/check_facade.f90 &
!            -o build/check_facade.exe
! All artifacts are written under build/ - the source tree stays clean.
! Posture: criteria (a)-(d) exercise the holder mechanism itself - the
!   REAL 'adiabatic' package registers through its reg module
!   (reg_adiabatic_all, the same row the aggregated registry executes)
!   and binds the name, which seeds the statistics holder and makes the
!   occ_get channel live; the members are not dispatched by that batch.
!   Criteria (e)-(g) additionally
!   assemble a FORCE CHAIN check-side: a local analytic branch container
!   is bound through container_bind_pes (the public assembly face), so
!   every (e)/(g) assertion runs on real force_eval calls; the default
!   binding posture of the tree is untouched (the driver binds nothing).
!   The bound container branches on a check-side fill selector: 1 = the
!   scalar harmonic channel, 2 = a two-surface diabatic channel (linear
!   H_ij + dh_ij + an additive scalar member), 3 = the one-surface
!   diabatic twin of branch 1 (same V and gradient delivered as 1x1
!   constructs, additive channel absent) - every branch is a local
!   procedure of this check; none of them ever enters src/. Criteria
!   (h)-(j) dispatch the REAL registered 'adiabatic' row through the two
!   electronic-phase entries on that same bound force chain (real calls,
!   not stub invocations). Loud-stop verification approach (same scheme
!   as the sibling checks): criteria (k)/(l) run as CHILD PROCESSES - the
!   check re-spawns itself with a mode argument (execute_command_line,
!   stderr redirected to build/); the child drives the real stop path and
!   must print a naming message on stderr and exit nonzero. The batch-4
!   real-container criteria ((m)/(n)) run as children for a second
!   reason: the parent's single PES slot is already bound check-side for
!   (e)-(j) (one container per process - the bind face loud-stops on a
!   second bind), so every different container gets a fresh process; the
!   (o) children run with cwd = build/ because the container parameter
!   file is run-directory material.
! Criteria:
!   (a) handle identity: two ws_handle() calls yield pointers to the SAME
!       target (pointer association comparison), and a field written
!       through one handle is visible through the other handle AND through
!       the ws_get snapshot channel (the alias is onto the holder itself,
!       not onto a copy)
!   (b) copy-out semantics: ws_get returns a snapshot - scribbling the
!       snapshot (a scalar component and an array component) leaves the
!       holder intact, and a fresh ws_get reproduces the holder's
!       pre-scribble content
!   (c) reset semantics: with the holder dirtied (repr set, scalar set,
!       every per-mode array allocated and filled), ws_reset returns it to
!       the DEFINED empty state - a subsequent ws_get shows repr = r_empty,
!       every per-mode array unallocated, the scalar zeroed (the one-step
!       survival contract: representation constructs live one integration
!       step, then die - not merely become unreadable)
!   (d) snapshot-style vs handle-style coexistence: the existing consumer
!       pattern of the scalar force branch (occ_get request -> ws_get
!       snapshot -> compose into the snapshot) runs unchanged and leaves
!       the HOLDER untouched (zero interference), while the handle-style
!       compose (the fill-site posture that writes the holder itself)
!       reads back through the snapshot channel without retroactively
!       perturbing an earlier snapshot. The regression half of this
!       criterion is external: check_propagator (the force_eval scalar
!       branch end-to-end on real forces) and check_driver (the full
!       four-phase loop) re-run green on the rewritten holder.
!       Historical note (d1): when (d) was written the production fill
!       site still composed into a ws_get snapshot; with the force_eval
!       realization the production fill went handle-first - the
!       snapshot-consumer pattern validated here remains a legal
!       consumption form (ws_get stays a public copy-out primitive), but
!       it no longer describes the fill site itself
!   (e) end-to-end fill chain through the handle: with the analytic
!       container bound through the public face, force_eval returns the
!       analytic force component-wise at two geometries (proportionality
!       to the gradient - the conversion constant cancels) and e_tot
!       decomposes as t_kin(p) + v_pot(q) via the composed conversion
!       measured by one calibration call (the eV->internal literal lives
!       only in force_facade - this check NEVER writes it); the container
!       consumed the occupation request (i_want = [1] observed inside the
!       container), and AFTER the call the holder - read back through a
!       fresh ws_handle() alias AND a ws_get snapshot - carries the
!       composed constructs of the returned geometry (repr, the request
!       field, g_add, v_add): the fill reaches the HOLDER itself and
!       SURVIVES the call (the one-step conveyance to the electronic
!       phases)
!   (f) conversion single-site: a recursive scan of every src/ .f90 file
!       (findstr on a Windows host, grep elsewhere) lists the files whose
!       text contains the conversion literal; each listed file is then
!       scanned line by line, and exactly ONE occurrence of the literal
!       sits at EXECUTABLE position (not inside a comment) in the whole
!       tree - the ev_kcal parameter declaration of force_facade.f90
!       (the unit-contract grep requirement, executable face)
!   (g) the multi-surface diabatic fold: with the container delivering a
!       two-surface linear H_ij + dh_ij + an additive scalar member,
!       force_eval returns the closed-form eigen-fold - the occupied
!       eigenvalue's gradient u^T (dH/dq) u with u the algebraic
!       eigenvector of the LOWER eigenvalue (the ascending occupation
!       row) plus the additive kernel, component-wise (rel 1e-12), and
!       e_tot = t_kin + (lambda_occ + v_add) via the measured conversion
!       (rel 1e-12); the delivered h_ij SURVIVES in the holder (the
!       facade diagonalizes a copy - the constructs stay for the
!       electronic phases). Degenerate equivalence: the one-surface diab
!       twin vs the scalar channel at the same state - force and e_tot
!       compared BITWISE (the attempt), falling back to rel 1e-15 with a
!       printed note
!   (h) electronic-phase dispatch reachable + no-op-safe: with the real
!       'adiabatic' package bound and the holder filled by a real
!       force_eval call (plus a superset dirty - the five arrays the
!       scalar fill does not deliver, allocated through the alias), the
!       pair elec_prop(sta, dt) / mqc_statistic(sta) REALLY DISPATCHES
!       through the table without crashing and leaves the physical state
!       bitwise unchanged (q/p/f/t, n_surf, the ontic containers a/rho
!       staying unallocated); elec_prop additionally leaves the holder
!       INTACT - content equality on g_add/v_add/i_want plus survival of
!       the dirtied arrays (the step-end reset belongs to the
!       mqc_statistic tail ALONE; an early reset fails here)
!   (i) step-end reset at the mqc_statistic tail: after the dirty holder
!       of (h) went through both dispatches, a snapshot reads back the
!       DEFINED empty state - repr = r_empty, every one of the seven
!       per-mode arrays gone (not merely the delivered subset), the
!       scalar zeroed; a second fill + dispatch cycle clears again
!   (j) statistics write-back channels: across the dispatch cycles
!       occ_get hands back [1] CONSTANT (the degenerate occupation the
!       match member leaves untouched) and n_hop_get 0 (the run-lifetime
!       accumulator holding its degenerate value across steps - the
!       step-end reset of (i) never touches the statistics holder, the
!       two holders have different lifetimes by statute)
!   (k) H-form match loud-stop (child): a check-side package demanding
!       h_form 'diab' bound against a container declaring h_deliver
!       'scalar' -> container_init loud-stops, stderr naming
!       container_init and BOTH sides of the failed subset-of match
!   (l) unbound-method posture (child): elec_prop / mqc_statistic called
!       before any method_bind loud-stop, stderr naming the entry and
!       method_bind (the skeleton's legislation: unbound = loud-stop,
!       never a silent skip)
!   (m) container registry path (child): the harmonic well bound through
!       its OWN drag-in row - the build-generated registry_gen module,
!       the exact registration call the driver runs - NOT the check-side
!       direct bind of (e). The child runs registry_gen_all (which
!       registers every method package AND binds the reference PES with
!       delivery declared 'scalar'), registry_gen_init (the assembly
!       seam, guard-skipped without staged input), container_init
!       (green) and then a real force_eval: the force components and the
!       e_tot decomposition of the default-constant well (K = 0.04
!       eV/A^2, the numeric anchor) at two geometries through the
!       measured conversion (1e-13 / 1e-12; the eV->internal literal
!       never written here)
!   (n) module direct bind (children): the morse and leps potentials of
!       the reference module, bound through container_bind_pes by
!       transparent check-side adapters (the direct-bind counterpart of
!       (m)'s registry path) and driven by a real force_eval: the
!       analytic gradient vs the CENTRAL DIFFERENCE of the module energy
!       (tolerance 1e-8, step h = 1e-6, component-wise - the analytic-
!       gradient criterion family), the facade fold reproducing the
!       module gradient through the measured conversion (1e-12), and
!       the e_tot decomposition at the same geometry with momenta
!   (o) box parameter file (children, cwd = build/): the container's
!       parameter channel reference_pars.txt read from the RUN directory:
!       the ABSENT file leaves the default constant and both optional
!       slots unbound (k = default, term = F, class = F); the GOOD file
!       (K_HARM/Q_TERM/E_SPLIT) stages the constant and provides both
!       optional slots (k, term = T, class = T - asserted on the child's
!       printed marker); the BAD files (unknown keyword, duplicate
!       keyword, unreadable value, non-positive constant) loud-stop
!       naming reference_load and carrying the file name + the keyword
!   (p) the leps diabatic full-matrix delivery (module face, in
!       process): leps_diab at a fixed geometry delivers H symmetric and
!       EVERY dh face symmetric (bitwise, both triangles); the analytic
!       faces (1,1)/(1,2) match the central difference of the delivered
!       H entries (1e-8, h = 1e-6); the eigen identity - eig_sym of the
!       delivered H yields lambda_1 = the scalar leps_vg energy (1e-12) -
!       and the eigenvector-fold gradient sum_ij u_i u_j (dH/dq)_ij
!       computed through BOTH delivered triangles = the scalar leps_vg
!       gradient (1e-12): the n_surf = 1 diab-vs-scalar degenerate
!       equivalence, made executable at the module face
!=====================================================================
program check_facade
   use consts,        only: e_conv
   use state,         only: state_t, state_create, state_destroy
   use config,        only: electronic
   use config_roster, only: roster
   use elec_ws,       only: ws_t, r_empty, r_scalar, r_diab
   use elec_facade,   only: ws_handle, ws_get, ws_reset, occ_get, n_hop_get, &
                            method_reg, method_bind, elec_prop, mqc_statistic
   use force_facade,  only: force_eval, container_bind_pes, container_init
   use reg_adiabatic, only: reg_adiabatic_all
   use linalg,        only: eig_sym
   use registry_gen,  only: registry_gen_all, registry_gen_init
   use reference_pot, only: morse_vg, leps_vg, leps_diab
   use reference_box, only: reference_load, reference_k, reference_term_on, &
                            reference_class_on
   implicit none
   type(ws_t) :: w, w2                     ! snapshot instances (ws_get dummies)
   type(ws_t), pointer :: p1, p2           ! handle aliases
   integer, allocatable :: occ_req(:)      ! the delivery request (occ_get)
   integer :: nfail = 0
   ! ---------- fill-chain section state ((e)-(g)) ----------
   type(state_t) :: sa                     ! the force-chain state (2 atoms / 6 dofs)
   integer :: nfail0                       ! per-criterion failure watermark
   integer :: m                            ! fold loop index
   integer :: fill_sel = 0                 ! the bound container's fill-mode selector
   integer :: seen_nwant = 0               ! the request the container observed (length)
   integer :: seen_want1 = 0               ! the request the container observed (first entry)
   integer :: nexec, nhit                  ! (f) executable-position hit counters
   integer :: v, ios, jj                   ! (f) scan unit / status / scratch
   character(len=512) :: path              ! (f) scan path buffer
   character(len=512) :: exec_path         ! (f) the file that carried the last hit
   character(len=512) :: cmd, ln           ! (f) scan command / line buffer
   logical :: bitexact                     ! (g) the twin's bitwise attempt
   real(8) :: cal, ea, eb, tk              ! the measured conversion + energy channels
   real(8) :: lam, un, u1, u2              ! (g) the closed-form eigenpair
   real(8) :: fs(6), es, et                ! (g) the twin comparison carriers
   real(8), allocatable :: gref(:)         ! (g) the manual-diagonalization reference fold
   ! harmonic test-potential parameters (per-dof, 6 dofs / 2 atoms of mass 1 and 4 amu)
   real(8), parameter :: keV(6)   = [ 0.95d0, 1.15d0, 0.85d0, 1.25d0, 0.75d0, 1.10d0 ] ! stiffness [eV/A^2]
   real(8), parameter :: qe_ctr(6) = [ 0.12d0,-0.18d0, 0.28d0,-0.14d0, 0.22d0,-0.06d0 ] ! centers [A]
   ! criterion-(e) geometries / momenta
   real(8), parameter :: ueA(6) = [ 0.30d0,-0.21d0, 0.18d0,-0.12d0, 0.24d0,-0.09d0 ] ! calibration offset [A]
   real(8), parameter :: ueB(6) = [-0.16d0, 0.26d0,-0.09d0, 0.31d0,-0.20d0, 0.14d0 ] ! second geometry [A]
   real(8), parameter :: pB(6)  = [ 1.05d0,-0.65d0, 0.48d0,-0.85d0, 0.55d0, 1.20d0 ] ! second-geometry momenta
   ! criterion-(g) diab branch: the linear two-surface H_ij (u = q - qe_ctr)
   real(8), parameter :: hb1 = -0.30d0, hb2 = 0.22d0, hb12 = 0.10d0 ! constant H offsets [eV]
   real(8), parameter :: a1(6) = [ 0.80d0,-0.50d0, 0.60d0,-0.40d0, 0.70d0,-0.30d0 ] ! dH11/du [eV/A]
   real(8), parameter :: a2(6) = [-0.40d0, 0.70d0,-0.30d0, 0.80d0,-0.50d0, 0.60d0 ] ! dH22/du [eV/A]
   real(8), parameter :: cc(6) = [ 0.25d0,-0.35d0, 0.15d0,-0.45d0, 0.30d0,-0.20d0 ] ! dH12/du [eV/A]
   real(8), parameter :: g_uni(6) = [ 0.05d0,-0.08d0, 0.06d0,-0.04d0, 0.07d0,-0.03d0 ] ! additive kernel [eV/A]
   real(8), parameter :: v_uni = -0.35d0   ! additive energy [eV]
   real(8), parameter :: ug(6) = [ 0.26d0,-0.19d0, 0.23d0,-0.11d0, 0.27d0,-0.15d0 ] ! (g) geometry [A]
   real(8), parameter :: pg(6) = [ 0.20d0,-0.13d0, 0.17d0,-0.10d0, 0.15d0, 0.22d0 ] ! (g) momenta
   ! ---------- dispatch-batch locals ((h)-(l)) ----------
   character(len=32) :: mode               ! the child spawn's mode argument
   type(state_t) :: st1                    ! the child-mode minimal state
   real(8) :: qref(6), pref(6), fref(6)    ! the state-invariance references
   real(8) :: tref                         ! (ditto - the time member)
   real(8), parameter :: dtv = 0.05d0      ! the dt channel value [10 fs] fed to elec_prop
                                           ! (the fixed staged step's stand-in; the degenerate
                                           ! member ignores it - the criterion needs the channel
                                           ! CARRIED, not consumed)

   ! ---------- child dispatch: one loud-stop case per spawn (criteria (k)/(l)) ----------
   ! the child drives the real stop path; if the chain ever RETURNS, the child
   ! prints the bug witness and exits 0 - the parent's exitstat check fails
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      if (trim(mode) == 'mismatch') then
         ! (k): a check-side package demanding the diab H form, bound against a
         ! container declaring only the scalar delivery - container_init's
         ! subset-of match must loud-stop naming both sides
         call method_reg('chk_diab', 'diab', noop_prop, noop_match)
         electronic%method = 'chk_diab'
         call container_bind_pes(check_pes, 'scalar')
         call container_init()
      else if (trim(mode) == 'unbound-prop') then
         ! (l): elec_prop before any method_bind - no row bound, loud-stop
         call state_create(st1, 1)
         call elec_prop(st1, dtv)
         call state_destroy(st1)
      else if (trim(mode) == 'unbound-match') then
         ! (l): mqc_statistic before any method_bind - same posture
         call state_create(st1, 1)
         call mqc_statistic(st1)
         call state_destroy(st1)
      else if (trim(mode) == 'mbind') then
         ! (m): the reference container through its own drag-in row (the
         ! registry_gen path the driver assembles) - the child terminates
         ! itself (marker on success, FAIL + stop 1 on any assertion)
         call child_mbind()
      else if (trim(mode) == 'nmorse' .or. trim(mode) == 'nleps') then
         ! (n): one module potential direct-bound per child (fresh process -
         ! one container slot per process)
         call child_npot(trim(mode) == 'nmorse', trim(mode))
      else if (trim(mode) == 'obox') then
         ! (o): the parameter-file child (cwd = build/ - the file the parent
         ! staged there); good posture prints the marker, bad files loud-stop
         call child_obox()
      end if
      print '(a)', 'child: NO LOUD-STOP (the child chain returned - bug)'
      stop 0
   end if

   ! the real 'adiabatic' package registers through its reg module (the same
   ! call the aggregated registry executes) + binds the name, which seeds the
   ! statistics holder - occ_get is live for criterion (d); batches (a)-(d)
   ! never dispatch the members, batch (h)-(j) dispatches them for real
   call reg_adiabatic_all()
   call method_bind('adiabatic')

   ! ---------- (a) handle identity + write-through observability ----------
   call ws_reset()
   p1 => ws_handle()
   p2 => ws_handle()
   allocate (p1%g_add(3))
   p1%g_add = [1.0d0, 2.0d0, 3.0d0]
   p1%v_add = 3.25d0
   call ws_get(w)
   if (.not. (associated(p1) .and. associated(p2) .and. associated(p1, p2)) .or. &
       p2%v_add /= 3.25d0 .or. any(p2%g_add /= [1.0d0, 2.0d0, 3.0d0]) .or. &
       w%v_add /= 3.25d0 .or. any(w%g_add /= [1.0d0, 2.0d0, 3.0d0])) then
      print '(a)', 'FAIL (a) handle identity: the two handles do not alias one '// &
                   'target, or a handle write is invisible to the other handle / ws_get'
      nfail = nfail + 1
   else
      print '(a)', '(a) handle identity: PASS (two ws_handle() calls alias one '// &
                   'target; the write through one is seen by the other and by ws_get)'
   end if

   ! ---------- (b) copy-out semantics: the snapshot is not the holder ----------
   w%v_add = -999.0d0                      ! scribble on the snapshot taken above
   w%g_add = -999.0d0                      ! (a scalar and an array component)
   call ws_get(w2)                         ! a fresh snapshot of the holder
   if (w2%v_add /= 3.25d0 .or. any(w2%g_add /= [1.0d0, 2.0d0, 3.0d0])) then
      print '(a)', 'FAIL (b) copy-out: writing the snapshot perturbed the holder'
      nfail = nfail + 1
   else
      print '(a)', '(b) copy-out: PASS (scribbling the snapshot leaves the holder '// &
                   'intact; a fresh ws_get reproduces the held content)'
   end if

   ! ---------- (c) reset semantics: full clear to the defined empty state ----------
   call ws_reset()                         ! clean slate so every allocation below is fresh
   p1 => ws_handle()
   p1%repr = r_diab
   p1%v_add = -777.0d0
   allocate (p1%i_want(2))
   p1%i_want = [1, 2]
   allocate (p1%g_add(3))
   p1%g_add = 1.0d0
   allocate (p1%h_ij(2, 2))
   p1%h_ij = 1.0d0
   allocate (p1%dh_ij(3, 2, 2))
   p1%dh_ij = 1.0d0
   allocate (p1%v_surf(2))
   p1%v_surf = 1.0d0
   allocate (p1%g_surf(3, 2))
   p1%g_surf = 1.0d0
   allocate (p1%d_ij(3, 2))
   p1%d_ij = 1.0d0
   call ws_reset()
   call ws_get(w)
   if (w%repr /= r_empty .or. w%v_add /= 0.0d0 .or. allocated(w%i_want) .or. &
       allocated(w%g_add) .or. allocated(w%h_ij) .or. allocated(w%dh_ij) .or. &
       allocated(w%v_surf) .or. allocated(w%g_surf) .or. allocated(w%d_ij)) then
      print '(a)', 'FAIL (c) reset: the holder did not return to the defined empty '// &
                   'state (repr, scalar, or a per-mode array survived)'
      nfail = nfail + 1
   else
      print '(a)', '(c) reset: PASS (repr = r_empty, all seven per-mode arrays '// &
                   'gone, scalar zeroed - the constructs die at step end)'
   end if

   ! ---------- (d) snapshot-style vs handle-style coexistence ----------
   call ws_reset()
   call occ_get(occ_req)
   ! the snapshot fill pattern of the scalar force branch, verbatim: request
   ! through occ_get, snapshot through ws_get, compose into the snapshot
   call ws_get(w)
   if (allocated(w%i_want)) deallocate (w%i_want)
   allocate (w%i_want(size(occ_req)))
   w%i_want = occ_req
   allocate (w%g_add(3))
   w%g_add = 0.5d0
   w%v_add = -1.5d0
   w%repr = r_scalar
   ! zero interference: the snapshot fill must leave the holder in its empty
   ! state (the pre-freeze consumer behavior, unchanged by the rewrite)
   p1 => ws_handle()
   if (size(occ_req) /= 1 .or. occ_req(1) /= 1 .or. &
       p1%repr /= r_empty .or. p1%v_add /= 0.0d0 .or. allocated(p1%i_want) .or. &
       allocated(p1%g_add)) then
      print '(a)', 'FAIL (d) coexistence: the occ channel broke, or the snapshot '// &
                   'fill trespassed on the holder'
      nfail = nfail + 1
   else
      print '(a)', '(d) coexistence part 1: PASS (occ request seeds [1]; the '// &
                   'snapshot fill leaves the holder empty - zero interference)'
   end if
   ! the handle-style compose (the fill-site posture onto the frozen holder):
   ! writes the holder itself, then reads back through the snapshot channel
   p1%repr = r_scalar
   allocate (p1%i_want(size(occ_req)))
   p1%i_want = occ_req
   allocate (p1%g_add(3))
   p1%g_add = 0.25d0
   p1%v_add = -2.5d0
   call ws_get(w2)
   if (p1%repr /= r_scalar .or. p1%v_add /= -2.5d0 .or. any(p1%g_add /= 0.25d0) .or. &
       w2%repr /= r_scalar .or. w2%v_add /= -2.5d0 .or. any(w2%g_add /= 0.25d0) .or. &
       any(w2%i_want /= 1) .or. &
       w%v_add /= -1.5d0 .or. any(w%g_add /= 0.5d0)) then
      print '(a)', 'FAIL (d) coexistence: the handle read-back broke, or the earlier '// &
                   'snapshot was retroactively perturbed'
      nfail = nfail + 1
   else
      print '(a)', '(d) coexistence part 2: PASS (the handle write reads back '// &
                   'through ws_get; the earlier snapshot keeps its own content)'
   end if
   call ws_reset()                         ! leave the holder in the empty state

   ! ---------- (e) end-to-end fill chain through the handle ----------
   ! assembly check-side: the roster masses, a 2-atom state, and the local branch
   ! container bound through the PUBLIC assembly face (one slot, bound once - the
   ! (g) fill branches ride the same binding through the fill_sel selector)
   allocate (roster%mass(2))
   roster%mass = [ 1.0d0, 4.0d0 ]
   call state_create(sa, 2)
   call container_bind_pes(check_pes)
   fill_sel = 1
   ! calibration: one force_eval at a known single-dof offset (p = 0) measures the
   ! composed conversion constant cal (the eV->internal product the facade applies;
   ! the literal lives only in force_facade - never written here)
   sa%q = qe_ctr + ueA
   sa%p = 0.0d0
   call force_eval(sa, ea)
   cal = -sa%f(1)/(keV(1)*ueA(1))
   nfail0 = nfail                ! the (e) water level: taken BEFORE the first pre-assertion
                                  ! so any prior fail() suppresses the PASS line below
   if (cal <= 0.0d0) then
      call fail('(e) measured conversion constant non-positive: '//trim(d2s(cal)))
   end if
   ! geometry A: the analytic force component-wise (proportionality - cal cancels
   ! into the same composition on both sides)
   do m = 1, 6
      if (abs(sa%f(m) + (keV(m)*(sa%q(m) - qe_ctr(m)))*cal) > &
          1.0d-13*max(abs((keV(m)*(sa%q(m) - qe_ctr(m)))*cal), 1.0d0)) then
         call fail('(e) analytic-force mismatch at the calibration geometry, dof '//trim(i2s(m)))
      end if
   end do
   ! geometry B: force + the e_tot decomposition via the measured conversion
   sa%q = qe_ctr + ueB
   sa%p = pB
   call force_eval(sa, eb)
   do m = 1, 6
      if (abs(sa%f(m) + (keV(m)*(sa%q(m) - qe_ctr(m)))*cal) > &
          1.0d-13*max(abs((keV(m)*(sa%q(m) - qe_ctr(m)))*cal), 1.0d0)) then
         call fail('(e) analytic-force mismatch, geometry B, dof '//trim(i2s(m)))
      end if
   end do
   tk = 0.0d0
   do m = 1, 2
      tk = tk + dot_product(sa%p(3*m - 2:3*m), sa%p(3*m - 2:3*m))/(2.0d0*roster%mass(m))
   end do
   tk = tk/e_conv
   if (abs(eb - (tk + 0.5d0*dot_product(keV, (sa%q - qe_ctr)*(sa%q - qe_ctr))*cal/e_conv)) > &
       1.0d-12*max(abs(eb), 1.0d0)) then
      call fail('(e) e_tot decomposition mismatch, geometry B: '//trim(d2s(eb)))
   end if
   ! the request handshake: the container must have consumed i_want = [1] (the held
   ! degenerate occupation seeded by method_bind)
   if (seen_nwant /= 1 .or. seen_want1 /= 1) then
      call fail('(e) the container did not observe the occupation request (i_want /= [1])')
   end if
   ! the holder read-back: AFTER the call a fresh alias AND a snapshot carry the
   ! composed constructs of the returned geometry - the fill reached the holder and
   ! survived the call (one-step conveyance); the references recompute the
   ! container's own expressions from sa%q, so equality is exact. Stepwise guards:
   ! size()/element queries run only behind their allocation checks (a broken fill
   ! must FAIL here, not crash). The PASS line is water-level gated: it prints only
   ! when NO (e) assertion above failed (nfail == nfail0, the watermark taken before
   ! the first (e) pre-assertion - a pre-assertion FAIL suppresses the PASS line)
   p1 => ws_handle()
   call ws_get(w2)
   bitexact = (p1%repr == r_scalar .and. allocated(p1%g_add) .and. allocated(p1%i_want))
   if (bitexact) bitexact = (size(p1%g_add) == 6 .and. size(p1%i_want) == 1 .and. &
       p1%i_want(1) == 1 .and. all(p1%g_add == keV*(sa%q - qe_ctr)) .and. &
       p1%v_add == 0.5d0*dot_product(keV, (sa%q - qe_ctr)*(sa%q - qe_ctr)))
   if (.not. bitexact) then
      call fail('(e) holder read-back: the fill did not reach the holder / survive the call')
   else if (w2%repr /= r_scalar .or. .not. allocated(w2%g_add) .or. .not. allocated(w2%i_want) .or. &
       any(w2%g_add /= p1%g_add) .or. w2%v_add /= p1%v_add .or. any(w2%i_want /= 1)) then
      call fail('(e) snapshot read-back: ws_get does not reproduce the held constructs')
   else if (nfail == nfail0) then
      print '(a)', '(e) fill chain through the handle: PASS (analytic force + e_tot at two'// &
         ' geometries; request consumed; the composed constructs sit in the holder after'// &
         ' the call - alias and snapshot agree)'
   end if

   ! ---------- (f) conversion single-site ----------
   ! list every src/ file whose text carries the conversion literal, then count the
   ! EXECUTABLE-position occurrences (comment text does not count); the unit contract
   ! demands exactly one - the facade's conversion-parameter declaration
   nfail0 = nfail
   if (windows_host()) then
      cmd = 'findstr /s /m /l /c:"23.0605" *.f90 > build/fcd_grep.txt 2>&1'
   else
      cmd = "grep -r -l -F --include='*.f90' 23.0605 . > build/fcd_grep.txt 2>&1"
   end if
   call execute_command_line(trim(cmd), wait=.true.)
   nexec = 0
   exec_path = ''
   open (newunit=v, file='build/fcd_grep.txt', status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call fail('(f) the conversion-literal file listing could not be produced')
   else
      do
         read (v, '(a)', iostat=ios) ln
         if (ios /= 0) exit
         path = adjustl(ln)
         jj = len_trim(path)
         if (jj > 0) then
            if (path(jj:jj) == achar(13)) then       ! strip a carried CR terminator
               path(jj:jj) = ' '
               path = adjustl(path)
            end if
         end if
         if (path(1:2) == './') path = adjustl(path(3:))
         if (len_trim(path) == 0) cycle
         call scan_exec_hits(trim(path), nhit)
         if (nhit > 0) then
            nexec = nexec + nhit
            exec_path = trim(path)
         end if
      end do
      close (v)
      if (nexec /= 1 .or. index(exec_path, 'force_facade.f90') == 0) then
         call fail('(f) conversion single-site violated: '//trim(i2s(nexec))// &
            ' executable-position occurrence(s) in src/, last in "'//trim(exec_path)//'"'// &
            ' (want exactly 1, in force_facade.f90)')
      else
         print '(a)', '(f) conversion single-site: PASS (exactly one executable-position'// &
            ' occurrence of the conversion literal in src/: '//trim(exec_path)//')'
      end if
   end if

   ! ---------- (g) the multi-surface diabatic fold ----------
   ! main: the two-surface linear H_ij with an additive scalar member - force_eval
   ! must return the closed-form eigen-fold of the occupied (lowest, ascending)
   ! eigenvalue plus the additive kernel, and the delivered h_ij must SURVIVE in
   ! the holder (the facade diagonalizes a copy)
   fill_sel = 2
   sa%q = qe_ctr + ug
   sa%p = pg
   call force_eval(sa, ea)
   nfail0 = nfail
   lam = 0.5d0*((hb1 + dot_product(a1, ug) + hb2 + dot_product(a2, ug)) &
                - sqrt((hb1 + dot_product(a1, ug) - hb2 - dot_product(a2, ug))**2 &
                       + 4.0d0*(hb12 + dot_product(cc, ug))**2))
   un = sqrt((hb12 + dot_product(cc, ug))**2 + (lam - (hb1 + dot_product(a1, ug)))**2)
   u1 = (hb12 + dot_product(cc, ug))/un
   u2 = (lam - (hb1 + dot_product(a1, ug)))/un
   allocate (gref(6))
   do m = 1, 6
      gref(m) = (u1*u1)*a1(m) + 2.0d0*(u1*u2)*cc(m) + (u2*u2)*a2(m) + g_uni(m)
   end do
   do m = 1, 6
      if (abs(sa%f(m) + gref(m)*cal) > 1.0d-12*max(abs(gref(m)*cal), 1.0d0)) then
         call fail('(g) eigen-fold force mismatch, dof '//trim(i2s(m))//': '// &
            trim(d2s(sa%f(m)))//' vs closed form '//trim(d2s(-gref(m)*cal)))
      end if
   end do
   tk = 0.0d0
   do m = 1, 2
      tk = tk + dot_product(sa%p(3*m - 2:3*m), sa%p(3*m - 2:3*m))/(2.0d0*roster%mass(m))
   end do
   tk = tk/e_conv
   if (abs(ea - (tk + (lam + v_uni)*cal/e_conv)) > 1.0d-12*max(abs(ea), 1.0d0)) then
      call fail('(g) eigen-fold e_tot mismatch: '//trim(d2s(ea))//' vs closed form '// &
         trim(d2s(tk + (lam + v_uni)*cal/e_conv)))
   end if
   p1 => ws_handle()
   bitexact = (p1%repr == r_diab .and. allocated(p1%h_ij) .and. allocated(p1%dh_ij))
   if (bitexact) bitexact = (size(p1%h_ij, 1) == 2 .and. size(p1%h_ij, 2) == 2 .and. &
       all(p1%h_ij == reshape([hb1 + dot_product(a1, ug), hb12 + dot_product(cc, ug), &
                               hb12 + dot_product(cc, ug), hb2 + dot_product(a2, ug)], &
                              [2, 2])))
   if (.not. bitexact) then
      call fail('(g) holder read-back: the delivered h_ij did not survive in the holder')
   else if (nfail == nfail0) then
      print '(a)', '(g) diab fold: PASS (two-surface eigen-fold + additive kernel = closed'// &
         ' form to 1e-12; e_tot = t_kin + (lambda_occ + v_add); h_ij survives in the holder)'
   end if
   deallocate (gref)
   ! twin: the one-surface diab delivery vs the scalar channel at the SAME state -
   ! bitwise attempt (identical arithmetic through the 1x1 eigenpair), 1e-15 fallback
   fill_sel = 1
   sa%q = qe_ctr + ueB
   sa%p = pB
   call force_eval(sa, es)
   fs = sa%f
   fill_sel = 3
   sa%q = qe_ctr + ueB
   sa%p = pB
   call force_eval(sa, et)
   bitexact = (all(sa%f == fs) .and. et == es)
   if (bitexact) then
      print '(a)', '(g) degenerate equivalence: PASS BITWISE (one-surface diab twin ='// &
         ' scalar channel: force and e_tot identical bit for bit)'
   else if (all(abs(sa%f - fs) <= 1.0d-15*max(abs(fs), 1.0d0)) .and. &
            abs(et - es) <= 1.0d-15*max(abs(es), 1.0d0)) then
      print '(a)', '(g) degenerate equivalence: PASS at the 1e-15 fallback (bitwise'// &
         ' attempt failed - arithmetic-order deviation, within tolerance)'
   else
      call fail('(g) degenerate equivalence: one-surface diab twin vs scalar channel'// &
         ' deviates beyond 1e-15 (force max dev '//trim(d2s(maxval(abs(sa%f - fs))))// &
         ', e_tot dev '//trim(d2s(abs(et - es)))//')')
   end if
   ! ---------- (h)/(i)/(j): the electronic-phase dispatch batch ----------
   ! cycle 1: the production fill (a real force_eval - the same chain the driver
   ! runs) + a superset dirty of the five arrays the scalar fill does not deliver,
   ! so the (i) full-clear assertion sees ALL seven arrays, not the delivered
   ! subset; the state-invariance references are taken at the dirty point
   fill_sel = 1
   sa%q = qe_ctr + ueB
   sa%p = pB
   call force_eval(sa, ea)
   p1 => ws_handle()
   allocate (p1%h_ij(2, 2))
   p1%h_ij = 7.0d0
   allocate (p1%dh_ij(6, 2, 2))
   p1%dh_ij = 7.0d0
   allocate (p1%v_surf(2))
   p1%v_surf = 7.0d0
   allocate (p1%g_surf(6, 2))
   p1%g_surf = 7.0d0
   allocate (p1%d_ij(6, 2))
   p1%d_ij = 7.0d0
   qref = sa%q
   pref = sa%p
   fref = sa%f
   tref = sa%t
   nfail0 = nfail                          ! the (h) water level
   call elec_prop(sa, dtv)                 ! REAL dispatch through the bound row
   if (.not. all(sa%q == qref) .or. .not. all(sa%p == pref) .or. .not. all(sa%f == fref) .or. &
       sa%t /= tref .or. sa%s%n_surf /= 1 .or. allocated(sa%s%a) .or. allocated(sa%s%rho)) then
      call fail('(h) elec_prop changed the physical state under the adiabatic no-op law')
   end if
   if (p1%repr /= r_scalar .or. .not. allocated(p1%g_add) .or. .not. allocated(p1%i_want) .or. &
       any(p1%i_want /= 1) .or. .not. all(p1%g_add == keV*(sa%q - qe_ctr)) .or. &
       p1%v_add /= 0.5d0*dot_product(keV, (sa%q - qe_ctr)*(sa%q - qe_ctr)) .or. &
       .not. allocated(p1%h_ij) .or. .not. all(p1%h_ij == 7.0d0) .or. &
       .not. allocated(p1%dh_ij) .or. .not. allocated(p1%v_surf) .or. &
       .not. allocated(p1%g_surf) .or. .not. allocated(p1%d_ij)) then
      call fail('(h) elec_prop disturbed the held ws (early reset or consumption - the'// &
                ' step-end reset belongs to the mqc_statistic tail alone)')
   end if
   call mqc_statistic(sa)                  ! REAL dispatch: the no-op match + the tail reset
   if (.not. all(sa%q == qref) .or. .not. all(sa%p == pref) .or. .not. all(sa%f == fref) .or. &
       sa%t /= tref .or. sa%s%n_surf /= 1 .or. allocated(sa%s%a) .or. allocated(sa%s%rho)) then
      call fail('(h) mqc_statistic changed the physical state under the adiabatic no-op law')
   else if (nfail == nfail0) then
      print '(a)', '(h) dispatch reachable + no-op-safe: PASS (elec_prop + mqc_statistic'// &
         ' through the real registered row: no crash, the state bitwise unchanged,'// &
         ' the holder intact across elec_prop)'
   end if
   ! (i) cycle 1: the dirty holder went through both dispatches - the snapshot
   ! must read back the defined empty state
   nfail0 = nfail                          ! the (i) water level
   call ws_get(w)
   if (w%repr /= r_empty .or. w%v_add /= 0.0d0 .or. allocated(w%i_want) .or. &
       allocated(w%g_add) .or. allocated(w%h_ij) .or. allocated(w%dh_ij) .or. &
       allocated(w%v_surf) .or. allocated(w%g_surf) .or. allocated(w%d_ij)) then
      call fail('(i) the mqc_statistic tail did not return the holder to the defined'// &
                ' empty state (repr / scalar / a per-mode array survived)')
   end if
   ! cycle 2: a fresh fill + dispatch pair on a different geometry - the holder
   ! dies again at the tail while the statistics channels (read after the cycle)
   ! keep their run-lifetime values
   fill_sel = 1
   sa%q = qe_ctr + ug
   sa%p = pg
   call force_eval(sa, ea)
   qref = sa%q
   pref = sa%p
   fref = sa%f
   call elec_prop(sa, dtv)
   call mqc_statistic(sa)
   call ws_get(w)
   if (w%repr /= r_empty .or. w%v_add /= 0.0d0 .or. allocated(w%i_want) .or. &
       allocated(w%g_add)) then
      call fail('(i) cycle 2: the tail reset did not fire on the refilled step')
   end if
   if (.not. all(sa%q == qref) .or. .not. all(sa%p == pref) .or. .not. all(sa%f == fref)) then
      call fail('(h) cycle 2: the dispatch pair changed the physical state')
   else if (nfail == nfail0) then
      print '(a)', '(i) step-end reset at the mqc_statistic tail: PASS (the dirty holder'// &
         ' reads back the defined empty state after the dispatch - all seven arrays'// &
         ' gone, scalar zeroed; clears again on the refilled step)'
   end if
   ! (j): the statistics write-back channels after two dispatch cycles - occ [1]
   ! (already read [1] pre-dispatch at criteria (d)/(e): constant across the
   ! dispatch boundary) and n_hop 0 (the accumulator survived two step-end
   ! resets untouched - the two holders have different lifetimes by statute)
   nfail0 = nfail                          ! the (j) water level
   call occ_get(occ_req)
   if (size(occ_req) /= 1 .or. occ_req(1) /= 1 .or. n_hop_get() /= 0) then
      call fail('(j) the statistics channels broke under dispatch (occ /= [1] or'// &
                ' n_hop /= 0 after the dispatch cycles)')
   else if (nfail == nfail0) then
      print '(a)', '(j) statistics write-back channels: PASS (occ = [1] constant across'// &
         ' the dispatch boundary, n_hop = 0 accumulated - untouched by the step-end'// &
         ' resets of (i))'
   end if

   ! ---------- (k)/(l): the loud-stop children ----------
   call expect_loudstop_at('container_init', 'mismatch', 'diab', 'scalar', &
      '(k) H-form mismatch loud-stop')
   call expect_loudstop_at('elec_prop', 'unbound-prop', 'method_bind', '', &
      '(l) unbound elec_prop loud-stop')
   call expect_loudstop_at('mqc_statistic', 'unbound-match', 'method_bind', '', &
      '(l) unbound mqc_statistic loud-stop')

   ! ---------- (m)-(p): the reference-container batch ----------
   call crit_m()      ! (m) the harm well through the container reg path (child)
   call crit_n()      ! (n) morse/leps module direct bind (children)
   call crit_o()      ! (o) box parameter file parsing (children, cwd = build/)
   call crit_p()      ! (p) leps diab full-matrix delivery + symmetry (module face)

   call ws_reset()                         ! leave the holder in the empty state
   call state_destroy(sa)

   print '(a,i0,a)', '=== check_facade: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! the local branch test container (bound once; the check-side fill_sel picks the
   ! delivery mode): 1 = the scalar harmonic channel, 2 = the two-surface linear
   ! diabatic channel with an additive scalar member, 3 = the one-surface diabatic
   ! twin of branch 1 (same V and gradient as the scalar branch, delivered as 1x1
   ! constructs; the additive channel stays absent). The gradient and energy in eV
   ! (the member unit contract). Branches 1 and 3 recompute u = q - qe_ctr and the
   ! V expression with IDENTICAL grouping, so the twin comparison is arithmetic-
   ! order clean
   subroutine check_pes(q, mass, ws)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused - pure PES branches)
      type(ws_t), intent(inout) :: ws  ! the holder (received through the facade's alias)
      real(8) :: u(6)
      u = q - qe_ctr
      seen_nwant = 0                    ! the request handshake log: what the container
      seen_want1 = 0                    ! observed in the pre-written request field
      if (allocated(ws%i_want)) then
         seen_nwant = size(ws%i_want)
         seen_want1 = ws%i_want(1)
      end if
      select case (fill_sel)
      case (1)                          ! scalar harmonic: V = 1/2 sum k u^2
         allocate (ws%g_add(size(q)))
         ws%g_add = keV*u                                 ! dV/dq [eV/A]
         ws%v_add = 0.5d0*dot_product(keV, u*u)           ! V [eV]
         ws%repr = r_scalar
      case (2)                          ! two-surface linear diabatic H + additive member
         allocate (ws%h_ij(2, 2), ws%dh_ij(size(q), 2, 2), ws%g_add(size(q)))
         ws%h_ij(1, 1) = hb1 + dot_product(a1, u)
         ws%h_ij(2, 2) = hb2 + dot_product(a2, u)
         ws%h_ij(1, 2) = hb12 + dot_product(cc, u)
         ws%h_ij(2, 1) = ws%h_ij(1, 2)
         ws%dh_ij(:, 1, 1) = a1                            ! dH_ij/du [eV/A] (symmetric)
         ws%dh_ij(:, 2, 2) = a2
         ws%dh_ij(:, 1, 2) = cc
         ws%dh_ij(:, 2, 1) = cc
         ws%g_add = g_uni                                  ! the additive member's kernel
         ws%v_add = v_uni                                  ! and energy [eV]
         ws%repr = r_diab
      case (3)                          ! one-surface diabatic twin of branch 1
         allocate (ws%h_ij(1, 1), ws%dh_ij(size(q), 1, 1))
         ws%dh_ij(:, 1, 1) = keV*u                         ! same expression as case (1)
         ws%h_ij(1, 1) = 0.5d0*dot_product(keV, u*u)       ! same expression as case (1)
         ws%repr = r_diab                                  ! (v_add stays at its zeroed default)
      end select
   end subroutine check_pes

   ! degenerate no-op method members for the (k) child only: the check-side
   ! 'chk_diab' package row demands an H form its container does not deliver -
   ! the members are never dispatched (the loud-stop fires at assembly), so
   ! no-op bodies honestly stand for "members exist" in the registration
   subroutine noop_prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine noop_prop

   subroutine noop_match(sta)
      type(state_t), intent(inout) :: sta
   end subroutine noop_match

   ! scan_exec_hits(path, nhit) - count this file's EXECUTABLE-position occurrences
   ! of the conversion literal: a hit inside comment text (a '!' at or before the
   ! hit position) does not count; a trailing comment after code counts (conservative
   ! - the single-site statute fails loudly rather than silently)
   subroutine scan_exec_hits(path, nhit)
      character(len=*), intent(in) :: path
      integer, intent(out) :: nhit
      integer :: uu, ios2, p, cq
      character(len=512) :: l2
      nhit = 0
      open (newunit=uu, file=path, status='old', action='read', iostat=ios2)
      if (ios2 /= 0) return
      do
         read (uu, '(a)', iostat=ios2) l2
         if (ios2 /= 0) exit
         p = index(l2, '23.0605')
         if (p > 0) then
            cq = index(l2, '!')
            if (cq == 0 .or. cq > p) nhit = nhit + 1
         end if
      end do
      close (uu)
   end subroutine scan_exec_hits

   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

   ! spawn the child for one loud-stop case and assert the stop: nonzero exit,
   ! stderr names the authority (the stopping entry) and carries the expected
   ! key / extra needle (an empty needle is skipped) - the sibling checks' scheme.
   ! The optional in_build spawns the child with cwd = build/ (run-directory
   ! material); the redirect target lands inside build/ either way
   subroutine expect_loudstop_at(who, cmode, key, needle, label, in_build)
      character(len=*), intent(in) :: who, cmode, key, needle, label
      logical, intent(in), optional :: in_build
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      logical :: ib
      integer :: stat, v, ios
      buf = ''
      ib = .false.
      if (present(in_build)) ib = in_build
      ! the shell redirect needs native separators on Windows (forward slashes make
      ! cmd misparse the redirect target)
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      if (ib) then
         ! cwd = build/ - the exe sits in that same directory; the EXPLICIT
         ! .-relative path spares the shell's executable search (cmd does not
         ! reliably search the cwd; '.' is never on a POSIX PATH)
         cmd = 'cd build && .'//sep//'check_facade.exe '//trim(cmode)//' 2>facade_err.txt'
      else
         cmd = 'build'//sep//'check_facade.exe '//trim(cmode)//' 2>build'//sep//'facade_err.txt'
      end if
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/facade_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) line
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(line)
         end do
         close (v)
      end if
      if (stat == 0) then
         print '(a)', 'FAIL '//trim(label)//': expected a loud-stop, the child exited 0'
         nfail = nfail + 1
      else if (index(buf, trim(who)) == 0) then
         print '(a)', 'FAIL '//trim(label)//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (len_trim(key) > 0 .and. index(buf, trim(key)) == 0) then
         print '(a)', 'FAIL '//trim(label)//': stderr missing the key "'//trim(key)//'"'
         nfail = nfail + 1
      else if (len_trim(needle) > 0 .and. index(buf, trim(needle)) == 0) then
         print '(a)', 'FAIL '//trim(label)//': stderr missing "'//trim(needle)//'"'
         nfail = nfail + 1
      else
         print '(a)', trim(label)//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop_at

   !------------------------------------------------------------------
   ! (m)/(n)/(o)/(p): the reference-container batch (criteria above)
   !------------------------------------------------------------------

   ! crit_m() - (m): the harm well through the container's own drag-in row
   subroutine crit_m()
      call expect_ok_child('mbind', 'mbind: OK', &
         '(m) harm through the container reg path (registry_gen drag-in)')
   end subroutine crit_m

   ! crit_n() - (n): the morse / leps module potentials direct-bound
   subroutine crit_n()
      call expect_ok_child('nmorse', 'nmorse: OK', '(n) morse module direct bind')
      call expect_ok_child('nleps', 'nleps: OK', '(n) leps module direct bind')
   end subroutine crit_n

   ! crit_o() - (o): the container parameter file (children with cwd = build/):
   !             the absent-file default posture, the good parse, the four
   !             named bad channels
   subroutine crit_o()
      character(len=80) :: lgood(3), lbad1(2), lbad2(2), lbad3(1), lbad4(1)
      ! absent file: default constant, both optional slots unbound
      call remove_file('build/reference_pars.txt')
      call expect_ok_child_in('obox', 'obox: k='//trim(d2s(4.0d-2))//' term=F class=F', &
         '(o) box parameter file absent: defaults, slots unbound')
      ! good file: the constant staged + both optional slots provided
      lgood = [ character(len=80) :: 'K_HARM = 0.123', 'Q_TERM = 9.5', 'E_SPLIT = 1.5' ]
      call write_pars('build/reference_pars.txt', lgood)
      call expect_ok_child_in('obox', 'obox: k='//trim(d2s(0.123d0))//' term=T class=T', &
         '(o) box parameter file good parse')
      ! bad files: every defect a named loud-stop
      lbad1 = [ character(len=80) :: 'K_HARM = 0.1', 'FOO = 3.0' ]
      call write_pars('build/reference_pars.txt', lbad1)
      call expect_loudstop_at('reference_load', 'obox', 'reference_pars.txt', 'FOO', &
         '(o) unknown-keyword loud-stop', in_build=.true.)
      lbad2 = [ character(len=80) :: 'K_HARM = 0.1', 'K_HARM = 0.2' ]
      call write_pars('build/reference_pars.txt', lbad2)
      call expect_loudstop_at('reference_load', 'obox', 'reference_pars.txt', 'twice', &
         '(o) duplicate-keyword loud-stop', in_build=.true.)
      lbad3 = [ character(len=80) :: 'K_HARM = abc' ]
      call write_pars('build/reference_pars.txt', lbad3)
      call expect_loudstop_at('reference_load', 'obox', 'reference_pars.txt', 'unreadable', &
         '(o) unreadable-value loud-stop', in_build=.true.)
      lbad4 = [ character(len=80) :: 'K_HARM = -0.5' ]
      call write_pars('build/reference_pars.txt', lbad4)
      call expect_loudstop_at('reference_load', 'obox', 'reference_pars.txt', 'positive', &
         '(o) non-positive-constant loud-stop', in_build=.true.)
      call remove_file('build/reference_pars.txt')   ! no stale file for later children
   end subroutine crit_o

   ! crit_p() - (p): the leps diabatic full-matrix delivery at the module face
   subroutine crit_p()
      real(8), parameter :: hcd = 1.0d-6   ! the central-difference step [A]
      real(8), parameter :: q3(9) = [ 0.10d0,-0.20d0, 0.15d0, 1.25d0, 0.30d0,-0.10d0, &
                                      0.05d0, 1.90d0, 0.40d0 ]  ! the (p) geometry [A]
      real(8) :: qh(9), ql(9)              ! the central-difference probes
      real(8) :: vsc, gsc(9)               ! the scalar reference (energy / gradient)
      real(8) :: hh(2, 2), dhd(9, 2, 2)    ! the diabatic delivery
      real(8) :: hhp(2, 2), hhm(2, 2)      ! the probe deliveries
      real(8) :: dhp(9, 2, 2), dhm(9, 2, 2)  ! (their gradient companions, scratch)
      real(8) :: hev(2, 2), evv(2)         ! the eig_sym copy (vectors by column)
      real(8) :: gfold(9), cd11, cd12
      integer :: i9
      nfail0 = nfail                       ! the (p) water level
      call leps_vg(3, q3, vsc, gsc)
      call leps_diab(3, q3, hh, dhd)
      ! symmetry: H symmetric and EVERY dh face symmetric, bitwise - the
      ! full-matrix delivery contract (both triangles)
      if (hh(1, 2) /= hh(2, 1) .or. any(dhd(:, 2, 1) /= dhd(:, 1, 2)) .or. &
          any(dhd(:, 2, 2) /= dhd(:, 1, 1))) then
         call fail('(p) delivery not symmetric (H or a dh face carries'// &
                   ' mismatched triangles)')
      end if
      ! the analytic faces vs the central difference of the delivered H entries
      do i9 = 1, 9
         qh = q3
         qh(i9) = q3(i9) + hcd
         ql = q3
         ql(i9) = q3(i9) - hcd
         call leps_diab(3, qh, hhp, dhp)
         call leps_diab(3, ql, hhm, dhm)
         cd11 = (hhp(1, 1) - hhm(1, 1))/(2.0d0*hcd)
         if (abs(dhd(i9, 1, 1) - cd11) > 1.0d-8*max(abs(dhd(i9, 1, 1)), 1.0d0)) then
            call fail('(p) dh face (1,1) vs central difference, dof '//trim(i2s(i9)))
         end if
         cd12 = (hhp(1, 2) - hhm(1, 2))/(2.0d0*hcd)
         if (abs(dhd(i9, 1, 2) - cd12) > 1.0d-8*max(abs(dhd(i9, 1, 2)), 1.0d0)) then
            call fail('(p) dh face (1,2) vs central difference, dof '//trim(i2s(i9)))
         end if
      end do
      ! the eigen identity: the delivered H diagonalizes to the scalar surface
      hev = hh
      call eig_sym(hev, evv)
      if (abs(evv(1) - vsc) > 1.0d-12*max(abs(vsc), 1.0d0)) then
         call fail('(p) eigen identity: lambda_1 /= the scalar leps energy ('// &
                   trim(d2s(evv(1)))//' vs '//trim(d2s(vsc))//')')
      end if
      ! the eigenvector fold through BOTH delivered triangles = the scalar gradient
      do i9 = 1, 9
         gfold(i9) = hev(1, 1)*hev(1, 1)*dhd(i9, 1, 1) &
                   + hev(1, 1)*hev(2, 1)*dhd(i9, 1, 2) &
                   + hev(2, 1)*hev(1, 1)*dhd(i9, 2, 1) &
                   + hev(2, 1)*hev(2, 1)*dhd(i9, 2, 2)
      end do
      if (any(abs(gfold - gsc) > 1.0d-12*max(1.0d0, abs(gsc)))) then
         call fail('(p) eigenvector fold /= the scalar leps gradient (max dev '// &
                   trim(d2s(maxval(abs(gfold - gsc))))//')')
      else if (nfail == nfail0) then
         print '(a)', '(p) leps diab full-matrix delivery: PASS (H and every dh face'// &
            ' symmetric bitwise; faces (1,1)/(1,2) = central difference to 1e-8;'// &
            ' lambda_1 = the scalar energy and the two-triangle eigenvector fold ='// &
            ' the scalar gradient, 1e-12)'
      end if
   end subroutine crit_p

   !------------------------------------------------------------------
   ! the batch-4 children (self-respawn modes; each terminates itself)
   !------------------------------------------------------------------

   ! child_mbind() - the (m) child: bind the reference container through the
   !                 build-generated registry (the driver's own registration
   !                 path), assemble, and force_eval the default-constant well
   subroutine child_mbind()
      real(8), parameter :: kk = 0.04d0    ! the container's default well constant [eV/A^2]
      type(state_t) :: sb
      real(8) :: cal, ea, eb, tk
      integer :: i
      ! the drag-in path exactly as the aggregated registry runs it: member
      ! init + the PES bind, the assembly seam (guard-skipped - no input
      ! staged), assembly, then a real force_eval
      call registry_gen_all()
      call registry_gen_init()
      call container_init()
      allocate (roster%mass(2))
      roster%mass = [ 1.0d0, 4.0d0 ]
      call state_create(sb, 2)
      ! calibration: one call at a known offset (p = 0) measures the composed
      ! conversion constant (the literal lives only in force_facade)
      sb%q = ueA
      sb%p = 0.0d0
      call force_eval(sb, ea)
      cal = -sb%f(1)/(kk*ueA(1))
      if (cal <= 0.0d0) then
         print '(a)', 'FAIL (m) child: measured conversion constant non-positive'
         stop 1
      end if
      do i = 1, 6
         if (abs(sb%f(i) + kk*sb%q(i)*cal) > 1.0d-13*max(abs(kk*sb%q(i)*cal), 1.0d0)) then
            print '(a)', 'FAIL (m) child: force mismatch at the calibration geometry,'// &
               ' dof '//trim(i2s(i))
            stop 1
         end if
      end do
      ! second geometry with momenta: force + the e_tot decomposition
      sb%q = ueB
      sb%p = pB
      call force_eval(sb, eb)
      do i = 1, 6
         if (abs(sb%f(i) + kk*sb%q(i)*cal) > 1.0d-13*max(abs(kk*sb%q(i)*cal), 1.0d0)) then
            print '(a)', 'FAIL (m) child: force mismatch, geometry B, dof '//trim(i2s(i))
            stop 1
         end if
      end do
      tk = 0.0d0
      do i = 1, 2
         tk = tk + dot_product(sb%p(3*i - 2:3*i), sb%p(3*i - 2:3*i))/(2.0d0*roster%mass(i))
      end do
      tk = tk/e_conv
      if (abs(eb - (tk + 0.5d0*kk*dot_product(sb%q, sb%q)*cal/e_conv)) > &
          1.0d-12*max(abs(eb), 1.0d0)) then
         print '(a)', 'FAIL (m) child: e_tot decomposition mismatch: '//trim(d2s(eb))
         stop 1
      end if
      print '(a)', 'mbind: OK'
      stop 0
   end subroutine child_mbind

   ! child_npot(use_morse, tag) - the (n) child: one module potential
   !                 direct-bound through a transparent adapter, assembled,
   !                 and force_eval'd - analytic gradient vs central
   !                 difference, facade fold, e_tot decomposition
   subroutine child_npot(use_morse, tag)
      logical, intent(in) :: use_morse
      character(len=*), intent(in) :: tag
      real(8), parameter :: hcd = 1.0d-6   ! the central-difference step [A]
      real(8), parameter :: qmo(6) = [ 0.0d0, 0.0d0, 0.0d0, 1.32d0, 0.21d0,-0.17d0 ] ! morse geometry [A]
      real(8), parameter :: qle(9) = [ 0.10d0,-0.20d0, 0.15d0, 1.25d0, 0.30d0,-0.10d0, &
                                       0.05d0, 1.90d0, 0.40d0 ]                    ! leps geometry [A]
      real(8), parameter :: pmo(6) = [ 0.5d0,-0.3d0, 0.4d0,-0.2d0, 0.6d0,-0.1d0 ]  ! momenta
      real(8), parameter :: ple(9) = [ 0.3d0,-0.2d0, 0.1d0, 0.4d0,-0.3d0, 0.2d0, &
                                       -0.1d0, 0.5d0,-0.4d0 ]
      type(state_t) :: sb
      real(8) :: qa(9), qh(9), ql(9), ga(9), gh(9), gl(9), gcdv(9)
      real(8) :: va, vh, vl, cal, ea, eb, tk
      integer :: nm, nd, i
      qa = 0.0d0
      ga = 0.0d0
      if (use_morse) then
         nm = 2
         qa(1:6) = qmo
         call container_bind_pes(morse_wrap)
      else
         nm = 3
         qa(1:9) = qle
         call container_bind_pes(leps_wrap)
      end if
      nd = 3*nm
      ! the degenerate 'adiabatic' method row must sit in the table before
      ! assembly (container_init binds it) - the direct-bind child registers
      ! the real package exactly as the aggregated registry would
      call reg_adiabatic_all()
      call container_init()
      allocate (roster%mass(nm))
      roster%mass = 1.0d0
      call state_create(sb, nm)
      ! the analytic reference from the module face
      if (use_morse) then
         call morse_vg(nm, qa, va, ga)
      else
         call leps_vg(nm, qa, va, ga)
      end if
      ! the facade fold: the measured conversion carries the module gradient
      sb%q(1:nd) = qa(1:nd)
      sb%p = 0.0d0
      call force_eval(sb, ea)
      cal = -dot_product(sb%f(1:nd), ga(1:nd))/dot_product(ga(1:nd), ga(1:nd))
      if (cal <= 0.0d0) then
         print '(a)', 'FAIL '//trim(tag)//' child: measured conversion constant non-positive'
         stop 1
      end if
      do i = 1, nd
         if (abs(sb%f(i) + ga(i)*cal) > 1.0d-12*max(abs(ga(i)*cal), 1.0d0)) then
            print '(a)', 'FAIL '//trim(tag)//' child: facade-fold force mismatch, dof '// &
               trim(i2s(i))
            stop 1
         end if
      end do
      ! the analytic gradient vs the central difference of the module energy
      do i = 1, nd
         qh = 0.0d0
         ql = 0.0d0
         qh(1:nd) = qa(1:nd)
         ql(1:nd) = qa(1:nd)
         qh(i) = qa(i) + hcd
         ql(i) = qa(i) - hcd
         if (use_morse) then
            call morse_vg(nm, qh, vh, gh)
            call morse_vg(nm, ql, vl, gl)
         else
            call leps_vg(nm, qh, vh, gh)
            call leps_vg(nm, ql, vl, gl)
         end if
         gcdv(i) = (vh - vl)/(2.0d0*hcd)
         if (abs(ga(i) - gcdv(i)) > 1.0d-8*max(abs(ga(i)), 1.0d0)) then
            print '(a)', 'FAIL '//trim(tag)//' child: analytic gradient vs central'// &
               ' difference, dof '//trim(i2s(i))//' ('//trim(d2s(ga(i)))//' vs '// &
               trim(d2s(gcdv(i)))//')'
            stop 1
         end if
      end do
      ! the e_tot decomposition at the same geometry with momenta
      if (use_morse) then
         sb%p(1:6) = pmo
      else
         sb%p(1:9) = ple
      end if
      call force_eval(sb, eb)
      tk = 0.0d0
      do i = 1, nm
         tk = tk + dot_product(sb%p(3*i - 2:3*i), sb%p(3*i - 2:3*i))/(2.0d0*roster%mass(i))
      end do
      tk = tk/e_conv
      if (abs(eb - (tk + va*cal/e_conv)) > 1.0d-12*max(abs(eb), 1.0d0)) then
         print '(a)', 'FAIL '//trim(tag)//' child: e_tot decomposition mismatch: '// &
            trim(d2s(eb))
         stop 1
      end if
      print '(a)', trim(tag)//': OK'
      stop 0
   end subroutine child_npot

   ! child_obox() - the (o) child: load the container parameter file from the
   !                run directory (cwd = build/) and print the staged posture
   subroutine child_obox()
      call reference_load()
      print '(a)', 'obox: k='//trim(d2s(reference_k()))//' term='// &
         merge('T', 'F', reference_term_on())//' class='// &
         merge('T', 'F', reference_class_on())
      stop 0
   end subroutine child_obox

   ! the transparent (n) adapters: the module potentials (v, g signature)
   ! delivered through the class-A container face (q, mass, ws) - the module
   ! material itself, bound check-side
   subroutine morse_wrap(q, mass, ws)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      type(ws_t), intent(inout) :: ws  ! workspace (filled: scalar mode)
      real(8) :: v                     ! the module energy return [eV]
      ws%repr = r_scalar
      if (allocated(ws%g_add)) deallocate (ws%g_add)
      allocate (ws%g_add(size(q)))
      call morse_vg(size(q)/3, q, v, ws%g_add)
      ws%v_add = v
   end subroutine morse_wrap

   subroutine leps_wrap(q, mass, ws)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      type(ws_t), intent(inout) :: ws  ! workspace (filled: scalar mode)
      real(8) :: v                     ! the module energy return [eV]
      ws%repr = r_scalar
      if (allocated(ws%g_add)) deallocate (ws%g_add)
      allocate (ws%g_add(size(q)))
      call leps_vg(size(q)/3, q, v, ws%g_add)
      ws%v_add = v
   end subroutine leps_wrap

   ! spawn the child for one success case (a fresh process, so its container
   ! bind is its own) and assert exit 0 + the expected stdout marker; the
   ! child's output (stdout, stderr merged) rides into any FAIL message
   subroutine expect_ok_child(cmode, marker, label)
      character(len=*), intent(in) :: cmode, marker, label
      character(len=512) :: cmd
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_facade.exe '//trim(cmode)//' >build'//sep// &
            'facade_ok.txt 2>&1'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp_local('build/facade_ok.txt', buf)
      if (stat /= 0) then
         print '(a)', 'FAIL '//trim(label)//': the child exited nonzero - '//trim(buf)
         nfail = nfail + 1
      else if (index(buf, trim(marker)) == 0) then
         print '(a)', 'FAIL '//trim(label)//': stdout missing "'//trim(marker)// &
            '" - '//trim(buf)
         nfail = nfail + 1
      else
         print '(a)', trim(label)//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_ok_child

   ! the cwd = build/ variant (the (o) children: the parameter file is
   ! run-directory material)
   subroutine expect_ok_child_in(cmode, marker, label)
      character(len=*), intent(in) :: cmode, marker, label
      character(len=512) :: cmd
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      ! cwd = build/ - the exe sits in that same directory; the EXPLICIT
      ! .-relative path spares the shell's executable search (cmd does not
      ! reliably search the cwd; '.' is never on a POSIX PATH)
      cmd = 'cd build && .'//sep//'check_facade.exe '//trim(cmode)//' >facade_ok.txt 2>&1'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp_local('build/facade_ok.txt', buf)
      if (stat /= 0) then
         print '(a)', 'FAIL '//trim(label)//': the child exited nonzero - '//trim(buf)
         nfail = nfail + 1
      else if (index(buf, trim(marker)) == 0) then
         print '(a)', 'FAIL '//trim(label)//': stdout missing "'//trim(marker)// &
            '" - '//trim(buf)
         nfail = nfail + 1
      else
         print '(a)', trim(label)//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_ok_child_in

   ! write one container parameter file (one KEYWORD = VALUE per record)
   subroutine write_pars(fname, lines)
      character(len=*), intent(in) :: fname, lines(:)
      integer :: v, j
      open (newunit=v, file=fname, status='replace', action='write')
      do j = 1, size(lines)
         write (v, '(a)') trim(lines(j))
      end do
      close (v)
   end subroutine write_pars

   ! remove one file (host-shaped; the Windows shell wants native separators)
   subroutine remove_file(path)
      character(len=*), intent(in) :: path
      character(len=512) :: pnat
      integer :: j
      if (windows_host()) then
         pnat = ''
         do j = 1, len_trim(path)
            if (path(j:j) == '/') then
               pnat(j:j) = '\'
            else
               pnat(j:j) = path(j:j)
            end if
         end do
         call execute_command_line('del "'//trim(pnat)//'" >nul 2>nul', wait=.true.)
      else
         call execute_command_line('rm -f "'//trim(path)//'"', wait=.true.)
      end if
   end subroutine remove_file

   ! slurp one text file into a single-line buffer (bound-guarded; a missing
   ! file leaves the buffer empty)
   subroutine slurp_local(fname, buf)
      character(len=*), intent(in) :: fname
      character(len=*), intent(inout) :: buf
      character(len=512) :: line
      integer :: v, ios
      buf = ''
      open (newunit=v, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(buf) + len_trim(line) + 1 < len(buf)) then
            buf = trim(buf)//' '//trim(line)
         end if
      end do
      close (v)
   end subroutine slurp_local

   ! fail(msg) - record one criterion failure
   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      print '(a)', 'FAIL '//trim(msg)
      nfail = nfail + 1
   end subroutine fail

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

end program check_facade
