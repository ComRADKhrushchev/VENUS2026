!=====================================================================
! check_facade.f90 - interface-mechanism criteria program (E27 contract):
!                    the statistics-holder seeding (d) + the force_eval
!                    delivery chain (e)-(f) + the electronic-phase
!                    dispatch batch (h)/(j)/(l) + the reference-container
!                    batch (m)-(p) + the acceptance batch (q)/(s):
!                    gradient self-consistency sweep and the old-vs-new
!                    bit-level force comparison (criteria defined FIRST
!                    in this header; this program executes them)
! E27 retirement note (2026-09-16 user ruling): the former (a)/(b)/(c)
!   ws-holder legs, the (g) interface eigen-fold, the (i) step-end ws
!   reset, the (k) H-form matchmaking named abort and the (r) interface-side
!   degenerate equivalence RETIRED with the ws interchange (elec_ws.f90
!   deleted; force_eval narrowed to the composed-kernel delivery; the
!   fold physics survives at the module face in (p)). The surviving
!   legs keep their historical letters for archival traceability.
! Compile (from src/; run make FIRST - the aggregated registry module is
! a build product; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree). STEP 1 builds the archived legacy PES
! chain read-only (criterion (s) links its two evaluation entries into
! this same executable - one compiler, one link, one double arithmetic;
! the objects and .mod files land under build/old_xcheck/, the archived
! tree itself is never written):
!   gfortran -c -Jbuild/old_xcheck -o build/old_xcheck/venus_params.o ../venus2026/venus_params.f90
!   gfortran -c -Ibuild/old_xcheck -Jbuild/old_xcheck -o build/old_xcheck/venus_data.o ../venus2026/venus_data.f90
!   gfortran -c -Jbuild/old_xcheck -o build/old_xcheck/input_parser.o ../venus2026/input_parser.f90
!   gfortran -c -Jbuild/old_xcheck -o build/old_xcheck/emt_nn_pes.o ../venus2026/src_TEST/emt_nn_pes.f90
!   gfortran -c -Ibuild/old_xcheck -Jbuild/old_xcheck -o build/old_xcheck/rst_pes.o ../venus2026/src_TEST/rst_pes.f90
!   gfortran -c -Ibuild/old_xcheck -Jbuild/old_xcheck -o build/old_xcheck/rst_slab_bvk.o ../venus2026/src_TEST/rst_slab_bvk.f90
!   gfortran -c -Ibuild/old_xcheck -Jbuild/old_xcheck -o build/old_xcheck/test_potentials.o ../venus2026/src_TEST/test_potentials.f90
!   gfortran -c -Ibuild/old_xcheck -Jbuild/old_xcheck -o build/old_xcheck/interface_TEST.o ../venus2026/src_TEST/interface_TEST.f90
! STEP 2 compiles the check (same flags as every sibling check - the
! legacy objects join the link as the final inputs):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 utils/rng.f90 &
!            globals/config_atoms.f90 globals/state.f90 &
!            interface/elec_interface.f90 methods/elec_adiabatic.f90 &
!            methods/reg_adiabatic.f90 utils/linalg.f90 utils/spectrum.f90 utils/hessian.f90 &
!            interface/force_interface.f90 interface/spectrum_interface.f90 &
!            interface/sampler.f90 interface/bath.f90 utils/geometry.f90 utils/geom_init.f90 dynamics/propagator.f90 &
!            methods/samp_boltzmann.f90 methods/samp_ebk.f90 methods/beam_laws.f90 methods/inc_surface.f90 methods/inc_pair.f90 methods/inc_single.f90 &
!            methods/samp_j.f90 methods/samp_rotation.f90 methods/samp_surface_oscillator.f90 &
!            methods/samp_normalmode.f90 methods/samp_barrier_excitation.f90 &
!            methods/bath_andersen.f90 methods/samp_glo_target.f90 &
!            methods/samp_thermalize.f90 methods/samp_wigner.f90 methods/seam_cast.f90 &
!            methods/reg_prop.f90 &
!            methods/reg_samp_boltzmann.f90 methods/reg_samp_ebk.f90 &
!            methods/reg_inc_surface.f90 methods/reg_inc_pair.f90 methods/reg_inc_single.f90 methods/reg_samp_j.f90 &
!            methods/reg_samp_rotation.f90 methods/reg_samp_surface_oscillator.f90 &
!            methods/reg_samp_normalmode.f90 methods/reg_samp_barrier_excitation.f90 &
!            methods/reg_bath_andersen.f90 methods/reg_samp_glo_target.f90 &
!            methods/reg_samp_thermalize.f90 methods/reg_samp_wigner.f90 &
!            output/final_state.f90 &
!            ../systems/reference/reference_pot.f90 &
!            ../systems/reference/reference_box.f90 &
!            interface/container_sched.f90 &
!            build/gfortran/registry_gen.f90 ../scripts/check_facade.f90 &
!            -Ibuild/old_xcheck build/old_xcheck/venus_params.o &
!            build/old_xcheck/venus_data.o build/old_xcheck/input_parser.o &
!            build/old_xcheck/emt_nn_pes.o build/old_xcheck/rst_pes.o &
!            build/old_xcheck/rst_slab_bvk.o build/old_xcheck/test_potentials.o &
!            build/old_xcheck/interface_TEST.o -o build/check_facade.exe
! All artifacts are written under build/ - the source tree stays clean.
! Posture: criterion (d) exercises the statistics-holder seeding - the
!   REAL 'adiabatic' package registers through its reg module
!   (reg_adiabatic_all, the same row the aggregated registry executes)
!   and binds the name, which seeds the holder and makes the occ_get /
!   n_hop_get channels live; the members are not dispatched by that
!   batch. Criterion (e) additionally
!   assembles a FORCE CHAIN check-side: a local analytic container
!   is bound through container_bind_pes (the public assembly face), so
!   every (e) assertion runs on real force_eval calls; the default
!   binding posture of the tree is untouched (the driver binds nothing).
!   The container delivers the composed scalar kernel (the E27
!   delivery contract: g and v fully assigned in eV) - a local
!   procedure of this check; none of it ever enters src/. Criterion
!   (h)/(j) dispatch the REAL registered 'adiabatic' row through the two
!   electronic-phase entries on that same bound force chain (real calls,
!   not stub invocations). Named-abort verification approach (same scheme
!   as the sibling checks): criterion (l) runs as CHILD PROCESSES - the
!   check re-spawns itself with a mode argument (execute_command_line,
!   stderr redirected to build/); the child drives the real stop path and
!   must print a naming message on stderr and exit nonzero. The batch-4
!   real-container criteria ((m)/(n)) run as children for a second
!   reason: the parent's single PES slot is already bound check-side for
!   (e)-(j) (one container per process - the bind face aborts with a named error on a
!   second bind), so every different container gets a fresh process; the
!   (o) children run with cwd = build/ because their obox_input.txt case
!   input is run-directory material (the container keys ride the
!   program input). The acceptance batch rides the same
!   respawn scheme: its (q)/(s) children each bind ONE potential per
!   fresh process (the single force slot), and its (s) children
!   additionally drive the legacy engine's PES entries in the SAME
!   process - the bit-level comparison is in-process by design (one
!   executable, one floating-point arithmetic; a cross-executable
!   comparison via printed tables could not carry the last bits).
!   Red-phase lesson (the fake-green first run): the children respawn
!   THIS executable at its own path - a broken red-phase copy must
!   occupy the spawn path itself (same -o target), or the child
!   criteria run the intact binary and pass green for nothing.
! Criteria:
!   (d) statistics-holder seeding: after the real 'adiabatic' package
!       registered and the name bound, occ_get hands back the seeded
!       degenerate occupation [1] and n_hop_get 0 (the run-lifetime
!       accumulator at its initial value) - the holders are live before
!       any dispatch, and (j) re-asserts them across the dispatch
!       boundary. [The former (a)/(b)/(c) ws-holder legs and the
!       ws-assertion halves of the old (d) RETIRED AT E27 - the holder
!       itself no longer exists]
!   (e) end-to-end delivery chain: with the analytic container bound
!       through the public face, force_eval returns the analytic force
!       component-wise at two geometries (proportionality to the
!       gradient - the conversion constant cancels) and e_tot decomposes
!       as t_kin(p) + v_pot(q) via the composed conversion measured by
!       one calibration call (the eV->internal literal lives only in
!       force_interface - this check NEVER writes it). [The i_want request
!       handshake and the holder read-back halves RETIRED AT E27]
!   (f) conversion single-site: a recursive scan of every src/ .f90 file
!       (findstr on a Windows host, grep elsewhere; paths under build/ are
!       build artifacts - stale check copies, generated registries - and
!       sit OUTSIDE the scan domain) lists the files whose
!       text contains the conversion literal; each listed file is then
!       scanned line by line, and exactly ONE occurrence of the literal
!       sits at EXECUTABLE position (not inside a comment) in the whole
!       tree - the ev_kcal parameter declaration of force_interface.f90
!       (the unit-contract grep requirement, executable face)
!   (h) electronic-phase dispatch reachable + no-op-safe: with the real
!       'adiabatic' package bound and the force chain live, the pair
!       elec_prop(sta, dt) / mqc_statistic(sta) REALLY DISPATCHES
!       through the table without crashing and leaves the physical state
!       bitwise unchanged (q/p/f/t, n_surf, the ontic containers a/rho
!       staying unallocated). [The held-ws-intact assertions RETIRED AT
!       E27; the step-end reset subject also RETIRED (former (i))]
!   (j) statistics write-back channels: across the dispatch cycles
!       occ_get hands back [1] CONSTANT (the degenerate occupation the
!       match member leaves untouched) and n_hop_get 0 (the run-lifetime
!       accumulator holding its degenerate value across steps)
!   (l) unbound-method posture (child): elec_prop / mqc_statistic called
!       before any method_bind named abort, stderr naming the entry and
!       method_bind (the skeleton's legislation: unbound = named abort,
!       never a silent skip). [The former (k) H-form mismatch child
!       RETIRED AT E27 - the matchmaking died with the ws interchange]
!   (m) container registry path (child): the harmonic well bound through
!       its OWN drag-in row - the build-generated registry_gen module,
!       the exact registration call the driver runs - NOT the check-side
!       direct bind of (e). The child runs registry_gen_all (which
!       registers every method package AND binds the reference PES),
!       registry_gen_init (the assembly
!       seam, guard-skipped without buffered input), container_init
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
!       gradient criterion family), the interface fold reproducing the
!       module gradient through the measured conversion (1e-12), and
!       the e_tot decomposition at the same geometry with momenta
!   (o) box container keys (children, cwd = build/): the container's
!       parameter channel = the declared input keys of obox_input.txt
!       (input_declare_keys -> read_input -> pull -> reference_load):
!       ABSENT container keys leave the default constant and both
!       optional slots unbound (k = default, term = F, class = F); the
!       GOOD keys (K_HARM/Q_TERM/E_SPLIT) stage the constant and provide
!       both optional slots (k, term = T, class = T - asserted on the
!       child's printed marker, after a clean input_audit_declared);
!       the BAD inputs (unknown key, duplicate key, unreadable value,
!       non-positive constant) named abort at their owning site
!       (read_input / pull / reference_load) naming the key
!   (p) the leps diabatic full-matrix delivery (module face, in
!       process): leps_diab at a fixed geometry delivers H symmetric and
!       EVERY dh face symmetric (bitwise, both triangles); the analytic
!       faces (1,1)/(1,2) match the central difference of the delivered
!       H entries (1e-8, h = 1e-6); the eigen identity - eig_sym of the
!       delivered H yields lambda_1 = the scalar leps_vg energy (1e-12) -
!       and the eigenvector-fold gradient sum_ij u_i u_j (dH/dq)_ij
!       computed through BOTH delivered triangles = the scalar leps_vg
!       gradient (1e-12): the n_surf = 1 diab-vs-scalar degenerate
!       equivalence, made executable at the module face. [This leg is
!       where the fold physics survives after the E27 retirement of the
!       interface-side (g)/(r)]
!   (q) gradient self-consistency, the full sweep (children, one
!       potential per process): EACH of the three reference potentials
!       x three constructed geometries (an asymmetric point, a
!       symmetric-degenerate point - the isosceles pair distances equal
!       BITWISE / the stationary pair at re / the near-origin well point
!       - and a stretched-or-wide point): the module analytic gradient
!       vs the CENTRAL DIFFERENCE of the module energy (h = 1e-6,
!       tolerance 1e-8 relative with floor 1 eV/A, component-wise - the
!       (n) criterion family extended to harm and to multi-geometry),
!       then through a real force_eval the interface fold reproducing the
!       module gradient via the measured conversion (1e-12) at EVERY
!       geometry, and the e_tot = T + V decomposition with momenta at
!       the first geometry (1e-12). The harm child rides the CONTAINER
!       REGISTRY PATH (the registry_gen drag-in with the default
!       constant well K = 0.04 eV/A^2); the morse / leps children ride
!       the transparent direct-bind adapters (the (n) posture)
!   (s) the old-vs-new bit-level force comparison (children, one
!       potential per process, the legacy PES chain linked read-only):
!       the SAME Q fed element-wise into the legacy engine's two
!       evaluation entries - POT0 (returns Vpot = v_ev x 23.0605
!       [kcal/mol]) and DPESHON (accumulates PDOT = PDOT - g x 23.0605
!       x 0.04184 [a.u.], the accumulator zeroed before each call) -
!       and into force_eval (same potential, same parameters: the
!       module defaults 4.746 eV / 1.401 A / 1.028 1/A / 0.164 for
!       morse and leps on both sides; for harm the legacy well constant
!       buffered through ITS OWN keyword channel at 0.04 eV/A^2 against
!       the container's default well), at p = 0 (T = 0 exactly, so
!       e_tot carries the potential channel alone): assert
!       max|dF| = max|dE| = 0 BITWISE - same double arithmetic, the
!       character-identical analytic forms, the same left-to-right
!       conversion composition (g x 23.0605 x 0.04184). A nonzero
!       difference FAILS and is then located per the disposition plan
!       (a port error is fixed in the port; a benign arithmetic-order
!       difference is recorded with its value and cause and falls back
!       to 1e-15 with the archive note). The legacy side is driven
!       through its own call convention - the keyword file buffered by
!       the parent, load_input, POTPRE, then per-geometry POT0 /
!       DPESHON on the live coordinate array; the geometry set is the
!       (q) sweep (asymmetric / symmetric-degenerate / stretched-wide
!       per potential)
!=====================================================================
program check_facade
   use consts,        only: e_conv
   use state,         only: state_t, state_create, state_destroy
   use config_atoms, only: list_atoms
   use elec_interface,   only: occ_get, n_hop_get, method_bind, &
                            elec_prop, mqc_statistic
   use force_interface,  only: force_eval, container_bind_pes, container_init
   use reg_adiabatic, only: reg_adiabatic_all
   use linalg,        only: eig_sym
   use registry_gen,  only: registry_gen_all, registry_gen_init, registry_gen_keys
   use reference_pot, only: harm_vg, morse_vg, leps_vg, leps_diab
   use reference_box, only: reference_params_t, reference_keys, reference_load, &
                           reference_k, reference_term_on, &
                           reference_class_on
   use input,        only: read_input, input_declare_keys, input_audit_declared, pull
   implicit none
   integer, allocatable :: occ_req(:)      ! the held occupation set (occ_get)
   integer :: nfail = 0
   ! ---------- delivery-chain section state ((e)) ----------
   type(state_t) :: sa                     ! the force-chain state (2 atoms / 6 dofs)
   integer :: nfail0                       ! per-criterion failure watermark
   integer :: m                            ! fold loop index
   integer :: nexec, nhit                  ! (f) executable-position hit counters
   integer :: v, ios, jj                   ! (f) scan unit / status / scratch
   character(len=512) :: path              ! (f) scan path buffer
   character(len=512) :: exec_path         ! (f) the file that carried the last hit
   character(len=512) :: cmd, ln           ! (f) scan command / line buffer
   real(8) :: cal, ea, eb, tk              ! the measured conversion + energy channels
   ! harmonic test-potential parameters (per-dof, 6 dofs / 2 atoms of mass 1 and 4 amu)
   real(8), parameter :: keV(6)   = [ 0.95d0, 1.15d0, 0.85d0, 1.25d0, 0.75d0, 1.10d0 ] ! stiffness [eV/A^2]
   real(8), parameter :: qe_ctr(6) = [ 0.12d0,-0.18d0, 0.28d0,-0.14d0, 0.22d0,-0.06d0 ] ! centers [A]
   ! criterion-(e) geometries / momenta
   real(8), parameter :: ueA(6) = [ 0.30d0,-0.21d0, 0.18d0,-0.12d0, 0.24d0,-0.09d0 ] ! calibration offset [A]
   real(8), parameter :: ueB(6) = [-0.16d0, 0.26d0,-0.09d0, 0.31d0,-0.20d0, 0.14d0 ] ! second geometry [A]
   real(8), parameter :: pB(6)  = [ 1.05d0,-0.65d0, 0.48d0,-0.85d0, 0.55d0, 1.20d0 ] ! second-geometry momenta
   ! ---------- dispatch-batch locals ((h)-(l)) ----------
   character(len=32) :: mode               ! the child spawn's mode argument
   type(state_t) :: st1                    ! the child-mode minimal state
   real(8) :: qref(6), pref(6), fref(6)    ! the state-invariance references
   real(8) :: tref                         ! (ditto - the time member)
   real(8), parameter :: dtv = 0.05d0      ! the dt channel value [10 fs] fed to elec_prop
                                           ! (the fixed buffered step's stand-in; the degenerate
                                           ! member ignores it - the criterion needs the channel
                                           ! CARRIED, not consumed)

   ! ---------- acceptance-batch sweep tables ((q)/(r)/(s)) ----------
   ! harm geometries [A]: asymmetric A / asymmetric B / near-stationary
   real(8), parameter :: qh_list(6, 3) = reshape([ 0.30d0,-0.21d0, 0.18d0,-0.12d0, 0.24d0,-0.09d0, &
                                                   -0.16d0, 0.26d0,-0.09d0, 0.31d0,-0.20d0, 0.14d0, &
                                                    0.03d0,-0.02d0, 0.025d0,-0.015d0, 0.02d0,-0.01d0 ], [6, 3])
   ! morse geometries [A]: asymmetric / stretched / stationary at re (the pair on z)
   real(8), parameter :: qm_list(6, 3) = reshape([ 0.0d0, 0.0d0, 0.0d0, 1.32d0, 0.21d0,-0.17d0, &
                                                   0.0d0, 0.0d0, 0.0d0, 2.31d0,-0.18d0, 0.25d0, &
                                                   0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.401d0 ], [6, 3])
   ! leps geometries [A]: asymmetric / isosceles (r13 = r23 by the SAME expression - the two
   ! branch distances equal bitwise) / wide near-dissociated (the two sheets closest)
   real(8), parameter :: ql_list(9, 3) = reshape([ 0.10d0,-0.20d0, 0.15d0, 1.25d0, 0.30d0,-0.10d0, &
                                                   0.05d0, 1.90d0, 0.40d0, &
                                                   0.0d0, 0.0d0, 0.0d0, 1.10d0, 0.0d0, 0.0d0, &
                                                   0.55d0, 0.80d0, 0.0d0, &
                                                   0.0d0, 0.0d0, 0.0d0, 4.10d0, 0.30d0,-0.20d0, &
                                                  -1.60d0, 3.90d0, 0.50d0 ], [9, 3])
   real(8), parameter :: pm2(6) = [ 0.5d0,-0.3d0, 0.4d0,-0.2d0, 0.6d0,-0.1d0 ]  ! 2-atom momenta
   real(8), parameter :: pm3(9) = [ 0.3d0,-0.2d0, 0.1d0, 0.4d0,-0.3d0, 0.2d0, &
                                   -0.1d0, 0.5d0,-0.4d0 ]                        ! 3-atom momenta

   ! ---------- child dispatch: one named abort case per spawn (criteria (k)/(l)) ----------
   ! the child drives the real stop path; if the chain ever RETURNS, the child
   ! prints the bug witness and exits 0 - the parent's exitstat check fails
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      if (trim(mode) == 'unbound-prop') then
         ! (l): elec_prop before any method_bind - no row bound, named abort
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
         ! buffered there); good posture prints the marker, bad files named abort
         call child_obox()
      else if (mode(1:2) == 'q-') then
         ! (q): the gradient-sweep child (one potential per process)
         call child_q(pot_kind(mode(3:)))
      else if (mode(1:2) == 's-') then
         ! (s): the old-vs-new bit-comparison child (one potential per process;
         ! the legacy PES entries driven in the same process)
         call child_s(pot_kind(mode(3:)))
      end if
      print '(a)', 'child: NO LOUD-STOP (the child chain returned - bug)'
      stop 0
   end if

   ! the real 'adiabatic' package registers through its reg module (the same
   ! call the aggregated registry executes) + binds the name, which seeds the
   ! statistics holder - occ_get / n_hop_get are live for criterion (d); the
   ! members are dispatched for real at batch (h)/(j)
   call reg_adiabatic_all()
   call method_bind('adiabatic')

   ! ---------- (d) statistics-holder seeding ----------
   ! [the former (a)/(b)/(c) ws-holder legs RETIRED AT E27 - the holder no
   !  longer exists; this leg keeps the surviving statistics half]
   call occ_get(occ_req)
   if (size(occ_req) /= 1 .or. occ_req(1) /= 1 .or. n_hop_get() /= 0) then
      print '(a)', 'FAIL (d) statistics-holder seeding: occ /= [1] or n_hop /= 0'// &
                   ' at the bind point'
      nfail = nfail + 1
   else
      print '(a)', '(d) statistics-holder seeding: PASS (occ seeded [1], n_hop 0 -'// &
                   ' the holders live before any dispatch)'
   end if

   ! ---------- (e) end-to-end delivery chain ----------
   ! assembly check-side: the atom list masses, a 2-atom state, and the local
   ! analytic container bound through the PUBLIC assembly face (one slot,
   ! bound once)
   allocate (list_atoms%mass(2))
   list_atoms%mass = [ 1.0d0, 4.0d0 ]
   call state_create(sa, 2)
   call container_bind_pes(check_pes)
   ! calibration: one force_eval at a known single-dof offset (p = 0) measures the
   ! composed conversion constant cal (the eV->internal product the interface applies;
   ! the literal lives only in force_interface - never written here)
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
      tk = tk + dot_product(sa%p(3*m - 2:3*m), sa%p(3*m - 2:3*m))/(2.0d0*list_atoms%mass(m))
   end do
   tk = tk/e_conv
   if (abs(eb - (tk + 0.5d0*dot_product(keV, (sa%q - qe_ctr)*(sa%q - qe_ctr))*cal/e_conv)) > &
       1.0d-12*max(abs(eb), 1.0d0)) then
      call fail('(e) e_tot decomposition mismatch, geometry B: '//trim(d2s(eb)))
   end if
   ! the PASS line is water-level gated: it prints only when NO (e) assertion
   ! above failed (nfail == nfail0, the watermark taken before the first (e)
   ! pre-assertion)
   if (nfail == nfail0) then
      print '(a)', '(e) delivery chain: PASS (analytic force + e_tot at two'// &
         ' geometries through the composed-kernel contract, the measured'// &
         ' conversion carrying both)'
   end if

   ! ---------- (f) conversion single-site ----------
   ! list every src/ file whose text carries the conversion literal, then count the
   ! EXECUTABLE-position occurrences (comment text does not count); the unit contract
   ! demands exactly one - the interface's conversion-parameter declaration
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
         ! the scan domain is the SOURCE tree: paths under build/ (stale
         ! check copies, generated material) are artifacts and are skipped
         if (index(path, 'build/') > 0 .or. index(path, 'build\') > 0) cycle
         call scan_exec_hits(trim(path), nhit)
         if (nhit > 0) then
            nexec = nexec + nhit
            exec_path = trim(path)
         end if
      end do
      close (v)
      if (nexec /= 1 .or. index(exec_path, 'force_interface.f90') == 0) then
         call fail('(f) conversion single-site violated: '//trim(i2s(nexec))// &
            ' executable-position occurrence(s) in src/, last in "'//trim(exec_path)//'"'// &
            ' (want exactly 1, in force_interface.f90)')
      else
         print '(a)', '(f) conversion single-site: PASS (exactly one executable-position'// &
            ' occurrence of the conversion literal in src/: '//trim(exec_path)//')'
      end if
   end if

   ! ---------- (h)/(j): the electronic-phase dispatch batch ----------
   ! [the former (g) interface eigen-fold + the (i) step-end ws reset RETIRED AT
   !  E27; the fold physics survives at the module face in (p)]
   ! cycle 1: a real force_eval (the same chain the driver runs); the
   ! state-invariance references are taken after it
   sa%q = qe_ctr + ueB
   sa%p = pB
   call force_eval(sa, ea)
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
   call mqc_statistic(sa)                  ! REAL dispatch: the no-op match
   if (.not. all(sa%q == qref) .or. .not. all(sa%p == pref) .or. .not. all(sa%f == fref) .or. &
       sa%t /= tref .or. sa%s%n_surf /= 1 .or. allocated(sa%s%a) .or. allocated(sa%s%rho)) then
      call fail('(h) mqc_statistic changed the physical state under the adiabatic no-op law')
   else if (nfail == nfail0) then
      print '(a)', '(h) dispatch reachable + no-op-safe: PASS (elec_prop + mqc_statistic'// &
         ' through the real registered row: no crash, the state bitwise unchanged)'
   end if
   ! cycle 2: a fresh force + dispatch pair on a different geometry
   nfail0 = nfail
   sa%q = qe_ctr + ueA
   sa%p = pB
   call force_eval(sa, ea)
   qref = sa%q
   pref = sa%p
   fref = sa%f
   call elec_prop(sa, dtv)
   call mqc_statistic(sa)
   if (.not. all(sa%q == qref) .or. .not. all(sa%p == pref) .or. .not. all(sa%f == fref)) then
      call fail('(h) cycle 2: the dispatch pair changed the physical state')
   else if (nfail == nfail0) then
      print '(a)', '(h) cycle 2: PASS (a fresh force + dispatch pair on a different'// &
         ' geometry - the no-op law holds across cycles)'
   end if
   ! (j): the statistics write-back channels after two dispatch cycles - occ [1]
   ! (constant across the dispatch boundary) and n_hop 0 (the run-lifetime
   ! accumulator)
   nfail0 = nfail                          ! the (j) water level
   call occ_get(occ_req)
   if (size(occ_req) /= 1 .or. occ_req(1) /= 1 .or. n_hop_get() /= 0) then
      call fail('(j) the statistics channels broke under dispatch (occ /= [1] or'// &
                ' n_hop /= 0 after the dispatch cycles)')
   else if (nfail == nfail0) then
      print '(a)', '(j) statistics write-back channels: PASS (occ = [1] constant across'// &
         ' the dispatch boundary, n_hop = 0 accumulated)'
   end if

   ! ---------- (l): the named abort children ----------
   call expect_loudstop_at('elec_prop', 'unbound-prop', 'method_bind', '', &
      '(l) unbound elec_prop loud-stop')
   call expect_loudstop_at('mqc_statistic', 'unbound-match', 'method_bind', '', &
      '(l) unbound mqc_statistic loud-stop')

   ! ---------- (m)-(p): the reference-container batch ----------
   call crit_m()      ! (m) the harm well through the container reg path (child)
   call crit_n()      ! (n) morse/leps module direct bind (children)
   call crit_o()      ! (o) box parameter file parsing (children, cwd = build/)
   call crit_p()      ! (p) leps diab full-matrix delivery + symmetry (module face)

   ! ---------- (q)/(s): the acceptance batch ----------
   call crit_q()      ! (q) gradient self-consistency sweep (children)
   call crit_s()      ! (s) old-vs-new bit-level force comparison (children)

   call state_destroy(sa)

   print '(a,i0,a)', '=== check_facade: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! the local analytic test container (bound once): the scalar harmonic
   ! channel - V = 1/2 sum k (q - qe_ctr)^2, g = keV*(q - qe_ctr) [eV/A]
   ! (the E27 delivery contract: g and v fully assigned, in eV)
   subroutine check_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused - a pure PES)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A]
      real(8), intent(out) :: v        ! V [eV]
      g = keV*(q - qe_ctr)                                 ! dV/dq [eV/A]
      v = 0.5d0*dot_product(keV, (q - qe_ctr)*(q - qe_ctr)) ! V [eV]
   end subroutine check_pes

   ! scan_exec_hits(path, nhit) - count this file's EXECUTABLE-position occurrences
   ! of the conversion literal: a hit inside comment text (a '!' at or before the
   ! hit position) does not count; a trailing comment after code counts (conservative
   ! - the single-site rule fails loudly rather than silently)
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

   ! spawn the child for one named abort case and assert the stop: nonzero exit,
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
         cmd = 'cd build && .'//sep//'check_facade.exe '//trim(cmode)//' 2>interface_err.txt'
      else
         cmd = 'build'//sep//'check_facade.exe '//trim(cmode)//' 2>build'//sep//'interface_err.txt'
      end if
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/interface_err.txt', status='old', action='read', iostat=ios)
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

   ! crit_o() - (o): the container's declared input keys (children with
   !             cwd = build/, reading obox_input.txt): the absent-keys
   !             default posture, the good parse, the four named bad
   !             channels (each aborting at its owning site)
   subroutine crit_o()
      character(len=80) :: req(6)
      req = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 5.0', 'B_MAX = 2.0', &
              'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = BOLTZMANN' ]
      ! absent container keys: default constant; under Tier 1 presence-binding
      ! the term/classify wrappers are bound and the found flags gate them off
      call write_pars('build/obox_input.txt', req)
      call expect_ok_child_in('obox', 'obox: k='//trim(d2s(4.0d-2))//' term=F class=F', &
         '(o) container keys absent: defaults, found-flag gates off')
      ! good keys: the constant buffered + both optional slots provided
      call write_pars('build/obox_input.txt', [req, [character(len=80) :: &
         'K_HARM = 0.123', 'Q_TERM = 9.5', 'E_SPLIT = 1.5']])
      call expect_ok_child_in('obox', 'obox: k='//trim(d2s(0.123d0))//' term=T class=T', &
         '(o) container keys good parse')
      ! bad inputs: every defect a named abort at its owning site
      call write_pars('build/obox_input.txt', [req, [character(len=80) :: 'FOO = 3.0']])
      call expect_loudstop_at('read_input', 'obox', 'unknown key FOO', 'FOO', &
         '(o) unknown-key abort', in_build=.true.)
      call write_pars('build/obox_input.txt', [req, [character(len=80) :: &
         'K_HARM = 0.1', 'K_HARM = 0.2']])
      call expect_loudstop_at('read_input', 'obox', 'duplicate key', 'K_HARM', &
         '(o) duplicate-key abort', in_build=.true.)
      call write_pars('build/obox_input.txt', [req, [character(len=80) :: 'K_HARM = abc']])
      call expect_loudstop_at('pull', 'obox', 'bad real', 'K_HARM', &
         '(o) unreadable-value abort', in_build=.true.)
      call write_pars('build/obox_input.txt', [req, [character(len=80) :: 'K_HARM = -0.5']])
      call expect_loudstop_at('reference_load', 'obox', 'K_HARM must be positive', 'positive', &
         '(o) non-positive-constant abort', in_build=.true.)
      call remove_file('build/obox_input.txt')   ! no stale file for later children
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
   ! (q)/(r)/(s): the acceptance batch (criteria above)
   !------------------------------------------------------------------

   ! crit_q() - (q): the gradient self-consistency sweep, one potential per
   !                child process (the harm child rides the container registry
   !                path with the default-constant well; the morse / leps
   !                children ride the transparent direct-bind adapters)
   subroutine crit_q()
      call expect_ok_child('q-harm', 'q-harm: OK', &
         '(q) harm gradient sweep, container path (analytic vs central difference)')
      call expect_ok_child('q-morse', 'q-morse: OK', &
         '(q) morse gradient sweep, direct bind (analytic vs central difference)')
      call expect_ok_child('q-leps', 'q-leps: OK', &
         '(q) leps gradient sweep, direct bind (analytic vs central difference)')
   end subroutine crit_q

   ! crit_s() - (s): the old-vs-new bit-level force comparison, one potential
   !                per child process; the parent stages the legacy keyword
   !                files first (the legacy engine's own PES-selection channel;
   !                the harm well constant rides that channel at 0.04 eV/A^2)
   subroutine crit_s()
      character(len=80) :: lharm(3), lmorse(1), lleps(1)
      lharm = [ character(len=80) :: 'TEST_PES = HARMONIC', &
                                      'HARM_K = 0.04, 0.04, 0.04', &
                                      'HARM_X0 = 0.0, 0.0, 0.0' ]
      lmorse = [ character(len=80) :: 'TEST_PES = MORSE' ]
      lleps = [ character(len=80) :: 'TEST_PES = LEPS' ]
      call write_pars('build/old_in_harm.txt', lharm)
      call write_pars('build/old_in_morse.txt', lmorse)
      call write_pars('build/old_in_leps.txt', lleps)
      call expect_ok_child('s-harm', 's-harm: OK', &
         '(s) old-vs-new bit-level, harm (same Q through both engines)')
      call expect_ok_child('s-morse', 's-morse: OK', &
         '(s) old-vs-new bit-level, morse (same Q through both engines)')
      call expect_ok_child('s-leps', 's-leps: OK', &
         '(s) old-vs-new bit-level, leps (same Q through both engines)')
   end subroutine crit_s

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
      ! buffered), assembly, then a real force_eval
      call registry_gen_keys()                    ! the slot handover runs first (the driver order)
      call registry_gen_all()
      call registry_gen_init()
      call container_init()
      allocate (list_atoms%mass(2))
      list_atoms%mass = [ 1.0d0, 4.0d0 ]
      call state_create(sb, 2)
      ! calibration: one call at a known offset (p = 0) measures the composed
      ! conversion constant (the literal lives only in force_interface)
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
         tk = tk + dot_product(sb%p(3*i - 2:3*i), sb%p(3*i - 2:3*i))/(2.0d0*list_atoms%mass(i))
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
   !                 difference, interface fold, e_tot decomposition
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
      allocate (list_atoms%mass(nm))
      list_atoms%mass = 1.0d0
      call state_create(sb, nm)
      ! the analytic reference from the module face
      if (use_morse) then
         call morse_vg(nm, qa, va, ga)
      else
         call leps_vg(nm, qa, va, ga)
      end if
      ! the interface fold: the measured conversion carries the module gradient
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
            print '(a)', 'FAIL '//trim(tag)//' child: interface-fold force mismatch, dof '// &
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
         tk = tk + dot_product(sb%p(3*i - 2:3*i), sb%p(3*i - 2:3*i))/(2.0d0*list_atoms%mass(i))
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

   ! child_obox() - the (o) child: declare the container keys, read the
   !                case input obox_input.txt from the run directory
   !                (cwd = build/), pull the buffered keys, audit the
   !                consumption, load the module state, print the posture
   subroutine child_obox()
      type(reference_params_t) :: p
      call input_declare_keys(reference_keys())
      call read_input('obox_input.txt')
      call pull('K_HARM',  p%k_harm)
      call pull('Q_TERM',  p%q_term,  p%q_term_found)
      call pull('E_SPLIT', p%e_split, p%e_split_found)
      call input_audit_declared()
      call reference_load(p)
      print '(a)', 'obox: k='//trim(d2s(reference_k()))//' term='// &
         merge('T', 'F', reference_term_on())//' class='// &
         merge('T', 'F', reference_class_on())
      stop 0
   end subroutine child_obox

   !------------------------------------------------------------------
   ! the acceptance-batch children (self-respawn modes; each terminates
   ! itself) + their shared new-side assembly helper
   !------------------------------------------------------------------

   ! bind_new_side(kind) - assemble the new-side force chain for one (q)/(s)
   !                       child: harm rides the container registry path (the
   !                       driver's own registration + assembly seam, default
   !                       well constant); morse / leps ride the transparent
   !                       direct-bind adapters (the (n) posture)
   subroutine bind_new_side(kind)
      integer, intent(in) :: kind            ! 1 = harm, 2 = morse, 3 = leps
      if (kind == 1) then
      call registry_gen_keys()                    ! the slot handover runs first (the driver order)
         call registry_gen_all()
         call registry_gen_init()
      else
         if (kind == 2) then
            call container_bind_pes(morse_wrap)
         else
            call container_bind_pes(leps_wrap)
         end if
         call reg_adiabatic_all()
      end if
      call container_init()
   end subroutine bind_new_side

   ! child_q(kind) - the (q) child: one reference potential swept over its
   !                 three constructed geometries - the module analytic
   !                 gradient vs the CENTRAL DIFFERENCE of the module energy
   !                 at EVERY geometry (h = 1e-6, tol 1e-8 rel with floor
   !                 1 eV/A), the interface fold through the measured conversion
   !                 (1e-12), and the e_tot decomposition with momenta at the
   !                 first geometry; one evidence line per geometry
   subroutine child_q(kind)
      integer, intent(in) :: kind            ! 1 = harm, 2 = morse, 3 = leps
      real(8), parameter :: hcd = 1.0d-6     ! the central-difference step [A]
      character(len=8) :: name
      type(state_t) :: sb
      real(8) :: qa(9), qh(9), ql(9), ga(9), gh(9), gl(9), gcdv(9)
      real(8) :: va, vh, vl, cal, ea, eb, tk, dcmax, dfmax
      integer :: nm, nd, i, j
      select case (kind)
      case (1)
         name = 'q-harm'
      case (2)
         name = 'q-morse'
      case (3)
         name = 'q-leps'
      case default
         print '(a)', 'FAIL child: unknown potential kind '//trim(i2s(kind))
         stop 1
      end select
      call bind_new_side(kind)
      if (kind == 3) then
         nm = 3
      else
         nm = 2
      end if
      nd = 3*nm
      allocate (list_atoms%mass(nm))
      list_atoms%mass = 1.0d0
      call state_create(sb, nm)
      cal = 0.0d0
      do j = 1, 3
         qa = 0.0d0
         ga = 0.0d0
         if (kind == 1) then
            qa(1:6) = qh_list(:, j)
            call harm_vg(nm, qa, va, ga)
         else if (kind == 2) then
            qa(1:6) = qm_list(:, j)
            call morse_vg(nm, qa, va, ga)
         else
            qa(1:9) = ql_list(:, j)
            call leps_vg(nm, qa, va, ga)
         end if
         ! the interface fold at this geometry; the first geometry's call also
         ! measures the composed conversion (one calibration, reused after)
         sb%q(1:nd) = qa(1:nd)
         sb%p = 0.0d0
         call force_eval(sb, ea)
         if (j == 1) then
            cal = -dot_product(sb%f(1:nd), ga(1:nd))/dot_product(ga(1:nd), ga(1:nd))
            if (cal <= 0.0d0) then
               print '(a)', 'FAIL '//trim(name)//' child: measured conversion'// &
                  ' constant non-positive'
               stop 1
            end if
         end if
         dfmax = 0.0d0
         do i = 1, nd
            dfmax = max(dfmax, abs(sb%f(i) + ga(i)*cal)/max(abs(ga(i)*cal), 1.0d0))
            if (abs(sb%f(i) + ga(i)*cal) > 1.0d-12*max(abs(ga(i)*cal), 1.0d0)) then
               print '(a)', 'FAIL '//trim(name)//' child: interface-fold force mismatch,'// &
                  ' geometry '//trim(i2s(j))//', dof '//trim(i2s(i))
               stop 1
            end if
         end do
         ! the analytic gradient vs the central difference of the module energy
         dcmax = 0.0d0
         do i = 1, nd
            qh = 0.0d0
            ql = 0.0d0
            qh(1:nd) = qa(1:nd)
            ql(1:nd) = qa(1:nd)
            qh(i) = qa(i) + hcd
            ql(i) = qa(i) - hcd
            if (kind == 1) then
               call harm_vg(nm, qh, vh, gh)
               call harm_vg(nm, ql, vl, gl)
            else if (kind == 2) then
               call morse_vg(nm, qh, vh, gh)
               call morse_vg(nm, ql, vl, gl)
            else
               call leps_vg(nm, qh, vh, gh)
               call leps_vg(nm, ql, vl, gl)
            end if
            gcdv(i) = (vh - vl)/(2.0d0*hcd)
            dcmax = max(dcmax, abs(ga(i) - gcdv(i))/max(abs(ga(i)), 1.0d0))
            if (abs(ga(i) - gcdv(i)) > 1.0d-8*max(abs(ga(i)), 1.0d0)) then
               print '(a)', 'FAIL '//trim(name)//' child: analytic gradient vs central'// &
                  ' difference, geometry '//trim(i2s(j))//', dof '//trim(i2s(i))//' ('// &
                  trim(d2s(ga(i)))//' vs '//trim(d2s(gcdv(i)))//')'
               stop 1
            end if
         end do
         print '(a,i0,a,es10.3,a,es10.3)', trim(name)//' geom ', j, &
            ': CD max rel dev ', dcmax, ', fold max rel dev ', dfmax
         ! the e_tot decomposition with momenta, first geometry only
         if (j == 1) then
            if (nd == 6) then
               sb%p(1:6) = pm2
            else
               sb%p(1:9) = pm3
            end if
            call force_eval(sb, eb)
            tk = 0.0d0
            do i = 1, nm
               tk = tk + dot_product(sb%p(3*i - 2:3*i), sb%p(3*i - 2:3*i))/(2.0d0*list_atoms%mass(i))
            end do
            tk = tk/e_conv
            if (abs(eb - (tk + va*cal/e_conv)) > 1.0d-12*max(abs(eb), 1.0d0)) then
               print '(a)', 'FAIL '//trim(name)//' child: e_tot decomposition mismatch: '// &
                  trim(d2s(eb))
               stop 1
            end if
         end if
      end do
      print '(a)', trim(name)//': OK'
      stop 0
   end subroutine child_q

   ! child_s(kind) - the (s) child: the SAME Q element-wise into the legacy
   !                 engine's two evaluation entries (driven in THIS process -
   !                 one executable, one double arithmetic; the bit level needs
   !                 it) and into force_eval (same potential, same parameters),
   !                 at p = 0 (T = 0 exactly, e_tot carries V alone): assert
   !                 max|dF| = max|dE| = 0 bitwise, the 1e-15 fallback with the
   !                 printed value and cause for a benign arithmetic-order
   !                 difference, a FAIL beyond that (the disposition plan)
   subroutine child_s(kind)
      use input_parser, only: load_input
      use venus_data, only: natoms, Q, PDOT, NSURF, NGLO
      integer, intent(in) :: kind            ! 1 = harm, 2 = morse, 3 = leps
      character(len=8) :: name
      character(len=24) :: oldfile
      type(state_t) :: sb
      real(8) :: qarr(1)                     ! the legacy dummy (it evaluates on Q)
      real(8) :: vpot, en, df, de
      integer :: nm, nd, j
      logical :: bitw
      select case (kind)
      case (1)
         name = 's-harm'
         oldfile = 'build/old_in_harm.txt'
      case (2)
         name = 's-morse'
         oldfile = 'build/old_in_morse.txt'
      case (3)
         name = 's-leps'
         oldfile = 'build/old_in_leps.txt'
      case default
         print '(a)', 'FAIL child: unknown potential kind '//trim(i2s(kind))
         stop 1
      end select
      ! the legacy side first, through its own call convention: the keyword
      ! file the parent buffered, load_input, POTPRE (the PES selection), then
      ! the live-array posture the entries evaluate on
      call load_input(trim(oldfile))
      call potpre()
      NSURF = 0
      NGLO = 0
      if (kind == 3) then
         nm = 3
      else
         nm = 2
      end if
      nd = 3*nm
      natoms = nm
      allocate (Q(nd), PDOT(nd))
      ! the new side: the same assembly posture as the (q) child
      call bind_new_side(kind)
      allocate (list_atoms%mass(nm))
      list_atoms%mass = 1.0d0
      call state_create(sb, nm)
      do j = 1, 3
         if (kind == 1) then
            sb%q(1:6) = qh_list(:, j)
         else if (kind == 2) then
            sb%q(1:6) = qm_list(:, j)
         else
            sb%q(1:9) = ql_list(:, j)
         end if
         sb%p = 0.0d0                        ! T = 0 exactly: e_tot carries V alone
         call force_eval(sb, en)
         Q(1:nd) = sb%q(1:nd)                ! the same Q, element-wise
         call POT0(nd, vpot, qarr)           ! Vpot = v_ev x 23.0605 [kcal/mol]
         PDOT = 0.0d0                        ! the accumulator zeroed: PDOT = -g x 23.0605 x C1
         call DPESHON(nd, qarr)
         bitw = all(sb%f(1:nd) == PDOT(1:nd)) .and. en == vpot
         df = maxval(abs(sb%f(1:nd) - PDOT(1:nd)))
         de = abs(en - vpot)
         print '(a,i0,a,l2,a,es12.4,a,es12.4)', trim(name)//' geom ', j, ': bitwise = ', bitw, &
            ', max|dF| = ', df, ', max|dE| = ', de
         if (.not. bitw) then
            if (all(abs(sb%f(1:nd) - PDOT(1:nd)) <= 1.0d-15*max(1.0d0, abs(PDOT(1:nd)))) .and. &
                de <= 1.0d-15*max(1.0d0, abs(vpot))) then
               print '(a)', trim(name)//' note: geometry '//trim(i2s(j))//' - bitwise'// &
                  ' attempt failed, deviation within 1e-15 (arithmetic-order difference,'// &
                  ' value above, cause recorded - the disposition plan applied)'
            else
               print '(a)', 'FAIL '//trim(name)//' child: geometry '//trim(i2s(j))// &
                  ' deviates beyond 1e-15 (max|dF| = '//trim(d2s(df))//', max|dE| = '// &
                  trim(d2s(de))//') - locate the difference source (port error vs'// &
                  ' arithmetic order) per the disposition plan'
               stop 1
            end if
         end if
      end do
      print '(a)', trim(name)//': OK'
      stop 0
   end subroutine child_s

   ! pot_kind(name) - the (q)/(s) child's potential selector from the mode tail
   function pot_kind(name)
      character(len=*), intent(in) :: name
      integer :: pot_kind
      select case (trim(name))
      case ('harm')
         pot_kind = 1
      case ('morse')
         pot_kind = 2
      case ('leps')
         pot_kind = 3
      case default
         pot_kind = 0
      end select
   end function pot_kind

   ! the transparent (n) adapters: the module potentials (v, g signature)
   ! delivered through the container force entry (q, mass, g, v) - the module
   ! material itself, bound check-side
   subroutine morse_wrap(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A]
      real(8), intent(out) :: v        ! V [eV]
      call morse_vg(size(q)/3, q, v, g)
   end subroutine morse_wrap

   subroutine leps_wrap(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A]
      real(8), intent(out) :: v        ! V [eV]
      call leps_vg(size(q)/3, q, v, g)
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
            'interface_ok.txt 2>&1'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp_local('build/interface_ok.txt', buf)
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
      cmd = 'cd build && .'//sep//'check_facade.exe '//trim(cmode)//' >interface_ok.txt 2>&1'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp_local('build/interface_ok.txt', buf)
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

   ! write one case-input file (one KEYWORD = VALUE per record)
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
