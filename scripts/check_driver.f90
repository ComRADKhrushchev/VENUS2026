!=====================================================================
! check_driver.f90 - criteria program for the whole driver (assembly ->
!                    evolution -> statistics; criteria defined FIRST in
!                    this header, this program executes them)
! Compile (from src/; run make FIRST - the driver exe and the generated registry
! module are build products, and make's prebuild step aggregates the dragged-in
! systems/ container into that module; -Jbuild keeps the .mod files out of the
! source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 utils/rng.f90 &
!            globals/state.f90 globals/config_atoms.f90 &
!            interface/elec_interface.f90 utils/linalg.f90 utils/spectrum.f90 &
!            utils/hessian.f90 interface/force_interface.f90 interface/spectrum_interface.f90 &
!            utils/geom_init.f90 dynamics/propagator.f90 interface/bath.f90 &
!            interface/sampler.f90 utils/geometry.f90 &
!            methods/samp_boltzmann.f90 methods/samp_ebk.f90 methods/beam_laws.f90 methods/inc_surface.f90 methods/inc_pair.f90 methods/inc_single.f90 &
!            methods/samp_j.f90 methods/samp_rotation.f90 methods/samp_surface_oscillator.f90 &
!            methods/samp_normalmode.f90 methods/samp_barrier_excitation.f90 &
!            methods/bath_andersen.f90 methods/samp_glo_target.f90 &
!            methods/samp_thermalize.f90 methods/samp_wigner.f90 methods/elec_adiabatic.f90 &
!            methods/seam_cast.f90 &
!            methods/reg_adiabatic.f90 methods/reg_prop.f90 &
!            methods/reg_samp_boltzmann.f90 methods/reg_samp_ebk.f90 &
!            methods/reg_inc_surface.f90 methods/reg_inc_pair.f90 methods/reg_inc_single.f90 methods/reg_samp_j.f90 &
!            methods/reg_samp_rotation.f90 methods/reg_samp_surface_oscillator.f90 &
!            methods/reg_samp_normalmode.f90 methods/reg_samp_barrier_excitation.f90 &
!            methods/reg_bath_andersen.f90 methods/reg_samp_glo_target.f90 &
!            methods/reg_samp_thermalize.f90 methods/reg_samp_wigner.f90 &
!            output/recorder.f90 output/final_state.f90 &
!            output/plot_hist.f90 &
!            ../systems/reference/reference_pot.f90 &
!            ../systems/reference/reference_box.f90 &
!            interface/container_sched.f90 &
!            ../systems/reference/reg_reference.f90 &
!            build/gfortran/registry_gen.f90 ../scripts/check_driver.f90 &
!            -o build/check_driver.exe
! Run (from src/): build\check_driver.exe
!
! Observables this check greps (the design choices, stated per the brief):
!   - the ASSEMBLY-SUMMARY line the driver prints at its assembly-phase end:
!     'assembly complete: trajectories=<n> atoms=<n> fragments=<n> paradigm=<n>
!     members=<n> init=ok columns=<n>' (columns = the selected recorder column
!     count through rec_n_sel(); the assembly-phase field set stays the same -
!     the count face named for what it counts).
!   - the STATISTICS-SUMMARY block the driver prints after the last trajectory:
!     'statistics: trajectories=<n> valid=<n>' / 'outcome counts: <name>=<n> ...'
!     or 'outcome counts: none (no classification bound - skipped)' /
!     'conservation drift: trajectories=<n> max|e_fin-e0|=<v> kcal/mol' /
!     'timing: cpu=<v> s'. The drift's data source is the driver's choice: the
!     recorder final-state row (rec_fin_val('e_tot') - force_eval at the final
!     geometry, a real reading), stated in the driver.
!   - the OUTPUT FILE 'trajectories.txt' (cwd-relative, the driver's registered
!     recording file): '# frame traj=<t> step=<s> t=<time>' headers + the
!     selected-column line 't_kin=.. v_pot=.. e_tot=..' + gated '# sec_qp'
!     flattened sections + one '# final:' row group. This check parses it.
! Criteria:
!   (a) END-TO-END THROUGH (the real driver exe, child process): a real input
!       file + a real two-fragment system folder (two diatomic .xyz files -
!       the generic hypothetical molecular system) + the dragged-in reference
!       container -> exit 0, stdout carries the assembly-summary
!       line with the expected field values (trajectories=3 atoms=4 fragments=2
!       paradigm=0 members=11 columns=3) AND the statistics-summary block, and
!       trajectories.txt exists and is non-empty
!   (b) END-TO-END ENERGY CONSERVATION (the whole chain's first conservation
!       criterion): in (a)'s trajectories.txt every frame's e_tot stays within
!       the stated tolerance of its trajectory's FIRST frame's e_tot (the head
!       force_eval e0 and the first frame's e_tot run the same arithmetic on
!       the same state, so first-frame constancy anchors e0 transitively).
!       Tolerance: the velocity-Verlet energy excursion for a harmonic mode
!       scales ~(omega*dt)^2 of the energy scale; with DT=0.05, K_HARM=0.04
!       eV/A^2 and unit-mass atoms omega ~ 0.20 per 10 fs, so (omega*dt)^2 ~
!       1.0e-4 of ~70 kcal/mol ~ 7e-3 kcal/mol; the criterion fixes
!       5.0e-2 kcal/mol (~7x that bound). The measured max prints in the PASS
!       line
!   (r) RADAU END-TO-END (the coverage gap this criterion closes: the
!       adaptive member's guards and cross-sequence carry had unit-level
!       coverage only, never through the real four-phase loop): the (a)
!       case material with INTEGRATOR=RADAU, ORDER=8 (ss = 1e-8) ->
!       exit 0, the two summary blocks, 2 trajectories x 41 frames
!       (safety-net-only, the (c) posture), and every frame's e_tot
!       within the radau tolerance of its trajectory's first frame -
!       the boot rule fires per trajectory (the fresh-state reset), the
!       B-carry crosses sequences, and the close-out state guard runs
!       on every trial, all on the real force chain. Tolerance
!       derivation (the (b) precedent's shape): RA15 controls the
!       per-sequence error at ss; the unit-level long-run anchor (Morse
!       well, order 8, hundreds of sequences) measures max|E-E0| ~ 1e-9
!       kcal/mol with its window fixed at 1e-7; this run is 40
!       sequences per trajectory on the smoother harmonic chain, so the
!       criterion fixes 1.0e-6 kcal/mol - 10x the unit window, ~4
!       decades tighter than the fixed-step tolerance 5e-2 (adaptive
!       control, not step-size luck). The measured max prints in the
!       PASS line
!   (c) FRAME STRUCTURE: frame headers carry traj/step/t with traj strictly
!       1..N, step strictly 1..M+1 in file order, t strictly increasing within
!       each trajectory; the column line's order is t_kin, v_pot, e_tot (the
!       recorder's selected-column legislation); frames per trajectory =
!       MAX_STEPS+1 (the MAX_STEPS-only run - this case binds no termination)
!   (d) MAX_STEPS SAFETY NET: a small MAX_STEPS input (5, one trajectory) ->
!       exactly 6 frames in trajectories.txt, the last at step 6 with t>0
!   (e) TERM-BOUND PATH: the container's parameter file carries Q_TERM (the
!       coordinate-threshold termination) -> the run terminates EARLY (frames
!       per trajectory < MAX_STEPS+1 and >= 2, the last frame's t below the
!       safety-net horizon MAX_STEPS*DT) - the container's termination judgment
!       observable through the frame count and time
!   (f) STATISTICS FIRST VERSION: the no-classification run prints the skipped
!       semantics line 'outcome counts: none (no classification bound -
!       skipped)' and the conservation-drift summary line; the
!       classification-bound run (E_SPLIT = 1.0 in the parameter file) prints
!       the outcome table over the REAL final-state values: E_REL = 200
!       kcal/mol leaves the pair's final relative radial energy two orders
!       above the split at the Q_TERM termination point (the weak well
!       K = 0.04 barely deflects the fly-by), so both trajectories classify
!       'scattering' -> 'outcome counts: scattering=2'. Expectation migrated
!       from the earlier 'capture=2', which read the zero-initialized fin
!       before the final-state quantity computation existed (e_rel = 0 sits
!       below any positive split); the real-value margin to the threshold is
!       two orders of magnitude, not a knife-edge
!   (g) SAMPLING CHAIN: (a)'s run drives samp_run through the real-input
!       DIST_SCHEME (the thermalize member on fragment 1, the rotation member
!       on fragment 2 - both on the real force chain) and the incident channel
!       for the two-fragment gas pair; observable: the FIRST frame's sec_qp
!       section carries non-zero momenta (the sampling really happened; the
!       case input's GWRITE_LEVEL=2 row is retained as documentation - it now
!       restates the bare-input default, not a departure from level 1)
!   (det) DETERMINISM: the same input run twice -> trajectories.txt
!       byte-identical (the rng single stream + the fixed arithmetic; compared
!       with fc /b on Windows, cmp -s elsewhere)
!   (E7) PHASE-RESTRICTION SOURCE ASSERTION (the method stated per the brief):
!       read src/driver.f90 as text, locate the integration-loop head (the
!       do-while line naming term_flag), scan to its end do, skip blank and
!       comment lines, collect every 'call <name>' target PLUS the
!       container_term occurrence (the phase-4 entry is a function assignment,
!       not a call statement), and require the collected multiset to equal
!       exactly {prop_step, elec_prop, mqc_statistic, rec_frame,
!       container_term} - nothing else inside the integration loop.
!       LIMIT: the scan stops at the FIRST 'end do' - a nested do inside
!       the integration loop would truncate the collected set (the
!       current loop body is unnested; the limit is declared here)
!   (h) assembly-order error named abort (child): a bad SYSTEM_DIR (the folder
!       does not exist) -> nonzero exit, stderr names sysdef_load - the
!       existing input/sysdef channel penetrates; the driver adds NO judgment
!       of its own
!   (i) size guard: UNREACHABLE through the driver's public face - the driver
!       derives state_create's natoms from size(list_atoms%mass) itself, so
!       size(sta%q) = 3*size(list_atoms%mass) holds structurally; the guard is
!       carried in the driver per the guard rule and names both sides if it
!       ever fires (the untriggerable-guard precedent: header note here, no
!       criterion body)
!   (j) seam criterion (the input-buffer -> member-params assembly channel):
!       (j1) round-trip, IN-PROCESS through the same public chain in the same
!            fixed assembly order (registry_gen_keys -> read_input -> rng_init
!            -> sysdef_load -> list_atoms_load -> state_create -> registry_gen_all
!            -> registry_gen_init -> input_audit_declared -> container_init ->
!            rec_select level entry):
!            the carrier bracket GLO_TARGET[T_GLO=350, GAMMA_GLO=0.04, W_GHOST=1.25] / DT=
!            0.05 / DIST_SCHEME=GLO_TARGET -> the member getters read the
!            exact values back (glo_target_t_glo / _gamma / _w), and
!            glo_target_dt_loaded() = DT (the dt-ordering contract); also
!            samp_n_dist() = 12 and rec_n_sel() = 3 through the same chain;
!            and container_term(sta) = .false. (no container keys buffered in
!            this case - the termination slot unbound, the legal
!            safety-net-only posture)
!       (j2) bad value (child, the real driver exe): T_GLO = -5 -> nonzero
!            exit, stderr names glo_target_init and the field t_glo
!   (k) unbound-classification authority (this check's OWN child mode
!       'eclass', self-respawn): the full assembly chain in a run directory
!       WITHOUT the container parameter file (termination and classification
!       unbound), then container_classify(fin) -> nonzero exit, stderr names
!       container_classify. (The unbound-PES posture has no criterion any
!       more: the verification container is dragged into the build on every
!       compile, so a container-free registry is no longer instantiable
!       through the real build - that posture was verified while the tree was
!       container-free; the classification slot keeps the fail-loud witness.)
!   (l) HISTOGRAM FILES (the statistics-phase data products): the term+classify
!       case's run directory carries 'hist_e_rel.dat' and 'hist_thta.dat' (the
!       driver's registered statistics output names), each 40 bin lines + one
!       TOTAL & P line; bin centers constant-spaced with the first at half the
!       spacing (the driver-side range rule [0, column max + 1) puts every row
!       in range by construction); the e_rel file counts BOTH trajectories
!       (TOTAL = 2, P = 2/2 = 1); the theta file's truth is fully derived: the
!       gas-phase run leaves the surface scattering angle bitwise 0, so every
!       row lands in bin 1 (bin-1 count 2, all other bins 0, spacing exactly
!       0.025 = 1.0/40)
!   (m) SUMMARY TRUTH: stdout carries the histogram line 'histograms: rows=2
!       files=hist_e_rel.dat hist_thta.dat' after the drift line; and the
!       PRINTED conservation drift matches the file-derived drift - per
!       trajectory, |final-row e_tot - first-frame e_tot| read back from
!       trajectories.txt, max over trajectories, within 1e-3 + 1e-3*|derived|
!       (the printed field is es12.4, the file's e_tot digits es24.16 - the
!       tolerance covers only the print width)
!   (n) OCCUPATION COLUMNS END TO END: OBSERVABLES = occ,n_hop on a
!       two-trajectory safety-net case -> the assembly line reads columns=5,
!       and EVERY frame column line and every final row in trajectories.txt
!       carries occ = 1.0 and n_hop = 0.0 - the degenerate adiabatic values
!       through the whole chain (input key -> selection union -> extractor ->
!       getter -> file)
!   (phg) RANGE GUARD (child 'phguard', self-respawn): plot_hist on a spec
!       with v_hi <= v_lo (kind 1) -> nonzero exit, stderr names plot_hist -
!       the bin-width choke point guards every caller
!   (q) NON-LINEAR FRAGMENT ANCHOR (child 'nlfrag', self-respawn; the
!       rotation-strip anchor the final-state criteria program deferred - its
!       constructed fragments were all linear): a real 3-atom non-collinear
!       fragment (the centered H3 triangle buffered in the case folder, real
!       list_atoms masses), pure rigid rotation about z (p_i = m_i*(omega x r_i),
!       omega = 0.07 z-hat): fin_eval must recover j_a = (0,0,I_zz*omega) in
!       hbar units, e_rot_a = (I_zz*omega^2/2) in kcal/mol, the internal
!       kinetic e_int_a(1) EQUAL to it (two computational paths to one number
!       - the cross-identity), the archive rotation STRIPPED (max|p_fin| <=
!       1e-10 against a p scale of ~0.06 - THE non-linear anchor; linear
!       fragments keep their rotation by design), q_fin the COM-shifted
!       triangle, e_int_a(2) the buffered well on a check-buffered constant
!       (calibrated through the bound container), e_cm bitwise 0 (a pure
!       rotation carries no net momentum) and the pair/B fields bitwise 0
!       (the single-fragment skip)
! Child-process scheme (the sibling checks' discipline): every named abort path
!   runs in a child - driver-exe children run with cwd = their case directory
!   (the driver reads ./input_qct.txt relative to its cwd; 'cd <case> && exe'
!   carries both shell families); the setup copies the build product to
!   build/drv_venus.exe once before spawning, PROBING BOTH COMPILER TREES
!   (build/gfortran/venus.e first, then build/ifx/venus.e; neither present ->
!   a NAMED setup failure - 'run make first' - instead of a silent copy
!   failure). The check's own child ('eclass') runs from src/ and reads its
!   input by build-relative path.
! Test materials live ONLY here and under build/ - the driver holds zero
!   system knowledge; every case folder is written fresh by make_cases().
!=====================================================================
program check_driver
   use control,       only: i_seed
   use consts,        only: dyn_ev, hbar_act
   use config,        only: reactants, observables, electronic, prop_cfg => propagator
   use config_atoms, only: list_atoms, list_atoms_load
   use state,         only: state_t, state_create
   use rng,           only: rng_init
   use sysdef_interface, only: sysdef_load
   use registry_gen,  only: registry_gen_all, registry_gen_init, registry_gen_keys
   use force_interface,  only: container_init, container_term, container_bind_pes, force_eval
   use sampler,       only: samp_n_dist
   use recorder,      only: rec_select, rec_n_sel
   use final_state,   only: fin_t, fin_eval, container_classify
   use hist,          only: hist_spec_t, hist_data_t, plot_hist, k_e_hist
   use elec_interface,   only: method_bind
   use reg_adiabatic, only: reg_adiabatic_all
   use reference_pot, only: harm_set_constant
   use reference_box, only: reference_pes
   use input,         only: read_input, input_audit_declared
   use samp_glo_target, only: glo_target_t_glo, glo_target_gamma, glo_target_w, &
                              glo_target_dt_loaded
   implicit none
   integer :: nfail = 0
   character(len=32) :: cmode
   character(len=2048) :: obuf_a = ''    ! case_a's stdout, kept for (f)
   character(len=2048) :: obuf_t = ''    ! case_t's stdout, kept for (f)

   ! ---------- child modes (self-respawn): the named abort paths that cannot run here ----------
   if (command_argument_count() >= 1) then
      call get_command_argument(1, cmode)
      select case (trim(cmode))
      case ('eclass')
         call child_eclass()
      case ('phguard')
         call child_phguard()
      case ('nlfrag')
         call child_nlfrag()
      case default
         print '(a)', 'check_driver: unknown child mode '//trim(cmode)
         stop 2
      end select
      print '(a)', 'check_driver: child mode '//trim(cmode)//' returned (it must stop)'
      stop 2
   end if

   ! ---------- case material under build/ ----------
   call make_cases()

   ! ---------- criteria ----------
   call crit_a()      ! (a) end-to-end through + the two summary blocks
   call crit_b()      ! (b) end-to-end energy conservation
   call crit_r()      ! (r) radau end to end (completion + adaptive conservation)
   call crit_c()      ! (c) frame structure (MAX_STEPS-only run)
   call crit_g()      ! (g) the sampling chain really ran
   call crit_d()      ! (d) MAX_STEPS safety net
   call crit_e()      ! (e) term-bound early termination (+ case_t stdout for (f))
   call crit_f()      ! (f) statistics first version (unbound + bound)
   call crit_det()    ! (det) determinism, same input twice
   call crit_e7()     ! (E7) the integration loop's direct-call set
   call crit_h()      ! (h) bad SYSTEM_DIR named abort
   call crit_j1()     ! (j1) seam round-trip + unbound container_term
   call crit_j2()     ! (j2) bad member value
   call crit_k()      ! (k) unbound-classification child
   call crit_l()      ! (l) the histogram files (existence/format/count truth)
   call crit_m()      ! (m) the summary lines (histogram line + drift vs file)
   call crit_n()      ! (n) occ/n_hop end to end through the real driver
   call crit_phg()    ! (phg) the plot_hist range guard (child)
   call crit_q()      ! (q) the non-linear fragment anchor (child)

   if (nfail == 0) then
      print '(a)', 'check_driver: ALL PASS'
   else
      print '(a,i0,a)', 'check_driver: ', nfail, ' FAILURE(S)'
      stop 1
   end if
contains
   !------------------------------------------------------------------
   ! child_eclass() - the (k) child: the full assembly chain on a system folder
   !                  WITHOUT container keys (termination and classification
   !                  unbound), then container_classify - the unbound call must
   !                  be the named-abort authority
   subroutine child_eclass()
      type(state_t) :: stc
      type(fin_t) :: fc
      character(len=32) :: nm
      call registry_gen_keys()
      call read_input('build/drv_case_a/input_eclass.txt')
      call rng_init(i_seed(1))
      call sysdef_load(reactants%system_dir)
      call list_atoms_load()
      call state_create(stc, size(list_atoms%mass), 1)   ! the degenerate adiabatic surface count
      call registry_gen_all()
      call registry_gen_init()     ! container keys absent -> slots stay unbound
      call input_audit_declared()
      call container_init()
      nm = container_classify(fc)  ! the unbound call - must named abort here
      print '(a)', 'child eclass: container_classify returned ', trim(nm)
   end subroutine child_eclass

   !------------------------------------------------------------------
   ! child_phguard() - the (phg) child: plot_hist on an inverted range
   !                  (v_hi < v_lo, kind 1) - the bin-width choke point must
   !                  named abort here
   subroutine child_phguard()
      type(hist_spec_t) :: sh
      type(hist_data_t) :: dh
      allocate (dh%e_rel(1))
      dh%e_rel = 1.0d0
      dh%n_traj = 1
      sh = hist_spec_t(k_e_hist, 'build/drv_phg.dat', 2.0d0, 1.0d0, 4, 0, .true.)
      call plot_hist(sh, dh)       ! the inverted range - must named abort here
   end subroutine child_phguard

   !------------------------------------------------------------------
   ! child_nlfrag() - the (q) child: the real 3-atom non-collinear list_atoms (the
   !                  centered H3 triangle of the case folder, real masses)
   !                  driven through the analysis chain, one pure rigid
   !                  rotation about z constructed on it, fin_eval, and every
   !                  anchor asserted (derivations: the header criterion (q))
   subroutine child_nlfrag()
      type(state_t) :: stc
      type(fin_t) :: fc
      real(8), parameter :: kk = 0.123d0   ! the check-buffered well constant [eV/A^2]
      real(8), parameter :: omg = 0.07d0   ! the constructed rotation rate [rad/(10 fs)]
      real(8) :: rc(3), m3(3), ii_zz, ee, cal, sq, pmax
      integer :: a, d
      call read_input('build/drv_case_n/input_nlfrag.txt')
      call rng_init(i_seed(1))
      call sysdef_load(reactants%system_dir)
      call list_atoms_load()
      ! the analysis chain: degenerate method row + the reference PES + the
      ! check's own well constant (the buffered-potential anchor rides the buffered
      ! value, not the package default)
      call reg_adiabatic_all()
      call method_bind('adiabatic')
      call container_bind_pes(reference_pes)
      call harm_set_constant(kk)
      call state_create(stc, size(list_atoms%mass), 1)   ! the degenerate adiabatic surface count
      ! conversion calibration: one force_eval at a single-offset geometry with
      ! p = 0 (V = kk/2 * |q|^2 eV) measures the kcal/mol-per-eV factor through
      ! the bound container - the conversion literal never lands here
      stc%q = 0.0d0
      stc%p = 0.0d0
      stc%f = 0.0d0
      stc%q(1) = 1.0d0
      call force_eval(stc, ee)
      cal = 2.0d0*ee/kk
      ! the constructed pure rotation on the real triangle: rc = COM(qz),
      ! I_zz about z, p_i = m_i*(omega x r_i) with omega = omg*z-hat (no COM
      ! velocity - the rotation is the whole motion)
      rc = 0.0d0
      m3 = list_atoms%mass(1:3)
      do a = 1, 3
         rc(:) = rc(:) + m3(a)*list_atoms%frag(1)%qz(3*a-2:3*a)
      end do
      rc = rc(:)/sum(m3)
      ii_zz = 0.0d0
      sq = 0.0d0
      do a = 1, 3
         stc%q(3*a-2:3*a) = list_atoms%frag(1)%qz(3*a-2:3*a)
         stc%p(3*a-2) = -m3(a)*omg*(stc%q(3*a-1) - rc(2))
         stc%p(3*a-1) =  m3(a)*omg*(stc%q(3*a-2) - rc(1))
         stc%p(3*a)   = 0.0d0
         ii_zz = ii_zz + m3(a)*((stc%q(3*a-2) - rc(1))**2 + (stc%q(3*a-1) - rc(2))**2)
      end do
      do d = 1, 9
         sq = sq + list_atoms%frag(1)%qz(d)**2
      end do
      stc%f = 0.0d0
      stc%t = 0.4d0
      call fin_eval(stc, fc)
      ! ---- anchors ----
      call chk_nl('(q) j_a(3)', fc%j_a(3), ii_zz*omg/hbar_act)
      call chk_z('(q) j_a(1)', fc%j_a(1))
      call chk_z('(q) j_a(2)', fc%j_a(2))
      call chk_nl('(q) j_mag_a', fc%j_mag_a, ii_zz*omg/hbar_act)
      call chk_nl('(q) e_rot_a [inertia-tensor path]', fc%e_rot_a, &
                  0.5d0*ii_zz*omg*omg*dyn_ev)
      call chk_nl('(q) e_int_a(1) [COM-frame path - the cross-identity]', fc%e_int_a(1), &
                  0.5d0*ii_zz*omg*omg*dyn_ev)
      call chk_nl('(q) e_int_a(2) [buffered well, check-buffered constant]', fc%e_int_a(2), &
                  cal*(kk/2.0d0)*sq)
      call chk_z('(q) e_cm [pure rotation carries no net momentum]', fc%e_cm)
      call chk_z('(q) e_rel [single fragment]', fc%e_rel)
      call chk_z('(q) oam_mag [single fragment]', fc%oam_mag)
      call chk_z('(q) b_fin [single fragment]', fc%b_fin)
      call chk_z('(q) ang(19) [gas: no surface pair]', fc%ang(19))
      call chk_z('(q) e_int_b(1) [no second fragment]', fc%e_int_b(1))
      call chk_nl('(q) t_life', fc%t_life, 0.4d0)   ! 4.0 x 10-fs units -> 0.4 fs (2026-10-06 rebase)
      ! the archive: rotation stripped (THE non-linear anchor) + COM-shifted rows
      pmax = 0.0d0
      do d = 1, 9
         pmax = max(pmax, abs(fc%p_fin(d)))
      end do
      if (pmax > 1.0d-10) then
         call fail('(q) the non-linear fragment archive must strip the rotation:'// &
                   ' max|p_fin| = '//d2s(pmax)//' (the p scale is ~0.06)')
      end if
      do a = 1, 3
         do d = 1, 3
            call chk_nl('(q) q_fin (COM-shifted triangle)', fc%q_fin(3*(a-1) + d), &
                        list_atoms%frag(1)%qz(3*(a-1) + d) - rc(d))
         end do
      end do
      if (nfail == 0) then
         print '(a)', 'nlfrag: non-linear fragment anchors hold (rotation stripped,'// &
            ' j/e_rot/kinetic cross-identity, buffered well, COM-shifted archive)'
         stop 0
      end if
      stop 1
   end subroutine child_nlfrag

   !------------------------------------------------------------------
   ! crit_a() - (a): the real driver exe end to end through the two-fragment
   !             case; the assembly line, the statistics block, the output file
   subroutine crit_a()
      character(len=2048) :: obuf, ebuf
      character(len=256) :: want
      character(len=256) :: fline
      integer :: stat, v, ios
      call run_driver('drv_case_a', obuf, ebuf, stat)
      obuf_a = obuf
      want = 'assembly complete: trajectories=3 atoms=4 fragments=2 paradigm=0'// &
             ' members=11 init=ok columns=3'
      if (stat /= 0) then
         call fail('(a): the driver child exited nonzero - stderr: '//trim(ebuf))
         return
      else if (index(obuf, trim(want)) == 0) then
         call fail('(a): stdout missing the assembly line "'//trim(want)//'" - got "'// &
                  trim(obuf)//'"')
         return
      end if
      if (index(obuf, 'statistics: trajectories=3 valid=3') == 0) then
         call fail('(a): stdout missing the statistics count line - got "'//trim(obuf)//'"')
         return
      end if
      if (index(obuf, 'conservation drift: trajectories=3') == 0 .or. &
          index(obuf, 'timing: cpu=') == 0) then
         call fail('(a): stdout missing the drift/timing lines - got "'//trim(obuf)//'"')
         return
      end if
      open (newunit=v, file='build/drv_case_a/trajectories.txt', status='old', &
            action='read', iostat=ios)
      if (ios /= 0) then
         call fail('(a): trajectories.txt was not created in the case directory')
         return
      end if
      read (v, '(a)', iostat=ios) fline
      close (v)
      if (ios /= 0 .or. len_trim(fline) == 0) then
         call fail('(a): trajectories.txt is empty')
         return
      end if
      print '(a)', '(a): PASS - the real driver ran 3 trajectories through the real'// &
           ' input, the two-fragment system and the dragged-in container; both'// &
           ' summary blocks and trajectories.txt verified'
   end subroutine crit_a

   !------------------------------------------------------------------
   ! crit_b() - (b): every frame's e_tot within tolerance of its trajectory's
   !             first frame's e_tot (the head e0 anchored transitively)
   subroutine crit_b()
      real(8), parameter :: tol = 5.0d-2   ! kcal/mol - the stated VV bound ~7e-3, ~7x margin
      integer, parameter :: max_traj = 8
      real(8) :: e_first(max_traj), dev_max(max_traj)
      integer :: ntr(max_traj), cur
      call scan_columns('build/drv_case_a/trajectories.txt', max_traj, ntr, e_first, &
                        dev_max, cur)
      if (cur < 1 .or. any(ntr(1:cur) == 0)) then
         call fail('(b): no frames parsed from trajectories.txt')
         return
      end if
      if (any(dev_max(1:cur) > tol)) then
         call fail('(b): energy drift exceeds the stated tolerance '//d2s(tol)// &
                  ' kcal/mol - per-trajectory max: '//trim(dev_str(dev_max, cur)))
         return
      end if
      print '(a)', '(b): PASS - e_tot constant across all frames of all '//i2s(cur)// &
           ' trajectories; max deviation from the first frame '//d2s(maxval(dev_max(1:cur)))// &
           ' kcal/mol (tolerance '//d2s(tol)//', the stated VV bound ~7e-3)'
   end subroutine crit_b

   !------------------------------------------------------------------
   ! crit_r() - (r): the radau end-to-end run - completion + conservation
   !             inside the stated adaptive tolerance (derivation: header (r))
   subroutine crit_r()
      character(len=2048) :: obuf, ebuf
      character(len=256) :: want
      integer :: stat, ntraj, nfr, last_s
      real(8) :: last_t
      real(8), parameter :: tol = 1.0d-6   ! kcal/mol - the stated RA15 bound (header (r))
      integer, parameter :: max_traj = 8
      real(8) :: e_first(max_traj), dev_max(max_traj)
      integer :: ntr(max_traj), cur
      call run_driver('drv_case_r', obuf, ebuf, stat)
      if (stat /= 0) then
         call fail('(r): the radau driver child exited nonzero - stderr: '//trim(ebuf))
         return
      end if
      want = 'assembly complete: trajectories=2 atoms=4 fragments=2 paradigm=0'// &
             ' members=11 init=ok columns=3'
      if (index(obuf, trim(want)) == 0) then
         call fail('(r): stdout missing the assembly line "'//trim(want)//'" - got "'// &
                  trim(obuf)//'"')
         return
      end if
      if (index(obuf, 'statistics: trajectories=2 valid=2') == 0 .or. &
          index(obuf, 'conservation drift: trajectories=2') == 0) then
         call fail('(r): stdout missing the statistics block - got "'//trim(obuf)//'"')
         return
      end if
      call traj_summary('build/drv_case_r/trajectories.txt', ntraj, nfr, last_s, last_t)
      if (ntraj /= 2 .or. nfr /= 41) then
         call fail('(r): expected 2 trajectories x 41 frames (MAX_STEPS+1, the'// &
                  ' safety-net-only posture), got '//i2s(ntraj)//' x last-count '//i2s(nfr))
         return
      end if
      call scan_columns('build/drv_case_r/trajectories.txt', max_traj, ntr, e_first, &
                        dev_max, cur)
      if (cur < 2 .or. any(ntr(1:cur) == 0)) then
         call fail('(r): no frames parsed from the radau trajectories.txt')
         return
      end if
      if (any(dev_max(1:cur) > tol)) then
         call fail('(r): radau energy drift exceeds the stated tolerance '//d2s(tol)// &
                  ' kcal/mol - per-trajectory max: '//trim(dev_str(dev_max, cur)))
         return
      end if
      print '(a)', '(r): PASS - the adaptive member ran the real four-phase loop'// &
           ' end to end (boot per trajectory, B-carry across sequences, state guard'// &
           ' on every trial); max e_tot deviation '//d2s(maxval(dev_max(1:cur)))// &
           ' kcal/mol (tolerance '//d2s(tol)//', the stated RA15 bound; the'// &
           ' fixed-step tolerance 5e-2 for contrast)'
   end subroutine crit_r

   !------------------------------------------------------------------
   ! crit_c() - (c): frame header monotonicity, column order, frames =
   !             MAX_STEPS+1 per trajectory (the MAX_STEPS-only run)
   subroutine crit_c()
      character(len=512) :: line
      integer :: v, ios, tj, sj, cur_t, cur_s, itraj
      integer :: nfr, bad
      real(8) :: tv, t_prev
      logical :: expect_col
      integer :: p1, p2, p3
      cur_t = 0
      cur_s = 0
      nfr = 0
      bad = 0
      itraj = 0
      expect_col = .false.
      open (newunit=v, file='build/drv_case_a/trajectories.txt', status='old', &
            action='read', iostat=ios)
      if (ios /= 0) then
         call fail('(c): trajectories.txt missing')
         return
      end if
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, '# frame traj=') == 1) then
            call parse_frame(line, tj, sj, tv)
            if (tj /= cur_t) then          ! a new trajectory begins
               if (cur_t > 0 .and. cur_s /= 41) bad = bad + 1   ! previous traj short/long
               itraj = itraj + 1
               if (tj /= cur_t + 1) bad = bad + 1
               cur_t = tj
               nfr = 1
            else
               nfr = nfr + 1
               if (sj /= cur_s + 1) bad = bad + 1
               if (tv <= t_prev) bad = bad + 1
            end if
            cur_s = sj
            t_prev = tv
            expect_col = .true.
         else if (expect_col .and. index(line, '=') > 0 .and. index(line, '#') == 0) then
            ! the column line right after a frame header
            p1 = index(line, 't_kin=')
            p2 = index(line, 'v_pot=')
            p3 = index(line, 'e_tot=')
            if (p1 == 0 .or. p2 == 0 .or. p3 == 0 .or. .not. (p1 < p2 .and. p2 < p3)) then
               bad = bad + 1
            end if
            expect_col = .false.
         end if
      end do
      close (v)
      if (cur_t > 0 .and. cur_s /= 41) bad = bad + 1    ! the last trajectory's count
      ! frames per trajectory = MAX_STEPS+1 = 41; itraj must be 3
      if (bad > 0) then
         call fail('(c): frame structure violated (monotonicity or column order)')
         return
      end if
      if (itraj /= 3 .or. nfr /= 41) then
         call fail('(c): expected 3 trajectories x 41 frames (MAX_STEPS+1, the'// &
                  ' safety-net-only run), got '//i2s(itraj)//' x last-count '//i2s(nfr))
         return
      end if
      print '(a)', '(c): PASS - traj/step/t monotone, column order t_kin<v_pot<e_tot,'// &
           ' 3 trajectories x 41 frames (MAX_STEPS+1, no termination bound)'
   end subroutine crit_c

   !------------------------------------------------------------------
   ! crit_g() - (g): the first frame's sec_qp section carries non-zero momenta
   subroutine crit_g()
      character(len=512) :: line
      integer :: v, ios, nat, k, pos
      real(8) :: vals(8), pnorm
      logical :: ok, in_sec, first_frame
      pnorm = 0.0d0
      nat = 0
      in_sec = .false.
      first_frame = .true.
      open (newunit=v, file='build/drv_case_a/trajectories.txt', status='old', &
            action='read', iostat=ios)
      if (ios /= 0) then
         call fail('(g): trajectories.txt missing')
         return
      end if
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, '# frame traj=') == 1 .and. .not. first_frame) exit
         if (index(line, '# sec_qp begin') == 1) then
            in_sec = .true.
            cycle
         end if
         if (index(line, '# sec_qp end') == 1) then
            in_sec = .false.
            first_frame = .false.
            cycle
         end if
         if (in_sec) then
            ! atom row '<index> <symb> <q(3)> <p(3)>' - skip the index and
            ! symbol tokens, then read the six reals; the last three are the
            ! momenta
            pos = 1
            call skip_tok(line, pos)
            call skip_tok(line, pos)
            do k = 1, 6
               call tok_real(line, pos, vals(k), ok)
               if (.not. ok) exit
            end do
            if (ok) then
               nat = nat + 1
               pnorm = pnorm + abs(vals(4)) + abs(vals(5)) + abs(vals(6))
            end if
         end if
      end do
      close (v)
      if (nat < 4) then
         call fail('(g): the first frame carries fewer than 4 sec_qp atom rows'// &
                  ' (got '//i2s(nat)//')')
         return
      end if
      if (pnorm <= 1.0d-6) then
         call fail('(g): the initial momenta are zero - the sampling did not happen')
         return
      end if
      print '(a)', '(g): PASS - DIST_SCHEME members ran on the real chain and the'// &
           ' incident channel assembled the pair; first-frame |p| sum = '//d2s(pnorm)
   end subroutine crit_g

   !------------------------------------------------------------------
   ! crit_d() - (d): MAX_STEPS=5, one trajectory -> exactly 6 frames
   subroutine crit_d()
      character(len=2048) :: obuf, ebuf
      integer :: stat, ntraj, nfr, last_s
      real(8) :: last_t
      call run_driver('drv_case_m', obuf, ebuf, stat)
      if (stat /= 0) then
         call fail('(d): the driver child exited nonzero - stderr: '//trim(ebuf))
         return
      end if
      call traj_summary('build/drv_case_m/trajectories.txt', ntraj, nfr, last_s, last_t)
      if (ntraj /= 1 .or. nfr /= 6) then
         call fail('(d): expected exactly 6 frames (1 + MAX_STEPS=5), got '//i2s(nfr))
         return
      end if
      if (last_s /= 6 .or. last_t <= 0.2d0) then
         call fail('(d): the last frame is not the safety-net step (step '//i2s(last_s)// &
                  ', t '//d2s(last_t)//')')
         return
      end if
      print '(a)', '(d): PASS - stopped exactly at the MAX_STEPS frame count (6),'// &
           ' last step 6 t '//d2s(last_t)
   end subroutine crit_d

   !------------------------------------------------------------------
   ! crit_e() - (e): the term-bound case terminates early (frame count below the
   !             safety-net horizon); also keeps case_t's stdout for (f)
   subroutine crit_e()
      character(len=2048) :: obuf, ebuf
      integer :: stat, ntraj, nfr, last_s
      real(8) :: last_t
      call run_driver('drv_case_t', obuf, ebuf, stat)
      obuf_t = obuf
      if (stat /= 0) then
         call fail('(e): the driver child exited nonzero - stderr: '//trim(ebuf))
         return
      end if
      call traj_summary('build/drv_case_t/trajectories.txt', ntraj, nfr, last_s, last_t)
      if (ntraj /= 2) then
         call fail('(e): expected 2 trajectories, got '//i2s(ntraj))
         return
      end if
      if (nfr < 2 .or. nfr >= 401) then
         call fail('(e): the term-bound run did not terminate early - frames '//i2s(nfr)// &
                  ' (MAX_STEPS=400, safety-net horizon 401)')
         return
      end if
      if (last_t >= 20.0d0) then
         call fail('(e): the last frame sits at the safety-net horizon (t '//d2s(last_t)//')')
         return
      end if
      print '(a)', '(e): PASS - the container termination fired mid-run: '//i2s(nfr)// &
           ' frames, last step '//i2s(last_s)//' t '//d2s(last_t)//' (horizon 401 frames)'
   end subroutine crit_e

   !------------------------------------------------------------------
   ! crit_f() - (f): the unbound-classification skipped line + the drift line in
   !             case_a's stdout; the bound outcome table in case_t's stdout
   subroutine crit_f()
      if (index(obuf_a, 'outcome counts: none (no classification bound - skipped)') == 0) then
         call fail('(f): the no-classification run must print the skipped-semantics'// &
                  ' outcome line - got "'//trim(obuf_a)//'"')
         return
      end if
      if (index(obuf_a, 'conservation drift: trajectories=3 max|e_fin-e0|=') == 0) then
         call fail('(f): the conservation-drift summary line missing - got "'// &
                  trim(obuf_a)//'"')
         return
      end if
      if (index(obuf_t, 'outcome counts: scattering=2') == 0) then
         call fail('(f): the classification-bound run must print the outcome table'// &
                  ' (scattering=2 on the real final-state values) - got "'// &
                  trim(obuf_t)//'"')
         return
      end if
      print '(a)', '(f): PASS - unbound run: skipped line + drift summary; bound run:'// &
           ' outcome table scattering=2 (the real fin values: the high-energy'// &
           ' pair stays above the split)'
   end subroutine crit_f

   !------------------------------------------------------------------
   ! crit_det() - (det): the same input twice -> trajectories.txt byte-identical
   subroutine crit_det()
      character(len=2048) :: obuf, ebuf
      character(len=512) :: cmd
      character(len=1) :: sep
      integer :: stat, statc
      logical :: ok1
      call run_driver('drv_case_a', obuf, ebuf, stat)
      if (stat /= 0) then
         call fail('(det): the first rerun exited nonzero - stderr: '//trim(ebuf))
         return
      end if
      if (windows_host()) then
         sep = '\'
         call execute_command_line('copy /y "build\drv_case_a\trajectories.txt"'// &
              ' "build\drv_case_a\trajectories_ref.txt" >nul 2>nul', wait=.true., &
              exitstat=statc)
      else
         sep = '/'
         call execute_command_line('cp -f build/drv_case_a/trajectories.txt'// &
              ' build/drv_case_a/trajectories_ref.txt', wait=.true., exitstat=statc)
      end if
      ok1 = (statc == 0)
      call run_driver('drv_case_a', obuf, ebuf, stat)
      if (stat /= 0) then
         call fail('(det): the second rerun exited nonzero - stderr: '//trim(ebuf))
         return
      end if
      if (windows_host()) then
         cmd = 'fc /b build'//sep//'drv_case_a'//sep//'trajectories.txt build'//sep// &
               'drv_case_a'//sep//'trajectories_ref.txt >nul 2>nul'
      else
         cmd = 'cmp -s build/drv_case_a/trajectories.txt build/drv_case_a/trajectories_ref.txt'
      end if
      call execute_command_line(trim(cmd), wait=.true., exitstat=statc)
      if (.not. ok1 .or. statc /= 0) then
         call fail('(det): the two runs differ (or the compare failed) - exitstat '// &
                  i2s(statc))
         return
      end if
      print '(a)', '(det): PASS - same input twice, trajectories.txt byte-identical'// &
           ' (rng single stream, deterministic end to end)'
   end subroutine crit_det

   !------------------------------------------------------------------
   ! crit_e7() - (E7): the source-level integration-loop call-set assertion
   !             (method per the header: call targets + the container_term
   !             occurrence, blank/comment lines skipped)
   subroutine crit_e7()
      character(len=512) :: line, ls
      character(len=64) :: names(16), expect(5)
      integer :: nn, v, ios, p, q, i, ltr
      logical :: inloop, closed
      expect = [ character(len=64) :: 'prop_step', 'elec_prop', 'mqc_statistic', &
                                       'rec_frame', 'container_term' ]
      nn = 0
      inloop = .false.
      closed = .false.
      open (newunit=v, file='driver.f90', status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call fail('(E7): src/driver.f90 not readable from the run directory')
         return
      end if
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         ls = adjustl(line)
         if (len_trim(ls) == 0 .or. ls(1:1) == '!') cycle
         if (.not. inloop) then
            if (index(ls, 'do while') /= 0 .and. index(ls, 'term_flag') /= 0) then
               inloop = .true.
            end if
            cycle
         end if
         if (index(ls, 'end do') == 1) then
            closed = .true.
            exit
         end if
         ltr = len_trim(ls)
         p = 1
         do
            q = index(ls(p:), 'call ')
            if (q == 0) exit
            p = p + q + 4
            do while (p <= ltr .and. ls(p:p) == ' ')
               p = p + 1
            end do
            q = p
            do while (q <= ltr .and. ls(q:q) /= '(' .and. ls(q:q) /= ' ')
               q = q + 1
            end do
            nn = nn + 1
            if (nn > size(names)) then
               call fail('(E7): more call sites than capacity inside the loop')
               return
            end if
            names(nn) = ls(p:q-1)
         end do
         if (index(ls, 'container_term') /= 0) then
            nn = nn + 1
            names(nn) = 'container_term'
         end if
      end do
      close (v)
      if (.not. inloop .or. .not. closed) then
         call fail('(E7): the integration loop (do-while naming term_flag .. end do)'// &
                  ' was not found in src/driver.f90')
         return
      end if
      if (nn /= 5) then
         call fail('(E7): the integration loop carries '//i2s(nn)// &
                  ' direct entries, expected exactly 5')
         return
      end if
      do i = 1, 5
         if (count(names(1:nn) == expect(i)) /= 1) then
            call fail('(E7): the entry '//trim(expect(i))//' does not appear exactly'// &
                     ' once in the integration loop')
            return
         end if
      end do
      print '(a)', '(E7): PASS - the integration loop calls exactly the five'// &
           ' legislated entries: prop_step, elec_prop, mqc_statistic, rec_frame,'// &
           ' container_term'
   end subroutine crit_e7
   !------------------------------------------------------------------
   ! crit_h() - (h): bad SYSTEM_DIR -> the sysdef channel's own named abort penetrates
   subroutine crit_h()
      character(len=2048) :: obuf, ebuf
      integer :: stat
      call run_driver('drv_case_b', obuf, ebuf, stat)
      if (stat == 0) then
         call fail('(h): expected a loud-stop, the driver child exited 0')
      else if (index(ebuf, 'sysdef_load') == 0) then
         call fail('(h): stderr does not name sysdef_load - got "'//trim(ebuf)//'"')
      else
         print '(a)', '(h): PASS ('//trim(ebuf)//')'
      end if
   end subroutine crit_h

   !------------------------------------------------------------------
   ! crit_j1() - (j1): the in-process seam round-trip (the legislated order) and
   !              the unbound container_term posture (no container keys here)
   subroutine crit_j1()
      type(state_t) :: st
      call registry_gen_keys()
      call read_input('build/drv_case_d/input_dpos.txt')
      call rng_init(i_seed(1))
      call sysdef_load(reactants%system_dir)
      call list_atoms_load()
      call state_create(st, size(list_atoms%mass), 1)   ! the degenerate adiabatic surface count
      call registry_gen_all()
      call registry_gen_init()
      call input_audit_declared()
      call container_init()
      call rec_select(observables%rec_level)
      if (len_trim(observables%obs_list) > 0) call rec_select(trim(observables%obs_list))
      if (glo_target_t_glo() /= 350.0d0 .or. glo_target_gamma() /= 0.04d0 .or. &
          glo_target_w() /= 1.25d0) then
         call fail('(j1): member getter round-trip mismatch - t_glo/gamma/w read '// &
                  d2s(glo_target_t_glo())//' '//d2s(glo_target_gamma())//' '// &
                  d2s(glo_target_w()))
      else
         print '(a)', '(j1): PASS - buffered member keys round-trip: t_glo=350 gamma=0.04 w=1.25'
      end if
      if (glo_target_dt_loaded() /= 0.05d0 .or. glo_target_dt_loaded() /= prop_cfg%dt) then
         call fail('(j1): dt-ordering mismatch - glo_target_dt_loaded() = '// &
                  d2s(glo_target_dt_loaded())//' but DT buffered '//d2s(prop_cfg%dt))
      else
         print '(a)', '(j1): PASS - the member init snapshotted the finalized dt (0.05)'
      end if
      if (samp_n_dist() /= 11) then
         call fail('(j1): member-table count expected 11 (10 scheme + one paradigm-armed incident member - thermostat left the family 2026-09-22), got '// &
                  i2s(samp_n_dist()))
      else if (rec_n_sel() /= 3) then
         call fail('(j1): selected column count expected 3 (the framework energy '// &
                  'columns at the default level), got '//i2s(rec_n_sel()))
      else
         print '(a)', '(j1): PASS - member table carries 12 rows; selected columns 3'
      end if
      if (container_term(st)) then
         call fail('(j1): unbound container_term must return .false.')
      else
         print '(a)', '(j1): PASS - container_term on the unbound slot returns .false.'
      end if
   end subroutine crit_j1

   !------------------------------------------------------------------
   ! crit_j2() - (j2): T_GLO = -5 through the real driver exe
   subroutine crit_j2()
      character(len=2048) :: obuf, ebuf
      integer :: stat
      call run_driver('drv_case_d', obuf, ebuf, stat)
      if (stat == 0) then
         call fail('(j2): expected a loud-stop, the driver child exited 0')
      else if (index(ebuf, 'glo_target_init') == 0) then
         call fail('(j2): stderr does not name glo_target_init - got "'//trim(ebuf)//'"')
      else if (index(ebuf, 't_glo') == 0) then
         call fail('(j2): stderr missing the field name t_glo - got "'//trim(ebuf)//'"')
      else
         print '(a)', '(j2): PASS ('//trim(ebuf)//')'
      end if
   end subroutine crit_j2

   !------------------------------------------------------------------
   ! crit_k() - (k): this check's own child runs the unbound-classification chain
   subroutine crit_k()
      character(len=512) :: cmd
      character(len=2048) :: ebuf
      character(len=1) :: sep
      integer :: stat
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_driver.exe eclass 2>build'//sep//'drv_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp('build/drv_err.txt', ebuf)
      if (stat == 0) then
         call fail('(k): expected a loud-stop, the eclass child exited 0')
      else if (index(ebuf, 'container_classify') == 0) then
         call fail('(k): stderr does not name container_classify - got "'//trim(ebuf)//'"')
      else
         print '(a)', '(k): PASS ('//trim(ebuf)//')'
      end if
   end subroutine crit_k

   !------------------------------------------------------------------
   ! crit_l() - (l): the two histogram files of the term+classify case -
   !             existence, 40+1 line shape, constant spacing, count truth
   subroutine crit_l()
      real(8) :: cen(50), cnt(50), tot, pp, dlt
      integer :: nbin, nv, k, bad
      ! ---- hist_e_rel.dat ----
      call read_hist('build/drv_case_t/hist_e_rel.dat', 50, cen, cnt, nbin, tot, pp, nv)
      if (nv == 0) then
         call fail('(l): hist_e_rel.dat was not created in the case directory')
         return
      end if
      if (nv /= 41 .or. nbin /= 40) then
         call fail('(l): hist_e_rel.dat expected 40 bin lines + 1 TOTAL line, got '// &
                  i2s(nv)//' lines / '//i2s(nbin)//' bins')
         return
      end if
      bad = 0
      dlt = cen(2) - cen(1)
      if (abs(cen(1) - dlt/2.0d0) > 1.0d-4) bad = bad + 1
      do k = 2, 39
         if (abs((cen(k+1) - cen(k)) - dlt) > 1.0d-4) bad = bad + 1
      end do
      if (abs(sum(cnt(1:40)) - 2.0d0) > 1.0d-6) bad = bad + 1
      if (abs(tot - 2.0d0) > 1.0d-6 .or. abs(pp - 1.0d0) > 1.0d-6) bad = bad + 1
      if (bad > 0) then
         call fail('(l): hist_e_rel.dat violates the shape/count truth (spacing,'// &
                  ' first center, TOTAL/P)')
         return
      end if
      ! ---- hist_thta.dat (the fully derived truth: bitwise-0 angles) ----
      call read_hist('build/drv_case_t/hist_thta.dat', 50, cen, cnt, nbin, tot, pp, nv)
      if (nv == 0) then
         call fail('(l): hist_thta.dat was not created in the case directory')
         return
      end if
      if (nv /= 41 .or. nbin /= 40) then
         call fail('(l): hist_thta.dat expected 40 bin lines + 1 TOTAL line, got '// &
                  i2s(nv)//' lines / '//i2s(nbin)//' bins')
         return
      end if
      bad = 0
      dlt = cen(2) - cen(1)
      if (abs(dlt - 0.025d0) > 1.0d-4 .or. abs(cen(1) - 0.0125d0) > 1.0d-4) bad = bad + 1
      if (abs(cnt(1) - 2.0d0) > 1.0d-6) bad = bad + 1
      do k = 2, 40
         if (abs(cnt(k)) > 1.0d-6) bad = bad + 1
      end do
      if (abs(tot - 2.0d0) > 1.0d-6 .or. abs(pp - 1.0d0) > 1.0d-6) bad = bad + 1
      if (bad > 0) then
         call fail('(l): hist_thta.dat violates the derived truth (bin-1 count 2,'// &
                  ' spacing 0.025, TOTAL/P)')
         return
      end if
      print '(a)', '(l): PASS - both histogram files exist with the 40+TOTAL shape;'// &
           ' e_rel counts both rows (P=1); theta bin 1 carries all rows (gas-phase'// &
           ' angle bitwise 0, spacing exactly 0.025)'
   end subroutine crit_l

   !------------------------------------------------------------------
   ! crit_m() - (m): the histogram summary line + the printed drift vs the
   !             file-derived drift
   subroutine crit_m()
      integer, parameter :: max_traj = 8
      real(8) :: e_first(max_traj), dev_max(max_traj), efin(max_traj)
      integer :: ntr(max_traj), cur, nfin, p, ios, t
      real(8) :: printed, derived
      if (index(obuf_t, 'histograms: rows=2 files=hist_e_rel.dat hist_thta.dat') == 0) then
         call fail('(m): stdout missing the histogram summary line - got "'// &
                  trim(obuf_t)//'"')
         return
      end if
      p = index(obuf_t, 'max|e_fin-e0|=')
      if (p == 0) then
         call fail('(m): stdout missing the drift line - got "'//trim(obuf_t)//'"')
         return
      end if
      read (obuf_t(p+14:), *, iostat=ios) printed
      if (ios /= 0) then
         call fail('(m): the printed drift value is unparseable - got "'// &
                  trim(obuf_t(p+14:p+25))//'"')
         return
      end if
      call scan_columns('build/drv_case_t/trajectories.txt', max_traj, ntr, e_first, &
                        dev_max, cur)
      call scan_finals('build/drv_case_t/trajectories.txt', max_traj, efin, nfin)
      if (cur < 2 .or. nfin < 2 .or. any(ntr(1:2) == 0)) then
         call fail('(m): could not parse first frames and final rows for both'// &
                  ' trajectories')
         return
      end if
      derived = 0.0d0
      do t = 1, cur
         derived = max(derived, abs(efin(t) - e_first(t)))
      end do
      if (abs(printed - derived) > 1.0d-3 + 1.0d-3*abs(derived)) then
         call fail('(m): the printed drift '//d2s(printed)//' disagrees with the'// &
                  ' file-derived '//d2s(derived)//' kcal/mol')
         return
      end if
      print '(a)', '(m): PASS - the histograms summary line present; printed drift '// &
           d2s(printed)//' = file-derived max|final e_tot - first e_tot| '//d2s(derived)
   end subroutine crit_m

   !------------------------------------------------------------------
   ! crit_n() - (n): occ/n_hop through the real driver, every frame line and
   !             final row (the degenerate adiabatic values end to end).
   !             The case input omits GWRITE_LEVEL, so under the bare-input default
   !             (level 2) every frame also carries a sec_qp section - the
   !             scan skips the section bodies (atom rows share the final
   !             row's no-# no-= shape but not its column count/semantics)
   subroutine crit_n()
      character(len=2048) :: obuf, ebuf
      character(len=512) :: line
      integer :: stat, v, ios, nfr, nfin, bad, tj
      real(8) :: tt, w1, w2, w3, vo, vn
      logical :: want_col, in_sec
      call run_driver('drv_case_o', obuf, ebuf, stat)
      if (stat /= 0) then
         call fail('(n): the driver child exited nonzero - stderr: '//trim(ebuf))
         return
      end if
      if (index(obuf, 'columns=5') == 0) then
         call fail('(n): the assembly line must read columns=5 (energy three + the'// &
                  ' two occupation columns) - got "'//trim(obuf)//'"')
         return
      end if
      nfr = 0
      nfin = 0
      bad = 0
      want_col = .false.
      in_sec = .false.
      open (newunit=v, file='build/drv_case_o/trajectories.txt', status='old', &
            action='read', iostat=ios)
      if (ios /= 0) then
         call fail('(n): trajectories.txt missing in the case directory')
         return
      end if
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, '# sec_qp begin') == 1) then
            in_sec = .true.
         else if (index(line, '# sec_qp end') == 1) then
            in_sec = .false.
         else if (in_sec) then
            cycle
         else if (want_col) then
            ! the frame's column line: occ/n_hop by name token
            want_col = .false.
            nfr = nfr + 1
            if (tok_val(line, 'occ', vo) .and. tok_val(line, 'n_hop', vn)) then
               if (abs(vo - 1.0d0) > 1.0d-9 .or. abs(vn) > 1.0d-9) bad = bad + 1
            else
               bad = bad + 1
            end if
         else if (index(line, '# frame traj=') == 1) then
            want_col = .true.
         else if (len_trim(line) > 0 .and. index(line, '#') == 0 .and. &
                  index(line, '=') == 0) then
            ! a final row: traj t t_kin v_pot e_tot occ n_hop (positional)
            nfin = nfin + 1
            read (line, *, iostat=ios) tj, tt, w1, w2, w3, vo, vn
            if (ios /= 0 .or. abs(vo - 1.0d0) > 1.0d-9 .or. abs(vn) > 1.0d-9) then
               bad = bad + 1
            end if
         end if
      end do
      close (v)
      if (nfr /= 12 .or. nfin /= 2) then
         call fail('(n): expected 12 frame column lines + 2 final rows, got '// &
                  i2s(nfr)//' + '//i2s(nfin))
         return
      end if
      if (bad > 0) then
         call fail('(n): '//i2s(bad)//' lines carry wrong occ/n_hop values (want'// &
                  ' 1.0 / 0.0 everywhere - the adiabatic degeneracy)')
         return
      end if
      print '(a)', '(n): PASS - occ=1.0 and n_hop=0.0 on all 12 frame lines and both'// &
           ' final rows: input key -> selection -> extractor -> getter -> file'
   end subroutine crit_n

   !------------------------------------------------------------------
   ! crit_phg() - (phg): this check's own child drives the inverted range
   subroutine crit_phg()
      character(len=512) :: cmd
      character(len=2048) :: ebuf
      character(len=1) :: sep
      integer :: stat
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_driver.exe phguard 2>build'//sep//'drv_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp('build/drv_err.txt', ebuf)
      if (stat == 0) then
         call fail('(phg): expected a loud-stop, the phguard child exited 0')
      else if (index(ebuf, 'plot_hist') == 0) then
         call fail('(phg): stderr does not name plot_hist - got "'//trim(ebuf)//'"')
      else
         print '(a)', '(phg): PASS ('//trim(ebuf)//')'
      end if
   end subroutine crit_phg

   !------------------------------------------------------------------
   ! crit_q() - (q): the non-linear fragment anchor child
   subroutine crit_q()
      character(len=512) :: cmd
      character(len=2048) :: obuf, ebuf
      character(len=1) :: sep
      integer :: stat
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_driver.exe nlfrag >build'//sep//'drv_qout.txt'// &
            ' 2>build'//sep//'drv_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp('build/drv_qout.txt', obuf)
      call slurp('build/drv_err.txt', ebuf)
      if (stat /= 0) then
         call fail('(q): the nlfrag child exited nonzero - stdout: '//trim(obuf)// &
                  ' stderr: '//trim(ebuf))
      else if (index(obuf, 'nlfrag: non-linear fragment anchors hold') == 0) then
         call fail('(q): the child output is missing the anchor marker - '//trim(obuf))
      else
         print '(a)', '(q): PASS - the 3-atom non-collinear fragment: rotation'// &
              ' stripped from the archive (p_fin ~ 0), j/e_rot/internal-kinetic'// &
              ' cross-identity, buffered well, COM-shifted rows'
      end if
   end subroutine crit_q

   !------------------------------------------------------------------
   ! run_driver(cdir, obuf, ebuf, stat) - spawn the driver exe with cwd = one case
   !                  directory (the driver reads ./input_qct.txt relative to its
   !                  cwd; 'cd <case> && exe' carries both shell families); stdout
   !                  and stderr land in the case dir and are slurped back
   subroutine run_driver(cdir, obuf, ebuf, stat)
      character(len=*), intent(in) :: cdir
      character(len=*), intent(out) :: obuf, ebuf
      integer, intent(out) :: stat
      character(len=1) :: sep
      character(len=512) :: cmd
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'cd build'//sep//trim(cdir)//' && ..'//sep//'drv_venus.exe >out.txt 2>err.txt'
      obuf = ''
      ebuf = ''
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp('build/'//trim(cdir)//'/out.txt', obuf)
      call slurp('build/'//trim(cdir)//'/err.txt', ebuf)
   end subroutine run_driver

   !------------------------------------------------------------------
   ! make_cases() - write every case folder under build/ (idempotent; the driver
   !                  exe copy gets a .exe extension because Windows cmd cannot
   !                  execute the build's venus.e name directly; BOTH compiler
   !                  trees are probed and a missing exe is a NAMED failure)
   subroutine make_cases()
      character(len=80) :: xyza(4), xyzb(4), xyz1(4), xyz3(5)
      character(len=80) :: la(12), lb(6), ld(6), ldp(6), lec(9)
      character(len=80) :: lm(9), lt(14)
      character(len=80) :: lr(13), lo(10), ln(7)
      logical :: have_g, have_i
      character(len=64) :: srcexe, srcw
      integer :: j
      xyza = [ character(len=80) :: '2', 'generic diatomic fragment a', &
               'H 0.0 0.0 0.0', 'H 0.7 0.0 0.0' ]
      xyzb = [ character(len=80) :: '2', 'generic diatomic fragment b', &
               'H 0.0 0.0 0.0', 'H 0.8 0.0 0.0' ]
      xyz1 = [ character(len=80) :: '2', 'diatomic test fragment', &
               'H 0.0 0.0 0.0', 'H 0.7 0.0 0.0' ]
      ! the non-collinear triangle for the (q) anchor: centered (COM at the
      ! origin for equal masses) and equilateral in the z = 0 plane, so the
      ! inertia tensor is diagonal (I_xy cancels pairwise) and a z-axis
      ! rotation exercises the full inertia-tensor / rotation-strip path
      xyz3 = [ character(len=80) :: '3', 'non-collinear triangle fragment', &
               'H 0.9 0.0 0.0', 'H -0.45 0.7794228634059948 0.0', &
               'H -0.45 -0.7794228634059948 0.0' ]
      ! (a): the end-to-end case - two-fragment gas-phase system, thermalize on
      ! fragment 1, rotation on fragment 2, the incident channel for the pair;
      ! GWRITE_LEVEL=2 opens the sec_qp flattened section (the decision-fixed
      ! bare-input default - the momenta criterion reads the
      ! section); NO container parameter file -> termination and
      ! classification unbound
      la = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', &
             'DIST_SCHEME = THERMALIZE[T_EQ=300.0, N_EQ=25, DT_EQ=0.05], ROTATION[J_ROT=3]', &
             'PROJECTILE = pair_a', 'TARGET = pair_b', 'N_TRAJ = 3', &
             'MAX_STEPS = 40', 'ISEED = 20260915', 'GWRITE_LEVEL = 2' ]
      ! (h): the bad system_dir case (the folder does not exist)
      lb = [ character(len=80) :: 'SYSTEM_DIR = no_such_system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION' ]
      ! (j2): the bad member value (driver child, cwd-relative system dir)
      ld = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', &
             'DIST_SCHEME = GLO_TARGET[T_GLO=-5.0, GAMMA_GLO=0.04, W_GHOST=1.25]' ]
      ! (j1): the positive round-trip (the check parent reads it from src/, so the
      ! system dir is build-relative)
      ldp = [ character(len=80) :: 'SYSTEM_DIR = build/drv_case_d/system', 'E_REL = 1.0', &
              'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.05', &
              'DIST_SCHEME = GLO_TARGET[T_GLO=350.0, GAMMA_GLO=0.04, W_GHOST=1.25]' ]
      ! (k): the eclass child input (build-relative; the child runs from src/,
      ! whose run directory carries no container parameter file; one member
      ! per fragment - a member instance carries one parameter set)
      lec = [ character(len=80) :: 'SYSTEM_DIR = build/drv_case_a/system', 'E_REL = 1.0', &
              'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.05', &
              'DIST_SCHEME = THERMALIZE[T_EQ=300.0, N_EQ=25, DT_EQ=0.05], ROTATION', &
              'PROJECTILE = pair_a', 'TARGET = pair_b', 'N_TRAJ = 1' ]
      ! (d): the MAX_STEPS safety-net case - one trajectory, five steps
      lm = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION', 'N_TRAJ = 1', &
             'MAX_STEPS = 5', 'ISEED = 20260915' ]
      ! (r): the radau end-to-end case - the (a) material under the adaptive
      ! member (ORDER=8 -> ss = 1e-8, the unit-level anchor's accuracy;
      ! MAX_STEPS=40 keeps the safety-net-only posture, 2 trajectories keep
      ! the run light at ~44 force calls per boot sequence)
      lr = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'INTEGRATOR = RADAU', 'ORDER = 8', &
             'DIST_SCHEME = THERMALIZE[T_EQ=300.0, N_EQ=25, DT_EQ=0.05], ROTATION[J_ROT=3]', &
             'PROJECTILE = pair_a', 'TARGET = pair_b', 'N_TRAJ = 2', 'MAX_STEPS = 40', &
             'ISEED = 20260915' ]
      ! (e)/(f): the term+classify bound case - a high collision energy so the
      ! amplitude climb crosses Q_TERM mid-run; MAX_STEPS=400 is the safety net
      lt = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 200.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', &
             'DIST_SCHEME = THERMALIZE[T_EQ=300.0, N_EQ=25, DT_EQ=0.05], ROTATION[J_ROT=3]', &
             'PROJECTILE = pair_a', 'TARGET = pair_b', 'N_TRAJ = 2', &
             'MAX_STEPS = 400', 'ISEED = 20260915', &
             'K_HARM = 0.04', 'Q_TERM = 9.5', 'E_SPLIT = 1.0' ]   ! container keys
      ! (n): the occupation-column case - the list entry appends occ/n_hop to
      ! the energy three (columns=5); a short safety-net run keeps the file
      ! small while exercising every frame line and final row
      lo = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION', 'N_TRAJ = 2', &
             'MAX_STEPS = 5', 'ISEED = 20260915', 'OBSERVABLES = occ,n_hop' ]
      ! (q): the non-linear fragment case - the child reads it from src/ with a
      ! build-relative system dir; DIST_SCHEME satisfies the required-key
      ! grammar but the child never samples (it binds the analysis chain itself)
      ln = [ character(len=80) :: 'SYSTEM_DIR = build/drv_case_n/system', 'E_REL = 1.0', &
             'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION', &
             'N_TRAJ = 1' ]
      ! wipe every case folder BEFORE rebuilding it: a case folder must hold
      ! exactly what THIS run writes - a stale file from an earlier round
      ! (a renamed system sample, a relocated parameter file) would ride the
      ! system-folder scan as an extra reactant and corrupt the case
      call wipe_dir('build/drv_case_a')
      call wipe_dir('build/drv_case_b')
      call wipe_dir('build/drv_case_d')
      call wipe_dir('build/drv_case_m')
      call wipe_dir('build/drv_case_r')
      call wipe_dir('build/drv_case_t')
      call wipe_dir('build/drv_case_o')
      call wipe_dir('build/drv_case_n')
      call make_dir('build/drv_case_a/system')
      call make_dir('build/drv_case_b')
      call make_dir('build/drv_case_d/system')
      call make_dir('build/drv_case_m/system')
      call make_dir('build/drv_case_r/system')
      call make_dir('build/drv_case_t/system')
      call make_dir('build/drv_case_o/system')
      call make_dir('build/drv_case_n/system')
      call write_sample('build/drv_case_a/system/pair_a.xyz', xyza)
      call write_sample('build/drv_case_a/system/pair_b.xyz', xyzb)
      call write_sample('build/drv_case_d/system/h2.xyz', xyz1)
      call write_sample('build/drv_case_m/system/h2.xyz', xyz1)
      call write_sample('build/drv_case_r/system/pair_a.xyz', xyza)
      call write_sample('build/drv_case_r/system/pair_b.xyz', xyzb)
      call write_sample('build/drv_case_t/system/pair_a.xyz', xyza)
      call write_sample('build/drv_case_t/system/pair_b.xyz', xyzb)
      call write_sample('build/drv_case_o/system/h2.xyz', xyz1)
      call write_sample('build/drv_case_n/system/h3.xyz', xyz3)
      call write_sample('build/drv_case_a/input_qct.txt', la)
      call write_sample('build/drv_case_a/input_eclass.txt', lec)
      call write_sample('build/drv_case_b/input_qct.txt', lb)
      call write_sample('build/drv_case_d/input_qct.txt', ld)
      call write_sample('build/drv_case_d/input_dpos.txt', ldp)
      call write_sample('build/drv_case_m/input_qct.txt', lm)
      call write_sample('build/drv_case_r/input_qct.txt', lr)
      call write_sample('build/drv_case_t/input_qct.txt', lt)
      call write_sample('build/drv_case_o/input_qct.txt', lo)
      call write_sample('build/drv_case_n/input_nlfrag.txt', ln)
      ! the driver exe: probe BOTH compiler trees, name the failure when neither
      ! exists (a silent copy failure is not diagnosable)
      inquire (file='build/gfortran/venus.e', exist=have_g)
      inquire (file='build/ifx/venus.e', exist=have_i)
      srcexe = ''
      if (have_g) then
         srcexe = 'build/gfortran/venus.e'
      else if (have_i) then
         srcexe = 'build/ifx/venus.e'
      else
         call fail('setup: no driver executable - run make first (looked for'// &
                   ' build/gfortran/venus.e and build/ifx/venus.e)')
         return
      end if
      if (windows_host()) then
         srcw = ''
         do j = 1, len_trim(srcexe)
            if (srcexe(j:j) == '/') then
               srcw(j:j) = '\'
            else
               srcw(j:j) = srcexe(j:j)
            end if
         end do
         call execute_command_line('copy /y "'//trim(srcw)//'" "build\drv_venus.exe"'// &
                                   ' >nul 2>nul', wait=.true.)
      else
         call execute_command_line('cp -f '//trim(srcexe)//' build/drv_venus.exe', &
                                   wait=.true.)
      end if
   end subroutine make_cases

   !------------------------------------------------------------------
   ! parsing helpers (the trajectories.txt readers the frame criteria share)
   !------------------------------------------------------------------

   ! parse_frame(line, tj, sj, tv) - read traj/step/t out of one frame header
   subroutine parse_frame(line, tj, sj, tv)
      character(len=*), intent(in) :: line
      integer, intent(out) :: tj, sj
      real(8), intent(out) :: tv
      integer :: p, ios
      tj = 0
      sj = 0
      tv = 0.0d0
      p = index(line, 'traj=')
      if (p > 0) read (line(p+5:), *, iostat=ios) tj
      p = index(line, 'step=')
      if (p > 0) read (line(p+5:), *, iostat=ios) sj
      p = index(line, ' t=')
      if (p > 0) read (line(p+3:), *, iostat=ios) tv
   end subroutine parse_frame

   ! col_val(line, name, v) - read the value of one '<name>=' token in a column line
   subroutine col_val(line, name, v)
      character(len=*), intent(in) :: line, name
      real(8), intent(out) :: v
      integer :: p, ios
      v = 0.0d0
      p = index(line, name//'=')
      if (p > 0) read (line(p+len_trim(name)+1:), *, iostat=ios) v
   end subroutine col_val

   ! scan_columns(fname, ...) - one pass: per trajectory, the first frame's e_tot and
   ! the max deviation of every later frame's e_tot from it
   subroutine scan_columns(fname, maxt, ntr, efirst, devmax, nseen)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: maxt
      integer, intent(out) :: ntr(maxt), nseen
      real(8), intent(out) :: efirst(maxt), devmax(maxt)
      character(len=512) :: line
      integer :: v, ios, tj, sj
      real(8) :: tv, ev
      logical :: expect_col
      nseen = 0
      ntr = 0
      efirst = 0.0d0
      devmax = 0.0d0
      tj = 0
      expect_col = .false.
      open (newunit=v, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, '# frame traj=') == 1) then
            call parse_frame(line, tj, sj, tv)
            if (tj /= nseen) then
               nseen = tj
               if (tj >= 1 .and. tj <= maxt) ntr(tj) = 0
            end if
            if (tj >= 1 .and. tj <= maxt) ntr(tj) = ntr(tj) + 1
            expect_col = .true.
         else if (expect_col .and. index(line, 'e_tot=') > 0 .and. &
                  index(line, '#') == 0) then
            call col_val(line, 'e_tot', ev)
            if (tj >= 1 .and. tj <= maxt) then
               if (ntr(tj) == 1) then
                  efirst(tj) = ev
               else
                  devmax(tj) = max(devmax(tj), abs(ev - efirst(tj)))
               end if
            end if
            expect_col = .false.
         end if
      end do
      close (v)
   end subroutine scan_columns

   ! traj_summary(fname, ...) - trajectories seen, the max per-trajectory frame count,
   ! and the final frame's step and time in the whole file
   subroutine traj_summary(fname, ntraj, nfr, last_s, last_t)
      character(len=*), intent(in) :: fname
      integer, intent(out) :: ntraj, nfr, last_s
      real(8), intent(out) :: last_t
      character(len=512) :: line
      integer :: v, ios, tj, sj, cur_t, cnt
      real(8) :: tv
      ntraj = 0
      nfr = 0
      last_s = 0
      last_t = 0.0d0
      cur_t = 0
      cnt = 0
      open (newunit=v, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, '# frame traj=') == 1) then
            call parse_frame(line, tj, sj, tv)
            if (tj /= cur_t) then
               cur_t = tj
               ntraj = ntraj + 1
               cnt = 1
            else
               cnt = cnt + 1
            end if
            nfr = max(nfr, cnt)
            last_s = sj
            last_t = tv
         end if
      end do
      close (v)
   end subroutine traj_summary

   ! read_hist(fname, ...) - one plain histogram file back: bin centers/counts
   ! (list-directed - format-agnostic) and the TOTAL & P line; a missing file
   ! reports nlines = 0 (the caller's existence check)
   subroutine read_hist(fname, cap, cen, cnt, nbin, tot, pp, nlines)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: cap
      real(8), intent(out) :: cen(cap), cnt(cap), tot, pp
      integer, intent(out) :: nbin, nlines
      character(len=256) :: line
      integer :: v, ios, ios2
      nbin = 0
      nlines = 0
      tot = -1.0d0
      pp = -1.0d0
      cen = 0.0d0
      cnt = 0.0d0
      open (newunit=v, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         nlines = nlines + 1
         if (index(line, 'TOTAL & P') == 1) then
            read (line(10:), *, iostat=ios2) tot, pp
            if (ios2 /= 0) then
               tot = -1.0d0
               pp = -1.0d0
            end if
         else if (nbin < cap) then
            nbin = nbin + 1
            read (line, *, iostat=ios2) cen(nbin), cnt(nbin)
            if (ios2 /= 0) then
               nbin = nbin - 1
               cen(nbin+1) = 0.0d0
               cnt(nbin+1) = 0.0d0
            end if
         end if
      end do
      close (v)
   end subroutine read_hist

   ! scan_finals(fname, ...) - per trajectory the e_tot of its final row (rows
   ! are the plain data lines: traj t <selected values>; e_tot is the third
   ! value under the default energy-three selection; frame column lines carry
   ! '=' and are skipped)
   subroutine scan_finals(fname, maxt, efin, nfin)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: maxt
      real(8), intent(out) :: efin(maxt)
      integer, intent(out) :: nfin
      character(len=512) :: line
      integer :: v, ios, tj
      real(8) :: tt, v1, v2, v3
      efin = 0.0d0
      nfin = 0
      open (newunit=v, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (v, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         if (index(line, '#') == 1) cycle
         if (index(line, '=') > 0) cycle
         read (line, *, iostat=ios) tj, tt, v1, v2, v3
         if (ios /= 0) cycle
         if (tj >= 1 .and. tj <= maxt) then
            efin(tj) = v3
            nfin = max(nfin, tj)
         end if
      end do
      close (v)
   end subroutine scan_finals

   ! tok_val(line, name, v) - read one '<name>=' token's value out of a frame
   ! column line (token-boundary guarded - a name must start after a blank or
   ! at the line head)
   logical function tok_val(line, name, v)
      character(len=*), intent(in) :: line, name
      real(8), intent(out) :: v
      integer :: p, ios
      tok_val = .false.
      v = 0.0d0
      p = index(line, trim(name)//'=')
      if (p == 0) return
      if (p > 1 .and. line(p-1:p-1) /= ' ') return
      read (line(p+len_trim(name)+1:), *, iostat=ios) v
      tok_val = (ios == 0)
   end function tok_val

   ! chk_nl(label, got, want) - relative anchor check, tol 1e-10 with floor 1.0
   ! (the final-state check's discipline: fixed short chains of double
   ! operations; a wrong formula deviates at the 1e-2..1e0 level)
   subroutine chk_nl(label, got, want)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got, want
      if (abs(got - want) > 1.0d-10*max(abs(want), 1.0d0)) then
         call fail(trim(label)//': got '//d2s(got)//', want '//d2s(want))
      end if
   end subroutine chk_nl

   ! chk_z(label, got) - bitwise-zero check (skip-by-construction contract)
   subroutine chk_z(label, got)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got
      if (got /= 0.0d0) then
         call fail(trim(label)//': must be bitwise zero, got '//d2s(got))
      end if
   end subroutine chk_z

   ! dev_str(d, n) - render a per-trajectory deviation list (the fail message)
   function dev_str(d, n) result(s)
      real(8), intent(in) :: d(:)
      integer, intent(in) :: n
      character(len=96) :: s
      character(len=16) :: tmp
      integer :: j
      s = ''
      do j = 1, n
         write (tmp, '(es12.4)') d(j)
         s = trim(s)//' '//trim(adjustl(tmp))
      end do
   end function dev_str

   ! skip_tok(line, pos) - advance past one whitespace-separated token (the
   ! sec_qp row's index and symbol fields are not reals)
   subroutine skip_tok(line, pos)
      character(len=*), intent(in) :: line
      integer, intent(inout) :: pos
      do while (pos <= len(line))
         if (line(pos:pos) /= ' ') exit
         pos = pos + 1
      end do
      do while (pos <= len(line))
         if (line(pos:pos) == ' ') exit
         pos = pos + 1
      end do
   end subroutine skip_tok

   ! tok_real(line, pos, val, ok) - read the next whitespace-separated real token
   subroutine tok_real(line, pos, val, ok)
      character(len=*), intent(in) :: line
      integer, intent(inout) :: pos
      real(8), intent(out) :: val
      logical, intent(out) :: ok
      integer :: a, b, ios
      ok = .false.
      val = 0.0d0
      a = pos
      do while (a <= len(line))
         if (line(a:a) /= ' ') exit
         a = a + 1
      end do
      if (a > len(line)) return
      b = a
      do while (b <= len(line))
         if (line(b:b) == ' ') exit
         b = b + 1
      end do
      read (line(a:b-1), *, iostat=ios) val
      ok = (ios == 0)
      pos = b
   end subroutine tok_real

   !------------------------------------------------------------------
   ! small helpers (the sibling checks' discipline)
   !------------------------------------------------------------------
   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      print '(a)', 'FAIL '//trim(msg)
      nfail = nfail + 1
   end subroutine fail

   ! write one sample file, one array element per record (materials live under build/)
   subroutine write_sample(fname, lines)
      character(len=*), intent(in) :: fname
      character(len=*), intent(in) :: lines(:)
      integer :: v, j
      open (newunit=v, file=fname, status='replace', action='write')
      do j = 1, size(lines)
         write (v, '(a)') trim(lines(j))
      end do
      close (v)
   end subroutine write_sample

   ! slurp one text file into a single-line buffer (bound-guarded; missing file
   ! leaves the buffer empty)
   subroutine slurp(fname, buf)
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
   end subroutine slurp

   ! make one sample folder (idempotent - a second run tolerates its own leftovers)
   subroutine make_dir(path)
      character(len=*), intent(in) :: path
      if (windows_host()) then
         call execute_command_line('mkdir "'//trim(path)//'" 2>nul', wait=.true.)
      else
         call execute_command_line('mkdir -p "'//trim(path)//'" 2>/dev/null', wait=.true.)
      end if
   end subroutine make_dir

   ! wipe one case folder wholesale (materials live ONLY under build/; the
   ! fresh rewrite is what makes make_cases robust against stale leftovers)
   subroutine wipe_dir(path)
      character(len=*), intent(in) :: path
      if (windows_host()) then
         call execute_command_line('rd /s /q "'//trim(path)//'" >nul 2>nul', wait=.true.)
      else
         call execute_command_line('rm -rf "'//trim(path)//'" 2>/dev/null', wait=.true.)
      end if
   end subroutine wipe_dir

   ! Windows-host probe (same runtime OS branch the sibling checks use)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

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
end program check_driver
