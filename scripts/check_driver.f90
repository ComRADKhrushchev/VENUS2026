!=====================================================================
! check_driver.f90 - criteria program for the whole driver (assembly ->
!                    evolution -> statistics; criteria defined FIRST in
!                    this header, this program executes them)
! Compile (from src/; run make FIRST - the driver exe and the generated registry
! module are build products, and make's prebuild step aggregates the dragged-in
! systems/ container into that module; -Jbuild keeps the .mod files out of the
! source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 utils/rng.f90 &
!            globals/state.f90 globals/config_roster.f90 facade/elec_ws.f90 &
!            facade/elec_facade.f90 facade/force_facade.f90 facade/sampler.f90 &
!            utils/linalg.f90 utils/geometry.f90 dynamics/propagator.f90 &
!            methods/dist_boltzmann.f90 methods/dist_ebk.f90 methods/dist_incident.f90 &
!            methods/dist_j.f90 methods/dist_rotation.f90 methods/dist_surface_oscillator.f90 &
!            methods/dist_normalmode.f90 methods/dist_barrier_excitation.f90 &
!            methods/dist_thermostat.f90 methods/dist_glo_target.f90 &
!            methods/dist_thermalize.f90 methods/elec_adiabatic.f90 &
!            methods/reg_adiabatic.f90 methods/reg_prop.f90 &
!            methods/reg_dist_boltzmann.f90 methods/reg_dist_ebk.f90 &
!            methods/reg_dist_incident.f90 methods/reg_dist_j.f90 &
!            methods/reg_dist_rotation.f90 methods/reg_dist_surface_oscillator.f90 &
!            methods/reg_dist_normalmode.f90 methods/reg_dist_barrier_excitation.f90 &
!            methods/reg_dist_thermostat.f90 methods/reg_dist_glo_target.f90 &
!            methods/reg_dist_thermalize.f90 output/recorder.f90 output/final_state.f90 &
!            ../systems/test_harmonic/harmonic_box.f90 &
!            ../systems/test_harmonic/reg_test_harmonic.f90 &
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
!       the generic hypothetical molecular system) + the dragged-in harmonic
!       verification container -> exit 0, stdout carries the assembly-summary
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
!   (c) FRAME STRUCTURE: frame headers carry traj/step/t with traj strictly
!       1..N, step strictly 1..M+1 in file order, t strictly increasing within
!       each trajectory; the column line's order is t_kin, v_pot, e_tot (the
!       recorder's selected-column legislation); frames per trajectory =
!       MAX_STEPS+1 (the safety-net-only run - this case binds no termination)
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
!       classification-bound run (E_SPLIT in the parameter file) prints a
!       non-empty outcome table 'outcome counts: <name>=<n>' (the channel is
!       what is verified - the final-state quantity computation behind fin is
!       a later landing, so the bound name reads the initialized defaults)
!   (g) SAMPLING CHAIN: (a)'s run drives samp_run through the real-input
!       DIST_SCHEME (the thermalize member on fragment 1, the rotation member
!       on fragment 2 - both on the real force chain) and the incident channel
!       for the two-fragment gas pair; observable: the FIRST frame's sec_qp
!       section carries non-zero momenta (the sampling really happened; the
!       case input sets GWRITE_LEVEL=2 to open the section - the input
!       per-key default is the minimal level 1)
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
!       container_term} - nothing else inside the integration loop
!   (h) assembly-order error loud-stop (child): a bad SYSTEM_DIR (the folder
!       does not exist) -> nonzero exit, stderr names sysdef_load - the
!       existing input/sysdef channel penetrates; the driver adds NO judgment
!       of its own
!   (i) size guard: UNREACHABLE through the driver's public face - the driver
!       derives state_create's natoms from size(roster%mass) itself, so
!       size(sta%q) = 3*size(roster%mass) holds structurally; the guard is
!       carried in the driver per the guard statute and names both sides if it
!       ever fires (the untriggerable-guard precedent: header note here, no
!       criterion body)
!   (j) seam criterion (the input-staging -> member-params assembly channel):
!       (j1) round-trip, IN-PROCESS through the same public chain in the same
!            legislated assembly order (read_input -> rng_init -> sysdef_load
!            -> roster_load -> state_create -> registry_gen_all ->
!            registry_gen_init -> container_init -> rec_select level entry):
!            the input keys T_GLO=350 / GAMMA_GLO=0.04 / W_GHOST=1.25 / DT=
!            0.05 / DIST_SCHEME=GLO_TARGET -> the member getters read the
!            exact values back (glo_target_t_glo / _gamma / _w), and
!            glo_target_dt_loaded() = DT (the dt-ordering contract); also
!            samp_n_dist() = 11 and rec_n_sel() = 3 through the same chain;
!            and container_term(sta) = .false. (the parameter file absent in
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
! Child-process scheme (the sibling checks' discipline): every loud-stop path
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
   use config,        only: reactants, observables, electronic, prop_cfg => propagator
   use config_roster, only: roster, roster_load
   use state,         only: state_t, state_create
   use rng,           only: rng_init
   use sysdef_facade, only: sysdef_load
   use registry_gen,  only: registry_gen_all, registry_gen_init
   use force_facade,  only: container_init, container_term
   use sampler,       only: samp_n_dist
   use recorder,      only: rec_select, rec_n_sel
   use final_state,   only: fin_t, container_classify
   use input,         only: read_input
   use dist_glo_target, only: glo_target_t_glo, glo_target_gamma, glo_target_w, &
                              glo_target_dt_loaded
   implicit none
   integer :: nfail = 0
   character(len=32) :: cmode
   character(len=2048) :: obuf_a = ''    ! case_a's stdout, kept for (f)
   character(len=2048) :: obuf_t = ''    ! case_t's stdout, kept for (f)

   ! ---------- child modes (self-respawn): the loud-stop paths that cannot run here ----------
   if (command_argument_count() >= 1) then
      call get_command_argument(1, cmode)
      select case (trim(cmode))
      case ('eclass')
         call child_eclass()
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
   call crit_c()      ! (c) frame structure (safety-net-only run)
   call crit_g()      ! (g) the sampling chain really ran
   call crit_d()      ! (d) MAX_STEPS safety net
   call crit_e()      ! (e) term-bound early termination (+ case_t stdout for (f))
   call crit_f()      ! (f) statistics first version (unbound + bound)
   call crit_det()    ! (det) determinism, same input twice
   call crit_e7()     ! (E7) the integration loop's direct-call set
   call crit_h()      ! (h) bad SYSTEM_DIR loud-stop
   call crit_j1()     ! (j1) seam round-trip + unbound container_term
   call crit_j2()     ! (j2) bad member value
   call crit_k()      ! (k) unbound-classification child

   if (nfail == 0) then
      print '(a)', 'check_driver: ALL PASS'
   else
      print '(a,i0,a)', 'check_driver: ', nfail, ' FAILURE(S)'
      stop 1
   end if
contains
   !------------------------------------------------------------------
   ! child_eclass() - the (k) child: the full assembly chain on a system folder
   !                  WITHOUT the container parameter file (termination and
   !                  classification unbound), then container_classify - the
   !                  unbound call must be the loud-stop authority
   subroutine child_eclass()
      type(state_t) :: stc
      type(fin_t) :: fc
      character(len=32) :: nm
      call read_input('build/drv_case_a/input_eclass.txt')
      call rng_init(i_seed(1))
      call sysdef_load(reactants%system_dir)
      call roster_load()
      call state_create(stc, size(roster%mass), max(electronic%n_surf, 1))
      call registry_gen_all()
      call registry_gen_init()     ! parameter file absent -> slots stay unbound
      call container_init()
      nm = container_classify(fc)  ! the unbound call - must loud-stop here
      print '(a)', 'child eclass: container_classify returned ', trim(nm)
   end subroutine child_eclass

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
   ! crit_c() - (c): frame header monotonicity, column order, frames =
   !             MAX_STEPS+1 per trajectory (the safety-net-only run)
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
      if (index(obuf_t, 'outcome counts: capture=2') == 0) then
         call fail('(f): the classification-bound run must print the outcome table'// &
                  ' (capture=2) - got "'//trim(obuf_t)//'"')
         return
      end if
      print '(a)', '(f): PASS - unbound run: skipped line + drift summary; bound run:'// &
           ' outcome table capture=2'
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
   ! crit_h() - (h): bad SYSTEM_DIR -> the sysdef channel's own loud-stop penetrates
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
   !              the unbound container_term posture (no parameter file here)
   subroutine crit_j1()
      type(state_t) :: st
      call read_input('build/drv_case_d/input_dpos.txt')
      call rng_init(i_seed(1))
      call sysdef_load(reactants%system_dir)
      call roster_load()
      call state_create(st, size(roster%mass), max(electronic%n_surf, 1))
      call registry_gen_all()
      call registry_gen_init()
      call container_init()
      call rec_select(observables%rec_level)
      if (len_trim(observables%obs_list) > 0) call rec_select(trim(observables%obs_list))
      if (glo_target_t_glo() /= 350.0d0 .or. glo_target_gamma() /= 0.04d0 .or. &
          glo_target_w() /= 1.25d0) then
         call fail('(j1): member getter round-trip mismatch - t_glo/gamma/w read '// &
                  d2s(glo_target_t_glo())//' '//d2s(glo_target_gamma())//' '// &
                  d2s(glo_target_w()))
      else
         print '(a)', '(j1): PASS - staged member keys round-trip: t_glo=350 gamma=0.04 w=1.25'
      end if
      if (glo_target_dt_loaded() /= 0.05d0 .or. glo_target_dt_loaded() /= prop_cfg%dt) then
         call fail('(j1): dt-ordering mismatch - glo_target_dt_loaded() = '// &
                  d2s(glo_target_dt_loaded())//' but DT staged '//d2s(prop_cfg%dt))
      else
         print '(a)', '(j1): PASS - the member init snapshotted the finalized dt (0.05)'
      end if
      if (samp_n_dist() /= 11) then
         call fail('(j1): member-table count expected 11 (10 scheme + incident), got '// &
                  i2s(samp_n_dist()))
      else if (rec_n_sel() /= 3) then
         call fail('(j1): selected column count expected 3 (the framework energy '// &
                  'columns at the default level), got '//i2s(rec_n_sel()))
      else
         print '(a)', '(j1): PASS - member table carries 11 rows; selected columns 3'
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
      character(len=80) :: xyza(4), xyzb(4), xyz1(4)
      character(len=80) :: la(14), lb(6), ld(9), ldp(9), lec(10)
      character(len=80) :: lm(9), lt(13), lpars(3)
      logical :: have_g, have_i
      character(len=64) :: srcexe, srcw
      integer :: j
      xyza = [ character(len=80) :: '2', 'generic diatomic fragment a', &
               'H 0.0 0.0 0.0', 'H 0.7 0.0 0.0' ]
      xyzb = [ character(len=80) :: '2', 'generic diatomic fragment b', &
               'H 0.0 0.0 0.0', 'H 0.8 0.0 0.0' ]
      xyz1 = [ character(len=80) :: '2', 'diatomic test fragment', &
               'H 0.0 0.0 0.0', 'H 0.7 0.0 0.0' ]
      ! (a): the end-to-end case - two-fragment gas-phase system, thermalize on
      ! fragment 1, rotation on fragment 2, the incident channel for the pair;
      ! GWRITE_LEVEL=2 opens the sec_qp flattened section (the input per-key
      ! default is the minimal level 1 - the momenta criterion reads the
      ! section); NO container parameter file -> termination and
      ! classification unbound
      la = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = THERMALIZE,ROTATION', &
             'T_EQ = 300.0', 'N_EQ = 25', 'DT_EQ = 0.05', 'J_ROT = 3', 'N_TRAJ = 3', &
             'MAX_STEPS = 40', 'ISEED = 20260915', 'GWRITE_LEVEL = 2' ]
      ! (h): the bad system_dir case (the folder does not exist)
      lb = [ character(len=80) :: 'SYSTEM_DIR = no_such_system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION' ]
      ! (j2): the bad member value (driver child, cwd-relative system dir)
      ld = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = GLO_TARGET', 'T_GLO = -5.0', &
             'GAMMA_GLO = 0.04', 'W_GHOST = 1.25' ]
      ! (j1): the positive round-trip (the check parent reads it from src/, so the
      ! system dir is build-relative)
      ldp = [ character(len=80) :: 'SYSTEM_DIR = build/drv_case_d/system', 'E_REL = 1.0', &
              'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = GLO_TARGET', &
              'T_GLO = 350.0', 'GAMMA_GLO = 0.04', 'W_GHOST = 1.25' ]
      ! (k): the eclass child input (build-relative; the child runs from src/,
      ! whose run directory carries no container parameter file; one member
      ! per fragment - a member instance carries one parameter set)
      lec = [ character(len=80) :: 'SYSTEM_DIR = build/drv_case_a/system', 'E_REL = 1.0', &
              'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.05', &
              'DIST_SCHEME = THERMALIZE,ROTATION', 'T_EQ = 300.0', 'N_EQ = 25', &
              'DT_EQ = 0.05', 'N_TRAJ = 1' ]
      ! (d): the MAX_STEPS safety-net case - one trajectory, five steps
      lm = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION', 'N_TRAJ = 1', &
             'MAX_STEPS = 5', 'ISEED = 20260915' ]
      ! (e)/(f): the term+classify bound case - a high collision energy so the
      ! amplitude climb crosses Q_TERM mid-run; MAX_STEPS=400 is the safety net
      lt = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 200.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = THERMALIZE,ROTATION', &
             'T_EQ = 300.0', 'N_EQ = 25', 'DT_EQ = 0.05', 'J_ROT = 3', 'N_TRAJ = 2', &
             'MAX_STEPS = 400', 'ISEED = 20260915' ]
      lpars = [ character(len=80) :: 'K_HARM = 0.04', 'Q_TERM = 9.5', 'E_SPLIT = 1.0' ]
      ! wipe every case folder BEFORE rebuilding it: a case folder must hold
      ! exactly what THIS run writes - a stale file from an earlier round
      ! (a renamed system sample, a relocated parameter file) would ride the
      ! system-folder scan as an extra reactant and corrupt the case
      call wipe_dir('build/drv_case_a')
      call wipe_dir('build/drv_case_b')
      call wipe_dir('build/drv_case_d')
      call wipe_dir('build/drv_case_m')
      call wipe_dir('build/drv_case_t')
      call make_dir('build/drv_case_a/system')
      call make_dir('build/drv_case_b')
      call make_dir('build/drv_case_d/system')
      call make_dir('build/drv_case_m/system')
      call make_dir('build/drv_case_t/system')
      call write_sample('build/drv_case_a/system/pair_a.xyz', xyza)
      call write_sample('build/drv_case_a/system/pair_b.xyz', xyzb)
      call write_sample('build/drv_case_d/system/h2.xyz', xyz1)
      call write_sample('build/drv_case_m/system/h2.xyz', xyz1)
      call write_sample('build/drv_case_t/system/pair_a.xyz', xyza)
      call write_sample('build/drv_case_t/system/pair_b.xyz', xyzb)
      ! the container parameter file is RUN material (cwd of the driver child),
      ! not system-definition material - the system-folder scan would reject
      ! its unknown extension
      call write_sample('build/drv_case_t/harmonic_pars.txt', lpars)
      call write_sample('build/drv_case_a/input_qct.txt', la)
      call write_sample('build/drv_case_a/input_eclass.txt', lec)
      call write_sample('build/drv_case_b/input_qct.txt', lb)
      call write_sample('build/drv_case_d/input_qct.txt', ld)
      call write_sample('build/drv_case_d/input_dpos.txt', ldp)
      call write_sample('build/drv_case_m/input_qct.txt', lm)
      call write_sample('build/drv_case_t/input_qct.txt', lt)
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
