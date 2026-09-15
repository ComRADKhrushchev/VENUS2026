!=====================================================================
! check_driver.f90 - criteria program for the driver's assembly phase (criteria
!                    defined FIRST in this header; this program executes them)
! Compile (from src/; run make FIRST - the driver exe and the generated registry
! module are build products; -Jbuild keeps the .mod files out of the source tree):
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
!            build/gfortran/registry_gen.f90 ../scripts/check_driver.f90 &
!            -o build/check_driver.exe
! Run (from src/): build\check_driver.exe
! Design choice, the ASSEMBLY-COMPLETION OBSERVABLE (stated per the brief): the
!   driver prints ONE assembly-summary line at the end of its assembly phase -
!     assembly complete: trajectories=<n> atoms=<n> fragments=<n> paradigm=<n>
!     members=<n> init=ok colgroups=<n>
!   - trajectories = control%n_traj; atoms = size(roster%mass); fragments =
!   size(roster%frag); paradigm = reactants%surface_model; members = the sampler
!   member-table row count through the NEW public getter samp_n_dist() (a
!   read-only count face; the table itself stays private); init = 'ok' - a
!   field reached only when the init aggregation returned clean (the fail-loud
!   posture makes the printed line itself the witness); colgroups = the selected
!   recorder column count through the NEW public getter rec_n_sel() (same
!   minimal-count-face shape). The criteria grep this line. The two count
!   getters are part of this observable design, not general face growth.
! Criteria:
!   (a) assembly chain through (the REAL driver exe, child process): a real
!       input file + a real system folder under build/ (build/drv_case_a: one
!       diatomic .xyz = the single-fragment gas-phase system, DIST_SCHEME=
!       ROTATION, N_TRAJ=3) -> exit 0, no loud-stop, stdout carries the
!       assembly-summary line with the expected field values (trajectories=3
!       atoms=2 fragments=1 paradigm=0 members=11 init=ok colgroups=3; members
!       = the 10 scheme members + the incident row; colgroups = 3 because the
!       recorder's lazy framework preregistration puts t_kin/v_pot/e_tot into
!       every preset level and the default rec_level=2 selects them - no
!       external columns exist tree-wide, so 3 is the true reading)
!   (b) assembly-order error loud-stop (child, the real driver exe): a bad
!       SYSTEM_DIR (the folder does not exist) -> nonzero exit, stderr names
!       sysdef_load - the existing input/sysdef channel penetrates; the driver
!       adds NO new judgment of its own
!   (c) size guard: UNREACHABLE through the driver's public face - the driver
!       derives state_create's natoms from size(roster%mass) itself, so
!       size(sta%q) = 3*size(roster%mass) holds structurally; the guard is
!       carried in the driver per the guard statute and names both sides if it
!       ever fires (the untriggerable-guard precedent: header note here, no
!       criterion body)
!   (d) seam criterion (the input-staging -> member-params assembly channel):
!       (d1) round-trip, IN-PROCESS through the same public chain in the same
!            legislated assembly order (read_input -> rng_init -> sysdef_load
!            -> roster_load -> state_create -> registry_gen_all ->
!            registry_gen_init -> container_init -> rec_select level entry):
!            the input keys T_GLO=350 / GAMMA_GLO=0.04 / W_GHOST=1.25 / DT=
!            0.05 / DIST_SCHEME=GLO_TARGET -> the member getters read the
!            exact values back (glo_target_t_glo / _gamma / _w), and
!            glo_target_dt_loaded() = DT (the dt-ordering contract: dt is
!            finalized before the member init that snapshots it); also
!            samp_n_dist() = 11 and rec_n_sel() = 3 through the same chain
!       (d2) bad value (child, the real driver exe): T_GLO = -5 -> nonzero
!            exit, stderr names glo_target_init and the field t_glo
!   (e) unbound-container posture (the idle chain - no systems/ container
!       dragged in anywhere):
!       (e1) covered by (a): the same child run binds no PES, the assembly
!            completes, exit 0 (assembly neither binds nor crashes)
!       (e2) force_eval first-call authority (this check's OWN child mode
!            'eforce', self-respawn): the full unbound assembly chain, then
!            force_eval -> nonzero exit, stderr names force_eval
!       (e3) container_term unbound = .false., unit-level through the public
!            face: after the in-process (d1) chain, container_term(sta)
!            returns .false. without stopping
! Child-process scheme (the sibling checks' discipline): every loud-stop path
!   runs in a child - driver-exe children run with cwd = their case directory
!   (the driver reads ./input_qct.txt relative to its cwd; 'cd <case> && exe'
!   carries both shell families); the build names the driver venus.e (no .exe
!   extension - Windows cmd cannot execute that name), so the setup copies the
!   build product to build/drv_venus.exe once before spawning. The check's own
!   child ('eforce') runs from src/ and reads its input by build-relative path.
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
   use force_facade,  only: container_init, force_eval, container_term
   use sampler,       only: samp_n_dist
   use recorder,      only: rec_select, rec_n_sel
   use input,         only: read_input
   use dist_glo_target, only: glo_target_t_glo, glo_target_gamma, glo_target_w, &
                              glo_target_dt_loaded
   implicit none
   integer :: nfail = 0
   character(len=32) :: cmode

   ! ---------- child modes (self-respawn): the loud-stop paths that cannot run here ----------
   if (command_argument_count() >= 1) then
      call get_command_argument(1, cmode)
      select case (trim(cmode))
      case ('eforce')
         call child_eforce()
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
   call crit_a()      ! (a) + (e1)
   call crit_b()      ! (b)
   call crit_d1()     ! (d1) + (e3)
   call crit_d2()     ! (d2)
   call crit_e2()     ! (e2)

   if (nfail == 0) then
      print '(a)', 'check_driver: ALL PASS'
   else
      print '(a,i0,a)', 'check_driver: ', nfail, ' FAILURE(S)'
      stop 1
   end if
contains
   !------------------------------------------------------------------
   ! child_eforce() - the (e2) child: the full unbound assembly chain (the legislated
   !                  order), then force_eval - the FIRST force call must be
   !                  the loud-stop authority (the class-A slot was never bound)
   subroutine child_eforce()
      type(state_t) :: stc
      real(8) :: e
      call read_input('build/drv_case_e/input_qct.txt')
      call rng_init(i_seed(1))
      call sysdef_load(reactants%system_dir)
      call roster_load()
      call state_create(stc, size(roster%mass), max(electronic%n_surf, 1))
      call registry_gen_all()
      call registry_gen_init()
      call container_init()      ! unbound class-A slot: the idle chain proceeds
      call force_eval(stc, e)    ! the first force call - must loud-stop here
      print '(a)', 'child eforce: force_eval returned (it must loud-stop)'
   end subroutine child_eforce

   !------------------------------------------------------------------
   ! crit_a() - (a) + (e1): the real driver exe through the good case
   subroutine crit_a()
      character(len=2048) :: obuf, ebuf
      character(len=256) :: want
      integer :: stat
      call run_driver('drv_case_a', obuf, ebuf, stat)
      want = 'assembly complete: trajectories=3 atoms=2 fragments=1 paradigm=0'// &
             ' members=11 init=ok colgroups=3'
      if (stat /= 0) then
         call fail('(a): the driver child exited nonzero - stderr: '//trim(ebuf))
      else if (index(obuf, trim(want)) == 0) then
         call fail('(a): stdout missing the summary line "'//trim(want)//'" - got "'// &
                  trim(obuf)//'"')
      else
         print '(a)', '(a)+(e1): PASS - the unbound assembly ran through the real driver;'// &
              ' summary line verified (no container bound anywhere)'
      end if
   end subroutine crit_a

   !------------------------------------------------------------------
   ! crit_b() - (b): bad SYSTEM_DIR -> the sysdef channel's own loud-stop penetrates
   subroutine crit_b()
      character(len=2048) :: obuf, ebuf
      integer :: stat
      call run_driver('drv_case_b', obuf, ebuf, stat)
      if (stat == 0) then
         call fail('(b): expected a loud-stop, the driver child exited 0')
      else if (index(ebuf, 'sysdef_load') == 0) then
         call fail('(b): stderr does not name sysdef_load - got "'//trim(ebuf)//'"')
      else
         print '(a)', '(b): PASS ('//trim(ebuf)//')'
      end if
   end subroutine crit_b

   !------------------------------------------------------------------
   ! crit_d1() - (d1) + (e3): the in-process seam round-trip (the legislated order) and
   !                          the unbound container_term posture
   subroutine crit_d1()
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
         call fail('(d1): member getter round-trip mismatch - t_glo/gamma/w read '// &
                  d2s(glo_target_t_glo())//' '//d2s(glo_target_gamma())//' '// &
                  d2s(glo_target_w()))
      else
         print '(a)', '(d1): PASS - staged member keys round-trip: t_glo=350 gamma=0.04 w=1.25'
      end if
      if (glo_target_dt_loaded() /= 0.05d0 .or. glo_target_dt_loaded() /= prop_cfg%dt) then
         call fail('(d1): dt-ordering mismatch - glo_target_dt_loaded() = '// &
                  d2s(glo_target_dt_loaded())//' but DT staged '//d2s(prop_cfg%dt))
      else
         print '(a)', '(d1): PASS - the member init snapshotted the finalized dt (0.05)'
      end if
      if (samp_n_dist() /= 11) then
         call fail('(d1): member-table count expected 11 (10 scheme + incident), got '// &
                  i2s(samp_n_dist()))
      else if (rec_n_sel() /= 3) then
         call fail('(d1): selected column count expected 3 (the framework energy '// &
                  'columns at the default level), got '//i2s(rec_n_sel()))
      else
         print '(a)', '(d1): PASS - member table carries 11 rows; selected columns 3'
      end if
      if (container_term(st)) then
         call fail('(e3): unbound container_term must return .false.')
      else
         print '(a)', '(e3): PASS - container_term on the unbound slot returns .false.'
      end if
   end subroutine crit_d1

   !------------------------------------------------------------------
   ! crit_d2() - (d2): T_GLO = -5 through the real driver exe
   subroutine crit_d2()
      character(len=2048) :: obuf, ebuf
      integer :: stat
      call run_driver('drv_case_d', obuf, ebuf, stat)
      if (stat == 0) then
         call fail('(d2): expected a loud-stop, the driver child exited 0')
      else if (index(ebuf, 'glo_target_init') == 0) then
         call fail('(d2): stderr does not name glo_target_init - got "'//trim(ebuf)//'"')
      else if (index(ebuf, 't_glo') == 0) then
         call fail('(d2): stderr missing the field name t_glo - got "'//trim(ebuf)//'"')
      else
         print '(a)', '(d2): PASS ('//trim(ebuf)//')'
      end if
   end subroutine crit_d2

   !------------------------------------------------------------------
   ! crit_e2() - (e2): this check's own child runs the unbound chain then force_eval
   subroutine crit_e2()
      character(len=512) :: cmd
      character(len=2048) :: ebuf
      character(len=1) :: sep
      integer :: stat
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_driver.exe eforce 2>build'//sep//'drv_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp('build/drv_err.txt', ebuf)
      if (stat == 0) then
         call fail('(e2): expected a loud-stop, the eforce child exited 0')
      else if (index(ebuf, 'force_eval') == 0) then
         call fail('(e2): stderr does not name force_eval - got "'//trim(ebuf)//'"')
      else
         print '(a)', '(e2): PASS ('//trim(ebuf)//')'
      end if
   end subroutine crit_e2

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
   !                  execute the build's venus.e name directly)
   subroutine make_cases()
      character(len=80) :: xyz(4), la(7), lb(6), ld(9), ldp(9), le(7)
      xyz = [ character(len=80) :: '2', 'diatomic test fragment', &
              'H 0.0 0.0 0.0', 'H 0.7 0.0 0.0' ]
      ! (a): the good case - single-fragment gas-phase system, rotation member
      la = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION', 'N_TRAJ = 3' ]
      ! (b): the bad system_dir case (the folder does not exist)
      lb = [ character(len=80) :: 'SYSTEM_DIR = no_such_system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION' ]
      ! (d2): the bad member value (driver child, cwd-relative system dir)
      ld = [ character(len=80) :: 'SYSTEM_DIR = system', 'E_REL = 1.0', 'B_MAX = 2.0', &
             'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = GLO_TARGET', 'T_GLO = -5.0', &
             'GAMMA_GLO = 0.04', 'W_GHOST = 1.25' ]
      ! (d1): the positive round-trip (the check parent reads it from src/, so the
      ! system dir is build-relative)
      ldp = [ character(len=80) :: 'SYSTEM_DIR = build/drv_case_d/system', 'E_REL = 1.0', &
              'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = GLO_TARGET', &
              'T_GLO = 350.0', 'GAMMA_GLO = 0.04', 'W_GHOST = 1.25' ]
      ! (e2): the eforce child input (reuses case a's system folder from src/)
      le = [ character(len=80) :: 'SYSTEM_DIR = build/drv_case_a/system', 'E_REL = 1.0', &
             'B_MAX = 2.0', 'R_SEP = 8.0', 'DT = 0.05', 'DIST_SCHEME = ROTATION', &
             'N_TRAJ = 1' ]
      call make_dir('build/drv_case_a/system')
      call make_dir('build/drv_case_b')
      call make_dir('build/drv_case_d/system')
      call make_dir('build/drv_case_e')
      call write_sample('build/drv_case_a/system/h2.xyz', xyz)
      call write_sample('build/drv_case_d/system/h2.xyz', xyz)
      call write_sample('build/drv_case_a/input_qct.txt', la)
      call write_sample('build/drv_case_b/input_qct.txt', lb)
      call write_sample('build/drv_case_d/input_qct.txt', ld)
      call write_sample('build/drv_case_d/input_dpos.txt', ldp)
      call write_sample('build/drv_case_e/input_qct.txt', le)
      if (windows_host()) then
         call execute_command_line('copy /y "build\gfortran\venus.e" "build\drv_venus.exe"'// &
                                   ' >nul 2>nul', wait=.true.)
      else
         call execute_command_line('cp -f build/gfortran/venus.e build/drv_venus.exe', &
                                   wait=.true.)
      end if
   end subroutine make_cases

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
