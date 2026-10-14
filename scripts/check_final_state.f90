!=====================================================================
! check_final_state.f90 - criteria program for the final-state quantity
!                         computation fin_eval (fragment internal energies /
!                         angular momenta, generalized COM relative quantities,
!                         scattering-angle family, the COM-frame archive) and
!                         its classification channel (criteria defined FIRST
!                         in this header; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files
! out of the source tree; the well constant is buffered check-side through the
! potential module's own setter - a NON-default value, so the buffered-potential
! anchors genuinely ride the bound container):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 utils/rng.f90 &
!            globals/state.f90 globals/config_atoms.f90 &
!            interface/elec_interface.f90 methods/elec_adiabatic.f90 &
!            utils/linalg.f90 interface/force_interface.f90 &
!            utils/geometry.f90 output/final_state.f90 &
!            ../systems/reference/reference_pot.f90 &
!            ../systems/reference/reference_box.f90 &
!            ../scripts/check_final_state.f90 -o build/check_final_state.exe
! Run (from src/): build\check_final_state.exe
!
! Constructed material (criterion (a)'s state, the anchor derivations):
!   masses [1,2 | 3,4]; fragment A = atoms (1,2) at (1,0,0)/(-0.5,0,0)
!   -> COM (0,0,0); fragment B = atoms (3,4) at (0,1,7)/(0,-0.75,7)
!   -> COM (0,0,7). COM velocities v_cm_A = (0.1,-0.05,0.05),
!   v_cm_B = (-0.2,0.15,0.02); rigid rotations omega_a = 0.4 z-hat,
!   omega_b = 0.25 x-hat (momenta p_i = m_i*(v_cm + omega x r_i)):
!     J_A = (0,0,0.6) = I_A*omega with I_A = 1*1^2+2*0.5^2 = 1.5;
!     E_rot_A = internal kinetic = 0.5*1.5*0.4^2 = 0.12 (internal units) - the
!       rotational energy (inertia-tensor path) and the internal kinetic
!       (COM-frame momentum path) are TWO computational paths to the same
!       number for a pure rigid rotation: a cross-identity anchor;
!     J_B = (1.3125,0,0) = I_B*omega with I_B = 3*1^2+4*0.75^2 = 5.25;
!     E_rot_B = internal kinetic = 0.5*5.25*0.25^2 = 0.1640625;
!     relative motion (gas mu = 3*7/10 = 2.1): q_r = (0,0,-7),
!     v_r = (0.3,-0.2,0.03), radial speed -0.03 -> e_rel = 2.1*9e-4/2,
!     e_rel_sq = 2.1*0.1309/2, oam = 2.1*(-1.4,-2.1,0),
!     b_fin = |q x v|/|v| = sqrt(6.37)/sqrt(0.1309) (mu-invariant);
!     e_cm = |(0.3,-0.15,0.15)+(-1.4,1.05,0.14)|^2/(2*10) = 2.1041/20;
!     ang(1) = angle(oam, j_a) = 90 deg, ang(2) = angle(oam, j_b) with
!       cos = -1.4/sqrt(6.37), ang(3) = angle(j_a, j_b) = 90 deg;
!     the equilibrium geometries qz_A = (0.9,0,-0.45,0 pattern),
!     qz_B = the current B pattern: the buffered internal potential of A keeps
!       A current and places B at qz_B + COM_B (= the current B sites:
!       sum |q|^2 = 100.8125), the buffered potential of B places A at
!       qz_A + COM_A (sum |q|^2 = 0.81+0.2025+99.5625 = 100.575) - the two
!       buffered anchors DIFFER, so the buffer really reads the fragment table;
!     the archive: q_fin = q - q_cm per fragment, p_fin = p - m*v_cm with the
!       rigid rotation subtracted first (the linear-fragment omega channel
!       carries 0 - the two-atom branch computes no omega - so p_fin RETAINS
!       the rotation: (0,0.4,0)/(0,-0.4,0) for A, (0,0,0.75)/(0,0,-0.75) for B);
!     t_life = t*10 = 37 at t = 3.7.
!   Surface leg (same state, surface_model = 1): mu = m_a = 3.0 ->
!   e_rel = 3*9e-4/2, e_rel_sq = 3*0.1309/2, oam = 3*(-1.4,-2.1,0);
!   the scattering angles from v_cm_A = (0.1,-0.05,0.05):
!   ang(19) = acos(0.05/sqrt(0.015)) [deg], ang(20) = 2*pi -
!   acos(0.1/sqrt(0.0125)) (the negative-y azimuth branch) [deg].
!   Criterion (b)'s states: two one-atom fragments (m = 2 at (0,0,5),
!   m = 3 at (0,0,-5), mu = 1.2): PARALLEL v_cm (0,0,0.2)/(0,0,-0.1) - the
!   whole relative motion is radial (e_rel = e_rel_sq = 1.2*0.09/2, oam = 0,
!   b_fin = 0, all internal fields skipped-zero for one-atom fragments);
!   PERPENDICULAR v_cm (0.2,0,0)/(0,0,0) - purely tangential (e_rel = 0
!   EXACTLY, e_rel_sq = 1.2*0.04/2, oam = (0,2.4,0), b_fin = |q x v|/|v| = 10
!   = the full separation: the two energy-splitting limits the classification
!   threshold consumes).
!   Criterion (c)'s state: ONE diatomic fragment (m = [1,2], the A geometry,
!   v_cm = (0.05,0,0), omega = 0.4 z-hat): j_a = (0,0,0.6), e_rot_a =
!   internal kinetic = 0.12, buffered potential = K/2*1.25, e_cm = 0.15^2/6,
!   and EVERY fragment-B / relative / angle field skipped-zero.
!
! Tolerances (derivations):
!   - scalars: rel 1e-10 (floor 1.0). Every field is a fixed short chain of
!     double-precision multiplies/adds over 2-4 atoms; per-operation roundoff
!     ~1e-15 rel, at most a few dozen operations per value - 1e-10 leaves
!     five decades of margin while a wrong formula (wrong mass factor, wrong
!     projection, wrong fold) deviates at the 1e-2..1e0 level.
!   - angles: abs 1e-9 deg (~1.7e-11 rad). The anchors re-derive the cosine
!     through a different arithmetic order than the implementation; the
!     constructed cosines sit away from the acos poles (|cos| <= 0.71), where
!     the condition number is O(1).
!   - bitwise zeros: quantities skipped by construction (one-atom fragments,
!     the not-yet-wired initial-value angle channels, the degenerate
!     relative-motion fields) must be EXACTLY 0.0 - the receiver's fields are
!     default-initialized and the skip must leave them untouched; any nonzero
!     bit means a wrong branch executed.
!   - the eV->kcal/mol conversion is CALIBRATED check-side (one force_eval at
!     a known single-offset geometry, p = 0): the conversion literal never
!     lands in this file (the unit-contract grep discipline).
!
! Criteria:
!   (a) constructed-state per-field anchors: the two-fragment state above ->
!       every fin field checked against its hand-derived anchor (internal
!       energies incl. the buffered potential via the measured conversion,
!       angular momenta with the hbar fold, relative quantities, e_cm, the
!       final-only scattering angles, the archive, t_life); then the SURFACE
!       leg re-runs the same state with surface_model = 1 -> the mu = m_a
!       relative quantities and the ang(19)/ang(20) scattering pair anchored
!       (gas leaves them exactly 0). The not-yet-wired channels (the
!       initial-value angles 4..18, the surface-normal angles 17..18) are
!       asserted bitwise-zero as their deferred-channel contract
!   (b) energy-splitting limits: the parallel one-atom-pair state (all
!       relative energy in the radial/translational channel: e_rel = e_rel_sq,
!       oam = 0 bitwise) vs the perpendicular state (e_rel = 0 bitwise, the
!       energy in the tangential channel, oam/b_fin maximal: b_fin = the full
!       separation) - the two limits the classification threshold consumes;
!       one-atom fragments leave every internal/angular field bitwise-zero
!       and the archive rows exactly zero
!   (c) the single-fragment degenerate path: the diatomic state above ->
!       fragment-A fields anchored, every fragment-B / relative-motion /
!       angle field bitwise-zero, e_cm and t_life still computed (they are
!       whole-state quantities, not pair quantities)
!   (d) the classification channel end to end (child 'cls', cwd =
!       build/fs_cls/ with the container parameter file buffered E_SPLIT = 0.5
!       - between the two limit values of criterion (b)'s e_rel): the child
!       binds the reference container's judgment through the public slot,
!       runs fin_eval on BOTH limit states and classifies through the REAL
!       bound entry -> 'scattering' for the parallel state (e_rel ~ 1.29
!       kcal/mol > 0.5) and 'capture' for the perpendicular state (e_rel = 0):
!       the dummy judgment consumes REAL fin values
!   (e) named abort faces (children): 'eun' - fin_eval before any list_atoms
!       assembly; 'esz' - a state whose size disagrees with the atom list atom
!       count; 'e3f' - a three-fragment list_atoms (the receiver carries two
!       fragment channels). Each child drives the real guard path; stderr
!       must name fin_eval and the exit must be nonzero
! Child-process scheme (the sibling checks' discipline): named abort paths and
!   the single-container-per-process classification case run as self-respawn
!   children (execute_command_line, output redirected under build/); a child
!   that RETURNS from a path expected to named abort prints the bug witness and
!   exits 0 - the parent's exitstat assertion fails.
! Registered deviations from the archived source semantics (documented, part
!   of this check's contract): the angular-magnitude skip threshold 1e-5
!   [hbar] and the cosine clamp [-1,1] are applied to ALL final-value angle
!   terms (the archived form clamped only its first term and carried no
!   final-magnitude threshold - the initial-value channels that governed
!   thresholds are not wired yet); the surface scattering pair rides the
!   fragment-A COM velocity uniformly (the archived surface branch
!   reconstructed a crossing velocity - that branch is excised); the linear
!   (two-atom) fragment keeps its rotation in the archive (the archived
!   rotation routine computed no omega for linear fragments - carried over
!   faithfully); the internal potential carries no reference-energy
!   subtraction yet (the reference lands with the assembly-phase channel).
!=====================================================================
program check_final_state
   use consts,        only: e_conv, hbar_code, dtor, two_pi
   use state,         only: state_t, state_create, state_destroy
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use elec_interface,   only: method_bind
   use force_interface,  only: force_eval, container_bind_pes
   use elec_adiabatic, only: adiabatic_register => member_register
   use reference_pot, only: harm_set_constant
   use reference_box, only: reference_params_t, reference_pes, reference_load, reference_classify
   use final_state,   only: fin_t, fin_eval, container_bind_classify, &
                            container_classify
   implicit none
   integer :: nfail = 0
   integer :: nfail0                     ! per-criterion failure watermark
   character(len=32) :: mode             ! the child spawn's mode argument
   ! the buffered well constant [eV/A^2] - deliberately NOT the container's
   ! default, so the buffered-potential anchors ride the check-buffered value
   real(8), parameter :: kk = 0.123d0

   ! ---------- child dispatch (one guard/classification case per spawn) ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('cls')
         call child_cls()
      case ('eun')
         call child_eun()
      case ('esz')
         call child_esz()
      case ('e3f')
         call child_e3f()
      case default
         print '(a)', 'check_final_state: unknown child mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the child chain returned - bug)'
      stop 0
   end if

   ! ---------- parent assembly: the real degenerate method row + the reference
   ! container bound through the public face + the check-buffered well constant ----------
   call adiabatic_register()
   call method_bind('adiabatic')
   call container_bind_pes(reference_pes)
   call harm_set_constant(kk)

   ! ---------- criteria ----------
   call crit_a()      ! (a) constructed-state anchors, gas + surface legs
   call crit_b()      ! (b) the energy-splitting limits (one-atom pair)
   call crit_c()      ! (c) the single-fragment degenerate path
   call crit_d()      ! (d) the classification channel end to end (child)
   call crit_e()      ! (e) the named abort faces (children)

   if (nfail == 0) then
      print '(a)', 'check_final_state: ALL PASS'
   else
      print '(a,i0,a)', 'check_final_state: ', nfail, ' FAILURE(S)'
      stop 1
   end if
contains
   !------------------------------------------------------------------
   ! crit_a() - (a): the two-fragment constructed state, gas leg then the
   !             surface leg (same state, surface_model = 1)
   subroutine crit_a()
      type(state_t) :: st
      type(fin_t) :: fin
      real(8) :: cal, e_tot
      real(8), parameter :: s637 = 6.37d0, s1309 = 0.1309d0
      nfail0 = nfail
      ! the atom list: 4 atoms, 2 fragments, NON-current equilibrium geometry
      ! for A (the buffer differentiator), current-equivalent for B
      call wipe_list_atoms()
      allocate (list_atoms%mass(4))
      list_atoms%mass = [ 1.0d0, 2.0d0, 3.0d0, 4.0d0 ]
      allocate (list_atoms%frag(2))
      list_atoms%frag(1)%nat = 2
      allocate (list_atoms%frag(1)%list(2))
      list_atoms%frag(1)%list = [ 1, 2 ]
      allocate (list_atoms%frag(1)%qz(6))
      list_atoms%frag(1)%qz = [ 0.9d0, 0.0d0, 0.0d0, -0.45d0, 0.0d0, 0.0d0 ]
      list_atoms%frag(1)%mass = 3.0d0
      list_atoms%frag(2)%nat = 2
      allocate (list_atoms%frag(2)%list(2))
      list_atoms%frag(2)%list = [ 3, 4 ]
      allocate (list_atoms%frag(2)%qz(6))
      list_atoms%frag(2)%qz = [ 0.0d0, 1.0d0, 0.0d0, 0.0d0, -0.75d0, 0.0d0 ]
      list_atoms%frag(2)%mass = 7.0d0
      list_atoms%i_proj = 1                   ! role rows: the surface quantities follow i_proj
      list_atoms%i_targ = 2
      ! the conversion calibration: one force_eval at a single-offset geometry
      ! with p = 0 (V = K/2*|q|^2 eV) measures the composed kcal/mol-per-eV
      ! factor through the bound container - the literal never lands here
      call state_create(st, 4)
      st%q = [ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
               0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 ]
      st%p = 0.0d0
      call force_eval(st, e_tot)
      cal = e_tot*2.0d0/kk
      call state_destroy(st)
      if (cal <= 0.0d0) then
         call fail('(a) measured conversion constant non-positive: '//d2s(cal))
         return
      end if
      ! rebuild the constructed state (the calibration consumed it)
      call state_create(st, 4)
      st%q = [ 1.0d0, 0.0d0, 0.0d0, -0.5d0, 0.0d0, 0.0d0, &
               0.0d0, 1.0d0, 7.0d0, 0.0d0, -0.75d0, 7.0d0 ]
      st%p = [ 0.1d0, 0.35d0, 0.05d0, 0.2d0, -0.5d0, 0.1d0, &
               -0.6d0, 0.45d0, 0.81d0, -0.8d0, 0.6d0, -0.67d0 ]
      st%t = 3.7d0
      reactants%surface_model = 0
      call fin_eval(st, fin)
      ! ---- gas-leg anchors ----
      call chk_r('(a) e_int_a(1) [internal kinetic, COM-frame path]', &
                 fin%e_int_a(1), 0.12d0/e_conv)
      call chk_r('(a) e_rot_a [inertia-tensor path - the cross-identity]', &
                 fin%e_rot_a, 0.12d0/e_conv)
      call chk_r('(a) j_a(1)', fin%j_a(1), 0.0d0/hbar_code)
      call chk_r('(a) j_a(2)', fin%j_a(2), 0.0d0/hbar_code)
      call chk_r('(a) j_a(3)', fin%j_a(3), 0.6d0/hbar_code)
      call chk_r('(a) j_mag_a', fin%j_mag_a, 0.6d0/hbar_code)
      call chk_r('(a) e_int_b(1)', fin%e_int_b(1), 0.1640625d0/e_conv)
      call chk_r('(a) e_rot_b', fin%e_rot_b, 0.1640625d0/e_conv)
      call chk_r('(a) j_b(1)', fin%j_b(1), 1.3125d0/hbar_code)
      call chk_r('(a) j_b(2)', fin%j_b(2), 0.0d0/hbar_code)
      call chk_r('(a) j_b(3)', fin%j_b(3), 0.0d0/hbar_code)
      call chk_r('(a) j_mag_b', fin%j_mag_b, 1.3125d0/hbar_code)
      call chk_r('(a) e_int_a(2) [buffered: B at its equilibrium sites]', &
                 fin%e_int_a(2), cal*(kk/2.0d0)*100.8125d0)
      call chk_r('(a) e_int_b(2) [buffered: A at its equilibrium sites]', &
                 fin%e_int_b(2), cal*(kk/2.0d0)*100.575d0)
      call chk_r('(a) e_int_a(3)', fin%e_int_a(3), 0.12d0/e_conv + cal*(kk/2.0d0)*100.8125d0)
      call chk_r('(a) e_int_b(3)', fin%e_int_b(3), 0.1640625d0/e_conv + cal*(kk/2.0d0)*100.575d0)
      call chk_r('(a) e_rel [radial projection]', fin%e_rel, 0.000945d0/e_conv)
      call chk_r('(a) e_rel_sq [full velocity]', fin%e_rel_sq, 0.137445d0/e_conv)
      call chk_r('(a) oam(1)', fin%oam(1), -2.94d0/hbar_code)
      call chk_r('(a) oam(2)', fin%oam(2), -4.41d0/hbar_code)
      call chk_z('(a) oam(3)', fin%oam(3))
      call chk_r('(a) oam_mag', fin%oam_mag, 2.1d0*sqrt(s637)/hbar_code)
      call chk_r('(a) b_fin', fin%b_fin, sqrt(s637)/sqrt(s1309))
      call chk_r('(a) e_cm', fin%e_cm, 0.105205d0/e_conv)
      call chk_deg('(a) ang(1) [oam . j_a]', fin%ang(1), 90.0d0)
      call chk_deg('(a) ang(2) [oam . j_b]', fin%ang(2), acos(-1.4d0/sqrt(s637))/dtor)
      call chk_deg('(a) ang(3) [j_a . j_b]', fin%ang(3), 90.0d0)
      call chk_arr_z('(a) ang(4..18) [initial-value channels, deferred]', fin%ang(4:18))
      call chk_z('(a) ang(19) [gas: no surface scattering pair]', fin%ang(19))
      call chk_z('(a) ang(20) [gas: no surface scattering pair]', fin%ang(20))
      call chk_r('(a) t_life', fin%t_life, 37.0d0)
      call chk_archive('(a)', fin, &
                       [ 1.0d0, 0.0d0, 0.0d0, -0.5d0, 0.0d0, 0.0d0, &
                         0.0d0, 1.0d0, 0.0d0, 0.0d0, -0.75d0, 0.0d0 ], &
                       [ 0.0d0, 0.4d0, 0.0d0, 0.0d0, -0.4d0, 0.0d0, &
                         0.0d0, 0.0d0, 0.75d0, 0.0d0, 0.0d0, -0.75d0 ])
      if (nfail == nfail0) then
         print '(a)', '(a) gas leg: PASS (every fin field on its hand-derived anchor;'// &
            ' buffered potentials differ per fragment table row; archive retains the'// &
            ' linear-fragment rotation; deferred channels bitwise-zero)'
      end if
      ! ---- surface leg: same state, mu = m_a, the scattering pair ----
      nfail0 = nfail
      reactants%surface_model = 1
      call fin_eval(st, fin)
      reactants%surface_model = 0
      call chk_r('(a/s) e_rel [mu = m_a]', fin%e_rel, 0.00135d0/e_conv)
      call chk_r('(a/s) e_rel_sq [mu = m_a]', fin%e_rel_sq, 0.19635d0/e_conv)
      call chk_r('(a/s) oam(1)', fin%oam(1), -4.2d0/hbar_code)
      call chk_r('(a/s) oam(2)', fin%oam(2), -6.3d0/hbar_code)
      call chk_r('(a/s) oam_mag', fin%oam_mag, 3.0d0*sqrt(s637)/hbar_code)
      call chk_r('(a/s) b_fin [mu-invariant]', fin%b_fin, sqrt(s637)/sqrt(s1309))
      call chk_deg('(a/s) ang(19) [scattering polar]', fin%ang(19), &
                   acos(0.05d0/sqrt(0.015d0))/dtor)
      call chk_deg('(a/s) ang(20) [scattering azimuth, negative-y branch]', fin%ang(20), &
                   (two_pi - acos(0.1d0/sqrt(0.0125d0)))/dtor)
      call chk_r('(a/s) e_int_a(2) [buffer untouched by the paradigm]', &
                 fin%e_int_a(2), cal*(kk/2.0d0)*100.8125d0)
      if (nfail == nfail0) then
         print '(a)', '(a) surface leg: PASS (mu = m_a relative quantities, the'// &
            ' scattering polar/azimuth pair from the fragment-A COM velocity;'// &
            ' the negative-y azimuth branch exercised)'
      end if
      ! ---- surface role-mirror leg: the projectile role on fragment B rescales
      ! the surface relative quantities by m_b/m_a = 7/3 (b_fin mu-invariant) ----
      nfail0 = nfail
      list_atoms%i_proj = 2
      reactants%surface_model = 1
      call fin_eval(st, fin)
      reactants%surface_model = 0
      list_atoms%i_proj = 1
      call chk_r('(a/r) e_rel [mu = m_b: the role row drives]', fin%e_rel, &
                 (7.0d0/3.0d0)*0.00135d0/e_conv)
      call chk_r('(a/r) e_rel_sq [mu = m_b]', fin%e_rel_sq, (7.0d0/3.0d0)*0.19635d0/e_conv)
      call chk_r('(a/r) oam_mag [mu = m_b]', fin%oam_mag, 7.0d0*sqrt(s637)/hbar_code)
      call chk_r('(a/r) b_fin [mu-invariant]', fin%b_fin, sqrt(s637)/sqrt(s1309))
      if (nfail == nfail0) then
         print '(a)', '(a) surface role-mirror leg: PASS (i_proj = 2: the surface'// &
            ' mu follows the resolved role row, b_fin unchanged)'
      end if
      call state_destroy(st)
   end subroutine crit_a

   !------------------------------------------------------------------
   ! crit_b() - (b): the energy-splitting limits on a one-atom pair
   subroutine crit_b()
      type(state_t) :: st
      type(fin_t) :: fp, fq
      nfail0 = nfail
      call wipe_list_atoms()
      allocate (list_atoms%mass(2))
      list_atoms%mass = [ 2.0d0, 3.0d0 ]
      allocate (list_atoms%frag(2))
      list_atoms%frag(1)%nat = 1
      allocate (list_atoms%frag(1)%list(1))
      list_atoms%frag(1)%list = [ 1 ]
      allocate (list_atoms%frag(1)%qz(3))
      list_atoms%frag(1)%qz = 0.0d0
      list_atoms%frag(1)%mass = 2.0d0
      list_atoms%frag(2)%nat = 1
      allocate (list_atoms%frag(2)%list(1))
      list_atoms%frag(2)%list = [ 2 ]
      allocate (list_atoms%frag(2)%qz(3))
      list_atoms%frag(2)%qz = 0.0d0
      list_atoms%frag(2)%mass = 3.0d0
      call state_create(st, 2)
      st%q = [ 0.0d0, 0.0d0, 5.0d0, 0.0d0, 0.0d0, -5.0d0 ]
      st%t = 1.1d0
      ! PARALLEL: the whole relative motion radial
      st%p = [ 0.0d0, 0.0d0, 0.4d0, 0.0d0, 0.0d0, -0.3d0 ]
      call fin_eval(st, fp)
      ! PERPENDICULAR: purely tangential
      st%p = [ 0.4d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 ]
      call fin_eval(st, fq)
      ! --- parallel anchors ---
      call chk_r('(b) par e_rel = e_rel_sq [all-radial limit]', &
                 fp%e_rel, 0.054d0/e_conv)
      call chk_r('(b) par e_rel_sq', fp%e_rel_sq, 0.054d0/e_conv)
      call chk_arr_z('(b) par e_int_a [one-atom skip]', fp%e_int_a)
      call chk_arr_z('(b) par e_int_b [one-atom skip]', fp%e_int_b)
      call chk_arr_z('(b) par j_a', fp%j_a)
      call chk_z('(b) par j_mag_a', fp%j_mag_a)
      call chk_z('(b) par e_rot_a', fp%e_rot_a)
      call chk_arr_z('(b) par j_b', fp%j_b)
      call chk_z('(b) par j_mag_b', fp%j_mag_b)
      call chk_z('(b) par e_rot_b', fp%e_rot_b)
      call chk_arr_z('(b) par oam', fp%oam)
      call chk_z('(b) par oam_mag', fp%oam_mag)
      call chk_z('(b) par b_fin', fp%b_fin)
      call chk_arr_z('(b) par ang (whole family)', fp%ang)
      call chk_arch_z('(b) par q_fin', fp%q_fin, 6)
      call chk_arch_z('(b) par p_fin', fp%p_fin, 6)
      call chk_r('(b) par e_cm', fp%e_cm, 0.001d0/e_conv)
      ! --- perpendicular anchors ---
      call chk_z('(b) perp e_rel [all-tangential limit]', fq%e_rel)
      call chk_r('(b) perp e_rel_sq', fq%e_rel_sq, 0.024d0/e_conv)
      call chk_r('(b) perp oam(2)', fq%oam(2), 2.4d0/hbar_code)
      call chk_z('(b) perp oam(1)', fq%oam(1))
      call chk_z('(b) perp oam(3)', fq%oam(3))
      call chk_r('(b) perp oam_mag', fq%oam_mag, 2.4d0/hbar_code)
      call chk_r('(b) perp b_fin [the full separation]', fq%b_fin, 10.0d0)
      call chk_arr_z('(b) perp e_int_a', fq%e_int_a)
      call chk_arr_z('(b) perp e_int_b', fq%e_int_b)
      call chk_arch_z('(b) perp q_fin', fq%q_fin, 6)
      call chk_arch_z('(b) perp p_fin', fq%p_fin, 6)
      call chk_r('(b) perp e_cm', fq%e_cm, 0.016d0/e_conv)
      if (nfail == nfail0) then
         print '(a)', '(b): PASS (parallel: e_rel = e_rel_sq, oam/b_fin/allskip fields'// &
            ' bitwise-zero; perpendicular: e_rel bitwise-zero, e_rel_sq/oam/b_fin'// &
            ' maximal - the two splitting limits the classification consumes)'
      end if
      call state_destroy(st)
   end subroutine crit_b

   !------------------------------------------------------------------
   ! crit_c() - (c): the single-fragment degenerate path
   subroutine crit_c()
      type(state_t) :: st
      type(fin_t) :: fin
      real(8) :: cal, e_tot
      nfail0 = nfail
      call wipe_list_atoms()
      allocate (list_atoms%mass(2))
      list_atoms%mass = [ 1.0d0, 2.0d0 ]
      allocate (list_atoms%frag(1))
      list_atoms%frag(1)%nat = 2
      allocate (list_atoms%frag(1)%list(2))
      list_atoms%frag(1)%list = [ 1, 2 ]
      allocate (list_atoms%frag(1)%qz(6))
      list_atoms%frag(1)%qz = [ 0.9d0, 0.0d0, 0.0d0, -0.45d0, 0.0d0, 0.0d0 ]
      list_atoms%frag(1)%mass = 3.0d0
      ! the conversion recalibrated for the 2-atom list_atoms (same construction)
      call state_create(st, 2)
      st%q = [ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 ]
      st%p = 0.0d0
      call force_eval(st, e_tot)
      cal = e_tot*2.0d0/kk
      call state_destroy(st)
      if (cal <= 0.0d0) then
         call fail('(c) measured conversion constant non-positive: '//d2s(cal))
         return
      end if
      call state_create(st, 2)
      st%q = [ 1.0d0, 0.0d0, 0.0d0, -0.5d0, 0.0d0, 0.0d0 ]
      st%p = [ 0.05d0, 0.4d0, 0.0d0, 0.1d0, -0.4d0, 0.0d0 ]
      st%t = 2.5d0
      call fin_eval(st, fin)
      call chk_r('(c) j_a(3)', fin%j_a(3), 0.6d0/hbar_code)
      call chk_r('(c) j_mag_a', fin%j_mag_a, 0.6d0/hbar_code)
      call chk_z('(c) j_a(1)', fin%j_a(1))
      call chk_z('(c) j_a(2)', fin%j_a(2))
      call chk_r('(c) e_rot_a', fin%e_rot_a, 0.12d0/e_conv)
      call chk_r('(c) e_int_a(1)', fin%e_int_a(1), 0.12d0/e_conv)
      call chk_r('(c) e_int_a(2) [buffered: no other fragment to place]', &
                 fin%e_int_a(2), cal*(kk/2.0d0)*1.25d0)
      call chk_r('(c) e_int_a(3)', fin%e_int_a(3), 0.12d0/e_conv + cal*(kk/2.0d0)*1.25d0)
      call chk_arr_z('(c) e_int_b [no second fragment]', fin%e_int_b)
      call chk_arr_z('(c) j_b', fin%j_b)
      call chk_z('(c) j_mag_b', fin%j_mag_b)
      call chk_z('(c) e_rot_b', fin%e_rot_b)
      call chk_z('(c) e_rel', fin%e_rel)
      call chk_z('(c) e_rel_sq', fin%e_rel_sq)
      call chk_arr_z('(c) oam', fin%oam)
      call chk_z('(c) oam_mag', fin%oam_mag)
      call chk_z('(c) b_fin', fin%b_fin)
      call chk_arr_z('(c) ang (whole family)', fin%ang)
      call chk_r('(c) e_cm [whole-state quantity]', fin%e_cm, 0.00375d0/e_conv)
      call chk_r('(c) t_life', fin%t_life, 25.0d0)
      call chk_archive('(c)', fin, &
                       [ 1.0d0, 0.0d0, 0.0d0, -0.5d0, 0.0d0, 0.0d0 ], &
                       [ 0.0d0, 0.4d0, 0.0d0, 0.0d0, -0.4d0, 0.0d0 ])
      if (nfail == nfail0) then
         print '(a)', '(c): PASS (fragment-A fields anchored; every pair/B field'// &
            ' bitwise-zero; e_cm and t_life still computed - whole-state quantities)'
      end if
      call state_destroy(st)
   end subroutine crit_c

   !------------------------------------------------------------------
   ! crit_d() - (d): the classification channel end to end (the 'cls' child,
   !             cwd = build/fs_cls/, E_SPLIT = 0.5 buffered in the parameter
   !             receiver between the two limit values)
   subroutine crit_d()
      call make_dir('build/fs_cls')
      call expect_ok_child_in('cls', 'cls: scattering', 'cls: capture', &
         '(d) classification channel end to end (real fin values consumed)')
      call wipe_dir('build/fs_cls')
   end subroutine crit_d

   !------------------------------------------------------------------
   ! crit_e() - (e): the named abort faces (children)
   subroutine crit_e()
      call expect_loudstop_at('fin_eval', 'eun', 'atom list', &
         '(e) unassembled-list_atoms loud-stop')
      call expect_loudstop_at('fin_eval', 'esz', 'mismatch', &
         '(e) state/list_atoms size-mismatch loud-stop')
      call expect_loudstop_at('fin_eval', 'e3f', 'fragment', &
         '(e) three-fragment loud-stop')
   end subroutine crit_e

   !------------------------------------------------------------------
   ! the children (self-respawn modes; each terminates itself)
   !------------------------------------------------------------------

   ! child_cls() - the (d) child: bind the reference container's judgment
   !               through the public slot, classify BOTH limit states
   subroutine child_cls()
      type(state_t) :: st
      type(fin_t) :: fin
      type(reference_params_t) :: pload            ! the buffered parameter set
      character(len=32) :: nm
      call adiabatic_register()
      call method_bind('adiabatic')
      call container_bind_pes(reference_pes)
      pload%e_split = 0.5d0
      pload%e_split_found = .true.
      call reference_load(pload)
      call container_bind_classify(reference_classify)
      call wipe_list_atoms()
      allocate (list_atoms%mass(2))
      list_atoms%mass = [ 2.0d0, 3.0d0 ]
      allocate (list_atoms%frag(2))
      list_atoms%frag(1)%nat = 1
      allocate (list_atoms%frag(1)%list(1))
      list_atoms%frag(1)%list = [ 1 ]
      allocate (list_atoms%frag(1)%qz(3))
      list_atoms%frag(1)%qz = 0.0d0
      list_atoms%frag(1)%mass = 2.0d0
      list_atoms%frag(2)%nat = 1
      allocate (list_atoms%frag(2)%list(1))
      list_atoms%frag(2)%list = [ 2 ]
      allocate (list_atoms%frag(2)%qz(3))
      list_atoms%frag(2)%qz = 0.0d0
      list_atoms%frag(2)%mass = 3.0d0
      call state_create(st, 2)
      st%q = [ 0.0d0, 0.0d0, 5.0d0, 0.0d0, 0.0d0, -5.0d0 ]
      st%t = 1.0d0
      ! the parallel state: e_rel ~ 1.29 kcal/mol > 0.5 -> scattering
      st%p = [ 0.0d0, 0.0d0, 0.4d0, 0.0d0, 0.0d0, -0.3d0 ]
      call fin_eval(st, fin)
      nm = container_classify(fin)
      if (trim(nm) /= 'scattering') then
         print '(a)', 'FAIL cls child: parallel state classified "'//trim(nm)// &
            '" (want scattering - e_rel above the split)'
         stop 1
      end if
      print '(a)', 'cls: scattering'
      ! the perpendicular state: e_rel = 0 -> capture
      st%p = [ 0.4d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 ]
      call fin_eval(st, fin)
      nm = container_classify(fin)
      if (trim(nm) /= 'capture') then
         print '(a)', 'FAIL cls child: perpendicular state classified "'//trim(nm)// &
            '" (want capture - e_rel below the split)'
         stop 1
      end if
      print '(a)', 'cls: capture'
      stop 0
   end subroutine child_cls

   ! child_eun() - the (e) child: fin_eval before any list_atoms assembly
   subroutine child_eun()
      type(state_t) :: st
      type(fin_t) :: fin
      call state_create(st, 1)
      call fin_eval(st, fin)        ! must named abort (list_atoms never loaded)
      call state_destroy(st)
   end subroutine child_eun

   ! child_esz() - the (e) child: a 3-atom state against a 4-atom list_atoms
   subroutine child_esz()
      type(state_t) :: st
      type(fin_t) :: fin
      call fill_two_frag_list_atoms()
      call state_create(st, 3)      ! 9 dofs vs the atom list's 12
      call fin_eval(st, fin)        ! must named abort (size mismatch)
      call state_destroy(st)
   end subroutine child_esz

   ! child_e3f() - the (e) child: a three-fragment list_atoms
   subroutine child_e3f()
      type(state_t) :: st
      type(fin_t) :: fin
      integer :: i
      call wipe_list_atoms()
      allocate (list_atoms%mass(3))
      list_atoms%mass = [ 1.0d0, 2.0d0, 3.0d0 ]
      allocate (list_atoms%frag(3))
      do i = 1, 3
         list_atoms%frag(i)%nat = 1
         allocate (list_atoms%frag(i)%list(1))
         list_atoms%frag(i)%list = [ i ]
         allocate (list_atoms%frag(i)%qz(3))
         list_atoms%frag(i)%qz = 0.0d0
         list_atoms%frag(i)%mass = list_atoms%mass(i)
      end do
      call state_create(st, 3)
      call fin_eval(st, fin)        ! must named abort (fragment count)
      call state_destroy(st)
   end subroutine child_e3f

   !------------------------------------------------------------------
   ! assertion helpers (watermark discipline: chk_* record failures only)
   !------------------------------------------------------------------

   ! chk_r - relative scalar check, tol 1e-10 with floor 1.0 (header derivation)
   subroutine chk_r(label, got, want)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got, want
      if (abs(got - want) > 1.0d-10*max(abs(want), 1.0d0)) then
         call fail(trim(label)//': got '//d2s(got)//', want '//d2s(want))
      end if
   end subroutine chk_r

   ! chk_deg - degree-valued check, abs 1e-9 (header derivation)
   subroutine chk_deg(label, got, want)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got, want
      if (abs(got - want) > 1.0d-9) then
         call fail(trim(label)//': got '//d2s(got)//' deg, want '//d2s(want)//' deg')
      end if
   end subroutine chk_deg

   ! chk_z - bitwise-zero scalar check (skip-by-construction contract)
   subroutine chk_z(label, got)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got
      if (got /= 0.0d0) then
         call fail(trim(label)//': must be bitwise zero, got '//d2s(got))
      end if
   end subroutine chk_z

   ! chk_arr_z - bitwise-zero check for a fixed-shape array member (the skip
   !              contract: any nonzero entry means a wrong branch executed)
   subroutine chk_arr_z(label, got)
      character(len=*), intent(in) :: label
      real(8), intent(in) :: got(:)
      integer :: i
      do i = 1, size(got)
         if (got(i) /= 0.0d0) then
            call fail(trim(label)//': entry '//i2s(i)//' must be bitwise zero, got '// &
                      d2s(got(i)))
            return
         end if
      end do
   end subroutine chk_arr_z

   ! chk_archive - the COM-frame archive vs its anchors (allocation-guarded:
   !                a missing archive FAILs, never crashes)
   subroutine chk_archive(label, fin, q_want, p_want)
      character(len=*), intent(in) :: label
      type(fin_t), intent(in) :: fin
      real(8), intent(in) :: q_want(:), p_want(:)
      integer :: i
      if (.not. allocated(fin%q_fin) .or. .not. allocated(fin%p_fin)) then
         call fail(trim(label)//' archive: q_fin/p_fin not allocated')
         return
      end if
      if (size(fin%q_fin) /= size(q_want) .or. size(fin%p_fin) /= size(p_want)) then
         call fail(trim(label)//' archive: wrong length '//i2s(size(fin%q_fin)))
         return
      end if
      do i = 1, size(q_want)
         if (abs(fin%q_fin(i) - q_want(i)) > 1.0d-10*max(abs(q_want(i)), 1.0d0)) then
            call fail(trim(label)//' archive q_fin('//i2s(i)//'): got '// &
                      d2s(fin%q_fin(i))//', want '//d2s(q_want(i)))
            return
         end if
         if (abs(fin%p_fin(i) - p_want(i)) > 1.0d-10*max(abs(p_want(i)), 1.0d0)) then
            call fail(trim(label)//' archive p_fin('//i2s(i)//'): got '// &
                      d2s(fin%p_fin(i))//', want '//d2s(p_want(i)))
            return
         end if
      end do
   end subroutine chk_archive

   ! chk_arch_z - the zero archive (one-atom fragments / the degenerate path):
   !              allocated at the full state length with every entry
   !              EXACTLY zero (allocation-guarded: missing = FAIL, no crash)
   subroutine chk_arch_z(label, arr, n)
      character(len=*), intent(in) :: label
      real(8), allocatable, intent(in) :: arr(:)
      integer, intent(in) :: n
      integer :: i
      if (.not. allocated(arr)) then
         call fail(trim(label)//': the archive array is not allocated')
         return
      end if
      if (size(arr) /= n) then
         call fail(trim(label)//': the archive array has wrong length '//i2s(size(arr)))
         return
      end if
      do i = 1, n
         if (arr(i) /= 0.0d0) then
            call fail(trim(label)//': archive entry '//i2s(i)//' must be bitwise'// &
                      ' zero, got '//d2s(arr(i)))
            return
         end if
      end do
   end subroutine chk_arch_z

   !------------------------------------------------------------------
   ! small helpers (the sibling checks' discipline)
   !------------------------------------------------------------------

   ! wipe_list_atoms() - drop any previously filled list_atoms ahead of a re-fill
   !                 (allocating an already-allocated array is a runtime error)
   subroutine wipe_list_atoms()
      integer :: i
      if (allocated(list_atoms%mass)) deallocate (list_atoms%mass)
      if (allocated(list_atoms%frag)) then
         do i = 1, size(list_atoms%frag)
            if (allocated(list_atoms%frag(i)%list)) deallocate (list_atoms%frag(i)%list)
            if (allocated(list_atoms%frag(i)%qz)) deallocate (list_atoms%frag(i)%qz)
         end do
         deallocate (list_atoms%frag)
      end if
      list_atoms%i_proj = 0
      list_atoms%i_targ = 0
   end subroutine wipe_list_atoms

   ! fill_two_frag_list_atoms() - the (a)-shaped 4-atom list_atoms (the 'esz' child)
   subroutine fill_two_frag_list_atoms()
      call wipe_list_atoms()
      allocate (list_atoms%mass(4))
      list_atoms%mass = [ 1.0d0, 2.0d0, 3.0d0, 4.0d0 ]
      allocate (list_atoms%frag(2))
      list_atoms%frag(1)%nat = 2
      allocate (list_atoms%frag(1)%list(2))
      list_atoms%frag(1)%list = [ 1, 2 ]
      allocate (list_atoms%frag(1)%qz(6))
      list_atoms%frag(1)%qz = [ 0.9d0, 0.0d0, 0.0d0, -0.45d0, 0.0d0, 0.0d0 ]
      list_atoms%frag(1)%mass = 3.0d0
      list_atoms%frag(2)%nat = 2
      allocate (list_atoms%frag(2)%list(2))
      list_atoms%frag(2)%list = [ 3, 4 ]
      allocate (list_atoms%frag(2)%qz(6))
      list_atoms%frag(2)%qz = [ 0.0d0, 1.0d0, 0.0d0, 0.0d0, -0.75d0, 0.0d0 ]
      list_atoms%frag(2)%mass = 7.0d0
   end subroutine fill_two_frag_list_atoms

   ! expect_loudstop_at - spawn one guard child and assert the stop: nonzero
   !                      exit, stderr names the authority and the key
   subroutine expect_loudstop_at(who, cmode, key, label)
      character(len=*), intent(in) :: who, cmode, key, label
      character(len=512) :: cmd
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_final_state.exe '//trim(cmode)// &
            ' 2>build'//sep//'fs_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp('build/fs_err.txt', buf)
      if (stat == 0) then
         call fail(trim(label)//': expected a loud-stop, the child exited 0')
      else if (index(buf, trim(who)) == 0) then
         call fail(trim(label)//': stderr does not name '//trim(who)//' - "'// &
                   trim(buf)//'"')
      else if (index(buf, trim(key)) == 0) then
         call fail(trim(label)//': stderr missing the key "'//trim(key)//'" - "'// &
                   trim(buf)//'"')
      else
         print '(a)', trim(label)//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop_at

   ! expect_ok_child_in - spawn one success child with cwd = build/fs_cls/
   !                      (the parameter file is run-directory material) and
   !                      assert exit 0 + both expected stdout markers
   subroutine expect_ok_child_in(cmode, marker1, marker2, label)
      character(len=*), intent(in) :: cmode, marker1, marker2, label
      character(len=512) :: cmd
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'cd build'//sep//'fs_cls && ..'//sep//'check_final_state.exe '// &
            trim(cmode)//' >ok.txt 2>&1'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      call slurp('build/fs_cls/ok.txt', buf)
      if (stat /= 0) then
         call fail(trim(label)//': the child exited nonzero - '//trim(buf))
      else if (index(buf, trim(marker1)) == 0 .or. index(buf, trim(marker2)) == 0) then
         call fail(trim(label)//': stdout missing "'//trim(marker1)//'" / "'// &
                  trim(marker2)//'" - '//trim(buf))
      else
         print '(a)', trim(label)//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_ok_child_in

   ! write_lines - write one text file, one array element per record
   subroutine write_lines(fname, lines)
      character(len=*), intent(in) :: fname, lines(:)
      integer :: v, j
      open (newunit=v, file=fname, status='replace', action='write')
      do j = 1, size(lines)
         write (v, '(a)') trim(lines(j))
      end do
      close (v)
   end subroutine write_lines

   ! slurp - one text file into a single-line buffer (missing file = empty)
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

   ! make_dir / wipe_dir - host-shaped folder helpers
   subroutine make_dir(path)
      character(len=*), intent(in) :: path
      if (windows_host()) then
         call execute_command_line('mkdir "'//trim(path)//'" 2>nul', wait=.true.)
      else
         call execute_command_line('mkdir -p "'//trim(path)//'" 2>/dev/null', wait=.true.)
      end if
   end subroutine make_dir

   subroutine wipe_dir(path)
      character(len=*), intent(in) :: path
      if (windows_host()) then
         call execute_command_line('rd /s /q "'//trim(path)//'" >nul 2>nul', wait=.true.)
      else
         call execute_command_line('rm -rf "'//trim(path)//'" 2>/dev/null', wait=.true.)
      end if
   end subroutine wipe_dir

   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      print '(a)', 'FAIL '//trim(msg)
      nfail = nfail + 1
   end subroutine fail

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
end program check_final_state
