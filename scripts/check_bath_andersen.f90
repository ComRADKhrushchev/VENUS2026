!=====================================================================
! check_bath_andersen.f90 - bath species criteria program (criteria
!                              defined FIRST in this header; executed below)
! Compile (from src/, linked against the build objects):
!   gfortran -O2 -Jbuild/gfortran -Ibuild/gfortran -o build/gfortran/check_bath_andersen.exe \
!            ../scripts/check_bath_andersen.f90 \
!            $(find build/gfortran -name '*.o' ! -name 'driver.o')
! Criteria:
!   (a) width table: sigma(i)^2 = m(i)*kb*T_bath for three distinct masses
!       (the whole-system scope - every atom, no fragment binding)
!   (b) certain-collision sweep (nu*dt >= 1) = the exact Maxwell reset: over
!       20000 sweeps on a one-atom state, the standardized component
!       z = p_x/sigma holds mean 0 and variance 1 within the 4/sqrt(N) and
!       4*sqrt(2/N) bands (the fresh-velocity-per-hit contract)
!   (c) collision probability: nu*dt = 0.1, sentinel momenta, one sweep per
!       trial over 4000 trials - the hit fraction lies within the 4-sigma
!       binomial band of 0.1; every HIT rewrites all three components and
!       every MISS keeps its sentinel bitwise
!   (d)/(e) inert loops: BATH = none, and ANDERSEN with N_BATH = 0, both
!       leave the state bitwise untouched (the archived-stream contract)
!   (f) named aborts (children): t_bath <= 0; nu_coll <= 0; bath_setup with
!       N_BATH > 0 and DT_BATH = 0 -> stderr names the caller and the field
!   (g) freeze statute: an armed bath on a fragment-only system (ownerless
!       set empty) leaves the state bitwise untouched
!   (h) separation prober: a decaying synthetic force - the ladder stops at
!       the first separation whose environment force stays under tolerance
!   (i) thermal footprint: a mixed system (one fragment + ownerless atoms) -
!       the fragment rows are restored bitwise while the ownerless atoms KEEP
!       their bathed q/p (the thermalized snapshot is the trajectory start)
!=====================================================================
program check_bath_andersen
   use bath_andersen, only: andersen_params_t, andersen_init, andersen_sweep, &
                              andersen_sigma, andersen_nat, andersen_t_bath, andersen_nu
   use geom_init,      only: geom_sep_probe
   use bath,            only: bath_run, bath_setup, bath_armed
   use bath_andersen, only: andersen_register => member_register
   use config,          only: reactants
   use config_atoms,    only: list_atoms
   use state,           only: state_t, state_create, state_destroy
   use rng,             only: rng_init
   use consts,          only: kb_act
   implicit none
   type(state_t) :: sta
   type(andersen_params_t) :: p
   real(8), parameter :: sentinel = -777.0d0
   real(8) :: z, s1, s2, expo, band, frac, d_out
   real(8) :: q2p(6)
   integer :: nfail = 0, nfail0, i, t, nhit
   logical :: mask3(3) = .true.
   character(len=256) :: mode
   integer, parameter :: ndraw = 20000, mtrial = 4000

   ! ---------- child dispatch: one bad case per spawn ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('bad-t')                     ! (f1): non-positive bath temperature
         call load_roster()
         p%t_bath = 0.0d0
         p%nu_coll = 0.1d0
         call andersen_init(p)
      case ('bad-nu')                    ! (f2): non-positive collision frequency
         call load_roster()
         p%t_bath = 300.0d0
         p%nu_coll = 0.0d0
         call andersen_init(p)
      case ('bad-pes')                   ! (f4): the interaction never decays
         q2p = [ 0.0d0, 0.0d0, 0.0d0, 3.0d0, 0.0d0, 0.0d0 ]
         call geom_sep_probe(flat_force, [ 1.0d0, 1.0d0 ], [ 1 ], q2p, 8.0d0, 1.0d-6, d_out)
      case ('bad-dt')                    ! (f3): N_BATH > 0 with DT_BATH = 0
         call andersen_register()
         call load_roster()
         p%t_bath = 300.0d0
         p%nu_coll = 0.1d0
         call andersen_init(p)
         reactants%bath = 2
         reactants%n_bath = 5
         reactants%dt_bath = 0.0d0
         call bath_setup()
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the species returned - bug)'
      stop 0
   end if

   call rng_init(20260922)

   ! ---------- (a) width table ----------
   nfail0 = nfail
   call load_roster()
   p%t_bath = 300.0d0
   p%nu_coll = 0.137d0
   call andersen_init(p)
   if (andersen_nat() /= 3 .or. andersen_t_bath() /= 300.0d0 .or. &
       andersen_nu() /= 0.137d0) then
      print '(a)', 'FAIL (a) width table: nat or parameter getters wrong'
      nfail = nfail + 1
   else
      do i = 1, 3
         expo = andersen_sigma(i)**2/(list_atoms%mass(i)*kb_act*300.0d0)
         if (abs(expo - 1.0d0) > 1.0d-12) then
            print '(a,i0,a)', 'FAIL (a) width table: atom ', i, ' sigma^2 != m*kb*T'
            nfail = nfail + 1
            exit
         end if
      end do
      if (nfail == nfail0) print '(a)', &
         '(a) width table: PASS (sigma^2 = m*kb*T for all three atoms; whole-system scope)'
   end if

   ! ---------- (b) certain-collision sweep = the exact Maxwell reset ----------
   nfail0 = nfail
   p%nu_coll = 10.0d0                            ! nu*dt = 10 >= 1: every atom collides
   call andersen_init(p)
   call state_create(sta, 3)
   s1 = 0.0d0
   s2 = 0.0d0
   do t = 1, ndraw
      call andersen_sweep(sta, 1.0d0, mask3)
      z = sta%p(1)/andersen_sigma(1)
      s1 = s1 + z
      s2 = s2 + z*z
   end do
   if (abs(s1/ndraw) > 4.0d0/sqrt(dble(ndraw)) .or. &
       abs(s2/ndraw - 1.0d0) > 4.0d0*sqrt(2.0d0/dble(ndraw))) then
      print '(a)', 'FAIL (b) Maxwell reset: standardized moments outside the 4-sigma bands'
      nfail = nfail + 1
   else
      print '(a)', '(b) Maxwell reset: PASS (z mean/var within 4-sigma bands over 20000 resets)'
   end if
   call state_destroy(sta)

   ! ---------- (c) collision probability + untouched contract ----------
   nfail0 = nfail
   p%nu_coll = 0.1d0
   call andersen_init(p)
   call state_create(sta, 3)
   nhit = 0
   do t = 1, mtrial
      sta%p(1:3) = sentinel
      call andersen_sweep(sta, 1.0d0, mask3)   ! nu*dt = 0.1
      if (sta%p(1) /= sentinel .or. sta%p(2) /= sentinel .or. sta%p(3) /= sentinel) then
         nhit = nhit + 1                          ! a hit rewrites ALL three components
         if (sta%p(1) == sentinel .or. sta%p(2) == sentinel .or. sta%p(3) == sentinel) then
            print '(a)', 'FAIL (c) collision trial: a hit left a component untouched'
            nfail = nfail + 1
            exit
         end if
      end if
   end do
   frac = dble(nhit)/dble(mtrial)
   band = 4.0d0*sqrt(0.1d0*0.9d0/dble(mtrial))
   if (nfail == nfail0 .and. abs(frac - 0.1d0) > band) then
      print '(a)', 'FAIL (c) collision probability: hit fraction outside the 4-sigma band'
      nfail = nfail + 1
   else if (nfail == nfail0) then
      print '(a)', '(c) collision probability: PASS (hit fraction in band; misses bitwise)'
   end if
   call state_destroy(sta)

   ! ---------- (d)/(e) inert loops: bitwise no-ops ----------
   nfail0 = nfail
   call andersen_register()                ! the species row the provision binds
   call state_create(sta, 3)
   sta%q = sentinel
   sta%p = sentinel
   reactants%bath = 1                            ! NONE
   call bath_setup()
   call bath_run(sta)
   reactants%bath = 2                            ! ANDERSEN, zero steps
   reactants%n_bath = 0
   reactants%dt_bath = 0.0d0
   call bath_setup()
   call bath_run(sta)
   if (any(sta%q /= sentinel) .or. any(sta%p /= sentinel) .or. bath_armed()) then
      print '(a)', 'FAIL (d/e) inert bath: an unarmed or zero-step loop touched the state'
      nfail = nfail + 1
   else
      print '(a)', '(d/e) inert bath: PASS (BATH=none and N_BATH=0 both bitwise no-ops)'
   end if
   call state_destroy(sta)

   ! ---------- (g) freeze statute: armed bath on a fragment-only system is
   ! a bitwise no-op (ownerless set empty - no seeding, no sweep draws; the
   ! fragments are save/placed/propagated/restored bitwise) ----------
   nfail0 = nfail
   call load_roster()
   allocate (list_atoms%frag(1))               ! all three atoms ONE fragment:
   list_atoms%frag(1)%nat = 3                  ! member-owned, frozen by statute
   allocate (list_atoms%frag(1)%list(3), list_atoms%frag(1)%qz(9))
   list_atoms%frag(1)%list = [ 1, 2, 3 ]
   list_atoms%frag(1)%qz = 0.0d0
   p%t_bath = 300.0d0
   p%nu_coll = 0.1d0
   call andersen_init(p)
   call bind_zero_force()   ! the species row is registered by the (d/e) leg
   call state_create(sta, 3)
   sta%q = sentinel
   sta%p = 2.0d0*sentinel
   reactants%bath = 2
   reactants%n_bath = 3
   reactants%dt_bath = 0.05d0
   reactants%r_sep = 8.0d0
   call bath_setup()
   call bath_run(sta)
   if (any(sta%q /= sentinel) .or. any(sta%p /= 2.0d0*sentinel) .or. .not. bath_armed()) then
      print '(a)', 'FAIL (g) freeze statute: the armed loop left a footprint on a member-owned fragment'
      nfail = nfail + 1
   else
      print '(a)', '(g) freeze statute: PASS (armed bath on a fragment-only system: bitwise no-op)'
   end if
   call state_destroy(sta)

   ! ---------- (i) thermal footprint: a mixed system (one nat=2 fragment,
   ! atom 3 ownerless) - the fragment rows are restored bitwise while the
   ! ownerless atom KEEPS its bathed q/p. nu*dt >= 1: the ownerless atom
   ! collides every sweep, so a correct loop ends with nonzero bathed
   ! momenta and a propagated geometry; a full-array restore instead
   ! discards the bath (ownerless back at the q0 seed, zero momenta).
   ! The fragment carries zero net momentum (the projection reference -
   ! a nonzero fragment net would leak into the ownerless rows) ----------
   nfail0 = nfail
   call load_roster()
   list_atoms%q0(7:9) = [ 5.0d0, 0.0d0, 0.0d0 ]  ! the ownerless atom's seed site
   deallocate (list_atoms%frag)
   allocate (list_atoms%frag(1))                ! atoms 1-2 ONE fragment, atom 3 ownerless
   list_atoms%frag(1)%nat = 2
   allocate (list_atoms%frag(1)%list(2), list_atoms%frag(1)%qz(6))
   list_atoms%frag(1)%list = [ 1, 2 ]
   list_atoms%frag(1)%qz = 0.0d0
   p%t_bath = 300.0d0
   p%nu_coll = 100.0d0                          ! nu*dt = 5 >= 1: certain collision every sweep
   call andersen_init(p)
   call state_create(sta, 3)
   sta%q = sentinel
   sta%p(1:3) = [ 1.0d0, 2.0d0, 3.0d0 ]         ! zero net fragment momentum
   sta%p(4:6) = [ -1.0d0, -2.0d0, -3.0d0 ]
   sta%p(7:9) = 3.0d0*sentinel                  ! overwritten by the ownerless seeding
   reactants%bath = 2
   reactants%n_bath = 3
   reactants%dt_bath = 0.05d0
   call bath_setup()
   call bath_run(sta)
   if (any(sta%q(1:6) /= sentinel) .or. any(sta%p(1:3) /= [ 1.0d0, 2.0d0, 3.0d0 ]) .or. &
       any(sta%p(4:6) /= [ -1.0d0, -2.0d0, -3.0d0 ])) then
      print '(a)', 'FAIL (i) freeze statute: the mixed-system loop left a footprint on the fragment'
      nfail = nfail + 1
   else
      print '(a)', '(i) freeze statute: PASS (fragment rows bitwise restored in the mixed system)'
   end if
   if (all(sta%p(7:9) == 0.0d0)) then
      print '(a)', 'FAIL (i) thermal footprint: the ownerless momenta are zero'// &
         ' (the bathed state was discarded by the unfreeze)'
      nfail = nfail + 1
   else
      print '(a)', '(i) thermal footprint: PASS (the ownerless atom keeps its bathed momenta)'
   end if
   if (all(sta%q(7:9) == list_atoms%q0(7:9))) then
      print '(a)', 'FAIL (i) thermal footprint: the ownerless geometry sits at the seed'// &
         ' (the bathed positions were discarded by the unfreeze)'
      nfail = nfail + 1
   else
      print '(a)', '(i) thermal footprint: PASS (the ownerless atom keeps its bathed geometry)'
   end if
   if (.not. bath_armed()) then
      print '(a)', 'FAIL (i) arming: the provision did not arm'
      nfail = nfail + 1
   end if
   call state_destroy(sta)

   ! ---------- (h) separation prober: decaying synthetic force, ladder from 8 ----------
   nfail0 = nfail
   d_out = -1.0d0
   ! atom 1 at origin (the fragment), atom 2 at x=3 (the environment); the
   ! callback puts a force F2x = -exp(-(x1 - x2)) on atom 2 only
   q2p = [ 0.0d0, 0.0d0, 0.0d0, 3.0d0, 0.0d0, 0.0d0 ]
   call geom_sep_probe(decay_force, [ 1.0d0, 1.0d0 ], [ 1 ], q2p, 8.0d0, 1.0d-6, d_out)
   ! ladder 8 -> 12 (exp(-9)=1.2e-4 > tol) -> 18 (exp(-15)=3e-7 <= tol): d*=18,
   ! fragment COM placed at [18, 0, 0]
   if (abs(d_out - 18.0d0) > 1.0d-12 .or. abs(q2p(1) - 18.0d0) > 1.0d-12 .or. &
       abs(q2p(2)) > 1.0d-12 .or. abs(q2p(3)) > 1.0d-12) then
      print '(a)', 'FAIL (h) separation prober: ladder distance or placement wrong'
      nfail = nfail + 1
   else
      print '(a)', '(h) separation prober: PASS (decaying force: d*=18 A, COM at [18,0,0])'
   end if

   ! ---------- (f) named aborts (children) ----------
   call expect_loudstop('bad-t',  'andersen_init', 't_bath',  '(f1) bad temperature loud-stop')
   call expect_loudstop('bad-nu', 'andersen_init', 'nu_coll', '(f2) bad frequency loud-stop')
   call expect_loudstop('bad-dt', 'bath_setup',     'dt_bath', '(f3) zero-step-size loud-stop')
   call expect_loudstop('bad-pes', 'geom_sep_probe', 'does not decay', '(f4) non-decaying PES loud-stop')

   print '(a,i0,a)', '=== check_bath_andersen: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! decaying synthetic force: F2x = -exp(-(x1 - x2)) on atom 2 (a pair
   ! interaction that vanishes with separation; atom 1 force-free)
   subroutine decay_force(q, f)
      real(8), intent(in) :: q(:)
      real(8), intent(out) :: f(:)
      f = 0.0d0
      f(4) = -exp(-abs(q(1) - q(4)))   ! bounded both ways (the probe's
      ! parking reference sits 1000 A aside - a one-sided exp would overflow)
   end subroutine decay_force

   ! non-decaying synthetic force: F2x = -1/(1 + (x1 - x2)) - fragment-
   ! dependent but only algebraically decaying (still ~5e-5 at the ladder
   ! cap ~17700 A), never under tol: the suspect PES)
   subroutine flat_force(q, f)
      real(8), intent(in) :: q(:)
      real(8), intent(out) :: f(:)
      f = 0.0d0
      f(4) = -1.0d0/(1.0d0 + (q(1) - q(4)))
   end subroutine flat_force

   ! bind a zero-force container so bath_run's propagation is free flight
   subroutine bind_zero_force()
      use force_interface, only: container_bind_pes, container_init
      use elec_interface, only: method_reg
      use propagator, only: prop_register => member_register
      use config, only: propagator
      call method_reg('adiabatic', noop_prop, noop_match)
      call container_bind_pes(zero_pes)
      call container_init()
      call prop_register()
      propagator%id = 1
      propagator%order = 1
   end subroutine bind_zero_force

   subroutine zero_pes(q, mass, g, v)
      real(8), intent(in) :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      g = 0.0d0
      v = 0.0d0
   end subroutine zero_pes

   subroutine noop_prop(sta, dt)
      use state, only: state_t
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine noop_prop

   subroutine noop_match(sta)
      use state, only: state_t
      type(state_t), intent(inout) :: sta
   end subroutine noop_match

   ! three-atom roster by direct field fill (masses only - the species scope)
   subroutine load_roster()
      if (allocated(list_atoms%mass)) deallocate (list_atoms%mass)
      allocate (list_atoms%mass(3))
      list_atoms%mass = [ 1.008d0, 12.011d0, 35.45d0 ]
      if (allocated(list_atoms%q0)) deallocate (list_atoms%q0)
      allocate (list_atoms%q0(9))
      list_atoms%q0 = 0.0d0
   end subroutine load_roster

   ! spawn the child for one bad case and assert the named abort
   subroutine expect_loudstop(cmode, who, key1, label)
      character(len=*), intent(in) :: cmode, who, key1, label
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'gfortran'//sep//'check_bath_andersen.exe '//trim(cmode)// &
            ' 2>build'//sep//'gfortran'//sep//'bath_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build'//sep//'gfortran'//sep//'bath_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) line
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(line)
         end do
         close (v)
      end if
      if (stat == 0) then
         print '(a)', 'FAIL '//label//': expected a loud-stop, the child exited 0'
         nfail = nfail + 1
      else if (index(buf, who) == 0) then
         print '(a)', 'FAIL '//label//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (index(buf, trim(key1)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key1)//'"'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop

   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_bath_andersen
