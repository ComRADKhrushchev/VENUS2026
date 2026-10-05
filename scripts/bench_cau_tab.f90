!=====================================================================
! bench_cau_tab.f90 - the preliminary validation bench: cau_tab vs
!   cau_peri in ONE process (module names do not collide)
! Design:
!   Two legs, argv selects the cau_tab table mode:
!     identity - the table is generated at run time by sampling the
!                cau_peri analytic chain at the three sites (E0 from
!                calc_POT_spec, h/U/J/VCP from cau_peri_hubbard_val,
!                h tilde = the average of h(1)/h(2)) - SAME physics,
!                different machinery -> the machinery benchmark
!     rst | current - the sparse channel tables from data/c_au/tables
!                -> the data-lineage influence (cau_peri analytic rides
!                along as the old-container reference)
!   Trajectory model: engine layout [S, Y, C], S at the origin, ghost Y
!   at (0,0,0.9), C at (site_x, site_y, z0) with a downward momentum
!   for the requested kinetic energy; velocity-Verlet dt = 0.5 fs with
!   the TDHF prop member when armed (the check_cau_tab vstep pattern);
!   fixed step count, no termination judgment (the z range is tracked).
!   Conditions: leg A (identity only) scatters from 4.5 A at 0.76/1.46
!   eV (the historic smoke energies; the analytic well is shallow
!   enough to stay in band); legs with the deep-well lineage tables
!   (rst well ~ -6 eV, current SCF-E0 well ~ -13 eV) use a bound
!   oscillation launched at 3.0 A with 0.15 eV - a scattering launch
!   would drive the inner turning point below the tables' data floors.
!   Metrics: |H_end - H_0| [kcal/mol], the C z range, max |q_tab - q_peri|
!   along the trajectory (leg A), the impurity occupations and e_elec
!   at the end (armed legs). Compile: like check_cau_tab minus the box.
!   Run from scripts/ (needs cau_peri's Analytic_Potential.txt staged).
!=====================================================================
program bench_cau_tab
   use cau_tab_tab,  only: tab_eval, sh_E0, sh_h, sh_U, sh_J, sh_VCP, si_top, si_bri, si_hol
   use cau_tab_pot,  only: cau_tab_load, ct_k => k_surf, ct_kg => k_ghost, ct_q0 => Q0_ghost
   use cau_tab_tdhf, only: cau_tab_tdhf_setup, cau_tab_pes, cau_tab_prop, &
                           cau_tab_col_nup, cau_tab_col_ndn, cau_tab_col_eelec, cau_tab_col_emf
   use cau_peri_pot, only: cau_peri_load, cau_peri_hubbard_val
   use cau_peri_tdhf, only: cau_peri_tdhf_setup, cau_peri_pes, cau_peri_prop, &
                            cau_peri_col_nup, cau_peri_col_ndn, cau_peri_col_eelec
   use surface_potentials, only: parm_nas, calc_pot_spec_and_grad, cp_k => k_surf, cp_kg => k_ghost, cp_q0 => Q0_ghost
   use state,        only: state_t, state_create, state_destroy
   implicit none
   real(8), parameter :: ev_kcal = 23.0605d0, e_conv = 0.04184d0
   real(8), parameter :: a_lat = 2.95d0
   real(8), parameter :: site_xy(3, 2) = reshape([0.0d0, 1.475d0, 1.475d0, &
                                                  0.0d0, 0.0d0, 0.8515917d0], [3, 2])
   real(8) :: au_mass = 196.97d0, c_mass = 12.01d0
   character(len=16) :: mode
   character(len=256) :: tabfile
   character(len=64) :: line
   real(8) :: ks
   real(8) :: mass3(3)
   integer, parameter :: ncond = 2
   character(len=2) :: seeds(ncond)
   real(8) :: ekins(ncond), z0s(ncond)
   integer :: sites(ncond), ic

   call get_command_argument(1, mode)
   if (len_trim(mode) == 0) mode = 'identity'

   ! ---- stage the engines ----
   call copy_file('../systems/inactive/cau_peri/Analytic_Potential.txt', 'Analytic_Potential.txt')
   call cau_peri_load()
   if (command_argument_count() >= 3) then
      call get_command_argument(3, line)
      if (trim(line) == 'nn') then
         call copy_file('../systems/inactive/cau_peri/nn_weights.txt', 'nn_weights.txt')
         call cau_peri_load('nn_weights.txt')            ! the NN branch reference
         print '(a)', 'reference branch: NN (nn_weights.txt)'
      end if
   end if
   if (trim(mode) == 'identity') then
      call write_identity_table('bench_tab.txt')
      tabfile = 'bench_tab.txt'
   else if (trim(mode) == 'rst') then
      tabfile = 'C:/Users/13614/Desktop/MAIN/01_WORK/data/c_au/tables/tab_rst_sparse.txt'
   else if (trim(mode) == 'current') then
      tabfile = 'C:/Users/13614/Desktop/MAIN/01_WORK/data/c_au/tables/tab_current_sparse.txt'
   else
      print '(a)', 'fatal: mode must be identity|rst|current'
      stop 1
   end if
   call cau_tab_load(trim(tabfile))
   ks = 1.0d0
   if (command_argument_count() >= 2) then
      call get_command_argument(2, line)
      read (line, *) ks
   end if
   ct_k = [ks, ks, ks]
   ct_kg = [ks, ks, ks]
   ct_q0 = [0.0d0, 0.0d0, 0.0d0]
   cp_k = [ks, ks, ks]
   cp_kg = [ks, ks, ks]
   cp_q0 = [0.0d0, 0.0d0, 0.0d0]
   if (ks /= 1.0d0) print '(a,f6.2)', 'DIAGNOSTIC spring scale K_SURF=K_GHOST = ', ks
   mass3 = [au_mass, au_mass, c_mass]

   print '(a,a,a,i0,a)', 'bench_cau_tab mode=', trim(mode), '  (both engines in-process; dt=0.5 fs)'

   ! ---- the HISTORICAL production conditions (smoke_cau_peri/real_
   !      conditions/input_{E1D,E3P}_k1.txt, retrieved 2026-10-01): single C
   !      launched from 6.0 A at normal incidence onto the top site
   !      (AIM 0,0); E1D = the 1D seed at 17.526 kcal/mol (0.76 eV), E3P =
   !      the GS(3P) seed at 33.6683 kcal/mol (1.46 eV); dt = 0.5 fs,
   !      2000-step cap, exit when the C rises back past z = 6.0 (Z_TERM);
   !      bath 20 orbitals / 7 eV / 300 K / MSBSTP = 4. ORDER: the adiabatic
   !      batch first - the armed flag is module state and an armed setup
   !      poisons every later adiabatic pass ----
   seeds = ['1D', 'GS']
   ekins = [0.76d0, 1.46d0]
   z0s = [6.0d0, 6.0d0]
   sites = [1, 1]
   do ic = 1, ncond
      call run_pair('u', ekins(ic), seeds(ic), z0s(ic), sites(ic), 3000, .false.)
   end do
   do ic = 1, ncond
      call run_pair('a', ekins(ic), seeds(ic), z0s(ic), sites(ic), 3000, .true.)
   end do

contains

   subroutine run_pair(tag, ekin, seed, z0, site, nstep, armed)
      character(len=*), intent(in) :: tag, seed
      real(8), intent(in) :: ekin, z0
      integer, intent(in) :: site, nstep
      logical, intent(in) :: armed
      type(state_t) :: st1, st2
      real(8) :: mp(9), mt(9), pz, dqmax, dqhalf
      real(8), allocatable :: qh1(:, :), qh2(:, :)
      integer :: k, j2
      pz = -sqrt(2.0d0*c_mass*ekin*e_conv)
      if (armed) then
         call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, seed, 4, 0, 0.05d0)
         call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, seed, 4, 1, 0, 0.05d0)
      end if
      allocate (qh1(9, 0:nstep), qh2(9, 0:nstep))
      ! engine 1: cau_peri analytic
      call state_create(st1, 3, 1)
      st1%q = 0.0d0
      st1%q(6) = 0.9d0
      st1%q(7:9) = [site_xy(site, 1), site_xy(site, 2), z0]
      st1%p = 0.0d0
      st1%p(9) = pz
      st1%t = 0.0d0
      call traj_peri(st1, nstep, armed, mp, qh1)
      ! engine 2: cau_tab (the mode's table)
      call state_create(st2, 3, 1)
      st2%q = 0.0d0
      st2%q(6) = 0.9d0
      st2%q(7:9) = [site_xy(site, 1), site_xy(site, 2), z0]
      st2%p = 0.0d0
      st2%p(9) = pz
      st2%t = 0.0d0
      call traj_tab(st2, nstep, armed, mt, qh2)
      ! the machinery-deviation metric: max |q_tab(t) - q_peri(t)| up to the
      ! earlier exit step, componentwise
      k = min(int(mp(2)), int(mt(2)))
      dqmax = 0.0d0
      dqhalf = 0.0d0
      do j2 = 0, k
         dqmax = max(dqmax, maxval(abs(qh2(:, j2) - qh1(:, j2))))
         if (j2 <= k/2) dqhalf = max(dqhalf, maxval(abs(qh2(:, j2) - qh1(:, j2))))
      end do
      print '(a,a,a,f5.2,a,a,a,i1,a,l1)', ' ', tag, ' E=', ekin, ' eV seed=', seed, &
         ' site=', site, ' armed=', armed
      print '(a,es11.3,a,f7.3,a,f7.3,a,f7.3,a,i0,a,l1)', '   peri : drift[kcal]=', mp(1), &
         '  z in [', mp(4), ',', mp(5), ']  z_end=', mp(6), ' steps=', int(mp(2)), &
         ' out=', mp(3) > 0.5d0
      print '(a,es11.3,a,f7.3,a,f7.3,a,f7.3,a,i0,a,l1)', '   tab  : drift[kcal]=', mt(1), &
         '  z in [', mt(4), ',', mt(5), ']  z_end=', mt(6), ' steps=', int(mt(2)), &
         ' out=', mt(3) > 0.5d0
      if (armed) then
         print '(a,f6.3,a,f6.3,a,f9.3)', '   peri : n_up=', mp(7), ' n_dn=', mp(8), &
            ' e_elec=', mp(9)
         print '(a,f6.3,a,f6.3,a,f9.3)', '   tab  : n_up=', mt(7), ' n_dn=', mt(8), &
            ' e_elec=', mt(9)
      end if
      print '(a,es9.2,a,es9.2)', '   |q_tab-q_peri|: half-traj max=', dqhalf, '  full max=', dqmax
      deallocate (qh1, qh2)
      call state_destroy(st1)
      call state_destroy(st2)
   end subroutine run_pair

   !------------------------------------------------------------------
   ! traj_peri / traj_tab - the two vstep loops; the tab loop also
   !                 measures max|q_tab - q_peri| against the peri
   !                 reference trajectory held in st_ref
   !------------------------------------------------------------------
   subroutine traj_peri(st, nstep, armed, m, qh)
      type(state_t), intent(inout) :: st
      integer, intent(in) :: nstep
      logical, intent(in) :: armed
      real(8), intent(out) :: m(9)
      real(8), intent(out) :: qh(9, 0:nstep)
      real(8) :: gg(9), vv, h0, h1
      logical :: passed
      integer :: k
      call cau_peri_pes(st%q, mass3, gg, vv)
      h0 = tot_e(st, vv)
      m(4) = st%q(9)
      m(5) = st%q(9)
      m(2) = nstep
      m(3) = 0.0d0
      passed = .false.
      qh(:, 0) = st%q
      do k = 1, nstep
         call vstep_peri(st, 0.05d0, armed)
         call cau_peri_pes(st%q, mass3, gg, vv)
         m(4) = min(m(4), st%q(9))
         m(5) = max(m(5), st%q(9))
         qh(:, k) = st%q
         if (st%q(9) < 5.9d0) passed = .true.
         if (passed .and. st%q(9) > 6.0d0) then
            m(2) = k
            m(3) = 1.0d0                     ! scattered out (past Z_TERM)
            exit
         end if
      end do
      h1 = tot_e(st, vv)
      m(1) = (h1 - h0)
      m(6) = st%q(9)
      m(7:9) = 0.0d0
      if (armed) then
         m(7) = cau_peri_col_nup(st)
         m(8) = cau_peri_col_ndn(st)
         m(9) = cau_peri_col_eelec(st)
      end if
   end subroutine traj_peri

   subroutine traj_tab(st, nstep, armed, m, qh)
      type(state_t), intent(inout) :: st
      integer, intent(in) :: nstep
      logical, intent(in) :: armed
      real(8), intent(out) :: m(9)
      real(8), intent(out) :: qh(9, 0:nstep)
      real(8) :: gg(9), vv, h0, h1
      logical :: passed
      integer :: k
      call cau_tab_pes(st%q, mass3, gg, vv)
      h0 = tot_e(st, vv)
      m(4) = st%q(9)
      m(5) = st%q(9)
      m(2) = nstep
      m(3) = 0.0d0
      passed = .false.
      qh(:, 0) = st%q
      do k = 1, nstep
         call vstep_tab(st, 0.05d0, armed)
         call cau_tab_pes(st%q, mass3, gg, vv)
         m(4) = min(m(4), st%q(9))
         m(5) = max(m(5), st%q(9))
         qh(:, k) = st%q
         if (st%q(9) < 5.9d0) passed = .true.
         if (passed .and. st%q(9) > 6.0d0) then
            m(2) = k
            m(3) = 1.0d0                     ! scattered out (past Z_TERM)
            exit
         end if
      end do
      h1 = tot_e(st, vv)
      m(1) = (h1 - h0)
      m(6) = st%q(9)
      m(7:9) = 0.0d0
      if (armed) then
         m(7) = cau_tab_col_nup(st)
         m(8) = cau_tab_col_ndn(st)
         m(9) = cau_tab_col_eelec(st)
      end if
   end subroutine traj_tab

   subroutine vstep_peri(st, dt, evolve)
      type(state_t), intent(inout) :: st
      real(8), intent(in) :: dt
      logical, intent(in) :: evolve
      real(8) :: gg(9), vv
      integer :: ia
      call cau_peri_pes(st%q, mass3, gg, vv)
      do ia = 1, 3
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      do ia = 1, 3
         st%q(3*ia - 2:3*ia) = st%q(3*ia - 2:3*ia) + st%p(3*ia - 2:3*ia)/mass3(ia)*dt
      end do
      call cau_peri_pes(st%q, mass3, gg, vv)
      do ia = 1, 3
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      if (evolve) then
         st%t = st%t + dt
         call cau_peri_prop(st, dt)
      end if
   end subroutine vstep_peri

   subroutine vstep_tab(st, dt, evolve)
      type(state_t), intent(inout) :: st
      real(8), intent(in) :: dt
      logical, intent(in) :: evolve
      real(8) :: gg(9), vv
      integer :: ia
      call cau_tab_pes(st%q, mass3, gg, vv)
      do ia = 1, 3
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      do ia = 1, 3
         st%q(3*ia - 2:3*ia) = st%q(3*ia - 2:3*ia) + st%p(3*ia - 2:3*ia)/mass3(ia)*dt
      end do
      call cau_tab_pes(st%q, mass3, gg, vv)
      do ia = 1, 3
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      if (evolve) then
         st%t = st%t + dt
         call cau_tab_prop(st, dt)
      end if
   end subroutine vstep_tab

   real(8) function tot_e(st, v)
      type(state_t), intent(in) :: st
      real(8), intent(in) :: v
      integer :: ia
      tot_e = 0.0d0
      do ia = 1, 3
         tot_e = tot_e + dot_product(st%p(3*ia - 2:3*ia), st%p(3*ia - 2:3*ia))/(2.0d0*mass3(ia))
      end do
      tot_e = tot_e/e_conv + v*ev_kcal
   end function tot_e

   !------------------------------------------------------------------
   ! write_identity_table(fname) - the channel identity table sampled
   !                 from the cau_peri analytic chain (z 1.0..6.0 @0.05 - deep
   !                 enough for the armed mean-field push)
   !------------------------------------------------------------------
   subroutine write_identity_table(fname)
      character(len=*), intent(in) :: fname
      character(len=8) :: shn(5) = [character(len=3) :: 'E0', 'h', 'U', 'J', 'VCP']
      character(len=8) :: sin_(3) = ['top   ', 'bridge', 'hollow']
      integer, parameter :: nz = 117
      real(8), parameter :: z0 = 0.2d0, dz = 0.05d0
      integer :: u, ish, isi, k
      real(8) :: z, e, gx, gy, gz
      real(8) :: qs(9), hv(2), uv, jv, vcv
      open (newunit=u, file=fname, status='replace', action='write')
      write (u, '(a)') '# bench identity table (sampled from the cau_peri analytic chain)'
      write (u, '(a,f8.4)') 'LATTICE ', a_lat
      do ish = 1, 5
         write (u, '(a,a)') 'SHEET ', trim(shn(ish))
         do isi = 1, 3
            write (u, '(a,a,f9.4,f9.4)') 'SITE ', trim(sin_(isi)), site_xy(isi, 1), site_xy(isi, 2)
            do k = 1, nz
               z = z0 + (k - 1)*dz
               if (ish == 1) then
                  call calc_pot_spec_and_grad(parm_nas, site_xy(isi, 1), site_xy(isi, 2), z, e, gx, gy, gz)
               else
                  qs = 0.0d0
                  qs(6) = 0.9d0
                  qs(7:9) = [site_xy(isi, 1), site_xy(isi, 2), z]
                  call cau_peri_hubbard_val(qs, hv, uv, jv, vcv)
                  select case (ish)
                  case (2); e = 0.5d0*(hv(1) + hv(2))
                  case (3); e = uv
                  case (4); e = jv
                  case (5); e = vcv
                  end select
               end if
               write (u, '(f10.5,1x,es24.16)') z, e
            end do
         end do
      end do
      close (u)
   end subroutine write_identity_table

   subroutine copy_file(from, to)
      character(len=*), intent(in) :: from, to
      integer :: ua, ub, ios
      character(len=256) :: line
      open (newunit=ua, file=from, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print '(a)', 'fatal: unreadable '//trim(from)
         stop 1
      end if
      open (newunit=ub, file=to, status='replace', action='write')
      do
         read (ua, '(a)', iostat=ios) line
         if (ios /= 0) exit
         write (ub, '(a)') trim(line)
      end do
      close (ua)
      close (ub)
   end subroutine copy_file
end program bench_cau_tab
