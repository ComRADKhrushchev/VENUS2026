!=====================================================================
! check_cau_tab.f90 - the table-PES container criterion program
!                         (criteria defined FIRST in this header)
! Criteria:
!   (t1) Akima math leg: a table built from f(z) = sin(1.3z) + 0.2z
!       (0.05 A grid, all six sheets) - midpoint evaluations vs the
!       exact value (1e-4) and slope (1e-3); the above-range hold
!       branch returns the last point with an exactly zero slope
!   (t4) loader leg: the generated table parses back with the written
!       lattice (2.95), the expected point count at every sheet x site
!   (t2c) table ground kernel vs central difference of E0 at the
!       GENERIC site (C above the anchor offset by (0.37, 0.21, 2.2125)):
!       <= 1e-5*max(1,|g|) eV/A (the table kernel is the exact
!       derivative of its own interpolant - the (d)-class tolerance,
!       not the analytic-engine FD allowance of check_cau_peri (c))
!   (t2dh/t2h) Hubbard kernels and the composed fold vs central
!       differences of the value deliveries: <= 1e-5*max(1,|g|)
!   (t3a) identity-table site agreement ON grid points: cau_tab vs
!       cau_peri ground/hubbard values at the three site geometries
!       (z = 2.20 on the 0.05 grid): <= 1e-9 (both evaluate the same
!       analytic sheet through the same site value; only Fourier
!       recombination rounding remains)
!   (t3b) identity-table site agreement OFF grid (z = 2.225, the
!       midpoint): <= 2e-3*max(1,|v|) - the measured Akima error at
!       the 0.05 grid, REPORTED (off-site lateral difference between
!       the energy-Fourier of the table and the parameter-Fourier of
!       the analytic engine is a modeling choice, not covered here -
!       both probes sit exactly ON the site positions)
!   (a)(b)(e)(f)(f2)(g)(i1)(i2)(i3)(h2)(o)(j)(k)(l) - the engine legs
!       ported verbatim from check_cau_peri (same thresholds): the
!       renamed TDHF copy is exercised in table mode end to end
! Identity table: generated at run time at the three sites on the z
!   grid [0.5, 6.0] step 0.05 - the E0 sheet samples the cau_peri
!   analytic ground sheet (calc_POT_spec), the h/U/J/VCP sheets sample
!   cau_peri_hubbard_val (h tilde = the average of h(1)/h(2), the
!   crystal field cancels) - written to tab_gen.txt, then loaded
!   through cau_tab_load: the loader, the Fourier lateral combination,
!   the channel delivery and the deliveries are all covered by the
!   comparison against cau_peri (the channel-direct scheme 2026-10-01).
! Compile (from scripts/, Git Bash, gfortran; both containers compile
!   side by side - module names do not collide):
!   gfortran -O2 -o check_cau_tab.exe \
!     ../src/globals/consts.f90 ../src/globals/state.f90 \
!     ../src/utils/linalg.f90 ../src/utils/rng.f90 \
!     ../systems/inactive/cau_peri/surface_potentials.f90 \
!     ../systems/inactive/cau_peri/cau_peri_pot.f90 \
!     ../systems/inactive/cau_tab/cau_tab_tab.f90 \
!     ../systems/inactive/cau_tab/cau_tab_pot.f90 \
!     ../systems/inactive/cau_tab/cau_tab_tdhf.f90 \
!     ../systems/inactive/cau_tab/cau_tab_box.f90 \
!     check_cau_tab.f90
!   ./check_cau_tab.exe   (run from scripts/)
!=====================================================================
program check_cau_tab
   use cau_tab_tab,  only: tab_eval, tab_loaded, tab_lattice, tab_info_out, &
                           sh_E0, sh_h, sh_U, sh_J, sh_VCP, si_top, si_bri, si_hol
   use cau_tab_pot,  only: cau_tab_load, cau_tab_ground, cau_tab_hubbard, cau_tab_hubbard_val, &
                           ct_k => k_surf, ct_kg => k_ghost, ct_q0 => Q0_ghost
   use cau_tab_box,  only: cau_tab_pes_r
   use cau_tab_tdhf, only: cau_tab_tdhf_setup, cau_tab_tdhf_armed, cau_tab_pes, cau_tab_prop, &
                           bath_out, rho_seed_out, rho_herm_out, hcore_out, &
                           kernel_caches_out, hub_vals_out, cau_tab_sync, &
                           cau_tab_col_nup, cau_tab_col_ndn, cau_tab_col_mag, &
                           cau_tab_col_eelec, cau_tab_col_emf, rho_ref_out
   use cau_peri_pot, only: cau_peri_load, cau_peri_ground, cau_peri_hubbard_val
   use surface_potentials, only: cp_k => k_surf, cp_kg => k_ghost, cp_q0 => Q0_ghost, &
                                 parm_nas, parm_coup, parm_1d, parm_3p, parm_mc, parm_pc, &
                                 calc_pot_spec_and_grad, calc_pot0_and_grad, calc_coup_and_grad
   use state,        only: state_t, state_create, state_destroy
   implicit none
   character(len=*), parameter :: shelf = '../systems/inactive/'
   real(8), parameter :: ev_kcal = 23.0605d0, e_conv = 0.04184d0
   real(8), parameter :: a_lat = 2.95d0
   real(8), parameter :: site_xy(3, 2) = reshape([0.0d0, 1.475d0, 1.475d0, &
                                                  0.0d0, 0.0d0, 0.8515917d0], [3, 2])
   integer, parameter :: nz_grid = 111        ! z grid [0.5, 6.0] step 0.05
   real(8), parameter :: z_grid0 = 0.5d0, dz_grid = 0.05d0
   integer :: natoms, ncoord, npass, nfail, dirn, i
   integer :: kk, ii, i2, is2
   real(8) :: offd
   complex(8), allocatable :: rref(:, :, :)
   real(8), allocatable :: q0(:), mass(:), g(:), gnum(:)
   real(8) :: qr(9), gr(9), mr(3), vr          ! (l) roster-order scratch
   real(8) :: e0, v, vp, vm, dx, s13
   real(8) :: au_mass = 196.97d0, c_mass = 12.01d0
   type(state_t) :: sta
   npass = 0
   nfail = 0
   s13 = 1.0d0/sqrt(3.0d0)

   ! ---- 0. engine init (cau_peri analytic branch = generator + comparator;
   !      cau_tab table branch = the examinee) + the 3-atom geometry ----
   call copy_file(trim(shelf)//'cau_peri/Analytic_Potential.txt', 'Analytic_Potential.txt')
   call cau_peri_load()
   natoms = 3
   ncoord = 9
   allocate (q0(ncoord), mass(natoms), g(ncoord), gnum(ncoord))
   q0 = 0.0d0
   q0(4:6) = [0.0d0, 0.0d0, 0.9d0]         ! ghost Y at a generic offset
   q0(7:9) = [0.37d0, 0.21d0, 2.2125d0]     ! incident C at the generic site, z at a
   ! table-segment MIDPOINT (2.2125 between the 2.20/2.25 knots) - the central
   ! differences stay inside one polynomial piece; probing ON a knot would
   ! straddle the C1 slope kink and add ~h*dE''/4 ~ 1e-4 to the FD
   mass(1) = au_mass                        ! surface atom S
   mass(2) = au_mass                        ! ghost Y
   mass(3) = c_mass                         ! incident C

   ! ---- (t1) Akima math leg on a known function table ----
   block
      real(8) :: ev, ed, fref, fpref, zhold
      integer :: k
      call write_fn_table('tab_fn.txt')
      call cau_tab_load('tab_fn.txt')
      do k = 34, 96, 31                     ! midpoints 2.175 / 3.725 / 5.275
         vp = z_grid0 + (k - 1)*dz_grid - 0.5d0*dz_grid
         call tab_eval(sh_J, si_top, vp, ev, ed)
         fref = sin(1.3d0*vp) + 0.2d0*vp
         fpref = 1.3d0*cos(1.3d0*vp) + 0.2d0
         call cmp_one('t1 val', ev, fref, 1.0d-4, npass, nfail)
         call cmp_one('t1 slo', ed, fpref, 1.0d-3, npass, nfail)
      end do
      ! the above-range hold branch: the last point value, an exactly zero slope
      zhold = z_grid0 + (nz_grid - 1)*dz_grid
      call tab_eval(sh_VCP, si_hol, 7.3d0, ev, ed)
      call cmp_one('t1 hold', ev, sin(1.3d0*zhold) + 0.2d0*zhold, 1.0d-9, npass, nfail)
      call check_c('t1 hold0', ed == 0.0d0, npass, nfail)
   end block

   ! ---- generate + load the identity table (the physics examinee) ----
   call write_gen_table('tab_gen.txt')
   call cau_tab_load('tab_gen.txt')

   ! ---- (t4) loader leg ----
   block
      real(8) :: lat
      integer :: npts(5, 3)
      call tab_info_out(lat, npts)
      call check_c('t4 lat', abs(lat - a_lat) < 1.0d-12, npass, nfail)
      call check_c('t4 npts', all(npts == nz_grid), npass, nfail)
      call check_c('t4 load', tab_loaded(), npass, nfail)
   end block

   ! ---- oscillator constants: both containers on the box defaults ----
   ct_k = [1.0d0, 1.0d0, 1.0d0]
   ct_kg = [1.0d0, 1.0d0, 1.0d0]
   ct_q0 = [0.0d0, 0.0d0, 0.0d0]
   cp_k = [1.0d0, 1.0d0, 1.0d0]
   cp_kg = [1.0d0, 1.0d0, 1.0d0]
   cp_q0 = [0.0d0, 0.0d0, 0.0d0]

   ! ---- (t2c) table ground kernel vs central difference ----
   dx = 1.0d-4
   do dirn = 1, 4                          ! C x/y/z + surface atom S x (dof 1)
      call cau_tab_ground(q0, e0, g)       ! the kernel AT q0 (the check_cau_peri
      ! (c) pattern let fd_kernel's last (-dx) call overwrite g and absorbed the
      ! shift in its 5e-3 engine allowance; the 1e-5 table tolerance cannot)
      call fd_kernel(dof(dirn), dx, gnum)
      call cmp_one(ctag('t2c', dirn), g(dof(dirn)), gnum(dof(dirn)), 1.0d-5, npass, nfail)
   end do

   ! ---- (t3a/t3b) identity-table site agreement vs cau_peri ----
   block
      real(8) :: qs(9), eg1, eg2, gg(9), h1(2), h2(2), u1, u2, j1, j2, vc1, vc2
      integer :: isi
      character(len=6) :: tag
      do isi = 1, 3
         ! (t3a) ON the grid: z = 2.20 (both sides evaluate the same sheet value)
         qs = 0.0d0
         qs(4:6) = [0.0d0, 0.0d0, 0.9d0]
         qs(7:9) = [site_xy(isi, 1), site_xy(isi, 2), 2.20d0]
         call cau_tab_ground(qs, eg1, gg)
         call cau_peri_ground(qs, eg2, gg)
         call cmp_one(site_tag('t3a E0', isi), eg1, eg2, 1.0d-9, npass, nfail)
         call cau_tab_hubbard_val(qs, h1, u1, j1, vc1)
         call cau_peri_hubbard_val(qs, h2, u2, j2, vc2)
         call cmp_one(site_tag('t3a h1', isi), h1(1), h2(1), 1.0d-9, npass, nfail)
         call cmp_one(site_tag('t3a h2', isi), h1(2), h2(2), 1.0d-9, npass, nfail)
         call cmp_one(site_tag('t3a U', isi), u1, u2, 1.0d-9, npass, nfail)
         call cmp_one(site_tag('t3a J', isi), j1, j2, 1.0d-9, npass, nfail)
         call cmp_one(site_tag('t3a V', isi), vc1, vc2, 1.0d-9, npass, nfail)
         ! (t3b) OFF the grid midpoint: the Akima error at the 0.05 grid, REPORTED
         qs(9) = 2.225d0
         call cau_tab_ground(qs, eg1, gg)
         call cau_peri_ground(qs, eg2, gg)
         call cmp_one(site_tag('t3b E0', isi), eg1, eg2, 2.0d-3, npass, nfail)
         call cau_tab_hubbard_val(qs, h1, u1, j1, vc1)
         call cau_peri_hubbard_val(qs, h2, u2, j2, vc2)
         call cmp_one(site_tag('t3b h1', isi), h1(1), h2(1), 2.0d-3, npass, nfail)
         call cmp_one(site_tag('t3b h2', isi), h1(2), h2(2), 2.0d-3, npass, nfail)
         call cmp_one(site_tag('t3b U', isi), u1, u2, 2.0d-3, npass, nfail)
         call cmp_one(site_tag('t3b J', isi), j1, j2, 2.0d-3, npass, nfail)
         call cmp_one(site_tag('t3b V', isi), vc1, vc2, 2.0d-3, npass, nfail)
         print '(a,i0,a,es9.2,a,es9.2)', 't3b report site ', isi, ': max|dE0| = ', abs(eg1 - eg2), &
            '  |dJ| = ', abs(j1 - j2)
      end do
   end block

   ! ---- (l) roster-order boundary wrapper: permutation equivalence ----
   call cau_tab_pes(q0, mass, g, v)               ! engine order, direct (unarmed)
   qr = [q0(7:9), q0(1:6)]                        ! roster order: [C, S, Y]
   mr = [mass(3), mass(1), mass(2)]
   call cau_tab_pes_r(qr, mr, gr, vr)
   if (vr == v .and. all(gr == [g(7:9), g(1:6)])) then
      npass = npass + 1
      print '(a)', '(l) roster wrapper: PASS (permutation equivalence, bitwise)'
   else
      nfail = nfail + 1
      print '(a)', 'FAIL (l) roster wrapper: permutation mismatch'
   end if

   ! ---- (i1) degenerate point: exactness + tight drift ----
   if (cau_tab_tdhf_armed()) then
      print '(a)', 'FAIL i1 order: the package armed before the degenerate-point leg'
      nfail = nfail + 1
   end if
   call crit_traj_unarmed(npass, nfail)

   ! ---- (a) + production arm ----
   call cau_tab_tdhf_setup(4, 7.0d0, 300.0d0, 'GS', 1, 1, 0, 0.05d0)
   call crit_gauss(npass, nfail, s13)
   call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 1, 0, 0.05d0)
   call crit_bath(npass, nfail)
   call crit_seed_gs(npass, nfail)

   ! ---- (e) Hcore structure after one prop call ----
   call state_create(sta, natoms, 1)
   sta%q = q0
   sta%t = 0.0d0
   call cau_tab_prop(sta, 0.05d0)
   call crit_hcore(npass, nfail)
   call state_destroy(sta)

   ! ---- (g) unitary evolution ----
   call crit_unitary(npass, nfail)

   ! ---- (f2) 1D seed ----
   call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, '1D', 4, 1, 0, 0.05d0)
   call crit_seed_1d(npass, nfail)
   call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 1, 0, 0.05d0)

   ! ---- (t2dh) Hubbard kernels vs FD of the values ----
   call fd_hub_all(npass, nfail, dx)

   ! ---- (t2h) the composed fold ----
   call cau_tab_pes(q0, mass, g, v)
   do dirn = 1, 4
      call dof_fd_v(dof(dirn), dx, vp, vm)
      call cmp_one(ctag('t2h', dirn), g(dof(dirn)), (vp - vm)/(2*dx), 1.0d-5, npass, nfail)
   end do

   ! ---- (h2) the composed fold with a NONDIAGONAL rho (ported verbatim:
   ! the coupling-kernel partner elements read exactly the REAL parts) ----
   call state_create(sta, natoms, 1)
   sta%q = q0
   do kk = 1, 100
      sta%t = kk*0.05d0
      call cau_tab_prop(sta, 0.05d0)
   end do
   call rho_ref_out(rref)
   offd = 0.0d0
   do is2 = 1, 2
      do ii = 1, 2
         do i2 = 3, size(rref, 1)
            offd = max(offd, abs(real(rref(ii, i2, is2))))
         end do
      end do
   end do
   if (offd < 1.0d-3) then
      print '(a,es9.2)', 'FAIL (h2) backdrop: max |Re rho(imp,bath)| = ', offd, &
         ' - the leg is blind (evolve longer or a closer geometry)'
      nfail = nfail + 1
   else
      call cau_tab_pes(sta%q, mass, g, v)
      do dirn = 1, 4
         call dof_fd_q(sta%q, dof(dirn), dx, vp, vm)
         call cmp_one(ctag('h2', dirn), g(dof(dirn)), (vp - vm)/(2*dx), 1.0d-5, npass, nfail)
      end do
      print '(a,es9.2)', '(h2) nondiagonal-rho fold backdrop: max |Re rho(imp,bath)| = ', offd
   end if
   call state_destroy(sta)

   ! ---- (i3) static-geometry mean-field theorem (msbstp=16) ----
   block
      real(8) :: eq(9), g9(9), vv, eb, ee2, mf1, mf2
      integer :: kk2
      call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 16, 1, 0, 0.05d0)
      eq = q0
      eq(9) = 2.2d0                          ! near-collision height
      sta%q = eq
      sta%t = 0.0d0
      call cau_tab_pes(eq, mass, g9, vv)     ! establish hcore for the baseline
      eb = 0.0d0
      do kk2 = 1, 300
         sta%t = kk2*0.05d0
         call cau_tab_prop(sta, 0.05d0)
      end do
      ee2 = cau_tab_col_eelec(sta)
      mf2 = cau_tab_col_emf(sta)
      ! re-evaluate the baseline: seed + the same hcore
      call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 16, 1, 0, 0.05d0)
      call cau_tab_pes(eq, mass, g9, vv)
      mf1 = cau_tab_col_emf(sta)
      eb = cau_tab_col_eelec(sta) + mf1
      if (abs(ee2 + mf2 - eb) > 0.05d0) then
         print '(a,es10.3)', 'FAIL (i3) static theorem: |dE| = ', abs(ee2 + mf2 - eb), &
            ' > 0.05 eV (300 steps, msbstp=16)'
         nfail = nfail + 1
      else
         print '(a,es10.3)', '(i3) static mean-field theorem: PASS (|dE| = ', &
            abs(ee2 + mf2 - eb), ' eV over 300 steps, msbstp=16)'
      end if
      call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 1, 0, 0.05d0)
   end block

   ! ---- (o) far-field freeze (ruling 2026-09-29), table mode ----
   block
      real(8) :: qf(9), gf(9), ef, hv(2), uv, jv, vcpv, hb(2), ub, jb, vcb
      real(8) :: dhf(2, 9), duf(9), djf(9), dvcf(9), e1, e2v
      qf = 0.0d0
      qf(6) = 0.9d0
      qf(1) = 0.15d0
      qf(9) = 7.0d0
      call cau_tab_hubbard(qf, hv, uv, jv, vcpv, dhf, duf, djf, dvcf)
      call cau_tab_ground(qf, ef, gf)
      ! (o1) the exact-zero set at 7 A: coupling kernel, every z kernel row,
      ! the C-z ground kernel
      if (maxval(abs(dvcf)) < 1.0d-12 .and. abs(dhf(1, 3)) + abs(dhf(1, 6)) + abs(dhf(1, 9)) < 1.0d-12 &
          .and. abs(dhf(2, 3)) + abs(dhf(2, 6)) + abs(dhf(2, 9)) < 1.0d-12 &
          .and. abs(duf(3)) + abs(duf(6)) + abs(duf(9)) < 1.0d-12 &
          .and. abs(djf(3)) + abs(djf(6)) + abs(djf(9)) < 1.0d-12 &
          .and. abs(gf(9)) < 1.0d-12) then
         npass = npass + 1
         print '(a)', '(o1) far-field exact zeros: PASS (dvcp=0, z-kernels=0 at 7 A)'
      else
         nfail = nfail + 1
         print '(a,es10.2,a,es10.2,a,es10.2)', 'FAIL (o1) far-field exact zeros: |dvcp|=', &
            maxval(abs(dvcf)), ' |dh_z|=', maxval(abs(dhf(:, [3, 6, 9]))), &
            ' |g0z|=', abs(gf(9))
      end if
      ! (o2) lateral residuals on the frozen sheet: the site-asymptote
      ! differences of the sampled data (bounded; REPORTED)
      if (maxval(abs(dhf)) < 2.0d0 .and. maxval(abs(duf)) < 2.0d0 .and. &
          maxval(abs(djf)) < 2.0d0 .and. maxval(abs(gf(7:9))) < 2.0d0) then
         npass = npass + 1
         print '(a)', '(o2) far-field lateral bound: PASS (site-asymptote artifact, bounded)'
         print '(4(a,es10.2))', ' |dh|=', maxval(abs(dhf)), ' |du|=', maxval(abs(duf)), &
            ' |dj|=', maxval(abs(djf)), ' |g0_C|=', maxval(abs(gf(7:9)))
      else
         nfail = nfail + 1
         print '(a,4(a,es10.2))', 'FAIL (o2) far-field lateral residual: ', &
            ' |dh|=', maxval(abs(dhf)), ' |du|=', maxval(abs(duf)), &
            ' |dj|=', maxval(abs(djf)), ' |g0_C|=', maxval(abs(gf(7:9)))
      end if
      qf(9) = 9.0d0
      call cau_tab_hubbard(qf, hb, ub, jb, vcb, dhf, duf, djf, dvcf)
      if (hb(1) == hv(1) .and. hb(2) == hv(2) .and. ub == uv .and. jb == jv .and. &
          vcb == vcpv) then
         npass = npass + 1
         print '(a)', '(o) far-field values: PASS (frozen - 7 A equals 9 A)'
      else
         print '(a)', 'FAIL (o) far-field values drift beyond the freeze line'
         nfail = nfail + 1
      end if
      qf(9) = 5.3d0
      call cau_tab_ground(qf, ef, gf)
      qf(9) = 5.3d0 + dx
      call cau_tab_ground(qf, e1, g)
      qf(9) = 5.3d0 - dx
      call cau_tab_ground(qf, e2v, g)
      call cmp_one('o ramp', gf(9), (e1 - e2v)/(2*dx), 1.0d-5, npass, nfail)
   end block

   ! ---- (i2) armed trajectory ----
   call crit_traj_armed(npass, nfail)

   ! ---- (j) sync publication + (k) recorder columns ----
   call crit_sync_cols(npass, nfail)

   print '(/,a,i0,a,i0,a)', 'check_cau_tab: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) call exit(1)

contains

   ! copy_file(from, to) - portable text copy (READAP reads
   !                 Analytic_Potential.txt from the run directory)
   subroutine copy_file(from, to)
      character(len=*), intent(in) :: from, to
      integer :: ua, ub, ios
      character(len=256) :: line
      open (newunit=ua, file=from, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print '(a)', 'FAIL copy_file: unreadable '//trim(from)
         call exit(1)
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

   !------------------------------------------------------------------
   ! gen_sheet_eval(ish, x, y, z, e) - one identity-channel sample: E0
   !                 from the cau_peri analytic ground sheet, h/U/J/VCP
   !                 from cau_peri_hubbard_val at the same geometry
   !                 (h tilde = the average of h(1)/h(2); gradients
   !                 discarded)
   !------------------------------------------------------------------
   subroutine gen_sheet_eval(ish, x, y, z, e)
      integer, intent(in) :: ish
      real(8), intent(in) :: x, y, z
      real(8), intent(out) :: e
      real(8) :: gx, gy, gz
      real(8) :: qs(9), hv(2), uv, jv, vcv
      select case (ish)
      case (sh_E0)
         call calc_pot_spec_and_grad(parm_nas, x, y, z, e, gx, gy, gz)
      case (sh_h:sh_VCP)
         qs = 0.0d0
         qs(6) = 0.9d0
         qs(7:9) = [x, y, z]
         call cau_peri_hubbard_val(qs, hv, uv, jv, vcv)
         select case (ish)
         case (sh_h);   e = 0.5d0*(hv(1) + hv(2))
         case (sh_U);   e = uv
         case (sh_J);   e = jv
         case (sh_VCP); e = vcv
         end select
      case default
         print '(a)', 'FAIL gen_sheet_eval: bad sheet index'
         call exit(1)
      end select
   end subroutine gen_sheet_eval

   !------------------------------------------------------------------
   ! write_gen_table(fname) - the identity table: sample the analytic
   !                 engine at the three sites on the z grid (the
   !                 generator runs AFTER cau_peri_load filled READAP)
   !------------------------------------------------------------------
   subroutine write_gen_table(fname)
      character(len=*), intent(in) :: fname
      character(len=8) :: shn(5) = [character(len=3) :: 'E0', 'h', 'U', 'J', 'VCP']
      character(len=8) :: sin_(3) = ['top   ', 'bridge', 'hollow']
      integer :: u, ish, isi, k
      real(8) :: z, e
      open (newunit=u, file=fname, status='replace', action='write')
      write (u, '(a)') '# cau_tab identity table (generated by check_cau_tab)'
      write (u, '(a,f8.4)') 'LATTICE ', a_lat
      do ish = 1, 5
         write (u, '(a,a)') 'SHEET ', trim(shn(ish))
         do isi = 1, 3
            write (u, '(a,a,f9.4,f9.4)') 'SITE ', trim(sin_(isi)), site_xy(isi, 1), site_xy(isi, 2)
            do k = 1, nz_grid
               z = z_grid0 + (k - 1)*dz_grid
               call gen_sheet_eval(ish, site_xy(isi, 1), site_xy(isi, 2), z, e)
               write (u, '(f10.5,1x,es24.16)') z, e
            end do
         end do
      end do
      close (u)
   end subroutine write_gen_table

   ! write_fn_table(fname) - the (t1) function table (all six sheets
   !                 carry f(z) = sin(1.3z) + 0.2z)
   subroutine write_fn_table(fname)
      character(len=*), intent(in) :: fname
      character(len=8) :: shn(5) = [character(len=3) :: 'E0', 'h', 'U', 'J', 'VCP']
      character(len=8) :: sin_(3) = ['top   ', 'bridge', 'hollow']
      integer :: u, ish, isi, k
      real(8) :: z
      open (newunit=u, file=fname, status='replace', action='write')
      write (u, '(a)') '# cau_tab function table (t1 leg)'
      write (u, '(a,f8.4)') 'LATTICE ', a_lat
      do ish = 1, 5
         write (u, '(a,a)') 'SHEET ', trim(shn(ish))
         do isi = 1, 3
            write (u, '(a,a,f9.4,f9.4)') 'SITE ', trim(sin_(isi)), site_xy(isi, 1), site_xy(isi, 2)
            do k = 1, nz_grid
               z = z_grid0 + (k - 1)*dz_grid
               write (u, '(f10.5,1x,es24.16)') z, sin(1.3d0*z) + 0.2d0*z
            end do
         end do
      end do
      close (u)
   end subroutine write_fn_table

   integer function dof(d)
      integer, intent(in) :: d
      if (d == 4) then
         dof = 1                            ! probe the surface atom S x (3-atom model)
      else
         dof = d
      end if
   end function dof

   function ctag(tag, dirn) result(s)
      character(len=*), intent(in) :: tag
      integer, intent(in) :: dirn
      character(len=8) :: s
      write (s, '(a,i0)') tag, dirn
   end function ctag

   function site_tag(base, isi) result(s)
      character(len=*), intent(in) :: base
      integer, intent(in) :: isi
      character(len=10) :: s
      write (s, '(a,a1,i0)') base, '-', isi
   end function site_tag

   subroutine crit_gauss(np, nf, s)
      integer, intent(inout) :: np, nf
      real(8), intent(in) :: s
      real(8) :: em(4), ez
      call bath_out(em, ez)
      call check_c('a ext+', abs(em(4) - 3.5d0*(1.0d0 + s)/2.0d0) < 1.0d-12, np, nf)
      call check_c('a ext-', abs(em(1) + 3.5d0*(1.0d0 + s)/2.0d0) < 1.0d-12, np, nf)
      call check_c('a inn+', abs(em(3) - 3.5d0*(1.0d0 - s)/2.0d0) < 1.0d-12, np, nf)
      call check_c('a inn-', abs(em(2) + 3.5d0*(1.0d0 - s)/2.0d0) < 1.0d-12, np, nf)
   end subroutine crit_gauss

   subroutine crit_bath(np, nf)
      integer, intent(inout) :: np, nf
      real(8) :: em(20), ez, s_lo
      integer :: kk
      call bath_out(em, ez)
      do kk = 1, 20
         call check_c('b sym', abs(em(kk) + em(21 - kk)) < 1.0d-14, np, nf)
         call check_c('b rng', abs(em(kk)) <= 3.5d0 + 1.0d-12, np, nf)
      end do
      s_lo = 0.0d0
      do kk = 1, 20
         if (em(kk) < 0.0d0) s_lo = s_lo + em(kk)
      end do
      call check_c('b zero', abs(ez - 2.0d0*s_lo) < 1.0d-12, np, nf)
   end subroutine crit_bath

   subroutine crit_seed_gs(np, nf)
      integer, intent(inout) :: np, nf
      real(8) :: occ(2), trb
      call rho_seed_out(occ, trb)
      call check_c('f GS up', abs(occ(1) - 2.0d0) < 1.0d-14, np, nf)
      call check_c('f GS dn', abs(occ(2)) < 1.0d-14, np, nf)
      call check_c('f bath', abs(trb - 2.0d0 - 20.0d0) < 5.0d-1, np, nf)
      call check_c('f edge', exp(-3.454d0/(1.9872198404d-3*300.0d0/23.0605d0)) < 1.0d-3, np, nf)
   end subroutine crit_seed_gs

   subroutine crit_seed_1d(np, nf)
      integer, intent(inout) :: np, nf
      real(8) :: occ(2), trb
      call rho_seed_out(occ, trb)
      call check_c('f2 up', abs(occ(1) - 1.0d0) < 1.0d-14, np, nf)
      call check_c('f2 dn', abs(occ(2) - 1.0d0) < 1.0d-14, np, nf)
   end subroutine crit_seed_1d

   subroutine crit_hcore(np, nf)
      integer, intent(inout) :: np, nf
      real(8) :: hc(22, 22, 2)
      real(8) :: symm, mir
      integer :: ii, jj, is
      call hcore_out(hc)
      symm = 0.0d0
      mir = 0.0d0
      do is = 1, 2
         do ii = 1, 22
            do jj = 1, 22
               symm = max(symm, abs(hc(ii, jj, is) - hc(jj, ii, is)))
            end do
         end do
         do ii = 1, 2
            do jj = 3, 12
               mir = max(mir, abs(hc(ii, jj, is) - hc(ii, jj + 10, is)))
            end do
         end do
         mir = max(mir, abs(hc(1, 2, is)))
      end do
      call check_c('e sym', symm < 1.0d-14, np, nf)
      call check_c('e mir', mir < 1.0d-14, np, nf)
   end subroutine crit_hcore

   subroutine crit_unitary(np, nf)
      integer, intent(inout) :: np, nf
      real(8) :: tr0, tr1, herm
      integer :: kk
      call state_create(sta, natoms, 1)
      sta%q = q0
      sta%t = 0.0d0
      call rho_herm_out(tr0, herm)
      do kk = 1, 100
         sta%t = kk*0.05d0
         call cau_tab_prop(sta, 0.05d0)
      end do
      call rho_herm_out(tr1, herm)
      call check_c('g trace', abs(tr1 - tr0) < 1.0d-9, np, nf)
      call check_c('g herm', herm < 1.0d-12, np, nf)
      call state_destroy(sta)
   end subroutine crit_unitary

   subroutine fd_hub_all(np, nf, h)
      integer, intent(inout) :: np, nf
      real(8), intent(in) :: h
      real(8) :: gg(ncoord), vv, mass1(natoms)
      real(8) :: dhv(2, ncoord), duv(ncoord), djv(ncoord), dvvv(ncoord)
      real(8) :: hv(2), uv, jv, vcpv
      real(8) :: hp(2), up, jp, vpp, hm(2), um, jm, vmm
      real(8), allocatable :: qp(:)
      integer :: d
      mass1 = 1.0d0
      allocate (qp(ncoord))
      call cau_tab_pes(q0, mass1, gg, vv)
      call kernel_caches_out(dhv, duv, djv, dvvv)
      do d = 1, 4
         qp = q0
         qp(dof(d)) = q0(dof(d)) + h
         call cau_tab_pes(qp, mass1, gg, vv)
         call hub_vals_out(hp, up, jp, vpp)
         qp(dof(d)) = q0(dof(d)) - h
         call cau_tab_pes(qp, mass1, gg, vv)
         call hub_vals_out(hm, um, jm, vmm)
         qp(dof(d)) = q0(dof(d))
         call cmp_one(ctag('t2dh', d), dhv(1, dof(d)), (hp(1) - hm(1))/(2*h), 1.0d-5, np, nf)
         call cmp_one(ctag('t2dU', d), duv(dof(d)), (up - um)/(2*h), 1.0d-5, np, nf)
         call cmp_one(ctag('t2dJ', d), djv(dof(d)), (jp - jm)/(2*h), 1.0d-5, np, nf)
         call cmp_one(ctag('t2dV', d), dvvv(dof(d)), (vpp - vmm)/(2*h), 1.0d-5, np, nf)
      end do
      deallocate (qp)
   end subroutine fd_hub_all

   subroutine dof_fd_v(k, h, vpp, vmm)
      integer, intent(in) :: k
      real(8), intent(in) :: h
      real(8), intent(out) :: vpp, vmm
      real(8) :: gg(ncoord)
      real(8), allocatable :: qp(:)
      allocate (qp(ncoord))
      qp = q0
      qp(k) = q0(k) + h
      call cau_tab_pes(qp, mass, gg, vpp)
      qp(k) = q0(k) - h
      call cau_tab_pes(qp, mass, gg, vmm)
      deallocate (qp)
   end subroutine dof_fd_v

   ! dof_fd_q(qb, k, h, vpp, vmm) - the FD pair of the armed potential at an
   !                 arbitrary base geometry qb (rho frozen in module state)
   subroutine dof_fd_q(qb, k, h, vpp, vmm)
      real(8), intent(in) :: qb(:), h
      integer, intent(in) :: k
      real(8), intent(out) :: vpp, vmm
      real(8) :: gg(ncoord)
      real(8), allocatable :: qp(:)
      allocate (qp(ncoord))
      qp = qb
      qp(k) = qb(k) + h
      call cau_tab_pes(qp, mass, gg, vpp)
      qp(k) = qb(k) - h
      call cau_tab_pes(qp, mass, gg, vmm)
      deallocate (qp)
   end subroutine dof_fd_q

   subroutine crit_traj_unarmed(np, nf)
      integer, intent(inout) :: np, nf
      real(8) :: h0, h1, vv, e0_ref, gg(ncoord), g0_ref(ncoord)
      integer :: step
      call state_create(sta, natoms, 1)
      sta%q = q0
      sta%t = 0.0d0
      sta%p = 0.0d0
      sta%p(3) = -c_mass*0.30d0
      call cau_tab_pes(sta%q, mass, gg, vv)
      call cau_tab_ground(sta%q, e0_ref, g0_ref)
      call check_c('i1 v==E0', abs(vv - e0_ref) < 1.0d-14, np, nf)
      call check_c('i1 g==g0', maxval(abs(gg - g0_ref)) < 1.0d-14, np, nf)
      h0 = tot_e(sta, vv)
      do step = 1, 300
         call vstep(sta, 0.05d0, .false.)
         call cau_tab_pes(sta%q, mass, gg, vv)
      end do
      h1 = tot_e(sta, vv)
      call check_c('i1 drift', abs(h1 - h0) < 0.05d0, np, nf)
      print '(a,es12.4)', 'i1 report: unarmed table drift [kcal/mol] = ', h1 - h0
      call state_destroy(sta)
   end subroutine crit_traj_unarmed

   subroutine crit_traj_armed(np, nf)
      integer, intent(inout) :: np, nf
      real(8) :: h0, h1, vv, gg(ncoord), tr0, tr1, herm
      integer :: step
      call state_create(sta, natoms, 1)
      sta%q = q0
      sta%t = 0.0d0
      sta%p = 0.0d0
      sta%p(3) = -c_mass*0.30d0
      call cau_tab_pes(sta%q, mass, gg, vv)
      h0 = tot_e(sta, vv)
      call rho_herm_out(tr0, herm)
      do step = 1, 300
         call vstep(sta, 0.05d0, .true.)
         call cau_tab_pes(sta%q, mass, gg, vv)
      end do
      h1 = tot_e(sta, vv)
      call rho_herm_out(tr1, herm)
      call check_c('i2 bound', abs(h1 - h0) < 200.0d0, np, nf)
      call check_c('i2 trace', abs(tr1 - tr0) < 1.0d-6, np, nf)
      print '(a,es12.4,a,es10.2)', 'i2 report: armed table drift [kcal/mol] = ', h1 - h0, &
         '  dtrace = ', tr1 - tr0
      call state_destroy(sta)
   end subroutine crit_traj_armed

   subroutine fd_kernel(k, dx, gn)
      integer, intent(in) :: k
      real(8), intent(in) :: dx
      real(8), intent(out) :: gn(:)
      real(8) :: ep, em
      real(8) :: gl(ncoord)
      real(8), allocatable :: qp(:)
      allocate (qp(ncoord))
      qp = q0
      qp(k) = q0(k) + dx
      call cau_tab_ground(qp, ep, gl)
      qp(k) = q0(k) - dx
      call cau_tab_ground(qp, em, gl)
      gn = 0.0d0
      gn(k) = (ep - em)/(2*dx)
      deallocate (qp)
   end subroutine fd_kernel

   subroutine vstep(st, dt, evolve)
      type(state_t), intent(inout) :: st
      real(8), intent(in) :: dt
      logical, intent(in) :: evolve
      real(8) :: gg2(ncoord), vv2
      integer :: ia
      call cau_tab_pes(st%q, mass, gg2, vv2)
      do ia = 1, natoms
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg2(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      do ia = 1, natoms
         st%q(3*ia - 2:3*ia) = st%q(3*ia - 2:3*ia) + st%p(3*ia - 2:3*ia)/mass(ia)*dt
      end do
      call cau_tab_pes(st%q, mass, gg2, vv2)
      do ia = 1, natoms
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg2(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      if (evolve) then
         st%t = st%t + dt
         call cau_tab_prop(st, dt)
      end if
   end subroutine vstep

   real(8) function tot_e(st, v)
      type(state_t), intent(in) :: st
      real(8), intent(in) :: v
      integer :: ia
      tot_e = 0.0d0
      do ia = 1, natoms
         tot_e = tot_e + dot_product(st%p(3*ia - 2:3*ia), st%p(3*ia - 2:3*ia))/(2.0d0*mass(ia))
      end do
      tot_e = tot_e/e_conv + v*ev_kcal
   end function tot_e

   !------------------------------------------------------------------
   ! (j)/(k) sync publication into sta%s%rho + the recorder columns
   !------------------------------------------------------------------
   subroutine crit_sync_cols(np, nf)
      integer, intent(inout) :: np, nf
      type(state_t) :: st
      complex(8), allocatable :: ref(:, :, :)
      real(8) :: vv, dmax
      integer :: i, j, is, n
      call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 1, 0, 0.05d0)   ! fresh GS seed
      call state_create(st, natoms, 1)
      st%q = q0
      st%t = 0.0d0
      st%p = 0.0d0
      st%p(3) = -c_mass*0.30d0
      call cau_tab_pes(st%q, mass, g, vv)                  ! armed fold at the seed
      call check_c('k nup seed', abs(cau_tab_col_nup(st) - 2.0d0) < 1.0d-14, np, nf)
      call check_c('k ndn seed', abs(cau_tab_col_ndn(st)) < 1.0d-14, np, nf)
      call check_c('k mag seed', abs(cau_tab_col_mag(st) - 2.0d0) < 1.0d-14, np, nf)
      call cau_tab_sync(st%s)
      call rho_ref_out(ref)
      n = size(ref, 1)
      call check_c('j shape', allocated(st%s%rho) .and. all(shape(st%s%rho) == [2*n, 2*n]), np, nf)
      call check_c('j nsurf', st%s%n_surf == 2*n, np, nf)
      dmax = 0.0d0
      do is = 1, 2
         do j = 1, n
            do i = 1, n
               dmax = max(dmax, abs(st%s%rho((is - 1)*n + i, (is - 1)*n + j) - ref(i, j, is)))
            end do
         end do
      end do
      call check_c('j copy', dmax == 0.0d0, np, nf)
      do i = 1, 50
         st%t = st%t + 0.05d0
         call cau_tab_prop(st, 0.05d0)
      end do
      call cau_tab_sync(st%s)
      call check_c('k nup range', cau_tab_col_nup(st) > -0.05d0 .and. cau_tab_col_nup(st) < 3.05d0, np, nf)
      call check_c('k ndn range', cau_tab_col_ndn(st) > -0.05d0 .and. cau_tab_col_ndn(st) < 3.05d0, np, nf)
      call check_c('k eelec fin', abs(cau_tab_col_eelec(st)) < 1.0d4, np, nf)
      call check_c('k emf fin', abs(cau_tab_col_emf(st)) < 1.0d4, np, nf)
      deallocate (ref)
      call state_destroy(st)
   end subroutine crit_sync_cols

   subroutine cmp_one(tag, ga, gb, tol, np, nf)
      character(len=*), intent(in) :: tag
      real(8), intent(in) :: ga, gb, tol
      integer, intent(inout) :: np, nf
      real(8) :: tv
      tv = tol*max(1.0d0, abs(ga), abs(gb))
      if (abs(ga - gb) <= tv) then
         np = np + 1
      else
         nf = nf + 1
         print '(a,a,a,es12.4,a,es12.4,a,es12.4)', 'FAIL ', tag, ' an=', ga, ' fd=', gb, ' tol=', tv
      end if
   end subroutine cmp_one

   subroutine check_c(tag, ok, np, nf)
      character(len=*), intent(in) :: tag
      logical, intent(in) :: ok
      integer, intent(inout) :: np, nf
      if (ok) then
         np = np + 1
      else
         nf = nf + 1
         print '(a,a)', 'FAIL ', tag
      end if
   end subroutine check_c
end program check_cau_tab
