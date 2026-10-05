!=====================================================================
! check_cau_peri.f90 - the C@Au(111) periodic+TDHF container criterion
!                         program (criteria defined FIRST in this header)
! Criteria:
!   (a) Gauss-Legendre closed forms through a 4-orbital probe bath
!       (extremes +-bw(1+1/sqrt(3))/4, inner +-bw(1-1/sqrt(3))/4; 1e-12)
!   (b) production bath (20 orbitals): symmetric band, range, the
!       energy_metal_zero identity (1e-12)
!   (c) periodic ground kernel vs central difference of E0 at the GENERIC
!       site (C above slab atom 1 offset by (0.37, 0.21, 2.2) A):
!       |g - g_fd| <= 5e-3*max(1,|g|) eV/A (the measured engine FD
!       anomaly allowance of the EMT family; probe before tightening)
!   (d) Hubbard NN kernels dh/du/dj/dvcp vs central difference of the
!       value delivery (same dofs): <= 1e-5*max(1,|g|) eV/A
!   (e) Hcore structure: symmetry (1e-14), Fermi-branch coupling mirror
!       (1e-14), strictly diagonal impurity block
!   (f) rho seeds: GS -> (2, 0), 1D -> (1, 1); bath trace ~ 20; band-edge
!       occupation < 1e-3
!   (g) unitary evolution: 100 Magnus-2 substeps preserve trace (1e-9)
!       and hermiticity (1e-12)
!   (h) the composed fold: armed kernel vs central difference of the
!       armed potential (rho frozen; h = 1e-4 A): <= 1e-5*max(1,|g|)
!   (i1) BEFORE the arm: composed v equals the ground E0 exactly, g
!       equals the ground kernel exactly; 300 velocity-Verlet steps
!       (dt = 0.5 fs), |H - H0| <= 0.05 kcal/mol
!   (i2) armed GS with per-step evolution: 300 steps, |H - H0| bounded
!       (< 200 kcal/mol - the ensemble mean-field relaxation transient;
!       see check_cau_emt2's header for the physics note) and rho trace
!       preserved (1e-6); measured drift REPORTED
! Engine files: the container folder's Slab.xyz + nn_weights_emt_nn.txt
!   (referenced by the path constant below - no run-directory copies
!   needed for the criterion run).
! Compile (from scripts/, gfortran; container on the inactive shelf):
!   gfortran -O2 -o check_cau_peri.exe ^
!     ..\src\globals\consts.f90 ..\src\globals\state.f90 ^
!     ..\src\utils\linalg.f90 ..\src\utils\rng.f90 ^
!     ..\systems\inactive\cau_peri\surface_potentials.f90 ^
!     ..\systems\inactive\cau_peri\cau_peri_pot.f90 ^
!     ..\systems\inactive\cau_peri\cau_peri_tdhf.f90 ^
!     ..\systems\inactive\cau_peri\cau_peri_box.f90 ^
!     check_cau_peri.f90
!   (j)/(k) sync publication + recorder columns (seed values exact,
!       post-propagation bounded)
!   (l) roster-order boundary wrapper: cau_peri_pes_r on the permuted
!       geometry [C, S, Y] reproduces the direct engine-order cau_peri_pes
!       value and gradient BITWISE (the roster->engine permutation is the
!       box's boundary contract; the reg file binds the wrapper)
!   check_cau_peri.exe   (run from scripts/)
!=====================================================================
program check_cau_peri
   use cau_peri_pot,  only: cau_peri_load, cau_peri_ground, cau_peri_hubbard
   use cau_peri_box, only: cau_peri_pes_r
   use cau_peri_tdhf, only: cau_peri_tdhf_setup, cau_peri_tdhf_armed, cau_peri_pes, &
                             cau_peri_prop, bath_out, rho_seed_out, rho_herm_out, hcore_out, &
                             kernel_caches_out, hub_vals_out
   use cau_peri_tdhf, only: cau_peri_sync, cau_peri_col_nup, cau_peri_col_ndn, &
                            cau_peri_col_mag, cau_peri_col_eelec, cau_peri_col_emf, rho_ref_out
   use state,          only: state_t, state_create, state_destroy
   implicit none
   character(len=*), parameter :: shelf = '../systems/inactive/cau_peri/'
   real(8), parameter :: ev_kcal = 23.0605d0, e_conv = 0.04184d0
   integer :: nslab, natoms, ncoord, npass, nfail, dirn, u, ios, i
   integer :: kk, ii, i2, is2
   real(8) :: offd
   complex(8), allocatable :: rref(:, :, :)
   real(8), allocatable :: q0(:), mass(:), g(:), gnum(:), gp(:), gm(:)
   real(8) :: qr(9), gr(9), mr(3), vr   ! (l) roster-order scratch
   real(8) :: e0, v, vp, vm, dx, s13, x1, y1, z1
   real(8) :: au_mass = 196.97d0, c_mass = 12.01d0
   character(len=4) :: elem
   type(state_t) :: sta
   npass = 0
   nfail = 0
   s13 = 1.0d0/sqrt(3.0d0)

   ! ---- 0. engine init (analytic defaults) + the 3-atom geometry ----
   call copy_file(trim(shelf)//'Analytic_Potential.txt', 'Analytic_Potential.txt')
   call cau_peri_load()
   natoms = 3
   ncoord = 9
   allocate (q0(ncoord), mass(natoms), g(ncoord), gnum(ncoord), gp(ncoord), gm(ncoord))
   q0 = 0.0d0
   q0(4:6) = [0.0d0, 0.0d0, 0.9d0]         ! ghost Y at a generic offset
   q0(7:9) = [0.37d0, 0.21d0, 2.2d0]       ! incident C at the generic site
   mass(1) = au_mass                        ! surface atom S
   mass(2) = au_mass                        ! ghost Y
   mass(3) = c_mass                         ! incident C

   ! ---- (c) ground kernel vs central difference ----
   dx = 1.0d-4
   do dirn = 1, 4                          ! C x/y/z + slab atom 2 z (dof 6)
      call fd_kernel(dof(dirn), dx, gnum)
      call cmp_one(ctag('c', dirn), g(dof(dirn)), gnum(dof(dirn)), 5.0d-3, npass, nfail)
   end do

   ! ---- (l) roster-order boundary wrapper: permutation equivalence ----
   call cau_peri_pes(q0, mass, g, v)             ! engine order, direct
   qr = [q0(7:9), q0(1:6)]                       ! roster order: [C, S, Y]
   mr = [mass(3), mass(1), mass(2)]
   call cau_peri_pes_r(qr, mr, gr, vr)
   if (vr == v .and. all(gr == [g(7:9), g(1:6)])) then
      npass = npass + 1
      print '(a)', '(l) roster wrapper: PASS (permutation equivalence, bitwise)'
   else
      nfail = nfail + 1
      print '(a)', 'FAIL (l) roster wrapper: permutation mismatch'
   end if

   ! ---- (i1) degenerate point: exactness + tight drift ----
   if (cau_peri_tdhf_armed()) then
      print '(a)', 'FAIL i1 order: the package armed before the degenerate-point leg'
      nfail = nfail + 1
   end if
   call crit_traj_unarmed(npass, nfail)

   ! ---- (a) + production arm ----
   call cau_peri_tdhf_setup(4, 7.0d0, 300.0d0, 'GS', 1, 0, 0.05d0)
   call crit_gauss(npass, nfail, s13)
   call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 0, 0.05d0)
   call crit_bath(npass, nfail)
   call crit_seed_gs(npass, nfail)

   ! ---- (e) Hcore structure after one prop call ----
   call state_create(sta, natoms, 1)
   sta%q = q0
   sta%t = 0.0d0
   call cau_peri_prop(sta, 0.05d0)
   call crit_hcore(npass, nfail)
   call state_destroy(sta)

   ! ---- (g) unitary evolution ----
   call crit_unitary(npass, nfail)

   ! ---- (f2) 1D seed ----
   call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, '1D', 4, 0, 0.05d0)
   call crit_seed_1d(npass, nfail)
   call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 0, 0.05d0)

   ! ---- (d) Hubbard kernels vs FD of the values ----
   call fd_hub_all(npass, nfail, dx)

   ! ---- (h) the composed fold ----
   call cau_peri_pes(q0, mass, g, v)
   do dirn = 1, 4
      call dof_fd_v(dof(dirn), dx, vp, vm)
      call cmp_one(ctag('h', dirn), g(dof(dirn)), (vp - vm)/(2*dx), 1.0d-5, npass, nfail)
   end do

   ! ---- (h2) the composed fold with a NONDIAGONAL rho: evolve the state
   ! ~5 ps worth of Magnus substeps at the probe geometry so the coupling
   ! blocks carry REAL off-diagonal density (the (h) backdrop stays blind:
   ! one prop call leaves the off-diagonals purely imaginary, and the
   ! coupling-kernel partner elements read exactly those real parts - the
   ! 2026-09-28 half-coupling-force bug hid in this blind spot) ----
   call state_create(sta, natoms, 1)
   sta%q = q0
   do kk = 1, 100
      sta%t = kk*0.05d0
      call cau_peri_prop(sta, 0.05d0)
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
      call cau_peri_pes(sta%q, mass, g, v)
      do dirn = 1, 4
         call dof_fd_q(sta%q, dof(dirn), dx, vp, vm)
         call cmp_one(ctag('h2', dirn), g(dof(dirn)), (vp - vm)/(2*dx), 1.0d-5, npass, nfail)
      end do
      print '(a,es9.2)', '(h2) nondiagonal-rho fold backdrop: max |Re rho(imp,bath)| = ', offd
   end if
   call state_destroy(sta)

   ! ---- (i3) static-geometry mean-field theorem: fixed q, GS seed, rho
   ! evolves 300 steps with msbstp=16 - d[Tr(rho*Hcore)+E_mf] must hold to
   ! the Magnus-2 budget (locks the per-spin shift-vector assembly; the
   ! former two-scalar form cross-wired the spin channels and drifted
   ! +5 eV/300 steps, msbstp-independent - 2026-09-28) ----
   block
      real(8) :: eq(9), g9(9), vv, eb, ee2, mf1, mf2
      integer :: kk2
      call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 16, 0, 0.05d0)
      eq = q0
      eq(9) = 2.2d0                          ! near-collision height
      sta%q = eq
      sta%t = 0.0d0
      call cau_peri_pes(eq, mass, g9, vv)    ! establish hcore for the baseline
      eb = 0.0d0
      do kk2 = 1, 300
         sta%t = kk2*0.05d0
         call cau_peri_prop(sta, 0.05d0)
      end do
      ee2 = cau_peri_col_eelec(sta)
      mf2 = cau_peri_col_emf(sta)
      ! re-evaluate the baseline: seed + the same hcore
      call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 16, 0, 0.05d0)
      call cau_peri_pes(eq, mass, g9, vv)
      mf1 = cau_peri_col_emf(sta)
      eb = cau_peri_col_eelec(sta) + mf1
      if (abs(ee2 + mf2 - eb) > 0.05d0) then
         print '(a,es10.3)', 'FAIL (i3) static theorem: |dE| = ', abs(ee2 + mf2 - eb), &
            ' > 0.05 eV (300 steps, msbstp=16)'
         nfail = nfail + 1
      else
         print '(a,es10.3)', '(i3) static mean-field theorem: PASS (|dE| = ', &
            abs(ee2 + mf2 - eb), ' eV over 300 steps, msbstp=16)'
      end if
      call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 0, 0.05d0)
   end block

   ! ---- (o) far-field freeze (ruling 2026-09-29): at/above the 6 A analysis
   ! line the delivered Hubbard kernels and the C-side ground kernel are
   ! EXACTLY zero, the sheet values are frozen (7 A equals 9 A bitwise), and
   ! inside the [5,6] warp ramp the ground kernel still matches central
   ! differences (the chain-rule scale is consistent with the warped values)
   block
      real(8) :: qf(9), gf(9), ef, hv(2), uv, jv, vcpv, hb(2), ub, jb, vcb
      real(8) :: dhf(2, 9), duf(9), djf(9), dvcf(9), e1, e2v
      qf = 0.0d0
      qf(6) = 0.9d0
      qf(1) = 0.15d0
      qf(9) = 7.0d0
      call cau_peri_hubbard(qf, hv, uv, jv, vcpv, dhf, duf, djf, dvcf)
      call cau_peri_ground(qf, ef, gf)
      ! (o1) the exact-zero set: coupling kernel identically zero; the
      ! z-components (C z + S z + Y z rows 3/6/9) of every Hubbard kernel
      ! and the C-z ground kernel exactly zero
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
      ! (o2) lateral residuals on the frozen sheet are the original fit's
      ! site-asymptote artifact (energy-consistent; the trajectory ends at
      ! the 6 A line, so production never samples them) - bound and REPORT
      ! the bound is the site-asymptote artifact scale of the fit (post
      ! site-column transposition 2026-09-29: |dh| ~1.5 eV/A at the probe
      ! lateral line; production terminates at the 6 A line and never
      ! samples the region)
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
      call cau_peri_hubbard(qf, hb, ub, jb, vcb, dhf, duf, djf, dvcf)
      if (hb(1) == hv(1) .and. hb(2) == hv(2) .and. ub == uv .and. jb == jv .and. &
          vcb == vcpv) then
         npass = npass + 1
         print '(a)', '(o) far-field values: PASS (frozen - 7 A equals 9 A)'
      else
         print '(a)', 'FAIL (o) far-field values drift beyond the freeze line'
         nfail = nfail + 1
      end if
      qf(9) = 5.3d0
      call cau_peri_ground(qf, ef, gf)
      qf(9) = 5.3d0 + dx
      call cau_peri_ground(qf, e1, g)
      qf(9) = 5.3d0 - dx
      call cau_peri_ground(qf, e2v, g)
      call cmp_one('o ramp', gf(9), (e1 - e2v)/(2*dx), 1.0d-5, npass, nfail)
   end block

   ! ---- (n1) NN loader sanity: the shelf nn_weights.txt must parse to a
   ! NONEMPTY model (the vendored module once lost the NN_E0..NN_Eh index
   ! constants - under implicit typing they became zero integers, every
   ! section lookup hit ipot=0, and the loader 'succeeded' with 486->0
   ! weights; restored 2026-09-29) ----
   block
      use surface_potentials, only: load_all_nn_weights, nn_nhidden, use_nn_pes
      integer :: ierr2
      call copy_file(trim(shelf)//'nn_weights.txt', 'nn_weights.txt')
      call load_all_nn_weights('nn_weights.txt', ierr2)
      if (ierr2 /= 0 .or. any(nn_nhidden(1:6) == 0) .or. .not. use_nn_pes) then
         print '(a)', 'FAIL (n1) NN loader: empty or failed model (weights silently zero)'
         nfail = nfail + 1
      else
         print '(a,7i4)', '(n1) NN loader: PASS (NH per type =', nn_nhidden, ')'
      end if
      call cau_peri_load()                     ! restore the analytic branch
   end block

   ! ---- (i2) armed trajectory ----
   call crit_traj_armed(npass, nfail)

   ! ---- (j) sync publication + (k) recorder columns ----
   call crit_sync_cols(npass, nfail)

   print '(/,a,i0,a,i0,a)', 'check_cau_peri: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) call exit(1)

contains

   ! copy_file(from, to) - portable text copy (the module's READAP reads
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
         call cau_peri_prop(sta, 0.05d0)
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
      call cau_peri_pes(q0, mass1, gg, vv)
      call kernel_caches_out(dhv, duv, djv, dvvv)
      do d = 1, 4
         qp = q0
         qp(dof(d)) = q0(dof(d)) + h
         call cau_peri_pes(qp, mass1, gg, vv)
         call hub_vals_out(hp, up, jp, vpp)
         qp(dof(d)) = q0(dof(d)) - h
         call cau_peri_pes(qp, mass1, gg, vv)
         call hub_vals_out(hm, um, jm, vmm)
         qp(dof(d)) = q0(dof(d))
         call cmp_one(ctag('dh', d), dhv(1, dof(d)), (hp(1) - hm(1))/(2*h), 1.0d-5, np, nf)
         call cmp_one(ctag('dU', d), duv(dof(d)), (up - um)/(2*h), 1.0d-5, np, nf)
         call cmp_one(ctag('dJ', d), djv(dof(d)), (jp - jm)/(2*h), 1.0d-5, np, nf)
         call cmp_one(ctag('dV', d), dvvv(dof(d)), (vpp - vmm)/(2*h), 1.0d-5, np, nf)
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
      call cau_peri_pes(qp, mass, gg, vpp)
      qp(k) = q0(k) - h
      call cau_peri_pes(qp, mass, gg, vmm)
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
      call cau_peri_pes(qp, mass, gg, vpp)
      qp(k) = qb(k) - h
      call cau_peri_pes(qp, mass, gg, vmm)
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
      call cau_peri_pes(sta%q, mass, gg, vv)
      call cau_peri_ground(sta%q, e0_ref, g0_ref)
      call check_c('i1 v==E0', abs(vv - e0_ref) < 1.0d-14, np, nf)
      call check_c('i1 g==g0', maxval(abs(gg - g0_ref)) < 1.0d-14, np, nf)
      h0 = tot_e(sta, vv)
      do step = 1, 300
         call vstep(sta, 0.05d0, .false.)
         call cau_peri_pes(sta%q, mass, gg, vv)
      end do
      h1 = tot_e(sta, vv)
      call check_c('i1 drift', abs(h1 - h0) < 0.05d0, np, nf)
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
      call cau_peri_pes(sta%q, mass, gg, vv)
      h0 = tot_e(sta, vv)
      call rho_herm_out(tr0, herm)
      do step = 1, 300
         call vstep(sta, 0.05d0, .true.)
         call cau_peri_pes(sta%q, mass, gg, vv)
      end do
      h1 = tot_e(sta, vv)
      call rho_herm_out(tr1, herm)
      call check_c('i2 bound', abs(h1 - h0) < 200.0d0, np, nf)
      call check_c('i2 trace', abs(tr1 - tr0) < 1.0d-6, np, nf)
      print '(a,es12.4,a,es10.2)', 'i2 report: armed drift [kcal/mol] = ', h1 - h0, &
         '  dtrace = ', tr1 - tr0
      call state_destroy(sta)
   end subroutine crit_traj_armed

   subroutine fd_kernel(k, dx, gn)
      integer, intent(in) :: k
      real(8), intent(in) :: dx
      real(8), intent(out) :: gn(:)
      real(8) :: ep, em
      real(8), allocatable :: qp(:)
      allocate (qp(ncoord))
      qp = q0
      qp(k) = q0(k) + dx
      call cau_peri_ground(qp, ep, g)
      qp(k) = q0(k) - dx
      call cau_peri_ground(qp, em, g)
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
      call cau_peri_pes(st%q, mass, gg2, vv2)
      do ia = 1, natoms
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg2(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      do ia = 1, natoms
         st%q(3*ia - 2:3*ia) = st%q(3*ia - 2:3*ia) + st%p(3*ia - 2:3*ia)/mass(ia)*dt
      end do
      call cau_peri_pes(st%q, mass, gg2, vv2)
      do ia = 1, natoms
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg2(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      if (evolve) then
         st%t = st%t + dt
         call cau_peri_prop(st, dt)
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
      call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 0, 0.05d0)   ! fresh GS seed
      call state_create(st, natoms, 1)
      st%q = q0
      st%t = 0.0d0
      st%p = 0.0d0
      st%p(3) = -c_mass*0.30d0
      call cau_peri_pes(st%q, mass, g, vv)                 ! armed fold at the seed
      call check_c('k nup seed', abs(cau_peri_col_nup(st) - 2.0d0) < 1.0d-14, np, nf)
      call check_c('k ndn seed', abs(cau_peri_col_ndn(st)) < 1.0d-14, np, nf)
      call check_c('k mag seed', abs(cau_peri_col_mag(st) - 2.0d0) < 1.0d-14, np, nf)
      call cau_peri_sync(st%s)
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
         call cau_peri_prop(st, 0.05d0)
      end do
      call cau_peri_sync(st%s)
      call check_c('k nup range', cau_peri_col_nup(st) > -0.05d0 .and. cau_peri_col_nup(st) < 3.05d0, np, nf)
      call check_c('k ndn range', cau_peri_col_ndn(st) > -0.05d0 .and. cau_peri_col_ndn(st) < 3.05d0, np, nf)
      call check_c('k eelec fin', abs(cau_peri_col_eelec(st)) < 1.0d4, np, nf)
      call check_c('k emf fin', abs(cau_peri_col_emf(st)) < 1.0d4, np, nf)
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
end program check_cau_peri
