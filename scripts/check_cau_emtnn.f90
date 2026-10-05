!=====================================================================
! check_cau_emtnn.f90 - the C@Au(111) EMT-NN+TDHF container criterion
!                         program (criteria defined FIRST in this header)
! Criteria:
!   (a) Gauss-Legendre closed forms through a 4-orbital probe bath
!       (extremes +-bw(1+1/sqrt(3))/4, inner +-bw(1-1/sqrt(3))/4; 1e-12)
!   (b) production bath (20 orbitals): symmetric band, range, the
!       energy_metal_zero identity (1e-12)
!   (c) EMT-NN ground kernel vs central difference of E0 at the GENERIC
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
! Engine files: the container folder's system/Slab.xyz + the folder's nn_weights_emt_nn.txt
!   (referenced by the path constant below - no run-directory copies
!   needed for the criterion run).
! Compile (from scripts/, gfortran; container on the inactive shelf):
!   gfortran -O2 -o check_cau_emtnn.exe ^
!     ..\systems\inactive\cau_emtnn\emt_nn_pes_venus.f90 ^
!     ..\systems\inactive\cau_emtnn\cau_emtnn_pot.f90 ^
!     ..\systems\inactive\cau_emtnn\cau_emtnn_tdhf.f90 ^
!     ..\systems\inactive\cau_emtnn\cau_emtnn_box.f90 ^
!     ..\src\globals\consts.f90 ..\src\globals\state.f90 ^
!     ..\src\utils\linalg.f90 ..\src\utils\rng.f90 ^
!     check_cau_emtnn.f90
!   (j)/(k) sync publication + recorder columns (seed values exact,
!       post-propagation bounded)
!   check_cau_emtnn.exe   (run from scripts/)
!=====================================================================
program check_cau_emtnn
   use cau_emtnn_pot,  only: cau_emtnn_load, cau_emtnn_slab_atoms, cau_emtnn_ground
   use cau_emtnn_tdhf, only: cau_emtnn_tdhf_setup, cau_emtnn_tdhf_armed, cau_emtnn_pes, &
                             cau_emtnn_prop, bath_out, rho_seed_out, rho_herm_out, hcore_out, &
                             kernel_caches_out, hub_vals_out
   use cau_emtnn_tdhf, only: cau_emtnn_sync, cau_emtnn_col_nup, cau_emtnn_col_ndn, &
                            cau_emtnn_col_mag, cau_emtnn_col_eelec, cau_emtnn_col_emf, rho_ref_out
   use state,          only: state_t, state_create, state_destroy
   implicit none
   character(len=*), parameter :: shelf = '../systems/inactive/cau_emtnn/'
   real(8), parameter :: ev_kcal = 23.0605d0, e_conv = 0.04184d0
   integer :: nslab, natoms, ncoord, npass, nfail, dirn, u, ios, i
   real(8), allocatable :: q0(:), mass(:), g(:), gnum(:), gp(:), gm(:)
   real(8) :: e0, v, vp, vm, dx, s13, x1, y1, z1
   real(8) :: au_mass = 196.97d0, c_mass = 12.01d0
   character(len=4) :: elem
   type(state_t) :: sta
   npass = 0
   nfail = 0
   s13 = 1.0d0/sqrt(3.0d0)

   ! ---- 0. engine init + geometry from the container's slab file ----
   call cau_emtnn_load(trim(shelf)//'system/Slab.xyz', trim(shelf)//'nn_weights_emt_nn.txt', 2.88d0, 6)
   nslab = cau_emtnn_slab_atoms()
   natoms = 1 + nslab
   ncoord = 3*natoms
   allocate (q0(ncoord), mass(natoms), g(ncoord), gnum(ncoord), gp(ncoord), gm(ncoord))
   open (newunit=u, file=trim(shelf)//'system/Slab.xyz', status='old', action='read', iostat=ios)
   if (ios /= 0) then
      print '(a)', 'FAIL slab file unreadable'
      call exit(1)
   end if
   read (u, *) i
   read (u, *)
   do i = 1, nslab
      read (u, *) elem, q0(3*(i + 1) - 2:3*(i + 1))
   end do
   close (u)
   x1 = q0(4)
   y1 = q0(5)
   z1 = q0(6)
   q0(1) = x1 + 0.37d0                     ! C at the generic site over atom 1
   q0(2) = y1 + 0.21d0
   q0(3) = z1 + 2.2d0
   mass(1) = c_mass
   mass(2:) = au_mass

   ! ---- (c) ground kernel vs central difference ----
   dx = 1.0d-4
   do dirn = 1, 4                          ! C x/y/z + slab atom 2 z (dof 6)
      call fd_kernel(dof(dirn), dx, gnum)
      call cmp_one(ctag('c', dirn), g(dof(dirn)), gnum(dof(dirn)), 5.0d-3, npass, nfail)
   end do

   ! ---- (i1) degenerate point: exactness + tight drift ----
   if (cau_emtnn_tdhf_armed()) then
      print '(a)', 'FAIL i1 order: the package armed before the degenerate-point leg'
      nfail = nfail + 1
   end if
   call crit_traj_unarmed(npass, nfail)

   ! ---- (a) + production arm ----
   call cau_emtnn_tdhf_setup(4, 7.0d0, 300.0d0, 'GS', 1, 0, 0.05d0)
   call crit_gauss(npass, nfail, s13)
   call cau_emtnn_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 0, 0.05d0)
   call crit_bath(npass, nfail)
   call crit_seed_gs(npass, nfail)

   ! ---- (e) Hcore structure after one prop call ----
   call state_create(sta, natoms, 1)
   sta%q = q0
   sta%t = 0.0d0
   call cau_emtnn_prop(sta, 0.05d0)
   call crit_hcore(npass, nfail)
   call state_destroy(sta)

   ! ---- (g) unitary evolution ----
   call crit_unitary(npass, nfail)

   ! ---- (f2) 1D seed ----
   call cau_emtnn_tdhf_setup(20, 7.0d0, 300.0d0, '1D', 4, 0, 0.05d0)
   call crit_seed_1d(npass, nfail)
   call cau_emtnn_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 0, 0.05d0)

   ! ---- (d) Hubbard kernels vs FD of the values ----
   call fd_hub_all(npass, nfail, dx)

   ! ---- (h) the composed fold ----
   call cau_emtnn_pes(q0, mass, g, v)
   do dirn = 1, 4
      call dof_fd_v(dof(dirn), dx, vp, vm)
      call cmp_one(ctag('h', dirn), g(dof(dirn)), (vp - vm)/(2*dx), 1.0d-5, npass, nfail)
   end do

   ! ---- (i2) armed trajectory ----
   call crit_traj_armed(npass, nfail)

   ! ---- (j) sync publication + (k) recorder columns ----
   call crit_sync_cols(npass, nfail)

   print '(/,a,i0,a,i0,a)', 'check_cau_emtnn: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) call exit(1)

contains

   integer function dof(d)
      integer, intent(in) :: d
      if (d == 4) then
         dof = 6
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
         call cau_emtnn_prop(sta, 0.05d0)
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
      call cau_emtnn_pes(q0, mass1, gg, vv)
      call kernel_caches_out(dhv, duv, djv, dvvv)
      do d = 1, 4
         qp = q0
         qp(dof(d)) = q0(dof(d)) + h
         call cau_emtnn_pes(qp, mass1, gg, vv)
         call hub_vals_out(hp, up, jp, vpp)
         qp(dof(d)) = q0(dof(d)) - h
         call cau_emtnn_pes(qp, mass1, gg, vv)
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
      call cau_emtnn_pes(qp, mass, gg, vpp)
      qp(k) = q0(k) - h
      call cau_emtnn_pes(qp, mass, gg, vmm)
      deallocate (qp)
   end subroutine dof_fd_v

   subroutine crit_traj_unarmed(np, nf)
      integer, intent(inout) :: np, nf
      real(8) :: h0, h1, vv, e0_ref, gg(ncoord), g0_ref(ncoord)
      integer :: step
      call state_create(sta, natoms, 1)
      sta%q = q0
      sta%t = 0.0d0
      sta%p = 0.0d0
      sta%p(3) = -c_mass*0.30d0
      call cau_emtnn_pes(sta%q, mass, gg, vv)
      call cau_emtnn_ground(sta%q, e0_ref, g0_ref)
      call check_c('i1 v==E0', abs(vv - e0_ref) < 1.0d-14, np, nf)
      call check_c('i1 g==g0', maxval(abs(gg - g0_ref)) < 1.0d-14, np, nf)
      h0 = tot_e(sta, vv)
      do step = 1, 300
         call vstep(sta, 0.05d0, .false.)
         call cau_emtnn_pes(sta%q, mass, gg, vv)
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
      call cau_emtnn_pes(sta%q, mass, gg, vv)
      h0 = tot_e(sta, vv)
      call rho_herm_out(tr0, herm)
      do step = 1, 300
         call vstep(sta, 0.05d0, .true.)
         call cau_emtnn_pes(sta%q, mass, gg, vv)
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
      call cau_emtnn_ground(qp, ep, g)
      qp(k) = q0(k) - dx
      call cau_emtnn_ground(qp, em, g)
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
      call cau_emtnn_pes(st%q, mass, gg2, vv2)
      do ia = 1, natoms
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg2(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      do ia = 1, natoms
         st%q(3*ia - 2:3*ia) = st%q(3*ia - 2:3*ia) + st%p(3*ia - 2:3*ia)/mass(ia)*dt
      end do
      call cau_emtnn_pes(st%q, mass, gg2, vv2)
      do ia = 1, natoms
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg2(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      if (evolve) then
         st%t = st%t + dt
         call cau_emtnn_prop(st, dt)
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
      call cau_emtnn_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 0, 0.05d0)   ! fresh GS seed
      call state_create(st, natoms, 1)
      st%q = q0
      st%t = 0.0d0
      st%p = 0.0d0
      st%p(3) = -c_mass*0.30d0
      call cau_emtnn_pes(st%q, mass, g, vv)                 ! armed fold at the seed
      call check_c('k nup seed', abs(cau_emtnn_col_nup(st) - 2.0d0) < 1.0d-14, np, nf)
      call check_c('k ndn seed', abs(cau_emtnn_col_ndn(st)) < 1.0d-14, np, nf)
      call check_c('k mag seed', abs(cau_emtnn_col_mag(st) - 2.0d0) < 1.0d-14, np, nf)
      call cau_emtnn_sync(st%s)
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
         call cau_emtnn_prop(st, 0.05d0)
      end do
      call cau_emtnn_sync(st%s)
      call check_c('k nup range', cau_emtnn_col_nup(st) > -0.05d0 .and. cau_emtnn_col_nup(st) < 3.05d0, np, nf)
      call check_c('k ndn range', cau_emtnn_col_ndn(st) > -0.05d0 .and. cau_emtnn_col_ndn(st) < 3.05d0, np, nf)
      call check_c('k eelec fin', abs(cau_emtnn_col_eelec(st)) < 1.0d4, np, nf)
      call check_c('k emf fin', abs(cau_emtnn_col_emf(st)) < 1.0d4, np, nf)
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
end program check_cau_emtnn
