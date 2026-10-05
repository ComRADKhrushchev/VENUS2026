!=====================================================================
! cau_rst_pot.f90 - C@Au(111) RST-Hubbard ground PES + Hubbard delivery
!   (pure physics layer; q in, eV out, zero program knowledge)
! Design:
!   Ground chain: Born-von Karman slab elasticity (12-nearest-neighbor
!   harmonic tensors from the BVK parameter line, built once over the
!   file slab) + the GS-state pairwise NN (C-slab 1D-NN + Born-Mayer
!   repulsive + smooth cutoff) + the weights-file baseline.
!   Hubbard chain (direct mode): h = NN_h (linear), U = softplus(NN_U),
!   J = softplus(NN_J) clamped to [0, 0.95U], Vcp = NN_Vcp (linear);
!   gradients through calc_rst_energy_grad chains (sigmoid for the
!   softplus pair, saturation gate on the J clamp).
! DEVIATIONS from the vendored interface (2026-09-21, all criterion-
!   verified by check_cau_rst legs (c)/(d)/(h)):
!   1. SLAB SYNC: the old interface passed the module slab to the NN
!      engine without syncing trajectory positions (the NN saw the
!      frozen file geometry; only the BvK term followed slab motion).
!      This container syncs slab%coords from q before every engine call
!      - the NN pairs see the movable slab.
!   2. GRADIENT CHAINS RESTORED: the old dH_orb shipped with debug
!      flags zero_dh = .true., zero_dU = .true. (Test E/F scaffolding)
!      and dW carried a stray 0.99 off-diagonal factor; its dV_orb
!      filled only the first 9 dofs (a 3-atom-interface copy-paste
!      residue). This container delivers the full consistent chains.
!   3. J-CLAMP GRADIENT: the clamp [0, 0.95U] saturates dJ inside the
!      clamped region (0 there), which the old code never expressed.
! Units: eV / A throughout; the pairwise-NN wrap is orthogonal
!   (box = slab span + 6 A padding, the engine's own convention).
!=====================================================================
module cau_rst_pot
   use rst_hubbard_pes
   implicit none
   private
   public :: cau_rst_load, cau_rst_slab_atoms
   public :: cau_rst_ground, cau_rst_hubbard, cau_rst_hubbard_val
   ! fit interface (probe_fit, 2026-09-25): expose the six loaded state
   ! vectors so an external LM driver can rewrite rst_hub(k)%w in place;
   ! slab / BvK / baseline module state stays as loaded by cau_rst_load
   public :: rst_hub, have_hub

   type(rst_slab_data) :: slab                 ! engine slab (synced from q per call)
   type(rst_state) :: rst_hub(6)               ! 1=h 2=U 3=J 4=(unused) 5=Vcp 6=GS
   logical :: have_hub(6) = .false.
   real(8) :: nn_baseline = 0.0d0
   logical :: loaded = .false.

   ! -- BvK slab elasticity state (built once over the file slab) --
   real(8) :: bvk_cutoff = 6.0d0, bvk_a = 2.95d0
   real(8) :: bvk_alpha = 0.5d0, bvk_beta = 0.3d0, bvk_gamma = 0.1d0
   real(8), allocatable :: eq_slab(:,:), bvktensor(:,:,:,:)
   integer, allocatable :: slab_pairs(:,:)
   integer :: n_slab = 0

contains

   !------------------------------------------------------------------
   ! dist_pbc(dx, dy) - orthogonal minimum image with the slab box
   !------------------------------------------------------------------
   subroutine dist_pbc(dx, dy)
      real(8), intent(inout) :: dx, dy
      dx = dx - slab%box_lx*nint(dx/slab%box_lx)
      dy = dy - slab%box_ly*nint(dy/slab%box_ly)
   end subroutine dist_pbc

   !------------------------------------------------------------------
   ! cau_rst_load(slab_file, weights_file, bvk_file) - slab + weights +
   !                 BVK parameters + one-time pair/tensor build
   !------------------------------------------------------------------
   subroutine cau_rst_load(slab_file, weights_file, bvk_file)
      character(len=*), intent(in) :: slab_file, weights_file
      character(len=*), intent(in), optional :: bvk_file
      integer :: u, ios, i, j, icnt
      real(8) :: rij(3), tensor(3,3)
      logical :: found
      character(len=256) :: line
      call load_slab_xyz(slab_file, slab)
      n_slab = slab%n_atoms
      call load_rst_weights_standalone(weights_file, rst_hub, have_hub, nn_baseline)
      ! BVK parameters: the "BVK:" line of Analytic_Potential.txt, defaults below
      found = .false.
      if (present(bvk_file)) then
         open (newunit=u, file=bvk_file, status='old', action='read', iostat=ios)
         if (ios == 0) then
            do
               read (u, '(a)', iostat=ios) line
               if (ios /= 0) exit
               if (index(line, 'BVK:') > 0) then
                  read (line(index(line, 'BVK:') + 4:), *, iostat=ios) &
                     bvk_cutoff, bvk_a, bvk_alpha, bvk_beta, bvk_gamma
                  if (ios == 0) found = .true.
                  exit
               end if
            end do
            close (u)
         end if
      end if
      if (.not. found) then
         bvk_cutoff = 6.0d0
         bvk_a = 2.95d0
         bvk_alpha = 0.5d0
         bvk_beta = 0.3d0
         bvk_gamma = 0.1d0
      end if
      ! pair list + tensors over the file geometry (equilibrium positions)
      allocate (slab_pairs(n_slab, 12), eq_slab(n_slab, 3), bvktensor(n_slab, 12, 3, 3))
      slab_pairs = 0
      eq_slab = slab%coords
      do i = 1, n_slab
         icnt = 0
         do j = i + 1, n_slab
            rij = eq_slab(i, :) - eq_slab(j, :)
            call dist_pbc(rij(1), rij(2))
            if (dot_product(rij, rij) > (bvk_a + 0.1d0)**2 .or. icnt >= 12) cycle
            icnt = icnt + 1
            slab_pairs(i, icnt) = j
            call build_tensor(i, j, bvktensor(i, icnt, :, :))
         end do
      end do
      loaded = .true.
   end subroutine cau_rst_load

   integer function cau_rst_slab_atoms()
      cau_rst_slab_atoms = n_slab
   end function cau_rst_slab_atoms

   !------------------------------------------------------------------
   ! build_tensor(i, j, t) - FCC bond-direction force-constant tensor
   !                 (the old BuildTensor logic, module-state parameters)
   !------------------------------------------------------------------
   subroutine build_tensor(i, j, t)
      integer, intent(in) :: i, j
      real(8), intent(out) :: t(3,3)
      real(8) :: rij(3), r_dir(3), a_base(3,3,6), fcc2slab(3,3), base_dir(3,6)
      real(8) :: a_, paral
      integer :: ib
      real(8) :: asqr2, asqr3, asqr6
      asqr2 = 1.0d0/sqrt(2.0d0)
      asqr3 = 1.0d0/sqrt(3.0d0)
      asqr6 = 1.0d0/sqrt(6.0d0)
      t = 0.0d0
      fcc2slab = reshape([ &
         asqr2, asqr2, 0.0d0, &
         asqr6, -asqr6, -2*asqr6, &
         asqr3, -asqr3, asqr3], [3,3], order=[2,1])
      rij = eq_slab(i, :) - eq_slab(j, :)
      call dist_pbc(rij(1), rij(2))
      a_ = sqrt(dot_product(rij, rij))
      r_dir = rij/a_
      r_dir = matmul(transpose(fcc2slab), r_dir)
      base_dir(:,1) = [0.0d0, asqr2, asqr2]
      a_base(:,:,1) = reshape([bvk_alpha, 0.0d0, 0.0d0, 0.0d0, bvk_beta, bvk_gamma, &
                               0.0d0, bvk_gamma, bvk_beta], [3,3])
      base_dir(:,2) = [0.0d0, asqr2, -asqr2]
      a_base(:,:,2) = reshape([bvk_alpha, 0.0d0, 0.0d0, 0.0d0, bvk_beta, -bvk_gamma, &
                               0.0d0, -bvk_gamma, bvk_beta], [3,3])
      base_dir(:,3) = [asqr2, 0.0d0, asqr2]
      a_base(:,:,3) = reshape([bvk_beta, 0.0d0, bvk_gamma, 0.0d0, bvk_alpha, 0.0d0, &
                               bvk_gamma, 0.0d0, bvk_beta], [3,3])
      base_dir(:,4) = [-asqr2, 0.0d0, asqr2]
      a_base(:,:,4) = reshape([bvk_beta, 0.0d0, -bvk_gamma, 0.0d0, bvk_alpha, 0.0d0, &
                               -bvk_gamma, 0.0d0, bvk_beta], [3,3])
      base_dir(:,5) = [asqr2, asqr2, 0.0d0]
      a_base(:,:,5) = reshape([bvk_beta, bvk_gamma, 0.0d0, bvk_gamma, bvk_beta, 0.0d0, &
                               0.0d0, 0.0d0, bvk_alpha], [3,3])
      base_dir(:,6) = [asqr2, -asqr2, 0.0d0]
      a_base(:,:,6) = reshape([bvk_beta, -bvk_gamma, 0.0d0, -bvk_gamma, bvk_beta, 0.0d0, &
                               0.0d0, 0.0d0, bvk_alpha], [3,3])
      do ib = 1, 6
         paral = dot_product(base_dir(:,ib), r_dir)
         if (abs(abs(paral) - 1) < 1e-3) then
            t = matmul(fcc2slab, matmul(a_base(:,:,ib), transpose(fcc2slab)))
            exit
         end if
      end do
   end subroutine build_tensor

   !------------------------------------------------------------------
   ! sync_slab(q) - copy trajectory slab positions into the engine slab
   !------------------------------------------------------------------
   subroutine sync_slab(q)
      real(8), intent(in) :: q(:)
      integer :: k
      do k = 1, n_slab
         slab%coords(1:3, k) = q(3*k + 1:3*k + 3)
      end do
   end subroutine sync_slab

   !------------------------------------------------------------------
   ! cau_rst_ground(q, e0, g0) - BvK elastic + GS pairwise NN + baseline
   !------------------------------------------------------------------
   subroutine cau_rst_ground(q, e0, g0)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: e0
      real(8), intent(out) :: g0(:)
      real(8) :: r_vec(3), r_eq(3), r_pt(3), t(3,3), e_nn
      real(8) :: dEdx, dEdy, dEdz
      real(8), allocatable :: dEdR(:,:)
      integer :: i, k, i1, i2, natoms
      natoms = size(q)/3
      if (.not. loaded) then
         call stop_pot('cau_rst_ground', 'the engine is not initialized (seam order error)')
      end if
      if (natoms /= 1 + n_slab) then
         call stop_pot('cau_rst_ground', 'the geometry size does not match the engine slab')
      end if
      call sync_slab(q)
      ! 1. BvK harmonic over the 12-neighbor pair list (built at eq geometry)
      e0 = 0.0d0
      g0 = 0.0d0
      do i = 1, n_slab
         do k = 1, 12
            i2 = slab_pairs(i, k)
            if (i2 == 0 .or. i2 == 1) cycle
            i1 = i + 1                          ! +1: atom 1 is the adsorbate
            r_vec = q(3*i1 - 2:3*i1) - q(3*i2 - 2:3*i2)
            r_eq = eq_slab(i, :) - eq_slab(i2, :)
            call dist_pbc(r_vec(1), r_vec(2))
            call dist_pbc(r_eq(1), r_eq(2))
            r_pt = r_vec - r_eq
            t = bvktensor(i, k, :, :)
            e0 = e0 + 0.5d0*dot_product(matmul(t, r_pt), r_pt)
            g0(3*i1 - 2:3*i1) = g0(3*i1 - 2:3*i1) + matmul(t, r_pt)
            g0(3*i2 - 2:3*i2) = g0(3*i2 - 2:3*i2) - matmul(t, r_pt)
         end do
      end do
      ! 2. GS-state pairwise NN (+ Born-Mayer repulsive + baseline)
      allocate (dEdR(3, n_slab))
      call calc_rst_energy_grad(q(1:3), slab, rst_hub(6), e_nn, dEdx, dEdy, dEdz, dEdR=dEdR)
      e0 = e0 + e_nn + nn_baseline
      g0(1:3) = g0(1:3) + [dEdx, dEdy, dEdz]
      do i = 1, n_slab
         g0(3*(i + 1) - 2:3*(i + 1)) = g0(3*(i + 1) - 2:3*(i + 1)) - dEdR(1:3, i)
      end do
      deallocate (dEdR)
   end subroutine cau_rst_ground

   !------------------------------------------------------------------
   ! cau_rst_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp) - direct-mode
   !                 delivery with the restored gradient chains
   !------------------------------------------------------------------
   subroutine cau_rst_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: h(2), u, j, vcp
      real(8), intent(out) :: dh(2, size(q)), du(size(q)), dj(size(q)), dvcp(size(q))
      real(8) :: h_raw, u_raw, j_raw, sig_u, sig_j
      real(8) :: e_t, dhx, dhy, dhz, dux, duy, duz, djx, djy, djz
      real(8), allocatable :: dr_h(:,:), dr_u(:,:), dr_j(:,:), dr_v(:,:)
      real(8) :: u_raw_at, j_raw_at
      integer :: k
      call sync_slab(q)
      ! 1. raw NN values + the softplus/clamp activations
      call calc_rst_energy(q(1:3), slab, rst_hub(1), h_raw)
      call calc_rst_energy(q(1:3), slab, rst_hub(2), u_raw)
      call calc_rst_energy(q(1:3), slab, rst_hub(3), j_raw)
      call calc_rst_energy(q(1:3), slab, rst_hub(5), vcp)
      u = softplus(u_raw)
      j = softplus(j_raw)
      j = min(max(j, 0.0d0), 0.95d0*u)
      h(1) = h_raw
      h(2) = h_raw
      ! 2. raw NN gradients + activation chains
      allocate (dr_h(3, n_slab), dr_u(3, n_slab), dr_j(3, n_slab), dr_v(3, n_slab))
      call calc_rst_energy_grad(q(1:3), slab, rst_hub(1), e_t, dhx, dhy, dhz, dEdR=dr_h)
      call calc_rst_energy_grad(q(1:3), slab, rst_hub(2), e_t, dux, duy, duz, dEdR=dr_u)
      call calc_rst_energy_grad(q(1:3), slab, rst_hub(3), e_t, djx, djy, djz, dEdR=dr_j)
      sig_u = sigmoid(u_raw)
      j_raw_at = softplus(j_raw)                ! pre-clamp activation
      if (j_raw_at <= 0.0d0 .or. j_raw_at >= 0.95d0*u) then
         sig_j = 0.0d0                          ! clamp saturation: dJ = 0
      else
         sig_j = sigmoid(j_raw)
      end if
      dh = 0.0d0
      du = 0.0d0
      dj = 0.0d0
      dvcp = 0.0d0
      dh(1, 1:3) = [dhx, dhy, dhz]
      dh(2, 1:3) = [dhx, dhy, dhz]
      du(1:3) = sig_u*[dux, duy, duz]
      dj(1:3) = sig_j*[djx, djy, djz]
      call calc_rst_energy_grad(q(1:3), slab, rst_hub(5), e_t, dhx, dhy, dhz, dEdR=dr_v)
      dvcp(1:3) = [dhx, dhy, dhz]
      do k = 1, n_slab
         dh(1, 3*k + 1:3*k + 3) = -dr_h(1:3, k)
         dh(2, 3*k + 1:3*k + 3) = -dr_h(1:3, k)
         du(3*k + 1:3*k + 3) = -sig_u*dr_u(1:3, k)
         dj(3*k + 1:3*k + 3) = -sig_j*dr_j(1:3, k)
         dvcp(3*k + 1:3*k + 3) = -dr_v(1:3, k)
      end do
      deallocate (dr_h, dr_u, dr_j, dr_v)
   end subroutine cau_rst_hubbard

   !------------------------------------------------------------------
   ! cau_rst_hubbard_val(q, h, u, j, vcp) - values-only delivery
   !------------------------------------------------------------------
   subroutine cau_rst_hubbard_val(q, h, u, j, vcp)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: h(2), u, j, vcp
      real(8) :: h_raw, u_raw, j_raw
      call sync_slab(q)
      call calc_rst_energy(q(1:3), slab, rst_hub(1), h_raw)
      call calc_rst_energy(q(1:3), slab, rst_hub(2), u_raw)
      call calc_rst_energy(q(1:3), slab, rst_hub(3), j_raw)
      call calc_rst_energy(q(1:3), slab, rst_hub(5), vcp)
      u = softplus(u_raw)
      j = softplus(j_raw)
      j = min(max(j, 0.0d0), 0.95d0*u)
      h(1) = h_raw
      h(2) = h_raw
   end subroutine cau_rst_hubbard_val

   ! stop_pot(who, msg) - the named-abort channel (stderr + STOP 1)
   subroutine stop_pot(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_rst_pot error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container physics layer is not usable)'
      stop 1
   end subroutine stop_pot
end module cau_rst_pot
