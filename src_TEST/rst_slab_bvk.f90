module rst_slab_bvk
!*******************************************************************************
! RST + BVK surface model for TEST_PES=RST (C / Au(111)).
!
! Ported from the 2D development repo's interface_RST.f90 — the slab_model
! module (BVK parameter read), PAIRATOMS (nearest-neighbour pair search),
! BuildTensor (Born-von Karman 3x3 force-constant tensors) and the slab-
! elastic halves of POT0/DPESHON — decoupled from venus_data/venus_params so
! it also compiles into a standalone smoke driver (scripts/smoke_rst_bvk.f90).
!
! Physics (identical to the 2D interface):
!   E_total[eV] = E_BVK(slab displacements vs eq_slab, PBC min-image)
!               + E_RST_GS(adatom vs the *equilibrium* rst_pes slab coords)
!               + EHF_baseline
!   Forces:     the adatom gets -grad E_RST_GS; slab atoms get BVK pair
!               elastic forces only. The 2D interface likewise applies no
!               RST force to slab atoms and keeps the RST pair geometry on
!               the equilibrium slab (rst_pes works on its own slab%coords,
!               which are never synced from the live coordinates).
!   BVK: nearest-neighbour force-constant tensors built once at setup from
!   the equilibrium slab (a_lat from the '# BVK:' line of the BVK file);
!   harmonic pair energy 0.5*u_ij*T_ij*u_ij, each pair counted once.
!
! Units: eV / angstrom throughout (interface_TEST converts to code units).
! Setup entry points are called once from init_test_potentials; per-step
! work goes through rst_total_vg.
!*******************************************************************************
   use rst_pes, only: rst_slab_data, rst_state, generate_fcc_slab, &
                      load_slab_xyz, load_rst_weights_standalone, &
                      calc_rst_energy, calc_rst_energy_grad
   implicit none
   private

   ! --- BVK parameters (Analytic_Potential.txt '# BVK:' line) ---------------
   real(8), save :: bvk_cutoff = 0.0d0    ! recorded, unused (2D convention)
   real(8), save :: bvk_alat   = 0.0d0    ! Au NN distance of the fitted slab
   real(8), save :: bvk_alpha  = 0.0d0
   real(8), save :: bvk_beta   = 0.0d0
   real(8), save :: bvk_gamma  = 0.0d0

   ! --- BVK pair bookkeeping (built once at setup) ---------------------------
   real(8),    save, allocatable :: eq_slab(:,:)      ! (nslab,3) equilibrium
   integer,    save, allocatable :: slab_pairs(:,:)   ! (nslab,12), 0 = empty
   real(8),    save, allocatable :: bvktensor(:,:,:,:) ! (nslab,12,3,3)

   ! --- RST pairwise NN weights + static slab container ----------------------
   type(rst_state),    save :: rst_diab(7)   ! 1-6 diabatic, 7 = GS (adiabatic)
   type(rst_slab_data), save :: rst_slab
   real(8), save :: nn_baseline = 0.0d0
   logical, save :: rst_loaded(7) = .false.

   ! --- Simulation box of the slab (min-image for the BVK sums) --------------
   real(8), save :: bvk_boxlx = 0.0d0, bvk_boxly = 0.0d0   ! angstrom
   real(8), save :: bvk_skew  = 0.0d0                       ! rad (pi/2 here)

   public :: rst_bvk_read_params, rst_bvk_load_weights
   public :: rst_bvk_generate_slab, rst_bvk_load_slab_file
   public :: rst_bvk_nslab, rst_bvk_box, rst_bvk_alat_value
   public :: rst_bvk_slab_coords, rst_bvk_energy_forces
   public :: rst_total_vg, rst_gs_energy, rst_gs_grad, rst_gs_available

contains

   !--- Read BVK parameters: '# BVK: cutoff a_lat alpha beta gamma' ------------
   ! 2D READSLAB convention, minus the silent default fallback: a missing
   ! '# BVK:' line stops loudly (the asset ships with the header).
   subroutine rst_bvk_read_params(filename)
      character(len=*), intent(in) :: filename
      character(len=512) :: line
      integer :: u, ios, pos
      logical :: found

      found = .false.
      open(newunit=u, file=trim(filename), status='old', iostat=ios)
      if (ios /= 0) then
         write(6,*) 'ERROR: cannot open RST BVK parameter file: ', trim(filename)
         stop
      end if
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         pos = index(line, 'BVK:')
         if (pos > 0) then
            read(line(pos+4:), *, iostat=ios) &
               bvk_cutoff, bvk_alat, bvk_alpha, bvk_beta, bvk_gamma
            if (ios == 0) then
               found = .true.
            else
               write(6,*) 'ERROR: bad "# BVK:" line in ', trim(filename)
               close(u); stop
            end if
            exit
         end if
      end do
      close(u)
      if (.not. found) then
         write(6,*) 'ERROR: no "# BVK: cutoff a_lat alpha beta gamma" line in ', &
                    trim(filename)
         stop
      end if

      write(6,'(A)') '  RST BVK slab parameters:'
      write(6,'(A,F10.4)') '    cutoff    = ', bvk_cutoff
      write(6,'(A,F10.4)') '    a_lat     = ', bvk_alat
      write(6,'(A,F10.4)') '    alpha_bvK = ', bvk_alpha
      write(6,'(A,F10.4)') '    beta_bvK  = ', bvk_beta
      write(6,'(A,F10.4)') '    gamma4    = ', bvk_gamma
   end subroutine rst_bvk_read_params

   !--- Load the 7-state RST pairwise NN weights ------------------------------
   subroutine rst_bvk_load_weights(filename)
      character(len=*), intent(in) :: filename
      call load_rst_weights_standalone(trim(filename), rst_diab, rst_loaded, &
                                       nn_baseline)
      write(6,'(A,ES12.4)') '  RST nn_baseline (EHF) = ', nn_baseline
   end subroutine rst_bvk_load_weights

   !--- Slab sources -----------------------------------------------------------
   ! GENERATE (default): ideal FCC(111) slab at the fitted a_lat, top layer
   ! at z=0 with a Top-site atom at the origin — geometrically self-consistent
   ! with the RST weights and the BVK a_lat, no extra asset file needed.
   ! FILE: external XYZ slab (e.g. the EMT Slab.xyz), auto box estimate.
   subroutine rst_bvk_generate_slab(a, n_layers, n_side)
      real(8), intent(in) :: a
      integer, intent(in) :: n_layers, n_side
      call generate_fcc_slab(a, n_layers, n_side, rst_slab)
      call rst_bvk_build_pairs()
   end subroutine rst_bvk_generate_slab

   subroutine rst_bvk_load_slab_file(filename)
      character(len=*), intent(in) :: filename
      call load_slab_xyz(trim(filename), rst_slab)
      call rst_bvk_build_pairs()
   end subroutine rst_bvk_load_slab_file

   !--- Nearest-neighbour pairs + BVK tensors on the equilibrium slab ----------
   ! 2D PAIRATOMS: pairs i<j with min-image distance <= a_lat + 0.1 A, at most
   ! 12 per atom; eq_slab is the equilibrium geometry (zero-energy reference).
   subroutine rst_bvk_build_pairs()
      integer :: i, j, icnt, nsl
      real(8) :: rij(3), dist2, tensor(3,3)

      nsl = rst_slab%n_atoms
      if (allocated(eq_slab)) deallocate(eq_slab)
      if (allocated(slab_pairs)) deallocate(slab_pairs)
      if (allocated(bvktensor)) deallocate(bvktensor)
      allocate(eq_slab(nsl, 3), slab_pairs(nsl, 12), bvktensor(nsl, 12, 3, 3))
      slab_pairs = 0

      bvk_boxlx = rst_slab%box_lx
      bvk_boxly = rst_slab%box_ly
      bvk_skew  = acos(-1.0d0)/2.0d0        ! rectangular slab box (2D: 90 deg)

      ! Fill the full equilibrium geometry first: the pair search below
      ! reads eq_slab(j,:) for j > i, which must not be uninitialised.
      do i = 1, nsl
         eq_slab(i,:) = rst_slab%coords(:,i)
      end do

      do i = 1, nsl
         icnt = 0
         do j = i+1, nsl
            call dist_pbc_local(eq_slab(i,:), eq_slab(j,:), rij)
            dist2 = dot_product(rij, rij)
            if (dist2 > (bvk_alat + 0.1d0)**2 .or. icnt >= 12) cycle
            icnt = icnt + 1
            slab_pairs(i, icnt) = j
            call build_tensor_local(i, j, tensor)
            bvktensor(i, icnt, :, :) = tensor
         end do
      end do
      write(6,'(A,I0,A)') '  BVK pair entries found: ', count(slab_pairs /= 0)
      write(6,'(A,I0,A,F8.3,A,F8.3,A)') '  BVK pairs built on ', nsl, &
         ' slab atoms; box ', bvk_boxlx, ' x ', bvk_boxly, ' A (skew 90 deg)'
   end subroutine rst_bvk_build_pairs

   !--- BVK 3x3 force-constant tensor for one pair (2D BuildTensor) ------------
   ! Rotates the 6 fcc basis tensors into the slab frame; the pair direction
   ! (min-image, in slab frame) selects the matching basis (|cos| ~ 1).
   subroutine build_tensor_local(index1, index2, tensor)
      integer, intent(in)  :: index1, index2
      real(8), intent(out) :: tensor(3,3)
      real(8) :: rij(3), r_dir(3), a_base(3,3,6), fcc2slab(3,3)
      real(8) :: base_dir(3,6), asqr2, asqr3, asqr6, paral, d
      integer :: i

      asqr2 = 1.0d0/dsqrt(2.0d0)
      asqr3 = 1.0d0/dsqrt(3.0d0)
      asqr6 = 1.0d0/dsqrt(6.0d0)
      tensor = 0.0d0
      fcc2slab = reshape([ &
         asqr2,  asqr2,  0.0d0, &
         asqr6, -asqr6, -2.0d0*asqr6, &
         asqr3, -asqr3,  asqr3], [3,3], order=[2,1])

      call dist_pbc_local(eq_slab(index1,:), eq_slab(index2,:), rij)
      d = sqrt(dot_product(rij, rij))
      r_dir = rij/d
      r_dir = matmul(transpose(fcc2slab), r_dir)

      base_dir(:,1) = [0.0d0, asqr2, asqr2]
      a_base(:,:,1) = reshape([ &
         bvk_alpha, 0.0d0,     0.0d0,    &
         0.0d0,     bvk_beta,  bvk_gamma,&
         0.0d0,     bvk_gamma, bvk_beta ], [3,3])

      base_dir(:,2) = [0.0d0, asqr2, -asqr2]
      a_base(:,:,2) = reshape([ &
         bvk_alpha, 0.0d0,      0.0d0,    &
         0.0d0,     bvk_beta,  -bvk_gamma,&
         0.0d0,    -bvk_gamma,  bvk_beta ], [3,3])

      base_dir(:,3) = [asqr2, 0.0d0, asqr2]
      a_base(:,:,3) = reshape([ &
         bvk_beta,  0.0d0,     bvk_gamma,&
         0.0d0,     bvk_alpha, 0.0d0,    &
         bvk_gamma, 0.0d0,     bvk_beta ], [3,3])

      base_dir(:,4) = [-asqr2, 0.0d0, asqr2]
      a_base(:,:,4) = reshape([ &
         bvk_beta,   0.0d0,     -bvk_gamma,&
         0.0d0,      bvk_alpha,  0.0d0,   &
         -bvk_gamma, 0.0d0,      bvk_beta ], [3,3])

      base_dir(:,5) = [asqr2, asqr2, 0.0d0]
      a_base(:,:,5) = reshape([ &
         bvk_beta,  bvk_gamma, 0.0d0,   &
         bvk_gamma, bvk_beta,  0.0d0,   &
         0.0d0,     0.0d0,     bvk_alpha], [3,3])

      base_dir(:,6) = [asqr2, -asqr2, 0.0d0]
      a_base(:,:,6) = reshape([ &
         bvk_beta,  -bvk_gamma, 0.0d0,   &
         -bvk_gamma, bvk_beta,  0.0d0,   &
         0.0d0,      0.0d0,     bvk_alpha], [3,3])

      do i = 1, 6
         paral = dot_product(base_dir(:,i), r_dir)
         if (abs(abs(paral) - 1.0d0) < 1.0d-3) then
            tensor = matmul(fcc2slab, matmul(a_base(:,:,i), transpose(fcc2slab)))
            exit
         end if
      end do
   end subroutine build_tensor_local

   !--- Min-image distance vector, skewed box, 9 images (2D Dist_pbc) ----------
   subroutine dist_pbc_local(qa, qb, rvec)
      real(8), intent(in)  :: qa(3), qb(3)
      real(8), intent(out) :: rvec(3)
      real(8) :: lat_x(3), lat_y(3), rr(3), d2(9), rimgs(3,9)
      integer :: k, kmin

      lat_x = (/ bvk_boxlx, 0.0d0, 0.0d0 /)
      lat_y = (/ bvk_boxly*cos(bvk_skew), bvk_boxly*sin(bvk_skew), 0.0d0 /)
      do k = 1, 3
         rimgs(k,1) = qb(k)
         rimgs(k,2) = qb(k) + lat_x(k)
         rimgs(k,3) = qb(k) - lat_x(k)
         rimgs(k,4) = qb(k) + lat_y(k)
         rimgs(k,5) = qb(k) - lat_y(k)
         rimgs(k,6) = qb(k) + lat_x(k) + lat_y(k)
         rimgs(k,7) = qb(k) + lat_x(k) - lat_y(k)
         rimgs(k,8) = qb(k) - lat_x(k) + lat_y(k)
         rimgs(k,9) = qb(k) - lat_x(k) - lat_y(k)
      end do
      do k = 1, 9
         rr = qa - rimgs(:,k)
         d2(k) = dot_product(rr, rr)
      end do
      kmin = minloc(d2, dim=1)
      rvec = qa - rimgs(:,kmin)
   end subroutine dist_pbc_local

   !--- Small accessors ---------------------------------------------------------
   integer function rst_bvk_nslab()
      rst_bvk_nslab = rst_slab%n_atoms
   end function rst_bvk_nslab

   subroutine rst_bvk_box(lx, ly)
      real(8), intent(out) :: lx, ly
      lx = bvk_boxlx
      ly = bvk_boxly
   end subroutine rst_bvk_box

   real(8) function rst_bvk_alat_value()
      rst_bvk_alat_value = bvk_alat
   end function rst_bvk_alat_value

   !--- Slab equilibrium coords, layout (3, nslab) -> QZB -----------------------
   subroutine rst_bvk_slab_coords(coords_out)
      real(8), intent(out) :: coords_out(3, rst_slab%n_atoms)
      coords_out = rst_slab%coords
   end subroutine rst_bvk_slab_coords

   logical function rst_gs_available()
      rst_gs_available = rst_loaded(7)
   end function rst_gs_available

   !--- BVK elastic energy + slab forces for live slab positions ----------------
   ! E_pair = 0.5*u*T*u with u = d_ij(live) - d_ij(eq), min-image; each pair
   ! counted once (2D POT0 / DPESHON pair loop, forces = -dE/dR).
   subroutine rst_bvk_energy_forces(coords, nsl, e_bvk, forces)
      integer, intent(in)  :: nsl
      real(8), intent(in)  :: coords(3, nsl)
      real(8), intent(out) :: e_bvk, forces(3, nsl)
      real(8) :: r_vec(3), r_eq(3), r_pt(3), tensor(3,3), f_vec(3)
      integer :: j, k, j2

      e_bvk  = 0.0d0
      forces = 0.0d0
      do j = 1, nsl
         do k = 1, 12
            j2 = slab_pairs(j, k)
            if (j2 == 0) cycle
            call dist_pbc_local(coords(:,j),  coords(:,j2),  r_vec)
            call dist_pbc_local(eq_slab(j,:), eq_slab(j2,:), r_eq)
            r_pt   = r_vec - r_eq
            tensor = bvktensor(j, k, :, :)
            e_bvk  = e_bvk + 0.5d0*dot_product(matmul(tensor, r_pt), r_pt)
            f_vec  = matmul(r_pt, tensor)      ! T symmetric: = T*u = -dE/dr_j
            forces(:,j)  = forces(:,j)  - f_vec
            forces(:,j2) = forces(:,j2) + f_vec
         end do
      end do
   end subroutine rst_bvk_energy_forces

   !--- RST GS adsorbate-slab energy (eV), equilibrium slab geometry -----------
   subroutine rst_gs_energy(xyz, e_nn)
      real(8), intent(in)  :: xyz(3)
      real(8), intent(out) :: e_nn
      if (.not. rst_loaded(7)) then
         write(6,*) 'ERROR: RST GS state not loaded (missing in weight file)'
         stop
      end if
      call calc_rst_energy(xyz, rst_slab, rst_diab(7), e_nn)
      e_nn = e_nn + nn_baseline
   end subroutine rst_gs_energy

   !--- RST GS gradient (eV/A), equilibrium slab geometry -----------------------
   subroutine rst_gs_grad(xyz, e_nn, g3)
      real(8), intent(in)  :: xyz(3)
      real(8), intent(out) :: e_nn, g3(3)
      if (.not. rst_loaded(7)) then
         write(6,*) 'ERROR: RST GS state not loaded (missing in weight file)'
         stop
      end if
      call calc_rst_energy_grad(xyz, rst_slab, rst_diab(7), e_nn, g3(1), g3(2), g3(3))
      e_nn = e_nn + nn_baseline
   end subroutine rst_gs_grad

   !--- Total PES for TEST_PES=RST ---------------------------------------------
   ! q(1:3) = C adatom, q(4:3*natom) = slab (live). Returns V[eV] and
   ! g = +dV/dq [eV/A] — the interface_TEST contract (g = -forces).
   subroutine rst_total_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      real(8) :: coords(3, natom-1), forces(3, natom-1), e_bvk, e_nn, g3(3)
      integer :: i, nsl

      nsl = natom - 1
      v = 0.0d0
      g = 0.0d0

      ! Slab elastic (live slab positions)
      do i = 1, nsl
         coords(1,i) = q(3*(i+1)-2)
         coords(2,i) = q(3*(i+1)-1)
         coords(3,i) = q(3*(i+1))
      end do
      call rst_bvk_energy_forces(coords, nsl, e_bvk, forces)
      v = v + e_bvk
      do i = 1, nsl
         g(3*(i+1)-2) = -forces(1,i)
         g(3*(i+1)-1) = -forces(2,i)
         g(3*(i+1))   = -forces(3,i)
      end do

      ! RST GS adsorbate term (equilibrium slab geometry, 2D convention)
      if (rst_loaded(7)) then
         call rst_gs_grad(q(1:3), e_nn, g3)
         v = v + e_nn
         g(1:3) = g3
      end if
   end subroutine rst_total_vg

end module rst_slab_bvk
