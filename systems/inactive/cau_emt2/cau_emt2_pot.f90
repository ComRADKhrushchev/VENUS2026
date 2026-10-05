!=====================================================================
! cau_emt2_pot.f90 - C@Au(111) EMT ground PES + Hubbard delivery (pure
!   physics layer of the container; q in, eV out, zero program knowledge)
! Design:
!   Wraps the vendored pes_emt_mod engine (2-species classical EMT +
!   Hubbard polynomial machinery, translated from the C@Au_TDHF/2D tree).
!   cau_emt2_ground delivers the ground-state energy E0[eV] and kernel
!   g0 = dE0/dq[eV/A] (the engine returns forces -dE/dq; sign verified by
!   check_cau_emt2 criterion (c)). Each call refreshes the module-cached
!   adsorbate embedding density ad_density and its gradient - the shared
!   state the Hubbard delivery reads one call later.
!   cau_emt2_hubbard delivers h, U, J, Vcp[eV] with kernels [eV/A];
!   cau_emt2_hubbard_val the values only (no density-gradient chain).
!   Both read the density cache - a cau_emt2_ground call at the same q
!   must come first (the ordering the old interface DPESHON -> H_Orb used).
!   Atom 1 = C adsorbate (NATOMA = 1, fixed), atoms 2..n = Au slab.
!=====================================================================
module cau_emt2_pot
   use pes_emt_mod
   implicit none
   private
   public :: cau_emt2_load_params, cau_emt2_set_box, cau_emt2_params_loaded
   public :: cau_emt2_box_out, cau_emt2_cuts_out
   public :: cau_emt2_ground, cau_emt2_hubbard, cau_emt2_hubbard_val

   type(emt_2spec_params) :: par            ! EMT + Hubbard polynomial parameter set
   logical :: loaded = .false.              ! parameter-load flag
   real(8) :: box_lx = 17.7d0               ! surface cell vector A length [A]
   real(8) :: box_ly = 17.7d0               ! surface cell vector B length [A]
   real(8) :: box_skew = 1.0471975511965976d0 ! cell angle [rad] (60 deg = hexagonal; the
                                              ! production Au(111) 6x6 supercell value)

contains

   !------------------------------------------------------------------
   ! cau_emt2_load_params(file) - read the EMT parameter file (14 EMT
   !                 keywords + 32 Hubbard polynomial coefficients, the
   !                 vendored read_emt_params grammar "KEY value"); a
   !                 missing/unreadable file is a named abort
   !------------------------------------------------------------------
   subroutine cau_emt2_load_params(file)
      character(len=*), intent(in) :: file  ! parameter file name (run-directory relative)
      integer :: ios
      call read_emt_params(file, par, ios)
      if (ios /= 0) then
         call stop_pot('cau_emt2_load_params', 'cannot read the EMT parameter file "'// &
                       trim(file)//'" (missing or malformed - the container is not usable)')
      end if
      loaded = .true.
   end subroutine cau_emt2_load_params

   !------------------------------------------------------------------
   ! cau_emt2_set_box(lx, ly, skew_deg) - set the periodic surface cell
   !                 (must match the slab lattice the atom list placed)
   !------------------------------------------------------------------
   subroutine cau_emt2_set_box(lx, ly, skew_deg)
      real(8), intent(in) :: lx, ly         ! cell vector lengths [A]
      real(8), intent(in) :: skew_deg       ! cell angle [deg]
      box_lx = lx
      box_ly = ly
      box_skew = skew_deg*(3.14159265358979323846d0/180.0d0)
   end subroutine cau_emt2_set_box

   logical function cau_emt2_params_loaded()
      cau_emt2_params_loaded = loaded .and. hubbard_poly_loaded
   end function cau_emt2_params_loaded

   !------------------------------------------------------------------
   ! cau_emt2_box_out(lx, ly, skew) / cau_emt2_cuts_out(rcut, acut) -
   !                 read-only hand-out of the cell / cutoff state
   !------------------------------------------------------------------
   subroutine cau_emt2_box_out(lx, ly, skew)
      real(8), intent(out) :: lx, ly, skew
      lx = box_lx
      ly = box_ly
      skew = box_skew
   end subroutine cau_emt2_box_out

   subroutine cau_emt2_cuts_out(rcut, acut)
      real(8), intent(out) :: rcut, acut
      rcut = par%rcut
      acut = par%acut
   end subroutine cau_emt2_cuts_out

   !------------------------------------------------------------------
   ! cau_emt2_ground(q, e0, g0) - ground PES: E0 and kernel dE0/dq
   !                 (also refreshes the density cache the Hubbard
   !                 delivery reads; forces sign verified by criterion (c))
   !------------------------------------------------------------------
   subroutine cau_emt2_ground(q, e0, g0)
      real(8), intent(in)  :: q(:)          ! coordinates [A], flattened 3*natoms (atom 1 = C)
      real(8), intent(out) :: e0            ! ground potential energy [eV]
      real(8), intent(out) :: g0(:)         ! kernel dE0/dq [eV/A] (size = size(q))
      integer :: natoms, i
      real(8), allocatable :: q2(:,:), f2(:,:)
      real(8) :: energy
      ! 1. guard + pack the flat vector into the engine's (3, natoms) layout
      natoms = size(q)/3
      if (.not. loaded) then
         call stop_pot('cau_emt2_ground', 'the EMT parameters are not loaded (seam order error)')
      end if
      allocate (q2(3, natoms), f2(3, natoms))
      do i = 1, natoms
         q2(1:3, i) = q(3*i - 2:3*i)
      end do
      ! 2. one engine call (caches ad_density + ad_density_grad module-side)
      call compute_emt_2species(natoms, 1, q2, par, energy, f2, &
                                box_lx=box_lx, box_ly=box_ly, box_skew=box_skew)
      ! 3. deliver the kernel: the engine returns forces -dE/dq
      do i = 1, natoms
         g0(3*i - 2:3*i) = -f2(1:3, i)
      end do
      e0 = energy
      deallocate (q2, f2)
   end subroutine cau_emt2_ground

   !------------------------------------------------------------------
   ! pack_geom(q, q_imp, q_slab) - split the flat vector into the
   !                 impurity position and the slab block
   !------------------------------------------------------------------
   subroutine pack_geom(q, q_imp, q_slab)
      real(8), intent(in) :: q(:)
      real(8), intent(out) :: q_imp(3)
      real(8), intent(out), allocatable :: q_slab(:,:)
      integer :: n_slab, i_slab, i3
      n_slab = size(q)/3 - 1
      allocate (q_slab(3, n_slab))
      q_imp = q(1:3)
      do i_slab = 1, n_slab
         i3 = 3*i_slab
         q_slab(1:3, i_slab) = q(i3 + 1:i3 + 3)
      end do
   end subroutine pack_geom

   !------------------------------------------------------------------
   ! pair_sums(q, ps) - the four polynomial pairwise sums (skipped when
   !                 the amplitude is zero - the old H_Orb skip)
   !------------------------------------------------------------------
   subroutine pair_sums(q, ps)
      real(8), intent(in) :: q(:)
      real(8), intent(out) :: ps(4)
      real(8) :: q_imp(3)
      real(8), allocatable :: q_slab(:,:)
      integer :: n_slab, k
      call pack_geom(q, q_imp, q_slab)
      n_slab = size(q_slab, 2)
      ps = 0.0d0
      do k = 1, 4
         if (abs(hubb_poly(k)%C_amp) > 1.0d-15) then
            ps(k) = compute_pairwise_sum(q_imp, q_slab, n_slab, hubb_poly(k)%D_dist, &
                                         par%rcut, par%acut, box_lx, box_ly, box_skew)
         end if
      end do
      deallocate (q_slab)
   end subroutine pair_sums

   !------------------------------------------------------------------
   ! cau_emt2_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp) - Hubbard
   !                 delivery with kernels: values [eV] (h per spatial
   !                 orbital - both equal in this series), kernels [eV/A]
   !------------------------------------------------------------------
   subroutine cau_emt2_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp)
      real(8), intent(in)  :: q(:)          ! coordinates [A]
      real(8), intent(out) :: h(2), u, j, vcp  ! Hubbard parameters [eV]
      real(8), intent(out) :: dh(2, size(q)), du(size(q)), dj(size(q)), dvcp(size(q)) ! kernels [eV/A]
      integer :: natoms, n_slab, ncoord
      real(8) :: q_imp(3)
      real(8), allocatable :: q_slab(:,:)
      real(8) :: ps(4)
      ! 1. guards + geometry pack
      natoms = size(q)/3
      n_slab = natoms - 1
      ncoord = 3*natoms
      if (.not. cau_emt2_params_loaded()) then
         call stop_pot('cau_emt2_hubbard', 'the Hubbard polynomials are not loaded (seam order error)')
      end if
      if (.not. allocated(ad_density_grad)) then
         call stop_pot('cau_emt2_hubbard', 'the density cache is empty (call cau_emt2_ground '// &
                       'at the same q first - an ordering error)')
      end if
      call pack_geom(q, q_imp, q_slab)
      ! 2. values (pairwise sums + polynomials on the cached density)
      call pair_sums(q, ps)
      h(1) = eval_hubbard_poly(ad_density, hubb_poly(1), ps(1))
      u    = eval_hubbard_poly(ad_density, hubb_poly(2), ps(2))
      j    = eval_hubbard_poly(ad_density, hubb_poly(3), ps(3))
      vcp  = eval_hubbard_poly(ad_density, hubb_poly(4), ps(4))
      h(2) = h(1)
      ! 3. kernels (density chain + pairwise chain, the vendored assembler)
      dh = 0.0d0
      du = 0.0d0
      dj = 0.0d0
      dvcp = 0.0d0
      call fill_hubbard_param_grad(dh(1, :), ncoord, 1, q_imp, q_slab, n_slab, par%rcut, par%acut, &
                                   1, box_lx, box_ly, box_skew)
      dh(2, :) = dh(1, :)
      call fill_hubbard_param_grad(du, ncoord, 2, q_imp, q_slab, n_slab, par%rcut, par%acut, &
                                   1, box_lx, box_ly, box_skew)
      call fill_hubbard_param_grad(dj, ncoord, 3, q_imp, q_slab, n_slab, par%rcut, par%acut, &
                                   1, box_lx, box_ly, box_skew)
      call fill_hubbard_param_grad(dvcp, ncoord, 4, q_imp, q_slab, n_slab, par%rcut, par%acut, &
                                   1, box_lx, box_ly, box_skew)
      deallocate (q_slab)
   end subroutine cau_emt2_hubbard

   !------------------------------------------------------------------
   ! cau_emt2_hubbard_val(q, h, u, j, vcp) - values-only Hubbard delivery
   !                 (no density-gradient chain; the prop path uses this)
   !------------------------------------------------------------------
   subroutine cau_emt2_hubbard_val(q, h, u, j, vcp)
      real(8), intent(in)  :: q(:)          ! coordinates [A]
      real(8), intent(out) :: h(2), u, j, vcp  ! Hubbard parameters [eV]
      real(8) :: ps(4)
      if (.not. cau_emt2_params_loaded()) then
         call stop_pot('cau_emt2_hubbard_val', 'the Hubbard polynomials are not loaded (seam order error)')
      end if
      call pair_sums(q, ps)
      h(1) = eval_hubbard_poly(ad_density, hubb_poly(1), ps(1))
      u    = eval_hubbard_poly(ad_density, hubb_poly(2), ps(2))
      j    = eval_hubbard_poly(ad_density, hubb_poly(3), ps(3))
      vcp  = eval_hubbard_poly(ad_density, hubb_poly(4), ps(4))
      h(2) = h(1)
   end subroutine cau_emt2_hubbard_val

   ! stop_pot(who, msg) - the named-abort channel (stderr + STOP 1)
   subroutine stop_pot(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_emt2_pot error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container physics layer is not usable)'
      stop 1
   end subroutine stop_pot
end module cau_emt2_pot
