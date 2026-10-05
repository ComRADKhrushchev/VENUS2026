!=====================================================================
! cau_emtnn_pot.f90 - C@Au(111) EMT-NN ground PES + Hubbard delivery
!   (pure physics layer; q in, eV out, zero program knowledge)
! Design:
!   Wraps the vendored emt_nn_pes_venus engine (EMT + pair-NN potential +
!   Hubbard NN for h/U/J/Vcp, transplanted from the C@Au_TDHF/2D fit
!   kernel). cau_emtnn_load runs the engine one-time initialization
!   (Slab.xyz + the weights file; the bare-slab baseline freezes inside
!   the engine). cau_emtnn_ground delivers E0[eV] (absolute - see the
!   engine's DEVIATION note) and kernel g0 = dE0/dq[eV/A] on ALL atoms
!   (movable slab; the engine returns forces -dE/dq; criterion (c)
!   verifies the sign). cau_emtnn_hubbard / _val deliver h (both spatial
!   orbitals equal in this series), U, J, Vcp[eV] and kernels [eV/A]
!   through the engine's NN chains. PBC inside the engine with the slab
!   cell from the slab file.
!=====================================================================
module cau_emtnn_pot
   use emt_nn_pes_venus
   implicit none
   private
   public :: cau_emtnn_load, cau_emtnn_slab_atoms
   public :: cau_emtnn_ground, cau_emtnn_hubbard, cau_emtnn_hubbard_val

   logical :: loaded = .false.              ! engine-init flag
   integer :: n_slab = 0                    ! engine slab atom count

contains

   !------------------------------------------------------------------
   ! cau_emtnn_load(slab_file, weights_file, a_lat, n_side) - engine
   !                 one-time initialization (frozen bare-slab baseline)
   !------------------------------------------------------------------
   subroutine cau_emtnn_load(slab_file, weights_file, a_lat, n_side)
      character(len=*), intent(in) :: slab_file, weights_file
      real(8), intent(in) :: a_lat
      integer, intent(in) :: n_side
      integer :: u, ios
      call emt_nn_init(slab_file, weights_file, a_lat, n_side)
      open (newunit=u, file=slab_file, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call stop_pot('cau_emtnn_load', 'cannot reopen the slab file "'//trim(slab_file)//'"')
      end if
      read (u, *, iostat=ios) n_slab
      close (u)
      if (ios /= 0 .or. n_slab < 1) then
         call stop_pot('cau_emtnn_load', 'the slab file "'//trim(slab_file)// &
                       '" has no atom count on line 1')
      end if
      loaded = .true.
   end subroutine cau_emtnn_load

   integer function cau_emtnn_slab_atoms()
      if (.not. loaded) then
         call stop_pot('cau_emtnn_slab_atoms', 'the engine is not initialized (seam order error)')
      end if
      cau_emtnn_slab_atoms = n_slab
   end function cau_emtnn_slab_atoms

   !------------------------------------------------------------------
   ! cau_emtnn_ground(q, e0, g0) - ground PES: E0 and kernel dE0/dq
   !------------------------------------------------------------------
   subroutine cau_emtnn_ground(q, e0, g0)
      real(8), intent(in)  :: q(:)          ! coordinates [A], 3*(1+n_slab)
      real(8), intent(out) :: e0            ! ground potential energy [eV]
      real(8), intent(out) :: g0(:)         ! kernel dE0/dq [eV/A]
      integer :: natoms, i
      real(8), allocatable :: q2(:,:), f2(:,:)
      natoms = size(q)/3
      if (.not. loaded) then
         call stop_pot('cau_emtnn_ground', 'the engine is not initialized (seam order error)')
      end if
      if (natoms /= 1 + n_slab) then
         call stop_pot('cau_emtnn_ground', 'the geometry size does not match the engine slab')
      end if
      allocate (q2(3, natoms), f2(3, natoms))
      do i = 1, natoms
         q2(1:3, i) = q(3*i - 2:3*i)
      end do
      call calc_emt_nn_energy(natoms, 1, q2, e0, f2)
      ! E0 semantics: the engine's interaction energy E_int = E_total minus
      ! the frozen bare-slab baseline - a CONSTANT offset (drift-invariant;
      ! the same reporting the old interface used)
      do i = 1, natoms
         g0(3*i - 2:3*i) = -f2(1:3, i)
      end do
      deallocate (q2, f2)
   end subroutine cau_emtnn_ground

   !------------------------------------------------------------------
   ! pack_geom(q, q_imp, q_slab) - impurity / slab split
   !------------------------------------------------------------------
   subroutine pack_geom(q, q_imp, q_slab)
      real(8), intent(in) :: q(:)
      real(8), intent(out) :: q_imp(3)
      real(8), intent(out), allocatable :: q_slab(:,:)
      integer :: n_sl, i_slab, i3
      n_sl = size(q)/3 - 1
      allocate (q_slab(3, n_sl))
      q_imp = q(1:3)
      do i_slab = 1, n_sl
         i3 = 3*i_slab
         q_slab(1:3, i_slab) = q(i3 + 1:i3 + 3)
      end do
   end subroutine pack_geom

   !------------------------------------------------------------------
   ! cau_emtnn_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp) - Hubbard
   !                 delivery with kernels [eV/A] (engine NN chains)
   !------------------------------------------------------------------
   subroutine cau_emtnn_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: h(2), u, j, vcp
      real(8), intent(out) :: dh(2, size(q)), du(size(q)), dj(size(q)), dvcp(size(q))
      integer :: ncoord
      real(8) :: q_imp(3)
      real(8), allocatable :: q_slab(:,:)
      if (.not. loaded) then
         call stop_pot('cau_emtnn_hubbard', 'the engine is not initialized (seam order error)')
      end if
      call pack_geom(q, q_imp, q_slab)
      ncoord = size(q)
      call calc_hubbard(q_imp, q_slab, n_slab, h(1), u, j, vcp)
      h(2) = h(1)
      du = 0.0d0
      dj = 0.0d0
      dvcp = 0.0d0
      call fill_hubbard_nn_grad(du, ncoord, 2, q_imp, q_slab, n_slab, 1)
      call fill_hubbard_nn_grad(dj, ncoord, 3, q_imp, q_slab, n_slab, 1)
      call fill_hubbard_nn_grad(dvcp, ncoord, 4, q_imp, q_slab, n_slab, 1)
      call fill_hubbard_nn_grad(dh(1, :), ncoord, 1, q_imp, q_slab, n_slab, 1)
      dh(2, :) = dh(1, :)
      deallocate (q_slab)
   end subroutine cau_emtnn_hubbard

   !------------------------------------------------------------------
   ! cau_emtnn_hubbard_val(q, h, u, j, vcp) - values-only delivery
   !------------------------------------------------------------------
   subroutine cau_emtnn_hubbard_val(q, h, u, j, vcp)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: h(2), u, j, vcp
      real(8) :: q_imp(3)
      real(8), allocatable :: q_slab(:,:)
      if (.not. loaded) then
         call stop_pot('cau_emtnn_hubbard_val', 'the engine is not initialized (seam order error)')
      end if
      call pack_geom(q, q_imp, q_slab)
      call calc_hubbard(q_imp, q_slab, n_slab, h(1), u, j, vcp)
      h(2) = h(1)
      deallocate (q_slab)
   end subroutine cau_emtnn_hubbard_val

   ! stop_pot(who, msg) - the named-abort channel (stderr + STOP 1)
   subroutine stop_pot(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_emtnn_pot error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container physics layer is not usable)'
      stop 1
   end subroutine stop_pot
end module cau_emtnn_pot
