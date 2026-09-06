!******************************************************************************
! smoke_rst_bvk — standalone acceptance driver for TEST_PES=RST.
!
! Links only rst_pes + rst_slab_bvk (no VENUS core). Run from the repo root:
!   ifx -r8 -double-size=64 -i8 -O2 -w \
!     src_TEST/rst_pes.f90 src_TEST/rst_slab_bvk.f90 scripts/smoke_rst_bvk.f90 \
!     -o smoke_rst_bvk.e
!   ./smoke_rst_bvk.e data/rst/nn_weights_rst.txt data/rst/Analytic_Potential.txt
!
! Checks:
!   1) RST GS adsorbate energy at the Top site (0,0,z), z = 1.5..3.0 A:
!      finite, well depth in the -1..-3 eV range, gradient nonzero;
!      analytic gradient vs central finite differences.
!   2) BVK elastic slab term: zero at the equilibrium slab; displacing one
!      atom by 0.1 A gives E > 0 and a restoring force pointing back to the
!      equilibrium position; force vs central finite differences.
!******************************************************************************
program smoke_rst_bvk
   use rst_slab_bvk
   implicit none
   integer, parameter :: nlay = 4, nside = 6
   real(8), parameter :: alat = 2.95d0
   real(8) :: e_nn, g3(3), g3_fd(3), h, fdmax
   real(8) :: coords(3, nlay*nside*nside), forces(3, nlay*nside*nside)
   real(8) :: e_bvk0, e_bvk1, f_fd, f_ana
   real(8) :: z, e_min, z_min, e_here, e_inf
   integer :: k, nsl, iz
   character(len=512) :: wfile, bfile
   real(8) :: xyz(3)

   call get_command_argument(1, wfile)
   call get_command_argument(2, bfile)
   if (len_trim(wfile) == 0) wfile = 'data/rst/nn_weights_rst.txt'
   if (len_trim(bfile) == 0) bfile = 'data/rst/Analytic_Potential.txt'

   call rst_bvk_read_params(trim(bfile))
   call rst_bvk_load_weights(trim(wfile))
   call rst_bvk_generate_slab(alat, nlay, nside)
   nsl = rst_bvk_nslab()
   write(*,*)

   !--- 1) RST GS at the Top site: scan z, analytic vs FD gradient -------------
   ! Energy zero: value at z = 25 A (every slab atom beyond the 6 A cutoff),
   ! i.e. the asymptotic E0_GS + EHF_baseline constant; E_ads = E(z) - E_ref
   ! is the adsorption energy. (The VENUS-side potential keeps the baseline
   ! offset — dynamics only involves differences, and POTENZ's VZERO shifts
   ! the initial geometry to zero.)
   call rst_gs_energy((/0.0d0, 0.0d0, 25.0d0/), e_inf)
   write(*,'(A)') '=== smoke 1: RST GS adsorbate at Top site (0,0,z) ==='
   write(*,'(A,ES12.4)') ' asymptotic E(z=25 A) = ', e_inf
   write(*,'(A)') '   z(A)   E_ads(eV)  dE/dz(eV/A)'
   e_min = 1.0d30
   do iz = 10, 45
      z = dble(iz)/10.0d0
      call rst_gs_grad((/0.0d0, 0.0d0, z/), e_here, g3)
      write(*,'(F6.2,ES13.4,ES15.4)') z, e_here - e_inf, g3(3)
      if (e_here < e_min) then
         e_min = e_here
         z_min = z
      end if
   end do
   write(*,'(A,F5.2,A,ES11.4,A)') ' well minimum: z = ', z_min, &
      ' A, E_ads = ', e_min - e_inf, ' eV'

   xyz = (/0.0d0, 0.0d0, 2.0d0/)
   call rst_gs_grad(xyz, e_nn, g3)
   h = 1.0d-4
   fdmax = 0.0d0
   do k = 1, 3
      call rst_gs_energy(xyz + h*unit_vec(k), e_nn)
      e_bvk0 = e_nn
      call rst_gs_energy(xyz - h*unit_vec(k), e_nn)
      g3_fd(k) = (e_bvk0 - e_nn)/(2.0d0*h)
      fdmax = max(fdmax, abs(g3_fd(k) - g3(k)))
   end do
   write(*,'(A)') '--- at z = 2.0 A: analytic vs finite-difference gradient ---'
   write(*,'(A,3ES14.5)') ' analytic : ', g3
   write(*,'(A,3ES14.5)') ' FD       : ', g3_fd
   write(*,'(A,ES10.3)')  ' max diff : ', fdmax
   write(*,'(A,L1)')      ' |grad| > 0              : ', sqrt(dot_product(g3,g3)) > 0.0d0
   write(*,'(A,ES9.2)')   ' |grad| (eV/A)           : ', sqrt(dot_product(g3,g3))

   !--- 2) BVK elastic slab term ----------------------------------------------
   write(*,*)
   write(*,'(A)') '=== smoke 2: BVK elastic slab ==='
   call rst_bvk_slab_coords(coords)
   call rst_bvk_energy_forces(coords, nsl, e_bvk0, forces)
   write(*,'(A,ES12.4)') ' E_BVK at equilibrium (eV)   : ', e_bvk0
   write(*,'(A,ES12.4)') ' |F| max at equilibrium      : ', maxval(abs(forces))

   coords(1,1) = coords(1,1) + 0.1d0     ! displace atom 1 by +0.1 A in x
   call rst_bvk_energy_forces(coords, nsl, e_bvk1, forces)
   f_ana = forces(1,1)                    ! analytic Fx at u = 0.1 A
   write(*,'(A,ES12.4)') ' E_BVK after 0.1 A bump (eV) : ', e_bvk1
   write(*,'(A,3ES12.4)')' F(atom 1) (eV/A)            : ', forces(:,1)
   write(*,'(A,L1)')     ' E_BVK > 0                   : ', e_bvk1 > 0.0d0
   write(*,'(A,L1)')     ' Fx(atom 1) restoring (<0)   : ', forces(1,1) < 0.0d0

   ! FD check of the atom-1 x-force: F_x = -dE/dx about u = 0.1 A
   coords(1,1) = coords(1,1) + h
   call rst_bvk_energy_forces(coords, nsl, e_nn, forces)
   e_bvk0 = e_nn
   coords(1,1) = coords(1,1) - 2.0d0*h
   call rst_bvk_energy_forces(coords, nsl, e_nn, forces)
   f_fd = -(e_bvk0 - e_nn)/(2.0d0*h)
   write(*,'(A)') '--- atom-1 x-force: analytic vs finite difference ---'
   write(*,'(A,ES14.5)') ' analytic : ', f_ana
   write(*,'(A,ES14.5)') ' FD       : ', f_fd
   write(*,'(A,ES10.3)') ' diff     : ', abs(f_ana - f_fd)

contains

   function unit_vec(k) result(u)
      integer, intent(in) :: k
      real(8) :: u(3)
      u = 0.0d0
      u(k) = 1.0d0
   end function unit_vec

end program smoke_rst_bvk
