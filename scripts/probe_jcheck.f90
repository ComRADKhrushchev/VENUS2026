!=====================================================================
! probe_jcheck.f90 - J(z) sign/structure and seed-ordering audit for cau_peri
!   (both the analytic branch and the NN branch)
! Design:
!   For z along the on-top normal, prints J(z), U(z) from the hubbard chain,
!   and the composed potential V(z) seen by EACH diabatic seed (3P = 'GS'
!   parallel-spin seed, 1D = antiparallel seed; their difference must equal
!   J(z) exactly - a one-line algebra identity at diagonal seeds). Run twice:
!   'ana' (analytic branch) and 'nn' (nn_weights.txt branch).
! Compile (scripts/):
!   gfortran -O2 -o probe_jcheck.exe ../src/globals/consts.f90
!     ../src/globals/state.f90 ../src/utils/linalg.f90 ../src/utils/rng.f90
!     ../systems/inactive/cau_peri/surface_potentials.f90
!     ../systems/inactive/cau_peri/cau_peri_pot.f90
!     ../systems/inactive/cau_peri/cau_peri_tdhf.f90
!     ../systems/inactive/cau_peri/cau_peri_box.f90 probe_jcheck.f90
!=====================================================================
program probe_jcheck
   use cau_peri_pot,  only: cau_peri_load
   use cau_peri_tdhf, only: cau_peri_tdhf_setup, cau_peri_pes, hub_vals_out
   implicit none
   real(8) :: q(9), mass(3), g(9), v3p, v1d, hv(2), u, j, vcp
   character(len=8) :: mode
   integer :: iz, k

   call get_command_argument(1, mode)
   if (len_trim(mode) == 0) mode = 'ana'
   if (trim(mode) == 'nn') then
      call cau_peri_load('nn_weights.txt')
   else
      call cau_peri_load()
   end if
   mass = [196.97d0, 196.97d0, 12.01d0]
   print '(a)', '# z  J(z)  U(z)  V_3Pseed(z)  V_1Dseed(z)  V1D-V3P  (eV)'
   do iz = 0, 70
      q = 0.0d0
      q(6) = 0.9d0
      q(9) = 1.0d0 + 0.1d0*iz
      call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 1, 0, 0.05d0)
      call cau_peri_pes(q, mass, g, v3p)
      call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, '1D', 1, 0, 0.05d0)
      call cau_peri_pes(q, mass, g, v1d)
      call hub_vals_out(hv, u, j, vcp)
      write (*, '(f6.2,5f12.4)') q(9), j, u, v3p, v1d, v1d - v3p
   end do
end program probe_jcheck
