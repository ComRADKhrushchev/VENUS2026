!=====================================================================
! hook_reference.f90 - PES-probe hook: the reference container's well
!                      face (systems/reference/reference_pot.f90,
!                      harm_vg), the same member code the production
!                      chain reaches through reference_box/reference_pes
! Note: the constant is buffered at the box default K = 0.04 eV/A^2
!   (reference_box k_default). The probe binds the potential face, so
!   the box's run-directory parameter file (reference_pars.txt) is NOT
!   read - stage another constant in pes_init below to scan a
!   different well.
! Compile (from scripts/):
!   gfortran -O2 -o probe_pes.exe ../systems/reference/reference_pot.f90 hook_reference.f90 probe_pes.f90
!=====================================================================
module pes_hook
   use reference_pot, only: harm_vg, harm_set_constant
   implicit none
   private
   public :: pes_init, pes_entry, pes_desc
   character(len=*), parameter :: pes_desc = 'reference harmonic well (harm_vg, K = 0.04 eV/A^2)'
contains
   subroutine pes_init()
      call harm_set_constant(0.04d0)   ! the box default (reference_box k_default)
   end subroutine pes_init

   subroutine pes_entry(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      call harm_vg(size(mass), q, v, g)
   end subroutine pes_entry
end module pes_hook
