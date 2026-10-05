!=====================================================================
! hook_h3_bkmp2.f90 - PES-probe hook: the H3 BKMP2 container's potential
!                     face (systems/inactive/h3_bkmp2/bkmp2_pot.f90,
!                     bkmp2_vg), the same member code the production
!                     chain reaches through bkmp2_box/bkmp2_pes; Å/eV in,
!                     eV/Å gradient out (the module owns the bohr/
!                     hartree conversions)
! Compile (from scripts/):
!   gfortran -O2 -o probe_pes.exe ../systems/inactive/h3_bkmp2/bkmp2_pot.f90 hook_h3_bkmp2.f90 probe_pes.f90
!=====================================================================
module pes_hook
   use bkmp2_pot, only: bkmp2_vg
   implicit none
   private
   public :: pes_init, pes_entry, pes_desc
   character(len=*), parameter :: pes_desc = 'H3 BKMP2 (bkmp2_vg)'
contains
   subroutine pes_init()
   end subroutine pes_init

   subroutine pes_entry(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      if (size(mass) /= 3) then
         write (0, '(a)') 'hook_h3_bkmp2: fatal: the H3 face needs exactly 3 atoms'
         stop 1
      end if
      call bkmp2_vg(3, q, v, g)
   end subroutine pes_entry
end module pes_hook
