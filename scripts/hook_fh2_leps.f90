!=====================================================================
! hook_fh2_leps.f90 - PES-probe hook: the F+H2 modified-LEPS container's
!                     potential face (systems/inactive/fh2_leps/
!                     fh2_pot.f90, fh2_vg), the same member code the
!                     production chain reaches through fh2_box/fh2_pes;
!                     parameters at the legislated literature defaults
!                     (the T-S set - fh2_set_defaults stages exactly the
!                     module's initial values; the call is kept explicit
!                     as self-documentation)
! Compile (from scripts/):
!   gfortran -O2 -o probe_pes.exe ../systems/inactive/fh2_leps/fh2_pot.f90 hook_fh2_leps.f90 probe_pes.f90
!=====================================================================
module pes_hook
   use fh2_pot, only: fh2_vg, fh2_set_defaults
   implicit none
   private
   public :: pes_init, pes_entry, pes_desc
   character(len=*), parameter :: pes_desc = 'F+H2 modified LEPS (fh2_vg, default literature set)'
contains
   subroutine pes_init()
      call fh2_set_defaults()
   end subroutine pes_init

   subroutine pes_entry(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      if (size(mass) /= 3) then
         write (0, '(a)') 'hook_fh2_leps: fatal: the F+H2 face needs exactly 3 atoms'
         stop 1
      end if
      call fh2_vg(3, q, v, g)
   end subroutine pes_entry
end module pes_hook
