!=====================================================================
! hook_h2ag111.f90 - PES-probe hook: the H2/Ag(111) NN container's
!                    potential face (systems/inactive/h2ag111PES/
!                    h2ag111_pot.f90 - the legacy pot0 energy +
!                    dpeshon gradient pair), the same member code the
!                    production chain reaches through h2ag111_box/
!                    h2ag111_pes. The molecule is the two H atoms in
!                    q; the slab is implicit in the network.
! Validity (the pot's own header): the network interpolates its
!   training region - keep scans inside the fitted domain (an intact
!   parallel molecule is usable to Z ~ 8 A; stretched or tilted far
!   geometries extrapolate garbage)
! Run material: pes_init reads weights-h2ag.txt + biases-h2ag.txt from
!   the CURRENT directory - copy them from the container folder beside
!   the run
! Compile (from scripts/):
!   gfortran -O2 -o probe_pes.exe ../systems/inactive/h2ag111PES/h2ag111_pot.f90 hook_h2ag111.f90 probe_pes.f90
!=====================================================================
module pes_hook
   implicit none
   private
   public :: pes_init, pes_entry, pes_desc
   character(len=*), parameter :: pes_desc = 'H2/Ag(111) NN PES (pot0 + dpeshon)'

   ! the legacy init (reads the weight/bias files from cwd) is exposed
   ! as the hook's pes_init directly - it is already the no-argument
   ! setup the contract asks for
   interface
      subroutine pot0(natoms, q, v)          ! the legacy face: energy [eV]
         real(8) :: q(*), v
         integer :: natoms
      end subroutine pot0
      subroutine dpeshon(natoms, q, pdot)    ! the legacy face: dV/dq [eV/A]
         real(8) :: q(*), pdot(*)
         integer :: natoms
      end subroutine dpeshon
      subroutine pes_init()
      end subroutine pes_init
   end interface
contains

   subroutine pes_entry(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      if (size(mass) /= 2) then
         write (0, '(a)') 'hook_h2ag111: fatal: the H2/Ag(111) face needs exactly 2 atoms'
         stop 1
      end if
      call pot0(2, q, v)
      call dpeshon(2, q, g)
   end subroutine pes_entry
end module pes_hook
