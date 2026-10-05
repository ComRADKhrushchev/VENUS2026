!=====================================================================
! reg_inc_single.f90 - register the 'incident_single' member (guarded:
!   the molecular paradigm with a single fragment)
! Design:
!   One samp_reg call fills the row when the run derived the molecular
!   paradigm with one fragment - the paradigm-guarded registration keeps
!   exactly ONE incident member armed per process. The member owns no
!   parameters, so this reg carries NO _init entry (the optional half of
!   the naming contract): there is nothing to load.
!=====================================================================
module reg_inc_single
   use sampler,       only: samp_reg
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use inc_single,    only: inc_single_draw, inc_single_realize
   implicit none
   private
   public :: reg_inc_single_all
contains
   !------------------------------------------------------------------
   ! reg_inc_single_all() - fill the member-table row (single-fragment
   !                        molecular paradigm only)
   !------------------------------------------------------------------
   subroutine reg_inc_single_all()
      if (.not. armed_here()) return
      call samp_reg('incident_single', inc_single_draw, inc_single_realize)
   end subroutine reg_inc_single_all

   ! armed_here() - the paradigm guard (single-fragment molecular; the atom
   !                list must be assembled - the mirror of the dispatch rule)
   logical function armed_here()
      armed_here = allocated(list_atoms%frag) .and. reactants%surface_model == 0 &
                   .and. size(list_atoms%frag) <= 1
   end function armed_here
end module reg_inc_single
