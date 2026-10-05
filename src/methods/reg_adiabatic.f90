!=====================================================================
! reg_adiabatic.f90 - register the 'adiabatic' electronic-method package
! Design:
!   One method_reg call fills the whole 'adiabatic' row; the reg file is
!   the atomic registration unit (a new method = its method file + one
!   thin reg file; zero edits to existing files).
!   The package is the all-no-op degenerate built-in: no amplitudes to
!   propagate, occ stays [1] constant, no sync member (the null pointer
!   is the legal 'no deliverable cache' form).
!=====================================================================
module reg_adiabatic
   use elec_interface,   only: method_reg
   use elec_adiabatic, only: ad_prop, ad_match   ! package members; their conformance
   ! to the abstract interfaces is checked at the method_reg call site
   implicit none
   private
   public :: reg_adiabatic_all
contains
   !------------------------------------------------------------------
   ! reg_adiabatic_all() - fill the 'adiabatic' package row (one call, one row)
   subroutine reg_adiabatic_all()
      call method_reg('adiabatic', ad_prop, ad_match)
   end subroutine reg_adiabatic_all
end module reg_adiabatic
