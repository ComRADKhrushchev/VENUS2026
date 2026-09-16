!=====================================================================
! reg_dist_ebk.f90 - register the 'ebk' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam with no staged keys: no input key stages
!   an ebk parameter (the spectroscopic constants have no production
!   source), so a selected member loud-stops at its own w_e check
!   (fail-loud); an unselected member skips.
!=====================================================================
module reg_dist_ebk
   use sampler,   only: samp_reg, samp_code
   use config,    only: reactants
   use dist_ebk,  only: ebk_sample, ebk_params_t, ebk_init
   implicit none
   private
   public :: reg_dist_ebk_all
   public :: reg_dist_ebk_init
contains
   !------------------------------------------------------------------
   ! reg_dist_ebk_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_ebk_all()
      call samp_reg('ebk', ebk_sample)
   end subroutine reg_dist_ebk_all

   !------------------------------------------------------------------
   ! reg_dist_ebk_init() - assembly seam: no staged keys, the member's own
   !                        init runs on field defaults (fail-loud)
   subroutine reg_dist_ebk_init()
      type(ebk_params_t) :: p            ! member-parameter receiver (field defaults;
                                         ! no input key stages an ebk parameter)
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('ebk')) == 0) return
      call ebk_init(p)
   end subroutine reg_dist_ebk_init
end module reg_dist_ebk
