!=====================================================================
! reg_dist_ebk.f90 - registration module (methods/): ONE samp_reg call
!                    registers the 'ebk' distribution member row - the
!                    reg file is the atomic registration unit (a new
!                    member = its method file + this one thin file;
!                    zero edits to existing files)
! Contract:
!   reg_dist_ebk_all() - register the ebk member's sampling routine
!     under its scheme word 'ebk' (the word the input grammar accepts
!     for DIST_SCHEME and the per-fragment dispatch resolves). Passing
!     ebk_sample as samp_reg's procedure dummy checks its conformance
!     to the dist_proc_i interface at compile time. Executed once in
!     the assembly phase, before any samp_run.
!   reg_dist_ebk_init() - the ASSEMBLY SEAM entry (this member's shape:
!     NO input key stages an ebk parameter - the spectroscopic constants
!     w_e/w_ex_e/b_rot have no production source yet, so the receiver
!     keeps its field defaults and a SELECTED ebk loud-stops at its own
!     w_e check - the fail-loud posture; the routine exists so the seam
!     is uniform across the member family). Guard: an UNSELECTED member
!     skips (no dist_scheme row carries its word). Three-party zero
!     knowledge: the driver knows no member, the input knows no member,
!     the member never touches the driver.
! Units: none (registration carries names and code references only).
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
   ! reg_dist_ebk_init() - the assembly seam: no staged keys, the member's
   !                        own init runs on field defaults (fail-loud)
   subroutine reg_dist_ebk_init()
      type(ebk_params_t) :: p            ! member-parameter receiver (field defaults;
                                         ! no input key stages an ebk parameter)
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('ebk')) == 0) return
      call ebk_init(p)
   end subroutine reg_dist_ebk_init
end module reg_dist_ebk
