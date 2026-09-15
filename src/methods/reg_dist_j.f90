!=====================================================================
! reg_dist_j.f90 - registration module (methods/): ONE samp_reg call
!                  registers the 'j' distribution member row - the reg
!                  file is the atomic registration unit (a new member =
!                  its method file + this one thin file; zero edits to
!                  existing files)
! Contract:
!   reg_dist_j_all() - register the j member's sampling routine under
!     its scheme word 'j' (the word the input grammar accepts for
!     DIST_SCHEME and the per-fragment dispatch resolves). Passing
!     j_sample as samp_reg's procedure dummy checks its conformance to
!     the dist_proc_i interface at compile time. Executed once in the
!     assembly phase, before any samp_run.
!   reg_dist_j_init() - the ASSEMBLY SEAM entry: load this member's
!     staged input parameters (the input framework's public staging
!     rows) into the member's params receiver - absent keys keep the
!     field defaults - then run the member's own init (validation +
!     carrier binding). Guard: an UNSELECTED member skips (no
!     dist_scheme row carries its word - member init's carrier battery
!     would otherwise loud-stop 'zero carriers' on every run that does
!     not use the member). Three-party zero knowledge: the driver knows
!     no member, the input knows no member, the member never touches
!     the driver.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_j
   use sampler, only: samp_reg, samp_code
   use config,  only: reactants
   use input,   only: staging_n_rows, staging_key, staging_val
   use dist_j,  only: j_sample, j_params_t, j_init
   implicit none
   private
   public :: reg_dist_j_all
   public :: reg_dist_j_init
contains
   !------------------------------------------------------------------
   ! reg_dist_j_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_j_all()
      call samp_reg('j', j_sample)
   end subroutine reg_dist_j_all

   !------------------------------------------------------------------
   ! reg_dist_j_init() - the assembly seam: staged member keys ->
   !                     params receiver -> the member's own init
   subroutine reg_dist_j_init()
      type(j_params_t) :: p            ! member-parameter receiver (field defaults)
      character(len=512) :: sval       ! one staged raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('j')) == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('T_ROT'); read (sval, *) p%t_rot
         end select
      end do
      call j_init(p)
   end subroutine reg_dist_j_init
end module reg_dist_j
