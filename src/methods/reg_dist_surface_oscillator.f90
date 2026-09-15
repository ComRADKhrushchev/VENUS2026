!=====================================================================
! reg_dist_surface_oscillator.f90 - registration module (methods/): ONE
!                          samp_reg call registers the
!                          'surface_oscillator' distribution member row
!                          - the reg file is the atomic registration unit
!                          (a new member = its method file + this one
!                          thin file; zero edits to existing files)
! Contract:
!   reg_dist_surface_oscillator_all() - register the surface-oscillator
!     member's sampling routine under its scheme word
!     'surface_oscillator' (the word the input grammar accepts for
!     DIST_SCHEME and the per-fragment dispatch resolves). Passing
!     surface_oscillator_sample as samp_reg's procedure dummy checks its
!     conformance to the dist_proc_i interface at compile time. Executed
!     once in the assembly phase, before any samp_run.
!   reg_dist_surface_oscillator_init() - the ASSEMBLY SEAM entry: load
!     this member's staged input parameters (the input framework's
!     public staging rows) into the member's params receiver - absent
!     keys keep the field defaults - then run the member's own init
!     (validation + carrier binding). Guard: an UNSELECTED member skips
!     (no dist_scheme row carries its word - member init's carrier
!     battery would otherwise loud-stop 'zero carriers' on every run
!     that does not use the member). Three-party zero knowledge: the
!     driver knows no member, the input knows no member, the member
!     never touches the driver.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_surface_oscillator
   use sampler,                 only: samp_reg, samp_code
   use config,                  only: reactants
   use input,                   only: staging_n_rows, staging_key, staging_val
   use dist_surface_oscillator, only: surface_oscillator_sample, &
                                      surf_osc_params_t, surface_oscillator_init
   implicit none
   private
   public :: reg_dist_surface_oscillator_all
   public :: reg_dist_surface_oscillator_init
contains
   !------------------------------------------------------------------
   ! reg_dist_surface_oscillator_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_surface_oscillator_all()
      call samp_reg('surface_oscillator', surface_oscillator_sample)
   end subroutine reg_dist_surface_oscillator_all

   !------------------------------------------------------------------
   ! reg_dist_surface_oscillator_init() - the assembly seam: staged member keys ->
   !                                       params receiver -> the member's own init
   subroutine reg_dist_surface_oscillator_init()
      type(surf_osc_params_t) :: p     ! member-parameter receiver (field defaults)
      character(len=512) :: sval       ! one staged raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('surface_oscillator')) == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('N_OSC');   read (sval, *) p%n_osc
         case ('T_OSC');   read (sval, *) p%t_osc
         case ('N_LEVEL'); read (sval, *) p%n_level
         case ('K_OSC');   read (sval, *) p%k_osc
         end select
      end do
      call surface_oscillator_init(p)
   end subroutine reg_dist_surface_oscillator_init
end module reg_dist_surface_oscillator
