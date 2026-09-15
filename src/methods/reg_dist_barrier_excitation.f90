!=====================================================================
! reg_dist_barrier_excitation.f90 - registration module (methods/): ONE
!                          samp_reg call registers the
!                          'barrier_excitation' distribution member row
!                          - the reg file is the atomic registration unit
!                          (a new member = its method file + this one
!                          thin file; zero edits to existing files)
! Contract:
!   reg_dist_barrier_excitation_all() - register the barrier-excitation
!     member's sampling routine under its scheme word
!     'barrier_excitation' (the word the input grammar accepts for
!     DIST_SCHEME and the per-fragment dispatch resolves). Passing
!     barrier_excitation_sample as samp_reg's procedure dummy checks its
!     conformance to the dist_proc_i interface at compile time. Executed
!     once in the assembly phase, before any samp_run.
!   reg_dist_barrier_excitation_init() - the ASSEMBLY SEAM entry: load
!     this member's staged input parameters (the input framework's
!     public staging rows) into the member's params receiver - absent
!     keys keep the field defaults - then run the member's own init
!     (validation + carrier binding). Guard: an UNSELECTED member skips
!     (no dist_scheme row carries its word - member init's carrier
!     battery would otherwise loud-stop 'zero carriers' on every run
!     that does not use the member). The mode-frequency and
!     mode-coefficient tables w_mode/c_mode are ASSEMBLY DATA, not input
!     keys - no production source stages them yet, so a selected
!     barrier_excitation loud-stops at its own unallocated-table check
!     (the fail-loud posture). Three-party zero knowledge: the driver
!     knows no member, the input knows no member, the member never
!     touches the driver.
! Units: none (registration carries names and code references only).
!=====================================================================
module reg_dist_barrier_excitation
   use sampler,                only: samp_reg, samp_code
   use config,                 only: reactants
   use input,                  only: staging_n_rows, staging_key, staging_val, staging_line
   use seam_cast,              only: seam_real, seam_int
   use dist_barrier_excitation, only: barrier_excitation_sample, &
                                      barrier_excitation_params_t, barrier_excitation_init
   implicit none
   private
   public :: reg_dist_barrier_excitation_all
   public :: reg_dist_barrier_excitation_init
contains
   !------------------------------------------------------------------
   ! reg_dist_barrier_excitation_all() - fill the member-table row (one call, one row)
   subroutine reg_dist_barrier_excitation_all()
      call samp_reg('barrier_excitation', barrier_excitation_sample)
   end subroutine reg_dist_barrier_excitation_all

   !------------------------------------------------------------------
   ! reg_dist_barrier_excitation_init() - the assembly seam: staged member keys ->
   !                                      params receiver -> the member's own init
   subroutine reg_dist_barrier_excitation_init()
      type(barrier_excitation_params_t) :: p   ! member-parameter receiver (field defaults)
      character(len=512) :: sval               ! one staged raw value
      integer :: i
      if (.not. allocated(reactants%dist_scheme)) return
      if (count(reactants%dist_scheme == samp_code('barrier_excitation')) == 0) return
      do i = 1, staging_n_rows()
         sval = staging_val(i)
         select case (trim(staging_key(i)))
         case ('E_STAB');  call seam_real(staging_key(i), sval, staging_line(i), p%e_stab)
         case ('N_E_BAR'); call seam_int(staging_key(i), sval, staging_line(i), p%n_e_bar)
         case ('E_BAR');   call seam_real(staging_key(i), sval, staging_line(i), p%e_bar)
         case ('T_BAR');   call seam_real(staging_key(i), sval, staging_line(i), p%t_bar)
         end select
      end do
      call barrier_excitation_init(p)
   end subroutine reg_dist_barrier_excitation_init
end module reg_dist_barrier_excitation
