!=====================================================================
! sampler.f90 - L2 interface-B' framework: general initial-state distribution sampler (distribution-member
!                registration + scheme dispatch - distributions -> (Q0,P0,s0); the SELECT integer dispatch
!                becomes a member table)
!                (see previous line)
! Reference: src_VENUS/SELECT.f90:71-87 (top gating - old task fork for read-in direct coordinates),
!       :89-142 (fragment placement + DH dissociation-energy reference - energy to the facade, geometry to
!       the roster, excised),
!       :144-186 (fragment A/B initialization dispatch - where the old integer-switch dispatch becomes a
!       member table),
!       :188-207 (collision-assembly dispatch: surface/gas two routes), :1300-1334 (setup_surface_collision),
!       :1337-1404 (setup_gas_collision - the two routes merge into the incident-channel member),
!       :1407-1491 (finalize_and_report - reporting to the recorder, excised);
!       distribution member roster (X2/X3): boltzmann/ebk/incident (batch-3 skeleton) + 9 more -> C5
!
! -- E14 sysdef hand-off (user ruling 2026-09-08; spec §6 data-flow tail / §9: sampler KEPT,
!    comment only - mechanism unchanged) --
!   Initial geometries and fragment data arrive from the sysdef hand-off: sysdef_load parses
!   the container folder's system-definition files and hands the equilibrium geometries /
!   fragment masses over BY POPULATING THE ROSTER (spec §6: "hand equilibrium geometries to
!   sampler (B' placement)"). The B' members (dist_boltzmann / dist_ebk / dist_incident)
!   read equilibrium geometries and fragment masses THROUGH THE ROSTER - the data itself is
!   the hand-off ("geometry through the roster" in samp_run step 4 already read this way).
!   Layering statute: NO module dependency on sysdef_facade - this module carries no use
!   statement on it and never calls it (the sampler never uses sysdef_facade; the roster is
!   the one and only geometry-knowledge channel into the B' family)
!   (see previous line)
!=====================================================================
module sampler
   use state, only: state_t
   implicit none
   private
   public :: samp_reg, samp_run

   ! Distribution member interface: members sample (q0,p0,s0) from their own distribution and write the
   ! components they own into the physical state
   abstract interface
      subroutine dist_proc_i(sta)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! initial state (the q/p/s components owned by the member are
         ! overwritten by sampling)
      end subroutine dist_proc_i
   end interface

   ! Distribution member table (B' registry family - a new member = a new file + one registration line, T2)
   integer, parameter :: max_dist = 16   ! member-table capacity [row]
   type :: dist_entry_t
      character(len=32) :: name = ''     ! member name (the scheme key of samp_run)
      procedure(dist_proc_i), pointer, nopass :: proc => null() ! member sampling routine
   end type dist_entry_t
   type(dist_entry_t) :: dist_tbl(max_dist) ! distribution member table (module-private; filled in C5)
   integer :: n_dist = 0                   ! registered row count [row]
contains
   !------------------------------------------------------------------
   ! samp_reg(name, proc) - register a distribution member into the member table
   ! Reference: docs/final_completion_plan.md §7 X2/X3 (SELECT as a member table - the 12-member roster)
   subroutine samp_reg(name, proc)
      character(len=*), intent(in) :: name  ! member name (e.g. 'boltzmann'/'ebk'/'incident')
      procedure(dist_proc_i) :: proc        ! member sampling routine
      ! Pseudocode:
      ! 1. Duplicate name or full table -> loud-stop
      ! 2. Register the row {name, proc} (executed in the assembly phase by the driver;
      !    each dist_* member registers one line)
   end subroutine samp_reg

   !------------------------------------------------------------------
   ! samp_run(scheme, sta) - dispatch sampling by scheme name: call the member to write initial-state
   ! components
   ! Reference: src_VENUS/SELECT.f90:144-186 (fragment dispatch), :188-207 (collision
   ! assembly dispatch)
   subroutine samp_run(scheme, sta)
      character(len=*), intent(in) :: scheme ! sampling scheme name (the config%reactants%dist_scheme key)
                                             ! the member name it denotes - old integer
                                             ! dispatch -> member table, X3;
      ! integer-code -> member-name mapping wired in C5)
      type(state_t), intent(inout) :: sta    ! initial state (q0/p0/s0 - overwritten by the member's sampling)
      ! Pseudocode:
      ! 1. Look up the scheme member in dist_tbl (unregistered -> loud-stop naming it)
      !    - the distribution member table is filled in C5 (currently a bare frame:
      !    batch-3 examples boltzmann/ebk/incident
      !    skeleton the first three; the remaining 9 members surface_oscillator/rotation/j/
      !    barrier_excitation/thermostat/glo_target/thermalize/normalmode/qmmicro -> C5 - E3)
      ! 2. Call the member proc(sta): sample the initial-state components it owns from its
      !    own distribution (member = code + its own data; the rng/linalg/geometry
      !    instruments are usable internally - preparation members calling evolution
      !    services internally is legal, boundary 4)
      ! 3. Multi-fragment systems call the members in order per config%reactants%
      !    dist_scheme (per-fragment distributions, frag-indexed; E15) + the incident-channel member
      !    (collision energy / impact parameter / orientation - the old setup_gas_collision
      !    :1337-1404 / setup_surface_collision :1300-1334 two routes merged)
      ! 4. The old SELECT master orchestration is excised: fragment placement and the DH
      !    dissociation-energy reference (:89-142) - energy goes through the facade,
      !    geometry through the roster; finalize_and_report (:1407-1491) - belongs to
      !    the recorder's statistics columns; the sampler only does distribution -> (q0,p0,s0)
      !    (see previous line)
   end subroutine samp_run
end module sampler
