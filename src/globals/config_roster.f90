!=====================================================================
! config_roster.f90 - representation layer: roster - the system-composition table (masses/symbols/
!                     fragment tables/unit cell),
!                     including assembly roster_load and clearing roster_clear
! Reference: venus_input.f90:248-256 (mass section) / :508-518 (fragment-A table assembly) /
!       :743-766 (fragment-B table,
!       B side finalized in phase C) / :996-1000 (unit-cell parameters); src_VENUS/SURF.f:58-83
!       (surface unit-cell derivation);
!       src_VENUS/SELECT.f90:713-715, 1028-1030 (equilibrium geometry QZ assembly);
!       venus_data.f90:349, 585-586 (unit-cell / lattice-constant declarations and defaults)
!=====================================================================
module config_roster
   implicit none
   private
   public :: roster_t, roster, roster_load, roster_clear

   ! Roster: the system-composition table - per-atom element symbol, mass, and fragment membership,
   ! plus fragment geometry and unit-cell parameters (K1 plan terminology; fields registered per K1
   ! plan §1.4)
   type :: roster_t
      real(8), allocatable :: mass(:)          ! per-atom masses [amu] (old W; rigid-surface atoms
                                               ! fixed with a huge mass ~1d30 - venus_input.f90:752)
      character(len=4), allocatable :: symb(:) ! per-atom element symbols [-] (new field; the old
                                               ! library stored masses only, symbols read straight
      ! from input - E14: the symbols arrive from the system-definition files, the sysdef
      ! facade products (spec §6))
      real(8) :: mass_a = 0.0d0                ! fragment-A total mass [amu] (old WTA(1) = sum W(LA))
      integer, allocatable :: list_a(:)        ! fragment-A atom index list [global atom index] (old LA(1,:))
      real(8), allocatable :: qz_a(:)          ! fragment-A equilibrium geometry (principal-axis coordinates) [Å]
      ! (old QZA(1,:), length 3*nat_a)
      integer :: nat_a = 0                     ! fragment-A atom count [count] (old NATOMA(1))
      real(8), allocatable :: q0(:)            ! whole-system equilibrium geometry [Å] (old QZ;
      ! assembled from fragment geometries by atom index - SELECT.f90:713-715)
      real(8) :: skew = 0.0d0                  ! unit-cell skew angle [rad] (old SKEW; input
      ! degrees -> radians, venus_input.f90:999)
      real(8) :: box_lx = 0.0d0                ! unit-cell x edge length [Å] (old BOXLX)
      real(8) :: box_ly = 0.0d0                ! unit-cell y edge length [Å] (old BOXLY)
      real(8) :: box_pad = 2.0d0               ! half-width of the random cluster drop box [Å]
      ! (old RND_BOX, venus_input.f90:307)
      real(8) :: a_lat = 0.0d0                ! surface lattice constant [Å] (system parameter - no
      ! generic default: the old-library Au value 2.94998 from venus_data.f90:585 is a system bias and
      ! is not carried over; surface paradigms must supply it via input/blackbox - roster_load
      ! validation intercepts a missing value)
      ! Left blank: fragment-B table (nat_b/list_b/mass_b/qz_b - old NATOMB/LB/WTB/QZB,
      ! including the surface
      ! four unit-cell atoms) and multi-fragment extension (K1 risk note: max_frag array
      ! growth with the roster) -
      ! finalized in phase C together with the surface members / incident channel
   end type

   type(roster_t) :: roster    ! roster instance (same instance pattern as the config module;
                               ! assembled by roster_load,
                               ! state%natoms = size(roster%mass) - wiring in phase C driver/input)
contains
   !------------------------------------------------------------------
   ! roster_load() - assemble and validate the roster (E14 spec §6/§9: data source = the sysdef
   !                 facade products - the system-definition files parsed by sysdef_load
   !                 (files -> sysdef products -> roster); was the input keywords
   !                 (ATOM_MASSES / fragment table / SKEW/BOXLX/BOXLY), which no longer carry
   !                 system composition. roster_load shrinks to allocation + algebra +
   !                 validation: per-field fill from the sysdef products, fragment mass sums,
   !                 equilibrium-geometry placement, unit-cell derivation. Ordering statute:
   !                 sysdef_load runs FIRST - the driver assembly choreography (Task 12); a
   !                 roster_load ahead of the sysdef products is an assembly-order error)
   ! Reference: venus_input.f90:248-256, 508-518, 743-766, 996-1000 (the OLD input-keyword
   !       source - superseded as source by spec §6; kept as citation for the assembly
   !       algebra); src_VENUS/SURF.f:58-83; src_VENUS/SELECT.f90:713-715, 1028-1030;
   !       src/facade/sysdef_facade.f90 (sysdef_load - the data source; its step-4 fill notes)
   subroutine roster_load()
      ! Pseudocode:
      ! 1. Mass table: allocate mass/symb by atom count; fill from the sysdef products (per-atom
      !    element symbols + masses, the symbol->mass lookup done in sysdef_load; the old input-
      !    keyword two-value broadcast rule :248-256 is superseded as data source)
      ! 2. Fragment-A table: list_a = 1..nat_a in order, mass_a = sum(mass(list_a)),
      !    qz_a reads the equilibrium geometry (a sysdef product; old :508-518 for the algebra)
      ! 3. Rigid-surface paradigm: surface atoms fixed with a huge mass ~1d30 (old :752);
      !    fragment-B table left blank (phase C; old :743-766 for reference only)
      ! 4. q0 assembly: place each fragment's equilibrium geometry into the whole-system
      !    geometry by atom index
      !    (old SELECT.f90:713-715 filled by A, :1028-1030 by B)
      ! 5. Unit cell: the periodic-box paradigm reads skew (degrees -> radians) / box_lx /
      !    box_ly from the sysdef products (the cell-carrying system-definition files; the old
      !    input keys :996-1000 superseded as source);
      !    the rigid-surface paradigm derives box_lx/box_ly/skew from the four surface-atom
      !    coordinates (law of cosines,
      !    old SURF.f:58-83; a degenerate unit cell reports an error and stops, old :74-79)
      !    - B table wired in phase C
      ! 6. Validation: nat_a <= atom count, list_a indices in range, all masses positive
      !    (illegal -> loud-stop); surface paradigm (config%reactants%surface_model) with a_lat <= 0
      !    -> loud-stop naming the missing lattice constant (a system parameter with no generic
      !    default - the 2026-09-08 audit removed the old-library Au bias 2.94998)
   end subroutine roster_load

   !------------------------------------------------------------------
   ! roster_clear() - clear the roster (deallocate all allocatable members and reset scalar defaults)
   ! Reference: venus_data.f90:585-586 (unit-cell / drop-box / lattice-constant defaults are
   ! the reset targets)
   subroutine roster_clear()
      ! Pseudocode:
      ! 1. Deallocate mass/symb/list_a/qz_a/q0
      ! 2. Reset scalars: mass_a=0, nat_a=0, skew=0, box_lx=0, box_ly=0,
      !    box_pad=2.0, a_lat=0 (generic defaults; a_lat has no system default - see the
      !    field statute above)
   end subroutine roster_clear
end module config_roster
