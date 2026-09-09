!=====================================================================
! sysdef_facade.f90 - L2 system-definition facade (E14 spec §6): one assembly-phase entry,
!                    sysdef_load(scan_dir) — folder discovery, parser dispatch, mass
!                    lookup, roster fill, paradigm derivation, sampler hand-off
! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §6 (granularity: one
!       file per reactant, file base name = reactant name; discovery: assembly-phase scan
!       of the blackbox folder for known extensions, no manifest layer; formats: POSCAR
!       (with cell) / .xyz (no cell) initially, a new format = a new parser member,
!       grammar versioning; derivation rule: cell present -> surface, absent -> molecular
!       - a generic formal rule lawfully in generic code; data flow: scan -> generic
!       parser (utils) -> element symbols against the generic mass table -> fill roster
!       (geometry / cell / masses) -> derive paradigm -> hand equilibrium geometries to
!       sampler (B' placement); ownership: file contents = blackbox material, parsers /
!       mass table / derivation rule = generic code), §3 (assembly-phase entries:
!       sysdef_load is consumed by the driver assembly phase - no evolution-loop
!       presence), §9 (affected-files rows: config surface_model derived-not-input;
!       config_roster roster_load data source = sysdef facade products);
!       src/utils/specio.f90 (the readers + mass table dispatched into);
!       src/globals/config_roster.f90 (roster_t - the fill-target fields);
!       venus2026/src_VENUS/SELECT.f90:1300-1334 setup_surface_collision / :1337-1404
!       setup_gas_collision (citation anchor: the two routes the derived paradigm selects
!       - the old paradigm fork was code, the new one is data)
!
! -- Data-material contract, not a fifth entry (spec §3/§6) --
!   The evolution loop's four entries (force_eval / elec_prop / mqc_statistic / predicate)
!   run per step; sysdef_load runs ONCE in the assembly phase and produces DATA (roster
!   rows + the derived paradigm code), not callable physics. The geometric definition of
!   the reactants is system knowledge - blackbox material in the E11/E12 sense: drag in a
!   folder and the system is fully defined (no manifest to edit, no registration lines).
!   (see previous line)
!
! -- Ownership split (spec §6 ownership row; E12/E13 audit line) --
!   Generic here (zero system knowledge): the extension scan, per-extension dispatch,
!   symbol->mass lookup, the roster-fill mechanics, the cell-presence paradigm rule. The
!   system author owns only the file CONTENTS (geometries, cells, elements). Nothing in
!   this module branches on system identity; the paradigm code it writes is a formal fact
!   of the files (cell present/absent), not a physics judgment.
!   (see previous line)
!
! -- Field notes (grammar lands C4; the field statutes are Task 11's) --
!   scan_dir is a character dummy of this entry; the config-side path field (system_dir -
!   the blackbox folder path as config data) is noted here in comments only. The paradigm
!   field reactants%surface_model is an integer code (config.f90: 0 = gas phase, 1/2 =
!   the two surface paradigms): the molecular/surface fork is derivable from file shape
!   and is written here (derived, not input - spec §9); the 1-vs-2 surface-paradigm
!   selection is NOT derivable from cell presence alone and arrives with further input
!   (C4 grammar; roster_load's rigid-surface vs periodic-box assembly branches there).
!   (see previous line)
!=====================================================================
module sysdef_facade
   use config,        only: reactants
   use config_roster, only: roster
   use specio,        only: read_poscar, read_xyz, elem_mass
   implicit none
   private
   public :: sysdef_load
contains
   !------------------------------------------------------------------
   ! sysdef_load(scan_dir) - assembly-phase materialization of the system definition (the
   !              spec §6 data flow end to end; consumed once by the driver assembly phase)
   ! Reference: spec §6 (data flow row); §3 (assembly-phase entries); §9 (config_roster
   !       data source = sysdef facade products)
   subroutine sysdef_load(scan_dir)
      character(len=*), intent(in) :: scan_dir  ! system-definition folder (the blackbox
                                                ! folder; the config-side path field
                                                ! system_dir reaches this dummy at the C4
                                                ! driver wiring - grammar lands C4)
      ! Pseudocode:
      ! 1. Discovery (spec §6: no manifest layer): scan scan_dir for known extensions
      !    (the exact name POSCAR and the extensions .poscar / .xyz, case-insensitive);
      !    file base name = reactant name (A.poscar, B.xyz, ...); multiple files per
      !    extension legal, one reactant each; zero known-extension files -> loud-stop
      !    (no system defined - an empty blackbox folder is an assembly error, not an
      !    empty system)
      ! 2. Per-file dispatch by extension: POSCAR/.poscar -> read_poscar (the cell is
      !    retained), .xyz -> read_xyz (no cell - by format definition). Each file opened
      !    here on a scratch unit, the reader called with the unit number, the unit closed
      !    (the readers take unit numbers, not names). A reactant name given in BOTH
      !    formats -> loud-stop naming it; a reader soft failure (ok=.false.) ->
      !    loud-stop echoing errmsg with the file name (the ok/errmsg channel exists for
      !    the readers' own testability; the assembly phase converts a failure into a hard
      !    stop - malformed blackbox material is unusable)
      !    (see previous line)
      ! 3. Masses: elem_mass(names(i)) per atom (symbol -> [amu]; an unknown symbol
      !    loud-stops inside elem_mass naming it)
      ! 4. Roster fill (spec §6 "fill roster (geometry / cell / masses)"; the fill target
      !    is config_roster's roster_t, assembly/validation staying in roster_load - the
      !    representation layer's own entry, its data source now the sysdef products,
      !    spec §9): per-atom roster%mass / roster%symb (files concatenated in name
      !    order); fragment tables one file per fragment - roster%frag(i) for the i-th file
      !    (the homogeneous frag_t row: list/nat/mass/qz; E15 arrayed, no A/B asymmetry -
      !    the surface rows' extra content is config_roster's own phase-C blank, field note
      !    there); equilibrium geometry roster%q0 (fragment geometries placed by atom
      !    index); unit cell roster%box_lx / roster%box_ly / roster%skew for cell-carrying
      !    files (+ a_lat is a system parameter with no generic default - roster_load
      !    validation intercepts a missing value). The exact conveyance grammar between
      !    the parsed material and roster_load (argument list vs staging tables) is
      !    finalized in C4 with the input grammar; the K2 statute is the direction:
      !    files -> sysdef products -> roster
      !    (see previous line)
      ! 5. Paradigm derivation (spec §6 derivation rule - generic and lawful, the same
      !    class as E13's retained n_surf>=2-with-r_scalar check): any file carries a
      !    cell -> reactants%surface_model = a surface paradigm code (1/2); none -> 0
      !    (gas phase / molecular)
      ! 6. Sampler hand-off (spec §6 "hand equilibrium geometries to sampler (B'
      !    placement)"): the equilibrium geometries + fragment masses are now ROSTER DATA
      !    - the B' members (dist_boltzmann / dist_ebk / dist_incident) read them through
      !    the roster; this module takes no dependency on sampler (no use statement - the
      !    hand-off is the data itself). The citation-anchor routes above
      !    (setup_surface_collision / setup_gas_collision) are what the derived paradigm
      !    now selects through data instead of a code fork
      !    (see previous line)
      ! 7. Validation loud-stop (formal): mixed surface/molecular definitions across
      !    files (some with cell, some without) -> loud-stop naming the conflict - the
      !    step-5 rule demands homogeneity, a half-periodic system is no paradigm the
      !    closed set defines
   end subroutine sysdef_load
end module sysdef_facade
