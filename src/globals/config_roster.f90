!=====================================================================
! config_roster.f90 - representation layer: roster - the system-composition table (masses/symbols/
!                     fragment tables/unit cell),
!                     including assembly roster_load and clearing roster_clear
! Reference: venus_input.f90:248-256 (mass section) / :508-518 (fragment-A table assembly) /
!       :743-766 (fragment-B table - reference for the phase-C surface rows) /
!       :996-1000 (unit-cell parameters); src_VENUS/SURF.f:58-83
!       (surface unit-cell derivation);
!       src_VENUS/SELECT.f90:713-715, 1028-1030 (equilibrium geometry QZ assembly);
!       venus_data.f90:349, 585-586 (unit-cell / lattice-constant declarations and defaults);
!       src/facade/sysdef_facade.f90 (the staging accessors - the S2 transfer source;
!       its header self-documents the roster_load usage contract)
!=====================================================================
module config_roster
   use config, only: reactants
   use sysdef_facade, only: staging_loaded, staging_frag_count, staging_frag_nat, &
                            staging_atom_symbols, staging_atom_masses, staging_atom_geometry, &
                            staging_cell_present, staging_unit_cell
   implicit none
   private
   public :: roster_t, frag_t, roster, roster_load, roster_clear

   ! Roster: the system-composition table - per-atom element symbol, mass, and fragment membership,
   ! plus fragment geometry and unit-cell parameters (K1 plan terminology; fields registered per K1
   ! plan §1.4). Fragment table symmetric and arrayed (E15, review-gate ruling 2026-09-09):
   ! one homogeneous frag_t row per fragment - no A/B asymmetry, no 2-fragment cap (the K1
   ! plan's max_frag risk-note extension direction landed early at the review gate)
   type :: frag_t
      integer :: nat = 0                   ! fragment atom count [count] (old NATOMA(i)/NATOMB)
      real(8) :: mass = 0.0d0              ! fragment total mass [amu] (old WTA(i)/WTB = sum mass(list))
      integer, allocatable :: list(:)      ! fragment atom index list [global atom index] (old LA(i,:)/LB)
      real(8), allocatable :: qz(:)        ! fragment equilibrium geometry (principal-axis coordinates)
                                           ! [Å] (old QZA(1,:)/QZB, length 3*nat)
   end type

   type :: roster_t
      real(8), allocatable :: mass(:)          ! per-atom masses [amu] (old W; rigid-surface atoms
                                               ! fixed with a huge mass ~1d30 - venus_input.f90:752)
      character(len=4), allocatable :: symb(:) ! per-atom element symbols [-] (new field; the old
                                               ! library stored masses only, symbols read straight
      ! from input - E14: the symbols arrive from the system-definition files, the sysdef
      ! facade products (spec §6))
      type(frag_t), allocatable :: frag(:)     ! fragment table - one homogeneous row per fragment
                                               ! (E15; frag(1) carries the old fragment-A semantics;
                                               ! the surface rows' extra content is the phase-C blank,
                                               ! field note at the type end)
      real(8), allocatable :: q0(:)            ! whole-system equilibrium geometry [Å] (old QZ;
      ! assembled from fragment geometries by atom index - SELECT.f90:713-715)
      real(8) :: skew = 0.0d0                  ! unit-cell skew angle [rad] (old SKEW; the degrees ->
      ! radians conversion happens at fill time - the old-library input key :999 carried degrees;
      ! E14: the value arrives from the sysdef products or the surface-atom derivation,
      ! roster_load step 5)
      real(8) :: box_lx = 0.0d0                ! unit-cell x edge length [Å] (old BOXLX)
      real(8) :: box_ly = 0.0d0                ! unit-cell y edge length [Å] (old BOXLY)
      real(8) :: box_pad = 2.0d0               ! half-width of the random cluster drop box [Å]
      ! (old RND_BOX, venus_input.f90:307)
      real(8) :: a_lat = 0.0d0                ! surface lattice constant [Å] (system parameter - no
      ! generic default: the old-library Au value 2.94998 from venus_data.f90:585 is a system bias and
      ! is not carried over; surface paradigms must supply it via input/container - roster_load
      ! validation intercepts a missing value)
      ! Left blank (content, not form - E15 settled the form): what the surface-paradigm fragment
      ! rows carry beyond the homogeneous shape (the old fragment-B table held the four surface
      ! unit-cell atoms, venus_input.f90:743-766) -
      ! finalized in phase C together with the surface members / incident channel
   end type

   type(roster_t) :: roster    ! roster instance (same instance pattern as the config module;
                               ! assembled by roster_load,
                               ! state dimensions = 3*size(roster%mass) (E16: no cached
                               ! natoms/ndof fields - size-derived) - wiring in phase C driver/input)
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
   !       src/facade/sysdef_facade.f90 (sysdef_load - the data source; the S2 staging
   !       accessors consumed below)
   ! Implemented flow (S2 transfer, ruling A - one-way pull from the sysdef staging):
   !   0. guard: staging_loaded() = .false. -> loud-stop (assembly-order error)
   !   1. per-atom transfer: symb/mass/geometry pulled from the staging getters (the
   !      symbol->mass lookup already happened in sysdef_load - not repeated here);
   !      internal consistency of the staged lengths checked defensively
   !   2. fragment tables: frag(i) = the contiguous atom-index run of fragment i (rows in
   !      the sysdef name order): nat = staging_frag_nat(i), list = offset+1 .. offset+nat,
   !      mass = sum(mass(list)) (old :508-518 fragment-A algebra, :743-766 the B-side
   !      reference), qz = the run's geometry packed 3*nat (old QZA/QZB fill)
   !   3. q0 assembly: the per-fragment qz runs concatenated by atom index = the packed
   !      staged geometry (old SELECT.f90:713-715 by A, :1028-1030 by B)
   !   4. unit cell (periodic-box route): the staged 3x3 cell (lattice vectors BY ROW,
   !      the specio convention) -> box_lx = |a1|,
   !      box_ly = |a2|, skew = acos(a1.a2/(|a1||a2|)) in RADIANS (the law of cosines,
   !      old SURF.f:58-83; the degrees->radians conversion note concerns the old input
   !      route - from cell vectors the angle arrives in radians directly); a zero-length
   !      in-plane vector is a degenerate cell -> loud-stop (old :74-79 statute); no cell
   !      staged -> the box scalars stay at their defaults
   !   5. validation (D4) on the assembled table: per fragment nat <= atom count, list
   !      indices in range, all masses positive, surface paradigm with a_lat <= 0 ->
   !      loud-stop naming the missing lattice constant (a system parameter with no
   !      generic default - the 2026-09-08 audit removed the old-library Au bias 2.94998)
   ! Ownership note: roster_load re-fills only the DERIVED members (mass/symb/frag/q0/
   ! skew/box_lx/box_ly); a_lat and box_pad are input-owned (the residual input keys,
   ! C4.4) and carried through untouched - a_lat is validated in step 5, never written
   ! Phase-C blanks (field notes, not fill): the rigid-surface huge-mass fixing (~1d30,
   ! old :752) and the surface-row extra content - the staging API records the cell as a
   ! single system-level fact without per-fragment provenance, so which rows are surface
   ! rows arrives with the phase-C surface members / incident channel
   subroutine roster_load()
      character(len=4), allocatable :: symb(:)
      real(8), allocatable :: amass(:), geom(:,:)
      real(8) :: cell(3,3), box_lx, box_ly, cos_skew
      character(len=256) :: msg
      integer :: n_atoms, n_frag, i, j, k, offset

      ! 0. ordering-statute guard (S2 usage contract)
      if (.not. staging_loaded()) then
         call stop_assembly('no staged sysdef products - sysdef_load must run first '// &
                            '(assembly-order error, the config_roster ordering statute)')
      end if

      ! 1. per-atom transfer (fresh copies from the staging getters)
      symb  = staging_atom_symbols()
      amass = staging_atom_masses()
      geom  = staging_atom_geometry()
      n_atoms = size(symb)
      if (n_atoms <= 0 .or. size(amass) /= n_atoms .or. &
          size(geom,1) /= 3 .or. size(geom,2) /= n_atoms) then
         call stop_assembly('staged per-atom products inconsistent (empty or mismatched lengths)')
      end if
      n_frag = staging_frag_count()
      if (n_frag <= 0) call stop_assembly('staged fragment count is zero')

      call reset_derived()
      allocate (roster%symb(n_atoms), roster%mass(n_atoms), roster%q0(3*n_atoms))
      allocate (roster%frag(n_frag))
      roster%symb = symb
      roster%mass = amass

      ! 2. fragment tables: one contiguous run per fragment (old :508-518 algebra)
      offset = 0
      do i = 1, n_frag
         associate (fr => roster%frag(i))
            fr%nat = staging_frag_nat(i)
            if (fr%nat < 0) then
               write (msg, '(a,i0,a,i0)') 'staged fragment ', i, ' has a negative atom count ', fr%nat
               call stop_assembly(trim(msg))
            end if
            allocate (fr%list(fr%nat), fr%qz(3*fr%nat))
            do j = 1, fr%nat
               k = offset + j
               fr%list(j) = k
               fr%qz(3*j-2:3*j) = geom(:,k)
            end do
            fr%mass = sum(roster%mass(fr%list))
            offset = offset + fr%nat
         end associate
      end do
      if (offset /= n_atoms) then
         write (msg, '(a,i0,a,i0)') 'staged fragment atom counts sum to ', offset, &
            ' but the staged atom count is ', n_atoms
         call stop_assembly(trim(msg))
      end if

      ! 3. q0 assembly: the whole-system equilibrium geometry packed by atom index
      do k = 1, n_atoms
         roster%q0(3*k-2:3*k) = geom(:,k)
      end do

      ! 4. unit cell (periodic-box route: law of cosines on the staged cell vectors)
      if (staging_cell_present()) then
         cell = staging_unit_cell()
         box_lx = sqrt(dot_product(cell(1,:), cell(1,:)))
         box_ly = sqrt(dot_product(cell(2,:), cell(2,:)))
         if (box_lx <= 0.0d0 .or. box_ly <= 0.0d0) then
            write (msg, '(a,es12.4,a,es12.4,a)') 'degenerate unit cell: box_lx=', box_lx, &
               ' box_ly=', box_ly, ' - zero-length in-plane cell vector in the staged unit cell'
            call stop_assembly(trim(msg))
         end if
         cos_skew = dot_product(cell(1,:), cell(2,:))/(box_lx*box_ly)
         cos_skew = max(-1.0d0, min(1.0d0, cos_skew))
         roster%box_lx = box_lx
         roster%box_ly = box_ly
         roster%skew = acos(cos_skew)
      end if

      ! 5. validation (D4) on the assembled table
      do i = 1, size(roster%frag)
         if (roster%frag(i)%nat > n_atoms) then
            write (msg, '(a,i0,a,i0,a,i0)') 'fragment ', i, ' atom count ', roster%frag(i)%nat, &
               ' exceeds the system atom count ', n_atoms
            call stop_assembly(trim(msg))
         end if
         do j = 1, roster%frag(i)%nat
            k = roster%frag(i)%list(j)
            if (k < 1 .or. k > n_atoms) then
               write (msg, '(a,i0,a,i0,a,i0)') 'fragment ', i, ' list entry ', j, &
                  ' is out of range: ', k
               call stop_assembly(trim(msg))
            end if
         end do
      end do
      do k = 1, n_atoms
         if (roster%mass(k) <= 0.0d0) then
            write (msg, '(a,i0,a,es12.4)') 'atom ', k, ' carries a non-positive mass ', roster%mass(k)
            call stop_assembly(trim(msg))
         end if
      end do
      if (reactants%surface_model /= 0 .and. roster%a_lat <= 0.0d0) then
         call stop_assembly('surface paradigm (surface_model=1) is missing the surface lattice '// &
                            'constant (a_lat <= 0) - a system parameter with no generic default; '// &
                            'supply it via input')
      end if
   end subroutine roster_load

   !------------------------------------------------------------------
   ! roster_clear() - clear the roster (deallocate all allocatable members and reset scalar defaults)
   ! Reference: venus_data.f90:585-586 (unit-cell / drop-box / lattice-constant defaults are
   ! the reset targets)
   ! Implemented: mass/symb/frag (each row's list/qz go with the array)/q0 deallocated;
   ! skew/box_lx/box_ly = 0, box_pad = 2.0, a_lat = 0 (generic defaults; a_lat has no
   ! system default - see the field statute above)
   subroutine roster_clear()
      if (allocated(roster%mass)) deallocate (roster%mass)
      if (allocated(roster%symb)) deallocate (roster%symb)
      if (allocated(roster%frag)) deallocate (roster%frag)
      if (allocated(roster%q0)) deallocate (roster%q0)
      roster%skew = 0.0d0
      roster%box_lx = 0.0d0
      roster%box_ly = 0.0d0
      roster%box_pad = 2.0d0
      roster%a_lat = 0.0d0
   end subroutine roster_clear

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! reset_derived() - drop only the DERIVED members ahead of a re-fill (roster_load's
   ! internal reset; the input-owned scalars a_lat/box_pad are NOT touched here - the
   ! full reset including them is roster_clear, the public entry)
   subroutine reset_derived()
      if (allocated(roster%mass)) deallocate (roster%mass)
      if (allocated(roster%symb)) deallocate (roster%symb)
      if (allocated(roster%frag)) deallocate (roster%frag)
      if (allocated(roster%q0)) deallocate (roster%q0)
      roster%skew = 0.0d0
      roster%box_lx = 0.0d0
      roster%box_ly = 0.0d0
   end subroutine reset_derived

   ! stop_assembly(msg) - the loud-stop channel (same semantics as sysdef_facade C4.2:
   ! name the caller, print the named message, STOP 1)
   subroutine stop_assembly(msg)
      character(len=*), intent(in) :: msg
      write (0, '(a)') 'roster_load: assembly error: '//trim(msg)
      write (0, '(a)') 'roster_load: loud-stop (the roster is unusable)'
      stop 1
   end subroutine stop_assembly

end module config_roster
