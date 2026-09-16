!=====================================================================
! config_roster.f90 - system-composition table (roster): per-atom symbol/mass,
!                     fragment table, equilibrium geometry, unit cell
! Design:
!   roster_load pulls one-way from the sysdef facade staging (sysdef_load must run
!   first - assembly-order error otherwise): per-atom transfer, contiguous fragment
!   runs, q0 packing, unit-cell derivation (law of cosines, radians), validation,
!   then the A_LAT/BOXPAD window reconcile - their single write site (molecular
!   paradigm forbids both keys, surface paradigm requires both).
!   a_lat/box_pad are input-owned window scalars; all other members are derived.
!=====================================================================
module config_roster
   use config, only: reactants
   use sysdef_facade, only: staging_loaded, staging_frag_count, staging_frag_nat, &
                            staging_atom_symbols, staging_atom_masses, staging_atom_geometry, &
                            staging_cell_present, staging_unit_cell
   use input, only: staging_n_rows, staging_key, staging_val, staging_line
   implicit none
   private
   public :: roster_t, frag_t, roster, roster_load, roster_clear

   ! Fragment table symmetric and arrayed: one homogeneous frag_t row per
   ! fragment - no A/B asymmetry, no 2-fragment cap
   type :: frag_t
      integer :: nat = 0                   ! fragment atom count [count]
      real(8) :: mass = 0.0d0              ! fragment total mass [amu] (= sum of member masses)
      integer, allocatable :: list(:)      ! fragment atom index list [global atom index]
      real(8), allocatable :: qz(:)        ! fragment equilibrium geometry (principal-axis coordinates)
                                           ! [Å], length 3*nat
   end type

   type :: roster_t
      real(8), allocatable :: mass(:)          ! per-atom masses [amu] (rigid-surface atoms fixed
                                               ! with a huge mass ~1d30)
      character(len=4), allocatable :: symb(:) ! per-atom element symbols [-] (arrive from the
                                               ! system-definition files via the sysdef facade
                                               ! products)
      type(frag_t), allocatable :: frag(:)     ! fragment table - one homogeneous row per fragment;
                                               ! the surface rows' extra content is still blank
      real(8), allocatable :: q0(:)            ! whole-system equilibrium geometry [Å] (assembled
                                                ! from the fragment geometries by atom index)
      real(8) :: skew = 0.0d0                  ! unit-cell skew angle [rad] (derived from the
                                                ! staged cell vectors, in radians)
      real(8) :: box_lx = 0.0d0                ! unit-cell x edge length [Å]
      real(8) :: box_ly = 0.0d0                ! unit-cell y edge length [Å]
      real(8) :: box_pad = 2.0d0               ! half-width of the random cluster drop box [Å]
                                                ! (input-owned through the BOXPAD key - a
                                                ! surface-window key, written only by the
                                                ! roster_load window reconcile, its single
                                                ! write site)
      real(8) :: a_lat = 0.0d0                ! surface lattice constant [Å] (system parameter with
                                               ! no generic default: surface paradigms must supply
                                               ! the A_LAT input key - the roster_load window
                                               ! reconcile intercepts a missing value and is its
                                               ! only write site)
      ! Left blank (content, not form): what the surface-paradigm fragment rows carry
      ! beyond the homogeneous shape - finalized together with the surface members /
      ! incident channel
   end type

   type(roster_t) :: roster    ! roster instance (same instance pattern as the config module;
                               ! assembled by roster_load; state dimensions =
                               ! 3*size(roster%mass) - no cached natoms/ndof fields,
                               ! size-derived)
contains
   !------------------------------------------------------------------
   ! roster_load() - assemble and validate the roster from the sysdef staging
   !                 products (one-way pull; sysdef_load must run first -
   !                 assembly-order error otherwise); flow steps numbered below
   !------------------------------------------------------------------
   subroutine roster_load()
      character(len=4), allocatable :: symb(:)
      real(8), allocatable :: amass(:), geom(:,:)
      real(8) :: cell(3,3), box_lx, box_ly, cos_skew
      character(len=256) :: msg
      integer :: n_atoms, n_frag, i, j, k, offset

      ! 0. ordering-statute guard (assembly-order error)
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

      ! 2. fragment tables: one contiguous run per fragment
      offset = 0
      do i = 1, n_frag
         associate (fr => roster%frag(i))
            fr%nat = staging_frag_nat(i)
            if (fr%nat < 0) then
               write (msg, '(a,i0,a,i0)') 'staged fragment ', i, ' has a negative atom count ', fr%nat
               call stop_assembly(trim(msg))
            end if
            ! bounds guard BEFORE any allocate/geom access: this fragment's contiguous
            ! run must fit inside the staged atom count (a too-large nat loud-stops
            ! named here - never an out-of-bounds geom read; fail-safe defensive order)
            if (offset + fr%nat > n_atoms) then
               write (msg, '(a,i0,a,i0,a,i0,a,i0)') 'fragment ', i, ' atom count ', fr%nat, &
                  ' exceeds the remaining staged atoms (', n_atoms - offset, ' left)'
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

      ! 5. validation on the assembled table
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

      ! 6. surface-window reconcile (the paradigm condition of A_LAT/BOXPAD - the
      !    single write site of the window scalars)
      call window_reconcile()
   end subroutine roster_load

   !------------------------------------------------------------------
   ! roster_clear() - deallocate all allocatable members and reset the scalars
   !                   to their generic defaults
   !------------------------------------------------------------------
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

   ! reset_derived() - drop the DERIVED members ahead of a re-fill; the window scalars
   !                    reset to generic defaults too, so a re-load cannot inherit the
   !                    previous load's window values (step 6 re-pulls from the staging)
   subroutine reset_derived()
      if (allocated(roster%mass)) deallocate (roster%mass)
      if (allocated(roster%symb)) deallocate (roster%symb)
      if (allocated(roster%frag)) deallocate (roster%frag)
      if (allocated(roster%q0)) deallocate (roster%q0)
      roster%skew = 0.0d0
      roster%box_lx = 0.0d0
      roster%box_ly = 0.0d0
      roster%box_pad = 2.0d0
      roster%a_lat = 0.0d0
   end subroutine reset_derived

   ! window_reconcile() - the A_LAT/BOXPAD paradigm condition at assembly time:
   ! molecular paradigm: neither key may appear; surface paradigm: both required,
   ! staged values pulled into the roster (the single write site of the window
   ! scalars). The staged-value pulls are named interceptions (key + source line +
   ! the quoted raw string on an unreadable value); unreadable values loud-stop in
   ! read_input's kind check first, guarded anyway per the fail-safe stance
   subroutine window_reconcile()
      integer :: k, ios
      real(8) :: v
      character(len=512) :: sval    ! staging-value buffer (an internal read needs a
                                    ! variable, not an accessor call, as its unit)
      character(len=640) :: msg     ! loud-stop message buffer (key + line + quoted value)
      if (reactants%surface_model == 0) then
         if (user_row('A_LAT') > 0) then
            call stop_assembly('key A_LAT is a surface-paradigm key, but this system derived '// &
                               'the molecular paradigm (surface_model = 0) - the molecular '// &
                               'paradigm carries no surface lattice constant')
         end if
         if (user_row('BOXPAD') > 0) then
            call stop_assembly('key BOXPAD is a surface-paradigm key, but this system derived '// &
                               'the molecular paradigm (surface_model = 0) - the molecular '// &
                               'paradigm has no periodic-box fill half-width')
         end if
         return
      end if
      k = user_row('A_LAT')
      if (k == 0) then
         call stop_assembly('surface paradigm (surface_model = 1) is missing the required key '// &
                            'A_LAT - the surface lattice constant has no generic default; '// &
                            'supply it in the input')
      end if
      sval = staging_val(k)
      read (sval, *, iostat=ios) v
      if (ios /= 0) then
         write (msg, '(a,i0,a,a,a)') 'key A_LAT (line ', staging_line(k), &
            '): bad real value "', trim(sval), '"'
         call stop_assembly(trim(msg))
      end if
      roster%a_lat = v
      k = user_row('BOXPAD')
      if (k == 0) then
         call stop_assembly('surface paradigm (surface_model = 1) is missing the required key '// &
                            'BOXPAD - the periodic-box fill half-width is required alongside '// &
                            'A_LAT; supply it in the input')
      end if
      sval = staging_val(k)
      read (sval, *, iostat=ios) v
      if (ios /= 0) then
         write (msg, '(a,i0,a,a,a)') 'key BOXPAD (line ', staging_line(k), &
            '): bad real value "', trim(sval), '"'
         call stop_assembly(trim(msg))
      end if
      roster%box_pad = v
   end subroutine window_reconcile

   ! user_row(key) - staging row of a user-written input key (0 = absent); the input
   !                  layer never injects a default for the window keys, so any
   !                  staged row is a user row
   integer function user_row(key)
      character(len=*), intent(in) :: key
      integer :: i
      user_row = 0
      do i = 1, staging_n_rows()
         if (trim(staging_key(i)) == trim(key)) then
            user_row = i
            return
         end if
      end do
   end function user_row

   ! stop_assembly(msg) - the loud-stop channel: print the named message, STOP 1
   subroutine stop_assembly(msg)
      character(len=*), intent(in) :: msg
      write (0, '(a)') 'roster_load: assembly error: '//trim(msg)
      write (0, '(a)') 'roster_load: loud-stop (the roster is unusable)'
      stop 1
   end subroutine stop_assembly

end module config_roster
