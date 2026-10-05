!=====================================================================
! config_atoms.f90 - system-composition table (list_atoms): per-atom symbol/mass,
!                     fragment table, equilibrium geometry, unit cell
! Design:
!   list_atoms_load pulls one-way from the sysdef interface buffer (sysdef_load must run
!   first - assembly-order error otherwise): per-atom transfer, contiguous fragment
!   runs, q0 packing, unit-cell derivation (law of cosines, radians), validation,
!   then the A_LAT/BOXPAD window reconcile - their single write site (molecular
!   paradigm forbids both keys, surface paradigm requires both) - and the
!   PROJECTILE/TARGET role reconcile (its single write site: gas pairs require
!   both keys naming distinct fragments, surface requires PROJECTILE alone,
!   single-fragment gas rejects either).
!   a_lat is the input-owned window scalar; all other members are derived.
!=====================================================================
module config_atoms
   use config, only: reactants
   use sysdef_interface, only: buffer_loaded, buffer_frag_count, buffer_frag_nat, &
                            buffer_frag_name, buffer_atom_symbols, buffer_atom_masses, &
                            buffer_atom_geometry, buffer_cell_present, buffer_unit_cell
   use input, only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   implicit none
   private
   public :: list_atoms_t, frag_t, list_atoms, list_atoms_load, list_atoms_clear

   ! Fragment table symmetric and arrayed: one homogeneous frag_t row per
   ! fragment - no A/B asymmetry, no 2-fragment cap
   type :: frag_t
      integer :: nat = 0                   ! fragment atom count [count]
      character(len=64) :: name = ''       ! fragment name (the system-definition file stem,
                                           ! transferred from the sysdef buffer; the
                                           ! PROJECTILE/TARGET role keys name fragments
                                           ! through it) [-]
      real(8) :: mass = 0.0d0              ! fragment total mass [amu] (= sum of member masses)
      integer, allocatable :: list(:)      ! fragment atom index list [global atom index]
      real(8), allocatable :: qz(:)        ! fragment equilibrium geometry (principal-axis coordinates)
                                           ! [Å], length 3*nat
   end type

   type :: list_atoms_t
      real(8), allocatable :: mass(:)          ! per-atom masses [amu] (rigid-surface atoms fixed
                                               ! with a huge mass ~1d30)
      character(len=4), allocatable :: symb(:) ! per-atom element symbols [-] (arrive from the
                                               ! system-definition files via the sysdef interface
                                               ! products)
      type(frag_t), allocatable :: frag(:)     ! fragment table - one homogeneous row per fragment;
                                               ! the surface rows' extra content is still blank
      real(8), allocatable :: q0(:)            ! whole-system equilibrium geometry [Å] (assembled
                                                ! from the fragment geometries by atom index)
      real(8) :: skew = 0.0d0                  ! unit-cell skew angle [rad] (derived from the
                                                ! buffered cell vectors, in radians)
      real(8) :: box_lx = 0.0d0                ! unit-cell x edge length [Å]
      real(8) :: box_ly = 0.0d0                ! unit-cell y edge length [Å]
      real(8) :: a_lat = 0.0d0                ! surface lattice constant [Å] (system parameter with
                                               ! no generic default: surface paradigms must supply
                                               ! the A_LAT input key - the list_atoms_load window
                                               ! reconcile intercepts a missing value and is its
                                               ! only write site)
      integer :: i_proj = 0                   ! projectile fragment row [-] (the role reconcile's
                                               ! single write: the PROJECTILE key names
                                               ! frag(i_proj); 0 = no collision roles in this run)
      integer :: i_targ = 0                   ! target fragment row [-] (the TARGET key, same
                                               ! write site; 0 = none - surface paradigm and
                                               ! single-fragment gas carry no target row)
      ! Left blank (content, not form): what the surface-paradigm fragment rows carry
      ! beyond the homogeneous shape - finalized together with the surface members /
      ! incident channel
   end type

   type(list_atoms_t) :: list_atoms    ! list_atoms instance (same instance pattern as the config module;
                               ! assembled by list_atoms_load; state dimensions =
                               ! 3*size(list_atoms%mass) - no cached natoms/ndof fields,
                               ! size-derived)
contains
   !------------------------------------------------------------------
   ! list_atoms_load() - assemble and validate the atom list from the sysdef buffer
   !                 products (one-way pull; sysdef_load must run first -
   !                 assembly-order error otherwise); flow steps numbered below
   !------------------------------------------------------------------
   subroutine list_atoms_load()
      character(len=4), allocatable :: symb(:)
      real(8), allocatable :: amass(:), geom(:,:)
      real(8) :: cell(3,3), box_lx, box_ly, cos_skew
      character(len=256) :: msg
      integer :: n_atoms, n_frag, i, j, k, offset

      ! 0. ordering-rule guard (assembly-order error)
      if (.not. buffer_loaded()) then
         call stop_assembly('no buffered sysdef products - sysdef_load must run first '// &
                            '(assembly-order error, the config_atoms ordering statute)')
      end if

      ! 1. per-atom transfer (fresh copies from the buffer getters)
      symb  = buffer_atom_symbols()
      amass = buffer_atom_masses()
      geom  = buffer_atom_geometry()
      n_atoms = size(symb)
      if (n_atoms <= 0 .or. size(amass) /= n_atoms .or. &
          size(geom,1) /= 3 .or. size(geom,2) /= n_atoms) then
         call stop_assembly('buffered per-atom products inconsistent (empty or mismatched lengths)')
      end if
      n_frag = buffer_frag_count()
      if (n_frag <= 0) call stop_assembly('buffered fragment count is zero')

      call reset_derived()
      allocate (list_atoms%symb(n_atoms), list_atoms%mass(n_atoms), list_atoms%q0(3*n_atoms))
      allocate (list_atoms%frag(n_frag))
      list_atoms%symb = symb
      list_atoms%mass = amass

      ! 2. fragment tables: one contiguous run per fragment
      offset = 0
      do i = 1, n_frag
         associate (fr => list_atoms%frag(i))
            fr%nat = buffer_frag_nat(i)
            fr%name = buffer_frag_name(i)
            if (fr%nat < 0) then
               write (msg, '(a,i0,a,i0)') 'buffered fragment ', i, ' has a negative atom count ', fr%nat
               call stop_assembly(trim(msg))
            end if
            ! bounds guard BEFORE any allocate/geom access: this fragment's contiguous
            ! run must fit inside the buffered atom count (a too-large nat aborts with a named error
            ! named here - never an out-of-bounds geom read; fail-safe defensive order)
            if (offset + fr%nat > n_atoms) then
               write (msg, '(a,i0,a,i0,a,i0,a,i0)') 'fragment ', i, ' atom count ', fr%nat, &
                  ' exceeds the remaining buffered atoms (', n_atoms - offset, ' left)'
               call stop_assembly(trim(msg))
            end if
            allocate (fr%list(fr%nat), fr%qz(3*fr%nat))
            do j = 1, fr%nat
               k = offset + j
               fr%list(j) = k
               fr%qz(3*j-2:3*j) = geom(:,k)
            end do
            fr%mass = sum(list_atoms%mass(fr%list))
            offset = offset + fr%nat
         end associate
      end do
      if (offset /= n_atoms) then
         write (msg, '(a,i0,a,i0)') 'buffered fragment atom counts sum to ', offset, &
            ' but the buffered atom count is ', n_atoms
         call stop_assembly(trim(msg))
      end if

      ! 3. q0 assembly: the whole-system equilibrium geometry packed by atom index
      do k = 1, n_atoms
         list_atoms%q0(3*k-2:3*k) = geom(:,k)
      end do

      ! 4. unit cell (periodic-box route: law of cosines on the buffered cell vectors)
      if (buffer_cell_present()) then
         cell = buffer_unit_cell()
         box_lx = sqrt(dot_product(cell(1,:), cell(1,:)))
         box_ly = sqrt(dot_product(cell(2,:), cell(2,:)))
         if (box_lx <= 0.0d0 .or. box_ly <= 0.0d0) then
            write (msg, '(a,es12.4,a,es12.4,a)') 'degenerate unit cell: box_lx=', box_lx, &
               ' box_ly=', box_ly, ' - zero-length in-plane cell vector in the buffered unit cell'
            call stop_assembly(trim(msg))
         end if
         cos_skew = dot_product(cell(1,:), cell(2,:))/(box_lx*box_ly)
         cos_skew = max(-1.0d0, min(1.0d0, cos_skew))
         list_atoms%box_lx = box_lx
         list_atoms%box_ly = box_ly
         list_atoms%skew = acos(cos_skew)
      end if

      ! 5. validation on the assembled table
      do i = 1, size(list_atoms%frag)
         if (list_atoms%frag(i)%nat > n_atoms) then
            write (msg, '(a,i0,a,i0,a,i0)') 'fragment ', i, ' atom count ', list_atoms%frag(i)%nat, &
               ' exceeds the system atom count ', n_atoms
            call stop_assembly(trim(msg))
         end if
         do j = 1, list_atoms%frag(i)%nat
            k = list_atoms%frag(i)%list(j)
            if (k < 1 .or. k > n_atoms) then
               write (msg, '(a,i0,a,i0,a,i0)') 'fragment ', i, ' list entry ', j, &
                  ' is out of range: ', k
               call stop_assembly(trim(msg))
            end if
         end do
      end do
      do k = 1, n_atoms
         if (list_atoms%mass(k) <= 0.0d0) then
            write (msg, '(a,i0,a,es12.4)') 'atom ', k, ' carries a non-positive mass ', list_atoms%mass(k)
            call stop_assembly(trim(msg))
         end if
      end do

      ! 6. surface-window reconcile (the paradigm condition of A_LAT/BOXPAD - the
      !    single write site of the window scalars)
      call window_reconcile()

      ! 7. role reconcile (the paradigm condition of PROJECTILE/TARGET - the single
      !    write site of the collision-role rows)
      call role_reconcile()
   end subroutine list_atoms_load

   !------------------------------------------------------------------
   ! list_atoms_clear() - deallocate all allocatable members and reset the scalars
   !                   to their generic defaults
   !------------------------------------------------------------------
   subroutine list_atoms_clear()
      if (allocated(list_atoms%mass)) deallocate (list_atoms%mass)
      if (allocated(list_atoms%symb)) deallocate (list_atoms%symb)
      if (allocated(list_atoms%frag)) deallocate (list_atoms%frag)
      if (allocated(list_atoms%q0)) deallocate (list_atoms%q0)
      list_atoms%skew = 0.0d0
      list_atoms%box_lx = 0.0d0
      list_atoms%box_ly = 0.0d0
      list_atoms%a_lat = 0.0d0
      list_atoms%i_proj = 0
      list_atoms%i_targ = 0
   end subroutine list_atoms_clear

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! reset_derived() - drop the DERIVED members ahead of a re-fill; the window scalars
   !                    reset to generic defaults too, so a re-load cannot inherit the
   !                    previous load's window values (step 6 re-pulls from the buffer)
   subroutine reset_derived()
      if (allocated(list_atoms%mass)) deallocate (list_atoms%mass)
      if (allocated(list_atoms%symb)) deallocate (list_atoms%symb)
      if (allocated(list_atoms%frag)) deallocate (list_atoms%frag)
      if (allocated(list_atoms%q0)) deallocate (list_atoms%q0)
      list_atoms%skew = 0.0d0
      list_atoms%box_lx = 0.0d0
      list_atoms%box_ly = 0.0d0
      list_atoms%a_lat = 0.0d0
      list_atoms%i_proj = 0
      list_atoms%i_targ = 0
   end subroutine reset_derived

   ! window_reconcile() - the A_LAT paradigm condition at assembly time:
   ! molecular paradigm: the key may not appear; surface paradigm: required,
   ! buffered values pulled into the atom list (the single write site of the window
   ! scalars). The buffered-value pulls are named interceptions (key + source line +
   ! the quoted raw string on an unreadable value); unreadable values abort with a named error in
   ! read_input's kind check first, guarded anyway per the fail-safe stance
   subroutine window_reconcile()
      integer :: k, ios
      real(8) :: v
      character(len=512) :: sval    ! buffer-value buffer (an internal read needs a
                                    ! variable, not an accessor call, as its unit)
      character(len=640) :: msg     ! named abort message buffer (key + line + quoted value)
      if (reactants%surface_model == 0) then
         if (user_row('A_LAT') > 0) then
            call stop_assembly('key A_LAT is a surface-paradigm key, but this system derived '// &
                               'the molecular paradigm (surface_model = 0) - the molecular '// &
                               'paradigm carries no surface lattice constant')
         end if
         return
      end if
      k = user_row('A_LAT')
      if (k == 0) then
         call stop_assembly('surface paradigm (surface_model = 1) is missing the required key '// &
                            'A_LAT - the surface lattice constant has no generic default; '// &
                            'supply it in the input')
      end if
      sval = buffer_val(k)
      read (sval, *, iostat=ios) v
      if (ios /= 0) then
         write (msg, '(a,i0,a,a,a)') 'key A_LAT (line ', buffer_line(k), &
            '): bad real value "', trim(sval), '"'
         call stop_assembly(trim(msg))
      end if
      list_atoms%a_lat = v
   end subroutine window_reconcile

   ! role_reconcile() - the PROJECTILE/TARGET paradigm condition at assembly time (the
   !                    single write site of the role rows). Molecular paradigm: one
   !                    fragment carries no collision roles and rejects either key; two
   !                    or more fragments REQUIRE both keys, naming distinct fragments
   !                    (no positional default exists). Surface paradigm: PROJECTILE
   !                    required (it names the one active fragment); TARGET is a
   !                    gas-pair key and rejects. Names match the fragment table
   !                    case-insensitively; a miss or an ambiguous fold aborts naming
   !                    every available fragment name
   subroutine role_reconcile()
      character(len=512) :: sval
      character(len=800) :: msg
      integer :: k, nf, ip, it
      nf = size(list_atoms%frag)
      list_atoms%i_proj = 0
      list_atoms%i_targ = 0
      if (reactants%surface_model == 0) then
         if (nf == 1) then
            if (user_row('PROJECTILE') > 0 .or. user_row('TARGET') > 0) then
               call stop_assembly('a single-fragment gas-phase system carries no collision '// &
                  'roles - keys PROJECTILE/TARGET name the fragments of a collision pair '// &
                  'and are not applicable here')
            end if
            return
         end if
         k = user_row('PROJECTILE')
         if (k == 0) then
            write (msg, '(a,i0,a)') 'molecular paradigm with ', nf, ' fragments is missing '// &
               'the required key PROJECTILE - the projectile role has no positional '// &
               'default; name it after one of the fragments: '//frag_names()
            call stop_assembly(trim(msg))
         end if
         sval = buffer_val(k)
         ip = name_row('PROJECTILE', sval, buffer_line(k))
         k = user_row('TARGET')
         if (k == 0) then
            write (msg, '(a,i0,a)') 'molecular paradigm with ', nf, ' fragments is missing '// &
               'the required key TARGET - the target role has no positional default; '// &
               'name it after one of the fragments: '//frag_names()
            call stop_assembly(trim(msg))
         end if
         sval = buffer_val(k)
         it = name_row('TARGET', sval, buffer_line(k))
         if (ip == it) then
            call stop_assembly('keys PROJECTILE and TARGET both name fragment "'// &
               trim(list_atoms%frag(ip)%name)//'" - the projectile and the target must be '// &
               'two distinct fragments; available names: '//frag_names())
         end if
         list_atoms%i_proj = ip
         list_atoms%i_targ = it
      else
         if (user_row('TARGET') > 0) then
            call stop_assembly('key TARGET is a gas-pair role key, but this system derived '// &
               'the surface paradigm - the surface is not an atom list fragment and the '// &
               'target role does not exist here')
         end if
         k = user_row('PROJECTILE')
         if (k == 0) then
            call stop_assembly('surface paradigm (surface_model = 1) is missing the required '// &
               'key PROJECTILE - it names the one active (incident) fragment; available '// &
               'names: '//frag_names())
         end if
         sval = buffer_val(k)
         list_atoms%i_proj = name_row('PROJECTILE', sval, buffer_line(k))
      end if
   end subroutine role_reconcile

   ! name_row(key, nm, ln) - list_atoms fragment row whose name matches nm (case-insensitive
   !                         trim fold; 0 is never returned on a miss - the miss is the
   !                         named abort below); an ambiguous fold (two fragments sharing
   !                         a case-insensitive name) likewise aborts
   integer function name_row(key, nm, ln)
      character(len=*), intent(in) :: key, nm
      integer, intent(in) :: ln
      character(len=800) :: msg
      integer :: i, n_hit
      name_row = 0
      n_hit = 0
      do i = 1, size(list_atoms%frag)
         if (ci_eq(list_atoms%frag(i)%name, nm)) then
            n_hit = n_hit + 1
            name_row = i
         end if
      end do
      if (n_hit == 0) then
         write (msg, '(a,i0,a)') 'key '//trim(key)//': no list_atoms fragment named "'// &
            trim(nm)//'" (line ', ln, '); available names: '//frag_names()
         call stop_assembly(trim(msg))
      else if (n_hit > 1) then
         call stop_assembly('key '//trim(key)//': fragment name "'//trim(nm)//'" matches '// &
            'more than one list_atoms fragment (case-insensitive fold) - the role keys '// &
            'cannot name fragments unambiguously; available names: '//frag_names())
      end if
   end function name_row

   ! ci_eq(a, b) - case-insensitive equality of two trimmed strings
   pure logical function ci_eq(a, b)
      character(len=*), intent(in) :: a, b
      integer :: i, ca, cb
      ci_eq = len_trim(a) == len_trim(b)
      if (.not. ci_eq) return
      do i = 1, len_trim(a)
         ca = iachar(a(i:i))
         cb = iachar(b(i:i))
         if (ca >= iachar('a') .and. ca <= iachar('z')) ca = ca - 32
         if (cb >= iachar('a') .and. cb <= iachar('z')) cb = cb - 32
         if (ca /= cb) then
            ci_eq = .false.
            return
         end if
      end do
   end function ci_eq

   ! frag_names() - the atom list fragment names as a comma list (abort-message payload)
   function frag_names() result(s)
      character(len=512) :: s
      integer :: i
      s = ''
      do i = 1, size(list_atoms%frag)
         if (i > 1) s = trim(s)//', '
         s = trim(s)//trim(list_atoms%frag(i)%name)
      end do
   end function frag_names

   ! user_row(key) - buffer row of a user-written input key (0 = absent); the input
   !                  layer never injects a default for the window keys, so any
   !                  buffered row is a user row
   integer function user_row(key)
      character(len=*), intent(in) :: key
      integer :: i
      user_row = 0
      do i = 1, buffer_n_rows()
         if (trim(buffer_key(i)) == trim(key)) then
            user_row = i
            return
         end if
      end do
   end function user_row

   ! stop_assembly(msg) - the named-abort channel: print the named message, STOP 1
   subroutine stop_assembly(msg)
      character(len=*), intent(in) :: msg
      write (0, '(a)') 'list_atoms_load: assembly error: '//trim(msg)
      write (0, '(a)') 'list_atoms_load: fatal (the atom list is unusable)'
      stop 1
   end subroutine stop_assembly

end module config_atoms
