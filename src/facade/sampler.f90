!=====================================================================
! sampler.f90 - initial-state distribution sampling framework (facade):
!               distribution-member registration + per-trajectory
!               orchestration - the per-fragment members write the
!               initial-state components they own, the incident-channel
!               member assembles the collision after them
! Contract:
!   samp_reg(name, proc) - register one distribution member row (assembly
!     phase; a member = one methods/ file + one registration line).
!     Duplicate name or full table -> loud-stop naming it.
!   samp_run(scheme, sta) - prepare one trajectory's initial state in one
!     call, in this order:
!     1. assembly guards: the roster fragment table must be assembled
!        (roster_load runs first) and the state sized to the roster
!        (size(q) = size(p) = 3*size(roster%mass)); violation -> loud-stop
!        naming the mismatch.
!     2. count reconcile BEFORE any dispatch: the per-fragment scheme array
!        reactants%dist_scheme must cover the roster exactly (its length
!        equals the roster fragment count) - mismatch -> loud-stop naming
!        BOTH values.
!     3. scheme lookup: the scheme argument names the incident-channel
!        member - the one orchestration member NOT selected by a
!        dist_scheme row (the collision-assembly channel has no config
!        field of its own, so the driver passes its member name, today
!        'incident'); unregistered -> loud-stop naming it.
!     4. per-fragment dispatch in roster fragment order: fragment i's member
!        name comes from the scheme-code table below (the integer codes of
!        reactants%dist_scheme denote member names - the same words the
!        input grammar accepts for DIST_SCHEME); an out-of-range code or an
!        unregistered name -> loud-stop naming it. Only sta is passed: the
!        fragment index and the distribution parameters are member-init
!        knowledge (member = code + its own data; sta is the one runtime
!        channel), and each member writes only the q/p/s components it
!        owns.
!   samp_code(name) - the scheme-code table row of a member word (0 =
!     absent); the single source a member init consults for its fragment
!     binding (which dist_scheme rows select it).
!     5. incident channel ONCE, AFTER all fragment members - dispatched
!        unless the system is single-fragment gas-phase (one fragment and
!        surface_model = 0: a lone gas-phase molecule has no collision
!        assembly, the relative-motion constructs would degenerate).
!   Layering statute: the sampler family carries no use of sysdef_facade -
!   the roster is the one and only geometry/composition channel into this
!   family (this module reads roster%frag / roster%mass and nothing else of
!   the system's composition).
! Units: none internal (distribution shapes and their units are member-side
! knowledge); table sizes in rows (capacity max_dist = 16).
!=====================================================================
module sampler
   use state, only: state_t
   use config, only: reactants
   use config_roster, only: roster
   implicit none
   private
   public :: samp_reg, samp_run, samp_code

   ! Distribution member interface: members sample (q0,p0,s0) from their own
   ! distribution and write the components they own into the physical state
   abstract interface
      subroutine dist_proc_i(sta)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! initial state (the q/p/s components owned by the member are
         ! overwritten by sampling)
      end subroutine dist_proc_i
   end interface

   ! Distribution member table (the registry family - a new member = a new
   ! methods/ file + one registration line)
   integer, parameter :: max_dist = 16   ! member-table capacity [row]
   type :: dist_entry_t
      character(len=32) :: name = ''     ! member name (the dispatch key)
      procedure(dist_proc_i), pointer, nopass :: proc => null() ! member sampling routine
   end type dist_entry_t
   type(dist_entry_t) :: dist_tbl(max_dist) ! distribution member table (module-private)
   integer :: n_dist = 0                   ! registered row count [row]

   ! Scheme-code table: the integer codes of reactants%dist_scheme denote
   ! member names (the same words the input grammar accepts for DIST_SCHEME);
   ! a new member adds its grammar word together with its row here, next to
   ! its registration line. The incident channel is NOT a per-fragment
   ! scheme word: it is selected by the samp_run scheme argument, so it
   ! carries no row here (its registration row in the member table above
   ! stays - the scheme-argument lookup path uses it)
   integer, parameter :: n_scheme_code = 7
   character(len=32), parameter :: scheme_code_name(n_scheme_code) = &
      [ character(len=32) :: 'boltzmann', 'ebk', 'rotation', 'j', &
                              'surface_oscillator', 'normalmode', 'barrier_excitation' ]
contains
   !------------------------------------------------------------------
   ! samp_code(name) - scheme-code table row of a member word (0 = absent):
   !                   the single source a member init consults to find which
   !                   dist_scheme rows select it (fragment binding)
   pure integer function samp_code(name)
      character(len=*), intent(in) :: name  ! member word (e.g. 'boltzmann')
      integer :: i
      samp_code = 0
      do i = 1, n_scheme_code
         if (trim(scheme_code_name(i)) == trim(name)) then
            samp_code = i
            exit
         end if
      end do
   end function samp_code

   !------------------------------------------------------------------
   ! samp_reg(name, proc) - register a distribution member into the member table
   subroutine samp_reg(name, proc)
      character(len=*), intent(in) :: name  ! member name (e.g. 'boltzmann'/'ebk'/'incident')
      procedure(dist_proc_i) :: proc        ! member sampling routine
      if (find_row(name) /= 0) then
         call stop_samp('samp_reg', 'distribution member "'//trim(name)// &
            '" is already registered (duplicate registration line)')
      end if
      if (n_dist >= max_dist) then
         call stop_samp('samp_reg', 'the distribution member table is full ('// &
            trim(i8_str(max_dist))//' rows) - member "'//trim(name)//'" not registered')
      end if
      n_dist = n_dist + 1
      dist_tbl(n_dist)%name = name
      dist_tbl(n_dist)%proc => proc
   end subroutine samp_reg

   !------------------------------------------------------------------
   ! samp_run(scheme, sta) - dispatch sampling: call the members to write the
   ! initial-state components (one call prepares one whole trajectory)
   subroutine samp_run(scheme, sta)
      character(len=*), intent(in) :: scheme ! incident-channel member name (the one
                                             ! orchestration member not selected by a
                                             ! dist_scheme row; the driver passes the
                                             ! incident member's name)
      type(state_t), intent(inout) :: sta    ! initial state (q0/p0/s0 - overwritten by
                                             ! the dispatched members, each in the
                                             ! components it owns)
      integer :: n_frag, nat, irow, frow, i

      ! 1. assembly guards (ordering statute + state sizing)
      if (.not. allocated(roster%frag)) then
         call stop_samp('samp_run', 'no roster fragment table - roster_load must run '// &
                        'first (assembly-order error)')
      end if
      n_frag = size(roster%frag)
      nat = 0
      if (allocated(roster%mass)) nat = size(roster%mass)
      if (size(sta%q) /= 3*nat .or. size(sta%p) /= 3*nat) then
         call stop_samp('samp_run', 'state size mismatch: size(sta%q) = '// &
            trim(i8_str(size(sta%q)))//' dof, size(sta%p) = '//trim(i8_str(size(sta%p)))// &
            ' dof but the roster carries '//trim(i8_str(nat))//' atoms ('//trim(i8_str(3*nat))// &
            ' dof) - state_create must follow the roster')
      end if

      ! 2. count reconcile BEFORE any dispatch
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_samp('samp_run', 'no per-fragment sampling scheme (reactants%'// &
                        'dist_scheme unallocated - the DIST_SCHEME provision is missing)')
      end if
      if (size(reactants%dist_scheme) /= n_frag) then
         call stop_samp('samp_run', 'dist_scheme length '// &
            trim(i8_str(size(reactants%dist_scheme)))//' does not match the roster fragment '// &
            'count '//trim(i8_str(n_frag))//' - the per-fragment schemes must cover the '// &
            'fragments exactly')
      end if

      ! 3. incident-channel member lookup (fail fast, before any dispatch)
      irow = find_row(scheme)
      if (irow == 0) then
         call stop_samp('samp_run', 'sampling scheme "'//trim(scheme)// &
                        '" is not a registered distribution member')
      end if

      ! 4. per-fragment dispatch in roster fragment order (sta is the one
      !    runtime channel; the fragment binding is member-init knowledge)
      do i = 1, n_frag
         if (reactants%dist_scheme(i) < 1 .or. reactants%dist_scheme(i) > n_scheme_code) then
            call stop_samp('samp_run', 'fragment '//trim(i8_str(i))//' carries scheme code '// &
               trim(i8_str(reactants%dist_scheme(i)))//' outside the scheme-code table (1..'// &
               trim(i8_str(n_scheme_code))//')')
         end if
         frow = find_row(scheme_code_name(reactants%dist_scheme(i)))
         if (frow == 0) then
            call stop_samp('samp_run', 'fragment '//trim(i8_str(i))//' scheme member "'// &
               trim(scheme_code_name(reactants%dist_scheme(i)))// &
               '" is not registered (incomplete member library assembly)')
         end if
         call dist_tbl(frow)%proc(sta)
      end do

      ! 5. incident channel once, after the fragment members (skipped only
      !    for the single-fragment gas-phase system)
      if (n_frag > 1 .or. reactants%surface_model /= 0) then
         call dist_tbl(irow)%proc(sta)
      end if
   end subroutine samp_run

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! find_row(name) - member-table row of name (0 = absent)
   function find_row(name) result(row)
      character(len=*), intent(in) :: name
      integer :: row, i
      row = 0
      do i = 1, n_dist
         if (trim(dist_tbl(i)%name) == trim(name)) then
            row = i
            return
         end if
      end do
   end function find_row

   ! i8_str(v) - integer to decimal string for named loud-stop messages
   ! (oversized local, fixed width - callers trim() when splicing into text)
   function i8_str(v) result(s)
      integer, intent(in) :: v
      character(len=32) :: s
      write (s, '(i0)') v
   end function i8_str

   ! stop_samp(where, msg) - the loud-stop channel (same semantics as the
   ! sibling facades: name the caller, print the named message, STOP 1)
   subroutine stop_samp(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': loud-stop (the initial state is not sampled)'
      stop 1
   end subroutine stop_samp

end module sampler
