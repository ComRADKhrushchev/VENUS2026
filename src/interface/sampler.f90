!=====================================================================
! sampler.f90 - initial-state distribution sampling interface:
!   distribution-member registration + per-trajectory two-phase dispatch
!   (draw phase -> realize phase)
! Design:
!   samp_reg(name, draw, realize) registers a member row (a new member = a
!   new methods/ file + one registration line). Every member splits into a
!   DRAW procedure (consumes the RNG, builds its module-private sample
!   record, touches no state) and a REALIZE procedure (writes q/p from the
!   record, draws no RNG). samp_run(scheme, sta) prepares one trajectory's
!   initial state: DRAW PHASE - each DISTINCT per-fragment member named by
!   reactants%dist_scheme (list_atoms order, dedup), then the orientation
!   stage's draws (the RANDOM_ORIENT option; one rot_rand_rmat per eligible
!   fragment - in the surface paradigm only the projectile fragment is
!   eligible: a substrate fragment's crystal frame must not rotate), then
!   the incident member's draw; REALIZE PHASE - the member realizations in
!   the same list_atoms order (the historical fused dispatch order - a
!   compatibility statute: reordering would change what the equilibration
!   members' runs see), then the orientation stage applies its drawn
!   matrices, then the incident member assembles the collision. The
!   incident channel splits per paradigm (incident_surface /
!   incident_pair / incident_single; samp_incident_word resolves which)
!   and is dispatched UNIFORMLY - the trivial single-reactant member is a
!   bitwise no-op (the historical skip-guard produced the identical
!   stream). Because realizations draw no RNG, the single draw phase
!   reproduces the historical fused stream exactly (both RANDOM_ORIENT
!   settings). The list_atoms is the only composition channel into this
!   family; units are member-side.
!=====================================================================
module sampler
   use state, only: state_t
   use config, only: reactants
   use config_atoms, only: list_atoms
   use geometry, only: rot_rand_rmat, rot_apply
   use bath, only: bath_run
   implicit none
   private
   public :: samp_reg, samp_run, samp_code, samp_n_dist, samp_incident_word

   ! Distribution member interface: the DRAW half consumes the RNG and builds
   ! the member-private sample record (no state access); the REALIZE half
   ! writes the owned q/p components from the record (no RNG draws)
   abstract interface
      subroutine dist_draw_i()
      end subroutine dist_draw_i
      subroutine dist_realize_i(sta)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! initial state (the q/p components owned by the
                                              ! member are overwritten from the drawn record)
      end subroutine dist_realize_i
   end interface

   ! Distribution member table (the registry family - a new member = a new
   ! methods/ file + one registration line)
   integer, parameter :: max_dist = 16   ! member-table capacity [row]
   type :: dist_entry_t
      character(len=32) :: name = ''     ! member name (the dispatch key)
      procedure(dist_draw_i), pointer, nopass :: draw => null()       ! draw routine
      procedure(dist_realize_i), pointer, nopass :: realize => null() ! realize routine
   end type dist_entry_t
   type(dist_entry_t) :: dist_tbl(max_dist) ! distribution member table (module-private)
   integer :: n_dist = 0                   ! registered row count [row]

   ! Scheme-code table: the integer codes of reactants%dist_scheme denote
   ! member names (the same words the input grammar accepts for DIST_SCHEME);
   ! a new member adds its grammar word together with its row here, next to
   ! its registration line. The incident channel is NOT a per-fragment
   ! scheme word: it splits per paradigm into incident_surface /
   ! incident_pair / incident_single (samp_incident_word resolves which),
   ! so those rows carry no scheme code here (their registration rows in
   ! the member table above stay - the scheme-argument lookup path uses them)
   integer, parameter :: n_scheme_code = 10
   character(len=32), parameter :: scheme_code_name(n_scheme_code) = &
      [ character(len=32) :: 'boltzmann', 'ebk', 'rotation', 'j', &
                              'surface_oscillator', 'normalmode', 'barrier_excitation', &
                              'glo_target', 'thermalize', 'wigner' ]

   integer, parameter :: max_frag = 64    ! fragment-count working ceiling [fragment]
                                          ! (the orientation-stage draw/apply arrays are
                                          ! sized by it; rosters stay far below)
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
   ! samp_n_dist() - registered member count (read-only public probe; the
   !                 assembly-complete summary's 'members=' witness)
   pure integer function samp_n_dist()
      samp_n_dist = n_dist
   end function samp_n_dist

   !------------------------------------------------------------------
   ! samp_reg(name, draw, realize) - register a distribution member into the
   !                                    member table
   subroutine samp_reg(name, draw, realize)
      character(len=*), intent(in) :: name ! member name (e.g. 'boltzmann'/'ebk'/'incident_surface')
      procedure(dist_draw_i) :: draw       ! member draw routine (RNG -> sample record)
      procedure(dist_realize_i) :: realize ! member realize routine (record -> q/p)
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
      dist_tbl(n_dist)%draw => draw
      dist_tbl(n_dist)%realize => realize
   end subroutine samp_reg

   !------------------------------------------------------------------
   ! samp_incident_word() - the paradigm-resolved incident member's name
   !                 (the driver passes it to samp_run; the regs register
   !                 exactly this one member per process - an unregistered
   !                 word fails loudly at the samp_run lookup)
   !------------------------------------------------------------------
   function samp_incident_word() result(word)
      character(len=32) :: word
      integer :: nf
      if (reactants%surface_model /= 0) then
         word = 'incident_surface'
      else
         nf = 0
         if (allocated(list_atoms%frag)) nf = size(list_atoms%frag)
         if (nf > 1) then
            word = 'incident_pair'
         else
            word = 'incident_single'
         end if
      end if
   end function samp_incident_word

   !------------------------------------------------------------------
   ! samp_run(scheme, sta) - dispatch sampling in two phases: DRAW (every
   !   member's random content first - the whole trajectory's RNG
   !   consumption happens here, roster/list_atoms order, then the
   !   orientation matrices, then the incident prescription), then REALIZE
   !   (state writes only, the same order; no RNG anywhere in this phase)
   subroutine samp_run(scheme, sta)
      character(len=*), intent(in) :: scheme ! incident-channel member name (the one
                                             ! orchestration member not selected by a
                                             ! dist_scheme row; the driver passes
                                             ! samp_incident_word() - the
                                             ! paradigm-resolved member's name)
      type(state_t), intent(inout) :: sta    ! initial state (q0/p0/s0 - overwritten by
                                             ! the dispatched members, each in the
                                             ! components it owns)
      integer :: n_frag, nat, irow, frow, i, n_called
      integer :: called(max_dist)               ! member rows already dispatched this
                                                ! call (dedup: one call per distinct member)
      logical :: orient_on, eligible(max_frag)  ! orientation-stage mask (per fragment)
      real(8) :: orient_r(3, 3, max_frag)       ! the drawn orientation matrices

      ! 1. assembly guards (ordering rule + state sizing)
      if (.not. allocated(list_atoms%frag)) then
         call stop_samp('samp_run', 'no list_atoms fragment table - list_atoms_load must run '// &
                        'first (assembly-order error)')
      end if
      n_frag = size(list_atoms%frag)
      nat = 0
      if (allocated(list_atoms%mass)) nat = size(list_atoms%mass)
      if (size(sta%q) /= 3*nat .or. size(sta%p) /= 3*nat) then
         call stop_samp('samp_run', 'state size mismatch: size(sta%q) = '// &
            trim(i8_str(size(sta%q)))//' dof, size(sta%p) = '//trim(i8_str(size(sta%p)))// &
            ' dof but the atom list carries '//trim(i8_str(nat))//' atoms ('//trim(i8_str(3*nat))// &
            ' dof) - state_create must follow the atom list')
      end if
      if (n_frag > max_frag) then
         call stop_samp('samp_run', 'fragment count '//trim(i8_str(n_frag))// &
            ' exceeds the sampler working capacity - raise max_frag')
      end if

      ! 2. count reconcile BEFORE any dispatch
      if (.not. allocated(reactants%dist_scheme)) then
         call stop_samp('samp_run', 'no per-fragment sampling scheme (reactants%'// &
                        'dist_scheme unallocated - the DIST_SCHEME provision is missing)')
      end if
      if (size(reactants%dist_scheme) /= n_frag) then
         call stop_samp('samp_run', 'dist_scheme length '// &
            trim(i8_str(size(reactants%dist_scheme)))//' does not match the atom list fragment '// &
            'count '//trim(i8_str(n_frag))//' - the per-fragment schemes must cover the '// &
            'fragments exactly')
      end if

      ! 3. incident-channel member lookup (fail fast, before any dispatch)
      irow = find_row(scheme)
      if (irow == 0) then
         call stop_samp('samp_run', 'sampling scheme "'//trim(scheme)// &
            '" is not a registered distribution member')
      end if

      ! 4. orientation-stage eligibility mask (computed once; both phases and
      !    the eligibility statute read the same mask): the RANDOM_ORIENT
      !    option orients every nat>=2 fragment whose scheme word is not
      !    'rotation' (that member already draws its orientation inside its
      !    fixed-J realization - a second draw would be distributionally
      !    neutral but wasteful); in the surface paradigm ONLY the projectile
      !    fragment is eligible - a substrate fragment carries the crystal
      !    frame, which the surface paradigm must not rotate
      orient_on = reactants%random_orient
      do i = 1, n_frag
         eligible(i) = orient_on .and. list_atoms%frag(i)%nat >= 2 .and. &
                       reactants%dist_scheme(i) /= samp_code('rotation') .and. &
                       (reactants%surface_model == 0 .or. i == list_atoms%i_proj)
      end do

      ! 5. DRAW PHASE (all RNG consumption, in the historical order):
      !    per-fragment members (each DISTINCT row once, at its first
      !    list_atoms occurrence - the member serves every fragment that
      !    selected it inside that one call), then the orientation
      !    matrices, then the incident prescription. Realizations draw no
      !    RNG, so this single phase reproduces the fused stream exactly
      n_called = 0
      called = 0
      do i = 1, n_frag
         ! a monoatomic fragment has no internal degrees of freedom - no
         ! distribution member is dispatched for it (its q/p ride the incident
         ! channel alone; the scheme word on it is the input's way of saying
         ! "this fragment carries no internal sampling")
         if (list_atoms%frag(i)%nat < 2) cycle
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
         if (any(called(1:n_called) == frow)) cycle
         n_called = n_called + 1
         called(n_called) = frow
         call dist_tbl(frow)%draw()
      end do
      do i = 1, n_frag
         if (eligible(i)) call rot_rand_rmat(orient_r(:,:,i))
      end do
      call dist_tbl(irow)%draw()

      ! 6. REALIZE PHASE (state writes only, the same list_atoms order - the
      !    historical fused dispatch order; a dynamics-last reordering would
      !    change what the equilibration members' runs see and stays a parked
      !    adjudication), then the orientation stage applies its drawn
      !    matrices, then the incident channel assembles the collision
      n_called = 0
      called = 0
      do i = 1, n_frag
         if (list_atoms%frag(i)%nat < 2) cycle   ! monoatomic: no internal DOFs
         frow = find_row(scheme_code_name(reactants%dist_scheme(i)))
         if (any(called(1:n_called) == frow)) cycle
         n_called = n_called + 1
         called(n_called) = frow
         call dist_tbl(frow)%realize(sta)
      end do

      ! 6b. the pre-evolution bath loop (the BATH provision): N_BATH steps of
      !     [propagation + one collision sweep] over the realized interiors,
      !     ending with one exact whole-system drift removal - BEFORE the
      !     orientation stage and the incident channel (the bath must not
      !     see the beam momenta; the incident placement after it re-places
      !     the COMs). This stage draws RNG INSIDE the loop (the collision
      !     trials depend on the propagated state) - the one sanctioned
      !     exception to the draw/realize split, owned by bath.f90; unarmed
      !     (BATH = NONE) it is a bitwise no-op
      call bath_run(sta)

      do i = 1, n_frag
         if (eligible(i)) then
            call rot_apply(sta, list_atoms%mass, list_atoms%frag(i)%list, orient_r(:,:,i))
         end if
      end do
      call dist_tbl(irow)%realize(sta)
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

   ! i8_str(v) - integer to decimal string for named abort messages
   ! (oversized local, fixed width - callers trim() when splicing into text)
   function i8_str(v) result(s)
      integer, intent(in) :: v
      character(len=32) :: s
      write (s, '(i0)') v
   end function i8_str

   ! stop_samp(where, msg) - the named-abort channel (same semantics as the
   ! sibling interfaces: name the caller, print the named message, STOP 1)
   subroutine stop_samp(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_samp

end module sampler
