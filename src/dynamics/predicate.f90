!=====================================================================
! predicate.f90 - L3 grammar interpreter: termination-predicate expression-tree parsing + evaluation
!                  (leaves: four built-in primitives + the open registered primitive set, E14 track 1;
!                  interior nodes: three boolean combinators, a closed set, E14 track 2 - system
!                  criteria demoted to data)
! Reference: src_VENUS/TEST.f:88-94 (prototype of the hardwired Z_HEIGHT criterion: rd=Q(3) height above the
!       surface + P(3)>0 outward velocity, topmost surface layer at z=0); src_VENUS/TEST.f:106-146
!       (COM_DISTANCE prototype: mass-weighted center-of-mass accumulation :107-130, center-of-mass distance
!       and receding-velocity dot product :133-143, RBAR/RMAX two-threshold test :144-145);
!       docs/final_completion_plan.md §5 N2 (grammar and interface spec - in force as amended by E14:
!       its closed-primitive-set wording is superseded, see the E14 note below); §3 B2
!       (the TEST general criteria registered -> grammar members of this module)
! Input grammar (TERMINATE= expression, a new keyword registered in param_mapping - N2 spec):
!   expression := primitive | AND(e1,e2) | OR(e1,e2) | NOT(e)
!   primitive := COM_DISTANCE(fragmentA, fragmentB, threshold-Ang[, receding-velocity threshold])  <- TEST.f:106-146 prototype
!           | Z_HEIGHT(fragment, threshold-Ang[, OUTWARD])     <- replaces the TEST.f:88-94 hardwired criterion
!           | MAX_STEPS(steps)                                 <- replaces the step-count magic number (always a safety net)
!           | BOUND(fragmentA, fragmentB, lower, upper)        <- bounded window (constrained scattering)
!           | PRIM('name')                                     <- registered system primitive (E14 track 1: 'name' must
!                                                                match a pred_reg_prim row; both-track composition e.g.
!                                                                AND(PRIM('desorb_z'), MAX_STEPS(500000)))
!   fragment := roster fragment name (A/B/...) or ATOM(i)
!
! -- E14 two-track composition (user ruling 2026-09-08; docs/superpowers/specs/2026-09-08-facade-e14-design.md §7) --
!   Track 1 (system material): each blackbox registers its termination primitives via
!   pred_reg_prim(name, proc) - e.g. "adsorbate height above surface" written by the system author;
!   system knowledge lives in the blackbox (the leaf set OPENS). Track 2 (program-side grammar):
!   boolean composition (AND/OR/NOT, closed set) and the MAX_STEPS safety net stay in predicate -
!   composition logic is program knowledge (pure Boolean algebra - the §3 closure argument survives
!   verbatim; the old N2 wording that read the whole grammar as closed is hereby revised).
!   pred_parse: the PRIM leaf matches the registered table at parse time (unknown name -> grammar
!   error naming it); pred_eval: k_prim nodes look up the table and call the member proc with sta.
!=====================================================================
module predicate
   use state,         only: state_t
   use config_roster, only: roster
   implicit none
   private
   public :: frag_ref_t, expr_t
   public :: pred_parse, pred_eval
   public :: com_dist, z_height, bound, steps
   public :: pred_and, pred_or, pred_not
   public :: pred_reg_prim
   public :: k_com_dist, k_z_height, k_bound, k_steps, k_and, k_or, k_not, k_prim

   ! Node-kind codes (values of expr_t%node_kind; E14 two-track: the boolean combinator codes are a
   ! closed set - a new combinator is a grammar-versioning matter; the primitive leaf side is open -
   ! registered system primitives enter as the k_prim leaf without touching this enumeration)
   integer, parameter :: k_com_dist = 1  ! built-in primitive: COM_DISTANCE (fragment-pair center-of-mass distance)
   integer, parameter :: k_z_height = 2  ! built-in primitive: Z_HEIGHT (fragment height above the surface)
   integer, parameter :: k_bound    = 3  ! built-in primitive: BOUND (bounded-window escape)
   integer, parameter :: k_steps    = 4  ! built-in primitive: MAX_STEPS (step-count safety net)
   integer, parameter :: k_and      = 5  ! boolean combinator: AND(e1,e2)
   integer, parameter :: k_or       = 6  ! boolean combinator: OR(e1,e2)
   integer, parameter :: k_not      = 7  ! boolean combinator: NOT(e)
   integer, parameter :: k_prim     = 8  ! registered primitive leaf: PRIM('name') (system termination
                                         ! primitive dispatched through the prim_tbl registry - E14 track 1)

   ! Fragment reference (in-tree representation of the grammar's "fragment": roster fragment name or single atom)
   type :: frag_ref_t
      integer :: i_frag = 0   ! roster fragment number [-] (A/B/... -> 1/2/...; 0 = unset)
      integer :: i_atom = 0   ! single-atom index [-] (ATOM(i) grammar; >0 overrides the fragment reference -
                              ! a single-atom criterion, generalizing the old TEST.f:89 hardwired atom 1)
   end type frag_ref_t

   ! Predicate expression node (recursive tree: primitives are leaves - built-in four plus the
   ! registered PRIM leaf (E14); boolean combinators are interior nodes)
   type :: expr_t
      integer :: node_kind = 0      ! node kind (k_com_dist...k_not, k_prim; 0 = empty / not assembled)
      type(frag_ref_t) :: frag_a   ! fragment reference A (used by com_dist/z_height/bound)
      type(frag_ref_t) :: frag_b   ! fragment reference B (used by com_dist/bound)
      real(8) :: th_lo    = 0.0d0   ! main/lower threshold [Å] (com_dist distance; z_height height; bound lower)
      real(8) :: th_hi    = 0.0d0   ! upper threshold [Å] (bound only)
      real(8) :: v_rec    = 0.0d0   ! receding-velocity threshold [Å/(10 fs)] (com_dist optional part)
      logical :: has_v_rec = .false. ! flag: receding-velocity threshold present (grammar optional part)
      logical :: outward  = .false.  ! off-surface velocity flag (z_height optional OUTWARD: requires p_z>0)
      integer :: n_max    = 0       ! step limit [step] (for steps)
      character(len=32) :: prim_name = ''  ! registered primitive name [-] (the PRIM('name') leaf key; read
                                           ! only when node_kind == k_prim - the prim_tbl lookup key)
      type(expr_t), allocatable :: sub(:)  ! subexpressions (and/or: 2; not: 1; primitives: none)
   end type expr_t

   ! Registered system-primitive member interface (E14 track 1): the whole system criterion lives
   ! behind this signature - the state in, a boolean out; the grammar conveys nothing else to it
   ! (parameters, if any, ride the member's own closure, not the expression)
   abstract interface
      logical function prim_proc_i(sta)
         import :: state_t
         type(state_t), intent(in) :: sta   ! physical state (read-only: the system criterion's own inputs)
      end function prim_proc_i
   end interface

   ! System-primitive registry (E14 track 1: the leaf set opens - system knowledge registers through
   ! pred_reg_prim, the same drag-in shape as pes_reg/method_reg; filled in the assembly phase by the
   ! system's reg_*.f90, read-only at parse/eval time)
   integer, parameter :: max_prim = 16     ! system-primitive table capacity [row] (registered leaves;
                                           ! the four built-ins are not rows of this table)
   type :: prim_entry_t
      character(len=32) :: name = ''      ! primitive name (the PRIM('name') grammar leaf key)
      procedure(prim_proc_i), pointer, nopass :: proc => null() ! system-criterion evaluation routine
   end type prim_entry_t
   type(prim_entry_t) :: prim_tbl(max_prim)  ! system-primitive registry (module-private)
   integer :: n_prim = 0                     ! registered primitive row count [row]
contains
   !------------------------------------------------------------------
   ! pred_reg_prim(name, proc) - register a SYSTEM termination primitive (spec §7 track 1:
   !                            system knowledge lives in the blackbox; the leaf set opens, the
   !                            boolean combinator set stays closed - N2 revision). Duplicate name
   !                            or full table -> loud-stop.
   ! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §7 (two-track composition:
   !       track 1 - each blackbox registers its termination primitives, e.g. "adsorbate height
   !       above surface" written by the system author; TERMINATE= composes leaves of both tracks
   !       on equal footing); docs/final_completion_plan.md §5 N2 (the superseded closed-grammar
   !       wording); E11 drag-in shape (systems/<system>/reg_<system>_all carries the pred_reg_prim
   !       lines - a new system primitive is a new file, zero edits here)
   ! pred_parse: grammar gains primitive := ... | PRIM('name') - the leaf matches the registered
   !   table at parse time (unknown name -> grammar error naming it)
   ! pred_eval: k_prim nodes look up the table and call the member proc with sta
   subroutine pred_reg_prim(name, proc)
      character(len=*), intent(in) :: name  ! primitive name (the PRIM('name') grammar leaf key; the
                                            ! expr_t%prim_name spelling it must match)
      procedure(prim_proc_i) :: proc         ! system-criterion evaluation routine (logical function of
                                            ! the state alone - the whole criterion lives behind it)
      ! Pseudocode:
      ! 1. Duplicate name or full table -> loud-stop
      ! 2. Register the row {name, proc} into prim_tbl / n_prim (executed in the assembly phase; the
      !    four built-in primitives stay built-in leaves alongside the table - not rows of it)
   end subroutine pred_reg_prim

   !------------------------------------------------------------------
   ! pred_parse(str, expr, ios, errmsg) - parse a TERMINATE expression string -> expression tree
   ! Reference: docs/final_completion_plan.md §5 N2 (grammar, as amended by E14 - PRIM leaf, see the
   !       module head); TEST.f:106-146/:88-94 (primitive semantics)
   subroutine pred_parse(str, expr, ios, errmsg)
      character(len=*), intent(in) :: str          ! expression string (the value section of TERMINATE=)
      type(expr_t), intent(out) :: expr            ! parsed expression tree (root node)
      integer, intent(out), optional :: ios        ! parse status [-] (0 = success; nonzero = grammar error -
                                                   ! the §5 N2 error channel; without it, errors loud-stop)
      character(len=*), intent(out), optional :: errmsg ! error description (paired with ios)
      ! Pseudocode:
      ! 1. Recursive-descent expansion per the grammar: expression := primitive | AND(e1,e2) | OR(e1,e2) | NOT(e)
      !    (paren matching and comma splitting; keywords COM_DISTANCE/Z_HEIGHT/MAX_STEPS/BOUND/PRIM/
      !    AND/OR/NOT/ATOM matched in uppercase)
      ! 2. Load primitive parameters into node fields: fragment references (roster names A/B/... -> i_frag;
      !    ATOM(i) -> i_atom), th_lo/th_hi [Å], optional receding-velocity threshold -> v_rec with
      !    has_v_rec=.true.,
      !    OUTWARD -> outward=.true., step count -> n_max
      ! 3. Registered leaf (E14 track 1): PRIM('name') -> node_kind=k_prim with the quoted name stored in
      !    prim_name; the leaf matches the registered table prim_tbl AT PARSE TIME - an unknown name is a
      !    grammar error naming it (quoting / argument-passing conventions are the C7 grammar deliverable)
      ! 4. Boolean combinators attach their subexpressions recursively to sub(:) (and/or binary, not unary)
      ! 5. Grammar error (unknown primitive / wrong argument count / illegal fragment reference /
      !    unregistered PRIM name) -> set the ios error bit + an errmsg description (or loud-stop); on
      !    success ios=0
      ! 6. The default expression used when no TERMINATE keyword is present (same-parameter version of the
      !    current TEST general criterion - N2.2, preserving pre-B2 behavior) is not assembled in this
      !    routine: the assembly layer (input/driver, C3) builds the same-parameter tree from the existing
      !    parameters (old RCRIT-style thresholds)
   end subroutine pred_parse

   !------------------------------------------------------------------
   ! pred_eval(sta, expr, i_step) - evaluate the expression tree in the current physical state -> termination boolean
   ! Reference: src_VENUS/TEST.f:106-146 (primitive-evaluation prototype); TEST.f:88-94 (z_height semantics)
   logical function pred_eval(sta, expr, i_step)
      type(state_t), intent(in) :: sta          ! physical state (read-only: q/p feed the evaluation)
      type(expr_t), intent(in) :: expr          ! expression tree (from pred_parse or the assembly layer)
      integer, intent(in), optional :: i_step   ! current step index [step] (safety-net channel of the steps
                                                ! primitive, passed in by the driver evolution loop; without it
                                                ! MAX_STEPS never fires)
      ! Result: the truth value of the expression in the current state (.true. = termination/event triggered)
      ! Pseudocode:
      ! 1. Empty tree (node_kind=0) -> false (nothing assembled, nothing blocked; default-expression
      !    assembly see pred_parse note 6)
      ! 2. Boolean-combinator nodes evaluate their subexpressions recursively, then call the combinator:
      !    k_and -> pred_and(eval(sub(1)), eval(sub(2)));
      !    k_or  -> pred_or(eval(sub(1)), eval(sub(2)));
      !    k_not -> pred_not(eval(sub(1))) (left-to-right recursion order = short-circuit evaluation order)
      ! 3. Built-in primitive nodes call one of the four primitives by kind (node fields are the dummy
      !    arguments):
      !    k_com_dist -> com_dist(sta, frag_a, frag_b, th_lo, v_rec [passed when has_v_rec])
      !    k_z_height -> z_height(sta, frag_a, th_lo, outward)
      !    k_bound    -> bound(sta, frag_a, frag_b, th_lo, th_hi)
      !    k_steps    -> steps(i_step, n_max)
      ! 4. Registered leaf (E14 track 1): k_prim nodes look up prim_name in prim_tbl -> call the member
      !    proc with sta -> the node value is the function result (not found -> loud-stop naming it - the
      !    parse-time check of pred_parse step 3 makes this a belt-and-suspenders formal stop)
   end function pred_eval

   !------------------------------------------------------------------
   ! com_dist(sta, ref_a, ref_b, r_th, v_rec) - primitive: fragment-pair center-of-mass distance
   !                                           (optional receding-velocity condition)
   ! Reference: src_VENUS/TEST.f:106-146 (center-of-mass accumulation :107-130, distance/receding dot
   !       product :133-143, two-threshold test :144-145 - how the old three-valued semantics dissolve, see note 5)
   logical function com_dist(sta, ref_a, ref_b, r_th, v_rec)
      type(state_t), intent(in) :: sta          ! physical state (read-only)
      type(frag_ref_t), intent(in) :: ref_a     ! fragment reference A (roster fragment or single atom)
      type(frag_ref_t), intent(in) :: ref_b     ! fragment reference B
      real(8), intent(in) :: r_th               ! distance threshold [Å] (single-threshold form of old RBAR/RMAX)
      real(8), intent(in), optional :: v_rec    ! receding-velocity threshold [Å/(10 fs)] (not checked if absent;
                                                 ! the old RMAX termination criterion is the v_rec=0 equivalent)
      ! Result: center-of-mass distance >= r_th, and (when v_rec is present) radial receding velocity > v_rec
      ! Pseudocode:
      ! 1. Expand fragment references to atom lists: i_atom>0 -> a one-atom pseudo-fragment; otherwise the
      !    roster fragment table (roster%list_a; the B-side table is deferred with the roster - phase C)
      ! 2. Per fragment compute the mass-weighted center of mass and its velocity: q_cm = sum(w.q)/sum(w),
      !    v_cm = sum(p)/sum(w) (per-atom accumulation, old :107-130; velocity = momentum aggregate over total mass)
      ! 3. r_cm = |q_cmA - q_cmB|; v_r = (q_cmA-q_cmB).(v_cmA-v_cmB)
      !    (dot product not divided by the distance, old :133-143)
      ! 4. Truth = r_cm >= r_th; when v_rec is present also require v_r > v_rec (AND semantics)
      ! 5. Mapping of the old dual-threshold dual-event semantics: the RBAR intermediate event (old :144) is
      !    this primitive without v_rec;
      !    RMAX terminate+receding (old :145) is the v_rec=0 version - the old three-valued semantics
      !    (NTEST 0/1/2) dissolve into boolean combinations and outcome classification (outcome_tbl, C7)
   end function com_dist

   !------------------------------------------------------------------
   ! z_height(sta, ref, z_th, outward) - primitive: fragment center-of-mass height above the surface
   !                                     (optional outward velocity)
   ! Reference: src_VENUS/TEST.f:88-94 (hardwired C-atom criterion: rd=Q(3), P(3)>0, NC>50 warm-up -
   !       this primitive generalizes it: hardwired atom 1 -> fragment reference)
   logical function z_height(sta, ref, z_th, outward)
      type(state_t), intent(in) :: sta           ! physical state (read-only)
      type(frag_ref_t), intent(in) :: ref        ! fragment reference (for surface systems: the adsorbate fragment)
      real(8), intent(in) :: z_th                ! height threshold [Å] (reference: topmost surface layer at z=0)
      logical, intent(in), optional :: outward   ! true = additionally require off-surface velocity p_z>0 (old P(3)>0)
      ! Result: fragment height z >= z_th, and (when outward) fragment p_z > 0
      ! Pseudocode:
      ! 1. Expand the fragment reference (as com_dist step 1); i_atom>0 takes that atom directly:
      !    z=q(3i), p_z=p(3i) (old :89 is exactly the single-atom form rd=Q(3))
      ! 2. Fragment center-of-mass height z_cm = sum(w.q_z)/sum(w), center-of-mass momentum p_z = sum(p_z)
      !    (mass-weighted, same convention as com_dist step 2)
      ! 3. Truth = z_cm >= z_th; when outward is present also require p_z > 0 (the AND semantics of old :91-92)
      ! 4. The old :90 NC>50 warm-up skip is not carried over - protection against spurious triggering at the
      !    initial position is borne by expression composition (e.g. AND(NOT(...), ...)) or by the driver
      !    skipping the first frame; finalized at C3 assembly
   end function z_height

   !------------------------------------------------------------------
   ! bound(sta, ref_a, ref_b, r_lo, r_hi) - primitive: bounded-window escape (constrained scattering)
   ! Reference: src_VENUS/TEST.f:133-143 (distance computation reuses com_dist steps 1-3); grammar N2 spec
   logical function bound(sta, ref_a, ref_b, r_lo, r_hi)
      type(state_t), intent(in) :: sta       ! physical state (read-only)
      type(frag_ref_t), intent(in) :: ref_a  ! fragment reference A
      type(frag_ref_t), intent(in) :: ref_b  ! fragment reference B
      real(8), intent(in) :: r_lo            ! window lower threshold [Å] (too-close boundary - recombination / entering the reaction zone)
      real(8), intent(in) :: r_hi            ! window upper threshold [Å] (too-far boundary - scattering exit)
      ! Result: true = center-of-mass distance escapes [r_lo, r_hi] (triggered on leaving the bounded region;
      ! false = stay inside the window and keep evolving)
      ! Pseudocode:
      ! 1. Reuse com_dist's fragment expansion and center-of-mass distance steps (steps 1-3) to get r_cm
      ! 2. Truth = (r_cm <= r_lo) or (r_cm >= r_hi)
      ! 3. Constrained-scattering usage: TERMINATE=BOUND(A,B,r_lo,r_hi) or combined with COM_DISTANCE to
      !    implement constrained-window sampling; sub-classification of escape events goes to the
      !    outcome-classification table (C7)
   end function bound

   !------------------------------------------------------------------
   ! steps(i_step, n_max) - primitive: step-count safety net (always the fallback, replacing the old magic number)
   ! Reference: docs/final_completion_plan.md §5 N2 (MAX_STEPS - always a safety net)
   logical function steps(i_step, n_max)
      integer, intent(in) :: i_step   ! current step index [step] (passed in by the driver evolution loop - pred_eval channel)
      integer, intent(in) :: n_max    ! step limit [step]
      ! Result: true = i_step >= n_max (terminate at the limit - the ever-present safety net against
      ! runaway trajectories)
      ! Pseudocode:
      ! 1. Truth = i_step >= n_max
      ! 2. Default assembly: with no TERMINATE keyword the default expression contains this primitive
      !    (N2.2 same-parameter version)
   end function steps

   !------------------------------------------------------------------
   ! pred_and(l_a, l_b) - boolean combinator AND(e1,e2): both subexpressions true
   ! Reference: docs/final_completion_plan.md §5 N2 (grammar: AND(e1,e2))
   logical function pred_and(l_a, l_b)
      logical, intent(in) :: l_a   ! truth value of subexpression e1
      logical, intent(in) :: l_b   ! truth value of subexpression e2
      ! Result = l_a .and. l_b
      ! Pseudocode:
      ! 1. Conjunction of the two subexpression truth values (evaluation order / short-circuiting is borne
      !    by pred_eval's recursion order)
   end function pred_and

   !------------------------------------------------------------------
   ! pred_or(l_a, l_b) - boolean combinator OR(e1,e2): either subexpression true
   ! Reference: docs/final_completion_plan.md §5 N2 (grammar: OR(e1,e2))
   logical function pred_or(l_a, l_b)
      logical, intent(in) :: l_a   ! truth value of subexpression e1
      logical, intent(in) :: l_b   ! truth value of subexpression e2
      ! Result = l_a .or. l_b
      ! Pseudocode:
      ! 1. Disjunction of the two subexpression truth values
   end function pred_or

   !------------------------------------------------------------------
   ! pred_not(l_a) - boolean combinator NOT(e): negation of the subexpression
   ! Reference: docs/final_completion_plan.md §5 N2 (grammar: NOT(e))
   logical function pred_not(l_a)
      logical, intent(in) :: l_a   ! truth value of subexpression e
      ! Result = .not. l_a
      ! Pseudocode:
      ! 1. Negation of the subexpression truth value
   end function pred_not
end module predicate
