!=====================================================================
! analysis/outcome.f90 - L4 analysis: outcome classification (outcome_tbl data-table driven -
!                          row-by-row predicate evaluation; the old duplicated outcome-count
!                          blocks turned into data)
! (see previous line)
! Reference: src_VENUS/VENUS.f90:252-284 (outcome summary printout: NONREACT/NREACT(IPATH)/
!       NLONG/NEFF - prototype of the statistics-phase count tables), :552-571 and :645-657
!       (duplicated
!       reactive/non-reactive count blocks - nearly identical, turned into table rows),
!       :559-565 (energy-conservation gate
!       DELM - drift-threshold predicate prototype), :669-685 (the long-time NLONG block -
!       lifetime-threshold predicate),
!       :669-673 (the -1 invalid placeholder of long-time final-state rows - data-column
!       semantics);
!       src_VENUS/TEST.f:106-146 (NPATH outcome-assignment criterion prototype - fragment COM
!       distance;
!       how the NTEST three-valued semantics dissolve, see predicate.f90 note 5);
!       config%event%outcome_tbl (K1 §1.4 event provisions - "data not code"; wired in C7)
! Excised: the old NPATH integer outcomes (1 = non-reactive, >1 = reactive path number) and the
! NTEST/NFINAL exit
!       semantics - all turned into predicate rows of the outcome table; outcome names are data,
!       never enumerated in code
!=====================================================================
module outcome
   use config,      only: event
   use final_state, only: fin_t
   implicit none
   private
   public :: outc_expr_t, outc_row_t, outc_tbl_t, outc_classify
   public :: o_e_rel, o_b_fin, o_n_a, o_n_b, o_e_drift, o_t_life, o_and, o_or, o_not

   ! Predicate-leaf observable keys (values of outc_expr_t%key; a closed set - new keys are a
   ! grammar-versioning matter)
   integer, parameter :: o_e_rel   = 1  ! final relative translational energy e_rel [kcal/mol] (scattering outcomes)
   integer, parameter :: o_b_fin   = 2  ! final impact parameter b_fin [Å]
   integer, parameter :: o_n_a     = 3  ! product-A atom count [count] (fragment-rearrangement outcomes - the old NPATH assignment)
   integer, parameter :: o_n_b     = 4  ! product-B atom count [count]
   integer, parameter :: o_e_drift = 5  ! energy-conservation drift H_max-H_min [kcal/mol]
                                        ! (old DELM gate
   ! VENUS.f90:559-565 - values come from trajectory
   ! aggregation of the recorder energy columns,
   ! channel finalized in C8)
   integer, parameter :: o_t_life  = 6  ! trajectory lifetime [fs] (old NLONG gate :669-685)
   integer, parameter :: o_and     = 7  ! combinator AND(e1,e2)
   integer, parameter :: o_or      = 8  ! combinator OR(e1,e2)
   integer, parameter :: o_not     = 9  ! combinator NOT(e)

   ! Outcome predicate expression node (recursive tree - same structure as predicate%expr_t,
   ! leaf keys drawn from the final-state quantity domain)
   type :: outc_expr_t
      integer :: key = 0           ! observable/combinator key (o_*)
      integer :: cmp = 0           ! comparison direction [-] (1 = >=, 2 = <=; unused for combinators)
      real(8) :: th  = 0.0d0       ! threshold [units of the key]
      type(outc_expr_t), allocatable :: sub(:)  ! subexpressions (and/or: 2, not: 1, leaves: none)
   end type outc_expr_t

   ! Outcome table row: outcome name + row condition (table order is priority - first true row
   ! wins)
   type :: outc_row_t
      character(len=32) :: name = ''  ! outcome name (e.g. 'nonreact'/'react'/'long' - data)
      type(outc_expr_t) :: cond       ! row-condition predicate
   end type outc_row_t

   type :: outc_tbl_t
      type(outc_row_t), allocatable :: row(:)  ! row table
      integer :: n_row = 0                       ! row count [row]
      character(len=32) :: name_def = 'unclassified'  ! all-false default name [-]
   end type outc_tbl_t
contains
   !------------------------------------------------------------------
   ! outc_classify(fin, tbl) - outcome classification: row-by-row predicate evaluation, the
   ! first true row's name wins
   ! Reference: VENUS.f90:552-571/:645-657/:669-685 (the three count blocks as data);
   !       TEST.f:106-146 (react/nonreact criterion prototype)
   function outc_classify(fin, tbl) result(name)
      type(fin_t), intent(in) :: fin        ! final-state quantities (fin_eval output - the
      ! leaf value domain)
      type(outc_tbl_t), intent(in) :: tbl   ! parsed outcome table (source =
                                             ! config%event%outcome_tbl - parsed once in
                                             ! the assembly phase; the string -> table parser
                                             ! follows the predicate grammar C4/C7, same
                                             ! family as predicate's pred_parse,
      ! no separate routine at this stage)
      character(len=32) :: name             ! result: outcome name [-]
      ! Pseudocode:
      ! 1. Per row r = 1..tbl%n_row (table order is priority): evaluate cond recursively -
      !    leaves fetch by key: o_e_rel -> fin%e_rel, o_b_fin -> fin%b_fin,
      !    o_n_a/o_n_b -> product atom counts (roster fragment counts - the atom-assignment
      !    analysis of fragment rearrangement
      !    finalized in C7, the old NPATH detected membership changes via TEST),
      !    o_e_drift -> conservation drift
      !    (recorder energy-column aggregation - C8 channel), o_t_life -> fin%t_life;
      !    compare against th per cmp (>=/<=)
      ! 2. Combinator keys recurse on subexpressions: o_and -> both true, o_or -> either true,
      !    o_not -> negation
      !    (evaluation order = short-circuit order, as in predicate%pred_eval note 2)
      ! 3. First true row -> name = tbl%row(r)%name, return (the old react/nonreact dichotomy
      !    is simply the minimal two-row table)
      !    (see previous line)
      ! 4. All false -> name = tbl%name_def (the long-time 'long' row can be an explicit row
      !    or the default row -
      !    counterpart of the old NLONG block :669-685)
   end function outc_classify
end module outcome
