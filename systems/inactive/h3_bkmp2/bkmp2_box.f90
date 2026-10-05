!=====================================================================
! bkmp2_box.f90 - the H+H2 BKMP2 container's concrete wiring: the
!                 aggregate force entry over bkmp2_pot, the parameter
!                 channel, the pair-separation termination and the
!                 exchange-reaction classification
! Design:
!   Atom order (the system-folder contract): atoms 1,2 = the H2 molecule,
!   atom 3 = the incident H. All three atoms are H - the BKMP2 surface is
!   permutation symmetric, the labeling matters only to the channel
!   classification below. Parameters: the declared container keys of
!   input_qct.txt (the reg file declares the vocabulary and pulls the
!   buffered rows into bkmp2_params_t at the assembly seam; absent keys
!   leave both optional slots unbound). Keys:
!     R_TERM = <v>  termination: any pair distance beyond [A] (absent -> unbound)
!     R_H2   = <v>  classification: H-H bond threshold [A] (absent -> unbound)
!   Both provided values must be positive (validated at bkmp2_load).
!   eV on the delivery channel (the eV-to-internal conversion belongs to
!   the force interface alone); R_TERM/R_H2 are run material, not system
!   definition. Termination = any pair separated (the scattering is over);
!   classification reads the COM-frame archive fin%q_fin (bond lengths are
!   translation invariant): the incident atom bonded to atom 1 or 2 names
!   the channel 'reaction' (exchange), else the molecule pair bonded names
!   it 'nonreactive', else 'dissociation' (no H-H pair bound).
!=====================================================================
module bkmp2_box
   use state,       only: state_t
   use final_state, only: fin_t
   use bkmp2_pot,   only: bkmp2_vg
   implicit none
   private
   public :: bkmp2_params_t, bkmp2_keys, bkmp2_load, bkmp2_pes, bkmp2_term, bkmp2_classify
   public :: bkmp2_term_on, bkmp2_class_on

   ! the parameter receiver: field defaults ARE the container defaults (an
   ! absent input key leaves its field untouched at the pull; the found
   ! flags carry the optional-slot provision semantics)
   type :: bkmp2_params_t
      real(8) :: r_term = 8.0d0          ! R_TERM [A]
      real(8) :: r_h2   = 1.2d0          ! R_H2 [A]
      logical :: r_term_found = .false.
      logical :: r_h2_found   = .false.
   end type bkmp2_params_t

   real(8) :: r_term = 8.0d0        ! termination pair-distance threshold [A]
   real(8) :: r_h2   = 1.2d0        ! classification H-H bond threshold [A]
   logical :: term_on  = .false.    ! R_TERM was provided (the slot may bind)
   logical :: class_on = .false.    ! R_H2 was provided (the slot may bind)
   logical :: loaded   = .false.    ! bkmp2_load already ran (double-load guard)
contains
   !------------------------------------------------------------------
   ! bkmp2_term_on() / bkmp2_class_on() - the optional-slot availability
   !                 flags (the drag-in consults them at the seam)
   !------------------------------------------------------------------
   logical function bkmp2_term_on()
      bkmp2_term_on = term_on
   end function bkmp2_term_on

   logical function bkmp2_class_on()
      bkmp2_class_on = class_on
   end function bkmp2_class_on

   !------------------------------------------------------------------
   ! bkmp2_keys() - the container's declared-key vocabulary (the reg file
   !                 hands it to input_declare_keys before the parse)
   !------------------------------------------------------------------
   function bkmp2_keys() result(keys)
      character(len=16) :: keys(2)
      keys = [character(len=16) :: 'R_TERM', 'R_H2']
   end function bkmp2_keys

   !------------------------------------------------------------------
   ! bkmp2_load(p) - take the pulled parameter set into the module state
   !                 (validations here: a non-positive threshold is not a
   !                 threshold). Runs once, at the assembly seam. The
   !                 BKMP2 surface constants are not run material (the
   !                 archived fit is the definition), so the keys carry
   !                 the two judgment thresholds only
   !------------------------------------------------------------------
   subroutine bkmp2_load(p)
      type(bkmp2_params_t), intent(in) :: p
      character(len=32) :: vs
      if (loaded) then
         call stop_box('bkmp2_load', 'the container parameters were already'// &
                       ' loaded (a second load - an assembly-choreography error, formal)')
      end if
      if (p%r_term_found .and. p%r_term <= 0.0d0) then
         write (vs, '(g0)') p%r_term
         call stop_box('bkmp2_load', 'R_TERM must be positive (got '//trim(vs)//')')
      end if
      if (p%r_h2_found .and. p%r_h2 <= 0.0d0) then
         write (vs, '(g0)') p%r_h2
         call stop_box('bkmp2_load', 'R_H2 must be positive (got '//trim(vs)//')')
      end if
      r_term = p%r_term
      r_h2 = p%r_h2
      term_on = p%r_term_found
      class_on = p%r_h2_found
      loaded = .true.
   end subroutine bkmp2_load

   !------------------------------------------------------------------
   ! bkmp2_pes(q, mass, g, v) - the container's aggregate force entry: the
   !                 BKMP2 surface over bkmp2_pot (H + H + H, three atoms).
   !                 The delivery contract fully assigns g and v; mass is
   !                 unused (a pure potential)
   !------------------------------------------------------------------
   subroutine bkmp2_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A] (all dof, flattened)
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A] (the interface folds F = -g)
      real(8), intent(out) :: v        ! V [eV]
      if (size(q) /= 9 .or. size(mass) /= 3) then
         call stop_box('bkmp2_pes', 'the H+H2 container needs exactly 3 atoms (9 dofs)'// &
                       ' - the system folder must define H, H, H in this order')
      end if
      call bkmp2_vg(3, q, v, g)
   end subroutine bkmp2_pes

   !------------------------------------------------------------------
   ! bkmp2_term(sta) - the termination judgment: any atom pair separated
   !                 beyond R_TERM (the scattering is over). Reached only
   !                 when the drag-in bound the slot
   !------------------------------------------------------------------
   logical function bkmp2_term(sta)
      type(state_t), intent(in) :: sta  ! physical state (read-only snapshot)
      real(8) :: r12, r13, r23
      bkmp2_term = .false.
      if (.not. term_on) return
      r12 = sqrt(dot_product(sta%q(4:6) - sta%q(1:3), sta%q(4:6) - sta%q(1:3)))
      r13 = sqrt(dot_product(sta%q(7:9) - sta%q(1:3), sta%q(7:9) - sta%q(1:3)))
      r23 = sqrt(dot_product(sta%q(7:9) - sta%q(4:6), sta%q(7:9) - sta%q(4:6)))
      bkmp2_term = (r12 > r_term) .or. (r13 > r_term) .or. (r23 > r_term)
   end function bkmp2_term

   !------------------------------------------------------------------
   ! bkmp2_classify(fin) - the channel judgment over the atom labeling
   !                 (COM-frame archive; bond lengths are invariant):
   !                 incident atom 3 bonded to atom 1 or 2 -> 'reaction'
   !                 (exchange); else molecule pair 1-2 bonded ->
   !                 'nonreactive'; else 'dissociation'
   !------------------------------------------------------------------
   function bkmp2_classify(fin) result(name)
      type(fin_t), intent(in) :: fin  ! final-state quantities (read-only)
      character(len=32) :: name       ! the outcome name [-]
      real(8) :: r12, r13, r23
      if (.not. allocated(fin%q_raw)) then
         call stop_box('bkmp2_classify', 'the raw final coordinate archive is not'// &
                       ' allocated (the classification reads fin%q_raw)')
      end if
      r12 = sqrt(dot_product(fin%q_raw(4:6) - fin%q_raw(1:3), fin%q_raw(4:6) - fin%q_raw(1:3)))
      r13 = sqrt(dot_product(fin%q_raw(7:9) - fin%q_raw(1:3), fin%q_raw(7:9) - fin%q_raw(1:3)))
      r23 = sqrt(dot_product(fin%q_raw(7:9) - fin%q_raw(4:6), fin%q_raw(7:9) - fin%q_raw(4:6)))
      ! atom labels: 1 = projectile H (from H.xyz), 2,3 = target H2 (from H2.xyz);
      ! r12/r13 = projectile-to-target distances, r23 = the target H2's own bond.
      ! Exchange (reaction) = the projectile ends up bonded to a target atom;
      ! nonreactive = the original target pair (2,3) stays intact, projectile free
      if (r12 < r_h2 .or. r13 < r_h2) then
         name = 'reaction'
      else if (r23 < r_h2) then
         name = 'nonreactive'
      else
         name = 'dissociation'
      end if
   end function bkmp2_classify

   !------------------------------------------------------------------
   ! stop_box(who, msg) - the fatal channel (name the caller, print the
   !                 named message, STOP 1)
   !------------------------------------------------------------------
   subroutine stop_box(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': container error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container is not usable)'
      stop 1
   end subroutine stop_box
end module bkmp2_box
