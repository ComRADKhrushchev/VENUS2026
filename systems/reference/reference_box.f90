!=====================================================================
! reference_box.f90 - isotropic harmonic well: the reference container's
!                     concrete Hamiltonian + its termination and
!                     classification judgments
! REFERENCE MATERIAL (pure mathematics, zero system knowledge): the
!   potential V(q) = K/2 * |q|^2 treats every Cartesian degree of freedom
!   alike - any system folder is interchangeable, and only the optional
!   parameter keys below change behavior. This is not a physical model of
!   anything; it exists so the whole program can be driven end to end
!   against a closed-form force law whose exact dynamics are known.
! Container slots served (the Tier 1 magic-name wrappers bind them):
!   PES           reference_pes - the composed kernel g = K*q [eV/A], the
!                 GRADIENT of the well (the interface folds the force
!                 F = -g), v = K/2*|q|^2 [eV] (both through
!                 harm_vg of reference_pot - this module carries the
!                 constant K and hands it to the member through
!                 harm_set_constant; the delivery unit contract:
!                 deliveries in eV - the eV-to-internal conversion
!                 belongs to the force interface alone, so K is carried in
!                 eV/A^2 here)
!   termination   reference_term - any |q_i| > Q_TERM (the coordinate-
!                 threshold form; bound only when Q_TERM was provided)
!   classification reference_classify - e_rel above/below E_SPLIT ->
!                 'scattering'/'capture' (bound only when E_SPLIT was
!                 provided; e_rel arrives in eV, so E_SPLIT is an eV
!                 number)
! Parameters: the declared container keys of input_qct.txt - K_HARM = <v>
!   force constant [eV/A^2] (default 0.04), Q_TERM = <v> termination
!   threshold [A] (absent -> unbound), E_SPLIT = <v> classification split
!   [eV] (absent -> unbound; follows fin%e_rel since the 2026-10-06 rebase). The Tier 1 load() wrapper declares the
!   vocabulary (keys) and pulls the buffered rows into reference_params_t at
!   the assembly seam; absent keys keep the field defaults, the found
!   flags carry the optional-slot provision semantics. The keys() wrapper
!   lazily stages the default constant at its first (pre-parse, unguarded)
!   call - the reg-era registration-time staging's anchor, so pre-seam
!   consumers see the computable default well (fix round R-C2); the seam
!   load later overwrites with the effective constant. A non-positive
!   K_HARM/Q_TERM aborts by name at reference_load (hand-written input
!   material - a typo must not silently run a different container).
! Units: see per-key list; energies in eV on the delivery channel, the
!   split in kcal/mol on the final-state channel.
!=====================================================================
module reference_box
   use state,       only: state_t
   use final_state, only: fin_t
   use input,       only: pull
   use reference_pot, only: harm_vg, harm_set_constant
   implicit none
   private
   public :: reference_params_t, reference_keys, reference_load, reference_pes
   public :: reference_term, reference_classify
   public :: reference_term_on, reference_class_on
   public :: reference_k, reference_k_default
   public :: pes, keys, load, term, classify   ! Tier 1 magic names (thin wrappers)

   ! the parameter receiver: field defaults ARE the container defaults (an
   ! absent input key leaves its field untouched at the pull; the found
   ! flags report key presence for the optional-slot semantics)
   type :: reference_params_t
      real(8) :: k_harm  = 0.04d0    ! force constant [eV/A^2] (K_HARM)
      real(8) :: q_term  = 0.0d0     ! termination threshold [A] (Q_TERM)
      logical :: q_term_found = .false.
      real(8) :: e_split = 0.0d0     ! classification split [eV] (E_SPLIT)
      logical :: e_split_found = .false.
   end type reference_params_t

   real(8), parameter :: k_default = 0.04d0 ! default force constant [eV/A^2]
   real(8) :: k_harm  = k_default ! force constant [eV/A^2]
   real(8) :: q_term  = 0.0d0     ! termination threshold [A]
   real(8) :: e_split = 0.0d0     ! classification split [eV]
   logical :: term_on  = .false.  ! Q_TERM was provided (the slot may bind)
   logical :: class_on = .false.  ! E_SPLIT was provided (the slot may bind)
   logical :: loaded   = .false.  ! reference_load already ran (double-load guard)
   logical :: staged_default = .false. ! the default well constant is already
   !                                    staged into the member (the keys()-time
   !                                    lazy guard + the pes() safety net)
contains
   !------------------------------------------------------------------
   ! reference_keys() - the container's declared-key vocabulary (the reg
   !                 file hands it to input_declare_keys before the parse)
   !------------------------------------------------------------------
   function reference_keys() result(keys)
      character(len=16) :: keys(3)
      keys = [character(len=16) :: 'K_HARM', 'Q_TERM', 'E_SPLIT']
   end function reference_keys

   !------------------------------------------------------------------
   ! reference_load(p) - take the pulled parameter set into the module
   !                 state (validations here: a non-positive K_HARM is
   !                 not a well, a non-positive Q_TERM is not a threshold).
   !                 Runs once, at the assembly seam. Whatever the
   !                 outcome, the member constant is handed to the
   !                 potential module before returning (the pulled value,
   !                 or the default when the key was absent)
   !------------------------------------------------------------------
   subroutine reference_load(p)
      type(reference_params_t), intent(in) :: p
      character(len=32) :: vs
      if (loaded) then
         call stop_box('reference_load', 'the container parameters were already'// &
                       ' loaded (a second load - an assembly-choreography error, formal)')
      end if
      if (p%k_harm <= 0.0d0) then
         write (vs, '(g0)') p%k_harm
         call stop_box('reference_load', 'K_HARM must be positive (got '//trim(vs)// &
                       ' eV/A^2 - a non-positive constant is not a well)')
      end if
      if (p%q_term_found .and. p%q_term <= 0.0d0) then
         write (vs, '(g0)') p%q_term
         call stop_box('reference_load', 'Q_TERM must be positive (got '//trim(vs)//' A)')
      end if
      k_harm   = p%k_harm
      q_term   = p%q_term
      e_split  = p%e_split
      term_on  = p%q_term_found
      class_on = p%e_split_found
      call harm_set_constant(k_harm)          ! the member init hand-off (always)
      staged_default = .true.                 ! the effective constant owns the member
      !                                         now - the lazy guard must never overwrite it
      loaded = .true.
   end subroutine reference_load

   !------------------------------------------------------------------
   ! reference_pes(q, mass, g, v) - the container's aggregate force entry: the
   !                 isotropic well, computed by the member harm_vg (which
   !                 carries the constant this module buffered through
   !                 harm_set_constant). The delivery contract fully assigns
   !                 g and v; mass is unused (a pure potential)
   !------------------------------------------------------------------
   subroutine reference_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A] (all dof, flattened)
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A] (the interface folds F = -g)
      real(8), intent(out) :: v        ! V [eV]
      call harm_vg(size(q)/3, q, v, g)
   end subroutine reference_pes

   !------------------------------------------------------------------
   ! reference_term(sta) - the termination judgment: any coordinate leaves
   !                 the threshold ball (bound only when Q_TERM buffered)
   !------------------------------------------------------------------
   logical function reference_term(sta)
      type(state_t), intent(in) :: sta
      reference_term = .false.
      if (.not. term_on) return
      reference_term = any(abs(sta%q) > q_term)
   end function reference_term

   !------------------------------------------------------------------
   ! reference_classify(fin) - the classification judgment: relative
   !                 translational energy above the split -> 'scattering',
   !                 otherwise 'capture' (read-only final-state quantities)
   !------------------------------------------------------------------
   function reference_classify(fin) result(name)
      type(fin_t), intent(in) :: fin  ! final-state quantities (read-only)
      character(len=32) :: name       ! the outcome name [-]
      if (fin%e_rel > e_split) then
         name = 'scattering'
      else
         name = 'capture'
      end if
   end function reference_classify

   !------------------------------------------------------------------
   ! reference_term_on() / reference_class_on() - the optional-slot
   !                 provision queries (the drag-in file consults them to
   !                 decide which slots to bind)
   !------------------------------------------------------------------
   pure logical function reference_term_on()
      reference_term_on = term_on
   end function reference_term_on

   pure logical function reference_class_on()
      reference_class_on = class_on
   end function reference_class_on

   !------------------------------------------------------------------
   ! reference_k() / reference_k_default() - the constant getters (the
   !                 assembly verification face: what the input buffered /
   !                 what the default is)
   !------------------------------------------------------------------
   pure function reference_k() result(k)
      real(8) :: k                          ! the buffered force constant [eV/A^2]
      k = k_harm
   end function reference_k

   pure function reference_k_default() result(k)
      real(8) :: k                          ! the default force constant [eV/A^2]
      k = k_default
   end function reference_k_default

   !------------------------------------------------------------------
   ! Tier 1 magic-name wrappers (the auto-weave contract: a reg-less
   !                 systems/ folder gets these five names woven into
   !                 registry_gen; each delegates to the public API above,
   !                 which the check programs and probe hooks keep using)
   !------------------------------------------------------------------
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      if (.not. staged_default) then    ! pre-seam consumers: the default well is
         call harm_set_constant(reference_k_default())   ! computable (the ruled
         staged_default = .true.        ! safety net; inert in woven flows - the
      end if                            ! keys() anchor already staged)
      call reference_pes(q, mass, g, v)
   end subroutine pes

   function keys() result(kk)
      character(len=16) :: kk(3)
      if (.not. staged_default) then    ! the earliest woven call (pre-parse,
         call harm_set_constant(reference_k_default())   ! unguarded): stage the
         staged_default = .true.        ! default well so pre-seam consumers see
      end if                            ! a computable constant - the reg-era
      kk = reference_keys()             ! registration-time staging's Tier 1
   end function keys                    ! anchor (fix round R-C2; a pes-anchored
                                        ! guard would flip the member constant
                                        ! inside a consumer's evaluation loop -
                                        ! see rc2-verbatim-facade.log; the seam
                                        ! load later overwrites with the
                                        ! effective constant)

   subroutine load()
      type(reference_params_t) :: p
      call pull('K_HARM',  p%k_harm)
      call pull('Q_TERM',  p%q_term,  p%q_term_found)
      call pull('E_SPLIT', p%e_split, p%e_split_found)
      call reference_load(p)
   end subroutine load

   logical function term(sta)
      type(state_t), intent(in) :: sta
      term = reference_term(sta)   ! direct call: the body self-gates on
      !                               term_on (an unbound run stays .false.)
   end function term

   function classify(fin) result(name)
      type(fin_t), intent(in) :: fin
      character(len=32) :: name
      if (reference_class_on()) then
         name = reference_classify(fin)
      else
         name = 'capture'         ! the no-E_SPLIT catch-all: the classify
      end if                      ! body's else-branch channel (no split
   end function classify          ! criterion can name 'scattering')

   !------------------------------------------------------------------
   ! stop_box(who, msg) - the container abort channel (name the caller,
   !                 print the named message, STOP 1 - the sibling
   !                 modules' semantics)
   !------------------------------------------------------------------
   subroutine stop_box(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': container parameter error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container is not usable)'
      stop 1
   end subroutine stop_box
end module reference_box
