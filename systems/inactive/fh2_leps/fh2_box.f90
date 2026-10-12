!=====================================================================
! fh2_box.f90 - the F+H2 container's concrete wiring: the aggregate
!               force entry over fh2_pot, the parameter channel, the
!               fragment-separation termination and the HF-channel
!               classification
! Design:
!   Atom order (the system-folder contract): atom 1 = F, atoms 2,3 = H.
!   Parameters: the declared container keys of input_qct.txt (the Tier 1
!   magic-name wrappers declare the vocabulary and pull the buffered rows
!   into fh2_params_t at the assembly seam; absent keys keep the legislated
!   defaults of fh2_pot). Keys:
!     D_FH/R0_FH/BETA_FH, D_HH/R0_HH/BETA_HH  [eV / A / 1/A]
!     DL0_HH/DL0_FH/AMOD                      [-] (the T-S Sato set;
!     all nine surface keys must appear together - a half-replaced
!     surface is a transcription error)
!     R_TERM = <v>  termination: any pair distance beyond [A] (absent -> unbound)
!     R_HF   = <v>  classification: HF formed when r(F,H) below [A] (absent -> unbound)
!   Every provided value must be positive (validated at fh2_load).
!   eV on the delivery channel (the eV-to-internal conversion belongs to
!   the force interface alone); R_TERM/R_HF are run material, not system
!   definition. Termination = any pair separated (scattering done);
!   classification reads the COM-frame archive fin%q_fin (bond lengths
!   are translation invariant): an F-H pair below R_HF names the channel
!   'reaction' (HF formed), else 'nonreactive'.
!=====================================================================
module fh2_box
   use state,       only: state_t
   use final_state, only: fin_t
   use input,       only: pull
   use fh2_pot,     only: fh2_vg, fh2_set_params, fh2_set_defaults
   implicit none
   private
   public :: fh2_params_t, fh2_keys, fh2_load, fh2_pes, fh2_term, fh2_classify
   public :: fh2_term_on, fh2_class_on
   public :: pes, keys, load, term, classify   ! Tier 1 magic names (thin wrappers)

   ! the parameter receiver: field defaults ARE the container defaults (an
   ! absent input key leaves its field untouched at the pull; the found
   ! flags carry the optional-slot and all-or-none surface-set semantics)
   type :: fh2_params_t
      real(8) :: r_term = 8.0d0          ! R_TERM [A]
      real(8) :: r_hf   = 1.6d0          ! R_HF [A]
      logical :: r_term_found = .false.
      logical :: r_hf_found   = .false.
      real(8) :: surf(9) = 0.0d0         ! D_FH/R0_FH/BETA_FH/D_HH/R0_HH/BETA_HH/DL0_HH/DL0_FH/AMOD
      logical :: surf_found(9) = .false.
   end type fh2_params_t

   real(8) :: r_term = 8.0d0        ! termination pair-distance threshold [A]
   real(8) :: r_hf   = 1.6d0        ! classification F-H bond threshold [A]
   logical :: term_on  = .false.    ! R_TERM was provided (the slot may bind)
   logical :: class_on = .false.    ! R_HF was provided (the slot may bind)
   logical :: loaded   = .false.    ! fh2_load already ran (double-load guard)
contains
   !------------------------------------------------------------------
   ! fh2_term_on() / fh2_class_on() - the optional-slot availability
   !                 flags (the drag-in consults them at the seam)
   !------------------------------------------------------------------
   logical function fh2_term_on()
      fh2_term_on = term_on
   end function fh2_term_on

   logical function fh2_class_on()
      fh2_class_on = class_on
   end function fh2_class_on

   !------------------------------------------------------------------
   ! fh2_keys() - the container's declared-key vocabulary (the reg file
   !                 hands it to input_declare_keys before the parse)
   !------------------------------------------------------------------
   function fh2_keys() result(keys)
      character(len=16) :: keys(11)
      keys = [character(len=16) :: 'R_TERM', 'R_HF', 'D_FH', 'R0_FH', 'BETA_FH', &
              'D_HH', 'R0_HH', 'BETA_HH', 'DL0_HH', 'DL0_FH', 'AMOD']
   end function fh2_keys

   !------------------------------------------------------------------
   ! fh2_load(p) - take the pulled parameter set into the module state:
   !                 every provided value must be positive, the nine
   !                 surface keys must appear together (a half-replaced
   !                 surface is a transcription error), then the surface
   !                 hand-off runs (the override set or the legislated
   !                 defaults). Runs once, at the assembly seam
   !------------------------------------------------------------------
   subroutine fh2_load(p)
      type(fh2_params_t), intent(in) :: p
      character(len=32) :: vs
      if (loaded) then
         call stop_box('fh2_load', 'the container parameters were already'// &
                       ' loaded (a second load - an assembly-choreography error, formal)')
      end if
      if (p%r_term_found .and. p%r_term <= 0.0d0) then
         write (vs, '(g0)') p%r_term
         call stop_box('fh2_load', 'R_TERM must be positive (got '//trim(vs)//')')
      end if
      if (p%r_hf_found .and. p%r_hf <= 0.0d0) then
         write (vs, '(g0)') p%r_hf
         call stop_box('fh2_load', 'R_HF must be positive (got '//trim(vs)//')')
      end if
      if (any(p%surf_found .and. p%surf <= 0.0d0)) then
         write (vs, '(g0)') minval(p%surf, mask=p%surf_found)
         call stop_box('fh2_load', 'a surface key carried a non-positive value (got '// &
                       trim(vs)//' - all of D_FH/R0_FH/BETA_FH/D_HH/R0_HH/BETA_HH/'// &
                       'DL0_HH/DL0_FH/AMOD must be positive)')
      end if
      if (any(p%surf_found) .and. .not. all(p%surf_found)) then
         call stop_box('fh2_load', 'a partial surface-parameter override '// &
                       '(all nine of D_FH/R0_FH/BETA_FH/D_HH/R0_HH/BETA_HH/'// &
                       'DL0_HH/DL0_FH/AMOD must appear together)')
      end if
      r_term = p%r_term
      r_hf = p%r_hf
      term_on = p%r_term_found
      class_on = p%r_hf_found
      if (all(p%surf_found)) then
         call fh2_set_params(p%surf(1), p%surf(2), p%surf(3), p%surf(4), &
                             p%surf(5), p%surf(6), p%surf(7), p%surf(8), p%surf(9))
      else
         call fh2_set_defaults()
      end if
      loaded = .true.
   end subroutine fh2_load

   !------------------------------------------------------------------
   ! fh2_pes(q, mass, g, v) - the container's aggregate force entry: the
   !                 modified-LEPS surface over fh2_pot (F + H + H, three
   !                 atoms). The delivery contract fully assigns g and v;
   !                 mass is unused (a pure potential)
   !------------------------------------------------------------------
   subroutine fh2_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A] (all dof, flattened)
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A] (the interface folds F = -g)
      real(8), intent(out) :: v        ! V [eV]
      if (size(q) /= 9 .or. size(mass) /= 3) then
         call stop_box('fh2_pes', 'the F+H2 container needs exactly 3 atoms (9 dofs)'// &
                       ' - the system folder must define F, H, H in this order')
      end if
      call fh2_vg(3, q, v, g)
   end subroutine fh2_pes

   !------------------------------------------------------------------
   ! fh2_term(sta) - the termination judgment: any atom pair separated
   !                 beyond R_TERM (the scattering is over). Reached only
   !                 when the drag-in bound the slot
   !------------------------------------------------------------------
   logical function fh2_term(sta)
      type(state_t), intent(in) :: sta  ! physical state (read-only snapshot)
      real(8) :: r12, r13, r23
      fh2_term = .false.
      if (.not. term_on) return
      r12 = sqrt(dot_product(sta%q(4:6) - sta%q(1:3), sta%q(4:6) - sta%q(1:3)))
      r13 = sqrt(dot_product(sta%q(7:9) - sta%q(1:3), sta%q(7:9) - sta%q(1:3)))
      r23 = sqrt(dot_product(sta%q(7:9) - sta%q(4:6), sta%q(7:9) - sta%q(4:6)))
      fh2_term = (r12 > r_term) .or. (r13 > r_term) .or. (r23 > r_term)
   end function fh2_term

   !------------------------------------------------------------------
   ! fh2_classify(fin) - the channel judgment: an F-H pair below R_HF
   !                 names 'reaction' (HF formed), else 'nonreactive'
   !                 (COM-frame archive; bond lengths are invariant)
   !------------------------------------------------------------------
   function fh2_classify(fin) result(name)
      type(fin_t), intent(in) :: fin  ! final-state quantities (read-only)
      character(len=32) :: name       ! the outcome name [-]
      real(8) :: r12, r13
      if (.not. allocated(fin%q_fin)) then
         call stop_box('fh2_classify', 'the final coordinate archive is not allocated'// &
                       ' (the classification reads fin%q_fin)')
      end if
      r12 = sqrt(dot_product(fin%q_fin(4:6) - fin%q_fin(1:3), fin%q_fin(4:6) - fin%q_fin(1:3)))
      r13 = sqrt(dot_product(fin%q_fin(7:9) - fin%q_fin(1:3), fin%q_fin(7:9) - fin%q_fin(1:3)))
      if (r12 < r_hf .or. r13 < r_hf) then
         name = 'reaction'
      else
         name = 'nonreactive'
      end if
   end function fh2_classify

   !------------------------------------------------------------------
   ! Tier 1 magic-name wrappers (the auto-weave contract: a reg-less
   !                 systems/ folder gets these five names woven into
   !                 registry_gen; each delegates to the public API above,
   !                 which the check programs and probe hooks keep using)
   !------------------------------------------------------------------
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      call fh2_pes(q, mass, g, v)
   end subroutine pes

   function keys() result(kk)
      character(len=16) :: kk(11)
      kk = fh2_keys()
   end function keys

   subroutine load()
      type(fh2_params_t) :: p
      call pull('R_TERM',  p%r_term,  p%r_term_found)
      call pull('R_HF',    p%r_hf,   p%r_hf_found)
      call pull('D_FH',    p%surf(1), p%surf_found(1))
      call pull('R0_FH',   p%surf(2), p%surf_found(2))
      call pull('BETA_FH', p%surf(3), p%surf_found(3))
      call pull('D_HH',    p%surf(4), p%surf_found(4))
      call pull('R0_HH',   p%surf(5), p%surf_found(5))
      call pull('BETA_HH', p%surf(6), p%surf_found(6))
      call pull('DL0_HH',  p%surf(7), p%surf_found(7))
      call pull('DL0_FH',  p%surf(8), p%surf_found(8))
      call pull('AMOD',    p%surf(9), p%surf_found(9))
      call fh2_load(p)
   end subroutine load

   logical function term(sta)
      type(state_t), intent(in) :: sta
      term = .false.
      if (fh2_term_on()) term = fh2_term(sta)
   end function term

   function classify(fin) result(name)
      type(fin_t), intent(in) :: fin
      character(len=32) :: name
      if (fh2_class_on()) then
         name = fh2_classify(fin)
      else
         name = 'nonreactive'   ! the no-R_HF catch-all: the classify body's
                                ! else-branch channel (no HF criterion can
                                ! name 'reaction')
      end if
   end function classify

   !------------------------------------------------------------------
   ! stop_box(who, msg) - the container abort channel (name the caller,
   !                 print the named message, STOP 1)
   !------------------------------------------------------------------
   subroutine stop_box(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': fh2 container error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container is not usable)'
      stop 1
   end subroutine stop_box
end module fh2_box
