!=====================================================================
! h2ag111_box.f90 - the H2/Ag(111) container's concrete wiring: the
!                   one-time network initialization (the run-directory
!                   data files), the aggregate force entry, the
!                   desorption termination and the molecular-vs-
!                   dissociated classification
! Design:
!   Two dynamic atoms only (the Ag(111) slab is rigid inside the PES):
!   q(1:3) = H1, q(4:6) = H2, z the height over the surface. Units are
!   native A/eV end to end (no conversion in this container). The
!   network parameters are RUN material: h2ag111_load() requires
!   weights-h2ag.txt / biases-h2ag.txt in the run directory (copy them
!   from the container folder) and calls pes_init exactly once.
!   Parameters: the declared container keys of input_qct.txt (the reg
!   file declares the vocabulary and pulls the buffered rows into
!   h2ag111_params_t at the assembly seam; absent keys leave both
!   optional slots unbound). Keys:
!     Z_TERM = <v>  termination: BOTH atom heights above [A] (the
!                  molecule left the surface region; absent -> unbound;
!                  keep it inside the network's fitted domain - the PES
!                  carries no asymptotic tail)
!     R_H2   = <v>  classification: H-H bond threshold [A] (absent -> unbound)
!   Both provided values must be positive (validated at h2ag111_load).
!   Termination = desorption (both z above Z_TERM); classification reads
!   the COM-frame archive fin%q_fin: r(H1,H2) below R_H2 names the
!   channel 'scattering' (molecule intact), else 'dissociation'
!   (dissociated chemisorption). eV on the delivery channel (the
!   eV-to-internal conversion belongs to the force interface alone).
!=====================================================================
module h2ag111_box
   use state,       only: state_t
   use final_state, only: fin_t
   implicit none
   private
   public :: h2ag111_params_t, h2ag111_keys, h2ag111_load, h2ag111_pes, &
             h2ag111_term, h2ag111_classify
   public :: h2ag111_term_on, h2ag111_class_on

   ! the parameter receiver: field defaults ARE the container defaults (an
   ! absent input key leaves its field untouched at the pull; the found
   ! flags carry the optional-slot provision semantics)
   type :: h2ag111_params_t
      real(8) :: z_term = 5.0d0          ! Z_TERM [A]
      real(8) :: r_h2   = 1.2d0          ! R_H2 [A]
      logical :: z_term_found = .false.
      logical :: r_h2_found   = .false.
   end type h2ag111_params_t

   ! the archived PES entries (this folder's h2ag111_pot.f90, external
   ! subroutines - assumed-size interfaces; units A / eV / eV-per-A)
   interface
      subroutine pot0(natoms, q, v)
         integer :: natoms
         real(8) :: q(*), v
      end subroutine pot0
      subroutine dpeshon(natoms, q, pdot)
         integer :: natoms
         real(8) :: q(*), pdot(*)
      end subroutine dpeshon
      subroutine pes_init()
      end subroutine pes_init
   end interface

   real(8) :: z_term = 5.0d0        ! termination height threshold [A]
   real(8) :: r_h2   = 1.2d0        ! classification H-H bond threshold [A]
   logical :: term_on  = .false.    ! Z_TERM was provided (the slot may bind)
   logical :: class_on = .false.    ! R_H2 was provided (the slot may bind)
   logical :: loaded   = .false.    ! h2ag111_load already ran (double-load guard)
contains
   !------------------------------------------------------------------
   ! h2ag111_term_on() / h2ag111_class_on() - the optional-slot
   !                 availability flags (the drag-in consults them)
   !------------------------------------------------------------------
   logical function h2ag111_term_on()
      h2ag111_term_on = term_on
   end function h2ag111_term_on

   logical function h2ag111_class_on()
      h2ag111_class_on = class_on
   end function h2ag111_class_on

   !------------------------------------------------------------------
   ! h2ag111_keys() - the container's declared-key vocabulary (the reg
   !                 file hands it to input_declare_keys before the parse)
   !------------------------------------------------------------------
   function h2ag111_keys() result(keys)
      character(len=16) :: keys(2)
      keys = [character(len=16) :: 'Z_TERM', 'R_H2']
   end function h2ag111_keys

   !------------------------------------------------------------------
   ! h2ag111_load(p) - the one-time container initialization at the
   !                 assembly seam: (1) require the two network data
   !                 files in the run directory, (2) call pes_init once,
   !                 (3) take the pulled judgment thresholds into the
   !                 module state. Every defect a named fatal
   !------------------------------------------------------------------
   subroutine h2ag111_load(p)
      type(h2ag111_params_t), intent(in) :: p
      logical :: exist
      character(len=32) :: vs
      if (loaded) then
         call stop_box('h2ag111_load', 'the container was already loaded'// &
                       ' (a second load - an assembly-choreography error, formal)')
      end if
      ! 1. the network parameters are run material (pes_init reads them
      !    cwd-relative by its archived contract)
      inquire (file='weights-h2ag.txt', exist=exist)
      if (.not. exist) then
         call stop_box('h2ag111_load', 'weights-h2ag.txt not found in the run'// &
                       ' directory (copy it from the h2ag111PES container folder)')
      end if
      inquire (file='biases-h2ag.txt', exist=exist)
      if (.not. exist) then
         call stop_box('h2ag111_load', 'biases-h2ag.txt not found in the run'// &
                       ' directory (copy it from the h2ag111PES container folder)')
      end if
      ! 2. the archived one-time initialization (reads the two files,
      !    allocates and fills the network)
      call pes_init()
      ! 3. the judgment thresholds (absent keys = both slots unbound)
      if (p%z_term_found .and. p%z_term <= 0.0d0) then
         write (vs, '(g0)') p%z_term
         call stop_box('h2ag111_load', 'Z_TERM must be positive (got '//trim(vs)//')')
      end if
      if (p%r_h2_found .and. p%r_h2 <= 0.0d0) then
         write (vs, '(g0)') p%r_h2
         call stop_box('h2ag111_load', 'R_H2 must be positive (got '//trim(vs)//')')
      end if
      z_term = p%z_term
      r_h2 = p%r_h2
      term_on = p%z_term_found
      class_on = p%r_h2_found
      loaded = .true.
   end subroutine h2ag111_load

   !------------------------------------------------------------------
   ! h2ag111_pes(q, mass, g, v) - the container's aggregate force entry:
   !                 the NN surface's energy (POT0) and the authors'
   !                 central-difference gradient (DPESHON). The delivery
   !                 contract fully assigns g and v; mass is unused (a
   !                 pure potential over a rigid surface)
   !------------------------------------------------------------------
   subroutine h2ag111_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A] (all dof, flattened)
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A] (the interface folds F = -g)
      real(8), intent(out) :: v        ! V [eV]
      if (size(q) /= 6 .or. size(mass) /= 2) then
         call stop_box('h2ag111_pes', 'the H2/Ag(111) container needs exactly'// &
                       ' 2 dynamic atoms (6 dofs) - the system folder must define'// &
                       ' the two H atoms only (the Ag(111) slab is rigid inside the PES)')
      end if
      call pot0(2, q, v)
      call dpeshon(2, q, g)
   end subroutine h2ag111_pes

   !------------------------------------------------------------------
   ! h2ag111_term(sta) - the termination judgment: BOTH atoms above
   !                 Z_TERM (the molecule desorbed - the scattering is
   !                 over). Reached only when the drag-in bound the slot
   !------------------------------------------------------------------
   logical function h2ag111_term(sta)
      type(state_t), intent(in) :: sta  ! physical state (read-only snapshot)
      h2ag111_term = .false.
      if (.not. term_on) return
      h2ag111_term = (sta%q(3) > z_term) .and. (sta%q(6) > z_term)
   end function h2ag111_term

   !------------------------------------------------------------------
   ! h2ag111_classify(fin) - the channel judgment: r(H1,H2) below R_H2
   !                 names 'scattering' (molecule intact), else
   !                 'dissociation' (COM-frame archive; bond lengths
   !                 are translation invariant)
   !------------------------------------------------------------------
   function h2ag111_classify(fin) result(name)
      type(fin_t), intent(in) :: fin  ! final-state quantities (read-only)
      character(len=32) :: name       ! the outcome name [-]
      real(8) :: rhh
      if (.not. allocated(fin%q_fin)) then
         call stop_box('h2ag111_classify', 'the final coordinate archive is not'// &
                       ' allocated (the classification reads fin%q_fin)')
      end if
      rhh = sqrt(dot_product(fin%q_fin(4:6) - fin%q_fin(1:3), fin%q_fin(4:6) - fin%q_fin(1:3)))
      if (rhh < r_h2) then
         name = 'scattering'
      else
         name = 'dissociation'
      end if
   end function h2ag111_classify

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
end module h2ag111_box
