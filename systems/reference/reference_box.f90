!=====================================================================
! reference_box.f90 - isotropic harmonic well: the reference container's
!                     concrete Hamiltonian + its termination and
!                     classification judgments
! REFERENCE MATERIAL (pure mathematics, zero system knowledge): the
!   potential V(q) = K/2 * |q|^2 treats every Cartesian degree of freedom
!   alike - any system folder is interchangeable, and only the optional
!   parameter file below changes behavior. This is not a physical model of
!   anything; it exists so the whole program can be driven end to end
!   against a closed-form force law whose exact dynamics are known.
! Container slots served (the drag-in file binds them):
!   PES           reference_pes - the composed kernel g = K*q [eV/A], the
!                 GRADIENT of the well (the facade folds the force
!                 F = -g), v = K/2*|q|^2 [eV] (both through
!                 harm_vg of reference_pot - this module carries the
!                 constant K and hands it to the member through
!                 harm_set_constant; the class-A unit contract:
!                 deliveries in eV - the eV-to-internal conversion
!                 belongs to the force facade alone, so K is carried in
!                 eV/A^2 here)
!   termination   reference_term - any |q_i| > Q_TERM (the coordinate-
!                 threshold form; bound only when Q_TERM was provided)
!   classification reference_classify - e_rel above/below E_SPLIT ->
!                 'scattering'/'capture' (bound only when E_SPLIT was
!                 provided; e_rel arrives in kcal/mol, so E_SPLIT is a
!                 kcal/mol number)
! Parameters: reference_pars.txt in the RUN directory (cwd-relative - the
!   same material location the input file reads from; the container
!   parameter file is RUN material, not system-definition material: the
!   system-folder scan treats every file there as a reactant definition
!   and loud-stops on any unrecognized extension, so a parameter file
!   cannot live beside the .xyz files). One KEYWORD = VALUE per line,
!   '#' starts a comment, blank lines ignored. An ABSENT file means all
!   defaults and both optional slots unbound (the safety-net-only,
!   no-classification run). Keys:
!     K_HARM  = <v>  force constant [eV/A^2]        (default 0.04)
!     Q_TERM  = <v>  termination threshold [A]      (absent -> unbound)
!     E_SPLIT = <v>  classification split [kcal/mol] (absent -> unbound)
!   A duplicated keyword, an unknown keyword, an unreadable value or a
!   non-positive K_HARM/Q_TERM loud-stops naming the file and the keyword
!   (the parameter file is hand-written material - a typo must not
!   silently run a different container).
! Units: see per-key list; energies in eV on the class-A channel, the
!   split in kcal/mol on the final-state channel.
!=====================================================================
module reference_box
   use state,       only: state_t
   use final_state, only: fin_t
   use reference_pot, only: harm_vg, harm_set_constant
   implicit none
   private
   public :: reference_load, reference_pes, reference_term, reference_classify
   public :: reference_term_on, reference_class_on
   public :: reference_k, reference_k_default

   real(8), parameter :: k_default = 0.04d0 ! default force constant [eV/A^2]
   real(8) :: k_harm  = k_default ! force constant [eV/A^2]
   real(8) :: q_term  = 0.0d0     ! termination threshold [A]
   real(8) :: e_split = 0.0d0     ! classification split [kcal/mol]
   logical :: term_on  = .false.  ! Q_TERM was provided (the slot may bind)
   logical :: class_on = .false.  ! E_SPLIT was provided (the slot may bind)
   logical :: loaded   = .false.  ! reference_load already ran (double-load guard)
contains
   !------------------------------------------------------------------
   ! reference_load() - read the container parameter file reference_pars.txt
   !                 from the run directory (absent = defaults, optional
   !                 slots unbound; present = parsed, every defect a named
   !                 loud-stop). Runs once, at the assembly seam. Whatever
   !                 the outcome, the member constant is handed to the
   !                 potential module before returning (the parsed value,
   !                 or the default when no file / no key)
   subroutine reference_load()
      character(len=512) :: path, line, key, val
      character(len=64)  :: fname
      logical :: exist
      integer :: u, ios, p
      logical :: got_k, got_q, got_e
      if (loaded) then
         call stop_box('reference_load', 'the container parameters were already'// &
                       ' loaded (a second load - an assembly-choreography error, formal)')
      end if
      loaded = .true.
      fname = 'reference_pars.txt'
      path = fname
      inquire (file=trim(path), exist=exist)
      if (.not. exist) then                   ! absent file: defaults, slots unbound
         call harm_set_constant(k_harm)
         return
      end if
      got_k = .false.
      got_q = .false.
      got_e = .false.
      open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call stop_box('reference_load', 'cannot open '//trim(path))
      end if
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         p = index(line, '#')
         if (p > 0) line = line(:p-1)       ! strip the comment tail
         p = index(line, '=')
         if (p == 0) then
            if (len_trim(adjustl(line)) > 0) then
               call stop_box('reference_load', 'line without "=" in '//trim(path)//': "'// &
                             trim(adjustl(line))//'" (grammar: KEYWORD = VALUE)')
            end if
            cycle                          ! blank line
         end if
         key = adjustl(line(:p-1))
         val = adjustl(line(p+1:))
         if (len_trim(key) == 0 .or. len_trim(val) == 0) then
            call stop_box('reference_load', 'empty keyword or value in '//trim(path)//': "'// &
                          trim(line)//'"')
         end if
         select case (trim(key))
         case ('K_HARM')
            if (got_k) call stop_dup('K_HARM')
            read (val, *, iostat=ios) k_harm
            if (ios /= 0) call stop_val('K_HARM', val)
            if (k_harm <= 0.0d0) then
               call stop_box('reference_load', 'K_HARM must be positive in '//trim(path)// &
                             ' (got '//trim(val)//' eV/A^2 - a non-positive constant is'// &
                             ' not a well)')
            end if
            got_k = .true.
         case ('Q_TERM')
            if (got_q) call stop_dup('Q_TERM')
            read (val, *, iostat=ios) q_term
            if (ios /= 0) call stop_val('Q_TERM', val)
            if (q_term <= 0.0d0) then
               call stop_box('reference_load', 'Q_TERM must be positive in '//trim(path)// &
                             ' (got '//trim(val)//' A)')
            end if
            got_q = .true.
            term_on = .true.
         case ('E_SPLIT')
            if (got_e) call stop_dup('E_SPLIT')
            read (val, *, iostat=ios) e_split
            if (ios /= 0) call stop_val('E_SPLIT', val)
            got_e = .true.
            class_on = .true.
         case default
            call stop_box('reference_load', 'unknown keyword "'//trim(key)//'" in '// &
                          trim(path)//' (the closed key set: K_HARM / Q_TERM / E_SPLIT)')
         end select
      end do
      close (u)
      call harm_set_constant(k_harm)          ! the member init hand-off (always)
   contains
      ! stop_dup(k) / stop_val(k, v) - the two message shapes of the parse loop
      subroutine stop_dup(k)
         character(len=*), intent(in) :: k
         call stop_box('reference_load', 'keyword '//trim(k)//' appears twice in '//trim(path))
      end subroutine stop_dup
      subroutine stop_val(k, v)
         character(len=*), intent(in) :: k, v
         call stop_box('reference_load', 'value of '//trim(k)//' unreadable in '//trim(path)// &
                       ': "'//trim(v)//'"')
      end subroutine stop_val
   end subroutine reference_load

   !------------------------------------------------------------------
   ! reference_pes(q, mass, g, v) - the container aggregate class-A entry: the
   !                 isotropic well, computed by the member harm_vg (which
   !                 carries the constant this module staged through
   !                 harm_set_constant). The delivery contract fully assigns
   !                 g and v; mass is unused (a pure potential)
   subroutine reference_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A] (all dof, flattened)
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A] (the facade folds F = -g)
      real(8), intent(out) :: v        ! V [eV]
      call harm_vg(size(q)/3, q, v, g)
   end subroutine reference_pes

   !------------------------------------------------------------------
   ! reference_term(sta) - the termination judgment: any coordinate leaves
   !                 the [-Q_TERM, Q_TERM] box (read-only state). Reached
   !                 only when the drag-in bound the slot
   logical function reference_term(sta)
      type(state_t), intent(in) :: sta  ! physical state (read-only snapshot)
      reference_term = .false.
      if (.not. term_on) return
      reference_term = any(abs(sta%q) > q_term)
   end function reference_term

   !------------------------------------------------------------------
   ! reference_classify(fin) - the classification judgment: relative
   !                 translational energy above the split -> 'scattering',
   !                 otherwise 'capture' (read-only final-state quantities)
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
   pure logical function reference_term_on()
      reference_term_on = term_on
   end function reference_term_on

   pure logical function reference_class_on()
      reference_class_on = class_on
   end function reference_class_on

   !------------------------------------------------------------------
   ! reference_k() / reference_k_default() - the constant getters (the
   !                 assembly verification face: what the file staged /
   !                 what the default is)
   pure function reference_k() result(k)
      real(8) :: k                          ! the staged force constant [eV/A^2]
      k = k_harm
   end function reference_k

   pure function reference_k_default() result(k)
      real(8) :: k                          ! the default force constant [eV/A^2]
      k = k_default
   end function reference_k_default

   !------------------------------------------------------------------
   ! stop_box(who, msg) - the loud-stop channel (name the caller, print the
   !                 named message, STOP 1 - the sibling modules' semantics)
   subroutine stop_box(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': container parameter error: '//trim(msg)
      write (0, '(a)') trim(who)//': loud-stop (the container is not usable)'
      stop 1
   end subroutine stop_box
end module reference_box
