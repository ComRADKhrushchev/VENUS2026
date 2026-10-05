!=====================================================================
! container_sched.f90 - the container assembly schedule: the nine
!   wiring slots, their lifecycle order, and the initialization-
!   assembly guard - written ONCE here so a container's reg file
!   carries content only (its wiring slots), never choreography
! Design:
!   container_sched_bind fills the slots from a container's reg file
!   (one handover; slots are single - one container per process, a
!   second fill of any slot aborts, mirroring container_bind_pes).
!   A null slot is a legal skip. container_sched_keys runs the key
!   vocabulary slot before the parse; container_sched_all runs the
!   unconditional slots at registration; container_sched_init gates on
!   buffered input (the guard, once) then runs the initialization and
!   output slots in fixed order. All conditional logic (key gates,
!   method selection) lives inside the container's own slot bodies -
!   this module knows time, not content.
!=====================================================================
module container_sched
   use config, only: reactants
   implicit none
   private
   public :: container_sched_bind, container_sched_keys, container_sched_all, &
             container_sched_init

   abstract interface
      subroutine slot_i               ! every slot body is argument-less
      end subroutine
   end interface

   ! The nine wiring slots: slot 1 input (keys / load), slot 2 regular
   ! (force / method / export), slot 3 initialization (term / arm),
   ! slot 4 output (classify / columns)
   type :: sched_slots_t
      procedure(slot_i), pointer, nopass :: keys     => null() ! declare the key vocabulary (before the parse)
      procedure(slot_i), pointer, nopass :: load     => null() ! pull the buffered keys + load (at the seam)
      procedure(slot_i), pointer, nopass :: force    => null() ! bind the force slot (registration)
      procedure(slot_i), pointer, nopass :: method   => null() ! register the method rows (registration)
      procedure(slot_i), pointer, nopass :: export   => null() ! registration-time exports (spectrum)
      procedure(slot_i), pointer, nopass :: term     => null() ! bind the termination slot (seam, key-gated)
      procedure(slot_i), pointer, nopass :: arm      => null() ! arm the method package (seam, selection-gated)
      procedure(slot_i), pointer, nopass :: classify => null() ! bind the classification slot (seam, key-gated)
      procedure(slot_i), pointer, nopass :: columns  => null() ! register the recorder columns (seam, selection-gated)
   end type
   type(sched_slots_t) :: w           ! the slots (module-private, run-lifetime)

contains
   !------------------------------------------------------------------
   ! container_sched_bind(...) - the container's single handover: fill
   !                 every present slot (the keyword list is the
   !                 author's full menu); an already-filled slot aborts
   !------------------------------------------------------------------
   subroutine container_sched_bind(keys, load, force, method, export, term, &
                                    arm, classify, columns)
      procedure(slot_i), optional :: keys, load, force, method, export, term, &
                                     arm, classify, columns ! the container's slot bodies
      if (present(keys))     call bind_one('keys',     w%keys,     keys)
      if (present(load))     call bind_one('load',     w%load,     load)
      if (present(force))    call bind_one('force',    w%force,    force)
      if (present(method))   call bind_one('method',   w%method,   method)
      if (present(export))   call bind_one('export',   w%export,   export)
      if (present(term))     call bind_one('term',     w%term,     term)
      if (present(arm))      call bind_one('arm',      w%arm,      arm)
      if (present(classify)) call bind_one('classify', w%classify, classify)
      if (present(columns))  call bind_one('columns',  w%columns,  columns)
   end subroutine container_sched_bind

   !------------------------------------------------------------------
   ! container_sched_keys() - run the key-vocabulary slot BEFORE the
   !                 parse (the grammar admits the container key names)
   !------------------------------------------------------------------
   subroutine container_sched_keys()
      if (associated(w%keys)) call w%keys()
   end subroutine container_sched_keys

   !------------------------------------------------------------------
   ! container_sched_all() - run the unconditional slots at
   !                 registration (nothing staged is required:
   !                 force slot, method rows, registration-time exports)
   !------------------------------------------------------------------
   subroutine container_sched_all()
      if (associated(w%force))  call w%force()
      if (associated(w%method)) call w%method()
      if (associated(w%export)) call w%export()
   end subroutine container_sched_all

   !------------------------------------------------------------------
   ! container_sched_init() - the initialization assembly: the guard
   !                 once (no buffered input -> no seam slot runs),
   !                 then the seam slots in fixed order
   !------------------------------------------------------------------
   subroutine container_sched_init()
      ! 1. the guard: nothing buffered -> skip (no load, no bind, no arm)
      if (len_trim(reactants%system_dir) == 0) return
      ! 2. the seam slots in fixed order: parameter load, termination
      !    slot, classification slot, method arming, recorder columns
      !    (column registration order = the slot body's call order)
      if (associated(w%load))     call w%load()
      if (associated(w%term))     call w%term()
      if (associated(w%classify)) call w%classify()
      if (associated(w%arm))      call w%arm()
      if (associated(w%columns))  call w%columns()
   end subroutine container_sched_init

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! bind_one(name, slot, proc) - fill one slot; an already-filled slot
   ! aborts (one container per process)
   subroutine bind_one(name, slot, proc)
      character(len=*), intent(in) :: name   ! the slot's name (abort payload)
      procedure(slot_i), pointer :: slot     ! the slot to fill
      procedure(slot_i) :: proc              ! the container's slot body
      if (associated(slot)) then
         call stop_cs('container_sched_bind', 'the "'//trim(name)//'" slot is already '// &
                      'filled (a second container handed its slots over - one container '// &
                      'per process)')
      end if
      slot => proc
   end subroutine bind_one

   ! stop_cs(who, msg) - the named-abort channel (name the caller, print the named
   ! message, STOP 1)
   subroutine stop_cs(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': schedule error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the schedule is not usable)'
      stop 1
   end subroutine stop_cs
end module container_sched
