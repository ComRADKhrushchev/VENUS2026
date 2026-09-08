! state.f90 - physical-state variables (Q/P/s/F) container: held by the driver, passed by argument
!             (K1 minimal set;
!              E9 multi-surface upgrade: elec_t gains electronic amplitudes and a hop counter -
!              the data slot of the nonadiabatic channel;
!              E10 representation statute: amp index semantics defined by the class-B method
!              member; representation constructs do not enter state)
module state
  implicit none
  private
  public :: elec_t, state_t, state_create, state_destroy
  type :: elec_t          ! electronic state s (runtime): adiabatic = degenerate point (single surface)
     integer :: n_surf = 1   ! number of potential-energy surfaces
     integer :: i_surf = 1   ! current surface index
     ! E10 representation statute (enacted in comments, field names untouched): amp only
     ! promises "a complex array whose |.|^2 can be read as
     ! populations"; the index semantics are defined by the class-B method member - surface
     ! representation i = surface number; single-electron representation
     ! i = flattened electron x state tensor product (fermionic antisymmetrization is the
     ! member's own affair); configuration representation i = configuration index.
                                       ! Representation constructs (H_ij/d_ij/V_i etc.) have no thing-in-itself status and do
     ! not enter state_t (E10 axiom C
     ! - only the facade workspace ws remains, see facade/elec_ws.f90; one-step lifetime)
     complex(8), allocatable :: amp(:) ! electronic amplitudes c_i(t) (|c_i|^2 = populations;
     ! adiabatic degenerate point = single element,
     ! identically 1 - keeping the degenerate point well-defined)
     integer :: n_hop = 0    ! accumulated hop count [count] (nonadiabatic observable; always 0 for adiabatic)
  end type
  type :: state_t
     integer :: natoms = 0   ! number of atoms (fixed at runtime from the roster)
     integer :: ndof   = 0   ! number of degrees of freedom = 3*natoms
     real(8), allocatable :: q(:) ! coordinates [Å]
     real(8), allocatable :: p(:) ! momenta [amu·Å/(10 fs)]
     real(8), allocatable :: f(:) ! forces [internal energy/Å]
     type(elec_t) :: s
     real(8) :: t = 0.0d0     ! trajectory time [10 fs] (held by the driver)
  end type
contains
  ! Minimal implementation: allocate q/p/f and the electronic amplitudes by natoms
  ! (ndof = 3*natoms; n_surf defaults to 1
  ! - the adiabatic degenerate point; on failure write to stderr and stop)
  subroutine state_create(st, natoms, n_surf)
     type(state_t), intent(out) :: st
     integer, intent(in) :: natoms
     integer, intent(in), optional :: n_surf ! number of surfaces [-] (default 1;
                                              ! multi-surface takes config%electronic%n_surf)
     integer :: n_s, stat
     n_s = 1
     if (present(n_surf)) n_s = n_surf
     if (n_s < 1) then
        write(0,'(a,i0)') 'state_create: illegal n_surf (<1): ', n_s
        stop 1
     end if
     st%natoms = natoms
     st%ndof   = 3*natoms
     allocate(st%q(st%ndof), st%p(st%ndof), st%f(st%ndof), st%s%amp(n_s), stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_create: allocation of q/p/f/amp failed, stat = ', stat
        stop 1
     end if
     st%q = 0.0d0
     st%p = 0.0d0
     st%f = 0.0d0
     st%s%n_surf = n_s
     st%s%i_surf = 1
     st%s%amp = (0.0d0, 0.0d0)
     st%s%amp(1) = (1.0d0, 0.0d0)  ! ground-surface population 1 (well-defined at the adiabatic degenerate point)
     st%s%n_hop = 0
  end subroutine state_create

  ! Minimal implementation: deallocate q/p/f/amp and reset the counters; on failure write
  ! to stderr and stop
  subroutine state_destroy(st)
     type(state_t), intent(inout) :: st
     integer :: stat
     deallocate(st%q, st%p, st%f, st%s%amp, stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_destroy: deallocation of q/p/f/amp failed, stat = ', stat
        stop 1
     end if
     st%natoms = 0
     st%ndof   = 0
     st%s%n_surf = 1
     st%s%i_surf = 1
     st%s%n_hop  = 0
  end subroutine state_destroy
end module
