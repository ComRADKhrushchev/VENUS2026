!=====================================================================
! state.f90 - physical-state container (Q/P/s/F/t), held by the driver and
!             passed by argument
! Design:
!   state_t is purely Q/P/s/F/t; dimensions are size-derived (ndof = size(q),
!   natoms = size(list_atoms%mass)) - no cached labels, no statistical fields (occ/n_hop
!   live in the elec_interface statistics holder). elec_t carries the dual ontic
!   containers a/rho, symmetric with no subordination prescribed (amplitude-ODE
!   methods propagate a, density-matrix methods rho, adiabatic allocates neither);
!   index semantics are defined by the electronic method member.
!=====================================================================
module state
  implicit none
  private
  public :: state_t, elec_t, state_create, state_destroy

  ! Electronic state s (runtime): purely the quantum-state ontic layer - the dual
  ! containers a/rho plus the dimension label n_surf. Index semantics of a/rho
  ! (surface number / electron x state / configuration cardinality / z = x+ip
  ! encoding) are defined by the electronic method member.
  type :: elec_t
     integer :: n_surf = 1   ! number of potential-energy surfaces (under the
                             ! single-electron representation read as the
                             ! configuration-space cardinality - semantics defined
                             ! by the electronic method member, field name unchanged)
     complex(8), allocatable :: a(:,:)  ! coefficient matrix (propagation container)
     complex(8), allocatable :: rho(:,:)! density matrix (ontic for rho-native methods;
                                        ! derived cache for a-native methods - the mapping
                                        ! form is method knowledge, refreshed by sync_proc)
  end type

  type :: state_t            ! purely physical: Q/P/s/F/t - dimensions size-derived,
                             ! statistical fields live in the elec_interface holder
     real(8), allocatable :: q(:) ! coordinates [Å] (size = ndof = 3*natoms)
     real(8), allocatable :: p(:) ! momenta [amu·Å/(10 fs)] (size = ndof)
     real(8), allocatable :: f(:) ! forces [internal energy/Å] (size = ndof)
     real(8) :: t = 0.0d0     ! trajectory time [10 fs] (held by the driver)
     type(elec_t) :: s        ! electronic ontic state
  end type
contains
  !------------------------------------------------------------------
  ! state_create(st, natoms [, n_surf]) - allocate q/p/f and seed the electronic
  !                                        ontic containers; stderr + stop on failure
  !------------------------------------------------------------------
  subroutine state_create(st, natoms, n_surf)
     type(state_t), intent(out) :: st
     integer, intent(in) :: natoms
     integer, intent(in), optional :: n_surf ! number of surfaces [-] (default 1;
                                              ! multi-surface takes config%electronic%n_surf)
     integer :: n_s, stat
     ! 1. validate n_surf
     n_s = 1
     if (present(n_surf)) n_s = n_surf
     if (n_s < 1) then
        write(0,'(a,i0)') 'state_create: illegal n_surf (<1): ', n_s
        stop 1
     end if
     ! 2. allocate and zero q/p/f
     allocate(st%q(3*natoms), st%p(3*natoms), st%f(3*natoms), stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_create: allocation of q/p/f failed, stat = ', stat
        stop 1
     end if
     st%q = 0.0d0
     st%p = 0.0d0
     st%f = 0.0d0
     ! 3. seed the electronic layer
     st%s%n_surf = n_s
     if (present(n_surf)) then ! degenerate 1x1 ontic containers (member-init
                               ! reallocates to the method's declared shape)
        allocate(st%s%a(1,1), st%s%rho(1,1), stat=stat)
        if (stat /= 0) then
           write(0,'(a,i0)') 'state_create: allocation of a/rho failed, stat = ', stat
           stop 1
        end if
        st%s%a = (0.0d0, 0.0d0)
        st%s%a(1,1) = (1.0d0, 0.0d0)   ! degenerate-point population 1 (well-defined at the
                                       ! adiabatic degenerate point)
        st%s%rho = (0.0d0, 0.0d0)
        st%s%rho(1,1) = (1.0d0, 0.0d0) ! symmetric degenerate seed - no subordination prescribed
     end if
  end subroutine state_create

  !------------------------------------------------------------------
  ! state_destroy(st) - deallocate q/p/f/a/rho and reset the label (each
  !                     deallocation guarded; stderr + stop on failure)
  !------------------------------------------------------------------
  subroutine state_destroy(st)
     type(state_t), intent(inout) :: st
     integer :: stat
     ! 1. deallocate q/p/f
     deallocate(st%q, st%p, st%f, stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_destroy: deallocation of q/p/f failed, stat = ', stat
        stop 1
     end if
     ! 2. deallocate a and rho independently (member init may have reallocated either)
     if (allocated(st%s%a)) then
        deallocate(st%s%a, stat=stat)
        if (stat /= 0) then
           write(0,'(a,i0)') 'state_destroy: deallocation of a failed, stat = ', stat
           stop 1
        end if
     end if
     if (allocated(st%s%rho)) then
        deallocate(st%s%rho, stat=stat)
        if (stat /= 0) then
           write(0,'(a,i0)') 'state_destroy: deallocation of rho failed, stat = ', stat
           stop 1
        end if
     end if
     ! 3. reset the label
     st%s%n_surf = 1
  end subroutine state_destroy
end module
