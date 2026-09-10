! state.f90 - physical-state variables (Q/P/s/F) container: held by the driver, passed by argument
!             (K1 minimal set;
!              E9 multi-surface upgrade: elec_t gains electronic amplitudes and a hop counter -
!              the data slot of the nonadiabatic channel;
!              E10 representation statute: amp index semantics defined by the class-B method
!              member; representation constructs do not enter state;
!              E14 dual-container revision (spec §4): amp upgraded to the coefficient matrix
!              a(:,:) with rho(:,:) added (symmetric ontic candidates - no subordination
!              legislated); i_surf generalized to the occupation vector occ(:), relocated
!              with n_hop to state_t top level - classical-side statistical fields;
!              E16 container-purity revision (review-gate ruling 2026-09-10): the cached
!              dimension labels natoms/ndof deleted (size-derived: ndof = size(q),
!              natoms = size(roster%mass)) and the statistical fields occ/n_hop relocated
!              to the elec_facade statistics holder (sampled constructs have no
!              thing-in-itself status — the E10 axiom family extended) - state_t is
!              purely Q/P/s/F/t)
module state
  implicit none
  private
  public :: state_t, elec_t, state_create, state_destroy

  ! E14 (user rulings 2026-09-08): (1) amp(:) upgraded to the coefficient matrix a(:,:);
  ! (2) rho(:,:) added — symmetric dual ontic containers, no subordination legislated
  ! (state legislates only the EXISTENCE of two allocatable complex matrices; which one a
  ! method propagates is its own declaration — amplitude-ODE methods on a, density-matrix-
  ! dynamics methods on rho, adiabatic/LZSH/MDEF allocate neither); (3) i_surf generalized
  ! to the occupation vector occ(:) — a CLASSICAL-SIDE STATISTICAL field (the discrete
  ! output of mqc_statistic's stochastic sampling of the populations; under the
  ! single-electron representation multiple states are occupied simultaneously, so a single
  ! surface index is not general) [E16: occ no longer lives in state_t — relocated to the
  ! elec_facade statistics holder; the generalization ruling itself stands]; (4) index-
  ! semantics statutes stay member-defined
  ! (E11/E12 route): a rows/columns = surface number / electron×state / configuration
  ! cardinality / mapping z = x+ip encoding — defined by the class-B method member.
  type :: elec_t          ! electronic state s (runtime): purely the quantum-state ontic layer -
                          ! the dual containers a/rho plus the dimension label n_surf (E14, spec §4)
     integer :: n_surf = 1   ! number of potential-energy surfaces (representation-statute extension:
                             ! under the single-electron representation read as the configuration-space
                             ! cardinality - semantics defined by the class-B method member, field
                             ! name unchanged)
     complex(8), allocatable :: a(:,:)  ! coefficient matrix (propagation container; upgraded
                                        ! from the complex vector amp)
     complex(8), allocatable :: rho(:,:)! density matrix (ontic for rho-native methods;
                                        ! derived cache for a-native methods — the mapping
                                        ! form is method knowledge, refreshed by sync_proc)
  end type

  type :: state_t            ! purely physical (E16): five members Q/P/s/F/t — the cached
                             ! dimension labels deleted (ndof = size(q), natoms =
                             ! size(roster%mass)) and the statistical fields relocated
                             ! to the elec_facade statistics holder (occ/n_hop: sampled
                             ! constructs, not things-in-themselves)
     real(8), allocatable :: q(:) ! coordinates [Å] (size = ndof = 3*natoms)
     real(8), allocatable :: p(:) ! momenta [amu·Å/(10 fs)] (size = ndof)
     real(8), allocatable :: f(:) ! forces [internal energy/Å] (size = ndof)
     real(8) :: t = 0.0d0     ! trajectory time [10 fs] (held by the driver)
     type(elec_t) :: s        ! electronic ontic state
  end type
contains
  ! Minimal implementation: allocate q/p/f by natoms (ndof = 3*natoms; E16: no dimension
  ! labels stored — consumers derive size(q); the degenerate occupation occ=[1] is no
  ! longer seeded here — it arrives with the method declaration at method_bind, the
  ! elec_facade statistics holder; a/rho seeded at the degenerate 1x1 shape when n_surf
  ! present - larger ontic containers (e.g. IESH n_e x n_s) are reallocated by the method
  ! member-init hook, blackbox_init step, spec §4 "Allocation timing"; n_surf defaults
  ! to 1 - the adiabatic degenerate point; on failure write to stderr and stop)
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
     allocate(st%q(3*natoms), st%p(3*natoms), st%f(3*natoms), stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_create: allocation of q/p/f failed, stat = ', stat
        stop 1
     end if
     st%q = 0.0d0
     st%p = 0.0d0
     st%f = 0.0d0
     st%s%n_surf = n_s
     if (present(n_surf)) then ! degenerate 1x1 ontic containers (E14 allocation timing:
                               ! member-init reallocates to the method's declared shape)
        allocate(st%s%a(1,1), st%s%rho(1,1), stat=stat)
        if (stat /= 0) then
           write(0,'(a,i0)') 'state_create: allocation of a/rho failed, stat = ', stat
           stop 1
        end if
        st%s%a = (0.0d0, 0.0d0)
        st%s%a(1,1) = (1.0d0, 0.0d0)   ! degenerate-point population 1 (well-defined at the
                                       ! adiabatic degenerate point; carries the old amp seed)
        st%s%rho = (0.0d0, 0.0d0)
        st%s%rho(1,1) = (1.0d0, 0.0d0) ! symmetric degenerate seed - no subordination legislated
     end if
  end subroutine state_create

  ! Minimal implementation: deallocate q/p/f/a/rho and reset the labels (a/rho are
  ! seeded together in state_create but independently reallocated by member init, E14 -
  ! each deallocation is guarded separately); on failure write to stderr and stop
  subroutine state_destroy(st)
     type(state_t), intent(inout) :: st
     integer :: stat
     deallocate(st%q, st%p, st%f, stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_destroy: deallocation of q/p/f failed, stat = ', stat
        stop 1
     end if
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
     st%s%n_surf = 1
  end subroutine state_destroy
end module
