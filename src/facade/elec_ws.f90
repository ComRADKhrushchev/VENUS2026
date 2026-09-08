!=====================================================================
! elec_ws.f90 - L2 facade-internal wiring contract: electronic workspace ws_t (the container of
!               representation constructs - filled by class-A members, consumed by class-B members;
!               one-step lifetime, rebuilt at every evaluation)
! Reference: docs/plans/K2.md E10 (user adjudication 2026-09-07: interface invariance / black-box
!       property / thing-in-itself principle - the three axioms); docs/design_spec_completeness.md §3.2-2 (two axes)
! Design axioms (enacted by E10):
!   A Interface invariance - class-A/class-B member signatures are invariant across systems: the signature
!     carries only the container ws_t, never shape arrays like h_ij/d_ij (surface count and representation
!     stay out of the signature);
!     (see previous line)
!   B Black-box property - class A + class B together form one system black box; outside the box
!     (driver/propagator/output) only (Q,P,s[,dt])->(F,E) is visible; ws never appears in a
!     facade_force signature outside the box;
!     (see previous line)
!     (see previous line)
!   C Thing-in-itself principle - representation constructs such as H_ij/d_ij/V_i have no
!     thing-in-itself status (they are "for us" quantities that change with the representation choice):
!     they do not enter state_t, do not leave the box, and live for one step (declared locally in
!     facade_force, dying inside the facade).
!     (see previous line)
!     (see previous line)
! This module is the dummy-type contract visible to members (a Fortran constraint: signature dummy types
! must be visible to the member modules) - dependency direction: force_facade and the
! interface_electronic/ members all use this module (L3->L2 downward dependency is legal); this module
! itself has zero dependencies (only implicit none).
! (see previous line)
! (see previous line)
!=====================================================================
module elec_ws
   implicit none
   private
   public :: ws_t
   public :: r_empty, r_scalar, r_diab, r_adiab

   ! ws fill-mode values (the repr tag - representation demoted from an interface fork to a data tag, E10)
   integer, parameter :: r_empty  = 0  ! empty / not assembled (facade resets before each evaluation)
   integer, parameter :: r_scalar = 1  ! single-surface scalar mode (g_add/v_add - plain additive members)
   integer, parameter :: r_diab   = 2  ! multi-surface diabatic mode (h_ij/dh_ij; natural for IESH)
   integer, parameter :: r_adiab  = 3  ! multi-surface adiabatic mode (v_surf/g_surf/d_ij)

   ! Electronic workspace: the facade-private container of intermediates (members read/write it via the
   ! signature; invisible outside the box)
   type :: ws_t
      integer :: repr = r_empty   ! fill mode (r_*: self-declared by the filling member; the facade and
                                  ! class-B members interpret the container by it - representation is data, not signature)
      integer :: i_want = 0       ! request field (written by the facade: a trimming request - in adiabatic
                                  ! mode, the i of the current-surface coupling row d_i,.; the member reads it.
                                  ! i_surf stays out of the class-A signature and is conveyed through this field)
      ! -- scalar-mode container --
      real(8), allocatable :: g_add(:)  ! force kernel -dV/dq [eV/Å] (an equivalent force kernel for stochastic laws)
      real(8)              :: v_add     ! energy contribution [eV]
      ! -- diab-mode container --
      real(8), allocatable :: h_ij(:,:)   ! diabatic Hamiltonian matrix [eV] (symmetric of order n_surf;
                                          ! the member delivers the full matrix; symmetry is the member's own affair)
      real(8), allocatable :: dh_ij(:,:,:) ! matrix-element gradients [eV/Å] (ndof x n_surf^2;
                                          ! compressed delivery for large state spaces is a grammar-versioning matter)
      ! -- adiab-mode container --
      real(8), allocatable :: v_surf(:)   ! per-surface energies [eV]
      real(8), allocatable :: g_surf(:,:) ! per-surface force kernels [eV/Å] (ndof x n_surf)
      real(8), allocatable :: d_ij(:,:)   ! current-surface derivative couplings <phi_i|grad_q|phi_j> [1/Å]
                                          ! (ndof x n_surf - delivered by the requested i_want row,
                                          ! O(ndof*n_surf) rather than all pairs; gap-weighted variants are
                                          ! left to the member/consumer - the contract accepts plain d)
   end type ws_t
contains
   ! This module is a pure type/constant contract with no routines (type lifecycle management belongs to
   ! force_facade: reset repr=r_empty before each evaluation, allocate/reuse containers per mode -
   ! finalized in C2 together with the facade concretization)
   ! (see previous line)
end module elec_ws
