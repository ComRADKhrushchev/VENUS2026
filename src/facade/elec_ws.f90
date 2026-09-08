!=====================================================================
! elec_ws.f90 - L2 facade-internal wiring contract: electronic workspace ws_t (the container of
!               representation constructs - filled at force evaluation, consumed by elec_prop /
!               mqc_statistic; one integration-step lifetime spanning the three entries, E14)
! Reference: docs/plans/K2.md E10 (user adjudication 2026-09-07: interface invariance / black-box
!       property / thing-in-itself principle - the three axioms);
!       docs/superpowers/specs/2026-09-08-facade-e14-design.md §3 (E14, 2026-09-08: four-entry
!       contract + ws holding - the B/C axioms, ws lifetime, and request channel revised);
!       docs/design_spec_completeness.md §3.2-2 (two axes)
! Design axioms (enacted by E10; B and C amended by E14):
!   A Interface invariance - class-A/class-B member signatures are invariant across systems: the signature
!     carries only the container ws_t, never shape arrays like h_ij/d_ij (surface count and representation
!     stay out of the signature);
!     (see previous line)
!   B Black-box property - class A + class B together form one system black box; outside the box
!     (driver/propagator/output) only the four signatures of spec §3 are visible (force_eval /
!     elec_prop / mqc_statistic / termination via predicate); ws never appears in any of them;
!     (see previous line)
!     (see previous line)
!   C Thing-in-itself principle - representation constructs such as H_ij/d_ij/V_i have no
!     thing-in-itself status (they are "for us" quantities that change with the representation choice):
!     they do not enter state_t, do not leave the box, and live for ONE INTEGRATION STEP (E14: extended
!     from one evaluation) — filled at force evaluation, consumed by elec_prop / mqc_statistic,
!     reset to r_empty at step end. Holder: a module-private instance of elec_facade (still
!     invisible outside the box; declared in the box, dies in the box).
!     (see previous line)
!     (see previous line)
! This module is the dummy-type contract visible to members (a Fortran constraint: signature dummy types
! must be visible to the member modules) - dependency direction: force_facade -> elec_facade -> elec_ws;
! the member modules (pes_*/elec_*/dist_*) use this module for type visibility only (L3->L2 downward
! dependency is legal); this module itself has zero dependencies (only implicit none).
! (see previous line)
! (see previous line)
!=====================================================================
module elec_ws
   implicit none
   private
   public :: ws_t
   public :: r_empty, r_scalar, r_diab, r_adiab

   ! ws fill-mode values (the repr tag - representation demoted from an interface fork to a data tag, E10)
   integer, parameter :: r_empty  = 0  ! empty / not assembled (reset at each step end - E14 one-step lifetime)
   integer, parameter :: r_scalar = 1  ! single-surface scalar mode (g_add/v_add - plain additive members)
   integer, parameter :: r_diab   = 2  ! multi-surface diabatic mode (h_ij/dh_ij; natural for IESH)
   integer, parameter :: r_adiab  = 3  ! multi-surface adiabatic mode (v_surf/g_surf/d_ij)

   ! Electronic workspace: the box-private container of intermediates (held as a module-private
   ! instance of elec_facade, E14; members read/write it through the signature or the facade-internal
   ! channel; invisible outside the box)
   type :: ws_t
      integer :: repr = r_empty   ! fill mode (r_*: self-declared by the filling member; the facade and
                                  ! class-B members interpret the container by it - representation is data, not signature)
      integer, allocatable :: i_want(:)  ! request field - the occupation-set request (the adiab-mode
                                  ! delivery request: the set of currently occupied single-particle
                                  ! states whose coupling rows are requested; length 1 = the
                                  ! single-active-surface degenerate form; written by force_eval from
                                  ! sta%occ, read by the member. The occupation set stays out of the
                                  ! class-A signature and is conveyed through this field)
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
      real(8), allocatable :: d_ij(:,:)   ! derivative couplings <phi_i|grad_q|phi_j> of the requested
                                          ! i_want rows [1/Å] (ndof x n_surf per row - one row per
                                          ! occupied state, E14; O(ndof*n_surf) per row rather than
                                          ! all pairs; gap-weighted variants are left to the
                                          ! member/consumer - the contract accepts plain d)
   end type ws_t
contains
   ! This module is a pure type/constant contract with no routines (type lifecycle management belongs
   ! to elec_facade, E14: holder of the module-private instance + repr=r_empty reset at step end;
   ! allocate/reuse containers per mode - finalized in C2 together with the facade concretization)
   ! (see previous line)
end module elec_ws
