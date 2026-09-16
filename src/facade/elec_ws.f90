!=====================================================================
! elec_ws.f90 - electronic workspace ws_t: the container of
!   representation constructs filled at force evaluation, consumed by
!   the electronic phases
! Design:
!   ws is box-private: never in state_t, never a public-signature
!   argument, one integration-step lifetime (filled by force_eval,
!   consumed by elec_prop / mqc_statistic, reset to r_empty at step
!   end). Member signatures carry only ws_t, never shape arrays like
!   h_ij/d_ij; representation is a data tag, not a signature fork.
!   Zero-dependency dummy-type contract; lifecycle belongs to elec_facade.
!=====================================================================
module elec_ws
   implicit none
   private
   public :: ws_t
   public :: r_empty, r_scalar, r_diab, r_adiab

   ! ws fill-mode values (the repr tag - representation as data, not a signature fork)
   integer, parameter :: r_empty  = 0  ! empty / not assembled (reset at each step end)
   integer, parameter :: r_scalar = 1  ! single-surface scalar mode (g_add/v_add - plain additive members)
   integer, parameter :: r_diab   = 2  ! multi-surface diabatic mode (h_ij/dh_ij; natural for IESH)
   integer, parameter :: r_adiab  = 3  ! multi-surface adiabatic mode (v_surf/g_surf/d_ij)

   ! Electronic workspace: the box-private container of intermediates (held as a
   ! module-private instance of elec_facade; invisible outside the box)
   type :: ws_t
      integer :: repr = r_empty   ! fill mode (r_*: self-declared by the filling member; the facade and
                                  ! class-B members interpret the container by it - representation is data, not signature)
      integer, allocatable :: i_want(:)  ! request field - the occupation-set request (the set of
                                  ! currently occupied single-particle states whose coupling rows
                                  ! are requested; length 1 = the single-active-surface degenerate
                                  ! form; written by force_eval from the held occ, read by the
                                  ! member; the occupation set stays out of the class-A signature)
      ! -- scalar-mode container --
      real(8), allocatable :: g_add(:)  ! force kernel -dV/dq [eV/Å] (an equivalent force kernel for stochastic laws)
      real(8)              :: v_add = 0.0d0  ! energy contribution [eV] (default-initialized:
                                             ! a construct-only delivery leaves it at zero by
                                             ! construction, not by luck)
      ! -- diab-mode container --
      real(8), allocatable :: h_ij(:,:)   ! diabatic Hamiltonian matrix [eV] (symmetric of order n_surf;
                                          ! the member delivers the full matrix; symmetry is the member's own affair)
      real(8), allocatable :: dh_ij(:,:,:) ! matrix-element gradients [eV/Å] (ndof x n_surf^2;
                                          ! FULL-MATRIX delivery: the facade's ordered-pair fold
                                          ! reads BOTH triangles - symmetry is not assumed on read)
      ! -- adiab-mode container --
      real(8), allocatable :: v_surf(:)   ! per-surface energies [eV]
      real(8), allocatable :: g_surf(:,:) ! per-surface force kernels [eV/Å] (ndof x n_surf)
      real(8), allocatable :: d_ij(:,:)   ! derivative couplings <phi_i|grad_q|phi_j> of the requested
                                          ! i_want rows [1/Å] (ndof x n_surf per row - one row per
                                          ! occupied state; gap-weighted variants are left to the
                                          ! member/consumer - the contract accepts plain d)
   end type ws_t
contains
   ! Pure type/constant contract with no routines (lifecycle management belongs to
   ! elec_facade: the holder + the ws_handle alias + the full-clear reset at step end)
end module elec_ws
