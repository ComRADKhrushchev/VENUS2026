!=====================================================================
! reg_samp_ebk.f90 - register the 'ebk' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the quantum numbers arrive as the
!   carrier's bracket parameters (N_VIB / N_ROT - field-default (0,0)
!   when absent, the ro-vibrational ground state; one row per carrier -
!   the 2026-10-13 bracket ruling realized per fragment, two EBK
!   fragments may sit on different quantum states), the spectroscopic
!   constants arrive as assembly data through the spectrum seam's
!   diatomic record PER CARRIER (MANUAL spectrum-table DIAT - the
!   Herzberg constants of the record, literature data with no generic
!   derivation), and the moment of inertia is derived from that
!   carrier's buffered equilibrium geometry (mu*r_eq^2, the dist_j
!   convention - never duplicated as an input). An unselected member
!   skips.
!=====================================================================
module reg_samp_ebk
   use sampler,        only: samp_reg, samp_code
   use config,         only: reactants
   use config_atoms,   only: list_atoms
   use input,          only: par_int
   use spectrum_interface, only: spectrum_frag_diatomic
   use samp_ebk,       only: ebk_draw, ebk_realize, ebk_params_t, ebk_init
   implicit none
   private
   public :: reg_samp_ebk_all
   public :: reg_samp_ebk_init
contains
   !------------------------------------------------------------------
   ! reg_samp_ebk_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_ebk_all()
      call samp_reg('ebk', ebk_draw, ebk_realize)
   end subroutine reg_samp_ebk_all

   !------------------------------------------------------------------
   ! reg_samp_ebk_init() - assembly seam: per-carrier bracket quantum
   !                        keys -> per-carrier constants pull ->
   !                        derived inertia -> member init
   !------------------------------------------------------------------
   subroutine reg_samp_ebk_init()
      type(ebk_params_t) :: p            ! member-parameter receiver (field defaults)
      real(8) :: m1, m2, req, dv(3)
      integer :: i, kc, code
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('ebk')
      if (count(reactants%dist_scheme == code) == 0) return
      allocate (p%car(count(reactants%dist_scheme == code)))
      kc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         kc = kc + 1
         ! 1. the quantum numbers (absent bracket names keep the (0,0) field
         !    defaults - the ro-vibrational ground state)
         call par_int('EBK', i, 'N_VIB', p%car(kc)%n_vib)
         call par_int('EBK', i, 'N_ROT', p%car(kc)%j_rot)
         ! 2. the constants: pulled at THIS carrier's composition
         call spectrum_frag_diatomic(i, p%car(kc)%w_e, p%car(kc)%w_ex_e, p%car(kc)%b_rot)
         ! 3. the moment of inertia from the same carrier's buffered
         !    geometry (mu*r_eq^2)
         associate (fr => list_atoms%frag(i))
            m1 = list_atoms%mass(fr%list(1))
            m2 = list_atoms%mass(fr%list(2))
            dv = fr%qz(4:6) - fr%qz(1:3)
            req = norm2(dv)
         end associate
         p%car(kc)%ai_rot = (m1*m2/(m1 + m2))*req*req
      end do
      call ebk_init(p)
   end subroutine reg_samp_ebk_init
end module reg_samp_ebk
