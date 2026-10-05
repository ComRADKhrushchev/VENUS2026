!=====================================================================
! reg_samp_ebk.f90 - register the 'ebk' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the quantum numbers arrive as buffered
!   member keys (N_VIB / N_ROT - field-default (0,0) when absent, the
!   ro-vibrational ground state; J_ROT stays the rotation member's key),
!   the spectroscopic constants arrive as assembly data through the
!   spectrum buffer's diatomic record (container-legislated, Herzberg
!   constants of the record - one shared set serves every carrier, the
!   equal-treatment statute; pulled at the first ebk carrier's
!   composition), and the moment of inertia is derived from that
!   carrier's buffered equilibrium geometry (mu*r_eq^2, the dist_j
!   convention - never duplicated as an input). An unselected member
!   skips.
!=====================================================================
module reg_samp_ebk
   use sampler,        only: samp_reg, samp_code
   use config,         only: reactants
   use config_atoms,   only: list_atoms
   use input,          only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,      only: seam_int
   use spectrum_export, only: spectrum_frag_diatomic
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
   ! reg_samp_ebk_init() - assembly seam: buffered quantum keys ->
   !                        constants pull at the first carrier ->
   !                        derived inertia -> member init
   !------------------------------------------------------------------
   subroutine reg_samp_ebk_init()
      type(ebk_params_t) :: p            ! member-parameter receiver (field defaults)
      character(len=512) :: sval         ! one buffered raw value
      real(8) :: m1, m2, req, dv(3)
      integer :: i, code, first
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('ebk')
      if (count(reactants%dist_scheme == code) == 0) return
      ! 1. the quantum numbers (absent keys keep the (0,0) field defaults)
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('N_VIB'); call seam_int(buffer_key(i), sval, buffer_line(i), p%n_vib)
         case ('N_ROT'); call seam_int(buffer_key(i), sval, buffer_line(i), p%j_rot)
         end select
      end do
      ! 2. the constants: pulled at the FIRST ebk carrier's composition
      !    (one shared set - equal treatment)
      first = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code) then
            first = i
            exit
         end if
      end do
      call spectrum_frag_diatomic(first, p%w_e, p%w_ex_e, p%b_rot)
      ! 3. the moment of inertia from the same carrier's buffered
      !    geometry (mu*r_eq^2)
      associate (fr => list_atoms%frag(first))
         m1 = list_atoms%mass(fr%list(1))
         m2 = list_atoms%mass(fr%list(2))
         dv = fr%qz(4:6) - fr%qz(1:3)
         req = norm2(dv)
      end associate
      p%ai_rot = (m1*m2/(m1 + m2))*req*req
      call ebk_init(p)
   end subroutine reg_samp_ebk_init
end module reg_samp_ebk
