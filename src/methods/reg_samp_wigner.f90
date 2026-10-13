!=====================================================================
! reg_samp_wigner.f90 - register the 'wigner' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the carrier's bracket parameters
!   (E_SCALE - one row per carrier, the 2026-10-13 bracket ruling; an
!   absent name keeps the member's field default) -> the mode-table
!   pull -> member init. The tables are assembly data, not input keys -
!   a selected member aborts with a named error at its own unallocated
!   check.
!=====================================================================
module reg_samp_wigner
   use sampler,        only: samp_reg, samp_code
   use config,         only: reactants
   use input,          only: par_real
   use spectrum_interface, only: spectrum_frag_w
   use samp_wigner,    only: wigner_draw, wigner_realize, wigner_params_t, wigner_init
   implicit none
   private
   public :: reg_samp_wigner_all
   public :: reg_samp_wigner_init
contains
   !------------------------------------------------------------------
   ! reg_samp_wigner_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_wigner_all()
      call samp_reg('wigner', wigner_draw, wigner_realize)
   end subroutine reg_samp_wigner_all

   !------------------------------------------------------------------
   ! reg_samp_wigner_init() - assembly seam: per-carrier bracket keys ->
   !                          the mode-table pull -> member init. The
   !                          tables are assembly data: pulled from the
   !                          SPECTRUM_SOURCE-selected producer for
   !                          EVERY carrier fragment, one ragged table
   !                          per carrier
   !------------------------------------------------------------------
   subroutine reg_samp_wigner_init()
      type(wigner_params_t) :: p     ! member-parameter receiver (field defaults)
      integer :: i, k, code, n_carry
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('wigner')
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) return
      allocate (p%e_scale(n_carry))
      p%e_scale = 1.0d0                  ! the member's field default (the ground
                                          ! state; an absent name keeps it)
      allocate (p%tbl(n_carry))
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         k = k + 1
         call par_real('WIGNER', i, 'E_SCALE', p%e_scale(k))
         call spectrum_frag_w(i, p%tbl(k)%w, p%tbl(k)%c)
      end do
      call wigner_init(p)
   end subroutine reg_samp_wigner_init
end module reg_samp_wigner
