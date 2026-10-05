!=====================================================================
! reg_samp_wigner.f90 - register the 'wigner' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: buffered input keys -> params receiver
!   (absent keys keep defaults) -> the mode-table pull -> member init.
!   An unselected member skips. The w_mode/c_mode tables are assembly
!   data, not input keys - a selected member aborts with a named error
!   at its own unallocated check.
!=====================================================================
module reg_samp_wigner
   use sampler,        only: samp_reg, samp_code
   use config,         only: reactants
   use input,          only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,      only: seam_real
   use spectrum_export, only: spectrum_frag_w
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
   ! reg_samp_wigner_init() - assembly seam: buffered keys -> the
   !                          mode-table pull -> member init. The tables
   !                          are assembly data: pulled from the spectrum
   !                          buffer for EVERY carrier fragment (a
   !                          composition-matched container export first,
   !                          the generic derivation otherwise), one
   !                          ragged table per carrier
   !------------------------------------------------------------------
   subroutine reg_samp_wigner_init()
      type(wigner_params_t) :: p     ! member-parameter receiver (field defaults)
      character(len=512) :: sval     ! one buffered raw value
      integer :: i, k, code, n_carry
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('wigner')
      n_carry = count(reactants%dist_scheme == code)
      if (n_carry == 0) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('E_SCALE'); call seam_real(buffer_key(i), sval, buffer_line(i), p%e_scale)
         end select
      end do
      allocate (p%tbl(n_carry))
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code) then
            k = k + 1
            call spectrum_frag_w(i, p%tbl(k)%w, p%tbl(k)%c)
         end if
      end do
      call wigner_init(p)
   end subroutine reg_samp_wigner_init
end module reg_samp_wigner
