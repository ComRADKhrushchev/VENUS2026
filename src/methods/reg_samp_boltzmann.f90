!=====================================================================
! reg_samp_boltzmann.f90 - register the 'boltzmann' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: buffered input keys -> params receiver
!   (absent keys keep defaults) -> one frequency-table pull per carrier
!   (spectrum buffer) -> the member's own init; an unselected member skips.
!   The tables are assembly data, not input keys - a selected member aborts
!   with a named error at its own unallocated check.
!=====================================================================
module reg_samp_boltzmann
   use sampler,        only: samp_reg, samp_code
   use config,         only: reactants
   use config_atoms,   only: list_atoms
   use input,          only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,      only: seam_real
   use spectrum_interface, only: spectrum_frag_w
   use samp_boltzmann, only: boltz_draw, boltz_realize, boltz_params_t, boltz_init
   implicit none
   private
   public :: reg_samp_boltzmann_all
   public :: reg_samp_boltzmann_init
contains
   !------------------------------------------------------------------
   ! reg_samp_boltzmann_all() - fill the member-table row (one call, one row)
   subroutine reg_samp_boltzmann_all()
      call samp_reg('boltzmann', boltz_draw, boltz_realize)
   end subroutine reg_samp_boltzmann_all

   !------------------------------------------------------------------
   ! reg_samp_boltzmann_init() - assembly seam: buffered keys -> the
   !                              frequency-table pull -> member init.
   !                              The tables are assembly data: pulled from
   !                              the SPECTRUM_SOURCE-selected producer
   !                              (the internal derivation on COMPUTE,
   !                              the spectrum-table file on MANUAL) for
   !                              EVERY carrier fragment, one ragged
   !                              table per carrier
   !------------------------------------------------------------------
   subroutine reg_samp_boltzmann_init()
      type(boltz_params_t) :: p          ! member-parameter receiver (field defaults)
      character(len=512) :: sval         ! one buffered raw value
      integer :: i, k, code, n_carry, n_word
      if (.not. allocated(reactants%dist_scheme)) return
      code = samp_code('boltzmann')
      n_word = count(reactants%dist_scheme == code)
      if (n_word == 0) return
      ! build tables only for polyatomic carriers; a monoatomic carrier rides
      ! the incident channel alone (BOLTZMANN on it = "no internal sampling")
      n_carry = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code .and. list_atoms%frag(i)%nat > 1) then
            n_carry = n_carry + 1
         end if
      end do
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('T_VIB_A'); call seam_real(buffer_key(i), sval, buffer_line(i), p%t_vib_a)
         case ('T_VIB_B'); call seam_real(buffer_key(i), sval, buffer_line(i), p%t_vib_b)
         end select
      end do
      allocate (p%tbl(n_carry))
      k = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) == code .and. list_atoms%frag(i)%nat > 1) then
            k = k + 1
            call spectrum_frag_w(i, p%tbl(k)%w, p%tbl(k)%c)
         end if
      end do
      call boltz_init(p)
   end subroutine reg_samp_boltzmann_init
end module reg_samp_boltzmann
