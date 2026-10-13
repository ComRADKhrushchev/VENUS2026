!=====================================================================
! reg_samp_boltzmann.f90 - register the 'boltzmann' distribution member
! Design:
!   One samp_reg call fills the member-table row under the scheme word
!   the input grammar accepts for DIST_SCHEME; the reg file is the
!   atomic registration unit (a new member = its method file + one thin
!   reg file; zero edits to existing files).
!   _init is the assembly seam: the carrier's bracket parameters
!   (T_VIB - one row per word-carrying fragment, the 2026-10-13 bracket
!   ruling; a bare word keeps the field default and the member's own
!   validation speaks) -> one frequency-table pull per polyatomic carrier
!   (spectrum buffer) -> the member's own init; an unselected member
!   skips. The tables are assembly data, not input keys - a selected
!   member aborts with a named error at its own unallocated check.
!=====================================================================
module reg_samp_boltzmann
   use sampler,        only: samp_reg, samp_code
   use config,         only: reactants
   use config_atoms,   only: list_atoms
   use input,          only: par_real
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
   ! reg_samp_boltzmann_init() - assembly seam: per-carrier bracket keys ->
   !                              the frequency-table pull -> member init.
   !                              Parameters index EVERY word-carrying
   !                              fragment (a monoatomic carrier's bracket is
   !                              consumed for the audit but carries no
   !                              table); the tables are assembly data:
   !                              pulled from the SPECTRUM_SOURCE-selected
   !                              producer for every POLYATOMIC carrier, one
   !                              ragged table per carrier
   !------------------------------------------------------------------
   subroutine reg_samp_boltzmann_init()
      type(boltz_params_t) :: p          ! member-parameter receiver (field defaults)
      integer :: i, k, wc, code, n_carry, n_word
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
      allocate (p%t_vib(n_word))
      p%t_vib = 0.0d0                     ! absent bracket name keeps the default;
                                           ! boltz_init's own validation rejects it
      wc = 0
      do i = 1, size(reactants%dist_scheme)
         if (reactants%dist_scheme(i) /= code) cycle
         wc = wc + 1
         call par_real('BOLTZMANN', i, 'T_VIB', p%t_vib(wc))
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
