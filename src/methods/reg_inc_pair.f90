!=====================================================================
! reg_inc_pair.f90 - register the 'incident_pair' member (guarded: the
!   molecular paradigm with two or more fragments)
! Design:
!   One samp_reg call fills the row when the run derived the molecular
!   paradigm and carries a collision pair (nf > 1; the member's own init
!   holds the exact nf == 2 statute) - the paradigm-guarded registration
!   keeps exactly ONE incident member armed per process. _init is the
!   assembly seam: buffered beam keys -> params receiver -> the member's
!   own init; the surface-only AIM_X/AIM_Y keys are rejected loudly when
!   they appear in a pair run (a drop-site key on a pair collision is a
!   configuration error, not a silently ignored row).
!=====================================================================
module reg_inc_pair
   use sampler,       only: samp_reg
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use input,         only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,     only: seam_real, seam_int
   use beam_laws,     only: beam_params_t
   use inc_pair,      only: inc_pair_draw, inc_pair_realize, inc_pair_init
   implicit none
   private
   public :: reg_inc_pair_all
   public :: reg_inc_pair_init
contains
   !------------------------------------------------------------------
   ! reg_inc_pair_all() - fill the member-table row (molecular pair only)
   !------------------------------------------------------------------
   subroutine reg_inc_pair_all()
      if (.not. armed_here()) return
      call samp_reg('incident_pair', inc_pair_draw, inc_pair_realize)
   end subroutine reg_inc_pair_all

   !------------------------------------------------------------------
   ! reg_inc_pair_init() - assembly seam: buffered keys -> receiver ->
   !                       member init (molecular pair only)
   !------------------------------------------------------------------
   subroutine reg_inc_pair_init()
      type(beam_params_t) :: p        ! parameter receiver (field defaults)
      character(len=512) :: sval      ! one buffered raw value
      integer :: i
      if (.not. armed_here()) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('N_E_REL');   call seam_int(buffer_key(i), sval, buffer_line(i), p%n_e_rel)
         case ('T_TRANS');   call seam_real(buffer_key(i), sval, buffer_line(i), p%t_trans)
         case ('V_WIDTH');   call seam_real(buffer_key(i), sval, buffer_line(i), p%v_width)
         case ('N_THTA');    call seam_int(buffer_key(i), sval, buffer_line(i), p%n_thta)
         case ('THTA_MAX');  call seam_real(buffer_key(i), sval, buffer_line(i), p%thta_max)
         case ('N_CHI');     call seam_int(buffer_key(i), sval, buffer_line(i), p%n_chi)
         case ('CHI');       call seam_real(buffer_key(i), sval, buffer_line(i), p%chi)
         case ('N_B');       call seam_int(buffer_key(i), sval, buffer_line(i), p%n_b)
         case ('N_AIM');     call seam_int(buffer_key(i), sval, buffer_line(i), p%n_aim)
         case ('AIM_X', 'AIM_Y')
            call stop_regp('reg_inc_pair_init', 'key '//trim(buffer_key(i))//' is a '// &
               'surface-paradigm drop-site key (N_AIM = 2 + AIM_X/AIM_Y) - the pair '// &
               'member drops at the origin or per N_B')
         end select
      end do
      call inc_pair_init(p)
   end subroutine reg_inc_pair_init

   ! armed_here() - the paradigm guard (molecular pair; the atom list must
   !                be assembled - the mirror of the dispatch rule)
   logical function armed_here()
      armed_here = allocated(list_atoms%frag) .and. reactants%surface_model == 0 &
                   .and. size(list_atoms%frag) > 1
   end function armed_here

   ! stop_regp(where, msg) - the named-abort channel (STOP 1)
   subroutine stop_regp(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      write (0, '(a)') trim(where)//': fatal (the initial state is not sampled)'
      stop 1
   end subroutine stop_regp
end module reg_inc_pair
