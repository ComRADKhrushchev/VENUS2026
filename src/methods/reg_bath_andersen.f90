!=====================================================================
! reg_bath_andersen.f90 - register the 'andersen' bath species
! Design:
!   One bath_reg call fills the species-table row under the word the BATH
!   selector accepts; the reg file is the atomic registration unit (a new
!   species = its method file + one thin reg file; zero edits elsewhere).
!   _init is the assembly seam: an armed BATH provision pulls the species
!   keys (T_BATH/NU_COLL) into the receiver and inits the species, then
!   provisions the loop shape via bath_setup (BATH = NONE skips both - the
!   species stays unloaded and the loop stays bitwise inert).
!=====================================================================
module reg_bath_andersen
   use bath,            only: bath_reg, bath_setup
   use config,          only: reactants
   use input,           only: buffer_n_rows, buffer_key, buffer_val, buffer_line
   use seam_cast,       only: seam_real
   use bath_andersen, only: andersen_sweep, andersen_params_t, andersen_init
   implicit none
   private
   public :: reg_bath_andersen_all
   public :: reg_bath_andersen_init
contains
   !------------------------------------------------------------------
   ! reg_bath_andersen_all() - fill the species-table row (one call, one row)
   subroutine reg_bath_andersen_all()
      call bath_reg('andersen', andersen_sweep)
   end subroutine reg_bath_andersen_all

   !------------------------------------------------------------------
   ! reg_bath_andersen_init() - assembly seam: BATH provision -> species
   !                              keys -> species init -> loop provisioning
   subroutine reg_bath_andersen_init()
      type(andersen_params_t) :: p   ! species-parameter receiver (field defaults)
      character(len=512) :: sval       ! one buffered raw value
      integer :: i
      if (reactants%bath /= 2) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('T_BATH');  call seam_real(buffer_key(i), sval, buffer_line(i), p%t_bath)
         case ('NU_COLL'); call seam_real(buffer_key(i), sval, buffer_line(i), p%nu_coll)
         end select
      end do
      call andersen_init(p)
      call bath_setup()
   end subroutine reg_bath_andersen_init
end module reg_bath_andersen
