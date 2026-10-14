!=====================================================================
! check_methods_weave.f90 - criterion: the methods fixed-name weave family
!   (w) legs: vocabulary single-sourcing - input's dist_tbl is the only
!       word<->code table; sampler's samp_code delegates to it
!   (added by later tasks: drill_methods_weave.sh owns the tree legs)
! Build (gfortran, from src/):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 &
!            globals/config.f90 utils/specio.f90 interface/sysdef_interface.f90 input.f90 &
!            globals/config_atoms.f90 globals/state.f90 utils/rng.f90 utils/geometry.f90 &
!            utils/linalg.f90 utils/hessian.f90 utils/geom_init.f90 utils/spectrum.f90 &
!            interface/elec_interface.f90 interface/force_interface.f90 dynamics/propagator.f90 &
!            interface/bath.f90 interface/sampler.f90 &
!            ../scripts/check_methods_weave.f90 -o build/check_methods_weave.exe
!=====================================================================
program check_methods_weave
   use input,   only: input_dist_code, input_dist_word
   use sampler, only: samp_code
   implicit none
   integer :: n_pass, n_fail
   n_pass = 0; n_fail = 0

   call eq_i('w1 dist_code uppercase', input_dist_code('BOLTZMANN'), 1)
   call eq_i('w2 dist_code lowercase', input_dist_code('boltzmann'), 1)
   call eq_i('w3 dist_code absent',    input_dist_code('nosuchword'), 0)
   call eq_s('w4 dist_word by code',   input_dist_word(10), 'wigner')
   call eq_s('w5 dist_word invalid',   input_dist_word(99), '')
   call eq_i('w6 samp_code delegates', samp_code('ebk'), 2)
   call eq_i('w7 samp_code absent',    samp_code('nosuchword'), 0)
   call eq_s('w8 dist_word lowercase form', input_dist_word(5), 'surface_oscillator')

   write (6, '(a,i0,a,i0)') 'check_methods_weave: pass=', n_pass, ' fail=', n_fail
   if (n_fail > 0) stop 1
contains
   subroutine eq_i(tag, got, want)
      character(len=*), intent(in) :: tag
      integer, intent(in) :: got, want
      if (got == want) then
         n_pass = n_pass + 1; write (6, '(a)') 'PASS: '//trim(tag)
      else
         n_fail = n_fail + 1; write (6, '(a,i0,a,i0)') 'FAIL: '//trim(tag)//' got=', got, ' want=', want
      end if
   end subroutine eq_i

   subroutine eq_s(tag, got, want)
      character(len=*), intent(in) :: tag, got, want
      if (trim(got) == trim(want)) then
         n_pass = n_pass + 1; write (6, '(a)') 'PASS: '//trim(tag)
      else
         n_fail = n_fail + 1; write (6, '(a)') 'FAIL: '//trim(tag)//' got=['//trim(got)//'] want=['//trim(want)//']'
      end if
   end subroutine eq_s
end program check_methods_weave
