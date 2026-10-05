! dbg_h_compare.f90 - h tilde(z) at the top site through three faces:
!   the NN branch (cau_peri + nn_weights.txt), the rst channel table
!   (cau_tab + tab_rst_sparse.txt), and the analytic branch. Temporary
!   diagnostic for the 2026-10-01 rst-vs-NN discrepancy hunt.
program dbg_h_compare
   use cau_tab_pot,  only: cau_tab_load, cau_tab_hubbard_val
   use cau_peri_pot, only: cau_peri_load, cau_peri_hubbard_val
   implicit none
   real(8) :: qs(9), h(2), u, j, vcp, hn(2), un, jn, vcn
   real(8) :: z
   integer :: i
   call copy_file('../systems/inactive/cau_peri/nn_weights.txt', 'nn_weights.txt')
   call copy_file('../systems/inactive/cau_peri/Analytic_Potential.txt', 'Analytic_Potential.txt')
   call cau_peri_load('nn_weights.txt')      ! the NN branch
   call cau_tab_load('C:/Users/13614/Desktop/MAIN/01_WORK/data/c_au/tables/tab_rst_sparse.txt')
   print '(a)', '   z      h_NN     h_tab   |dh/dz|_tab(eV/A)   U_NN     U_tab    J_NN     J_tab'
   do i = 0, 45
      z = 1.5d0 + i*0.1d0
      qs = 0.0d0
      qs(6) = 0.9d0
      qs(7:9) = [0.0d0, 0.0d0, z]            ! directly above the anchor = top
      call cau_peri_hubbard_val(qs, hn, un, jn, vcn)
      call cau_tab_hubbard_val(qs, h, u, j, vcp)
      print '(f6.2, 2f9.3, a, 2f9.3, 2f9.3)', z, 0.5d0*(hn(1) + hn(2)), 0.5d0*(h(1) + h(2)), &
         '              ', un, u, jn, j
   end do

contains
   subroutine copy_file(from, to)
      character(len=*), intent(in) :: from, to
      integer :: ua, ub, ios
      character(len=256) :: line
      open (newunit=ua, file=from, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print '(a)', 'fatal: unreadable '//trim(from)
         stop 1
      end if
      open (newunit=ub, file=to, status='replace', action='write')
      do
         read (ua, '(a)', iostat=ios) line
         if (ios /= 0) exit
         write (ub, '(a)') trim(line)
      end do
      close (ua); close (ub)
   end subroutine copy_file
end program dbg_h_compare
