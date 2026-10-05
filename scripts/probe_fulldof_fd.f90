!=====================================================================
! probe_fulldof_fd.f90 - supplementary full-dof FD probe of the armed
!   composed fold (E_elec + E_mf + E0 vs g): check_cau_peri legs
!   (h)/(h2) probe only dofs 1-3 (the S atom; dof(4)=1 remap); this
!   probe sweeps ALL NINE dofs including the incident-C side (7-9) and
!   the ghost-Y side (4-6), with a nondiagonal rho (100 prop calls at
!   the probe geometry). rho frozen: only cau_peri_pes runs.
!=====================================================================
program probe_fulldof_fd
   use cau_peri_pot,  only: cau_peri_load
   use cau_peri_tdhf, only: cau_peri_tdhf_setup, cau_peri_pes, cau_peri_prop
   use state,         only: state_t, state_create, state_destroy
   implicit none
   character(len=*), parameter :: shelf = '../systems/inactive/cau_peri/'
   real(8), parameter :: ev_kcal = 23.0605d0
   real(8) :: q0(9), mass(3), g(9), v, vp, vm, dx, rel, tol
   integer :: k, npass, nfail
   type(state_t) :: sta
   npass = 0
   nfail = 0
   dx = 1.0d-4
   call copy_file(trim(shelf)//'Analytic_Potential.txt', 'Analytic_Potential.txt')
   call cau_peri_load()
   mass = [196.97d0, 196.97d0, 12.01d0]
   q0 = 0.0d0
   q0(4:6) = [0.0d0, 0.0d0, 0.9d0]
   q0(7:9) = [0.37d0, 0.21d0, 2.2d0]
   ! arm + evolve rho to a nondiagonal state at the probe geometry
   call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', 4, 0, 0.05d0)
   call state_create(sta, 3, 1)
   sta%q = q0
   do k = 1, 100
      sta%t = k*0.05d0
      call cau_peri_prop(sta, 0.05d0)
   end do
   call state_destroy(sta)
   ! FD sweep over all nine dofs (rho frozen in module state)
   do k = 1, 9
      call cau_peri_pes(q0, mass, g, v)
      call fd_v(q0, k, +dx, vp)
      call fd_v(q0, k, -dx, vm)
      rel = abs(g(k) - (vp - vm)/(2*dx))/max(1.0d0, abs(g(k)))
      tol = 1.0d-5
      if (rel <= tol) then
         npass = npass + 1
         print '(a,i0,a,es9.2)', 'dof ', k, ' PASS  rel=', rel
      else
         nfail = nfail + 1
         print '(a,i0,a,es9.2,a,es12.5,a,es12.5)', 'dof ', k, ' FAIL  rel=', rel, &
            '  an=', g(k), '  fd=', (vp - vm)/(2*dx)
      end if
   end do
   print '(a,i0,a,i0)', 'probe_fulldof_fd: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) call exit(1)

contains

   subroutine fd_v(qb, k, h, vv)
      real(8), intent(in) :: qb(:), h
      integer, intent(in) :: k
      real(8), intent(out) :: vv
      real(8) :: qp(9), gg(9)
      qp = qb
      qp(k) = qb(k) + h
      call cau_peri_pes(qp, mass, gg, vv)
   end subroutine fd_v

   subroutine copy_file(from, to)
      character(len=*), intent(in) :: from, to
      integer :: ua, ub, ios
      character(len=256) :: line
      open (newunit=ua, file=from, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print '(a)', 'unreadable '//trim(from)
         call exit(1)
      end if
      open (newunit=ub, file=to, status='replace', action='write')
      do
         read (ua, '(a)', iostat=ios) line
         if (ios /= 0) exit
         write (ub, '(a)') trim(line)
      end do
      close (ua); close (ub)
   end subroutine copy_file
end program probe_fulldof_fd
