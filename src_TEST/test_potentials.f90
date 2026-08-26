module test_potentials
!*******************************************************************************
! TEST PES suite — system-independent potentials for VENUS core validation
!
! Three analytic potentials, selected at runtime via input keyword
!   TEST_PES = HARMONIC | MORSE | LEPS     (default HARMONIC)
!
! Unit convention (the PES interface contract, see interface_TEST.f90):
! positions in bohr (Q), potential in eV, gradient in eV/bohr.
! interface_TEST converts eV -> kcal/mol (x23.0605) in POT0 and accumulates
! -grad*23.0605*C1 into PDOT.
!
! Analytic references (verified in cases/):
!   HARMONIC: omega_i = sqrt(k_i/m)              — cases/integrator_matrix
!   MORSE:    E(n) = we(n+1/2) - wexe(n+1/2)^2    — cases/ebk_fixed_nj, cases/mb_thermal
!   LEPS:     collinear H3 barrier ~9.8 kcal/mol  — cases/barrier_saddle
!*******************************************************************************
   use input_parser
   implicit none
   private
   public :: init_test_potentials, test_pot_v, test_pot_vg, test_pes_name

   integer, parameter :: PES_HARMONIC = 1, PES_MORSE = 2, PES_LEPS = 3
   integer, save :: test_pes_id = PES_HARMONIC

   ! HARMONIC: independent anisotropic well on every atom
   !   V = sum_i 0.5*[kx*(x_i-x0)^2 + ky*(y_i-y0)^2 + kz*(z_i-z0)^2]
   real(8), save :: h_k(3) = (/ 1.0d0, 1.0d0, 1.0d0 /)   ! eV/bohr^2
   real(8), save :: h_x0(3) = 0.0d0                       ! bohr

   ! MORSE: diatomic pair potential on atoms 1-2
   !   V = De*(em-1)^2 - De, em = exp(-a(r-re))  => V(re) = -De, V(inf) = 0
   !   zero of energy at dissociation; well depth De (standard spectroscopy
   !   convention: E_n measured from well bottom = we(n+1/2) - wexe(n+1/2)^2)
   real(8), save :: m_de = 4.746d0      ! eV    (H2)
   real(8), save :: m_re = 1.401d0      ! bohr  (H2)
   real(8), save :: m_a  = 1.028d0      ! 1/bohr(H2)

   ! LEPS: triatomic on atoms 1-2-3, classic H3 parameters (CGM formulation)
   !   pairs p=1:(1,2) 2:(2,3) 3:(1,3)
   !   Q_p(r) = De/2*[(3+d)e^-2a(r-re) - (2+2d)e^-a(r-re)]        (triplet)
   !   A_p(r) = De/2*[(3-d)e^-2a(r-re) + (2-2d)e^-a(r-re)] + De*d  (singlet)
   !   V = (Q1+Q2+Q3)/2 - sqrt((A1-Q1)^2 + (A2-Q2)^2 + (A3-Q3)^2)/2
   real(8), save :: l_de    = 4.746d0   ! eV
   real(8), save :: l_re    = 1.401d0   ! bohr
   real(8), save :: l_a     = 1.028d0   ! 1/bohr
   real(8), save :: l_delta = 0.164d0   ! Sato parameter

   character(len=16), save :: test_pes_name = 'HARMONIC'

contains

   subroutine init_test_potentials()
      character(len=64) :: s
      test_pes_name = upper(get_str('TEST_PES', 'HARMONIC'))
      select case (trim(test_pes_name))
      case ('HARMONIC'); test_pes_id = PES_HARMONIC
      case ('MORSE');    test_pes_id = PES_MORSE
      case ('LEPS');     test_pes_id = PES_LEPS
      case default
         write(6,*) 'ERROR: TEST_PES must be HARMONIC, MORSE, or LEPS (got ', &
                    trim(test_pes_name), ')'
         stop
      end select

      if (has_keyword('HARM_K')) call get_real_arr('HARM_K', h_k, 3)
      if (has_keyword('HARM_X0')) call get_real_arr('HARM_X0', h_x0, 3)

      m_de = get_real('MORSE_DE', m_de)
      m_re = get_real('MORSE_RE', m_re)
      m_a  = get_real('MORSE_A',  m_a)

      l_de    = get_real('LEPS_DE',    l_de)
      l_re    = get_real('LEPS_RE',    l_re)
      l_a     = get_real('LEPS_A',     l_a)
      l_delta = get_real('LEPS_DELTA', l_delta)

      write(6,*) '=== TEST PES : ', trim(test_pes_name), ' ==='
      select case (test_pes_id)
      case (PES_HARMONIC)
         write(6,'(A,3F12.6)') ' k (eV/bohr^2) :', h_k
         write(6,'(A,3F12.6)') ' x0 (bohr)     :', h_x0
      case (PES_MORSE)
         write(6,'(A,F12.6)') ' De (eV)   :', m_de
         write(6,'(A,F12.6)') ' re (bohr) :', m_re
         write(6,'(A,F12.6)') ' a (1/bohr):', m_a
      case (PES_LEPS)
         write(6,'(A,F12.6)') ' De (eV)    :', l_de
         write(6,'(A,F12.6)') ' re (bohr)  :', l_re
         write(6,'(A,F12.6)') ' a (1/bohr) :', l_a
         write(6,'(A,F12.6)') ' Sato delta :', l_delta
      end select
   end subroutine init_test_potentials

   !--- energy only (eV) --------------------------------------------------------
   subroutine test_pot_v(natom, q, v)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v
      real(8) :: g(3*natom)
      call test_pot_vg(natom, q, v, g)
   end subroutine test_pot_v

   !--- energy (eV) + gradient (eV/bohr); g = dV/dq -----------------------------
   subroutine test_pot_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      g = 0.0d0
      select case (test_pes_id)
      case (PES_HARMONIC); call harm_vg(natom, q, v, g)
      case (PES_MORSE);    call morse_vg(natom, q, v, g)
      case (PES_LEPS);     call leps_vg(natom, q, v, g)
      end select
   end subroutine test_pot_vg

   !--- anisotropic harmonic well on each atom -----------------------------------
   subroutine harm_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      integer :: i, k
      v = 0.0d0
      do i = 1, natom
         do k = 1, 3
            v = v + 0.5d0*h_k(k)*(q(3*(i-1)+k) - h_x0(k))**2
            g(3*(i-1)+k) = h_k(k)*(q(3*(i-1)+k) - h_x0(k))
         end do
      end do
   end subroutine harm_vg

   !--- Morse pair on atoms 1-2 ---------------------------------------------------
   subroutine morse_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      real(8) :: d(3), r, em, dvdr
      d = q(4:6) - q(1:3)
      r = sqrt(dot_product(d, d))
      em = exp(-m_a*(r - m_re))
      v = m_de*(em - 1.0d0)**2 - m_de          ! V(re) = -De, V(inf) = 0
      ! V = De*(em-1)^2 - De  => dV/dr = 2*De*(em-1)*(-a*em)
      dvdr = -2.0d0*m_de*(em - 1.0d0)*m_a*em
      g(1:3) = -dvdr*d/r                        ! atom 1
      g(4:6) =  dvdr*d/r                        ! atom 2
   end subroutine morse_vg

   !--- LEPS triatomic on atoms 1-2-3 ---------------------------------------------
   subroutine leps_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      real(8) :: r(3), u(3,3)        ! pair distances; u(p,k) = d r_p / d(cart k of atom B-A scheme below)
      real(8) :: qq(3), aa(3), dqq(3), daa(3)
      real(8) :: s2, rt, dvdq(3), dvda(3), dvdr(3)
      integer :: p, k
      integer :: pa(3), pb(3)
      data pa /1, 2, 1/  ! pair atoms (A<B): (1,2),(2,3),(1,3)
      data pb /2, 3, 3/

      do p = 1, 3
         call pair_dist(q(3*pa(p)-2:3*pa(p)), q(3*pb(p)-2:3*pb(p)), r(p))
         call leps_branch(r(p), qq(p), aa(p), dqq(p), daa(p))
      end do

      s2 = sum((aa - qq)**2)
      rt = sqrt(s2)
      v = 0.5d0*sum(qq) - 0.5d0*rt
      ! d V/d Q_p = 1/2 + (A_p-Q_p)/(2 rt);  d V/d A_p = -(A_p-Q_p)/(2 rt)
      do p = 1, 3
         dvdq(p) = 0.5d0 + (aa(p) - qq(p))/(2.0d0*rt)
         dvda(p) = -(aa(p) - qq(p))/(2.0d0*rt)
         dvdr(p) = dvdq(p)*dqq(p) + dvda(p)*daa(p)
      end do

      g = 0.0d0
      do p = 1, 3
         ! r_p = |x_B - x_A|: dr/dx_A = -u, dr/dx_B = +u
         call pair_unit(q(3*pa(p)-2:3*pa(p)), q(3*pb(p)-2:3*pb(p)), u(p,:))
         do k = 1, 3
            g(3*(pa(p)-1)+k) = g(3*(pa(p)-1)+k) - dvdr(p)*u(p,k)
            g(3*(pb(p)-1)+k) = g(3*(pb(p)-1)+k) + dvdr(p)*u(p,k)
         end do
      end do
   end subroutine leps_vg

   subroutine pair_dist(a, b, r)
      real(8), intent(in)  :: a(3), b(3)
      real(8), intent(out) :: r
      real(8) :: d(3)
      d = b - a
      r = sqrt(dot_product(d, d))
   end subroutine pair_dist

   subroutine pair_unit(a, b, u)
      real(8), intent(in)  :: a(3), b(3)
      real(8), intent(out) :: u(3)
      real(8) :: d(3)
      d = b - a
      u = d/sqrt(dot_product(d, d))
   end subroutine pair_unit

   ! Q (triplet) and A (singlet) branch energies and r-derivatives
   subroutine leps_branch(r, qe, ae, dqe, dae)
      real(8), intent(in)  :: r
      real(8), intent(out) :: qe, ae, dqe, dae
      real(8) :: em, dem
      em  = exp(-l_a*(r - l_re))
      dem = -l_a*em
      qe = 0.5d0*l_de*((3.0d0 + l_delta)*em*em - (2.0d0 + 2.0d0*l_delta)*em)
      ae = 0.5d0*l_de*((3.0d0 - l_delta)*em*em + (2.0d0 - 2.0d0*l_delta)*em) &
           + l_delta*l_de
      dqe = 0.5d0*l_de*(2.0d0*(3.0d0 + l_delta)*em*dem &
                        - (2.0d0 + 2.0d0*l_delta)*dem)
      dae = 0.5d0*l_de*(2.0d0*(3.0d0 - l_delta)*em*dem &
                        + (2.0d0 - 2.0d0*l_delta)*dem)
   end subroutine leps_branch

   pure function upper(s) result(u)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: u
      integer :: i, ic
      u = s
      do i = 1, len(s)
         ic = iachar(s(i:i))
         if (ic >= 97 .and. ic <= 122) u(i:i) = achar(ic - 32)
      end do
   end function upper

end module test_potentials
