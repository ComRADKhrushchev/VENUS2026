module test_potentials
!*******************************************************************************
! TEST PES suite — system-independent potentials for VENUS core validation
!
! Potentials, selected at runtime via input keyword
!   TEST_PES = HARMONIC | MORSE | LEPS | EMT-NN     (default HARMONIC)
!
! Unit convention (the PES interface contract, see interface_TEST.f90):
! positions in angstrom (Q), potential in eV, gradient in eV/angstrom.
! interface_TEST converts eV -> kcal/mol (x23.0605) in POT0 and accumulates
! -grad*23.0605*C1 into PDOT.
!
! Analytic references (verified in cases/):
!   HARMONIC: omega_i = sqrt(k_i/m)              — cases/integrator_matrix
!   MORSE:    E(n) = we(n+1/2) - wexe(n+1/2)^2    — cases/ebk_fixed_nj, cases/mb_thermal
!   LEPS:     collinear H3 barrier ~9.8 kcal/mol  — cases/barrier_saddle
!
! EMT-NN: C/Au(111) EMT+neural-network surface PES (engine transplanted
! verbatim into src_TEST/emt_nn_pes.f90, module emt_nn_pes_venus).
!   Atom order: atom 1 = C adsorbate, atoms 2..NATOMS = Au slab (order
!   must match the slab file). Slab-file/weight-file paths default to
!   the CWD and are keyword-overridable:
!     EMTNN_SLAB_FILE    (default Slab.xyz)
!     EMTNN_WEIGHTS_FILE (default nn_weights_emt_nn.txt)
!     EMTNN_ALAT         (default 2.88 A)
!     EMTNN_NSIDE        (default 6)
!   Reference assets: data/emt_nn/.
!*******************************************************************************
   use input_parser
   use venus_params, only: C4
   use venus_data
   use emt_nn_pes_venus, only: emt_nn_init, calc_emt_nn_energy, emt_nn_box, &
                               emt_nn_natoms, emt_nn_rcut, emt_nn_acut
   implicit none
   private
   public :: init_test_potentials, test_pot_v, test_pot_vg, test_pes_name

   integer, parameter :: PES_HARMONIC = 1, PES_MORSE = 2, PES_LEPS = 3, &
                         PES_EMTNN = 4
   integer, save :: test_pes_id = PES_HARMONIC

   ! HARMONIC: independent anisotropic well on every atom
   !   V = sum_i 0.5*[kx*(x_i-x0)^2 + ky*(y_i-y0)^2 + kz*(z_i-z0)^2]
   real(8), save :: h_k(3) = (/ 1.0d0, 1.0d0, 1.0d0 /)   ! eV/angstrom^2
   real(8), save :: h_x0(3) = 0.0d0                       ! angstrom

   ! MORSE: diatomic pair potential on atoms 1-2
   !   V = De*(em-1)^2 - De, em = exp(-a(r-re))  => V(re) = -De, V(inf) = 0
   !   zero of energy at dissociation; well depth De (standard spectroscopy
   !   convention: E_n measured from well bottom = we(n+1/2) - wexe(n+1/2)^2)
   real(8), save :: m_de = 4.746d0      ! eV    (H2)
   real(8), save :: m_re = 1.401d0      ! angstrom  (H2)
   real(8), save :: m_a  = 1.028d0      ! 1/angstrom (H2)

   ! LEPS: triatomic on atoms 1-2-3, classic H3 parameters (CGM formulation)
   !   pairs p=1:(1,2) 2:(2,3) 3:(1,3)
   !   Q_p(r) = De/2*[(3+d)e^-2a(r-re) - (2+2d)e^-a(r-re)]        (triplet)
   !   A_p(r) = De/2*[(3-d)e^-2a(r-re) + (2-2d)e^-a(r-re)] + De*d  (singlet)
   !   V = (Q1+Q2+Q3)/2 - sqrt((A1-Q1)^2 + (A2-Q2)^2 + (A3-Q3)^2)/2
   real(8), save :: l_de    = 4.746d0   ! eV
   real(8), save :: l_re    = 1.401d0   ! angstrom
   real(8), save :: l_a     = 1.028d0   ! 1/angstrom
   real(8), save :: l_delta = 0.164d0   ! Sato parameter

   character(len=16), save :: test_pes_name = 'HARMONIC'

contains

   subroutine init_test_potentials()
      character(len=64) :: s
      character(len=256) :: slab_file, weights_file
      character(len=2) :: elem
      real(8) :: a_lat, blx, bly, bskew, cx, cy, cz
      integer :: n_side, ios, u_slab, n_hdr, i_sl
      logical :: ex
      test_pes_name = upper(get_str('TEST_PES', 'HARMONIC'))
      select case (trim(test_pes_name))
      case ('HARMONIC'); test_pes_id = PES_HARMONIC
      case ('MORSE');    test_pes_id = PES_MORSE
      case ('LEPS');     test_pes_id = PES_LEPS
      case ('EMT-NN');   test_pes_id = PES_EMTNN
      case default
         write(6,*) 'ERROR: TEST_PES must be HARMONIC, MORSE, LEPS, or EMT-NN (got ', &
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
         write(6,'(A,3F12.6)') ' k (eV/A^2) :', h_k
         write(6,'(A,3F12.6)') ' x0 (A)     :', h_x0
      case (PES_MORSE)
         write(6,'(A,F12.6)') ' De (eV)   :', m_de
         write(6,'(A,F12.6)') ' re (A)    :', m_re
         write(6,'(A,F12.6)') ' a (1/A)   :', m_a
      case (PES_LEPS)
         write(6,'(A,F12.6)') ' De (eV)    :', l_de
         write(6,'(A,F12.6)') ' re (A)     :', l_re
         write(6,'(A,F12.6)') ' a (1/A)    :', l_a
         write(6,'(A,F12.6)') ' Sato delta :', l_delta
      case (PES_EMTNN)
         ! Engine assets: paths default to the CWD (loud-stop on missing),
         ! keyword-overridable so case dirs can point at data/emt_nn/.
         slab_file    = trim(get_str('EMTNN_SLAB_FILE',    'Slab.xyz'))
         weights_file = trim(get_str('EMTNN_WEIGHTS_FILE', 'nn_weights_emt_nn.txt'))
         a_lat        = get_real('EMTNN_ALAT', 2.88d0)
         n_side       = get_int ('EMTNN_NSIDE', 6)
         write(6,'(A,A)')    ' slab file    :', trim(slab_file)
         write(6,'(A,A)')    ' weight file  :', trim(weights_file)
         write(6,'(A,F12.6)') ' a_lat (A)   :', a_lat
         write(6,'(A,I12)')   ' n_side      :', n_side
         ! Loads slab baseline + weights; any missing asset/section stops inside.
         call emt_nn_init(slab_file, weights_file, a_lat, n_side)
         write(6,'(A,I0,A)') ' slab atoms  :', emt_nn_natoms(), ' (NATOMS must be 1 + this)'
         write(6,'(A,F12.6)') ' rcut (A)    :', emt_nn_rcut()
         write(6,'(A,F12.6)') ' acut (1/A)  :', emt_nn_acut()

         ! Engine cell overrides the input cell (2D interface_EMT_NN POTPRE
         ! convention): the PES is only self-consistent with the slab cell it
         ! was fitted on, and SURF.f impact sampling / GWRITE PBC use these.
         call emt_nn_box(blx, bly, bskew)
         BOXLX = blx
         BOXLY = bly
         SKEW  = bskew
         write(6,'(A,F8.4,A,F8.4,A,F6.1,A)') &
            ' engine cell : ', BOXLX, ' x ', BOXLY, ' A, skew ', SKEW/C4, &
            ' deg (input BOXLX/BOXLY/SKEW overridden by engine geometry)'

         ! Atom-count cross-check (loud-stop, 2D POTPRE convention).
         ! POTPRE runs before venus_input reads the NATOMB keyword, so fill
         ! it here the same way the 2D interface does.
         NATOMB(1) = NATOMS - NATOMA(1)
         if (NATOMA(1) /= 1) then
            write(6,*) 'ERROR: TEST_PES=EMT-NN requires NATOMA = 1 (got ', &
                       NATOMA(1), '): atom 1 = adsorbate, 2..NATOMS = slab'
            stop
         end if
         if (NATOMB(1) /= emt_nn_natoms()) then
            write(6,*) 'ERROR: NATOMB = ', NATOMB(1), &
                       ' but slab file has ', emt_nn_natoms(), ' atoms'
            stop
         end if

         ! Slab geometry into QZB (2D interface_EMT_NN POTPRE pattern):
         ! SELECT.setup_b_coords copies QZB into the global Q array, and
         ! calc_emt_nn_energy syncs the Q slab columns into the engine on
         ! every call. Without this the engine receives a zero slab (all Au
         ! at the origin) and the EMT energy diverges (~-6.6 keV).
         inquire(file=slab_file, exist=ex)
         if (.not. ex) then
            write(6,*) 'ERROR: slab file not found on re-read: ', trim(slab_file)
            stop
         end if
         open(newunit=u_slab, file=slab_file, status='old', iostat=ios)
         if (ios /= 0) then
            write(6,*) 'ERROR: cannot open slab file: ', trim(slab_file)
            stop
         end if
         read(u_slab, *, iostat=ios) n_hdr
         if (ios /= 0 .or. n_hdr /= NATOMB(1)) then
            write(6,*) 'ERROR: slab file header declares ', n_hdr, &
                       ' atoms, expected ', NATOMB(1)
            close(u_slab)
            stop
         end if
         read(u_slab, *, iostat=ios)   ! comment line
         if (ios /= 0) then
            write(6,*) 'ERROR: cannot read slab file comment line'
            close(u_slab)
            stop
         end if
         do i_sl = 1, NATOMB(1)
            read(u_slab, *, iostat=ios) elem, cx, cy, cz
            if (ios /= 0) then
               write(6,*) 'ERROR: short read in slab file at atom ', i_sl
               close(u_slab)
               stop
            end if
            QZB(1, 3*i_sl-2) = cx
            QZB(1, 3*i_sl-1) = cy
            QZB(1, 3*i_sl)   = cz
         end do
         close(u_slab)
         write(6,'(A,I0,A)') ' Slab loaded into QZB: ', NATOMB(1), ' atoms'
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

   !--- energy (eV) + gradient (eV/angstrom); g = dV/dq -------------------------
   subroutine test_pot_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      g = 0.0d0
      select case (test_pes_id)
      case (PES_HARMONIC); call harm_vg(natom, q, v, g)
      case (PES_MORSE);    call morse_vg(natom, q, v, g)
      case (PES_LEPS);     call leps_vg(natom, q, v, g)
      case (PES_EMTNN);    call emtnn_vg(natom, q, v, g)
      end select
   end subroutine test_pot_vg

   !--- EMT-NN: atom 1 = C adsorbate, atoms 2..n = Au slab -----------------------
   ! Engine returns E[eV] and forces(3,n) [eV/A] with f = -dE/dR, so the
   ! gradient contract (g = +dV/dq) is g = -forces.
   subroutine emtnn_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      real(8) :: q_emt(3, natom)
      real(8) :: forces(3, natom)
      integer :: k
      do k = 1, natom
         q_emt(1, k) = q(3*(k-1)+1)
         q_emt(2, k) = q(3*(k-1)+2)
         q_emt(3, k) = q(3*(k-1)+3)
      end do
      call calc_emt_nn_energy(natom, 1, q_emt, v, forces)
      do k = 1, natom
         g(3*(k-1)+1) = -forces(1, k)
         g(3*(k-1)+2) = -forces(2, k)
         g(3*(k-1)+3) = -forces(3, k)
      end do
   end subroutine emtnn_vg

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
