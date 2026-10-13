!=====================================================================
! check_samp_normalmode.f90 - normalmode distribution-member criteria
!                            program (criteria defined FIRST in this
!                            header; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 &
!            globals/config_atoms.f90 globals/state.f90 &
!            utils/linalg.f90 utils/hessian.f90 utils/geom_init.f90 utils/spectrum.f90 interface/elec_interface.f90 interface/force_interface.f90 dynamics/propagator.f90 interface/bath.f90 &
!            utils/rng.f90 utils/geometry.f90 interface/sampler.f90 methods/samp_normalmode.f90 &
!            ../scripts/check_samp_normalmode.f90 -o build/check_samp_normalmode.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter rule: the member is driven by directly
! constructed parameters and a directly filled list_atoms (a three-atom
! carrier, masses 1/2/3 amu at a non-symmetric geometry). Each
! carrier's mode table is SYNTHETIC: six orthonormal mass-weighted columns spanning the
! complement of the three translations (built here by Gram-Schmidt from
! the canonical basis), with arbitrary positive distinct frequencies
! w = [0.8, 1.1, 1.4, 1.9, 2.3, 3.1] rad/(10 fs) - the member's contract
! is the measure on a GIVEN table (frequency derivation from a Hessian
! is assembly-phase mechanism, not member code). Expected distributions
! are computed in this program from theory, not from the member.
! Mass-weighting convention (the member's, restated as the anchor; the
! member seeds the carrier about the COM origin, the family posture):
!   qz_com = qz - COM(qz);  x = sqrt(m)*(q - qz_com),  y = p/sqrt(m)
!   Q = C^T x,  P = C^T y;  conversely q = qz_com + (C Q)/sqrt(m),
!   p = sqrt(m)*(C P).  Mode energy E_i = 0.5*(P_i^2 + w_i^2*Q_i^2).
! Named-abort verification approach (same scheme as the sibling checks):
! every named-abort criterion runs as a CHILD PROCESS - the check re-spawns
! itself with a mode argument (execute_command_line); the child drives
! the REAL member entry into the stop and must print a naming message on
! stderr and exit nonzero. The parent asserts the exit status and the
! stderr content.
! Criteria (carrier: nat = 3, s = n_mode = 6, e_vib = 0.54205243 eV
! (= 12.5 kcal/mol / 23.0605); mode-1 marginal of the sphere split is
! Beta(1/2, 5/2)):
!   (a) total energy exact: 2000 draws; the energy re-expanded from the
!       WRITTEN state (Q = C^T sqrt(m)*(q - qz_com), P = C^T p/sqrt(m),
!       then sum 0.5*(P_i^2 + w_i^2*Q_i^2)) equals e_tot to 1e-10
!       relative every draw
!   (b) zero total momentum AND subspace containment: |sum p| <=
!       1e-9*sqrt(2*M*e_tot) (M the carrier mass); the residual of
!       sqrt(m)*(q - qz_com) and of p/sqrt(m) after projection onto
!       span(C) has norm <= 1e-10 of the vector norm (only mode-subspace
!       content is written)
!   (c) COM at the origin: the COM of the written q is the origin
!       (COM-shifted seed; the displacement is translation-free) to
!       1e-10 Angstrom every draw
!   (d) mode-energy marginal: E_1/e_tot ~ Beta(1/2, 5/2); equal-CDF
!       10-bin chi-square with closed-form quantiles (x = sin^2(theta),
!       CDF(theta) = (3t/4 + sin2t/2 + sin4t/16)/(3*pi/8) solved by
!       bisection), N = 50000, df = 9, critical 16.92 (alpha = 0.05);
!       mean = 1/6 within 5 sigma (Var = 5/144, sigma of the mean =
!       sqrt(Var/N) = 8.3e-4)
!   (e) phase equipartition: mean kinetic energy = e_tot/2 within
!       5 sigma with sigma = e_tot*sqrt(3/(64*N)) (Var(KE) = 3*e_tot^2/64:
!       KE/e_tot = sum_k u_k^2 sin^2(phi_k) with independent uniform phi,
!       Var(sin^2) = 1/8 and E[sum_k u_k^4] = 3/(s+2) = 3/8 for s = 6 -
!       the SUM of the fourth moments, not the single-mode 1/16)
!   (f) non-positive total energy named abort (child): e_vib = 0 ->
!       stderr names normalmode_init and e_vib
!   (g) zero-mode table named abort (child): n_mode = 0 (both tables
!       empty) -> stderr names normalmode_init and the mode count
!   (h) table-length mismatch named abort (child): tbl(1)%c 6 columns vs
!       tbl(1)%w 7 entries -> stderr names normalmode_init and the mismatch
!   (i) table-row shape named abort (child): tbl(1)%c 8 rows vs 3*nat = 9
!       -> stderr names normalmode_init and the shape
!   (j) non-positive frequency named abort (child): tbl(1)%w(2) = 0 ->
!       stderr names normalmode_init and the frequency
!   (k) non-orthonormal table named abort (child): tbl(1)%c column 2
!       replaced by column 1 -> stderr names normalmode_init and
!       orthonormality
!   (l) translation-violating table named abort (child): tbl(1)%c column
!       1 tilted by a translation vector (still unit, still orthogonal to
!       the other columns) -> stderr names normalmode_init and the
!       translation
!   (m) no-fragment-binding named abort (child): stderr names
!       normalmode_init and the normalmode scheme
!   (n) table-count mismatch named abort (child): TWO fragments select
!       the scheme but the parameter set carries ONE table -> stderr
!       names normalmode_init and the carrier count
!   (o) sample-before-init named abort (child): stderr names
!       normalmode_sample
!   (p) dual-carrier equal treatment (positive): two fragments, two
!       tables (the second a copy of the first), one init and one
!       sampling loop - each carrier's atom block carries E =
!       e_vib*ev_dyn exactly (both carriers' q/p are written)
!=====================================================================
program check_samp_normalmode
   use consts,        only: ev_dyn
   use rng,           only: rng_init
   use config,        only: reactants
   use config_atoms, only: list_atoms, list_atoms_clear
   use state,         only: state_t, state_create, state_destroy
   use samp_normalmode, only: normalmode_params_t, normalmode_init, normalmode_sample
   implicit none
   type(state_t) :: sta
   type(normalmode_params_t) :: par
   integer, parameter :: ndraw_exact = 2000
   integer, parameter :: ndraw_stat = 50000
   integer, parameter :: nbin = 10
   integer, parameter :: s_mode = 6
   integer :: nfail = 0, draw, ib, k, it
   integer :: freq(nbin)
   real(8) :: e_vib, e_tot, amean, emean, sig, chi2
   real(8) :: qz(9), mass3(3), com_new(3), psc
   real(8) :: w_tab(s_mode), c_tab(9, s_mode)
   real(8) :: x_vec(9), y_vec(9), qmode(s_mode), pmode(s_mode)
   real(8) :: e_sum, ke, xres, yres, psum(3), kmean
   real(8) :: x_edge(nbin - 1), lo, hi, mid, gmax, e1_share
   character(len=32) :: mode
   logical :: ok_all

   allocate (par%e_vib(1))   ! single-carrier check face
   par%e_vib(1) = 0.43364194d0
   e_vib = 0.54205243d0   ! [eV] = 12.5 kcal/mol / 23.0605
   e_tot = e_vib*ev_dyn
   mass3 = [ 1.0d0, 2.0d0, 3.0d0 ]
   psc = sqrt(2.0d0*sum(mass3)*e_tot)          ! momentum scale for the (b) tolerance
   w_tab = [ 0.8d0, 1.1d0, 1.4d0, 1.9d0, 2.3d0, 3.1d0 ]
   qz = non_geom()
   call build_table(qz, mass3, c_tab)          ! synthetic orthonormal table (self-checked)

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 3)
      call fill_list_atoms(3, qz, mass3)
      par%e_vib(1) = e_vib
      select case (trim(mode))
      case ('bad-evib')                     ! (f): e_vib = 0
         call set_scheme([6])
         par%e_vib(1) = 0.0d0
         call table_to_par(par, 1, c_tab, w_tab)
         call normalmode_init(par)
      case ('bad-nmode')                    ! (g): both tables empty
         call set_scheme([6])
         allocate (par%tbl(1))
         allocate (par%tbl(1)%w(0), par%tbl(1)%c(9, 0))
         call normalmode_init(par)
      case ('bad-mismatch')                 ! (h): 6 columns vs 7 frequencies
         call set_scheme([6])
         allocate (par%tbl(1))
         allocate (par%tbl(1)%w(7), par%tbl(1)%c(9, 6))
         par%tbl(1)%w = [ w_tab, 2.9d0 ]
         par%tbl(1)%c = c_tab
         call normalmode_init(par)
      case ('bad-shape')                    ! (i): 8 rows vs 3*nat = 9
         call set_scheme([6])
         allocate (par%tbl(1))
         allocate (par%tbl(1)%w(s_mode), par%tbl(1)%c(8, s_mode))
         par%tbl(1)%w = w_tab
         par%tbl(1)%c(1:8, :) = c_tab(1:8, :)
         call normalmode_init(par)
      case ('bad-omega')                    ! (j): tbl(1)%w(2) = 0
         call set_scheme([6])
         call table_to_par(par, 1, c_tab, w_tab)
         par%tbl(1)%w(2) = 0.0d0
         call normalmode_init(par)
      case ('bad-orth')                     ! (k): column 2 = column 1
         call set_scheme([6])
         call table_to_par(par, 1, c_tab, w_tab)
         par%tbl(1)%c(:, 2) = par%tbl(1)%c(:, 1)
         call normalmode_init(par)
      case ('bad-trans')                    ! (l): translation admixture
         call set_scheme([6])
         call tilt_column(par, c_tab, w_tab)
         call normalmode_init(par)
      case ('bad-bind0')                    ! (m): no carrier
         call set_scheme([1])
         call table_to_par(par, 1, c_tab, w_tab)
         call normalmode_init(par)
      case ('bad-bind2')                    ! (n): two carriers, one table
         call fill_two(qz, mass3)
         call set_scheme([6, 6])
         call table_to_par(par, 1, c_tab, w_tab)
         par%e_vib = [ e_vib, e_vib ]       ! two rows: the TABLE count is the mismatch
         call normalmode_init(par)
      case ('bad-noinit')                   ! (o): sample without init
         call set_scheme([6])
         call normalmode_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- (a)/(b)/(c) exact identities ----------
   call fill_list_atoms(3, qz, mass3)
   call set_scheme([6])
   call state_create(sta, 3)
   call table_to_par(par, 1, c_tab, w_tab)
   par%e_vib(1) = e_vib
   call rng_init(104729)
   call normalmode_init(par)
   ok_all = .true.
   do draw = 1, ndraw_exact
      call normalmode_sample(sta)
      call reexpand(sta, qz, mass3, c_tab, x_vec, y_vec, qmode, pmode, com_new)
      e_sum = 0.0d0
      do k = 1, s_mode
         e_sum = e_sum + 0.5d0*(pmode(k)**2 + w_tab(k)**2*qmode(k)**2)
      end do
      if (abs(e_sum - e_tot) > 1.0d-10*e_tot) ok_all = .false.                 ! (a)
      psum = sta%p(1:3) + sta%p(4:6) + sta%p(7:9)                              ! (b)
      if (norm2(psum) > 1.0d-9*psc) ok_all = .false.
      xres = norm2(x_vec - matmul(c_tab, qmode))
      yres = norm2(y_vec - matmul(c_tab, pmode))
      if (xres > 1.0d-10*max(norm2(x_vec), 1.0d-12)) ok_all = .false.
      if (yres > 1.0d-10*max(norm2(y_vec), 1.0d-12)) ok_all = .false.
      if (norm2(com_new) > 1.0d-10) ok_all = .false.                           ! (c)
   end do
   if (ok_all) then
      print '(a,es12.4,a)', '(a)-(c) energy/momentum/COM: PASS (E = e_vib*ev_dyn = ', e_tot, &
         ' internal exact 1e-10; |sum p| ~ 0, COM at origin, state in span(C))'
   else
      print '(a)', 'FAIL (a)-(c): energy identity, momentum, COM or subspace containment broke'
      nfail = nfail + 1
   end if
   call state_destroy(sta)

   ! ---------- (d)/(e) statistical criteria (fresh seeding) ----------
   call fill_list_atoms(3, qz, mass3)
   call set_scheme([6])
   call state_create(sta, 3)
   call table_to_par(par, 1, c_tab, w_tab)
   call rng_init(15485863)
   call normalmode_init(par)

   ! Beta(1/2,5/2) equal-CDF bin edges by bisection on the closed form
   gmax = 0.375d0*acos(-1.0d0)                 ! G(pi/2) = 3*pi/8
   do k = 1, nbin - 1
      lo = 0.0d0
      hi = 0.5d0*acos(-1.0d0)
      do it = 1, 80
         mid = 0.5d0*(lo + hi)
         if (g_theta(mid) < dble(k)*gmax/dble(nbin)) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x_edge(k) = sin(0.5d0*(lo + hi))**2
   end do

   freq = 0
   amean = 0.0d0
   kmean = 0.0d0
   do draw = 1, ndraw_stat
      call normalmode_sample(sta)
      call reexpand(sta, qz, mass3, c_tab, x_vec, y_vec, qmode, pmode, com_new)
      e1_share = 0.5d0*(pmode(1)**2 + w_tab(1)**2*qmode(1)**2)/e_tot
      amean = amean + e1_share
      ib = 1
      do k = 1, nbin - 1
         if (e1_share > x_edge(k)) ib = ib + 1
      end do
      freq(ib) = freq(ib) + 1
      do k = 1, s_mode
         kmean = kmean + 0.5d0*pmode(k)**2
      end do
   end do
   amean = amean/dble(ndraw_stat)
   kmean = kmean/dble(ndraw_stat)
   chi2 = 0.0d0
   do ib = 1, nbin
      chi2 = chi2 + (dble(freq(ib)) - ndraw_stat/dble(nbin))**2/(ndraw_stat/dble(nbin))
   end do
   emean = 1.0d0/6.0d0
   sig = sqrt((5.0d0/144.0d0)/dble(ndraw_stat))
   if (chi2 > 16.92d0) then
      print '(a,f9.3,a)', 'FAIL (d) Beta chi-square: chi2 = ', chi2, ' > 16.92 (df=9, alpha=0.05)'
      nfail = nfail + 1
   else if (abs(amean - emean) > 5.0d0*sig) then
      print '(a,f9.6,a,f9.6)', 'FAIL (d) Beta mean: sample ', amean, ' vs 1/6 = ', emean
      nfail = nfail + 1
   else
      print '(a,f7.3,a,f8.5,a)', '(d) mode-1 energy share: PASS (chi2 = ', chi2, &
         ' <= 16.92, df = 9; mean ', amean, ' vs Beta(1/2,5/2) mean 1/6; N = 50000)'
   end if
   sig = e_tot*sqrt(3.0d0/(64.0d0*dble(ndraw_stat)))
   if (abs(kmean - 0.5d0*e_tot) > 5.0d0*sig) then
      print '(a,es12.4,a,es12.4)', 'FAIL (e) equipartition: mean KE ', kmean, ' vs e_tot/2 = ', 0.5d0*e_tot
      nfail = nfail + 1
   else
      print '(a,es11.4,a)', '(e) phase equipartition: PASS (mean KE = ', kmean, &
         ' = e_tot/2 within 5 sigma; N = 50000)'
   end if
   call state_destroy(sta)

   ! ---------- (p) dual-carrier equal treatment: two fragments, two
   ! tables (the second a copy of the first), one init and one sampling
   ! loop - each carrier's atom block carries the full energy ----------
   call fill_two(qz, mass3)
   call set_scheme([6, 6])
   call state_create(sta, 6)
   call table_to_par(par, 2, c_tab, w_tab)
   par%e_vib = [ e_vib, e_vib ]      ! one energy, two rows
   call rng_init(611953)
   call normalmode_init(par)
   ok_all = .true.
   do draw = 1, ndraw_exact
      call normalmode_sample(sta)
      do k = 0, 1                                 ! both carriers' atom blocks
         if (abs(frag_energy(sta, qz, mass3, c_tab, 9*k) - e_tot) > 1.0d-10*e_tot) ok_all = .false.
      end do
   end do
   if (ok_all) then
      print '(a,es12.4,a)', '(p) dual-carrier equal treatment: PASS (each carrier block carries E = ', &
         e_tot, ' = e_vib*ev_dyn exactly; both fragments are written)'
   else
      print '(a)', 'FAIL (p) dual-carrier: a fragment block does not carry the full energy'
      nfail = nfail + 1
   end if
   call state_destroy(sta)

   ! ---------- (f)-(o) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-evib',     'normalmode_init', 'e_vib',       '(f) non-positive energy loud-stop')
   call expect_loudstop('bad-nmode',    'normalmode_init', 'mode count',  '(g) zero-mode-table loud-stop')
   call expect_loudstop('bad-mismatch', 'normalmode_init', 'mismatch',    '(h) table-length mismatch loud-stop')
   call expect_loudstop('bad-shape',    'normalmode_init', 'shape',       '(i) table-row shape loud-stop')
   call expect_loudstop('bad-omega',    'normalmode_init', 'frequency',   '(j) non-positive frequency loud-stop')
   call expect_loudstop('bad-orth',     'normalmode_init', 'orthonormal', '(k) non-orthonormal table loud-stop')
   call expect_loudstop('bad-trans',    'normalmode_init', 'translation', '(l) translation-violating table loud-stop')
   call expect_loudstop('bad-bind0',    'normalmode_init', 'normalmode scheme', '(m) no-carrier loud-stop')
   call expect_loudstop('bad-bind2',    'normalmode_init', 'carriers',    '(n) table-count/carrier-count mismatch loud-stop')
   call expect_loudstop('bad-noinit',   'normalmode_draw', 'normalmode_init', '(o) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_samp_normalmode: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! closed-form kernel of the Beta(1/2,5/2) CDF (monotone on [0,pi/2])
   real(8) function g_theta(th)
      real(8), intent(in) :: th
      g_theta = 0.75d0*th + 0.5d0*sin(2.0d0*th) + 0.0625d0*sin(4.0d0*th)
   end function g_theta

   ! bent three-atom geometry (non-principal, non-symmetric)
   function non_geom() result(g)
      real(8) :: g(9)
      g = 0.0d0
      g(1) = -0.7d0
      g(4) = 0.5d0
      g(8) = 0.45d0
      g(9) = 0.2d0
   end function non_geom

   ! synthetic mode table: Gram-Schmidt the canonical basis against the
   ! three mass-weighted translation directions; keep the first six
   ! orthonormal vectors of the complement (self-checked before use)
   subroutine build_table(g, mss, ctab)
      real(8), intent(in) :: g(9), mss(3)
      real(8), intent(out) :: ctab(9, s_mode)
      real(8) :: kept(9, 9), trans(9, 3), v(9), nrm, sm
      integer :: i, k, nkept, d
      trans = 0.0d0
      sm = sqrt(sum(mss))
      do d = 1, 3
         do i = 1, 3
            trans(3*(i - 1) + d, d) = sqrt(mss(i))
         end do
         trans(:, d) = trans(:, d)/sm
      end do
      kept(:, 1:3) = trans
      nkept = 3
      do k = 1, 9
         v = 0.0d0
         v(k) = 1.0d0
         do i = 1, nkept
            v = v - dot_product(v, kept(:, i))*kept(:, i)
         end do
         nrm = norm2(v)
         if (nrm > 1.0d-6) then
            nkept = nkept + 1
            kept(:, nkept) = v/nrm
         end if
         if (nkept == 9) exit
      end do
      ctab = kept(:, 4:9)
      do i = 1, s_mode
         do k = 1, s_mode
            if (abs(dot_product(ctab(:, i), ctab(:, k)) - merge(1.0d0, 0.0d0, i == k)) > 1.0d-10) &
               stop 'fixture: mode table not orthonormal'
         end do
         do d = 1, 3
            if (abs(dot_product(ctab(:, i), trans(:, d))) > 1.0d-10) &
               stop 'fixture: mode table not translation-orthogonal'
         end do
      end do
   end subroutine build_table

   ! re-expand the WRITTEN state into mass-weighted mode coordinates (the
   ! seed reference is the COM-shifted equilibrium geometry)
   subroutine reexpand(st, g, mss, ctab, xv, yv, qm, pm, cm)
      type(state_t), intent(in) :: st
      real(8), intent(in) :: g(9), mss(3), ctab(9, s_mode)
      real(8), intent(out) :: xv(9), yv(9), qm(s_mode), pm(s_mode), cm(3)
      real(8) :: qzc(9), cg(3)
      integer :: i
      cg = 0.0d0
      do i = 1, 3
         cg = cg + mss(i)*g(3*i-2:3*i)
      end do
      cg = cg/sum(mss)
      cm = 0.0d0
      do i = 1, 3
         cm = cm + mss(i)*st%q(3*i-2:3*i)
      end do
      cm = cm/sum(mss)
      do i = 1, 3
         qzc(3*i-2:3*i) = g(3*i-2:3*i) - cg
         xv(3*i-2:3*i) = sqrt(mss(i))*(st%q(3*i-2:3*i) - qzc(3*i-2:3*i))
         yv(3*i-2:3*i) = st%p(3*i-2:3*i)/sqrt(mss(i))
      end do
      qm = matmul(transpose(ctab), xv)
      pm = matmul(transpose(ctab), yv)
   end subroutine reexpand

   ! copy the synthetic table into the member parameters, one table per
   ! carrier (the second and further tables are copies of the same
   ! synthetic fixture); gfortran forbids allocating a parent object and
   ! its components in one statement, so the outer array is allocated
   ! first and the sibling components in a second statement
   subroutine table_to_par(p, n, ctab, wtab)
      type(normalmode_params_t), intent(inout) :: p
      integer, intent(in) :: n
      real(8), intent(in) :: ctab(9, s_mode), wtab(s_mode)
      integer :: j
      if (allocated(p%tbl)) deallocate (p%tbl)
      allocate (p%tbl(n))
      do j = 1, n
         allocate (p%tbl(j)%w(s_mode), p%tbl(j)%c(9, s_mode))
         p%tbl(j)%w = wtab
         p%tbl(j)%c = ctab
      end do
   end subroutine table_to_par

   ! column 1 tilted by a translation admixture (unit, still orthogonal to
   ! columns 2..6 - only the translation guard can catch it)
   subroutine tilt_column(p, ctab, wtab)
      type(normalmode_params_t), intent(inout) :: p
      real(8), intent(in) :: ctab(9, s_mode), wtab(s_mode)
      real(8) :: t1(9), sm, nrm
      integer :: i
      call table_to_par(p, 1, ctab, wtab)
      sm = sqrt(sum(mass3))
      t1 = 0.0d0
      do i = 1, 3
         t1(3*(i - 1) + 1) = sqrt(mass3(i))/sm
      end do
      p%tbl(1)%c(:, 1) = ctab(:, 1) + 0.1d0*t1
      nrm = norm2(p%tbl(1)%c(:, 1))
      p%tbl(1)%c(:, 1) = p%tbl(1)%c(:, 1)/nrm
   end subroutine tilt_column

   ! total mode energy re-expanded from one carrier's atom block (the (a)
   ! identity applied per carrier; off = 3*(first atom - 1))
   real(8) function frag_energy(st, g, mss, ctab, off)
      type(state_t), intent(in) :: st
      real(8), intent(in) :: g(9), mss(3), ctab(9, s_mode)
      integer, intent(in) :: off
      real(8) :: qzc(9), cg(3), xv(9), yv(9), qm(s_mode), pm(s_mode)
      integer :: i
      cg = 0.0d0
      do i = 1, 3
         cg = cg + mss(i)*g(3*i-2:3*i)
      end do
      cg = cg/sum(mss)
      do i = 1, 3
         qzc(3*i-2:3*i) = g(3*i-2:3*i) - cg
         xv(3*i-2:3*i) = sqrt(mss(i))*(st%q(off+3*i-2:off+3*i) - qzc(3*i-2:3*i))
         yv(3*i-2:3*i) = st%p(off+3*i-2:off+3*i)/sqrt(mss(i))
      end do
      qm = matmul(transpose(ctab), xv)
      pm = matmul(transpose(ctab), yv)
      frag_energy = 0.0d0
      do i = 1, s_mode
         frag_energy = frag_energy + 0.5d0*(pm(i)**2 + w_tab(i)**2*qm(i)**2)
      end do
   end function frag_energy

   ! mini-list_atoms by direct field fill: one fragment
   subroutine fill_list_atoms(nat, g, mss)
      integer, intent(in) :: nat
      real(8), intent(in) :: g(:), mss(:)
      integer :: k
      call list_atoms_clear()
      allocate (list_atoms%mass(nat))
      list_atoms%mass = mss
      allocate (list_atoms%frag(1))
      list_atoms%frag(1)%nat = nat
      list_atoms%frag(1)%mass = sum(mss)
      allocate (list_atoms%frag(1)%list(nat), list_atoms%frag(1)%qz(3*nat))
      list_atoms%frag(1)%list = [ (k, k = 1, nat) ]
      list_atoms%frag(1)%qz = g
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine fill_list_atoms

   subroutine fill_two(g, mss)
      real(8), intent(in) :: g(:), mss(:)
      integer :: k
      call list_atoms_clear()
      allocate (list_atoms%mass(6))
      list_atoms%mass = [ mss, mss ]
      allocate (list_atoms%frag(2))
      do k = 1, 2
         list_atoms%frag(k)%nat = 3
         list_atoms%frag(k)%mass = sum(mss)
         allocate (list_atoms%frag(k)%list(3), list_atoms%frag(k)%qz(9))
         list_atoms%frag(k)%list = [ 3*k - 2, 3*k - 1, 3*k ]
         list_atoms%frag(k)%qz = g
      end do
      reactants%surface_model = 0
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
   end subroutine fill_two

   ! preset the per-fragment scheme codes
   subroutine set_scheme(codes)
      integer, intent(in) :: codes(:)
      if (allocated(reactants%dist_scheme)) deallocate (reactants%dist_scheme)
      allocate (reactants%dist_scheme(size(codes)))
      reactants%dist_scheme = codes
   end subroutine set_scheme

   ! spawn the child for one bad case and assert the named abort
   subroutine expect_loudstop(cmode, who, key, label)
      character(len=*), intent(in) :: cmode, who, key, label
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_samp_normalmode.exe '//trim(cmode)// &
            ' 2>build'//sep//'normalmode_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/normalmode_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) line
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(line)
         end do
         close (v)
      end if
      if (stat == 0) then
         print '(a)', 'FAIL '//label//': expected a loud-stop, the child exited 0'
         nfail = nfail + 1
      else if (index(buf, trim(who)) == 0) then
         print '(a)', 'FAIL '//label//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (index(buf, trim(key)) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key)//'"'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop

   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_samp_normalmode
