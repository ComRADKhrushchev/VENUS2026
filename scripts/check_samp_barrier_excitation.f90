!=====================================================================
! check_samp_barrier_excitation.f90 - barrier_excitation distribution-
!                                     member criteria program (criteria
!                                     defined FIRST in this header; this
!                                     program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 &
!            globals/config_atoms.f90 globals/state.f90 &
!            utils/rng.f90 utils/geometry.f90 interface/sampler.f90 methods/samp_barrier_excitation.f90 &
!            ../scripts/check_samp_barrier_excitation.f90 -o build/check_samp_barrier_excitation.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter rule: the member is driven by directly
! constructed parameters and a directly filled list_atoms (a three-atom
! carrier, masses 1/2/3 amu at a non-symmetric geometry). Each
! carrier's saddle data is SYNTHETIC: one mass-weighted unit
! reaction-coordinate vector orthogonal to the three translations,
! plus five stable-mode columns
! (Gram-Schmidt of the canonical basis against the translations and the
! reaction coordinate) with arbitrary positive distinct frequencies
! w = [0.9, 1.3, 1.7, 2.2, 2.8] rad/(10 fs) - the member's contract is
! the measure on GIVEN saddle data (deriving it from a Hessian at a
! saddle is assembly-phase mechanism, not member code). Expected
! distributions are computed in this program from theory, not from the
! member.
! Mass-weighting convention (the member's, restated as the anchor; the
! member seeds the carrier about the COM origin):
!   x = sqrt(m)*(q - qz_com),  y = p/sqrt(m);  stable modes Q = C^T x,
!   P = C^T y;  the reaction-coordinate observables are d_rc = c_rc . x
!   (displacement along the rc, zero by construction) and p_rc =
!   c_rc . y (the barrier momentum, +/- sqrt(2*e_conv*e_bar) in mode 0,
!   sqrt(2*kb*t_bar*Exp(1)) in mode 1).
! Named-abort verification approach (same scheme as the sibling checks):
! every named-abort criterion runs as a CHILD PROCESS - the check re-spawns
! itself with a mode argument (execute_command_line); the child drives
! the REAL member entry into the stop and must print a naming message on
! stderr and exit nonzero. The parent asserts the exit status and the
! stderr content.
! Criteria (carrier: nat = 3, s = n_mode = 5, e_stab = 8 kcal/mol,
! e_bar = 5 kcal/mol (mode 0), t_bar = 900 K (mode 1); the stable-share
! marginal of the sphere split is Beta(1/2, 2)):
!   (a) exact battery, mode 0, 2000 draws: stable energy
!       sum_k 0.5*(P_k^2 + w_k^2*Q_k^2) = e_stab*e_conv to 1e-10;
!       |c_rc . y| = sqrt(2*e_conv*e_bar) to 1e-12; c_rc . x = 0 to
!       1e-12; |sum p| <= 1e-9*sqrt(2*M*e_tot_scale); COM at the origin
!       to 1e-10; x and y lie in span([C, c_rc]) (residual 1e-10)
!   (b) sign symmetry, mode 0, N = 20000: the fraction of positive
!       c_rc . y is 1/2 within 5 sigma (sigma = 0.5/sqrt(N) = 3.5e-3)
!   (c) thermal barrier momentum, mode 1, N = 50000: z = (c_rc . y)^2/2
!       ~ kb*t_bar*Exp(1); equal-CDF 10-bin chi-square (F(z) =
!       1 - exp(-z/(kb*t_bar))), df = 9, critical 16.92 (alpha = 0.05);
!       mean = kb*t_bar within 5 sigma (sigma = kb*t_bar/sqrt(N));
!       SIGN-INDEPENDENCE (defect regression): the per-sign mean of z
!       agrees between the positive and negative halves within
!       5*(2*kb*t_bar/sqrt(N)) - a magnitude-sign correlation fails this
!   (d) stable-share marginal, mode 1 draws: E_1/e_stab ~ Beta(1/2, 2);
!       equal-CDF 10-bin chi-square with closed-form quantiles
!       (CDF(s) = 1.5*s - 0.5*s^3, s = sqrt(x), solved by bisection),
!       df = 9, critical 16.92; mean = 1/5 within 5 sigma
!       (Var = 8/175, sigma of the mean = sqrt(Var/N))
!   (e) negative stabilization energy named abort (child): e_stab = -1 ->
!       stderr names barrier_excitation_init and e_stab
!   (f) bad energy-law selector named abort (child): n_e_bar = 2 ->
!       stderr names barrier_excitation_init and n_e_bar
!   (g) negative barrier energy named abort (child): mode 0, e_bar = -2 ->
!       stderr names barrier_excitation_init and e_bar
!   (h) non-positive barrier temperature named abort (child): mode 1,
!       t_bar = 0 -> stderr names barrier_excitation_init and t_bar
!   (i) non-positive stable frequency named abort (child): tbl(1)%w(3) = 0
!       -> stderr names barrier_excitation_init and the frequency
!   (j) non-orthonormal stable table named abort (child): tbl(1)%c column
!       2 replaced by column 1 -> stderr names barrier_excitation_init
!       and orthonormality
!   (k) non-unit reaction coordinate named abort (child): tbl(1)%rc
!       scaled by 1.1 -> stderr names barrier_excitation_init and
!       "not unit"
!   (l) translation-violating reaction coordinate named abort (child):
!       tbl(1)%rc tilted by a translation admixture -> stderr names
!       barrier_excitation_init and the translation
!   (m) stable-mode/reaction-coordinate overlap named abort (child):
!       tbl(1)%rc mixed with stable column 1 -> stderr names
!       barrier_excitation_init and the overlap
!   (n) stabilization energy without stable modes named abort (child):
!       e_stab = 6 with an empty stable table -> stderr names
!       barrier_excitation_init and the stable support
!   (o) table-length mismatch named abort (child): tbl(1)%c 5 columns vs
!       tbl(1)%w 6 entries -> stderr names barrier_excitation_init and
!       the mismatch
!   (p) reaction-coordinate length named abort (child): tbl(1)%rc of
!       length 8 -> stderr names barrier_excitation_init and the length
!   (q) no-fragment-binding named abort (child): stderr names
!       barrier_excitation_init and the barrier_excitation scheme
!   (r) table-count mismatch named abort (child): TWO fragments select
!       the scheme but the parameter set carries ONE table -> stderr
!       names barrier_excitation_init and the carrier count
!   (s) sample-before-init named abort (child): stderr names
!       barrier_excitation_sample
!   (t) dual-carrier equal treatment (positive): two fragments, two
!       saddle tables (the second a copy of the first), one init - each
!       carrier's rc momentum observable c_rc . y = +/- sqrt(2*e_conv*e_bar)
!=====================================================================
program check_samp_barrier_excitation
   use consts,        only: kb_code, e_conv
   use rng,           only: rng_init
   use config,        only: reactants
   use config_atoms, only: list_atoms, list_atoms_clear
   use state,         only: state_t, state_create, state_destroy
   use samp_barrier_excitation, only: barrier_excitation_params_t, barrier_excitation_init, &
                                      barrier_excitation_sample
   implicit none
   type(state_t) :: sta
   type(barrier_excitation_params_t) :: par
   integer, parameter :: ndraw_exact = 2000
   integer, parameter :: ndraw_sign = 20000
   integer, parameter :: ndraw_stat = 50000
   integer, parameter :: nbin = 10
   integer, parameter :: s_mode = 5
   integer :: nfail = 0, draw, ib, k, it
   integer :: freq(nbin)
   integer :: n_plus, n_minus
   real(8) :: e_stab, e_bar, pbar0, e_tot_s, psc
   real(8) :: qz(9), mass3(3), com_new(3)
   real(8) :: w_tab(s_mode), c_tab(9, s_mode), rc_tab(9)
   real(8) :: x_vec(9), y_vec(9), qmode(s_mode), pmode(s_mode)
   real(8) :: d_rc, p_rc, e_sum, z_val, amean, emean, sig, chi2, cdf
   real(8) :: zmean_p, zmean_m, e1_share, s_edge(nbin - 1)
   real(8) :: lo, hi, mid, sval
   character(len=32) :: mode
   logical :: ok_all

   e_stab = 8.0d0
   e_bar = 5.0d0
   e_tot_s = e_stab*e_conv
   pbar0 = sqrt(2.0d0*e_bar*e_conv)
   mass3 = [ 1.0d0, 2.0d0, 3.0d0 ]
   psc = sqrt(2.0d0*sum(mass3)*(e_tot_s + pbar0**2/2.0d0))
   w_tab = [ 0.9d0, 1.3d0, 1.7d0, 2.2d0, 2.8d0 ]
   qz = non_geom()
   call build_saddle(qz, mass3, rc_tab, c_tab)   ! synthetic saddle data (self-checked)

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call state_create(sta, 3)
      call fill_list_atoms(3, qz, mass3)
      par%e_stab = e_stab
      par%e_bar = e_bar
      par%t_bar = 900.0d0
      select case (trim(mode))
      case ('bad-estab')                  ! (e): e_stab = -1
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         par%e_stab = -1.0d0
         call barrier_excitation_init(par)
      case ('bad-nbar')                   ! (f): n_e_bar = 2
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         par%n_e_bar = 2
         call barrier_excitation_init(par)
      case ('bad-ebar')                   ! (g): mode 0, e_bar = -2
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         par%e_bar = -2.0d0
         call barrier_excitation_init(par)
      case ('bad-tbar')                   ! (h): mode 1, t_bar = 0
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         par%n_e_bar = 1
         par%t_bar = 0.0d0
         call barrier_excitation_init(par)
      case ('bad-omega')                  ! (i): tbl(1)%w(3) = 0
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         par%tbl(1)%w(3) = 0.0d0
         call barrier_excitation_init(par)
      case ('bad-orth')                   ! (j): column 2 = column 1
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         par%tbl(1)%c(:, 2) = par%tbl(1)%c(:, 1)
         call barrier_excitation_init(par)
      case ('bad-rcnorm')                 ! (k): tbl(1)%rc scaled by 1.1
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         par%tbl(1)%rc = 1.1d0*par%tbl(1)%rc
         call barrier_excitation_init(par)
      case ('bad-rctrans')                ! (l): translation admixture in c_rc
         call set_scheme([7])
         call tilt_rc(par, rc_tab, c_tab, w_tab)
         call barrier_excitation_init(par)
      case ('bad-rcperp')                 ! (m): c_rc mixed with stable column 1
         call set_scheme([7])
         call mix_rc(par, rc_tab, c_tab, w_tab)
         call barrier_excitation_init(par)
      case ('bad-nostab')                 ! (n): e_stab > 0, empty stable table
         call set_scheme([7])
         allocate (par%tbl(1))
         allocate (par%tbl(1)%w(0), par%tbl(1)%c(9, 0), par%tbl(1)%rc(9))
         par%tbl(1)%rc = rc_tab
         par%e_stab = 6.0d0
         call barrier_excitation_init(par)
      case ('bad-mismatch')               ! (o): 5 columns vs 6 frequencies
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         deallocate (par%tbl(1)%w)
         allocate (par%tbl(1)%w(6))
         par%tbl(1)%w = [ w_tab, 3.3d0 ]
         call barrier_excitation_init(par)
      case ('bad-rclen')                  ! (p): tbl(1)%rc of length 8
         call set_scheme([7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         deallocate (par%tbl(1)%rc)
         allocate (par%tbl(1)%rc(8))
         par%tbl(1)%rc = rc_tab(1:8)
         call barrier_excitation_init(par)
      case ('bad-bind0')                  ! (q): no carrier
         call set_scheme([1])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         call barrier_excitation_init(par)
      case ('bad-bind2')                  ! (r): two carriers, one table
         call fill_two(qz, mass3)
         call set_scheme([7, 7])
         call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
         call barrier_excitation_init(par)
      case ('bad-noinit')                 ! (s): sample without init
         call set_scheme([7])
         call barrier_excitation_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- (a) exact battery, mode 0 ----------
   call fill_list_atoms(3, qz, mass3)
   call set_scheme([7])
   call state_create(sta, 3)
   call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
   par%e_stab = e_stab
   par%e_bar = e_bar
   call rng_init(104729)
   call barrier_excitation_init(par)
   ok_all = .true.
   do draw = 1, ndraw_exact
      call barrier_excitation_sample(sta)
      call reexpand(sta, qz, mass3, c_tab, rc_tab, x_vec, y_vec, qmode, pmode, com_new)
      e_sum = 0.0d0
      do k = 1, s_mode
         e_sum = e_sum + 0.5d0*(pmode(k)**2 + w_tab(k)**2*qmode(k)**2)
      end do
      if (abs(e_sum - e_tot_s) > 1.0d-10*e_tot_s) ok_all = .false.         ! stable energy
      p_rc = dot_product(rc_tab, y_vec)
      d_rc = dot_product(rc_tab, x_vec)
      if (abs(abs(p_rc) - pbar0) > 1.0d-12*pbar0) ok_all = .false.         ! barrier momentum
      if (abs(d_rc) > 1.0d-12) ok_all = .false.                            ! zero rc displacement
      if (norm2(y_vec - matmul(c_tab, pmode) - p_rc*rc_tab) > 1.0d-10*max(norm2(y_vec), 1.0d-12)) &
         ok_all = .false.                                                   ! y in span([C, rc])
      if (norm2(x_vec - matmul(c_tab, qmode)) > 1.0d-10*max(norm2(x_vec), 1.0d-12)) &
         ok_all = .false.                                                   ! x in span(C)
      if (norm2(sta%p(1:3) + sta%p(4:6) + sta%p(7:9)) > 1.0d-9*psc) ok_all = .false.
      if (norm2(com_new) > 1.0d-10) ok_all = .false.
   end do
   if (ok_all) then
      print '(a,es11.4,a)', '(a) mode-0 exact battery: PASS (E_stab = ', e_tot_s, &
         ' internal, |p_rc| = sqrt(2*e_conv*e_bar), d_rc = 0; P, COM, span identities)'
   else
      print '(a)', 'FAIL (a) mode-0 exact battery: an energy/momentum/displacement identity broke'
      nfail = nfail + 1
   end if
   call state_destroy(sta)

   ! ---------- (b) sign symmetry, mode 0 (fresh seeding) ----------
   call fill_list_atoms(3, qz, mass3)
   call set_scheme([7])
   call state_create(sta, 3)
   call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
   call rng_init(32452843)
   call barrier_excitation_init(par)
   n_plus = 0
   do draw = 1, ndraw_sign
      call barrier_excitation_sample(sta)
      call reexpand(sta, qz, mass3, c_tab, rc_tab, x_vec, y_vec, qmode, pmode, com_new)
      if (dot_product(rc_tab, y_vec) > 0.0d0) n_plus = n_plus + 1
   end do
   n_minus = ndraw_sign - n_plus
   sig = 0.5d0/sqrt(dble(ndraw_sign))
   if (abs(dble(n_plus)/dble(ndraw_sign) - 0.5d0) > 5.0d0*sig) then
      print '(a,i0,a,i0)', 'FAIL (b) sign symmetry: ', n_plus, ' positive vs ', n_minus, ' negative'
      nfail = nfail + 1
   else
      print '(a,i0,a,i0,a)', '(b) sign symmetry: PASS (', n_plus, ' positive / ', n_minus, &
         ' negative, N = 20000, 5 sigma band)'
   end if
   call state_destroy(sta)

   ! ---------- (c)/(d) thermal barrier momentum + stable share ----------
   call fill_list_atoms(3, qz, mass3)
   call set_scheme([7])
   call state_create(sta, 3)
   call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
   par%n_e_bar = 1
   par%t_bar = 900.0d0
   call rng_init(15485863)
   call barrier_excitation_init(par)

   ! Beta(1/2,2) equal-CDF bin edges: CDF(s) = 1.5*s - 0.5*s^3, s = sqrt(x)
   do k = 1, nbin - 1
      lo = 0.0d0
      hi = 1.0d0
      do it = 1, 80
         mid = 0.5d0*(lo + hi)
         if (cdf_beta2(mid) < dble(k)/dble(nbin)) then
            lo = mid
         else
            hi = mid
         end if
      end do
      s_edge(k) = 0.5d0*(lo + hi)
   end do

   freq = 0
   amean = 0.0d0
   zmean_p = 0.0d0
   zmean_m = 0.0d0
   n_plus = 0
   n_minus = 0
   do draw = 1, ndraw_stat
      call barrier_excitation_sample(sta)
      call reexpand(sta, qz, mass3, c_tab, rc_tab, x_vec, y_vec, qmode, pmode, com_new)
      p_rc = dot_product(rc_tab, y_vec)
      z_val = 0.5d0*p_rc*p_rc
      amean = amean + z_val
      if (p_rc > 0.0d0) then
         n_plus = n_plus + 1
         zmean_p = zmean_p + z_val
      else
         n_minus = n_minus + 1
         zmean_m = zmean_m + z_val
      end if
      cdf = 1.0d0 - exp(-z_val/(kb_code*900.0d0))
      ib = min(nbin, max(1, int(cdf*dble(nbin)) + 1))
      freq(ib) = freq(ib) + 1
   end do
   amean = amean/dble(ndraw_stat)
   chi2 = 0.0d0
   do ib = 1, nbin
      chi2 = chi2 + (dble(freq(ib)) - ndraw_stat/dble(nbin))**2/(ndraw_stat/dble(nbin))
   end do
   emean = kb_code*900.0d0
   sig = emean/sqrt(dble(ndraw_stat))
   zmean_p = zmean_p/dble(n_plus)
   zmean_m = zmean_m/dble(n_minus)
   if (chi2 > 16.92d0) then
      print '(a,f9.3,a)', 'FAIL (c) thermal chi-square: chi2 = ', chi2, ' > 16.92 (df=9, alpha=0.05)'
      nfail = nfail + 1
   else if (abs(amean - emean) > 5.0d0*sig) then
      print '(a,es12.4,a,es12.4)', 'FAIL (c) thermal mean: sample ', amean, ' vs kb*T = ', emean
      nfail = nfail + 1
   else if (abs(zmean_p - zmean_m) > 5.0d0*2.0d0*emean/sqrt(dble(ndraw_stat))) then
      print '(a,2es12.4)', 'FAIL (c) sign independence: per-sign means ', zmean_p, zmean_m
      nfail = nfail + 1
   else
      print '(a,f7.3,a,es11.4,a)', '(c) thermal barrier momentum: PASS (chi2 = ', chi2, &
         ' <= 16.92, df = 9; mean ', amean, ' = kb*T; per-sign means agree; N = 50000)'
   end if
   call state_destroy(sta)

   ! (d) stable share, separate pass for clean binning
   call fill_list_atoms(3, qz, mass3)
   call set_scheme([7])
   call state_create(sta, 3)
   call saddle_to_par(par, 1, rc_tab, c_tab, w_tab)
   par%n_e_bar = 1
   par%t_bar = 900.0d0
   call rng_init(49979687)
   call barrier_excitation_init(par)
   freq = 0
   amean = 0.0d0
   do draw = 1, ndraw_stat
      call barrier_excitation_sample(sta)
      call reexpand(sta, qz, mass3, c_tab, rc_tab, x_vec, y_vec, qmode, pmode, com_new)
      e1_share = 0.5d0*(pmode(1)**2 + w_tab(1)**2*qmode(1)**2)/e_tot_s
      amean = amean + e1_share
      sval = sqrt(max(e1_share, 0.0d0))
      ib = 1
      do k = 1, nbin - 1
         if (sval > s_edge(k)) ib = ib + 1
      end do
      freq(ib) = freq(ib) + 1
   end do
   amean = amean/dble(ndraw_stat)
   chi2 = 0.0d0
   do ib = 1, nbin
      chi2 = chi2 + (dble(freq(ib)) - ndraw_stat/dble(nbin))**2/(ndraw_stat/dble(nbin))
   end do
   emean = 0.2d0
   sig = sqrt((8.0d0/175.0d0)/dble(ndraw_stat))
   if (chi2 > 16.92d0 .or. abs(amean - emean) > 5.0d0*sig) then
      print '(a,f9.3,a,f9.6)', 'FAIL (d) stable share: chi2 = ', chi2, ' mean = ', amean
      nfail = nfail + 1
   else
      print '(a,f7.3,a,f8.5,a)', '(d) stable-mode-1 energy share: PASS (chi2 = ', chi2, &
         ' <= 16.92, df = 9; mean ', amean, ' vs Beta(1/2,2) mean 1/5; N = 50000)'
   end if
   call state_destroy(sta)

   ! ---------- (t) dual-carrier equal treatment: two fragments, two
   ! saddle tables (the second a copy of the first), one init - each
   ! carrier's rc momentum observable carries the barrier kick ------
   call fill_two(qz, mass3)
   call set_scheme([7, 7])
   call state_create(sta, 6)
   call saddle_to_par(par, 2, rc_tab, c_tab, w_tab)
   par%e_stab = e_stab
   par%n_e_bar = 0
   par%e_bar = e_bar
   call rng_init(611953)
   call barrier_excitation_init(par)
   ok_all = .true.
   do draw = 1, ndraw_exact
      call barrier_excitation_sample(sta)
      do k = 0, 1                                 ! both carriers' atom blocks
         if (abs(abs(frag_rc_p(sta, mass3, 9*k, rc_tab)) - pbar0) > 1.0d-12*pbar0) ok_all = .false.
      end do
   end do
   if (ok_all) then
      print '(a,es12.4,a)', '(t) dual-carrier equal treatment: PASS (|c_rc . y| = ', pbar0, &
         ' = sqrt(2*e_conv*e_bar) on both carriers)'
   else
      print '(a)', 'FAIL (t) dual-carrier: a carrier block missed the barrier momentum'
      nfail = nfail + 1
   end if
   call state_destroy(sta)

   ! ---------- (e)-(s) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-estab',   'barrier_excitation_init', 'e_stab',  '(e) negative stabilization energy loud-stop')
   call expect_loudstop('bad-nbar',    'barrier_excitation_init', 'n_e_bar', '(f) bad energy-law selector loud-stop')
   call expect_loudstop('bad-ebar',    'barrier_excitation_init', 'e_bar',   '(g) negative barrier energy loud-stop')
   call expect_loudstop('bad-tbar',    'barrier_excitation_init', 't_bar',   '(h) non-positive barrier temperature loud-stop')
   call expect_loudstop('bad-omega',   'barrier_excitation_init', 'frequency', '(i) non-positive stable frequency loud-stop')
   call expect_loudstop('bad-orth',    'barrier_excitation_init', 'orthonormal', '(j) non-orthonormal stable table loud-stop')
   call expect_loudstop('bad-rcnorm',  'barrier_excitation_init', 'not unit', '(k) non-unit reaction coordinate loud-stop')
   call expect_loudstop('bad-rctrans', 'barrier_excitation_init', 'translation', '(l) translation-violating rc loud-stop')
   call expect_loudstop('bad-rcperp',  'barrier_excitation_init', 'overlap', '(m) stable-mode/rc overlap loud-stop')
   call expect_loudstop('bad-nostab',  'barrier_excitation_init', 'stable',  '(n) energy-without-support loud-stop')
   call expect_loudstop('bad-mismatch','barrier_excitation_init', 'mismatch', '(o) table-length mismatch loud-stop')
   call expect_loudstop('bad-rclen',   'barrier_excitation_init', 'length',  '(p) reaction-coordinate length loud-stop')
   call expect_loudstop('bad-bind0',   'barrier_excitation_init', 'barrier_excitation scheme', '(q) no-carrier loud-stop')
   call expect_loudstop('bad-bind2',   'barrier_excitation_init', 'carriers', '(r) table-count/carrier-count mismatch loud-stop')
   call expect_loudstop('bad-noinit',  'barrier_excitation_draw', 'barrier_excitation_init', '(s) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_samp_barrier_excitation: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! closed-form Beta(1/2,2) CDF in s = sqrt(x)
   real(8) function cdf_beta2(sv)
      real(8), intent(in) :: sv
      cdf_beta2 = 1.5d0*sv - 0.5d0*sv**3
   end function cdf_beta2

   ! bent three-atom geometry (non-principal, non-symmetric)
   function non_geom() result(g)
      real(8) :: g(9)
      g = 0.0d0
      g(1) = -0.7d0
      g(4) = 0.5d0
      g(8) = 0.45d0
      g(9) = 0.2d0
   end function non_geom

   ! synthetic saddle data: unit reaction coordinate orthogonal to the
   ! translations, then five stable columns Gram-Schmidted against the
   ! translations and the reaction coordinate (self-checked before use)
   subroutine build_saddle(g, mss, rct, ctab)
      real(8), intent(in) :: g(9), mss(3)
      real(8), intent(out) :: rct(9), ctab(9, s_mode)
      real(8) :: trans(9, 3), kept(9, 9), v(9), nrm, sm
      integer :: i, k, nkept, d
      trans = 0.0d0
      sm = sqrt(sum(mss))
      do d = 1, 3
         do i = 1, 3
            trans(3*(i - 1) + d, d) = sqrt(mss(i))
         end do
         trans(:, d) = trans(:, d)/sm
      end do
      ! reaction coordinate: canonical e_5 stripped of translations
      v = 0.0d0
      v(5) = 1.0d0
      do d = 1, 3
         v = v - dot_product(v, trans(:, d))*trans(:, d)
      end do
      rct = v/norm2(v)
      kept(:, 1:3) = trans
      kept(:, 4) = rct
      nkept = 4
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
      ctab = kept(:, 5:9)
      ! fixture sanity
      if (abs(norm2(rct) - 1.0d0) > 1.0d-12) stop 'fixture: rc not unit'
      do i = 1, s_mode
         if (abs(dot_product(ctab(:, i), rct)) > 1.0d-10) stop 'fixture: stable column not rc-orthogonal'
         do d = 1, 3
            if (abs(dot_product(ctab(:, i), trans(:, d))) > 1.0d-10) &
               stop 'fixture: stable column not translation-orthogonal'
         end do
         do k = 1, s_mode
            if (abs(dot_product(ctab(:, i), ctab(:, k)) - merge(1.0d0, 0.0d0, i == k)) > 1.0d-10) &
               stop 'fixture: stable table not orthonormal'
         end do
      end do
   end subroutine build_saddle

   ! re-expand the WRITTEN state into mass-weighted saddle coordinates
   subroutine reexpand(st, g, mss, ctab, rct, xv, yv, qm, pm, cm)
      type(state_t), intent(in) :: st
      real(8), intent(in) :: g(9), mss(3), ctab(9, s_mode), rct(9)
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

   ! copy the synthetic saddle data into the member parameters, one table
   ! per carrier (the second and further tables are copies of the same
   ! synthetic fixture); gfortran forbids allocating a parent object and
   ! its components in one statement, so the outer array is allocated
   ! first and the sibling components in a second statement
   subroutine saddle_to_par(p, n, rct, ctab, wtab)
      type(barrier_excitation_params_t), intent(inout) :: p
      integer, intent(in) :: n
      real(8), intent(in) :: rct(9), ctab(9, s_mode), wtab(s_mode)
      integer :: j
      if (allocated(p%tbl)) deallocate (p%tbl)
      allocate (p%tbl(n))
      do j = 1, n
         allocate (p%tbl(j)%w(s_mode), p%tbl(j)%c(9, s_mode), p%tbl(j)%rc(9))
         p%tbl(j)%w = wtab
         p%tbl(j)%c = ctab
         p%tbl(j)%rc = rct
      end do
   end subroutine saddle_to_par

   ! c_rc tilted by a translation admixture (still unit - only the
   ! translation guard can catch it)
   subroutine tilt_rc(p, rct, ctab, wtab)
      type(barrier_excitation_params_t), intent(inout) :: p
      real(8), intent(in) :: rct(9), ctab(9, s_mode), wtab(s_mode)
      real(8) :: t1(9), sm
      integer :: i
      call saddle_to_par(p, 1, rct, ctab, wtab)
      sm = sqrt(sum(mass3))
      t1 = 0.0d0
      do i = 1, 3
         t1(3*(i - 1) + 1) = sqrt(mass3(i))/sm
      end do
      p%tbl(1)%rc = rct + 0.15d0*t1
      p%tbl(1)%rc = p%tbl(1)%rc/norm2(p%tbl(1)%rc)
   end subroutine tilt_rc

   ! c_rc mixed with stable column 1, renormalized (violates the rc-stable
   ! orthogonality only)
   subroutine mix_rc(p, rct, ctab, wtab)
      type(barrier_excitation_params_t), intent(inout) :: p
      real(8), intent(in) :: rct(9), ctab(9, s_mode), wtab(s_mode)
      call saddle_to_par(p, 1, rct, ctab, wtab)
      p%tbl(1)%rc = rct + 0.2d0*ctab(:, 1)
      p%tbl(1)%rc = p%tbl(1)%rc/norm2(p%tbl(1)%rc)
   end subroutine mix_rc

   ! rc momentum observable of one carrier's atom block (the (a) c_rc . y
   ! identity applied per carrier; off = 3*(first atom - 1))
   real(8) function frag_rc_p(st, mss, off, rct)
      type(state_t), intent(in) :: st
      real(8), intent(in) :: mss(3), rct(9)
      integer, intent(in) :: off
      real(8) :: yv(9)
      integer :: i
      do i = 1, 3
         yv(3*i-2:3*i) = st%p(off+3*i-2:off+3*i)/sqrt(mss(i))
      end do
      frag_rc_p = dot_product(rct, yv)
   end function frag_rc_p

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
      cmd = 'build'//sep//'check_samp_barrier_excitation.exe '//trim(cmode)// &
            ' 2>build'//sep//'barrier_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/barrier_err.txt', status='old', action='read', iostat=ios)
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

end program check_samp_barrier_excitation
