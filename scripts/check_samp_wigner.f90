!=====================================================================
! check_samp_wigner.f90 - wigner distribution-member criteria program
!                          (criteria defined FIRST in this header; this
!                          program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 interface/sysdef_interface.f90 input.f90 &
!            globals/config_atoms.f90 globals/state.f90 &
!            utils/linalg.f90 utils/hessian.f90 utils/geom_init.f90 utils/spectrum.f90 interface/elec_interface.f90 interface/force_interface.f90 dynamics/propagator.f90 interface/bath.f90 &
!            utils/rng.f90 utils/geometry.f90 interface/sampler.f90 methods/samp_wigner.f90 &
!            ../scripts/check_samp_wigner.f90 -o build/check_samp_wigner.exe
! All artifacts are written under build/ - the source tree stays clean.
! Synthetic-parameter rule: the member is driven by directly constructed
! parameters and a directly filled list_atoms (a three-atom carrier, masses
! 1/2/3 amu at a non-symmetric geometry). The mode table is SYNTHETIC:
! six orthonormal mass-weighted columns spanning the complement of the
! three translations (built here by Gram-Schmidt from the canonical
! basis), with frequencies w = [0.8, 1.1, 1.4, 1.9, 2.3, 3.1] rad/(10 fs)
! - the member's contract is the measure on a GIVEN table (frequency
! derivation from a Hessian is assembly-phase mechanism, not member code).
! Expected distributions are computed here from theory, not from the
! member output. Theory anchor (the ground-state Wigner law, energies in
! dynamical units - hbar_evfs*ev_dyn, the member's own fold):
!   per mode E = (hbar/2)*w; sigma_p = sqrt(E), sigma_q = sqrt(E)/w;
!   Q, P independent Gaussians N(0, sigma^2); per-mode energy
!   E*(g1^2+g2^2)/2 exponential with mean E, so <T> = <V> = E/2,
!   <Q^2> = E/w^2, <P^2> = E. Numbers for this table:
!   E = [0.02540, 0.03493, 0.04446, 0.06033, 0.07303, 0.09844],
!   sum E = 0.33659, sum E^2 = 0.022506 (internal units).
! Mass-weighting convention (the member's, restated as the anchor):
!   x = sqrt(m)*(q - qz_com),  y = p/sqrt(m),  Q = C^T x,  P = C^T y;
!   conversely q = qz_com + (C Q)/sqrt(m), p = sqrt(m)*(C P).
! Named-abort verification approach (same scheme as the sibling checks):
! every named-abort criterion runs as a CHILD PROCESS - the check
! re-spawns itself with a mode argument (execute_command_line); the
! child drives the REAL member entry into the stop and must print a
! naming message on stderr and exit nonzero. The parent asserts the
! exit status and the stderr content.
! Criteria (carrier: nat = 3, s = n_mode = 6, e_scale = 1; N = 200000
! unless noted):
!   (a) mean total mode energy = sum E_k within 5 sigma of the
!       exponential-sum law (Var per draw = sum E_k^2, so
!       5 sigma = 5*sqrt(sum E_k^2/N) = 1.68e-3)
!   (b) per-mode variances: <Q_k^2> = E_k/w_k^2 AND <P_k^2> = E_k for
!       every mode, each within 5 sigma = 5*sigma^2*sqrt(2/N)
!       (variance of a Gaussian square is 2 sigma^4)
!   (c) q-p decorrelation: |corr(Q_k, P_k)| <= 5/sqrt(N) = 0.0112 for
!       every mode (coordinate and momentum are independent Gaussians)
!   (d) virial in the mean: |<T_k> - <V_k>| <= 5*E_k/sqrt(N) for every
!       mode (Var(T - V) = E_k^2, T and V independent)
!   (e) Gaussian shape on mode 1: |skewness| <= 5*sqrt(6/N) = 0.0274
!       and |kurtosis - 3| <= 5*sqrt(24/N) = 0.0548
!   (f) structural exacts, every draw: net momentum <=
!       1e-9*sqrt(2*M*sum E) (M = carrier mass); the COM of the written
!       q sits at the ORIGIN (COM-origin output) to 1e-10 Angstrom; the
!       residual of sqrt(m)(q - qz_com)
!       and of p/sqrt(m) after projection onto span(C) is <= 1e-10 of
!       the vector norm (only mode-subspace content is written)
!   (g) e_scale plumbing: e_scale = 0.5 -> <Q_1^2> = 0.5*E_1/w_1^2 =
!       0.019846 within 5 sigma (the shared scale reaches every mode)
!   (h) dual-carrier equal treatment: two fragments (6-mode and 3-mode
!       ragged tables) at ONE e_scale -> each carrier's mode-1 <Q^2>
!       matches its own E_1/w_1^2 (carrier 2: w = [0.9, 1.6, 2.7],
!       <Q_1^2> = 0.035282) - the carriers do not cross-talk
!   (i) non-positive e_scale named abort (child): e_scale = 0 ->
!       stderr names wigner_init and the field
!   (j) no-fragment-binding named abort (child): dist_scheme carries no
!       wigner row -> stderr names wigner_init and the member word
!   (k) unallocated-table-set named abort (child): tbl unallocated ->
!       stderr names wigner_init and the mismatch
!   (l) table-count mismatch named abort (child): TWO fragments carry
!       the wigner code but the parameter set carries ONE table ->
!       stderr names wigner_init and the carrier count
!   (m) inner-table unallocated named abort (child): the table set
!       carries an unallocated eigenvector component -> stderr names
!       wigner_init and the unallocated table
!   (n) non-positive frequency named abort (child): a 0 rad/(10 fs)
!       table entry -> stderr names wigner_init and the frequency
!   (o) non-orthonormal table named abort (child): one eigenvector
!       entry perturbed by 1e-3 -> stderr names wigner_init and the
!       orthonormality
!   (p) sample-before-init named abort (child): wigner_sample called
!       with no wigner_init -> stderr names wigner_draw (the two-phase guard site)
!=====================================================================
program check_samp_wigner
   use consts,         only: hbar_evfs, ev_dyn
   use rng,            only: rng_init
   use config,         only: reactants
   use config_atoms,  only: list_atoms, list_atoms_clear
   use state,          only: state_t, state_create, state_destroy
   use samp_wigner,    only: wigner_params_t, wigner_init, wigner_sample
   implicit none
   type(state_t) :: sta
   type(wigner_params_t) :: par
   integer, parameter :: ndraw = 200000
   integer, parameter :: nmode1 = 6                 ! carrier-1 mode count
   integer, parameter :: nmode2 = 3                 ! carrier-2 mode count (dual leg)
   integer :: nfail = 0, nfail0, draw, k
   integer, parameter :: code_w = 10                ! the wigner scheme code
   real(8), parameter :: w1(nmode1) = [ 0.8d0, 1.1d0, 1.4d0, 1.9d0, 2.3d0, 3.1d0 ]  ! rad/(10 fs)
   real(8), parameter :: w2(nmode2) = [ 0.9d0, 1.6d0, 2.7d0 ]                      ! rad/(10 fs)
   real(8), parameter :: mss1(3) = [ 1.0d0, 2.0d0, 3.0d0 ]   ! carrier-1 masses [amu]
   real(8), parameter :: mss2(2) = [ 2.0d0, 5.0d0 ]          ! carrier-2 masses [amu]
   real(8), parameter :: g1(9) = [ 0.0d0, 0.1d0, -0.2d0, 1.5d0, 0.3d0, -0.1d0, -0.7d0, 1.9d0, 0.5d0 ]
   real(8), parameter :: g2(6) = [ 0.0d0, 0.0d0, 0.0d0, 1.1d0, 0.4d0, 0.2d0 ]
   real(8) :: e1(nmode1), e2(nmode2), sig_q2_1(nmode1), sig_q2_2(nmode2)
   real(8), allocatable :: ctab1(:,:), ctab2(:,:)
   real(8), allocatable :: qm(:), pm(:)
   real(8) :: esum, esq_sum, emean, etol
   real(8) :: q2(nmode1), p2(nmode1), qp(nmode1), tbar(nmode1), vbar(nmode1)
   real(8) :: m1, mm2, mm3, mm4, mu2, mu3, mu4, skw, krt, corr, tol
   real(8) :: pmom_max, dcom_max, resq_max, resp_max, x(9), y(9), vres
   real(8) :: pmom(3), dcom(3)
   integer :: list1(3) = [ 1, 2, 3 ], list2(2) = [ 4, 5 ]
   character(len=32) :: mode

   ! ---------- child dispatch: one bad case per spawn, must named abort ----------
   allocate (par%e_scale(1))   ! single-carrier check face
   par%e_scale(1) = 1.0d0
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      call fill_list_atoms(g1, mss1)
      call state_create(sta, 3)
      select case (trim(mode))
      case ('bad-escale')                   ! (i): e_scale = 0
         call set_scheme([code_w])
         call build_table(mss1, nmode1, ctab1)
         call table_to_par(1.0d0, w1, ctab1, par)
         par%e_scale(1) = 0.0d0
         call wigner_init(par)
      case ('bad-bind0')                    ! (j): no fragment carries the code
         call set_scheme([1])
         call build_table(mss1, nmode1, ctab1)
         call table_to_par(1.0d0, w1, ctab1, par)
         call wigner_init(par)
      case ('bad-notbl')                    ! (k): unallocated table set
         call set_scheme([code_w])
         par%e_scale(1) = 1.0d0
         call wigner_init(par)
      case ('bad-tblcount')                 ! (l): two carriers, one table
         call fill_two(g1, mss1, g2, mss2)
         call state_destroy(sta)
         call state_create(sta, 5)
         call set_scheme([code_w, code_w])
         call build_table(mss1, nmode1, ctab1)
         call table_to_par(1.0d0, w1, ctab1, par)
         par%e_scale = [ 1.0d0, 1.0d0 ]     ! two rows: the TABLE count is the mismatch
         call wigner_init(par)
      case ('bad-empty')                    ! (m): unallocated eigenvector component
         call set_scheme([code_w])
         allocate (par%tbl(1))
         par%e_scale(1) = 1.0d0
         allocate (par%tbl(1)%w(nmode1))
         par%tbl(1)%w = w1
         call wigner_init(par)
      case ('bad-freq')                     ! (n): a zero frequency entry
         call set_scheme([code_w])
         call build_table(mss1, nmode1, ctab1)
         call table_to_par(1.0d0, w1, ctab1, par)
         par%tbl(1)%w(2) = 0.0d0
         call wigner_init(par)
      case ('bad-orth')                     ! (o): one entry perturbed by 1e-3
         call set_scheme([code_w])
         call build_table(mss1, nmode1, ctab1)
         call table_to_par(1.0d0, w1, ctab1, par)
         par%tbl(1)%c(1, 1) = par%tbl(1)%c(1, 1) + 1.0d-3
         call wigner_init(par)
      case ('bad-noinit')                   ! (p): sample without init
         call wigner_sample(sta)
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (the member returned - bug)'
      stop 0
   end if

   ! ---------- theory constants (the anchor, independent of the member) ----------
   e1 = hbar_evfs*ev_dyn*w1/2.0d0      ! [amu*Ang^2/fs^2] per mode
   e2 = hbar_evfs*ev_dyn*w2/2.0d0
   sig_q2_1 = e1/w1**2
   sig_q2_2 = e2/w2**2
   esum = sum(e1)
   esq_sum = sum(e1**2)
   call build_table(mss1, nmode1, ctab1)
   call build_table(mss2, nmode2, ctab2)

   ! ---------- (a)-(f) single-carrier statistics ----------
   call fill_list_atoms(g1, mss1)
   call set_scheme([code_w])
   call state_create(sta, 3)
   call table_to_par(1.0d0, w1, ctab1, par)
   call rng_init(20260921)
   call wigner_init(par)
   allocate (qm(nmode1), pm(nmode1))
   esum = 0.0d0
   q2 = 0.0d0; p2 = 0.0d0; qp = 0.0d0; tbar = 0.0d0; vbar = 0.0d0
   m1 = 0.0d0; mm2 = 0.0d0; mm3 = 0.0d0; mm4 = 0.0d0
   pmom_max = 0.0d0; dcom_max = 0.0d0; resq_max = 0.0d0; resp_max = 0.0d0
   do draw = 1, ndraw
      call wigner_sample(sta)
      call expand_modes(sta%q, sta%p, list1, mss1, g1, ctab1, qm, pm)
      esum = esum + sum(0.5d0*(pm**2 + (w1*qm)**2))
      q2 = q2 + qm**2
      p2 = p2 + pm**2
      qp = qp + qm*pm
      tbar = tbar + 0.5d0*pm**2
      vbar = vbar + 0.5d0*(w1*qm)**2
      m1 = m1 + qm(1); mm2 = mm2 + qm(1)**2; mm3 = mm3 + qm(1)**3; mm4 = mm4 + qm(1)**4

      ! structural exacts of this draw: net momentum, COM, subspace containment
      pmom = 0.0d0
      dcom = 0.0d0
      do k = 1, 3
         pmom = pmom + sta%p(3*list1(k)-2:3*list1(k))
         dcom = dcom + mss1(k)*sta%q(3*list1(k)-2:3*list1(k))
      end do
      pmom_max = max(pmom_max, norm2(pmom))
      dcom_max = max(dcom_max, norm2(dcom))
      x = mass_weight(sta%q, list1, mss1, g1)
      y = p_over_sqrtm(sta%p, list1, mss1)
      vres = norm2(x - matmul(ctab1, qm))
      if (norm2(x) > 0.0d0) resq_max = max(resq_max, vres/norm2(x))
      vres = norm2(y - matmul(ctab1, pm))
      if (norm2(y) > 0.0d0) resp_max = max(resp_max, vres/norm2(y))
   end do
   esum = esum/dble(ndraw)
   q2 = q2/dble(ndraw); p2 = p2/dble(ndraw); qp = qp/dble(ndraw)
   tbar = tbar/dble(ndraw); vbar = vbar/dble(ndraw)
   m1 = m1/dble(ndraw); mm2 = mm2/dble(ndraw); mm3 = mm3/dble(ndraw); mm4 = mm4/dble(ndraw)

   ! (a) mean total mode energy
   etol = 5.0d0*sqrt(esq_sum/dble(ndraw))
   nfail0 = nfail
   if (abs(esum - sum(e1)) > etol) then
      print '(a,f9.6,a,f9.6,a,es9.2)', 'FAIL (a) mean total energy: sample ', esum, &
         ' vs theory ', sum(e1), ' outside 5 sigma ', etol
      nfail = nfail + 1
   else
      print '(a,f9.6,a,f9.6,a,es9.2)', '(a) mean total mode energy: PASS (sample ', esum, &
         ' vs theory ', sum(e1), ', 5 sigma = ', etol, ')'
   end if

   ! (b) per-mode variances (Q and P)
   do k = 1, nmode1
      tol = 5.0d0*sig_q2_1(k)*sqrt(2.0d0/dble(ndraw))
      if (abs(q2(k) - sig_q2_1(k)) > tol) then
         print '(a,i0,a,f9.6,a,f9.6)', 'FAIL (b) Q variance mode ', k, ': sample ', q2(k), &
            ' vs theory ', sig_q2_1(k)
         nfail = nfail + 1
      end if
      tol = 5.0d0*e1(k)*sqrt(2.0d0/dble(ndraw))
      if (abs(p2(k) - e1(k)) > tol) then
         print '(a,i0,a,f9.6,a,f9.6)', 'FAIL (b) P variance mode ', k, ': sample ', p2(k), &
            ' vs theory ', e1(k)
         nfail = nfail + 1
      end if
   end do
   if (nfail == nfail0) print '(a)', '(b) per-mode variances: PASS (<Q_k^2> = E_k/w_k^2 and <P_k^2> = E_k, all 6 modes)'

   ! (c) q-p decorrelation
   nfail0 = nfail
   do k = 1, nmode1
      corr = qp(k)/sqrt(q2(k)*p2(k))
      if (abs(corr) > 5.0d0/sqrt(dble(ndraw))) then
         print '(a,i0,a,f9.6)', 'FAIL (c) corr(Q,P) mode ', k, ': ', corr
         nfail = nfail + 1
      end if
   end do
   if (nfail == nfail0) print '(a)', '(c) q-p decorrelation: PASS (|corr| <= 0.0112, all 6 modes)'

   ! (d) virial in the mean
   nfail0 = nfail
   do k = 1, nmode1
      tol = 5.0d0*e1(k)/sqrt(dble(ndraw))
      if (abs(tbar(k) - vbar(k)) > tol) then
         print '(a,i0,a,f10.6,a,f10.6)', 'FAIL (d) virial mode ', k, ': <T> = ', tbar(k), &
            ' vs <V> = ', vbar(k)
         nfail = nfail + 1
      end if
   end do
   if (nfail == nfail0) print '(a)', '(d) virial in the mean: PASS (<T_k> = <V_k> = E_k/2, all 6 modes)'

   ! (e) Gaussian shape on mode 1 (central moments from raw moments)
   mu2 = mm2 - m1**2
   mu3 = mm3 - 3.0d0*m1*mu2 - m1**3
   mu4 = mm4 - 4.0d0*m1*mm3 + 6.0d0*m1**2*mu2 + 3.0d0*m1**4
   skw = mu3/max(mu2, 1.0d-300)**1.5d0
   krt = mu4/max(mu2, 1.0d-300)**2
   nfail0 = nfail
   if (abs(skw) > 5.0d0*sqrt(6.0d0/dble(ndraw))) then
      print '(a,f9.6,a)', 'FAIL (e) skewness: ', skw, ' outside 5 sigma 0.0274'
      nfail = nfail + 1
   else if (abs(krt - 3.0d0) > 5.0d0*sqrt(24.0d0/dble(ndraw))) then
      print '(a,f9.6,a)', 'FAIL (e) kurtosis: ', krt, ' outside 3 +- 0.0548'
      nfail = nfail + 1
   else
      print '(a,f8.4,a,f8.4,a)', '(e) Gaussian shape mode 1: PASS (skewness ', skw, &
         ', kurtosis ', krt, '; N(0,sigma^2) shape)'
   end if

   ! (f) structural exacts
   nfail0 = nfail
   if (pmom_max > 1.0d-9*sqrt(2.0d0*sum(mss1)*sum(e1))) then
      print '(a,es9.2)', 'FAIL (f) net momentum: max ', pmom_max
      nfail = nfail + 1
   else if (dcom_max > 1.0d-10) then
      print '(a,es9.2)', 'FAIL (f) COM drift: max ', dcom_max
      nfail = nfail + 1
   else if (resq_max > 1.0d-10 .or. resp_max > 1.0d-10) then
      print '(a,2es9.2)', 'FAIL (f) subspace containment: residuals ', resq_max, resp_max
      nfail = nfail + 1
   else
      print '(a,es9.2,a,es9.2)', '(f) structural exacts: PASS (max |sum p| = ', pmom_max, &
         ', max COM drift = ', dcom_max, ', containment <= 1e-10 every draw)'
   end if
   call state_destroy(sta)

   ! ---------- (g) e_scale plumbing: 0.5 halves every width ----------
   call fill_list_atoms(g1, mss1)
   call set_scheme([code_w])
   call state_create(sta, 3)
   call table_to_par(0.5d0, w1, ctab1, par)
   call wigner_init(par)
   q2 = 0.0d0
   do draw = 1, ndraw
      call wigner_sample(sta)
      call expand_modes(sta%q, sta%p, list1, mss1, g1, ctab1, qm, pm)
      q2 = q2 + qm**2
   end do
   q2 = q2/dble(ndraw)
   nfail0 = nfail
   do k = 1, nmode1
      tol = 5.0d0*(0.5d0*sig_q2_1(k))*sqrt(2.0d0/dble(ndraw))
      if (abs(q2(k) - 0.5d0*sig_q2_1(k)) > tol) then
         print '(a,i0,a,f9.6,a,f9.6)', 'FAIL (g) e_scale variance mode ', k, ': sample ', q2(k), &
            ' vs theory ', 0.5d0*sig_q2_1(k)
         nfail = nfail + 1
      end if
   end do
   if (nfail == nfail0) print '(a)', '(g) e_scale plumbing: PASS (e_scale = 0.5 halves <Q^2> on every mode)'
   call state_destroy(sta)

   ! ---------- (h) dual-carrier equal treatment ----------
   call fill_two(g1, mss1, g2, mss2)
   call set_scheme([code_w, code_w])
   call state_create(sta, 5)
   if (allocated(par%tbl)) deallocate (par%tbl)
   allocate (par%tbl(2))
   par%e_scale = [ 1.0d0, 1.0d0 ]        ! one scale, two rows
   allocate (par%tbl(1)%w(nmode1), par%tbl(1)%c(9, nmode1))
   par%tbl(1)%w = w1
   par%tbl(1)%c = ctab1
   allocate (par%tbl(2)%w(nmode2), par%tbl(2)%c(6, nmode2))
   par%tbl(2)%w = w2
   par%tbl(2)%c = ctab2
   call wigner_init(par)
   q2 = 0.0d0
   do draw = 1, ndraw
      call wigner_sample(sta)
      call expand_modes(sta%q, sta%p, list1, mss1, g1, ctab1, qm, pm)
      q2(1) = q2(1) + qm(1)**2
      call expand_modes(sta%q, sta%p, list2, mss2, g2, ctab2, qm, pm)
      q2(2) = q2(2) + qm(1)**2
   end do
   q2(1:2) = q2(1:2)/dble(ndraw)
   nfail0 = nfail
   if (abs(q2(1) - sig_q2_1(1)) > 5.0d0*sig_q2_1(1)*sqrt(2.0d0/dble(ndraw)) .or. &
       abs(q2(2) - sig_q2_2(1)) > 5.0d0*sig_q2_2(1)*sqrt(2.0d0/dble(ndraw))) then
      print '(a,2f9.6,a,2f9.6)', 'FAIL (h) dual-carrier variances: samples ', q2(1:2), &
         ' vs theory ', [ sig_q2_1(1), sig_q2_2(1) ]
      nfail = nfail + 1
   else
      print '(a,2f9.6,a)', '(h) dual-carrier equal treatment: PASS (mode-1 <Q^2> = ', q2(1:2), &
         ' - each carrier follows its own table at the one e_scale)'
   end if
   call state_destroy(sta)

   ! ---------- (i)-(p) named-abort criteria (child processes) ----------
   call expect_loudstop('bad-escale',  'wigner_init',   'e_scale',     '(i) non-positive e_scale loud-stop')
   call expect_loudstop('bad-bind0',   'wigner_init',   'wigner',      '(j) no-fragment-binding loud-stop')
   call expect_loudstop('bad-notbl',   'wigner_init',   'mismatch',    '(k) unallocated-table-set loud-stop')
   call expect_loudstop('bad-tblcount','wigner_init',   'carriers',    '(l) table-count/carrier-count mismatch loud-stop')
   call expect_loudstop('bad-empty',   'wigner_init',   'unallocated', '(m) inner-table unallocated loud-stop')
   call expect_loudstop('bad-freq',    'wigner_init',   'frequen',     '(n) non-positive frequency loud-stop')
   call expect_loudstop('bad-orth',    'wigner_init',   'orthonormal', '(o) non-orthonormal table loud-stop')
   call expect_loudstop('bad-noinit',  'wigner_draw', 'wigner_init', '(p) sample-before-init loud-stop')

   print '(a,i0,a)', '=== check_samp_wigner: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! one-carrier parameter set from the shared table (scale kept explicit)
   subroutine table_to_par(scale, wtab, ctab, p)
      real(8), intent(in) :: scale
      real(8), intent(in) :: wtab(:), ctab(:,:)
      type(wigner_params_t), intent(out) :: p
      allocate (p%e_scale(1))
      p%e_scale(1) = scale
      allocate (p%tbl(1))
      allocate (p%tbl(1)%w(size(wtab)), p%tbl(1)%c(size(ctab, 1), size(ctab, 2)))
      p%tbl(1)%w = wtab
      p%tbl(1)%c = ctab
   end subroutine table_to_par

   ! SYNTHETIC orthonormal mass-weighted columns spanning the complement of
   ! the three translations (Gram-Schmidt from the canonical basis)
   subroutine build_table(mss, ncol, ctab)
      real(8), intent(in) :: mss(:)
      integer, intent(in) :: ncol
      real(8), allocatable, intent(out) :: ctab(:,:)
      integer :: nd, i, a, kk, nkeep
      real(8) :: t(3*size(mss), 3), v(3*size(mss)), nrm
      nd = 3*size(mss)
      t = 0.0d0
      do i = 1, 3
         do a = 1, size(mss)
            t(3*(a-1)+i, i) = sqrt(mss(a))
         end do
         t(:, i) = t(:, i)/norm2(t(:, i))
      end do
      allocate (ctab(nd, ncol))
      nkeep = 0
      do kk = 1, nd
         v = 0.0d0
         v(kk) = 1.0d0
         do i = 1, 3
            v = v - dot_product(t(:, i), v)*t(:, i)
         end do
         do i = 1, nkeep
            v = v - dot_product(ctab(:, i), v)*ctab(:, i)
         end do
         nrm = norm2(v)
         if (nrm > 0.5d0) then
            nkeep = nkeep + 1
            ctab(:, nkeep) = v/nrm
            if (nkeep == ncol) exit
         end if
      end do
   end subroutine build_table

   ! mode expansion of the written state about the member's own seeding
   ! reference (the COM-shifted seed): Q = C^T sqrt(m)(q - (qz - com)),
   ! P = C^T p/sqrt(m) (the sample is the displacement, not the seed)
   subroutine expand_modes(q, p, list, mss, qz, ctab, qm, pm)
      real(8), intent(in) :: q(:), p(:), mss(:), qz(:), ctab(:,:)
      integer, intent(in) :: list(:)
      real(8), intent(out) :: qm(:), pm(:)
      integer :: a, k
      real(8) :: com(3), x(3*size(mss)), y(3*size(mss))
      com = 0.0d0
      do a = 1, size(mss)
         com = com + mss(a)*qz(3*a-2:3*a)
      end do
      com = com/sum(mss)
      do a = 1, size(mss)
         x(3*a-2:3*a) = sqrt(mss(a))*(q(3*list(a)-2:3*list(a)) - (qz(3*a-2:3*a) - com))
         y(3*a-2:3*a) = p(3*list(a)-2:3*list(a))/sqrt(mss(a))
      end do
      do k = 1, size(qm)
         qm(k) = dot_product(ctab(:, k), x)
         pm(k) = dot_product(ctab(:, k), y)
      end do
   end subroutine expand_modes

   ! x = sqrt(m)*(q - (qz - com)) over one fragment's atoms (containment probe)
   function mass_weight(q, list, mss, qz) result(x)
      real(8), intent(in) :: q(:), mss(:), qz(:)
      integer, intent(in) :: list(:)
      real(8) :: x(3*size(mss)), com(3)
      integer :: a
      com = 0.0d0
      do a = 1, size(mss)
         com = com + mss(a)*qz(3*a-2:3*a)
      end do
      com = com/sum(mss)
      do a = 1, size(mss)
         x(3*a-2:3*a) = sqrt(mss(a))*(q(3*list(a)-2:3*list(a)) - (qz(3*a-2:3*a) - com))
      end do
   end function mass_weight

   ! y = p/sqrt(m) over one fragment's atoms (containment probe)
   function p_over_sqrtm(p, list, mss) result(y)
      real(8), intent(in) :: p(:), mss(:)
      integer, intent(in) :: list(:)
      real(8) :: y(3*size(mss))
      integer :: a
      do a = 1, size(mss)
         y(3*a-2:3*a) = p(3*list(a)-2:3*list(a))/sqrt(mss(a))
      end do
   end function p_over_sqrtm

   ! single-fragment list_atoms by direct field fill
   subroutine fill_list_atoms(g, mss)
      real(8), intent(in) :: g(:), mss(:)
      integer :: nat
      nat = size(mss)
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

   ! two-fragment list_atoms by direct field fill (atoms 1..3 and 4..5)
   subroutine fill_two(ga, mssa, gb, mssb)
      real(8), intent(in) :: ga(:), mssa(:), gb(:), mssb(:)
      call list_atoms_clear()
      allocate (list_atoms%mass(size(mssa) + size(mssb)))
      list_atoms%mass(1:size(mssa)) = mssa
      list_atoms%mass(size(mssa)+1:) = mssb
      allocate (list_atoms%frag(2))
      list_atoms%frag(1)%nat = size(mssa)
      list_atoms%frag(1)%mass = sum(mssa)
      allocate (list_atoms%frag(1)%list(size(mssa)), list_atoms%frag(1)%qz(3*size(mssa)))
      list_atoms%frag(1)%list = [ (k, k = 1, size(mssa)) ]
      list_atoms%frag(1)%qz = ga
      list_atoms%frag(2)%nat = size(mssb)
      list_atoms%frag(2)%mass = sum(mssb)
      allocate (list_atoms%frag(2)%list(size(mssb)), list_atoms%frag(2)%qz(3*size(mssb)))
      list_atoms%frag(2)%list = [ (k, k = size(mssa)+1, size(mssa)+size(mssb)) ]
      list_atoms%frag(2)%qz = gb
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

   ! spawn the child for one bad case and assert the named abort: nonzero exit,
   ! stderr names the authority and carries the expected key
   subroutine expect_loudstop(cmode, who, key, label)
      character(len=*), intent(in) :: cmode, who, key, label
      character(len=512) :: cmd, line
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      ! the shell redirect needs native separators on Windows
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_samp_wigner.exe '//trim(cmode)//' 2>build'//sep//'wigner_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/wigner_err.txt', status='old', action='read', iostat=ios)
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

   ! Windows-host probe (same runtime OS branch the sibling checks use)
   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

end program check_samp_wigner
