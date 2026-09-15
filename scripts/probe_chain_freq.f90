!=====================================================================
! probe_chain_freq.f90 - chain-level frequency-truth probe: real
!                         container forces -> central-difference Hessian
!                         -> eig_sym -> spectrum -> the boltzmann /
!                         normalmode members sampling on the DERIVED
!                         (real-frequency) tables (criteria defined FIRST
!                         in this header; this program executes them).
!                         Companion run material of
!                         cases/case_chain_frequencies (the case README
!                         documents the anchors and the wiring finding).
! Compile (from src/, against the built gfortran object set - driver.o
! excluded, the driver is the only other program unit; objects live in
! per-source subfolders of build/gfortran, .mod files at its root):
!   gfortran -O2 -Ibuild/gfortran -o build/gfortran/probe_chain_freq.exe \
!     ../scripts/probe_chain_freq.f90 \
!     $(find build/gfortran -name '*.o' ! -name 'driver.o')
! Run (cwd = the leg's run directory - reads input_qct.txt, system/,
! reference_pars.txt from cwd, exactly like the driver; the leg word
! selects which member's chain is demonstrated):
!   cd cases/case_chain_frequencies/boltz && \
!     ../../../src/build/gfortran/probe_chain_freq.exe boltz
!   cd cases/case_chain_frequencies/normalmode && \
!     ../../../src/build/gfortran/probe_chain_freq.exe normalmode
! Approach: the REAL assembly chain first (the driver's assembly steps
! 1-8: read_input -> rng -> sysdef -> roster -> state -> registry ->
! seam -> container_init), so the container is live and force_eval
! delivers the true composed forces of the bound reference well. The
! frequency derivation then runs driver-form (the assembly-time
! auto-export channel is mechanism-without-wiring - the case README
! registers the finding): a force wrapper conforming to hessian's
! force_proc feeds force_eval (the container chain) into hessian_cd at
! the roster run geometry; the mass-weighted Hessian goes through
! linalg's eig_sym; the derived spectrum feeds the members' PUBLIC
! faces (params -> member init -> samp_run dispatch -> member sample).
! Expected values are computed in this program from the analytic
! formulas of the constructed harmonic system (theory is the anchor,
! not the member output). The probe writes NO files - every evidence
! line lands on stdout. Exit 0 = all anchors green; exit 1 = any
! failure.
! The harmonic analytic contracts (hand side of every anchor):
!   well            V = K_HARM/2*sum|q_i|^2, K_HARM [eV/A^2] (the
!                   container parameter file, read from cwd); F =
!                   -grad V through force_eval lands in internal units
!                   per A with the single legislated conversion pair
!                   ev_kcal*e_conv = 23.0605*0.04184 = 0.96485132
!                   internal per eV (program constants, restated here
!                   as the analytic formula's own literals)
!   Cartesian H     K_int = K_HARM*0.96485132 [internal/A^2] on every
!                   diagonal; coordinate i's force never reads
!                   coordinate j, so every off-diagonal is the
!                   difference of two bit-identical force components =
!                   exactly 0.0
!   mass-weighting  hmw(i,j) = H(i,j)/sqrt(m_i*m_j) (mass replicated
!                   per component); eigenvalues lambda in 1/(10 fs)^2,
!                   omega = sqrt(lambda) in rad/(10 fs)
!   spectrum        omega_a = sqrt(K_int/m_a), 3-fold per atom; the
!                   wavenumber table handed to boltzmann is the FULL
!                   3N spectrum (the well pins every Cartesian dof -
!                   nothing is free), nu = omega/(2*pi*c) with
!                   c = 2.99792458e-4 cm per 10 fs (nu = omega*530.884)
!   projection      the normalmode member's contract demands columns
!                   orthogonal to the three mass-weighted translations
!                   u_d (components sqrt(m_a)/sqrt(sum m) along axis d),
!                   so the eig_sym input there is the Hessian projected
!                   onto the translation complement,
!                   hp = P*hmw*P with P = I - sum_d u_d u_d^T (3
!                   near-null eigenvalues + 3N-3 modes); with EQUAL
!                   masses the projected spectrum stays analytically
!                   uniform omega = sqrt(K_int/m) (every vector is an
!                   eigenvector of a multiple of the identity)
!   member geometry the normalmode member seeds about the COM origin:
!                   qz_com = qz - COM(qz) (qz = the roster fragment
!                   seed), x = sqrt(m)*(q - qz_com), y = p/sqrt(m),
!                   Q = C^T x, P = C^T y, mode energy
!                   E_i = 0.5*(P_i^2 + w_i^2*Q_i^2)
! Criteria (boltz leg, HCl carrier - two distinct masses sharpen the
! spectrum anchor; N = 200000 draws, one seeded draw loop):
!   (F1) Hessian structure: every diagonal within 1e-10 RELATIVE of
!        K_int; every off-diagonal exactly 0.0 (bitwise)
!   (F2) full spectrum: each of the 3N sorted omega within 1e-10
!        relative of the sorted analytic set {sqrt(K_int/m_a)}
!        (HCl: 3x omega_Cl = 0.032996 + 3x omega_H = 0.195673,
!        nu = 17.517 / 103.879 cm^-1 at K_HARM = 0.04)
!   (B1) geometric chi-square on the highest-nu mode (mode 6, an H
!        mode): buckets n = 0..4 plus the >=5 tail lumped (6 buckets),
!        theory q = exp(-1.43878*nu_analytic/T) = exp(-1.4371) =
!        0.23762 at T = 104.0 K, chi-square df = 5, critical 11.07
!        (alpha = 0.05); PASS requires chi2 <= 11.07
!   (B2) per-mode means: every mode's sample mean within 5 sigma of
!        its own q_k/(1-q_k) (sigma = sqrt(q/(1-q)^2/N); theory means
!        Cl 3.6465 (5 sigma 0.0460), H 0.31163 (5 sigma 7.15e-3) - the
!        two degenerate groups must not cross-talk)
!   (B3) quantum-number integrity: every drawn n of every mode is a
!        non-negative integer (getter returns integers by type)
! Criteria (normalmode leg, H3 carrier - three EQUAL masses make the
! projected spectrum analytically uniform; s = 3N-3 = 6 modes,
! e_vib = 12.5 kcal/mol -> e_tot = 12.5*e_conv = 0.523 internal):
!   (F3) projected spectrum + table contract: the 3 null eigenvalues
!        <= 1e-10 absolute; each mode omega within 1e-10 relative of
!        sqrt(K_int/m); the eigenvector table orthonormal (C^T*C - I)
!        to 1e-12; every column orthogonal to every mass-weighted
!        translation u_d to 1e-12; the eigen-residual
!        ||hp*c_k - lambda_k*c_k|| <= 1e-10 of the spectral scale
!        lambda = K_int/m (the member's own init battery re-validates
!        the same table at its 1e-8 tolerances - passing normalmode_init
!        is itself chain evidence)
!   (N1) exact identities (2000 draws, the member's measure on the
!        GIVEN derived table): the energy re-expanded from the WRITTEN
!        state (Q = C^T sqrt(m)*(q - qz_com), P = C^T p/sqrt(m), then
!        sum 0.5*(P^2 + w^2*Q^2)) equals e_tot to 1e-10 relative;
!        |sum p| <= 1e-9*sqrt(2*M*e_tot) (M the carrier mass); the
!        residuals of sqrt(m)*(q - qz_com) and of p/sqrt(m) after
!        projection onto span(C) <= 1e-10 of the vector norm; the COM
!        of the written q at the origin to 1e-10 Angstrom
!   (N2) mode-1 energy share ~ Beta(1/2,5/2): equal-CDF 10-bin
!        chi-square with the closed-form quantiles (x = sin^2(theta),
!        CDF kernel G(theta) = 3*theta/4 + sin(2*theta)/2 +
!        sin(4*theta)/16, G(pi/2) = 3*pi/8, edges by 80-step bisection),
!        N = 50000, df = 9, critical 16.92 (alpha = 0.05); mean = 1/6
!        within 5 sigma (sigma = sqrt(5/144/N) = 8.33e-4)
!   (N3) phase equipartition: mean kinetic energy = e_tot/2 within
!        5 sigma, sigma = e_tot*sqrt(3/(64*N)) (Var(KE) = 3*e_tot^2/64)
! Determinism: rng_init(i_seed(1)) reseeds before each statistical
! pass (the boltz loop; the normalmode N1 pass; the N2/N3 pass reseeds
! with i_seed(1)+1) - double runs are byte-identical on stdout.
!=====================================================================
module chain_probe_support
   use state,        only: state_t, state_create
   use force_facade, only: force_eval
   implicit none
   private
   public :: probe_force_setup, probe_force
   type(state_t), save :: sta_h    ! the force chain's scratch state
contains
   !------------------------------------------------------------------
   ! probe_force_setup(natoms) - create the scratch state the wrapper
   !                             rides on (one per probe process)
   subroutine probe_force_setup(natoms)
      integer, intent(in) :: natoms
      call state_create(sta_h, natoms)
   end subroutine probe_force_setup

   !------------------------------------------------------------------
   ! probe_force(q, f) - the force_proc wrapper: the REAL container
   !                     chain (force_eval: request -> container ->
   !                     fold -> the single eV-to-internal conversion)
   subroutine probe_force(q, f)
      real(8), intent(in)  :: q(:)   ! position vector [A]
      real(8), intent(out) :: f(:)   ! forces at q [internal/A]
      real(8) :: e_dummy
      sta_h%q(1:size(q)) = q
      call force_eval(sta_h, e_dummy)
      f = sta_h%f(1:size(f))
   end subroutine probe_force
end module chain_probe_support

program probe_chain_freq
   use consts,        only: e_conv, pi, two_pi
   use rng,           only: rng_init
   use control,       only: i_seed
   use input,         only: read_input, staging_n_rows, staging_key, staging_val, staging_line
   use config,        only: reactants
   use config_roster, only: roster, roster_load
   use state,         only: state_t, state_create, state_destroy
   use sysdef_facade, only: sysdef_load
   use registry_gen,  only: registry_gen_all, registry_gen_init
   use force_facade,  only: container_init
   use sampler,       only: samp_code, samp_run
   use seam_cast,     only: seam_real
   use hessian,       only: hessian_cd
   use linalg,        only: eig_sym
   use reference_box, only: reference_k
   use dist_boltzmann,  only: boltz_params_t, boltz_init, boltz_last_n
   use dist_normalmode, only: normalmode_params_t, normalmode_init
   use chain_probe_support, only: probe_force_setup, probe_force
   implicit none

   ! analytic-side constants (the formulas' own literals; ev_kcal and
   ! e_conv are the program's legislated pair, restated here)
   real(8), parameter :: ev_kcal_lit = 23.0605d0      ! kcal/mol per eV (the facade's literal)
   real(8), parameter :: hc_k_lit    = 1.43878d0      ! h*c/k [cm*K] (the member's literal)
   real(8), parameter :: c_light     = 2.99792458d-4 ! speed of light [cm per 10 fs]
   real(8), parameter :: w_per_rad   = 1.0d0/(two_pi*c_light) ! cm^-1 per rad/(10 fs)

   ! sample sizes (the anchors' statistical bands are quoted for these)
   integer, parameter :: ndraw_boltz = 200000
   integer, parameter :: ndraw_exact = 2000
   integer, parameter :: ndraw_stat  = 50000
   integer, parameter :: nbin = 10
   integer, parameter :: s_mode = 6    ! 3N-3 modes of the H3 carrier

   type(state_t) :: sta
   character(len=32) :: leg
   character(len=512) :: sval
   integer :: nfail = 0, nfail0, i, j, k, d, nat, n, draw, ib, it, nneg
   real(8) :: k_int, t_vib, e_vib, e_tot, psc
   real(8), allocatable :: sqm(:), hess(:,:), hmw(:,:), aa(:,:), lam(:), omega(:), nu(:)
   real(8), allocatable :: lam_ana(:), omega_ana(:)
   real(8) :: worst_rel, worst_abs, chi2, q1, qexp(6)
   integer :: freq(6)

   ! ---- leg selection ----
   if (command_argument_count() /= 1) then
      print '(a)', 'probe_chain_freq: usage: probe_chain_freq.exe <boltz|normalmode>'
      stop 2
   end if
   call get_command_argument(1, leg)
   if (trim(leg) /= 'boltz' .and. trim(leg) /= 'normalmode') then
      print '(a)', 'probe_chain_freq: unknown leg "'//trim(leg)//'" (boltz | normalmode)'
      stop 2
   end if

   ! ---- the REAL assembly chain (the driver's assembly steps 1-8) ----
   call read_input('input_qct.txt')
   call rng_init(i_seed(1))
   call sysdef_load(reactants%system_dir)
   call roster_load()
   nat = size(roster%mass)
   n = 3*nat
   call state_create(sta, nat)
   call registry_gen_all()
   call registry_gen_init()          ! DIST_SCHEME = ROTATION seam (the harmless carrier)
   call container_init()
   call probe_force_setup(nat)
   if (size(roster%frag) /= 1) then
      print '(a,i0)', 'probe_chain_freq: fixture expects exactly one roster fragment, got ', &
         size(roster%frag)
      stop 2
   end if
   if (trim(leg) == 'boltz' .and. nat /= 2) stop 'probe_chain_freq: boltz fixture expects 2 atoms'
   if (trim(leg) == 'normalmode' .and. nat /= 3) stop 'probe_chain_freq: normalmode fixture expects 3 atoms'

   k_int = reference_k()*ev_kcal_lit*e_conv    ! analytic Cartesian Hessian [internal/A^2]
   print '(a,a,a,i0,a)', 'probe: leg = ', trim(leg), '  (roster atoms = ', nat, &
      ')  - real container forces through force_eval'
   print '(a,es12.6,a)', 'probe: K_HARM staged -> analytic K_int = K*23.0605*0.04184 = ', k_int, &
      ' internal/A^2'

   ! ---- frequency derivation: real forces -> Hessian -> mass-weighting ----
   allocate (sqm(n), hess(n,n), hmw(n,n), aa(n,n), lam(n), omega(n), nu(n), lam_ana(n), omega_ana(n))
   call hessian_cd(probe_force, roster%q0, hess)     ! default step (the honest assembly site)

   ! (F1) Hessian structure: diagonal vs K_int, off-diagonal bitwise zero
   worst_rel = 0.0d0
   do i = 1, n
      worst_rel = max(worst_rel, abs(hess(i,i) - k_int)/k_int)
   end do
   worst_abs = 0.0d0
   do j = 1, n
      do i = 1, n
         if (i /= j) worst_abs = max(worst_abs, abs(hess(i,j)))
      end do
   end do
   nfail0 = nfail
   print '(a,es10.3,a,es10.3)', 'probe: (F1) worst diagonal rel = ', worst_rel, &
      '  worst off-diagonal abs = ', worst_abs
   if (worst_rel >= 1.0d-10) then
      print '(a)', 'FAIL (F1): diagonal outside 1e-10 relative of K_int'
      nfail = nfail + 1
   else if (worst_abs /= 0.0d0) then
      print '(a)', 'FAIL (F1): off-diagonal not bitwise zero'
      nfail = nfail + 1
   else
      print '(a)', '(F1) Hessian structure: PASS (diag = K_int within 1e-10 rel, off-diag exactly 0.0;'
      print '(a)', '     central difference at the default step on the container chain)'
   end if

   ! mass-weighted Hessian (mass replicated per component)
   do i = 1, n
      sqm(i) = sqrt(roster%mass((i + 2)/3))
   end do
   hmw = hess/(spread(sqm, 2, n)*spread(sqm, 1, n))

   ! analytic spectrum: 3 copies of K_int/m_a per atom, ascending
   i = 0
   do j = 1, nat
      do d = 1, 3
         i = i + 1
         lam_ana(i) = k_int/roster%mass(j)
      end do
   end do
   call sort_real(lam_ana)
   omega_ana = sqrt(lam_ana)

   if (trim(leg) == 'boltz') then
      ! ============ boltz leg: full spectrum -> geometric sampling ============
      aa = hmw
      call eig_sym(aa, lam)
      omega = sqrt(lam)                            ! ascending by construction
      nu = omega*w_per_rad                         ! the derived wavenumber table [cm^-1]

      ! (F2) per-entry spectrum anchor
      nfail0 = nfail
      worst_rel = 0.0d0
      do k = 1, n
         worst_rel = max(worst_rel, abs(omega(k) - omega_ana(k))/omega_ana(k))
      end do
      print '(a,es10.3)', 'probe: (F2) worst omega rel = ', worst_rel
      do k = 1, n
         print '(a,i2,a,f10.7,a,f10.4,a)', 'probe:   mode ', k, ': omega = ', omega(k), &
            ' rad/(10 fs) -> nu = ', nu(k), ' cm^-1'
      end do
      if (worst_rel >= 1.0d-10) then
         print '(a)', 'FAIL (F2): derived spectrum outside 1e-10 relative of the analytic set'
         nfail = nfail + 1
      else
         print '(a,f10.7,a,f10.7,a)', '(F2) full spectrum: PASS (3x ', omega(n), ' + 3x ', omega(1), &
            ' = analytic {sqrt(K_int/m_a)} within 1e-10 rel)'
      end if

      ! staged T_VIB_A (the seam's own read pattern)
      t_vib = 0.0d0
      do i = 1, staging_n_rows()
         if (trim(staging_key(i)) == 'T_VIB_A') then
            sval = staging_val(i)
            call seam_real(staging_key(i), sval, staging_line(i), t_vib)
         end if
      end do
      if (t_vib <= 0.0d0) stop 'probe_chain_freq: T_VIB_A not staged (boltz leg needs it)'

      ! the member's public face: derived table -> init -> samp_run dispatch
      block
         type(boltz_params_t) :: par
         real(8) :: amean_m(n), emean_m(n), sig_m(n)
         integer :: kk, nb
         par%t_vib_a = t_vib
         allocate (par%w_mode_a(n))
         par%w_mode_a = nu
         par%n_mode_a = n
         deallocate (reactants%dist_scheme)
         allocate (reactants%dist_scheme(1))
         reactants%dist_scheme = samp_code('boltzmann')
         call boltz_init(par)          ! member validation + fragment binding
         call rng_init(i_seed(1))      ! seeded statistical pass

         ! one draw loop feeds (B1)/(B2)/(B3)
         amean_m = 0.0d0
         freq = 0
         nneg = 0
         do draw = 1, ndraw_boltz
            call samp_run('incident', sta)
            do kk = 1, n
               ib = boltz_last_n(kk)
               amean_m(kk) = amean_m(kk) + dble(ib)
               if (ib < 0) nneg = nneg + 1            ! (B3) witness
            end do
            nb = boltz_last_n(n)                      ! the highest-nu mode (an H mode)
            if (nb < 5) then
               freq(nb + 1) = freq(nb + 1) + 1
            else
               freq(6) = freq(6) + 1                  ! >=5 tail lumped into bucket 6
            end if
         end do
         amean_m = amean_m/dble(ndraw_boltz)

         ! (B1) geometric chi-square, theory q from the ANALYTIC nu of mode 6
         q1 = exp(-hc_k_lit*omega_ana(n)*w_per_rad/t_vib)
         qexp(1) = 1.0d0 - q1
         do i = 2, 5
            qexp(i) = (1.0d0 - q1)*q1**dble(i - 1)
         end do
         qexp(6) = q1**5.0d0
         chi2 = 0.0d0
         do i = 1, 6
            chi2 = chi2 + (dble(freq(i)) - ndraw_boltz*qexp(i))**2/(ndraw_boltz*qexp(i))
         end do
         nfail0 = nfail
         print '(a,f9.3,a,f9.6,a)', 'probe: (B1) chi2 = ', chi2, '  (theory q = ', q1, ')'
         if (chi2 > 11.07d0) then
            print '(a)', 'FAIL (B1): chi-square above 11.07 (df=5, alpha=0.05)'
            nfail = nfail + 1
         else
            print '(a,f7.3,a,i0,a)', '(B1) geometric chi-square: PASS (chi2 = ', chi2, &
               ' <= 11.07, df = 5, N = ', ndraw_boltz, ')'
         end if

         ! (B2) per-mode means vs each mode's own q_k/(1-q_k)
         do kk = 1, n
            q1 = exp(-hc_k_lit*omega_ana(kk)*w_per_rad/t_vib)
            emean_m(kk) = q1/(1.0d0 - q1)
            sig_m(kk) = sqrt(q1/(1.0d0 - q1)**2/dble(ndraw_boltz))
         end do
         do kk = 1, n
            if (abs(amean_m(kk) - emean_m(kk)) > 5.0d0*sig_m(kk)) then
               print '(a,i2,a,f9.5,a,f9.5,a,es9.2)', 'FAIL (B2) mode ', kk, ': mean ', amean_m(kk), &
                  ' vs theory ', emean_m(kk), ' outside 5 sigma = ', 5.0d0*sig_m(kk)
               nfail = nfail + 1
            end if
         end do
         if (nfail == nfail0) then
            print '(a,3f8.4,a,3f8.4)', '(B2) per-mode means: PASS (Cl group ', amean_m(1:3), &
               ' / H group ', amean_m(4:6), '; each within its own 5 sigma)'
            print '(a,f9.5,a,f9.5)', 'probe:   theory means: Cl = ', emean_m(1), '  H = ', emean_m(n)
            print '(a,es9.2,a,es9.2)', 'probe:   5 sigma:    Cl = ', 5.0d0*sig_m(1), '  H = ', 5.0d0*sig_m(n)
         end if

         ! (B3) quantum-number integrity
         if (nneg > 0) then
            print '(a,i0,a)', 'FAIL (B3): ', nneg, ' negative quantum number(s) drawn'
            nfail = nfail + 1
         else
            print '(a)', '(B3) quantum-number integrity: PASS (all drawn n >= 0, integers by getter type)'
         end if
      end block

   else
      ! ======== normalmode leg: projected spectrum -> microcanonical ========
      block
         type(normalmode_params_t) :: par
         real(8) :: trans(9, 3), hp(9, 9), c_tab(9, s_mode), w_tab(s_mode)
         real(8) :: qz(9), xv(9), yv(9), qm(s_mode), pm(s_mode)
         real(8) :: com3(3), psum(3), res(9), tu(3, 9), uu(3, 3)
         real(8) :: m_sum, lam_ana_uni, e_sum, xres, yres, rmax, dev
         real(8) :: amean, emean, sig, kmean, e1s, x_edge(nbin - 1), lo, hi, mid, gmax
         integer :: binct(nbin)
         logical :: ok_all

         ! translation projector: orthonormal mass-weighted translations u_d
         m_sum = sum(roster%mass)
         trans = 0.0d0
         do d = 1, 3
            do i = 1, nat
               trans(3*(i - 1) + d, d) = sqrt(roster%mass(i))/sqrt(m_sum)
            end do
         end do
         ! hp = P*hmw*P with P = I - U*U^T, expanded as
         ! hmw - U*(U^T*hmw) - (hmw*U)*U^T + U*(U^T*hmw*U)*U^T (then symmetrized)
         tu = matmul(transpose(trans), hmw)          ! 3x9 = U^T*hmw
         uu = matmul(tu, trans)                      ! 3x3 = U^T*hmw*U
         hp = hmw - matmul(trans, tu) - matmul(transpose(tu), transpose(trans)) &
                + matmul(trans, matmul(uu, transpose(trans)))
         hp = 0.5d0*(hp + transpose(hp))
         aa = hp
         call eig_sym(aa, lam)                       ! lam ascending: 3 nulls + 6 modes
         w_tab = sqrt(lam(4:9))
         c_tab = aa(:, 4:9)                          ! eigenvectors by column

         ! (F3) projected spectrum + table contract
         lam_ana_uni = k_int/roster%mass(1)          ! equal masses -> uniform analytic lambda
         nfail0 = nfail
         worst_abs = maxval(abs(lam(1:3)))
         worst_rel = 0.0d0
         do k = 1, s_mode
            worst_rel = max(worst_rel, abs(w_tab(k)**2 - lam_ana_uni)/lam_ana_uni)
         end do
         dev = 0.0d0
         do k = 1, s_mode
            do i = 1, s_mode
               dev = max(dev, abs(dot_product(c_tab(:, k), c_tab(:, i)) - merge(1.0d0, 0.0d0, k == i)))
            end do
         end do
         rmax = 0.0d0
         do d = 1, 3
            do k = 1, s_mode
               rmax = max(rmax, abs(dot_product(trans(:, d), c_tab(:, k))))
            end do
         end do
         xres = 0.0d0
         do k = 1, s_mode
            res = matmul(hp, c_tab(:, k)) - (w_tab(k)**2)*c_tab(:, k)
            xres = max(xres, norm2(res))
         end do
         print '(a,es10.3,a,es10.3,a)', 'probe: (F3) null |lambda| max = ', worst_abs, &
            '  mode lambda rel = ', worst_rel, ' (vs K_int/m)'
         print '(a,es10.3,a,es10.3,a,es10.3)', 'probe: (F3) orthonormality dev = ', dev, &
            '  translation overlap = ', rmax, '  eigen-residual = ', xres
         do k = 1, s_mode
            print '(a,i2,a,f10.7,a)', 'probe:   mode ', k, ': omega = ', w_tab(k), &
               ' rad/(10 fs) (projected spectrum)'
         end do
         if (worst_abs > 1.0d-10 .or. worst_rel >= 1.0d-10 .or. dev > 1.0d-12 .or. &
             rmax > 1.0d-12 .or. xres > 1.0d-10*lam_ana_uni) then
            print '(a)', 'FAIL (F3): projected spectrum or table contract out of tolerance'
            nfail = nfail + 1
         else
            print '(a,f10.7,a)', '(F3) projected spectrum + table: PASS (3 nulls <= 1e-10, 6 modes = ', &
               sqrt(lam_ana_uni), ' = sqrt(K_int/m) within 1e-10 rel; C orthonormal 1e-12)'
            print '(a)', '     translation-orthogonal 1e-12, eigen-residual within 1e-10 of the scale'
         end if

         ! staged E_VIB (the seam's own read pattern)
         e_vib = 0.0d0
         do i = 1, staging_n_rows()
            if (trim(staging_key(i)) == 'E_VIB') then
               sval = staging_val(i)
               call seam_real(staging_key(i), sval, staging_line(i), e_vib)
            end if
         end do
         if (e_vib <= 0.0d0) stop 'probe_chain_freq: E_VIB not staged (normalmode leg needs it)'
         e_tot = e_vib*e_conv
         psc = sqrt(2.0d0*m_sum*e_tot)

         ! the member's public face: derived table -> init (its own 1e-8 battery
         ! re-validates the DERIVED table) -> samp_run dispatch
         par%e_vib = e_vib
         allocate (par%w_mode(s_mode), par%c_mode(9, s_mode))
         par%w_mode = w_tab
         par%c_mode = c_tab
         deallocate (reactants%dist_scheme)
         allocate (reactants%dist_scheme(1))
         reactants%dist_scheme = samp_code('normalmode')
         call normalmode_init(par)
         print '(a)', 'probe:   normalmode_init PASSED on the derived table (member 1e-8 battery green)'

         ! (N1) exact identities (2000 draws)
         qz = roster%frag(1)%qz
         call rng_init(i_seed(1))
         ok_all = .true.
         do draw = 1, ndraw_exact
            call samp_run('incident', sta)
            call reexpand(sta, qz, c_tab, xv, yv, qm, pm, com3)
            e_sum = 0.0d0
            do k = 1, s_mode
               e_sum = e_sum + 0.5d0*(pm(k)**2 + w_tab(k)**2*qm(k)**2)
            end do
            if (abs(e_sum - e_tot) > 1.0d-10*e_tot) ok_all = .false.
            psum = sta%p(1:3) + sta%p(4:6) + sta%p(7:9)
            if (norm2(psum) > 1.0d-9*psc) ok_all = .false.
            call project_out(xv, c_tab, res)
            xres = norm2(res)
            call project_out(yv, c_tab, res)
            yres = norm2(res)
            if (xres > 1.0d-10*max(norm2(xv), 1.0d-12)) ok_all = .false.
            if (yres > 1.0d-10*max(norm2(yv), 1.0d-12)) ok_all = .false.
            if (norm2(com3) > 1.0d-10) ok_all = .false.
         end do
         if (ok_all) then
            print '(a,es12.4,a)', '(N1) exact identities: PASS (E = e_vib*e_conv = ', e_tot, &
               ' internal to 1e-10; |sum p| <= 1e-9*sqrt(2*M*E); state in span(C); COM at origin;'
            print '(a,i0,a)', '     N = ', ndraw_exact, ' draws)'
         else
            print '(a)', 'FAIL (N1): energy identity, momentum, COM or subspace containment broke'
            nfail = nfail + 1
         end if

         ! (N2)/(N3) statistical pass (fresh seeding)
         ! Beta(1/2,5/2) equal-CDF bin edges by bisection on the closed form
         gmax = 0.375d0*pi                    ! G(pi/2) = 3*pi/8
         do k = 1, nbin - 1
            lo = 0.0d0
            hi = 0.5d0*pi
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
         call rng_init(i_seed(1) + 1)
         binct = 0
         amean = 0.0d0
         kmean = 0.0d0
         do draw = 1, ndraw_stat
            call samp_run('incident', sta)
            call reexpand(sta, qz, c_tab, xv, yv, qm, pm, com3)
            e1s = 0.5d0*(pm(1)**2 + w_tab(1)**2*qm(1)**2)/e_tot
            amean = amean + e1s
            ib = 1
            do k = 1, nbin - 1
               if (e1s > x_edge(k)) ib = ib + 1
            end do
            binct(ib) = binct(ib) + 1
            do k = 1, s_mode
               kmean = kmean + 0.5d0*pm(k)**2
            end do
         end do
         amean = amean/dble(ndraw_stat)
         kmean = kmean/dble(ndraw_stat)
         chi2 = 0.0d0
         do ib = 1, nbin
            chi2 = chi2 + (dble(binct(ib)) - ndraw_stat/dble(nbin))**2/(ndraw_stat/dble(nbin))
         end do
         nfail0 = nfail
         emean = 1.0d0/6.0d0
         sig = sqrt((5.0d0/144.0d0)/dble(ndraw_stat))
         print '(a,f9.3,a,f9.6,a)', 'probe: (N2) chi2 = ', chi2, '  mean share = ', amean, &
            '  (theory mean 1/6)'
         if (chi2 > 16.92d0) then
            print '(a)', 'FAIL (N2): Beta chi-square above 16.92 (df=9, alpha=0.05)'
            nfail = nfail + 1
         else if (abs(amean - emean) > 5.0d0*sig) then
            print '(a,f9.6,a,f9.6)', 'FAIL (N2): Beta mean ', amean, ' vs 1/6 = ', emean
            nfail = nfail + 1
         else
            print '(a,f7.3,a,f8.5,a,i0,a)', '(N2) mode-1 energy share: PASS (chi2 = ', chi2, &
               ' <= 16.92, df = 9; mean ', amean, ' vs Beta(1/2,5/2) mean 1/6; N = ', ndraw_stat, ')'
         end if
         sig = e_tot*sqrt(3.0d0/(64.0d0*dble(ndraw_stat)))
         print '(a,es11.4,a,es11.4)', 'probe: (N3) mean KE = ', kmean, '  vs e_tot/2 = ', 0.5d0*e_tot
         if (abs(kmean - 0.5d0*e_tot) > 5.0d0*sig) then
            print '(a)', 'FAIL (N3): equipartition mean KE outside 5 sigma of e_tot/2'
            nfail = nfail + 1
         else
            print '(a,es11.4,a,i0,a)', '(N3) phase equipartition: PASS (mean KE = ', kmean, &
               ' = e_tot/2 within 5 sigma; N = ', ndraw_stat, ')'
         end if
      end block
   end if

   print '(a,i0,a)', '=== probe_chain_freq: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! closed-form kernel of the Beta(1/2,5/2) CDF (monotone on [0,pi/2])
   real(8) function g_theta(th)
      real(8), intent(in) :: th
      g_theta = 0.75d0*th + 0.5d0*sin(2.0d0*th) + 0.0625d0*sin(4.0d0*th)
   end function g_theta

   ! re-expand the WRITTEN state into mass-weighted mode coordinates (the
   ! member's convention restated: seed reference = COM-shifted qz)
   subroutine reexpand(st, gz, ctab, xv, yv, qm, pm, cm)
      type(state_t), intent(in) :: st
      real(8), intent(in)  :: gz(9), ctab(9, s_mode)
      real(8), intent(out) :: xv(9), yv(9), qm(s_mode), pm(s_mode), cm(3)
      real(8) :: qzc(9), cgv(3)
      integer :: a
      cgv = 0.0d0
      do a = 1, 3
         cgv = cgv + roster%mass(a)*gz(3*a-2:3*a)
      end do
      cgv = cgv/sum(roster%mass)
      cm = 0.0d0
      do a = 1, 3
         cm = cm + roster%mass(a)*st%q(3*a-2:3*a)
      end do
      cm = cm/sum(roster%mass)
      do a = 1, 3
         qzc(3*a-2:3*a) = gz(3*a-2:3*a) - cgv
         xv(3*a-2:3*a) = sqrt(roster%mass(a))*(st%q(3*a-2:3*a) - qzc(3*a-2:3*a))
         yv(3*a-2:3*a) = st%p(3*a-2:3*a)/sqrt(roster%mass(a))
      end do
      qm = matmul(transpose(ctab), xv)
      pm = matmul(transpose(ctab), yv)
   end subroutine reexpand

   ! span(C) projection residual: v - C*(C^T v)
   subroutine project_out(v, ctab, r)
      real(8), intent(in)  :: v(9), ctab(9, s_mode)
      real(8), intent(out) :: r(9)
      r = v - matmul(ctab, matmul(transpose(ctab), v))
   end subroutine project_out

   ! ascending insertion sort (n <= 9 - clarity over speed)
   subroutine sort_real(v)
      real(8), intent(inout) :: v(:)
      real(8) :: tmp
      integer :: a, b
      do a = 2, size(v)
         tmp = v(a)
         b = a - 1
         do while (b >= 1 .and. v(b) > tmp)
            v(b + 1) = v(b)
            b = b - 1
         end do
         v(b + 1) = tmp
      end do
   end subroutine sort_real

end program probe_chain_freq
