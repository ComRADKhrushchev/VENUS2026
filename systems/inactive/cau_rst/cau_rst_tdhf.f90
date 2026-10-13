!=====================================================================
! cau_rst_tdhf.f90 - the RST-Hubbard container's composed force entry +
!   the TDHF electronic method package members (per-container copy of
!   the cau_emt2 template per the 2026-09-21 ruling: the TDHF method
!   stays in the container folder, one copy per C@Au container)
! Design:
!   Two-site Anderson/Hubbard impurity (N_sptl = 2 spatial orbitals, both
!   fed from the scalar EMT density-polynomial delivery) coupled to a
!   band-width bath discretized by Gauss-Legendre quadrature (N_orb/2
!   nodes -> 2*(N_orb/2) bath levels, symmetric +/- band_width/2, Fermi
!   level at band center). Ported from C@Au_TDHF/2D TDHF_U.f90 (the
!   2025.12.31 mean-field scheme; the old file's NAC/FSSH pieces are not
!   used by it and are not ported).
!   Electronic phase: exp(i*E_int*dt/hbar_evfs) with E_int [eV] and dt
!   [fs] (the 2026-10-06 eV-A-amu-fs base, consts.f90) - numerically the
!   same phase as the old chain E[eV]*23.0605*e_conv*dt/hbar_code.
!   cau_rst_pes (the force-slot contract): V = E0 + E_elec + E_mf -
!   energy_metal_zero [eV], the matched kernel with the rho-weighted
!   Hcore-kernel trace + the U/J mean-field chain [eV/A]. Unarmed (the
!   adiabatic degenerate point) it reduces to the pure EMT ground PES.
!   Per-trajectory rho lifecycle: the driver offers no trajectory-head
!   signal to containers (a real contract gap this prototype surfaces),
!   so the seed - a q-independent diabatic occupation - is re-applied at
!   (1) the container's own true termination judgment (the common end)
!   and (2) the prop member's time-backwards guard (t reset = new
!   trajectory). The first force of a max_steps-exhausted trajectory
!   still sees the previous rho - a documented residual, verified only
!   to be harmless by check_cau_rst criterion (i).
!=====================================================================
module cau_rst_tdhf
   use consts,       only: hbar_evfs, kb_ev, pi
   use cau_rst_pot
   use linalg,       only: eig_sym
   use rng,          only: rng_u
   use state,        only: state_t, elec_t
   implicit none
   private
   public :: cau_rst_tdhf_setup, cau_rst_tdhf_armed, cau_rst_tdhf_end_traj
   public :: cau_rst_diag_out
   public :: bath_out, rho_seed_out, rho_herm_out, hcore_out          ! criterion introspection
   public :: kernel_caches_out, hub_vals_out                          ! (check_cau_rst legs)
   public :: cau_rst_pes                        ! the container's composed force entry
   public :: cau_rst_prop, cau_rst_match        ! the 'tdhf_hub' method-package members
   public :: cau_rst_sync                        ! the package's container-publish member
   public :: cau_rst_col_nup, cau_rst_col_ndn, cau_rst_col_mag, &  ! recorder columns
             cau_rst_col_eelec, cau_rst_col_emf                    ! (obs_list-selected)
   public :: rho_ref_out                          ! rho reference copy (criterion leg)

   ! -- conversion constants (the legislated pair; delivery stays in eV) --
   ! (ev_kcal/e_conv pair retired with the 2026-10-06 eV base: phases ride
   !  hbar_evfs [eV*fs], bath kT rides kb_ev [eV/K])

   ! -- bath / model dimensions --
   integer :: n_orb = 20                         ! spinless bath orbital count (even)
   integer :: n_sptl = 2                         ! impurity spatial orbitals (fixed by the model)
   integer :: n_siz = 0                          ! total single-particle dimension n_orb + n_sptl
   real(8) :: band_width = 7.0d0                 ! metal bath width [eV]
   real(8) :: t_elec = 300.0d0                   ! electronic temperature [K] (bath Fermi-Dirac)
   real(8) :: mu_ev = 0.0d0                      ! chemical potential [eV] (band center - set_occ(.,0))
   character(len=2) :: start_state = 'GS'        ! initial impurity occupation ('GS' | '1D')
   integer :: msbstp = 4                         ! electronic substeps per prop call
   integer :: perturb_mode = 0                   ! h degeneracy break: 0=off 1=per Hcore build 2=once
   real(8) :: perturb_amp = 0.05d0               ! break half-width [eV]
   logical :: perturb_done = .false.             ! mode-2 once flag

   ! -- quadrature + bath --
   real(8), allocatable :: e_levels(:), w_levels(:)  ! Gauss-Legendre nodes/weights on [-1,1]
   real(8), allocatable :: e_metal(:)                ! bath orbital energies [eV], symmetric band
   real(8) :: energy_metal_zero = 0.0d0              ! bath zero-energy reference [eV]

   ! -- Hamiltonian / electronic state (module-private, run-lifetime) --
   real(8), allocatable :: hcore(:, :, :), fock(:, :, :)   ! (n_siz, n_siz, 2)
   real(8), allocatable :: e_hf(:, :), c_hf(:, :, :)       ! Fock eigenpairs (per spin)
   real(8) :: h_vec(2) = 0.0d0, vg_vec(2) = 0.0d0          ! impurity level / coupling [eV]
   real(8) :: u_val = 0.0d0, j_val = 0.0d0                 ! Hubbard U / exchange J [eV]
   complex(8), allocatable :: rhom(:, :, :)                ! density matrix, diabatic rep (n,n,2)
   logical :: armed = .false.                              ! package-arm flag (the seam sets it)
   real(8) :: t_seen = -1.0d30                             ! last prop-member trajectory time [fs]

   ! -- cached delivery (one ground+Hubbard evaluation per geometry) --
   real(8) :: e0_cache = 0.0d0
   real(8), allocatable :: g0_cache(:), du_cache(:), dj_cache(:), dvcp_cache(:)
   real(8), allocatable :: dh_cache(:, :)      ! (2, ncoord) - per-orbital h kernels

   ! -- diagnostics (read back by the criterion program / future recorder columns) --
   real(8) :: diag_e_elec = 0.0d0, diag_e_mf = 0.0d0
   real(8) :: diag_imp_occ(2) = 0.0d0, diag_mag_mom = 0.0d0

contains

   !------------------------------------------------------------------
   ! gauss_legendre(n, x, w) - nodes/weights of the Gauss-Legendre rule
   !                 on [-1, 1] (Newton iteration on the Legendre
   !                 recurrence; checked against closed forms by
   !                 criterion (a))
   !------------------------------------------------------------------
   subroutine gauss_legendre(n, x, w)
      integer, intent(in)  :: n
      real(8), intent(out) :: x(n), w(n)
      integer :: i, it, k
      real(8) :: z, z1, pm1, p0, p1, dp
      do i = 1, (n + 1)/2
         ! 1. initial guess (the Chebyshev-like node seed)
         z = cos(pi*(real(i, 8) - 0.25d0)/(real(n, 8) + 0.5d0))
         ! 2. Newton on P_n(z) = 0 through the three-term recurrence
         do it = 1, 100
            pm1 = 1.0d0
            p0 = z
            do k = 2, n
               p1 = ((2.0d0*real(k, 8) - 1.0d0)*z*p0 - (real(k, 8) - 1.0d0)*pm1)/real(k, 8)
               pm1 = p0
               p0 = p1
            end do
            dp = real(n, 8)*(z*p0 - pm1)/(z*z - 1.0d0)
            z1 = z
            z = z1 - p0/dp
            if (abs(z - z1) < 1.0d-15) exit
         end do
         ! 3. weight + symmetric mirror (ascending node order)
         x(i) = -z
         x(n + 1 - i) = z
         w(i) = 2.0d0/((1.0d0 - z*z)*dp*dp)
         w(n + 1 - i) = w(i)
      end do
   end subroutine gauss_legendre

   !------------------------------------------------------------------
   ! cau_rst_tdhf_setup(norb, bw, telec, start, nsub, pmode, pamp) -
   !                 one-time bath/model initialization; arms the package
   !------------------------------------------------------------------
   subroutine cau_rst_tdhf_setup(norb, bw, telec, start, nsub, pmode, pamp)
      integer, intent(in) :: norb, nsub, pmode
      real(8), intent(in) :: bw, telec, pamp
      character(len=*), intent(in) :: start
      integer :: i, nh
      ! 0. re-entry hygiene (a second setup - the criterion program's
      !    small-bath probe - deallocates the previous model first)
      if (allocated(e_levels)) then
         deallocate (e_levels, w_levels, e_metal, hcore, fock, e_hf, c_hf, rhom)
      end if
      ! 1. validate + adopt the model parameters
      if (mod(norb, 2) /= 0 .or. norb < 2) then
         call stop_td('cau_rst_tdhf_setup', 'N_ORB must be even and >= 2 (the bath comes in +/- pairs)')
      end if
      n_orb = norb
      n_siz = n_orb + n_sptl
      band_width = bw
      t_elec = telec
      start_state = start
      msbstp = max(1, nsub)
      perturb_mode = pmode
      perturb_amp = pamp
      perturb_done = .false.
      ! 2. quadrature + bath energies (the old Init_IESH mapping verbatim:
      !    node i fills level nh+1-i below the Fermi level and nh+i above)
      nh = n_orb/2
      allocate (e_levels(nh), w_levels(nh), e_metal(n_orb))
      call gauss_legendre(nh, e_levels, w_levels)
      do i = 1, nh
         e_metal(nh + 1 - i) = -band_width/2.0d0*(1.0d0 + e_levels(i))/2.0d0
         e_metal(nh + i) = band_width/2.0d0*(1.0d0 + e_levels(i))/2.0d0
      end do
      energy_metal_zero = 0.0d0
      do i = 1, nh
         energy_metal_zero = energy_metal_zero - band_width/2.0d0*(1.0d0 + e_levels(i))/2.0d0
      end do
      energy_metal_zero = 2.0d0*energy_metal_zero
      ! 3. allocate the electronic state and seed it (diabatic constants)
      allocate (hcore(n_siz, n_siz, 2), fock(n_siz, n_siz, 2))
      allocate (e_hf(n_siz, 2), c_hf(n_siz, n_siz, 2), rhom(n_siz, n_siz, 2))
      call seed_rho()
      armed = .true.
   end subroutine cau_rst_tdhf_setup

   logical function cau_rst_tdhf_armed()
      cau_rst_tdhf_armed = armed
   end function cau_rst_tdhf_armed

   !------------------------------------------------------------------
   ! cau_rst_tdhf_end_traj() - trajectory-end rho reseed (called by the
   !                 container's own true termination judgment; the seed
   !                 is q-independent and idempotent)
   !------------------------------------------------------------------
   subroutine cau_rst_tdhf_end_traj()
      if (armed) call seed_rho()
   end subroutine cau_rst_tdhf_end_traj

   !------------------------------------------------------------------
   ! seed_rho() - the q-independent initial diabatic occupation (the old
   !                 set_occ: impurity per start state + Fermi-Dirac bath)
   !------------------------------------------------------------------
   subroutine seed_rho()
      integer :: i, is
      real(8) :: kT, prob
      rhom = (0.0d0, 0.0d0)
      ! 1. impurity occupation (the old convention: 'GS' = both impurity
      !    electrons in spin 1; '1D' = one electron in each spin channel)
      if (start_state == 'GS') then
         rhom(1, 1, 1) = (1.0d0, 0.0d0)
         rhom(2, 2, 1) = (1.0d0, 0.0d0)
      else if (start_state == '1D') then
         rhom(1, 1, 1) = (1.0d0, 0.0d0)
         rhom(2, 2, 2) = (1.0d0, 0.0d0)
      else
         call stop_td('seed_rho', 'START_U0 must be GS or 1D (unknown impurity start state "'// &
                      trim(start_state)//'")')
      end if
      ! 2. bath: Fermi-Dirac at mu_ev with kT = kb_ev*T [eV]
      kT = kb_ev*t_elec
      do is = 1, 2
         do i = n_sptl + 1, n_siz
            if (kT > 1.0d-8) then
               prob = 1.0d0/(1.0d0 + exp((e_metal(i - n_sptl) - mu_ev)/kT))
            else
               if (e_metal(i - n_sptl) < mu_ev) then
                  prob = 1.0d0
               else
                  prob = 0.0d0
               end if
            end if
            rhom(i, i, is) = cmplx(prob, 0.0d0, 8)
         end do
      end do
      ! 3. reset the trajectory-time memory (the prop guard re-learns it)
      t_seen = -1.0d30
   end subroutine seed_rho

   !------------------------------------------------------------------
   ! hubbard_eval(q, want_grad) - one ground + Hubbard delivery at q,
   !                 cached module-side (values always; kernels on request)
   !------------------------------------------------------------------
   subroutine hubbard_eval(q, want_grad)
      real(8), intent(in) :: q(:)
      logical, intent(in) :: want_grad
      integer :: ncoord
      if (.not. allocated(g0_cache)) then
         ncoord = size(q)
         allocate (g0_cache(ncoord), dh_cache(2, ncoord), du_cache(ncoord), &
                   dj_cache(ncoord), dvcp_cache(ncoord))
      end if
      call cau_rst_ground(q, e0_cache, g0_cache)
      if (want_grad) then
         call cau_rst_hubbard(q, h_vec, u_val, j_val, vg_vec(1), &
                               dh_cache, du_cache, dj_cache, dvcp_cache)
         vg_vec(2) = vg_vec(1)
      else
         call cau_rst_hubbard_val(q, h_vec, u_val, j_val, vg_vec(1))
         vg_vec(2) = vg_vec(1)
      end if
   end subroutine hubbard_eval

   !------------------------------------------------------------------
   ! assemble_hcore() - fill hcore from the cached h_vec/vg_vec delivery
   !                 (impurity on-sites, bath diagonals, coupling blocks
   !                 sqrt(w/2)*Vg - the old buildHcore body, in eV)
   !------------------------------------------------------------------
   subroutine assemble_hcore()
      integer :: i, ii, is, nh
      nh = n_orb/2
      hcore = 0.0d0
      do ii = 1, n_sptl
         hcore(ii, ii, 1) = h_vec(ii)
         hcore(ii, ii, 2) = h_vec(ii)
      end do
      do is = 1, 2
         do i = 1, nh
            hcore(n_sptl + i, n_sptl + i, is) = e_metal(i)
            hcore(n_sptl + nh + i, n_sptl + nh + i, is) = e_metal(nh + i)
            do ii = 1, n_sptl
               hcore(ii, n_sptl + i, is) = sqrt(w_levels(i)/2.0d0)*vg_vec(ii)
               hcore(n_sptl + i, ii, is) = hcore(ii, n_sptl + i, is)
               hcore(ii, n_sptl + nh + i, is) = hcore(ii, n_sptl + i, is)
               hcore(n_sptl + nh + i, ii, is) = hcore(ii, n_sptl + i, is)
            end do
         end do
      end do
   end subroutine assemble_hcore

   !------------------------------------------------------------------
   ! build_hcore_prop(q) - the prop-path Hcore: delivery at q + the keyed
   !                 h degeneracy break (the old perturbation; the rng
   !                 draw is the container's internal stochastic law)
   !------------------------------------------------------------------
   subroutine build_hcore_prop(q)
      real(8), intent(in) :: q(:)
      call hubbard_eval(q, .false.)
      if (perturb_mode == 1 .or. (perturb_mode == 2 .and. .not. perturb_done)) then
         h_vec(1) = h_vec(1) + (rng_u() - 0.5d0)*2.0d0*perturb_amp
         h_vec(2) = h_vec(2) - (rng_u() - 0.5d0)*2.0d0*perturb_amp
         perturb_done = .true.
      end if
      call assemble_hcore()
   end subroutine build_hcore_prop

   !------------------------------------------------------------------
   ! build_fock() - mean-field Fock assembly from the CURRENT rhom +
   !                 diagonalization per spin (the old buildFock: W/J
   !                 tensors from the scalar U/J delivery, eig_sym)
   !------------------------------------------------------------------
   subroutine build_fock()
      integer :: i, is
      real(8) :: nv(2, 2), vsm(2, 2)
      real(8) :: wmat(2, 2), jmat(2, 2), wmj(2, 2)
      ! 1. impurity occupation vector (|Re rho_ii|) + interaction tensors
      call occ_vec(nv)
      wmat = u_val
      jmat(1, 1) = u_val
      jmat(2, 2) = u_val
      jmat(1, 2) = j_val
      jmat(2, 1) = j_val
      wmj = wmat - jmat
      ! 2. per-spin shift vectors (the original vs_hf(:,1)/(:,2)): spin 1
      !    carries W*n_2 + (W-J)*n_1 on BOTH orbitals, spin 2 mirrored -
      !    PER-SPIN VECTORS (the former two-scalar form took component 1 of
      !    the spin-1 vector and component 2 of the spin-2 vector and applied
      !    the pair to both spins - cross-wiring the channels whenever the
      !    spin occupations differ; F was then no longer dE_mf/dn and the
      !    mean-field energy theorem failed structurally, 2026-09-28)
      vsm(:, 1) = matmul(wmat, nv(:, 2)) + matmul(wmj, nv(:, 1))
      vsm(:, 2) = matmul(wmat, nv(:, 1)) + matmul(wmj, nv(:, 2))
      ! 3. Fock = Hcore + shift on the impurity diagonal; diagonalize
      do is = 1, 2
         fock(:, :, is) = hcore(:, :, is)
         do i = 1, n_sptl
            fock(i, i, is) = fock(i, i, is) + vsm(i, is)
         end do
         c_hf(:, :, is) = fock(:, :, is)
         call eig_sym(c_hf(:, :, is), e_hf(:, is))
      end do
   end subroutine build_fock

   !------------------------------------------------------------------
   ! evolve_rho(dt) - one propagation of rhom with the CURRENT Fock
   !                 eigenpairs (rho -> C e^{iEdt} (C^T rho C) e^{-iEdt} C^T;
   !                 phase exp(-i*E[eV]*dt[fs]/hbar_evfs[eV*fs]) - the
   !                 2026-10-06 eV-A-amu-fs base)
   !------------------------------------------------------------------
   subroutine evolve_rho(dt)
      real(8), intent(in) :: dt                ! substep length [fs]
      integer :: n, is, i
      complex(8), allocatable :: cc(:, :), tmp(:, :), rho_adi(:, :)
      complex(8) :: ph(size(rhom, 1))
      n = size(rhom, 1)
      allocate (cc(n, n), tmp(n, n), rho_adi(n, n))
      do is = 1, 2
         ! 1. rotate to the adiabatic representation
         cc = cmplx(c_hf(:, :, is), 0.0d0, 8)
         rho_adi = matmul(conjg(transpose(cc)), matmul(rhom(:, :, is), cc))
         ! 2. free propagation phases on the diagonal
         do i = 1, n
            ph(i) = exp(cmplx(0.0d0, 1.0d0, 8)*e_hf(i, is)*dt/hbar_evfs)
         end do
         do i = 1, n
            rho_adi(i, :) = rho_adi(i, :)*ph(i)
            rho_adi(:, i) = rho_adi(:, i)*conjg(ph(i))
         end do
         ! 3. rotate back
         tmp = matmul(cc, matmul(rho_adi, conjg(transpose(cc))))
         rhom(:, :, is) = tmp
      end do
      deallocate (cc, tmp, rho_adi)
   end subroutine evolve_rho

   !------------------------------------------------------------------
   ! occ_vec(nv) - impurity occupation snapshot [|Re rho_ii|]
   !------------------------------------------------------------------
   subroutine occ_vec(nv)
      real(8), intent(out) :: nv(2, 2)
      integer :: i
      do i = 1, n_sptl
         nv(i, 1) = abs(real(rhom(i, i, 1)))
         nv(i, 2) = abs(real(rhom(i, i, 2)))
      end do
   end subroutine occ_vec

   !------------------------------------------------------------------
   ! mf_energy(nv) - the mean-field two-electron energy (the old E_mf):
   !                 n1.W.n2 + 1/2 n1.(W-J).n1 + 1/2 n2.(W-J).n2 [eV]
   !------------------------------------------------------------------
   function mf_energy(nv) result(e)
      real(8), intent(in) :: nv(2, 2)
      real(8) :: e
      real(8) :: wmat(2, 2), jmat(2, 2), wmj(2, 2)
      wmat = u_val
      jmat(1, 1) = u_val
      jmat(2, 2) = u_val
      jmat(1, 2) = j_val
      jmat(2, 1) = j_val
      wmj = wmat - jmat
      e = dot_product(nv(:, 1), matmul(wmat, nv(:, 2))) &
        + 0.5d0*dot_product(nv(:, 1), matmul(wmj, nv(:, 1))) &
        + 0.5d0*dot_product(nv(:, 2), matmul(wmj, nv(:, 2)))
   end function mf_energy

   !------------------------------------------------------------------
   ! mf_kernel(k, nv, out) - dE_mf/dq_k [eV/A] through the U/J kernels
   !                 (dW = dU everywhere; dJ = dU diagonal, dJ off-diag)
   !------------------------------------------------------------------
   subroutine mf_kernel(k, nv, out)
      integer, intent(in) :: k
      real(8), intent(in) :: nv(2, 2)
      real(8), intent(out) :: out
      real(8) :: dw(2, 2), dj(2, 2), dwmj(2, 2)
      dw = du_cache(k)
      dj(1, 1) = du_cache(k)
      dj(2, 2) = du_cache(k)
      dj(1, 2) = dj_cache(k)
      dj(2, 1) = dj_cache(k)
      dwmj = dw - dj
      out = dot_product(nv(:, 1), matmul(dw, nv(:, 2))) &
          + 0.5d0*dot_product(nv(:, 1), matmul(dwmj, nv(:, 1))) &
          + 0.5d0*dot_product(nv(:, 2), matmul(dwmj, nv(:, 2)))
   end subroutine mf_kernel

   !------------------------------------------------------------------
   ! cau_rst_pes(q, mass, g, v) - the container's composed force entry
   !                 (the force-slot contract: g = dV/dq[eV/A], v[eV],
   !                 both fully assigned; the adiabatic degenerate point
   !                 delivers the pure EMT ground PES)
   !------------------------------------------------------------------
   subroutine cau_rst_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      integer :: k, i, is, ii, nh
      real(8) :: nv(2, 2), e_elec, e_mf, fk, tr, fmf
      ! 1. unarmed degenerate point: pure EMT ground PES (matches the old
      !    CALTYP = -1 force exactly; the Hubbard fold never activates)
      if (.not. armed) then
         call hubbard_eval(q, .false.)
         g = g0_cache
         v = e0_cache
         return
      end if
      ! 2. armed fold: full delivery (ground + Hubbard kernels) at q
      call hubbard_eval(q, .true.)
      call occ_vec(nv)
      ! 3. Hcore at the unperturbed delivery (the pes/prop pair must come
      !    from ONE potential function; the keyed perturbation lives only
      !    on the prop path, as in the old code)
      call assemble_hcore()
      e_elec = 0.0d0
      do is = 1, 2
         do i = 1, n_siz
            do k = 1, n_siz
               e_elec = e_elec + real(rhom(i, k, is))*hcore(k, i, is)
            end do
         end do
      end do
      e_mf = mf_energy(nv)
      ! 4. composed potential: V = E0 + E_elec + E_mf - energy_metal_zero [eV]
      v = e0_cache + e_elec + e_mf - energy_metal_zero
      ! 5. composed kernel: g = dE0/dq + sum_is Tr(rho.dHcore/dq) + dE_mf/dq
      !    (impurity diagonal dh for both spatial orbitals; coupling blocks
      !    dVcp*sqrt(w/2) for both impurity rows and both Fermi branches)
      nh = n_orb/2
      g = g0_cache
      do k = 1, size(q)
         fk = 0.0d0
         do is = 1, 2
            do ii = 1, n_sptl
               fk = fk + real(rhom(ii, ii, is))*dh_cache(ii, k)
            end do
            do i = 1, nh
               tr = sqrt(w_levels(i)/2.0d0)*dvcp_cache(k)
               do ii = 1, n_sptl
                  ! both symmetric elements of each coupling block
                  ! contribute (the full-trace convention: Grad(bath,imp) =
                  ! Grad(imp,bath) - the partner element was missing here,
                  ! halving the coupling force vs the energy's derivative)
                  fk = fk + real(rhom(ii, n_sptl + i, is))*tr
                  fk = fk + real(rhom(n_sptl + i, ii, is))*tr
                  fk = fk + real(rhom(ii, n_sptl + nh + i, is))*tr
                  fk = fk + real(rhom(n_sptl + nh + i, ii, is))*tr
               end do
            end do
         end do
         call mf_kernel(k, nv, fmf)
         g(k) = g(k) + fk + fmf
      end do
      ! 6. diagnostics
      diag_e_elec = e_elec
      diag_e_mf = e_mf
      diag_imp_occ(1) = nv(1, 1) + nv(2, 1)
      diag_imp_occ(2) = nv(1, 2) + nv(2, 2)
      diag_mag_mom = diag_imp_occ(1) - diag_imp_occ(2)
   end subroutine cau_rst_pes

   !------------------------------------------------------------------
   ! cau_rst_prop(sta, dt) - the 'tdhf_hub' propagation member: per call
   !                 msbstp Magnus-2 predictor/corrector substeps at the
   !                 fixed endpoint geometry sta%q (the old TDHF_Evolution
   !                 inner loop); includes the time-backwards trajectory
   !                 guard (t reset = a new trajectory -> reseed rho)
   !------------------------------------------------------------------
   subroutine cau_rst_prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
      integer :: i, is
      real(8) :: dtq
      complex(8), allocatable :: rho_save(:, :, :)
      real(8), allocatable :: fock_t(:, :, :)
      if (.not. armed) return                      ! degenerate member (not selected)
      ! 1. trajectory guard: the trajectory time went backwards (samp_run
      !    reset it) - a new trajectory, reseed the diabatic occupation
      if (sta%t < t_seen - 0.5d0*max(dt, 1.0d-12)) call seed_rho()
      t_seen = sta%t
      ! 2. fixed Hcore for the whole call (endpoint geometry)
      call build_hcore_prop(sta%q)
      allocate (rho_save(n_siz, n_siz, 2), fock_t(n_siz, n_siz, 2))
      dtq = dt/real(msbstp, 8)
      ! 3. the old Magnus-2 sequence per substep: t-Fock -> evolve guess ->
      !    guess-Fock -> average -> re-diagonalize -> evolve the original
      do i = 1, msbstp
         rho_save = rhom
         call build_fock()
         fock_t = fock
         call evolve_rho(dtq)
         call build_fock()
         rhom = rho_save
         fock = 0.5d0*(fock_t + fock)
         do is = 1, 2
            c_hf(:, :, is) = fock(:, :, is)
            call eig_sym(c_hf(:, :, is), e_hf(:, is))
         end do
         call evolve_rho(dtq)
      end do
      deallocate (rho_save, fock_t)
   end subroutine cau_rst_prop

   !------------------------------------------------------------------
   ! cau_rst_match(sta) - the 'tdhf_hub' matching member: the mean-field
   !                 scheme has no hops (occ stays the seeded [1], n_hop
   !                 stays 0, no rng consumption) - a pure no-op member
   !------------------------------------------------------------------
   subroutine cau_rst_match(sta)
      type(state_t), intent(inout) :: sta
      ! the Ehrenfest/TDHF statistic: nothing to sample, nothing to write
      ! back (the impurity observables travel through cau_rst_diag_out)
   end subroutine cau_rst_match

   !------------------------------------------------------------------
   ! cau_rst_diag_out(...) - read back the last pes-call diagnostics
   !                 (criterion program / future recorder columns)
   !------------------------------------------------------------------
   subroutine cau_rst_diag_out(e_elec, e_mf, imp_up, imp_dn, mag)
      real(8), intent(out) :: e_elec, e_mf, imp_up, imp_dn, mag
      e_elec = diag_e_elec
      e_mf = diag_e_mf
      imp_up = diag_imp_occ(1)
      imp_dn = diag_imp_occ(2)
      mag = diag_mag_mom
   end subroutine cau_rst_diag_out

   !------------------------------------------------------------------
   ! the criterion-introspection outs (check_cau_rst legs (a)(b)(e)(f)(g)
   !                 (d)(h) - read-only views of the module state)
   !------------------------------------------------------------------
   subroutine bath_out(em, ez)
      real(8), intent(out) :: em(:), ez
      em = e_metal
      ez = energy_metal_zero
   end subroutine bath_out

   subroutine rho_seed_out(occ, tr)
      real(8), intent(out) :: occ(2), tr
      integer :: i, is
      do is = 1, 2
         occ(is) = real(rhom(1, 1, is)) + real(rhom(2, 2, is))
      end do
      tr = 0.0d0
      do is = 1, 2
         do i = 1, n_siz
            tr = tr + real(rhom(i, i, is))
         end do
      end do
   end subroutine rho_seed_out

   subroutine rho_herm_out(tr, herm)
      real(8), intent(out) :: tr, herm
      integer :: i, is
      tr = 0.0d0
      herm = 0.0d0
      do is = 1, 2
         do i = 1, n_siz
            tr = tr + real(rhom(i, i, is))
            herm = max(herm, abs(rhom(i, i, is) - conjg(rhom(i, i, is))))
         end do
         herm = max(herm, maxval(abs(rhom(:, :, is) - conjg(transpose(rhom(:, :, is))))))
      end do
   end subroutine rho_herm_out

   subroutine hcore_out(hc)
      real(8), intent(out) :: hc(:, :, :)
      hc = hcore
   end subroutine hcore_out

   subroutine kernel_caches_out(dh, du, dj, dvcp)
      real(8), intent(out) :: dh(:,:), du(:), dj(:), dvcp(:)
      if (.not. allocated(dh_cache)) then
         dh = 0.0d0
         du = 0.0d0
         dj = 0.0d0
         dvcp = 0.0d0
         return
      end if
      dh = dh_cache
      du = du_cache
      dj = dj_cache
      dvcp = dvcp_cache
   end subroutine kernel_caches_out

   subroutine hub_vals_out(h, u, j, vcp)
      real(8), intent(out) :: h(2), u, j, vcp
      h = h_vec
      u = u_val
      j = j_val
      vcp = vg_vec(1)
   end subroutine hub_vals_out

   !------------------------------------------------------------------
   ! cau_rst_sync(s) - the package's container-publish member: copy the
   !                 module rho into sta%s%rho as spin blocks ((is-1)*n+1
   !                 .. is*n on both axes, spin 1 then 2 - the member's
   !                 declared index semantics; the elec_prop /
   !                 mqc_statistic sync tails call this every step, so
   !                 the recorder sees the post-propagation state)
   !------------------------------------------------------------------
   subroutine cau_rst_sync(s)
      type(elec_t), intent(inout) :: s
      integer :: is, n
      if (.not. armed) return                      ! degenerate posture: nothing published
      n = n_siz
      if (.not. allocated(s%rho) .or. size(s%rho, 1) /= 2*n .or. size(s%rho, 2) /= 2*n) then
         if (allocated(s%rho)) deallocate (s%rho)
         allocate (s%rho(2*n, 2*n))
      end if
      s%rho = (0.0d0, 0.0d0)
      do is = 1, 2
         s%rho((is - 1)*n + 1:is*n, (is - 1)*n + 1:is*n) = rhom(:, :, is)
      end do
      s%n_surf = 2*n                               ! configuration cardinality (spin-orbitals)
   end subroutine cau_rst_sync

   !------------------------------------------------------------------
   ! the recorder columns (col_val_i conformance checked at the
   ! rec_reg_col call site; each reads the LIVE module state - fresher
   ! than the pes-time caches by half a step at rec_frame time)
   !------------------------------------------------------------------
   real(8) function cau_rst_col_nup(sta)
      type(state_t), intent(in) :: sta
      cau_rst_col_nup = real(rhom(1, 1, 1)) + real(rhom(2, 2, 1))
   end function cau_rst_col_nup

   real(8) function cau_rst_col_ndn(sta)
      type(state_t), intent(in) :: sta
      cau_rst_col_ndn = real(rhom(1, 1, 2)) + real(rhom(2, 2, 2))
   end function cau_rst_col_ndn

   real(8) function cau_rst_col_mag(sta)
      type(state_t), intent(in) :: sta
      cau_rst_col_mag = cau_rst_col_nup(sta) - cau_rst_col_ndn(sta)
   end function cau_rst_col_mag

   real(8) function cau_rst_col_eelec(sta)
      type(state_t), intent(in) :: sta
      integer :: i, k, is
      cau_rst_col_eelec = 0.0d0
      if (.not. armed) return
      do is = 1, 2
         do i = 1, n_siz
            do k = 1, n_siz
               cau_rst_col_eelec = cau_rst_col_eelec + real(rhom(i, k, is))*hcore(k, i, is)
            end do
         end do
      end do
   end function cau_rst_col_eelec

   real(8) function cau_rst_col_emf(sta)
      type(state_t), intent(in) :: sta
      real(8) :: nv(2, 2)
      cau_rst_col_emf = 0.0d0
      if (.not. armed) return
      call occ_vec(nv)
      cau_rst_col_emf = mf_energy(nv)
   end function cau_rst_col_emf

   !------------------------------------------------------------------
   ! rho_ref_out(r) - reference copy of the module rho (criterion leg)
   !------------------------------------------------------------------
   subroutine rho_ref_out(r)
      complex(8), intent(out), allocatable :: r(:, :, :)
      allocate (r(n_siz, n_siz, 2))
      r = rhom
   end subroutine rho_ref_out

   ! stop_td(who, msg) - the named-abort channel (stderr + STOP 1)
   subroutine stop_td(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_rst_tdhf error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container evolution layer is not usable)'
      stop 1
   end subroutine stop_td
end module cau_rst_tdhf
