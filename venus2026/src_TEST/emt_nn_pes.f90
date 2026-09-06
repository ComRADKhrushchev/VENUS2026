! ================================================================
! PORT NOTE (venus2026_stage, 2026-09-01)
!   Verbatim copy of ../2D/src_NN/emt_nn_pes_venus.f90 (module
!   `emt_nn_pes_venus`, name kept for 1:1 traceability with the fit
!   side) into the TEST-PES framework. Exposed to the runtime keyword
!   TEST_PES=EMT-NN via src_TEST/test_potentials.f90; the
!   POTPRE/POT0/DPESHON contract is served by src_TEST/interface_TEST.f90
!   (V[eV] -> kcal/mol x23.0605, forces[eV/A] -> PDOT x23.0605*C1).
!   Data assets (Slab.xyz, nn_weights_emt_nn.txt, emt_2spec_params.txt)
!   live in data/emt_nn/; runtime paths are keyword-overridable
!   (EMTNN_SLAB_FILE / EMTNN_WEIGHTS_FILE, default = CWD).
!   Any numerical change here MUST be re-synced with the fit-side
!   original — the body below is untouched transplant.
! ================================================================
! emt_nn_pes.f90 (module emt_nn_pes_venus) — EMT-NN PES engine
!
! Direct transplant of the AnalyticModel fit_emt_nn kernel
! (pot_emt.f90 compute_emt_energy + Hubbard NN), kept 1:1 in sync
! with the fitting program so fitted weights load unchanged.
!
!   E_total = Ecoh2 + Ecoh1 - V22 - V_NN_eff + 0.5*(vref2+vref1) + E0_nn
!   V_NN_eff = sum_j theta(r_j) * NN_pair(r_j)
!   E_int    = E_total - E_slab_bare_total   (bare baseline frozen at init)
!   h/U/J/Vcp = NN_k(s_ad) + sum_j theta_soft(r_j)*NN_pair_k(r_j)
!
!   Ecoh1: analytic EMT path by default; if the weight file carries a
!   ## Ecoh_NN section it REPLACES Ecoh1 = NN_ecoh(s_ad) (4b e0fit family,
!   use_ecoh_nn=true on the fit side). emtnn_use_ecoh_nn flips on at init.
!
! !! CONSTANTS ARE LOCKED TO THE FIT SIDE !!
!   acut = log(10)/(rr - rcut)   — ln(10) convention.
!   (pes_emt_mod.f90 uses smooth_factor=1e4 => acut = ln(10^4), a
!    STEEPER cutoff. The fit parameters were fitted under ln(10);
!   never copy cutoff constants from pes_emt_mod into this engine.)
!
! Public API (interface_EMT_2spec style):
!   emt_nn_init(slab_file, weights_file, a_lat, n_side)  ! failure => STOP
!   calc_emt_nn_energy(natoms, natoma, Q_emt, E, forces)  ! E[eV], forces(3,natoms)[eV/A]
!   calc_hubbard(Q_imp, Q_slab, n_slab, h, U, J, Vcp)     ! eV, current slab coords
!   fill_hubbard_nn_grad(grad, ncoord, ipoly, Q_imp, Q_slab, n_slab, natoma)
!   calc_hubbard_grad_fd(...)          ! FD fallback behind use_hubb_grad_fd
!   emt_nn_loaded / emt_nn_hubbard_loaded / emt_nn_a_lat / emt_nn_box /
!   emt_nn_natoms / emt_nn_rcut / emt_nn_acut / emt_nn_s_ad
!
! Units: E in eV, forces in eV/Angstrom with f = -dE/dR.
! Slab coordinates are taken from the Q arrays at every call
! (movable slab); the bare-slab baseline is frozen at init.
! Single-adsorbate kernel: natoma must be 1.
! ================================================================
module emt_nn_pes_venus
    implicit none
    private

    ! ================================================================
    !  Types (verbatim from fit_emt_nn/pot_emt.f90)
    ! ================================================================
    type :: rst_slab_data
        integer :: n_atoms = 0, n_surf = 0
        real(8), allocatable :: coords(:,:)
        integer, allocatable :: layer_id(:)
        real(8) :: box_lx = 36.0d0, box_ly = 36.0d0, box_skew = 1.5707963267948966d0
    end type rst_slab_data

    type :: emt_params
        real(8) :: eta2   = 0.0d0   ! eta2 : electron density decay (1/A)
        real(8) :: n0     = 0.0d0   ! n0   : equilibrium density (1/A^3)
        real(8) :: e0     = 0.0d0   ! E0   : cohesive energy scale (eV)
        real(8) :: lambda = 0.0d0   ! lambda : well-width parameter (1/A)
        real(8) :: v0     = 0.0d0   ! V0   : pair-interaction strength (eV)
        real(8) :: kappa  = 0.0d0   ! kappa: pair-interaction decay (1/A)
        real(8) :: s0     = 0.0d0   ! s0   : neutral-sphere radius (A)
    end type emt_params

    type :: emt_component
        real(8) :: E_coh_slab = 0.0d0
        real(8) :: V_ref_slab = 0.0d0
        real(8) :: E_coh_ad   = 0.0d0
        real(8) :: V_ref_ad   = 0.0d0
        real(8) :: V_AuAu     = 0.0d0
        real(8) :: V_AuAd     = 0.0d0
    end type emt_component

    ! ================================================================
    !  Universal constants (fit side, verbatim)
    ! ================================================================
    real(8), parameter :: PI = 3.14159265358979323846d0
    real(8), parameter :: sr3 = 1.7320508075688772d0
    real(8), parameter :: sr2o3 = 0.8164965809277260d0
    real(8), parameter :: beta = (16.0d0*acos(-1.0d0)/3.0d0)**(1.0d0/3.0d0) / sqrt(2.0d0)
    real(8), parameter :: bohr2ang = 0.52917721092d0
    real(8), parameter :: sqrt2 = sqrt(2.0d0)
    real(8), parameter :: sqrt3 = sqrt(3.0d0)
    real(8), parameter :: twelfth = 1.0d0 / 12.0d0
    real(8), parameter :: tolerance = 1.0d-30
    integer, parameter :: n_emt_param = 7

    ! Gentle cutoff for Hubbard pairwise NN (fit: acut_soft = 0.3)
    real(8), parameter :: acut_soft = 0.3d0

    ! ================================================================
    !  Module state (fit side, verbatim where applicable)
    ! ================================================================
    type(rst_slab_data), save :: slab

    real(8), save :: E_slab_bare_total   = 0.0d0
    real(8), save :: E_coh_bare_sum      = 0.0d0
    real(8), save :: V_ref_bare_raw      = 0.0d0
    real(8), save :: V_AuAu_bare_raw     = 0.0d0

    ! Au substrate parameters (JCP 143, 124708, Table II)
    type(emt_params), save :: Au = emt_params( &
        eta2=3.197d0, n0=0.042d0, e0=-3.826d0, &
        lambda=4.182d0, v0=0.348d0, kappa=3.249d0, s0=1.642d0)

    real(8), save :: betas0_1, betas0_2
    real(8), save :: betaeta2_1, betaeta2_2
    real(8), save :: kappadbeta_1, kappadbeta_2
    real(8), save :: rnn_1(3), rnn_2(3)

    real(8), save :: rcut, acut
    real(8), save :: x_2(3)
    real(8), save :: igamma1_2, igamma2_2

    ! NN architecture (fit defaults; header line overrides on load)
    integer, save :: emtnn_nh_pair = 8
    real(8), save :: emtnn_r_scale_pair = 3.0d0
    integer, save :: emtnn_nh_hubb = 6
    real(8), save :: emtnn_hubb_scale = 1.0d0
    integer, save :: emtnn_nh_hubb_pair = 3
    real(8), save :: emtnn_r_scale_hubb_pair = 3.0d0
    integer, save :: emtnn_nh_ecoh = 8
    real(8), save :: emtnn_ecoh_scale = 1.0d0
    logical, save :: emtnn_use_ecoh_nn = .false.
    logical, save :: ecoh_weights_loaded = .false.

    ! V_rep: near-field exponential repulsion, dynamics-only fix (2026-08-25).
    ! The pair NN extrapolates spuriously attractive for r < ~1.8 A, giving
    ! Bridge/HCP/FCC subsurface artifact wells (z ~ 0.4-0.8) deeper than the
    ! physical adsorption well -> trajectories tunnel into the slab. Add
    ! V_rep(r) = A*exp(-(r-r0)/rho) per C-Au pair, cosine-tapered to EXACTLY
    ! zero at rc_vrep, so all r >= rc_vrep geometries stay bit-identical to
    ! the fit engine (parity preserved; validated by scripts/vrep_validate.py:
    ! all 4 sites PASS, global min back at physical z 1.95-2.25, well shift 0).
    logical,  save :: emtnn_vrep_on = .true.
    real(8), parameter :: emtnn_vrep_a   = 0.7d0    ! eV
    real(8), parameter :: emtnn_vrep_r0  = 2.0d0    ! A
    real(8), parameter :: emtnn_vrep_rho = 0.094d0  ! A
    real(8), parameter :: emtnn_vrep_rc  = 2.15d0   ! A, exact-zero beyond

    ! GS weight vector workspace (7 EMT + pair NN + E0 [+ ecoh NN])
    real(8), allocatable, save :: fit_weights(:)
    logical, save :: have_fit_weights = .false.

    ! Extracted GS parameters (init pulls them from fit_weights once)
    real(8), save :: ad7_w(7) = 0.0d0
    real(8), save :: E0_w = 0.0d0
    real(8), allocatable, save :: pair_w(:)
    real(8), allocatable, save :: ecoh_w(:)   ! Ecoh NN weights (allocated only if file has ## Ecoh_NN)

    ! Hubbard NN weights: hubb_nn_weights(:,1:4) = h, U, J, Vcp
    real(8), allocatable, save :: hubb_nn_weights(:,:)
    real(8), allocatable, save :: hubb_pair_nn_weights(:,:)

    ! Cache of last-computed Hubbard parameters
    real(8), save :: emt_cached_h = 0.0d0
    real(8), save :: emt_cached_U = 0.0d0
    real(8), save :: emt_cached_J = 0.0d0
    real(8), save :: emt_cached_V = 0.0d0
    real(8), save :: emt_cached_s_ad = 0.0d0

    ! Hard-cutoff range (fit: cutoff_range 1.5)
    real(8), save :: emt_cutoff_range = 1.5d0

    ! Engine bookkeeping
    character(len=256), save :: slab_file_val = 'Slab.xyz'
    logical, save :: engine_loaded = .false.
    logical, save :: hubbard_loaded = .false.
    logical, save :: hubbard_section_seen = .false.
    real(8), save :: a_lat_val = 2.88d0
    integer, save :: n_side_val = 6
    logical, save :: use_hubb_grad_fd = .false.

    public :: emt_nn_init
    public :: calc_emt_nn_energy, calc_hubbard
    public :: fill_hubbard_nn_grad, calc_hubbard_grad_fd, use_hubb_grad_fd
    public :: emt_nn_loaded, emt_nn_hubbard_loaded
    public :: emt_nn_a_lat, emt_nn_box, emt_nn_natoms
    public :: emt_nn_rcut, emt_nn_acut, emt_nn_s_ad
    public :: emt_component

contains

    ! ================================================================
    !  eval_1d_nn — verbatim from fit_emt_nn/nn_engine.f90
    !
    !  weights layout: w1(1:nh), b1(nh+1:2*nh), w2(2*nh+1:3*nh), b2(3*nh+1)
    !  Analytical dV/dz includes the 1/z_scale factor.
    ! ================================================================
    subroutine eval_1d_nn(z, weights, z_scale, V, dVdz)
        real(8), intent(in)  :: z
        real(8), intent(in)  :: weights(:)
        real(8), intent(in)  :: z_scale
        real(8), intent(out) :: V, dVdz

        integer :: i, nh
        real(8) :: w1, b1, w2, h, zs, dh_dz

        nh = (size(weights) - 1) / 3
        zs = z / z_scale

        V = weights(3*nh + 1)
        dVdz = 0.0d0

        do i = 1, nh
            w1 = weights(i)
            b1 = weights(nh + i)
            w2 = weights(2*nh + i)

            h = tanh(w1 * zs + b1)
            V = V + w2 * h
            dh_dz = (1.0d0 - h * h) * w1 / z_scale
            dVdz = dVdz + w2 * dh_dz
        end do
    end subroutine eval_1d_nn

    ! ================================================================
    !  pbc_wrap — verbatim from fit_emt_nn/cutoff_utils.f90
    !
    !  Cell: L = [box_lx, box_ly*cos(skew)] ; skew = pi/3 hexagonal
    ! ================================================================
    subroutine pbc_wrap(dx, dy, box_lx, box_ly, skew)
        real(8), intent(inout) :: dx, dy
        real(8), intent(in) :: box_lx, box_ly, skew

        real(8) :: cskew, sskew, inv_11, inv_12, inv_22, fx, fy

        cskew = cos(skew)
        sskew = sin(skew)

        if (abs(sskew) < 1.0d-12) then
            ! Orthogonal cell — Cartesian fallback
            dx = dx - box_lx * nint(dx / box_lx)
            dy = dy - box_ly * nint(dy / box_ly)
        else
            ! Fractional-coordinate wrapping
            inv_11 = 1.0d0 / box_lx
            inv_12 = -cskew / (box_lx * sskew)
            inv_22 = 1.0d0 / (box_ly * sskew)
            fx = inv_11 * dx + inv_12 * dy
            fy = inv_22 * dy
            fx = fx - nint(fx)
            fy = fy - nint(fy)
            dx = box_lx * fx + box_ly * cskew * fy
            dy = box_ly * sskew * fy
        end if
    end subroutine pbc_wrap

    ! ================================================================
    !  EMT kernel — verbatim transplant of fit_emt_nn/pot_emt.f90
    !                compute_emt_energy (full analytical gradient)
    !
    !  Sign conventions follow the fit original EXACTLY:
    !    f_array = dEcoh_array - dV_array + 0.5*dvref_array
    !  where f_array stores FORCE (= -dE/dr); the returned
    !  dEdx/dEdy/dEdz are ENERGY gradients (= -force).
    !
    !  ad_params(1:7) = (eta2, n0, E0, lambda, V0, kappa, s0)
    !  Uses the module slab (coords synced by the public wrappers).
    ! ================================================================
    subroutine compute_emt_energy(xyz, ad_params, E_int, dEdx, dEdy, dEdz, &
                                  f_slab, components, pair_nn_weights, E0_nn, ecoh_nn_weights)
        real(8), intent(in)  :: xyz(3)           ! adsorbate position (A)
        real(8), intent(in)  :: ad_params(7)     ! adsorbate EMT parameters
        real(8), intent(out) :: E_int            ! interaction energy (eV)
        real(8), intent(out), optional :: dEdx, dEdy, dEdz  ! energy gradient (ad)
        real(8), intent(out), optional :: f_slab(:,:)       ! forces on slab atoms (3,na)
        type(emt_component), intent(out), optional :: components
        real(8), intent(in),  optional :: pair_nn_weights(:) ! pair NN weights (3*nh+1)
        real(8), intent(in),  optional :: E0_nn              ! NN energy offset
        real(8), intent(in),  optional :: ecoh_nn_weights(:) ! Ecoh NN weights

        ! --- species parameters unpacked ---
        real(8) :: eta2_1, n0_1, e0_1, lambda_1, v0_1, kappa_1, s0_1
        real(8) :: eta2_2, n0_2, e0_2, lambda_2, v0_2, kappa_2, s0_2

        ! --- cross-species mixing ---
        real(8) :: chi_21, chi_12

        ! --- cutoff & gamma ---
        real(8) :: x_1(3), r3temp(3)
        real(8) :: igamma1_1, igamma2_1

        ! --- sigma accumulators ---
        real(8), allocatable :: sigma_22(:), sigma_21(:)
        real(8) :: sigma_12

        ! --- pair repulsion accumulators ---
        real(8) :: V_22, V_NN_eff

        ! --- neutral-sphere radius ---
        real(8), allocatable :: s_2(:)
        real(8) :: s_1

        ! --- cohesive energy ---
        real(8) :: Ecoh_1, Ecoh_2

        ! --- Ecoh NN locals ---
        logical  :: use_ec
        real(8)  :: s_ad_ec, dEcoh_ds, coef_ec

        ! --- exp(-lambda*s) workspace ---
        real(8), allocatable :: exp_lambda_s2(:)
        real(8) :: exp_lambda_s1

        ! --- reference pair potential ---
        real(8) :: vref_1, vref_2

        ! --- total energy ---
        real(8) :: E_total

        ! --- gradient arrays (original sign convention: -d/dr) ---
        real(8), allocatable :: dsigma_22(:,:,:)    ! (3, na, na)
        real(8), allocatable :: dsigma_21_2(:,:)    ! (3, na)
        real(8), allocatable :: dsigma_21_1(:,:)    ! (3, na)
        real(8), allocatable :: dsigma_12_1(:)      ! (3)
        real(8), allocatable :: dsigma_12_2(:,:)    ! (3, na)

        real(8), allocatable :: dV_22_2(:,:)        ! (3, na)
        real(8), allocatable :: dV_NN_eff_2(:,:)    ! (3, na)
        real(8), allocatable :: dV_NN_eff_1(:)      ! (3)

        real(8) :: V_nn, dV_nn_dr, dtheta_nn

        ! --- V_rep near-field repulsion (see module header comment) ---
        real(8) :: V_rep_tot, vrep_w, dvrep
        real(8) :: dV_rep_1(3)
        real(8), allocatable :: dV_rep_2(:,:)    ! (3, na)

        real(8), allocatable :: ds_2_2(:,:,:)   ! (3, na, na)
        real(8), allocatable :: ds_2_1(:,:)     ! (3, na)
        real(8), allocatable :: ds_1_2(:,:)     ! (3, na)
        real(8), allocatable :: ds_1_1(:)       ! (3)

        real(8), allocatable :: dEcoh_2_2(:,:)  ! (3, na)
        real(8), allocatable :: dEcoh_1_2(:,:)  ! (3, na)
        real(8), allocatable :: dEcoh_2_1(:)    ! (3)
        real(8), allocatable :: dEcoh_1_1(:)    ! (3)

        real(8), allocatable :: dvref_2_2(:,:)  ! (3, na)
        real(8), allocatable :: dvref_1_2(:,:)  ! (3, na)
        real(8), allocatable :: dvref_2_1(:)    ! (3)
        real(8), allocatable :: dvref_1_1(:)    ! (3)

        ! --- force accumulators (original convention: f = -dE/dr) ---
        real(8), allocatable :: f_ad(:)          ! (3)
        real(8), allocatable :: f_slab_loc(:,:)  ! (3, na)

        ! --- loop / workspace ---
        real(8) :: dx, dy, dz, r, theta, rtemp, rtemp1
        real(8) :: dtheta(3), vec(3)
        integer  :: i, j, na
        logical  :: need_grad

        na = slab%n_atoms
        need_grad = present(dEdx) .or. present(dEdy) .or. present(dEdz) .or. present(f_slab)

        ! --- Phase 1: unpack parameters ---
        eta2_1   = ad_params(1)
        n0_1     = ad_params(2)
        e0_1     = ad_params(3)
        lambda_1 = ad_params(4)
        v0_1     = ad_params(5)
        kappa_1  = ad_params(6)
        s0_1     = ad_params(7)

        eta2_2   = Au%eta2
        n0_2     = Au%n0
        e0_2     = Au%e0
        lambda_2 = Au%lambda
        v0_2     = Au%v0
        kappa_2  = Au%kappa
        s0_2     = Au%s0

        betas0_1     = beta * s0_1
        betaeta2_1   = beta * eta2_1
        kappadbeta_1 = kappa_1 / beta

        rnn_1(1) = betas0_1
        rnn_1(2) = rnn_1(1) * sqrt2
        rnn_1(3) = rnn_1(1) * sqrt3

        chi_21 = n0_1 / n0_2 * exp(0.5d0/bohr2ang * (s0_2 - s0_1))
        chi_12 = 1.0d0 / chi_21

        ! --- Phase 2: adsorbate cutoff & gamma ---
        x_1 = [1.0d0, 0.5d0, 2.0d0] / (1.0d0 + exp(acut*(rnn_1 - rcut)))

        r3temp = rnn_1 - betas0_1
        igamma1_1 = 1.0d0 / sum(x_1 * exp(-eta2_1 * r3temp))
        igamma2_1 = 1.0d0 / sum(x_1 * exp(-kappadbeta_1 * r3temp))

        ! --- Phase 4a: allocate ---
        allocate(sigma_22(na), sigma_21(na), s_2(na), exp_lambda_s2(na))

        if (need_grad) then
            allocate(dsigma_22(3, na, na))
            allocate(dsigma_21_2(3, na), dsigma_21_1(3, na))
            allocate(dsigma_12_1(3), dsigma_12_2(3, na))
            allocate(dV_22_2(3, na))
            allocate(dV_NN_eff_2(3, na), dV_NN_eff_1(3))
            if (emtnn_vrep_on) allocate(dV_rep_2(3, na))
            allocate(ds_2_2(3, na, na), ds_1_1(3))
            allocate(ds_2_1(3, na), ds_1_2(3, na))
            allocate(dEcoh_2_2(3, na), dEcoh_2_1(3))
            allocate(dEcoh_1_2(3, na), dEcoh_1_1(3))
            allocate(dvref_2_2(3, na), dvref_2_1(3))
            allocate(dvref_1_2(3, na), dvref_1_1(3))
            allocate(f_ad(3), f_slab_loc(3, na))
        end if

        ! --- Phase 4b: initialize accumulators ---
        sigma_22 = 0.0d0; sigma_21 = 0.0d0; sigma_12 = 0.0d0
        V_22 = 0.0d0; V_NN_eff = 0.0d0
        V_rep_tot = 0.0d0
        if (need_grad) then
            dV_rep_1 = 0.0d0
            if (emtnn_vrep_on) dV_rep_2 = 0.0d0
        end if

        if (need_grad) then
            dsigma_22   = 0.0d0
            dsigma_21_2 = 0.0d0; dsigma_21_1 = 0.0d0
            dsigma_12_1 = 0.0d0; dsigma_12_2 = 0.0d0
            dV_22_2 = 0.0d0
            dV_NN_eff_2 = 0.0d0; dV_NN_eff_1 = 0.0d0
            ds_2_2  = 0.0d0; ds_1_1  = 0.0d0
            ds_2_1  = 0.0d0; ds_1_2  = 0.0d0
            dEcoh_2_2 = 0.0d0; dEcoh_2_1 = 0.0d0
            dEcoh_1_2 = 0.0d0; dEcoh_1_1 = 0.0d0
            dvref_2_2 = 0.0d0; dvref_2_1 = 0.0d0
            dvref_1_2 = 0.0d0; dvref_1_1 = 0.0d0
            f_ad = 0.0d0; f_slab_loc = 0.0d0
        end if

        ! --- Phase 4c: Au-Au pairs ---
        do i = 1, na - 1
            do j = i + 1, na
                dx = slab%coords(1,i) - slab%coords(1,j)
                dy = slab%coords(2,i) - slab%coords(2,j)
                call pbc_wrap(dx, dy, slab%box_lx, slab%box_ly, slab%box_skew)
                dz = slab%coords(3,i) - slab%coords(3,j)
                r  = sqrt(dx*dx + dy*dy + dz*dz)

                if (r > emt_cutoff_range * rcut) cycle

                if (r > 1.0d-12) then
                    vec = [dx, dy, dz] / r
                else
                    vec = [0.0d0, 0.0d0, 0.0d0]
                end if
                rtemp  = exp(acut*(r - rcut))
                theta  = 1.0d0 / (1.0d0 + rtemp)
                rtemp1 = acut * rtemp * theta

                ! sigma_22
                rtemp = theta * exp(-eta2_2 * (r - betas0_2))
                sigma_22(i) = sigma_22(i) + rtemp
                sigma_22(j) = sigma_22(j) + rtemp

                if (need_grad) then
                    dtheta = (eta2_2 + rtemp1) * rtemp * vec
                    dsigma_22(:,i,i) = dsigma_22(:,i,i) - dtheta
                    dsigma_22(:,j,j) = dsigma_22(:,j,j) + dtheta
                    dsigma_22(:,j,i) = dtheta
                    dsigma_22(:,i,j) = -dtheta
                end if

                ! V_22
                rtemp = theta * exp(-kappadbeta_2 * (r - betas0_2))
                V_22 = V_22 + rtemp

                if (need_grad) then
                    dtheta = (kappadbeta_2 + rtemp1) * rtemp * vec
                    dV_22_2(:,i) = dV_22_2(:,i) + dtheta
                    dV_22_2(:,j) = dV_22_2(:,j) - dtheta
                end if
            end do
        end do

        ! --- Phase 4d: Au-Ad cross pairs ---
        do j = 1, na
            dx = xyz(1) - slab%coords(1,j)
            dy = xyz(2) - slab%coords(2,j)
            call pbc_wrap(dx, dy, slab%box_lx, slab%box_ly, slab%box_skew)
            dz = xyz(3) - slab%coords(3,j)
            r  = sqrt(dx*dx + dy*dy + dz*dz)

            if (r > emt_cutoff_range * rcut) cycle

            if (r > 1.0d-12) then
                vec = [dx, dy, dz] / r
            else
                vec = [0.0d0, 0.0d0, 0.0d0]
            end if
            rtemp  = exp(acut*(r - rcut))
            theta  = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut * rtemp * theta

            ! sigma_21
            rtemp = theta * exp(-eta2_1 * (r - betas0_1))
            sigma_21(j) = sigma_21(j) + rtemp
            if (need_grad) then
                dtheta = (eta2_1 + rtemp1) * rtemp * vec
                dsigma_21_2(:,j) = dsigma_21_2(:,j) + dtheta
                dsigma_21_1(:,j) = -dtheta
            end if

            ! sigma_12
            rtemp = theta * exp(-eta2_2 * (r - betas0_2))
            sigma_12 = sigma_12 + rtemp
            if (need_grad) then
                dtheta = (eta2_2 + rtemp1) * rtemp * vec
                dsigma_12_1(:) = dsigma_12_1(:) - dtheta
                dsigma_12_2(:,j) = dtheta
            end if

            ! V_NN_eff = sum theta(r_j)*NN_pair(r_j), gamma absorbed into NN
            call eval_1d_nn(r, pair_nn_weights, emtnn_r_scale_pair, V_nn, dV_nn_dr)

            V_NN_eff = V_NN_eff + theta * V_nn
            if (need_grad) then
                ! d/dR[theta*NN] = theta*(rtemp1*NN - dNN/dr)*vec
                dtheta_nn = theta * (rtemp1 * V_nn - dV_nn_dr)
                dV_NN_eff_1(:) = dV_NN_eff_1(:) + dtheta_nn * vec(:)
                dV_NN_eff_2(:,j) = dV_NN_eff_2(:,j) - dtheta_nn * vec(:)
            end if

            ! V_rep near-field repulsion (exact zero for r >= rc, so the
            ! well region and all fit-parity geometries are untouched)
            if (emtnn_vrep_on .and. r < emtnn_vrep_rc) then
                rtemp = emtnn_vrep_a * exp(-(r - emtnn_vrep_r0) / emtnn_vrep_rho)
                if (r > emtnn_vrep_r0) then
                    ! cosine taper 1 -> 0 on (r0, rc): w, dw/dr
                    rtemp1 = pi / (emtnn_vrep_rc - emtnn_vrep_r0)
                    vrep_w = 0.5d0 * (1.0d0 + cos(rtemp1 * (r - emtnn_vrep_r0)))
                    dvrep  = rtemp * (-vrep_w / emtnn_vrep_rho &
                              - 0.5d0 * rtemp1 * sin(rtemp1 * (r - emtnn_vrep_r0)))
                else
                    vrep_w = 1.0d0
                    dvrep  = -rtemp / emtnn_vrep_rho
                end if
                V_rep_tot = V_rep_tot + rtemp * vrep_w
                if (need_grad) then
                    dV_rep_1(:)   = dV_rep_1(:) + dvrep * vec(:)
                    dV_rep_2(:,j) = dV_rep_2(:,j) - dvrep * vec(:)
                end if
            end if
        end do

        ! --- Phase 5: gamma scaling ---
        sigma_22 = sigma_22 * igamma1_2
        V_22     = V_22 * igamma2_2 * v0_2
        sigma_21 = sigma_21 * igamma1_2
        sigma_12 = sigma_12 * igamma1_1

        if (need_grad) then
            dsigma_22   = dsigma_22   * igamma1_2
            dsigma_21_2 = dsigma_21_2 * igamma1_2
            dsigma_21_1 = dsigma_21_1 * igamma1_2
            dsigma_12_1 = dsigma_12_1 * igamma1_1
            dsigma_12_2 = dsigma_12_2 * igamma1_1
            dV_22_2     = dV_22_2     * igamma2_2 * v0_2
        end if

        ! --- Phase 6: neutral-sphere radius s ---
        s_1 = max(tolerance, chi_12 * sigma_12)
        do i = 1, na
            s_2(i) = max(tolerance, sigma_22(i) + chi_21 * sigma_21(i))
        end do

        if (need_grad) then
            ds_2_2 = -dsigma_22
            ds_1_1 = 0.0d0

            do i = 1, na
                ds_2_2(:, i, i) = ds_2_2(:, i, i) - chi_21 * dsigma_21_2(:, i)
                do j = 1, na
                    ds_2_2(:, i, j) = ds_2_2(:, i, j) / (betaeta2_2 * s_2(j))
                end do
                ds_1_2(:, i) = -chi_12 * dsigma_12_2(:, i) / (betaeta2_1 * s_1)
            end do

            ds_1_1(:) = ds_1_1(:) - chi_12 * dsigma_12_1(:)
            ds_1_1(:) = ds_1_1(:) / (betaeta2_1 * s_1)
            do j = 1, na
                ds_2_1(:, j) = -chi_21 * dsigma_21_1(:, j) / (betaeta2_2 * s_2(j))
            end do
        end if

        ! s = -ln(sigma_total/12) / (beta*eta2)
        s_1 = -log(s_1 * twelfth) / betaeta2_1
        do i = 1, na
            s_2(i) = -log(s_2(i) * twelfth) / betaeta2_2
        end do

        ! --- Phase 7: cohesive energy ---
        Ecoh_1 = 0.0d0
        Ecoh_2 = 0.0d0

        use_ec = present(ecoh_nn_weights) .and. emtnn_use_ecoh_nn
        if (use_ec) then
            s_ad_ec = -log(max(sigma_12, tolerance) * twelfth) / betaeta2_1
            call eval_1d_nn(s_ad_ec, ecoh_nn_weights, emtnn_ecoh_scale, Ecoh_1, dEcoh_ds)
            exp_lambda_s1 = 0.0d0
        else
            exp_lambda_s1 = exp(-lambda_1 * s_1)
            Ecoh_1 = ((1.0d0 + lambda_1 * s_1) * exp_lambda_s1)
            Ecoh_1 = Ecoh_1 * e0_1
        end if

        do i = 1, na
            exp_lambda_s2(i) = exp(-lambda_2 * s_2(i))
            Ecoh_2 = Ecoh_2 + ((1.0d0 + lambda_2 * s_2(i)) * exp_lambda_s2(i) - 1.0d0)
        end do

        Ecoh_2 = Ecoh_2 * e0_2

        if (need_grad) then
            do i = 1, na
                do j = 1, na
                    dEcoh_2_2(:,i) = dEcoh_2_2(:,i) + s_2(j) * exp_lambda_s2(j) * ds_2_2(:,i,j)
                end do
            end do
            do j = 1, na
                dEcoh_2_1(:) = dEcoh_2_1(:) + s_2(j) * exp_lambda_s2(j) * ds_2_1(:,j)
            end do
            dEcoh_2_2 = dEcoh_2_2 * e0_2 * lambda_2 * lambda_2
            dEcoh_2_1 = dEcoh_2_1 * e0_2 * lambda_2 * lambda_2

            if (use_ec) then
                coef_ec = dEcoh_ds * (-1.0d0/(betaeta2_1 * max(sigma_12, tolerance)))
                dEcoh_1_1(:) = -coef_ec * dsigma_12_1(:)
                do i = 1, na
                    dEcoh_1_2(:,i) = -coef_ec * dsigma_12_2(:,i)
                end do
            else
                do i = 1, na
                    dEcoh_1_2(:,i) = dEcoh_1_2(:,i) + s_1 * exp_lambda_s1 * ds_1_2(:,i)
                end do
                dEcoh_1_1(:) = dEcoh_1_1(:) + s_1 * exp_lambda_s1 * ds_1_1(:)
                dEcoh_1_2 = dEcoh_1_2 * e0_1 * lambda_1 * lambda_1
                dEcoh_1_1 = dEcoh_1_1 * e0_1 * lambda_1 * lambda_1
            end if
        end if

        ! --- Phase 8: reference pair potential ---
        vref_1 = 0.0d0
        vref_2 = 0.0d0

        do i = 1, na
            rtemp = exp(-kappa_2 * s_2(i))
            vref_2 = vref_2 + rtemp
            if (need_grad) then
                do j = 1, na
                    dvref_2_2(:,j) = dvref_2_2(:,j) + rtemp * ds_2_2(:,j,i)
                end do
                dvref_2_1(:) = dvref_2_1(:) + rtemp * ds_2_1(:,i)
            end if
        end do

        rtemp = exp(-kappa_1 * s_1)
        vref_1 = vref_1 + rtemp
        if (need_grad) then
            do j = 1, na
                dvref_1_2(:,j) = dvref_1_2(:,j) + rtemp * ds_1_2(:,j)
            end do
            dvref_1_1(:) = dvref_1_1(:) + rtemp * ds_1_1(:)
        end if

        rtemp     = 12.0d0 * v0_2
        vref_2    = vref_2    * rtemp
        if (need_grad) then
            dvref_2_2 = dvref_2_2 * rtemp * kappa_2
            dvref_2_1 = dvref_2_1 * rtemp * kappa_2
        end if

        rtemp     = 12.0d0 * v0_1
        vref_1    = vref_1    * rtemp
        if (need_grad) then
            dvref_1_2 = dvref_1_2 * rtemp * kappa_1
            dvref_1_1 = dvref_1_1 * rtemp * kappa_1
        end if

        ! --- Phase 9: total energy assembly ---
        E_total = Ecoh_2 + Ecoh_1 - V_22 - V_NN_eff &
                + 0.5d0 * (vref_2 + vref_1) + V_rep_tot
        if (present(E0_nn)) E_total = E_total + E0_nn

        ! --- Phase 10: force assembly ---
        if (need_grad) then
            do i = 1, na
                f_slab_loc(:,i) = dEcoh_2_2(:,i) + dEcoh_1_2(:,i) - dV_22_2(:,i) &
                    - dV_NN_eff_2(:,i) + 0.5d0 * (dvref_2_2(:,i) + dvref_1_2(:,i))
                if (emtnn_vrep_on) f_slab_loc(:,i) = f_slab_loc(:,i) - dV_rep_2(:,i)
            end do
            f_ad(:) = dEcoh_2_1(:) + dEcoh_1_1(:) &
                - dV_NN_eff_1(:) + 0.5d0 * (dvref_2_1(:) + dvref_1_1(:))
            if (emtnn_vrep_on) f_ad(:) = f_ad(:) - dV_rep_1(:)
        end if

        ! --- Phase 11: interaction energy, outputs ---
        E_int = E_total - E_slab_bare_total

        if (need_grad) then
            if (present(dEdx)) dEdx = -f_ad(1)
            if (present(dEdy)) dEdy = -f_ad(2)
            if (present(dEdz)) dEdz = -f_ad(3)
            if (present(f_slab)) then
                if (size(f_slab, 1) >= 3 .and. size(f_slab, 2) >= na) then
                    f_slab(1:3, 1:na) = f_slab_loc(1:3, 1:na)
                end if
            end if
        end if

        if (present(components)) then
            components%E_coh_slab = Ecoh_2 - E_coh_bare_sum
            components%V_ref_slab = 0.5d0 * (vref_2 - V_ref_bare_raw)
            components%E_coh_ad   = Ecoh_1
            components%V_ref_ad   = 0.5d0 * vref_1
            components%V_AuAu     = V_22 - V_AuAu_bare_raw
            components%V_AuAd     = -V_NN_eff + V_rep_tot
        end if

        ! --- deallocate ---
        deallocate(sigma_22, sigma_21, s_2, exp_lambda_s2)
        if (need_grad) then
            deallocate(dsigma_22)
            deallocate(dsigma_21_2, dsigma_21_1, dsigma_12_1, dsigma_12_2)
            deallocate(dV_22_2, dV_NN_eff_2, dV_NN_eff_1)
            if (emtnn_vrep_on) deallocate(dV_rep_2)
            deallocate(ds_2_2, ds_1_1, ds_2_1, ds_1_2)
            deallocate(dEcoh_2_2, dEcoh_1_2, dEcoh_2_1, dEcoh_1_1)
            deallocate(dvref_2_2, dvref_1_2, dvref_2_1, dvref_1_1)
            deallocate(f_ad, f_slab_loc)
        end if
    end subroutine compute_emt_energy

    ! ================================================================
    !  compute_sigma_12 — verbatim from fit pot_emt.f90
    ! ================================================================
    function compute_sigma_12(xyz, eta2_ad, s0_ad) result(sigma_12)
        real(8), intent(in) :: xyz(3)
        real(8), intent(in) :: eta2_ad, s0_ad
        real(8) :: sigma_12

        real(8) :: betas0_1, rnn_1(3), x_1(3), r3temp(3)
        real(8) :: igamma1_1
        real(8) :: dx, dy, dz, r, rtemp, theta
        integer :: j, na

        na = slab%n_atoms
        if (na <= 0) then
            sigma_12 = 0.0d0
            return
        end if

        betas0_1 = beta * s0_ad
        rnn_1(1) = betas0_1
        rnn_1(2) = rnn_1(1) * sqrt2
        rnn_1(3) = rnn_1(1) * sqrt3

        x_1 = [1.0d0, 0.5d0, 2.0d0] / (1.0d0 + exp(acut*(rnn_1 - rcut)))

        r3temp = rnn_1 - betas0_1
        igamma1_1 = 1.0d0 / sum(x_1 * exp(-eta2_ad * r3temp))

        sigma_12 = 0.0d0
        do j = 1, na
            dx = xyz(1) - slab%coords(1,j)
            dy = xyz(2) - slab%coords(2,j)
            call pbc_wrap(dx, dy, slab%box_lx, slab%box_ly, slab%box_skew)
            dz = xyz(3) - slab%coords(3,j)
            r  = sqrt(dx*dx + dy*dy + dz*dz)

            if (r > emt_cutoff_range * rcut) cycle

            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)

            rtemp = theta * exp(-Au%eta2 * (r - betas0_2))
            sigma_12 = sigma_12 + rtemp
        end do

        sigma_12 = sigma_12 * igamma1_1
    end function compute_sigma_12

    ! ================================================================
    !  compute_hubbard_s_ad — verbatim from fit pot_emt.f90
    !  (NO chi_12 factor — matches the fit convention)
    ! ================================================================
    function compute_hubbard_s_ad(xyz, eta2_ad, s0_ad) result(s_ad)
        real(8), intent(in) :: xyz(3)
        real(8), intent(in) :: eta2_ad, s0_ad
        real(8) :: s_ad

        real(8) :: sigma_12

        sigma_12 = compute_sigma_12(xyz, eta2_ad, s0_ad)
        s_ad = -log(max(sigma_12, tolerance) * twelfth) / (beta * eta2_ad)
    end function compute_hubbard_s_ad

    ! ================================================================
    !  compute_hubbard_all — verbatim from fit pot_emt.f90
    ! ================================================================
    subroutine compute_hubbard_all(xyz, ad_params, h, U, J, Vcp)
        real(8), intent(in) :: xyz(3)
        real(8), intent(in) :: ad_params(7)
        real(8), intent(out) :: h, U, J, Vcp

        real(8) :: s_ad, d_dummy
        real(8) :: h_s, U_s, J_s, Vcp_s
        real(8) :: pair_h, pair_U, pair_J, pair_V
        real(8) :: dx, dy, dz, r, rtemp, theta, V_nn_pair
        integer :: islab, na

        if (.not. allocated(hubb_nn_weights)) then
            allocate(hubb_nn_weights(3*emtnn_nh_hubb+1, 4))
            hubb_nn_weights = 0.0d0
        end if

        ! 1. EMT neutral-sphere coordinate
        s_ad = compute_hubbard_s_ad(xyz, ad_params(1), ad_params(7))

        ! 2. s_ad NN for each Hubbard parameter
        call eval_1d_nn(s_ad, hubb_nn_weights(:,1), emtnn_hubb_scale, h_s, d_dummy)
        call eval_1d_nn(s_ad, hubb_nn_weights(:,2), emtnn_hubb_scale, U_s, d_dummy)
        call eval_1d_nn(s_ad, hubb_nn_weights(:,3), emtnn_hubb_scale, J_s, d_dummy)
        call eval_1d_nn(s_ad, hubb_nn_weights(:,4), emtnn_hubb_scale, Vcp_s, d_dummy)

        ! 3. Pairwise NN correction — sum over slab atoms
        pair_h = 0.0d0; pair_U = 0.0d0; pair_J = 0.0d0; pair_V = 0.0d0
        if (allocated(hubb_pair_nn_weights)) then
            na = slab%n_atoms
            do islab = 1, na
                dx = xyz(1) - slab%coords(1,islab)
                dy = xyz(2) - slab%coords(2,islab)
                call pbc_wrap(dx, dy, slab%box_lx, slab%box_ly, slab%box_skew)
                dz = xyz(3) - slab%coords(3,islab)
                r  = sqrt(dx*dx + dy*dy + dz*dz)

                ! Gentle cutoff (acut_soft) — long-range tail matters
                rtemp = exp(acut_soft*(r - rcut))
                theta = 1.0d0 / (1.0d0 + rtemp)

                if (theta > 1.0d-15) then
                    call eval_1d_nn(r, hubb_pair_nn_weights(:,1), &
                        emtnn_r_scale_hubb_pair, V_nn_pair, d_dummy)
                    pair_h = pair_h + theta * V_nn_pair

                    call eval_1d_nn(r, hubb_pair_nn_weights(:,2), &
                        emtnn_r_scale_hubb_pair, V_nn_pair, d_dummy)
                    pair_U = pair_U + theta * V_nn_pair

                    call eval_1d_nn(r, hubb_pair_nn_weights(:,3), &
                        emtnn_r_scale_hubb_pair, V_nn_pair, d_dummy)
                    pair_J = pair_J + theta * V_nn_pair

                    call eval_1d_nn(r, hubb_pair_nn_weights(:,4), &
                        emtnn_r_scale_hubb_pair, V_nn_pair, d_dummy)
                    pair_V = pair_V + theta * V_nn_pair
                end if
            end do
        end if

        ! 4. Combine s_ad NN + pairwise NN
        h   = h_s   + pair_h
        U   = U_s   + pair_U
        J   = J_s   + pair_J
        Vcp = Vcp_s + pair_V

        ! 5. Cache results
        emt_cached_h    = h
        emt_cached_U    = U
        emt_cached_J    = J
        emt_cached_V    = Vcp
        emt_cached_s_ad = s_ad
    end subroutine compute_hubbard_all

    ! ================================================================
    !  load_slab_xyz — verbatim from fit pot_emt.f90
    !  (layer_id classification hardcodes 2.95 — no physical effect,
    !   layer_id has no consumers in the engine)
    ! ================================================================
    subroutine load_slab_xyz(filename, slab_out, a_lat, n_side)
        character(len=*),intent(in)::filename; type(rst_slab_data),intent(out)::slab_out
        real(8),intent(in),optional::a_lat; integer,intent(in),optional::n_side
        integer::unit,ios,n,i; real(8)::x,y,z,z_max,d_layer_est; character(len=2)::elem; character(len=256)::sline
        open(newunit=unit,file=trim(filename),status='old',iostat=ios)
        if(ios/=0)then; write(*,*) "ERROR: Cannot open slab file: ",trim(filename); stop; end if
        read(unit,*,iostat=ios)n; read(unit,'(A)',iostat=ios)sline
        slab_out%n_atoms=n
        if(allocated(slab_out%coords))deallocate(slab_out%coords)
        if(allocated(slab_out%layer_id))deallocate(slab_out%layer_id)
        allocate(slab_out%coords(3,n)); allocate(slab_out%layer_id(n))
        do i=1,n; read(unit,*,iostat=ios)elem,x,y,z; slab_out%coords(1,i)=x; slab_out%coords(2,i)=y; slab_out%coords(3,i)=z; end do
        close(unit)
        if(present(a_lat).and.present(n_side))then; slab_out%box_lx=n_side*a_lat; slab_out%box_ly=n_side*a_lat; slab_out%box_skew=PI/3.0d0; end if
        z_max=maxval(slab_out%coords(3,:)); d_layer_est=2.95d0*sr2o3; slab_out%n_surf=0
        do i=1,n
            if(slab_out%coords(3,i)>z_max-0.5d0*d_layer_est)then; slab_out%layer_id(i)=1; slab_out%n_surf=slab_out%n_surf+1
            else; slab_out%layer_id(i)=2; end if
        end do
        write(*,"(A,A,A,I0,A,I0)")"  Slab loaded: ",trim(filename),": ",slab_out%n_atoms," atoms, ",slab_out%n_surf," surface"
    end subroutine load_slab_xyz

    ! ================================================================
    !  init_slab_baseline — fit emt_init Phase B1-B7 (verbatim),
    !  always loading the slab from file. The bare-slab energy is
    !  FROZEN here (movable-slab dynamics subtract the initial
    !  baseline, i.e. E_int is relative to the pristine slab).
    ! ================================================================
    subroutine init_slab_baseline()
        real(8) :: r3temp(3)
        real(8) :: rr
        real(8) :: rtemp, theta
        real(8), allocatable :: sigma_22(:), s_2(:)
        real(8) :: V_22, Ecoh_2, vref_2
        real(8) :: dx, dy, dz, r
        integer :: i, j, na

        write(*,"(A)") "EMT: Initializing bare slab baseline (fit_emt_nn convention)"

        call load_slab_xyz(slab_file_val, slab, a_lat_val, n_side_val)
        na = slab%n_atoms
        write(*,*) "EMT slab initialized: ", na, " atoms"

        ! Phase B1 — FCC geometry & Au derived quantities
        betas0_2     = beta * Au%s0
        betaeta2_2   = beta * Au%eta2
        kappadbeta_2 = Au%kappa / beta

        rnn_2(1) = betas0_2
        rnn_2(2) = rnn_2(1) * sqrt2
        rnn_2(3) = rnn_2(1) * sqrt3

        rcut = betas0_2 * sqrt3
        rr          = 4.0d0 * rcut / (sqrt3 + 2.0d0)
        acut = log(10.0d0) / (rr - rcut)

        x_2 = [1.0d0, 0.5d0, 2.0d0] / (1.0d0 + exp(acut*(rnn_2 - rcut)))

        r3temp = rnn_2 - betas0_2
        igamma1_2 = 1.0d0 / sum(x_2 * exp(-Au%eta2 * r3temp))
        igamma2_2 = 1.0d0 / sum(x_2 * exp(-kappadbeta_2 * r3temp))

        write(*,"(A,F8.4,A,F10.6)")      "EMT: rcut = ", rcut, "  acut = ", acut
        write(*,"(A,F8.6,A,F8.6)")       "EMT: igamma1_2 = ", igamma1_2, "  igamma2_2 = ", igamma2_2

        ! Phase B2 — Au-Au sigma and pair accumulation
        allocate(sigma_22(na), s_2(na))
        sigma_22 = 0.0d0
        V_22     = 0.0d0

        do i = 1, na - 1
            do j = i + 1, na
                dx = slab%coords(1,i) - slab%coords(1,j)
                dy = slab%coords(2,i) - slab%coords(2,j)
                call pbc_wrap(dx, dy, slab%box_lx, slab%box_ly, slab%box_skew)
                dz = slab%coords(3,i) - slab%coords(3,j)
                r  = sqrt(dx*dx + dy*dy + dz*dz)

                if (r > emt_cutoff_range * rcut) cycle

                rtemp  = exp(acut*(r - rcut))
                theta  = 1.0d0 / (1.0d0 + rtemp)

                rtemp = theta * exp(-Au%eta2 * (r - betas0_2))
                sigma_22(i) = sigma_22(i) + rtemp
                sigma_22(j) = sigma_22(j) + rtemp

                rtemp = theta * exp(-kappadbeta_2 * (r - betas0_2))
                V_22 = V_22 + rtemp
            end do
        end do

        ! Phase B3 — gamma scaling
        sigma_22 = sigma_22 * igamma1_2
        V_22     = V_22 * igamma2_2 * Au%v0

        ! Phase B4 — neutral-sphere radius
        do i = 1, na
            s_2(i) = max(tolerance, sigma_22(i))
            s_2(i) = -log(s_2(i) * twelfth) / betaeta2_2
        end do

        ! Phase B5 — cohesive energy
        Ecoh_2 = 0.0d0
        do i = 1, na
            Ecoh_2 = Ecoh_2 + ((1.0d0 + Au%lambda * s_2(i)) &
                          * exp(-Au%lambda * s_2(i)) - 1.0d0)
        end do
        Ecoh_2 = Ecoh_2 * Au%e0

        ! Phase B6 — reference pair potential
        vref_2 = 0.0d0
        do i = 1, na
            rtemp = exp(-Au%kappa * s_2(i))
            vref_2 = vref_2 + rtemp
        end do
        vref_2 = vref_2 * 12.0d0 * Au%v0

        ! Phase B7 — bare slab total energy
        E_slab_bare_total = Ecoh_2 - V_22 + 0.5d0 * vref_2

        E_coh_bare_sum  = Ecoh_2
        V_ref_bare_raw  = vref_2
        V_AuAu_bare_raw = V_22

        deallocate(sigma_22, s_2)

        write(*,"(A,F12.6)") "EMT: E_slab_bare_total = ", E_slab_bare_total
        write(*,"(A,F12.6)") "EMT:   E_coh_2         = ", Ecoh_2
        write(*,"(A,F12.6)") "EMT:   V_22            = ", V_22
        write(*,"(A,F12.6)") "EMT:   0.5*vref_2      = ", 0.5d0*vref_2
    end subroutine init_slab_baseline

    ! ================================================================
    !  Weight loading — adapted from fit pot_emt.f90 emt_load_weights
    !
    !  KEY FIX vs the fit version: the fit loader ends its 'NH_pair:'
    !  branch with 'cycle', so when the production header has all three
    !  keys on ONE line ("# NH_pair: 8 NH_hubb: 8 NH_hubb_pair: 3")
    !  NH_hubb / NH_hubb_pair are never read from the file — the fit
    !  relies on fit_config.txt to preset them. There is no config in
    !  dynamics, so this version scans EVERY key on the same line.
    !  Missing file / missing GS / missing Hubbard section => STOP.
    ! ================================================================
    subroutine emt_load_weights(filename)
        character(len=*),intent(in)::filename
        character(200)::line; integer::unit,ios,idx,i
        integer :: n_nn, nread, n_ec
        real(8), allocatable :: tmpw(:)
        logical :: found_gs
        n_nn = 3*emtnn_nh_pair + 1
        n_ec = 0
        if (emtnn_use_ecoh_nn) n_ec = 3*emtnn_nh_ecoh + 1
        have_fit_weights = .false.
        ecoh_weights_loaded = .false.
        hubbard_section_seen = .false.
        if (allocated(fit_weights)) deallocate(fit_weights)
        allocate(fit_weights(7+n_nn+1+3*emtnn_nh_ecoh+1+7)); fit_weights = 0.0d0
        found_gs = .false.
        open(newunit=unit,file=filename,status='old',iostat=ios)
        if(ios/=0)then
            write(*,*) "ERROR: EMT-NN weight file not found: ",trim(filename)
            write(*,*) "       The dynamics engine refuses to run with a zero potential."
            stop
        end if
        do
            read(unit,'(A)',iostat=ios)line; if(ios/=0)exit; line=adjustl(line)
            if(len_trim(line)==0)cycle
            ! Parse NH info — NO early cycle: all keys may share one line
            idx=index(line,'NH_pair:')
            if(idx>0)then
                read(line(idx+8:),*,iostat=ios) emtnn_nh_pair
                if(ios==0) then
                    n_nn = 3*emtnn_nh_pair + 1
                    deallocate(fit_weights)
                    allocate(fit_weights(7+n_nn+1+3*emtnn_nh_ecoh+1+7)); fit_weights = 0.0d0
                end if
            end if
            idx=index(line,'NH_hubb_pair:')
            if(idx>0)then
                read(line(idx+13:),*,iostat=ios) emtnn_nh_hubb_pair
                if (ios==0) then
                    if (allocated(hubb_pair_nn_weights)) deallocate(hubb_pair_nn_weights)
                    allocate(hubb_pair_nn_weights(3*emtnn_nh_hubb_pair+1, 4))
                    hubb_pair_nn_weights = 0.0d0
                end if
            end if
            idx=index(line,'NH_hubb:')
            if(idx>0)then
                read(line(idx+8:),*,iostat=ios) emtnn_nh_hubb
                if (ios==0) then
                    if (allocated(hubb_nn_weights)) deallocate(hubb_nn_weights)
                    allocate(hubb_nn_weights(3*emtnn_nh_hubb+1, 4))
                    hubb_nn_weights = 0.0d0
                end if
            end if
            idx=index(line,'NH_ecoh:')
            if(idx>0)then
                read(line(idx+8:),*,iostat=ios) emtnn_nh_ecoh
                if (ios==0) then
                    n_ec = 3*emtnn_nh_ecoh + 1
                    if (size(fit_weights) < 7+n_nn+1+n_ec+7) then
                        allocate(tmpw(7+n_nn+1+n_ec+7)); tmpw = 0.0d0
                        tmpw(1:size(fit_weights)) = fit_weights
                        call move_alloc(tmpw, fit_weights)
                    end if
                end if
            end if
            ! Parse state weights — GS only
            if(index(line,'# STATE:')==1)then
                read(unit,'(A)',iostat=ios)line; if(ios/=0)exit
                line=adjustl(line)
                if(index(line,'## EMT_params')>0 .or. index(line,'# EMT_params')>0)then
                    read(unit,'(A)',iostat=ios)line
                    read(line,*,iostat=ios) fit_weights(1:7)
                    read(unit,'(A)',iostat=ios)line  ! ## Pair_NN
                    nread = 0
                    do while (nread < n_nn)
                        read(unit,'(A)',iostat=ios)line
                        if (ios /= 0) exit
                        line = adjustl(line)
                        if (len_trim(line) == 0) exit
                        if (line(1:1) == '#') then
                            backspace(unit); exit
                        end if
                        do while (len_trim(line) > 0 .and. nread < n_nn)
                            idx = index(trim(line), ' ')
                            if (idx == 0) idx = len_trim(line) + 1
                            nread = nread + 1
                            read(line(1:idx-1), *, iostat=ios) fit_weights(7+nread)
                            line = adjustl(line(idx:))
                        end do
                    end do
                    if (nread < n_nn) then
                        write(*,*) "ERROR: Pair_NN section truncated: ",nread," of ",n_nn
                        stop
                    end if
                    read(unit,'(A)',iostat=ios)line  ! # E0
                    read(unit,'(A)',iostat=ios)line  ! E0 value
                    read(line,*,iostat=ios) fit_weights(8+n_nn)
                    found_gs = .true.
                    write(*,"(A,7F9.4)") "  Loaded EMT-NN weights: ", fit_weights(1:7)
                else
                    ! Legacy 7-param format — rejected (obsolete architecture)
                    write(*,*) "ERROR: legacy 7-param weight format detected — obsolete."
                    write(*,*) "       Provide a fit_emt_nn GS weight file (## EMT_params + Pair_NN)."
                    stop
                end if
                cycle
            end if
            ! Ecoh NN weights (engine keeps the analytic Ecoh path)
            if (index(line, '## Ecoh_NN') > 0) then
                n_ec = 3*emtnn_nh_ecoh + 1
                if (size(fit_weights) < 8+n_nn+n_ec) then
                    allocate(tmpw(7+n_nn+1+n_ec+7)); tmpw = 0.0d0
                    tmpw(1:size(fit_weights)) = fit_weights
                    call move_alloc(tmpw, fit_weights)
                end if
                nread = 0
                do while (nread < n_ec)
                    read(unit,'(A)',iostat=ios)line
                    if (ios /= 0) exit
                    line = adjustl(line)
                    if (len_trim(line) == 0) exit
                    if (line(1:1) == '#') then
                        backspace(unit); exit
                    end if
                    do while (len_trim(line) > 0 .and. nread < n_ec)
                        idx = index(trim(line), ' ')
                        if (idx == 0) idx = len_trim(line) + 1
                        nread = nread + 1
                        read(line(1:idx-1), *, iostat=ios) fit_weights(8+n_nn+nread)
                        line = adjustl(line(idx:))
                    end do
                end do
                if (nread < n_ec) then
                    write(*,*) "ERROR: Ecoh_NN section truncated:", nread, "of", n_ec
                    stop
                end if
                ecoh_weights_loaded = .true.
                write(*,"(A,I0,A)") "  Loaded Ecoh NN weights (", n_ec, &
                    " values) — will be activated at engine init"
                cycle
            end if
            ! Hubbard NN weights
            if (index(line, '# HUBBARD_NN') > 0) then
                hubbard_section_seen = .true.
                call load_hubbard_nn_section(unit); cycle
            end if
        end do; close(unit)
        if (found_gs) have_fit_weights = .true.
        if (.not. have_fit_weights) then
            write(*,*) "ERROR: no '# STATE: GS' section in ",trim(filename)
            stop
        end if
    end subroutine emt_load_weights

    ! ================================================================
    !  load_hubbard_nn_section — fit version + cross-check of the
    !  in-section '# nh_hubb:' line against the header value.
    ! ================================================================
    subroutine load_hubbard_nn_section(unit)
        integer, intent(in) :: unit
        character(200) :: line
        integer :: ios, ip, nw, nw_pair, nread, idx, nh_check
        nw = 3*emtnn_nh_hubb + 1
        if (.not. allocated(hubb_nn_weights)) then
            allocate(hubb_nn_weights(nw, 4))
        else if (size(hubb_nn_weights,1) /= nw) then
            deallocate(hubb_nn_weights)
            allocate(hubb_nn_weights(nw, 4))
        end if
        hubb_nn_weights = 0.0d0

        ! In-section "# nh_hubb: N" — cross-check against the header value
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) return
        line = adjustl(line)
        idx = index(line, 'nh_hubb:')
        if (idx > 0) then
            read(line(idx+8:), *, iostat=ios) nh_check
            if (ios == 0 .and. nh_check /= emtnn_nh_hubb) then
                write(*,*) "ERROR: HUBBARD_NN section nh_hubb = ",nh_check, &
                           " but header says ",emtnn_nh_hubb
                write(*,*) "       (The fit loader skips this line; the dynamics", &
                           " loader cross-checks it.)"
                stop
            end if
        else
            backspace(unit)
        end if

        ! --- Read s_ad NN weights (4 blocks) ---
        do ip = 1, 4
            read(unit, '(A)', iostat=ios) line  ! skip comment header
            if (ios /= 0) return
            nread = 0
            do while (nread < nw)
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(line)
                if (len_trim(line) == 0) exit
                if (line(1:1) == '#' .and. nread > 0) then
                    backspace(unit); exit
                end if
                do while (len_trim(line) > 0 .and. nread < nw)
                    idx = index(trim(line), ' ')
                    if (idx == 0) idx = len_trim(line) + 1
                    nread = nread + 1
                    read(line(1:idx-1), *, iostat=ios) hubb_nn_weights(nread, ip)
                    line = adjustl(line(idx:))
                end do
            end do
            if (nread < nw) then
                write(*,*) "ERROR: Hubbard s_ad NN block ",ip," truncated: ",nread," of ",nw
                stop
            end if
        end do

        ! --- Try to read pairwise NN weights (4 blocks, optional) ---
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) then
            write(*,*) "  Loaded Hubbard s_ad NN weights"
            return
        end if
        line = adjustl(line)
        if (line(1:1) == '#' .and. index(line, 'pairwise') > 0) then
            nw_pair = 3*emtnn_nh_hubb_pair + 1
            if (.not. allocated(hubb_pair_nn_weights)) then
                allocate(hubb_pair_nn_weights(nw_pair, 4))
            else if (size(hubb_pair_nn_weights,1) /= nw_pair) then
                deallocate(hubb_pair_nn_weights)
                allocate(hubb_pair_nn_weights(nw_pair, 4))
            end if
            hubb_pair_nn_weights = 0.0d0

            do ip = 1, 4
                if (ip > 1) then
                    read(unit, '(A)', iostat=ios) line  ! skip comment
                    if (ios /= 0) return
                end if
                nread = 0
                do while (nread < nw_pair)
                    read(unit, '(A)', iostat=ios) line
                    if (ios /= 0) exit
                    line = adjustl(line)
                    if (len_trim(line) == 0) exit
                    if (line(1:1) == '#' .and. nread > 0) then
                        backspace(unit); exit
                    end if
                    do while (len_trim(line) > 0 .and. nread < nw_pair)
                        idx = index(trim(line), ' ')
                        if (idx == 0) idx = len_trim(line) + 1
                        nread = nread + 1
                        read(line(1:idx-1), *, iostat=ios) hubb_pair_nn_weights(nread, ip)
                        line = adjustl(line(idx:))
                    end do
                end do
                if (nread < nw_pair) then
                    write(*,*) "ERROR: Hubbard pairwise NN block ",ip," truncated: ",nread," of ",nw_pair
                    stop
                end if
            end do
            write(*,*) "  Loaded Hubbard NN weights (s_ad + pairwise)"
        else
            write(*,*) "  Loaded Hubbard s_ad NN weights"
        end if
    end subroutine load_hubbard_nn_section

    ! ================================================================
    !  Public API — initialization
    ! ================================================================
    subroutine emt_nn_init(slab_file, weights_file, a_lat, n_side)
        character(len=*), intent(in) :: slab_file, weights_file
        real(8), intent(in) :: a_lat
        integer, intent(in) :: n_side
        integer :: n_nn

        slab_file_val = slab_file
        a_lat_val  = a_lat
        n_side_val = n_side

        call init_slab_baseline()
        call emt_load_weights(weights_file)

        if (.not. hubbard_section_seen) then
            write(*,*) "ERROR: no '# HUBBARD_NN' section in ",trim(weights_file)
            write(*,*) "       h/U/J/Vcp come from the Hubbard NN — the dynamics"
            write(*,*) "       cannot run with silent zero electronic structure."
            stop
        end if
        hubbard_loaded = .true.

        ! Extract GS parameters once (fit emt_eval layout)
        n_nn = 3*emtnn_nh_pair + 1
        if (size(fit_weights) < 8+n_nn) then
            write(*,*) "ERROR: fit_weights shorter than 7+NN+1"
            stop
        end if
        ad7_w(1:7) = fit_weights(1:7)
        if (allocated(pair_w)) deallocate(pair_w)
        allocate(pair_w(n_nn))
        pair_w(1:n_nn) = fit_weights(8:7+n_nn)
        E0_w = fit_weights(8+n_nn)

        ! Ecoh NN: extracted only when the weight file carried a ## Ecoh_NN
        ! section (ecoh_weights_loaded set by emt_load_weights). Without it
        ! the engine keeps the analytic Ecoh path (emtnn_use_ecoh_nn = F).
        if (ecoh_weights_loaded) then
            if (allocated(ecoh_w)) deallocate(ecoh_w)
            allocate(ecoh_w(3*emtnn_nh_ecoh + 1))
            ecoh_w = fit_weights(8+n_nn+1 : 8+n_nn+3*emtnn_nh_ecoh+1)
            emtnn_use_ecoh_nn = .true.
            write(*,"(A,I0,A)") "EMT-NN: Ecoh NN ACTIVE (", 3*emtnn_nh_ecoh+1, &
                " weights) — replaces analytic Ecoh(adsorbate) path"
        end if

        engine_loaded = .true.

        write(*,"(A,I0,A,I0,A,I0,A,I0)") &
            "EMT-NN: NH_pair = ",emtnn_nh_pair, &
            "  NH_hubb = ",emtnn_nh_hubb, &
            "  NH_hubb_pair = ",emtnn_nh_hubb_pair, &
            "  NH_ecoh = ",emtnn_nh_ecoh
        write(*,"(A,F8.4,A,F10.6,A)") "EMT-NN engine ready: rcut = ",rcut, &
            "  acut = ",acut,"  (ln10 family ~6.1 — must NOT be ~25)"
        write(*,"(A,ES14.6)") "EMT-NN: E0_w (NN offset) = ", E0_w
        write(*,"(A)") "EMT-NN: energies E_int[eV] = E_total - bare slab (frozen baseline)"
    end subroutine emt_nn_init

    ! ================================================================
    !  Public API — total energy + forces
    !
    !  Q_emt(3,natoms): column 1 = adsorbate, columns 2..natoms = slab.
    !  The slab columns are synced into the engine slab every call
    !  (movable slab). forces(3,natoms) in eV/A, f = -dE/dR.
    ! ================================================================
    subroutine calc_emt_nn_energy(natoms, natoma, Q_emt, E, forces, components)
        integer, intent(in) :: natoms, natoma
        real(8), intent(in) :: Q_emt(3,natoms)
        real(8), intent(out) :: E
        real(8), intent(out), optional :: forces(3,natoms)
        type(emt_component), intent(out), optional :: components

        real(8) :: dEdx, dEdy, dEdz
        integer :: k, na

        if (.not. engine_loaded) then
            write(*,*) "ERROR calc_emt_nn_energy: engine not initialized"; stop
        end if
        if (natoma /= 1) then
            write(*,*) "ERROR calc_emt_nn_energy: kernel supports a single adsorbate,", &
                       " natoma = ",natoma; stop
        end if
        na = slab%n_atoms
        if (natoms /= 1 + na) then
            write(*,*) "ERROR calc_emt_nn_energy: natoms = ",natoms, &
                       " but engine slab has ",na," atoms"; stop
        end if

        ! Sync movable slab coordinates into the engine slab
        do k = 1, na
            slab%coords(1,k) = Q_emt(1,k+1)
            slab%coords(2,k) = Q_emt(2,k+1)
            slab%coords(3,k) = Q_emt(3,k+1)
        end do

        if (present(forces)) then
            if (allocated(ecoh_w)) then
                call compute_emt_energy(Q_emt(1:3,1), ad7_w, E, &
                                        dEdx=dEdx, dEdy=dEdy, dEdz=dEdz, &
                                        f_slab=forces(1:3,2:natoms), &
                                        components=components, &
                                        pair_nn_weights=pair_w, E0_nn=E0_w, &
                                        ecoh_nn_weights=ecoh_w)
            else
                call compute_emt_energy(Q_emt(1:3,1), ad7_w, E, &
                                        dEdx=dEdx, dEdy=dEdy, dEdz=dEdz, &
                                        f_slab=forces(1:3,2:natoms), &
                                        components=components, &
                                        pair_nn_weights=pair_w, E0_nn=E0_w)
            end if
            forces(1,1) = -dEdx
            forces(2,1) = -dEdy
            forces(3,1) = -dEdz
        else
            if (allocated(ecoh_w)) then
                call compute_emt_energy(Q_emt(1:3,1), ad7_w, E, &
                                        components=components, &
                                        pair_nn_weights=pair_w, E0_nn=E0_w, &
                                        ecoh_nn_weights=ecoh_w)
            else
                call compute_emt_energy(Q_emt(1:3,1), ad7_w, E, &
                                        components=components, &
                                        pair_nn_weights=pair_w, E0_nn=E0_w)
            end if
        end if
    end subroutine calc_emt_nn_energy

    ! ================================================================
    !  Public API — Hubbard parameters at the current geometry
    ! ================================================================
    subroutine calc_hubbard(Q_imp, Q_slab, n_slab, h, U, J, Vcp)
        real(8), intent(in) :: Q_imp(3)
        real(8), intent(in) :: Q_slab(3,n_slab)
        integer, intent(in) :: n_slab
        real(8), intent(out) :: h, U, J, Vcp
        integer :: k

        if (.not. hubbard_loaded) then
            write(*,*) "ERROR calc_hubbard: Hubbard NN not loaded"; stop
        end if
        if (n_slab /= slab%n_atoms) then
            write(*,*) "ERROR calc_hubbard: n_slab = ",n_slab, &
                       " but engine slab has ",slab%n_atoms; stop
        end if

        do k = 1, n_slab
            slab%coords(1,k) = Q_slab(1,k)
            slab%coords(2,k) = Q_slab(2,k)
            slab%coords(3,k) = Q_slab(3,k)
        end do

        call compute_hubbard_all(Q_imp, ad7_w, h, U, J, Vcp)
    end subroutine calc_hubbard

    ! ================================================================
    !  fill_hubbard_nn_grad — ANALYTICAL gradient of the Hubbard-NN
    !  parameters (new code; the fit has no gradient path for these).
    !
    !    param_k = NN_k(s_ad) + sum_j theta_soft(r_j)*NNpair_k(r_j)
    !
    !  Chain rule (fit conventions, pot_emt.f90:894-906 + 1057-1130):
    !    s_ad  = -ln(max(sigma12,tol)/12)/(beta*eta2_ad)
    !          => ds_ad/dsigma = -1/(beta*eta2_ad*sigma)
    !    sigma12 = igamma1_1 * sum_j theta(r_j)*exp(-Au%eta2*(r_j-betas0_2))
    !          => d/dr(theta*e) = -theta*e*(Au%eta2 + acut*(1-theta))
    !             with the SAME r > cutoff_range*rcut skip as the value path
    !    theta_soft' = -acut_soft*theta_soft*(1-theta_soft)
    !    f'(r) = theta_soft'*NN + theta_soft*dNN/dr   (pairwise term)
    !
    !  grad layout (pes_emt_mod fill_hubbard_param_grad style):
    !    grad(1:3)              = d(param)/d(Q_imp)
    !    grad(3+3*(j-1)+k)      = d(param)/d(Q_slab(j,k))
    !  ipoly: 1=h, 2=U, 3=J, 4=Vcp
    ! ================================================================
    subroutine fill_hubbard_nn_grad(grad, ncoord, ipoly, Q_imp, Q_slab, n_slab, natoma)
        real(8), intent(out) :: grad(ncoord)
        integer, intent(in) :: ncoord, ipoly, n_slab, natoma
        real(8), intent(in) :: Q_imp(3)
        real(8), intent(in) :: Q_slab(3,n_slab)

        real(8) :: eta2_ad, s0_ad, betas0_1_l, rnn_1(3), x_1(3), r3temp(3), igamma1_1_l
        real(8) :: dx, dy, dz, r, rtemp, theta, eth, dsig_dr
        real(8) :: sigma_sum, sigma_12_loc, s_ad_loc, Vnn_s, dNN_ds, coef_sig
        real(8) :: theta_s, dtheta_s, V_pair, dV_pair, fpair
        real(8) :: vec(3)
        real(8), allocatable :: gsig(:)
        integer :: j, na

        grad = 0.0d0
        if (.not. hubbard_loaded) then
            write(*,*) "ERROR fill_hubbard_nn_grad: Hubbard NN not loaded"; return
        end if
        if (natoma /= 1) then
            write(*,*) "ERROR fill_hubbard_nn_grad: natoma must be 1"; stop
        end if
        if (ncoord /= 3*(1+n_slab)) then
            write(*,*) "ERROR fill_hubbard_nn_grad: ncoord = ",ncoord, &
                       " expected ",3*(1+n_slab); stop
        end if
        if (ipoly < 1 .or. ipoly > 4) then
            write(*,*) "ERROR fill_hubbard_nn_grad: ipoly = ",ipoly; stop
        end if
        na = n_slab
        if (na /= slab%n_atoms) then
            write(*,*) "ERROR fill_hubbard_nn_grad: slab size mismatch"; stop
        end if

        ! Sync slab coordinates (value path uses the same geometry)
        do j = 1, na
            slab%coords(1,j) = Q_slab(1,j)
            slab%coords(2,j) = Q_slab(2,j)
            slab%coords(3,j) = Q_slab(3,j)
        end do

        ! Adsorbate-dependent gamma (same as compute_sigma_12)
        eta2_ad = ad7_w(1)
        s0_ad   = ad7_w(7)
        betas0_1_l = beta * s0_ad
        rnn_1(1) = betas0_1_l
        rnn_1(2) = rnn_1(1) * sqrt2
        rnn_1(3) = rnn_1(1) * sqrt3
        x_1 = [1.0d0, 0.5d0, 2.0d0] / (1.0d0 + exp(acut*(rnn_1 - rcut)))
        r3temp = rnn_1 - betas0_1_l
        igamma1_1_l = 1.0d0 / sum(x_1 * exp(-eta2_ad * r3temp))

        allocate(gsig(ncoord))
        gsig = 0.0d0
        sigma_sum = 0.0d0

        do j = 1, na
            dx = Q_imp(1) - slab%coords(1,j)
            dy = Q_imp(2) - slab%coords(2,j)
            call pbc_wrap(dx, dy, slab%box_lx, slab%box_ly, slab%box_skew)
            dz = Q_imp(3) - slab%coords(3,j)
            r  = sqrt(dx*dx + dy*dy + dz*dz)
            if (r > 1.0d-12) then
                vec = [dx, dy, dz] / r
            else
                vec = [0.0d0, 0.0d0, 0.0d0]
            end if

            ! --- pairwise NN term (acut_soft cutoff, no skip, guard) ---
            rtemp = exp(acut_soft*(r - rcut))
            theta_s = 1.0d0 / (1.0d0 + rtemp)
            if (theta_s > 1.0d-15 .and. allocated(hubb_pair_nn_weights)) then
                call eval_1d_nn(r, hubb_pair_nn_weights(:,ipoly), &
                                emtnn_r_scale_hubb_pair, V_pair, dV_pair)
                dtheta_s = -acut_soft * theta_s * (1.0d0 - theta_s)
                fpair = dtheta_s * V_pair + theta_s * dV_pair
                grad(1:3) = grad(1:3) + fpair * vec
                grad(3+3*(j-1)+1) = grad(3+3*(j-1)+1) - fpair * vec(1)
                grad(3+3*(j-1)+2) = grad(3+3*(j-1)+2) - fpair * vec(2)
                grad(3+3*(j-1)+3) = grad(3+3*(j-1)+3) - fpair * vec(3)
            end if

            ! --- sigma_12 term (hard skip, matches value path) ---
            if (r > emt_cutoff_range * rcut) cycle
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            eth = theta * exp(-Au%eta2 * (r - betas0_2))
            sigma_sum = sigma_sum + eth
            dsig_dr = -eth * (Au%eta2 + acut * (1.0d0 - theta))
            gsig(1:3) = gsig(1:3) + dsig_dr * vec
            gsig(3+3*(j-1)+1) = gsig(3+3*(j-1)+1) - dsig_dr * vec(1)
            gsig(3+3*(j-1)+2) = gsig(3+3*(j-1)+2) - dsig_dr * vec(2)
            gsig(3+3*(j-1)+3) = gsig(3+3*(j-1)+3) - dsig_dr * vec(3)
        end do

        ! --- chain rule through s_ad into the s_ad-NN ---
        sigma_12_loc = sigma_sum * igamma1_1_l
        s_ad_loc = -log(max(sigma_12_loc, tolerance) * twelfth) / (beta * eta2_ad)
        call eval_1d_nn(s_ad_loc, hubb_nn_weights(:,ipoly), emtnn_hubb_scale, &
                        Vnn_s, dNN_ds)
        coef_sig = igamma1_1_l * dNN_ds * &
                   (-1.0d0 / (beta * eta2_ad * max(sigma_12_loc, tolerance)))
        grad = grad + coef_sig * gsig

        deallocate(gsig)
    end subroutine fill_hubbard_nn_grad

    ! ================================================================
    !  calc_hubbard_grad_fd — FD fallback for fill_hubbard_nn_grad
    !  (same signature; flip use_hubb_grad_fd = .true. to use it)
    ! ================================================================
    subroutine calc_hubbard_grad_fd(grad, ncoord, ipoly, Q_imp, Q_slab, n_slab, natoma)
        real(8), intent(out) :: grad(ncoord)
        integer, intent(in) :: ncoord, ipoly, n_slab, natoma
        real(8), intent(in) :: Q_imp(3)
        real(8), intent(in) :: Q_slab(3,n_slab)

        real(8), parameter :: dx = 1.0d-4
        real(8) :: Qp(3), Sp(3,n_slab), val_p, val_m, hd, Ud, Jd, Vd
        integer :: i, k, c

        grad = 0.0d0
        if (ncoord /= 3*(1+n_slab)) then
            write(*,*) "ERROR calc_hubbard_grad_fd: ncoord mismatch"; stop
        end if

        ! impurity components
        do k = 1, 3
            Qp = Q_imp; Qp(k) = Q_imp(k) + dx
            call calc_hubbard(Qp, Q_slab, n_slab, hd, Ud, Jd, Vd)
            val_p = pick(ipoly, hd, Ud, Jd, Vd)
            Qp = Q_imp; Qp(k) = Q_imp(k) - dx
            call calc_hubbard(Qp, Q_slab, n_slab, hd, Ud, Jd, Vd)
            val_m = pick(ipoly, hd, Ud, Jd, Vd)
            grad(k) = (val_p - val_m) / (2.0d0 * dx)
        end do

        ! slab components
        do i = 1, n_slab
            do k = 1, 3
                Sp = Q_slab
                Sp(k,i) = Q_slab(k,i) + dx
                call calc_hubbard(Q_imp, Sp, n_slab, hd, Ud, Jd, Vd)
                val_p = pick(ipoly, hd, Ud, Jd, Vd)
                Sp = Q_slab
                Sp(k,i) = Q_slab(k,i) - dx
                call calc_hubbard(Q_imp, Sp, n_slab, hd, Ud, Jd, Vd)
                val_m = pick(ipoly, hd, Ud, Jd, Vd)
                c = 3 + 3*(i-1) + k
                grad(c) = (val_p - val_m) / (2.0d0 * dx)
            end do
        end do

    contains
        function pick(ii, a, b, c3, d) result(v)
            integer, intent(in) :: ii
            real(8), intent(in) :: a, b, c3, d
            real(8) :: v
            select case(ii)
            case(1); v = a
            case(2); v = b
            case(3); v = c3
            case default; v = d
            end select
        end function pick
    end subroutine calc_hubbard_grad_fd

    ! ================================================================
    !  Public queries
    ! ================================================================
    logical function emt_nn_loaded()
        emt_nn_loaded = engine_loaded
    end function emt_nn_loaded

    logical function emt_nn_hubbard_loaded()
        emt_nn_hubbard_loaded = hubbard_loaded
    end function emt_nn_hubbard_loaded

    real(8) function emt_nn_a_lat()
        emt_nn_a_lat = a_lat_val
    end function emt_nn_a_lat

    subroutine emt_nn_box(blx, bly, bskew)
        real(8), intent(out) :: blx, bly, bskew
        blx = slab%box_lx
        bly = slab%box_ly
        bskew = slab%box_skew
    end subroutine emt_nn_box

    integer function emt_nn_natoms()
        emt_nn_natoms = slab%n_atoms
    end function emt_nn_natoms

    real(8) function emt_nn_rcut()
        emt_nn_rcut = rcut
    end function emt_nn_rcut

    real(8) function emt_nn_acut()
        emt_nn_acut = acut
    end function emt_nn_acut

    real(8) function emt_nn_s_ad()
        emt_nn_s_ad = emt_cached_s_ad
    end function emt_nn_s_ad

end module emt_nn_pes_venus
