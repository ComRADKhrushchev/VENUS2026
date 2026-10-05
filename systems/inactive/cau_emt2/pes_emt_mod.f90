!############################################################################
! pes_emt_mod.f90 — 2-Species Classical EMT Potential Module for VENUS
!
! Translated from md_tian2 pes_emt_mod.f90 (Goettingen group)
!   compute_emt_2species: lines 483-1027 of the original
!
! Strictly preserves the original formulas — only adapts the data
! interface from md_tian2 (universe_mod) to VENUS (Q/forces arrays).
!
! DEVIATION (VENUS2026 systems/cau_emt2, 2026-09-21): the PBC minimum
! image at the four wrap sites was the fractional-rectangle wrap
! (f - nint(f)), which is NOT a true minimum image for a skewed cell —
! pairs in the fractional corner region reconstruct to up to the
! half-diagonal (e.g. 15.33 A instead of the true 8.85 A for the 17.7 A
! hexagonal cell), giving discontinuous and partly wrong Hubbard
! pairwise sums. All four sites now use min_image_true (a 9-neighbor
! fractional-shift shortest reconstruction). Verified by
! scripts/check_cau_emt2.f90 criteria (d)/(h); the original formulas
! elsewhere are untouched.
!
! Energy formula (original line 976):
!   E = Ecoh_1 + Ecoh_2 - V_11 - V_22
!       - 0.5*(V_12 + V_21 - vref_1 - vref_2)
!
! 14 parameters: 7 EMT params × 2 species
!   eta2, n0, e0, lambda, v0, kappa, s0
!############################################################################

module pes_emt_mod
    implicit none
    private

    ! ================================================================
    ! Mathematical / physical constants
    ! ================================================================
    real(8), parameter :: pi_val   = 3.14159265358979323846d0
    real(8), parameter :: sqrt2    = 1.4142135623730951d0
    real(8), parameter :: sqrt3    = 1.7320508075688772d0
    real(8), parameter :: twelfth  = 1.0d0 / 12.0d0
    real(8), parameter :: tolerance = 1.0d-30
    real(8), parameter :: bohr2ang  = 0.52917721067d0
    real(8), parameter :: cutoff_factor = 1.5d0
    real(8), parameter :: smooth_factor = 1.0d4  ! controls FD cutoff sharpness

    ! FCC geometry factors (original: select_beta, select_nneighs)
    real(8), parameter :: beta_val = (16.0d0 * pi_val / 3.0d0)**(1.0d0/3.0d0) / sqrt2
    real(8), parameter :: nneighs(3) = [12.0d0, 6.0d0, 24.0d0]

    ! ================================================================
    ! 2-species EMT parameter type
    ! ================================================================
    type, public :: emt_2spec_params
        ! Raw parameters: index 1 = species A (adsorbate), 2 = species B (slab)
        real(8) :: eta2(2)   = 0.0d0   ! orbital exponent (1/Ang)
        real(8) :: n0(2)     = 0.0d0   ! reference density (1/Ang^3)
        real(8) :: e0(2)     = 0.0d0   ! cohesive energy (eV)
        real(8) :: lambda(2) = 0.0d0   ! embedding exponent (1/Ang)
        real(8) :: v0(2)     = 0.0d0   ! pair potential amplitude (eV)
        real(8) :: kappa(2)  = 0.0d0   ! pair potential exponent (1/Ang)
        real(8) :: s0(2)     = 0.0d0   ! neutral-sphere radius (Ang)
        real(8) :: energy_shift = 0.0d0 ! fixed energy shift (eV), added to total E
        ! Derived quantities (computed by read_emt_params)
        real(8) :: betas0(2)     = 0.0d0
        real(8) :: betaeta2(2)   = 0.0d0
        real(8) :: kappadbeta(2) = 0.0d0
        real(8) :: chi(2,2)      = 0.0d0
        real(8) :: rcut           = 0.0d0
        real(8) :: acut           = 0.0d0
        real(8) :: igamma1(2)     = 0.0d0
        real(8) :: igamma2(2)     = 0.0d0
    end type emt_2spec_params

    ! ================================================================
    ! Hubbard polynomial parameter type
    !
    ! p = a + [b*dens^(2/3) + c*dens^(4/3) + d*dens^2]*exp(A*(dens-B))
    !       + C * sum_j[theta(r_j)/(r_j - D)]
    ! ================================================================
    type :: hubbard_poly_params
        real(8) :: a       = 0.0d0   ! constant offset
        real(8) :: b       = 0.0d0   ! dens^(2/3) coefficient
        real(8) :: c       = 0.0d0   ! dens^(4/3) coefficient
        real(8) :: d       = 0.0d0   ! dens^2 coefficient
        real(8) :: A_exp   = 0.0d0   ! exponential amplitude
        real(8) :: B_shift = 0.0d0   ! exponential density shift
        real(8) :: C_amp   = 0.0d0   ! pairwise correction amplitude
        real(8) :: D_dist  = 0.0d0   ! pairwise distance offset (Ang)
    end type hubbard_poly_params

    ! --- Module-level cached Hubbard state (public) ---
    type(hubbard_poly_params) :: hubb_poly(4)    ! h, U, J, Vcp
    logical  :: hubbard_poly_loaded = .false.
    real(8)  :: ad_density = 0.0d0      ! cached density at impurity
    real(8), allocatable :: ad_density_grad(:,:) ! (3,natoms) density gradient
    real(8) :: emt_pair_sum(4) = 0.0d0           ! cached pairwise sums (h,U,J,Vcp)

    ! Public API
    public :: read_emt_params
    public :: compute_emt_2species
    public :: hubbard_poly_params
    public :: hubb_poly, hubbard_poly_loaded
    public :: ad_density, ad_density_grad
    public :: emt_pair_sum
    public :: eval_hubbard_poly, dpolyddens
    public :: compute_pairwise_sum, compute_pairwise_grad
    public :: fill_hubbard_param_grad
    public :: fd_cutoff

contains

    ! ================================================================
    ! min_image_true — true 2D minimum image for a (possibly skewed)
    ! cell, replacing the fractional-rectangle wrap (see the DEVIATION
    ! note in the file header): explore the 9 neighbor fractional
    ! shifts of the wrapped offsets and keep the shortest reconstructed
    ! Cartesian vector. dz is left unchanged (no z periodicity).
    ! ================================================================
    pure subroutine min_image_true(dx, dy, blx, bly, cskew, sskew, inv11, inv12, inv22)
        real(8), intent(inout) :: dx, dy
        real(8), intent(in) :: blx, bly, cskew, sskew, inv11, inv12, inv22
        real(8) :: fx0, fy0, fx, fy, tx, ty, d2, d2min
        integer :: sx, sy
        fx0 = inv11*dx + inv12*dy
        fy0 = inv22*dy
        fx0 = fx0 - nint(fx0)
        fy0 = fy0 - nint(fy0)
        d2min = huge(1.0d0)
        do sx = -1, 1
            do sy = -1, 1
                fx = fx0 + sx
                fy = fy0 + sy
                tx = blx*fx + bly*cskew*fy
                ty = bly*sskew*fy
                d2 = tx*tx + ty*ty
                if (d2 < d2min) then
                    d2min = d2
                    dx = tx
                    dy = ty
                end if
            end do
        end do
    end subroutine min_image_true

    ! ================================================================
    ! read_emt_params — read 14 EMT parameters + energy shift from file
    !
    ! File format: one keyword=value per line
    ! Keywords: ETA2_A, N0_A, E0_A, LAMBDA_A, V0_A, KAPPA_A, S0_A
    !           ETA2_B, N0_B, E0_B, LAMBDA_B, V0_B, KAPPA_B, S0_B
    !           ESHIFT
    !
    ! Computes all derived quantities (betas0, chi, rcut, igamma, etc.)
    ! ================================================================
    subroutine read_emt_params(filename, params, ios_out)
        character(len=*), intent(in)    :: filename
        type(emt_2spec_params), intent(out) :: params
        integer, intent(out) :: ios_out

        integer :: unit, ios, i, k
        real(8) :: rr
        real(8) :: rnn(3), x_val(3), r3temp(3)
        character(len=200) :: line
        character(len=50)  :: key
        real(8) :: val

        ! Open file
        open(newunit=unit, file=trim(filename), status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) "WARNING: Cannot open EMT params file: ", trim(filename)
            write(*,*) "  Using default Au parameters for both species"
            ios_out = ios
            call set_default_params(params)
            return
        end if

        ! Read keyword-value pairs
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle    ! skip empty lines
            if (line(1:1) == '#') cycle        ! skip comments
            if (line(1:1) == '!') cycle

            ! Split into key and value
            i = index(line, ' ')
            if (i == 0) cycle
            key = adjustl(line(1:i-1))
            read(line(i:), *, iostat=ios) val
            if (ios /= 0) cycle

            call lower_case(key)

            select case (trim(key))
            case ('eta2_a');    params%eta2(1)   = val
            case ('n0_a');      params%n0(1)     = val
            case ('e0_a');      params%e0(1)     = val
            case ('lambda_a');  params%lambda(1) = val
            case ('v0_a');      params%v0(1)     = val
            case ('kappa_a');   params%kappa(1)  = val
            case ('s0_a');      params%s0(1)     = val
            case ('eta2_b');    params%eta2(2)   = val
            case ('n0_b');      params%n0(2)     = val
            case ('e0_b');      params%e0(2)     = val
            case ('lambda_b');  params%lambda(2) = val
            case ('v0_b');      params%v0(2)     = val
            case ('kappa_b');   params%kappa(2)  = val
            case ('s0_b');      params%s0(2)     = val
            case ('eshift');    params%energy_shift = val
            ! Hubbard polynomial parameters (optional)
            case ('h_a');       hubb_poly(1)%a       = val; hubbard_poly_loaded = .true.
            case ('h_b');       hubb_poly(1)%b       = val
            case ('h_c');       hubb_poly(1)%c       = val
            case ('h_d');       hubb_poly(1)%d       = val
            case ('h_aexp');    hubb_poly(1)%A_exp   = val
            case ('h_bshift');  hubb_poly(1)%B_shift = val
            case ('h_camp');    hubb_poly(1)%C_amp   = val
            case ('h_ddist');   hubb_poly(1)%D_dist  = val
            case ('u_a');       hubb_poly(2)%a       = val; hubbard_poly_loaded = .true.
            case ('u_b');       hubb_poly(2)%b       = val
            case ('u_c');       hubb_poly(2)%c       = val
            case ('u_d');       hubb_poly(2)%d       = val
            case ('u_aexp');    hubb_poly(2)%A_exp   = val
            case ('u_bshift');  hubb_poly(2)%B_shift = val
            case ('u_camp');    hubb_poly(2)%C_amp   = val
            case ('u_ddist');   hubb_poly(2)%D_dist  = val
            case ('j_a');       hubb_poly(3)%a       = val; hubbard_poly_loaded = .true.
            case ('j_b');       hubb_poly(3)%b       = val
            case ('j_c');       hubb_poly(3)%c       = val
            case ('j_d');       hubb_poly(3)%d       = val
            case ('j_aexp');    hubb_poly(3)%A_exp   = val
            case ('j_bshift');  hubb_poly(3)%B_shift = val
            case ('j_camp');    hubb_poly(3)%C_amp   = val
            case ('j_ddist');   hubb_poly(3)%D_dist  = val
            case ('vcp_a');     hubb_poly(4)%a       = val; hubbard_poly_loaded = .true.
            case ('vcp_b');     hubb_poly(4)%b       = val
            case ('vcp_c');     hubb_poly(4)%c       = val
            case ('vcp_d');     hubb_poly(4)%d       = val
            case ('vcp_aexp');  hubb_poly(4)%A_exp   = val
            case ('vcp_bshift');hubb_poly(4)%B_shift = val
            case ('vcp_camp');  hubb_poly(4)%C_amp   = val
            case ('vcp_ddist'); hubb_poly(4)%D_dist  = val
            end select
        end do
        close(unit)
        ios_out = 0

        call compute_derived(params)
        call print_params(params)

    contains
        ! Lowercase conversion helper
        subroutine lower_case(str)
            character(len=*), intent(inout) :: str
            integer :: j, ic
            do j = 1, len_trim(str)
                ic = ichar(str(j:j))
                if (ic >= 65 .and. ic <= 90) str(j:j) = char(ic + 32)
            end do
        end subroutine lower_case
    end subroutine read_emt_params

    ! ================================================================
    ! set_default_params — Au default values when no file present
    ! ================================================================
    subroutine set_default_params(params)
        type(emt_2spec_params), intent(out) :: params
        params%eta2   = [3.197d0, 3.197d0]
        params%n0     = [0.042d0, 0.042d0]
        params%e0     = [-3.826d0, -3.826d0]
        params%lambda = [4.182d0, 4.182d0]
        params%v0     = [0.348d0, 0.348d0]
        params%kappa  = [3.249d0, 3.249d0]
        params%s0     = [1.642d0, 1.642d0]
        call compute_derived(params)
    end subroutine set_default_params

    ! ================================================================
    ! compute_derived — compute all derived quantities from raw params
    ! (original lines 537-592)
    ! ================================================================
    subroutine compute_derived(params)
        type(emt_2spec_params), intent(inout) :: params
        real(8) :: rr, rnn(3), x_val(3), r3temp(3)
        integer :: i, k

        ! beta * s0  (line 537-538)
        params%betas0(:) = beta_val * params%s0(:)

        ! beta * eta2  (line 541-542)
        params%betaeta2(:) = beta_val * params%eta2(:)

        ! kappa / beta  (line 547-548)
        params%kappadbeta(:) = params%kappa(:) / beta_val

        ! chi coupling parameters (lines 552-553)
        params%chi(2,1) = params%n0(1) / params%n0(2) * &
            exp(0.5d0 / bohr2ang * (params%s0(2) - params%s0(1)))
        params%chi(1,2) = 1.0d0 / params%chi(2,1)
        params%chi(1,1) = 1.0d0
        params%chi(2,2) = 1.0d0

        ! Cutoff: rcut = max(betas0) * sqrt3  (line 574)
        params%rcut = maxval(params%betas0) * sqrt3

        ! rr = 4*rcut / (sqrt3+2)  (line 576)
        rr = 4.0d0 * params%rcut / (sqrt3 + 2.0d0)
        params%acut = log(smooth_factor) / (rr - params%rcut)

        ! Gamma1, Gamma2 for each species (lines 578-592)
        do i = 1, 2
            ! NNN distances
            rnn(1) = params%betas0(i)
            rnn(2) = rnn(1) * sqrt2
            rnn(3) = rnn(1) * sqrt3

            ! x factors (using unified fd_cutoff)
            do k = 1, 3
                x_val(k) = nneighs(k) * twelfth &
                    * fd_cutoff(rnn(k), params%rcut, params%acut)
                r3temp(k) = rnn(k) - params%betas0(i)
            end do

            ! igamma1
            params%igamma1(i) = 1.0d0 / &
                sum(x_val * exp(-params%eta2(i) * r3temp))

            ! igamma2
            params%igamma2(i) = 1.0d0 / &
                sum(x_val * exp(-params%kappadbeta(i) * r3temp))
        end do
    end subroutine compute_derived

    ! ================================================================
    ! print_params — diagnostic output
    ! ================================================================
    subroutine print_params(params)
        type(emt_2spec_params), intent(in) :: params
        write(*, '(A)') "EMT 2-species parameters loaded:"
        write(*, '(A,7F10.4)') "  Species A: ", params%eta2(1), params%n0(1), &
            params%e0(1), params%lambda(1), params%v0(1), params%kappa(1), params%s0(1)
        write(*, '(A,7F10.4)') "  Species B: ", params%eta2(2), params%n0(2), &
            params%e0(2), params%lambda(2), params%v0(2), params%kappa(2), params%s0(2)
        write(*, '(A,2F15.8)') "  chi_AB, chi_BA: ", params%chi(1,2), params%chi(2,1)
        write(*, '(A,F10.4,A,F10.4)') "  rcut = ", params%rcut, "  acut = ", params%acut
        write(*, '(A,F12.6)') "  energy_shift = ", params%energy_shift
        if (hubbard_poly_loaded) then
            write(*, '(A)') "Hubbard polynomial parameters loaded:"
            write(*, '(A,8ES12.4)') "  h:   ", hubb_poly(1)%a, hubb_poly(1)%b, &
                hubb_poly(1)%c, hubb_poly(1)%d, hubb_poly(1)%A_exp, &
                hubb_poly(1)%B_shift, hubb_poly(1)%C_amp, hubb_poly(1)%D_dist
            write(*, '(A,8ES12.4)') "  U:   ", hubb_poly(2)%a, hubb_poly(2)%b, &
                hubb_poly(2)%c, hubb_poly(2)%d, hubb_poly(2)%A_exp, &
                hubb_poly(2)%B_shift, hubb_poly(2)%C_amp, hubb_poly(2)%D_dist
            write(*, '(A,8ES12.4)') "  J:   ", hubb_poly(3)%a, hubb_poly(3)%b, &
                hubb_poly(3)%c, hubb_poly(3)%d, hubb_poly(3)%A_exp, &
                hubb_poly(3)%B_shift, hubb_poly(3)%C_amp, hubb_poly(3)%D_dist
            write(*, '(A,8ES12.4)') "  Vcp: ", hubb_poly(4)%a, hubb_poly(4)%b, &
                hubb_poly(4)%c, hubb_poly(4)%d, hubb_poly(4)%A_exp, &
                hubb_poly(4)%B_shift, hubb_poly(4)%C_amp, hubb_poly(4)%D_dist
        end if
    end subroutine print_params

    ! ================================================================
    ! fd_cutoff — Fermi-Dirac smooth cutoff function
    !
    ! theta(r) = 1 / (1 + exp(acut * (r - rcut)))
    !
    ! The smoothing sharpness is controlled by acut:
    !   acut = log(smooth_factor) / (rr - rcut)
    ! where smooth_factor (default 10000) controls how sharply
    ! the cutoff transitions at r = rcut.
    !
    ! Optional output: dtheta_dr = -acut * theta * (1 - theta)
    ! ================================================================
    real(8) function fd_cutoff(r, rcut, acut, dtheta_dr)
        real(8), intent(in) :: r, rcut, acut
        real(8), intent(out), optional :: dtheta_dr
        real(8) :: exp_term

        exp_term = exp(acut * (r - rcut))
        fd_cutoff = 1.0d0 / (1.0d0 + exp_term)
        if (present(dtheta_dr)) then
            dtheta_dr = -acut * fd_cutoff * (1.0d0 - fd_cutoff)
        end if
    end function fd_cutoff

    ! ================================================================
    ! eval_hubbard_poly — evaluate Hubbard polynomial at density
    !
    ! p = a + [b*dens^(2/3) + c*dens^(4/3) + d*dens^2]*exp(A*(dens-B))
    !       + C * pair_sum
    !
    ! Reference: pot_emt.f90 lines 999-1015
    ! ================================================================
    pure real(8) function eval_hubbard_poly(density, poly, pair_sum)
        real(8), intent(in) :: density
        type(hubbard_poly_params), intent(in) :: poly
        real(8), intent(in) :: pair_sum
        real(8) :: dens23, dens43, dens2, poly_body, dens_clamp

        dens_clamp = max(density, 0.0d0)
        dens23 = dens_clamp ** (2.0d0 / 3.0d0)
        dens43 = dens_clamp ** (4.0d0 / 3.0d0)
        dens2  = dens_clamp * dens_clamp
        poly_body = poly%b * dens23 + poly%c * dens43 + poly%d * dens2
        
        eval_hubbard_poly = poly%a + poly_body * exp(poly%A_exp * (dens_clamp - poly%B_shift))+ poly%C_amp * pair_sum
    end function eval_hubbard_poly

    ! ================================================================
    ! dpolyddens — derivative of Hubbard polynomial w.r.t. density
    !
    ! dp/d(dens) = [b*(2/3)dens^(-1/3) + c*(4/3)dens^(1/3) + 2d*dens
    !               + A*(b*dens^(2/3) + c*dens^(4/3) + d*dens^2)]
    !            * exp(A*(dens-B))
    ! ================================================================
    pure real(8) function dpolyddens(density, poly)
        real(8), intent(in) :: density
        type(hubbard_poly_params), intent(in) :: poly
        real(8) :: dens23, dens43, dens2, body, dbody, dens_clamp

        dens_clamp = max(density, 1.0d-10)  ! avoid dens^(-1/3) divergence
        dens23 = dens_clamp ** (2.0d0 / 3.0d0)
        dens43 = dens_clamp ** (4.0d0 / 3.0d0)
        dens2  = dens_clamp * dens_clamp
        body   = poly%b * dens23 + poly%c * dens43 + poly%d * dens2
        dbody  = poly%b * (2.0d0/3.0d0) * dens23 / dens_clamp &
               + poly%c * (4.0d0/3.0d0) * dens43 / dens_clamp &
               + poly%d * 2.0d0 * dens_clamp
        dpolyddens = (dbody + poly%A_exp * body) &
                   * exp(poly%A_exp * (dens_clamp - poly%B_shift))
    end function dpolyddens

    ! ================================================================
    ! compute_pairwise_sum — pairwise correction term
    !
    ! pair_sum = sum_j theta(r_j) / max(r_j - D_dist, 1e-8)
    ! theta(r) = 1 / (1 + exp(acut*(r - rcut)))
    !
    ! Reference: pot_emt.f90 lines 1023-1050
    ! ================================================================
    real(8) function compute_pairwise_sum(Q_imp, Q_slab, n_slab, &
                                          D_dist, rcut, acut, &
                                          box_lx, box_ly, box_skew)
        real(8), intent(in) :: Q_imp(3)
        real(8), intent(in) :: Q_slab(3, n_slab)
        integer, intent(in) :: n_slab
        real(8), intent(in) :: D_dist, rcut, acut
        real(8), intent(in), optional :: box_lx, box_ly, box_skew
        real(8) :: dx, dy, dz, r, theta, blx, bly, bskew, fx, fy
        real(8) :: cskew, sskew, inv_11, inv_12, inv_22
        real(8), parameter :: acut_gentle = 0.3d0
        integer :: j
        logical :: use_pbc

        compute_pairwise_sum = 0.0d0
        use_pbc = present(box_lx) .and. present(box_ly) .and. present(box_skew)
        if (use_pbc) then
            blx = box_lx; bly = box_ly; bskew = box_skew
            cskew = cos(bskew); sskew = sin(bskew)
            if (abs(sskew) < 1.0d-12) then
                inv_11 = 1.0d0 / blx; inv_12 = 0.0d0; inv_22 = 1.0d0 / bly
            else
                inv_11 = 1.0d0 / blx
                inv_12 = -cskew / (blx * sskew)
                inv_22 = 1.0d0 / (bly * sskew)
            end if
        end if

        do j = 1, n_slab
            dx = Q_imp(1) - Q_slab(1, j)
            dy = Q_imp(2) - Q_slab(2, j)
            dz = Q_imp(3) - Q_slab(3, j)
            if (use_pbc) then
                call min_image_true(dx, dy, blx, bly, cskew, sskew, inv_11, inv_12, inv_22)
            end if
            r = sqrt(dx*dx + dy*dy + dz*dz)
            ! Gentle cutoff matching pot_emt.f90: no distance skip, acut=0.3
            theta = fd_cutoff(r, rcut, acut_gentle)
            compute_pairwise_sum = compute_pairwise_sum &
                + theta / max(r - D_dist, 1.0d-8)
        end do
    end function compute_pairwise_sum

    ! ================================================================
    ! compute_pairwise_grad — gradient of pairwise correction term
    !
    ! d/dQ_C [theta/(r-D)] = [dtheta/dr*(r-D) - theta]/(r-D)^2 * vec
    ! d/dQ_slab_j = -d/dQ_C (Newton's 3rd law)
    ! dtheta/dr = -acut * theta * (1 - theta)
    !
    ! Output: dpair_dQ(3, 1+n_slab) where column 1 = impurity, columns 2: = slab
    ! ================================================================
    subroutine compute_pairwise_grad(Q_imp, Q_slab, n_slab, &
                                     D_dist, rcut, acut, dpair_dQ, &
                                     box_lx, box_ly, box_skew)
        real(8), intent(in) :: Q_imp(3)
        real(8), intent(in) :: Q_slab(3, n_slab)
        integer, intent(in) :: n_slab
        real(8), intent(in) :: D_dist, rcut, acut
        real(8), intent(out) :: dpair_dQ(3, 1+n_slab)
        real(8), intent(in), optional :: box_lx, box_ly, box_skew
        real(8) :: dx, dy, dz, r, theta, dtheta_dr, r_minus_D, inv_denom
        real(8) :: dpair_dr, vec(3)
        real(8), parameter :: acut_gentle = 0.3d0
        real(8) :: blx, bly, bskew, fx, fy, cskew, sskew, inv_11, inv_12, inv_22
        integer :: j
        logical :: use_pbc

        dpair_dQ = 0.0d0
        use_pbc = present(box_lx) .and. present(box_ly) .and. present(box_skew)
        if (use_pbc) then
            blx = box_lx; bly = box_ly; bskew = box_skew
            cskew = cos(bskew); sskew = sin(bskew)
            if (abs(sskew) < 1.0d-12) then
                inv_11 = 1.0d0 / blx; inv_12 = 0.0d0; inv_22 = 1.0d0 / bly
            else
                inv_11 = 1.0d0 / blx
                inv_12 = -cskew / (blx * sskew)
                inv_22 = 1.0d0 / (bly * sskew)
            end if
        end if

        do j = 1, n_slab
            dx = Q_imp(1) - Q_slab(1, j)
            dy = Q_imp(2) - Q_slab(2, j)
            dz = Q_imp(3) - Q_slab(3, j)
            if (use_pbc) then
                call min_image_true(dx, dy, blx, bly, cskew, sskew, inv_11, inv_12, inv_22)
            end if
            r = sqrt(dx*dx + dy*dy + dz*dz)
            vec = [dx/r, dy/r, dz/r]
            ! Gentle cutoff matching pot_emt.f90
            theta = fd_cutoff(r, rcut, acut_gentle, dtheta_dr)
            r_minus_D = max(r - D_dist, 1.0d-8)
            inv_denom = 1.0d0 / (r_minus_D * r_minus_D)
            dpair_dr = (dtheta_dr * r_minus_D - theta) * inv_denom
            ! Impurity (column 1)
            dpair_dQ(:, 1) = dpair_dQ(:, 1) + dpair_dr * vec
            ! Slab atom j (column 1+j): Newton's 3rd law
            dpair_dQ(:, 1+j) = -dpair_dr * vec
        end do
    end subroutine compute_pairwise_grad

    ! ================================================================
    ! fill_hubbard_param_grad — gradient of one Hubbard parameter
    !
    ! d(param)/dQ_k = dpolyddens * ad_density_grad(dir,i)
    !               + C_amp * d(pair_sum)/dQ_k
    !
    ! This consolidates the repeated gradient-assembly logic from
    ! dH_orb (h, U, J) and dV_orb (Vcp).
    ! ================================================================
    subroutine fill_hubbard_param_grad(grad, ncoord, ipoly, &
                                        Q_imp, Q_slab, n_slab, rcut, acut, &
                                        natoma, box_lx, box_ly, box_skew)
        real(8), intent(inout) :: grad(ncoord)   ! output gradient (eV/Ang)
        integer, intent(in) :: ncoord            ! = 3*NATOMS
        integer, intent(in) :: ipoly             ! 1=h, 2=U, 3=J, 4=Vcp
        real(8), intent(in) :: Q_imp(3)          ! impurity position
        real(8), intent(in) :: Q_slab(3, n_slab) ! slab positions
        integer, intent(in) :: n_slab, natoma
        real(8), intent(in) :: rcut, acut
        real(8), intent(in) :: box_lx, box_ly, box_skew

        real(8) :: dpolyd_val
        real(8), allocatable :: dpair_dQ(:,:)
        integer :: i, i3, i_slab, dir

        ! Density contribution: dpolyddens * ad_density_grad
        dpolyd_val = dpolyddens(ad_density, hubb_poly(ipoly))
        do i = 1, ncoord/3
            i3 = 3*(i-1)
            do dir = 1, 3
                grad(i3+dir) = dpolyd_val * ad_density_grad(dir, i)
            end do
        end do

        ! Pairwise contribution (only if C_amp != 0)
        if (abs(hubb_poly(ipoly)%C_amp) > 1.0d-15) then
            allocate(dpair_dQ(3, 1+n_slab))
            call compute_pairwise_grad(Q_imp, Q_slab, n_slab, &
                hubb_poly(ipoly)%D_dist, rcut, acut, dpair_dQ, &
                box_lx=box_lx, box_ly=box_ly, box_skew=box_skew)
            ! Impurity (column 1 → VENUS indices 1:3)
            do dir = 1, 3
                grad(dir) = grad(dir) + hubb_poly(ipoly)%C_amp * dpair_dQ(dir, 1)
            end do
            ! Slab atoms (columns 2: → VENUS indices after NATOMA)
            do i_slab = 1, n_slab
                i3 = 3 * (natoma + i_slab - 1)
                do dir = 1, 3
                    grad(i3+dir) = grad(i3+dir) &
                        + hubb_poly(ipoly)%C_amp * dpair_dQ(dir, 1+i_slab)
                end do
            end do
            deallocate(dpair_dQ)
        end if
    end subroutine fill_hubbard_param_grad

    ! ================================================================
    ! compute_emt_2species — 2-species classical EMT energy + forces
    !
    ! Strict translation of compute_emt_2species from VenusNAMD
    ! (original lines 483-1027).
    !
    ! Arguments:
    !   natoms   : total number of atoms
    !   natoma   : number of species-A (adsorbate) atoms
    !   Q        : Cartesian coordinates (3, natoms) in Angstrom
    !   params   : EMT parameters (with derived quantities)
    !   Energy   : output total energy (eV)
    !   forces   : output forces (3, natoms) in eV/Angstrom
    ! ================================================================
    subroutine compute_emt_2species(natoms, natoma, Q, params, Energy, forces, &
                                     debug_terms, force_terms, box_lx, box_ly, box_skew)
        integer, intent(in) :: natoms, natoma
        real(8), intent(in) :: Q(3, natoms)
        type(emt_2spec_params), intent(in) :: params
        real(8), intent(out) :: Energy
        real(8), intent(out) :: forces(3, natoms)
        real(8), intent(out), optional :: debug_terms(8)
        real(8), intent(out), optional :: force_terms(8, 3, natoms)
        real(8), intent(in), optional :: box_lx, box_ly, box_skew

        integer :: i, j, k, idx_i, idx_j
        real(8) :: r, theta, rtemp, rtemp1, rtemp2, alpha, expterm, cutoff_dist, rmin, rmin_nopbc
        real(8) :: vec(3), dtheta(3), dx, dy, dz
        real(8) :: blx, bly, bskew
        real(8) :: fx, fy, inv_11, inv_12, inv_22, cskew, sskew
        integer :: imin, imin_nopbc
        logical :: use_pbc

        ! --- Per-atom scalar accumulators (original lines 502-509) ---
        real(8) :: sigma_22(natoms), sigma_11(natoms)
        real(8) :: sigma_21(natoms), sigma_12(natoms)
        real(8) :: V_22, V_11, V_21, V_12
        real(8) :: Ecoh_1, Ecoh_2, vref_1, vref_2
        real(8) :: s_1(natoms), s_2(natoms)
        real(8) :: exp_lambda_s1(natoms), exp_lambda_s2(natoms)

        ! --- Per-atom vector accumulators (original lines 510-517) ---
        real(8) :: dsigma_21_2(3, natoms), dsigma_12_1(3, natoms)
        real(8) :: dEcoh_2_2(3, natoms), dEcoh_2_1(3, natoms)
        real(8) :: dEcoh_1_2(3, natoms), dEcoh_1_1(3, natoms)
        real(8) :: dV_22_2(3, natoms), dV_11_1(3, natoms)
        real(8) :: dV_21_2(3, natoms), dV_21_1(3, natoms)
        real(8) :: dV_12_1(3, natoms), dV_12_2(3, natoms)
        real(8) :: dvref_2_2(3, natoms), dvref_2_1(3, natoms)
        real(8) :: dvref_1_2(3, natoms), dvref_1_1(3, natoms)

        ! --- 4-index arrays: allocatable to avoid stack overflow ---
        real(8), allocatable :: dsigma_22(:,:,:), dsigma_21_1(:,:,:)
        real(8), allocatable :: dsigma_11(:,:,:), dsigma_12_2(:,:,:)
        real(8), allocatable :: ds_2_2(:,:,:), ds_1_1(:,:,:)
        real(8), allocatable :: ds_2_1(:,:,:), ds_1_2(:,:,:)
        integer :: alloc_err

        ! Allocate large arrays on heap
        allocate(dsigma_22(3, natoms, natoms), dsigma_21_1(3, natoms, natoms), &
                 dsigma_11(3, natoms, natoms), dsigma_12_2(3, natoms, natoms), &
                 ds_2_2(3, natoms, natoms), ds_1_1(3, natoms, natoms), &
                 ds_2_1(3, natoms, natoms), ds_1_2(3, natoms, natoms), &
                 stat=alloc_err)
        if (alloc_err /= 0) then
            write(*,*) "ERROR: Failed to allocate 4D arrays in compute_emt_2species"
            Energy = 0.0d0; forces = 0.0d0; return
        end if

        ! ============================================================
        ! Initialize PBC — general 2D lattice vectors (hexagonal or orthogonal)
        !   L = [BOXLX  BOXLY*cos(SKEW)]
        !       [  0    BOXLY*sin(SKEW)]
        ! Precompute coefficients for Cartesian->fractional transform
        ! in the inner pair loop (no divisions, no trig calls per pair).
        ! SKEW=90° → cos=0, sin=1 → bitwise identical to simple Cartesian wrap.
        ! ============================================================
        use_pbc = present(box_lx) .and. present(box_ly) .and. present(box_skew)
        if (use_pbc) then
            blx = box_lx; bly = box_ly; bskew = box_skew
            cskew = cos(bskew)
            sskew = sin(bskew)
            if (abs(sskew) < 1.0d-12) then
                ! Degenerate cell (collinear vectors) — fall back to orthogonal
                inv_11 = 1.0d0 / blx
                inv_12 = 0.0d0
                inv_22 = 1.0d0 / bly
                cskew = 0.0d0
            else
                inv_11 = 1.0d0 / blx
                inv_12 = -cskew / (blx * sskew)
                inv_22 = 1.0d0 / (bly * sskew)
            end if
        end if

        ! ============================================================
        ! Initialize accumulators (original lines 603-644)
        ! ============================================================
        sigma_22 = 0.0d0; sigma_11 = 0.0d0
        sigma_12 = 0.0d0; sigma_21 = 0.0d0
        dsigma_22 = 0.0d0; dsigma_11 = 0.0d0
        dsigma_21_2 = 0.0d0; dsigma_12_1 = 0.0d0
        dsigma_12_2 = 0.0d0; dsigma_21_1 = 0.0d0

        s_1 = 0.0d0; s_2 = 0.0d0
        ds_1_1 = 0.0d0; ds_1_2 = 0.0d0
        ds_2_1 = 0.0d0; ds_2_2 = 0.0d0

        Ecoh_1 = 0.0d0; Ecoh_2 = 0.0d0
        dEcoh_1_1 = 0.0d0; dEcoh_1_2 = 0.0d0
        dEcoh_2_1 = 0.0d0; dEcoh_2_2 = 0.0d0

        V_22 = 0.0d0; V_11 = 0.0d0
        V_21 = 0.0d0; V_12 = 0.0d0
        dV_22_2 = 0.0d0; dV_11_1 = 0.0d0
        dV_21_2 = 0.0d0; dV_21_1 = 0.0d0
        dV_12_2 = 0.0d0; dV_12_1 = 0.0d0

        vref_2 = 0.0d0; vref_1 = 0.0d0
        dvref_2_2 = 0.0d0; dvref_2_1 = 0.0d0
        dvref_1_2 = 0.0d0; dvref_1_1 = 0.0d0

        cutoff_dist = cutoff_factor * params%rcut

        ! ============================================================
        ! Sigma and Pair-wise Contributions (original lines 646-763)
        ! Loop over all atom pairs
        ! ============================================================
        do i = 1, natoms - 1
            idx_i = type_idx(i, natoma)

            do j = i + 1, natoms
                idx_j = type_idx(j, natoma)

                ! Compute distance vector i -> j (original convention: r_i - r_j)
                dx = Q(1, i) - Q(1, j)
                dy = Q(2, i) - Q(2, j)
                dz = Q(3, i) - Q(3, j)
                if (use_pbc) then
                    ! 2D lattice-vector minimum image (no z-wrapping for surface)
                    call min_image_true(dx, dy, blx, bly, cskew, sskew, inv_11, inv_12, inv_22)
                end if
                r = sqrt(dx*dx + dy*dy + dz*dz)

                if (r < tolerance) cycle
                if (r > cutoff_dist) cycle

                vec(1) = dx / r
                vec(2) = dy / r
                vec(3) = dz / r

                ! Cut-off function (Fermi-Dirac, unified fd_cutoff)
                theta = fd_cutoff(r, params%rcut, params%acut)
                ! rtemp1 = acut * exp_term * theta = acut * (1 - theta)
                rtemp1 = params%acut * (1.0d0 - theta)

                if (idx_i == 2 .and. idx_j == 2) then
                    ! --- Both type 2 (slab-slab) (original lines 657-687) ---

                    ! sigma_22 contribution
                    rtemp = theta * exp(-params%eta2(2) * (r - params%betas0(2)))
                    sigma_22(i) = sigma_22(i) + rtemp
                    sigma_22(j) = sigma_22(j) + rtemp

                    dtheta = (params%eta2(2) + rtemp1) * rtemp * vec
                    dsigma_22(:, i, i) = dsigma_22(:, i, i) - dtheta
                    dsigma_22(:, j, j) = dsigma_22(:, j, j) + dtheta
                    dsigma_22(:, j, i) = dtheta
                    dsigma_22(:, i, j) = -dtheta

                    ! V_22 contribution
                    rtemp = theta * exp(-params%kappadbeta(2) * (r - params%betas0(2)))
                    V_22 = V_22 + rtemp

                    dtheta = (params%kappadbeta(2) + rtemp1) * rtemp * vec
                    dV_22_2(:, i) = dV_22_2(:, i) + dtheta
                    dV_22_2(:, j) = dV_22_2(:, j) - dtheta

                else if (idx_i == 1 .and. idx_j == 1) then
                    ! --- Both type 1 (ad-ad) (original lines 690-718) ---

                    ! sigma_11 contribution
                    rtemp = theta * exp(-params%eta2(1) * (r - params%betas0(1)))
                    sigma_11(i) = sigma_11(i) + rtemp
                    sigma_11(j) = sigma_11(j) + rtemp

                    dtheta = (params%eta2(1) + rtemp1) * rtemp * vec
                    dsigma_11(:, i, i) = dsigma_11(:, i, i) - dtheta
                    dsigma_11(:, j, j) = dsigma_11(:, j, j) + dtheta
                    dsigma_11(:, j, i) = dtheta
                    dsigma_11(:, i, j) = -dtheta

                    ! V_11 contribution
                    rtemp = theta * exp(-params%kappadbeta(1) * (r - params%betas0(1)))
                    V_11 = V_11 + rtemp

                    dtheta = (params%kappadbeta(1) + rtemp1) * rtemp * vec
                    dV_11_1(:, i) = dV_11_1(:, i) + dtheta
                    dV_11_1(:, j) = dV_11_1(:, j) - dtheta

                else if (idx_i /= idx_j) then
                    ! --- Cross-species (original lines 721-761) ---

                    ! sigma_21: type1 contribution TO type2 density
                    rtemp = theta * exp(-params%eta2(1) * (r - params%betas0(1)))
                    sigma_21(j) = sigma_21(j) + rtemp
                    dtheta = (params%eta2(1) + rtemp1) * rtemp * vec
                    dsigma_21_2(:, j) = dsigma_21_2(:, j) + dtheta
                    dsigma_21_1(:, i, j) = -dtheta

                    ! sigma_12: type2 contribution TO type1 density
                    rtemp = theta * exp(-params%eta2(2) * (r - params%betas0(2)))
                    sigma_12(i) = sigma_12(i) + rtemp
                    dtheta = (params%eta2(2) + rtemp1) * rtemp * vec
                    dsigma_12_1(:, i) = dsigma_12_1(:, i) - dtheta
                    dsigma_12_2(:, j, i) = dtheta

                    ! V_21
                    rtemp = theta * exp(-params%kappadbeta(1) * (r - params%betas0(1)))
                    V_21 = V_21 + rtemp
                    dtheta = (params%kappadbeta(1) + rtemp1) * rtemp * vec
                    dV_21_1(:, i) = dV_21_1(:, i) + dtheta
                    dV_21_2(:, j) = dV_21_2(:, j) - dtheta

                    ! V_12
                    rtemp = theta * exp(-params%kappadbeta(2) * (r - params%betas0(2)))
                    V_12 = V_12 + rtemp
                    dtheta = (params%kappadbeta(2) + rtemp1) * rtemp * vec
                    dV_12_1(:, i) = dV_12_1(:, i) + dtheta
                    dV_12_2(:, j) = dV_12_2(:, j) - dtheta

                end if
            end do
        end do

        ! ============================================================
        ! Divide by cut-off scaling factors (original lines 767-788)
        ! ============================================================
        sigma_22 = sigma_22 * params%igamma1(2)
        V_22     = V_22     * params%igamma2(2) * params%v0(2)
        sigma_11 = sigma_11 * params%igamma1(1)
        V_11     = V_11     * params%igamma2(1) * params%v0(1)
        sigma_21 = sigma_21 * params%igamma1(2)
        sigma_12 = sigma_12 * params%igamma1(1)
        V_21     = V_21     * params%igamma2(2) * params%v0(2) * params%chi(2,1)
        V_12     = V_12     * params%igamma2(1) * params%v0(1) * params%chi(1,2)

        dsigma_22   = dsigma_22   * params%igamma1(2)
        dsigma_11   = dsigma_11   * params%igamma1(1)
        dsigma_21_2 = dsigma_21_2 * params%igamma1(2)
        dsigma_21_1 = dsigma_21_1 * params%igamma1(2)
        dsigma_12_1 = dsigma_12_1 * params%igamma1(1)
        dsigma_12_2 = dsigma_12_2 * params%igamma1(1)
        dV_22_2     = dV_22_2     * params%igamma2(2) * params%v0(2)
        dV_11_1     = dV_11_1     * params%igamma2(1) * params%v0(1)
        dV_21_2     = dV_21_2     * params%igamma2(2) * params%v0(2) * params%chi(2,1)
        dV_21_1     = dV_21_1     * params%igamma2(2) * params%v0(2) * params%chi(2,1)
        dV_12_2     = dV_12_2     * params%igamma2(1) * params%v0(1) * params%chi(1,2)
        dV_12_1     = dV_12_1     * params%igamma2(1) * params%v0(1) * params%chi(1,2)

        ! ============================================================
        ! Neutral Sphere Radius (original lines 792-845)
        ! ============================================================

        ! Combined sigma: s_raw = sigma_self + chi_cross * sigma_cross
        do i = 1, natoms
            if (type_idx(i, natoma) == 1) then
                s_1(i) = sigma_11(i) + params%chi(1,2) * sigma_12(i)
            else
                s_2(i) = sigma_22(i) + params%chi(2,1) * sigma_21(i)
            end if
        end do

        ! Initialize ds = -dsigma (full array copy)
        ds_1_1 = -dsigma_11
        ds_2_2 = -dsigma_22

        ! Add cross-species ds contributions and divide by betaeta2 * s
        do i = 1, natoms
            if (type_idx(i, natoma) == 2) then
                ! Type 2 (slab)
                ds_2_2(:, i, i) = ds_2_2(:, i, i) - params%chi(2,1) * dsigma_21_2(:, i)
                do j = 1, natoms
                    if (type_idx(j, natoma) == 2) then
                        if (s_2(j) > tolerance) ds_2_2(:, i, j) = ds_2_2(:, i, j) / (params%betaeta2(2) * s_2(j))
                    else if (type_idx(j, natoma) == 1) then
                        if (s_1(j) > tolerance) ds_1_2(:, i, j) = -params%chi(1,2) * dsigma_12_2(:, i, j) / (params%betaeta2(1) * s_1(j))
                    end if
                end do
            else
                ! Type 1 (adsorbate)
                ds_1_1(:, i, i) = ds_1_1(:, i, i) - params%chi(1,2) * dsigma_12_1(:, i)
                do j = 1, natoms
                    if (type_idx(j, natoma) == 1) then
                        if (s_1(j) > tolerance) ds_1_1(:, i, j) = ds_1_1(:, i, j) / (params%betaeta2(1) * s_1(j))
                    else if (type_idx(j, natoma) == 2) then
                        if (s_2(j) > tolerance) ds_2_1(:, i, j) = -params%chi(2,1) * dsigma_21_1(:, i, j) / (params%betaeta2(2) * s_2(j))
                    end if
                end do
            end if
        end do

        ! Final neutral-sphere radius (original lines 844-845)
        do i = 1, natoms
            if (type_idx(i, natoma) == 1) then
                s_1(i) = -log(max(tolerance, s_1(i)) * twelfth) / params%betaeta2(1)
            else
                s_2(i) = -log(max(tolerance, s_2(i)) * twelfth) / params%betaeta2(2)
            end if
        end do

        ! ============================================================
        ! Cache Hubbard density & gradient at impurity (type-1 atom)
        !
        ! Uses sigma_12(1) — already igamma-scaled, NO chi factor
        ! (matches reference pot_emt.f90 compute_hubbard_density).
        ! d(dens)/dQ derived from existing dsigma_12_1/dsigma_12_2
        ! which are already igamma-scaled at lines 496-497.
        ! ============================================================
        if (hubbard_poly_loaded) then
            if (.not. allocated(ad_density_grad)) allocate(ad_density_grad(3, natoms))
            ad_density_grad = 0.0d0

            ! Use actual sigma (no clamp) for correct power-law behaviour
            rtemp = sigma_12(1)
            rtemp1 = params%betaeta2(1)
            rtemp2 = 0.5d0 / bohr2ang - rtemp1  ! eta_eff
            alpha = -rtemp2 / rtemp1              ! = 1 - 0.5/(bohr*betaeta2)
            ad_density = params%n0(1) * (rtemp * twelfth)**alpha

            ! d(dens)/dQ = n0 * alpha * (twelfth)^alpha * sigma^(alpha-1) * dsigma/dQ
            ! This avoids the 1/sigma singularity of the original formula.
            if (rtemp > 0.0d0) then
                rtemp = params%n0(1) * alpha * (twelfth**alpha) * rtemp**(alpha - 1.0d0)
            else
                rtemp = 0.0d0
            end if
            do i = 1, natoms
                if (type_idx(i, natoma) == 1) then
                    ! dsigma_12_1(:,i) already accumulated over all type-2 j
                    ad_density_grad(:, i) = rtemp * dsigma_12_1(:, i)
                end if
            end do
            ! For slab (type-2) atoms: sum dsigma_12_2(:, j, i) over type-1 i
            do j = 1, natoms
                if (type_idx(j, natoma) == 2) then
                    do i = 1, natoms
                        if (type_idx(i, natoma) == 1) then
                            ad_density_grad(:, j) = ad_density_grad(:, j) &
                                + rtemp * dsigma_12_2(:, j, i)
                        end if
                    end do
                end if
            end do
        end if

        ! ============================================================
        ! Cohesive Function (original lines 849-914)
        ! ============================================================
        rtemp1 = 0.5d0 / bohr2ang - params%betaeta2(1)
        rtemp2 = 0.5d0 / bohr2ang - params%betaeta2(2)

        do i = 1, natoms
            if (type_idx(i, natoma) == 1) then
                exp_lambda_s1(i) = exp(-params%lambda(1) * s_1(i))
                Ecoh_1 = Ecoh_1 + &
                    ((1.0d0 + params%lambda(1) * s_1(i)) * exp_lambda_s1(i) - 1.0d0)
            else
                exp_lambda_s2(i) = exp(-params%lambda(2) * s_2(i))
                Ecoh_2 = Ecoh_2 + &
                    ((1.0d0 + params%lambda(2) * s_2(i)) * exp_lambda_s2(i) - 1.0d0)
            end if
        end do

        Ecoh_1 = Ecoh_1 * params%e0(1)
        Ecoh_2 = Ecoh_2 * params%e0(2)

        ! dEcoh (original lines 871-914)
        do i = 1, natoms
            if (type_idx(i, natoma) == 2) then
                do j = 1, natoms
                    if (type_idx(j, natoma) == 2) then
                        dEcoh_2_2(:, i) = dEcoh_2_2(:, i) + &
                            s_2(j) * exp_lambda_s2(j) * ds_2_2(:, i, j)
                    else if (type_idx(j, natoma) == 1) then
                        dEcoh_1_2(:, i) = dEcoh_1_2(:, i) + &
                            s_1(j) * exp_lambda_s1(j) * ds_1_2(:, i, j)
                    end if
                end do
            else
                do j = 1, natoms
                    if (type_idx(j, natoma) == 2) then
                        dEcoh_2_1(:, i) = dEcoh_2_1(:, i) + &
                            s_2(j) * exp_lambda_s2(j) * ds_2_1(:, i, j)
                    else if (type_idx(j, natoma) == 1) then
                        dEcoh_1_1(:, i) = dEcoh_1_1(:, i) + &
                            s_1(j) * exp_lambda_s1(j) * ds_1_1(:, i, j)
                    end if
                end do
            end if
        end do

        dEcoh_2_2 = dEcoh_2_2 * params%e0(2) * params%lambda(2) * params%lambda(2)
        dEcoh_1_2 = dEcoh_1_2 * params%e0(1) * params%lambda(1) * params%lambda(1)
        dEcoh_2_1 = dEcoh_2_1 * params%e0(2) * params%lambda(2) * params%lambda(2)
        dEcoh_1_1 = dEcoh_1_1 * params%e0(1) * params%lambda(1) * params%lambda(1)

        ! ============================================================
        ! Reference Pair Potential (original lines 918-972)
        ! ============================================================
        do i = 1, natoms
            if (type_idx(i, natoma) == 2) then
                expterm = exp(-params%kappa(2) * s_2(i))
                vref_2 = vref_2 + expterm

                do j = 1, natoms
                    if (type_idx(j, natoma) == 2) then
                        dvref_2_2(:, j) = dvref_2_2(:, j) + &
                            expterm * ds_2_2(:, j, i)
                    else if (type_idx(j, natoma) == 1) then
                        dvref_2_1(:, j) = dvref_2_1(:, j) + &
                            expterm * ds_2_1(:, j, i)
                    end if
                end do
            else
                expterm = exp(-params%kappa(1) * s_1(i))
                vref_1 = vref_1 + expterm

                do j = 1, natoms
                    if (type_idx(j, natoma) == 1) then
                        dvref_1_1(:, j) = dvref_1_1(:, j) + &
                            expterm * ds_1_1(:, j, i)
                    else if (type_idx(j, natoma) == 2) then
                        dvref_1_2(:, j) = dvref_1_2(:, j) + &
                            expterm * ds_1_2(:, j, i)
                    end if
                end do
            end if
        end do

        rtemp     = 12.0d0 * params%v0(2)
        vref_2    = vref_2 * rtemp
        dvref_2_2 = dvref_2_2 * rtemp * params%kappa(2)
        dvref_2_1 = dvref_2_1 * rtemp * params%kappa(2)

        rtemp     = 12.0d0 * params%v0(1)
        vref_1    = vref_1 * rtemp
        dvref_1_2 = dvref_1_2 * rtemp * params%kappa(1)
        dvref_1_1 = dvref_1_1 * rtemp * params%kappa(1)

        ! ============================================================
        ! TOTAL ENERGY (original line 976) + fixed energy shift
        ! ============================================================
        Energy = Ecoh_1 + Ecoh_2 - V_11 - V_22 &
               - 0.5d0 * (V_12 + V_21 - vref_1 - vref_2) &
               + params%energy_shift

        ! Optional debug output of energy components
        if (present(debug_terms)) then
            debug_terms(1) = Ecoh_1
            debug_terms(2) = Ecoh_2
            debug_terms(3) = V_11
            debug_terms(4) = V_22
            debug_terms(5) = V_12
            debug_terms(6) = V_21
            debug_terms(7) = vref_1
            debug_terms(8) = vref_2
        end if

        ! ============================================================
        ! TOTAL FORCES (original lines 1001-1011)
        ! ============================================================
        forces = 0.0d0

        do i = 1, natoms
            if (type_idx(i, natoma) == 2) then
                ! Type 2 (slab)
                forces(:, i) = dEcoh_2_2(:, i) + dEcoh_1_2(:, i) &
                             - dV_22_2(:, i) &
                             - 0.5d0 * (dV_21_2(:, i) + dV_12_2(:, i) &
                                      - dvref_2_2(:, i) - dvref_1_2(:, i))
            else
                ! Type 1 (adsorbate)
                forces(:, i) = dEcoh_2_1(:, i) + dEcoh_1_1(:, i) &
                             - dV_11_1(:, i) &
                             - 0.5d0 * (dV_21_1(:, i) + dV_12_1(:, i) &
                                      - dvref_2_1(:, i) - dvref_1_1(:, i))
            end if
        end do

        ! Optional: pack force components + debug intermediates
        if (present(force_terms)) then
            ! Debug: dump sigma/s/ds for slab atom 2
            if (natoms >= 2 .and. type_idx(2, natoma) == 2) then
                write(6, '(A)') '--- DEBUG sigma/s/ds for slab atom 2 ---'
                write(6, '(A,2ES14.5)') '  sigma_22(2), sigma_21(2)=', sigma_22(2), sigma_21(2)
                write(6, '(A,ES14.5)') '  s_2(2) [sigma_raw]=', s_2(2)
                write(6, '(A,ES14.5)') '  s_2(2) [after log]=', -log(max(tolerance, sigma_22(2) + params%chi(2,1)*sigma_21(2))*twelfth)/params%betaeta2(2)
                write(6, '(A,3ES14.5)') '  dsigma_22(:,2,2)=', dsigma_22(:,2,2)
                write(6, '(A,3ES14.5)') '  dsigma_21_2(:,2)=', dsigma_21_2(:,2)
                write(6, '(A,3ES14.5)') '  ds_2_2(:,2,2)=', ds_2_2(:,2,2)
                write(6, '(A,3ES14.5)') '  ds_2_2(:,2,3)=', ds_2_2(:,2,3)
                write(6, '(A,3ES14.5)') '  dsigma_22(:,2,3)=', dsigma_22(:,2,3)
                write(6, '(A,ES14.5)') '  Ecoh2 contrib: s₂*exp(-λs₂)=', s_2(2)*exp(-params%lambda(2)*s_2(2))
                write(6, '(A,ES14.5)') '  dEcoh factor: e0*λ²=', params%e0(2)*params%lambda(2)*params%lambda(2)
                write(6, '(A,ES14.5)') '  vref factor: 12*v0*κ=', 12.0d0*params%v0(2)*params%kappa(2)
                write(6, '(A)') '-----------------------------------------'
            end if
            force_terms = 0.0d0
            do i = 1, natoms
                if (type_idx(i, natoma) == 2) then
                    ! Slab atom
                    force_terms(1, :, i) = dEcoh_2_2(:, i)
                    force_terms(2, :, i) = dEcoh_1_2(:, i)
                    force_terms(3, :, i) = -dV_22_2(:, i)
                    force_terms(4, :, i) = -0.5d0 * (dV_21_2(:, i) + dV_12_2(:, i))
                    force_terms(5, :, i) = 0.5d0 * (dvref_2_2(:, i) + dvref_1_2(:, i))
                    ! 6-8 are zero for slab
                else
                    ! Adatom
                    force_terms(6, :, i) = dEcoh_2_1(:, i)
                    force_terms(7, :, i) = dEcoh_1_1(:, i)
                    force_terms(8, :, i) = -dV_11_1(:, i)  - 0.5d0 * (dV_21_1(:, i) + dV_12_1(:, i) - dvref_2_1(:, i) - dvref_1_1(:, i))
                    ! 1-5 are zero for adatom
                end if
            end do
        end if

        ! Deallocate large arrays
        deallocate(dsigma_22, dsigma_21_1, dsigma_11, dsigma_12_2, &
                   ds_2_2, ds_1_1, ds_2_1, ds_1_2)

    end subroutine compute_emt_2species

    ! ================================================================
    ! type_idx — helper: map atom index to species type (1 or 2)
    ! ================================================================
    pure integer function type_idx(i, natoma)
        integer, intent(in) :: i, natoma
        if (i <= natoma) then
            type_idx = 1   ! species A (adsorbate)
        else
            type_idx = 2   ! species B (slab)
        end if
    end function type_idx

    ! ================================================================
    ! pbc_min_image_2d — 2D minimum image with general lattice vectors
    !
    ! Cell defined by BOXLX, BOXLY, SKEW (radians):
    !   Lx = (box_lx, 0)
    !   Ly = (box_ly * cos(skew), box_ly * sin(skew))
    !
    ! Transforms to fractional coordinates, wraps to [-0.5, 0.5],
    ! transforms back. With SKEW=90° (pi/2) this is bitwise-identical
    ! to the old Cartesian nint wrapping.
    !
    ! dz is left unchanged (no z-periodicity for surface slab).
    ! ================================================================
    subroutine pbc_min_image_2d(dx, dy, box_lx, box_ly, skew)
        real(8), intent(inout) :: dx, dy
        real(8), intent(in) :: box_lx, box_ly, skew
        real(8) :: sin_skew, cos_skew, fx, fy

        sin_skew = sin(skew)
        cos_skew = cos(skew)

        ! Guard against degenerate cells (nearly collinear lattice vectors)
        if (abs(sin_skew) < 1.0d-12) then
            dx = dx - box_lx * nint(dx / box_lx)
            dy = dy - box_ly * nint(dy / box_ly)
            return
        end if

        ! 1. Cartesian -> fractional coordinates:  f = inv(L) * dr
        fx = dx / box_lx - dy * cos_skew / (box_lx * sin_skew)
        fy = dy / (box_ly * sin_skew)

        ! 2-3. True minimum image (9-neighbor fractional shifts; see the
        !      DEVIATION note in the file header) - fractional rectangle
        !      wrap first, then the shortest of the 9 reconstructions
        call min_image_true(dx, dy, box_lx, box_ly, cos_skew, sin_skew, &
                            1.0d0/box_lx, -cos_skew/(box_lx*sin_skew), 1.0d0/(box_ly*sin_skew))

    end subroutine pbc_min_image_2d

end module pes_emt_mod
