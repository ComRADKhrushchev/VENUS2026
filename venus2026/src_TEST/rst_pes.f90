! ================================================================
! rst_pes.f90 — Zero-dependency RST pairwise NN PES module
!
! Provides: FCC(111) slab generation, pairwise 1D-NN interatomic
!           potential with distance input, energy as sum over
!           adsorbate–slab atom pairs within cutoff, site-organized
!           data file I/O, and standalone NN weight file I/O.
!
! This module has NO dependencies on other project modules.
! It can be compiled standalone for Venus MD interface use.
! ================================================================
module rst_pes
    implicit none
    private

    ! --- Constants ---
    real(8), parameter :: pi = 3.14159265358979323846d0
    real(8), parameter :: sqrt3 = 1.7320508075688772d0
    real(8), parameter :: sqrt2_over_3 = 0.8164965809277260d0   ! sqrt(2/3)
    real(8), parameter :: default_a = 2.95d0                     ! Au NN distance
    real(8), parameter :: default_r_cut = 6.0d0
    real(8), parameter :: default_d_cut = 0.5d0
    real(8), parameter :: default_r_scale = 3.0d0

    ! --- Slab data type ---
    type :: rst_slab_data
        integer :: n_atoms = 0
        integer :: n_surf = 0               ! number of surface-layer atoms
        real(8), allocatable :: coords(:,:) ! (3, n_atoms)
        integer, allocatable :: layer_id(:) ! 1 = surface, 2 = inner
        real(8) :: box_lx = 36.0d0
        real(8) :: box_ly = 36.0d0
    end type rst_slab_data

    ! --- Single-state pairwise NN type ---
    type :: rst_state
        integer :: nh_surf = 0              ! hidden neurons, surface layer
        integer :: nh_inner = 0             ! hidden neurons, inner layer
        real(8) :: r_cut = default_r_cut    ! cutoff radius (Angstrom)
        real(8) :: d_cut = default_d_cut    ! smoothing width (Angstrom)
        real(8) :: r_scale = default_r_scale ! distance scaling factor
        real(8) :: A_rep = 0.0d0            ! repulsive core amplitude (eV)
        real(8) :: B_rep = 0.0d0            ! repulsive core exponent (1/Angstrom)
        real(8), allocatable :: w(:)        ! flat weights
    end type rst_state

    ! --- Per-site DFT data type ---
    type :: rst_site_data
        character(len=8) :: site_name = ""
        real(8) :: x = 0.0d0, y = 0.0d0
        integer :: npts = 0
        real(8), allocatable :: z(:)        ! z values (Angstrom)
        real(8), allocatable :: wgt(:)      ! weights
        real(8), allocatable :: E(:)        ! DFT energies (eV)
    end type rst_site_data

    ! --- State name labels ---
    character(len=3), parameter :: rst_state_names(7) = &
        ["C+ ", "1D ", "3P ", "AN ", "Vcp", "E_h", "GS "]

    ! --- Site names for 4 high-symmetry sites ---
    character(len=8), parameter :: site_names(4) = &
        ["Top     ", "Bridge  ", "HCP     ", "FCC     "]

    public :: rst_slab_data, rst_state, rst_site_data
    public :: rst_state_names, site_names
    public :: generate_fcc_slab, load_slab_xyz, write_slab_xyz
    public :: get_site_xy, get_site_index
    public :: pbc_min_image, smooth_cutoff
    public :: eval_1d_nn_dist, eval_repulsive_core
    public :: calc_rst_energy, calc_rst_energy_grad
    public :: read_rst_data, deallocate_rst_sites
    public :: load_rst_weights_standalone, save_rst_weights_standalone

contains

    ! ================================================================
    ! Generate FCC(111) slab
    !
    ! 2D hexagonal lattice basis: a1 = (a, 0), a2 = (a/2, a*sqrt3/2)
    ! ABC stacking: layer offset cycles through (0,0), (a/2, a/(2*s3)), (0, a/s3)
    ! Interlayer spacing: d = a * sqrt(2/3) ≈ 2.408 Angstrom
    ! ================================================================
    subroutine generate_fcc_slab(a, n_layers, n_side, slab)
        real(8), intent(in) :: a
        integer, intent(in) :: n_layers, n_side
        type(rst_slab_data), intent(out) :: slab

        real(8) :: d_layer, ox, oy, x, y, z, z_max
        integer :: ilay, ltype, i, j, idx, n_per_layer

        d_layer = a * sqrt2_over_3
        n_per_layer = n_side * n_side
        slab%n_atoms = n_layers * n_per_layer

        if (allocated(slab%coords)) deallocate(slab%coords)
        if (allocated(slab%layer_id)) deallocate(slab%layer_id)
        allocate(slab%coords(3, slab%n_atoms))
        allocate(slab%layer_id(slab%n_atoms))

        ! Set box dimensions to contain the rhombic atom arrangement
        slab%box_lx = n_side * a + a     ! pad by one lattice spacing
        slab%box_ly = n_side * a * sqrt3 / 2.0d0 + a

        idx = 0
        do ilay = 0, n_layers - 1
            ltype = mod(ilay, 3)
            select case(ltype)
            case(0)  ! A-type: no offset
                ox = 0.0d0
                oy = 0.0d0
            case(1)  ! B-type: FCC hollow offset
                ox = a / 2.0d0
                oy = a / (2.0d0 * sqrt3)
            case(2)  ! C-type: HCP hollow offset
                ox = 0.0d0
                oy = a / sqrt3
            end select
            z = -dble(ilay) * d_layer

            do i = 0, n_side - 1
                do j = 0, n_side - 1
                    idx = idx + 1
                    x = dble(i) * a + dble(j) * a / 2.0d0 + ox
                    y = dble(j) * a * sqrt3 / 2.0d0 + oy

                    ! Shift so a surface-layer atom sits at (0,0) = Top site
                    x = x - (n_side / 2) * a * 1.5d0
                    y = y - (n_side / 2) * a * sqrt3 / 2.0d0

                    slab%coords(1, idx) = x
                    slab%coords(2, idx) = y
                    slab%coords(3, idx) = z
                end do
            end do
        end do

        ! Auto-classify layers: highest-z layer = surface, rest = inner
        z_max = maxval(slab%coords(3, :))
        slab%n_surf = 0
        do i = 1, slab%n_atoms
            if (slab%coords(3, i) > z_max - 0.5d0 * d_layer) then
                slab%layer_id(i) = 1   ! surface
                slab%n_surf = slab%n_surf + 1
            else
                slab%layer_id(i) = 2   ! inner
            end if
        end do

        write(*, "(A,I0,A,I0,A,F8.3,A,F8.3,A,F6.3)") &
            "  FCC slab generated: ", slab%n_atoms, " atoms (", &
            slab%n_surf, " surface), box ", slab%box_lx, " x ", &
            slab%box_ly, " Ang^2, d_layer = ", d_layer
    end subroutine generate_fcc_slab

    ! ================================================================
    ! Load slab from XYZ file and auto-classify layers
    ! Expects standard XYZ: line1=n_atoms, line2=comment, then Au x y z
    ! ================================================================
    subroutine load_slab_xyz(filename, slab)
        character(len=*), intent(in) :: filename
        type(rst_slab_data), intent(out) :: slab

        integer :: unit, ios, n, i
        real(8) :: x, y, z, z_max, d_layer_est
        character(len=2) :: elem
        character(len=200) :: line

        open(newunit=unit, file=trim(filename), status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Cannot open slab file: ", trim(filename)
            stop
        end if

        read(unit, *, iostat=ios) n
        if (ios /= 0) then
            write(*,*) "ERROR: Failed to read atom count from ", trim(filename)
            close(unit); stop
        end if
        read(unit, '(A)', iostat=ios) line   ! skip comment line

        slab%n_atoms = n
        if (allocated(slab%coords)) deallocate(slab%coords)
        if (allocated(slab%layer_id)) deallocate(slab%layer_id)
        allocate(slab%coords(3, n))
        allocate(slab%layer_id(n))

        do i = 1, n
            read(unit, *, iostat=ios) elem, x, y, z
            if (ios /= 0) then
                write(*,*) "ERROR: Failed to read atom ", i, " from ", trim(filename)
                close(unit); stop
            end if
            slab%coords(1, i) = x
            slab%coords(2, i) = y
            slab%coords(3, i) = z
        end do
        close(unit)

        ! Estimate box dimensions from coordinate extent + 2*A margin
        slab%box_lx = maxval(slab%coords(1, :)) - minval(slab%coords(1, :)) + 6.0d0
        slab%box_ly = maxval(slab%coords(2, :)) - minval(slab%coords(2, :)) + 6.0d0

        ! Auto-classify layers: find z-gaps, topmost layer = surface
        z_max = maxval(slab%coords(3, :))
        d_layer_est = default_a * sqrt2_over_3   ! ~2.408 A for Au
        slab%n_surf = 0
        do i = 1, n
            if (slab%coords(3, i) > z_max - 0.5d0 * d_layer_est) then
                slab%layer_id(i) = 1
                slab%n_surf = slab%n_surf + 1
            else
                slab%layer_id(i) = 2
            end if
        end do

        write(*, "(A,A,A,I0,A,I0,A)") "  Slab loaded from ", trim(filename), &
            ": ", slab%n_atoms, " atoms, ", slab%n_surf, " surface"
    end subroutine load_slab_xyz

    ! ================================================================
    ! Write slab to standard XYZ file
    ! ================================================================
    subroutine write_slab_xyz(filename, slab)
        character(len=*), intent(in) :: filename
        type(rst_slab_data), intent(in) :: slab
        integer :: unit, i
        open(newunit=unit, file=filename, status='replace')
        write(unit, '(I0)') slab%n_atoms
        write(unit, '(A)') "# Au(111) slab generated by AnalyticModel"
        do i = 1, slab%n_atoms
            write(unit, '(A,3F14.6)') "Au", slab%coords(1,i), slab%coords(2,i), slab%coords(3,i)
        end do
        close(unit)
        write(*,*) "Slab XYZ written: ", trim(filename), " (", slab%n_atoms, " atoms)"
    end subroutine write_slab_xyz

    ! ================================================================
    ! Get (x,y) coordinates for 4 high-symmetry sites on Au(111)
    !
    ! site 1 = Top     → (0, 0)
    ! site 2 = Bridge  → (a/2, 0)
    ! site 3 = FCC     → (a/2, a/(2*sqrt3))   hollow with 2nd-layer atom below
    ! site 4 = HCP     → (0, a/sqrt3)          hollow with 3rd-layer atom below
    ! ================================================================
    subroutine get_site_xy(isite, a, x, y, name)
        integer, intent(in) :: isite
        real(8), intent(in) :: a
        real(8), intent(out) :: x, y
        character(len=*), intent(out) :: name

        select case(isite)
        case(1)  ! Top
            x = 0.0d0
            y = 0.0d0
            name = "Top"
        case(2)  ! Bridge
            x = a / 2.0d0
            y = 0.0d0
            name = "Bridge"
        case(3)  ! HCP hollow (B-layer atom directly below)
            x = a / 2.0d0
            y = a / (2.0d0 * sqrt3)
            name = "HCP"
        case(4)  ! FCC hollow (no B-layer atom below, C-layer below)
            x = 0.0d0
            y = a / sqrt3
            name = "FCC"
        case default
            x = 0.0d0; y = 0.0d0; name = "Unknown"
        end select
    end subroutine get_site_xy

    ! ================================================================
    ! Look up site index from name string (case-insensitive prefix)
    ! Returns 0 if not found
    ! ================================================================
    integer function get_site_index(name)
        character(len=*), intent(in) :: name
        integer :: i
        get_site_index = 0
        do i = 1, 4
            if (index(name, trim(site_names(i))) > 0 .or. &
                index(name, trim(adjustl(site_names(i)))) > 0) then
                get_site_index = i
                return
            end if
        end do
        ! Also try case-insensitive
        do i = 1, 4
            if (len_trim(name) >= len_trim(site_names(i))) then
                ! Simple prefix match
            end if
        end do
    end function get_site_index

    ! ================================================================
    ! Minimum image convention for (dx, dy) in rectangular box
    ! ================================================================
    subroutine pbc_min_image(dx, dy, box_lx, box_ly)
        real(8), intent(inout) :: dx, dy
        real(8), intent(in) :: box_lx, box_ly

        dx = dx - box_lx * nint(dx / box_lx)
        dy = dy - box_ly * nint(dy / box_ly)
    end subroutine pbc_min_image

    ! ================================================================
    ! Smooth cutoff function (cosine taper)
    !
    ! f(r) = 1                            for r <= r_cut - d_cut
    !      = 0.5*(1 + cos(pi*t))          for r_cut - d_cut < r < r_cut
    !      = 0                            for r >= r_cut
    ! where t = (r - (r_cut - d_cut)) / d_cut
    ! ================================================================
    subroutine smooth_cutoff(r, r_cut, d_cut, f, dfdr)
        real(8), intent(in) :: r, r_cut, d_cut
        real(8), intent(out) :: f, dfdr

        real(8) :: t, r_inner

        r_inner = r_cut - d_cut

        if (r <= r_inner) then
            f = 1.0d0
            dfdr = 0.0d0
        else if (r >= r_cut) then
            f = 0.0d0
            dfdr = 0.0d0
        else
            t = (r - r_inner) / d_cut
            f = 0.5d0 * (1.0d0 + cos(pi * t))
            dfdr = -0.5d0 * pi / d_cut * sin(pi * t)
        end if
    end subroutine smooth_cutoff

    ! ================================================================
    ! 1D Neural Network: V(r) with analytical dV/dr
    !
    ! Input: distance r (scalar), weights(:) size = 3*nh + 1
    ! Layout: w1(1:nh), b1(nh+1:2*nh), w2(2*nh+1:3*nh), b2(3*nh+1)
    ! ================================================================
    subroutine eval_1d_nn_dist(r, weights, nh, r_scale, V, dVdr)
        real(8), intent(in) :: r
        real(8), intent(in) :: weights(:)
        integer, intent(in) :: nh
        real(8), intent(in) :: r_scale
        real(8), intent(out) :: V, dVdr

        integer :: i
        real(8) :: w1, b1, w2, h, rs, dh_dr

        rs = r / r_scale

        V = weights(3 * nh + 1)     ! b2 (output bias)
        dVdr = 0.0d0

        do i = 1, nh
            w1 = weights(i)
            b1 = weights(nh + i)
            w2 = weights(2 * nh + i)

            h = tanh(w1 * rs + b1)
            V = V + w2 * h
            dh_dr = (1.0d0 - h * h) * w1 / r_scale
            dVdr = dVdr + w2 * dh_dr
        end do
    end subroutine eval_1d_nn_dist

    ! ================================================================
    ! Repulsive core for short-range Pauli repulsion
    !
    ! V_rep(r) = A_rep * exp(-B_rep * r)
    ! dV_rep/dr = -B_rep * V_rep
    !
    ! Returns 0 if A_rep <= 0 (no repulsive core).
    ! ================================================================
    subroutine eval_repulsive_core(r, A_rep, B_rep, V_rep, dVdr_rep)
        real(8), intent(in) :: r, A_rep, B_rep
        real(8), intent(out) :: V_rep, dVdr_rep

        real(8) :: exp_br

        if (A_rep > 0.0d0) then
            exp_br = exp(-B_rep * r)
            V_rep = A_rep * exp_br
            dVdr_rep = -B_rep * V_rep
        else
            V_rep = 0.0d0
            dVdr_rep = 0.0d0
        end if
    end subroutine eval_repulsive_core

    ! ================================================================
    ! Total RST energy: sum over adsorbate-slab atom pairs
    !
    ! E(xyz) = E0 + sum_i [ f_cut(r_i) * V_layer(r_i) ]
    ! ================================================================
    subroutine calc_rst_energy(xyz, slab, state, Energy)
        real(8), intent(in) :: xyz(3)
        type(rst_slab_data), intent(in) :: slab
        type(rst_state), intent(in) :: state
        real(8), intent(out) :: Energy

        real(8) :: dEdx, dEdy, dEdz

        call calc_rst_energy_grad(xyz, slab, state, Energy, dEdx, dEdy, dEdz)
    end subroutine calc_rst_energy

    ! ================================================================
    ! Total RST energy + Cartesian gradient
    !
    ! dE/dx_ad = sum_i [ (dV/dr * f_cut + V * df_cut/dr) * (x_ad - X_i) / r_i ]
    ! ================================================================
    subroutine calc_rst_energy_grad(xyz, slab, state, Energy, dEdx, dEdy, dEdz)
        real(8), intent(in) :: xyz(3)
        type(rst_slab_data), intent(in) :: slab
        type(rst_state), intent(in) :: state
        real(8), intent(out) :: Energy, dEdx, dEdy, dEdz

        integer :: i, n_surf_params, n_inner_params
        real(8) :: dx, dy, dz, r, f, dfdr, V, dVdr, dE_dr, E0
        real(8) :: V_rep, dVdr_rep
        real(8) :: w_surf_vals(3 * state%nh_surf + 1)
        real(8) :: w_inner_vals(3 * state%nh_inner + 1)

        n_surf_params = 3 * state%nh_surf + 1
        n_inner_params = 3 * state%nh_inner + 1

        ! Extract layer-specific weight blocks
        if (state%nh_surf > 0) then
            w_surf_vals = state%w(1 : n_surf_params)
        end if
        if (state%nh_inner > 0) then
            w_inner_vals = state%w(n_surf_params + 1 : n_surf_params + n_inner_params)
        end if
        E0 = state%w(n_surf_params + n_inner_params + 1)

        Energy = E0
        dEdx = 0.0d0
        dEdy = 0.0d0
        dEdz = 0.0d0

        do i = 1, slab%n_atoms
            dx = xyz(1) - slab%coords(1, i)
            dy = xyz(2) - slab%coords(2, i)
            dz = xyz(3) - slab%coords(3, i)

            call pbc_min_image(dx, dy, slab%box_lx, slab%box_ly)

            r = sqrt(dx*dx + dy*dy + dz*dz)
            if (r > state%r_cut) cycle

            call smooth_cutoff(r, state%r_cut, state%d_cut, f, dfdr)

            ! Evaluate NN for this pair based on slab atom layer
            if (slab%layer_id(i) == 1 .and. state%nh_surf > 0) then
                call eval_1d_nn_dist(r, w_surf_vals, state%nh_surf, &
                                     state%r_scale, V, dVdr)
            else if (slab%layer_id(i) == 2 .and. state%nh_inner > 0) then
                call eval_1d_nn_dist(r, w_inner_vals, state%nh_inner, &
                                     state%r_scale, V, dVdr)
            else
                V = 0.0d0
                dVdr = 0.0d0
            end if

            ! Add repulsive core if enabled (A_rep > 0)
            if (state%A_rep > 0.0d0) then
                call eval_repulsive_core(r, state%A_rep, state%B_rep, &
                                         V_rep, dVdr_rep)
                V = V + V_rep
                dVdr = dVdr + dVdr_rep
            end if

            Energy = Energy + V * f

            ! Gradient: chain rule dE/dr * dr/d{xyz}
            if (r > 1.0d-12) then
                dE_dr = dVdr * f + V * dfdr
                dEdx = dEdx + dE_dr * dx / r
                dEdy = dEdy + dE_dr * dy / r
                dEdz = dEdz + dE_dr * dz / r
            end if
        end do
    end subroutine calc_rst_energy_grad

    ! ================================================================
    ! Read RST site-organized data file
    !
    ! Format:
    !   # SITE: <name>
    !   # x: <val>  y: <val>
    !   z  weight  energy
    !   ...
    ! ================================================================
    subroutine read_rst_data(filename, sites, total_pts)
        character(len=*), intent(in) :: filename
        type(rst_site_data), intent(out), allocatable :: sites(:)
        integer, intent(out) :: total_pts

        integer :: unit, ios, isite, i, nlines
        integer :: count_sites(4)   ! line count per site
        character(len=200) :: line
        character(len=64) :: key, value
        real(8) :: a_lat, z_val, w_val, e_val, x_tmp, y_tmp
        integer :: space_idx

        a_lat = default_a

        ! --- First pass: count lines per site ---
        allocate(sites(4))
        do isite = 1, 4
            sites(isite)%site_name = trim(site_names(isite))
            sites(isite)%npts = 0
            sites(isite)%x = 0.0d0
            sites(isite)%y = 0.0d0
        end do
        count_sites = 0
        isite = 0

        open(newunit=unit, file=trim(filename), status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Cannot open RST data file: ", trim(filename)
            stop
        end if

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle

            if (line(1:1) == '#') then
                ! Check for site header
                if (index(line, 'SITE:') > 0) then
                    ! Find which site
                    do isite = 1, 4
                        if (index(line, trim(site_names(isite))) > 0) then
                            exit
                        end if
                    end do
                    if (isite > 4) isite = 0  ! unknown site
                end if
                ! Check for x: y: coordinates
                if (isite >= 1 .and. isite <= 4) then
                    if (index(line, 'x:') > 0 .or. index(line, 'X:') > 0) then
                        space_idx = index(line, 'x:')
                        if (space_idx == 0) space_idx = index(line, 'X:')
                        read(line(space_idx+2:), *, iostat=ios) x_tmp
                        if (ios == 0) sites(isite)%x = x_tmp
                    end if
                    if (index(line, 'y:') > 0 .or. index(line, 'Y:') > 0) then
                        space_idx = index(line, 'y:')
                        if (space_idx == 0) space_idx = index(line, 'Y:')
                        read(line(space_idx+2:), *, iostat=ios) y_tmp
                        if (ios == 0) sites(isite)%y = y_tmp
                    end if
                end if
                cycle
            end if

            ! Data line: z weight energy
            if (isite >= 1 .and. isite <= 4) then
                read(line, *, iostat=ios) z_val, w_val, e_val
                if (ios == 0) count_sites(isite) = count_sites(isite) + 1
            end if
        end do
        close(unit)

        ! --- Allocate and fill ---
        total_pts = 0
        do isite = 1, 4
            sites(isite)%npts = count_sites(isite)
            if (count_sites(isite) > 0) then
                allocate(sites(isite)%z(count_sites(isite)))
                allocate(sites(isite)%wgt(count_sites(isite)))
                allocate(sites(isite)%E(count_sites(isite)))
            end if
            total_pts = total_pts + count_sites(isite)
        end do

        ! Fill default (x,y) for sites that weren't specified in the file
        do isite = 1, 4
            if (abs(sites(isite)%x) < 1.0d-10 .and. &
                abs(sites(isite)%y) < 1.0d-10 .and. isite > 1) then
                call get_site_xy(isite, a_lat, sites(isite)%x, &
                                 sites(isite)%y, sites(isite)%site_name)
            end if
        end do
        ! Top site explicit (0,0) is correct by default

        ! --- Second pass: read data ---
        count_sites = 0
        isite = 0

        open(newunit=unit, file=trim(filename), status='old', iostat=ios)

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle

            if (line(1:1) == '#') then
                if (index(line, 'SITE:') > 0) then
                    do isite = 1, 4
                        if (index(line, trim(site_names(isite))) > 0) exit
                    end do
                    if (isite > 4) isite = 0
                end if
                cycle
            end if

            if (isite >= 1 .and. isite <= 4) then
                read(line, *, iostat=ios) z_val, w_val, e_val
                if (ios == 0) then
                    count_sites(isite) = count_sites(isite) + 1
                    i = count_sites(isite)
                    sites(isite)%z(i) = z_val
                    sites(isite)%wgt(i) = w_val
                    sites(isite)%E(i) = e_val
                end if
            end if
        end do
        close(unit)

        write(*, "(A,A,A)") "  RST data loaded from ", trim(filename), ":"
        do isite = 1, 4
            if (sites(isite)%npts > 0) then
                write(*, "(A,A,A,F8.4,A,F8.4,A,I0,A)") "    ", &
                    trim(sites(isite)%site_name), " (", &
                    sites(isite)%x, ",", sites(isite)%y, "): ", &
                    sites(isite)%npts, " points"
            end if
        end do
        write(*, "(A,I0)") "    Total: ", total_pts
    end subroutine read_rst_data

    ! ================================================================
    ! Deallocate rst_site_data array
    ! ================================================================
    subroutine deallocate_rst_sites(sites)
        type(rst_site_data), intent(inout), allocatable :: sites(:)
        integer :: i
        if (.not. allocated(sites)) return
        do i = 1, size(sites)
            if (allocated(sites(i)%z)) deallocate(sites(i)%z)
            if (allocated(sites(i)%wgt)) deallocate(sites(i)%wgt)
            if (allocated(sites(i)%E)) deallocate(sites(i)%E)
        end do
        deallocate(sites)
    end subroutine deallocate_rst_sites

    ! ================================================================
    ! Standalone RST NN weight file loader
    !
    ! File format:
    !   # RST Pairwise NN PES Weights
    !   # r_scale=...  d_cut=...  r_cut=...  tanh activation
    !   # STATE: <name>
    !   # NH_surf: <n>  NH_inner: <m>
    !   ## Surface_Layer
    !   <3*n+1 values>
    !   ## Inner_Layer
    !   <3*m+1 values>
    !   # E0
    !   <1 value>
    !   ... repeat for each state ...
    !   # EHF_baseline <value>  (optional, trailing line)
    !
    ! Returns: rst(7) — indices 1-6 are diabatic, 7 is GS
    !          have_rst(7) — .true. if state was found
    !          baseline — EHF_baseline_shift value (0.0d0 if not found)
    ! ================================================================
    subroutine load_rst_weights_standalone(filename, rst, have_rst, baseline)
        character(len=*), intent(in) :: filename
        type(rst_state), intent(out) :: rst(7)
        logical, intent(out) :: have_rst(7)
        real(8), intent(out) :: baseline

        character(len=200) :: line
        integer :: unit, ios, idx, ns, ni, n_surf_params, n_inner_params
        integer :: i, is, n_tot, idx_a, idx_b
        character(len=3) :: sname
        real(8) :: r_scale_val, r_cut_val, d_cut_val, A_val, B_val

        have_rst = .false.
        baseline = 0.0d0
        r_scale_val = default_r_scale
        r_cut_val = default_r_cut
        d_cut_val = default_d_cut

        open(newunit=unit, file=filename, status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Cannot open RST NN weight file: ", trim(filename)
            stop
        end if

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)

            ! Parse global scale/cutoff parameters from comment lines
            idx = index(line, 'r_scale=')
            if (idx > 0) read(line(idx+8:), *, iostat=ios) r_scale_val
            idx = index(line, 'r_cut=')
            if (idx > 0) read(line(idx+6:), *, iostat=ios) r_cut_val
            idx = index(line, 'd_cut=')
            if (idx > 0) read(line(idx+6:), *, iostat=ios) d_cut_val

            idx = index(line, '# STATE:')
            if (idx == 0) cycle
            if (idx > 1) cycle   ! only at start of line

            sname = adjustl(line(idx+8:))
            do is = 1, 7
                if (trim(sname) == trim(rst_state_names(is))) exit
            end do
            if (is > 7) then
                write(*,*) "WARNING: Unknown state in RST file: ", trim(sname)
                cycle
            end if

            ! Read NH_surf and NH_inner
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            idx = index(line, 'NH_surf:')
            if (idx == 0) then
                write(*,*) "ERROR: Missing NH_surf for state ", trim(sname)
                close(unit); stop
            end if
            ! Parse "NH_surf: N  NH_inner: M" format
            read(line(idx+8:), *, iostat=ios) ns
            if (ios /= 0) then
                write(*,*) "ERROR: Bad NH_surf for state ", trim(sname)
                close(unit); stop
            end if
            idx = index(line, 'NH_inner:')
            if (idx > 0) then
                read(line(idx+9:), *, iostat=ios) ni
            else
                ni = 0
            end if
            if (ios /= 0 .or. ns <= 0) then
                write(*,*) "ERROR: Bad NH for state ", trim(sname)
                close(unit); stop
            end if

            n_surf_params = 3 * ns + 1
            n_inner_params = 3 * ni + 1
            n_tot = n_surf_params + n_inner_params + 1

            rst(is)%nh_surf = ns
            rst(is)%nh_inner = ni
            rst(is)%r_scale = r_scale_val
            rst(is)%r_cut = r_cut_val
            rst(is)%d_cut = d_cut_val
            if (allocated(rst(is)%w)) deallocate(rst(is)%w)
            allocate(rst(is)%w(n_tot))

            ! Skip comment lines; parse # Repulsive: A=… B=… if present
            do
                read(unit, '(A)') line
                if (index(line, '## Surface_Layer') > 0) exit
                if (index(line, '# Repulsive:') > 0) then
                    idx_a = index(line, 'A=')
                    idx_b = index(line, 'B=')
                    if (idx_a > 0) read(line(idx_a+2:), *, iostat=ios) A_val
                    if (idx_b > 0) read(line(idx_b+2:), *, iostat=ios) B_val
                    if (ios == 0) then
                        rst(is)%A_rep = A_val
                        rst(is)%B_rep = B_val
                    end if
                    cycle
                end if
                if (line(1:1) /= '#') then
                    backspace(unit)
                    exit
                end if
            end do
            do i = 1, n_surf_params
                read(unit, *, iostat=ios) rst(is)%w(i)
                if (ios /= 0) then
                    write(*,*) "ERROR reading Surface_Layer of ", trim(sname), &
                        " i=", i, " ios=", ios
                    close(unit); stop
                end if
            end do

            ! Read Inner_Layer block (may be absent if ni==0)
            if (ni > 0) then
                ! Skip comment lines until we find ## Inner_Layer header
                do
                    read(unit, '(A)') line
                    if (index(line, '## Inner_Layer') > 0) exit
                    if (line(1:1) /= '#') then
                        backspace(unit)
                        exit
                    end if
                end do
                do i = 1, n_inner_params
                    read(unit, *, iostat=ios) rst(is)%w(n_surf_params + i)
                    if (ios /= 0) then
                        write(*,*) "ERROR reading Inner_Layer of ", trim(sname)
                        close(unit); stop
                    end if
                end do
            end if

            ! Read E0 — skip any comment lines until we find "# E0"
            do
                read(unit, '(A)') line
                if (index(line, '# E0') > 0 .or. index(line, '#E0') > 0) exit
                if (line(1:1) /= '#') then
                    ! Direct numeric value without # E0 header
                    backspace(unit)
                    exit
                end if
            end do
            read(unit, *, iostat=ios) rst(is)%w(n_tot)
            if (ios /= 0) then
                write(*,*) "ERROR reading E0 of ", trim(sname)
                close(unit); stop
            end if

            have_rst(is) = .true.
            if (rst(is)%A_rep > 0.0d0) then
                write(*,"(A,A,A,I0,A,I0,A,I0,A,F6.1,A,F5.2,A)") "  Loaded ", trim(sname), &
                    " (NH_surf=", ns, ", NH_inner=", ni, &
                    ", ", n_tot, " params, rep=(A=", rst(is)%A_rep, &
                    " B=", rst(is)%B_rep, "))"
            else
                write(*,"(A,A,A,I0,A,I0,A,I0,A)") "  Loaded ", trim(sname), &
                    " (NH_surf=", ns, ", NH_inner=", ni, &
                    ", ", n_tot, " params)"
            end if
        end do
        close(unit)

        ! --- Second pass: scan for # EHF_baseline line ---
        open(newunit=unit, file=filename, status='old', iostat=ios)
        if (ios == 0) then
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                idx = index(line, 'EHF_baseline')
                if (idx > 0) then
                    read(line(idx+12:), *, iostat=ios) baseline
                    if (ios == 0) then
                        write(*,"(A,F12.6)") "  Loaded EHF_baseline = ", baseline
                    end if
                    exit
                end if
            end do
            close(unit)
        end if

        write(*,"(A,A,A)") "RST NN weights loaded from ", trim(filename)
    end subroutine load_rst_weights_standalone

    ! ================================================================
    ! Standalone RST NN weight file writer (full rewrite)
    ! ================================================================
    subroutine save_rst_weights_standalone(filename, rst, have_rst, baseline)
        character(len=*), intent(in) :: filename
        type(rst_state), intent(in) :: rst(7)
        logical, intent(in) :: have_rst(7)
        real(8), intent(in) :: baseline

        integer :: unit, is, i, n_surf_params, n_inner_params, n_tot

        open(newunit=unit, file=filename, status='replace')
        write(unit, '(A)') "# RST Pairwise NN PES Weights"
        write(unit, '(A,F5.2,A,F5.2,A,F5.2,A)') &
            "# r_scale=", default_r_scale, "  d_cut=", default_d_cut, &
            "  r_cut=", default_r_cut, "  tanh activation"
        write(unit, '(A)') "# Layout per state: # STATE: <name> / " &
            // "# NH_surf: <n>  NH_inner: <m> / " &
            // "## Surface_Layer / ## Inner_Layer / # E0"

        do is = 1, 7
            if (.not. have_rst(is)) cycle
            if (.not. allocated(rst(is)%w)) cycle

            n_surf_params = 3 * rst(is)%nh_surf + 1
            n_inner_params = 3 * rst(is)%nh_inner + 1
            n_tot = n_surf_params + n_inner_params + 1

            write(unit, *)
            write(unit, '(A,A)') "# STATE: ", trim(rst_state_names(is))
            write(unit, '(A,I0,A,I0)') "# NH_surf: ", rst(is)%nh_surf, &
                "  NH_inner: ", rst(is)%nh_inner

            write(unit, '(A)') "## Surface_Layer"
            do i = 1, n_surf_params
                write(unit, '(ES24.16)') rst(is)%w(i)
            end do

            if (rst(is)%nh_inner > 0) then
                write(unit, '(A)') "## Inner_Layer"
                do i = 1, n_inner_params
                    write(unit, '(ES24.16)') rst(is)%w(n_surf_params + i)
                end do
            end if

            write(unit, '(A)') "# E0"
            write(unit, '(ES24.16)') rst(is)%w(n_tot)
        end do

        if (abs(baseline) > 1.0d-12) then
            write(unit, '(A,ES24.16)') "# EHF_baseline ", baseline
        end if

        close(unit)
        write(*,*) "RST NN weights saved to ", trim(filename)
    end subroutine save_rst_weights_standalone

end module rst_pes
