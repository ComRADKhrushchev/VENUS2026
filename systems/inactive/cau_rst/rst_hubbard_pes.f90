!=====================================================================
! rst_hubbard_pes.f90 - vendored RST pairwise-NN engine (Hubbard direct
!   mode), extracted VERBATIM from C@Au_TDHF/2D src_NN/interface_RST_Hubbard.f90
!   (module rst_hubbard_pes). Zero external dependencies.
! DEVIATIONS (VENUS2026 systems/cau_rst, 2026-09-21): none in this file -
!   the interface-level deviations (slab syncing, restored gradient
!   chains, debug-flag removal) live in cau_rst_pot.f90 and are listed
!   in its header.
!=====================================================================
module rst_hubbard_pes
    implicit none

    ! --- math constants ---
    real(8), parameter :: pi_val = 3.14159265358979323846d0    ! renamed to avoid Venus COMMON/CONSTN/PI conflict
    real(8), parameter :: sqrt3 = 1.7320508075688772d0
    real(8), parameter :: sqrt2_over_3 = 0.8164965809277260d0
    real(8), parameter :: cons = 2.95d0

    private :: pi_val  ! prevent name clash with Venus COMMON blocks

    ! --- state name mapping (Hubbard direct) ---
    character(len=3), parameter :: state_names(6) = &
        ["h  ", "U  ", "J  ", "AN ", "Vcp", "GS "]

    ! ===== Types =====

    type :: rst_slab_data
        integer :: n_atoms = 0
        integer :: n_surf  = 0
        real(8), allocatable :: coords(:,:)
        integer, allocatable :: layer_id(:)
        real(8) :: box_lx = 36.0d0
        real(8) :: box_ly = 36.0d0
    end type rst_slab_data

    type :: rst_state
        integer  :: nh_surf  = 0
        integer  :: nh_inner = 0
        real(8)  :: r_cut   = 6.0d0
        real(8)  :: d_cut   = 0.5d0
        real(8)  :: r_scale = 3.0d0
        real(8)  :: rep_A   = 0.0d0
        real(8)  :: rep_B   = 3.0d0
        logical  :: use_unified = .false.
        real(8), allocatable :: w(:)
    end type rst_state

    public :: rst_state, rst_slab_data
    public :: calc_rst_energy, calc_rst_energy_grad
    public :: load_rst_weights_standalone, save_rst_weights_standalone
    public :: load_slab_xyz

contains

!==============================================================================
! 1. Math utilities
!==============================================================================

    pure real(8) function softplus(x)
        real(8), intent(in) :: x
        if (x > 50.0d0) then
            softplus = x
        else if (x < -50.0d0) then
            softplus = 0.0d0
        else
            softplus = log(1.0d0 + exp(x))
        end if
    end function softplus

    pure real(8) function sigmoid(x)
        real(8), intent(in) :: x
        if (x > 50.0d0) then
            sigmoid = 1.0d0
        else if (x < -50.0d0) then
            sigmoid = 0.0d0
        else
            sigmoid = 1.0d0 / (1.0d0 + exp(-x))
        end if
    end function sigmoid

    subroutine smooth_cutoff(r, r_cut, d_cut, f, dfdr)
        real(8), intent(in)  :: r, r_cut, d_cut
        real(8), intent(out) :: f, dfdr
        real(8) :: t, r_inner
        r_inner = r_cut - d_cut
        if (r <= r_inner) then
            f = 1.0d0; dfdr = 0.0d0
        else if (r >= r_cut) then
            f = 0.0d0; dfdr = 0.0d0
        else
            t = (r - r_inner) / d_cut
            f = 0.5d0 * (1.0d0 + cos(pi_val * t))
            dfdr = -0.5d0 * pi_val / d_cut * sin(pi_val * t)
        end if
    end subroutine smooth_cutoff

    subroutine pbc_min_image(dx, dy, box_lx, box_ly)
        real(8), intent(inout) :: dx, dy
        real(8), intent(in)    :: box_lx, box_ly
        dx = dx - box_lx * nint(dx / box_lx)
        dy = dy - box_ly * nint(dy / box_ly)
    end subroutine pbc_min_image

!==============================================================================
! 2. Born-Mayer repulsive potential
!==============================================================================

    subroutine born_mayer(r, A, B, V, dVdr)
        real(8), intent(in)  :: r, A, B
        real(8), intent(out) :: V, dVdr
        real(8) :: er
        if (abs(A) < 1.0d-15) then
            V = 0.0d0; dVdr = 0.0d0; return
        end if
        er = exp(-B * r)
        V = A * er
        dVdr = -A * B * er
    end subroutine born_mayer

    subroutine born_mayer_grad(r, A, B, V, dVdr, d2Vdr2)
        real(8), intent(in)  :: r, A, B
        real(8), intent(out) :: V, dVdr, d2Vdr2
        real(8) :: er
        if (abs(A) < 1.0d-15) then
            V = 0.0d0; dVdr = 0.0d0; d2Vdr2 = 0.0d0; return
        end if
        er = exp(-B * r)
        V = A * er
        dVdr = -A * B * er
        d2Vdr2 = A * B * B * er
    end subroutine born_mayer_grad

!==============================================================================
! 3. 1D pairwise NN evaluation
!==============================================================================

    subroutine eval_1d_nn_dist(r, weights, nh, r_scale, V, dVdr, d2Vdr2)
        real(8), intent(in)  :: r
        real(8), intent(in)  :: weights(:)
        integer, intent(in)  :: nh
        real(8), intent(in)  :: r_scale
        real(8), intent(out) :: V, dVdr
        real(8), intent(out), optional :: d2Vdr2
        integer :: i
        real(8) :: w1, b1, w2, h, rs, dh_dr
        rs = r / r_scale
        V = weights(3 * nh + 1)
        dVdr = 0.0d0
        if (present(d2Vdr2)) d2Vdr2 = 0.0d0
        do i = 1, nh
            w1 = weights(i)
            b1 = weights(nh + i)
            w2 = weights(2 * nh + i)
            h = tanh(w1 * rs + b1)
            V = V + w2 * h
            dh_dr = (1.0d0 - h * h) * w1 / r_scale
            dVdr = dVdr + w2 * dh_dr
            if (present(d2Vdr2)) then
                d2Vdr2 = d2Vdr2 - 2.0d0 * h * (1.0d0 - h * h) * (w1 / r_scale)**2 * w2
            end if
        end do
    end subroutine eval_1d_nn_dist

!==============================================================================
! 4. Pair potential composition (NN + Born-Mayer)
!==============================================================================

    subroutine eval_repulsive(r, A, B, V_rep, dV_rep, d2V_rep)
        real(8), intent(in)  :: r, A, B
        real(8), intent(out) :: V_rep, dV_rep
        real(8), intent(out), optional :: d2V_rep
        real(8) :: d2
        if (present(d2V_rep)) then
            call born_mayer_grad(r, A, B, V_rep, dV_rep, d2V_rep)
        else
            call born_mayer_grad(r, A, B, V_rep, dV_rep, d2)
        end if
    end subroutine eval_repulsive

    subroutine eval_pair_potential(r, weights, nh, r_scale, rep_A, rep_B, V, dVdr, d2Vdr2)
        real(8), intent(in)  :: r, weights(:), r_scale, rep_A, rep_B
        integer, intent(in)  :: nh
        real(8), intent(out) :: V, dVdr
        real(8), intent(out), optional :: d2Vdr2
        real(8) :: V_nn, dV_nn, d2V_nn, V_rep, dV_rep, d2V_rep
        call eval_1d_nn_dist(r, weights, nh, r_scale, V_nn, dV_nn, d2V_nn)
        call eval_repulsive(r, rep_A, rep_B, V_rep, dV_rep, d2V_rep)
        V = V_nn + V_rep
        dVdr = dV_nn + dV_rep
        if (present(d2Vdr2)) d2Vdr2 = d2V_nn + d2V_rep
    end subroutine eval_pair_potential

!==============================================================================
! 5. RST energy & gradient
!==============================================================================

    subroutine calc_rst_energy(xyz, slab, state, Energy, no_repulsive)
        real(8), intent(in)              :: xyz(3)
        type(rst_slab_data), intent(in)  :: slab
        type(rst_state), intent(in)      :: state
        real(8), intent(out)             :: Energy
        logical, intent(in), optional    :: no_repulsive
        real(8) :: dEdx, dEdy, dEdz
        call calc_rst_energy_grad(xyz, slab, state, Energy, dEdx, dEdy, dEdz, no_repulsive)
    end subroutine calc_rst_energy

    subroutine calc_rst_energy_grad(xyz, slab, state, Energy, dEdx, dEdy, dEdz, no_repulsive, dEdR)
        real(8), intent(in)              :: xyz(3)
        type(rst_slab_data), intent(in)  :: slab
        type(rst_state), intent(in)      :: state
        real(8), intent(out)             :: Energy, dEdx, dEdy, dEdz
        logical, intent(in), optional    :: no_repulsive
        real(8), intent(out), optional   :: dEdR(3, slab%n_atoms)
        logical  :: skip_repulsive
        integer  :: i, n_surf_params, n_inner_params
        real(8)  :: dx, dy, dz, r, f, dfdr, V, dVdr, dE_dr, E0
        real(8)  :: w_surf_vals(3 * state%nh_surf + 1)
        real(8)  :: w_inner_vals(3 * state%nh_inner + 1)

        skip_repulsive = .false.
        if (present(no_repulsive)) skip_repulsive = no_repulsive

        n_surf_params  = 3 * state%nh_surf  + 1
        n_inner_params = 3 * state%nh_inner + 1
        if (state%nh_surf  > 0) w_surf_vals  = state%w(1:n_surf_params)
        if (state%nh_inner > 0) w_inner_vals = state%w(n_surf_params+1 : n_surf_params+n_inner_params)
        E0 = state%w(n_surf_params + n_inner_params + 1)

        Energy = E0
        dEdx = 0.0d0; dEdy = 0.0d0; dEdz = 0.0d0
        if (present(dEdR)) dEdR = 0.0d0

        do i = 1, slab%n_atoms
            dx = xyz(1) - slab%coords(1, i)
            dy = xyz(2) - slab%coords(2, i)
            dz = xyz(3) - slab%coords(3, i)
            call pbc_min_image(dx, dy, slab%box_lx, slab%box_ly)
            r = sqrt(dx*dx + dy*dy + dz*dz)
            if (r > state%r_cut) cycle
            call smooth_cutoff(r, state%r_cut, state%d_cut, f, dfdr)

            if (slab%layer_id(i) == 1 .and. state%nh_surf > 0) then
                if (skip_repulsive) then
                    call eval_1d_nn_dist(r, w_surf_vals, state%nh_surf, state%r_scale, V, dVdr)
                else
                    call eval_pair_potential(r, w_surf_vals, state%nh_surf, &
                        state%r_scale, state%rep_A, state%rep_B, V, dVdr)
                end if
            else if (slab%layer_id(i) == 2) then
                if (state%use_unified .and. state%nh_surf > 0) then
                    if (skip_repulsive) then
                        call eval_1d_nn_dist(r, w_surf_vals, state%nh_surf, state%r_scale, V, dVdr)
                    else
                        call eval_pair_potential(r, w_surf_vals, state%nh_surf, &
                            state%r_scale, state%rep_A, state%rep_B, V, dVdr)
                    end if
                else if (state%nh_inner > 0) then
                    if (skip_repulsive) then
                        call eval_1d_nn_dist(r, w_inner_vals, state%nh_inner, state%r_scale, V, dVdr)
                    else
                        call eval_pair_potential(r, w_inner_vals, state%nh_inner, &
                            state%r_scale, state%rep_A, state%rep_B, V, dVdr)
                    end if
                else
                    V = 0.0d0; dVdr = 0.0d0
                end if
            else
                V = 0.0d0; dVdr = 0.0d0
            end if

            Energy = Energy + V * f
            if (r > 1.0d-12) then
                dE_dr = dVdr * f + V * dfdr
                dEdx = dEdx + dE_dr * dx / r
                dEdy = dEdy + dE_dr * dy / r
                dEdz = dEdz + dE_dr * dz / r
                if (present(dEdR)) then
                    dEdR(1, i) = dEdR(1, i) + dE_dr * dx / r
                    dEdR(2, i) = dEdR(2, i) + dE_dr * dy / r
                    dEdR(3, i) = dEdR(3, i) + dE_dr * dz / r
                end if
            end if
        end do
    end subroutine calc_rst_energy_grad

!==============================================================================
! 6. Weight I/O
!==============================================================================

    subroutine parse_repulsive_line(line, rep_A, rep_B, use_unified)
        character(len=*), intent(in)    :: line
        real(8),          intent(inout) :: rep_A, rep_B
        logical,          intent(inout) :: use_unified
        integer :: idx_a, idx_b, idx_u, ios
        character(len=20) :: ustr
        idx_a = index(line, 'A=')
        if (idx_a > 0) read(line(idx_a+2:), *, iostat=ios) rep_A
        idx_b = index(line, 'B=')
        if (idx_b > 0) read(line(idx_b+2:), *, iostat=ios) rep_B
        idx_u = index(line, 'unified=')
        if (idx_u > 0) then
            read(line(idx_u+8:), *, iostat=ios) ustr
            if (ios == 0) use_unified = (trim(ustr) == 'true' .or. trim(ustr) == 'True')
        end if
    end subroutine parse_repulsive_line

    subroutine load_rst_weights_standalone(filename, rst, have_rst, baseline)
        character(len=*), intent(in)  :: filename
        type(rst_state),  intent(out) :: rst(6)
        logical,          intent(out) :: have_rst(6)
        real(8),          intent(out) :: baseline
        character(len=200) :: line
        integer :: unit, ios, idx, ns, ni, n_surf_params, n_inner_params, i, is, n_tot
        character(len=3) :: sname
        real(8) :: r_scale_val, r_cut_val, d_cut_val

        have_rst = .false.; baseline = 0.0d0
        r_scale_val = 3.0d0; r_cut_val = 6.0d0; d_cut_val = 0.5d0

        open(newunit=unit, file=filename, status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Cannot open RST NN weight file: ", trim(filename); stop
        end if

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            idx = index(line, 'r_scale=')
            if (idx > 0) read(line(idx+8:), *, iostat=ios) r_scale_val
            idx = index(line, 'r_cut=')
            if (idx > 0) read(line(idx+6:), *, iostat=ios) r_cut_val
            idx = index(line, 'd_cut=')
            if (idx > 0) read(line(idx+6:), *, iostat=ios) d_cut_val
            idx = index(line, '# STATE:')
            if (idx /= 1) cycle

            sname = adjustl(line(idx+8:))
            do is = 1, 6
                if (trim(sname) == trim(state_names(is))) exit
            end do
            if (is > 6) cycle

            read(unit, '(A)', iostat=ios) line
            idx = index(line, 'NH_surf:')
            if (idx == 0) then
                write(*,*) "ERROR: Missing NH_surf for state ", trim(sname); close(unit); stop
            end if
            read(line(idx+8:), *, iostat=ios) ns
            idx = index(line, 'NH_inner:')
            if (idx > 0) then; read(line(idx+9:), *, iostat=ios) ni
            else; ni = 0; end if

            n_surf_params = 3*ns + 1; n_inner_params = 3*ni + 1
            n_tot = n_surf_params + n_inner_params + 1
            rst(is)%nh_surf = ns; rst(is)%nh_inner = ni
            rst(is)%r_scale = r_scale_val; rst(is)%r_cut = r_cut_val; rst(is)%d_cut = d_cut_val
            if (allocated(rst(is)%w)) deallocate(rst(is)%w)
            allocate(rst(is)%w(n_tot))

            ! Skip optional per-state r_cut/d_cut/r_scale comment line (= or :)
            read(unit, '(A)') line
            do while (index(line, 'r_cut') > 0 .or. index(line, 'd_cut') > 0 .or. index(line, 'r_scale') > 0)
                read(unit, '(A)') line
            end do
            if (index(line, 'Repulsive:') > 0 .and. index(line, '#') > 0) then
                call parse_repulsive_line(line, rst(is)%rep_A, rst(is)%rep_B, rst(is)%use_unified)
                read(unit, '(A)') line
            end if
            do i = 1, n_surf_params
                read(unit, *, iostat=ios) rst(is)%w(i)
            end do
            if (ni > 0) then
                read(unit, '(A)') line
                do i = 1, n_inner_params
                    read(unit, *, iostat=ios) rst(is)%w(n_surf_params + i)
                end do
            end if
            read(unit, '(A)') line
            read(unit, *, iostat=ios) rst(is)%w(n_tot)
            have_rst(is) = .true.
            write(*,"(A,A,A,I0,A,I0,A,I0,A)") "  Loaded ", trim(sname), &
                " (NH_surf=", ns, ", NH_inner=", ni, ", ", n_tot, " params)"
        end do
        close(unit)

        open(newunit=unit, file=filename, status='old', iostat=ios)
        if (ios == 0) then
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                idx = index(line, 'EHF_baseline')
                if (idx > 0) then
                    read(line(idx+12:), *, iostat=ios) baseline
                    if (ios == 0) write(*,"(A,F12.6)") "  Loaded EHF_baseline = ", baseline
                    exit
                end if
            end do
            close(unit)
        end if
        write(*,"(A,A,A)") "RST NN weights loaded from ", trim(filename)
    end subroutine load_rst_weights_standalone

    subroutine save_rst_weights_standalone(filename, rst, have_rst, baseline)
        character(len=*), intent(in) :: filename
        type(rst_state),  intent(in) :: rst(6)
        logical,          intent(in) :: have_rst(6)
        real(8),          intent(in) :: baseline
        integer :: unit, is, i, n_surf_params, n_inner_params, n_tot

        open(newunit=unit, file=filename, status='replace')
        write(unit, '(A)') "# RST Pairwise NN PES Weights"
        write(unit, '(A)') "# r_scale=3.0  d_cut=0.5  r_cut=6.0  tanh activation"
        write(unit, '(A)') "# HUBBARD_MODE: direct"
        write(unit, '(A)') "# Layout per state: # STATE: <name> / " &
            // "# NH_surf: <n>  NH_inner: <m> / " &
            // "## Surface_Layer / ## Inner_Layer / # E0"
        do is = 1, 6
            if (.not. have_rst(is)) cycle
            if (.not. allocated(rst(is)%w)) cycle
            n_surf_params  = 3*rst(is)%nh_surf  + 1
            n_inner_params = 3*rst(is)%nh_inner + 1
            n_tot = n_surf_params + n_inner_params + 1
            write(unit, *)
            write(unit, '(A,A)') "# STATE: ", trim(state_names(is))
            write(unit, '(A,I0,A,I0)') "# NH_surf: ", rst(is)%nh_surf, "  NH_inner: ", rst(is)%nh_inner
            if (rst(is)%rep_A > 1.0d-12 .or. rst(is)%use_unified) then
                if (rst(is)%use_unified) then
                    write(unit, '(A,F10.4,A,F8.4,A)') "# Repulsive: A=", rst(is)%rep_A, &
                        " B=", rst(is)%rep_B, " unified=true"
                else
                    write(unit, '(A,F10.4,A,F8.4,A)') "# Repulsive: A=", rst(is)%rep_A, &
                        " B=", rst(is)%rep_B, " unified=false"
                end if
            end if
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
        if (abs(baseline) > 1.0d-12) write(unit, '(A,ES24.16)') "# EHF_baseline ", baseline
        close(unit)
        write(*,*) "RST NN weights saved to ", trim(filename)
    end subroutine save_rst_weights_standalone

!==============================================================================
! 7. Slab I/O
!==============================================================================

    subroutine load_slab_xyz(filename, slab)
        character(len=*), intent(in)  :: filename
        type(rst_slab_data), intent(out) :: slab
        integer :: unit, ios, n, i
        real(8) :: x, y, z, z_max, d_layer_est
        character(len=2) :: elem
        character(len=200) :: line

        open(newunit=unit, file=trim(filename), status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Cannot open slab file: ", trim(filename); stop
        end if
        read(unit, *, iostat=ios) n
        if (ios /= 0) then
            write(*,*) "ERROR: Failed to read atom count from ", trim(filename)
            close(unit); stop
        end if
        read(unit, '(A)', iostat=ios) line

        slab%n_atoms = n
        if (allocated(slab%coords))  deallocate(slab%coords)
        if (allocated(slab%layer_id)) deallocate(slab%layer_id)
        allocate(slab%coords(3, n))
        allocate(slab%layer_id(n))

        do i = 1, n
            read(unit, *, iostat=ios) elem, x, y, z
            if (ios /= 0) then
                write(*,*) "ERROR: Failed to read atom ", i; close(unit); stop
            end if
            slab%coords(1, i) = x
            slab%coords(2, i) = y
            slab%coords(3, i) = z
        end do
        close(unit)

        slab%box_lx = maxval(slab%coords(1, :)) - minval(slab%coords(1, :)) + 6.0d0
        slab%box_ly = maxval(slab%coords(2, :)) - minval(slab%coords(2, :)) + 6.0d0

        z_max = maxval(slab%coords(3, :))
        d_layer_est = cons * sqrt2_over_3
        slab%n_surf = 0
        do i = 1, n
            if (slab%coords(3, i) > z_max - 0.5d0 * d_layer_est) then
                slab%layer_id(i) = 1
                slab%n_surf = slab%n_surf + 1
            else
                slab%layer_id(i) = 2
            end if
        end do

        write(*,"(A,A,A,I0,A,I0,A)") "  Slab loaded from ", trim(filename), &
            ": ", slab%n_atoms, " atoms, ", slab%n_surf, " surface"
    end subroutine load_slab_xyz

end module rst_hubbard_pes
