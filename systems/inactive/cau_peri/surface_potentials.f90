!=====================================================================
! surface_potentials.f90 - vendored 2D analytic PES module (periodic
!   series), extracted from C@Au_TDHF/2D src_NN/surface_potentials.f90.
! GRADIENT FIXES (VENUS2026 systems/cau_peri, 2026-09-21): the analytic
!   gradient code of Calc_POT_spec_and_Grad was NOT the gradient of its
!   own energy (FD-probed: z-component sign flipped whenever the image
!   plane p5 sat above z, in-plane components off by the p5 chain) -
!   the rho = |z| - p5 convention needs d(rho)/dz = sign(z) and
!   d(rho)/dp5 = -1 with no sign factor. Fixed at the two sites flagged
!   "DEVIATION FIX" below. Criterion (c) of check_cau_peri verifies.
! QUARANTINE (VENUS2026 systems/cau_peri, 2026-09-21): the two host
!   `use` lines (venus_params / venus_data) are dropped and their only
!   consumed symbols - NDA3 (workspace dimension) and NN_POT_TYPES -
!   are defined locally with the host values (NDA = 1000 -> NDA3 = 3000;
!   NN types = 7); the save_prev_matrices GWRITE-debug tail (the sole
!   IESHParas consumer) is removed. Everything else is verbatim.
!=====================================================================
!==============================================================================
!  surface_potentials — 解析势能面模块
!
!  为 C@Au(111) TDHF-FSSH 提供所有解析势能函数及其梯度。
!  势能面依赖三个坐标:
!    (u,v) — 表面六角晶格坐标 (reduced coordinates on Au(111))
!    z    — 垂直于表面的高度
!
!  (u,v) 与 Cartesian (x,y) 的转换:
!    u = x - y/√3
!    v = 2*y/√3
!
!  势能函数:
!    Calc_POT0_and_Grad   — Morse 势 (用于构造 Hubbard U, J, h_til)
!    Calc_coup_and_Grad   — 耦合 V_k(z) = A/(1 + B*exp(C*z))
!    Calc_POT_spec_and_Grad — 扩展 Rydberg 势 (基态 PES)
!==============================================================================
module surface_potentials
    implicit real*8 (a-h,o-z)

    ! quarantine locals (host values: NDA = 1000 -> NDA3 = 3000; NN types = 7)
    integer, parameter :: NDA3 = 3000
    integer, parameter :: NN_POT_TYPES = 7
    ! the NN type indices (host constants in the original venus_params.f90;
    ! never vendored - under implicit typing they silently became implicit
    ! integers of value 0, which zeroed every section lookup and made the
    ! loader return an empty NN model with use_nn_pes = .true. - restored
    ! verbatim from the original numbering, 2026-09-29)
    integer, parameter :: NN_MAX_HIDDEN = 20
    integer, parameter :: NN_E0 = 1, NN_1D = 2, NN_3P = 3, NN_mC = 4
    integer, parameter :: NN_pC = 5, NN_COUP = 6, NN_Eh = 7

    ! =========================================================================
    ! 物理常数与晶格参数
    ! =========================================================================
    real*8, parameter :: pii = 3.14159265358979323846d0
    real*8, parameter :: sq3 = 1.7320508075688772d0
    real*8 :: cons = 2.95d0          ! Au(111) 晶格常数 (bohr)

    ! =========================================================================
    ! 势能参数数组 (从 Analytic_Potential.txt 读入)
    ! =========================================================================
    ! parm_nas(1:16)  — 基态扩展 Rydberg 势: 5 个 effective par + E_0
    ! parm_coup(1:10) — 耦合 V_k: 3 组 effective par (A,B,C) + ...
    ! parm_1D(1:10)   — ¹D 态 Morse: (D_top,D_bri,D_hollow; a_top,a_bri,a_hollow;
    !                    z0_top,z0_bri,z0_hollow; E_0)
    ! parm_3P(1:10)   — ³P 态 Morse
    ! parm_mC(1:10)   — 负离子态 (−C) Morse
    ! parm_pC(1:10)   — 正离子态 (+C) Morse
    real*8 :: parm_nas(16), parm_coup(10)
    real*8 :: parm_1D(10), parm_3P(10), parm_mC(10), parm_pC(10)

    ! =========================================================================
    ! 表面谐振子力常数 (kcal/mol/bohr²)
    ! =========================================================================
    real*8 :: k_surf(3)     ! 表面振子 (atom1-atom2耦合)
    real*8 :: k_ghost(3)    ! ghost自谐振 (atom2自洽势阱)
    real*8 :: Q0_ghost(3)   ! ghost平衡位置 (Å)

    ! =========================================================================
    ! 受力记录数组 (用于能量守恒 debug)
    ! =========================================================================
    real(8), allocatable :: force_e0(:), force_core(:), force_mfv(:)
    real(8), allocatable :: pforce_e0(:), pforce_core(:), pforce_mfv(:)
    real*8 :: work_e0, work_core, work_mf
    real*8 :: work_e0_sub, work_core_sub, work_mf_sub  ! substep work accumulators
    real*8 :: pE_e0, pE_core, pE_mf,p_H
    ! 10-term decomposition: previous-step matrices (for GWRITE [3]-[6])
    real*8, allocatable :: old_Hcore(:,:,:), old_W(:,:,:), old_J(:,:,:)
    real*8, allocatable :: old_rhom(:,:,:)
    real*8, allocatable :: cur_Hcore(:,:,:), cur_W(:,:,:), cur_J(:,:,:)
    real*8, allocatable :: cur_rhom(:,:,:)
    ! Energy drift correction: errors from rho-evolution in previous step
    real*8 :: E_err_core_prev, E_err_mf_prev, E_err_total_prev
    real*8 :: dQ_prev(NDA3)  ! displacement from previous step

    ! =========================================================================
    ! 晶场劈裂 — 供 dH_orb 使用的梯度 (kcal/mol/bohr)
    ! =========================================================================
    real*8 :: cf_ddx_saved, cf_ddy_saved, cf_ddz_saved

    ! =========================================================================
    ! Neural Network 势能
    ! =========================================================================
    ! Per-type NN architecture metadata (set by load_all_nn_weights)
    integer :: nn_nhidden(NN_POT_TYPES)       ! hidden neurons per type (NH)
    integer :: nn_nw(NN_POT_TYPES)            ! weights per site (3*NH+1) or per type for RST
    integer :: nn_total(NN_POT_TYPES)         ! total weights per type (3*nw or nw for RST)
    integer :: nn_offset(NN_POT_TYPES)        ! start index (1-based) in nn_weights_all
    real*8 :: nn_e0_offset(NN_POT_TYPES) = 0.0d0  ! per-type E0 energy offset (eV)
    real*8, allocatable :: nn_weights_all(:)  ! flat weight array, allocated at load time
    logical :: use_nn_pes = .false.           ! .true. = global NN PES enabled
    logical :: use_eh_correction = .false.   ! .true. = add E_h to orbital energies
    logical :: use_rst_nn = .false.           ! .true. = RST pairwise format (Surf+Inner layers)
    integer :: nn_nhidden_surf(NN_POT_TYPES)  ! RST: NH_surf per type
    integer :: nn_nhidden_inner(NN_POT_TYPES) ! RST: NH_inner per type

contains

!==============================================================================
!  effective_par — 表面参数 Fourier 插值 (标量版)
!
!  在 Au(111) 六角晶格上，用 top / bridge / hollow 三个高对称点的
!  已知参数值，通过二阶 Fourier 级数插值出任意 (u,v) 点的有效参数值。
!
!  Fourier 基函数:
!    f1(u,v) = cos(2πu/a) + cos(2πv/a) + cos(2π(u+v)/a)
!    f2(u,v) = cos(4πu/a) + cos(4πv/a) + cos(4π(u+v)/a)
!    where a = cons = 2.95 bohr (lattice constant)
!
!  p(u,v) = c0 + c1*f1 + c2*f2
!
!  c0, c1, c2 由三个已知点确定:
!    Top site (u=0,v=0):        f1=3,  f2=3  → p_top  = c0 + 3c1 + 3c2
!    Bridge site:               已知值          → p_bri
!    Hollow site:               已知值          → p_hollow = p_cen
!
!  解得:
!    c1 = (p_top - p_bri) / 4
!    c0 = (p_top + 2*p_cen) / 3
!    c2 = (p_top - 3*c1 - c0) / 3
!==============================================================================
    subroutine effective_par(u, v, p_top, p_bri, p_cen, p_m)
        real(8), intent(in)  :: p_top, p_bri, p_cen  ! top/bridge/hollow 已知值
        real(8), intent(in)  :: u, v                    ! 表面坐标
        real(8), intent(out) :: p_m                     ! 插值结果
        real(8) :: fir, snd, c0, c1, c2

        ! 计算 Fourier 基函数在 (u,v) 处的值
        fir = cos(2.0d0*pii*u/cons) + cos(2.0d0*pii*v/cons) &
            + cos(2.0d0*pii*(u+v)/cons)
        snd = cos(4.0d0*pii*u/cons) + cos(4.0d0*pii*v/cons) &
            + cos(4.0d0*pii*(u+v)/cons)

        ! 由三个高对称点值计算 Fourier 系数
        c1 = (p_top - p_bri) / 4.0d0
        c0 = (p_top + 2.0d0 * p_cen) / 3.0d0
        c2 = (p_top - 3.0d0 * c1 - c0) / 3.0d0

        ! 插值
        p_m = c0 + c1 * fir + c2 * snd
    end subroutine

!==============================================================================
!  effective_par_grad — 表面参数 Fourier 插值 (含 Cartesian 梯度)
!
!  与 effective_par 相同的插值逻辑，同时计算 dp/du, dp/dv 并转换为
!  Cartesian 坐标下的 dp/dx, dp/dy。
!
!  链式法则:
!    dp/dx = dp/du * du/dx + dp/dv * dv/dx
!    dp/dy = dp/du * du/dy + dp/dv * dv/dy
!
!  由 u = x - y/√3, v = 2y/√3:
!    du/dx = 1,   du/dy = -1/√3
!    dv/dx = 0,   dv/dy =  2/√3
!
!  因此:
!    dp/dx = dp/du
!    dp/dy = (-1/√3)*dp/du + (2/√3)*dp/dv
!==============================================================================
    subroutine effective_par_grad(u, v, p_top, p_bri, p_cen, p_m, dp_dx, dp_dy)
        real*8, intent(in)  :: p_top, p_bri, p_cen  ! top/bridge/hollow 已知值
        real*8, intent(in)  :: u, v                    ! 表面坐标
        real*8, intent(out) :: p_m                     ! 插值结果
        real*8, intent(out) :: dp_dx, dp_dy            ! Cartesian 梯度

        real*8 :: c0, c1, c2
        real*8 :: arg_u, arg_v, arg_uv
        real*8 :: fir, snd
        real*8 :: dfir_du, dfir_dv, dsnd_du, dsnd_dv
        real*8 :: dfir_dx, dfir_dy, dsnd_dx, dsnd_dy

        ! --- Fourier 系数 (同 effective_par) ---
        c1 = (p_top - p_bri) / 4.0d0
        c0 = (p_top + 2.0d0 * p_cen) / 3.0d0
        c2 = (p_top - 3.0d0 * c1 - c0) / 3.0d0

        ! --- 基函数在 (u,v) 处的值 ---
        arg_u  = 2.0d0 * pii * u / cons
        arg_v  = 2.0d0 * pii * v / cons
        arg_uv = 2.0d0 * pii * (u + v) / cons

        fir = cos(arg_u) + cos(arg_v) + cos(arg_uv)
        snd = cos(2.0d0*arg_u) + cos(2.0d0*arg_v) + cos(2.0d0*arg_uv)

        p_m = c0 + c1 * fir + c2 * snd

        ! --- 基函数对 (u,v) 的偏导数 ---
        ! df1/du = -(2π/a)*[sin(arg_u) + sin(arg_uv)]
        ! df1/dv = -(2π/a)*[sin(arg_v) + sin(arg_uv)]
        dfir_du = -(2.0d0 * pii / cons) * (sin(arg_u) + sin(arg_uv))
        dfir_dv = -(2.0d0 * pii / cons) * (sin(arg_v) + sin(arg_uv))

        ! df2/du = -(4π/a)*[sin(2*arg_u) + sin(2*arg_uv)]
        ! df2/dv = -(4π/a)*[sin(2*arg_v) + sin(2*arg_uv)]
        dsnd_du = -(4.0d0 * pii / cons) * (sin(2.0d0*arg_u) + sin(2.0d0*arg_uv))
        dsnd_dv = -(4.0d0 * pii / cons) * (sin(2.0d0*arg_v) + sin(2.0d0*arg_uv))

        ! --- 转换为 Cartesian 梯度 ---
        ! dp/dx = dp/du * du/dx + dp/dv * dv/dx = dp/du * 1 + dp/dv * 0 = dp/du
        dfir_dx = dfir_du
        dsnd_dx = dsnd_du

        ! dp/dy = dp/du * (-1/√3) + dp/dv * (2/√3)
        dfir_dy = -dfir_du / sq3 + 2.0d0 * dfir_dv / sq3
        dsnd_dy = -dsnd_du / sq3 + 2.0d0 * dsnd_dv / sq3

        dp_dx = c1 * dfir_dx + c2 * dsnd_dx
        dp_dy = c1 * dfir_dy + c2 * dsnd_dy
    end subroutine effective_par_grad

!==============================================================================
!  Calc_crystal_field — 晶体场劈裂 (2px/2py 轨道能级差)
!
!  在 Au(111) 六角晶格上通过 effective_par_grad 进行 Fourier 插值，
!  计算 C 原子 2px/2py 轨道的能级劈裂 (kcal/mol)。
!
!  高对称点取值:
!    Top site (C₃v):    delta = cf_top    (简并, 默认 0)
!    Bridge site (C₂v): delta = cf_bri    (劈裂最大)
!    Hollow site (C₃v): delta = cf_hcp    (简并, 默认 0)
!
!  输出:
!    delta_cf             — 劈裂能 (eV), = h₂ − h₁（参数按 kcal/mol 撰写，
!                           交付处乘 1/23.0605 换算——2026-09-28 单位错误裁决：
!                           原实现把 kcal/mol 参数当 eV 直接交付，劈裂放大 23 倍）
!    ddelta_dx,dy,dz      — Cartesian 梯度 eV/A（同一换算）
!
!  调用后模块变量 cf_ddx_saved / cf_ddy_saved 被更新，供 dH_orb 使用。
!==============================================================================
subroutine Calc_crystal_field(x, y, z, delta_cf, ddelta_dx, ddelta_dy, ddelta_dz)
    real*8, intent(in)  :: x, y, z
    real*8, intent(out) :: delta_cf
    real*8, intent(out) :: ddelta_dx, ddelta_dy, ddelta_dz

    real*8, parameter :: cf_top   = 0.0d0
    real*8, parameter :: cf_bri   = 5.0d0
    real*8, parameter :: cf_hcp   = 0.0d0
    real*8, parameter :: cf_const = 0.01d0
    real*8, parameter :: z_switch = 3.0d0
    real*8, parameter :: z_width  = 1.0d0

    real*8 :: u_loc, v_loc, delta_2d, d2d_dx, d2d_dy
    real*8 :: fz, dfz_dz, tanh_arg, sech2
    real*8 :: kcal_to_ev

    u_loc = x - y / sq3
    v_loc = 2.0d0 * y / sq3

    call effective_par_grad(u_loc, v_loc, cf_top, cf_bri, cf_hcp, &
                            delta_2d, d2d_dx, d2d_dy)

    ! z-decay: fz = 0.5*(1 - tanh((z-z0)/w)),  1 near surface, 0 at large z
    tanh_arg = (z - z_switch) / z_width
    fz       = 0.5d0 * (1.0d0 - tanh(tanh_arg))
    sech2    = 1.0d0 / cosh(tanh_arg)
    sech2    = sech2 * sech2
    dfz_dz   = -0.5d0 / z_width * sech2

    ! kcal/mol 参数 -> eV 交付（立法换算对 23.0605 的倒数；cf_bri=5.0 kcal
    ! = 0.2168 eV 才是本意——用户裁决 2026-09-28，此前未换算致劈裂 23 倍过量）
    kcal_to_ev = 1.0d0 / 23.0605d0
    delta_cf   = (delta_2d + cf_const) * fz * kcal_to_ev
    ddelta_dx  = d2d_dx * fz * kcal_to_ev
    ddelta_dy  = d2d_dy * fz * kcal_to_ev
    ddelta_dz  = (delta_2d + cf_const) * dfz_dz * kcal_to_ev

    cf_ddx_saved = ddelta_dx
    cf_ddy_saved = ddelta_dy
    cf_ddz_saved = ddelta_dz
end subroutine Calc_crystal_field

!==============================================================================
!  Calc_POT0_and_Grad — 普通 Morse 势 (用于 1D/3P/mC/pC 轨道)
!
!  形式: V(z; D,a,z0) = D * [x² − 2x] + E₀
!        其中 x = exp(−a·(z − z0))
!
!  Morse 参数 (D, a, z0) 随 (u,v) 位置通过 effective_par_grad 插值:
!    p_m(1) = D(u,v)   — 阱深     [eV]
!    p_m(2) = a(u,v)   — 宽度参数 [1/bohr]
!    p_m(3) = z0(u,v)  — 平衡位置 [bohr]
!    E_0               — 能量零点 [eV] (parm(10), 不随 (u,v) 变化)
!
!  梯度公式 (链式法则求导):
!    dV/dD       = x² − 2x
!    dV/da       = −2·D·x·(x−1)·(z−z0)
!    dV/dz0      =  2·D·a·x·(x−1)
!    dV/dz       = −2·D·a·x·(x−1)
!
!  (x,y) 方向梯度通过 dE_dp(i) × dp_dx/dp_dy 累积得到。
!==============================================================================
    subroutine Calc_POT0_and_Grad(parm, x, y, z, Energy, dEdx, dEdy, dEdz)
        real*8, intent(in)  :: parm(:)  ! [D_i(3), a_i(3), z0_i(3), E_0]  (10)
        real*8, intent(in)  :: x, y, z
        real*8, intent(out) :: Energy, dEdx, dEdy, dEdz

        real*8 :: u, v
        real*8 :: p_m(3)        ! D(u,v), a(u,v), z0(u,v) at this (x,y)
        real*8 :: dp_dx(3), dp_dy(3)  ! ∂{D,a,z0}/∂x, ∂{D,a,z0}/∂y
        real*8 :: dE_dp(3)      ! ∂V/∂D, ∂V/∂a, ∂V/∂z0
        real*8 :: x_morse       ! exp(-a*(z-z0))
        integer :: i

        ! --- 转换到表面坐标系 (u,v) ---
        u = x - y / sq3
        v = 2.0d0 * y / sq3

        ! --- 对 Morse 的三个参数 (D, a, z0) 分别做 Fourier 插值 ---
        ! parm(1:3)   = D  at top / bridge / hollow
        ! parm(4:6)   = a  at top / bridge / hollow
        ! parm(7:9)   = z0 at top / bridge / hollow
        do i = 1, 3
            call effective_par_grad(u, v, &
                parm(3*i-2), parm(3*i-1), parm(3*i), &
                p_m(i), dp_dx(i), dp_dy(i))
        end do

        E_0 = parm(10)   ! 能量零点, 与 (u,v) 无关

        ! --- 普通 Morse 势 (无 soft cap) ---
        ! V(z) = D·(x² − 2x) + E₀,  x = exp(−a·(z−z0))
        x_morse = exp(-p_m(2) * (z - p_m(3)))

        Energy = p_m(1) * (x_morse*x_morse - 2.0d0 * x_morse) + E_0

        ! --- dV/dD = x² − 2x ---
        dE_dp(1) = x_morse*x_morse - 2.0d0 * x_morse

        ! --- dV/da = D·(2x−2)·dx/da = D·(2x−2)·(−(z−z0)·x)
        !          = −2·D·x·(x−1)·(z−z0) ---
        dE_dp(2) = -2.0d0 * p_m(1) * x_morse * (x_morse - 1.0d0) * (z - p_m(3))

        ! --- dV/dz0 = D·(2x−2)·dx/dz0 = D·(2x−2)·(a·x)
        !           = 2·D·a·x·(x−1) ---
        dE_dp(3) = 2.0d0 * p_m(1) * p_m(2) * x_morse * (x_morse - 1.0d0)

        ! --- dV/dz = D·(2x−2)·dx/dz = D·(2x−2)·(−a·x)
        !          = −2·D·a·x·(x−1) ---
        dEdz = -2.0d0 * p_m(1) * p_m(2) * x_morse * (x_morse - 1.0d0)

        ! --- (x,y) 梯度 = Σ ∂V/∂p_i · ∂p_i/∂(x,y) ---
        dEdx = dE_dp(1)*dp_dx(1) + dE_dp(2)*dp_dx(2) + dE_dp(3)*dp_dx(3)
        dEdy = dE_dp(1)*dp_dy(1) + dE_dp(2)*dp_dy(2) + dE_dp(3)*dp_dy(3)
    end subroutine

!==============================================================================
!  Calc_coup_and_Grad — 耦合势 V_k(z) 及其梯度
!
!  形式: V(z; A,B,C) = A / [1 + B·exp(C·z)]
!
!  物理意义: 描述吸附粒子与金属表面之间的电子耦合强度随距离的衰减。
!  当 z → +∞ (远表面), exp(C·z) → ∞ (C>0), V → 0  (耦合消失)
!  当 z → −∞ (近表面), exp(C·z) → 0,  V → A  (最大耦合)
!
!  参数 (A, B, C) 随 (u,v) 通过 effective_par_grad 插值:
!    p_m(1) = A(u,v)  — 振幅 [eV]
!    p_m(2) = B(u,v)  — 无量纲前置因子
!    p_m(3) = C(u,v)  — 衰减长度倒数 [1/bohr]
!
!  梯度公式:
!    令 denom = 1 + B·exp(C·z)
!    dV/dA = 1 / denom
!    dV/dB = −A·exp(C·z) / denom²
!    dV/dC = −A·B·z·exp(C·z) / denom²
!    dV/dz = −A·B·C·exp(C·z) / denom²
!==============================================================================
    subroutine Calc_coup_and_Grad(parm, x, y, z, Energy, dEdx, dEdy, dEdz)
        real*8, intent(in)  :: parm(:)  ! [A1,A2,A3, B1,B2,B3, C1,C2,C3, ...]
        real*8, intent(in)  :: x, y, z
        real*8, intent(out) :: Energy, dEdx, dEdy, dEdz

        real*8 :: u, v, denom
        real*8 :: p_m(3)        ! A(u,v), B(u,v), C(u,v) at this (x,y)
        real*8 :: dp_dx(3), dp_dy(3)
        real*8 :: dE_dp(3)      ! ∂V/∂A, ∂V/∂B, ∂V/∂C
        real*8 :: expCz          ! exp(C·z)
        integer :: i

        ! --- 表面坐标系 ---
        u = x - y / sq3
        v = 2.0d0 * y / sq3

        ! --- Fourier 插值 A, B, C ---
        ! parm(1:3)   = A  at top/bridge/hollow
        ! parm(4:6)   = B  at top/bridge/hollow
        ! parm(7:9)   = C  at top/bridge/hollow
        do i = 1, 3
            call effective_par_grad(u, v, &
                parm(3*i-2), parm(3*i-1), parm(3*i), &
                p_m(i), dp_dx(i), dp_dy(i))
        end do

        expCz = exp(p_m(3) * z)
        denom = 1.0d0 + p_m(2) * expCz

        ! --- 势能值 ---
        Energy = p_m(1) / denom

        ! --- ∂V/∂A = 1 / denom ---
        dE_dp(1) = 1.0d0 / denom

        ! --- ∂V/∂B = −A·exp(Cz) / denom² ---
        dE_dp(2) = -p_m(1) * expCz / (denom * denom)

        ! --- ∂V/∂C = −A·B·z·exp(Cz) / denom² ---
        dE_dp(3) = -p_m(1) * p_m(2) * z * expCz / (denom * denom)

        ! --- ∂V/∂z = −A·B·C·exp(Cz) / denom² ---
        dEdz = -p_m(1) * p_m(2) * p_m(3) * expCz / (denom * denom)

        ! --- Cartesian (x,y) 梯度 ---
        dEdx = dE_dp(1)*dp_dx(1) + dE_dp(2)*dp_dx(2) + dE_dp(3)*dp_dx(3)
        dEdy = dE_dp(1)*dp_dy(1) + dE_dp(2)*dp_dy(2) + dE_dp(3)*dp_dy(3)
    end subroutine

!==============================================================================
!  Calc_POT_spec_and_Grad — 扩展 Rydberg 势 (基态 PES)
!
!  形式: V(z; p1..p5) = −p₁ · poly(ρ) · exp(−p₂·ρ) + E₀
!        其中 ρ = |z − p₅|  (原子到表面像平面的距离)
!            poly(ρ) = 1 + p₂·ρ + p₃·ρ² + p₄·ρ³
!
!  物理意义: 描述吸附原子(或分子)与金属表面之间的 physisorption 相互作用。
!  p₅ 是表面像平面 (image plane) 的位置，随 (u,v) 变化。
!  用 abs(z−p₅) 确保 ρ ≥ 0 (势能在表面对称)。
!
!  参数:
!    p_m(1:3) — 标准 Fourier 插值参数 (top/bridge/hollow)
!    p_m(4)   — 由 parm(11:13) 单独插值的参数
!    p_m(5)   — 由 parm(14:16) 单独插值的"像平面"位置 z_im(u,v)
!    E_0      — 能量零点 (parm(10))
!
!  梯度:
!    dV/dρ = −p₁·exp(−p₂·ρ)·[poly'(ρ) − p₂·poly(ρ)]
!    poly'(ρ) = p₂ + 2·p₃·ρ + 3·p₄·ρ²
!    dV/dz = dV/dρ · sign(z − p₅)
!    dV/dp₅ = −dV/dρ · sign(z − p₅) = −dV/dz
!==============================================================================
    subroutine Calc_POT_spec_and_Grad(parm, x, y, z, Energy, dEdx, dEdy, dEdz)
        real*8, intent(in)  :: parm(:)  ! 16 parameters
        real*8, intent(in)  :: x, y, z
        real*8, intent(out) :: Energy, dEdx, dEdy, dEdz

        real*8 :: u, v
        real*8 :: p_m(5)       ! 5 effective parameters at (u,v)
        real*8 :: dp_dx(5), dp_dy(5)
        real*8 :: dE_dp(5)     ! ∂V/∂p_i

        real*8 :: E_0          ! 能量零点
        real*8 :: rho          ! |z − z_im|
        real*8 :: poly         ! 1 + p2·ρ + p3·ρ² + p4·ρ³
        real*8 :: exp_term     ! exp(−p2·ρ)
        real*8 :: dE_drho      ! ∂V/∂ρ
        real*8 :: dpoly_drho   ! poly'(ρ)
        integer :: i

        ! --- 表面坐标系 ---
        u = x - y / sq3
        v = 2.0d0 * y / sq3

        E_0 = parm(10)

        ! --- 前三组参数: parm(1:3), parm(4:6), parm(7:9) ---
        ! 标准 top/bridge/hollow 插值
        do i = 1, 3
            call effective_par_grad(u, v, &
                parm(3*i-2), parm(3*i-1), parm(3*i), &
                p_m(i), dp_dx(i), dp_dy(i))
        end do

        ! --- 第四组参数: parm(11:13)  top/bridge/hollow ---
        call effective_par_grad(u, v, &
            parm(11), parm(12), parm(13), &
            p_m(4), dp_dx(4), dp_dy(4))

        ! --- 第五组参数: parm(14:16)  "像平面"位置 top/bridge/hollow ---
        call effective_par_grad(u, v, &
            parm(14), parm(15), parm(16), &
            p_m(5), dp_dx(5), dp_dy(5))

        ! =====================================================================
        ! 扩展 Rydberg 势:
        !   V(ρ) = −p₁ · (1 + p₂ρ + p₃ρ² + p₄ρ³) · exp(−p₂·ρ) + E₀
        !   ρ = |z − p₅|
        ! =====================================================================

        rho = abs(z) - p_m(5)

        ! 多项式部分: 1 + p₂·ρ + p₃·ρ² + p₄·ρ³
        poly = 1.0d0 + p_m(2)*rho + p_m(3)*rho**2 + p_m(4)*rho**3

        ! 指数衰减因子
        exp_term = exp(-p_m(2) * rho)

        ! 势能值
        Energy = -p_m(1) * poly * exp_term + E_0

        ! --- 解析梯度 ---
        ! poly'(ρ) = p₂ + 2·p₃·ρ + 3·p₄·ρ²
        dpoly_drho = p_m(2) + 2.0d0*p_m(3)*rho + 3.0d0*p_m(4)*rho**2

        ! dV/dρ = −p₁·exp(−p₂ρ)·[poly'(ρ) − p₂·poly(ρ)]
        dE_drho = -p_m(1) * exp_term * (dpoly_drho - p_m(2)*poly)

        ! --- dV/dp_i (i=1..5) ---
        ! ∂V/∂p₁ = −poly · exp(−p₂ρ)
        dE_dp(1) = -poly * exp_term

        ! ∂V/∂p₂ = −p₁·ρ·exp(−p₂ρ)·(1 − poly)
        ! (包含 poly 中的线性项和 exp 对 p₂ 的导数)
        dE_dp(2) = -p_m(1) * rho * exp_term * (1.0d0 - poly)

        ! ∂V/∂p₃ = −p₁·ρ²·exp(−p₂ρ)
        dE_dp(3) = -p_m(1) * rho**2 * exp_term

        ! ∂V/∂p₄ = −p₁·ρ³·exp(−p₂ρ)
        dE_dp(4) = -p_m(1) * rho**3 * exp_term

        ! DEVIATION FIX (VENUS2026 systems/cau_peri, 2026-09-21): the energy
        ! uses rho = |z| - p5, so d(rho)/dp5 = -1 with NO sign factor and
        ! d(rho)/dz = sign(z) - the original sign(z - p5) factors made the
        ! gradient inconsistent with the energy (FD-probed; see
        ! check_cau_peri criterion (c))
        dE_dp(5) = -dE_drho

        ! ∂V/∂z = dV/dρ · dρ/dz = dE_drho · sign(z)   (rho = |z| - p5)
        dEdz = dE_drho * sign(1.0d0, z)

        ! --- Cartesian (x,y) 梯度 = Σ ∂V/∂p_i · ∂p_i/∂(x,y) ---
        dEdx = dot_product(dE_dp, dp_dx)
        dEdy = dot_product(dE_dp, dp_dy)
    end subroutine

!==============================================================================
!  READAP — 从 Analytic_Potential.txt 读取势能参数
!
!  文件格式:
!    Line 3:  C+    D0 D1 D2  a1 a2 a3  z1 z2 z3  e0
!    Line 4:  1D    D0 D1 D2  a1 a2 a3  z1 z2 z3  e0  (¹D 态 Morse)
!    Line 5:  3P    D0 D1 D2  a1 a2 a3  z1 z2 z3  e0  (³P 态 Morse)
!    Line 6:  AN    D0 D1 D2  a1 a2 a3  z1 z2 z3  e0  (负离子 −C 态 Morse)
!    Line 7:  NA    p1..p10                           (基态 Rydberg 第1行)
!    Line 8:  NA2   p11..p16                          (基态 Rydberg 第2行)
!    Line 10: CP    A1 A2 A3  B1 B2 B3  C1 C2 C3     (耦合参数)
!==============================================================================
    subroutine READAP
        character(len=2) :: dum
        open(181, file = 'Analytic_Potential.txt')
        read(181, *)           ! 跳过标题行
        read(181, *)           ! 跳过列标题行
        read(181, *) dum, parm_pC    ! +C 态 (正离子态) Morse 参数
        read(181, *) dum, parm_1D    ! ¹D 态 Morse 参数
        read(181, *) dum, parm_3P    ! ³P 态 Morse 参数
        read(181, *) dum, parm_mC    ! −C 态 (负离子态) Morse 参数
        read(181, *) dum, parm_nas(1:10)   ! 基态 Rydberg: 前10参数
        read(181, *) dum, parm_nas(11:16)  ! 基态 Rydberg: 后6参数
        read(181, *)           ! 跳过标题行
        read(181, *) dum, parm_coup   ! 耦合 V_k 参数
        close(181)
        ! 位点列换位（用户裁决 2026-09-29）：拟合管线输出的 top/hollow 两列是
        ! 互换的（NA 行声明 D_top=6.14 > D_bri=5.42 > D_hol=4.30，物理层级应为
        ! hollow 最深 / top 最浅）。装载后把每个位点三元组的 1<->3 对调，使
        ! 几何 top 点收到携带 top 物理的那一列
        call swap_site_cols(parm_pC)
        call swap_site_cols(parm_1D)
        call swap_site_cols(parm_3P)
        call swap_site_cols(parm_mC)
        call swap_site_cols(parm_coup)
        call swap_site_cols_nas(parm_nas)
    end subroutine

    subroutine swap_site_cols(p)        ! 10 槽 Morse/耦合族: 组 (1:3),(4:6),(7:9)
        real*8, intent(inout) :: p(:)
        integer :: b
        real*8 :: t
        do b = 1, 9, 3
            t = p(b); p(b) = p(b+2); p(b+2) = t
        end do
    end subroutine

    subroutine swap_site_cols_nas(p)    ! 16 槽 Rydberg: 组 1,4,7,11,14（槽 10 为 E_0）
        real*8, intent(inout) :: p(:)
        integer :: b
        real*8 :: t
        do b = 1, 9, 3
            t = p(b); p(b) = p(b+2); p(b+2) = t
        end do
        do b = 11, 14, 3
            t = p(b); p(b) = p(b+2); p(b+2) = t
        end do
    end subroutine

!==============================================================================
!  fMORSE — 标准 Morse 势函数 (已弃用, 保留供参考)
!  dMORSE — Morse 势对 r 的导数 (已弃用, 保留供参考)
!
!  V(r) = D·[exp(−2a(r−re)) − 2·exp(−a(r−re))]
!  dV/dr = D·[−2a·exp(−2a(r−re)) + 2a·exp(−a(r−re))]
!==============================================================================
    real*8 function fMORSE(D, a, re, r)
        real*8, intent(in) :: D, a, re, r
        fMORSE = D * (exp(-2.0d0*a*(r-re)) - 2.0d0*exp(-a*(r-re)))
    end function

    real*8 function dMORSE(D, a, re, r)
        real*8, intent(in) :: D, a, re, r
        dMORSE = D * (-2.0d0*a*exp(-2.0d0*a*(r-re)) + 2.0d0*a*exp(-a*(r-re)))
    end function

!==============================================================================
!  eval_1d_nn — 单隐藏层 NN: z → V(z), dV/dz
!
!  权重布局: w1(1:NH), b1(NH+1:2*NH), w2(2*NH+1:3*NH), b2(3*NH+1)
!  V  = b2 + Σ w2_i * tanh(w1_i * z/z_scale + b1_i)
!  dV = Σ w2_i * (1 - h_i^2) * w1_i / z_scale
!==============================================================================
    subroutine eval_1d_nn(z, weights, nhidden, V, dVdz)
        real*8, intent(in)  :: z
        integer, intent(in) :: nhidden
        real*8, intent(in)  :: weights(:)   ! size = 3*nhidden + 1
        real*8, intent(out) :: V, dVdz
        real*8, parameter   :: z_scale = 3.0d0
        integer :: i, nw
        real*8 :: w1, b1, w2, h, zs

        nw = 3 * nhidden + 1
        zs = z / z_scale
        V = weights(nw)   ! b2
        dVdz = 0.0d0

        do i = 1, nhidden
            w1 = weights(i)
            b1 = weights(nhidden + i)
            w2 = weights(2 * nhidden + i)
            h  = tanh(w1 * zs + b1)
            V   = V   + w2 * h
            dVdz = dVdz + w2 * (1.0d0 - h * h) * w1 / z_scale
        end do
    end subroutine eval_1d_nn

!==============================================================================
!  Calc_NN_and_Grad — 通用 NN 势能评估 (任意势能类型) + Cartesian 梯度
!
!  对给定势能类型 ipot (NN_E0..NN_COUP), 从 nn_weights_all 提取对应权重,
!  计算 3 个位点的 NN(z) 值及其梯度, 再通过 effective_par_grad 插值。
!
!  z-安全机制 (2026-05-26):
!  当 z < z_safe_center (0.35 bohr) 时, 通过 tanh 开关函数平滑过渡到
!  安全势能: 对 PES 态使用排斥墙, 对耦合态趋于零。
!  这防止了 NN 在训练域外 (z < 0) 产生非物理深阱 (~ -68 eV)。
!==============================================================================
    subroutine Calc_NN_and_Grad(ipot, x, y, z, Energy, dEdx, dEdy, dEdz)
        integer, intent(in)  :: ipot
        real*8, intent(in)  :: x, y, z
        real*8, intent(out) :: Energy, dEdx, dEdy, dEdz
        real*8 :: u, v
        real*8 :: V_top, V_bri, V_hol
        real*8 :: dVdz_top, dVdz_bri, dVdz_hol
        real*8 :: dummy_dx, dummy_dy
        real*8 :: V_surf, dVdz_surf, V_inner, dVdz_inner
        integer :: nh, nw, base, nw_surf, nw_inner

        if (.not. allocated(nn_weights_all)) then
            Energy = 0.0d0; dEdx = 0.0d0; dEdy = 0.0d0; dEdz = 0.0d0
            return
        end if

        if (use_rst_nn) then
            nw_surf = 3 * nn_nhidden_surf(ipot) + 1
            nw_inner = 3 * nn_nhidden_inner(ipot) + 1
            base = nn_offset(ipot)

            call eval_1d_nn(z, nn_weights_all(base:base+nw_surf-1), &
                nn_nhidden_surf(ipot), V_surf, dVdz_surf)
            call eval_1d_nn(z, nn_weights_all(base+nw_surf:base+nw_surf+nw_inner-1), &
                nn_nhidden_inner(ipot), V_inner, dVdz_inner)

            Energy = V_surf + V_inner + nn_e0_offset(ipot)
            dEdz = dVdz_surf + dVdz_inner
            dEdx = 0.0d0
            dEdy = 0.0d0
            return
        end if

        nh   = nn_nhidden(ipot)
        nw   = nn_nw(ipot)
        base = nn_offset(ipot)

        u = x - y / sq3
        v = 2.0d0 * y / sq3

        ! 位点块换位（用户裁决 2026-09-29）：权重文件的每类型位点块相对物理
        ! 是 (hollow, bridge, top) 序——第三块读作 top、第一块读作 hollow
        call eval_1d_nn(z, nn_weights_all(base+2*nw : base+3*nw-1), nh, V_top, dVdz_top)
        call eval_1d_nn(z, nn_weights_all(base+nw : base+2*nw-1), nh, V_bri, dVdz_bri)
        call eval_1d_nn(z, nn_weights_all(base : base+nw-1), nh, V_hol, dVdz_hol)

        call effective_par_grad(u, v, V_top, V_bri, V_hol, Energy, dEdx, dEdy)
        Energy = Energy + nn_e0_offset(ipot)
        ! dEdx, dEdy remain in eV/A — caller converts to kcal/mol via ×23.0605

        call effective_par_grad(u, v, dVdz_top, dVdz_bri, dVdz_hol, dEdz, dummy_dx, dummy_dy)
        ! dEdz remains in eV/A — caller converts to kcal/mol via ×23.0605
    end subroutine Calc_NN_and_Grad

!==============================================================================
!  Calc_NN_POT_and_Grad — NN 基态势能 + Cartesian 梯度 (兼容包装器)
!==============================================================================
    subroutine Calc_NN_POT_and_Grad(x, y, z, Energy, dEdx, dEdy, dEdz)
        real*8, intent(in)  :: x, y, z
        real*8, intent(out) :: Energy, dEdx, dEdy, dEdz
        call Calc_NN_and_Grad(NN_E0, x, y, z, Energy, dEdx, dEdy, dEdz)
    end subroutine Calc_NN_POT_and_Grad

!==============================================================================
!  load_all_nn_weights — 从主文件加载全部 6 类势能的 NN 权重
!
!  文件格式:
!    # 注释行
!    # SECTION: <name>  NH=<int>   ← 识别势能类型 + 隐藏神经元数
!    ## Top / ## Bridge / ## Hollow  ← 位点标签
!    每段 (3*NH+1) 个实数
!
!  两遍扫描: 第一遍收集 NH → 分配 → 第二遍读权重
!==============================================================================
    subroutine load_all_nn_weights(filename, ierror)
        character(len=*), intent(in)  :: filename
        integer, intent(out) :: ierror
        character(len=200) :: line, section_name
        integer :: i, j, unit, ios, ipot, nh, nw, total_all, idx
        integer :: sections_found, nh_pos, pending_ipot, sec_pos, state_pos
        integer :: nh_surf, nh_inner, nh_pos2, col1, col2
        integer :: nw_surf, nw_inner, base
        logical :: file_exists, pending_nh

        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            write(*,*) 'load_all_nn_weights: file not found: ', trim(filename)
            write(*,*) '  Falling back to analytic potentials for all types.'
            use_nn_pes = .false.
            ierror = 1
            return
        end if

        open(newunit=unit, file=trim(filename), status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'load_all_nn_weights: cannot open ', trim(filename)
            use_nn_pes = .false.
            ierror = 1
            return
        end if

        ! ---- Pass 1: scan SECTION/STATE headers for NH values ----
        nn_nhidden = 0
        sections_found = 0
        pending_nh = .false.
        pending_ipot = 0

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) /= '#') cycle

            ! ---- Handle pending NH= from previous STATE: line ----
            if (pending_nh) then
                nh_pos = index(line, 'NH:')
                if (nh_pos > 0) then
                    read(line(nh_pos+3:), *, iostat=ios) nh
                    if (ios /= 0 .or. nh < 1 .or. nh > NN_MAX_HIDDEN) then
                        write(*,*) 'load_all_nn_weights: bad NH in # NH: line'
                        use_nn_pes = .false.
                        close(unit)
                        ierror = 1
                        return
                    end if
                    ipot = pending_ipot
                    nn_nhidden(ipot) = nh
                    nn_nw(ipot)     = 3 * nh + 1
                    nn_total(ipot)  = 3 * nn_nw(ipot)
                    sections_found = sections_found + 1
                    pending_nh = .false.
                    pending_ipot = 0
                else
                    ! Check for RST format: # NH_surf: <n>  NH_inner: <m>
                    nh_pos = index(line, 'NH_surf:')
                    if (nh_pos > 0) then
                        nh_pos2 = index(line, 'NH_inner:')
                        if (nh_pos2 == 0) then
                            write(*,*) 'load_all_nn_weights: expected NH_inner: in RST line: ', trim(line)
                            use_nn_pes = .false.
                            close(unit)
                            ierror = 1
                            return
                        end if
                        col1 = index(line(nh_pos:), ':') + nh_pos - 1
                        col2 = index(line(nh_pos2:), ':') + nh_pos2 - 1
                        read(line(col1+1:nh_pos2-1), *, iostat=ios) nh_surf
                        if (ios == 0) read(line(col2+1:), *, iostat=ios) nh_inner
                        if (ios /= 0 .or. nh_surf < 1 .or. nh_surf > NN_MAX_HIDDEN &
                            .or. nh_inner < 1 .or. nh_inner > NN_MAX_HIDDEN) then
                            write(*,*) 'load_all_nn_weights: bad NH_surf/NH_inner in line: ', trim(line)
                            use_nn_pes = .false.
                            close(unit)
                            ierror = 1
                            return
                        end if
                        use_rst_nn = .true.
                        ipot = pending_ipot
                        nn_nhidden_surf(ipot) = nh_surf
                        nn_nhidden_inner(ipot) = nh_inner
                        nn_nhidden(ipot) = nh_surf  ! sentinel for section check
                        nn_nw(ipot) = (3*nh_surf + 1) + (3*nh_inner + 1)
                        nn_total(ipot) = nn_nw(ipot)
                        sections_found = sections_found + 1
                        pending_nh = .false.
                        pending_ipot = 0
                    else
                        write(*,*) 'load_all_nn_weights: expected # NH: or # NH_surf: after STATE, got: ', trim(line)
                        use_nn_pes = .false.
                        close(unit)
                        ierror = 1
                        return
                    end if
                end if
                cycle
            end if

            ! ---- Detect section/state header ----
            sec_pos = index(line, 'SECTION:')
            state_pos = index(line, 'STATE:')
            ! Filter out embedded references: # text ... # STATE: ...
            if (sec_pos > 0 .and. len_trim(line(2:sec_pos-1)) > 0) sec_pos = 0
            if (state_pos > 0 .and. len_trim(line(2:state_pos-1)) > 0) state_pos = 0
            if (sec_pos == 0 .and. state_pos == 0) cycle

            if (sec_pos > 0) then
                section_name = adjustl(line(sec_pos+8:))
            else
                section_name = adjustl(line(state_pos+6:))
            end if

            ! Extract NH= value if present inline
            nh_pos = index(section_name, 'NH=')
            if (nh_pos > 0) then
                read(section_name(nh_pos+3:), *, iostat=ios) nh
                if (ios /= 0 .or. nh < 1 .or. nh > NN_MAX_HIDDEN) then
                    write(*,*) 'load_all_nn_weights: bad NH in section: ', trim(section_name)
                    use_nn_pes = .false.
                    close(unit)
                    ierror = 1
                    return
                end if
                section_name = adjustl(section_name(1:nh_pos-1))
            end if

            ! Map section name to potential type index
            select case (trim(section_name))
            case ('E0', 'GS');    ipot = NN_E0
            case ('1D');          ipot = NN_1D
            case ('3P');          ipot = NN_3P
            case ('mC', 'AN');    ipot = NN_mC
            case ('pC', 'C+');    ipot = NN_pC
            case ('COUP', 'Vcp'); ipot = NN_COUP
            case ('E_h', 'EH');  ipot = NN_Eh
            case default
                write(*,*) 'load_all_nn_weights: WARNING: unknown section: ', trim(section_name)
                ipot = 0
            end select

            if (ipot > 0) then
                if (nh_pos > 0) then
                    ! Old format: NH= was inline, store immediately
                    nn_nhidden(ipot) = nh
                    nn_nw(ipot)     = 3 * nh + 1
                    nn_total(ipot)  = 3 * nn_nw(ipot)
                    sections_found = sections_found + 1
                else
                    ! New format: NH on next # NH: line
                    pending_nh = .true.
                    pending_ipot = ipot
                end if
            end if
        end do

        ! All-or-nothing: all 6 required sections must be present (NN_Eh is optional)
        if (use_rst_nn) then
            if (any(nn_nhidden_surf(1:NN_Eh-1) == 0)) then
                write(*,*) 'load_all_nn_weights: missing required sections (RST)'
                write(*,*) '  Falling back to analytic potentials for all types.'
                use_nn_pes = .false.
                close(unit)
                ierror = 1
                return
            end if
        else
            if (any(nn_nhidden(1:NN_Eh-1) == 0)) then
                write(*,*) 'load_all_nn_weights: missing required sections'
                write(*,*) '  Falling back to analytic potentials for all types.'
                use_nn_pes = .false.
                close(unit)
                ierror = 1
                return
            end if
        end if
        ! E_h is optional — zero out if not found
        if (use_rst_nn) then
            if (nn_nhidden_surf(NN_Eh) == 0) then
                nn_nw(NN_Eh) = 0
                nn_total(NN_Eh) = 0
            end if
        else
            if (nn_nhidden(NN_Eh) == 0) then
                nn_nw(NN_Eh) = 0
                nn_total(NN_Eh) = 0
            end if
        end if

        ! Compute offsets and total size
        nn_offset(1) = 1
        do ipot = 2, NN_POT_TYPES
            nn_offset(ipot) = nn_offset(ipot-1) + nn_total(ipot-1)
        end do
        total_all = nn_offset(NN_POT_TYPES) + nn_total(NN_POT_TYPES) - 1

        ! Allocate and initialise
        if (allocated(nn_weights_all)) deallocate(nn_weights_all)
        allocate(nn_weights_all(total_all))
        nn_weights_all = 0.0d0

        ! ---- Pass 2: rewind and read actual weights ----
        rewind(unit)
        ipot = 0

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            line = adjustl(line)
            if (len_trim(line) == 0) cycle

            ! Detect section header (SECTION: or STATE:)
            if (line(1:1) == '#') then
                sec_pos = index(line, 'SECTION:')
                state_pos = index(line, 'STATE:')
                ! Filter out embedded references: # text ... # STATE: ...
                if (sec_pos > 0 .and. len_trim(line(2:sec_pos-1)) > 0) sec_pos = 0
                if (state_pos > 0 .and. len_trim(line(2:state_pos-1)) > 0) state_pos = 0
                if (sec_pos > 0 .or. state_pos > 0) then
                    if (sec_pos > 0) then
                        section_name = adjustl(line(sec_pos+8:))
                    else
                        section_name = adjustl(line(state_pos+6:))
                    end if
                    nh_pos = index(section_name, 'NH=')
                    if (nh_pos > 0) section_name = adjustl(section_name(1:nh_pos-1))
                    select case (trim(section_name))
                    case ('E0', 'GS');    ipot = NN_E0
                    case ('1D');          ipot = NN_1D
                    case ('3P');          ipot = NN_3P
                    case ('mC', 'AN');    ipot = NN_mC
                    case ('pC', 'C+');    ipot = NN_pC
                    case ('COUP', 'Vcp'); ipot = NN_COUP
                    case ('E_h', 'EH');  ipot = NN_Eh
                    case default;         ipot = 0
                    end select
                    cycle
                end if
            end if

            if (line(1:1) == '#') cycle
            if (line(1:2) == '##') cycle
            if (ipot == 0) cycle

            if (use_rst_nn) then
                ! RST: read Surface_Layer weights, skip ## Inner_Layer, read Inner_Layer weights
                nw_surf = 3 * nn_nhidden_surf(ipot) + 1
                nw_inner = 3 * nn_nhidden_inner(ipot) + 1
                base = nn_offset(ipot)

                ! Read Surface_Layer weights
                backspace(unit)
                do i = 1, nw_surf
                    read(unit, *, iostat=ios) nn_weights_all(base + i - 1)
                    if (ios /= 0) then
                        write(*,*) 'load_all_nn_weights: ERROR reading RST surf weights, ipot=', ipot
                        use_nn_pes = .false.
                        close(unit)
                        ierror = 1
                        return
                    end if
                end do

                ! Skip ## Inner_Layer marker
                read(unit, '(A)', iostat=ios) line

                ! Read Inner_Layer weights
                do i = 1, nw_inner
                    read(unit, *, iostat=ios) nn_weights_all(base + nw_surf + i - 1)
                    if (ios /= 0) then
                        write(*,*) 'load_all_nn_weights: ERROR reading RST inner weights, ipot=', ipot
                        use_nn_pes = .false.
                        close(unit)
                        ierror = 1
                        return
                    end if
                end do
            else
                nw = nn_nw(ipot)
                backspace(unit)
                do j = 1, 3
                    do i = 1, nw
                        idx = nn_offset(ipot) + (j-1)*nw + i - 1
                        read(unit, *, iostat=ios) nn_weights_all(idx)
                        if (ios /= 0) then
                            write(*,*) 'load_all_nn_weights: ERROR reading weights, ipot=', ipot, ' site=', j
                            use_nn_pes = .false.
                            close(unit)
                            ierror = 1
                            return
                        end if
                    end do
                    if (j < 3) read(unit, '(A)', iostat=ios) line  ! skip ## site label
                end do
            end if
            ! Try to read E0 offset after Hollow weights
            read(unit, '(A)', iostat=ios) line
            if (ios == 0) then
                line = adjustl(line)
                if (index(line, '# E0') > 0 .and. &
                    index(line, 'STATE:') == 0 .and. &
                    index(line, 'SECTION:') == 0 .and. &
                    index(line, 'EHF') == 0) then
                    read(unit, *, iostat=ios) nn_e0_offset(ipot)
                    if (ios /= 0) nn_e0_offset(ipot) = 0.0d0
                else
                    backspace(unit)
                end if
            end if
            ipot = 0
        end do

        close(unit)

        use_nn_pes = .true.
        ierror = 0
        use_eh_correction = (nn_nhidden(NN_Eh) > 0)
        if (use_rst_nn) then
            write(*,*) 'RST NN weights loaded: all ', NN_POT_TYPES, &
                ' types from ', trim(filename)
            write(*,*) '  Total weights allocated: ', total_all
            do ipot = 1, NN_POT_TYPES
                write(*,*) '  Type', ipot, ': NH_surf=', nn_nhidden_surf(ipot), &
                    '  NH_inner=', nn_nhidden_inner(ipot), &
                    '  NW/total=', nn_total(ipot)
            end do
        else
            write(*,*) 'NN weights loaded: all ', NN_POT_TYPES, &
                ' types from ', trim(filename)
            write(*,*) '  Total weights allocated: ', total_all
            if (use_eh_correction) then
                write(*,*) '  E_h correction enabled (NN type 7)'
            end if
            do ipot = 1, NN_POT_TYPES
                write(*,*) '  Type', ipot, ': NH=', nn_nhidden(ipot), &
                    '  NW/site=', nn_nw(ipot), '  total=', nn_total(ipot)
            end do
        end if
    end subroutine load_all_nn_weights

end module
