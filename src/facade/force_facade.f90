!=====================================================================
! force_facade.f90 — L2 接口 A 门面：(Q,P,s)→(F,E)。内部 Σ甲类力成员（可叠加）
!                    × 单选乙类电子律成员；单位契约履行点（×23.0605×e_conv 仅此处）
! 参考: venus2026/DVDQ.f:10-22（PDOT 清零 + DPESHON 分派——加性合成的旧形态）；
!       src_TEST/interface_TEST.f90:12-16（POTPRE 成员初始化钩子）、:18-30
!       （POT0 能量换算 ×23.0605）、:32-44（DPESHON 梯度累加
!       PDOT(i)=PDOT(i)−g(i)*23.0605*C1——新门面换算与之逐字符同形，N1.1
!       bit 级可比）；ENERGY.f:54-68（T=Σp²/2m/C1、V=POT0 口径、H=T+V）；
!       docs/final_completion_plan.md §5 N1（注册/装配/能力级协商表规格——直接生效）
! 成员单位契约：甲类成员交付 eV 单位力核/能量（eV 为交换单位——沿旧
! interface_TEST 的 v/g 链）；eV→内部单位换算集中本门面（T4 grep 要求 23.0605
! 只出现于本文件——N1 成员签名注的"F_add code 单位"据此读作"经门面换算后"）
! 成员表 C2 填充（pes_* 六文件延后——K2 计划 E3；主势+加性随机律成员）
!
! —— E9 多面升级（用户裁决 2026-09-07：立即升级 + 双签名并存 + 一包全升）——
! 甲类新增两条多面签名（与单面标量签名并存，成员按数据源天然形态自选）：
!   pes_diab_i  透热表象：交付 H 矩阵 + 矩阵元梯度（IESH 模型哈密顿量天然此形；
!               绝热视图经 utils/linalg 对角化现算）
!   pes_adiab_i 绝热表象：交付各面 V_i/力 + 当前面导数耦合 d_i,·（仅当前面行
!               ——FSSH 跃迁只消费当前面耦合，O(ndof·n_surf) 而非全对）
! 乙类签名扩 dt（跃迁概率 g_ij(dt) 依赖步长；radau 序列边界"一步"无良定义
! ——协商表以积分器轴排除该组合，见 triad_tbl）。能力级协商数据化：合法组合
! 从写死单行改为 triad_tbl 数据表（加行=加数据不改代码）。
!=====================================================================
module force_facade
   use consts,        only: e_conv
   use state,         only: state_t, elec_t
   use config,        only: electronic, propagator
   use config_roster, only: roster
   implicit none
   private
   public :: facade_reg_pes, facade_reg_pes_diab, facade_reg_pes_adiab
   public :: facade_reg_elec, facade_init, facade_force, facade_info

   ! eV→kcal/mol 换算（单位契约字面值——T4：全库仅本文件出现）
   real(8), parameter :: ev_kcal = 23.0605d0

   ! 甲类成员（可叠加）接口：Q,W → 力核+能量贡献（eV 交换单位；换算在门面边界）
   abstract interface
      subroutine pes_proc_i(q, mass, g_add, v_add)
         real(8), intent(in)  :: q(:)     ! 坐标 [Å]（全自由度展平）
         real(8), intent(in)  :: mass(:)  ! 每原子质量 [amu]（随机律/质量律成员用；
                                          ! 纯势成员可忽略）
         real(8), intent(out) :: g_add(:) ! 力核 −dV/dq（随机律为等效力核）[eV/Å]
         real(8), intent(out) :: v_add    ! 能量贡献 [eV]
      end subroutine pes_proc_i
      ! 注：随机律成员内部调 rng（N1——签名承认随机律；单流由 rng 担保）
   end interface

   ! 甲类多面成员接口·透热表象（E9 双签名之一）：Q,W → H 矩阵 + 矩阵元梯度
   abstract interface
      subroutine pes_diab_i(q, mass, h_ij, dh_ij)
         real(8), intent(in)  :: q(:)      ! 坐标 [Å]（全自由度展平）
         real(8), intent(in)  :: mass(:)   ! 每原子质量 [amu]
         real(8), intent(out) :: h_ij(:,:)   ! 透热哈密顿矩阵 H_ij(Q) [eV]
                                            ! （n_surf×n_surf 对称；成员交付全矩阵）
         real(8), intent(out) :: dh_ij(:,:,:) ! 矩阵元梯度 ∂H_ij/∂q [eV/Å]
                                            ! （ndof×n_surf×n_surf；对称性成员自理。
                                            ! 大态空间（IESH 上百态）压缩交付
                                            ! 属文法版本化事项——消费者入场时定形）
      end subroutine pes_diab_i
   end interface

   ! 甲类多面成员接口·绝热表象（E9 双签名之二）：Q,W,当前面 → 各面量 + 当前面耦合
   abstract interface
      subroutine pes_adiab_i(q, mass, i_surf, v_surf, g_surf, d_ij)
         real(8), intent(in)  :: q(:)       ! 坐标 [Å]（全自由度展平）
         real(8), intent(in)  :: mass(:)    ! 每原子质量 [amu]
         integer, intent(in)  :: i_surf     ! 当前面索引 [—]（跃迁后随新面重算——
                                            ! FSSH hop 后力重算位）
         real(8), intent(out) :: v_surf(:)  ! 各面能量 V_i(Q) [eV]（n_surf——
                                            ! 本征值天然全给）
         real(8), intent(out) :: g_surf(:,:) ! 各面力核 −dV_i/dq [eV/Å]
                                            ! （ndof×n_surf）
         real(8), intent(out) :: d_ij(:,:)  ! 当前面导数耦合 d_i,j = <φ_i|∇_q|φ_j>
                                            ! [1/Å]（ndof×n_surf——仅当前面行，
                                            ! O(ndof·n_surf) 而非全对 n_surf²/2；
                                            ! 能隙加权变体（NAC×能隙）交成员/
                                            ! 消费者自理，契约只认裸 d）
      end subroutine pes_adiab_i
   end interface

   ! 乙类成员（单选）接口：电子传播律——电子态演化 + 核力修饰（绝热=两退化点）
   abstract interface
      subroutine elec_proc_i(s, f, dt)
         import :: elec_t
         type(elec_t), intent(inout) :: s  ! 电子态（律成员按律演化 amp/i_surf/
                                           ! n_hop；绝热：不变）
         real(8), intent(inout) :: f(:)    ! 核力 [内部能量/Å]（律成员修饰——面
                                           ! 选择/平均/切换与动能重整；绝热：恒等）
         real(8), intent(in)    :: dt      ! 当前实际步长 [10 fs]（跃迁概率
                                           ! g_ij(dt) 所需——E9；绝热成员忽略）
      end subroutine elec_proc_i
   end interface

   ! 注册表（甲类可叠加表 + 乙类单选槽；装配相位登记，演化期只读）
   integer, parameter :: max_pes = 16     ! 甲类成员表容量 [行]（主势+加性随机律；C2 填充）
   integer, parameter :: k_pes_scalar = 1 ! 甲类成员类别：单面标量（g_add/v_add——现行）
   integer, parameter :: k_pes_diab   = 2 ! 甲类成员类别：多面透热（h_ij/dh_ij——E9）
   integer, parameter :: k_pes_adiab  = 3 ! 甲类成员类别：多面绝热（v_surf/g_surf/d_ij——E9）
   type :: pes_entry_t
      character(len=32) :: name = ''      ! 成员名（config%electronic%pes_main 装配键）
      integer :: kind = k_pes_scalar      ! 成员类别（k_pes_*——按类别调对应槽）
      procedure(pes_proc_i),  pointer, nopass :: proc  => null() ! 标量求值例程
      procedure(pes_diab_i),  pointer, nopass :: diab  => null() ! 透热求值例程（E9）
      procedure(pes_adiab_i), pointer, nopass :: adiab => null() ! 绝热求值例程（E9）
   end type pes_entry_t
   type(pes_entry_t) :: pes_tbl(max_pes)  ! 甲类注册表（模块私有）
   integer :: n_pes = 0                   ! 已登记甲类行数 [行]
   character(len=32) :: elec_name = ''    ! 乙类成员名（模块私有；单选槽；空=未装配）
   procedure(elec_proc_i), pointer :: elec_proc ! 乙类例程（模块私有；装配相位
                                                ! facade_init 挂槽；nopass 仅组件属性）

   ! 能力级合法组合表（E9 数据化——原写死单行改表数据；加行=加数据不改代码）
   ! 四轴：H 能力级 n_surf 区间 × 主势类别 × 电子律 × 积分器（积分器轴以三布尔
   ! 承载）。radau 列在多面行恒假：自适应序列边界"一步"无良定义，跃迁概率
   ! g_ij(dt) 失据——非绝热只配定步长积分器（verlet/symple）。
   integer, parameter :: m_adiabatic = 1  ! 电子律编号：绝热退化点（config%electronic%method 键）
   integer, parameter :: m_tsh      = 2   ! 电子律编号：TSH/FSSH（将来 elec_tsh 注册键）
   integer, parameter :: m_iesh     = 3   ! 电子律编号：IESH（将来 elec_iesh 注册键）
   integer, parameter :: max_triad = 4    ! 合法组合表容量 [行]
   integer, parameter :: n_surf_unbounded = huge(1) ! 多面行能力级上限（无界意）
   type :: triad_row_t
      integer :: n_surf_lo = 1            ! H 能力级下限（n_surf ≥）
      integer :: n_surf_hi = 1            ! H 能力级上限（单面行 1..1；多面行 ..无界）
      integer :: pes_kind  = k_pes_scalar ! 主势成员类别（k_pes_*）
      integer :: method    = m_adiabatic  ! 电子律编号（m_*）
      logical :: ok_verlet = .true.       ! 积分器 1（verlet/beeman）放行位
      logical :: ok_symple = .true.       ! 积分器 2（symple）放行位
      logical :: ok_radau  = .true.       ! 积分器 3（radau）放行位（多面行假）
   end type triad_row_t
   type(triad_row_t), parameter :: triad_tbl(max_triad) = (/ &
      ! 行 1（现行）：单面标量 × 绝热 × 三积分器全放行
      triad_row_t(1, 1, k_pes_scalar, m_adiabatic, .true., .true., .true.), &
      ! 行 2：多面透热 × TSH × 定步长积分器
      triad_row_t(2, n_surf_unbounded, k_pes_diab,  m_tsh, .true., .true., .false.), &
      ! 行 3：多面绝热 × TSH × 定步长积分器
      triad_row_t(2, n_surf_unbounded, k_pes_adiab, m_tsh, .true., .true., .false.), &
      ! 行 4：多面透热 × IESH × 定步长积分器（iesh 透热天然——模型哈密顿量；
      !       adiab×iesh 不放行）
      triad_row_t(2, n_surf_unbounded, k_pes_diab,  m_iesh, .true., .true., .false.) &
      /)
   ! 注：行 2-4 的 method（tsh/iesh）当前无乙类成员注册——协商表放行在先，
   ! facade_init 乙类装配（按 method 查注册槽）loud-stop 兜底拦截，双闸门
contains
   !------------------------------------------------------------------
   ! facade_reg_pes(name, proc) — 登记甲类力成员·单面标量（可叠加：主势 + 随机律/约束势）
   ! 参考: docs/final_completion_plan.md §5 N1（注册规格）；T2（新成员=新文件+一行注册）
   subroutine facade_reg_pes(name, proc)
      character(len=*), intent(in) :: name  ! 成员名（装配键——config%electronic%pes_main）
      procedure(pes_proc_i) :: proc         ! 成员求值例程
      ! 伪代码：
      ! 1. 重名或表满 → loud-stop
      ! 2. 登记行 {name, kind=k_pes_scalar, proc}（装配相位执行——driver；
      !    pes_* 成员 C2 引入后各注册一行）
   end subroutine facade_reg_pes

   !------------------------------------------------------------------
   ! facade_reg_pes_diab(name, proc) — 登记甲类力成员·多面透热（E9 双签名之一：
   !                                主势专用——多面行不参与加性叠加，见 facade_init）
   ! 参考: docs/final_completion_plan.md §5 N1（注册规格）；E9 裁决（双签名并存）
   subroutine facade_reg_pes_diab(name, proc)
      character(len=*), intent(in) :: name  ! 成员名（装配键——config%electronic%pes_main）
      procedure(pes_diab_i) :: proc         ! 成员透热求值例程（H 矩阵 + 梯度）
      ! 伪代码：
      ! 1. 重名或表满 → loud-stop
      ! 2. 登记行 {name, kind=k_pes_diab, diab}（装配相位执行——driver；
      !    多面 pes 成员随 elec_tsh 消费者引入）
   end subroutine facade_reg_pes_diab

   !------------------------------------------------------------------
   ! facade_reg_pes_adiab(name, proc) — 登记甲类力成员·多面绝热（E9 双签名之二：
   !                                各面 V_i/力 + 当前面导数耦合）
   ! 参考: docs/final_completion_plan.md §5 N1（注册规格）；E9 裁决（双签名并存）
   subroutine facade_reg_pes_adiab(name, proc)
      character(len=*), intent(in) :: name  ! 成员名（装配键——config%electronic%pes_main）
      procedure(pes_adiab_i) :: proc        ! 成员绝热求值例程
      ! 伪代码：
      ! 1. 重名或表满 → loud-stop
      ! 2. 登记行 {name, kind=k_pes_adiab, adiab}（装配相位执行——driver）
   end subroutine facade_reg_pes_adiab

   !------------------------------------------------------------------
   ! facade_reg_elec(name, proc) — 登记乙类电子律成员（单选）
   ! 参考: docs/final_completion_plan.md §5 N1（乙类单选规格）
   subroutine facade_reg_elec(name, proc)
      character(len=*), intent(in) :: name  ! 成员名（config%electronic%method 映射键）
      procedure(elec_proc_i) :: proc        ! 成员例程
      ! 伪代码：
      ! 1. 槽已占 → loud-stop（乙类单选——重复注册即装配错）
      ! 2. 登记槽 {name, proc}（elec_adiabatic 退化成员由装配相位注册）
   end subroutine facade_reg_elec

   !------------------------------------------------------------------
   ! facade_init() — 成员装配 + 能力级组合校验（非法组合 loud-stop）
   ! 参考: docs/final_completion_plan.md §5 N1（协商表规格）；E9（表数据化——
   !       triad_tbl 四行：现行单面绝热 + 多面 tsh/iesh 预放行）；
   !       src_TEST/interface_TEST.f90:12-16（POTPRE 型成员初始化钩子——D3
   !       ghost 参数装配位）；docs/design_spec_completeness.md §3.2-2
   !       （能力级与方法的合法性做类型检查，非法组合显式拒绝并指明冲突的两轴）
   subroutine facade_init()
      ! 伪代码：
      ! 1. 甲类装配：按 config%electronic%pes_main 从 pes_tbl 固化主势成员行
      !    + 加性成员行（加性成员表——config%electronic 留白字段，C2 随旧开关
      !    NGLO/NFC 映射填充，文法向后兼容）；装配表空（无主势成员）→ loud-stop。
      !    多面体系（n_surf ≥ 2）：主势行 kind 必须为 k_pes_diab/k_pes_adiab
      !    （标量主势配多面 → loud-stop 指明冲突）；加性成员限标量类（随机律/
      !    约束势与电子态无关——只进核力通道，H 矩阵通道唯主势承载，见
      !    facade_force 步 4'）
      ! 2. 乙类装配：按 config%electronic%method 选律成员（1=绝热 →
      !    elec_adiabatic 退化成员；2/3=tsh/iesh 将来注册）；未登记 → loud-stop
      !    （协商表放行在先、注册槽空此后拦——双闸门，triad_tbl 注）
      ! 3. 成员初始化钩子（POTPRE 型）：逐装配成员调其 init（参数经成员自带
      !    params——D3；钩子槽的接线形态 C2 定形）
      ! 4. 能力级组合校验（E9 查表化）：四元组
      !    (n_surf ∈ 行区间) × (主势 kind = 行 pes_kind) × (method = 行 method)
      !    × (积分器 id 放行位——config%propagator%id 映射 ok_verlet/ok_symple/
      !    ok_radau) 逐行匹配 triad_tbl；无命中行 → 打印冲突各轴现值并 STOP
      !    （规范 §3.2-2 类型检查——指明两轴；典型拦截：多面 × radau——
      !    跃迁概率 g_ij(dt) 在自适应序列边界失据）
   end subroutine facade_init

   !------------------------------------------------------------------
   ! facade_force(sta, e_tot[, dt]) — (Q,P,s)→(F,E)：Σ成员力核 → 换算 → 电子律修饰；
   !                            E = T(P) + ΣV；多面路径 E9（H 矩阵通道唯主势）
   ! 参考: venus2026/DVDQ.f:10-22（清零+分派）；src_TEST/interface_TEST.f90:32-44
   !       （同形累加）、:18-30（能量口径）；ENERGY.f:54-68（T/H 报告口径）
   subroutine facade_force(sta, e_tot, dt)
      type(state_t), intent(inout) :: sta  ! 物理态（q/p/s 读入；f 写出——Σ合成力）
      real(8), intent(out) :: e_tot        ! 总能 H = T(P)+ΣV [kcal/mol]（N1 口径；
                                           ! 多面路径 V=当前面/律成员口径，步 5'）
      real(8), intent(in), optional :: dt  ! 当前实际步长 [10 fs]（跃迁概率
                                           ! g_ij(dt) 通道——E9；缺省取
                                           ! config%propagator%dt——定步长语义；
                                           ! 绝热律忽略之）
      ! 伪代码：
      ! 1. 力核工作区 g(:) 清零、v_sum = 0（旧 DVDQ.f:12-14 PDOT 清零位）；
      !    dt_eff = 缺省 config%propagator%dt（prop_cfg 别名）或传入值
      ! 2. 单面路径（n_surf = 1——现行）：逐装配甲类成员：
      !    proc(sta%q, roster%mass, g_add, v_add)；g(:) += g_add、v_sum += v_add
      !    （eV 累加；随机律成员内部调 rng）
      ! 3. 单位契约履行点（T4——全库唯一）：
      !    sta%f 预清零后逐元素 sta%f(i) = sta%f(i) − g(i)·ev_kcal·e_conv
      !    （与旧 DPESHON "PDOT(i)=PDOT(i)-g(i)*23.0605*C1" 逐字符同形——N1.1
      !    bit 级可比；e_conv 即旧 C1）
      ! 4. 电子律修饰：elec_proc(sta%s, sta%f, dt_eff)（绝热成员恒等——单面表示下
      !    复合塌缩为纯 PES，规范 §3.2-2；E9 后签名带 dt）
      ! 5. 能量：v_pot = v_sum·ev_kcal [kcal/mol]（旧 POT0 :26 口径）；
      !    t_kin = Σ p²/(2·roster%mass)/e_conv [kcal/mol]（逐原子折叠——
      !    旧 ENERGY.f:58-66 口径）；e_tot = t_kin + v_pot（旧 :68 H=T+V）
      ! ————————————————————————————————————————————————
      ! 2'. 多面路径（n_surf ≥ 2，E9——协商表放行后）：
      !    a. 主势按表象分派：kind=diab → diab(sta%q, roster%mass, h_ij, dh_ij)
      !       （透热矩阵与梯度，eV 通道；绝热视图消费方经 utils/linalg 对角化
      !       现算）；kind=adiab → adiab(sta%q, roster%mass, sta%s%i_surf,
      !       v_surf, g_surf, d_ij)（当前面行耦合——hop 后随新面重算）
      !    b. 加性标量成员照常叠加 g(:)（随机律/约束势与电子态无关——只进
      !       核力通道；H 矩阵通道 h_ij/dh_ij 唯主势承载）
      !    c. 换算同步 3（eV→内部单位，契约唯此点）
      !    d. elec_proc(sta%s, sta%f, dt_eff)：非绝热律在此完成 amp 传播、
      !       跃迁判定（rng 单流）与 i_surf 切换后的力修饰/动能重整——
      !       消费细节随 elec_tsh/elec_iesh 成员定形（C 阶段）
      !    e. 能量口径：t_kin 同步 5；v_pot = 当前面 V（adiab：v_surf(i_surf)；
      !       diab：布居加权或当前面——律成员口径，消费者定形）；
      !       e_tot = t_kin + v_pot·ev_kcal
   end subroutine facade_force

   !------------------------------------------------------------------
   ! facade_info() — 诊断输出：成员清单与能力级组合
   ! 参考: docs/final_completion_plan.md §5 N1（facade_info——成员清单/能力级）
   subroutine facade_info()
      ! 伪代码：
      ! 1. 打印甲类注册/装配成员清单（名 + 主势/加性标注 + 类别
      !    scalar/diab/adiab）与乙类成员名
      ! 2. 打印当前能力级组合：H 表示 n_surf × 主势类别 × 电子律 × 积分器
      !    （config 三源）+ triad_tbl 命中行号（未装配则示空）
      ! 3. 打印单位契约常量（ev_kcal/e_conv）与装配状态——排查装配错用
   end subroutine facade_info
end module force_facade
