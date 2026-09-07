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
! —— E10 工作区修订（用户裁决 2026-09-07，三公理立法——修正 E9）——
!   A 接口不变性：abstract interface 对体系不变——签名只传容器 ws_t（facade/
!     elec_ws.f90），不传 h_ij/d_ij 等形态数组；表象（透热/绝热）从签名分叉
!     降级为 ws%repr 数据标记。E9 的 pes_diab_i/pes_adiab_i 双签名与
!     pes_adiab_i 的 i_surf 入参据此撤销（表象与电子态知识不进甲类签名）。
!   B 黑箱性：甲类+乙类合为一个 system black box——箱外只见
!     (Q,P,s[,dt])→(F,E)（facade_force 签名零变化）；ws 是 facade_force
!     局部变量，门面内声明、门面内消亡。
!   C 物自体原则：H_ij/d_ij/V_i 等表象构件不具有物自体地位——不进 state_t、
!     不出箱外、生命周期一步（每步求值复位重建）。
! E9 保留资产：elec_t 振幅/跃迁计数、乙类 dt 通道、协商表数据化 triad_tbl
!（pes_kind 列语义改为 ws 填充模式 repr；radau 多面行恒假约束不变——
!   自适应序列边界"一步"无良定义，跃迁概率 g_ij(dt) 失据）。
!=====================================================================
module force_facade
   use consts,        only: e_conv
   use state,         only: state_t, elec_t
   use config,        only: electronic, propagator
   use config_roster, only: roster
   use elec_ws,       only: ws_t, r_scalar, r_diab, r_adiab
   implicit none
   private
   public :: facade_reg_pes
   public :: facade_reg_elec, facade_init, facade_force, facade_info

   ! eV→kcal/mol 换算（单位契约字面值——T4：全库仅本文件出现）
   real(8), parameter :: ev_kcal = 23.0605d0

   ! 甲类成员（可叠加）接口·E10 归一：(q,mass) → ws（体系/表象双不变——
   ! 单面标量与多面透热/绝热统一为 ws 填充模式；形态数组不进签名）
   abstract interface
      subroutine pes_proc_i(q, mass, ws)
         import :: ws_t
         real(8), intent(in)  :: q(:)     ! 坐标 [Å]（全自由度展平）
         real(8), intent(in)  :: mass(:)  ! 每原子质量 [amu]（随机律/质量律成员用；
                                          ! 纯势成员可忽略）
         type(ws_t), intent(inout) :: ws  ! 工作区（inout：门面预写请求字段
                                          ! i_want，成员读之并按自身模式填充
                                          ! 交付字段——scalar 填 g_add/v_add；
                                          ! diab 填 h_ij/dh_ij；adiab 填
                                          ! v_surf/g_surf/d_ij〔i_want 行〕；
                                          ! 见 elec_ws.f90 容器契约）
      end subroutine pes_proc_i
      ! 注：随机律成员内部调 rng（N1——签名承认随机律；单流由 rng 担保）
   end interface

   ! 乙类成员（单选）接口：电子传播律——电子态演化 + 核力修饰（绝热=两退化点）
   abstract interface
      subroutine elec_proc_i(s, f, dt, ws)
         import :: elec_t, ws_t
         type(elec_t), intent(inout) :: s  ! 电子态（律成员按律演化 amp/i_surf/
                                           ! n_hop；绝热：不变）
         real(8), intent(inout) :: f(:)    ! 核力 [内部能量/Å]（律成员修饰——面
                                           ! 选择/平均/切换与动能重整；绝热：恒等）
         real(8), intent(in)    :: dt      ! 当前实际步长 [10 fs]（跃迁概率
                                           ! g_ij(dt) 所需——E9；绝热成员忽略）
         type(ws_t), intent(in)  :: ws     ! 工作区（只读消费：按 ws%repr 解释
                                           ! 容器——箱内接线，E10；成员不支持
                                           ! 的 repr → loud-stop，如 adiab×iesh）
      end subroutine elec_proc_i
   end interface

   ! 注册表（甲类可叠加表 + 乙类单选槽；装配相位登记，演化期只读）
   integer, parameter :: max_pes = 16     ! 甲类成员表容量 [行]（主势+加性随机律；C2 填充）
   type :: pes_entry_t
      character(len=32) :: name = ''      ! 成员名（config%electronic%pes_main 装配键）
      integer :: kind = r_scalar          ! 成员 ws 填充模式（r_*——登记数据，
                                          ! E10：表象是数据不是接口分叉）
      procedure(pes_proc_i), pointer, nopass :: proc => null() ! 成员求值例程
                                          ! （单签名——scalar/diab/adiab 统一）
   end type pes_entry_t
   type(pes_entry_t) :: pes_tbl(max_pes)  ! 甲类注册表（模块私有）
   integer :: n_pes = 0                   ! 已登记甲类行数 [行]
   character(len=32) :: elec_name = ''    ! 乙类成员名（模块私有；单选槽；空=未装配）
   procedure(elec_proc_i), pointer :: elec_proc ! 乙类例程（模块私有；装配相位
                                                ! facade_init 挂槽；nopass 仅组件属性）

   ! 能力级合法组合表（E9 数据化；E10 pes_kind 列语义改 ws 填充模式——
   ! 加行=加数据不改代码）。四轴：H 能力级 n_surf 区间 × 主势 repr 模式 ×
   ! 电子律 × 积分器（积分器轴以三布尔承载）。radau 列在多面行恒假：自适应
   ! 序列边界"一步"无良定义，跃迁概率 g_ij(dt) 失据——非绝热只配定步长
   ! 积分器（verlet/symple）。
   integer, parameter :: m_adiabatic = 1  ! 电子律编号：绝热退化点（config%electronic%method 键）
   integer, parameter :: m_tsh      = 2   ! 电子律编号：TSH/FSSH（将来 elec_tsh 注册键）
   integer, parameter :: m_iesh     = 3   ! 电子律编号：IESH（将来 elec_iesh 注册键）
   integer, parameter :: max_triad = 4    ! 合法组合表容量 [行]
   integer, parameter :: n_surf_unbounded = huge(1) ! 多面行能力级上限（无界意）
   type :: triad_row_t
      integer :: n_surf_lo = 1            ! H 能力级下限（n_surf ≥）
      integer :: n_surf_hi = 1            ! H 能力级上限（单面行 1..1；多面行 ..无界）
      integer :: pes_kind  = r_scalar     ! 主势 ws 填充模式（r_*——E10 repr 语义）
      integer :: method    = m_adiabatic  ! 电子律编号（m_*）
      logical :: ok_verlet = .true.       ! 积分器 1（verlet/beeman）放行位
      logical :: ok_symple = .true.       ! 积分器 2（symple）放行位
      logical :: ok_radau  = .true.       ! 积分器 3（radau）放行位（多面行假）
   end type triad_row_t
   type(triad_row_t), parameter :: triad_tbl(max_triad) = (/ &
      ! 行 1（现行）：单面 scalar 模式 × 绝热 × 三积分器全放行
      triad_row_t(1, 1, r_scalar, m_adiabatic, .true., .true., .true.), &
      ! 行 2：多面 diab 模式 × TSH × 定步长积分器
      triad_row_t(2, n_surf_unbounded, r_diab,  m_tsh, .true., .true., .false.), &
      ! 行 3：多面 adiab 模式 × TSH × 定步长积分器
      triad_row_t(2, n_surf_unbounded, r_adiab, m_tsh, .true., .true., .false.), &
      ! 行 4：多面 diab 模式 × IESH × 定步长积分器（iesh 消费 diab 容器天然——模型哈密顿量；
      !       adiab×iesh 不放行）
      triad_row_t(2, n_surf_unbounded, r_diab,  m_iesh, .true., .true., .false.) &
      /)
   ! 注：行 2-4 的 method（tsh/iesh）当前无乙类成员注册——协商表放行在先，
   ! facade_init 乙类装配（按 method 查注册槽）loud-stop 兜底拦截，双闸门
contains
   !------------------------------------------------------------------
   ! facade_reg_pes(name, proc[, repr]) — 登记甲类力成员（可叠加：主势 + 随机律/
   !                  约束势；E10 归一单入口——表象经 repr 数据标记，不再是注册分叉）
   ! 参考: docs/final_completion_plan.md §5 N1（注册规格）；T2（新成员=新文件+
   !       一行注册）；E10 裁决（单签名完全归一——E9 双签名撤销）
   subroutine facade_reg_pes(name, proc, repr)
      character(len=*), intent(in) :: name  ! 成员名（装配键——config%electronic%pes_main）
      procedure(pes_proc_i) :: proc         ! 成员求值例程（单签名——scalar/diab/
                                            ! adiab 统一经 ws 容器交付）
      integer, intent(in), optional :: repr ! ws 填充模式 r_*（缺省 r_scalar——
                                            ! 六单面成员注册行零触及；多面主势
                                            ! 注册行带 r_diab/r_adiab）
      ! 伪代码：
      ! 1. 重名或表满 → loud-stop
      ! 2. 登记行 {name, kind=repr（present 判缺省 r_scalar）, proc}（装配相位
      !    执行——driver；多面 pes 成员随 elec_tsh 消费者引入）
   end subroutine facade_reg_pes

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
      !    多面体系（n_surf ≥ 2）：主势行 ws 填充模式须为 r_diab/r_adiab
      !    （r_scalar 主势配多面 → loud-stop 指明冲突）；加性成员限 r_scalar
      !    （随机律/约束势与电子态无关——只进核力通道，多面构件通道唯主势
      !    承载，见 facade_force 步 4'）
      ! 2. 乙类装配：按 config%electronic%method 选律成员（1=绝热 →
      !    elec_adiabatic 退化成员；2/3=tsh/iesh 将来注册）；未登记 → loud-stop
      !    （协商表放行在先、注册槽空此后拦——双闸门，triad_tbl 注）
      ! 3. 成员初始化钩子（POTPRE 型）：逐装配成员调其 init（参数经成员自带
      !    params——D3；钩子槽的接线形态 C2 定形）
      ! 4. 能力级组合校验（E9 查表化）：四元组
      !    (n_surf ∈ 行区间) × (主势 ws 模式 = 行 pes_kind) × (method = 行 method)
      !    × (积分器 id 放行位——config%propagator%id 映射 ok_verlet/ok_symple/
      !    ok_radau) 逐行匹配 triad_tbl；无命中行 → 打印冲突各轴现值并 STOP
      !    （规范 §3.2-2 类型检查——指明两轴；典型拦截：多面 × radau——
      !    跃迁概率 g_ij(dt) 在自适应序列边界失据）
   end subroutine facade_init

   !------------------------------------------------------------------
   ! facade_force(sta, e_tot[, dt]) — (Q,P,s)→(F,E)：成员经 ws 交付 → 门面合成 →
   !                            换算 → 电子律消费；E = T(P) + V（E10 统一路径）
   ! 参考: venus2026/DVDQ.f:10-22（清零+分派）；src_TEST/interface_TEST.f90:32-44
   !       （同形累加）、:18-30（能量口径）；ENERGY.f:54-68（T/H 报告口径）
   subroutine facade_force(sta, e_tot, dt)
      type(state_t), intent(inout) :: sta  ! 物理态（q/p/s 读入；f 写出——Σ合成力）
      real(8), intent(out) :: e_tot        ! 总能 H = T(P)+V [kcal/mol]（N1 口径；
                                           ! 多面模式 V=当前面/律成员口径，见 4'）
      real(8), intent(in), optional :: dt  ! 当前实际步长 [10 fs]（跃迁概率
                                           ! g_ij(dt) 通道——E9；缺省取
                                           ! config%propagator%dt——定步长语义；
                                           ! 绝热律忽略之）
      ! 门面局部工作区（E10 公理 B/C：箱内声明、箱内消亡、一步生命周期——
      ! 对箱外（driver/propagator/output）不可见，facade_force 签名零触及）
      type(ws_t) :: ws                     ! 电子工作区（主势多面构件 + scalar
                                           ! 成员贡献容器；成员填充/律消费共道）
      ! 伪代码：
      ! 1. 力核工作区 g(:) 清零、v_sum = 0（旧 DVDQ.f:12-14 PDOT 清零位）；
      !    ws%repr = r_empty（每步复位——公理 C 一步生命周期）；dt_eff = 缺省
      !    config%propagator%dt（prop_cfg 别名）或传入值
      ! 2. 主势求值（装配主势成员——单签名统一调用）：ws%i_want = sta%s%i_surf
      !    （adiab 模式当前面耦合行请求——i_surf 不进甲类签名，经请求字段
      !    传达，E10）；proc(sta%q, roster%mass, ws)——成员按自身模式填充 ws
      !    （scalar：g_add/v_add；diab：h_ij/dh_ij；adiab：v_surf/g_surf/
      !    d_ij〔i_want 行〕）
      ! 3. scalar 贡献并入 g(:)/v_sum：主势为 scalar 模式时并入其 g_add/v_add；
      !    逐加性成员（限 r_scalar 模式行——随机律/约束势与电子态无关，只进
      !    核力通道）：proc(sta%q, roster%mass, ws) 后并入（eV 累加；成员间
      !    门面复位 ws——交付不叠加、并入由门面显式执行；随机律成员内部调
      !    rng；多面构件通道唯主势承载）
      ! 4. 单位契约履行点（T4——全库唯一）：
      !    sta%f 预清零后逐元素 sta%f(i) = sta%f(i) − g(i)·ev_kcal·e_conv
      !    （与旧 DPESHON "PDOT(i)=PDOT(i)-g(i)*23.0605*C1" 逐字符同形——N1.1
      !    bit 级可比；e_conv 即旧 C1；多面模式 g 的构成见 4'）
      ! 5. 电子律消费：elec_proc(sta%s, sta%f, dt_eff, ws)（只读消费——按
      !    ws%repr 解释容器，成员不支持的 repr → loud-stop。绝热成员恒等并
      !    忽略 ws——单面表示下复合塌缩为纯 PES，规范 §3.2-2；非绝热律在此
      !    完成 amp 传播、跃迁判定（rng 单流）与 i_surf 切换后的力修饰/动能
      !    重整——消费细节随 elec_tsh/elec_iesh 成员定形（C 阶段））
      ! 6. 能量：v_pot = v_sum·ev_kcal [kcal/mol]（旧 POT0 :26 口径）；
      !    t_kin = Σ p²/(2·roster%mass)/e_conv [kcal/mol]（逐原子折叠——
      !    旧 ENERGY.f:58-66 口径）；e_tot = t_kin + v_pot（旧 :68 H=T+V）
      ! ————————————————————————————————————————————————
      ! 4'. 多面模式分支（主势行模式 = r_diab/r_adiab——步 3-6 的构成置换）：
      !    a. 力核构成：r_diab → ws%h_ij 经 utils/linalg 对角化现算绝热视图
      !       （本征值=各面 V、本征矢=表象变换 U），当前面绝热力核取变换后
      !       当前行，跃迁概率 g_ij(dt) 所需 dh_ij·U 由乙类成员自取（门面不
      !       预构）；r_adiab → g = ws%g_surf(:, sta%s%i_surf)（当前面行——
      !       hop 后乙类成员更新 i_surf，随新面行重算属下一求值步）
      !    b. v_pot 口径：adiab → v_surf(i_surf)·ev_kcal；diab → 布居加权或
      !       当前面（乙类律成员口径，消费者定形）；scalar 加性贡献照常并入；
      !       e_tot = t_kin + v_pot
      !    c. radau 注记：多面行协商表恒假（triad_tbl）——自适应子步无"一步"
      !       dt 语义，本路径只见定步长（propagator 头注记同源）
   end subroutine facade_force

   !------------------------------------------------------------------
   ! facade_info() — 诊断输出：成员清单与能力级组合
   ! 参考: docs/final_completion_plan.md §5 N1（facade_info——成员清单/能力级）
   subroutine facade_info()
      ! 伪代码：
      ! 1. 打印甲类注册/装配成员清单（名 + 主势/加性标注 + ws 填充模式
      !    r_*）与乙类成员名
      ! 2. 打印当前能力级组合：H 能力级 n_surf × 主势 ws 模式 × 电子律 × 积分器
      !    （config 三源）+ triad_tbl 命中行号（未装配则示空）
      ! 3. 打印单位契约常量（ev_kcal/e_conv）与装配状态——排查装配错用
   end subroutine facade_info
end module force_facade
