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
!=====================================================================
module force_facade
   use consts,        only: e_conv
   use state,         only: state_t, elec_t
   use config,        only: electronic, propagator
   use config_roster, only: roster
   implicit none
   private
   public :: facade_reg_pes, facade_reg_elec, facade_init, facade_force, facade_info

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

   ! 乙类成员（单选）接口：电子传播律——电子态演化 + 核力修饰（绝热=两退化点）
   abstract interface
      subroutine elec_proc_i(s, f)
         import :: elec_t
         type(elec_t), intent(inout) :: s  ! 电子态（律成员按律演化；绝热：不变）
         real(8), intent(inout) :: f(:)    ! 核力 [内部能量/Å]（律成员修饰；绝热：恒等）
      end subroutine elec_proc_i
   end interface

   ! 注册表（甲类可叠加表 + 乙类单选槽；装配相位登记，演化期只读）
   integer, parameter :: max_pes = 16     ! 甲类成员表容量 [行]（主势+加性随机律；C2 填充）
   type :: pes_entry_t
      character(len=32) :: name = ''      ! 成员名（config%electronic%pes_main 装配键）
      procedure(pes_proc_i), pointer, nopass :: proc => null() ! 成员求值例程
   end type pes_entry_t
   type(pes_entry_t) :: pes_tbl(max_pes)  ! 甲类注册表（模块私有）
   integer :: n_pes = 0                   ! 已登记甲类行数 [行]
   character(len=32) :: elec_name = ''    ! 乙类成员名（模块私有；单选槽；空=未装配）
   procedure(elec_proc_i), pointer :: elec_proc ! 乙类例程（模块私有；装配相位
                                                ! facade_init 挂槽；nopass 仅组件属性）
contains
   !------------------------------------------------------------------
   ! facade_reg_pes(name, proc) — 登记甲类力成员（可叠加：主势 + 随机律/约束势）
   ! 参考: docs/final_completion_plan.md §5 N1（注册规格）；T2（新成员=新文件+一行注册）
   subroutine facade_reg_pes(name, proc)
      character(len=*), intent(in) :: name  ! 成员名（装配键——config%electronic%pes_main）
      procedure(pes_proc_i) :: proc         ! 成员求值例程
      ! 伪代码：
      ! 1. 重名或表满 → loud-stop
      ! 2. 登记行 {name, proc}（装配相位执行——driver；pes_* 成员 C2 引入后
      !    各注册一行）
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
   ! facade_init() — 成员装配 + 能力级三元组校验（非法组合 loud-stop）
   ! 参考: docs/final_completion_plan.md §5 N1（协商表：合法三元组当前仅
   !       (single_surface, adiabatic, {verlet,symple,radau})）；
   !       src_TEST/interface_TEST.f90:12-16（POTPRE 型成员初始化钩子——D3
   !       ghost 参数装配位）；docs/design_spec_completeness.md §3.2-2
   !       （能力级与方法的合法性做类型检查，非法组合显式拒绝并指明冲突的两轴）
   subroutine facade_init()
      ! 伪代码：
      ! 1. 甲类装配：按 config%electronic%pes_main 从 pes_tbl 固化主势成员行
      !    + 加性成员行（加性成员表——config%electronic 留白字段，C2 随旧开关
      !    NGLO/NFC 映射填充，文法向后兼容）；装配表空（无主势成员）→ loud-stop
      ! 2. 乙类装配：按 config%electronic%method 选律成员（1=绝热 →
      !    elec_adiabatic 退化成员）；未登记 → loud-stop
      ! 3. 成员初始化钩子（POTPRE 型）：逐装配成员调其 init（参数经成员自带
      !    params——D3；钩子槽的接线形态 C2 定形）
      ! 4. 能力级三元组校验（协商表数据化）：
      !    (H 表示能力级 = config%electronic%n_surf) × (电子律 = config%electronic%method)
      !    × (积分器 = config%propagator%id ∈ {verlet, symple, radau})
      !    ——当前合法仅 (单面, 绝热, {verlet|symple|radau})；
      !    非法组合打印冲突两轴并 STOP（规范 §3.2-2 类型检查——指明两轴）
   end subroutine facade_init

   !------------------------------------------------------------------
   ! facade_force(sta, e_tot) — (Q,P,s)→(F,E)：Σ成员力核 → 换算 → 电子律修饰；
   !                            E = T(P) + ΣV
   ! 参考: venus2026/DVDQ.f:10-22（清零+分派）；src_TEST/interface_TEST.f90:32-44
   !       （同形累加）、:18-30（能量口径）；ENERGY.f:54-68（T/H 报告口径）
   subroutine facade_force(sta, e_tot)
      type(state_t), intent(inout) :: sta  ! 物理态（q/p/s 读入；f 写出——Σ合成力）
      real(8), intent(out) :: e_tot        ! 总能 H = T(P)+ΣV [kcal/mol]（N1 口径）
      ! 伪代码：
      ! 1. 力核工作区 g(:) 清零、v_sum = 0（旧 DVDQ.f:12-14 PDOT 清零位）
      ! 2. 逐装配甲类成员：proc(sta%q, roster%mass, g_add, v_add)；
      !    g(:) += g_add、v_sum += v_add（eV 累加；随机律成员内部调 rng）
      ! 3. 单位契约履行点（T4——全库唯一）：
      !    sta%f 预清零后逐元素 sta%f(i) = sta%f(i) − g(i)·ev_kcal·e_conv
      !    （与旧 DPESHON "PDOT(i)=PDOT(i)-g(i)*23.0605*C1" 逐字符同形——N1.1
      !    bit 级可比；e_conv 即旧 C1）
      ! 4. 电子律修饰：elec_proc(sta%s, sta%f)（绝热成员恒等——单面表示下
      !    复合塌缩为纯 PES，规范 §3.2-2）
      ! 5. 能量：v_pot = v_sum·ev_kcal [kcal/mol]（旧 POT0 :26 口径）；
      !    t_kin = Σ p²/(2·roster%mass)/e_conv [kcal/mol]（逐原子折叠——
      !    旧 ENERGY.f:58-66 口径）；e_tot = t_kin + v_pot（旧 :68 H=T+V）
   end subroutine facade_force

   !------------------------------------------------------------------
   ! facade_info() — 诊断输出：成员清单与能力级三元组
   ! 参考: docs/final_completion_plan.md §5 N1（facade_info——成员清单/能力级）
   subroutine facade_info()
      ! 伪代码：
      ! 1. 打印甲类注册/装配成员清单（名 + 主势/加性标注）与乙类成员名
      ! 2. 打印当前能力级三元组：H 表示 × 电子律 × 积分器（config 三源）
      ! 3. 打印单位契约常量（ev_kcal/e_conv）与装配状态——排查装配错用
   end subroutine facade_info
end module force_facade
