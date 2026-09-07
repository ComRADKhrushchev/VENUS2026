!=====================================================================
! recorder.f90 — L4 记录：轨迹帧/终态的列式记录器（列 = 名字 + 抽取器的
!                注册机制；预设级或名单选列——体系知识经 roster/名单参数化）
! 参考: src_VENUS/GWRITE.f90:25-31（帧头：轨迹号/步号/时间 + T/V/H 能量行）、
!       :49-60（fort.1NNN 分子文件：步头 + E0(eV)/T(eV)/H(eV) 行 + 指定原子
!       Q/V 三元组）、:106-117（全原子 Q/P 展平段——GWRITE_LEVEL≥2）、
!       :122-127（逐原子力段——GWRITE_LEVEL≥5）、:133-143（NFQP 全坐标
!       stdout 段）、:145（FLUSH 检查点冲刷）；
!       src_VENUS/HIST.f90:8-24（record_all 终态记录列：Ntz/Nb/Site/RX0/RY0/
!       Zmin/Vfinal/Erel/Theta/Ephi——rec_final 列组来源；格位列经 surfcell
!       :19——lattice 传参化归 plot_hist 模块）；
!       src_VENUS/VENUS.f90:405/:450/:459（GWRITE 调用位——检查点间隔 NIP 的
!       旧编排，新版经 propagator 检查点回调）、:572（CLOSE(1000+NTZ) 单元
!       收尾——编排归 driver）；
!       src_TEST/interface_TEST.f90:18-30（POT0 能量口径——能量列抽取器口径）
! 剥离警示（GWRITE 四处体系知识——参考勿带入，K2 计划取材注）：
!   (1) 原子标签写死：'Q(C)='/'V(C)=' 以原子 1 为体系指定原子（:59-60）
!       → 观测对象经名册片段/名单指定；
!   (2) 质量阈元素名：W(i)<50 → 'C '/'Au ' 标签（:110-114）→ 标签取
!       roster%symb；
!   (3) nearest-Au 诊断：指定原子到最近 Au 的距离/坐标/分量扫描（:68-99）
!       → 最近邻诊断列参数化（目标原子与候选集经名册指定，C8 定形）；
!   (4) PBC 盒写死：由 QZB(1,·) 极值 + 6.0 Å padding 推 box（:69-87）→
!       盒参数经 roster 晶胞字段（box_lx/box_ly/skew）；
!   另：E0/23.0605 eV 显示换算（:57）不在本模块散写——eV 口径归单位契约
!   （force_facade 的 ev_kcal 唯一点，T4 同源）；旧随机数回显段（:30）与
!   旧温度回显段（:33-45）属旧任务分支——不迁
!=====================================================================
module recorder
   use state,         only: state_t
   use config,        only: observables
   use config_roster, only: roster
   use force_facade,  only: facade_force
   implicit none
   private
   public :: col_val_i, rec_reg_col, rec_select, rec_frame, rec_final

   ! 列抽取器接口：由物理态抽取一个标量观测值（测量类消费演化服务——能量列
   ! 内部经 facade_force 求势能，方向合法；规范 §2.2 测量类定义）
   abstract interface
      real(8) function col_val_i(sta)
         import :: state_t
         type(state_t), intent(in) :: sta   ! 物理态（只读快照）
      end function col_val_i
   end interface

   ! 列注册表（列 = 名字 + 抽取器；装配期注册，记录期只读）
   integer, parameter :: max_col = 64        ! 注册表容量 [列]
   type :: col_entry_t
      character(len=16) :: name = ''         ! 列名（obs_list 名单键；1-16 字符）
      procedure(col_val_i), pointer, nopass :: val => null() ! 抽取器
   end type col_entry_t
   type(col_entry_t) :: col_tbl(max_col)     ! 列注册表（模块私有）
   integer :: n_col = 0                      ! 已注册列数 [列]
   integer :: sel_idx(max_col) = 0           ! 已选列的注册表索引 [—]
   integer :: n_sel = 0                      ! 已选列数 [列]（rec_select 装配）
   logical :: sec_qp  = .false.              ! 全原子 Q/P 展平段门（rec_level≥2）
   logical :: sec_frc = .false.              ! 逐原子力段门（rec_level≥5）

   ! 选列泛型两入口：预设级 or 名单（GWRITE_LEVEL 四级语义 GWRITE.f90:11-15
   ! 与 config%observables%rec_level/obs_list 两键——K1 §1.4 观测设置槽）
   interface rec_select
      module procedure rec_sel_lvl   ! 按预设级（rec_level 键）
      module procedure rec_sel_lst   ! 按名单（obs_list 键——逗号分隔列名）
   end interface rec_select
contains
   !------------------------------------------------------------------
   ! rec_reg_col(name, proc) — 注册观测列（名字 + 抽取器各一行）
   ! 参考: docs/final_completion_plan.md §5（列注册机制——C4 记录器）；
   !       interface_TEST.f90:18-30（能量列口径——框架预注册列的参照）
   subroutine rec_reg_col(name, proc)
      character(len=*), intent(in) :: name  ! 列名（obs_list 名单键）
      procedure(col_val_i) :: proc          ! 列抽取器
      ! 伪代码：
      ! 1. 重名或表满 → loud-stop（指名列名）
      ! 2. 登记行 {name, proc}（装配相位执行——driver/成员 init 各注册一行）
      ! 3. 框架预注册列（本模块装配时先调）：t_kin（Σp²/2·roster%mass/e_conv，
      !    旧 T 口径 ENERGY.f:58-66）、v_pot/e_tot（经 facade_force(sta, e)
      !    ——observables%e_kin/e_pot/e_tot 开关各对应一列，旧 GWRITE.f90:31
      !    的 T/V/H 行）；成员/分析量列（e_rel、z_min、n_bnc 等）随成员与
      !    analysis 各自注册——C8 接线
   end subroutine rec_reg_col

   !------------------------------------------------------------------
   ! rec_select(level / list) — 选列泛型（两入口见模块声明部 interface 块）
   !
   ! rec_level 预设列组数据表（config%observables%rec_level 键——旧
   ! GWRITE_LEVEL 语义 GWRITE.f90:11-15 的参数化重述）：
   !   级 0 = 最小：帧头（轨迹/步/时间）+ t_kin/v_pot/e_tot 能量行（旧 :28/:31）
   !   级 1 = +结构诊断列（最近邻距离/坐标等——旧 :63-99 诊断段；体系知识
   !          参数化后入列，剥离警示 (3)(4)）
   !   级 2 = +全原子 Q/P 展平段 sec_qp（默认——旧 :106-117，向后兼容级；
   !          标签取 roster%symb——剥离警示 (2)）
   !   级 5 = +逐原子力段 sec_frc（旧 :122-127）

   !------------------------------------------------------------------
   ! rec_sel_lvl(level) — 按预设级选列（置 sel_idx 与两段门）
   ! 参考: GWRITE.f90:11-15（四级语义）
   subroutine rec_sel_lvl(level)
      integer, intent(in) :: level   ! 预设级 [—]（0/1/2/5；缺省调用传
                                     ! observables%rec_level）
      ! 伪代码：
      ! 1. 依上方数据表置 sel_idx（级 n 含级 <n 的全部列）与段门
      !    sec_qp/sec_frc（级 ≥2/≥5 开）
      ! 2. 表中列未注册 → loud-stop 指名（装配次序错——注册先于选列）
   end subroutine rec_sel_lvl

   !------------------------------------------------------------------
   ! rec_sel_lst(list) — 按名单选列（逗号分隔列名逐名查注册表）
   ! 参考: config%observables%obs_list（K1 §1.4——名单键）
   subroutine rec_sel_lst(list)
      character(len=*), intent(in) :: list  ! 列名名单（逗号分隔；
                                            ! 名单含保留段名 'all_qp'/'all_frc'
                                            ! 时开对应展平段——C4 定形）
      ! 伪代码：
      ! 1. 逗号分段逐名查 col_tbl（未注册 loud-stop 指名）
      ! 2. 命中行索引依序入 sel_idx、计 n_sel；段名按保留字开段门
   end subroutine rec_sel_lst

   !------------------------------------------------------------------
   ! rec_frame(sta, unit) — 写一检查点帧（逐已选列取值 + 门控展平段）
   ! 参考: GWRITE.f90:25-31（帧头/能量行）、:49-60（分子文件步段）、
   !       :106-127（两展平段）、:145（flush）
   subroutine rec_frame(sta, unit)
      type(state_t), intent(in) :: sta   ! 物理态（只读快照——检查点回调传入）
      integer, intent(in) :: unit        ! 输出单元号 [—]（driver 编排——旧
                                         ! 6/8/1000+NTZ 族，单元分配 C4 定形）
      ! 伪代码：
      ! 1. 帧头：轨迹号 traj_i/步号 step_i/时间 sta%t（旧 :27-28——序号通道
      !    经闭包宿主或可选哑元传入，骨架签名按计划两哑元——C4 定形）
      ! 2. 逐已选列 j = 1..n_sel：val = col_tbl(sel_idx(j))%val(sta)，
      !    依列格式化写一行/一节（能量行口径旧 :31）
      ! 3. sec_qp 开：全原子逐行 q/p 六元（标签 roster%symb——剥离警示 (2)；
      !    旧 :109-115）
      ! 4. sec_frc 开：逐原子力三元（旧 :122-127；sta%f）
      ! 5. flush(unit)（旧 :145）
   end subroutine rec_frame

   !------------------------------------------------------------------
   ! rec_final(sta, unit) — 写轨迹终态行（record_all 列组——统计相位消费的
   !                        逐轨迹记录）
   ! 参考: HIST.f90:8-24（record_all：终态行各列 + 首行表头）、:18（表头列序
   !       Ntz/Nb/Site/RX0/RY0/Zmin/Vfinal/Erel/Theta/Ephi）
   subroutine rec_final(sta, unit)
      type(state_t), intent(in) :: sta   ! 物理态（终态快照——只读）
      integer, intent(in) :: unit        ! 输出单元号 [—]
      ! 伪代码：
      ! 1. traj_i = 1 时写表头行（旧 HIST.f90:18——列名序）
      ! 2. 终态行各列取值：e_rel/thta/z_min/v_fin/n_bnc/x0/y0/site——终态量
      !    经 analysis/final_state 的 fin_t（列抽取器注册制——C8 接线；
      !    旧 record_all 直读全局 AErelS 等的成员化）；格位列经 plot_hist
      !    的格位归属工具（lattice 参数化——surfcell HIST.f90:313-373 的
      !    传参化形态，模块间接线 C8）
      ! 3. flush(unit)；本例程只写行——单元开闭与文件名编排归 driver
      !    （旧 VENUS.f90:572 CLOSE 位）
   end subroutine rec_final
end module recorder
