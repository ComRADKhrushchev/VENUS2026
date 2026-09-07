!=====================================================================
! sampler.f90 — L2 接口 B′ 框架：初态分布通用抽样器（分布成员注册 + 方案分发
!                ——分布→(Q0,P0,s0)；SELECT 整数分发改成员表）
! 参考: src_VENUS/SELECT.f90:71-87（顶部门控——读入直配坐标的旧任务分流）、
!       :89-142（片段位置 + DH 分离能参考——能量归 facade、几何归 roster，剥除）、
!       :144-186（片段 A/B 初始化分派——旧整数开关分发的成员表化位）、
!       :188-207（碰撞装配分派：表面/气相两路）、:1300-1334（setup_surface_collision）、
!       :1337-1404（setup_gas_collision——两路归并入射通道成员）、
!       :1407-1491（finalize_and_report——报告归 recorder，剥除）；
!       分布成员名册（X2/X3）：boltzmann/ebk/incident（批 3 建骨架）+ 余 9 个 → C5
!=====================================================================
module sampler
   use state, only: state_t
   implicit none
   private
   public :: samp_reg, samp_run

   ! 分布成员接口：按自带分布抽 (q0,p0,s0) 的成员负责分量写入物理态
   abstract interface
      subroutine dist_proc_i(sta)
         import :: state_t
         type(state_t), intent(inout) :: sta  ! 初态（成员负责的 q/p/s 分量被抽样改写）
      end subroutine dist_proc_i
   end interface

   ! 分布成员表（B′ 注册族——新增成员 = 新文件 + 一行注册，T2）
   integer, parameter :: max_dist = 16   ! 成员表容量 [行]
   type :: dist_entry_t
      character(len=32) :: name = ''     ! 成员名（方案键——samp_run 的 scheme）
      procedure(dist_proc_i), pointer, nopass :: proc => null() ! 成员抽样例程
   end type dist_entry_t
   type(dist_entry_t) :: dist_tbl(max_dist) ! 分布成员表（模块私有；C5 填充）
   integer :: n_dist = 0                   ! 已登记行数 [行]
contains
   !------------------------------------------------------------------
   ! samp_reg(name, proc) — 登记分布成员进成员表
   ! 参考: docs/final_completion_plan.md §7 X2/X3（SELECT 成员表化——12 成员名册）
   subroutine samp_reg(name, proc)
      character(len=*), intent(in) :: name  ! 成员名（如 'boltzmann'/'ebk'/'incident'）
      procedure(dist_proc_i) :: proc        ! 成员抽样例程
      ! 伪代码：
      ! 1. 重名或表满 → loud-stop
      ! 2. 登记行 {name, proc}（装配相位执行——driver；dist_* 成员各注册一行）
   end subroutine samp_reg

   !------------------------------------------------------------------
   ! samp_run(scheme, sta) — 按方案名分发抽样：调成员写初态分量
   ! 参考: src_VENUS/SELECT.f90:144-186（片段分派）、:188-207（碰撞装配分派）
   subroutine samp_run(scheme, sta)
      character(len=*), intent(in) :: scheme ! 抽样方案名（config%reactants%dist_scheme
                                             ! 所指成员名——旧整数分发→成员表，X3；
                                             ! 整数编号→成员名映射 C5 接线定形）
      type(state_t), intent(inout) :: sta    ! 初态（q0/p0/s0——被成员抽样改写）
      ! 伪代码：
      ! 1. 查分布成员表 dist_tbl 得 scheme 成员（未登记 loud-stop 指名）
      !    ——分布成员表 C5 填充（现仅框架：批 3 三示例 boltzmann/ebk/incident
      !    建骨架，余 9 个成员 surface_oscillator/rotation/j/barrier_excitation/
      !    thermostat/glo_target/thermalize/normalmode/qmmicro → C5——E3）
      ! 2. 调成员 proc(sta)：按自带分布抽其负责的初态分量（成员 = 代码 + 自带
      !    数据；内部可用 rng/linalg/geometry 仪器——制备成员内部调演化服务
      !    合法，边界 4）
      ! 3. 多片段体系依 config%reactants%dist_scheme 依序调各成员（片段 A/B
      !    分布）+ 入射通道成员（碰撞能/碰撞参数/取向——旧 setup_gas_collision
      !    :1337-1404 / setup_surface_collision :1300-1334 两路归并）
      ! 4. 旧 SELECT 总编排剥除：片段位置与 DH 分离能参考（:89-142）——能量经
      !    facade、几何经 roster；finalize_and_report（:1407-1491）——归 recorder
      !    统计列；sampler 只承担 分布→(q0,p0,s0)
   end subroutine samp_run
end module sampler
