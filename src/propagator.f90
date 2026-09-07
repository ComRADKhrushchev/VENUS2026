!=====================================================================
! propagator.f90 — L3 封闭族：核传播子（verlet/symple/radau 三成员注册 +
!                   统一回调步进——积分器本体无体系知识）
! 参考: src_VENUS/VERLET.f90:31-79（velocity Verlet：力前算 :34-41、位置全步
!       :52-56、力后算+半步 :60-67）；:90-170（Beeman 三阶：位置预测 :98-103、
!       速度预测 :105-109、力暂存与重算 :111-125、速度校正 :163-168）；
!       src_VENUS/SYMPLE.f:14-86（Schlier 4/6/8 阶系数装载）、:88-100（交错
!       漂移-踢击主循环）；src_VENUS/RADAU.f:90-94（Radau 节点 hh）、:138-175
!       （w/u 权与三角系数 c/d/ra）、:183-195（首序列 6 迭代）、:197-311（子步
!       预测-差分修正）、:313-336（重启判据与序列推进）、:390-396（自适应步长
!       tp=dir·tcorr·(ss/hv)^(1/9)）；RADOUT :438-498（内嵌编排——剥除为 hooks
!       回调，G4）；
!       剥离清单（参考时勿带入）：VERLET.f90:44-50/:71-79/:151-161（ghost 朗之万）、
!       :130-148（GLO/LDFA 内嵌力）、:179-270（FRICFORCE/DIAG2 链）——随机律一律经
!       force_facade 成员汇总（P-D1/D4）；'CALLING RADOUT' 调试打印（RADAU.f:112/
!       :380/:428）与 PRINFO 注释残留不迁（G4）
!=====================================================================
module propagator
   use state,         only: state_t
   use config,        only: prop_cfg => propagator   ! 别名规避与模块同名冲突
                                                      ! （prop_cfg 即 config 实例 propagator）
   use config_roster, only: roster
   use force_facade,  only: facade_force
   implicit none
   private
   public :: hooks_t, prop_reg, prop_step
   public :: verlet_step, symple_step, radau_step

   ! 检查点回调接口（RADOUT 编排剥除位——G4：输出→recorder、终止→predicate）
   abstract interface
      subroutine hook_out_i(sta)             ! 输出钩子：一检查点的帧写出
         import :: state_t
         type(state_t), intent(in) :: sta    ! 物理态（只读快照）
      end subroutine hook_out_i
      logical function hook_term_i(sta)      ! 终止谓词求值点：真=终止本轨迹
         import :: state_t
         type(state_t), intent(in) :: sta    ! 物理态（只读快照）
      end function hook_term_i
   end interface

   ! 演化回调束（driver 装配；空指针=该类回调不挂——旧 RADOUT 两回调位的数据化）
   type :: hooks_t
      procedure(hook_out_i),  pointer, nopass :: out  => null() ! 输出钩子（旧 RADOUT→GWRITE 位）
      procedure(hook_term_i), pointer, nopass :: term => null() ! 终止谓词点（旧 RADOUT→TEST 位）
   end type hooks_t

   ! 积分器成员接口：推进一个步长（verlet/symple）或一个自适应序列（radau）
   abstract interface
      subroutine integ_proc_i(sta, hooks)
         import :: state_t, hooks_t
         type(state_t), intent(inout) :: sta   ! 物理态（q/p/t 被推进；f 为力工作区）
         type(hooks_t), intent(in)    :: hooks ! 检查点回调束（序列成员在边界调）
      end subroutine integ_proc_i
   end interface

   ! 封闭族注册表（三成员装配期各登记一行；表容量含余量）
   integer, parameter :: max_integ = 8        ! 注册表容量 [行]
   type :: integ_entry_t
      integer :: id = 0                       ! 积分器编号（config%propagator%id 的键）
      procedure(integ_proc_i), pointer, nopass :: proc => null() ! 成员步进例程
   end type integ_entry_t
   type(integ_entry_t) :: reg_tbl(max_integ)  ! 注册表（模块私有）
   integer :: n_reg = 0                       ! 已登记行数 [行]
contains
   !------------------------------------------------------------------
   ! prop_reg(id, proc) — 登记积分器成员进封闭族注册表
   ! 参考: docs/final_completion_plan.md §7 X4（积分器统一协议 + 三成员注册）
   subroutine prop_reg(id, proc)
      integer, intent(in) :: id           ! 积分器编号 [—]（1=verlet/2=symple/3=radau——
                                          ! config%propagator%id 的键）
      procedure(integ_proc_i) :: proc     ! 成员步进例程
      ! 伪代码：
      ! 1. 重复 id 或表满 → loud-stop（封闭族仅三成员，装配期各调一次）
      ! 2. 登记行 {id, proc}（driver 装配相位执行——E7）
   end subroutine prop_reg

   !------------------------------------------------------------------
   ! prop_step(sta, hooks) — 步进一分派：按 config%propagator%id 调已注册成员
   ! 参考: src_VENUS/RADAU.f:438-498（RADOUT 编排剥除→hooks 回调的旧形态）
   subroutine prop_step(sta, hooks)
      type(state_t), intent(inout) :: sta   ! 物理态（经成员推进一个步长/序列）
      type(hooks_t), intent(in)    :: hooks ! 检查点回调束（透传成员）
      ! 伪代码：
      ! 1. 按 prop_cfg%id 查注册表 reg_tbl（未注册 loud-stop 指名 id）
      ! 2. 调成员 proc(sta, hooks)：推 sta 前进 prop_cfg%dt
      !    （radau=一个自适应序列；成员内部按需多次调 facade_force）
      ! 3. 检查点协议（成员内履行，先 out 后 term）：
      !    hooks%out 挂了才调——一帧写出（recorder 接线）；
      !    hooks%term 挂了才调——真=终止本轨迹（终止旗标由 driver 侧闭包置位：
      !    hooks 为 nopass 指针，闭包宿主在 driver 持旗标——E7 演化环的终止通道）
   end subroutine prop_step

   !------------------------------------------------------------------
   ! verlet_step(sta, hooks) — velocity Verlet（order=1）/ Beeman 三阶（order=2）单步
   ! 参考: src_VENUS/VERLET.f90:31-79（VV）；:90-170（Beeman）
   subroutine verlet_step(sta, hooks)
      type(state_t), intent(inout) :: sta   ! 物理态（q/p/t 推进；f 力工作区）
      type(hooks_t), intent(in)    :: hooks ! 检查点回调束（步末调）
      ! 局部工作区留 C6 声明（旧 FT2/QTP1 力与坐标暂存——力差分是 Beeman 分支状态）
      ! 伪代码：
      ! 1. dt = prop_cfg%dt、mass = roster%mass、order = prop_cfg%order
      ! 2. order=1（velocity Verlet，旧 :31-79）：
      ! 2a. facade_force(sta, e)（步首力 F(t)——旧 :34-41 力位）
      ! 2b. p += f·dt/2（半步动量，旧 :39-41）
      ! 2c. q += (p/m)·dt（位置全步，旧 :52-56）
      ! 2d. facade_force(sta, e)（步末力 F(t+dt)，旧 :60-62 力位）
      ! 2e. p += f·dt/2（半步收尾，旧 :65-67）
      ! 3. order=2（Beeman 三阶，旧 :90-170）：
      ! 3a. 位置预测：q += (p/m)·dt + dt²·(4f(t)−f(t−dt))/(6m)（旧 :98-103）
      ! 3b. 动量预测：p += dt·(3f(t)−f(t−dt))/2（旧 :105-109）
      ! 3c. 力暂存 f(t)→f(t−dt)（旧 :111-114）；facade_force(sta, e)（力重算，旧 :123-125）
      ! 3d. 动量校正：p = m·(q−q_prev)/dt + dt·(2f(t+dt)+f(t))/(6m)
      !     （旧 :163-168；q_prev=3a 前的 q，旧 QTP1 暂存）
      ! 4. 旧 ghost 朗之万（:44-50/:71-79/:151-161）与 GLO/LDFA 内嵌力
      !    （:130-148、FRICFORCE :179-270）不迁——随机律经 facade_force 成员
      !    汇总（P-D1/D4），积分器本体无体系知识
      ! 5. sta%t += dt；步末检查点回调（先 out 后 term）
   end subroutine verlet_step

   !------------------------------------------------------------------
   ! symple_step(sta, hooks) — 4/6/8 阶辛几何单步（Schlier 系数交错漂移-踢击）
   ! 参考: src_VENUS/SYMPLE.f:14-86（系数装载）、:88-100（主循环）
   subroutine symple_step(sta, hooks)
      type(state_t), intent(inout) :: sta   ! 物理态（q/p/t 推进；f 力工作区）
      type(hooks_t), intent(in)    :: hooks ! 检查点回调束（步末调）
      ! 系数缓存留 C6 声明（旧 A(0:34) + FIRST 旗标——模块私有静态）
      ! 伪代码：
      ! 1. 首次调用按 order（4/6/8，config%propagator%order）装载 Schlier 系数
      !    a(0:2N)（N=5/8/17；对称后半镜像复制，旧 :14-86 的 FIRST 块）
      ! 2. 交错循环 i=0..N-1（旧 :88-97）：
      ! 2a. facade_force(sta, e)（踢击前力——旧 :94-96 力位）
      ! 2b. p += dt·a(2i)·f（动量踢击，旧 :91）
      ! 2c. q += dt·a(2i+1)·p/m（位置漂移，旧 :92）
      ! 3. 末拍：facade_force(sta, e)；p += dt·a(2N)·f（旧 :98-100）
      ! 4. sta%t += dt；步末检查点回调（先 out 后 term）
   end subroutine symple_step

   !------------------------------------------------------------------
   ! radau_step(sta, hooks) — Everhart RA15 一个 Gauss-Radau 序列（自适应步长）
   ! 参考: src_VENUS/RADAU.f:90-94（节点 hh）、:138-175（w/u 权与 c/d/ra）、
   !       :183-195（首序列 6 迭代 + 步首力）、:197-311（子步预测-差分修正）、
   !       :313-336（重启判据与序列推进）、:390-396（下一序列步长）；
   !       RADOUT :438-498（编排剥除→hooks，G4）
   subroutine radau_step(sta, hooks)
      type(state_t), intent(inout) :: sta   ! 物理态（q/p/t 推进；f 力工作区）
      type(hooks_t), intent(in)    :: hooks ! 检查点回调束（序列边界调）
      ! 常数与 B/G/E 系数阵缓存留 C6 声明（模块私有静态——序列间携带）
      ! 伪代码：
      ! 1. 首次调用装载常数：Radau 节点 hh(2:8)（旧 :92-94）、w/u 权（旧 :138-143）、
      !    三角形系数 c/d/ra（旧 :153-175）
      ! 2. 速度折叠 v = p/m（旧 :103-107）；首序列迭代数 6、后续序列 2
      !    （旧 :185、:419）
      ! 3. 序列内子步循环（迭代 × 7 子步，旧 :197-311）：
      ! 3a. 由 B 系数七阶展开预测子步坐标 yy（旧 :220-245）
      ! 3b. facade_force 于子步坐标求值（旧 :250-260 力位；f_j = f/m 折叠）
      ! 3c. 差分修正：g 系数逐阶差商 + B 系数三角累加（旧 :262-311）
      ! 4. 序列收尾：q/p 推进一序列长 tt（旧 :329-336）；末项 hv 估计
      !    （旧 :313-317）；重启判据 tp/tt≤1 则弃序列重启（旧 :319-327）
      ! 5. 自适应步长：tp = tcorr·(ss/hv)^(1/9)（旧 :390；指数 1/9 旧 :74；
      !    ss = 10^(−order) 精度参数，方向 dir 随 tf 符号——旧 :133-135）
      ! 6. 步长夹逼 dt_min/dt_max——prop_cfg 留白字段，C6 引入
      !    （K1 计划 §1.3 留白注；引入前以 ss 精度单控）
      ! 7. 序列边界检查点回调（先 out 后 term；旧 RADOUT :449-478 的
      !    GWRITE/TEST 调用位——ENERGY_1/FINAL/GFINAL 统计编排不迁，归统计
      !    相位）；旧 NIP 间隔输出（旧 :99-101 prtime）由检查点间隔参数承担——C6 定形
      ! 8. sta%t += 序列长 tt；term 真即收束本步返 driver
   end subroutine radau_step
end module propagator
