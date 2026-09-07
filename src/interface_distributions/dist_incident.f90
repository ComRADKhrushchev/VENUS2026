!=====================================================================
! dist_incident.f90 — L3 分布成员示例（B′ 注册族）：入射通道分布——碰撞能/
!                     碰撞参数/取向的规定与抽样（示例 3/3；表面/气相两路归并）
! 参考: src_VENUS/SURF.f:120-187（取向角规定：极角 P(θ)∝sinθ 抽样
!       cosθ=1−U(1−cosθ_max) :133-141、方位角 [0,2π) 均匀 :142-149、落点
!       晶胞内均匀 :153-169、瞄准圆盘 b=b_max·√U :171-187——手册 Eq. V.19/
!       V.24-27）、:211-217（簇盒随机落点）、:258-325（碰撞能三规定：热布居
!       GAMA(2) :261-266、固定 :278-292、密度加权高斯拒绝抽样 :295-325）、
!       :331-338（入射方向单位矢——指向表面）、:343-374（相对速度按质量折叠
!       进片段动量）、:451-458（DENS=v³·exp(−((v−v0)/α)²)）、
!       src_VENUS/SELECT.f90:1300-1334/:1337-1404（表面/气相两路装配归并依据
!       ——sampler.f90 已引）
! 注册：本成员进 sampler 分布成员表（samp_reg 一行——装配相位执行；C5 接线）
!=====================================================================
module dist_incident
   implicit none
   private
   public :: inc_params_t, inc_init

   ! 与 config 分工：碰撞能 e_rel [kcal/mol]、最大碰撞参数 b_max [Å]、初始
   ! 间距 r_sep [Å] 是体系级入射规定——已在 config%reactants（K1 §1.4），
   ! 成员参数不重复存；本成员只携带分布形状与取向规定（模式开关 + 形状参数），
   ! sample 定形时（C5）按模式从 config%reactants 取 e_rel/b_max 分配
   type :: inc_params_t
      integer :: n_e_rel  = 1      ! 碰撞能规定模式 [—]（旧 NREL：0=热布居、
                                   ! 1=固定值、2=密度加权高斯）
      real(8) :: t_trans  = 0.0d0  ! 平动温度 [K]（旧 TRANS——模式 0 热布居：
                                   ! E_rel=R·GAMA(2)·T 换算，SURF.f:261-266）
      real(8) :: v_width  = 0.0d0  ! 速度展宽 α（模式 2 密度加权 DENS 的宽度，
                                   ! 内部速度单位 [Å/(10 fs)]；旧 ALPHAA）
      integer :: n_thta   = 0      ! 极角规定模式 [—]（固定束流 θ=θ_max vs
                                   ! P(θ)∝sinθ 抽样——SURF.f:133-141）
      real(8) :: thta_max = 0.0d0  ! 极角上限 θ_max [rad]（旧 THTA——输入度→
                                   ! 弧度经 consts 的 dtor，文法 C4）
      integer :: n_chi    = 0      ! 方位角规定模式 [—]（固定 χ vs [0,2π)
                                   ! 均匀——SURF.f:142-149）
      real(8) :: chi      = 0.0d0  ! 方位角 χ [rad]（旧 CHI——度→弧度同上）
      integer :: n_b      = 0      ! 瞄准模式 [—]（固定 b=b_max vs b=b_max·√U
                                   ! 圆盘均匀——b² 均布 [0,b_max²]，手册
                                   ! Eq. V.19、SURF.f:171-187）
      integer :: n_aim    = 0      ! 落点模式 [—]（表面晶胞内均匀 vs 簇盒随机
                                   ! ——SURF.f:153-169/:211-217；盒半宽归
                                   ! roster%box_pad、晶胞边/斜角归 roster）
      ! 留白：Gauss 密度峰值缓存（旧 PVMAX 型）与拒绝抽样上界（旧 5·v1）、
      ! 束流/热型研究的规定组合预设——C5 成员具体化时定形
   end type
contains
   !------------------------------------------------------------------
   ! inc_init(p) — 初始化成员（校验模式组合、预计算取向常量）
   ! 参考: src_VENUS/SURF.f:120-150（取向角规定——见模块头）
   subroutine inc_init(p)
      type(inc_params_t), intent(in) :: p ! 成员参数（模式开关/形状参数）
      ! 伪代码：
      ! 1. 校验模式组合合法性（非法 loud-stop）：n_e_rel ∈ {0,1,2} 且模式 0
      !    需 t_trans > 0、模式 2 需 v_width > 0；n_thta/n_chi/n_b/n_aim
      !    取值合法；取向规定仅表面/簇构型有意义（构型范式来自
      !    config%reactants%surface_model——装配相位核对）
      ! 2. 预计算取向常量：sin/cos(θ_max) 与 sin/cos(χ)；固定束流入射方向
      !    单位矢 u = (−sinθ·cosχ, −sinθ·sinχ, −|cosθ|)（指向表面，
      !    SURF.f:331-338——速度沿 −z 保证入射）
   end subroutine inc_init
   ! —— sample 方向（C5 定形）——落点 (rx0,ry0) 按 n_aim 规定抽样（晶胞内
   !   均匀含斜角修正/簇盒 ±半宽）；瞄准 b = b_max·√U（或固定 b_max）、方位
   !   φ 均匀 [0,2π)；E_rel 按模式分配（0: GAMA(2) 热布居、1: 定值取
   !   config%reactants%e_rel、2: DENS=v³·exp(−((v−v0)/α)²) 拒绝抽样）；
   !   相对速度沿 u 折算动量按质量分写两片段动量（表面/气相两路同式归并）
   ! sample 例程：C5 成员具体化时定形（本阶段不声明）
   ! 注册行注记：samp_reg 登记本成员（名 'incident'）进 sampler 分布成员表
   ! ——C5 接线
end module dist_incident
