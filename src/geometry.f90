!=====================================================================
! geometry.f90 — L5 仪器：几何变换（取向旋转/质心分离/角速度与角动量/对齐步进/求积节点）
! 参考: src_VENUS/ROTATE.f90:19-54（欧拉随机取向）；src_VENUS/ROTATEJKM.f:21-108（(J,M) 取向）；
!       src_VENUS/ArbitraryAxisRotation.f:18-75（任意轴旋转）；src_VENUS/CENMAS.f90:16-50（质心分离）；
!       src_VENUS/ANGVEL.f:9-20（角速度扣除）；src_VENUS/ROTN.f90:18-104（角动量/转动能）；
!       src_TEST/interface_TEST.f90:76-108（刚体对齐步进 Z/Y）；src_VENUS/GLPAR.f:23-69（Gauss-Legendre）
!=====================================================================
module geometry
   use consts, only: pi, two_pi
   use state,  only: state_t
   implicit none
   private
   public :: rot_euler, rot_jkm, rot_axis, cenmas, angvel, amom, align_step, gl_nodes

   ! 质心系工作区（派生量不入 state——K1 决策 D2；旧 QQ/PP 公共块内部化为模块私有；
   ! 尺寸随片段原子数增长，生命周期 C 阶段定形）：
   real(8), allocatable :: q_cm_frm(:)  ! 质心系坐标 [Å]（旧 QQ；cenmas 建立，rot_*/angvel/amom 消费）
   real(8), allocatable :: p_cm_frm(:)  ! 质心系动量 [amu·Å/(10 fs)]（旧 PP）
contains
   !------------------------------------------------------------------
   ! rot_euler(st, mass, list) — 片段绕质心按随机欧拉角整体旋转（坐标与动量同变换）
   ! 参考: src_VENUS/ROTATE.f90:19-54
   subroutine rot_euler(st, mass, list)
      type(state_t), intent(inout) :: st  ! 物理态（q/p 的 list 内分量被旋转改写）
      real(8), intent(in) :: mass(:)      ! 每原子质量 [amu]（roster%mass，按全局原子序索引）
      integer, intent(in) :: list(:)      ! 片段原子索引表 [全局原子序]（旧 LL）
      ! 伪代码：
      ! 1. 片段内原子做质心分离（复用 cenmas 步骤，得质心系工作区）
      ! 2. 抽 u 三次（rng_u）：phi = two_pi*u、costht = 2*u-1、chi = two_pi*u；
      !    thta = acos(costht)（旧 :19-25）
      ! 3. 由 (phi,thta,chi) 展开 3×3 旋转矩阵九元（旧 :26-39）
      ! 4. 质心系坐标与动量逐原子右乘该矩阵，写回 st%q/st%p 并同步工作区（旧 :40-54）
   end subroutine rot_euler

   !------------------------------------------------------------------
   ! rot_jkm(st, mass, list, j_vec, j_rot, m_proj) — 按 (J,M) 取向旋转片段
   !（J 旋到 z 轴 + 投影角 + 绕 J 随机方位 + 绕法向随机）
   ! 参考: src_VENUS/ROTATEJKM.f:21-108；src_VENUS/ROTATEY.f:13-29；
   !       src_VENUS/ROTATEZ.f:13-27；src_VENUS/ArbitraryAxisRotation.f:63-72
   subroutine rot_jkm(st, mass, list, j_vec, j_rot, m_proj)
      type(state_t), intent(inout) :: st  ! 物理态（q/p 的 list 内分量被旋转改写）
      real(8), intent(in) :: mass(:)      ! 每原子质量 [amu]
      integer, intent(in) :: list(:)      ! 片段原子索引表 [全局原子序]
      real(8), intent(in) :: j_vec(3)     ! 片段角动量矢量 [amu·Å²/(10 fs)]（旧 AM(1:3)，amom 算出）
      integer, intent(in) :: j_rot        ! 转动量子数 J [无量纲]（旧 JROT）
      integer, intent(in) :: m_proj       ! 投影量子数 M [无量纲]（旧首哑元传 M 值）
      ! 伪代码：
      ! 1. 片段内原子做质心分离（同 rot_euler 前置）
      ! 2. 绕 z 旋转把 j_vec 旋入 xz 平面（旧 :21-41）
      ! 3. xz 平面内旋转把 j_vec 旋到 z 轴（旧 :43-62）
      ! 4. 绕 y 旋转使 j 与面法向（z）成角
      !    acos(m_proj/sqrt(j_rot*(j_rot+1)))（旧 :64-70，绕 y 步进同 ROTATEY）
      ! 5. 绕当前 j 轴（过原点、方向 (sin(thta),0,cos(thta))）随机方位
      !    phi = two_pi*u 旋转：逐原子坐标与动量调 rot_axis
      !    （旧 :72-103；旧版以度制传角，新版统一弧度）
      ! 6. 绕 z 轴均匀随机旋转使终分布对法向柱对称（旧 :105-108，绕 z 步进同 ROTATEZ）
   end subroutine rot_jkm

   !------------------------------------------------------------------
   ! rot_axis(pt_a, pt_b, pt_in, ang, pt_out) — 点绕任意轴线（pt_a→pt_b）旋转（Murray 闭式）
   ! 参考: src_VENUS/ArbitraryAxisRotation.f:18-75（闭式展开 :63-72；轴向单位化 :32-35）
   subroutine rot_axis(pt_a, pt_b, pt_in, ang, pt_out)
      real(8), intent(in)  :: pt_a(3)    ! 轴线起点坐标 [Å]（旧 p1）
      real(8), intent(in)  :: pt_b(3)    ! 轴线终点坐标 [Å]（旧 p2；轴线方向 = pt_b - pt_a）
      real(8), intent(in)  :: pt_in(3)   ! 旋转前坐标 [Å]（旧 x1）
      real(8), intent(in)  :: ang        ! 旋转角 [rad]（逆时针；旧 th 以度计，新版统一弧度）
      real(8), intent(out) :: pt_out(3)  ! 旋转后坐标 [Å]（旧 x2）
      ! 伪代码：
      ! 1. 轴向单位向量 (u,v,w) = (pt_b - pt_a)/|pt_b - pt_a|（旧 :32-35）
      ! 2. 按 Murray 闭式三分量展开写 pt_out（cos(ang)/sin(ang) 组合，
      !    点到轴的投影项 + 旋转项，旧 :63-72）
   end subroutine rot_axis

   !------------------------------------------------------------------
   ! cenmas(st, mass, list, q_cm, v_cm) — 质心位置/速度计算 + 质心系坐标动量分离
   ! 参考: src_VENUS/CENMAS.f90:16-50
   subroutine cenmas(st, mass, list, q_cm, v_cm)
      type(state_t), intent(in) :: st    ! 物理态（只读；list 内原子参与求和）
      real(8), intent(in) :: mass(:)     ! 每原子质量 [amu]
      integer, intent(in) :: list(:)     ! 片段原子索引表 [全局原子序]
      real(8), intent(out) :: q_cm(3)    ! 质心坐标 [Å]
      real(8), intent(out) :: v_cm(3)    ! 质心速度 [Å/(10 fs)]（旧名 VCM 实为速度 = Σp/Σm）
      ! 伪代码：
      ! 1. v_cm = Σp(list)/Σm、q_cm = Σ(m·q)(list)/Σm（旧 :21-37）
      ! 2. 质心系工作区写分离量：q_cm_frm = q - q_cm、p_cm_frm = p - m·v_cm
      !    （旧 :39-50；工作区供 rot_euler/rot_jkm/angvel/amom 复用——旧 QQ/PP）
   end subroutine cenmas

   !------------------------------------------------------------------
   ! angvel(st, mass, list, omega) — 从动量中扣除刚体转动贡献 m·(omega×r)
   ! 参考: src_VENUS/ANGVEL.f:9-20
   subroutine angvel(st, mass, list, omega)
      type(state_t), intent(inout) :: st  ! 物理态（st%p 的 list 内分量被扣除改写）
      real(8), intent(in) :: mass(:)      ! 每原子质量 [amu]
      integer, intent(in) :: list(:)      ! 片段原子索引表 [全局原子序]
      real(8), intent(in) :: omega(3)     ! 角速度 [rad/(10 fs)]（旧 WX/WY/WZ，amom 回传）
      ! 伪代码：
      ! 1. 前置：cenmas 已建立质心系工作区（旧 QQ）
      ! 2. 逐原子 st%p -= mass·(omega × r_cm)（三分量展开，旧 :14-16）
      ! 3. 质心系动量工作区同步为扣除后的 p（p_cm_frm = st%p，旧 :17-19）
   end subroutine angvel

   !------------------------------------------------------------------
   ! amom(st, mass, list, j_vec, e_rot, omega) — 角动量/转动能（含惯量张量求逆）
   ! 参考: src_VENUS/ROTN.f90:18-104
   subroutine amom(st, mass, list, j_vec, e_rot, omega)
      type(state_t), intent(in) :: st          ! 物理态（只读；经质心系工作区取 r、p）
      real(8), intent(in) :: mass(:)           ! 每原子质量 [amu]
      integer, intent(in) :: list(:)           ! 片段原子索引表 [全局原子序]
      real(8), intent(out) :: j_vec(3)         ! 角动量矢量 [amu·Å²/(10 fs)]（旧 AM(1:3)）
      real(8), intent(out) :: e_rot            ! 转动能 [kcal/mol]（旧 EROT = 内部单位/e_conv）
      real(8), intent(out), optional :: omega(3) ! 可选回传角速度 [rad/(10 fs)]（旧 WX/WY/WZ，angvel 消费）
      ! 伪代码：
      ! 1. 前置：cenmas 已建立质心系工作区
      ! 2. j_vec = Σ r×p（质心系三分量，旧 :25-33）；单原子片段 j_vec/e_rot 全零（旧 :18-24）
      ! 3. 双原子片段：e_rot = |J|²/(2·I)，I = Σ m·r²（线性特例，旧 :35-47）
      ! 4. 多原子：惯量张量六分量（旧 :52-69）→ 行列式 → 伴随法求逆（旧 :70-81）
      ! 5. omega = I⁻¹·j_vec（旧 :83-87）；e_rot = (omega·j_vec)/2/e_conv（旧 :104）
      ! 6. 行列式奇异（|det|<0.01）退回线性特例公式（旧 :88-102）
   end subroutine amom

   !------------------------------------------------------------------
   ! align_step(st, i_axis, ang, list) — 刚体对齐单步（绕 z 或绕 y、绕原点；坐标动量同变换）
   ! 参考: src_TEST/interface_TEST.f90:76-90（Z 步）；:94-108（Y 步）
   subroutine align_step(st, i_axis, ang, list)
      type(state_t), intent(inout) :: st  ! 物理态（q/p 的 list 内分量被旋转改写）
      integer, intent(in) :: i_axis       ! 转轴（3=绕 z，旧 ALIGN_ZSTEP；2=绕 y，旧 ALIGN_YSTEP）
      real(8), intent(in) :: ang          ! 绕 z：转角 phi [rad]；绕 y：目标方向余弦 cos(thta) [无量纲]
      integer, intent(in) :: list(:)      ! 对齐原子索引表（旧固定为片段前两原子——键轴）
      ! 伪代码：
      ! 1. 绕 z 步：c=cos(ang)、s=sin(ang)，(x,y)→(x·c+y·s, -x·s+y·c)
      !    （旧 :81-89 对坐标与动量平行施加）
      ! 2. 绕 y 步：c=ang、s=sqrt(1-c²)，(x,z)→(x·c-z·s, x·s+z·c)
      !    （旧 :98-107 对坐标与动量平行施加）
   end subroutine align_step

   !------------------------------------------------------------------
   ! gl_nodes(x_lo, x_hi, node, wts) — Gauss-Legendre 求积节点与权重（对称折半 + 牛顿精化）
   ! 参考: src_VENUS/GLPAR.f:23-50（节点/权重映射）；:55-69（勒让德多项式递推 PL）
   subroutine gl_nodes(x_lo, x_hi, node, wts)
      real(8), intent(in)  :: x_lo       ! 积分下限 [被积变量单位]
      real(8), intent(in)  :: x_hi       ! 积分上限 [被积变量单位]
      real(8), intent(out) :: node(:)    ! 求积节点 [被积变量单位]（n = size(node)）
      real(8), intent(out) :: wts(:)     ! 求积权重 [被积变量单位]（与 node 对齐）
      ! 伪代码：
      ! 1. n = size(node)；m = (n+1)/2（对称折半，旧 :23）
      ! 2. 各根初值 x0 = cos(pi*(i-0.25)/(n+0.5))（旧 :28）
      ! 3. 牛顿-拉夫逊精化 x0 ← x0 - P_n(x0)/dP_n(x0)，P_n 三项递推、
      !    dP_n = n*(x0*P_n - P_{n-1})/(x0²-1)，至 |Δ|<1e-13（旧 :32-35 与 PL :60-68）
      ! 4. 映射到 [x_lo,x_hi]：中点 ± 半宽·x0；权重 = 半宽/((1-x0²)·dP_n²)；
      !    对称半边镜像填充另一半（旧 :43-50）
   end subroutine gl_nodes
end module geometry
