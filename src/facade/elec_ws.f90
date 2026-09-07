!=====================================================================
! elec_ws.f90 — L2 门面内部接线契约：电子工作区 ws_t（表象构件的容器——
!               甲类成员填充、乙类成员消费；一步生命周期，每步求值重建）
! 参考: docs/plans/K2.md E10（用户裁决 2026-09-07：接口不变性/黑箱性/
!       物自体原则三公理）；docs/design_spec_completeness.md §3.2-2（两轴）
! 设计公理（E10 立法）：
!   A 接口不变性——甲/乙类成员签名对体系不变：签名只传容器 ws_t，
!     不传 h_ij/d_ij 等形态数组（面数、表象不进签名）；
!   B 黑箱性——甲类+乙类合为一个 system black box；箱外（driver/
!     propagator/output）只见 (Q,P,s[,dt])→(F,E)；ws 不出现于 facade_force
!     箱外签名；
!   C 物自体原则——H_ij/d_ij/V_i 等表象构件不具有物自体地位（随表象选择
!     而变的"为我们"之量）：不进 state_t、不出箱外、生命周期一步
!     （facade_force 局部声明、门面内消亡）。
! 本模块为成员可见的哑元类型契约（Fortran 约束：签名哑元类型必须对
! 成员模块可见）——依赖方向：force_facade 与 interface_electronic/ 成员
! 均 use 本模块（L3→L2 向下依赖合法）；本模块零依赖（仅 implicit none）。
!=====================================================================
module elec_ws
   implicit none
   private
   public :: ws_t
   public :: r_empty, r_scalar, r_diab, r_adiab

   ! ws 填充模式值（repr 标记——表象从接口分叉降级为数据标记，E10）
   integer, parameter :: r_empty  = 0  ! 空/未装配（门面每步求值前复位）
   integer, parameter :: r_scalar = 1  ! 单面标量模式（g_add/v_add——六单面成员）
   integer, parameter :: r_diab   = 2  ! 多面透热模式（h_ij/dh_ij；IESH 天然）
   integer, parameter :: r_adiab  = 3  ! 多面绝热模式（v_surf/g_surf/d_ij）

   ! 电子工作区：门面私有中间量容器（成员经签名读写；对箱外不可见）
   type :: ws_t
      integer :: repr = r_empty   ! 填充模式（r_*：成员填充时自报；门面与
                                  ! 乙类成员按此解释容器——表象是数据不是签名）
      integer :: i_want = 0       ! 请求字段（门面写入：裁剪请求——绝热模式
                                  ! 当前面耦合行 d_i,· 的 i；成员读。
                                  ! i_surf 不进甲类签名，经此字段传达）
      ! —— scalar 模式容器 ——
      real(8), allocatable :: g_add(:)  ! 力核 −dV/dq [eV/Å]（随机律为等效力核）
      real(8)              :: v_add     ! 能量贡献 [eV]
      ! —— diab 模式容器 ——
      real(8), allocatable :: h_ij(:,:)   ! 透热哈密顿矩阵 [eV]（n_surf 阶对称；
                                          ! 成员交付全矩阵；对称性成员自理）
      real(8), allocatable :: dh_ij(:,:,:) ! 矩阵元梯度 [eV/Å]（ndof×n_surf²；
                                          ! 大态空间压缩交付属文法版本化事项）
      ! —— adiab 模式容器 ——
      real(8), allocatable :: v_surf(:)   ! 各面能量 [eV]
      real(8), allocatable :: g_surf(:,:) ! 各面力核 [eV/Å]（ndof×n_surf）
      real(8), allocatable :: d_ij(:,:)   ! 当前面导数耦合 <φ_i|∇_q|φ_j> [1/Å]
                                          ! （ndof×n_surf——按 i_want 行交付，
                                          ! O(ndof·n_surf) 而非全对；能隙加权
                                          ! 变体交成员/消费者自理，契约认裸 d）
   end type ws_t
contains
   ! 本模块纯类型/常数契约，无例程（类型生命周期管理归 force_facade：
   ! 每步求值复位 repr=r_empty、按模式分配/复用容器——C2 随门面具体化定形）
end module elec_ws
