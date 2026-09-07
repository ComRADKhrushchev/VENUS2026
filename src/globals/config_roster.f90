!=====================================================================
! config_roster.f90 — 表象层：原子清单（roster）——体系构成表（质量/符号/片段表/晶胞），
!                     含装配 roster_load 与清空 roster_clear
! 参考: venus_input.f90:248-256（质量段）/:508-518（片段 A 表装配）/:743-766（片段 B 表，
!       B 侧 C 阶段定形）/:996-1000（晶胞参数）；src_VENUS/SURF.f:58-83（表面晶胞推导）；
!       src_VENUS/SELECT.f90:713-715、1028-1030（平衡几何 QZ 装配）；
!       venus_data.f90:349、585-586（晶胞/晶格常数声明与默认）
!=====================================================================
module config_roster
   implicit none
   private
   public :: roster_t, roster, roster_load, roster_clear

   ! 原子清单（roster）：体系构成表——每原子的元素符号、质量、片段归属，
   ! 及片段几何与晶胞参数（K1 计划术语约定；字段按 K1 计划 §1.4 登记名）
   type :: roster_t
      real(8), allocatable :: mass(:)          ! 每原子质量 [amu]（旧 W；刚性表面原子
                                               ! 置大数 ~1d30 固定——venus_input.f90:752）
      character(len=4), allocatable :: symb(:) ! 每原子元素符号 [—]（新增字段；旧库只存
                                               ! 质量，符号直读输入——文法 C4 定形）
      real(8) :: mass_a = 0.0d0                ! 片段 A 总质量 [amu]（旧 WTA(1) = Σ W(LA)）
      integer, allocatable :: list_a(:)        ! 片段 A 原子索引表 [全局原子序]（旧 LA(1,:)）
      real(8), allocatable :: qz_a(:)          ! 片段 A 平衡几何（主轴系坐标）[Å]
                                               ! （旧 QZA(1,:)，长 3*nat_a）
      integer :: nat_a = 0                     ! 片段 A 原子数 [个]（旧 NATOMA(1)）
      real(8), allocatable :: q0(:)            ! 全体系平衡几何 [Å]（旧 QZ；由片段几何
                                               ! 按原子索引装配——SELECT.f90:713-715）
      real(8) :: skew = 0.0d0                  ! 晶胞斜角 [rad]（旧 SKEW；输入度→弧度，
                                               ! venus_input.f90:999）
      real(8) :: box_lx = 0.0d0                ! 晶胞 x 边长 [Å]（旧 BOXLX）
      real(8) :: box_ly = 0.0d0                ! 晶胞 y 边长 [Å]（旧 BOXLY）
      real(8) :: box_pad = 2.0d0               ! 簇构型随机落点半宽 [Å]（旧 RND_BOX，
                                               ! venus_input.f90:307）
      real(8) :: a_lat = 2.94998d0             ! 表面晶格常数 [Å]（旧 auro，
                                               ! venus_data.f90:585 默认值）
      ! 留白：B 片段表（nat_b/list_b/mass_b/qz_b——旧 NATOMB/LB/WTB/QZB，含表面
      ! 晶胞四原子）与多片段扩展（K1 风险注：max_frag 随原子清单扩数组）——
      ! C 阶段随表面成员/入射通道定形
   end type

   type(roster_t) :: roster    ! 原子清单实例（同 config 模块的实例模式；roster_load 装配，
                               ! state%natoms = size(roster%mass)——接线 C 阶段 driver/input）
contains
   !------------------------------------------------------------------
   ! roster_load() — 装配与校验原子清单（已解析输入数据 → 质量表/片段表/平衡几何/晶胞）
   ! 参考: venus_input.f90:248-256、508-518、743-766、996-1000；
   !       src_VENUS/SURF.f:58-83；src_VENUS/SELECT.f90:713-715、1028-1030
   subroutine roster_load()
      ! 伪代码：
      ! 1. 质量表：按原子数分配 mass/symb；缺省广播规则（只给 2 值则第 2 值
      !    广播到其余原子，旧 :248-256）；符号填充留 C4
      ! 2. 片段 A 表：list_a = 1..nat_a 顺排、mass_a = Σ mass(list_a)、
      !    qz_a 读入平衡几何（旧 :508-518）
      ! 3. 表面刚性范式：表面原子质量置大数 ~1d30 固定（旧 :752）；
      !    B 片段表留白（C 阶段，旧 :743-766 仅作对照）
      ! 4. q0 装配：各片段平衡几何按原子索引落位成全体系几何
      !    （旧 SELECT.f90:713-715 由 A 填、:1028-1030 由 B 填）
      ! 5. 晶胞：周期盒范式直读 skew（度→弧度）/box_lx/box_ly（旧 :996-1000）；
      !    刚性表面范式由表面四原子坐标推 box_lx/box_ly/skew（余弦定理，
      !    旧 SURF.f:58-83；退化晶胞报错停机，旧 :74-79）——B 表 C 阶段接入
      ! 6. 校验：nat_a ≤ 原子数、list_a 索引不越界、质量全正（非法 loud-stop）
   end subroutine roster_load

   !------------------------------------------------------------------
   ! roster_clear() — 清空原子清单（释放全部 allocatable 成员并复位标量默认）
   ! 参考: venus_data.f90:585-586（晶胞/落点/晶格常数默认值即复位目标）
   subroutine roster_clear()
      ! 伪代码：
      ! 1. 释放 mass/symb/list_a/qz_a/q0
      ! 2. 标量复位：mass_a=0、nat_a=0、skew=0、box_lx=0、box_ly=0、
      !    box_pad=2.0、a_lat=2.94998（旧默认值，venus_data.f90:585-586）
   end subroutine roster_clear
end module config_roster
