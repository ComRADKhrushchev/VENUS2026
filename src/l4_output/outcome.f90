!=====================================================================
! analysis/outcome.f90 — L4 分析：结局分类（outcome_tbl 数据表驱动——逐行
!                          谓词求值；旧双份结局计数块的数据化）
! 参考: src_VENUS/VENUS.f90:252-284（结局汇总打印：NONREACT/NREACT(IPATH)/
!       NLONG/NEFF——统计相位计数表原型）、:552-571 与 :645-657（双份
!       反应/非反应计数块——几乎全同，数据化为表行）、:559-565（能量守恒门
!       DELM——漂移阈谓词原型）、:669-685（长时程 NLONG 块——寿命阈谓词）、
!       :669-673（长时程终态行的 -1 无效位——数据列语义）；
!       src_VENUS/TEST.f:106-146（NPATH 结局归属判据原型——片段质心距离；
!       NTEST 三值语义消解见 predicate.f90 注 5）；
!       config%event%outcome_tbl（K1 §1.4 事件规定——"数据非代码"；C7 接入）
! 剥离：旧 NPATH 整数结局（1=非反应、>1=反应路径号）与 NTEST/NFINAL 出口
!       语义——全部数据化为结局表的谓词行；结局名即数据，不枚举进代码
!=====================================================================
module outcome
   use config,      only: event
   use final_state, only: fin_t
   implicit none
   private
   public :: outc_expr_t, outc_row_t, outc_tbl_t, outc_classify
   public :: o_e_rel, o_b_fin, o_n_a, o_n_b, o_e_drift, o_t_life, o_and, o_or, o_not

   ! 谓词叶子观测键（outc_expr_t%key 取值；封闭集——新键属文法版本化事项）
   integer, parameter :: o_e_rel   = 1  ! 终态相对平动能 e_rel [kcal/mol]（散射结局）
   integer, parameter :: o_b_fin   = 2  ! 终态碰撞参数 b_fin [Å]
   integer, parameter :: o_n_a     = 3  ! 产物 A 原子数 [个]（片段重组结局——旧 NPATH 归属）
   integer, parameter :: o_n_b     = 4  ! 产物 B 原子数 [个]
   integer, parameter :: o_e_drift = 5  ! 能量守恒漂移 H_max−H_min [kcal/mol]（旧 DELM 门
                                        ! VENUS.f90:559-565——值来自 recorder 能量列轨迹聚合，
                                        ! 通道 C8 定形）
   integer, parameter :: o_t_life  = 6  ! 轨迹寿命 [fs]（旧 NLONG 门 :669-685）
   integer, parameter :: o_and     = 7  ! 组合 AND(e1,e2)
   integer, parameter :: o_or      = 8  ! 组合 OR(e1,e2)
   integer, parameter :: o_not     = 9  ! 组合 NOT(e)

   ! 结局谓词表达式节点（递归树——结构同 predicate%expr_t，叶子键改终态量域）
   type :: outc_expr_t
      integer :: key = 0           ! 观测键/组合键（o_*）
      integer :: cmp = 0           ! 比较方向 [—]（1=≥、2=≤；组合键不用）
      real(8) :: th  = 0.0d0       ! 阈值 [键的单位]
      type(outc_expr_t), allocatable :: sub(:)  ! 子式（and/or=2 个、not=1 个、叶子无）
   end type outc_expr_t

   ! 结局表行：结局名 + 行条件（表序即优先级——首真行中选）
   type :: outc_row_t
      character(len=32) :: name = ''  ! 结局名（如 'nonreact'/'react'/'long'——数据）
      type(outc_expr_t) :: cond       ! 行条件谓词
   end type outc_row_t

   type :: outc_tbl_t
      type(outc_row_t), allocatable :: row(:)  ! 行表
      integer :: n_row = 0                       ! 行数 [行]
      character(len=32) :: name_def = 'unclassified'  ! 全假默认名 [—]
   end type outc_tbl_t
contains
   !------------------------------------------------------------------
   ! outc_classify(fin, tbl) — 结局分类：逐行谓词求值，首真行名中选
   ! 参考: VENUS.f90:552-571/:645-657/:669-685（三类计数块的数据化）；
   !       TEST.f:106-146（react/nonreact 判据原型）
   function outc_classify(fin, tbl) result(name)
      type(fin_t), intent(in) :: fin        ! 终态量（fin_eval 产物——叶子取值域）
      type(outc_tbl_t), intent(in) :: tbl   ! 已解析结局表（表源 =
                                             ! config%event%outcome_tbl——装配
                                             ! 相位解析一次；串→表的解析器随
                                             ! 谓词文法 C4/C7 定形，与
                                             ! predicate 的 pred_parse 同族，
                                             ! 本阶段不另立例程）
      character(len=32) :: name             ! 函数值：结局名 [—]
      ! 伪代码：
      ! 1. 逐行 r = 1..tbl%n_row（表序即优先级）：递归求值 cond——
      !    叶子按 key 取值：o_e_rel→fin%e_rel、o_b_fin→fin%b_fin、
      !    o_n_a/o_n_b→产物原子数（名册片段计数——片段重组的原子归属分析
      !    C7 定形，旧 NPATH 由 TEST 检测成员变化）、o_e_drift→守恒漂移
      !    （recorder 能量列聚合——C8 通道）、o_t_life→fin%t_life；
      !    按 cmp（≥/≤）与 th 比较
      ! 2. 组合键递归子式：o_and→两真、o_or→任一真、o_not→取反
      !    （求值序=短路序，同 predicate%pred_eval 注 2）
      ! 3. 首个真行 → name = tbl%row(r)%name，返回（旧 react/nonreact 两分
      !    即两行表的最小例）
      ! 4. 全假 → name = tbl%name_def（长时程 'long' 行可作显式行或默认行——
      !    旧 NLONG 块 :669-685 对应）
   end function outc_classify
end module outcome
