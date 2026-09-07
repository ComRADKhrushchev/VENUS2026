!=====================================================================
! predicate.f90 — L3 封闭族：终止谓词文法解释器（表达式树解析 + 求值；
!                  四原语 + 三布尔组合——体系判据退为数据）
! 参考: src_VENUS/TEST.f:88-94（Z_HEIGHT 写死判据原型：rd=Q(3) 离面高度 +
!       P(3)>0 外向速度，表面顶层 z=0）；src_VENUS/TEST.f:106-146（COM_DISTANCE
!       原型：质量加权质心累加 :107-130、质心距与退行速度点积 :133-143、
!       RBAR/RMAX 双阈判 :144-145）；docs/final_completion_plan.md §5 N2（文法与
!       接口规格——直接生效）；§3 B2（TEST 通用判据登记 → 本模块文法成员）
! 输入文法（TERMINATE= 表达式，param_mapping 注册的新关键字——N2 规格）：
!   表达式 := 原式 | AND(e1,e2) | OR(e1,e2) | NOT(e)
!   原式   := COM_DISTANCE(片段A, 片段B, 阈值Å[, 退行速度阈])  ← TEST.f:106-146 原型
!           | Z_HEIGHT(片段, 阈值Å[, OUTWARD])                 ← 取代 TEST.f:88-94 写死判据
!           | MAX_STEPS(步数)                                  ← 取代步数魔数（恒为安全网）
!           | BOUND(片段A, 片段B, 下阈, 上阈)                   ← 有界窗（约束散射）
!   片段   := 名册片段名（A/B/…）或 ATOM(i)
!=====================================================================
module predicate
   use state,         only: state_t
   use config_roster, only: roster
   implicit none
   private
   public :: frag_ref_t, expr_t
   public :: pred_parse, pred_eval
   public :: com_dist, z_height, bound, steps
   public :: pred_and, pred_or, pred_not
   public :: k_com_dist, k_z_height, k_bound, k_steps, k_and, k_or, k_not

   ! 节点类别编号（expr_t%node_kind 取值；封闭集——新增原式属文法版本化事项）
   integer, parameter :: k_com_dist = 1  ! 原式：COM_DISTANCE（片段对质心距离）
   integer, parameter :: k_z_height = 2  ! 原式：Z_HEIGHT（片段离面高度）
   integer, parameter :: k_bound    = 3  ! 原式：BOUND（有界窗越窗）
   integer, parameter :: k_steps    = 4  ! 原式：MAX_STEPS（步数安全网）
   integer, parameter :: k_and      = 5  ! 布尔组合：AND(e1,e2)
   integer, parameter :: k_or       = 6  ! 布尔组合：OR(e1,e2)
   integer, parameter :: k_not      = 7  ! 布尔组合：NOT(e)

   ! 片段引用（文法"片段"的树内表示：名册片段名或单原子）
   type :: frag_ref_t
      integer :: i_frag = 0   ! 名册片段号 [—]（A/B/… → 1/2/…；0=未设）
      integer :: i_atom = 0   ! 单原子序 [—]（ATOM(i) 文法；>0 时覆盖片段引用——
                              ! 单原子判据，旧 TEST.f:89 写死 atom 1 的通用化）
   end type frag_ref_t

   ! 谓词表达式节点（递归树：原式为叶节点，布尔组合为内节点）
   type :: expr_t
      integer :: node_kind = 0      ! 节点类别（k_com_dist…k_not；0=空/未装配）
      type(frag_ref_t) :: frag_a   ! 片段引用 A（com_dist/z_height/bound 用）
      type(frag_ref_t) :: frag_b   ! 片段引用 B（com_dist/bound 用）
      real(8) :: th_lo    = 0.0d0   ! 主阈/下阈 [Å]（com_dist 距离阈；z_height 高度阈；bound 下阈）
      real(8) :: th_hi    = 0.0d0   ! 上阈 [Å]（bound 专用）
      real(8) :: v_rec    = 0.0d0   ! 退行速度阈 [Å/(10 fs)]（com_dist 可选段）
      logical :: has_v_rec = .false. ! 退行速度阈在场旗标（文法可选段有无）
      logical :: outward  = .false.  ! 离面速度旗标（z_height 可选 OUTWARD：要求 p_z>0）
      integer :: n_max    = 0       ! 步数上限 [步]（steps 用）
      type(expr_t), allocatable :: sub(:)  ! 子表达式（and/or=2 个；not=1 个；原式无）
   end type expr_t
contains
   !------------------------------------------------------------------
   ! pred_parse(str, expr, ios, errmsg) — 解析 TERMINATE 表达式串 → 表达式树
   ! 参考: docs/final_completion_plan.md §5 N2（文法）；TEST.f:106-146/:88-94（原式语义）
   subroutine pred_parse(str, expr, ios, errmsg)
      character(len=*), intent(in) :: str          ! 表达式串（TERMINATE= 的值段）
      type(expr_t), intent(out) :: expr            ! 解析所得表达式树（根节点）
      integer, intent(out), optional :: ios        ! 解析状态 [—]（0=成功；非 0=文法错——
                                                   ! §5 N2 错误通道；缺省时错即 loud-stop）
      character(len=*), intent(out), optional :: errmsg ! 错误描述（与 ios 配套）
      ! 伪代码：
      ! 1. 递归下降按文法展开：表达式 := 原式 | AND(e1,e2) | OR(e1,e2) | NOT(e)
      !    （括号配对与逗号分段；关键字 COM_DISTANCE/Z_HEIGHT/MAX_STEPS/BOUND/
      !    AND/OR/NOT/ATOM 大写匹配）
      ! 2. 原式参数装节点字段：片段引用（名册片段名 A/B/… → i_frag；ATOM(i) →
      !    i_atom）、th_lo/th_hi [Å]、可选退行速度阈 → v_rec 且 has_v_rec=.true.、
      !    OUTWARD → outward=.true.、步数 → n_max
      ! 3. 布尔组合子式递归挂 sub(:)（and/or 两元、not 一元）
      ! 4. 文法错（未知原式/参数个数不符/片段引用不合法）→ ios 置错位 +
      !    errmsg 描述（或 loud-stop）；成功 ios=0
      ! 5. 无 TERMINATE 关键字时的默认表达式（现行 TEST 通用判据同参版本——N2.2，
      !    保证 B2 前行为不变）不在本例程装配：由装配层（input/driver，C3）按
      !    既有参数（旧 RCRIT 类阈值）构造同参树
   end subroutine pred_parse

   !------------------------------------------------------------------
   ! pred_eval(sta, expr, i_step) — 在当前物理态求值表达式树 → 终止布尔
   ! 参考: src_VENUS/TEST.f:106-146（原式求值原型）；TEST.f:88-94（z_height 语义）
   logical function pred_eval(sta, expr, i_step)
      type(state_t), intent(in) :: sta          ! 物理态（只读：q/p 参与求值）
      type(expr_t), intent(in) :: expr          ! 表达式树（pred_parse 或装配层产物）
      integer, intent(in), optional :: i_step   ! 当前步序 [步]（steps 原语的安全网
                                                ! 通道，driver 演化环传入；缺省时
                                                ! MAX_STEPS 不触发）
      ! 函数值：表达式在当前态的真值 [.true.=触发终止/事件]
      ! 伪代码：
      ! 1. 空树（node_kind=0）→ 假（未装配不拦；默认表达式装配见 pred_parse 注 5）
      ! 2. 布尔组合节点递归求值子式后调组合例程：
      !    k_and → pred_and(eval(sub(1)), eval(sub(2)))；
      !    k_or  → pred_or(eval(sub(1)), eval(sub(2)))；
      !    k_not → pred_not(eval(sub(1)))（先左后右的递归序=短路求值序）
      ! 3. 原式节点按类别调四原式之一（节点字段即哑元）：
      !    k_com_dist → com_dist(sta, frag_a, frag_b, th_lo, v_rec〔has_v_rec 时传〕)
      !    k_z_height → z_height(sta, frag_a, th_lo, outward)
      !    k_bound    → bound(sta, frag_a, frag_b, th_lo, th_hi)
      !    k_steps    → steps(i_step, n_max)
   end function pred_eval

   !------------------------------------------------------------------
   ! com_dist(sta, ref_a, ref_b, r_th, v_rec) — 原式：两片段质心距离（可选退行速度条件）
   ! 参考: src_VENUS/TEST.f:106-146（质心累加 :107-130、质心距/退行点积 :133-143、
   !       双阈判 :144-145——旧三值语义的消解见注 5）
   logical function com_dist(sta, ref_a, ref_b, r_th, v_rec)
      type(state_t), intent(in) :: sta          ! 物理态（只读）
      type(frag_ref_t), intent(in) :: ref_a     ! 片段引用 A（名册片段或单原子）
      type(frag_ref_t), intent(in) :: ref_b     ! 片段引用 B
      real(8), intent(in) :: r_th               ! 距离阈 [Å]（旧 RBAR/RMAX 的单阈化）
      real(8), intent(in), optional :: v_rec    ! 退行速度阈 [Å/(10 fs)]（缺省不检查；
                                                 ! 旧 RMAX 终止判据以 v_rec=0 等价）
      ! 函数值：质心距 ≥ r_th 且（有 v_rec 时）径向退行速度 > v_rec
      ! 伪代码：
      ! 1. 片段引用展开为原子表：i_atom>0 → 单原子伪片段；否则名册片段表
      !    （roster%list_a；B 侧片段表随 roster 留白——C 阶段定形）
      ! 2. 两片段各算质量加权质心与质心速度：q_cm = Σ(w·q)/Σw、
      !    v_cm = Σp/Σw（逐原子累加，旧 :107-130；速度=动量聚合除总质量）
      ! 3. r_cm = |q_cmA − q_cmB|；v_r = (q_cmA−q_cmB)·(v_cmA−v_cmB)
      !    （点积不除距离，旧 :133-143）
      ! 4. 真值 = r_cm ≥ r_th；v_rec 在场时再要求 v_r > v_rec（AND 语义）
      ! 5. 旧双阈双事件语义映射：RBAR 中间事件（旧 :144）= 本原式无 v_rec 版；
      !    RMAX 终止+退行（旧 :145）= 带 v_rec=0 版——旧三值语义（NTEST 0/1/2）
      !    消解为布尔组合与结局分类（outcome_tbl，C7）
   end function com_dist

   !------------------------------------------------------------------
   ! z_height(sta, ref, z_th, outward) — 原式：片段质心离面高度（可选外向速度）
   ! 参考: src_VENUS/TEST.f:88-94（C 原子写死判据：rd=Q(3)、P(3)>0、NC>50 预热——
   !       本原式即其通用化：写死原子 1 → 片段引用）
   logical function z_height(sta, ref, z_th, outward)
      type(state_t), intent(in) :: sta           ! 物理态（只读）
      type(frag_ref_t), intent(in) :: ref        ! 片段引用（表面体系=吸附物种片段）
      real(8), intent(in) :: z_th                ! 高度阈 [Å]（基准=表面顶层层合 z=0）
      logical, intent(in), optional :: outward   ! 真=另要求离面速度 p_z>0（旧 P(3)>0）
      ! 函数值：片段高度 z ≥ z_th 且（outward 时）片段 p_z > 0
      ! 伪代码：
      ! 1. 片段引用展开（同 com_dist 步 1）；i_atom>0 直取该原子
      !    z=q(3i)、p_z=p(3i)（旧 :89 即单原子形态 rd=Q(3)）
      ! 2. 片段质心高度 z_cm = Σ(w·q_z)/Σw、质心动量 p_z = Σp_z
      !    （质量加权，口径同 com_dist 步 2）
      ! 3. 真值 = z_cm ≥ z_th；outward 在场时再要求 p_z > 0（旧 :91-92 的 AND 语义）
      ! 4. 旧 :90 的 NC>50 预热跳过不迁——防初位误触发保护由表达式组合（如
      !    AND(NOT(…), …)）或 driver 初帧免检承担，C3 装配时定形
   end function z_height

   !------------------------------------------------------------------
   ! bound(sta, ref_a, ref_b, r_lo, r_hi) — 原式：有界窗越窗（约束散射）
   ! 参考: src_VENUS/TEST.f:133-143（距离计算复用 com_dist 步 1-3）；文法 N2 规格
   logical function bound(sta, ref_a, ref_b, r_lo, r_hi)
      type(state_t), intent(in) :: sta       ! 物理态（只读）
      type(frag_ref_t), intent(in) :: ref_a  ! 片段引用 A
      type(frag_ref_t), intent(in) :: ref_b  ! 片段引用 B
      real(8), intent(in) :: r_lo            ! 窗下阈 [Å]（过近界——复合/进入反应区）
      real(8), intent(in) :: r_hi            ! 窗上阈 [Å]（过远界——散射离场）
      ! 函数值：真 = 质心距越出 [r_lo, r_hi]（离开有界区即触发；假 = 窗内续演化）
      ! 伪代码：
      ! 1. 复用 com_dist 的片段展开与质心距离步骤（步 1-3），得 r_cm
      ! 2. 真值 = (r_cm ≤ r_lo) 或 (r_cm ≥ r_hi)
      ! 3. 约束散射用法：TERMINATE=BOUND(A,B,r_lo,r_hi) 或与 COM_DISTANCE 组合
      !    实现约束窗采样；越窗的事件类别细分交结局分类表（C7）
   end function bound

   !------------------------------------------------------------------
   ! steps(i_step, n_max) — 原式：步数安全网（恒为兜底，取代旧步数魔数）
   ! 参考: docs/final_completion_plan.md §5 N2（MAX_STEPS——恒为安全网）
   logical function steps(i_step, n_max)
      integer, intent(in) :: i_step   ! 当前步序 [步]（driver 演化环传入——pred_eval 通道）
      integer, intent(in) :: n_max    ! 步数上限 [步]
      ! 函数值：真 = i_step ≥ n_max（达上限即终止——防失止轨迹的恒在安全网）
      ! 伪代码：
      ! 1. 真值 = i_step ≥ n_max
      ! 2. 默认装配：无 TERMINATE 关键字时默认表达式含本原式（N2.2 同参版本）
   end function steps

   !------------------------------------------------------------------
   ! pred_and(l_a, l_b) — 布尔组合 AND(e1,e2)：两子式皆真
   ! 参考: docs/final_completion_plan.md §5 N2（文法：AND(e1,e2)）
   logical function pred_and(l_a, l_b)
      logical, intent(in) :: l_a   ! 子式 e1 真值
      logical, intent(in) :: l_b   ! 子式 e2 真值
      ! 函数值 = l_a .and. l_b
      ! 伪代码：
      ! 1. 两子式真值合取（求值序/短路由 pred_eval 递归序承担）
   end function pred_and

   !------------------------------------------------------------------
   ! pred_or(l_a, l_b) — 布尔组合 OR(e1,e2)：任一子式真
   ! 参考: docs/final_completion_plan.md §5 N2（文法：OR(e1,e2)）
   logical function pred_or(l_a, l_b)
      logical, intent(in) :: l_a   ! 子式 e1 真值
      logical, intent(in) :: l_b   ! 子式 e2 真值
      ! 函数值 = l_a .or. l_b
      ! 伪代码：
      ! 1. 两子式真值析取
   end function pred_or

   !------------------------------------------------------------------
   ! pred_not(l_a) — 布尔组合 NOT(e)：子式取反
   ! 参考: docs/final_completion_plan.md §5 N2（文法：NOT(e)）
   logical function pred_not(l_a)
      logical, intent(in) :: l_a   ! 子式 e 真值
      ! 函数值 = .not. l_a
      ! 伪代码：
      ! 1. 子式真值取反
   end function pred_not
end module predicate
