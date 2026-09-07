!=====================================================================
! dist_ebk.f90 — L3 分布成员示例（B′ 注册族）：EBK 半经典态分布——固定
!                 (n,J) 量子数的双原子初态成员（B′ 第二实现原型；示例 2/3）
! 参考: src_VENUS/INITEBK.f90:1-38（能级初值与频率步长 :15-16、转动能矩
!       AL=√(J(J+1))·C7 :17-18、折合质量 :19、固定点迭代 :24-33、转折点
!       内缩 :35-36、动量尺度 PTEST :37）、src_VENUS/FINLNJ.f:1-104（作用量
!       积分的 RKR 反演——有效势 VEFF=V−DH+离心项 :48、转折点扫描 :56-81、
!       Gauss-Legendre 求积 :85-95）、src_VENUS/SELECT.f90:337-351（固定量子
!       数分支 TRV_A<0 与能级初值 ENJA=(n+½)·ν·cm2cal·C1 :347）、
!       src_VENUS/HOMOQP.f90（r/pr → 双原子坐标动量装配——SELECT.f90:373 调用）、
!       cases/ebk_fixed_nj/README.md:18-25（能级公式 hc[ωe(n+½)−ωexe(n+½)²]
!       +hc·Be·J(J+1) 与 Morse 谱常数对照）、cases/ebk_fixed_nj/input_qct.txt:36-39
!       （固定量子数/转动惯量输入位）
! 注册：本成员进 sampler 分布成员表（samp_reg 一行——装配相位执行；C5 接线）
!=====================================================================
module dist_ebk
   implicit none
   private
   public :: ebk_params_t, ebk_init

   ! 成员自带参数（E2；按旧 INITEBK 及其调用链实际读到的输入列举，留白允许
   ! ——C 阶段定形。折合质量不存——由 roster 质量表算（旧 W(LBOND)，
   ! INITEBK.f90:19）；分离能不存——旧 SELECT.f90:121 DH=V 由门面求值，
   ! 能量归 facade）
   type :: ebk_params_t
      integer :: n_vib  = 0        ! 固定振动量子数 n [—]（旧 NNA——TRV_A<0
                                   ! 规定，input_qct.txt:37）
      integer :: j_rot  = 0        ! 固定转动量子数 J [—]（旧 JA——同上 :38）
      real(8) :: w_e    = 0.0d0    ! 谱常数 ωe [cm⁻¹]（Morse 谐频项——能级公式
                                   ! (n+½) 一次项；对照值 2326.4）
      real(8) :: w_ex_e = 0.0d0    ! 非谐谱常数 ωexe [cm⁻¹]（能级公式 (n+½)² 项；
                                   ! 对照值 35.35——由势参数导出，C5 定形）
      real(8) :: b_rot  = 0.0d0    ! 转动常数 Be [cm⁻¹]（刚性转子项 J(J+1)）
      real(8) :: ai_rot = 0.0d0    ! 转动惯量 [amu·Å²]（旧 AIA——装配/取向用；
                                   ! 对照值 0.98925，input_qct.txt:39）
      ! 留白：转折点初值 r_min/r_max [Å]、固定点迭代容差（旧 1e-6）与上限
      ! （旧 1000 次，INITEBK.f90:31-32）、求积节点数（旧 50）——C5 定形
   end type
contains
   !------------------------------------------------------------------
   ! ebk_init(p) — 初始化成员（校验量子数/谱常数、预计算能级与转动能矩）
   ! 参考: src_VENUS/INITEBK.f90:15-19（HNU 步长、AL=√(J(J+1))·C7、折合质量）
   subroutine ebk_init(p)
      type(ebk_params_t), intent(in) :: p ! 成员参数（量子数/谱常数）
      ! 伪代码：
      ! 1. 校验参数合法性（非法 loud-stop）：n_vib ≥ 0、j_rot ≥ 0、w_e > 0、
      !    w_ex_e ≥ 0、b_rot ≥ 0、ai_rot > 0
      ! 2. 预计算态能级初值 [kcal/mol]：
      !    e_nj = (w_e·(n+½) − w_ex_e·(n+½)² + b_rot·J(J+1))·wvn2e
      !    （cm⁻¹→kcal/mol 经 consts；旧初值 (n+½)·ν·cm2cal·C1，
      !    SELECT.f90:347——能级公式对照 cases/ebk_fixed_nj/README.md:20）
      ! 3. 预计算转动能矩 al = √(J(J+1))·hbar_code [ħ 单位]（旧 AL=AM·C7，
      !    INITEBK.f90:17-18）与频率步长 hnu = e_nj/(n+½)（旧 HNU :16）
   end subroutine ebk_init
   ! —— sample 方向（C5 定形）——固定点迭代 e_nj ← e_nj + (n−an)·hnu 至
   !   |n−an| ≤ 1e-6（an 为 FINLNJ 型作用量积分——WKB 密度 p(r)∝1/|p_r|，
   !   转折点 [r_min,r_max] 内拒绝抽样核间距 r、径向动量 pr=√(2μ·summ) 符号
   !   对半翻转，再经 HOMOQP 型装配写初态分量）
   ! sample 例程：C5 成员具体化时定形（本阶段不声明）
   ! 注册行注记：samp_reg 登记本成员（名 'ebk'）进 sampler 分布成员表——C5 接线
end module dist_ebk
