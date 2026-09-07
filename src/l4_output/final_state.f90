!=====================================================================
! analysis/final_state.f90 — L4 分析：终态量计算（广义 COM 相对能量/角动量/
!                            散射角族——名册片段表数据驱动，气相/表面统一）
! 参考: src_VENUS/FINAL.f:100-109（片段 A 角动量/转动能——CENMAS+ROTN，
!       /C7 折 ℏ）、:113-129（片段内动能 Σp²/2m/C1 与内势能 V−DELH——
!       DVDQ_1+ENERGY_1 于片段几何）、:140-144/:295-301（双原子 n/j 反演
!       FINLNJ）、:149-174（角速度扣除 angvel + QFINAL/PFINAL 终态存档）、
!       :356-372（全体系质心平动能 ETCM）、:376-385（折合质量两式：气相
!       WTA·WTB/(WTA+WTB)、表面 WTA）、:387-410（相对运动：质心差/速度差 →
!       EREL/ERELSQ/OAM/BF）、:432-531（气相散射角 ANG(1..16)——VI/VF、
!       LI/LF、JAI/JBI 两两点积 acos）、:538-586（表面角 ANG(17..18)——面法向
!       由 NN1-NN3 三原子叉积）、:598-637（表面入射/散射角 ANG(17..20)）；
!       src_VENUS/GFINAL.f:64-75（寿命 TLIFE 与能量守恒漂移打印——漂移归
!       recorder 能量列）、:199-202（AErelS/AThtaS 终态记录——归 recorder，
!       边界注）
! 剥离警示（EREL 写死覆盖——K2 计划取材注）：
!   GFINAL.f:78-81 以 C 原子动能 NORM2(P(7:9))²/2W(3)/C1 覆写 EREL（1D 模型
!   体系知识）——不迁：e_rel 一律按 FINAL.f:387-401 广义质心相对式；
!   同源不迁：FINAL.f:616-625 NZDOWN==1 分支写死 B 速度取 P(7:9)/W(3)——
!   表面终态速度通道经名册片段指定
!=====================================================================
module final_state
   use consts,        only: hbar_code, e_conv, dtor
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use geometry,      only: cenmas, amom, angvel
   use force_facade,  only: facade_force
   implicit none
   private
   public :: fin_t, fin_eval

   ! 终态量承接类型（fin_eval 出口；字段=旧全局终态量的成员化）
   type :: fin_t
      real(8) :: e_int_a(3) = 0.0d0   ! 片段 A 内能 [kcal/mol]：(1)内动能 (2)内势能−参考 (3)总（旧 EA(1..3)）
      real(8) :: e_int_b(3) = 0.0d0   ! 片段 B 内能 [kcal/mol]（旧 EB(1..3)；B 表 roster 留白——C 阶段）
      real(8) :: j_a(3)     = 0.0d0   ! 片段 A 角动量矢量 [ℏ]（旧 AMA(1:3)，/C7 折算）
      real(8) :: j_mag_a    = 0.0d0   ! 片段 A 角动量模 [ℏ]（旧 AMA(4)）
      real(8) :: e_rot_a    = 0.0d0   ! 片段 A 转动能 [kcal/mol]（旧 EROTA）
      real(8) :: j_b(3)     = 0.0d0   ! 片段 B 角动量矢量 [ℏ]（旧 AMB(1:3)）
      real(8) :: j_mag_b    = 0.0d0   ! 片段 B 角动量模 [ℏ]（旧 AMB(4)）
      real(8) :: e_rot_b    = 0.0d0   ! 片段 B 转动能 [kcal/mol]（旧 EROTB）
      real(8) :: e_rel      = 0.0d0   ! 相对平动能（径向投影）[kcal/mol]（旧 EREL——FINAL.f:394-401 广义式）
      real(8) :: e_rel_sq   = 0.0d0   ! 相对平动能（全速度）[kcal/mol]（旧 ERELSQ）
      real(8) :: oam(3)     = 0.0d0   ! 轨道角动量矢量 [ℏ]（旧 OAM(1:3)，/C7 折算）
      real(8) :: oam_mag    = 0.0d0   ! 轨道角动量模 [ℏ]（旧 OAM(4)）
      real(8) :: b_fin      = 0.0d0   ! 终态碰撞参数 [Å]（旧 BF = OAM(4)/μ/|v_rel|）
      real(8) :: e_cm       = 0.0d0   ! 全体系质心平动能 [kcal/mol]（旧 ETCM）
      real(8) :: ang(20)    = 0.0d0   ! 散射角族 [°]（旧 ANG(1..20)：1-3 LF·JAF/JBF·JAF/JAF·JBF，4-7 VI·VF/LI·LF/LI·JAF/LI·JBF，8-16 JAI/JBI 组合，17-18 表面入射 θ/χ，19-20 表面散射 θ/χ）
      real(8) :: t_life     = 0.0d0   ! 轨迹寿命 [fs]（旧 TLIFE = TIME·NC·10——GFINAL.f:64；新版 sta%t×10）
      real(8), allocatable :: q_fin(:) ! 终态质心系坐标存档 [Å]（旧 QFINAL——简正模分析用）
      real(8), allocatable :: p_fin(:) ! 终态质心系动量存档 [amu·Å/(10 fs)]（旧 PFINAL）
      ! 留白：双原子 n/j 反演（旧 AN/AJ/BN/BJ——FINLNJ，随 ebk 成员定形）；
      !       初值参考量（旧 VI/OAMI/AMAI/AMBI——sampler 装配时快照，传交通道
      !       C 阶段定形）；转动能抽样涨落 SDA/SDB（旧 NCROT 重复抽样平均——
      !       归统计相位）；片段 B 表（roster 留白同源）
   end type fin_t
contains
   !------------------------------------------------------------------
   ! fin_eval(sta, fin) — 终态量计算：片段内能/角动量 + 广义 COM 相对量 +
   !                      散射角族（名册片段表驱动，无体系分支）
   ! 参考: FINAL.f:100-174、:356-410、:432-531、:538-637（分段见模块头）
   subroutine fin_eval(sta, fin)
      type(state_t), intent(in) :: sta    ! 物理态（终态快照——只读；片段几何
                                          ! 暂存技巧见步 2，不改 sta 本体）
      type(fin_t), intent(out) :: fin     ! 终态量承接（本例程唯一出口）
      ! 伪代码：
      ! 1. 片段 A 属性（旧 :100-109）：cenmas(sta, roster%mass, roster%list_a,
      !    q_cm_a, v_cm_a) 建质心系；amom(…, j_a, e_rot_a, omega_a)——j_a 与
      !    j_mag_a 折 ℏ（除 hbar_code，旧 /C7）；单原子片段全零（旧 :24-33）
      ! 2. 片段 A 内能（旧 :113-129）：内动能 = 质心系 Σp²/(2·mass)/e_conv；
      !    内势能 = facade_force 于片段平衡排布坐标求 v（旧 :49-62 的 QZB+QCMB
      !    暂存求值技巧——新版经门面于局部暂存态求值，不写回 sta）减参考
      !    分离能（旧 V−DELH；DELH 能量基准归装配相位——C 阶段定形）；
      !    e_int_a(3) = (1)+(2)
      ! 3. 片段 B 同式（旧 :264-301；B 表 roster 留白——C 阶段接入）
      ! 4. e_cm（旧 :356-372）：全原子 Σp 三分量平方 / Σmass /2 /e_conv 求和
      ! 5. 相对运动（旧 :376-410）：折合质量 μ 依名册片段表——两片段气相式
      !    μ = m_a·m_b/(m_a+m_b)、表面式（单活性片段）μ = m_a（旧 :376-385
      !    两分支统一为片段表数据驱动，reactants%surface_model 只定通道）；
      !    q_r = q_cm_a − q_cm_b、v_r 同；e_rel = μ·(v_r·q̂_r)²/2/e_conv（径向
      !    投影，旧 :394-401）；e_rel_sq = μ·|v_r|²/2/e_conv；oam = μ·(q_r×v_r)
      !    折 ℏ（旧 :403-410）；b_fin = |oam|/μ/|v_r|（旧 :407）
      ! 6. 气相散射角 ang(1..16)（旧 :432-531）：初/终相对速度矢量（VI 与 v_r）
      !    与轨道/转动角动量矢量（OAMI/OAM、AMAI/AMBI 与 j_a/j_b）两两点积
      !    acos、除 dtor 转度（旧 /C4）；模 < 1e-5 的角动量跳项（旧 :460/:482
      !    门槛）；点积夹逼 [−1,1]（旧 :435-436）；初值量 VI/OAMI/AMAI/AMBI
      !    来自初态快照（留白通道）
      ! 7. 表面角 ang(17..20)（旧 :598-637）：入射角由初速 VI（θ/χ 两式
      !    :600-612）、散射角由活性片段质心速度 v_cm_a（:614-637；NZDOWN
      !    分支写死不迁——剥离警示见模块头）；极角退化位（|v_xy|≤1e-6 →
      !    χ=0，旧 :604-611/:629-636）；面法向基准（旧 :543-555 由参考三原子
      !    叉积——数据化经名册/晶胞字段，C 阶段定形）
      ! 8. 终态存档（旧 :149-174）：angvel 扣除刚体转动后的质心系坐标/动量
      !    → q_fin/p_fin（供简正模分析——density/成员消费）
      ! 9. t_life = sta%t·10 [fs]（内部时间 10 fs → fs，旧 TLIFE 口径）
      ! 10. 边界：GFINAL 的记录职责（AErelS/AThtaS :199-202）归 recorder；
      !     EREL 写死覆盖（GFINAL.f:78-81）不迁——e_rel 保持步 5 广义式
   end subroutine fin_eval
end module final_state
