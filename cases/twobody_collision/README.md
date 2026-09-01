# twobody_collision — 双反应物碰撞组装（MB / MB）

> 运行：`cd cases/twobody_collision && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

验证双反应物组装（`INIT_SAMPLING_A/B=MB`）在 H₂(A)+H(B) / LEPS 势上的正确性
（NT=200）：碰撞参数 b² 均匀、Euler 取向球面均匀、固定相对平动能 Erel 精确闭合、
系统质心静止与刚体键长固定。

## 2. 理论与体系

体系为 TEST 势 LEPS 三原子：A 刚性 H₂（键长 0.7414 Å），B 单原子 H，全为 1.008 u。
流程：A 双原子 MB（TRVA=0 → 内部动量恒零、冷刚性转子；位置取刚性平衡键，
`ROTATE` 随机取向）→ `setup_gas_collision`——碰撞参数 SB=BMAX·√U（故 b² 在
[0,BMAX²] 上均匀）；B 平移到 (0, SB, √(S²−SB²))；相对速度 √(2·SEREL/μ) 按质量比
分配 VELA/VELB，总质心严格静止。Euler 取向：RPHI=2π·R1、cosθ=2R2−1、RCHI=2π·R3
三次独立均匀抽取，ZYZ 旋转矩阵组装——球面均匀。全链跑 RAND0(ISEED)，双跑逐位确定。

## 3. 方法与流程

1. 读取输入并逐轨迹打印 IMPACT PARAMETER（stdout，results.txt；第 1 条 0.482 Å）。
2. 初始化随机数（ISEED=20260832），NT=200 条轨迹抽取 SB=BMAX·√U 与 Euler 取向；
   `setup_gas_collision` 组装 A/B 坐标与按质量比分配的速度。
3. 每轨迹传播 NS=5 步并打印相空间（fort.1001-1200 每文件一条轨迹）。

## 4. 核心验证

见 `fig_long_trajectory.png`。

> 图注：左 = 原子 1 到最近原子的距离 r_min(t)（NS=5 短轨迹内即 A 片段键伴距离，
> 0.74-0.78 近恒定——碰撞远未发生），右 = 原子 1 动量模 |P|(t)（fort.1001-1010，
> 10 条轨迹）。判定依据：SB²/BMAX² KS D=0.0345（p=0.96，200 行 IMPACT PARAMETER，
> SB ∈ [0.241, 4.992] Å）；cos(θ) KS D=0.0512（p=0.65）；max|Erel−5.0|=0.0001
> kcal/mol；键长偏差 ≤1.1×10⁻⁵ Å，总动量 |ΣP|≤1.0×10⁻⁵（质心静止）。

> 散射判据（已修复）：此前 140/200 轨迹跑满 NS、坐标漂移至 1e6 Å，是 TEST.f 气相
> 终止判据缺失所致；现为通用质心距离判据 RCM_AB≥RMAX=8 Å + 外向质心速度
> （NTEST=2）。重跑验证：NON-REACTIVE 200/200，TOO LONG 归零。图亦为散射终止
> NT=20/NS=300000 重跑重画：r_AB(t) 在 ~8 Å（RMAX）处截止——轨迹在此判定为已
> 散射而终止，不再无限漂移。
