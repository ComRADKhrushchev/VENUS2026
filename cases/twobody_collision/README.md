# twobody_collision — 双反应物碰撞组装（MB / MB）

> 运行：`cd cases/twobody_collision && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

本 case 验证双反应物组装（`INIT_SAMPLING_A/B=MB`）在 H₂(A)+H(B) / LEPS 势上的正确性
（NT=200）：碰撞参数 b² 均匀、Euler 取向球面均匀、固定相对平动能 Erel 精确闭合、
系统质心静止与刚体键长固定。

## 2. 理论与体系

体系为 TEST 势 LEPS 三原子：A 刚性 H₂（键长 0.7414 Å），B 单原子 H，全为 1.008 u。
LEPS 参数取自 input_qct.txt（De=4.746 eV、re=1.401 bohr、a=1.028 bohr⁻¹、Δ=0.164；
源码 `src_TEST/test_potentials.f90:44-47` 定死单位为 eV/bohr/bohr⁻¹，input 注释
"[Å][Å⁻¹]"有误，数值一致仅单位标注错误）。机制链：A 双原子 MB（TRVA=0 → 内部动量恒零、冷刚性转子；位置取刚性
平衡键；`ROTATE` 随机取向）→ `setup_gas_collision`——碰撞参数 SB=BMAX·√U（故 b² 在
[0,BMAX²] 上均匀）；B 平移到 (0, SB, √(S²−SB²))；`NREL=1` 时 SEREL=EREL·C1（C1=0.04184
为 kcal/mol→a.u. 换算，`venus_params.f90:60-80`），相对速度 √(2·SEREL/μ) 按质量比分配
VELA/VELB，总质心严格静止。Euler 取向：RPHI=2π·R1、cosθ=2R2−1、RCHI=2π·R3 三次独立
均匀抽取，ZYZ 旋转矩阵组装——球面均匀。全链跑 RAND0(ISEED)，双跑逐位确定。

## 3. 方法与流程

1. 程序读取输入并逐轨迹回显 IMPACT PARAMETER（stdout，results.txt；第 1 条 0.482 Å）。
2. 初始化随机数（ISEED=20260832），对 NT=200 条轨迹抽取 SB=BMAX·√U 与 Euler 取向。
3. `setup_gas_collision` 组装 A/B 坐标与按质量比分配的速度（stdout，fort.1001 起可见
   A 内部动量恒零、B 位于 (0, SB, √(S²−SB²))）。
4. 每轨迹传播 NS=5 步并打印相空间（fort.1001-1200 每文件一条轨迹）。

## 4. 核心验证

**判据**：b² 应在 [0,BMAX²] 均匀、cos(θ) 取向应均匀（KS 检验），固定相对能应精确闭合。
**实测**（results.txt，200 行 IMPACT PARAMETER，SB ∈ [0.241, 4.992] Å）：
SB²/BMAX² KS D=0.0345（p=0.96）；cos(θ) KS D=0.0512（p=0.65）；max|Erel−5.0|=
0.0001 kcal/mol（标准差 8.9×10⁻¹⁶）。**结论**：碰撞参数与取向采样机制按设计工作。
佐证：键长偏差 ≤1.1×10⁻⁵ Å，总动量 |ΣP|≤1.0×10⁻⁵（质心静止）。TRVA=0 时 H₂ 为
冷刚性转子（无振/转能）；温度模式 Erel 采样（NREL=0）为统计性分支，不在本 case 范围。

> 图：`fig_trajectory.png` —— 10 条入库轨迹的相空间（fort.1001-1010 活数据）：
> 左为原子 1 到最近原子的距离 r_min(t)（NS=5 短轨迹内即 A 片段键伴距离，0.74-0.78 bohr
> 近恒定——碰撞远未发生），右为原子 1 动量模 |P|(t)。
