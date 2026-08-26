# 双反应物碰撞（MB / MB）— 碰撞参数、Euler 取向与固定相对能

> 运行：`cd cases/twobody_collision && ../../venus_test.e`

## 验证目的

验证双反应物组装（`INIT_SAMPLING_A/B=MB`）在 H₂(A)+H(B) / LEPS 势上的正确性
（NT=200）：碰撞参数 b² 均匀、Euler 取向球面均匀、固定相对平动能 Erel 精确闭合、
系统质心静止、刚体键长固定。

## 原理

体系：TEST 势 LEPS；A 为刚性 H₂（键长 0.7414 Å），B 为单原子 H。

机制链：A 双原子 MB（TRVA=0 → 内部动量恒零；位置取刚性平衡键；`ROTATE` 随机
取向）→ `setup_gas_collision`——碰撞参数 SB=BMAX·√U（故 b² 在 [0,BMAX²] 上
均匀）；B 平移到 (0, SB, √(S²−SB²))；`NREL=1` 时 SEREL=EREL·C1，相对速度
√(2·SEREL/μ) 按质量比分配 VELA/VELB（总质心严格静止）。

Euler 取向：RPHI=2π·R1、cosθ=2R2−1、RCHI=2π·R3（三次独立均匀抽取），ZYZ 旋转
矩阵组装——球面均匀。全链跑 RAND0(ISEED)（GASDEV 输出乘 DESKET=0 不触系统
种子）→ 双跑逐位确定。

## 预期与结果

- **Erel 精确闭合**：max|Erel−5.0|=0.0001 kcal/mol（标准差 8.9×10⁻¹⁶）✓
- **b² 与取向均匀**：SB²/BMAX² KS D=0.0345（p=0.96）；cos(θ) KS D=0.0512
  （p=0.65）✓
- **几何与动量约束**：键长偏差 ≤1.1×10⁻⁵ Å；总动量 |ΣP|≤1.0×10⁻⁵ ✓

补注：TRVA=0 时 H₂ 为冷刚性转子（无振/转能）；温度模式 Erel 采样（NREL=0）
为统计性分支，不在本 case 范围。

> 图：`fig_trajectory.png` —— 10 条入库轨迹的相空间（fort.1001-1010 活数据）：
> 左=原子 1 到最近原子的距离 r_min(t)（NS=5 短轨迹内即 A 片段键伴距离，
> 0.74-0.78 bohr 近恒定——碰撞远未发生），右=原子 1 动量模 |P|(t)。GWRITE 能量
> 打印为陈旧值（见 morse_bootstrap 登记），本图只用 fort 相空间。
