# 微正则简正模采样 — 固定内能的模能量分配

> 运行：`cd cases/microcanonical_normalmodes && ../../venus_test.e`
>（ISEED 固定，双跑逐位一致）。

## 1. 目的

验证微正则简正模采样（INIT_SAMPLING_A=MICROCANONICAL，即 NACTA=2）在 LEPS 等边
H₃（3×1.008 u）体系上的统计正确性：模间能量分配闭合至固定内能 ENMT_A=1.0
kcal/mol（NT=1000）；同时锁定祖传缺陷 F19（振动角动量剥离）的确定性签名。

## 2. 理论与体系

TEST 势 LEPS（test_potentials.f90:44-48；Sato Δ=0.164）。简正模频谱（fort.26
实测）：三个平动模 <1 cm⁻¹，E′ 二重简并对 1013.5 cm⁻¹（切向），A₁′ 呼吸模
1803.8 cm⁻¹（径向）。流程：`select_polyatomic_a`→`NMODE`（频谱与简正向量）→
序贯能量抽样（各模份额 p_i=rest·(1−u^(1/(NMBAR−i)))，数学上等价均匀
Dirichlet(1,1,1)）→`INITQP`（相位组装 + 振动角动量剥离；能量闭合缩放 P 与
(Q−QZ) 同乘 √(HSCALE/H)，直至 |H−HSCALE|/HSCALE≤10⁻³，NSCALE≤50）。

## 3. 方法与流程

1. 读入 QZA_EQ 等边构型，NMODE 做简正模分析并打印频谱（fort.26 + stdout 频谱表）。
2. 序贯抽样把 ENMT_A 分配到各振动模，INITQP 组装相空间并做振动角动量剥离（F19）。
3. 能量闭合缩放至完整非谐 LEPS 势（stdout 逐轨迹 `INTERNAL ENERGY` 迭代收敛）；
   CENMAS 去质心动量后输出（fort.8/fort.1001），fort.999 汇总。

## 4. 核心验证

见 `fig_sampling_stats.png`。

> 图注：左 = 各轨迹初始内能（0.99983-1.00054 kcal/mol，红线 = 目标 ENMT_A=1.0），
> 右 = 逐轨迹缩放收敛过程（首猜偏离 → √(HSCALE/H) 同乘收敛；stdout `INTERNAL
> ENERGY`，10 条轨迹）。判定依据：1000/1000 轨迹末态内能相对 ENMT_A 偏差 ≤10⁻³
>（缩放循环验收门，最大值恰命中门限）；频谱对照 1013.47/1013.51/1803.83 cm⁻¹
>（fort.26，容差 ±1 cm⁻¹）。F19 签名：呼吸模能量份额均值 0.2008 ∈ [0.17,0.24]
>（无偏应为 0.333），逐位可复现（INITQP.f:131-137 角动量剥离非线性，登记不修）。
