# 微正则简正模采样 — 固定内能的模能量分配

> 运行：`cd cases/microcanonical_normalmodes && ../../venus_test.e`
>（ISEED 固定，双跑逐位一致）。

## 1. 目的

本 case 验证微正则简正模采样（INIT_SAMPLING_A=MICROCANONICAL，即
NACTA=2）在 LEPS 等边 H₃ 体系上的统计正确性：简正模频谱锚点、
模间能量分配、以及能量闭合至固定内能 ENMT_A；同时锁定祖传缺陷
F19（振动角动量剥离）的确定性签名。判据来源为微正则系综理论：
固定总内能下各简正模能量份额应服从均匀 Dirichlet(1,1,1) 分布。

## 2. 理论与体系

体系为等边 H₃（3×1.008 u）气相 TEST 势 LEPS。势参数（源码
src_TEST/test_potentials.f90:44-48；input 注释的 [kcal/mol][Å][Å⁻¹]
单位有误，源码为 eV/bohr 制）：

| 参数 | 值 | 单位 | 出处 |
|---|---|---|---|
| De | 4.746 | eV | LEPS_DE，test_potentials.f90:44 段 |
| re | 1.401 | bohr | LEPS_RE |
| a | 1.028 | bohr⁻¹ | LEPS_A |
| Sato Δ | 0.164 | 无量纲 | LEPS_DELTA，test_potentials.f90:48 |
| ENMT_A | 1.0 | kcal/mol | input_qct.txt |
| NT | 1000 | 条 | input_qct.txt |

简正模频谱（fort.26 频谱块实测）：三个平动模 <1 cm⁻¹，E′ 二重
简并对 1013.5 cm⁻¹（切向），A₁′ 呼吸模 1803.8 cm⁻¹（径向）。
机制链：`select_polyatomic_a`→`NMODE`（频谱与简正向量）→序贯能量
抽样（各模份额 p_i=rest·(1−u^(1/(NMBAR−i)))，数学上等价均匀
Dirichlet(1,1,1)）→`INITQP`（相位组装 + 振动角动量剥离；能量闭合
缩放 P 与 (Q−QZ) 同乘 √(HSCALE/H)，直至 |H−HSCALE|/HSCALE≤10⁻³，
NSCALE≤50）。换算常数 C1=0.04184（kcal/mol→code，
venus_params.f90:69）。

## 3. 方法与流程

1. 读入 QZA_EQ 等边构型，NMODE 做简正模分析，打印频谱与向量
   （fort.26 频谱块 + stdout 频谱表，见 results.txt）。
2. 序贯抽样把 ENMT_A 分配到各振动模，INITQP 组装相空间并做
   振动角动量剥离（F19，见 §4）。
3. 能量闭合缩放至完整非谐 LEPS 势，stdout 逐轨迹打印
   `INTERNAL ENERGY` 迭代收敛（results.txt 示例回显）。
4. CENMAS 去质心动量后入库（fort.8/fort.1001），fort.999 汇总。

## 4. 核心验证：能量闭合

判据：每条轨迹末态内能相对 ENMT_A 偏差 ≤10⁻³（缩放循环验收门）。
实测（NT=1000 stdout 末次 `INTERNAL ENERGY`，results.txt）：
1000/1000 轨迹达标（最大值恰命中门限），初始内能散布
0.99983-1.00054 kcal/mol 后收敛于 1.0。结论：微正则能量闭合通过。
佐证：频谱锚点 1013.47/1013.51/1803.83 cm⁻¹（fort.26，容差
±1 cm⁻¹）；F19 签名——呼吸模能量份额均值 0.2008 ∈ [0.17,0.24]
（无偏应为 0.333），E′ 对合计 0.7992，逐位可复现（祖传行为：
INITQP.f:131-137 的角动量剥离非线性，登记不修，凡依赖此采样的
统计结论均须携带该偏置）；质心动量去除 max|Σp|=1.0×10⁻⁵。

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量
>（stdout `INTERNAL ENERGY` 回显）：左 = 各轨迹初始内能
>（0.99983-1.00054 kcal/mol，红线 = 目标 ENMT_A=1.0），右 = 逐轨迹
> 缩放收敛过程（首猜偏离 → √(HSCALE/H) 同乘收敛）。
