# 正交采样（ORTHANT）— 固定内能 + 对称陀螺初始条件

> 运行：`cd cases/orthant_sampling && ../../venus_test.e`（LEPS 等边 H₃，NT=1000
> 基线，输入 input_qct.txt，ISEED=20260822 固定）。

## 1. 目的

验证正交采样链（`INIT_SAMPLING_A=ORTHANT`，手册 V.5-8）在 LEPS 等边 H₃ 上的语义：
从全局极小出发、固定内能与固定对称陀螺态 (J,K)=(2,0) 采样，核心考察能量闭合的
精确性，兼查转动能量解析值与相空间符号对称性。

## 2. 理论与体系

体系为 TEST 势 LEPS 等边 H₃（3×1.008 u；De=4.746 eV、re=1.401 Å、a=1.028 Å⁻¹、
Sato Δ=0.164）。全局极小为等边三角形，边长 s=1.759910 Å（QZA_EQ），主惯量
I_perp=m·s²/2（二重简并）与 I_axis=m·s²，实测 1.56107/3.12214 amu·Å²。固定态下
转动能解析确定：EROTT = J(J+1)·ℏ²/(2·I_perp)。能量闭合目标 HSCALE =
EROTT+ENMT_A = 10.185254 kcal/mol（ENMT_A=10.0 输入 + EROTT 解析值 0.185254）：|H−HSCALE|/HSCALE ≥ 10⁻³ 时 P 与 (Q−QZ) 同乘
√(HSCALE/H)，NSCALE ≤ 50。**PSCALE_A=5 说明**：动能主导的初猜，PSCALE=1 时软
E′ 模可致缩放循环过冲停机——原版 VENUS 行为，非缺陷。

## 3. 方法与流程

1. 读入 QZA_EQ，在极小处做力常数分析并对角化（`fort.26` 频谱块）。
2. `select_polyatomic_a` 走正卦限单位向量采样：6N 维相空间逐步条件抽取单位随机
   向量，按分量比例分配动量与坐标位移，每分量 50% 概率翻号；QMAX/QMIN 逐坐标
   ±0.1 Å 步进至 V−V(QZ) ≥ ENMT。
3. 去质心速度；两步差值角动量插值（J′=J−Js，ω=I⁻¹J′）赋 (J,K)=(2,0)；能量闭合
   缩放至 HSCALE（stdout 打印迭代过程）；轨迹短演化输出（`fort.1001`–`1200`，
   NS=10），汇总统计（`fort.999`）。

## 4. 核心验证

见 `fig_sampling_stats.png`。

> 图注：左 = 各轨迹初始内能（10.185–10.195 kcal/mol，虚线 = HSCALE =
> EROTT+ENMT_A = 10.185），右 = 逐轨迹缩放收敛过程（stdout `INTERNAL ENERGY`，
> 10 条轨迹）。判定依据：1000 条轨迹末态内能相对 HSCALE 最大相对偏差 9.99×10⁻⁴
> ≤ 10⁻³（程序闭合门限）；EROTT = 0.1852540 vs 理论对照值 0.1852524 kcal/mol
>（相对偏差 8.6×10⁻⁶）；9 个动量分量正号率 ∈ [0.479, 0.518] ⊂ 50%±3σ [0.44, 0.56]。
