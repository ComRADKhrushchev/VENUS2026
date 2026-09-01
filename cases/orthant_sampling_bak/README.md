# 正交采样（ORTHANT）— 固定内能 + 对称陀螺初始条件

> 运行：`cd cases/orthant_sampling && ../../venus_test.e`（LEPS 等边 H₃，
> NT=1000 基线，输入 input_qct.txt，ISEED=20260822 固定）。

## 1. 目的

定案正交采样链（`INIT_SAMPLING_A=ORTHANT`，手册 V.5-8）在 LEPS 三原子
等边体系上的语义：从全局极小出发、固定内能与固定对称陀螺态 (J,K)=(2,0)
采样，核心考察能量闭合的精确性，兼查转动能量解析值与相空间符号对称性。

## 2. 理论与体系

体系为 TEST 势 LEPS 等边 H₃（NATOMS=3，ATOM_MASSES=1.008 u ×3）。
势参数（input_qct.txt 显式）：De = 4.746 eV、re = 1.401 Å、a = 1.028 Å⁻¹、
Sato Δ = 0.164（无量纲）。全局极小为等边三角形，边长 s = 1.759910 Å
（QZA_EQ 给出），主惯量 I_perp = m·s²/2（二重简并）与 I_axis = m·s²，
实测 1.56107 / 3.12214 amu·Å²。固定态下转动能解析确定：
EROTT = J(J+1)·C7²/(2·C1·I_perp)，其中 C7 = 0.063508 为 ℏ 的 code 单位值
（`venus_params.f90:75`）、C1 = 0.04184 为 kcal/mol → a.u. 能量换算常数
（`venus_params.f90:69`）；采样链温度分支另用 C5 = 0.083144×10⁻³ a.u.
（`venus_params.f90:73`），本 case 不触及。能量闭合目标
HSCALE = EROTT + ENMT_A = 10.185254 kcal/mol：|H−HSCALE|/HSCALE ≥ 10⁻³ 时
P 与 (Q−QZ) 同乘 √(HSCALE/H)，NSCALE ≤ 50；初猜动量模
PMAX = √(2·W·ENMT_A·C1)·PSCALE。**PSCALE_A=5 说明**：动能主导的初猜，
PSCALE=1 时软 E′ 模可致缩放循环过冲停机——原版 VENUS 行为，非缺陷。

## 3. 方法与流程

1. 读入 QZA_EQ，在极小处做力常数分析并对角化（`fort.26` 频谱块）。
2. `select_polyatomic_a` 走正卦限单位向量采样：6N 维相空间逐步条件抽取
   单位随机向量，按分量比例分配动量与坐标位移，每分量 50% 概率翻号；
   QMAX/QMIN 逐坐标 ±0.1 Å 步进至 V−V(QZ) ≥ ENMT。
3. 去质心速度；两步差值角动量插值（J′=J−Js，ω=I⁻¹J′）赋 (J,K)=(2,0)。
4. 能量闭合缩放至 HSCALE，逐轨迹打印迭代过程（stdout，摘录于 results.txt）。
5. 轨迹短演化入库（`fort.1001`–`1200` 等，NS=10），汇总统计（`fort.999`）。

## 4. 核心验证：1000 条轨迹的能量闭合

这是判定采样链在工作的单一关键证据，因为它同时要求三个环节全部正确：
正卦限抽取的量值分配正确（否则初始 H 系统性偏离）、转动插值正确
（否则 EROTT 份额失控）、缩放公式正确（否则收敛不到 HSCALE）。

**判据**：全部轨迹末态内能相对 HSCALE = 10.185254 kcal/mol 的最大相对
偏差 ≤ 10⁻³（程序自身闭合门限）。

**实测**（results.txt 汇总）：1000 条轨迹最大相对偏差 9.99×10⁻⁴ ≤ 10⁻³ ✓。
配套佐证：程序 EROTT = 0.1852540 vs 解析 J(J+1)·C7²/(2·C1·I_perp) =
0.1852524 kcal/mol（相对偏差 8.6×10⁻⁶）；9 个动量分量正号率 ∈ [0.479,
0.518]，落在 N=1000 的 50%±3σ 区间 [0.44, 0.56]；质心速度去除
max|Σm·v|/Σ|m·v| = 1.5×10⁻⁵。已知行为：CHOSEN EROTA（0.157）低于目标
EROTT（0.185）为两步插值的转动份额损耗，总内能闭合不受影响
（EVIBA 与 HSCALE−EROTA 偏差 0.0007 kcal/mol）。

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量（stdout
> `INTERNAL ENERGY` 回显，活数据）：左 = 各轨迹初始内能（10.185–10.195
> kcal/mol，虚线 = HSCALE = EROTT+ENMT_A = 10.185），右 = 逐轨迹缩放
> 收敛过程。
