# 正交采样（ORTHANT）— 固定内能 + 对称陀螺初始条件

> 运行：`cd cases/orthant_sampling && ../../venus_test.e`

## 验证目的

验证正交采样链（`INIT_SAMPLING_A=ORTHANT`）在 LEPS 三原子等边 H₃ 类体系上的统计
正确性（NT=1000）：从 LEPS 全局极小出发，以固定内能 `ENMT_A=10 kcal/mol` 与固定
对称陀螺态 (J,K)=(2,0) 采样，考察能量闭合、转动能量解析值与相空间符号对称性。

## 原理

体系：TEST 势 LEPS（Sato Δ=0.164，H-H 对 De=4.746 eV、re=1.401 Å、a=1.028 Å⁻¹）。
全局极小为等边三角形（边长 s=1.7599 Å），主惯量 I_perp=m·s²/2（二重简并）、
I_axis=m·s²（实测 1.56107 / 3.12214 amu·Å²）。

机制链：`select_polyatomic_a` → 正卦限单位向量采样（6N 维相空间内逐步条件抽取
单位随机向量，按分量比例分配动量与坐标位移，每分量 50% 概率翻号；QMAX/QMIN 逐坐标
±0.1 Å 步进至 V−V(QZ)≥ENMT；PMAX=√(2·W·ENMT·C1)·PSCALE）→ 去质心速度 →
两步差值角动量插值（J′=J−Js，ω=I⁻¹J′）→ 能量闭合缩放。

- **转动**：`NROT_A=2` 选对称陀螺，固定 (J,K)=(2,0)，故
  EROTT=J(J+1)·C7²/(2·C1·I_perp) 解析确定（C7=0.063508 为 ℏ 的 code 单位值，
  C1=0.04184 为 kcal/mol→code 能量换算常数）。
- **能量闭合**：|H−HSCALE|/HSCALE ≥ 10⁻³ 时 P 与 (Q−QZ) 同乘 √(HSCALE/H)，
  NSCALE≤50；HSCALE=EROTT+ENMT_A。

（`PSCALE_A=5` 使初猜动能为主导；PSCALE=1 时软 E′ 模可致缩放循环过冲停机——
原版 VENUS 行为，非缺陷。）

## 预期与结果

- **能量闭合**：1000 条轨迹末态内能相对 HSCALE=10.185254 kcal/mol 的最大相对偏差
  9.99×10⁻⁴ ≤ 10⁻³ ✓
- **转动能量解析闭合**：程序 EROTT=0.1852540 vs 解析 0.1852524 kcal/mol
  （相对偏差 8.6×10⁻⁶）✓
- **符号对称**：9 个动量分量的正号率 ∈ [0.479, 0.518]，落在 50%±3σ 区间
  [0.44, 0.56]（N=1000）✓

补注：`CHOSEN EROTA`（0.157）低于目标 EROTT（0.185）——两步插值的转动份额损耗，
为已知行为；总内能闭合不受影响（EVIBA 与 HSCALE−EROTA 偏差 0.0007 kcal/mol）。
质心速度去除 max|Σm·v|/Σ|m·v|=1.5×10⁻⁵。

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量（stdout
> `INTERNAL ENERGY` 回显，活数据）：左=各轨迹初始内能（10.185-10.195 kcal/mol，
> 虚线=HSCALE=EROTT+ENMT_A=10.185——闭合含固定转动份额），右=逐轨迹缩放
> 收敛过程。
