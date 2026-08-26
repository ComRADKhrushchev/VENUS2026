# 积分器矩阵 — 阶数收敛扫描 + 长时间能量漂移

> 运行：`cd cases/integrator_matrix && ../../venus_test.e`（VERLET 10 ps 基线；
> 完整 8 入口 × DT 矩阵由内部扫描驱动，结论见下）

## 验证目的

以简谐振子解析解为参照，验证 8 个积分器入口的阶数收敛性与长时间能量漂移，
并以 step-0 全矩阵逐位一致作确定性锚。

## 原理

解析参照：x(t)=x₀·cos(ωt)+(v₀/ω)·sin(ωt)，ω²=k_code/m，
k_code=23.0605·C1·k=0.965173·k（单位 (10 fs)⁻²）。

- **阶数收敛**：总时长固定，DT={0.04, 0.02, …, 0.00125} 扫描，log-log 斜率
  vs 解析终态给出有效阶数。
- **长时间漂移**：DT=0.01，10 ps（NS=100000）；能量由 fort 坐标/动量重算
  E=ΣP²/2m+Σ½k_code·(r−X0)²（GWRITE 的 E0/T/H 打印在 ADIABATIC TEST 构建下
  为死值，不用）。

关键发现：**RK4 关键字为不可达入口**（LLL=0 被 RADAU 拒绝，RK 仅作 ADAMS
起步子）；**ADAMS 损坏**（有效 1 阶、漂移最差）；SYMPLECTIC-4 系数本身仅 2 阶；
RADAU-ADAPTIVE 为容差驱动，只入漂移不入扫描。

## 预期与结果

- **阶数**：VERLET/BEEMAN/SYMPLECTIC-4 斜率≈2；SYMPLECTIC-6/8 与 RADAU-FIXED
  在最小 DT 达机器量化平台（更高阶）；ADAMS 斜率≈1 = 损坏 ✓
- **漂移排序**（10 ps）：RADAU-ADAPTIVE(1.0e-06) < SYMPLECTIC-6/8/RADAU-FIXED
  (1.6e-06) < BEEMAN(7.1e-05) < VERLET(2.1e-04) < SYMPLECTIC-4(2.8e-04)
  < ADAMS(1.3e+00，损坏) ✓
- **确定性锚**：step-0 全矩阵逐位一致（随机数流与 INTEGRATOR/DT/NS/NIP 解耦）✓

补注：RK4 与 ADAMS 已列为禁用组合；`NIP` 须为 7 的倍数（ADAMS 宏步 7 子步，
MOD(NC,NIP) 输出门）。

> 图：`fig_drift_matrix.png` —— 8 积分器入口 10 ps 能量漂移对比（log 轴柱状，
> 数据取 results.txt 漂移表）：RADAU-ADAPTIVE 1.0×10⁻⁶ 至 ADAMS 1.3×10⁰ eV
> 跨六个量级；红色 = ADAMS（损坏，已移除）。
