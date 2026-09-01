# 积分器矩阵 — 阶数收敛扫描与长时间能量漂移

> 运行：`cd cases/integrator_matrix && ../../venus_test.e`（本输入为 VERLET
> 10 ps 基线；完整 8 入口 × DT 矩阵由 sweep.py 驱动，结论见 §4）。

## 1. 目的

以简谐振子解析解为参照，判定 8 个积分器入口的阶数收敛性与长时间能量漂移，
并以 step-0 全矩阵逐位一致作为确定性判定。

## 2. 理论与体系

体系为双原子（2×1.008 u）气相 TEST 势 HARMONIC：每原子独立各向异性简谐阱
（k_x/k_y/k_z = 1/4/9，x0=0；test_potentials.f90:25-26），产生三个非简并简正
频率。基线 DT=0.01（10 fs）、NS=100000（input_qct.txt）。ω²=23.0605·C1·k/m=
0.965173·k/m（(10 fs)⁻²）；解析参照 x(t)=x₀cos(ωt)+(v₀/ω)sin(ωt)。

## 3. 方法与流程

1. sweep.py 生成 8 入口 × DT 矩阵；本目录 input_qct.txt 仅承载 VERLET 10 ps 基线
   （Tier B 双跑门禁判定）。
2. 阶数收敛扫描：总时长固定 T=1.0 code，DT 自 0.04 逐次折半，终态误差对解析解做
   log-log 拟合得有效阶数（results.txt 阶数收敛片段）。
3. 长时间漂移：DT=0.01 跑 10 ps，能量由 fort.8/fort.1001 坐标动量按
   E=ΣP²/2m+Σ½k_code(r−X0)² 重算（GWRITE 的 E0/T/H 在 ADIABATIC TEST 构建下为
   死值，不可用）；漂移排序与门禁结论记于 results.txt。

## 4. 核心验证

见 `fig_drift_matrix.png`。

> 图注：8 积分器入口 10 ps 能量漂移对比（log 轴柱状，数据取 results.txt 漂移表），
> 红色 = ADAMS（损坏）。判定依据：VERLET 的 DT 由 0.04 折至 0.005 时 err_pos 由
> 2.068×10⁻³ 降至 2.804×10⁻⁵（斜率 2.06 ≈ 2，二阶）；SYMPLECTIC-6/8 与
> RADAU-FIXED 达机器量化平台；ADAMS 斜率仅 1.01（有效一阶，损坏）；RK4 关键字为
> 不可达入口（LLL=0 被 RADAU 拒绝）。10 ps 漂移排序：RADAU-ADAPTIVE 1.0×10⁻⁶ <
> SYMPLECTIC-6/8/RADAU-FIXED 1.6×10⁻⁶ < BEEMAN 7.1×10⁻⁵ < VERLET 2.1×10⁻⁴ <
> SYMPLECTIC-4 2.8×10⁻⁴ < ADAMS 1.3×10⁰；step-0 全矩阵逐位一致。
