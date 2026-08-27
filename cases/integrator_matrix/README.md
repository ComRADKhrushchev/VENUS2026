# 积分器矩阵 — 阶数收敛扫描与长时间能量漂移

> 运行：`cd cases/integrator_matrix && ../../venus_test.e`（本输入为 VERLET
> 10 ps 基线；完整 8 入口 × DT 矩阵由 sweep.py 驱动，结论见 §4）。

## 1. 目的

本 case 以简谐振子解析解为参照，判定 8 个积分器入口的阶数收敛性与
长时间能量漂移，并以 step-0 全矩阵逐位一致作为确定性锚。判据来源为
数值积分的一般收敛阶理论与保守系统长时间能量有界性。

## 2. 理论与体系

体系为双原子（各 1.008 u）气相 TEST 势 HARMONIC：每原子独立各向
异性简谐阱（三方向力常数不同，产生三个非简并简正频率）。势参数
（源码 src_TEST/test_potentials.f90:25-26；input 注释的 [eV/Å²][Å]
单位有误，以源码 bohr 制为准）：

| 参数 | 值 | 单位 | 出处 |
|---|---|---|---|
| k_x, k_y, k_z | 1, 4, 9 | eV/bohr² | HARM_K，test_potentials.f90:25 |
| 平衡位置 x0 | 0, 0, 0 | bohr | HARM_X0，test_potentials.f90:26 |
| 原子质量 | 1.008 × 2 | u | ATOM_MASSES |
| 基线 DT / NS | 0.01 / 100000 | 10 fs / 步 | input_qct.txt |

单位链：V=eV×23.0605 kcal/mol，再乘 C1 进 code 能量，C1=0.04184
（kcal/mol→code，venus_params.f90:69），故 ω²=23.0605·C1·k/m=
0.965173·k/m（(10 fs)⁻²）。解析参照 x(t)=x₀cos(ωt)+(v₀/ω)sin(ωt)。

## 3. 方法与流程

1. sweep.py 生成 8 入口 × DT 矩阵；本目录 input_qct.txt 仅承载
   VERLET 10 ps 基线（Tier B 双跑门禁锚）。
2. 阶数收敛扫描：总时长固定 T=1.0 code，DT 自 0.04 逐次折半，
   终态误差对解析解做 log-log 拟合得有效阶数（results.txt 阶数
   收敛片段）。
3. 长时间漂移：DT=0.01 跑 10 ps，能量由 fort.8/fort.1001 的坐标
   与动量按 E=ΣP²/2m+Σ½k_code(r−X0)² 重算（GWRITE 的 E0/T/H 在
   ADIABATIC TEST 构建下为死值，不可用）。
4. 汇总与门禁：漂移排序与逐位一致结论记于 results.txt；fort.999
   为单跑轨迹计数。

## 4. 核心验证：积分阶数收敛

判据：二阶方法（VERLET/BEEMAN）log-log 斜率≈2，高阶方法
（SYMPLECTIC-6/8、RADAU-FIXED）在最小 DT 达机器量化平台。实测
（results.txt 阶数收敛片段）：VERLET 的 DT 由 0.04 折至 0.005 时
err_pos 由 2.068×10⁻³ 降至 2.804×10⁻⁵，斜率 2.06；SYMPLECTIC-6/8
与 RADAU-FIXED 达量化平台；ADAMS 斜率仅 1.01。结论：阶数判定通过；
ADAMS 为有效一阶（损坏），RK4 关键字为不可达入口（LLL=0 被 RADAU
拒绝，RK 仅作 ADAMS 起步子）。佐证：10 ps 漂移排序 RADAU-ADAPTIVE
1.0×10⁻⁶ < SYMPLECTIC-6/8/RADAU-FIXED 1.6×10⁻⁶ < BEEMAN 7.1×10⁻⁵
< VERLET 2.1×10⁻⁴ < SYMPLECTIC-4 2.8×10⁻⁴ < ADAMS 1.3×10⁰
（results.txt 漂移表），且 step-0 全矩阵逐位一致（随机数流与
INTEGRATOR/DT/NS/NIP 解耦）。

> 图：`fig_drift_matrix.png` —— 8 积分器入口 10 ps 能量漂移对比
>（log 轴柱状，数据取 results.txt 漂移表），红色 = ADAMS（损坏）。
