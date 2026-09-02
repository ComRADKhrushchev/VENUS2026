# Morse 自举 — TEST 势体系无关性的首证

> 运行：`cd cases/morse_bootstrap && ../../venus_test.e`（MORSE H₂ 单轨迹基线，
> 输入 input_qct.txt，ISEED=12345 固定）。

## 1. 目的

验证 TEST 势构建在不读取任何体系数据文件时的正确性：MORSE 势 H₂ 双原子单轨迹
（GAS-PHASE、ADIABATIC、VERLET、NS=5000），考察能量量值的理论对照、相空间演化
与固定随机种子下的运行确定性。

## 2. 理论与体系

体系为内置 MORSE H₂ 双原子（NATOMS=2，1.008 u ×2）。势参数取程序内置默认
（`src_TEST/test_potentials.f90`）：De = 4.746 eV、re = 1.401 Å、a = 1.028 Å⁻¹；
势能及梯度全部解析给出。Morse 振动态闭式解 E(n)=ωe(n+½)−ωexe(n+½)²，
ωe=a·√(2De/μ)，为初始势能提供解析参照。初始条件为相对平动能 0（纯势能起点）：
原子初始间距 6.00 Å，总能量即初始势能。

## 3. 方法与流程

1. 程序以纯势能起点组装初始构型（stdout 横幅，摘录于 results.txt）。
2. VERLET 主循环积分 5000 步、DT=0.01（10 fs），逐步打印能量；相空间逐步落盘
   （`fort.1001`）并汇总（`fort.999`、`fort.666`）。
3. 事后由 fort 坐标/动量重算能量漂移（方法见 `cases/integrator_matrix`）并核对
   双跑确定性。

## 4. 核心验证


见 `fig_potential_curve.png`。

> 图注：左 = 原子最近距离 r_min(t)，右 = 动量模 |P|(t)；两量反相振荡为振子物理
> 签名（fort.1001 单轨迹）。判定依据：初始势能 −1.576872 kcal/mol 与理论对照值
> −1.5769（=−0.0684 eV）一致；同 DT 下 10 ps 漂移 2.1×10⁻⁴ eV，与
> integrator_matrix 水平相当；ISEED=12345 双跑全部产物逐位一致。
