# Morse 势函数正确性的验证

> 运行：`cd cases/morse_bootstrap && ../../venus_test.e`（MORSE H₂ 单轨迹基线，
> 输入 input_qct.txt，ISEED=12345 固定）。

## 1. 目的

验证 TEST 势构建在不读取任何体系数据文件时的正确性

## 2. 理论与体系

体系为内置 MORSE H₂ 双原子（NATOMS=2，1.008 u ×2）。势参数取程序内置默认
（`src_TEST/test_potentials.f90`）：De = 4.746 eV、re = 1.401 Å、a = 1.028 Å⁻¹；
势能及梯度全部解析给出。Morse 振动态闭式解 E(n)=ωe(n+½)−ωexe(n+½)²，
ωe=a·√(2De/μ)，为初始势能提供解析参照。

## 3. 方法与流程
-

## 4. 核心验证
见 `fig_potential_curve.png`。

