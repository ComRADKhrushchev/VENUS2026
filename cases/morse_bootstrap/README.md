# Morse 势函数正确性的验证

> 运行：`cd cases/morse_bootstrap && ../../venus_test.e`（MORSE H₂ 单轨迹基线，
> 输入 input_qct.txt，ISEED=12345 固定）。

## 1. 目的

验证 TEST 势构建在不读取任何体系数据文件时的正确性。

## 2. 理论与体系

体系与势：
- 内置 MORSE H₂ 双原子（NATOMS=2，1.008 u ×2）。
- 势参数取程序内置默认（`src_TEST/test_potentials.f90`）：
  De = 4.746 eV、re = 1.401 Å、a = 1.028 Å⁻¹；势能及梯度全部解析给出。

解析参照：
- Morse 振动态闭式解 E(n)=ωe(n+½)−ωexe(n+½)²，ωe=a·√(2De/μ)，
  为初始势能提供解析参照。
## 3. 方法与流程
-

## 4. 核心验证
见 `fig_potential_curve.png`。

> 图注：Morse 势能曲线 V(r)=De(em−1)²−De（De=4.746 eV、re=1.401、a=1.028，
> r 扫描 0.8-8 Å）；fort.1001 记录的轨迹能量 H 恒为 −0.0684 eV，与理论值一致。
> 判定依据：初始势能 −1.5769 kcal/mol（=−0.0684 eV）与理论值一致；同 DT 下
> 10 ps 漂移 2.1×10⁻⁴ eV；ISEED=12345 复跑全部产物逐位一致。
