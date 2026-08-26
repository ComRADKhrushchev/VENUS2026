# Morse 自举 — TEST 势体系无关性的首证

> 运行：`cd cases/morse_bootstrap && ../../venus_test.e`

## 验证目的

验证 TEST 势构建在不读取任何体系数据文件时的正确性：MORSE 势 H₂ 双原子单轨迹
（GAS-PHASE、ADIABATIC、VERLET，NS=5000），考察能量守恒与固定随机种子下的运行确定性。

## 原理

`TEST_PES=MORSE` 由程序内置 H₂ 参数驱动（De=4.746 eV、re=1.401 Å、a=1.028 Å⁻¹；
输入未显式指定 `MORSE_*` 关键字时取默认值），势能与梯度全部由
`src_TEST/test_potentials.f90` 解析给出，不依赖任何外部势文件。Morse 振动态有闭式解
E(n)=ωe(n+½)−ωexe(n+½)²、ωe=a·√(2De/μ)（μ 为约化质量），为能量量值提供解析锚。

初始条件：相对平移能 0 kcal/mol——纯势能起点，总能量等于初始势能
（−1.5769 kcal/mol = −0.0684 eV，与 MORSE 解析锚一致）。速度 Verlet 的能量漂移由
fort 坐标/动量重算验证（方法与结论见 `cases/integrator_matrix`：同 DT 下
10 ps 漂移 2.1×10⁻⁴ eV）。固定 `ISEED` 时全部输出逐位确定。

**输出行为登记**：GWRITE 的 E0/T/H(eV)（fort 头部）与逐步 KINETIC/TOTAL 打印
（stdout）在 VERLET 主循环路径下不刷新——VENUS.f90 积分后直接调 GWRITE、不调
ENERGY_1，故这些打印全程恒为采样末次 ENERGY_1 的陈旧值，**打印恒定不构成守恒
证据**；守恒验证以 fort 相空间重算为准。

## 预期与结果

- **相空间演化**：r_min(t) 在 0.73↔6.18 bohr 间振荡、|P|(t) 反相（r 峰对应 P 谷）
  ——纯势能起点的振子动力学，见下图 ✓
- **打印陈旧值行为**：见上文登记（恒定打印 = 采样期末值，非守恒证据）✓
- **确定性**：固定 `ISEED=12345` 双跑全部产物逐位一致 ✓

> 图：`fig_trajectory.png` —— 单轨迹相空间（fort.1001 活数据）：左=原子最近距离
> r_min(t)，右=动量模 |P|(t)；两量反相振荡为振子物理签名。
