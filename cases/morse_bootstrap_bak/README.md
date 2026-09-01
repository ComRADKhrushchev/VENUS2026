# Morse 自举 — TEST 势体系无关性的首证

> 运行：`cd cases/morse_bootstrap && ../../venus_test.e`（MORSE H₂ 单轨迹基线，
> 输入 input_qct.txt，ISEED=12345 固定）。

## 1. 目的

验证 TEST 势构建在不读取任何体系数据文件时的正确性：MORSE 势 H₂ 双原子
单轨迹（GAS-PHASE、ADIABATIC、VERLET、NS=5000），考察能量量值的解析锚、
相空间演化与固定随机种子下的运行确定性。

## 2. 理论与体系

体系为内置 MORSE H₂ 双原子（NATOMS=2，ATOM_MASSES=1.008 u ×2，`TEST_PES=MORSE`）。
势参数取程序内置默认（`src_TEST/test_potentials.f90`）：De = 4.746 eV、
re = 1.401 Å、a = 1.028 Å⁻¹；势能及其梯度全部解析给出，不依赖外部势文件。
Morse 振动态有闭式解 E(n)=ωe(n+½)−ωexe(n+½)²，其中 ωe=a·√(2De/μ)
（μ 为约化质量），为初始势能提供解析参照。初始条件为相对平动能 0 kcal/mol
（`ELEC_METHOD=ADIABATIC` 路径下纯势能起点）：原子初始间距 6.00 Å，
总能量即初始势能。能量单位换算常数 C1 = 0.04184（kcal/mol → a.u.，
`venus_params.f90:69`），本卡与各采样卡共用；采样链中另引用
C5 = 0.083144×10⁻³ a.u.（Boltzmann 温度常数，`venus_params.f90:73`）与
C7 = 0.063508（ℏ 的 code 单位值，`venus_params.f90:75`），此处一并登记备引。

## 3. 方法与流程

1. 程序以纯势能起点组装初始构型（回显见 stdout 横幅，摘录于 results.txt）。
2. VERLET 主循环积分 5000 步、DT=0.01（单位 10 fs），逐步打印
   KINETIC/POTENTIAL/TOTAL ENERGY（stdout）。
3. 相空间逐步落盘（`fort.1001`，坐标动量 + E0/T/H(eV) 头）与采样末能量
   汇总（`fort.999`、`fort.666`）。
4. 事后由 fort 坐标/动量重算能量漂移（方法见 `cases/integrator_matrix`）
   并核对双跑确定性。

**输出行为登记**：GWRITE 的 E0/T/H(eV)（fort 头部）与 stdout 逐步
KINETIC/POTENTIAL/TOTAL 打印在 VERLET 主循环下不刷新——VENUS.f90 积分后
直接调 GWRITE、不重调 ENERGY_1（对照 VENUS.f90:403 与 ：448），故这些打印
全程恒为采样末次 ENERGY_1 的陈旧值（kcal/mol），打印恒定不构成守恒证据；
守恒验证以 fort 相空间重算为准（其他 case 引用本条）。

## 4. 核心验证：纯势能起点的振子相空间演化

这是判定 TEST 势链在工作的单一关键证据，因为它同时要求三个环节全部正确：
MORSE 势参数解析正确（否则势能量值偏离锚点）、梯度正确（否则轨迹不构成
振子）、积分器正确（否则能量漂移破坏振荡结构）。

**判据**：r_min(t) 呈全幅振荡、|P|(t) 反相（r 峰对应 P 谷）；初始势能
−1.5769 kcal/mol（=−0.0684 eV）与 MORSE 解析锚一致；fort 重算 10 ps
漂移与 integrator_matrix 同 DT 水平相当。

**实测**（results.txt + fort.1001）：初始势能 −1.576872 kcal/mol；r_min 由
6.1768 bohr 起、末步 5.6324 bohr，全程在 0.73↔6.18 bohr 间振荡、|P| 反相
（fig_trajectory.png）；同 DT 下 10 ps 漂移 2.1×10⁻⁴ eV
（cases/integrator_matrix）；ISEED=12345 双跑全部产物逐位一致。结论：通过。
配套佐证：fort.666 回显势能 ±0.0317052 kcal/mol 对（等距双点）；
fort.999 计 1 条轨迹、耗时 0.203 s。E0/T/H(eV) 打印全程纹丝不动
（−0.06838）而坐标动量持续演化，即登记的陈旧值行为的直接证据。

> 图：`fig_trajectory.png` —— 单轨迹相空间（fort.1001 活数据）：左 = 原子
> 最近距离 r_min(t)，右 = 动量模 |P|(t)；两量反相振荡为振子物理签名。
