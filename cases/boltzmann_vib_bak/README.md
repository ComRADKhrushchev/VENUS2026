# boltzmann_vib — 玻尔兹曼振动采样

> 运行：`cd cases/boltzmann_vib && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

本 case 验证玻尔兹曼振动采样（`INIT_SAMPLING_A=BOLTZMANN-VIB`，TVIB_A=2000 K，NT=1000）：
每个简正模的量子数应服从几何分布，且采样平均能量与温度定义一致。

## 2. 理论与体系

体系为等边 H₃（NATOMS=3，全为 1.008 u），TEST 势 LEPS（`src_TEST/test_potentials.f90`）。
LEPS 参数取自 input_qct.txt：De=4.746 eV、re=1.401 bohr、a=1.028 bohr⁻¹、Sato Δ=0.164
（源码 `src_TEST/test_potentials.f90:44-47` 定死单位为 eV/bohr/bohr⁻¹；input 注释
"[Å][Å⁻¹]"有误，数值一致仅单位标注错误）。等边极小处谱为
简并 E′ 对 1013.5 cm⁻¹ ×2 与呼吸模 1803.8 cm⁻¹（stdout 频谱表，见 results.txt）。
量子数抽样在 `THRMAN` 中实现：GAMA=−ln(U) 服从指数分布，n=trunc(GAMA/DUM)，
DUM=1.43878·ν/T=hν/k_BT，故 P(n)=(1−q)qⁿ，q=exp(−hν/k_BT)。

T=2000 K 时 q(E′)=0.482（双峰简并）、q(呼吸)=0.273，⟨n⟩=q/(1−q)=0.93/0.93/0.38。
单轨迹 EVIBA=Σ(n_i+½)ω_i（谐波能级和），解析期望 E[EVIBA]=12.815 kcal/mol
（kcal/mol→a.u. 换算 C1=0.04184，`venus_params.f90:60-80`）；随后能量闭合缩放至完整非谐 LEPS 势。

## 3. 方法与流程

1. 程序读取输入并回显采样模式与 VIBRATIONAL TEMPERATURE=2000.0（stdout，results.txt）。
2. 初始化随机数（ISEED=20260827），对 NT=1000 条轨迹由 `THRMAN` 逐模抽取几何分布量子数。
3. 按简正模组合初始坐标与动量，做能量闭合缩放（stdout `INTERNAL ENERGY` 回显）。
4. 每轨迹传播 NS=10 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

**判据**：三模量子数直方图应通过几何分布 KS 检验（Dmax ≤ 3/√N=0.0949），且各 q 值与温度定义一致。
**实测**（results.txt 直方图 n=0..5，模1 498/253/131/68/24/13，模2 512/251/131/62/27/9，
模3 741/192/54/9/3/1）：三模 KS 统计 Dmax=0.0197/0.0101/0.0142，全部远低于门限。
**结论**：几何分布抽样机制按设计工作。佐证：EVIBA 均值 12.702 vs 解析期望 12.815 kcal/mol
（相对偏差 0.0088，门限 2%），且逐轨迹 EVIBA=Σ(n_i+½)ω_i 成立（最大相对偏差 0.0010）。

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量（stdout `INTERNAL ENERGY` 回显，
> 即缩放后 ESEL）：左为各轨迹初始内能（5.5-30.9 kcal/mol 波动，逐模几何分布抽样加缩放闭合的直接体现），
> 右为逐轨迹缩放收敛过程。
