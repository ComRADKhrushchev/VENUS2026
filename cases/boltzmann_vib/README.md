# boltzmann_vib — 玻尔兹曼振动采样

> 运行：`cd cases/boltzmann_vib && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

验证玻尔兹曼振动采样（`INIT_SAMPLING_A=BOLTZMANN-VIB`，TVIB_A=2000 K，NT=1000）：
每个简正模量子数服从几何分布，且采样平均能量与温度定义一致。

## 2. 理论与体系

体系为等边 H₃（NATOMS=3，全为 1.008 u），TEST 势 LEPS；等边极小处谱为简并 E′ 对
1013.5 cm⁻¹ ×2 与呼吸模 1803.8 cm⁻¹（stdout 频谱表，results.txt）。抽样在 `THRMAN`
实现：GAMA=−ln(U) 服从指数分布，n=trunc(GAMA/DUM)，DUM=1.43878·ν/T=hν/k_BT，
故 P(n)=(1−q)qⁿ，q=exp(−hν/k_BT)。T=2000 K 时 q(E′)=0.482（双峰简并）、
q(呼吸)=0.273，⟨n⟩=q/(1−q)=0.93/0.93/0.38；理论期望 E[EVIBA]=12.815 kcal/mol。

## 3. 方法与流程

1. 读取输入并打印采样模式与 VIBRATIONAL TEMPERATURE=2000.0（stdout，results.txt）。
2. 初始化随机数（ISEED=20260827），NT=1000 条轨迹由 `THRMAN` 逐模抽取几何分布
   量子数，按简正模组合坐标动量并做能量闭合缩放（stdout `INTERNAL ENERGY`）。
3. 每轨迹传播 NS=10 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

见 `fig_sampling_stats.png`。

> 图注：左 = 各轨迹初始内能（5.5-30.9 kcal/mol，逐模几何分布抽样加缩放闭合的
> 直接体现），右 = 逐轨迹缩放收敛过程（stdout `INTERNAL ENERGY`，10 条轨迹）。
> 判定依据：三模量子数直方图 KS Dmax=0.0197/0.0101/0.0142，全部远低于门限
> 3/√N=0.0949；EVIBA 均值 12.702 vs 理论值 12.815 kcal/mol（相对偏差 0.0088，
> 门限 2%），逐轨迹 EVIBA=Σ(n_i+½)ω_i 成立（最大相对偏差 0.0010）。
