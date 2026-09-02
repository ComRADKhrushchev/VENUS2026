# 热 n,J 采样 — 量子热振动/转动态分布

> 运行：`cd cases/mb_thermal && ../../venus_test.e`（ISEED 固定，双跑逐位一致）。

## 1. 目的

验证量子热振转采样在双原子分子上的统计正确性：振动量子数服从几何分布、
转动量子数服从带核自旋权重的 Boltzmann 分布、(n,J) 到能量的映射逐位确定。

## 2. 理论与体系

TEST 势 MORSE（test_potentials.f90:36-38）；TRV_A=TROT_A=2000 K、NT=1000
（input_qct.txt）；转动惯量 AIA=0.98925 u·Å²（venus_input 以 QZA_EQ 重算覆盖）。
抽样公式（THRMAN.f、PROBJ.f）：振动 n=trunc(GAMA(1)/DUM)，GAMA(1)=−ln(U) 服从
指数分布、DUM=ℏω/k_BT，即 P(n)=(1−q)qⁿ、q=exp(−hν/k_BT)；转动 B=48.5085/(2·AIA·T)，
P(J)∝WGT(J)(2J+1)exp(−B·J(J+1))，H₂ 正/仲核自旋权重 WGT=奇 J 0.75/偶 J 0.25。
T=2000 K 时 q=0.1876、B=0.012259。

## 3. 方法与流程

1. `THRMAN` 按几何分布逐轨迹抽取 n（stdout `NNA =`，1000 条，results.txt）；
   `PROBJ` 以拒绝采样抽取 J（stdout `JA =`）。
2. `INITEBK` 按 (n,J) 态做 EBK 固定点采样，得核间距与径向动量 PR=√(2μ·SUMM(r))，
   每轨迹打印 ENJA/DUM/RMIN；fort.8/fort.1001 输出逐步轨迹，fort.999 汇总。

## 4. 核心验证

见 `fig_sampling_stats.png`。

> 图注：左面板为振动量子数 n 的经验直方图（NT=1000）与几何分布理论线
> （q=0.1876），右面板为转动量子数 J 的分布与核自旋加权理论线
> P(J)∝WGT(2J+1)e^{−BJ(J+1)}（正/仲权重 0.75/0.25）。
> 判定依据：n 分布 KS Dmax=0.0176 ≤ 3/√N=0.0949（直方图 n=0..5 为
> 830/141/24/5/0/0）；J 分布 KS Dmax=0.0129，偶/奇 J 计数 287/713 与权重吻合；
> (n,J)→ENJA 逐位确定，范围 [0.1381, 1.2230] ⊂ (0, De)。
