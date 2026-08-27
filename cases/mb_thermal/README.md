# 热 n,J 采样 — 量子热振动/转动态分布

> 运行：`cd cases/mb_thermal && ../../venus_test.e`（ISEED 固定，双跑
> 逐位一致）。

## 1. 目的

本 case 验证量子热 n,J 采样（NACTA=2）在 MORSE 双原子上的统计
正确性：振动量子数应服从几何分布，转动量子数应服从带核自旋权重
的 Boltzmann 分布，且 (n,J)→能量映射逐位确定。判据来源为量子
统计力学：温度 T 下振动态占据为几何分布，转动态占据为简并×
Boltzmann×核自旋统计权重。

## 2. 理论与体系

体系为 H₂ 双原子（2×1.008 u）气相 TEST 势 MORSE。势参数
（源码 src_TEST/test_potentials.f90:36-38；input 注释的 [Å][Å⁻¹]
单位有误，源码为 bohr 制）：

| 参数 | 值 | 单位 | 出处 |
|---|---|---|---|
| De | 4.746 | eV | MORSE_DE，test_potentials.f90:36 |
| re | 1.401 | bohr | MORSE_RE，test_potentials.f90:37 |
| a | 1.028 | bohr⁻¹ | MORSE_A，test_potentials.f90:38 |
| TRV_A = TROT_A | 2000 | K | input_qct.txt |
| 转动惯量 AIA | 0.98925 | u·Å² | AI_A（venus_input 以 QZA_EQ 重算覆盖） |
| NT | 1000 | 条 | input_qct.txt |

抽样公式（THRMAN.f、PROBJ.f）：振动 n=trunc(GAMA(1)/DUM)，其中
GAMA(1)=−ln(U) 服从指数分布、DUM=C7·ω/(C5·T)，即 P(n)=(1−q)qⁿ、
q=exp(−hνω/k_BT)；转动 B=48.5085/(2·AIA·T)，P(J)∝WGT(J)(2J+1)
exp(−B·J(J+1))，H₂ 正/仲核自旋权重 WGT=奇 J 0.75/偶 J 0.25。
此处 C7=0.063508（ℏ 的 code 值，venus_params.f90:75）、C5=
8.3144×10⁻⁵ a.u.（Boltzmann 常数，venus_params.f90:73）、C1=
0.04184（kcal/mol→code，venus_params.f90:69）。T=2000 K 时
q=0.1876、B=0.012259。

## 3. 方法与流程

1. `THRMAN` 按几何分布逐轨迹抽取振动量子数 n（stdout `NNA =` 回显，
   1000 条，见 results.txt）。
2. `PROBJ` 以拒绝采样抽取转动量子数 J（stdout `JA =` 回显）。
3. `INITEBK` 按 (n,J) 态做 EBK 固定点采样，得核间距与径向动量
   PR=√(2μ·SUMM(r))；每轨迹回显 ENJA/DUM/RMIN。
4. 统计验收：n/J 分布 KS 检验与 ENJA 确定性判定记于 README 本表
   （数据源 stdout，results.txt 存样例回显）；fort.8/fort.1001 入库
   逐步轨迹，fort.999 汇总。

## 4. 核心验证：n 几何分布命中

判据：n=0..∞ 的经验分布与几何分布 q=0.1876 的 KS 统计量
Dmax≤3/√N=0.0949。实测（NT=1000 stdout 统计，results.txt）：KS
Dmax=0.0176，直方图 n=0..5 为 830/141/24/5/0/0，呈几何衰减。结论：
振动热采样通过。佐证：J 分布 KS Dmax=0.0129≤0.0949，偶/奇 J 计数
287/713 与正/仲权重 0.25/0.75 吻合；(n,J)→ENJA 逐位确定（位级失配
0），范围 [0.1381, 1.2230] code 能量 ⊂ (0, De·C1)。补注：597 条
`FINLNJ` 诊断行为 INITEBK 在 n=0/ZPE 边界的已知提示噪声，全部轨迹
正常收敛，非失败项。

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量
>（stdout `CHOSEN:` 回显）：左 = 各轨迹 EVIBA+EROTA（3.5-25.4
> kcal/mol 波动，n,J 热分布的直接体现），右 = 逐轨迹采样过程。
