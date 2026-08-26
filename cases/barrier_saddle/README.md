# 势垒鞍点起始 — BARRIER 模式简正模-反应坐标完整语义

> 运行：`cd cases/barrier_saddle && ../../venus_test.e`（NBAR=1 基线）
> 变体：`cp <变体文件> input_qct.txt && ../../venus_test.e`（差异见文末表）

## 验证目的

定案 `TASK=BARRIER`（NSELT=3）势垒激发模式的完整语义：鞍点起始几何、
简正模-反应坐标重排、NBAR=1/2/3 各分支与拒绝路径。

## 原理

体系：LEPS 3H 共线（NLINA=1；DE=4.746 eV、RE=1.401 Å、A=1.028 Å⁻¹、Δ=0.164、
m=1.008 u×3）；鞍点为指标-2 共线驻点 `QZA_EQ=(±1.6039, 0, 0; 0,0,0)`，
V_s=−129.6692 kcal/mol。

模重排（NSELT=3 专属）：NMA=3N−(6−NLINA)=4、NMBAR=3、IBARR=1——采样模
={0.5134 cm⁻¹ 转动伪模, 1225.26 反对称伸缩, 1791.12 对称伸缩}，
反应坐标模=弯曲虚模（打印 |ν|=449.18 cm⁻¹）。

NBAR 分支：NBAR=1 → PBAR=√(2·C1·EBAR)（反应坐标能量精确 ==EBAR）；
NBAR=2 → PBAR=√(−2·C5·TBAR·ln(1−RAND))（指数/Boltzmann 分布）；
NBAR=3 → Beyer-Swinehart 微正则（见 `cases/fixed_energy_rc`）。

## 预期与结果

- **谱与鞍点解析锚**：四模频率 449.18/449.18/1225.26/1791.12 cm⁻¹（容差 1），
  V(INITQP)=V_s=−129.6692 kcal/mol ✓
- **NBAR=1 精确 + NBAR=2 指数分布**：RC 能量==EBAR 精确；TBAR=2000 K 千轨迹
  均值 3.87 vs kT=3.97，KS D=0.0277<0.0515 ✓
- **退化/拒绝分支**：NBAR=0 退化 RC=0；NATOMB≠0 拒绝；NSELT=3+MB 强制终止
  （F17 修复后 loud STOP）✓

补注：哈密顿量极差 0（NS+1=201 步块）；模重排语义确认（伪转动模取大量子数、
量子数每轨迹重抽）；IJDIR 分支幅值 |q|=PBAR·|C|₂=0.907507。

> 图：`fig_trajectory.png` —— 入库轨迹相空间（fort.1001-1002 活数据）：左=原子 1
> 到最近原子距离 r_min(t)——NS=200 步内 1.47→2.45→1.47 bohr 一次往返（鞍点起始、
> EBAR 反应坐标推力下的演化），右=|P|(t)。GWRITE 能量打印为陈旧值（见
> morse_bootstrap 登记）。

| 变体 | 差异 |
|---|---|
| input_nbar2.txt | NBAR=2, TBAR=2000 K, NT=1000/NS=50（指数分布） |
| input_ijdir.txt | IJDIR=1 IDIR=1 JDIR=3, NT=100/NS=50 |
| input_ijdir_nofix.txt | IJDIR=1 无 IDIR/JDIR（F16 越界读探针） |
| input_nbar0.txt | NBAR=0 退化（无 EBAR 行） |
| input_natomb.txt | NATOMA=1/NATOMB=2（拒绝断言） |
| input_mb.txt | INIT_SAMPLING_A=MB（F17 强制终止断言） |
