# ebk_fixed_nj — 固定 n,J 的 EBK 量子态采样

> 运行：`cd cases/ebk_fixed_nj && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

验证 EBK 半经典量化在双原子分子上固定振转量子态采样的正确性：固定 (n,J)
下态能量不变、核间距服从 WKB 概率密度、振转能量与解析值闭合。

## 2. 理论与体系

采样流程：
- `NACTA=0`（双原子 MB/EBK 分支）+ `TRV_A<0`（固定量子数标志）→ `INITEBK` 采样。
- r 按 WKB 经典密度采样：p(r)∝1/|p_r(r)|。
- 径向动量 PR=√(2μ·SUMM(r))（SUMM 为有效势），符号以 50% 概率翻转。
- 经典可达域由 EBK 转折点 [RMIN, RMAX] 限定。

能级（含非谐修正与刚性转子项）：

  EVIBA = hc[ωe(n+½) − ωexe(n+½)²] + hc·Be·J(J+1)

参数与理论对照：

- MORSE 谱常数 ωe=2326.4 cm⁻¹、ωexe=35.35 cm⁻¹。
- n=3、J=2 时理论对照值 EVIBA=22.3347 kcal/mol。

## 3. 方法与流程

1. 读取输入并打印 TRV/TROT/N/J 量子数与频谱表（stdout，results.txt）。
2. 每轨迹传播 NS=10 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

见 `fig_sampling_stats.png`。

> 图注：三个面板分别为 EVIBA（kcal/mol）、EROTA（kcal/mol）的轨迹分布与

> 初始键长 r 的概率密度（NT=999 有效样本）。EVIBA 恒为 22.3134（固定量子数

> n=3 的直接体现）；r 分布由 WKB 密度 p(r)∝1/|p_r(r)| 支配。

> 判定依据：ENJA/RMIN/RMAX 跨样本相对展宽全为 0（态不变性）；r 分布 KS

> Dmax=0.0254 ≤ 3/√N=0.0949，可达域 [1.0431, 1.9820] 与解析转折点一致；EVIBA

> 均值 22.3134 vs 理论对照值 22.3347 kcal/mol（相对偏差 9.6×10⁻⁴，门限 0.5%）。
