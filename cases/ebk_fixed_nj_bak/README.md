# ebk_fixed_nj — 固定 n,J 的 EBK 量子态采样

> 运行：`cd cases/ebk_fixed_nj && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

本 case 验证 EBK（Einstein-Brillouin-Keller）固定量子态采样在 MORSE 双原子上的正确性：
固定振动量子数 n=3、转动量子数 J=2，采样 999 条有效轨迹（NT=1000），考察 EBK 态不变性、
核间距分布的 WKB 密度命中与振转能量解析闭合。

## 2. 理论与体系

体系为 TEST 势 MORSE 双原子 H₂（NATOMS=2，各 1.008 u，约化质量 μ=0.5040 u）。
MORSE 参数取自 input_qct.txt：De=4.746 eV、re=1.401 bohr、a=1.028 bohr⁻¹
（源码 `src_TEST/test_potentials.f90:35-37` 定死单位为 eV/bohr/bohr⁻¹；input 注释
"[Å][Å⁻¹]"有误，数值一致仅单位标注错误）。采样链：`NACTA=0`（双原子 MB/EBK 分支）+ `TRV_A<0`（固定量子数标志）
→ `INITEBK` 逐轨迹采样。采样按 WKB 经典密度进行：p(r)∝1/|p_r(r)|，
径向动量 PR=√(2μ·SUMM(r))（SUMM 为有效势），符号以 50% 概率翻转；经典可达域由 EBK 能量
转折点 [RMIN, RMAX] 限定。EBK 能级（含非谐修正与刚性转子项）：

  EVIBA = hc[ωe(n+½) − ωexe(n+½)²] + hc·Be·J(J+1)

MORSE 谱常数 ωe=2326.4 cm⁻¹、ωexe=35.35 cm⁻¹；n=3、J=2 时解析 EVIBA=22.3347 kcal/mol
（kcal/mol→a.u. 换算 C1=0.04184，`venus_params.f90:60-80`）。

## 3. 方法与流程

1. 程序读取输入并回显 TRV/TROT/N/J 量子数与频谱表（stdout，results.txt）。
2. 初始化随机数（ISEED=20260820），对 NT=1000 条轨迹由 `INITEBK` 按 WKB 密度采样 r、
   组装径向/转动动量，回显 CHOSEN EROTA/EVIBA（results.txt）。
3. 每轨迹传播 NS=10 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

**判据**：固定 n,J 下 EBK 能量与转折点跨样本应为常数（态不变性），核间距分布应命中
WKB 密度（KS Dmax ≤ 3/√N=0.0949）。**实测**：ENJA/RMIN/RMAX 跨 999 个采样相对展宽全为 0；
r 分布 KS 统计 Dmax=0.0254，实测可达域 [1.0431, 1.9820] 与解析转折点 [1.0430, 1.9820] 一致；
单样本 CHOSEN EVIBA=22.313 kcal/mol（results.txt 第 1 条轨迹）。**结论**：固定量子态 EBK
采样按设计工作。佐证：程序 EVIBA 均值 22.3134 vs 解析 22.3347 kcal/mol（相对偏差
9.6×10⁻⁴，门限 0.5%）；PR 符号正号率 0.530（N=999 时 3σ≈±0.016），符合 50% 随机翻转语义。

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量（stdout `CHOSEN:` 回显）：
> 左为各轨迹 EVIBA+EROTA（22.46-22.84 kcal/mol，EVIBA 部分恒为 22.3134，固定量子数语义的直接体现），
> 右为逐轨迹采样过程。
