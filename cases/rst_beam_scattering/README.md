# rst_beam_scattering — RST C/Au(111) 束流散射多轨迹

> 运行：`cd cases/rst_beam_scattering && ../../venus_test.e`（主输入文件 `input_qct.txt`，NT=20 存档态）

## 1. 目的

在 RST 势下以 Au(111) 板块复现束流散射链：法向入射（θ=0°）、瞄准点在斜胞内
均匀抽取（NRNDXY=1）、固定相对能 Erel=14.528 kcal/mol，NT=20 条轨迹检验
C 原子的吸附动力学（z(t) 下降-吸附振荡、r_min(t) 物理范围、能量守恒）。

## 2. 理论与体系

RST 势：C-Au 成对 1D-NN（tanh 激活，GS 态 NH_surf=12/NH_inner=8，63 参数，
nn_weights_rst.txt）+ cosine 平滑截断（r_cut=6、d_cut=0.5）+ 排斥核 A=10/B=3，
吸附质受力；板块为 Born-von Karman 弹性 Au(111) slab（Analytic_Potential.txt：
cutoff 6.0 / a_lat 2.95 / alpha 0.5 / beta 0.3 / gamma 0.1，GENERATE 生成
Au₁₄₄，BOXLX=20.65、BOXLY=18.2786、SKEW=90°）。QZA_EQ=0,0,-3.5 → 初始
C 高度 z=2.5 Å，vz<0 朝表面。NS=30000（DT=0.01=0.1 fs → 3000 fs）、NIP=50、
ISEED=20260821、GWRITE_LEVEL=1（fort.10NN 头含 E0/T/H、Q(C)、r_min）。

## 3. 方法与流程

1. 先以 smoke 检验通过（顶位吸附阱 z=2.90 Å / E_ads=−2.81 eV；解析梯度 vs FD
   最大差 5×10⁻⁹；BVK 弹性检验解析=FD 10⁻¹³）后正式跑 NT=20（run_full.log）。
2. plot_fig.py 逐行 split 解析 fort.1001-1020 叠绘 20 条 z(t) 与 r_min(t)
   → fig_scattering.png（dpi=150，ps 轴）。
3. 终止分类：块数=601（=NS/NIP+1）即跑满 → 吸附；<601 即散射提前终止。

## 4. 核心验证

见 `fig_scattering.png`。

> 图注：上 = 20 条轨迹 C 高度 z(t)（红虚线=顶位阱 z_eq=2.90 Å，图角注 E_ads=−2.81 eV
> vs 入射能 0.63 eV，阱深/入射能≈4.5× → 捕获必然），下 = C–Au 最近距离 r_min(t)。
> 判定依据：吸附 20/20（全部 601 块跑满，无散射提前终止）；z_min 全域
> 1.334-2.193 Å，末态 z 2.10-3.71 Å；全程 H 极差 0.0000 eV（微正则守恒）；
> r_min 全程 ≥1.45 Å，物理合理。
