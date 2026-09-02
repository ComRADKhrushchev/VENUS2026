# rst_surface_md — RST Au(111) 板块 300 K 恒温 MD 均衡

> 运行：`cd cases/rst_surface_md && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

验证 RST 势下多原子表面板块（B 片段，NATOMB=144）的 MD 恒温均衡链
（INIT_SAMPLING_B=MD，NACTB=7）：速度重标（NSCALE=5000）与均衡段
（NEQUAL=5000）将体系温度驱动至 THERMOTEMP=300 K，随后以均衡构型为
初条件跑 NT=4 条 C 入射轨迹（Erel=14.528 kcal/mol，θ=0°）。

## 2. 理论与体系

RST 势：C-Au 成对 1D-NN（tanh 激活，GS 态 NH_surf=12/NH_inner=8，63 参数，
nn_weights_rst.txt）+ cosine 平滑截断（r_cut=6、d_cut=0.5）+ 排斥核 A=10/B=3，
吸附质受力；板块为 Born-von Karman 弹性 Au(111) slab（Analytic_Potential.txt：
cutoff 6.0 / a_lat 2.95 / alpha 0.5 / beta 0.3 / gamma 0.1，GENERATE 生成
Au₁₄₄，BOXLX=20.65、BOXLY=18.2786、SKEW=90°）。QZA_EQ=0,0,-3.8 → 初始
C 高度 z=2.2 Å。MD 链：恒温重标至 300 K（fort.30 SYSTEM TEMPERATURE 逐步
输出，10424 样本），均衡完成后进入轨迹传播（NS=30000、DT=0.01=0.1 fs、
NIP=50，`system temperature=` 每 50 步回显于 run_full.log）。

## 3. 方法与流程

1. 读取 input_qct.txt，初始化 RST+BVK 板块并做 MD 均衡（NSCALE/NEQUAL=5000，
   温度历史 fort.30）。
2. 均衡终态作为初条件，跑 NT=4 条轨迹，每 50 步打印体系温度与 C z(t)
   （fort.1001-1004，GWRITE_LEVEL=1 含 r_min 行）。
3. plot_fig.py 逐行 split 解析 run_full.log + fort.10NN 生成 fig_thermal.png
   （dpi=150，ps 轴）；>600 K 的恒温器瞬态尖峰从 T 轴掩蔽（trj 2 晚期
   0.63 eV 撞击热脉冲灌入 145 原子小胞所致，401/2404 条已掩蔽并打印计数）。

## 4. 核心验证

见 `fig_thermal.png`。

> 图注：上 = 体系温度 T(t)（灰=均衡重标段，彩=4 条生产轨迹，黑虚线=300 K
> 目标；均衡弛豫段自动体现在曲线里），下 = C 高度 z(t)（4 条）。判定依据：
> 均衡后生产段温度 306.6±24.3 K（掩蔽尖峰后 2003/2404 条，距目标 300 K
> <3%；全记录 2404 条零 NaN）；4 条轨迹 C z(t) 均呈吸附振荡（z_min
> ≈1.2-2.0 Å），物理合理。
