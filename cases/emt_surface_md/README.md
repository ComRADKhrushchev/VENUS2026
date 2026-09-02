# emt_surface_md — EMT-NN Au(111) 板块 300 K 恒温 MD 均衡

> 运行：`cd cases/emt_surface_md && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

验证 EMT-NN 势下多原子表面板块（B 片段，NATOMB=144）的 MD 恒温均衡链
（INIT_SAMPLING_B=MD，NACTB=7）：速度重标（NSCALE=5000）与均衡段
（NEQUAL=5000）将体系温度驱动至 THERMOTEMP=300 K 附近，随后以均衡构型为
初条件跑 NT=4 条 C 入射轨迹（Erel=14.528 kcal/mol，θ=0°）。

## 2. 理论与体系

体系为 C(A) + Au₁₄₄ 六层 12×12 斜胞板块(B)（Slab.xyz，SKEW=60°、
BOXLX/BOXLY=17.7 Å），EMT-NN 势（EMTNN_ALAT=2.88 Å，权重
nn_weights_emt_nn.txt）。QZA_EQ=0,0,-3.8 → 初始 C 高度 z=2.3 Å（S=6.0）。
MD 链：恒温重标至 300 K（fort.30 SYSTEM TEMPERATURE 逐步输出），均衡完成后
进入轨迹传播（NS=30000、DT=0.01=0.1 fs、NIP=50，`system temperature=` 每 50 步
回显于 run.log）。

## 3. 方法与流程

1. 读取 input_qct.txt，初始化 EMT-NN 板块并做 MD 均衡（NSCALE/NEQUAL=5000，
   温度历史 fort.30，共 10424 样本）。
2. 均衡终态作为初条件，跑 NT=4 条轨迹，每 50 步打印体系温度与 C z(t)
   （fort.1001-1004，GWRITE_LEVEL=1 含 r_min 行）。
3. plot_fig.py 解析 fort.30 + run.log + fort.10NN 生成 fig_thermal.png。

## 4. 核心验证

见 `fig_thermal.png`。

> 图注：上 = 体系温度 T(t)（蓝=均衡段 fort.30，彩=4 条生产轨迹 run.log，黑虚线
> =300 K 目标），下 = C 高度 z(t)（4 条）。判定依据：均衡段 T 均值
> 295.3±20.6 K（后半段 296.2 K，距目标 <2%）；生产段 T 均值 306.7±21.0 K
>（C 撞击沉积 14.528 kcal/mol 导致 ~7 K 温升，物理合理）；全部坐标无 NaN，
> C z(t) 均呈吸附振荡（z_min 1.32-1.78 Å）。
