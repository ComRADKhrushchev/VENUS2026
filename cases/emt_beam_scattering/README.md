# emt_beam_scattering — EMT-NN C/Au(111) 束流散射多轨迹

> 运行：`cd cases/emt_beam_scattering && ../../venus_test.e`（主输入文件 `input_qct.txt`，NT=20 存档态）

## 1. 目的

在 EMT-NN 势下以真实 Au(111) 板块复现束流散射链：法向入射（θ=0°）、瞄准点
在斜胞内均匀抽取（NRNDXY=1）、固定相对能 Erel=14.528 kcal/mol，NT=20 条轨迹
检验 C 原子的散射/吸附动力学（z(t) 下降-反弹/吸附、r_min(t) 物理范围、无 NaN）。

## 2. 理论与体系

体系为 C(A) + Au₁₄₄ 板块(B)（12×12 六层斜胞，SKEW=60°、BOXLX/BOXLY=17.7 Å，
周期边界）。EMT-NN 势：EMTNN_ALAT=2.88 Å，板块 Slab.xyz，权重
nn_weights_emt_nn.txt（均取自 ../../data/emt_nn/）。QZA_EQ=0,0,-3.5 叠加束流
分离 S=6.0 → 初始 C 高度 z=2.5 Å，vz<0 朝表面。NS=30000（DT=0.01=0.1 fs →
3000 fs）、NIP=50、ISEED=20260821、GWRITE_LEVEL=1（fort.10NN 头含 E0/T/H、
Q(C)、Au_nearest/r_min）。终止判据：NC>50 后 z≥RMAX=6 且离面 → 散射完成；
跑满 NS → 吸附。

## 3. 方法与流程

1. 先以 NT=2 冒烟（80.5 s）验证无 NaN、z(t) 物理行为（两轨迹均下降-反弹-吸附，
   z_min=1.90/1.30 Å），通过后改 NT=20 正式跑（806.7 s，运行日志
   run_production_nt20.log）。
2. plot_fig.py 解析 fort.1001-1020 叠绘 20 条 z(t) 与 r_min(t)
   → fig_scattering.png。
3. 终止分类：轨迹块数=601（=NS/NIP+1）即跑满 → 吸附；<601 即散射提前终止。

## 4. 核心验证

见 `fig_scattering.png`。

> 图注：上 = 20 条轨迹 C 高度 z(t)，下 = C–Au 最近距离 r_min(t)（ps 轴）。
> 判定依据：20/20 轨迹无 NaN；全部 601 块跑满 → 吸附 20、散射反弹 0
>（末态 z=1.34-3.06 Å，远低于 RMAX=6）；z_min 全域 1.131-1.899 Å（均值
> 1.320 Å），r_min 全程 ≥1.13 Å，与 C-Au 平衡键长（~1.9-2.1 Å）及入射
> 动能压缩深度物理一致；z(t) 普遍呈下降-反弹-振荡吸附行为。
