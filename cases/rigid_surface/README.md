# rigid_surface — 表面振子（GLO）束流散射

> 运行：`cd cases/rigid_surface && ../../venus_test.e`（NSURF=2 + NGLO=1
> 表面振子模型，输入 input_qct.txt）。

## 1. 目的

验证表面振子（surface-oscillator / GLO）分支：表面 Au 原子经弹簧连接
ghost 耗散原子（广义 Langevin 振子，Tully GLE），气体 H 原子束流入射。
检验束流组装、振子热初条件采样（GLOSELECT）与振子-气体耦合动力学。

## 2. 理论与体系

体系与布局：
- H（气体，片段 A）+ 表面 Au + GLO ghost，NATOMB=0，
  NATOMS = NATOMA + 2（venus_input 的 GLO 布局检查）。
- 质量分配：H 1.008 u、表面 Au 与 ghost 各 196.97 u。

GLO 模型（`src_VENUS/GLO.f90`）：
- 表面原子-ghost 弹簧频率 WS1 = ghost 频率 WG1 = 52.1 cm⁻¹（AU 输入，
  GLOINIT 内 ×413.4 换算）；ghost 摩擦 FCG；采样温度 TVIB_B。
- GLOSELECT 按 TVIB_B 高斯采样振子位移/动量；ghost 受摩擦 + 白噪声
  （广义 Langevin 耗散）。
- 弹簧力完全由 GLO 方程描述；TEST 势（HARMONIC）对振子原子不加阱，
  仅描述气体原子。

## 3. 方法与流程

1. 运行 NT=20、NS=20000（2 ps）轨迹，输出 fort.1001-1020。
2. plot_fig.py 绘制 H 的 z(t) 叠图与表面 Au 振幅。

## 4. 核心验证

见 `fig_long_trajectory.png`。

> 图注：20 条轨迹的 H z(t)（入射高度 8.0 Å）与表面 Au 的 z 振荡。
> 判定依据：全部轨迹数值有限；表面 Au 振子被激发（振幅量级与 52 cm⁻¹
> 弱弹簧及入射能转移自洽）；HARMONIC 测试势对 H 无近距排斥，H 穿越
> z=0 平面到达镜像位置（测试势固有属性，非缺陷；真实表面散射用 RST 势）。
