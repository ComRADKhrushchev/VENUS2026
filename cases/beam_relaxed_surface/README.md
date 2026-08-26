# 松弛表面束流几何（NSURF=1）— 束流入射几何组装

> 运行：`cd cases/beam_relaxed_surface && ../../venus_test.e`

## 验证目的

验证 RELAXED-SURFACE（NSURF=1）束流几何链：固定入射角 θ=30°、固定方位 χ=0、
瞄准点斜胞均匀、隐含固定碰撞参数 b=S·sinθ=4.0 Å、固定相对能闭合、表面静止——
H₂(A) + 3 原子表面 / HARMONIC 势（NT=200）。

## 原理

体系：TEST 势 HARMONIC；A 刚性 H₂（键长 0.7414 Å），B 3 原子表面（`QZB_EQ`
构成 5×5 矩形胞，SKEW=90°）。

机制链：`SURF.f`——斜胞由 BOXLX/BOXLY/SKEW 构造；`NRNDXY=1` 时瞄准点在斜胞内
均匀抽取；A 质心偏移 (S·sinθ·cosχ, S·sinθ·sinχ, S)（S=8.0 Å 恒定入射高度）；
速度方向 U=(−sinθ·cosχ, −sinθ·sinχ, −|cosθ|) 恒朝表面；`NREL=1` 时
VELA=√(2·SEREL/WTA)、VELB=0；`NZDOWN=0`：A 接束流、B 静止。全链跑
RAND0(ISEED) → 双跑逐位确定。

登记限制：HARMONIC 势无原子间表面力 → B 片段 Hessian 退化，NMODE 诊断表含 NaN
（已知 HARMONIC 伪迹，不传播到轨迹物理）；`NZDOWN=1` 在 TEST 构建不可达
（桩停机）。

## 预期与结果

- **入射几何精确**：A 高度 z=8.0（偏差 0）；θ=29.9993°；速度方向与
  (−sin30°, 0, −cos30°) 偏差 ≤1.1×10⁻⁵；b=S·sin30°=4.0 精确（4.4×10⁻¹⁶）✓
- **瞄准点均匀**：RX0/BOXLX KS D=0.0369（p=0.94）、RY0/BOXLY KS D=0.0538
  （p=0.59）✓
- **Erel 闭合与表面静止**：相对能 5.0 kcal/mol（偏差 0）；B 动量 <10⁻³ ✓

> 图：`fig_trajectory.png` —— 10 条入库轨迹的相空间（fort.1001-1010 活数据）：
> 左=原子 1 到最近原子距离 r_min(t)（NS=5 短轨迹内即 A 片段键伴距离，0.74 bohr
> 恒定——入射束离表面尚远），右=原子 1 动量模 |P|(t)。GWRITE 能量打印为陈旧值
> （见 morse_bootstrap 登记），本图只用 fort 相空间。
