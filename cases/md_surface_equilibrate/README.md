# 表面 MD 恒温均衡（INIT_SAMPLING_B=MD）— 恒温表面采样

> 运行：`cd cases/md_surface_equilibrate && ../../venus_test.e`
> 注：恒温分支经 GASDEV 调用 Fortran 内建 `random_number`（系统种子），双跑非逐位；
> 验收以统计量为准。

## 验证目的

验证 MD 恒温表面采样（`INIT_SAMPLING_B=MD`，NACTB=7）在 FULL-SURFACE 板块上的统计
正确性（NT=100，`THERMOTEMP=300 K`）：均衡结束时 B 片段温度分布应匹配理论形式
300·χ²(9)/9（3 个 B 原子，3N=9 个自由度）。

## 原理

体系：TEST 势 HARMONIC（每原子独立简谐阱，k=1 eV/Å²，平衡位置原点）；FULL-SURFACE
预设下 A 为单原子、B 为 3 原子板块。

机制链：`md_equilibrate_b`——A 坐标提升 100 bohr 冻结；B 动量 Maxwell-Boltzmann
采样（P=GASDEV·DESKET·√W，DESKET=√(k_B·T·C1)，C1=0.04184），去总动量；速度重标
（NSCALE 阶段）后每 100 步 Andersen 重采样——均衡末端 B 动量为标准正态 MB 样本。
温度 ∝ ΣP²（3N 个独立高斯动量的平方和），故 T ~ 300·χ²(9)/9：理论均值 300 K、
标准差 300·√(2/9)=141.4 K。均衡结束 A 恢复到 QZA_EQ（z=11.0 bohr）。

## 预期与结果

- **温度分布命中**：100 条轨迹均衡末端温度 vs χ²(9)——KS D=0.084（p=0.45），
  接受分布假设 ✓
- **矩匹配**：均值 288.0 K（理论 300，偏差 −4%，门限 ±15%）；标准差 119.6 K
  （理论 141.4，偏差 −15%，门限 ±25%）✓
- **A 冻结-恢复语义**：均衡期 A 固定于提升位置，结束恢复至 z=11.0 bohr ✓

> 图：`fig_equilibration.png` —— 温度采样（stdout `system temperature` 行，
> 活数据，100 轨迹 ×4 打印）：左=各轨迹 step-0 温度散点（79-639 K、均值 288.0，
> 点线=恒温目标 300 K——χ²(9) 形态的散布），右=全部温度点云 vs 目标线。
