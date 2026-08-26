# 表面 MD 零温分支（THERMOTEMP=0）— 确定性路径

> 运行：`cd cases/md_surface_zeroT && ../../venus_test.e`

## 验证目的

验证 NACT=7（`INIT_SAMPLING_B=MD`）的零温分支：`THERMOTEMP≤0.01` 时 B 片段动量
置零、早期返回——该路径不调用 GASDEV、不消耗系统随机种子，是 NACT=7 链中唯一
双跑逐位确定的锚；同时验证 A 片段不被提升冻结。

## 原理

体系：TEST 势 HARMONIC（同 md_surface_equilibrate；A 单原子 + B 3 原子板块）。
零温守卫（`THERMOTEMP≤0.01`）下 `md_equilibrate_b` 将 B 动量置零并早期返回：
无「EQUALIBRATION ... IS NOW OVER」横幅、`THERMO` 恒温浴不执行、A 不提升
（从 MB 初值自由飞行）。因不经 GASDEV（系统种子），全部输出逐位确定。

## 预期与结果

- **B 冻结**：B 原子坐标与动量最大模 |q,p|=0.00e+00（6 步），温度回读全 0.0 ✓
- **A 自由飞行**：A z 跨步变化 0.01316 bohr（非冻结），z0=11.0000 ✓
- **逐位确定**：双跑全部产物逐位一致（仅时间戳行差异）✓

> 图：`fig_trajectory.png` —— 入库轨迹相空间（fort.1001 活数据）：左=A 到最近
> 表面原子距离 r_min(t)——零温下 B 冻结，A 自 z=11.19 bohr 起自由下落（NS=5 短程
> 内降至 11.17），右=|P|(t)（下落加速增长）。GWRITE 能量打印为陈旧值（见
> morse_bootstrap 登记）。
