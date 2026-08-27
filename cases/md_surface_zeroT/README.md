# 表面 MD 零温分支 — 确定性路径

> 运行：`cd cases/md_surface_zeroT && ../../venus_test.e`（零温分支
> 不调用 GASDEV，双跑逐位一致）。

## 1. 目的

本 case 定案 NACT=7（INIT_SAMPLING_B=MD）的零温分支语义：
THERMOTEMP≤0.01 时 B 片段动量置零并早期返回，该路径不调用
GASDEV、不消耗系统随机数，是 NACT=7 链中唯一双跑逐位确定的锚；
同时验证 A 片段不经抬升冻结、保持自由飞行。判据来源为零温极限
的物理要求（无热运动）与确定性判据（同输入同输出）。

## 2. 理论与体系

体系为 4 原子（1 A + 3 B，各 1.008 u）FULL-SURFACE TEST 势
HARMONIC（每原子独立简谐阱），与 md_surface_equilibrate 同构。
势参数（源码 src_TEST/test_potentials.f90:25-26）：

| 参数 | 值 | 单位 | 出处 |
|---|---|---|---|
| k_x = k_y = k_z | 1 | eV/bohr² | HARM_K，test_potentials.f90:25 |
| 平衡位置 x0 | 0, 0, 0 | bohr | HARM_X0，test_potentials.f90:26 |
| THERMOTEMP | 0.0 | K（≤0.01 触发零温分支） | input_qct.txt |
| NT / NS | 1 / 5 | 条 / 步 | input_qct.txt |

机制链（SELECT.f90 md_equilibrate_b 入口守卫）：THERMOTEMP≤0.01
时仅将 B 各原子动量 P 置零后 return——无抬升冻结（QTEMP+100 逻辑
仅在恒温分支）、无 `EQUALIBRATION ... IS NOW OVER` 横幅、恒温浴
THERMO 不执行。B 冻结于阱底、A 从 MB 初值自由飞行；能量换算常数
C1=0.04184（kcal/mol→code，venus_params.f90:69），此处仅经 GWRITE
温度回读间接出现。

## 3. 方法与流程

1. setup_b_coords 将 QZB_EQ 板块构型置入 Q；A 按 MB 采样初值
   （fort.8 首段坐标动量）。
2. md_equilibrate_b 零温守卫命中：B 动量置零、早期返回（stdout
   `THERMOTEMP = 0.0` 回显与全零 `system temperature=`，见
   results.txt）。
3. A 自由飞行 NS=5 步；GWRITE.f90:51-63 每步打印 B 温度（全 0.0）。
4. fort.1001 入库逐步相空间（H 恒 461.80825 eV），fort.999 汇总
   单轨迹统计。

## 4. 核心验证：B 冻结且路径逐位确定

判据：B 原子坐标与动量全步为零（温度回读全 0.0），且双跑全部
产物逐位一致。实测（fort.1001 六个输出步 + stdout，results.txt）：
B 三原子坐标与动量逐字段 0.00e+00，`system temperature=` 六次全为
0.0；双跑仅时间戳行差异。结论：零温分支语义与确定性均通过。佐证：
A 的 z 跨步变化 0.01316 bohr（fort.1001 step-0 的 11.00000 与
step-5 的 10.98684 之差，非冻结自由飞行，z0=11.0000 bohr）。

> 图：`fig_trajectory.png` —— 入库轨迹相空间（fort.1001 活数据）：
> 左 = A 到最近表面原子距离 r_min(t)（零温下 B 冻结，A 自 11.194
> bohr 起自由下落），右 = |P|(t)（下落加速增长）。
