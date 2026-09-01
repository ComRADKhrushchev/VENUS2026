# 表面 MD 恒温均衡 — 恒温板块采样的温度分布

> 运行：`cd cases/md_surface_equilibrate && ../../venus_test.e`（恒温
> 分支经 GASDEV 调用 Fortran 内建 random_number 系统种子，双跑非逐位，
> 验收以统计量为准）。

## 1. 目的

本 case 验证 MD 恒温表面采样（INIT_SAMPLING_B=MD，即 NACTB=7）在
FULL-SURFACE 板块上的统计正确性：恒温均衡结束时 B 片段温度分布应
匹配理论形式 300·χ²(9)/9（3 个 B 原子、3N=9 个自由度）。判据来源
为统计力学：3N 个独立高斯动量平方和服从 χ² 分布。

## 2. 理论与体系

体系为 4 原子（1 A + 3 B，各 1.008 u）FULL-SURFACE TEST 势
HARMONIC：每原子独立简谐阱。势参数（源码
src_TEST/test_potentials.f90:25-26）：

| 参数 | 值 | 单位 | 出处 |
|---|---|---|---|
| k_x = k_y = k_z | 1 | eV/bohr² | HARM_K，test_potentials.f90:25 |
| 平衡位置 x0 | 0, 0, 0 | bohr | HARM_X0，test_potentials.f90:26 |
| THERMOTEMP | 300 | K | input_qct.txt |
| NSCALE / NEQUAL | 200 / 0 | 步 | input_qct.txt |
| NT / NS | 100 / 3 | 条 / 步 | input_qct.txt |

机制链（SELECT.f90 md_equilibrate_b）：A 坐标存 QTEMP 并抬升
100 bohr 冻结；B 动量按 Maxwell-Boltzmann 采样
P=GASDEV·DESKET·√W，DESKET=√(0.00198717·T·C1)，其中
C1=0.04184（kcal/mol→code，venus_params.f90:69）；去总动量后
速度重标（THERMO(0,1)），每 100 步 Andersen 重采样——保存前最后
一步操作是 Andersen 重采样，故均衡末端 B 动量为标准正态 MB 样本，
温度 T~300·χ²(9)/9：均值 300 K、标准差 300·√(2/9)=141.4 K。

## 3. 方法与流程

1. setup_b_coords 将 QZB_EQ 板块构型置入 Q（B 居于原点阱底）。
2. md_equilibrate_b：A 抬升冻结、B 动量 MB 采样 + 恒温均衡
   （stdout `DOING MD SAMPLING` 与 `EQUALIBRATION ... IS NOW OVER`
   横幅，见 results.txt）。
3. restore_a_positions 恢复 A 至真实位置（z=11.0 bohr，无 +100
   残留）。
4. GWRITE.f90:51-63 每动力学步打印 B 片段 `system temperature=`
   （stdout，每轨迹 4 段，首段即均衡末端温度）；fort.8/fort.1001
   入库轨迹，fort.52 记录散射布局。

## 4. 核心验证：温度分布命中

判据：100 条轨迹均衡末端温度经验分布与 300·χ²(9)/9 的 KS 统计量
Dmax≤0.3（=3/√N）。实测（NT=100 stdout 首段温度统计，
results.txt）：KS D=0.084（p=0.45），接受分布假设。结论：恒温
采样统计正确。佐证：均值 288.0 K（理论 300，−4%，门限 [255,345]）、
标准差 119.6 K（理论 141.4，−15%，门限 [106,177]）；A 片段均衡期
固定于抬升位置、结束恢复至 z=11.0 bohr（fort.1001 step-0）。

> 图：`fig_equilibration.png` —— 温度采样（stdout `system
> temperature` 行，100 轨迹 ×4 打印）：左 = 各轨迹 step-0 温度散点
>（79-639 K、均值 288.0，点线 = 恒温目标 300 K），右 = 全部温度
> 点云 vs 目标线。
