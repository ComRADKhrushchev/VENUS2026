# fixed_energy_rc — 固定能量含反应坐标的微正则采样

> 运行：`cd cases/fixed_energy_rc && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

验证固定能量含反应坐标的微正则采样（`INIT_SAMPLING_A=FIXED-ENERGY`，`NBAR=3`
自动触发 NACTA=6）在等边 LEPS H₃ 上的正确性（NT=1000，ENMT_A=20.0 kcal/mol）：
反应坐标（RC）能量的微正则态分配与总能量闭合。

## 2. 理论与体系

流程：`NBAR=3` → `SELECT.f90` 强制 NACTA=6 → Beyer-Swinehart 网格态计数
（`SELECT.f90:749-801`，STEP=0.1·ν₉ + idnint 舍入）→ RC 能量按微正则态分配 →
`INITQP` 组装。总能量 ENMT_A 分解为振动能 EVIBA 与 RC 能，能量闭合缩放保证
EVIBA+RC 守恒。注：共线构型存在 0.5134 cm⁻¹ 弯曲低频模式会使态计数步长过小而
溢出停机（登记缺陷 F22，等边构型规避）；本链此前因 ENMT_A 读入分支漏 NACTA=6
（缺陷 F21，已一行修复于 `venus_input.f90`）不可用。

## 3. 方法与流程

1. 读取输入并打印 FIXED ENERGY SAMPLING INCLUDING RC 与每轨迹 RC 能量（stdout，
   results.txt；第 1 条 E_RC=2.7598 kcal/mol）。
2. 初始化随机数（ISEED=20260829），NT=1000 条轨迹做网格态计数，抽取 13 个微正则
   态与 RC 态；`INITQP` 组装坐标动量并能量闭合缩放（EVIBA=ENMT_A−E_RC）。
3. 每轨迹传播 NS=5 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

见 `fig_sampling_stats.png`。

> 图注：左 = 各轨迹 EVIBA（9.2-18.5 kcal/mol 波动，E_RC 逐轨迹微正则重分配的直接
> 体现，EVIBA+E_RC 闭合 20.0），右 = 逐轨迹缩放收敛过程（stdout `INTERNAL
> ENERGY`，10 条轨迹）。判定依据：1000 条轨迹能量闭合最大偏差 0.00094（门限 1%）；
> 13 个微正则态 KS（n8 Dmax=0.0077、n9 Dmax=0.0177，≤3/√N=0.0949），联合 χ²=6.2
>（df=12，门限 24.05），匹配网格 PMF 而非均匀分布。
