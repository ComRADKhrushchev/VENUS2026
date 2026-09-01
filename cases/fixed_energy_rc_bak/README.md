# fixed_energy_rc — 固定能量含反应坐标的微正则采样

> 运行：`cd cases/fixed_energy_rc && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

本 case 验证固定能量含反应坐标的微正则采样（`INIT_SAMPLING_A=FIXED-ENERGY`，`NBAR=3`
自动触发 NACTA=6）在等边 LEPS H₃ 上的正确性（NT=1000，ENMT_A=20.0 kcal/mol）：
反应坐标（RC）能量的微正则态分配与总能量闭合。

## 2. 理论与体系

体系为 TEST 势 LEPS 等边 H₃ 极小（NLINA=0，无弯曲伪模；QZA_EQ 边长 1.759910 Å，
取自 input_qct.txt）。LEPS 参数单位为 eV/bohr/bohr⁻¹（`src_TEST/test_potentials.f90`）。
机制链：`NBAR=3` → `SELECT.f90` 强制 NACTA=6 → Beyer-Swinehart 网格态计数
（`SELECT.f90:749-801`，STEP=0.1·ν₉ + idnint 舍入）→ RC 能量按微正则态分配 →
`INITQP` 组装。总能量 ENMT_A 分解为振动能 EVIBA 与 RC 能，能量闭合缩放保证
EVIBA+RC 守恒；kcal/mol→a.u. 换算 C1=0.04184（`venus_params.f90:60-80`）。

核心科学事实：程序实现的是 Beyer-Swinehart 网格态计数而非理想均匀分布，观测精确
匹配网格 PMF（联合 χ²=6.2/12 df；若按均匀假设 χ²=36.4 被拒绝）；idnint 舍入是固有的
文档化近似。共线构型存在 0.5134 cm⁻¹ 弯曲伪模会使 STEP=0.05134 → NSTEP 溢出停机
（登记缺陷 F22，等边构型规避）；本链此前因 ENMT_A 读入分支漏 NACTA=6（缺陷 F21，
已一行修复于 `venus_input.f90`）不可用。

## 3. 方法与流程

1. 程序读取输入并回显 FIXED ENERGY SAMPLING INCLUDING RC 与每轨迹 RC 能量
   （stdout，results.txt；第 1 条 E_RC=2.7598 kcal/mol）。
2. 初始化随机数（ISEED=20260829），对 NT=1000 条轨迹做 Beyer-Swinehart 网格态计数，
   抽取 13 个微正则态与 RC 态。
3. `INITQP` 组装坐标与动量，能量闭合缩放（stdout `INTERNAL ENERGY` 回显 EVIBA=ENMT_A−E_RC）。
4. 每轨迹传播 NS=5 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

**判据**：总能量应闭合——EVIBA+RC 相对 ENMT_A 的偏差应低于 1%；微正则态分配应匹配
网格 PMF 而非均匀分布。**实测**（results.txt 直方图 n8 模 0:278 1:216 2:217 3:139 4:74 5:76、
n9 模 0:482 1:321 2:131 3:66）：1000 条轨迹能量闭合最大偏差 0.00094；13 个微正则态 KS
（n8 Dmax=0.0077、n9 Dmax=0.0177，≤ 3/√N=0.0949），联合 χ²=6.2（df=12，门限 24.05）。
**结论**：网格态计数与 RC 分配机制按设计工作。佐证：RC 态匹配 1000/1000，
RC 能量 ∈ (0.500, 15.972) kcal/mol ⊂ (0, ENMT_A)。

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量（stdout `INTERNAL ENERGY`
> 回显，活数据）：左为各轨迹 EVIBA（9.2-18.5 kcal/mol 波动，E_RC 逐轨迹微正则重分配的
> 直接体现，EVIBA+E_RC 闭合 20.0），右为逐轨迹缩放收敛过程。
