# 固定能量含反应坐标（FIXED-ENERGY, NBAR=3）— 微正则反应坐标采样

> 运行：`cd cases/fixed_energy_rc && ../../venus_test.e`

## 验证目的

验证固定能量含反应坐标的微正则采样（`INIT_SAMPLING_A=FIXED-ENERGY`，`NBAR=3`
自动触发 NACTA=6）在等边 LEPS H₃ 上的正确性（NT=1000，`ENMT_A=20.0 kcal/mol`）：
反应坐标能量的微正则态分配与总能量闭合。

## 原理

体系：TEST 势 LEPS 等边 H₃ 极小（NLINA=0）。几何选择的关键性：共线构型存在
0.5134 cm⁻¹ 弯曲伪模，会使网格步长 STEP=0.05134 → NSTEP 溢出停机——登记为缺陷
**F22**（不修，等边构型规避）。

机制链：`NBAR=3` → `SELECT.f90` 强制 NACTA=6 → Beyer-Swinehart 网格态计数
（`SELECT.f90:749-801`，STEP=0.1·ν₉ + idnint 舍入）→ 反应坐标（RC）能量按微正则
态分配 → `INITQP` 组装。总能量 ENMT_A 分解为振动能 EVIBA 与 RC 能，能量闭合缩放
保证 EVIBA+RC 守恒。

核心科学事实：程序实现的是 Beyer-Swinehart **网格态计数**而非理想均匀分布——
观测精确匹配网格 PMF（联合 χ²=6.2/12 df；若按均匀假设 χ²=36.4 被拒绝）。
idnint 舍入是固有的文档化近似，非缺陷。

（本链此前因 **F21** 缺陷不可用：ENMT_A 读入分支漏了 NACTA=6 → ENMTA=0 →
零能守卫 STOP；`venus_input.f90` 一行修复后启用。）

## 预期与结果

- **能量闭合**：1000 条轨迹 EVIBA+RC 相对 ENMT_A 最大偏差 0.00094（门限 1%）✓
- **网格态精确匹配**：13 个微正则态 KS（n8 Dmax=0.0077、n9 Dmax=0.0177，
  ≤ 0.0949），联合 χ²=6.2（df=12，门限 24.05）匹配网格 PMF；RC 态匹配
  1000/1000 ✓
- **RC 能量有界**：RC ∈ (0.500, 15.972) kcal/mol ⊂ (0, ENMT_A) ✓

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量（stdout
> `INTERNAL ENERGY` 回显即 EVIBA=ENMT_A−E_RC，活数据）：左=各轨迹 EVIBA
> （9.2-18.5 kcal/mol 波动——E_RC 逐轨迹微正则重分配的直接体现，EVIBA+E_RC
> 闭合 20.0），右=逐轨迹缩放收敛过程。
