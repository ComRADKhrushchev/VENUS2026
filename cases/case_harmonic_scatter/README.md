# case_harmonic_scatter — 谐振容器双片段散射（verlet 基准）

> 运行：`cd cases/case_harmonic_scatter && ../../src/build/gfortran/venus.e`
> （驱动器读 cwd 的 `input_qct.txt`；跑 gfortran 树——ifx 目录门环境备忘沿
> docs/plans/V.md 验收 4，本 case 不承载 ifx 复验。）
> 重跑产物：`trajectories.txt`、`hist_e_rel.dat`、`hist_thta.dat`（全文件不入库，
> 摘录见第 4 节与 `results.txt`）。

## 1. 验证目的

新树 case 集的守恒/确定性/统计面基准：真程序（input → 组装 → 采样 → 传播 →
记录 → 统计的全四相循环）在谐振井容器上端到端跑绿。与
scripts/check_driver.f90 的 (b)/(det)/(f)/(l) 判据同源——判据在 build/ 下一次性
构造并随判据程序跑，本 case 在 cases/ 下以可重复重跑的运行材料落档同一套
推理（判据先行：本 README 第 3 节预期锚先写，后跑，第 4 节对照入档）。

## 2. 原理

体系与容器：
- 两个 H2 型双原子片段：`system/frag_a.xyz`（弹丸，键长 0.7 A）、
  `system/frag_b.xyz`（靶，键长 0.8 A）；全 H，质量 1.008 u。
- 势能唯一来自拖入的 reference 容器：各向同性谐振井
  V = K/2 · Σ|q_i|²，K = K_HARM = 0.04 eV/A²（`reference_pars.txt`），
  作用于全部 12 个笛卡尔自由度；片段之间无相互作用——体系即 12 个同频
  简谐模，每原子频率 ω = √(K/m_H) ≈ 0.196 /(10 fs)
  （换算 1 eV/u = 0.9648 A²/(10 fs)²）。

初始条件（采样链，片段内成员先写内部、incident 通道后组装对）：
- DIST_SCHEME = THERMALIZE,ROTATION：片段 1（弹丸）Maxwell 种子 300 K +
  25 步 DT_EQ=0.05 平衡；片段 2（靶）转动 J = 3。
- incident 通道：E_REL = 200 kcal/mol（固定碰撞能）；N_B = 1 碰撞盘
  b² 均匀于 [0, B_MAX²]；R_SEP = 8 A 定点投放（靶质心在原点、弹丸质心在
  (b, 0, 8)）；θ=0 顺 z 束。ISEED = 20260915 定死随机流。

动力学与终止：
- verlet（velocity-Verlet），DT = 0.05（时间单位 10 fs，即 0.5 fs/步）。
- Q_TERM = 9.5 A：任一坐标越界即终止（散射完成）；MAX_STEPS = 400 安全网。
- 几何图景（推导用）：初始井势 ≈ K/2·2·(b²+64)·换算 ≈ 66–68 kcal/mol；
  每原子总能 ~3.5 eV → 振幅 10–13 A，两片段的原子都能越过 |q|=9.5，
  终止在 ~10² 步（远早于 400 步地平线 t=20）。
- 分类：E_SPLIT = 1.0 kcal/mol，终态相对径向能量 e_rel 高于此 → scattering。

## 3. 预期锚（先写后跑——判据先行）

**A1 能量守恒（帧级）**：trajectories.txt 中每帧 e_tot 与本轨迹首帧 e_tot
之差 ≤ **2.0e-1 kcal/mol**；stdout 打印行 `conservation drift: ... max|e_fin-e0|`
同界。
推导（沿 check_driver (b) 的方法，以 case 尺度重推）：velocity-Verlet 对简谐模
的能量偏移按 (ω·dt)² 比例于模能量且不长期积累（阴影哈密顿量守恒）。本 case
ω ≈ 0.196 /(10 fs)、dt = 0.05 → (ω·dt)² ≈ 9.6e-5。能标
e0 = E_REL(200) + 初始井势(≈67) + 内能（300 K 平衡 + J=3 转动，≲10）
≈ 270 kcal/mol，推导包络取 ≤ 300。偏移界 ≈ 9.6e-5 × 300 ≈ 3e-2 kcal/mol；
锚定 2.0e-1（约 7 倍裕度——(b) 判据 5e-2 对 7e-3 的同一裕度纪律移到 case 尺度，
能标放大约 4 倍、锚随之放宽 4 倍）。

**A2 确定性**：同输入双跑，trajectories.txt 逐字节一致（`fc /b`；rng 单流 +
固定算术，双跑差异即缺陷）。

**A3 统计面**（8 条轨迹全部完成、结局表/直方图产出且格式核对）：
- 退出码 0；assembly 行
  `assembly complete: trajectories=8 atoms=4 fragments=2 paradigm=0 members=11 init=ok columns=3`；
- statistics 行 `statistics: trajectories=8 valid=8`；
- 结局表 `outcome counts: scattering=8`。推导：终止时刻（首个 |q|>9.5）井势
  ≈ 1.2e2 kcal/mol 量级，e0 ≈ 268 → 平动动能 ≳ 1e2 kcal/mol 且以相对径向分量
  为主 → e_rel ~ 1e2 kcal/mol，高出 E_SPLIT = 1.0 约两个量级（弱阱 K=0.04
  基本不偏折飞越对；裕度两个量级，非刀锋判断）→ 8 条全 scattering；
- 中途终止：每轨迹帧数 < MAX_STEPS+1 = 401 且末帧 t < 20（10 fs 单位）；
- 直方图文件 `hist_e_rel.dat` / `hist_thta.dat`：各 40 个 bin 行 + 1 行
  `TOTAL & P`；bin 中心等间距、首中心 = 半个间距；e_rel 侧计数和 = 8、
  P = 8/8 = 1.0。

**A4 thta 退化文件说明（C8 移交的用户可见注意）**：气相跑没有表面对，
终态散射角通道 thta 逐位为 0 属预期——`hist_thta.dat` 是退化文件：8 行
全部落在 bin 1（bin-1 计数 8、其余 0），间距恰 0.025 = 1.0/40（值域规则
[0, max+1) = [0,1) 均分 40 bin）。该文件的存在与形状均为预期行为，非缺陷；
是否按表面范式门控产出属候选处置（C8 停放区登记），与本 case 判定无关。

## 4. 结果对照（跑后入档，2026-09-15，gfortran 树 venus.e）

stdout 摘录（全文见 results.txt）：

```
assembly complete: trajectories=8 atoms=4 fragments=2 paradigm=0 members=11 init=ok columns=3
statistics: trajectories=8 valid=8
outcome counts: scattering=8
conservation drift: trajectories=8 max|e_fin-e0|=  1.2648E-03 kcal/mol
histograms: rows=8 files=hist_e_rel.dat hist_thta.dat
```

| 锚 | 结果 | 判定 |
|---|---|---|
| A1 帧级守恒 ≤ 2.0e-1 | 8 轨迹最坏帧偏移 1.2648e-03 kcal/mol（打印行同值）；每轨迹 e0 在 241.0–263.6 kcal/mol（推导包络 ≤ 300 内）；帧时刻按 0.05 步进、traj/step/t 单调 | 绿（低于锚两个量级，与推导界 ~3e-2 同量级） |
| A2 确定性 | 同输入双跑 trajectories.txt 逐字节一致（cmp）；stdout 除 timing 行外亦一致 | 绿 |
| A3 统计面 | 退出码 0；assembly 行逐字命中；trajectories=8 valid=8；outcome counts: scattering=8；每轨迹帧数 98–117（< 401）、末帧 t 4.85–5.80（< 20，Q_TERM 中途终止）；hist_e_rel.dat 40 bin 行 + TOTAL & P 行、计数和 8、P = 1.0、bin 中心等距（打印精度内最大偏差 1e-6）、首中心 = 半间距；非零 bin 为 32/35/36/40（各计 2，e_rel ≈ 113–145 kcal/mol——高出 E_SPLIT = 1.0 两个量级，与推导一致） | 绿 |
| A4 thta 退化 | hist_thta.dat 40 bin 行 + TOTAL & P 行；bin 1 计数恰 8、其余 bin 全 0；间距恰 0.025 = 1.0/40；P = 1.0——与预期逐项一致 | 绿（预期内退化文件） |

