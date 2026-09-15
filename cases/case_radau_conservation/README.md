# case_radau_conservation — 谐振容器自适应守恒（RADAU 基准）

> 运行：`cd cases/case_radau_conservation && ../../src/build/gfortran/venus.e`
> （驱动器读 cwd 的 `input_qct.txt`；跑 gfortran 树——ifx 目录门环境备忘沿
> docs/plans/V.md 验收 4，本 case 不承载 ifx 复验。）
> 重跑产物：`trajectories.txt`、`hist_e_rel.dat`、`hist_thta.dat`（全文件不入库，
> 摘录见第 4 节与 `results.txt`）。

## 1. 验证目的

同一体系在自适应积分器下的守恒基准：与 case_harmonic_scatter 完全同一
体系、同一容器、同一随机种子，只把积分器换成 INTEGRATOR = RADAU
（ORDER = 8）——自适应成员（RA15 Gauss-Radau 序列：每步一条自适应序列、
boot 每轨迹重置、B 系数跨序列携带）在真程序全四相循环里的守恒面以 case
形落档。与 scripts/check_driver.f90 的 (r) 判据同源——判据在 build/ 下随
判据程序跑，本 case 在 cases/ 下可重复重跑。

## 2. 原理

体系、容器、初始条件与 case_harmonic_scatter 逐字相同（见其 README 第 2 节；
两 case 的 system/ 与 reference_pars.txt 内容一致，input 只差 INTEGRATOR/ORDER
与 TITLE）。同 ISEED → 同一随机流（Maxwell 种子等采样抽取逐位相同），但
thermalize 的平衡运行走注册传播器成员（内部把步长换成 DT_EQ 跑 N_EQ 步/
序列）——radau 下平衡是自适应序列、覆盖的物理时间远长于 verlet 的定点步，
平衡后状态不同。两 case 的守恒锚因此各自独立锚定（同系综对照，非同初态
对照）；初态一致性的实情见第 4 节修正记录。

动力学：
- RADAU：Everhart RA15，8 点左 Gauss-Radau 配置；一次 prop_step = 一条
  自适应序列，序列长从 DT = 0.05（10 fs 单位）播种；局部误差控制在
  ss = 10^(-ORDER) = 1e-8；DT_MIN/DT_MAX 不设（0 = 不钳制，合法默认）。
- 终止与分类同 harmonic case：Q_TERM = 9.5 A 中途终止、E_SPLIT = 1.0
  → scattering。

## 3. 预期锚（先写后跑——判据先行）

**R1 自适应守恒（帧级）**：trajectories.txt 中每帧 e_tot 与本轨迹首帧 e_tot
之差 ≤ **1.0e-6 kcal/mol**；stdout 打印行 `conservation drift: ... max|e_fin-e0|`
同界。
推导（沿 check_driver (r) 的推导，以 case 尺度重推）：RA15 把每序列局部
误差控制在 ss = 1e-8；单元级长程锚（Morse 阱、order 8、数百序列）实测
max|E-E0| ~ 1e-9 kcal/mol、判据窗口 1e-7。本 case 每轨迹 ~1e2 条序列
（终止 ~10² 步）且在更光滑的谐振链上；相对 (r) 判据的 40 序列情形，序列数
约 3 倍、能标约 4 倍（控制器为绝对误差界，能标放大不稀释 ss，但保守计入）。
锚定 1.0e-6 kcal/mol = 单元窗口 10 倍，比 verlet case 的锚 2.0e-1 严 5 个
量级——自适应误差控制，而非步长运气。

**R2 完成面与分类一致**：
- 退出码 0；assembly 行含 `trajectories=8 atoms=4 fragments=2 paradigm=0
  members=11 init=ok columns=3`；statistics 行 `statistics: trajectories=8
  valid=8`；
- 结局表 `outcome counts: scattering=8`——分类判断的是终态物理量 e_rel
  （高出 E_SPLIT 两个量级，推导同 harmonic case A3），与积分器无关；
- 每轨迹帧数 < 401 且末帧 t < 20（Q_TERM 中途终止，同一物理）。

**R3 步长自适应行为观察通道**：步径路径无任何打印（propagator 立法：
no debug prints live in the stepping path），但帧头时刻 t 是步长信息的观察
通道——radau 一次 prop_step = 一条自适应序列，sta%t 按接受序列长推进
（DT 只播种首序列），故相邻帧的 t 之差 = 各条接受序列的长度。预期形态：
首差 = DT = 0.05，此后控制器在误差远低于 ss 的光滑段按稳定比率放大序列长
（几何增长），直至终止；帧数因此远少于同物理时长下定步长积分器的帧数。
（a priori 版本曾误判帧时刻按固定 DT 推进而省略本通道——跑后修正，见
第 4 节；R1/R2 数值锚不受影响。）

## 4. 结果对照（跑后入档，2026-09-15，gfortran 树 venus.e）

stdout 摘录（全文见 results.txt）：

```
assembly complete: trajectories=8 atoms=4 fragments=2 paradigm=0 members=11 init=ok columns=3
statistics: trajectories=8 valid=8
outcome counts: scattering=8
conservation drift: trajectories=8 max|e_fin-e0|=  5.6843E-14 kcal/mol
histograms: rows=8 files=hist_e_rel.dat hist_thta.dat
```

| 锚 | 结果 | 判定 |
|---|---|---|
| R1 帧级自适应守恒 ≤ 1.0e-6 | 8 轨迹最坏帧偏移 1.1369e-13 kcal/mol（打印行 e_fin 口径 5.6843E-14，traj 1；帧级最坏在 traj 8）；e0 范围 184.7–275.7 kcal/mol | 绿（低于锚约 7 个量级——自适应控制在谐振链上近机器精度） |
| R2 完成面与分类 | 退出码 0；assembly 行逐字命中；trajectories=8 valid=8；scattering=8；每轨迹帧数 12–13（< 401）、末帧 t 4.94–6.96（< 20，Q_TERM 中途终止） | 绿 |
| R3 步长通道 | 通道确认存在：帧 t 差序列 0.0500, 0.0700, 0.0980, 0.1372, …, 2.0248（traj 1）——首差恰为 DT 种子，此后每序列 ×1.4 几何增长，12 序列到达 t=6.96；步径路径无打印（立法不变） | 绿（通道写明） |

跑后修正记录（a priori 描述错、锚判定不受影响的两处，诚实入档）：
1. 第 2 节原声称两 case 初态逐位一致——错。thermalize 平衡走注册传播器
   成员：verlet 25 步 ×0.05 = 1.25 时间单位 vs radau 25 条自适应序列
   （按 ×1.4 增长合计约 5.6e2 时间单位），平衡后状态迥异，e0 逐轨迹不同
   （traj 1 的 e0 两 case 相近纯属巧合）；Maxwell 种子抽取本身逐位相同。
2. R3 原判"帧时刻按固定 DT 推进、通道省略"——错。DT 只播种首序列，
   sta%t 按接受序列长推进，通道即上表所写；已改写第 3 节 R3。

