# case_surface_window — 表面窗口端到端（POSCAR 胞 + A_LAT/BOXPAD + verlet）

> 运行：`cd cases/case_surface_window && ../../src/build/gfortran/venus.e`
> （驱动器读 cwd 的 `input_qct.txt`；跑 gfortran 树——ifx 目录门环境备忘沿
> docs/plans/V.md 验收 4，本 case 不承载 ifx 复验。）
> 重跑产物：`trajectories.txt`、`hist_e_rel.dat`、`hist_thta.dat`（全文件不入库，
> 摘录见第 4 节与 `results.txt`）。

## 1. 验证目的

表面范式全链首例：真程序在**派生**表面范式（POSCAR 含胞 → surface_model=1，
E14 规约"任一文件带胞 → 表面类"）下端到端跑绿——POSCAR 解析（含**速度块跳过**
路径，C5.6 首件变体的 case 形覆盖）、胞字段派生（box_lx/box_ly/skew）、
表面窗口键 reconcile（A_LAT/BOXPAD 表面范式必键、装配期拉取，E23）、单片段
弹丸入射（胞内投放 N_AIM=1——胞字段的动力学消费面）、谐振容器传播与统计面。
与 check_input (o) 族的 input→sysdef→roster 表面链判据同源——判据在 build/ 下
一次性构造并随判据程序跑，本 case 在 cases/ 下以可重复重跑的运行材料落档
全链版本（判据先行：本 README 第 3 节预期锚先写，后跑，第 4 节对照入档）。

统计面对偶：气相跑（case_harmonic_scatter）的 thta 直方图是退化文件、e_rel
有真值；本表面跑恰好倒过来——thta 有真值（表面散射角通道 ang(19) 单片段即可
算）、e_rel 退化为恒 0（相对运动块仅双片段计算）。两 case 合起来把两个通道
的"真值/退化"两面各走一遍。

## 2. 原理

体系与范式：

- `system/` 恰一个文件 `projectile.poscar`，承载全部体系定义（胞 + 分子片段）。
  **表面范式单片段规约**（incident 成员契约）：表面范式的 roster 恰承载一个
  片段——弹丸；表面不是 roster 片段（胞是系统级单例事实，单胞法令）。故分子
  片段搭载 POSCAR 文件本身（一个文件 = 一个反应物 = 一个片段）；在 system/ 再
  加一个 .xyz 会把片段数翻倍 → 装配期具名停机（inc_init 拒绝表面范式片段数
  ≠1，第 4 节 A6 探针取证）。
- 胞（Cartesian 坐标模式，标度因子 1.0）：
  `a1 = 2·A_LAT·x̂ = (6,0,0)`、`a2 = 2·A_LAT·ŷ = (0,6,0)`、`a3 = (0,0,20)`
  （真空垫层——roster 只取 |a1|、|a2| 与 a1·a2 夹角，a3 长度不入任何字段）。
  A_LAT = 3.0（输入键）→ roster 派生 box_lx = box_ly = 6.0 = 2·A_LAT
  （构造式，十进制精确）、skew = π/2（sin(skew)=1 > 0，N_AIM=1 的合法性前提
  由构造保证）。
- 分子片段：H2 型双原子，2 个 H（质量 1.008 u），键长 0.7 A 沿 x 放置、质心在
  原点（±0.35, 0, 0）。
- POSCAR 带**速度块**（坐标块后每原子一行，值 +0.50/-0.50 x̂）：specio 读过
  不消费（速度是初值材料、归采样器所有）。本跑初态动量全部来自采样链，帧动量
  不含速度块值——解析干净完成（run 不停机）即该路径的覆盖证词。
- 势能唯一来自 reference 容器：各向同性谐振井 V = K/2·Σ|q|²，K = K_HARM =
  0.04 eV/A²，作用于弹丸全部 6 个笛卡尔自由度；片段无内部势。

初始条件（采样链：片段成员先写内部，incident 通道后组装）：

- DIST_SCHEME = THERMALIZE（单片段一行）：Maxwell 种子 300 K + 精确漂移移除
  （质心动量投影清零）+ 25 步 DT_EQ=0.05 平衡（走注册传播器 verlet）。
- incident 通道表面分支：E_REL = 200 kcal/mol（表面 mu_eff = 弹丸质量
  2.016 u → v = 2.881 A/(10 fs)）；N_B = 1 碰撞盘 b² 均匀于 [0, B_MAX²]、盘
  方位 φ 均匀（B_MAX = 2.0）；N_AIM = 1 **胞内投放**——瞄准点 s1·a1 + s2·a2
  （s 均匀 [0,1)，即 x0,y0 ∈ [0, 6.0)——roster 胞字段 box_lx/box_ly/skew 的
  动力学消费面）；θ=0 顺 z 固定束（N_THTA=0 默认，无横向后撤）；R_SEP = 8 A
  投放高（弹丸质心落在 z = 8.0 恰）；表面取反冲——弹丸质心速度 = +v·u_inc。
  ISEED = 20260916 定死随机流。

动力学与终止：

- verlet（velocity-Verlet，INTEGRATOR 默认），DT = 0.05（时间单位 10 fs）。
- Q_TERM = 9.5 A：任一坐标越界即终止（出射完成）；MAX_STEPS = 400 安全网。
- **无 E_SPLIT**（reference_pars.txt 不含该键）：分类槽不绑——单片段跑的
  相对运动量 e_rel 恒 0（fin_eval 相对运动块仅双片段计算），E_SPLIT 判决会
  退化成常数 'capture'；合法无分类跑（E19 卫兵跳过路径）是本 case 的既定
  姿态，统计面打 `outcome counts: none`。
- 几何图景（推导用）：弹丸质心自 z=8 以 v=2.881 下行；井中每原子简谐
  ω = sqrt(conv·K/m_H) ≈ 0.1957 /(10 fs)（换算 1 eV/u = 0.9648 A²/(10 fs)²
  乘于分子——同 case_harmonic_scatter 第 2 节；方向勘误见第 4 节），z 振幅
  A_z = sqrt(8² + (v/ω)²) ≈ 16.8 →
  越过原点后在 |z| > 9.5 出射，约 1.1e2 步（远早于 400 步地平线）；横向
  x0,y0 ∈ [0, 8] 各自小幅振荡（振幅 ≤ ~8.4 < 9.5，不触终止）。出射时
  v_z ≈ -2.7 主导、横向分量 ∝ ω·(x0,y0)·sin(ωt_e) → 表面散射角
  thta = acos(v_z/|v|) ∈ [142°, 180°]（下行半球，v_z < 0），随投放点/盘位
  散布约 38°——thta 直方图非退化的推导基础。

## 3. 预期锚（先写后跑——判据先行）

**A1 装配链与统计面行**：退出码 0；assembly 行逐字
`assembly complete: trajectories=8 atoms=2 fragments=1 paradigm=1 members=11 init=ok columns=3`
（paradigm=1 = POSCAR 胞派生表面范式；fragments=1 = 单文件单片段；两窗口键
缺一即 reconcile 具名停机——**run 走完装配本身即证 A_LAT/BOXPAD 拉取通过**）；
statistics 行 `statistics: trajectories=8 valid=8`；结局行
`outcome counts: none (no classification bound - skipped)`（无分类合法跑）；
drift 行与直方图行存在（A2/A4 各自锚定）。

**A2 能量守恒（帧级）**：trajectories.txt 中每帧 e_tot 与本轨迹首帧 e_tot 之差
≤ **2.0e-1 kcal/mol**；stdout 打印行 `conservation drift: ... max|e_fin-e0|`
同界。
推导（沿 check_driver (b) / case_harmonic_scatter A1 的方法，本体系尺度重推）：
velocity-Verlet 对简谐模的能量偏移按 (ω·dt)² 比例于模能量且不长期积累（阴影
哈密顿量守恒）。本 case ω ≈ 0.1957 /(10 fs)（井 K=0.04 eV/A² 对 H；换算
方向勘误见第 4 节）、dt = 0.05 → (ω·dt)² ≈ 9.6e-5。能标 e0 = E_REL(200) +
初始井势（质心在 (x0,y0,8)，x0²+y0² ≤ (6+2)² → 两原子 Σ|q|² ≤ 2·(64+64)，
井势 ≤ K/2·256 = 5.12 eV ≈ 118 kcal/mol）+ 平衡内能（300 K 种子三内禀
自由度，~1 kcal/mol 量级）→ 包络取 ≤ 400 kcal/mol。偏移界 ≈ 9.6e-5 × 400
≈ 3.8e-2（~4e-2）kcal/mol；锚定 2.0e-1（约 5 倍裕度——(b) 判据 5e-2 对 7e-3 的同一
裕度纪律；V.1 同锚值，包络按本 case 尺度重推写明）。

**A3 POSCAR/A_LAT 几何自洽（构造式 + 可观察面）**：
- 构造式即锚：胞矢量 a1 = 2·A_LAT·x̂、a2 = 2·A_LAT·ŷ（A_LAT = 3.0 为输入
  键，胞边 6.0 十进制精确）→ roster%box_lx = box_ly = 2·A_LAT、skew = π/2。
- 可观察面（trajectories.txt 首帧）：弹丸质心 = (x0 + b·cosφ, y0 + b·sinφ,
  8.0)——z 恰等于 R_SEP（com_shift 精确平移）；x、y ∈ [-2.0, 8.0]（胞内投放
  [0,6) + 碰撞盘 ±B_MAX）；8 轨迹逐条成立。
- BOXPAD 注记：BOXPAD = 2.0 装配期落地 roster%box_pad（reconcile 拉取，A1
  的完成即证）；其下游消费面 = incident 通道 N_AIM=2 模式（簇盒投放
  [-box_pad, box_pad]²）——本 case 未行使（N_AIM=1 行使胞字段面），不造
  消费面。

**A4 thta 非退化（与 case_harmonic_scatter A4 对偶）+ e_rel 退化说明**：
- hist_thta.dat 非退化：≥ 3 个非零 bin 且无任何 bin 计数 = 8；全部 thta 真值
  ∈ (90°, 180°]（下行半球——出射 v_z < 0 的推导面）；40 bin 行 + 1 行
  `TOTAL & P`；计数和 8、P = 1.0。
- 独立复核：末帧 sec_qp 动量 → 弹丸质心速度 v_c = Σp/M →
  thta' = acos(v_z/|v_c|)·(180/π) 逐轨迹重算，与直方图一致（以重算值按
  [0, max+1)·40 bin 规则独立重分 bin，逐 bin 计数一致）。
- hist_e_rel.dat 退化说明（预期内，非缺陷——与气相跑 thta 退化文件对偶）：
  fin_eval 相对运动块仅双片段计算 → 表面单片段 e_rel 逐位 0 → 值域
  [0, max+1) = [0,1) 全落 bin 1（bin-1 计数恰 8、其余 0、间距 0.025 = 1/40、
  P = 1.0）。与 C8 停放区"直方图产出按 surface_model 门控候选"同源登记，
  与本 case 判定无关。

**A5 确定性**：同输入双跑，trajectories.txt 逐字节一致（rng 单流 + 固定算术，
双跑差异即缺陷）。

**A6 组成法令探针（负锚，跑后取证）**：在 system/ 临时加入第二个文件
（frag.xyz，双原子弹丸几何）→ 装配期具名停机：stderr 含 `inc_init` 与
"exactly one projectile fragment"；退出码 ≠ 0。证表面范式 roster 组成法令
——分子片段搭载 POSCAR（单文件）的设计依据。

## 4. 结果对照（跑后入档，2026-09-15，gfortran 树 venus.e）

stdout 摘录（全文见 results.txt）：

```
assembly complete: trajectories=8 atoms=2 fragments=1 paradigm=1 members=11 init=ok columns=3
statistics: trajectories=8 valid=8
outcome counts: none (no classification bound - skipped)
conservation drift: trajectories=8 max|e_fin-e0|=  7.4535E-04 kcal/mol
histograms: rows=8 files=hist_e_rel.dat hist_thta.dat
```

| 锚 | 结果 | 判定 |
|---|---|---|
| A1 装配链/统计面 | 退出码 0；assembly 行逐字命中（paradigm=1、fragments=1、两窗口键过 reconcile）；valid=8；outcome counts: none（无分类合法跑）；drift/直方图行在 | 绿 |
| A2 帧级守恒 ≤ 2.0e-1 | 8 轨迹最坏帧偏移 1.9419e-03 kcal/mol（打印行 7.4535E-04，终帧族）；e0 ∈ 207.8–265.9（包络 ≤ 400 内）；帧数 113–129（< 401）、末帧 t 5.60–6.40（< 20，Q_TERM 中途终止）；traj/step/t 单调 | 绿（低于锚两个量级，与推导界 ~4e-2 同量级） |
| A3 POSCAR/A_LAT 自洽 | 构造式 6.0 = 2×3.0 十进制精确；首帧弹丸质心 z = 8.0000000000（恰 R_SEP）×8；x ∈ [1.64, 4.79]、y ∈ [-0.07, 6.29] 全落 [-2, 8] 包络（胞内投放 [0,6) + 盘 ±2）；BOXPAD 落地注记按第 3 节 | 绿 |
| A4 thta 非退化 | thta 真值（末帧 sec_qp 动量独立重算 v_c = Σp/M → acos(v_z/|v_c|)）146.67–171.99°，全落下行半球 (90°, 180°]；hist_thta.dat 非零 bin 34/38/39/40（计数 1/2/3/2，4 个非零 bin、最大 bin 计数 3 < 8），与独立重分 bin **逐 bin 一致**；40 bin 行 + TOTAL & P、计数和 8、P = 1.0、间距 4.324647 = (max+1)/40、首中心 = 半间距；hist_e_rel.dat 退化如预期（bin 1 计数恰 8、其余 0、间距恰 0.025、P = 1.0）——表面单片段 e_rel 恒 0，与 case_harmonic_scatter A4 的 thta 退化文件对偶 | 绿（thta 非退化 + e_rel 退化均为预期形态） |
| A5 确定性 | 同输入双跑 trajectories.txt / hist_thta.dat 逐字节一致（cmp）；stdout 除 timing 行外一致；stderr 均空 | 绿 |
| A6 组成法令探针 | system/ 临时加 frag.xyz → exit 1，stderr：`inc_init: a surface-paradigm incident channel carries exactly one projectile fragment (the surface is not a roster fragment); fragment count =2` + loud-stop 行；移除后复跑绿 | 绿（负锚命中——分子片段搭载 POSCAR 的设计依据实证） |

跑后认识修正（a priori 描述 vs 实测，诚实入档；两处均不触及任何锚的成立）：

1. 第 2 节几何图景预估"约 1.1e2 步"：实测 113–129 帧（含热扰动展宽），
   量级吻合，无修正必要；出射时刻随轨迹不同（热噪声使 z 出射步散布），
   与 V.1 的帧数散布同源。
2. thta 推导区间 [142°, 180°]：实测 [146.67°, 171.99°] 落在推导区间内
   （上界 180° 需投放点恰在原点，8 次抽样未及；非锚项，无影响）。

勘误（V.5 收口批，2026-09-16——V.3 审查移交；锚值与全部判定不受影响）：
第 2/3 节的 ω 换算方向原写反——sqrt(K/(conv·m_H)) = 0.2028 应为
**sqrt(conv·K/m_H) ≈ 0.1957**（conv 乘于分子，同 V.1/case_b_layer 口径）。
推导链随之更正：(ω·dt)² 1.03e-4 → 9.6e-5（入界窄 7%）、A_z 16.3 → 16.8、
井势包络子项 ~1.8e2 → **~118**（K/2·2·(64+64) = 5.12 eV 复算，原值偏松）、
偏移界 ~4e-2 → 3.8e-2；能标包络 ≤ 400 与锚 2.0e-1 均不动
（200 + 118 + ~1 = 319 ≤ 400，裕度反而放宽）。

