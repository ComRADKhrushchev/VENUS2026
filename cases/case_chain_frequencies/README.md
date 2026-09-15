# case_chain_frequencies — 链级频率真值（Hessian→频谱→采样 后验）

> 运行（两腿各在己目录，探针读 cwd 的 `input_qct.txt`/`system/`/`reference_pars.txt`，
> 与驱动器同姿态）：
> `cd cases/case_chain_frequencies/boltz && ../../../src/build/gfortran/probe_chain_freq.exe boltz`
> `cd cases/case_chain_frequencies/normalmode && ../../../src/build/gfortran/probe_chain_freq.exe normalmode`
> 探针编译（src/ 下，沿 probe_glo 配方）：
> `gfortran -O2 -Ibuild/gfortran -o build/gfortran/probe_chain_freq.exe ../scripts/probe_chain_freq.f90 $(find build/gfortran -name '*.o' ! -name 'driver.o')`
> 探针不写文件——全部证据行落 stdout（入档 `results.txt`，见第 4 节）；
> 退出码 0 = 全锚绿。判据主体双写在探针头注（scripts/probe_chain_freq.f90）与本 README
> 第 3 节，先写后跑。

## 1. 验证目的

C5-R1 的链级另一半：**Hessian→eig_sym→频谱→boltzmann/normalmode 以真频率采样**
的链级真值后验（C5 期判据为合成频率注入版，本 case 换真容器真力导出的频率表）；
同时是 dens_states/eig_out 的前置生产者演示（C8 移交项的兑现）。附带任务：
**查明装配期频率导出通道的接线状态并如实定性**（第 2 节）。

## 2. 原理

**通道查明（先行，结论：机制存在但程序内未接线 = 停放，非缺陷）**，四条证据：
- reg_dist_boltzmann / reg_dist_normalmode 头注与代码逐字记录：频率表
  w_mode_a / w_mode+c_mode 是 "ASSEMBLY DATA, not an input key - **no production
  source stages it yet**"——两 seam 只暂存 T_VIB_A/T_VIB_B 与 E_VIB 三个标量键，
  表本身无任何生产者暂存；选中成员而无表 → 成员自己的 fail-loud 停机
  （check_dist_boltzmann (f) 已证零谱停机路径）。
- src/ 全树 `hessian_cd` 零生产调用者（仅 scripts/ 判据程序与本 case 探针调用）。
- driver.f90 尾注同族记录：dens_states / eig_out 消费频谱组与本征对表——
  "No such producer exists in this tree"（同一停放家族）。
- 两成员头注明言 "frequency derivation is assembly-time mechanism, not member
  code"——成员只存表、消费表。

定性：**停放（机制空转件）而非缺陷**——fail-loud 姿态本身是设计（无表选中即具名停机，
不落静默零表）；接线上缺口是"装配期自动导出通道"（成员 init 经容器真力算
Hessian→EIGN→填表）尚未落地，登记为候选后续件。

**演示形态（按简报：驱动程序式）**：`scripts/probe_chain_freq.f90` 小驱动。
- **真装配链先行**（驱动器装配步 1-8：read_input → rng_init → sysdef_load →
  roster_load → state_create → registry_gen_all → registry_gen_init →
  container_init），容器活置，force_eval 交付绑定 reference 井的合成真力。
- **频率导出**：force_proc 形式的力包装（探针 scratch 态上跑 force_eval）喂
  `hessian_cd`（默认步长 h = 1e-4，取 roster 运行几何处）→ 质量加权
  hmw(i,j) = H(i,j)/√(m_i·m_j) → `linalg.eig_sym`（EIGN 逐字移植，升序）→
  ω = √λ → ν = ω/(2πc)。
- **成员公开面**：导出表填成员 params → 成员 init（自身校验电池在**导出表**上复跑）
  → samp_run 分发（探针在装配后改写 reactants%dist_scheme 为目标成员码——
  输入文件里 DIST_SCHEME = ROTATION+J_ROT=3 是无害占位成员，两 seam 的
  未选中守卫因此跳过 boltzmann/normalmode 的装配期 init，绕开零谱停机）。
  T_VIB_A / E_VIB 从暂存行用 seam 自身读法（seam_real）取出。

**谐振解析契约（所有锚的手算侧）**：
- 井 V = K_HARM/2·Σ|q_i|²，K_HARM = 0.04 eV/Å²；力过 force_eval 落内能制：
  K_int = 0.04×23.0605×0.04184 = 0.04×0.96485132 = **0.0385941 内能/Å²**
  （程序立法常数对；CODATA 值 0.9648535 相对差 2.3e-6——观察项，见第 3 节注）。
- 笛卡尔 H：对角全为 K_int；坐标 i 的力不读坐标 j（井各向同性、逐分量独立），
  非对角是两个位级相同力分量之差 = **恰 0.0**（中心差分对二次势精确）。
- 谱：ω_a = √(K_int/m_a) 每原子三重；ν = ω·530.8837（c = 2.99792458e-4 cm/(10 fs)）。
- **表交接约定**：boltzmann 收全 3N 谱（井钉住全部笛卡尔自由度——无自由平动）；
  normalmode 收平动投影谱（成员契约要求列与质量加权平动正交）：
  P = I − Σ_d u_d u_dᵀ，hp = P·hmw·P → 3 近零 + (3N−3) 模。**等质量（H3）使
  投影谱解析均匀** ω = √(K_int/m)（单位阵倍数的任何向量都是特征向量）——
  解析锚在投影后存活；不等质量（HCl）投影谱不再均匀，故 HCl 只走 boltz 腿。
- 两腿载体：boltz = HCl（m_H = 1.008, m_Cl = 35.45 u；ω_Cl = 0.032996、
  ω_H = 0.195673 rad/(10 fs)，ν = 17.517 / 103.879 cm⁻¹——两个可分辨解析频率
  锐化谱锚）；normalmode = H3 弯曲（三等质量，s = 3N−3 = 6 模）。
- normalmode 成员几何约定（再展开用）：qz_com = qz − COM(qz)，
  x = √m·(q−qz_com)，y = p/√m，Q = Cᵀx，P = Cᵀy，模能量 E_i = ½(P_i²+ω_i²Q_i²)。

## 3. 预期锚（先写后跑——判据先行）

常数侧：K_int = 0.0385941；ω_Cl = 0.032996、ω_H = 0.195673 rad/(10 fs)；
ν_Cl = 17.517、ν_H = 103.879 cm⁻¹。统计段 N 与临界值如下，容差沿
check_hessian / check_dist_boltzmann / check_dist_normalmode 判据形式换真频率参数。

**F1 Hessian 结构（两腿共用）**：对角 |H(i,i)−K_int|/K_int ≤ **1e-10**（相对）；
非对角 == **0.0**（按位）。沿 check_hessian (a)（1e-8 相对/1e-9 绝对）收紧：
中心差分对二次势精确，残差只来自 (q±h) 舍入（~1e-12 相对），1e-10 留 ~30 倍裕度。

**boltz 腿（HCl，N = 200000，每统计段前 rng_init(ISEED) 复种）**：

- **F2 全谱**：eig_sym 升序 ω 逐个对解析集 {3×ω_Cl, 3×ω_H} 相对 ≤ **1e-10**。
- **B1 几何 chi-square**（最高 ν 模 = 模 6，H 模）：桶 n = 0..4 外加 ≥5 尾
  （6 桶），理论 q = exp(−1.43878·ν_解析/T) = exp(−1.4371) = **0.23761**
  （T = 104.0 K；期望计数 152475/36229/8607/2046/486/152，全 ≥ 5——chi-square
  有效）；χ² df = 5，临界 **11.07**（α = 0.05）。理论侧用解析 ν（C5 法令：
  理论是锚，非成员输出）。
- **B2 逐模均值**：每模样本均值落各自 q_k/(1−q_k) 的 **5σ** 内
  （σ = √(q/(1−q)²/N)）；理论均值 Cl **3.6467**（5σ = **0.0460**）、
  H **0.3117**（5σ = **7.15e-3**）——两简并组不得串扰。
- **B3 量子数完整性**：全部抽样 n ≥ 0（getter 整数类型）。

**normalmode 腿（H3，s = 6，E_VIB = 12.5 kcal/mol → e_tot = 12.5·e_conv =
0.523 内能；统计段前 rng_init 复种）**：

- **F3 投影谱 + 表契约**：3 个近零本征值 |λ| ≤ **1e-10**（绝对）；6 模 ω 对
  √(K_int/m) = 0.195673 相对 ≤ **1e-10**；本征矢表 CᵀC−I ≤ **1e-12**；每列与
  质量加权平动 u_d 正交 ≤ **1e-12**；本征残差 ‖hp·c_k−λ_k·c_k‖ ≤ **1e-10×(K_int/m)**
  （谱标 0.0382878）。成员自身 1e-8 电池在导出表上复验（normalmode_init 绿）
  本身即链证据。
- **N1 精确恒等式**（2000 抽样）：再展开能量 Σ½(P²+ω²Q²) = e_tot 相对
  **1e-10**；|Σp| ≤ **1e-9·√(2·M·e_tot)**（M = 载体质量和 3.024 u）；
  √m·(q−qz_com) 与 p/√m 的 span(C) 投影残差 ≤ **1e-10**×向量范数；
  写出 q 的 COM 在原点 ≤ **1e-10 Å**。
- **N2 模 1 能量份额 ~ Beta(1/2,5/2)**：等 CDF 10 桶 χ² df = 9，临界 **16.92**
  （α = 0.05；闭形式分位 x = sin²θ，G(θ) = 3θ/4 + sin2θ/2 + sin4θ/16，
  G(π/2) = 3π/8，80 步二分求边），N = 50000；均值 **1/6** 于 **5σ**
  （σ = √(5/144/N) = 8.33e-4）。
- **N3 相位均衡**：平均动能 = e_tot/2 于 **5σ**（σ = e_tot·√(3/(64N))；
  Var(KE) = 3·e_tot²/64——u⁴ 矩之和 3/(s+2) = 3/8，非单模 1/16）。

**D 确定性（两腿）**：同输入双跑 stdout 逐字节一致（每统计段前复种 rng）。

**常数注（观察项，非缺陷）**：解析侧必须用程序立法换算对
23.0605×0.04184 = 0.96485132 内能/eV（force_facade 唯一换算点的两个字面量），
非 CODATA 0.9648535（相对差 2.3e-6）——否则 F1 的 1e-10 锚会假失败。
立法舍入常数姿态入任务报告观察项。

## 4. 结果对照（跑后入档）

正式跑 2026-09-16，gfortran 树（`src/build/gfortran/probe_chain_freq.exe`，两腿各在己目录），
退出码均 0；stdout 全文入档 `results.txt`。**全锚绿**：

| 锚 | 预期（第 3 节） | 实测 | 判定 |
|---|---|---|---|
| F1（boltz） | 对角 ≤1e-10 相对；非对角 ==0.0 | 最坏对角相对 2.718e-13；非对角按位 0.000E+00 | 绿 |
| F2 | 3×ω_Cl + 3×ω_H 各 ≤1e-10 相对 | 最坏 1.359e-13；ω = 0.0329953/0.1956726，ν = 17.5167/103.8794 cm⁻¹ | 绿 |
| B1 | χ² ≤ 11.07（理论 q = 0.23761） | χ² = 5.853（q = 0.237613） | 绿 |
| B2 | Cl 3.6467（5σ 0.0460）、H 0.3117（5σ 7.15e-3） | Cl 3.6490/3.6547/3.6400；H 0.3132/0.3118/0.3147——各落己 5σ | 绿 |
| B3 | 全部 n ≥ 0 | 零负数 | 绿 |
| F1（normalmode） | 同上 | 最坏对角相对 7.213e-13；非对角按位 0 | 绿 |
| F3 | 零本征值 ≤1e-10；ω ≤1e-10 相对；正交 ≤1e-12×2；残差 ≤1e-10×谱标 | 零 1.388e-17；λ 相对 5.408e-13；正交 9.992e-16；平动交叠 1.974e-16；残差 3.717e-16；6 模 ω 全 0.1956726 | 绿 |
| 成员电池 | normalmode_init 在导出表上复验 | 绿（成员自身 1e-8 容差全过——链证据） | 绿 |
| N1 | 能量/动量/span/COM 四恒等式（2000 抽样） | 全过（e_tot = 0.5230 内能精确重展开） | 绿 |
| N2 | χ² ≤ 16.92；均值 1/6 于 5σ（8.33e-4） | χ² = 4.305；均值 0.166586（差 8.1e-5） | 绿 |
| N3 | 均动能 = e_tot/2 于 5σ（2.53e-3） | 0.26064 vs 0.26150 | 绿 |
| D | 双跑 stdout 逐字节一致 | 两腿 cmp 均无差异 | 绿 |

常数注实证：F1 以立法对 0.96485132 通过（2.7e-13/7.2e-13，恰为中心差分舍入水平）——
若误用 CODATA 值将假失败于 2.3e-6，与第 3 节预期一致。

诚实记录（判据程序自身缺陷，非主体程序缺陷）：normalmode 腿首版探针的再展开例程有一处
索引笔误——COM 三向量被按原子下标切片（原子 2/3 越界读栈内存：首原子正确、后两原子被
栈垃圾污染成恒等式破坏与偶发 NaN）。`-fcheck=all` 定位，正式跑前修复；成员
dist_normalmode 自身恒等式（动量/COM/动量侧 span）全程完美，无辜。此缺陷只存在于
scripts/ 探针，主体 src/ 不受影响。

