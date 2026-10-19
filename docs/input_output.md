# 程序输入输出说明（2026/10/14）


基本格式：
- `KEYWORD = VALUE`；
- 可通过`#`/`!` 整行注释；在行尾写 `,` 可以续行；
- 关键字名大小写不敏感；一键仅允许出现一次，未知/重复键指名中止。  

单位制：
- 能量 eV，长度 Å，时间 fs，温度 K，角度 rad。
  
键表的定义参见 `src/input.f90`, 如有冲突，以源代码为准；

程序仍在持续更新调试中，本说明可能会随实际情况调整；如遇到问题，欢迎反馈；

## 程序如何读取输入文件

- 输入文件名暂时固定为 `input_qct.txt`，程序在**运行目录**下直接读取（`cd` 到案例目录再运行）；


## 反应物构型文件（SYSTEM_DIR）

- `SYSTEM_DIR` 关键字指向一个文件夹，**每个反应物片段一个文件**（`.xyz` 或 `.poscar`）

- **文件名去扩展名 = 片段名**——`PROJECTILE` / `TARGET` 按此名称指定片段；

- `.xyz` 格式：原子数行 + 注释行 + 逐原子 `元素 x y z` 行（行尾可选三个
  `T/F` 约束旗标，见下）；
- `.poscar` 格式：标准的 VASP POSCAR 文件；
- **自由度约束（DOF mask，2026-10-16 立法）**：POSCAR 的 Selective dynamics
  节与坐标行尾三个 `T/F` 旗标（以及 `.xyz` 行尾同形三列）按坐标解析——
  `T`=自由、`F`=约束。力对所有坐标照常计算（VASP 同款），积分器只推进
  自由分量；约束方向动量在采样出口清零、恒为零。**三方向全 F 的原子为
  ghost（背景）原子**：不参与采样/统计/热浴域，可被表面 `TARGET` 放置。
  速记：元素名 `G` 前缀（如 `GAu`）≡ `Au` + `F F F`（`Ga`/`Gd`/`Ge` 等
  真实元素名整体优先查表，不受影响）。


## 产物声明文件（SYSTEM_DIR/products/）

- 配合可选键 `PRODUCT_NMA = T` 使用：统计段启用末态简正模分析器（对每条轨迹的
  末态做参考构型最优旋转拟合 → 模投影 → 振动量子数），并输出通道判词；
- `SYSTEM_DIR/products/` 子文件夹中**每个产物一个 `.xyz` / `.poscar`**，文件名
  （去扩展名）= 产物单位名；无该子文件夹也合法——分析单位退为反应物片段
  （弹性/解离通道计数照常）；
- 产物为**物种语义**：不绑定固定原子序号，逐轨迹在元素相容的系统原子子集中
  取拟合 RMS 最小者为实现（H+H₂ 交换反应自动覆盖）；产物名与片段名撞车 = 中止
  （交换场景请用可区分命名，如 `H2_exch.xyz`）；
- 产物振动谱走 `SPECTRUM_SOURCE = COMPUTE` 通道（声明构型探针 + 全系统
  Hessian + 衰变门；产物谱暂无 MANUAL 通道）；
- 输出（运行目录）：每分析单位一个 `nma_<单位名>.dat`（列 `traj t drift atoms
  fit_rms e_kin_vib e_pot_vib e_har n_1..n_m`，数据面无阈值——drift 列留给后处理
  自滤）+ `nma_label.dat`（每轨迹一行通道标签）；stdout 终局汇总块内增
  `nma <单位>` 汇总标签行与 `channels` 判词计数标签行（与 classify 的
  `outcomes` 计数标签行并存）；
- 表面范式（`SCATTERING = SURFACE`）与 `PRODUCT_NMA = T` 可共存（2026-10-15 批
  选择规则，替代旧互斥句）：分析单位仅收 PROJECTILE 弹丸片段（nat≥2）与产物
  声明，衬底载体片段被排除、永不拉谱（有衬底被排除时派生期 stdout 注记一行）；
  表面+无 PROJECTILE 角色 = 装配期中止。ER 型带出产物由产物声明的物种语义天然
  支持；单原子弹丸+无产物声明 = 零单位合法形态：不写任何 `nma_*.dat` /
  `nma_label.dat`，判词面不产生。


## 常用关键字
- 基本计算流程控制，一般必须写；

| 键 | 默认 | 含义 |
|---|---|---|
| `TASK` | `TRAJECTORY` | 任务类型（目前仅此一个合法输入） |
| `N_TRAJ` | `1` | 轨迹条数 |
| `INTEGRATOR` | `VERLET` | 积分器：`VERLET` / `SYMPLECTIC` / `RADAU` |
| `ISEED` | `1` | 随机流种子 |
| `MAX_STEPS` | `100000` | 每轨迹步数上限（安全网；无终止判据时跑满为止） |

## 必需关键字
- 无论什么体系都必须写的关键字，没有默认值；

| 键 | 含义 |
|---|---|
| `DT` | 积分时间步长 [fs] |
| `SYSTEM_DIR` | 反应物构型文件夹：每片段一个 `.xyz` / `.poscar`，文件名 = 片段名|
| `E_REL` | 碰撞能 [eV] |
| `R_SEP` | 初始质心间距（投放距离）[Å] |
| `DIST_SCHEME` | 每个片段的采样方案，逗号分隔，**词数 = 片段数**；关键字可带括号参数：`<方案名>[参数=值, ...]`（见下文采样方案自带参数节） |

## 入射通道参数
- 弹丸的入射方式：控制碰撞能、方向与落点；

| 键 | 含义 |
|---|---|
| `N_E_REL` | 碰撞能模式：0=热分布 / 1=定值(给定`E_REL`) / 2=密度加权 |
| `T_TRANS` | 平动温度 [K]（N_E_REL=0 时须 >0） |
| `V_WIDTH` | 速度展宽 [Å/fs]（N_E_REL=2 时须 >0） |
| `N_THTA` | 极角模式：0=固定 / 1=按 sinθ 采样 |
| `THTA_MAX` | 极角界/固定值 [rad]（0..π/2） |
| `N_CHI` | 方位角模式：0=固定 / 1=均匀采样 |
| `CHI` | 固定方位角 [rad]（N_CHI=0 时用） |
| `N_B` | 瞄准模式：0=固定 b=B_MAX / 1=均匀采样 |
| `N_AIM` | 表面散射落点模式：0=原点 / 1=胞内均匀 / 2=定点 |
| `AIM_X` / `AIM_Y` | 定点落点 x/y [Å] |


## 条件关键字
- 填写方式与散射类型（`SCATTERING`）相关的关键字

| 键 | 含义 | 表面散射填写 | 气相散射填写 |
|---|---|---|---|
| `SCATTERING` | 决定体系是表面散射还是气相散射 | `SURFACE` |`GAS` |
| `CELL` | 晶胞三行晶矢 [Å]| 如 `6 0 0, 0 6 0, 0 0 20` | 不支持 |
| `PROJECTILE` | 作为弹丸的分子片段的名称 | 分子片段名称 | 与 `TARGET` 不同名的片段 |
| `TARGET` | 靶片段名（2026-10-16 起表面为**可选**的显式衬底片段——入射通道将其原子复位到文件晶格位，热浴域与统计排除以它为锚；不写则维持隐式表面参考系） | 可选：衬底片段名 | 与`PROJECTILE`不同名的片段 |
| `B_MAX` | 碰撞参数最大值 [Å] | > 0的数值 |不支持（表面晶格随机入射）|

## 可选关键字
- 不写即取默认值；

| 键 | 默认 | 含义 |
|---|---|---|
| `TITLE` | （空） | 输出中的标题行 |
| `ELEC_METHOD` | `ADIABATIC` | 非绝热方法的注册名，参见自定义PES接入说明 |
| `GWRITE_LEVEL` | `2` | 记录级别，限 {0,1,2,5}：≥2 记每步 q/p，≥5 额外写入受力 |
| `OBSERVABLES` | （空） | 让程序额外输出的物理量的名称 |
| `RANDOM_ORIENT` | `F` | 是否随机化各片段的初始构型 |
| `PRODUCT_NMA` | `F` | 是否启用末态简正模分析器（统计段逐轨迹分析 + 通道判词；产物声明见"产物声明文件"节） |
| `NMA_FIT_TOL` | `0.0` | 判词拟合阈值 [Å]；0 = 逐单位自动取 0.4 × 参考构型最短原子间距（只影响判词贴标签，不影响 fit_rms 数据列） |
| `ORDER` | `1` | 积分器阶：verlet 1=velocity-Verlet、2=三阶 Beeman；symplectic 组合阶 4/6/8 |
| `BATH` | `NONE` | 是否在轨迹开始前做热浴：`NONE` / `ANDERSEN` / `NOSE_HOOVER` / `THERMALIZE`（2026-10-16 迁入：一次性 Maxwell 播种+自由平衡跑——`seed` 形态，环内无温控维持；作用域=衬底域，见 TARGET）** |
| `N_BATH` | `0` | 浴循环步数 |
| `DT_BATH` | `0.0` | 浴循环步长 [fs] |
| `T_BATH` | — | 浴温度 [K]（`BATH=ANDERSEN`/`NOSE_HOOVER` 时须 >0） |
| `NU_COLL` | — | Andersen 碰撞频率 [1/fs]（`BATH=ANDERSEN` 时须 >0） |
| `TAU_NOSE` | — | Nosé–Hoover 链周期 [fs]（`BATH=NOSE_HOOVER` 时须 >0；Nosé 处方 Q_k=f_k·kb·T·(τ/2π)²，取近体系振动周期） |
| `SPECTRUM_SOURCE` | `COMPUTE` | 振动谱（部分采样方式需要）的计算方式：`COMPUTE` =  Hessian 直接计算；`MANUAL` = 手动输入参数 |
| `SPECTRUM_FILE` | （空） | 谱表文件路径（须 `SPECTRUM_SOURCE=MANUAL`；） |

**注**：`NOSE_HOOVER`（2026-10-15 入册）为双链 Nosé–Hoover 恒温浴，step 接管形态
（浴环内自带 Liouville VV 分裂，不依赖注册积分器）；逐轨迹链动量按正则分布重播
+ 首步无主原子 Maxwell 播种（环内抽签豁免内）；判据 check_bath_nose_hoover。

## 采样方案自带参数
- 各采样方法的参数随 DIST_SCHEME 处关键字的子括号给出，形式为：
  `DIST_SCHEME = <方案名称>[参数=值, 参数=值], <方案名称>[...], ...`

| 方案词 | 括号参数 | 含义 |
|---|---|---|
| `BOLTZMANN` | `T_VIB` | Boltzmann 振动温度 [K] |
| `EBK` | `N_VIB` / `N_ROT` | 振动/转动量子数 |
| `ROTATION` | `J_ROT` | 固定转动角动量量子数 J |
| `J` | `T_ROT` | 转动温度 [K]（ |
| `NORMALMODE` | `E_VIB` | 振动能 [eV] |
| `BARRIER_EXCITATION` | `E_STAB` | 鞍点稳定模分能 [eV] |
| | `N_E_BAR` | 反应坐标能量模式：0=定值 / 1=热（缺省 0） |
| | `E_BAR` | 反应坐标能量 [eV]（≥0） |
| | `T_BAR` | 反应坐标温度 [K] |
| `GLO_TARGET` | `T_GLO` | GLO模型的温度 [K] |
| | `GAMMA_GLO` | GLO 摩擦系数 |
| | `W_GHOST` | ghost 振子频率 [1/fs] |
| `WIGNER` | `E_SCALE` | Wigner 能量标度（零点能份额，须 >0） |
| `EROT` | `T_ROT` | 非线性片段热转动温度 [K]（三主轴各 ½kbT 经典热抽，须 >0；线性转子用 `J`+`ROTATION`） |
| `NONE` | （无） | 该片段不采内态：单原子片段的诚实用词（替代历史上的哑成员词），亦可用于不采内态的多原子载体（如冻结基底）；无成员、无参数、派发跳过，`NONE[...]` 带参数会被审计中止 |

**注**：本部分仅有EBK、Boltzmann两个方法进行过测试，目前不保证其余方案实际有效，且后续可能会进行更多调整；

示例（两个 EBK 片段各处一个量子态）：
```
DIST_SCHEME = EBK[N_VIB=0, N_ROT=2], EBK[N_VIB=3, N_ROT=1]
```

示例（单原子弹丸片段用 NONE，靶分子 EBK 采样）：
```
DIST_SCHEME = NONE, EBK[N_VIB=0, N_ROT=2]
```

## 其他输入

可由自定义势能面额外声明，随
input_qct.txt 一起写，接入
自定义势能面见 [custom_pes.md](custom_pes.md)。

## 运行输出

在运行目录（即 `input_qct.txt` 所在处）产生下列文件：

| 文件 | 类别 | 内容 |
|---|---|---|
| `trajectories.txt` | 轨迹记录 | 文件头三行图例 + 逐帧块 + 每轨迹一条末态行 |
| `hist_e_rel.dat` / `hist_thta.dat` | 统计直方图 | 末态相对平动能 / 散射角分布 |
| stdout | 运行汇总 | 开场横幅 + 装配读数；终局重框汇总块（统计 / 通道计数 / 守恒漂移 / 计时） |
| stderr | 报错 | 正常运行为空；出错为报错框（`!` 线框：类别 / 文件 / 站名 + 指名中止消息行） |

以下样例均取自 [cases/case_harmonic_scatter](../cases/case_harmonic_scatter)
的实际运行产物。

### trajectories.txt —— 文件头图例 + 帧块（每步一块）

文件开头三行图例（单位与行形自描述；观测列名随 `OBSERVABLES` 选择动态生成，
`level <n>` 随 `GWRITE_LEVEL` 落档），随后逐帧块：

```text
# VENUS2026 trajectory record (units: eV | A | amu | fs; p amu*A/fs; f eV/A)
# each frame: '# frame traj=<id> step=<n> t=<t [fs]>' + observable row: t_kin=<v [eV]> v_pot=<v [eV]> e_tot=<v [eV]> (level 2)
# '# final: traj t t_kin v_pot e_tot' header + one final row per trajectory (t [fs], values [eV])
# frame traj=1 step=1 t=  0.0000000000000000E+00
t_kin=8.8210396325053004E+00 v_pot=2.6086196116948859E+00 e_tot=1.1429659244200186E+01
# sec_qp begin [q A | p amu*A/fs]
1 H     -1.6914187273464129E-01  1.0720118906927656E+00  7.8311562014828864E+00  1.7402814056400123E-02  1.4291265564038527E-02 -2.1870780307566229E-01
2 H      … （每原子一行：序号 元素 qx qy qz px py pz）
# sec_qp end
```

- 数值列默认 `t_kin/v_pot/e_tot` 三列； `OBSERVABLES ` 选用后
  追加到本行。

- 节标 begin 行带单位后缀（`# sec_qp begin [q A | p amu*A/fs]`、
  `# sec_frc begin [f eV/A]`——`GWRITE_LEVEL>=5` 时写入受力节）；end 行朴素无后缀。

- 单位：t [fs]、能量 [eV]、坐标 [Å]、动量 [amu·Å/fs]、受力 [eV/Å]。

### trajectories.txt —— 末态行（每轨迹一行）

```text
# final: traj t t_kin v_pot e_tot
1   5.5500000000000000E+01   6.6910000611638925E+00   4.7387101558770137E+00   1.1429710217040906E+01
```

### stdout（开场横幅 + 装配块；终局汇总块）

```text
=====================================================================
 VENUS2026  adiabatic QCT
=====================================================================
 assembly     trajectories=8 atoms=4 fragments=2 paradigm=0
              members=11 init=ok columns=3
=====================================================================
 VENUS2026  run summary
=====================================================================
 statistics   trajectories=8 valid=8
 outcomes      scattering=8
 drift        max|e_fin-e0|=  6.1156E-05 eV
 histograms   rows=8 files=hist_e_rel.dat hist_thta.dat
---------------------------------------------------------------------
 timing       cpu=    0.12 s
=====================================================================
```

- 重线 `=`×69 开合、轻线 `-`×69 分段；12 字符标签列 + 1 空格（内容第 15 列起，
  续行缩进至值列——装配块的 members 行即续行形态）；数值 token 形态逐字保留。
- `outcomes` 标签行即 `classify` 返回的通道名的去向（容器限定时标签为
  `outcomes (container)`；无分类绑定时为
  `outcomes     none (no classification bound - skipped)`）。
- 启用末态简正模分析器时，`histograms` 行后接 `nma <单位>` 汇总标签行与
  `channels` 判词计数标签行（见"产物声明文件"节）。

### 直方图 .dat

```text
    0.097960    0.000000
    0.293881    0.000000
    0.489802    0.000000
```

两列 = 箱中心、计数；等宽分箱（[lo, hi)），行数即箱数。
