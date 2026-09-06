# VENUS2026 程序完备化最终计划

> **执行者须知**：本计划按"阶段-任务-步骤"组织，步骤用复选框（`- [ ]`）跟踪。
> 每任务以"验收命令通过 + git 提交"结束；提交 message 引用 audit 条目名（如
> `P-B1: GWRITE 观测量退为配置（audit: GWRITE）`）。执行前先读
> `docs/design_spec_completeness.md`（裁决标准）与 `venus2026/CLAUDE.md`（构建/单位约定）。
> **展示协议**（用户裁决 2026-09-06）：一切方案/骨架/实现的展示一律"先大纲、
> 经确认、再扩充"——每次只给当前层大纲，经用户确认后才展开下一层，禁止一次
> 性倾倒全部具体内容。
> **执行策略**（用户裁决同日）：骨架优先——Phase K 先生成全部占位文件（全局
> 变量真实定义 + 每例程输入/输出签名与伪代码），结构审核确认后逐文件具体化。
> **独立新程序**（用户裁决同日二次修订）：目标是交付一份独立的新程序；旧代码
> （venus2026/）仅作纯参考（算法与公式出处），全程只读——不修改、不清理、
> 不迁移、不退役。
> **分阶段计划文件**（用户裁决同日）：每执行一个小 phase（如 K1、C2），先提交
> 该 phase 的单独计划文件 `docs/plans/<phase>.md`（目标/交付物/文件清单/步骤/
> 验收），经确认后再执行；本计划只保留序列、映射与量尺。

状态：**执行基准**（本文档为唯一执行依据）
日期：2026-09-06（同日删除 remediation_plan.md 与 audit_brief.md——前者的纯化步骤
已全部并入本文，后者为已完成的审查委托书；避免双头执行与文档混乱）
依据链：design_spec_completeness.md（规范）→ completeness_audit.md/.json（质料清单 99 条目）
→ optimal_structure.md（目标结构与质料映射，**本计划基准**）

**目标：** 交付一份**独立的新程序**——目标结构即 optimal_structure.md §1 的
L0-L5 分层，可变代码集中于接口 A/B′ 成员与封闭族注册行。旧代码（venus2026/）
仅作纯参考。**验证体系重建**（用户裁决 2026-09-06）：旧 10 个 case 全部舍弃
作基准；新验证 case 与判据在架构（K 骨架）落成后随 C 组件逐个拟定，V 阶段
汇编成集并全量执行。

**架构：** 独立新程序 + 骨架优先（2026-09-06 用户裁决，两次修订）。新程序自成
一棵源码树（含自己的构建）；旧库 venus2026/ 只读参考。Phase K 生成全部 L0-L5
占位骨架（全局变量真实定义 + 每例程签名/伪代码，新树独立编译），经审核门后
Phase C 逐文件具体化（按 §3 映射回旧库参考算法，组件级新旧对照验证），最后
Phase V 执行重建的验证 case 体系（判据覆盖率检查 + 全量运行）。每小 phase
执行前先出单独计划文件。

**技术栈：** Fortran（ifx，oneAPI 2025.2 验证）、make（单入口 Makefile）、git。
测试方法适配说明：本项目的"测试环"= 新验证 case（最小输入 + 解析/统计锚点 +
守恒律判据，随 C 组件拟定）+ 临时同量参照（旧程序可算同一物理量时对照，非
基准），不以单元测试框架为载体。

**Spec：** docs/design_spec_completeness.md；结构映射细则 docs/optimal_structure.md §2。

## 全局约束（每任务隐含遵守）

1. 内部单位制 Å/amu/(10 fs)，换算常数 C1=0.04184、C7=0.063508、23.0605（eV→kcal/mol）、
   C8=1.987e-3；换算只发生在 force_facade 边界（单位契约履行点）。
2. 编译 `make`（venus2026/ 下；Windows 在 Git Bash 执行，链接器路径变量覆盖见 CLAUDE.md）；
   flags `-r8 -double-size=64 -i8 -O2 -w`。
3. 输入只认 `KEYWORD=VALUE` 文法（input_parser）；新关键字一律走
   param_mapping 注册表模式，旧整数关键字回退保留。
4. 演化期不调用操作类（相位限制，规范 2.5）；制备成员内部调演化服务合法（边界 4）。
5. 单一随机流：新程序自 C1 起即单流（rng.f90，Knuth RAND0 族，ISEED 状态），
   全库无第二流。
6. 改模块文件后删除陈旧 `.mod` 再编（CLAUDE.md 惯例）。
7. 行号引用以 2026-09-04 audit 为准；执行时若行号漂移，以符号名 grep 重新定位。
8. C 各步以"最小验证输入 + 解析/统计判据"收口（判据先行：先定判据再写实现）；
   V/S 全链路执行新 case 集；旧程序仅作临时同量参照，不作门槛。

---

## 0. 基线核实（2026-09-06，与 audit 的差异）

执行前已核实，以下为既成事实，计划据此调整：

| 核实项 | 现状 | 计划含义 |
|---|---|---|
| `src_VENUS/GLO.f`（未编译滞留副本） | **已不在树** | P-A1 对应条目自然完成 |
| `ROTATEX_b/ROTATEZ_b` | **已不在树** | 同上 |
| GASDEV.f/VFDATE.f/POTEN.f/EBOND.f/THERMBATH.f/WEBOND/WLBOND/WENMOD/MICROCI/VOLPSCONE/LMODE/LMEXCT | 仍在 | P-A1/A2 照删 |
| gasdev 住 interface_TEST.f90:112（函数） | 仍在 | P-C1 照做 |
| GWRITE :59/95-99/110-114、TEST :88-94、VERLET :47,74,137 ghost/GLO 内嵌 | 仍在，行号与 audit 一致 | P-B/P-D 照做 |
| Makefile MKL 三库链接 | 仍在 | P-A5 照做 |
| `venus2026/.git` | **不存在** | Phase 0 必须先建版本控制 |
| cases/ 数量 | **10 个**（boltzmann_vib、ebk_fixed_nj、fixed_normalmode_qnums、integrator_matrix、leps_potential、mb_thermal、rigid_surface、rst_beam_scattering、rst_surface_md、twobody_collision） | **全部舍弃作验证基准**（用户裁决 2026-09-06：旧 case 是"旧行为保持"用例，判定依据参差且烙有 F15/双随机流等缺陷）；目录留旧库作历史参考；新 case 体系待 K 骨架后拟定。CLAUDE.md"16 个"过时，P0 修正 |
| docs/completeness_audit.json | 在 | S2 复审对象 |

## 1. 目标终态（验收态定义）

结构（optimal_structure §1 摘录，文件名即归位目标）：

```
L0 控制核        driver.f90          装配→演化→统计的调用序列，无物理
L1 物自体容器    state.f90           Q, P, s, F, t（可分配，维数=名册内容参数）
L1 表象层        config.f90          场景、名册、分布、谓词、结局分类、观测量配置
L2 接口 A 门面   force_facade.f90    (Q,P,s)→(F,E)；内部 Σ力贡献成员 × 单选电子传播律
L2 接口 B′框架   sampler.f90         分布→(Q0,P0,s0) 通用抽样器
L3 注册族成员    pes_*.f90 / dist_*.f90
L3 封闭族        propagator.f90 / predicate.f90
L4 测量          recorder.f90        轨迹序列化 + 谓词求值
L4 后处理        analysis/           结局分类执行、统计、绘图（可离线）
L5 仪器          rng.f90, linalg.f90, geometry.f90
```

终态判据（T1-T7，收尾 S 逐条核验）：

- **T1 文件结构**：上表文件全部在场，Makefile SRC 按新依赖序排列；
- **T2 可变代码位置**：新增一个 PES/分布成员 = 新文件 + 一行注册，既有代码零触及
  （以"新增 pes_leps 变体试件"演练验证）；
- **T3 双实现矩阵**：A 轴 HARMONIC+MORSE+LEPS+EMT-NN+RST（≥2 ✓）；B′ 轴
  THRMAN+INITEBK+dist_surface_oscillator（≥3）；构型层面分子 vs 表面两范式
  （新验证 case 各设代表）；
- **T4 单位契约**：grep 全库，23.0605 与 C1 力换算只出现在 force_facade 边界；
- **T5 单一随机流**：全库无 `RANDOM_NUMBER` 残留（rng.f90 独此一流）；
- **T6 三相位**：driver 只含装配→演化→统计调用序列；演化环内只有
  facade/propagator/recorder/predicate；
- **T7 完备性复审**：对新程序全例程按规范三问复审，产出新 completeness 审计
  （md+json）；除不变类/成员类外无残余。旧 audit 保留为历史参考。

## 2. 执行序列（2026-09-06 用户裁决：独立新程序 + 骨架优先 + 分阶段计划文件）

```
Phase 0 准备（git + 新程序树立项）
   0.1 git init（工作区根；旧库自此只读）
   0.2 新程序树立项（src/：Makefile 骨架 + README + analysis/ 占位）
   0.3 旧库 CLAUDE.md case 计数修正
   ↓
Phase K 骨架生成（新树内全部 L0-L5 文件占位）
   K1 全局变量真实定义（state.f90 / config.f90 / 常数；
       venus_data 约 40 表象块的逐块去向分配）
   K2 每例程占位：输入/输出签名 + 伪代码 + 参考出处注释（旧库 文件:行）
   K3 新树独立编译通过（桩主程序；编译器验证接口一致性）
   ↓ 【审核门：只展示骨架大纲（文件树 + 每文件例程一行清单），
        逐层确认后放行；伪代码逐文件按展示协议再展开】
Phase C 具体化（逐文件；每文件先大纲、确认后扩充；组件级新旧对照）
   C1 L5 仪器：rng（单一随机流）/ linalg / geometry   ← 参考 RAND0 族/gasdev/EIGN/ROTATE 族
   C2 force_facade + elec_adiabatic + pes_* 成员       ← 参考 interface_TEST/emt_nn/rst/VERLET 内嵌力
   C3 predicate 文法 + 终止判据装配                   ← 参考 TEST.f:106-146
   C4 recorder + config_roster + 输入文法             ← 参考 GWRITE/HIST（剥离 audit 所列违规）
   C5 sampler + dist_* 成员表                         ← 参考 SELECT/THRMAN/INITEBK/SURF/…
   C6 propagator 封闭族（verlet/symple/radau）        ← 参考 VERLET/SYMPLE/RADAU（剥编排）
   C7 driver + 结局分类数据                           ← 参考 VENUS.f90 三相位编排
   C8 analysis/ 离线族                                ← 参考 FINAL/GFINAL/EIGOUT/plot_*
   ↓
Phase V 全链路验证（C 阶段各组件验证输入与判据汇编成正式 case 集；
   判据覆盖率检查 + 全量执行 + 同 ISEED 复现性抽验）
   ↓
S 收尾（新程序完备性复审 / 文档收口 / spec §4 回填）
```

要点：
- **旧库只读**：venus2026/ 不修改、不清理死代码、不退役（原 C0/T 两类旧树操作
  废弃）；audit 的 文件:行 引用仅作参考定位。
- **分阶段计划文件**：每执行一个小 phase，先提交 `docs/plans/<phase>.md`
  （P0/K1/K2/K3/C1…C8/V/S），经确认后执行；本计划保留序列、映射与量尺。
- **验证体系重建**（用户裁决 2026-09-06）：旧 10 case 全部舍弃作基准；各 C 步
  拟定自己的最小验证输入 + 判据（解析/统计锚点 > 守恒律/分布检验；判据设计可
  借鉴旧 case README 的锚点思想，不沿用其输入与输出基线），V 阶段汇编成正式
  case 集。临时同量参照：旧程序能算同一物理量时可作对照（非基准、非门槛）。
- **输入文法**：以旧 KEYWORD=VALUE 为基础扩展构造子新关键字
  （TERMINATE/OBSERVABLES/ELEMENTS/…，文法定案在 C4）；临时对照旧程序时
  可复用旧输入。
- **骨架漂移禁令**：具体化中发现签名不适用，回审核门出示差异、确认后修订。
- C 序即依赖序：仪器→门面→谓词/记录器→抽样→传播→driver→后处理；
  C2 需 C1（成员内部调 rng），C7 需 C2-C6（driver 只编排）。
- §5（原 N1-N3）= 骨架中三个关键文件的接口规格；§6（原 P-A…P-G）= 移植时的
  剥离/操作要点（audit 违规即"参考时勿带入"清单）；§7（原 X 轨）= 新树归位核对。

## 3. 落点映射表（旧库关注点 → 新树落点）

（2026-09-06 二次修订：独立新程序语境下，本表语义 = 各关注点在新树的落位；
"原落点"列保留历史出处。旧库 venus2026/ 只读，下列目标文件均指新树。）

| 步骤 | 原落点 | 最终落点 |
|---|---|---|
| B1 GWRITE 观测量配置化 | GWRITE 内部表 | recorder.f90 列注册（N3 产出的机制） |
| B2 TEST 通用判据登记 | TEST.f 内部 | predicate.f90 文法成员（N2 产出） |
| C1 gasdev/GAMA 迁移 | "src_VENUS 独立仪器模块" | rng.f90（L5，与 RAND0/RAND1/RANDST 同居） |
| D1/D2/D4 随机律成员化 | 接口 A 成员（位置未定） | pes_glo_spring / pes_langevin_ghost / pes_ldfa，注册进 force_facade.f90 |
| D3 ghost 参数装配 | "成员初始化钩子" | facade 成员 init 协议（POTPRE 型） |
| E1 结局块合并+分类 | VENUS.f90 单一路径 | driver 调 analysis 结局分类执行器 + △6 分类数据 |
| E2 观测积累 | "测量记录器模块" | recorder.f90 终态列 |
| E3 体系分支 | "装配层" | 构型装配器（roster，N3 地基） |
| G1/G3 反写表象改交付 | "输入文件或构型装配器" | 构型装配器 + config 文法 |
| G4 RADOUT 回调 | "向控制核回调" | driver 的输出钩子 + predicate 求值点 |

冲突裁决（两文档不一致处，以基准 optimal_structure 为准；2026-09-06 二次修订后
裁决 1/2 自动消解——旧树不动，新树直接按基准映射落位，无删/留之争）：

1. **ENMODE**：remediation A2 删除（MPLOT 死）；optimal_structure 判 ◐ 随
   dist_normalmode 放置。裁决：**不随 A2 删**；其 MPLOT 死调用链删除，X3 时归位为
   dist_normalmode 的装配期分析工具；届时若全库确无消费者，按 audit 原裁决删除并记录。
2. **FRICFORCE/DIAG2**：A2 删（死入口）与 D4 复活（LDFA 成员）并存。裁决：A2 中
   从 VERLET.f90 **删除**（连同 DIAG2/LINEQ/LUDCMP/LUBKSB/FEHLER），物理内容在
   D4 以 pes_ldfa 成员形态在门面重生（素材经 git 历史可取）。
3. **case 数量**：16 → 10（见 §0）。
4. **harmonic_sampling**：audit 判"接入 B′ 或删除"。裁决：接入——X2 注册为
   dist_surface_oscillator（B′ 第三实现），SELECT 成员表加行。

编号对照：optimal_structure §5 所称 R0→R7 即原 remediation_plan（已删除）的
P-A…P-G+收尾；本文沿用 P-X 编号并新增 N-X/X-X/S-X。

---

## 4. Phase 0 准备（版本控制与新程序树立项）

**Files:** 工作区根 `.git`（新建）、`docs/plans/`（分阶段计划文件目录）、
新程序树 `src/`（命名已在 P0 计划文件定案）。

- [x] **0.1 git 初始化**：工作区根 `git init && git add -A && git commit`（docs/、
  venus2026/ 参考、README 入管；旧库自此只读）。可选接远程（push 前征得用户同意）。
- [x] **0.2 新程序树立项**：`src/Makefile`（ifx，flags 同旧库约定，Windows
  链接器遮蔽处理段移植旧库平台段）+ `src/README.md`（L0-L5 一览，占位说明
  K 阶段生成）+ `src/analysis/README.md`（L4 用途）。
- [x] **0.3 旧库文档数字修正**：venus2026/CLAUDE.md"16 个示例 case"→10（避免
  参考误导；此为旧库唯一允许的改动）。
- [x] 提交 `P0: 版本控制与新程序树立项`。执行细则见 docs/plans/P0.md（v2）。

## 5. 骨架接口规格：facade / predicate / recorder+roster（原乙轨 N1-N3）

（2026-09-06 修订：N1-N3 并入 Phase K，作骨架中三个关键文件的接口规格；
check_* 验证程序挪至对应 C 步骤验收。本节签名与文法定义对 K 直接生效。）

### N1 force_facade.f90 + elec_adiabatic（△1+△5）

**Files:**
- Create: `src/force_facade.f90`（新树 L2）
- Create: `src/elec_adiabatic.f90`
- Test: `scripts/check_facade.f90`（组件级对照，独立编译）

**Interfaces（产出，后续任务依赖）：**
```fortran
module force_facade
  public :: facade_register_pes(name, pes_proc)      ! 甲类成员注册（可叠加）
  public :: facade_register_elec(name, elec_proc)    ! 乙类成员注册（单选）
  public :: facade_init(config)                      ! 组装成员表+能力级协商，非法组合 loud-stop
  public :: facade_force(Q, P, s, F, E)              ! (Q,P,s)→(F,E)；F=Σ成员贡献；E=T(P)+ΣV
  public :: facade_info                             ! 成员清单/能力级（诊断输出）
end module
! 甲类成员签名：Q,W → F_add(code 单位), E_add；随机律成员内部调 rng（签名承认随机律）
! 乙类成员签名(退化点)：elec_adiabatic: ds/dt≡0，力规则恒等——显式存在约 20 行
```
能力级协商表（facade_init 内数据）：合法三元组当前仅
`(single_surface, adiabatic, {verlet,symple,radau})`；非法组合打印冲突两轴并 STOP。

- [ ] **N1.1** 写 force_facade.f90：成员注册表（name→procedure 数组）、加性合成循环、
  能力级协商。换算表达式与 DPESHON 完全同形（`-grad×23.0605×C1` 累加），保证 bit 级可比。
- [ ] **N1.2** 写 elec_adiabatic.f90（退化成员：ds/dt≡0、力规则恒等）。
- [ ] **N1.3** 写五个甲类成员文件：pes_reference（HARMONIC/MORSE/LEPS 体系无关
  参考势，参考 test_potentials.f90）、pes_emt_nn、pes_rst（算法自旧库对应文件
  移植，数据外置与路径关键字沿用），注册进门面。
- [ ] **N1.4** 写 scripts/check_facade.f90：临时构造输入跑旧程序取参照轨迹
  （≥100 帧；或对参考势用解析梯度手算），同 Q 下分别经 `DVDQ_1→DPESHON`（旧路）与
  `facade_force`（新路）求 F/E，逐元素 diff。
- [ ] **N1.5 验收**：`make` 通过；check_facade 全帧 max|ΔF|=0、max|ΔE|=0（bit 级，
  不为零则查换算次序并修至为零）。提交 `N1: force_facade 门面+退化电子成员（旁路验证 bit 级一致）`。

### N2 predicate.f90 文法解释器（△2）

**Files:**
- Create: `src/predicate.f90`（新树 L3 封闭族）
- Test: `scripts/check_predicate.f90`

**Interfaces（产出）：**
```fortran
module predicate
  public :: predicate_parse(str, expr, ios, errmsg)  ! 输入文法表达式 → expr 树
  public :: predicate_eval(Q, P, W, roster, expr) → logical
end module
```
输入文法（新关键字 `TERMINATE =`，param_mapping 注册）：
```
表达式 := 原式 | AND(e1,e2) | OR(e1,e2) | NOT(e)
原式   := COM_DISTANCE(片段A, 片段B, 阈值Å[, 退行速度阈])   ← TEST.f:106-146 原型
        | Z_HEIGHT(片段, 阈值Å[, OUTWARD])                  ← 取代 TEST.f:88-94 写死判据
        | MAX_STEPS(步数)                                   ← 取代 NC 魔数（恒为安全网）
        | BOUND(片段A, 片段B, 下阈, 上阈)                    ← 有界窗（约束散射）
片段   := 名册片段名（A/B/…）或 ATOM(i)
```

- [ ] **N2.1** 写文法数据结构（原式/布尔组合节点）+ 解析器 + 求值器（COM 距离与
  退行速度算法逐行取自 TEST.f:106-146，保证数值同源）。
- [ ] **N2.2** 无 `TERMINATE` 关键字时的默认表达式 = 现行 TEST 通用判据同参版本
  （RCRIT 等既有参数映射进表达式），保证 B2 前行为不变。
- [ ] **N2.3** 写 check_predicate.f90：对临时生成的旧程序轨迹回放（或构造
  状态序列），断言每步布尔值与
  TEST.f 现行输出一致。
- [ ] **N2.4 验收**：`make` 通过；check_predicate 全轨迹一致。提交 `N2: 谓词文法解释器`。

### N3 recorder.f90 列注册 + roster 名册（△4+△3 地基）

**Files:**
- Create: `src/recorder.f90`、`src/config_roster.f90`（新树）
- Test: 组件级对照在 C4 接线时进行；本件自测走编译+冒烟

**Interfaces（产出）：**
```fortran
module recorder
  public :: recorder_register_col(name, col_proc)    ! 列 = 名字 + 抽取器
  public :: recorder_select(cols_spec)                ! OBSERVABLES= 名单或 GWRITE_LEVEL 预设组
  public :: recorder_frame(Q, P, F, unit)             ! 写一帧所选列
end module
module config_roster
  public :: atom_spec(symbol, mass, fragment)         ! 名册元素结构
  public :: roster_load(ELEMENTS= 列表) / roster_symbol(i) → 'C '/'Au '（原格式宽度）
end module
```
新关键字：`ELEMENTS = C, AU(111), …`（顺序对应原子序，缺省回退 `AT<i>` 标签）；
`OBSERVABLES = 列名,…`（与 GWRITE_LEVEL 预设组并存，后者引用列组表）。

- [ ] **N3.1** 写 recorder 列注册机制（列组预设表数据化：level 0/1/2/3 → 列名单）。
- [ ] **N3.2** 写 config_roster：元素符号表（符号→2 字符标签；质量仍走现行 W 装配）、
  ELEMENTS 关键字接入。
- [ ] **N3.3 验收**：`make` 通过；全 case 回归不变（新模块尚无人消费）。
  提交 `N3: recorder 列注册与名册数据结构`。

## 6. 质料操作清单：原纯化主链（P-A…P-G）

（2026-09-06 二次修订：旧库只读，本节步骤均不在旧库执行。其角色 = **移植参考
时的剥离/操作要点**——audit 所列违规点即"参考时勿带入新树"的清单，行号与
条目名用于定位参考代码。P-A 死代码清除整节作废：旧库不清理，新树无死代码
可清（audit A8 清单 = 新树不必实现的功能清单）。各任务内嵌的"验收"行是旧计划
的旧树验收，已被 §2/§9 的新验收方式（组件判据验证 + V 全链路新 case 体系）
取代，仅步骤要点有效。）

### P-A 死代码清除（零物理风险）

**Files:** 删 A1/A2 所列文件；改 VENUS.f90、venus_data.f90、VERLET.f90、GLO.f90、Makefile。

- [ ] **A1** 删文件：`src_VENUS/{GASDEV.f, VFDATE.f, POTEN.f, EBOND.f, THERMBATH.f,
  WEBOND.f, WLBOND.f, WENMOD.f}`（GLO.f、ROTATEX_b/ROTATEZ_b 已不在树；核实 Makefile
  SRC 无残留引用）。
- [ ] **A2** 删死入口及依赖链：MICROCI.f、VOLPSCONE.f、LMODE.f、LMEXCT.f（EBOND 已在
  A1）；ENMODE.f 的 MPLOT 调用链删除（文件本体保留待 X3 裁决，见 §3 冲突 1），
  venus_input.f90:692 的 MPLOT=0 写死关键字一并清；
  SURF.f 内 PBC；VERLET.f90 内 FRICFORCE+DIAG2+LINEQ+LUDCMP+LUBKSB+FEHLER（物理素材
  D4 重生）；GLO.f90 内 GLOEQU、DISTRIBUTE。
- [ ] **A3** VENUS.f90 死块：NCHKP 分支（:182-232）、RKUTTA 系数（:86-94）、idplt
  直方图死块（:286-299，含 Qr1/Qr2 位点常数与 SYSTEM 壳调用）、NCVIB 死关键字。
- [ ] **A4** venus_data 死块：MORSEB、STRETB、CUBEB、COORS、GAUSS2/GAUHES、GDD、
  HIST_ATOM_NAD、U0/U1/V01、RKUTTA 块。
- [ ] **A5** Makefile：删 MKL 三库链接（全库无调用）与已删文件 SRC 条目。
- [ ] **验收**：`make` 通过；`bash scripts/regress.sh` 10 case 全 PASS（bit 级）。
  提交 `P-A: 死代码清除（audit A8 清单）`。

### P-B 测量层纯化（观测量选择退为配置）

**Files:** 改 GWRITE.f90、TEST.f、GFINAL.f、FINAL.f、POTENZ.f、HIST.f90；接线 recorder/predicate/roster。

- [ ] **B1** GWRITE 剥体系知识：删 atom1=C 特写行（:59-60）、质量阈元素名（:110-114，
  改 roster_symbol）、nearest-Au 诊断（:95-99，移 analysis 作可选离线工具）、+6.0 Å
  PBC 盒（:78-79，改配置）、:43 内联 kB 改用 venus_params 常数；GWRITE_LEVEL 各层
  观测量集合改走 recorder 列组表。
- [ ] **B2** TEST.f：删活的 C/Au 判据（:88-94）与两段注释体系知识（:54-78、:149-192）；
  终止判据改 predicate.f90 装配（无 TERMINATE 关键字走 N2.2 默认表达式）。受影响
  case（rst_beam_scattering、rst_surface_md、rigid_surface 等 NSURF>0 者）在
  input_qct.txt 增加与现行判据等价的 TERMINATE 行（阈值取现值；等价翻译以
  TEST.f 现文为准）。
- [ ] **B3** GFINAL.f:80 删写死原子 3 的 EREL 覆盖，采用 FINAL 广义 COM 相对能量；
  同步删 :66-69 current_U0 未初始化打印。
- [ ] **B4** 删 fort.777/666 的 U0/U1/V01 打印（FINAL.f:131-135、POTENZ.f:76-82，
  全库无写入者）。
- [ ] **B5** record_all/HIST：auro=2.94998 与六方格子（PI/3、top/brg/hcp/fcc）改
  slab 名册参数（新关键字 `SURFACE_LATTICE = HEX(auro)|FCC(a)`，入 config）。
- [ ] **B6**（增补，audit FINAL 行）删 FINAL.f:616-620 NZDOWN 写死原子 3 死分支
  （输入层 loud-stop 桩保留）。
- [ ] **验收**：regress 全 PASS——rst 类 case 轨迹除元素标签来源变化外数值一致；
  GWRITE 元素名与 ELEMENTS 名册一致（rst case 输入补 ELEMENTS 行后 diff 干净）。
  提交 `P-B: 测量层纯化（audit: GWRITE/TEST/GFINAL/FINAL/record_all）`。

### P-C 随机流统一（可复现性）

**Files:** 建 rng.f90；改 interface_TEST.f90、GAMA.f、VERLET.f90、VENUS.f90、GLO.f90、THERMO.f 等消费者。

- [ ] **C1** rng.f90（L5）：RAND0/RAND1/RANDST 原样迁入；gasdev 从
  interface_TEST.f90:112 迁入并改走主流（内部均匀源 RANDOM_NUMBER→RAND1，Box-Muller
  与缓存态保留）；GAMA 一并迁入（IDUM 形参恢复实效）。
- [ ] **C2** 全库唯一流核验：grep `RANDOM_NUMBER` 仅 rng.f90 内部转换点允许一处
  （若全改 RAND1 则零处）；所有消费者（VERLET 噪声、GLOSELECT、THERMO Andersen、
  QMMICRO 等）经 rng.f90。
- [ ] **C3** VENUS.f90:184,231 RANDST 双调用修复为读入后单次（F15 关闭）。
- [ ] **验收**：含噪声/恒温 case（rigid_surface、mb_thermal、rst_surface_md——执行时
  以各输入 NFC/NGLO/NTHERM 复核名单）同 ISEED 两次运行 bit 级一致；其余 case 与
  baseline bit 级一致。**预期基线更换**：上述三 case 基线快照更新，README 标注
  "P-C 换流，旧流不可复现"。提交 `P-C: 单一随机流（audit A2；F15 关闭）`。

### P-D 噪声律入接口 A（成员注册进 force_facade）

**前置**：N1 完成；P-B/P-C 已完成（常量先参数化、换流在前）。

**Files:** 建 pes_glo_spring.f90、pes_langevin_ghost.f90、pes_ldfa.f90；改 VERLET.f90、ENERGY.f、venus_input.f90、force_facade.f90（注册行）。

- [ ] **D1** VERLET.f90 的 ghost Langevin（:44-50,71-79,153-161）与 GLO 弹簧+随机+
  摩擦（:137-148）抽为成员：pes_langevin_ghost（随机律）、pes_glo_spring（弹簧+
  随机+摩擦）；原子索引经成员 ctx 参数化（不再写死 atom2/NATOMA+1/+2）；成员内部
  调 rng。注册行入 facade_init（由 NGLO/NFC 开关映射为成员装配，文法向后兼容）。
- [ ] **D2** ENERGY.f:40-47 GLO 弹簧势能删，能量侧经 facade 汇总（E=T(P)+ΣV，
  kcal/mol 报告口径不变：T/C1 后同标）。
- [ ] **D3** venus_input.f90:616-621 ghost 参数装配（W(2) 写死、413.4 换算）移入
  成员 init 钩子（facade_init 时调）。
- [ ] **D4** pes_ldfa.f90：以 FRICFORCE 物理（NFC=1/2 两律）重生为随机律第三成员。
- [ ] **D5** VERLET 残留 CALTYP 分支（:34-36,60-62,123-125）清除。
- [ ] **验收**：rigid_surface 能量分布与 P-C 后基线 bit 级一致（同流同算术序）；
  `grep -n "NGLO\|ghost" VERLET.f90` 零命中——积分器本体无体系知识。
  提交 `P-D: 噪声律成员化入接口 A（audit: VERLET/ENERGY/FRICFORCE）`。

### P-E 控制核清理（VENUS.f90）

**Files:** 改 VENUS.f90；接线 recorder/config_roster/predicate。

- [ ] **E1** 结局处理+统计块双份复制（:506-576 与 :581-662）合并为单一结局路径。
- [ ] **E2** 观测量积累（BOUNCE/z_pos_min 写死原子 3 :399-400,454-458；AE* 数组
  :544-548,631-635,671-677）移入 recorder 终态列；原子索引退为配置。
- [ ] **E3** 体系分支（GLO 布局 :322-335、表面对齐 :340、THTA 特例 :344-346）移交
  构型装配器（roster/装配期）。
- [ ] **E4** 空检查语句（:462-464）清除。
- [ ] **E5**（=△6）结局分类 if 链数据化：`OUTCOME = 名字: 谓词表达式` 行 × 计数桶
  （复用 N2 文法）；VENUS.f90 结局路径只查表计数。
- [ ] **验收**：regress 全 PASS；VENUS.f90 内 grep 无晶格常数/元素名/写死原子序；
  只剩装配→演化→统计编排 + 少量循环变量。提交 `P-E: 控制核纯化（audit: VENUS）`。

### P-F 全局状态分区（venus_data）

- [ ] **F1** venus_data 以注释显式三分区：物自体区（Q/P/PDOT/QDOT/Q_old——标注
  →state）/表象区（→config）/仪器区（→rng），每块行注释写明 X 轨去向。
- [ ] **F2** NN_POT_TYPES 索引（venus_params.f90:88-98）迁入 emt_nn 成员文件；
  auro 已随 B5 迁出，核实。
- [ ] **F3** 死块清理与 A4 联动复核。
- [ ] **验收**：`make` 通过 + regress 全 PASS。提交 `P-F: 全局状态分区标注`。

### P-G 输入文法与构型装配收尾

- [ ] **G1** test_potentials.f90 成员反写表象（BOXLX/BOXLY/SKEW/NATOMB/QZB，
  :152-163,260-293）改为经输入/构型装配器交付；PES 成员只留 Q→V/F 与数据读取。
- [ ] **G2** JMAXCALC.f:35-48 的 H2/D2 核自旋统计（质量判别）退为输入数据
  （`NUCLEAR_SPIN_WEIGHTS =` 表或数据文件）。
- [ ] **G3** BOXLX/BOXLY 默认几何（venus_input.f90:997-998）与 NATOMB 写死域
  （:703-704）改为必填（缺失 loud-stop 指名关键字）或名册推导；SURF.f:48-50 的
  SKEW=π/4 默认魔数同步必填化；venus_input.f90:752 的 W=1D30 刚性质量改为装配器
  按 RIGID 模型语义显式交付。
- [ ] **G4** RADAU.f 的 RADOUT 内嵌编排（:438-498）改为向控制核回调（输出钩子+
  谓词求值点，driver 侧接收）；'CALLING RADOUT'（:112,380,428）及 SELECT/SURF/
  ROTen/THERMO/RADAU/FINLNJ 调试打印残留一并清；PRINFO 约 70 行注释残留清；
  ROTEN :101-105 调试打印清；ROTATEJKM 尾部激光注释块清。
- [ ] **验收**：regress 全 PASS；EMT-NN/RST case 数值不变（成员内部合法性保留）。
  提交 `P-G: 输入文法与构型装配收尾（audit: test_potentials/JMAXCALC/RADAU）`。

## 7. 归位核对清单：原汇合轨 X

（2026-09-06 二次修订：X1 并入 K1（state/config 原生定义）；X2/X3/X4/X5/X6
并入 C7/C5/C6/C1/C8；X7 并入 P0/C7（新树构建与 EXTENDING）。本节保留作
V 验证前的新树完整性逐条核对。）

- [ ] **X1 state/config 实拆**：按 F1 标注执行。两步提交：先抽 state.f90
  （物自体 5 组+维数参数；全库 use 点机械替换），regress；再收拢 config.f90
  （表象块），regress。venus_params 的 C1..C8/PI 保留（单位契约，L1 常数）。
- [ ] **X2 driver 定型**：VENUS.f90 → driver.f90（三相位调用序列）；harmonic_sampling
  注册为 dist_surface_oscillator（B′ 第三实现，激活现状死编译件）。
- [ ] **X3 sampler 成员表**：SELECT 的 NACTA/NACTB 整数分发改 dist_* 成员表注册。
  名册：dist_boltzmann(THRMAN)、dist_ebk(INITEBK)、dist_surface_oscillator、
  dist_rotation(ROTEN)、dist_j(PROBJ+JMAXCALC)、dist_incident(SURF+DENS)、
  dist_barrier_excitation(BAREXC)、dist_thermostat(THERMO)、dist_glo_target(GLOSELECT)、
  dist_thermalize(md_equilibrate)、dist_normalmode(INITQP/HOMOQP/ORThAN +
  NMODE/FMTRX/FGMTRX[±ENMODE 按裁决]）、dist_qmmicro(QMMICRO)。顺带：SELECT
  :1139-1143 md_equilibrate 的 +100 Å 魔数参数化（默认值不变，行为不变）。
- [ ] **X4 propagator.f90**：积分器统一 init/step/输出回调协议；verlet/symple/radau
  三成员注册（G4 已剥编排）。
- [ ] **X5 L5 归位**：linalg.f90(EIGN)、geometry.f90(ROTATE 族/ROTATENO/ROTATEJKM/
  ROTATEX/Y/Z/CENMAS/ANGVEL/ROTN/ALIGN_ZSTEP/ALIGN_YSTEP/rotaline/GLPAR+PL)；
  rng.f90 已在（C1）；CPUSEC 入 rng 同级仪器位（VFDATE 已弃）。
- [ ] **X6 analysis/ 离线化**：FINAL/GFINAL 广义 COM/EREL/OAM 部分（结局分类执行器，
  消费 E5 数据）、EIGOUT、DENQ+SUMQ（DEL=1 cm⁻¹ 顺带配置化）、HIST 的 plot_* 五
  函数（格子常数传参化）、nearest-Au 诊断（B1 移来）、EIGOUT 的 IP=77/fort.77
  通道复用拆分；conserv_check、smoke_rst_bvk 纳入构建（scripts 告别游离命令行）。
- [ ] **X7 构建与文档**：Makefile SRC 按 L5→L1→L2→L3→L0 依赖序重排；
  EXTENDING.md 重写（扩展点=pes_* 成员/dist_* 成员/谓词原语/结局分类行/积分器注册）。
- [ ] **每任务验收**：`make` + regress 全 PASS（bit 级）。X7 后跑 T2 演练：新增
  `pes_leps` 变体成员（如改参数的 LEPSII）只动新文件+注册行，diff 确认既有文件零触及。

## 8. 收尾 S

- [ ] **S1** 新验证 case 体系全量执行收口 + 各 case results.txt 生成
  （新程序产出，含判据与容差记录）。
- [ ] **S2** 新程序完备性复审：按规范三问全例程审查，产出新 completeness 审计
  （md+json）；旧 audit 保留为历史参考。
- [ ] **S3** 文档收口：新树 CLAUDE.md/README/EXTENDING.md（扩展点=pes_* 成员/
  dist_* 成员/谓词原语/结局分类行/积分器注册）；design_spec_completeness.md §4
  留白填写（接口 A/B′ 形式签名 + 退化点良定义性检验 + 双实现验证记录：A 轴
  RST+EMT-NN+GLO 随机成员+参考势、B′ 轴 THRMAN+INITEBK+surface_oscillator、
  构型 twobody vs rst、电子轴绝热退化+TDHF-FSSH 留位）；optimal_structure.md
  状态改"已实现（新程序）"；根 README 更新布局（src/ 新程序 + venus2026/ 参考）。
- [ ] **S4** Notion "VENUS-TDHF 通用化"（阶段 P1 架构）任务对账：本计划各阶段
  勾选状态同步。

## 9. 验收量尺与回归矩阵

**验证矩阵**（2026-09-06 三次修订：旧 case 舍弃，验证体系重建）：

| 时点 | 范围 | 判据层级 |
|---|---|---|
| C1-C8（逐组件） | 每组件最小验证输入 + 判据（先定判据后写实现） | 解析锚点 > 统计分布（KS 等）> 守恒律/闭合量 |
| V | C 各步验证输入汇编为正式 case 集；判据覆盖率检查（每注册族成员、每封闭族成员至少一个 case 覆盖）；全量执行 + 同 ISEED 复现性抽验 | 同上 + 复现性（同 ISEED 两次运行 bit 级一致） |
| S | 全量收口 + results.txt | 同上 |

**结构量尺**（S 阶段逐条）：T1-T7（§1）；不变量三问抽检（规范 2.6）——抽
driver/force_facade/recorder/predicate/sampler 各一例程，问"是否只经手 Q/P/s/F、
是否按契约索服务、内部有无体系分支/常数"。

## 10. 风险与回退

| 风险 | 缓解 |
|---|---|
| 参考代码缺陷传染（audit 违规点带入新树） | §6 剥离清单逐条对照；check_* 组件级对照截获数值偏离 |
| 从零重写的行为漂移（同算法不同实现顺序→浮点差异） | 判据以解析/统计/守恒为主（对实现顺序不敏感）；需要逐位精度处（如换算表达式同形）在组件计划中显式声明；旧程序临时同量对照辅助归因 |
| 骨架与实际不符（C 中发现签名不适用） | 骨架漂移禁令：回审核门出示差异，经确认后修订骨架，禁止静默改签名 |
| 新旧混淆（venus2026/ 参考 vs src/ 新树） | 旧库只读；新树自 P0 起独立构建；新树每文件头注明参考出处（旧库 文件:行） |
| 验证判据覆盖不足（旧 case 舍弃后的新风险） | 判据先行原则：C 各步先定判据再写实现；V 阶段判据覆盖率检查（注册族/封闭族成员逐一覆盖）；缺口即补 case |
| C5 sampler 成员最多、面最大 | 成员表化对照 NACTA/NACTB 分发表逐支核验；integrator_matrix 等 case 覆盖多分支 |
| Windows/Linux 双平台构建 | 新树 Makefile 双平台自 K3 起持续可编；V 末双平台各跑一次（如可行） |
| 回退 | 每小 phase 独立计划文件 + 独立提交组；`git revert` 单 phase 粒度恢复 |

## 11. 交付物清单（完成标志）

1. 独立新程序全树（src/，L0-L5，T1）+ 新验证 case 体系全量通过（V/S1）；
2. 新程序完备性审计报告（S2）；
3. 文档收口全套（S3）：spec §4 填写、新树 CLAUDE/README/EXTENDING、
   optimal_structure 状态更新、根 README 布局；
4. 过程留痕：docs/plans/ 每小 phase 一份计划文件 + git 提交组，全程可回溯。
