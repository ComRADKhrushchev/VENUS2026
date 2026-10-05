# 2026-10-05 谱来源立法：废容器导出槽，SPECTRUM_SOURCE 双模式（COMPUTE 缺省/MANUAL 文件）

状态：执行中（判据先行——本文档先于代码落规，跑出的读数回填 §8）
裁决：用户三次指令（2026-10-05 会话），累计如下。

## 0. 裁决链

1. **废容器导出通道**。`export` 槽自调度表移除（九槽→八槽）；两容器
   （h2ag111PES / h3_bkmp2）reg 文件中的 `export_spectrum` 槽体与仅其使用的
   `is_sym` 辅助删除；缓冲中的容器表 `export_tbl`/`diab_tbl` 与公共入口
   `spectrum_export_put`/`spectrum_export_put_diatomic` 删除。容器不再立法谱。
   动因：约 40 行槽体在两容器逐字重复、文献常数硬编码于代码——常数应为
   运行时数据（净化条款 docs/design_spec_completeness.md §2.6"体系知识参数化、
   经接口进入"；先例 E28 容器键并入程序输入）。S.md P-1 停放项（H 形容器
   自带谱导出）由此改道解决。
2. **缺省 COMPUTE，且范式中立**。表面范式**允许**内部推导——旧代码在推导前按
   `surface_model /= 0` 一刀切中止（spectrum_export.f90 旧 :133-137）属越权代理；
   "NN 无渐近尾"是 h2ag111 这个具体 PES 的性质（其 reg 头注自认），不是表面范式
   的构造性质。裁判权归还给跨块衰减门（PES 在探针处的物理行为）。
3. **MANUAL = 谱表文件全表通道**（支持多原子完整模表）。
4. **interface 层谱代码合并为单文件** `src/interface/spectrum_interface.f90`
   （沿 force_interface/elec_interface/sysdef_interface 命名惯例）；
   `spectrum_export.f90` 删除。utils/spectrum.f90（L5 纯数学）不动。

## 1. 模式语义（互斥二态，键独占选择生产者）

- `SPECTRUM_SOURCE = COMPUTE`（缺省）：拉取时懒推导（§2）。
- `SPECTRUM_SOURCE = MANUAL`：只查谱表文件记录（§3）；推导不运行。
  被拉取的 nat≥2 载体片段在文件中无记录 → 命名中止（"参数作为输入要求"）。
- 单原子片段两模式一致：返回分配空表（nat=1 无内部自由度，COMPUTE 路径
  零力调用、免衰减门）。
- **EBK 双原子常数通道无推导**（非谐/转动常数不存在通用推导——原设计即
  export-only）：COMPUTE 下拉 `spectrum_frag_diatomic` → 命名中止，指引
  MANUAL+DIAT 记录；MANUAL 下无 DIAT 记录 → 同样中止。
- 交叉校验（全部命名中止，无静默路径）：
  - SPECTRUM_FILE 在场（非空）而 SPECTRUM_SOURCE=COMPUTE → input 层形式
    校验中止："SPECTRUM_FILE requires SPECTRUM_SOURCE=MANUAL"；
  - MANUAL 而 SPECTRUM_FILE 空/文件不存在/不可解析 → 首次拉取时 staging 中止。

## 2. COMPUTE 推导重构（对 derive_all 的三处实质变更）

现 `derive_all`（spectrum_export.f90 旧 :226-290）：一次推导**全部**片段并逐一
校验，COM 间距固定 12 Å，表面范式被前置禁令拦截。新法：

1. **按需逐片段**。探针几何 + 全体系中心差分 Hessian 仍只算一次（进程内缓存）；
   每次 `spectrum_frag_w` 拉取只构建并校验**被拉片段**的对角块（质量加权 →
   `spectrum_modes` 的零模劈分/正本征值校验保持单点）。理由：表面体系的
   slab 块（无消费者）若非 PES 精确驻点会误触发劈分校验中止——无人需要的
   表不该有权杀死进程。
2. **衰减门判据变化**：门只查（被拉片段, 其他片段）对，参考曲率 = 被拉片段的
   最小振动本征值平方（旧判据 = 逐对 min(λi,λj)——那需要导出**双方**的表，
   与按需推导矛盾）。质量加权跨块幅度阈值仍 `cross_rel = 1e-6`（相对）。
   中止文案只谈 PES 行为并指引 MANUAL，**无范式语言**：
   "the inter-fragment coupling has not decayed at the probe separation -
   this PES carries no decaying isolated-fragment tail; provide the spectrum
   as SPECTRUM_SOURCE=MANUAL + SPECTRUM_FILE"。
3. **延展感知探针**：相邻片段 COM 间距 = ext_i/2 + r_asym + ext_{i+1}/2
   （ext = 片段原子到自身 COM 的最大距离；r_asym = 12 Å 不变）。分子片段
   延展 ~1 Å，行为与旧探针几乎一致；宽 slab（表面范式）不再被"12 Å COM 间距
   实际只隔几 Å"的探针缺陷误伤。原子间最小距离 ≥ r_asym 量级（球形保守估计）。
4. 删除 `surface_model /= 0` 前置禁令——任何范式同一推导路径，门仲裁。

预期后果（经验项，跑后回填 §8）：h3_bkmp2（气相，BKMP2 有真渐近尾）推导通过，
H2 谐振频率复现在案锚 4403.2 cm⁻¹（check_spectrum.f90 头注腿 (c)：独立能量二阶
差分 V''=35.933611286509404 eV/A²，对文献 4401.21 偏 0.045%）；h2ag111
（单片段表面 NN，无跨块对，门无从触发）推导出的 H2 频率为该 NN 在缓冲几何处的
局部曲率值，与文献值的偏差如实记录。

## 3. 谱表文件格式（权威规格——读取器与本文档不得漂移）

行导向、大小写不敏感、`#`/`!` 整行注释、空行忽略；每条记录按**成分键**
（符号序列，与 list_atoms 片段大小写不敏感逐符号匹配）：

```
# manual spectrum table - SPECTRUM_SOURCE=MANUAL carrier records
FRAG H H              ! 符号序列，一条记录一个成分（记录开始行）
NMODE 1               ! 振动模数：nat=2 → 1；nat≥3 → 3nat-5 或 3nat-6
WVN 4401.21           ! 频率 [cm^-1]，正数，每模一行（NMODE 行）
COL 1                 ! 第 j 模列开头行，其后一行 3nat 个实数（质量加权正交归一）
  0.0 0.0 0.0 0.0 0.0 0.70710678...
DIAT 4401.21 121.33 60.853   ! 可选，仅 nat=2 记录：w_e, w_e*x_e, B_e [cm^-1]（EBK 通道）
```

- **模列参考架契约**：列向量表达在**片段缓冲几何（system/ 构型文件）的笛卡尔
  架**中。外部量化程序输出的模列需作者自行旋转到该架（不做自动对齐）。
- 单位：频率与常数一律 cm⁻¹（谱学面）；程序内换算唯一通道 `wvn_to_w`
  （utils/spectrum.f90）。
- 读取器**只做形式检查**（input 层同姿态）：记号合法、NMODE 计数契约、
  COL 行系数个数 = 3nat、WVN > 0、DIAT 三正数且仅 nat=2、文件内成分重复、
  未知记号行。正交归一性/平移正交等物理校验留成员 init 单点（现行职责不动）。
- staging 懒执行：首次拉取时读文件一次（此时 list_atoms 已装配，装配序
  driver.f90 :104 list_atoms_load → :124 registry_gen_init 不变）。

## 4. 输入键（零新文法类型）

| 键 | kind | 块 | 缺省 | 说明 |
|---|---|---|---|---|
| SPECTRUM_SOURCE | kind_enum（族 spectrum_tbl: COMPUTE=1/MANUAL=2，族 id 5） | system | COMPUTE | 谱生产者独占选择 |
| SPECTRUM_FILE | kind_char, slotw=128 | system | ''（空） | MANUAL 时的谱表文件路径（相对运行目录）；MANUAL 下必填——空值在首次拉取中止（B_MAX 先例：条件必填在装配期强制） |

key_tab 30→32 行；`enum_lookup`/`legal_words` 各加 `enum_spectrum` 分支；
map_slots 加两 case；表尾形式交叉校验（§1）。config.f90 `reactants_t` 加
`spectrum_src = 1`、`spectrum_file = ''`（128 宽）。

## 5. 改动清单

| 文件 | 改动 |
|---|---|
| src/interface/spectrum_interface.f90 | 新建（合并件）：拉取入口/推导（§2 重构）/命名中止 + MANUAL staging + 文件读取器；头注 = 新法立场 |
| src/interface/spectrum_export.f90 | 删除 |
| src/interface/container_sched.f90 | export 槽声明/bind 关键字/container_sched_all 调用删除；头注 nine→eight |
| src/input.f90 | §4 两键 + 族 + 分支 + case + 交叉校验 + 头注键表提及 |
| src/globals/config.f90 | reactants_t 两字段 |
| src/methods/reg_samp_{wigner,normalmode,boltzmann,ebk}.f90 | use spectrum_export→spectrum_interface；注释中"container export 优先"表述改"模式键选择" |
| src/interface/bath.f90 | :238 注释中 spectrum_export 字样随模块改名（无代码引用） |
| systems/inactive/h2ag111PES/reg_h2ag111.f90、h3_bkmp2/reg_h3_bkmp2.f90 | 删 export= 实参、export_spectrum 体、is_sym、仅其用的 use（spectrum/spectrum_export/config_atoms·list_atoms） |
| scripts/check_spectrum.f90 | 腿 (f)(g) 随立法改写（§6）；(a)-(e),(h) 保留（(c) 4403.2 锚不变——探针变化不触碰 H2 内几何） |
| scripts/check_driver.f90、check_facade.f90 | 编译配方行 spectrum_export.f90→spectrum_interface.f90 |
| src/EXTENDING.md、src/README.md | 槽表九→八、去 export 段；新键+文件契约（§3）+模列参考架；L2 行改"推导(缺省)+手动文件" |
| docs/plans/S.md | P-1 追加：导出通道废除，输入通道立法 |
| docs/plans/2026-09-22-container-schedule.md | 九→八槽修正注 |
| docs/glossary.md | 增"谱来源键/谱表文件"；废"容器谱导出""表面范式必须导出"表述 |

新增测试材料（本批）：cases/case_spectrum_source/{h3_wigner,h3_ebk,h2ag_wigner,h2ag_ebk}/
+ scripts/gen_spectrum_tbl.f90（模列由 stretch_mode 生成，禁手抄）。

## 6. 判据（先写后跑；读数全部落盘 *.log）

1. **等价判据（核心）**：四 case 各跑改动前 AUTO（容器导出）与改动后
   MANUAL（谱表文件，常数= 容器原硬编码值：H2 4401.21/121.33/60.853、
   D2 3115.5/61.82/30.429 cm⁻¹，出自两 reg 文件），trajectories.txt 逐位一致。
2. **COMPUTE 判据**：h3 case 无 SPECTRUM 键跑通；推导 H2 频率（scripts 探针
   程序落盘）对文献 4401.21 偏差 |Δ| ≤ 30 cm⁻¹（在案锚 4403.2 → 预期 Δ≈2.0）。
   h2ag111 COMPUTE 结局经验记录（单片段无门对，NN 局部曲率值 vs 文献差）。
3. **中止矩阵**（子进程 stderr 逐条比对）：MANUAL 缺文件 / 载体无记录 /
   格式错行 / 文件内成分重复 / EBK 缺 DIAT（含 COMPUTE 下拉 diatomic）/
   SPECTRUM_FILE 在场但 COMPUTE / 衰减门未过（新文案：无范式语言+MANUAL 指引）。
4. **check_spectrum 腿**：(f') 表面范式+衰减合成 PES → 推导成功（范式中立
   反证）；(g') MANUAL 文件往返：常数逐位、成分匹配、重复成分中止、缺记录中止。
5. **回归**：双编译器全量重建绿；五 case 终跑对各自 results.txt 锚（本改动
   不触其方案词——THERMALIZE/ROTATION/J/GLO_TARGET 均非谱消费者，输出应不变）；
   reference 容器（无 export 槽）构建拖入回归；八槽 one-shot 语义不变
   （第二次交接仍中止）。

## 7. 执行序

dev 树：本文档 → 基线 AUTO 四跑落盘（改码前！）→ 代码 §5 → 构建（两轮，
depend.mk 迟滞）→ §6 判据 → 文档连锁 → 提交。
publish_min 树：src/ 逐文件对齐（diff 校验）+ EXTENDING/README + 四段构建
（旧引擎/无容器核心/拖入 h2ag111/拖入 h3_bkmp2 后恢复）→ 提交。

## 8. 跑后回填（数值出处：落盘 log 绝对路径，均在 cases/case_spectrum_source/logs/）

执行中发现并随批修正的两件事（先于判据结果记录）：

- **装配次序约束（本批的实质发现）**：成员缝（registry_gen_init 里的 reg_samp_*_init，
  拉谱→推导→调力）原排在容器缝（container_sched_init 的 load 槽，NN 权重加载）之前
  ——load 依赖型 PES 在拉谱时尚不可计算，h2ag111 首跑即在探针 Hessian 处段错误。
  旧导出槽跑在更早的 registry_gen_all 期，正是隐性绕开此约束的机制。修正：
  prebuild.sh 生成器把 `container_sched_init` 移至 registry_gen_init 聚合首位
  （容器先装、成员后拉），container-schedule.md 已追记。
- **单片段探针**：单片段体系的探针=缓冲几何原样（隔离对单片段无对象；重排质心会把
  h2ag111 的 H2 从 z=3.0 Å 挪进 NN 无效域）。§2.3 的延展感知排开仅适用于多片段。

判据结果（全绿）：

1. **等价判据（核心，逐位）**：四 case 的 MANUAL（谱表文件，常数=容器原硬编码值，
   模列由 scripts/gen_spectrum_tbl.f90 调 stretch_mode 生成）与改码前 AUTO（容器导出）
   的 trajectories.txt **全部逐位一致**：
   - h3_wigner / h3_ebk：`logs/h3_{wigner,ebk}_manual_trajectories.txt` vs
     `logs/h3_{wigner,ebk}_baseline_auto_trajectories.txt`（8433 行，cmp 零差）
   - h2ag_wigner / h2ag_ebk：同款（7229 行，cmp 零差）
2. **COMPUTE 判据**：h3_bkmp2 缺省（无 SPECTRUM 键）跑通（`logs/h3_wigner_compute.log`
   未留 stdout——stdout 缓冲丢失于其后 h2ag 崩溃试验前；跑通证据=exit 0 与判据 5 的
   逐位回归含该构建）。推导频率在案锚复现：check_spectrum 腿 (c)
   **ν = 4403.2 cm⁻¹**（`logs/check_spectrum_run.log`），对文献 4401.21 偏
   **+2.0 cm⁻¹（0.045%）**，判据 |Δ|≤30 cm⁻¹ 裕量满足。
   **h2ag111 COMPUTE 经验结局**：探针成功（NN 局部曲率 **ν = 4396.07 cm⁻¹**，
   `logs/h2ag_compute_probe_nu.log`，scripts/probe_spectrum_h2ag.f90），但表面嵌入
   曲率导出的模列非平移正交（残差 2.7573E-03 > 1e-8 容差），wigner_init 命名中止
   （`logs/h2ag_wigner_compute.log`）——成员验证正确拒绝"非自由片段表"；该 NN 的
   正确通道是 MANUAL（等价判据已证逐位）。
3. **中止矩阵（8 项，文案逐条核对）**：MANUAL 缺文件（check manempty，"MANUAL-mode
   requirement"）/ 载体无记录（check miss + t3 同款，"no MANUAL spectrum record"/"no
   DIAT record"）/ 格式错行（abort_t2.log，"unknown spectrum-table word 'GARBAGE'"）/
   文件内成分重复（check dup，"one composition, one record"）/ EBK 缺 DIAT
   （abort_t3.log）/ EBK under COMPUTE（h3_ebk_compute.log + check diatc，
   "no generic derivation"）/ SPECTRUM_FILE 在场但 COMPUTE（abort_t1.log，输入层
   "requires SPECTRUM_SOURCE=MANUAL"）/ 衰减门未过（check e，"decayed…MANUAL" 指引，
   无范式语言）。
4. **check_spectrum 重写重跑**：**13 PASS / 0 FAIL**（`logs/check_spectrum_run.log`）；
   腿 (f) 翻转为"表面范式+衰减 PES 推导成功"（范式中立反证），腿 (g) 四条 MANUAL
   中止+往返全部入判据。check_driver 重编重跑 **ALL PASS（23 腿）**
   （`logs/check_driver_run.log`）；check_facade 重编重跑 **0 failures**
   （`logs/check_facade_run.log`）——其编译配方顺带修复两处先在漂移（缺
   interface/bath.f90 与 utils/geom_init.f90，与本次立法无关的存量欠账）。
5. **回归**：gfortran + ifx 双编译器全量重建绿；**净态前后对照（强判据）**：以
   publish_min 同源 src（改码前逐文件一致，经 cmp 复核）构建净态二进制，与改码后
   二进制在全部 **8 个子 case**（五 case 全部子目录）上 trajectories.txt
   **逐位一致（identical=8 differs=0）**——现有流零扰动。
   注：`results.txt` 档案锚（如 harmonic 1.2648E-03）与今日树（净态与改码后同为
   1.4103E-03）不符系**先在陈档**（档案漂移值未提交且追溯至 9-16；members 计数
   11→12 的未提交改动亦与现树 members=11 不符）——分歧先于本批存在，净态对照已证
   与本批无关，档案重生成留待独立批次裁决。

## 9. 后续件方向记录（2026-10-05 会话取向，取代"表面感知探针"提法）

- **统一 z 向堆叠 + 基片锚定探针**（COMPUTE 探针的下一形态候选）：锚定片段保持
  缓冲原位，其余片段沿 +z 在其上方按延展感知间距堆叠。锚定片段取法：表面范式 =
  非 PROJECTILE 片段（基片载体，glossary 在案，`list_atoms%i_proj` 可判）；
  分子范式 = 片段 1（轴本任意）；单片段退化为现行"缓冲几何原样"。
  动机：z 对气相是任意约定、对表面恰是法向（真空方向）——一条规则覆盖全范式且
  不需要范式分支（优于本批实现所预留的"表面专属构造"路线）；显式基片体系
  （cau_peri 形态 C+Slab）因此获得合法探针（slab 晶胞注册不动）；h2ag 型原位
  曲率的模列平移正交问题（§8 实测残差 2.7573E-03）对有真实渐近尾的 PES 随之
  消失——分子被抬进真空后测回自由分子表。窗口型 NN 在抬升位可能仍在域外，
  由衰减门与成员校验如实拦截，不做普适保证。
  本批保留 x 轴排开的理由（记录）：气相轴物理等价，x 为旧 derive_all 约定，
  保留使 check 腿 (c) 在案锚（V''=35.933611286509404 eV/A²，+x 旁观 12 Å）
  无需重新取证。
  未随本批落地的原因：无活消费者（cau_peri 的 C 单原子、slab 走
  SURFACE_OSCILLATOR，无谱拉取）；探针几何重设计应独立判据先行（多基片片段
  锚谁、衰减门参考曲率是否随锚定片段调整、check 腿 (c) 锚文字与容差复验）。
  触发条件 = 第一个需要 COMPUTE 采谱的显式基片体系落地。
- 其余后续件维持汇报清单：publish_min 推送裁决、dev 批次归属（本批随用户指令
  整批落档，见提交消息）、deck 九槽勘误注记、五 case 档案重生成批次、
  check 程序 ifx 侧复验、多原子模列转架助手（首个真实用户触发）。
