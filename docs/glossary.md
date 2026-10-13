# 术语与代号词汇表（唯一权威）

状态：v1，2026-09-17 随术语规范化专项建立（方案与执行记录：docs/plans/2026-09-17-terminology-normalization.md）。
地位：行文称谓与代号命名空间的唯一登记处。新代号入册两要件：**定义处在册 + 跨文档引用带限定**。每阶段收尾固定做一次术语同步（新增/退役词条入册）。
语言规约：中文行文用中文规范形，英文行文用英文规范形（AGENTS.md "match the surrounding language"）；同一概念每语言一个规范形。

---

## 1. 称谓规范表（2026-09-17 用户决定）

### 1.1 替换定名（十一条）

| 概念 | 规范中文 | 规范英文 | 弃用形 |
|---|---|---|---|
| systems/<sys>/ 文件夹实体 | **物理系统容器**（重复处可简"系统容器"） | physical-system container | 单字"容器"（行文不再用）、blackbox/黑盒（史层） |
| 接口族模块类称（原名 facade 族，2026-09-22 决定改名） | **接口**（力接口/电子接口/体系定义接口/采样接口） | 以模块名直称：force_interface / elec_interface / sysdef_interface / sampler（目录 interface/） | 门面、facade（史层原词不改） |
| 分发行→成员 init 的通道 | **初始化组装机构** | initialization assembly | 装配缝 |
| 容器交付的受力计算（container_pes_i 全赋值 g,v） | **受力（计算/交付）** | force evaluation | 组合核、composed force kernel |
| roster 类称（标识符 2026-09-22 决定改名，覆盖同表"标识符 config_roster 不变"旧注） | **（具体）列表**：体系组成列表/片段列表 | atom list；标识符 list_atoms / list_atoms_t / list_atoms_load，模块 config_atoms | 名册、roster（史层原词不改） |
| 解析后暂存待组装期消费的中间存储（原 staging，2026-09-22 决定改名） | **缓冲（区）** | buffer；标识符 buffer_*（input 行缓冲 / sysdef 交付缓冲） | staging / staged（史层原词不改） |
| driver 积分循环 | **演化循环** | evolution loop | 演化环、"四相位体/五入口"计数名 |
| 按决定固定取值的换算/物理常数 | **物理常数** | physical constants | 立法常数 |
| 项目大阶段 | **阶段** | stage / phase | 相（"C5 相"→"C5 阶段"；程序三阶段本就称阶段） |
| 同初值验证性重跑 | **两次运行** | two runs | 双跑 |
| 用户/架构裁定 | **决定** | decision | 裁决、裁定、重裁（史层原词不改） |
| 骨架→真实现 | **实现** | implement | 成真 |
| 弹丸/靶的显式输入声明（PROJECTILE/TARGET，2026-09-21 决定） | **角色键**（弹丸键/靶键） | role keys | "片段 1=弹丸"位置约定（史层；现行无位置缺省，多片段体系必填） |
| 容器参数的输入通道（E28，2026-09-22 决定） | **容器键**（声明键） | declared container keys | 运行目录 *_pars.txt 旁通道（史层，已退役）；现行：box 键表 + reg `declare_keys` 声明 + `load_params` pull + 消费审计 |
| 容器 reg 的通用时序载体（2026-09-22 决定；2026-10-05 八槽修订） | **装配时刻表**（八接线槽） | container schedule；模块 container_sched，容器侧唯一入口 reg_<sys>_slots，四个接线区块 declare_keys/load_params·bind_force/register_method·bind_term/arm_method·bind_classify/register_columns（注册期谱导出槽已废除——2026-10-05 谱来源立法） | 容器 reg 自持 _keys/_all/_init 三入口（史层形态，已退役）；九槽形态（2026-09-22 至 2026-10-05，含 export_spectrum 槽，史层） |
| 谱生产者的输入选择键（2026-10-05 立法） | **谱来源键** | SPECTRUM_SOURCE = COMPUTE（缺省，内部 Hessian 推导，范式中立）/ MANUAL（谱表文件唯一生产者） | 容器导出通道（注册期 export_spectrum 槽 + spectrum_export 缓冲，史层，已废除）；EBK 双原子常数无推导通道，仅 MANUAL DIAT 记录 |
| MANUAL 模式的谱数据文件（2026-10-05 立法） | **谱表文件** | spectrum table；SPECTRUM_FILE 指向，行导向格式（FRAG 成分键/NMODE/WVN[cm⁻¹]/COL 模列/DIAT 双原子常数），列表达在片段缓冲几何笛卡尔架 | 权威规格 docs/plans/2026-10-05-spectrum-source.md §3；生成器范本 scripts/gen_spectrum_tbl.f90 |
| 同方案词选中的片段（成员视角） | **载体（片段）** | carrier (fragment) | "单载体规则"（2026-09-21 决定后同词多片段合法） |
| 一份成员参数服务全部载体 | **平等处理** | equal treatment | 按侧面（A/B 面）参数差异化——登记为后续件 |
| 采样后对片段整体的各向同性取向抽取（RANDOM_ORIENT=T 时启用，2026-09-22 决定） | **取向随机化** | orientation randomization (stage) | rotation 成员内部的同型操作称固定 J 的取向实现（其载体被普适阶段跳过） |
| 分布成员的两段式接口：draw 消耗随机流写成员私有记录、realize 写 q/p 不耗随机流（2026-09-22 方案 A） | **抽签/落地** | draw/realize | 融合包装 xxx_sample 保留为判据程序入口；落地顺序=名册序是兼容立法 |
| 入射通道三成员 | incident channel splits per assembly paradigm: inc_surface (incident_surface) / inc_pair (incident_pair) / inc_single (incident_single); samp_incident_word() resolves, paradigm-guarded registration keeps members=11 | 入射通道按范式拆分 | inc_surface/inc_pair/inc_single | 2026-09-28 裁决 docs/plans/2026-09-28-incident-split-surface-paradigm.md |
| 定点落点 | N_AIM=2 + AIM_X/AIM_Y [Angstrom, POSCAR frame] - the surface fixed-incidence option (b offset stacks on top; B_MAX=0 = on-site beam under N_B=0) | 定点落点 | AIM_X/AIM_Y | surface-only mode; criterion check_inc_surface (v) |
| 基片片段 | substrate carrier: a surface-paradigm fragment other than the PROJECTILE fragment - its own DIST_SCHEME member samples it; the beam assembly never re-places or rotates it (crystal-frame exemption) | 基片片段 | substrate carrier | the relaxed single-fragment statute (2026-09-28); cau_peri C.xyz + Slab.poscar is the in-repo form |
| samp_ 前缀族 | the per-fragment initial-state distribution family file/module prefix (renamed from dist_ 2026-09-28; registration words unchanged) | 采样分布成员前缀 | samp_*/reg_samp_*/check_samp_* | the incident family uses inc_*/beam_laws |
| 演化循环前的碰撞平衡机制（BATH=ANDERSEN 时启用，thermostat 摘出 dist 家族 2026-09-22） | **浴循环** | bath loop | species=bath 物种（模块 bath_andersen，Andersen 重置）；循环含环内抽签的明文豁免 |
| BATH 激活即隐含：无主原子（基板）从 list_atoms%q0 就位；有方案的片段 q/p 逐位冻结并经 geom_sep_probe 摆到无相互作用距离（兼作 PES 健康检查） | **浴前置播种 / 冻结立法** | bath seeding / freeze statute | 探针度量=片段对环境受力的扰动（全局场在差分中消去） |
| Tier 1 纯势容器自动接线识别的十个 public 过程名（pes/keys/load/term/classify + prop/match/sync/init/method_name；术语定名与方法包初始化名 arm→init 均 2026-10-13 用户决定——名称须自释初始化义） | **预留名** | reserved names | 魔数名（史层不改：裁决记录 docs/plans/2026-10-12-auto-pes-container.md 原词保留）；arm/武装（史层：arm→init 更名前的槽位与预留名，docs/plans 与既有裁决记录原词不改） |

### 1.2 已清除的特殊说法（用常规陈述代替）

| 已清除 | 常规表达 |
|---|---|
| 停放（停放项/停放区） | "已登记、暂不处理（的事项）"；P-1..P-7 行文称"暂缓事项 P-x"（S.md 表格原文不改） |
| 空转 | 无操作（no-op） |
| 组合根 | 陈述式："物理系统容器是体系专属代码（H 定义、体系专属方法包）的组合点"；英文 composition point |
| 锚 / 锚点 | R8 文献语境=**参考值**（锚点表→参考值表、读图锚→读图取值）；case 语境=**判据（编号）** |
| loud-stop | 中文=**报错即停（带明确错误信息）**；英文注释=named abort / abort with a named error；**运行期消息前缀=fatal:**（批 C 已改 66 处） |
| 具名停机 / 点名 | 随 loud-stop 清除："错误消息指明责任键" |
| class-A / class-B（设计期成员分类代号，2026-09-18 用户决定全退役） | 受力侧成员=**受力成员**（force members，Interface A 行）；电子侧成员=**电子方法成员**（electronic method member）。源码注释与现行文档不再使用；史层（docs/plans、final_completion_plan 伪代码）原词不改 |

保留的常规词（不清除）：骨架（skeleton）、装配（assembly）、判据（criterion）、成员（member）、注册（registry）、安全网（safety net，指 MAX_STEPS 兜底时）。

### 1.3 流程动词约定

| 规范形 | 用法 | 分工 |
|---|---|---|
| 决定 | ruling/decision | 用户/架构裁定事项 |
| 退役 | retire | 已落地机制的退出（撤编并入此） |
| 撤销 | withdraw | 计划项未落地即取消 |
| 实现 | implement | 骨架→真实现 |

### 1.4 英文注释层

清除真自造词（loud-stop→named abort、composition root→composition point 陈述、idle-chain posture→no container bound、safety-net-only run→MAX_STEPS-only run、legislated→fixed by decision、statute→rule）；保留业界常规词（container、facade、member、registry、seam——seam 与标识符 seam_cast 同族）。
**leg/腿**（2026-09-24 用户指示）：行文层不再使用，用平实表述（"计算方案/该组计算/曲线"）；
代码与判据层同名标识符（probe_chain_freq 的 `leg`、判据条目标签）仍按"定义处保留"不动，
全面清理待后续授权批。

## 2. 遗产词汇面（不改名，只立档）

| 词 | 含义 | 备注 |
|---|---|---|
| GWRITE_LEVEL | 记录详细级（0/1/2/5），裸缺省 2 | 旧库 GWRITE 例程名的键名遗产 |
| sec_qp / sec_frc | 记录器分节：section q,p / section forces | recorder.f90 已补 "sec_ = section" 展开注 |
| input_qct.txt | 硬编码输入名 | QCT = quasi-classical trajectory（driver.f90 已补展开注） |
| thta | theta 的旧库拼写变体 | 输出 hist_thta.dat 沿用 |
| MQC | mixed quantum-classical | elec_facade.f90 已补展开注 |
| wvn | wavenumber | consts.f90 已补展开注 |
| e_conv / ev_kcal / e2wvn / wvn2e / hbar_code / kb_code | 物理常数标识符 | consts.f90；**守护句：取值按既有决定固定（固定换算对 23.0605×0.04184=0.96485132），检查程序逐字使用，不得按 CODATA 重推** |
| venus.e / results.txt / stdout_* | 可执行名 / 回归对账基准 / case 存档 | stdout_official/rerun/run + stderr_official/rerun |
| expect_loudstop | 判据程序的子进程断言助手（标识符） | 断言退出码≠0 + stderr 含责任入口名与责任键名；不依赖 fatal 前缀 |
| fort.* | 旧库产物通道名 | 仅 venus2026/（只读旧树）使用 |

| BOXPAD / 簇盒投放 | 已退役（E29，2026-09-22）：表面窗口键 BOXPAD 与 incident 通道 N_AIM=2 簇盒投放同批删除；表面窗口键唯余 A_LAT | retired key |

## 3. 代号命名空间登记

### 3.1 现存各族（属主与引用限定式）

规则：**每个代号族有唯一属主文档；族内裸用仅限属主文档；跨文档引用必须带限定。**

| 族 | 属主 | 限定式示例 |
|---|---|---|
| E1–E28（决定号） | docs/plans/K2.md（E1–E26）+ S.md（E27）+ T.md（E28） | 裸用特权族；与 P-E 任务 E1–E5 歧义时写"P-E 任务 E1" |
| 阶段 P0/K1–K3/C1–C8/V/S/R8 | docs/plans 文件名 | "C5 阶段" |
| 总案 §5 接口规格 N1–N3 | final_completion_plan §5 | 总案 N2 |
| 审计登记注记 N1–N7 | completeness_audit_src §2 | audit N1 |
| 审计架构发现 A1–A8 | completeness_audit §2 | audit A2 |
| 暂缓事项 P-1..P-7 | S.md 终局表 | S 暂缓 P-3 |
| V 阶段发现 V-1..V-6 | V.md 发现清单 | V-4 |
| 结构判据 T1–T7 | final_completion_plan §1 | 判据 T4（第二个裸用特权族） |
| 计划内决定 R/S/D 系 | 各 plan 自身 | C5-R1、C6 决定 H、C2-S3；fh2 决定点引用写 R8:D5 |
| 判据条目号 (a)–(s) | 各 check 程序头注 | check_facade (p)；裸用仅限该程序与其 results 存档 |
| case 判据编号 | 各 cases/*/README（已加范围声明行） | case_b_layer T2；README 内部可裸用 |
| 旧库常数 C1..C8 | K1 §1.1 对照表 | venus_params C1 |
| R8 子系（R8-1…、R8-A…、M-R8a…、[P1]–[P9]、A0–A2/B1–B3/O1–O10、T-三分容差） | R8.md / fh2_leps_anchors.md（已加范围声明） | 参考值 A1(R8-1) |
| 旧例程专名（POT0/DPESHON/GWRITE…） | completeness_audit / optimal_structure 映射表 | 仅作出处标签；check_facade (s) 旧链编译清单照用 |

### 3.2 未来规则

1. 全局决定续 E28+ 还是封线——**待决定**；若续号，E28+ 统一记入新台账 docs/plans/E_ledger.md（阶段文件只留链接）。
2. 新阶段/新顶层文档用描述性文件名（日期-主题式；先例：本专项与 superpowers plans），不再启用新字母。
3. 计划内决定一律 `<plan>-` 前缀。
4. 判据条目号一程序内不复用已退役字母（check_facade 的洞保留——"保留历史字母"）；新增子条目用数字。
5. case 判据编号选字母约定（成员首字母）：T=thermalize、H=thermostat、G=glo_target、F=频率结构、B=boltzmann、N=normalmode。
6. 撇号：文档语境 U+2032（B′），代码/标识符语境 ASCII（B'）——**字符规范本体待决定**（现状已按此统一）。
7. 源码注释规约：零决定号（E 号）/零任务阶段标号；判据头注自含（"criteria defined FIRST in this header"）。现存例外（在案登记）：input.f90 与 elec_facade.f90 的 5 处字符串字面量（E21/E22/E13，报错消息正文）；check_* 头注/体内的判据出处注（E20/E21/E22/E26/E27、E7 条目号等）——判据程序头注作为"判据先行的定义处"允许引决定号，生产代码 src/ 注释已零决定号。

## 4. 孤码登记表（定义已佚，含义在案）

| 孤码 | 含义（抢救自上下文） | 出处线索 |
|---|---|---|
| F14 | RANDST 双调用缺陷 | C5.md 附录、final plan §6；总表随 remediation_plan.md（2026-09-06 删）佚失 |
| F15 | 双随机流缺陷 | 同上 |
| F19 | 振动角动量剥离 | 同上 |
| F20 | LOCAL-MODE 无参数源 | 同上 |
| F22 | 共线 FIXED-ENERGY 溢出 | 同上 |
| F23 | MICROCI 死循环 | 同上 |
| I-2 | E14 重建终审集成分项 | C2.md 四处引用；定义文件不在树 |
| △1–△6 / X1–X7 / R0–R7 / P-A..P-G | 迁移期任务/缺件/步骤族 | final plan / optimal_structure 可考，仅历史引用 |

## 5. grep 审计命令集（Git Bash，仓库根执行）

```bash
# 1) 生产代码注释层零决定号（允许清单外的 E 号=0）
python .zcode/fterm.py survey "E2" "src/*.f90" "src/*/*.f90" | grep ":C:"
# 2) 已清除特殊说法零回潮（活文档层与注释层；史层目录除外）
grep -rn "停放\|空转\|组合根\|立法\|双跑\|裁决\|门面\|名册\|演化环\|具名停机" \
  src/CLAUDE.md src/EXTENDING.md src/README.md src/methods/README.md \
  README.md AGENTS.md cases/*/README.md docs/plans/R8.md docs/benchmarks/
python - <<'EOF'
import sys, glob, io
sys.path.insert(0, ".zcode")
from fterm import regions
for f in sorted(set(glob.glob("src/*.f90")+glob.glob("src/*/*.f90")+glob.glob("scripts/*.f90"))):
    for ln, line in enumerate(io.open(f, encoding="utf-8", errors="replace").read().split("\n"), 1):
        if any(k == "comment" and "loud-stop" in t for k, t in regions(line)):
            print(f"{f}:{ln}")
EOF
# 3) 文档引 API 名存在性抽查（防幽灵名；输出应为空）
grep -rn "pes_reg" src/ README.md AGENTS.md
# 4) B′ 字符域检查（.md 文档应只用 U+2032；代码语境 ASCII）
grep -rn "B' " src/*.md docs/*.md README.md AGENTS.md   # 应无命中
# 5) 注释剥离等价（大批注释改动后；应报 0）
#    见 docs/plans/2026-09-17-terminology-normalization.md §7 批 B 验证脚本
```

注：审计依赖 `.zcode/fterm.py`（会话工作件，未入库）。将其正式化时移入 scripts/ 并以 Fortran 判据程序规约撰写头注——**待决定**。

## 6. 执行记录摘要（2026-09-17）

- 批 A/A2（纯文档）：12 个 .md 修缮与称谓替换；中文活文档旧称谓归零；英文活文档（AGENTS/根 README/src README/methods README/spec §4.2）同步。
- 批 B/B2（纯注释）：E27 注释×4 改写、八件判据头注统一、缩写展开×4（QCT/sec/MQC/wvn）、注释层自造词清除 352 行（loud-stop 0 残留；legislated/statute/LEGISLATED 同清）；**注释剥离等价：58 个 tracked .f90 代码层零差异**。
- 批 C（字符串）：66 处运行期消息 loud-stop→fatal；**终验：代码层全部差异恰为该替换，零意外**。
- 验证：gfortran+ifx 双编译绿；五 case 回归（换树纪律）trajectories/hist/stdout(非 timing) 逐字节一致 ×2 轮（批 B 后、批 C 后）；probe_glo 14 PASS、probe_chain_freq 两腿 0 FAIL；错误路径抽查（双容器装配）消息前缀=fatal。
- 遗留在案：① scripts 判据程序 PASS/FAIL 标签字符串中 "loud-stop" 字样（如 "(d) non-positive temperature loud-stop"）未改——超出授权批 C 范围（批 C 限 src），待指示；② check_* 头注/体内的判据出处决定号（E20/E21/E22/E26/E27 等）按"判据定义处"保留；③ leg/腿 全家待 §8-7 决定。
