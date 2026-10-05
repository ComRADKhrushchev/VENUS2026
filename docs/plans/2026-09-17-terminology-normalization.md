# 术语与代号规范化方案（v2，含 2026-09-17 用户决定）

状态：v2。2026-09-17 用户决定（称谓替换十一条 + 特殊说法清除）已并入 §2；同日用户授权全部拟修缮内容，批 A–D **已执行完毕**（见 §10 执行记录）；§8 待决定事项仍开放。
性质：文档治理专项，非阶段任务；本文件命名即示范"新顶层文档用描述性名、不再启用新字母代号"（先例：docs/superpowers/plans/2026-09-08-facade-e14-rebuild.md）。
依据：全库代号盘点（附录 A，关键坐标已逐条抽查证实）；AGENTS.md 既有规约；提交 3badaa0 注释规约先例（注释剥离等价验证）。
修订记录：v1 2026-09-17 初稿；v2 同日并入用户决定，重写 §2、调整 §5/§7/§8。

---

## 0. 问题陈述

源码与说明文档中对代码功能的说明，累计形成 20 余族文字代号（决定号、阶段字母、判据条目号、case 判据编号、轴代号、输入键遗产名……）与 60 余个自造术语。四类实际问题：

1. **字母撞车**：C/E/N/A/B/T/R/S/D/F 十个字母均有 2–6 个互不相干语义（附录 A.3）。"C1" 可以是阶段、旧库常数或净化任务；check_driver 的条目号 "(E7)" 与决定号 E7 同形。
2. **孤码不可解引用**：F14/F15/F19/F20/F22/F23（旧树缺陷登记）总表随 remediation_plan.md 删除而佚失，C5.md 照引；"I-2 carry"（C2.md 四处）与 "E 线"（src/EXTENDING.md:47）定义不可考。
3. **规约违反与陈旧文档**：E27 决定把编号写回了 4 处源码注释（违反 3badaa0"注释零 E 号"规约）；src/README.md:62 引用不存在的 `pes_reg`；src/methods/README.md 状态行过期；8 件判据程序头注保留 C1.x/C4.x 旧式引用。
4. **自造称谓不完整、不确切**：门面、组合核、名册、装配缝等草案期随意规定流传至今；停放、空转、组合根、锚、loud-stop 等特殊说法完全可用常规陈述替代。

## 1. 分层模型：史 / 活 / 码

| 层 | 范围 | 处置 |
|---|---|---|
| **史（冻结层）** | docs/plans/*.md 全部历史文本（含 K2.md E10–E14 块的 blackbox/ws 原词——E17 决定的 grep 基准）、completeness_audit.md、superpowers 归档、results 存档 | 不改写、不重编号、不翻译；孤码仅在词汇表登记，不动原文 |
| **活（活文档层）** | AGENTS.md、根/src README、src/methods/README、src/output/README、src/CLAUDE.md、src/EXTENDING.md、design_spec §4.2 现行结论区、cases/*/README、docs/benchmarks/*、R8.md（在案进行中） | 规范化主战场：按 §2 替换表与 §3 命名空间规则修订 |
| **码（源码注释层）** | src/**/*.f90 注释、scripts/check_*.f90 头注 | 注释零 E 号/零任务阶段标号；判据头注自含；英文注释中的自造词按 §2.3 清除 |

代码标识符（sec_qp、GWRITE_LEVEL、thta、input_qct.txt、elec_facade、config_roster、seam_cast、expect_loudstop 等）**不改名**——旧库兼容面与回归纪律优先；词汇表建立"行文称谓 ↔ 标识符"映射。

## 2. 称谓规范（2026-09-17 用户决定）

### 2.1 替换表（已裁，十一条）

| 旧称谓 | 规范称谓 | 英文行文 | 说明 |
|---|---|---|---|
| 容器 | **物理系统容器** | physical-system container；重复处可简"系统容器" | 不再单用"容器"一字；史层原词不改 |
| 门面 | **接口** | 以模块名直称（force_facade / elec_facade / sysdef_facade / sampler） | 四接口中文：力接口/电子接口/体系定义接口/采样接口。歧义注记见 §2.5-1 |
| 装配缝 | **初始化组装机构** | initialization assembly | 指分发行→成员 init 的通道；标识符 seam_cast 等不变 |
| 组合核 | **受力** | force evaluation（受力计算） | 指物理系统容器交付的受力计算（container_pes_i 对 g、v 全赋值）；行文多作"受力计算/受力交付" |
| 名册 | **（具体）列表** | 按语境：体系组成列表 / 片段列表 | 不再用类称"名册"；标识符 config_roster 不变 |
| 演化环 | **演化循环** | evolution loop | driver 的积分循环；需要时直接列五个调用入口，不再计数命名（"四相位体/五入口"争议随术语取消而消解） |
| 立法常数 | **物理常数** | physical constants | 守护句见 §2.5-2——"固定取值、逐字使用"的原规定不变 |
| 相 | **阶段** | stage | "C5 相"→"C5 阶段"；与程序三阶段（装配→演化→统计）靠限定词区分 |
| 双跑 | **两次运行** | two runs | "双跑逐字节一致"→"同初值两次运行输出逐字节一致" |
| 裁决 | **决定** | decision | 用户裁决→用户决定；"裁定/重裁"一并改"决定"；E 台账行文称"E 编号决定记录" |
| 成真 | **实现** | implement | 骨架→实现直接写"实现" |

### 2.2 特殊说法清除表（已裁五项 + 同精神拟清除七项）

**已裁（用户点名，用常规陈述代替）：**

| 特殊说法 | 常规表达 |
|---|---|
| 停放（停放项/停放区） | "已登记、暂不处理（的事项）"；P-1..P-7 行文称"暂缓事项 P-x"；史层表格不改 |
| 空转 | 无操作（no-op）／不做任何处理 |
| 组合根 | 不作称谓；需要时陈述："物理系统容器是体系专属代码（哈密顿量定义、体系专属方法包）的组合点" |
| 锚（锚点） | R8 文献"锚点表"→**参考值表**（读图锚→读图取值）；case 的"锚"→**判据编号**（"锚表"→"判据表"） |
| loud-stop | 中文：**报错即停（带明确错误信息）**；英文注释/消息前缀：常规措辞（建议 `fatal:`，见 §5 批 C） |

**拟清除（同精神引申，待决定 §8-7）：**

| 特殊说法 | 常规表达 |
|---|---|
| 腿（leg） | 判据条目（编号不变，(a)(b)…照用） |
| 移交 | 交接（清单） |
| 收口 | （该阶段）验收完成 |
| 后窗 | 后续扩展点 |
| 拖入（拖入件） | "新文件加入目录即被 wildcard 自动编译"；reg_*.f90 称"薄注册文件" |
| 零触碰 | 不修改既有文件（的扩展方式） |
| 点名／具名停机 | 随 loud-stop 一并清除；"点名责任键"→"错误消息指明责任键" |

保留的常规词（不清除）：骨架（skeleton，常规）、装配（assembly，常规）、判据（criterion，常规）、层 L0–L5、成员（member）、注册（registry/register）。

### 2.3 英文注释层的清除范围（拟，待决定 §8-8）

英文源码注释只清除**真自造词**，常规业界词保留：

- 清除：loud-stop（73 处注释）、composition root（1 处）、idle-chain posture、safety-net-only、leg（改 criterion item）。
- 保留：container、facade（仅当指模块时建议直称模块名）、member、registry、assembly seam→改 initialization assembly（8 处文档 + 注释零星）。

### 2.4 遗产词汇面（立档不改名）

| 词 | 含义 | 处置 |
|---|---|---|
| GWRITE_LEVEL | 记录详细级（0/1/2/5），裸缺省 2（E26） | 键名照用；词汇表注明源于旧库 GWRITE 例程 |
| sec_qp / sec_frc | 记录器分节：section q,p / section forces | recorder.f90 补一行"sec=section"展开 |
| input_qct.txt | 硬编码输入名 | 注释补 QCT = quasi-classical trajectory 展开 |
| thta | theta 的旧库拼写变体 | 输出 hist_thta.dat 沿用 |
| MQC / wvn | mixed quantum-classical / wavenumber | 定义处注释补展开 |
| e_conv / ev_kcal / e2wvn / wvn2e / hbar_code / kb_code | 物理常数标识符 | 照用；词汇表映射 |
| venus.e / results.txt / stdout_* | 可执行名 / 回归对账基准 / case 存档 | stdout_* 命名约定入词汇表 |

### 2.5 歧义与守护注记

1. **"接口"的三处撞名**：① 本表新义＝四个 *_facade/sampler 模块的类称；② Fortran `interface` 声明（container_pes_i 等抽象接口）——行文称"interface 声明/抽象接口"；③ design_spec 历史词"接口 A / 接口 B′"（指成员族）——活文档改称"**A 成员族 / B′ 成员族**"（史层 §3.3 原文不动），"A 轴/B′ 轴"随之称"A 实现轴/B′ 实现轴"。
2. **物理常数守护句**（原"立法"规定的实质内容保留，随 §8-6 决定）："本库物理常数取值按既有决定固定（如 23.0605、0.04184，立法对 0.96485132），检查程序中逐字使用，不得按 CODATA 重推。"——防止"物理常数"引发按权威值重推的误读。
3. **"受力"用法**：作名词指 container_pes_i 交付的受力计算结果；独立出现时建议写"受力计算/受力交付"。
4. **"阶段"双语境**：项目阶段（C5 阶段）与程序阶段（演化阶段）靠限定词区分，无歧义风险记录在词汇表。

## 3. 代号命名空间规则

### 3.1 现存各族登记（属主与引用限定式）

规则：**每个代号族有唯一属主文档；族内裸用仅限属主文档；跨文档引用必须带限定。**

| 族 | 属主 | 限定式示例 |
|---|---|---|
| E1–E27（决定号） | K2.md（E1–E26，E19 乱序于末）+ S.md（E27） | 裸用特权族；与 P-E 任务 E1–E5 歧义时写"P-E 任务 E1" |
| 阶段 P0/K1–K3/C1–C8/V/S/R8 | docs/plans 文件名 | "C5 阶段" |
| 总案 §5 接口规格 N1–N3 | final_completion_plan §5 | 总案 N2 |
| 审计登记注记 N1–N7 | completeness_audit_src §2 | audit N1 |
| 审计架构发现 A1–A8 | completeness_audit §2 | audit A2 |
| 暂缓事项 P-1..P-7（原"停放项"） | S.md 终局表 | S 暂缓 P-3 |
| V 阶段发现 V-1..V-6 | V.md 发现清单 | V-4 |
| 结构判据 T1–T7 | final_completion_plan §1 | 判据 T4（第二个裸用特权族） |
| 计划内决定 R/S/D 系 | 各 plan 自身 | C5-R1、C6 决定 H、C2-S3；fh2 锚点表 D 系引用写 R8:D5 |
| 判据条目号 (a)–(s) | 各 check 程序头注 | check_facade (p)；裸用仅限该程序与其 results 存档 |
| case 判据编号 A/R/T/F/B/N… | 各 cases/*/README | case_b_layer T2；README 内部可裸用 |
| 旧库常数 C1..C8 | K1 §1.1 对照表 | venus_params C1 |
| R8 子系（R8-1…、R8-A…、M-R8a…、[P1]–[P9]、A0–A2/B1–B3/O1–O10、T-三分容差） | R8.md / fh2_leps_anchors.md | 参考值 A1(R8-1)；fh2_leps_anchors.md 文件头补范围声明"本文 D/参考值/O 系仅限本文" |
| 旧例程专名（POT0/DPESHON/GWRITE…） | completeness_audit / optimal_structure 映射表 | 仅作出处标签；check_facade (s) 腿旧链编译清单照用 |

### 3.2 未来新增规则

1. **全局决定**：续 E28+（待决定 §8-1）；编号一次性、不复用、不重排。E28+ 统一记入新台账 docs/plans/E_ledger.md，阶段文件只留链接。
2. **新阶段/新顶层文档**：描述性文件名（日期-主题式），不再启用新字母。
3. **计划内决定**：一律 `<plan>-` 前缀（C5-R1 先例规则化）。
4. **判据条目号**：一程序内不复用已退役字母（check_facade 的洞保留，头注已声明"保留历史字母"）；新增子条目用数字。
5. **case 判据编号选字母**：成员首字母约定（T=thermalize、H=thermostat、G=glo_target、F=频率结构、B=boltzmann、N=normalmode）成文入词汇表。
6. **撇号字符规范**：文档语境统一 U+2032（B′），代码/标识符语境 ASCII（B'）（待决定 §8-2）。
7. **中英双写**：维持 AGENTS.md"match the surrounding language"；词汇表立双向对照，同一概念每语言一个规范形。

## 4. 孤码登记表（入词汇表，不动史层原文）

| 孤码 | 已知含义（抢救自上下文） | 状态 |
|---|---|---|
| F14 | RANDST 双调用缺陷 | 定义总表随 remediation_plan.md（2026-09-06 删）佚失；含义在案（C5 附录/final plan §6） |
| F15 | 双随机流缺陷 | 同上 |
| F19 | 振动角动量剥离 | 同上 |
| F20 | LOCAL-MODE 无参数源 | 同上 |
| F22 | 共线 FIXED-ENERGY 溢出 | 同上 |
| F23 | MICROCI 死循环 | 同上 |
| I-2 | E14 重建终审集成分项（C2.md 四处引用） | 定义已佚 |
| △1–△6 / X1–X7 / R0–R7 / P-A..P-G（含内层任务号） | 迁移期任务/缺件/步骤族 | 含义在 final plan/optimal_structure 可考，仅历史引用 |
| "E 线" | = E 编号决定记录 | EXTENDING.md:47 孤词，修缮项 4 改写 |

## 5. 修缮清单

**批 A——活文档缺陷修缮（纯 .md）：**

1. 新建 docs/glossary.md（§2 全表 + §3 命名空间登记 + §4 孤码表 + §2.4 遗产词汇面 + grep 审计命令集）。
2. src/README.md:62 —— `pes_reg` 幽灵名 → `container_bind_pes`。
3. src/methods/README.md:20-21,86 —— 状态行重写（dist_* 11 员已齐、只差 qmmicro；"phase C"措辞清理）。
4. src/EXTENDING.md:47 —— "E 线裁决记录" → "E 编号决定记录（docs/plans/K2.md、S.md）"。
5. B′ 写法统一：src/methods/README.md:21,67 的 ASCII B' → B′（U+2032）。
6. AGENTS.md / 根 README 退役词核查：ws / h_form / h_deliver 每处出现均带 E27 指针。
7. cases/*/README 判据表头与 fh2_leps_anchors.md 文件头补范围声明。

**批 A2——活文档称谓替换（§2.1/§2.2 已裁项）：**

8. 中文活文档（src/CLAUDE.md、src/EXTENDING.md、cases/*/README.md、docs/plans/R8.md、docs/benchmarks/*.md、根 README）按替换表改写，约 200 处（规模实测：锚 104、容器 39、双跑 17、立法 14、裁决 13、停放 7、门面 3、装配缝 3、组合根 2、空转 1、成真 1）。
9. 英文活文档（AGENTS.md、根/src README、methods/README、spec §4.2）：legislated constants→physical constants（附守护句）、composition root→组合点陈述、E27 段措辞同步。

**批 B——源码注释规约（零代码字符变动）：**

10. force_facade.f90:7,156 与 elec_facade.f90:10,196 —— 4 处 E27 注释改写为无编号现行结论（字符串字面量例外清单保留：input.f90:182,282,317,368 之 E21/E22、elec_facade.f90:212 之 E13）。
11. 判据头注风格统一：8 件旧式（check_rng / check_linalg / check_geometry / check_specio / check_sysdef / check_roster / check_input / check_recorder 的 C1.x/C4.x 引用）改为"criteria defined FIRST in this header"自含式。
12. 缩写展开补注（各一行）：input_qct.txt=QCT、sec=section、MQC、wvn。

**批 B2——源码注释英文自造词清除（纯注释，待 §8-8 决定范围）：**

13. loud-stop 注释 73 处 → 常规措辞（"abort with a named error message"）；composition root、idle-chain posture、safety-net-only、leg 同步清除。

**批 C——运行期字符串替换（触碰代码，单列，待单独授权）：**

14. src 内 65 处 STOP/write 消息前缀 `loud-stop (...)` → 建议统一 `fatal: (...)`；**expect_loudstop 断言的子串（如 'thermostat_init'、't_bath'）不动**，故判据程序预期不需改；验证腿见 §7。

**批 D——治理落章：**

15. AGENTS.md Conventions 增"术语与代号"节（替换表要点 + 命名空间规则 + grep 审计命令）；阶段收尾模板补"术语同步"步。

**相邻但不并入（另行）：** src/build/gfortran/sys/well2 陈旧演练对象属 build 卫生，`make clean` 即除。

## 6. 治理条款

1. **词汇表（docs/glossary.md）为术语唯一权威**；新代号入册两要件：定义处在册 + 跨文档引用带限定。
2. 阶段收尾固定新增"术语同步"步。
3. **grep 审计腿**（命令集入词汇表附录，择要入 AGENTS.md）：
   - 源码注释 E 号残留 = 0（允许清单：5 处字符串字面量）；
   - 文档引 API 名的存在性抽查（防 pes_reg 类幽灵）；
   - B′ 字符域检查（.md 用 U+2032；代码用 ASCII）；
   - 已清除特殊说法零回潮（停放/空转/组合根/锚/loud-stop 在活文档层与注释层 grep = 0，史层目录除外）。
4. 大批注释改动沿用注释剥离等价验证（3badaa0 先例脚本）。

## 7. 执行批次与验证

| 批 | 内容 | 性质 | 验证腿 |
|---|---|---|---|
| A | 修缮项 1–7 | 纯 .md | git diff 仅涉 .md；grep 审计零命中 |
| A2 | 修缮项 8–9 | 纯 .md | 同上；替换前后逐词计数对账（旧词活文档层 = 0） |
| B | 修缮项 10–12 | 纯注释 | 注释剥离后与 HEAD 逐字符一致；`make` + `make FC=ifx` 双编译绿；五 case 回归 results.txt 逐字节一致；check_* 全绿 |
| B2 | 修缮项 13 | 纯注释 | 同批 B |
| C | 修缮项 14 | **代码触碰** | 全部 check_* 重跑全绿（重点 child-process 报错即停条目族：check_dist_thermostat 等）；五 case 回归 results.txt 逐字节一致（正常路径不含错误消息）；错误路径 stderr 抽样人工比对；双编译绿 |
| D | 修缮项 15 | AGENTS.md/模板 | 复跑 grep 审计 |

时序约束：R8 执行期间有不修改共享代码的纪律——批 B/B2/C 涉 src/ 与 scripts/，是否等 R8 收口或作显式例外，见 §8-4。

## 8. 待决定事项

**已由 2026-09-17 决定解决**：演化循环命名（原 §8.3）；动词合并议题整体被"决定/实现/两次运行"替换表取代（原 §8.4）；"停放/停驶"术语取消、纪律改陈述式（原 §8.5）。

**仍开放：**

1. E 序列未来：续 E28+ 并新立 docs/plans/E_ledger.md 统一记账（建议）？还是封线改"日期+属主文档"引用？
2. B′ 字符规范：文档 U+2032 / 代码 ASCII（建议，成文既成事实）？还是全面 ASCII？
3. 词汇表居所：docs/glossary.md（建议）？
4. 时序：批 B/B2/C 等 R8 收口后执行，还是作为显式例外现在并入？
5. 批 C 授权：65 处运行期字符串前缀替换（建议 `fatal:`）是否照准？——这是唯一触碰代码的批次。
6. 物理常数守护句（§2.5-2）是否照准？
7. 拟清除表（§2.2 下半：腿/移交/收口/后窗/拖入/零触碰/点名）是否照准？
8. 英文注释层清除范围（§2.3：仅真自造词；container/facade 等常规词保留）是否照准？

## 9. 明确不做清单

- 不重编号任何历史台账（E / 阶段 / 判据条目 / case 编号 / V- / P- / 判据 T 系）——grep 基准与法定可追溯性优先。
- 不改代码标识符，含旧库兼容词汇面（sec_qp、GWRITE_LEVEL、thta、input_qct.txt、results.txt、venus.e）与机制名（elec_facade、config_roster、seam_cast、expect_loudstop、scheme_code_name）。
- 不改写 docs/plans 历史块（K2 的 blackbox/ws 原词是 E17 决定的 grep 基准；F 系照引不动；"停放项"表格原样）。
- 不翻译、不统一史层文档；替换只约束活文档层与源码注释层。
- 不借本专项顺手做功能或结构调整。

---

## 10. 执行记录（2026-09-17，用户授权后落地）

授权原文：「授权同意拟修缮的内容，待决定事项之后再说」——§5 十五项全数执行；受 §8 缓决约束的部分按保守口径（拟清除表下半与 leg/腿 未动）。

**批 A/A2（纯文档，12 个 .md）**：修缮项 1–9 全部落地。docs/glossary.md 建立（称谓表+命名空间登记+孤码表+审计命令）；pes_reg 幽灵名、methods/README 过期状态行、"E 线"孤词、B′ 分裂、范围声明（五 case README + fh2 参考值表头）、退役词指针核查完毕；中文活文档称谓替换两轮（脚本+全角标点补刀），旧称谓（停放/空转/组合根/锚/立法/双跑/裁决/门面/组合核/装配缝/名册/演化环/具名停机）在活文档层 grep 归零；英文活文档（AGENTS/根 README/src README/methods README/spec §4.2）同步（legislated→physical/fixed-by-decision、composition root→composition point、composed force kernel→force evaluation、"four-phase step"→演化循环五入口、ws%repr 陈述句清除）。执行中发现并顺带处理：EXTENDING"同_folder"错字、"全四相循环"→"全流程循环"（cases）、radau/chain_freq 的"propagator 立法/立法不变"→"明文规定/规定不变"。

**批 B/B2（纯注释，零代码字符变动）**：E27 注释×4（force_facade×2、elec_facade×2）改写为无决定号现行结论；八件判据头注统一为"criteria defined FIRST in this header"自含式（rng/linalg/geometry/specio/sysdef/roster/input/recorder），体内 C4.x/C1.x 引用同步清除；缩写展开×4（driver 的 QCT、recorder 的 sec、elec_facade 的 MQC、consts 的 wvn）；consts 头注 LEGISLATED 措辞、legislated/statute 家族（34 行）一并清除；注释层自造词替换共 352 行，loud-stop 注释层归零（src+scripts）。**验证：注释剥离等价（区段感知，剔除纯注释行）——58 个 tracked .f90 代码层零差异。**执行中两处脚本损伤（consts 丢注释标记、recorder 丢 sec_qp 声明行）被等价验证当场抓获并修复——验证腿有效性的一次实证。

**批 C（字符串，66 处）**：src 运行期消息前缀 loud-stop→fatal（66 行，较预估 65 多 1）；expect_loudstop 断言依赖的责任入口名/责任键名子串未动。**终验：对 HEAD 施以同变换后与工作版逐字符一致（零意外差异）；错误路径抽查（双容器装配）消息前缀=fatal。**

**批 D**：AGENTS.md Conventions 增"Terminology"节（规范称谓+退役词+引用限定+收尾术语同步步）；E 号注释规约补登记两类在案例外（check 头注判据出处注、5 处字符串字面量）。

**验证腿（全数通过）**：
1. gfortran + ifx 双编译器全量重建绿（批 B 后、批 C 后、还树后共三轮）。
2. 五 case 回归（按 R8 换树纪律：撤 fh2_leps 树重建后跑）两轮（批 B 后、批 C 后）：六运行目录 trajectories.txt / hist_*.dat / stdout（滤 timing 行）与基线**逐字节一致**，stderr 全空；probe_glo 14 PASS；probe_chain_freq 两腿 exit 0、0 FAIL。
3. 全判据套件重编重跑：**27 件全绿**（0 failure(s)，无任何 FAIL 行；含 check_facade 旧链位级 (s) 腿）+ 还树后 check_fh2_leps **18 PASS/0 FAIL**。
4. grep 审计五条全过（生产注释层 E 号=0、活文档旧称谓=0、幽灵名=0、B′ 域=0、规则自指命中除外）。

**遗留在案（未做，各有原因）**：
- scripts 判据程序 PASS/FAIL 标签字符串中的 "loud-stop" 字样（如 '(d) non-positive temperature loud-stop'）——批 C 授权范围限 src/，scripts 标签属判据程序自身输出文案，待指示是否同改。
- leg/腿 全家（§8-7 缓决；probe_chain_freq 有同名标识符 `leg`，盲改有风险）。
- check_* 头注/体内的判据出处决定号（E20/E21/E22/E26/E27、E7 条目号）——按"判据先行的定义处"定位保留，已登记为规约例外。
- 拟清除表下半（腿/移交/收口/后窗/拖入/零触碰/点名）与 §8 其余事项一并待用户决定。

**过程注记**：判据套件首轮 check_facade 构建失败系本批执行脚本的命令拼接 bug（三元优先级致命令重复、双 -J），修正后绿——非主体代码问题；.zcode/fterm.py（区段感知工具）与判据日志存 .zcode/（会话工作区，未入库），正式化事项已登记 glossary §5 注。

---

## 附录 A：代号族全景盘点（2026-09-17，v1 存档）

### A.1 族清单（族 → 例 → 含义 → 定义处 → 现状）

| 族 | 例 | 含义 | 定义处 | 现状 |
|---|---|---|---|---|
| E1–E27 | E27 | 全局决定台账（E9 起升格用户决定） | K2.md（E1–E26，E19 乱序于末）、S.md（E27） | 活跃；E27 后停更 |
| 阶段 | P0/K1–K3/C1–C8/V/S/R8 | 执行序大阶段（C3 缺号=E18 撤销留空位） | final_completion_plan §2 + docs/plans 文件名 | P0..S 完成；R8 进行中 |
| 阶段内步骤 | C2.1、V.5、S4.1、R8-A..D | 步骤复选框编号 | 各 plan | 活跃（R8） |
| 判据 T1–T7 | T4 | 七条终态结构判据 | final plan §1；核验在 S.md | 活跃 |
| 总案 N1–N3 | N2 | 骨架接口规格（N2 已被 E18 撤编） | final plan §5 | 历史，仍被 P 节当出处引 |
| audit N1–N7 | N1 | 完成度复审登记注记 | completeness_audit_src §2 | 活跃 |
| 暂缓事项 P-1..P-7 | P-3 | 停放项终局表（原名"停放项"） | S.md | 活跃 |
| P-A..P-G（内层 A1–A5…G1–G4） | P-C | 旧净化主线（移植剔除点清单） | final plan §6 | 历史 |
| X1–X7 | X3 | 收敛轨放置核对 | final plan §7 | 历史 |
| R0–R7 | R3 | 迁移策略步骤 | optimal_structure §5 | 历史 |
| R8 系 | R8-1、R8-L、M-R8a | 文献基准族与里程碑 | final plan §9、R8.md | 活跃 |
| 计划内决定 | R-a..R-d、R1–R5、决定 H–O、D1–Dn、S1–Sn | 各 plan 范围/设计决定 | C2/C5/C6/C7 等 | 活跃引用（R5/R1 跨文件被引） |
| E14 重建 Task 1–14 | T11 | superpowers 重建任务号 | superpowers plans/sdd | 历史；T11 被 C4.md 引 |
| 判据条目号 (a)–(s) | (p) | 判据程序子判据（原"腿"） | 各 check 头注 | 活跃；check_facade 七条目退役留洞；数件跳号（normalmode 缺 n、barrier 缺 r、surface 缺 l、input 缺 o5） |
| case 判据编号 | b_layer T2、chain_freq F1 | case README 判据编号（原"锚"） | 各 cases/*/README | 活跃 |
| 轴 | A 轴 / B′ 轴 / 电子轴 | 双实现矩阵两维 | final plan T3、spec §3.3 | 活跃 |
| class-A/B | class-B | A 族接口内两类成员 | K2 E9–E14 | 半退役（注释仍用） |
| Interface B（新义） | Interface B | spec §4.2 把 elec_facade 称电子轴接口 | spec §4.2 | 与 B′ 并存，B 字母三义 |
| E10 公理 A/B/C | 公理 B | ws 三公理 | K2 E10 | 历史 |
| branch 4′ | 4′a | 接口伪代码分支 | K2 E10/E12 | 历史 |
| 处置符号 | ● ◐ ✕ △、△1–△6 | 结构映射判定 | optimal_structure §2/§3 | 历史 |
| V-1..V-6 | V-4 | V 阶段发现清单 | V.md | 活跃（V-4 立法对被 CLAUDE.md 引） |
| 常数 C1..C8 | C1=0.04184 | 旧库 venus_params 换算常数 | final plan 全局约束 1、K1 §1.1 | 活跃（立法对） |
| 输入键 | GWRITE_LEVEL、ELEC_METHOD… | 键盘语法 | C4 S3 等 | 活跃；退役键作文史存 |
| 记录分节 | sec_qp / sec_frc | recorder 门 | recorder.f90 | 活跃 |
| ws 家族 | ws、h_form、r_diab | 已退役互换词汇 | K2/C2/S | 史层残留 + 活文档退役叙述 |
| 旧例程专名 | POT0、GWRITE、THRMAN | 出处标签族 | audit/optimal_structure | 历史；check_facade (s) 旧链编译照用 |
| L0–L5 | L2 | 层号 | spec/src README | 活跃（目录名不带层号） |
| fh2 子宇宙 | [P1]–[P9]、A0–A2、D1–D10、T-三分 | 文献/参考值组/决定点/容差类 | fh2_leps_anchors.md | 活跃（R8-1） |
| ①②③④ | ① | 逐点列举 | K1/K2/V/S/R8 | 活跃 |

### A.2 术语量级（v1 存档；替换规模实测见 §5 批 A2）

判据/判据程序 153 处、容器 81、装配 55、立法系 56、裁决 51、loud-stop 139、成真 53、停放 49、腿 32、移交 35、骨架 37、缝 20–22——中英双写按文档语言切换是制度化现状（AGENTS.md 规约），此前无对照表。

### A.3 撞车矩阵（字母 → 语义数 → 各语义）

| 字母 | 语义数 | 语义清单 |
|---|---|---|
| C | 4 | 阶段 C1–C8；旧库常数 C1..C8；P-C 任务 C1–C3；C5-R1 类计划决定 |
| E | 3 | 决定 E1–E27；P-E 任务 E1–E5；check_driver 条目号 (E7) |
| N | 3 | 总案 N1–N3；audit N1–N7；chain_freq 编号 N1–N3 |
| A | 6+ | A 族接口；A 轴；audit A1–A8；P-A 任务 A1–A5；case 编号 A1–A6；fh2 参考值组 A0–A2；E10 公理 A |
| B | 5 | B′ 族；Interface B（电子轴）；class-B；B 层；roster 片段 B |
| T | 4 | 判据 T1–T7；b_layer 编号 T1–T3；E14 重建 Task 1–14；R8 容差 T-三分 |
| R | 6 | R0–R7；R8；R-a..R-d；C5 R1–R5；radau 编号 R1–R4；radau B 系数 R 系判据 |
| S | 3 | 各 plan 范围裁定 S1–Sn；S 阶段及其任务 S1–S4；S2=新审计 |
| D | 2 | 各 plan 设计决定 D1–Dn；fh2 决定点 D1–D10 |
| F | 2 | 旧库缺陷 F14…（孤码）；chain_freq 编号 F1–F3 |
| V | 2 | V 阶段；V-1..V-6 发现清单 |

### A.4 已证实的具体缺陷坐标

- src/README.md:62 `pes_reg` 幽灵（实为 container_bind_pes）。
- src/facade/force_facade.f90:7,156；src/facade/elec_facade.f90:10,196 —— 注释中 E27（违 3badaa0 规约；字符串字面量例外：input.f90:182,282,317,368 与 elec_facade.f90:212）。
- src/methods/README.md:20-21,86 状态行过期（"remaining 9 arrive C5"，实际 11 员齐）。
- src/EXTENDING.md:47 "E 线" 孤词。
- B′ 分裂：src/EXTENDING.md:31（U+2032）vs src/methods/README.md:21,67（ASCII）。
- 判据头注旧式 8 件：check_rng / check_linalg / check_geometry / check_specio / check_sysdef / check_roster / check_input / check_recorder。
- "四相位体"（K2/C7）vs "五入口"（audit_src/CLAUDE.md）计数名并存——随"演化循环"命名消解。
- K2.md 台账 E19 乱序排于 E26 后；E27 定义散在 S.md。
- loud-stop 分布：src 138 处 = 65 处运行期字符串 + 73 处注释；判据程序以 expect_loudstop 子进程断言 stderr（断言子串为责任键名，非前缀）。
