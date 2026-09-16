# S - 收尾相(新程序完成度复审 + 文档收尾)

状态: **已完成**(2026-09-16;执行中落裁 E27 并当日落地,见下节)
日期: 2026-09-16
依据: docs/final_completion_plan.md §8(Closing S)+ §1 T1-T7 + §9(S 行判据);
V.md →S 移交区(2026-09-16,能力面→证据对照表 + 停放项终局四条 + 既有登记)。

## E27 裁决(2026-09-16 用户裁决,S 相执行中落裁——容器组合根)

**动因**:用户以自定义 Hubbard 哈密顿量传播方法压测架构,暴露双轴设计(方法居 methods/、
体系居 systems/)迫使"体系→方法"通用语(ws 数据词汇表)存在,其代价链为:ws schema →
force_eval 的 repr switch → 每新增交付形态改共享代码 → 容器格式收紧 + 新体系两处注册。
E10(ws 容器)/E13(h_form·h_deliver 撮合)/E14(方法轴)的相关部分被本裁覆盖。

**裁决内容**:
1. **容器 = 组合根**:H 定义与体系特有的非绝热演化方法同居 `systems/<sys>/` 一个文件夹
   (共享私有模块,无通用语);方法包的 reg 文件可居容器文件夹(prebuild 通配符本就聚合
   `systems/*/reg_*.f90`,一处注册)。方法换用仍走 ELEC_METHOD 名字键(单表纯名字查找,
   E13 纯形式权威不变; adiabatic 退化行仍为缺席缺省,键可命名库内行或容器行——保住
   "同一体系跑纯绝热参考"的实验面)。
2. **force_eval 契约收窄为普适组合核交付**:容器交付已组合的活性力核 g[eV/A] 与势能 v[eV]
   (container_pes_i(q, mass, g, v),全赋值契约);换算唯一点(ev_kcal·e_conv, T4)与
   E=T+V 折叠留在门面,**select case 整体退役**。原 r_diab 分支的对角化+占据折叠物理不删
   ——在 check_facade (p) 腿的模块面(leps_diab + eig_sym 折叠等价)存活;H 形容器的通用
   折叠助手为候选后续件(git 历史可取)。
3. **elec_ws 退役**(文件删除):repr 标签/构造通道/一步生命周期/传送包装问题全部消解;
   "有状态成员 wrapper 终形"停放项(P-6 之一)随之**关闭**。elec_facade 保留四入口签名契约
   (prop(sta,dt)/match(sta)/sync(s))、occ/n_hop 运行期持有者(新增成员写回入口
   occ_set/n_hop_set)、退化缺省行;method_reg 签名去 h_form,h_form/h_deliver 撮合退役
   (container_init 只剩名字绑定 + 闲链姿态)。
4. **methods/ 章程收窄**:通用方法包(adiabatic 退化行等)与组件库(amp_*/fssh_g_ij/
   hop_decide/momentum_rescale 等 level-2 目录)为主;体系专属方法包按 1 居容器文件夹。

**与提案措辞的偏差(明记)**:ELEC_METHOD 键未按"容器自带即缺省"实现,改为统一的单表
名字查找——"容器缺省"方案使体系无法显式选择绝热参考跑(键缺席与显式 adiabatic 不可区分,
需引入键在场跟踪);统一名字键保住方法-换键实验面且零优先级分支,一处注册不受影响。

**落地与验证(2026-09-16)**:生产 5 文件(force_facade 重写/elec_facade 重写/elec_ws 删除/
reference_box·reg_reference·reg_adiabatic 适配);判据程序 6 个适配(check_facade 手术:
(a)(b)(c)(g)(i)(k)(r) 七腿随主体退役并在头注明注,(d)(e)(f)(h)(j)(l)(m)(n)(o)(p)(q)(s)
 存活——(p) 腿承接折叠物理;五小件签名/配方机械适配)。验证四腿:①gfortran+ifx 双编译
器全量重建绿;②五 case 终跑与 V 相入库 results.txt 逐字节一致(算术路径未变:标量通道
g=v 直达,门面折叠表达式同一);③六判据程序重编重跑全绿——check_facade 0 失败,(s) 旧链
位级腿三势仍 max|dF|=max|dE|=0,(p) 折叠等价 1e-12;④T4(23.0605 仅 force_facade)与
T5(全库零 RANDOM_NUMBER)复验不动。S 相后续件按 E27 后树执行:T2 演练用新签名重做、
S2 审计审 E27 后树、S3 文档批书写新契约、停放项表 P-6 的 wrapper 条目关闭。


## 目标

按总案 §8 完成收尾四件(S1 案例体系收口 / S2 完成度复审 / S3 文档收尾 / S4 外部对账),
逐项核验 T1-T7 结构判据与不变式三问抽查,V 移交停放项按各处登记终局清点。
本相原则上**零源码改动**(纯复审/文档/验证);任一 T 判据或复审发现红项,
按项目惯例单独立案裁决后再动代码,不在本相内顺手修。

## 现状事实(开工前已核实)

- 五 case 的 results.txt 均已由 V 相生成入库;V.5 终跑抽验 harmonic/surface 两 case 逐字节复现。
- src/README.md 在案且为现行层表(已含 E12/E14/E18/E19 调整);src/CLAUDE.md、EXTENDING.md 不存在。
- spec §4.2(正式签名 + 退化点良定 + 双实现记录)仍为空白待回填;optimal_structure.md 状态行仍为
  "design draft";根 README.md 的 layout 仍指向 venus2026/ 为程序源(过时,需更新为 src/ 新树 + 旧树参考)。
- 工作区仅一个未追踪文件 AGENTS.md(本会话 /init 产物,随 S3 文档批一并提交);
  src/ 根散落 ad-hoc 编译残留 .obj/.mod(均未被 git 追踪,.gitignore 已覆盖)——收尾卫生项,直接删除。
- 总案 T3 字面材料名(RST/EMT-NN)属完整形态资产、不在本仓库(根 README 在案)——T3 按其判据本意
  (A 轴 ≥2 实现 / B′ 轴 ≥3 / 配置面双范式)以本仓库现实载体核验,字面名标注为历史目标。

## 任务分解

### S1 案例体系收口(总案 §8 S1)

- [x] S1.1 五 case 全量终跑核验(2026-09-16,E27 后二次执行):三直跑 case 证据行与入库
  results.txt 逐字节一致(harmonic/radau 原生一致;surface 存档为 LF 归一化,剥 CR 后一致);
  b_layer 三子跑一致 + probe_glo 14 PASS;chain_freq 两腿探针空白归一口径一致(存档手工折行),
  退出码全 0、stderr 空(chain_freq 的 STOP 0 为设计路径)。**覆盖方式:全量,未用降级路径。**
- [x] S1.2 判据与容差记录完备性:五 case README §3 预期锚/§4 对照齐全(判据先行在案);
  能力面→证据对照表八行全有载体。

### S2 新程序完成度复审(总案 §8 S2;= T7)

- [x] S2.1/S2.2 全树 50 文件 234 例程三问复审完成 → **docs/completeness_audit_src.{md,json}**
  (零红项;三问 (c) 机械 grep 全树零命中;登记注记七条各有在案登记)。
- [x] S2.3 不变式深审抽查四处(driver/force_facade/recorder/sampler)——审计报告 §3。

### S3 文档收尾(总案 §8 S3)

- [x] S3.1 src/CLAUDE.md(新树构建/运行/判据/契约要点,承接旧树 CLAUDE.md 角色)。
- [x] S3.2 src/EXTENDING.md(五类扩展点按 E27 组合根契约 + 落地演练清单)。
- [x] S3.3 src/README.md 收尾态(Build 区重写为 K3 现实;utils 补 hessian;dynamics 行
  predicate 残留修正;挂 CLAUDE/EXTENDING 衔接)。
- [x] S3.4 design_spec_completeness.md §4.2 回填(接口 A/B/B′ 落地签名 + 退化点良定
  三条 + 双实现记录按本仓库现实载体,RST/EMT-NN 字面名历史目标注记)。
- [x] S3.5 optimal_structure.md 状态行 → implemented(含裁决链漂移注记);根 README.md
  全文重写(src/ 新程序 + systems/ 组合根 + venus2026/ 只读参考 + 状态段)。
- [x] S3.6 提交文档批(含 AGENTS.md 的 E27 段落微调)。

### S4 外部对账(总案 §8 S4)

- [ ] S4.1 Notion "VENUS-TDHF 泛化" 任务勾选态同步——**用户侧条目**,S 收口提醒。

### T1-T7 结构判据逐项核验(证据入本件验收区)

- [x] T1 文件结构:50 源文件清点与 src/README 层表相符;漂移注记三条(predicate 撤 E18/
  analysis 并入 output/utils 的 hessian 待 README 补行——S3.3 一并修)。
- [x] T2 零触碰扩展演练(E27 新签名二次执行):systems/well2/(pot+reg 两新文件)拖入 →
  第二次 make 全绿(登记过的 depend.mk 单轮迟滞)→ registry_gen 聚合 reg_well2(use+call)、
  depend.mk 三边 → **git diff tracked 零改动**(仅新增未追踪文件夹)→ 撤销后重建绿、树干净。
- [x] T3 双实现矩阵(按判据本意以现实载体):A 轴 = reference 三势(harmonic/morse/leps,
  check_facade (q) 解析 vs 中心差分 + (s) 旧链位级,E27 后仍全 0)≥2 ✓;B′ 轴 = 11 个 dist
  成员(五 case + 11 个 check_dist_* 判据程序覆盖)≥3 ✓;配置面双范式 = 气相四 case
  (paradigm=0)+ 表面 case(POSCAR 胞,paradigm=1)✓;电子轴 = adiabatic 退化成员(完整形态
  非绝热轴不在本仓库——总案 T3 字面名 RST/EMT-NN 同此,历史目标注记)。
- [x] T4 单位契约:23.0605 与 ev_kcal 仅 force_facade(定义+注释);组合换算式
  `-g*ev_kcal*e_conv` 唯一(facade fold 点);e_conv 定义在 consts,其余消费均为 I/O 侧
  kcal/mol↔内部暂存(成员参数缝/报告列)——力边界换算唯一,T4 判据满足。
- [x] T5 单随机流:全库 grep RANDOM_NUMBER 零命中( rng 全流 RAND0 族)。
- [x] T6 三相:driver 审读——装配 10 步→演化(轨迹环×步进环:prop_step/elec_prop/
  mqc_statistic/rec_frame/container_term)→统计(fin_eval 单趟+分类计数+直方图+漂移);
  演化环内仅门面/传播器/记录器入口;E18 漂移注记(谓词位=container_term 容器槽)。
- [x] T7 = S2(新审计产出在案)。

### 停放项终局清点(逐条裁决;E27 后复审)

| # | 登记处 | 内容 | S 相裁决 | 去向 |
|---|---|---|---|---|
| P-1 | V-5/chain_freq §2 | 装配期频率导出通道未接线(机制存在零生产接线,四证据) | **维持停放 → 候选后续件** | E27 后形态明确:H 形容器可自带 Hessian→谱导出(容器文件夹私有通道);触发条件=首个多面容器落地 |
| P-2 | V 停放②/results_check_facade_ifx.txt | ifx 下 check_facade (g)/(r) 两红腿 | **关闭(E27 消解)** | 判据腿随 ws 构造通道退役而消失;ifx 存档留历史参考 |
| P-3 | C8→V 第 5 条 | 采样器→统计相初值快照通道(VI/OAMI/AMAI/AMBI/SDA/SDB 生产者,解锁 fin%ang(4..18)) | **维持停放** | 随 C8 触发条件;E27 不影响(recorder 读 holder 的 getter 面未动) |
| P-4 | V 停放④ | check_facade (g) FAIL 行不打印 got-vs-expected(scripts/ 改进注记) | **关闭(E27 消解)** | (g) 腿对象已退役;其余容差判据均打印的通例仍在 |
| P-5 | C8 终审停放区 | hist_thta 气相门控候选 / col_occ 守卫 / fin_t 初值参考量(=P-3 族) / SDA/SDB(=P-3 族) / plot_hist allocated 守卫 / lattice%q_ref 占位 | **六条全维持停放** | 均为设计变更/守卫候选,走流程触发;无一条因 E27 失效 |
| P-6 | C2 停放区 | (r) 1e-15 性状 / fssh·iesh 后窗件 / r_adiab 构造通道 / 有状态成员 wrapper 终形 | (r) **关闭(对象腿随 E27 退役)**;wrapper 终形 **关闭(E27 消解:ws 传送退役,成员数据走容器文件夹私有模块,写回入口 occ_set/n_hop_set 已落地)**;fssh·iesh 后窗件 **维持→形态改写**(E27 后此类方法包居容器文件夹,待第一消费者);r_adiab 构造通道 **维持**(旧树材料) | fssh/iesh 触发条件不变(非绝热轴开工) |
| P-7 | b_layer README | NU_COLL 携带无消费面注记 | **维持(如实注记)** | 无消费面即无行使面;首个消费者出现时处置 |

### 卫生项

- [ ] 删除 src/ 根 ad-hoc 编译残留(.obj/.mod,均未追踪);确认 make clean 后全量重建双编译器绿。

## 验收标准

1. `make`(gfortran)与 `make FC=ifx` 双编译链接绿(卫生项后全量重建)。
2. S1 终跑核验按声明的覆盖方式全绿;能力面→证据对照表无空行。
3. docs/completeness_audit_src.{md,json} 产出:src/ 全例程三问结论齐,无"体系知识漏入"红项
  (如有红项如实记录并单独立案,不阻塞文档收尾)。
4. T1-T7 逐项 PASS(或如实标注漂移依据),证据行入本件验收区。
5. S3 文档六件齐;git 提交组清晰(建议两批:复审+判据批、文档批)。
6. 停放项终局表覆盖 P-1..P-7 无遗漏。

## 风险

- 统计 case 终跑时长(S1.1 已备降级路径,覆盖方式如实标注)。
- S2 复审若发现红项(体系知识漏入等)——如实入审计报告并单独立案,本相不改源码。
- 总案字面判据与现实漂移(T1/T3)——以裁决链(E12/E14/E18/E19)为依据注记,不擅自改总案原文。


## 验收证据(S 收口,2026-09-16)

1. **双编译器**:make clean 后 gfortran + ifx 全量重建绿(卫生项 + E27 后复建共两轮)。
2. **五 case 回归**:E27 前后各一轮终跑,均与 V 相入库 results.txt 逐字节/空白归一
   一致(E27 轮见 S1.1 证据行)——算术路径未变的直接证据。
3. **判据程序**:六个受 E27 影响件重编重跑全绿(check_facade 0 失败、(s) 三势位级
   max|dF|=max|dE|=0、(p) 折叠等价 1e-12;driver/final_state ALL PASS;
   thermalize/recorder/propagator 0 失败);未受影响件不涉改动面。
4. **T1-T7**:逐项 PASS(T3 电子轴为绝热退化点,完整形态非绝热轴不在本仓库——判据
   本意与现实载体注记在案);T2 演练两度执行(ws 时代 + E27 新签名),零触碰结论同。
5. **S2/T7**:docs/completeness_audit_src.{md,json}——50 文件 234 例程,零红项,
   零未分类;旧审计保留历史参考。
6. **停放项终局表**:P-1..P-7 逐条裁决(关闭四:wrapper 终形、(r) 性状、(g)/(r) 红腿
   两条目——皆随 E27 消解;维持七;P-1 升格候选后续件且 E27 后形态明确)。
7. **S4(用户侧)**:Notion "VENUS-TDHF 泛化"(Phase P1)勾选态待用户同步——本件提醒。

## 交付物清单

- docs/completeness_audit_src.{md,json}(S2/T7)
- src/CLAUDE.md、src/EXTENDING.md、src/README.md(收尾态)
- docs/design_spec_completeness.md §4.2 回填、docs/optimal_structure.md 状态、README.md
- docs/plans/S.md(本件:E27 裁决记录 + T/S 证据 + 停放项终局表)
- AGENTS.md(工作区指导,会话 /init 产物 + E27 微调)
