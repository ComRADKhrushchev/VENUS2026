# S - 收尾相(新程序完成度复审 + 文档收尾)

状态: **执行中**(用户确认 2026-09-16 开工;同日 S 相执行中落裁 E27,见下节)
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

- [ ] S1.1 五 case 全量终跑核验(gfortran 树):stdout 证据行与入库 results.txt 逐字节一致
  (timing 行除外),退出码 0、stderr 空;b_layer/chain_freq 两统计 case 若时长不可承受,
  依 V.5 抽验模式各补一 case,其余以各自 V 相双跑判据在案作记录性收口(执行时如实标注覆盖方式)。
- [ ] S1.2 判据与容差记录完备性清点:每 case README/results.txt 载明判据、锚与容差;
  能力面→证据对照表(V 移交)逐行有载体、无空行。

### S2 新程序完成度复审(总案 §8 S7;= T7)

- [ ] S2.1 例程清单:src/ 全树(含 systems/)模块与例程逐个登记(脚本辅助 + 人工过)。
- [ ] S2.2 三问复审:每例程过 spec 三问(只处理 Q/P/s/F 吗 / 按契约请求服务吗 /
  内含体系分支或常数吗),结论入新审计 docs/completeness_audit_src.{md,json}
  (旧 docs/completeness_audit.{md,json} 保留为历史参考,新件头注声明取代关系)。
- [ ] S2.3 不变式三问抽查:driver / force_facade / recorder / sampler 各抽一例程深审
  (原案的 predicate 抽查位按 E18 撤销,以 force_facade 容器终止槽位替补)。

### S3 文档收尾(总案 §8 S3)

- [ ] S3.1 新树指导文件:src/CLAUDE.md(构建/运行/判据程序约定,承接旧树 CLAUDE.md 的角色;
  根 AGENTS.md 已有的内容引用不复制)。
- [ ] S3.2 src/EXTENDING.md:扩展点按现实落地写——pes_* 体系成员(systems/)、dist_* 成员
  (methods/ + 薄 reg_*)、终止原语(pred_reg_prim 两轨)、结果分类行、积分器注册(reg_prop 闭族);
  各附"新增一件 = 新文件 + 薄 reg、零触碰现有文件"的 drag-in 演练步骤。
- [ ] S3.3 src/README.md 收尾态更新(状态行 + 与 EXTENDING/CLAUDE 的衔接)。
- [ ] S3.4 design_spec_completeness.md §4.2 回填:接口 A/B′ 现行正式签名(force_eval /
  elec_prop / mqc_statistic + container_term;samp_reg/dist_proc_i)按落地形态书写;
  退化点良定检查记录(adiabatic 全 no-op 退化成员 / 表面范式单片段跳过);双实现验证记录
  按现实载体(A 轴三参考势位级 + 双编译器位级;B′ 轴五 case + probe 驱动覆盖的成员面;
  配置面气相/表面双范式;电子轴 adiabatic 退化)。
- [ ] S3.5 optimal_structure.md 状态行 → "implemented (new program)";根 README.md layout 更新
  (src/ 新程序 + venus2026/ 只读参考;移除过时的迁移注记)。
- [ ] S3.6 提交文档批(含根 AGENTS.md)。

### S4 外部对账(总案 §8 S4)

- [ ] S4.1 Notion "VENUS-TDHF 泛化" 任务勾选态同步——外部系统,无程序内通道;
  本件登记为**用户侧条目**,S 收口时提醒执行。

### T1-T7 结构判据逐项核验(总案 §1;证据入本件验收区)

- [ ] T1 文件结构:对照 src/README.md 现行层表逐文件在场(总案 §1 表的字面名按 E12/E14/E18/E19
  漂移注记:predicate 併入容器槽、analysis/ 并入 output/、L3 成员分居 methods/+systems/)。
- [ ] T2 零触碰扩展演练:复制一份参考势成员为变体(如 pes_leps 变体)拖入,make 全绿,
  git diff 证明除新增文件外零改动;演练件随后撤销,证据(命令 + diff 摘要)入本件。
- [ ] T3 双实现矩阵:按上节"现状事实"第三条的姿态清点现实载体,逐轴记录覆盖证据。
- [ ] T4 单位契约:全库 grep 23.0605 / 0.04184 —— 仅 force_facade 换算点 + consts 立法定义点。
- [ ] T5 单随机流:全库 grep RANDOM_NUMBER —— 仅 rng.f90 内部转换点。
- [ ] T6 三相:driver 审读 + grep——装配→演化→统计调用序,演化循环内仅 facade/propagator/recorder。
- [ ] T7 = S2(新审计产出即为 T7 证据)。

### 停放项终局清点(V 移交;逐条裁决:维持停放/升级后续件/关闭)

- [ ] P-1 装配期频率导出通道未接线(V-5;四证据在案)。
- [ ] P-2 ifx 下 check_facade 两红腿((g) 1 ulp 构造性不等价;(r) eig_sym 算术序 1.29e-15 既有停放族)。
- [ ] P-3 采样器→统计相初值快照通道(VI/OAMI/AMAI/AMBI/SDA/SDB 生产者,解锁 fin%ang(4..18))。
- [ ] P-4 scripts/ 改进注记:check_facade (g) FAIL 行不打印 got-vs-expected。
- [ ] P-5 C8 终审停放区六条(hist_thta 气相门控候选 / col_occ 守卫 / fin_t 初值参考量 / SDA/SDB /
  plot_hist allocated 守卫 / lattice%q_ref)。
- [ ] P-6 C2 停放区五条((r) 1e-15 性状 / fssh / iesh 后窗件 / r_adiab 构造通道 / 有状态成员 wrapper 终形)。
- [ ] P-7 b_layer NU_COLL 携带无消费面注记。
- 清点产物:本件"停放项终局表"——每条一行:登记处 / 触发条件 / S 相裁决(维持/升级/关闭)/ 去向。

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
