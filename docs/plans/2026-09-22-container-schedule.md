# 2026-09-22 容器装配时刻表 + 接口族/列表/buffer 改名（用户决定）

状态：已落地（本文件即决定记录；执行与验证证据见文末）。

## 1. 决定链（用户三项决定，同日连续）

1. **容器 reg 接线区块模板**：容器 reg 文件按接线功能分四个区块——输入接线（词表声明 +
   键拉取装载）、常规接线（力槽、方法行）、初始化接线（终止槽、方法武装）、输出接线
   （分类槽、记录列、注册期导出）。
2. **时刻表上收通用代码**：三入口的调用次序、守卫、空槽跳过不再由每个容器 reg 自持，
   单源化到通用模块 `interface/container_sched.f90`；容器侧只剩"内容 + 一次交接"。
3. **改名包**（覆盖既有登记条目）：
   - `facade/` → `interface/`；模块 force_facade / elec_facade / sysdef_facade →
     force_interface / elec_interface / sysdef_interface（规范中文"接口"不变，英文载体
     与中文对齐；词汇表"标识符不变"旧注就此覆盖）；
   - roster 族 → list_atoms 族：实例/类型 roster / roster_t → list_atoms /
     list_atoms_t，roster_load / roster_clear → list_atoms_load / list_atoms_clear，
     模块 config_roster → config_atoms（文件随改）；
   - staging 族 → buffer 族：staging_*（input 行缓冲 / sysdef 交付缓冲共 15 个标识符）
     → buffer_*；注释散文 staged/staging → buffered/buffer；模块 spectrum_stage →
     spectrum_export，入口 spectrum_export_stage → spectrum_export_put；
   - 交接动词 wire（提案暂名一）→ zone（提案暂名二，用户同日撤销——"随口说的"）→ 终形
     **按实际功能命名区块**：declare_keys / load_params / bind_force /
     register_method / export_spectrum / bind_term / arm_method / bind_classify / register_columns，
     签名接口 slot_i（槽位语义承仓库既有"力槽/终止槽"词汇），交接入口
     container_sched_bind，容器侧唯一公共入口 reg_<sys>_slots。

## 2. 设计要点

- **container_sched（九槽 + 守卫，写一次）**：无参抽象接口 slot_i；九个过程指针槽
  （nopass，缺省 null，槽即文档）；`container_sched_bind` 一次交接（关键字清单即容器
  接线菜单；槽单占——二次填充命名中止，与 container_bind_pes 同语义，一容器一进程）。
  生命周期：`container_sched_keys`（解析前跑 declare_keys）；`container_sched_all`
  （注册期跑 bind_force → register_method → export_spectrum）；`container_sched_init`
  （初始化组装机构：守卫一次——system_dir 空 → 无组装区运行；然后固定次序 load_params →
  bind_term → bind_classify → arm_method → register_columns）。
- **容器侧形态**：唯一公共入口 `reg_<sys>_slots()`（私有功能名区体按四个接线区块排布——
  输入 / 常规 / 初始化 / 输出）。
  守卫、时序、空槽跳过全部离开容器；键门控/方法选择门控留在区体内。
- **prebuild 契约**：methods reg 保持 `<mod>_all`（+可选 `<mod>_init`）；systems reg
  探测 `<mod>_slots`；聚合三入口名不变（driver 零改动）。registry_gen_keys =
  [容器 slots 交接] + [container_sched_keys]；registry_gen_all = [methods _all] +
  [container_sched_all]；registry_gen_init = [methods _init] + [container_sched_init]。
- **调用顺序保真**（逐容器核对）：cau 系 pulls→load→term→setup→五列 = load→term→
  arm→columns；四甲型 load→term→classify；h2ag111 力槽→谱导出；记录列五行注册次序
  原位保留（列顺序=注册顺序）。

## 3. 落地清单

- 新文件：`src/interface/container_sched.f90`。
- 重写 8 个容器 reg（reference + fh2_leps + h3_bkmp2 + h2ag111 + cau_emt2/emtnn/peri/rst，
  接线区块模板）。
- 重写 `scripts/prebuild.sh` 容器探测段与聚合发射；改名 sweep 覆盖 src/scripts/systems/
  cases 全部 .f90 与 .sh（133 文件）+ 活文档（AGENTS/src README/EXTENDING/CLAUDE/根
  README/cases README）。
- 判据程序适配：check_driver / check_facade 配方补 interface/container_sched.f90；
  check_facade 两子进程与 scripts/probe_chain_freq、cases/case_b_layer/probe_glo 补
  `registry_gen_keys()` 前置（driver 时序步 1，slots 交接在解析前完成）；9 个
  check_dist_* 配方补 utils/geometry.f90（取向随机化批次遗留的陈旧配方，顺手修复）；
  check_final_state (e) 期望词随消息对齐（'atom list'）。
- 档案随迁（改名/前批决定的直接后果，逐条在案）：四 case results.txt 的
  members=11→12（前批 dist_wigner 既有事实，check_driver 已按 12 校准而案例档案未随，
  本批随迁）；probe_glo 两行 "roster mass"→"list_atoms mass"、case_surface_window 一行
  abort 措辞、case_chain_frequencies 两行探针腿行——均为用户批准的标识符改名在输出
  字面量上的直接映射。
- 不改名清单（引用锚保真）：scripts/check_*.f90 文件名（check_facade / check_roster
  等为跨文档引用锚，内容已随 sweep）；venus2026/ 旧树；docs/plans/ 史层文档原词。

## 4. 验证四腿（全绿）

1. **双编译器**：make clean 后 gfortran、ifx 全量重建链接成功（含 depend.mk 单轮迟滞
   的第二次 ifx make）。
2. **判据程序全量重编重跑**：check_driver ALL PASS（端到端含新聚合链）；check_facade
   0 失败（(s) 位级比对仍 max|dF|=max|dE|=0）；check_input / check_roster /
   check_sampler / check_final_state / check_dist_* 十二件 / check_specio /
   check_sysdef / check_recorder / check_propagator / check_spectrum /
   check_geometry / check_hessian / check_linalg / check_rng 全绿；cau 四件
   94–99 PASS 0 FAIL。
3. **五 case 终跑**：八个子目录运行全 exit 0；全部锚行（assembly/statistics/outcome/
   drift/histograms）与入库 results.txt 一致（drift 值逐位相同——数值链未动）；
   probe_glo 输出与档案逐字节一致（--strip-trailing-cr 口径，档案 CRLF）；probe_chain_freq
   两腿 0 失败、腿行与档案一致。
4. **grep 审计**：src/scripts/systems/cases 内 roster/facade/staging 残留仅剩四个受保护
   判据名/短语（check_facade、check_roster 文件名锚与两处历史短语）；词汇表同步
   （接口/列表/buffer/装配时刻表四行更新或新增，含"标识符不变"旧注的覆盖说明）。

## 5. 已知边界

- slot_i 为无参接口，签名检查只到"无参"粒度（填错槽编译器不拦——与既有三入口形态
  风险持平；槽名即文档）。
- 第九槽（注册期谱导出）现仅 h2ag111 使用；新增接线类型 = 改 container_sched
  （封闭槽集，与积分器闭族同姿态）。

---

## 修订（2026-10-05 谱来源立法）

- **第九槽（注册期谱导出 export_spectrum）废除**：九槽 → **八槽**。谱的生产者改为
  SPECTRUM_SOURCE 键选择——内部 Hessian 推导（COMPUTE，缺省，范式中立）或
  MANUAL 谱表文件（interface/spectrum_interface.f90，合并原 spectrum_export 模块）。
  `container_sched_bind` 的 `export` 哑元删除；本文其余处对 export_spectrum/
  第九槽/九槽的描述为 2026-09-22 至 2026-10-05 间的史层形态。
- **registry_gen_init 内部次序**：`container_sched_init`（容器缝：参数 load）移至
  聚合首位——load 依赖型 PES（如 NN）必须在任何探其力的成员缝（谱推导）之前
  可计算。装配时刻表的四区块划分与守卫语义不变。
  见 docs/plans/2026-10-05-spectrum-source.md §2。
