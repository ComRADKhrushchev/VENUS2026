# T 相——E28 容器参数通道并入程序输入（2026-09-22，用户裁决"按形态A实现"）

> 动因：用户对"一次运行两个输入文件"（input_qct.txt + *_pars.txt）不满。前此多轮讨论确立
> 方案：容器参数改走程序输入的**声明容器键**通道，消费侧用形态 A（一行一键 pull）。
> 本文件是该裁决的记录与验收证据；架构裁决依据 docs/design_spec_completeness.md §4.2
> 之 E28 括注。

## 1. 裁决内容（E28）

1. **运行目录参数文件通道退役**（`reference_pars.txt` / `cau_*_pars.txt` / `fh2_pars.txt` /
   `bkmp2_pars.txt` / `h2ag111_pars.txt` 全部撤销）。容器参数 = input_qct.txt 的声明容器键。
2. **准入保持解析期封闭文法**（E21/E22 立场不变）：input.f90 新增运行期声明键集
   `input_declare_keys`——键合法 ⇔ 程序键 ∪ 成员派发键 ∪ 声明键；表外键照旧解析期命名中止。
   声明键不得与程序键/成员键撞名（声明即中止；cau_emtnn 的 `A_LAT` 撞表面窗口程序键，改名
   `SLAB_A_LAT`——本批唯一键名变更）。
3. **装配序新增纯词表阶段**：prebuild 探针聚合各 reg 的可选 `<reg>_keys()` 入口为
   `registry_gen_keys()`（探针手法与 `_init` 同款）；driver 在 read_input **之前**调用它，
   `_all`/`_init` 时序不变（h2ag111 的谱导出在 `_all` 期读 roster，位置不动）。
4. **消费侧（形态 A）**：`pull(key, target [, found])` 泛型（real/int/char/real(8)(3)）居
   input.f90——缺键保留目标字段缺省，found 哈出可选槽存在性（term_on/class_on 语义），
   坏值命名中止（键+行号）；未声明键的 pull 与双 pull 均为容器作者错误命名中止。
5. **消费审计**：driver 在 `registry_gen_init` 后调 `input_audit_declared`——staged 而无人
   消费的声明键命名中止（兜住"声明表与 pull 字面量漂移"与"SYSTEM_DIR 未暂存而键已写"两类
   失配）。
6. **box 形态**：`*_params_t`（字段缺省=容器缺省+found 旗）+ `*_keys()` 键表 + `*_load(p)`
   收进模块态并做值校验；各 box 的私有文件解析器（split_kv/rdreal/stop_dup 族）删除。
   reg 文件三个可选入口：`_keys()`（声明）/`_all()`（绑槽+方法行，不变）/`_init()`
   （guard→pull→load→绑可选槽）。

## 2. 落地清单

- **src/**：input.f90（声明集/准入/pull/审计/容量+非空校验）；driver.f90（装配序
  registry_gen_keys→read_input→…→registry_gen_all→registry_gen_init→audit→container_init，
  步骤号与契约注释同步）；scripts/prebuild.sh（probe_keys + registry_gen_keys 发射）。
  Makefile 零改动（keys 聚合居 registry_gen 同一文件）。
- **容器（8）**：reference、cau_emt2、cau_emtnn、cau_peri、cau_rst、fh2_leps、h2ag111PES、
  h3_bkmp2——box 重写参数通道、reg 增 `_keys()` 与 pull 化 `_init()`。
- **cases（8 运行目录）**：pars 键并入 input_qct.txt（带"容器键"分节注释），参数文件删除，
  README 提及处同步。
- **判据**：check_input 新增 (r1)–(r5)（声明→准入→typed pull→审计干净；staged-未消费/
  未声明 pull/双 pull/空值四中止腿）；check_facade (o) 腿重写为 obox_input.txt 通道
  （缺键缺省/好键 staged/坏输入四路各归其主：read_input/pull/reference_load）；
  check_final_state (d) 参数构造化；check_driver (j1)/(k)/(e)/(r) 装配序镜像+审计+case_t
  键并入。随批吸收两处**既有漂移**（非本批所致，树中 wigner 成员使 members=11→12：
  (a)/(r)/(j1) 三处断言与 stdout 期望对齐 12；判据配方注释补 dist_wigner/spectrum_stage/
  spectrum/hessian 模块行）；check_cau_rst 配方注释的行接损伤（cau_rstst_hubbard_pes.f90）
  顺手修复。
- **文档**：EXTENDING.md 容器节（参数通道条款+reg 形态）；AGENTS.md 架构规则区 E28 段；
  design_spec §4.2 E28 括注；glossary E 编号行与"容器键"术语；本文件。

## 3. 验收证据（2026-09-22）

1. **双编译器**：`make`（gfortran）与 `make FC=ifx` 全量构建绿，0 error。
2. **拖入冒烟（7/7）**：七个 inactive 容器逐个拖入 systems/ 顶层，第二次 make 编译链接全绿
   （首次失败为 depend.mk 单轮迟滞，房规在案），随后拖回恢复。
3. **判据全绿**：check_input 0 fail（含 (r1)–(r5)）；check_facade 0 fail（(o) 六腿新通道）；
   check_driver ALL PASS；check_final_state ALL PASS；check_cau_emt2/emtnn/peri/rst =
   99+94+94+94 PASS / 0 FAIL（与基线完全一致）。
4. **八运行目录终跑**：全部 rc=0、stderr 空；有基线产物的六目录 trajectories.txt 与两直方图
   **逐字节一致**（数值链零扰动）。
5. **stdout 对档**：与入库 results.txt 的差异仅 members=11→12（上述既有漂移）与 timing 行
   （档案自注非锚）；其余判据行逐字一致。results.txt 档案未改动——漂移登记于此，
   重入档属相收口事务。

## 4. 交付物

- 生产：src/input.f90、src/driver.f90、scripts/prebuild.sh、systems/reference/ 三件、
  systems/inactive/ 七容器 box+reg 对；
- cases：八 input_qct.txt（含容器键）、删八参数文件、README 六处提及；
- 判据：check_input/check_facade/check_driver/check_final_state/check_cau_rst 五件适配，
  check_input 新五腿；
- 文档：EXTENDING.md、AGENTS.md、docs/design_spec_completeness.md、docs/glossary.md、
  docs/plans/T.md（本件）。

## 5. 状态

已完成（用户裁决 2026-09-22"按形态A实现"；执行同日）。悬挂：results.txt 档案的
members=11→12 重入档；未提交工作树中与本批并行的 A+B/wigner 批（其判据配方缺口已随批补齐）。

## 6. 增补（同日第二批）——cau_peri 判据下标裁决 + 两项核实

用户三项（2026-09-22）：

1. **文献核实（h2ag111 刚性）**：归档引擎（h2ag111.f：两 H 原子 6 自由度，q(3)/q(6) 为
   "相对 Ag(111) 表面的高度"）即 Jiang & Guo 2014 JCP 之 PIP-NN 方法对 H2+Ag(111) DFT
   数据的拟合（symmetrize→npoly=13 不变式→NN getpot 结构吻合）——6 维拟合以刚性板为
   构造性前提；板弛豫/声子修正属该领域后续独立课题，不在此 PES 内。结论：容器
   system/ 仅含 H2.poscar 是正确的（表面非 roster 片段，表面范式 TARGET 非法）。
2. **system/ 缺构型行为实测**（断链前 venus.e）：空目录 → 命名中止 "no system-definition
   files in the container folder"（exit 1）；目录缺失 → "folder not found"（exit 1）。
   无静默路径。
3. **cau_peri 判据下标裁决落地**（用户裁定"改判据下标"）：
   - roster 采家族命名排序约定 [C, 锚 S, ghost Y]（C.xyz 先于 Slab.xyz）；引擎内部布局
     仍 [S, Y, C]（周期核以 q(1:3) 为主原子）。
   - box 新增边界置换包装器 `cau_peri_pes_r`（q/mass 置换→引擎→梯度回置）与
     `cau_peri_prop_r`（局部几何副本置换；rho 生命周期在引擎模块态，sta 原样传递）；
     match/sync/col_* 不读 q，直连不动。
   - 判据下标：`cau_peri_term` 由 q(9)−q(3) 改为 **q(3)−q(6)**（roster 序：原子 1=C、
     原子 2=锚）。reg 的 zone_force/zone_method 改绑包装器。
   - 反应物文件补齐：system/C.xyz（C 于原点）+ system/Slab.xyz（Au 锚原点 + Au ghost
     名义 0.9 Å 振子偏移，待产品化调参——注释在案）。
   - 判据新增 (l)：roster 序包装器对 permuted 几何复现引擎序直接调用的 v 与 g **逐位
     一致**。check_cau_peri 重编重跑 **95 PASS / 0 FAIL**（基线 94+新腿）。

**并行会话注记**：本增补执行期间另一会话正在树上做 interface 层改名（facade→interface）
与 sampler draw 接口重构，`make` 处于中段红态；reg_cau_peri 的编译验证（依赖被改模块）
待该批落地后补跑，此为已知悬挂项。cau_peri 判据编译不依赖被改模块，验证有效。

## 7. 增补（同日第三批）——E29 BOXPAD 键与簇盒投放退役（用户裁决"团簇碰撞功能不需要"）

路线 (a)：键与能力同生共死，无悬空。改动：

1. **src/input.f90**：键表 31→30（BOXPAD 行删除；A_LAT 成为唯一表面窗口键）；
   map_slots 分组例去 BOXPAD。
2. **src/globals/config_atoms.f90**：list_atoms%box_pad 成员、两处复位、窗口对账的
   BOXPAD 双向条款（分子范式拒斥 + 表面范式必填/拉取）全删——window_reconcile 收窄为
   A_LAT 单条款。
3. **src/methods/dist_incident.f90**：投放点模式表收窄为 0（原点）/1（胞内均匀）；
   n_aim 合法域 0..1；box_pad 正值守卫、pad_c 通道、case(2) 抽取分支删除。
4. **判据**：check_input (o3)/(o5)/(o7) 三腿退役、余腿重排 (o1)-(o4)、模板与断言同步；
   check_roster 模板/断言/注释六处；check_dist_incident (i) 簇盒分布腿与 (s) 零半宽
   中止腿退役（bad-pad 子模式删）。
5. **cases**：case_surface_window 输入删 BOXPAD 行、头注与 README 五处改写（含
   退役注记条目改写）。
6. **数值面**：N_AIM=1 下 box_pad 从不进算术——case 重跑逐字节一致可验。

### E29 验收证据

- 双编译器构建绿；
- check_input / check_roster / check_dist_incident 重编重跑全绿；
- case_surface_window 重跑：rc=0、trajectories.txt 与两直方图逐字节一致、stdout
  除 timing 外逐字一致（装配行 members= 为 E28 增补已注记之既有漂移口径）。

## 8. 增补（2026-09-24）——probe_seq 坐标序列 PES 评估小工具（C-Au 四容器通用）

用户需求：读一个 C 原子坐标序列，逐点输出基态能量与位置依赖 Hubbard 参数。形态沿
probe_pes 的 hook 思想，全部落在 scripts/，零触碰 src/、systems/、cases/：

1. **scripts/probe_seq.f90**：~150 行通用驱动。序列文件每行一个 C 位置 `x y z` [Å]
   （`#` 注释/空行跳过，坏行命名中止带行号）；逐点调 hook 输出
   `i x y z V h(1) h(2) U J Vcp` 表 + V 最小/最大摘要（仅 stdout）。
2. **hook 契约扩展**（translation unit `pes_hook`）：pes_desc / pes_init（容器加载
   链）/ pes_geom(nat,names,mass,q0,ic)（引擎序基准构型，slab 取容器自带
   system/Slab.xyz + C 下标）/ pes_entry（ground 面）/ pes_hub（hubbard_val 面）。
3. **四个 hook**：hook_cau_emt2（scripts/ 锚点参数副本 + box 缺省胞
   17.7×17.7/60°）、hook_cau_emtnn（box 缺省 2.88/6 + 货架三件）、hook_cau_rst
   （slab/weights/BvK 货架路径）、hook_cau_peri（拷 Analytic_Potential.txt 到运行
   目录——READAP 固定文件名；直写引擎公共变量 k_surf/k_ghost/Q0_ghost = box 缺省；
   引擎序 [S,Y,C]，ic=3）。
4. **配方**（各链 引擎+pot+hook+probe_seq+specio 五件，pot 层不依赖 consts/
   state/linalg/rng）：见 probe_seq.f90 头部四行。

验收：四 exe 编译绿（-Wall 自检：新文件唯一告警为契约形参 mass 未使用，预期内）；
四容器 top 位下降序列逐点全列有限、沿路径平滑、物理形态合理（peri 的 h(1)≠h(2)
晶场劈裂可见，emt2/peri 的 Vcp 随高度指数衰减，emtnn 相互作用能远高趋 0，emtnn 于
z=2.0 有 −14.3 eV 吸附阱）；缺序列文件/坏行/空序列均 fatal 带 STOP 1。工具级新增
（scripts/ 五文件），不入判据体系；产物 exe 用后即删（配方一行重建）。
   - 作图验证（同日）：scripts/seqscan/ 编排脚本从各容器自带 slab 算 top/bridge/hollow
高对称点（peri 无晶格，取名义 d=2.88 的 +x 横向偏移），z=1.0..6.0 每 0.05 一点，每容器一张
六面板能力曲线图（V/h(1)/h(2)/U/J/Vcp）；四图经视觉验收通过（含一轮 ylim 裁剪缺陷修正：
面板取值域须三曲线汇总后统一定限）。物理面：远高渐近跨位形一致（emt2 32.59 / rst 969.61 /
emtnn ≈0），h(1)=h(2) 于 emt2/emtnn/rst 简并、peri 可见晶场劈裂。
   - 迁移（同日，用户裁决"C-Au 工具归 C-Au 目录"）：整套套件（probe_seq.f90 + 四 hook +
seqscan 编排与图 + emt_2spec_params.txt 副本 + README）移出本仓库至 MAIN/3/probe_seq/，
shelf/配方/REPO 路径常量改指 ../../VENUS2026；迁移后四容器编译+冒烟+全量曲线复跑数值
逐位一致。scripts/ 下本套件五文件不再保留（emt_2spec_params.txt 仓库原版仍留 scripts/ 供
check_cau_emt2 配方用）。

## 9. 增补（2026-09-30）——cau_emtnn 权重文件勘误替换（用户裁决"做一下权重文件替换"）

- **发现**：probe_seq 六面板曲线复检（本日，套件 exe 先按引擎现状重编——旧 exe 落后引擎
  源码两小时）暴露 cau_emtnn 的 Hubbard 参数曲线全面病态：J 随高度从 +1.36 衰减、z≈5.6
  过零、远处饱和 -0.49 eV；Vcp 远处收敛 -3.24 而非 0；h/U 停在 -26.4/+15.5 错平台
  （训练目标：J 为 Hund 耦合全程 ≈+1.23 平直、Vcp 衰减归零、h@6.25=-15.9、U@6.25=+10.3）。
  在训练网格 311 原坐标逐点对照：h RMSE 22.2 eV / U 18.1 / Vcp 4.5 / J 0.55。
- **根因**：容器 2026-09-22 建立时拷入的是旧项目 2026-08-20 04:21 的过期快照（与
  `MAIN/3/C@Au_TDHF/2D/_rmse_check/nn_weights_emt_nn.txt` MD5 逐位一致），早于当日 06:09
  的 Hubbard s_ad 分支仿射重标尺修复（au-placeholder 标尺 η₂=3.197/s₀=1.642 → 本文件 EMT
  标尺 η₂=5.0967/s₀=0.7160，`_affine_fix_hubbard.py`）。正确终稿为
  `MAIN/3/AnalyticModel/fit_emt_nn/final_weights.txt`（08-20 06:10，含 Ecoh_NN 25 权重与
  尾部仿射注释，09:51 出全套终验图）；当日 04-06 时段多套候选权重（4b_ckpt16 / e0fit_full /
  pre-affine）在旧项目流转，04:21 那份恰是未修复版。勘误过程同时修正本档早前一处误判：
  "旧引擎同权重对照分岔→疑移植缺陷"不成立——`Hubbard_at_data.dat` 出自 e0fit_full 权重
  （`fit_config.txt` 指定），并非容器那份；容器引擎移植本身经 final_weights 验证无罪。
- **替换**：`systems/inactive/cau_emtnn/nn_weights_emt_nn.txt` ← final_weights.txt 整文件
  覆盖。加载器兼容性预核（多键 NH 头行 / `## Ecoh_NN` 节 / 密度块注释盲跳 / `pairwise`
  关键词检测 / `# nh_hubb:` 交叉核对均匹配）。旧权重存档于上述 _rmse_check 路径（md5 核实）。
- **验证四腿**：①加载横幅全链——final 参数集 5.0967…、Ecoh NN ACTIVE（25 权重）、Hubbard
  s_ad+pairwise 双块、NH_ecoh=8；②probe_seq 三位形复扫——J 全程 +1.21~+1.28（z=2~6；10 Å
  训练域外外推 +1.32，有界）、Vcp 物理窗归零（z=4~6 处 |≤0.26|）、h@6=-15.68 / U@6=+10.3
  （目标 -15.6/+10.3）、h(1)=h(2) 严格简并、全列有限、z∈[2.5,8] 最大二阶差分 ≤0.13 eV，
  h 通道对训练点 RMSE 0.14 eV；③`check_cau_emtnn` 重编重跑 94 PASS / 0 FAIL（i2 armed
  漂移 -6.5 kcal/mol，界 200）；④seqscan 四容器六面板图刷新（其余三容器权重未动，输出
  等价）。地面随新权重改变：V 井 -15.0~-14.7@z1.7-2.1 → -13.7/-12.5/-11.7@z2.0-2.2，
  远渐近 -0.39 → +2.88 eV（E0_w=9.28 计入）——属预期（旧 V 亦为旧地面产物）。
- **生产流不受影响**：cau_emtnn 居 inactive 货架（build glob 两层之下），reference 容器
  与五 case 零触碰，无需重跑入库 results。
- **后续候选件**：`AnalyticModel/fit_emt_nn/nn_weights_gsfit_full.txt`（08-30，基动能面
  对 DFT 参考再拟合；HUBBARD_NN 八块与 final_weights 逐字节相同，仅 EMT_params/Pair_NN/
  Ecoh_NN 不同）——地面若需采纳为独立裁决事项。
