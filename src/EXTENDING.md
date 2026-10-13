# EXTENDING.md — 新程序扩展点(E27 后契约)

接一个自定义势能面（含可选非绝热方法包）的**最短路径是 Tier 1**：一个文件夹、
若干 .f90、零 reg 文件。需要更复杂的接线（多方法包、自定义记录列、复杂参数体）
时升级 Tier 2（§1b 的 reg 三件套——reg 在场则自动接线静默退出，pot 侧预留名
函数全部可被 reg 直接复用，升级零迁移）。

## 1a. Tier 1 纯势容器——十预留名词表（封闭集合）

`systems/<sys>/` 文件夹**不含 reg_\*.f90** 即 Tier 1 容器。文件夹内 .f90 暴露
下列预留名（名字即约定，签名即既有接口，编译期核对；裁决记录
docs/plans/2026-10-12-auto-pes-container.md）：

**势能面组**

| 预留名 | 形 | 必备性 | 说明 |
|---|---|---|---|
| `pes(q, mass, g, v)` | subroutine | 恰一个 | E27 力交付：q[A]/mass[amu] 进（q/g 一维展平 size=3×nat），g[eV/Å]/v[eV] 全赋值出 |
| `keys()` | function | 可选 | `character(len=16)` 数组容器键词表（parse 前声明） |
| `load()` | subroutine | 可选 | 无参自持 `pull`；**亦即 seam 期初始化钩子**——读参数文件（NN 权重模式）/分配在此，先于一切力消费者与 init |
| `term(sta)` | function | 可选 | 终止判断（键缺席门控用 found 旗标自实现；未绑 = 跑满 MAX_STEPS） |
| `classify(fin)` | function | 可选 | 通道分类（返回通道名） |

**方法包组**（自定义非绝热方法）

| 预留名 | 形 | 必备性 | 说明 |
|---|---|---|---|
| `prop(sta, dt)` | subroutine | 与 match 成对 | 电子传播成员（`elec_prop_i` 同型） |
| `match(sta)` | subroutine | 与 prop 成对 | 统计匹配成员（跳变、occ/n_hop 写回） |
| `sync(s)` | subroutine | 可选随对 | a↔rho 表象换算尾巴 |
| `init()` | subroutine | 可选随对 | 选中后初始化（门控归生成代码；参数自持走 keys/load） |
| `method_name()` | function | 随对必备 | 返回 `character(len=32)` **小写**注册名；`ELEC_METHOD = <名>` 选用 |

规则一句话：势能面组在场即绑；方法组注册无条件、初始化有门；reg 在场则十名全忽略
（reg 权威）；名字须 public；prebuild 对 pes 0/≥2、任一名 ≥2、形不符、prop 与
match 不成对、method_name 缺对、sync/init 孤立、双自动文件夹各有一条指名中止。

现成范本：`systems/inactive/auto_sample/`（pot.f90 / io.f90 / tdhf_stub.f90；
判据 `scripts/check_auto_pes.f90` + 演练 `scripts/drill_auto_pes.sh`）。
激活 = 文件夹搬回 systems/ 顶层（一次一容器；与 reference 并存会在装配期
STOP 1——测试时把 reference 暂存 systems/inactive/）。

## 1b. Tier 2 全容器——reg 三件套（体系专属代码的组合点）

一个文件夹 = 一个体系,自足材料:

```
systems/mysystem/
  my_pes.f90        体系物理: H 定义/势能/力律/随机律, eV 交付
  my_box.f90        (可选) 参数接收器(键表 + params_t) + 终止/分类判断
  reg_mysystem.f90  拖入件(接线区块模板, 2026-09-22): 唯一公共入口
                    reg_mysystem_slots() 一次交接填满 interface/container_sched
                    的八个槽(2026-10-05: 注册期谱导出槽废除); 区体全部私有, 按
                    declare_keys/load_params(输入), bind_force/register_method
                    (常规), bind_term/init_method(初始化),
                    bind_classify/register_columns(输出)分组——空槽合法跳过,
                    时序与守卫归 container_sched
  system/           反应物构型(2026-09-22 裁决): 一反应物一文件(.xyz 无晶胞 /
                    .poscar 带晶胞), 只放构型——扫描对未知扩展名即停(文件夹根的
                    .f90/引擎参数文件不进扫描); SYSTEM_DIR 指向此子文件夹
```

容器 reg 骨架(照抄即得, 现成范本 systems/reference/reg_reference.f90 与
systems/inactive/cau_rst/reg_cau_rst.f90):

```fortran
subroutine reg_mysystem_slots()          ! 唯一公共入口: 一次交接
   call container_sched_bind(keys    =declare_keys, &     ! 声明键词表(解析前)
                             load    =load_params, &     ! pull 暂存键 -> load
                             force   =bind_force, &    ! container_bind_pes(my_pes)
                             method  =register_method, &   ! method_reg(体系专属方法包, 可无)
                             term    =bind_term, &     ! 键门控 container_bind_term
                             init    =init_method, &     ! 方法包初始化(ELEC_METHOD 门控)
                             classify=bind_classify, & ! 键门控 container_bind_classify
                             columns =register_columns)    ! rec_reg_col(选择门控)
end subroutine
! ---- 以下为私有功能名区体: 每体只做自己的接线 + 自己的门控 ----
```

- **构型文件随容器**(2026-09-22 裁决): 体系专属容器的反应物构型放容器文件夹的
  `system/` 子文件夹——扫描纯净性所迫(逐文件判据的即停不容代码/参数文件混入),
  同时"拖入一个文件夹,体系即完整定义"成立。名字排序即 list_atoms 序:C+Slab → C 先
  (cau 引擎布局: 原子 1 = 入射 C,其后 slab),F+H2 → F 先(fh2 聚合入口要求
  F,H,H),均天然吻合。cau_emtnn/cau_rst 的 `system/Slab.xyz` 双重角色: 引擎按
  `SLAB_FILE` 键读它,sysdef 把它作靶片段定义(判据程序路径常量已随迁)。
  参考容器(纯数学)不带构型——五案例的案例侧 system/ 照旧是构型来源。
- **力交付契约(E27)**:`subroutine my_pes(q, mass, g, v)`——q[A]/mass[amu] 进,
  组合完成的受力 `g = dV/dq`[eV/A] 与势能 `v`[eV] **全赋值**出;换算归力接口,勿自换。
  多面体系自折活性面(对角化+占据折叠;通用助手登记为候选后续件,git 历史可取
  force_interface 退役前实现)。
- 随机律成员在容器内部和(g_add 语义已内化),rng 单流直用。
- 体系专属**电子方法包**(如自定义 Hubbard 传播)同文件夹放置:prop/match 成员 +
  `reg_*` 里一行 `method_reg('my_tsh', my_prop, my_match)`;成员与 H 模块共享私有数据
  (E27:无互换工作区);输入键 `ELEC_METHOD = my_tsh` 选用(绝热参考跑仍是
  `ELEC_METHOD = adiabatic` 或缺省)。
- 参数经程序输入的**容器键**(E28, 运行目录参数文件通道退役): box 持键表 `my_keys()` 与
参数接收器 `my_params_t`(字段缺省即容器缺省, found 哈尔存在性); reg 的 `declare_keys` 区体
声明词表(`input_declare_keys`——prebuild 聚合为 `registry_gen_keys`, 驱动器在解析前调用,
封闭文法由此放行这些键名), `load_params` 区体一行一键 `pull` 进接收器(缺键保留缺省), 交
`my_load(p)` 收进模块态并做值校验; `input_audit_declared` 兜底"buffered 而无人消费"。
容器键不得与程序键/成员键撞名(声明即中止)。参照 systems/reference/reference_box.f90。
- **势能面扫描成像**(`scripts/probe_pes.f90`,独立工具,判据程序同族,不进构建树):
  通用扫描引擎经一个 `pes_hook` 薄适配文件绑定容器的势能面(判据程序同款 pot 层入口,
  与生产链同一条数值代码路径);按 `pes_scan.txt`(KEYWORD = VALUE)沿 1~2 根轴
  (BOND 键长 / CART 绝对坐标)扫 V 与投影解析斜率,中心差分对照自证,
  `plot_pes.py`(python+matplotlib)出 PNG;基几何/质量可取自容器 `system/` 的构型
  值。现成 hook:reference / fh2_leps / h3_bkmp2 / h2ag111(编译配方见 probe 头部;
  cau_* 平均场容器绑 `cau_<sys>_ground` 基态面,文件集见 check_cau_* 配方)。
- **容器不立法谱**(2026-10-05 谱来源立法):容器不再携带注册期谱导出槽——谱的
  生产者是程序内部的 Hessian 推导(SPECTRUM_SOURCE=COMPUTE,缺省)或运行输入的
  谱表文件(SPECTRUM_SOURCE=MANUAL + SPECTRUM_FILE),见 §6。

## 2. 初始态分布成员(methods/samp_*.f90,B′ 族)

新成员 = `methods/samp_mydist.f90`(init/draw/realize 两段式 + 报错即停条目族)+ 薄
`methods/reg_samp_mydist.f90`(`samp_reg` 一行 + `_init` 初始化组装机构消费成员键);
`sampler` 的 scheme 码表加一词(scheme_code_name——这是族表数据非逻辑分支)。
(2026-09-28: 族前缀 dist_→samp_ 更名;入射通道不在本族——它按组装范式拆为
inc_surface/inc_pair/inc_single 三成员 + beam_laws 共享组件,见
docs/plans/2026-09-28-incident-split-surface-paradigm.md。)
零体系数据:质量/几何取 list_atoms,参数走输入键。
成员按多载体章程实现(2026-09-21 决定):同方案词可选多个片段——init 绑定载体
列表(仅零载体终止),sample 一次调用循环全部载体,**一份参数平等服务每个载体**;
装配数据类字段(频谱/模式表)逐载体一表(reg 缝按载体拉取)。

## 3. 通用方法包与组件(methods/,库章程)

- 通用包:elec_*.f90 + reg(行进方法表,任何体系可经 ELEC_METHOD 选用)。
- 组件(level 2):amp_*(振幅积分器)/ fssh_g_ij / hop_decide / momentum_rescale /
  dom_reset 计划位——容器或成员经 `use` 组合,物理单拷贝。

## 4. 积分器(闭族,methods/reg_prop.f90)

verlet/symple/radau 固定成员制;新增积分器 = dynamics/propagator.f90 内实现 +
reg_prop 注册行(闭族非开放注册面——成员资格是决定事项,见 E 编号决定记录
docs/plans/K2.md、S.md)。

## 5. 记录列与后处理

- 记录列:output/recorder.f90 列注册表加行(rec_reg_col;列提取器纯形式)。
- 结果分类:物理系统容器侧 container_bind_classify(E19 槽)。
- 直方图类:output/plot_hist.f90 五类族(数据驱动)。

## 6. 谱来源(SPECTRUM_SOURCE,2026-10-05 立法)

采样成员(wigner/normalmode/boltzmann 拉模表, ebk 拉双原子常数)消费的谱由
**一个键独占选择的生产者**提供(模块 interface/spectrum_interface.f90;容器导出
通道已废除):

- `SPECTRUM_SOURCE = COMPUTE`(缺省):内部 Hessian 推导,范式中立——表面范式
  与气相同一推导路径,裁判是跨块衰减门(PES 在探针处的物理行为),不是范式标签。
  单片段体系的探针=缓冲几何原样;多片段按延展感知间距沿 x 排开(原子间最小距离
  ≥ r_asym=12 Å 量级)。EBK 双原子常数通道无推导(非谐/转动常数是文献数据),
  COMPUTE 下拉即命名中止并指引 MANUAL。
- `SPECTRUM_SOURCE = MANUAL`:`SPECTRUM_FILE` 谱表文件是唯一生产者——被采样
  载体在文件中必须有记录(缺记录/空文件/文件带错即命名中止,"参数作为输入要求")。
  `SPECTRUM_FILE` 在场而 SPECTRUM_SOURCE=COMPUTE 是输入层形式错误(无静默忽略)。

谱表文件格式(权威规格:docs/plans/2026-10-05-spectrum-source.md §3):行导向、
大小写不敏感、`#`/`!` 整行注释;每条记录 `FRAG`(成分键,与 list_atoms 片段
逐符号大小写不敏感匹配)/`NMODE`(nat=2→1;nat≥3→3nat-5 或 3nat-6)/
`WVN`(频率[cm⁻¹],每模一行,严格升序)/`COL j`(其后一行恰 3nat 个质量加权
正交归一系数)/可选 `DIAT w_e w_e·x_e B_e`(仅 nat=2,EBK 通道)。

**模列参考架契约**:列向量表达在**片段缓冲几何(容器 system/ 构型文件)的笛卡尔
架**中;外部量化程序输出的模列须作者自行旋转到该架(不做自动对齐)。生成器范本
`scripts/gen_spectrum_tbl.f90`(双原子记录,列出自 stretch_mode 本体)。

## 演练清单(新扩展件落地后)

1. `cd src && make`(如新增 reg:再 make 一次)→ 双编译器绿;
2. 相关 scripts/check_* 重编重跑;
3. case 回归走值域判据 `scripts/check_ev_regression.sh`(守恒漂移报告、
   初帧能量锚对 baseline 旧制换算 ≤1e-3 相对差、直方图分布比对);
   基线两代:baseline_kcal_v1(旧制,kcal/mol 门面)与 baseline_ev_v2(eV 制)。
   逐字节等价门已随 2026-10-06 单位改制(eV/Å/amu/fs)废止——浮点舍入路径
   随单位制改变,位级复现自此不再成立,回归以值域判据为准;
4. T4/T5 grep 复验(23.0605/0.04184/349.755/0.063508/0.083144 在 src/ 零命中;
   RANDOM_NUMBER 零命中)。
