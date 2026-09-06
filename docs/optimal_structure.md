# VENUS2026 最优程序结构与质料映射

状态：设计草案（自完备性规范演绎，不拘泥既有代码组织）
日期：2026-09-05
依据：design_spec_completeness.md（规范）、completeness_audit.md（质料清单 99 条目）
方法：先从标准演绎目标结构，再将现有子例程作为质料放置。放置结果分四类：
可放置（原样或轻改）、需改造后放置、无位置（出局）、结构缺质料（须新写）。

---

## 1. 目标结构（自规范演绎）

### 1.1 分层

```
L0 控制核        driver.f90          装配→演化→统计的调用序列，无物理
L1 物自体容器    state.f90           Q, P, s, F, t（可分配，维数=名册内容参数）
L1 表象层        config.f90          输入文法产物：场景、名册、分布、谓词、结局分类、观测量配置
L2 接口 A 门面   force_facade.f90    (Q,P,s) → (F,E)；内部组装 Σ力贡献成员 × 单选电子传播律
L2 接口 B′框架   sampler.f90         分布 → (Q0,P0,s0) 通用抽样器
L3 注册族成员    pes_*.f90 / dist_*.f90    可换实现体
L3 封闭族        propagator.f90 / predicate.f90    有限选择 + 文法解释器
L4 测量          recorder.f90        轨迹序列化 + 谓词求值（演化期驻留的全部）
L4 后处理        analysis/           结局分类执行、统计、绘图（可离线）
L5 仪器          rng.f90, linalg.f90, geometry.f90  随机流、对角化、旋转
```

要点：
- 每层只向下依赖，L3 成员互不依赖（注册族力学）；
- L2 门面是控制核看到的全部物理——换 PES、换方法、换场景都不触及 L0/L4；
- L5 仪器独立成层，纠正"gasdev 住在 PES 文件里"式的层次倒置；
- venus_data 单一全局块解体：物自体进 L1 state，表象进 L1 config，
  仪器状态（RNG 流）进 L5。

### 1.2 各层形式定义（对应规范条款）

| 层 | 规范身份 | 签名 | 变更性 |
|---|---|---|---|
| driver | 控制核 | 编排 | 不变 |
| state | 物自体 | 数据 | 不变（形式） |
| config | 表象 | 数据 | 不变（文法单调扩展） |
| force_facade | 接口 A 门面 | (Q,P,s)→(F,E) | 不变（框架代码） |
| pes_* | 接口 A 甲类成员 | Q→(V,∇V)，可叠加 | 注册开放 |
| elec_* | 接口 A 乙类成员 | (H,s)→(ds/dt,力规则)，单选 | 注册开放（绝热=退化成员） |
| dist_* | 接口 B′ 成员 | 表象→分布 | 注册开放 |
| sampler | 操作类框架 | 分布→(Q0,P0) | 不变 |
| propagator | 封闭族 | 力泛函→步进 | 封闭（加成员=注册行） |
| predicate | 封闭文法 | 状态→布尔 | 不变（解释器） |
| recorder | 测量 | 物自体流→记录 | 不变 |
| analysis | 后处理 | 记录→统计 | 追加律 |

---

## 2. 质料映射总表

判定符号：●=可放置（原样/轻改）；◐=需改造后放置；✕=无位置出局；
△=结构缺质料（须新写，括号给可借鉴的现有素材）。
出处行号省略，详见 completeness_audit.md。

### L0 控制核

| 质料 | 放置 | 说明 |
|---|---|---|
| VENUS.f90 主循环骨架 | ◐ | 剥离 E1-E4 清单所列混入物后，剩余编排即 driver 雏形 |
| （E1 结局块双份复制） | ◐ | 合并后成为统计相位入口 |

### L1 state / config

| 质料 | 放置 | 说明 |
|---|---|---|
| venus_data 的 Q/Q_old/P/PDOT/QDOT 块 | ◐ | 迁入 state.f90，其余约 40 个表象块全部迁 config |
| venus_params 的 C1..C8/PI | ● | 单位契约保留（L1 常数） |
| venus_params 的 NN_POT_TYPES 索引 | ◐ | 迁入 pes_emt_nn 成员内部 |
| input_parser / param_mapping / presets | ● | config 的文法层原样保留（五张映射表=封闭族注册行） |
| venus_input 读取骨架 | ◐ | 几何/名册/ghost 装配段移交构型装配器（见 △3） |
| △ roster 构型装配器 | △ | 缺。质料：venus_input 几何段 + test_potentials 反写段（BOXLX/SKEW/NATOMB/QZB）+ SURF 表面准备 + FIXROTDATM 对齐 |

### L2 接口 A

| 质料 | 放置 | 说明 |
|---|---|---|
| DVDQ_1 力分派 | ◐ | 其"清零 PDOT→调 DPESHON"逻辑并入 force_facade |
| interface_TEST 的 POTPRE/POT0/DPESHON | ◐ | POTPRE 型初始化钩子成为门面协议；×23.0605、×C1 换算留边界（单位契约履行点） |
| test_potentials: HARMONIC/MORSE/LEPS | ● | pes_reference 成员（体系无关参考势，双实现检验用） |
| emt_nn_pes_venus | ● | pes_emt_nn 成员（自包含+数据外置，已是正确形态） |
| rst_pes + rst_slab_bvk | ● | pes_rst 成员（含 PBC，自己的投影逻辑） |
| （GLO 弹簧力，现嵌 VERLET） | ◐ | pes_glo_spring 成员 |
| （ghost Langevin，现嵌 VERLET） | ◐ | pes_langevin_ghost 成员（随机律） |
| FRICFORCE+NFC | ◐ | pes_ldfa 成员（复活） |
| ENERGY_1 的 NSURF/NGLO 分派与 GLO 势能 | ◐ | 能量侧并入 facade（E = T(P) + ΣV 由门面汇总） |
| △ force_facade 框架 | △ | 缺。质料：DVDQ+interface_TEST 的分派壳；新写内容=成员表数据结构、加性合成、能力级协商 |
| △ elec_adiabatic 退化成员 | △ | 缺（trivial：ds/dt≡0、力规则恒等，但须显式存在——退化点良定义） |

### L2 接口 B′

| 质料 | 放置 | 说明 |
|---|---|---|
| SELECT 抽样器框架 | ◐ | NACTA/NACTB 整数分发改为 dist_* 成员表后即为 sampler |
| THRMAN | ● | dist_boltzmann（双实现之一） |
| INITEBK | ● | dist_ebk（双实现之二） |
| harmonic_sampling | ◐ | dist_surface_oscillator（注册即活，代替删除） |
| ROTEN | ● | dist_rotation（NROT 三律） |
| PROBJ + JMAXCALC | ◐ | dist_j；H2/D2 核自旋表数据化后并入 |
| SURF（入射通道）+ DENS | ◐ | dist_incident（相对运动分布） |
| BAREXC | ● | dist_barrier_excitation |
| THERMO | ◐ | dist_thermostat（Andersen/标温两律；仅装配期用） |
| ROTATE/ROTATENO/ROTATEJKM/ROTATEX/Y/Z | ● | L5 geometry（取向实现，被各 dist 调用） |
| HOMOQP/INITQP/ORThAN/QMMICRO | ◐ | 简正模坐标装配并入 dist_normalmode / dist_qmmicro；NACTA=4 死分支随迁弃 |
| GLOSELECT | ● | dist_glo_target（GLO 靶制备） |
| md_equilibrate（SELECT 内） | ● | dist_thermalize（成员内部调演化服务，2.3 边界 4 允许） |
| NMODE/FMTRX/FGMTRX/ENMODE | ◐ | 装配期分析工具，随 dist_normalmode 放置（GAUSS 死分支弃） |

### L3 封闭族

| 质料 | 放置 | 说明 |
|---|---|---|
| VERLET | ◐ | 剥 GLO/Langevin 内嵌后=封闭传播子 verlet |
| SYMPLE / RADAU | ◐ | RADAU 剥 RADOUT 编排（改回调）后并入 |
| TEST 的通用 COM 判据段 | ◐ | 谓词文法的 and/or/com_distance 原语素材 |
| △ predicate.f90 文法解释器 | △ | 缺。质料：TEST:106-146；新写=文法数据结构（阈值比较/布尔组合/有界窗）+求值器 |
| △ 结局分类配置 | △ | 缺（现为 VENUS.f90 结局块内的写死 if 链）。质料：E1 合并后的分类逻辑数据化 |

### L4 测量 / 后处理

| 质料 | 放置 | 说明 |
|---|---|---|
| GWRITE | ◐ | 剥四处体系知识后=recorder（序列化器）；GWRITE_LEVEL 保留为配置入口 |
| record_all/HIST 记录列 | ◐ | 并入 recorder 的记录表（列选择配置化） |
| FINAL/GFINAL | ◐ | 迁 analysis/；GFINAL:80 EREL 覆盖弃；广义 COM/EREL/OAM 部分即结局分类执行器素材 |
| EIGOUT/DENQ+SUMQ | ● | analysis（态密度/本征输出） |
| HIST 的 plot_* 五函数 | ◐ | analysis 离线绘图（surfcell/Dist_cell 通用实现保留，格子常数改传参） |
| conserv_check / smoke_rst_bvk | ● | analysis 验收工具原样保留 |

### L5 仪器

| 质料 | 放置 | 说明 |
|---|---|---|
| RAND0/RAND1/RANDST | ● | rng.f90 主流（ISEED3 状态） |
| gasdev/GAMA | ◐ | 迁入 rng.f90 并改走主流（单一随机流） |
| EIGN | ● | linalg.f90 |
| GLPAR+PL / rotaline | ● | 数学工具 |
| CPUSEC | ● | 计时（VFDATE 弃——SYSTEM 壳调用） |
| CENMAS/ANGVEL/ROTN | ● | geometry/kinematics |
| ALIGN_ZSTEP/ALIGN_YSTEP | ● | geometry |

### 无位置出局（与 audit 冗余清单一致）

RK4 系数与死入口、ADAMS（已 loud-stop）、MICROCI+VOLPSCONE、LMODE+LMEXCT+EBOND、
ENMODE（MPLOT 死）、POTEN、POTENZ 幽灵打印、THERMBATH、GLOEQU、DISTRIBUTE、PBC、
LINEQ/LUDCMP/LUBKSB/FEHLER、WEBOND/WLBOND/WENMOD、ROTATEX_b/ROTATEZ_b、
grandom/ran1/gasdev_0、GLO.f 滞留副本、VENUS.f90 三死块、venus_data 死块群、
MKL 链接、debug 打印群（SELECT/SURF/ROTen/THERMO/RADAU/FINLNJ 残留）。

---

## 3. 缺口清单（结构需要、库里没有的质料）

按必要性排序：

1. **force_facade 框架**——门面是接口 A 的本体，现库只有分派雏形（DVDQ+POTPRE）。
   新写量：成员注册表（name→procedure 数组）、加性合成循环、能力级协商
   （H 表示 × 传播律 × 积分器兼容三元组检查）。约 200-300 行。
2. **predicate 文法解释器**——终止谓词从硬编码 if 变为数据驱动求值。
   质料：TEST:106-146。新写量：文法解析（输入文件中的谓词表达式）+求值器。
   约 150 行。
3. **roster 构型装配器**——名册/几何/元素表的统一装配点。质料分散在
   venus_input/SURF/test_potentials 三处。新写量：名册数据结构 + 元素表。
   约 150 行。
4. **recorder 观测量配置表**——GWRITE_LEVEL 从"写死的层"变为"用户可选列"。
   新写量：列注册机制。约 50-100 行。
5. **elec_adiabatic 退化成员**——trivial 但必须显式存在（约 20 行），
   为将来 TDHF-FSSH 作为乙类成员回归留位。
6. **结局分类数据**——E1 合并后把 if 链改查表（约 50 行）。

## 4. 结构自检（对照规范）

- 相位限制：driver 三相位调用序列，演化期仅 facade/propagator/recorder/predicate；
- 因子化：任一 pes_* 与 elec_*（现在只有退化成员）自由配对，无捆绑例程；
- 双实现在场：A 轴 HARMONIC+MORSE+LEPS+EMT-NN+RST（≥2），B′ 轴
  THRMAN+INITEBK（+surface_oscillator 注册后 ≥3），构型层面分子/表面两范式
  （twobody vs rst_* cases）；
- 单位契约：换算只在 facade 边界（×23.0605、×C1），state 内部恒 code 单位；
- 单一随机流：rng.f90，ISEED 全流可复现。

## 5. 迁移策略含义

质料映射显示：出局项与 R0 清单重合；◐ 项大多对应 R1-R6 各步；△ 缺口
是纯新增（不动旧码，可与纯化并行开发、旁路验证）。建议实施顺序不变
（R0→R7），但 △1/△2 可提前动工——它们是新结构的地基，且不依赖任何
纯化步骤的完成。
