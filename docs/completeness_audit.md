# VENUS2026 子例程完备性审查报告

- 日期：2026-09-04
- 裁决标准：`docs/design_spec_completeness.md`（物自体/表象二分、四类纯粹子例程、控制核、接口 A/B′、封闭层、三值判定）
- 审查方式：只读静态审查。逐例程通读根源码 10 文件、src_VENUS 61 文件、src_TEST 6 文件、scripts 2 工具，并对全部疑似死入口做了全库调用点核对（grep 交叉验证）。
- 机器可读总表：`docs/completeness_audit.json`（99 条目）。
- 判定用语：**不变类** / **接口或封闭族成员** / **冗余候选或待纯化**。
- 类别用语：演化（物自体→物自体）、测量（物自体→表象）、操作（表象→物自体）、后处理（表象→表象）、控制核、基设（表象层基础设施）、仪器、数学。

---

## 0. 与已知校准点的对照（先行结论）

| 校准点 | 审查结论 |
|---|---|
| GWRITE.f90 历史上硬编码 atom3=C | **未清**。现行形态演化为 **atom1=C**（`GWRITE.f90:59-60` 的 `Q(C)=`/`V(C)=`），并新增两处体系知识：`W(i)<50 amu → 'C '` 否则 `'Au '` 的质量阈元素标注（:110-114）、`Au_nearest`/`Au_xyz` 最近邻诊断（:95-99）。违规从"写死原子 3"变为"写死原子 1 + 质量阈元素名"，性质未变。 |
| RK4 死入口、ADAMS 损坏 → 冗余候选 | **确认**。RK4：`param_mapping.f90` 仍接受 `RK4`（→INTEGRATOR=1, LLL=0），但 `venus_input.f90:363-365` 拒绝 LLL=0 → 读入即停的死入口；`VENUS.f90:86-94` 的 RKUTTA 系数 RAA1..RC2 计算后全库无消费者。ADAMS：INTEGRATOR=0 在 `venus_input.f90:339-343` 显式 loud-stop——已从"损坏"升级为"显式拒绝"，出局结论成立。 |
| TEST.f 终止判据 → 应为封闭谓词文法成员 | **确认**。`TEST.f:88-94` 有**活的**体系特异分支（NSURF>0 时以写死的 atom1 高度 Q(3) 判散射、NC>50 魔数、P(3)>0 外向速度——C/Au(111) 判据）；:54-78（F+H2O）与 :149-192（CH4/Ni(111)）为两段注释掉的死体系知识；:106-146 的通用 COM 距离+退行速度判据是封闭文法成员的正确原型。 |
| 编译期无 PES 变体（TEST_PES 运行时选择）→ 已符合接口 A | **确认**。`interface_TEST.f90` 头注释明示单构建多势；HARMONIC/MORSE/LEPS/EMT-NN/RST 五成员经 `TEST_PES` 关键字运行时注册，无编译期变体增 proliferation。 |

四项校准点结论与审查一致，无需要说明的分歧。

---

## 1. 总表（全例程覆盖）

判定列缩写：不=不变类；族=接口/封闭族成员；冗=冗余候选/待纯化。违规点给文件:行号。

### 1.1 根源码（10 文件）

| 例程 | 文件 | 四类归属 | 判定 | 违规点 / 备注 |
|---|---|---|---|---|
| VENUS (program) | VENUS.f90 | 控制核 | 冗 | :86-94 RKUTTA 系数死计算；:182-232 NCHKP 死分支（无输入来源）；:184,231 RANDST 双调用（F15）；:286-299 直方图死块（idplt=0 写死）+ Qr1=2.94998/Qr2=2.55476 位点常数 + SYSTEM 壳调用；:322-335 GLO 布局分支；:340 表面对齐调用；:399-400,454-458 BOUNCE/z_pos_min 观测量积累写死原子 3；:462-464 空检查语句；:544-548,631-635,671-677 AE* 观测记录内嵌；:506-576 与 :581-662 结局块整段双份复制 |
| venus_params | venus_params.f90 | 基设（单位契约/常数） | 冗 | :88-98 NN_POT_TYPES=7 及 E0/1D/3P/离子态索引——C+Au 系 NN 势知识写死在常数模块；C1..C8/PI 为正当单位契约 |
| venus_data + allocate_venus_data | venus_data.f90 | 全局状态 | 冗 | 混居：物自体（Q/P/PDOT/QDOT/Q_old）与约 40 个表象块同模块、无注释区分；:585 auro=2.94998（Au 晶格常数）；:382-384 HIST_ATOM_NAD 无写入者；:389-391,417-418 MORSEB/STRETB 死数据；CUBEB/COORS/GDD/GAUSS2/U0-U1-V01/RKUTTA 等死块 |
| read_venus_input | venus_input.f90 | 操作（表象读取装配） | 冗 | :616-621 ghost Langevin 装配写死 W(2) 与 413.4 换算；:703-704 NATOMB 域写死 4/10；:997-998 BOXLX=32.5/BOXLY=30.7 默认几何；:752 W=1D30 刚性质量；:692 MPLOT=0 写死（判死整条模式布居通路）；:392 NCVIB 死关键字。读取文法本身合格 |
| input_parser | input_parser.f90 | 基设（表象文法） | 不 | 纯文法（关键字/续行/大小写/构造子），旧格式 loud-stop |
| param_mapping | param_mapping.f90 | 基设（枚举注册） | 不 | 五张字符串→整数映射表=封闭族注册行；已移除项 loud-stop |
| presets | presets.f90 | 基设（固化预设） | 不 | 五预设 set_default 数据化——"体系规定从代码搬进数据"的正确方向 |
| sample_surface_oscillator | harmonic_sampling.f90 | 操作（B′ 成员素材） | 冗 | 全库无调用点（编译但死）；应接入接口 B′ 或删除 |
| ENERGY_1 | ENERGY.f | 测量（T/V/H 泛函） | 冗 | :16-49 NSURF/NGLO 分派；:40-47 GLO 弹簧势能内嵌（写死 NATOMA+1/+2）——接口 A 内容越界；POT_SH 注释死调用 |
| DVDQ_1 | DVDQ.f | 演化（力评估分派） | 不 | 清零 PDOT→调 DPESHON；接口 A 的唯一消费点，PES 当黑盒 |

### 1.2 src_TEST（6 文件）

| 例程 | 文件 | 四类归属 | 判定 | 违规点 / 备注 |
|---|---|---|---|---|
| POTPRE | interface_TEST.f90 | 演化（接口 A 注册） | 不 | 运行时 TEST_PES 单构建多势——符合接口 A |
| POT0 | interface_TEST.f90 | 演化（接口 A 边界） | 不 | eV→kcal/mol（×23.0605）换算集中于边界（单位契约履行）；副作用写 E0 供诊断（轻微） |
| DPESHON | interface_TEST.f90 | 演化（接口 A 边界） | 不 | −grad×23.0605×C1 累加 PDOT；签名承认加性贡献 |
| FIXROTDATM | interface_TEST.f90 | 操作（取向对齐） | 族 | :58-72 写死原子 1-2 轴向 z；参数化后为封闭成员 |
| ALIGN_ZSTEP / ALIGN_YSTEP | interface_TEST.f90 | 操作（运动学） | 不 | 刚体旋转步进 |
| gasdev | interface_TEST.f90 | 仪器（正态随机律） | 冗 | :112-136 核心仪器住在 PES 接口文件（层次倒置，靠 Makefile 顺序敏感存活）；走 RANDOM_NUMBER 流（与 RAND0 流并存，ISEED 不可复现）。活代码 |
| test_potentials（模块全族） | test_potentials.f90 | 演化（接口 A 成员族） | 族 | :152-158,287-293 成员反写表象（BOXLX/BOXLY/SKEW）；:163,260 成员写名册 NATOMB；:210-213,274-281 成员写装配几何 QZB。HARMONIC/MORSE/LEPS 为体系无关参考成员；EMT-NN/RST 携带 C/Au(111) 属成员内部合法 |
| emt_nn_pes_venus | emt_nn_pes.f90 | 演化（接口 A 成员） | 族 | 自包含引擎（不 use venus_data），数据外置 data/emt_nn/ 且路径关键字可覆盖——成员携带外部数据的正确形态 |
| rst_pes | rst_pes.f90 | 演化（接口 A 成员辅助） | 族 | 零依赖模块；FCC(111) slab 生成/权重 IO |
| rst_slab_bvk | rst_slab_bvk.f90 | 演化（接口 A 成员） | 族 | BVK 弹性 slab + RST GS；rst_total_vg 交付 V 与梯度；含 PBC |

### 1.3 src_VENUS（61 文件）

| 例程 | 文件 | 四类归属 | 判定 | 违规点 / 备注 |
|---|---|---|---|---|
| ANGVEL | ANGVEL.f | 操作（运动学） | 不 | 角速度扣除 |
| rotaline | ArbitraryAxisRotation.f | 数学 | 不 | 任意轴旋转 |
| BAREXC | BAREXC.f90 | 操作（势垒激发成员） | 族 | NSELT=3 反应坐标动量；IRC 方向经配置 |
| CENMAS | CENMAS.f90 | 操作（运动学） | 不 | COM 分离；QQ/PP 工作区 |
| CPUSEC | CPUSEC.f | 仪器 | 不 | 计时 |
| DENQ / SUMQ | DENQ.f | 后处理（态密度） | 不 | Beyer-Swinehart；DEL=1 cm⁻¹ 写死（可配置化） |
| EBOND | EBOND.f | 测量（局域模能量） | 冗 | 调用者仅 MPLOT 死分支与 LMEXCT（NACTA=4 已拒）→死入口；依赖无写入者的 MORSEB |
| EIGN | EIGN.f | 数学 | 不 | Givens-Householder 对角化 |
| EIGOUT | EIGOUT.f | 测量（输出） | 不 | IP=77 与 fort.77 复用（轻微） |
| ENMODE | ENMODE.f | 测量（简正模能量） | 冗 | 唯一调用在 MPLOT==1 块，MPLOT 被写死 0→死入口 |
| FGMTRX | FGMTRX.f | 操作（简正模分析） | 不 | :26-40 Gaussian Hessian 死分支（GAUHES 无写入者） |
| FINAL | FINAL.f | 测量（终态性质） | 冗 | :131-135 fort.777 写 U0/U1/V01（全库无写入者，打印未定义值）；:616-620 NZDOWN 分支写死原子 3（P(7:9)/W(3)）。广义 COM/EREL/OAM/ANG1-20 部分合格 |
| FINLNJ | FINLNJ.f | 测量（RKR 反演 n,J） | 不 | :60-63 debug 赋值残留；经 ENERGY_1 黑盒取势 |
| FMTRX | FMTRX.f | 操作（力常数差分） | 不 | :25-60 Gaussian 分支死（GAUHES 无写入者、nmo 未初始化） |
| GAMA | GAMA.f | 仪器（Gamma 分布） | 不 | :9 已改走 RANDOM_NUMBER 内在流——与 RAND0 流分裂 |
| grandom | GASDEV.f | 仪器 | 冗 | 无调用者（与 interface_TEST 的 gasdev 重复实现） |
| ran1 | GASDEV.f | 仪器 | 冗 | 无调用者 |
| gasdev_0 | GASDEV.f | 仪器 | 冗 | :48-71 函数名与返回名不符——损坏死代码 |
| GLOINIT | GLO.f90 | 操作（GLO 单位换算） | 族 | 413.4 换算 |
| GLOSELECT | GLO.f90 | 操作（GLO 靶制备成员） | 族 | 接口 B′ 侧成员 |
| GLOEQU | GLO.f90 | 操作 | 冗 | 无调用者（被 md_equilibrate 取代） |
| DISTRIBUTE | GLO.f90 | 后处理 | 冗 | 无调用者 |
| （整文件） | GLO.f | — | 冗 | 未列入 Makefile SRC 的 GLO.f90 近重复文件，滞留源树 |
| GLPAR (+PL) | GLPAR.f | 数学 | 不 | Gauss-Legendre 求积 |
| GWRITE | GWRITE.f90 | 测量（轨迹序列化） | 冗 | 见 §0 校准点：:59-60 atom1=C；:95-99 Au_nearest；:110-114 质量阈元素名；:78-79 QZB 极值+6.0 Å PBC 盒；:43 内联 kB；GWRITE_LEVEL 分层是配置（好）但各层观测量集合写死——不完备测量子例程，应纯化 |
| record_all | HIST.f90 | 测量（终态记录） | 冗 | :5-6 自声明"仅限 C/Au"；:19 auro+PI/3 六方格子写死；:18-22 记录列含 NAMD 死量 AEphS |
| plot_Erel_hist 等 5 个 | HIST.f90 | 后处理（直方图） | 冗 | :165,340-343 FCC(111) 位点分类写死；唯一入口是 VENUS.f90 idplt=0 死开关→现为死代码；区间/箱数为参数（好） |
| surfcell / Dist_cell | HIST.f90 | 后处理（几何） | 冗 | 实现通用（le/skew 传参），调用方写死 auro 与 PI/3 |
| HOMOQP | HOMOQP.f90 | 操作（双原子装配） | 不 | 取向分发经 NTHTA 配置 |
| INITEBK | INITEBK.f90 | 操作（EBK 成员） | 族 | 接口 B′ 双实现素材（Boltzmann+EBK）之一 |
| INITQP | INITQP.f90 | 操作（简正模装配） | 不 | :46,227-233 NACTA=4 死分支；能量标度重采样循环 |
| JMAXCALC | JMAXCALC.f | 操作（J 分布上界） | 冗 | :35-48 按质量 1.008/2.014 写死 H2/D2 核自旋统计——体系知识应退为输入数据 |
| LMEXCT | LMEXCT.f | 操作（局域模） | 冗 | F20：输入即拒；依赖无写入者的 MORSEB |
| LMODE | LMODE.f | 操作（局域模） | 冗 | 同上，两处调用点均死 |
| MICROCI | MICROCI.f | 操作（CI 微正则） | 冗 | F23：死循环，输入即拒 |
| VOLPSCONE | VOLPSCONE.f | 操作 | 冗 | 仅 MICROCI 调用→传递性死代码 |
| NMODE | NMODE.f90 | 操作（分析驱动） | 不 | 坐标/梯度保存-恢复包裹 FMTRX |
| ORTHAN | ORTHAN.f90 | 操作（正交象限成员） | 族 | 装配期用 TEST 谓词作几何界（初始化性测量，相位合法） |
| POTEN | POTEN.f | 测量（取向扫描） | 冗 | 无调用者（旧任务删除后遗留） |
| POTENZ | POTENZ.f | 操作（能量参考） | 不 | :76-82 fort.666 写无写入者的 U0/U1/V01 |
| PRINFO | PRINFO.f | 基设（回显） | 不 | 约 70 行注释掉的旧读取残留 |
| PROBJ | PROBJ.f90 | 操作（J 抽样） | 不 | 核自旋权重经 WGT（由 JMAXCALC 填） |
| RADAU | RADAU.f | 演化（RADAU15 封闭传播子） | 族 | :438-498 RADOUT 内嵌编排（调 GWRITE/TEST/FINAL/GFINAL+维护 NAST 状态机——编排越位）；:112,380,428 'CALLING RADOUT' 调试打印。算法本体不变 |
| RADOUT | RADAU.f | 测量+编排 | 冗 | 应改为向控制核回调的钩子 |
| RAND0 / RAND1 / RANDST | RAND0.f / RAND1.f / RANDST.f | 仪器（Knuth 流） | 不 | 可复现主流；F15 双调用在调用方 VENUS.f90 |
| ROTATE / ROTATENO | ROTATE.f90 | 操作（取向） | 不 / 族 | ROTATENO 为场致取向分布成员（注释残留 NO/CO 语境，实现通用） |
| ROTATEX / ROTATEY / ROTATEZ | 同名 .f | 操作（运动学） | 不 | |
| ROTATEX_b / ROTATEZ_b | ROTATEX.f / ROTATEZ.f | 操作（运动学） | 冗 | 无调用者（2024.11 定角变体） |
| ROTATEJKM / ROTATEJM | 同名 .f | 操作（(J,M) 取向） | 不 | ROTATEJKM 尾部激光偏振注释块（死知识备注） |
| ROTEN | ROTEN.f90 | 操作（转动抽样成员） | 族 | :101-105 调试打印+一致性 stop 残留；NROT=0/1/2 三律 |
| ROTN | ROTN.f90 | 测量（角动量/转动能） | 不 | |
| SELECT | SELECT.f90 | 操作（抽样器框架+制备律分发） | 不 | :156-160,362-400,417,1321-1326 调试打印残留；:1139-1143 MD 平衡化 +100 Å 魔数。NACTA/NACTB 整数分发=封闭族注册行；md_equilibrate 装配期调演化服务（2.3 边界 4 允许） |
| SURF | SURF.f | 操作（入射通道成员） | 族 | :191 write(88) 调试；:48-50 SKEW=π/4 魔数；NSURF=1/2 双实现检验在场 |
| PBC | SURF.f | 演化（PBC 折影） | 冗 | 唯一调用（ENERGY.f:17）被注释→死入口；PBC 实际由 rst_slab_bvk/GWRITE 各自实现 |
| DENS | SURF.f | 操作（NREL=2 速度分布） | 族 | 活路径 |
| SYMPLE | SYMPLE.f | 演化（辛积分器封闭传播子） | 族 | :94-96 CALTYP 分支残迹（绝热版恒真）；Schlier 系数 |
| TEST | TEST.f | 测量（终止谓词） | 冗 | 见 §0 校准点：:88-94 活的 C/Au 判据；:54-78/:149-192 两段注释体系知识；:106-146 通用 COM 判据为文法成员原型 |
| THERMO | THERMO.f | 操作（恒温成员） | 族 | :56,115-122 调试输出残留；info=0/1 两律（标度/Andersen） |
| THERMBATH | THERMBATH.f | 操作（恒温浴） | 冗 | 无调用者；NTHERMB 写死 0、NTHMID 无写入者 |
| THRMAN | THRMAN.f | 操作（Boltzmann 量子数成员） | 族 | 接口 B′ 双实现素材之一 |
| VFDATE | VFDATE.f | 仪器 | 冗 | 无调用者；经 SYSTEM(date) 壳调用 |
| VERLET | VERLET.f90 | 演化（VV/Beeman 封闭传播子） | 冗 | :44-50,71-79 ghost Langevin 写死原子 2（k=3+i）；:137-148 GLO 弹簧+随机+摩擦内嵌 Beeman 分支——接口 A 力贡献越界进封闭层；:153-161 同；:34-36,60-62,123-125 CALTYP 残迹。核传播本体不变 |
| FRICFORCE | VERLET.f90 | 演化（LDFA 摩擦+噪声候选成员） | 冗 | :130-133 调用点被注释→死入口；NFC=1/2 随机律成员现成素材 |
| DIAG2 | VERLET.f90 | 数学 | 冗 | 仅被死入口 FRICFORCE 调用 |
| LINEQ / LUDCMP / LUBKSB / FEHLER | VERLET.f90 | 数学 | 冗 | 全库无调用者；FEHLER 为空壳 |
| WEBOND / WLBOND / WENMOD | 同名 .f | 后处理（布居写出） | 冗 | 无调用者 |
| QMMICRO | QMMICRO.f | 操作（量子微正则成员） | 族 | NACTA/B=8 活路径 |

### 1.4 scripts（独立工具，不在主构建）

| 例程 | 文件 | 归属 | 判定 | 备注 |
|---|---|---|---|---|
| conserv_check | scripts/conserv_check.f90 | 后处理 | 不 | 离线逐帧能量重构；暴露主程序演化环内不刷新 T/V/H 的接口缺口 |
| smoke_rst_bvk | scripts/smoke_rst_bvk.f | 后处理 | 不 | RST/BVK 成员的解析/数值梯度验收——接口 A 双实现验证程序 |

**统计**：99 条目——不变类 39、接口/封闭族成员 19、冗余候选/待纯化 41。

---

## 2. 架构级发现

### A1 全局状态混居（venus_data）

物自体区（Q、Q_old、P、PDOT、QDOT，:75-84）与约 40 个表象块（打印旗标、名册/几何 QZA/QZB/LA/LB、GLO 参数、恒温参数、GLOP ghost、统计数组 HIST_ATOM 等）同模块存放，无任何注释区分两类——违反"混居须报"判据。另有：`auro=2.94998`（Au 最近邻距离，:585，被 HIST/record_all 活用）；十余个无写入者或无消费者的死块（MORSEB、STRETB、CUBEB、COORS、GAUSS2/GAUHES、GDD、HIST_ATOM_NAD、U0/U1/V01、RKUTTA 块）。

### A2 双随机数流并存（可复现性缺口）

- 可复现主流：RAND0/RAND1/RANDST（Knuth 64 位乘同余+游混，ISEED3 状态）——SELECT/SURF/PROBJ/INITQP 等抽样走此流。
- 编译器内在流：`gasdev`（interface_TEST.f90:112）、`GAMA`（GAMA.f:9）、`grandom` 走 `RANDOM_NUMBER`——VERLET 的 Langevin 噪声、GLOSELECT、Maxwell-Boltzmann 动量、THERMO 的 Andersen 重抽样均在此流上。
- 后果：ISEED 无法控制第二条流，同一 ISEED 下含噪声/恒温的轨迹不可复现；且 **gasdev 这个核心仪器物理上住在 src_TEST/interface_TEST.f90（PES 接口文件）**，src_VENUS 核心反向依赖 src_TEST 提供仪器（层次倒置，仅靠 Makefile 顺序敏感成立）。

### A3 接口 A 符合性

- **符合**：无编译期 PES 变体（运行时 TEST_PES）；单位换算集中于 interface_TEST 边界（×23.0605、×C1，符合单位契约"换算集中于接口边界"）；五成员一行 case 注册；EMT-NN/RST 成员自包含、数据外置、路径可覆盖。
- **违反**：
  1. GLO 弹簧力与随机/摩擦力内嵌 `VERLET.f90:137-148`，GLO 弹簧势能内嵌 `ENERGY.f:40-47`，ghost Langevin 内嵌 `VERLET.f90:44-50,71-79,153-161`——按规范 3.2"噪声律归力贡献族（要求接口 A 签名承认随机律）"，这些应是接口 A 成员，且原子索引（atom2、NATOMA+1/+2）写死；
  2. `test_potentials.f90` 成员反向写表象与名册（BOXLX/BOXLY/SKEW 覆盖输入晶胞、NATOMB(1)、QZB slab 几何）——成员越出 Q→V/F 签名，几何/名册应经输入文件或构型装配器；
  3. FRICFORCE（LDFA 摩擦+噪声）是现成的随机律成员素材，但调用点被注释成死入口。

### A4 相位限制检验（演化期无操作类）

**满足**。演化环（VENUS.f90 INTEGRATION_LOOP 与积分器内部）只调用：DVDQ_1（演化）、TEST/GWRITE/FINAL/GFINAL/ENERGY_1/record_all（测量，含终止谓词）、积分器。操作类（SELECT/SURF/GLOSELECT/ROTATE* 等）全部位于装配相位。装配期 `SELECT→NMODE/DVDQ/ENERGY`、`md_equilibrate→VERLET/SYMPLE/THERMO` 属"制备成员内部调用演化服务"，规范 2.3 边界 4 明示允许。唯一越位是 `RADOUT`（RADAU.f:438-498）在传播子内部复制控制核编排（调 FINAL/GFINAL、维护 NAST 状态机、自行 return）——属**编排越位**而非相位违规。

### A5 编译期耦合

- Makefile 单次 ifx 调用按序编译（顺序敏感）；gasdev 的层次倒置（见 A2）依赖此顺序。
- `src_VENUS/GLO.f` 未列入 SRC——不编译的近重复文件滞留源树。
- MKL 三库（Windows 显式路径/Linux -qmkl）被链接，但全部 Fortran 源无任何 MKL/LAPACK 调用——**死构建依赖**（对角化全走自带 EIGN/DIAG2）。
- scripts 两工具游离于构建系统外（自带编译命令行）。

### A6 控制核不完备（VENUS.f90）

控制核混入：体系分支（GLO 布局 :322-335、表面对齐 :340、THTA 特例 :344-346）；观测量选择与积累（BOUNCE/z_pos_min 写死原子 3、AE* 数组记录）；结局处理+统计块整段双份复制（:506-576 与 :581-662，含 EROTA 平均与反应计数）；三块死代码（NCHKP 分支、RKUTTA 系数、idplt 直方图块）。按分解定理，控制核应只剩"装配→演化→统计"的调用序列。

### A7 测量层不完备（观测量选择未退为配置）

- GWRITE：atom1=C、质量阈元素名、nearest-Au、+6.0 Å PBC 盒（校准点详情见 §0）。
- TEST：活的 C/Au 高度判据（体系分支）。
- record_all/HIST：Au(111) 六方格子（auro、PI/3、top/brg/hcp/fcc）。
- **GFINAL.f:80**：`EREL=NORM2(P(7:9))**2/2/W(3)/C1`——用写死原子 3 的动能覆盖 FINAL 算出的广义相对能量（注释自认"Erel for 1D model"），是结局分类数据被体系知识污染的最重一处。
- FINAL.f:131-135 与 POTENZ.f:76-82：fort.777/666 写 U0/U1/V01——全库无写入者的 NAMD 遗留，打印未定义值。
- GFINAL.f:66-69：current_U0 无定义（隐式局部未初始化即打印）。

### A8 冗余总量（约 26 处）

死文件：GLO.f。死入口/死代码：harmonic_sampling（无调用）、VFDATE、POTEN、EBOND、ENMODE（MPLOT 死）、WEBOND、WLBOND、WENMOD、MICROCI+VOLPSCONE（F23）、LMODE+LMEXCT（F20）、GLOEQU、DISTRIBUTE、THERMBATH、PBC、FRICFORCE+DIAG2、LINEQ+LUDCMP+LUBKSB、FEHLER、ROTATEX_b、ROTATEZ_b、grandom、ran1、gasdev_0、RKUTTA 系数、NCHKP 分支、idplt 直方图块、NCVIB 死关键字、MKL 链接、venus_data 死块群。F20/F23/F14 等已在输入层 loud-stop（防误用），但实现体仍滞留编译。

---

## 3. 纯化建议清单（现状 → 目标位置）

1. **TEST.f:88-94 表面 C 原子判据** → 封闭谓词文法成员：`z_height(原子, 阈值, 外向速度)` 与 `com_distance(片段对, 阈值, 退行速度)` 的布尔组合；体系判据经输入文件装配（:106-146 通用判据已是原型）；两段注释体系知识删除。
2. **GWRITE 观测量写死（:59-60,95-99,110-114）** → 序列化器只写 Q/P/F 数组与单位标签；元素名由输入名册提供；C 原子/nearest-Au 诊断退为后处理配置。
3. **GFINAL.f:80 EREL 覆盖（写死原子 3）** → 删除覆盖、采用 FINAL 广义 COM 相对能量；确需自由度子集时由结局分类配置指定。
4. **VERLET.f90:44-50,137-161 与 ENERGY.f:40-47 的 GLO/ghost 力** → 接口 A 力贡献族成员（GLO-spring 成员、Langevin-ghost 成员，签名承认随机律与摩擦）；积分器/能量例程回到纯核传播/纯泛函。FRICFORCE 可顺带复活为 LDFA 成员。
5. **venus_input.f90:616-621 ghost 参数装配（W(2) 写死）** → 移入 GLO/Langevin 成员自己的初始化钩子（POTPRE 型）。
6. **venus_data.f90:585 auro 与 HIST 六方格子（PI/3、top/brg/hcp/fcc）** → 表象输入参数/表面格子类型配置（slab 名册数据）。
7. **venus_params.f90:88-98 NN 势类型索引** → NN PES 成员文件内部常量。
8. **JMAXCALC.f:35-48 H2/D2 核自旋统计（质量判别）** → 输入数据（核自旋权重表）。
9. **venus_input.f90:997-998 BOXLX/BOXLY 默认几何与 :703-704 NATOMB 写死域** → 必填输入或由 slab 名册推导；校验泛化为名册一致性检查。
10. **VENUS.f90 观测量积累与双份结局块** → 测量记录器模块+单一结局路径；BOUNCE/z_pos_min 的原子索引退为配置；死块（NCHKP/RKUTTA/idplt）删除。
11. **gasdev 位置与随机流** → 移至 src_VENUS 独立仪器模块；全库统一单一随机流（或全走 RAND0，或全走内在 RNG），恢复 ISEED 全流可复现；GAMA 的 IDUM 形参恢复实效。
12. **test_potentials.f90:152-163,260-293 成员写 BOXLX/NATOMB/QZB** → 构型装配器（名册/几何经输入文件交付），PES 成员只保留 Q→V/F 映射与数据文件读取。
13. **RADAU.f RADOUT 内嵌编排** → 传播子向控制核回调（输出钩子+终止谓词求值点），编排归还控制核；三处 'CALLING RADOUT' 调试打印删除。
14. **约 26 处死代码（A8 清单）** → 删除；harmonic_sampling 若保留则作为接口 B′ 成员正式注册（表面振子 Boltzmann 律），与 THRMAN/INITEBK 构成三实现。
15. **venus_data 分区** → 至少以注释划分物自体区/表象区/仪器区，死块随删除清理；MKL 链接去除（或改用其 LAPACK 替换 EIGN，二选一，消除"链接而无消费者"状态）。

---

## 4. 结语

本库在**接口 A 的注册形态**（运行时 TEST_PES、单位换算集中于边界、成员自包含+数据外置）与**输入文法**（keyword=value+构造子映射+预设数据化）上已达到规范要求的封闭形态；接口 B′ 的双实现素材（THRMAN=Boltzmann、INITEBK=EBK）与靶制备成员（GLOSELECT、md_equilibrate、SURF 入射通道）齐备。主要欠缺集中在三处：(1) **噪声律未入接口 A**（GLO/ghost/LDFA 力写死在积分器与能量例程内）；(2) **测量层观测量选择未退为配置**（GWRITE/TEST/HIST/GFINAL 四处体系知识，含 GFINAL.f:80 对广义 EREL 的写死覆盖）；(3) **控制核与全局状态的混居与死代码**（约 26 处冗余、双随机流、双份结局块）。按第 3 节清单纯化后，可变代码将集中于接口 A/B′ 成员与封闭族注册行，与规范 2.3 的可变代码类分布结论对齐。
