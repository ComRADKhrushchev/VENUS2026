# V 计划：全链验证跑（数值后验 + case 集建立）

> 依据：docs/final_completion_plan.md §2（V/S 行：new case set over the full
> chain；旧程序只作 ad-hoc 同量参考非门）；全局约束 8（判据先行）；C5-R1
> （B 层/链级数值后验——法定）；C2/C7/C8 后验清单（→V 各条）；GWRITE_LEVEL
> 用户裁决仍悬（随 V 收口再催一次）。
> 用户 2026-09-15 放行执行。
状态：V 已完成（2026-09-16）
日期：2026-09-15（确认）/ 2026-09-16（收口）

**目标：** 完工主体放上**真力链**跑数值后验：C5 期机制先行的 B 层三员与
频率依赖成员现在有真力可验；建立新树 **case 集**（cases/ 惯例沿旧库：
README/input/expected）；守恒判据家族全链全跑。

**范围裁定（拟）：**
- S1 **case 集建立**（cases/ 新目录，沿旧库三件套惯例：README.md（验证
  目的/原理/预期与结果）+ input_qct.txt + system/ 文件夹 + results 记录）：
  - case_harmonic_scatter：分子两片段 + 谐振容器 + verlet——能量守恒/
    统计面基准（与 check_driver 判据呼应但以 case 形落档）；
  - case_radau_conservation：同体系 INTEGRATOR=RADAU——自适应守恒；
  - case_surface_window：表面体系（POSCAR + A_LAT/BOXPAD + 分子片段）
    ——表面窗口键集端到端（气相判据 gfortran 树承载）；
  - case_b_layer：B 层三员数值后验的 case 化（thermostat/thermalize/
    glo_target 各一小节，判据见 S2）；
  - case_chain_frequencies：Hessian→频谱→boltzmann/normalmode 真值
    （S3 的 case 化）。
- S2 **B 层数值后验**（C5-R1 法定）：真力链下——thermalize（Maxwell 种子
  + MD 平衡：平衡后动能趋近目标温度的统计带——容差与样本量写明）；
  thermostat（速度重置的平衡温度维持）；glo_target（FDT 装载的数值一致
  性：阻尼-噪声关系在链上的实现核对——解析后验而非长程统计）。
  - 450× 备忘（V.1 移交，task-1 报告）：radau 下 thermalize 平衡覆盖约
    450× verlet 的物理时长（自适应步长序列增长所致）——成员语义即"走
    注册传播器"，属行为差异非缺陷；B 层平衡统计以 verlet 为标定基
    （case 用 verlet，与标定一致），radau 复用需按其步长序列重标定。
- S3 **链级频率真值**（C5-R1 法定 + dens/eig 前置演示）：reference 容器
  真力 → Hessian → eig_sym → 频谱（谐振体系的解析频率对照——位级/
  容差写明）→ boltzmann/normalmode 以真频率采样（几何分布/微正则锚复验
  ——C5 判据的合成频率版换成真频率版）。
- S4 **守恒家族全跑 + V 报告**：全 case 能量/动量守恒判据汇总表；ifx
  旧链位级复验（可选——环境允许则做）；thta 退化文件用户说明（C8 移交）
  入 case README；GWRITE_LEVEL 裁决催办注记。

**设计裁决（拟）：**
- D1 判据先行：每 case 的预期锚先写进 README（解析/统计带），跑后
  results 对照——不绿即 BLOCKED 修主体（V 是验证不是修补场：发现的
  主体缺陷回相应模块的 fix 流程，记 V 发现清单）。
- D2 E12：case 材料是体系定义（用户侧）——体系知识进 case 文件夹
  合法（这正是 container 材料的家）。
- D3 R5 沿用（case README 是文档非源码——旧库引用亦不禁，但保持简洁）。

**步骤：**
- [x] **V.1 case_harmonic_scatter + case_radau_conservation**：两基准 case
  建立 + 跑 + results 入档（提交 560b2ae，全锚绿）；
- [x] **V.2 case_surface_window**：表面窗口端到端 case（提交 c8f0ae2，
  全锚绿；ω 换算方向勘误随 V.5 批入档——见验收证据）；
- [x] **V.3 B 层数值后验**（case_b_layer）：三员真力链判据（提交 18dddf9，
  19 锚腿 18 绿 1 红腿定性为试跑标定失误非主体缺陷——见发现清单 V-1）；
- [x] **V.4 链级频率真值**（case_chain_frequencies）：S3 全链（提交 c2c5025，
  全锚绿；通道停放定性——见发现清单 V-5）；
- [x] **V.5 V 报告 + 守恒家族汇总**：S4；提交（本件——验收证据与 →S 移交
  区如下）。

**验收：**
1. 五 case 建立，README 预期锚先行，results 全绿入档；
2. B 层三员与频率链数值后验落地（C5-R1 关闭）；
3. 守恒家族汇总表 + V 发现清单（主体缺陷若有——回修流程）；
4. 双编译器绿（case 跑 gfortran 树——ifx 目录门环境备忘沿档）。

**验收证据（V.5 收口，2026-09-16）：**

五 case 收口（b_layer 一条红腿经 AR(1) 链代数定性为标定失误非主体缺陷
——V-1；其余锚全绿；判据先行：README §3 锚先写后跑、§4 对照入档）+
B 层三员数值后验关闭 + 链级频率真值关闭（C5-R1 两半齐）——**主体缺陷 0**。

守恒家族汇总表（能量 = 帧级守恒锚 vs 实测；动量 = 既有零总动量类/
动量面判据援引；全 case 确定性判据双跑逐字节一致）：

| case（积分器·体系） | 能量守恒锚（帧级 kcal/mol） | 实测·打印行 max\|e_fin−e0\| | 实测·帧级最坏 | 动量面（既有判据援引） |
|---|---|---|---|---|
| harmonic_scatter（verlet·气相双片段） | ≤ 2.0e-1 | 1.2648E-03 | 1.2648e-03 | —（井力下总动量不守恒，无此判据面） |
| radau_conservation（radau·同体系） | ≤ 1.0e-6 | 5.6843E-14 | 1.1369e-13 | —（同上；完成/分类面积分器无关 R2） |
| surface_window（verlet·表面单片段） | ≤ 2.0e-1 | 7.4535E-04 | 1.9419e-03 | thta 真值 = 末帧动量独立重算逐轨迹一致（A4） |
| b_layer（verlet ×3 跑） | ≤ 5.0e-3 | 2.9746E-06 / 1.7462E-05 / 0.0000E+00 | 1.617e-04 / 1.7462e-05 / 逐位 0 | thermalize 精确 COM 投影（χ²₉ 恰 12−3 自由度，T1 带心 0.954 证实）+ thermostat 首帧动量独立重算最坏 4.0e-16（H2） |
| chain_frequencies（采样链·不传播） | N1 能量重展开 ≤ 1e-10 相对 | e_tot = 0.5230 精确重展开 | 同左 | N1 \|Σp\| ≤ 1e-9·√(2·M·e_tot) 全过（零总动量恒等式） |

（glo 静态跑全帧 t_kin = v_pot = e_tot = 逐位 0——"ghost 不触真态"的位级
守恒证词；全 case 同输入双跑逐字节一致 ×8 跑在案。）

**ifx 旧链位级复验（S4 可选项——环境允许，已做）**：按 check_facade 头注
STEP 1 同法以 ifx 编旧链 8 件（`-r8 -double-size=64 -i8 -O2 -w`——旧库
自家 Makefile 同旗标）入 `build/old_xcheck_ifx/`；check_facade.f90 以 ifx
编译（新旧两 .mod 树并读）并链现有 build/ifx 树全部对象（除 driver.o）+
旧链 8 件 → `build/check_facade_ifx.exe`（产物全落 build/ 不入源树）。
**(s) 旧-新位级力对照 ifx 版全绿**：harm/morse/leps 三势 × 三几何
bitwise = T、max|dF| = max|dE| = 0.0000E+00——四配置补全（{gfortran 旧,
gfortran 新} C2 在案 + {ifx 旧, ifx 新} 本次），单编译器内位级等价两
编译器各证一次（跨编译器位级等价不承诺，C2 既定）。全电池 31 PASS /
2 FAIL，两红腿均非 (s) 且非主体：(g) holder read-back 的 h_ij 严格 ==
比较在 ifx 下差 1 ulp——检查侧两路计算**构造性**不等价（fill 的
u = (qe_ctr+ug)−qe_ctr 对 ug 六分量中两分量差 1 ulp，Python 复算证实；
gfortran 的 dot_product 恰好吸收之），-O0 复跑仍在（非优化序）；容差腿
（1e-12）与全部结构/生命周期事实（(e)/(p)/(h)/(i)/(j)、(g) 孪生位级、
(r) 1x1 位级）ifx 全绿。(r) 2x2 肢 = C2 停放区既有永久性状（eig_sym
折叠算术序），ifx 下力偏差 1.2854e-15（刚过 1e-15 门槛）、e_tot 逐位 0
——同族非缺陷。判据行摘录入档：docs/plans/results_check_facade_ifx.txt
（两条 FAIL + (s) 位级块及孪生/同族行 + 汇总行，字面未改；build/ 原件
不入库）。

**终跑核验（抽验两 case）**：case_harmonic_scatter 与 case_surface_window
各复跑一次（gfortran 树 venus.e，2026-09-16）——stdout 五证据行
（assembly/statistics/outcome/drift/histograms）与入库 results.txt 逐字节
一致（timing 行除外），退出码 0、stderr 空；全部 case 的确定性已由各自
双跑判据证（V.1-V.4 在案），抽验确认入库状态可复现。

**V 发现清单（终局）——主体缺陷：0**：

| # | 类别 | 内容 | 处置 |
|---|---|---|---|
| V-1 | 标定教训（非缺陷） | b_layer T3 总体比 0.5385 出锚 [0.495, 0.520]：试跑 N=8 只见近非混合 AR(1) q 链头部（a=−0.9964、混合窗 139），携带包络欠采样；链代数预测 0.553 包络内。教训：统计腿试跑长度须超过链混合时间（复用需 N≥139 或显式携带模型） | case_b_layer §4 在案；无回修项 |
| V-2 | 观察项（注册语义，非缺陷） | 驱动器新轨迹复位只清 sta%t（成员只写拥有分量）——单片段气相单员跑 q 跨轨迹链式携带；散射跑 incident 每轨迹重写 q 不可见。"每轨迹独立热化"语义属设计变更走流程 | case_b_layer §2 定量在案 |
| V-3 | 文档勘误（已闭环） | case_surface_window ω 换算方向（sqrt(K/(conv·m))=0.2028 → sqrt(conv·K/m)≈0.1957，入界窄 7%）+ 井势包络子项（~1.8e2 → ~118）+ case_b_layer T2b 谷相位 π/2→3π/2——三处 V.5 勘误批已修，锚值与判定不受影响 | 两 README 勘误注在案 |
| V-4 | 常数观察项（非缺陷） | 解析侧必须用立法换算对 23.0605×0.04184 = 0.96485132 内能/eV（force_facade 唯一换算点两字面量；CODATA 0.9648535 相对差 2.3e-6）——否则 1e-10 级锚假失败。后续判据的常数姿态沿此 | chain_freq §3 常数注 + §4 实证 |
| V-5 | 通道停放（非缺陷） | 装配期频率导出通道：机制存在（hessian_cd + 成员表）但程序内未接线（四证据：两 seam 头注/全树零生产调用者/driver 尾注/成员 fail-loud）；缺口 = 装配期自动导出（成员 init 经容器真力算 Hessian→EIGN→填表） | chain_freq §2 查明在案；候选后续件 |
| V-6 | 探针缺陷（scripts/ 范围，非主体） | probe_chain_freq 首版再展开例程 COM 索引笔误（越界读栈）——-fcheck=all 定位、正式跑前修复；成员 dist_normalmode 无辜 | chain_freq §4 诚实记录 |

**GWRITE_LEVEL 催办（第二次，随 V 收口）**：每键缺省 '1' 的缺省语义
裁决/文档调和（C7 移交用户）**仍悬**——S 期开工前宜裁决（缺省上调或
维持并在文档写明；裁决权在用户）。

**收尾移交（→S，2026-09-16）：**
- **能力面→证据对照表**（S 完工性复审的 V 相输入）：

| 能力面 | 验证载体 |
|---|---|
| 全四相循环（input→装配→采样→传播→记录→统计） | harmonic_scatter A1/A3（assembly 行逐字 + 散射分类面） |
| verlet 帧级守恒 + 确定性 | harmonic A1/A2（1.2648e-03 对锚 2.0e-1；双跑逐字节） |
| radau 自适应守恒 + 步长序列通道 | radau R1/R3（1.1e-13 对锚 1.0e-6；首差 = DT、×1.4 几何增长） |
| 表面范式端到端（POSCAR 胞/A_LAT/BOXPAD/单片段弹丸/速度块跳过/无分类合法跑） | surface A1-A6（paradigm=1 装配 + thta 真值面 + e_rel 退化对偶面 + 组成法令负锚） |
| B 层三员真力链（thermalize/thermostat/glo_target，C5-R1） | b_layer T/H/G 系（χ² 精确带 0.954/0.956 + FDT 位级 14/14） |
| 链级频率真值（Hessian→谱→boltzmann/normalmode，C5-R1） | chain_freq F/B/N 系（谱位级 1e-13 级 + 统计带全绿） |
| 换算常数立法（force_facade 唯一换算点） | chain_freq 常数注（立法对 0.96485132 实证）+ 各 case e0 复算 |
| 旧链位级等价（参考同量对照，非门） | check_facade (s) 双编译器（gfortran C2 + ifx V.5 四配置） |

- **停放项终局**：V 相新增/确认四条——① 装配期频率导出通道未接线
  （V-5，候选后续件，chain_freq §2 四证据）；② ifx 下 check_facade 两
  红腿均为判据脆性/既有停放族（(g) 严格 == 的 1 ulp 构造性不等价；
  (r) 2x2 eig_sym 算术序 1.29e-15——C2 停放区既有条目的 ifx 量值）；
  ③ C8→V 第 5 条采样器→统计相初值快照通道（初值参考量 VI/OAMI/AMAI/
  AMBI 与转动能重采样均值 SDA/SDB 的生产者；落地即解锁 fin%ang(4..18)）
  V 期未落地——后验无可行使面，随其 C8 触发条件继续停放，→S 可见；
  ④ scripts/ 改进注记：check_facade (g) FAIL 行不打印 got-vs-expected
  （其余容差判据均打印）——候选后续件，暂不改代码。
  既有登记不动：C8 终审停放区（hist_thta 气相门控候选/col_occ 守卫/
  fin_t 初值参考量/SDA/SDB/plot_hist allocated 守卫/lattice%q_ref）、
  C2 停放区（(r) 1e-15 性状/fssh/iesh 后窗件/r_adiab 构造通道/有状态
  成员 wrapper 终形）、b_layer NU_COLL 携带无消费面注记——S 复审按
  各处登记清点。
- **GWRITE_LEVEL 仍悬用户**（上文催办）。

**风险：**
- B 层统计判据的容差标定（样本量 vs 运行时长）——先小样本试跑再定
  锚，写明；
- 表面 case 的 POSCAR 样件与 A_LAT 一致性（几何自洽）。
