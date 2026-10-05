# C@Au TDHF 容器化研究报告（四系列全套落地）

日期 2026-09-21 · 请求：研究 C@Au_TDHF/2D 的 EMT_2spec / EMT_NN / periodic / RST 四系
PES + TDHF_U 动力学能否建成 VENUS2026 的一系列 systems/ 容器，并尝试落地。
**同日续批：四系列容器一次建成**（cau_emt2/emtnn/rst/peri，判据
88+84+84+84 PASS / 0 FAIL，gfortran+ifx 双编译零错；见 §九）。

## 一、总结论

**可行**。E27 裁决（容器=组合根，体系专属电子方法包同居容器文件夹）正是为此形态
设计的——"自定义 Hubbard 哈密顿量压测"即本体系。本次以 **EMT_2spec 系列落地原型容器
`cau_emt2`**（PES 腿 + Hubbard 交付 + TDHF 电子方法包 + 终止判断），判据程序
`scripts/check_cau_emt2.f90` **84 PASS / 0 FAIL**，gfortran 与 ifx 双编译零错。
容器现收存于 `systems/inactive/cau_emt2/`（构建 glob 之外；移回 `systems/` 即激活）。

## 二、源码勘察与映射

源树 `C:\Users\13614\Desktop\MAIN\3\C@Au_TDHF\2D`（主程序 VENUS.f90 + IESHParas +
TDHF_U + src_NN 平铺接口）。全部系列共用旧 VENUS 钩子约定，映射到 E27 契约：

| 旧钩子 | E27 契约位 |
|---|---|
| POTPRE（读参数/slab） | reg _init 装配缝（参数文件通道，`<sys>_pars.txt`） |
| POT0(NDUM,V) + DPESHON | `container_pes_i(q,mass,g,v)` 组合力条目（eV 全赋值） |
| H_Orb/V_orb/dH_orb/dV_orb | 容器 pot 层的 Hubbard 交付（h/U/J/Vcp + 核），由组合 pes 与方法包共享 |
| Init_IESH/set_occ/buildHcore/buildFock/evolve_rhom/TDHF_Evolution | 容器文件夹电子方法包 `method_reg('tdhf_hub', prop, match)`（ELEC_METHOD 选用） |
| CALTYP=−1 绝热参考 | ELEC_METHOD=adiabatic → 包不武装 → pes 退化为纯 EMT（判据 i1 精确验证） |
| Calc_Energy/Calc_force 折叠 | 组合 pes 内折叠（v = E0+E_elec+E_mf−E_metal_zero；判据 h 精确验证 v/g 配对） |

四系列 + 第五变体的移植要点：

| 系列 | 引擎 | 依赖 | 运行材料 | 备注 |
|---|---|---|---|---|
| EMT_2spec | pes_emt_mod（零依赖，本次已移植） | 无 | emt_2spec_params.txt（EMT14+Hubbard32）+ Slab.xyz | 纯解析；可动 slab；**原型载体** |
| EMT_NN | emt_nn_pes_venus（自足模块） | 无 | nn_weights_emt_nn.txt + Slab.xyz | h/U/J/Vcp 走 NN 段；梯度有解析装配 |
| periodic | surface_potentials | venus_params/venus_data/input_parser（需隔离改写） | nn_weights.txt（可选） | 2D 解析晶场 + Morse/Rydberg；杂质量坐标 (u,v,z) |
| RST | rst_pes（零依赖） | 无 | nn_weights_rst.txt（7 态权重）+ Slab.xyz | 两两 NN + BvK 板弹性；多透热态交付 |
| RST_Hubbard | rst_hubbard_pes（自足内联） | 仅 IESHParas | NN 权重 + slab | h/U(softplus)/J(softplus)/Vcp/GS 直供 |

**共性**：TDHF 引擎（浴离散、Fock、Magnus-2、均场折叠）对四系列是同一套，给定
{h, U, J, Vcp + 各自梯度} 即通用。**用户裁决（2026-09-21）：TDHF 方法留在容器文件夹
内，不升格 methods/**——升格需要"体系→方法"通用交付接口（E27 已退役的通用语），且
ρ 生命周期/组合折叠与容器私有模块耦合。代价是各 C@Au 容器各持一份引擎拷贝
（registry 族成员不互相依赖）；接受为有意识的取舍。

## 三、原型容器构成（systems/inactive/cau_emt2/）

- `pes_emt_mod.f90` — 原样引入的 EMT 引擎（含一处明记修复，见 §4.1）
- `cau_emt2_pot.f90` — 纯物理层：E0+核、h/U/J/Vcp 值+核（eV 交付，密度缓存共享）
- `cau_emt2_tdhf.f90` — TDHF 引擎 + 组合 pes（武装态折叠均场项）+ 'tdhf_hub' 包成员
  （prop = Magnus-2 子步；match = Ehrenfest 无跳跃 no-op）+ ρ 生命周期
- `cau_emt2_box.f90` — cau_emt2_pars.txt 参数通道（KEYWORD = VALUE）+ Z_TERM 散射终止
- `reg_cau_emt2.f90` — container_bind_pes + method_reg('tdhf_hub') + 装配缝武装

常数链：相位 = E[eV]×23.0605×e_conv×dt/hbar_code（consts 立法对 0.96485132/0.063508，
即旧 C1/C7 链，dt 为 10 fs 内部单位）；kT = r_kcal×T/23.0605；μ=带心（旧 set_occ(.,0)）。
单位姿态：容器内全程 eV，换算只在 force_facade（T4 不变量保持）。

## 四、实质发现（三件）

### 4.1 vendored 引擎 PBC 最小镜像缺陷（已修，文件头明记）

`compute_pairwise_sum/grad`、EMT 主环、`pbc_min_image_2d` 四处包装用的是分数坐标
`(f−nint(f))` 矩形包装——对 60° 斜胞**不是真最小镜像**：分数角区原子重建距离可达
半对角（17.7 Å 六角胞实测 15.33 Å vs 真最小 8.85 Å）。后果：Hubbard 两两和在镜像
翻转面处跳变（探针实测 h 值跳 0.0117 eV，FD 发散至 −58 eV/Å）且角区原子系统错值
（gentle 0.3 截断下 θ(8.85)≈0.25，贡献不可忽略）。对合位形（C 恰在 6×6 胞顶位）时
ii=3/jj=3 列原子恰落边界 → 症状最烈。E0 主环因硬截断(1.5·rcut)先于角区贡献而幸免，
仅余 1.86e-4 eV/Å 的边界 nint 不对称力。**修复**：新增 `min_image_true`（9 邻域分数
移位取最短重建），四点统一替换。修复后对称位形 ps(x±1e-6) 逐位相等；**此修复使容器
不再逐位复现旧工程轨迹**（旧结果内嵌此缺陷），判据按物理正确性（FD/对称性）验证。
注意：EMT 主环角区错距对 θ(≈7.7 Å, acut=24.93)≈e^-64 无能量影响，修复对 E0/F0 仅
消除边界不对称。

### 4.2 引擎能量输出的 ∝h FD 反常（解析核为准，容差据此定）

E0 的中心差分与解析核之差精确线性于 h（系数 ~15 eV/Å²，z 探针；h=1e-6 时两者 9 位
一致）——引擎能量输出含 ~1e-11 eV 级非解析残余，解析核正确。判据 (c) 容差按实测
5e-3·max(1,|g|)，头部注明。旧工程自己的 fd_check 以 `if (fd_count >= 0) return` 禁用
——佐证其力/能自洽性从未严格验证。

### 4.3 ensemble-TDHF 均场无严格总能量守恒（旧工程亦然）

固定几何下 GS 种子的 E_elec+E_mf 从 −51.6 eV 弛豫至 +1.0 eV 后稳定（物理弛豫瞬态，
非积分器错误；ρ 迹/厄米性全程守恒）。冻结 ρ 的武装折叠守恒至 6e-3 kcal/mol（判据
B1 腿）——v/g 配对正确；含演化的轨迹漂移 ~100 kcal/mol 与 dt、耦合强度均不敏感，
为此均场方案固有。旧工程 IESHParas 的 `energy_feedback`（速度重标定）正是对此的
生产补丁，未移植。判据 (i2) 取有界性界（<200 kcal）+ 迹守恒 + 实测值报告，头部如实
记载。后续若需严格守恒：能量重标定补丁或更细的核-电交错（见 §6 缺口）。

## 五、判据证据（scripts/check_cau_emt2.f90，84 PASS / 0 FAIL）

(a) Gauss-Legendre 闭式（4 轨道探针浴）；(b) 生产浴对称带+零能恒等式；(c) EMT 核
vs FD（一般位形，含引擎反常容差说明）；(d) h/U/J/Vcp 核 vs FD（一般位形，1e-5）；
(e) Hcore 对称/分支镜像/杂质块对角；(f) GS/1D 种子占据 + FD 浴；(g) 100 子步幺正
（迹/厄米）；(h) **组合折叠 v/g 配对**（1e-5，核心判据）；(i1) 绝热简并点精确退化
+ 300 步漂移 <0.05 kcal；(i2) 武装轨迹有界 + 迹守恒。EMT/H 参数锚点 = 生产
emt_2spec_params.txt 全集内置。

## 六、契约缺口与新决定事项候选

1. **轨迹头容器钩子缺失**（本次实证）：封闭五入口演化环（prop_step→elec_prop→
   mqc_statistic→rec_frame→container_term）没有任何"新轨迹开始"信号可达 ρ 承载容器；
   组合力条目只收 q。原型以兜底实现（终止真值时重播种 + prop 成员时间倒退守卫），
   max_steps 终止轨迹的首力仍见旧 ρ。候选后续件：samp_run 头部的容器复位钩子。
2. **多容器选择停决事项——运营约定已立（2026-09-21 用户裁决），机制裁决仍停放**：
   用户裁决"所有非 reference 容器平时一律收纳于 systems/inactive/，需要测哪个才
   把哪个取出"（inactive 比 systems/<sys>/ 深一层，构建 glob 不聚合）。据此
   fh2_leps / h3_bkmp2 / h2ag111PES / cau_emt2 四容器已归位 inactive/，根目录
   containers/ 双轨收存库撤销，venus.e 恢复 reference 单容器绿态（五 case 冒烟
   漂移与会话基线逐位一致）。取出测试时的并存 double-bind（one-container-per-
   process + 通配符全拖入）仍为机制层停决事项：判据程序路径独立编译不受影响；
   venus.e 实跑目标容器需测试期间将 reference 一并收入 inactive（换位），或裁决
   选择机制（输入驱动绑定守卫 / 每容器一可执行 / 条件化 *_all）。系列化（多容器
   同树常驻）前必须机制裁决。
3. occ_set/n_hop_set 统计通道对本方案（Ehrenfest 无跳跃）保持简并 [1]/0；电荷/自旋
   观测（杂质占据、磁矩）现经 `cau_emt2_diag_out` 内省出口，记录列化（rec_reg_col）
   与终态分类（container_bind_classify，需按轨迹缓存 ρ 观测）为自然后续件。

## 七、系列化路线图（其余三系列）

1. EMT_NN：引擎自足，容器 = 引擎 + pot(NN 交付) + box(权重文件+slab) + reg；
   TDHF 方法包按裁决随容器文件夹各持一份（从 cau_emt2_tdhf.f90 复制改接交付层）。
2. RST / RST_Hubbard：rst_pes 零依赖可直接引入；多透热态交付（7 态）在组合 pes 内
   按 ρ 折叠——需要把 RST 的态→(h,U,J,Vcp) 链并入交付层。
3. periodic：需先解除 surface_potentials 对 venus_params/venus_data/input_parser 的
   耦合（隔离改写量最大，建议最后做）。
4. 每容器判据程序沿用 check_cau_emt2 形态（FD/对称/幺正/简并点/轨迹五族）。

## 八、验证腿汇总

- 判据 84 PASS / 0 FAIL（gfortran，scratch 目录运行，锚点内置）；
- 容器五文件 gfortran + ifx(-O2) 编译零错（ifx 2025.2 全路径调用）；
- 容器激活态入树 `make` 绿（registry_gen/depend.mk 正确聚合 cau_emt2，双容器
  double-bind 行为与立法一致实测）；
- 容器收存后 `make` 绿；**未做** ifx 全树 make 与五 case 回归——src/ 树本批零改动，
  且 systems/ 顶层三容器并存为报告期用户态（见 §6.2），全树构建/回归在该裁决后做
  才有意义。


## 九、四系列全套落地（2026-09-21 续批）

四容器同居 `systems/inactive/`（收存区裁决），TDHF 引擎按裁决随各容器夹持拷贝
（`cau_<n>_tdhf.f90`，自 cau_emt2 模板机械改名；模板已推广 h(2)/dh(2,:) 缓存以
承载 periodic 的晶场劈裂交付）。全部数据件随容器夹（h2ag111 先例）。

| 容器 | 系列 | 判据 | 引擎要点与明记偏离 |
|---|---|---|---|
| cau_emt2 | EMT_2spec | 88/0 | PBC 真最小镜像修复（§4.1）；h(2) 推广回归后 88 |
| cau_emtnn | EMT_NN | 84/0 | pbc_wrap 同款修复；E0 为相对冻结裸板的相互作用能（常数平移语义，头注明记） |
| cau_rst | RST_Hubbard | 84/0 | ①slab 逐调用同步（旧接口 NN 见冻结文件几何）②梯度链修复（旧 dH_orb 调试旗标 zero_dh/zero_dU=.true. 与 dW 的 0.99 残迹；旧 dV_orb 只填 9 自由度）③J 钳位饱和梯度（新） |
| cau_peri | periodic | 84/0 | 隔离改写（去 venus_params/venus_data 宿主耦合，NDA3/NN_POT_TYPES 本地化，剔除 IESHParas debug 尾）；**引擎梯度符号修复**：rho=\|z\|-p5 约定下 d(rho)/dz=sign(z)、d(rho)/dp5=-1（原 sign(z-p5) 因子使面内梯度偏 -8/sqrt(3) 倍、z 分量翻号，FD 探针定位，修复后六函数 FD 一致） |

过程发现（判据方法论）：emtnn 的 slab 力"失配"最终定位为**判据程序自身的
slab 文件读取 bug**（首原子漏存）——引擎无误，slab 全动力学完好（用户异议
"slab 受力是 EMT 关键部分"成立，冻结板方案已撤销）。教训入档：FD 判据的
几何构造必须逐原子核对。

验证腿：四判据程序（FD/对称/幺正/简并点/轨迹五族，86 项均值）gfortran 全绿；
四容器物理栈 ifx 2025.2 按依赖序编译零错（16 文件）；判据程序运行材料
（Slab.xyz/权重/Analytic_Potential.txt/emt 参数锚点）全部自足（容器夹或内置）。


## 十、电子态发布与容器观测列（A+B 落地，2026-09-21 用户批准）

**B（发布）**：四容器的方法包现以 `method_reg('tdhf_hub', prop, match, sync)` 注册
sync 成员——把模块私有 rhom 按自旋块拷入 `sta%s%rho`（(2·n_siz)² 块对角，成员
声明下标语义，s%n_surf 同步为 2·n_siz）；elec_prop / mqc_statistic 的 sync 尾
每步自动调用，记录器看到的总是演化后电子态。1×1 退化种子到方法形状的重分配
即 state_destroy 注释预期的"成员初始化可重分配"路径。

**A（容器列）**：各 reg 缝在 tdhf_hub 选中时 `rec_reg_col` 注册五列——
`n_imp_up` / `n_imp_dn` / `mag_mom`（杂质占据与磁矩，读活的 rhom 对角）与
`e_elec_ev` / `e_mf_ev`（活 ρ × 缓存 Hcore 的 Tr 与均场能）。外部列为
obs_list-only（lvl_min = huge）：运行经 `OBS_LIST = n_imp_up, n_imp_dn, mag_mom,
e_elec_ev, e_mf_ev` 选用；提取器在 rec_frame 时刻执行，读到的 ρ 比 pes 缓存
新鲜半步（prop 在 rec 之前）。未选 tdhf_hub 的绝热跑不注册这些列。

**判据**：(j) 发布一致性（形状/n_surf/逐位拷贝）与 (k) 列值（GS 种子
2/0/2 精确 + 演化后有界 + 能量有限）四容器入档；总分 99+94+94+94 PASS/0 FAIL；
ifx 四栈复验零错。

（§6.3 预登记的"记录列化"后续件就此关闭；剩余：container_bind_classify 终态
分类需按轨迹缓存 ρ 观测——观测列已把数据送到 rec_final 侧，分类件可复用同
一发布通道。）


## 十一、第五容器 cau_tab：DFT 数据点表 + 1D Akima（2026-09-30 裁决）

**定位**：散射敏感性实验容器——把 DFT 数据点本身变成实验旋钮（手工调整表数据
的局部特征，重跑散射，测其对总体散射结果的影响）。PES 唯一来源是 `TABLE_FILE`
指定的表文件（六片 GS/1D/3P/mC/pC/VCP × 三位点 top/bridge/hollow × z 数据点，
格式立法见容器夹 `TABLE_FORMAT.md`），每位点 1D Akima 插值（C1 局部，稀疏点无
过冲），横向三系数傅里叶与 NN 路径同构（周期取表头 LATTICE 行）。越界政策：
低于首点 fatal 点名（不外推物理），高于末点末值保持+零梯度（与 6 Å 冻结线
自洽）。装载回显把每次运行的表形态写进 stdout 作出处。

**构成**：`cau_tab_tab.f90`（表装载/自检/Akima）、`cau_tab_pot.f90`（物理层，
warp/远场衰减/核装配自 cau_peri_pot 移植不变）、`cau_tab_tdhf.f90`
（cau_peri_tdhf 机械改名拷贝——与基线同引擎，A/B 对比唯一差异是 PES 源；
同步义务记头注）、`cau_tab_box.f90`（12 声明键 = cau_peri 13 键去 USE_NN/
NN_FILE 加 TABLE_FILE）、`reg_cau_tab.f90`、`system/` 两文件拷贝。

**与 cau_peri 的物理差异（明记）**：位点处两种构造重合（傅里叶在位点精确再现
位点值）；离点后表容器做"能量级傅里叶"、解析引擎做"参数级傅里叶"，一般不同
——建模选择非缺陷。判据 t3 在位点对拍锁定装配正确性。

**判据**（`scripts/check_cau_tab.f90`，判据程序同编译两容器）：(t1) Akima 数学
腿；(t4) 装载回读；(t2c/t2dh/t2h) 表核-值 FD 对拍 1e-5 档（探针 z 取段中点
2.2125——表格点上跨节点中差会吃 C1 斜率折痕 ~1e-4）；(t3a/t3b) 恒等表（解析
引擎三位点采样生成）vs cau_peri 在/离格点对拍（1e-9 / 2e-3，实测中点误差
|dE0|≤1.3e-4 eV、|dJ|≤2.3e-5 eV @0.05 Å 网格）；引擎腿 (a)(b)(e)(f)(f2)(g)
(i1)(i2)(i3)(h2)(o)(j)(k)(l) 自 check_cau_peri 移植同阈。**150 PASS / 0 FAIL，
gfortran+ifx 双编译同绿**；i1 无武装漂移 −0.019 kcal（C1 力 0.5 fs 步长下守恒
良好）、i2 武装漂移 −8.68 kcal（有界，家族量级）。

裁决记录：docs/plans/2026-09-30-cau-tab-container.md（六条：另立容器/数据后
供/三位点傅里叶/Akima/TDHF 副本随行/无修改通道直改表副本）。
