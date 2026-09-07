# interface_electronic/ — 电子层面扩展成员（两轴同家）

规范 §3.2-2 电子层面的两条独立轴的可换成员同居此目录，以文件前缀区分：

| 轴 | 前缀 | 接口 | 现有 / 计划 |
|---|---|---|---|
| H 场提供方·单面标量 | pes_* | 甲类：facade_reg_pes（可叠加贡献族）。成员交付 **eV 力核**（g_add[eV/Å], v_add[eV]）；×23.0605×C1 换算集中在 facade（单位契约，主计划 §5 N1） | C2 六件：reference / emt_nn / rst / glo_spring / langevin_ghost / ldfa |
| H 场提供方·多面透热 | pes_* | 甲类：facade_reg_pes_diab（主势专用——多面行不做加性叠加）。交付 H 矩阵 + 矩阵元梯度（eV、eV/Å）；绝热视图消费方经 utils/linalg 对角化现算 | E9 预留；随 elec_tsh/elec_iesh 消费者引入 |
| H 场提供方·多面绝热 | pes_* | 甲类：facade_reg_pes_adiab（主势专用）。交付各面 V_i/力核 + 当前面导数耦合行 d_i,·（eV、eV/Å、1/Å——O(ndof·n_surf) 而非全对） | 同上 |
| 电子传播子（量子-经典方法） | elec_* | 乙类：facade_reg_elec（单选）。签名 elec_proc_i(s, f, dt)——E9 加步长（跃迁概率 g_ij(dt) 所需；绝热忽略）。绝热 = 退化点（ds/dt≡0、力规则恒等），接口在退化点良定义 | elec_adiabatic（先行）；IESH / FSSH 将来同住 |

**两轴独立**（任意 PES × 任意方法合法配对，笛卡尔积结构）由代码机制保证，
与目录编组无关：模块互相独立（pes_* 与 elec_* 零引用）、两个注册表
（reg_pes / reg_elec 分立）、facade_init 能力级协商（非法组合 loud-stop
并指名冲突两轴）。

**E9 多面升级（2026-09-07）**：甲类双签名并存（透热/绝热，成员按数据源天然
形态自选）+ 乙类签名扩 dt + 能力级协商数据化 triad_tbl（四行：单面绝热现行
+ 多面 tsh/iesh 预放行；**radau 在多面行恒假**——自适应序列边界"一步"无
良定义，跃迁概率 g_ij(dt) 失据，非绝热只配定步长积分器）。

接入 = 新文件 + 装配相位一行注册，既有代码零触及（T2 判据演练点）。

参考质料（旧库只读）：interface_TEST.f90 / EMT-NN / RST / VERLET.f90 内嵌力。
