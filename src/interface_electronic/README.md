# interface_electronic/ — 电子层面扩展成员（两轴同家）

规范 §3.2-2 电子层面的两条独立轴的可换成员同居此目录，以文件前缀区分：

| 轴 | 前缀 | 接口 | 现有 / 计划 |
|---|---|---|---|
| H 场提供方（可叠加贡献族） | pes_* | 甲类：facade_reg_pes（单入口，可选 repr=r_scalar〔缺省〕/r_diab/r_adiab——E10 归一：表象是登记数据不是注册分叉）。签名 pes_proc_i(q, mass, ws)——成员按自身 ws 填充模式交付：scalar 填 g_add[eV/Å]/v_add[eV]；diab 填 H 矩阵 + 矩阵元梯度（eV、eV/Å）；adiab 填各面 V_i/力核 + 当前面导数耦合行 d_i,·（经 ws%i_want 请求——i_surf 不进签名）。×23.0605×C1 换算集中在 facade（单位契约，主计划 §5 N1）；多面主势不做加性叠加 | C2 六件：reference / emt_nn / rst / glo_spring / langevin_ghost / ldfa；多面成员随 elec_tsh/elec_iesh 消费者引入 |
| 电子传播子（量子-经典方法） | elec_* | 乙类：facade_reg_elec（单选）。签名 elec_proc_i(s, f, dt, ws)——E9 加步长（跃迁概率 g_ij(dt) 所需）、E10 加工作区（只读消费：按 ws%repr 解释容器，成员不支持的 repr → loud-stop）。绝热 = 退化点（ds/dt≡0、力规则恒等、忽略 dt 与 ws），接口在退化点良定义 | elec_adiabatic（先行）；IESH / FSSH 将来同住 |

**两轴独立**（任意 PES × 任意方法合法配对，笛卡尔积结构）由代码机制保证，
与目录编组无关：模块互相独立（pes_* 与 elec_* 零引用）、两个注册表
（reg_pes / reg_elec 分立）、facade_init 能力级协商（非法组合 loud-stop
并指名冲突两轴）。

**E9 多面升级 → E10 工作区修订（2026-09-07）**：E9 引入甲类双签名（透热/
绝热）+ 乙类签名扩 dt + 协商表数据化 triad_tbl（四行：单面绝热现行 + 多面
tsh/iesh 预放行）。E10 公理化修正（三公理：接口不变性/黑箱性/物自体原则）
——甲类签名完全归一 pes_proc_i(q, mass, ws)：表象构件（H_ij/d_ij/V_i）
收纳于 ws_t 容器（facade/elec_ws.f90——表象从签名分叉降级为 ws%repr 数据
标记），不进 state_t、不出门面箱外、生命周期一步；双注册例程撤销
（facade_reg_pes 单入口带可选 repr）。**radau 在多面行恒假**不变（自适应
序列边界"一步"无良定义，跃迁概率 g_ij(dt) 失据，非绝热只配定步长积分器）。

接入 = 新文件 + 装配相位一行注册，既有代码零触及（T2 判据演练点）。

参考质料（旧库只读）：interface_TEST.f90 / EMT-NN / RST / VERLET.f90 内嵌力。
