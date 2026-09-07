# interface_electronic/ — 电子层面扩展成员（两轴同家）

规范 §3.2-2 电子层面的两条独立轴的可换成员同居此目录，以文件前缀区分：

| 轴 | 前缀 | 接口 | 现有 / 计划 |
|---|---|---|---|
| H 场提供方（势能面的推广） | pes_* | 甲类：facade_reg_pes（可叠加贡献族）。成员交付 **eV 力核**（g_add[eV/Å], v_add[eV]）；×23.0605×C1 换算集中在 facade（单位契约，主计划 §5 N1） | C2 六件：reference / emt_nn / rst / glo_spring / langevin_ghost / ldfa |
| 电子传播子（量子-经典方法） | elec_* | 乙类：facade_reg_elec（单选）。绝热 = 退化点（ds/dt≡0、力规则恒等），接口在退化点良定义 | elec_adiabatic（先行）；IESH / FSSH 将来同住 |

**两轴独立**（任意 PES × 任意方法合法配对，笛卡尔积结构）由代码机制保证，
与目录编组无关：模块互相独立（pes_* 与 elec_* 零引用）、两个注册表
（reg_pes / reg_elec 分立）、facade_init 能力级协商（非法组合 loud-stop
并指名冲突两轴）。

接入 = 新文件 + 装配相位一行注册，既有代码零触及（T2 判据演练点）。

参考质料（旧库只读）：interface_TEST.f90 / EMT-NN / RST / VERLET.f90 内嵌力。
