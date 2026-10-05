# systems/inactive/ —— 容器收存区（build glob 之外）

**用户裁决（2026-09-21）：所有非 reference 容器平时一律收纳于本目录；需要测哪个
才把哪个取出。**

本目录比 `systems/<sys>/` 深一层，`src/Makefile` 的 `systems/*/*.f90` 通配符不会
聚合到这里——收存于此的容器**不进入构建**。`systems/` 顶层常驻 `reference/`
（五 case 的生产容器）；测试某容器时把该文件夹移到 `systems/` 顶层（拖入即安装，
E27"文件夹在场即安装"立法不变）→ `make` 两轮 → 测试 → 移回本目录。

注意：venus.e 立法为 one-container-per-process——目标容器与 reference 并存时
`registry_gen_all` 装配期即 double-bind STOP 1。**判据程序（scripts/check_*）独立
编译容器，不受此限**；若需用 venus.e 实跑目标容器，测试期间 reference 一并收入
本目录（换位），或待多容器选择机制裁决（docs/plans/S.md 停放项表；本案分析见
docs/cau_container_study.md §6.2）。

（沿革：本目录前身为根目录 `containers/` 收存库 + 原型停放区双轨，2026-09-21
裁决后合并为一。）

收存清单（C@Au TDHF 四系列全套，2026-09-21 一次建成）：

- cau_emt2/ — C@Au(111) EMT_2spec 系列（判据 scripts/check_cau_emt2.f90
  88 PASS/0 FAIL；EMT 解析 + Hubbard 密度多项式交付；含 vendored 引擎
  PBC 最小镜像修复，见文件头 DEVIATION）
- cau_emtnn/ — C@Au(111) EMT-NN 系列（判据 scripts/check_cau_emtnn.f90
  84 PASS/0 FAIL；EMT+pairNN 势 + Hubbard NN 交付；引擎 pbc_wrap 同款
  最小镜像修复；判据曾暴露的 slab 力失配系判据几何 bug，引擎无误）
- cau_rst/ — C@Au(111) RST-Hubbard 直连系列（判据 scripts/check_cau_rst.f90
  84 PASS/0 FAIL；BvK 板弹性 + 两两 NN + softplus 交付链；三处明记
  偏离：slab 同步/梯度链修复（旧接口调试旗标 zero_dh/zero_dU 与截断
  dV_orb）/J 钳位梯度；数据件 Slab.xyz + nn_weights_rst_hubbard.txt +
  Analytic_Potential.txt 随容器夹）
- cau_peri/ — C@Au(111) periodic 2D 解析系列（判据 scripts/check_cau_peri.f90
  84 PASS/0 FAIL；三原子模型 + 晶场劈裂 h 交付；vendored 引擎两处
  梯度符号修复（rho=|z|-p5 约定）+ 隔离改写（去 venus 宿主耦合），
  见 surface_potentials.f90 头注；Analytic_Potential.txt 随容器夹）
- （四容器共同能力：`ELEC_METHOD = tdhf_hub` 时方法包经 sync 成员把电子态发布进
sta%s%rho（自旋块 (2·n_siz)²），并注册五个 obs_list-only 记录列——OBS_LIST =
n_imp_up, n_imp_dn, mag_mom, e_elec_ev, e_mf_ev；C@Au 研究记录与四系列移植
要点：docs/cau_container_study.md）
- fh2_leps/ — F+H2 改型 LEPS 容器（判据 scripts/check_fh2_leps.f90；自会话早期
  状态逐字恢复，恢复后判据全量复跑通过）
- h3_bkmp2/ — H+H2 BKMP2 容器（判据 scripts/check_h3_bkmp2.f90；势为存档固定
  格式 systems/bkmp2.f 的自由格式逐位等价转录）
- h2ag111PES/ — H2/Ag(111) 神经网络容器（判据 scripts/check_h2ag111.f90；原分发
  包四件 + 容器三件同居；运行目录需携带 weights/biases-h2ag.txt 副本；谱导出：
  气相 H2 文献常数 w_e=4401.21 cm^-1，经 spectrum_stage 导出通道）
