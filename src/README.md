# src/ — VENUS2026 完备化新程序

VENUS 化学动力学程序的完备化重构：绝热 QCT，L0-L5 分层（层概念见下表，
目录名不带层标号）。

**唯一执行依据：[docs/final_completion_plan.md](../docs/final_completion_plan.md)**
（目标终态 §1、执行序列 §2、验收量尺 §9）。架构裁决标准见
docs/design_spec_completeness.md，目标结构与质料映射见 docs/optimal_structure.md。
旧库 venus2026/ 只读参考（算法与公式出处），不迁移、不修改。

## 布局（入口在根；interface_ 前缀 = 可扩展注册族成员）

```
src/
  driver.f90  input.f90                          入口：主控编排（装配→演化→统计）+ 输入文法
  globals/     consts  state  control  config  config_roster
                                                  变量与规定（全体 use 的根基）
  facade/      force_facade  sampler             物理总入口：力门面 + 抽样框架
  dynamics/    predicate  propagator             封闭族（谓词文法、三积分器；成员固定）
  interface_electronic/    elec_adiabatic        电子层面两轴成员家（乙类先行；
                                                  甲类 pes_* C2 落此——见该目录 README）
  interface_distributions/ dist_boltzmann  dist_ebk  dist_incident
                                                  初态分布成员（B′ 扩展；示例 3 个，余者 C5）
  output/      recorder  final_state  outcome  density  plot_hist
                                                  在线记录（列注册）+ 离线后处理
  utils/       rng  linalg  geometry             仪器（单一随机流、对角化、旋转）
  Makefile  README.md
```

| 层 | 目录 | 职责 |
|---|---|---|
| L0 控制核 | 根 driver.f90 | 装配→演化→统计的调用序列，无物理 |
| L1 变量与规定 | globals/（input 留根） | 物理态容器、控制流、固化导出规定、原子清单 |
| L2 门面与框架 | facade/ | (Q,P,s)→(F,E) 门面（E9 甲类三签名：单面标量/多面透热/多面绝热；协商表 triad_tbl 数据化）；抽样框架 |
| L3 封闭族 | dynamics/ | 谓词文法、传播子（固定选择，非扩展点） |
| L3 注册族成员 | interface_*/ | 两目录三族：electronic（甲类 pes_* H 场提供方 + 乙类 elec_* 电子传播律——规范 §3.2-2 电子层面两轴同家）/ distributions（B′ 初态分布）——新文件 + 一行注册即接入 |
| L4 记录与后处理 | output/ | 轨迹序列化（列注册）、终态/结局/态密度/绘图（可离线） |
| L5 仪器 | utils/ | 随机流（单一）、对角化、旋转 |

分层要点：每层只向下依赖，注册族成员互不依赖（成员内部可调 utils 仪器）；
facade 是控制核看到的全部物理；单位换算只发生在 force_facade 边界；
全库单一随机流（rng.f90）。

## 构建（P0 骨架）

```bash
make            # 单次 ifx 调用编译+链接，产出 venus.e（SRC 填入后生效）
make clean      # 删除 venus.e 与 .mod/.obj/.pdb/.ilk
```

- 编译器：Intel `ifx`（oneAPI 2025.2 验证）；flags `-r8 -double-size=64 -i8 -O2 -w`。
- `SRC` 暂为空列表，K3 阶段（新树独立编译）按依赖序填入
  （utils → globals → facade → dynamics → interface_* → output → input → driver）。
- 新程序无 MKL 依赖；Windows 链接器遮蔽处理（`-Qlocation,link`、`LIB` 导出、
  `MSVC_BASE`/`WINKIT_LIB`/`INTEL_BASE` 覆盖机制）移植自旧库 Makefile 平台段，
  用法见 venus2026/CLAUDE.md「Build」节。
