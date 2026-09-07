# src/ — VENUS2026 完备化新程序

VENUS 化学动力学程序的完备化重构：绝热 QCT，L0-L5 分层。

**唯一执行依据：[docs/final_completion_plan.md](../docs/final_completion_plan.md)**
（目标终态 §1、执行序列 §2、验收量尺 §9）。架构裁决标准见
docs/design_spec_completeness.md，目标结构与质料映射见 docs/optimal_structure.md。
旧库 venus2026/ 只读参考（算法与公式出处），不迁移、不修改。

## 布局（最底层 L1 变量与规定在根目录；其余按层入子目录）

```
src/
  consts.f90  state.f90  control.f90  config.f90    ← 根：全局变量与规定（全体 use 的根基）
  config_roster.f90  input.f90                       ← 根：原子清单 + 输入文法
  l0_driver/   driver.f90                            主控（装配→演化→统计，无物理）
  l2_facade/   force_facade.f90  elec_adiabatic.f90  sampler.f90
                                                     力门面（Σ成员+能力级）、退化电子成员、抽样框架
  l3_family/   predicate.f90  propagator.f90         封闭族（谓词文法、三积分器）
               dist_*.f90                            初态分布成员（示例 3 个；余者 C5）
               （pes_*.f90 力成员 C2 落此）
  l4_output/   recorder.f90  final_state.f90  outcome.f90  density.f90  plot_hist.f90
                                                     记录（列注册）与后处理（离线）
  l5_tools/    rng.f90  linalg.f90  geometry.f90    仪器（单一随机流、对角化、旋转）
  Makefile  README.md
```

| 层 | 内容 | 职责 |
|---|---|---|
| L0 控制核 | l0_driver/ | 装配→演化→统计的调用序列，无物理 |
| L1 变量与规定（根目录） | consts/state/control/config/config_roster/input | 物理态容器、控制流、固化导出规定、原子清单、输入文法 |
| L2 门面与框架 | l2_facade/ | (Q,P,s)→(F,E) 门面；退化电子成员；抽样框架 |
| L3 封闭族+注册成员 | l3_family/ | 谓词文法、传播子；pes_*/dist_* 可换成员（新文件+一行注册即接入） |
| L4 记录与后处理 | l4_output/ | 轨迹序列化（列注册）、终态/结局/态密度/绘图（可离线） |
| L5 仪器 | l5_tools/ | 随机流（单一）、对角化、旋转 |

分层要点：每层只向下依赖，L3 成员互不依赖（注册族力学）；L2 门面是控制核
看到的全部物理；单位换算只发生在 force_facade 边界；全库单一随机流（rng.f90）。

## 构建（P0 骨架）

```bash
make            # 单次 ifx 调用编译+链接，产出 venus.e（SRC 填入后生效）
make clean      # 删除 venus.e 与 .mod/.obj/.pdb/.ilk
```

- 编译器：Intel `ifx`（oneAPI 2025.2 验证）；flags `-r8 -double-size=64 -i8 -O2 -w`。
- `SRC` 暂为空列表，K3 阶段（新树独立编译）按 L5→L1→L2→L3→L0 依赖序填入
  （根目录 L1 六件在 L5 之后、子目录之前）。
- 新程序无 MKL 依赖；Windows 链接器遮蔽处理（`-Qlocation,link`、`LIB` 导出、
  `MSVC_BASE`/`WINKIT_LIB`/`INTEL_BASE` 覆盖机制）移植自旧库 Makefile 平台段，
  用法见 venus2026/CLAUDE.md「Build」节。
