# src/ — VENUS2026 完备化新程序

VENUS 化学动力学程序的完备化重构：绝热 QCT，L0-L5 分层。

**唯一执行依据：[docs/final_completion_plan.md](../docs/final_completion_plan.md)**
（目标终态 §1、执行序列 §2、验收量尺 §9）。架构裁决标准见
docs/design_spec_completeness.md，目标结构与质料映射见 docs/optimal_structure.md。
旧库 venus2026/ 只读参考（算法与公式出处），不迁移、不修改。

## 布局（L0-L5 文件一览）

| 层 | 文件 | 职责 | 状态 |
|---|---|---|---|
| L0 控制核 | driver.f90 | 装配→演化→统计的调用序列，无物理 | K 阶段生成骨架 |
| L1 物自体容器 | state.f90 | Q, P, s, F, t（可分配，维数=名册内容参数） | K 阶段生成骨架 |
| L1 表象层 | config.f90 | 场景、名册、分布、谓词、结局分类、观测量配置 | K 阶段生成骨架 |
| L2 接口 A 门面 | force_facade.f90 | (Q,P,s)→(F,E)；内部 Σ力贡献成员 × 单选电子传播律 | K 阶段生成骨架 |
| L2 接口 A 乙类退化成员 | elec_adiabatic.f90 | ds/dt≡0、力规则恒等（显式存在的退化点） | K 阶段生成骨架 |
| L2 接口 B′ 框架 | sampler.f90 | 分布→(Q0,P0,s0) 通用抽样器 | K 阶段生成骨架 |
| L3 注册族成员 | pes_*.f90 / dist_*.f90 | 可换实现体（新文件 + 一行注册即接入） | K 阶段生成骨架 |
| L3 封闭族 | propagator.f90 / predicate.f90 | 力泛函→步进（有限选择）；状态→布尔（文法解释器） | K 阶段生成骨架 |
| L4 测量 | recorder.f90 | 轨迹序列化 + 谓词求值 | K 阶段生成骨架 |
| L4 后处理 | analysis/ | 结局分类执行、统计、绘图（可离线），见 [analysis/README.md](analysis/README.md) | K 阶段生成骨架 |
| L5 仪器 | rng.f90, linalg.f90, geometry.f90 | 随机流（单一）、对角化、旋转 | K 阶段生成骨架 |

分层要点：每层只向下依赖，L3 成员互不依赖（注册族力学）；L2 门面是控制核
看到的全部物理；单位换算只发生在 force_facade 边界；全库单一随机流（rng.f90）。

## 构建（P0 骨架）

```bash
make            # 单次 ifx 调用编译+链接，产出 venus.e（SRC 填入后生效）
make clean      # 删除 venus.e 与 .mod/.obj/.pdb/.ilk
```

- 编译器：Intel `ifx`（oneAPI 2025.2 验证）；flags `-r8 -double-size=64 -i8 -O2 -w`。
- `SRC` 暂为空列表，K3 阶段（新树独立编译）按 L5→L1→L2→L3→L0 依赖序填入。
- 新程序无 MKL 依赖；Windows 链接器遮蔽处理（`-Qlocation,link`、`LIB` 导出、
  `MSVC_BASE`/`WINKIT_LIB`/`INTEL_BASE` 覆盖机制）移植自旧库 Makefile 平台段，
  用法见 venus2026/CLAUDE.md「Build」节。
