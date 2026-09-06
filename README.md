# VENUS2026 项目

VENUS 程序精简与重构的独立工作项目。

## 布局

```
venus2026/    程序源码（自 C@Au_TDHF/venus2026_stage 迁出，git 管理文件
              223 个，工作产物 fort.*/stdout* 未迁移；原仓库含 .git 历史，
              远程 github.com/ComRADKhrushchev/VENUS2026.git）
docs/         设计文档
  design_spec_completeness.md   程序完备性设计规范（推演文档）
```

## 当前阶段任务

精简 VENUS2026：削除冗余功能。架构裁决依据
docs/design_spec_completeness.md——四类纯粹子例程 + 控制核的分解、
不可约接口 A（力贡献门面，内部双注册轴）与 B′（初态分布族）、
封闭选择族（核传播子、终止谓词文法、结局分类数据）。

## 注意

- venus2026/ 目前无 .git（迁移时仅复制了受管文件）；需要历史时回
  C@Au_TDHF/venus2026_stage 取，或在原仓库继续用 git 操作后同步
- 构建与运行约定见 venus2026/CLAUDE.md（Makefile 单入口、ifx、
  内部单位 Å/amu/10 fs、KEYWORD=VALUE 输入）
