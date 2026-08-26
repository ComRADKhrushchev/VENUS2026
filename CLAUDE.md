# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VENUS2026 — VENUS 化学动力学 Fortran 程序的最小可运行形态（体系无关 TEST 势构建）。
只含经典动力学（`ELEC_METHOD=ADIABATIC`）；非绝热方法（TDHF-FSSH/IESH/MDEF）的
入口为 loud-stop 桩。C+Au 体系资产、NN PES、替代 PES 接口均不在本仓库
（完整形态见内部开发仓库）。

## Build

唯一构建入口是根目录 `Makefile`（无其他构建脚本）：

```bash
make            # 单次 ifx 调用编译+链接，产出 venus_test.e
make clean      # 删除可执行文件与 .mod/.obj/.pdb/.ilk
```

- 编译器：Intel `ifx`（oneAPI 2025.2 验证）。Linux 需先 `source setvars.sh`。
-  flags：`-r8 -double-size=64 -i8 -O2 -w`（默认实数/整数均 64 位）。
- MKL：Windows 经显式 lib 路径（`INTEL_BASE` 变量可覆盖）；Linux 经 `-qmkl=sequential`。
- **Windows 链接器遮蔽**：Git Bash 的 `/usr/bin/link` 遮蔽 MSVC link.exe，
  Makefile 用 `-Qlocation,link,"$(MSVC_BASE)/bin/Hostx64/x64"` 绕过并导出 `LIB`；
  版本路径不匹配时 `make MSVC_BASE=... WINKIT_LIB=... INTEL_BASE=...` 覆盖。
- 源文件按模块依赖顺序列于 `SRC`（单次调用，顺序敏感）。

## Run

```bash
cd cases/<name> && ../../venus_test.e
```

- 主输入**硬编码**读运行目录下的 `input_qct.txt`（`venus_input.f90:191`），
  不经 stdin、无命令行覆盖。
- `TASK=READ-QP`/`NORMAL-MODE` 的坐标（和动量）从 stdin 顺序读入。
- 程序在运行目录写 `fort.*` 产物（已被 .gitignore 覆盖）。

## Directory Structure

```
├── 根源码：VENUS.f90（主程序）、venus_params（常数）、venus_data（全局状态）、
│   venus_input（输入+初始坐标）、input_parser、param_mapping、presets、
│   harmonic_sampling、ENERGY.f、DVDQ.f
├── src_VENUS/   引擎子程序库（~60 个 .f/.f90）：
│   积分器 VERLET/RADAU/SYMPLE；初条件 SELECT/INITEBK/INITQP/GLPAR/
│   ORTHAN/THRMAN/PROBJ/QMMICRO/MICROCI；旋转 ROTATE*；分析 FINAL/GFINAL/EIGOUT；
│   恒温 THERMO/THERMBATH/GASDEV；浴 GLO；输出 writer TEST.f/GWRITE.f90/HIST.f90
├── src_TEST/    test_potentials.f90（HARMONIC/MORSE/LEPS 解析势+梯度）
│   + interface_TEST.f90（POTPRE 接口，运行时 PES 由 TEST_PES 关键字选择）
├── docs/        EXTENDING.md（接入自定义 PES/采样/积分器/电子方法的扩展指南）
└── cases/       16 个示例 case（README.md + input_qct.txt + results.txt + 验证图）
```

## Architecture

### Module dependency order（编译顺序敏感）

```
venus_params（纯常数）
  └─ venus_data（全局轨迹状态；allocate_venus_data(n) 在读 NATOMS 后调用一次）
       ├─ input_parser → param_mapping → presets → harmonic_sampling
       ├─ test_potentials（src_TEST）
       ├─ venus_input（读入+预设应用+关键字解析）
       └─ interface_TEST（POTPRE）
```

### Internal Unit System（改动代码必须遵守）

| 量 | 单位 | 量 | 单位 |
|---|---|---|---|
| 长度 Q | Å | 能量 | amu·Å²/(10 fs)² |
| 质量 W | amu | 动量 P | amu·Å/(10 fs)（P=W·v） |
| 时间 DT | 10 fs | 力 PDOT | code 能量/Å |

换算常数：`C1`=0.04184（kcal/mol→code）、`C7`=0.063508（ℏ code 单位）、
23.0605（eV→kcal/mol）、`C8`=1.987×10⁻³（k_B kcal/mol/K）。

关键规则：
1. POT0 返回 kcal/mol；动力学用 code 单位——力常数 k(kcal/mol/Å²) 须乘 C1；
2. 采样能量（kcal/mol）乘 C1 后才可算振幅/动量；
3. P²/(2m) 给 code 能量，与 kcal/mol 比较前须除 C1；
4. ENERGY.f 中 `T=T/C1` 后 H=T+V（同为 kcal/mol）。

## Input Format

仅支持 `KEYWORD=VALUE` 格式（`input_parser.f90`；旧顺序格式已删除，使用即报错停机）。
注释 `#`/`!`，续行 `,` 结尾。字符串关键字大小写不敏感；旧整数关键字回退。

**行内注释陷阱**：解析器只跳过整行注释。字符串关键字带行内注释会把注释吞进值里
导致匹配失败；数值关键字带行内注释安全（含 UTF-8 中文注释）。

## TEST PES 与 loud-stop 桩

`TEST_PES=HARMONIC|MORSE|LEPS` 运行时选择；参数 `HARM_K`/`MORSE_DE/RE/A`/
`LEPS_DE/RE/A/DELTA`，缺省取内置 H₂/H₃ 值。不读任何外部势文件。

本构建的 loud-stop 桩（调用即打印消息并 STOP，不静默出错物理）：
`ELEC_METHOD=TDHF-FSSH/IESH/MDEF`、`NZDOWN=1`、`TASK=NORMAL-MODE`（F14）。

## Known Defects（登记，不修）

| 编号 | 缺陷 | 触发 |
|---|---|---|
| F14 | NSELT=-1 死代码入口 | `TASK=NORMAL-MODE` → loud STOP |
| F15 | RANDST 双调用 | NSELT=2/3 用户 ISEED 无效（序列由种子 1 决定） |
| F17 | NSELT=3+MB 曾静默冻结 | 已修为 loud STOP |
| F19 | 振动角动量剥离 | 微正则简正模：呼吸模能量 → E′ 对（0.20 vs 0.33） |
| F20 | LOCAL-MODE 参数无来源 | `INIT_SAMPLING_*=LOCAL-MODE` → NaN+挂起 |
| F22 | 共线伪模 NSTEP 溢出 | 共线构型 + FIXED-ENERGY（NBAR=3）→ 停机 |
| F23 | CI-QM-MICRO 死循环 | `INIT_SAMPLING_*=CI-QM-MICRO` |
| F24/D1-D4 | 刚性表面四重缺陷 | `SURFACE_MODEL=RIGID` → 越界读 → 全轨迹 NaN；旧 `NSURF=2` 重映射为 RELAXED（D1） |
| — | RK4 不可达 / ADAMS 损坏 | 有效 1 阶、漂移最差；禁用 |

缺陷签名的可复现证据：`cases/rigid_surface_defect`、`cases/microcanonical_normalmodes`、
`cases/integrator_matrix`。

## cases/ 约定

每 case 一个语义名目录：`README.md`（验证目的/原理/预期与结果三节）+
`input_qct.txt`（可直接运行）+ `results.txt`（程序输出摘录）+ 个别坐标/变体文件。
新增 case 遵循同一结构；统计类断言写明分布、门限与样本量；缺陷登记类 case
的「预期」即损坏签名。

## Key Conventions

- 固定格式 `.f` 与自由格式 `.f90` 混排；`IMPLICIT DOUBLE PRECISION (A-H,O-Z)`。
- 全局状态全在 `venus_data` 模块；切换分支或改模块文件后删除陈旧 `.mod`。
- 旧 `SIZES` 文件已被 `venus_params.f90` 取代（本仓库不含）。
- NCHKP checkpoint 功能已裁剪（死关键字，输入无副作用）。
