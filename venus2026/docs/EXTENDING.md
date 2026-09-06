# 扩展指南：接入自定义势能面与动力学方法

本指南面向要在 VENUS2026 上接入自己的势能面（PES）、初始条件采样方法或
积分器的开发者。发布仓库自带一套体系无关的 TEST 势（HARMONIC / MORSE / LEPS），
它同时就是**接入自定义组件的活模板**——本文所有示例均指向仓库内可直接阅读
的源码。

## 1. 扩展点地图

一条轨迹的生命周期与扩展点位置：

```
input_qct.txt
   │  venus_input.f90（关键字解析；解析完成后 CALL POTPRE ←── 【PES 初始化钩子】
   ▼
SELECT.f90（初始条件采样，按 NACTA/NACTB 分派 ←──────────【采样方法扩展点】
   │   各采样例程：INITEBK / INITQP / ORTHAN / THRMAN / QMMICRO / …
   ▼
VENUS.f90 主循环（按 INTEGRATOR 分派 ←───────────────────【积分器扩展点】
   │   INTEGRATOR=1 → RADAU    INTEGRATOR=2 → SYMPLE    INTEGRATOR=3 → VERLET
   │   每步调 DVDQ.f → CALL DPESHON ←───────────────────【PES 梯度入口】
   │   采样/输出时调 ENERGY.f（ENERGY_1）→ CALL POT0 ←───【PES 能量入口】
   ▼
fort.* 输出（相空间）+ stdout（轨迹统计）
```

| 想接入什么 | 改哪里 |
|---|---|
| 自定义势能面 | 实现 `POTPRE` / `POT0` / `DPESHON` 三符号契约（§2） |
| 自定义初始条件采样 | `param_mapping.f90` 的 `map_init_sampling` + `SELECT.f90` 分派（§3） |
| 自定义积分器 | `VENUS.f90` 主循环 + `param_mapping.f90` 的 `map_integrator`（§4.1） |
| 非绝热电子方法 | `map_elec_method` 的 loud-stop 桩即预留扩展点（§4.2） |

## 2. 接入自定义势能面

### 2.1 三符号契约

PES 以三个 Fortran 外部子程序接入核心（样板：[src_TEST/interface_TEST.f90](../src_TEST/interface_TEST.f90)）：

| 符号 | 签名要点 | 调用方 / 时机 |
|---|---|---|
| `POTPRE()` | 无参 | `venus_input.f90`——输入解析完成后调用一次，做势初始化（读参数、打印头部） |
| `POT0(NDUM, Vpot, Qarr)` | 返回 `Vpot` [kcal/mol] | `ENERGY.f`（ENERGY_1）、`VENUS.f90`——能量求值 |
| `DPESHON(NDUM, Qarr)` | 把梯度累加进全局 `PDOT` | `DVDQ.f`——每积分步的力 |
| `FIXROTDATM(I)` | 表面键旋转 | 仅 `NZDOWN=1` 表面构建需要；TEST 构建为 loud-stop 桩 |
| `GASDEV()` | 正态随机数 [函数] | 恒温浴 THERMO/GLO 调用；Box-Muller 实现在 interface 层（interface_TEST.f90 内），`src_VENUS/GASDEV.f` 另提供等价的 `grandom()` |

两个约定务必注意：

1. **POT0 / DPESHON 都在 live 全局坐标 `Q` 上求值**（`venus_data` 模块），
   哑元 `Qarr` 只是接口占位——照抄 interface_TEST 的写法即可；
2. **DPESHON 是累加**（`PDOT(i) = PDOT(i) - g(i)*23.0605*C1`），不是赋值。

### 2.2 单位契约（最容易出错的地方）

核心的内部单位制（详见 [CLAUDE.md](../CLAUDE.md)）：质量 amu、时间步 DT 以
10 fs 计、动量 amu·Å/(10 fs)、能量 code 单位 = amu·Å²/(10 fs)²。
PES 接口层的换算责任划分：

| 你的势例程 | 收到 | 必须返回 | 换算（在 interface 层做） |
|---|---|---|---|
| `v(q)` | 坐标 | **eV** | `POT0`：`Vpot = v_ev * 23.0605` → kcal/mol |
| `vg(q)` | 坐标 | **eV/坐标单位** 的 v 与梯度 | `DPESHON`：`PDOT -= g * 23.0605 * C1` |

常数：`23.0605`（eV→kcal/mol）、`C1 = 0.04184`（kcal/mol→code 能量），均在
`venus_params`。

**坐标语义由你的势定义**。TEST 势把 Q 的数值按 bohr 解释（`m_re=1.401` bohr 等
参数与 `QZA_EQ` 输入自洽）；换成别的构建时同样只需保持"势参数、几何输入、梯度
量纲"三者自洽。判断标准：`cases/` 任一 case 的 fort 输出回读能量应与 stdout
打印一致（对照法见 `cases/leps_potential/results.txt` 的首末步摘录）。

### 2.3 实现步骤（照抄 TEST 势样板）

1. **写势模块**（模板 [src_TEST/test_potentials.f90](../src_TEST/test_potentials.f90)）：
   - 只 `public` 两个求值例程 `v` 与 `vg`（分离便于只求能量时省梯度）；
   - 势参数经输入关键字读取（`get_real('MORSE_DE', m_de)` 模式），缺省值内置；
   - 解析锚写进注释（TEST 势的 Morse 能级公式即例）——供将来写验证 case 用。
2. **写 interface 三符号**（模板 [src_TEST/interface_TEST.f90](../src_TEST/interface_TEST.f90)，
   通常 <50 行）。替换路径：直接修改 interface_TEST 并换 `use` 你的模块；
   并存路径：新建 `interface_X.f90` 并在 Makefile 中替换 `src_TEST/interface_TEST.f90`
   一行。
3. **登记 Makefile**：`SRC` 列表按模块依赖顺序单次编译——你的势模块放在
   `venus_input.f90` **之前**（POTPRE 要被它调用），interface 文件放它**之后**。
4. **重建并验证**：`make && cd cases/<name> && ../../venus_test.e`。
   建议先复制 `cases/leps_potential` 做最小冒烟（2 原子、5000 步、秒级）。

### 2.4 PES 验证建议

- **解析锚优先**：能写出闭式解的量（频率、能级、鞍点）先对（TEST 势每个 case
  的「预期与结果」小节全是这种锚）；
- **守恒验证用 fort 重算**，不要信能量打印（原因见 §5.1）：
  `cases/integrator_matrix` 给出了从 fort.1001 坐标/动量重算 E 的方法与
  10 ps 漂移门限的先例；
- **逐位确定性**：固定 `ISEED` 复跑 diff 全部产物，是回归验证的最强手段。

## 3. 接入自定义初始条件采样

关键字流程：`INIT_SAMPLING_A=XXX`（input_qct.txt）→
[param_mapping.f90](../param_mapping.f90) `map_init_sampling` 映射为整数
`NACTA` → [src_VENUS/SELECT.f90](../src_VENUS/SELECT.f90) 按 `NACTA` 分派到
采样例程（`NACTB`/片段 B 同理）。

现役映射（`map_init_sampling`）：`MB`=0、`ORTHANT`=1、`MICROCANONICAL`=2、
`NORMAL-MODE`=3、`LOCAL-MODE`=4、`BOLTZMANN-VIB`=5、`FIXED-ENERGY`=6、
`MD`=7、`QM-MICRO`=8、`CI-QM-MICRO`=9。

接入两条路径：

- **替换既有分支**：直接改对应例程（如 `INITQP.f90`、`ORTHAN.f90`）。
  这些例程均为现代化后的风格（implicit none、自由格式、无 GOTO），
  本身就是新代码的风格样板；
- **新增 NACTA 值**：`map_init_sampling` 加一行映射 → `SELECT.f90` 加分派
  分支 → 新例程文件登记进 Makefile `SRC`。

采样代码的三条铁律（违反即静默错误物理）：

1. 采样能量是 kcal/mol，**乘 `C1` 后**才能算振幅/动量；
2. `P²/(2m)` 给的是 code 能量，与 kcal/mol 比较前**先除 `C1`**；
3. 每轨迹动量组装完必须去质心动量（`CENMAS`），`cases/` 全部采样 case 的
   「预期与结果」都带 max|Σp| 检验。

## 4. 接入自定义积分器与电子方法

### 4.1 积分器

关键字流程：`INTEGRATOR=XXX` → `map_integrator` → 整数 `INTEGRATOR` +
子模式 `LLL` → `VENUS.f90` 主循环分派（1=RADAU、2=SYMPLE、3=VERLET 族）。

新增积分器三步：新例程（模板 `src_VENUS/VERLET.f90` / `RADAU.f` / `SYMPLE.f`）
→ `VENUS.f90` 主循环加 `ELSEIF (INTEGRATOR.EQ.N)` 分支 → `map_integrator`
加映射。注意：

- 主循环在 `MOD(NC,NIP).EQ.0` 时调 `GWRITE` 打印/落盘——新积分器推进 `NC`
  的方式必须与该门兼容；
- 阶数与漂移验证照抄 `cases/integrator_matrix`（DT 扫描 log-log 斜率 +
  fort 重算 10 ps 漂移）；历史上 ADAMS 例程因名义 6 阶实际 1 阶被移除
  （该 case 的漂移柱状图 `fig_drift_matrix.png` 即验收证据形态）。

### 4.2 非绝热电子方法（预留扩展点）

本构建只含经典动力学：`ELEC_METHOD` 只接受 `ADIABATIC`；
`TDHF-FSSH` / `IESH` / `MDEF` 在 `map_elec_method`（param_mapping.f90）
**loud-stop**——打印明确错误并停机，不静默回退到经典。

这是刻意的裁剪（完整非绝热实现见内部开发仓库）。在此接入非绝热方法需要：

1. `map_elec_method` 放行你的方法名并映射 `CALTYP`；
2. 势接口从三符号契约扩展为**多态势能**（`POT0` 需返回所选态能量与跃迁
   驱动，`DPESHON` 需按当前活性态给力）；
3. `VENUS.f90` 主循环加分支决策（每步或每次跃迁后重置力缓存）。

## 5. 验证注意事项

### 5.1 GWRITE 能量打印是陈旧值（重要）

`VENUS.f90` 的 VERLET/SYMPLE 分支积分后直接调 `GWRITE`、**不调 `ENERGY_1`**，
因此：

- fort 头部的 `E0/T/H(eV)` 与 stdout 的逐行 `KINETIC/POTENTIAL/TOTAL ENERGY`
  在整个轨迹期**恒为采样末次 `ENERGY_1` 的值**；
- 打印恒定**不构成能量守恒的证据**（直接证据：
  `cases/leps_potential/results.txt` 首末步摘录——坐标已演化、H(eV) 纹丝不动）；
- 守恒验证一律从 fort 相空间重算（方法与门限先例：
  `cases/integrator_matrix`）。

采样阶段的能量回显（`INTERNAL ENERGY =`、`CHOSEN: EROTA/EVIBA`、
`system temperature=`）是**活数据**——采样缩放循环以 `ENERGY_1` 判敛，
这些行每个新 case 都能直接用（各 case 验证图的形态见
`fig_sampling_energy.png` 系列）。

### 5.2 loud-stop 哲学

不支持的关键字组合一律打印原因并 `STOP`（`map_elec_method`、`FIXROTDATM`
桩、`TASK=NORMAL-MODE` 均如此）。接入自定义组件时保持这一约定：**宁可 loudly
失败，不可静默出错物理**。

### 5.3 为新组件建 case

每个新组件配一个语义名目录 `cases/<name>/`：`README.md`
（验证目的/原理/预期与结果三节，统计断言写明分布、门限与样本量）+
`input_qct.txt`（可直接运行）+ `results.txt`（输出摘录）+ 验证图 +
轨迹类 case 入库 10 条 `fort.1001`-`fort.1010` 供检验。全部 16 个现有
case 均为此形态，任取一个照抄即可。

## 6. 快速清单

接入自定义 PES：

- [ ] 势模块：public `v`/`vg`，参数走输入关键字，解析锚写注释
- [ ] interface 三符号 `POTPRE`/`POT0`/`DPESHON`（eV→kcal/mol→code 两级换算）
- [ ] Makefile `SRC` 按依赖顺序登记（势模块在 venus_input 前、interface 在后）
- [ ] `make` 通过；冒烟 case 复跑逐位一致
- [ ] 验证 case：解析锚 + fort 重算守恒 + 图表

接入采样 / 积分器 / 电子方法：

- [ ] `param_mapping.f90` 加映射（非法值 loud-stop）
- [ ] `SELECT.f90` / `VENUS.f90` 加分派分支
- [ ] 新例程 implicit none、自由格式、无 GOTO（样板：`INITQP.f90`）
- [ ] 单位铁律：kcal/mol ×C1 才算动量；比较前 ÷C1；去质心动量
- [ ] 验证 case：分布 KS/χ² 断言 + 逐位确定性
