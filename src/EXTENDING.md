# EXTENDING.md — 新程序扩展点(E27 后契约)

每类扩展 = **新文件 + 薄 reg 文件**,拖入即编译(Makefile wildcard + prebuild 聚合),
**零触碰现有文件**(T2 演练在案:docs/plans/S.md 验收区)。登记迟滞:新增 reg 后第二次
make 即全绿(depend.mk 单轮迟滞,Makefile 头注在案)。

## 1. 体系容器(systems/<system>/)——组合根

一个文件夹 = 一个体系,自足材料:

```
systems/mysystem/
  my_pes.f90        体系物理: H 定义/势能/力律/随机律, eV 交付
  my_box.f90        (可选) 参数文件通道 + 终止/分类判断
  reg_mysystem.f90  拖入件: reg_mysystem_all() 绑 container_bind_pes(my_pes)
                    [+ container_bind_term / container_bind_classify]
                    [可选 _init(): 参数装载缝, 与 reg_dist_* 同形]
```

- **力交付契约(E27)**:`subroutine my_pes(q, mass, g, v)`——q[A]/mass[amu] 进,
  组合好的活性力核 `g = dV/dq`[eV/A] 与势能 `v`[eV] **全赋值**出;换算归门面,勿自换。
  多面体系自折活性面(对角化+占据折叠;通用助手登记为候选后续件,git 历史可取
  force_facade 退役前实现)。
- 随机律成员在容器内部和(g_add 语义已内化),rng 单流直用。
- 体系专属**电子方法包**(如自定义 Hubbard 传播)同_folder 放置:prop/match 成员 +
  `reg_*` 里一行 `method_reg('my_tsh', my_prop, my_match)`;成员与 H 模块共享私有数据
  (E27:无互换工作区);输入键 `ELEC_METHOD = my_tsh` 选用(绝热参考跑仍是
  `ELEC_METHOD = adiabatic` 或缺省)。
- 参数经运行目录参数文件(参照 systems/reference/reference_box.f90 的通道形态)。

## 2. 初始态分布成员(methods/dist_*.f90,B′ 族)

新成员 = `methods/dist_mydist.f90`(init/sample + 具名停机族)+ 薄
`methods/reg_dist_mydist.f90`(`samp_reg` 一行 + `_init` 装配缝消费成员键);
`sampler` 的 scheme 码表加一词(scheme_code_name——这是族表数据非逻辑分支)。
零体系数据:质量/几何取 roster,参数走输入键。

## 3. 通用方法包与组件(methods/,库章程)

- 通用包:elec_*.f90 + reg(行进方法表,任何体系可经 ELEC_METHOD 选用)。
- 组件(level 2):amp_*(振幅积分器)/ fssh_g_ij / hop_decide / momentum_rescale /
  dom_reset 计划位——容器或成员经 `use` 组合,物理单拷贝。

## 4. 积分器(闭族,methods/reg_prop.f90)

verlet/symple/radau 固定成员制;新增积分器 = dynamics/propagator.f90 内实现 +
reg_prop 注册行(闭族非开放注册面——成员资格是裁决事项,见 E 线裁决记录)。

## 5. 记录列与后处理

- 记录列:output/recorder.f90 列注册表加行(rec_reg_col;列提取器纯形式)。
- 结果分类:容器侧 container_bind_classify(E19 槽)。
- 直方图类:output/plot_hist.f90 五类族(数据驱动)。

## 演练清单(新扩展件落地后)

1. `cd src && make`(如新增 reg:再 make 一次)→ 双编译器绿;
2. 相关 scripts/check_* 重编重跑;
3. 五 case 回归对 results.txt(逐字节);
4. T4/T5 grep 复验(23.0605 仅 force_facade;RANDOM_NUMBER 零命中)。
