# VENUS2026 — QCT 化学动力学程序（最小可编译核心）

绝热 QCT（准经典轨迹）化学动力学程序的完成式重写：L0–L5 分层的 `src/` 新树。
本仓库当前只发布**最小可编译核心**（`src/` 全部源码 + `scripts/prebuild.sh` 注册聚合脚本）；
案例（cases/）、判据程序、研究文档（docs/）与物理系统容器材料（systems/）暂缓上传。

## 构建

```bash
cd src
make              # gfortran（默认）-> build/gfortran/venus.e
make FC=ifx       # Intel ifx      -> build/ifx/venus.e
make clean        # 删除整个 build/ 子树
```

- 单一入口 Makefile：wildcard 源码拾取（`src/*.f90`、`src/*/*.f90`、`systems/*/*.f90`）
  + 每次构建运行 `scripts/prebuild.sh` 生成 `build/<FC>/registry_gen.f90`（聚合注册模块——
  注册是构建产物，driver 零手写注册行）与 `depend.mk`。
- 无容器拖入时同样可编译可运行（退化绝热参考面仍绑定）。

## 架构速览

```
src/
  driver.f90 / input.f90      L0 控制内核：装配 → 演化 → 统计，零物理
  globals/                    L1 变量与供给（state/control/config/config_atoms/consts）
  interface/                  L2 接口与框架（force/elec/sysdef/sampler/container_sched/bath）
  dynamics/                   L3 闭族积分器（verlet · symple · radau，固定成员制）
  methods/                    L3 方法库（elec_* / samp_* / inc_* / beam_laws + reg_*）
  output/                     L4 记录与后处理（recorder/final_state/density/plot_hist）
  utils/                      L5 仪器（rng 单流/linalg/geometry/specio/hessian）
```

详见 `src/README.md`（分层规则）、`src/EXTENDING.md`（扩展点）、`src/CLAUDE.md`（构建/判据要点）。
