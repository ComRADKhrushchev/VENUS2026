# VENUS2026 — QCT 化学动力学程序（最小可编译核心）

绝热 QCT（准经典轨迹）化学动力学程序的完成式重写：L0–L5 分层的 `src/` 新树。
本仓库当前发布**最小可编译核心 + 两个已验证体系的容器 + 旧版 VENUS05 标定树**；
案例（cases/）、判据程序、研究文档（docs/）与其余容器材料暂缓上传。

## 已验证体系的容器（systems/inactive/）

两套物理系统容器已完成对新程序与旧引擎的**逐位交叉验证**（能量与梯度多几何点
比对，全部 BITWISE IDENTICAL），平时收纳于 `systems/inactive/`（不进构建）：

- `h3_bkmp2/` — H+H2 / BKMP2 势能面容器（f90 转录 `bkmp2_pot.f90` 与
  VENUS05 官方引擎 12 几何点逐位一致；含 `system/` 反应物构型 H + H2）
- `h2ag111PES/` — H2/Ag(111) 神经网络势容器（6 表面几何点逐位一致；
  运行需 `weights-h2ag.txt` / `biases-h2ag.txt` 随运行目录，已随容器携带）

激活方式：把容器文件夹移到 `systems/` 顶层 → `make` 两轮（depend.mk 单轮迟滞）。
注意 one-container-per-process 立法：目标容器与其他容器并存时装配期 double-bind
即停；一次只激活一个。

## 旧版 VENUS05 标定树（benchmark_oldvenus/）

从 `VENUS_TSH_HO2` 剥离 EANN/NN、TSH、AmberFF 与表面支路后的干净 VENUS05
血统树，即上述两体系交叉验证的旧引擎一侧（`sys/bkmp2_engine.f`、`sys/h2ag_engine.f`）：

```bash
cd benchmark_oldvenus && make    # gfortran -std=legacy -> venus.e
```

`work/`、`work_h2ag/` 为冒烟运行目录（输入 + 通道文件）；134MB 的扫能运行数据
（h2ag_sweep/）不入库。详见其 `README.md`。


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
