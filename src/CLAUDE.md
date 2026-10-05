# CLAUDE.md — src/ 新树工作指导

本件承接旧树 `venus2026/CLAUDE.md` 的角色(该件继续描述旧参考树,只读)。工作区级通用
指导见根 `AGENTS.md`;架构权威见 `docs/design_spec_completeness.md` 与 `src/README.md`。

## 构建与运行

```bash
cd src
make              # gfortran(缺省)→ build/gfortran/venus.e
make FC=ifx       # ifx(项目约定旗标 -r8 -double-size=64 -i8 -O2 -w)→ build/ifx/venus.e
make clean        # 删整个 build/(两编译器子树)
```

- 单一入口 Makefile;wildcard 拖入 + `scripts/prebuild.sh` 每 make 再生
  `build/<FC>/registry_gen.f90`(聚合注册)与 `depend.mk`(对象级依赖)。
  **新增 reg 文件后需第二次 make**(登记在案的残余限制:depend.mk 单轮迟滞)。
- 产物只落 `build/<FC>/`;源树绝不留 .mod/.obj(K2 陈旧 mod 遮蔽教训——"Unexpected EOF")。
- Windows ifx 的链接器遮蔽/模块目录问题已在 Makefile 平台段处理(MSVC link.exe 路径、
  `-module:` 冒号形态),无需手工干预。
- 运行:`cd cases/<name> && ../../src/build/gfortran/venus.e`——驱动器硬读 cwd 的
  `input_qct.txt`;产物 `trajectories.txt` / `hist_*.dat`(gitignored)。

## 判据体系(scripts/)

`scripts/check_*.f90` 为组件级判据程序(头注载编译配方与判据先行声明),探针
(`probe_glo.exe` / `probe_chain_freq.exe`)链接 build/gfortran 对象集(除 driver.o)。
改接口/成员后:相关 check 重编重跑 + 双编译器重建 + 五 case 回归对 `results.txt`。

## 契约要点(改代码前必读)

- **E27 物理系统容器**:物理系统容器(systems/<sys>/)交付组合完成的受力计算
  `container_pes_i(q,mass,g,v)`(eV,全赋值);力接口唯一换算点 `ev_kcal*e_conv`(T4);
  ws 互换已退役(决定记录:docs/plans/S.md E27)。
- 方法包(电子方法)可居 methods/(通用)或物理系统容器文件夹(体系专属,与 H 定义共享私有模块);
  注册一律 `method_reg(name, prop, match [, sync])`,装配 `method_bind` 纯名字查找。
- 成员统计写回:`occ_set(occ)` / `n_hop_set(n)`;读取 getter 供 recorder。
- 内部单位 Å/amu/(10 fs);物理常数定于 `globals/consts.f90`(取值按既有决定固定),
  解析侧须用固定换算对(23.0605×0.04184=0.96485132,V-4 在案),不用 CODATA 重推。
- 单随机流(rng.f90);演化阶段随机消费只在 match 成员。
- 演化循环五入口固定:prop_step → elec_prop → mqc_statistic → rec_frame → container_term。
- 输入 KEYWORD=VALUE 独立语法(E21/E22);成员键经暂存行由各 reg 文件的 `_init`
  初始化组装机构消费。
