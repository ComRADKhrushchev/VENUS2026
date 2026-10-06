# baseline_kcal_v1 —— 旧单位制（kcal/mol 门面 / 10fs 时间基）值域基线

**用途**：单位制改造（eV–Å–amu–fs，2026-10-06 裁决）前的对照基线。
改造后的验收判据（守恒 / 初帧能量锚 / 分布比对）以本目录为旧制参照；
逐字节等价门自本次改造起废止（详见 src/EXTENDING.md 演练清单）。

## 生成条件（出处）

- 代码状态：commit `0b9faaf`（工作树干净，无未提交改动）
- 构建：`cd src && make clean && make`（gfortran，单容器两轮装配，
  容器按 systems/inactive/README.md 拖入程序激活：
  h3 两 case 配 h3_bkmp2，h2ag 两 case 配 h2ag111PES，测毕移回）
- 运行：各 case 目录内 `../../../src/build/gfortran/venus.e`，
  运行日志 `logs/p0_run_<case>.log`
- 日期：2026-10-06

## 各 case 实际输入

- h3_ebk：`input_used.txt` = `input_manual_qct.txt` 变体
  （EBK 双原子常数通道需 SPECTRUM_SOURCE=MANUAL + DIAT 记录；
  仓库默认 `input_qct.txt` 的 COMPUTE 模式对双原子常数无通用推导，会 STOP 1）
- h3_wigner / h2ag_ebk / h2ag_wigner：`input_used.txt` = 该 case 的 `input_qct.txt` 原样

## 产物单位（旧制）

- trajectories.txt：t 列 [10 fs]（内部计数原样），能量列 [kcal/mol]，
  q [Å]，p [amu·Å/(10 fs)]；终态行 = traj, t, t_kin, v_pot, e_tot
- hist_e_rel.dat / hist_thta.dat：横轴物理量同上单位约定
