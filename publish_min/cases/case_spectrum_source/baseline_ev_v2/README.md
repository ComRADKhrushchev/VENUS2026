# baseline_ev_v2 —— eV/Å/amu/fs 制基线（2026-10-06 单位改制定版产物）

**生成条件**：单位改制工作树（C1 提交态），gfortran 构建，容器按
systems/inactive/ 拖入程序逐个激活（h3_bkmp2 / h2ag111PES），各 case 目录内
`../../../src/build/gfortran/venus.e` 运行。判据门 4/4 全绿
（logs/check_ev_regression_gfortran.log（副本见本目录 logs/）：A2 初帧锚 ≤6.7e-6、A3 结构、A4 时间轴）。

**产物单位（新制）**：t 列 [fs]，能量列 [eV]，q [Å]，p [amu·Å/fs]；
终态行 = traj, t, t_kin, v_pot, e_tot。

**与 baseline_kcal_v1 的关系**：v1 是改制前旧制产物（kcal/mol、10fs 计数），
作为值域判据的旧制参照永久保留；v2 是改制后的现行基线。后续回归对 v2 判据
（A2 锚对 v1 的 ÷23.0605 比对仅改制验收用一次）。
