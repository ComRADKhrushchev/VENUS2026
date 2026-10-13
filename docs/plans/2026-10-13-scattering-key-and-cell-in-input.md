# 2026-10-13 散射类型主控键 + 晶胞入输入 + A_LAT 退役 —— 裁决记录

## 裁决（用户，2026-10-13）

1. **A_LAT 键移除**：依托 a_lat 的格点统计（top/bridge/hcp/fcc 落位分组直方图、
   位置图折叠）是对特定（刚性、规则）表面的设计；对一般表面——特别是可动
   slab——完全多余。此类统计交给后处理程序。
2. **散射类型主控键**：范式（表面散射 / 气相散射）不再由构型文件是否带胞推导
   （.poscar ⇒ 表面的旧规则退役），改由输入键决定。
3. **晶胞入输入**：体系晶胞写进 input_qct.txt，构型文件只携带几何（.poscar 的
   胞仅作该文件 Direct 坐标的参考架；单帧规则保留：至多一个带胞文件）。

## 落地语义

- `SCATTERING`：枚举 `GAS`（默认）/ `SURFACE`，块二；槽写
  `reactants%surface_model`（GAS→0，SURFACE→1）。
- `CELL`：条件键（kind_cell），值为三行晶矢（三逗号组各三实数，行存 [Å]）；
  槽写 `reactants%cell` + `cell_present`。装配期窗口（config_atoms
  window_reconcile，先于盒派生）：SURFACE 缺 CELL → 指名中止；GAS 写 CELL →
  指名中止。表面盒（box_lx/box_ly/skew，余弦定理）由 CELL 派生——**值权威**，
  不与文件胞交叉校验。N_AIM=1 落点域与 AIM_X/Y 坐标架随之取 CELL 派生盒。
- `A_LAT`：键、config/list_atoms 字段、inc_surface 公度门、plot_hist 的
  k_pos/k_site/cell_fold/lat_t 格点统计机械整体删除（k_bnc 反弹分组保留，
  直方图类收缩为 1/2/3）。驱动侧 lattice 交接（driver.f90）同删。
- sysdef：删除范式派生（merge(1,0,any_cell)）与 buffer_cell_present/
  buffer_unit_cell 访问器；文件胞不再进入体系层。
- 默认 GAS 保持未写键的存档流不变（旧分子案例零扰动）。

## 改动面

src：input.f90（键表 27→28：+SCATTERING/+CELL、−A_LAT；kind_posreal→kind_cell；
parse_cell）、config.f90（reactants%cell/cell_present、surface_model 注释）、
sysdef_interface.f90（去派生/访问器）、config_atoms.f90（window_reconcile 重建、
盒源切 CELL、a_lat 字段删）、inc_surface.f90（公度门删）、plot_hist.f90 +
driver.f90（格点统计删）。
判据：check_input.f90（(o) 窗判据重写为 SCATTERING/CELL、v-cell 载体、计数）、
check_roster.f90（(b)/(e)/(f) 改 CELL 驱动、子链补 read_input）、
check_inc_surface.f90（(k)/(l) 腿随公度门退役；三处速度断言单位迁移
e_conv→dyn_ev——该判据是 2026-10-06 改基的漏网件，另发现 scripts/ 下陈旧
.mod 使其对旧制常数静默编译，已清除）。
案例：surface_window 与 spectrum_source 的 h2ag_ebk/h2ag_wigner 输入改写
（SCATTERING/CELL 取代 A_LAT；CELL 值 = 原 POSCAR 胞原值）。
GUI：keyspec（+SCATTERING/+CELL、−A_LAT、KIND_CELL 解析）、validate（范式改
读 SCATTERING；A_LAT 规则→CELL 必填/禁写；公度规则删）、sysdir（surface_model
推导删，扫描胞降为信息性）、preview（胞/落点框取 CELL 键）、input_page
（SCATTERING 下拉+CELL 输入；范式呈现随主控键刷新）。

## 验证（读数落盘 /tmp/paradigm_cell.log 与 gui/logs/）

- gfortran + ifx 双编译器构建绿；
- check_input / check_roster / check_inc_surface：0 失败；
- GUI 套件：93 passed 2 deselected；
- 数值中立性：case_surface_window 重跑轨迹与改前**逐字节一致**；h2ag_ebk/
  h2ag_wigner 活置 h2ag111PES（容器架规程，跑后还原）重跑轨迹与
  baseline_ev_v2 **逐字节一致**、对 baseline_kcal_v1 值域门 ALL GATES PASS；
  气相抽查 case_harmonic_scatter 逐字节一致。
