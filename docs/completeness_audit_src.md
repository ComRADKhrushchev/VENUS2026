# VENUS2026 新程序(src/)例程完成度复审报告

- 日期: 2026-09-16(S 相 S2/T7;E27 容器组合根重裁当日——本审计覆盖 E27 后的树)
- 裁决标准: `docs/design_spec_completeness.md`(物自体/表象二分、四类纯子程序、控制内核、
  接口 A/B′、闭合选择族、三问复审)
- 审计方法: 只读静态复审。全树 50 源文件(driver/globals/facade/dynamics/methods/output/utils/
  systems)逐文件头注与例程清单核对;三问之(c)(内含体系分支或常数)另做机械化全树 grep
  (元素名/硬编码原子序号模式/旧库魔数 2.94998/413.4/1D30——**全树零命中**);
  动态证据援引 V 相五 case 与 scripts/ 判据程序(check_* 全绿记录)。
- 机器可读主表: `docs/completeness_audit_src.json`(与本 md 同源生成)。
- 结论词汇: 每例程三问各 ✓/✗;总体 **PASS**(零红项)/ **PASS-注**(带登记注记)。
- 类别词汇: 内核(L0)/ 供给(L1)/ 门面与框架(L2)/ 闭族成员(L3 传播器)/ 方法库成员(L3 dist_*)/
  电子方法包(L3 elec_*)/ 注册材料(reg_*)/ 缝工具(seam_cast)/ 测量(L4 记录)/ 后处理(L4)/
  仪器(L5)/ 体系材料(systems/——体系知识是其本职,三问之(c)按定义豁免,审其零程序知识)。

---

## 0. 结论先行

| 审计面 | 结论 |
|---|---|
| 三问 (c) 机械 grep(元素名/原子序号/旧库魔数) | **全树零命中**——通用代码零体系知识渗漏 |
| 三问 (b)(只用下层、按契约请求服务) | 全树依赖单向(成员→utils,门面→L1,驱动→门面);分层 grep 无上溯(use 面经 depend.mk 全量核对) |
| 体系材料隔离 | systems/reference/ 自足(pure-math 头注在案),经薄 reg 进装配,零程序知识 |
| 红项 | **0**;登记注记 7 条(下表 §2),全部有在案登记与去向 |
| 与旧审计(docs/completeness_audit.md,99 条)的关系 | 旧审计为旧树(venus2026/)历史参考,保留不动;本报告为 src/ 新树的 T7 交付 |

## 1. 主表(逐文件逐例程)

三问列 = (a) 只处理 Q/P/s/F 或己方供给物 ✓ (b) 按契约请求服务 ✓ (c) 内含体系分支/常数 ✓;
"✓✓✓"为全过;注记以 §2 编号引用。

### 1.1 内核与入口(L0)

| 文件 | 例程 | 类 | (a)(b)(c) | 注 |
|---|---|---|---|---|
| src/driver.f90 | venus | 内核 | ✓✓✓ | 装配→演化→统计三段纯编排;演化环内仅 facade/propagator/recorder 五入口(T6);登记注 N1/N2/N3 |
| src/input.f90 | read_input, parse_file, collect_input, check_required, inject_defaults, legal_words, line_tag, map_slots, parse_elist, parse_enum, parse_int, parse_level, parse_log, parse_posreal, parse_real, reset_staging, staging_key, staging_val, i2s | 供给(根入口) | ✓✓✓ | 独立语法无回退轨(E21);成员键暂存原始串供装配缝;两块槽映射;零物理 |

### 1.2 供给层(L1)

| 文件 | 例程 | 类 | (a)(b)(c) | 注 |
|---|---|---|---|---|
| src/globals/state.f90 | state_create, state_destroy | 供给 | ✓✓✓ | 纯 Q/P/s/F/t 容器,量纲按 roster 内容参数推导 |
| src/globals/config.f90 | (模块单例: propagator/electronic/reactants/observables 四组) | 供给 | ✓✓✓ | 表象层配置,零物理 |
| src/globals/control.f90 | (模块单例: n_traj/i_seed/max_steps) | 供给 | ✓✓✓ | 编排变量,零物理 |
| src/globals/consts.f90 | (立法常数参数) | 供给 | ✓✓✓ | C1=0.04184 等唯一发生点;换算消费仅 facade + I/O 缝(T4 在案) |
| src/globals/config_roster.f90 | roster_load, roster_clear, reset_derived, window_reconcile, stop_assembly | 供给 | ✓✓✓ | roster 单向取自 sysdef 暂存(装配序法);Q-TERM/E_SPLIT 表面窗键 reconcile(E23) |

### 1.3 门面与框架(L2)

| 文件 | 例程 | 类 | (a)(b)(c) | 注 |
|---|---|---|---|---|
| src/facade/force_facade.f90 | force_eval | 门面 | ✓✓✓ | E27 组合核交付契约;换算唯一点 ev_kcal·e_conv(f90:30 定义/fold 点) |
| | container_pes_i(接口), container_bind_pes, container_bind_term, container_term, container_init, stop_ff | 门面 | ✓✓✓ | 纯形式权威(E13/E27):签名/单位/槽位/名字,无物理合法性表 |
| src/facade/elec_facade.f90 | elec_prop, mqc_statistic | 门面 | ✓✓✓ | 四入口之二,表查派发 + sync 尾 |
| | method_reg, method_bind | 门面 | ✓✓✓ | 单表纯名字查找;行可来自 methods/ 或 systems/(E27 组合根) |
| | occ_get, occ_set, n_hop_get, n_hop_set | 门面 | ✓✓✓ | 运行期统计持有者读写对(E27 新增写回入口,形式校验) |
| | elec_prop_i, mqc_match_i, sync_proc_i(接口), stop_ele | 门面 | ✓✓✓ | 成员契约文档位,conformance 在注册点机检 |
| src/facade/sampler.f90 | samp_run | 门面 | ✓✓✓ | 逐片段派发 + incident 编排位后置一次(单片段气相跳过) |
| | samp_reg, samp_code, samp_n_dist, dist_proc_i(接口), find_row, i8_str, stop_samp | 门面 | ✓✓✓ | B′ 成员表;E13 纯形式 |
| src/facade/sysdef_facade.f90 | sysdef_load, list_folder, classify_entry, collect_system, reset_staging, staging_*(6), | 门面 | ✓✓✓ | 体系定义发现/解析派发/质量表/暂存;一文件一反应物;表面范式派生 |

### 1.4 闭合族(L3 传播器)

| 文件 | 例程 | 类 | (a)(b)(c) | 注 |
|---|---|---|---|---|
| src/dynamics/propagator.f90 | prop_step, prop_reg | 闭族 | ✓✓✓ | 派发 + 注册;固定成员制 |
| | verlet_step, vv_advance, symple_step, symple_load | 闭族成员 | ✓✓✓ | 力只经 force_eval(头注在案) |
| | radau_step, radau_load, radau_boot, rad_clamp, rad_forget, beeman_forget | 闭族成员 | ✓✓✓ | RA15 自适应;boot 每轨迹重置/ B 系数携带(V 相 R 系判据) |
| | integ_proc_i(接口), stop_prop, stop_propi | 闭族 | ✓✓✓ | |

### 1.5 方法库——电子方法包(L3)

| 文件 | 例程 | 类 | (a)(b)(c) | 注 |
|---|---|---|---|---|
| src/methods/elec_adiabatic.f90 | ad_prop, ad_match | 电子方法包 | ✓✓✓ | 全 no-op 退化点(occ 恒 [1],零 rng) |
| src/methods/reg_adiabatic.f90 | reg_adiabatic_all | 注册材料 | ✓✓✓ | 一行注册 |

### 1.6 方法库——B′ 分布成员(L3;例程组: init/sample/公共面 + stop_* 具名停机族)

| 文件 | 例程 | 类 | (a)(b)(c) | 注 |
|---|---|---|---|---|
| src/methods/dist_boltzmann.f90 | boltz_init, boltz_sample + stop_*(4) | 方法库成员 | ✓✓✓ | 频率表为装配数据非输入键(P-1 家族:选中无表 fail-loud) |
| src/methods/dist_ebk.f90 | ebk_init, ebk_sample, ebk_morse, ebk_r_turn, quantize_step, turn_bounds, ebk_last + stop_*(4) | 方法库成员 | ✓✓✓ | Morse WKB 精确谱推导;零体系数据 |
| src/methods/dist_incident.f90 | inc_init, inc_sample, com_shift, add_mom + stop_*(3) | 方法库成员 | ✓✓✓ | 全族唯一编排位成员(专用选择通道 + 后置次序,samp_run 保证);体系级处方读 config 不入成员参 |
| src/methods/dist_rotation.f90 | rotation_init, rotation_sample + stop_*(4) | 方法库成员 | ✓✓✓ | |
| src/methods/dist_j.f90 | j_init, j_sample + stop_*(4) | 方法库成员 | ✓✓✓ | |
| src/methods/dist_surface_oscillator.f90 | surface_oscillator_init, surface_oscillator_sample + stop_*(4) | 方法库成员 | ✓✓✓ | 模型处方(SO 模型),质量/几何取 roster,零体系数据(S 相复议在案) |
| src/methods/dist_normalmode.f90 | normalmode_init, normalmode_sample + stop_*(4) | 方法库成员 | ✓✓✓ | 模式表为装配数据(P-1 家族) |
| src/methods/dist_barrier_excitation.f90 | barrier_excitation_init, barrier_excitation_sample + stop_*(4) | 方法库成员 | ✓✓✓ | |
| src/methods/dist_thermostat.f90 | thermostat_init, thermostat_sample + stop_*(3) | 方法库成员 | ✓✓✓ | |
| src/methods/dist_glo_target.f90 | glo_target_init, glo_target_sample + stop_*(3) | 方法库成员 | ✓✓✓ | FDT 代数装载;幽灵变量经 getter 交付(成员自有通道) |
| src/methods/dist_thermalize.f90 | thermalize_init, thermalize_sample + stop_*(4) | 方法库成员 | ✓✓✓ | 制备成员内部调演化服务(边界 4 合法:步长存取复原) |
| src/methods/seam_cast.f90 | seam_real, seam_int, seam_tag | 缝工具 | ✓✓✓ | 纯形式转换 |
| src/methods/reg_prop.f90 | reg_prop_all | 注册材料 | ✓✓✓ | 闭族积分器注册 |
| src/methods/reg_dist_*.f90 (×10) | *_all, *_init | 注册材料 | ✓✓✓ | 原子注册单元 + 装配缝(暂存键→成员 init;未选中守卫跳过) |

### 1.7 测量与后处理(L4)

| 文件 | 例程 | 类 | (a)(b)(c) | 注 |
|---|---|---|---|---|
| src/output/recorder.f90 | rec_reg_col, rec_sel_lvl, rec_sel_lst, rec_frame, rec_final, rec_fin_val, fw_reg, reg_row, check_ready, stop_rec | 测量 | ✓✓✓ | 列注册表 + 级/名单并集选列;骨架列无通用生产者(P-3 家族,不聚合不造值) |
| src/output/final_state.f90 | fin_eval, classify_i(接口), container_bind_classify, container_classify | 后处理 | ✓✓✓ | 广义 COM 相对量/散射角族;分类为容器槽(E19) |
| src/output/density.f90 | dens_states, eig_out, sum_count, swap_triangles | 后处理 | ✓✓✓ | Beyer-Swinehart;频谱组/本征对表无生产者(P-1 家族,驱动器尾注在案) |
| src/output/plot_hist.f90 | plot_hist, column_hist, bounce_hists, pos_map, site_hists, cell_fold | 后处理 | ✓✓✓ | 五直方图类;晶格参数化(roster 字段) |

### 1.8 仪器(L5)

| 文件 | 例程 | 类 | (a)(b)(c) | 注 |
|---|---|---|---|---|
| src/utils/rng.f90 | rng_init (+内部均匀/正态/gamma 变换) | 仪器 | ✓✓✓ | 全库唯一随机流(T5 零 RANDOM_NUMBER) |
| src/utils/linalg.f90 | eig_sym | 仪器 | ✓✓✓ | EIGN 逐字移植(C1 判据位级) |
| src/utils/geometry.f90 | cenmas, amom, angvel, rot_axis, rot_euler, rot_jkm, align_step, gl_nodes | 仪器 | ✓✓✓ | 无状态 (state,mass,list) 三元组例程 |
| src/utils/specio.f90 | read_poscar, read_xyz, elem_mass, split_tokens, fail | 仪器 | ✓✓✓ | 元素质量表为通用供给(非体系知识;头注在案) |
| src/utils/hessian.f90 | hessian_cd, force_proc(接口), require_finite | 仪器 | ✓✓✓ | 纯形式二阶差分;生产消费者零(P-1 家族四证据之一) |

### 1.9 体系材料(systems/;三问之 (c) 按定义豁免——审其反向:零程序知识)

| 文件 | 例程 | 类 | (c')零程序知识 | 注 |
|---|---|---|---|---|
| systems/reference/reference_pot.f90 | harm_vg, morse_vg, leps_vg, leps_diab, harm_set_constant, pair_dist, pair_unit | 体系材料 | ✓ | pure-math 头注在案;eV 交付契约;文献参考常数为解析形默认值(参数通道可覆盖) |
| systems/reference/reference_box.f90 | reference_load, reference_pes, reference_term, reference_classify, reference_term_on, reference_class_on, stop_box, stop_dup, stop_val | 体系材料 | ✓ | 参数文件通道 + 终止/分类判断(容器槽本职) |
| systems/reference/reg_reference.f90 | reg_reference_all, reg_reference_init | 体系材料 | ✓ | 拖入件;绑定经公开门面 |

## 2. 登记注记(全部有在案登记与去向,非红项)

| # | 注记 | 登记处/去向 |
|---|---|---|
| N1 | driver 的 dens_states/eig_out 不接线(无频谱组/本征对表生产者) | P-1 停放(装配期频率导出通道候选后续件) |
| N2 | driver 新轨迹复位只清 sta%t;单片段气相单员跑 q 跨轨迹携带(V-2 观察项) | case_b_layer §2 定量在案;语义变更走流程 |
| N3 | driver 统计相 n_valid 占位语义(无有效性判断,每轨迹计有效);lattice%q_ref = 0 占位 | C8 终审停放区(P-5);fin%ang(4..18) 消费面随 P-3 |
| N4 | recorder 骨架列(VI/OAMI/AMAI/AMBI/SDA/SDB 等)无通用生产者,不聚合不造值 | P-3 停放(采样器→统计相初值快照通道) |
| N5 | dist boltzmann/normalmode 的频率/模式表是装配数据,生产侧导出通道未接线(选中无表 fail-loud) | P-1 停放(与 N1 同族) |
| N6 | hessian_cd 生产消费者为零(判据程序与探针调用在案) | P-1 四证据之一 |
| N7 | 体系材料豁免句:systems/ 的本职即体系知识;其反向审计(零程序知识:不 use 门面以外的程序单元、不自建输出通道、报错走具名停机)全过 | 本报告 §1.9 |

## 3. 不变式三问深审抽查(S2.3;四处各一例程)

| 抽查位 | 例程 | 深审结论 |
|---|---|---|
| driver | venus(统计相聚合段) | 只经 fin_eval/rec_fin_val/plot_hist 消费;直读 roster 派生晶格字段属编排数据传递;无势能/常数/体系分支 ✓ |
| force_facade | force_eval | E27 后纯两步:容器调用 + 唯一换算 + T+V 折叠;零物理分支(switch 已退役);质量表只用于动能折叠 ✓ |
| recorder | col_tkin/rec_final 通道(t_kin 折叠 :407-415, occ/n_hop 列 :441-456) | 能量列经 e_conv 一次换算(I/O 侧,非力边界);occ/n_hop 经持有者 getter 只读;零体系知识 ✓ |
| sampler | samp_run | 派发纯数据驱动(scheme 码表);incident 后置与单片段守卫是编排契约非物理;成员绑定知识在成员 init ✓ |

## 4. 与总案 T 判据的对接

T7(本审计)= 全例程三问结论齐,零红项;T4/T5 机械 grep 在案(§0);
T6(driver 三相)由 §1.1 与深审抽查覆盖;T1/T2/T3 见 docs/plans/S.md 验收区。
