# 显式弹丸/靶声明（PROJECTILE/TARGET）—— 三项决定与两阶段落地（已完成）

> 状态：**两阶段均已落地并过完整验证四腿（2026-09-21）**。阶段一＝角色键；阶段二＝
> 同方案词多片段放开（派发去重 + 十成员多载体化）。用户决定三条见下；平行批次
> （频谱链 + 多容器，后者收拢至 systems/inactive/）已合入后实施阶段二。

## 一、用户裁决（2026-09-21）

1. **键形态：双键显式**——`PROJECTILE = <片段名>` 与 `TARGET = <片段名>` 各自点名，
   交叉校验；片段名 = 体系定义文件词干（sysdef 暂存名），大小写不敏感匹配。
2. **缺省行为：多片段体系必填，无位置回退**——遗漏/拼错即命名终止（破坏性变更：
   全部多片段/表面输入必须补键）。
3. **范围：一并解除同方案词双片段终止**——`DIST_SCHEME = THERMALIZE,THERMALIZE`
   合法，一份成员参数平等服务其全部载体片段（"不区分具体分子"）。

必填矩阵（范式条件，仿 A_LAT/BOXPAD 先例，裁决于 roster_load 装配期）：

| 范式 | PROJECTILE | TARGET |
|---|---|---|
| 气相双片段 | 必填 | 必填（与 PROJECTILE 互异） |
| 表面范式 | 必填（命名唯一片段） | 出现即终止（表面不是 roster 片段） |
| 单片段气相 | 出现即终止 | 出现即终止 |

背景动因：原硬规则"roster 片段 1 = 弹丸、片段 2 = 靶"（文件名字典序承载身份，
重命名文件会静默改变谁是弹丸）——本件将角色承载从位置改为显式声明。

## 二、阶段一交付物（已入库工作区）

- `src/input.f90`：key_tab 25→27 行（PROJECTILE/TARGET，kind_char，slotw=64，
  cond=.true.——不注缺省、不做无条件必填、无槽位写入）；map_slots 表面窗口键
  分支扩含两新键。
- `src/globals/config_roster.f90`：`frag_t` 新增 `name`（len=64，建表时从
  staging_frag_name 填充）；`roster_t` 新增 `i_proj/i_targ`（0=未指派）；
  roster_load 第 7 步 **role_reconcile**（必填矩阵 + 名字解析 + 互异校验 +
  自证式报错——未命中时列出全部可用片段名与行号）；roster_clear/reset_derived
  同步清零。
- `src/methods/dist_incident.f90`：绑定块改读 roster%i_proj/i_targ（表面/气相
  两分支的 at_proj/at_targ/m1/m2/mu_eff/w1_frac/w2_frac）；行域守卫；报错文案
  去位置硬编码。
- `src/output/final_state.f90`：表面 μ = frag(i_proj)%mass；ang(19)/ang(20) 用
  vcms(i_proj,:)；表面角色未解析守卫；A/B 通道标签维持位置性（明记非本件范围）。
- 案例补键：case_harmonic_scatter 与 case_radau_conservation 加
  `PROJECTILE = frag_a` / `TARGET = frag_b`；case_surface_window 加
  `PROJECTILE = projectile`；b_layer 三子跑与 chain_frequencies 两子跑为单片段
  气相，按裁定不补。

判据程序：check_input 新增 (q1)–(q8) 八腿（正链解析/缺键/未知名/同名/单片段
误用/表面误用 TARGET/表面缺 PROJECTILE，w-pr- 子分派 + window_pair 材料）；
check_roster 双链断言名字与角色行 + (g) 清零检查；check_dist_incident 三 setup
填角色行 + 新增 (t) 角色镜像腿（固定模下互换角色：放置互换、相对运动不变）；
check_final_state crit_a 填角色 + wipe_roster 清角色 + 新增表面角色镜像腿
（i_proj=2 → 表面量按 m_b/m_a 缩放，b_fin 不变）；check_driver 四个两片段材料
（la/lec/lr/lt）补键 + 头注编译配方补 spectrum/hessian 模块。

## 三、阶段一验证证据（2026-09-21，当时树状态：reference + 频谱链在途文件）

1. gfortran 全库构建链接绿；
2. check_input / check_roster / check_dist_incident / check_final_state 四程序
   重编重跑全绿（0 failures，含全部新腿）；
3. 五案例回归：七个 venus.e 跑（含补键三案例）trajectories.txt 与入库锚
   **逐字节一致**；chain_frequencies 两探针腿 0 failures；
4. T4/T5 对本批零暴露（无新换算常数、RANDOM_NUMBER 零命中）。

未跑腿（被并发阻断，阶段二验证轮统一补）：ifx 双编译器腿；check_driver
重编重跑（材料已就绪，编译被 registry_gen 中新出现的 reg_cau_emt2 打断）。

## 四、阶段二交付物（已落地）

1. **sampler 派发去重**（src/facade/sampler.f90）：samp_run 第 4 步改为"每个不同
   成员行恰调用一次"（roster 首现位置）；同词双选 = 一次调用服务两个载体。
2. **十分布成员多载体化**：载体列表绑定（零载体终止保留、>1 终止全部删除）、
   派生态按载体行、sample 循环载体、随机流按载体序；**一份参数平等服务全部载体**。
   - 装配数据类字段逐载体一表：boltzmann `w_mode_c(:)`（reg 缝逐载体
     `spectrum_frag_wvn` 拉取——与频谱链批次的载体模型融合）、normalmode
     `mode_tbl_t tbl(:)`（reg 缝逐载体 `spectrum_frag_w`）、barrier_excitation
     `saddle_tbl_t tbl(:)`（仍无生产侧拉取链，选即在未分配检查终止——原姿态）。
   - getter 统一**可选载体参（缺省=载体 1）**：thermalize/j/ebk/glo_target/
     thermostat/surface_oscillator 的既有单载体调用面零改动可编译（probe_glo
     未动）；boltz_last_n(k, m) 双参（probe_chain_freq 已适配）。
   - thermalize 语义：逐载体 Maxwell 播种+漂移扣除（载体序）后**一次**共享热化跑。
3. **判据**：check_sampler 新增 (l) 去重腿（同词双选=单调用，stand_both 多载体
   形态）；十个成员判据各加"双载体平等处理"正腿（thermalize (k)/rotation (l)/
   j (j)/ebk (l)/thermostat (g)/glo_target (k)/surface_oscillator (m)/
   boltzmann (l)/normalmode (p)/barrier (t)）；六条旧"双选终止"负腿就地转为
   合法新语义正腿（头注注明 retired with multi-carrier support）；
   check_dist_boltzmann 的 (g)/(j) 改为内表未分配/表数-载体数失配；
   probe_chain_freq 参数构造与 getter 适配。

## 五、验证四腿（2026-09-21 终扫，全绿）

1. 双编译器全量重建：`make` 与 `make FC=ifx` 均 0 错误。
2. 18 个判据程序按各自头注配方重编重跑全 0 失败（含 check_driver ALL PASS、
   check_final_state ALL PASS；check_recorder/check_propagator 无关抽查零附带损伤）。
3. 五案例回归：六个 venus.e 跑（补键三案例 + 单片段三案例）trajectories.txt 与
   入库锚**逐字节一致**；chain_frequencies 两探针腿 0 失败（参数形状适配后）。
4. T4/T5：23.0605 仅 force_facade、生产树 RANDOM_NUMBER 零命中（本批零新换算
   点、零新随机流）。

单载体路径位级保持由第 3 腿逐字节锚直接证明（去重与多载体化不改变任何现有
案例的随机流消耗次序）。

## 六、并发注记（历史，留档）

- 阶段一执行期间另一会话并行落盘频谱链（utils/spectrum.f90 +
  facade/spectrum_stage.f90 + reg_dist_boltzmann/normalmode 改造）与多容器
  （fh2_leps/h2ag111PES/h3_bkmp2/cau_emt2，后收拢至 systems/inactive/——
  prebuild 的 systems/*/*.f90 单层通配天然排除）。两批文件集零重叠；阶段二
  在其合入后实施，boltzmann/normalmode 的多载体改造以其终态为准融合。
- 阶段一验证轮中 ifx 腿与 check_driver 重跑曾因此推迟，并入第五节终扫补齐。
