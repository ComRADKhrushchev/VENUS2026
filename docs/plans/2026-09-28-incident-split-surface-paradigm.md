# 2026-09-28 入射通道按范式拆分 + cau_peri 表面范式化 + dist_→samp_ 全家族更名

用户裁决（2026-09-28，会话批准的实施方案）：

## 动因

cau_peri（C@Au 周期 + TDHF 容器）的声子对比冒烟暴露两点：其一，该体系的物理形态是
表面散射（入射位点 = 胞内坐标），却因容器 system/ 用两个 .xyz 片段而落在分子范式，
横向几何只能借碰撞盘参数 b 表达；其二，dist_incident 一个成员内嵌双范式分支
（surf_case），单片段法令（"表面不是 list_atoms 片段"）使带动力学基片原子的表面
体系无法在表面范式下装配。

## 裁决四条

1. **入射通道按组装范式拆三成员**（前缀 `inc_`）：
   - `inc_surface`（词 `incident_surface`）：表面范式——恰一个 PROJECTILE 弹丸片段，
     其余片段为基片（**单片段法令放宽**：旧"nf==1"改为"nf≥1 + i_proj 点名"；
     case_surface_window 的单片段形态仍合法，cau_peri 的 C.xyz + Slab.poscar 双片段
     形态成为新在案形态）；弹丸 COM 速度全额指派（w1=1，表面=动量汇）。
   - `inc_pair`（词 `incident_pair`）：双反应物气体对——现行分子分支原样迁入
     （nf==2、双角色、mu_eff 分速 w1/w2、(t) 角色镜像）。
   - `inc_single`（词 `incident_single`）：平凡单反应物成员——draw/realize 皆 no-op
     （契约文档化）；单片段气相从"跳过守卫"改为"统一派发平凡成员"。
2. **dist_ 前缀全家族更名 samp_**：10 个逐片段分布成员 + reg + check 程序
   （注册词 'boltzmann' 等与 dist_tbl 码表不动）；入射家族用 inc_ 前缀。
3. **表面入射定点选项**：`N_AIM = 2` + `AIM_X`/`AIM_Y`（笛卡尔 Å，POSCAR 坐标系，
   与 b 偏移、R_SEP 后撒同坐标系可叠加）；键入 dispatch 表主 = incident_surface。
4. **cau_peri 容器表面范式化**：Slab.xyz → Slab.poscar（六角胞对齐引擎周期
   cons=2.95，surface_potentials.f90:47）；跑法 PROJECTILE=C（无 TARGET）、
   A_LAT=2.95、DIST_SCHEME=BOLTZMANN,SURFACE_OSCILLATOR。

## 随附裁决三条

- **按范式守卫注册**（保 `members=11` 计数与五 case 档案稳定）：reg_inc_surface_all
  仅 surface_model/=0 注册；reg_inc_pair_all 仅 model=0 且 nf>1；reg_inc_single_all
  仅 model=0 且 nf<=1——每次进程恰武装一个入射成员。
- **RANDOM_ORIENT 豁免基片片段**：取向资格掩码在表面范式下限 i_proj（弹丸）——
  表面晶格取向不可转。
- **表面在位束 B_MAX=0 合法化**（限 N_B=0 固定瞄准）：定点/在位研究需要 b≡0；
  气体对成员保持 b_max>0 旧法令。

## 共享组件

`methods/beam_laws.f90`（level-2 组件，无 reg）：束流法则四层抽取（取向/方位/瞄准/
能量）+ 参数验证 + com_shift/set_com_vel 放置助手，`beam_t` 类型化实例——
成员间不互相依赖的立法以组件形态满足（EXTENDING.md level-2 章）。

## 回归口径

- gfortran（档案编译器）逐字节：b_layer/chain_frequencies 与在档 results.txt 一致；
  harmonic_scatter/radau/surface_window 与 2026-09-27 set_com_vel 修复后基线一致
  （三 case 档案随该批重入档）。
- members=11 计数不变（守卫注册）；抽样序（逐片段成员→取向矩阵→入射）不变；
  inc_single 不耗 RNG 不写状态——单片段气相流位级不变。

## 记录

- 实施批：本文件同日（会话执行）；判据拆分 check_inc_surface / check_inc_pair /
  check_inc_single（旧 check_dist_incident 拆解）+ check_samp_* 更名族。

---

## 落地与验证（2026-09-28 执行批）

- 落地：beam_laws 组件 + inc_surface/inc_pair/inc_single 三成员 + 守卫式 reg 三件
  （reg_inc_single 无 _init——平凡成员无参数）；sampler 统一派发 + samp_incident_word()
  + 取向基片豁免；driver 接线；input dispatch 表（35 行，AIM_X/AIM_Y 入表）；
  dist_*→samp_* 全家族更名（10 成员 + 10 reg + 10 check + probe/facade 配方）；
  cau_peri system/Slab.xyz→Slab.poscar（六角胞 a=2.95）；Makefile SRC_TREE 过滤
  build/ 整目录（并发会话探针源混入 wildcard 的修复，K3 法令推论）。
- 判据：check_dist_incident 拆三（check_inc_surface/pair/single，新腿 (v) 定点/
  (w) 基片放宽/在位束/A_LAT 整数比/取向豁免）；check_sampler 增 (p) 词解析 (q)
  基片豁免、改 (b)(o) 统一派发语义；顺带修 8 个在案陈红腿（两相化后 sample→draw
  守卫名过时：j/rotation/boltz/glo/normalmode/surface_osc/wigner/barrier 的
  sample-before-init 腿）。check_facade 引用不存在的 input_parser 模块——预存断裂
  （并发会话在途），本批跳过登记。
- 回归：gfortran+ifx 双净重建绿；七 case 腿六腿逐字节一致（surface_window/radau/
  b_layer×3/chain_freq×2），harmonic_scatter ULP 级差异（868 行观测最大相对
  1.98e-15 on v_pot；stdout 五证据行逐字节一致）——判读为采样链重组的编译器收缩
  噪声，随三 case 重入档批（09-27 set_com_vel 批已在流程）收口；members=11 计数
  全程稳定（守卫注册）。
- 冒烟（cau_peri 表面范式，本地激活 dance）：定点腿 AIM=(0,0)+b=0.15 在位顶位
  ——C 落 (0.15,0,6.0) 全额束流（表面动量汇，w=1），slab 保持 POSCAR 框架不动，
  RC=0 漂移 6.4e-4；胞内均匀腿落六角胞随机位；TDHF 武装腿 2/2 有效、占据列活性
  （GS 种子精确）、漂移 57 kcal（瞬态在族）。容器已复位 reference 并复验
  surface_window 位级一致。产物 smoke_cau_peri/surface_smoke/。


---

## 追加裁决：表面范式输入禁止书写 B_MAX（问题提出 2026-09-28，立法 2026-10-01）

**问题**（用户，2026-09-28）：周期表面的入射位点本应胞内随机，为什么结果里会有
b 依赖？冲击参数盘是气相碰撞概念，表面横向落点应由落点模式独占表述。

**落实**（2026-10-01）：`inc_surface_init` 增加输入缓冲扫描——表面范式下任何
**真实书写**的 B_MAX 行（`buffer_line > 0`，排除框架的缺省注入行）触发命名中止：
"the lateral placement belongs to the N_AIM drop-point modes (origin / unit cell /
fixed site); the collision-disk radius is a gas-phase pair concept"。`input.f90`
的 B_MAX 键注释同步（PAIR 通道专用）。N_AIM 的三种落点模式（原点 / 胞内均匀 /
固定点）不变。

**端到端验证**：写了 B_MAX=2.0 的 case_surface_window 输入 → rc=1 命名中止 ✓；
移除后重跑 8/8 有效（守恒 6.7e-4 kcal/mol）✓；check_inc_surface 全绿（其腿直
操参数记录、不经输入缓冲，语义不受影响）✓。

**连带修正**：case_surface_window 原输入为"N_AIM=1 胞内均匀 + B_MAX=2.0/N_B=1
碰撞盘"的双重展宽（首帧 y 曾达 6.29 出胞）——正是本裁决针对的形态；输入修正
（去 B_MAX、N_B=0）并全量重出 results.txt（纯胞内投放：y ∈ [0.30, 5.65] ⊂ 胞，
双跑逐字节一致）。云端 cau_tab 系综批次输入同步清理。
