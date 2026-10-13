# 2026-10-13 裁定：DIST_SCHEME 逐片段括号参数（方案词参数的唯一合法通道）

## 背景与裁定

用户提出：希望在 `DIST_SCHEME` 的值里以 `方法名[参数A, 参数B]` 的形式就地给出
各片段的采样参数（例：两个 EBK 片段各处一个量子态）。可行性评估后用户裁定
三点口径：

1. **不写任何回退**——写了括号的片段自洽，缺的参数就是缺（走成员字段默认）；
2. **括号形式是唯一合法的参数通道**——21 个方案词成员键（T_VIB/N_VIB/N_ROT/
   J_ROT/T_ROT/N_OSC/T_OSC/N_LEVEL/K_OSC/E_VIB/E_STAB/N_E_BAR/E_BAR/T_BAR/
   T_GLO/GAMMA_GLO/W_GHOST/T_EQ/N_EQ/DT_EQ/E_SCALE）从分派表退役，顶层写它们
   即"未知键"中止；
3. 通道级成员键（入射束流 11 键 + 浴物种 T_BATH/NU_COLL）不随片段走，仍为
   顶层键（dispatch_tbl 34→13 行）。

由此实现了 2026-09-21 等权条款登记的"逐侧参数差异化"挂起项——按片段逐载体化。

## 语法（input.f90 parse_elist 重写）

- 值按**括号深度 0 的逗号**切分片段（括号内逗号是该载体参数的分隔符）；
  整值括号不配对先行指名中止；
- 片段项 = `词` 或 `词[名=值, 名=值]`：词查 dist 词表；括号前无词、空括号
  `[]`、非 `名=值` 项、空名/空值、名字符集外（限字母数字下划线）、超宽名、
  同括号重名、闭括号后残留文本、括号嵌套——均为指名 DIST_SCHEME + 行号的
  形式错误；大小写折叠（词与参数名折大写）；
- 解析只管形式；**参数名是否属于该词、值域如何，归成员 init**（分层不变）。

## 交付缝（与容器键通道对称）

- input 模块暂存逐片段参数对（dist_pars），公共入口 `par_int/par_real(word,
  ifrag, name, target [, found])` 供成员装配缝取参：缺席名保持成员默认
  （found=F——无回退的语义正是"默认仍在成员里"）；错误指名词+片段+参数名；
- `input_audit_dist_par`（driver 装配后调用，紧随 input_audit_declared）：
  暂存而无人消费的参数对指名中止——成员不认识的参数名不会被无声吞掉
  （判据 check_input (s5)）。

## 成员侧逐载体化

十个 dist 词成员的 params_t 标量全部改为逐载体数组（boltzmann t_vib(:)、
ebk car(:) 行（含逐载体谱常数拉取与 mu·r_eq² 转动惯量——不再"首载体一套
常数大家共用"）、rotation j_rot(:)（|L| 逐载体）、j t_rot(:)、
surface_oscillator 四参数、normalmode e_vib(:)、barrier_excitation 四参数、
glo_target 三参数（getter 增可选载体行）、thermalize t_eq 逐载体 +
N_EQ/DT_EQ 逐载体但**须一致**（一次共享平衡跑，不一致指名中止）、
wigner e_scale(:)）；BOLTZMANN 的温度行覆盖全部词载体（单原子词载体的括号
被消费与校验但不建模表）。RNG 消费序不变——参数不进随机流。

## 值中立性

- 8 案例（含 4 个 manual 变体）共 15 个运行重写为括号形式后 stdout 逐字节
  一致（logs/bracket_bytecompare.txt；唯一差异 CPU 计时行）；
- case_chain_frequencies/boltz 的顶层 `T_VIB=104` 与 normalmode 的
  `E_VIB=0.54205243` 是旧语法下**仅探针消费**的行（venus 运行从未读过——
  方案词是 ROTATION 占位），随本裁定改为探针自持常数（probe_chain_freq
  同步完成改基欠账清偿，w_per_rad/c_light 按 fs 折算，两腿 0 failures）。

## 顺带清偿的欠账（本批必要前置）

- **判据族单位改基漏账**（23349b7 记账只列了 check_driver/final_state/facade，
  漏了整个 check_samp_* 族与 probe_chain_freq）：十个 check_samp_* 的
  hbar_code/kb_code/e_conv → hbar_act/hbar_evfs/kb_act/kb_ev/wvn2ev·ev_dyn
  按成员单位逐处换算，kcal 时代字面量（42/8/5/12.5 kcal/mol）按 23.0605
  折 eV（与案例换算同链）；全部 0 failures；
- **判据配方闭包欠账**：八个 check_samp_* 头注释配方缺 sampler 的依赖闭包
  （bath/propagator/force/elec/spectrum/hessian/geom_init/linalg——自 bath
  并入 sampler 起就编不过），补齐；
- **check_driver**（23349b7 记账 deferred）：单位常数 5 处、t_life ×10 锚、
  我批的 plot_hist 面变动（lat_t 退役）与 electronic%n_surf 退役四处——
  全部清偿，23 腿 ALL PASS；
- 内嵌输入数组括号化 + 数组定长声明对齐。

## 判据与验证

- check_input：(s) 族 17 条新腿（s1/s2 暂存与取参、s3 十个退役键拦截、
  s4 六种畸形括号、s5 消费审计、s6 空值、s7 成员侧数值转型）——先红
  （logs/bracket_red_stepA.log 17 失败 + 步骤B编译红）后绿
  （logs/bracket_green_input.log 0 failures）；
- check_samp_* 十程序 0 failures（logs/bracket_green_samp_checks.log）；
- check_driver ALL PASS（logs/bracket_green_check_driver.log）；
- probe_chain_freq 两腿 0 failures；drill_auto_pes 全绿
  （logs/bracket_green_drill_auto_pes.log）；
- GUI：keyspec 镜像 13 行 + SCHEME_PARAMS 括号词表 + parse/format_dist_scheme；
  validate 逐片段规则；input_page 方案行内联参数字段；pytest 89 passed
  （余 smoke/pes_scan 两项系后台 drill 拖动 auto_sample 的环境性失败）。

## 文档

- docs/input_keys.md："采样方案自带参数"改为括号语法节（逐词参数表 + 缺省值
  + 双 EBK 片段示例）；
- AGENTS.md：等权条款 bullet 改写（逐载体实现 + 指向本记录）。
