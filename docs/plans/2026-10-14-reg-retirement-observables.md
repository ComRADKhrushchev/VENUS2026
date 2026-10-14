# 2026-10-14 reg_ 型容器支持退役 + 第 11 预留名 observables —— 裁决记录

## 裁决（用户，2026-10-14）

1. **第 11 预留名定名 `observables`**：自定义记录列槽（即 reg 的 columns 槽）接入
   Tier 1 自动接线，名称与输入键 OBSERVABLES 对齐（选列的键名即此）。
2. **reg_ 型支持以静默忽略方式解除**：prebuild 不再识别 reg_*.f90 的"reg 权威"
   地位——有无 reg 一律走预留名扫描；reg 文件作为普通源文件照编（模块无人调用，
   无害陪葬）；无弃用警告、无 loud-stop。
3. **五容器包装法转换的验收锚 = 判据全绿 + cau_tab 逐位****：五判据转换前后
   全绿且读数不变；cau_tab 用归档冒烟逐位比对。不为其余四容器另建冒烟案例。

背景方向（用户，2026-10-13）："今后应当解除对 reg_ 型的支持，毕竟已经没有
自定义的必要了"。本批即该方向的落地。

## 事实基础（2026-10-14 会话核实，均可复验）

- 五 cau reg 槽位清一色六槽：keys/force/method/term/init/columns；前五个均有
  预留名对应，唯 columns 缺——预留名词表扩展仅需一名。
- `export_spectrum` 槽全库零使用（systems/ 下无任何容器引用）——无需第 12 名；
  谱来源走采样成员（EBK/WIGNER），不经容器导出。
- Tier 2 的另两条在册存在理由均未使用：一文件夹多方法包（五家均为单包
  tdhf_hub）、复杂参数体（13 个容器键走 keys/load 预留名完全够）。
- 五判据（check_cau_tab/peri/emt2/emtnn/rst）编译配方只编 box/pot/tdhf/tab，
  不含 reg——reg 退役对判据编译零冲击；包装法只加薄包装，判据读数应不变。
- 容器现状：Tier 1 = reference（活性）+ fh2_leps/h2ag111PES/h3_bkmp2 +
  auto_sample 样例；Tier 2 = cau 五连（2026-10-13 07fb8ce 已完成常数同步）。

## 设计（2026-10-14 用户批准）

### ① 词表扩展——第 11 预留名 observables

- 签名同 keys()/load() 惯例：无参子程序；体内连调
  `rec_reg_col('列名', 函数)`（与 reg 的 register_columns 槽同形）。
- prebuild：扫描表（subroutine 名集）+ 重复/形态校验 + 中止文案 + use 行
  织入序（pes/keys/load/term/classify + prop/match/sync/init/method_name +
  observables 十一名）+ 十一名封闭词表声明一处扩展。
- auto_sample 样例增补一个 observables 观察列作织入正例；drill 织入腿
  随批 +1 断言（use auto_sample_met, only: ..., observables 归属模块视
  样例落点而定）。

### ② 五容器包装法转换（模板 = tier1-convert 四容器已验证同款）

- box 保留全部公有旧名（判据程序与 probe 钩子零改动），另暴露薄包装：
  - 纯势四名：pes/keys/load/term（term 包装直调体，体自带短路则同
    reference 先例）；
  - 方法包五名：prop/match/sync/method_name/init（cau 家族方法包在
    *_tdhf 模块——auto_sample 样例即此形态的织入正例；init 包装内
    ELEC_METHOD 门控照抄 reg 现体）；
  - 新 observables：包装体内转调 reg 现在的 register_columns 内容
    （rec_reg_col 逐行照抄）。
- keys/load 包装承载各家 13 个容器键的 input_declare_keys + pull
  （逐字照抄 reg 的 declare_keys/load_params 体，同 reference 转换先例）。
- 删五份 reg_*.f90。

### ③ reg 通道静默忽略

- prebuild 删除"文件夹含 reg_*.f90 则 reg 权威、不做织入"分支——扫描
  一律走预留名。
- reg_*.f90 落入 Makefile 源文件通配照编；其模块（reg_<sys>_slots）无人
  调用，为无害陪葬（转换后的容器已删 reg；仅用户手中未转换旧文件夹受此
  语义影响）。
- 旧 reg 容器拖入 systems/ 的自然失败模式 = 预留在场的既有 named abort
  （"wires nothing: expose exactly one 'pes'"）；该中止文案中"或写
  reg_*.f90 全容器"提法随批删除（reg 已非合法通道）。

### ④ 周边同步

- src/EXTENDING.md：删 reg 骨架节（§1b），Tier 2 描述改史层注记
  （"2026-10-14 退役"）。
- AGENTS.md：容器架条目并存句收窄（auto+reg 分支成死代码，删）；
  预留名词表十一名。
- docs/glossary.md：Tier 2 全容器条目加退役史层标记；预留名条目补
  observables。
- GUI（gui/core/pes_scan.py 等）：删 Tier2 提示分支与判词；扫描本就只认
  预留名，机制零改动。
- scripts/drill_auto_pes.sh：auto+reg 语义腿收尾；新增 observables 织入
  断言。

### ⑤ 验收（裁决锚）

1. 五判据转换前后全绿且读数不变（150/103/99/94/94 PASS 基线，
   logs/ 前后对照落盘）。
2. cau_tab 归档冒烟逐位（publish_min/cases/_smoke_cau_tab_ev 归档 +
   迁移输入，trajectories.txt 逐位一致；observables 列选入的健全性跑
   rc=0 columns 对应增长）。
3. 正位树 drill 30/0（含新 observables 织入腿）。
4. 双编译器构建 rc=0 零错误。
5. GUI 套件全绿（pes_scan 词表改动涉及）。
6. 五案例回归（reference 未动，预期全绿——常设门）。

### ⑥ 明确不动

- methods/ 家族 reg_*（reg_samp_*/reg_inc_*/reg_bath_*/reg_adiabatic 等
  方法注册机制，与容器 reg 是两回事）完全不动。
- publish_min/ 旧树与其内未跟踪副本不动（其使命随本批完成而终结的退役
  另议）。

## 开放边界（实施时注意）

- 静默忽略下 reg_*.f90 照编：若某旧 reg 引用了已被后续批删除的接口，
  编译会断——这是编译期自然暴露，不为旧 reg 保接口（用户裁决的代价，
  在案）。
- observables 织入的调用时点 = 输出装配组（与 reg 的 register_columns
  同位）；无 observables 的容器合法缺省（槽空跳过）。

## 执行注记（2026-10-14，实施计划 docs/plans/2026-10-14-reg-retirement-plan.md 五任务全数完成）

- **T1**（b0fba9e + 勘误 df2bd0b）：prebuild 四处名集+形态 case+auto_slots
  columns 绑定+auto_columns 发射；样例 io.f90 observables 正例。勘误起因：
  钻头退出陷阱第二次抹掉未提交样例编辑（同 2fb9be5 先例）——**教训固化：
  样例编辑必须先提交再跑钻头**。
- **T2**（84994b1）：cau_tab_wires.f90 十名包装（wires 文件形态，判据配方零
  触碰）；判据 150/0 与基线逐字一致；冒烟逐位复现归档；columns=11。
- **T3**（b714d02）：四兄弟 wires（pull 行自 reg 逐字抽取；peri 续行符一字
  之漏经编译级验证当场捕获）；四判据 103/99/94/94 与基线逐字一致；四家
  wires 语法编译全过（模块链 -J 临时 .mod + 主构建 .mod）。
- **T4**（becd12b）：prebuild 删 reg 让位分支；wires-nothing 文案去 reg 提法
  （prebuild+GUI 镜像两处）；EXTENDING/AGENTS/glossary/GUI 同步；静默忽略
  语义实证（旧 reg 容器拖入→新文案中止，logs/regret_silent_ignore.log）；
  GUI 套件 137/0。
- **T5**（终局门）：八案例回归 8×rc=0（logs/regret_cases.log，首跑五例 127
  为执行路径笔误——嵌套案例需 ../../../，非程序问题）；双编译器 gfortran/
  ifx 双 rc=0 零错误（logs/regret_make_gf/ifx.log）；钻头 30 实腿全绿
  （logs/regret_drill_final.log）——zero-touch 唯一红腿归因**外来并行编辑**
  （docs/custom_pes.md + docs/input_keys.md，非本批文件，归属并行文档会话，
  未代处置）；干净树 31/0 留待该会话落盘后随批复验。
