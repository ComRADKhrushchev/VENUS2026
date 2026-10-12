# 2026-10-12 纯势容器（Tier 1 auto-PES）——零 reg 自动接线设计 + 裁决记录

## 0. 动因与裁决链

用户对现行自定义 PES 接入的三点不满（2026-10-12）：

1. reg 注册方式对新用户过于繁琐（pot + box + reg 三件套，真实最小例 F+H2 LEPS 约 460 行，
   其中 reg 88 行 90% 与 reference 逐字相同）；
2. 本质上仍要写一个格式特定的势能面子例程，谈不上"直接接入"；
3. 可注册项八九个（container_sched 八槽 + 各注册表）太多。

本轮用户裁决（逐条问答，2026-10-12）：

- **D1 目标形态**：遗留 Fortran 源码 + Python/脚本势两类；**Python 桥接后置**，
  本批只做纯 Fortran 接口的精简（外部势通道另案）。
- **D2 激活机制不动**：维持 2026-09-21 裁决（systems/ 顶层一次一容器，搬文件夹即激活）。
- **D3 入口契约仅全量**：Tier 1 只接受 `(q, mass, g, v)` 全赋值入口；
  "能量 only + 框架内算差分"的提案被否。
- **D4 方案一 + 魔数名轻量槽**：零 reg 自动接线当选；且 auto 检测除必备 `pes` 外
  再扫可选魔数名（`keys`/`load`/`term`/`classify`），接受检测矩阵复杂化换取
  Tier 1 能力更全（用户在"三层处理 vs 魔数名轻量槽"中选后者）。
- **D5 非绝热方法包入词表**（同日追加问询）：`prop`/`match`/`sync`/`arm` 四名入表
  （配对约束见 §1）；注册名**不由模块名/文件名派生**（两条派生规则均被否），由
  第十个魔数名——显式定名函数 `method_name()`（作者数据，返回
  `character(len=32)`）提供；生成代码在注册行与 arm 门两处调同一函数，一致性由
  构造保证。自定义记录列（`rec_reg_col` 多参注册词汇）仍留 Tier 2。

设计分节（契约 / 机械 / 缺省与升级 / 文档判据）均已逐节获用户同意，D5 增量
（四名 + 定名函数）另经专项问询获批（2026-10-12）。

## 1. Tier 1 纯势容器——十名词表（封闭集合）

一个 systems/ 顶层文件夹，**无 `reg_*.f90`**，即纯势容器候选。作者面 = 文件夹内
一个（或多个）.f90，暴露下列**魔数名**中的若干个（名字即约定，签名即 E14/E19/E27
既有接口，编译期核对）。词表分两组：

**势能面组**

| 魔数名 | 形 | 必备性 | 签名要点 | 织入槽位 |
|---|---|---|---|---|
| `pes` | subroutine | **恰一个** | `(q, mass, g, v)`：q[A]/mass[amu] 进（q/g 一维展平，size=3×nat——`container_pes_i` 实际形态，force_interface.f90:49），g[eV/A]/v[eV] 全赋值出（E27 力交付契约原样） | force（registration） |
| `keys` | function | 可选 0/1 | `character(len=16)` 数组（`reference_keys` 同型），容器键词表 | keys（parse 前） |
| `load` | subroutine | 可选 0/1 | 无参；模块**自持 `pull`**（`use input, only: pull`），缺键保缺省；**亦即势能面侧的 seam 期初始化钩子**——pull 之后（或无键）就地读参数文件/分配（NN 权重模式：`keys` 声明文件名键 + `load` 读入模块态，cau_emtnn/cau_tab 同款压扁版；无键 `load` = 约定路径的纯初始化，亦合法）。时序保证：load 先于 arm/一切力消费者（container_sched_init 固定序 + 装配立法注记） | load（seam） |
| `term` | function | 可选 0/1 | `term_i` 同型：`(sta) → logical` | term（seam） |
| `classify` | function | 可选 0/1 | `classify_i` 同型：`(fin) → name(len=32)` | classify（seam） |

**方法包组**（自定义非绝热方法，E27 容器方法包的压扁版）

| 魔数名 | 形 | 必备性 | 签名要点 | 织入 |
|---|---|---|---|---|
| `prop` | subroutine | 与 `match` **成对**可选 | `elec_prop_i` 同型：`(sta, dt)` 电子传播成员 | `method_reg` 第二参（registration） |
| `match` | subroutine | 与 `prop` 成对可选 | `mqc_match_i` 同型：`(sta)` 统计匹配成员（跳变、`occ_set`/`n_hop_set` 写回） | `method_reg` 第三参 |
| `sync` | subroutine | 可选，须随对 | `sync_proc_i` 同型：`(s)`——a↔rho 表象换算尾巴 | `method_reg` 第四参 |
| `arm` | subroutine | 可选，须随对 | 无参；运行期初始化（分配/参数落位，参数自持走 `keys`/`load` 通道） | arm 槽（seam），**门控归生成代码** |
| `method_name` | function | **随对必备** | 返回 `character(len=32)` **小写**名——注册名（作者数据；input_qct.txt 的 `ELEC_METHOD` 值。输入层将键值小写化后查表，input.f90——返回非小写名即查表落空，故定为契约） | 注册行与 arm 门**两处同源调用** |

规则：

- **词表封闭**：不再增名。边界新表述：十个名字全部是**零参槽位体或固定签名
  成员**；唯一的注册词汇是方法名，且由魔数函数显式提供。带多参注册词汇的扩展
  （`rec_reg_col(名, 提取器)`、一个文件夹注册多个方法包）不入 Tier 1——那是
  Tier 2 全容器的事。
- **势能面组在场即绑（无条件）**；键缺席门控（如 R_TERM 缺席 → 跑满 MAX_STEPS）
  由模块内部 found 旗标实现（范本示例见 §3）。
- **方法包组注册无条件、武装有门**：`method_reg` 行总是注册（名字查表语义，与
  methods/ 全员注册同理）；生成代码的 arm 门
  `if (trim(electronic%method) == method_name()) call arm()` 与注册行同源——
  两处调同一 `method_name()`，不可能不一致。方法组在场而跑绝热参考
  （ELEC_METHOD 缺省/adiabatic）→ 行惰性、arm 不调（绝热参考面白送）。
  occ/n_hop 统计由 elec_interface 既有持有者承接（occ_set/n_hop_set 通道）。
- **名字可居文件夹内任意模块**（检测按名找归属模块；通常全在一个模块；方法组
  各名亦可跨模块，唯注册名取 `method_name()` 返回值与归属无关）。
- **公共性**：入口须 public——织入体 `use <mod>, only: <名>` 编译期核对，
  private 即编译错（文档明示预期错误形态）。
- **reg 权威**：文件夹有 `reg_*.f90` → 十名全部忽略，一切照旧（自动接线静默退出）。

### 有名中止（prebuild 期，全部指名）

1. 无 reg 且 `pes` 计数 0 —— "systems/<X> 已安装（文件夹在场）但未接线：暴露一个
   `pes` 入口，或写 `reg_<X>.f90` 走全容器"（比今天"编译过、首次力调用才死"更早更响）；
2. 任一名计数 ≥2 —— 指名两处定义；
3. 名命中但形不符（如 `subroutine keys`、`subroutine method_name`）—— 指名
   （把安静失败变有声）；
4. 自动文件夹 ≥2 —— 一进程一容器（2026-09-21 裁决）的构建期前置；
5. 自动容器与 reg 容器（cau 族）并存激活 → 既有排程第二填 STOP 1 原样生效
   （测试自动容器时暂存该 reg 容器；两自动文件夹归第 4 条构建期中止——D6 裁决后
   reference 自身即自动容器）；
6. `prop` XOR `match` —— 写半个传播器不算包，指名；
7. 对在场而无 `method_name`，或 `method_name` 在场而无对 —— 指名；
8. `sync` / `arm` 在场而无 prop+match 对 —— 指名。

### 零新增监管（既有兜底照管一致性问题）

- `load` 无 `keys`：合法——无键可拉时作**无输入依赖的纯 seam 期初始化**（读约定
  路径文件/分配/常量落位）；封闭语法本就不放行未声明键，pull 落空保缺省；
- `keys` 无 `load`：既有 `input_audit_declared` 兜底中止（buffered 而无人消费）；
- 签名不符：编译期接口核对即错；
- 撤掉 `pes`：prebuild 中止 1；
- `method_name` 返回值与既有方法名撞名（如 'adiabatic' 或 methods/ 包名）：
  `method_reg` 既有 duplicate 名中止（运行期兜底，elec_interface.f90:170）；
- 名字字面超 32 字符：编译期错（result 声明 `len=32`，字面超长即报）。

## 2. 构建机械

### prebuild.sh（zmods 收集后新增 auto 扫描与织入，约 +50 行 bash）

- 对 systems/ 顶层每个含 .f90 且无 `reg_*.f90` 的文件夹：按名扫
  `^[[:space:]]*(subroutine|function)[[:space:]]+<名>[[:space:]]*(`（大小写不敏感），
  计数与形核对如 §1；归属模块 = 名定义行之前最近的 module 语句（awk 扫描，
  复用既有 mod_owner 语义）。
- 命中恰一文件夹时，织入**本就生成的 registry_gen.f90**（无新源文件、无新构建产物）：

```fortran
module registry_gen
   use container_sched, only: container_sched_keys, container_sched_all, &
                              container_sched_init
   use force_interface, only: container_bind_pes, container_bind_term   ! 按在场织入
   use final_state,     only: container_bind_classify                   ! classify 在场时
   use input,           only: input_declare_keys                        ! keys 在场时
   use elec_interface,  only: method_reg                                ! 方法组在场时
   use config,          only: electronic                                ! arm 在场时
   use mysys_pot,  only: pes, keys, load, term, classify                ! 按检测到的名列表；
                                                                         ! 名跨多个模块时逐模块各一行 use
   use mysys_tdhf, only: prop, match, sync, arm, method_name            ! 方法组同规
   ...
contains
   subroutine auto_slots()                       ! 织入的自动手（只填在场槽）
      call container_sched_bind(force    = auto_force, &
                                [method  = auto_method,] &              ! 方法组在场时
                                [arm     = auto_arm,] &
                                [keys    = auto_keys,] &                 ! 按在场
                                [load    = auto_load,] &
                                [term    = auto_term,] &
                                [classify = auto_classify])
   end subroutine auto_slots
   subroutine auto_force();  call container_bind_pes(pes);          end subroutine
   subroutine auto_keys();   call input_declare_keys(keys());       end subroutine
   subroutine auto_load();   call load()                            end subroutine   ! 模块自持
   subroutine auto_term();   call container_bind_term(term)         end subroutine
   subroutine auto_classify(); call container_bind_classify(classify) end subroutine
   subroutine auto_method()                                  ! 注册行（无条件）
      call method_reg(method_name(), prop, match [, sync])
   end subroutine auto_method
   subroutine auto_arm()                                      ! 武装（门控在此，同源派生）
      if (trim(electronic%method) == method_name()) call arm()
   end subroutine auto_arm

   subroutine registry_gen_keys()
      call auto_slots()                           ! 一行，置于 zmods 调用之前
      ...（既有内容不动）
```

- 生命周期零新增：keys 前于 parse、force 于 registration、seam 归既有守卫——
  排程（2026-09-22 裁决：时间归排程、内容归槽体）原封不动，auto_slots 只是
  "一只只填在场槽的手"。

### Makefile（两行 + 头注改写）

`registry_gen.o` 的静态前提从 `$(OBJ_REG)` 扩为 `$(OBJ_REG) $(OBJ_SYS)`：
闭掉"生成文件 use 非 reg 模块"这条新边在 depend.mk 单轮迟滞窗口里的首编竞态
（拖入即 make，**首编即绿**——实现时实测核验并在判据中固化；重建代价为一个
微小文件的偶发重编）。头注同步改写（"reg-object prerequisites" →
"systems/methods object prerequisites (auto-wired use edges included)"）。

depend.mk 侧零改动：registry_gen 的 use 边（用户模块）走既有 mod_owner →
emit_deps 机制自动覆盖。

## 3. 缺省语义与升级路径

### Tier 1 缺省（全部是既有"空槽合法跳过"形态，零新物理）

| 绑定项 | 缺省 | 依据 |
|---|---|---|
| 终止 | 跑满 `MAX_STEPS`（程序键） | `container_term` 未绑定恒 false（force_interface.f90:130-131） |
| 分类 | 无分类运行 | final_state 既有无绑定形态 |
| 输入键 | 无容器键，语法封闭如常；PES 参数=源内常数或经 `keys`/`load` 声明拉取 | E28 通道的压扁版 |
| 电子方法 | 方法组缺席 → 绝热参考（ELEC_METHOD 缺省）；在场未选中 → 注册行惰性、arm 不调（既有查表语义） | method/arm 槽 |
| 谱 | SPECTRUM_SOURCE=COMPUTE 照常（力槽绑定即可算） | 2026-10-05 谱立法 |
| 几何 | SYSTEM_DIR 输入键指案例侧 system/，或容器自带 system/ 子文件夹（两条既有通道） | 2026-09-22 构型裁决 |

### 作者面范本（一个文件即全部；完整可编译版随实现入库为判据材料）

```fortran
! systems/mysys/pot.f90 —— Tier 1 纯势容器
module mysys_pot
   use input, only: pull                  ! load 需要；纯势（无键）可省
   implicit none
   private
   public :: pes, keys, load, term        ! 暴露什么，接什么
   real(8) :: k_harm = 5.0d0, r_term = 12.0d0   ! 源内常数即缺省
   logical  :: r_term_found = .false.
contains
   subroutine pes(q, mass, g, v)                          ! 必备（q/g 一维展平，size=3*nat）
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      ! ...物理本体...
   end subroutine pes
   function keys() result(kk)                             ! 可选
      character(len=16) :: kk(2)
      kk = [character(len=16) :: 'K_HARM', 'R_TERM']
   end function keys
   subroutine load()                                      ! 可选（模块自持 pull）
      call pull('K_HARM', k_harm)
      call pull('R_TERM', r_term, r_term_found)
   end subroutine load
   logical function term(sta)                             ! 可选（键缺席→false）
      use state, only: state_t
      type(state_t), intent(in) :: sta
      term = .false.
      if (r_term_found) term = .false.                    ! 实际判据由作者给
   end function term
end module
```

非绝热方法包（方法组；同文件夹另一模块，需要时才写）：

```fortran
! systems/mysys/tdhf.f90 —— Tier 1 方法包（方法组五名）
module mysys_tdhf
   use state, only: state_t, elec_t
   implicit none
   private
   public :: prop, match, sync, arm, method_name
contains
   function method_name() result(s)      ! 随对必备：注册名（ELEC_METHOD 的值）
      character(len=32) :: s
      s = 'mysys_tdhf'                   ! 任意名字，作者定；两处调用同源
   end function method_name
   subroutine prop(sta, dt)              ! 成对必备：电子传播（elec_prop_i 同型）
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
      ! ...振幅/密度矩阵演化本体...
   end subroutine prop
   subroutine match(sta)                 ! 成对必备：统计匹配（mqc_match_i 同型；
      type(state_t), intent(inout) :: sta  !  跳变与 occ_set/n_hop_set 写回在此）
   end subroutine match
   subroutine sync(s)                    ! 可选：a↔rho 表象换算尾巴（sync_proc_i 同型）
      type(elec_t), intent(inout) :: s
   end subroutine sync
   subroutine arm()                      ! 可选：选中后的运行期初始化（门控归生成代码）
      ! ...分配/参数落位（参数经 keys/load 通道已入模块态）...
   end subroutine arm
end module
```

### 升级 = 写 reg，pot 十名全部直接复用

哪天需要自定义记录列 / 一个文件夹多个方法包 / 更复杂的参数体了：同文件夹加
`reg_<sys>.f90`（+ 可选 box）——prebuild 见 reg 即静默退出自动接线、reg 权威
接管编排。**pot 一字不动**：reg 区体直接 `use mysys_pot, only: pes, keys, load, term`
与 `use mysys_tdhf, only: prop, match, sync, arm, method_name` 引用同一批函数
（`container_bind_pes(pes)`、`method_reg(method_name(), prop, match, sync)`……）。
所谓 box 在 Tier 1 里被压扁成了 pot 模块自己的 `keys`/`load` 两函数。

## 4. 非目标与登记后续件

- **Python/脚本势桥接**（D1 后置）：子进程管道协议 / i-PI socket 方案本轮只留缝
  （十名词表封闭、外部势另立容器形态另案裁决）。
- **能量 only 入口 + 框架差分**：D3 否决，不做。
- **通用距离终止键**（TERM_R 类，程序级、范式感知——气对按碎片间距/表面按弹射体
  高度）：真正的新裁决（触碰 E14/E19 容器属位），登记为后续件候选，不混本批。
- **激活机制改革**（SYSTEM=xxx 输入键选择、多容器共存常驻编译）：D2 明确不动。

## 5. 判据与验收

1. **判据程序** `scripts/check_auto_pes.f90` + 演练脚本（判据族同款）：
   - 正面：范本容器拖入 systems/ → **首次 make 即绿**（首编即绿核验点）→
     `build/<FC>/registry_gen.f90` 织入内容断言（use 行、auto_slots、调用行）→
     MAX_STEPS 安全网跑通一条；`keys`/`load`/`term`（键在场与缺席两分支）/
     `classify` 路径各覆盖一条；方法组两条——`ELEC_METHOD = <method_name()>` 值
     选中跑通一条（arm 被调的见证），ELEC_METHOD 缺省（adiabatic）下方法行惰性、
     arm 未调的见证各一；
   - 负面：pes 0 / pes ≥2 / 任一名 ≥2 / 双自动文件夹 / 形不符（subroutine keys），
     以及配对四条——`prop` 无 `match`、对无 `method_name`、`method_name` 无对、
     `sync`/`arm` 无对——与 `method_name` 返回 'adiabatic' 撞名（method_reg 既有
     duplicate 中止），中止文案逐条断言；
   - T2 式零触碰演练：拖入 → 判据绿 → 拖出重建 → git tracked 零改动。
2. **回归**：五案例 results.txt 逐字节一致；双编译器全量重建绿；T4/T5 grep 复验
   （23.0605 仅 force_interface；RANDOM_NUMBER 零命中）。
3. **文档落点**：EXTENDING.md §1 重写为 Tier 1 领篇（十名词表 + 两组范本 + 升级
   注记），现行 8 槽三件套降为 Tier 2 全容器；README.md / src/CLAUDE.md / AGENTS.md
   注册段同步增补；本文件即裁决记录。

## 6. 状态

- 2026-10-12：设计定稿，用户逐节批准（D1-D4 裁决链 + 四节设计）；同日 D5 增量
  批准（方法包四名 + `method_name()` 定名函数，词表五名 → 十名），文档随批修订。
  待实现。
- 2026-10-12：实现落地（实施计划 docs/plans/2026-10-12-auto-pes-implementation.md）——
  prebuild 织入 + Makefile 静态前提 + auto_sample fixture + check/drill 判据 + 文档四件。
  收尾回归：双编译器干净重建绿、check_auto_pes ALL PASS、drill 干净树 20/20 全绿；
  六案例零扰动——surface_window 与 spectrum_source 四腿对在档基线逐字节一致（后者按
  容器架规程活置 h2ag111PES/h3_bkmp2 跑，trajectories 逐位同），b_layer/
  chain_frequencies/harmonic/radau 对 results.txt 的分歧均为先在陈档（members 12→11
  迁移与 2026-10-05 谱来源批 §8.5 在案漂移；本批零 src 源码触碰、registry_gen 对
  基线逐字节同、当日读数与前树在案读数逐项一致，档案重生成留待独立批次）；T4/T5
  代码域零命中（全域 grep 18/1 命中均为先在停架容器字面量与文档自引，与前树计数
  零差）。状态：已实现。
- 2026-10-12：**四容器同步批（D6 裁决）**——可转尽转落地（提交 edac2d1 / 43e5d11 /
  28b8a4d / a76a4e6 / 9132bf1；实施记录 .superpowers/sdd/2026-10-12-tier1-conversion/）：
  - **转换形态（包装法，四容器同模板）**：reference / fh2_leps / h2ag111PES / h3_bkmp2
    的 box 保留全部现行公有名（判据程序与 probe 钩子零改动），另暴露五个魔数名薄包装
    pes/keys/load/term/classify，删 `reg_<sys>.f90` 落入 prebuild 自动织入（无 reg 的
    systems/ 文件夹 = Tier 1）；拉键序列以原 reg 为准逐字照抄。
  - **D6（用户裁决）**：可转尽转——四容器全转 Tier 1；两个自动（零 reg）文件夹并存
    维持 prebuild 中止（"two auto containers are installed"，drill 共存腿随迁）。
  - **R-C1（用户裁决，结局行翻转）**：Tier 1 在场即绑（auto_classify 无条件
    container_bind_classify，§1 语义）——六个无 E_SPLIT 的 reference 生产腿
    （b_layer thermalize/andersen/glo、chain_frequencies boltz/normalmode、
    surface_window）的统计行从 `outcome counts: none (no classification bound -
    skipped)` 翻为兜底通道计数 `capture=N`（32/32/4/1/1/8 逐腿）——与谱来源腿
    dissociation=4 同款语义；验收锚相应收窄（trajectories 逐字节 + stdout 除
    timing 行、该六腿再除 outcome counts 行）。
  - **R-C2（用户裁决，惰性缺省预置）**：reg 时代 bind_force 注册期的
    harm_set_constant(0.04) 预置，织入流以 keys() 包装首调（parse 前无卫兵调用点）
    惰性预置复刻；pes() 内保留守卫为手持消费者安全网；reference_load 缝期置有效
    常数后置 staged_default=.true.（防缺省覆写输入覆写的有效 K_HARM）。
  - **等价证据**：六案例 pre/post 收窄锚全等——trajectories 逐字节 12/12（b_layer
    ×3/chain×2/harmonic/radau/surface_window/谱来源四腿）+ h3_ebk COMPUTE 中止腿
    stdout+stderr 两侧同款（合计 14 项逐字节比对，Task 4 口径）；收窄 stdout 全等
    （harmonic/radau 的 scattering=8 与谱来源四腿 dissociation=4 全行保留）。判据族
    全绿（check_driver/check_facade/check_final_state/check_auto_pes/
    check_fh2_leps/check_h2ag111/check_h2ag111_pes/check_h3_bkmp2），drill 30/30，
    双编译器干净重建绿；cau 五连保持 Tier 2（rec_reg_col 记录列所迫）仅文档标注。
    状态：已实现（含文档收尾批）。
