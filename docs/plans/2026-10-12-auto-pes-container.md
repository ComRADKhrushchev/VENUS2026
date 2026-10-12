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

设计分节（契约 / 机械 / 缺省与升级 / 文档判据）均已逐节获用户同意（2026-10-12）。

## 1. Tier 1 纯势容器——五名词表（封闭集合）

一个 systems/ 顶层文件夹，**无 `reg_*.f90`**，即纯势容器候选。作者面 = 文件夹内
一个（或多个）.f90，暴露下列**魔数名**中的若干个（名字即约定，签名即 E14/E19/E27
既有接口，编译期核对）：

| 魔数名 | 形 | 必备性 | 签名要点 | 织入槽位 |
|---|---|---|---|---|
| `pes` | subroutine | **恰一个** | `(q, mass, g, v)`：q[A]/mass[amu] 进，g[eV/A]/v[eV] 全赋值出（E27 力交付契约原样） | force（registration） |
| `keys` | function | 可选 0/1 | `character(len=16)` 数组（`reference_keys` 同型），容器键词表 | keys（parse 前） |
| `load` | subroutine | 可选 0/1 | 无参；模块**自持 `pull`**（`use input, only: pull`），缺键保缺省 | load（seam） |
| `term` | function | 可选 0/1 | `term_i` 同型：`(sta) → logical` | term（seam） |
| `classify` | function | 可选 0/1 | `classify_i` 同型：`(fin) → name(len=32)` | classify（seam） |

规则：

- **词表封闭**：不再增名。带注册词汇的扩展（`method_reg(name,prop,match)`、
  `rec_reg_col`、arm 武装、方法包）不入 Tier 1——那是 Tier 2 全容器的事。
  边界原则一句话：**Tier 1 = 五个无参槽位体**。
- **在场即绑（无条件）**；键缺席门控（如 R_TERM 缺席 → 跑满 MAX_STEPS）由模块
  内部 found 旗标实现（范本示例见 §3）。
- **名字可居文件夹内任意模块**（检测按名找归属模块；通常全在一个模块）。
- **公共性**：入口须 public——织入体 `use <mod>, only: <名>` 编译期核对，
  private 即编译错（文档明示预期错误形态）。
- **reg 权威**：文件夹有 `reg_*.f90` → 五名全部忽略，一切照旧（自动接线静默退出）。

### 有名中止（prebuild 期，全部指名）

1. 无 reg 且 `pes` 计数 0 —— "systems/<X> 已安装（文件夹在场）但未接线：暴露一个
   `pes` 入口，或写 `reg_<X>.f90` 走全容器"（比今天"编译过、首次力调用才死"更早更响）；
2. 任一名计数 ≥2 —— 指名两处定义；
3. 名命中但形不符（如 `subroutine keys`）—— 指名（把安静失败变有声）；
4. 自动文件夹 ≥2 —— 一进程一容器（2026-09-21 裁决）的构建期前置；
5. 自动容器与 reference（或其他 reg 容器）并存激活 → 既有排程第二填 STOP 1
   原样生效（测试自动容器仍需暂存 reference，机制不变）。

### 零新增监管（既有兜底照管一致性问题）

- `load` 无 `keys`：合法空转（封闭语法本就不放行未声明键，pull 落空保缺省）；
- `keys` 无 `load`：既有 `input_audit_declared` 兜底中止（buffered 而无人消费）；
- 签名不符：编译期接口核对即错；
- 撤掉 `pes`：prebuild 中止 1。

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
   use mysys_pot, only: pes, keys, load, term, classify                 ! 按检测到的名列表；
                                                                         ! 名跨多个模块时逐模块各一行 use
   ...
contains
   subroutine auto_slots()                       ! 织入的纯势手（只填在场槽）
      call container_sched_bind(force    = auto_force, &
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
| 电子方法 | 无，绝热参考（ELEC_METHOD 缺省） | method/arm 槽空 |
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
   subroutine pes(q, mass, g, v)                          ! 必备
      real(8), intent(in)  :: q(:,:), mass(:)
      real(8), intent(out) :: g(:,:), v
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

### 升级 = 写 reg，pot 五函数全部直接复用

哪天需要方法包/记录列/更复杂的参数体了：同文件夹加 `reg_<sys>.f90`（+ 可选
box）——prebuild 见 reg 即静默退出自动接线、reg 权威接管编排。**pot 一字不动**：
reg 区体直接 `use mysys_pot, only: pes, keys, load, term` 引用同一批函数
（`container_bind_pes(pes)`、`input_declare_keys(keys())`……）。所谓 box 在
Tier 1 里被压扁成了 pot 模块自己的 `keys`/`load` 两函数。

## 4. 非目标与登记后续件

- **Python/脚本势桥接**（D1 后置）：子进程管道协议 / i-PI socket 方案本轮只留缝
  （五名词表封闭、外部势另立容器形态另案裁决）。
- **能量 only 入口 + 框架差分**：D3 否决，不做。
- **通用距离终止键**（TERM_R 类，程序级、范式感知——气对按碎片间距/表面按弹射体
  高度）：真正的新裁决（触碰 E14/E19 容器属位），登记为后续件候选，不混本批。
- **激活机制改革**（SYSTEM=xxx 输入键选择、多容器共存常驻编译）：D2 明确不动。

## 5. 判据与验收

1. **判据程序** `scripts/check_auto_pes.f90` + 演练脚本（判据族同款）：
   - 正面：范本容器拖入 systems/ → **首次 make 即绿**（首编即绿核验点）→
     `build/<FC>/registry_gen.f90` 织入内容断言（use 行、auto_slots、调用行）→
     MAX_STEPS 安全网跑通一条；`keys`/`load`/`term`（键在场与缺席两分支）/
     `classify` 路径各覆盖一条；
   - 负面：pes 0 / pes ≥2 / 任一名 ≥2 / 双自动文件夹 / 形不符（subroutine keys），
     中止文案逐条断言；
   - T2 式零触碰演练：拖入 → 判据绿 → 拖出重建 → git tracked 零改动。
2. **回归**：五案例 results.txt 逐字节一致；双编译器全量重建绿；T4/T5 grep 复验
   （23.0605 仅 force_interface；RANDOM_NUMBER 零命中）。
3. **文档落点**：EXTENDING.md §1 重写为 Tier 1 领篇（五名词表 + 范本 + 升级注记），
   现行 8 槽三件套降为 Tier 2 全容器；README.md / src/CLAUDE.md / AGENTS.md
   注册段同步增补；本文件即裁决记录。

## 6. 状态

- 2026-10-12：设计定稿，用户逐节批准（D1-D4 裁决链 + 四节设计）。待实现。
