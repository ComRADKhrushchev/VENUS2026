# 四容器 Tier 1 格式同步——实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** reference / fh2_leps / h2ag111PES / h3_bkmp2 四容器转为 Tier 1 零 reg 魔数名形态（box 加五个薄包装、保留原公有名、删 reg 文件）；cau 五连保持 Tier 2（记录列所迫）并文档标注；并存语义维持构建期中止（用户裁决 D6，drill 重排）；六案例输出 pre/post 逐字节等价。

**Architecture:** 低风险包装法——各 box 模块**保留全部现行公有入口**（判据程序 check_\*、probe 钩子 hook_\*、probe_spectrum_h2ag 全部按旧名工作，零改动），**新增**五个魔数名薄包装（pes/keys/load/term/classify），reg 文件删除使文件夹落入 prebuild 自动织入。唯一连锁：check_driver.f90 与 check_facade.f90 的编译配方含 reg_reference.f90 行，随 reference 转换适配。

**Tech Stack:** Fortran 2008（gfortran 缺省 / ifx 第二编译器）、bash（prebuild 既有机械）。

**Spec:** docs/plans/2026-10-12-auto-pes-container.md（十词表与 Tier 1 语义）；本批用户裁决 D6（可转尽转四容器；双自动文件夹维持 prebuild 中止）。勘察底稿：.superpowers/sdd/2026-10-12-tier1-conversion/survey.md（若不存在，勘察结论已内嵌本计划各任务的入口清单）。

## Global Constraints

- **六案例 pre/post 逐字节等价**是 reference 转换的验收锚（Trajectories+stdout 除 timing 行；对照对象是转换前同树实测，非陈档 results.txt）。
- 各容器**判据程序照旧绿**（check_fh2_leps / check_h2ag111 / check_h2ag111_pes / check_h3_bkmp2 / check_final_state / probe_spectrum_h2ag 不改不重编即绿；check_driver / check_facade 仅动配方中被删文件行）。
- 键缺席门控语义 = Tier 1 设计语义（在场即绑、found 旗标模块内门控）：term 包装 `term=.false.; if(<sys>_term_on()) term=<sys>_term(sta)`；classify 包装门控关闭时返回该容器分类体的兜底通道名（读各 classify 体确定；无兜底则 'long' 并在提交信息注记）——键在场流（全部入库证据）行为不变。
- reference 的 `harm_set_constant` 预置：reg 的 registration 期预置随 reg 消失；`reference_load` 本就在每次 load 时预置生效常数（reference_box.f90:108）——seam 先于一切力消费者，等价。
- 停架容器转换后仍停架（systems/inactive/ 两层深盲区）；验证用"暂存 reference → 拖入 → make → 织入 grep → 归架"循环，任何退出形态还原。
- T4/T5 grep、双编译器、E 号禁令、中文详叙提交——沿前批全部纪律。
- 工作区：worktree `.worktrees/tier1-convert`（分支 tier1-convert，基点 master@e6b45fa）。

---

### Task 1: fh2_leps 转换

**Files:**
- Modify: `systems/inactive/fh2_leps/fh2_box.f90`（加五包装）
- Delete: `systems/inactive/fh2_leps/reg_fh2_leps.f90`

**Interfaces:**
- Consumes: 现行 fh2_box 公有 API（全部保留）；reg 的 11 个 pull（R_TERM/R_HF + 9 个 surf 键）移入 load()。
- Produces: 模块 fh2_box 新增 public: pes/keys/load/term/classify（magic 名）。

- [ ] **Step 1: box 加五包装**（public 行与 contains 区追加；拉键序列照抄 reg_fh2_leps.f90:52-62）

```fortran
   public :: pes, keys, load, term, classify   ! Tier 1 magic names (thin wrappers)

   ! ---- Tier 1 magic-name wrappers (auto-weave contract; the public API above stays) ----
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      call fh2_pes(q, mass, g, v)
   end subroutine pes
   function keys() result(kk)
      character(len=16) :: kk(11)
      kk = fh2_keys()
   end function keys
   subroutine load()
      type(fh2_params_t) :: p
      call pull('R_TERM',  p%r_term,  p%r_term_found)
      call pull('R_HF',    p%r_hf,   p%r_hf_found)
      call pull('D_FH',    p%surf(1), p%surf_found(1))
      call pull('R0_FH',   p%surf(2), p%surf_found(2))
      call pull('BETA_FH', p%surf(3), p%surf_found(3))
      call pull('D_HH',    p%surf(4), p%surf_found(4))
      call pull('R0_HH',   p%surf(5), p%surf_found(5))
      call pull('BETA_HH', p%surf(6), p%surf_found(6))
      call pull('DL0_HH',  p%surf(7), p%surf_found(7))
      call pull('DL0_FH',  p%surf(8), p%surf_found(8))
      call pull('AMOD',    p%surf(9), p%surf_found(9))
      call fh2_load(p)
   end subroutine load
   logical function term(sta)
      type(state_t), intent(in) :: sta
      term = .false.
      if (fh2_term_on()) term = fh2_term(sta)
   end function term
   function classify(fin) result(name)
      type(fin_t), intent(in) :: fin
      character(len=32) :: name
      if (fh2_class_on()) then
         name = fh2_classify(fin)
      else
         name = '<fh2 兜底通道——读 fh2_classify 体确定, 无则 long>'
      end if
   end function classify
```

注意：box 需补 `use input, only: pull`（若未有）；state_t/fin_t 已在 use 域。classify 兜底通道由实现者读 fh2_classify（fh2_box.f90:161）后填定并在报告注明。

- [ ] **Step 2: 删 reg 文件**：`git rm systems/inactive/fh2_leps/reg_fh2_leps.f90`
- [ ] **Step 3: 验证**
  - check_fh2_leps 重编重跑（pot 未动，配方不变）→ 全绿；
  - 拖入循环：暂存 reference（mv systems/reference systems/inactive/reference）→ mv systems/inactive/fh2_leps systems/fh2_leps → cd src && make clean && make（EXIT=0）→ 织入 grep：`use fh2_box, only: pes, keys, load, term, classify` 与 `container_sched_bind(keys = auto_keys, &` 命中 → 归架还原 → make 绿；
  - `bash scripts/prebuild.sh . build/gfortran`（还原态）确认 fh2 不在织入（inactive 盲区）。
- [ ] **Step 4: Commit**（中文详叙：包装法/拉键迁移/门控语义/验证三腿）

### Task 2: h2ag111PES 转换

同 Task 1 模式。`h2ag111_box.f90` 加五包装：load() 拉 Z_TERM/R_H2（照抄 reg_h2ag111.f90:60-62）后调 h2ag111_load(p)（文件检查+pes_init 原地不动——NN 权重挂载位本就在这）；term/classify 门控包装（兜底通道读 h2ag111_classify 体确定）。删 reg_h2ag111.f90。验证：check_h2ag111 + check_h2ag111_pes 重编重跑绿 + probe_spectrum_h2ag 重编重跑绿（旧名 API 未动）+ 拖入构建织入 grep。h2ag111PES 文件夹名含大写——文件夹名任意，不影响。

### Task 3: h3_bkmp2 转换

同模式。`bkmp2_box.f90` 加五包装：load() 拉 R_TERM/R_H2（照抄 reg_h3_bkmp2.f90:58-60）后调 bkmp2_load(p)；term/classify 门控包装（classify 读 fin%q_raw 的既有体不动，兜底通道读体确定）。删 reg_h3_bkmp2.f90。验证：check_h3_bkmp2 重编重跑绿 + 拖入构建织入 grep。

### Task 4: reference 转换 + 判据配方适配 + 六案例 pre/post

**Files:**
- Modify: `systems/reference/reference_box.f90`（加五包装）
- Delete: `systems/reference/reg_reference.f90`
- Modify: `scripts/check_driver.f90`、`scripts/check_facade.f90`（配方中 reg_reference 行的处理——见 Step 3）

- [ ] **Step 1: PRE 快照**——六案例在**当前**二进制下各跑一次，trajectories.txt 与 stdout（除 timing 行）存 `.superpowers/sdd/2026-10-12-tier1-conversion/pre_case_<name>/`。
- [ ] **Step 2: box 加五包装 + 删 reg**——load() 拉 K_HARM/Q_TERM/E_SPLIT（照抄 reg_reference.f90:56-60）后调 reference_load(p)；term 包装直调 reference_term（其体已自带 term_on 短路——reference_box.f90:131，勘察确认）；classify 门控包装（reference_class_on 关闭时兜底通道读 reference_classify 体确定）。`use input, only: pull` 补入。git rm reg_reference.f90。
- [ ] **Step 3: 判据配方适配**——先 grep 两判据对 `registry_gen\|reg_reference_slots` 的**使用**（不止编译行）：若仅编译未调用 → 删配方中 reg_reference.f90 与 registry_gen.f90 行；若有调用（如 check_driver 走真实装配序）→ 以最小改动改为直接绑定序（细节在实现报告注明理由）。check_final_state 配方不含 reg（勘察确认）不动。
- [ ] **Step 4: 重建 + POST 比对**——make clean && make（reference 现为唯一自动容器，织入生效；EXIT=0）→ 六案例重跑 → 与 Step 1 快照逐字节 cmp（trajectories 全等 + stdout 除 timing 全等）。任何一字节差异 = 阻断，停下诊断。
- [ ] **Step 5: 判据族**——check_driver / check_facade / check_final_state / check_auto_pes 重编重跑绿。
- [ ] **Step 6: Commit**（含 pre/post 逐字节等价证据链）。

### Task 5: drill 重排（D6 并存语义）

**Files:** Modify: `scripts/drill_auto_pes.sh`

- 编排改为：leg 0 暂存 reference（mv systems/reference systems/inactive/reference）→ leg 1 make clean + 拖入 auto_sample + 单次 make（原首编腿语义保留）→ …… 运行腿照旧（reference 已暂存）→ 原"运行期 STOP 1 并存腿"（腿 3）**替换**为：拖入态下把 reference 搬回 → 直接跑 `bash ../scripts/prebuild.sh . build/gfortran` 断言 exit≠1... 期望 exit=1 且 stderr 含 `two auto containers are installed`（构建期中止，D6）→ 再暂存 reference 继续后续腿 → 负面腿族/零触碰照旧（restore_tree 已含 reference 归位）。腿数维持 30（腿 3 语义替换）；PATS 与断言串相应更新；重跑全绿。
- Commit。

### Task 6: 文档 + 全量回归收尾

- EXTENDING.md 1b：reg 骨架范本指针改 `systems/inactive/cau_rst/reg_cau_rst.f90`（唯一存续 reg 族范本）；1a 的"现成范本"段补三物理容器已 Tier 1 化一句；§1b 补 cau 五连"Tier 2（自定义记录列 rec_reg_col 所迫）"标注。
- systems/inactive/README.md：fh2_leps/h2ag111PES/h3_bkmp2 条目改"Tier 1 魔数名（无 reg）"；cau 五连标注 Tier 2 原因。
- src/README.md 的 reg 结构段补一句"或 Tier 1 零 reg（EXTENDING §1a）"。
- 回归：T4/T5 grep（src/+systems/，代码域零命中）；双编译器干净重建；六案例对 Task 4 pre 快照复比；drill 全绿；check 全族绿。
- docs/plans/2026-10-12-auto-pes-container.md §6 追加"四容器同步批（D6 裁决）"段。
- 终批提交。

---

## Self-Review

1. **覆盖**：四容器（T1-T4）+ D6 drill（T5）+ 文档与回归（T6）= 用户裁决 D6 全范围；cau 族仅文档标注（T6）。无缺口。
2. **占位符**：三处"兜底通道读体确定"是有意留白（实现者读代码后填定并注明——非 TBD 式含糊，是"由被读代码决定"的明确指令）；其余代码块完整。
3. **一致性**：包装法统一五容器形态（含 auto_sample 先例）；拉键序列逐一取自勘察的 reg 行号；pre/post 锚只属 Task 4（其余容器无生产流）。
