# reg_ 型容器支持退役 + 第 11 预留名 observables —— 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 解除容器 reg_*.f90 形态的支持（静默忽略式），以第 11 预留名 `observables` 补齐 Tier 1 的自定义记录列能力，并把 cau 五连转为 Tier 1 包装法形态。

**Architecture:** prebuild 的预留名扫描扩至十一名（新增 observables 子程序名，织入到排程既有 columns 槽）；每容器新增独立 `*_wires.f90` 承载全部薄包装（既有文件与判据配方零触碰）；删除 prebuild 的"reg 让位"分支实现静默忽略。

**Tech Stack:** bash(prebuild.sh) / Fortran(systems/*) / pytest(gui/) / 判据程序(gfortran standalone)

**Spec:** docs/plans/2026-10-14-reg-retirement-observables.md（裁决记录——本计划从该文件论证，执行者两者都读）

## Global Constraints

- 三项用户裁决（2026-10-14）：第 11 预留名定名 `observables`；reg 解除=静默忽略（无警告无 loud-stop，reg 文件照编为无害陪葬）；验收锚=五判据前后全绿且读数不变 + cau_tab 归档逐位。
- 判据基线（2026-10-13 07fb8ce 批落盘）：check_cau_tab 150/0、cau_peri 103/0、cau_emt2 99/0、cau_emtnn 94/0、cau_rst 94/0（logs/sync_check_cau_*.log）。
- methods/ 家族 reg_*（reg_samp_*/reg_inc_*/reg_bath_*/reg_adiabatic）完全不动——那是方法注册机制，与容器 reg 是两回事。
- 冒烟正例：publish_min/cases/_smoke_cau_tab_ev/（归档 trajectories.txt + pes_table.txt + system/）+ 迁移输入 build/cau_smoke_mig/input_qct.txt。
- 所有验证读数落盘 logs/（用户纪律 2026-10-04）。
- 激活测试的架规：暂存 reference（mv systems/reference systems/inactive/reference）、测毕恢复、make 两次、恢复后重建。

---

### Task 1: prebuild 第 11 预留名 observables + auto_sample 样例 + drill 断言

**Files:**
- Modify: `scripts/prebuild.sh:118-119`（扫描名集）、`:142-143`（名循环与形态）、`:167` 前无改动、`:206`（use 行收集序）、`:229-272`（auto_slots args 与包装发射）
- Modify: `systems/inactive/auto_sample/io.f90`（加 observables 正例）
- Modify: `scripts/drill_auto_pes.sh:61-71`（PATS 织入断言）

**Interfaces:**
- Produces: 预留名词表十一名（pes/keys/load/term/classify + prop/match/sync/init/method_name + observables）；织入发射新增 `columns = auto_columns` 关键字形与 `auto_columns()` 包装（体内 `call observables()`）。
- 消费: 排程既有 columns 槽（container_sched_bind 的 columns 虚参，reg 时代即存在，零改动）。

- [ ] **Step 1: drill 断言先行（红）**

`scripts/drill_auto_pes.sh` 的 PATS 块（约 61-71 行）改两处：

```
use auto_sample_io, only: load, classify
```
改为
```
use auto_sample_io, only: load, classify, observables
```
并在 `force = auto_force, &` 行后新增一行 PAT：
```
columns = auto_columns
```

- [ ] **Step 2: 样例加 observables 正例**

`systems/inactive/auto_sample/io.f90`：public 列表加 `observables`；模块 use 加 `use recorder, only: rec_reg_col`（io.f90 是依赖隔离层，加 recorder 不破坏其独立编译性——drill 整树编译，recorder 在树内）；contains 尾部加：

```fortran
   !------------------------------------------------------------------
   ! observables() - the reserved-name column registration (recorder
   !                 columns the run opts into via OBSERVABLES)
   !------------------------------------------------------------------
   subroutine observables()
      call rec_reg_col('stub_count', stub_count)
   end subroutine observables

   real(8) function stub_count(sta)
      type(state_t), intent(in) :: sta
      stub_count = 1.0d0
   end function stub_count
```

（io.f90 若尚未 use state 则补 `use state, only: state_t`。）

- [ ] **Step 3: prebuild 三处名集 + 织入发射**

`scripts/prebuild.sh`：
1. 行 118：`split("pes load prop match sync init", sa)` → `split("pes load prop match sync init observables", sa)`
2. 行 142：`for nm in pes keys load term classify prop match sync init method_name; do` → 尾部加 ` observables`
3. 行 143：`case "$nm" in pes|load|prop|match|sync|init) want=s;;` → `pes|load|prop|match|sync|init|observables)`
4. 行 206（use 行收集循环）与 142 同步加 ` observables`
5. auto_slots args 段（行 229-236），在 classify 行后加：
```bash
      [ -n "${auto_mod[observables]:-}" ] && args="$args columns"
```
6. 包装发射段（行 270-272 classify 包装后）加：
```bash
      [ -n "${auto_mod[observables]:-}" ] && { echo "   subroutine auto_columns()"
         echo "      call observables()"
         echo "   end subroutine auto_columns"; }
```

- [ ] **Step 4: 内环快验（免整 drill）**

```bash
cd src && bash ../scripts/prebuild.sh . build/gfortran && grep -n "observables\|columns = auto_columns" build/gfortran/registry_gen.f90 | head -5
```
预期：`use auto_sample_io, only: load, classify, observables` 与 `columns = auto_columns` 各至少 1 命中。（此时 auto_sample 在 inactive 架，prebuild 不扫——把样例临时拖到 systems/ 顶层再跑，验毕归架；reference 需同时暂存避免双 auto 中止。）

- [ ] **Step 5: 全量 drill（绿）**

```bash
bash scripts/drill_auto_pes.sh > logs/regret_drill_t1.log 2>&1; echo rc=$?
```
预期：`drill summary: pass=31 fail=0`（原 30 腿 + 新 columns 织入断言 1；use 行修改不增腿数）。

- [ ] **Step 6: Commit**

```bash
git add scripts/prebuild.sh scripts/drill_auto_pes.sh systems/inactive/auto_sample/io.f90
git commit -m "第 11 预留名 observables(2026-10-14 裁决定名): prebuild 扫描/校验/use 行收集/织入发射四处扩名, auto_slots 绑排程既有 columns 槽(auto_columns 包装直调作者 observables 体); auto_sample 样例 io.f90 增 observables 正例(rec_reg_col stub_count); drill PATS +1 断言(columns = auto_columns)与 use 行扩名。TDD: drill 断言先行红→实现后 31/0(logs/regret_drill_t1.log)"
```

---

### Task 2: cau_tab 转 Tier 1（wires 文件 + 删 reg）

**Files:**
- Create: `systems/inactive/cau_tab/cau_tab_wires.f90`
- Delete: `systems/inactive/cau_tab/reg_cau_tab.f90`

**Interfaces:**
- Consumes: Task 1 的 observables 预留名；cau_tab_box/cau_tab_tdhf 既有公有名（全部保持不动——判据与 probe 钩子零改动）。
- Produces: cau_tab 文件夹以十名（无 classify）+ observables 共十一名中的十个形态被自动接线（pes/keys/load/term/prop/match/sync/method_name/init/observables）。

**关键设计（比 reference 先例更彻底的零触碰）**：薄包装**不进 box**，独立 `cau_tab_wires.f90` 承载——box 若加 `use input/recorder/config` 会断五判据配方（配方不含这些模块）。wires 文件只被构建通配与 prebuild 扫描消费，判据配方显式列文件、不含它。

- [ ] **Step 1: 写 cau_tab_wires.f90（内容=reg 七体逐字搬运）**

```fortran
!=====================================================================
! cau_tab_wires.f90 - the Tier 1 reserved-name wrappers (the folder's
!   auto-weave contract: prebuild weaves these names into registry_gen
!   when the folder sits reg-less at systems/ top level; each delegates
!   to the public API the check programs and probe hooks keep using).
!   Replaces reg_cau_tab.f90 (the reg form retired 2026-10-14,
!   docs/plans/2026-10-14-reg-retirement-observables.md).
! Units: none
!=====================================================================
module cau_tab_wires
   use state,          only: state_t, elec_t
   use config,         only: electronic
   use input,          only: input_declare_keys, pull
   use recorder,       only: rec_reg_col
   use cau_tab_box,    only: cau_tab_params_t, cau_tab_keys, cau_tab_load, &
                             cau_tab_term, cau_tab_term_on, &
                             cau_tab_pes_r, cau_tab_prop_r, &
                             cau_pars_start, cau_pars_norb, cau_pars_bw, cau_pars_telec, &
                             cau_pars_msbstp, cau_pars_nsbstp, cau_pars_pmode, cau_pars_pamp
   use cau_tab_tdhf,   only: cau_tab_match, cau_tab_tdhf_setup, cau_tab_sync, &
                             cau_tab_col_nup, cau_tab_col_ndn, cau_tab_col_mag, &
                             cau_tab_col_eelec, cau_tab_col_emf, &
                             cau_tab_col_eimp, cau_tab_col_ebath, cau_tab_col_ehyb
   implicit none
   private
   public :: pes, keys, load, term
   public :: prop, match, sync, method_name, init, observables
contains
   !------------------------------------------------------------------
   ! pure-PES four
   !------------------------------------------------------------------
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      call cau_tab_pes_r(q, mass, g, v)   ! the roster-order boundary wrapper
   end subroutine pes

   function keys() result(kk)
      character(len=16) :: kk(13)
      kk = cau_tab_keys()
   end function keys

   subroutine load()
      type(cau_tab_params_t) :: p
      call pull('START_U0',     p%start)
      call pull('N_ORB',        p%norb)
      call pull('BAND_WIDTH',   p%bw)
      call pull('TELECTRO',     p%telec)
      call pull('MSBSTP',       p%msbstp)
      call pull('NSBSTP',       p%nsbstp)
      call pull('PERTURB_MODE', p%pmode)
      call pull('PERTURB_AMP',  p%pamp)
      call pull('K_SURF',       p%ksurf)
      call pull('K_GHOST',      p%kghost)
      call pull('Q0_GHOST',     p%q0ghost)
      call pull('TABLE_FILE',   p%table_file)
      call pull('Z_TERM',       p%z_term, p%z_term_found)
      call cau_tab_load(p)
   end subroutine load

   logical function term(sta)
      type(state_t), intent(in) :: sta
      if (.not. cau_tab_term_on()) then   ! the reg-era bind gate, carried
         term = .false.                   ! by the wrapper (the body does
         return                           ! not self-gate; an unbound run
      end if                              ! must answer .false.)
      term = cau_tab_term(sta)
   end function term

   !------------------------------------------------------------------
   ! the method package (tdhf_hub)
   !------------------------------------------------------------------
   subroutine prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
      call cau_tab_prop_r(sta, dt)
   end subroutine prop

   subroutine match(sta)
      type(state_t), intent(inout) :: sta
      call cau_tab_match(sta)
   end subroutine match

   subroutine sync(s)
      type(elec_t), intent(inout) :: s
      call cau_tab_sync(s)
   end subroutine sync

   function method_name() result(s)
      character(len=32) :: s
      s = 'tdhf_hub'
   end function method_name

   subroutine init()
      if (trim(electronic%method) == 'tdhf_hub') then
         call cau_tab_tdhf_setup(cau_pars_norb, cau_pars_bw, cau_pars_telec, &
                                 cau_pars_start, cau_pars_msbstp, cau_pars_nsbstp, &
                                 cau_pars_pmode, cau_pars_pamp)
      end if
   end subroutine init

   !------------------------------------------------------------------
   ! the record columns (the 11th reserved name, 2026-10-14)
   !------------------------------------------------------------------
   subroutine observables()
      if (trim(electronic%method) == 'tdhf_hub') then
         call rec_reg_col('n_imp_up', cau_tab_col_nup)
         call rec_reg_col('n_imp_dn', cau_tab_col_ndn)
         call rec_reg_col('mag_mom', cau_tab_col_mag)
         call rec_reg_col('e_elec_ev', cau_tab_col_eelec)
         call rec_reg_col('e_mf_ev', cau_tab_col_emf)
         call rec_reg_col('e_imp_ev', cau_tab_col_eimp)
         call rec_reg_col('e_bath_ev', cau_tab_col_ebath)
         call rec_reg_col('e_hyb_ev', cau_tab_col_ehyb)
      end if
   end subroutine observables
end module cau_tab_wires
```

**签名核对（写前必查）**：`cau_tab_keys()` 返回的字符长度与元素数以上例 `character(len=16) :: kk(13)` 为占位——执行时以 cau_tab_box.f90 的实际声明为准照抄；`cau_tab_prop_r`/`cau_tab_match` 的哑元签名以各自模块声明为准。

- [ ] **Step 2: 删 reg**

```bash
git rm systems/inactive/cau_tab/reg_cau_tab.f90
```

- [ ] **Step 3: 判据前后对照**

```bash
cd scripts && gfortran -O2 -o check_cau_tab.exe ../src/globals/consts.f90 ../src/globals/state.f90 ../src/utils/linalg.f90 ../src/utils/rng.f90 ../systems/inactive/cau_peri/surface_potentials.f90 ../systems/inactive/cau_peri/cau_peri_pot.f90 ../systems/inactive/cau_tab/cau_tab_tab.f90 ../systems/inactive/cau_tab/cau_tab_pot.f90 ../systems/inactive/cau_tab/cau_tab_tdhf.f90 ../systems/inactive/cau_tab/cau_tab_box.f90 check_cau_tab.f90 && ./check_cau_tab.exe > ../logs/regret_check_cau_tab.log 2>&1; echo rc=$?
diff ../logs/sync_check_cau_tab.log ../logs/regret_check_cau_tab.log && echo IDENTICAL
```
预期：`150 PASS, 0 FAIL` 且与基线逐字一致（判据不编译 wires/reg，恒等性由包装法结构保证）。

- [ ] **Step 4: 激活冒烟（架规）+ observables 全列跑**

```bash
cp -r systems/inactive/cau_tab systems/cau_tab && mv systems/reference systems/inactive/reference
cd src && make > /tmp/m1.log 2>&1; make > /tmp/m2.log 2>&1; echo "make2 rc=$?"
cd .. && RUN=build/regret_tab && rm -rf $RUN && mkdir -p $RUN && \
  cp publish_min/cases/_smoke_cau_tab_ev/pes_table.txt $RUN/ && \
  cp -r publish_min/cases/_smoke_cau_tab_ev/system $RUN/ && \
  cp build/cau_smoke_mig/input_qct.txt $RUN/ && cd $RUN && \
  ../../src/build/gfortran/venus.e > stdout.txt 2> stderr.txt; echo rc=$?
cd ../.. && cmp publish_min/cases/_smoke_cau_tab_ev/trajectories.txt build/regret_tab/trajectories.txt && echo BITWISE
```
再跑列健全性：复制 build/cau_sync_cols/input_qct.txt（OBSERVABLES 含 e_imp_ev/e_bath_ev/e_hyb_ev 的那版）到新 run 目录同法跑，预期 rc=0、`columns=11`。
两份 stdout 落盘 `logs/regret_tab_smoke.log`、`logs/regret_tab_cols.log`。
恢复：`rm -rf systems/cau_tab && mv systems/inactive/reference systems/reference && cd src && make`。

- [ ] **Step 5: Commit**

```bash
git add systems/inactive/cau_tab && git commit -m "cau_tab 转 Tier 1(包装法, wires 文件形态): 新增 cau_tab_wires.f90 承载十名薄包装(pes/keys/load/term+prop/match/sync/method_name/init+observables)——独立文件而非进 box, 保五判据配方与既有文件零触碰(box 加 use input/recorder/config 会断判据); term/init/observables 包装承接 reg 时代的 ELEC_METHOD/term_on 门控; 删 reg_cau_tab.f90。验证: check_cau_tab 150/0 与 07fb8ce 基线逐字一致(logs/regret_check_cau_tab.log); 激活冒烟 trajectories 与归档逐位(logs/regret_tab_smoke.log); observables 全列 rc=0 columns=11(logs/regret_tab_cols.log)"
```

---

### Task 3: 四兄弟转 Tier 1（同模板）

**Files:**
- Create: `systems/inactive/cau_peri/cau_peri_wires.f90`、`cau_emt2/cau_emt2_wires.f90`、`cau_emtnn/cau_emtnn_wires.f90`、`cau_rst/cau_rst_wires.f90`
- Delete: 对应四份 `reg_*.f90`

**Interfaces:** 同 Task 2（模板即 cau_tab_wires.f90）。

- [ ] **Step 1: 按差异表生成四份 wires**

模板 = Task 2 的 cau_tab_wires.f90，逐容器替换以下要素（全部从各自 reg 文件逐字搬运，reg 删除前它是唯一权威源）：

| 要素 | cau_peri | cau_emt2 | cau_emtnn | cau_rst |
|---|---|---|---|---|
| 模块前缀 | cau_peri | cau_emt2 | cau_emtnn | cau_rst |
| pes 目标 | cau_peri_pes_r | cau_emt2_pes | cau_emtnn_pes | cau_rst_pes |
| prop 目标 | cau_peri_prop_r | cau_emt2_prop | cau_emtnn_prop | cau_rst_prop |
| load 的 pull 集 | 13 键: START_U0/N_ORB/BAND_WIDTH/TELECTRO/MSBSTP/PERTURB_MODE/PERTURB_AMP/K_SURF/K_GHOST/Q0_GHOST/USE_NN/NN_FILE/Z_TERM(+found) | 12 键: …PERTURB_AMP/EMT_FILE/BOXLX/BOXLY/SKEW/Z_TERM(+found) | 12 键: …PERTURB_AMP/SLAB_FILE/WEIGHTS_FILE/SLAB_A_LAT/N_SIDE/Z_TERM(+found) | 11 键: …PERTURB_AMP/SLAB_FILE/WEIGHTS_FILE/BVK_FILE/Z_TERM(+found) |
| init 体 setup 实参 | 照抄各自 reg 的 init_method 体（cau_pars_* 列表各家不同） | 同左 | 同左 | 同左 |
| observables 体 | 6 列(n_imp_up/n_imp_dn/mag_mom/e_elec_ev/e_mf_ev) | 6 列同左 | 6 列同左 | 6 列同左 |
| method_name | 'tdhf_hub' | 'tdhf_hub' | 'tdhf_hub' | 'tdhf_hub' |

use 段按各家实际引用的目标名+cau_pars_* 集从 reg 的 use 块照抄（前缀替换）。

- [ ] **Step 2: 删四份 reg + 判据前后对照**

```bash
git rm systems/inactive/cau_peri/reg_cau_peri.f90 systems/inactive/cau_emt2/reg_cau_emt2.f90 systems/inactive/cau_emtnn/reg_cau_emtnn.f90 systems/inactive/cau_rst/reg_cau_rst.f90
```
对四家各跑判据（配方同 07fb8ce 批的 logs/sync_* 对应命令），落盘 `logs/regret_check_cau_{peri,emt2,emtnn,rst}.log`，各与基线 diff 预期逐字一致（103/99/94/94 PASS）。

- [ ] **Step 3: Commit**

```bash
git add systems/inactive/cau_peri systems/inactive/cau_emt2 systems/inactive/cau_emtnn systems/inactive/cau_rst && git commit -m "cau 四兄弟转 Tier 1(同 cau_tab 模板 wires 形态): peri/emt2/emtnn/rst 各增 *_wires.f90 十名薄包装, 删四份 reg; 差异要素(pes/prop 目标两派/各家 pull 集与 init 实参/6 列集)逐字搬运自各自 reg。验证: 四判据与 07fb8ce 基线逐字一致(logs/regret_check_cau_*.log)"
```

---

### Task 4: prebuild reg 静默忽略 + 中止文案 + 周边同步

**Files:**
- Modify: `scripts/prebuild.sh:100`（删 reg 让位行）、`:167`（中止文案）、`:188-189`（注释）
- Modify: `src/EXTENDING.md`（删 reg 骨架节，Tier 2 史层注记）
- Modify: `AGENTS.md`（并存句收窄、词表十一名、容器架条目）
- Modify: `docs/glossary.md`（Tier 2 条目退役标记、预留名条目补 observables）
- Modify: `gui/core/pes_scan.py:45,156`（删 Tier2 分支与注记）、`gui/pages/pes_page.py:243`（删显式注册容器 info 行）
- Modify: `gui/tests/test_pes_scan.py:65-68`（test_reg_folder_is_tier2 随语义改）

**Interfaces:** 无代码接口变化——纯退役与文案。

- [ ] **Step 1: prebuild 删 reg 让位分支**

删除行 100：`ls "$d"reg_*.f90 >/dev/null 2>&1 && continue`。
行 167 中止文案：`... is installed (folder present, no reg_*.f90) but wires nothing: expose exactly one 'pes' entry, or write a reg_*.f90 full container` → `... is installed (folder present) but wires nothing: expose exactly one 'pes' entry (the reserved names)`（**保留前缀** `wires nothing: expose exactly one 'pes'`——drill negative-1 grep 的是这个固定子串）。
行 188-189 注释中 "a reg-less systems/ folder" → "a systems/ folder"。

- [ ] **Step 2: 静默忽略语义验证**

把一份**未转换**的旧 reg 容器文件夹临时拖入 systems/（用 publish_min/systems/cau_tab——它仍是 reg 形态）：
```bash
cp -r publish_min/systems/cau_tab systems/cau_tab_oldreg
(cd src && bash ../scripts/prebuild.sh . build/gfortran) > /tmp/silent.log 2>&1; echo rc=$?
grep "wires nothing" /tmp/silent.log | head -2
```
预期：rc≠0 且中止文案为预留名版本（旧 reg 容器无预留名 → 既有 named abort，静默忽略的失败模式）。验毕 `rm -rf systems/cau_tab_oldreg`。

- [ ] **Step 3: 文档与 GUI 同步**

- EXTENDING.md：删 reg 骨架节（§1b），其位置放一行史层注记（"reg_*.f90 全容器形态 2026-10-14 退役——裁决记录 docs/plans/2026-10-14-reg-retirement-observables.md"）。
- AGENTS.md：①命名契约段 "a systems/ (container) reg hosts `<mod>_slots` only" 条目改史层注记；②容器架条目 "auto+reg 并存仍装配期 STOP 1" 分支删（reg 已不存在）；③预留名句 "ten reserved names" → "eleven reserved names (… + observables)"。
- glossary.md：Tier 2 全容器行加 "（2026-10-14 退役，史层）"；预留名行的十名列举补 observables。
- pes_scan.py:45 注释与 :156 的 reg 分支删（is_tier1 语义: 含 reg 不再是特殊形态——reg 文件按普通文件对待，文件夹身份由预留名决定）；pes_page.py:243 的 info 行删；test_pes_scan.py 的 test_reg_folder_is_tier2 改为断言新语义（reg 文件在场不再改变判词）或删除该测试并在提交说明。

- [ ] **Step 4: GUI 套件 + drill**

```bash
cd gui && QT_QPA_PLATFORM=offscreen python -m pytest tests -q > logs/pytest_regret_t4.log 2>&1; tail -2 logs/pytest_regret_t4.log
cd .. && bash scripts/drill_auto_pes.sh > logs/regret_drill_t4.log 2>&1; echo rc=$?
```
预期：套件全绿；drill 31/0。

- [ ] **Step 5: Commit**

```bash
git add scripts/prebuild.sh src/EXTENDING.md AGENTS.md docs/glossary.md gui/core/pes_scan.py gui/pages/pes_page.py gui/tests/test_pes_scan.py && git commit -m "reg_ 型支持退役(2026-10-14 用户裁决: 静默忽略式): prebuild 删 reg 让位分支(行 100)——有无 reg 一律走预留名扫描, reg_*.f90 落源文件通配照编为无害陪葬; wires-nothing 中止文案去 reg 提法(保留 drill grep 子串); EXTENDING 删 reg 骨架节+史层注记; AGENTS 并存句收窄+词表十一名; glossary Tier 2 退役标记; GUI pes_scan/pes_page 删显式注册分支与判词。验证: 旧 reg 文件夹拖入→预留名中止实证(/tmp/silent.log 复跑于 logs/regret_silent_ignore.log); GUI 套件全绿(logs/pytest_regret_t4.log); drill 31/0(logs/regret_drill_t4.log)"
```

---

### Task 5: 终局验收门

**Files:** 无新改动——全库终验。

- [ ] **Step 1: 五案例回归**

```bash
for c in case_harmonic_scatter case_radau_conservation case_b_layer/thermalize case_b_layer/andersen case_b_layer/glo case_chain_frequencies/boltz case_chain_frequencies/normalmode case_surface_window; do (cd cases/$c && ../../src/build/gfortran/venus.e > /dev/null 2>&1; echo "$c rc=$?"); done | tee logs/regret_cases.log
```
预期：八行 rc=0（reference 未动，全绿预期）。

- [ ] **Step 2: 双编译器 + 终局 drill**

```bash
cd src && make > ../logs/regret_make_gf.log 2>&1; echo gf=$?
make FC=ifx > ../logs/regret_make_ifx.log 2>&1; echo ifx=$?
cd .. && bash scripts/drill_auto_pes.sh > logs/regret_drill_final.log 2>&1; echo rc=$?
```
预期：gf=0、ifx=0、drill 31/0。

- [ ] **Step 3: 终局提交（裁决记录补执行注记）**

在 docs/plans/2026-10-14-reg-retirement-observables.md 尾部追加执行摘要段（各任务提交哈希+验证读数索引），单独提交。

---

## Self-Review 记录

- 规格覆盖：spec §①→Task 1；§②→Task 2/3；§③→Task 4 Step 1-2；§④→Task 4 Step 3；§⑤→各任务验证步+Task 5；§⑥→Global Constraints。无缺口。
- 占位符扫描：Task 3 差异表的 "照抄各自 reg" 指向磁盘上的权威源文件（reg 删除前搬运），要素值已内联；Task 2 的签名核对注记是核对指令非占位。
- 类型一致性：observables(子程序)/keys(函数)/term(逻辑函数)/method_name(字符函数) 与 prebuild 名集分类（sa/fk）一致；wires 模块名与文件名一致。
