# 纯势容器（Tier 1 auto-PES）实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 落地零 reg 自动接线——systems/ 下无 reg_\*.f90 的文件夹凭十魔数名（pes/keys/load/term/classify + prop/match/sync/arm/method_name）由 prebuild 织入 registry_gen.f90，新用户一个 .f90 文件接入自定义 PES 与非绝热方法包。

**Architecture:** 全部机械在 scripts/prebuild.sh（检测 + 八条有名中止 + 织入本就生成的 registry_gen.f90）与 src/Makefile（registry_gen.o 静态前提扩含 OBJ_SYS，闭首编竞态）两处；fixture 容器停架 systems/inactive/auto_sample/ 作判据程序与演练的载体；判据 = check_auto_pes.f90（物理数学面，dep-light）+ drill_auto_pes.sh（接线与生命周期面，真实二进制）；文档四件同步；全案例逐字节回归。

**Tech Stack:** Fortran 2008（gfortran 15.2 缺省 / ifx 第二编译器）、GNU bash + awk（prebuild 既有机械）、GNU make。

**Spec:** docs/plans/2026-10-12-auto-pes-container.md（含裁决链 D1-D5 与勘误批：pes 一维展平 q(:)，size=3×nat）——本计划从该文档论证，执行者两份都读。

## Global Constraints

- 构建唯一入口 `cd src && make`（gfortran 缺省）/ `make FC=ifx`；产物只落 `build/<FC>/`；绝不在 src/ 里编译。
- **无自动容器在场时 prebuild 输出必须与改动前逐字节相同**（六案例逐字节回归的根基）。
- T5：全库（src/ 与 systems/）零 `RANDOM_NUMBER`；T4：字面 `23.0605` 只许 force_interface.f90 命中。fixture 不引入新换算与新随机流。
- fixture 停架 `systems/inactive/auto_sample/`（两层深 = 构建盲区；激活 = 搬回 systems/ 顶层；一次一容器裁决不变）。
- Fortran 源风格：小写、文件头一行用途 + Design note ≤12 行、子程序头 1-2 行 `---` 包裹；新增源码不含 E 号叙事（那是 docs/plans 的事）。
- 提交信息：中文详叙（证据 + 验证腿）；每个 Task 一个提交。
- 案例回归锚：stdout 除 `timing:` 行外与入库 `results.txt` 一致（README A2 口径；若个别案例 results.txt 为摘录，以该案例 README 第 4 节口径为准）。
- 现有六案例：case_b_layer / case_chain_frequencies / case_harmonic_scatter / case_radau_conservation / case_spectrum_source / case_surface_window。
- 工作区注意：仓库有大量**未入库的 venus2026/ 旧树删除**（不属本批，勿动勿提交，零触碰检查须排除）；每次提交只 add 本批文件。
- 激活态的 fixture 是未跟踪路径（git 只认 inactive 停架位）——负面腿的还原规程 = `rm -rf` 激活副本 + `git checkout -- systems/inactive/auto_sample`，详见 Task 2 Step 5。

---

### Task 1: fixture 容器（systems/inactive/auto_sample/）

**Files:**
- Create: `systems/inactive/auto_sample/pot.f90`（pes/keys/term + 参数态 + setter）
- Create: `systems/inactive/auto_sample/io.f90`（load/classify——input/final_state 依赖隔离层）
- Create: `systems/inactive/auto_sample/tdhf_stub.f90`（方法组五名 + arm 见证文件）
- Create: `systems/inactive/auto_sample/system/frag_a.xyz`、`system/frag_b.xyz`（逐字复制 cases/case_harmonic_scatter/system/ 同名文件）

**Interfaces:**
- Produces（后续任务依赖）: 模块 `auto_sample_pot`（public: pes/keys/term/auto_sample_set_params）、模块 `auto_sample_io`（public: load/classify）、模块 `auto_sample_met`（public: prop/match/sync/arm/method_name，`method_name()=='auto_tdhf'`）；K_SPRING 缺省 5.0 eV/Å²；R_TERM 缺席 = MAX_STEPS-only。

**拆文件动机**：pot.f90 只依赖 state（term 要 state_t）——判据程序 Task 5 单编它不拖 input/final_state 链；io.f90 独居 input/final_state 依赖；load 经 pot 的公共 setter `auto_sample_set_params` 落位参数（E28 box 模式的压扁；setter 是 fixture 内部 API，cau_pars_\* 同例）。

- [ ] **Step 1: 写 pot.f90**

```fortran
!=====================================================================
! pot.f90 - the auto_sample fixture's PES group: origin-spring well
!   (Tier 1 pure-PES container example; pes/keys/term magic names)
! Design:
!   V = 0.5*k*sum(q_i^2), g = k*q - every atom springs to the origin,
!   geometry-agnostic (any nat), closed-form checkable. k stages from
!   K_SPRING (default 5.0 eV/A^2) via auto_sample_set_params (io.f90's
!   load hands it over); r_term/r_found stage from R_TERM (absent key =
!   the MAX_STEPS-only run, term stays false).
! Units: q[A], g[eV/A], v[eV]
!=====================================================================
module auto_sample_pot
   use state, only: state_t
   implicit none
   private
   public :: pes, keys, term
   public :: auto_sample_set_params   ! fixture wiring: io.f90 load's handover
   real(8) :: k_spring = 5.0d0        ! spring constant [eV/A^2] (K_SPRING)
   real(8) :: r_term = 0.0d0          ! termination radius [A] (R_TERM)
   logical :: r_term_found = .false.  ! R_TERM presence flag [-]
contains
   !------------------------------------------------------------------
   ! keys() - the container key vocabulary (pure data, pre-parse)
   !------------------------------------------------------------------
   function keys() result(kk)
      character(len=16) :: kk(2)
      kk = [character(len=16) :: 'K_SPRING', 'R_TERM']
   end function keys

   !------------------------------------------------------------------
   ! pes(q, mass, g, v) - the force entry (E27 delivery, flat dof)
   !------------------------------------------------------------------
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)      ! coordinates [A] (3*nat, flattened)
      real(8), intent(in)  :: mass(:)   ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)      ! dV/dq [eV/A]
      real(8), intent(out) :: v         ! V [eV]
      v = 0.5d0*k_spring*dot_product(q, q)
      g = k_spring*q
   end subroutine pes

   !------------------------------------------------------------------
   ! term(sta) - termination: any coordinate leaves the R_TERM ball
   !             (key absent -> always false = MAX_STEPS-only)
   !------------------------------------------------------------------
   logical function term(sta)
      type(state_t), intent(in) :: sta
      term = .false.
      if (r_term_found) term = any(abs(sta%q) > r_term)
   end function term

   !------------------------------------------------------------------
   ! auto_sample_set_params(k, r, found) - the load handover (io.f90)
   !------------------------------------------------------------------
   subroutine auto_sample_set_params(k, r, found)
      real(8), intent(in) :: k, r
      logical, intent(in) :: found
      k_spring = k
      r_term = r
      r_term_found = found
   end subroutine auto_sample_set_params
end module auto_sample_pot
```

- [ ] **Step 2: 写 io.f90**

```fortran
!=====================================================================
! io.f90 - the auto_sample fixture's input/output group: load/classify
!   magic names (the input/final_state dependency layer)
! Design:
!   load pulls K_SPRING/R_TERM buffered keys and hands them to pot's
!   module state via auto_sample_set_params (absent key keeps the
!   default; found carries the R_TERM optional-slot semantics).
!   classify is the wiring-exercise stub: constant 'bound' (fin unread)
!   - the drill asserts it end-to-end through the real binary.
! Units: none (judgments carry code references only)
!=====================================================================
module auto_sample_io
   use input,       only: pull
   use final_state, only: fin_t
   use auto_sample_pot, only: auto_sample_set_params
   implicit none
   private
   public :: load, classify
contains
   !------------------------------------------------------------------
   ! load() - pull the buffered container keys, stage pot's state
   !------------------------------------------------------------------
   subroutine load()
      real(8) :: k, r
      logical :: found
      k = 5.0d0                          ! field defaults ARE the container defaults
      r = 0.0d0
      found = .false.
      call pull('K_SPRING', k)
      call pull('R_TERM', r, found)
      call auto_sample_set_params(k, r, found)
   end subroutine load

   !------------------------------------------------------------------
   ! classify(fin) - wiring-exercise stub: the bound channel (fin unread)
   !------------------------------------------------------------------
   function classify(fin) result(name)
      type(fin_t), intent(in) :: fin     ! final-state quantities (unused)
      character(len=32) :: name
      name = 'bound'
   end function classify
end module auto_sample_io
```

- [ ] **Step 3: 写 tdhf_stub.f90**

```fortran
!=====================================================================
! tdhf_stub.f90 - the auto_sample fixture's method group: a no-op
!   package named auto_tdhf (prop/match/sync/arm/method_name)
! Design:
!   The no-op propagation and matching forms are the legal degenerate
!   members (elec_prop_i lists no-op; match without a hop writes
!   nothing - occ stays the seeded [1]). arm writes the witness file
!   auto_arm_witness.txt into the run directory - the drill's
!   selection-gate evidence (present iff the run armed the package).
! Units: none
!=====================================================================
module auto_sample_met
   use state, only: state_t, elec_t
   implicit none
   private
   public :: prop, match, sync, arm, method_name
   logical :: armed = .false.            ! witness latch [-]
contains
   !------------------------------------------------------------------
   ! method_name() - the registration name (lowercase, author data)
   !------------------------------------------------------------------
   function method_name() result(s)
      character(len=32) :: s
      s = 'auto_tdhf'
   end function method_name

   !------------------------------------------------------------------
   ! prop(sta, dt) - no-op propagation (legal degenerate form)
   !------------------------------------------------------------------
   subroutine prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine prop

   !------------------------------------------------------------------
   ! match(sta) - no-op matching (no hop, occ stays the seeded [1])
   !------------------------------------------------------------------
   subroutine match(sta)
      type(state_t), intent(inout) :: sta
   end subroutine match

   !------------------------------------------------------------------
   ! sync(s) - no-op container sync (no deliverable cache)
   !------------------------------------------------------------------
   subroutine sync(s)
      type(elec_t), intent(inout) :: s
   end subroutine sync

   !------------------------------------------------------------------
   ! arm() - the selected-package witness (writes the drill's file once)
   !------------------------------------------------------------------
   subroutine arm()
      if (armed) return
      armed = .true.
      open (20, file='auto_arm_witness.txt', status='replace', action='write')
      write (20, '(a)') 'auto_tdhf armed'
      close (20)
   end subroutine arm
end module auto_sample_met
```

- [ ] **Step 4: 复制构型文件并做独立语法编译验证**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026
mkdir -p systems/inactive/auto_sample/system
cp cases/case_harmonic_scatter/system/frag_a.xyz cases/case_harmonic_scatter/system/frag_b.xyz \
   systems/inactive/auto_sample/system/
mkdir -p build_tmp_auto && cd build_tmp_auto
gfortran -c -J. -o state.o ../src/globals/state.f90
gfortran -c -J. -I. -o pot.o ../systems/inactive/auto_sample/pot.f90
gfortran -c -J. -I. -o met.o ../systems/inactive/auto_sample/tdhf_stub.f90
cd .. && rm -rf build_tmp_auto
```
Expected: 三对象编译绿（io.f90 不在此验——它引 input/final_state，整树构建里编，Task 3 验证）。

- [ ] **Step 5: Commit**

```bash
git add systems/inactive/auto_sample/
git commit -m "auto_sample fixture 容器入库(Tier 1 纯势容器范例+判据载体): pot.f90(pes/keys/term 魔数名——原点弹簧井 V=0.5k·q², K_SPRING 缺省 5.0, R_TERM 缺席=MAX_STEPS; auto_sample_set_params 为 io 落位 setter, cau_pars_* 同例)/io.f90(load/classify——input/final_state 依赖隔离层, classify 为接线练习常量 'bound')/tdhf_stub.f90(prop/match/sync/method_name 合法无操作形 + arm 写 auto_arm_witness.txt 见证文件)/system/ 两片段构型逐字复制 case_harmonic_scatter。独立语法编译 state/pot/met 三绿(io 需整树链, Task 3)"
```

---

### Task 2: prebuild.sh 自动扫描与织入

**Files:**
- Modify: `scripts/prebuild.sh`（zmods 收集之后插入扫描段；registry_gen emission 块内三处插入）

**Interfaces:**
- Consumes: Task 1 的 fixture（三模块名）。
- Produces: 无自动容器时输出逐字节不变；fixture 激活时织入 auto_slots 手；八条 prebuild 有名中止（exit 1 + stderr 指名）。第 5 条中止（双容器 STOP 1）是运行期行为，归 Task 4。

- [ ] **Step 1: 基线快照（无自动容器时输出不变的前置锚）**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026/src
bash ../scripts/prebuild.sh . build/gfortran
cp build/gfortran/registry_gen.f90 /tmp/registry_gen.baseline
```

- [ ] **Step 2: 插入自动扫描段（zmods 的 `done` 之后、`{` emission 块之前）**

```bash
# ---------- auto-PES container scan (2026-10-12, docs/plans/2026-10-12-auto-pes-container.md) ----------
# A systems/* folder WITHOUT reg_*.f90 is a pure-PES container (Tier 1): the
# ten magic names (pes/keys/load/term/classify + prop/match/sync/arm/method_name)
# its modules expose get woven into registry_gen below. A folder WITH reg_*.f90
# is untouched (the explicit reg is authoritative - Tier 2).
auto_fatal() { # $1 = message (the design's named aborts, prebuild side)
   echo "prebuild.sh: auto-PES error: $1" >&2
   exit 1
}

# candidate folder: one level under systems/, holding .f90, no reg_*.f90
auto_dir=""
for d in "$SRC/../systems"/*/; do
   [ -d "$d" ] || continue
   ls "$d"*.f90 >/dev/null 2>&1 || continue
   ls "$d"reg_*.f90 >/dev/null 2>&1 && continue
   if [ -n "$auto_dir" ]; then
      auto_fatal "two auto containers are installed ('$auto_dir' and '$d') - one container per process"
   fi
   auto_dir="$d"
done

# magic-name scan: emits "kind name module" per definition (kind s/f).
# Comment tails are stripped and the keyword may sit behind a type/purity
# prefix - "logical function term(sta)" is the idiomatic declaration form
# (reference_box writes it; a line-start anchor would silently MISS it).
# "end subroutine pes" never matches: the name must carry the argument paren.
# A wrong-kind hit (subroutine keys / function pes) emits its actual kind -
# the caller's form check fatals on it.
auto_scan_file() { # $1 = source path
   awk '
      BEGIN { split("pes load prop match sync arm", sa); for (i in sa) sk[sa[i]] = 1
              split("keys term classify method_name", fa); for (i in fa) fk[fa[i]] = 1 }
      {
         line = $0
         sub(/!.*/, "", line)                       # strip the comment tail
         sub(/^[[:space:]]+/, "", line)             # normalize the lead
         if (line ~ /^module[[:space:]]+[a-zA-Z_]/) {
            split(line, mw, /[[:space:]]+/)
            if (mw[2] != "procedure") mod = mw[2]
         }
         if (match(line, /(subroutine|function)[[:space:]]+[a-zA-Z_]/)) {
            if (RSTART > 1 && substr(line, 1, RSTART - 1) ~ /end[[:space:]]*$/) next
            split(substr(line, RSTART), w, /[[:space:]]+/)
            kw = w[1]; nm = tolower(w[2]); sub(/\(.*/, "", nm)
            if (nm in sk || nm in fk) print (kw == "subroutine" ? "s" : "f"), nm, mod
         }
      }' "$1"
}

declare -A auto_mod                     # magic name -> owning module
auto_pair=no
if [ -n "$auto_dir" ]; then
   defs=$(for f in "$auto_dir"*.f90; do auto_scan_file "$f"; done)
   for nm in pes keys load term classify prop match sync arm method_name; do
      case "$nm" in pes|load|prop|match|sync|arm) want=s;; *) want=f;; esac
      got=$(printf '%s\n' "$defs" | awk -v n="$nm" '$2 == n {print}')
      cnt=$(printf '%s\n' "$got" | grep -c .) || true
      [ "$cnt" -eq 0 ] && continue
      [ "$cnt" -ge 2 ] && auto_fatal "magic name '$nm' is defined $cnt times in $auto_dir (exactly one allowed)"
      read -r kind _ mod _ <<<"$(printf '%s\n' "$got")"
      [ "$kind" != "$want" ] && auto_fatal "magic name '$nm' must be a $([ "$want" = s ] && echo subroutine || echo function), the other form was found in $auto_dir"
      [ -z "$mod" ] && auto_fatal "magic name '$nm' sits outside any module in $auto_dir (the weave needs a module owner)"
      auto_mod[$nm]=$mod
   done
   n_pm=0
   [ -n "${auto_mod[prop]}" ] && n_pm=$((n_pm+1))
   [ -n "${auto_mod[match]}" ] && n_pm=$((n_pm+1))
   if [ "$n_pm" -eq 1 ]; then
      auto_fatal "prop and match must be defined together in $auto_dir (a half package is not a method)"
   fi
   [ "$n_pm" -eq 2 ] && auto_pair=yes
   if [ "$auto_pair" = yes ]; then
      [ -z "${auto_mod[method_name]}" ] && auto_fatal "prop+match are defined but method_name is missing in $auto_dir (the registration name is author data)"
   else
      [ -n "${auto_mod[sync]}" ] && auto_fatal "sync is defined but prop+match are not (a lone sync member is not a package) in $auto_dir"
      [ -n "${auto_mod[arm]}" ] && auto_fatal "arm is defined but prop+match are not (a lone arm is not a package) in $auto_dir"
      [ -n "${auto_mod[method_name]}" ] && auto_fatal "method_name is defined but prop+match are not in $auto_dir"
   fi
   [ -z "${auto_mod[pes]}" ] && auto_fatal "$auto_dir is installed (folder present, no reg_*.f90) but wires nothing: expose exactly one 'pes' entry, or write a reg_*.f90 full container"
fi
```

（`read -r kind _ mod _`：got 行形如 `s pes auto_sample_pot`——第二字段是名、第三字段是归属模块。）

- [ ] **Step 3: 织入 emission 块（`{ ... } > "$BLD/registry_gen.f90.tmp"` 内三处插入）**

插入点 A——`use container_sched` 两行之后、`for m in $mods` 循环之前：

```bash
   if [ -n "$auto_dir" ]; then
      [ -n "${auto_mod[term]}" ] && echo "   use force_interface, only: container_bind_pes, container_bind_term"
      [ -z "${auto_mod[term]}" ] && echo "   use force_interface, only: container_bind_pes"
      [ -n "${auto_mod[classify]}" ] && echo "   use final_state,     only: container_bind_classify"
      [ -n "${auto_mod[keys]}" ] && echo "   use input,           only: input_declare_keys"
      if [ "$auto_pair" = yes ]; then
         echo "   use elec_interface,  only: method_reg"
         [ -n "${auto_mod[arm]}" ] && echo "   use config,          only: electronic"
      fi
      for m in $(printf '%s\n' "${auto_mod[@]}" | sort -u); do
         names=""; first=1
         for nm in pes keys load term classify prop match sync arm method_name; do
            if [ "${auto_mod[$nm]}" = "$m" ]; then
               if [ "$first" -eq 1 ]; then names="$nm"; first=0; else names="$names, $nm"; fi
            fi
         done
         echo "   use $m, only: $names"
      done
   fi
```

插入点 B——`echo "contains"` 之后、`echo "   subroutine registry_gen_keys()"` 之前：

```bash
   if [ -n "$auto_dir" ]; then
      args=""
      [ -n "${auto_mod[keys]}" ]     && args="$args keys"
      [ -n "${auto_mod[load]}" ]     && args="$args load"
      args="$args force"
      [ "$auto_pair" = yes ]         && args="$args method"
      [ -n "${auto_mod[term]}" ]     && args="$args term"
      [ -n "${auto_mod[arm]}" ]      && args="$args arm"
      [ -n "${auto_mod[classify]}" ] && args="$args classify"
      echo "   subroutine auto_slots()"
      ia=0
      for a in $args; do
         ia=$((ia+1))
         if [ "$ia" -eq 1 ]; then
            printf '      call container_sched_bind(%s' "$a"
         else
            printf ', &\n                                %s' "$a"
         fi
      done
      echo ")"
      echo "   end subroutine auto_slots"
      [ -n "${auto_mod[keys]}" ] && { echo "   subroutine auto_keys()"
         echo "      call input_declare_keys(keys())"
         echo "   end subroutine auto_keys"; }
      [ -n "${auto_mod[load]}" ] && { echo "   subroutine auto_load()"
         echo "      call load()"
         echo "   end subroutine auto_load"; }
      echo "   subroutine auto_force()"
      echo "      call container_bind_pes(pes)"
      echo "   end subroutine auto_force"
      if [ "$auto_pair" = yes ]; then
         echo "   subroutine auto_method()"
         if [ -n "${auto_mod[sync]}" ]; then echo "      call method_reg(method_name(), prop, match, sync)"
         else echo "      call method_reg(method_name(), prop, match)"; fi
         echo "   end subroutine auto_method"
      fi
      [ -n "${auto_mod[term]}" ] && { echo "   subroutine auto_term()"
         echo "      call container_bind_term(term)"
         echo "   end subroutine auto_term"; }
      [ -n "${auto_mod[arm]}" ] && { echo "   subroutine auto_arm()"
         echo "      if (trim(electronic%method) == method_name()) call arm()"
         echo "   end subroutine auto_arm"; }
      [ -n "${auto_mod[classify]}" ] && { echo "   subroutine auto_classify()"
         echo "      call container_bind_classify(classify)"
         echo "   end subroutine auto_classify"; }
   fi
```

插入点 C——`registry_gen_keys()` 体内、zmods 调用循环之前：

```bash
   [ -n "$auto_dir" ] && echo "      call auto_slots()"
```

emission 头注（文件头 `! registry_gen.f90 - BUILD-GENERATED...` 注记块末）追加两行：

```bash
   echo "! An auto container's weave (auto_slots) enters here when a reg-less"
   echo "! systems/ folder is present - Tier 1, docs/plans/2026-10-12-auto-pes-container.md."
```

- [ ] **Step 4: 正面验证（fixture 激活 → 织入内容断言 → 基线回归）**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026
mv systems/inactive/auto_sample systems/auto_sample
(cd src && bash ../scripts/prebuild.sh . build/gfortran)
for pat in \
   "use auto_sample_pot, only: pes, keys, term" \
   "use auto_sample_io, only: load, classify" \
   "use auto_sample_met, only: prop, match, sync, arm, method_name" \
   "use elec_interface,  only: method_reg" \
   "use config,          only: electronic" \
   "call method_reg(method_name(), prop, match, sync)" \
   "if (trim(electronic%method) == method_name()) call arm()" \
   "call auto_slots()"; do
   grep -qF "$pat" src/build/gfortran/registry_gen.f90 && echo "WOVEN-OK: $pat" || echo "WOVEN-FAIL: $pat"
done
mv systems/auto_sample systems/inactive/auto_sample
(cd src && bash ../scripts/prebuild.sh . build/gfortran)
cmp src/build/gfortran/registry_gen.f90 /tmp/registry_gen.baseline && echo BASELINE-IDENTICAL
```
Expected: 八行 `WOVEN-OK` + `BASELINE-IDENTICAL`。

- [ ] **Step 5: 负面验证（prebuild 中止逐条）**

每条腿的固定规程（激活副本是未跟踪路径，还原 = 整删 + checkout 停架位）：

```bash
# 腿模板（cwd = 仓库根）：
mv systems/inactive/auto_sample systems/auto_sample
# <按表构造违规>
(cd src && bash ../scripts/prebuild.sh . build/gfortrun 2>/dev/null; true)   # 占位行，见下——实际跑 build/gfortran
```
实际模板（每腿执行）：

```bash
mv systems/inactive/auto_sample systems/auto_sample
# <构造违规（见下表）>
( cd src && bash ../scripts/prebuild.sh . build/gfortran ) > /tmp/neg.out 2>&1
echo "exit=$?"; grep -o "auto-PES error: .*" /tmp/neg.out | head -1
rm -rf systems/auto_sample
git checkout -- systems/inactive/auto_sample
```

| # | 构造（作用于激活副本 systems/auto_sample/） | 期望 stderr 含 |
|---|---|---|
| 1 | `mv systems/auto_sample/pot.f90 /tmp/pot_keep.f90` | `wires nothing: expose exactly one 'pes'` |
| 2 | `cat >> systems/auto_sample/pot.f90 <<'EOT'`<br/>`module dup_mod`<br/>`contains`<br/>`subroutine pes(q, mass, g, v)`<br/>`real(8), intent(in) :: q(:), mass(:)`<br/>`real(8), intent(out) :: g(:), v`<br/>`end subroutine`<br/>`end module`<br/>`EOT` | `'pes' is defined 2 times` |
| 3 | pot.f90 内 `function keys` 改 `subroutine keys`（sed 或手编） | `must be a function, the other form` |
| 4 | `cp -r systems/auto_sample systems/auto_two` | `two auto containers are installed` |
| 6 | tdhf_stub.f90 内 `subroutine match` 改名 `subroutine match2`（含 public 行） | `prop and match must be defined together` |
| 7 | tdhf_stub.f90 内 `function method_name` 改名 `method_label`（含 public 行） | `method_name is missing` |
| 8 | tdhf_stub.f90 内注释掉 prop 与 match 两个子程序（保留 sync/arm/method_name） | `sync is defined but prop+match are not` |
| 8b | 同上但只保留 method_name（注释 prop/match/sync/arm） | `method_name is defined but prop+match are not` |

Expected：每腿 exit=1 且 grep 命中期望文案；`git checkout` 后停架位无 diff（`git status --porcelain systems/` 仅显预存无关项）。

- [ ] **Step 6: Commit**

```bash
git add scripts/prebuild.sh
git commit -m "prebuild 自动扫描与织入(Tier 1): 无 reg 的 systems/* 文件夹十魔数名检测(锚定行首 subroutine|function 定义, 大小写不敏感, 名-形核对, 归属=定义前最近 module 语句)+八条有名中止(exit 1 指名: pes 0/名≥2/形不符/双自动文件夹/prop XOR match/对无 method_name/method_name 无对/sync|arm 无对)+registry_gen 织入(interface use 行按在场/逐模块 use 名单按固定名序逗号列/auto_slots 手与七个 auto_* 体/registry_gen_keys 内 call auto_slots() 先于 zmods); 验证: 织入八断言全中, 无自动容器输出与基线 cmp 逐字节一致, 负面八腿逐条指名中止"
```

---

### Task 3: Makefile 静态前提 + 首编即绿

**Files:**
- Modify: `src/Makefile`（registry_gen.o 规则及其头注，现约 94-105 行）

**Interfaces:**
- Produces: `$(BUILD)/registry_gen.o: ... $(OBJ_REG) $(OBJ_SYS)`——拖入即单次 make 首编即绿（织入代码的 Fortran 编译期验证同腿完成）。

- [ ] **Step 1: 改规则与头注**

把规则上方的解释注（"the generated unit gets an EXPLICIT rule..."段）末尾补两句、规则行加 `$(OBJ_SYS)`：

```make
# ...（原注释前三句保留）... - a static rule closes that window. The set covers
# reg objects AND all systems objects: the auto-PES weave (2026-10-12) makes
# registry_gen `use` NON-reg systems modules (a Tier 1 container's pot/io/met),
# and those edges must be ordered on the FIRST make after the drag, not the
# second. Ruling: docs/plans/2026-10-12-auto-pes-container.md §2.
SRC_REGM := $(wildcard methods/reg_*.f90)
SRC_REGS := $(wildcard ../systems/*/reg_*.f90)
OBJ_REG  := $(patsubst %.f90,$(BUILD)/%.o,$(SRC_REGM)) \
            $(patsubst ../systems/%.f90,$(BUILD)/sys/%.o,$(SRC_REGS))
$(BUILD)/registry_gen.o: $(BUILD)/registry_gen.f90 $(OBJ_REG) $(OBJ_SYS) | $(BUILD)
	$(FC) $(FFLAGS) $< -o $@
```

- [ ] **Step 2: 首编即绿验证（干净树 + 拖入 + 单次 make）**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026
mv systems/inactive/auto_sample systems/auto_sample
(cd src && make clean && make) > /tmp/first_make.log 2>&1
echo "EXIT=$?"; tail -3 /tmp/first_make.log
ls src/build/gfortran/venus.e
```
Expected: `EXIT=0`、venus.e 存在、log 无 "Can't open module"/"File not found"（单次 make 全绿 = 织入代码通过编译器的 use/签名核对）。

- [ ] **Step 3: 还原并提交**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026
mv systems/auto_sample systems/inactive/auto_sample
(cd src && make) > /tmp/restore_make.log 2>&1; echo "EXIT=$?"
git add src/Makefile
git commit -m "Makefile: registry_gen.o 静态前提扩含 OBJ_SYS(Tier 1 首编即绿): 自动织入使生成单元 use 非 reg 的 systems 模块(拖入容器的 pot/io/met), 该边在 depend.mk 单轮迟滞窗口内需静态排序; 验证——拖入后 make clean+单次 make 全绿(织入代码过编译器 use/签名核对), 还原停架后重建绿"
```

---

### Task 4: 演练脚本 scripts/drill_auto_pes.sh

**Files:**
- Create: `scripts/drill_auto_pes.sh`（可执行）

**Interfaces:**
- Consumes: Task 1-3 全部。
- Produces: 一键演练——正面（首编绿/织入断言/绝热缺省跑/方法跑 arm 见证/R_TERM 两分支/双跑确定性）+ 运行期负面（双容器 STOP 1）+ 零触碰收尾；任何退出形态都还原树（trap）。

- [ ] **Step 1: 写脚本**

```bash
#!/usr/bin/env bash
# drill_auto_pes.sh - the Tier 1 auto-PES container drill (criterion family):
#   positive legs (first-make green, woven content, adiabatic default run with
#   the arm witness ABSENT, ELEC_METHOD = auto_tdhf run with the witness
#   PRESENT, R_TERM both branches, double-run determinism) + the runtime
#   coexistence abort (reference + auto = STOP 1) + the zero-touch close.
# Leaves the tree as found (trap restores); exits nonzero on any failed leg.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/src"; BLD="$SRC/build/gfortran"
FIX="$ROOT/systems/auto_sample"; SHELF="$ROOT/systems/inactive/auto_sample"
REF="$ROOT/systems/reference"; REFSHELF="$ROOT/systems/inactive/reference"
DRILL="$BLD/drill_case"
pass=0; fail=0
ok()  { echo "PASS: $1"; pass=$((pass+1)); }
bad() { echo "FAIL: $1"; fail=$((fail+1)); }

restore_tree() {
   [ -d "$FIX" ] && mv "$FIX" "$SHELF"
   [ -d "$REFSHELF" ] && mv "$REFSHELF" "$REF"
   rm -rf "$DRILL"
}
trap restore_tree EXIT

cd "$SRC"

# ---- leg 1: drag in, FIRST make green ----
mv "$SHELF" "$FIX"
if make >/tmp/drill_make.log 2>&1; then ok "first make green (drag-in, no second make)"
else bad "first make"; tail -20 /tmp/drill_make.log; exit 1; fi

# ---- leg 2: woven content ----
while IFS= read -r pat; do
   grep -qF "$pat" "$BLD/registry_gen.f90" && ok "woven: $pat" || bad "woven: $pat"
done <<'PATS'
use auto_sample_pot, only: pes, keys, term
use auto_sample_io, only: load, classify
use auto_sample_met, only: prop, match, sync, arm, method_name
call method_reg(method_name(), prop, match, sync)
if (trim(electronic%method) == method_name()) call arm()
call auto_slots()
PATS

# ---- leg 3: runtime coexistence abort (reference + auto = STOP 1 at assembly;
#      registry_gen_keys runs before the input read, so no input file needed) ----
"$BLD/venus.e" >/dev/null 2>/tmp/drill_coexist.err
if grep -q "slot is already filled" /tmp/drill_coexist.err; then ok "coexistence STOP 1"
else bad "coexistence STOP 1"; cat /tmp/drill_coexist.err; fi

# ---- leg 4: shelve reference, rebuild, scratch-case runs ----
mv "$REF" "$REFSHELF"
make >/dev/null 2>&1 || { bad "shelved-reference rebuild"; exit 1; }
mkdir -p "$DRILL/system"
cp "$FIX/system/frag_a.xyz" "$FIX/system/frag_b.xyz" "$DRILL/system/"
cat > "$DRILL/input_qct.txt" <<'EOF'
# auto_sample drill case - Tier 1 auto wiring
SYSTEM_DIR = system
TITLE = auto sample drill - origin-spring well
E_REL = 50.0
B_MAX = 1.0
R_SEP = 6.0
DT = 0.05
DIST_SCHEME = THERMALIZE,ROTATION
PROJECTILE = frag_a
TARGET = frag_b
T_EQ = 100.0
N_EQ = 5
DT_EQ = 0.05
J_ROT = 3
N_B = 1
N_TRAJ = 2
MAX_STEPS = 50
ISEED = 20261012
GWRITE_LEVEL = 2
# ---- auto_sample container keys ----
K_SPRING = 0.05
EOF
cd "$DRILL"
rm -f auto_arm_witness.txt
if "$BLD/venus.e" >stdout_a.txt 2>stderr_a.txt; then ok "adiabatic run (default ELEC_METHOD)"
else bad "adiabatic run"; cat stderr_a.txt; fi
[ ! -f auto_arm_witness.txt ] && ok "arm NOT called under adiabatic default" || bad "arm witness present under adiabatic"
grep -q "outcome counts: bound=2" stdout_a.txt && ok "classify stub ran (bound=2)" || bad "classify outcome line"
# determinism: run again, trajectories byte-identical
cp trajectories.txt /tmp/drill_traj_first.txt
"$BLD/venus.e" >/dev/null 2>&1
cmp -s /tmp/drill_traj_first.txt trajectories.txt && ok "double-run determinism" || bad "determinism"
# R_TERM branch: present key -> immediate termination, still classified
printf 'R_TERM = 2.0\n' >> input_qct.txt
if "$BLD/venus.e" >stdout_b.txt 2>stderr_b.txt; then ok "R_TERM run"
else bad "R_TERM run"; cat stderr_b.txt; fi
grep -q "outcome counts: bound=2" stdout_b.txt && ok "R_TERM branch classified" || bad "R_TERM outcome"
# method branch: select the package -> arm witness appears
sed -i '/^K_SPRING = 0.05/a ELEC_METHOD = auto_tdhf' input_qct.txt
rm -f auto_arm_witness.txt
if "$BLD/venus.e" >stdout_c.txt 2>stderr_c.txt; then ok "auto_tdhf run"
else bad "auto_tdhf run"; cat stderr_c.txt; fi
[ -f auto_arm_witness.txt ] && ok "arm witness written under ELEC_METHOD = auto_tdhf" || bad "arm witness missing"

# ---- leg 5: zero-touch close ----
cd "$SRC"; restore_tree; trap - EXIT
make >/dev/null 2>&1 || { bad "restored rebuild"; exit 1; }
dirty=$(git -C "$ROOT" status --porcelain | grep -v 'venus2026/' | grep -v '^??')
if [ -z "$dirty" ]; then ok "zero-touch (tracked tree clean, pre-existing deletions excluded)"
else bad "zero-touch"; echo "$dirty" | head; fi
echo "drill summary: pass=$pass fail=$fail"
[ "$fail" -eq 0 ]
```

- [ ] **Step 2: 跑演练至全绿**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026
chmod +x scripts/drill_auto_pes.sh && bash scripts/drill_auto_pes.sh 2>&1 | tee /tmp/drill_run.log
```
Expected: 末行 `drill summary: pass=17 fail=0`，退出码 0。失败腿的修复属本任务交付（fixture/织入/输入键三方都可能要动），修完重跑至全绿。

- [ ] **Step 3: Commit**

```bash
git add scripts/drill_auto_pes.sh
git commit -m "drill_auto_pes.sh 演练脚本入库(Tier 1 判据族): 十七腿——拖入首编即绿/织入内容六断言/双容器并存运行期 STOP 1/暂存 reference 后演练案例(绝热缺省跑 arm 不调+classify 常量出线, 双跑 trajectories 逐字节确定, R_TERM=2.0 早终止分支, ELEC_METHOD=auto_tdhf arm 见证文件在场)/零触碰收尾(tracked 净, 预存旧树删除显式排除); trap 保证任何退出形态还原树"
```

---

### Task 5: 判据程序 scripts/check_auto_pes.f90

**Files:**
- Create: `scripts/check_auto_pes.f90`

**Interfaces:**
- Consumes: `auto_sample_pot`（pes/keys）、`auto_sample_met`（method_name）、`state` 模块——零 input/final_state 链（load/classify/arm 的接线面归 Task 4 演练）。
- Produces: 四判据 + 方法名判据全 PASS 的独立程序。

- [ ] **Step 1: 写程序**

```fortran
!=====================================================================
! check_auto_pes.f90 - the auto_sample fixture's criterion program
!   (criteria defined FIRST in this header; this program executes them)
! Criteria:
!   (a) pes vs closed form over a battery (nat=1..4, k = the staged
!       default 5.0 eV/A^2): v = 0.5*k*sum(q^2), g = k*q,
!       max |d| <= 1e-12 (exact arithmetic)
!   (b) pes vs central-difference gradient (h = 1e-5 A, nat=3 fixed
!       battery): max |g - g_fd| <= 1e-8 eV/A (truncation ~ k*h^2/6)
!   (c) keys() is exactly ['K_SPRING', 'R_TERM'] (count + strings)
!   (d) full assignment: g/v prefilled with huge() come back finite
!   (e) method_name() == 'auto_tdhf' and is lowercase
! Compile (cwd = scripts/):
!   mkdir -p build && gfortran -c -Jbuild -o build/state.o ../src/globals/state.f90
!   gfortran -c -Jbuild -Ibuild -o build/auto_sample_pot.o ../systems/inactive/auto_sample/pot.f90
!   gfortran -c -Jbuild -Ibuild -o build/auto_sample_met.o ../systems/inactive/auto_sample/tdhf_stub.f90
!   gfortran -c -Jbuild -Ibuild -o build/check_auto_pes.o check_auto_pes.f90
!   gfortran -o build/check_auto_pes.exe build/state.o build/auto_sample_pot.o build/auto_sample_met.o build/check_auto_pes.o
!=====================================================================
program check_auto_pes
   use auto_sample_pot, only: pes, keys
   use auto_sample_met,  only: method_name
   implicit none
   integer :: nfail
   nfail = 0
   call crit_a; call crit_b; call crit_c; call crit_d; call crit_e
   if (nfail == 0) then
      print '(a)', 'check_auto_pes: ALL PASS'
   else
      print '(a,i0,a)', 'check_auto_pes: ', nfail, ' FAILURES'
      stop 1
   end if
contains
   subroutine crit_a
      integer, parameter :: nats(4) = [1, 2, 3, 4]
      real(8), parameter :: k_ref = 5.0d0   ! the staged default (criterion constant)
      real(8), allocatable :: q(:), m(:), g(:)
      real(8) :: v, dmax
      integer :: it, i
      dmax = 0.0d0
      do it = 1, size(nats)
         allocate (q(3*nats(it)), m(nats(it)), g(3*nats(it)))
         q = [(0.37d0*real(i, 8) + 0.11d0*real(it, 8), i = 1, 3*nats(it))]
         m = real(nats(it), 8)
         call pes(q, m, g, v)
         dmax = max(dmax, abs(v - 0.5d0*k_ref*dot_product(q, q)), &
                    maxval(abs(g - k_ref*q)))
         deallocate (q, m, g)
      end do
      if (dmax <= 1.0d-12) then
         print '(a,es12.3)', 'PASS: (a) closed form, max|d| = ', dmax
      else
         nfail = nfail + 1; print '(a,es12.3)', 'FAIL: (a) closed form, max|d| = ', dmax
      end if
   end subroutine crit_a
   subroutine crit_b
      real(8), parameter :: h = 1.0d-5
      real(8) :: q(9), m(3), g(9), gp(9), gd(9), gfd(9), v, v_hi, v_lo, dmax
      integer :: i
      q = [0.51d0, -1.32d0, 0.74d0, 2.05d0, -0.18d0, 1.44d0, -0.92d0, 0.33d0, 1.11d0]
      m = 1.0d0
      call pes(q, m, g, v)
      do i = 1, 9
         q(i) = q(i) + h;  call pes(q, m, gp, v);  v_hi = v
         q(i) = q(i) - 2*h; call pes(q, m, gd, v);  v_lo = v
         q(i) = q(i) + h
         gfd(i) = -(v_hi - v_lo)/(2*h)          ! g = dV/dq
      end do
      dmax = maxval(abs(g - gfd))
      if (dmax <= 1.0d-8) then
         print '(a,es12.3)', 'PASS: (b) central difference, max|d| = ', dmax
      else
         nfail = nfail + 1; print '(a,es12.3)', 'FAIL: (b) central difference, max|d| = ', dmax
      end if
   end subroutine crit_b
   subroutine crit_c
      character(len=16) :: kk(2)
      kk = keys()
      if (trim(kk(1)) == 'K_SPRING' .and. trim(kk(2)) == 'R_TERM') then
         print '(a)', 'PASS: (c) keys() vocabulary exact'
      else
         nfail = nfail + 1; print '(a)', 'FAIL: (c) keys() vocabulary'
      end if
   end subroutine crit_c
   subroutine crit_d
      real(8) :: q(3), m(1), g(3), v
      q = [0.1d0, 0.2d0, 0.3d0]; m = 1.0d0
      g = huge(1.0d0); v = huge(1.0d0)         ! garbage prefill
      call pes(q, m, g, v)
      if (all(abs(g) < huge(1.0d0)) .and. abs(v) < huge(1.0d0)) then
         print '(a)', 'PASS: (d) full assignment (garbage overwritten)'
      else
         nfail = nfail + 1; print '(a)', 'FAIL: (d) full assignment'
      end if
   end subroutine crit_d
   subroutine crit_e
      character(len=32) :: s, lows
      s = method_name()
      lows = s
      call to_lower(lows)
      if (trim(s) == 'auto_tdhf' .and. s == lows) then
         print '(a)', 'PASS: (e) method_name() == ''auto_tdhf'' (lowercase)'
      else
         nfail = nfail + 1; print '(a)', 'FAIL: (e) method_name()'
      end if
   end subroutine crit_e
   subroutine to_lower(str)
      character(len=*), intent(inout) :: str
      integer :: i, ic
      do i = 1, len_trim(str)
         ic = iachar(str(i:i))
         if (ic >= iachar('A') .and. ic <= iachar('Z')) str(i:i) = achar(ic + 32)
      end do
   end subroutine to_lower
end program check_auto_pes
```

- [ ] **Step 2: 编译并跑**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026/scripts
mkdir -p build
gfortran -c -Jbuild -o build/state.o ../src/globals/state.f90
gfortran -c -Jbuild -Ibuild -o build/auto_sample_pot.o ../systems/inactive/auto_sample/pot.f90
gfortran -c -Jbuild -Ibuild -o build/auto_sample_met.o ../systems/inactive/auto_sample/tdhf_stub.f90
gfortran -c -Jbuild -Ibuild -o build/check_auto_pes.o check_auto_pes.f90
gfortran -o build/check_auto_pes.exe build/state.o build/auto_sample_pot.o build/auto_sample_met.o build/check_auto_pes.o
./build/check_auto_pes.exe | tee check_auto_pes.out
```
Expected: 六行 PASS + `check_auto_pes: ALL PASS`。

- [ ] **Step 3: Commit**

```bash
git add scripts/check_auto_pes.f90 scripts/check_auto_pes.out
git commit -m "check_auto_pes.f90 判据程序入库(Tier 1 物理数学面): 五判据——(a) 弹簧井闭式 nat=1..4 电池 max|d|≤1e-12 (b) 中心差分梯度 h=1e-5 ≤1e-8 eV/A (c) keys() 词表恰为 K_SPRING/R_TERM (d) 全赋值(huge 预填被覆写) (e) method_name()=='auto_tdhf' 且小写; 单编三对象(state/pot/met)零 input 链——load/classify/arm 接线面归 drill_auto_pes.sh; 运行 ALL PASS 读数落盘并入库 check_auto_pes.out"
```

---

### Task 6: 文档四件（EXTENDING 领篇重构 + 三处注册段增补）

**Files:**
- Modify: `src/EXTENDING.md`（§1 拆 1a/1b：Tier 1 十名词表领篇，现三件套降为 Tier 2）
- Modify: `AGENTS.md`（Registration 段增补）
- Modify: `src/CLAUDE.md`（构建/判据要点增补）
- Modify: `README.md`（状态段增补）

**Interfaces:**
- Consumes: Task 1-5 已入库的机制、fixture 与判据名。

- [ ] **Step 1: EXTENDING.md §1 重构**

文件头三行导语（"每类扩展 = 新文件 + 薄 reg 文件..."）替换为：

```markdown
# EXTENDING.md — 新程序扩展点(E27 后契约)

接一个自定义势能面（含可选非绝热方法包）的**最短路径是 Tier 1**：一个文件夹、
若干 .f90、零 reg 文件。需要更复杂的接线（多方法包、自定义记录列、复杂参数体）
时升级 Tier 2（§1b 的 reg 三件套——reg 在场则自动接线静默退出，pot 侧魔数名
函数全部可被 reg 直接复用，升级零迁移）。
```

原 `## 1. 物理系统容器（systems/<system>/）——体系专属代码的组合点` 整段之前插入新小节 `## 1a. Tier 1 纯势容器——十魔数名词表（封闭集合）`，内容：

```markdown
## 1a. Tier 1 纯势容器——十魔数名词表（封闭集合）

`systems/<sys>/` 文件夹**不含 reg_\*.f90** 即 Tier 1 容器。文件夹内 .f90 暴露
下列魔数名（名字即约定，签名即既有接口，编译期核对；裁决记录
docs/plans/2026-10-12-auto-pes-container.md）：

**势能面组**

| 魔数名 | 形 | 必备性 | 说明 |
|---|---|---|---|
| `pes(q, mass, g, v)` | subroutine | 恰一个 | E27 力交付：q[A]/mass[amu] 进（q/g 一维展平 size=3×nat），g[eV/Å]/v[eV] 全赋值出 |
| `keys()` | function | 可选 | `character(len=16)` 数组容器键词表（parse 前声明） |
| `load()` | subroutine | 可选 | 无参自持 `pull`；**亦即 seam 期初始化钩子**——读参数文件（NN 权重模式）/分配在此，先于一切力消费者与 arm |
| `term(sta)` | function | 可选 | 终止判断（键缺席门控用 found 旗标自实现；未绑 = 跑满 MAX_STEPS） |
| `classify(fin)` | function | 可选 | 通道分类（返回通道名） |

**方法包组**（自定义非绝热方法）

| 魔数名 | 形 | 必备性 | 说明 |
|---|---|---|---|
| `prop(sta, dt)` | subroutine | 与 match 成对 | 电子传播成员（`elec_prop_i` 同型） |
| `match(sta)` | subroutine | 与 prop 成对 | 统计匹配成员（跳变、occ/n_hop 写回） |
| `sync(s)` | subroutine | 可选随对 | a↔rho 表象换算尾巴 |
| `arm()` | subroutine | 可选随对 | 选中后初始化（门控归生成代码；参数自持走 keys/load） |
| `method_name()` | function | 随对必备 | 返回 `character(len=32)` **小写**注册名；`ELEC_METHOD = <名>` 选用 |

规则一句话：势能面组在场即绑；方法组注册无条件、武装有门；reg 在场则十名全忽略
（reg 权威）；名字须 public；prebuild 对 pes 0/≥2、任一名 ≥2、形不符、prop 与
match 不成对、method_name 缺对、sync/arm 孤立、双自动文件夹各有一条指名中止。

现成范本：`systems/inactive/auto_sample/`（pot.f90 / io.f90 / tdhf_stub.f90；
判据 `scripts/check_auto_pes.f90` + 演练 `scripts/drill_auto_pes.sh`）。
激活 = 文件夹搬回 systems/ 顶层（一次一容器；与 reference 并存会在装配期
STOP 1——测试时把 reference 暂存 systems/inactive/）。
```

原 §1 标题改为 `## 1b. Tier 2 全容器——reg 三件套（体系专属代码的组合点）`，正文原样保留。

- [ ] **Step 2: 三处注册段增补**

AGENTS.md 的 "Registration is a BUILD ARTIFACT" 段末追加：

```markdown
- Tier 1 pure-PES containers (2026-10-12): a systems/* folder WITHOUT reg_*.f90 is
  auto-wired by prebuild from ten magic names (pes/keys/load/term/classify +
  prop/match/sync/arm/method_name) - no reg file at all; a reg_*.f90 present disables
  the weave (reg authoritative). Ruling: docs/plans/2026-10-12-auto-pes-container.md;
  fixture systems/inactive/auto_sample; criterion check_auto_pes.f90 + drill_auto_pes.sh.
```

src/CLAUDE.md 的判据/构建要点列表追加一行：

```markdown
- Tier 1 纯势容器：无 reg 文件夹由 prebuild 按十魔数名自动接线（判据
  scripts/check_auto_pes.f90 + scripts/drill_auto_pes.sh；裁决
  docs/plans/2026-10-12-auto-pes-container.md）。
```

README.md 状态段追加一句：

```markdown
自定义 PES 接入新增 Tier 1 快通道——单文件魔数名容器零 reg 接入
（docs/plans/2026-10-12-auto-pes-container.md）。
```

- [ ] **Step 3: 一致性检查与提交**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026
grep -c "Tier 1" src/EXTENDING.md AGENTS.md src/CLAUDE.md README.md   # 期望四行均 ≥1
grep -n "五名" src/EXTENDING.md                                       # 期望无输出
git add src/EXTENDING.md AGENTS.md src/CLAUDE.md README.md
git commit -m "文档四件 Tier 1 增补批: EXTENDING.md 头部导语与 §1 拆 1a/1b——1a Tier 1 十名词表领篇(两组表格+规则一句话+auto_sample 范本与判据指引+激活/暂存说明), 1b reg 三件套原样降级; AGENTS.md 注册段/CLAUDE.md 判据要点/README 状态段各增 Tier 1 一句(裁决出处+fixture+判据名)"
```

---

### Task 7: 回归与收尾

**Files:**
- Modify: `docs/plans/2026-10-12-auto-pes-container.md`（§6 状态 → 已实现）

**Interfaces:**
- Consumes: 全部前序任务（fixture 停架态：reference 独占 systems/）。

- [ ] **Step 1: 双编译器干净重建**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026/src
make clean && make > /tmp/reg_gfortran.log 2>&1; echo "GFORTRAN_EXIT=$?"
make FC=ifx > /tmp/reg_ifx.log 2>&1; echo "IFX_EXIT=$?"
```
Expected: 双 EXIT=0。

- [ ] **Step 2: 六案例逐字节回归**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026
for c in case_b_layer case_chain_frequencies case_harmonic_scatter case_radau_conservation case_spectrum_source case_surface_window; do
  (cd "cases/$c" && ../../src/build/gfortran/venus.e > /tmp/reg_stdout.txt 2>/tmp/reg_stderr.txt
   grep -v '^timing:' /tmp/reg_stdout.txt | grep -v '^#' > /tmp/reg_now.txt
   grep -v '^timing:' results.txt | grep -v '^#' > /tmp/reg_arch.txt
   if cmp -s /tmp/reg_now.txt /tmp/reg_arch.txt; then echo "PASS $c"; else echo "FAIL $c (口径差异见该案例 README §4)"; fi)
done
```
Expected: 六行 PASS。若个别 FAIL：先查该案例 README 第 4 节的入库口径（results.txt 可能是摘录）——口径一致的差异才是真回归；如实记录每一案例的比对结论。

- [ ] **Step 3: T4/T5 复验**

```bash
cd C:/Users/13614/Desktop/MAIN/01_WORK/dev/venus2026/VENUS2026
grep -rn "23.0605" src/ systems/ | grep -v force_interface.f90    # 期望无输出
grep -rn "RANDOM_NUMBER" src/ systems/                            # 期望无输出
```

- [ ] **Step 4: 判据族终跑与状态收尾**

```bash
cd scripts && ./build/check_auto_pes.exe | tail -1        # ALL PASS
cd .. && bash scripts/drill_auto_pes.sh 2>&1 | tail -1    # drill summary: pass=17 fail=0
```
编辑 docs/plans/2026-10-12-auto-pes-container.md §6 追加：

```markdown
- 2026-10-XX：实现落地（实施计划 docs/plans/2026-10-12-auto-pes-implementation.md）——
  prebuild 织入 + Makefile 静态前提 + auto_sample fixture + check/drill 判据 + 文档四件；
  六案例回归绿、双编译器绿、T4/T5 复验零命中、drill 全绿。状态：已实现。
```

- [ ] **Step 5: 终批提交**

```bash
git add docs/plans/2026-10-12-auto-pes-container.md
git commit -m "纯势容器(Tier 1)实现收尾: 设计文档状态→已实现; 回归证据——双编译器干净重建绿/六案例 stdout(除 timing/注释行)对 results.txt 逐字节一致/T4(23.0605 仅 force_interface)T5(RANDOM_NUMBER 零命中)复验/check_auto_pes ALL PASS/drill_auto_pes 十七腿全绿"
```

---

## Self-Review（计划自审记录）

1. **Spec 覆盖**：设计 §1 词表/八中止 → Task 2（七条 prebuild 中止 + STOP 1 归 Task 4 leg 3）；§2 机械（织入/Makefile 两行/首编即绿）→ Task 2/3；§3 缺省语义与范本 → Task 1 fixture + Task 4 运行腿（MAX_STEPS 分支/R_TERM 分支/绝热惰性/arm 见证）；§5 判据（正负用例/零触碰/双编译器/六案例/T4T5）→ Task 4/5/7；§5.3 文档落点 → Task 6；D5 方法组 → Task 1 tdhf_stub + Task 2 配对中止 + Task 4 arm 双向见证。无缺口。
2. **占位符扫描**：无 TBD/TODO；全部代码块为终态（Task 5 crit_b 为单段正确 FD 循环）。
3. **名称/类型一致性**：模块名 auto_sample_pot/io/met 于 Task 1/2/4/5 一致；魔数名与设计 §1 表一致；`auto_sample_set_params` 定义（pot）与调用（io）成对；判据程序不引用 load/classify（拆文件动机自洽）；drill 断言串与 Task 2 Step 4 织入断言同源。
