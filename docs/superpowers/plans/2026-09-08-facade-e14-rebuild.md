# E14 Facade Rebuild Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild the facade layer as the four-entry blackbox contract (force_eval / elec_prop / mqc_statistic / termination), rebuild interface_* as the flat methods/ library, and land all E14 state/roster/propagator/driver revisions.

**Architecture:** Delete-and-rewrite of the comment-legislated skeleton (the project is in the K2 skeleton stage — every unit is a pseudocode-contract file; no executable bodies). Contract per spec §2-8; dependency direction force_facade → elec_facade → elec_ws; ws held as a module-private instance in elec_facade.

**Tech Stack:** Fortran 2008 free-form, comment-legislated skeletons; gfortran 15.2 + Intel ifx 2025.2 verification.

**Spec:** `docs/superpowers/specs/2026-09-08-facade-e14-design.md` — the plan argues from the spec; executors read both. Section references below ("spec §N") point into it.

## Global Constraints

- **Language form:** every unit is a skeleton — module head + `! Reference:` citations + signed procedures whose bodies are `! Pseudocode:` comment blocks. Match the existing style exactly (see `src/facade/elec_ws.f90` for module-head style, `src/facade/force_facade.f90` for registry + pseudocode style). All comments in English.
- **Compile verification (run after every task, from repo root, PowerShell):**
```powershell
$src = "src"; $tmp = "build\e14check"
Remove-Item -Recurse -Force $tmp -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force $tmp | Out-Null
$base = @('consts','state','control','config','rng','linalg','geometry','specio','elec_ws','config_roster','elec_facade')
$tail = @('elec_adiabatic','reg_adiabatic','force_facade','sysdef_facade','predicate','propagator','sampler','dist_boltzmann','dist_ebk','dist_incident','recorder','final_state','outcome','density','plot_hist','input','driver')
$files = @(); foreach ($m in ($base + $tail)) { $files += (Get-ChildItem -Recurse "$src\*.f90" | Where-Object { $_.BaseName -eq $m }).FullName }
Copy-Item $files -Destination $tmp
Push-Location $tmp; $ok = $true
foreach ($f in ($base + $tail)) { gfortran -c -J. "$f.f90"; if ($LASTEXITCODE -ne 0) { $ok = $false; break } }
Pop-Location; if ($ok) { "COMPILE_OK" } else { "COMPILE_FAIL" }
```
  Expected: `COMPILE_OK`. (Tasks 1-3 add files progressively — until a unit exists, drop it from `$base`/`$tail` in this command. From Task 5 on, `Get-ChildItem -Recurse` already picks up files from `src\methods\`.)
- **ifx verification (Task 14 only):** same order, `ifx -c -module . <file>.f90` with flags `-r8 -double-size=64 -i8 -O2 -w` in a fresh directory (the .mod formats are incompatible — never mix compilers in one directory).
- **Commit style:** one commit per task; message prefix `K2 收口(预审门): E14 …`; end with `Co-Authored-By: Claude Code <noreply@anthropic.com>`.
- **Do not touch** `venus2026/` (read-only old library) or `docs/plans/K1.md`, `P0.md`.

---

### Task 1: state.f90 — dual containers + occupation vector

**Files:**
- Modify: `src/globals/state.f90` (whole file rewrite in place, keeping existing header Reference lines)

**Interfaces:**
- Produces: `type(elec_t)` with fields `n_surf, a(:,:), rho(:,:)` (all `complex(8), allocatable` for the matrices); `type(state_t)` with `q/p/f/t`, `s`, `occ(:)` (integer, allocatable), `n_hop`; `state_create(sta, natoms [, n_surf])`, `state_destroy(sta)` unchanged signatures.
- All downstream units (elec_facade, force_facade, recorder, final_state, propagator, predicate, sampler, dist_*) consume these names.

- [ ] **Step 1: Rewrite the file to this skeleton (keep the existing `! Reference:` header block verbatim; extend the header with the E14 ruling note)**

```fortran
module state
   implicit none
   private
   public :: state_t, elec_t, state_create, state_destroy

   ! E14 (user rulings 2026-09-08): (1) amp(:) upgraded to the coefficient matrix a(:,:);
   ! (2) rho(:,:) added — symmetric dual ontic containers, no subordination legislated
   ! (state legislates only the EXISTENCE of two allocatable complex matrices; which one a
   ! method propagates is its own declaration — amplitude-ODE methods on a, density-matrix-
   ! dynamics methods on rho, adiabatic/LZSH/MDEF allocate neither); (3) i_surf generalized
   ! to the occupation vector occ(:) — a CLASSICAL-SIDE STATISTICAL field (the discrete
   ! output of mqc_statistic's stochastic sampling of the populations; under the
   ! single-electron representation multiple states are occupied simultaneously, so a single
   ! surface index is not general); (4) index-semantics statutes stay member-defined
   ! (E11/E12 route): a rows/columns = surface number / electron×state / configuration
   ! cardinality / mapping z = x+ip encoding — defined by the class-B method member.
   type :: elec_t
      integer :: n_surf = 1              ! (keep existing E12 semantics statute comment)
      complex(8), allocatable :: a(:,:)  ! coefficient matrix (propagation container; upgraded
                                         ! from the complex vector amp)
      complex(8), allocatable :: rho(:,:)! density matrix (ontic for rho-native methods;
                                         ! derived cache for a-native methods — the mapping
                                         ! form is method knowledge, refreshed by sync_proc)
   end type elec_t

   type :: state_t
      integer :: natoms = 0              ! (keep existing comments)
      integer :: ndof   = 0
      real(8), allocatable :: q(:)       ! [Å]
      real(8), allocatable :: p(:)       ! [amu·Å/(10 fs)]
      real(8), allocatable :: f(:)       ! [internal energy/Å]
      real(8) :: t = 0.0d0
      type(elec_t) :: s                  ! electronic ontic state
      integer, allocatable :: occ(:)     ! occupation vector — classical-side statistical
                                         ! field (spec §4): occ(k) = state of electron k
                                         ! (independent-electron methods); length 1
                                         ! (single-active-surface / configuration number);
                                         ! [1] constant (adiabatic). Written by
                                         ! mqc_statistic (a hop = one component change);
                                         ! read by force_eval (adiab-mode delivery-request
                                         ! source). Dimension n_e set at member init,
                                         ! reallocated together with a.
      integer :: n_hop = 0               ! hop counter — statistical observable (recorder-only
                                         ! read; relocated from elec_t, spec §4)
   end type state_t
contains
   ! state_create(sta, natoms [, n_surf]): allocates degenerate shapes — q/p/f(3*natoms),
   ! occ = [1], a/rho 1x1 when n_surf present. Pseudocode notes: larger ontic containers
   ! (e.g. IESH n_e x n_s) are reallocated by the method member-init hook (blackbox_init
   ! step, spec §4 "Allocation timing").
   ! state_destroy(sta): deallocates everything including occ/a/rho.
   ! (keep the existing pseudocode bodies of both routines, adjusting the allocation list)
end module state
```

- [ ] **Step 2: Compile** (drop not-yet-existing units from the verification list: `specio, elec_facade` out; `elec_adiabatic/force_facade/...` still compile against the old field set only if untouched — `force_facade` references `elec_t` by name only, so it still compiles; run the Global Constraints command with `$base`/`$tail` trimmed to existing files). Expected: `COMPILE_OK`.

- [ ] **Step 3: Commit**

```powershell
git add src/globals/state.f90
git commit -m "K2 收口(预审门): E14 T1 state 修订——amp→a(:,:) 双对称容器, i_surf→occ(:) 占据向量经典侧统计字段"
```

---

### Task 2: elec_ws.f90 — one-step lifetime + occupation-set request

**Files:**
- Modify: `src/facade/elec_ws.f90` (header axioms + ws_t field comments; no field changes except the request field)

**Interfaces:**
- Produces: `ws_t` with `repr`, request field generalized, `g_add/v_add/h_ij/dh_ij/v_surf/g_surf/d_ij` unchanged; `r_empty/r_scalar/r_diab/r_adiab` unchanged.

- [ ] **Step 1: Revise the header and field comments** — key legislated wording to land (replacing the E10 axiom wording where superseded):

  - Axiom B (black-box property): outside the box the four signatures of spec §3 are visible (`force_eval / elec_prop / mqc_statistic` / termination via predicate); ws never appears in any of them.
  - Axiom C (thing-in-itself): ws lifetime = ONE INTEGRATION STEP (E14: extended from one evaluation) — filled at force evaluation, consumed by elec_prop / mqc_statistic, reset to `r_empty` at step end. Holder: module-private instance of elec_facade (still invisible outside the box; declared in the box, dies in the box).
  - Request field: `i_want` generalized to the occupation-set request — `integer, allocatable :: i_want(:)` (the adiab-mode delivery request: the set of currently occupied single-particle states whose coupling rows are requested; length 1 = the single-active-surface degenerate form; written by force_eval from `sta%occ`, read by the member). Change the field declaration accordingly.
  - Keep the dependency-direction note (members use this module for type visibility; zero dependencies) and the "pure type/constant contract" closing note, updated: lifecycle management belongs to elec_facade (holder + reset), allocation/reuse per mode finalized in C2.

- [ ] **Step 2: Compile** (same trimmed list as Task 1). Expected: `COMPILE_OK`.

- [ ] **Step 3: Commit**

```powershell
git add src/facade/elec_ws.f90
git commit -m "K2 收口(预审门): E14 T2 elec_ws 修订——ws 生命期一步(跨三入口), 请求通道泛化为占据集 i_want(:)"
```

---

### Task 3: elec_facade.f90 — electronic phases + method-package registry + ws holder

**Files:**
- Create: `src/facade/elec_facade.f90`

**Interfaces:**
- Consumes: `state` (state_t/elec_t), `config` (electronic), `elec_ws` (ws_t, r_*).
- Produces: `elec_prop(sta, dt)`, `mqc_statistic(sta)`, `method_reg(name, prop_proc, match_proc [, sync_proc])`, `method_bind(name)` (assembly lookup called by blackbox_init), module-private `ws` instance with the accessor pair `ws_get()`/`ws_reset()` (subroutines, to avoid pointer plumbing), plus the abstract interfaces `elec_prop_i`, `mqc_match_i`, `sync_proc_i` that method-package members implement.

- [ ] **Step 1: Write the file to this skeleton** (full Reference/pseudocode style per `src/facade/force_facade.f90`; comment content argues from spec §3 entries 2-3, §4 sync contract, §5 registry):

```fortran
!=====================================================================
! elec_facade.f90 - L2 electronic-phase facade: entries 2-3 of the E14 four-entry contract
!   (elec_prop: §4.1 layer 2 electronic propagation; mqc_statistic: §4.1 layer 3
!   state-statistics matching) + the method-package registry + the ws instance holder.
! Reference: docs/superpowers/specs/2026-09-08-facade-e14-design.md §3 (entries 2-3,
!   ws holding), §4 (sync tail calls), §5 (method package registration);
!   docs/design_spec_completeness.md §4.1 (three-layer decomposition — layer semantics);
!   old-library propagation sites venus2026/src_VENUS/*.f (citation pass during C-phase
!   implementation — no skeleton-time citation).
!=====================================================================
module elec_facade
   use state,   only: state_t, elec_t
   use config,  only: electronic
   use elec_ws, only: ws_t, r_empty
   implicit none
   private
   public :: elec_prop, mqc_statistic, method_reg, method_bind, ws_get, ws_reset

   ! Method-package member interfaces (spec §5 — one registration fills the whole package;
   ! the registering author knows the pairing is legal — no combination table, E13 carried)
   abstract interface
      subroutine elec_prop_i(sta, dt)      ! layer 2: propagation member
         import :: state_t
         type(state_t), intent(inout) :: sta  ! reads ws (via ws_get) + sta%q; writes the
         ! method's ontic container (sta%s%a or sta%s%rho)
         real(8), intent(in) :: dt             ! [10 fs] (hop-probability channel — E9/E14)
      end subroutine elec_prop_i
      subroutine mqc_match_i(sta)          ! layer 3: matching member
         import :: state_t
         type(state_t), intent(inout) :: sta  ! writes sta%occ / sta%p (optionally sta%f via
         ! a force_eval self-call — intra-box, legal; optionally the ontic container for
         ! decay-of-mixing resets). Sole evolution-phase rng consumption point (§4.1).
      end subroutine mqc_match_i
      subroutine sync_proc_i(s)            ! container synchronization (spec §4)
         import :: elec_t
         type(elec_t), intent(inout) :: s  ! brings the method's declared deliverable
         ! containers up to date; direction/form (a→rho, rho→a, none) is method knowledge
      end subroutine sync_proc_i
   end interface

   ! Method-package table (single-select by name at assembly; registered via reg_*.f90 —
   ! spec §5 two-level registration; capacity 16)
   integer, parameter :: max_method = 16
   type :: method_entry_t
      character(len=32) :: name = ''
      procedure(elec_prop_i), pointer, nopass :: prop  => null()
      procedure(mqc_match_i), pointer, nopass :: match => null()
      procedure(sync_proc_i), pointer, nopass :: sync  => null()  ! null = no cache need (legal)
   end type method_entry_t
   type(method_entry_t) :: method_tbl(max_method)
   integer :: n_method = 0
   integer :: i_assembled = 0          ! row bound at blackbox_init from config%electronic%method

   ! The ws instance — E10 axiom C / E14: module-private holder, one-integration-step
   ! lifetime, invisible outside the box (never a dummy of any public signature)
   type(ws_t) :: ws

contains
   ! ws_get(w) — hand the instance to force_facade (fill site); members never call this
   ! (they receive ws through their signature — wait: prop/match members read ws how?
   ! spec §3: they read the held ws. Conveyance: elec_prop/mqc_statistic fetch the module
   ! instance and pass it into the member through an internal wrapper signature carrying
   ! ws_t — the member-facing abstract interfaces above take only (sta[,dt]); the wrapper
   ! layer is finalized in C2 together with the member signature freeze; the skeleton
   ! legislates: members read ws through the facade-internal channel, never as a public
   ! signature argument — E10 axiom A intact).
   ! ws_reset() — repr = r_empty (each step end; force_eval head re-request)
   ! method_reg(name, prop, match [, sync]) — fills one package row; duplicate name loud-stop
   !   (formal); the row is the atomic registration unit (spec §5)
   ! method_bind(name) — assembly-time name lookup into method_tbl (config%electronic%method);
   !   binds i_assembled; not found / table empty -> loud-stop (formal check; called from
   !   blackbox_init — the init choreography stays in force_facade as the single assembly
   !   entry point)
   ! elec_prop(sta, dt) — table dispatch to the assembled row's prop member; tail: sync
   !   call if attached (pure-form)
   ! mqc_statistic(sta) — table dispatch to the assembled row's match member; tail: sync
   !   call if attached (pure-form)
end module elec_facade
```

  Expand every `!` pseudocode note above into full `! Pseudocode:` blocks in the delivered file (steps per spec §3 entries 2-3: dispatch, member call, sync tail, loud-stops on unassembled table).

- [ ] **Step 2: Compile** (add `elec_facade` to the verification list after `elec_ws`; `force_facade` not yet rewritten — it does not use elec_facade yet, so order it before/after freely; keep trimmed list green). Expected: `COMPILE_OK`.

- [ ] **Step 3: Commit**

```powershell
git add src/facade/elec_facade.f90
git commit -m "K2 收口(预审门): E14 T3 elec_facade 新建——elec_prop/mqc_statistic 入口+方法包注册表+ws 实例持有"
```

---

### Task 4: force_facade.f90 — delete and rewrite (entry 1 + class-A registry)

**Files:**
- Delete + Create: `src/facade/force_facade.f90` (replace whole content)
- Modify: `src/dynamics/propagator.f90:30` and `src/output/recorder.f90:47` and `src/output/final_state.f90:41` — the `use force_facade, only: facade_force` lines become `use force_facade, only: force_eval` (comment references to `facade_force(` in those files also renamed; full pseudocode revisions of those files happen in Tasks 8/10-11 — here only the name).

**Interfaces:**
- Consumes: `elec_facade` (ws_get/ws_reset, method-table query for blackbox_init validation), `elec_ws`, `state`, `config`, `config_roster`, `consts`.
- Produces: `force_eval(sta, e_tot)`, `pes_reg(name, proc [, repr])`, `blackbox_init()`, `blackbox_info()`; `pes_proc_i` abstract interface `(q, mass, ws)` unchanged (E10 single signature).

- [ ] **Step 1: Rewrite force_facade.f90 to this skeleton** (carry over verbatim from the current file: the venus2026/DVDQ.f + interface_TEST.f90 Reference block, the ev_kcal unit-contract note, the pes_proc_i interface with its comments, the pes_entry_t table, the E10/E13/E14 header rulings — the E13 triad withdrawal note stays; add the E14 header note: three things the blackbox machine does are now three DRIVER-orchestrated entries, facade_force-the-orchestrator is withdrawn):

```fortran
module force_facade
   use consts,        only: e_conv
   use state,         only: state_t
   use config,        only: electronic
   use config_roster, only: roster
   use elec_ws,       only: ws_t, r_scalar, r_diab, r_adiab, r_empty
   use elec_facade,   only: ws_get, ws_reset, method_bind
   implicit none
   private
   public :: force_eval, pes_reg, blackbox_init, blackbox_info

   real(8), parameter :: ev_kcal = 23.0605d0   ! (keep T4 unique-point note)

   ! pes_proc_i abstract interface — CARRIED VERBATIM from the current file (E10 single
   ! signature, comments unchanged except: the request-field sentence now reads i_want(:)
   ! occupation set)

   ! pes_entry_t / pes_tbl / n_pes — CARRIED VERBATIM

contains
   ! pes_reg(name, proc [, repr]) — carried from facade_reg_pes (rename only; pseudocode
   !   unchanged: duplicate/full loud-stop; one registration line per member, T2)
   ! blackbox_init() — assembly + purely formal validation spanning BOTH facades (spec §3
   !   assembly entries):
   !   1. class-A assembly from pes_tbl per config%electronic%pes_main (carried pseudocode:
   !      empty table loud-stop; n_surf>=2 with r_scalar main PES loud-stop; additive rows
   !      r_scalar)
   !   2. method-slot assembly: call method_bind(config%electronic%method) (name lookup —
   !      formal; not found loud-stop). The table lives in elec_facade; the init
   !   choreography stays here as the single assembly entry point
   !   3. member-init hooks (POTPRE type — carried)
   !   4. observable-vs-container formal check (spec §4): an observable column demanding a
   !      container no assembled package member delivers -> loud-stop naming the conflict
   ! force_eval(sta, e_tot) — entry 1 (spec §3):
   !   1. zero g(:), v_sum; ws_reset(); dt removed from the signature (E9 channel withdrawn)
   !   2. main-PES evaluation: request occupation set — ws%i_want(:) = sta%occ (adiab-mode
   !      delivery request, generalized from the scalar; length 1 = degenerate);
   !      proc(sta%q, roster%mass, ws) via ws_get instance
   !   3. fold scalar contributions (carried pseudocode — member reset between calls)
   !   4. unit-contract fulfillment (carried verbatim: sta%f(i) = sta%f(i) - g(i)*ev_kcal*e_conv)
   !   5. NO electron-law consumption here (E14: withdrawn to elec_prop/mqc_statistic) —
   !      the ws of the current geometry is left in place for the later phases
   !   6. energy (carried: v_pot/t_kin/e_tot conventions)
   !   4'. multi-surface branch (carried; r_adiab takes the occupation-set rows
   !       g_surf(:, occ) — current-occupation force kernel; the U-delivery note carried,
   !       consumer re-anchored to elec_prop/mqc_statistic)
   ! blackbox_info() — carried (prints member lists, method name, assembly state)
end module force_facade
```

- [ ] **Step 2: Fix the three `use` lines** in propagator.f90 / recorder.f90 / final_state.f90 (`facade_force` → `force_eval`; also rename in-file comment call sites `facade_force(sta, e` → `force_eval(sta, e`). Do not otherwise touch those files.

- [ ] **Step 3: Compile** (full list minus not-yet-existing `sysdef_facade`, `specio`; `elec_adiabatic/dist_*` still at their old paths — the recursive Get-ChildItem picks them up wherever they live). Expected: `COMPILE_OK`.

- [ ] **Step 4: Commit**

```powershell
git add src/facade/force_facade.f90 src/dynamics/propagator.f90 src/output/recorder.f90 src/output/final_state.f90
git commit -m "K2 收口(预审门): E14 T4 force_facade 重写——force_eval 入口+blackbox_init 跨门面装配, 电子相位剥离, 三处 use 行改名"
```

---

### Task 5: methods/ library — directory rebuild + elec_adiabatic package

**Files:**
- Create: `src/methods/` directory
- Move: `src/interface_electronic/elec_adiabatic.f90` → `src/methods/elec_adiabatic.f90` (rewrite module to the method-package form)
- Create: `src/methods/reg_adiabatic.f90` (thin reg module)
- Move: `src/interface_electronic/README.md` → `src/methods/README.md` (rewrite for the library role; fold in interface_distributions content in Task 6)
- Delete: `src/interface_electronic/` (after emptying)

**Interfaces:**
- Consumes: `state` (elec_t/state_t), `elec_ws` (ws_t), `elec_facade` (elec_prop_i/mqc_match_i abstract interfaces, method_reg).
- Produces: `elec_adiabatic` module with degenerate members conforming to `elec_prop_i`/`mqc_match_i` (the adiabatic package: prop = no-op, match = no-op, no sync); `reg_adiabatic` module with `reg_adiabatic_all()` calling `method_reg('adiabatic', ad_prop, ad_match)`.

- [ ] **Step 1: Rewrite elec_adiabatic.f90 as the degenerate method package** — module head cites spec §5 (library), §5 two-exposure-levels note; members:

```fortran
module elec_adiabatic
   use state,       only: state_t
   use elec_facade, only: elec_prop_i, mqc_match_i   ! contract documentation only (the
   ! members conform to these interfaces; conformance is checked at the method_reg call
   ! site in reg_adiabatic.f90 — E14 member-signature freeze lands in C2)
   implicit none
   private
   public :: ad_prop, ad_match
   ! (the old ad_dsdt/ad_force_rule adapter content is superseded: E14 splits the
   !  electronic phases out of force composition — the adiabatic package is the all-no-op
   !  degenerate point; keep an E9-note that the degenerate-equivalence verification
   !  criterion (C2) now compares four-entry orchestration with pure PES)
contains
   subroutine ad_prop(sta, dt)   ! elec_prop_i: no-op (no amplitudes under adiabatic)
      type(state_t), intent(inout) :: sta; real(8), intent(in) :: dt
      ! Pseudocode: 1. return (no-op; ignores ws and dt — E13 wording carried)
   end subroutine ad_prop
   subroutine ad_match(sta)      ! mqc_match_i: no-op (occ stays [1]; no rng use)
      type(state_t), intent(inout) :: sta
      ! Pseudocode: 1. return (no-op — single surface, nothing to sample)
   end subroutine ad_match
end module elec_adiabatic
```

  (Expand pseudocode/comment statutes in delivery; unused-dummy notes per compiler: add `associate` nothing — comment-only bodies are the project norm.)

- [ ] **Step 2: Write reg_adiabatic.f90**

```fortran
! reg_adiabatic.f90 — thin library-level registration module (spec §5: the reg file is
! the atomic registration unit; prebuild aggregates methods/*/reg_*.f90 + systems/*/reg_*.f90
! — K3 scope; the K2 skeleton wires it from the driver pseudocode transitionally)
module reg_adiabatic
   use elec_facade,  only: method_reg
   use elec_adiabatic, only: ad_prop, ad_match
   implicit none
   private
   public :: reg_adiabatic_all
contains
   subroutine reg_adiabatic_all()   ! one call fills the whole package (spec §5)
      ! Pseudocode: 1. call method_reg('adiabatic', ad_prop, ad_match)
   end subroutine reg_adiabatic_all
end module reg_adiabatic
```

- [ ] **Step 3: Write methods/README.md** — the library charter: flat built-in method library; file prefixes self-describe the family (`elec_*` method packages, `dist_*` distribution members, `amp_*` shared amplitude integrators phase C); two exposure levels (whole packages via config key / reusable components via use — spec §5 "Two exposure levels" verbatim in summary); drag-in rule (new method = new file + thin reg_*.f90; Makefile wildcard picks it up; prebuild aggregates).
  Delete `src/interface_electronic/` (verify empty first: `Get-ChildItem src\interface_electronic` returns nothing).

- [ ] **Step 4: Compile** (`elec_adiabatic` and `reg_adiabatic` in the list — `reg_adiabatic` compiles but nothing uses it yet; place it after `elec_adiabatic`). Expected: `COMPILE_OK`.

- [ ] **Step 5: Commit**

```powershell
git add src/methods/ ; git add -A src/interface_electronic/
git commit -m "K2 收口(预审门): E14 T5 methods/ 库重建——elec_adiabatic 方法包化+reg_adiabatic 薄注册+库章程 README"
```

---

### Task 6: methods/ library — dist_* migration

**Files:**
- Move: `src/interface_distributions/dist_boltzmann.f90`, `dist_ebk.f90`, `dist_incident.f90` → `src/methods/` (module content unchanged — signatures untouched; only the header gains one line: "Library home: methods/ (E14) — a dist_* member is a pre-stored generic distribution method; its samp_reg registration line lands at C5 together with its sample routine (K2 E3 slimming ruling); see methods/README.md")
- Modify: `src/methods/README.md` (created in Task 5)
- Delete: `src/interface_distributions/` (verify empty first)

**Interfaces:**
- Consumes: nothing new (the modules' internals are untouched).
- Produces: unchanged `dist_boltzmann/dist_ebk/dist_incident` modules at their new home.

**NOTE — no reg file for dist_* at K2:** `samp_reg` takes a `procedure(dist_proc_i)` — a `dist_proc_i(sta)` sample routine. Per the K2 E3 slimming ruling the three dist skeletons carry only `*_init(params)` (sample routines land at C5), so a `reg_dist_seed` module calling `samp_reg` would not typecheck today. Each dist member's thin `reg_dist_*.f90` therefore lands at C5 together with its sample routine (drag-in rule unchanged). `reg_adiabatic.f90` is the sole reg file of the K2 library seed.

- [ ] **Step 1: Move the three files** (git mv or copy+delete), adding the one-line header note to each.

- [ ] **Step 2: Extend methods/README.md** — add the dist family section: `dist_*` members (pre-stored generic initial-state distributions, parameters via config — E12 wording carried); registration timing note (samp_reg lines + thin reg files land at C5 with the sample routines, per K2 E3); the remaining 9 dist members listed as C5 arrivals.

- [ ] **Step 3: Compile** (dist files at their new path — the recursive Get-ChildItem picks them up from `src\methods\`). Expected: `COMPILE_OK`.

- [ ] **Step 4: Commit**

```powershell
git add -A src/methods/ src/interface_distributions/
git commit -m "K2 收口(预审门): E14 T6 dist 三件迁入 methods/——interface_distributions 撤编, 注册推迟至 C5(E3 裁定)"
```

---

### Task 7: specio.f90 + sysdef_facade.f90 — system-definition files

**Files:**
- Create: `src/utils/specio.f90` (format readers + element mass table)
- Create: `src/facade/sysdef_facade.f90` (discovery / dispatch / roster fill / paradigm derivation / sampler hand-off)

**Interfaces:**
- Consumes: `config_roster` (roster), `config` (reactants%surface_model), `state` (nothing — operates on roster data).
- Produces: `read_poscar(u, title, cell, names, xyz, ok [, errmsg])`, `read_xyz(u, names, xyz, ok [, errmsg])`, `elem_mass(symbol) -> real(8)` [amu]; `sysdef_load(scan_dir)`.

- [ ] **Step 1: Write specio.f90 skeleton**

```fortran
!=====================================================================
! specio.f90 - L4 system-definition file readers (generic, zero system knowledge — E14
!   spec §6: parsers/mass table/derivation rule are generic code; file CONTENTS are
!   blackbox material). Formats: POSCAR (VASP, cell + fractional/cartesian coords) and
!   .xyz (no cell); new formats = new reader members (grammar versioning).
! Reference: spec §6 ownership table; POSCAR/.xyz public format definitions.
!=====================================================================
module specio
   implicit none
   private
   public :: read_poscar, read_xyz, elem_mass
contains
   ! read_poscar(u, title, cell, names, xyz, ok [, errmsg]) — parse one POSCAR: comment
   !   line, scale factor (applied to cell AND coords per VASP semantics), 3 lattice
   !   vectors, element-name line, counts line, Direct/Cartesian tag, coordinates.
   !   Pseudocode: line-by-line read; malformed -> ok=.false. + errmsg
   ! read_xyz(u, names, xyz, ok [, errmsg]) — parse one .xyz frame: count line, comment
   !   line, "Symbol x y z" rows. (No cell — by format definition.)
   ! elem_mass(symbol) — periodic-table lookup [amu] (generic physical-constant table;
   !   unknown symbol -> loud-stop naming it)
end module specio
```

- [ ] **Step 2: Write sysdef_facade.f90 skeleton**

```fortran
!=====================================================================
! sysdef_facade.f90 - L2 system-definition facade (E14 spec §6): discovery, parsing
!   dispatch, roster fill, paradigm derivation, sampler hand-off. NOT an evolution-loop
!   entry — a data-material contract consumed once in the assembly phase.
! Reference: spec §6 (granularity/discovery/formats/derivation/ownership/data flow);
!   venus2026/src_VENUS/SELECT.f90:1300-1334 setup_surface_collision / :1337-1404
!   setup_gas_collision (the two routes the derived paradigm selects — citation anchor).
!=====================================================================
module sysdef_facade
   use config,        only: reactants
   use config_roster, only: roster
   use specio
   implicit none
   private
   public :: sysdef_load
contains
   ! sysdef_load(scan_dir) — spec §6 data flow, pseudocode:
   ! 1. Scan scan_dir (the blackbox folder, config-side path) for known extensions
   !    (POSCAR/.poscar/.xyz); file base name = reactant name; multiple files per
   !    extension legal (one reactant each); zero files -> loud-stop (no system defined)
   ! 2. Per file dispatch: .poscar/POSCAR -> read_poscar (cell retained); .xyz -> read_xyz
   !    (no cell). A reactant given in both formats -> loud-stop naming it
   ! 3. Masses: elem_mass(names) per atom
   ! 4. Fill roster: fragment tables (one fragment per file — reactant name), equilibrium
   !    geometries, unit cell (surface systems)
   ! 5. Paradigm derivation (generic formal rule, spec §6): any cell present ->
   !    reactants%surface_model = surface; none -> molecular
   ! 6. Hand-off note: equilibrium geometries + fragment masses are now roster data —
   !    sampler B' members read them through the roster (no module dependency on sampler)
   ! 7. Validation loud-stops: mixed surface/molecular definitions across files
   !    (some with cell, some without) -> loud-stop naming the conflict
end module sysdef_facade
```

- [ ] **Step 3: Compile** (add `specio` after `geometry`, `sysdef_facade` after `force_facade`). Expected: `COMPILE_OK`.

- [ ] **Step 4: Commit**

```powershell
git add src/utils/specio.f90 src/facade/sysdef_facade.f90
git commit -m "K2 收口(预审门): E14 T7 sysdef 门面+specio 解析器——体系定义文件契约(发现/解析/填充/范式派生/交接)"
```

---

### Task 8: propagator.f90 — hooks withdrawal + endpoint-force contract

**Files:**
- Modify: `src/dynamics/propagator.f90` (whole-file revision)

**Interfaces:**
- Consumes: `state`, `config`, `config_roster`, `force_facade` (`force_eval` — renamed in Task 4).
- Produces: `prop_reg(id, proc)`, `prop_step(sta)` — single argument; `integ_proc_i(sta)`; `hooks_t`/`hook_out_i`/`hook_term_i` DELETED.

- [ ] **Step 1: Revise the file** — delete the hooks block (`hook_out_i`, `hook_term_i`, `hooks_t`, the second prop_step dummy, all hooks pseudocode); `integ_proc_i(sta)` single-dummy; prop_step pseudocode becomes: dispatch per `config%propagator%id`, member advances sta by `prop_cfg%dt`, **endpoint-force contract in the header** (spec §3: "prop_step returns only after the endpoint force has been evaluated — ws corresponds to the returned sta%q; under MQC the VV right-end force cannot be folded into the next step's left end; the amplitude propagation of phase 2 needs the endpoint H"); frame output / termination move to the driver (E7 revision note). Header gains the E14 ruling note (E7 "only direct call is prop_step" withdrawn; propagator = pure integrator; radau fixed-step note on multi-surface paths carried, E13).

- [ ] **Step 2: Compile**. Expected: `COMPILE_OK`.

- [ ] **Step 3: Commit**

```powershell
git add src/dynamics/propagator.f90
git commit -m "K2 收口(预审门): E14 T8 propagator 回归纯积分器——hooks 撤编, 终点力契约立法"
```

---

### Task 9: predicate.f90 — primitive registry (two-track composition)

**Files:**
- Modify: `src/dynamics/predicate.f90`

**Interfaces:**
- Consumes: `state`, `config_roster`.
- Produces: existing grammar (`pred_parse`/`pred_eval`/four primitives/three combinators/`k_*` codes) + `pred_reg_prim(name, proc)`, new code `k_prim = 8`, `prim_proc_i` abstract interface (`logical function prim_proc_i(sta)`), leaf-node fields `prim_name` on `expr_t`.

- [ ] **Step 1: Add to the module** (rest of the file unchanged): the abstract interface, a primitive table (`max_prim = 16`, `prim_entry_t{name, procedure(prim_proc_i), pointer, nopass :: proc}`, `n_prim`), `k_prim = 8` parameter, `prim_name` character field on `expr_t` (used when `node_kind == k_prim`), and:

```fortran
   ! pred_reg_prim(name, proc) — register a SYSTEM termination primitive (spec §7 track 1:
   ! system knowledge lives in the blackbox; the leaf set opens, the boolean combinator
   ! set stays closed — N2 revision). Duplicate name or full table -> loud-stop.
   ! pred_parse: grammar gains primitive := ... | PRIM('name') — the leaf matches the
   !   registered table at parse time (unknown name -> grammar error naming it)
   ! pred_eval: k_prim nodes look up the table and call the member proc with sta
```

  Expand pseudocode blocks accordingly; update the module header grammar block and the N2 citation with the E14 two-track note.

- [ ] **Step 2: Compile**. Expected: `COMPILE_OK`.

- [ ] **Step 3: Commit**

```powershell
git add src/dynamics/predicate.f90
git commit -m "K2 收口(预审门): E14 T9 predicate 双轨组合——原语开放注册集+PRIM 叶节点, 布尔组合子闭集不变"
```

---

### Task 10: recorder.f90 + final_state.f90 — rho observables + rename sweep

**Files:**
- Modify: `src/output/recorder.f90`, `src/output/final_state.f90`

**Interfaces:**
- Consumes: `state` (occ, rho), `force_facade` (`force_eval`).
- Produces: unchanged public faces (`rec_select/rec_frame/rec_final`, `fin_t/fin_eval`).

- [ ] **Step 1: recorder.f90** — comment statutes: population/coherence observable columns read `sta%s%rho` diagonals/off-diagonals and never touch the semantics of `a` (spec §4); the `n_hop` column reads `sta%n_hop`; the `i_surf`-era column comments become occupation-vector comments (`sta%occ`); energy columns via `force_eval(sta, e)` (renamed in Task 4 — verify no stragglers).

- [ ] **Step 2: final_state.f90** — same sweep (`facade_force` → `force_eval` in remaining comments; `sta%s%i_surf`-era references → `sta%occ`).

- [ ] **Step 3: Compile**. Expected: `COMPILE_OK`.

- [ ] **Step 4: Commit**

```powershell
git add src/output/recorder.f90 src/output/final_state.f90
git commit -m "K2 收口(预审门): E14 T10 recorder/final_state 观测列改读 rho/occ——force_eval 名字清扫"
```

---

### Task 11: config / config_roster / input / sampler — derivation and grammar notes

**Files:**
- Modify: `src/globals/config.f90`, `src/globals/config_roster.f90`, `src/input.f90`, `src/facade/sampler.f90`

**Interfaces:**
- Produces: no signature changes — comment statutes only.

- [ ] **Step 1: config.f90** — `electronic%method` comment: name-string lookup into the methods/ table (`'adiabatic'` seed; `'fssh'/'iesh'` phase C); `reactants%surface_model` comment: DERIVED by sysdef_load from cell presence (no longer an input key — spec §6); note the blackbox-folder path field (name `system_dir`, grammar finalized in C4).

- [ ] **Step 2: config_roster.f90** — `roster_load` comment: data source = sysdef facade products (was input keywords); the E12 a_lat statute carried (0.0 default + loud-stop).

- [ ] **Step 3: input.f90** — grammar notes: TERMINATE gains the `PRIM('name')` leaf (registered primitives); `surface_model` key withdrawn from the grammar (derived); `method`/`pes_main` remain name strings.

- [ ] **Step 4: sampler.f90** — header comment: initial geometries / fragment data arrive from the sysdef hand-off (spec §9: sampler kept, comment only — B' members read equilibrium geometries through the roster; no module dependency on sysdef_facade).

- [ ] **Step 5: Compile**. Expected: `COMPILE_OK`.

- [ ] **Step 6: Commit**

```powershell
git add src/globals/config.f90 src/globals/config_roster.f90 src/input.f90 src/facade/sampler.f90
git commit -m "K2 收口(预审门): E14 T11 config/roster/input/sampler 注释立法——method 名串查表, surface_model 派生, PRIM 语法, sysdef 交接"
```

---

### Task 12: driver.f90 — four-phase orchestration

**Files:**
- Modify: `src/driver.f90`

**Interfaces:**
- Consumes: everything produced by Tasks 3-9 (`force_eval`, `elec_prop`, `mqc_statistic`, `sysdef_load`, `blackbox_init`, no `hooks_t`).
- Produces: the program skeleton (no library consumers beyond it).

- [ ] **Step 1: Revise the pseudocode** per spec §8 (keep the excision-list header block and the three-phase structure; changes):
  - use list: drop `hooks_t`; `force_facade` now `force_eval, pes_reg, blackbox_init`; add `elec_facade: elec_prop, mqc_statistic`, `sysdef_facade: sysdef_load`, `reg_adiabatic` (transitional registration form until K3 prebuild; note E11; dist registration absent — C5 per Task 6 note).
  - Phase 1: insert `sysdef_load(system_dir)` after `roster_load` (order note: sysdef fills the roster, so roster_load shrinks to allocation/algebra — pseudocode note); registration step cites `reg_adiabatic_all()` + `pes_reg` lines; `blackbox_init()`; delete the hooks wiring step; `pred_parse` unchanged.
  - Phase 2 loop body per spec §8 verbatim: `prop_step(sta)` → `elec_prop(sta, dt_eff)` → `mqc_statistic(sta)` → `rec_frame(sta, u_out)` → `term_flag = pred_eval(sta, expr_term, step_i)`; delete `hooks` locals; add `dt_eff` local (from `prop_cfg%dt` — fixed-step source note).
  - E7 phase-restriction paragraph rewritten: direct calls = prop_step + facade phase entries + rec_frame + pred_eval.

- [ ] **Step 2: Compile** (full tree). Expected: `COMPILE_OK`.

- [ ] **Step 3: Commit**

```powershell
git add src/driver.f90
git commit -m "K2 收口(预审门): E14 T12 driver 四相位编排——prop_step→elec_prop→mqc_statistic→rec_frame→pred_eval, hooks 撤编"
```

---

### Task 13: documentation — READMEs, K2.md E14 entry, master plan

**Files:**
- Modify: `src/README.md`, `docs/plans/K2.md`, `docs/final_completion_plan.md`
- (methods/README.md already created in Task 5.)

**Interfaces:** none (documentation).

- [ ] **Step 1: src/README.md** — file tree updated (facade 5 files, methods/, interface_* gone); the call-graph note updated to the four-phase loop.

- [ ] **Step 2: K2.md** — append the canonical **E14** ruling entry after E13 (style: match the E10-E13 entries — ruling basis in the user's terms, the withdrawal/supersession ledger from spec §"Ruling revision ledger", the open points carried). Also update K2.md's "File tree" section if it lists interface_* paths.

- [ ] **Step 3: final_completion_plan.md §5 N1** — supersession note: E14 four-entry contract replaces the single-entry orchestration; the negotiation/registration wording points to the spec.

- [ ] **Step 4: Compile sanity** (docs only — rerun the gfortran verification anyway to confirm the tree is still green). Expected: `COMPILE_OK`.

- [ ] **Step 5: Commit**

```powershell
git add src/README.md docs/plans/K2.md docs/final_completion_plan.md
git commit -m "K2 收口(预审门): E14 T13 文档——K2.md E14 裁决条目, 主计划 §5 N1 supersession, README 树更新"
```

---

### Task 14: closure verification — dual compiler + whole-tree audit

**Files:** none modified (verification only; fixes land as follow-up edits within this task if needed).

- [ ] **Step 1: gfortran full-tree compile** (Global Constraints command, complete 28-unit list). Expected: `COMPILE_OK`.

- [ ] **Step 2: ifx full-tree compile** — fresh directory, same order, `ifx -r8 -double-size=64 -i8 -O2 -w -c -module . <f>.f90`. Expected: exit 0 (ifx warning #6843 on empty intent(out) bodies is known-benign — progress.md precedent).

- [ ] **Step 3: Whole-tree comment-consistency audit** — grep sweeps, all must return zero NEW hits (historical citations in Reference lines are legitimate; body/comment references are not):
  - `facade_force|facade_reg_pes|facade_reg_elec|facade_init|facade_info|facade_elec_prop|facade_mqc` outside `docs/` historical records
  - `amp\(` field references in src/ (the vector is gone — only `a(:,:)` remains)
  - `i_surf` in src/ (only historical statute citations naming the withdrawn field remain; no live field uses)
  - `hooks_t|hook_out_i|hook_term_i` in src/
  - `interface_electronic|interface_distributions` in src/ and docs/plans/K2.md's live tree section
  - `use force_facade, only: facade_force` anywhere
- [ ] **Step 4: spec-coverage crosscheck** — walk spec §"Ruling revision ledger" (9 rows) + §2-§9 and confirm each landed in Tasks 1-13; list any gap and fix it in this task.

- [ ] **Step 5: Final commit (if fixes landed) + closure report** in the conversation: units compiled, both compilers, audit sweep results, task-to-spec mapping.

```powershell
git add -A src/ docs/
git commit -m "K2 收口(预审门): E14 T14 终验——gfortran+ifx 双绿, 全树注释一致性审查通过, E14 收口"
```

---

## Post-plan notes for executors

- **Style calibration:** before Task 1, read `src/facade/elec_ws.f90` and the current `src/facade/force_facade.f90` in full — they are the gold standard for header rulings, `! Reference:` citations, `(see previous line)` continuation convention, and `! Pseudocode:` bodies.
- **Naming:** entry names are physical (force_eval/elec_prop/mqc_statistic/blackbox_init/blackbox_info/sysdef_load/pes_reg/method_reg/pred_reg_prim); module names keep *_facade; library files keep elec_*/dist_* prefixes; reg files are `reg_<package>_all` modules.
- **The skeleton norm:** comment-only bodies with signed dummies; a dummy that a pseudocode routine does not use gets a trailing comment (see `ad_prop` in Task 5), not `associate` tricks.
