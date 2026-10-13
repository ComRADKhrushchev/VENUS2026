# 自定义势能面接入说明（Tier 1 纯势容器）

最短路径：`systems/` 下一个文件夹、若干 `.f90`、**零 reg 文件**——构建系统
（scripts/prebuild.sh）按十个**预留名**（reserved names）自动接线，`make`
一次即完成。完整裁决记录
[docs/plans/2026-10-12-auto-pes-container.md](plans/2026-10-12-auto-pes-container.md)
（其中原称"魔数名"）；全部扩展点（Tier 2 全容器、采样成员、积分器等）见
[src/EXTENDING.md](../src/EXTENDING.md)。

## 预留名是什么

Tier 1 容器不写注册文件。prebuild 扫描文件夹内 `.f90` 的 public 过程名，按
名字接线：叫 `pes` 的接受力槽，叫 `keys` 的接键表声明，以此类推——这十个
名字即**预留名**（词表封闭；名字即约定，签名即既有接口，编译期核对）。
名字可居文件夹内任意模块，入口须 `public`（private 即编译错）。

文件夹内放了 `reg_*.f90` 即 Tier 2 全容器，预留名全部失效（reg 权威）。
反向升级零迁移：预留名函数可被 reg 直接复用。

## 文件夹地址

- 容器：`systems/<体系名>/`（如 `systems/mysys/`），`.f90` 放文件夹根。
- 构型（如随容器走）：`systems/<体系名>/system/`，每片段一个 `.xyz` 或
  `.poscar`（文件名干 = 片段名，名序即原子布局；文件只带几何，体系胞由输入 CELL 键给）。案例也可自带 case 侧
  `system/` 由 `SYSTEM_DIR` 指向——两条通道任选。
- 激活 = 文件夹位于 `systems/` 顶层（搬入即激活）；不用时停回
  `systems/inactive/`。一次只激活一个容器；与 reference 并存会在装配期中止
  （测自动容器时把 reference 暂存 `systems/inactive/`）。

## 基本规范

1. 力交付契约：`q`[Å] / `mass`[amu] 进，`g = dV/dq`[eV/Å] 与 `v`[eV]
   **全赋值**出；单位换算归力接口，勿自换。
2. 参数两条路：源内常数，或容器键（`keys` 声明 + `load` 拉取；键名不得与
   程序键/成员键撞名，声明即中止）。
3. 构建：`cd src && make`——prebuild 自动检测并织入（拖入即 make，首编即绿）。
   出错为指名中止：`pes` 缺失或 ≥2、任一名 ≥2、形不符（如 `subroutine keys`）、
   双自动文件夹、方法组不成对，均指名报错。
4. 范本：`systems/inactive/auto_sample/`（pot.f90 / io.f90 / tdhf_stub.f90；
   判据 `scripts/check_auto_pes.f90` + 演练 `scripts/drill_auto_pes.sh`）——
   拷走改物理本体即可。

下两节按预留名逐名组织：每名一小节，先一行说明（形 + 什么时候写），后跟
一块示例代码。示例贯穿同一个谐振子小例（弹簧常数 `K_HARM`、终止球
`R_TERM`，与范本 auto_sample 同构），各节代码可直接拼成一个能跑的容器；
方法包组示例为合法最简形（无操作骨架），物理本体随需替换。

## 势能面组（默认只写这组）

### `pes` —— 力交付入口（唯一必写，恰一个）

```fortran
subroutine pes(q, mass, g, v)
   real(8), intent(in)  :: q(:), mass(:)   ! q/mass 一维展平, size = 3×nat
   real(8), intent(out) :: g(:), v         ! g = dV/dq [eV/Å], v = V [eV], 全赋值
   v = 0.5d0*k_harm*dot_product(q, q)      ! 例: 各原子向原点的弹簧井
   g = k_harm*q
end subroutine pes
```

### `keys` —— 容器键词表声明（有容器键时写）

形：function 返回 `character(len=16)` 数组；解析前声明。

```fortran
function keys() result(kk)
   character(len=16) :: kk(2)
   kk = [character(len=16) :: 'K_HARM', 'R_TERM']
end function keys
```

### `load` —— 参数拉取 / 初始化钩子（有容器键或需初始化时写）

形：无参 subroutine，模块自持 `use input, only: pull`。读参数文件（如 NN
权重）、分配等初始化也在此（seam 期，先于一切力调用）；缺键保留缺省，
`found` 带出键在场。

```fortran
subroutine load()
   use input, only: pull
   call pull('K_HARM', k_harm)              ! 缺键: k_harm 保留模块缺省值
   call pull('R_TERM', r_term, r_found)     ! r_found: R_TERM 是否在场
end subroutine load
```

### `term` —— 终止判据（可选）

形：function 返回 `logical`。未绑 = 跑满 `MAX_STEPS`；键缺席门控用 `found`
旗标自实现。

```fortran
logical function term(sta)
   use state, only: state_t
   type(state_t), intent(in) :: sta
   term = .false.
   if (r_found) term = any(abs(sta%q) > r_term)   ! R_TERM 在场才启用
end function term
```

### `classify` —— 末态分类（可选）

形：function 返回 `character(len=32)` 通道名。

```fortran
function classify(fin) result(name)
   use final_state, only: fin_t
   type(fin_t), intent(in) :: fin
   character(len=32) :: name
   name = 'bound'                                  ! 例: 恒返 bound (fin 未读)
end function classify
```

## 非绝热方法包组（默认不写）

只在做表面跳跃类非绝热动力学时写（同文件夹另一模块）。注册无条件、初始化
有门：不设 `ELEC_METHOD` 时方法行惰性、`init` 不调，绝热参考面照常可用。

### `prop` —— 电子传播成员（与 `match` 成对）

```fortran
subroutine prop(sta, dt)
   use state, only: state_t
   type(state_t), intent(inout) :: sta
   real(8), intent(in) :: dt
   ! ……电子态传播……
end subroutine prop
```

### `match` —— 统计匹配成员（与 `prop` 成对）

跳变判断，occ/n_hop 写回。

```fortran
subroutine match(sta)
   use state, only: state_t
   type(state_t), intent(inout) :: sta
   ! ……跳变判断; occ/n_hop 写回……
end subroutine match
```

### `sync` —— 表象换算尾巴（可选随对）

a↔rho 表象换算。

```fortran
subroutine sync(s)
   use state, only: elec_t
   type(elec_t), intent(inout) :: s
   ! ……a ↔ rho 表象换算……
end subroutine sync
```

### `init` —— 选中后的运行期初始化（可选随对）

分配/参数落位；是否调用门控归生成代码（未被 `ELEC_METHOD` 选中即不调）。

```fortran
subroutine init()
   ! ……分配/参数落位……
end subroutine init
```

### `method_name` —— 注册名（随对必备）

返回 `character(len=32)` **小写**名；`ELEC_METHOD = <名>` 选用。

```fortran
function method_name() result(s)
   character(len=32) :: s
   s = 'my_tdhf'
end function method_name
```

## 升级注记

哪天需要自定义记录列、一个文件夹多个方法包或更复杂的参数体——同文件夹加
`reg_*.f90` 即升级 Tier 2（见 src/EXTENDING.md §1b），上述函数全部直接
复用，一字不动。
