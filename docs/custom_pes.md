# 自定义势能面接入说明(2026/10/14)

> 本文件与 [input_output.md](input_output.md) 并列为本程序现行仅有的两份正式使用说明。

- 路径：`<程序根目录>/systems/` 内；
- 基本形式：数量、名称不限的若干 `.f90`，只需要包含特定名称的子例程即可在编译时由`prebuild.sh`自动识别；详情见下。
- 子例程的参数形式必须与示例一致，如`pes`必须为 `pes(q,mass,g,v)`



## 势能面与体系

### `pes` —— 提供能量/力（唯一必须）

```fortran
subroutine pes(q, mass, g, v)
   real(8), intent(in)  :: q(:), mass(:)   ! q/mass 一维展平, size = 3×nat
   real(8), intent(out) :: g(:), v         ! g = dV/dq [eV/Å], v = V [eV]

   v = 0.5d0*k_harm*dot_product(q, q)      ! 例: 各原子向原点的弹簧井
   g = k_harm*q

   !....[可以在这里接入任意其他自定义势能面]

end subroutine pes
```
### `term` —— 终止判据（多数情况下必须）

形：function 输入物理状态`sta`，返回 `logical`;
- 轨迹会在满足规定的条件(`term == .true.`)时终止。如不提供，则每条轨迹都仅会跑到`MAX_STEPS`时终止；

```fortran
logical function term(sta)
   use state, only: state_t
   type(state_t), intent(in) :: sta
   term = .false.
   if (r_found) term = any(abs(sta%q) > r_term)   ! R_TERM 在场才启用
end function term
```

### `classify` —— 末态分类（可选）

形：function 输入末态`fin`，返回 `character(len=32)` 通道名。

```fortran
function classify(fin) result(name)
   use final_state, only: fin_t
   type(fin_t), intent(in) :: fin
   character(len=32) :: name
   name = 'bound'                                  ! 例: 恒返 bound (fin 未读)
end function classify
```

### 数据结构参考

两处入参都是程序内部结构体：`state_t` 定义在
[src/globals/state.f90](../src/globals/state.f90)，`fin_t` 定义在
[src/output/final_state.f90](../src/output/final_state.f90)。
```fortran
type :: state_t                       ! 纯物理态 Q/P/s/F/t
   real(8), allocatable :: q(:)       ! 坐标 [Å]（size = ndof = 3×natoms）
   real(8), allocatable :: p(:)       ! 动量 [amu·Å/fs]
   real(8), allocatable :: f(:)       ! 力 [amu·Å/fs²]
   real(8) :: t                       ! 轨迹时间 [fs]（driver 持有）
   type(elec_t) :: s                  ! 电子本征态（定义见非绝热方法一节）
end type state_t
```

`sta` 在动力学中持续演化，`fin`则是轨迹终止时装配的末态量：

```fortran
type :: fin_t                         ! 末态量
   real(8) :: e_int_a(3), e_int_b(3)  ! 两片段内能 [eV]：内动能/对参考势能/总量
   real(8) :: j_a(3), j_b(3)          ! 两片段角动量矢量 [ħ]
   real(8) :: j_mag_a, j_mag_b        ! 两片段角动量模 [ħ]
   real(8) :: e_rot_a, e_rot_b        ! 两片段转动能 [eV]
   real(8) :: e_rel, e_rel_sq         ! 相对平动能 [eV]：径向投影 / 全速度
   real(8) :: oam(3), oam_mag         ! 轨道角动量矢量 / 模 [ħ]
   real(8) :: b_fin                   ! 末态碰撞参数 [Å]
   real(8) :: e_cm                    ! 整体质心平动能 [eV]
   real(8) :: ang(20)                 ! 散射角族 [deg]（1-3 末值对、4-16 初值组
                                      !   合、17-18 表面入射、19-20 表面散射；
                                      !   其中 4-18 依赖初值快照、暂无生产者）
   real(8) :: t_life                  ! 轨迹寿命 [fs]（= sta%t）
   real(8), allocatable :: q_fin(:)   ! 末态 COM 系坐标归档 [Å]（简正分析用）
   real(8), allocatable :: p_fin(:)   ! 末态 COM 系动量归档 [amu·Å/fs]
   real(8), allocatable :: q_raw(:)   ! 末态原始系坐标 [Å]（全原子，轨迹原样——
                                      !   通道分类等片段间距消费方读这个）
end type fin_t
```

## 输入输出

### `keys` —— 自定义输入关键字声明

形：function 返回 `character(len=16)` 数组；解析前声明。

```fortran
function keys() result(kk)
   character(len=16) :: kk(2)
   kk = [character(len=16) :: 'K_HARM', 'R_TERM']
end function keys
```

### `load` ——初始化（需初始化时写）

形：无参 subroutine，自持 `use input, only: pull`。
- 在这里填入keys中的自定义关键字，它们会在程序读取input时被识别；
- 读参数文件（如 NN
权重）、分配等初始化也在此（会先于一切力调用）

```fortran
subroutine load()
   use input, only: pull
   call pull('K_HARM', k_harm)              ! 现在在input中可以规定简谐振子的力常数
   call pull('R_TERM', r_term, r_found)     ! 

   ! 在这里可以额外进行读取参数、数组分配等等功能

end subroutine load
```





## 自定义非绝热方法（实验性）

只在做自定义非绝热动力学时写，仍在测试中；


### `prop` —— 电子传播成员（与 `match` 成对）

```fortran
subroutine prop(sta, dt)
   use state, only: state_t
   type(state_t), intent(inout) :: sta
   real(8), intent(in) :: dt
   ! ……电子态传播，具体算法由内置哈密顿量和方法决定……
end subroutine prop
```

### `match` —— 准经典统计匹配（与 `prop` 成对）

非绝热方法中使得经典参数和量子演化匹配的操作；某种广义的PES跳变；

```fortran
subroutine match(sta)
   use state, only: state_t
   type(state_t), intent(inout) :: sta
   ! ……跳变判断，操作等等……
end subroutine match
```

### `sync` —— 容器同步（可选随对）

把方法自持的电子态发布到外部交付容器 `s`（即 `sta%s`）：方法可以把
电子态留在包内私有模块里演化（如 cau TDHF 的密度矩阵），接口在每步
`prop` 与 `match` 之后各调一次 `sync`，把声明交付的数据（`s%a`/`s%rho`）
更新到位

```fortran
subroutine sync(s)
   use state, only: elec_t
   type(elec_t), intent(inout) :: s   ! 发布: 包内电子态 -> s%a / s%rho
   ! ……把方法自持的电子态拷入交付容器……
end subroutine sync
```

### `init` —— 方法初始化

非绝热方法本身的数组分配/参数读取等等；

```fortran
subroutine init()
   ! ……分配/参数落位……
end subroutine init
```

### `method_name` —— 方法注册名

返回 `character(len=32)` **小写**名；在`ELEC_METHOD = <名>` 中选用。

```fortran
function method_name() result(s)
   character(len=32) :: s
   s = 'my_tdhf'
end function method_name
```
### 数据结构参考

(程序对电子状态的规定，可能用到)

```fortran
type :: elec_t                        ! 电子态容器（sta 的 s 分量）
   integer :: n_surf                  ! 势能面数——单体表象下读作组态空间基数，
                                      !   语义由电子方法成员定义（字段名不变）
   complex(8), allocatable :: a(:,:)  ! 振幅系数矩阵（传播容器）
   complex(8), allocatable :: rho(:,:)! 密度矩阵——rho 本征方法的本体；a 本征方
                                      !   法的派生缓存（由 sync 成员刷新）
end type elec_t
```

