# VENUS 现代化改进对照

以旧版 VENUS 为参照，展示本仓库的现代化改进结果；

## 1. COMMON 块 → 模块

**旧版**：主程序声明约 48 个 COMMON 块，靠 INCLUDE 'SIZES' 的编译期
定长数组支撑，任何子例程都可静默修改全局状态：

```fortran
      PROGRAM VENUS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA+4)
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      ...（共约 48 个）
```

**新版**：全部状态收进两个模块——`venus_params`（编译期常数）与
`venus_data`（轨迹状态，数组改为运行时 ALLOCATABLE，按实际原子数分配）：

```fortran
program VENUS
  use input_parser, only: get_real
  use venus_params
  use venus_data
  use venus_input, only: read_venus_input, NMA, ENU, ...
```

```fortran
! venus_data.f90 — 全局轨迹状态（替代全部 COMMON 块）
module venus_data
  use venus_params
  implicit none
  public
  integer :: natoms = 0          ! 运行时维度（原 NDA 编译期定长）
  integer :: i3n   = 0
  real(8) :: T, V, H, TIME       ! 原 COMMON/PRLIST/
  integer :: NTZ, NT
  real(8), allocatable :: Q(:), P(:), PDOT(:), QDOT(:)   ! 原 COMMON/QPDOT//PQDOT/
  real(8), allocatable :: W(:)
  ...
end module
```

收益：显式依赖（`use` 可追溯）、数组按需分配（无 NDA 越界风险）、
编译器可检查接口。

## 2. GOTO 网络 → 结构化命名循环

**旧版**：主循环用数字标签 + GOTO 织成网（0VENUS.f 计 35 处 GOTO）：

```fortran
      IF(NCHKP.EQ.0)GOTO 223
      ...
      IF (NC.EQ.NS.AND.NCHKP.EQ.1) THEN
         NTZ=NTZ-1
         GOTO 451
      ENDIF
      IF (NSELT.EQ.0.AND.NCHKP.EQ.-1) GOTO 425
      IF (NSELT.EQ.1) GOTO 447
      IF (NSELT.EQ.0) GOTO 400
      GOTO 400
  223 IF (NSELT.EQ.2.OR.NSELT.EQ.3) CALL RANDST(ISEED)
  451 NTZ=NTZ+1
      ...
  400 NC=NC+1
      ...
            GOTO 451          ! 下一轨迹
      ...
            GOTO 400          ! 下一积分步
```

**新版**：双层命名循环 + 布尔标志，控制流一眼可读（原 GOTO 位置以
注释留痕便于对照）：

```fortran
  TRAJECTORY_LOOP: DO
     NTZ = NTZ + 1
     ! Original 451: check for completion
     IF (NTZ > NT) THEN
        ...（全局统计输出）
        STOP
     ENDIF
     ...
     trajectory_done = .false.
     INTEGRATION_LOOP: DO WHILE (.NOT. trajectory_done)
        ! Original 400: NC=NC+1
        NC = NC + 1
        IF (INTEGRATOR.EQ.1) THEN
           ...
           ! Original GOTO 451: next trajectory
           trajectory_done = .true.
           CYCLE INTEGRATION_LOOP
        ELSEIF (INTEGRATOR.EQ.2) THEN
           CALL SYMPLE(LLL,NSELECT)
           ...
        ENDIF
     END DO INTEGRATION_LOOP
  END DO TRAJECTORY_LOOP
```

## 3. 顺序 READ(5,\*) 输入 → 关键字解析

**旧版**：输入按行序强耦合——增删一个参数全链重排，且靠魔法数字
（NSELT=2/NACTA=5/NSURF=1）表达语义：

```fortran
      READ(5,*)NSELT,NSURF,NTHTA
      IF (NSURF.NE.3) NTHTA=-1
      ...
      READ(5,*)NATOMS,NFC,NGLO
      ...
      READ(5,*)NOB,BMAX
      ...
      READ(5,*)NPATHS
      DO I=2,M
         READ(5,*)RMAX(I),RBAR(I),NATOMA(I),NATOMB(I),DELH(I)
```

**新版**：`KEYWORD = VALUE` 逐行解析（input_parser + param_mapping +
presets），语义化关键字带默认值与合法性检查：

```fortran
    if (has_keyword('TASK')) then
       call map_task(trim(get_str('TASK','TRAJECTORY')), NSELT, ios, errmsg)
       if (ios /= 0) then
          write(6,*) trim(errmsg)     ! loud-stop：显式报错而非静默错读
          stop
       end if
    end if
```

输入文件从"必须背行序"变为自解释的（节选自 cases/rst_beam_scattering）：

```text
TEST_PES = RST
MODEL = FULL-SURFACE
EREL = 14.528, 6.5      # 入射能量 [kcal/mol]，高度 [Å]（SURF.f ACZ=S）
ISEED = 20260821        # 随机种子
NT = 20                 # 轨迹数
```

## 4. 固定格式 .f → 自由格式 .f90

**旧版**（72 列限制、列 6 续行符、C 列注释）：

```fortran
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
C         SELECT INITIAL CONDITIONS
C
      CALL SELECT
```

**新版**（自由格式，主程序与重构过的核心子程序均为 .f90）：

```fortran
  real(8), allocatable :: QZA(:,:), QZB(:,:)   ! 按片段×坐标运行时分配
  ...
  CALL SELECT
```

## 5. 移除的死代码与损坏分支

- ADAMSM 积分器（名义 6 阶实测 1 阶，能量漂移最差）及其启动链
  RUNGEK/PARTI；NMA.f（929 行零调用）
- MAINLZ 调度死岛（GPATH/MPATH/MPATHO/STATPT，非动力学功能）与
  NSELT=1 分支
- 半开入口（LOCAL-MODE/CI-QM-MICRO 关键字曾被接受但实现缺失）改为
  loud-stop 显式拒绝

## 6. 行为保持的验证方式

所有删除/重构以"删前删后轨迹逐位一致"为验收（同一输入 diff 全部
输出文件），即上述改动不改变任何数值行为——现代化是纯结构性的。

---

*对照基准：`OLD VENUS VERSIONS/VENUS_TSH_HO2/src_VENUS/0VENUS.f`
（2135 行）。本文各代码块为两侧实际摘录，非示意。*
