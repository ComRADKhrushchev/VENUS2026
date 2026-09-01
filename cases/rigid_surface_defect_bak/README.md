# rigid_surface_defect — 刚性表面（NSURF=2）缺陷登记

> 运行：`cd cases/rigid_surface_defect && ../../venus_test.e`（主输入文件 `input_qct.txt`）
> 注：越界读非确定 → 双跑不逐位（差异本身即缺陷证据）。

## 1. 目的

本 case 为缺陷登记：锁定 NSURF=2（RIGID 表面）在当前 fork 下的实际损坏行为——四重缺陷
D1–D4。任何未来的 NSURF=2 修复都必须让本卡的损坏签名消失，再按刚性表面五坐标截面语义重建。

## 2. 理论与体系

体系为 TEST 势 HARMONIC；A 单原子气体（入射高度 8.0 Å），B 4 原子表面，NT=20、NS=5。
四重缺陷（源码出处如下）：

| 缺陷 | 内容 | 出处 |
|---|---|---|
| D1 预设语义错误 | `MODEL=RIGID-SURFACE` 写 NSURF=2，但无 `SURFACE_MODEL` 时旧值重映射 map_old_nsurf(2)=1 → 实际 NSURF=1；本输入显式 `SURFACE_MODEL=RIGID` 绕开 D1 | venus_input.f90 |
| D2 B 坐标忽略 | 仅 NSURF=0 时读 QZB_EQ → NSURF=2 下 B 原子 Q=0；SELECT 跳过 B 采样；W(B)=10³⁰（刚性无限质量） | venus_input.f90:774 / SELECT.f90:164 |
| D3 胞定义越界读 | 读 Q(3·(NATOMS+1)..3·(NATOMS+4)) 越出分配 → 零填充邻堆 → EDGE 全零 → SKEW=NaN → RX0=NaN → 整条轨迹 NaN | SURF.f:30-50 |
| D4 五坐标采样未实现 | 刚性表面应有的 b+φ₁ 圆盘瞄准、θ sin-θ CDF、φ₂ 取向、s 分离均未实现；NSURF=2 与 NSURF=1 共享斜胞入射链 | SURF.f |

## 3. 方法与流程

1. 程序读取输入并打印横幅 'NSURF = 2' + 'RIGID SURFACE'（stdout，results.txt）。
2. 初始化随机数（ISEED=20260834），对 NT=20 条轨迹走 NSURF=2 斜胞入射链（fort.88 记录
   瞄准点 RX0/RY0）。
3. 组装 A 束流初始条件，B 原子坐标动量置零（D2 零化）。
4. 每轨迹传播 NS=5 步并打印相空间（fort.1001-1020，全部含 NaN）。

## 4. 核心验证

**判据（损坏签名，非正确性）**：D3 越界读应使 SKEW=NaN → 瞄准点 RX0=NaN 传播到全部轨迹。
**实测**（fort.88 全部 20 行）：RX0='NaN'、RY0≈1e-155 垃圾值（非用户 5×5 胞）；fort.1001
step 0 即 H(eV)=NaN、Q(C) x 分量 NaN，全部 20 条 fort.1NNN 含 NaN；stdout 计 240 处 NaN。
**结论**：NaN 传播链与 D3 越界读签名被锁定。佐证：B 原子 q=p=(0,0,0)（fort.1001 逐行
Au 全零，D2）；与 NSURF 无关的束流部分仍精确（A z=8.0、速度方向 (−sin30°, 0, −cos30°)）；
双跑差异仅限 RY0 垃圾值与 vy ±0.0 符号噪声。本 case 是 F24（RIGID → NaN）的登记锚。
