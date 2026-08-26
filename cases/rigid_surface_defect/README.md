# 刚性表面缺陷登记（NSURF=2）— 当前损坏行为锁定

> 运行：`cd cases/rigid_surface_defect && ../../venus_test.e`
> 注：越界读非确定 → 双跑不逐位（差异本身即缺陷证据）。

## 验证目的

**缺陷登记 case**：锁定 NSURF=2（RIGID 表面）在当前 fork 下的实际（损坏）行为——
四重缺陷 D1–D4。任何未来的 NSURF=2 修复都必须让本卡断言**失败**，再按刚性表面
五坐标截面语义重建。

## 原理

体系：TEST 势 HARMONIC；A 单原子气体（高度 8.0 Å），B 4 原子表面。

- **D1 预设语义错误**：`MODEL=RIGID-SURFACE` 写 NSURF=2，但无 `SURFACE_MODEL`
  关键字时旧值重映射 map_old_nsurf(2)=1 → 实际 NSURF=1；正确触发须显式
  `SURFACE_MODEL=RIGID`（本输入以此绕开 D1，锁定 NSURF=2 真行为）。
- **D2 B 坐标忽略**：`venus_input.f90:774` 仅在 NSURF=0 时读 QZB_EQ → NSURF=2
  下 B 原子 Q=0；`SELECT.f90:164` 又跳过 B 采样；W(B)=10³⁰（刚性无限质量）。
- **D3 胞定义越界读**：`SURF.f:30-50` 读取 Q(3·(NATOMS+1)..3·(NATOMS+4)) 越出
  分配 → 零填充邻堆 → EDGE 全零 → SKEW=NaN → RX0=NaN → 整条轨迹 NaN。
- **D4 五坐标采样未实现**：刚性表面应有的 b+φ₁ 圆盘瞄准、θ sin-θ CDF、φ₂ 取向、
  s 分离均未实现；NSURF=2 与 NSURF=1 共享斜胞入射链。

## 预期与结果（损坏签名，非正确性）

- **NaN 传播链**：fort.88 RX0='NaN'（×20）、RY0≈10⁻¹⁵⁵ 垃圾值；全部 20 条
  fort.1NNN 轨迹文件含 NaN ✓
- **B 原子零化**：B q=p=(0,0,0)（QZB_EQ 被忽略 + B 采样被跳过）✓
- **缺陷局域性**：与 NSURF 无关的束流部分仍精确——A z=8.0、速度方向
  (−sin30°, 0, −cos30°) ✓

补注：双跑差异仅限 RY0 垃圾值与 vy 的 ±0.0 符号噪声——越界读非确定 + NaN
污染，本身即 D3 证据。本 case 是 **F24**（RIGID → NaN）的登记锚。
