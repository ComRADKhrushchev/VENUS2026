# rigid_surface_defect — 刚性表面（NSURF=2）缺陷修复回归验证

> 关联：真实势下的表面散射/均衡现已由 `cases/emt_beam_scattering`（EMT-NN 束流散射）
> 与 `cases/emt_surface_md`（EMT-NN 300 K 恒温均衡）覆盖。

> 运行：`cd cases/rigid_surface_defect && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

本 case 为缺陷修复回归验证：NSURF=2（RIGID 表面）的四重缺陷 D1–D4 中 D1–D3 已修复，
本卡复跑确认损坏签名消失；D4 部分实现，剩余范围见下表。

## 2. 理论与体系

体系为 TEST 势 HARMONIC；A 单原子气体（入射高度 8.0 Å），B 4 原子表面，NT=20、
NS=5。四重缺陷及修复状态：

| 缺陷 | 内容 | 出处 | 状态 |
|---|---|---|---|
| D1 预设语义错误 | `MODEL=RIGID-SURFACE` 无 `SURFACE_MODEL` 时 map_old_nsurf(2)=1 → 实际 NSURF=1；本输入显式 `SURFACE_MODEL=RIGID` 绕开 | venus_input.f90 | **已修复**：旧值重映射按 NSURF 原义转发 |
| D2 B 坐标忽略 | 仅 NSURF=0 读 QZB_EQ → NSURF=2 下 B 原子 Q=0、W(B)=10³⁰ | venus_input.f90:774 / SELECT.f90:164 | **已修复**：NSURF=2 亦读入 B 表面坐标 |
| D3 胞定义越界读 | 读 Q(3·(NATOMS+1)..3·(NATOMS+4)) 越界 → SKEW=NaN → 全轨迹 NaN | SURF.f:30-50 | **已修复**：按 5×5 方格胞显式组装，越界读消除 |
| D4 五坐标采样未实现 | b+φ₁ 圆盘瞄准、θ sin-θ CDF、φ₂ 取向、s 分离未实现 | SURF.f | **部分实现**：斜胞入射链可用；五坐标截面语义采样仍缺 |

## 3. 方法与验证

1. 复跑同输入（ISEED=20260834，NT=20，NS=5），检查 fort.88 瞄准点与 fort.1001-1020
   相空间轨迹。
2. **修复前签名（历史）**：fort.88 全 20 行 RX0=NaN、RY0≈1e-155 垃圾值；fort.1001
   step 0 即 H=NaN。
3. **修复后（2026-09-01）**：fort.88 RX0/RY0 全部有限（首行 4.64e-2/0.750），
   fort.1001 无 NaN，H 守恒极差 0.0000，束流 A z=8.0 精确。
4. 见 `fig_long_trajectory.png`：修复后入射原子 z(t) 自 8.0 Å 起正常下降（fort.1001，
   NS=5 原始数据）——正常行为，替代旧 NaN 证据图。

> 附注：HARMONIC 测试势下方格表面无近距排斥，A 长时程可穿过表面平面（测试势固有
> 属性，非缺陷）；真实刚性表面散射应使用 EMT-NN 势（见 `cases/emt_beam_scattering`）。
