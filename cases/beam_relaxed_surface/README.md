# beam_relaxed_surface — 松弛表面束流入射几何组装

> 运行：`cd cases/beam_relaxed_surface && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

本 case 验证 RELAXED-SURFACE（NSURF=1）束流几何组装链：固定入射角 θ、固定方位角 χ、
瞄准点在斜胞内均匀抽取、隐含固定碰撞参数 b=S·sinθ、固定相对能与表面静止初始条件，
全部一次运行中闭合。

## 2. 理论与体系

体系为 H₂(A) + 3 原子表面(B)，共 5 个 H 原子，TEST 势 HARMONIC（NT=200）。A 为刚性 H₂，
键长 0.7414 Å；B 由 `QZB_EQ` 构成 5×5 矩形胞，SKEW=90°。关键几何参数：

| 参数 | 值 | 出处 |
|---|---|---|
| 入射角 θ | 30.0°（固定） | input_qct.txt / SURF.f |
| 方位角 χ | 0.0°（固定） | input_qct.txt / SURF.f |
| 入射高度 S | 8.0 Å | input_qct.txt |
| 隐含碰撞参数 b | S·sinθ = 4.0 Å | SURF.f 几何构造 |
| 相对平动能 SEREL | 5.0 kcal/mol | input_qct.txt |

几何机制在 `SURF.f` 中实现：斜胞由 BOXLX/BOXLY/SKEW 构造，`NRNDXY=1` 时瞄准点在斜胞内均匀抽取；
A 质心偏移 (S·sinθ·cosχ, S·sinθ·sinχ, S)，速度方向 U=(−sinθ·cosχ, −sinθ·sinχ, −|cosθ|) 恒朝表面；
`NREL=1` 时 VELA=√(2·SEREL/WTA)、VELB=0；`NZDOWN=0` 表示 A 接束流、B 静止。
其中 kcal/mol→a.u. 换算用 C1=0.04184（`venus_params.f90:60-80`，即 1 kcal/mol = 0.04184 a.u.）。

## 3. 方法与流程

1. 程序读取 input_qct.txt 并回显束流几何参数（stdout，见 results.txt）。
2. 初始化随机数（ISEED=20260833），对 NT=200 条轨迹循环抽取瞄准点（RX0、RY0，回显于 results.txt）。
3. 组装 A 片段初始坐标与动量（质心偏移与速度方向按第 2 节公式），B 片段置静止（fort.1001-1200 逐步打印）。
4. 每条轨迹运行短轨迹传播（NS=5 步），打印相空间（fort.1001 起每文件一条轨迹）。
5. 末尾做正则模分析诊断（fort.26）。

## 4. 核心验证

**判据**：隐含碰撞参数 b 应精确等于 S·sinθ=4.0 Å，且 A 高度 z=8.0、速度方向与
(−sin30°, 0, −cos30°) 的偏差应在数值噪声量级；瞄准点应通过均匀性 KS 检验。
**实测**：results.txt 与 fort.1001 记录 z=8.00000（step 0，偏差 0）、θ=29.9993°、
b=4.0（偏差 4.4×10⁻¹⁶）、速度方向偏差 ≤1.1×10⁻⁵；RX0/BOXLX 的 KS D=0.0369（p=0.94）、
RY0/BOXLY 的 KS D=0.0538（p=0.59）。**结论**：束流几何链按设计工作。
佐证：相对能 5.0 kcal/mol 偏差 0、B 动量 <10⁻³，且全链经 RAND0(ISEED) 双跑逐位确定。
已知 HARMONIC 势伪迹：B 片段 Hessian 退化，fort.26 正则模诊断表含 NaN，不传播到轨迹物理。

> 图：`fig_trajectory.png` —— 10 条入库轨迹的相空间（fort.1001-1010 活数据）：
> 左为原子 1 到最近原子距离 r_min(t)（NS=5 短轨迹内为 A 片段键伴距离，0.7414 Å 恒定），
> 右为原子 1 动量模 |P|(t)。
