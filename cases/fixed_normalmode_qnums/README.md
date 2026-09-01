# fixed_normalmode_qnums — 固定简正模量子数采样

> 运行：`cd cases/fixed_normalmode_qnums && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

验证固定简正模量子数采样（`INIT_SAMPLING_A=NORMAL-MODE`）在 LEPS 等边 H₃ 极小
（NLINA=0，QZA_EQ 边长 1.759910 Å）上的正确性（NT=1000）：固定 `ANQ_A=(1,2,3)`，
考察量子数打印、谐波能级能量闭合与逐轨迹能量集中度（能量固定、相位随机）。

## 2. 理论与体系

谱为 E′ 二重简并对 1013.5 cm⁻¹ 与 A₁′ 呼吸模 1803.8 cm⁻¹（stdout 频谱表）。
流程与微正则简正模采样相同（`NMODE` → `INITQP`），仅能量分配不同：每模取谐波
能级 (n_i+½)ℏω_i，谐波能级和

  E = [1.5×1013.47 + 2.5×1013.51 + 3.5×1803.83] / 349.755 = 29.64 kcal/mol

（349.755 为 cm⁻¹→kcal/mol 换算因子）。随后能量闭合缩放至完整非谐 LEPS 势；
随机相位由 `INITQP` 逐轨迹重抽——能量固定、相位随机，是固定量子数语义的关键。
关联缺陷：`NACTA=4`（LOCAL-MODE）的 Morse 参数在 fork 中无输入来源，W(0) 越界 →
NaN 与死循环，登记为缺陷 F20（不修）。

## 3. 方法与流程

1. 读取输入并打印 NORMAL MODE SAMPLING 与量子数表（stdout，results.txt）。
2. 初始化随机数（ISEED=20260826），NT=1000 条轨迹按 ANQ_A=(1,2,3) 取谐波能级分配
   能量，`INITQP` 抽取随机相位并组装坐标动量，能量闭合缩放。
3. 每轨迹传播 NS=10 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

见 `fig_sampling_stats.png`。

> 图注：左 = 各轨迹初始内能（29.61-29.66 kcal/mol，极窄展宽——固定量子数下能量
> 唯一、相位是唯一自由度），右 = 逐轨迹缩放收敛过程（stdout `INTERNAL ENERGY`，
> 10 条轨迹）。判定依据：ANQ 行每轨迹 1.00 2.00 3.00 与输入一致；EVIBA 均值
> 29.636 vs 理论对照值 29.642 kcal/mol（相对偏差 0.0010，门限 0.5%）；1000 条
> 轨迹 EVIBA 展宽 0.059 kcal/mol（< 0.1，相位为唯一自由度）。
