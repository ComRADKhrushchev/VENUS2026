# qm_micro — 量子微正则简正模采样

> 运行：`cd cases/qm_micro && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

验证量子微正则简正模采样（`INIT_SAMPLING_A=QM-MICRO`，NACTA=8）在等边 LEPS H₃
上的正确性（NT=100，ENMT_A=8.38 kcal/mol）：量子数守恒、能量不可达模拒绝、
简并对均匀取态与能量闭合。

## 2. 理论与体系

体系为 TEST 势 LEPS 等边 H₃（NLINA=0；E′ 对 1013.5 cm⁻¹，呼吸模 1803.8 cm⁻¹）。
机制在 `QMMICRO`（手册外扩充，前 fork 谱系继承）实现：可用量子能量
DUM1=(ENMT_A−ZPE)·CAL2CM（cm⁻¹）；高频模以均匀候选经 Beyer-Swinehart 态密度拒绝
（DENQ 比值裁决）；末 2 模在 ±100 cm⁻¹ 容差网格上枚举均匀取态。本 case 态选择：
DUM1=1015.7 cm⁻¹ 时呼吸模不可达（占据将留负剩余能量，DENQ 拒绝，物理正确）；
E′ 简并对恰接纳 (1,0)/(0,1) 两态均匀。

**登记隐患（本 case 未触发）**：若 DUM1 使网格枚举 SCOUNT=0，`QMMICRO.f:94-106`
留 ANQ(1)/ANQ(2) 未赋值——其他 ENMT_A 值下的潜在隐患，保留不改。

## 3. 方法与流程

1. 读取输入并打印 QUANTUM MICROCANONICAL NORMAL MODE SAMPLING 与 HSCALE=8.380
   （stdout，results.txt）。
2. 初始化随机数（ISEED=20260830），NT=100 条轨迹由 `QMMICRO` 逐模抽取量子数（每
   轨迹一行 ANQ，共 100 行）；`INITQP` 组装坐标动量并能量闭合缩放。
3. 每轨迹传播 NS=3 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

见 `fig_sampling_stats.png`。

> 图注：左 = 各轨迹初始内能（8.376-8.386 kcal/mol，红线 = 目标 ENMT_A=8.38——
> E′ 简并对取态 (1,0)/(0,1) 的能量简并），右 = 逐轨迹缩放收敛过程（如轨迹 1：
> 9.56→8.30→8.39；stdout `INTERNAL ENERGY`，10 条轨迹）。判定依据：Σ(n3)=0
>（呼吸模全部被拒）、Σ(n1+n2)=100=NT；(1,0,0) 占比 57/100（理论对照值 0.5）；
> max|E_vib−ENMT_A|=2.07 cm⁻¹（远小于网格容差 ±100 cm⁻¹）。
