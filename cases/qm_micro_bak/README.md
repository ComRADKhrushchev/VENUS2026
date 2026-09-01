# qm_micro — 量子微正则简正模采样

> 运行：`cd cases/qm_micro && ../../venus_test.e`（主输入文件 `input_qct.txt`）

## 1. 目的

本 case 验证量子微正则简正模采样（`INIT_SAMPLING_A=QM-MICRO`，NACTA=8）在等边 LEPS H₃
上的正确性（NT=100，ENMT_A=8.38 kcal/mol）：量子数守恒、能量不可达模拒绝、简并对均匀取态
与能量闭合。

## 2. 理论与体系

体系为 TEST 势 LEPS 等边 H₃（NLINA=0；E′ 对 1013.5 cm⁻¹，呼吸模 1803.8 cm⁻¹）。
LEPS 参数单位为 eV/bohr/bohr⁻¹（`src_TEST/test_potentials.f90`）。机制链在 `QMMICRO`
（手册外扩充，前 fork 谱系继承）中实现：可用量子能量 DUM1=(ENMT_A−ZPE)·CAL2CM（cm⁻¹；
kcal/mol→a.u. 换算 C1=0.04184，`venus_params.f90:60-80`）；高频模以均匀候选经
Beyer-Swinehart 态密度拒绝（DENQ 比值裁决）；末 2 模在 ±100 cm⁻¹ 容差网格上枚举均匀取态。
本 case 的态选择：DUM1=1015.7 cm⁻¹ 时呼吸模（1803.8 cm⁻¹）不可达——占据将留负剩余能量，
DENQ 拒绝（物理正确）；E′ 简并对恰接纳 (1,0)/(0,1) 两态均匀。全链跑 RAND0(ISEED)，
无系统种子依赖，双跑逐位确定。

登记隐患（本 case 未触发）：若 DUM1 使网格枚举 SCOUNT=0，`QMMICRO.f:94-106` 留
ANQ(1)/ANQ(2) 未赋值——其他 ENMT_A 值下的潜在隐患，保留不改。

## 3. 方法与流程

1. 程序读取输入并回显 QUANTUM MICROCANONICAL NORMAL MODE SAMPLING 与 HSCALE=8.380
   （stdout，results.txt）。
2. 初始化随机数（ISEED=20260830），对 NT=100 条轨迹由 `QMMICRO` 逐模抽取量子数
   （每轨迹一行 ANQ 回显，共 100 行，results.txt）。
3. `INITQP` 按量子数组装坐标动量，能量闭合缩放（stdout `INTERNAL ENERGY` 回显）。
4. 每轨迹传播 NS=3 步并打印相空间（fort.1001 起每文件一条轨迹）。

## 4. 核心验证

**判据**：呼吸模不可达时应全部被拒绝（Σ(n3)=0），简并对取态应在两态间均匀
（占比 ∈ [0.3, 0.7]），总能量闭合应远小于网格容差 ±100 cm⁻¹。**实测**（results.txt，
100 行 ANQ 回显）：Σ(n1+n2)=100=NT、Σ(n3)=0；(1,0,0) 占比 57/100（理论 0.5）；
max|E_vib−ENMT_A|=2.07 cm⁻¹（E_vib 2928.9 vs 2930.9 cm⁻¹）。**结论**：态密度拒绝与
简并对均匀取态机制按设计工作。佐证：第二能量点（ENMT_A=11.19 kcal/mol）网格枚举三态
均匀（74/62/64，χ²=1.24，p=0.54）；双跑逐位一致（106/107 fort 文件一致，stdout 仅
耗时行差异）。

> 图：`fig_sampling_energy.png` —— 10 条入库轨迹的采样能量（stdout `INTERNAL ENERGY`
> 回显，活数据）：左为各轨迹初始内能（8.376-8.386 kcal/mol，红线为目标 ENMT_A=8.38——
> E′ 简并对取态 (1,0)/(0,1) 的能量简并），右为逐轨迹缩放收敛过程（如轨迹 1：9.56→8.30→8.39）。
