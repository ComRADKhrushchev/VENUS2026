# benchmark_oldvenus —— 干净的 VENUS05 血统标定树

从 `Codes/OLD VENUS VERSIONS/VENUS_TSH_HO2` 剥离 EANN/NN、TSH、AmberFF 与表面
(GLO/摩擦) 支路后重建的最小树，专用于将来的标定计算（首个目标：H+H2/BKMP2
转动激发截面，与本仓库新程序的 QCT 对比）。原树未动一字。

## 布局

```
benchmark_oldvenus/
├── Makefile          # gfortran；make → venus.e
├── SIZES             # 尺寸参数（原树 src_VENUS/SIZES 原样）
├── src/              # 72 个库存机械件（.f，原样复制，零改动）
│   └── VERLET.f      # 原树 0VERLET.f 换名（CZZ 的 VERLET.f90 已弃）
├── sys/              # 系统侧五缝 + 官方引擎
│   ├── PES_bkmp2.f90       # POTPRE/POT0/DPESHON 用户势能缝（BKMP2 占位）
│   ├── bkmp2_engine.f      # 官方 surface950621 引擎（本仓库 systems/bkmp2.f 副本）
│   ├── TEST_sys.f          # 通道判定（原 src_NN/TEST_3ATOMS.f，读 ./channel）
│   ├── GWRITE_sys.f        # 逐轨迹打印（原 src_NN/GWRITE_3ATOMS.f）
│   └── friction_stub.f90   # FRICTION/FRICFORCE/readff/swap 桩（NGLO=0 时不可达）
└── work/             # H+H2 冒烟运行目录（input.inp + channel + smoke_output.txt）
```

## 切割清单（相对原树）

- 删 `src_NN/` 全部（EANN 势链、TSH、friction、AmberFF/V_MFF/Zmatrix、readinput…）
- 删 TSH 死件 `POT_SH.f`、死 main `0VENUS.f`、CZZ 改版 `VERLET.f90`
- 删 FFpara、channel 可执行、work_dir 的 fort.* 与 .mod/.o 陈旧产物
- 保留 VENUS05-Manual.pdf 于原树（未随本树分发；输入语义详见手册）
- 五个活缝由 sys/ 重实现：POTPRE/POT0/DPESHON（势能）、TEST（通道）、GWRITE（打印）；
  FRICTION/FRICFORCE/readff/swap 为噪声桩

## 缝契约（sys/PES_bkmp2.f90 头注同文）

- `POT0(NATOMS, V)`：读 COMMON/QPDOT 的 Q [Å]，返回 V 内部单位（kcal/mol×C1）
- `DPESHON(NATOMS)`：填 `PDOT = −∂V/∂q` [内部单位/Å]（**物理力非梯度**，
  0VERLET.f 头注 `F(T)=−DV/DQ=PDOT`；PDOT 直接推进坐标）
- 内部单位与主程序一致：C1=0.04184（kcal/mol→内部）、C7=0.063508（ℏ）、Å/amu/(10 fs)

## 冒烟结果（work/，2026-09-28）

H+H2（EBK v=0,j=0，NREL=1 固定 8.07 kcal/mol，B_MAX=6，S=10 Å，verlet 0.1 fs）：

- 初始 E_tot = **−95.22 kcal/mol** = V(H2,re)=−109.5 + ZPE 6.2 + E_rel 8.07 ✓
  （BKMP2 引擎 + 单位换算 + EBK 初态 + 束流装配全链正确）
- 双轨迹终止（~130 fs）、分类路径 1、FINLNJ 末态 **N′=0.34/0.01，J′=0.00**（连续值，
  j(j+1) 反演约定，与本仓库 convA 同式）
- verlet 能量振荡 ±3 kcal/mol（0.1 fs 步长下 H2 零点振动采样所致）；
  **INTEGRATOR=4 (RADAU) 守恒达 3×10⁻⁸ kcal/mol** 但 TIME 参数语义不同（作为精度/
  种子），10000 循环内不终止——生产用需查手册调参

## 已知问题（不阻塞使用）

1. 汇总计数器全零（"NON-REACTIVE/REACTIVE/EFFECTIVE"）——变体的主循环计数分支
   对本配置失效；逐轨迹 GFINAL 输出完整正确，截面统计按轨迹行后处理即可
2. FINLNJ 对远距"产物对"打印 59 条 WARNING（中间态分析噪声，库存行为）
3. work/ 每次运行产生 MFFinput/fort.*（机器 OPEN 语句残留，无害）

## 输入行映射（work/input.inp 逐行注释；比原树模板多三行）

原树 work_dir/input.inp 模板**滞后于代码**：缺 num_res 行与 CZZ 2024/3/8 的
NNAtomnum/NNAtoms/NFFAtoms 三行（零项列表读也会吞一行——三行都要给）。本树
input.inp 为核对过的行序；NREL 语义与模板注释相反：**NREL=1 才是固定平动能
(kcal/mol)，NREL=0 是温度采样**。channel 文件第 I 行原子对必须与路径 I 的定义
一致（路径 1=非反应对 (2,3)，2=(1,2)，3=(1,3)）。

## 下一步（σ(0→2) 基准）

同 E、同 b² 均匀、同 B_MAX=6 下跑 10⁴ 轨迹/能点，从逐轨迹 "PRODUCT IS A DIATOM:
N=/J=" 行取 j′ 分箱（round 即 convA），σ = πB_MAX²·N/N_tot，与本仓库 rotex3 对比。
