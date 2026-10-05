# 2026-09-30 cau_tab 容器立项：DFT 数据点表 + 1D Akima 插值驱动的散射敏感性实验容器

## 动因

用户要测"某个局部特征对总体散射结果的影响"：把 DFT 数据点本身变成实验旋钮——手工
调整表数据的局部特征，重跑散射，对比电子通道转换、能量分配、漂移口径。为此另立
功能单一的新容器 `cau_tab`：PES 唯一来源是表文件（六片 × 三位点 × z 数据点 +
每位点 1D Akima 插值），无解析式、无 NN、无中间拟合层。cau_peri 不动（基线容器）。

## 裁决记录（用户裁决 2026-09-30，六条）

1. **另立新容器**：不改 cau_peri；新容器表驱动唯一路径，居 `systems/inactive/cau_tab/`
   架上。
2. **数据谱系**：用户另行提供新 DFT 表；本批交付格式规范（容器文件夹
   `TABLE_FORMAT.md`）+ 装载器 + 判据。格式以旧树 rst 六态文件为样本。
3. **横向处理**：三位点（top/bridge/hollow）+ 三系数傅里叶，与 cau_peri 的 NN 路径
   同构；hollow 取 FCC/HCP 均值还是单取，由数据制备定并在表头注明。
4. **插值元**：Akima（C1 局部，稀疏点无过冲，无需求解线代）。
5. **TDHF 副本随容器携带**（家族先例：四 cau 容器各持引擎拷贝）：电荷转移/自旋通道
   可观，与 cau_peri 基线同引擎同键同缺省——A/B 散射对比中唯一差异是 PES 数据源。
   `ELEC_METHOD = adiabatic` 参考面照留（绝热退化为纯 GS 片基态）。
6. **无声明式修改通道**：直改表副本。表文件即实验旋钮（`TABLE_FILE` 键指向运行目录
   的手编副本）；装载器自检 + 装载回显承担手编表防错与出处记录。

## 容器构成

| 文件 | 内容 |
|---|---|
| `cau_tab_tab.f90` | 表存储/装载/自检 + Akima 建斜率与求值（二分定段 + 三次 Hermite，值+导数）。越界：低于首点 fatal 点名 sheet/site/z；高于末点末值保持+零梯度。装载回显进 stdout |
| `cau_tab_pot.f90` | 物理层：自持三系数傅里叶横向对（周期取表头 LATTICE）；`state_val` 缝（每片三位点 tab_eval → 傅里叶合成，z 梯度由位点斜率同一傅里叶合成）；warp_z 冻结线、Vcp 远场衰减、hub_chain、ground/hubbard 交付、核装配自 cau_peri_pot 移植（代数不变）；晶场劈裂解析自带（2026-09-28 kcal→eV 单位裁决随行） |
| `cau_tab_tdhf.f90` | cau_peri_tdhf.f90 机械改名拷贝（含 2026-09-28/29 两轮引擎修复的现形态）；与 cau_peri 引擎同步义务记于头注 |
| `cau_tab_box.f90` | 12 声明键（cau_peri 的 13 键去 USE_NN/NN_FILE 加 TABLE_FILE）；roster↔engine 置换包装与终止判断照抄 |
| `reg_cau_tab.f90` | `reg_cau_tab_slots()` 单一交接，四布线组同构 reg_cau_peri |
| `system/` | C.xyz + Slab.poscar（自 cau_peri 复制；反应物定义居容器文件夹） |
| `TABLE_FORMAT.md` | 表格式规范 + 手编指南 + 越界政策 |

## 物理形态说明（与 cau_peri 的差异，明记）

- **位点处两种横向构造重合**：cau_peri 解析引擎做"位点参数傅里叶插值 + z 向解析式"，
  cau_tab 做"位点能量 Akima + 能量傅里叶插值"。在三个位点上两者给出同一函数（傅里叶
  在位点精确再现位点值）；**离开位点，二者一般不同**（参数级 vs 能量级插值的固有差
  异）。这是表容器的物理建模选择，不是缺陷；判据 (t3) 在位点上对拍、(t3b) 报告离点
  差异尺度。
- 判据程序可用解析引擎在三位点采样生成"恒等表"，验证全链装配（装载器→傅里叶→
  hub_chain→交付）与 cau_peri 在位点几何一致。

## 判据（scripts/check_cau_tab.f90）

- (t1) Akima 数学腿：已知光滑函数建表，格点中点求值与导数 vs 精确值。
- (t2) FD 核-值对拍腿：解析引擎生成的恒等表装载后，ground/hubbard 核 vs 中差
  （1e-5·max(1,|g|)，同 check_cau_peri (c)(d) 口径；表核是自身插值函数的精确导数，
  故阈值取 (d) 档而非 (c) 的引擎异常容差）。
- (t3a) 位点对拍（在格点上）：cau_tab vs cau_peri 的 ground/hubbard_val 在三个位点
  几何一致（插值函数过数据点 → 容差 1e-10）。
- (t3b) 离点对拍（格点中点）：同上，容差按插值误差定标（2e-3·max(1,|v|)），
  measured 上报。
- (t4) 装载腿：表结构自检 + tab_info_out 回读（周期/位点数）。
- 引擎腿移植自 check_cau_peri：(a)(b)(e)(f)(f2)(g)(i1)(i2)(i3)(h2)(o)(j)(k)(l)，
  阈值不变（TDHF 拷贝等价性 + 表模式下的引擎行为）。
- (n1) NN 装载腿不移植（本容器无 NN）。

## 工作流与限制

- **A/B 对比受一容器一进程 statute 约束**：cau_peri 与 cau_tab 不能同时居 `systems/`
  顶层；对比 = 换文件夹 + make 两次，input_qct.txt 逐键相同。判据程序侧不受此限
  （独立编译两容器，模块名不冲突）。
- 稀疏网格（每态 10~32 点量级）对局部特征的形状控制有限；需受控光滑扰动时在表里
  加行加点即可。
- 晶场劈裂 δcf 保持解析（无 DFT 数据），与 cau_peri 同源同值。

## 验证腿

1. check_cau_tab gfortran 编译运行全绿（新腿 + 移植腿）。
2. 容器源 gfortran + ifx 双编译零错（家族四容器先例）。
3. check_cau_peri 复跑基线不动（cau_peri 零触碰的证据）。
4. 主构建树零影响（inactive 架，glob 之外）；激活 = 移入 systems/ 顶层 + make 两次。

## 后续候选（待裁决）

- 离线平滑预处理工具（旧树 interp_data 的 Reinsch 平滑 → 均匀网格表）：若 Akima
  力折痕在长轨迹能量守恒上超预期，升格此路径。
- 用户 DFT 表到位后：首套真实数据表 + 冒烟散射批（对比 cau_peri 解析基线）。


---

## 追加裁决（用户 2026-10-01）：通道直读 + GS 退参考 + E0 基底

**动因**：四透热态片各自插值再经 hub_chain 相减（J = E₁D−E₃P 等），每片的插值误差
与 C1 折痕（抽稀后各片节点位置还不同）叠加进通道，可致 J/U/h 振荡。

**裁决三条**：
1. **通道直读**：表携带五通道 E0/h/U/J/VCP，各自 Akima+傅里叶插值直接交付；
   hub_chain 四态加减退役（`channel_eval` 替位）。梯度核为各通道自身插值的精确
   导数，无跨片算术。晶场劈裂 δcf 仍解析（h(1/2) = h ± δcf）。
2. **GS 数据退参考**：GS 的 DFT 数据不参与计算；力学基底 = E0 通道。
3. **数据源**：rst 表 E0 ← rst/E0_Data_GS_rst.txt（SCF 链 E0 导出；其文件头
   a=2.88 字样为导出工具陈旧模板，站点坐标是 a=2.95 自帧）；current 表通道 ←
   current/Data_{h,U,J,V}_aug.txt，E0 ← scf_nas/E0_Data_GS.txt（**同日续裁**：
   用户指认 SCF 链导出为 current 线 E0 源，maker 以 argv[5] 覆盖参数接入；
   数据性质：z < 1.25 Å 深点为 aug 外推噪声（hcp 振荡至 −16 eV，运行期下限
   以下不可达），覆盖带内无正排斥壁，压穿 1.25 Å 的轨迹 fatal）。

**落地**：cau_tab_tab 词表 6→5；cau_tab_pot 换 channel_eval；TDHF 拷贝零改动
（入口签名不变）；判据生成器改通道采样（E0 ← 解析基态片，h/U/J/V ← cau_peri
hubbard_val，h̃ = h₁h₂ 均值），150 PASS/0 FAIL gfortran+ifx 复验；两表重建，
恒等表实测通道直读收益：t3b |dJ| 4.4e-7（片加减时代 1.6e-5，降 36 倍）。


## 追加裁决 III（用户 2026-10-01）：E₀ 定案 + NN 差接受 + energy_feedback 弃用

1. **E₀ = SCF 链 E₀**（两表数据源现状即定案；NN 的 GS 拟合基态片为过时口径）。
2. NN 与 rst 通道数据的 h̃/U 系统差（~0.4–0.9 eV，J 一致）**接受为谱系差**，不重训。
3. **energy_feedback 不复刻**：历史程序硬编码开启的每步速度重标定总能量钳制
   （旧树 CLAUDE.md 在案），按现行裁决视为过时——venus2026 未携带即最终态。
   附记（同日实测）：无钳制下 NN 分支 TDHF 在真实生产条件同样发生平均场捕获
   （z 至 −1.5 Å、漂移 −88 kcal）；rst 表修复顶位伪影后正常散射（漂移
   ≤0.12 kcal）。裸动力学口径下各谱系的中场 h̃ 梯度决定捕获强度——这正是
   敏感性实验容器要测的量。
