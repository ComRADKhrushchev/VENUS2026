# V 计划：全链验证跑（数值后验 + case 集建立）

> 依据：docs/final_completion_plan.md §2（V/S 行：new case set over the full
> chain；旧程序只作 ad-hoc 同量参考非门）；全局约束 8（判据先行）；C5-R1
> （B 层/链级数值后验——法定）；C2/C7/C8 后验清单（→V 各条）；GWRITE_LEVEL
> 用户裁决仍悬（随 V 收口再催一次）。
> 用户 2026-09-15 放行执行。
状态：已确认（2026-09-15）
日期：2026-09-15

**目标：** 完工主体放上**真力链**跑数值后验：C5 期机制先行的 B 层三员与
频率依赖成员现在有真力可验；建立新树 **case 集**（cases/ 惯例沿旧库：
README/input/expected）；守恒判据家族全链全跑。

**范围裁定（拟）：**
- S1 **case 集建立**（cases/ 新目录，沿旧库三件套惯例：README.md（验证
  目的/原理/预期与结果）+ input_qct.txt + system/ 文件夹 + results 记录）：
  - case_harmonic_scatter：分子两片段 + 谐振容器 + verlet——能量守恒/
    统计面基准（与 check_driver 判据呼应但以 case 形落档）；
  - case_radau_conservation：同体系 INTEGRATOR=RADAU——自适应守恒；
  - case_surface_window：表面体系（POSCAR + A_LAT/BOXPAD + 分子片段）
    ——表面窗口键集端到端（气相判据 gfortran 树承载）；
  - case_b_layer：B 层三员数值后验的 case 化（thermostat/thermalize/
    glo_target 各一小节，判据见 S2）；
  - case_chain_frequencies：Hessian→频谱→boltzmann/normalmode 真值
    （S3 的 case 化）。
- S2 **B 层数值后验**（C5-R1 法定）：真力链下——thermalize（Maxwell 种子
  + MD 平衡：平衡后动能趋近目标温度的统计带——容差与样本量写明）；
  thermostat（速度重置的平衡温度维持）；glo_target（FDT 装载的数值一致
  性：阻尼-噪声关系在链上的实现核对——解析后验而非长程统计）。
  - 450× 备忘（V.1 移交，task-1 报告）：radau 下 thermalize 平衡覆盖约
    450× verlet 的物理时长（自适应步长序列增长所致）——成员语义即"走
    注册传播器"，属行为差异非缺陷；B 层平衡统计以 verlet 为标定基
    （case 用 verlet，与标定一致），radau 复用需按其步长序列重标定。
- S3 **链级频率真值**（C5-R1 法定 + dens/eig 前置演示）：reference 容器
  真力 → Hessian → eig_sym → 频谱（谐振体系的解析频率对照——位级/
  容差写明）→ boltzmann/normalmode 以真频率采样（几何分布/微正则锚复验
  ——C5 判据的合成频率版换成真频率版）。
- S4 **守恒家族全跑 + V 报告**：全 case 能量/动量守恒判据汇总表；ifx
  旧链位级复验（可选——环境允许则做）；thta 退化文件用户说明（C8 移交）
  入 case README；GWRITE_LEVEL 裁决催办注记。

**设计裁决（拟）：**
- D1 判据先行：每 case 的预期锚先写进 README（解析/统计带），跑后
  results 对照——不绿即 BLOCKED 修主体（V 是验证不是修补场：发现的
  主体缺陷回相应模块的 fix 流程，记 V 发现清单）。
- D2 E12：case 材料是体系定义（用户侧）——体系知识进 case 文件夹
  合法（这正是 container 材料的家）。
- D3 R5 沿用（case README 是文档非源码——旧库引用亦不禁，但保持简洁）。

**步骤：**
- [ ] **V.1 case_harmonic_scatter + case_radau_conservation**：两基准 case
  建立 + 跑 + results 入档；
- [ ] **V.2 case_surface_window**：表面窗口端到端 case；
- [ ] **V.3 B 层数值后验**（case_b_layer）：三员真力链判据；
- [ ] **V.4 链级频率真值**（case_chain_frequencies）：S3 全链；
- [ ] **V.5 V 报告 + 守恒家族汇总**：S4；提交。

**验收：**
1. 五 case 建立，README 预期锚先行，results 全绿入档；
2. B 层三员与频率链数值后验落地（C5-R1 关闭）；
3. 守恒家族汇总表 + V 发现清单（主体缺陷若有——回修流程）；
4. 双编译器绿（case 跑 gfortran 树——ifx 目录门环境备忘沿档）。

**风险：**
- B 层统计判据的容差标定（样本量 vs 运行时长）——先小样本试跑再定
  锚，写明；
- 表面 case 的 POSCAR 样件与 A_LAT 一致性（几何自洽）。
