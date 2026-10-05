# cau_peri 容器能量守恒解析证书（conservation certificate）

日期：2026-10-04
对象：`systems/inactive/cau_peri/`（cau_peri_tdhf.f90 @ 09-29 14:21、cau_peri_pot.f90 @ 09-29 19:28、
cau_peri_box.f90 @ 09-28 11:33、reg_cau_peri.f90）+ 驱动链（src/driver.f90、src/dynamics/propagator.f90、
src/output/recorder.f90）。
方法：手工解析审计——把"由 F=ma 推出机械能守恒"的功-能定理推导在 TDHF 容器上逐项重做，
核查四个解析条件（C1 厄米性 / C2 变分一致性 / C3 力一致性 / C4 无未记账显式时间依赖），
对每个违反给出解析漂移率公式与可观测特征。**不涉及任何数值传播。**

---

## 0. 符号与监测量

- ρ ≡ rhom(i,k,is)：二自旋单粒子密度矩阵（diabatic 表象），n_siz = n_orb + 2。
- h ≡ hcore：实对称 one-body 矩阵（杂质对角 h_vec(1/2) = h̃ ± cf，bath 对角 e_metal（q 无关），
  耦合块 sqrt(w/2)·Vg，双 Fermi 分支等耦）。
- D ≡ diag(vsm)：Fock 的杂质对角位移，vsm(:,1) = W·n(:,2) + (W−J)·n(:,1)，自旋 2 镜像；
  W = U·ones(2,2)，J = [[U, j],[j, U]]。
- 监测能量（cau_peri_pes 第 4 步 + recorder 的 e_tot）：
  **V(q, ρ) = E0(q) + Σ_is Tr(ρ_is h_is(q)) + E_mf(n) − E_metal_zero**，
  E_mf = n₁ᵀWn₂ + ½n₁ᵀ(W−J)n₁ + ½n₂ᵀ(W−J)n₂，n = ρ 的杂质对角（occ_vec）。
- 主恒等式：若传播生成元 F̃ 满足 iℏρ̇ = [F̃(ρ), ρ] 且 F̃ = δV_elec/δρ，则
  dV_elec/dt = (1/iℏ)·Tr([F, F̃]·ρ) —— F̃ = F 时恒为零；
  若不等，**解析漂移率 = (1/iℏ)·Tr([F, F̃]ρ)**，只含瞬时量，单快照可算。

## 1. 逐条件审计结论（PERTURB_MODE = 0 时）

### C1 厄米性 — 通过
- hcore 构造对称（assemble_hcore 对每个耦合块填入镜像元素）；Fock = hcore + 实对角位移，实对称。
- evolve_rho 是严格酉共轭（ρ → C e^{iεt}(C†ρC)e^{−iεt}C†，行乘 ph、列乘 conjg(ph) 成对），
  Hermiticity 与迹在精确算术下守恒（rho_herm_out 是该条件的运行时探针）。

### C2 变分一致性（F̃ = δV_elec/δρ）— 通过（2026-09-28 修复后）
- vsm = ∂E_mf/∂n 逐项核对：∂/∂n₁[W·n₂ + ½(W−J)n₁ + ½(W−J)n₂] = Wn₂ + (W−J)n₁ ✓（自旋 2 镜像 ✓）。
  每自旋向量形式正确——旧的二标量交叉接线正是 C2 违反（F ≠ dE_mf/dn），已修复。
- 监测能量用 Tr(ρh)（one-body）+ E_mf（不是 Tr(F̃ρ)）——正确的变分记账形式，无双计数。
- occ_vec 的 |Re ρ_ii|：初值 PSD（FD 占据 ∈ [0,1]）+ 酉演化保谱 ⟹ ρ(t) 恒 PSD ⟹ 对角 ≥ 0，
  abs() 在整条轨迹上恒为恒等映射，解析无害（不产生 kink）。
- Gauss-Legendre 权重 w 只以 sqrt(w/2) 进入 h 的耦合块，q 无关，不参与变分导数。✓

### C3 力一致性（g = −∂V/∂q，Hellmann–Feynman）— 通过（2026-09-28/29 修复后）
g = g0 + Σ_is [Σ_ii ρ_ii·∂h_ii/∂q + Σ_i Σ_ii 2Re(ρ_ii,bath_i)·sqrt(w_i/2)·∂Vg/∂q（双分支）] + ∂E_mf/∂q|_n
与 V(q,ρ) 的显式导数逐项一致：
- E0/g0（cau_peri_ground）：配对面 ±de、双谐振子核，精确。
- 杂质对角：dh_cache(1/2, k) 分别对应 h(1/2) = h̃ ± cf（晶体场劈裂进了力）✓。
- 耦合块：双对称元素都计入（09-29 修复：旧代码只算一半，耦合力减半）；vg_vec(2)=vg_vec(1)
  且 dvcp 同源，双杂质行共用 tr ✓。bath 对角 q 无关 ⟹ 无贡献 ✓。
- E_mf 核（mf_kernel）：dW = dU·ones、dJ 对角 = dU、off-diag = dj，与 E_mf 的张量结构一致，
  只取显式 q 依赖（n 冻结）——正是 HF 力的正确形态 ✓。
- Pulay 项：bath 基（离散能级）与 q 无关 ⟹ Pulay 力解析恒为零——"放弃 Pulay 项"是
  精确操作而非近似，有解析证书。
- 势能层（cau_peri_pot）：warp_z 的 dz 链式法则核对无误；Vcp 远场衰减的梯度
  （dv·(1−sw) − v·dsw，用未衰减值）自洽；冻结区（z≥6 Å）值与梯度同时归零 ✓。
- 界面重排：cau_peri_pes_r 与 cau_peri_prop_r 用同一 [C,S,Y]→[S,Y,C] 置换 ⟹ pes 与 prop
  在同一几何、同一 h 上工作 ✓（"同一势函数"约束成立）。

### C4 无未记账显式时间依赖 — 通过（mode 0；mode ≥ 1 违反，见 §2.1）
- 新容器无随机浴（ghost Y 是确定性谐振子，非 GASDEV 郎之万）；bath 参数、外场均无。
- 唯一随机源是 build_hcore_prop 的 rng_u（PERTURB_MODE ≥ 1 时）。

## 2. 残余解析违反项（全部给出公式与特征）

### 2.1 PERTURB_MODE ≥ 1（默认 0 = 关闭；已有代码内裁决记录）
prop 路径给 h_vec 加随机劈裂 Δ = diag(δ₁, −δ₂)（±perturb_amp，每 prop 调用重抽），pes 路径不加
（"pes/prop 对必须来自同一势函数"）⟹ 传播生成元 F̃ = F + Δ ≠ δV/δρ。解析漂移率：

  **Ė = (1/iℏ)·Tr([F, Δ]·ρ)** —— [F,Δ] 仅在杂质-bath 耦合块非零，漂移由杂化相干 Re(ρ_imp,bath) 驱动。

- mode 1（每次调用重抽）：能量做随机游走（每步跳变 ~ pamp·Δn），**任何步长都救不了**——
  这就是代码注释"non-energy-conserving stochastic law"的解析实质。
- mode 2（只抽一次）：V = E′ − Σρ_ii·δ_i，E′ 守恒 ⟹ V 只做有界振荡（幅度 ≤ ~2·pamp·2 ≈ 0.2 eV），
  无长期漂移。

### 2.2 Beeman 历史项的 ρ 滞后（唯一的积分器级解析缝隙；量级极小）
驱动时序（driver.f90:190-191）：`prop_step`（核步内所有 force_eval——含端点力——都读冻结的 ρ(t)）
→ `elec_prop`（端点几何上更新 ρ → ρ(t+dt)）→ `rec_frame`（e_tot 用**新** ρ 重算 V）。
记账核对：ΔV(q_{t+dt}, ρ_t→ρ_{t+dt}) 的电子部分被固定 q 幺正传播守恒 ⟹ 端点力用旧 ρ 不产生
显式失衡。**唯一混入旧 ρ 的是 Beeman 的 F(t−dt) 历史缓冲**（f_prev 上一步暂存，携带 ρ_{t−dt}）：

  **ΔE_step ≈ (dt²/12)·[∂F/∂ρ·ρ̇]·v**（AM-3 正确子中 −1/12 权重项的 ρ 滞后）。

特征：正比 dt（非 dt²）、**对 msbstp 完全不敏感**、符号由电子力演化与速度的相关性决定。
量级估计（dt=0.05×10fs、v~40 Å/ps、Ḟ_ρ~1–100 eV/Å/ps）：**~10⁻³–10⁻¹ eV/ps** —— 比观测的
−0.8 eV/ps 小 1–3 个量级，不足以解释主漂移，但在 msbstp 扫描的平台期应可辨。
判别实验：同 dt 下 INTEGRATOR order=1（VV，无历史项，解析上无此缝隙）vs order=2（Beeman）A/B。
注意 symple/radau 的全部踢力都用 ρ(t)（一致冻结），与 VV 同属"平衡"类，无此缝隙。

### 2.3 Magnus-2 子步截断（数值项，非解析违反）
显式中点预测-校正（fock 平均后对角化、从 rho_save 重演化）在固定 q 下守恒 Tr(ρF̄)，与监测 E 的
差为 O(dtq³)/子步 ⟹ 全局 O((dt/msbstp)²)——即回归测得的 msbstp⁻² 分量，属可收敛数值误差。

### 2.4 radau + elec_prop 的 dt 失配（潜在物理错误，非守恒错误）
driver 把固定的 `prop_cfg%dt` 传给 elec_prop，而 radau 的实际序列长 tt 是自适应的。
tdhf_hub 若配 radau，电子步长系统性错误（dt ≠ tt）。**不破坏能量记账**（任意步长的幺正传播
照样守恒固定 q 电子能），但传播的物理相位错。当前产线（verlet/beeman，固定 dt）不受影响；
建议在 elec_interface 或 driver 加一条命名守卫（radau 与电子方法互斥）。

## 3. 证书边界（未覆盖的暴露面）

1. **NN 路径的 V/∇V 一致性**（USE_NN=1 的产线输入）：probe_fulldof_fd 只 FD 了解析路径
   （cau_peri_load() 无参调用）。calc_nn_and_grad / NN 劈裂表 / calc_crystal_field 的 NN 分支
   未做逐项解析或 FD 核对。若 NN 梯度与其能量不是严格导数对，就是一个 O(1) 的 C3 违反——
   建议把 probe_fulldof_fd 加 NN 变体（load('nn_weights.txt')）跑一遍 9-DOF FD。
   rst vs cur 谱系的漂移差异若与此相关，该测试可直接定位。
2. 子系统能量（E3P 驻留段等通道量）：解析上只有 E_tot = T + V(q,ρ) 被守恒证书覆盖；
   通道能量的"渗漏"需用等大反号的转移项对账（电子↔核），对不上才是记账缺陷。
3. 首步 max_steps 耗尽后"第一个力看到上一条轨迹的 ρ"（头部注释已记录的残余）：影响 e0 基准
   一帧，有界，非长期项。

## 4. 单快照功率审计（免长传播的终判实验，按需执行）

在任意一条已有轨迹的相邻两帧上：
1. Ė_obs = [E_tot(t+dt) − E_tot(t)]/dt（从 trajectories.txt 直接取）；
2. 电子残差 = V_elec(q_{t+dt}, ρ_{t+dt}) − V_elec(q_{t+dt}, ρ_t)（重放一帧即可，不需要重跑）；
3. 冲量失配 = v̄·(dt/12)·[F(q_{t−dt},ρ_t) − F(q_{t−dt},ρ_{t−dt})]（Beeman 历史缝隙项）；
4. 显式时间依赖项 = ⟨∂Ĥ/∂t⟩（mode 0 时恒零）。
哪一项的大小与符号对上 Ė_obs，缺陷/物理就在哪一项。全程零新轨迹成本。

## 5. 总裁定

**PERTURB_MODE = 0 下，cau_peri 容器（tdhf_hub 方法包 + 组合 pes + verlet/vv/beeman/symple 驱动）
在解析层面满足能量守恒的全部四个条件**：传播生成元恰为监测能量的变分导数（09-28 修复）、
核力为监测势的精确 Hellmann–Feynman 导数（09-28/29 修复 + Pulay 恒零）、生成元厄米、
无未记账随机/外场项。残余的解析效应只有 Beeman 历史项 ρ 滞后（~10⁻³–10⁻¹ eV/ps，dt 线性、
msbstp 平坦）与 Magnus-2 截断（msbstp⁻²、可收敛）。**若在 mode 0 下仍观测到显著超出这两项
量级的 msbstp 平坦漂移，来源必在证书边界之外：NN 路径的 V/∇V 失配（§3.1，最可疑）或
所监测量为子系统通道能量（§3.2，物理交换而非缺陷）。**
