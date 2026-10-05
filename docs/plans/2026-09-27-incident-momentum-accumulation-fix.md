# 修复记录——incident 通道动量累加泄漏（2026-09-27）

> 动因：h3_bkmp2 转动激发生产（rotex, 10 能点 × 5000 轨）能量守恒检查发现 4996/5000
> 条轨迹末态能量爆炸（中位数 +2.7e5 kcal/mol，正常 ≈ -103）。帧级取证定位为初始态
> 缺陷而非积分缺陷；根因在 dist_incident 的动量写入语义。

## 1. 根因

`dist_incident.f90` 的 `add_mom`（inc_realize 的动量助手）把片段 COM 速度**累加**到
sta%p 上（`p += m*vel`），而 driver 的轨迹循环从不复位 sta（只复位 sta%t）。泄漏路径：

- 多原子片段：每条轨迹被其分布成员 realize 整体覆写 q/p（成员契约零净动量）——
  累加等价于赋值，无泄漏（五个在档 case、h2ag111/D2 腿全部属此类，不受影响）。
- **单原子片段**（nat<2，采样器按 2026-09-22 单原子跳过立法不派发内禀成员）：
  没有任何成员拥有其动量——弹丸动量逐轨迹叠加弹丸速度一次。

数值指纹（run_0, E=0.05 eV, p_beam=0.2546 amu·Å/(10 fs)）：

```
e0(N) = -103.5 + 0.775·N²  kcal/mol        （0.775 = p²/2m 弹丸动能）
N=6:  -75.3 (实测)  N=88:  +5782 (实测)  N=1638: +2.07e6 (实测, 预测 2.08e6)
```

伴随效应：终止步 ~6000/N（速度随 N 增长→更早到 R_TERM）、档案帧数调和衰减
（Σ6000/N ≈ 57466 与实测帧总数吻合）。冒烟（4 轨）泄漏仅 0.775×16≈12 kcal/mol
——不可见，故漏网。

## 2. 修复（赋值语义）

`add_mom` → `set_com_vel`：先减去片段当前 COM 速度再加目标速度
（`p += m*(vel - v_com_current)`）。多原子载体进来时 v_com=0，与旧 add 逐位相同；
单原子弹丸的跨轨迹残留被擦除、按处方赋值。头注两处同步（"momenta additive" →
"COM velocities assigned"）。改动仅 `src/methods/dist_incident.f90`。

## 3. 验证

1. **判据**：check_dist_incident 新增 (u) 腿——未复位态上二次装配必须复现新鲜态的
   位置与 COM 动量（旧 add 下必 FAIL）；全套重编重跑 **0 failure**（含新腿）。
   check_driver 重编重跑 **ALL PASS**。
2. **双编译器**：`make`（gfortran）与 `make FC=ifx` 全绿。
3. **五 case 回归**：八运行目录重跑——b_layer 三目录 + chain 两目录 trajectories.txt
   与直方图**逐字节一致**；三个 THERMALIZE case（harmonic_scatter / radau_conservation /
   surface_window）traj≥2 改变、守恒锚全部仍满足（漂移 1.4e-3 / 1.1e-13 / 5.4e-4，
   界 2.0e-1 / 1.0e-6 / 2.0e-1）。差异机制：THERMALIZE 平衡化从上一轨迹末态几何出发，
   N_EQ 步内片段在阱中获得 COM 振荡速度，旧 add 把束流叠加其上（污染处方 E_REL），
   新 set 按处方赋值（保留内部振动）。traj=1 逐位相同（阱底静止出发，COM 恒零）。
   **重入档登记**：三 case 的 results.txt 档案数值属旧缺陷态，按 E28 先例暂不改动，
   待相收口重入档。
4. **集群重跑**（修复后 venus.e，rotex2）：10 能点 × 5000 轨全 nonreactive，
   全轨迹漂移 ≤ 5.5e-2 kcal/mol；e0 逐轨迹分布 [-102.27, -101.89]（残差为初构型
   相互作用随 b/位相的物理变化；旧泄漏下应 ~+31000）。

## 4. 作废与影响面

- **rotex 第一轮全部作废**（分类计数、能量扫描俱不可信——N≳12 后每条轨迹实跑
  能量被平方级污染，远超势垒）。
- h2ag111 验证结果不受影响：H2/D2 弹丸均为双原子（EBK 整体覆写），路径上无泄漏
  （判据 3 的逐字节腿即其直接证据）。
- 转动激发验证改由 rotex2（无取向）+ rotex3（RANDOM_ORIENT=T 取向平均）承担。
