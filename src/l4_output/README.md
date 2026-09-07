# analysis/ — L4 离线后处理族

L4 后处理层：对已记录的轨迹/观测数据做**结局分类执行、统计、绘图**——
可离线运行，不参与演化环。

边界约束：**不触物理态**（不修改 Q/P/s/F，不调 facade/propagator），
只做**表象→表象**的变换与归约（记录→统计）。演化期驻留的全部测量
（轨迹序列化 + 谓词求值）在 L4 recorder.f90，不在本目录。

## 计划成员（K 阶段生成骨架，C8 具体化）

| 成员 | 用途 | 参考质料（旧库，只读） |
|---|---|---|
| 结局分类执行器 | 消费 driver 的 `OUTCOME = 名字: 谓词表达式` 分类数据，逐轨迹判类计数；广义 COM/EREL/OAM 计算 | FINAL.f / GFINAL.f |
| 态密度/本征输出 | DENQ+SUMQ、EIGOUT（DEL=1 cm⁻¹ 配置化；fort.77 通道复用拆分） | DENQ.f / EIGOUT.f |
| 离线绘图 | HIST 的 plot_* 五函数（surfcell/Dist_cell 通用实现保留，格子常数传参化） | HIST.f90 |
| 诊断工具 | nearest-atom 诊断（自 recorder 剥离的可选离线工具）；conserv_check、smoke_rst_bvk 验收工具纳入构建 | GWRITE.f90 / scripts/ |

执行依据：[docs/final_completion_plan.md](../../docs/final_completion_plan.md) §2（C8）、
§7（X6）。
