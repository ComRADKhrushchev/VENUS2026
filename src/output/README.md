# output/ — L4 记录与离线后处理

- **在线（演化环内）**：recorder.f90——轨迹序列化（rec_reg_col 列注册 +
  rec_select/rec_frame/rec_final）。演化环内合法调用四类之一
  （facade/propagator/recorder/predicate）。
- **离线（不参与演化环）**：final_state.f90（广义 COM 相对能量/角动量/终态量）、
  outcome.f90（结局分类执行——查 event%outcome_tbl 数据表）、density.f90
  （Beyer-Swinehart 态密度 + 本征输出）、plot_hist.f90（五类直方图，模块名 hist）。

边界约束：离线四件**不触物理态**（不修改 Q/P/s/F，不调 facade/propagator），
只做表象→表象的变换与归约（记录→统计）。

参考质料（旧库只读）：GWRITE.f90 / FINAL.f / GFINAL.f / VENUS.f90 结局块 /
DENQ.f / EIGOUT.f / HIST.f90。
