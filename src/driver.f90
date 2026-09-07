!=====================================================================
! driver.f90 — L0 主控：program venus——三相位编排（装配 → 演化 → 统计）；
!               物理态持有者（state_t 参数传递——K1 计划 D1）
! 参考: venus2026/VENUS.f90:71-115（初始化数组清零 + 读入段——装配相位原型）、
!       :157-167（积分参数装配——propagator 槽的旧位）、:182-232（随机流
!       初始化 RANDST :184/:231——rng_init 位）、:239-705（TRAJECTORY_LOOP ×
!       INTEGRATION_LOOP 两环结构——演化相位骨架原型）、:316-349（初值装载
!       与 SELECT 调用 :349——samp_run 位）、:402-413（初力/首帧 GWRITE/
!       QSAVE-PSAVE 暂存——轨迹环首原型）、:419-465（积分环 INTEGRATOR 三
!       分支分派 :423-465——prop_step 统一）、:505-576 与 :581-662（NFINAL/
!       NTEST 两出口的结局块）、:669-699（长时程出口）、:243-301（完成段：
!       结局汇总与直方图五调用 :286-299——统计相位原型）
! 剥离清单（参考时勿带入，K2 计划取材注——体系分支/双份结局块/死块全剥离）：
!   - RKUTTA 系数死块 :86-94（RKUTTA 在 K1 §2 死块清单——永不进入）
!   - 断点续算分支 :182-232（NCHKP——control 留白，将来按需引入）
!   - 体系分支 :322-346（NSURF==2·NGLO 表面四原子装载/FIXROTDATM 对齐/
!     THTA 度弧分叉——构型范式数据化 config%reactants%surface_model，
!     细节归 sampler/roster 成员）
!   - 双份结局块 :505-576 与 :581-662（两入口的 FINAL/GFINAL/record_all/
!     反应计数几乎逐行全同——归并为本骨架统计相位单次逐轨迹处理）
!   - 反弹/最小高追踪 :453-464（BOUNCE/z_pos_min——体系诊断量，归 recorder
!     列注册）
!   - 模布居追踪 :363-395/:468-502（MPLOT/ENMODE/EBOND——归 analysis 留白 C8）
!   - 转动能重复抽样平均 :511-542/:598-629（NCROT ERAVA/ERAVB 双份——
!     fin_t 留白/统计相位）
!   - 能量守恒追踪与门 :451-452/:461-464/:559-565（HTMAX-HTMIN——归 recorder
!     能量列聚合 + outcome 表行）
!   - STOP/计时双份打印 :243-301（999/6 双通道逐行重复——归并统计段一次）
!=====================================================================
program venus
   use control,        only: task, n_traj, title, i_seed
   use config,         only: event, observables, reactants
   use config_roster,  only: roster, roster_load
   use state,          only: state_t, state_create, state_destroy
   use rng,            only: rng_init
   use predicate,      only: expr_t, pred_parse, pred_eval
   use propagator,     only: hooks_t, prop_reg, prop_step, &
                             verlet_step, symple_step, radau_step
   use force_facade,   only: facade_reg_pes, facade_reg_elec, &
                             facade_init, facade_force
   use elec_adiabatic, only: ad_dsdt, ad_force_rule
   use sampler,        only: samp_reg, samp_run
   use recorder,       only: rec_select, rec_frame, rec_final
   use final_state,    only: fin_t, fin_eval
   use outcome,        only: outc_tbl_t, outc_classify
   use density,        only: dens_states, eig_out
   use hist,           only: hist_spec_t, hist_data_t, lat_t, plot_hist
   use input,          only: read_input
   implicit none

   ! —— 编排计数与旗标（局部变量——E7；旧 NTZ/NC/NX 全局编排量收归于此）——
   integer :: traj_i        ! 轨迹序 [条]（旧 NTZ；1..n_traj）
   integer :: step_i        ! 步序 [步]（旧 NC——pred_eval 的 MAX_STEPS 通道）
   logical :: term_flag     ! 本轨迹终止旗标（hooks%term 闭包宿主持有——E7
                            ! 演化环终止通道，propagator.f90 注 3 同源）
   real(8) :: e0            ! 轨迹初态总能 [kcal/mol]（facade_force 出口——
                            ! 守恒漂移基准，旧 E0）
   real(8) :: t_cpu(2)      ! 计时对 [s]（旧 START_TIME/FINISH_TIME——统计段打印）
   integer :: u_out         ! 记录单元号 [—]（帧/终态行共用；旧 6/8/1000+NTZ
                            ! 族的单元编排 C4 定形）
   character(len=32) :: outc_name  ! 本轨迹结局名 [—]（统计相位——outc_classify 出口）
   ! —— 相位间承接的数据结构 ——
   type(state_t) :: sta                      ! 物理态（driver 持有——K1 D1；
                                             ! 演化相位独占写、统计相位只读）
   type(state_t), allocatable :: sta_fin(:)  ! 各轨迹终态快照 [条]（演化→统计
                                             ! 承接；旧 QSAVE/PSAVE 暂存位同源)
   type(hooks_t) :: hooks                    ! 检查点回调束（rec_frame/pred_eval 接线）
   type(expr_t)  :: expr_term                ! 终止谓词树（装配相位 pred_parse 产物）
   type(outc_tbl_t) :: outc_tbl              ! 结局分类表（outcome_tbl 解析产物——C7）
   type(fin_t) :: fin                        ! 终态量（统计相位逐轨迹承接）
   type(hist_spec_t) :: spec_hist            ! 直方图规格（统计相位装配——C8）
   type(hist_data_t) :: data_hist            ! 轨迹行表（recorder 终态行聚合——C8）
   type(lat_t) :: lattice                    ! 格子常数（roster 晶胞映射——统计相位）
   integer, allocatable :: cnt_outc(:)       ! 各结局计数 [条]（与 outc_tbl 行对齐；
                                             ! 旧 NONREACT/NREACT/NLONG 计数组）

   ! 伪代码（三相位——E7；规范 §2.5 相位限制表：装配=操作类+初始化性测量、
   ! 演化=演化类+测量类、统计=后处理类）：
   !
   ! ===== 相位 1：装配 =====
   ! 1. read_input('input_qct.txt')——输入框架（KEYWORD=VALUE → control/config
   !    槽/名册暂存/成员参数；文法 C4；旧 VENUS.f90:115 read_venus_input 位）
   ! 2. rng_init(i_seed 折算)——随机流建立（旧 :184/:231 RANDST(ISEED)；
   !    control%i_seed 8 字→i_seed 折算映射 C4——input 步 5）
   ! 3. roster_load()——原子清单装配（质量表/片段 A 表/平衡几何/晶胞；
   !    旧 :316-335 初值装载位的名册化）
   ! 4. state_create(sta, size(roster%mass), electronic%n_surf)——物理态分配
   !    （natoms 由名册定；n_surf 传 config 值——电子振幅按面数分配，缺省 1
   !    绝热退化点；E9）
   ! 5. 成员注册（各一行——注册族装配期约定）：
   !    prop_reg(1, verlet_step)、prop_reg(2, symple_step)、
   !    prop_reg(3, radau_step)（封闭三成员——X4，统一旧 :423-465 三分支）；
   !    facade_reg_elec('adiabatic', …)——乙类退化成员（ad_dsdt/ad_force_rule
   !    到槽的适配接线 C2）；facade_reg_pes(主势名, …)——pes_* 成员 C2 填充；
   !    多面主势成员（E9 双签名）：facade_reg_pes_diab(名, …)/
   !    facade_reg_pes_adiab(名, …)——随 elec_tsh/iesh 消费者引入（注册行
   !    同为一行，T2 不变）；samp_reg(分布成员名, …)——dist_* 成员 C5 接线
   !    （boltzmann/ebk/incident 批 3 骨架 + 余 9 个）
   ! 6. facade_init()——成员装配 + 能力级三元组校验（非法组合 loud-stop）
   ! 7. rec_select(observables%rec_level)（或 rec_select(observables%obs_list)
   !    名单入口）——记录列组选择（注册先于选列）
   ! 8. pred_parse(event%term_expr, expr_term)——终止树装配；无 TERMINATE 键
   !    时默认表达式（同参版——C3 装配层构造，predicate.f90 注 5）；
   !    event%outcome_tbl → outc_tbl 解析（C7 文法定形件）；spec_hist/
   !    lattice 装配（直方图规格与 roster 晶胞映射——C8）
   ! 9. hooks 接线（E7 终止通道——propagator.f90 注 3）：hooks%out => 输出
   !    闭包（宿主内部调 rec_frame(sta, u_out)——nopass 指针 + 宿主持单元号，
   !    闭包宿主例程形态 C3 定形）；hooks%term => 终止闭包（内部调
   !    pred_eval(sta, expr_term, step_i)，真则置 term_flag——旗标本程序持有）
   ! 10. allocate(sta_fin(n_traj)、cnt_outc(outc_tbl%n_row)、data_hist 列——
   !     相位间承接结构装配；计时起点 t_cpu(1)
   !
   ! ===== 相位 2：演化（轨迹环 × 积分环）=====
   ! 【相位限制声明（E7——规范 §2.5）】：积分环（下方内环）内只出现
   !    facade/propagator/recorder/predicate 四类调用——直接调用仅
   !    prop_step（propagator 类）；facade_force（facade 类）由传播子成员
   !    内部求值；rec_frame（recorder 类）经 hooks%out 检查点回调；
   !    pred_eval（predicate 类）经 hooks%term 闭包。环内无 sampler/
   !    analysis/input/roster 调用；编排计数 traj_i/step_i 为局部变量。
   ! 轨迹环：do traj_i = 1, n_traj（旧 TRAJECTORY_LOOP :239）
   ! 1. samp_run(scheme, sta)——逐轨迹初态制备（操作类调用——规范 §2.3
   !    边界 4：装配相位的实际调用含分布成员，制备成员属装配性调用；
   !    位于轨迹环首、积分环外；旧 :349 SELECT 位。scheme 取
   !    reactants%dist_scheme 映射的成员名——C5 接线）
   ! 2. facade_force(sta, e0)——初力与初态总能（旧 :402-403 DVDQ_1/ENERGY_1）
   ! 3. rec_frame(sta, u_out)——首帧（旧 :405 GWRITE 位）
   ! 4. step_i = 0、term_flag = .false.（编排复位）
   ! 5. 积分环：do while (.not. term_flag)（旧 INTEGRATION_LOOP :419）
   ! 5a. step_i = step_i + 1（编排计数——MAX_STEPS 安全网通道）
   ! 5b. prop_step(sta, hooks)——推进一检查点（成员内部：facade_force
   !     多次 + 检查点回调先 out（rec_frame）后 term（pred_eval→term_flag）；
   !     旧 :423-465 三分支 + :450/:459 间隔输出的统一位）
   ! 6. 环尾（积分环外、轨迹环尾）：sta_fin(traj_i) = sta——终态快照存档
   !    （旧 :411-413 QSAVE/PSAVE 暂存位同源）；rec_final(sta, u_out)——
   !    终态行（测量类；recorder 终态行含体系诊断列——bounce/z_min 等，
   !    剥离清单第 5 条的落位）
   !    （旧双份结局块 :505-576/:581-662 的 FINAL/GFINAL/record_all/反应计数
   !    与长时程出口 :669-699——全部延后统计相位处理，见剥离清单）
   !
   ! ===== 相位 3：统计（后处理类——outcome/analysis 调用）=====
   ! 1. 逐轨迹：do traj_i = 1, n_traj
   ! 1a. fin_eval(sta_fin(traj_i), fin)——终态量（旧 FINAL 广义部分）
   ! 1b. outc_name = outc_classify(fin, outc_tbl)；cnt_outc 命中行 +1
   !     （旧 :552-571/:645-657 反应计数与 :669-685 长时程计数——数据化）
   ! 1c. data_hist 聚合（fin/recorder 终态行 → e_rel/thta/x0/y0/z_min/
   !     v_fin/n_bnc/site 列——旧 GFINAL.f:199-202 与 record_all 的装配位，
   !     C8 接线）
   ! 2. analysis 调用：plot_hist(spec_hist, data_hist, lattice)（五类直方图
   !    ——旧 :286-299 plot_* 五调用位；HISTO 目录编排在此）；需要时
   !    dens_states/eig_out（态密度/本征输出）
   ! 3. 汇总打印：结局计数表/有效轨迹数/能量守恒漂移汇总/计时（旧
   !    :243-284 双通道归并一次；守恒门 verdict 交 outcome 表行）
   ! 4. state_destroy(sta)、释放 sta_fin/cnt_outc/data_hist——收尾
end program venus
