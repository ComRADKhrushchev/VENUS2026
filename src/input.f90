!=====================================================================
! input.f90 — L0 输入：KEYWORD=VALUE 框架（逐行解析 → control/config 槽映射
!             与成员参数分发；文法定案 C4——骨架先立框架与映射位）
! 参考: venus2026/input_parser.f90:67-152（load_input：空行与 '#'/'!' 注释
!       跳过 :107-108、'=' 划分键/值 :125-130、值段行尾 ',' 续行吸收
!       :134-145、键大写归一 :170、首数据行 '=' 有无判文法 :111-120）、
!       :196-239（extract_nth 逗号分段取第 n 值）、:268-458（get_int/
!       get_real/get_str/get_real_arr/get_int_arr/has_keyword/
!       count_keyword_values/get_*_element 取数族）、:489-502（set_default——
!       缺省键只在用户未设时注入）；
!       venus2026/param_mapping.f90:34-160（枚举映射表族：map_elec_method
!       :34-55、map_surface_model :60-75、map_task :80-101、map_integrator
!       :106-121+、map_init_sampling、map_old_nsurf——字符串枚举→整数编号，
!       非法值指名停机）；
!       venus2026/presets.f90:33-100（apply_preset：MODEL 键→五预设分发、
!       预设体=set_default 缺省键表）；
!       venus2026/venus_input.f90:187-204（装载 + MODEL 预设 + 旧序贯文法
!       拒绝位 :200-204）、:209-231（TITLE×2/ELEC_METHOD 新键与 CALTYP 旧键
!       回退双轨 :214-222/NATOMS/GWRITE_LEVEL）、:248-256（ATOM_MASSES 质量
!       表 + 两值广播规则）、:403-413（INIT_SAMPLING_A 新键与 NACTA 旧整数
!       回退双轨——成员参数分发原型）、:495-506（NACTA 分支打印——分发位
!       雏形）、:607-681（TRV_A 等成员参数逐分支读取）、:996-1000（SKEW/
!       BOXLX/BOXLY 晶胞段）；
!       文法定案 C4（K2 计划 §1 批 4 表注）——新键集、TERMINATE=/OUTCOME=
!       表文法与旧整数键兼容表 C4 定稿；本骨架立框架与映射位，不锁键集
!=====================================================================
module input
   use control, only: task, n_traj, title, i_seed
   use config,  only: propagator, electronic, reactants, event, observables
   implicit none
   private
   public :: read_input

   ! 解析暂存（KEYWORD=VALUE 行表——模块私有；roster_load 与成员 init 的
   ! 消费接口 C4 定形）
   type :: kv_row_t
      character(len=32)  :: key = ''   ! 关键字（大写归一；旧 MAX_KEY_LEN 缩编）
      character(len=256) :: val = ''   ! 值段（逗号分段；续行已并回）
   end type kv_row_t
   type(kv_row_t), allocatable :: kv_tbl(:)  ! 行表（模块私有）
   integer :: n_kv = 0                          ! 行数 [行]
contains
   !------------------------------------------------------------------
   ! read_input(file) — 读入与映射输入（解析 → 预设 → 枚举 → 槽映射 →
   !                   名册暂存 → 成员参数分发 → 校验）
   ! 参考: venus_input.f90:187-269（读入主流程原型）；input_parser.f90:67-152
   subroutine read_input(file)
      character(len=*), intent(in) :: file  ! 输入文件名（如 'input_qct.txt'——
                                            ! 旧 venus_input.f90:189 写死名，
                                            ! 新版哑元传入）
      ! 伪代码：
      ! 1. 逐行解析入 kv_tbl：空行与注释行（'#'/'!' 首字符）跳过；'=' 划分
      !    键/值两段；值段行尾 ',' 续行吸收（并回下一行）；键大写归一
      !    （旧 input_parser.f90:97-148；行号计入错误通道）
      ! 2. 文法检测与旧整数关键字回退：首数据行无 '=' → 旧序贯文法——按旧
      !    整数关键字流回退读取（旧 venus_input.f90:200-204 为硬拒绝；新版
      !    回退读法或转换指引 C4 裁决定形——骨架立"检测 + 回退位"槽，
      !    不写死策略）
      ! 3. 预设展开：MODEL 键 → 预设分发（五预设——气相/弛豫表面/刚性表面/
      !    全表面/GLO 表面），预设体 = 缺省键注入（不覆盖用户已设键——
      !    set_default 语义 input_parser.f90:489-502；旧 presets.f90:33-100）
      ! 4. 枚举映射（映射表族思想——param_mapping.f90）：字符串枚举→整数
      !    编号，非法值指名停机——电子律→electronic%method（旧
      !    map_elec_method）、构型范式→reactants%surface_model、任务→task、
      !    积分器+阶→propagator%id/order（旧 map_integrator）、抽样方案→
      !    reactants%dist_scheme（旧 map_init_sampling）；无对应新枚举键时
      !    旧整数键直读兼容（CALTYP/NSURF/INTEGRATOR/NACTA·NACTB——旧
      !    venus_input.f90:214-222/:403-413 双轨回退位）
      ! 5. 控制槽映射：task、n_traj（旧 NT）、title(2)（旧 TITLE1/2）、
      !    i_seed 折算 8 字（旧 ISEED 标量→rng 种子折算映射 C4——rng_init
      !    消费）
      ! 6. config 槽映射：
      !    propagator ← id/order/dt（旧 INTEGRATOR/LLL/HTMIN 族）；
      !    electronic ← method/n_surf/pes_main（旧 CALTYP/NSURF/TEST_PES）；
      !    reactants ← n_frag/surface_model/dist_scheme(2)/e_rel/b_max/r_sep
      !    （旧 NACTA·B/NSURF/SEREL/BMAX/S）；
      !    event ← term_expr/outcome_tbl（TERMINATE=/OUTCOME= 新文法键——
      !    predicate 模块头文法；C3/C7 接入）；
      !    observables ← rec_level/obs_list/e_kin/e_pot/e_tot（旧
      !    GWRITE_LEVEL→rec_level）
      ! 7. 名册原始段暂存：ATOM_MASSES 质量表 + 两值广播规则（只给 2 值则
      !    广播到其余原子——旧 venus_input.f90:248-256）、片段 A 表与平衡
      !    几何（旧 :508-518）、晶胞参数 SKEW/BOXLX/BOXLY（旧 :996-1000）；
      !    暂存数据由 roster_load 消费（消费接口 C4 定形——roster 不引
      !    本模块，方向单向）
      ! 8. 成员参数分发：dist_*/pes_* 各成员自带参数键 → 各成员 init 分发
      !    （成员 = 代码 + 自带数据；旧 TRV_A/TROT_A/AI_A 等分散于 NACTA
      !    各分支 venus_input.f90:607-681——新版集中分发表，C2/C5 随成员
      !    接线逐行填充）
      ! 9. 校验：必填键缺位/数值域违例/枚举非法 → loud-stop 指名键与行号
      !    （错误通道口径同 pred_parse 的 ios/errmsg）
   end subroutine read_input
end module input
