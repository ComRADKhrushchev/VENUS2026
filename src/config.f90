! config.f90 —— 表象层配置：可观测量 + 固化操作导出的四类输入规定（type 封装 + 留白）
module config
  implicit none
  private
  public :: propagator_t, electronic_t, reactants_t, event_t, observables_t
  public :: propagator, electronic, reactants, event, observables

  ! ===== 固化操作导出的输入规定（规范 §3.1-3.2）=====

  ! 核传播子（规范 §3.2-1：功能有界，封闭族固化）
  type :: propagator_t
     integer :: id    = 1     ! 积分器编号（verlet/symple/radau 封闭三成员）
     integer :: order = 0     ! 阶数/变体参数
     real(8) :: dt    = 0.01d0 ! 步长 [10 fs]
     ! 留白：自适应步长界 dt_min/dt_max（RADAU 具体化时引入，C6）
  end type

  ! 电子层面（规范 §3.2-2：两独立轴——H 表示能力级 × 电子传播律）
  type :: electronic_t
     integer :: method = 1       ! 电子传播律（1=绝热退化点；非绝热预留）
     integer :: n_surf  = 1      ! H 表示能力级（单面；多面/耦合留待非绝热）
     character(len=32) :: pes_main = ''  ! 主势成员名（H 场提供方，接口 A）
     ! 留白：加性力成员表（GLO 弹簧/朗之万 ghost/LDFA——C2 随门面引入）；
     !       介电/外场等能力级扩展
  end type

  ! 反应物规定（规范 §3.2-3：构型规定 + 初态统计分布，皆输入数据）
  type :: reactants_t
     integer :: n_frag = 1           ! 片段数
     integer :: surface_model = 0    ! 构型范式（0=气相；1/2=表面两范式）
     integer :: dist_scheme(2) = 0   ! 各片段抽样方案（指向 dist_* 成员）
     real(8) :: e_rel  = 0.0d0       ! 碰撞能（相对运动分布）[kcal/mol]
     real(8) :: b_max  = 0.0d0       ! 最大碰撞参数 [Å]
     real(8) :: r_sep  = 0.0d0       ! 初始间距 [Å]
     ! 构型细节（名册/几何/晶胞）由 config_roster 承载（K2 骨架）；
     ! 留白：取向规定、各片段分布参数——由 dist_* 成员自带（成员=代码+数据）
  end type

  ! 事件规定（规范 §3.2-3：终止谓词 + 结局分类，数据非代码）
  type :: event_t
     character(len=256)  :: term_expr   = '' ! 终止谓词表达式（C3 接入）
     character(len=1024) :: outcome_tbl = '' ! 结局分类表（C7 接入）
     ! 留白：事件类别扩展（约束窗/多事件耦合——谓词文法定形时补）
  end type

  ! ===== 可观测量（通用能量项与记录设置）=====
  type :: observables_t
     integer :: rec_level = 2          ! 记录详细度（预设列组）
     character(len=256) :: obs_list = '' ! 观测列名单（C4 记录器列注册）
     logical :: e_kin = .true.         ! 动能 T（通用能量项）
     logical :: e_pot = .true.         ! 势能 V（通用能量项）
     logical :: e_tot = .true.         ! 总能 H（通用能量项）
     ! 留白：更多观测量随记录器列注册扩展（角动量、局域量等，C4/C8）
  end type

  type(propagator_t) :: propagator
  type(electronic_t) :: electronic
  type(reactants_t)  :: reactants
  type(event_t)      :: event
  type(observables_t) :: observables
end module
