! control.f90 —— 控制流必须变量（独立模块）：主控编排输入，不含物理规定
module control
  ! 控制流输入：主控编排所需，不含物理规定（规范 §2.2：控制核无物理知识）
  implicit none
  public
  integer :: task    = 0   ! 任务方案（轨迹计算/读入坐标/…，输入 TASK 映射）
  integer :: n_traj  = 1   ! 轨迹条数（旧 NT）
  character(len=80) :: title(2) = ''  ! 运行标题
  integer :: i_seed(8) = 0 ! 随机数种子（8 字，rng 初始化用）
  ! 留白：断点续算、并行/进程管理等运行控制——将来按需引入（规范 §1 范围条件）
end module
