! state.f90 —— 物理态变量（Q/P/s/F）容器：driver 持有、参数传递（K1 最小集；
!              E9 多面升级：elec_t 扩电子振幅与跃迁计数——非绝热通道的数据位；
!              E10 表象立法：amp 索引语义由乙类方法成员定义，表象构件不入 state）
module state
  implicit none
  private
  public :: elec_t, state_t, state_create, state_destroy
  type :: elec_t          ! 电子态 s（运行时）：绝热 = 退化点（单面）
     integer :: n_surf = 1   ! 势能面数
     integer :: i_surf = 1   ! 当前面索引
     ! E10 表象立法（注释立法，字段名不动）：amp 只承诺"复数组 + |·|² 可读作
     ! 布居"，索引语义由乙类方法成员定义——面表象 i=面号；单电子表象
     ! i=电子×态张量积展平（费米子反称化随成员自理）；组态表象 i=组态编号。
     ! 表象构件（H_ij/d_ij/V_i 等）不具物自体地位，不入 state_t（E10 公理 C
     ! ——只留门面工作区 ws，见 facade/elec_ws.f90；一步生命周期）
     complex(8), allocatable :: amp(:) ! 电子振幅 c_i(t)（|c_i|²=布居；
                                       ! 绝热退化点=单元素恒 1——退化点良定义保持）
     integer :: n_hop = 0    ! 累计跃迁次数 [次]（非绝热观测量；绝热恒 0）
  end type
  type :: state_t
     integer :: natoms = 0   ! 原子数（运行时由原子清单确定）
     integer :: ndof   = 0   ! 自由度数 = 3*natoms
     real(8), allocatable :: q(:) ! 坐标 [Å]
     real(8), allocatable :: p(:) ! 动量 [amu·Å/(10 fs)]
     real(8), allocatable :: f(:) ! 力 [内部能量/Å]
     type(elec_t) :: s
     real(8) :: t = 0.0d0     ! 轨迹时间 [10 fs]（主控程序持有）
  end type
contains
  ! 最小实现：按 natoms 分配 q/p/f 与电子振幅（ndof = 3*natoms；n_surf 缺省 1
  ! ——绝热退化点；失败写 stderr 并停机）
  subroutine state_create(st, natoms, n_surf)
     type(state_t), intent(out) :: st
     integer, intent(in) :: natoms
     integer, intent(in), optional :: n_surf ! 势能面数 [—]（缺省 1；
                                              ! 多面取 config%electronic%n_surf）
     integer :: n_s, stat
     n_s = 1
     if (present(n_surf)) n_s = n_surf
     if (n_s < 1) then
        write(0,'(a,i0)') 'state_create: n_surf 非法 (<1): ', n_s
        stop 1
     end if
     st%natoms = natoms
     st%ndof   = 3*natoms
     allocate(st%q(st%ndof), st%p(st%ndof), st%f(st%ndof), st%s%amp(n_s), stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_create: q/p/f/amp 分配失败, stat = ', stat
        stop 1
     end if
     st%q = 0.0d0
     st%p = 0.0d0
     st%f = 0.0d0
     st%s%n_surf = n_s
     st%s%i_surf = 1
     st%s%amp = (0.0d0, 0.0d0)
     st%s%amp(1) = (1.0d0, 0.0d0)  ! 基面布居 1（绝热退化点良定义）
     st%s%n_hop = 0
  end subroutine state_create

  ! 最小实现：释放 q/p/f/amp 并清零计数，失败写 stderr 并停机
  subroutine state_destroy(st)
     type(state_t), intent(inout) :: st
     integer :: stat
     deallocate(st%q, st%p, st%f, st%s%amp, stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_destroy: q/p/f/amp 释放失败, stat = ', stat
        stop 1
     end if
     st%natoms = 0
     st%ndof   = 0
     st%s%n_surf = 1
     st%s%i_surf = 1
     st%s%n_hop  = 0
  end subroutine state_destroy
end module
