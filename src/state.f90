! state.f90 —— 物理态变量（Q/P/s/F）容器：driver 持有、参数传递（K1 最小集）
module state
  implicit none
  private
  public :: elec_t, state_t, state_create, state_destroy
  type :: elec_t          ! 电子态 s（运行时）：绝热 = 退化点（单面）
     integer :: n_surf = 1   ! 势能面数
     integer :: i_surf = 1   ! 当前面索引
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
  ! 最小实现：按 natoms 分配 q/p/f（ndof = 3*natoms），失败写 stderr 并停机
  subroutine state_create(st, natoms)
     type(state_t), intent(out) :: st
     integer, intent(in) :: natoms
     integer :: stat
     st%natoms = natoms
     st%ndof   = 3*natoms
     allocate(st%q(st%ndof), st%p(st%ndof), st%f(st%ndof), stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_create: q/p/f 分配失败, stat = ', stat
        stop 1
     end if
     st%q = 0.0d0
     st%p = 0.0d0
     st%f = 0.0d0
  end subroutine state_create

  ! 最小实现：释放 q/p/f 并清零计数，失败写 stderr 并停机
  subroutine state_destroy(st)
     type(state_t), intent(inout) :: st
     integer :: stat
     deallocate(st%q, st%p, st%f, stat=stat)
     if (stat /= 0) then
        write(0,'(a,i0)') 'state_destroy: q/p/f 释放失败, stat = ', stat
        stop 1
     end if
     st%natoms = 0
     st%ndof   = 0
  end subroutine state_destroy
end module
