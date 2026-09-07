! consts.f90 —— 单位契约：全程序共享的换算常数与数学常数最小集（结构上限随消费者自带）
module consts
  implicit none
  private
  public :: e_conv, r_kcal, kb_code, hbar_code, pi, half_pi, two_pi, dtor, e2wvn, wvn2e

  ! 依 §1.1 表格逐项（行末括注旧名；pi 族先声明以供 dtor 表达 π/180）：
  real(8), parameter :: e_conv    = 0.04184d0                ! kcal/mol→内部能量单位（旧名 C1）
  real(8), parameter :: r_kcal    = 1.9872198404d-3          ! 气体常数 [kcal/mol/K]（旧名 C8）
  real(8), parameter :: kb_code   = 0.083144d-3              ! 玻尔兹曼常数内部单位组合（旧名 C5）
  real(8), parameter :: hbar_code = 0.063508d0               ! ℏ 的内部单位值（旧名 C7）
  real(8), parameter :: pi        = 3.14159265358979323846d0 ! 圆周率 π（旧名 PI）
  real(8), parameter :: half_pi   = pi/2.0d0                 ! 圆周率族：π/2（旧名 HALFPI）
  real(8), parameter :: two_pi    = 2.0d0*pi                 ! 圆周率族：2π（旧名 TWOPI）
  real(8), parameter :: dtor      = pi/180.0d0               ! π/180，度→弧度（旧名 C4）
  real(8), parameter :: e2wvn     = 349.755d0                ! kcal/mol↔cm⁻¹：正向换算（旧名 cal2cm）
  real(8), parameter :: wvn2e     = 1.0d0/e2wvn              ! kcal/mol↔cm⁻¹：e2wvn 之倒数（旧名 cm2cal）
end module consts
