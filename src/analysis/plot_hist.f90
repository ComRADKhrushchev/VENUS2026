!=====================================================================
! analysis/plot_hist.f90 — L4 分析：五类直方图/落点图（归一参数化；格子常数
!                           全部传参——无体系常数）
!（模块名 hist：Fortran 不许模块与内含例程同名，例程名 plot_hist 为计划
!  §1 批 4 表定名——模块让名取域名）
! 参考: src_VENUS/HIST.f90:26-68（plot_Erel_hist：箱心 lo+(i−½)d :48-53、
!       [lo,hi) 计箱 int((v−lo)/d)+1 :56-60、TOTAL & P 归一行 :66）、
!       :70-93（plot_Erel_impact_pos：晶胞折叠落点 + (x,y,e_rel,z_min) 四列
!       写出 :87-92）、:95-131（plot_theta_hist——同形角度版）、:134-230
!       （plot_site_histograms：格位分组双直方图——top/brg/hcp/fcc 四组、
!       每组 TOTAL & P :219-221、n_bin≤1000 护栏 :167-172）、:232-311
!       （plot_Erel_by_bounce：反弹次数分组 + merge_at 合并阈 :265-274、
!       ≥merge_at 并末档 :283-296、n_bnc<0 跳过 :284）、:313-373（surfcell：
!       Wigner-Seitz 型最近格位归属——晶胞四角/五桥/二心候选点距 minloc
!       :340-372）、:375-404（Dist_cell：晶胞 PBC 折叠——格矢分数坐标取整
!       回叠 :393-403）；
!       src_VENUS/VENUS.f90:286-299（五函数调用位 + HISTO 目录编排——归
!       统计相位 driver）
! 剥离警示（体系知识写死——新版全传参，K2 计划取材注）：
!   (1) auro 与 PI/3 晶格常数写死（HIST.f90:19、:89、:188——Au(111) 体系）
!       → lat_t 传参（与 roster 晶胞字段同源——a_lat/skew）；
!   (2) 格位名 'top'/'brg'/'hcp'/'fcc' 与三角格四角/桥/心几何写死（:165、
!       :340-345）→ 格位候选集数据化（C8 定形；骨架保留四名对照）
!=====================================================================
module hist
   implicit none
   private
   public :: hist_spec_t, hist_data_t, lat_t, plot_hist
   public :: k_e_hist, k_t_hist, k_pos, k_site, k_bnc

   ! 五类规格键（spec%kind 取值——旧五个 plot_* 函数的参数化归一）
   integer, parameter :: k_e_hist = 1  ! 类 1：能量直方图（旧 plot_Erel_hist）
   integer, parameter :: k_t_hist = 2  ! 类 2：角度直方图（旧 plot_theta_hist）
   integer, parameter :: k_pos    = 3  ! 类 3：落点-量图（旧 plot_Erel_impact_pos）
   integer, parameter :: k_site   = 4  ! 类 4：格位分组直方图（旧 plot_site_histograms）
   integer, parameter :: k_bnc    = 5  ! 类 5：按反弹次数分组直方图（旧 plot_Erel_by_bounce）

   ! 直方图规格（量程/箱数/归一——全部输入数据）
   type :: hist_spec_t
      integer :: kind = 0                  ! 五类键（k_*）
      character(len=64) :: file_out = ''   ! 输出文件名（旧 'HISTO/…'——目录
                                           ! 编排归 driver；分组类派生名
                                           ! = file_out 前缀 + 组名）
      real(8) :: v_lo = 0.0d0              ! 量程下限 [量的单位]（旧 Erel_min/theta_min）
      real(8) :: v_hi = 0.0d0              ! 量程上限（旧 Erel_max/theta_max）
      integer :: n_bin = 0                 ! 箱数 [个]（≤1000 护栏，旧 :167-172/:265-268）
      integer :: merge_at = 0              ! 类 5 合并阈 [次]（反弹次数 ≥ 此值并末档）
      logical :: norm_p = .false.          ! 末行 TOTAL & P 归一（旧 :66/:220）
   end type hist_spec_t

   ! 轨迹行表（recorder 终态行的聚合——旧 AErelS/AInitX 等全局数组的成员化）
   type :: hist_data_t
      integer :: n_traj = 0                    ! 轨迹行数 [行]（旧 NT）
      real(8), allocatable :: e_rel(:)         ! 终态相对能列 [kcal/mol]（旧 AErelS）
      real(8), allocatable :: thta(:)          ! 散射角列 [°]（旧 AThtaS = ANG(19)）
      real(8), allocatable :: x0(:)            ! 初始落点 x 列 [Å]（旧 AInitX = RX0）
      real(8), allocatable :: y0(:)            ! 初始落点 y 列 [Å]（旧 AInitY = RY0）
      real(8), allocatable :: z_min(:)         ! 最小高度列 [Å]（旧 AZminS）
      real(8), allocatable :: v_fin(:)         ! 终态势能列 [kcal/mol]（旧 AVfinS）
      integer, allocatable :: n_bnc(:)         ! 反弹次数列 [次]（旧 NBncS；−1=无效位）
      character(len=3), allocatable :: site(:) ! 格位列 [—]（旧 site——record_all :19）
   end type hist_data_t

   ! 格子常数（全部传参——无体系常数；与 roster 晶胞字段同源）
   type :: lat_t
      real(8) :: a_lat = 0.0d0     ! 晶格常数 [Å]（旧 auro 写死——剥离警示 (1)）
      real(8) :: skew = 0.0d0      ! 晶胞斜角 [rad]（旧 PI/3 写死——剥离警示 (1)）
      real(8) :: q_ref(2) = 0.0d0  ! 参考格点坐标 [Å]（旧 Qref——晶胞原点；
                                   ! 旧调用位 Qr2 = 2.55476 写死 VENUS.f90:290）
   end type lat_t
contains
   !------------------------------------------------------------------
   ! plot_hist(spec, data, lattice) — 五类直方图/落点图分发（归一参数化）
   ! 参考: HIST.f90:26-311（五函数）、:313-404（surfcell/Dist_cell 两工具）；
   !       格位候选集数据化前，类 4/格位列按四名对照执行——C8 定形
   subroutine plot_hist(spec, data, lattice)
      type(hist_spec_t), intent(in) :: spec      ! 直方图规格（类键/量程/箱数/归一）
      type(hist_data_t), intent(in) :: data      ! 轨迹行表（recorder 聚合）
      type(lat_t), intent(in) :: lattice         ! 格子常数（传参——无体系常数）
      ! 伪代码：
      ! 1. 箱几何：d = (v_hi − v_lo)/n_bin；箱心 c_i = v_lo + (i−0.5)·d
      !    （旧 :48-53/:175-180）；护栏：n_bin>1000 截断、merge_at<1 抬 1
      !    （旧 :167-172/:265-273）
      ! 2. 按 spec%kind 分发：
      ! 2a. k_e_hist / k_t_hist（旧 :50-67/:115-130）：取 data%e_rel 或
      !     data%thta，落 [v_lo,v_hi) 计箱 idx = int((v−v_lo)/d)+1（越界不计，
      !     旧 :56-60）；逐箱写 (c_i, count)；norm_p 时末行
      !     'TOTAL & P'（Σcount、Σcount/n_traj——旧 :66）
      ! 2b. k_pos（旧 :87-92）：逐轨迹落点晶胞 PBC 折叠（同 2c 首步）后写
      !     (x, y, e_rel, z_min) 四列（旧 write(98,'(4F12.6)') :90）
      ! 2c. k_site（旧 :185-229）：逐轨迹先折叠到参考晶胞——格矢
      !     e_x=(a_lat,0)、e_y=(a_lat·cos(skew), a_lat·sin(skew))，分数坐标
      !     a1/a2 取整回叠（旧 Dist_cell :393-403）；再最近格位归属——晶胞
      !     四角/五桥边心/两面心候选点距 minloc，映射组名（旧 surfcell
      !     :340-372；候选集数据化 C8）；按格位分组各自计箱写派生名文件；
      !     norm_p 每组 TOTAL & P（除该组轨迹数——旧 :219-221）
      ! 2d. k_bnc（旧 :282-309）：n_bnc<0 跳过（旧 :284）；n_bnc ≥ merge_at
      !     并末档（旧 :285-289）；逐档写直方图文件（档名
      !     _bounce_{0..merge_at−1} 与 _bounce_ge{merge_at}——旧 :298-309）
      ! 3. 输出单元 open(file_out, status='replace') … close（旧 :62-67 位）；
      !    HISTO 目录编排归 driver 统计相位（旧 VENUS.f90:293）
   end subroutine plot_hist
end module hist
