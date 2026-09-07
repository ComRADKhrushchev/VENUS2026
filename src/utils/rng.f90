!=====================================================================
! rng.f90 — L5 仪器：均匀随机流（Knuth 乘同余 + 游混）与分布变换（正态/伽马）
! 参考: src_VENUS/RANDST.f:47-105（流初始化）；src_VENUS/RAND0.f:19-23（游混取数）；
!       src_VENUS/RAND1.f:59-106（乘同余步进，基 256 八位数字）；
!       src_TEST/interface_TEST.f90:112-136（gasdev 正态对）；
!       src_VENUS/GAMA.f:7-13（伽马/指数）
!=====================================================================
module rng
   implicit none
   private
   public :: rng_init, rng_u, rng_gauss, rng_gamma

   ! 流状态（模块内部化——K1 计划 §2：旧 RANLST/ISEED3/IBFCTR 收归本模块私有）：
   integer, parameter :: dig_base = 256     ! 数字进制 2^8（旧 IBFCTR）
   integer :: seed_dig(8) = 0               ! 流种子 8 位基 256 数字（旧 ISEED3）
   real(8) :: tbl_shuf(100) = 0.0d0         ! 游混表 100 项（旧 RANLST）
   ! Box-Muller 缓存对状态（旧 gasdev 的 iset/gset）：
   integer :: i_cache = 0                   ! 缓存旗标（0=无缓存，1=有缓存；旧 iset）
   real(8) :: g_cache = 0.0d0               ! 缓存的第二个正态偏离（旧 gset）
contains
   !------------------------------------------------------------------
   ! rng_init(i_seed) — 由十进制种子建立新随机流（基 256 换算 + 演练 + 填游混表）
   ! 参考: src_VENUS/RANDST.f:47-105
   subroutine rng_init(i_seed)
      integer, intent(in) :: i_seed   ! 十进制种子（0 ≤ i_seed ≤ 2^31-1；旧标量 ISEED，
                                      ! control%i_seed 的折算映射 C4 输入定形）
      ! 伪代码：
      ! 1. i_seed 逐位除以 dig_base 折成 8 位数字存入 seed_dig（旧 :52-56）
      ! 2. seed_dig 翻倍并按位归一（≥dig_base 则减、向上进位），保证为偶数（旧 :65-76）
      ! 3. seed_dig(1) 加一并归一，保证流种子为奇数（旧 :81-90）
      ! 4. 连续两轮各抽 100 个数填满游混表 tbl_shuf（演练生成器；两轮同循环
      !    不同果，旧 :100-105）
   end subroutine rng_init

   !------------------------------------------------------------------
   ! rng_u() — 取下一个 [0,1) 均匀偏离（游混取数 + 乘同余步进补位）
   ! 参考: src_VENUS/RAND0.f:19-23；src_VENUS/RAND1.f:66-106
   real(8) function rng_u()
      ! 函数值：[0,1) 均匀随机数 [无量纲]
      ! 伪代码：
      ! 1. j = int(99*tbl_shuf(100)) + 1（游混指针，旧 RAND0.f:19）
      ! 2. 函数值 = tbl_shuf(100)；tbl_shuf(100) = tbl_shuf(j)（旧 :20-21）
      ! 3. tbl_shuf(j) = 乘同余步进值（旧 :22 调 RAND1）：
      ! 3a. 种子 8 位数字与乘数 8 位数字（Knuth 乘数 6364136223846793005 的
      !     基 256 表示 (/45,127,149,76,45,244,81,88/)，旧 RAND1.f:59）
      !     逐位卷积乘加，逢 dig_base 归一进位（中间 16 位工作数组，旧 :66-94）
      ! 3b. 低 8 位回写 seed_dig（旧 :96-98）
      ! 3c. 8 位数字按 256 进制折成浮点，再乘 1/dig_base 得 [0,1)（旧 :102-106）
   end function rng_u

   !------------------------------------------------------------------
   ! rng_gauss() — 标准正态偏离 N(0,1)（Box-Muller 极坐标形式 + 缓存对）
   ! 参考: src_TEST/interface_TEST.f90:112-136（gasdev）
   real(8) function rng_gauss()
      ! 函数值：标准正态偏离 [无量纲]
      ! 伪代码：
      ! 1. 缓存旗标为 1：函数值 = g_cache，清旗标返回（旧 :131-134）
      ! 2. 否则拒绝采样：抽 u 对 v=2*rng_u()-1，至 r2=v1^2+v2^2 落入 (0,1)
      !    （单位圆内均匀投点，旧 :118-125）
      ! 3. fac = sqrt(-2*ln(r2)/r2)（旧 :127）
      ! 4. g_cache = v1*fac 置旗标；函数值 = v2*fac（一对缓存，旧 :128-130）
   end function rng_gauss

   !------------------------------------------------------------------
   ! rng_gamma(order) — 整数阶伽马分布偏离（order=1 退化为单位指数分布）
   ! 参考: src_VENUS/GAMA.f:7-13
   real(8) function rng_gamma(order)
      integer, intent(in) :: order    ! 伽马分布阶数（正整数；1=指数分布，旧哑元 IA）
      ! 函数值：Gamma(order,1) 偏离 [无量纲]
      ! 伪代码：
      ! 1. dum = order 个 rng_u() 连乘（旧 :8-11；旧版混用内部 random_number，
      !    新版统一走 rng_u 流，保持流可复现）
      ! 2. 函数值 = -log(dum)（旧 :12）
   end function rng_gamma
end module rng
