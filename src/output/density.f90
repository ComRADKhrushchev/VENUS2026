!=====================================================================
! analysis/density.f90 — L4 分析：量子态密度（Beyer-Swinehart 直接计数）与
!                         本征值/向量格式化输出
! 参考: src_VENUS/DENQ.f:10-18（DENQ：两探针能级 E±50 各求累计态数，密度 =
!       (SUM1−SUM2)/100 中心差分）、:25-96（SUMQ 直接计数：能量格宽 DEL=1.0
!       :40、+20 格余量 :50-51、各模有效整数频 NEFIN=int(W/DEL+0.499999)
!       :57-59、T 数组加法序列 T(j+ne_i)+=T(j) :63-77、累计化 T(i)+=T(i−1)
!       :81-83、E/DEL 线性内插 :86-89、负探针早退 :44-47、T 分配/释放 :53/:93）；
!       src_VENUS/EIGOUT.f:1-39（本征输出：上下三角互转取行向量 :13-19、
!       6 列块写本征值 EIG 与逐行向量分量 :21-31、写后转置复原 :33-38）
!=====================================================================
module density
   implicit none
   private
   public :: dens_states, eig_out
contains
   !------------------------------------------------------------------
   ! dens_states(w_mode, e_probe, rho) — 态密度（Beyer-Swinehart 直接计数 +
   !                                      中心差分）
   ! 参考: DENQ.f:10-18（差分壳）；DENQ.f:25-96（SUMQ 直接计数——内联为步 1）
   subroutine dens_states(w_mode, e_probe, rho)
      real(8), intent(in) :: w_mode(:)  ! 简正模频率表 [能量格]（旧 WK——以格宽
                                        ! DEL 为单位的有效频率；cm⁻¹→格 折算经
                                        ! consts 于装配相位完成——C 阶段定形）
      real(8), intent(in) :: e_probe    ! 探针能量 [能量格]（旧 EN）
      real(8), intent(out) :: rho       ! 态密度 dN/dE [1/格]（旧 RHO）
      ! 伪代码：
      ! 1. 累计态数 N(E)（Beyer-Swinehart 直接计数，旧 SUMQ :25-96）：
      ! 1a. n_grid = int(e_probe/DEL + 0.499999) + 20（DEL=1.0 格；+20 格
      !     余量，旧 :40-51；e_probe<0 → N=0 早退，旧 :44-47）
      ! 1b. 逐模有效整数频 ne_i = int(w_mode(i)/DEL + 0.499999)（旧 :57-59）
      ! 1c. 工作数组 t(1..n_grid)：t(1)=1、其余清零（旧 :63-66）
      ! 1d. 逐模加法序列：t(j+ne_i) += t(j)，j = 1..n_grid−ne_i
      !     （逐模单遍前向累加——Beyer-Swinehart 核心，旧 :70-77）
      ! 1e. 累计化：t(i) += t(i−1)（逐格前缀和——N(E) 化，旧 :81-83）
      ! 1f. N(E) = t 于 e_probe/DEL 的线性内插（旧 :86-89）；
      !     t 释放（旧 :93）
      ! 2. rho = (N(e_probe+50) − N(e_probe−50))/100（±50 格中心差分除
      !     100 格窗宽，旧 DENQ.f:10-18；步 1 对两探针各执行一次）
   end subroutine dens_states

   !------------------------------------------------------------------
   ! eig_out(a, eig, unit) — 本征值/向量格式化输出（6 列块；写后复原）
   ! 参考: EIGOUT.f:13-39
   subroutine eig_out(a, eig, unit)
      real(8), intent(inout) :: a(:,:)  ! 本征向量阵（入口：向量按列——linalg
                                         ! eig_sym 约定；写前临时转置改按行，
                                         ! 写后复原原地）
      real(8), intent(in) :: eig(:)     ! 本征值（升序，与 a 列对齐；旧全局
                                         ! EIG 改哑元传入）
      integer, intent(in) :: unit       ! 输出单元号 [—]（旧哑元 IP）
      ! 伪代码：
      ! 1. a 上下三角互换（列存向量临时改行存——旧 :13-19）
      ! 2. 6 列块循环（旧 :21-31）：块头列号行；本征值行（eig(块列)）；
      !    逐行 i 写行号 + 向量分量 a(i, 块列)
      ! 3. a 再次上下三角互换复原（旧 :33-38）
   end subroutine eig_out
end module density
