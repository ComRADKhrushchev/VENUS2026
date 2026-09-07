!=====================================================================
! elec_adiabatic.f90 — L3 乙类退化成员：绝热电子律（ds/dt≡0、力规则恒等——
!                      单面表示下复合塌缩为纯 PES，接口在退化点良定义）
! 参考: docs/design_spec_completeness.md §3.2-2（两独立轴/退化点良定义/
!       绝热：F=−∇V_α——消费单面表示时复合塌缩为纯 PES）；
!       docs/final_completion_plan.md §5 N1（乙类成员退化点——显式存在约 20 行）
! 注册：facade_reg_elec(name, proc) 由装配相位执行（driver）；本退化成员两
! 例程到乙类接口槽的适配接线 C2 定形（E9 乙类签名扩 dt、E10 再扩 ws——
! elec_proc_i(s, f, dt, ws)；适配器忽略 dt 与 ws：绝热律无跃迁判定、
! 不消费多面构件——退化点对整条通道恒等）
!=====================================================================
module elec_adiabatic
   use state, only: elec_t
   implicit none
   private
   public :: ad_dsdt, ad_force_rule
contains
   !------------------------------------------------------------------
   ! ad_dsdt(s) — 电子态传播规则：绝热退化点 ds/dt ≡ 0（电子态冻结当前面）
   ! 参考: docs/design_spec_completeness.md §3.2-2（绝热是退化点）
   subroutine ad_dsdt(s)
      type(elec_t), intent(inout) :: s   ! 电子态（绝热律下不变——单面无耦合无跃迁）
      ! 伪代码：
      ! 1. 恒等退化：s 保持当前面 i_surf 不变（ds/dt ≡ 0——传播子对 s 为空推进）
   end subroutine ad_dsdt

   !------------------------------------------------------------------
   ! ad_force_rule(f) — 电子律力修饰：绝热退化点恒等（Σ单面势力即最终核力）
   ! 参考: docs/design_spec_completeness.md §3.2-2（绝热：F=−∇V_α）
   subroutine ad_force_rule(f)
      real(8), intent(inout) :: f(:)   ! 核力 [内部能量/Å]（进出恒等）
      ! 伪代码：
      ! 1. 恒等退化：f 原样通过（无面权重平均、无面切换修饰）
   end subroutine ad_force_rule
end module elec_adiabatic
