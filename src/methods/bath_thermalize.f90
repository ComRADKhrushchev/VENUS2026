!=====================================================================
! bath_thermalize.f90 - bath species: the one-shot Maxwell-seed + free
!                       equilibration run (the THERMALIZE migration,
!                       2026-10-16 user ruling)
! Design:
!   The retired dist-registry member samp_thermalize was the ONLY scheme
!   word that consumed the PES and the whole-state geometry (its realize
!   ran prop_step x n_eq) - a bath-domain mechanism by engineering fact,
!   re-homed here as the third bath registration form (seed): the loop
!   head Maxwell-seeds the bath-domain atoms once, the loop body is a bare
!   free run (no in-loop coupling - the physical boundary the 2026-10-15
!   ruling recorded: a one-shot equilibration, not a reservoir), and
!   bath_run's tail exact drift removal closes it. Parameters ride the
!   BATH key family verbatim: BATH = THERMALIZE + N_BATH + DT_BATH +
!   T_BATH (zero new keys; the per-carrier multi-temperature capability
!   retired with the sampling identity - one T_BATH for the domain).
!   The seeding kernel (per-atom Gaussian over free components only +
!   the domain's mass-weighted COM-momentum projection) is the retired
!   member's draw core, per-atom over the substrate domain.
! Units: T [K]; sigma_p [amu*A/fs] = sqrt(m*kb*T); p [amu*A/fs]; dt [fs]
!=====================================================================
module bath_thermalize
   use consts,        only: kb_act
   use config,        only: reactants
   use config_atoms,  only: list_atoms
   use input,         only: input_bath_code, buffer_n_rows, buffer_key, buffer_val, buffer_line
   use state,         only: state_t
   use rng,           only: rng_gauss
   use seam_cast,     only: seam_real
   use bath,          only: bath_reg, bath_setup
   use diag,          only: diag_fatal, diag_cat_sampler
   implicit none
   private
   public :: member_register, member_init
   public :: thermalize_seed, thermalize_t_bath

   real(8) :: t_bath_c = -1.0d0     ! provisioned temperature [K] (negative = unset)
   logical :: inited = .false.
contains

   !------------------------------------------------------------------
   ! thermalize_seed(sta, mask) - the loop-head one-shot seeding: per-atom
   !                               Maxwell over the FREE components of the
   !                               masked atoms, then one exact mass-weighted
   !                               COM-momentum projection over the same
   !                               domain (the retired draw core, per-atom)
   !------------------------------------------------------------------
   subroutine thermalize_seed(sta, mask)
      type(state_t), intent(inout) :: sta
      logical, intent(in) :: mask(:)
      integer :: i, j
      real(8) :: m, sig, wt, p_tot(3)
      if (t_bath_c <= 0.0d0) then
         call stop_thx('thermalize_seed', 'no provisioned temperature (member_init '// &
            'must run first - T_BATH > 0 under BATH = THERMALIZE)')
      end if
      do i = 1, size(mask)
         if (.not. mask(i)) cycle
         m = list_atoms%mass(i)
         sig = sqrt(m*kb_act*t_bath_c)   ! kb_act: the action-context
                                             ! Boltzmann (p is amu*A/fs - the
                                             ! sigma^2 = m*kb_act*T face rule,
                                             ! consts.f90 face table)
         do j = 1, 3
            if (allocated(list_atoms%dof)) then
               if (i <= size(list_atoms%dof, 2)) then
                  if (.not. list_atoms%dof(j, i)) then
                     sta%p(3*i - 3 + j) = 0.0d0   ! constrained component: stays zero
                     cycle
                  end if
               end if
            end if
            sta%p(3*i - 3 + j) = sig*rng_gauss()
         end do
      end do
      ! one exact COM-momentum projection over the domain (mass-weighted):
      ! p_i <- p_i - (m_i/M)*P - the retired member's per-carrier drift
      ! removal, domain-wide (bath_run's tail projection then re-projects -
      ! idempotent, exact arithmetic both times)
      wt = 0.0d0
      p_tot = 0.0d0
      do i = 1, size(mask)
         if (.not. mask(i)) cycle
         wt = wt + list_atoms%mass(i)
         p_tot = p_tot + sta%p(3*i - 2:3*i)
      end do
      if (wt > 0.0d0) then
         do i = 1, size(mask)
            if (.not. mask(i)) cycle
            sta%p(3*i - 2:3*i) = sta%p(3*i - 2:3*i) - (list_atoms%mass(i)/wt)*p_tot
         end do
      end if
   end subroutine thermalize_seed

   !------------------------------------------------------------------
   ! member_register - the fixed-name weave entry (prebuild binds it into
   !                   registry_gen_all; the seed form: loop-head seeding,
   !                   bare-run loop body)
   !------------------------------------------------------------------
   subroutine member_register()
      call bath_reg('thermalize', seed=thermalize_seed)
   end subroutine member_register

   !------------------------------------------------------------------
   ! member_init - the fixed-name init seam: under BATH = THERMALIZE pull
   !               T_BATH (> 0 required - the seed temperature) and
   !               provision the loop shape through bath_setup
   !------------------------------------------------------------------
   subroutine member_init()
      character(len=512) :: sval
      integer :: i
      if (reactants%bath /= input_bath_code('THERMALIZE')) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('T_BATH'); call seam_real(buffer_key(i), sval, buffer_line(i), t_bath_c)
         end select
      end do
      if (t_bath_c <= 0.0d0) then
         call stop_thx('member_init', 'BATH = THERMALIZE requires T_BATH > 0 '// &
            '(the Maxwell-seed temperature), got a non-positive value')
      end if
      inited = .true.
      call bath_setup()
   end subroutine member_init

   pure real(8) function thermalize_t_bath()
      thermalize_t_bath = t_bath_c
   end function thermalize_t_bath

   subroutine stop_thx(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') 'bath_thermalize ('//trim(who)//'): '//trim(msg)
      call diag_fatal(diag_cat_sampler, 'bath_thermalize.f90', trim(who), trim(msg))
   end subroutine stop_thx

end module bath_thermalize
