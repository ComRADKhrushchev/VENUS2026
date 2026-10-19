!=====================================================================
! bath_nose_hoover.f90 - bath species: Nosé–Hoover two-chain thermostat
!   - a SELF-CONTAINED bath step (the step-takeover form: the chain must
!   interleave with the kick/drift beats, which an opaque prop_step
!   cannot offer, so the species carries its own Liouville VV splitting
!   - ANT manual p.71: NH cannot ride a plain/position-Verlet/Beeman
!   integrator, only the Liouville approach to VV). One bath step:
!   [chain exact ops dt/2] [p *= exp(-v1*dt/2)] [kick dt/2, drift, force,
!   kick dt/2] [p *= exp(-v1*dt/2)] [chain exact ops dt/2]; chain masses
!   follow the Nosé prescription Q_k = f_k*kb*T*(tau_nose/2*pi)^2 (ANT A2;
!   tau_nose = the chain oscillation period, chosen near the system's
!   typical vibrational period or the bath collision period)
! Design:
!   The chain itself is DETERMINISTIC (no random draws in the step); the
!   per-trajectory RESET seeds the chain from its CANONICAL equilibrium -
!   eta positions zeroed, chain momenta drawn p_eta_k ~ N(0, Q_k*kb*T) -
!   and the FIRST step of the trajectory seeds the masked atoms' momenta
!   Maxwell at T_BATH: a pure chain can only RESCALE momenta, so a rest
!   start (the ownerless seeding p = 0 at the equilibrium geometry) is a
!   fixed point with no injection mechanism (Andersen's velocity redraws
!   inject; NH does not) - the thermal seed is the standard NH
!   initialization, and both draws ride the bath loop's sanctioned
!   in-loop RNG exception. Equilibrium seeds additionally give every
!   trajectory an independent canonical start (a chain zeroed to rest
!   would land the ownerless atoms on the SAME endpoint every
!   trajectory - no sampling power). A roster with no ownerless atoms
!   (a fragment-only system under the freeze statute) draws NOTHING and
!   the step early-returns - the armed run stays bitwise a no-op, the
!   Andersen leg's freeze-twin contract. The whole STATE is advanced by
!   the step (frozen fragments included - they are bitwise restored by
!   the loop's unfreeze anyway, mirroring the sweep species' whole-state
!   prop_step); the thermostat and the seed touch the masked atoms only.
!   Units: t_bath [K]; tau_nose [fs]; G [amu*Ang^2/fs^2]; Q
!   [amu*Ang^2/fs^2]; p_eta [amu*Ang^2/fs]; dt [fs].
!=====================================================================
module bath_nose_hoover
   use consts,        only: kb_act, pi
   use rng,           only: rng_gauss
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use force_interface, only: force_eval
   use bath,          only: bath_reg, bath_setup
   use input,         only: buffer_n_rows, buffer_key, buffer_val, buffer_line, &
                            input_bath_code
   use seam_cast,     only: seam_real
   use diag,          only: diag_fatal, diag_cat_sampler
   implicit none
   private
   public :: nose_hoover_params_t, nose_hoover_init, nose_hoover_step, &
             nose_hoover_reset, nose_hoover_n_reset, nose_hoover_chain, &
             member_register, member_init

   ! Species-owned parameters
   type :: nose_hoover_params_t
      real(8) :: t_bath   = 0.0d0   ! bath temperature [K]
      real(8) :: tau_nose = 0.0d0   ! chain oscillation period [fs] (the Nosé
                                    ! prescription: the thermostat period near the
                                    ! system's vibrational period)
   end type nose_hoover_params_t

   ! Derived species state (module-private; rebuilt by every init; the chain
   ! variables are the STATEFUL part - re-seeded per trajectory by reset)
   integer :: nat_c = 0              ! atom count [-]
   integer :: nf_c = 0               ! bathed degrees of freedom 3*N_ownerless [-]
                                     ! (the roster's ownerless count, process-fixed)
   real(8) :: t_bath_c = 0.0d0       ! bath temperature [K]
   real(8) :: tau_nose_c = 0.0d0     ! chain oscillation period [fs]
   real(8) :: q1_c = 0.0d0           ! chain mass 1 [amu*Ang^2/fs^2]
   real(8) :: q2_c = 0.0d0           ! chain mass 2 [amu*Ang^2/fs^2]
   logical :: inited = .false.       ! init completed [flag]
   real(8) :: eta1 = 0.0d0           ! chain position 1 [-]
   real(8) :: eta2 = 0.0d0           ! chain position 2 [-]
   real(8) :: p_eta1 = 0.0d0         ! chain momentum 1 [amu*Ang^2/fs]
   real(8) :: p_eta2 = 0.0d0         ! chain momentum 2 [amu*Ang^2/fs]
   integer :: n_reset = 0            ! reset call count [-] (the criterion's
                                     ! once-per-trajectory witness)
   logical :: first_step = .true.    ! the trajectory's first step pending
                                     ! (the Maxwell seed of the masked momenta)
contains
   !------------------------------------------------------------------
   ! nose_hoover_init(p) - validate the bath parameters, derive the bathed
   !                       degrees of freedom (the roster's ownerless atoms)
   !                       and the chain masses, zero the chain
   !------------------------------------------------------------------
   subroutine nose_hoover_init(p)
      type(nose_hoover_params_t), intent(in) :: p ! species parameters
      integer :: i, fi, fa, n_own

      ! 1. parameter validation (fail-loud, naming field + value)
      if (p%t_bath <= 0.0d0) then
         call stop_nhz('nose_hoover_init', 'bath temperature t_bath [K] must be positive, got', &
            p%t_bath)
      end if
      if (p%tau_nose <= 0.0d0) then
         call stop_nhz('nose_hoover_init', 'chain period tau_nose [fs] must be positive '// &
            '(the Nosé prescription: near the system vibrational period), got', p%tau_nose)
      end if

      ! 2. atom-list presence (the step walks it) + the ownerless count
      !    (the bath scope: atoms no roster fragment claims)
      if (.not. allocated(list_atoms%mass)) then
         call stop_nhx('nose_hoover_init', 'no atom mass table - the atom list must be '// &
            'assembled before species init (assembly-order error)')
      end if
      nat_c = size(list_atoms%mass)
      n_own = nat_c
      if (allocated(list_atoms%frag)) then
         do i = 1, nat_c
            do fi = 1, size(list_atoms%frag)
               do fa = 1, list_atoms%frag(fi)%nat
                  if (list_atoms%frag(fi)%list(fa) == i) then
                     n_own = n_own - 1
                     exit
                  end if
               end do
            end do
         end do
      end if
      nf_c = 3*n_own
      t_bath_c = p%t_bath
      tau_nose_c = p%tau_nose

      ! 3. chain masses by the Nosé prescription (ANT A2:
      !    Q = f*kb*T*(tau/2pi)^2 - one chain variable per scale)
      q1_c = dble(max(nf_c, 1))*kb_act*t_bath_c*(tau_nose_c/(2.0d0*pi))**2
      q2_c = kb_act*t_bath_c*(tau_nose_c/(2.0d0*pi))**2
      eta1 = 0.0d0
      eta2 = 0.0d0
      p_eta1 = 0.0d0
      p_eta2 = 0.0d0
      first_step = .true.
      inited = .true.
   end subroutine nose_hoover_init

   !------------------------------------------------------------------
   ! nose_hoover_reset() - the per-trajectory re-seed (bath.f90's loop-head
   !                       hook): eta positions zeroed, chain momenta drawn
   !                       from their canonical equilibrium
   !                       p_eta_k ~ N(0, Q_k*kb*T), the first-step seed
   !                       armed. A roster with no ownerless atoms draws
   !                       NOTHING (the armed run stays bitwise a no-op -
   !                       the freeze-twin contract)
   !------------------------------------------------------------------
   subroutine nose_hoover_reset()
      n_reset = n_reset + 1
      first_step = .true.
      eta1 = 0.0d0
      eta2 = 0.0d0
      if (nf_c <= 0) return
      p_eta1 = sqrt(q1_c*kb_act*t_bath_c)*rng_gauss()
      p_eta2 = sqrt(q2_c*kb_act*t_bath_c)*rng_gauss()
   end subroutine nose_hoover_reset

   !------------------------------------------------------------------
   ! nose_hoover_step(sta, dt, mask) - ONE self-contained thermostat step:
   !                       [chain dt/2] [scale dt/2] [VV kick/drift/kick]
   !                       [scale dt/2] [chain dt/2]. The whole state is
   !                       advanced (the loop's unfreeze restores the frozen
   !                       fragments bitwise); the thermostat touches the
   !                       masked atoms only. Zero RNG (the first step's
   !                       Maxwell seed aside); empty scope returns untouched
   !------------------------------------------------------------------
   subroutine nose_hoover_step(sta, dt, mask)
      type(state_t), intent(inout) :: sta ! physical state (advanced one bath step)
      real(8), intent(in) :: dt           ! the bath step size [fs]
      logical, intent(in) :: mask(:)      ! per-atom thermostat scope [flag]
      real(8) :: scl, edum, k_m
      integer :: i, i3, n_m
      if (.not. inited) then
         call stop_nhx('nose_hoover_step', 'step called before nose_hoover_init - the '// &
            'species carries no parameter set')
      end if
      if (size(sta%p) < 3*nat_c .or. size(sta%q) < 3*nat_c) then
         call stop_nhx('nose_hoover_step', 'the state vector ends before the last '// &
            'atom - state_create must follow the atom list')
      end if
      if (size(mask) < nat_c) then
         call stop_nhx('nose_hoover_step', 'the thermostat scope mask ends before the last '// &
            'atom - the bath stage must build it over the whole atom list')
      end if
      n_m = count(mask(1:nat_c))
      if (n_m == 0) return                 ! empty scope: bitwise no-op
      if (3*n_m /= nf_c) then
         call stop_nhx('nose_hoover_step', 'the thermostat scope count disagrees with the '// &
            'init-time ownerless count (roster drift between init and the loop)')
      end if

      ! 0. the trajectory's first step seeds the masked momenta Maxwell at
      !    T_BATH (the chain can only rescale - a rest start is a fixed
      !    point with no injection; the thermal seed is the standard NH
      !    initialization, and the draw rides the in-loop exception)
      if (first_step) then
         do i = 1, nat_c
            if (.not. mask(i)) cycle
            i3 = 3*(i - 1)
            sta%p(i3+1) = sqrt(list_atoms%mass(i)*kb_act*t_bath_c)*rng_gauss()
            sta%p(i3+2) = sqrt(list_atoms%mass(i)*kb_act*t_bath_c)*rng_gauss()
            sta%p(i3+3) = sqrt(list_atoms%mass(i)*kb_act*t_bath_c)*rng_gauss()
         end do
         first_step = .false.
      end if

      ! 1. the first thermostat half: chain exact ops, then the scale factor
      call chain_ops(sta, mask, 0.5d0*dt)
      scl = exp(-(p_eta1/q1_c)*0.5d0*dt)
      do i = 1, nat_c
         if (.not. mask(i)) cycle
         i3 = 3*(i - 1)
         sta%p(i3+1) = sta%p(i3+1)*scl
         sta%p(i3+2) = sta%p(i3+2)*scl
         sta%p(i3+3) = sta%p(i3+3)*scl
      end do

      ! 2. the velocity-Verlet nucleus over the WHOLE state (kick dt/2,
      !    drift, force, kick dt/2 - the Lioulli splitting the chain needs)
      call force_eval(sta, edum)
      do i3 = 1, size(sta%p)
         sta%p(i3) = sta%p(i3) + 0.5d0*dt*sta%f(i3)
      end do
      do i3 = 1, size(sta%q)
         sta%q(i3) = sta%q(i3) + dt*sta%p(i3)/list_atoms%mass((i3+2)/3)
      end do
      call force_eval(sta, edum)
      do i3 = 1, size(sta%p)
         sta%p(i3) = sta%p(i3) + 0.5d0*dt*sta%f(i3)
      end do

      ! 3. the second thermostat half: the same chain-driven scale factor,
      !    then the chain exact ops at the new momenta
      do i = 1, nat_c
         if (.not. mask(i)) cycle
         i3 = 3*(i - 1)
         sta%p(i3+1) = sta%p(i3+1)*scl
         sta%p(i3+2) = sta%p(i3+2)*scl
         sta%p(i3+3) = sta%p(i3+3)*scl
      end do
      call chain_ops(sta, mask, 0.5d0*dt)
   end subroutine nose_hoover_step

   !------------------------------------------------------------------
   ! chain_ops(sta, mask, t) - the exactly-solved chain advance over t at a
   !                           fixed G1 (the masked kinetic energy): the
   !                           p_eta2 half impulses bracket the exact
   !                           p_eta1 drag solution; the eta positions ride
   !------------------------------------------------------------------
   subroutine chain_ops(sta, mask, t)
      type(state_t), intent(in) :: sta
      logical, intent(in) :: mask(:)
      real(8), intent(in) :: t
      real(8) :: k_m, g1, g2, xx, efac
      integer :: i, i3
      k_m = 0.0d0
      do i = 1, nat_c
         if (.not. mask(i)) cycle
         i3 = 3*(i - 1)
         k_m = k_m + (sta%p(i3+1)**2 + sta%p(i3+2)**2 + sta%p(i3+3)**2)/ &
               (2.0d0*list_atoms%mass(i))
      end do
      g1 = 2.0d0*k_m - dble(nf_c)*kb_act*t_bath_c
      g2 = p_eta1*p_eta1/q1_c - kb_act*t_bath_c
      p_eta2 = p_eta2 + 0.5d0*t*g2           ! half impulse at the old p1
      if (abs(p_eta2) > 1.0d-300) then
         xx = p_eta2*t/q2_c                  ! exact p1 under the p2 drag
         efac = exp(-xx)
         p_eta1 = p_eta1*efac + g1*(q2_c/p_eta2)*(1.0d0 - efac)
      else
         p_eta1 = p_eta1 + g1*t              ! p2 ~ 0: the drag vanishes
      end if
      g2 = p_eta1*p_eta1/q1_c - kb_act*t_bath_c
      p_eta2 = p_eta2 + 0.5d0*t*g2           ! half impulse at the new p1
      eta1 = eta1 + t*(p_eta1/q1_c)
      eta2 = eta2 + t*(p_eta2/q2_c)
   end subroutine chain_ops

   !------------------------------------------------------------------
   ! getters (copies of the derived species state - the criterion's anchors)
   !------------------------------------------------------------------
   pure integer function nose_hoover_n_reset()
      nose_hoover_n_reset = n_reset
   end function nose_hoover_n_reset

   pure subroutine nose_hoover_chain(e1, e2, pe1, pe2)
      real(8), intent(out) :: e1, e2, pe1, pe2
      e1 = eta1
      e2 = eta2
      pe1 = p_eta1
      pe2 = p_eta2
   end subroutine nose_hoover_chain

   !------------------------------------------------------------------
   ! private named-abort helpers (name the caller; value variants; abort via diag_fatal)
   subroutine stop_nhx(where, msg)
      character(len=*), intent(in) :: where, msg
      write (0, '(a)') trim(where)//': '//trim(msg)
      call diag_fatal(diag_cat_sampler, 'bath_nose_hoover.f90', trim(where), &
         'the bath loop is not usable')
   end subroutine stop_nhx

   subroutine stop_nhz(where, msg, val)
      character(len=*), intent(in) :: where, msg
      real(8), intent(in) :: val
      write (0, '(a,a,es12.4)') trim(where)//': ', trim(msg), val
      call diag_fatal(diag_cat_sampler, 'bath_nose_hoover.f90', trim(where), &
         'the bath loop is not usable')
   end subroutine stop_nhz

   !------------------------------------------------------------------
   ! member_register - the fixed-name weave entry (prebuild binds it into
   !                   registry_gen_all); the step-takeover form + the
   !                   per-trajectory re-seed ride the registration
   !------------------------------------------------------------------
   subroutine member_register()
      call bath_reg('nose_hoover', step=nose_hoover_step, reset=nose_hoover_reset)
   end subroutine member_register

   !------------------------------------------------------------------
   ! member_init - the fixed-name init seam (registry_gen_init calls it):
   !               an armed NOSE_HOOVER provision pulls T_BATH/TAU_NOSE
   !               and provisions the loop shape (the species carries its
   !               own Liouville VV splitting - no registered-propagator
   !               dependency to guard)
   !------------------------------------------------------------------
   subroutine member_init()
      type(nose_hoover_params_t) :: p   ! species-parameter receiver (field defaults)
      character(len=512) :: sval        ! one buffered raw value
      integer :: i
      if (reactants%bath /= input_bath_code('NOSE_HOOVER')) return
      do i = 1, buffer_n_rows()
         sval = buffer_val(i)
         select case (trim(buffer_key(i)))
         case ('T_BATH');   call seam_real(buffer_key(i), sval, buffer_line(i), p%t_bath)
         case ('TAU_NOSE'); call seam_real(buffer_key(i), sval, buffer_line(i), p%tau_nose)
         end select
      end do
      call nose_hoover_init(p)
      call bath_setup()
   end subroutine member_init

end module bath_nose_hoover
