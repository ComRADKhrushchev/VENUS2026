!=====================================================================
! cau_peri_pot.f90 - C@Au(111) periodic 2D-analytic ground PES +
!   Hubbard delivery (pure physics layer; q in, eV out)
! Design:
!   Three-atom model of the vendored surface_potentials engine
!   (quarantined; see its header): atom 1 = surface atom S, atom 2 =
!   ghost oscillator Y, atom 3 = incident C. The C-S pair rides the
!   (u,v,z) Rydberg/Morse analytic surface (Fourier interpolation over
!   top/bridge/hollow parameters from READAP defaults or the NN path);
!   the S-Y surface oscillator and the Y self-harmonic close the model.
!   Hubbard chain: J = E_1D - E_3P, U = E_mC - E_1D - E_3P + E_pC,
!   h_til = (E_1D - E_pC) - U from the four state Morse surfaces,
!   with the CRYSTAL-FIELD SPLIT h(1/2) = h_til +/- delta_cf(u,v,z);
!   Vcp from the coupling surface. Kernels chain the analytic state
!   gradients onto atoms 3/1 (+/- d/dr_rel) and the crystal-field
!   gradient into the split. All eV / A.
!=====================================================================
module cau_peri_pot
   use surface_potentials
   implicit none
   private
   public :: cau_peri_load, cau_peri_ground, cau_peri_hubbard, cau_peri_hubbard_val

   logical :: loaded = .false.

contains

   !------------------------------------------------------------------
   ! warp_z(z, z_eff, dz) - the far-field evaluation-height warp (one
   !                 definition, shared by state_val and the crystal-field
   !                 call): at/above the 6 A freeze line every sheet is
   !                 evaluated at exactly z_eff = 6 A (all z gradients zero,
   !                 values hold their asymptote); C1 smoothstep ramp on
   !                 [5, 6]; below 5 A z_eff = z, dz = 1 (bit-identical)
   !------------------------------------------------------------------
   subroutine warp_z(z, z_eff, dz)
      real(8), intent(in)  :: z
      real(8), intent(out) :: z_eff, dz
      real(8) :: t, sw
      real(8), parameter :: z_freeze = 6.0d0, z_knee = 5.0d0
      z_eff = z
      dz = 1.0d0
      if (z > z_knee) then
         t = min((z - z_knee)/(z_freeze - z_knee), 1.0d0)
         sw = 3.0d0*t*t - 2.0d0*t*t*t
         z_eff = z - sw*(z - z_freeze)
         dz = 1.0d0 - sw - 6.0d0*t*(1.0d0 - t)/(z_freeze - z_knee)*(z - z_freeze)
      end if
   end subroutine warp_z

   !------------------------------------------------------------------
   ! cau_peri_load(nn_weights_file) - READAP defaults (analytic mode)
   !                 or the NN weight file (use_nn_pes)
   !------------------------------------------------------------------
   subroutine cau_peri_load(nn_file)
      character(len=*), intent(in), optional :: nn_file
      integer :: ierr
      if (present(nn_file)) then
         call load_all_nn_weights(nn_file, ierr)
         if (ierr /= 0) then
            call stop_pot('cau_peri_load', 'cannot load the NN weight file "'//trim(nn_file)//'"')
         end if
      else
         call readap()
      end if
      loaded = .true.
   end subroutine cau_peri_load

   !------------------------------------------------------------------
   ! state_val(ipot_r, r, e, de) - one analytic state value + gradient
   !                 (NN path when use_nn_pes, else the Morse/Rydberg
   !                 family; ipot_r selects the parm set)
   !------------------------------------------------------------------
   subroutine state_val(ipot_r, r, e, de)
      integer, intent(in) :: ipot_r
      real(8), intent(in) :: r(3)
      real(8), intent(out) :: e, de(3)
      real(8) :: z_eff, dz, t, sw
      ! far-field freeze (ruling 2026-09-29): the vendored sheets' extrapolation
      ! beyond ~5.5 A is unphysical (NN tanh tail, analytic far-gradient flip);
      ! the evaluation height is warped to z_eff so that at/above the 6 A
      ! analysis line every sheet is frozen at its 6 A geometry - all delivered
      ! gradients are exactly zero there, values hold their asymptotes. C1
      ! smoothstep ramp on [5, 6]; below 5 A the call is untouched (z_eff = z,
      ! dz = 1). Full record: docs/plans/2026-09-28-coupling-kernel-fix.md.
      real(8), parameter :: z_freeze = 6.0d0, z_knee = 5.0d0
      z_eff = r(3)
      dz = 1.0d0
      if (r(3) > z_knee) then
         t = min((r(3) - z_knee)/(z_freeze - z_knee), 1.0d0)
         sw = 3.0d0*t*t - 2.0d0*t*t*t
         z_eff = r(3) - sw*(r(3) - z_freeze)
         dz = 1.0d0 - sw - 6.0d0*t*(1.0d0 - t)/(z_freeze - z_knee)*(r(3) - z_freeze)
      end if
      if (use_nn_pes) then
         ! the analytic branch encodes the sheets as 0=ground,1=1D,2=3P,
         ! 3=mC,4=pC,5=coupling; the NN branch numbers the SAME sheets
         ! NN_E0=1..NN_COUP=6 - the +1 map. (A partial 0->NN_E0 merge once
         ! left 1..5 reading one sheet off: every Hubbard channel wrong by
         ! one network, "Vcp" was the pC surface - multi-eV, non-decaying,
         ! charging the impurity through the whole approach band - 2026-09-29)
         if (ipot_r < 0 .or. ipot_r > 5) then
            call stop_pot('state_val', 'internal state index out of range')
         end if
         call calc_nn_and_grad(ipot_r + 1, r(1), r(2), z_eff, e, de(1), de(2), de(3))
         de(3) = de(3)*dz
      else if (ipot_r == 0) then
         call calc_pot_spec_and_grad(parm_nas, r(1), r(2), z_eff, e, de(1), de(2), de(3))
         de(3) = de(3)*dz
      else if (ipot_r == 5) then
         call calc_coup_and_grad(parm_coup, r(1), r(2), z_eff, e, de(1), de(2), de(3))
         de(3) = de(3)*dz
      else if (ipot_r == 1) then
         call calc_pot0_and_grad(parm_1d, r(1), r(2), z_eff, e, de(1), de(2), de(3))
         de(3) = de(3)*dz
      else if (ipot_r == 2) then
         call calc_pot0_and_grad(parm_3p, r(1), r(2), z_eff, e, de(1), de(2), de(3))
         de(3) = de(3)*dz
      else if (ipot_r == 3) then
         call calc_pot0_and_grad(parm_mc, r(1), r(2), z_eff, e, de(1), de(2), de(3))
         de(3) = de(3)*dz
      else if (ipot_r == 4) then
         call calc_pot0_and_grad(parm_pC, r(1), r(2), z_eff, e, de(1), de(2), de(3))
         de(3) = de(3)*dz
      else
         call stop_pot('state_val', 'internal state index out of range')
      end if
   end subroutine state_val

   !------------------------------------------------------------------
   ! cau_peri_ground(q, e0, g0) - pair surface + oscillator harmonics
   !------------------------------------------------------------------
   subroutine cau_peri_ground(q, e0, g0)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: e0
      real(8), intent(out) :: g0(:)
      real(8) :: r(3), e, de(3), dsy(3), dy(3)
      integer :: k
      if (.not. loaded) then
         call stop_pot('cau_peri_ground', 'the engine is not initialized (seam order error)')
      end if
      if (size(q) /= 9) then
         call stop_pot('cau_peri_ground', 'the periodic model holds exactly 3 atoms (9 dofs)')
      end if
      r = q(7:9) - q(1:3)
      call state_val(0, r, e, de)
      e0 = e + 0.5d0*sum(k_surf(:)*(q(4:6) - q(1:3))**2) &
              + 0.5d0*sum(k_ghost(:)*(q(4:6) - Q0_ghost(:))**2)
      ! kernels: pair on C(+)/S(-); oscillators on S/Y
      g0 = 0.0d0
      g0(7:9) = g0(7:9) + de
      g0(1:3) = g0(1:3) - de
      g0(1:3) = g0(1:3) - k_surf(:)*(q(4:6) - q(1:3))
      g0(4:6) = g0(4:6) + k_surf(:)*(q(4:6) - q(1:3)) + k_ghost(:)*(q(4:6) - Q0_ghost(:))
   end subroutine cau_peri_ground

   !------------------------------------------------------------------
   ! hub_chain(r, h_til, u, j, vcp, dh_t, du, dj, dv) - the state
   !                 combination chains (values + r_rel kernels)
   !------------------------------------------------------------------
   subroutine hub_chain(r, h_til, u, j, vcp, dh_t, du, dj, dv)
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: h_til, u, j, vcp
      real(8), intent(out) :: dh_t(3), du(3), dj(3), dv(3)
      real(8) :: e1, e3, em, ep, ec
      real(8) :: d1(3), d3(3), dm(3), dp(3), dc(3)
      call state_val(1, r, e1, d1)
      call state_val(2, r, e3, d3)
      call state_val(3, r, em, dm)
      call state_val(4, r, ep, dp)
      call state_val(5, r, ec, dc)
      j = e1 - e3
      u = em - e1 - e3 + ep
      h_til = (e1 - ep) - u
      vcp = ec
      dj = d1 - d3
      du = dm - d1 - d3 + dp
      dh_t = (d1 - dp) - du
      dv = dc
   end subroutine hub_chain

   !------------------------------------------------------------------
   ! cau_peri_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp) - delivery
   !                 with the crystal-field split and kernels; the far
   !                 field beyond the 6 A line: h/U/J hold their frozen
   !                 sheet values with zero z gradients, Vcp DECAYS TO
   !                 ZERO (the physical asymptote - a frozen nonzero
   !                 coupling keeps the impurity-bath hybridization and
   !                 the charge churn alive at any distance)
   !------------------------------------------------------------------
   subroutine cau_peri_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: h(2), u, j, vcp
      real(8), intent(out) :: dh(2, 9), du(9), dj(9), dvcp(9)
      real(8) :: r(3), h_t, dh_t(3), d_u(3), d_j(3), d_v(3)
      real(8) :: cf, dcf(3), z_eff, dz, t, sw, dsw
      real(8), parameter :: z_freeze = 6.0d0, z_knee = 5.0d0
      integer :: k
      if (.not. loaded) then
         call stop_pot('cau_peri_hubbard', 'the engine is not initialized (seam order error)')
      end if
      r = q(7:9) - q(1:3)
      call hub_chain(r, h_t, u, j, vcp, dh_t, d_u, d_j, d_v)
      ! crystal field on the warped height (the frozen-zone cf is the 6 A
      ! sheet value; its z kernel rides dz)
      call warp_z(r(3), z_eff, dz)
      call calc_crystal_field(r(1), r(2), z_eff, cf, dcf(1), dcf(2), dcf(3))
      dcf(3) = dcf(3)*dz
      ! Vcp far-field decay: v_used = v(z_eff)*(1-S(z)); beyond the freeze
      ! line the coupling is exactly zero (gradient included). d_v(3) arrives
      ! from hub_chain already carrying the warp chain rule (dz factor)
      if (r(3) > z_knee) then
         t = min((r(3) - z_knee)/(z_freeze - z_knee), 1.0d0)
         sw = 3.0d0*t*t - 2.0d0*t*t*t
         dsw = 6.0d0*t*(1.0d0 - t)/(z_freeze - z_knee)
         d_v(3) = d_v(3)*(1.0d0 - sw) - vcp*dsw
         d_v(1) = d_v(1)*(1.0d0 - sw)
         d_v(2) = d_v(2)*(1.0d0 - sw)
         vcp = vcp*(1.0d0 - sw)
      end if
      h(1) = h_t + cf
      h(2) = h_t - cf
      ! kernels: chains onto C(+)/S(-) with the split sign on cf
      dh = 0.0d0
      du = 0.0d0
      dj = 0.0d0
      dvcp = 0.0d0
      do k = 1, 3
         dh(1, k + 6) = dh_t(k) + dcf(k)
         dh(1, k) = -(dh_t(k) + dcf(k))
         dh(2, k + 6) = dh_t(k) - dcf(k)
         dh(2, k) = -(dh_t(k) - dcf(k))
         du(k + 6) = d_u(k)
         du(k) = -d_u(k)
         dj(k + 6) = d_j(k)
         dj(k) = -d_j(k)
         dvcp(k + 6) = d_v(k)
         dvcp(k) = -d_v(k)
      end do
   end subroutine cau_peri_hubbard

   !------------------------------------------------------------------
   ! cau_peri_hubbard_val(q, h, u, j, vcp) - values-only delivery
   !------------------------------------------------------------------
   subroutine cau_peri_hubbard_val(q, h, u, j, vcp)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: h(2), u, j, vcp
      real(8) :: r(3), h_t, dh_t(3), d_u(3), d_j(3), d_v(3), cf, dcf(3)
      real(8) :: z_eff, dz, t, sw
      real(8), parameter :: z_freeze = 6.0d0, z_knee = 5.0d0
      r = q(7:9) - q(1:3)
      call hub_chain(r, h_t, u, j, vcp, dh_t, d_u, d_j, d_v)
      call warp_z(r(3), z_eff, dz)
      call calc_crystal_field(r(1), r(2), z_eff, cf, dcf(1), dcf(2), dcf(3))
      ! the same Vcp far-field decay as the kernel delivery
      if (r(3) > z_knee) then
         t = min((r(3) - z_knee)/(z_freeze - z_knee), 1.0d0)
         sw = 3.0d0*t*t - 2.0d0*t*t*t
         vcp = vcp*(1.0d0 - sw)
      end if
      h(1) = h_t + cf
      h(2) = h_t - cf
   end subroutine cau_peri_hubbard_val

   ! stop_pot(who, msg) - the named-abort channel (stderr + STOP 1)
   subroutine stop_pot(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_peri_pot error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container physics layer is not usable)'
      stop 1
   end subroutine stop_pot
end module cau_peri_pot
