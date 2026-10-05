!=====================================================================
! cau_tab_pot.f90 - the table-PES container's physics layer: the
!   DFT-point sheets + oscillator harmonics and the Hubbard delivery
!   (pure physics; q in, eV out)
! Design:
!   Three-atom model identical to cau_peri (atom 1 = surface atom S,
!   atom 2 = ghost oscillator Y, atom 3 = incident C; the C-S pair
!   rides r = q_C - q_S). CHANNEL-DIRECT delivery (2026-10-01 user
!   ruling): the five channels (E0 / h / U / J / VCP) are each read
!   from their own table sheet - per-site 1D Akima at the warped
!   height + the three-site Fourier lateral combination (values ->
!   lateral gradient, site slopes -> dz gradient; the period rides the
!   LATTICE row). The former four-diabat subtraction (J = E_1D-E_3P,
!   U, h_til) is RETIRED: it stacked each sheet's interpolation error
!   and C1 kinks into the channels; the GS DFT data is reference-only
!   and the mechanical base is the E0 sheet. The crystal-field split
!   h(1/2) = h +/- delta_cf stays analytic (kcal-parameter convention,
!   the 2026-09-28 unit ruling). All eV / A.
!   Ported from cau_peri_pot.f90 (2026-09-30): warp_z freeze line, the
!   Vcp far-field decay, the kernel assembly.
!=====================================================================
module cau_tab_pot
   use consts,       only: pi
   use cau_tab_tab
   implicit none
   private
   public :: cau_tab_load, cau_tab_ground, cau_tab_hubbard, cau_tab_hubbard_val
   public :: k_surf, k_ghost, Q0_ghost           ! oscillator constants (the box assigns)

   real(8), parameter :: sq3 = 1.7320508075688772d0
   real(8) :: k_surf(3) = [1.0d0, 1.0d0, 1.0d0]     ! S-Y pair spring [eV/A^2]
   real(8) :: k_ghost(3) = [1.0d0, 1.0d0, 1.0d0]    ! ghost self spring [eV/A^2]
   real(8) :: Q0_ghost(3) = [0.0d0, 0.0d0, 0.0d0]   ! ghost equilibrium [A]
   logical :: loaded = .false.

contains

   !------------------------------------------------------------------
   ! eff_par(u, v, p_top, p_bri, p_cen, p_m, dp_dx, dp_dy) - the
   !                 three-site Fourier lateral combination + Cartesian
   !                 gradient (the effective_par_grad form; the period
   !                 rides the table's LATTICE row)
   !------------------------------------------------------------------
   subroutine eff_par(u, v, p_top, p_bri, p_cen, p_m, dp_dx, dp_dy)
      real(8), intent(in)  :: p_top, p_bri, p_cen
      real(8), intent(in)  :: u, v
      real(8), intent(out) :: p_m, dp_dx, dp_dy
      real(8) :: c0, c1, c2
      real(8) :: arg_u, arg_v, arg_uv
      real(8) :: fir, snd
      real(8) :: dfir_du, dfir_dv, dsnd_du, dsnd_dv
      real(8) :: cons
      cons = tab_lattice()
      ! 1. Fourier coefficients from the three site values
      c1 = (p_top - p_bri)/4.0d0
      c0 = (p_top + 2.0d0*p_cen)/3.0d0
      c2 = (p_top - 3.0d0*c1 - c0)/3.0d0
      ! 2. basis at (u,v) + first/second harmonic derivatives
      arg_u = 2.0d0*pi*u/cons
      arg_v = 2.0d0*pi*v/cons
      arg_uv = 2.0d0*pi*(u + v)/cons
      fir = cos(arg_u) + cos(arg_v) + cos(arg_uv)
      snd = cos(2.0d0*arg_u) + cos(2.0d0*arg_v) + cos(2.0d0*arg_uv)
      p_m = c0 + c1*fir + c2*snd
      dfir_du = -(2.0d0*pi/cons)*(sin(arg_u) + sin(arg_uv))
      dfir_dv = -(2.0d0*pi/cons)*(sin(arg_v) + sin(arg_uv))
      dsnd_du = -(4.0d0*pi/cons)*(sin(2.0d0*arg_u) + sin(2.0d0*arg_uv))
      dsnd_dv = -(4.0d0*pi/cons)*(sin(2.0d0*arg_v) + sin(2.0d0*arg_uv))
      ! 3. u = x - y/sqrt(3), v = 2y/sqrt(3) -> du/dx = 1, dv/dx = 0,
      !    du/dy = -1/sqrt(3), dv/dy = 2/sqrt(3)
      dp_dx = c1*dfir_du + c2*dsnd_du
      dp_dy = c1*(-dfir_du/sq3 + 2.0d0*dfir_dv/sq3) + c2*(-dsnd_du/sq3 + 2.0d0*dsnd_dv/sq3)
   end subroutine eff_par

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
   ! cau_tab_load(table_file) - load the DFT-point table (the box hands
   !                 the TABLE_FILE value; re-loadable per tab_load)
   !------------------------------------------------------------------
   subroutine cau_tab_load(tab_file)
      character(len=*), intent(in) :: tab_file
      call tab_load(tab_file)
      loaded = .true.
   end subroutine cau_tab_load

   !------------------------------------------------------------------
   ! state_val(isheet, r, e, de) - one sheet value + gradient at the
   !                 pair coordinate r = q_C - q_S: per-site Akima at
   !                 the warped height, Fourier lateral combination of
   !                 the site values (lateral gradient) and of the site
   !                 slopes (dz gradient, chained with the warp scale)
   !------------------------------------------------------------------
   subroutine state_val(isheet, r, e, de)
      integer, intent(in) :: isheet
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: e, de(3)
      real(8) :: z_eff, dz, u, v, dd
      real(8) :: vt, vb, vh, dt_, db, dh
      call warp_z(r(3), z_eff, dz)
      u = r(1) - r(2)/sq3
      v = 2.0d0*r(2)/sq3
      call tab_eval(isheet, si_top, z_eff, vt, dt_)
      call tab_eval(isheet, si_bri, z_eff, vb, db)
      call tab_eval(isheet, si_hol, z_eff, vh, dh)
      call eff_par(u, v, vt, vb, vh, e, de(1), de(2))
      call eff_par(u, v, dt_, db, dh, de(3), dd, dd)
      de(3) = de(3)*dz
   end subroutine state_val

   !------------------------------------------------------------------
   ! cau_tab_ground(q, e0, g0) - GS sheet + oscillator harmonics
   !------------------------------------------------------------------
   subroutine cau_tab_ground(q, e0, g0)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: e0
      real(8), intent(out) :: g0(:)
      real(8) :: r(3), e, de(3)
      if (.not. loaded) then
         call stop_pot('cau_tab_ground', 'the engine is not initialized (seam order error)')
      end if
      if (size(q) /= 9) then
         call stop_pot('cau_tab_ground', 'the periodic model holds exactly 3 atoms (9 dofs)')
      end if
      r = q(7:9) - q(1:3)
      call state_val(sh_E0, r, e, de)
      e0 = e + 0.5d0*sum(k_surf(:)*(q(4:6) - q(1:3))**2) &
              + 0.5d0*sum(k_ghost(:)*(q(4:6) - Q0_ghost(:))**2)
      ! kernels: pair on C(+)/S(-); oscillators on S/Y
      g0 = 0.0d0
      g0(7:9) = g0(7:9) + de
      g0(1:3) = g0(1:3) - de
      g0(1:3) = g0(1:3) - k_surf(:)*(q(4:6) - q(1:3))
      g0(4:6) = g0(4:6) + k_surf(:)*(q(4:6) - q(1:3)) + k_ghost(:)*(q(4:6) - Q0_ghost(:))
   end subroutine cau_tab_ground

   !------------------------------------------------------------------
   ! channel_eval(r, h_til, u, j, vcp, dh_t, du, dj, dv) - the four
   !                 electronic channels, each its OWN sheet evaluation
   !                 (value + r_rel kernel; no cross-sheet arithmetic -
   !                 the 2026-10-01 channel-direct ruling)
   !------------------------------------------------------------------
   subroutine channel_eval(r, h_til, uu, j, vcp, dh_t, du, dj, dv)
      real(8), intent(in)  :: r(3)
      real(8), intent(out) :: h_til, uu, j, vcp
      real(8), intent(out) :: dh_t(3), du(3), dj(3), dv(3)
      call state_val(sh_h, r, h_til, dh_t)
      call state_val(sh_U, r, uu, du)
      call state_val(sh_J, r, j, dj)
      call state_val(sh_VCP, r, vcp, dv)
   end subroutine channel_eval

   !------------------------------------------------------------------
   ! cau_tab_hubbard(q, h, u, j, vcp, dh, du, dj, dvcp) - delivery
   !                 with the crystal-field split and kernels; the far
   !                 field beyond the 6 A line: h/U/J hold their frozen
   !                 sheet values with zero z gradients, Vcp DECAYS TO
   !                 ZERO (the physical asymptote)
   !------------------------------------------------------------------
   subroutine cau_tab_hubbard(q, h, uu, j, vcp, dh, du, dj, dvcp)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: h(2), uu, j, vcp
      real(8), intent(out) :: dh(2, 9), du(9), dj(9), dvcp(9)
      real(8) :: r(3), h_t, dh_t(3), d_u(3), d_j(3), d_v(3)
      real(8) :: cf, dcf(3), z_eff, dz, t, sw, dsw
      real(8), parameter :: z_freeze = 6.0d0, z_knee = 5.0d0
      integer :: k
      if (.not. loaded) then
         call stop_pot('cau_tab_hubbard', 'the engine is not initialized (seam order error)')
      end if
      r = q(7:9) - q(1:3)
      call channel_eval(r, h_t, uu, j, vcp, dh_t, d_u, d_j, d_v)
      ! crystal field on the warped height (the frozen-zone cf is the 6 A
      ! sheet value; its z kernel rides dz)
      call warp_z(r(3), z_eff, dz)
      call calc_crystal_field(r(1), r(2), z_eff, cf, dcf(1), dcf(2), dcf(3))
      dcf(3) = dcf(3)*dz
      ! Vcp far-field decay: v_used = v(z_eff)*(1-S(z)); beyond the freeze
      ! line the coupling is exactly zero (gradient included). d_v(3) arrives
      ! from channel_eval already carrying the warp chain rule (dz factor)
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
   end subroutine cau_tab_hubbard

   !------------------------------------------------------------------
   ! cau_tab_hubbard_val(q, h, u, j, vcp) - values-only delivery
   !------------------------------------------------------------------
   subroutine cau_tab_hubbard_val(q, h, uu, j, vcp)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: h(2), uu, j, vcp
      real(8) :: r(3), h_t, dh_t(3), d_u(3), d_j(3), d_v(3), cf, dcf(3)
      real(8) :: z_eff, dz, t, sw
      real(8), parameter :: z_freeze = 6.0d0, z_knee = 5.0d0
      r = q(7:9) - q(1:3)
      call channel_eval(r, h_t, uu, j, vcp, dh_t, d_u, d_j, d_v)
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
   end subroutine cau_tab_hubbard_val

   !------------------------------------------------------------------
   ! calc_crystal_field(x, y, z, delta_cf, ddelta_dx, ddelta_dy,
   !                 ddelta_dz) - the analytic 2px/2py site split with
   !                 the tanh z-decay, kcal-parameter convention, eV
   !                 delivery (the 2026-09-28 unit ruling: the kcal
   !                 parameters are converted by the legislated pair's
   !                 inverse at delivery). Verbatim values from
   !                 surface_potentials (cau_peri); the lateral period
   !                 rides the table's LATTICE row
   !------------------------------------------------------------------
   subroutine calc_crystal_field(x, y, z, delta_cf, ddelta_dx, ddelta_dy, ddelta_dz)
      real(8), intent(in)  :: x, y, z
      real(8), intent(out) :: delta_cf, ddelta_dx, ddelta_dy, ddelta_dz
      real(8), parameter :: cf_top = 0.0d0
      real(8), parameter :: cf_bri = 5.0d0
      real(8), parameter :: cf_hcp = 0.0d0
      real(8), parameter :: cf_const = 0.01d0
      real(8), parameter :: z_switch = 3.0d0
      real(8), parameter :: z_width = 1.0d0
      real(8) :: u_loc, v_loc, delta_2d, d2d_dx, d2d_dy
      real(8) :: fz, dfz_dz, tanh_arg, sech2
      real(8) :: kcal_to_ev
      u_loc = x - y/sq3
      v_loc = 2.0d0*y/sq3
      call eff_par(u_loc, v_loc, cf_top, cf_bri, cf_hcp, delta_2d, d2d_dx, d2d_dy)
      ! z-decay: fz = 0.5*(1 - tanh((z-z0)/w)), 1 near surface, 0 at large z
      tanh_arg = (z - z_switch)/z_width
      fz = 0.5d0*(1.0d0 - tanh(tanh_arg))
      sech2 = 1.0d0/cosh(tanh_arg)
      sech2 = sech2*sech2
      dfz_dz = -0.5d0/z_width*sech2
      kcal_to_ev = 1.0d0/23.0605d0
      delta_cf = (delta_2d + cf_const)*fz*kcal_to_ev
      ddelta_dx = d2d_dx*fz*kcal_to_ev
      ddelta_dy = d2d_dy*fz*kcal_to_ev
      ddelta_dz = (delta_2d + cf_const)*dfz_dz*kcal_to_ev
   end subroutine calc_crystal_field

   ! stop_pot(who, msg) - the named-abort channel (stderr + STOP 1)
   subroutine stop_pot(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': cau_tab_pot error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the container physics layer is not usable)'
      stop 1
   end subroutine stop_pot
end module cau_tab_pot
