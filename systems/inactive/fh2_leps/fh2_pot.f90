!=====================================================================
! fh2_pot.f90 - the F+H2 modified-LEPS potential (Takayanagi-Sato
!               surface): analytic V and gradient dV/dq in eV, pure
!               mathematics with zero program knowledge
! Design:
!   The extended-LEPS assembly is legislated verbatim from Muckerman,
!   JCP 56, 2997 (1972) Eqs. (1)(2a)(2b): with e_i = exp(-beta_i(r_i-r0_i)),
!   Q'_i = D_i/(4(1+Dl_i))*[(3+Dl_i)e^2 - (2+6*Dl_i)e],
!   A'_i = D_i/(4(1+Dl_i))*[(1+3Dl_i)e^2 - (6+2Dl_i)e],
!   V = sum(Q') - sqrt(A1^2+A2^2+A3^2-A1A2-A2A3-A3A1). The pair sum
!   Q'+A' = D(e^2-2e) exactly (a lone pair degenerates to plain Morse,
!   independent of the Sato parameter) - this invariant is the internal
!   fingerprint of the correct transcription, machine-checked by
!   check_fh2_leps together with four Muckerman saddle barriers (0.90/
!   0.89/1.43/2.55 kcal, agreement 0.02 kcal) and the T-S bent saddle
!   (0.18 kcal vs the published 8 meV). The T-S modification is the
!   angle-modulated Sato parameters of Takayanagi & Sato, CPL 144, 191
!   (1988): Dl_HH = Dl0_HH + a*sin^2(phi), Dl_HF likewise, phi = angle
!   between the H-H axis and the F-to-bond-midpoint vector. Legislated
!   constants (defaults, overridable via the box parameter channel):
!   Morse base [P8] Brown et al., JCP 82, 188 (1985) Table X - FH pair
!   De = 141.196 kcal/mol = 6.1228507621 eV, r0 = 0.917 A, beta =
!   2.2187 1/A; HH pair De = 109.449 kcal/mol = 4.7461676893 eV,
!   r0 = 0.7419 A, beta = 1.942 1/A (consistent with [P6] Table I);
!   Sato set Dl0_HH = 0.173, Dl0_FH = 0.165, a = 0.090 ([P2]).
!   Units: eV and Angstrom throughout; NO unit conversion happens here.
!=====================================================================
module fh2_pot
   implicit none
   private
   public :: fh2_vg, fh2_set_params, fh2_set_defaults

   ! Legislated defaults (the literature set above; the kcal->eV folding
   ! happened once at transcription time with the interface's legislated
   ! factor - a documentation fact, not a runtime conversion)
   real(8), parameter :: d_fh_d  = 6.1228507621d0   ! FH Morse De [eV] (141.196 kcal/mol)
   real(8), parameter :: r0_fh_d = 0.917d0          ! FH Morse re [A]
   real(8), parameter :: b_fh_d  = 2.2187d0         ! FH Morse beta [1/A]
   real(8), parameter :: d_hh_d  = 4.7461676893d0   ! HH Morse De [eV] (109.449 kcal/mol)
   real(8), parameter :: r0_hh_d = 0.7419d0         ! HH Morse re [A]
   real(8), parameter :: b_hh_d  = 1.942d0          ! HH Morse beta [1/A]
   real(8), parameter :: dl0_hh_d = 0.173d0         ! Sato base, HH pair [-] ([P2])
   real(8), parameter :: dl0_fh_d = 0.165d0         ! Sato base, FH pairs [-] ([P2])
   real(8), parameter :: amod_d   = 0.090d0         ! sin^2(phi) modulation [-] ([P2])

   real(8), save :: d_fh = d_fh_d, r0_fh = r0_fh_d, b_fh = b_fh_d
   real(8), save :: d_hh = d_hh_d, r0_hh = r0_hh_d, b_hh = b_hh_d
   real(8), save :: dl0_hh = dl0_hh_d, dl0_fh = dl0_fh_d, amod = amod_d

   real(8), parameter :: tiny_s = 1.0d-30    ! sqrt guard of the London form

contains
   !------------------------------------------------------------------
   ! fh2_set_defaults() - re-stage the legislated literature set (the
   !                 anchor battery of the criterion program switches
   !                 parameter sets and restores between them)
   !------------------------------------------------------------------
   subroutine fh2_set_defaults()
      d_fh = d_fh_d
      r0_fh = r0_fh_d
      b_fh = b_fh_d
      d_hh = d_hh_d
      r0_hh = r0_hh_d
      b_hh = b_hh_d
      dl0_hh = dl0_hh_d
      dl0_fh = dl0_fh_d
      amod = amod_d
   end subroutine fh2_set_defaults

   !------------------------------------------------------------------
   ! fh2_set_params(...) - the box parameter channel write: replace the
   !                 legislated defaults (validated by the caller)
   !------------------------------------------------------------------
   subroutine fh2_set_params(dfh, rfh, bfh, dhh, rhh, bhh, dlhh, dlfh, am)
      real(8), intent(in) :: dfh, rfh, bfh, dhh, rhh, bhh, dlhh, dlfh, am
      d_fh = dfh
      r0_fh = rfh
      b_fh = bfh
      d_hh = dhh
      r0_hh = rhh
      b_hh = bhh
      dl0_hh = dlhh
      dl0_fh = dlfh
      amod = am
   end subroutine fh2_set_params

   !------------------------------------------------------------------
   ! fh2_vg(natom, q, v, g) - the full surface: energy [eV] + gradient
   !                 dV/dq [eV/A] for F (atom 1) + H (atom 2) + H (atom 3);
   !                 the Sato parameters ride on sin^2(phi), so the
   !                 gradient carries the chain rule through the angle
   !------------------------------------------------------------------
   subroutine fh2_vg(natom, q, v, g)
      integer, intent(in)  :: natom
      real(8), intent(in)  :: q(3*natom)
      real(8), intent(out) :: v, g(3*natom)
      real(8) :: xf(3), xh1(3), xh2(3)          ! atom positions
      real(8) :: r1, r2, r3                     ! pair distances: F-H1, F-H2, H1-H2
      real(8) :: e1, e2, e3                     ! Morse exponentials
      real(8) :: q1, q2, q3, a1, a2, a3         ! branches Q' and A'
      real(8) :: dq1, dq2, dq3, da1, da2, da3   ! d/dr branches (fixed Dl)
      real(8) :: dqa1, dqa2, dqa3, daa1, daa2, daa3 ! d/dDl branches
      real(8) :: dl1, dl2, dl3                  ! Sato parameters (modulated)
      real(8) :: s, rs                            ! London form and its root
      real(8) :: dv1, dv2, dv3                    ! dV/dr_i at fixed angle
      real(8) :: dvt, dvdlt(3)                    ! dV/d(sin^2 phi), dV/dDl_i
      real(8) :: u1(3), u2(3), u3(3)              ! pair unit vectors (atom A -> B)
      real(8) :: nhh(3), w(3), wm                 ! HH axis unit, F->midpoint vector
      real(8) :: cph, tph                         ! cos(phi), sin^2(phi)
      real(8) :: dc(3, 3), dt(3, 3), gph(3, 3)    ! per-atom d(cos)/dq, d(sin^2)/dq, angle gradient
      real(8) :: what(3)                          ! unit F->midpoint vector
      ! 1. Geometry: pairs (F,H1) (F,H2) (H1,H2); phi between the H1-H2
      !    axis and the F-to-bond-midpoint vector
      xf = q(1:3)
      xh1 = q(4:6)
      xh2 = q(7:9)
      r1 = sqrt(dot_product(xh1 - xf, xh1 - xf))
      r2 = sqrt(dot_product(xh2 - xf, xh2 - xf))
      r3 = sqrt(dot_product(xh2 - xh1, xh2 - xh1))
      u1 = (xh1 - xf)/r1
      u2 = (xh2 - xf)/r2
      u3 = (xh2 - xh1)/r3
      nhh = u3
      w = 0.5d0*(xh1 + xh2) - xf
      wm = sqrt(dot_product(w, w))
      cph = dot_product(nhh, w)/wm
      tph = 1.0d0 - cph*cph
      dl1 = dl0_fh + amod*tph                  ! FH pairs share the modulated value
      dl2 = dl1
      dl3 = dl0_hh + amod*tph
      ! 2. Branches and their derivatives (pair i: D, r0, beta per kind)
      call branch(d_fh, r0_fh, b_fh, r1, dl1, q1, a1, dq1, da1, dqa1, daa1)
      call branch(d_fh, r0_fh, b_fh, r2, dl2, q2, a2, dq2, da2, dqa2, daa2)
      call branch(d_hh, r0_hh, b_hh, r3, dl3, q3, a3, dq3, da3, dqa3, daa3)
      ! 3. Assembly V = sum(Q') - sqrt(S); dV/dr_i at fixed angle
      s = a1*a1 + a2*a2 + a3*a3 - a1*a2 - a2*a3 - a3*a1
      v = q1 + q2 + q3 - sqrt(max(s, tiny_s))
      rs = sqrt(max(s, tiny_s))
      dv1 = dq1 - da1*(2.0d0*a1 - a2 - a3)/(2.0d0*rs)
      dv2 = dq2 - da2*(2.0d0*a2 - a1 - a3)/(2.0d0*rs)
      dv3 = dq3 - da3*(2.0d0*a3 - a1 - a2)/(2.0d0*rs)
      ! 4. Angle channel: dV/dDl_i then dV/d(sin^2 phi) = a_mod * sum
      dvdlt(1) = dqa1 - daa1*(2.0d0*a1 - a2 - a3)/(2.0d0*rs)
      dvdlt(2) = dqa2 - daa2*(2.0d0*a2 - a1 - a3)/(2.0d0*rs)
      dvdlt(3) = dqa3 - daa3*(2.0d0*a3 - a1 - a2)/(2.0d0*rs)
      dvt = amod*(dvdlt(1) + dvdlt(2) + dvdlt(3))
      ! 5. d(cos phi)/dq per atom, written as two orthogonal channels:
      !    axis channel (what - c*n)/r3 and radial channel (n - c*what)/m;
      !    F carries -radial, H1/H2 carry -/+ axis and +radial/2 (the
      !    three atom gradients sum to zero - translation invariance;
      !    both channels vanish in the collinear limit)
      what = w/wm
      dc(1, :) = -(nhh - cph*what)/wm
      dc(2, :) = -(what - cph*nhh)/r3 + 0.5d0*(nhh - cph*what)/wm
      dc(3, :) = (what - cph*nhh)/r3 + 0.5d0*(nhh - cph*what)/wm
      dt = -2.0d0*cph*dc
      gph = dvt*dt
      ! 6. Total gradient: pair-distance chain rule + angle channel
      g(1:3) = dv1*(-u1) + dv2*(-u2) + gph(1, :)
      g(4:6) = dv1*u1 + dv3*(-u3) + gph(2, :)
      g(7:9) = dv2*u2 + dv3*u3 + gph(3, :)
   end subroutine fh2_vg

   !------------------------------------------------------------------
   ! branch(...) - one pair: Q', A' and their d/dr, d/dDl at fixed r-set
   !------------------------------------------------------------------
   subroutine branch(dd, r0, bb, r, dl, qv, av, dq, da, dqdl, dadl)
      real(8), intent(in)  :: dd, r0, bb, r, dl
      real(8), intent(out) :: qv, av, dq, da, dqdl, dadl
      real(8) :: e, e2
      e = exp(-bb*(r - r0))
      e2 = e*e
      qv = dd/(4.0d0*(1.0d0 + dl))*((3.0d0 + dl)*e2 - (2.0d0 + 6.0d0*dl)*e)
      av = dd/(4.0d0*(1.0d0 + dl))*((1.0d0 + 3.0d0*dl)*e2 - (6.0d0 + 2.0d0*dl)*e)
      dq = dd*bb*e/(4.0d0*(1.0d0 + dl))*((2.0d0 + 6.0d0*dl) - 2.0d0*(3.0d0 + dl)*e)
      da = dd*bb*e/(4.0d0*(1.0d0 + dl))*((6.0d0 + 2.0d0*dl) - 2.0d0*(1.0d0 + 3.0d0*dl)*e)
      ! d/dDl: Q'+A' is Dl-independent (plain Morse) - the pair invariant
      dqdl = -dd*(e2 + 2.0d0*e)/(2.0d0*(1.0d0 + dl)**2)
      dadl = dd*(e2 + 2.0d0*e)/(2.0d0*(1.0d0 + dl)**2)
   end subroutine branch
end module fh2_pot
