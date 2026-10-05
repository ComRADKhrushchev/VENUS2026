!=====================================================================
! check_fh2_leps.f90 - the F+H2 modified-LEPS container criterion
!                      program (criteria defined FIRST in this header;
!                      this program executes them)
! Criteria (against docs/benchmarks/fh2_leps_anchors.md):
!   (a) analytic gradient vs central difference over a geometry
!       battery (bent/collinear saddles, product valley, bent random,
!       near-asymptote, stretched H-F-H): max |g_an - g_fd| <= 5e-7 eV/A
!       (h = 5e-6 A: the H-F-H exchange region is curvature-stiff and
!       the FD truncation error scales h^2 - verified 16x drop for 4x h)
!   (b) the pair invariant Q'+A' = Morse: a lone pair degenerates to
!       plain Morse independent of the Sato parameter - F far away at
!       two orientations, H2 swept; and the HF product channel:
!       |V - D(e^2-2e)| <= 1e-12 eV
!   (c) translation invariance of V (1e-13 eV) and zero net force
!       (1e-12 eV/A); rotation invariance of V (1e-10 eV)
!   (d) the anchor battery (kcal/mol via the fixed 23.0605, a
!       criterion-driver-local constant - the T4 unit contract governs
!       the library, not the checks):
!         d1 exothermicity 31.747 +- 0.01 kcal (141.196 - 109.449)
!         d2 Muckerman Surface 2 (Dl 0.163/0.163, amod 0): collinear
!            saddle (rFH 1.630, rHH 0.755) barrier 0.893 +- 0.03 kcal
!            (published 0.89, [P6] Table I)
!         d3 M5 (Dl 0.106/0.167, amod 0): saddle (1.5413, 0.7618)
!            barrier 1.064 +- 0.05 kcal (published saddle geometry
!            1.54/0.76, barrier in the 0.9-1.8 kcal literature range)
!         d4 T-S fixed set: bent saddle F(1.71258,1.22425),
!            rHH 0.75874 -> barrier 0.1794 +- 0.02 kcal (published
!            "negligible, 8 meV at a bent configuration", [P1] p.228)
!         d5 T-S collinear saddle (1.6587, 0.7537) barrier 0.7781
!            +- 0.03 kcal; and bent below collinear by > 0.4 kcal
!            (the bent transition state, [P2]/[P3])
!   (e) parameter setter round-trip: re-setting the fixed set
!       changes V by <= 1e-15 eV
! Compile (from src/ after make; modules land under build/):
!   gfortran -c -Jbuild -o build/fh2_pot.o ../systems/inactive/fh2_leps/fh2_pot.f90
!   gfortran -c -Ibuild -Jbuild -o build/check_fh2_leps.o ../scripts/check_fh2_leps.f90
!   gfortran -o build/check_fh2_leps.exe build/fh2_pot.o build/check_fh2_leps.o
!   ./build/check_fh2_leps.exe
!=====================================================================
program check_fh2_leps
   use fh2_pot, only: fh2_vg, fh2_set_params, fh2_set_defaults
   implicit none
   integer, parameter :: ntest = 6
   real(8) :: q(9), v, g(9), v2, g2(9), vfd, gnum(9), dq(9)
   real(8) :: vas, kcal, exo, bar, tol, hh
   integer :: it, k, npass, nfail, ios
   npass = 0
   nfail = 0
   kcal = 23.0605d0            ! the criterion-driver-local kcal->eV factor (checks only)
   hh = 5.0d-6                 ! central-difference step [A] (see header note)

   ! ---- (a) analytic gradient vs central difference ----
   do it = 1, ntest
      select case (it)
      case (1)   ! T-S bent saddle
         q = [1.71258d0, 1.22425d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.75874d0, 0.0d0, 0.0d0]
      case (2)   ! collinear saddle region
         q = [0.0d0, 0.0d0, 0.0d0, 1.6587d0, 0.0d0, 0.0d0, 2.4124d0, 0.0d0, 0.0d0]
      case (3)   ! product valley (HF formed, H leaving)
         q = [0.0d0, 0.0d0, 0.0d0, 0.917d0, 0.03d0, 0.0d0, 6.0d0, 1.0d0, -0.5d0]
      case (4)   ! generic bent
         q = [0.4d0, -1.1d0, 0.7d0, 1.3d0, 0.2d0, 0.9d0, 1.9d0, -0.6d0, 0.1d0]
      case (5)   ! near-asymptote
         q = [15.0d0, 2.0d0, -1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.7419d0, 0.0d0, 0.0d0]
      case (6)   ! stretched exchange region
         q = [0.0d0, 0.0d0, 0.0d0, 1.2d0, 0.0d0, 0.0d0, 1.2d0, 0.05d0, 0.0d0]
      end select
      call fh2_vg(3, q, v, g)
      do k = 1, 9
         dq = 0.0d0
         dq(k) = hh
         call fh2_vg(3, q + dq, v2, g2)
         vfd = v2
         dq(k) = -hh
         call fh2_vg(3, q + dq, v2, g2)
         vfd = (vfd - v2)/(2.0d0*hh)
         gnum(k) = vfd
      end do
      tol = 5.0d-7
      if (maxval(abs(g - gnum)) <= tol) then
         npass = npass + 1
         print '(a,i1,a,es9.2)', '(a) gradient geometry ', it, '  PASS  max|dg| = ', maxval(abs(g - gnum))
      else
         nfail = nfail + 1
         print '(a,i1,a,es9.2)', '(a) gradient geometry ', it, '  FAIL  max|dg| = ', maxval(abs(g - gnum))
      end if
   end do

   ! ---- (b) the lone-pair Morse invariant ----
   call fh2_vg(3, [20.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.74d0, 0.0d0, 0.0d0], v, g)
   v2 = 4.7461676893d0*(exp(-1.942d0*(0.74d0 - 0.7419d0))**2 - 2.0d0*exp(-1.942d0*(0.74d0 - 0.7419d0)))
   if (abs(v - v2) <= 1.0d-12) then
      npass = npass + 1
      print '(a,es9.2)', '(b) H2 lone pair (collinear)   PASS  |dV| = ', abs(v - v2)
   else
      nfail = nfail + 1
      print '(a,es9.2)', '(b) H2 lone pair (collinear)   FAIL  |dV| = ', abs(v - v2)
   end if
   call fh2_vg(3, [20.0d0, 7.0d0, 3.0d0, 0.0d0, 0.0d0, 0.0d0, 0.80d0, 0.0d0, 0.0d0], v, g)
   v2 = 4.7461676893d0*(exp(-1.942d0*(0.80d0 - 0.7419d0))**2 - 2.0d0*exp(-1.942d0*(0.80d0 - 0.7419d0)))
   if (abs(v - v2) <= 1.0d-12) then
      npass = npass + 1
      print '(a)', '(b) H2 lone pair (bent F)      PASS'
   else
      nfail = nfail + 1
      print '(a,es9.2)', '(b) H2 lone pair (bent F)      FAIL  |dV| = ', abs(v - v2)
   end if
   call fh2_vg(3, [0.0d0, 0.0d0, 0.0d0, 0.90d0, 0.0d0, 0.0d0, 20.0d0, 5.0d0, -2.0d0], v, g)
   v2 = 6.1228507621d0*(exp(-2.2187d0*(0.90d0 - 0.917d0))**2 - 2.0d0*exp(-2.2187d0*(0.90d0 - 0.917d0)))
   if (abs(v - v2) <= 1.0d-12) then
      npass = npass + 1
      print '(a)', '(b) HF product channel         PASS'
   else
      nfail = nfail + 1
      print '(a,es9.2)', '(b) HF product channel         FAIL  |dV| = ', abs(v - v2)
   end if

   ! ---- (c) translation / rotation invariance, zero net force ----
   q = [0.4d0, -1.1d0, 0.7d0, 1.3d0, 0.2d0, 0.9d0, 1.9d0, -0.6d0, 0.1d0]
   call fh2_vg(3, q, v, g)
   call fh2_vg(3, q + [1.7d0, -0.3d0, 4.2d0, 1.7d0, -0.3d0, 4.2d0, 1.7d0, -0.3d0, 4.2d0], v2, g2)
   tol = abs(v - v2)
   v2 = max(abs(g(1) + g(4) + g(7)), abs(g(2) + g(5) + g(8)), abs(g(3) + g(6) + g(9)))
   if (tol <= 1.0d-13 .and. v2 <= 1.0d-12) then
      npass = npass + 1
      print '(a,es9.2)', '(c) translation / zero force  PASS  |dV|,|sum g| = ', max(tol, v2)
   else
      nfail = nfail + 1
      print '(a,es9.2,es9.2)', '(c) translation / zero force  FAIL  |dV|,|sum g| = ', tol, v2
   end if
   ! rotate ~31.7 deg about (1,1,1)/sqrt(3) (Rodrigues), r . q
   block
      real(8) :: ax(3), cth, sth, rmat(3, 3)
      ax = [1.0d0, 1.0d0, 1.0d0]/sqrt(3.0d0)
      cth = cos(0.553d0)
      sth = sin(0.553d0)
      rmat = reshape([cth + ax(1)*ax(1)*(1 - cth), ax(2)*ax(1)*(1 - cth) + ax(3)*sth, ax(3)*ax(1)*(1 - cth) - ax(2)*sth, &
                      ax(1)*ax(2)*(1 - cth) - ax(3)*sth, cth + ax(2)*ax(2)*(1 - cth), ax(3)*ax(2)*(1 - cth) + ax(1)*sth, &
                      ax(1)*ax(3)*(1 - cth) + ax(2)*sth, ax(2)*ax(3)*(1 - cth) - ax(1)*sth, cth + ax(3)*ax(3)*(1 - cth)], [3, 3])
      call fh2_vg(3, [matmul(rmat, q(1:3)), matmul(rmat, q(4:6)), matmul(rmat, q(7:9))], v2, g2)
      if (abs(v - v2) <= 1.0d-10) then
         npass = npass + 1
         print '(a,es9.2)', '(c) rotation invariance       PASS  |dV| = ', abs(v - v2)
      else
         nfail = nfail + 1
         print '(a,es9.2)', '(c) rotation invariance       FAIL  |dV| = ', abs(v - v2)
      end if
   end block

   ! ---- (d) the anchor battery ----
   ! reactant asymptote (F far, collinear): V = -D_HH exactly (pair invariant)
   call fh2_vg(3, [50.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.7419d0, 0.0d0, 0.0d0], vas, g)
   ! d1: exothermicity (HF + H channel)
   call fh2_vg(3, [0.0d0, 0.0d0, 0.0d0, 0.917d0, 0.0d0, 0.0d0, 50.0d0, 0.0d0, 0.0d0], v, g)
   exo = (vas - v)*kcal
   if (abs(exo - 31.747d0) <= 0.01d0) then
      npass = npass + 1
      print '(a,f8.4)', '(d1) exothermicity            PASS  ', exo
   else
      nfail = nfail + 1
      print '(a,f8.4)', '(d1) exothermicity            FAIL  ', exo
   end if
   ! d2: Muckerman Surface 2 (Dl 0.163/0.163, amod 0)
   call fh2_set_params(6.1228507621d0, 0.917d0, 2.2187d0, 4.7461676893d0, 0.7419d0, 1.942d0, 0.163d0, 0.163d0, 0.0d0)
   call fh2_vg(3, [0.0d0, 0.0d0, 0.0d0, 1.630d0, 0.0d0, 0.0d0, 2.385d0, 0.0d0, 0.0d0], v, g)
   call fh2_vg(3, [50.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.7419d0, 0.0d0, 0.0d0], vas, g)
   bar = (v - vas)*kcal
   if (abs(bar - 0.893d0) <= 0.03d0) then
      npass = npass + 1
      print '(a,f8.4)', '(d2) Surface-2 barrier        PASS  ', bar
   else
      nfail = nfail + 1
      print '(a,f8.4)', '(d2) Surface-2 barrier        FAIL  ', bar
   end if
   ! d3: M5 (Dl 0.106/0.167, amod 0)
   call fh2_set_params(6.1228507621d0, 0.917d0, 2.2187d0, 4.7461676893d0, 0.7419d0, 1.942d0, 0.106d0, 0.167d0, 0.0d0)
   call fh2_vg(3, [0.0d0, 0.0d0, 0.0d0, 1.5413d0, 0.0d0, 0.0d0, 2.3031d0, 0.0d0, 0.0d0], v, g)
   call fh2_vg(3, [50.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.7419d0, 0.0d0, 0.0d0], vas, g)
   bar = (v - vas)*kcal
   if (abs(bar - 1.064d0) <= 0.05d0) then
      npass = npass + 1
      print '(a,f8.4)', '(d3) M5 barrier               PASS  ', bar
   else
      nfail = nfail + 1
      print '(a,f8.4)', '(d3) M5 barrier               FAIL  ', bar
   end if
   ! d4/d5: the T-S fixed set (defaults)
   call fh2_set_defaults()
   call fh2_vg(3, [1.71258d0, 1.22425d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.75874d0, 0.0d0, 0.0d0], v, g)
   call fh2_vg(3, [50.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.7419d0, 0.0d0, 0.0d0], vas, g)
   bar = (v - vas)*kcal
   if (abs(bar - 0.1794d0) <= 0.02d0) then
      npass = npass + 1
      print '(a,f8.4,a)', '(d4) T-S bent barrier         PASS  ', bar, ' kcal (published 8 meV)'
   else
      nfail = nfail + 1
      print '(a,f8.4)', '(d4) T-S bent barrier         FAIL  ', bar
   end if
   call fh2_vg(3, [0.0d0, 0.0d0, 0.0d0, 1.6587d0, 0.0d0, 0.0d0, 2.4124d0, 0.0d0, 0.0d0], v, g)
   bar = (v - vas)*kcal
   if (abs(bar - 0.7781d0) <= 0.03d0) then
      npass = npass + 1
      print '(a,f8.4)', '(d5) T-S collinear barrier    PASS  ', bar
   else
      nfail = nfail + 1
      print '(a,f8.4)', '(d5) T-S collinear barrier    FAIL  ', bar
   end if
   if (bar - 0.1794d0 > 0.4d0) then
      npass = npass + 1
      print '(a)', '(d5) bent below collinear      PASS'
   else
      nfail = nfail + 1
      print '(a)', '(d5) bent below collinear      FAIL'
   end if

   ! ---- (e) setter round-trip ----
   call fh2_vg(3, [0.4d0, -1.1d0, 0.7d0, 1.3d0, 0.2d0, 0.9d0, 1.9d0, -0.6d0, 0.1d0], v, g)
   call fh2_set_params(6.1228507621d0, 0.917d0, 2.2187d0, 4.7461676893d0, 0.7419d0, 1.942d0, 0.173d0, 0.165d0, 0.090d0)
   call fh2_vg(3, [0.4d0, -1.1d0, 0.7d0, 1.3d0, 0.2d0, 0.9d0, 1.9d0, -0.6d0, 0.1d0], v2, g2)
   if (abs(v - v2) <= 1.0d-15) then
      npass = npass + 1
      print '(a)', '(e) setter round-trip         PASS'
   else
      nfail = nfail + 1
      print '(a,es9.2)', '(e) setter round-trip         FAIL  |dV| = ', abs(v - v2)
   end if

   print '(/,a,i0,a,i0,a)', 'check_fh2_leps: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) then
      call exit(1)
   end if
end program check_fh2_leps
