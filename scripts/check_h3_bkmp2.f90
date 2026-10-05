!=====================================================================
! check_h3_bkmp2.f90 - the H+H2 BKMP2 container criterion program
!                       (criteria defined FIRST in this header; this
!                       program executes them)
! Criteria:
!   (a) transcription battery: eleven geometries (pair-limit, collinear
!       saddle, equilateral/scalene compact, near-linear,
!       near-asymptote, isosceles - covering every surface branch
!       including the compact corrections) evaluated through bkmp2_vg in
!       A/eV against the values recorded from the ARCHIVED fixed-form
!       original systems/bkmp2.f (same compiler, same battery; the
!       original and the transcription were verified bit-identical):
!       |V - V_ref*ha_ev| <= 1e-9 eV, |g - g_ref*(ha_ev/bohr_a)| <=
!       1e-7 eV/A (headroom over the observed 0 for flag variation)
!   (b) analytic gradient vs central difference over three geometries
!       (h = 1e-5 A): max |g_an - g_fd| <= 2e-6 eV/A
!   (c) permutation symmetry (the H3 surface): atoms (1,2,3) -> (3,2,1)
!       leaves V (1e-12 eV) and permutes g (1e-10 eV/A)
!   (d) translation invariance of V (1e-12 eV), zero net force
!       (1e-10 eV/A); rotation invariance of V (1e-10 eV)
!   (e) the collinear saddle anchor: barrier above the H2+H asymptote
!       at r(HH) = 1.757 bohr = 9.611 +- 0.02 kcal/mol (published
!       9.61 kcal/mol; kcal via the criterion-driver-local 23.0605 -
!       the T4 unit contract governs the library, not the checks)
! Recorded anchors (systems/bkmp2.f, legacy units Eh and Eh/bohr):
!   geo i: (r12, r23, r13) [bohr]; V [Eh]; g1, g2, g3 [Eh/bohr]
!   H2 pair minimum: r = 1.4014719756078895 bohr, V = -0.17449577089710522 Eh
! Compile (from src/ after make; modules land under build/):
!   gfortran -c -Jbuild -o build/bkmp2_pot.o ../systems/inactive/h3_bkmp2/bkmp2_pot.f90
!   gfortran -c -Ibuild -Jbuild -o build/check_h3_bkmp2.o ../scripts/check_h3_bkmp2.f90
!   gfortran -o build/check_h3_bkmp2.exe build/bkmp2_pot.o build/check_h3_bkmp2.o
!   ./build/check_h3_bkmp2.exe
!=====================================================================
program check_h3_bkmp2
   use bkmp2_pot, only: bkmp2_vg
   implicit none
   integer, parameter :: ng = 11
   real(8), parameter :: bohr_a = 0.529177210903d0   ! 1 bohr [A] (the module's literal, restated)
   real(8), parameter :: ha_ev  = 27.211386245988d0  ! 1 hartree [eV] (likewise)
   real(8) :: trip(3, ng), vref(ng), gref(9, ng)
   real(8) :: q(9), qs(9), v, g(9), v2, g2(9), vfd, dq(9), gnum(9), kcal
   integer :: it, k, npass, nfail
   npass = 0
   nfail = 0
   kcal = 23.0605d0            ! the criterion-driver-local kcal->eV factor (checks only)

   ! ---- the battery: (r12, r23, r13) [bohr] ----
   trip(:, 1) = [1.401d0, 100.0d0, 100.0d0]    ! H2 near equilibrium, far H
   trip(:, 2) = [1.757d0, 1.757d0, 3.514d0]    ! collinear symmetric saddle
   trip(:, 3) = [1.4d0, 1.4d0, 1.4d0]          ! equilateral compact
   trip(:, 4) = [0.9d0, 1.4d0, 1.2d0]          ! compact scalene
   trip(:, 5) = [1.8d0, 2.7d0, 1.9d0]          ! generic bent scalene
   trip(:, 6) = [4.0d0, 8.0d0, 5.0d0]          ! near-asymptote
   trip(:, 7) = [1.2d0, 3.0d0, 2.5d0]          ! one compact pair
   trip(:, 8) = [2.0d0, 2.0d0, 3.9d0]          ! near-linear
   trip(:, 9) = [0.7d0, 1.3d0, 1.3d0]          ! deep compact
   trip(:, 10) = [2.3d0, 1.1d0, 2.9d0]         ! mixed compact scalene
   trip(:, 11) = [1.5489d0, 1.5489d0, 2.6d0]   ! isosceles
   ! recorded energies [Eh] from the archived original
   vref = [ -1.7449572968772770d-01, -1.5917600012106270d-01, -2.0887270859178741d-02, &
             1.0345050633156369d-01, -1.3147395124917829d-01, -1.6763182729284447d-02, &
            -1.4271434230899468d-01, -1.5039616787382096d-01,  2.6325206924979261d-01, &
            -1.2261720721672316d-01, -1.3380046072171428d-01 ]
   ! recorded gradients [Eh/bohr], components (g1x,g1y,g1z, g2x,..., g3z)
   gref = 0.0d0
   gref(1, 1) = 1.7467217954735917d-04; gref(2, 1) = -9.7732619056172556d-14
   gref(4, 1) = -1.7467217954735917d-04; gref(5, 1) = -9.7732619056172569d-14
   gref(7, 1) = -2.2876966251409342d-28; gref(8, 1) = 1.9546523811234514d-13
   gref(1, 2) = -8.2747010632377993d-11; gref(7, 2) = 8.2747010632377993d-11
   gref(1, 3) = 1.2886567399382123d-01; gref(2, 3) = 7.4400631569635237d-02
   gref(4, 3) = -1.2886567399382123d-01; gref(5, 3) = 7.4400631569635237d-02
   gref(8, 3) = -1.4880126313927047d-01
   gref(1, 4) = 5.3022870851349324d-01; gref(2, 4) = 1.9481747606397756d-01
   gref(4, 4) = -5.7422037654814018d-01; gref(5, 4) = 1.1327726664591210d-01
   gref(7, 4) = 4.3991668034646908d-02; gref(8, 4) = -3.0809474270988968d-01
   gref(1, 5) = -3.6495604051172595d-02; gref(2, 5) = -2.4912657145978918d-02
   gref(4, 5) = 7.5905887072310295d-03; gref(5, 5) = 3.0095720362580045d-02
   gref(7, 5) = 2.8905015343941565d-02; gref(8, 5) = -5.1830632166011276d-03
   gref(1, 6) = -2.2070531151451064d-02; gref(2, 6) = -5.3822335944636818d-04
   gref(4, 6) = 2.2477547227058792d-02; gref(5, 6) = -1.7107751222437650d-05
   gref(7, 6) = -4.0701607560772842d-04; gref(8, 6) = 5.5533111066880588d-04
   gref(1, 7) = 1.0679613449936212d-01; gref(2, 7) = 1.2672623160631401d-02
   gref(4, 7) = -1.2002567379325890d-01; gref(5, 7) = 1.4525313658577968d-02
   gref(7, 7) = 1.3229539293896791d-02; gref(8, 7) = -2.7197936819209369d-02
   gref(1, 8) = -2.6668170653630113d-02; gref(2, 8) = 2.9081397627080201d-03
   gref(4, 8) = 3.8935772927253184d-03; gref(5, 8) = -1.7084405143124724d-02
   gref(7, 8) = 2.2774593360904795d-02; gref(8, 8) = 1.4176265380416704d-02
   gref(1, 9) = 1.2825447862386348d+00; gref(2, 9) = 1.4866266735407482d-01
   gref(4, 9) = -1.2825447862386348d+00; gref(5, 9) = 1.4866266735407482d-01
   gref(7, 9) = -1.3877787807814457d-17; gref(8, 9) = -2.9732533470814965d-01
   gref(1, 10) = 3.6750922480466706d-02; gref(2, 10) = 7.3429731064002944d-03
   gref(4, 10) = 5.7951393696542682d-02; gref(5, 10) = 1.8430891827214260d-01
   gref(7, 10) = -9.4702316177009388d-02; gref(8, 10) = -1.9165189137854288d-01
   gref(1, 11) = 5.3311512917634696d-02; gref(2, 11) = 2.3373551977674729d-02
   gref(4, 11) = -1.0183639469704860d-02; gref(5, 11) = 1.5721548399219964d-02
   gref(7, 11) = -4.3127873447929836d-02; gref(8, 11) = -3.9095100376894690d-02

   ! ---- (a) transcription battery ----
   do it = 1, ng
      call place(trip(1, it), trip(2, it), trip(3, it), q)
      call bkmp2_vg(3, q, v, g)
      v2 = abs(v - vref(it)*ha_ev)
      g2 = abs(g - gref(:, it)*(ha_ev/bohr_a))
      if (v2 <= 1.0d-9 .and. maxval(g2) <= 1.0d-7) then
         npass = npass + 1
         print '(a,i2,a,es9.2,a,es9.2)', '(a) battery geometry ', it, '  PASS  |dV| = ', v2, &
            '  max|dg| = ', maxval(g2)
      else
         nfail = nfail + 1
         print '(a,i2,a,es9.2,a,es9.2)', '(a) battery geometry ', it, '  FAIL  |dV| = ', v2, &
            '  max|dg| = ', maxval(g2)
      end if
   end do

   ! ---- (b) analytic gradient vs central difference ----
   do it = 1, 3
      select case (it)
      case (1)   ! generic bent
         q = [0.4d0, -1.1d0, 0.7d0, 1.3d0, 0.2d0, 0.9d0, 1.9d0, -0.6d0, 0.1d0]
      case (2)   ! near the collinear saddle
         q = [0.0d0, 0.0d0, 0.0d0, 0.9299d0, 0.0d0, 0.0d0, 1.8598d0, 0.02d0, 0.0d0]
      case (3)   ! compact
         q = [0.3d0, 0.1d0, 0.0d0, 0.77d0, 0.05d0, 0.0d0, 1.0d0, 0.5d0, 0.2d0]
      end select
      call bkmp2_vg(3, q, v, g)
      do k = 1, 9
         dq = 0.0d0
         dq(k) = 1.0d-5
         call bkmp2_vg(3, q + dq, v2, g2)
         vfd = v2
         dq(k) = -1.0d-5
         call bkmp2_vg(3, q + dq, v2, g2)
         gnum(k) = (vfd - v2)/(2.0d-5)
      end do
      v2 = maxval(abs(g - gnum))
      if (v2 <= 2.0d-6) then
         npass = npass + 1
         print '(a,i1,a,es9.2)', '(b) finite difference ', it, '     PASS  max|dg| = ', v2
      else
         nfail = nfail + 1
         print '(a,i1,a,es9.2)', '(b) finite difference ', it, '     FAIL  max|dg| = ', v2
      end if
   end do

   ! ---- (c) permutation symmetry: (1,2,3) -> (3,2,1) ----
   q = [0.4d0, -1.1d0, 0.7d0, 1.3d0, 0.2d0, 0.9d0, 1.9d0, -0.6d0, 0.1d0]
   call bkmp2_vg(3, q, v, g)
   qs = [q(7:9), q(4:6), q(1:3)]
   call bkmp2_vg(3, qs, v2, g2)
   v2 = abs(v - v2)
   g2 = abs([g2(7:9), g2(4:6), g2(1:3)] - g)   ! atom map (1,2,3)->(3,2,1): blocks 1<->3 swap
   if (v2 <= 1.0d-12 .and. maxval(g2) <= 1.0d-10) then
      npass = npass + 1
      print '(a,es9.2)', '(c) permutation symmetry    PASS  |dV|,max|dg| = ', max(v2, maxval(g2))
   else
      nfail = nfail + 1
      print '(a,es9.2)', '(c) permutation symmetry    FAIL  |dV|,max|dg| = ', max(v2, maxval(g2))
   end if

   ! ---- (d) translation / zero force / rotation ----
   call bkmp2_vg(3, q + [1.7d0, -0.3d0, 4.2d0, 1.7d0, -0.3d0, 4.2d0, 1.7d0, -0.3d0, 4.2d0], v2, g2)
   vfd = abs(v - v2)
   v2 = max(abs(g(1) + g(4) + g(7)), abs(g(2) + g(5) + g(8)), abs(g(3) + g(6) + g(9)))
   if (vfd <= 1.0d-12 .and. v2 <= 1.0d-10) then
      npass = npass + 1
      print '(a,es9.2)', '(d) translation / zero force PASS  |dV|,|sum g| = ', max(vfd, v2)
   else
      nfail = nfail + 1
      print '(a,es9.2,es9.2)', '(d) translation / zero force FAIL  |dV|,|sum g| = ', vfd, v2
   end if
   ! rotate ~31.7 deg about (1,1,1)/sqrt(3) (Rodrigues)
   block
      real(8) :: ax(3), cth, sth, rmat(3, 3)
      ax = [1.0d0, 1.0d0, 1.0d0]/sqrt(3.0d0)
      cth = cos(0.553d0)
      sth = sin(0.553d0)
      rmat = reshape([cth + ax(1)*ax(1)*(1 - cth), ax(2)*ax(1)*(1 - cth) + ax(3)*sth, ax(3)*ax(1)*(1 - cth) - ax(2)*sth, &
                      ax(1)*ax(2)*(1 - cth) - ax(3)*sth, cth + ax(2)*ax(2)*(1 - cth), ax(3)*ax(2)*(1 - cth) + ax(1)*sth, &
                      ax(1)*ax(3)*(1 - cth) + ax(2)*sth, ax(2)*ax(3)*(1 - cth) - ax(1)*sth, cth + ax(3)*ax(3)*(1 - cth)], [3, 3])
      call bkmp2_vg(3, [matmul(rmat, q(1:3)), matmul(rmat, q(4:6)), matmul(rmat, q(7:9))], v2, g2)
      if (abs(v - v2) <= 1.0d-10) then
         npass = npass + 1
         print '(a,es9.2)', '(d) rotation invariance     PASS  |dV| = ', abs(v - v2)
      else
         nfail = nfail + 1
         print '(a,es9.2)', '(d) rotation invariance     FAIL  |dV| = ', abs(v - v2)
      end if
   end block

   ! ---- (e) the collinear saddle anchor ----
   call place(1.757d0, 1.757d0, 3.514d0, q)
   call bkmp2_vg(3, q, v, g)
   call place(1.4014719756078895d0, 100.0d0, 100.0d0, qs)   ! the H2 pair minimum, far H
   call bkmp2_vg(3, qs, v2, g2)
   v = (v - v2)*kcal                                          ! barrier above H2 + H [kcal/mol]
   if (abs(v - 9.611d0) <= 0.02d0) then
      npass = npass + 1
      print '(a,f8.4)', '(e) collinear barrier        PASS  ', v
   else
      nfail = nfail + 1
      print '(a,f8.4)', '(e) collinear barrier        FAIL  ', v
   end if

   print '(/,a,i0,a,i0,a)', 'check_h3_bkmp2: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) then
      call exit(1)
   end if
contains
   ! place(a, b, c, qq) - atoms realizing (r12, r23, r13), in A
   subroutine place(a, b, c, qq)
      real(8), intent(in) :: a, b, c
      real(8), intent(out) :: qq(9)
      real(8) :: x3, y3
      qq = 0.0d0
      qq(4) = a*bohr_a                       ! atom 2 on the x axis at r12
      x3 = (a*a + c*c - b*b)/(2.0d0*a)       ! atom 3 satisfies r13, r23
      y3 = sqrt(max(c*c - x3*x3, 0.0d0))
      qq(7) = x3*bohr_a
      qq(8) = y3*bohr_a
   end subroutine place
end program check_h3_bkmp2
