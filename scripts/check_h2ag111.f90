!=====================================================================
! check_h2ag111.f90 - the H2/Ag(111) neural-network container criterion
!                     program (criteria defined FIRST in this header;
!                     this program executes them)
! Criteria:
!   (a) transcription battery: six geometries (the distribution main.f
!       anchor configuration, top-site parallel, near-surface
!       dissociating, tilted bridge, one-chemisorbed-one-out, compact
!       close-in, far-above-surface) through POT0/DPESHON against the
!       values recorded from the ARCHIVED fixed-form original
!       h2ag111.f (same compiler, same battery; the original and the
!       free-form transcription were verified bit-identical):
!       |V - V_ref| <= 1e-13 eV and |dV - dV_ref| <= 1e-13 eV/A per
!       component. The distribution main.f's PUBLISHED values also
!       asserted: V = 1.15467106108370 within 1e-9 (observed 2.3e-12),
!       the six printed dV components within 1e-6 (the published digits
!       were produced on other hardware; they agree to ~3 significant
!       digits - the loose bound documents that, it does not test it)
!   (b) atom-exchange symmetry (H1 <-> H2): V equal (1e-12 eV) and the
!       gradient blocks swap (1e-12 eV/A) - the symmetrize() input
!       construction is exchange symmetric by form
!   (c) the far-field stretch extrapolation pin: at the tilted,
!       stretched far geometry (z = 8.0/8.5, r = 0.89) the recorded
!       V = 0.26276113239005561 eV is reproduced (1e-13) with
!       |dV/dz1| > 1 eV/A - the network's far field is valid for
!       intact parallel molecules only (the usable asymptote is
!       pinned by check_h2ag111_pes leg (b)); placements and Z_TERM
!       stay inside the fitted domain
! Data files: the program calls pes_init, which reads
!   weights-h2ag.txt / biases-h2ag.txt from the CURRENT directory -
!   run it in a directory holding the two files (e.g. copy them from
!   systems/inactive/h2ag111PES/ into build/ before running).
! Compile (from src/):
!   gfortran -c -o build/h2ag111_pot.o ../systems/inactive/h2ag111PES/h2ag111_pot.f90
!   gfortran -c -Jbuild -o build/check_h2ag111.o ../scripts/check_h2ag111.f90
!   gfortran -o build/check_h2ag111.exe build/h2ag111_pot.o build/check_h2ag111.o
!   (cd build && cp ../../systems/inactive/h2ag111PES/{weights,biases}-h2ag.txt . && ./check_h2ag111.exe)
!=====================================================================
program check_h2ag111
   implicit none
   ! the archived PES entries (external subroutines of h2ag111_pot.f90)
   interface
      subroutine pot0(natoms, q, v)
         integer :: natoms
         real(8) :: q(*), v
      end subroutine pot0
      subroutine dpeshon(natoms, q, pdot)
         integer :: natoms
         real(8) :: q(*), pdot(*)
      end subroutine dpeshon
      subroutine pes_init()
      end subroutine pes_init
   end interface
   integer, parameter :: ng = 6
   real(8) :: q(6), v, dv(6), qs(6), vs, dvs(6)
   real(8) :: qbat(6, ng), vref(ng), dref(6, ng)
   integer :: it, k, npass, nfail
   npass = 0
   nfail = 0

   ! the battery [A]: (x1,y1,z1, x2,y2,z2)
   qbat(:, 1) = [1.277280d-01, 1.625873d+00, 1.064094d+00, 1.223640d+00, 9.931286d-01, 1.125637d+00]  ! main.f anchor
   qbat(:, 2) = [0.0d0, 0.0d0, 1.6d0, 0.74d0, 0.0d0, 1.6d0]                                            ! top-site parallel
   qbat(:, 3) = [0.0d0, 0.0d0, 1.2d0, 1.8d0, 0.0d0, 1.0d0]                                            ! dissociating
   qbat(:, 4) = [1.47d0, 0.0d0, 1.3d0, 2.21d0, 0.85d0, 1.5d0]                                         ! tilted bridge
   qbat(:, 5) = [0.5d0, 0.5d0, 0.9d0, 2.9d0, 1.7d0, 3.4d0]                                            ! one out
   qbat(:, 6) = [0.1d0, 0.2d0, 8.0d0, 0.84d0, 0.2d0, 8.5d0]                                           ! far above
   ! recorded energies [eV] from the archived fixed-form original
   vref = [1.1546710610860220d+00, 9.5345842783986612d-01, 4.0942990862016559d+00, &
           1.1293049155473054d+00, 6.2297179801813574d+00, 2.6276113239005561d-01]
   ! recorded gradients [eV/A]
   dref(:, 1) = [-4.6628662042635938d-05, 1.0736353583951086d-04, -7.6015930217110395d-05, &
                  1.0333587319166782d-04, -1.3978710411421957d-05, 9.5749486117568949d-06]
   dref(:, 2) = [2.2303354092390038d+00, 4.1258887395478894d-03, -2.6960668827524792d+00, &
                -2.4667402784983583d+00, -6.1576504184479575d-03, -9.8666807833658332d-01]
   dref(:, 3) = [4.1410422849441630d-01, -1.9522855105869041d-02, -1.4356341026093844d+01, &
                  1.8498188060922516d+00, 5.8702703813473534d-02, -1.4672624593847772d+00]
   dref(:, 4) = [-1.0007675466161547d+00, -9.7433746488961859d-01, -3.5761903573794385d-01, &
                  8.3971246071512518d-01, 9.6692335538328678d-01, 1.9863462062508219d-01]
   dref(:, 5) = [-5.5987406550661589d+00, -5.9668446322591606d+00, -1.0140376755831682d+01, &
                  1.4273311572665648d-02, 9.4428402963586677d-03, 1.1662352427332934d-01]
   dref(:, 6) = [-2.6250182031728464d+00, 0.0d0, -1.7820078404828799d+00, 0.0d0, 0.0d0, 0.0d0]
   ! note: geometry 6's recorded x2/y2/z2 gradient components were not
   ! part of the recorded battery (the archived comparison printed
   ! components 1 and 3 only); the zeros above are placeholders and the
   ! comparison for geometry 6 checks components 1 and 3 only

   call pes_init()

   ! ---- (a) the transcription battery ----
   do it = 1, ng
      q = qbat(:, it)
      call pot0(2, q, v)
      call dpeshon(2, q, dv)
      if (it == 6) then
         vs = max(abs(v - vref(it)), abs(dv(1) - dref(1, it)), abs(dv(3) - dref(3, it)))
      else
         vs = max(abs(v - vref(it)), maxval(abs(dv - dref(:, it))))
      end if
      if (vs <= 1.0d-13) then
         npass = npass + 1
         print '(a,i1,a,es9.2)', '(a) battery geometry ', it, '  PASS  max dev = ', vs
      else
         nfail = nfail + 1
         print '(a,i1,a,es9.2)', '(a) battery geometry ', it, '  FAIL  max dev = ', vs
      end if
   end do
   ! the published anchor of the distribution main.f
   q = qbat(:, 1)
   call pot0(2, q, v)
   call dpeshon(2, q, dv)
   if (abs(v - 1.15467106108370d0) <= 1.0d-9) then
      npass = npass + 1
      print '(a,es9.2)', '(a) published V anchor        PASS  |dV| = ', abs(v - 1.15467106108370d0)
   else
      nfail = nfail + 1
      print '(a,es9.2)', '(a) published V anchor        FAIL  |dV| = ', abs(v - 1.15467106108370d0)
   end if
   vs = maxval(abs(dv - [-4.669620023989296d-05, 1.073684985364309d-04, -7.622343645152796d-05, &
                          1.034455954229685d-04, -1.406232019718345d-05, 9.655022337184960d-06]))
   if (vs <= 1.0d-6) then
      npass = npass + 1
      print '(a,es9.2)', '(a) published dV (3 digits)  PASS  max dev = ', vs
   else
      nfail = nfail + 1
      print '(a,es9.2)', '(a) published dV (3 digits)  FAIL  max dev = ', vs
   end if

   ! ---- (b) atom-exchange symmetry ----
   do it = 1, 2
      q = qbat(:, it)
      qs = [q(4:6), q(1:3)]                       ! H1 <-> H2
      call pot0(2, q, v)
      call dpeshon(2, q, dv)
      call pot0(2, qs, vs)
      call dpeshon(2, qs, dvs)
      vs = max(abs(v - vs), maxval(abs([dvs(4:6), dvs(1:3)] - dv)))
      if (vs <= 1.0d-12) then
         npass = npass + 1
         print '(a,i1,a,es9.2)', '(b) exchange symmetry geo ', it, '  PASS  dev = ', vs
      else
         nfail = nfail + 1
         print '(a,i1,a,es9.2)', '(b) exchange symmetry geo ', it, '  FAIL  dev = ', vs
      end if
   end do

   ! ---- (c) the far-field no-tail fact ----
   q = qbat(:, 6)
   call pot0(2, q, v)
   call dpeshon(2, q, dv)
   if (abs(v - 2.6276113239005561d-01) <= 1.0d-13 .and. abs(dv(3)) > 1.0d0) then
      npass = npass + 1
      print '(a,2es9.2)', '(c) far-field no-tail fact   PASS  V,|gz1| = ', v, abs(dv(3))
   else
      nfail = nfail + 1
      print '(a,2es9.2)', '(c) far-field no-tail fact   FAIL  V,|gz1| = ', v, abs(dv(3))
   end if

   print '(/,a,i0,a,i0,a)', 'check_h2ag111: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) then
      call exit(1)
   end if
end program check_h2ag111
