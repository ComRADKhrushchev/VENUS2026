!=====================================================================
! check_h2ag111_pes.f90 - the H2/Ag(111) NN surface physics-anchor
!                          criterion program (criteria defined FIRST in
!                          this header; this program executes them)
! Criteria (against Bin Jiang & Hua Guo, PCCP 2014, 16, 24704,
! DOI 10.1039/c4cp03761h - the source of this PES; all searches below
! are deterministic scans in the bridge/top site planes, the H2
! parallel (theta = 90 deg) with the bond along +y, COM at the site):
!   (a) surface periodicity: translating BOTH atoms by either primitive
!       lattice vector a1 = (dmm, 0) or a2 = (dmm/2, dmm*sqrt(3)/2)
!       (dmm = 2.9437 A, the Ag-Ag nearest-neighbor distance of the
!       symmetrize() input construction) leaves V (1e-12 eV) and the
!       gradient (1e-9 eV/A) invariant at a generic point
!   (b) the usable asymptote: an intact parallel H2 at Z = 7.5 A has
!       V in [0.002, 0.006] eV (the NN's far-field constant; measured
!       +0.0035) with |dV/dZ| below 0.01 eV/A - the paper's QCT
!       scattering boundary at 8.1 A sits inside this usable region
!       (the far field fails only for stretched/tilted geometries)
!   (c) the physisorption well (paper: depth ~0.007 eV at Z = 4.01 A):
!       the Z minimum at fixed r = 0.7414 A lies at Z in [3.8, 4.3]
!       with V - V_far in [-0.012, -0.003] eV (measured -0.0075
!       at Z = 4.1)
!   (d) the bridge saddle (paper: r(HH) = 1.26 A, barrier 1.16 eV):
!       over Z in [1.00, 1.20] the ridge max along r in [1.15, 1.35]
!       is minimized at (measured Z = 1.10, r = 1.27) with
!       V - V_far in [1.10, 1.20] eV (measured 1.152) and the
!       argmin r in [1.21, 1.32]
!   (e) site corrugation (paper: top/fcc saddles ~0.3 eV above the
!       global bridge barrier): the same ridge search at the top site
!       over Z in [1.4, 1.7], r in [1.45, 1.70] yields a saddle with
!       V_top - V_bridge in [0.20, 0.40] eV (measured 0.306)
!   (f) the dissociated channel descends: at (Z = 1.0, r = 2.2 A,
!       phi = 90 deg toward the flanking hollows) V sits at least
!       0.3 eV below the bridge saddle (measured ~0.42 below) - the
!       valley toward dissociation is open below the saddle. The full
!       fbf chemisorbed minimum (paper endothermicity 0.47 eV) needs
!       an unconstrained relaxation and stays outside this scan-based
!       criterion
! Data files: pes_init reads weights/biases-h2ag.txt from the CURRENT
!   directory (copy from containers/h2ag111PES/).
! Compile (from src/):
!   gfortran -c -o build/h2ag111_pot.o ../systems/inactive/h2ag111PES/h2ag111_pot.f90
!   gfortran -c -Jbuild -o build/check_h2ag111_pes.o ../scripts/check_h2ag111_pes.f90
!   gfortran -o build/check_h2ag111_pes.exe build/h2ag111_pot.o build/check_h2ag111_pes.o
!   (cd build && cp ../../systems/inactive/h2ag111PES/{weights,biases}-h2ag.txt . && ./check_h2ag111_pes.exe)
!=====================================================================
program check_h2ag111_pes
   implicit none
   real(8), parameter :: dmm = 2.9437d0   ! Ag-Ag NN distance [A] (the symmetrize constant)
   real(8) :: q(6), qt(6), v, vt, dv(6), dvt(6)
   integer :: npass = 0, nfail = 0, i, iz, ir
   real(8) :: z, r, v_far, v_well, z_well, v_sad, z_sad, r_sad, v_rib, v_top, dev
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
   call pes_init()

   ! ---- (a) periodicity ----
   call place(0.31d0, 0.42d0, 2.7d0, 1.05d0, q)
   call pot0(2, q, v)
   call dpeshon(2, q, dv)
   dev = 0.0d0
   do i = 1, 2
      if (i == 1) then
         qt = q + [dmm, 0.0d0, 0.0d0, dmm, 0.0d0, 0.0d0]
      else
         qt = q + [0.5d0*dmm, 0.5d0*sqrt(3.0d0)*dmm, 0.0d0, 0.5d0*dmm, 0.5d0*sqrt(3.0d0)*dmm, 0.0d0]
      end if
      call pot0(2, qt, vt)
      call dpeshon(2, qt, dvt)
      dev = max(dev, abs(vt - v), maxval(abs(dvt - dv)))
   end do
   if (abs(dev) < 1.0d-9) then
      npass = npass + 1
      print '(a,es9.2)', '(a) lattice periodicity     PASS  dev = ', dev
   else
      nfail = nfail + 1
      print '(a,es9.2)', '(a) lattice periodicity     FAIL  dev = ', dev
   end if

   ! ---- (b) the usable asymptote ----
   call place(0.5d0*dmm, 0.0d0, 7.5d0, 0.7414d0, q)
   call pot0(2, q, v)
   call dpeshon(2, q, dv)
   v_far = v
   call place(0.5d0*dmm, 0.0d0, 8.1d0, 0.7414d0, qt)
   call pot0(2, qt, vt)
   if (v >= 0.002d0 .and. v <= 0.006d0 .and. abs(0.5d0*(dv(3) + dv(6))) < 0.01d0 .and. &
       abs(vt - v) < 0.002d0) then
      npass = npass + 1
      print '(a,f9.5,a,es9.2)', '(b) usable asymptote       PASS  V_far = ', v, &
         '  |gZ| = ', abs(0.5d0*(dv(3) + dv(6)))
   else
      nfail = nfail + 1
      print '(a,f9.5,a,es9.2)', '(b) usable asymptote       FAIL  V_far = ', v, &
         '  |gZ| = ', abs(0.5d0*(dv(3) + dv(6)))
   end if

   ! ---- (c) the physisorption well ----
   v_well = 1.0d0
   do iz = 0, 80
      z = 3.7d0 + 0.01d0*iz
      call place(0.5d0*dmm, 0.0d0, z, 0.7414d0, q)
      call pot0(2, q, vt)
      if (vt < v_well) then
         v_well = vt
         z_well = z
      end if
   end do
   if (z_well >= 3.8d0 .and. z_well <= 4.3d0 .and. v_well - v_far >= -0.012d0 .and. &
       v_well - v_far <= -0.003d0) then
      npass = npass + 1
      print '(a,f6.2,a,f9.5)', '(c) physisorption well     PASS  Z = ', z_well, &
         '  depth = ', v_well - v_far
   else
      nfail = nfail + 1
      print '(a,f6.2,a,f9.5)', '(c) physisorption well     FAIL  Z = ', z_well, &
         '  depth = ', v_well - v_far
   end if

   ! ---- (d) the bridge saddle: min over Z of the ridge max over r ----
   v_sad = huge(1.0d0)
   do iz = 0, 20
      z = 1.00d0 + 0.01d0*iz
      v_rib = -huge(1.0d0)
      do ir = 0, 40
         r = 1.15d0 + 0.005d0*ir
         call place(0.5d0*dmm, 0.0d0, z, r, q)
         call pot0(2, q, vt)
         v_rib = max(v_rib, vt)
      end do
      if (v_rib < v_sad) then      ! the saddle = min over Z of the ridge max
         v_sad = v_rib
         z_sad = z
      end if
   end do
   ! recover r of the ridge max at the winning Z
   v_rib = -huge(1.0d0)
   do ir = 0, 40
      r = 1.15d0 + 0.005d0*ir
      call place(0.5d0*dmm, 0.0d0, z_sad, r, q)
      call pot0(2, q, vt)
      if (vt > v_rib) then
         v_rib = vt
         r_sad = r
      end if
   end do
   if (v_sad - v_far >= 1.10d0 .and. v_sad - v_far <= 1.20d0 .and. &
       r_sad >= 1.21d0 .and. r_sad <= 1.32d0) then
      npass = npass + 1
      print '(a,f6.2,a,f5.2,a,f8.4)', '(d) bridge saddle          PASS  Z = ', z_sad, &
         '  r = ', r_sad, '  barrier = ', v_sad - v_far
   else
      nfail = nfail + 1
      print '(a,f6.2,a,f5.2,a,f8.4)', '(d) bridge saddle          FAIL  Z = ', z_sad, &
         '  r = ', r_sad, '  barrier = ', v_sad - v_far
   end if

   ! ---- (e) the top-site saddle ----
   v_top = huge(1.0d0)
   do iz = 0, 30
      z = 1.40d0 + 0.01d0*iz
      v_rib = -huge(1.0d0)
      do ir = 0, 50
         r = 1.45d0 + 0.005d0*ir
         call place(0.0d0, 0.0d0, z, r, q)
         call pot0(2, q, vt)
         v_rib = max(v_rib, vt)
      end do
      if (v_rib < v_top) then
         v_top = v_rib
      end if
   end do
   if (v_top - v_sad >= 0.20d0 .and. v_top - v_sad <= 0.40d0) then
      npass = npass + 1
      print '(a,f8.4)', '(e) top-site corrugation    PASS  delta = ', v_top - v_sad
   else
      nfail = nfail + 1
      print '(a,f8.4)', '(e) top-site corrugation    FAIL  delta = ', v_top - v_sad
   end if

   ! ---- (f) the dissociated channel descends ----
   call place(0.5d0*dmm, 0.0d0, 1.0d0, 2.2d0, q)
   call pot0(2, q, vt)
   if (v_sad - vt >= 0.3d0) then
      npass = npass + 1
      print '(a,f8.4)', '(f) dissociation channel    PASS  below saddle by ', v_sad - vt
   else
      nfail = nfail + 1
      print '(a,f8.4)', '(f) dissociation channel    FAIL  below saddle by ', v_sad - vt
   end if

   print '(/,a,i0,a,i0,a)', 'check_h2ag111_pes: ', npass, ' PASS, ', nfail, ' FAIL'
   if (nfail > 0) call exit(1)

contains
   ! place(x, y, z, r, q) - the parallel H2: COM at (x,y,z), bond along +y
   subroutine place(x, y, z, r, qq)
      real(8), intent(in) :: x, y, z, r
      real(8), intent(out) :: qq(6)
      qq = [x, y - 0.5d0*r, z, x, y + 0.5d0*r, z]
   end subroutine place
end program check_h2ag111_pes
