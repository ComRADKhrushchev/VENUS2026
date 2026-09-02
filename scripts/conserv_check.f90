!******************************************************************************
! conserv_check — offline per-frame total-energy reconstruction for
! TEST_PES=RST runs.
!
! Why: the engine refreshes T/V/H only at trajectory start/end
! (ENERGY_1 is not called inside the integration loop), so the per-step
! E0/T/H lines in fort.1001 are start-of-trajectory values, not live
! energies. This tool rebuilds the true per-frame total energy:
!   T  from the printed momenta (fort.1001 GWRITE_LEVEL>=2 atom lines,
!      columns 5-7 = P in amu*A/(10 fs)),
!   V  by calling rst_total_vg — the very same code path the dynamics
!      uses (BVK live-slab elastic + RST GS adatom, eV) — on the printed
!      coordinates,
! and reports the drift of H = T + V against the first frame.
!
! Input: fort.1001 converted to plain numbers, one line per atom,
! 6 floats each (x y z px py pz), NATOMS lines per frame, frames in
! output order. From the run directory:
!   grep -E "^(C|Au) " fort.1001 | awk '{printf "%.8f %.8f %.8f %.8f %.8f %.8f\n", $2,$3,$4,$5,$6,$7}' > frames.txt
!
! Build & run (repo root):
!   ifx -r8 -double-size=64 -i8 -O2 -w \
!     src_TEST/rst_pes.f90 src_TEST/rst_slab_bvk.f90 scripts/conserv_check.f90 \
!     -o conserv_check.e
!   ./conserv_check.e data/rst/nn_weights_rst.txt data/rst/Analytic_Potential.txt \
!     frames.txt 145 12.01 196.97 50
!
! Slab is regenerated with the default GENERATE parameters (a=2.95 A,
! 4 layers, 6x6) — must match the slab the run used. The constant
! VZERO/EHF_baseline offsets cancel in the drift.
!******************************************************************************
program conserv_check
   use rst_slab_bvk
   implicit none
   integer :: natom, nip, i, u, ios, iframe
   real(8) :: m1, m2, am, t_kin, v, h, h0, drift, maxabs, step
   real(8), parameter :: C1 = 0.04184d0, EVKCAL = 23.0605d0
   real(8), allocatable :: q(:), p(:), g(:)
   character(len=512) :: wfile, bfile, ffile, cbuf

   call get_command_argument(1, wfile)
   call get_command_argument(2, bfile)
   call get_command_argument(3, ffile)
   call get_command_argument(4, cbuf); read(cbuf,*) natom
   call get_command_argument(5, cbuf); read(cbuf,*) m1
   call get_command_argument(6, cbuf); read(cbuf,*) m2
   cbuf = '50'
   if (command_argument_count() >= 7) then
      call get_command_argument(7, cbuf)
   end if
   read(cbuf,*) nip

   call rst_bvk_read_params(trim(bfile))
   call rst_bvk_load_weights(trim(wfile))
   call rst_bvk_generate_slab(2.95d0, 4, 6)
   if (rst_bvk_nslab() /= natom - 1) then
      write(6,*) 'ERROR: regenerated slab has ', rst_bvk_nslab(), &
                 ' atoms but NATOMS-1 = ', natom - 1
      stop
   end if

   allocate(q(3*natom), p(3*natom), g(3*natom))
   open(newunit=u, file=trim(ffile), status='old', iostat=ios)
   if (ios /= 0) then
      write(6,*) 'ERROR: cannot open frames file: ', trim(ffile)
      stop
   end if

   h0 = 0.0d0
   maxabs = 0.0d0
   iframe = 0
   write(*,'(A)') '  frame    step       T+V (eV)      drift (eV)'
   do
      do i = 1, natom
         read(u,*,iostat=ios) q(3*i-2), q(3*i-1), q(3*i), &
                              p(3*i-2), p(3*i-1), p(3*i)
         if (ios /= 0) exit
      end do
      if (ios /= 0) exit
      iframe = iframe + 1

      t_kin = 0.0d0
      do i = 1, natom
         am = m2
         if (i == 1) am = m1
         t_kin = t_kin + (p(3*i-2)**2 + p(3*i-1)**2 + p(3*i)**2)/(2.0d0*am)
      end do
      t_kin = t_kin/(C1*EVKCAL)          ! code energy -> eV

      call rst_total_vg(natom, q, v, g)  ! eV, incl. EHF baseline
      h = t_kin + v
      if (iframe == 1) h0 = h
      drift = h - h0
      if (abs(drift) > maxabs) maxabs = abs(drift)
      step = dble((iframe - 1)*nip)
      write(*,'(I7,F10.1,F14.8,ES13.3)') iframe, step, h, drift
   end do
   write(*,'(A,I0,A,ES10.3)') ' frames: ', iframe, &
      '  max |H drift| (eV) = ', maxabs
end program conserv_check
