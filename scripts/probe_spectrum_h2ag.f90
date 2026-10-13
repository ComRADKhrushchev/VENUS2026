!=====================================================================
! probe_spectrum_h2ag.f90 - one-shot diagnostic: the H2/Ag(111) NN's
!   LOCAL curvature frequency at the case's buffered placement
!   (COMPUTE-mode probe value of the 2026-10-05 spectrum-source test
!   batch - docs/plans/2026-10-05-spectrum-source.md §8; empirical
!   record, not a criterion program)
! Usage (from src/, with the h2ag111PES container dragged in and built;
!   run from a directory holding weights-h2ag.txt / biases-h2ag.txt):
!   ../scripts/probe_spectrum_h2ag.exe
! Compile:
!   gfortran -O2 -Jbuild/gfortran -o ../scripts/probe_spectrum_h2ag.exe &
!     ../scripts/probe_spectrum_h2ag.f90 build/gfortran/sys/h2ag111PES/h2ag111_box.o &
!     build/gfortran/sys/h2ag111PES/h2ag111_pot.o build/gfortran/interface/force_interface.o &
!     build/gfortran/utils/hessian.o build/gfortran/utils/spectrum.o &
!     build/gfortran/utils/linalg.o build/gfortran/globals/consts.o &
!     build/gfortran/globals/config.o build/gfortran/globals/config_atoms.o &
!     build/gfortran/utils/specio.o
!=====================================================================
program probe_spectrum_h2ag
   use h2ag111_box,      only: h2ag111_params_t, h2ag111_load, h2ag111_pes
   use force_interface,  only: container_bind_pes, container_probe_force
   use input,            only: read_input
   use config,           only: reactants
   use sysdef_interface, only: sysdef_load
   use config_atoms,     only: list_atoms_load
   use hessian,          only: hessian_cd
   use spectrum,         only: spectrum_modes, w_to_wvn
   use specio,           only: elem_mass
   implicit none
   type(h2ag111_params_t) :: p
   real(8) :: qz(6), hess(6, 6), mass(2)
   real(8), allocatable :: w(:), c(:,:)
   call read_input('input_qct.txt')           ! the grammar stages SCATTERING/CELL (the surface window)
   call sysdef_load(reactants%system_dir)     ! stage the case's system definition
   call list_atoms_load()                     ! assemble the atom list (2 atoms)
   call h2ag111_load(p)                       ! reads weights/biases from cwd
   call container_bind_pes(h2ag111_pes)
   qz = [0.0d0, -0.3707d0, 3.0d0, 0.0d0, 0.3707d0, 3.0d0]   ! H2.poscar buffered geometry
   mass = [elem_mass('H'), elem_mass('H')]
   call hessian_cd(wrap, qz, hess)
   call spectrum_modes(hess, mass, qz, w, c)
   write (*, '(a,es20.10,a)') 'h2ag111 NN local-curvature nu = ', w_to_wvn(w(1)), ' cm^-1'
   write (*, '(a)') 'literature gas-phase H2 w_e = 4401.21 cm^-1 (Herzberg)'
contains
   subroutine wrap(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      call container_probe_force(q, f)
   end subroutine wrap
end program probe_spectrum_h2ag
