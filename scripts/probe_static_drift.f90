!=====================================================================
! probe_static_drift.f90 - fixed-geometry TDHF energy-conservation probe
!   (the A-group discriminator of the 2026-09-28 energy-debug plan)
! Design:
!   Holds q frozen (no nuclear motion), evolves rho by the member's own
!   Magnus-2 substep loop, and records d[Tr(rho*Hcore)+E_mf] per step.
!   By the mean-field theorem this MUST stay zero (Tr(rho_dot F) = 0 for
!   the unitary flow) up to the Magnus discretization error - a secular
!   drift that does not collapse with MSBSTP is a structural F-vs-E_mf
!   inconsistency; one that scales ~ 1/msbstp^2 is discretization.
!   Command line: z_near[1|0] msbstp nsteps dt  (defaults 1 4 300 0.05)
! Compile (scripts/, standalone against the container like check_cau_peri):
!   gfortran -O2 -o probe_static_drift.exe ../src/globals/consts.f90
!     ../src/globals/state.f90 ../src/utils/linalg.f90 ../src/utils/rng.f90
!     ../systems/inactive/cau_peri/surface_potentials.f90
!     ../systems/inactive/cau_peri/cau_peri_pot.f90
!     ../systems/inactive/cau_peri/cau_peri_tdhf.f90
!     ../systems/inactive/cau_peri/cau_peri_box.f90 probe_static_drift.f90
!=====================================================================
program probe_static_drift
   use consts,       only: e_conv
   use rng,          only: rng_init
   use state,        only: state_t, state_create, state_destroy
   use cau_peri_pot, only: cau_peri_load
   use cau_peri_tdhf, only: cau_peri_tdhf_setup, cau_peri_prop, cau_peri_pes, cau_peri_col_eelec, &
                            cau_peri_col_emf, rho_ref_out
   implicit none
   type(state_t) :: sta
   complex(8), allocatable :: r(:, :, :)
   real(8) :: q(9), mass(3), dt
   real(8) :: e0, emax, esum, e, eprev, drho, drhomax, eel0, emf0
   integer :: msbstp, nsteps, k, narg, u
   character(len=16) :: arg
   logical :: near

   call cau_peri_load()
   if (command_argument_count() >= 5) then
      call cau_peri_load('nn_weights.txt')
   end if
   near = .true.; msbstp = 4; nsteps = 300; dt = 0.05d0
   narg = command_argument_count()
   if (narg >= 1) then; call get_command_argument(1, arg); read (arg, *) k; near = (k /= 0); end if
   if (narg >= 2) then; call get_command_argument(2, arg); read (arg, *) msbstp; end if
   if (narg >= 3) then; call get_command_argument(3, arg); read (arg, *) nsteps; end if
   if (narg >= 4) then; call get_command_argument(4, arg); read (arg, *) dt; end if

   call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, 'GS', msbstp, 0, 0.05d0)
   mass = 196.97d0; mass(3) = 12.01d0
   q = 0.0d0
   q(6) = 0.9d0                               ! ghost Y at its offset
   if (near) then
      q(9) = 2.2d0                            ! engine layout [S, Y, C]: C near-collision
   else
      q(9) = 6.0d0
   end if
   call state_create(sta, 3, 1)
   sta%q = q
   sta%t = 0.0d0
   call rng_init(20260928)                    ! deterministic (PERTURB_MODE=0 draws nothing)
   block
      real(8) :: g9(9), vv
      call cau_peri_pes(q, mass, g9, vv)      ! establish hcore BEFORE the E baseline
   end block

   open (newunit=u, file='static_drift.txt', status='replace')
   write (u, '(a)') '# step t E_sum dE_step rmin dE_elec dE_mf'
   eprev = cau_peri_col_eelec(sta) + cau_peri_col_emf(sta)
   eel0 = cau_peri_col_eelec(sta)
   emf0 = cau_peri_col_emf(sta)
   e0 = eprev
   call rho_ref_out(r)
   emax = 0.0d0; esum = 0.0d0; drhomax = 0.0d0
   do k = 1, nsteps
      sta%t = k*dt
      call cau_peri_prop(sta, dt)
      e = cau_peri_col_eelec(sta) + cau_peri_col_emf(sta)
      drho = 0.0d0
      block
         complex(8), allocatable :: r2(:, :, :)
         integer :: i, j, s2
         call rho_ref_out(r2)
         do s2 = 1, 2
            do i = 1, size(r, 1)
               do j = 1, size(r, 1)
                  drho = max(drho, abs(r2(i, j, s2) - r(i, j, s2)))
               end do
            end do
         end do
         call move_alloc(r2, r)
      end block
      drhomax = max(drhomax, drho)
      emax = max(emax, abs(e - e0))
      esum = e - e0
      block
         real(8) :: rmin
         integer :: i3, j3, s3
         rmin = 1.0d0
         do s3 = 1, 2
            do i3 = 1, 2
               rmin = min(rmin, real(r(i3, i3, s3)))
            end do
         end do
         write (u, '(i6,f9.3,f14.8,es12.3,f10.5,2f14.8)') k, sta%t, e, e - eprev, rmin, &
         cau_peri_col_eelec(sta) - eel0, cau_peri_col_emf(sta) - emf0
      end block
      eprev = e
   end do
   close (u)
   print '(a,l1,a,i0,a,i0,a,f6.3,a)', 'probe_static_drift: near=', near, ' msbstp=', msbstp, &
      ' nsteps=', nsteps, ' dt=', dt, ' [10fs]'
   print '(a,es12.4,a)', '  net dE = ', esum, ' eV; max |dE| along = '
   print '(a,es12.4)', '  max |dE| = ', emax
   print '(a,es12.4)', '  max per-step |drho| = ', drhomax
   call state_destroy(sta)
end program probe_static_drift
