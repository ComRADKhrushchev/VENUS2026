!=====================================================================
! ens_cau_tab.f90 - the cau_tab scattering ensemble: outgoing kinetic
!   energy distributions for the two DFT lineages, one process (the
!   table loader is re-loadable: run lineage A, swap the table, run B)
! Design:
!   Initial conditions mirror the production recipes
!   (smoke_cau_peri/real_conditions/input_{E1D,E3P}_k1.txt): C launched
!   from 6.0 A at normal incidence aimed at the top site with the
!   impact parameter drawn on the b <= 0.3 A disc (B_MAX), beam energy
!   17.526 kcal/mol (0.76 eV, 1D seed) or 33.6683 kcal/mol (1.46 eV,
!   GS/3P seed); the surface oscillator (S-Y pair, k = 1 eV/A^2) is
!   thermally sampled at T_OSC = 300 K (Y displacement sigma_x =
!   sqrt(kT/k) around its 0.9 A offset, momenta Maxwell, both atoms)
!   - the DIST_SCHEME = BOLTZMANN,SURFACE_OSCILLATOR spirit; RNG = the
!   repo stream, ISEED = 20260928.
!   Integration: velocity-Verlet dt = 0.5 fs, TDHF Magnus-2 (20-orb
!   bath, 7 eV, 300 K, MSBSTP=4), the rho lifecycle reseeds itself at
!   each trajectory head (the prop time-guard). Exit when the C rises
!   past z = 6.0 after entering (Z_TERM); 4000-step cap (2 ps).
!   Per-trajectory record (ens_<lineage>_<cond>.txt):
!     i  status(1=out,0=trapped)  steps  Ekin_out[eV]  z_min  z_end  n_up  n_dn
!   plus per-leg summary lines. Usage:
!     ens_cau_tab <table_rst> <table_current> <N_per_leg> [iseed] [legsel]
!   legsel: 0 = all four legs (default); 1..4 runs ONE leg as its own
!   process (server parallelization: four processes, one leg each, each
!   with a fresh RNG stream from the same ISEED).
!=====================================================================
program ens_cau_tab
   use cau_tab_pot,  only: cau_tab_load, ct_k => k_surf, ct_kg => k_ghost, ct_q0 => Q0_ghost
   use cau_tab_tdhf, only: cau_tab_tdhf_setup, cau_tab_pes, cau_tab_prop, &
                           cau_tab_col_nup, cau_tab_col_ndn, cau_tab_col_eelec
   use rng,          only: rng_init, rng_u, rng_gauss
   use state,        only: state_t, state_create, state_destroy
   implicit none
   real(8), parameter :: ev_kcal = 23.0605d0, e_conv = 0.04184d0
   real(8) :: au_mass = 196.97d0, c_mass = 12.01d0
   real(8), parameter :: b_max = 0.30d0, z_launch = 6.0d0, z_term = 6.0d0
   real(8), parameter :: kt_ev = 1.9872198404d-3*300.0d0/ev_kcal   ! kT at 300 K [eV]
   integer, parameter :: nstep_cap = 4000
   real(8) :: mass3(3)
   character(len=256) :: tab_rst, tab_cur, arg
   integer :: ntraj, iseed, ilin, icond, i, legsel
   character(len=8) :: linname(2), condname(2)
   character(len=2) :: seeds(2)
   real(8) :: ekins(2)

   ntraj = 200
   iseed = 20260928
   legsel = 0
   call get_command_argument(1, tab_rst)
   call get_command_argument(2, tab_cur)
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg)
      read (arg, *) ntraj
   end if
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg)
      read (arg, *) iseed
   end if
   if (command_argument_count() >= 5) then
      call get_command_argument(5, arg)
      read (arg, *) legsel
   end if
   mass3 = [au_mass, au_mass, c_mass]
   linname = ['rst    ', 'current']
   condname = ['E1D_0.76', 'E3P_1.46']
   seeds = ['1D', 'GS']
   ekins = [17.526d0/ev_kcal, 33.6683d0/ev_kcal]        ! eV
   call rng_init(iseed)
   ct_k = [1.0d0, 1.0d0, 1.0d0]
   ct_kg = [1.0d0, 1.0d0, 1.0d0]
   ct_q0 = [0.0d0, 0.0d0, 0.0d0]

   i = 0
   do ilin = 1, 2
      if (ilin == 1) call cau_tab_load(trim(tab_rst))
      if (ilin == 2) call cau_tab_load(trim(tab_cur))
      do icond = 1, 2
         i = i + 1
         if (legsel == 0 .or. legsel == i) then
            if (legsel > 0) call rng_init(iseed)   ! one leg = own stream from ISEED
            call run_leg(linname(ilin), condname(icond), seeds(icond), ekins(icond))
         end if
      end do
   end do

contains

   subroutine run_leg(lin, cond, seed, ekin)
      character(len=*), intent(in) :: lin, cond, seed
      real(8), intent(in) :: ekin
      character(len=64) :: fname
      type(state_t) :: st
      real(8) :: gg(9), vv, pz0, b, phi, sx, sp_au, eout, zmin, zend
      real(8) :: nup, ndn
      integer :: k, nout, u, status, kused
      logical :: passed
      pz0 = -sqrt(2.0d0*c_mass*ekin*e_conv)
      sx = sqrt(kt_ev/1.0d0)                  ! Y thermal displacement sigma [A] (k=1 eV/A^2)
      sp_au = sqrt(au_mass*kt_ev*e_conv)      ! Au thermal momentum sigma [amu A / 10 fs]
      call cau_tab_tdhf_setup(20, 7.0d0, 300.0d0, seed, 4, 1, 0, 0.05d0)
      write (fname, '(a,a,a,a,a)') 'ens_', trim(lin), '_', trim(cond), '.txt'
      open (newunit=u, file=fname, status='replace', action='write')
      nout = 0
      do i = 1, ntraj
         ! 1. sample: impact parameter on the disc, oscillator thermal state
         b = b_max*sqrt(rng_u())
         phi = 2.0d0*3.14159265358979323846d0*rng_u()
         call state_create(st, 3, 1)
         st%q = 0.0d0
         st%q(4:6) = [0.0d0, 0.0d0, 0.9d0] + sx*[rng_gauss(), rng_gauss(), rng_gauss()]
         st%q(7:9) = [b*cos(phi), b*sin(phi), z_launch]
         st%p = 0.0d0
         st%p(1:3) = sp_au*[rng_gauss(), rng_gauss(), rng_gauss()]   ! anchor S thermal
         st%p(4:6) = sp_au*[rng_gauss(), rng_gauss(), rng_gauss()]   ! ghost Y thermal
         st%p(9) = pz0
         st%t = 0.0d0
         ! 2. integrate (the prop time-guard reseeds rho at each head)
         call cau_tab_pes(st%q, mass3, gg, vv)
         zmin = st%q(9)
         zend = st%q(9)
         status = 0
         kused = nstep_cap
         passed = .false.
         do k = 1, nstep_cap
            call vstep(st, 0.05d0)
            call cau_tab_pes(st%q, mass3, gg, vv)
            zmin = min(zmin, st%q(9))
            zend = st%q(9)
            if (st%q(9) < 5.9d0) passed = .true.
            if (passed .and. st%q(9) > z_term) then
               status = 1
               kused = k
               exit
            end if
         end do
         eout = dot_product(st%p(7:9), st%p(7:9))/(2.0d0*c_mass)/e_conv   ! C KE [eV]
         nup = cau_tab_col_nup(st)
         ndn = cau_tab_col_ndn(st)
         write (u, '(i6,i3,i6,3f10.4,2f8.3)') i, status, kused, eout, zmin, zend, nup, ndn
         if (status == 1) nout = nout + 1
         call state_destroy(st)
      end do
      close (u)
      print '(a,a,a,a,a,i0,a,i0)', 'leg ', trim(lin), ' ', trim(cond), ': scattered ', nout, &
         '/', ntraj
   end subroutine run_leg

   subroutine vstep(st, dt)
      type(state_t), intent(inout) :: st
      real(8), intent(in) :: dt
      real(8) :: gg(9), vv
      integer :: ia
      call cau_tab_pes(st%q, mass3, gg, vv)
      do ia = 1, 3
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      do ia = 1, 3
         st%q(3*ia - 2:3*ia) = st%q(3*ia - 2:3*ia) + st%p(3*ia - 2:3*ia)/mass3(ia)*dt
      end do
      call cau_tab_pes(st%q, mass3, gg, vv)
      do ia = 1, 3
         st%p(3*ia - 2:3*ia) = st%p(3*ia - 2:3*ia) - gg(3*ia - 2:3*ia)*ev_kcal*e_conv*0.5d0*dt
      end do
      st%t = st%t + dt
      call cau_tab_prop(st, dt)
   end subroutine vstep

end program ens_cau_tab
