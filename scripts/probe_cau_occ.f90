!=====================================================================
! probe_cau_occ.f90 - per-step spin-orbital occupation dump for the cau_peri
!   TDHF member (the four impurity components + bath trace, all trajectories)
! Design:
!   Mirrors the driver's assembly (steps 1-10) and evolution choreography
!   (samp_run -> force_eval -> [prop_step -> elec_prop -> mqc_statistic ->
!   dump -> container_term]) with the recorder replaced by an occ_dump.txt
!   table: traj, step, t[10fs], C height above the anchor [A], the four
!   spin-orbital impurity occupations (o1_up, o2_up, o1_dn, o2_dn = Re
!   rho_11s / rho_22s) and the bath trace per spin. Reads input_qct.txt
!   from the cwd exactly like the driver (run it in a leg's run directory;
!   the cau_peri container must be the active build).
! Compile (from src/, against the danced tree's object set):
!   gfortran -O2 -Jbuild/gfortran -c ../scripts/probe_cau_occ.f90 -o build/ck_occ.o
!   gfortran build/ck_occ.o <tree objects, driver.o excluded> -o build/probe_cau_occ.exe
!=====================================================================
program probe_cau_occ
   use consts,            only: e_conv
   use rng,               only: rng_init
   use control,           only: n_traj, i_seed, max_steps
   use config,            only: reactants, electronic, observables, prop_cfg => propagator
   use config_atoms,      only: list_atoms, list_atoms_load
   use state,             only: state_t, state_create, state_destroy
   use sysdef_interface,  only: sysdef_load
   use input,             only: read_input, input_audit_declared
   use registry_gen,      only: registry_gen_all, registry_gen_init, registry_gen_keys
   use sampler,           only: samp_run, samp_incident_word
   use force_interface,   only: container_init, force_eval, container_term
   use propagator,        only: prop_step
   use elec_interface,    only: elec_prop, mqc_statistic
   use cau_peri_tdhf,     only: rho_ref_out
   implicit none
   type(state_t) :: sta
   complex(8), allocatable :: r(:, :, :)
   real(8) :: e0, dt_eff, t
   integer :: traj_i, step_i, u, ul, ios, ns, i
   logical :: term_flag
   character(len=4) :: seq_mode

   ! 1.-10. the driver's assembly, verbatim order (no recorder selection:
   !        the probe owns its own dump)
   call registry_gen_keys()
   call read_input('input_qct.txt')
   call rng_init(i_seed(1))
   call sysdef_load(reactants%system_dir)
   call list_atoms_load()
   call state_create(sta, size(list_atoms%mass), max(electronic%n_surf, 1))
   call registry_gen_all()
   call registry_gen_init()
   call input_audit_declared()
   call container_init()
   dt_eff = prop_cfg%dt

   open (newunit=u, file='occ_dump.txt', status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write (0, '(a)') 'probe_cau_occ: cannot open occ_dump.txt in the run directory'
      stop 1
   end if
   write (u, '(a)') '# traj step t10fs h_A o1_up o2_up o1_dn o2_dn bath_up bath_dn'
   open (newunit=ul, file='occ_levels.txt', status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write (0, '(a)') 'probe_cau_occ: cannot open occ_levels.txt in the run directory'
      stop 1
   end if
   write (ul, '(a)') '# traj step e_elec e_mf r12_up r12_dn | bath diagonals up then dn (n_orb)'

   seq_mode = 'a'
   if (command_argument_count() > 0) call get_command_argument(1, seq_mode)
   do traj_i = 1, n_traj
      sta%t = 0.0d0
      call samp_run(samp_incident_word(), sta)
      call force_eval(sta, e0)
      call rho_ref_out(r)
      call dump_levels(traj_i, 0)
      ns = size(r, 1)
      call force_eval(sta, e0)
      write (u, '(i6,i6,f8.2,f7.3,f12.4,4f9.5,2f10.5)') traj_i, 0, sta%t, &
         sta%q(3) - sta%q(6), e0, real(r(1, 1, 1)), real(r(2, 2, 1)), &
         real(r(1, 1, 2)), real(r(2, 2, 2)), &
         sum([(real(r(i, i, 1)), i=3, ns)]), sum([(real(r(i, i, 2)), i=3, ns)])
      step_i = 0
      term_flag = .false.
      do while (.not. term_flag .and. step_i < max_steps)
         step_i = step_i + 1
         if (seq_mode == 'c') then
            ! symmetric split: half electronic, full nuclear, half electronic
            call elec_prop(sta, 0.5d0*dt_eff)
            call prop_step(sta)
            call elec_prop(sta, 0.5d0*dt_eff)
         else
            call prop_step(sta)
            call elec_prop(sta, dt_eff)
         end if
         call mqc_statistic(sta)
         call rho_ref_out(r)
         call force_eval(sta, e0)
         write (u, '(i6,i6,f8.2,f7.3,f12.4,4f9.5,2f10.5)') traj_i, step_i, sta%t, &
            sta%q(3) - sta%q(6), e0, real(r(1, 1, 1)), real(r(2, 2, 1)), &
            real(r(1, 1, 2)), real(r(2, 2, 2)), &
            sum([(real(r(i, i, 1)), i=3, ns)]), sum([(real(r(i, i, 2)), i=3, ns)])
         term_flag = container_term(sta)
         call dump_levels(traj_i, step_i)
      end do
      write (*, '(a,i0,a,i0)') 'probe_cau_occ: trajectory ', traj_i, ' ended at step ', step_i
   end do
   close (u)
   close (ul)
   call state_destroy(sta)
   stop 0
contains

   subroutine dump_levels(tr, st)
      use cau_peri_tdhf, only: cau_peri_diag_out, hub_vals_out
      integer, intent(in) :: tr, st
      real(8) :: ee, em, iu, id, mu2
      real(8) :: hv(2), uv, jv, vgv
      call cau_peri_diag_out(ee, em, iu, id, mu2)
      call hub_vals_out(hv, uv, jv, vgv)
      write (ul, '(i6,i6,9f11.5)', advance='no') tr, st, ee, em, &
         real(r(1, 2, 1)), real(r(1, 2, 2)), hv(1), uv, jv, vgv
      do i = 1, ns
         if (i > 2) write (ul, '(f9.5)', advance='no') real(r(i, i, 1))
      end do
      do i = 1, ns
         if (i > 2) write (ul, '(f9.5)', advance='no') real(r(i, i, 2))
      end do
      write (ul, '(a)') ''
   end subroutine dump_levels
end program probe_cau_occ
