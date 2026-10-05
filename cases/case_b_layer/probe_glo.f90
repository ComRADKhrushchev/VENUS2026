!=====================================================================
! probe_glo.f90 - glo_target ON-CHAIN loading-consistency criteria
!                 program (criteria defined FIRST in this header; this
!                 program executes them). Companion run material of
!                 cases/case_b_layer/glo (the case README documents the
!                 two-stage calibration and the anchors).
! Compile (from src/, against the built gfortran object set - driver.o
! excluded, the driver is the only other program unit; objects live in
! per-source subfolders of build/gfortran, .mod files at its root):
!   gfortran -O2 -Ibuild/gfortran -o build/gfortran/probe_glo.exe \
!     ../cases/case_b_layer/probe_glo.f90 \
!     $(find build/gfortran -name '*.o' ! -name 'driver.o')
! Run (cwd = the glo run directory - reads input_qct.txt, system/,
! reference_pars.txt from cwd, exactly like the driver):
!   cd cases/case_b_layer/glo && ../../../src/build/gfortran/probe_glo.exe
! Approach: unlike the mechanism-level check (synthetic parameters, the
! input grammar out of the loop), this probe drives the REAL assembly
! chain - read_input -> sysdef_load -> list_atoms_load -> registry_gen_all ->
! registry_gen_init (the driver's assembly steps 1,3,4,6,7; no rng and no
! state needed, nothing is sampled) - and then compares the member's
! loaded bath data against HAND values: the input-file literals and the
! IUPAC masses written below (a route independent of the atom list's
! symbol->mass parse). Exit 0 = all checks green; exit 1 = any failure
! (each failure prints check tag, got-vs-hand values, relative error).
! Criteria:
!   (a) parameter loading through the real chain: the t_glo/gamma/w
!       getters return the input literals T_GLO/GAMMA_GLO/W_GHOST
!       (buffered by read_input, pulled by the reg_dist_glo_target seam),
!       and dt_loaded returns the propagator step DT the member was
!       calibrated against - relative 1e-12 each
!   (b) carrier binding: the bound fragment is fragment 1, its atom
!       count is 2 (H+Cl), the ghost-state size is 3*nat = 6 (one ghost
!       oscillator per carrier atom, isotropic three components)
!   (c) list_atoms-mass route: list_atoms%mass of the carrier atoms equals the
!       hand IUPAC masses 1.008 / 35.45 amu (the xyz symbol ->
!       elem_mass parse chain), relative 1e-9
!   (d) fluctuation-dissipation algebra ON-CHAIN: for every carrier
!       atom sigma_noise(i)^2 = 2*gamma*m_i*kb*t_glo/dt with the HAND
!       masses and input literals (unequal masses make the per-atom
!       check sharp: H 0.644093, Cl 3.819673) - relative 1e-12
!   (e) stationary-width algebra ON-CHAIN: sigma_p(i)^2 = m_i*kb*t_glo
!       and sigma_q(i)^2 = kb*t_glo/(m_i*w_ghost^2) with the same hand
!       masses - relative 1e-12
!=====================================================================
program probe_glo
   use consts,        only: kb_code
   use input,         only: read_input
   use config,        only: reactants
   use config_atoms, only: list_atoms, list_atoms_load
   use sysdef_interface, only: sysdef_load
   use registry_gen,  only: registry_gen_all, registry_gen_init, registry_gen_keys
   use dist_glo_target, only: glo_target_t_glo, glo_target_gamma, glo_target_w, &
                              glo_target_dt_loaded, glo_target_frag, glo_target_nat, &
                              glo_target_ghost_size, glo_target_sigma_noise, &
                              glo_target_sigma_p, glo_target_sigma_q
   implicit none
   ! hand side: the input-file literals (glo/input_qct.txt) and the IUPAC
   ! masses (the specio table values, written here as literals - a route
   ! independent of the atom list's symbol->mass parse)
   real(8), parameter :: t_hand = 550.0d0      ! T_GLO [K]
   real(8), parameter :: g_hand = 0.09d0       ! GAMMA_GLO [1/(10 fs)]
   real(8), parameter :: w_hand = 0.0021d0     ! W_GHOST [1/(10 fs)]
   real(8), parameter :: dt_hand = 0.02d0      ! DT [10 fs]
   real(8), parameter :: m_hand(2) = [1.008d0, 35.45d0]  ! H, Cl [amu]
   integer :: i, nfail
   real(8) :: got, hand, rel
   nfail = 0

   ! the real loading chain (driver assembly steps 1,3,4,6,7; nothing sampled)
   call registry_gen_keys()                    ! the slot handover (the driver runs keys first)
   call read_input('input_qct.txt')
   call sysdef_load(reactants%system_dir)
   call list_atoms_load()
   call registry_gen_all()
   call registry_gen_init()

   ! (a) parameter loading through the real chain
   call chk_val('(a) t_glo getter vs T_GLO literal', glo_target_t_glo(), t_hand, 1.0d-12, nfail)
   call chk_val('(a) gamma getter vs GAMMA_GLO literal', glo_target_gamma(), g_hand, 1.0d-12, nfail)
   call chk_val('(a) w getter vs W_GHOST literal', glo_target_w(), w_hand, 1.0d-12, nfail)
   call chk_val('(a) dt_loaded vs DT literal', glo_target_dt_loaded(), dt_hand, 1.0d-12, nfail)

   ! (b) carrier binding
   if (glo_target_frag() /= 1) then
      write (6, '(a,i0)') 'FAIL (b) bound fragment expected 1, got ', glo_target_frag()
      nfail = nfail + 1
   else
      write (6, '(a)') 'PASS (b) bound fragment = 1'
   end if
   if (glo_target_nat() /= 2 .or. glo_target_ghost_size() /= 6) then
      write (6, '(a,2i0)') 'FAIL (b) nat/ghost_size expected 2/6, got ', &
         glo_target_nat(), glo_target_ghost_size()
      nfail = nfail + 1
   else
      write (6, '(a)') 'PASS (b) nat = 2 (H+Cl), ghost_size = 6'
   end if

   do i = 1, glo_target_nat()
      ! (c) list_atoms-mass route vs hand IUPAC masses
      got = list_atoms%mass(list_atoms%frag(glo_target_frag())%list(i))
      call chk_val('(c) list_atoms mass vs IUPAC hand, atom', got, m_hand(i), 1.0d-9, nfail, i)
      ! (d) fluctuation-dissipation algebra on-chain
      hand = sqrt(2.0d0*g_hand*m_hand(i)*kb_code*t_hand/dt_hand)
      call chk_val('(d) sigma_noise vs sqrt(2*gamma*m*kb*T/dt), atom', &
         glo_target_sigma_noise(i), hand, 1.0d-12, nfail, i)
      ! (e) stationary-width algebra on-chain
      hand = sqrt(m_hand(i)*kb_code*t_hand)
      call chk_val('(e) sigma_p vs sqrt(m*kb*T), atom', &
         glo_target_sigma_p(i), hand, 1.0d-12, nfail, i)
      hand = sqrt(kb_code*t_hand/(m_hand(i)*w_hand**2))
      call chk_val('(e) sigma_q vs sqrt(kb*T/(m*w^2)), atom', &
         glo_target_sigma_q(i), hand, 1.0d-12, nfail, i)
   end do

   if (nfail > 0) then
      write (6, '(a,i0)') 'probe_glo: FAILED checks = ', nfail
      stop 1
   end if
   write (6, '(a)') 'probe_glo: all checks passed (loading chain consistent with FDT algebra)'

contains

   ! value comparator with optional atom tag: relative-tolerance check,
   ! PASS/FAIL line with got-vs-hand and relative error
   subroutine chk_val(tag, got_v, hand_v, tol, nf, atom)
      character(len=*), intent(in) :: tag   ! check tag (printed)
      real(8), intent(in) :: got_v          ! member/getter value
      real(8), intent(in) :: hand_v         ! hand-computed value
      real(8), intent(in) :: tol            ! relative tolerance
      integer, intent(inout) :: nf          ! failure counter
      integer, intent(in), optional :: atom ! carrier-atom tag
      rel = abs(got_v - hand_v)/abs(hand_v)
      if (rel <= tol) then
         if (present(atom)) then
            write (6, '(a,a,i0,a,es9.2)') 'PASS ', trim(tag), atom, ' rel = ', rel
         else
            write (6, '(a,a,a,es9.2)') 'PASS ', trim(tag), ' rel = ', rel
         end if
      else
         if (present(atom)) then
            write (6, '(a,a,i0,a,2es16.8,a,es9.2)') 'FAIL ', trim(tag), atom, &
               ' got/hand = ', got_v, hand_v, ' rel = ', rel
         else
            write (6, '(a,a,a,2es16.8,a,es9.2)') 'FAIL ', trim(tag), &
               ' got/hand = ', got_v, hand_v, ' rel = ', rel
         end if
         nf = nf + 1
      end if
   end subroutine chk_val

end program probe_glo
