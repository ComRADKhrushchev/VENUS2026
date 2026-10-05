!=====================================================================
! probe_hub_audit.f90 - per-channel Hubbard energy/force audit for cau_peri
!   (verifications 2+3 of the 2026-09-29 request)
! Design:
!   mode 'fd'     - fixed near geometry, GS seed evolved 100 steps (real
!                   off-diagonal rho), then for each channel (impurity
!                   on-site h, impurity-bath coupling, mean-field U/J)
!                   compare the CONTAINER'S kernel pieces (assembled from
!                   kernel_caches_out / hcore_out / hub_vals_out) against
!                   central differences of that channel's energy part.
!   mode 'budget' - synthetic 3P beam trajectory (kick-drift-kick verlet
!                   on the composed force + the member's rho evolution),
!                   per step accumulating each force channel's work and
!                   each energy part's change; the closing table is the
!                   additive decomposition of the total drift.
! Compile (scripts/, standalone against the container):
!   gfortran -O2 -o probe_hub_audit.exe ../src/globals/consts.f90
!     ../src/globals/state.f90 ../src/utils/linalg.f90 ../src/utils/rng.f90
!     ../systems/inactive/cau_peri/surface_potentials.f90
!     ../systems/inactive/cau_peri/cau_peri_pot.f90
!     ../systems/inactive/cau_peri/cau_peri_tdhf.f90
!     ../systems/inactive/cau_peri/cau_peri_box.f90 probe_hub_audit.f90
!=====================================================================
program probe_hub_audit
   use consts,        only: e_conv
   use rng,           only: rng_init
   use state,         only: state_t, state_create, state_destroy
   use cau_peri_pot,  only: cau_peri_load, cau_peri_ground, cau_peri_hubbard_val
   use surface_potentials, only: k_surf, k_ghost, q0_ghost
   use cau_peri_tdhf, only: cau_peri_tdhf_setup, cau_peri_prop, cau_peri_pes, &
                            cau_peri_col_eelec, cau_peri_col_emf, &
                            cau_peri_col_nup, cau_peri_col_ndn, cau_peri_col_mag, &
                            rho_ref_out, hcore_out, kernel_caches_out, hub_vals_out
   implicit none
   type(state_t) :: sta
   complex(8), allocatable :: r(:, :, :)
   character(len=16) :: mode
   real(8), parameter :: dx = 1.0d-4
   integer :: ud

   character(len=4) :: seedw
   real(8) :: dt_arg, pt_arg, beam_ov, site_x, site_y, kfac
   integer :: ms_arg
   character(len=8) :: a3
   call get_command_argument(3, a3)
   if (trim(a3) == 'nn') then
      call cau_peri_load('nn_weights.txt')
   else
      call cau_peri_load()
   end if
   call get_command_argument(1, mode)
   if (len_trim(mode) == 0) mode = 'fd'
   seedw = 'GS'
   if (trim(mode) == 'drift') call get_command_argument(2, seedw)
   ! resolution override (drift mode): arg4 = dt [10 fs], arg5 = msbstp,
   ! arg6 = relative perturbation of the beam momentum (realization-noise test)
   dt_arg = 0.05d0
   ms_arg = 4
   pt_arg = 0.0d0
   block
      character(len=32) :: sarg
      if (command_argument_count() >= 4) then
         call get_command_argument(4, sarg); read (sarg, *) dt_arg
      end if
      if (command_argument_count() >= 5) then
         call get_command_argument(5, sarg); read (sarg, *) ms_arg
      end if
      if (command_argument_count() >= 6) then
         call get_command_argument(6, sarg); read (sarg, *) pt_arg
      end if
      beam_ov = 0.0d0
      if (command_argument_count() >= 7) then
         call get_command_argument(7, sarg); read (sarg, *) beam_ov
      end if
      site_x = 0.0d0
      site_y = 0.0d0
      if (command_argument_count() >= 8) then
         call get_command_argument(8, sarg); read (sarg, *) site_x
      end if
      if (command_argument_count() >= 9) then
         call get_command_argument(9, sarg); read (sarg, *) site_y
      end if
      kfac = 1.0d0
      if (command_argument_count() >= 10) then
         call get_command_argument(10, sarg); read (sarg, *) kfac
      end if
   end block
   ! substrate-oscillator constants: the probe bypasses the box's parameter
   ! load, and the vendored module declares k without initializers - set the
   ! PRODUCTION defaults (K_SURF/K_GHOST 1 eV/A^2, Q0_GHOST origin) explicitly,
   ! then apply the phonon A/B factor (arg10; 0 = phonon off)
   k_surf = [1.0d0, 1.0d0, 1.0d0]
   k_ghost = [1.0d0, 1.0d0, 1.0d0]
   q0_ghost = [0.0d0, 0.0d0, 0.0d0]
   k_surf = k_surf*kfac
   k_ghost = k_ghost*kfac
   print '(a,f6.3)', 'substrate oscillator k factor = ', kfac
   call cau_peri_tdhf_setup(20, 7.0d0, 300.0d0, seedw, ms_arg, 0, 0.05d0)
   call rng_init(20260929)
   if (trim(mode) == 'fd') then
      call audit_fd()
   else if (trim(mode) == 'cons') then
      call audit_cons()
   else if (trim(mode) == 'e0chk') then
      call e0chk()
   else if (trim(mode) == 'drift') then
      call drift_run(seedw, dt_arg, ms_arg, pt_arg, beam_ov, site_x, site_y)
   else
      call budget()
   end if

contains

   ! drift_run - synthetic beam trajectory dumping the full per-frame
   !                 bookkeeping (t, h, e_tot) for the drift-structure analysis;
   !                 dt/msbstp overridable for resolution scans (fixed 200 [10
   !                 fs] total time = 2 ps regardless of dt)
   subroutine drift_run(seedw, dt, ms, pt, beam_ov, site_x, site_y)
      character(len=4), intent(in) :: seedw
      real(8), intent(in) :: dt, pt, beam_ov, site_x, site_y
      integer, intent(in) :: ms
      real(8) :: q(9), p(9), m3(3), mdof(9), g(9), v, tkin, e0t, e_beam
      real(8) :: e0v, g0(9)
      real(8), allocatable :: bath0(:, :)
      character(len=4) :: seedf
      character(len=96) :: fname
      integer :: k, step, nsteps, i2, kk
      m3 = [196.97d0, 196.97d0, 12.01d0]
      do k = 1, 3
         mdof(3*k - 2:3*k) = m3(k)
      end do
      q = 0.0d0
      q(6) = 0.0d0                          ! Y at the ghost equilibrium (origin) - the
                                            ! oscillator starts relaxed, no stored strain
      q(1) = 0.15d0; q(9) = 6.0d0
      q(7) = q(1) + site_x; q(8) = site_y
      p = 0.0d0
      e_beam = 17.526d0
      if (seedw == 'GS' .or. seedw == '1DH') e_beam = 33.6683d0
      if (beam_ov > 0.0d0) e_beam = beam_ov
      p(9) = -sqrt(2.0d0*(e_beam*e_conv)*12.01d0)*(1.0d0 + pt)
      call state_create(sta, 3, 1)
      sta%q = q
      sta%t = 0.0d0
      seedf = seedw
      if (seedf == '1DH') seedf = '1D'
      nsteps = nint(200.0d0/dt)
      write (fname, '(a,a,a,i0,a,f6.4,a,es8.1,a,f5.2,a,f5.2,a,f6.1,a,f3.1,a)') 'drift_', trim(seedf), '_m', ms, '_d', dt, &
         '_p', pt, '_s', site_x, '_', site_y, '_e', beam_ov, '_k', kfac, '.txt'
      open (newunit=ud, file=fname, status='replace')
      write (ud, '(a)') '# step t10fs h_A t_kin V e_tot n_up n_dn E0 E_elec E_mf x_rel y_rel e_osc t_c'
      call cau_peri_pes(q, m3, g, v)
      call cau_peri_ground(q, e0v, g0)
      call rho_ref_out(r)
      allocate (bath0(size(r, 1) - 2, 2))
      bath0 = 0.0d0
      do kk = 1, 2
         do i2 = 3, size(r, 1)
            bath0(i2 - 2, kk) = real(r(i2, i2, kk))
         end do
      end do
      tkin = 0.5d0*sum(p*p/mdof)/0.96485d0
      e0t = tkin + v
      write (ud, '(i6,f9.3,f8.3,f12.5,f13.5,f13.5,2f8.4,3f13.5,2f9.3,f10.4,f10.4)') 0, 0.0d0, q(9) - q(3), &
         tkin, v, e0t, col_occ(1), col_occ(2), e0v, cau_peri_col_eelec(sta), &
         cau_peri_col_emf(sta), q(7) - q(1), q(8) - q(2), &
         0.5d0*sum(k_surf*(q(4:6) - q(1:3))**2) + 0.5d0*sum(k_ghost*(q(4:6) - q0_ghost)**2), &
         0.5d0*sum(p(7:9)**2)/12.01d0/0.96485d0
      do step = 1, nsteps
         p = p - 0.5d0*dt*g
         q = q + dt*p/mdof
         call cau_peri_pes(q, m3, g, v)
         p = p - 0.5d0*dt*g
         sta%q = q
         sta%t = step*dt
         call cau_peri_prop(sta, dt)
         call rho_ref_out(r)
         call cau_peri_ground(q, e0v, g0)
         tkin = 0.5d0*sum(p*p/mdof)/0.96485d0
         write (ud, '(i6,f9.3,f8.3,f12.5,f13.5,f13.5,2f8.4,3f13.5,2f9.3,f10.4,f10.4)') step, step*dt, q(9) - q(3), &
            tkin, v, tkin + v, col_occ(1), col_occ(2), e0v, cau_peri_col_eelec(sta), &
            cau_peri_col_emf(sta), q(7) - q(1), q(8) - q(2), &
            0.5d0*sum(k_surf*(q(4:6) - q(1:3))**2) + 0.5d0*sum(k_ghost*(q(4:6) - q0_ghost)**2), &
            0.5d0*sum(p(7:9)**2)/12.01d0/0.96485d0
      end do
      close (ud)
      print '(a,f12.5,a,f6.3,a,i4,a)', trim(seedw)//' drift: final e_tot - e0 = ', (tkin + v) - e0t, &
         '  [dt=', dt, ' msbstp=', ms, ']'
      ! bath dump: final bath-level occupations vs the initial FD seed - the
      ! electron-hole pair content left in the Fermi sea (eps from hcore diag)
      block
         complex(8), allocatable :: rr(:, :, :)
         real(8) :: hc2(22, 22, 2), dn, eh_e, eh_n
         integer :: ib, isb, nb
         call rho_ref_out(rr)
         call hcore_out(hc2)
         nb = min(size(rr, 1), 22)
         eh_e = 0.0d0
         eh_n = 0.0d0
         do isb = 1, 2
            do ib = 3, nb
               dn = real(rr(ib, ib, isb)) - bath0(ib - 2, isb)
               if (abs(dn) > 2.0d-3) then
                  write (*, '(a,i2,a,i1,a,f8.4,a,f8.3)') '  bath level ', ib - 2, ' spin', isb, &
                     ': dn=', dn, '  eps-eV=', hc2(ib, ib, isb)
               end if
               eh_e = eh_e + dn*hc2(ib, ib, isb)
               eh_n = eh_n + abs(dn)
            end do
         end do
         write (*, '(a,f9.4,a,f9.4)') '  bath: sum(dn*eps) =', eh_e, '   sum|dn| =', eh_n
         write (*, '(a,4f8.4)') '  impurity n1up/n2up/n1dn/n2dn =', &
            real(rr(1, 1, 1)), real(rr(2, 2, 1)), real(rr(1, 1, 2)), real(rr(2, 2, 2))
         write (*, '(a,2f8.4)') '  config weight  3P(parallel) / 1D(paired) =', &
            real(rr(1, 1, 1))*real(rr(2, 2, 1)) + real(rr(1, 1, 2))*real(rr(2, 2, 2)), &
            real(rr(1, 1, 1))*real(rr(2, 2, 2)) + real(rr(1, 1, 2))*real(rr(2, 2, 1))
      end block
      call state_destroy(sta)
   end subroutine drift_run

   real(8) function col_occ(ispin)
      integer, intent(in) :: ispin
      col_occ = real(r(1, 1, ispin)) + real(r(2, 2, ispin))
   end function col_occ

   ! channel energies at arbitrary q with the frozen module rho
   subroutine parts_at(q, e_h, e_coup, e_mf)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: e_h, e_coup, e_mf
      real(8) :: g(9), v, hcm(22, 22, 2)
      integer :: is, i, b
      call cau_peri_pes(q, [196.97d0, 196.97d0, 12.01d0], g, v)
      call hcore_out(hcm)
      e_h = 0.0d0
      e_coup = 0.0d0
      do is = 1, 2
         do i = 1, 2
            e_h = e_h + real(r(i, i, is))*hcm(i, i, is)
         end do
         do b = 3, 22
            do i = 1, 2
               e_coup = e_coup + 2.0d0*real(r(i, b, is))*hcm(i, b, is)
            end do
         end do
      end do
      e_mf = cau_peri_col_emf(sta)
   end subroutine parts_at

   ! the container's per-channel kernel pieces at the CURRENT pes caches
   subroutine kernels_now(gh, gcp, gmf, nvnow)
      real(8), intent(out) :: gh(9), gcp(9), gmf(9), nvnow(2, 2)
      real(8) :: dh(2, 9), du(9), dj(9), dvcp(9)
      real(8) :: hcm(22, 22, 2), h2v(2), uv, jv, vgv
      integer :: k, is, i, b
      call kernel_caches_out(dh, du, dj, dvcp)
      call hcore_out(hcm)
      call hub_vals_out(h2v, uv, jv, vgv)
      do is = 1, 2
         do i = 1, 2
            nvnow(i, is) = real(r(i, i, is))
         end do
      end do
      gh = 0.0d0; gcp = 0.0d0; gmf = 0.0d0
      do k = 1, 9
         do is = 1, 2
            do i = 1, 2
               gh(k) = gh(k) + nvnow(i, is)*dh(i, k)
            end do
            if (abs(vgv) > 1.0d-12) then
               do b = 3, 22
                  do i = 1, 2
                     gcp(k) = gcp(k) + 2.0d0*real(r(i, b, is))*(hcm(1, b, is)/vgv)*dvcp(k)
                  end do
               end do
            end if
         end do
         ! E_mf = U*S1*S2 + (U-J)*P, S = per-spin sums, P = nv11*nv21 + nv12*nv22
         ! (occupations are q-independent: dE_mf/dq = du*S1*S2 + (du-dj)*P)
         gmf(k) = du(k)*(nvnow(1, 1) + nvnow(2, 1))*(nvnow(1, 2) + nvnow(2, 2)) &
                + (du(k) - dj(k))*(nvnow(1, 1)*nvnow(2, 1) + nvnow(1, 2)*nvnow(2, 2))
      end do
   end subroutine kernels_now

   subroutine audit_fd()
      real(8) :: q(9), g(9), v
      real(8) :: e_h, e_cp, e_mf, ehp, ecp_p, emfp, ehm, ecp_m, emfm
      real(8) :: gh(9), gcp(9), gmf(9), nvnow(2, 2)
      integer :: k, kk
      q = 0.0d0
      q(6) = 0.9d0
      q(1) = 0.37d0; q(2) = 0.21d0; q(9) = 2.2d0
      call state_create(sta, 3, 1)
      sta%q = q
      do kk = 1, 100
         sta%t = kk*0.05d0
         call cau_peri_prop(sta, 0.05d0)
      end do
      call rho_ref_out(r)
      call cau_peri_pes(q, [196.97d0, 196.97d0, 12.01d0], g, v)
      call kernels_now(gh, gcp, gmf, nvnow)
      print '(a)', '== per-channel FD audit (evolved rho; dx=1e-4 A) =='
      print '(a)', 'k | ana_h        fd_h          | ana_coup     fd_coup       | ana_mf       fd_mf'
      do k = 1, 9
         q(k) = q(k) + dx
         call parts_at(q, ehp, ecp_p, emfp)
         q(k) = q(k) - 2*dx
         call parts_at(q, ehm, ecp_m, emfm)
         q(k) = q(k) + dx
         write (*, '(i1,1x,3(es12.4,1x,es12.4,4x))') k, gh(k), (ehp - ehm)/(2*dx), &
               gcp(k), (ecp_p - ecp_m)/(2*dx), gmf(k), (emfp - emfm)/(2*dx)
      end do
      call state_destroy(sta)
   end subroutine audit_fd

   ! audit_cons - full composed V-vs-g FD consistency scan along the C-S
   !                 approach coordinate at the probe trajectory's lateral
   !                 geometry (C-S offset x = -0.15, y = 0); a localized
   !                 mismatch band is the signature of an assembly seam
   subroutine audit_cons()
      real(8) :: q(9), g(9), v, v_p, v_m, fd, err, werr, hs, gmax
      real(8) :: e0v, g0(9), ehp, emfm
      integer :: k, kk, ih
      real(8) :: m3(3) = [196.97d0, 196.97d0, 12.01d0]
      q = 0.0d0
      q(6) = 0.9d0
      q(1) = 0.15d0; q(9) = 2.2d0
      call state_create(sta, 3, 1)
      sta%q = q
      do kk = 1, 100
         sta%t = kk*0.05d0
         call cau_peri_prop(sta, 0.05d0)
      end do
      print '(a)', '== composed V-vs-g consistency scan (dx=1e-4 A, evolved rho) =='
      print '(a)', 'h[A]    max_k |dV/dq_k(FD) - g_k| [eV/A]'
      gmax = 0.0d0
      do ih = 0, 96
         hs = 2.0d0 + 0.05d0*ih
         q = 0.0d0
         q(6) = 0.9d0
         q(1) = 0.15d0
         q(9) = hs
         call cau_peri_pes(q, m3, g, v)
         werr = 0.0d0
         do k = 1, 9
            q(k) = q(k) + dx
            call cau_peri_pes(q, m3, g, v_p)
            q(k) = q(k) - 2*dx
            call cau_peri_pes(q, m3, g, v_m)
            q(k) = q(k) + dx
            fd = (v_p - v_m)/(2*dx)
            err = abs(fd - g(k))
            if (err > werr) werr = err
         end do
         if (werr > gmax) gmax = werr
         if (werr > 1.0d-6) write (*, '(f6.2, 2x, es12.4)') hs, werr
      end do
      print '(a, es12.4)', 'scan max mismatch [eV/A] = ', gmax
      ! per-component breakdown at h = 3.2 (the early-event height)
      q = 0.0d0
      q(6) = 0.9d0
      q(1) = 0.15d0
      q(9) = 3.2d0
      call cau_peri_pes(q, m3, g, v)
      print '(a)', 'h=3.2 per-component fd - g [eV/A]:'
      do k = 1, 9
         q(k) = q(k) + dx
         call cau_peri_pes(q, m3, g, v_p)
         q(k) = q(k) - 2*dx
         call cau_peri_pes(q, m3, g, v_m)
         q(k) = q(k) + dx
         write (*, '(i2, 2x, 2es14.6)') k, (v_p - v_m)/(2*dx) - g(k), g(k)
      end do
      ! E0-only channel: FD of cau_peri_ground vs its delivered g0
      q = 0.0d0
      q(6) = 0.9d0
      q(1) = 0.15d0
      q(9) = 3.2d0
      call cau_peri_ground(q, e0v, g0)
      print '(a)', 'h=3.2 E0-only per-component fd - g0 [eV/A]:'
      do k = 1, 9
         q(k) = q(k) + dx
         call cau_peri_ground(q, e0v, g0)
         ehp = e0v
         q(k) = q(k) - 2*dx
         call cau_peri_ground(q, e0v, g0)
         emfm = e0v
         q(k) = q(k) + dx
         write (*, '(i2, 2x, 2es14.6)') k, (ehp - emfm)/(2*dx) - g0(k), g0(k)
      end do
      ! purity test: is E0(q) / pes(q) a pure function of q at frozen rho?
      ! interleave two geometries and re-visit - a changed value = call-history
      ! dependence in the vendored evaluator
      block
         real(8) :: qa(9), qb(9), va1, va2, vb1, vb2, ga(9), gb(9)
         real(8) :: ea1, ea2, eb1, eb2
         qa = 0.0d0; qa(6) = 0.9d0; qa(1) = 0.15d0; qa(9) = 3.2d0
         qb = qa; qb(9) = 5.0d0
         call cau_peri_ground(qa, ea1, ga)
         call cau_peri_ground(qb, eb1, gb)
         call cau_peri_ground(qa, ea2, ga)
         print '(a)', 'E0 purity: qa / qb / qa-revisit:'
         print '(3es20.12)', ea1, eb1, ea2
         call cau_peri_pes(qa, m3, ga, va1)
         call cau_peri_pes(qb, m3, gb, vb1)
         call cau_peri_pes(qa, m3, ga, va2)
         print '(a)', 'pes V purity (qa, qb, qa-revisit):'
         print '(3es20.12)', va1, vb1, va2
         call cau_peri_ground(qa, ea1, ga)
         print '(a)', 'E0(qa) after pes calls vs original:'
         print '(2es20.12)', ea1, ea2
      end block
      ! short replay: expose which quantity carries the early-frame wiggles
      block
         real(8) :: q2(9), p2(9), g2(9), vv, ee0, gg0(9), md2(9)
         integer :: st, kk
         do kk = 1, 3
            md2(3*kk - 2:3*kk) = m3(kk)
         end do
         q2 = 0.0d0
         q2(6) = 0.9d0
         q2(1) = 0.15d0; q2(9) = 6.0d0
         p2 = 0.0d0
         p2(9) = -sqrt(2.0d0*(17.526d0*e_conv)*12.01d0)
         call cau_peri_pes(q2, m3, g2, vv)
         print '(a)', 'replay: step  x_C    y_C    z_C     E0       g0z      V        col_eelec'
         do st = 1, 30
            p2 = p2 - 0.5d0*0.05d0*g2
            q2 = q2 + 0.05d0*p2/md2
            call cau_peri_pes(q2, m3, g2, vv)
            p2 = p2 - 0.5d0*0.05d0*g2
            sta%q = q2
            sta%t = st*0.05d0
            call cau_peri_prop(sta, 0.05d0)
            call cau_peri_ground(q2, ee0, gg0)
            write (*, '(i4, 2f7.3, f7.3, 2f9.4, f10.4, f10.4)') st, q2(7), q2(8), q2(9), &
               ee0, gg0(9), vv, cau_peri_col_eelec(sta)
         end do
      end block
      call state_destroy(sta)
   end subroutine audit_cons

   ! e0chk - E0 context test: same geometries evaluated fresh after setup vs
   !                 after a pes/prop warmup, plus the z-tail force scan
   subroutine e0chk()
      real(8) :: q(9), g(9), v, e, gz
      real(8) :: m3(3) = [196.97d0, 196.97d0, 12.01d0]
      integer :: k
      call state_create(sta, 3, 1)
      print '(a)', 'fresh (right after setup):'
      do k = 6, 2, -2
         q = 0.0d0; q(6) = 0.9d0; q(1) = 0.15d0; q(9) = real(k, 8)
         call cau_peri_ground(q, e, g)
         call cau_peri_pes(q, m3, g, v)
         write (*, '(a,f5.1,a,f12.5,a,f10.5,a,f12.5)') 'h=', real(k), '  E0=', e, '  g0z=', g(9), '  V=', v
      end do
      ! warmup: 100 prop steps at h=2.2 (the cons-mode warmup path)
      q = 0.0d0; q(6) = 0.9d0; q(1) = 0.15d0; q(9) = 2.2d0
      sta%q = q
      do k = 1, 100
         sta%t = k*0.05d0
         call cau_peri_prop(sta, 0.05d0)
      end do
      print '(a)', 'after 100 prop steps at h=2.2:'
      do k = 6, 2, -2
         q = 0.0d0; q(6) = 0.9d0; q(1) = 0.15d0; q(9) = real(k, 8)
         call cau_peri_ground(q, e, g)
         call cau_peri_pes(q, m3, g, v)
         write (*, '(a,f5.1,a,f12.5,a,f10.5,a,f12.5)') 'h=', real(k), '  E0=', e, '  g0z=', g(9), '  V=', v
      end do
      ! lateral landscape scan (y=0.10 line, z=2.2): E0 + Hubbard delivery
      block
         real(8) :: ql(9), gl(9), el, vl, hl2(2), ul, jl, vcl, xx
         integer :: ix
         do ix = 0, 10
            xx = 0.295d0*ix
            ql = 0.0d0
            ql(6) = 0.9d0
            ql(1) = 0.15d0
            ql(7) = ql(1) + xx
            ql(8) = 0.10d0
            ql(9) = 2.2d0
            call cau_peri_ground(ql, el, gl)
            call cau_peri_hubbard_val(ql, hl2, ul, jl, vcl)
            write (*, '(a,7f10.4)') 'land ', xx, el, hl2(1), ul, jl, vcl
         end do
      end block
      ! site-attribution check: the true high-symmetry sites must be lateral
      ! stationary points (g0x=g0y~0), lattice-periodic (x+a reproduces the
      ! same site), while a generic point carries a real lateral gradient
      block
         real(8) :: q4(9), g4(9), e4, zs, s1, s2
         integer :: is4, k4
         real(8) :: sx(4), sy(4)
         sx = [0.0d0, 1.475d0, 1.475d0, 0.7d0]
         sy = [0.0d0, 0.0d0, 0.8515d0, 0.4d0]
         do k4 = 1, 2
            zs = merge(2.5d0, 3.2d0, k4 == 1)
            print '(a,f4.1,a)', 'site check z=', zs, '   x      y       E0      g0x     g0y     g0z'
            do is4 = 1, 4
               q4 = 0.0d0
               q4(6) = 0.9d0
               q4(1) = 0.15d0
               q4(7) = q4(1) + sx(is4)
               q4(8) = sy(is4)
               q4(9) = zs
               call cau_peri_ground(q4, e4, g4)
               write (*, '(2f7.3, f10.4, 3f9.4)') sx(is4), sy(is4), e4, g4(7), g4(8), g4(9)
            end do
            q4 = 0.0d0
            q4(6) = 0.9d0
            q4(1) = 0.15d0
            q4(7) = q4(1) + 2.95d0
            q4(8) = 0.0d0
            q4(9) = zs
            call cau_peri_ground(q4, e4, g4)
            print '(a,f10.4,a,2f9.4)', '  top+(a,0): E0=', e4, '  g0xy=', g4(7), g4(8)
         end do
      end block
      ! per-site well-depth scan: min_z E0 at each site vs its far asymptote
      block
         real(8) :: q5(9), g5(9), e5, zs5, mn(3), mzt(3), asym(3)
         integer :: iz5, is5
         real(8) :: sx5(3), sy5(3)
         sx5 = [0.0d0, 1.475d0, 1.475d0]
         sy5 = [0.0d0, 0.0d0, 0.8515d0]
         do is5 = 1, 3
            q5 = 0.0d0
            q5(6) = 0.9d0
            q5(1) = 0.15d0
            q5(7) = q5(1) + sx5(is5)
            q5(8) = sy5(is5)
            q5(9) = 7.0d0
            call cau_peri_ground(q5, asym(is5), g5)
            mn(is5) = 1.0d30
            do iz5 = 12, 60
               zs5 = real(iz5, 8)/10.0d0
               q5(9) = zs5
               call cau_peri_ground(q5, e5, g5)
               if (e5 < mn(is5)) then
                  mn(is5) = e5
                  mzt(is5) = zs5
               end if
            end do
         end do
         print '(a)', 'per-site well:  site        z_min   E0_min    asympt   depth[eV]'
         print '(a,f5.2,f9.4,f10.4,f9.4)', '  top     ', mzt(1), mn(1), asym(1), mn(1) - asym(1)
         print '(a,f5.2,f9.4,f10.4,f9.4)', '  bridge  ', mzt(2), mn(2), asym(2), mn(2) - asym(2)
         print '(a,f5.2,f9.4,f10.4,f9.4)', '  hollow  ', mzt(3), mn(3), asym(3), mn(3) - asym(3)
      end block
      ! fixed-geometry prop test across the approach band: h/U/J/Vcp
      ! delivery + does the impurity occupation run away (charging)
      block
         real(8) :: q3(9), g3(9), v3, h3(2), u3, j3, vc3, o0, o1, hh
         integer :: k3, ih3
         do ih3 = 1, 6
            hh = 4.5d0 + 0.5d0*ih3
            q3 = 0.0d0
            q3(6) = 0.9d0
            q3(1) = 0.15d0
            q3(9) = hh
            call cau_peri_pes(q3, m3, g3, v3)
            call hub_vals_out(h3, u3, j3, vc3)
            o0 = cau_peri_col_nup(sta) + cau_peri_col_ndn(sta)
            do k3 = 1, 100
               sta%q = q3
               sta%t = k3*0.05d0
               call cau_peri_prop(sta, 0.05d0)
            end do
            o1 = cau_peri_col_nup(sta) + cau_peri_col_ndn(sta)
            print '(a,f5.1,a,5f9.3,a,f7.3,a,f7.3)', 'h=', hh, ' h1/h2/U/J/Vcp =', h3(1), h3(2), u3, j3, vc3, &
               '  occ:', o0, '->', o1
         end do
      end block
      print '(a)', 'z-tail scan (x=-0.15,y=0): z[A]  E0[eV]  g0z[eV/A]  gz_total[eV/A]'
      do k = 100, 5, -5
         gz = real(k, 8)/10.0d0
         q = 0.0d0; q(6) = 0.9d0; q(1) = 0.15d0; q(9) = gz
         call cau_peri_ground(q, e, g)
         call cau_peri_pes(q, m3, g, v)
         write (*, '(f6.2, f12.5, 2f11.5)') gz, e, g(9)
      end do
      call state_destroy(sta)
   end subroutine e0chk

   subroutine budget()
      real(8) :: q(9), qp(9), p(9), g(9), v, m3(3), mdof(9)
      real(8) :: e_h, e_cp, e_mf, e_h0, e_cp0, e_mf0, e0v, g0(9), e0_start
      real(8) :: gh(9), gcp(9), gmf(9), nvnow(2, 2)
      real(8) :: wh, wcp, wmf, w0, tkin, dt
      real(8) :: gh_p(9), gcp_p(9), gmf_p(9), first(9)
      integer :: k, step
      dt = 0.05d0
      m3 = [196.97d0, 196.97d0, 12.01d0]
      do k = 1, 3
         mdof(3*k - 2:3*k) = m3(k)
      end do
      q = 0.0d0
      q(6) = 0.0d0                          ! Y at the ghost equilibrium (origin) - the
                                            ! oscillator starts relaxed, no stored strain
      q(1) = 0.15d0; q(9) = 6.0d0
      q(7) = q(1) + site_x; q(8) = site_y
      p = 0.0d0
      p(9) = -sqrt(2.0d0*(33.6683d0*e_conv)*12.01d0)   ! 1.46 eV on C (engine atom 3)
      call state_create(sta, 3, 1)
      sta%q = q
      sta%t = 0.0d0
      call cau_peri_pes(q, m3, g, v)                    ! initial force (kick 1)
      call rho_ref_out(r)
      call parts_at(q, e_h0, e_cp0, e_mf0)
      call cau_peri_ground(q, e0v, g0)
      e0_start = e0v
      wh = 0.0d0; wcp = 0.0d0; wmf = 0.0d0; w0 = 0.0d0
      gh_p = 0.0d0; gcp_p = 0.0d0; gmf_p = 0.0d0; first = 0.0d0
      ! channel forces at the INITIAL point (trapezoid anchor)
      call cau_peri_pes(q, m3, g, v)
      call rho_ref_out(r)
      call kernels_now(gh_p, gcp_p, gmf_p, nvnow)
      first = g - gh_p - gcp_p - gmf_p
      do step = 1, 4000
         p = p - 0.5d0*dt*g
         qp = q
         q = q + dt*p/mdof
         call cau_peri_pes(q, m3, g, v)                 ! endpoint composed force
         p = p - 0.5d0*dt*g
         sta%q = q
         sta%t = step*dt
         call rho_ref_out(r)
         call kernels_now(gh, gcp, gmf, nvnow)          ! SAME rho as the pes force
         do k = 1, 9
            wh = wh - 0.5d0*(gh(k) + gh_p(k))*(q(k) - qp(k))
            wcp = wcp - 0.5d0*(gcp(k) + gcp_p(k))*(q(k) - qp(k))
            wmf = wmf - 0.5d0*(gmf(k) + gmf_p(k))*(q(k) - qp(k))
            w0 = w0 - 0.5d0*((g(k) - gh(k) - gcp(k) - gmf(k)) + first(k))*(q(k) - qp(k))
         end do
         gh_p = gh; gcp_p = gcp; gmf_p = gmf; first = g - gh - gcp - gmf
         if (q(9) - q(3) > 8.0d0) exit            ! scattered: budget to the exit
         call cau_peri_prop(sta, dt)                    ! rho evolves AFTER the work
      end do
      ! final table (energies in eV; works converted /e_conv; per-channel
      ! mismatch = W + dE, zero = the channel conserves; the sum of the
      ! mismatches reproduces the total drift additively)
      call parts_at(q, e_h, e_cp, e_mf)
      call cau_peri_ground(q, e0v, g0)
      tkin = 0.5d0*sum(p*p/mdof)/0.96485d0
      print '(a)', '== per-channel budget over the synthetic 3P trajectory (eV) =='
      print '(a,f14.4)', 'C height at end = ', q(9) - q(3)
      print '(a)', 'channel        W(force)      dE(book)      mismatch (W+dE)'
      print '(a,f14.4,f14.4,f14.4)', 'E0 ground    ', w0, e0v - e0_start, &
         w0 + (e0v - e0_start)
      print '(a,f14.4,f14.4,f14.4)', 'h impurity   ', wh, e_h - e_h0, &
         wh + (e_h - e_h0)
      print '(a,f14.4,f14.4,f14.4)', 'coupling     ', wcp, e_cp - e_cp0, &
         wcp + (e_cp - e_cp0)
      print '(a,f14.4,f14.4,f14.4)', 'mean-field   ', wmf, e_mf - e_mf0, &
         wmf + (e_mf - e_mf0)
      print '(a,f14.4)', 'T_final - T_initial = ', tkin - 1.46d0
      call state_destroy(sta)
   end subroutine budget

end program probe_hub_audit
