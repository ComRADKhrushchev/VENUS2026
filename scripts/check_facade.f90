!=====================================================================
! check_facade.f90 - facade-mechanism criteria program, batches 1-2:
!                    the frozen ws holder association (a)-(d) + the
!                    force_eval fill chain (e)-(g) (criteria defined
!                    FIRST in this header; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod
! files out of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 &
!            globals/config_roster.f90 globals/state.f90 facade/elec_ws.f90 &
!            facade/elec_facade.f90 utils/linalg.f90 facade/force_facade.f90 &
!            ../scripts/check_facade.f90 -o build/check_facade.exe
! All artifacts are written under build/ - the source tree stays clean.
! Posture: criteria (a)-(d) exercise the holder mechanism itself - a
!   degenerate no-op method row registers + binds the name, which seeds
!   the statistics holder and makes the occ_get channel live; the members
!   are never dispatched by that batch. Criteria (e)-(g) additionally
!   assemble a FORCE CHAIN check-side: a local analytic branch container
!   is bound through container_bind_pes (the public assembly face), so
!   every (e)/(g) assertion runs on real force_eval calls; the default
!   binding posture of the tree is untouched (the driver binds nothing).
!   The bound container branches on a check-side fill selector: 1 = the
!   scalar harmonic channel, 2 = a two-surface diabatic channel (linear
!   H_ij + dh_ij + an additive scalar member), 3 = the one-surface
!   diabatic twin of branch 1 (same V and gradient delivered as 1x1
!   constructs, additive channel absent) - every branch is a local
!   procedure of this check; none of them ever enters src/. No loud-stop
!   criteria here, hence no child-process machinery.
! Criteria:
!   (a) handle identity: two ws_handle() calls yield pointers to the SAME
!       target (pointer association comparison), and a field written
!       through one handle is visible through the other handle AND through
!       the ws_get snapshot channel (the alias is onto the holder itself,
!       not onto a copy)
!   (b) copy-out semantics: ws_get returns a snapshot - scribbling the
!       snapshot (a scalar component and an array component) leaves the
!       holder intact, and a fresh ws_get reproduces the holder's
!       pre-scribble content
!   (c) reset semantics: with the holder dirtied (repr set, scalar set,
!       every per-mode array allocated and filled), ws_reset returns it to
!       the DEFINED empty state - a subsequent ws_get shows repr = r_empty,
!       every per-mode array unallocated, the scalar zeroed (the one-step
!       survival contract: representation constructs live one integration
!       step, then die - not merely become unreadable)
!   (d) snapshot-style vs handle-style coexistence: the existing consumer
!       pattern of the scalar force branch (occ_get request -> ws_get
!       snapshot -> compose into the snapshot) runs unchanged and leaves
!       the HOLDER untouched (zero interference), while the handle-style
!       compose (the fill-site posture that writes the holder itself)
!       reads back through the snapshot channel without retroactively
!       perturbing an earlier snapshot. The regression half of this
!       criterion is external: check_propagator (the force_eval scalar
!       branch end-to-end on real forces) and check_driver (the full
!       four-phase loop) re-run green on the rewritten holder.
!       Historical note (d1): when (d) was written the production fill
!       site still composed into a ws_get snapshot; with the force_eval
!       realization the production fill went handle-first - the
!       snapshot-consumer pattern validated here remains a legal
!       consumption form (ws_get stays a public copy-out primitive), but
!       it no longer describes the fill site itself
!   (e) end-to-end fill chain through the handle: with the analytic
!       container bound through the public face, force_eval returns the
!       analytic force component-wise at two geometries (proportionality
!       to the gradient - the conversion constant cancels) and e_tot
!       decomposes as t_kin(p) + v_pot(q) via the composed conversion
!       measured by one calibration call (the eV->internal literal lives
!       only in force_facade - this check NEVER writes it); the container
!       consumed the occupation request (i_want = [1] observed inside the
!       container), and AFTER the call the holder - read back through a
!       fresh ws_handle() alias AND a ws_get snapshot - carries the
!       composed constructs of the returned geometry (repr, the request
!       field, g_add, v_add): the fill reaches the HOLDER itself and
!       SURVIVES the call (the one-step conveyance to the electronic
!       phases)
!   (f) conversion single-site: a recursive scan of every src/ .f90 file
!       (findstr on a Windows host, grep elsewhere) lists the files whose
!       text contains the conversion literal; each listed file is then
!       scanned line by line, and exactly ONE occurrence of the literal
!       sits at EXECUTABLE position (not inside a comment) in the whole
!       tree - the ev_kcal parameter declaration of force_facade.f90
!       (the unit-contract grep requirement, executable face)
!   (g) the multi-surface diabatic fold: with the container delivering a
!       two-surface linear H_ij + dh_ij + an additive scalar member,
!       force_eval returns the closed-form eigen-fold - the occupied
!       eigenvalue's gradient u^T (dH/dq) u with u the algebraic
!       eigenvector of the LOWER eigenvalue (the ascending occupation
!       row) plus the additive kernel, component-wise (rel 1e-12), and
!       e_tot = t_kin + (lambda_occ + v_add) via the measured conversion
!       (rel 1e-12); the delivered h_ij SURVIVES in the holder (the
!       facade diagonalizes a copy - the constructs stay for the
!       electronic phases). Degenerate equivalence: the one-surface diab
!       twin vs the scalar channel at the same state - force and e_tot
!       compared BITWISE (the attempt), falling back to rel 1e-15 with a
!       printed note
!=====================================================================
program check_facade
   use consts,        only: e_conv
   use state,         only: state_t, state_create, state_destroy
   use config_roster, only: roster
   use elec_ws,       only: ws_t, r_empty, r_scalar, r_diab
   use elec_facade,   only: ws_handle, ws_get, ws_reset, occ_get, &
                            method_reg, method_bind
   use force_facade,  only: force_eval, container_bind_pes
   implicit none
   type(ws_t) :: w, w2                     ! snapshot instances (ws_get dummies)
   type(ws_t), pointer :: p1, p2           ! handle aliases
   integer, allocatable :: occ_req(:)      ! the delivery request (occ_get)
   integer :: nfail = 0
   ! ---------- fill-chain section state ((e)-(g)) ----------
   type(state_t) :: sa                     ! the force-chain state (2 atoms / 6 dofs)
   integer :: nfail0                       ! per-criterion failure watermark
   integer :: m                            ! fold loop index
   integer :: fill_sel = 0                 ! the bound container's fill-mode selector
   integer :: seen_nwant = 0               ! the request the container observed (length)
   integer :: seen_want1 = 0               ! the request the container observed (first entry)
   integer :: nexec, nhit                  ! (f) executable-position hit counters
   integer :: v, ios, jj                   ! (f) scan unit / status / scratch
   character(len=512) :: path              ! (f) scan path buffer
   character(len=512) :: exec_path         ! (f) the file that carried the last hit
   character(len=512) :: cmd, ln           ! (f) scan command / line buffer
   logical :: bitexact                     ! (g) the twin's bitwise attempt
   real(8) :: cal, ea, eb, tk              ! the measured conversion + energy channels
   real(8) :: lam, un, u1, u2              ! (g) the closed-form eigenpair
   real(8) :: fs(6), es, et                ! (g) the twin comparison carriers
   real(8), allocatable :: gref(:)         ! (g) the manual-diagonalization reference fold
   ! harmonic test-potential parameters (per-dof, 6 dofs / 2 atoms of mass 1 and 4 amu)
   real(8), parameter :: keV(6)   = [ 0.95d0, 1.15d0, 0.85d0, 1.25d0, 0.75d0, 1.10d0 ] ! stiffness [eV/A^2]
   real(8), parameter :: qe_ctr(6) = [ 0.12d0,-0.18d0, 0.28d0,-0.14d0, 0.22d0,-0.06d0 ] ! centers [A]
   ! criterion-(e) geometries / momenta
   real(8), parameter :: ueA(6) = [ 0.30d0,-0.21d0, 0.18d0,-0.12d0, 0.24d0,-0.09d0 ] ! calibration offset [A]
   real(8), parameter :: ueB(6) = [-0.16d0, 0.26d0,-0.09d0, 0.31d0,-0.20d0, 0.14d0 ] ! second geometry [A]
   real(8), parameter :: pB(6)  = [ 1.05d0,-0.65d0, 0.48d0,-0.85d0, 0.55d0, 1.20d0 ] ! second-geometry momenta
   ! criterion-(g) diab branch: the linear two-surface H_ij (u = q - qe_ctr)
   real(8), parameter :: hb1 = -0.30d0, hb2 = 0.22d0, hb12 = 0.10d0 ! constant H offsets [eV]
   real(8), parameter :: a1(6) = [ 0.80d0,-0.50d0, 0.60d0,-0.40d0, 0.70d0,-0.30d0 ] ! dH11/du [eV/A]
   real(8), parameter :: a2(6) = [-0.40d0, 0.70d0,-0.30d0, 0.80d0,-0.50d0, 0.60d0 ] ! dH22/du [eV/A]
   real(8), parameter :: cc(6) = [ 0.25d0,-0.35d0, 0.15d0,-0.45d0, 0.30d0,-0.20d0 ] ! dH12/du [eV/A]
   real(8), parameter :: g_uni(6) = [ 0.05d0,-0.08d0, 0.06d0,-0.04d0, 0.07d0,-0.03d0 ] ! additive kernel [eV/A]
   real(8), parameter :: v_uni = -0.35d0   ! additive energy [eV]
   real(8), parameter :: ug(6) = [ 0.26d0,-0.19d0, 0.23d0,-0.11d0, 0.27d0,-0.15d0 ] ! (g) geometry [A]
   real(8), parameter :: pg(6) = [ 0.20d0,-0.13d0, 0.17d0,-0.10d0, 0.15d0, 0.22d0 ] ! (g) momenta

   ! degenerate no-op method row (the sibling checks' posture): registers +
   ! binds the name and seeds the statistics holder - occ_get is live for
   ! criterion (d); the members are never dispatched by this batch
   call method_reg('adiabatic', '', noop_prop, noop_match)
   call method_bind('adiabatic')

   ! ---------- (a) handle identity + write-through observability ----------
   call ws_reset()
   p1 => ws_handle()
   p2 => ws_handle()
   allocate (p1%g_add(3))
   p1%g_add = [1.0d0, 2.0d0, 3.0d0]
   p1%v_add = 3.25d0
   call ws_get(w)
   if (.not. (associated(p1) .and. associated(p2) .and. associated(p1, p2)) .or. &
       p2%v_add /= 3.25d0 .or. any(p2%g_add /= [1.0d0, 2.0d0, 3.0d0]) .or. &
       w%v_add /= 3.25d0 .or. any(w%g_add /= [1.0d0, 2.0d0, 3.0d0])) then
      print '(a)', 'FAIL (a) handle identity: the two handles do not alias one '// &
                   'target, or a handle write is invisible to the other handle / ws_get'
      nfail = nfail + 1
   else
      print '(a)', '(a) handle identity: PASS (two ws_handle() calls alias one '// &
                   'target; the write through one is seen by the other and by ws_get)'
   end if

   ! ---------- (b) copy-out semantics: the snapshot is not the holder ----------
   w%v_add = -999.0d0                      ! scribble on the snapshot taken above
   w%g_add = -999.0d0                      ! (a scalar and an array component)
   call ws_get(w2)                         ! a fresh snapshot of the holder
   if (w2%v_add /= 3.25d0 .or. any(w2%g_add /= [1.0d0, 2.0d0, 3.0d0])) then
      print '(a)', 'FAIL (b) copy-out: writing the snapshot perturbed the holder'
      nfail = nfail + 1
   else
      print '(a)', '(b) copy-out: PASS (scribbling the snapshot leaves the holder '// &
                   'intact; a fresh ws_get reproduces the held content)'
   end if

   ! ---------- (c) reset semantics: full clear to the defined empty state ----------
   call ws_reset()                         ! clean slate so every allocation below is fresh
   p1 => ws_handle()
   p1%repr = r_diab
   p1%v_add = -777.0d0
   allocate (p1%i_want(2))
   p1%i_want = [1, 2]
   allocate (p1%g_add(3))
   p1%g_add = 1.0d0
   allocate (p1%h_ij(2, 2))
   p1%h_ij = 1.0d0
   allocate (p1%dh_ij(3, 2, 2))
   p1%dh_ij = 1.0d0
   allocate (p1%v_surf(2))
   p1%v_surf = 1.0d0
   allocate (p1%g_surf(3, 2))
   p1%g_surf = 1.0d0
   allocate (p1%d_ij(3, 2))
   p1%d_ij = 1.0d0
   call ws_reset()
   call ws_get(w)
   if (w%repr /= r_empty .or. w%v_add /= 0.0d0 .or. allocated(w%i_want) .or. &
       allocated(w%g_add) .or. allocated(w%h_ij) .or. allocated(w%dh_ij) .or. &
       allocated(w%v_surf) .or. allocated(w%g_surf) .or. allocated(w%d_ij)) then
      print '(a)', 'FAIL (c) reset: the holder did not return to the defined empty '// &
                   'state (repr, scalar, or a per-mode array survived)'
      nfail = nfail + 1
   else
      print '(a)', '(c) reset: PASS (repr = r_empty, all seven per-mode arrays '// &
                   'gone, scalar zeroed - the constructs die at step end)'
   end if

   ! ---------- (d) snapshot-style vs handle-style coexistence ----------
   call ws_reset()
   call occ_get(occ_req)
   ! the snapshot fill pattern of the scalar force branch, verbatim: request
   ! through occ_get, snapshot through ws_get, compose into the snapshot
   call ws_get(w)
   if (allocated(w%i_want)) deallocate (w%i_want)
   allocate (w%i_want(size(occ_req)))
   w%i_want = occ_req
   allocate (w%g_add(3))
   w%g_add = 0.5d0
   w%v_add = -1.5d0
   w%repr = r_scalar
   ! zero interference: the snapshot fill must leave the holder in its empty
   ! state (the pre-freeze consumer behavior, unchanged by the rewrite)
   p1 => ws_handle()
   if (size(occ_req) /= 1 .or. occ_req(1) /= 1 .or. &
       p1%repr /= r_empty .or. p1%v_add /= 0.0d0 .or. allocated(p1%i_want) .or. &
       allocated(p1%g_add)) then
      print '(a)', 'FAIL (d) coexistence: the occ channel broke, or the snapshot '// &
                   'fill trespassed on the holder'
      nfail = nfail + 1
   else
      print '(a)', '(d) coexistence part 1: PASS (occ request seeds [1]; the '// &
                   'snapshot fill leaves the holder empty - zero interference)'
   end if
   ! the handle-style compose (the fill-site posture onto the frozen holder):
   ! writes the holder itself, then reads back through the snapshot channel
   p1%repr = r_scalar
   allocate (p1%i_want(size(occ_req)))
   p1%i_want = occ_req
   allocate (p1%g_add(3))
   p1%g_add = 0.25d0
   p1%v_add = -2.5d0
   call ws_get(w2)
   if (p1%repr /= r_scalar .or. p1%v_add /= -2.5d0 .or. any(p1%g_add /= 0.25d0) .or. &
       w2%repr /= r_scalar .or. w2%v_add /= -2.5d0 .or. any(w2%g_add /= 0.25d0) .or. &
       any(w2%i_want /= 1) .or. &
       w%v_add /= -1.5d0 .or. any(w%g_add /= 0.5d0)) then
      print '(a)', 'FAIL (d) coexistence: the handle read-back broke, or the earlier '// &
                   'snapshot was retroactively perturbed'
      nfail = nfail + 1
   else
      print '(a)', '(d) coexistence part 2: PASS (the handle write reads back '// &
                   'through ws_get; the earlier snapshot keeps its own content)'
   end if
   call ws_reset()                         ! leave the holder in the empty state

   ! ---------- (e) end-to-end fill chain through the handle ----------
   ! assembly check-side: the roster masses, a 2-atom state, and the local branch
   ! container bound through the PUBLIC assembly face (one slot, bound once - the
   ! (g) fill branches ride the same binding through the fill_sel selector)
   allocate (roster%mass(2))
   roster%mass = [ 1.0d0, 4.0d0 ]
   call state_create(sa, 2)
   call container_bind_pes(check_pes)
   fill_sel = 1
   ! calibration: one force_eval at a known single-dof offset (p = 0) measures the
   ! composed conversion constant cal (the eV->internal product the facade applies;
   ! the literal lives only in force_facade - never written here)
   sa%q = qe_ctr + ueA
   sa%p = 0.0d0
   call force_eval(sa, ea)
   cal = -sa%f(1)/(keV(1)*ueA(1))
   if (cal <= 0.0d0) then
      call fail('(e) measured conversion constant non-positive: '//trim(d2s(cal)))
   end if
   ! geometry A: the analytic force component-wise (proportionality - cal cancels
   ! into the same composition on both sides)
   do m = 1, 6
      if (abs(sa%f(m) + (keV(m)*(sa%q(m) - qe_ctr(m)))*cal) > &
          1.0d-13*max(abs((keV(m)*(sa%q(m) - qe_ctr(m)))*cal), 1.0d0)) then
         call fail('(e) analytic-force mismatch at the calibration geometry, dof '//trim(i2s(m)))
      end if
   end do
   ! geometry B: force + the e_tot decomposition via the measured conversion
   sa%q = qe_ctr + ueB
   sa%p = pB
   call force_eval(sa, eb)
   do m = 1, 6
      if (abs(sa%f(m) + (keV(m)*(sa%q(m) - qe_ctr(m)))*cal) > &
          1.0d-13*max(abs((keV(m)*(sa%q(m) - qe_ctr(m)))*cal), 1.0d0)) then
         call fail('(e) analytic-force mismatch, geometry B, dof '//trim(i2s(m)))
      end if
   end do
   tk = 0.0d0
   do m = 1, 2
      tk = tk + dot_product(sa%p(3*m - 2:3*m), sa%p(3*m - 2:3*m))/(2.0d0*roster%mass(m))
   end do
   tk = tk/e_conv
   if (abs(eb - (tk + 0.5d0*dot_product(keV, (sa%q - qe_ctr)*(sa%q - qe_ctr))*cal/e_conv)) > &
       1.0d-12*max(abs(eb), 1.0d0)) then
      call fail('(e) e_tot decomposition mismatch, geometry B: '//trim(d2s(eb)))
   end if
   ! the request handshake: the container must have consumed i_want = [1] (the held
   ! degenerate occupation seeded by method_bind)
   if (seen_nwant /= 1 .or. seen_want1 /= 1) then
      call fail('(e) the container did not observe the occupation request (i_want /= [1])')
   end if
   ! the holder read-back: AFTER the call a fresh alias AND a snapshot carry the
   ! composed constructs of the returned geometry - the fill reached the holder and
   ! survived the call (one-step conveyance); the references recompute the
   ! container's own expressions from sa%q, so equality is exact. Stepwise guards:
   ! size()/element queries run only behind their allocation checks (a broken fill
   ! must FAIL here, not crash). The PASS line is water-level gated: it prints only
   ! when NO (e) assertion above failed (nfail == nfail0, the watermark taken below)
   nfail0 = nfail
   p1 => ws_handle()
   call ws_get(w2)
   bitexact = (p1%repr == r_scalar .and. allocated(p1%g_add) .and. allocated(p1%i_want))
   if (bitexact) bitexact = (size(p1%g_add) == 6 .and. size(p1%i_want) == 1 .and. &
       p1%i_want(1) == 1 .and. all(p1%g_add == keV*(sa%q - qe_ctr)) .and. &
       p1%v_add == 0.5d0*dot_product(keV, (sa%q - qe_ctr)*(sa%q - qe_ctr)))
   if (.not. bitexact) then
      call fail('(e) holder read-back: the fill did not reach the holder / survive the call')
   else if (w2%repr /= r_scalar .or. .not. allocated(w2%g_add) .or. .not. allocated(w2%i_want) .or. &
       any(w2%g_add /= p1%g_add) .or. w2%v_add /= p1%v_add .or. any(w2%i_want /= 1)) then
      call fail('(e) snapshot read-back: ws_get does not reproduce the held constructs')
   else if (nfail == nfail0) then
      print '(a)', '(e) fill chain through the handle: PASS (analytic force + e_tot at two'// &
         ' geometries; request consumed; the composed constructs sit in the holder after'// &
         ' the call - alias and snapshot agree)'
   end if

   ! ---------- (f) conversion single-site ----------
   ! list every src/ file whose text carries the conversion literal, then count the
   ! EXECUTABLE-position occurrences (comment text does not count); the unit contract
   ! demands exactly one - the facade's conversion-parameter declaration
   nfail0 = nfail
   if (windows_host()) then
      cmd = 'findstr /s /m /l /c:"23.0605" *.f90 > build/fcd_grep.txt 2>&1'
   else
      cmd = "grep -r -l -F --include='*.f90' 23.0605 . > build/fcd_grep.txt 2>&1"
   end if
   call execute_command_line(trim(cmd), wait=.true.)
   nexec = 0
   exec_path = ''
   open (newunit=v, file='build/fcd_grep.txt', status='old', action='read', iostat=ios)
   if (ios /= 0) then
      call fail('(f) the conversion-literal file listing could not be produced')
   else
      do
         read (v, '(a)', iostat=ios) ln
         if (ios /= 0) exit
         path = adjustl(ln)
         jj = len_trim(path)
         if (jj > 0) then
            if (path(jj:jj) == achar(13)) then       ! strip a carried CR terminator
               path(jj:jj) = ' '
               path = adjustl(path)
            end if
         end if
         if (path(1:2) == './') path = adjustl(path(3:))
         if (len_trim(path) == 0) cycle
         call scan_exec_hits(trim(path), nhit)
         if (nhit > 0) then
            nexec = nexec + nhit
            exec_path = trim(path)
         end if
      end do
      close (v)
      if (nexec /= 1 .or. index(exec_path, 'force_facade.f90') == 0) then
         call fail('(f) conversion single-site violated: '//trim(i2s(nexec))// &
            ' executable-position occurrence(s) in src/, last in "'//trim(exec_path)//'"'// &
            ' (want exactly 1, in force_facade.f90)')
      else
         print '(a)', '(f) conversion single-site: PASS (exactly one executable-position'// &
            ' occurrence of the conversion literal in src/: '//trim(exec_path)//')'
      end if
   end if

   ! ---------- (g) the multi-surface diabatic fold ----------
   ! main: the two-surface linear H_ij with an additive scalar member - force_eval
   ! must return the closed-form eigen-fold of the occupied (lowest, ascending)
   ! eigenvalue plus the additive kernel, and the delivered h_ij must SURVIVE in
   ! the holder (the facade diagonalizes a copy)
   fill_sel = 2
   sa%q = qe_ctr + ug
   sa%p = pg
   call force_eval(sa, ea)
   nfail0 = nfail
   lam = 0.5d0*((hb1 + dot_product(a1, ug) + hb2 + dot_product(a2, ug)) &
                - sqrt((hb1 + dot_product(a1, ug) - hb2 - dot_product(a2, ug))**2 &
                       + 4.0d0*(hb12 + dot_product(cc, ug))**2))
   un = sqrt((hb12 + dot_product(cc, ug))**2 + (lam - (hb1 + dot_product(a1, ug)))**2)
   u1 = (hb12 + dot_product(cc, ug))/un
   u2 = (lam - (hb1 + dot_product(a1, ug)))/un
   allocate (gref(6))
   do m = 1, 6
      gref(m) = (u1*u1)*a1(m) + 2.0d0*(u1*u2)*cc(m) + (u2*u2)*a2(m) + g_uni(m)
   end do
   do m = 1, 6
      if (abs(sa%f(m) + gref(m)*cal) > 1.0d-12*max(abs(gref(m)*cal), 1.0d0)) then
         call fail('(g) eigen-fold force mismatch, dof '//trim(i2s(m))//': '// &
            trim(d2s(sa%f(m)))//' vs closed form '//trim(d2s(-gref(m)*cal)))
      end if
   end do
   tk = 0.0d0
   do m = 1, 2
      tk = tk + dot_product(sa%p(3*m - 2:3*m), sa%p(3*m - 2:3*m))/(2.0d0*roster%mass(m))
   end do
   tk = tk/e_conv
   if (abs(ea - (tk + (lam + v_uni)*cal/e_conv)) > 1.0d-12*max(abs(ea), 1.0d0)) then
      call fail('(g) eigen-fold e_tot mismatch: '//trim(d2s(ea))//' vs closed form '// &
         trim(d2s(tk + (lam + v_uni)*cal/e_conv)))
   end if
   p1 => ws_handle()
   bitexact = (p1%repr == r_diab .and. allocated(p1%h_ij) .and. allocated(p1%dh_ij))
   if (bitexact) bitexact = (size(p1%h_ij, 1) == 2 .and. size(p1%h_ij, 2) == 2 .and. &
       all(p1%h_ij == reshape([hb1 + dot_product(a1, ug), hb12 + dot_product(cc, ug), &
                               hb12 + dot_product(cc, ug), hb2 + dot_product(a2, ug)], &
                              [2, 2])))
   if (.not. bitexact) then
      call fail('(g) holder read-back: the delivered h_ij did not survive in the holder')
   else if (nfail == nfail0) then
      print '(a)', '(g) diab fold: PASS (two-surface eigen-fold + additive kernel = closed'// &
         ' form to 1e-12; e_tot = t_kin + (lambda_occ + v_add); h_ij survives in the holder)'
   end if
   deallocate (gref)
   ! twin: the one-surface diab delivery vs the scalar channel at the SAME state -
   ! bitwise attempt (identical arithmetic through the 1x1 eigenpair), 1e-15 fallback
   fill_sel = 1
   sa%q = qe_ctr + ueB
   sa%p = pB
   call force_eval(sa, es)
   fs = sa%f
   fill_sel = 3
   sa%q = qe_ctr + ueB
   sa%p = pB
   call force_eval(sa, et)
   bitexact = (all(sa%f == fs) .and. et == es)
   if (bitexact) then
      print '(a)', '(g) degenerate equivalence: PASS BITWISE (one-surface diab twin ='// &
         ' scalar channel: force and e_tot identical bit for bit)'
   else if (all(abs(sa%f - fs) <= 1.0d-15*max(abs(fs), 1.0d0)) .and. &
            abs(et - es) <= 1.0d-15*max(abs(es), 1.0d0)) then
      print '(a)', '(g) degenerate equivalence: PASS at the 1e-15 fallback (bitwise'// &
         ' attempt failed - arithmetic-order deviation, within tolerance)'
   else
      call fail('(g) degenerate equivalence: one-surface diab twin vs scalar channel'// &
         ' deviates beyond 1e-15 (force max dev '//trim(d2s(maxval(abs(sa%f - fs))))// &
         ', e_tot dev '//trim(d2s(abs(et - es)))//')')
   end if
   call ws_reset()                         ! leave the holder in the empty state
   call state_destroy(sa)

   print '(a,i0,a)', '=== check_facade: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! the local branch test container (bound once; the check-side fill_sel picks the
   ! delivery mode): 1 = the scalar harmonic channel, 2 = the two-surface linear
   ! diabatic channel with an additive scalar member, 3 = the one-surface diabatic
   ! twin of branch 1 (same V and gradient as the scalar branch, delivered as 1x1
   ! constructs; the additive channel stays absent). The gradient and energy in eV
   ! (the member unit contract). Branches 1 and 3 recompute u = q - qe_ctr and the
   ! V expression with IDENTICAL grouping, so the twin comparison is arithmetic-
   ! order clean
   subroutine check_pes(q, mass, ws)
      real(8), intent(in)  :: q(:)     ! coordinates [A]
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused - pure PES branches)
      type(ws_t), intent(inout) :: ws  ! the holder (received through the facade's alias)
      real(8) :: u(6)
      u = q - qe_ctr
      seen_nwant = 0                    ! the request handshake log: what the container
      seen_want1 = 0                    ! observed in the pre-written request field
      if (allocated(ws%i_want)) then
         seen_nwant = size(ws%i_want)
         seen_want1 = ws%i_want(1)
      end if
      select case (fill_sel)
      case (1)                          ! scalar harmonic: V = 1/2 sum k u^2
         allocate (ws%g_add(size(q)))
         ws%g_add = keV*u                                 ! dV/dq [eV/A]
         ws%v_add = 0.5d0*dot_product(keV, u*u)           ! V [eV]
         ws%repr = r_scalar
      case (2)                          ! two-surface linear diabatic H + additive member
         allocate (ws%h_ij(2, 2), ws%dh_ij(size(q), 2, 2), ws%g_add(size(q)))
         ws%h_ij(1, 1) = hb1 + dot_product(a1, u)
         ws%h_ij(2, 2) = hb2 + dot_product(a2, u)
         ws%h_ij(1, 2) = hb12 + dot_product(cc, u)
         ws%h_ij(2, 1) = ws%h_ij(1, 2)
         ws%dh_ij(:, 1, 1) = a1                            ! dH_ij/du [eV/A] (symmetric)
         ws%dh_ij(:, 2, 2) = a2
         ws%dh_ij(:, 1, 2) = cc
         ws%dh_ij(:, 2, 1) = cc
         ws%g_add = g_uni                                  ! the additive member's kernel
         ws%v_add = v_uni                                  ! and energy [eV]
         ws%repr = r_diab
      case (3)                          ! one-surface diabatic twin of branch 1
         allocate (ws%h_ij(1, 1), ws%dh_ij(size(q), 1, 1))
         ws%dh_ij(:, 1, 1) = keV*u                         ! same expression as case (1)
         ws%h_ij(1, 1) = 0.5d0*dot_product(keV, u*u)       ! same expression as case (1)
         ws%repr = r_diab                                  ! (v_add stays at its zeroed default)
      end select
   end subroutine check_pes

   ! degenerate no-op method members (the all-noop package this check registers
   ! so the method-name lookup resolves and the statistics holder is seeded)
   subroutine noop_prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine noop_prop

   subroutine noop_match(sta)
      type(state_t), intent(inout) :: sta
   end subroutine noop_match

   ! scan_exec_hits(path, nhit) - count this file's EXECUTABLE-position occurrences
   ! of the conversion literal: a hit inside comment text (a '!' at or before the
   ! hit position) does not count; a trailing comment after code counts (conservative
   ! - the single-site statute fails loudly rather than silently)
   subroutine scan_exec_hits(path, nhit)
      character(len=*), intent(in) :: path
      integer, intent(out) :: nhit
      integer :: uu, ios2, p, cq
      character(len=512) :: l2
      nhit = 0
      open (newunit=uu, file=path, status='old', action='read', iostat=ios2)
      if (ios2 /= 0) return
      do
         read (uu, '(a)', iostat=ios2) l2
         if (ios2 /= 0) exit
         p = index(l2, '23.0605')
         if (p > 0) then
            cq = index(l2, '!')
            if (cq == 0 .or. cq > p) nhit = nhit + 1
         end if
      end do
      close (uu)
   end subroutine scan_exec_hits

   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

   ! fail(msg) - record one criterion failure
   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      print '(a)', 'FAIL '//trim(msg)
      nfail = nfail + 1
   end subroutine fail

   ! small int/real to string helpers for failure messages
   function i2s(i)
      integer, intent(in) :: i
      character(len=12) :: i2s
      write (i2s, '(i0)') i
   end function i2s
   function d2s(x)
      real(8), intent(in) :: x
      character(len=24) :: d2s
      write (d2s, '(es12.4)') x
   end function d2s

end program check_facade
