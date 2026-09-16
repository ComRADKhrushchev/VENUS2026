!=====================================================================
! check_recorder.f90 - C4.5 acceptance criteria program (criteria defined FIRST in
!                      docs/plans/C4.md; this program executes them)
! Compile (from src/; modules in use-order first, -Jbuild keeps the .mod files out
! of the source tree):
!   gfortran -Jbuild globals/consts.f90 globals/control.f90 globals/config.f90 &
!            utils/specio.f90 facade/sysdef_facade.f90 input.f90 globals/state.f90 &
!            facade/elec_facade.f90 globals/config_roster.f90 &
!            utils/linalg.f90 facade/force_facade.f90 &
!            output/recorder.f90 ../scripts/check_recorder.f90 -o build/check_recorder.exe
! Sample output files are written under build/ and read back line-by-line - the
! source tree stays clean.
! S1 boundary (C4.md): the energy columns t_kin/v_pot/e_tot pass through force_eval,
!   now REAL on a check-bound zero-force container (degenerate method row +
!   container_bind_pes + container_init, assembled check-side) - their VALUES stay
!   unasserted here (zero potential: only presence / order / header, the column
!   MECHANISM). Numeric assertions cover the constructible quantities: frame-header
!   step/trajectory/time and the sec_qp q/p six-tuple rows (hand-built state_t). The
!   force section is asserted for GATING + header/line count only (the E20 posture).
! Placeholder statute -> reality: the occ / n_hop columns are REAL - the extractors
!   read the elec_facade statistics holder through the getters (occ_get /
!   n_hop_get), seeded at the method bind. Under the degenerate adiabatic row
!   the holder carries occ = [1] constant and n_hop = 0, so criterion (p)
!   asserts exactly those values end to end (file digits + the final-row
!   read-back) - the numeric channel for any surface-hopping method later.
! Loud-stop verification approach (same scheme as check_sysdef/check_roster): the
!   loud-stop criteria run as CHILD PROCESSES - the check re-spawns itself with a
!   mode argument (execute_command_line); the child triggers the real loud-stop path
!   and must print a naming message on stderr and exit nonzero.
! Criteria (C4.md C4.5 + D3/N3.1):
!   (a) registration + level 0/1 minimal group: select level 0 -> one frame carries
!       the frame header (traj=1 step=1 t exact) + the column line naming exactly
!       t_kin/v_pot/e_tot in registry order; no section markers; level 1 selects the
!       SAME column set (the old level-1 diagnostic block was excised system
!       knowledge - no generic level-1 columns exist)
!   (b) sec_qp gate + content: level 2 -> '# sec_qp begin/end' markers with exactly
!       natoms six-tuple lines between them; per-atom index/symbol/3q/3p numerically
!       exact against the hand-built state (constructible quantities); labels from
!       roster%symb (excision fix (2)); sec_frc absent at level 2 (positional proof)
!   (c) sec_frc gate: level 5 -> '# sec_frc begin/end' markers with exactly natoms
!       lines between them (gating + header only - force VALUES unasserted, E20)
!   (d) list entry + union semantics: level 1 then list 'n_hop,T_KIN,occ' ->
!       selection = the preset three (registry order) with n_hop/occ APPENDED in
!       list order (case-insensitive lookup, duplicate t_kin not re-listed)
!   (e) reserved section names: level 1 + list 'all_qp' opens the qp section; a
!       further list 'all_frc' opens the force section - gates union across the two
!       selection entries while the column set stays the preset three
!   (f) frame-header numeric channel: successive rec_frame calls advance step
!       1,2,3...; t equals sta%t per frame (0.7 / 1.4); after rec_final the counters
!       advance - the next frame reads traj=2 step=1; the selection/gates persist
!       across the rec_final boundary (only the counters reset)
!   (g) final-state line: rec_final at traj=1 writes the name header
!       ('# final: traj t <selected column names>') + one row (traj=1, t exact,
!       values parse); a second rec_final (traj=2) writes NO second header
!   (h) duplicate-name registration loud-stop (child)
!   (i) registry-capacity loud-stop (child): overflowing the 64-column registry
!       loud-stops naming the overflow
!   (j) unknown-column-name loud-stop in the list (child)
!   (k) empty-list-item loud-stop (child - no silent swallow)
!   (l) illegal level loud-stop (child: rec_select(3) - the domain is {0,1,2,5},
!       the C4.4 GWRITE_LEVEL set)
!   (m) rec_frame before any rec_select loud-stop (child - assembly-order error)
!   (n) state-larger-than-roster loud-stop (child, fix round 1): a 3-atom state against
!       a 2-atom roster hits the named size-mismatch stop - the equality guard fires
!       BEFORE col_tkin's roster%mass walk could read out of bounds
!   (o) final-row read-back (the statistics-phase value channel): rec_fin_val(name)
!       hands back the LAST written final row's named column value - the analytic
!       e_tot fold to machine precision, the case-insensitive name lookup, and the
!       e_tot = t_kin identity on the zero container; two children cover the named
!       loud-stops (a name outside the selection; a read before any rec_final)
!   (p) occupation columns on the real holder (the column-activation value
!       channel): level 1 + the list 'occ,n_hop' -> the frame column line and
!       the final row carry occ = 1.0 and n_hop = 0.0 exactly (the degenerate
!       adiabatic holder values through registration -> selection -> write
!       order), and rec_fin_val reads the same pair back - the registration /
!       selection / write-order consistency of the two columns with the energy
!       family, plus the whole-file trailing check (moved here from (f) when
!       the write sequence grew the third trajectory)
!=====================================================================
program check_recorder
   use recorder,     only: col_val_i, rec_reg_col, rec_select, rec_frame, rec_final, &
                           rec_fin_val
   use consts,       only: e_conv
   use state,        only: state_t, state_create, state_destroy
   use config_roster, only: roster
   use elec_facade,  only: method_reg
   use force_facade, only: container_bind_pes, container_init
   implicit none
   integer :: nfail = 0, nfail0
   character(len=256) :: mode
   character(len=512), allocatable :: lines(:)
   character(len=512) :: line
   character(len=16) :: names(8)
   real(8) :: vals(8), r6(6), tt, tol
   integer :: nlines, lp, nv, i, j, u, ios, traj, step
   type(state_t) :: sta
   ! hand-built constructible state (2 atoms): q/p/f exact short decimals, roster
   ! masses/symbols set directly (the roster instance is public and input-independent
   ! here - no sysdef chain needed for the recording mechanism)
   real(8), parameter :: q2(6) = [ 0.5d0,-1.0d0, 2.0d0, 0.25d0, 1.25d0,-0.75d0 ]
   real(8), parameter :: p2(6) = [ 1.0d0, 2.0d0, 3.0d0,-0.5d0, 0.5d0, 4.0d0 ]
   real(8) :: tk_ana                  ! the analytic kinetic fold [kcal/mol] (criterion o)
   tol = 1.0d-13

   ! ---------- child dispatch: one bad case per spawn, must loud-stop ----------
   if (command_argument_count() > 0) then
      call get_command_argument(1, mode)
      select case (trim(mode))
      case ('c-dup')                      ! framework t_kin preregistered -> duplicate
         call rec_reg_col('t_kin', col_probe)
      case ('c-full')                     ! 5 framework + 61 probes > 64 capacity
         do i = 1, 61
            write (mode, '(a,i0)') 'probe', i
            call rec_reg_col(trim(mode), col_probe)
         end do
      case ('c-unk')                      ! unregistered name in the list
         call rec_select('no_such_col')
      case ('c-empty')                    ! empty item in the list (no silent swallow)
         call rec_select('t_kin,,v_pot')
      case ('c-lvl')                      ! level outside the {0,1,2,5} domain
         call rec_select(3)
      case ('c-noselect')                 ! frame before any selection
         call state_create(sta, 1)
         call rec_frame(sta, 11)
      case ('c-big')                      ! state LARGER than the roster (fix round 1:
         allocate (roster%mass(2), roster%symb(2)) ! the old one-way guard let this read
         roster%mass = [ 1.5d0, 3.0d0 ]          ! roster%mass out of bounds in col_tkin)
         roster%symb = [ character(4) :: 'H', 'He' ]
         call rec_select(0)
         call state_create(sta, 3)
         call rec_frame(sta, 11)
      case ('c-nofin')                    ! statistics read before ANY rec_final ran
         print *, rec_fin_val('t_kin')    ! no row stashed yet - must loud-stop
      case ('c-badfin')                   ! a name outside the selection (the row
         allocate (roster%mass(2), roster%symb(2)) ! carries the selection only)
         roster%mass = [ 1.5d0, 3.0d0 ]
         roster%symb = [ character(4) :: 'H', 'He' ]
         call method_reg('adiabatic', noop_prop, noop_match)
         call container_bind_pes(zero_pes)
         call container_init()            ! the energy columns ride force_eval
         call rec_select(0)
         call state_create(sta, 2)
         open (newunit=u, file='build/rec_finchild.txt', status='replace', action='write')
         call rec_final(sta, u)           ! a good row first (the stash goes live)
         print *, rec_fin_val('zzz')      ! unselected name - must loud-stop
      case default
         print '(a)', 'child: unknown mode '//trim(mode)
         stop 2
      end select
      print '(a)', 'child: NO LOUD-STOP (returned - bug)'
      stop 0
   end if

   ! ---------- hand-built roster + physical state (all under the parent only) ----------
   allocate (roster%mass(2), roster%symb(2))
   roster%mass = [ 1.5d0, 3.0d0 ]
   roster%symb = [ character(4) :: 'H', 'He' ]
   call state_create(sta, 2)
   sta%q = q2
   sta%p = p2
   sta%f = [ 0.125d0,-0.25d0, 0.375d0, 1.5d0,-2.5d0, 3.5d0 ]

   ! the force chain the energy columns ride (real since the scalar force path
   ! landed): degenerate method row + check-local ZERO-force container + the
   ! single assembly entry - v_pot lands 0 on it, t_kin the analytic kinetic fold.
   ! Assembled CHECK-SIDE; the tree's default binding posture stays untouched
   call method_reg('adiabatic', noop_prop, noop_match)
   call container_bind_pes(zero_pes)
   call container_init()

   open (newunit=u, file='build/rec_frames.txt', status='replace', action='write')

   ! frame 1: level 0 (minimal preset group)
   sta%t = 0.7d0
   call rec_select(0)
   call rec_frame(sta, u)
   ! frame 2: level 1 (same generic column set - no level-1 columns exist)
   call rec_select(1)
   call rec_frame(sta, u)
   ! frame 3: level 2 (+ sec_qp)
   call rec_select(2)
   call rec_frame(sta, u)
   ! frame 4: level 5 (+ sec_qp + sec_frc)
   call rec_select(5)
   call rec_frame(sta, u)
   ! frame 5: level 1 + list (union; case-insensitive; duplicate not re-listed)
   call rec_select(1)
   call rec_select('n_hop,T_KIN,occ')
   call rec_frame(sta, u)
   ! frame 6: level 1 + reserved name all_qp (gate opens, columns stay the preset 3)
   call rec_select(1)
   call rec_select('all_qp')
   call rec_frame(sta, u)
   ! frame 7: further reserved name all_frc (both gates open)
   call rec_select('all_frc')
   call rec_frame(sta, u)
   ! final row of trajectory 1 (writes the name header once), then a frame of
   ! trajectory 2 (counters advanced), then the trajectory-2 final row (no header)
   call rec_final(sta, u)
   sta%t = 1.4d0
   call rec_frame(sta, u)
   call rec_final(sta, u)
   ! trajectory 3: the occupation columns live - level 1 + the occ/n_hop list
   ! (the same union the driver assembles from the OBSERVABLES key); one frame
   ! and one final row carry the degenerate adiabatic holder values
   call rec_select(1)
   call rec_select('occ,n_hop')
   sta%t = 2.1d0
   call rec_frame(sta, u)
   call rec_final(sta, u)
   close (u)

   ! ---------- read the whole file back ----------
   allocate (lines(400))
   nlines = 0
   open (newunit=u, file='build/rec_frames.txt', status='old', action='read')
   do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      nlines = nlines + 1
      lines(nlines) = line
   end do
   close (u)
   lp = 1

   ! ---------- (a) registration + level 0/1 minimal group ----------
   ! (absence of sections is proven POSITIONALLY: the line after each column line is
   !  the next frame header - any section marker in between fails expect_frame)
   nfail0 = nfail
   call expect_frame(1, 1, 0.7d0, tol)
   call expect_cols([character(16) :: 't_kin','v_pot ','e_tot '])
   call expect_frame(1, 2, 0.7d0, tol)
   call expect_cols([character(16) :: 't_kin','v_pot ','e_tot '])
   if (nfail == nfail0) print '(a)', &
      '(a) registration + level 0/1: PASS (energy three in registry order, header exact, no sections)'

   ! ---------- (b) sec_qp gate + content ----------
   nfail0 = nfail
   call expect_frame(1, 3, 0.7d0, tol)
   call expect_cols([character(16) :: 't_kin','v_pot ','e_tot '])
   call expect_pref('# sec_qp begin')
   call expect_qp_atom(1, 'H', q2(1:3), p2(1:3), tol)
   call expect_qp_atom(2, 'He', q2(4:6), p2(4:6), tol)
   call expect_pref('# sec_qp end')
   if (nfail == nfail0) print '(a)', &
      '(b) sec_qp: PASS (level-2 gate open, 2 six-tuple rows numerically exact, symb labels)'

   ! ---------- (c) sec_frc gate ----------
   nfail0 = nfail
   call expect_frame(1, 4, 0.7d0, tol)
   call expect_cols([character(16) :: 't_kin','v_pot ','e_tot '])
   call expect_pref('# sec_qp begin')
   lp = lp + 2                              ! two qp rows (already proven numerically)
   call expect_pref('# sec_qp end')
   call expect_pref('# sec_frc begin')
   lp = lp + 2                              ! two force rows - VALUES unasserted (E20)
   call expect_pref('# sec_frc end')
   if (nfail == nfail0) print '(a)', &
      '(c) sec_frc: PASS (level-5 gate open, marker pair + natoms lines, values unasserted per S1)'

   ! ---------- (d) list entry + union semantics ----------
   nfail0 = nfail
   call expect_frame(1, 5, 0.7d0, tol)
   call expect_cols([character(16) :: 't_kin','v_pot ','e_tot ','n_hop ','occ   '])
   if (nfail == nfail0) print '(a)', &
      '(d) list entry + union: PASS (preset 3 + n_hop/occ appended in list order, no duplicate)'

   ! ---------- (e) reserved section names ----------
   nfail0 = nfail
   call expect_frame(1, 6, 0.7d0, tol)
   call expect_cols([character(16) :: 't_kin','v_pot ','e_tot '])
   call expect_pref('# sec_qp begin')
   lp = lp + 2
   call expect_pref('# sec_qp end')
   ! sec_frc must stay closed here: the next line is frame 7's header (positional)
   call expect_frame(1, 7, 0.7d0, tol)
   call expect_cols([character(16) :: 't_kin','v_pot ','e_tot '])
   call expect_pref('# sec_qp begin')
   lp = lp + 2
   call expect_pref('# sec_qp end')
   call expect_pref('# sec_frc begin')
   lp = lp + 2
   call expect_pref('# sec_frc end')
   if (nfail == nfail0) print '(a)', &
      '(e) reserved names: PASS (all_qp/all_frc open the gates, column set stays the preset 3)'

   ! ---------- (g) final-state line ----------
   nfail0 = nfail
   call expect_final(1, 0.7d0, 3, tol)
   if (nfail == nfail0) print '(a)', &
      '(g) final line: PASS (name header once, row traj=1 t exact, values parse)'

   ! ---------- (f) frame-header numeric channel across trajectories ----------
   ! the gates/selection of frame 7 carry into trajectory 2 (module state persists
   ! across rec_final - documented semantics: only the counters reset)
   nfail0 = nfail
   call expect_frame(2, 1, 1.4d0, tol)
   call expect_cols([character(16) :: 't_kin','v_pot ','e_tot '])
   call expect_pref('# sec_qp begin')
   lp = lp + 2
   call expect_pref('# sec_qp end')
   call expect_pref('# sec_frc begin')
   lp = lp + 2
   call expect_pref('# sec_frc end')
   call expect_final_row_only(2, 1.4d0, 3, tol)   ! selection carries into traj 2; NO header
   if (nfail == nfail0) print '(a)', &
      '(f) header channel: PASS (step 1..7, t exact per frame, traj/step advance past rec_final)'

   ! ---------- (o) final-row read-back (the statistics-phase value channel) ----------
   ! the module stash holds the LAST final row (trajectory 3's - the (p) write
   ! above overwrote it); the state never changed q/p, so its values equal
   ! trajectory 2's parsed row digit for digit, and the analytic fold still
   ! names them - the stash is checked against all three
   nfail0 = nfail
   tk_ana = (dot_product(p2(1:3), p2(1:3))/(2.0d0*roster%mass(1)) + &
             dot_product(p2(4:6), p2(4:6))/(2.0d0*roster%mass(2)))/e_conv
   if (abs(rec_fin_val('e_tot') - tk_ana) > 1.0d-10) then
      call fail('(o) rec_fin_val(e_tot) off the analytic fold: '//d2s(rec_fin_val('e_tot'))// &
                ' vs '//d2s(tk_ana)//' kcal/mol')
   end if
   if (abs(rec_fin_val('T_KIN') - tk_ana) > 1.0d-10) then
      call fail('(o) case-insensitive lookup failed: rec_fin_val(T_KIN) = '// &
                d2s(rec_fin_val('T_KIN'))//' vs t_kin analytic '//d2s(tk_ana))
   end if
   if (abs(rec_fin_val('v_pot')) > 1.0d-10) then
      call fail('(o) v_pot must vanish on the zero container: '//d2s(rec_fin_val('v_pot')))
   end if
   i = 0
   do j = 1, nv
      if (trim(names(j)) == 'e_tot') i = j
   end do
   if (i == 0) then
      call fail('(o) e_tot missing from the parsed final-row names')
   else if (abs(rec_fin_val('e_tot') - vals(i)) > 1.0d-9) then
      call fail('(o) the stash disagrees with the written row: '//d2s(rec_fin_val('e_tot'))// &
                ' vs file '//d2s(vals(i)))
   end if
   if (nfail == nfail0) print '(a)', &
      '(o) final-row read-back: PASS (e_tot = t_kin = analytic fold to 1e-10, case-insensitive'// &
      ' lookup, stash = written row)'

   ! ---------- (p) occupation columns on the real holder (value channel) ----------
   ! trajectory 3's frame: the level-1 base call REPLACED the selection (the
   ! section gates closed again), then the occ/n_hop list appended - the column
   ! line carries the five tokens in family order with the two degenerate
   ! holder values EXACT; the final row carries the same pair positionally
   ! (rows are positional per the contract - the once-written header still
   ! names trajectory 1's three columns); the stash reads them back by name;
   ! the whole-file trailing check lives here now that the write sequence
   ! grew past (f)'s horizon
   nfail0 = nfail
   call expect_frame(3, 1, 2.1d0, tol)
   call expect_occ_line()
   call expect_final_row_only(3, 2.1d0, 5, tol)
   if (vals(4) /= 1.0d0 .or. vals(5) /= 0.0d0) then
      call fail('(p) the final row carries occ='//d2s(vals(4))//' n_hop='//d2s(vals(5))// &
                ' (want exactly 1.0 / 0.0)')
   end if
   if (rec_fin_val('occ') /= 1.0d0) then
      call fail('(p) rec_fin_val(occ) = '//d2s(rec_fin_val('occ'))//', want exactly 1.0')
   end if
   if (rec_fin_val('n_hop') /= 0.0d0) then
      call fail('(p) rec_fin_val(n_hop) = '//d2s(rec_fin_val('n_hop'))//', want exactly 0.0')
   end if
   if (lp <= nlines) then
      call fail('(p) unexpected trailing lines in the frame file (cursor at line '// &
                trim(i2s(lp))//' of '//trim(i2s(nlines))//')')
   end if
   if (nfail == nfail0) print '(a)', &
      '(p) occupation columns: PASS (column line + final row carry occ = 1.0 and'// &
      ' n_hop = 0.0 exactly in family order; rec_fin_val reads the pair back;'// &
      ' no trailing lines)'

   ! ---------- (h)+ loud-stop criteria (child processes) ----------
   call expect_loudstop('c-dup',     'duplicate',            'rec_reg_col', &
                        '(h) duplicate-name registration loud-stop')
   call expect_loudstop('c-full',    'full',                 'rec_reg_col', &
                        '(i) registry-capacity loud-stop')
   call expect_loudstop('c-unk',     'no_such_col',          'rec_sel_lst', &
                        '(j) unknown-column-name loud-stop')
   call expect_loudstop('c-empty',   'empty item',           'rec_sel_lst', &
                        '(k) empty-list-item loud-stop')
   call expect_loudstop('c-lvl',     'level',                'rec_sel_lvl', &
                        '(l) illegal-level loud-stop (domain 0/1/2/5)')
   call expect_loudstop('c-noselect','rec_select',           'rec_frame', &
                        '(m) frame-before-selection loud-stop')
   call expect_loudstop('c-big',     'size mismatch',        'rec_frame', &
                        '(n) state-larger-than-roster loud-stop (fix round 1)')
   call expect_loudstop('c-nofin',   'no final-state row',   'rec_fin_val', &
                        '(o1) read-before-any-rec_final loud-stop')
   call expect_loudstop('c-badfin',  'zzz',                  'rec_fin_val', &
                        '(o2) name-outside-the-selection loud-stop')

   call state_destroy(sta)
   print '(a,i0,a)', '=== check_recorder: ', nfail, ' failure(s) ==='
   if (nfail > 0) stop 1
   stop 0

contains

   ! a probe extractor conforming to col_val_i (the registration-mechanism test dummy)
   real(8) function col_probe(sta) result(v)
      type(state_t), intent(in) :: sta
      v = 0.0d0
   end function col_probe

   ! the check-local zero container: empty kernel and energy (eV member units)
   ! - v_pot/e_tot ride the zero potential, the mechanism-level posture
   subroutine zero_pes(q, mass, g, v)
      real(8), intent(in)  :: q(:)     ! coordinates [A] (unused - zero force)
      real(8), intent(in)  :: mass(:)  ! per-atom masses [amu] (unused)
      real(8), intent(out) :: g(:)     ! the composed kernel dV/dq [eV/A] (zero)
      real(8), intent(out) :: v        ! V [eV] (zero)
      g = 0.0d0
      v = 0.0d0
   end subroutine zero_pes

   ! degenerate no-op method members (the row container_init resolves)
   subroutine noop_prop(sta, dt)
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine noop_prop

   subroutine noop_match(sta)
      type(state_t), intent(inout) :: sta
   end subroutine noop_match

   ! fail(msg) - record one criterion failure
   subroutine fail(msg)
      character(len=*), intent(in) :: msg
      print '(a)', 'FAIL '//trim(msg)
      nfail = nfail + 1
   end subroutine fail

   ! fetch the current line (advancing the cursor) - a missing line is a failure
   logical function next_line(l)
      character(len=*), intent(out) :: l
      next_line = .false.
      l = ''
      if (lp > nlines) then
         call fail('unexpected end of the frame file (cursor at line '//trim(i2s(lp))//')')
         return
      end if
      l = lines(lp)
      lp = lp + 1
      next_line = .true.
   end function next_line

   ! assert the current line is a frame header with the exact traj/step/t values
   subroutine expect_frame(wtraj, wstep, wt, tol)
      integer, intent(in) :: wtraj, wstep
      real(8), intent(in) :: wt, tol
      character(len=512) :: l
      integer :: p1, p2, p3, ios
      if (.not. next_line(l)) return
      p1 = index(l, 'traj=')
      p2 = index(l, ' step=')
      p3 = index(l, ' t=')
      if (index(l, '# frame') == 0 .or. p1 == 0 .or. p2 == 0 .or. p3 == 0) then
         call fail('not a frame header line: '//trim(l))
         return
      end if
      read (l(p1+5:p2-1), *, iostat=ios) traj
      if (ios /= 0) call fail('unparseable traj field: '//trim(l))
      read (l(p2+6:p3-1), *, iostat=ios) step
      if (ios /= 0) call fail('unparseable step field: '//trim(l))
      read (l(p3+3:), *, iostat=ios) tt
      if (ios /= 0) call fail('unparseable t field: '//trim(l))
      if (traj /= wtraj .or. step /= wstep) then
         call fail('frame header traj/step mismatch: got '//trim(i2s(traj))//'/'//trim(i2s(step))// &
                   ' want '//trim(i2s(wtraj))//'/'//trim(i2s(wstep)))
      end if
      if (abs(tt - wt) > tol) then
         call fail('frame header t mismatch: '//trim(d2s(tt))//' vs '//trim(d2s(wt)))
      end if
   end subroutine expect_frame

   ! assert the current line is the column line with exactly the given names, in order
   subroutine expect_cols(want)
      character(len=*), intent(in) :: want(:)
      character(len=512) :: l, rest, tok
      character(len=16) :: got(8)
      integer :: n, i, ie, ios
      real(8) :: v
      if (.not. next_line(l)) return
      rest = adjustl(l)
      n = 0
      do while (len_trim(rest) > 0 .and. n < 8)
         i = index(rest, ' ')
         if (i == 0) i = len_trim(rest) + 1
         tok = rest(:i-1)
         rest = adjustl(rest(i:))
         ie = index(tok, '=')
         if (ie == 0) then
            call fail('column token without "'//'=": '//trim(tok))
            return
         end if
         n = n + 1
         got(n) = tok(:ie-1)
         read (tok(ie+1:), *, iostat=ios) v
         if (ios /= 0) then
            call fail('column value unparseable: '//trim(tok))
            return
         end if
         if (v - v /= 0.0d0) call fail('column value not finite: '//trim(tok))  ! NaN probe
      end do
      if (size(want) /= n) then
         call fail('column count mismatch: got '//trim(i2s(n))//' want '//trim(i2s(size(want))))
         return
      end if
      do i = 1, n
         if (got(i) /= want(i)) then
            call fail('column order mismatch at '//trim(i2s(i))//': got "'//trim(got(i))// &
                      '" want "'//trim(want(i))//'"')
            return
         end if
      end do
   end subroutine expect_cols

   ! assert the current line is a sec_qp atom row with exact index/symbol/q/p
   subroutine expect_qp_atom(wi, wsym, wq, wp, tol)
      integer, intent(in) :: wi
      character(len=*), intent(in) :: wsym
      real(8), intent(in) :: wq(3), wp(3), tol
      character(len=512) :: l
      character(len=4) :: sym
      integer :: iat, ios, k
      if (.not. next_line(l)) return
      read (l, *, iostat=ios) iat, sym, r6
      if (ios /= 0) then
         call fail('unparseable qp row: '//trim(l))
         return
      end if
      if (iat /= wi .or. trim(sym) /= trim(wsym)) then
         call fail('qp row index/symbol mismatch: '//trim(i2s(iat))//' '//trim(sym)// &
                   ' want '//trim(i2s(wi))//' '//trim(wsym))
         return
      end if
      do k = 1, 3
         if (abs(r6(k) - wq(k)) > tol .or. abs(r6(k+3) - wp(k)) > tol) then
            call fail('qp row value mismatch at component '//trim(i2s(k))//': '//trim(l))
            return
         end if
      end do
   end subroutine expect_qp_atom

   ! assert the current line starts with the given prefix
   subroutine expect_pref(pref)
      character(len=*), intent(in) :: pref
      character(len=512) :: l
      if (.not. next_line(l)) return
      if (index(l, trim(pref)) /= 1) call fail('expected "'//trim(pref)//'" got: '//trim(l))
   end subroutine expect_pref

   ! assert the final-state block: name header (once) + row with exact traj/t and
   ! nval parseable values; consumes both lines
   subroutine expect_final(wtraj, wt, nval, tol)
      integer, intent(in) :: wtraj, nval
      real(8), intent(in) :: wt, tol
      character(len=512) :: l
      integer :: ios
      if (.not. next_line(l)) return
      if (index(l, '# final: traj t') /= 1) then
         call fail('missing final-state name header: '//trim(l))
         return
      end if
      call parse_final_header(l, names, nv)
      if (nv /= nval) call fail('final header column count mismatch: '//trim(i2s(nv))// &
                                ' want '//trim(i2s(nval)))
      if (.not. next_line(l)) return
      read (l, *, iostat=ios) traj, tt, vals(1:max(nval,1))
      if (ios /= 0) then
         call fail('unparseable final row: '//trim(l))
         return
      end if
      if (traj /= wtraj .or. abs(tt - wt) > tol) then
         call fail('final row traj/t mismatch: '//trim(i2s(traj))//' '//trim(d2s(tt)))
      end if
   end subroutine expect_final

   ! assert a final-state ROW only (no header repeats from trajectory 2 on)
   subroutine expect_final_row_only(wtraj, wt, nval, tol)
      integer, intent(in) :: wtraj, nval
      real(8), intent(in) :: wt, tol
      character(len=512) :: l
      integer :: ios
      if (.not. next_line(l)) return
      if (index(l, '#') == 1) then
         call fail('unexpected final-state header repeat: '//trim(l))
         return
      end if
      read (l, *, iostat=ios) traj, tt, vals(1:max(nval,1))
      if (ios /= 0) then
         call fail('unparseable final row: '//trim(l))
         return
      end if
      if (traj /= wtraj .or. abs(tt - wt) > tol) then
         call fail('final row traj/t mismatch: '//trim(i2s(traj))//' '//trim(d2s(tt)))
      end if
   end subroutine expect_final_row_only

   ! assert the current line is the occupation column line: exactly five
   ! name=value tokens in family order (the energy three, then occ/n_hop in
   ! list order - the registration/selection/write-order consistency), with
   ! occ and n_hop carrying the degenerate holder values EXACTLY (the (p)
   ! value channel; the energy values stay unasserted per the S1 boundary)
   subroutine expect_occ_line()
      character(len=512) :: l, rest, tok
      character(len=16) :: got(8)
      real(8) :: wv(8), v
      integer :: n, i, ie, ios
      if (.not. next_line(l)) return
      rest = adjustl(l)
      n = 0
      do while (len_trim(rest) > 0 .and. n < 8)
         i = index(rest, ' ')
         if (i == 0) i = len_trim(rest) + 1
         tok = rest(:i-1)
         rest = adjustl(rest(i:))
         ie = index(tok, '=')
         if (ie == 0) then
            call fail('(p) column token without "=": '//trim(tok))
            return
         end if
         n = n + 1
         got(n) = tok(:ie-1)
         read (tok(ie+1:), *, iostat=ios) v
         if (ios /= 0) then
            call fail('(p) column value unparseable: '//trim(tok))
            return
         end if
         wv(n) = v
      end do
      if (n /= 5) then
         call fail('(p) expected the 5-token column line (energy three + occ/n_hop),'// &
                   ' got '//trim(i2s(n))//' tokens: '//trim(l))
         return
      end if
      if (trim(got(1)) /= 't_kin' .or. trim(got(2)) /= 'v_pot' .or. &
          trim(got(3)) /= 'e_tot' .or. trim(got(4)) /= 'occ' .or. &
          trim(got(5)) /= 'n_hop') then
         call fail('(p) column order mismatch: got "'//trim(got(1))//' '//trim(got(2))//' '// &
                   trim(got(3))//' '//trim(got(4))//' '//trim(got(5))// &
                   '" want "t_kin v_pot e_tot occ n_hop"')
         return
      end if
      if (wv(4) /= 1.0d0) then
         call fail('(p) the column line occ = '//d2s(wv(4))//', want exactly 1.0')
      end if
      if (wv(5) /= 0.0d0) then
         call fail('(p) the column line n_hop = '//d2s(wv(5))//', want exactly 0.0')
      end if
   end subroutine expect_occ_line

   ! parse '# final: traj t <names...>' into the trailing names
   subroutine parse_final_header(l, nm, n)
      character(len=*), intent(in) :: l
      character(len=16), intent(out) :: nm(8)
      integer, intent(out) :: n
      character(len=512) :: rest, tok
      integer :: i, pfx
      n = 0
      nm = ''
      pfx = index(l, '# final: traj t')
      if (pfx == 0) return
      rest = adjustl(l(pfx + len('# final: traj t'):))
      do while (len_trim(rest) > 0 .and. n < 8)
         i = index(rest, ' ')
         if (i == 0) i = len_trim(rest) + 1
         n = n + 1
         nm(n) = rest(:i-1)
         rest = adjustl(rest(i:))
      end do
   end subroutine parse_final_header

   ! spawn the child for one bad case and assert the loud-stop (check_roster scheme)
   subroutine expect_loudstop(cmode, key, who, label)
      character(len=*), intent(in) :: cmode, key, who, label
      character(len=512) :: cmd, ln
      character(len=1024) :: buf
      character(len=1) :: sep
      integer :: stat, v, ios
      buf = ''
      if (windows_host()) then
         sep = '\'
      else
         sep = '/'
      end if
      cmd = 'build'//sep//'check_recorder.exe '//trim(cmode)//' 2>build'//sep//'rec_err.txt'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      open (newunit=v, file='build/rec_err.txt', status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (v, '(a)', iostat=ios) ln
            if (ios /= 0) exit
            buf = trim(buf)//' '//trim(ln)
         end do
         close (v)
      end if
      if (stat == 0) then
         print '(a)', 'FAIL '//label//': expected a loud-stop, the child exited 0'
         nfail = nfail + 1
      else if (index(buf, who) == 0) then
         print '(a)', 'FAIL '//label//': stderr does not name '//trim(who)
         nfail = nfail + 1
      else if (index(buf, key) == 0) then
         print '(a)', 'FAIL '//label//': stderr missing the key "'//trim(key)//'"'
         nfail = nfail + 1
      else
         print '(a)', label//': PASS ('//trim(buf)//')'
      end if
   end subroutine expect_loudstop

   logical function windows_host()
      character(len=32) :: osname
      call get_environment_variable('OS', osname)
      windows_host = (trim(osname) == 'Windows_NT')
   end function windows_host

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

end program check_recorder
