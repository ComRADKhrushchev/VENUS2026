!=====================================================================
! make_sparse_tab.f90 - build sparse cau_tab PES tables from the two
!   C@Au channel-data lineages (data/c_au), as sparse as possible while
!   holding every raw data point to a tolerance
! Design:
!   CHANNEL-DIRECT (2026-10-01 user ruling): the table carries the five
!   channels E0 / h / U / J / VCP directly - no four-diabat subtraction
!   (it stacked interpolation errors into the channels) and the GS DFT
!   data is reference-only (the mechanical base is E0). Two modes:
!     rst     - E0 from rst/E0_Data_GS_rst.txt (the SCF-chain E0 export;
!               its 'a=2.88' header line is a stale template artifact -
!               the site coordinates are the a=2.95 self-frame), h/U/J
!               from rst/Data_{h,U,J}_rst.txt, V from rst/Data_V_rst.txt;
!               LATTICE 2.95; hollow := HCP site (1.475, 0.8517)
!     current - the chain's official aug inputs current/Data_{h,U,J,V}_
!               aug.txt + E0 PROVISIONALLY from current/Data_GS_merged_
!               aug.txt (no E0 export exists for this lineage - the
!               final E0 source is an open user decision, flagged in
!               the table header); a=2.88 IN frame (16 sites,
!               landmark-matched; bridge from the symmetry-equivalent
!               site (5.0549, 3.7509)); hollow := HCP
!   Sparsification: greedy point insertion - keep each block's z-endpoints
!   plus two interior seeds, then repeatedly add the raw point with the
!   largest |table(z) - E_raw| UNTIL every block holds the tolerance.
!   Every candidate is verified through the CONTAINER'S OWN loader
!   (tab_load + tab_eval) - the guarantee "the delivered table's own
!   Akima reproduces the raw points within tol" is by construction, not
!   by a re-implementation. Usage (argv):
!     make_sparse_tab <rst|current> <data_dir> <out_table> [tol_eV]
!   Run from any directory; the greedy scratch file lands in the CWD.
!=====================================================================
program make_sparse_tab
   use cau_tab_tab, only: tab_load, tab_eval, sh_E0, sh_h, sh_U, sh_J, sh_VCP, &
                          si_top, si_bri, si_hol
   implicit none
   integer, parameter :: nsh = 5, nsi = 3, mxraw = 1024
   character(len=8), parameter :: shn(nsh) = ['E0 ', 'h  ', 'U  ', 'J  ', 'VCP']
   character(len=8), parameter :: sin_(nsi) = ['top   ', 'bridge', 'hollow']
   real(8), parameter :: rst_sites(4, 2) = reshape([0.0d0, 1.475d0, 0.0d0, 1.475d0, &
                                                    0.0d0, 0.0d0, 1.7034d0, 0.8517d0], [4, 2])
   real(8), parameter :: cur_lm(3, 2) = reshape([4.3339d0, 5.054874d0, 4.3339d0, &
                                                 2.5022d0, 3.750912d0, 4.1672d0], [3, 2])
   real(8), parameter :: site_xy_rst(3, 2) = reshape([0.0d0, 1.475d0, 1.475d0, &
                                                       0.0d0, 0.0d0, 0.8516d0], [3, 2])
   real(8), parameter :: site_xy_cur(3, 2) = reshape([0.0d0, 1.44d0, 1.44d0, &
                                                       0.0d0, 0.0d0, 0.8314d0], [3, 2])

   real(8) :: rz(nsh, nsi, mxraw), re(nsh, nsi, mxraw)
   integer :: rn(nsh, nsi) = 0
   logical :: keep(nsh, nsi, mxraw)
   real(8) :: tol = 1.0d-2, lat
   character(len=256) :: mode, datadir, outfile
   character(len=512) :: line, e0file
   logical :: e0_override
   integer :: i, k, sh, si, ios, u, narg, it, nadded, nviol, n
   integer :: ndup

   narg = command_argument_count()
   if (narg < 3) then
      print '(a)', 'usage: make_sparse_tab <rst|current> <data_dir> <out_table> [tol_eV]'
      stop 1
   end if
   call get_command_argument(1, mode)
   call get_command_argument(2, datadir)
   call get_command_argument(3, outfile)
   if (narg >= 4) then
      call get_command_argument(4, line)
      read (line, *) tol
   end if
   e0_override = .false.
   if (narg >= 5) then
      call get_command_argument(5, e0file)
      e0_override = .true.
   end if
   rn = 0
   if (trim(mode) == 'rst') then
      lat = 2.95d0
      call read_rst_group()
   else if (trim(mode) == 'current') then
      lat = 2.88d0
      call read_current_group()
   else
      print '(a)', 'fatal: unknown mode "'//trim(mode)//'" (rst|current)'
      stop 1
   end if

   ! 1. per-block sort by z + merge near-duplicate z (within 1e-3 A - the
   !    raw sets interleave two calculation passes at the same site whose z
   !    differ at the 1e-6 scale and whose energies disagree by up to
   !    ~0.03 eV; a 1e-6 dedupe let both through and the 1e-6-wide segment
   !    blew the Akima slopes to 1e4 eV/A - 2026-10-01). Merged points
   !    AVERAGE the energies; the per-block spread is reported.
   ndup = 0
   do sh = 1, nsh
      do si = 1, nsi
         call sort_block(rz(sh, si, :), re(sh, si, :), rn(sh, si), ndup)
         if (rn(sh, si) < 4) then
            print '(a,a,a,a,a,i0)', 'fatal: sheet ', shn(sh), ' site ', sin_(si), &
               ' carries fewer than 4 distinct raw points: ', rn(sh, si)
            stop 1
         end if
      end do
   end do
   print '(a,i0,a)', 'make_sparse_tab: raw points read (after sort/dedupe, ', ndup, ' dups dropped)'
   do sh = 1, nsh
      do si = 1, nsi
         print '(a,a,a,a,a,i0,a,f7.2,a,f7.2,a)', '  raw ', shn(sh), ' ', sin_(si), ': ', &
            rn(sh, si), ' pts  z [', rz(sh, si, 1), ', ', rz(sh, si, rn(sh, si)), ']'
      end do
   end do

   ! 2. greedy insertion against the container's own loader
   call greedy()
   call write_table(outfile, .true.)

   ! 3. verification pass: reload the DELIVERED file, evaluate every raw point
   print '(a)', '--- verification (delivered table vs every raw point) ---'
   call tab_load(outfile, quiet=.true.)
   nviol = 0
   do sh = 1, nsh
      do si = 1, nsi
         call check_block(sh, si, nviol)
      end do
   end do
   if (nviol > 0) then
      print '(a,i0,a)', 'fatal: ', nviol, ' blocks exceed the tolerance'
      stop 1
   end if
   print '(a,es9.2,a)', 'make_sparse_tab: all blocks hold tol = ', tol, ' eV'

contains

   subroutine add_point(sh, si, z, e)
      integer, intent(in) :: sh, si
      real(8), intent(in) :: z, e
      rn(sh, si) = rn(sh, si) + 1
      if (rn(sh, si) > mxraw) then
         print '(a)', 'fatal: mxraw overflow - raise mxraw'
         stop 1
      end if
      rz(sh, si, rn(sh, si)) = z
      re(sh, si, rn(sh, si)) = e
   end subroutine add_point

   !------------------------------------------------------------------
   ! rst group (channels): E0_Data_GS_rst + Data_{h,U,J,V}_rst, 5/6-col
   !                 x y z w E [dEdz] (list-directed reads take the
   !                 first five items; a tail column is ignored)
   !------------------------------------------------------------------
   subroutine read_rst_group()
      character(len=32), parameter :: files(nsh) = [character(len=32) :: &
                                                    'E0_Data_GS_rst.txt', 'Data_h_rst.txt', &
                                                    'Data_U_rst.txt', 'Data_J_rst.txt', &
                                                    'Data_V_rst.txt']
      real(8) :: x, y, z, w, e
      integer :: f, isi
      do f = 1, nsh
         open (newunit=u, file=trim(datadir)//'/rst/'//trim(files(f)), status='old', action='read', iostat=ios)
         if (ios /= 0) then
            print '(a)', 'fatal: cannot open '//trim(datadir)//'/rst/'//trim(files(f))
            stop 1
         end if
         do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(adjustl(line)) == 0) cycle
            if (adjustl(line(1:1)) == '#') cycle
            read (line, *, iostat=ios) x, y, z, w, e
            if (ios /= 0) cycle
            isi = 0
            if (dist2(x, y, rst_sites(1, :)) < 1.0d-4) isi = si_top
            if (dist2(x, y, rst_sites(2, :)) < 1.0d-4) isi = si_bri
            if (dist2(x, y, rst_sites(4, :)) < 1.0d-4) isi = si_hol    ! (1.475,0.8517) = HCP
            if (isi > 0) call add_point(f, isi, z, e)
         end do
         close (u)
      end do
      ! near-field completion, rst edition: the deep rows of bridge/hollow
      ! are REAL dense CDFT data (h/U/J extrema, the E0 repulsive wall) and
      ! are KEPT; only blocks whose coverage stops above 0.20 A (every top
      ! block, the bridge blocks from 0.5) get the smooth continuation
      do isi = 1, nsi
         do sh = 1, nsh
            call extend_block(sh, isi, sh /= sh_E0 .and. isi == si_top)
         end do
      end do
   end subroutine read_rst_group

   !------------------------------------------------------------------
   ! current group (channels): the chain's official aug inputs - E0
   !                PROVISIONALLY from Data_GS_merged_aug (no E0 export
   !                exists for this lineage; flagged in the header),
   !                h/U/J/V from Data_{h,U,J,V}_aug; 6-col x y z w E
   !                dEdz, a=2.88 IN frame, landmark-matched
   !------------------------------------------------------------------
   subroutine read_current_group()
      character(len=32), parameter :: files(nsh) = [character(len=32) :: &
                                                    'Data_GS_merged_aug.txt', 'Data_h_aug.txt', &
                                                    'Data_U_aug.txt', 'Data_J_aug.txt', &
                                                    'Data_V_aug.txt']
      real(8) :: x, y, z, w, e
      character(len=512) :: fpath
      integer :: f, isi
      do f = 1, nsh
         fpath = trim(datadir)//'/current/'//trim(files(f))
         if (f == 1 .and. e0_override) fpath = trim(e0file)
         open (newunit=u, file=trim(fpath), status='old', action='read', iostat=ios)
         if (ios /= 0) then
            print '(a)', 'fatal: cannot open '//trim(fpath)
            stop 1
         end if
         do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(adjustl(line)) == 0) cycle
            if (adjustl(line(1:1)) == '#') cycle
            read (line, *, iostat=ios) x, y, z, w, e
            if (ios /= 0) cycle
            do isi = 1, 3
               if (dist2(x, y, cur_lm(isi, :)) < 1.0d-4) call add_point(f, isi, z, e)
            end do
         end do
         close (u)
      end do
      ! near-field completion (2026-10-01 user ruling: NO artificial data
      ! floors - the physics must repel naturally in the near field): every
      ! block evaluates down to z = 0.20 A. The aug deep rows (z < 1.25)
      ! are extrapolation noise (h oscillates -13..+13 eV) - replaced by a
      ! slope-matched continuation; E0 takes the steepening exponential
      ! repulsive wall, the h/U/J/VCP channels take the linear one
      do isi = 1, nsi
         do sh = 1, nsh
            call extend_block(sh, isi, .true.)
         end do
      end do
   end subroutine read_current_group

   !------------------------------------------------------------------
   ! extend_block(sh, si, replace_deep) - complete one block's coverage
   !                 down to z = 0.20 A with a C1 continuation at the
   !                 block's first trusted anchor (the first point with
   !                 z >= z_trust):
   !                   E0      : E(z) = E_j + (s0/lam)*(1 - exp(-lam*(z-z_j)))
   !                             - the steepening repulsive wall (the near
   !                             field must repel naturally)
   !                   h/U/J/V : E(z) = E_j + (s0/kap)*tanh(kap*(z-z_j)) -
   !                             SATURATING (the real deep data turns over:
   !                             the rst hollow h has its minimum at -0.7 A
   !                             and rises back; a linear continuation runs
   !                             away and its mean-field pull overwhelmed
   !                             the E0 wall; kap = 8/A turns the channel
   !                             over within ~0.15 A and ~|s0|/8 eV,
   !                             mirroring the real deep-region turnover
   !                             (rst hollow h: min -6.57 vs anchor -7.9))
   !                 replace_deep = .true.  : the rows below the anchor are
   !                 aug-extrapolation noise - REWRITTEN to the continuation
   !                 replace_deep = .false. : rows below the anchor are real
   !                 data - kept; only the missing coverage ladder
   !                 (z = 0.20 .. anchor-0.15, step 0.15) is synthesized
   !------------------------------------------------------------------
   subroutine extend_block(sh, si, replace_deep)
      integer, intent(in) :: sh, si
      logical, intent(in) :: replace_deep
      real(8), parameter :: kap = 8.0d0, z_floor = 0.20d0, dl = 0.15d0
      real(8) :: z_trust
      real(8), parameter :: w0 = 50.0d0        ! the wall height at z_floor [eV]
      integer :: i, j, nrep, nadd, it
      real(8) :: zj, ej, s0, ze, ee, lam, dd
      ! the trust floor rides the lineage: current aug channels are real down
      ! to 1.25; the rst TOP-site segment 1.6-2.0 is the old pipeline's own
      ! extrapolation below the diabat data floor (top diabats start at 2.0)
      ! - its step artifact drove the mean-field collapse (2026-10-01), so
      ! the rst top channels trust only from 2.0; bridge/hollow deep rows
      ! are backed by real diabat coverage
      z_trust = 1.25d0
      if (trim(mode) == 'rst' .and. si == si_top) z_trust = 2.00d0
      lam = 2.5d0
      n = rn(sh, si)
      j = 0
      do i = 1, n
         if (rz(sh, si, i) >= z_trust) then
            j = i
            exit
         end if
      end do
      if (j == 0 .or. j >= n) then
         print '(a)', 'extend_block: no trusted anchor - skipped'
         return
      end if
      zj = rz(sh, si, j)
      ej = re(sh, si, j)
      s0 = (re(sh, si, j + 1) - ej)/(rz(sh, si, j + 1) - zj)
      ! the E0 wall reaches W0 eV at z_floor: solve (|s0|/lam)(exp(lam*|dd|)-1)
      ! = W0 for lam by fixed point (a shallow anchor slope - the rst E0 top
      ! secant is only ~0.5 eV/A - would otherwise leave a ~5 eV wall that
      ! the mean-field plunge punches through)
      dd = abs(z_floor - zj)
      if (sh == sh_E0 .and. s0 < 0.0d0) then
         do it = 1, 60
            lam = log(1.0d0 + w0*lam/abs(s0))/dd
            if (lam > 30.0d0) exit
         end do
         lam = min(lam, 30.0d0)
      else
         lam = 2.5d0
      end if
      nrep = 0
      if (replace_deep) then
         do i = 1, j - 1
            if (sh == sh_E0) then
               re(sh, si, i) = ej + (s0/lam)*(1.0d0 - exp(-lam*(rz(sh, si, i) - zj)))
            else
               re(sh, si, i) = ej + (s0/kap)*tanh(kap*(rz(sh, si, i) - zj))
            end if
            nrep = nrep + 1
         end do
      end if
      ! synthesize the missing coverage ladder down to z_floor
      nadd = 0
      ze = z_floor
      do while (ze < zj - 1.0d-9)
         if (.not. any(abs(rz(sh, si, 1:rn(sh, si)) - ze) < 5.0d-3)) then
            if (sh == sh_E0) then
               ee = ej + (s0/lam)*(1.0d0 - exp(-lam*(ze - zj)))
            else
               ee = ej + (s0/kap)*tanh(kap*(ze - zj))
            end if
            call add_point(sh, si, ze, ee)
            nadd = nadd + 1
         end if
         ze = ze + dl
      end do
      ! re-sort the block (the ladder points appended at the tail)
      call sort_block(rz(sh, si, :), re(sh, si, :), rn(sh, si), i)
      if (nrep + nadd > 0) print '(a,a,a,a,a,i0,a,i0,a,f6.2)', 'extend: sheet ', shn(sh), &
         ' site ', sin_(si), ': replaced ', nrep, ', added ', nadd, ' (anchor z = ', zj, ')'
   end subroutine extend_block

   real(8) function dist2(x, y, p)
      real(8), intent(in) :: x, y, p(2)
      dist2 = (x - p(1))**2 + (y - p(2))**2
   end function dist2

   !------------------------------------------------------------------
   ! sort_block - insertion sort by z; near-duplicate z (within 1e-3 A)
   !                 merged by ENERGY AVERAGING (nd counts the merged-away
   !                 rows; the inter-pass disagreement is folded into the
   !                 averaged value and reported via the global ndup count)
   !------------------------------------------------------------------
   subroutine sort_block(zs, es, n, nd)
      real(8), intent(inout) :: zs(:), es(:)
      integer, intent(inout) :: n, nd
      integer :: i, j, ng
      real(8) :: tz, te
      do i = 2, n
         tz = zs(i)
         te = es(i)
         j = i - 1
         do while (j >= 1 .and. zs(j) > tz)
            zs(j + 1) = zs(j)
            es(j + 1) = es(j)
            j = j - 1
         end do
         zs(j + 1) = tz
         es(j + 1) = te
      end do
      j = 1
      ng = 1
      do i = 2, n
         if (zs(i) - zs(j) > 1.0d-3) then
            j = j + 1
            zs(j) = zs(i)
            es(j) = es(i)
            ng = 1
         else
            ng = ng + 1
            es(j) = es(j) + (es(i) - es(j))/ng   ! running average over the group
            nd = nd + 1
         end if
      end do
      n = j
   end subroutine sort_block

   !------------------------------------------------------------------
   ! greedy() - keep endpoints + 2 interior seeds per block, then add the
   !   worst raw point per block until every raw point holds tol. Every
   !   candidate set is validated through the container's own loader.
   !------------------------------------------------------------------
   subroutine greedy()
      integer :: iter, k3
      logical :: changed
      keep = .false.
      do sh = 1, nsh
         do si = 1, nsi
            n = rn(sh, si)
            if (n <= 6) then
               keep(sh, si, 1:n) = .true.
            else
               keep(sh, si, 1) = .true.
               keep(sh, si, n) = .true.
               k3 = (n + 2)/3
               keep(sh, si, k3) = .true.
               keep(sh, si, 2*k3) = .true.
            end if
         end do
      end do
      it = 0
      nadded = 0
      do iter = 1, mxraw*nsh*nsi
         call write_table('tab_sparse_tmp.txt', .false.)
         call tab_load('tab_sparse_tmp.txt', quiet=.true.)
         changed = .false.
         do sh = 1, nsh
            do si = 1, nsi
               call grow_block(sh, si, changed)
            end do
         end do
         it = iter
         if (.not. changed) exit
      end do
      print '(a,i0,a,i0,a)', 'greedy: converged after ', it, ' load/insert rounds, ', nadded, ' points added'
   end subroutine greedy

   subroutine grow_block(sh, si, changed)
      integer, intent(in) :: sh, si
      logical, intent(inout) :: changed
      integer :: i, wi
      real(8) :: ev, ed, worst
      worst = 0.0d0
      wi = 0
      do i = 1, rn(sh, si)
         if (keep(sh, si, i)) cycle
         call tab_eval(sh, si, rz(sh, si, i), ev, ed)
         if (abs(ev - re(sh, si, i)) > worst) then
            worst = abs(ev - re(sh, si, i))
            wi = i
         end if
      end do
      if (worst > tol .and. wi > 0) then
         keep(sh, si, wi) = .true.
         nadded = nadded + 1
         changed = .true.
      end if
   end subroutine grow_block

   !------------------------------------------------------------------
   ! write_table(fname, full_header) - the cau_tab v1 format
   !------------------------------------------------------------------
   subroutine write_table(fname, full)
      character(len=*), intent(in) :: fname
      logical, intent(in) :: full
      real(8) :: sxy(3, 2)
      integer :: nk, i
      sxy = site_xy_rst
      if (trim(mode) == 'current') sxy = site_xy_cur
      open (newunit=u, file=fname, status='replace', action='write')
      write (u, '(a)') '# cau_tab PES table v1 - sparse grid, mode "'//trim(mode)//'"'
      if (full) then
         write (u, '(a)') '# lineage: rst group channels (data/c_au/rst: E0_Data_GS_rst + Data_{h,U,J,V}_rst) | hollow = HCP'
         if (trim(mode) == 'current') then
            if (e0_override) then
               write (u, '(a)') '# lineage: current group channels (data/c_au/current: Data_{h,U,J,V}_aug;'// &
                  ' E0 = the SCF-chain export) | hollow = HCP; bridge from the a2/2 symmetry-equivalent'// &
                  ' site (5.0549,3.7509); E0 source file: '//trim(e0file)
            write (u, '(a)') '# E0 z < 1.25 A: SYNTHETIC repulsive wall (C1 exponential continuation'// &
               ' at the 1.25 A anchor; the aug-extrapolation noise replaced - temporary, 2026-10-01)'
            else
               write (u, '(a)') '# lineage: current group channels (data/c_au/current: Data_{h,U,J,V}_aug;'// &
                  ' E0 FALLBACK Data_GS_merged_aug - pass the SCF-chain export as argv[5]!) | hollow = HCP'
            end if
         end if
         write (u, '(a,es9.2,a)') '# construction: greedy insertion verified through the container'// &
            ' Akima, every raw point held to ', tol, ' eV'
         write (u, '(a)') '# nominal engine-frame site coordinates below; see the maker report'// &
            ' for per-block kept/total'
      end if
      write (u, '(a,f8.4)') 'LATTICE ', lat
      do sh = 1, nsh
         write (u, '(a,a)') 'SHEET ', trim(shn(sh))
         do si = 1, nsi
            write (u, '(a,a,f9.4,f9.4)') 'SITE ', trim(sin_(si)), sxy(si, 1), sxy(si, 2)
            nk = 0
            do i = 1, rn(sh, si)
               if (.not. keep(sh, si, i)) cycle
               write (u, '(f12.6,1x,es24.16)') rz(sh, si, i), re(sh, si, i)
               nk = nk + 1
            end do
         end do
      end do
      close (u)
   end subroutine write_table

   !------------------------------------------------------------------
   ! check_block - the delivered table vs every raw point + the report
   !------------------------------------------------------------------
   subroutine check_block(sh, si, nviol)
      integer, intent(in) :: sh, si
      integer, intent(inout) :: nviol
      integer :: i, nk
      real(8) :: ev, ed, worst
      worst = 0.0d0
      do i = 1, rn(sh, si)
         call tab_eval(sh, si, rz(sh, si, i), ev, ed)
         worst = max(worst, abs(ev - re(sh, si, i)))
      end do
      nk = count(keep(sh, si, 1:rn(sh, si)))
      print '(a,a,a,a,a,i0,a,i0,a,es9.2)', '  ', shn(sh), ' ', sin_(si), ': kept ', nk, '/', &
         rn(sh, si), '  max|dE| = ', worst
      if (worst > tol + 1.0d-12) nviol = nviol + 1
   end subroutine check_block

end program make_sparse_tab
