# -*- coding: utf-8 -*-
p = r"C:\Users\13614\Desktop\MAIN\01_WORK\dev\venus2026\VENUS2026\scripts\make_sparse_tab.f90"
s = open(p, encoding='utf-8').read()

s = s.replace("""! make_sparse_tab.f90 - build sparse cau_tab PES tables from the two
!   C@Au raw-data lineages (data/c_au), as sparse as possible while
!   holding every raw DFT point to a tolerance
! Design:
!   Two modes ('rst' | 'current') read the group's raw files into the
!   six-sheet x three-site blocks (sheet map: GS / 1D / 3P / mC / pC /
!   VCP; sites top / bridge / hollow):
!     rst     - rst/Data_{GS,S,T,A,C,V}_rst.txt, 5-col x y z w E,
!               a=2.95 frame; hollow := HCP site (1.475, 0.8517)
!     current - cluster/pes_{S,T,A,C}* (z E label) for the four diabats,
!               cluster/Data_GS_merged.txt (landmark-matched; bridge from
!               the a2/2 symmetry-equivalent site (5.0549, 3.7509)),
!               cluster/Data_V.txt (CDFT coupling scan, z_rel heights);
!               a=2.88 frame; hollow := HCP""",
"""! make_sparse_tab.f90 - build sparse cau_tab PES tables from the two
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
!               site (5.0549, 3.7509)); hollow := HCP""")

s = s.replace("""   integer, parameter :: nsh = 6, nsi = 3, mxraw = 512
   character(len=8), parameter :: shn(nsh) = ['GS ', '1D ', '3P ', 'mC ', 'pC ', 'VCP']""",
"""   integer, parameter :: nsh = 5, nsi = 3, mxraw = 1024
   character(len=8), parameter :: shn(nsh) = ['E0 ', 'h  ', 'U  ', 'J  ', 'VCP']""")

s = s.replace("""   use cau_tab_tab, only: tab_load, tab_eval, sh_GS, sh_1D, sh_3P, sh_mC, sh_pC, sh_VCP, &
                          si_top, si_bri, si_hol""",
"""   use cau_tab_tab, only: tab_load, tab_eval, sh_E0, sh_h, sh_U, sh_J, sh_VCP, &
                          si_top, si_bri, si_hol""")

old_rst = s[s.index("   !------------------------------------------------------------------\n   ! rst group:"):s.index("   real(8) function dist2")]
new_rst = """   !------------------------------------------------------------------
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
   end subroutine read_rst_group

"""
s = s.replace(old_rst, new_rst)

old_cur = s[s.index("   !------------------------------------------------------------------\n   ! current group:"):s.index("   real(8) function dist2")]
new_cur = """   !------------------------------------------------------------------
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
      integer :: f, isi
      do f = 1, nsh
         open (newunit=u, file=trim(datadir)//'/current/'//trim(files(f)), status='old', action='read', iostat=ios)
         if (ios /= 0) then
            print '(a)', 'fatal: cannot open '//trim(datadir)//'/current/'//trim(files(f))
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
   end subroutine read_current_group

"""
s = s.replace(old_cur, new_cur)

s = s.replace("""         write (u, '(a)') '# lineage: rst group (a=2.95 CDFT, data/c_au/rst) | hollow = HCP site'
         if (trim(mode) == 'current') write (u, '(a)') &
            '# lineage: current group (cluster TSO diabats + GS_merged + CDFT V, data/c_au/cluster) |'// &
            ' hollow = HCP; bridge from the a2/2 symmetry-equivalent GS site (5.0549,3.7509)'""",
"""         write (u, '(a)') '# lineage: rst group channels (data/c_au/rst: E0_Data_GS_rst + Data_{h,U,J,V}_rst) | hollow = HCP'
         if (trim(mode) == 'current') write (u, '(a)') &
            '# lineage: current group channels (data/c_au/current: Data_{h,U,J,V}_aug; E0 PROVISIONALLY'// &
            ' Data_GS_merged_aug - the final E0 source is an open decision) | hollow = HCP;'// &
            ' bridge from the a2/2 symmetry-equivalent site (5.0549,3.7509)'""")
open(p, 'w', encoding='utf-8').write(s)
print("maker patched OK")
