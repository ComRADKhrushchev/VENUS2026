# -*- coding: utf-8 -*-
p = r"C:\Users\13614\Desktop\MAIN\01_WORK\dev\venus2026\VENUS2026\scripts\make_sparse_tab.f90"
s = open(p, encoding='utf-8').read()
assert "read_current_group" not in s, "current group already present"
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

   real(8) function dist2"""
s = s.replace("   real(8) function dist2", new_cur, 1)
open(p, 'w', encoding='utf-8').write(s)
print("current group re-inserted")
