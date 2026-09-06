!###############################################################################
    !# Subroutines for recording and plotting histograms of impact results
    !# CZZ, 2025.11.30, with thanks to #Github Copilot# for suggestions
    
    !#WARNING：For C atom scatter from Au surface only, DON'T FORGET!!!!
    !#If you use for other systems, please modify accordingly!!!
!###############################################################################
subroutine record_all
      use venus_data
      use venus_params
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     COMMON/HIST_ATOM/AErelS(99999),NBncS(99999),AInitX(99999),&  -- now in venus_data/venus_params
!     AInitY(99999),AZminS(99999),AVfinS(99999),AThtaS(99999)
!     COMMON/HIST_ATOM_NAD/AEphS(99999),AEelS(99999),AEu0S(99999),NCrstS(99999)  -- now in venus_data/venus_params
!     COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX  -- now in venus_data/venus_params
!     COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI  -- now in venus_data/venus_params
    character(len=3)::site
    if(NTZ==1) write(52,'(3A6,10A18)') 'Ntz','Nb','Site','RX0','RY0','Zmin','Vfinal','Erel','Theta', 'Ephi','Eel/f','U0/f'
        call surfcell((/0.0d0,0.0d0/),(/AInitX(NTZ),AInitY(NTZ)/),auro,PI/3,site)
        if(caltyp<0) then
            write(52,'(2I5,A6,7F18.7)')NTZ,NBncS(NTZ),site,AInitX(NTZ),AInitY(NTZ),AZminS(NTZ),AVfinS(NTZ),AErelS(NTZ),AThtaS(NTZ),AEphS(NTZ)
        end if
            
end subroutine
!###############################################################################
subroutine plot_Erel_hist(Erel_min,Erel_max, Nbin)
      use venus_data
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    integer i,Nbin,bin_index
    real*8 Erel_min,Erel_max,dErel
    real*8 Erel_bin(1000),Erel_count(1000)
    character(len=40)::filename
!     COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX  -- now in venus_data/venus_params
!     COMMON/HIST_ATOM/AErelS(99999),NBncS(99999),AInitX(99999),&  -- now in venus_data/venus_params
!        AInitY(99999),AZminS(99999),AVfinS(99999),AThtaS(99999)
!     COMMON/HIST_ATOM_NAD/AEphS(99999),AEelS(99999),AEU0S(99999),NCrstS(99999)  -- now in venus_data/venus_params
    !#-------------------------------------------------------------------
    !# Plot histogram of Erel from all impacts
    
    !# Inputs/~
    !     Erel_min,Erel_max: Energy range
    !     Nbin: Number of bins
    
    !# O. CZZ, 2025.11.30
    !#-------------------------------------------------------------------
    filename = 'HISTO/Erel_hist.dat'
    
    dErel = (Erel_max - Erel_min)/Nbin
    
    do i=1,Nbin
        Erel_bin(i) = Erel_min + (i-0.5d0)*dErel
        Erel_count(i) = 0.d0
    end do
    
    do i=1,NT
        if (AErelS(i) .ge. Erel_min .and. AErelS(i) .lt. Erel_max) then
           bin_index = int((AErelS(i) - Erel_min)/dErel) + 1
           Erel_count(bin_index) = Erel_count(bin_index) + 1.d0
        end if
    end do
    
    open(unit=99,file=filename,status='replace')
    do i=1,Nbin
        write(99,'(2F12.6)') Erel_bin(i), Erel_count(i)
    end do
    write(99,'(A9,F14.8,F12.6)')'TOTAL & P:',SUM(Erel_count(1:Nbin)),SUM(Erel_count(1:Nbin))/dble(NT)
    close(99)
    end subroutine
!###############################################################################
subroutine plot_Erel_impact_pos()
      use venus_data
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    integer i
        real*8 ppp_Q(2)
    character(len=40)::filename
!     COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX  -- now in venus_data/venus_params
!     COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI  -- now in venus_data/venus_params
!     COMMON/HIST_ATOM/AErelS(99999),NBncS(99999),AInitX(99999),&  -- now in venus_data/venus_params
!        AInitY(99999),AZminS(99999),AVfinS(99999),AThtaS(99999)
    !#-------------------------------------------------------------------
    !# Output E_rel and impact XY positions for surface plotting
    !# Writes file 'Erel_pos.dat' with columns: X Y E_rel
    
    !# O. CZZ, 2025.12.01
    !#-------------------------------------------------------------------
    filename = 'HISTO/Erel_pos.dat'
    open(unit=98,file=filename,status='replace')
    do i=1,NT
        call Dist_cell((/AInitX(i),AInitY(i)/),(/0.0d0,0.0d0/),auro,PI/3,ppp_Q)
        write(98,'(4F12.6)')ppp_Q, AErelS(i),AZminS(i)
    end do
    close(98)
end subroutine
!###############################################################################
subroutine plot_theta_hist(theta_min,theta_max,Nbin)
      use venus_data
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    integer i,Nbin,bin_index
    real*8 theta_min,theta_max,dtheta
    real*8 theta_bin(1000),theta_count(1000)
    character(len=64)::filename
!     COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX  -- now in venus_data/venus_params
!     COMMON/HIST_ATOM/AErelS(99999),NBncS(99999),AInitX(99999),&  -- now in venus_data/venus_params
!        AInitY(99999),AZminS(99999),AVfinS(99999),AThtaS(99999)
    !#-------------------------------------------------------------------
    !# Plot histogram of scattering angle from all impacts
    !#
    !# Inputs/~
    !#     theta_min, theta_max: Angle range (degrees)
    !#     Nbin: Number of bins
    !#
    !# O. CZZ, 2025.12.01
    !#-------------------------------------------------------------------
    filename = 'HISTO/Theta_hist.dat'
    dtheta = (theta_max - theta_min)/Nbin
    do i=1,Nbin
        theta_bin(i) = theta_min + (i-0.5d0)*dtheta
        theta_count(i) = 0.d0
    end do
    do i=1,NT
        if (AThtaS(i) .ge. theta_min .and. AThtaS(i) .lt. theta_max) then
           bin_index = int((AThtaS(i) - theta_min)/dtheta) + 1
           theta_count(bin_index) = theta_count(bin_index) + 1.d0
        end if
    end do
    open(unit=99,file=filename,status='replace')
    do i=1,Nbin
        write(99,'(2F12.6)') theta_bin(i), theta_count(i)
    end do
    close(99)
    end subroutine
    
!###############################################################################
subroutine plot_site_histograms(Erel_min,Erel_max,NEbin,theta_min,theta_max,NTbin,Qrefx,Qrefy)
      use venus_data
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    !#-------------------------------------------------------------------
    !# Produce separate histograms (E_rel and theta) for each site:
    !#   sites columns order: 1=top,2=brg,3=hcp,4=fcc
    !# Writes files:
    !#   Erel_hist_top.dat, Erel_hist_brg.dat, Erel_hist_hcp.dat, Erel_hist_fcc.dat
    !#   Theta_hist_top.dat, ...
    !#
    !# Inputs/~
    !#   Erel_min,Erel_max, NEbin    : energy range and bins
    !#   theta_min,theta_max,NTbin   : angle range and bins (same units as AThtaS)
    !#   le, skew, Qrefx,Qrefy       : parameters passed to `surfcell`
    !#
    !# O. CZZ, 2025.12.01
    !#-------------------------------------------------------------------
    integer NEbin,NTbin,countall(4)
    integer i, isite, bin, idx, unitE, unitT
    integer countsE(4,1000), countsT(4,1000)
    real*8 Erel_min,Erel_max,theta_min,theta_max
    real*8 dErel, dtheta
    real*8 Ecent(1000), Thetacent(1000)
    real*8 Q_i(2), Q_ref(2)
    character(len=3) :: site
    character(len=64) :: fname
    character(len=4), dimension(4) :: sname
!     COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI  -- now in venus_data/venus_params
!     COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX  -- now in venus_data/venus_params
!     COMMON/HIST_ATOM/AErelS(99999),NBncS(99999),AInitX(99999),&  -- now in venus_data/venus_params
!        AInitY(99999),AZminS(99999),AVfinS(99999),AThtaS(99999)
    sname = (/'top','brg','hcp','fcc'/)
    ! guard no-more-than-1000 bins
    if (NEbin .gt. 1000) then
       NEbin = 1000
    end if
    if (NTbin .gt. 1000) then
       NTbin = 1000
    end if
    dErel = (Erel_max - Erel_min)/dble(NEbin)
    dtheta = (theta_max - theta_min)/dble(NTbin)
    do bin = 1, NEbin
        Ecent(bin) = Erel_min + (dble(bin)-0.5d0)*dErel
    end do
    do bin = 1, NTbin
        Thetacent(bin) = theta_min + (dble(bin)-0.5d0)*dtheta
    end do
    countsE = 0
    countsT = 0
    Q_ref(1) = Qrefx
    Q_ref(2) = Qrefy
    do i = 1, NT
        Q_i(1) = AInitX(i)
        Q_i(2) = AInitY(i)
        call surfcell(Q_ref, Q_i, auro, PI/3, site)
        if (site == 'top') then
           idx = 1
        else if (site == 'brg') then
           idx = 2
        else if (site == 'hcp') then
           idx = 3
        else
           idx = 4
        end if
        countall(idx) = countall(idx) + 1
        
        ! energy histogram for this site
        if (AErelS(i) .ge. Erel_min .and. AErelS(i) .lt. Erel_max) then
           bin = int((AErelS(i) - Erel_min)/dErel) + 1
           if (bin .ge. 1 .and. bin .le. NEbin) countsE(idx,bin) = countsE(idx,bin) + 1
        end if
        ! theta histogram for this site
        if (AThtaS(i) .ge. theta_min .and. AThtaS(i) .lt. theta_max) then
           bin = int((AThtaS(i) - theta_min)/dtheta) + 1
           if (bin .ge. 1 .and. bin .le. NTbin) countsT(idx,bin) = countsT(idx,bin) + 1
        end if
    end do
    ! write per-site files
    do isite = 1, 4
        fname = 'HISTO/Erel_hist_'//trim(sname(isite))//'.dat'
        unitE = 10 + isite
        open(unit=unitE, file=fname, status='replace')
        do bin = 1, NEbin
            write(unitE,'(2F12.6)') Ecent(bin), dble(countsE(isite,bin))
        end do
        sum1 = SUM(countsE(isite,1:NEbin))
        write(unitE,'(A9,F14.8,F12.6)')'TOTAL & P:',sum1,sum1/countall(isite)
        close(unitE)
        fname = 'HISTO/Theta_hist_'//trim(sname(isite))//'.dat'
        unitT = 20 + isite
        open(unit=unitT, file=fname, status='replace')
        do bin = 1, NTbin
            write(unitT,'(2F12.6)') Thetacent(bin), dble(countsT(isite,bin))
        end do
        close(unitT)
    end do
    end subroutine
!###############################################################################  
subroutine plot_Erel_by_bounce(Erel_min,Erel_max,Nbin,merge_at)
      use venus_data
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    integer, parameter :: MAXBIN = 1000, MAXB = 200
    integer i, bin, bounce_idx, nbins_bounce, bin_idx
    integer counts( MAXB, MAXBIN )
    integer Nbin, merge_at
    integer nb
    real*8 Erel_min, Erel_max, dErel
    real*8 Ecenter(MAXBIN)
    character(len=64) :: fname
    integer unit
!     COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX  -- now in venus_data/venus_params
!     COMMON/HIST_ATOM/AErelS(99999),NBncS(99999),AInitX(99999),&  -- now in venus_data/venus_params
!        AInitY(99999),AZminS(99999),AVfinS(99999),AThtaS(99999)
    !#-------------------------------------------------------------------
    !# 按反弹次数统计 E_rel 分布（每个反弹次数一个直方图）
    !#
    !# 输入：
    !#   Erel_min, Erel_max : 能量范围
    !#   Nbin               : 每个直方图的箱数（<=1000）
    !#   merge_at           : 合并阈值；反弹次数 >= merge_at 的轨迹合并到最后一格
    !#
    !# 输出：
    !#   在工作目录生成文件：
    !#     Erel_hist_bounce_0.dat, ..., Erel_hist_bounce_{merge_at-1}.dat,
    !#     Erel_hist_bounce_ge{merge_at}.dat
    !#   每个文件两列：E_center  count
    !#
    !# 用法示例： call plot_Erel_by_bounce(0.d0,10.d0,100,10)
    !# O. CZZ, 2025.12.01
    !#-------------------------------------------------------------------
    ! 参数检查和界限
    if (Nbin .gt. MAXBIN) then
       write(*,*) 'Nbin too large, limit to ', MAXBIN
       Nbin = MAXBIN
    end if
    if (merge_at .lt. 1) merge_at = 1
    if (merge_at .gt. MAXB-1) then
       write(*,*) 'merge_at too large, limit to ', MAXB-1
       merge_at = MAXB-1
    end if
    nbins_bounce = merge_at + 1   ! 1..merge_at 对应 bounce=0..merge_at-1，第 nbins_bounce 对应 >=merge_at
    ! 初始化箱中心与计数
    dErel = (Erel_max - Erel_min)/dble(Nbin)
    do bin = 1, Nbin
        Ecenter(bin) = Erel_min + (dble(bin)-0.5d0)*dErel
    end do
    counts = 0
    ! 遍历所有轨迹，按照反弹次数分类并计数能量箱
    do i = 1, NT
        nb = NBncS(i)
        if (nb .lt. 0) cycle
        if (nb .ge. merge_at) then
           bounce_idx = nbins_bounce
        else
           bounce_idx = nb + 1
        end if
        if (AErelS(i) .ge. Erel_min .and. AErelS(i) .lt. Erel_max) then
           bin_idx = int((AErelS(i) - Erel_min)/dErel) + 1
           if (bin_idx .ge. 1 .and. bin_idx .le. Nbin) then
              counts(bounce_idx, bin_idx) = counts(bounce_idx, bin_idx) + 1
           end if
        end if
    end do
    ! 写出每个反弹次数对应的文件
    do bounce_idx = 1, nbins_bounce
        if (bounce_idx .lt. nbins_bounce) then
           write(fname,'(A,I0,A)') 'HISTO/Erel_hist_bounce_', bounce_idx-1, '.dat'
        else
           write(fname,'(A,I0,A)') 'HISTO/Erel_hist_bounce_ge', merge_at, '.dat'
        end if
        unit = 100 + bounce_idx
        open(unit=unit, file=fname, status='replace')
        do bin = 1, Nbin
            write(unit,'(2F12.6)') Ecenter(bin), dble(counts(bounce_idx,bin))
        end do
        close(unit)
    end do
end subroutine
!###############################################################################    
subroutine surfcell(Q_ref,Q_i,le,skew_d,site)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    !#-------------------------------------------------------------------
    !# Attribute random impact position to a certain surf. site
    !# Like Wigner-Setz Cell
    ! Inputs/~
    !     Q_ref: Any surface atomic position,
    !     Q_i: Impact position randomly selected
    !     le,skew_d: Surface cell length, angle
    ! Outputs/~
    !     Site: Attributed site, fcc/hcp/brg/top
    !O. 2025.12, CZZ
    !#-------------------------------------------------------------------
    real*8 Q_i(2),Q_ref(2),le,skew_d,PBC_Q(2)
    real*8 cell(2,4),edge_cent(2,5),fcc(2),hcp(2)
    real*8 dist_all(4+5+1+1) !# top+bridge+fcc+hcp
    integer i,min_index
    character(len=3)::site
    !#-------------------------------
    !#        Au3——c4——Au4
    !#       / \  fcc  /
    !#      c2   c5   c3
    !#     / hcp  \  /
    !#    Au1——c1——Au2
    !#  (Qref)
    !#------------------------------- 
    
    cell(:,1) = (/0d0,0d0/)
    cell(:,2) = (/le, 0d0/)
    cell(:,3) = (/le * cos(skew_d), le * sin(skew_d)/)
    cell(:,4) = (/le * (1.d0 + cos(skew_d)), le * sin(skew_d)/)
    fcc = (cell(:,2) + cell(:,4) + cell(:,3))/3
    hcp = (cell(:,1) + cell(:,3) + cell(:,2))/3
    call Dist_cell(Q_i,Q_ref,le,skew_d,PBC_Q)
    
    !top sites
    do i=1,4
        dist_all(i) = norm2(PBC_Q - cell(:,i))
    end do
    !bridge sites
    do i=1,4
        edge_cent(:,i) = (cell(:,i) + cell(:,mod(i,4)+1))/2
        dist_all(4+i) = norm2(PBC_Q - edge_cent(:,i))
    end do
        edge_cent(:,5) = (cell(:,1) + cell(:,3))/2
        dist_all(9) = norm2(PBC_Q - edge_cent(:,5))        
    dist_all(10) = norm2(PBC_Q - hcp) !hcp site
    dist_all(11) = norm2(PBC_Q - fcc) !fcc site
    
    ! Find the minimum distance
    min_index = minloc(dist_all,dim=1)
    if (min_index .le. 4) then
       site = 'top'
    else if (min_index .le. 9) then
       site = 'brg'
    else if (min_index .eq. 10) then
       site = 'hcp'
    else
       site = 'fcc'
    end if
end subroutine 
!###############################################################################
subroutine  Dist_cell(Q_i,Q_ref,le,skew_d,vec_r)
    implicit real*8 (a-h,o-z)
!     COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI  -- now in venus_data/venus_params
    !#-------------------------------------------------------------------
    !# PBC position for surface cells
    ! Inputs/~
    !     Q_i(2): Initial impact position
    !     len,skew: Surface cell length, angle
    !     Q_ref(2): Reference surface atomic position
     
    ! Outputs/~
    !    r(2): Distance vector between Q_i and Au1
     
    !CZZ, 2025.12
    !#-------------------------------------------------------------------
	real*8 Q_i(2),Q_ref(2),Q_op(2),vec_r(2),le,skew_d,a1,a2
	real*8 latvec_x(2),latvec_y(2),dist_vec(9)
	integer :: n1, n2
        latvec_x = (/le,0d0/)
        latvec_y = (/le * cos(SKEW_D),le * sin(SKEW_D)/)
        a2 = (Q_i(2)-Q_ref(2))/(le * sin(SKEW_D))
        a1 = (Q_i(1)-Q_ref(1) - a2 * le * cos(SKEW_D))/le
        
        if(a1 > 0) n1 = int(a1,4)
        if(a1 <= 0) n1 = int(a1,4) - 1
        if(a2 > 0) n2 = int(a2,4)
        if(a2 <= 0) n2 = int(a2,4) - 1
        Q_op = n1*latvec_x + n2*latvec_y + Q_ref
		vec_r =  Q_i - Q_op             
end subroutine
