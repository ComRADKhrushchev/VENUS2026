      SUBROUTINE GWRITE
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer :: i_nearest
!
!     Write per-step trajectory diagnostics to output files.
!     Called at each output step (every NIP integration steps)
!     and at trajectory start/end.
!
!     GWRITE_LEVEL controls verbosity:
!       0 = minimal:  C atom Q/V, total T/V/H
!       1 = medium:   + PBC nearest-Au structural diagnostic
!       2 = high:     + all-atom Q/P (DEFAULT, backward compatible)
!       5 = maximum:  + per-atom total forces
!
    1 FORMAT(2X,'THE CYCLE COUNT IS:',I14,16X,'TIME:',F12.3)
    2 FORMAT(2X,'KINETIC ENERGY: ',1PE17.9,'    POTENTIAL ENERGY: ',E17.9/2X,'TOTAL ENERGY:   ',E17.9)
    3 FORMAT(19X,'Q',37X,'P')
    4 FORMAT(F11.6,2F12.6,2X,3F12.6)
    9 FORMAT(1X,'XXXXXXXXXXXXXXXXXXXXXXXX TRAJECTORY NUMBER ',I4,' XXXXXXXXXXXXXXXXXXXXXXXXX')
   10 FORMAT(2X,'THE CURRENT RANDOM NUMBER IS: ',8I4,' BASE 256')
   17 FORMAT(1X,3F11.7)
!
      CNC=NC
      TI=TIME*CNC
      WRITE(6,9)NTZ
      WRITE(6,1)NC,TI

      IF (NSELT.EQ.2.OR.NSELT.EQ.3) WRITE(6,10)(ISEED0(9-I),I=1,8)
      WRITE(6,2)T,V,H

      if (nsel.eq.1) then
         n=natomb(1)-nrgd
         tb=0.d0
         do i=1,n
           j=lb(1,i)
           j3=3*j
           j2=j3-1
           j1=j2-1
           tb=tb+(p(j1)**2+p(j2)**2+p(j3)**2)/w(j)
         enddo
         tempinit=tb/(3.0d0*dble(n)*0.00198717d0*c1)
         write(6,*)'system temperature=',tempinit
      endif
!
!     Per-trajectory graphics / electronic diagnostic file (fort.1NNN)
!
      IF (NCOOR.EQ.1.AND.NSELT.NE.-2) THEN
         ! --- Always written: trajectory banner (fort.8) ---
         WRITE(8,9)NTZ
         WRITE(8,17)(Q(I),I=1,NATOMS*3)

         ! --- Level 0 (always): step header + total energies + C atom basic ---
         write(1000+NTZ,'(A,I8,A,F10.3,A)')'--- step',NC,' t(fs)=',NC*0.1,' ---'
         write(1000+NTZ,'(A6,F14.5,A6,F14.5,A6,F14.5)') &
            'E0(eV)',E0/23.0605,' T(eV)',T/23.0605,' H(eV)',H/23.0605

         write(1000+NTZ,'(A,3F14.5)')'Q(C)=',Q(1), Q(2), Q(3)
         write(1000+NTZ,'(A,3F14.5)')'V(C)=',P(1)/W(1), P(2)/W(1), P(3)/W(1)

         ! ============================================================
         ! Level >= 1: forces on C, electronic populations, Hubbard, PBC
         ! ============================================================
         IF (GWRITE_LEVEL >= 1) THEN


            ! --- PBC: distance from C to nearest Au (structural diagnostic) ---
            xmin_s = QZB(1, 1); xmax_s = QZB(1, 1)
            ymin_s = QZB(1, 2); ymax_s = QZB(1, 2)
            do j = 2, NATOMB(1)
                xj = QZB(1, 3*j-2); yj = QZB(1, 3*j-1)
                if (xj .lt. xmin_s) xmin_s = xj
                if (xj .gt. xmax_s) xmax_s = xj
                if (yj .lt. ymin_s) ymin_s = yj
                if (yj .gt. ymax_s) ymax_s = yj
            end do
            box_sx = xmax_s - xmin_s + 6.0d0
            box_sy = ymax_s - ymin_s + 6.0d0
            r_min = 1.0d30
            i_nearest = 0
            do i = 2, NATOMS
                dx_i = Q(1) - Q(3*i-2)
                dy_i = Q(2) - Q(3*i-1)
                dz_i = Q(3) - Q(3*i)
                dx_i = dx_i - box_sx * nint(dx_i / box_sx)
                dy_i = dy_i - box_sy * nint(dy_i / box_sy)
                r_i = sqrt(dx_i*dx_i + dy_i*dy_i + dz_i*dz_i)
                if (r_i .lt. r_min) then
                    r_min = r_i
                    i_nearest = i
                    dx_min = dx_i; dy_min = dy_i; dz_min = dz_i
                end if
            end do
            write(1000+NTZ,'(A12,I5,A6,F8.4)')'Au_nearest=',i_nearest, &
                ' r_min', r_min
            write(1000+NTZ,'(A12,3F10.4)')'Au_xyz(Q)= ', &
                Q(3*i_nearest-2), Q(3*i_nearest-1), Q(3*i_nearest)
            write(1000+NTZ,'(A12,3F10.4)')'sep_dx,dy,dz=', dx_min, dy_min, dz_min

         END IF  ! GWRITE_LEVEL >= 1

         ! ============================================================
         ! Level >= 2: all atoms Q/P and F_tot (slab coordinates/forces)
         ! ============================================================
         IF (GWRITE_LEVEL >=2) THEN

            ! Q and P for all atoms
            do i=1,NATOMS
                if (W(i) < 50.0d0) then
                    write(1000+NTZ,'(A3,6f14.5)')'C  ',Q(3*i-2:3*i),P(3*i-2:3*i)
                else
                    write(1000+NTZ,'(A3,6f14.5)')'Au ',Q(3*i-2:3*i),P(3*i-2:3*i)
                end if
            end do

         END IF  ! GWRITE_LEVEL >= 2
         
         ! ============================================================
         ! Level >= 5: Forces
         ! ============================================================
        IF (GWRITE_LEVEL >= 5) THEN
            ! Per-atom total forces (all atoms)
            do i=1,NATOMS
                write(1000+NTZ,'(A,I4,A,3E14.5)')'F_tot(',i,')= ',PDOT(3*i-2:3*i)
            end do
        end IF 

         
        
      ENDIF
!
!     Write coordinates and momenta to stdout (fort.6)
!
      IF (NFQP.NE.0) THEN
         WRITE(6,3)
         J=1
         DO iL=1,NATOMA(1)
            M=J+2
            WRITE(6,4)(Q(I),I=J,M),(P(I),I=J,M)
            J=J+3
         ENDDO
      ENDIF

      CALL FLUSH(6)
      RETURN
      END
