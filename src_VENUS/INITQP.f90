subroutine INITQP(WW, A, C, AM, WT, EINT, EROTS, AI, EROT, PHASE, N, NM, NPHASE, IPHASE)
  use venus_params
  use venus_data, only: AI_H => AI, L => LL
  use venus_data
  use venus_input, only: ENU, EDELTU, ENL, EDELTL
  implicit none
  real(8),  intent(in)    :: WW(NDA3), A(NDA3), C(NDA3yf,NDA3yf)
  real(8),  intent(inout) :: AM(4), EINT, EROT
  real(8),  intent(in)    :: WT, EROTS, AI(3), PHASE(5)
  integer,  intent(in)    :: N, NM, NPHASE, IPHASE(5)
  real(8)  :: QCM(3), VCM(3), COOR(NDA3), DCOOR(NDA3)
  real(8)  :: ESEL, DUM1, DUM2, DUM3, DUM, RAND, DDD, SDUM
  integer  :: NUMP, I, J, K, II, JJ, K3, IFLAG, ISEED
  logical  :: matched
  real(8), external :: RAND0
!
!         INITIALIZE COORDINATES AND MOMENTA FROM NORMAL MODE
!         PARAMETERS (FREQUENCY, AMPLITUDE,...)
!
5 format('  A-B INTERACTION ENERGY WHEN ENTERING INITQP =', &
        1pe18.9,' KCAL/MOL'/)
15 format(15x,'INTERNAL ENERGY =',1pe18.9,' KCAL/MOL')
25 format(/,10x,'CHOSEN:  EROT =',f7.3,' KCAL/MOL',/, &
        10x,'JX,JY,JZ =',3d13.5,' H-BAR',/)
35 format(9x,'VIBRATIONAL ANGULAR MOMENTUM =',1pe18.9,' H-BAR')
45 format(/,15x,'ONLY KINETIC ENERGY FOR FIRST ',i3,' MODES')
!
  ESEL = EINT
!
!         CALCULATE THE TOTAL ENERGY, WHICH IS THE REFERENCE ENERGY
!         (EZERO) WITH RESPECT TO ADDING EINT.
!
!      CALL DVDQ        ! by bin, 2016/10/02
  call ENERGY_1
  EZERO = H
  write(6,5) EZERO
!
!         SET COUNTER FOR NUMBER OF SCALING ATTEMPTS
!
  NSCALE = 0
!
!         SET IFLAG WHICH IS USED FOR LOCAL MODE EXCITATION BETWEEN
!         ATOMS NONI AND NONJ
!
  IFLAG = 0
  if (NACTA == 4) call LMODE(1, ENU, EDELTU, ENL, EDELTL)
!
!         STORE THE ANGULAR MOMENTUM VECTOR FROM SELECT
!
  DUM1 = AM(1)
  DUM2 = AM(2)
  DUM3 = AM(3)
!
!         CALCULATE NORMAL MODE COORDINATES AND VELOCITIES
!
!
!-->    modified by Bin, 09/14/2015
  resample: do
     NSCALE = 0
!-->    end
!         SET COOR(I) AND DCOOR(I) FOR MODES WHICH ONLY RECEIVE
!         KINETIC ENERGY SPECIFIED BY NUMP
!
     NUMP = 0
     if (NUMP >= 1) then
        do I = 1, NUMP
           COOR(I) = 0.D0
           DCOOR(I) = -WW(I)*A(I)
        end do
     end if
!
!         FOR OTHER MODES
!
!         TO SET THE VIBRATIONAL ANGULAR MOMENTUM FOR DEGENERATE BENDS,
!         THE PHASE OF THE SECOND MODE OF THE DEGENERATE PAIR IS SHIFTED
!         FROM THAT OF THE FIRST MODE BY THE ANGLE "PHASE".
!
     do I = 1+NUMP, NM
        matched = .false.
        if (NPHASE > 0) then
           do J = 1, NPHASE
              if (IPHASE(J) == I) then
                 DUM = DUM+PHASE(J)
                 matched = .true.
                 exit
              end if
           end do
        end if
        if (.not. matched) then
           RAND = RAND0(ISEED)
           DUM = TWOPI*RAND
        end if
        COOR(I) = A(I)*cos(DUM)
        DCOOR(I) = -WW(I)*A(I)*sin(DUM)
     end do
!
!         TRANSFORM FROM NORMAL MODE TO CARTESIAN COORDINATES AND VELOCIT
!
     do II = 1, N
        do K = 1, 3
           JJ = 3*II+1-K
           J = 3*L(II)+1-K
           Q(J) = 0.0D0
           P(J) = 0.0D0
           do I = 1, NM
              Q(J) = Q(J)+C(JJ,I)*COOR(I)
              P(J) = P(J)+C(JJ,I)*DCOOR(I)
           end do
           P(J) = P(J)*W(L(II))
           Q(J) = Q(J)+QZ(J)
        end do
     end do
!
!         CALCULATE CENTER OF MASS COORDINATES QQ AND MOMENTA PP
!
     inner: do
        call CENMAS(WT, QCM, VCM, N)
!
!         MOVE PP ARRAY TO P ARRAY AND QQ ARRAY TO Q ARRAY
!
        do I = 1, N
           J = 3*L(I)+1
           do K = 1, 3
              Q(J-K) = QQ(J-K)
              P(J-K) = PP(J-K)
           end do
        end do
!
!         ADD ANGULAR MOMENTUM VECTOR FROM SELECT TO THE MOLECULE.
!         CALCULATE THE REQUIRED ANGULAR VELOCITY AND ADD IT TO THE
!         MOLECULE.
!
!         IF NPHASE > 0 THERE IS VIBRATIONAL ANGULAR MOMENTUM ABOUT THE
!         X-AXIS OF A LINEAR MOLECULE.  THIS IS ADDED BY SETTING THE
!         QUANTUM NUMBERS AND PHASE FOR THE DEGENERATE BENDS.  THIS
!         ANGULAR MOMENTUM IS NOT SPURIOUS AND SHOULD NOT BE SUBTRACTED.
!
        call ROTN(AM, EROT, N)
        AM(1) = DUM1-AM(1)
        AM(2) = DUM2-AM(2)
        AM(3) = DUM3-AM(3)
        NAM = 1
        call ROTN(AM, EROT, N)
        NAM = 0
        WX = -WX
        WY = -WY
        WZ = -WZ
        if (NPHASE > 0) then
           WX = 0.0D0
           DUM = AM(1)/C7
           write(6,35) DUM
        end if

        call ANGVEL(N)
!
!         SCALE COORDINATES AND MOMENTA TO FIT THE TOTAL ENERGY
!         THE INITIAL CONDITION IS ACCEPTED IF THE CALCULATED AND
!         DESIRED ENERGY AGREE TO WITHIN 0.1 PER-CENT.
!
!         EINT IS THE SUM OF THE SELECTED INTERNAL VIBRATIONAL AND
!         ROTATIONAL ENERGIES.  DO NOT SCALE IF THE SELECTED VIBRATIONAL
!         ENERGY IS ZERO.
!
        if (IFLAG /= 1) then
           DDD = EINT-EROTS
           if (DDD /= 0.0D0) then
!-->    modified by Bin, 4/29/2014
!-->    shift z coorindates up in order to calculate the potential energy
              if (NSURF > 0) then
                 J = N
                 do I = 1, J
                    K3 = 3*I
                    Q(K3) = Q(K3)+ZASYM
                 end do
              end if
!-->    end Bin's changes
!            CALL DVDQ  ! by bin, 2016/10/02
              call ENERGY_1
!-->    modified by Bin, 4/29/2014
!-->    shift z coorindates back after the potential energy
              if (NSURF > 0) then
                 J = N
                 do I = 1, J
                    K3 = 3*I
                    Q(K3) = Q(K3)-ZASYM
                 end do
              end if
!-->    end Bin's changes

              ESEL = H-EZERO
              SDUM = abs(EINT-ESEL)/EINT
!
!         TEST TO SEE IF CALCULATED ENERGY IS OUT OF RANGE FOR ACCURATE
!         SCALING. IF SO, ONLY ADD KINETIC ENERGY TO AN ADDITIONAL MODE
!
!.........Modified by Bin 12/18/2013
!.........The following four lines were commented by Bin in order to match the Venus96's results of Ar-H2O collision
!.........Their influence on the entire trajectory calculations is not clear
!            IF (SDUM.GE.0.1D0 .AND. NUMP.LT.NM) THEN
!               NUMP=NUMP+1
!               GOTO 48
!            ENDIF
!.........End

              write(6,15) ESEL
              if (SDUM >= 0.001D0) then
                 NSCALE = NSCALE+1
!               IF (NSCALE.GT.50) STOP
!.........Modified by Bin 09/14/2015
                 if (NSCALE > 50) then
                    write(6,*) 'IMPROPER RANDOM NUMBER, USE THE NEXT ONE'
                    cycle resample
                 end if
!.........end
                 SDUM = sqrt(EINT/ESEL)
                 do I = 1, N
                    J = 3*L(I)+1
                    do K = 1, 3
                       P(J-K) = P(J-K)*SDUM
                       Q(J-K) = (Q(J-K)-QZ(J-K))*SDUM+QZ(J-K)
                    end do
                 end do
                 cycle inner
              end if
              if (NUMP > 0) write(6,45) NUMP
           end if
           if (NACTA /= 4) exit inner
        end if
!
!         CHOOSE CONDITIONS FOR LOCAL MODE(NACT=4)
!
        call LMEXCT
        IFLAG = 1
        if (JFLAG == 0) cycle inner
!
!         CALCULATE THE TOTAL ENERGY
!
!-->    modified by Bin, 4/29/2014
!-->    shift z coorindates up in order to calculate the potential energy
        if (NSURF > 0) then
           J = N
           do I = 1, J
              K3 = 3*I
              Q(K3) = Q(K3)+ZASYM
           end do
        end if
!-->    end Bin's changes

!      CALL DVDQ        ! by bin, 2016/10/02
        call ENERGY_1

!-->    modified by Bin, 4/29/2014
!-->    shift z coorindates back after the potential energy
        if (NSURF > 0) then
           J = N
           do I = 1, J
              K3 = 3*I
              Q(K3) = Q(K3)-ZASYM
           end do
        end if
!-->    end Bin's changes

        ESEL = H-EZERO
        exit inner
     end do inner
     exit resample
  end do resample
!
!         CALCULATE THE ROTATIONAL ENERGY
!
  call ROTN(AM, EROT, N)
!
  EINT = ESEL

!-->..Modified by Bin 4/24/2014
!      IF (N.EQ.NATOMS) RETURN
!      IF ((NSURF.NE.0).AND.(N.EQ.NATOMB(1))) RETURN

  if (N == NATOMS .and. NSURF /= 3) return
!-->..end
!
!         RANDOMLY ROTATE THE MOLECULE ABOUT ITS CENTER OF MASS
!         BY EULER'S ANGLES.
!         CENTER OF MASS COORDINATES QQ AND MOMENTA PP ARE PASSED FROM
!         SUBROUTINES CENMAS AND ANGVEL THROUGH COMMON BLOCK WASTE.
!
!-->..Modified by Bin 4/24/2014
!-->..Adapt this to control a specific orientation of initial state
!-->
  if (NTHTA >= 0) then
     if (NTHTA > JROTA) then
        write(6,*) 'NTHTA CAN NOT EXCEED JROTA'
        stop
     end if
     call ROTATEJKM(NTHTA, JROTA, AM, N)
  else
     call ROTATE(N)
  end if
!
!         RECALCULATE THE ROTATIONAL ENERGY AND ANGULAR MOMENTUM
!
  call ROTN(AM, EROT, N)
!
!         CALCULATE THE ROTATIONAL ENERGY
!
  DUM1 = AM(1)/C7
  DUM2 = AM(2)/C7
  DUM3 = AM(3)/C7
!
  write(6,25) EROT, DUM1, DUM2, DUM3
  write(6,15) ESEL
  write(6,*)
end subroutine INITQP
