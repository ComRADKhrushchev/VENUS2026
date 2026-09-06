SUBROUTINE SELECT
  use venus_params
  use venus_data, only: IA_H => IA, KB_H => KB, IB_H => IB, &
                       L => LL, A => A_MPATH
  use venus_data
  implicit double precision (a-h,o-z)
  parameter (NSTEP=5000)

  ! SELECT INITIAL CONDITIONS FOR COORDINATES AND MOMENTA

  dimension QMAXA(NDA3), QMINA(NDA3), PMAXA(NDA), QMAXB(NDA3), &
            QMINB(NDA3), PMAXB(NDA), ANNA(NDA3), ANNB(NDA3)
  dimension ENLOW(NDA3), ENMOD(NDA3), SUM(0:NSTEP)
  dimension QCM(3), VCM(3)
  dimension SAVE_Q(NDA3)
  logical MOVE_A
  save NMBAR, NMA, NMB
  save WWASTORE, WWBSTORE
  ! F10 fix (ported from isolated copy 7eb55ae): ENJA/ENJB were host-implicit
  ! locals in the original monolithic SELECT. After the split into internal
  ! subroutines, names never used in the host body and not in venus_data
  ! became per-subroutine implicit locals, so the EBK energy computed in
  ! select_diatom_a/select_diatom_b no longer reached finalize_and_report
  ! (printed as uninitialized NaN). Declaring them in the host restores host
  ! association for every internal subroutine. (ETAI/ERAI/ETBI/ERBI already
  ! live in venus_data.)
  double precision ENJA, ENJB
  ! F12 defense (ported from isolated copy 7eb55ae): EBK turning points /
  ! test momentum / angular momentum computed on trajectory 1 must survive
  ! into trajectories 2+ (the NSFLAG=1 fast path skips INITEBK). These were
  ! per-subroutine implicit locals whose static residue was not reliable
  ! across calls — promote to host scope next to ENJA/ENJB.
  double precision RMINA, RMAXA, RMASSA, PTESTA, ALA
  double precision RMINB, RMAXB, RMASSB, PTESTB, ALB
  save ENJA, ENJB
  save RMINA, RMAXA, RMASSA, PTESTA, ALA
  save RMINB, RMAXB, RMASSB, PTESTB, ALB

 45 FORMAT(//5X,'SELECT:NORMAL MODE QUANTUM NUMBERS')
 46 FORMAT(5X,10F10.2)
100 FORMAT(5X,'DIATOM A FREQUENCY =',F7.1,' CM-1, AND ENERGY =', &
           F7.2,' KCAL/MOL'/)
106 FORMAT(5X,'SELECT:   JXA,JYA,JZA=',1P3D13.5,' H-BAR'/)
117 FORMAT(5X,'REACTANT A')
124 FORMAT(/5X,'SELECT:    EROTA =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
128 FORMAT(5X,'DIATOM B FREQUENCY =',F7.1,' CM-1, AND ENERGY =', &
           F7.2,' KCAL/MOL'/)
136 FORMAT(5X,'SELECT:   JXB,JYB,JZB=',1P3D13.5,' H-BAR'/)
158 FORMAT(//5X,'REACTANT B')
162 FORMAT(/15X,'IMPACT PARAMETER=',F7.3,' A')
163 FORMAT(/5X,'SELECT:   EROTB =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
196 FORMAT(/5X,'RELATIVE TRANSLATIONAL ENERGY SELECTED: ',F7.2, &
           ' KCAL/MOL'/)
198 FORMAT(/5X,'CHOSEN:   LX,LY,LZ =',1P3D13.5,' H-BAR')
200 FORMAT(/5X,'CHOSEN:   EROT =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
206 FORMAT(/5X,'CHOSEN:   EROTA =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
207 FORMAT(/5X,'CHOSEN:   EROTB =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
208 FORMAT(/5X,'CHOSEN:   EVIBA =',F7.3,' KCAL/MOL'/)
209 FORMAT(/5X,'CHOSEN:   EVIBB =',F7.3,' KCAL/MOL'/)
201 FORMAT('THE ZERO POINT ENERGY CANNOT BE LARGER THAN', &
           ' THE AVAILABLE ENERGY')
202 FORMAT('ZERO POINT ENERGY ',F7.2)
203 FORMAT('AVAILABLE ENERGY ',F7.2)
205 FORMAT('INCREASE THE VALUE OF PARAMETER NSTEP TO AT LEAST ',I8)

  !============================================================
  !  Top-level gate: only NSELT=2 (TRAJECTORY) or NSELT=3
  !  (BARRIER) enter the main initialization body.
  !============================================================
  if (NSELT /= 2 .and. NSELT /= 3) then
     read(5,*) (Q(I), I=1, I3N)
     if (NSELT == 0) then
        read(5,*) (P(I), I=1, I3N)
     else
        do I = 1, I3N
           P(I) = 0.0d0
        end do
     end if
     call DVDQ_1
     call ENERGY_1
     return
  end if

  !============================================================
  !  Set up fragment positions and compute energy reference DH
  !============================================================
  call setup_fragment_positions

  ! Energy reference for separated reactants.
  ! Move the smaller fragment far away so the A-B interaction is zero.
  ! Uses LA/LB arrays for a general fragment separation.
  N_MOVE = NATOMA(1)
  MOVE_A = .true.
  if (NATOMB(1) < N_MOVE) then
     N_MOVE = NATOMB(1)
     MOVE_A = .false.
  end if

  do I = 1, N_MOVE
     if (MOVE_A) then
        IA = LA(1,I)
     else
        IA = LB(1,I)
     end if
     J3 = 3*IA
     K3 = 3*I
     SAVE_Q(K3-2) = Q(J3-2)
     SAVE_Q(K3-1) = Q(J3-1)
     SAVE_Q(K3)   = Q(J3)
     Q(J3-2) = 100.0d0
     Q(J3-1) = 100.0d0
     Q(J3)   = 100.0d0
  end do

  call ENERGY_1
  DH = V

  if (.not. MOVE_A .and. is_diatom_a()) then
     ! Fragment B was moved and diatom A EBK follows.
     ! Pass saved coords to select_diatom_a (restores at line 308).
     Q7_SAVE = SAVE_Q(1)
     Q8_SAVE = SAVE_Q(2)
     Q9_SAVE = SAVE_Q(3)
  else
     do I = 1, N_MOVE
        if (MOVE_A) then
           IA = LA(1,I)
        else
           IA = LB(1,I)
        end if
        J3 = 3*IA
        K3 = 3*I
        Q(J3-2) = SAVE_Q(K3-2)
        Q(J3-1) = SAVE_Q(K3-1)
        Q(J3)   = SAVE_Q(K3)
     end do
  end if

  !============================================================
  !  Fragment A initial conditions
  !============================================================
  if (is_diatom_a()) then
     call select_diatom_a(Q7_SAVE, Q8_SAVE, Q9_SAVE, DH)
  else
     call select_polyatomic_a(DH)
  end if

  !============================================================
  !  Fragment B initial conditions
  !============================================================
  if (NTZ > 1) then
     do kk = 1, 9
        if (P(kk) /= P(kk)) write(6,*) 'D_D:P126 k=', kk
     end do
  end if

  if (NSURF /= 2) then

     if (NATOMB(1) == 0) then
        if (NATOMA(1) /= 2) then
           continue  ! fall through to epilogue
        else
           call DVDQ_1
           call ENERGY_1
        end if
     else if (NATOMB(1) == 1) then
        ! B is a single atom — zero its coordinates
        J3 = 3*LB(1,1)
        J2 = J3 - 1
        J1 = J2 - 1
        Q(J1) = 0.0d0
        Q(J2) = 0.0d0
        Q(J3) = 0.0d0
     else if (is_diatom_b()) then
        call select_diatom_b(DH)
        call restore_a_positions
     else
        call select_polyatomic_b(DH)
     end if

  end if

  !============================================================
  !  Collision setup: gas-phase or surface
  !============================================================
  if (NTZ > 1) then
     do kk = 1, NI
        if (Q(kk) /= Q(kk)) write(6,*) 'D_E:N_Q160 k=', kk
        if (P(kk) /= P(kk)) write(6,*) 'D_E:N_P160 k=', kk
     end do
  end if

  if (NATOMB(1) == 0 .and. NSURF == 0) then
     ! F11 fix (ported from isolated copy 7eb55ae): old code GOTO 999 — a
     ! single-fragment gas-phase run has no bimolecular collision setup
     ! (SEREL/impact-parameter assembly divides by WTB=0, producing 0/0 =
     ! NaN in every momentum).
  else if (NSURF /= 0) then
     call setup_surface_collision
  else
     call setup_gas_collision
  end if

  !============================================================
  !  Epilogue: normal-mode analysis for products, diagnostics
  !============================================================
  call finalize_and_report

  NSFLAG = 1
  return

contains
  !==================================================================
  !  INTERNAL SUBROUTINES
  !==================================================================

  !------------------------------------------------------------
  logical function is_diatom_a()
    implicit double precision (a-h,o-z)
    is_diatom_a = (NATOMA(1) == 2 .and. NLINA /= 0)
  end function

  !------------------------------------------------------------
  logical function is_diatom_b()
    implicit double precision (a-h,o-z)
    is_diatom_b = (NATOMB(1) == 2 .and. NLINB /= 0)
  end function

  !------------------------------------------------------------
  subroutine setup_fragment_positions
    implicit double precision (a-h,o-z)
117 FORMAT(5X,'REACTANT A')

    ! If A is a single atom, zero its coordinates and momenta
    if (NATOMA(1) <= 1 .and. NATOMB(1) <= 2) then
       do K = 1, 3
          J = 3*LA(1,1) - 3 + K
          Q(K) = 0.0d0
          P(K) = 0.0d0
          QQ(K) = 0.0d0
          PP(K) = 0.0d0
       end do
       return
    end if

    if (NATOMA(1) <= 1 .and. NATOMB(1) > 2) then
       ! Full slab: place single adatom at QZA position
       do K = 1, 3
          Q(K) = QZA(1, K)
          P(K) = 0.0d0
          QQ(K) = QZA(1, K)
          PP(K) = 0.0d0
       end do
       return
    end if

    if (NSURF > 0) then
       ! Surface model: place fragment A at equilibrium z + ZASYM
       J = NATOMA(1)
       do I = 1, J
          K3 = 3*I
          Q(K3) = QZA(1,K3) + ZASYM
       end do
    else
       ! Gas phase: displace fragment B by ZASYM
       write(6,117)
       N = NATOMB(1)
       if (N /= 0) then
          do I = 1, N
             J3 = 3*LB(1,I)
             K3 = 3*I
             Q(J3) = QZB(1,K3) + ZASYM
          end do
       end if
    end if
  end subroutine setup_fragment_positions

  !------------------------------------------------------------
  subroutine select_diatom_a(Q7_SAVE, Q8_SAVE, Q9_SAVE, DH)
    implicit double precision (a-h,o-z)
    real*8, intent(in) :: Q7_SAVE, Q8_SAVE, Q9_SAVE, DH
100 FORMAT(5X,'DIATOM A FREQUENCY =',F7.1,' CM-1, AND ENERGY =', &
              F7.2,' KCAL/MOL'/)
106 FORMAT(5X,'SELECT:   JXA,JYA,JZA=',1P3D13.5,' H-BAR'/)

    LBOND(1) = LA(1,1)
    LBOND(2) = LA(1,2)

    ! Normal-mode analysis if requested
    if (NTZ == 1) then
       N = NATOMA(1)
       write(26,*) 'NORMAL MODES FOR FRAGMENT A IN PATH ', 1
       write(26,*)
       call NMODE(N, 0)
       DUM = EIG(6)
       ! F12 defense (ported from isolated copy 7eb55ae): store the diatom
       ! frequency for ALL paths (the old code stored WWASTORE only inside
       ! the temperature branch, so the fixed-n,J path had no NTZ>1
       ! recovery for DUM).
       WWASTORE = EIG(6)*C6
    end if

    if (TRVA >= 0.0d0) then

       if (NACTA == 0) then
          ! Maxwell-Boltzmann sampling for surface oscillator
          call maxwell_boltzmann_a
          return
       end if

       if (NSFLAG /= 1) then
          write(6,*) 'N AND J CHOSEN BASED ON TEMPERATURE'
          write(6,*) 'TRVA = ', TRVA
          ! Calculate vibrational quantum number N
          if (NTZ == 1) then
             WWA(1) = EIG(6)*C6
             WWASTORE = WWA(1)
          else
             WWA(1) = WWASTORE
             DUM = WWA(1)/C6
          end if
          NMBAR = 1
          call THRMAN(WWA, ANNA, TRVA, NMBAR)
          NNA = nint(ANNA(1))
          write(6,*) 'NNA = ', NNA
          ! Calculate rotational quantum number J
          WD1 = W(LBOND(1))
          WD2 = W(LBOND(2))
          call PROBJ(TROTA, AIA, ISEED, JA)
          write(6,*) 'JA = ', JA
       end if
    else
       write(6,*) 'N AND J USED AS INPUT'
    end if

    if (NSFLAG /= 1) then
       ! F12 defense (ported from isolated copy 7eb55ae): on trajectories
       ! after the first, DUM (=EIG(6)) was only set inside the NTZ==1
       ! block; recover it from WWASTORE like the temperature branch does,
       ! instead of using a stale local.
       if (NTZ > 1) DUM = WWASTORE/C6
       ENJA = (dble(NNA) + 0.5d0)*DUM*CM2CAL*C1
       call INITEBK(NNA, JA, RMINA, RMAXA, DH, RMASSA, ENJA, PTESTA, ALA)
       SDUM = ENJA/C1
       write(6,100) DUM, SDUM
    end if

    ! Restore fragment B coordinates after EBK quantization
    Q(7) = Q7_SAVE
    Q(8) = Q8_SAVE
    Q(9) = Q9_SAVE

    ! Rejection-sampling loop for diatom A bond distance
    call diatom_rejection_loop(RMINA, RMAXA, RMASSA, ALA, ENJA, &
                               PTESTA, DH, LBOND(1), LBOND(2), R, PR, SUMM, .true.)

    if (NTZ > 1) then
       write(6,*) 'D_A:HQP R=', R, ' PR=', PR, ' ALA=', ALA, &
                  ' RM=', RMASSA, ' AI=', AI
       write(6,*) 'D_A:ENJA=', ENJA, ' DUM=', DUM, ' RMIN=', RMINA
       write(6,*) 'D_A:RMAX=', RMAXA, ' SUMM=', SUMM
       do kk = 1, 9
          if (Q(kk) /= Q(kk)) write(6,*) 'D_A:N_QbH k=', kk
          if (P(kk) /= P(kk)) write(6,*) 'D_A:N_PbH k=', kk
       end do
    end if

    call HOMOQP(R, PR, ALA, AMA, RMASSA, AI)

    if (NTZ > 1) then
       write(6,*) 'D_B:HQPdone AI=', AI
       do kk = 1, 9
          if (P(kk) /= P(kk)) write(6,*) 'D_B:N_PaH k=', kk
       end do
    end if

    AMAI(1) = AMA(1)/C7
    AMAI(2) = AMA(2)/C7
    AMAI(3) = AMA(3)/C7
    AMAI(4) = AMA(4)/C7
    call ROTN(AMA, EROTA, 2)

    if (NTZ > 1) then
       do kk = 1, 9
          if (P(kk) /= P(kk)) write(6,*) 'D_C:N_PaR k=', kk
       end do
    end if

    write(6,106) AMAI(1), AMAI(2), AMAI(3)

    if (NTZ > 1) then
       do kk = 1, 9
          if (P(kk) /= P(kk)) write(6,*) 'D_C2:N_SEL k=', kk
       end do
    end if
  end subroutine select_diatom_a

  !------------------------------------------------------------
  subroutine maxwell_boltzmann_a
    implicit double precision (a-h,o-z)
    DESKET = sqrt(0.00198717d0*TRVA*C1)
    WT = WTA(1)
    do I = 1, NATOMA(1)
       J3 = 3*LA(1,I)
       J2 = J3 - 1
       J1 = J2 - 1
       J = LA(1,I)
       P(J1) = GASDEV*DESKET*sqrt(W(J))
       P(J2) = GASDEV*DESKET*sqrt(W(J))
       P(J3) = GASDEV*DESKET*sqrt(W(J))
    end do
    write(6,*) 'DBG-MBA TRVA=',TRVA,' DESKET=',DESKET,' P=',P(1:3)

    do I = 1, NATOMA(1)
       L(I) = LA(1,I)
    end do
    call CENMAS(WT, QCM, VCM, NATOMA(1))

    do I = 1, NATOMA(1)
       J = 3*LA(1,I)
       do K = 0, 2
          Q(J-K) = QQ(J-K)
          PP(J-K) = P(J-K)
       end do
    end do

    if (NZDOWN /= 1) call ROTATE(NATOMA(1))

    ! Remove velocities along z (replaced by translational energy later)
    do I = 1, NATOMA(1)
       J = 3*LA(1,I)
       P(J) = 0.0d0
       PP(J) = 0.0d0
    end do
  end subroutine maxwell_boltzmann_a

  !------------------------------------------------------------
  subroutine diatom_rejection_loop(RMIN, RMAX, RMASS, AL, ENJ, &
                                    PTEST, DH, LB1, LB2, R_out, PR_out, SUMM_out, &
                                    shift_z)
    implicit double precision (a-h,o-z)
    real*8, intent(in)  :: RMIN, RMAX, RMASS, AL, ENJ, PTEST, DH
    integer, intent(in) :: LB1, LB2
    real*8, intent(out) :: R_out, PR_out, SUMM_out
    logical, intent(in) :: shift_z
    real*8 :: DUM2, VDUM

    DUM2 = AL**2/(2.0d0*RMASS)

    do
       RAND = RAND0(ISEED)
       R_out = RMIN + (RMAX - RMIN)*RAND
       Q(3*LB1)     = -0.5d0*R_out
       Q(3*LB1 - 1) = 0.0d0
       Q(3*LB1 - 2) = 0.0d0
       Q(3*LB2)     =  0.5d0*R_out
       Q(3*LB2 - 1) = 0.0d0
       Q(3*LB2 - 2) = 0.0d0

       ! Shift z-coordinates up for surface potential evaluation (A only)
       if (shift_z .and. NSURF > 0) then
          J = NATOMA(1)
          do I = 1, J
             K3 = 3*I
             Q(K3) = Q(K3) + ZASYM
          end do
       else
          call DVDQ_1
       end if

       call ENERGY_1

       ! Shift z-coordinates back
       if (shift_z .and. NSURF > 0) then
          J = NATOMA(1)
          do I = 1, J
             K3 = 3*I
             Q(K3) = Q(K3) - ZASYM
          end do
       end if

       VDUM = (V - DH)*C1
       SUMM = ENJ - DUM2/R_out**2 - VDUM
       if (SUMM <= 0.0d0) then
          SUMM = 0.0d0
          PR_out = 0.0d0
       else
          PR_out = sqrt(2.0d0*RMASS*SUMM)
       end if

       RAND = RAND0(ISEED)
       if (PR_out > 0.0d0) then
          if (PTEST/PR_out >= RAND) exit
       else
          exit
       end if
    end do

    SUMM_out = SUMM
    RAND = RAND0(ISEED)
    if (RAND < 0.5d0) PR_out = -PR_out
  end subroutine diatom_rejection_loop

  !------------------------------------------------------------
  subroutine select_polyatomic_a(DH)
    implicit double precision (a-h,o-z)
    real*8, intent(in) :: DH
    real*8 :: DUMNAC6  ! saved across calls via host
201 FORMAT('THE ZERO POINT ENERGY CANNOT BE LARGER THAN', &
           ' THE AVAILABLE ENERGY')
202 FORMAT('ZERO POINT ENERGY ',F7.2)
203 FORMAT('AVAILABLE ENERGY ',F7.2)
205 FORMAT('INCREASE THE VALUE OF PARAMETER NSTEP TO AT LEAST ',I8)

    if (NACTA == 0) then
       call maxwell_boltzmann_poly_a
       return
    end if

    if (NACTA == 1) then
       call rot_and_init_a(ENMTA, DUMNAC6)
       return
    end if

    ! If NSFLAG=1, skip normal-mode analysis for non-quantum cases
    if (NSFLAG == 1) then
       if (NACTA == 5 .or. NACTA == 8 .or. NACTA == 9) then
          call normal_mode_quantum_a
       end if
       call rot_and_init_a(ENMTA, DUMNAC6)
       return
    end if

    ! Full normal-mode setup for polyatomic A
    N = NATOMA(1)
    K = 3*N
    M = 6 - NLINA
    NMA = K - M

    if (NSELT == 3) then
       NMBAR = NMA - 1
       IBARR = 1
    else
       NMBAR = NMA
       IBARR = 0
    end if

    write(26,*) 'NORMAL MODES FOR FRAGMENT A IN PATH ', 1
    write(26,*)
    call NMODE(N, 0)

    do I = 1, NMBAR
       WWA(I) = EIG(I+M+IBARR)*C6
       do J = 1, K
          CA(J,I) = A(I+M+IBARR, J)
       end do
    end do
    if (NSELT == 3) then
       WWA(NMA) = EIG(1)*C6
       do J = 1, K
          CA(J,NMA) = A(1, J)
       end do
    end if

    ! NACTA=6: Beyer-Swinehardt microcanonical setup
    if (NACTA == 6) then
       NBAR = 1
       DUM = ENMTA*CAL2CM
       ZPE = 0.0d0
       do I = 1, NMBAR
          ENMOD(I) = WWA(I)/C6
          ZPE = ZPE + WWA(I)/C6
          ENLOW(I) = ENMOD(I)*0.5d0*CM2CAL*C1
       end do
       ZPE = ZPE*0.5d0
       DUM = DUM - ZPE
       STEP = ENMOD(1)*0.1d0
       if (DUM/STEP >= dble(NSTEP)) then
          write(6,205) idnint(DUM/STEP) + 1
          stop
       end if
       if (ZPE/CAL2CM > ENMTA) then
          write(6,201)
          write(6,202) ZPE*CM2CAL
          write(6,203) ENMTA
          stop
       end if
       DUMNAC6 = DUM
    end if

    if (NACTA == 2 .or. NACTA == 6) then
       call rot_and_init_a(ENMTA, DUMNAC6)
       return
    end if

    ! Normal-mode quantum number selection
    call normal_mode_quantum_a

    call rot_and_init_a(ENMTA, DUMNAC6)
  end subroutine select_polyatomic_a

  !------------------------------------------------------------
  subroutine maxwell_boltzmann_poly_a
    implicit double precision (a-h,o-z)
    DESKET = sqrt(0.00198717d0*TVIBA*C1)
    WT = WTA(1)
    do I = 1, NATOMA(1)
       L(I) = LA(1,I)
       J3 = 3*LA(1,I)
       J2 = J3 - 1
       J1 = J2 - 1
       J = LA(1,I)
       P(J1) = GASDEV*DESKET*sqrt(W(J))
       P(J2) = GASDEV*DESKET*sqrt(W(J))
       P(J3) = GASDEV*DESKET*sqrt(W(J))
    end do

    call CENMAS(WT, QCM, VCM, NATOMA(1))

    do I = 1, NATOMA(1)
       J = 3*LA(1,I)
       do K = 0, 2
          Q(J-K) = QQ(J-K)
          PP(J-K) = P(J-K)
       end do
    end do

    call ROTATE(NATOMA(1))

    do I = 1, NATOMA(1)
       J = 3*LA(1,I)
       P(J) = 0.0d0
       PP(J) = 0.0d0
    end do
  end subroutine maxwell_boltzmann_poly_a

  !------------------------------------------------------------
  subroutine normal_mode_quantum_a
    implicit double precision (a-h,o-z)
 45 FORMAT(//5X,'SELECT:NORMAL MODE QUANTUM NUMBERS')
 46 FORMAT(5X,10F10.2)

    if (NACTA == 5) then
       call THRMAN(WWA, ANQA, TVIBA, NMBAR)
       write(6,45)
       write(6,46) (ANQA(I), I=1, NMBAR)
       write(9,46) (ANQA(I), I=1, NMBAR)
    end if

    if (NACTA == 9) call MICROCI(ECONE, EQNM)

    if (NACTA == 8 .or. NACTA == 9) then
       if (NACTA == 8) then
          EQNM = ENMTA
          call QMMICRO(WWA, ANQA, EQNM, NMA)
          write(6,45)
          write(6,46) (ANQA(I), I=1, NMA)
          write(9,46) (ANQA(I), I=1, NMA)
       else if (NACTA == 9) then
          call QMMICRO(WWA, ANQA, EQNM, NM3N8)
          write(6,45)
          write(6,46) (ANQA(I), I=1, NM3N8)
          write(9,46) (ANQA(I), I=1, NM3N8)
       end if
    end if

    ! Calculate normal-mode energies and amplitudes
    ENMDUM = 0.0d0
    do I = 1, NMBAR
       DUM = (ANQA(I) + 0.5d0)*WWA(I)/(C6*CAL2CM)
       ENMDUM = ENMDUM + DUM
       DUM = DUM*C1
       AMPA(I) = sqrt(2.0d0*DUM)/WWA(I)
    end do
    if (NACTA == 3 .or. NACTA == 4 .or. NACTA == 5) ENMTA = ENMDUM
    if (NACTA == 8 .or. NACTA == 9) ENM89 = ENMDUM
  end subroutine normal_mode_quantum_a

  !------------------------------------------------------------
  subroutine rot_and_init_a(ENMTA_in, DUMNAC6)
    implicit double precision (a-h,o-z)
    real*8, intent(inout) :: ENMTA_in
    real*8, intent(in) :: DUMNAC6
124 FORMAT(/5X,'SELECT:    EROTA =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)

    call ROTEN(AMA, AI, TROTA, EROTA, NROTA, NLINA, JROTA, KROTA)

    AMAI(1) = AMA(1)/C7
    AMAI(2) = AMA(2)/C7
    AMAI(3) = AMA(3)/C7
    AMAI(4) = sqrt(AMAI(1)**2 + AMAI(2)**2 + AMAI(3)**2)
    write(6,124) EROTA, AMAI(1), AMAI(2), AMAI(3)

    ETAI = EROTA + ENMTA_in
    if (NACTA == 8 .or. NACTA == 9) ETAI = EROTA + ENM89

    ! Set up atomic indices and equilibrium coordinates
    N = NATOMA(1)
    do I = 1, N
       L(I) = LA(1,I)
       J3 = 3*L(I)
       J2 = J3 - 1
       J1 = J2 - 1
       K3 = 3*I
       K2 = K3 - 1
       K1 = K2 - 1
       QZ(J1) = QZA(1,K1)
       QZ(J2) = QZA(1,K2)
       QZ(J3) = QZA(1,K3)
    end do
    DUM1 = WTA(1)

    ! Orthant sampling
    if (NACTA == 1) then
       call ORTHAN(AMA, DUM1, ENMTA_in, ETAI, QMAXA, QMINA, PMAXA, &
                   PSCALA, ERAI, N)
       return
    end if

    ! Classical microcanonical normal-mode energies
    if (NACTA == 2) then
       DUM = ENMTA_in*C1
       NN = NMBAR - 1
       do I = 1, NN
          RAND = RAND0(ISEED)
          SDUM = 1.0d0/dble(NMBAR - I)
          SDUM = DUM*(1.0d0 - RAND**SDUM)
          DUM = DUM - SDUM
          AMPA(I) = sqrt(2.0d0*SDUM)/WWA(I)
       end do
       AMPA(NMBAR) = sqrt(2.0d0*DUM)/WWA(NMBAR)
    end if

    ! Beyer-Swinehardt microcanonical sampling (NACTA=6)
    if (NACTA == 6) then
       DUM = DUMNAC6
       do II = 1, NMBAR
          IMAXST = idint(DUM/ENMOD(II))
          if (II /= NMBAR) then
             STEP = ENMOD(II+1)*0.1d0
             IEND = idnint(DUM/STEP)
             do J = 0, IEND
                SUM(J) = 1.0d0
             end do
             do J = II+1, NMBAR
                ISTART = idnint(ENMOD(J)/STEP)
                do K = ISTART, IEND
                   SUM(K) = SUM(K) + SUM(K-ISTART)
                end do
             end do
             PROB = 0.0d0
             do J = IMAXST, 0, -1
                ENDUM = DUM - dble(J)*ENMOD(II)
                IDDUM = idnint(ENDUM/STEP)
                PROB = PROB + SUM(IDDUM)
             end do
          else
             STEP = ENMOD(II)
             PROB = 0.0d0
             do J = IMAXST, 0, -1
                ENDUM = DUM - dble(J)*ENMOD(II)
                IDDUM = idnint(ENDUM/STEP)
                SUM(IDDUM) = 1.0d0
                PROB = PROB + SUM(IDDUM)
             end do
          end if
          ! Find the selected quantum number by rejection
          RAND = RAND0(ISEED)*PROB
          if (RAND > PROB) RAND = PROB
          PROB = 0.0d0
          do J = IMAXST, 0, -1
             ENDUM = DUM - dble(J)*ENMOD(II)
             IDDUM = idnint(ENDUM/STEP)
             PROB = SUM(IDDUM) + PROB
             if (RAND <= PROB) then
                SDUM = ENMOD(II)*dble(J)
                DUM = DUM - SDUM
                SDUM = SDUM*CM2CAL*C1
                AMPA(II) = sqrt(2.0d0*(SDUM + ENLOW(II)))/WWA(II)
                exit
             end if
          end do
       end do
       EBAR = DUM*CM2CAL
       ETAI = ETAI - EBAR
    end if

    call INITQP(WWA, AMPA, CA, AMA, DUM1, ETAI, EROTA, AI, ERAI, &
                PHASEA, N, NMBAR, NPHASA, IPHASA)

    if (NSELT == 3) call BAREXC(DUM1, CA, AMA, ERAI, N, NMA)
  end subroutine rot_and_init_a

  !==================================================================
  !  Fragment B internal subroutines
  !==================================================================

  !------------------------------------------------------------
  subroutine select_diatom_b(DH)
    implicit double precision (a-h,o-z)
    real*8, intent(in) :: DH
128 FORMAT(5X,'DIATOM B FREQUENCY =',F7.1,' CM-1, AND ENERGY =', &
              F7.2,' KCAL/MOL'/)
136 FORMAT(5X,'SELECT:   JXB,JYB,JZB=',1P3D13.5,' H-BAR'/)
158 FORMAT(//5X,'REACTANT B')

    write(6,158)
    ! Save fragment A coords
    N = NATOMA(1)
    do I = 1, N
       J3 = 3*LA(1,I)
       J2 = J3 - 1
       J1 = J2 - 1
       K3 = 3*I
       K2 = K3 - 1
       K1 = K2 - 1
       Q(J1) = QZA(1,K1)
       Q(J2) = QZA(1,K2)
       Q(J3) = QZA(1,K3) + ZASYM
       P(J1) = 0.0d0
       P(J2) = 0.0d0
       P(J3) = 0.0d0
    end do

    ! Reset Q array to QZ for B
    N = NATOMB(1)
    do I = 1, N
       J3 = 3*LB(1,I)
       K3 = 3*I
       Q(J3) = QZB(1,K3)
    end do

    LBOND(1) = LB(1,1)
    LBOND(2) = LB(1,2)

    if (NSFLAG /= 1) then
       if (NTZ == 1) then
          I = NATOMA(1)
          N = NATOMB(1)
          write(26,*) 'NORMAL MODES FOR FRAGMENT B IN PATH ', 1
          write(26,*)
          call NMODE(N, I)
          DUM = EIG(6)
       end if

       if (TRVB >= 0.0d0) then
          write(6,*) 'N AND J CHOSEN BASED ON TEMPERATURE'
          write(6,*) 'TRVB = ', TRVB
          if (NTZ == 1) then
             WWB(1) = EIG(6)*C6
             WWBSTORE = WWB(1)
          else
             WWB(1) = WWBSTORE
             DUM = WWB(1)/C6
          end if
          NMBAR = 1
          call THRMAN(WWB, ANNB, TRVB, NMBAR)
          NNB = nint(ANNB(1))
          write(6,*) 'NNB = ', NNB
          WD1 = W(LBOND(1))
          WD2 = W(LBOND(2))
          call PROBJ(TROTB, AIB, ISEED, JB)
          write(6,*) 'JB = ', JB
       else
          write(6,*) 'N AND J USED AS INPUT'
       end if

       ENJB = (dble(NNB) + 0.5d0)*DUM*CM2CAL*C1
       call INITEBK(NNB, JB, RMINB, RMAXB, DH, RMASSB, ENJB, PTESTB, ALB)
       SDUM = ENJB/C1
       write(6,128) DUM, SDUM
    end if

    ! Rejection-sampling loop for diatom B bond distance
    call diatom_rejection_loop(RMINB, RMAXB, RMASSB, ALB, ENJB, &
                               PTESTB, DH, LBOND(1), LBOND(2), R, PR, SUMM, .false.)

    call HOMOQP(R, PR, ALB, AMB, RMASSB, BI)
    AMBI(1) = AMB(1)/C7
    AMBI(2) = AMB(2)/C7
    AMBI(3) = AMB(3)/C7
    AMBI(4) = AMB(4)/C7
    call ROTN(AMB, EROTB, 2)
    write(6,136) AMBI(1), AMBI(2), AMBI(3)
  end subroutine select_diatom_b

  !------------------------------------------------------------
  subroutine select_polyatomic_b(DH)
    implicit double precision (a-h,o-z)
    real*8, intent(in) :: DH
 45 FORMAT(//5X,'SELECT:NORMAL MODE QUANTUM NUMBERS')
 46 FORMAT(5X,10F10.2)
158 FORMAT(//5X,'REACTANT B')

    write(6,158)
    ! Save fragment A coords and reset to QZ
    N = NATOMA(1)
    do I = 1, N
       J3 = 3*LA(1,I)
       J2 = J3 - 1
       J1 = J2 - 1
       K3 = 3*I
       K2 = K3 - 1
       K1 = K2 - 1
       Q(J1) = QZA(1,K1)
       Q(J2) = QZA(1,K2)
       if (NATOMB(1) <= 2) then
          Q(J3) = QZA(1,K3) + ZASYM
       else
          Q(J3) = QZA(1,K3)
       end if
       P(J1) = 0.0d0
       P(J2) = 0.0d0
       P(J3) = 0.0d0
       QQ(J1) = Q(J1)
       QQ(J2) = Q(J2)
       QQ(J3) = Q(J3)
       PP(J1) = 0.0d0
       PP(J2) = 0.0d0
       PP(J3) = 0.0d0
    end do

    ! Reset Q for B
    N = NATOMB(1)
    do I = 1, N
       J3 = 3*LB(1,I)
       K3 = 3*I
       Q(J3) = QZB(1,K3)
    end do

    ! NACTB=7: MD equilibration — skip B selection, go to MD
    if (NATOMB(1) > 2 .and. NACTB == 7) then
       call setup_b_coords
       call md_equilibrate_b
       call restore_a_positions
       return
    end if

    if (NACTB == 1) then
       call rot_and_init_b(ENMTB)
       return
    end if

    if (NACTB == 7) then
       call setup_b_coords
       call md_equilibrate_b
       call restore_a_positions
       return
    end if

    ! Normal-mode analysis for B
    if (NSFLAG == 1) then
       if (NACTB /= 5 .and. NACTB /= 8) then
          call rot_and_init_b(ENMTB)
          return
       end if
    end if

    N = NATOMB(1)
    K = 3*N
    M = 6 - NLINB
    NMB = K - M
    I = NATOMA(1)
    write(26,*) 'NORMAL MODES FOR FRAGMENT B IN PATH ', 1
    write(26,*)
    call NMODE(N, I)

    do I = 1, NMB
       WWB(I) = EIG(I+M)*C6
       do J = 1, K
          CB(J,I) = A(I+M, J)
       end do
    end do

    if (NACTB == 2) then
       call rot_and_init_b(ENMTB)
       return
    end if

    ! Quantum number selection
    if (NACTB == 5) then
       call THRMAN(WWB, ANQB, TVIBB, NMB)
       write(6,45)
       write(6,46) (ANQB(I), I=1, NMB)
       write(9,46) (ANQB(I), I=1, NMB)
    end if

    if (NACTB == 8) then
       EQNM = ENMTB
       call QMMICRO(WWB, ANQB, EQNM, NMB)
       write(6,45)
       write(6,46) (ANQB(I), I=1, NMB)
       write(9,46) (ANQB(I), I=1, NMB)
    end if

    ! Calculate normal-mode energies and amplitudes
    ENMDUM = 0.0d0
    do I = 1, NMB
       DUM = (ANQB(I) + 0.5d0)*WWB(I)/(C6*CAL2CM)
       ENMDUM = ENMDUM + DUM
       DUM = DUM*C1
       AMPB(I) = sqrt(2.0d0*DUM)/WWB(I)
    end do
    if (NACTB == 3 .or. NACTB == 4 .or. NACTB == 5) ENMTB = ENMDUM
    if (NACTB == 8) ENM89 = ENMDUM

    call rot_and_init_b(ENMTB)
  end subroutine select_polyatomic_b

  !------------------------------------------------------------
  subroutine setup_b_coords
    implicit double precision (a-h,o-z)
    N = NATOMB(1)
    do I = 1, N
       L(I) = LB(1,I)
       J3 = 3*L(I)
       J2 = J3 - 1
       J1 = J2 - 1
       K3 = 3*I
       K2 = K3 - 1
       K1 = K2 - 1
       QZ(J1) = QZB(1,K1)
       QZ(J2) = QZB(1,K2)
       QZ(J3) = QZB(1,K3)
       Q(J1)  = QZB(1,K1)
       Q(J2)  = QZB(1,K2)
       Q(J3)  = QZB(1,K3)
    end do
  end subroutine setup_b_coords

  !------------------------------------------------------------
  subroutine rot_and_init_b(ENMTB_in)
    implicit double precision (a-h,o-z)
    real*8, intent(inout) :: ENMTB_in
163 FORMAT(/5X,'SELECT:   EROTB =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)

    call ROTEN(AMB, BI, TROTB, EROTB, NROTB, NLINB, JROTB, KROTB)

    AMBI(1) = AMB(1)/C7
    AMBI(2) = AMB(2)/C7
    AMBI(3) = AMB(3)/C7
    AMBI(4) = sqrt(AMBI(1)**2 + AMBI(2)**2 + AMBI(3)**2)
    write(6,163) EROTB, AMBI(1), AMBI(2), AMBI(3)

    ETBI = EROTB + ENMTB_in
    if (NACTB == 8) ETBI = EROTB + ENM89

    call setup_b_coords
    if (NACTB == 7) return

    DUM1 = WTB(1)

    if (NACTB == 1) then
       N = NATOMB(1)
       call ORTHAN(AMB, DUM1, ENMTB_in, ETBI, QMAXB, QMINB, PMAXB, &
                   PSCALB, ERBI, N)
       return
    end if

    if (NACTB == 2) then
       DUM = ENMTB_in*C1
       NN = NMB - 1
       do I = 1, NN
          RAND = RAND0(ISEED)
          SDUM = 1.0d0/dble(NMB - I)
          SDUM = DUM*(1.0d0 - RAND**SDUM)
          DUM = DUM - SDUM
          AMPB(I) = sqrt(2.0d0*SDUM)/WWB(I)
       end do
       AMPB(NMB) = sqrt(2.0d0*DUM)/WWB(NMB)
    end if

    N = NATOMB(1)
    call INITQP(WWB, AMPB, CB, AMB, DUM1, ETBI, EROTB, BI, ERBI, &
                PHASEB, N, NMB, NPHASB, IPHASB)
  end subroutine rot_and_init_b

  !------------------------------------------------------------
  subroutine restore_a_positions
    implicit double precision (a-h,o-z)
    N = NATOMA(1)
    do I = 1, N
       do K = 1, 3
          J = 3*LA(1,I) - 3 + K
          Q(J) = QQ(J)
          P(J) = PP(J)
       end do
    end do
  end subroutine restore_a_positions

  !------------------------------------------------------------
  subroutine md_equilibrate_b
    implicit double precision (a-h,o-z)

    if (NACTB /= 7) return

    TELEC = THERMOTEMP

    ! Store coordinates and momenta of reactant A first
    if (THERMOTEMP <= 0.01d0) then
       N = NATOMB(1) - NRGD
       do I = 1, N
          J3 = 3*LB(1,I)
          J2 = J3 - 1
          J1 = J2 - 1
          J = LB(1,I)
          P(J1) = 0.0d0
          P(J2) = 0.0d0
          P(J3) = 0.0d0
       end do
       return
    end if

    N = NATOMA(1)
    do I = 1, N
       J3 = 3*LA(1,I)
       J2 = J3 - 1
       J1 = J2 - 1
       QTEMP(J3) = Q(J3)
       QTEMP(J2) = Q(J2)
       QTEMP(J1) = Q(J1)
       P(J3) = 0.0d0
       P(J2) = 0.0d0
       P(J1) = 0.0d0
    end do

    ! Shift QTEMP for fragment A far above the slab (z+100 A).
    ! The existing freeze-at-QTEMP logic (lines 1194/1234) will then
    ! keep fragment A at a safe distance during equilibration, avoiding
    ! spurious large repulsive forces on slab atoms.
    ! The final restore subtracts this offset to recover the true position.
    do I = 1, NATOMA(1)
       J3 = 3*LA(1,I)
       QTEMP(J3) = QTEMP(J3) + 100.0d0
       Q(J3)     = QTEMP(J3)              ! also for NTZ>1 DVDQ_1
    end do

    NCOORORG = NCOOR
    NCOOR = 0
    N = NATOMB(1) - NRGD

    TIMEORG = TIME
    ATIMEORG = ATIME
    TIME = 0.05d0
    ATIME = TIME
    if (INTEGRATOR == 1) then
       write(*,*) 'INTEGRATOR=1 IS NOT SUPPORTED IN NACTB=7'
       write(*,*) 'SEE SELECT.f90 FOR DETAILS'
       stop
    end if

    NC = 0

    if (NTZ == 1) then
       DESKET = sqrt(0.00198717d0*THERMOTEMP*C1)
       do I = 1, N
          J3 = 3*LB(1,I)
          J2 = J3 - 1
          J1 = J2 - 1
          J = LB(1,I)
          P(J1) = GASDEV*DESKET*sqrt(W(J))
          P(J2) = GASDEV*DESKET*sqrt(W(J))
          P(J3) = GASDEV*DESKET*sqrt(W(J))
       end do

       ! Ensure total linear momentum is zero
       SUMX = 0.0d0
       SUMY = 0.0d0
       SUMZ = 0.0d0
       do I = 1, N
          J3 = 3*LB(1,I)
          J2 = J3 - 1
          J1 = J2 - 1
          SUMX = SUMX + P(J1)
          SUMY = SUMY + P(J2)
          SUMZ = SUMZ + P(J3)
       end do
       SUMX = SUMX/dble(N)
       SUMY = SUMY/dble(N)
       SUMZ = SUMZ/dble(N)
       do I = 1, N
          J3 = 3*LB(1,I)
          J2 = J3 - 1
          J1 = J2 - 1
          P(J1) = P(J1) - SUMX
          P(J2) = P(J2) - SUMY
          P(J3) = P(J3) - SUMZ
       end do

       if (NSCALE >= 0) call THERMO(0, 1)
       NSE = NSCALE + NEQUAL
    else
       do I = 1, N
          J3 = 3*LB(1,I)
          J2 = J3 - 1
          J1 = J2 - 1
          Q(J3) = QTEMP(J3)
          Q(J2) = QTEMP(J2)
          Q(J1) = QTEMP(J1)
          P(J3) = PTEMP(J3)
          P(J2) = PTEMP(J2)
          P(J1) = PTEMP(J1)
       end do
       call DVDQ_1
       RAND = RAND0(ISEED)
       NSE = int(RAND*100) + 100
    end if

    if (NSE > 0) then
       if (INTEGRATOR == 2) then
          call DVDQ_1
       end if

       ! Equilibration: run dynamics with thermostat
       do I = 1, NSE
          NC = NC + 1
          if (INTEGRATOR == 2) then
             call SYMPLE(LLL, 0)
          else if (INTEGRATOR == 3) then
             call VERLET(LLL, TELEC, 0)
          end if
          do IA = 1, NATOMA(1)
             J3 = 3*LA(1,IA)
             J2 = J3 - 1
             J1 = J2 - 1
             Q(J3) = QTEMP(J3)
             Q(J2) = QTEMP(J2)
             Q(J1) = QTEMP(J1)
             P(J3) = 0.0d0
             P(J2) = 0.0d0
             P(J1) = 0.0d0
          end do
          call ENERGY_1
          if (NSEL == 1) then
             TB = 0.0d0
             do II = 1, N
                J = LB(1,II)
                J3 = 3*J
                J2 = J3 - 1
                J1 = J2 - 1
                TB = TB + (P(J1)**2 + P(J2)**2 + P(J3)**2)/W(J)
             end do
             TEMPINIT = TB/(3.0d0*dble(N)*0.00198717d0*C1)
             write(30,*) 'SYSTEM TEMPERATURE=', TEMPINIT
          end if
          if (NC <= (NSCALE+NEQUAL) .and. mod(NC,100) == 0) then
             call THERMO(NC, 1)
          end if
       end do
    end if

    ! Save coordinates and momenta for the surface
    do I = 1, N
       J3 = 3*LB(1,I)
       J2 = J3 - 1
       J1 = J2 - 1
       PTEMP(J3) = P(J3)
       PTEMP(J2) = P(J2)
       PTEMP(J1) = P(J1)
       QTEMP(J3) = Q(J3)
       QTEMP(J2) = Q(J2)
       QTEMP(J1) = Q(J1)
    end do

    ! Restore reactant A coordinates and momenta
    N = NATOMA(1)
    do I = 1, N
       J3 = 3*LA(1,I)
       J2 = J3 - 1
       J1 = J2 - 1
       P(J3) = PP(J3)
       P(J2) = PP(J2)
       P(J1) = PP(J1)
       Q(J3) = QTEMP(J3) - 100.0d0
       Q(J2) = QTEMP(J2)
       Q(J1) = QTEMP(J1)
    end do

    write(6,*) 'EQUALIBRATION FOR NACTB EQ 7 IS NOW OVER'

    TIME = TIMEORG
    ATIME = ATIMEORG
    NFINAL = 0
    NCOOR = NCOORORG
    NC = 0
  end subroutine md_equilibrate_b

  !==================================================================
  !  Collision setup and finalization
  !==================================================================

  !------------------------------------------------------------
  subroutine setup_surface_collision
    implicit double precision (a-h,o-z)

    call SURF(NSURF)

    if (NSURF == 1 .and. NATOMB(1) <= 2) then
       surf_z = Q(3*LA(1,1))
       do I = 1, NATOMB(1)
          J = 3*LB(1,I)
          Q(J) = Q(J) + surf_z
       end do
       do I = 1, NATOMA(1)
          K3 = 3*I
          Q(K3) = Q(K3) - ZASYM
       end do
       do I = 1, NATOMB(1)
          J = 3*LB(1,I)
          Q(J) = Q(J) - ZASYM
       end do
    end if

    if (NTZ > 1) then
       do kk = 1, NI
          if (Q(kk) /= Q(kk)) write(6,*) 'D_F:N_QS k=', kk
          if (P(kk) /= P(kk)) write(6,*) 'D_F:N_PS k=', kk
       end do
    end if

    if (NSURF == 2 .and. NGLO /= 0) then
       call GLOSELECT(TVIBB)
    end if

    call DVDQ_1
    call ENERGY_1
  end subroutine setup_surface_collision

  !------------------------------------------------------------
  subroutine setup_gas_collision
    implicit double precision (a-h,o-z)
162 FORMAT(/15X,'IMPACT PARAMETER=',F7.3,' A')
198 FORMAT(/5X,'CHOSEN:   LX,LY,LZ =',1P3D13.5,' H-BAR')

    if (NRNDXY == 1) then
       RAND = RAND0(ISEED)
       RX0 = RND_BOX*(2.0d0*RAND - 1.0d0)
       RAND = RAND0(ISEED)
       RY0 = RND_BOX*(2.0d0*RAND - 1.0d0)
       N = NATOMA(1)
       do I = 1, N
          J = 3*LA(1,I)
          Q(J-2) = Q(J-2) + RX0
          Q(J-1) = Q(J-1) + RY0
       end do
    end if

    SB = BMAX
    if (NOB /= 1) then
       RAND = RAND0(ISEED)
       SB = BMAX*sqrt(RAND)
    end if
    write(6,162) SB

    DUM1 = sqrt(S*S - SB*SB)
    N = NATOMB(1)
    do I = 1, N
       J = 3*LB(1,I)
       Q(J)   = Q(J)   + DUM1
       Q(J-1) = Q(J-1) + SB
    end do

    if (NREL == 0) then
       DUM = GAMA(2, ISEED)
       SEREL = 0.00198717d0*DUM*TRANS
       SEREL = SEREL*C1
    end if

    WT = WTA(1) + WTB(1)
    SDUM = WTA(1)*WTB(1)/WT
    DUM = sqrt(2.0d0*SEREL/SDUM)
    VELA = DUM*WTB(1)/WT
    VELB = VELA - DUM
    N = NATOMA(1)
    do I = 1, N
       J = 3*LA(1,I)
       P(J) = P(J) + VELA*W(LA(1,I))
    end do
    N = NATOMB(1)
    do I = 1, N
       J = 3*LB(1,I)
       P(J) = P(J) + VELB*W(LB(1,I))
    end do

    call DVDQ_1
    call ENERGY_1

    VI(1) = 0.0d0
    VI(2) = 0.0d0
    VI(3) = DUM
    VI(4) = DUM
    OAMI(1) = -SB*DUM*SDUM/C7
    OAMI(2) = 0.0d0
    OAMI(3) = 0.0d0
    OAMI(4) = abs(OAMI(1))
    write(6,198) (OAMI(I), I=1,3)
  end subroutine setup_gas_collision

  !------------------------------------------------------------
  subroutine finalize_and_report
    implicit double precision (a-h,o-z)
196 FORMAT(/5X,'RELATIVE TRANSLATIONAL ENERGY SELECTED: ',F7.2, &
           ' KCAL/MOL'/)
206 FORMAT(/5X,'CHOSEN:   EROTA =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
207 FORMAT(/5X,'CHOSEN:   EROTB =',F7.3,' KCAL/MOL'/ &
           15X,'JX,JY,JZ =',1P3D13.5,' H-BAR'/)
208 FORMAT(/5X,'CHOSEN:   EVIBA =',F7.3,' KCAL/MOL'/)
209 FORMAT(/5X,'CHOSEN:   EVIBB =',F7.3,' KCAL/MOL'/)

    ! Normal-mode analysis for each product channel (gas phase only)
    if (NSURF == 0) then
       if (NSFLAG /= 1) then
          if (NTZ == 1) then
             N = NATOMA(1) + NATOMB(1)
             do J = 1, 3*N
                QTEMP(J) = Q(J)
             end do

             do I = 2, NPATHS+1
                N = NATOMA(I)
                do K = 1, 3*N
                   Q(K) = QZA(I,K)
                end do
                M = NATOMB(I)
                do K = 1, 3*M
                   Q(3*N+K) = QZB(I,K)
                end do
                do K = 1, N
                   K3 = 3*K
                   Q(K3) = QZA(I,K3) + ZASYM
                end do

                if (N >= 2) then
                   write(26,*) 'NORMAL MODES FOR FRAGMENT A IN PATH ', I
                   write(26,*)
                   call NMODE(N, 0)
                end if
                if (M >= 2) then
                   write(26,*) 'NORMAL MODES FOR FRAGMENT B IN PATH ', I
                   write(26,*)
                   call NMODE(M, N)
                end if
             end do

             N = NATOMA(1) + NATOMB(1)
             do J = 1, 3*N
                Q(J) = QTEMP(J)
             end do
          end if
       end if
    end if

    ! Write final rotational angular momenta
    if (NATOMA(1) > 1) then
       AMAI(1) = AMA(1)/C7
       AMAI(2) = AMA(2)/C7
       AMAI(3) = AMA(3)/C7
       AMAI(4) = sqrt(AMAI(1)**2 + AMAI(2)**2 + AMAI(3)**2)
       if (NATOMA(1) == 2) then
          write(6,206) EROTA, AMAI(1), AMAI(2), AMAI(3)
          write(6,208) ENJA/C1
       else
          write(6,206) ERAI, AMAI(1), AMAI(2), AMAI(3)
          write(6,208) ETAI - ERAI
       end if
    end if

    if (NATOMB(1) > 1) then
       AMBI(1) = AMB(1)/C7
       AMBI(2) = AMB(2)/C7
       AMBI(3) = AMB(3)/C7
       AMBI(4) = sqrt(AMBI(1)**2 + AMBI(2)**2 + AMBI(3)**2)
       if (NATOMA(1) == 2) then
          write(6,207) EROTB, AMBI(1), AMBI(2), AMBI(3)
          write(6,209) ENJB/C1
       else
          write(6,207) ERBI, AMBI(1), AMBI(2), AMBI(3)
          write(6,209) ETBI - ERBI
       end if
    end if

    write(6,196) SEREL/C1
  end subroutine finalize_and_report

END SUBROUTINE SELECT
