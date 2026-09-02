!
!    VELOCITY VERLET (1) AND BEEMAN'S ALGORITHM (2) INTEGRATOR
!    CODED BY KI SONG, 5/23/05
!    CODED BY BIN JIANG, 2/15/2017
!
!    VELOCITY VERLET
!      X(T+DT)=X(T)+V(T)DT+F(T)/2M*DT**2
!      V(T+DT)=V(T)+(F(T+DT)+FT)/2M*DT
!
!    BEEMAN
!      X(T+DT)=X(T)-V(T)*DT+(4*F(T)-F(T-DT))*DT**2/M/6
!      V(T)=(X(T+DT)-X(T)/DT+(2*F(T+DT)+F(T))*DT/6
!
!    F(T)=-DV/DQ=PDOT
!    VENUS PUT PDOT LIKE THIS. SYMPLE.F USES THIS EXPRESSION.
!
      SUBROUTINE VERLET(NVERLET,TELEC,NSELECT)
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FT(NDA3),QTP1(NDA3),QTM1(NDA3),QTP2(NDA3),QTM2(NDA3)
      DIMENSION PC(NDA3),FT2(NDA3)
      DIMENSION Q_end_save(NDA3),Q_mid(NDA3)
      DIMENSION F0_e0(NDA3),F0_core(NDA3),F0_mfv(NDA3)
      DIMENSION Fmid_e0(NDA3),Fmid_core(NDA3),Fmid_mfv(NDA3)
      SAVE PC,QTM1,QTM2,QTP1,FT2
!
!   VELOCITY VERLET
!
        Q_old(1:NI) = Q(1:NI)
      IF (NC.EQ.1.OR.NVERLET.EQ.1) THEN
!
! ---- Compute force F(t) at start-of-step position (self-contained) ----
          if(CALTYP <0 .OR. NSELECT == 0) then
              call DVDQ_1
          end if


          DO K=1,NI
              P(K)=P(K)+PDOT(K)*ATIME*0.5D0  ! First half-step: v(t+dt/2)
          ENDDO

          ! Ghost anchor Langevin half-step 1 (VV)
          if (FCG > 0.0d0) then
             do i = 1, 3
                k = 3 + i
                F_lang_save(i) = -FCG * P(k) + GASDEV() * GSW_ghost
                P(k) = P(k) + F_lang_save(i) * ATIME * 0.5d0
             end do
          end if

          DO K=1,NI
              KK=(K+2)/3
              FT2(K)=PDOT(K)
              Q(K)=Q(K)+P(K)*ATIME/W(KK)      ! Position update: x(t+dt)
          ENDDO


! ---- Compute force F(t+dt) at end-of-step position ----
          if(CALTYP <0 .OR. NSELECT == 0) then
              call DVDQ_1
          end if
          

          DO K=1,NI
               P(K)=P(K)+PDOT(K)*ATIME*0.5D0  ! Second half-step: v(t+dt)
          ENDDO

          ! Ghost anchor Langevin half-step 2 + work (VV)
          dW_lang = 0.0d0
          if (FCG > 0.0d0) then
             do i = 1, 3
                k = 3 + i
                F_lang2 = -FCG * P(k) + GASDEV() * GSW_ghost
                P(k) = P(k) + F_lang2 * ATIME * 0.5d0
                dW_lang = dW_lang + 0.5d0*(F_lang_save(i)+F_lang2)*(Q(k)-Q_old(k))
             end do
             dW_lang = dW_lang / C1
          end if



!        ENDIF
!
!*********************************************************************
!   BEEMAN'S THIRD ORDER PROCEDUCE, BEEMAN, J. COMPUT. PHYS. 20, 130
!   (1976), MODIFIED BY TULLY ET AL., J. CHEM. PHYS. 71, 1630 (1979)
!*********************************************************************
!
      ELSE 

!   IF NC>1 THEN DO BEEMAN'S THIRD ORDER PROPAGATION, FT2 AND QTP1 ARE
!   SAVED AS THE FORCES AT THE (N-1)TH STEP AND THE COORDINATES AT THE
!   NTH STEP
!        IF (NGLO.NE.0) THEN

!        ENDIF
!...update position
          DO K=1,NI
            KK=(K+2)/3
            QTP1(K)=Q(K)
            Q(K)=Q(K)+P(K)*ATIME/W(KK)+ATIME**2*(4*PDOT(K)-FT2(K))/W(KK)/6
          ENDDO

!...predicted velocity (use FT2=F(t-dt) BEFORE overwriting)
          DO K=1,NI
            KK=(K+2)/3
            P(K)=P(K)+0.5D0*ATIME*(3*PDOT(K)-FT2(K))
          ENDDO

!...save current forces for next step
          DO K=1,NI
            FT2(K)=PDOT(K)
          ENDDO
        

!...UPDATE VELOCITY (MOMENTUM)

         !####################################################################
         !---------------------------------------------------------------------
         ! This is where positions are updated before calculating forces
         ! All force calculations MUST be done after this point
        if(CALTYP <0 .OR. NSELECT == 0) then
            call DVDQ_1
        end if
        !---------------------------------------------------------------------


!.....FOR LDFA MODEL
          IF (NFC.NE.0) THEN
  !           CALL FRICTION(NFC)
   !          CALL FRICFORCE(NFC,TELEC)
          END IF
!.....END 

!.....FOR GLO MODEL
        IF (NGLO.NE.0) THEN
          DO I=1,3
            J=3*(NATOMA(1))+I
            K=3*(NATOMA(1)+1)+I
            I1=NATOMA(1)+1
            I2=NATOMA(1)+2
            PDOT(J)=PDOT(J)-2D0*WS2(I)*W(I1)*Q(J)+WGS2(I)*W(I1)*Q(K)
            PDOT(K)=PDOT(K)-2D0*WG2(I)*W(I2)*Q(K)+WGS2(I)*W(I2)*Q(J)
            GN(K)=GASDEV()*GSW
            PDOT(K)=PDOT(K)-FCG*P(K)+GN(K)
          ENDDO
        ENDIF
!.....END

!.....Ghost anchor Langevin force (Beeman)
        dW_lang = 0.0d0
        if (FCG > 0.0d0) then
           do i = 1, 3
              k = 3 + i
              F_lang = -FCG * P(k) + GASDEV() * GSW_ghost
              PDOT(k) = PDOT(k) + F_lang
              dW_lang = dW_lang + F_lang * (Q(k) - Q_old(k))
           end do
           dW_lang = dW_lang / C1
        end if

!...corrected velocity

        DO K=1,NI
           KK=(K+2)/3
          P(K)=W(KK)*(Q(K)-QTP1(K))/ATIME+(2*PDOT(K)+FT2(K))*ATIME/6
        ENDDO
        
      ENDIF     ! BY BIN, 2017/2/16



      RETURN
      END

!-->FRICTION COEFFICIENTS DEPENDING ON MEAN ELECTRON RADIUS AND THE ATOM TYPE
!-->REFERENCE: PRL, 100, 116102 (2008).
        SUBROUTINE FRICFORCE(NFC,TELEC)
        use venus_params
        use venus_data, only: NFC_H => NFC
        use venus_data
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PDOTT(NDA3),VEL(NDA3),FCOEFT(NATOMA(1)*3,NATOMA(1)*3),FDIAG(NATOMA(1)*3),GRAN(NATOMA(1)*3),VELT(NDA3),FCOEFP(NATOMA(1)*3,NATOMA(1)*3),GRANT(NATOMA(1)*3)
        BOLTZ=0.00198717D0*C1
        FACT=DSQRT(2D0*BOLTZ*TELEC/ATIME)

        IF (NFC.EQ.0) RETURN
        DO K=1,NATOMA(1)*3
          I=(K+2)/3
          VEL(K)=P(K)/W(I)
        ENDDO

        IF (NFC.EQ.1) THEN
!  USE LDFA, THE UNIT OF FRICTION COEFFICIENT IS MASS/TIME
          DO I=1,NATOMA(1)
            J1=I*3-2
            J2=I*3-1
            J3=I*3
            G1=GASDEV
            G2=GASDEV
            G3=GASDEV
            PDOT(J1)=PDOT(J1)-FCOEF(J1,J1)*VEL(J1)+G1*DSQRT(FCOEF(J1,J1))*FACT
            PDOT(J2)=PDOT(J2)-FCOEF(J2,J2)*VEL(J2)+G2*DSQRT(FCOEF(J2,J2))*FACT
            PDOT(J3)=PDOT(J3)-FCOEF(J3,J3)*VEL(J3)+G3*DSQRT(FCOEF(J3,J3))*FACT
          ENDDO
          RETURN
        ENDIF
 
        IF (NFC.EQ.2) THEN
!....WITH ELECTROINC TEMPERATURE, WE HAVE TO FIRST DIAGONALIZE THE
!....FRICTION TENSOR MATRIX, GENERATE THE RANDOM FORCE IN THIS
!....EIGENSPACE, THEN TRANSFORM THE RESULTING FORCE BACK TO THE
!....CARTESIAN SPACE
          PDOTT=0D0
          VELT=0D0
          GRANT=0D0
          DO I=1,NATOMA(1)*3
            IF (FCOEF(I,I).LE.0D0) FCOEF(I,I)=0D0
          ENDDO

          DO K=1,NATOMA(1)*3
            DO J=1,NATOMA(1)*3
              FCOEFT(J,K)=FCOEF(J,K)
            ENDDO
          ENDDO
         
          CALL DIAG2(NATOMA(1)*3,NATOMA(1)*3,FDIAG,FCOEFT)
          
          DO I=1,NATOMA(1)*3
            IF (FDIAG(I).LE.0D0) FDIAG(I)=0D0
            GRAN(I)=GASDEV*DSQRT(FDIAG(I))*FACT  
          ENDDO

          FCOEFP=TRANSPOSE(FCOEFT)

          DO K=1,NATOMA(1)*3
            DO J=1,NATOMA(1)*3
              PDOTT(K)=PDOTT(K)+FCOEFT(J,K)*PDOT(J)
              VELT(K)=VELT(K)+FCOEFT(J,K)*VEL(J)
              GRANT(K)=GRANT(K)+FCOEFT(J,K)*GRAN(J)
            ENDDO
          ENDDO

          DO I=1,NATOMA(1)*3
             PDOTT(I)=PDOTT(I)-FDIAG(I)*VELT(I)+GRANT(I)
          ENDDO

          DO K=1,NATOMA(1)*3
              PDOT(K)=0D0
            DO J=1,NATOMA(1)*3
              PDOT(K)=PDOT(K)+FCOEFP(J,K)*PDOTT(J)
            ENDDO
          ENDDO
!....END

!....WITHOUT ELECTROINC TEMPERATURE
!          PDOTT=PDOT
!          DO I=1,NATOMA(1)*3
!            IF (FCOEF(I,I).LE.0D0) FCOEF(I,I)=0D0
!          ENDDO
!          DO K=1,NATOMA(1)*3
!            DO J=1,NATOMA(1)*3
!              PDOTT(J)=PDOTT(J)-FCOEF(J,K)*VEL(K)
!            ENDDO
!          ENDDO
!          PDOT=PDOTT
!....END
        ENDIF
        END SUBROUTINE FRICFORCE

!----------------------------------------------------------------------
      SUBROUTINE DIAG2(M,N,D,X)
!----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXDIM=1000)
      PARAMETER (EPS=2.5D-16,DINF=2.3D-308,TOL=DINF/EPS)
!
!      COMPUTATION OF ALL EIGENVALUES AND EIGENVECTORS OF A REAL
!      SYMMETRIC MATRIX BY THE METHOD OF QR TRANSFORMATIONS.
!      IF THE EUCLIDEAN NORM OF THE ROWS VARIES   S T R O N G L Y
!      MOST ACCURATE RESULTS MAY BE OBTAINED BY PERMUTING ROWS AND
!      COLUMNS TO GIVE AN ARRANGEMENT WITH INCREASING NORMS OF ROWS.
!
!      TWO MACHINE CONSTANTS MUST BE ADJUSTED APPROPRIATELY,
!      EPS = MINIMUM OF ALL X SUCH THAT 1+X IS GREATER THAN 1 ON THE
!            COMPUTER,
!      TOL = INF / EPS  WITH INF = MINIMUM OF ALL POSITIVE X REPRESEN-
!            TABLE WITHIN THE COMPUTER.
!      A DIMENSION STATEMENT E(160) MAY ALSO BE CHANGED APPROPRIATELY.
!
!      INPUT
!
!      (M)   NOT LARGER THAN 160,  CORRESPONDING VALUE OF THE ACTUAL
!            DIMENSION STATEMENT A(M,M), D(M), X(M,M),
!      (N)   NOT LARGER THAN (M), ORDER OF THE MATRIX,
!      (A)   THE MATRIX TO BE DIAGONALIZED, ITS LOWER TRIANGLE HAS TO
!            BE GIVEN AS  ((A(I,J), J=1,I), I=1,N),
!.....
!.....THE MATRIX #A# HAS BEEN REMOVED FROM THE PROCEDURE
!.....THE MATRIX #X# HAS TO BE PUT UP BY A PREDECESSOR ROUTINE C.....
!
!      OUTPUT
!
!      (D)   COMPONENTS D(1), ..., D(N) HOLD THE COMPUTED EIGENVALUES
!            IN ASCENDING SEQUENCE. THE REMAINING COMPONENTS OF (D) ARE
!            UNCHANGED,
!      (X)   THE COMPUTED EIGENVECTOR CORRESPONDING TO THE J-TH EIGEN-
!            VALUE IS STORED AS COLUMN (X(I,J), I=1,N). THE EIGENVECTORS
!            ARE NORMALIZED AND ORTHOGONAL TO WORKING ACCURACY. THE
!            REMAINING ENTRIES OF (X) ARE UNCHANGED.
!
!      ARRAY (A) IS UNALTERED. HOWEVER, THE ACTUAL PARAMETERS
!      CORRESPONDING TO (A) AND (X)  MAY BE IDENTICAL, ''OVERWRITING''
!      THE EIGENVECTORS ON (A).
!
!      LEIBNIZ-RECHENZENTRUM, MUNICH 1965
!
      DIMENSION   D(M), X(M,M)
      DIMENSION   E(MAXDIM)
!
!     CORRECT ADJUSTMENT FOR IBM 360/91 DOUBLE PRECISION
!
      IF(M.GT.MAXDIM) THEN
        WRITE(6,10) M,MAXDIM
10      FORMAT(' DIMENSION TOO LARGE IN DIAG2:',2I8)
        CALL FEHLER
      END IF
!     CALL ACCNT('DIAG2',1)
!     EPS=DLAMCH('E')
!     TOL=DLAMCH('U')/EPS
!
      IF(N.EQ.1) GO TO 400
      DO 11 I=1,N
      D(I)=0
11    E(I)=0
!
!     HOUSEHOLDER'S REDUCTION
!
      DO 150 I=N,2,-1
      L=I-2
      H=0.0D0
      G=X(I,I-1)
      IF(L.LE.0) GOTO 140
      DO 30 K=1,L
   30 H=H+X(I,K)*X(I,K)
      S=H+G*G
      IF(S.LT.TOL) THEN
        H=0.0D0
      ELSE IF(H.GT.0) THEN
        L=L+1
        F=G
        G=DSQRT(S)
        IF(F.GT.0) G=-G
        H=S-F*G
        X(I,I-1)=F-G
        F=0.0D0
!
        DO 110 J=1,L
        X(J,I)=X(I,J)/H
        S=0.0D0
        DO 80 K=1,J
   80   S=S+X(J,K)*X(I,K)
        J1=J+1
        IF(J1.GT.L) GO TO 100
        DO 90 K=J1,L
   90   S=S+X(K,J)*X(I,K)
  100   E(J)=S/H
  110   F=F+S*X(J,I)
!
        F=F/(2.D0*H)
!
        DO 120 J=1,L
  120   E(J)=E(J)-F*X(I,J)
!
        DO 130 J=1,L
        F=X(I,J)
        S=E(J)
        DO 130 K=1,J
  130   X(J,K)=X(J,K)-F*E(K)-X(I,K)*S
!
      END IF
  140 D(I)=H
  150 E(I-1)=G
!
!     ACCUMULATION OF TRANSFORMATION MATRICES
!
  160 D(1)=X(1,1)
      X(1,1)=1.0D0
      DO 220 I=2,N
      L=I-1
      IF(D(I)) 200,200,170
  170 DO 190 J=1,L
      S=0.0D0
      DO 180 K=1,L
  180 S=S+X(I,K)*X(K,J)
      DO 190 K=1,L
  190 X(K,J)=X(K,J)-S*X(K,I)
  200 D(I)=X(I,I)
      X(I,I)=1.0D0
  210 DO 220 J=1,L
      X(I,J)=0.0D0
  220 X(J,I)=0.0D0
!
!     DIAGONALIZATION OF THE TRIDIAGONAL MATRIX
!
      B=0.0
      F=0.0
      E(N)=0.0D0
!
      DO 340 L=1,N
      H=EPS*(DABS(D(L))+DABS(E(L)))
      IF (H.GT.B) B=H
!
!     TEST FOR SPLITTING
!
      DO 240 J=L,N
      IF (DABS(E(J)).LE.B) GOTO 250
  240 CONTINUE
!
!     TEST FOR CONVERGENCE
!
  250 IF(J.EQ.L) GO TO 340
!
!     SHIFT FROM UPPER 2*2 MINOR
!
  260 P=(D(L+1)-D(L))*0.5D0
      R=DSQRT(P*P+E(L)*E(L))
      IF(P.LT.0) THEN
        P=P+R
      ELSE
        P=P-R
      END IF
  290 H=D(L)+P
      DO 300 I=L,N
  300 D(I)=D(I)-H
      F=F+H
!
!     QR TRANSFORMATION
!
      P=D(J)
      C=1.0D0
      S=0.0D0
!
!     SIMULATION OF LOOP DO 330 I=J-1,L,(-1)
!
      J1=J-1
      DO 330 I=J1,L,-1
      G=C*E(I)
      H=C*P
!
!     PROTECTION AGAINST UNDERFLOW OF EXPONENTS
!
      IF (DABS(P).LT.DABS(E(I))) GOTO 310
      C=E(I)/P
      R=DSQRT(C*C+1.0D0)
      E(I+1)=S*P*R
      S=C/R
      C=1.0D0/R
      GO TO 320
  310 C=P/E(I)
      R=DSQRT(C*C+1.0D0)
      E(I+1)=S*E(I)*R
      S=1.0D0/R
      C=C/R
  320 P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
      DO 330 K=1,N
      H=X(K,I+1)
      X(K,I+1)=X(K,I)*S+H*C
  330 X(K,I)=X(K,I)*C-H*S
!
      E(L)=S*P
      D(L)=C*P
      IF (DABS(E(L)).GT.B) GO TO 260
!
!     CONVERGENCE
!
  340 D(L)=D(L)+F
!
!     ORDERING OF EIGENVALUES
!
      NI=N-1
  350 DO 380 I=1,NI
      K=I
      P=D(I)
      J1=I+1
      DO 360 J=J1,N
      IF(D(J).GE.P) GOTO 360
      K=J
      P=D(J)
  360 CONTINUE
      IF (K.EQ.I) GOTO 380
      D(K) =D(I)
      D(I)=P
      DO 370 J=1,N
      P=X(J,I)
      X(J,I)=X(J,K)
  370 X(J,K)=P
  380 CONTINUE
!
!     FIXING OF SIGN
!
      DO 385 I=1,N
      PM=0
      DO 386 J=1,N
      IF(PM.GT.DABS(X(J,I))) GOTO 386
      PM =DABS(X(J,I))
      K=J
  386 CONTINUE
      IF(X(K,I).GE.0) GOTO 385
      DO 387 J=1,N
  387 X(J,I)=-X(J,I)
  385 CONTINUE
  390 GO TO 410
!
!     SPECIAL TREATMENT OF CASE N = 1
!
  400 D(1)=X(1,1)
      X(1,1)=1.0D0
  410 CONTINUE
!     CALL ACCNT(' ',2)
      RETURN
      END
!----------------------------------------------------------------------
      SUBROUTINE LINEQ(Y,NY,B,NB,X,NX,N,M,INDX,V)
!----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(NY,1),B(NB,1),X(NX,1)
      DIMENSION INDX(N),V(N)
!
!.....LU DECOMPOSITION
!
      CALL LUDCMP(Y,N,NY,INDX,V)
!
      DO 50 J=1,M
      DO 10 I=1,N
10    X(I,J)=B(I,J)
!
!.....NOW SUBSTITUTE AND BACKSUBSTITUTE FOR EACH RHS
!
      CALL LUBKSB(Y,N,NY,INDX,X(1,J))
50    CONTINUE
      RETURN
      END
      SUBROUTINE LUDCMP(A,N,NP,INDX,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP),INDX(N),V(N)
      DATA THR/1.D-15/
!
!.....OBTAIN IMPLICIT SCALING INFORMATION
!
      DO 20 I=1,N
      AMAX=0
        DO 10 J=1,N
        IF(ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
  10    CONTINUE
      IF(AMAX.LT.THR) THEN
        WRITE(6,*) 'AMAX',AMAX
        STOP 'MATRIX SINGULAR'
      ENDIF
      V(I)=1.D0/AMAX
  20  CONTINUE
!
!.....LOOP OVER COLUMS
!
      DO 100 J=1,N
!
!.....FIRST PART FOR I<J
!
        DO 40 I=1,J-1
        SUM=A(I,J)
          DO 30 K=1,I-1
  30      SUM=SUM-A(I,K)*A(K,J)
  40    A(I,J)=SUM
!
!.....SECOND PART FOR I>=J
!
        AMAX=0
        DO 60 I=J,N
        SUM=A(I,J)
          DO 50 K=1,J-1
  50      SUM=SUM-A(I,K)*A(K,J)
        A(I,J)=SUM
        DUM=V(I)*ABS(SUM)
        IF(DUM.GE.AMAX) THEN
          AMAX=DUM
          IMAX=I
        END IF
  60    CONTINUE
!
!.....INTERCHANGE ROWS IF NECESSARY
!
        IF(J.NE.IMAX) THEN
          DO 70 K=1,N
          DUM=A(IMAX,K)
          A(IMAX,K)=A(J,K)
  70      A(J,K)=DUM
          V(IMAX)=V(J)
        END IF
        INDX(J)=IMAX
!
!.....NOW DIVIDE BY THE PIVOT ELEMENT
!
        IF(J.NE.N) THEN
          IF(ABS(A(J,J)).LT.THR) THEN
           WRITE(6,*) 'A(J,J)',A(J,J)
           STOP 'MATRIX SINGULAR'
          ENDIF
          DUM=1.0/A(J,J)
            DO 80 I=J+1,N
  80        A(I,J)=A(I,J)*DUM
        END IF
 100  CONTINUE
      IF(ABS(A(N,N)).LT.THR) THEN
       WRITE(6,*) 'A(N,N)',A(N,N)
       STOP 'MATRIX SINGULAR'
      ENDIF
      RETURN
      END
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,N),INDX(N),B(N)
!
!.....IFIRST WILL BE THE FIRST NONVANISHINNG ELEMENT OF B
!.....FORWARD SUBSTITUTION
!
      IFIRST=0
      DO 20 I=1,N
      LL=INDX(I)
      SUM=B(LL)
      B(LL)=B(I)
      IF(IFIRST.NE.0) THEN
        DO 10 J=1,I-1
!       DO 10 J=II,I-1
  10    SUM=SUM-A(I,J)*B(J)
      ELSE IF(SUM.NE.0) THEN
        IFIRST=I
      END IF
      B(I)=SUM
  20  CONTINUE
!
!.....NOW DO BACKSUBSTITUTION
!
      DO 40 I=N,1,-1
      SUM=B(I)
        DO 30 J=I+1,N
  30    SUM=SUM-A(I,J)*B(J)
  40  B(I)=SUM/A(I,I)
      RETURN
      END
!----------------------------------------------------------------------
      SUBROUTINE FEHLER
      END
