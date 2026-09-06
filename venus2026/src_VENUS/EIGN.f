      SUBROUTINE EIGN(A,VEC,NN,RHO)
      use venus_params
      use venus_data, only: EIG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         DIAGONALIZE A MATRIX A, OF WHICH ONLY LOWER TRIANGLE IS USED
C         AND DESTROYED, USING THE GIVENS-HOUSHOLDER ALGORITHM.
C         EIGENVALUES ARE RETURNED IN ALGEBRAIC ASCENDING ORDER IN ARRAY
C         EIG THE EIGENVECTORS ARE RETURNED IN VEC.   
C
C           PARAMETERS PASSED                   
C           RHO IS THE UPPER LIMIT FOR OFF-DIAGONAL     
C           NN IS THE SIZE OF THE MATRIX              
C
      DIMENSION WK(NDA3),BETASQ(NDA3),GAMMA(NDA3),BETA(NDA3)
      DIMENSION A(NN,NN),VEC(NN,NN)
      DIMENSION PK(NDA3),QK(NDA3),IPOSV(NDA3),IVPOS(NDA3),IORD(NDA3)
      RHOSQ=RHO*RHO
      N=NN
      IF (N.EQ.0) RETURN
      SHIFT = 0.0D0
      N1=N-1
      N2=N-2
      GAMMA(1)=A(1,1)
      IF (N2.NE.0) THEN 
         IF (N2.LT.0) GOTO 280
         DO NR=1,N2
            B=A(NR+1,NR)
            S=0.0D0
            DO I=NR,N2
               S=S+A(I+2,NR)**2
            ENDDO
C
C         PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION
C
            A(NR+1,NR)=0.0D0
            IF (S.GT.0) THEN
               S=S+B*B
               SGN=+1.0D0
               IF (B.LT.0) SGN = -1.D0
               SQRTS=SQRT(S)
               D=SGN/(SQRTS+SQRTS)
               TEMP=SQRT(0.5D0+B*D)
               WK(NR)=TEMP
               A(NR+1,NR)=TEMP
               D=D/TEMP
               B=-SGN*SQRTS
C
C         D IS FACTOR OF PROPORTIONALITY. COMPUTE AND SAVE W VECTOR
C
               DO I=NR,N2
                  TEMP=D*A(I+2,NR)
                  WK(I+1)=TEMP
                  A(I+2,NR)=TEMP
               ENDDO
C
C         PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR.
C         SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K).
C
               WTAW = 0.0D0
               DO I=NR,N1
                  SUM = 0.0D0
                  DO J=NR,I
                     SUM=SUM+A(I+1,J+1)*WK(J)
                  ENDDO
                  I1=I+1
                  IF ((N1-I1).GE.0) THEN
                     DO J=I1,N1
                        SUM=SUM+A(J+1,I+1)*WK(J)
                     ENDDO
                  ENDIF
                  PK(I)=SUM
                  WTAW=WTAW+SUM*WK(I)
               ENDDO
C
C         P VECTOR AND SCALAR K NOW STORED. NEXT COMPUTE Q VECTOR
C         AND FORM PAP MATRIX.
C
               DO I=NR,N1
                  QK(I)=PK(I)-WTAW*WK(I)
               ENDDO
               DO J=NR,N1
                  QJ=QK(J)
                  WJ=WK(J)
                  DO I=J,N1
                     A(I+1,J+1)=A(I+1,J+1)-2.*(WK(I)*QJ+WJ*QK(I))
                  ENDDO
               ENDDO
  250          BETA(NR)=B
            ENDIF
            BETASQ(NR)=B*B
            GAMMA(NR+1)=A(NR+1,NR+1)
         ENDDO
      ENDIF
      B=A(N,N-1)
      BETA(N-1)=B
      BETASQ(N-1)=B*B
      GAMMA(N)=A(N,N)
  280 BETASQ(N)=0.
C
C         ADJOIN AN IDENTITY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS
C
      DO I=1,N
         DO J=1,N
            VEC(I,J)=0.0D0
         ENDDO
         VEC(I,I)=1.0D0
      ENDDO
      M=N
      SUM=0.0D0
      NPAS=1
      GOTO 400
  310 SUM=SUM+SHIFT
      COSA=1.
      G=GAMMA(1)-SHIFT
      PP=G
      PPBS=PP*PP+BETASQ(1)
      PPBR=SQRT(PPBS)
      DO J=1,M
         COSAP=COSA
         IF (PPBS.EQ.0) THEN
            SINA = 0.0D0
            SINA2=0.0D0
            COSA=1.0D0
         ELSE
            SINA=BETA(J)/PPBR
            SINA2=BETASQ(J)/PPBS
            COSA=PP/PPBR
C
C         POSTMULTIPLY IDENTITY BY P-TRANSPOSE
C
            NT=J+NPAS
            IF (NT.GT.N) NT=N
            DO I=1,NT
               TEMP=COSA*VEC(J,I)+SINA*VEC(J+1,I)
               VEC(J+1,I)=-SINA*VEC(J,I)+COSA*VEC(J+1,I)
               VEC(J,I)=TEMP
            ENDDO
         ENDIF
         DIA=GAMMA(J+1)-SHIFT
         U=SINA2*(G+DIA)
         GAMMA(J)=G+U
         G=DIA-U
         PP=DIA*COSA-SINA*COSAP*BETA(J)
         IF (J.EQ.M) THEN
            BETA(J)=SINA*PP
            BETASQ(J)=SINA2*PP*PP
            GOTO 340
         ENDIF
         PPBS=PP*PP+BETASQ(J+1)
         PPBR=SQRT(PPBS)
         BETA(J)=SINA*PPBR
         BETASQ(J)=SINA2*PPBS
      ENDDO
  340 GAMMA(M+1)=G
C
C         TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT
C
      NPAS=NPAS+1
      IF (BETASQ(M).GT.RHOSQ) GOTO 410
  390 EIG(M+1)=GAMMA(M+1)+SUM
  400 BETA(M)=0.0D0
      BETASQ(M)=0.0D0
      M=M-1
      IF (M.EQ.0) GOTO 430
      IF (BETASQ(M).LE.RHOSQ) GOTO 390
  410 CONTINUE 
C
C         TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN
C         VALUE AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT
C
      A2=GAMMA(M+1)
      R2=0.5D0*A2
      R1=0.5D0*GAMMA(M)
      R12=R1+R2
      DIF=R1-R2
      TEMP=SQRT(DIF*DIF+BETASQ(M))
      R1=R12+TEMP
      R2=R12-TEMP
      DIF=ABS(A2-R1)-ABS(A2-R2)
      IF (DIF.GE.0) THEN
         SHIFT=R2
      ELSE
         SHIFT=R1
      ENDIF
      GOTO 310
  430 EIG(1)=GAMMA(1)+SUM
C
C         INITIALIZE AUXILIARY TABLES REQUIRED FOR
C         REARANGING THE VECTORS
C
      DO J=1,N
         IPOSV(J)=J
         IVPOS(J)=J
         IORD(J) = J
      ENDDO
C
C         USE A TRANSPOSITON SORT TO ORDER THE EIGENVALUES
C
      M=N
      GOTO 470
  450 DO J=1,M
         IF (EIG(J).GT.EIG(J+1)) THEN
            TEMP=EIG(J)
            EIG(J)=EIG(J+1)
            EIG(J+1)=TEMP
            ITEMP=IORD(J)
            IORD(J)=IORD(J+1)
            IORD(J+1)=ITEMP
         ENDIF
      ENDDO
  470 M=M-1
      IF (M.NE.0) GOTO 450
      IF (N1.NE.0) THEN
         DO IL=1,N1
            NV=IORD(IL)
            NP=IPOSV(NV)
            IF (NP.NE.IL) THEN
               LV=IVPOS(IL)
               IVPOS(NP)=LV
               IPOSV(LV)=NP
               DO I=1,N
                  TEMP=VEC(IL,I)
                  VEC(IL,I)=VEC(NP,I)
                  VEC(NP,I) = TEMP
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
C         BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX
C
      DO NRR=1,N
         K=N1
  510    K=K-1
         IF (K.GT.0) THEN
            SUM = 0.0
            DO I=K,N1
               SUM=SUM+VEC(NRR,I+1)*A(I+1,K)
            ENDDO
            SUM=SUM+SUM
            DO I=K,N1
               VEC(NRR,I+1)=VEC(NRR,I+1)-SUM*A(I+1,K)
            ENDDO
            GOTO 510
         ENDIF
      ENDDO
      RETURN
      END
