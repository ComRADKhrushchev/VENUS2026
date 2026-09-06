      SUBROUTINE QMMICRO(WW,ANQ,EQNM,NS)

      use venus_params
      use venus_data, only: NS_H => NS
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EQNMM(NS)
      DIMENSION WW(NDA3), ANQ(NDA3)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PSTATE
      INTEGER NS, FLAGRHOTEST, TCOUNT, SCOUNT
      INTEGER  MODE1MAX, MODE2MAX, FLAGCONIC
      LOGICAL  FLAGENV1, FLAGENV2, FLAGENV3

      ZPE=0.0D0
      DO I=1,NS
         ZPE=ZPE+0.5D0*C7*WW(I)
         WW(I)=WW(I)/C6
      ENDDO
      EQNMNOZPE=(EQNM-(ZPE/C1))*CAL2CM

C
C
      FLAGENV1=.TRUE.
      DO WHILE(FLAGENV1)

         FLAGENV2=.FALSE.
         DUM1=EQNMNOZPE

         DO I=1,NS-2

            IF (FLAGENV2) THEN
               EXIT
            ENDIF

            J=NS-I+1
            FLAGENV3=.FALSE.
            TCOUNT=0

            DO WHILE(.NOT.FLAGENV3)
               TCOUNT=TCOUNT+1
               DUM2=DUM1/WW(J)
               RAND=RAND0(ISEED)
               DUM3=DUM2*RAND
 
               ANQ(J)=IDINT(DUM3)
 
               EQNMM(J)=DBLE(ANQ(J))*WW(J)
               DUM2=DUM1
               CALL DENQ(WW,DUM2,NS-I,DUM3)
 
               DUM2=DUM1-EQNMM(J)
               CALL DENQ(WW,DUM2,NS-I,DUM5)
               RAND=RAND0(ISEED)
               IF (DUM5/DUM3 .GT. RAND) THEN
                  FLAGENV3=.TRUE.
               ELSE IF (TCOUNT .GT. 100) THEN
                  FLAGENV2=.TRUE.
                  EXIT
               ENDIF
            ENDDO

            DUM1=DUM1-EQNMM(J)

         ENDDO

         IF (I .GT. NS-2) FLAGENV1=.FALSE.

      ENDDO
C
C     CHOOSING THE ENERGIES FOR REMAINING 2 MODES
C     WITHIN A ERROR OF 100CM-1 (0.01196 IU)
C
      DUM1=DUM1
      DUM2=(DUM1+100.D0)/WW(1)
      MODE1MAX=IDINT(DUM2)
      DUM2=(DUM1+100.D0)/WW(2)
      MODE2MAX=IDINT(DUM2)

      ALLOCATE(PSTATE(MODE1MAX+1,MODE2MAX+1))
      PSTATE(:,:)=0
      SCOUNT=0

      DO I=1,MODE1MAX+1
         DO J=1,MODE2MAX+1
            DUM2=(DBLE(I-1)*WW(1))+(DBLE(J-1)*WW(2))
            DUM3=DUM1-DUM2
            IF (ABS(DUM3) .LT. 100.D0) THEN
               SCOUNT=SCOUNT+1
               PSTATE(I,J)=SCOUNT
            ENDIF
         ENDDO
      ENDDO

      RAND=RAND0(ISEED)

      DUM2=RAND*DBLE(SCOUNT)
      SCOUNT=IDINT(DUM2)+1

      DO I=1,MODE1MAX+1
         DO J=1,MODE2MAX+1
            IF (PSTATE(I,J) .EQ. SCOUNT) THEN
               ANQ(1)=I-1
               ANQ(2)=J-1
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(PSTATE)

      DO I=1,NS
         WW(I)=WW(I)*C6
      ENDDO

      RETURN
      END SUBROUTINE QMMICRO
C
C
C
