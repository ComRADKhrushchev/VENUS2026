      SUBROUTINE PRINFO(NFQP, NCOOR)
      use venus_params
      use venus_data, only: NFQP_H => NFQP, NCOOR_H => NCOOR
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer, intent(in) :: NFQP, NCOOR
C
C    INFORMATION TO BE PRINTED
C
  816 FORMAT('   NFQP=',I2,'  NCOOR= ',I2/)
C
      WRITE(6,816)NFQP,NCOOR
C
C         NFQP=0, DO NOT PRINT Q AND P ARRAYS
C         NCOOR=0, DO NOT WRITE COORDINATES INTO UNIT 8
C
!....removed by bin, 2016/7/30
c      READ(5,*)NFR,NUMR
c      WRITE(6,817)NFR,NUMR
c      IF (NFR.NE.0) THEN
c         READ(5,*)(JR(I),KR(I),I=1,NUMR)
c         WRITE(6,*)
c         DO I=1,NUMR
c            WRITE(6,818)JR(I),KR(I)
c         ENDDO
c         WRITE(6,*)
c      ENDIF
cC
c      READ(5,*)NFB,NUMB
c      WRITE(6,819)NFB,NUMB
c      IF (NFB.NE.0) THEN
c         WRITE(6,*)
c         READ(5,*)(KB(I),IB(I),MB(I),I=1,NUMB)
c         WRITE(6,829)(KB(I),IB(I),MB(I),I=1,NUMB)
c         WRITE(6,*)
c      ENDIF
cC
c      READ(5,*)NFA,NUMA
c      WRITE(6,822)NFA,NUMA
c      IF (NFA.NE.0) THEN
c         READ(5,*)(IA(I),I=1,NUMA)
c         WRITE(6,823)
c         WRITE(6,821)(IA(I),I=1,NUMA)
c         WRITE(6,*)
c      ENDIF
cC
c      READ(5,*)NFTAU,NUMTAU
c      WRITE(6,824)NFTAU,NUMTAU
c      IF (NFTAU.NE.0) THEN
c         READ(5,*)(ITAU(I),I=1,NUMTAU)
c         WRITE(6,825)
c         WRITE(6,821)(ITAU(I),I=1,NUMTAU)
c      ENDIF
cC
c      READ(5,*)NFTET,NUMTET
c      WRITE(6,947)NFTET,NUMTET
c      IF (NFTET.NE.0) THEN
c         READ(5,*)(ITET(I),I=1,NUMTET)
c         WRITE(6,948)
c         WRITE(6,821)(ITET(I),I=1,NUMTET)
c      ENDIF
cC
c      READ(5,*) NFDH,NUMDH
c      WRITE(6,953) NFDH,NUMDH
c      IF (NFDH.NE.0) THEN
c         READ(5,*)(IDH(I),I=1,NUMDH)
c         WRITE(6,954)
c         WRITE(6,821)(IDH(I),I=1,NUMDH)
c      ENDIF
c      WRITE(6,*)
cC
c      READ(5,*)NFHT,NUMHT
c      WRITE(6,990)NFHT,NUMHT
c      IF (NFHT.NE.0) THEN
c         READ(5,*)(IHT(I),I=1,NUMHT)
c         WRITE(6,992)
c         WRITE(6,994)(IHT(I),I=1,NUMHT)
c      ENDIF
c      WRITE(6,*)
c....end bin
cC
      RETURN
      END
