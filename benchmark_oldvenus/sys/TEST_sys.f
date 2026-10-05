      SUBROUTINE TEST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         CHECK FOR INTERMEDIATE AND FINAL EVENTS
C
C  Kyoyeon 06/07/10
C      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
C     *NPATH,NAST,NABL(NDP),NABM(NDP)
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/PSN2/PESN2,GA,RA,RB
      COMMON/FORCES/NATOMS,I3N,NFC,NGLO,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/WASTE/QQ(NDA3),PP(NDA3),WX,WY,WZ,L(NDA),NAM
      COMMON/TESTIN/VRELO,INTST
      COMMON/TESTSN2/GAO,NSAD,NCBA,NCAB,IBAR
      COMMON/FINALB/EROTA,EROTB,EA(3),EB(3),AMA(4),AMB(4),AN,AJ,BN,BJ,
     *OAM(4),EREL,ERELSQ,ETCM,BF,SDA,SDB,DELH(NDP),ANG(NDG),NFINAL
      COMMON/VMAXB/QVMAX(NDA3),PVMAX(NDA3),VMAX,NCVMAX      
c      COMMON/SELTB/QZ(nda3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      DIMENSION QCMA(3),VCMA(3),QCMB(3),VCMB(3),QR(3),VR(3)
      CHARACTER*10 TYPE
      CHARACTER*15 COMP
      REAL*8 D_DIS3(3),XXX(3),YYY(3),ZZZ(3),XX0,YY0,ZZ0
      INTEGER*4 I_PA3(3,2),II,ICHANNEL,I1,I2,J1,J2
      LOGICAL CH_ONCE
      SAVE CH_ONCE, I_PA3, ICHANNEL
      DATA CH_ONCE/.FALSE./
C
 900  FORMAT(4X,' TURNING POINT #  ','  CYCLE ','  RCM(A)',
     *       '    EA     ','    EB     ','    EROTA  ','    EROTB  ',
     *       '    JA     ','    JB     ','    L      ')
 903  FORMAT(8X,A10,I4,I8,F8.3,1P7D11.4)
 905  FORMAT(/5X,'$$$$BARRIER CROSSING NUMBER$$$$ ',I6,
     &       '  AT CYCLE',I8)
 906  FORMAT(5X,'$$$$BARRIER CROSSING FROM B TO A $$$$')
 907  FORMAT(5X,'$$$$BARRIER CROSSING FROM A TO B $$$$')
 935  FORMAT(7X,'RA= ',F7.3,3X,'RB= ',F7.3,3X,'GA= ',F7.3)
 910  FORMAT(4X,' TURNING POINT #  ',5X,' COMPLEX ',6X,
     *'  CYCLE ','  RCM(A)',
     *       '    EA     ','    EB     ','    EROTA  ','    EROTB  ',
     *       '    JA     ','    JB     ','    L      ')
 913  FORMAT(8X,A10,I4,2X,A17,I8,F8.3,1P7D11.4)
C
C         TEST FOR MAXIMUM IN POTENTIAL ENERGY
C
      CALL ENERGY      
      IF (V.GT.VMAX) THEN
        VMAX=V
        NCVMAX=NC
        DO I=1,NDA3
        QVMAX(I)=Q(I)
        PVMAX(I)=P(I)
        ENDDO
      ENDIF
c
c         TEST FOR REACHING RBAR(I) OR RMAX(I)
c
      NTEST=0
      M=NPATHS+1
      DO I=1,M
         NPATH=I

!#########################################################
! Added by Zexing Qu, 2023.7.13
! For 3 atoms reactions,

       IF (.NOT. CH_ONCE) THEN
          OPEN(23,FILE='channel')
          READ(23,*)ICHANNEL
          I_PA3=0
          DO II=1,ICHANNEL
            READ(23,*)I_PA3(II,1),I_PA3(II,2)
          ENDDO
          CLOSE(23)
          CH_ONCE=.TRUE.
       ENDIF
       
       D_DIS3=0.0D0
       XXX=0.0D0
       YYY=0.0D0
       ZZZ=0.0D0
       DO II=1,ICHANNEL
         I1=I_PA3(II,1)
         I2=I_PA3(II,2)
         J1=3*(I1-1)+1
         J2=3*(I2-1)+1
         XXX(II)=(W(I1)*Q(J1  )+W(I2)*Q(J2  ))/(W(I1)+W(I2))
         YYY(II)=(W(I1)*Q(J1+1)+W(I2)*Q(J2+1))/(W(I1)+W(I2))
         ZZZ(II)=(W(I1)*Q(J1+2)+W(I2)*Q(J2+2))/(W(I1)+W(I2))
       ENDDO
       
       DO II=1,ICHANNEL
         DO I1=1,3
           IF(I1/=I_PA3(II,1) .AND. I1/=I_PA3(II,2)) THEN
             J1=I1
           ENDIF
         ENDDO
         J2=3*(J1-1)+1
         XX0=Q(J2)
         YY0=Q(J2+1)
         ZZ0=Q(J2+2)
         D_DIS3(II)=(XXX(II)-XX0)**2+(YYY(II)-YY0)**2+(ZZZ(II)-ZZ0)**2
         D_DIS3(II)=DSQRT(D_DIS3(II))
       ENDDO

       if(I.eq.1)then
         if(D_DIS3(1).ge.RBAR(I))NTEST=1
         if(D_DIS3(1).ge.RMAX(I))NTEST=2
       endif

       if(I.eq.2)then
         if(D_DIS3(2).ge.RBAR(I))NTEST=1
         if(D_DIS3(2).ge.RMAX(I))NTEST=2
       endif

       if(I.eq.3)then
         if(D_DIS3(3).ge.RBAR(I))NTEST=1
         if(D_DIS3(3).ge.RMAX(I))NTEST=2
       endif

! --> end of 3 atoms reactions
!##########################################################

        IF (NTEST.GT.0) GOTO 3
      ENDDO
      NPATH=1
    3 RETURN
      END
