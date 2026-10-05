      SUBROUTINE GWRITE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         WRITE RELEVANT INFORMATION IN OUTPUT FILE
C
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/PARRAY/KR(300),JR(300),KB(300),MB(300),IB(300),IA(300),
     *ITAU(300),ITET(300),IDH(300),IHT(300)
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      COMMON/PRFLAG/NFQP,NCOOR,NFR,NUMR,NFB,NUMB,NFA,NUMA,NFTAU,NUMTAU,
     *NFTET,NUMTET,NFDH,NUMDH,NFHT,NUMHT
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/COORS/R(NDA*(NDA+1)/2),THETA(ND03),ALPHA(ND04),CTAU(ND06),
     *GR(ND08,5),TT(ND09,6),DANG(ND13I)
      COMMON/FORCES/NATOMS,I3N,NFC,NGLO,NA,NLJ,NTAU,NEXP,NGHOST,NTET,
     *NVRR,NVRT,NVTT,NANG,NAXT,NSN2,NRYD,NHFD,NLEPSA,NLEPSB,NDMBE,
     *NRAX,NONB,NMO,NCRCO6
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/RANCOM/RANLST(100),ISEED3(8),IBFCTR
      COMMON/TESTIN/VRELO,INTST
      COMMON/TABLEB/TABLE(42*NDA)
      COMMON/GPATHB/WM(NDA3),TEMP(NDP),AI1D(5),AAI(2),BBI(2),SYMM(5),
     *SYMA,SYMB,GTEMP(NDP),NFLAG(NDP),N1DR,N2DR
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      COMMON/VECTB/VI(4),OAMI(4),AMAI(4),AMBI(4),ETAI,ERAI,ETBI,ERBI
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),
     *QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
C  Kyoyeon 11/25/09
C      common/vrscal/nsel,nscale,nequal,thermotemp,nrgd
      common/vrscal/thermotemp,nsel,nscale,nequal,nrgd

!-->..Added by Bin 04/24/2014
      COMMON/HDIAG/HTMIN,HTMAX,ZASYM,DELM
!-->..End

      DIMENSION HT(300)
C
    1 FORMAT(2X,'THE CYCLE COUNT IS:',I14,16X,'TIME:',F12.3)
    2 FORMAT(2X,'KINETIC ENERGY: ',1PE17.9,'    POTENTIAL ENERGY: ',
     *E17.9/2X,'TOTAL ENERGY:   ',E17.9)
    3 FORMAT(19X,'Q',37X,'P')
    4 FORMAT(F11.6,2F12.6,2X,3F12.6)
    5 FORMAT(4X,' ATOMS  ',5X,'BOND LENGTH(A)')
    6 FORMAT(2X,2I4,10X,F7.3)
    7 FORMAT(5X,'ATOMS',7X,'ANGLE (DEGREES)')
    8 FORMAT(2X,3I3,8X,F8.3)
    9 FORMAT(1X,'XXXXXXXXXXXXXXXXXXXXXXXX TRAJECTORY NUMBER ',I4,
     *' XXXXXXXXXXXXXXXXXXXXXXXXX')
   10 FORMAT(2X,'THE CURRENT RANDOM NUMBER IS: ',8I4,' BASE 256')
   11 FORMAT(2X,'NUMBER',5X,'ALPHA(DEGREES)')
   12 FORMAT(4X,I3,12X,F8.3)
   13 FORMAT(2X,'NUMBER',5X,'TAU(DEGREES)')
   15 FORMAT(2X,'NUMBER',5X,'THETA(TETRAHEDRAL,DEGREES)')
   16 FORMAT(2X,'NUMBER',5X,'DIHEDRAL(DEGREES)')
   17 FORMAT(1X,3F11.7)
   20 FORMAT(5X,'HEIGHT FROM SURFACE FOR ATOMS:',20I4)
   22 FORMAT(2X,'NTZ:',I5,2X,'NC:',I8,2X,'HEIGHT:',20F9.3)
C
C         WRITE RELEVANT INFORMATION IN CHECKPOINT FILE
C
!   removed by bin, 2016/7/30
C      OPEN(50,FORM='UNFORMATTED')
C      REWIND(50)
C      WRITE(50)Q,P,QDOT,PDOT,TABLE,VRELO,RANLST,GTEMP,NFLAG,ISEED0,
C     *ISEED3,NX,NC,NTZ,INTST,NAST,IBFCTR,
C     *VI,OAMI,AMAI,AMBI,ETAI,ERAI,ETBI,ERBI
C      CLOSE(50)
!   end bin
C
C         WRITE TRAJECTORY INFORMATION
C
      CNC=NC
      TI=TIME*CNC
      WRITE(6,9)NTZ
      WRITE(6,1)NC,TI

      IF (NSELT.EQ.2.OR.NSELT.EQ.3) WRITE(6,10)(ISEED0(9-I),I=1,8)
      WRITE(6,2)T,V,H

!-->      Added by Bin, 12/23/2013
C      if(H.lt.HTMIN) HTMIN=H
C      if(H.gt.HTMAX) HTMAX=H
!-->      End

C
C         WRITE COORDINATES FOR GRAPHICS
C
      IF (NCOOR.EQ.1.AND.NSELT.NE.-2) THEN
         WRITE(8,9)NTZ
         WRITE(8,17)(Q(I),I=1,NATOMS*3)
!-->    added by Bin, 4/24/2014
        write(1000+NTZ,'(i6)')NATOMA(1)+NATOMB(1)
        write(1000+NTZ,'(f17.9,f10.4,f17.9)')V,NC*0.1,H
        do iatom=1,1
        write(1000+NTZ,'(f6.1,6f10.5)')W(iatom),Q(3*iatom-2:3*iatom),
     &P(3*iatom-2:3*iatom)/W(iatom)
        enddo
        do iatom=2,2
        write(1000+NTZ,'(f6.1,6f10.5)')W(iatom),Q(3*iatom-2:3*iatom),
     &P(3*iatom-2:3*iatom)/W(iatom)
        enddo
        do iatom=3,3
        write(1000+NTZ,'(f6.1,6f10.5)')W(iatom),Q(3*iatom-2:3*iatom),
     &P(3*iatom-2:3*iatom)/W(iatom)
        enddo

!-->    end Bin's changes

      ENDIF
C
C         WRITE COORDINATES AND MOMENTA 
C
      IF (NFQP.NE.0) THEN
         WRITE(6,3)
         J=1
         DO L=1,NATOMS
            M=J+2
            WRITE(6,4)(Q(I),I=J,M),(P(I),I=J,M)
            J=J+3
         ENDDO
      ENDIF
!!!    removed by Bin, 2016/7/30
cC
cC         CALCULATE AND WRITE POSSIBLE INTERATOMIC DISTANCES
cC
c      IF (NFR.NE.0) THEN
c         WRITE(6,5)
c         DO I=1,NUMR
c            J3=3*JR(I)
c            J2=J3-1
c            J1=J2-1
c            K3=3*KR(I)
c            K2=K3-1
c            K1=K2-1
c            T1=Q(K1)-Q(J1)
c            T2=Q(K2)-Q(J2)
c            T3=Q(K3)-Q(J3)
c            RR=SQRT(T1*T1+T2*T2+T3*T3)
c            WRITE(6,6)JR(I),KR(I),RR
c         ENDDO
c      ENDIF
cC
cC         CALCULATE AND WRITE POSSIBLE ANGLES 
cC
c      IF (NFB.NE.0) THEN
c         WRITE(6,7)
c         DO I=1,NUMB
c            K3=3*KB(I)
c            K2=K3-1
c            K1=K2-1
c            M3=3*MB(I)
c            M2=M3-1
c            M1=M2-1
c            I3=3*IB(I)
c            I2=I3-1
c            I1=I2-1
c            T1=Q(I1)-Q(M1)
c            T2=Q(I2)-Q(M2)
c            T3=Q(I3)-Q(M3)
c            T4=Q(K1)-Q(M1)
c            T5=Q(K2)-Q(M2)
c            T6=Q(K3)-Q(M3)
c            R1=SQRT(T1*T1+T2*T2+T3*T3)
c            R2=SQRT(T4*T4+T5*T5+T6*T6)
c            CTHETA=(T1*T4+T2*T5+T3*T6)/R1/R2
c            IF (CTHETA.GT. 1.00D0) CTHETA= 1.00D0
c            IF (CTHETA.LT.-1.00D0) CTHETA=-1.00D0
c            DUM=ACOS(CTHETA)/C4
c            WRITE(6,8)KB(I),MB(I),IB(I),DUM
c         ENDDO
c      ENDIF
cC
cC         WRITE ALPHA, TAU, TETRAHEDRAL AND DIHEDRAL ANGLES
cC         THESE ARE POTENTIAL FUNCTION - DEPENDENT
cC
c      IF (NFA.NE.0) THEN
c         WRITE(6,11)
c         DO I=1,NUMA
c            DUM=ALPHA(IA(I))/C4
c            WRITE(6,12)IA(I),DUM
c         ENDDO
c      ENDIF
c      IF (NFTAU.NE.0) THEN
c         WRITE(6,13)
c         DO I=1,NUMTAU
c            DUM=ACOS(CTAU(ITAU(I)))/C4
c            WRITE(6,12)ITAU(I),DUM
c         ENDDO
c      ENDIF
c      IF (NFTET.NE.0) THEN
c         WRITE(6,15)
c         DO I=1,NUMTET
c            DUM=TT(1,ITET(I))/C4
c            WRITE(6,12)ITET(I),DUM
c         ENDDO
c      ENDIF
c      IF (NFDH.NE.0) THEN
c         WRITE(6,16)
c         DO I=1,NUMDH
c            DUM=DANG(IDH(I))/C4
c            WRITE(6,12) IDH(I),DUM
c         ENDDO
c      ENDIF
cC
cC          CALCULATE POSITION HEIGHTS ABOVE THE SURFACE FOR SOME ATOMS
cC          SAVE DATA IN FORT.17
cC
c      IF (NFHT.NE.0) THEN
cc         WRITE(17,20)(IHT(I),I=1,NUMHT)
c         DO I=1,NUMHT
c            CALL HEIGHT(IHT(I),HT(I))
c         ENDDO
c         WRITE(17,22)NTZ,NC,(HT(I),I=1,NUMHT)           
c      ENDIF
c
c          calculate system temperature
c
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
C
!!!   end bin
      CALL FLUSH(6)
C
      RETURN
      END
