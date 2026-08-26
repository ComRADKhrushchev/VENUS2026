      SUBROUTINE TEST
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         CHECK FOR INTERMEDIATE AND FINAL EVENTS
C
C  Kyoyeon 06/07/10
C      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
C     *NPATH,NAST,NABL(NDP),NABM(NDP)
c      COMMON/SELTB/QZ(nda3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
      DIMENSION QCMA(3),VCMA(3),QCMB(3),VCMB(3),QR(3),VR(3)
      CHARACTER*10 TYPE
      CHARACTER*15 COMP
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
      !call ENERGY_1
      !IF (V.GT.VMAX) THEN
        !VMAX=V
        !NCVMAX=NC
        !DO I=1,NDA3
        !QVMAX(I)=Q(I)
        !PVMAX(I)=P(I)
        !ENDDO
      !ENDIF
c
c         TEST FOR REACHING RBAR(I) OR RMAX(I)
c
      NTEST=0
      M=NPATHS+1
      DO I=1,M
         NPATH=I

!-->   Added by Bin 12/18/2013
!-->   The original defination of each path in Venus is not complete and not convinient for equivalent channels
!-->   We use self-defined criteria which is system dependent for each path

c!-->   for F+H2O reaction only 
c       Rfh1=DSQRT((Q(1)-Q(4))**2+(Q(2)-Q(5))**2+(Q(3)-Q(6))**2)
c       Rfh2=DSQRT((Q(1)-Q(10))**2+(Q(2)-Q(11))**2+(Q(3)-Q(12))**2)
c       Rfo=DSQRT((Q(1)-Q(7))**2+(Q(2)-Q(8))**2+(Q(3)-Q(9))**2)
c
c       Roh1=DSQRT((Q(4)-Q(7))**2+(Q(5)-Q(8))**2+(Q(6)-Q(9))**2)
c       Roh2=DSQRT((Q(10)-Q(7))**2+(Q(11)-Q(8))**2+(Q(12)-Q(9))**2)
c
c       if(I.eq.1)then
c         if(Rfh1.ge.RBAR(I).and.Rfh2.ge.RBAR(I)
c     $ .and.Rfo.ge.RBAR(I))NTEST=1
c         if(Rfh1.ge.RMAX(I).and.Rfh2.ge.RMAX(I)
c     $ .and.Rfo.ge.RMAX(I))NTEST=2
c       endif
c
c       if(I.eq.2)then
c         if(Roh1.gt.Roh2.and.Roh1.ge.RBAR(I))NTEST=1
c         if(Roh1.gt.Roh2.and.Roh1.ge.RMAX(I))NTEST=2
c       endif
c            
c       if(I.eq.3)then
c         if(Roh2.gt.Roh1.and.Roh2.ge.RBAR(I))NTEST=1
c         if(Roh2.gt.Roh1.and.Roh2.ge.RMAX(I))NTEST=2
c       endif
!-->   end for F+H2O reaction 

c!-->   for C/Au(111) atomic scattering with slab
c       Check if C adatom (atom 1) has scattered away from the surface.
c       Surface top layer is at z=0 (from Slab.xyz). C atom z = Q(3).
c       Skip first 50 integration steps to avoid triggering on initial position.
c       Also require outward velocity (P(3) > 0) for non-reactive detection.
       if (I.eq.1) then
         rd = Q(3)  ! C atom height above surface (z=0)
         if (NC > 50) then
           if (rd.ge.RBAR(I) .and. P(3).gt.0.0d0) NTEST=1
           if (rd.ge.RMAX(I) .and. P(3).gt.0.0d0) NTEST=2
         end if
       endif

c!-->   end for H2/Cu(111)

!-->   for CH4/Ni(111) dissociative chemisorption only

!       zch4=(W(1)*Q(3)+W(2)*Q(6)+W(3)*Q(9)+W(4)*Q(12)+W(5)*Q(15))
!     &/(W(1)+W(2)+W(3)+W(4)+W(5))
!
!       rch1=DSQRT((Q(13)-Q(1))**2+(Q(14)-Q(2))**2+(Q(15)-Q(3))**2)
!       rch2=DSQRT((Q(13)-Q(4))**2+(Q(14)-Q(5))**2+(Q(15)-Q(6))**2)
!       rch3=DSQRT((Q(13)-Q(7))**2+(Q(14)-Q(8))**2+(Q(15)-Q(9))**2)
!       rch4=DSQRT((Q(13)-Q(10))**2+(Q(14)-Q(11))**2+(Q(15)-Q(12))**2)
!
!       if(I.eq.1)then
!         if(zch4.ge.RBAR(I)) NTEST=1
!         if(zch4.ge.RMAX(I).and.P(5).ge.0d0) NTEST=2
!       endif
!
!       if(I.eq.2)then
!         if(rch1.ge.rch2.and.rch1.ge.rch3
!     $ .and.rch1.ge.rch4.and.rch1.ge.RBAR(I))NTEST=1
!         if(rch1.ge.rch1.and.rch1.ge.rch3
!     $ .and.rch1.ge.rch4.and.rch1.ge.RMAX(I))NTEST=2
!       endif
!
!       if(I.eq.3)then
!         if(rch2.ge.rch1.and.rch2.ge.rch3
!     $ .and.rch2.ge.rch4.and.rch2.ge.RBAR(I))NTEST=1
!         if(rch2.ge.rch1.and.rch2.ge.rch3
!     $ .and.rch2.ge.rch4.and.rch2.ge.RMAX(I))NTEST=2
!       endif
!
!       if(I.eq.4)then
!         if(rch3.ge.rch1.and.rch3.ge.rch2
!     $ .and.rch3.ge.rch4.and.rch3.ge.RBAR(I))NTEST=1
!         if(rch3.ge.rch1.and.rch3.ge.rch2
!     $ .and.rch3.ge.rch4.and.rch3.ge.RMAX(I))NTEST=2
!       endif
!
!       if(I.eq.5)then
!         if(rch4.ge.rch1.and.rch4.ge.rch2
!     $ .and.rch4.ge.rch3.and.rch4.ge.RBAR(I))NTEST=1
!         if(rch4.ge.rch1.and.rch4.ge.rch2
!     $ .and.rch4.ge.rch3.and.rch4.ge.RMAX(I))NTEST=2
!       endif

!-->   end for CH4/Ni(111) dissociative chemisorption 

        IF (NTEST.GT.0) GOTO 3
      ENDDO
      NPATH=1
    3 RETURN
      END
