      SUBROUTINE TEST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         H2/Ag(111) termination: BOTH H atoms above Z_TERM = 8.1 A
C         (the molecule desorbed - matches the venus2026 Z_TERM statute
C         of the matching protocol).  Step-cap termination (NC = NX)
C         stays in the main loop.  fort.78 carries a diagnostic trace
C         (every 200 cycles) for the smoke validation.
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
     *NABL(NDP),NABM(NDP),NPATH,NAST
      ZTERM = 8.1D0
      CALL ENERGY
      NTEST = 0
      IF (Q(3).GT.ZTERM .AND. Q(6).GT.ZTERM) NTEST = 2
      IF (MOD(NC,50).EQ.0) WRITE(78,'(I8,1P4E14.5)') NC, Q(3), Q(6), V,
     &   DSQRT((Q(1)-Q(4))**2+(Q(2)-Q(5))**2+(Q(3)-Q(6))**2)
      RETURN
      END
