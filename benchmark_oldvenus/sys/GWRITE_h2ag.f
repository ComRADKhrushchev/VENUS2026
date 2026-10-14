      SUBROUTINE GWRITE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         H2/Ag(111) per-call record for post-processing: cycle, r(HH),
C         the two atom heights, current potential.  One row per call on
C         fort.77; the analysis takes the max-NC row of each trajectory
C         (rows group by NC dropping).
C
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      RHH = DSQRT((Q(1)-Q(4))**2+(Q(2)-Q(5))**2+(Q(3)-Q(6))**2)
      OPEN(97,FILE='traj_rows.txt',POSITION='APPEND')
      WRITE(97,'(I8,F12.6,1P3E16.7)') NC, RHH, Q(3), Q(6), V
      RETURN
      END
