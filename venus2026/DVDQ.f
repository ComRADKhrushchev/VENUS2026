      SUBROUTINE DVDQ_1
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         CALCULATE POTENTIAL ENERGY PARTIAL DERIVATIVES WITH
C         RESPECT TO COORDINATES (PDOT)
C     All data formerly in COMMON blocks — now in venus_data / venus_params
C
C         ZERO PDOT'S
C
      DO I=1,I3N
         PDOT(I)=0.0D0
      ENDDO

!     ADAPTED BY BIN 8/6/2016
C     NOW: CALCULATE PARTIALS FOR USER CUSTOMED POTENTIAL ENERGY
      IF (NGLO.EQ.0) THEN
         CALL DPESHON(NATOMS)
      ELSE
         CALL DPESHON(NATOMS-1)
      ENDIF
!     END
C


      RETURN
      END
