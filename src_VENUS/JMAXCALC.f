	SUBROUTINE JMAXCALC(T,AI,JMAX,AJPEAK)
	use venus_data, only: WD1, WD2, WGT
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     COMMON/WNS/WGT — now in venus_data

c       T: Kelvin
c	AI: Amu.A2
	PLIMIT=0.001
	IFIRST=0
	JTOTAL=500
	JMAX=0

C       H-bar**2/KB=48.5085 (FOR AI (IN AMU-A2) AND T IN KELVIN)
	H2KB=48.5085
	IF (T.LE.0.0D0) T = DABS(T)
	B=H2KB/(2.0D0*AI*T)
	TC=H2KB/2.0d0/AI

	IF(WD1.EQ.WD2)THEN
	  NSYM=2        !Bin 2020/2/29
	  WRITE(6,*)'DIATOM IS HOMONUCLEAR'
	ELSE
	  NSYM=1
	  WRITE(6,*)'DIATOM IS HETERONUCLEAR'
	ENDIF

	Q=T/(DBLE(NSYM)*TC)*(1.0D0+1.0D0/3.0D0*TC/T+
     *    1.0D0/15.0D0*(TC/T)**2+4.0D0/315.0D0*(TC/T)**3)

!       Bin, 2020/2/29
	WGT=1.0         ! Bin, nuclear spin factor, 2020/2/29
	DO K=1,JTOTAL+1
	   J=K-1
!       Set up the nuclear spin statics for H2 and D2
	   IF (WD1.EQ.WD2.AND.ABS(WD1-1.008D0).LE.0.01D0) THEN
	      IF (J/2*2.EQ.J) THEN
	         WGT(J)=0.25
	      ELSE
	         WGT(J)=0.75
	      ENDIF
	   ENDIF
	   IF (WD1.EQ.WD2.AND.ABS(WD1-2.014D0).LE.0.01D0) THEN
	      IF (J/2*2.EQ.J) THEN
	         WGT(J)=0.666667
	      ELSE
	         WGT(J)=0.333363
	      ENDIF
	   ENDIF
	   PJ=WGT(J)*(2*J+1)*EXP(-J*(J+1)*B)/Q
	     IF((IFIRST.EQ.0).AND.(PJ.LT.PLIMIT))THEN
	       JMAX=J
	       IFIRST=1
	     ENDIF
	ENDDO
!       end

!-->    modified by Bin, 2020/2/29
	AJPEAK=SQRT(T/2.0d0/TC)-0.5	!fractional peak position of J
!-->    end

	RETURN
	END

