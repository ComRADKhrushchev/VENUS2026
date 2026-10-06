        SUBROUTINE JMAXCALC(T,AI,JMAX,JPEAK)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
        COMMON/WNJ/WD1,WD2

c       T: Kelvin
c	AI: Amu.A2
	PLIMIT=0.001
	IFIRST=0
        JTOTAL=50

C       H**2/KB=48.5085 (FOR AI (IN AMU-A2) AND T IN KELVIN)
        H2KB=48.5085
        B=H2KB/(2.0D0*AI*T)
	TC=H2KB/2.0d0/AI

        IF(WD1.EQ.WD2)THEN
          NSYM=1
          WRITE(6,*)'DIATOM IS HOMONUCLEAR'
        ELSE
          NSYM=1
          WRITE(6,*)'DIATOM IS HETERONUCLEAR'
        ENDIF

        Q=T/(DBLE(NSYM)*TC)*(1.0D0+1.0D0/3.0D0*TC/T+
     *    1.0D0/15.0D0*(TC/T)**2+4.0D0/315.0D0*(TC/T)**3) 

	DO K=1,JTOTAL+1
           J=K-1
	   PJ=(2*J+1)*EXP(-J*(J+1)*B)/Q
             IF((IFIRST.EQ.0).AND.(PJ.LT.PLIMIT))THEN 
               JMAX=J
               IFIRST=1
             ENDIF	
	ENDDO
!-->    modified by Bin, 5/14/2014
	JPEAK=NINT(SQRT(T/2.0d0/TC))
!-->    end

	RETURN
	END



