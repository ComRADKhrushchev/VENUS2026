	SUBROUTINE  PROBJ(T,AI,ISEED,JFINAL)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)

        B=48.5085/(2*AI*T)
        CALL JMAXCALC(T,AI,JMAX,JPEAK) 
        PJMPQ=(2*JPEAK+1)*EXP(-JPEAK*(JPEAK+1)*B)    
10      DUM=RAND0(ISEED)
!-->    modified by Bin, 5/14/2014
        J = NINT(DUM*DBLE(JMAX))        
!-->    end
        PJQ=(2*J+1)*EXP(-J*(J+1)*B)
        PF=PJQ/PJMPQ
        PCOMP=RAND0(ISEED)
        IF(PF.LT.PCOMP) GOTO 10
         JFINAL = J
   
        RETURN
        END



