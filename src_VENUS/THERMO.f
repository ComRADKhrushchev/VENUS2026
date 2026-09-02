      SUBROUTINE THERMO(IPR,info)
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      COMMON/VRSCAL/NSEL,NSCALE,NEQUAL,THERMOTEMP,NRGD
C
      integer   info    !!!  added by Qinghua Liu 2019-9-15
      real  average_temp
      SVELSQ=0.00D0
C       
        average_temp=0.0d0
C    MAKE SURE THAT THE COORDINATES ARE NOT WRITTEN DURING MINIZATION
C
c        info=1
c      write(300,*)info

      N=NATOMB(1)-NRGD
C	
C    CALCULATE VELSQ.  RECALL VNEW = VOLD *SQRT(AHEAT/VELSQ)
C
	DO I=1,N
	  J3 = 3 * LB(1,I)
	  J2 = J3 - 1
	  J1 = J2 - 1
          J = LB(1,I)
          SVELSQ =SVELSQ+(P(J1)**2+P(J2)**2+P(J3)**2)/W(J)
	END DO
C
C   REFERENCE:  "MOLECULAR DYNAMICS SIMULATION" BY JIM HAILE  P.458
C   K = 1.38066 * 10(-23) J/K = 0.00198624 KCAL/MOL K
C   FACTEMP IS AHEAT ABOVE
C
        If(info.eq.0) THEN  !!! added by Qinghua Liu 2019-9-15
          FACTEMP =3.0D0*DBLE(N)*0.00198717D0*THERMOTEMP*C1
          FACTOR = SQRT( FACTEMP / SVELSQ )
          TEMPINITS=SVELSQ/(3.0 * DBLE(N) * 0.00198717D0 * C1)
          IF (MOD(IPR,100).EQ.0.AND.NSEL.EQ.1) THEN
             WRITE(6,'(1X,A,F12.6)')'SYSTEM TEMP',
     *             TEMPINITS
          ENDIF
          DO I=1,N
            J3 = 3 * LB(1,I)
            J2 = J3 - 1
            J1 = J2 - 1
            P( J1 ) = P( J1 ) * FACTOR
            P( J2 ) = P( J2 ) * FACTOR
            P( J3 ) = P( J3 ) * FACTOR
          END DO
        
        ENDif  !!! added by Qinghua Liu 2019-9-15


!!!  added by Qinghua Liu for Andersen thermostat  2019-9-15 
        If(info.eq.1) then
          DESKET=sqrt(0.00198717D0*THERMOTEMP*C1)
          write(300,*)"target temperature= ",thermotemp
           DO I = 1,N
              J3=3*LB(1,I)
              J2=J3-1
              J1=J2-1
              J=LB(1,I)
              debug = gasdev()
              P(J1)=GASDEV()*DESKET*SQRT(W(J))
              P(J2)=GASDEV()*DESKET*SQRT(W(J))
              P(J3)=GASDEV()*DESKET*SQRT(W(J))
           ENDDO
         Endif 
!!! added by Qinghua Liu for Andersen thermostat  2019-9-15



C
C    MAKE SURE THAT THE TOTAL LINEAR MOMENTUM EQUALS ZERO
C
cc changed by qinghua Liu 2020-1-1
!        SUMX = 0.00
!        SUMY = 0.00
!        SUMZ = 0.00
!        DO I=1,N
!          J3 = 3 * LB(1,I)
!          J2 = J3 - 1
!          J1 = J2 - 1
!          SUMX = SUMX + P( J1 )
!          SUMY = SUMY + P( J2 )
!          SUMZ = SUMZ + P( J3 )
!        END DO
!        SUMX = SUMX / REAL( N )
!        SUMY = SUMY / REAL( N )
!        SUMZ = SUMZ / REAL( N )
!        DO I=1,N
!          J3 = 3 * LB(1,I)
!          J2 = J3 - 1
!          J1 = J2 - 1
!          P(J1) = P(J1) - SUMX
!          P(J2) = P(J2) - SUMY
!          P(J3) = P(J3) - SUMZ
!        END DO
cc ended by qinghua Liu 2020-1-1
C
C     CALCULATE THE TEMPERATURE
C
      SVELSQ=0.00
      DO I=1,N
         J3 = 3 * LB(1,I)
         J2 = J3 - 1
         J1 = J2 - 1                                                    
         J = LB(1,I)
         SVELSQ = SVELSQ+(P(J1)**2+P(J2)**2+P(J3)**2)/W(J)
      END DO
C
C  COMPUT THE TEMPERATURE AFTER THE RESCALING. 
C
      TEMPINITS=SVELSQ/(3.0 * DBLE(N) * 0.00198717D0 * C1)
c        write(*,*)NC,NSCALE
      IF (MOD(IPR,100).EQ.0.AND.NSEL.EQ.1.AND.IPR.gt.NSCALE) THEN
      average_temp=average_temp+tempinits
      write(6,'(A30,f12.6,A30,f12.6)')"target temperature=", 
     &          thermotemp, "temperature resample  = ",  TEMPINITS
      ENDIF        
      if(IPR.eq.(NSCALE+NEQUAL))then
         write(6,*)"average_temperature= ", average_temp/NT
      endif
c      IF (MOD(IPR,100).EQ.0.AND.NSEL.EQ.1) THEN
c         WRITE(6,'(1X,A,F12.6)')'SYSTEM TEMP AFTER RESCALING',
c     *        TEMPINITS
c      ENDIF
 999  RETURN
      END

