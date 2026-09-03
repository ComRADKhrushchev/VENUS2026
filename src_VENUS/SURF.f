      SUBROUTINE SURF(NSURF)
      use venus_params
      use venus_data, only: NSURF_H => NSURF, PVMAX_H => PVMAX
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         BEAM TYPE AND THERMAL TYPE SURFACE SAMPLING 
C
C      COMMON/SURFB/NN1,NN2,NN3,NN4,THTA,PHI,NCHI,CHI,RX0,RY0,RZ0,
C     *NTHET,THET,NPHI1,PHI1,NPHI2,PHI2

C  ADDED BY BIN 2016/10/1
      DIMENSION EDGE(3,4)
C  ADDED BY ZXY 20170620
C
  10  FORMAT(/15X,'IMPACT PARAMETER=',F7.3,' A')
  12  FORMAT(/5X,'BEAM-TYPE STUDY: THETA=',F8.3,'  PHI=',F8.3,
     *'  KAI=',F8.3,' (DEG)'/)
  13  FORMAT(/5X,'GAS-SURFACE SCATTERING: THETA=',F8.3,'  PHI=',F8.3,
     *' (DEG)'/)
  14  FORMAT(/5X,'THERMAL-TYPE STUDY: THETA=',F8.3,'  PHI1=',F8.3,
     *'  PHI2=',F8.3,' (DEG)'/)
  16  FORMAT(/5X,'RELATIVE TRANSLATIONAL ENERGY SELECTED: ',F7.2,
     *' KCAL/MOL'/)
  18  FORMAT(/5X,'CHOSEN:   LX,LY,LZ =',1P3D13.5,' H-BAR')
C
C     REFERENCE Z-AXIS VECTOR IS ORTHOGONAL TO THE SURFACE PLANE
C
!   ADDED BY BIN, 2016/7/30
      IF (NSURF.EQ.2) THEN
C        F24/D3 fix: the unit cell (NN1 the origin, NN2/NN3/NN4 in
C        anticlockwise sequence, NN1-NN3 defining the reference plane) is
C        given by the first four atoms of fragment B -- the rigid slab
C        itself. The old indices NATOMS+1..NATOMS+4 lie outside the
C        initialised coordinate range (3*NATOMS), so the cell read zeros
C        and SKEW became NaN on every trajectory.
         IF (NGLO.EQ.0 .AND. NATOMB(1).LT.4) THEN
            WRITE(6,*)'ERROR: NSURF=2 (rigid cell) NEEDS >= 4 SURFACE',
     &                ' ATOMS'
            WRITE(6,*)'       NATOMB(1)=',NATOMB(1)
            STOP
         ENDIF
         IF (NGLO.NE.0) THEN
C           Surface-oscillator model (GLO): a single mobile surface atom
C           coupled by springs to a dissipative ghost - no unit cell.
C           Set a degenerate 1x1 reference (SKEW=pi/4 avoids the
C           tan(pi/2) singularity in the NRNDXY aiming formulas).
            BOXLX=0.0D0
            BOXLY=0.0D0
            SKEW=0.7853981633974483D0
            DO I=1,3
               DO J=1,4
                  EDGE(J,I)=0.0D0
               ENDDO
            ENDDO
            GOTO 330
         ENDIF
         NN1=LB(1,1)
         NN2=LB(1,2)
         NN3=LB(1,3)
         NN4=LB(1,4)
         WRITE(6,*)'SURFACE UNIT CELL IS DEFINED BY ATOMS'
     &,NN1,NN2,NN3,NN4
         EDGE(1:3,1)=Q(NN1*3-2:NN1*3)
         EDGE(1:3,2)=Q(NN2*3-2:NN2*3)
         EDGE(1:3,3)=Q(NN3*3-2:NN3*3)
         EDGE(1:3,4)=Q(NN4*3-2:NN4*3)
         BOXLX=DSQRT((Q(3*NN2-2)-Q(3*NN1-2))**2+(Q(3*NN2-1)-Q(3*NN1-1))
     &**2+(Q(3*NN2)-Q(3*NN1))**2)
         BOXLY=DSQRT((Q(3*NN4-2)-Q(3*NN1-2))**2+(Q(3*NN4-1)-Q(3*NN1-1))
     &**2+(Q(3*NN4)-Q(3*NN1))**2)
         DNN24=DSQRT((Q(3*NN4-2)-Q(3*NN2-2))**2+(Q(3*NN4-1)-Q(3*NN2-1))
     &**2+(Q(3*NN4)-Q(3*NN2))**2)
         IF (BOXLX.LE.0.0D0 .OR. BOXLY.LE.0.0D0) THEN
            WRITE(6,*)'ERROR: DEGENERATE RIGID-SURFACE CELL'
            WRITE(6,*)'       ATOMS ',NN1,NN2,NN3,NN4
            WRITE(6,*)'       BOXLX=',BOXLX,' BOXLY=',BOXLY
            STOP
         ENDIF
         CSKEW=(BOXLX**2+BOXLY**2-DNN24**2)/(2D0*BOXLX*BOXLY)
         IF (CSKEW.GE.1D0) CSKEW=1D0
         IF (CSKEW.LE.-1D0) CSKEW=-1D0
         SKEW=DACOS(CSKEW)
      ELSE
         SSKEW=DSIN(SKEW)
         CSKEW=DCOS(SKEW)
         WRITE(6,*)'SURFACE PRAMARY CELL IS DEFINED AS A PERIODIC BOX'
         EDGE(1,1)=0D0
         EDGE(2,1)=0D0
         EDGE(3,1)=0D0
         EDGE(1,2)=BOXLX
         EDGE(2,2)=0D0
         EDGE(3,2)=0D0
         EDGE(1,3)=BOXLX+BOXLY*CSKEW
         EDGE(2,3)=BOXLY*SSKEW
         EDGE(3,3)=0D0
         EDGE(1,4)=BOXLY*CSKEW
         EDGE(2,4)=BOXLY*SSKEW
         EDGE(3,4)=0D0
        do i=1,4
         WRITE(6,'(3F15.8)')EDGE(1:3,i)
        enddo
      ENDIF
  330 CONTINUE
!   END

!-->    ADDED BY BIN AT 9/30/2016
C        NSURF=2  REACTANT A IS A GASEOUS MOLECULE, REACTANT B IS THE
C                    SURFACE, INITIAL CONDITIONS OF A ARE SAMPLED AS
C                    USUAL AND SURFACE B IS MOVABLE WHICH IS SAMPLED VIA
C                    NACT=7
C
C        NSURF=3  REACTANTS INCLUDE A SURFACE, AND INITIAL CONDITIONS
C                    ARE SELECTED BY CHOOSING A RANDOM IMPACT POINT
C                    WITHIN A SURFACE UNIT CELL.  
C
C        NN1, NN2, NN3, AND NN4 DEFINE THE UNIT CELL, NN1 IS THE ORIGIN
C        NN2, NN3, AND NN4 IS IN THE ANTICLOCKWISE SEQUENCE.
C
      IF (NSURF.GT.0) THEN
C
C        SELECT THETA, AND CHI ANGLES 
C        DETERMINATE THE INCIDENT ANGLES
C
         THTA0=THTA
         RAND=RAND0(ISEED)
         CHI0=CHI
         IF (NCHI.NE.1) THEN
            RAND=RAND0(ISEED)
            CHI0=CHI*RAND
         ENDIF
         PI2=2.0D0*DACOS(0.0D0)
C        F24/D4 fix (VENUS05 manual Sec. K, Eq. V.24-V.25): polar incidence
C        angle for the rigid surface. NTHTA>=-1 keeps the beam behaviour
C        THETA=THTA (fixed); NTHTA=-2 samples THETA in [0,THTA=theta_max]
C        with P(theta)~sin(theta), i.e.
C        cos(theta)=1-R*(1-cos(theta_max)).
         IF (NSURF.EQ.2 .AND. NTHTA.EQ.-2) THEN
            RAND=RAND0(ISEED)
            THTA0=DACOS(1.0D0-RAND*(1.0D0-DCOS(THTA0)))
         ENDIF
C        F24/D4 fix (manual Eq. V.27): azimuth PHI2 of the incident
C        velocity. NCHI=1 keeps the fixed CHI; NCHI=2 samples PHI2 uniform
C        in [0,2pi). The legacy NCHI<>1 branch above is kept untouched for
C        NSURF=1.
         IF (NSURF.EQ.2 .AND. NCHI.EQ.2) THEN
            RAND=RAND0(ISEED)
            CHI0=PI2*RAND
         ENDIF
         WRITE(6,13)THTA0/C4,CHI0/C4

C
C        RANDOM SAMPLING IMPACT POINT POSITION (RX0,RY0,RZ0) IN THE UNIT CELL
C
C        FIND THE SKEW ANGLE FOR THE SURFACE FROM THE NN1, NN2 AND NN4 COORDINATES
C
C        THE EDGE OF UNIT CELL WILL DEPEND ON THE SKEW ANGLE
C        
!#MODIFIED by CZZ 2025.9.16
      if(NRNDXY == 1) then
         RAND=RAND0(ISEED)
         YMIN=EDGE(2,1)
         YMAX=EDGE(2,4)
         RY0=YMIN+(YMAX-YMIN)*RAND
         RAND=RAND0(ISEED)
         XMIN=RY0/DTAN(SKEW) !-->DEPENDS ON SKEW ANGLE
         XMAX=XMIN+BOXLX
         RX0=XMIN+(XMAX-XMIN)*RAND
      end if

C        F24/D4 fix (manual Eq. V.19/V.27): disk aiming of fragment A
C        around the origin: impact parameter b fixed (NOB=1) or
C        b=BMAX*sqrt(R) (b^2 uniform on [0,BMAX^2], Eq. V.19), azimuth
C        PHI1 uniform in [0,2pi). Applied for NSURF=2 when BMAX>0;
C        BMAX=0 (default) aims at the origin and consumes no random
C        numbers, so the beam-type stream is unchanged.
         IF (NSURF.EQ.2 .AND. BMAX.GT.0.0D0) THEN
            IF (NOB.EQ.1) THEN
               BDUM=BMAX
            ELSE
               BDUM=BMAX*DSQRT(RAND0(ISEED))
            ENDIF
            RAND=RAND0(ISEED)
            PHI1=PI2*RAND
            RX0=RX0+BDUM*DCOS(PHI1)
            RY0=RY0+BDUM*DSIN(PHI1)
         ENDIF

C        THE SURFACE IS IN XY PLANE
         RZ0=0D0
         WRITE(88,*)RX0,RY0
C
C        DETERMINATE RELATIVE POSITION OF THE MASS CENTER OF FRAGMENT 
C        A (ACX,ACY,ACZ) TO THE AIMING POINT IN THE REFERENCE FRAME
C
         CSTHETA=COS(THTA0)
         SNTHETA=SIN(THTA0)
         CSCHI=COS(CHI0)
         SNCHI=SIN(CHI0)
         DUM2=S*SNTHETA
         ACX=DUM2*CSCHI
         ACY=DUM2*SNCHI
!        For surface scattering: initial height always = S (RMAX)
!        measured perpendicular to surface. Incident angle only
!        affects lateral position (ACX, ACY) and velocity direction.
         ACZ=S
!-->    END OF BIN'S CHANGES

      ENDIF

C     For cluster (NSURF=0): randomize impact site within a box if NRNDXY=1
      IF (NSURF.EQ.0 .AND. NRNDXY.EQ.1) THEN
         RAND=RAND0(ISEED)
         RX0 = RND_BOX * (2.0D0*RAND - 1.0D0)
         RAND=RAND0(ISEED)
         RY0 = RND_BOX * (2.0D0*RAND - 1.0D0)
      ENDIF

C
C     DETERMINATE COORDINATES OF THE MASS CENTER OF FRAGMENT A IN THE
C     REFERENCE FRAME
C
      ACX=ACX+RX0
      ACY=ACY+RY0
      ACZ=ACZ+RZ0

      IF (NZDOWN.EQ.1) THEN
!        Surface oscillator model: offset the gas-phase species.
!        With NATOMB>0 the gas is fragment B; in the GLO layout
!        (NATOMB=0, A=gas + surface + ghost) it is fragment A itself.
         IF (NATOMB(1).GT.0) THEN
            N=NATOMB(1)
            DO I=1,N
               J=3*LB(1,I)
               Q(J)=Q(J)+ACZ
               Q(J-1)=Q(J-1)+ACY
               Q(J-2)=Q(J-2)+ACX
            ENDDO
         ELSE
            N=NATOMA(1)
            DO I=1,N
               J=3*LA(1,I)
               Q(J)=Q(J)+ACZ
               Q(J-1)=Q(J-1)+ACY
               Q(J-2)=Q(J-2)+ACX
            ENDDO
         ENDIF
      ELSE
         N=NATOMA(1)
         DO I=1,N
            J=3*LA(1,I)
            Q(J)=Q(J)+ACZ
            Q(J-1)=Q(J-1)+ACY
            Q(J-2)=Q(J-2)+ACX
         ENDDO
      END IF

C     IF NREL = 0
C     CHOOSE RELATIVE ENERGY FROM BOLTZMANN DISTRIBUTION
C
      IF (NREL.EQ.0) THEN
         DUM=GAMA(2,ISEED)
         SEREL=0.00198717D0*DUM*TRANS
C         WRITE(6,16)SEREL
         SEREL=SEREL*C1
      ENDIF
C
C     ADD RELATIVE TRANSLATIONAL ENERGY
C     RELATIVE VELOCITY VECTOR IS FROM BEAM CENTER TO THE AIMING POINT
C

!!!    modified by Bin at 2016/11/24, add velocity of the gaseous molecule only
C      WT=WTA(1)+WTB(1)
C      SDUM=WTA(1)*WTB(1)/WT
C      DUM=SQRT(2.0D0*SEREL/SDUM)
C      VELA=DUM*WTB(1)/WT
C      VELB=VELA-DUM
      IF(NREL.EQ.1)THEN
        IF (NZDOWN.EQ.1) THEN
C         Gas species mass: fragment B, or fragment A in the GLO
C         NATOMB=0 layout (where WTB=0 would divide by zero).
          WGAS=WTB(1)
          IF (NATOMB(1).LE.0) WGAS=WTA(1)
          DUM=SQRT(2.0D0*SEREL/WGAS)
          VELA=0D0
          VELB=DUM
        ELSE
          DUM=SQRT(2.0D0*SEREL/WTA(1))
          VELA=DUM
          VELB=0D0
        END IF
      ENDIF

!!    Added by zxy at 2017/06/20, sample velocity by density-weighted Gaussian distribution
      IF(NREL.EQ.2)THEN
        SEREL=TRANS*C1
        IF (NZDOWN.EQ.1) THEN
          WGAS=WTB(1)
          IF (NATOMB(1).LE.0) WGAS=WTA(1)
          VELA0=SQRT(2.0D0*SEREL/WGAS)
        ELSE
          VELA0=SQRT(2.0D0*SEREL/WTA(1))
        END IF
        VELA1=0.5D0*VELA0+SQRT(-1.5D0*ALPHAA**2+
     &  0.25D0*VELA0**2)
        VELAMAX=5.0D0*VELA1
        PVMAX=DENS(VELA0,VELA1)
  22    CONTINUE
        RAND=RAND0(ISEED)
        DUM=RAND*VELAMAX
        PVEL=DENS(VELA0,DUM)
        PVEL=PVEL/PVMAX
        RAND=RAND0(ISEED)
        IF(RAND.GT.PVEL) GOTO 22
c        write(11111,*)DUM,PVEL*PVMAX  !!zxywrite
        IF (NZDOWN.EQ.1) THEN
          SEREL=0.5D0*WTB(1)*DUM*DUM
          VELA=0D0
          VELB=DUM
        ELSE
          SEREL=0.5D0*WTA(1)*DUM*DUM
          VELA=DUM
          VELB=0D0
        END IF
      ENDIF

C     
C     UNIT RELATIVE VELOCITY VECTOR IN THE SURFACE REFERENCE FRAME
C
!-->    MODIFIED BY BIN AT 1/23/2014
      IF (NSURF.GT.0) THEN
         ! For surface scattering, ensure atom starts above surface and
         ! velocity always points toward surface (negative z direction)
         UIX=-SNTHETA*CSCHI
         UIY=-SNTHETA*SNCHI
         UIZ=-ABS(CSTHETA)
!-->    END BIN'S CHANGE
      ENDIF
      
C
C     CONTRIBUTION OF RELATIVE VELOCITY TO A AND B MOMENTUM
C
      IF (NZDOWN.EQ.1) THEN
!        Surface oscillator model: incident velocity to the gas species
!        (fragment B, or fragment A in the GLO NATOMB=0 layout).
         IF (NATOMB(1).GT.0) THEN
            N=NATOMB(1)
            DO I=1,N
               J=3*LB(1,I)
               PMOD=W(LB(1,I))*VELB
               P(J)=P(J)+PMOD*UIZ
               P(J-1)=P(J-1)+PMOD*UIY
               P(J-2)=P(J-2)+PMOD*UIX
            ENDDO
         ELSE
            N=NATOMA(1)
            DO I=1,N
               J=3*LA(1,I)
               PMOD=W(LA(1,I))*VELB
               P(J)=P(J)+PMOD*UIZ
               P(J-1)=P(J-1)+PMOD*UIY
               P(J-2)=P(J-2)+PMOD*UIX
            ENDDO
         ENDIF
      ELSE
         N=NATOMA(1)
         DO I=1,N
            J=3*LA(1,I)
            PMOD=W(LA(1,I))*VELA
            P(J)=P(J)+PMOD*UIZ
            P(J-1)=P(J-1)+PMOD*UIY
            P(J-2)=P(J-2)+PMOD*UIX
         ENDDO
      END IF
!-->    MODIFIED BY BIN AT 4/24/2014
!-->    NO INITIAL MOMENTUM ON SURFACE, IT IS FIXED

!-->    END BIN'S CHANGE
      CALL DVDQ_1
      CALL ENERGY_1
C
C     SAVE INITIAL RELATIVE VELOCITY AND ORBITAL ANGULAR MOMENTUM
C
      VI(1)=DUM*UIX
      VI(2)=DUM*UIY
      VI(3)=DUM*UIZ
      VI(4)=DUM
CC
CC     SAVE INITIAL ORBITAL ANGULAR MOMENTUM (NOTE THAT MASS CENTER OF B
CC     IS THE ORIGIN OF PROGRAM'S COORDINAT SYSTEM): 
CC               ->      ->    ->           ->
CC               OAM = (QCMA -QCMB) X (SDUM*VI)
CC
C      OAMI(1)=SDUM*DUM*(ACY*UIZ-ACZ*UIY)/C7
C      OAMI(2)=SDUM*DUM*(ACZ*UIX-ACX*UIZ)/C7
C      OAMI(3)=SDUM*DUM*(ACX*UIY-ACY*UIX)/C7
C      OAMI(4)=SQRT(OAMI(1)*OAMI(1)+OAMI(2)*OAMI(2)+
C     *       OAMI(3)*OAMI(3))
C      WRITE(6,18)OAMI(1),OAMI(2),OAMI(3)
CC
      RETURN
      END

C    MODIFIED FROM "COMPUTER SIMULATION OF LIQUIDS" BY ALLEN AND
C    TILDESLEY

      SUBROUTINE PBC (NATOMS)
      use venus_data, only: NATOMS_H => NATOMS
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C    *******************************************************************
C    ** PERIODIC BOUNDARY CONDITIONS FOR AN ARBYRARY RHOMBOID BOX.    **
C    **                                                               **
C    ** PERIODIC CORRECTIONS ARE APPLIED IN TWO DIMENSIONS X, Y.      **
C    ** IN MOST APPLICATIONS THE MOLECULES WILL BE CONFINED IN THE    **
C    ** Z DIRECTION BY REAL WALLS RATHER THAN BY PERIODIC BOUNDARIES. **
C    ** THE BOX IS CENTRED AT THE ORIGIN. THE X AXIS LIES ALONG THE   **
C    ** SIDE OF THE CELL                                              **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    RT3                SQRT(3.0) TO MACHINE ACCURACY      **
C    ** REAL    RT32               SQRT(3.0)/2.0                      **
C    ** REAL    RRT3               1.0/SQRT(3.0)                      **
C    ** REAL    RRT32              2.0/SQRT(3.0)                      **
C    *******************************************************************


C    *******************************************************************

        CSKEW=DCOS(SKEW)
        SSKEW=DSIN(SKEW)
        TSKEW=DTAN(SKEW)
        DO I=1,NATOMS
           J1=I*3-2
           J2=I*3-1
           J3=I*3
           UX1=Q(J1)-Q(J2)/TSKEW
           UY1=Q(J2)/SSKEW
           UX2=BOXLX*ANINT(UX1/BOXLX)
           UY2=BOXLY*ANINT(UY1/BOXLY)
           Q(J1) = Q(J1) - UX2 - UY2*SSKEW
           Q(J2) = Q(J2) - UY2*CSKEW
        ENDDO

        RETURN
        END

!!      Added by zxy, 20170620       
        FUNCTION DENS(X0,DUM)
        use venus_data
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        DENS=DUM*DUM*DUM*EXP(-(((DUM-X0)/ALPHAA)**2))

        RETURN
        END


