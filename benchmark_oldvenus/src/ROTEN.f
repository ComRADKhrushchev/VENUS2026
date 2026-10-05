!--> adapted by Bin, 2016/7/30
      SUBROUTINE ROTEN(AM,AI,TROT,EROT,NROT,NLIN,JROT,KROT)
!    end
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C         SELECT ANGULAR MOMENTUM AND ROTATIONAL ENERGY
C
      COMMON/CONSTN/C1,C2,C3,C4,C5,C6,C7,PI,HALFPI,TWOPI
      COMMON/SELTB/QZ(NDA3),NSELT,NSFLAG,NACTA,NACTB,NLINA,NLINB,NSURF
 
C
C         THE ARRAY SIZES OF THE RELEVANT CALLING ARGUEMENTS ARE
C         DEFINED IN COMMON BLOCKS IN SUBROUTINE SELECT.
C
C                   ROTEN:             SELECT:
C
C                   AM(4)              AMA(4)  AMB(4)
C                   AI(3)              AI(3)   BI(3)
C
      DIMENSION AM(4),AI(3)
C
C         NROT=0 , CHOOSE ROTATIONAL ENERGY FROM A THERMAL
C                  DISTRIBUTION BY ASSUMING A SYMMETRIC TOP.
C                  FARADAY DISCUSSIONS 55, 93(1973).
C                  IF Z IS THE SYMMETRY AXIS, THEN IX=IY.
C                  IF X IS THE SYMMETRY AXIS, THEN IY=IZ.
C         NROT=1 , ROTATIONAL ENERGY ABOUT EACH AXIS EQUALS RT/2.
!--> added by Bin, 2014/6/18
C         NROT=2 , SAMPLE ROTATIONAL ENERGY FOR A SYMMETRIC TOP (J,K)
C                  IF Z IS THE SYMMETRY AXIS, THEN IX=IY.
C                  IF X IS THE SYMMETRY AXIS, THEN IY=IZ.
!--> end
C
C         NLIN=0 , MOLECULE IS NONLINEAR
C         NLIN=1 , MOLECULE IS LINEAR
C         NOTE: , A LINEAR MOLECULE MUST LIE ALONG THE X-AXIS.
C                 (THE PROGRAM DOES IT AUTOMATICALLY NOW)
C         NOTE: , LZMAX**2 = 10*2*AI(1)*C5*TROT
C
      EROT=0.0D0
      AM(1)=0.0D0
      TEMP=C5*TROT
!--> added by Bin, 2014/6/18
      IF (NROT.EQ.2) THEN
              IF(KROT.GT.JROT) THEN
                      STOP 'KROT MUST BE SMALLER THAN JROT'
              ENDIF
              IF(JROT.LE.0.AND.KROT.LE.0) THEN
                      AM=0D0
                      EROT=0D0
                      RETURN
              ENDIF
              IF(NLIN.EQ.1) THEN
                      STOP 'NROT=2 IS ONLY FOR SYMMETRIC TOP'
              ENDIF
            DUM1=ABS(AI(1)-AI(2))
            DUM2=ABS(AI(2)-AI(3))
            IF (DUM1.LE.DUM2) THEN
               JDUM=3           !-->Primary axis is z axis
            DUM=DSQRT(AI(1)*AI(2))
            ELSE
               JDUM=1           !-->Primary axis is x axis
            DUM=DSQRT(AI(2)*AI(3))
            ENDIF
!-->total rotational energy, erot=B*J*(J+1)-(A-B)*K^2
        EROT=JROT*(JROT+1)/2D0/DUM
     &+KROT**2/2d0/AI(JDUM)-KROT**2/2d0/DUM
        EROTT=EROT*C7*C7/C1
!-->total rotational angular momentum J^=dsqrt(J(J+1))*HBAR
        AMJ=SQRT(DBLE(JROT*(JROT+1)))
        ALJ=AMJ*C7
!-->rotational angular momemtum about z axis
        AM(JDUM)=ALJ*KROT/AMJ
        RAND=RAND0(ISEED)
        IF (RAND.GT.0.5D0) AM(JDUM)=-AM(JDUM)
 
!-->rotational angular momemtum about x and y axis
         IF (JDUM.EQ.1) THEN
            DUM=SQRT(ALJ**2-AM(1)**2)
            RAND=RAND0(ISEED)
!-->randomly sample the angular momentum on x and y axis
            AM(2)=DUM*DSIN(TWOPI*RAND)
            AM(3)=DUM*DCOS(TWOPI*RAND)
         ELSEIF (JDUM.EQ.3) THEN
            DUM=SQRT(ALJ**2-AM(3)**2)
            RAND=RAND0(ISEED)
!-->randomly sample the angular momentum on x and y axis
            AM(1)=DUM*DSIN(TWOPI*RAND)
            AM(2)=DUM*DCOS(TWOPI*RAND)
         ENDIF
         EROT=(AM(1)**2/AI(1)+AM(2)**2/AI(2)+AM(3)**2/AI(3))/2.0D0/C1
         WRITE(*,*)'ROTATIONAL ENERGY BY J AND K ', EROTT
         WRITE(*,*)'ROTATIONAL ENERGY BY AM1, AM2, AND AM3 ', EROT
         IF (ABS(EROTT-EROT).GE.1d-5) THEN
                 WRITE(*,*)'ERROR IN ROTEN.f'
                 STOP
         ENDIF
!--> end
      ELSEIF (NROT.EQ.1) THEN
         DO I=NLIN+1,3
            AM(I)=SQRT(AI(I)*TEMP)
            RAND=RAND0(ISEED)
            IF (RAND.LT.0.5D0) AM(I)=-AM(I)
         ENDDO
         EROT=DBLE(3-NLIN)*TEMP/2.0D0/C1
      ELSEIF (NROT.EQ.0) THEN
         IF (NLIN.EQ.0) THEN
            DUM1=ABS(AI(1)-AI(2))
            DUM2=ABS(AI(2)-AI(3))
            IF (DUM1.LE.DUM2) THEN
               JDUM=3           !-->Primary axis is z axis
            ELSE
               JDUM=1           !-->Primary axis is x axis
            ENDIF
            ALZMAX=SQRT(20.0D0*AI(JDUM)*TEMP)
    2       CONTINUE
            RAND=RAND0(ISEED)
            AM(JDUM)=RAND*ALZMAX
            DUM=EXP(-AM(JDUM)**2/2.0D0/AI(JDUM)/TEMP)
            RAND=RAND0(ISEED)
            IF (RAND.GT.DUM) GOTO 2
            RAND=RAND0(ISEED)
            IF (RAND.GT.0.5D0) AM(JDUM)=-AM(JDUM)
            EROT=AM(JDUM)**2/AI(JDUM)
         ENDIF
         IF (NLIN.EQ.1.OR.JDUM.EQ.1) THEN
            RAND=RAND0(ISEED)
            DUM=SQRT(AI(2)*AI(3))
            AL=SQRT(AM(1)**2-2.0D0*DUM*TEMP*LOG(1.0D0-RAND))
            DUM=SQRT(AL**2-AM(1)**2)
            RAND=RAND0(ISEED)
            AM(2)=DUM*SIN(TWOPI*RAND)
            AM(3)=DUM*COS(TWOPI*RAND)
            EROT=(AM(2)**2/AI(2)+AM(3)**2/AI(3)+EROT)/2.0D0/C1
         ELSE
            RAND=RAND0(ISEED)
            DUM=SQRT(AI(1)*AI(2))
            AL=SQRT(AM(3)**2-2.0D0*DUM*TEMP*LOG(1.0D0-RAND))
            DUM=SQRT(AL**2-AM(3)**2)
            RAND=RAND0(ISEED)
            AM(1)=DUM*SIN(TWOPI*RAND)
            AM(2)=DUM*COS(TWOPI*RAND)
            EROT=(AM(1)**2/AI(1)+AM(2)**2/AI(2)+EROT)/2.0D0/C1
         ENDIF
      ENDIF
      RETURN
C
      END
