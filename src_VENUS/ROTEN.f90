!--> adapted by Bin, 2016/7/30
subroutine ROTEN(AM, AI, TROT, EROT, NROT, NLIN, JROT, KROT)
!    end
  use venus_params
  use venus_data, only: AI_H => AI, TEMP_H => TEMP
  use venus_data
  implicit none
  real(8),  intent(inout) :: AM(4), EROT
  real(8),  intent(in)    :: AI(3), TROT
  integer,  intent(in)    :: NROT, NLIN, JROT
  integer,  intent(inout) :: KROT
  real(8)  :: EROTT, DUM1, DUM2, DUM, AMJ, ALJ, RAND, ALZMAX, AL
  real(8)  :: TEMP
  integer  :: I, JDUM, ISEED
  real(8), external :: RAND0
!
!         SELECT ANGULAR MOMENTUM AND ROTATIONAL ENERGY
!
!
!         THE ARRAY SIZES OF THE RELEVANT CALLING ARGUEMENTS ARE
!         DEFINED IN COMMON BLOCKS IN SUBROUTINE SELECT.
!
!                   ROTEN:             SELECT:
!
!                   AM(4)              AMA(4)  AMB(4)
!                   AI(3)              AI(3)   BI(3)
!
!
!         NROT=0 , CHOOSE ROTATIONAL ENERGY FROM A THERMAL
!                  DISTRIBUTION BY ASSUMING A SYMMETRIC TOP.
!                  FARADAY DISCUSSIONS 55, 93(1973).
!                  IF Z IS THE SYMMETRY AXIS, THEN IX=IY.
!                  IF X IS THE SYMMETRY AXIS, THEN IY=IZ.
!         NROT=1 , ROTATIONAL ENERGY ABOUT EACH AXIS EQUALS RT/2.
!--> added by Bin, 2014/6/18
!         NROT=2 , SAMPLE ROTATIONAL ENERGY FOR A SYMMETRIC TOP (J,K)
!                  IF Z IS THE SYMMETRY AXIS, THEN IX=IY.
!                  IF X IS THE SYMMETRY AXIS, THEN IY=IZ.
!--> end
!
!         NLIN=0 , MOLECULE IS NONLINEAR
!         NLIN=1 , MOLECULE IS LINEAR
!         NOTE: , A LINEAR MOLECULE MUST LIE ALONG THE X-AXIS.
!                 (THE PROGRAM DOES IT AUTOMATICALLY NOW)
!         NOTE: , LZMAX**2 = 10*2*AI(1)*C5*TROT
!
  EROT = 0.0D0
  AM(1) = 0.0D0
  TEMP = C5*TROT
!--> added by Bin, 2014/6/18
  if (NROT == 2) then
     if (KROT > JROT) then
        stop 'KROT MUST BE SMALLER THAN JROT'
     end if
     if (JROT <= 0 .and. KROT <= 0) then
        AM = 0D0
        EROT = 0D0
        return
     end if
     if (NLIN == 1) KROT = 0

     DUM1 = abs(AI(1)-AI(2))
     DUM2 = abs(AI(2)-AI(3))
     if (DUM1 <= DUM2) then
        JDUM = 3           !-->Primary axis is z axis
        DUM = dsqrt(AI(1)*AI(2))
     else
        JDUM = 1           !-->Primary axis is x axis
        DUM = dsqrt(AI(2)*AI(3))
     end if
!-->total rotational energy, erot=B*J*(J+1)-(A-B)*K^2
     EROT = JROT*(JROT+1)/2D0/DUM+KROT**2/2d0/AI(JDUM)-KROT**2/2d0/DUM
     EROTT = EROT*C7*C7/C1
!-->total rotational angular momentum J^=dsqrt(J(J+1))*HBAR
     AMJ = sqrt(dble(JROT*(JROT+1)))
     ALJ = AMJ*C7
!-->rotational angular momemtum about z axis
     AM(JDUM) = ALJ*KROT/AMJ
     RAND = RAND0(ISEED)
     if (RAND > 0.5D0) AM(JDUM) = -AM(JDUM)

!-->rotational angular momemtum about x and y axis
     if (JDUM == 1) then
        DUM = sqrt(ALJ**2-AM(1)**2)
        RAND = RAND0(ISEED)
!-->randomly sample the angular momentum on y and z axis
        AM(2) = DUM*dsin(TWOPI*RAND)
        AM(3) = DUM*dcos(TWOPI*RAND)
     else if (JDUM == 3) then
        DUM = sqrt(ALJ**2-AM(3)**2)
        RAND = RAND0(ISEED)
!-->randomly sample the angular momentum on x and y axis
        AM(1) = DUM*dsin(TWOPI*RAND)
        AM(2) = DUM*dcos(TWOPI*RAND)
     end if
     if (NLIN == 1) then
        EROT = EROTT
     else
        EROT = (AM(1)**2/AI(1)+AM(2)**2/AI(2)+AM(3)**2/AI(3))/2.0D0/C1
     end if
     write(*,*) 'ROTATIONAL ENERGY BY J AND K ', EROTT
     write(*,*) 'ROTATIONAL ENERGY BY AM1, AM2, AND AM3 ', EROT
     if (abs(EROTT-EROT) >= 1d-5) then
        write(*,*) 'ERROR IN ROTEN.f'
        stop
     end if
!--> end
  else if (NROT == 1) then
     do I = NLIN+1, 3
        AM(I) = sqrt(AI(I)*TEMP)
        RAND = RAND0(ISEED)
        if (RAND < 0.5D0) AM(I) = -AM(I)
     end do
     EROT = dble(3-NLIN)*TEMP/2.0D0/C1
  else if (NROT == 0) then
     if (NLIN == 0) then
        DUM1 = abs(AI(1)-AI(2))
        DUM2 = abs(AI(2)-AI(3))
        if (DUM1 <= DUM2) then
           JDUM = 3           !-->Primary axis is z axis
        else
           JDUM = 1           !-->Primary axis is x axis
        end if
        ALZMAX = sqrt(20.0D0*AI(JDUM)*TEMP)
        do
           RAND = RAND0(ISEED)
           AM(JDUM) = RAND*ALZMAX
           DUM = exp(-AM(JDUM)**2/2.0D0/AI(JDUM)/TEMP)
           RAND = RAND0(ISEED)
           if (RAND <= DUM) exit
        end do
        RAND = RAND0(ISEED)
        if (RAND > 0.5D0) AM(JDUM) = -AM(JDUM)
        EROT = AM(JDUM)**2/AI(JDUM)
     end if
     if (NLIN == 1 .or. JDUM == 1) then
        RAND = RAND0(ISEED)
        DUM = sqrt(AI(2)*AI(3))
        AL = sqrt(AM(1)**2-2.0D0*DUM*TEMP*log(1.0D0-RAND))
        DUM = sqrt(AL**2-AM(1)**2)
        RAND = RAND0(ISEED)
        AM(2) = DUM*sin(TWOPI*RAND)
        AM(3) = DUM*cos(TWOPI*RAND)
        EROT = (AM(2)**2/AI(2)+AM(3)**2/AI(3)+EROT)/2.0D0/C1
     else
        RAND = RAND0(ISEED)
        DUM = sqrt(AI(1)*AI(2))
        AL = sqrt(AM(3)**2-2.0D0*DUM*TEMP*log(1.0D0-RAND))
        DUM = sqrt(AL**2-AM(3)**2)
        RAND = RAND0(ISEED)
        AM(1) = DUM*sin(TWOPI*RAND)
        AM(2) = DUM*cos(TWOPI*RAND)
        EROT = (AM(1)**2/AI(1)+AM(2)**2/AI(2)+EROT)/2.0D0/C1
     end if
  end if
!
end subroutine ROTEN
