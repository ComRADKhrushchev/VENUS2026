subroutine HOMOQP(R, PR, AL, AM, RMASS, A)
  use venus_params
  use venus_data
  implicit none
  real(8),  intent(in)    :: R, PR, AL, RMASS
  real(8),  intent(inout) :: AM(4)
  real(8),  intent(out)   :: A(3)
  real(8)  :: QCM(3), VCM(3), VREL, WT, VELA, VELB, RAND, DUM, EROT
  integer  :: I, J3, J2, J1, J, K, ISEED
  real(8), external :: RAND0
!
!         INITIALIZE COORDINATES AND MOMENTA FOR A DIATOM
!
!--> added by Bin, 2014/6/18
!
!         SET CARTESIAN COORDINATES AND MOMENTA
!
  do I = 1, 2
     J3 = 3*LL(I)
     J2 = J3-1
     J1 = J2-1
     Q(J1) = 0.0D0
     Q(J2) = 0.0D0
     Q(J3) = 0.0D0
     P(J2) = 0.0D0
     P(J3) = 0.0D0
  end do
  Q(J1) = R
  VREL = PR/RMASS
  WT = W(LBOND(1))+W(LBOND(2))
  VELA = VREL*W(LBOND(2))/WT
  VELB = VELA-VREL
  J1 = 3*LBOND(1)-2
  P(J1) = W(LBOND(1))*VELA
  J1 = 3*LBOND(2)-2
  P(J1) = W(LBOND(2))*VELB
!
!         SET INERTIA ARRAYS.  CHOOSE Y AND Z ANGULAR MOMENTUM COMPONENTS
!
  A(1) = 1.0D+20
  A(2) = RMASS*R**2
  A(3) = A(2)
  RAND = RAND0(ISEED)
  DUM = TWOPI*RAND
  AM(1) = 0.0D0
  AM(2) = AL*sin(DUM)
  AM(3) = AL*cos(DUM)
!--> added by Bin, 2014/6/25
  if (NTHTA >= 0) then
     AM(3) = AL
     AM(2) = 0.0D0
  end if
!--> end
!
!         CALCULATE CENTER OF MASS COORDIANTES QQ AND MOMENTA PP
!
  call CENMAS(WT, QCM, VCM, 2)
!
!         MOVE PP ARRAY TO P ARRAY AND QQ ARRAY TO Q ARRAY
!
  do I = 1, 2
     J = 3*LL(I)+1
     do K = 1, 3
        Q(J-K) = QQ(J-K)
        P(J-K) = PP(J-K)
     end do
  end do
!
!         ADD THE ANGULAR MOMENTUM VECTOR.  CALCULATE THE REQUIRED
!         ANGULAR VELOCITY AND ADD IT TO THE DIATOM, WHICH LIES
!         ALONG THE X-AXIS.
!
  WX = 0.0D0
  WY = -AM(2)/A(2)
  WZ = -AM(3)/A(3)
  call ANGVEL(2)

!-->    modified by Bin, 1/24/2014
  if (NATOMS == 2 .and. NSURF == 0) return
!-->    end

!-->    modified by Bin, 6/24/2014
  if (NTHTA >= 0) then
     if (NTHTA > JA) then
        stop 'NTHTA (Mj) CAN NOT EXCEED JA'
     else
!
!         ROTATE A DIATOM (N=2) TO MODEL A SPECIFIC (J,M) STATE BY A VECTOR MODEL
!
        call ROTATEJM(NTHTA, JA, 2)
     end if
  else
!
!         RANDOMLY ROTATE THE DIATOM ABOUT ITS CENTER OF MASS BY
!         EULER'S ANGLES.  CENTER OF MASS COORDINATES QQ AND MOMENTA
!         PP ARE PASSED FROM SUBROUTINES CENMAS AND ANGVEL THROUGH
!         COMMON BLOCK WASTE.
!
!     modifed by bin, 2020/4/3
     if (NTHTA == -1) then
        call ROTATE(2)
     else if (NTHTA == -2 .or. NTHTA == -3) then
        call ROTATENO(2, NTHTA)
     else
        stop 'NTHTA IS NOT DEFINED'
     end if
!     end 2020/4/3
  end if

!--> end
!
!         CALCULATE ANGULAR MOMENTUM AND COMPONENTS.
!
  call ROTN(AM, EROT, 2)
end subroutine HOMOQP
