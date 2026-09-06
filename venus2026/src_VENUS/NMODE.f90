subroutine NMODE(NATOM, NDIS)
  use venus_params
  use venus_data
  implicit none
  integer, intent(in) :: NATOM, NDIS
  real(8) :: X(NDA3), GZS(NDA3)
  integer :: I
!
!         DRIVER FOR NORMAL MODE ANALYSIS
!         DOES NOT ALTER EITHER COORDINATES NOR THE ENERGY GRADIENT.
!
!
!         WRITE RELEVANT INFORMATION IN CHECKPOINT FILE
!         (NORMAL MODE ANALYSIS OR REACTION PATH FOLLOWING)
!
!   removed this part by bin, 2016/7/30
!      IF (NSELT.LT.0) THEN
!         OPEN(50,FORM='UNFORMATTED')
!         REWIND(50)
!         WRITE(50)Q,P,QDOT,PDOT,TABLE,VRELO,RANLST,GTEMP,NFLAG,
!     *            ISEED0,ISEED3,NX,NC,NTZ,INTST,NAST,IBFCTR,
!     *            VI,OAMI,AMAI,AMBI,ETAI,ERAI,ETBI,ERBI
!         CLOSE(50)
!      ENDIF
!   end
!
!         SAVE COORDINATES AND GRADIENT
!
  I3N = 3*NATOM
  do I = 1, I3N
     X(I) = Q(I+3*NDIS)
     GZS(I) = PDOT(I+3*NDIS)
  end do
  !write(6,"(A10I5)") 'NDIS:',NDIS    NDIS=0
!
  call FMTRX(NATOM, NDIS, I3N)
!
  if (NSELT == -1) return
!
!         RESTORE COORDINATES AND GRADIENT
!
  do I = 1, I3N
     Q(I+3*NDIS) = X(I)
     PDOT(I+3*NDIS) = GZS(I)
  end do
end subroutine NMODE
