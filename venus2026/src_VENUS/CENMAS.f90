subroutine CENMAS(WT, QCM, VCM, N)
  use venus_params
  use venus_data
  implicit none
  real(8),  intent(in)  :: WT
  real(8),  intent(out) :: QCM(3), VCM(3)
  integer,  intent(in)  :: N
  integer :: I, J, J3, J2, J1
!
!         CALCULATE THE CENTER OF MASS MOMENTA AND COORDINATES
!
!
!         CENTER OF MASS COORDINATES AND MOMENTA ARE STORED IN
!         ARRAYS QQ AND PP
!
  do I = 1, 3
     VCM(I) = 0.0D0
     QCM(I) = 0.0D0
  end do
!
  do I = 1, N
     J = LL(I)
     J3 = 3*J
     J2 = J3-1
     J1 = J2-1
     VCM(1) = VCM(1)+P(J1)
     VCM(2) = VCM(2)+P(J2)
     VCM(3) = VCM(3)+P(J3)
     QCM(1) = QCM(1)+W(J)*Q(J1)
     QCM(2) = QCM(2)+W(J)*Q(J2)
     QCM(3) = QCM(3)+W(J)*Q(J3)
  end do
!
  do I = 1, 3
     VCM(I) = VCM(I)/WT
     QCM(I) = QCM(I)/WT
  end do
!
  do I = 1, N
     J = LL(I)
     J3 = 3*J
     J2 = J3-1
     J1 = J2-1
     PP(J1) = P(J1)-W(J)*VCM(1)
     PP(J2) = P(J2)-W(J)*VCM(2)
     PP(J3) = P(J3)-W(J)*VCM(3)
     QQ(J1) = Q(J1)-QCM(1)
     QQ(J2) = Q(J2)-QCM(2)
     QQ(J3) = Q(J3)-QCM(3)
  end do
end subroutine CENMAS
