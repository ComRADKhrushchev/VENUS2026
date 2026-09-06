      SUBROUTINE ROTATEY(CSTHTA,N)
      use venus_params
      use venus_data
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         RANDOMLY ROTATE A MOLECULE ABOUT ITS CENTER OF MASS BY
C         EULER'S ANGLES.
C
C
C         CENTER OF MASS COORDINATES(QQ) AND MOMENTA(PP) ARE PASSED FROM
C         SUBROUTINES CENMAS AND ANGVEL THROUGH COMMON BLOCK WASTE.
C 
      RTHTA=ACOS(CSTHTA)    
      SNTHTA=SIN(RTHTA)
      DO I=1,N
      J=3*LL(I)
         Q(J-2)=QQ(J-2)*CSTHTA+QQ(J)*SNTHTA
         Q(J-1)=QQ(J-1)
         Q(J)=-QQ(J-2)*SNTHTA+QQ(J)*CSTHTA
         QQ(J-2)=Q(J-2)
         QQ(J-1)=Q(J-1)
         QQ(J)=Q(J)
         P(J-2)=PP(J-2)*CSTHTA+PP(J)*SNTHTA
         P(J-1)=PP(J-1)
         P(J)=-PP(J-2)*SNTHTA+PP(J)*CSTHTA
         PP(J-2)=P(J-2)
         PP(J-1)=P(J-1)
         PP(J)=P(J)
      ENDDO
      RETURN
      END
