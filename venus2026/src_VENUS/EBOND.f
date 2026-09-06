      SUBROUTINE EBOND(H,T,R,I)
      use venus_params
      use venus_data, only: D, B, RMZ, W, N2J, N2K, Q, P, C1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         CALCULATE LOCAL MODE (MORSE OSCILLATOR) ENERGIES 
C
C
      DE=D(I)
      BETA=B(I)
      R0=RMZ(I)
      WA=W(N2J(I))
      WB=W(N2K(I))
      WT=WA+WB
      WJ=WB/WT
      WI=WA/WT
      UMASS=WA*WB/WT
      IZ=3*N2J(I)
      JZ=3*N2K(I)
      IY=IZ-1
      IX=IY-1
      JY=JZ-1
      JX=JY-1
      R=SQRT((Q(IX)-Q(JX))**2 + (Q(IY)-Q(JY))**2 + (Q(IZ)-Q(JZ))**2)
      T=((WJ*P(IX)-WI*P(JX))**2 + (WJ*P(IY)-WI*P(JY))**2 +
     *                            (WJ*P(IZ)-WI*P(JZ))**2)/2.0D0/UMASS
      V=EXP(-BETA*(R-R0))
      V=DE*(V**2-2*V+1.0D0)
      T=T/C1
      V=V/C1
      H= V + T
      RETURN
      END
