      SUBROUTINE FMTRX(NATOM,NDIS,I3NS)
      use venus_params
      use venus_data, only: Q, PDOT, FCOEF, P, QDOT, W,
     *                     EIG, NTZ, NSELT, DIM, DG,
     *                     A => A_MPATH, DA => DA_MPATH, B => B_MPATH,
     *                     GAUHES, HINC, NPTS, IARB, NACTA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical ngflag
C
C         EVALUATE THE FORCE CONSTANT MATRIX BY DIFFERENCING
C         THE GRADIENT OF THE POTENTIAL ENERGY FUNCTION.
C
      DIMENSION DIST(2),GZ(NDA3),FA(NDA3yf,NDA3yf)
      DATA PT5/0.5D0/
  100 FORMAT(10X,'***** CALCULATION OF FORCE CONSTANT  *****')
  102 FORMAT(///)
  103 FORMAT(15X,' FORCE CONSTANT MATRIX'/)
  104 FORMAT(/10X,'SYMMETRIZED CARTESIAN FORCE CONSTANT MATRIX'/)
c
c   This is to determine whether the gaussian will calculate the force
c   constant matrix or use numerical second derivative
c
      if(nmo.lt.0.or.(nmo.gt.0.and.ndis.eq.0.and.natom.eq.nmo))then
        ngflag=.true.
      else
        ngflag=.false.
      endif
C
C         SET PRINT FILE FOR EIGOUT
C
      IP=77
      DIST(1)=HINC
      DIST(2)=-HINC
      DO I=1,I3NS
         GZ(I)=-PDOT(I+3*NDIS)
      ENDDO
C
C         INITIALIZE SOME ARRAYS
C
      WRITE(ip,100)
C
C              GETTING FORCE MATRIX FROM GAUSSIAN 
c  For NMO < 0, this is entirely from Gaussian Hessian matrix.
c  Otherwise use the difference method to get the second derivative.
C
      IF (ngflag) THEN
         WRITE(6,*)'            GETTING FORCE MATRIX FROM GAUSSIAN'
         IARB=1                 ! FLAG FOR FREQUENCY CALCULATION
         CALL DVDQ_1
         IARB=0
         I3A=3*NDIS
         KI=I3A*(I3A+1)/2
         DO K=1,I3NS
            KI=KI+I3A
            DO I=1,K
              KI=KI+1
              FA(K,I)=GAUHES(KI)
              FA(I,K)=FA(K,I)
            ENDDO
         ENDDO
      ELSE
C
C             GETTING FORCE MATRIX FROM ANALYTICAL POTENTIALS
C
         IF (NPTS.GT.2.or.npts.le.0) NPTS=2
         DO I=1,I3NS
            DO J=1,NPTS
               Q(I+3*NDIS)=Q(I+3*NDIS)+DIST(J)
C
C              DISPLACE COORDINATE AND CALCULATE GRADIENT INTO DG
C
               CALL DVDQ_1
               Q(I+3*NDIS)=Q(I+3*NDIS)-DIST(J)

               DO K=1,I3NS
                  DG(K,J)=-PDOT(K+3*NDIS)
               ENDDO
            ENDDO
C
C              TWO POINT DIFFERENCE FORMULA
C
            IF (NPTS.NE.1) THEN
               DO K=1,I3NS
                  FA(K,I)=(DG(K,1)-DG(K,2))*PT5/HINC
               ENDDO
            ELSE
C
C              SIMPLE ONE POINT FORMULA
C
               DO K=1,I3NS
                  FA(K,I)=(DG(K,1)-GZ(K))/HINC
               ENDDO
            ENDIF
         ENDDO
C
      ENDIF
C

      DO I=1,I3NS
         EIG(I)=0.0D0
      ENDDO
      WRITE(ip,102)
      WRITE(ip,103)
      CALL EIGOUT(FA,I3NS,IP)
      DO I=1,I3NS
         DO J=1,I
            IF (ngflag) THEN
               DU=FA(I,J)
               DL=0.0D0
            ELSE
               DU=PT5*(FA(I,J)+FA(J,I))
               DL=PT5*(FA(I,J)-FA(J,I))
            ENDIF
            A(I,J)=DU
            A(J,I)=DL
         ENDDO
      ENDDO
      DO I=1,I3NS
         A(I,I)=FA(I,I)
      ENDDO
      WRITE(ip,102)
      WRITE(ip,104)
      CALL EIGOUT(A,I3NS,IP)
C
C         CALCULATE ARRAY DIM(150) USED FOR MASS-WEIGHTING
C         CONVERT TO MASS WEIGHTED COORDINATES AND CALCULATE THE NORMAL
C         MODES AND THE SPECTROSCOPIC FREQUENCIES.
C
      K=0
      DO I=1,NATOM
         DO J=1,3
            K=K+1
            DIM(K)=1.D0/SQRT(W(I+NDIS))
         ENDDO
      ENDDO
      CALL FGMTRX(I3NS)
      RETURN
      END
