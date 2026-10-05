!
!    VELOCITY VERLET (1) AND BEEMAN'S ALGORITHM (2) INTEGRATOR
!    CODED BY KI SONG, 5/23/05
!    CODED BY BIN JIANG, 2/15/2017
!
!    VELOCITY VERLET
!      X(T+DT)=X(T)+V(T)DT+F(T)/2M*DT**2
!      V(T+DT)=V(T)+(F(T+DT)+FT)/2M*DT
!
!    BEEMAN
!      X(T+DT)=X(T)-V(T)*DT+(4*F(T)-F(T-DT))*DT**2/M/6
!      V(T)=(X(T+DT)-X(T)/DT+(2*F(T+DT)+F(T))*DT/6
!
!    F(T)=-DV/DQ=PDOT
!    VENUS PUT PDOT LIKE THIS. SYMPLE.F USES THIS EXPRESSION.
!
      SUBROUTINE VERLET(NVERLET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON/QPDOT/Q(NDA3),PDOT(NDA3),FCOEF(NDA3,NDA3)
      COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)
      COMMON/GLOP/WS1(3),WG1(3),WS2(3),WG2(3),WGS1(3),WGS2(3),WEFF(3),GSW,FCG,COEFA,COEFB,GN(NDA3)
      COMMON/FORCES/NATOMS,I3N,NFC,NGLO
      COMMON/FRAGB/WTA(NDP),WTB(NDP),LA(NDP,NDA),LB(NDP,NDA),QZA(NDP,NDA3),QZB(NDP,NDA3),NATOMA(NDP),NATOMB(NDP)
      COMMON/INTEGR/ATIME,NI,NID
      COMMON/PRLIST/T,V,H,TIME,NTZ,NT,ISEED0(8),NC,NX
      COMMON/VRSCAL/THERMOTEMP,NSEL,NSCALE,NEQUAL,NRGD

      COMMON/DIABATS/V_Coup,t_normal(NDA3),itranstype,iPES,nPES,IFmin   !耦合 跃迁方向向量 跃迁类型(0=LZ,1=ZN) 势能面指标 是否为极小值(bool)
	  COMMON/T13/Qt3(480),PDOT1t3(480),PDOT2t3(480),Qt1(480),PDOT1t1(480),PDOT2t1(480),E1t1,E2t1,E1t3,E2t3  !存储前后两步的绝热力，透热力，能量
	  
	!  		COMMON/NNPESi/NNAtomnum,NNAtoms(10),NFFAtoms(10) !2024/3/7 CZZ
      DIMENSION FT(NDA3),QTP1(NDA3),QTM1(NDA3),QTP2(NDA3),QTM2(NDA3)
      DIMENSION QC(NDA3),FT2(NDA3)
      SAVE QC,QTM1,QTM2,QTP1,FT2
	  real*8:: H_bf,H_aft,OHval(4),Qtemp(3)
!########################################
! Added by Zexing Qu 2023.7.10
      COMMON/SH/SH_ENG
!########################################
!
!   VELOCITY VERLET
!
      IF (NC.EQ.1.OR.NVERLET.EQ.1) THEN
!--> CZZ 2024/3/7 判断哪个H离O最近，将它SWAP到H3

	CALL SWAP

	!--> End CZZ
	
!   IF NC=1 THEN DO NORMAL VERLET AND STORE PDOT TO FT2

!.....ADDED BY BIN FOR LANGEVIN EQUATIONS BY BBK INTEGRATOR OF VELOCITY VERLET FORMULATION
!.....ALGORITHM TAKEN FROM REFERENCE OF JENSEN AND FARAGO, MOLECULAR PHYSICS, 111, 983-991 (2013)
        IF (NGLO.NE.0) THEN
          DO K=1,NI-3
            KK=(K+2)/3
            FT2(K)=PDOT(K)
            Q(K)=Q(K)+P(K)*ATIME/W(KK)+PDOT(K)/(2.D0*W(KK))*ATIME**2
          ENDDO

          DO K=NI-2,NI
            KK=(K+2)/3
            GN(K)=GASDEV(IDUM)*GSW
            FT2(K)=PDOT(K)
            FT(K)=PDOT(K)-FCG*P(K)+GN(K)/ATIME
            Q(K)=Q(K)+COEFB*ATIME*P(K)/W(KK)+COEFB*ATIME**2*0.5D0*PDOT(K)/W(KK)+COEFB*ATIME*0.5D0*GN(K)/W(KK)
          ENDDO
          
        ELSE       ! FOR STANDARD VELOCITY VERLET
			
            Qt1 = Q  !#######ADDED BY CZZ
			
          DO K=1,NI
            KK=(K+2)/3
            FT2(K)=PDOT(K)
            Q(K)=Q(K)+P(K)*ATIME/W(KK)+PDOT(K)/(2.D0*W(KK))*ATIME**2
          ENDDO
	
		!	CALL Zmatrix(ZMat) !#### Added by CZZ, Calculate Z matrix !!!!!!Revist this!
		!	CALL Force_MFF(ZMat) !####Molecular Force Field
        ENDIF

!...UP DATE VELOCITY (MOMENTUM)
        CALL DVDQ

!.....  FOR LDFA MODEL
        IF (NFC.NE.0) THEN
           CALL FRICTION(NFC)
           CALL FRICFORCE(NFC)
        END IF
!.....  END 

!.....  FOR GLO MODEL
        IF (NGLO.NE.0) THEN
           DO I=1,3
             J=3*(NATOMA(1))+I
             K=3*(NATOMA(1)+1)+I
             I1=NATOMA(1)+1
             I2=NATOMA(1)+2
             PDOT(J)=PDOT(J)-2D0*WS2(I)*W(I1)*Q(J)+WGS2(I)*W(I1)*Q(K)
             PDOT(K)=PDOT(K)-2D0*WG2(I)*W(I2)*Q(K)+WGS2(I)*W(I2)*Q(J)
           ENDDO
!
!   P= MV SO P(T+DT)=P(T)+(F(T+DT)+F(T))*DT/2
!
          DO K=1,NI-3
            P(K)=P(K)+(FT2(K)+PDOT(K))*ATIME*0.5D0
          ENDDO

          DO K=NI-2,NI
            P(K)=COEFA*P(K)+ATIME*0.5D0*(COEFA*FT2(K)+PDOT(K))+COEFB*GN(K)
          ENDDO

			IF (NVERLET.NE.1) THEN
            DO K=NI-2,NI
              FT2(K)=FT(K)
              GN(K)=GASDEV(IDUM)*GSW/ATIME
              PDOT(K)=PDOT(K)-FCG*P(K)+GN(K)
            ENDDO
          ENDIF

        ELSE       ! FOR STANDARD VELOCITY VERLET，真
!
!   P= MV SO P(T+DT)=P(T)+(F(T+DT)+F(T))*DT/2
!
          DO K=1,NI
            P(K)=P(K)+(FT2(K)+PDOT(K))*ATIME*0.5D0
          ENDDO
		  
        
        DO K=1,NI     !计算Qt3
            KK=(K+2)/3
            Qt3(K)=Q(K)+P(K)*ATIME/W(KK)+PDOT(K)/(2.D0*W(KK))*ATIME**2
        ENDDO

	!#### Added by CZZ,for surface hopping	 
        !write(*,*)'!!!!!!!!N1.5',nPES 
		  Call TSH_Algorithm 
	!####

        ENDIF
!
!*********************************************************************
!   BEEMAN'S THIRD ORDER PROCEDUCE, BEEMAN, J. COMPUT. PHYS. 20, 130
!   (1976), MODIFIED BY TULLY ET AL., J. CHEM. PHYS. 71, 1630 (1979)
!********************************************************************
!
      ELSE 

!   IF NC>1 THEN DO BEEMAN'S THIRD ORDER PROPAGATION, FT2 AND QTP1 ARE
!   SAVED AS THE FORCES AT THE (N-1)TH STEP AND THE COORDINATES AT THE
!   NTH STEP
        IF (NGLO.NE.0) THEN

          DO K=NI-2,NI
            KK=(K+2)/3
            QC(K)=P(K)+0.5D0*ATIME*(3*PDOT(K)-FT2(K))
          ENDDO

        ENDIF

          DO K=1,NI
            KK=(K+2)/3
            QTP1(K)=Q(K)
            Q(K)=Q(K)+P(K)*ATIME/W(KK)+ATIME**2*(4*PDOT(K)-FT2(K))/W(KK)/6
            FT2(K)=PDOT(K)
          ENDDO

!...UPDATE VELOCITY (MOMENTUM)
          CALL DVDQ

!.....FOR LDFA MODEL
          IF (NFC.NE.0) THEN
             CALL FRICTION(NFC)
             CALL FRICFORCE(NFC)
          END IF
!.....END 

!.....FOR GLO MODEL
        IF (NGLO.NE.0) THEN
          DO I=1,3
            J=3*(NATOMA(1))+I
            K=3*(NATOMA(1)+1)+I
            I1=NATOMA(1)+1
            I2=NATOMA(1)+2
            PDOT(J)=PDOT(J)-2D0*WS2(I)*W(I1)*Q(J)+WGS2(I)*W(I1)*Q(K)
            PDOT(K)=PDOT(K)-2D0*WG2(I)*W(I2)*Q(K)+WGS2(I)*W(I2)*Q(J)
            GN(K)=GASDEV(IDUM)*GSW/ATIME
            PDOT(K)=PDOT(K)-FCG*QC(K)+GN(K)
          ENDDO
        ENDIF
!.....END 

           DO K=1,NI
            KK=(K+2)/3
            P(K)=W(KK)*(Q(K)-QTP1(K))/ATIME+(2*PDOT(K)+FT2(K))*ATIME/6
          ENDDO

      ENDIF     ! BY BIN, 2017/2/16
!############################
! ADDED by Zexing Qu 2023/7/10
      IF(DABS(SH_ENG)>0.0D0) THEN
!       P(:)=P(:)*scalingfactor(SH_ENG)
      ENDIF
!#############################

      RETURN
      END
	  
	  
