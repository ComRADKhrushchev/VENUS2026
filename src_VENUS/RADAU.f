      SUBROUTINE RADAU(TF,XL,LL,NIP,NITER)
c------------------------ki 
c  NITER is included for the time being
c******************************************************************
c
c RADAU (RA15)
c A 15th order integrator that uses Gauss-Radau spacings.
c   This means that it does 7 force calculations and on the
c   first time through it does 6 stes of 7.
c
c Reference : Edgar Everhart, 'an efficient integrator that uses
c             Gauss-Radau spacings'. In 'Dynamics of comets:
c             their origin and evolution', 185-202,
c             A. Carusi and G.B. Valsecci (eds.),
c             D. Reidel publishing company 1985.
c
c                                                N. Markovic, 1992
c
c******************************************************************
      use venus_params
      use venus_data, only: D_H => D, V_H => V, RA_H => RA, W_H => W,
     *                     VW => W,
     *                     LL_H => LL, LA_H => LA, LB_H => LB,
     *                     TT_H => TT, Q_H => Q, TEMP_H => TEMP
      use venus_data
      implicit double precision(a-h,o-z)
C******************************************************************
C         A 15TH ORDER INTEGRATOR THAT USES GAUSS-RADAU SPACINGS.
C         THIS MEANS THAT EACH SEQUENCE IS SPACED INTO 7 SUBSTEPS,
C         THE POSITION AND VELOCITY ARE ACCURATE TO THE 15TH ORDER.
C         THE FIRST SEQUENCE IS ITERNATED 6 TIMES AND LATER SEQUENCES
C         ARE ITERNATED 2 TIMES.
C         IN OTHER WORDS, IT DOES 7 FORCE CALCULATIONS FOR EACH
C         ITERATION OF SEQUENCE AND ON THE FIRST SEQUENCE (6
C         ITERATIONS) THROUGH IT DOES 1+6*7, WHILE ON THE LATER
C         SEQUENCES (2 ITERATIONS) THROUGH IT DOES 1+2*7.
C         NEW SEQUENCE WILL BE CONTINOUSLY GENERATED AS LONG AS
C         TOTAL INTEGRATION LENGTH IS LESS THAN THE GIVEN TIME
C         INTERVAL (TF).
C
C         ADAPTED FROM SUBROUTINE RA15 IN REF:
C             EDGAR EVERHART, 'AN EFFICIENT INTEGRATOR THAT USES
C             GAUSS-RADAU SPACINGS'. IN 'DYNAMICS OF COMETS:
C             THEIR ORIGIN AND EVOLUTION', 185-202,
C             A. CARUSI AND G.B. VALSECCI (EDS.),
C             D. REIDEL PUBLISHING COMPANY 1985.
C         WITH SLIGHT MODIFICATION SEQUENTLY BY N. MARKOVIC(1992),
C         K. BOLTON(1997) AND G. LI(1999)
C
C         NOTE THAT NAMES OF SOME VARIABLES IN RA15 HAVE BEEN CHANGED
C         BECAUSE THEY ARE ALREADY DEFINED AS OTHER VARIABLES IN VENUS.
C         THESE VARIABLES AND THEIR CHANGES ARE: NI->NITER, NV->NI,
C         Q->QT, X->Q, W->WC, WW->WTMP, H->HH AND T->TT.
C
C         TF IS THE LENGTH OF INTEGRATION (NEGATIVE WHEN INTEGRATING
C         BACKWARD). LL CONTROLS ACCURACY. THUS IF SS=10.**(-LL)
C         CONTROLS THE SIZE OF THE LAST TERM IN A SERIES. TRY LL=8 AND
C         WORK UP OR DOWN FROM THERE. HOWEVER, IF LL.LT.0, THEN XL IS
C         THE CONSTANT SEQUENCE SIZE USED. A NON-ZERO XL SETS THE SIZE
C         OF THE FIRST SEQUENCE REGARDLESS OF LL'S SIGN. NIP IS THE
C         NUMBER OF SEQUENCES BETWEEN INTERMEDIATE PRINTOUT.
C         TYPE OF DIFFERENTIAL EQUATION:
C                 NCLASS=1  :  Y' = F(Y,T)
C                 NCLASS=-2 :  Y''= F(Y,T)
C                 NCLASS=2  :  Y''= F(Y',Y,T)
C************************************************************************
c
c     nda = the no of atoms
c     These are commented out for the inclusion of sizes file-------ki
c     parameter (nda=1000)
c     parameter (nda3=nda3)
c
      real*4 tval,pw
      parameter (pw=1./9.)
c
C      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
C     *NPATH,NAST
c
c
      dimension cc(21),d(21),ra(21),hh(8),w(7),u(7),nw(8)
      dimension f1(nda3),fj(nda3),yy(nda3),zz(nda3),
     2          bb(7,nda3),gg(7,nda3),ee(7,nda3),bd(7,nda3),
     3          v(nda3)
      dimension xxsav(nda3),xx(nda3)
c
      logical npq,nsf,nper,ncl,nes
c
      save nw,zero,half,one,sr,hh
c
      data nw/0,0,1,3,6,10,15,21/
      data zero,half,one,sr/0.0d0,0.5d0,1.0d0,1.4d0/
      data hh/       0.0d0, .05626256053692215d0, .18024069173689236d0,
     *.35262471711316964d0, .54715362633055538d0, .73421017721541053d0,
     *.88532094683909577d0, .97752061356128750d0/
c
c*****************  initialization  *************************
      nclass=-2
c     printing every nip steps
      prstep=atime*dble(nip)
c      prstep=xl*dble(nip)/20.0d0
      prtime=prstep
c      write(6,*)'AS we enter integ'
      do 223 i=1,nv
         kn=1+(i-1)/3
         v(i)=p(i)/vw(kn)
c        write(6,*)i,kn,xx(i),v(i),vw(kn),p(i)
 223  continue
 666  format(i4,4f10.4)
c     call the RADOUT at the beginning of the trajectory
      
      if(ll.ge.0) then
         write(*,*)'CALLING RADOUT 0'
         tmstart=0.d0
         call RADOUT(tmstart,tf)
      endif

c      write(6,*)'IN RADAU 1',ll,tf,xl
c************************************************************
c
c
c     tcorr=10.d0**(-112.d0/9)
      tcorr=1.0d0
c
      nper=.false.
      nsf=.false.
      ncl=nclass.eq.1
      npq=nclass.lt.2
c y'=f(y,t) ncl=.true.  y''=f(y,t) ncl=.false.  y''=f(y',y,t) ncl=.false.
c nclass=1  npq=.true.  nclass=-2  npq=.true.   nclass=2      npq=.false.
c nsf is .false. on starting sequence, otherwise .true.
c nper is .true. only on last sequence of the integration
c nes is .true. only if ll is negative. then the sequence size is xl.
      dir=one
      if(tf.lt.zero) dir=-one
      nes=ll.lt.0
      xl=dabs(xl)*dir
c     pw=1./9.
      do 14 n=2,8
      ww=n+n*n
      if(ncl) ww=n
      w(n-1)=one/ww
      ww=n
   14 u(n-1)=one/ww
      do 22 k=1,nv
      if(ncl) v(k)=zero
      do 22 l=1,7
C     E NOT ZEROED IN EVERHART
      ee(l,k)=zero
      bd(l,k)=zero
   22 bb(l,k)=zero
      w1=half
      if(ncl) w1=one
      cc(1)=-hh(2)
      d(1)=hh(2)
      ra(1)=one/(hh(3)-hh(2))
      la=1
      lc=1
      do 73 k=3,7
      lb=la
      la=lc+1
      lc=nw(k+1)
      cc(la)=-hh(k)*cc(lb)
      cc(lc)=cc(la-1)-hh(k)
      d(la)=hh(2)*d(lb)
      d(lc)=-cc(lc)
      ra(la)=one/(hh(k+1)-hh(2))
      ra(lc)=one/(hh(k+1)-hh(k))
      if(k.eq.3) goto 73
      do 72 l=4,k
      ld=la+l-3
      le=lb+l-4
      cc(ld)=cc(le)-hh(k)*cc(le+1)
      d(ld)=d(le)+hh(l-1)*d(le+1)
   72 ra(ld)=one/(hh(k+1)-hh(l-1))
   73 continue
      ss=10.**(-ll)
      tp=0.1d0*dir
      if(xl.ne.zero) tp=xl
C     THIS LINE NOT INCLUDED IN E (below)
      if(nes) tp=xl
c*******  start of first sequence  **************************
c
 4000 ns=0
      nf=0
      ni=6
c
c***  initialization  tm=collision time
c***************************************************************
      tm=zero
c      write(6,*)'CALLING FROM RADAU 1'
      CALL DVDQ_1
      do  i=1,nv
         kn=1+(i-1)/3
         f1(i)=pdot(i)/vw(kn)
      enddo
c*** start of all sequences but the first  *********************
  722 nf=nf+1
  591 do 58 k=1,nv
      gg(1,k)=bb(1,k)+d(1)*bb(2,k)+
     *  d(2)*bb(3,k)+d(4)*bb(4,k)+d(7)*bb(5,k)+d(11)*bb(6,k)+d(16)*
     *  bb(7,k)
      gg(2,k)=            bb(2,k)+
     *  d(3)*bb(3,k)+d(5)*bb(4,k)+d(8)*bb(5,k)+d(12)*bb(6,k)+d(17)*
     *  bb(7,k)
      gg(3,k)=bb(3,k)+d(6)*bb(4,k)+d(9)*bb(5,k)+d(13)*bb(6,k)+d(18)*
     *  bb(7,k)
      gg(4,k)=            bb(4,k)+d(10)*bb(5,k)+d(14)*bb(6,k)+d(19)*
     *  bb(7,k)
      gg(5,k)=                         bb(5,k)+d(15)*bb(6,k)+d(20)*
     *  bb(7,k)
      gg(6,k)=                                      bb(6,k)+d(21)*
     *  bb(7,k)
   58 gg(7,k)=                                              bb(7,k)
      tt=tp
      t2=tt*tt
      if(ncl) t2=tt
c     dtval=dabs(tt*1.d+16)
      dtval=dabs(tt*1.d0)
      tval=real(dtval)
      do 175 m=1,ni
      do 174 j=2,8
      jd=j-1
c     jdm=j-2
      s=hh(j)
      q=s
      if(ncl) q=one
C     W() AND U() ARE SWOPPED AROUND IN E
      do 130 k=1,nv
      a=w(3)*bb(3,k)+s*(w(4)*bb(4,k)+s*(w(5)*bb(5,k)+s*(w(6)*bb(6,k)+
     *   s*w(7)*bb(7,k))))
      yy(k)=xx(k)+q*(tt*v(k)+t2*s*(f1(k)*w1+s*(w(1)*bb(1,k)+
     *  s*(w(2)*bb(2,k)+s*a))))
      if(npq) goto 130
      a=u(3)*bb(3,k)+s*(u(4)*bb(4,k)+s*(u(5)*bb(5,k)+s*(u(6)*bb(6,k)+
     *    s*u(7)*bb(7,k))))
      zz(k)=v(k)+s*tt*(f1(k)+s*(u(1)*bb(1,k)+s*(u(2)*bb(2,k)+s*a)))
Ckim      a=u(3)*bb(3,k)+s*(u(4)*bb(4,k)+s*(u(5)*bb(5,k)+s*(u(6)*bb(6,k)+
c     *   s*u(7)*bb(7,k))))
c      yy(k)=xx(k)+q*(tt*v(k)+t2*s*(f1(k)*w1+s*(u(1)*bb(1,k)+
c     *  s*(u(2)*bb(2,k)+s*a))))
c      if(npq) goto 130
c      a=w(3)*bb(3,k)+s*(w(4)*bb(4,k)+s*(w(5)*bb(5,k)+s*(w(6)*bb(6,k)+
c     *    s*w(7)*bb(7,k))))
c      zz(k)=v(k)+s*tt*(f1(k)+s*(w(1)*bb(1,k)+s*(w(2)*bb(2,k)+s*a)))
  130 continue

c      write(6,*)'IN RADAU 2',tf,tp,tm

c      write(6,*)'CALLING FROM RADAU 2'
      do i=1,nv
         xxsav(i)=xx(i)
         xx(i)=yy(i)
      enddo

      CALL DVDQ_1
      do  i=1,nv
         xx(i)=xxsav(i)
         kn=1+(i-1)/3
         fj(i)=pdot(i)/vw(kn)
      enddo
      nf=nf+1
      do 171 k=1,nv
      temp=gg(jd,k)
      gk=(fj(k)-f1(k))/s
      goto (102,102,103,104,105,106,107,108),j
  102 gg(1,k)=gk
      goto 160
  103 gg(2,k)=(gk-gg(1,k))*ra(1)
      goto 160
  104 gg(3,k)=((gk-gg(1,k))*ra(2)-gg(2,k))*ra(3)
      goto 160
  105 gg(4,k)=(((gk-gg(1,k))*ra(4)-gg(2,k))*ra(5)-gg(3,k))*ra(6)
      goto 160
  106 gg(5,k)=((((gk-gg(1,k))*ra(7)-gg(2,k))*ra(8)-gg(3,k))*ra(9)-
     *     gg(4,k))*ra(10)
      goto 160
  107 gg(6,k)=(((((gk-gg(1,k))*ra(11)-gg(2,k))*ra(12)-gg(3,k))*ra(13)-
     *     gg(4,k))*ra(14)-gg(5,k))*ra(15)
      goto 160
  108 gg(7,k)=((((((gk-gg(1,k))*ra(16)-gg(2,k))*ra(17)-gg(3,k))*ra(18)-
     *     gg(4,k))*ra(19)-gg(5,k))*ra(20)-gg(6,k))*ra(21)
  160 temp=gg(jd,k)-temp
      bb(jd,k)=bb(jd,k)+temp
      goto (171,171,203,204,205,206,207,208),j
  203 bb(1,k)=bb(1,k)+cc(1)*temp
      goto 171
  204 bb(1,k)=bb(1,k)+cc(2)*temp
      bb(2,k)=bb(2,k)+cc(3)*temp
      goto 171
  205 bb(1,k)=bb(1,k)+cc(4)*temp
      bb(2,k)=bb(2,k)+cc(5)*temp
      bb(3,k)=bb(3,k)+cc(6)*temp
      goto 171
  206 bb(1,k)=bb(1,k)+cc(7)*temp
      bb(2,k)=bb(2,k)+cc(8)*temp
      bb(3,k)=bb(3,k)+cc(9)*temp
      bb(4,k)=bb(4,k)+cc(10)*temp
      goto 171
  207 bb(1,k)=bb(1,k)+cc(11)*temp
      bb(2,k)=bb(2,k)+cc(12)*temp
      bb(3,k)=bb(3,k)+cc(13)*temp
      bb(4,k)=bb(4,k)+cc(14)*temp
      bb(5,k)=bb(5,k)+cc(15)*temp
      goto 171
  208 bb(1,k)=bb(1,k)+cc(16)*temp
      bb(2,k)=bb(2,k)+cc(17)*temp
      bb(3,k)=bb(3,k)+cc(18)*temp
      bb(4,k)=bb(4,k)+cc(19)*temp
      bb(5,k)=bb(5,k)+cc(20)*temp
      bb(6,k)=bb(6,k)+cc(21)*temp
  171 continue
  174 continue
      if(nes.or.m.lt.ni) goto 175
      hv=zero
      do 635 k=1,nv
  635 hv=dmax1(hv,dabs(bb(7,k)))
      hv=hv*w(7)/tval**7
  175 continue
      if(nsf) goto 180
      if(.not.nes) tp=tcorr*(ss/hv)**pw*dir
      if(nes) tp=xl
      if(nes) goto 170
      if(tp/tt.gt.one) goto 170
      tp=.8d0*tp
c      write(6,8)tm,tt,tp
c    8 format(' restart. tm,tt,tp',3f15.5)
      goto 4000
  170 nsf=.true.
  180 do 35 k=1,nv
      xx(k)=xx(k)+v(k)*tt+t2*(f1(k)*w1+bb(1,k)*w(1)+bb(2,k)*
     *w(2)+bb(3,k)*w(3)+bb(4,k)*w(4)+bb(5,k)*w(5)+bb(6,k)*w(6)+
     * bb(7,k)*w(7))
      if(ncl) goto 35
      v(k)=v(k)+tt*(f1(k)+bb(1,k)*u(1)+bb(2,k)*u(2)+bb(3,k)*u(3)
     *    +bb(4,k)*u(4)+bb(5,k)*u(5)+bb(6,k)*u(6)+bb(7,k)*u(7))
   35 continue
      tm=tm+tt
c      write(6,*)'CALLING FROM RADAU 3'
      CALL DVDQ_1
      do 555 i=1,nv
         kn=1+(i-1)/3
         f1(i)=pdot(i)/vw(kn)
 555  continue
 637  continue
c
c***************  test section   ******************************
c
c
c
c
c
c=================================================================
c========  place new tests, calculation of energies etc.  ========
c          in this region.
c                          coordinates:  x
c                          velocities:   v
c                          distances:    r
c                          potential e:  epot
c
c
      do  i=1,nv
         kn=1+(i-1)/3
         p(i)=v(i)*vw(kn)
      enddo
      if((ll.ge.0).and.(dabs(tm).ge.prtime-1.d-8)) then
         prtime=prtime+prstep
         call RADOUT(tm,tf)
         if(tm.eq.tf.or.ntest.eq.2) return
      endif
c=================================================================
c=================================================================
c
  636 ns=ns+1
      if(.not.nper) goto 78
      do i=1,nv
         kn=1+(i-1)/3
         p(i)=v(i)*vw(kn)
      enddo
      if(ll.ge.0) then
         write(*,*)'CALLING RADOUT 2'
         call RADOUT(tm,tf)
         if(tm.eq.tf.or.ntest.eq.2.or.ntest.eq.2) return
      endif

      return
   78 continue
c****************  end of test section  ******************************
      nf=nf+1
      if(nes) goto 341
      Tp=dir*tcorr*(ss/hv)**pw
      if(tp/tt.gt.sr) tp=tt*sr
  341 if(nes) tp=xl
C     THE LINE BELOW IS  A BIT DIFFERENT TO E
      if(dir*(tf-tm-tp).gt.xl) goto 77
      tp=tf-tm
      nper=.true.
   77 q=tp/tt
      do 39 k=1,nv
      if(ns.eq.1) goto 31
      do 20 j=1,7
   20 bd(j,k)=bb(j,k)-ee(j,k)
   31 ee(1,k)=      q*(bb(1,k)+2.d0*bb(2,k)+3.d0*bb(3,k)+
     *           4.d0*bb(4,k)+5.d0*bb(5,k)+6.d0*bb(6,k)+7.d0*
     *   bb(7,k))
      ee(2,k)=               q**2*(bb(2,k)+3.d0*bb(3,k)+
     *           6.d0*bb(4,k)+10.d0*bb(5,k)+15.d0*bb(6,k)+21.d0*
     *   bb(7,k))
      ee(3,k)=                             q**3*(bb(3,k)+
     *           4.d0*bb(4,k)+10.d0*bb(5,k)+20.d0*bb(6,k)+35.d0*
     *   bb(7,k))
      ee(4,k)=  q**4*(bb(4,k)+5.d0*bb(5,k)+15.d0*bb(6,k)+35.d0*
     *  bb(7,k))
      ee(5,k)=              q**5*(bb(5,k)+6.d0*bb(6,k)+21.d0*
     *   bb(7,k))
      ee(6,k)=                       q**6*(bb(6,k)+7.d0*bb(7,k))
      ee(7,k)=                                      q**7*bb(7,k)
      do 39 l=1,7
   39 bb(l,k)=ee(l,k)+bd(l,k)
      ni=2
      goto 722
  638 tm=tm+tt
c
      do i=1,nv
         kn=1+(i-1)/3
         p(i)=v(i)*vw(kn)
      enddo
      if(ll.ge.0) then
         write(*,*)'CALLING RADOUT 3'
         call RADOUT(tm,tf)
         if(tm.eq.tf.or.ntest.eq.2) return
      endif
c
      return
      end
c
c  radout.f
c
      SUBROUTINE RADOUT(tme,tf)
      implicit double precision (a-h,o-z)
c     PARAMETER(ND1=1000,NDP=10)
C      COMMON/TESTB/RMAX(NDP),RBAR(NDP),NTEST,NPATHS,NABJ(NDP),NABK(NDP),
C     *NPATH,NAST
c
c     this is where output will be written when using the variable
c     time step integrator
c
c     modified to use venus GWRITE and FINAL ------------------ki
c
      write(6,*)'TRAJECTORY NUMBER ',ntz
      time=tme
      write(6,*)'The time is ',time*10.d0,' fs'
c     write(6,*)'Total energy ',h
c     write(6,*)'Potential energy ',v
c     write(6,*)'Kinetic energy ',t
      call GWRITE
      call test
c     write(6,*)' in radout,ntest,nast,npath=',ntest,nast,npath
      IF (NTEST.EQ.1.AND.NAST.EQ.0) GOTO 406
      IF (NTEST.EQ.2.AND.NAST.EQ.1) GOTO 410
      IF (NTEST.EQ.2.AND.NAST.EQ.0) GOTO 410
      IF (NTEST.EQ.0.AND.NAST.NE.0) NAST=0
      IF (NTEST.EQ.1.AND.NAST.EQ.2) NAST=1
      nc=nc+1
      RETURN
  406 NAST=1
      WRITE(6,*)
  910 FORMAT(10X,'REACTION OCCURRED FOR PATH',I3)
      WRITE(6,910)NPATH
      CALL ENERGY_1
      CALL GWRITE
      RETURN
  410 CONTINUE
      CALL DVDQ_1
      CALL ENERGY_1
      CALL GWRITE
  414 CONTINUE
      CALL FINAL
      CALL GFINAL
c     if(ntest.eq.2) then
c        write(6,*)'THE CRITERIA TO STOP THE TRAJ ARE FULFILLED'
c        tme=tf
c        do i=1,nv,3
c           write(6,*)q(i),q(i+1),q(i+2)
c        enddo
c        do i=1,nv,3
c           write(6,*)p(i),p(i+1),p(i+2)
c        enddo
c     endif
c     if(tme.eq.0.0d0) then
c        do i=1,nv,3
c           write(6,*)q(i),q(i+1),q(i+2)
c        enddo
c        do i=1,nv,3
c           write(6,*)p(i),p(i+1),p(i+2)
c        enddo
c     endif
      return
      END
