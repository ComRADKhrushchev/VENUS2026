
      function grandom() result(zz) 
      ! generate a random normal deviate using the Box-Muller method 
      real*8              :: zz 
      real*8              :: dist,xx(2) 
      real*8, parameter   :: two = 2d0, minus_two = -2d0 
      do 
         call random_number(xx) 
         xx   = two*xx - 1d0 
         dist = xx(1)**2 + xx(2)**2 
         if (dist < 1d0) exit 
      end do 
      zz = xx(1)*sqrt(minus_two*log(dist)/dist) 
      end function grandom 
c
c	normal distribution random number generator taken
c	from 'numerical recipe in fortran 77'
c	Copyright (C) 1986-1992 Numerical Recipes Software
c
      FUNCTION ran1(idum)  
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV  
      REAL ran1,AM,EPS,RNMX  
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,  
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)  
      INTEGER j,k,iv(NTAB),iy  
      SAVE iv,iy  
      DATA iv /NTAB*0/, iy /0/  
      if (idum.le.0.or.iy.eq.0) then  
        idum=max(-idum,1)  
        do 11 j=NTAB+8,1,-1  
          k=idum/IQ  
          idum=IA*(idum-k*IQ)-IR*k  
          if (idum.lt.0) idum=idum+IM  
          if (j.le.NTAB) iv(j)=idum  
11      continue  
        iy=iv(1)  
      endif  
      k=idum/IQ  
      idum=IA*(idum-k*IQ)-IR*k  
      if (idum.lt.0) idum=idum+IM  
      j=1+iy/NDIV  
      iy=iv(j)  
      iv(j)=idum  
      ran1=min(AM*iy,RNMX)  
      return  
      END  
 
      FUNCTION gasdev_0(idum)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER idum  
      REAL*8 gasdev  
CU    USES ran1  
      INTEGER iset  
      REAL fac,gset,rsq,v1,v2,ran1  
      SAVE iset,gset  
      DATA iset/0/  
      if (iset.eq.0) then  
1       v1=2.*ran1(idum)-1.  
        v2=2.*ran1(idum)-1.  
        rsq=v1**2+v2**2  
        if(rsq.ge.1..or.rsq.eq.0.)goto 1  
        fac=sqrt(-2.*log(rsq)/rsq)  
        gset=v1*fac  
        gasdev=v2*fac  
        iset=1  
      else  
        gasdev=gset  
        iset=0  
      endif  
      return  
      END  
