!!!     set natoms=2 in your main program
!!!     Q(1:3) for x, y, and z coordiantes of H1 
!!!     Q(4:6) for x, y, and z coordiantes of H1 
!!!     V is the potential energy
        subroutine POT0(NATOMS,Q,V)
        implicit real*8 (a-h,o-z)
        parameter(npoly=13)
        parameter(rbohr=0.5291772d0)
        dimension Q(natoms*3)
        dimension cart(3,natoms),scoor(npoly)
        do i=1,natoms
        cart(1,i)=Q((i-1)*3+1)
        cart(2,i)=Q((i-1)*3+2)
        cart(3,i)=Q((i-1)*3+3)
        enddo
        call symmetrize(natoms,npoly,cart,scoor)
        call getpot(scoor,vpot)
        V=vpot

        return
        end subroutine POT0

        subroutine DPESHON(NATOMS,Q,PDOT)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        dimension DV(NATOMS*3,2),PDOT(NATOMS*3),Q(NATOMS*3)
        HINC=1.0d-4

        DO I=1,NATOMS*3
         Q(I)=Q(I)+HINC
         CALL pot0(NATOMS,Q,V)
         Q(I)=Q(I)-HINC
         DV(I,1)=V
         Q(I)=Q(I)-HINC
         CALL pot0(NATOMS,Q,V)
         Q(I)=Q(I)+HINC
         DV(I,2)=V
         PDOT(I)=(DV(I,1)-DV(I,2))/2.d0/HINC
        ENDDO
        RETURN
        end


!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
        module nnparam
        implicit none
        integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
        integer nscale
        integer, allocatable::nodes(:)
        real*8, allocatable::weight(:,:,:),bias(:,:)
        real*8, allocatable::pdel(:),pavg(:)
        end module nnparam

!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weights.txt'
!-->  biases saved in 'biases.txt'
!-->  one has to call this subroutine once and only once before calling the getpot() subroutine
        subroutine pes_init
        use nnparam
        implicit none
        integer i,ihid,iwe,inode1,inode2,ilay1,ilay2
        character f1*80
        open(111,file='weights-h2ag.txt')
        open(222,file='biases-h2ag.txt')
        read(111,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2
        allocate(nodes(nlayer),pdel(nscale),pavg(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(111,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
        nodemax=max(nodemax,nodes(i))
        enddo
       allocate(weight(nodemax,nodemax,2:nlayer),bias(nodemax,2:nlayer))
        read(111,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
!-->....see function tranfun() for details
        read(111,*)(pdel(i),i=1,nscale)
        read(111,*)(pavg(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2)
        iwe=iwe+1
        enddo
        iwe=iwe+1
        enddo
        enddo
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2)
        read(111,*)weight(inode2,inode1,ilay1)
        enddo
        read(222,*)bias(inode1,ilay1)
        enddo
        enddo
        write(*,*)'read all parameters done'
        close(111)
        close(222)
        end subroutine pes_init

        subroutine getpot(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
        do i=1,ninput
        y(i,1)=(x(i)-pavg(i))/pdel(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=bias(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weight(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=bias(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weight(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdel(nscale)+pavg(nscale)
        return
        end

        function tranfun(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun,x
c    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
c    ifunc=2, transfer function is log sigmoid function, 'logsig'
c    ifunc=3, transfer function is pure linear function, 'purelin'. It is imposed to the output layer by default
        if (ifunc.eq.1) then
        tranfun=tanh(x)
        else if (ifunc.eq.2) then
        tranfun=1d0/(1d0+exp(-x))
        else if (ifunc.eq.3) then
        tranfun=x
        endif
        return
        end

        subroutine symmetrize(natom,npoly,cart,scoor)
        implicit real*8(a-h,o-z)
        dimension cart(3,natom),scoor(npoly),coor(100)
        pi=dacos(-1d0)
        dmm=2.9437d0
        dsq3=dsqrt(3d0)
        alpha=2*pi/dmm
        gama=1d0
        xh1=cart(1,1)
        yh1=cart(2,1)
        zh1=cart(3,1)
        xh2=cart(1,2)
        yh2=cart(2,2)
        zh2=cart(3,2)
        xcom=(xh1+xh2)*0.5d0
        ycom=(yh1+yh2)*0.5d0
        zcom=(zh1+zh2)*0.5d0
        r=dsqrt((xh1-xh2)**2+(yh1-yh2)**2+(zh1-zh2)**2)
        coor(1)=dcos(2*alpha*yh1/dsq3)+2d0*dcos(alpha*xh1)
     &*dcos(alpha*yh1/dsq3)
        coor(2)=dsin(2*alpha*yh1/dsq3)-2d0*dcos(alpha*xh1)
     &*dsin(alpha*yh1/dsq3)
        coor(3)=dexp(-gama*zh1)
        coor(4)=dcos(2*alpha*yh2/dsq3)+2d0*dcos(alpha*xh2)
     &*dcos(alpha*yh2/dsq3)
        coor(5)=dsin(2*alpha*yh2/dsq3)-2d0*dcos(alpha*xh2)
     &*dsin(alpha*yh2/dsq3)
        coor(6)=dexp(-gama*zh2)
        coor(7)=dexp(-gama*r)
        scoor(1)=coor(1)+coor(4)
        scoor(2)=coor(2)+coor(5)
        scoor(3)=coor(3)+coor(6)
        scoor(4)=coor(1)*coor(4)
        scoor(5)=coor(2)*coor(5)
        scoor(6)=coor(3)*coor(6)
        scoor(7)=coor(1)*coor(2)+coor(4)*coor(5)
        scoor(8)=coor(1)*coor(3)+coor(4)*coor(6)
        scoor(9)=coor(2)*coor(3)+coor(5)*coor(6)
        scoor(10)=coor(7)
        coor(8)=dcos(2*alpha*ycom/dsq3)+2d0*dcos(alpha*xcom)
     &*dcos(alpha*ycom/dsq3)
        coor(9)=dsin(2*alpha*ycom/dsq3)-2d0*dcos(alpha*xcom)
     &*dsin(alpha*ycom/dsq3)
        coor(10)=dexp(-gama*zcom)
        scoor(11)=coor(8)
        scoor(12)=coor(9)
        scoor(13)=coor(10)
 
        return
        end subroutine symmetrize


