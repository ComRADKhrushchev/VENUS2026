c        implicit real*8 (a-h,o-z)
c        real*8 p1(3),p2(3),x1(3),x2(3)
c        p1(1)=6d0
c        p1(2)=-2d0
c        p1(3)=0d0
c        p2(1)=12d0
c        p2(2)=8d0
c        p2(3)=0d0
c        do the=0d0,360d0,10d0
c        x1(1)=10d0
c        x1(2)=6d0
c        x1(3)=0d0
c        call rotaline(p1,p2,x1,x2,the)
c        write(*,*)x2(1),x2(2),x2(3)
c        enddo
c        end

       subroutine rotaline(p1,p2,x1,x2,th)
!-->...rotate about an arbitrary line
!-->...p1(x,y,z) and p2(x,y,z) are two points define this line (vector:p1->p2)
!-->...x1(x,y,z) is a point before rotation, x2 (x,y,z) is a point after rotation about p1->p2 line
!-->...th is the rotation angle (anti-clockwise )
!-->...reference: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
        implicit real*8 (a-h,o-z)
        real*8 p1(3),p2(3),x1(3),x2(3)
        real*8 trmt(4,4)
        pi=dacos(-1d0)
        the=th*pi/180d0
        a=p1(1)
        b=p1(2)
        c=p1(3)
        dist=dsqrt((p2(1)-p1(1))**2+(p2(2)-p1(2))**2+(p2(3)-p1(3))**2)
        u=(p2(1)-p1(1))/dist
        v=(p2(2)-p1(2))/dist
        w=(p2(3)-p1(3))/dist

c!-->    obtain the transformation matrix first
c        trmt(1,1)=u*u+(v*v+w*w)*dcos(the)
c        trmt(2,1)=u*v*(1-dcos(the))-w*dsin(the)
c        trmt(3,1)=u*w*(1-dcos(the))+v*dsin(the)
c        trmt(4,1)=(a*(v*v+w*w)-u*(b*v+c*w))*(1-dcos(the))+(b*w-c*v)
c     &*dsin(the)
c        trmt(1,2)=u*v*(1-dcos(the))+w*dsin(the)
c        trmt(2,2)=v*v+(u*u+w*w)*dcos(the)
c        trmt(3,2)=v*w*(1-dcos(the))-u*dsin(the)
c        trmt(4,2)=(b*(u*u+w*w)-v*(a*u+c*w))*(1-dcos(the))+(c*u-a*w)
c     &*dsin(the)
c        trmt(1,3)=u*w*(1-dcos(the))-v*dsin(the)
c        trmt(2,3)=v*w*(1-dcos(the))+u*dsin(the)
c        trmt(3,3)=w*w+(u*u+v*v)*dcos(the)
c        trmt(4,3)=(c*(u*u+v*v)-w*(a*u+b*v))*(1-dcos(the))+(a*v-b*u)
c     &*dsin(the)
c        trmt(1,4)=0 
c        trmt(2,4)=0 
c        trmt(3,4)=0 
c        trmt(4,4)=1 
c
c!-->    obtain the coordinates after rotation then
c        x2(1)=x1(1)*trmt(1,1)+x1(2)*trmt(2,1)+x1(3)*trmt(3,1)+trmt(4,1)
c        x2(2)=x1(1)*trmt(1,2)+x1(2)*trmt(2,2)+x1(3)*trmt(3,2)+trmt(4,2)
c        x2(3)=x1(1)*trmt(1,3)+x1(2)*trmt(2,3)+x1(3)*trmt(3,3)+trmt(4,3)

!-->    obtain the coordinates after rotation directly
        x=x1(1)
        y=x1(2)
        z=x1(3)
        x2(1)=(a*(v*v+w*w)-u*(b*v+c*w-u*x-v*y-w*z))*(1d0-dcos(the))
     &+x*dcos(the)+(-c*v+b*w-w*y+v*z)*dsin(the)
        x2(2)=(b*(u*u+w*w)-v*(a*u+c*w-u*x-v*y-w*z))*(1d0-dcos(the))
     &+y*dcos(the)+(c*u-a*w+w*x-u*z)*dsin(the)
        x2(3)=(c*(u*u+v*v)-w*(a*u+b*v-u*x-v*y-w*z))*(1d0-dcos(the))
     &+z*dcos(the)+(-b*u+a*v-v*x+u*y)*dsin(the)

        return
        end subroutine
