        implicit real*8(a-h,o-z)
        parameter (natoms=2)
        dimension q(natoms*3),dv(natoms*3)

!.......call this subroutine once and just once for initializing the PES (reading parameters)
        call pes_init

        q(1:3)=(/    1.277280E-01,    1.625873E+00,    1.064094E+00/) !the coordinate of H1
        q(4:6)=(/    1.223640E+00,    9.931286E-01,    1.125637E+00/) !the coordinate of H2
!.......call the subroutine for a given geometry q (in cartesian coordiantes, Angstrom), which yields the energy in eV
              ! q(1:2) and q(4:5) are the x,y coordinate of H1 and H2. q(3) and q(6) are the height of the H1 and H2 relative to the Ag(111) surface.
        call pot0(natoms,q,v)  !call the subroutine to obtain potential energy v in eV.
        write(*,*)v
!   1.15467106108370    (your output should be very close to this one) 

!.......for numerical gradients with two center finite difference, unit is eV/Angstrom
        call dpeshon(natoms,q,dv) 
        write(*,*)dv

! -4.669620023989296E-005  1.073684985364309E-004 -7.622343645152796E-005
!  1.034455954229685E-004 -1.406232019718345E-005  9.655022337184960E-006
!       (your output should be very close to these ones)
        end
