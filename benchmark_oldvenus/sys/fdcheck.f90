program fdcheck
   implicit double precision (a-h, o-z)
   parameter (NDA = 160, NDA3 = NDA*3)
   common/QPDOT/Q(NDA3), PDOT(NDA3), FCOEF(NDA3, NDA3)
   common/CONSTN/C1, C2, C3, C4, C5, C6, C7, PI, HALFPI, TWOPI
   dimension q0(9), fd(9), rr(3), dvv(3)
   data q0 / 0.2d0, -0.1d0, 3.0d0,  0.9d0, 0.05d0, 0.4d0,  -0.3d0, 0.5d0, -0.8d0 /
   bohr = 0.529177210903d0
   rr(1) = sqrt((q0(4)-q0(1))**2+(q0(5)-q0(2))**2+(q0(6)-q0(3))**2)/bohr
   rr(2) = sqrt((q0(7)-q0(4))**2+(q0(8)-q0(5))**2+(q0(9)-q0(6))**2)/bohr
   rr(3) = sqrt((q0(7)-q0(1))**2+(q0(8)-q0(2))**2+(q0(9)-q0(3))**2)/bohr
   write (*, '("pairs bohr:", 3f9.4)') rr
   id = 1
   call bkmp2(rr, vv, dvv, id)
   write (*, '("engine V [Eh] =", f14.8, " dV [Eh/bohr]:", 3f12.6)') vv, dvv
   do k = 1, 9
      Q(k) = q0(k)
   end do
   C1 = 0.04184d0
   call POT0(3, v0)
   write (*, '("POT0 V [internal] =", f14.8)') v0
   call probe_pairs
   call DPESHON(3)
   write (*, '("DPESHON PDOT:", 3f13.6)') PDOT(1), PDOT(2), PDOT(3)
   dd = 1.0d-5
   call POT0(3, v0)
   do k = 1, 9
      Q(k) = q0(k) + dd
      call POT0(3, v1)
      Q(k) = q0(k) - dd
      call POT0(3, v2)
      Q(k) = q0(k)
      fd(k) = -(v1 - v2)/(2.0d0*dd)
   end do
   call DPESHON(3)
   write (*, '(a)') "  k      FD(-dV/dq)     DPESHON        ratio"
   do k = 1, 9
      write (*, '(i3, 2f14.6, f10.5)') k, fd(k), PDOT(k), PDOT(k)/fd(k)
   end do
end program fdcheck

subroutine probe_pairs
   implicit double precision (a-h, o-z)
   parameter (NDA = 160, NDA3 = NDA*3)
   dimension r(3)
   common/QPDOT/Q(NDA3), PDOT(NDA3), FCOEF(NDA3, NDA3)
   common/CONSTN/C1, C2, C3, C4, C5, C6, C7, PI, HALFPI, TWOPI
   call h3_pairs(r)
   write (*, '("probe h3_pairs r:", 3f9.4)') r
end subroutine probe_pairs
