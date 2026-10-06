!=====================================================================
! linalg.f90 - symmetric-matrix eigendecomposition (Givens-Householder +
!              shifted QL)
! Design:
!   eig_sym(a, d) reads a by the lower triangle, destroys it, and returns
!   the eigenvalues algebraically ascending in d with the eigenvectors
!   written back into a BY COLUMN. The convergence threshold is an
!   internal default (rho = 1e-13).
!=====================================================================
module linalg
   implicit none
   private
   public :: eig_sym
contains
   !------------------------------------------------------------------
   ! eig_sym(a, d) - symmetric-matrix eigenvalues (ascending) and eigenvectors
   !                 (a destroyed, vectors written back by column)
   !------------------------------------------------------------------
   subroutine eig_sym(a, d)
      real(8), intent(inout) :: a(:,:)  ! symmetric matrix n x n (lower triangle on input,
                                        ! destroyed; on exit eigenvectors BY COLUMN,
                                        ! ascending to match d)
      real(8), intent(out)   :: d(:)    ! eigenvalues (n of them, algebraically ascending)
      real(8), parameter :: rho = 1.0d-13    ! convergence threshold
      integer :: n, n1, n2, nr, i, j, k, m, npas, il, nv, np, lv, nt, itemp, i1
      logical :: entry_400
      integer, allocatable :: iposv(:), ivpos(:), iord(:)
      real(8), allocatable :: wk(:), betasq(:), gam(:), beta(:), pk(:), qk(:), eig(:), vec(:,:)
      real(8) :: rhosq, shift, b, s, sgn, sqrts, dd, temp, wtaw, sum, qj, wj, cosa, cosap, &
                 sina, sina2, g, pp, ppbs, ppbr, dia, u, r1, r2, r12, dif, a2
      n = size(d)
      if (n == 0) return
      if (n == 1) then
         ! the trivial eigenpair: a 1x1 matrix is its own eigen-decomposition (the
         ! general path below reads a(n,n-1)/beta(n-1) - out of range at n = 1)
         d(1) = a(1,1)
         a(1,1) = 1.0d0
         return
      end if
      allocate (wk(n), betasq(n), gam(n), beta(n), pk(n), qk(n), eig(n), vec(n,n), &
                iposv(n), ivpos(n), iord(n))
      rhosq = rho*rho
      n1 = n - 1
      n2 = n - 2
      gam(1) = a(1,1)
      ! ---- Givens-Householder tridiagonalization ----
      if (n2 > 0) then
         do nr = 1, n2
            b = a(nr+1,nr)
            s = 0.0d0
            do i = nr, n2
               s = s + a(i+2,nr)**2
            end do
            ! prepare for possible bypass of the transformation
            a(nr+1,nr) = 0.0d0
            if (s > 0.0d0) then
               s = s + b*b
               sgn = 1.0d0
               if (b < 0.0d0) sgn = -1.0d0
               sqrts = sqrt(s)
               dd = sgn/(sqrts + sqrts)
               temp = sqrt(0.5d0 + b*dd)
               wk(nr) = temp
               a(nr+1,nr) = temp
               dd = dd/temp
               b = -sgn*sqrts
               ! dd is the factor of proportionality; compute and save the W vector
               do i = nr, n2
                  temp = dd*a(i+2,nr)
                  wk(i+1) = temp
                  a(i+2,nr) = temp
               end do
               ! premultiply W by A to obtain P; accumulate the scalar k = W^T A W
               wtaw = 0.0d0
               do i = nr, n1
                  sum = 0.0d0
                  do j = nr, i
                     sum = sum + a(i+1,j+1)*wk(j)
                  end do
                  i1 = i + 1
                  if (n1 - i1 >= 0) then
                     do j = i1, n1
                        sum = sum + a(j+1,i+1)*wk(j)
                     end do
                  end if
                  pk(i) = sum
                  wtaw = wtaw + sum*wk(i)
               end do
               ! Q vector; form PAP
               do i = nr, n1
                  qk(i) = pk(i) - wtaw*wk(i)
               end do
               do j = nr, n1
                  qj = qk(j)
                  wj = wk(j)
                  do i = j, n1
                     a(i+1,j+1) = a(i+1,j+1) - 2.0d0*(wk(i)*qj + wj*qk(i))
                  end do
               end do
            end if
            beta(nr) = b          ! written also on the S=0 bypass path
            betasq(nr) = b*b
            gam(nr+1) = a(nr+1,nr+1)
         end do
      end if
      b = a(n,n-1)
      beta(n-1) = b
      betasq(n-1) = b*b
      gam(n) = a(n,n)
      betasq(n) = 0.0d0
      ! ---- adjoint an identity matrix; shifted QL iteration ----
      do i = 1, n
         do j = 1, n
            vec(i,j) = 0.0d0
         end do
         vec(i,i) = 1.0d0
      end do
      m = n
      sum = 0.0d0
      npas = 1
      ! Control-flow note: the loop is entered at the deflate block FIRST (m = n deflated
      ! without a value taken - eig(n+1) does not exist); entry_400 = .true. means arrived
      ! at the deflate block, .false. means arrived at the convergence tail of a sweep
      entry_400 = .true.
      ql_outer: do
         if (entry_400) then
            ! ---- deflate one position ----
            beta(m) = 0.0d0
            betasq(m) = 0.0d0
            m = m - 1
            if (m == 0) exit ql_outer
            if (betasq(m) <= rhosq) then
               eig(m+1) = gam(m+1) + sum          ! converged
               entry_400 = .true.
               cycle ql_outer
            end if
         else
            ! ---- convergence test after a sweep ----
            npas = npas + 1
            if (betasq(m) <= rhosq) then
               eig(m+1) = gam(m+1) + sum          ! converged (fall-through)
               entry_400 = .true.
               cycle ql_outer
            end if
         end if
         ! ---- shift = the corner-2x2 root nearest gam(m+1) ----
         a2 = gam(m+1)
         r2 = 0.5d0*a2
         r1 = 0.5d0*gam(m)
         r12 = r1 + r2
         dif = r1 - r2
         temp = sqrt(dif*dif + betasq(m))
         r1 = r12 + temp
         r2 = r12 - temp
         dif = abs(a2 - r1) - abs(a2 - r2)
         if (dif >= 0.0d0) then
            shift = r2
         else
            shift = r1
         end if
         ! ---- one shifted QL sweep over j = 1..m ----
         sum = sum + shift
         cosa = 1.0d0
         g = gam(1) - shift
         pp = g
         ppbs = pp*pp + betasq(1)
         ppbr = sqrt(ppbs)
         do j = 1, m
            cosap = cosa
            if (ppbs == 0.0d0) then
               sina = 0.0d0
               sina2 = 0.0d0
               cosa = 1.0d0
            else
               sina = beta(j)/ppbr
               sina2 = betasq(j)/ppbs
               cosa = pp/ppbr
               ! postmultiply the identity by P-transpose
               nt = j + npas
               if (nt > n) nt = n
               do i = 1, nt
                  temp = cosa*vec(j,i) + sina*vec(j+1,i)
                  vec(j+1,i) = -sina*vec(j,i) + cosa*vec(j+1,i)
                  vec(j,i) = temp
               end do
            end if
            dia = gam(j+1) - shift
            u = sina2*(g + dia)
            gam(j) = g + u
            g = dia - u
            pp = dia*cosa - sina*cosap*beta(j)
            if (j == m) then
               beta(j) = sina*pp
               betasq(j) = sina2*pp*pp
               exit
            end if
            ppbs = pp*pp + betasq(j+1)
            ppbr = sqrt(ppbs)
            beta(j) = sina*ppbr
            betasq(j) = sina2*ppbs
         end do
         gam(m+1) = g                                ! sweep tail
         entry_400 = .false.
         cycle ql_outer
      end do ql_outer
      eig(1) = gam(1) + sum
      ! ---- transposition sort of the eigenvalues + vector reordering ----
      do j = 1, n
         iposv(j) = j
         ivpos(j) = j
         iord(j) = j
      end do
      m = n
      do
         m = m - 1
         if (m == 0) exit
         do j = 1, m
            if (eig(j) > eig(j+1)) then
               temp = eig(j)
               eig(j) = eig(j+1)
               eig(j+1) = temp
               itemp = iord(j)
               iord(j) = iord(j+1)
               iord(j+1) = itemp
            end if
         end do
      end do
      if (n1 /= 0) then
         do il = 1, n1
            nv = iord(il)
            np = iposv(nv)
            if (np /= il) then
               lv = ivpos(il)
               ivpos(np) = lv
               iposv(lv) = np
               do i = 1, n
                  temp = vec(il,i)
                  vec(il,i) = vec(np,i)
                  vec(np,i) = temp
               end do
            end if
         end do
      end if
      ! ---- back-transform: accumulate the Householder transforms ----
      do i = 1, n
         k = n1
         do
            k = k - 1
            if (k <= 0) exit
            sum = 0.0d0
            do j = k, n1
               sum = sum + vec(i,j+1)*a(j+1,k)
            end do
            sum = sum + sum
            do j = k, n1
               vec(i,j+1) = vec(i,j+1) - sum*a(j+1,k)
            end do
         end do
      end do
      ! ---- write back: eigenvalues to d, vectors BY COLUMN to a (vec rows = vectors) ----
      d(1:n) = eig(1:n)
      do j = 1, n
         do i = 1, n
            a(i,j) = vec(j,i)
         end do
      end do
   end subroutine eig_sym
end module linalg
