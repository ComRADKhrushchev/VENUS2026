!=====================================================================
! geom_init.f90 - geometry-initialization prober: separate a fragment's
!   COM to a PES-validated non-interacting distance along +x
! Design:
!   Stateless utility over a force-probe callback (the hessian_cd pattern):
!   geom_sep_probe places the fragment's atoms with their COM at trial
!   distances d on a ladder (x1.5 growth from d0) and measures the
!   fragment's PERTURBATION OF THE ENVIRONMENT - the probe force on every
!   NON-fragment atom with the fragment at [d, 0, 0] minus the same with
!   the fragment parked far aside. Global container fields (a harmonic
!   well acting on every atom) cancel in the difference, so the measure is
!   exactly the fragment-environment interaction. The first d whose max
!   perturbation component falls under tol wins; the placement stays
!   applied and d returns. A system with no other atoms passes trivially
!   (nothing to protect; the fragment stays where it is). A named abort
!   fires when the ladder cap is reached without decay - the PES does not
!   separate, the suspect is the surface itself (the cheap PES sanity
!   check: an interaction that never decays with distance is a wrong or
!   non-physical delivery). Units: q [Angstrom]; f [internal] (numerically
!   ~ eV/Angstrom); tol [internal].
!=====================================================================
module geom_init
   implicit none
   private
   public :: geom_sep_probe

   ! force-probe callback (the container probe wrapped caller-side; the
   ! same interface hessian_cd consumes)
   abstract interface
      subroutine force_proc_i(q, f)
         real(8), intent(in) :: q(:)
         real(8), intent(out) :: f(:)
      end subroutine force_proc_i
   end interface

   real(8), parameter :: grow = 1.5d0     ! ladder growth factor [-]
   integer, parameter :: max_step = 20    ! ladder cap [step] (d0*grow^19 ~ 2200x d0)
contains
   !------------------------------------------------------------------
   ! geom_sep_probe(fproc, mass, list, sta, d0, tol, d_out) - place the
   !   list-fragment at the first ladder distance whose interaction force
   !   (whole-configuration probe minus the isolated-fragment probe, max
   !   component over the fragment atoms) stays under tol; the placement
   !   is applied to sta (internal geometry preserved, COM at [d, 0, 0])
   !------------------------------------------------------------------
   subroutine geom_sep_probe(fproc, mass, list, q, d0, tol, d_out)
      procedure(force_proc_i) :: fproc       ! the force probe (whole system)
      real(8), intent(in) :: mass(:)         ! per-atom masses [amu]
      integer, intent(in) :: list(:)         ! the fragment's atom list [global atom index]
      real(8), intent(inout) :: q(:)         ! whole-system coordinates [Angstrom] (the
                                             ! fragment is translated in place)
      real(8), intent(in) :: d0              ! ladder start distance [Angstrom]
      real(8), intent(in) :: tol             ! interaction-force tolerance [internal ~ eV/A]
      real(8), intent(out) :: d_out          ! the validated distance [Angstrom]

      real(8), allocatable :: f_ref(:), f_try(:)
      real(8) :: com(3), shift(3), d, fmax
      integer :: i, a, n3, nat, nother
      nat = size(mass)
      n3 = 3*nat
      nother = nat - size(list)
      if (nother == 0) then
         d_out = d0                        ! nothing to protect - trivial pass
         return
      end if
      allocate (f_ref(n3), f_try(n3))

      ! the reference: environment force with the fragment parked far aside
      call shift_frag(q, list, [ -1.0d3, 0.0d0, 0.0d0 ])
      call fproc(q, f_ref)
      call shift_frag(q, list, [ 1.0d3, 0.0d0, 0.0d0 ])

      ! ladder scan: fragment COM at [d, 0, 0], measure the perturbation of
      ! every NON-fragment atom's force
      d = d0
      do i = 1, max_step
         com = frag_com(q, mass, list)
         shift = [ d, 0.0d0, 0.0d0 ] - com
         call shift_frag(q, list, shift)
         call fproc(q, f_try)
         fmax = 0.0d0
         do a = 1, nat
            if (any(list == a)) cycle
            fmax = max(fmax, maxval(abs(f_try(3*a-2:3*a) - f_ref(3*a-2:3*a))))
         end do
         if (fmax <= tol) then
            d_out = d
            deallocate (f_ref, f_try)
            return
         end if
         ! revert and grow (the internal geometry never changed; only the
         ! translation is rewound)
         call shift_frag(q, list, -shift)
         d = d*grow
      end do
      call stop_geoinit(d0, tol)
   end subroutine geom_sep_probe

   ! frag_com - the fragment's mass-weighted COM from the flat q array
   function frag_com(q, mass, list) result(com)
      real(8), intent(in) :: q(:), mass(:)
      integer, intent(in) :: list(:)
      real(8) :: com(3), wt
      integer :: i, i3
      wt = sum(mass(list))
      com = 0.0d0
      do i = 1, size(list)
         i3 = 3*(list(i) - 1)
         com = com + mass(list(i))*q(i3+1:i3+3)
      end do
      com = com/wt
   end function frag_com

   ! shift_frag - translate the list atoms by shift
   subroutine shift_frag(q, list, shift)
      real(8), intent(inout) :: q(:)
      integer, intent(in) :: list(:)
      real(8), intent(in) :: shift(3)
      integer :: i, i3
      do i = 1, size(list)
         i3 = 3*(list(i) - 1)
         q(i3+1:i3+3) = q(i3+1:i3+3) + shift
      end do
   end subroutine shift_frag

   ! stop_geoinit - the named abort: the interaction never decayed
   subroutine stop_geoinit(d0, tol)
      real(8), intent(in) :: d0, tol
      write (0, '(a,es12.4,a,es12.4,a)') 'geom_sep_probe: the interaction force does not '// &
         'decay within the ladder (start ', d0, ' A, growth 1.5x, 20 steps) under the '// &
         'tolerance ', tol, ' - the suspect is the PES itself (a non-decaying '// &
         'separation dependence is a wrong or non-physical delivery)'
      write (0, '(a)') 'geom_sep_probe: fatal (the geometry initialization is not usable)'
      stop 1
   end subroutine stop_geoinit
end module geom_init
