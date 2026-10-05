!=====================================================================
! spectrum_export.f90 - assembly-time fragment-spectrum buffer: the
!                      container exports first, the generic derivation
!                      fills the rest (lazy, pull-driven)
! Design:
!   The mode tables the distribution members consume are assembly
!   data. Two producers, in priority order: (1) a container export
!   buffered by its reg file at registry_gen_all time (the container
!   folder's private channel - the only producer for PES shapes with
!   no decaying inter-fragment tail, e.g. surface neural networks);
!   (2) the generic derivation, run ONCE on first pull from a member
!   seam: fragments probe at their buffered internal geometry with COMs
!   separated along x by r_asym, one central-difference Hessian of the
!   bound container force (container_probe_force - the fold stays at
!   its unique site), per-fragment blocks through spectrum_modes, and
!   a mass-weighted inter-fragment cross-block gate that named-aborts
!   when the interaction has not decayed (the PES cannot be probed in
!   isolation - the container must export). Surface paradigms carry
!   no generic probe at all (the isolated fragment lies outside the
!   PES domain by construction) and require an export. A monoatomic
!   fragment derives an allocated-empty table without any force call.
!   A separate diatomic-constants record (put_diatomic / frag_diatomic)
!   carries the EBK member's spectroscopic constants - export-only
!   (no generic derivation exists for anharmonic/rotational constants).
!=====================================================================
module spectrum_export
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use force_interface,  only: container_probe_force
   use hessian,       only: hessian_cd, force_proc
   use spectrum,      only: spectrum_modes, w_to_wvn
   implicit none
   private
   public :: spectrum_export_put, spectrum_frag_w, spectrum_frag_wvn
   public :: spectrum_export_put_diatomic, spectrum_frag_diatomic

   integer, parameter :: export_cap = 8     ! buffered-export capacity [table]
   real(8), parameter :: r_asym = 12.0d0    ! inter-fragment COM separation at the probe [A]
   real(8), parameter :: cross_rel = 1.0d-6 ! cross-block gate vs mode curvature [rel]

   ! one mode table: frequencies [rad/(10 fs)] ascending + mass-weighted
   ! orthonormal columns (3nat x n_mode)
   type :: tab_t
      real(8), allocatable :: w(:), c(:,:)
   end type

   ! a container export: the composition key + the table (columns
   ! already aligned to the fragment's buffered geometry)
   type :: export_t
      character(len=4), allocatable :: symb(:)
      type(tab_t) :: tab
   end type

   ! a container's diatomic spectroscopic constants (the EBK member's
   ! channel - Herzberg-type constants of the record, container-
   ! legislated; the mode-table export carries only the harmonic w, so
   ! the anharmonic and rotational constants ride this separate record)
   type :: diab_t
      character(len=4), allocatable :: symb(:)
      real(8) :: w_e = 0.0d0      ! vibrational constant w_e [cm^-1]
      real(8) :: w_ex_e = 0.0d0   ! anharmonic constant w_e*x_e [cm^-1]
      real(8) :: b_rot = 0.0d0    ! rotational constant B_e [cm^-1]
   end type

   type(export_t), save :: export_tbl(export_cap)
   integer, save :: n_export = 0
   type(diab_t), save :: diab_tbl(export_cap)
   integer, save :: n_diab = 0
   type(tab_t), allocatable, save :: deriv_tbl(:)
   logical, save :: derived = .false.
contains
   !------------------------------------------------------------------
   ! spectrum_export_put(symb, w, c) - a container stages its exported
   !                 table (runs at registry_gen_all time, list_atoms in
   !                 place; the composition key matches list_atoms fragments
   !                 case-insensitively by symbol sequence; a duplicate
   !                 composition is a named abort)
   !------------------------------------------------------------------
   subroutine spectrum_export_put(symb, w, c)
      character(len=4), intent(in) :: symb(:)
      real(8), intent(in)          :: w(:), c(:,:)
      integer :: i
      if (size(c, 2) /= size(w)) then
         call stop_stage('spectrum_export_put', 'export table shape mismatch: column count'// &
                         ' does not match the frequency count')
      end if
      do i = 1, n_export
         if (symb_seq_eq(symb, export_tbl(i)%symb)) then
            call stop_stage('spectrum_export_put', 'two exports claim the same atom'// &
                            ' composition - one composition, one export')
         end if
      end do
      if (n_export >= export_cap) then
         call stop_stage('spectrum_export_put', 'export capacity exceeded')
      end if
      n_export = n_export + 1
      allocate (export_tbl(n_export)%symb(size(symb)))
      export_tbl(n_export)%symb = symb
      allocate (export_tbl(n_export)%tab%w(size(w)), export_tbl(n_export)%tab%c(size(c, 1), size(c, 2)))
      export_tbl(n_export)%tab%w = w
      export_tbl(n_export)%tab%c = c
   end subroutine spectrum_export_put

   !------------------------------------------------------------------
   ! spectrum_frag_w(frag, w, c) - the member-seam pull: fresh copies of
   !                 fragment frag's table (export match first, else the
   !                 generic derivation). w [rad/(10 fs)], c (3nat x n)
   !------------------------------------------------------------------
   subroutine spectrum_frag_w(frag, w, c)
      integer, intent(in) :: frag
      real(8), allocatable, intent(out) :: w(:), c(:,:)
      integer :: e
      if (.not. allocated(list_atoms%frag)) then
         call stop_stage('spectrum_frag_w', 'no list_atoms fragment table - the spectrum seam'// &
                         ' runs after list_atoms assembly (assembly-order error)')
      end if
      if (frag < 1 .or. frag > size(list_atoms%frag)) then
         call stop_stage('spectrum_frag_w', 'fragment index out of the atom list range')
      end if
      ! 1. container export (composition match)
      do e = 1, n_export
         if (symb_seq_eq(list_atoms%symb(list_atoms%frag(frag)%list), export_tbl(e)%symb)) then
            if (size(export_tbl(e)%tab%c, 1) /= 3*list_atoms%frag(frag)%nat) then
               call stop_stage('spectrum_frag_w', 'the buffered export matches the composition'// &
                               ' but not the atom count of the fragment')
            end if
            w = export_tbl(e)%tab%w                ! (re)allocation on assignment
            c = export_tbl(e)%tab%c
            return
         end if
      end do
      ! 2. surface paradigms have no generic isolated-fragment probe
      if (reactants%surface_model /= 0) then
         call stop_stage('spectrum_frag_w', 'the surface paradigm has no generic'// &
                         ' isolated-fragment spectrum probe - the container must stage'// &
                         ' its own spectrum export')
      end if
      ! 3. generic derivation, once per process
      if (.not. derived) call derive_all()
      w = deriv_tbl(frag)%w
      c = deriv_tbl(frag)%c
   end subroutine spectrum_frag_w

   !------------------------------------------------------------------
   ! spectrum_frag_wvn(frag, nu) - the same pull in wavenumbers [cm^-1]
   !                 (the boltzmann member's table contract)
   !------------------------------------------------------------------
   subroutine spectrum_frag_wvn(frag, nu)
      integer, intent(in) :: frag
      real(8), allocatable, intent(out) :: nu(:)
      real(8), allocatable :: wt(:), ct(:,:)
      integer :: k
      call spectrum_frag_w(frag, wt, ct)
      allocate (nu(size(wt)))
      do k = 1, size(wt)
         nu(k) = w_to_wvn(wt(k))
      end do
   end subroutine spectrum_frag_wvn

   !------------------------------------------------------------------
   ! spectrum_export_put_diatomic(symb, w_e, w_ex_e, b_rot) - a container
   !                 stages its diatomic spectroscopic constants (runs at
   !                 registry_gen_all time; composition-keyed like the mode
   !                 export; one composition, one record)
   !------------------------------------------------------------------
   subroutine spectrum_export_put_diatomic(symb, w_e, w_ex_e, b_rot)
      character(len=4), intent(in) :: symb(:)
      real(8), intent(in) :: w_e, w_ex_e, b_rot
      integer :: i
      do i = 1, n_diab
         if (symb_seq_eq(symb, diab_tbl(i)%symb)) then
            call stop_stage('spectrum_export_put_diatomic', 'two diatomic-constant'// &
                            ' records claim the same atom composition - one'// &
                            ' composition, one record')
         end if
      end do
      if (n_diab >= export_cap) then
         call stop_stage('spectrum_export_put_diatomic', 'record capacity exceeded')
      end if
      n_diab = n_diab + 1
      allocate (diab_tbl(n_diab)%symb(size(symb)))
      diab_tbl(n_diab)%symb = symb
      diab_tbl(n_diab)%w_e = w_e
      diab_tbl(n_diab)%w_ex_e = w_ex_e
      diab_tbl(n_diab)%b_rot = b_rot
   end subroutine spectrum_export_put_diatomic

   !------------------------------------------------------------------
   ! spectrum_frag_diatomic(frag, w_e, w_ex_e, b_rot) - the EBK member's
   !                 pull: fragment frag's composition-matched constants.
   !                 Export-only channel (no generic derivation exists
   !                 for anharmonic/rotational constants) - a miss is a
   !                 named abort telling the container to legislate them
   !------------------------------------------------------------------
   subroutine spectrum_frag_diatomic(frag, w_e, w_ex_e, b_rot)
      integer, intent(in) :: frag
      real(8), intent(out) :: w_e, w_ex_e, b_rot
      integer :: e
      if (.not. allocated(list_atoms%frag)) then
         call stop_stage('spectrum_frag_diatomic', 'no list_atoms fragment table - the'// &
                         ' spectrum seam runs after list_atoms assembly (assembly-order error)')
      end if
      if (frag < 1 .or. frag > size(list_atoms%frag)) then
         call stop_stage('spectrum_frag_diatomic', 'fragment index out of the atom list range')
      end if
      do e = 1, n_diab
         if (symb_seq_eq(list_atoms%symb(list_atoms%frag(frag)%list), diab_tbl(e)%symb)) then
            w_e = diab_tbl(e)%w_e
            w_ex_e = diab_tbl(e)%w_ex_e
            b_rot = diab_tbl(e)%b_rot
            return
         end if
      end do
      call stop_stage('spectrum_frag_diatomic', 'no diatomic constants staged for this'// &
                      ' fragment composition - the container must legislate and export'// &
                      ' them (spectrum_export_put_diatomic)')
   end subroutine spectrum_frag_diatomic

   !------------------------------------------------------------------
   ! derive_all() - the generic derivation over every list_atoms fragment:
   !                 probe geometry (buffered internal geometry, COMs
   !                 along x at r_asym spacing), one Hessian of the
   !                 bound container force, per-fragment mode tables,
   !                 cross-block decay gate
   !------------------------------------------------------------------
   subroutine derive_all()
      real(8) :: q_probe(3*size(list_atoms%mass)), hess(3*size(list_atoms%mass), 3*size(list_atoms%mass))
      real(8), allocatable :: hf(:,:), massf(:)
      integer :: i, j, a, b, nf, n3
      real(8) :: com(3), lam_min, cross, mca, mcb
      nf = size(list_atoms%frag)
      n3 = 3*size(list_atoms%mass)
      ! 1. probe geometry: each fragment at its buffered internal geometry,
      !    COM on the x axis at r_asym spacing
      do i = 1, nf
         associate (fr => list_atoms%frag(i))
            com = 0.0d0
            do j = 1, fr%nat
               com = com + list_atoms%mass(fr%list(j))*fr%qz(3*j-2:3*j)
            end do
            com = com/fr%mass
            do j = 1, fr%nat
               q_probe(3*fr%list(j)-2:3*fr%list(j)) = fr%qz(3*j-2:3*j) - com &
                                                       + [real(8):: (i-1)*r_asym, 0.0d0, 0.0d0]
            end do
         end associate
      end do
      ! 2. one Hessian of the bound container force (internal-per-A fold)
      call hessian_cd(probe_wrap, q_probe, hess)
      ! 3. per-fragment tables from the diagonal blocks
      allocate (deriv_tbl(nf))
      do i = 1, nf
         associate (fr => list_atoms%frag(i))
            allocate (hf(3*fr%nat, 3*fr%nat), massf(fr%nat))
            massf = list_atoms%mass(fr%list)
            do a = 1, fr%nat              ! the FULL fragment block (pair
               do b = 1, fr%nat           ! couplings inside the fragment ride
                  hf(3*a-2:3*a, 3*b-2:3*b) = hess(3*fr%list(a)-2:3*fr%list(a), &
                                                  3*fr%list(b)-2:3*fr%list(b))
               end do
            end do
            call spectrum_modes(hf, massf, fr%qz, deriv_tbl(i)%w, deriv_tbl(i)%c)
            deallocate (hf, massf)
         end associate
      end do
      ! 4. cross-block decay gate (mass-weighted): a non-decayed
      !    inter-fragment interaction fails loudly - the container's
      !    PES cannot be probed in isolation and must export
      do i = 1, nf
         do j = i + 1, nf
            if (size(deriv_tbl(i)%w) == 0 .or. size(deriv_tbl(j)%w) == 0) cycle
            lam_min = min(deriv_tbl(i)%w(1), deriv_tbl(j)%w(1))**2
            cross = 0.0d0
            do a = 1, 3*list_atoms%frag(i)%nat
               do b = 1, 3*list_atoms%frag(j)%nat
                  mca = sqrt(list_atoms%mass(list_atoms%frag(i)%list((a + 2)/3)))
                  mcb = sqrt(list_atoms%mass(list_atoms%frag(j)%list((b + 2)/3)))
                  cross = max(cross, abs(hess(3*list_atoms%frag(i)%list((a + 2)/3) - 3 + a, &
                                              3*list_atoms%frag(j)%list((b + 2)/3) - 3 + b))/(mca*mcb))
               end do
            end do
            if (cross > cross_rel*lam_min) then
               call stop_stage('derive_all', 'the inter-fragment coupling has not decayed at'// &
                               ' the probe separation - this PES carries no decaying tail and'// &
                               ' the container must stage its own spectrum export')
            end if
         end do
      end do
      derived = .true.
   end subroutine derive_all

   ! probe_wrap(q, f) - the bound container through the interface probe
   subroutine probe_wrap(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      call container_probe_force(q, f)
   end subroutine probe_wrap

   ! symb_seq_eq(a, b) - case-insensitive symbol-sequence equality
   pure logical function symb_seq_eq(a, b)
      character(len=4), intent(in) :: a(:), b(:)
      integer :: i
      symb_seq_eq = (size(a) == size(b))
      if (.not. symb_seq_eq) return
      do i = 1, size(a)
         if (ci4(a(i)) /= ci4(b(i))) then
            symb_seq_eq = .false.
            return
         end if
      end do
   end function symb_seq_eq

   ! ci4(s) - ASCII lowercase of a 4-character symbol
   pure function ci4(s) result(t)
      character(len=4), intent(in) :: s
      character(len=4) :: t
      integer :: i, c
      t = s
      do i = 1, len(t)
         c = iachar(t(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) t(i:i) = char(c + 32)
      end do
   end function ci4

   ! stop_stage(who, msg) - the named-abort channel
   subroutine stop_stage(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': spectrum buffer error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the spectrum table is unusable)'
      stop 1
   end subroutine stop_stage
end module spectrum_export
