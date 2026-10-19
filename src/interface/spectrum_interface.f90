!=====================================================================
! spectrum_interface.f90 - the member-seam spectrum producer: the mode
!   tables (and the EBK diatomic constants) the distribution members
!   pull at assembly come from exactly ONE producer, selected by the
!   SPECTRUM_SOURCE key (2026-10-05 spectrum-source legislation,
!   docs/plans/2026-10-05-spectrum-source.md; the container export
!   slot is withdrawn - spectrum_export_put/_put_diatomic and the
!   registration-time export channel no longer exist)
! Design:
!   COMPUTE (default, paradigm-neutral): the internal derivation, lazy
!   on first pull. The probe geometry places every fragment at its
!   buffered internal geometry with COMs on the x axis, consecutive
!   fragments separated by ext/2 + r_asym + ext/2 (ext = the
!   fragment's largest atom distance from its COM - a molecular pair
!   sits as before; a wide substrate still gets atom-scale
!   isolation). ONE central-difference Hessian of the bound container
!   force (container_probe_force - the unit fold stays at its unique
!   site) is cached process-wide; each pull builds and validates ONLY
!   the pulled fragment's block through spectrum_modes (an unconsumed
!   block - e.g. a substrate nobody samples - has no right to abort
!   the process). A mass-weighted cross-block decay gate runs per
!   (pulled, other) pair against the pulled fragment's smallest
!   vibrational eigenvalue: a PES whose inter-fragment interaction has
!   not decayed at the probe carries no isolated-fragment spectrum and
!   named-aborts pointing at MANUAL. Paradigm labels carry no vote
!   here - the old surface blanket ban is withdrawn (a surface PES
!   with a decaying tail derives like any other; a non-decaying one,
!   e.g. an interaction-region NN, fails the gate as a fact about THAT
!   PES, not about the paradigm).
!   The product-species pull spectrum_prod_w mirrors the COMPUTE probe
!   with its own cache: products land on element-matched system slots
!   (the container reads the slots' element identity - a species-level
!   probe must not park an atom on another element's slot), every
!   declared product at its declared internal geometry plus the
!   leftover system atoms, isolated the same ext/2 + r_asym + ext/2
!   way; ONE Hessian of the same bound container force serves every
!   product pull (V1 is COMPUTE-only - the product channel has no
!   MANUAL face; SPECTRUM_SOURCE gates fragments alone).
!   MANUAL: the spectrum-table file (SPECTRUM_FILE) is the only
!   producer - composition-keyed records of frequencies [cm^-1, the
!   spectroscopy face, converted here through wvn_to_w - the unique
!   conversion channel] plus mass-weighted orthonormal columns
!   expressed in the fragment's buffered-geometry Cartesian frame
!   (external tables must be rotated there by the author - no
!   automatic alignment), and optional diatomic constants (the EBK
!   member's channel; no generic derivation exists for
!   anharmonic/rotational constants, so EBK under COMPUTE aborts with
!   the same MANUAL pointer). The reader carries formal checks only
!   (counts, positivity, ascending order, one record per composition);
!   orthonormality and translation-orthogonality stay member-init
!   concerns. A monoatomic fragment needs no record and no force call
!   on either producer (an allocated-empty table).
!=====================================================================
module spectrum_interface
   use config,           only: reactants
   use config_atoms,     only: list_atoms, list_products
   use force_interface,  only: container_probe_force
   use hessian,          only: hessian_cd, force_proc
   use spectrum,         only: spectrum_modes, w_to_wvn, wvn_to_w
   use diag,             only: diag_fatal, diag_cat_output
   implicit none
   private
   public :: spectrum_frag_w, spectrum_frag_wvn, spectrum_frag_diatomic
   public :: spectrum_prod_w

   integer, parameter :: manual_cap = 8      ! manual-record capacity [record]
   integer, parameter :: sym_cap = 32        ! record composition capacity [atom]
   integer, parameter :: tok_cap = 96        ! one table line's token capacity [token]
   real(8), parameter :: r_asym = 12.0d0     ! atom-scale inter-fragment isolation at the probe [A]
   real(8), parameter :: cross_rel = 1.0d-6  ! cross-block gate vs mode curvature [rel]

   ! one mode table: frequencies [rad/fs] ascending + mass-weighted
   ! orthonormal columns (3nat x n_mode)
   type :: tab_t
      real(8), allocatable :: w(:), c(:,:)
   end type tab_t

   ! one manual file record: the composition key + the table + the
   ! optional diatomic constants of the record
   type :: rec_t
      character(len=4), allocatable :: symb(:)
      type(tab_t) :: tab
      logical :: has_diatic = .false.
      real(8) :: w_e = 0.0d0       ! vibrational constant w_e [cm^-1]
      real(8) :: w_ex_e = 0.0d0    ! anharmonic constant w_e*x_e [cm^-1]
      real(8) :: b_rot = 0.0d0     ! rotational constant B_e [cm^-1]
   end type rec_t

   ! derivation cache: the probe geometry + ONE full-system Hessian
   ! (process-wide), plus per-fragment tables built lazily on first pull
   real(8), allocatable, save :: q_probe(:), hess_c(:,:)
   type(tab_t), allocatable, save :: deriv_tbl(:)
   logical, allocatable, save :: deriv_done(:)
   logical, save :: probed = .false.

   ! product-species derivation cache: the product probe geometry + ONE
   ! full-system Hessian (its own state - the product placement is another
   ! probe, never mixed with the fragment cache above), the element-matched
   ! slot table (CSR: product i's system slots are
   ! prod_slot(prod_off(i):prod_off(i+1)-1)), and the per-product tables
   ! built lazily
   real(8), allocatable, save :: q_pprobe(:), hess_p(:,:)
   integer, allocatable, save :: prod_off(:), prod_slot(:)
   type(tab_t), allocatable, save :: prod_tbl(:)
   logical, allocatable, save :: prod_done(:)
   logical, save :: pprobed = .false.

   ! manual records, loaded once on the first MANUAL pull
   type(rec_t), save :: manual_tbl(manual_cap)
   integer, save :: n_manual = 0
   logical, save :: manual_loaded = .false.
contains
   !------------------------------------------------------------------
   ! spectrum_frag_w(frag, w, c) - the member-seam pull: fresh copies
   !                 of fragment frag's mode table from the selected
   !                 producer. w [rad/fs], c (3nat x n)
   !------------------------------------------------------------------
   subroutine spectrum_frag_w(frag, w, c)
      integer, intent(in) :: frag
      real(8), allocatable, intent(out) :: w(:), c(:,:)
      if (.not. allocated(list_atoms%frag)) then
         call stop_si('spectrum_frag_w', 'no list_atoms fragment table - the spectrum seam'// &
                      ' runs after list_atoms assembly (assembly-order error)')
      end if
      if (frag < 1 .or. frag > size(list_atoms%frag)) then
         call stop_si('spectrum_frag_w', 'fragment index out of the atom list range')
      end if
      if (reactants%spectrum_src == 2) then
         if (list_atoms%frag(frag)%nat <= 1) then      ! monoatomic: no record, no staging
            allocate (w(0), c(3, 0))
            return
         end if
         call manual_stage()
         call manual_table(frag, w, c)
      else
         call compute_table(frag, w, c)
      end if
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
   ! spectrum_frag_diatomic(frag, w_e, w_ex_e, b_rot) - the EBK member's
   !                 pull: fragment frag's constants from a MANUAL DIAT
   !                 record. No generic derivation exists for
   !                 anharmonic/rotational constants (they are
   !                 literature data) - COMPUTE aborts with the MANUAL
   !                 pointer, a missing record likewise
   !------------------------------------------------------------------
   subroutine spectrum_frag_diatomic(frag, w_e, w_ex_e, b_rot)
      integer, intent(in) :: frag
      real(8), intent(out) :: w_e, w_ex_e, b_rot
      integer :: k
      if (.not. allocated(list_atoms%frag)) then
         call stop_si('spectrum_frag_diatomic', 'no list_atoms fragment table - the'// &
                      ' spectrum seam runs after list_atoms assembly (assembly-order error)')
      end if
      if (frag < 1 .or. frag > size(list_atoms%frag)) then
         call stop_si('spectrum_frag_diatomic', 'fragment index out of the atom list range')
      end if
      if (reactants%spectrum_src /= 2) then
         call stop_si('spectrum_frag_diatomic', 'the diatomic-constants channel has no'// &
                      ' generic derivation (anharmonic/rotational constants are literature'// &
                      ' data) - provide SPECTRUM_SOURCE=MANUAL with a DIAT record in'// &
                      ' SPECTRUM_FILE')
      end if
      call manual_stage()
      do k = 1, n_manual
         if (symb_seq_eq(list_atoms%symb(list_atoms%frag(frag)%list), manual_tbl(k)%symb)) then
            if (.not. manual_tbl(k)%has_diatic) then
               call stop_si('spectrum_frag_diatomic', 'the matched MANUAL record carries no'// &
                            ' DIAT record - the EBK member needs (w_e, w_e*x_e, B_e) in'// &
                            ' SPECTRUM_FILE')
            end if
            w_e = manual_tbl(k)%w_e
            w_ex_e = manual_tbl(k)%w_ex_e
            b_rot = manual_tbl(k)%b_rot
            return
         end if
      end do
      call stop_si('spectrum_frag_diatomic', 'no MANUAL spectrum record for this fragment'// &
                   ' composition in SPECTRUM_FILE - MANUAL makes the file a requirement'// &
                   ' for every sampled carrier')
   end subroutine spectrum_frag_diatomic

   !==================================================================
   ! MANUAL producer
   !==================================================================

   ! manual_stage() - load the spectrum-table file once (first MANUAL pull)
   subroutine manual_stage()
      if (manual_loaded) return
      manual_loaded = .true.
      if (len_trim(reactants%spectrum_file) == 0) then
         call stop_si('manual_stage', 'SPECTRUM_SOURCE=MANUAL but SPECTRUM_FILE is empty -'// &
                      ' the spectrum-table file is a MANUAL-mode requirement')
      end if
      call read_table_file(trim(reactants%spectrum_file))
   end subroutine manual_stage

   ! manual_table(frag, w, c) - the composition-matched record's table
   subroutine manual_table(frag, w, c)
      integer, intent(in) :: frag
      real(8), allocatable, intent(out) :: w(:), c(:,:)
      integer :: k
      do k = 1, n_manual
         if (symb_seq_eq(list_atoms%symb(list_atoms%frag(frag)%list), manual_tbl(k)%symb)) then
            if (size(manual_tbl(k)%tab%c, 1) /= 3*list_atoms%frag(frag)%nat) then
               call stop_si('manual_table', 'the MANUAL record matches the composition'// &
                            ' but not the atom count of the fragment')
            end if
            w = manual_tbl(k)%tab%w                ! (re)allocation on assignment
            c = manual_tbl(k)%tab%c
            return
         end if
      end do
      call stop_si('manual_table', 'no MANUAL spectrum record for this fragment'// &
                   ' composition in SPECTRUM_FILE='''//trim(reactants%spectrum_file)// &
                   ''' - MANUAL makes the file a requirement for every sampled carrier')
   end subroutine manual_table

   !------------------------------------------------------------------
   ! read_table_file(path) - the spectrum-table reader (formal checks
   !                 only, input-layer posture; orthonormality stays a
   !                 member-init concern). Line-oriented, case-insensitive,
   !                 '#'/'!' whole-line comments, blank lines ignored:
   !                   FRAG s1 s2 ...      record start (composition key; nat>=2)
   !                   NMODE n             1 for nat=2; 3nat-5 or 3nat-6 for nat>=3
   !                   WVN x               one positive frequency [cm^-1] per line,
   !                                      strictly ascending, NMODE lines
   !                   COL j               column j header; the NEXT content line
   !                                      carries exactly 3nat reals
   !                   DIAT we wx be       optional, nat=2 only, once per record
   !                 (the authoritative format spec lives in
   !                 docs/plans/2026-10-05-spectrum-source.md §3)
   !------------------------------------------------------------------
   subroutine read_table_file(path)
      character(len=*), intent(in) :: path
      character(len=1024) :: line
      character(len=32) :: toks(tok_cap)
      character(len=4) :: syms(sym_cap)
      real(8), allocatable :: wvn(:), cmat(:,:)
      real(8) :: d_we, d_wx, d_be
      integer :: u, ios, ln, nt, k, nat, nmod, iwvn, icol, j, ierr
      logical :: ex, in_rec, nmod_set, diat_set

      inquire (file=trim(path), exist=ex)
      if (.not. ex) then
         call stop_si('read_table_file', 'SPECTRUM_FILE '''//trim(path)//''' does not'// &
                      ' exist (the path resolves from the run directory)')
      end if
      open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call stop_si('read_table_file', 'cannot open SPECTRUM_FILE '''//trim(path)//'''')
      end if
      ln = 0
      in_rec = .false.
      nmod = 0
      nmod_set = .false.
      diat_set = .false.
      iwvn = 0
      icol = 0
      nat = 0
      d_we = 0.0d0
      d_wx = 0.0d0
      d_be = 0.0d0
      wvn = [real(8)::]
      cmat = reshape([real(8)::], [1, 0])
      do
         call next_content(u, line, ln, ios)
         if (ios /= 0) exit
         call split_tokens(line, toks, nt)
         select case (to_upper(toks(1)))
         case ('FRAG')
            if (in_rec) then
               call commit_record()
            end if
            nat = nt - 1
            if (nat < 2 .or. nat > sym_cap) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': a FRAG record'// &
                            ' carries '//trim(i2s(nat))//' symbols (2..'//trim(i2s(sym_cap))// &
                            ' expected - a monoatomic fragment needs no record)')
            end if
            do k = 1, nat
               if (len_trim(toks(k + 1)) > 4) then
                  call stop_si('read_table_file', 'line '//trim(i2s(ln))//': symbol '''// &
                               trim(toks(k + 1))//''' exceeds the 4-character width')
               end if
               syms(k) = toks(k + 1)
            end do
            in_rec = .true.
            nmod_set = .false.
            diat_set = .false.
            iwvn = 0
            icol = 0
            deallocate (wvn, cmat, stat=ierr)
            allocate (wvn(0), cmat(1, 0))
         case ('NMODE')
            if (.not. in_rec .or. nmod_set) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': NMODE needs an'// &
                            ' open record without an NMODE row')
            end if
            if (nt /= 2) call stop_si('read_table_file', 'line '//trim(i2s(ln))// &
                                      ': NMODE takes exactly one integer')
            nmod = parse_int(toks(2), ierr)
            if (ierr /= 0) call stop_si('read_table_file', 'line '//trim(i2s(ln))// &
                                        ': bad NMODE integer '''//trim(toks(2))//'''')
            if (nat == 2 .and. nmod /= 1) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': a diatomic'// &
                            ' record carries exactly 1 mode, got '//trim(i2s(nmod)))
            end if
            if (nat >= 3 .and. nmod /= 3*nat - 5 .and. nmod /= 3*nat - 6) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': NMODE '// &
                            trim(i2s(nmod))//' violates the count contract (nat = '// &
                            trim(i2s(nat))//': 3nat-5 or 3nat-6)')
            end if
            nmod_set = .true.
            if (allocated(wvn)) deallocate (wvn)
            if (allocated(cmat)) deallocate (cmat)
            allocate (wvn(nmod), cmat(3*nat, nmod))
            cmat = 0.0d0
         case ('WVN')
            if (.not. in_rec .or. .not. nmod_set .or. iwvn >= nmod) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': a WVN row'// &
                            ' needs an open record with unread frequency slots')
            end if
            if (nt /= 2) call stop_si('read_table_file', 'line '//trim(i2s(ln))// &
                                      ': WVN takes exactly one positive real')
            call parse_real(toks(2), wvn(iwvn + 1), ierr)
            if (ierr /= 0 .or. wvn(iwvn + 1) <= 0.0d0) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': bad WVN value '''// &
                            trim(toks(2))//''' (a positive wavenumber [cm^-1])')
            end if
            if (iwvn >= 1 .and. wvn(iwvn + 1) <= wvn(iwvn)) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': WVN rows must'// &
                            ' be strictly ascending (the derive-path contract)')
            end if
            iwvn = iwvn + 1
         case ('COL')
            if (.not. in_rec .or. .not. nmod_set .or. iwvn /= nmod .or. icol >= nmod) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': a COL row needs'// &
                            ' its record''s frequencies complete and an open column slot')
            end if
            j = parse_int(toks(2), ierr)
            if (ierr /= 0 .or. nt /= 2 .or. j /= icol + 1) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': COL expects the'// &
                            ' next column index '//trim(i2s(icol + 1)))
            end if
            call next_content(u, line, ln, ios)
            if (ios /= 0) then
               call stop_si('read_table_file', 'COL '//trim(i2s(j))//' runs past the end'// &
                            ' of the file (the coefficient line is missing)')
            end if
            call split_tokens(line, toks, nt)
            if (nt /= 3*nat) then
               call stop_si('read_table_file', 'the COL '//trim(i2s(j))//' coefficient'// &
                            ' line carries '//trim(i2s(nt))//' values (exactly '// &
                            trim(i2s(3*nat))//' expected)')
            end if
            do k = 1, 3*nat
               call parse_real(toks(k), cmat(k, j), ierr)
               if (ierr /= 0) then
                  call stop_si('read_table_file', 'the COL '//trim(i2s(j))//' line''s'// &
                               ' value '''//trim(toks(k))//''' is not a real')
               end if
            end do
            icol = icol + 1
         case ('DIAT')
            if (.not. in_rec .or. .not. nmod_set .or. diat_set) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': DIAT needs an'// &
                            ' open record without a DIAT row')
            end if
            if (nat /= 2) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': DIAT is a'// &
                            ' diatomic-record row (this record carries '//trim(i2s(nat))// &
                            ' atoms)')
            end if
            if (nt /= 4) call stop_si('read_table_file', 'line '//trim(i2s(ln))// &
                                      ': DIAT takes exactly three reals (w_e, w_e*x_e, B_e)')
            call parse_real(toks(2), d_we, ierr)
            if (ierr == 0) call parse_real(toks(3), d_wx, ierr)
            if (ierr == 0) call parse_real(toks(4), d_be, ierr)
            if (ierr /= 0 .or. d_we <= 0.0d0 .or. d_wx <= 0.0d0 .or. d_be < 0.0d0) then
               call stop_si('read_table_file', 'line '//trim(i2s(ln))//': bad DIAT values'// &
                            ' (w_e > 0, w_e*x_e > 0, B_e >= 0 [cm^-1])')
            end if
            diat_set = .true.
         case default
            call stop_si('read_table_file', 'line '//trim(i2s(ln))//': unknown'// &
                         ' spectrum-table word '''//trim(toks(1))//'''')
         end select
      end do
      if (in_rec) call commit_record()
      close (u)

   contains

      ! commit_record() - validate completeness and stage the open record
      subroutine commit_record()
         integer :: r, m
         if (.not. nmod_set .or. iwvn /= nmod .or. icol /= nmod) then
            call stop_si('read_table_file', 'the '''//trim(syms(1))//''' record ends'// &
                         ' incomplete (frequencies/columns must cover NMODE = '// &
                         trim(i2s(nmod))//')')
         end if
         do r = 1, n_manual
            if (symb_seq_eq(syms(:nat), manual_tbl(r)%symb)) then
               call stop_si('read_table_file', 'two MANUAL records claim the same atom'// &
                            ' composition - one composition, one record')
            end if
         end do
         if (n_manual >= manual_cap) then
            call stop_si('read_table_file', 'manual-record capacity exceeded')
         end if
         n_manual = n_manual + 1
         allocate (manual_tbl(n_manual)%symb(nat))
         manual_tbl(n_manual)%symb = syms(:nat)
         allocate (manual_tbl(n_manual)%tab%w(nmod), manual_tbl(n_manual)%tab%c(3*nat, nmod))
         do m = 1, nmod
            manual_tbl(n_manual)%tab%w(m) = wvn_to_w(wvn(m))   ! cm^-1 -> rad/fs,
         end do                                                  ! the unique channel
         manual_tbl(n_manual)%tab%c = cmat
         manual_tbl(n_manual)%has_diatic = diat_set
         manual_tbl(n_manual)%w_e = d_we
         manual_tbl(n_manual)%w_ex_e = d_wx
         manual_tbl(n_manual)%b_rot = d_be
      end subroutine commit_record

   end subroutine read_table_file

   !==================================================================
   ! COMPUTE producer
   !==================================================================

   ! compute_table(frag, w, c) - the pulled fragment's derived table
   subroutine compute_table(frag, w, c)
      integer, intent(in) :: frag
      real(8), allocatable, intent(out) :: w(:), c(:,:)
      real(8), allocatable :: hf(:,:), massf(:)
      integer :: a, b
      if (list_atoms%frag(frag)%nat <= 1) then       ! monoatomic: no internal dof,
         allocate (w(0), c(3, 0))                    ! no force call, no gate
         return
      end if
      if (.not. probed) call build_probe()
      if (.not. deriv_done(frag)) then
         associate (fr => list_atoms%frag(frag))
            allocate (hf(3*fr%nat, 3*fr%nat), massf(fr%nat))
            massf = list_atoms%mass(fr%list)
            do a = 1, fr%nat              ! the FULL fragment block (pair
               do b = 1, fr%nat           ! couplings inside the fragment ride)
                  hf(3*a-2:3*a, 3*b-2:3*b) = hess_c(3*fr%list(a)-2:3*fr%list(a), &
                                                  3*fr%list(b)-2:3*fr%list(b))
               end do
            end do
            call spectrum_modes(hf, massf, fr%qz, deriv_tbl(frag)%w, deriv_tbl(frag)%c)
         end associate
         call decay_gate(frag)
         deriv_done(frag) = .true.
      end if
      w = deriv_tbl(frag)%w
      c = deriv_tbl(frag)%c
   end subroutine compute_table

   !------------------------------------------------------------------
   ! build_probe() - the probe geometry (every fragment at its buffered
   !                 internal geometry, COMs on the x axis, consecutive
   !                 fragments separated by ext/2 + r_asym + ext/2) and
   !                 ONE central-difference Hessian of the bound
   !                 container force, cached process-wide. A
   !                 single-fragment system needs no isolation - its
   !                 buffered geometry IS the probe (re-placing the COM
   !                 would serve nothing and would break surface
   !                 placement conventions, e.g. a molecule's height
   !                 above an implicit substrate)
   !------------------------------------------------------------------
   subroutine build_probe()
      integer :: nf, n3, i, j
      real(8), allocatable :: xpos(:)
      real(8) :: com(3), d
      real(8), allocatable :: ext(:), coms(:, :)
      nf = size(list_atoms%frag)
      n3 = 3*size(list_atoms%mass)
      allocate (q_probe(n3), hess_c(n3, n3), deriv_tbl(nf), deriv_done(nf))
      deriv_done = .false.
      if (nf == 1) then                 ! nothing to isolate from: the
         do j = 1, list_atoms%frag(1)%nat   ! buffered geometry as-is
            q_probe(3*list_atoms%frag(1)%list(j)-2:3*list_atoms%frag(1)%list(j)) = &
               list_atoms%frag(1)%qz(3*j-2:3*j)
         end do
      else
         allocate (ext(nf), coms(3, nf), xpos(nf))
         ! 1. per-fragment COM and extent (largest atom distance from the COM)
         do i = 1, nf
            com = 0.0d0
            do j = 1, list_atoms%frag(i)%nat
               com = com + list_atoms%mass(list_atoms%frag(i)%list(j))* &
                           list_atoms%frag(i)%qz(3*j-2:3*j)
            end do
            com = com/list_atoms%frag(i)%mass
            coms(:, i) = com
            ext(i) = 0.0d0
            do j = 1, list_atoms%frag(i)%nat
               d = norm2(list_atoms%frag(i)%qz(3*j-2:3*j) - com)
               ext(i) = max(ext(i), d)
            end do
         end do
         ! 2. place fragments along x with atom-scale isolation
         xpos(1) = 0.0d0
         do i = 2, nf
            xpos(i) = xpos(i - 1) + ext(i - 1)/2 + r_asym + ext(i)/2
         end do
         do i = 1, nf
            do j = 1, list_atoms%frag(i)%nat
               q_probe(3*list_atoms%frag(i)%list(j)-2:3*list_atoms%frag(i)%list(j)) = &
                  list_atoms%frag(i)%qz(3*j-2:3*j) - coms(:, i) + [xpos(i), 0.0d0, 0.0d0]
            end do
         end do
         deallocate (ext, coms, xpos)
      end if
      ! 3. one Hessian of the bound container force (internal-per-A fold)
      call hessian_cd(probe_wrap, q_probe, hess_c)
      probed = .true.
   end subroutine build_probe

   !------------------------------------------------------------------
   ! decay_gate(frag) - the mass-weighted cross-block gate: the pulled
   !                 fragment's Hessian block is its isolated spectrum
   !                 only if the interaction with every other fragment
   !                 has decayed at the probe. The reference curvature
   !                 is the PULLED fragment's smallest vibrational
   !                 eigenvalue (the pair-minimum of the withdrawn
   !                 all-fragments form would force deriving the other
   !                 side's table too - against the per-fragment pull).
   !                 The scan is per atom pair over the WHOLE 3x3
   !                 block, same form as decay_gate_prod
   !------------------------------------------------------------------
   subroutine decay_gate(frag)
      integer, intent(in) :: frag
      integer :: j, a, b, ca, cb, ga, gb, nf
      real(8) :: lam, cross, mca, mcb
      nf = size(list_atoms%frag)
      lam = deriv_tbl(frag)%w(1)**2
      do j = 1, nf
         if (j == frag) cycle
         if (list_atoms%frag(j)%nat <= 1) cycle   ! a monoatomic partner carries
                                                  ! no mode scale (the pair rides)
         cross = 0.0d0
         do a = 1, list_atoms%frag(frag)%nat
            ga = list_atoms%frag(frag)%list(a)
            mca = sqrt(list_atoms%mass(ga))
            do b = 1, list_atoms%frag(j)%nat
               gb = list_atoms%frag(j)%list(b)
               mcb = sqrt(list_atoms%mass(gb))
               do ca = 1, 3
                  do cb = 1, 3
                     cross = max(cross, abs(hess_c(3*ga - 3 + ca, 3*gb - 3 + cb))/(mca*mcb))
                  end do
               end do
            end do
         end do
         if (cross > cross_rel*lam) then
            call stop_si('decay_gate', 'the inter-fragment coupling has not decayed at'// &
                         ' the probe separation - this PES carries no decaying'// &
                         ' isolated-fragment tail; provide the spectrum as'// &
                         ' SPECTRUM_SOURCE=MANUAL + SPECTRUM_FILE')
         end if
      end do
   end subroutine decay_gate

   !==================================================================
   ! COMPUTE producer, product species
   !==================================================================

   !------------------------------------------------------------------
   ! spectrum_prod_w(iprod, w, c) - the product-species pull: declared
   !                 product iprod's mode table (COMPUTE always - V1 has
   !                 no MANUAL face for products). w [rad/fs] ascending,
   !                 c (3nat x n) mass-weighted orthonormal in the
   !                 declared geometry's Cartesian frame
   !------------------------------------------------------------------
   subroutine spectrum_prod_w(iprod, w, c)
      integer, intent(in) :: iprod
      real(8), allocatable, intent(out) :: w(:), c(:,:)
      real(8), allocatable :: hf(:,:), massf(:)
      integer :: a, b, ka, kb
      if (.not. allocated(list_products)) then
         call stop_si('spectrum_prod_w', 'no product table - the product spectrum'// &
                      ' seam runs after list_products assembly (assembly-order error)')
      end if
      if (iprod < 1 .or. iprod > size(list_products)) then
         call stop_si('spectrum_prod_w', 'product index out of the declared-product'// &
                      ' range')
      end if
      if (list_products(iprod)%nat <= 1) then       ! monoatomic: no internal dof,
         allocate (w(0), c(3, 0))                   ! no force call, no gate
         return
      end if
      if (.not. pprobed) call build_prod_probe()
      if (.not. prod_done(iprod)) then
         associate (pr => list_products(iprod))
            allocate (hf(3*pr%nat, 3*pr%nat), massf(pr%nat))
            massf = pr%mass
            do a = 1, pr%nat              ! the FULL product block (pair
               do b = 1, pr%nat           ! couplings inside the product ride)
                  ka = prod_slot(prod_off(iprod)+a-1)
                  kb = prod_slot(prod_off(iprod)+b-1)
                  hf(3*a-2:3*a, 3*b-2:3*b) = hess_p(3*ka-2:3*ka, 3*kb-2:3*kb)
               end do
            end do
            call spectrum_modes(hf, massf, pr%qz, prod_tbl(iprod)%w, prod_tbl(iprod)%c)
         end associate
         call decay_gate_prod(iprod)
         prod_done(iprod) = .true.
      end if
      w = prod_tbl(iprod)%w
      c = prod_tbl(iprod)%c
   end subroutine spectrum_prod_w

   !------------------------------------------------------------------
   ! build_prod_probe() - the product probe geometry and ONE
   !                 central-difference Hessian of the bound container
   !                 force, cached process-wide. Slot assignment is
   !                 element-matched: each product atom (declaration
   !                 order) takes the lowest unoccupied system slot of
   !                 the same element - the container reads the slots'
   !                 element identity (container_pes receives the
   !                 list_atoms slots), so a species-level probe must
   !                 land products on same-element atoms; an atom with
   !                 no unclaimed same-element slot names the product
   !                 and the element (an element-level rearrangement
   !                 check). The unmatched system atoms keep their slots
   !                 and close the line as one leftover group at its q0
   !                 arrangement. Groups sit at their internal
   !                 geometries, COMs on the x axis, separated
   !                 ext/2 + r_asym + ext/2; a lone product covering the
   !                 whole system needs no isolation - its declared
   !                 geometry IS the probe (mirrors build_probe's nf==1:
   !                 re-placing the COM would serve nothing and would
   !                 hand a translation-dependent container other than
   !                 the declared placement)
   !------------------------------------------------------------------
   subroutine build_prod_probe()
      integer :: np, n, n3, ng, i, j, k, ks, ntot, nrest
      integer, allocatable :: rslot(:)
      real(8), allocatable :: rcoord(:), rmass(:), ext(:), coms(:,:), xpos(:)
      logical, allocatable :: taken(:)
      np = size(list_products)
      if (.not. allocated(list_atoms%mass) .or. .not. allocated(list_atoms%symb)) then
         call stop_si('build_prod_probe', 'no atom list - the product probe needs'// &
                      ' the system slots (assembly-order error)')
      end if
      n = size(list_atoms%mass)
      n3 = 3*n
      ! 1. element-matched slot assignment (CSR: product i's slots are
      !    prod_slot(prod_off(i):prod_off(i+1)-1)); greedy lowest-slot per
      !    atom in declaration order - homonuclear symmetry makes the
      !    greedy choice unique and physically sufficient
      ntot = 0
      do i = 1, np
         ntot = ntot + list_products(i)%nat
      end do
      allocate (prod_off(np+1), prod_slot(ntot), taken(n))
      taken = .false.
      prod_off(1) = 1
      do i = 1, np
         do j = 1, list_products(i)%nat
            k = 0
            do ks = 1, n
               if (.not. taken(ks) .and. &
                   ci4(list_atoms%symb(ks)) == ci4(list_products(i)%symb(j))) then
                  k = ks
                  exit
               end if
            end do
            if (k == 0) then
               call stop_si('build_prod_probe', 'product '''//trim(list_products(i)%name)// &
                            ''' declares an atom of element '''//trim(list_products(i)%symb(j))// &
                            ''' but no unclaimed same-element system slot remains -'// &
                            ' products must be element rearrangements of the system'// &
                            ' atoms')
            end if
            prod_slot(prod_off(i)+j-1) = k
            taken(k) = .true.
         end do
         prod_off(i+1) = prod_off(i) + list_products(i)%nat
      end do
      ! 2. the leftover slots (ascending); their q0 arrangement is the
      !    closing group
      nrest = n - ntot
      allocate (rslot(nrest))
      ks = 0
      do k = 1, n
         if (.not. taken(k)) then
            ks = ks + 1
            rslot(ks) = k
         end if
      end do
      ng = np
      if (nrest > 0) ng = np + 1
      allocate (q_pprobe(n3), hess_p(n3, n3), prod_tbl(np), prod_done(np))
      prod_done = .false.
      ! 3. placement: groups at their internal geometries, COMs on the x
      !    axis with atom-scale isolation; the lone full-coverage product
      !    keeps its declared placement as-is
      if (ng == 1) then
         do j = 1, list_products(1)%nat
            k = prod_slot(prod_off(1)+j-1)
            q_pprobe(3*k-2:3*k) = list_products(1)%qz(3*j-2:3*j)
         end do
      else
         allocate (ext(ng), coms(3, ng), xpos(ng))
         do i = 1, np
            call grp_stat(list_products(i)%qz, list_products(i)%mass, coms(:, i), ext(i))
         end do
         if (nrest > 0) then
            allocate (rcoord(3*nrest), rmass(nrest))
            do ks = 1, nrest
               rcoord(3*ks-2:3*ks) = list_atoms%q0(3*rslot(ks)-2:3*rslot(ks))
               rmass(ks) = list_atoms%mass(rslot(ks))
            end do
            call grp_stat(rcoord, rmass, coms(:, ng), ext(ng))
         end if
         xpos(1) = 0.0d0
         do i = 2, ng
            xpos(i) = xpos(i-1) + ext(i-1)/2 + r_asym + ext(i)/2
         end do
         do i = 1, np
            do j = 1, list_products(i)%nat
               k = prod_slot(prod_off(i)+j-1)
               q_pprobe(3*k-2:3*k) = list_products(i)%qz(3*j-2:3*j) - coms(:, i) + &
                                     [xpos(i), 0.0d0, 0.0d0]
            end do
         end do
         do ks = 1, nrest
            k = rslot(ks)
            q_pprobe(3*k-2:3*k) = rcoord(3*ks-2:3*ks) - coms(:, ng) + [xpos(ng), 0.0d0, 0.0d0]
         end do
         deallocate (ext, coms, xpos)
         if (nrest > 0) deallocate (rcoord, rmass)
      end if
      deallocate (taken, rslot)
      ! 4. one Hessian of the bound container force (internal-per-A fold)
      call hessian_cd(probe_wrap, q_pprobe, hess_p)
      pprobed = .true.
   end subroutine build_prod_probe

   !------------------------------------------------------------------
   ! grp_stat(qz, mass, com, ext) - one probe group's placement facts:
   !                 the mass-weighted COM and the largest atom distance
   !                 from it (the isolation extent)
   !------------------------------------------------------------------
   subroutine grp_stat(qz, mass, com, ext)
      real(8), intent(in) :: qz(:), mass(:)
      real(8), intent(out) :: com(3), ext
      real(8) :: wsum, d
      integer :: j
      wsum = sum(mass)
      com = 0.0d0
      do j = 1, size(mass)
         com = com + mass(j)*qz(3*j-2:3*j)
      end do
      com = com/wsum
      ext = 0.0d0
      do j = 1, size(mass)
         d = norm2(qz(3*j-2:3*j) - com)
         ext = max(ext, d)
      end do
   end subroutine grp_stat

   !------------------------------------------------------------------
   ! decay_gate_prod(iprod) - the product version of the cross-block
   !                 gate: product iprod's Hessian block is its isolated
   !                 spectrum only if the interaction with EVERY atom
   !                 outside the product (another product's atoms and
   !                 the leftover atoms alike) has decayed at the probe.
   !                 The reference curvature is the PULLED product's
   !                 smallest vibrational eigenvalue (no partner skip -
   !                 the outside atoms are gated one by one, monoatomic
   !                 or not: the reference never reads the partner)
   !------------------------------------------------------------------
   subroutine decay_gate_prod(iprod)
      integer, intent(in) :: iprod
      integer :: a, ca, cb, ia, k, n
      real(8) :: lam, cross, mca, mcb
      n = size(list_atoms%mass)
      lam = prod_tbl(iprod)%w(1)**2
      do k = 1, n
         if (any(prod_slot(prod_off(iprod):prod_off(iprod+1)-1) == k)) cycle
         cross = 0.0d0
         mcb = sqrt(list_atoms%mass(k))
         do a = 1, list_products(iprod)%nat
            ia = prod_slot(prod_off(iprod)+a-1)
            mca = sqrt(list_atoms%mass(ia))
            do ca = 1, 3
               do cb = 1, 3
                  cross = max(cross, abs(hess_p(3*ia - 3 + ca, 3*k - 3 + cb))/(mca*mcb))
               end do
            end do
         end do
         if (cross > cross_rel*lam) then
            call stop_si('decay_gate_prod', 'the product-surroundings coupling has'// &
                         ' not decayed at the probe separation - this PES carries'// &
                         ' no decaying isolated-product tail (the product channel'// &
                         ' has no MANUAL producer in V1)')
         end if
      end do
   end subroutine decay_gate_prod

   ! probe_wrap(q, f) - the bound container through the interface probe
   subroutine probe_wrap(q, f)
      real(8), intent(in)  :: q(:)
      real(8), intent(out) :: f(:)
      call container_probe_force(q, f)
   end subroutine probe_wrap

   !==================================================================
   ! string helpers (reader mechanics - no system knowledge)
   !==================================================================

   ! next_content(u, line, ln, ios) - the next content line (blank and
   ! '#','!' whole-line comments skipped)
   subroutine next_content(u, line, ln, ios)
      integer, intent(in) :: u
      character(len=*), intent(out) :: line
      integer, intent(inout) :: ln
      integer, intent(out) :: ios
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) return
         ln = ln + 1
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle
         return
      end do
   end subroutine next_content

   ! split_tokens(line, toks, nt) - whitespace split of a content line
   subroutine split_tokens(line, toks, nt)
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: toks(:)
      integer, intent(out) :: nt
      integer :: i, n, ts
      nt = 0
      n = len_trim(line)
      i = 1
      do while (i <= n)
         do while (i <= n .and. (line(i:i) == ' ' .or. line(i:i) == achar(9)))
            i = i + 1
         end do
         if (i > n) exit
         ts = i
         do while (i <= n .and. line(i:i) /= ' ' .and. line(i:i) /= achar(9))
            i = i + 1
         end do
         nt = nt + 1
         if (nt > size(toks)) then
            call stop_si('split_tokens', 'a spectrum-table line carries more tokens'// &
                         ' than the reader capacity')
         end if
         toks(nt) = line(ts:i-1)
      end do
      if (nt == 0) then
         call stop_si('split_tokens', 'an empty content line reached the tokenizer')
      end if
   end subroutine split_tokens

   ! parse_int(s, ierr) / parse_real(s, val, ierr) - token casts
   integer function parse_int(s, ierr)
      character(len=*), intent(in) :: s
      integer, intent(out) :: ierr
      read (s, *, iostat=ierr) parse_int
   end function parse_int

   subroutine parse_real(s, val, ierr)
      character(len=*), intent(in) :: s
      real(8), intent(out) :: val
      integer, intent(out) :: ierr
      read (s, *, iostat=ierr) val
   end subroutine parse_real

   ! to_upper(s) - ASCII uppercase of one token (word selection only)
   pure function to_upper(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, c
      t = s
      do i = 1, len_trim(t)
         c = iachar(t(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) t(i:i) = char(c - 32)
      end do
   end function to_upper

   ! i2s(n) - decimal string of an integer
   pure function i2s(n) result(s)
      integer, intent(in) :: n
      character(len=16) :: s
      write (s, '(i0)') n
   end function i2s

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

   ! stop_si(who, msg) - the named-abort channel
   subroutine stop_si(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': spectrum error: '//trim(msg)
      call diag_fatal(diag_cat_output, 'spectrum_interface.f90', trim(who), &
         'the spectrum table is unusable')
   end subroutine stop_si
end module spectrum_interface
