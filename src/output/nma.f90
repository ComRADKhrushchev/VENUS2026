!=====================================================================
! nma.f90 - final-state normal-mode analyzer + generic channel labeler
! Design:
!   Two phases: nma_derive (driver assembly 9c) collects and validates
!   the units, pulls the mode tables, sets thresholds - every data-face
!   abort fires there, before any trajectory burns; nma_setup (the
!   statistics head) opens the files. Units = reactant fragments (fixed
!   atom sets = the elastic channel) UNION declared products (species
!   semantics: per trajectory the best-fitting element-compatible atom
!   subset). Surface paradigm (2026-10-15 ruling): fragment units = the
!   PROJECTILE fragment alone (substrate carriers never units, never
!   pulled); 'reactants' = all reactant units formed - zero reactant
!   units (monoatomic projectile) is vacuously true, a bare reflection
!   labels 'reactants'. Per unit per trajectory: COM strip -> quat_fit
!   -> rotate p -> mode projections (the inverse of samp_normalmode's
!   realize folds) -> n_k = E_k/(hbar*w) - 1/2, one row per trajectory
!   (no gate - the drift column filters). Labels: formed = best rms <
!   tol AND subset != any fragment set; else reactants/unassigned,
!   stems '+'; tol (default 0.4*min pair distance) is the label face's
!   alone.
!=====================================================================
module nma
   use consts,             only: ev_dyn, hbar_evfs
   use state,              only: state_t
   use config,             only: reactants
   use config_atoms,       only: list_atoms, list_products
   use fitrot,             only: quat_fit, eckart_rot, rot3
   use spectrum_interface, only: spectrum_frag_w, spectrum_prod_w
   use spectrum,           only: w_to_wvn
   use diag,               only: diag_fatal, diag_cat_output, diag_kv, diag_eol
   implicit none
   private
   public :: nma_arm_from_input, nma_derive, nma_setup, nma_traj, nma_summary, nma_armed

   integer, parameter :: subset_cap = 1000   ! product subset enumeration cap [subset]
   integer, parameter :: label_cap = 64      ! distinct channel-label capacity [label]
   character(len=*), parameter :: row_fmt = &
      '(i0,1x,es24.16,1x,es24.16,1x,a,1x,es24.16,3(1x,es24.16),999(1x,es24.16))'

   ! One analysis unit: a reactant fragment (identity semantics - the atom
   ! set is fixed) or a declared product (species semantics - the atom set
   ! resolves per trajectory as the best-fitting element-compatible subset)
   type :: unit_t
      character(len=64) :: name = ''             ! unit name (fragment/product file stem) [-]
      integer :: nat = 0                         ! unit atom count [count]
      logical :: is_product = .false.            ! species vs identity semantics [-]
      integer :: frag = 0                        ! reactant units: the list_atoms fragment row [-]
      integer :: iprod = 0                       ! product units: the list_products row [-]
      character(len=4), allocatable :: symb(:)   ! unit atom symbols [-]
      integer, allocatable :: list(:)            ! reactant units: the fixed atom set; product: unset
      real(8), allocatable :: qz(:)              ! COM-centered reference [A], 3*nat
      real(8), allocatable :: mass(:)            ! unit atom masses [amu]
      real(8), allocatable :: w(:)               ! mode frequencies [rad/fs]
      real(8), allocatable :: c(:,:)             ! mass-weighted columns (3nat x n_mode) [-]
      integer :: n_sub = 0                       ! candidate atom-set count [subset]
      integer, allocatable :: subs(:,:)          ! candidate atom sets (nat x n_sub), ascending
      real(8) :: tol = 0.0d0                     ! label threshold [A]
      integer :: u = 0                           ! nma_<name>.dat unit [-]
      integer :: n_rows = 0                      ! written rows [count]
      real(8) :: e_sum = 0.0d0, rms_max = 0.0d0, eck_max = 0.0d0
   end type unit_t
   type(unit_t), allocatable :: units(:)         ! the analysis-unit table []
   logical :: armed = .false., derived = .false., set_up = .false.
   real(8) :: fit_tol_in = 0.0d0                 ! NMA_FIT_TOL override [A] (0 = per-unit default)
   integer :: u_label = 0                        ! nma_label.dat unit [-]
   character(len=64) :: outc_names(label_cap) = ''  ! channel-label counting table [-]
   integer :: outc_cnt(label_cap) = 0
   integer :: n_outc = 0

   ! One product element pool: the system atoms of one element plus how
   ! many the product declaration draws from it
   type :: pool_t
      character(len=4) :: sym = ''
      integer, allocatable :: idx(:)
      integer :: k = 0
   end type pool_t
contains
   !------------------------------------------------------------------
   ! nma_arm_from_input(on, fit_tol) - driver-side arming: store the two
   !                      key values (fit_tol 0 = the per-unit auto default)
   !------------------------------------------------------------------
   subroutine nma_arm_from_input(on, fit_tol)
      logical, intent(in) :: on
      real(8), intent(in) :: fit_tol
      if (fit_tol < 0.0d0) then
         write (0, '(a,es12.4)') 'PRODUCT_NMA: NMA_FIT_TOL must be >= 0 A (0 = the'// &
            ' per-unit 0.4*r_min default), got ', fit_tol
         call diag_fatal(diag_cat_output, 'nma.f90', 'PRODUCT_NMA', &
            'the nma arm guard')
      end if
      armed = on
      fit_tol_in = fit_tol
   end subroutine nma_arm_from_input

   !------------------------------------------------------------------
   ! nma_armed() - the armed flag (the driver's stdout-naming condition)
   !------------------------------------------------------------------
   logical function nma_armed()
      nma_armed = armed
   end function nma_armed

   !------------------------------------------------------------------
   ! nma_derive() - the assembly-time derivation (driver step 9c): build
   !               the analysis-unit table (reactant fragments nat>=2 -
   !               under the surface paradigm the PROJECTILE fragment
   !               alone, substrate carriers excluded - UNION declared
   !               products nat>=2), guard it, pull the mode tables, set
   !               the label thresholds. Zero units (all monoatomic) is
   !               legal: no table, the note moves to setup
   !------------------------------------------------------------------
   subroutine nma_derive()
      type(pool_t), allocatable :: pools(:)
      real(8), allocatable :: wt(:), ct(:,:)
      real(8) :: qcm(3)
      integer :: i, n_unit, np, n_excl
      if (.not. armed) return
      if (derived) call stop_nma('nma_derive', 'derive already ran (one'// &
         ' analysis-unit table per process)')
      ! 1. collect the units: fragments nat>=2 first, then products nat>=2.
      !    Surface paradigm: the PROJECTILE fragment alone is a unit
      !    candidate - every other fragment is a substrate carrier,
      !    never a unit, never pulled (the sampler-side laziness mirror);
      !    no role (i_proj = 0) is an assembly mismatch
      if (reactants%surface_model /= 0 .and. list_atoms%i_proj == 0) then
         call stop_nma('nma_derive', 'surface paradigm (surface_model /= 0)'// &
            ' but list_atoms carries no PROJECTILE role (i_proj = 0) - the'// &
            ' surface paradigm legally has exactly one PROJECTILE fragment'// &
            ' (the incident-split ruling); the analysis units cannot be picked')
      end if
      n_unit = 0
      n_excl = 0
      do i = 1, size(list_atoms%frag)
         if (frag_is_unit(i)) then
            n_unit = n_unit + 1
         else if (reactants%surface_model /= 0 .and. i /= list_atoms%i_proj .and. &
                  list_atoms%frag(i)%nat >= 2) then
            n_excl = n_excl + 1            ! a substrate carrier that gas
         end if                            ! would have collected - excluded
      end do
      if (n_excl > 0) then
         call diag_kv('nma')
         write (*, '(a)') 'surface paradigm - substrate fragment(s)'// &
            ' excluded from analysis units'
      end if
      np = 0
      if (allocated(list_products)) np = size(list_products)
      do i = 1, np
         if (list_products(i)%nat >= 2) n_unit = n_unit + 1
      end do
      if (n_unit == 0) then
         derived = .true.
         return
      end if
      allocate (units(n_unit))
      n_unit = 0
      do i = 1, size(list_atoms%frag)
         if (.not. frag_is_unit(i)) cycle
         n_unit = n_unit + 1
         units(n_unit)%is_product = .false.
         units(n_unit)%frag = i
         units(n_unit)%name = list_atoms%frag(i)%name
         units(n_unit)%nat = list_atoms%frag(i)%nat
         units(n_unit)%list = list_atoms%frag(i)%list
         units(n_unit)%symb = list_atoms%symb(list_atoms%frag(i)%list)
         units(n_unit)%mass = list_atoms%mass(list_atoms%frag(i)%list)
         units(n_unit)%qz = list_atoms%frag(i)%qz
      end do
      do i = 1, np
         if (list_products(i)%nat < 2) cycle
         n_unit = n_unit + 1
         units(n_unit)%is_product = .true.
         units(n_unit)%iprod = i
         units(n_unit)%name = list_products(i)%name
         units(n_unit)%nat = list_products(i)%nat
         units(n_unit)%symb = list_products(i)%symb
         units(n_unit)%mass = list_products(i)%mass
         units(n_unit)%qz = list_products(i)%qz
      end do
      ! 2. COM-center every reference (the fit target frame)
      do i = 1, n_unit
         call center_ref(units(i)%qz, units(i)%mass, qcm)
      end do
      ! 3. product guards + subset expansion BEFORE any mode pull (the
      !    abort paths bind no container force)
      do i = 1, n_unit
         if (.not. units(i)%is_product) cycle
         call collide_check(i)
         call build_pools(i, pools)
         call multiset_check(i, pools)
         call cap_check(i, pools)
         call expand_subsets(i, pools)
      end do
      ! 4. reactant units: the single fixed candidate set
      do i = 1, n_unit
         if (units(i)%is_product) cycle
         units(i)%n_sub = 1
         allocate (units(i)%subs(units(i)%nat, 1))
         units(i)%subs(:, 1) = units(i)%list
      end do
      ! 5. mode pulls (fragments through the fragment seam, products
      !    through the product seam; nat>=2 only - no empty tables)
      do i = 1, n_unit
         if (units(i)%is_product) then
            call spectrum_prod_w(units(i)%iprod, wt, ct)
         else
            call spectrum_frag_w(units(i)%frag, wt, ct)
         end if
         if (size(wt) < 1) call stop_nma('nma_derive', 'unit '''//trim(units(i)%name)// &
            ''' pulled an empty mode table - a nat>=2 unit needs vibrational support')
         units(i)%w = wt
         units(i)%c = ct
      end do
      ! 6. label thresholds (NMA_FIT_TOL override or 0.4*min pair distance)
      do i = 1, n_unit
         if (fit_tol_in > 0.0d0) then
            units(i)%tol = fit_tol_in
         else
            units(i)%tol = 0.4d0*min_pair_dist(units(i)%qz)
         end if
      end do
      derived = .true.
   end subroutine nma_derive

   !------------------------------------------------------------------
   ! nma_setup() - the statistics-head phase: open the files (one
   !               nma_<name>.dat per unit + the shared label file). The
   !               derivation must already have run - the driver's
   !               assembly ordering is the contract - so opening has no
   !               failure path left. Zero units is legal: a stdout note,
   !               no files
   !------------------------------------------------------------------
   subroutine nma_setup()
      real(8), allocatable :: wt(:)
      integer :: i, k
      character(len=1024) :: s
      if (.not. armed) return
      if (.not. derived) call stop_nma('nma_setup', 'nma_derive must run first'// &
         ' (driver assembly ordering is the contract)')
      if (set_up) call stop_nma('nma_setup', 'setup already ran (one analysis-unit'// &
         ' table per process)')
      if (.not. allocated(units)) then
         call diag_kv('nma')
         write (*, '(a)') 'no analysis units (all monoatomic)'
         set_up = .true.
         return
      end if
      ! 1. files: one nma_<name>.dat per unit + the shared label file
      do i = 1, size(units)
         open (newunit=units(i)%u, file='nma_'//trim(units(i)%name)//'.dat', &
               status='replace', action='write')
         write (units(i)%u, '(a,a,a,i0,a,i0)') '# nma unit=', trim(units(i)%name), &
            ' nat=', units(i)%nat, ' modes=', size(units(i)%w)
         s = '# traj t drift atoms fit_rms e_kin_vib e_pot_vib e_har'
         do k = 1, size(units(i)%w)
            s = trim(s)//' n_'//i2s(k)
         end do
         write (units(i)%u, '(a)') trim(s)
         wt = [(w_to_wvn(units(i)%w(k)), k = 1, size(units(i)%w))]
         write (units(i)%u, '(a,999(1x,es24.16))') '# wavenumbers[cm-1]', wt
      end do
      open (newunit=u_label, file='nma_label.dat', status='replace', action='write')
      set_up = .true.
   end subroutine nma_setup

   !------------------------------------------------------------------
   ! nma_traj(i, sta, drift) - one trajectory's final state through every
   !               unit: best subset -> COM strip -> quat_fit -> rotated p
   !               -> mode projections -> one row per unit + one label row
   !------------------------------------------------------------------
   subroutine nma_traj(i, sta, drift)
      integer, intent(in) :: i
      type(state_t), intent(in) :: sta
      real(8), intent(in) :: drift
      real(8), allocatable :: qu(:,:), pu(:,:), qz3(:,:), qr(:,:), xw(:), pw(:), nv(:)
      real(8), allocatable :: rms_row(:)
      integer, allocatable :: sidx(:)
      real(8) :: u(3,3), rms, rms_best, eck, qk, pk, ekk, epk, ekin, epot, ehar
      integer :: iu, j, k, s, sbest
      character(len=512) :: atstr
      character(len=256) :: label
      logical :: all_formed
      if (.not. armed) return
      if (.not. derived .or. .not. set_up) call stop_nma('nma_traj', &
         'derive and setup must run first (the derive/open double flag)')
      if (.not. allocated(units)) return      ! zero-units run: no files, no labels
      ! per-row label inputs - allocated only past the allocated() guard
      ! (automatic bounds would evaluate size(units) at entry, UB on the
      ! zero-units degenerate path)
      allocate (rms_row(size(units)), sidx(size(units)))
      do iu = 1, size(units)
         associate (un => units(iu))
            qz3 = reshape(un%qz, [3, un%nat])
            ! 1. subset selection: products scan their element-compatible
            !    candidates for the smallest fit rms; reactant units carry
            !    the one fixed set
            rms_best = huge(1.0d0)
            sbest = 1
            do s = 1, un%n_sub
               call gather_unit(sta, un%subs(:, s), qu, pu)
               call strip_com(qu, pu, un%mass)
               call quat_fit(qz3, qu, un%mass, u, rms)
               if (rms < rms_best) then
                  rms_best = rms
                  sbest = s
               end if
            end do
            ! 2. the winning subset's fit frame (re-gather - the scan kept
            !    only the rms)
            call gather_unit(sta, un%subs(:, sbest), qu, pu)
            call strip_com(qu, pu, un%mass)
            call quat_fit(qz3, qu, un%mass, u, rms)
            qr = matmul(u, qu)
            do j = 1, un%nat
               call rot3(u, pu(:, j))
            end do
            eck = eckart_rot(qz3, qr, un%mass)
            ! 3. projections (samp_normalmode's realize folds, inverted):
            !    Q_k = sum sqrt(m) dq . c_k, P_k = sum p . c_k / sqrt(m)
            allocate (xw(3*un%nat), pw(3*un%nat), nv(size(un%w)))
            do j = 1, un%nat
               xw(3*j-2:3*j) = sqrt(un%mass(j))*(qr(:, j) - qz3(:, j))
               pw(3*j-2:3*j) = pu(:, j)/sqrt(un%mass(j))
            end do
            ekin = 0.0d0
            epot = 0.0d0
            do k = 1, size(un%w)
               qk = dot_product(xw, un%c(:, k))
               pk = dot_product(pw, un%c(:, k))
               ekk = 0.5d0*pk*pk/ev_dyn
               epk = 0.5d0*un%w(k)**2*qk*qk/ev_dyn
               nv(k) = (ekk + epk)/(hbar_evfs*un%w(k)) - 0.5d0
               ekin = ekin + ekk
               epot = epot + epk
            end do
            ehar = ekin + epot
            ! 4. the row (atoms column = space-joined global indices)
            atstr = ''
            do j = 1, un%nat
               atstr = trim(atstr)//' '//i2s(un%subs(j, sbest))
            end do
            write (un%u, row_fmt) i, sta%t, drift, adjustl(trim(atstr)), rms, ekin, epot, ehar, nv
            flush (un%u)
            ! 5. accumulators + this row's label inputs
            un%n_rows = un%n_rows + 1
            un%e_sum = un%e_sum + ehar
            un%rms_max = max(un%rms_max, rms)
            un%eck_max = max(un%eck_max, eck)
            rms_row(iu) = rms
            sidx(iu) = sbest
            deallocate (xw, pw, nv)
         end associate
      end do
      ! 6. the channel label: formed product stems joined '+' (formed =
      !    best rms < tol AND best subset != any fragment atom set - the
      !    elastic guard); no product label -> 'reactants' when every
      !    reactant fragment is formed (monoatomic fragments count as
      !    always formed - they are not units), else 'unassigned'
      label = ''
      do iu = 1, size(units)
         if (.not. units(iu)%is_product) cycle
         if (rms_row(iu) < units(iu)%tol .and. &
             .not. subset_is_frag(units(iu)%subs(:, sidx(iu)))) then
            if (len_trim(label) > 0) label = trim(label)//'+'
            label = trim(label)//trim(units(iu)%name)
         end if
      end do
      if (len_trim(label) == 0) then
         all_formed = .true.
         do iu = 1, size(units)
            if (.not. units(iu)%is_product) then
               if (rms_row(iu) >= units(iu)%tol) all_formed = .false.
            end if
         end do
         if (all_formed) then
            label = 'reactants'
         else
            label = 'unassigned'
         end if
      end if
      write (u_label, '(i0,1x,a)') i, trim(label)
      flush (u_label)
      call count_label(label)
   end subroutine nma_traj

   !------------------------------------------------------------------
   ! nma_summary() - per-unit stdout lines + the channel counts line,
   !                  then close every file
   !------------------------------------------------------------------
   subroutine nma_summary()
      integer :: iu
      real(8) :: eavg
      if (.not. armed .or. .not. set_up) return
      if (.not. allocated(units)) return
      do iu = 1, size(units)
         eavg = 0.0d0
         if (units(iu)%n_rows > 0) eavg = units(iu)%e_sum/units(iu)%n_rows
         call diag_kv('nma '//trim(units(iu)%name))
         write (*, '(a,i0,a,es12.4,a,es12.4,a,es12.4,a)') 'n=', units(iu)%n_rows, &
            ' E_vib=', eavg, ' eV max_fit_rms=', units(iu)%rms_max, &
            ' A max_eckart=', units(iu)%eck_max, ' amu*A^2'
         close (units(iu)%u)
         units(iu)%u = 0
      end do
      call diag_kv('channels')
      do iu = 1, n_outc
         write (*, '(a,i0)', advance='no') ' '//trim(outc_names(iu))//'=', outc_cnt(iu)
      end do
      call diag_eol()
      close (u_label)
      u_label = 0
   end subroutine nma_summary

   !==================================================================
   ! derive helpers
   !==================================================================

   !------------------------------------------------------------------
   ! frag_is_unit(i) - does fragment i become a reactant analysis unit?
   !               Gas paradigm: every nat>=2 fragment. Surface paradigm:
   !               the PROJECTILE fragment alone (nat>=2) - substrate
   !               carriers never become units, never pull a spectrum
   !               (the sampler-side laziness mirror)
   !------------------------------------------------------------------
   logical function frag_is_unit(i)
      integer, intent(in) :: i
      frag_is_unit = (list_atoms%frag(i)%nat >= 2)
      if (reactants%surface_model /= 0) then
         frag_is_unit = (i == list_atoms%i_proj .and. list_atoms%frag(i)%nat >= 2)
      end if
   end function frag_is_unit

   !------------------------------------------------------------------
   ! center_ref(qz, mass, qcm) - subtract the mass-weighted COM from a
   !               flat reference in place
   !------------------------------------------------------------------
   subroutine center_ref(qz, mass, qcm)
      real(8), intent(inout) :: qz(:)
      real(8), intent(in) :: mass(:)
      real(8), intent(out) :: qcm(3)
      integer :: j, nat
      nat = size(mass)
      qcm = 0.0d0
      do j = 1, nat
         qcm = qcm + mass(j)*qz(3*j-2:3*j)
      end do
      qcm = qcm/sum(mass)
      do j = 1, nat
         qz(3*j-2:3*j) = qz(3*j-2:3*j) - qcm
      end do
   end subroutine center_ref

   !------------------------------------------------------------------
   ! collide_check(iu) - product name vs the fragment-name namespace
   !               (case-insensitive) + the reserved label-file stem
   !------------------------------------------------------------------
   subroutine collide_check(iu)
      integer, intent(in) :: iu
      integer :: i
      if (ci_eq(units(iu)%name, 'label')) then
         call stop_nma('nma_derive', 'product name '''//trim(units(iu)%name)// &
            ''' collides with the reserved label-file stem (nma_label.dat)')
      end if
      do i = 1, size(list_atoms%frag)
         if (ci_eq(units(iu)%name, list_atoms%frag(i)%name)) then
            call stop_nma('nma_derive', 'product name '''//trim(units(iu)%name)// &
               ''' collides with reactant fragment '''//trim(list_atoms%frag(i)%name)// &
               ''' (case-insensitive) - the nma_<name>.dat namespace needs distinct names')
         end if
      end do
   end subroutine collide_check

   !------------------------------------------------------------------
   ! multiset_check(iu, pools) - every product element must be covered by
   !               the system's element multiset (absent element / count
   !               exceeded -> named abort naming the product)
   !------------------------------------------------------------------
   subroutine multiset_check(iu, pools)
      integer, intent(in) :: iu
      type(pool_t), intent(in) :: pools(:)
      character(len=512) :: msg
      integer :: j
      do j = 1, size(pools)
         if (size(pools(j)%idx) < pools(j)%k) then
            write (msg, '(a,a,a,a,a,i0,a,i0,a)') 'product ''', trim(units(iu)%name), &
               ''' declares element ', trim(pools(j)%sym), ' x', pools(j)%k, &
               ' but the system carries x', size(pools(j)%idx), &
               ' - products must be element rearrangements of the system atoms'
            call stop_nma('nma_derive', trim(msg))
         end if
      end do
   end subroutine multiset_check

   !------------------------------------------------------------------
   ! cap_check(iu, pools) - the element-pool product of binomials must
   !               stay under subset_cap (saturating count, named abort
   !               naming the product)
   !------------------------------------------------------------------
   subroutine cap_check(iu, pools)
      integer, intent(in) :: iu
      type(pool_t), intent(in) :: pools(:)
      character(len=512) :: msg
      integer :: j, ch, cnt
      cnt = 1
      do j = 1, size(pools)
         ch = choose_sat(size(pools(j)%idx), pools(j)%k)
         if (ch >= huge(0)) then
            cnt = huge(0)
            exit
         end if
         if (cnt > huge(0)/max(ch, 1)) then
            cnt = huge(0)
            exit
         end if
         cnt = cnt*ch
      end do
      if (cnt > subset_cap) then
         write (msg, '(a,a,a,i0,a,i0,a)') 'product ''', trim(units(iu)%name), &
            ''' needs ', cnt, ' element-compatible atom subsets - over the enumeration cap ', &
            subset_cap, ' (declare narrower products or a more specific system)'
         call stop_nma('nma_derive', trim(msg))
      end if
   end subroutine cap_check

   !------------------------------------------------------------------
   ! expand_subsets(iu, pools) - enumerate every element-compatible atom
   !               set (per-pool ascending combinations, sorted columns);
   !               cap_check has bounded the count
   !------------------------------------------------------------------
   subroutine expand_subsets(iu, pools)
      integer, intent(in) :: iu
      type(pool_t), intent(in) :: pools(:)
      integer :: cur(units(iu)%nat), nsub, cnt, j
      cnt = 1
      do j = 1, size(pools)
         cnt = cnt*choose_sat(size(pools(j)%idx), pools(j)%k)
      end do
      allocate (units(iu)%subs(units(iu)%nat, cnt))
      nsub = 0
      call walk_pools(pools, 1, cur, 0, units(iu)%subs, nsub)
      units(iu)%n_sub = nsub
      if (nsub /= cnt) call stop_nma('nma_derive', 'internal: subset expansion count'// &
         ' mismatch for product '''//trim(units(iu)%name)//'''')
   end subroutine expand_subsets

   !------------------------------------------------------------------
   ! build_pools(iu, pools) - the product's element pools: per distinct
   !               product element, the ascending system atom indices of
   !               that element and the declared draw count
   !------------------------------------------------------------------
   subroutine build_pools(iu, pools)
      integer, intent(in) :: iu
      type(pool_t), allocatable, intent(out) :: pools(:)
      logical :: dup
      integer :: j, k, m, npool, n
      n = size(list_atoms%symb)
      npool = 0
      do j = 1, units(iu)%nat
         dup = .false.
         do k = 1, j - 1
            if (ci4(units(iu)%symb(k)) == ci4(units(iu)%symb(j))) dup = .true.
         end do
         if (.not. dup) npool = npool + 1
      end do
      allocate (pools(npool))
      npool = 0
      do j = 1, units(iu)%nat
         dup = .false.
         do k = 1, j - 1
            if (ci4(units(iu)%symb(k)) == ci4(units(iu)%symb(j))) dup = .true.
         end do
         if (dup) cycle
         npool = npool + 1
         pools(npool)%sym = units(iu)%symb(j)
         m = 0
         do k = 1, units(iu)%nat
            if (ci4(units(iu)%symb(k)) == ci4(units(iu)%symb(j))) m = m + 1
         end do
         pools(npool)%k = m
         m = 0
         do k = 1, n
            if (ci4(list_atoms%symb(k)) == ci4(units(iu)%symb(j))) m = m + 1
         end do
         allocate (pools(npool)%idx(m))
         m = 0
         do k = 1, n
            if (ci4(list_atoms%symb(k)) == ci4(units(iu)%symb(j))) then
               m = m + 1
               pools(npool)%idx(m) = k
            end if
         end do
      end do
   end subroutine build_pools

   !------------------------------------------------------------------
   ! walk_pools / walk_choose - the two-level recursion: per pool choose
   !               k ascending indices, then descend into the next pool;
   !               a full assignment is sorted and stored as one column
   !------------------------------------------------------------------
   recursive subroutine walk_pools(pools, ip, cur, ncur, subs, nsub)
      type(pool_t), intent(in) :: pools(:)
      integer, intent(in) :: ip, ncur
      integer, intent(inout) :: cur(:), subs(:,:), nsub
      if (ip > size(pools)) then
         nsub = nsub + 1
         call sort_ins(cur, ncur)
         subs(1:ncur, nsub) = cur(1:ncur)
         return
      end if
      call walk_choose(pools(ip)%idx, pools(ip)%k, 1, 0, cur, ncur, pools, ip, subs, nsub)
   end subroutine walk_pools

   recursive subroutine walk_choose(idx, k, start, taken, cur, ncur, pools, ip, subs, nsub)
      integer, intent(in) :: idx(:), k, start, taken, ncur, ip
      integer, intent(inout) :: cur(:), subs(:,:), nsub
      type(pool_t), intent(in) :: pools(:)
      integer :: j
      if (taken == k) then
         call walk_pools(pools, ip + 1, cur, ncur + k, subs, nsub)
         return
      end if
      do j = start, size(idx) - (k - taken) + 1
         cur(ncur + taken + 1) = idx(j)
         call walk_choose(idx, k, j + 1, taken + 1, cur, ncur, pools, ip, subs, nsub)
      end do
   end subroutine walk_choose

   !==================================================================
   ! traj helpers
   !==================================================================

   !------------------------------------------------------------------
   ! gather_unit(sta, idx, qu, pu) - one atom set's Cartesian pair out of
   !               the raw final state
   !------------------------------------------------------------------
   subroutine gather_unit(sta, idx, qu, pu)
      type(state_t), intent(in) :: sta
      integer, intent(in) :: idx(:)
      real(8), allocatable, intent(out) :: qu(:,:), pu(:,:)
      integer :: j, g
      allocate (qu(3, size(idx)), pu(3, size(idx)))
      do j = 1, size(idx)
         g = idx(j)
         qu(:, j) = sta%q(3*g-2:3*g)
         pu(:, j) = sta%p(3*g-2:3*g)
      end do
   end subroutine gather_unit

   !------------------------------------------------------------------
   ! com_of(q3, mass, qcm) / pcom_of(p3, mass, vcm) - the mass-weighted
   !               COM and the COM velocity of one atom set (module
   !               private; cenmas eats state_t + list - unusable here)
   !------------------------------------------------------------------
   pure subroutine com_of(q3, mass, qcm)
      real(8), intent(in) :: q3(:,:), mass(:)
      real(8), intent(out) :: qcm(3)
      integer :: j
      qcm = 0.0d0
      do j = 1, size(mass)
         qcm = qcm + mass(j)*q3(:, j)
      end do
      qcm = qcm/sum(mass)
   end subroutine com_of

   pure subroutine pcom_of(p3, mass, vcm)
      real(8), intent(in) :: p3(:,:), mass(:)
      real(8), intent(out) :: vcm(3)
      real(8) :: ptot(3)
      integer :: j
      ptot = 0.0d0
      do j = 1, size(mass)
         ptot = ptot + p3(:, j)
      end do
      vcm = ptot/sum(mass)
   end subroutine pcom_of

   !------------------------------------------------------------------
   ! strip_com(qu, pu, mass) - translate to the unit COM frame and remove
   !               the unit's net momentum
   !------------------------------------------------------------------
   subroutine strip_com(qu, pu, mass)
      real(8), intent(inout) :: qu(:,:), pu(:,:)
      real(8), intent(in) :: mass(:)
      real(8) :: qcm(3), vcm(3)
      integer :: j
      call com_of(qu, mass, qcm)
      call pcom_of(pu, mass, vcm)
      do j = 1, size(mass)
         qu(:, j) = qu(:, j) - qcm
         pu(:, j) = pu(:, j) - mass(j)*vcm
      end do
   end subroutine strip_com

   !------------------------------------------------------------------
   ! subset_is_frag(sub) - does an atom set equal any reactant fragment's
   !               atom set (the elastic guard; both sides ascending)
   !------------------------------------------------------------------
   logical function subset_is_frag(sub)
      integer, intent(in) :: sub(:)
      integer :: i
      subset_is_frag = .false.
      do i = 1, size(list_atoms%frag)
         associate (fl => list_atoms%frag(i)%list)
            if (size(fl) == size(sub)) then
               if (all(fl == sub)) then
                  subset_is_frag = .true.
                  return
               end if
            end if
         end associate
      end do
   end function subset_is_frag

   !------------------------------------------------------------------
   ! count_label(label) - the channel counting table (insertion order)
   !------------------------------------------------------------------
   subroutine count_label(label)
      character(len=*), intent(in) :: label
      integer :: j
      do j = 1, n_outc
         if (trim(outc_names(j)) == trim(label)) then
            outc_cnt(j) = outc_cnt(j) + 1
            return
         end if
      end do
      if (n_outc >= label_cap) then
         call stop_nma('nma_traj', 'distinct channel labels exceed the counting-table'// &
            ' capacity '//i2s(label_cap))
      end if
      n_outc = n_outc + 1
      outc_names(n_outc) = label
      outc_cnt(n_outc) = 1
   end subroutine count_label

   !==================================================================
   ! small shared helpers
   !==================================================================

   ! min_pair_dist(qz) - the reference's smallest interatomic distance [A]
   pure function min_pair_dist(qz) result(dmin)
      real(8), intent(in) :: qz(:)
      real(8) :: dmin, d
      integer :: a, b, nat
      nat = size(qz)/3
      dmin = huge(1.0d0)
      do a = 1, nat - 1
         do b = a + 1, nat
            d = norm2(qz(3*a-2:3*a) - qz(3*b-2:3*b))
            dmin = min(dmin, d)
         end do
      end do
   end function min_pair_dist

   ! choose_sat(n, k) - C(n,k) exact, saturated at huge(0) (the
   ! multiplicative form stays integral at every step)
   pure function choose_sat(n, k) result(r)
      integer, intent(in) :: n, k
      integer :: j, r
      if (k < 0 .or. k > n) then
         r = 0
         return
      end if
      r = 1
      do j = 1, k
         if (r > huge(0)/(n - k + j)) then
            r = huge(0)
            return
         end if
         r = r*(n - k + j)/j
      end do
   end function choose_sat

   ! sort_ins(a, n) - in-place ascending insertion sort
   pure subroutine sort_ins(a, n)
      integer, intent(inout) :: a(:)
      integer, intent(in) :: n
      integer :: i, j, t
      do i = 2, n
         t = a(i)
         j = i - 1
         do while (j >= 1 .and. a(j) > t)
            a(j + 1) = a(j)
            j = j - 1
         end do
         a(j + 1) = t
      end do
   end subroutine sort_ins

   ! i2s(n) - decimal string of an integer
   pure function i2s(n) result(s)
      integer, intent(in) :: n
      character(len=16) :: s
      write (s, '(i0)') n
   end function i2s

   ! ci_eq(a, b) - case-insensitive equality of two trimmed names
   pure logical function ci_eq(a, b)
      character(len=*), intent(in) :: a, b
      integer :: i, ca, cb
      ci_eq = len_trim(a) == len_trim(b)
      if (.not. ci_eq) return
      do i = 1, len_trim(a)
         ca = iachar(a(i:i))
         cb = iachar(b(i:i))
         if (ca >= iachar('a') .and. ca <= iachar('z')) ca = ca - 32
         if (cb >= iachar('a') .and. cb <= iachar('z')) cb = cb - 32
         if (ca /= cb) then
            ci_eq = .false.
            return
         end if
      end do
   end function ci_eq

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

   ! stop_nma(who, msg) - the named-abort channel
   subroutine stop_nma(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': '//trim(msg)
      call diag_fatal(diag_cat_output, 'nma.f90', trim(who), &
         'the nma analyzer is unusable')
   end subroutine stop_nma

end module nma
