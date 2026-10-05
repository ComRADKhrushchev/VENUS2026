!=====================================================================
! recorder.f90 - column-style recorder of trajectory frames and final-state rows
! Design:
!   A column is a registry row (name + extractor); selection is preset level (the
!   base) plus an appended name list (union). Levels: 0/1 = frame header + energy
!   line; 2 = + all-atom Q/P section (the bare-input default); 5 = + per-atom
!   force section. Energy columns report kcal/mol; electronic columns read the
!   elec_interface statistics getters only.
!=====================================================================
module recorder
   use state,         only: state_t
   use config,        only: observables
   use config_atoms, only: list_atoms
   use force_interface,  only: force_eval
   use elec_interface,   only: occ_get, n_hop_get
   use consts,        only: e_conv
   implicit none
   private
   public :: col_val_i, rec_reg_col, rec_select, rec_frame, rec_final, rec_n_sel
   public :: rec_fin_val                                            ! final-row value getter
                                                                    ! (one named column of the
                                                                    ! last written final row)

   ! Column-extractor interface: one scalar observable from a physical-state snapshot
   abstract interface
      real(8) function col_val_i(sta)
         import :: state_t
         type(state_t), intent(in) :: sta   ! physical state (read-only snapshot)
      end function col_val_i
   end interface

   ! Column registry (column = name + extractor; filled at assembly, read-only while
   ! recording)
   integer, parameter :: max_col = 64        ! registry capacity [column]
   type :: col_entry_t
      character(len=16) :: name = ''         ! column name (obs_list list key; 1-16 characters)
      integer :: lvl_min = huge(1)           ! minimum preset level including
                                             ! this column (huge = no preset group - obs_list-only)
      procedure(col_val_i), pointer, nopass :: val => null() ! extractor
   end type col_entry_t
   type(col_entry_t) :: col_tbl(max_col)     ! column registry (module-private)
   integer :: n_col = 0                      ! registered column count [column]
   integer :: sel_idx(max_col) = 0           ! registry indices of the selected columns [-]
   integer :: n_sel = 0                      ! selected column count [column] (rec_select assembly)
   logical :: sec_qp  = .false.              ! all-atom Q/P flattened section (sec_ = section;
                                                   !    gate at rec_level>=2)
   logical :: sec_frc = .false.              ! per-atom force section (gate at rec_level>=5)
   logical :: sel_done = .false.             ! rec_select ran (assembly-order guard for the
                                             ! write entries)
   logical :: fw_done = .false.              ! framework preregistration flag (lazy, first use)
   integer :: traj_i = 1                     ! trajectory counter [-] (frame header / final row)
   integer :: step_i = 0                     ! per-trajectory frame counter [-]
   real(8) :: fin_val(max_col) = 0.0d0       ! last written final row's column values
                                             ! (the rec_fin_val read-back stash)
   integer :: fin_n = 0                      ! the stashed row's column count [column]
   logical :: fin_done = .false.             ! a final row has been written (the getter's
                                             ! assembly-order guard)

   ! Section-gate levels: rec_sel_lvl derives both gates from these constants
   integer, parameter :: lvl_qp  = 2         ! sec_qp opens at this preset level
   integer, parameter :: lvl_frc = 5         ! sec_frc opens at this preset level
   ! the closed preset-level domain
   integer, parameter :: lvl_domain(4) = [ 0, 1, 2, 5 ]

   ! Column-selection generic: preset level (rec_level key) or name list (obs_list key)
   interface rec_select
      module procedure rec_sel_lvl   ! by preset level (the rec_level key)
      module procedure rec_sel_lst   ! by list (the obs_list key - comma-separated column names)
   end interface rec_select
contains
   !------------------------------------------------------------------
   ! rec_reg_col(name, proc) - register an external observable column (enters with
   ! lvl_min = huge, i.e. obs_list-only; duplicate name / full table / overlong name
   ! aborts with a named error)
   !------------------------------------------------------------------
   subroutine rec_reg_col(name, proc)
      character(len=*), intent(in) :: name  ! column name (obs_list list key)
      procedure(col_val_i) :: proc          ! column extractor
      call fw_reg()
      call reg_row(name, proc, huge(1))
   end subroutine rec_reg_col

   !------------------------------------------------------------------
   ! rec_select(level / list) - column-selection generic: the preset-level entry SETS the
   ! selection and both section gates (base); the list entry APPENDS named columns and may
   ! open 'all_qp'/'all_frc' - the recorded set is the union of both entries
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   ! rec_sel_lvl(level) - select columns by preset level (single scan over lvl_min;
   ! replaces any previous selection, sets both section gates)
   !------------------------------------------------------------------
   subroutine rec_sel_lvl(level)
      integer, intent(in) :: level   ! preset level [-] (0/1/2/5; the default call passes
                                     ! observables%rec_level)
      integer :: j
      ! 1. lazy framework preregistration + closed level-domain validation
      call fw_reg()
      if (all(level /= lvl_domain)) then
         call stop_rec('rec_sel_lvl', 'unknown preset level - the recording-level set is '// &
                       '0/1/2/5 (GWRITE_LEVEL semantics, N3.1)')
      end if
      ! 2. single scan over lvl_min collects the preset group (replaces any selection)
      n_sel = 0
      do j = 1, n_col
         if (col_tbl(j)%lvl_min <= level) then
            n_sel = n_sel + 1
            sel_idx(n_sel) = j
         end if
      end do
      ! 3. set both section gates by the level
      sec_qp  = (level >= lvl_qp)
      sec_frc = (level >= lvl_frc)
      sel_done = .true.
   end subroutine rec_sel_lvl

   !------------------------------------------------------------------
   ! rec_sel_lst(list) - append comma-separated column names to the current selection
   ! (case-insensitive lookup; the reserved names 'all_qp'/'all_frc' open the sections)
   !------------------------------------------------------------------
   subroutine rec_sel_lst(list)
      character(len=*), intent(in) :: list  ! column-name list (comma-separated;
                                            ! the reserved section names 'all_qp'/'all_frc'
                                            ! open the corresponding flattened section)
      character(len=len(list)) :: rest
      character(len=16) :: tok
      integer :: ic, j
      call fw_reg()
      if (len_trim(list) == 0) then
         sel_done = .true.                   ! empty list key = no columns to add (legal no-op)
         return
      end if
      ! 1. one loop pass per comma-separated token
      rest = adjustl(list)
      do while (len_trim(rest) > 0)
         ic = index(trim(rest), ',')
         if (ic == 0) then
            tok = adjustl(trim(rest))
            rest = ''
         else
            tok = adjustl(rest(:ic-1))
            rest = adjustl(rest(ic+1:))
         end if
         if (len_trim(tok) == 0) then
            call stop_rec('rec_sel_lst', 'empty item in the column list (obs_list grammar: '// &
                          'comma-separated names, no empty items)')
         end if
         ! 2. dispatch: reserved section gates vs registry lookup, append if absent
         if (to_upper(trim(tok)) == 'ALL_QP') then
            sec_qp = .true.
         else if (to_upper(trim(tok)) == 'ALL_FRC') then
            sec_frc = .true.
         else
            j = col_lookup(tok)
            if (j == 0) then
               call stop_rec('rec_sel_lst', 'unknown column name "'//trim(tok)// &
                             '" (not in the recorder column registry)')
            end if
            if (.not. any(sel_idx(1:n_sel) == j)) then
               if (n_sel >= max_col) then
                  call stop_rec('rec_sel_lst', 'selection table full ('// &
                                'too many selected columns)')
               end if
               n_sel = n_sel + 1
               sel_idx(n_sel) = j
            end if
         end if
      end do
      sel_done = .true.
   end subroutine rec_sel_lst

   !------------------------------------------------------------------
   ! rec_frame(sta, unit) - write one checkpoint frame: selected column values plus the
   ! gated flattened sections
   ! Frame layout: '# frame traj= step= t=' header; then one 'name=value' line (selection
   ! order); '# sec_qp begin'/per-atom q,p rows/'# sec_qp end' gated at level>=2;
   ! '# sec_frc begin'/per-atom force rows/'# sec_frc end' gated at level>=5. Atom row =
   ! '<index> <symb> <3 reals> [<3 reals>]'; t = sta%t [10 fs] as held, energies
   ! [kcal/mol]; flush(unit) only at the final row
   !------------------------------------------------------------------
   subroutine rec_frame(sta, unit)
      type(state_t), intent(in) :: sta   ! physical state (read-only post-statistics snapshot)
      integer, intent(in) :: unit        ! output unit number [-] (driver owns open/close)
      integer :: j, k, nat
      ! 1. guards + frame counters
      call fw_reg()
      call check_ready(sta, 'rec_frame')
      nat = size(list_atoms%mass)
      step_i = step_i + 1
      ! 2. frame header + one name=value line of selected column values
      write (unit, '(a,i0,a,i0,a,es24.16)') '# frame traj=', traj_i, ' step=', step_i, ' t=', sta%t
      if (n_sel > 0) then
         do j = 1, n_sel
            if (j == 1) then
               write (unit, '(a)', advance='no') trim(trim(col_tbl(sel_idx(j))%name)//'='// &
                  num_str(col_tbl(sel_idx(j))%val(sta)))
            else
               write (unit, '(a)', advance='no') trim(' '//trim(col_tbl(sel_idx(j))%name)//'='// &
                  num_str(col_tbl(sel_idx(j))%val(sta)))
            end if
         end do
         write (unit, '(a)') ''
      end if
      ! 3. all-atom Q/P flattened section (gate: preset level >= 2)
      if (sec_qp) then
         write (unit, '(a)') '# sec_qp begin'
         do k = 1, nat
            write (unit, '(i0,1x,a,1x,6es24.16)') k, list_atoms%symb(k), &
               sta%q(3*k-2:3*k), sta%p(3*k-2:3*k)
         end do
         write (unit, '(a)') '# sec_qp end'
      end if
      ! 4. per-atom force section (gate: preset level >= 5)
      if (sec_frc) then
         write (unit, '(a)') '# sec_frc begin'
         do k = 1, nat
            write (unit, '(i0,1x,a,1x,3es24.16)') k, list_atoms%symb(k), sta%f(3*k-2:3*k)
         end do
         write (unit, '(a)') '# sec_frc end'
      end if
      ! 5. no per-frame flush: on high-latency object stores a per-frame
      !    fsync dominates the wall clock (measured ~100x trajectory
      !    throughput loss); durability of mid-trajectory frames is not a
      !    production requirement - the final-row flush below carries the
      !    record. Byte-neutral for any completed archive.
   end subroutine rec_frame

   !------------------------------------------------------------------
   ! rec_final(sta, unit) - write the trajectory final-state row: '# final: traj t
   ! <column names>' header once (first trajectory), then one row per trajectory =
   ! traj, t, selected column values. Values are computed ONCE into the fin_val stash
   ! shared with rec_fin_val; afterwards traj_i advances and step_i resets
   !------------------------------------------------------------------
   subroutine rec_final(sta, unit)
      type(state_t), intent(in) :: sta   ! physical state (final-state snapshot - read-only)
      integer, intent(in) :: unit        ! output unit number [-]
      integer :: j
      ! 1. guards
      call fw_reg()
      call check_ready(sta, 'rec_final')
      ! 2. evaluate the selected columns once into the stash
      do j = 1, n_sel
         fin_val(j) = col_tbl(sel_idx(j))%val(sta)
      end do
      ! 3. one-time name header, then the data row
      if (traj_i == 1) then
         write (unit, '(a)', advance='no') '# final: traj t'
         do j = 1, n_sel
            write (unit, '(a)', advance='no') ' '//trim(col_tbl(sel_idx(j))%name)
         end do
         write (unit, '(a)') ''
      end if
      write (unit, '(i0,1x,es24.16)', advance='no') traj_i, sta%t
      do j = 1, n_sel
         write (unit, '(1x,es24.16)', advance='no') fin_val(j)
      end do
      write (unit, '(a)') ''
      flush (unit)
      ! 4. publish the stash and advance the frame/final counters
      fin_n = n_sel
      fin_done = .true.
      traj_i = traj_i + 1
      step_i = 0
   end subroutine rec_final

   !------------------------------------------------------------------
   ! rec_fin_val(name) - read one named column's value back out of the LAST written
   ! final-state row (call between a rec_final and the next; the name must be among
   ! the selected columns)
   !------------------------------------------------------------------
   function rec_fin_val(name) result(v)
      character(len=*), intent(in) :: name  ! column name (must be among the SELECTED
                                            ! columns - the row carries those only)
      real(8) :: v                          ! the stashed value [the column's unit]
      integer :: j, row
      if (.not. fin_done) then
         call stop_rec('rec_fin_val', 'no final-state row written yet (rec_final must'// &
                       ' run first - a statistics read before any recording)')
      end if
      row = 0
      do j = 1, n_sel
         if (to_upper(trim(col_tbl(sel_idx(j))%name)) == to_upper(trim(name))) row = j
      end do
      if (row == 0) then
         call stop_rec('rec_fin_val', 'column "'//trim(name)//'" is not among the selected'// &
                       ' columns (the final row carries the selection only)')
      end if
      v = fin_val(row)
   end function rec_fin_val

   !------------------------------------------------------------------
   ! private helpers
   !------------------------------------------------------------------

   ! fw_reg() - lazy framework preregistration (runs once, on first registry use): the
   ! energy three at preset level 0 + the occ/n_hop placeholders (lvl_min = huge)
   subroutine fw_reg()
      if (fw_done) return
      fw_done = .true.
      call reg_row('t_kin', col_tkin, 0)
      call reg_row('v_pot', col_vpot, 0)
      call reg_row('e_tot', col_etot, 0)
      call reg_row('occ',   col_occ,  huge(1))
      call reg_row('n_hop', col_nhop, huge(1))
   end subroutine fw_reg

   ! reg_row(name, proc, lvl) - append one registry row (shared by the public entry and the
   ! framework preregistration); duplicate name / full table / overlong name aborts with a named error
   subroutine reg_row(name, proc, lvl)
      character(len=*), intent(in) :: name
      procedure(col_val_i) :: proc
      integer, intent(in) :: lvl
      if (len_trim(name) == 0 .or. len_trim(name) > 16) then
         call stop_rec('rec_reg_col', 'illegal column name length ("'//trim(name)// &
                       '"; 1-16 characters)')
      end if
      if (col_lookup(name) /= 0) then
         call stop_rec('rec_reg_col', 'duplicate column name "'//trim(name)// &
                       '" (already in the recorder column registry)')
      end if
      if (n_col >= max_col) then
         call stop_rec('rec_reg_col', 'column registry full ('// &
                       'too many registered columns; max_col exceeded)')
      end if
      n_col = n_col + 1
      col_tbl(n_col)%name = name
      col_tbl(n_col)%lvl_min = lvl
      col_tbl(n_col)%val => proc
   end subroutine reg_row

   ! col_lookup(name) - case-insensitive name lookup; 0 = not found
   integer function col_lookup(name)
      character(len=*), intent(in) :: name
      integer :: j
      col_lookup = 0
      do j = 1, n_col
         if (to_upper(trim(col_tbl(j)%name)) == to_upper(trim(name))) then
            col_lookup = j
            return
         end if
      end do
   end function col_lookup

   ! check_ready(sta, who) - write-path guards: selection made, list_atoms assembled, and the
   ! state an EXACT size fit to the atom list in both directions (an oversized state would
   ! read list_atoms%mass out of bounds in col_tkin)
   subroutine check_ready(sta, who)
      type(state_t), intent(in) :: sta
      character(len=*), intent(in) :: who
      integer :: nat
      if (.not. sel_done) then
         call stop_rec(who, 'called before any rec_select (assembly-order error: '// &
                       'column selection precedes recording)')
      end if
      if (.not. allocated(list_atoms%mass) .or. .not. allocated(list_atoms%symb)) then
         call stop_rec(who, 'list_atoms not assembled (list_atoms_load must run before recording)')
      end if
      nat = size(list_atoms%mass)
      if (size(sta%q) /= 3*nat .or. size(sta%p) /= 3*nat .or. size(sta%f) /= 3*nat) then
         call stop_rec(who, 'state/list_atoms size mismatch (state carries '//trim(i8_str(size(sta%q)))// &
                       ' q-dofs, the atom list carries '//trim(i8_str(nat))//' atoms - exact fit required)')
      end if
   end subroutine check_ready

   ! i8_str(i) - integer rendered as a string (the size-mismatch named abort names both sides)
   pure function i8_str(i) result(s)
      integer, intent(in) :: i
      character(len=12) :: s
      write (s, '(i0)') i
      s = trim(adjustl(s))
   end function i8_str

   ! stop_rec(who, msg) - the named-abort channel (same semantics as config_atoms /
   ! sysdef_interface: name the caller, print the named message, STOP 1)
   subroutine stop_rec(who, msg)
      character(len=*), intent(in) :: who, msg
      write (0, '(a)') trim(who)//': recording error: '//trim(msg)
      write (0, '(a)') trim(who)//': fatal (the recording is unusable)'
      stop 1
   end subroutine stop_rec

   ! num_str(v) - a real rendered as a blank-free ES string (keeps the column line's
   ! name=value tokens single-space-separated)
   pure function num_str(v) result(s)
      real(8), intent(in) :: v
      character(len=32) :: s
      write (s, '(es24.16)') v
      s = trim(adjustl(s))
   end function num_str

   ! to_upper(s) - ASCII uppercase of a string (the column-name lookup is
   ! case-insensitive)
   pure function to_upper(s)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: to_upper
      integer :: i, ic
      do i = 1, len(s)
         ic = iachar(s(i:i))
         if (ic >= iachar('a') .and. ic <= iachar('z')) then
            to_upper(i:i) = achar(ic - 32)
         else
            to_upper(i:i) = s(i:i)
         end if
      end do
   end function to_upper

   !------------------------------------------------------------------
   ! framework column extractors (the preregistered set)
   !------------------------------------------------------------------

   ! t_kin - kinetic energy T = sum_atoms p.p/(2 m) / e_conv [kcal/mol] (list_atoms masses)
   real(8) function col_tkin(sta) result(v)
      type(state_t), intent(in) :: sta
      integer :: k
      v = 0.0d0
      do k = 1, size(sta%p)/3
         v = v + dot_product(sta%p(3*k-2:3*k), sta%p(3*k-2:3*k))/(2.0d0*list_atoms%mass(k))
      end do
      v = v/e_conv
   end function col_tkin

   ! v_pot - potential energy via force_eval [kcal/mol]: V = H - T on a private copy of
   ! the snapshot (force_eval takes an inout state; the snapshot stays read-only)
   real(8) function col_vpot(sta) result(v)
      type(state_t), intent(in) :: sta
      type(state_t) :: wrk
      real(8) :: e
      wrk = sta
      e = 0.0d0
      call force_eval(wrk, e)
      v = e - col_tkin(sta)
   end function col_vpot

   ! e_tot - total energy H = T(P) + V via force_eval [kcal/mol] (same private-copy
   ! pattern as col_vpot)
   real(8) function col_etot(sta) result(v)
      type(state_t), intent(in) :: sta
      type(state_t) :: wrk
      wrk = sta
      v = 0.0d0
      call force_eval(wrk, v)
   end function col_etot

   ! occ - occupation column: the ACTIVE-STATE INDEX of the held occupation set, read via
   ! occ_get (seeded [1] at the method bind, so the adiabatic row records 1.0 on every
   ! line; no whole-vector scalar summary exists yet)
   real(8) function col_occ(sta) result(v)
      type(state_t), intent(in) :: sta
      integer, allocatable :: occ(:)
      call occ_get(occ)
      v = dble(occ(1))
      deallocate (occ)
   end function col_occ

   ! n_hop - hop-count column: the held hop counter via n_hop_get (a pure observable -
   ! the recorder only reads; 0 under the adiabatic row, incremented at hop detection
   ! inside the hopping member)
   real(8) function col_nhop(sta) result(v)
      type(state_t), intent(in) :: sta
      v = dble(n_hop_get())
   end function col_nhop

   !------------------------------------------------------------------
   ! rec_n_sel() - selected column count (read-only public probe)
   !------------------------------------------------------------------
   pure integer function rec_n_sel()
      rec_n_sel = n_sel
   end function rec_n_sel
end module recorder
