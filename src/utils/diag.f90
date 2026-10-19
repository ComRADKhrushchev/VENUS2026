module diag
!------------------------------------------------------------------
! Console diagnostics single entry: the stderr fatal frame (fixed form -
! category tag + file + site + the byte-preserved legacy line) and the
! stdout layout helpers (heavy banner / light rule / label column).
! Design note: pure output module - touches no RNG and no state; every
! fatal site collapses here so the frame format has one truth point
! (2026-10-15 ruling, docs/plans/2026-10-15-stdout-stderr-format.md).
!------------------------------------------------------------------
   implicit none
   private
   public :: diag_fatal, diag_banner, diag_heavy, diag_rule, diag_kv, diag_kv_cont, diag_eol
   public :: diag_cat_input, diag_cat_system, diag_cat_assembly, diag_cat_sampler, &
             diag_cat_dynamics, diag_cat_output, diag_cat_internal
   integer, parameter :: diag_cat_input    = 1   ! input parse/audit
   integer, parameter :: diag_cat_system   = 2   ! container/PES/seam/list_atoms
   integer, parameter :: diag_cat_assembly = 3   ! registry/schedule wiring
   integer, parameter :: diag_cat_sampler  = 4   ! distribution/bath/incident members
   integer, parameter :: diag_cat_dynamics = 5   ! propagation/force evaluation
   integer, parameter :: diag_cat_output   = 6   ! statistics/products output
   integer, parameter :: diag_cat_internal = 7   ! defensive invariants
   integer, parameter :: frame_w = 68            ! stderr frame line width [char]
   integer, parameter :: rule_w  = 69            ! stdout rule width [char]
   integer, parameter :: label_w = 12            ! stdout label-column width [char]
contains
   !------------------------------------------------------------------
   ! diag_fatal(cat, file, site, msg) - the stderr frame, then stop 1; the
   !   third frame line reconstructs the legacy 'site: fatal (msg)' byte for
   !   byte so every substring assertion in the check family survives.
   !------------------------------------------------------------------
   subroutine diag_fatal(cat, file, site, msg)
      integer, intent(in) :: cat
      character(len=*), intent(in) :: file, site, msg
      write (0, '(a)') repeat('!', frame_w)
      write (0, '(a,a,a,a,a,a,a)') '!!!  FATAL  [', trim(cat_name(cat)), ']  ', &
         trim(file), ' :: ', trim(site)
      write (0, '(a,a,a,a,a)') '!!!  ', trim(site), ': fatal (', trim(msg), ')'
      write (0, '(a)') repeat('!', frame_w)
      stop 1
   end subroutine diag_fatal

   ! diag_banner(title) - heavy rule + title + heavy rule (block open)
   subroutine diag_banner(title)
      character(len=*), intent(in) :: title
      call diag_heavy()
      write (*, '(a)') ' '//trim(title)
      call diag_heavy()
   end subroutine diag_banner

   ! diag_heavy() / diag_rule() - the single heavy / light rule lines
   subroutine diag_heavy()
      write (*, '(a)') repeat('=', rule_w)
   end subroutine diag_heavy

   subroutine diag_rule()
      write (*, '(a)') repeat('-', rule_w)
   end subroutine diag_rule

   ! diag_kv(label) - the label column, no advance (caller writes the values)
   subroutine diag_kv(label)
      character(len=*), intent(in) :: label
      character(len=256) :: lab
      lab = adjustl(label)
      if (len_trim(lab) < label_w) then
         write (*, '(a,a,a)', advance='no') ' ', lab(:label_w), ' '
      else
         write (*, '(a,a,a)', advance='no') ' ', trim(lab), ' '
      end if
   end subroutine diag_kv

   ! diag_kv_cont() - continuation indent to the value column, no advance
   subroutine diag_kv_cont()
      write (*, '(a)', advance='no') repeat(' ', label_w + 2)
   end subroutine diag_kv_cont

   ! diag_eol() - end a no-advance line
   subroutine diag_eol()
      write (*, '(a)') ''
   end subroutine diag_eol

   ! cat_name(cat) - the category tag word
   pure function cat_name(cat) result(s)
      integer, intent(in) :: cat
      character(len=9) :: s
      select case (cat)
      case (diag_cat_input);    s = 'INPUT'
      case (diag_cat_system);   s = 'SYSTEM'
      case (diag_cat_assembly); s = 'ASSEMBLY'
      case (diag_cat_sampler);  s = 'SAMPLER'
      case (diag_cat_dynamics); s = 'DYNAMICS'
      case (diag_cat_output);   s = 'OUTPUT'
      case default;             s = 'INTERNAL'
      end select
   end function cat_name
end module diag
