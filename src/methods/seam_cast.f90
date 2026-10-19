!=====================================================================
! seam_cast.f90 - buffered-string -> number cast for the assembly seam
! Design:
!   Pure form: key name, raw string and source line in, converted
!   number out - no member word, no driver symbol, no buffer accessor
!   (the caller hands what it already hoisted).
!   An unreadable value (e.g. a hand typo like T_GLO = 3O0) is a named
!   named abort here, never a bare runtime read error dying mid-assembly;
!   the message shape (key + line + the unreadable value, quoted)
!   matches the input framework's own value parsers.
!=====================================================================
module seam_cast
   use diag, only: diag_fatal, diag_cat_system
   implicit none
   private
   public :: seam_real, seam_int
contains
   !------------------------------------------------------------------
   ! seam_real(key, sval, ln, rval) - cast one buffered value to real(8);
   !              an unreadable string aborts with a named error (diag_fatal)
   subroutine seam_real(key, sval, ln, rval)
      character(len=*), intent(in) :: key, sval
      integer, intent(in) :: ln
      real(8), intent(out) :: rval
      integer :: ios
      read (sval, *, iostat=ios) rval
      if (ios == 0) return
      write (0, '(a)') 'seam_cast: key '//trim(key)//trim(seam_tag(ln))// &
                       ': bad real value "'//trim(sval)//'"'
      call diag_fatal(diag_cat_system, 'seam_cast.f90', 'seam_cast', &
         'the member parameter is not loaded')
   end subroutine seam_real

   !------------------------------------------------------------------
   ! seam_int(key, sval, ln, ival) - cast one buffered value to integer;
   !              an unreadable string aborts with a named error (diag_fatal)
   subroutine seam_int(key, sval, ln, ival)
      character(len=*), intent(in) :: key, sval
      integer, intent(in) :: ln
      integer, intent(out) :: ival
      integer :: ios
      read (sval, *, iostat=ios) ival
      if (ios == 0) return
      write (0, '(a)') 'seam_cast: key '//trim(key)//trim(seam_tag(ln))// &
                       ': bad integer value "'//trim(sval)//'"'
      call diag_fatal(diag_cat_system, 'seam_cast.f90', 'seam_cast', &
         'the member parameter is not loaded')
   end subroutine seam_int

   ! seam_tag(ln) - the source-line suffix of a message (a buffered key arrives
   !                from a data line, ln >= 1; ln <= 0 falls back to the
   !                default tag)
   function seam_tag(ln) result(tag)
      integer, intent(in) :: ln
      character(len=24) :: tag
      if (ln > 0) then
         write (tag, '(a,i0,a)') ' (line ', ln, ')'
      else
         tag = ' (default)'
      end if
   end function seam_tag
end module seam_cast
