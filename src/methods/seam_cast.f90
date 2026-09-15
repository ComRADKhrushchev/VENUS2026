!=====================================================================
! seam_cast.f90 - staged-string -> number cast for the assembly seam
!                 (methods/): the one shared conversion the reg
!                 member-parameter loaders read their staged raw values
!                 through - an unreadable value (a hand typo like
!                 T_GLO = 3O0) is a NAMED loud-stop here, never a bare
!                 runtime read error dying mid-assembly
! Placement (the zero-knowledge audit): the cast is PURE FORM - key
!   name, raw string and source line in, converted number out; no
!   member word, no driver symbol, no staging accessor (the caller
!   hands what it already hoisted). It lives beside the seam it serves,
!   so the input framework keeps its posture (value semantics are
!   member-side knowledge) and every reg file stays a thin loader.
! Message shape: the input framework's own value parsers - key + line
!   + the unreadable value, quoted - so a bad value reads identically
!   wherever it is caught.
!=====================================================================
module seam_cast
   implicit none
   private
   public :: seam_real, seam_int
contains
   !------------------------------------------------------------------
   ! seam_real(key, sval, ln, rval) - cast one staged value to real(8);
   !              an unreadable string loud-stops named (STOP 1)
   subroutine seam_real(key, sval, ln, rval)
      character(len=*), intent(in) :: key, sval
      integer, intent(in) :: ln
      real(8), intent(out) :: rval
      integer :: ios
      read (sval, *, iostat=ios) rval
      if (ios == 0) return
      write (0, '(a)') 'seam_cast: key '//trim(key)//trim(seam_tag(ln))// &
                       ': bad real value "'//trim(sval)//'"'
      write (0, '(a)') 'seam_cast: loud-stop (the member parameter is not loaded)'
      stop 1
   end subroutine seam_real

   !------------------------------------------------------------------
   ! seam_int(key, sval, ln, ival) - cast one staged value to integer;
   !              an unreadable string loud-stops named (STOP 1)
   subroutine seam_int(key, sval, ln, ival)
      character(len=*), intent(in) :: key, sval
      integer, intent(in) :: ln
      integer, intent(out) :: ival
      integer :: ios
      read (sval, *, iostat=ios) ival
      if (ios == 0) return
      write (0, '(a)') 'seam_cast: key '//trim(key)//trim(seam_tag(ln))// &
                       ': bad integer value "'//trim(sval)//'"'
      write (0, '(a)') 'seam_cast: loud-stop (the member parameter is not loaded)'
      stop 1
   end subroutine seam_int

   ! seam_tag(ln) - the source-line suffix of a message (a staged member key
   !                always arrives from a data line, ln >= 1; the ln <= 0 arm
   !                is the defensive default-tag shape)
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
