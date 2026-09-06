module input_parser
  implicit none
  private

  integer, parameter :: MAX_KEY_LEN = 64
  integer, parameter :: MAX_LINE_LEN = 2048
  integer, parameter :: MAX_VALUES = 5000

  type :: TokenPair
     character(len=MAX_KEY_LEN) :: keyword = ''
     character(len=MAX_LINE_LEN) :: values_str = ''
     type(TokenPair), pointer :: next => null()
  end type TokenPair

  type :: ParameterDict
     type(TokenPair), pointer :: head => null()
     integer :: n_entries = 0
  end type ParameterDict

  type(ParameterDict), save :: store
  logical, save :: format_is_new = .false.
  logical, save :: input_loaded = .false.
  integer, save :: line_number = 0

  ! Public interface
  public :: load_input, is_keyword_format
  public :: get_int, get_real, get_str
  public :: get_int_arr, get_real_arr
  public :: get_real_element, get_int_element
  public :: has_keyword, count_keyword_values
  public :: set_default
  public :: MAX_LINE_LEN

  interface get_int
     module procedure get_int_scalar
  end interface

  interface get_real
     module procedure get_real_scalar
  end interface

  interface get_str
     module procedure get_str_scalar
  end interface

contains

  ! ---------------------------------------------------------------------------
  ! Uppercase a string (in-place)
  ! ---------------------------------------------------------------------------
  function to_upper(str) result(upper)
    character(len=*), intent(in) :: str
    character(len=len_trim(str)) :: upper
    integer :: i, ch
    upper = trim(adjustl(str))
    do i = 1, len_trim(upper)
       ch = iachar(upper(i:i))
       if (ch >= iachar('a') .and. ch <= iachar('z')) then
          upper(i:i) = achar(ch - 32)
       end if
    end do
  end function to_upper

  ! ---------------------------------------------------------------------------
  ! Load input file and detect format
  ! ---------------------------------------------------------------------------
  subroutine load_input(filename, unit_num)
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: unit_num
    character(len=MAX_LINE_LEN) :: line, kw, vals
    integer :: u, ios, eqpos
    logical :: file_exists, first_data_line
    type(TokenPair), pointer :: current, newpair

    u = 10
    if (present(unit_num)) u = unit_num

    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
       write(0,*) 'INPUT_PARSER: ERROR: file not found: ', trim(filename)
       stop
    end if

    open(u, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(0,*) 'INPUT_PARSER: ERROR: cannot open: ', trim(filename)
       stop
    end if

    ! Initialise store
    if (associated(store%head)) call destroy_store()
    store%n_entries = 0
    format_is_new = .false.
    first_data_line = .true.
    line_number = 0

    do
       read(u, '(A)', iostat=ios) line
       if (ios < 0) exit  ! EOF
       if (ios > 0) then
          write(0,*) 'INPUT_PARSER: ERROR reading line ', line_number + 1
          stop
       end if
       line_number = line_number + 1

       ! Skip blank lines and comments
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#' .or. line(1:1) == '!') cycle

       ! Detect format on first data line
       if (first_data_line) then
          first_data_line = .false.
          eqpos = scan(line, '=')
          if (eqpos > 0) then
             format_is_new = .true.
          else
             format_is_new = .false.
             exit  ! old format — caller will re-read sequentially
          end if
       end if

       if (.not. format_is_new) exit

       ! Parse keyword = value line
       eqpos = scan(line, '=')
       if (eqpos == 0) then
          write(0,*) 'INPUT_PARSER: WARNING: line ', line_number, &
               ' has no "=" — skipping: ', trim(line)
          cycle
       end if

       kw = adjustl(line(1:eqpos-1))
       vals = adjustl(line(eqpos+1:))

       ! Handle continuation: if vals ends with ',' absorb next line(s)
       do while (len_trim(vals) > 0)
          if (vals(len_trim(vals):len_trim(vals)) == ',') then
             read(u, '(A)', iostat=ios) line
             if (ios /= 0) exit
             line_number = line_number + 1
             vals = trim(vals) // ' ' // trim(adjustl(line))
          else
             exit
          end if
       end do

       call append_pair(kw, vals)
    end do

    close(u)
    input_loaded = .true.
  end subroutine load_input

  ! ---------------------------------------------------------------------------
  ! Check if input file is in new keyword format
  ! ---------------------------------------------------------------------------
  function is_keyword_format() result(flag)
    logical :: flag
    flag = format_is_new
  end function is_keyword_format

  ! ---------------------------------------------------------------------------
  ! Append a keyword-value pair to the linked list
  ! ---------------------------------------------------------------------------
  subroutine append_pair(keyword, values_str)
    character(len=*), intent(in) :: keyword, values_str
    type(TokenPair), pointer :: newpair

    allocate(newpair)
    newpair%keyword = to_upper(keyword)
    newpair%values_str = trim(values_str)
    newpair%next => store%head
    store%head => newpair
    store%n_entries = store%n_entries + 1
  end subroutine append_pair

  ! ---------------------------------------------------------------------------
  ! Find a keyword in the store (case-insensitive); return pointer or null()
  ! ---------------------------------------------------------------------------
  function find_keyword(key) result(ptr)
    character(len=*), intent(in) :: key
    type(TokenPair), pointer :: ptr
    character(len=MAX_KEY_LEN) :: upper_key

    upper_key = to_upper(key)
    ptr => store%head
    do while (associated(ptr))
       if (trim(ptr%keyword) == trim(upper_key)) return
       ptr => ptr%next
    end do
  end function find_keyword

  ! ---------------------------------------------------------------------------
  ! Extract the n-th value from a comma-separated string
  ! ---------------------------------------------------------------------------
  subroutine extract_nth(values_str, nth, value_out, found)
    character(len=*), intent(in) :: values_str
    integer, intent(in) :: nth
    character(len=MAX_LINE_LEN), intent(out) :: value_out
    logical, intent(out) :: found
    integer :: start_pos, comma_pos, count, i, strlen
    character(len=MAX_LINE_LEN) :: temp

    value_out = ''
    found = .false.
    count = 0
    start_pos = 1
    strlen = len_trim(values_str)
    temp = ''

    do while (start_pos <= strlen)
       ! Skip leading spaces
       do while (start_pos <= strlen)
          if (values_str(start_pos:start_pos) /= ' ') exit
          start_pos = start_pos + 1
       end do
       if (start_pos > strlen) exit

       comma_pos = index(values_str(start_pos:), ',')
       if (comma_pos == 0) then
          temp = values_str(start_pos:strlen)
          count = count + 1
          if (count == nth) then
             value_out = trim(adjustl(temp))
             found = .true.
          end if
          exit
       else
          temp = values_str(start_pos:start_pos+comma_pos-2)
          count = count + 1
          if (count == nth) then
             value_out = trim(adjustl(temp))
             found = .true.
             return
          end if
          start_pos = start_pos + comma_pos
       end if
    end do
  end subroutine extract_nth

  ! ---------------------------------------------------------------------------
  ! Count comma-separated values
  ! ---------------------------------------------------------------------------
  function count_values(values_str) result(n)
    character(len=*), intent(in) :: values_str
    integer :: n, pos, strlen, i
    logical :: in_value

    n = 0
    strlen = len_trim(values_str)
    if (strlen == 0) return

    in_value = .false.
    do i = 1, strlen
       if (values_str(i:i) == ',') then
          n = n + 1
          in_value = .false.
       else if (values_str(i:i) /= ' ') then
          in_value = .true.
       end if
    end do
    if (in_value) n = n + 1
  end function count_values

  ! ---------------------------------------------------------------------------
  ! Get scalar integer
  ! ---------------------------------------------------------------------------
  function get_int_scalar(key, default) result(val)
    character(len=*), intent(in) :: key
    integer, intent(in) :: default
    integer :: val
    type(TokenPair), pointer :: ptr
    character(len=MAX_LINE_LEN) :: value_str
    logical :: found
    integer :: ios

    ptr => find_keyword(key)
    if (.not. associated(ptr)) then
       val = default
       return
    end if

    call extract_nth(ptr%values_str, 1, value_str, found)
    if (.not. found) then
       val = default
       return
    end if

    read(value_str, *, iostat=ios) val
    if (ios /= 0) val = default
  end function get_int_scalar

  ! ---------------------------------------------------------------------------
  ! Get scalar double precision
  ! ---------------------------------------------------------------------------
  function get_real_scalar(key, default) result(val)
    character(len=*), intent(in) :: key
    double precision, intent(in) :: default
    double precision :: val
    type(TokenPair), pointer :: ptr
    character(len=MAX_LINE_LEN) :: value_str
    logical :: found
    integer :: ios

    ptr => find_keyword(key)
    if (.not. associated(ptr)) then
       val = default
       return
    end if

    call extract_nth(ptr%values_str, 1, value_str, found)
    if (.not. found) then
       val = default
       return
    end if

    read(value_str, *, iostat=ios) val
    if (ios /= 0) val = default
  end function get_real_scalar

  ! ---------------------------------------------------------------------------
  ! Get n-th element from a keyword's values (1-based)
  ! ---------------------------------------------------------------------------
  function get_real_element(key, nth, default) result(val)
    character(len=*), intent(in) :: key
    integer, intent(in) :: nth
    double precision, intent(in) :: default
    double precision :: val
    type(TokenPair), pointer :: ptr
    character(len=MAX_LINE_LEN) :: value_str
    logical :: found
    integer :: ios

    ptr => find_keyword(key)
    if (.not. associated(ptr)) then
       val = default
       return
    end if

    call extract_nth(ptr%values_str, nth, value_str, found)
    if (.not. found) then
       val = default
       return
    end if

    read(value_str, *, iostat=ios) val
    if (ios /= 0) val = default
  end function get_real_element

  ! ---------------------------------------------------------------------------
  ! Get string
  ! ---------------------------------------------------------------------------
  function get_str_scalar(key, default) result(val)
    character(len=*), intent(in) :: key, default
    character(len=MAX_LINE_LEN) :: val
    type(TokenPair), pointer :: ptr
    character(len=MAX_LINE_LEN) :: value_str
    logical :: found

    ptr => find_keyword(key)
    if (.not. associated(ptr)) then
       val = default
       return
    end if

    val = trim(ptr%values_str)
    if (len_trim(val) == 0) val = default
  end function get_str_scalar

  ! ---------------------------------------------------------------------------
  ! Get real array — fill arr(1:n) from keyword values
  ! ---------------------------------------------------------------------------
  subroutine get_real_arr(key, arr, n)
    character(len=*), intent(in) :: key
    integer, intent(in) :: n
    double precision, intent(out) :: arr(n)
    type(TokenPair), pointer :: ptr
    character(len=MAX_LINE_LEN) :: value_str
    logical :: found
    integer :: i, ios

    ptr => find_keyword(key)
    if (.not. associated(ptr)) then
       arr(:) = 0.0d0
       return
    end if

    do i = 1, n
       call extract_nth(ptr%values_str, i, value_str, found)
       if (.not. found) then
          arr(i) = 0.0d0
       else
          read(value_str, *, iostat=ios) arr(i)
          if (ios /= 0) arr(i) = 0.0d0
       end if
    end do
  end subroutine get_real_arr

  ! ---------------------------------------------------------------------------
  ! Get integer array — fill arr(1:n) from keyword values
  ! ---------------------------------------------------------------------------
  subroutine get_int_arr(key, arr, n)
    character(len=*), intent(in) :: key
    integer, intent(in) :: n
    integer, intent(out) :: arr(n)
    type(TokenPair), pointer :: ptr
    character(len=MAX_LINE_LEN) :: value_str
    logical :: found
    integer :: i, ios

    ptr => find_keyword(key)
    if (.not. associated(ptr)) then
       arr(:) = 0
       return
    end if

    do i = 1, n
       call extract_nth(ptr%values_str, i, value_str, found)
       if (.not. found) then
          arr(i) = 0
       else
          read(value_str, *, iostat=ios) arr(i)
          if (ios /= 0) arr(i) = 0
       end if
    end do
  end subroutine get_int_arr

  ! ---------------------------------------------------------------------------
  ! Clean up linked list
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------
  ! Get n-th integer element from a keyword's values (1-based)
  ! ---------------------------------------------------------------------------
  function get_int_element(key, nth, default) result(val)
    character(len=*), intent(in) :: key
    integer, intent(in) :: nth
    integer, intent(in) :: default
    integer :: val
    type(TokenPair), pointer :: ptr
    character(len=MAX_LINE_LEN) :: value_str
    logical :: found
    integer :: ios

    ptr => find_keyword(key)
    if (.not. associated(ptr)) then
       val = default
       return
    end if

    call extract_nth(ptr%values_str, nth, value_str, found)
    if (.not. found) then
       val = default
       return
    end if

    read(value_str, *, iostat=ios) val
    if (ios /= 0) val = default
  end function get_int_element

  ! ---------------------------------------------------------------------------
  ! Check if a keyword exists in the store
  ! ---------------------------------------------------------------------------
  function has_keyword(key) result(found)
    character(len=*), intent(in) :: key
    logical :: found
    type(TokenPair), pointer :: ptr
    ptr => find_keyword(key)
    found = associated(ptr)
  end function has_keyword

  ! ---------------------------------------------------------------------------
  ! Return number of comma-separated values for a given keyword
  ! ---------------------------------------------------------------------------
  function count_keyword_values(key) result(n)
    character(len=*), intent(in) :: key
    integer :: n
    type(TokenPair), pointer :: ptr
    ptr => find_keyword(key)
    if (.not. associated(ptr)) then
       n = 0
    else
       n = count_values(ptr%values_str)
    end if
  end function count_keyword_values

  ! ---------------------------------------------------------------------------
  ! Add a default value for a keyword — only if not already set by user
  ! ---------------------------------------------------------------------------
  subroutine set_default(key, value_str)
    character(len=*), intent(in) :: key, value_str
    type(TokenPair), pointer :: ptr, newpair

    ptr => find_keyword(key)
    if (associated(ptr)) return   ! user already set this key, do not override

    allocate(newpair)
    newpair%keyword = trim(adjustl(key))
    newpair%values_str = trim(adjustl(value_str))
    newpair%next => store%head
    store%head => newpair
    store%n_entries = store%n_entries + 1
  end subroutine set_default

  ! ---------------------------------------------------------------------------
  ! Clean up linked list
  ! ---------------------------------------------------------------------------
  subroutine destroy_store()
    type(TokenPair), pointer :: current, nextpair

    current => store%head
    do while (associated(current))
       nextpair => current%next
       deallocate(current)
       current => nextpair
    end do
    store%head => null()
    store%n_entries = 0
  end subroutine destroy_store

end module input_parser
