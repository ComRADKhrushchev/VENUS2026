module param_mapping
  implicit none
  private

  public :: map_elec_method, map_surface_model, map_task
  public :: map_integrator, map_init_sampling
  public :: map_old_nsurf

  integer, parameter :: MAX_ERR_LEN = 512

contains

  ! ---------------------------------------------------------------------------
  ! Case-insensitive string comparison
  ! ---------------------------------------------------------------------------
  function str_eq(a, b) result(ok)
    character(len=*), intent(in) :: a, b
    logical :: ok
    integer :: i, ca, cb, la, lb
    la = len_trim(a); lb = len_trim(b)
    if (la /= lb) then; ok = .false.; return; end if
    do i = 1, la
       ca = iachar(a(i:i)); cb = iachar(b(i:i))
       if (ca >= iachar('a') .and. ca <= iachar('z')) ca = ca - 32
       if (cb >= iachar('a') .and. cb <= iachar('z')) cb = cb - 32
       if (ca /= cb) then; ok = .false.; return; end if
    end do
    ok = .true.
  end function str_eq

  ! ---------------------------------------------------------------------------
  ! ELEC_METHOD  string -> CALTYP integer
  ! ---------------------------------------------------------------------------
  subroutine map_elec_method(str, caltyp, iostat, errmsg)
    character(len=*), intent(in)  :: str
    integer,          intent(out) :: caltyp, iostat
    character(len=*), intent(out) :: errmsg

    if (str_eq(str, 'ADIABATIC')) then
       caltyp = -1; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'TDHF-FSSH') .or. str_eq(str, 'IESH') &
            .or. str_eq(str, 'MDEF')) then
       ! Non-adiabatic methods removed from this minimal TEST build — the
       ! full implementations live in the internal dev repository. Reaching
       ! here means the caller asked for physics this build cannot provide;
       ! stop loudly rather than silently fall back to classical dynamics.
       caltyp = -1; iostat = 1
       errmsg = 'ERROR: ELEC_METHOD='//trim(str)//' is a non-adiabatic method ' &
                //'that has been removed from this build. Only ADIABATIC ' &
                //'is supported.'
    else
       caltyp = -1; iostat = 1
       errmsg = 'ERROR: ELEC_METHOD must be ADIABATIC (TDHF-FSSH/IESH/MDEF are removed)'
    end if
  end subroutine map_elec_method

  ! ---------------------------------------------------------------------------
  ! SURFACE_MODEL  string -> NSURF integer (new numbering: 0=NONE,1=RELAXED,2=RIGID)
  ! ---------------------------------------------------------------------------
  subroutine map_surface_model(str, nsurf, iostat, errmsg)
    character(len=*), intent(in)  :: str
    integer,          intent(out) :: nsurf, iostat
    character(len=*), intent(out) :: errmsg

    if (str_eq(str, 'NONE')) then
       nsurf = 0; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'RELAXED')) then
       nsurf = 1; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'RIGID')) then
       nsurf = 2; iostat = 0; errmsg = ''
    else
       nsurf = 0; iostat = 1
       errmsg = 'ERROR: SURFACE_MODEL must be NONE, RELAXED, or RIGID'
    end if
  end subroutine map_surface_model

  ! ---------------------------------------------------------------------------
  ! TASK  string -> NSELT integer
  ! ---------------------------------------------------------------------------
  subroutine map_task(str, nselt, iostat, errmsg)
    character(len=*), intent(in)  :: str
    integer,          intent(out) :: nselt, iostat
    character(len=*), intent(out) :: errmsg

    if (str_eq(str, 'NORMAL-MODE')) then
       nselt = -1; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'READ-QP')) then
       nselt = 0; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'TRAJECTORY')) then
       nselt = 2; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'BARRIER')) then
       nselt = 3; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'PES-SCAN')) then
       nselt = 4; iostat = 0; errmsg = ''
    else
       nselt = 0; iostat = 1
       ! 已移除的旧值：STATIONARY-POINT(-3)/REACTION-PATH(-2)/MIN-ENERGY(1)
       errmsg = 'ERROR: TASK must be NORMAL-MODE, READ-QP, TRAJECTORY, &
                &BARRIER, or PES-SCAN'
    end if
  end subroutine map_task

  ! ---------------------------------------------------------------------------
  ! INTEGRATOR  string -> INTEGRATOR + LLL integers
  ! ---------------------------------------------------------------------------
  subroutine map_integrator(str, integrator, lll, iostat, errmsg)
    character(len=*), intent(in)  :: str
    integer,          intent(out) :: integrator, lll, iostat
    character(len=*), intent(out) :: errmsg

    if (str_eq(str, 'VERLET')) then
       integrator = 3; lll = 1; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'BEEMAN')) then
       integrator = 3; lll = 0; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'RK4')) then
       integrator = 1; lll = 0; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'RADAU-FIXED')) then
       integrator = 1; lll = -1; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'RADAU-ADAPTIVE')) then
       integrator = 1; lll = 6; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'SYMPLECTIC-4')) then
       integrator = 2; lll = 4; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'SYMPLECTIC-6')) then
       integrator = 2; lll = 6; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'SYMPLECTIC-8')) then
       integrator = 2; lll = 8; iostat = 0; errmsg = ''
    else
       integrator = -1; lll = -1; iostat = 1
       errmsg = 'ERROR: INTEGRATOR must be VERLET, BEEMAN, RK4, RADAU-FIXED, &
                &RADAU-ADAPTIVE, SYMPLECTIC-4, SYMPLECTIC-6, or SYMPLECTIC-8'
    end if
  end subroutine map_integrator

  ! ---------------------------------------------------------------------------
  ! INIT_SAMPLING_A / INIT_SAMPLING_B  string -> NACTA / NACTB integer
  ! ---------------------------------------------------------------------------
  subroutine map_init_sampling(str, nact, iostat, errmsg)
    character(len=*), intent(in)  :: str
    integer,          intent(out) :: nact, iostat
    character(len=*), intent(out) :: errmsg

    if (str_eq(str, 'MB')) then
       nact = 0; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'ORTHANT')) then
       nact = 1; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'MICROCANONICAL')) then
       nact = 2; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'NORMAL-MODE')) then
       nact = 3; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'LOCAL-MODE')) then
       nact = 4; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'BOLTZMANN-VIB')) then
       nact = 5; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'FIXED-ENERGY')) then
       nact = 6; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'MD')) then
       nact = 7; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'QM-MICRO')) then
       nact = 8; iostat = 0; errmsg = ''
    elseif (str_eq(str, 'CI-QM-MICRO')) then
       nact = 9; iostat = 0; errmsg = ''
    else
       nact = 0; iostat = 1
       errmsg = 'ERROR: INIT_SAMPLING must be MB, ORTHANT, MICROCANONICAL, &
                &NORMAL-MODE, LOCAL-MODE, BOLTZMANN-VIB, FIXED-ENERGY, &
                &MD, QM-MICRO, or CI-QM-MICRO'
    end if
  end subroutine map_init_sampling

  ! ---------------------------------------------------------------------------
  ! Remap old NSURF integer values to new numbering: 2->1, 3->2
  ! ---------------------------------------------------------------------------
  pure integer function map_old_nsurf(old_nsurf)
    integer, intent(in) :: old_nsurf
    select case (old_nsurf)
    case (2)
       map_old_nsurf = 1   ! old RELAXED -> new RELAXED (merged with ER)
    case (3)
       map_old_nsurf = 2   ! old RIGID -> new RIGID
    case default
       map_old_nsurf = old_nsurf
    end select
  end function map_old_nsurf

end module param_mapping
