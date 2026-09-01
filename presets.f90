module presets
  use input_parser, only: set_default
  implicit none
  private

  public :: apply_preset

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
  ! Dispatch preset by model name (case-insensitive)
  ! ---------------------------------------------------------------------------
  subroutine apply_preset(model_name, iostat, errmsg)
    character(len=*), intent(in)  :: model_name
    integer,          intent(out) :: iostat
    character(len=*), intent(out) :: errmsg

    if (str_eq(model_name, 'GAS-PHASE')) then
       call preset_gas_phase()
       iostat = 0; errmsg = ''
    elseif (str_eq(model_name, 'RELAXED-SURFACE')) then
       call preset_relaxed_surface()
       iostat = 0; errmsg = ''
    elseif (str_eq(model_name, 'RIGID-SURFACE')) then
       call preset_rigid_surface()
       iostat = 0; errmsg = ''
    elseif (str_eq(model_name, 'FULL-SURFACE')) then
       call preset_full_surface()
       iostat = 0; errmsg = ''
    elseif (str_eq(model_name, 'GLO-SURFACE')) then
       call preset_glo_surface()
       iostat = 0; errmsg = ''
    else
       iostat = 1
       errmsg = 'ERROR: MODEL must be GAS-PHASE, RELAXED-SURFACE, RIGID-SURFACE, FULL-SURFACE, or GLO-SURFACE'
    end if
  end subroutine apply_preset

  ! ---------------------------------------------------------------------------
  ! GAS-PHASE: no surface, gas-phase collision
  ! Sets behavioral params only; descriptive params (NATOMS, masses, QZA, etc.)
  ! must be specified by user.
  ! ---------------------------------------------------------------------------
  subroutine preset_gas_phase()
    call set_default('NSURF', '0')
    call set_default('NSELT', '2')
    call set_default('NACTA', '0')
    call set_default('NACTB', '0')
    call set_default('NTHTA', '-1')
    call set_default('NREL', '1')
    call set_default('INTEGRATOR', '3')
    call set_default('LLL', '1')
    call set_default('NIP', '1')
  end subroutine preset_gas_phase

  ! ---------------------------------------------------------------------------
  ! RELAXED-SURFACE: harmonic oscillator surface (3-atom model)
  ! Fragment A = surface oscillator (2 atoms: surface Au + anchor Au)
  ! Fragment B = incident atom
  ! NACTA=5: Boltzmann vibrational for the surface oscillator
  ! NACTB=0: Maxwell-Boltzmann for incident atom
  ! ---------------------------------------------------------------------------
  subroutine preset_relaxed_surface()
    call set_default('NSURF', '1')
    call set_default('NZDOWN', '1')
    call set_default('NSELT', '2')
    call set_default('NACTA', '5')
    call set_default('NACTB', '0')
    call set_default('NTHTA', '-1')
    call set_default('NRNDXY', '1')
    call set_default('NCHI', '1')
    call set_default('CHI', '0.0')
    call set_default('INTEGRATOR', '3')
    call set_default('LLL', '1')
    call set_default('NREL', '1')
    call set_default('NIP', '1')
    call set_default('NCROT', '200')
    call set_default('NCVIB', '0')
    call set_default('TRV_A', '300.0')
    call set_default('RMAX', '5.0')
    call set_default('RBAR', '3.5')
    call set_default('NOB', '0')
    call set_default('BMAX', '0.0')
    call set_default('NPATHS', '0')
  end subroutine preset_relaxed_surface

  ! ---------------------------------------------------------------------------
  ! RIGID-SURFACE: rigid surface atoms with infinite mass
  ! ---------------------------------------------------------------------------
  subroutine preset_rigid_surface()
    ! D1 fix (F24): declare the surface via the SURFACE_MODEL string path.
    ! The old NSURF='2' default alone went through the integer fallback,
    ! where map_old_nsurf interprets 2 in the OLD numbering (2=RELAXED) and
    ! remaps it to 1 — the RIGID-SURFACE preset silently yielded a relaxed
    ! surface. SURFACE_MODEL takes precedence over the NSURF integer in
    ! venus_input, so the preset now yields NSURF=2 as intended.
    call set_default('SURFACE_MODEL', 'RIGID')
    call set_default('NSURF', '2')
    call set_default('NSELT', '2')
    call set_default('NACTA', '0')
    call set_default('NACTB', '0')
    call set_default('NTHTA', '-1')
    call set_default('NRNDXY', '1')
    call set_default('INTEGRATOR', '3')
    call set_default('LLL', '1')
    call set_default('NREL', '1')
    call set_default('NIP', '1')
    call set_default('NCROT', '200')
    call set_default('RMAX', '5.0')
    call set_default('RBAR', '3.5')
  end subroutine preset_rigid_surface

  ! ---------------------------------------------------------------------------
  ! FULL-SURFACE: multi-atom surface slab with MD sampling
  ! ---------------------------------------------------------------------------
  subroutine preset_full_surface()
    call set_default('NSURF', '1')
    call set_default('NSELT', '2')
    call set_default('NACTA', '0')
    call set_default('NACTB', '7')
    call set_default('NTHTA', '-1')
    call set_default('NRNDXY', '1')
    call set_default('INTEGRATOR', '3')
    call set_default('LLL', '1')
    call set_default('NREL', '1')
    call set_default('NIP', '1')
    call set_default('NCROT', '200')
    call set_default('RMAX', '5.0')
    call set_default('RBAR', '3.5')
    call set_default('THERMOTEMP', '300.0')
  end subroutine preset_full_surface

  ! ---------------------------------------------------------------------------
  ! GLO-SURFACE: like RELAXED-SURFACE but atom 2 is a ghost Langevin oscillator
  ! with finite mass. Set FCG > 0 to enable Langevin friction + noise on atom 2.
  ! FCG = 0 (default) degenerates to RELAXED-SURFACE behaviour.
  ! ---------------------------------------------------------------------------
  subroutine preset_glo_surface()
    call preset_relaxed_surface()
    call set_default('FCG', '0.0')
  end subroutine preset_glo_surface

end module presets
