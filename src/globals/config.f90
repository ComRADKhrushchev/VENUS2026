!=====================================================================
! config.f90 - representation-layer configuration: input-provision groups
!              (propagator/electronic/reactants) and observables
! Design:
!   Type-encapsulated runtime configuration, one module singleton per group.
!   Field comments carry the input-key contracts; unspecified extensions are
!   left blank by design.
!=====================================================================
module config
  implicit none
  private
  public :: propagator_t, electronic_t, reactants_t, observables_t
  public :: propagator, electronic, reactants, observables

  ! ===== Input provisions =====

  ! Nuclear propagator (bounded functionality, frozen into a closed family)
  type :: propagator_t
     integer :: id    = 1     ! integrator code (the closed three-member set verlet/symple/radau)
     integer :: order = 0     ! order/variant parameter
     real(8) :: dt    = 0.01d0 ! step [10 fs]
     real(8) :: dt_min = 0.0d0 ! adaptive sequence-length lower bound [10 fs] (magnitude
                               ! bound; 0 = unclamped - consumed by the adaptive member only)
     real(8) :: dt_max = 0.0d0 ! adaptive sequence-length upper bound [10 fs] (magnitude
                               ! bound; 0 = unclamped - consumed by the adaptive member only)
  end type

  ! Electronic layer (two independent axes - H representation capability x
  ! electron propagation law)
  type :: electronic_t
     character(len=32) :: method = 'adiabatic'  ! electron propagation law package name (assembly
                                                ! key - a name-string lookup into the methods/
                                                ! library registry; 'adiabatic' = the all-no-op
                                                ! degenerate seed package; the package binds by
                                                ! this key, reusable components by use)
     integer :: n_surf  = 1      ! H representation capability (single surface; multi-surface/coupling waits for nonadiabatic)
     character(len=32) :: pes_main = ''  ! main-PES member name (H-field provider, interface A)
     ! Left blank: additive force member table (GLO springs / Langevin ghost / LDFA);
     ! dielectric / external-field and other capability extensions
  end type

  ! Reactant provisions (configuration provisions + initial-state statistical
  ! distributions; composition arrives from the container folder and surface_model
  ! is derived)
  type :: reactants_t
     integer :: n_frag = 1           ! number of fragments (= size(list_atoms%frag) after assembly)
     character(len=128) :: system_dir = '' ! container folder path scanned by sysdef_load
     integer :: surface_model = 0    ! configuration paradigm (0 = gas phase / molecular; 1/2 = the
                                     ! two surface paradigms) - DERIVED, not an input key:
                                     ! sysdef_load writes it by the cell-presence rule (any
                                     ! system-definition file carries a cell -> surface
                                     ! paradigm; none -> 0); the 1-vs-2 selection is not
                                     ! derivable from cell presence alone
     integer, allocatable :: dist_scheme(:) ! per-fragment sampling scheme (pointing to dist_* members;
                                            ! arrayed together with list_atoms%frag - allocated at
                                            ! input/assembly, size = n_frag)
     real(8) :: e_rel  = 0.0d0       ! collision energy (relative-motion distribution) [kcal/mol]
     real(8) :: b_max  = 0.0d0       ! maximum impact parameter [Å]
     real(8) :: r_sep  = 0.0d0       ! initial separation [Å]
     logical :: random_orient = .false. ! universal orientation stage switch (RANDOM_ORIENT;
                                       ! off = the archived five-case stream: the initial state
                                       ! keeps whatever orientation its member wrote; on = after
                                       ! the per-fragment members and before the incident
                                       ! channel, every nat>=2 fragment not served by the
                                       ! 'rotation' word gets one isotropic rot_euler pass
                                       ! about its COM - q and p together)
     integer :: bath = 0              ! bath provision selector (BATH enum code; 1 = none,
                                       ! 2 = Andersen velocity-reset species; 0 = the unset
                                       ! default, treated as none) - arms the pre-evolution
                                       ! bath loop inside samp_run (members realize -> bath
                                       ! loop -> orientation -> incident)
     integer :: n_bath = 0            ! bath-loop step count (N_BATH; 0 = the armed bath runs
                                       ! no steps - no collisions, no propagation)
     real(8) :: dt_bath = 0.0d0      ! bath-loop step size [10 fs] (DT_BATH; must be positive
                                       ! when N_BATH > 0)
     integer :: spectrum_src = 1     ! spectrum producer selector (SPECTRUM_SOURCE enum;
                                       ! 1 = COMPUTE - the internal Hessian derivation, lazy on
                                       ! first pull at the member seam; 2 = MANUAL - the
                                       ! spectrum-table file is the only producer; the
                                       ! container export channel is withdrawn - 2026-10-05
                                       ! spectrum-source legislation, docs/plans/)
     character(len=128) :: spectrum_file = '' ! spectrum-table file (SPECTRUM_FILE; MANUAL's
                                       ! carrier records - composition-keyed mode tables
                                       ! [cm^-1 face] + optional diatomic constants; a staged
                                       ! path under COMPUTE is a formal input error)

     ! Configuration details (list_atoms/geometry/unit cell) are carried by config_atoms;
     ! Left blank: per-fragment distribution parameters - owned by the dist_* members
     ! (the orientation prescription is the sampler stage above, not a member parameter;
     ! the bath species parameters T_BATH/NU_COLL are bath_andersen's seam keys)
  end type

  ! (the event provisions type is WITHDRAWN - termination is the container's own
  ! judgment (container_bind_term) and classification likewise (container_bind_classify,
  ! final_state); no event grammar key remains)

  ! ===== Observables (generic energy terms and recording settings) =====
  type :: observables_t
     integer :: rec_level = 2          ! recording verbosity (preset column groups)
     character(len=256) :: obs_list = '' ! observable column list (recorder column registry)
     logical :: e_kin = .true.         ! kinetic energy T (generic energy term)
     logical :: e_pot = .true.         ! potential energy V (generic energy term)
     logical :: e_tot = .true.         ! total energy H (generic energy term)
     ! Left blank: more observables extended via recorder column registration (angular
     ! momenta, local quantities, etc.)
  end type

  type(propagator_t) :: propagator
  type(electronic_t) :: electronic
  type(reactants_t)  :: reactants
  type(observables_t) :: observables
end module
