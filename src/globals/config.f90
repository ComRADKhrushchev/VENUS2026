! config.f90 - representation-layer configuration: observables + the four input-provision groups
!              derived by freezing (type-encapsulated + blanks left)
module config
  implicit none
  private
  public :: propagator_t, electronic_t, reactants_t, observables_t
  public :: propagator, electronic, reactants, observables

  ! ===== Input provisions derived by freezing (spec §3.1-3.2) =====

  ! Nuclear propagator (spec §3.2-1: bounded functionality, frozen into a closed family)
  type :: propagator_t
     integer :: id    = 1     ! integrator code (the closed three-member set verlet/symple/radau)
     integer :: order = 0     ! order/variant parameter
     real(8) :: dt    = 0.01d0 ! step [10 fs]
     real(8) :: dt_min = 0.0d0 ! adaptive sequence-length lower bound [10 fs] (magnitude
                               ! bound; 0 = unclamped - consumed by the adaptive member only)
     real(8) :: dt_max = 0.0d0 ! adaptive sequence-length upper bound [10 fs] (magnitude
                               ! bound; 0 = unclamped - consumed by the adaptive member only)
  end type

  ! Electronic layer (spec §3.2-2: two independent axes - H representation capability x
  ! electron propagation law)
  type :: electronic_t
     character(len=32) :: method = 'adiabatic'  ! electron propagation law PACKAGE name (assembly key,
                                                ! same form as pes_main - E13; E14: a name-string
                                                ! lookup into the methods/ library registry - the
                                                ! method_tbl of elec_facade, filled by the thin
                                                ! reg_*.f90 files (spec §5); 'adiabatic' = the
                                                ! all-no-op degenerate seed package (K2);
                                                ! 'ehrenfest'/'fssh'/'iesh' = phase-C arrivals;
                                                ! two exposure levels (spec §5): the whole package
                                                ! binds by this config key (level 1), reusable
                                                ! components are consumed by use (level 2))
     integer :: n_surf  = 1      ! H representation capability (single surface; multi-surface/coupling waits for nonadiabatic)
     character(len=32) :: pes_main = ''  ! main-PES member name (H-field provider, interface A)
     ! Left blank: additive force member table (GLO springs / Langevin ghost / LDFA -
     ! introduced with the facade in C2);
     ! dielectric / external-field and other capability extensions
  end type

  ! Reactant provisions (spec §3.2-3: configuration provisions + initial-state statistical
  ! distributions; E14: the configuration half is no longer wholly input data - composition
  ! arrives from the container folder and surface_model is derived, spec §6)
  type :: reactants_t
     integer :: n_frag = 1           ! number of fragments (= size(roster%frag) after assembly)
     character(len=128) :: system_dir = '' ! container folder path scanned by sysdef_load (E14 spec §6;
                                           ! grammar + default value finalized in C4)
     integer :: surface_model = 0    ! configuration paradigm (0 = gas phase / molecular; 1/2 = the two
                                     ! surface paradigms). E14: DERIVED, not an input key - sysdef_load
                                     ! writes it by the spec §6 cell-presence rule (any system-
                                     ! definition file carries a cell -> surface paradigm (1/2); none
                                     ! -> 0); the 1-vs-2 surface-paradigm selection is not derivable
                                     ! from cell presence alone and arrives with further input (C4)
     integer, allocatable :: dist_scheme(:) ! per-fragment sampling scheme (pointing to dist_* members;
                                            ! E15: arrayed together with roster%frag - allocated at
                                            ! input/assembly, size = n_frag; the K1 max_frag risk-note
                                            ! extension direction landed early at the review gate)
     real(8) :: e_rel  = 0.0d0       ! collision energy (relative-motion distribution) [kcal/mol]
     real(8) :: b_max  = 0.0d0       ! maximum impact parameter [Å]
     real(8) :: r_sep  = 0.0d0       ! initial separation [Å]
     ! Configuration details (roster/geometry/unit cell) are carried by config_roster
     ! (K2 skeleton);
     ! Left blank: orientation prescriptions and per-fragment distribution parameters -
     ! owned by the dist_* members (member = code + data)
  end type

  ! (E18/E19: the event provisions type is WITHDRAWN - termination is the container's own
  ! judgment (container_bind_term, force_facade) and classification likewise
  ! (container_bind_classify, final_state); neither event grammar key remains)

  ! ===== Observables (generic energy terms and recording settings) =====
  type :: observables_t
     integer :: rec_level = 2          ! recording verbosity (preset column groups)
     character(len=256) :: obs_list = '' ! observable column list (C4 recorder column registry)
     logical :: e_kin = .true.         ! kinetic energy T (generic energy term)
     logical :: e_pot = .true.         ! potential energy V (generic energy term)
     logical :: e_tot = .true.         ! total energy H (generic energy term)
     ! Left blank: more observables extended via recorder column registration (angular
     ! momenta, local quantities, etc., C4/C8)
  end type

  type(propagator_t) :: propagator
  type(electronic_t) :: electronic
  type(reactants_t)  :: reactants
  type(observables_t) :: observables
end module
