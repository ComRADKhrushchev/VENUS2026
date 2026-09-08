! config.f90 - representation-layer configuration: observables + the four input-provision groups
!              derived by freezing (type-encapsulated + blanks left)
module config
  implicit none
  private
  public :: propagator_t, electronic_t, reactants_t, event_t, observables_t
  public :: propagator, electronic, reactants, event, observables

  ! ===== Input provisions derived by freezing (spec §3.1-3.2) =====

  ! Nuclear propagator (spec §3.2-1: bounded functionality, frozen into a closed family)
  type :: propagator_t
     integer :: id    = 1     ! integrator code (the closed three-member set verlet/symple/radau)
     integer :: order = 0     ! order/variant parameter
     real(8) :: dt    = 0.01d0 ! step [10 fs]
     ! Left blank: adaptive step bounds dt_min/dt_max (introduced when RADAU is concretized, C6)
  end type

  ! Electronic layer (spec §3.2-2: two independent axes - H representation capability x
  ! electron propagation law)
  type :: electronic_t
     integer :: method = 1       ! electron propagation law (1 = adiabatic degenerate point; nonadiabatic reserved)
     integer :: n_surf  = 1      ! H representation capability (single surface; multi-surface/coupling waits for nonadiabatic)
     character(len=32) :: pes_main = ''  ! main-PES member name (H-field provider, interface A)
     ! Left blank: additive force member table (GLO springs / Langevin ghost / LDFA -
     ! introduced with the facade in C2);
     ! dielectric / external-field and other capability extensions
  end type

  ! Reactant provisions (spec §3.2-3: configuration provisions + initial-state statistical
  ! distributions, both as input data)
  type :: reactants_t
     integer :: n_frag = 1           ! number of fragments
     integer :: surface_model = 0    ! configuration paradigm (0 = gas phase; 1/2 = the two surface paradigms)
     integer :: dist_scheme(2) = 0   ! per-fragment sampling scheme (pointing to dist_* members)
     real(8) :: e_rel  = 0.0d0       ! collision energy (relative-motion distribution) [kcal/mol]
     real(8) :: b_max  = 0.0d0       ! maximum impact parameter [Å]
     real(8) :: r_sep  = 0.0d0       ! initial separation [Å]
     ! Configuration details (roster/geometry/unit cell) are carried by config_roster
     ! (K2 skeleton);
     ! Left blank: orientation prescriptions and per-fragment distribution parameters -
     ! owned by the dist_* members (member = code + data)
  end type

  ! Event provisions (spec §3.2-3: termination predicate + outcome classification, data not code)
  type :: event_t
     character(len=256)  :: term_expr   = '' ! termination-predicate expression (wired in C3)
     character(len=1024) :: outcome_tbl = '' ! outcome-classification table (wired in C7)
     ! Left blank: event-class extensions (constrained windows / coupled multi-events - added
     ! when the predicate grammar is finalized)
  end type

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
  type(event_t)      :: event
  type(observables_t) :: observables
end module
