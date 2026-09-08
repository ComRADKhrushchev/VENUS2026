! control.f90 - control-flow must variables (standalone module): driver-orchestration input,
!              no physics settings
module control
  ! Control-flow input: what the control kernel needs for orchestration, no physics settings
  ! (spec §2.2: the control kernel holds no physics knowledge)
  implicit none
  public
  integer :: task    = 0   ! task scheme (trajectory run / read-in coordinates / ...; maps the TASK input)
  integer :: n_traj  = 1   ! number of trajectories (old NT)
  character(len=80) :: title(2) = ''  ! run title
  integer :: i_seed(8) = 0 ! random seed (8 words; used to initialize rng)
  ! Left blank: checkpoint restart, parallelism/process management and other run controls -
  ! introduced later if needed (spec §1 scope conditions)
end module
