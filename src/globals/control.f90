!=====================================================================
! control.f90 - control-flow run variables (driver orchestration, no physics)
! Design:
!   What the control kernel needs for orchestration - the kernel holds no
!   physics knowledge. Left blank: checkpoint restart, parallelism/process
!   management and other run controls - introduced later if needed.
!=====================================================================
module control
  implicit none
  public
  integer :: task    = 0   ! task scheme (trajectory run / read-in coordinates / ...; maps the TASK input)
  integer :: n_traj  = 1   ! number of trajectories
 integer :: max_steps = 100000 ! safety-net step cap per trajectory [step] (the ever-present
                               ! MAX_STEPS net - the only termination judgment never
                               ! delegated to the container; grammar key MAX_STEPS=)
  character(len=80) :: title(2) = ''  ! run title
  integer :: i_seed(8) = 0 ! random seed (8 words; used to initialize rng)
end module
