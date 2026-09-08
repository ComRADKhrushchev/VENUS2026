!=====================================================================
! dist_incident.f90 - L3 example distribution member (B' registry family): incident-channel
!                     distribution - collision energy /
!                     impact parameter / orientation prescription and sampling (example 3/3;
!                     surface/gas two routes merged)
! Reference: src_VENUS/SURF.f:120-187 (orientation-angle prescription: polar angle sampled with
!       P(theta) ~ sin theta,
!       cos(theta) = 1-U(1-cos(theta_max)) :133-141, azimuth uniform on [0,2pi) :142-149,
!       drop point
!       uniform within the unit cell :153-169, aiming disk b = b_max*sqrt(U) :171-187 -
!       manual Eqs. V.19/
!       V.24-27), :211-217 (random cluster-box drop point), :258-325 (three collision-energy
!       prescriptions: thermal
!       GAMA(2) :261-266, fixed :278-292, density-weighted Gaussian rejection sampling :295-325),
!       :331-338 (incident-direction unit vector - pointing at the surface), :343-374 (relative
!       velocity folded by mass
!       into fragment momenta), :451-458 (DENS = v^3*exp(-((v-v0)/alpha)^2)),
!       src_VENUS/SELECT.f90:1300-1334/:1337-1404 (the basis for merging the surface/gas two
!       assembly routes
!       - already cited in sampler.f90)
! Registration: this member enters the sampler distribution member table (one samp_reg line -
! executed in the assembly phase; C5 wiring)
! Library home: methods/ (E14) — a dist_* member is a pre-stored generic distribution method;
! its samp_reg registration line lands at C5 together with its sample routine (K2 E3 slimming
! ruling); see methods/README.md
!=====================================================================
module dist_incident
   implicit none
   private
   public :: inc_params_t, inc_init

   ! Division of labor with config: the collision energy e_rel [kcal/mol], maximum impact
   ! parameter b_max [Å], and initial
   ! separation r_sep [Å] are system-level incident prescriptions - already in
   ! config%reactants (K1 §1.4);
   ! member parameters do not duplicate them; this member carries only the distribution shape
   ! and orientation prescription (mode switches + shape parameters),
   ! and when sample is finalized (C5) it takes e_rel/b_max from config%reactants per mode
   type :: inc_params_t
      integer :: n_e_rel  = 1      ! collision-energy prescription mode [-] (old NREL: 0 = thermal
                                   ! population, 1 = fixed value, 2 = density-weighted Gaussian)
      real(8) :: t_trans  = 0.0d0  ! translational temperature [K] (old TRANS - mode-0 thermal
                                   ! population: E_rel = R*GAMA(2)*T conversion, SURF.f:261-266)
      real(8) :: v_width  = 0.0d0  ! velocity width alpha (the width of the mode-2
                                   ! density-weighted DENS,
      ! in internal velocity units [Å/(10 fs)]; old ALPHAA)
      integer :: n_thta   = 0      ! polar-angle prescription mode [-] (fixed beam theta =
      ! theta_max vs
                                   ! P(theta) ~ sin theta sampling - SURF.f:133-141)
      real(8) :: thta_max = 0.0d0  ! polar-angle upper bound theta_max [rad] (old THTA - input
                                   ! degrees -> radians via consts dtor, grammar C4)
      integer :: n_chi    = 0      ! azimuth prescription mode [-] (fixed chi vs uniform on
      ! [0,2pi) - SURF.f:142-149)
      real(8) :: chi      = 0.0d0  ! azimuth chi [rad] (old CHI - degrees -> radians as above)
      integer :: n_b      = 0      ! aiming mode [-] (fixed b = b_max vs b = b_max*sqrt(U)
      ! disk-uniform - b^2 uniform on [0,b_max^2], manual
                                   ! Eq. V.19, SURF.f:171-187)
      integer :: n_aim    = 0      ! drop-point mode [-] (uniform within the surface unit cell
      ! vs random cluster box
      ! - SURF.f:153-169/:211-217; the box half-width belongs to
      ! roster%box_pad, unit-cell edges/skew to the roster)
      ! Left blank: Gaussian density-peak cache (old PVMAX style) and rejection-sampling cap
      ! (old 5*v1),
      ! beam/thermal research prescription-combination presets - finalized when the C5 members
      ! are concretized
   end type
contains
   !------------------------------------------------------------------
   ! inc_init(p) - initialize the member (validate the mode combination, precompute the
   ! orientation constants)
   ! Reference: src_VENUS/SURF.f:120-150 (orientation-angle prescription - see the module
   ! header)
   subroutine inc_init(p)
      type(inc_params_t), intent(in) :: p ! member parameters (mode switches / shape parameters)
      ! Pseudocode:
      ! 1. Validate the mode combination (illegal -> loud-stop): n_e_rel in {0,1,2} with mode 0
      !    requiring t_trans > 0 and mode 2 requiring v_width > 0; n_thta/n_chi/n_b/n_aim
      !    values legal; orientation prescriptions are meaningful only for surface/cluster
      !    configurations (the paradigm comes from
      !    config%reactants%surface_model - checked in the assembly phase)
      ! 2. Precompute orientation constants: sin/cos(theta_max) and sin/cos(chi); the fixed-beam
      !    incident-direction
      !    unit vector u = (-sin(theta)*cos(chi), -sin(theta)*sin(chi), -|cos(theta)|)
      !    (pointing at the surface,
      !    SURF.f:331-338 - velocity along -z guarantees incidence)
   end subroutine inc_init
   ! -- sample direction (finalized in C5) -- drop point (rx0,ry0) sampled per the n_aim
   !   prescription (uniform within the
   !   unit cell with skew correction / cluster box +/- half-width); aiming b = b_max*sqrt(U)
   !   (or fixed b_max), azimuth
   !   phi uniform on [0,2pi); E_rel assigned per mode (0: GAMA(2) thermal population, 1: the
   !   fixed value from
   !   config%reactants%e_rel, 2: DENS = v^3*exp(-((v-v0)/alpha)^2) rejection sampling);
   !   the relative velocity along u becomes momenta split by mass across the two fragments
   !   (surface/gas routes merged into one formula)
   ! sample routines: finalized when the C5 members are concretized (not declared at this stage)
   ! Registration-line note: samp_reg registers this member (name 'incident') in the sampler
   ! distribution member table
   ! - C5 wiring
end module dist_incident
