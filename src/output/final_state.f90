!=====================================================================
! analysis/final_state.f90 - L4 analysis: final-state quantities (generalized COM relative
!                            energies / angular momenta /
!                            scattering-angle family - roster fragment-table data-driven,
!                            gas/surface unified)
! Reference: src_VENUS/FINAL.f:100-109 (fragment-A angular momentum / rotational energy -
!       CENMAS+ROTN,
!       folded to hbar units by /C7), :113-129 (fragment internal kinetic energy sum p^2/2m/C1
!       and internal potential V-DELH -
!       DVDQ_1+ENERGY_1 at the fragment geometry), :140-144/:295-301 (diatomic n/j inversion
!       FINLNJ), :149-174 (angular-velocity subtraction angvel + QFINAL/PFINAL final-state
!       archive),
!       :356-372 (whole-system center-of-mass translational energy ETCM), :376-385 (two
!       reduced-mass formulas: gas-phase
!       WTA*WTB/(WTA+WTB), surface WTA), :387-410 (relative motion: center-of-mass difference
!       / velocity difference ->
!       EREL/ERELSQ/OAM/BF), :432-531 (gas-phase scattering angles ANG(1..16) - VI/VF,
!       LI/LF, JAI/JBI pairwise dot products then acos), :538-586 (surface angles ANG(17..18) -
!       - the surface normal from
!       the NN1-NN3 three-atom cross product), :598-637 (surface incident/scattering angles
!       ANG(17..20));
!       src_VENUS/GFINAL.f:64-75 (lifetime TLIFE and energy-conservation drift printout - the
!       drift belongs to the
!       recorder energy columns), :199-202 (AErelS/AThtaS final-state records - to the recorder,
!       boundary note)
! Excision warning (the hardwired EREL overwrite - K2 plan sourcing note):
!   GFINAL.f:78-81 overwrites EREL with the C-atom kinetic energy NORM2(P(7:9))^2/2W(3)/C1 (1D
!   model
!   system knowledge) - not carried over: e_rel always uses the FINAL.f:387-401 generalized
!   COM relative form;
!   Same source, also not carried over: the FINAL.f:616-625 NZDOWN==1 branch hardwires the B
!   velocity as P(7:9)/W(3) -
!   the surface final-state velocity channel is specified via roster fragments
!=====================================================================
module final_state
   use consts,        only: hbar_code, e_conv, dtor
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use geometry,      only: cenmas, amom, angvel
   use force_facade,  only: force_eval
   implicit none
   private
   public :: fin_t, fin_eval

   ! Final-state quantity receiver type (fin_eval output; fields = the old global final-state
   ! quantities turned into members)
   type :: fin_t
      real(8) :: e_int_a(3) = 0.0d0   ! fragment-A internal energy [kcal/mol]: (1) internal kinetic (2) internal potential - reference (3) total (old EA(1..3))
      real(8) :: e_int_b(3) = 0.0d0   ! fragment-B internal energy [kcal/mol] (old EB(1..3); B table deferred with the roster - phase C)
      real(8) :: j_a(3)     = 0.0d0   ! fragment-A angular momentum vector [hbar] (old AMA(1:3), /C7 conversion)
      real(8) :: j_mag_a    = 0.0d0   ! fragment-A angular momentum magnitude [hbar] (old AMA(4))
      real(8) :: e_rot_a    = 0.0d0   ! fragment-A rotational energy [kcal/mol] (old EROTA)
      real(8) :: j_b(3)     = 0.0d0   ! fragment-B angular momentum vector [hbar] (old AMB(1:3))
      real(8) :: j_mag_b    = 0.0d0   ! fragment-B angular momentum magnitude [hbar] (old AMB(4))
      real(8) :: e_rot_b    = 0.0d0   ! fragment-B rotational energy [kcal/mol] (old EROTB)
      real(8) :: e_rel      = 0.0d0   ! relative translational energy (radial projection) [kcal/mol] (old EREL - FINAL.f:394-401 generalized form)
      real(8) :: e_rel_sq   = 0.0d0   ! relative translational energy (full velocity) [kcal/mol] (old ERELSQ)
      real(8) :: oam(3)     = 0.0d0   ! orbital angular momentum vector [hbar] (old OAM(1:3), /C7 conversion)
      real(8) :: oam_mag    = 0.0d0   ! orbital angular momentum magnitude [hbar] (old OAM(4))
      real(8) :: b_fin      = 0.0d0   ! final impact parameter [Å] (old BF = OAM(4)/mu/|v_rel|)
      real(8) :: e_cm       = 0.0d0   ! whole-system center-of-mass translational energy [kcal/mol] (old ETCM)
      real(8) :: ang(20)    = 0.0d0   ! scattering-angle family [deg] (old ANG(1..20): 1-3 LF.JAF/JBF.JAF/JAF.JBF, 4-7 VI.VF/LI.LF/LI.JAF/LI.JBF, 8-16 JAI/JBI combinations, 17-18 surface incident theta/chi, 19-20 surface scattering theta/chi)
      real(8) :: t_life     = 0.0d0   ! trajectory lifetime [fs] (old TLIFE = TIME*NC*10 - GFINAL.f:64; new version sta%t*10)
      real(8), allocatable :: q_fin(:) ! final COM-frame coordinate archive [Å] (old QFINAL - for normal-mode analysis)
      real(8), allocatable :: p_fin(:) ! final COM-frame momentum archive [amu·Å/(10 fs)] (old PFINAL)
      ! Left blank: diatomic n/j inversion (old AN/AJ/BN/BJ - FINLNJ, finalized with the ebk
      ! member);
      ! initial-value reference quantities (old VI/OAMI/AMAI/AMBI - snapshotted when the sampler
      ! assembles, the hand-off channel
      ! finalized in phase C); rotational-energy sampling fluctuations SDA/SDB (old NCROT
      ! resampling average -
      ! to the statistics phase); fragment-B table (same-source roster blank)
   end type fin_t
contains
   !------------------------------------------------------------------
   ! fin_eval(sta, fin) - final-state quantity computation: fragment internal energies /
   ! angular momenta + generalized COM relative quantities +
   ! scattering-angle family (roster fragment-table driven, no system branches)
   ! Reference: FINAL.f:100-174, :356-410, :432-531, :538-637 (sections per the module header)
   subroutine fin_eval(sta, fin)
      type(state_t), intent(in) :: sta    ! physical state (final-state snapshot - read-only;
                                          ! the fragment-geometry staging trick is in step 2,
      ! never modifying sta itself)
      type(fin_t), intent(out) :: fin     ! final-state quantity receiver (the sole output of this routine)
      ! Pseudocode:
      ! 1. Fragment-A properties (old :100-109): cenmas(sta, roster%mass, roster%list_a,
      !    q_cm_a, v_cm_a) builds the COM frame; amom(..., j_a, e_rot_a, omega_a) - j_a and
      !    j_mag_a folded to hbar units (divide hbar_code, old /C7); one-atom fragments all zero (old :24-33)
      ! 2. Fragment-A internal energy (old :113-129): internal kinetic = COM-frame
      !    sum p^2/(2*mass)/e_conv;
      !    internal potential = v from force_eval at the fragment's equilibrium-arrangement
      !    coordinates (the old :49-62 QZB+QCMB
      !    staging-evaluation trick - the new version evaluates via the facade on a locally
      !    staged state carrying the trajectory-final sta%occ - the adiab-mode delivery
      !    request source of force_eval, spec §4 (the old scalar surface request generalized
      !    to the occupation vector; length 1 = the degenerate form), never writing back to
      !    sta) minus the reference
      !    dissociation energy (old V-DELH; the DELH energy reference belongs to the assembly
      !    phase - finalized in phase C);
      !    e_int_a(3) = (1)+(2)
      ! 3. Fragment-B same form (old :264-301; B table deferred with the roster - phase C)
      ! 4. e_cm (old :356-372): summed over all atoms, sum of the three squared p components /
      !    sum(mass) / 2 / e_conv
      ! 5. Relative motion (old :376-410): reduced mass mu per the roster fragment table -
      !    two-fragment gas-phase form
      !    mu = m_a*m_b/(m_a+m_b), surface form (single active fragment) mu = m_a (old :376-385
      !    two branches unified as fragment-table data-driven, reactants%surface_model only picks
      !    the channel);
      !    q_r = q_cm_a - q_cm_b, v_r likewise; e_rel = mu*(v_r.q̂_r)^2/2/e_conv (radial
      !    projection, old :394-401); e_rel_sq = mu*|v_r|^2/2/e_conv; oam = mu*(q_r x v_r)
      !    folded to hbar (old :403-410); b_fin = |oam|/mu/|v_r| (old :407)
      ! 6. Gas-phase scattering angles ang(1..16) (old :432-531): the initial/final relative
      ! velocity vectors (VI and v_r)
      !    against the orbital/rotational angular momentum vectors (OAMI/OAM, AMAI/AMBI and
      !    j_a/j_b), pairwise dot products
      !    then acos, divided by dtor to degrees (old /C4); angular momenta with magnitude < 1e-5
      !    skip the term (old :460/:482
      !    thresholds); dot products clamped to [-1,1] (old :435-436); the initial values
      !    VI/OAMI/AMAI/AMBI
      !    come from the initial-state snapshot (deferred channel)
      ! 7. Surface angles ang(17..20) (old :598-637): incident angles from the initial velocity
      !    VI (theta/chi two forms :600-612),
      !    scattering angles from the active-fragment COM velocity v_cm_a (:614-637; the NZDOWN
      !    hardwired branch not carried over - excision warning in the module header); polar
      !    degenerate positions (|v_xy| <= 1e-6 ->
      !    chi=0, old :604-611/:629-636); surface-normal reference (old :543-555 from the
      !    reference three-atom
      !    cross product - as data via roster/unit-cell fields, finalized in phase C)
      ! 8. Final-state archive (old :149-174): COM-frame coordinates/momenta after angvel
      !    removes rigid-body rotation
      !    -> q_fin/p_fin (for normal-mode analysis - consumed by density/members)
      ! 9. t_life = sta%t*10 [fs] (internal time 10 fs -> fs, the old TLIFE convention)
      ! 10. Boundary: GFINAL's recording duties (AErelS/AThtaS :199-202) belong to the recorder;
      !     the hardwired EREL overwrite (GFINAL.f:78-81) is not carried over - e_rel keeps the
      !     step-5 generalized form
   end subroutine fin_eval
end module final_state
