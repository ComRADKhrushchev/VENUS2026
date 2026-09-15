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
!
! -- E19 outcome-classification ownership (review-gate ruling 2026-09-13: the outcome
!    expression family of output/outcome.f90 withdrawn - the L4-review symmetry ruling
!    promised by E18) --
!   Same family, same verdict as E18: classification criteria (what counts as reactive /
!   scattering / adsorbed / long-lived) are problem knowledge, equally system-specific -
!   the container registers ONE classification function, container_bind_classify(proc),
!   private contract classify_i(fin) -> name, called by the driver statistics phase
!   through container_classify(fin). outcome.f90 is DELETED (the expression trees / the
!   o_* leaf keys / the combinators / the outcome_tbl grammar key all withdrawn with it);
!   the generic per-name counting stays in the driver. THIS module hosts the slot - the
!   home of the contract's argument type fin_t (hosting it in force_facade would close a
!   cycle: final_state use force_facade (fin_eval's staged evaluation) -> force_facade
!   use final_state). Unbound slot + a call = loud-stop (classification is a
!   statistics-phase contract - unlike termination there is no safety-net-only form; not
!   calling it is the legal no-classification run).
!=====================================================================
module final_state
   use consts,        only: hbar_code, e_conv, dtor, two_pi
   use state,         only: state_t
   use config,        only: reactants
   use config_roster, only: roster
   use geometry,      only: cenmas, amom, angvel
   use force_facade,  only: force_eval
   implicit none
   private
   public :: fin_t, fin_eval
   public :: container_bind_classify, container_classify, container_classify_bound
                                                          ! E19: the outcome-classification
                                                          ! slot + its statistics-phase entry
                                                          ! + the bound query (the driver-side
                                                          ! skip guard)

   ! Final-state quantity receiver type (fin_eval output; fields = the old global final-state
   ! quantities turned into members)
   type :: fin_t
      real(8) :: e_int_a(3) = 0.0d0   ! fragment-A internal energy [kcal/mol]: (1) internal kinetic (2) internal potential - reference (3) total (old EA(1..3))
      real(8) :: e_int_b(3) = 0.0d0   ! fragment-B internal energy [kcal/mol] (old EB(1..3); the roster%frag(2) row, E15)
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
      ! initial-value reference quantities (old VI/OAMI/AMAI/AMBI - deferred: no producer
      ! exists in this tree, and with them 15 of the 20 angle-family fields stay unreachable
      ! (ang(4..18) all consume an initial value - the final-value angles ang(1..3)/ang(19..20)
      ! alone are computed); the trigger is a sampler->statistics initial-value snapshot
      ! channel); rotational-energy sampling fluctuations SDA/SDB (old NCROT resampling
      ! average - deferred with the same channel: the resampling averages have no producer);
      ! the *_b fields' content: the roster%frag(2) row (E15 arrayed)
   end type fin_t

   ! The outcome-classification slot - E19: the container's own classification judgment
   ! (problem knowledge; the outcome expression family withdrawn). Contract: read-only fin,
   ! the outcome name out. Hosted HERE - the home of the argument type fin_t (see the E19
   ! header note for the cycle argument). Unbound slot + a call = loud-stop (the legal
   ! no-classification run simply never calls the entry)
   abstract interface
      function classify_i(fin) result(name)
         import :: fin_t
         type(fin_t), intent(in) :: fin    ! final-state quantities (read-only)
         character(len=32) :: name         ! the outcome name [-] (e.g. 'nonreact'/'react'/
                                           ! 'long' - data chosen by the container author)
      end function classify_i
   end interface
   procedure(classify_i), pointer :: container_classify_proc => null()
contains
   !------------------------------------------------------------------
   ! fin_eval(sta, fin) - final-state quantity computation: fragment internal energies /
   ! angular momenta + generalized COM relative quantities +
   ! scattering-angle family (roster fragment-table driven, no system branches)
   ! Reference: FINAL.f:100-174, :356-410, :432-531, :538-637 (sections per the module header)
   subroutine fin_eval(sta, fin)
      type(state_t), intent(in) :: sta    ! physical state (final-state snapshot - read-only;
                                          ! the rotation subtraction and the staged-geometry
      ! evaluation ride local copies, never modifying sta itself)
      type(fin_t), intent(out) :: fin     ! final-state quantity receiver (the sole output of this routine)
      type(state_t) :: wrk, stg           ! working copy (rotation subtraction) / staged state
      real(8) :: qcms(2, 3), vcms(2, 3)   ! per-fragment COM position / velocity
      real(8) :: q_cm(3), v_cm(3)         ! one fragment's COM (the cenmas returns)
      real(8) :: ptot(3)                  ! the whole-system linear momentum (per component)
      real(8) :: j_vec(3), e_rot, omg(3)  ! one fragment's angular momentum / energy / omega
      real(8) :: tkin, e_pot              ! internal kinetic / staged potential [kcal/mol]
      real(8) :: qr(3), vr(3)             ! relative COM position / velocity
      real(8) :: qr_mag, vr_mag, v_rad    ! their magnitudes + the radial speed
      real(8) :: mu, m_a, m_b             ! reduced / fragment masses [amu]
      real(8) :: oam_v(3), oam_m          ! orbital angular momentum (physical / magnitude)
      real(8) :: cs, v_xy                 ! an angle cosine / the transverse speed
      integer :: nf, i, j, k, a, d        ! fragment / atom / component loop indices
      ! ---- guards: the fragment table must be assembled, agree with the state, and fit
      ! the two fragment channels of the receiver (a single fragment is the legal
      ! degenerate run: every pair quantity stays zero)
      if (.not. allocated(roster%mass) .or. .not. allocated(roster%frag)) then
         write (0, '(a)') 'fin_eval: the roster was never assembled (no fragment table - '// &
                          'the final-state quantities need the fragment masses and lists)'
         write (0, '(a)') 'fin_eval: loud-stop (the final-state quantities are not computable)'
         stop 1
      end if
      if (size(sta%q) /= 3*size(roster%mass)) then
         write (0, '(a,i0,a,i0)') 'fin_eval: state/roster size mismatch: size(q) = ', size(sta%q), &
            ', 3*roster atoms = ', 3*size(roster%mass)
         write (0, '(a)') 'fin_eval: loud-stop (the state and the fragment table disagree)'
         stop 1
      end if
      nf = size(roster%frag)
      if (nf < 1 .or. nf > 2) then
         write (0, '(a,i0)') 'fin_eval: illegal fragment count in the roster: ', nf
         write (0, '(a)') 'fin_eval: loud-stop (the receiver carries two fragment channels)'
         stop 1
      end if
      ! ---- per-fragment properties + the COM-frame archive ----
      ! Each fragment in turn: build its COM frame, take the angular momentum /
      ! rotational energy / internal kinetic, strip the rigid-body rotation from the
      ! working momenta, then archive the fragment's rows relative to its own COM.
      ! A one-atom fragment carries no internal structure: every internal quantity
      ! stays zero and its archive rows stay zero (the lone atom sits at its own COM).
      wrk = sta
      allocate (fin%q_fin(size(sta%q)), fin%p_fin(size(sta%q)))
      fin%q_fin = 0.0d0
      fin%p_fin = 0.0d0
      do i = 1, nf
         call cenmas(wrk, roster%mass, roster%frag(i)%list, q_cm, v_cm)
         qcms(i, :) = q_cm
         vcms(i, :) = v_cm
         if (roster%frag(i)%nat < 2) cycle
         omg = 0.0d0                         ! a linear fragment computes no omega
         call amom(wrk, roster%mass, roster%frag(i)%list, j_vec, e_rot, omg)
         tkin = 0.0d0
         do k = 1, roster%frag(i)%nat
            a = roster%frag(i)%list(k)
            do d = 1, 3
               tkin = tkin + (wrk%p(3*(a - 1) + d) - roster%mass(a)*v_cm(d))**2 &
                              /(2.0d0*roster%mass(a))
            end do
         end do
         call angvel(wrk, roster%mass, roster%frag(i)%list, omg)
         do k = 1, roster%frag(i)%nat
            a = roster%frag(i)%list(k)
            do d = 1, 3
               fin%q_fin(3*(a - 1) + d) = wrk%q(3*(a - 1) + d) - q_cm(d)
               fin%p_fin(3*(a - 1) + d) = wrk%p(3*(a - 1) + d) - roster%mass(a)*v_cm(d)
            end do
         end do
         if (i == 1) then
            fin%e_int_a(1) = tkin/e_conv
            fin%j_a = j_vec/hbar_code
            fin%j_mag_a = sqrt(dot_product(j_vec, j_vec))/hbar_code
            fin%e_rot_a = e_rot
         else
            fin%e_int_b(1) = tkin/e_conv
            fin%j_b = j_vec/hbar_code
            fin%j_mag_b = sqrt(dot_product(j_vec, j_vec))/hbar_code
            fin%e_rot_b = e_rot
         end if
      end do
      ! ---- internal potentials: the staged-geometry evaluation ----
      ! The internal potential of fragment i is the potential energy of the FULL
      ! configuration that keeps fragment i at its current geometry and places every
      ! OTHER fragment at its tabulated equilibrium arrangement about that other
      ! fragment's own center of mass - the inter-fragment posture is left as it is,
      ! the other fragments' internal distortion is removed. The evaluation rides a
      ! zero-momentum local copy (the total energy then carries the potential alone),
      ! and the physical state is never touched. No reference energy is subtracted
      ! yet (the dissociation reference arrives with the assembly seam). Same family,
      ! registered with it - the staging-posture protocol: the old gas-phase branch
      ! staged the fragments at +ZASYM full separation before taking V, this
      ! implementation keeps the current inter-fragment posture, so e_int(2) carries
      ! the finite-separation residual inter-fragment potential (and the surface
      ! branch holds a reference-translation difference); the internal-energy
      ! reference protocol lands with the assembly seam.
      do i = 1, nf
         if (roster%frag(i)%nat < 2) cycle
         stg = sta
         stg%p = 0.0d0
         do j = 1, nf
            if (j == i) cycle
            do k = 1, roster%frag(j)%nat
               a = roster%frag(j)%list(k)
               do d = 1, 3
                  stg%q(3*(a - 1) + d) = roster%frag(j)%qz(3*(k - 1) + d) + qcms(j, d)
               end do
            end do
         end do
         call force_eval(stg, e_pot)
         if (i == 1) then
            fin%e_int_a(2) = e_pot
         else
            fin%e_int_b(2) = e_pot
         end if
      end do
      fin%e_int_a(3) = fin%e_int_a(1) + fin%e_int_a(2)
      fin%e_int_b(3) = fin%e_int_b(1) + fin%e_int_b(2)
      ! ---- whole-system COM translational energy (a whole-state quantity: computed
      ! for every fragment count, single-fragment runs included): the square of the
      ! per-component momentum SUMS - the COM drifts as one body of mass sum(m);
      ! the sum of squared momenta would instead re-count the per-fragment internal
      ! kinetic channels
      do d = 1, 3
         ptot(d) = sum(sta%p(d::3))
      end do
      fin%e_cm = dot_product(ptot, ptot)/(2.0d0*sum(roster%mass))/e_conv
      ! ---- relative motion of the fragment pair (two fragments only)
      if (nf == 2) then
         m_a = roster%frag(1)%mass
         m_b = roster%frag(2)%mass
         if (reactants%surface_model == 0) then
            mu = m_a*m_b/(m_a + m_b)         ! the gas-phase pair form
         else
            mu = m_a                         ! the surface form: the active fragment alone
         end if
         qr = qcms(1, :) - qcms(2, :)
         vr = vcms(1, :) - vcms(2, :)
         qr_mag = sqrt(dot_product(qr, qr))
         vr_mag = sqrt(dot_product(vr, vr))
         if (qr_mag >= 1.0d-12) then
            v_rad = dot_product(vr, qr)/qr_mag
         else
            v_rad = 0.0d0                    ! coincident centers: no radial direction
         end if
         fin%e_rel = mu*v_rad*v_rad/2.0d0/e_conv
         fin%e_rel_sq = mu*dot_product(vr, vr)/2.0d0/e_conv
         oam_v(1) = mu*(qr(2)*vr(3) - qr(3)*vr(2))
         oam_v(2) = mu*(qr(3)*vr(1) - qr(1)*vr(3))
         oam_v(3) = mu*(qr(1)*vr(2) - qr(2)*vr(1))
         oam_m = sqrt(dot_product(oam_v, oam_v))
         fin%oam = oam_v/hbar_code
         fin%oam_mag = oam_m/hbar_code
         if (vr_mag >= 1.0d-12) fin%b_fin = oam_m/(mu*vr_mag)
      end if
      ! ---- final-value scattering angles (the family that needs no initial values;
      ! the initial-value channels of the angle family stay at their deferred zero)
      if (nf == 2 .and. roster%frag(1)%nat /= 1 .and. roster%frag(2)%nat /= 1) then
         if (fin%oam_mag >= 1.0d-5 .and. fin%j_mag_a >= 1.0d-5) then
            cs = dot_product(fin%oam, fin%j_a)/(fin%oam_mag*fin%j_mag_a)
            fin%ang(1) = acos(max(-1.0d0, min(1.0d0, cs)))/dtor
         end if
         if (fin%oam_mag >= 1.0d-5 .and. fin%j_mag_b >= 1.0d-5) then
            cs = dot_product(fin%oam, fin%j_b)/(fin%oam_mag*fin%j_mag_b)
            fin%ang(2) = acos(max(-1.0d0, min(1.0d0, cs)))/dtor
         end if
         if (fin%j_mag_a >= 1.0d-5 .and. fin%j_mag_b >= 1.0d-5) then
            cs = dot_product(fin%j_a, fin%j_b)/(fin%j_mag_a*fin%j_mag_b)
            fin%ang(3) = acos(max(-1.0d0, min(1.0d0, cs)))/dtor
         end if
      end if
      ! ---- surface scattering pair (the polar/equatorial angles of the active
      ! fragment's COM velocity; the incident pair needs the initial velocity and
      ! stays deferred with it)
      if (reactants%surface_model /= 0) then
         vr_mag = sqrt(dot_product(vcms(1, :), vcms(1, :)))
         if (vr_mag >= 1.0d-12) then
            cs = vcms(1, 3)/vr_mag
            fin%ang(19) = acos(max(-1.0d0, min(1.0d0, cs)))/dtor
         end if
         v_xy = sqrt(vcms(1, 1)**2 + vcms(1, 2)**2)
         if (v_xy > 1.0d-6) then
            cs = vcms(1, 1)/v_xy
            if (vcms(1, 2) >= 0.0d0) then
               fin%ang(20) = acos(max(-1.0d0, min(1.0d0, cs)))/dtor
            else
               fin%ang(20) = (two_pi - acos(max(-1.0d0, min(1.0d0, cs))))/dtor
            end if
         end if
      end if
      ! ---- trajectory lifetime (internal time unit 10 fs -> fs)
      fin%t_life = sta%t*10.0d0
   end subroutine fin_eval
   !------------------------------------------------------------------
   ! container_bind_classify(proc) - E19 single-slot assembly injection: bind the
   !                 container's outcome-classification judgment (the caller: the
   !                 container's reg module, alongside container_bind_pes / bind_term)
   ! Reference: this file's E19 header note; E17(4)/E18 (problem knowledge, same family)
   subroutine container_bind_classify(proc)
      procedure(classify_i) :: proc       ! the container classification judgment
      if (associated(container_classify_proc)) then
         write (0, '(a)') 'container_bind_classify: the classification slot is already '// &
                        'bound (container_bind_classify ran twice - an assembly-'// &
                        'choreography error, formal)'
         write (0, '(a)') 'container_bind_classify: loud-stop (the container is unusable)'
         stop 1
      end if
      container_classify_proc => proc      ! executed once in the assembly phase by the
                                           ! container's reg module
   end subroutine container_bind_classify

   !------------------------------------------------------------------
   ! container_classify(fin) - statistics-phase entry (E19): the outcome name of one
   !                 trajectory's final state (the driver statistics step 1b; supersedes
   !                 outc_classify(fin, tbl))
   ! Reference: this file's E19 header note; VENUS.f90:252-284 (the old outcome summary
   !       printout - the generic per-name counting stays in the driver)
   function container_classify(fin) result(name)
      type(fin_t), intent(in) :: fin      ! final-state quantities (read-only)
      character(len=32) :: name           ! the outcome name [-]
      ! unbound slot + a call = loud-stop (a classification call without a bound
      ! judgment is an assembly error; the legal no-classification run guards with
      ! container_classify_bound() and simply never calls this entry)
      if (.not. associated(container_classify_proc)) then
         write (0, '(a)') 'container_classify: the classification slot is not bound '// &
                        '(container_bind_classify never ran - the legal no-classification '// &
                        'run guards with container_classify_bound and never calls here)'
         write (0, '(a)') 'container_classify: loud-stop (the container is unusable)'
         stop 1
      end if
      name = container_classify_proc(fin)
   end function container_classify

   !------------------------------------------------------------------
   ! container_classify_bound() - E19 bound query (public, read-only): has a
   !                 classification judgment been bound? The driver-side skip
   !                 guard for the legal no-classification run - the statistics
   !                 phase calls container_classify only when this is .true.
   pure logical function container_classify_bound()
      container_classify_bound = associated(container_classify_proc)
   end function container_classify_bound
end module final_state
