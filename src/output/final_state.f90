!=====================================================================
! final_state.f90 - final-state analysis quantities of one trajectory's end state
! Design:
!   Computes fragment internal energies / angular momenta, generalized COM relative
!   quantities, and the scattering-angle family, all driven by the atom list fragment
!   table (gas/surface unified, no system branches). Surface-specific quantities
!   (the surface mu and the ang(19)/ang(20) pair) follow the resolved projectile
!   role list_atoms%i_proj; the A/B channel bookkeeping stays positional. Angular
!   momenta in hbar, energies in kcal/mol, angles in degrees. Also hosts the
!   container's single outcome-classification slot (bind + statistics-phase call +
!   bound query).
!=====================================================================
module final_state
   use consts,        only: hbar_code, e_conv, dtor, two_pi
   use state,         only: state_t
   use config,        only: reactants
   use config_atoms, only: list_atoms
   use geometry,      only: cenmas, amom, angvel
   use force_interface,  only: force_eval
   implicit none
   private
   public :: fin_t, fin_eval
   public :: container_bind_classify, container_classify, container_classify_bound
                                                          ! the outcome-classification slot:
                                                          ! bind + statistics-phase call +
                                                          ! bound query (driver skip guard)

   ! Final-state quantity receiver type (fin_eval output)
   type :: fin_t
      real(8) :: e_int_a(3) = 0.0d0   ! fragment-A internal energy [kcal/mol]: (1) internal kinetic
                                      ! (2) internal potential - reference (3) total
      real(8) :: e_int_b(3) = 0.0d0   ! fragment-B internal energy [kcal/mol] (same layout; the atom list%frag(2) row)
      real(8) :: j_a(3)     = 0.0d0   ! fragment-A angular momentum vector [hbar]
      real(8) :: j_mag_a    = 0.0d0   ! fragment-A angular momentum magnitude [hbar]
      real(8) :: e_rot_a    = 0.0d0   ! fragment-A rotational energy [kcal/mol]
      real(8) :: j_b(3)     = 0.0d0   ! fragment-B angular momentum vector [hbar]
      real(8) :: j_mag_b    = 0.0d0   ! fragment-B angular momentum magnitude [hbar]
      real(8) :: e_rot_b    = 0.0d0   ! fragment-B rotational energy [kcal/mol]
      real(8) :: e_rel      = 0.0d0   ! relative translational energy (radial projection) [kcal/mol]
      real(8) :: e_rel_sq   = 0.0d0   ! relative translational energy (full velocity) [kcal/mol]
      real(8) :: oam(3)     = 0.0d0   ! orbital angular momentum vector [hbar]
      real(8) :: oam_mag    = 0.0d0   ! orbital angular momentum magnitude [hbar]
      real(8) :: b_fin      = 0.0d0   ! final impact parameter [Å] (= oam_mag/mu/|v_rel|)
      real(8) :: e_cm       = 0.0d0   ! whole-system center-of-mass translational energy [kcal/mol]
      real(8) :: ang(20)    = 0.0d0   ! scattering-angle family [deg] (1-3 final-value pairs,
                                      ! 4-16 initial-value-dependent combinations, 17-18 surface
                                      ! incident, 19-20 surface scattering)
      real(8) :: t_life     = 0.0d0   ! trajectory lifetime [fs] (= sta%t*10)
      real(8), allocatable :: q_fin(:) ! final COM-frame coordinate archive [Å] (for normal-mode analysis)
      real(8), allocatable :: p_fin(:) ! final COM-frame momentum archive [amu·Å/(10 fs)]
      real(8), allocatable :: q_raw(:) ! final raw system-frame coordinates [Å] (all atoms,
                                      ! as the trajectory left them - the channel-classification
                                      ! and any inter-fragment distance consumer reads THIS,
                                      ! not q_fin which is per-fragment COM-relative)
      ! Left blank: diatomic n/j inversion (finalized with the ebk member); initial-value
      ! reference quantities and rotational-energy sampling fluctuations (deferred: no
      ! producer exists yet; ang(4..18) all consume an initial value, so they stay
      ! unreachable until a sampler->statistics initial-value snapshot channel exists)
   end type fin_t

   ! The outcome-classification slot: an optional module-private procedure
   ! pointer holding the container's own classification judgment (problem
   ! knowledge). Contract: read-only fin in, the outcome name out; hosted
   ! here as the home of fin_t (hosting it in force_interface would close a
   ! use cycle). Left null and called = named abort (the legal
   ! no-classification run never calls)
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
   ! fin_eval(sta, fin) - compute all final-state quantities: fragment internal
   ! energies / angular momenta, generalized COM relative quantities, and the
   ! scattering-angle family (list_atoms fragment-table driven, no system branches)
   !------------------------------------------------------------------
   subroutine fin_eval(sta, fin)
      type(state_t), intent(in) :: sta    ! physical state (read-only; the rotation subtraction
      ! and the buffered-geometry evaluation ride local copies, never modifying sta)
      type(fin_t), intent(out) :: fin     ! final-state quantity receiver (the sole output)
      type(state_t) :: wrk, stg           ! working copy (rotation subtraction) / buffered state
      real(8) :: qcms(2, 3), vcms(2, 3)   ! per-fragment COM position / velocity
      real(8) :: q_cm(3), v_cm(3)         ! one fragment's COM (the cenmas returns)
      real(8) :: ptot(3)                  ! the whole-system linear momentum (per component)
      real(8) :: j_vec(3), e_rot, omg(3)  ! one fragment's angular momentum / energy / omega
      real(8) :: tkin, e_pot              ! internal kinetic / buffered potential [kcal/mol]
      real(8) :: qr(3), vr(3)             ! relative COM position / velocity
      real(8) :: qr_mag, vr_mag, v_rad    ! their magnitudes + the radial speed
      real(8) :: mu, m_a, m_b             ! reduced / fragment masses [amu]
      real(8) :: oam_v(3), oam_m          ! orbital angular momentum (physical / magnitude)
      real(8) :: cs, v_xy                 ! an angle cosine / the transverse speed
      integer :: nf, i, j, k, a, d        ! fragment / atom / component loop indices
      ! 1. guards: the fragment table must be assembled, agree with the state, and fit
      !    the two fragment channels of the receiver (a single fragment is the legal
      !    degenerate run: every pair quantity stays zero)
      if (.not. allocated(list_atoms%mass) .or. .not. allocated(list_atoms%frag)) then
         write (0, '(a)') 'fin_eval: the atom list was never assembled (no fragment table - '// &
                          'the final-state quantities need the fragment masses and lists)'
         write (0, '(a)') 'fin_eval: fatal (the final-state quantities are not computable)'
         stop 1
      end if
      if (size(sta%q) /= 3*size(list_atoms%mass)) then
         write (0, '(a,i0,a,i0)') 'fin_eval: state/list_atoms size mismatch: size(q) = ', size(sta%q), &
            ', 3*list_atoms atoms = ', 3*size(list_atoms%mass)
         write (0, '(a)') 'fin_eval: fatal (the state and the fragment table disagree)'
         stop 1
      end if
      nf = size(list_atoms%frag)
      if (nf < 1 .or. nf > 2) then
         write (0, '(a,i0)') 'fin_eval: illegal fragment count in the atom list: ', nf
         write (0, '(a)') 'fin_eval: fatal (the receiver carries two fragment channels)'
         stop 1
      end if
      if (reactants%surface_model /= 0 .and. (list_atoms%i_proj < 1 .or. list_atoms%i_proj > nf)) then
         write (0, '(a)') 'fin_eval: surface quantities need the resolved projectile role '// &
                          '(list_atoms%i_proj unset - the list_atoms_load role reconcile must run first)'
         write (0, '(a)') 'fin_eval: fatal (the final-state quantities are not computable)'
         stop 1
      end if
      ! 2. per-fragment properties + the COM-frame archive
      !    Each fragment in turn: build its COM frame, take the angular momentum /
      !    rotational energy / internal kinetic, strip the rigid-body rotation from
      !    the working momenta, then archive the fragment's rows relative to its own
      !    COM. A one-atom fragment carries no internal structure: every internal
      !    quantity stays zero and its archive rows stay zero.
      wrk = sta
      allocate (fin%q_fin(size(sta%q)), fin%p_fin(size(sta%q)))
      fin%q_fin = 0.0d0
      fin%p_fin = 0.0d0
      allocate (fin%q_raw(size(sta%q)))
      fin%q_raw = sta%q
      do i = 1, nf
         call cenmas(wrk, list_atoms%mass, list_atoms%frag(i)%list, q_cm, v_cm)
         qcms(i, :) = q_cm
         vcms(i, :) = v_cm
         if (list_atoms%frag(i)%nat < 2) cycle
         omg = 0.0d0                         ! a linear fragment computes no omega
         call amom(wrk, list_atoms%mass, list_atoms%frag(i)%list, j_vec, e_rot, omg)
         tkin = 0.0d0
         do k = 1, list_atoms%frag(i)%nat
            a = list_atoms%frag(i)%list(k)
            do d = 1, 3
               tkin = tkin + (wrk%p(3*(a - 1) + d) - list_atoms%mass(a)*v_cm(d))**2 &
                              /(2.0d0*list_atoms%mass(a))
            end do
         end do
         call angvel(wrk, list_atoms%mass, list_atoms%frag(i)%list, omg)
         do k = 1, list_atoms%frag(i)%nat
            a = list_atoms%frag(i)%list(k)
            do d = 1, 3
               fin%q_fin(3*(a - 1) + d) = wrk%q(3*(a - 1) + d) - q_cm(d)
               fin%p_fin(3*(a - 1) + d) = wrk%p(3*(a - 1) + d) - list_atoms%mass(a)*v_cm(d)
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
      ! 3. internal potentials: the buffered-geometry evaluation
      !    Fragment i's internal potential = the potential energy of the FULL
      !    configuration that keeps fragment i at its current geometry and places
      !    every other fragment at its tabulated equilibrium arrangement about its
      !    own COM (inter-fragment posture kept, the other fragments' internal
      !    distortion removed), taken on a zero-momentum local copy; sta is never
      !    touched. No reference energy is subtracted yet, and the current posture
      !    (not full separation) is kept, so e_int(2) carries the finite-separation
      !    residual inter-fragment potential.
      do i = 1, nf
         if (list_atoms%frag(i)%nat < 2) cycle
         stg = sta
         stg%p = 0.0d0
         do j = 1, nf
            if (j == i) cycle
            do k = 1, list_atoms%frag(j)%nat
               a = list_atoms%frag(j)%list(k)
               do d = 1, 3
                  stg%q(3*(a - 1) + d) = list_atoms%frag(j)%qz(3*(k - 1) + d) + qcms(j, d)
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
      ! 4. whole-system COM translational energy (computed for every fragment count,
      !    single-fragment runs included): the square of the per-component momentum
      !    SUMS - the COM drifts as one body of mass sum(m); the sum of squared
      !    momenta would re-count the per-fragment internal kinetic channels
      do d = 1, 3
         ptot(d) = sum(sta%p(d::3))
      end do
      fin%e_cm = dot_product(ptot, ptot)/(2.0d0*sum(list_atoms%mass))/e_conv
      ! 5. relative motion of the fragment pair (two fragments only)
      if (nf == 2) then
         m_a = list_atoms%frag(1)%mass
         m_b = list_atoms%frag(2)%mass
         if (reactants%surface_model == 0) then
            mu = m_a*m_b/(m_a + m_b)         ! the gas-phase pair form (role-symmetric)
         else
            mu = list_atoms%frag(list_atoms%i_proj)%mass   ! the surface form: the projectile fragment alone
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
      ! 6. final-value scattering angles (the family that needs no initial values; the
      !    initial-value channels of the angle family stay at their deferred zero)
      if (nf == 2 .and. list_atoms%frag(1)%nat /= 1 .and. list_atoms%frag(2)%nat /= 1) then
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
      ! 7. surface scattering pair (the polar/equatorial angles of the projectile
      !    fragment's COM velocity, following the resolved role row; the incident pair
      !    needs the initial velocity and stays deferred with it)
      if (reactants%surface_model /= 0) then
         vr_mag = sqrt(dot_product(vcms(list_atoms%i_proj, :), vcms(list_atoms%i_proj, :)))
         if (vr_mag >= 1.0d-12) then
            cs = vcms(list_atoms%i_proj, 3)/vr_mag
            fin%ang(19) = acos(max(-1.0d0, min(1.0d0, cs)))/dtor
         end if
         v_xy = sqrt(vcms(list_atoms%i_proj, 1)**2 + vcms(list_atoms%i_proj, 2)**2)
         if (v_xy > 1.0d-6) then
            cs = vcms(list_atoms%i_proj, 1)/v_xy
            if (vcms(list_atoms%i_proj, 2) >= 0.0d0) then
               fin%ang(20) = acos(max(-1.0d0, min(1.0d0, cs)))/dtor
            else
               fin%ang(20) = (two_pi - acos(max(-1.0d0, min(1.0d0, cs))))/dtor
            end if
         end if
      end if
      ! 8. trajectory lifetime (internal time unit 10 fs -> fs)
      fin%t_life = sta%t*10.0d0
   end subroutine fin_eval
   !------------------------------------------------------------------
   ! container_bind_classify(proc) - fill the classification slot: store
   !                 the container's outcome-naming judgment proc into the
   !                 module-private procedure pointer
   !                 container_classify_proc; the statistics phase then
   !                 calls it through container_classify. Optional - a run
   !                 without classification guards with
   !                 container_classify_bound and never calls. Executed
   !                 once, at assembly, by the container's reg module; a
   !                 second call aborts
   !------------------------------------------------------------------
   subroutine container_bind_classify(proc)
      procedure(classify_i) :: proc       ! the container classification judgment
      if (associated(container_classify_proc)) then
         write (0, '(a)') 'container_bind_classify: the classification slot is already '// &
                        'bound (container_bind_classify ran twice - an assembly-'// &
                        'choreography error, formal)'
         write (0, '(a)') 'container_bind_classify: fatal (the container is unusable)'
         stop 1
      end if
      container_classify_proc => proc      ! executed once in the assembly phase
   end subroutine container_bind_classify

   !------------------------------------------------------------------
   ! container_classify(fin) - statistics-phase entry: the outcome name of one
   ! trajectory's final state (unbound slot + a call aborts with a named error; the legal
   ! no-classification run guards with container_classify_bound and never calls)
   !------------------------------------------------------------------
   function container_classify(fin) result(name)
      type(fin_t), intent(in) :: fin      ! final-state quantities (read-only)
      character(len=32) :: name           ! the outcome name [-]
      if (.not. associated(container_classify_proc)) then
         write (0, '(a)') 'container_classify: the classification slot is not bound '// &
                        '(container_bind_classify never ran - the legal no-classification '// &
                        'run guards with container_classify_bound and never calls here)'
         write (0, '(a)') 'container_classify: fatal (the container is unusable)'
         stop 1
      end if
      name = container_classify_proc(fin)
   end function container_classify

   !------------------------------------------------------------------
   ! container_classify_bound() - has a classification judgment been bound?
   ! (read-only driver-side skip guard: call container_classify only when .true.)
   !------------------------------------------------------------------
   pure logical function container_classify_bound()
      container_classify_bound = associated(container_classify_proc)
   end function container_classify_bound
end module final_state
