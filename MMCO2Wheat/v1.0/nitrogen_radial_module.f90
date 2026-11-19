! ======================================================================== !
! NITROGEN_RADIAL_MODULE  (M&M-CO2Wheat v1.0)
!
! Driver record 189 n_radial_uptake. Nitrate and ammonium uptake at the ROOT
! SURFACE, per soil layer, around one representative root (Barber-Cushman):
!   b dC/dt = (1/r) d/dr [ r (D_l theta f_l) dC/dr + r v0 r0/r C ]   (radial
!   diffusion + mass flow toward the root), with transporter uptake on the root
!   and root-hair surfaces J = Imax C_s/(Km + C_s) at the LOCAL solution
!   concentration C_s. Implicit finite volume on a log grid, hourly steps,
!   Picard iteration on the kinetics, with N's own transport constants.
!
! The alternative it replaces when on (record 189 = 0): the layer supply
!   MIN(pool, active + mass flow + diffusion), where the three are ADDED, so ions
!   that diffuse toward the root were counted as absorbed without passing a
!   transporter (about 70 % of the uptake under record 187), and the active
!   Michaelis-Menten term used the BULK concentration.
!
! Constants (no new fitted value):
!  - D_l: Li & Gregory (1974) Geochim Cosmochim Acta 38:703, Table 1 p.705,
!    25 C: NO3- 19.0e-6, NH4+ 19.8e-6 cm2 s-1; temperature by Stokes-Einstein
!    with water viscosity.
!  - impedance f_l = 1.1 (theta - theta_th), Olesen et al. 2001 Soil Sci
!    166:89, p.94 Eq.(7a/b), theta_th by p.93 Eq.(5) from clay, silt, rho_b
!    (22 soils). Barraclough & Tinker 1981 (1.58 theta - 0.17, two loams) was
!    rejected: it gives f_l = 0 at the sandy WDF soil's field water content.
!  - buffer power: NO3- is not adsorbed, b = theta; NH4+ b = theta + rho_b Kd,
!    with the deck's measured bulk density and Kd_NH4 (driver_N.csv).
!  - root radius and hairs: wheat, Gahoonia et al. (1997) Plant Soil 191:181,
!    Table 2.
!  - Imax: the model's existing per-gram-root Vmax (Goyal & Huffaker 1986)
!    times the existing f_Nstatus(NNI) and a_root, spread over the root PLUS
!    hair surface per gram, so the measured per-gram capacity is conserved.
!  - Km: the existing Km_c_no3 and Km_nh4 (Goyal & Huffaker 1986).
!  - specific root length: the existing nud_SRL (Chen et al. 2020).
!
! Coupling to the layer pools: the radial profile keeps its shape between days
! and is rescaled each day so that its content equals the layer's ion pool
! (leaching, mineralization, nitrification and fertilizer act on the pool).
! The solve returns the day's POTENTIAL uptake; when the caller takes only a
! fraction of it (demand gate or stock cap), nrad_commit blends the profile
! back linearly, which conserves mass exactly because b is constant in C.
! ======================================================================== !
MODULE NITROGEN_RADIAL_MODULE
  USE TYPES, ONLY: wp
  USE SHARED, ONLY: nud_SRL, nud_csif_slope, nud_theta_th
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nrad_layer_day, nrad_commit, nrad_lv_on, r0   ! r0 also read by rec 81 mode 4 (main)

  REAL(wp), PARAMETER :: PI = 3.14159265358979_wp
  REAL(wp), PARAMETER :: MN = 14.0067_wp            ! g N mol-1
  INTEGER,  PARAMETER :: NR = 25, NLM = 40, NION = 2 ! radial cells, layers, ions (1 NO3, 2 NH4)
  REAL(wp), PARAMETER :: D_NO3_25 = 1.90e-9_wp      ! m2 s-1, Li & Gregory 1974 Table 1
  REAL(wp), PARAMETER :: D_NH4_25 = 1.98e-9_wp      ! m2 s-1, Li & Gregory 1974 Table 1
  REAL(wp), PARAMETER :: r0  = 8.0e-5_wp            ! m, root radius, Gahoonia 1997 Table 2
  REAL(wp), PARAMETER :: r_h = 6.0e-6_wp            ! m, hair radius, Gahoonia 1997 p.181
  REAL(wp), PARAMETER :: l_h = 8.0e-4_wp            ! m, hair length, Gahoonia 1997 Table 2
  REAL(wp), PARAMETER :: n_h = 3.0e4_wp             ! hairs m-1 root, Gahoonia 1997 Table 2
  REAL(wp), PARAMETER :: nrad_lv_on = 10.0_wp       ! m m-3, a layer joins once its roots exceed this (numerical; r1 = 18 cm)

  LOGICAL  :: act(NLM) = .FALSE.
  REAL(wp) :: Cr(NR, NION, NLM) = 0.0_wp            ! mol m-3 solution, profile after today's potential solve
  REAL(wp) :: C0(NR, NION, NLM) = 0.0_wp            ! profile before the solve (for the commit blend)
  REAL(wp) :: frL(0:NR, NLM) = 0.0_wp, rcL(NR, NLM) = 0.0_wp, LvL(NLM) = 0.0_wp

CONTAINS

  SUBROUTINE build_grid(Lv, fr, rc)
    REAL(wp), INTENT(IN)  :: Lv
    REAL(wp), INTENT(OUT) :: fr(0:NR), rc(NR)
    REAL(wp) :: r1, q
    INTEGER  :: i
    r1 = MAX(1.0_wp / SQRT(PI * MAX(Lv, 1.0_wp)), 2.0_wp * r0)
    q  = (r1 / r0)**(1.0_wp / REAL(NR, wp))
    fr(0) = r0
    DO i = 1, NR
      fr(i) = r0 * q**i
      rc(i) = SQRT(fr(i-1) * fr(i))
    END DO
    fr(NR) = r1
  END SUBROUTINE build_grid

  ! Content of one ion in a layer's domain, g N m-2 ground.
  REAL(wp) FUNCTION content(C, fr, b, Lv, dz)
    REAL(wp), INTENT(IN) :: C(NR), fr(0:NR), b, Lv, dz
    INTEGER :: i
    content = 0.0_wp
    DO i = 1, NR
      content = content + PI * (fr(i)**2 - fr(i-1)**2) * b * C(i)
    END DO
    content = content * Lv * dz * MN
  END FUNCTION content

  ! ---------------------------------------------------------------------- !
  ! One layer, one day. Inputs: k layer index; th volumetric water; dz layer
  ! thickness (m); T_soil (C); Lv root length density (m m-3); water_mm root
  ! water extraction of the layer (mm d-1); pool_no3/pool_nh4 the layer's ion
  ! pools (g N m-2); Imax_g_no3/Imax_g_nh4 the transporter capacity per gram of
  ! root DM today (g N g-1 d-1, already times f_Nstatus and a_root); Km_no3,
  ! Km_nh4 (mol m-3); b_nh4 the ammonium buffer power. Outputs: U_no3, U_nh4 the
  ! potential uptake of the day (g N m-2). A layer without roots returns 0.
  ! ---------------------------------------------------------------------- !
  SUBROUTINE nrad_layer_day(k, th, dz, T_soil, Lv, water_mm, pool_no3, pool_nh4, &
                            Imax_g_no3, Imax_g_nh4, Km_no3, Km_nh4, b_nh4, U_no3, U_nh4)
    INTEGER,  INTENT(IN)  :: k
    REAL(wp), INTENT(IN)  :: th, dz, T_soil, Lv, water_mm, pool_no3, pool_nh4
    REAL(wp), INTENT(IN)  :: Imax_g_no3, Imax_g_nh4, Km_no3, Km_nh4, b_nh4
    REAL(wp), INTENT(OUT) :: U_no3, U_nh4
    REAL(wp) :: fr(0:NR), rc(NR), frn(0:NR), rcn(NR), b(NION), pool(NION), Dl(NION), Km(NION), Ig(NION)
    REAL(wp) :: tot, fl, TK, eta, eta298, v0, w, A_per_m, Imax, U(NION), Cmean, vol, fnew
    INTEGER  :: j, i
    U_no3 = 0.0_wp ; U_nh4 = 0.0_wp
    IF (k < 1 .OR. k > NLM) RETURN
    IF (Lv < nrad_lv_on) THEN
      act(k) = .FALSE.
      RETURN
    END IF
    b    = [MAX(th, 1.0e-3_wp), MAX(b_nh4, 1.0e-3_wp)]
    pool = [MAX(pool_no3, 0.0_wp), MAX(pool_nh4, 0.0_wp)]
    Km   = [Km_no3, Km_nh4]
    Ig   = [MAX(Imax_g_no3, 0.0_wp), MAX(Imax_g_nh4, 0.0_wp)]
    IF (.NOT. act(k)) THEN
      ! a layer joins with an undepleted (uniform) profile
      CALL build_grid(Lv, fr, rc)
      frL(:,k) = fr ; rcL(:,k) = rc ; LvL(k) = Lv
      DO j = 1, NION
        Cr(:, j, k) = 1.0_wp
      END DO
      act(k) = .TRUE.
    END IF
    fr = frL(:,k) ; rc = rcL(:,k)
    ! root growth: new root length starts in soil at the domain mean
    IF (ABS(Lv - LvL(k)) > 1.0e-9_wp * LvL(k)) THEN
      fnew = MAX(0.0_wp, (Lv - LvL(k)) / Lv)
      vol = PI * (fr(NR)**2 - fr(0)**2)
      DO j = 1, NION
        Cmean = 0.0_wp
        DO i = 1, NR
          Cmean = Cmean + PI * (fr(i)**2 - fr(i-1)**2) * Cr(i, j, k)
        END DO
        Cmean = Cmean / vol
        Cr(:, j, k) = (1.0_wp - fnew) * Cr(:, j, k) + fnew * Cmean
      END DO
      CALL build_grid(Lv, frn, rcn)
      fr = frn ; rc = rcn
      frL(:,k) = fr ; rcL(:,k) = rc ; LvL(k) = Lv
    END IF
    ! rescale each profile so the domain holds exactly the layer pool
    DO j = 1, NION
      tot = content(Cr(:, j, k), fr, b(j), Lv, dz)
      IF (tot > 1.0e-30_wp) THEN
        Cr(:, j, k) = Cr(:, j, k) * (pool(j) / tot)
      ELSE
        Cr(:, j, k) = pool(j) / (MAX(content([(1.0_wp, i = 1, NR)], fr, b(j), Lv, dz), 1.0e-30_wp))
      END IF
      C0(:, j, k) = Cr(:, j, k)
    END DO
    ! transport coefficients today
    TK = T_soil + 273.15_wp
    eta    = 2.414e-5_wp * 10.0_wp**(247.8_wp / (TK - 140.0_wp))
    eta298 = 2.414e-5_wp * 10.0_wp**(247.8_wp / (298.15_wp - 140.0_wp))
    ! Impedance: Olesen et al. 2001 Soil Sci 166:89-96 p.94 Eq.(7a)/(7b), f_l = 1.1 (theta - theta_th),
    ! theta_th from the deck's clay, silt and bulk density (p.93 Eq.5; 22 soils, sands to clays).
    ! Barraclough & Tinker 1981 (1.58 theta - 0.17) was fit on two loams only and returns f_l = 0
    ! at theta <= 0.11, i.e. no diffusion in a sandy soil at field capacity.
    fl = nud_csif_slope * MAX(0.0_wp, th - nud_theta_th)
    Dl = [D_NO3_25, D_NH4_25] * (TK / 298.15_wp) * (eta298 / eta) * th * fl
    v0 = (MAX(water_mm, 0.0_wp) * 1.0e-3_wp / 86400.0_wp) / (2.0_wp * PI * r0 * Lv * dz)   ! m s-1 at the root surface
    w  = 2.0_wp * PI * r0 * v0                                                             ! m2 s-1 per m root
    A_per_m = 2.0_wp * PI * r0 + n_h * 2.0_wp * PI * r_h * l_h                              ! m2 absorbing surface per m root
    DO j = 1, NION
      ! mol m-2 s-1 of absorbing surface: the per-gram capacity over the surface per gram
      Imax = Ig(j) / MN / 86400.0_wp / (nud_SRL * A_per_m)
      CALL radial_ion_day(Cr(:, j, k), fr, rc, b(j), Dl(j), w, Imax, Km(j), U(j))
      U(j) = MIN(U(j) * Lv * dz * MN, pool(j))                                             ! g N m-2
    END DO
    U_no3 = U(1) ; U_nh4 = U(2)
  END SUBROUTINE nrad_layer_day

  ! Take the fraction f (0..1) of today's potential uptake of layer k: the
  ! profile becomes C0 + f (C_after - C0), so its content falls by exactly f U.
  SUBROUTINE nrad_commit(k, f)
    INTEGER,  INTENT(IN) :: k
    REAL(wp), INTENT(IN) :: f
    IF (k < 1 .OR. k > NLM) RETURN
    IF (.NOT. act(k)) RETURN
    Cr(:, :, k) = C0(:, :, k) + MAX(0.0_wp, MIN(1.0_wp, f)) * (Cr(:, :, k) - C0(:, :, k))
  END SUBROUTINE nrad_commit

  ! Implicit radial solve for one ion over one day, 24 hourly steps.
  ! U returned in mol per m root.
  SUBROUTINE radial_ion_day(C, fr, rc, b, Dsoil, w, Imax, Km, U)
    REAL(wp), INTENT(INOUT) :: C(NR)
    REAL(wp), INTENT(IN)    :: fr(0:NR), rc(NR), b, Dsoil, w, Imax, Km
    REAL(wp), INTENT(OUT)   :: U
    REAL(wp) :: a(0:NR), Hh(NR), vol(NR), lo(NR), di(NR), up(NR), rh(NR), Cs(NR), Ucoef(NR), dt, U0, ahair
    INTEGER  :: i, it, ih
    U = 0.0_wp
    DO i = 1, NR - 1
      a(i) = 2.0_wp * PI * fr(i) * Dsoil / (rc(i+1) - rc(i))
    END DO
    a(0) = 0.0_wp ; a(NR) = 0.0_wp
    DO i = 1, NR
      vol(i) = PI * (fr(i)**2 - fr(i-1)**2)
      ahair  = MAX(0.0_wp, MIN(fr(i), r0 + l_h) - MAX(fr(i-1), r0))
      Hh(i)  = n_h * 2.0_wp * PI * r_h * ahair
    END DO
    dt = 3600.0_wp
    DO ih = 1, 24
      Cs = C
      DO it = 1, 6
        DO i = 1, NR
          Ucoef(i) = Imax / (Km + MAX(Cs(i), 0.0_wp))
        END DO
        DO i = 1, NR
          di(i) = vol(i) * b / dt + Hh(i) * Ucoef(i)
          rh(i) = vol(i) * b / dt * C(i)
          lo(i) = 0.0_wp ; up(i) = 0.0_wp
          IF (i < NR) THEN
            di(i) = di(i) + a(i)
            up(i) = -(a(i) + w)
          END IF
          IF (i > 1) THEN
            di(i) = di(i) + a(i-1) + w
            lo(i) = -a(i-1)
          ELSE
            di(i) = di(i) + 2.0_wp * PI * r0 * Ucoef(1)
          END IF
        END DO
        CALL tridiag(lo, di, up, rh, Cs)
        Cs = MAX(Cs, 0.0_wp)
      END DO
      U0 = 2.0_wp * PI * r0 * Ucoef(1) * Cs(1)
      DO i = 1, NR
        U0 = U0 + Hh(i) * Ucoef(i) * Cs(i)
      END DO
      U = U + U0 * dt
      C = Cs
    END DO
  END SUBROUTINE radial_ion_day

  SUBROUTINE tridiag(a, b, c, d, x)
    REAL(wp), INTENT(IN)  :: a(NR), b(NR), c(NR), d(NR)
    REAL(wp), INTENT(OUT) :: x(NR)
    REAL(wp) :: cp(NR), dp(NR), m
    INTEGER :: i
    cp(1) = c(1) / b(1) ; dp(1) = d(1) / b(1)
    DO i = 2, NR
      m = b(i) - a(i) * cp(i-1)
      cp(i) = c(i) / m
      dp(i) = (d(i) - a(i) * dp(i-1)) / m
    END DO
    x(NR) = dp(NR)
    DO i = NR - 1, 1, -1
      x(i) = dp(i) - cp(i) * x(i+1)
    END DO
  END SUBROUTINE tridiag

END MODULE NITROGEN_RADIAL_MODULE
