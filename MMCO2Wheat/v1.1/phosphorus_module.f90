! ======================================================================== !
! PHOSPHORUS_MODULE  (M&M-CO2Wheat v1.1, development branch v1.1-P-module)
!
! A mechanistic soil-plant phosphorus module. Design document:
!   v1.1_P_module/DESIGN_P_module_mechanistic_v1.1.md (Google Drive).
! Parameters, each with source, page and status: input/driver_P.csv.
!
! Principles (design section 0):
!  - Every flux is a conservation law, a transport law, a sorption equilibrium,
!    a transporter kinetic law or a stoichiometric requirement.
!  - There is no P stress factor. Growth is limited only because a process
!    cannot get the P atoms it needs (construction payment), because leaf
!    cytosolic Pi sets the triose-phosphate utilization capacity, and (phase B)
!    because root cytosolic Pi gates the root hydraulic conductance.
!  - General: no site-specific constant lives here. Soil properties are inputs.
!  - OFF MEANS IDENTICAL: when input/driver_P.csv is absent or its switch is 0,
!    p_on stays .FALSE., no routine below changes any v1.0 state, and every
!    v1.0 deck reproduces bit for bit.
!
! Scope of this first code version (phase A), stated so that nothing is
! claimed that is not coded:
!  - one root-zone soil compartment (as MINERAL_MODULE), resolved radially
!    around a representative root (Barber-Cushman, implicit finite volume,
!    nonlinear Langmuir buffer, Picard iteration, root hairs, hourly steps);
!  - fast sorption (Langmuir with the Barrow electrostatic pH term), slow
!    diffusive penetration (first eigenmode, Barrow 1983), broadcast
!    fertilizer, net organic-P mineralization as a site input;
!  - Pht1 transporter abundance from synthesis repressed by cytosolic Pi;
!  - plant pools: structural (phospholipid + DNA), ribosomal RNA (growth-rate
!    hypothesis), cytosolic Pi, vacuolar Pi; leaf resorption at senescence;
!    grain P demand reported;
!  - the P payment gate in front of the N payment gate (growth_module);
!  - TPU from leaf cytosolic Pi (assim.f90, through SHARED p_tpu_scale).
!  Not yet coded (phase B, see design): fertilizer band compartment, DCPD
!  precipitation, phosphatase route, vertical P transport, the tiller-site
!  window, primordium scaling, root hydraulic conductance.
! ======================================================================== !
MODULE PHOSPHORUS_MODULE
  USE TYPES, ONLY: wp
  USE SHARED, ONLY: p_module_on, p_tpu_scale, nud_SRL, &
                    n4_sln_target, csl_LMA, n4_stem_nmax, npg_build_R   ! the N gate's tissue N targets (protoplasm scale)
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: p_on, read_p_params, p_soil_day, p_gate, p_plant_end_day, &
            p_write_daily, p_close, p_pL, p_pS, p_pR, p_temp_today, p_report, &
            p_pL_last, p_pS_last

  LOGICAL  :: p_on = .FALSE.
  REAL(wp) :: p_temp_today = 20.0_wp        ! daily mean temperature for the growth gate, deg C (set by the caller)
  REAL(wp) :: p_pL = 1.0_wp, p_pS = 1.0_wp, p_pR = 1.0_wp   ! today's P-paid fractions
  REAL(wp) :: p_pL_last = 1.0_wp, p_pS_last = 1.0_wp        ! paid fractions of the last gate day (read by leaf appearance and tiller buds)

  ! ---- physical constants ---------------------------------------------- !
  REAL(wp), PARAMETER :: MP   = 30.973762_wp      ! g P mol-1
  REAL(wp), PARAMETER :: RGAS = 8.314462618_wp    ! J mol-1 K-1
  REAL(wp), PARAMETER :: FARA = 96485.33212_wp    ! C mol-1
  REAL(wp), PARAMETER :: PI   = 3.14159265358979_wp
  REAL(wp), PARAMETER :: PROT_PER_N = 6.25_wp     ! g protein g-1 N (Jones factor)

  INTEGER, PARAMETER :: NR = 25                   ! radial cells

  ! ---- parameters (read from input/driver_P.csv; defaults are the sourced
  !      values of PARAMETERS_P_module.csv and are overwritten by the file) - !
  REAL(wp) :: pKa2_298  = 7.20_wp       ! -
  REAL(wp) :: dH2       = 3600.0_wp     ! J mol-1
  REAL(wp) :: D_H2PO4   = 8.46e-10_wp   ! m2 s-1, Li & Gregory 1974 Table 1 p.707
  REAL(wp) :: D_HPO4    = 7.34e-10_wp   ! m2 s-1, Li & Gregory 1974 Table 1 p.707
  ! soil (site inputs)
  REAL(wp) :: soil_pH   = 6.5_wp
  REAL(wp) :: rho_b     = 1300.0_wp     ! kg m-3
  REAL(wp) :: phi_tot   = 0.50_wp       ! total porosity
  REAL(wp) :: S_max     = 5.0e-3_wp     ! mol P kg-1 (site isotherm)
  REAL(wp) :: K_0       = 50.0_wp       ! m3 mol-1 (site isotherm, at pH_0)
  REAL(wp) :: pH_0      = 6.5_wp        ! pH at which the isotherm was measured
  REAL(wp) :: z_eff     = 1.5_wp        ! effective charge of the adsorbing ion
  REAL(wp) :: m_psi     = 0.0_wp        ! V per pH unit; 0 = no pH correction (site input)
  REAL(wp) :: dH_a      = 0.0_wp        ! J mol-1, sorption enthalpy (site input)
  REAL(wp) :: k_slow    = 0.0_wp        ! d-1, pi^2 D_s/a_p^2 (site input; 0 = off)
  REAL(wp) :: beta_slow = 0.0_wp        ! interior/surface partition (site input)
  REAL(wp) :: soiltest_mgkg = 15.0_wp   ! soil-test P, mg P kg-1 (site input)
  REAL(wp) :: eps_test  = 1.0_wp        ! recovery of fast-sorbed P by the soil test (site input)
  REAL(wp) :: Pmin_g_m2_d = 0.0_wp      ! net organic-P mineralization, g P m-2 d-1 (site input)
  INTEGER, PARAMETER :: NFERT = 4
  INTEGER  :: fert_day(NFERT) = -1      ! model day index iD of each broadcast application
  REAL(wp) :: fert_gP_m2(NFERT) = 0.0_wp
  ! root uptake
  REAL(wp) :: r0        = 8.0e-5_wp     ! m, root radius; Gahoonia et al. 1997 Table 2 (diam. 0.15-0.17 mm)
  REAL(wp) :: r_h       = 6.0e-6_wp     ! m, hair radius; Gahoonia 1997 p.181 (diam. 12 um)
  REAL(wp) :: l_h       = 8.0e-4_wp     ! m, hair length; Gahoonia 1997 Table 2 (0.49-1.27 mm)
  REAL(wp) :: n_h       = 3.0e4_wp      ! hairs m-1 root; Gahoonia 1997 Table 2 (24-38 mm-1), Bole 1973 (45-60 mm-1)
  REAL(wp) :: Imax_ind_g = 9.5_wp       ! umol h-1 g-1 root DM, fully induced; Otani et al. 1996 Table 1 (cereals 5.9-13.1)
  REAL(wp) :: K_m       = 0.017_wp      ! mol m-3; Otani 1996 Table 1 (cereals 13.7-21.3 uM)
  REAL(wp) :: C_min     = 1.8e-3_wp     ! mol m-3; Otani 1996 Table 1 (cereals 1.70-1.95 uM)
  REAL(wp) :: Q10_upt   = 2.0_wp        ! provisional
  REAL(wp) :: T_upt_ref = 25.0_wp       ! deg C of the Otani measurements
  REAL(wp) :: induct_R  = 3.0_wp        ! induced/replete Imax ratio (provisional)
  REAL(wp) :: k_deg_tr  = 1.0_wp        ! d-1, transporter turnover (provisional)
  ! plant
  REAL(wp) :: c_cyt_star = 6.5_wp       ! mol m-3 cytosol; Lee et al. 1990 (maize roots); conflict with Pratt 2009 flagged
  REAL(wp) :: c_cyt_floor_frac = 0.1_wp ! lowest cytosolic Pi as a fraction of c* (provisional)
  REAL(wp) :: v_cyt     = 4.0e-7_wp     ! m3 cytosol g-1 DM (provisional: 4 g water g-1 DM x 10 % cytosol)
  REAL(wp) :: K_TPT     = 1.0_wp        ! mol m-3; order of the feedback-limited stromal Pi, Sharkey & Vanderveer 1989 (provisional)
  REAL(wp) :: pP_struct_L = 5.5e-4_wp   ! g P g-1 DM in phospholipid + DNA of new lamina (provisional, Veneklaas 2012 Fig. 1)
  REAL(wp) :: pP_struct_S = 2.0e-4_wp   ! g P g-1 DM, stem (provisional)
  REAL(wp) :: pP_struct_R = 3.0e-4_wp   ! g P g-1 DM, root (provisional)
  REAL(wp) :: fP_RNA    = 0.091_wp      ! g P g-1 RNA; Raven 2013 Table 1
  ! Composition of NEW tissue (design 6.2, revised 6 Nov 2026). Veneklaas et al.
  ! 2012 New Phytol 195:306, p.309: nucleic-acid P 0.3-2 mg P g-1 DW, >= 85 % RNA;
  ! pool ratio lipid P : nucleic-acid P : Pi : P-esters = 1.3 : 2.0 : 1.4 : 1.
  REAL(wp) :: pP_NA_new = 1.0e-3_wp     ! g P g-1 new DM in nucleic acids (Veneklaas 2012 range 0.3-2 mg g-1; mid value)
  REAL(wp) :: pP_lip    = 6.5e-4_wp     ! g P g-1 new DM in phospholipid = 1.3/2.0 x pP_NA_new (Veneklaas 2012 ratio)
  REAL(wp) :: pP_ester  = 5.0e-4_wp     ! g P g-1 new DM in P-esters = 1/2.0 x pP_NA_new (Veneklaas 2012 ratio)
  REAL(wp) :: f_DNA     = 0.15_wp       ! DNA share of nucleic acids (Veneklaas 2012: RNA >= 85 %)
  REAL(wp) :: f_RNA_mature = 0.25_wp    ! RNA kept in mature tissue (Veneklaas 2012: cytosolic rRNA drops < 75 %)
  REAL(wp) :: k_RNA     = 0.1_wp        ! d-1, RNA decay toward the mature level (provisional)
  REAL(wp) :: kps_20    = 17.5_wp       ! NOT USED since 6 Nov 2026: the maximum rate gave ~10x too little RNA against Veneklaas 2012 / Sulpice 2014; kept for the record
  REAL(wp) :: Q10_ps    = 2.0_wp        ! Raven 2013 p.5
  REAL(wp) :: eta_P     = 0.5_wp        ! leaf P resorption; Veneklaas 2012 p.313 ("at least 50 %")
  REAL(wp) :: grainP_conc = 3.5e-3_wp   ! g P g-1 grain DM (provisional)
  REAL(wp) :: seedP_g_m2 = 0.04_wp      ! g P m-2 in the sown seed (site input)

  ! ---- state --------------------------------------------------------------- !
  LOGICAL  :: initialized = .FALSE.
  REAL(wp) :: Cr(NR) = 0.0_wp          ! mol m-3 solution, radial profile
  REAL(wp) :: Ss(NR) = 0.0_wp          ! mol kg-1, slow (interior) pool per cell
  REAL(wp) :: fr(0:NR) = 0.0_wp        ! face radii (m)
  REAL(wp) :: rc(NR) = 0.0_wp          ! cell-centre radii (m)
  REAL(wp) :: Lv_prev = 0.0_wp         ! m m-3
  REAL(wp) :: depth_prev = 0.0_wp      ! m
  REAL(wp) :: I_rel = 1.0_wp           ! transporter abundance relative to replete
  REAL(wp) :: Psoil_free = 0.0_wp      ! g P m-2: soil P below the current root depth / before roots (well-mixed)
  ! plant P (g P m-2)
  REAL(wp) :: Pstr_L = 0.0_wp, Pstr_S = 0.0_wp, Pstr_R = 0.0_wp
  REAL(wp) :: Prna = 0.0_wp, Pcyt = 0.0_wp, Pvac = 0.0_wp
  REAL(wp) :: Pgrain = 0.0_wp, Pdead = 0.0_wp, Pgrain_unmet = 0.0_wp
  REAL(wp) :: Prna_req_today = 0.0_wp
  REAL(wp) :: c_cyt = 6.5_wp
  REAL(wp) :: W_live_today = 0.0_wp
  REAL(wp) :: W_dead_prev = -1.0_wp, W_grain_prev = 0.0_wp
  ! ledgers
  REAL(wp) :: U_day = 0.0_wp, U_cum = 0.0_wp, fert_cum = 0.0_wp, min_cum = 0.0_wp
  REAL(wp) :: Psoil_start = 0.0_wp, Pplant_start = 0.0_wp
  REAL(wp) :: soil_resid_max = 0.0_wp, plant_resid_max = 0.0_wp
  REAL(wp) :: paid_today = 0.0_wp, req_today = 0.0_wp
  INTEGER  :: unit_p = -1

CONTAINS

  ! ======================================================================= !
  ! Input
  ! ======================================================================= !
  SUBROUTINE read_p_params(fname)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER :: u, ios, k
    CHARACTER(LEN=512) :: line
    CHARACTER(LEN=64)  :: key
    REAL(wp) :: val
    LOGICAL  :: ex
    INQUIRE (FILE=fname, EXIST=ex)
    IF (.NOT. ex) RETURN                       ! no file: module off, v1.0 identical
    OPEN (NEWUNIT=u, FILE=fname, STATUS='old', ACTION='read')
    READ (u, '(A)', IOSTAT=ios) line          ! header
    DO
      READ (u, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(line) == 0 .OR. line(1:1) == '#') CYCLE
      k = INDEX(line, ',')
      IF (k <= 1) CYCLE
      key = ADJUSTL(line(1:k-1))
      READ (line(k+1:), *, IOSTAT=ios) val
      IF (ios /= 0) CYCLE
      CALL set_param(TRIM(key), val)
    END DO
    CLOSE (u)
    p_module_on = p_on
    IF (p_on) WRITE (*,'(A)') ' P module: ON (input/driver_P.csv)'
  END SUBROUTINE read_p_params

  SUBROUTINE set_param(key, v)
    CHARACTER(LEN=*), INTENT(IN) :: key
    REAL(wp), INTENT(IN) :: v
    INTEGER :: i
    SELECT CASE (key)
    CASE ('p_module_switch'); p_on = (NINT(v) >= 1)
    CASE ('pKa2_298');   pKa2_298 = v
    CASE ('dH2');        dH2 = v
    CASE ('D_H2PO4');    D_H2PO4 = v
    CASE ('D_HPO4');     D_HPO4 = v
    CASE ('soil_pH');    soil_pH = v
    CASE ('rho_b');      rho_b = v
    CASE ('phi_tot');    phi_tot = v
    CASE ('S_max');      S_max = v
    CASE ('K_0');        K_0 = v
    CASE ('pH_0');       pH_0 = v
    CASE ('z_eff');      z_eff = v
    CASE ('m_psi');      m_psi = v
    CASE ('dH_a');       dH_a = v
    CASE ('k_slow');     k_slow = v
    CASE ('beta_slow');  beta_slow = v
    CASE ('soiltest_mgkg'); soiltest_mgkg = v
    CASE ('eps_test');   eps_test = v
    CASE ('Pmin_g_m2_d'); Pmin_g_m2_d = v
    CASE ('r0');         r0 = v
    CASE ('r_h');        r_h = v
    CASE ('l_h');        l_h = v
    CASE ('n_h');        n_h = v
    CASE ('Imax_ind_g'); Imax_ind_g = v
    CASE ('K_m');        K_m = v
    CASE ('C_min');      C_min = v
    CASE ('Q10_upt');    Q10_upt = v
    CASE ('T_upt_ref');  T_upt_ref = v
    CASE ('induct_R');   induct_R = MAX(1.0_wp, v)
    CASE ('k_deg_tr');   k_deg_tr = v
    CASE ('c_cyt_star'); c_cyt_star = v
    CASE ('c_cyt_floor_frac'); c_cyt_floor_frac = v
    CASE ('v_cyt');      v_cyt = v
    CASE ('K_TPT');      K_TPT = v
    CASE ('pP_struct_L'); pP_struct_L = v
    CASE ('pP_struct_S'); pP_struct_S = v
    CASE ('pP_struct_R'); pP_struct_R = v
    CASE ('fP_RNA');     fP_RNA = v
    CASE ('k_RNA'); k_RNA = v
    CASE ('f_RNA_mature'); f_RNA_mature = v
    CASE ('f_DNA'); f_DNA = v
    CASE ('pP_ester'); pP_ester = v
    CASE ('pP_lip'); pP_lip = v
    CASE ('pP_NA_new'); pP_NA_new = v
    CASE ('kps_20');     kps_20 = v
    CASE ('Q10_ps');     Q10_ps = v
    CASE ('eta_P');      eta_P = v
    CASE ('grainP_conc'); grainP_conc = v
    CASE ('seedP_g_m2'); seedP_g_m2 = v
    CASE DEFAULT
      DO i = 1, NFERT
        IF (key == 'fert_day_'//CHAR(48+i))   fert_day(i)   = NINT(v)
        IF (key == 'fert_gP_m2_'//CHAR(48+i)) fert_gP_m2(i) = v
      END DO
    END SELECT
  END SUBROUTINE set_param

  ! ======================================================================= !
  ! Soil chemistry
  ! ======================================================================= !
  PURE REAL(wp) FUNCTION Ka_T(TK)
    REAL(wp), INTENT(IN) :: TK
    REAL(wp) :: psi
    psi  = -m_psi * (soil_pH - pH_0)          ! change of the plane potential from the isotherm pH (V)
    Ka_T = K_0 * EXP(-dH_a / RGAS * (1.0_wp/TK - 1.0_wp/298.15_wp)) * EXP(-z_eff * FARA * psi / (RGAS * TK))
  END FUNCTION Ka_T

  PURE REAL(wp) FUNCTION Sf_of(C, Ka)
    REAL(wp), INTENT(IN) :: C, Ka
    Sf_of = S_max * Ka * MAX(C, 0.0_wp) / (1.0_wp + Ka * MAX(C, 0.0_wp))
  END FUNCTION Sf_of

  PURE REAL(wp) FUNCTION Q_of(C, th, Ka)        ! mol P m-3 soil in solution + fast sites
    REAL(wp), INTENT(IN) :: C, th, Ka
    Q_of = th * MAX(C, 0.0_wp) + rho_b * Sf_of(C, Ka)
  END FUNCTION Q_of

  PURE REAL(wp) FUNCTION b_of(C, th, Ka)        ! buffer power dQ/dC
    REAL(wp), INTENT(IN) :: C, th, Ka
    b_of = th + rho_b * S_max * Ka / (1.0_wp + Ka * MAX(C, 0.0_wp))**2
  END FUNCTION b_of

  ! Invert Q(C) = Qt (monotone): closed form of the quadratic.
  PURE REAL(wp) FUNCTION C_of_Q(Qt, th, Ka)
    REAL(wp), INTENT(IN) :: Qt, th, Ka
    REAL(wp) :: a, b, c
    IF (Qt <= 0.0_wp) THEN
      C_of_Q = 0.0_wp ; RETURN
    END IF
    ! th*Ka*C^2 + (th + rho_b*S_max*Ka - Qt*Ka) C - Qt = 0
    a = th * Ka ; b = th + rho_b * S_max * Ka - Qt * Ka ; c = -Qt
    C_of_Q = (-b + SQRT(b*b - 4.0_wp*a*c)) / (2.0_wp*a)
  END FUNCTION C_of_Q

  PURE REAL(wp) FUNCTION alpha1(TK)            ! H2PO4- fraction
    REAL(wp), INTENT(IN) :: TK
    REAL(wp) :: pKa
    pKa = pKa2_298 + dH2 / (RGAS * LOG(10.0_wp)) * (1.0_wp/TK - 1.0_wp/298.15_wp)
    alpha1 = 1.0_wp / (1.0_wp + 10.0_wp**(soil_pH - pKa))
  END FUNCTION alpha1

  ! Liquid diffusion, species-weighted, Stokes-Einstein in T (water viscosity
  ! from the Vogel equation), times Millington-Quirk tortuosity (1961).
  PURE REAL(wp) FUNCTION De_of(th, TK)
    REAL(wp), INTENT(IN) :: th, TK
    REAL(wp) :: a1, DL, eta, eta298
    a1  = alpha1(TK)
    eta    = 2.414e-5_wp * 10.0_wp**(247.8_wp / (TK - 140.0_wp))
    eta298 = 2.414e-5_wp * 10.0_wp**(247.8_wp / (298.15_wp - 140.0_wp))
    DL  = (a1 * D_H2PO4 + (1.0_wp - a1) * D_HPO4) * (TK / 298.15_wp) * (eta298 / eta)
    De_of = DL * MAX(th, 1.0e-3_wp)**(10.0_wp/3.0_wp) / phi_tot**2
  END FUNCTION De_of

  ! ======================================================================= !
  ! Radial grid for root length density Lv (m m-3)
  ! ======================================================================= !
  SUBROUTINE build_grid(Lv)
    REAL(wp), INTENT(IN) :: Lv
    REAL(wp) :: r1, q
    INTEGER :: i
    r1 = 1.0_wp / SQRT(PI * MAX(Lv, 1.0_wp))
    r1 = MAX(r1, 2.0_wp * r0)
    q  = (r1 / r0)**(1.0_wp / REAL(NR, wp))
    fr(0) = r0
    DO i = 1, NR
      fr(i) = r0 * q**i
      rc(i) = SQRT(fr(i-1) * fr(i))
    END DO
    fr(NR) = r1
  END SUBROUTINE build_grid

  REAL(wp) FUNCTION cell_vol(i)                 ! m3 per m root
    INTEGER, INTENT(IN) :: i
    cell_vol = PI * (fr(i)**2 - fr(i-1)**2)
  END FUNCTION cell_vol

  ! Soil P in the radial domain, g P m-2 ground.
  REAL(wp) FUNCTION soil_P_domain(th, Ka, Lv, depth)
    REAL(wp), INTENT(IN) :: th, Ka, Lv, depth
    INTEGER :: i
    soil_P_domain = 0.0_wp
    DO i = 1, NR
      soil_P_domain = soil_P_domain + cell_vol(i) * (Q_of(Cr(i), th, Ka) + rho_b * Ss(i))
    END DO
    soil_P_domain = soil_P_domain * Lv * depth * MP
  END FUNCTION soil_P_domain

  ! ======================================================================= !
  ! First day: soil from the soil test, plant from the seed
  ! ======================================================================= !
  SUBROUTINE p_initialize(th, TK, Lv, depth, W_live)
    REAL(wp), INTENT(IN) :: th, TK, Lv, depth, W_live
    REAL(wp) :: Ka, Sf0, C0
    Ka  = Ka_T(TK)
    ! soil_test = eps * Sf(C0) * MP * 1000 (mg kg-1); invert the Langmuir.
    Sf0 = MIN(0.999_wp * S_max, soiltest_mgkg / (eps_test * MP * 1000.0_wp))
    C0  = Sf0 / (Ka * (S_max - Sf0))
    CALL build_grid(Lv)
    Cr = C0
    Ss = beta_slow * Sf0                       ! slow reaction at equilibrium
    Lv_prev = Lv ; depth_prev = depth
    I_rel = 1.0_wp
    ! seed P: cytosol at its set point in the emerged tissue, the rest vacuolar
    Pcyt = MIN(seedP_g_m2, c_cyt_star * v_cyt * W_live * MP)
    Pvac = seedP_g_m2 - Pcyt
    c_cyt = c_cyt_star
    Psoil_start  = soil_P_domain(th, Ka, Lv, depth)
    Pplant_start = seedP_g_m2
    initialized = .TRUE.
    OPEN (NEWUNIT=unit_p, FILE='output/daily_P.csv', STATUS='replace', ACTION='write')
    WRITE (unit_p,'(A)') 'day,C_bulk_uM,C_root_surface_uM,Sf_mmol_kg,Ss_mmol_kg,buffer_power,Lv_m_m3,' // &
      'I_rel,uptake_gP_m2_d,Pvac_gP_m2,Pcyt_gP_m2,c_cyt_mM,Prna_gP_m2,Pstr_gP_m2,Pgrain_gP_m2,' // &
      'Pdead_gP_m2,p_paid_L,p_paid_S,p_paid_R,tpu_scale,P_req_gP_m2_d,P_paid_gP_m2_d,' // &
      'soil_resid_gP_m2,plant_resid_gP_m2'
  END SUBROUTINE p_initialize

  ! ======================================================================= !
  ! Daily soil step and root uptake (call once per day before growth)
  !   iD        model day index
  !   trans_mm  actual transpiration, mm d-1
  !   theta     root-zone volumetric water content
  !   depth     root depth, m
  !   W_root    root dry mass, g m-2
  !   T_soil    root-zone temperature, deg C
  !   W_live    living plant dry mass, g m-2 (for the seed-P initialization)
  ! ======================================================================= !
  SUBROUTINE p_soil_day(iD, trans_mm, theta, depth, W_root, T_soil, W_live)
    INTEGER,  INTENT(IN) :: iD
    REAL(wp), INTENT(IN) :: trans_mm, theta, depth, W_root, T_soil, W_live
    REAL(wp) :: TK, Ka, th, Lv, dz, De, v0, w, Imax, Ucoef(NR), Q0(NR), Cs(NR), b(NR)
    REAL(wp) :: lo(NR), di(NR), up(NR), rh(NR), Hh(NR), a(0:NR), ahair, U0, dt, add, fnew
    REAL(wp) :: Qtot_old, Qtot_new, Sfi, dSs, before, after, Cm
    INTEGER  :: i, it, ih, k
    IF (.NOT. p_on) RETURN
    TK = T_soil + 273.15_wp
    th = MAX(0.02_wp, theta)
    dz = MAX(0.05_wp, depth)
    Lv = MAX(1.0e3_wp, W_root * nud_SRL / dz)   ! m root m-3 soil (floor: emergence)
    Ka = Ka_T(TK)
    IF (.NOT. initialized) CALL p_initialize(th, TK, Lv, dz, W_live)
    U_day = 0.0_wp
    W_live_today = W_live
    p_pL = 1.0_wp ; p_pS = 1.0_wp ; p_pR = 1.0_wp

    before = soil_P_domain(th, Ka, Lv_prev, depth_prev)

    ! --- root growth: new root length starts in undepleted soil (design 4.4);
    !     the deepening root zone brings new soil at the initial state of the
    !     domain mean, conserved by construction. -------------------------- !
    IF (ABS(Lv - Lv_prev) > 1.0e-9_wp * Lv_prev .OR. ABS(dz - depth_prev) > 1.0e-12_wp) THEN
      CALL remap_profile(th, Ka, Lv, dz)
    END IF
    before = soil_P_domain(th, Ka, Lv, dz)

    ! --- inputs: broadcast fertilizer and net mineralization, spread uniformly
    add = Pmin_g_m2_d
    DO k = 1, NFERT
      IF (fert_day(k) == iD) add = add + fert_gP_m2(k)
      IF (fert_day(k) == iD) fert_cum = fert_cum + fert_gP_m2(k)
    END DO
    min_cum = min_cum + Pmin_g_m2_d
    IF (add > 0.0_wp) THEN
      Cm = add / (MP * dz)                      ! mol P m-3 soil
      DO i = 1, NR
        Cr(i) = C_of_Q(Q_of(Cr(i), th, Ka) + Cm, th, Ka)
      END DO
    END IF

    ! --- slow penetration (Barrow 1983), daily, per cell -------------------- !
    IF (k_slow > 0.0_wp) THEN
      DO i = 1, NR
        Sfi = Sf_of(Cr(i), Ka)
        dSs = (Ss(i) + k_slow * beta_slow * Sfi) / (1.0_wp + k_slow) - Ss(i)
        Ss(i) = Ss(i) + dSs
        Cr(i) = C_of_Q(Q_of(Cr(i), th, Ka) - rho_b * dSs, th, Ka)
      END DO
    END IF

    ! --- transporter abundance: synthesis repressed by cytosolic Pi -------- !
    ! dI/dt = k_deg [R K/(K + c) - I],  K = c*/(R-1), so I = 1 at c = c* and
    ! I = R at c = 0 (fully induced). Implicit daily step.
    IF (induct_R > 1.0_wp) THEN
      I_rel = (I_rel + k_deg_tr * induct_R * (c_cyt_star/(induct_R-1.0_wp)) / &
              (c_cyt_star/(induct_R-1.0_wp) + c_cyt)) / (1.0_wp + k_deg_tr)
    END IF

    ! --- radial Barber-Cushman, implicit, hourly ------------------------------ !
    De = De_of(th, TK)
    ! water influx at the root surface, m s-1 (daily mean)
    v0 = (MAX(trans_mm, 0.0_wp) * 1.0e-3_wp / 86400.0_wp) / (2.0_wp * PI * r0 * Lv * dz)
    w  = 2.0_wp * PI * r0 * v0                   ! m2 s-1 per m root
    ! Imax per root surface, mol m-2 s-1: Otani's fully induced value per g
    ! root DM, divided by the root surface per g DM (SRL x 2 pi r0), scaled to
    ! the replete level by 1/R and back up by today's transporter abundance.
    Imax = Imax_ind_g * 1.0e-6_wp / 3600.0_wp / (nud_SRL * 2.0_wp * PI * r0) &
           / induct_R * I_rel * Q10_upt**((T_soil - T_upt_ref) / 10.0_wp)
    DO i = 1, NR - 1
      a(i) = 2.0_wp * PI * fr(i) * De / (rc(i+1) - rc(i))
    END DO
    a(0) = 0.0_wp ; a(NR) = 0.0_wp
    DO i = 1, NR                                  ! hair surface per m root in cell i
      ahair = MAX(0.0_wp, MIN(fr(i), r0 + l_h) - MAX(fr(i-1), r0))
      Hh(i)  = n_h * 2.0_wp * PI * r_h * ahair
    END DO
    dt = 3600.0_wp
    DO ih = 1, 24
      Q0 = [(Q_of(Cr(i), th, Ka), i = 1, NR)]
      Cs = Cr
      DO it = 1, 6                                ! Picard on b(C) and on the uptake kinetics
        DO i = 1, NR
          b(i) = b_of(Cs(i), th, Ka)
          IF (Cs(i) > C_min) THEN
            Ucoef(i) = Imax / (K_m + Cs(i) - C_min)
          ELSE
            Ucoef(i) = 0.0_wp
          END IF
        END DO
        DO i = 1, NR
          di(i) = cell_vol(i) * b(i) / dt + Hh(i) * Ucoef(i)
          rh(i) = cell_vol(i) / dt * (Q0(i) - Q_of(Cs(i), th, Ka) + b(i) * Cs(i)) + Hh(i) * Ucoef(i) * C_min
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
            rh(i) = rh(i) + 2.0_wp * PI * r0 * Ucoef(1) * C_min
          END IF
        END DO
        CALL tridiag(lo, di, up, rh, Cs)
        Cs = MAX(Cs, 0.0_wp)
      END DO
      ! uptake this hour (mol per m root), exactly the flux the solver removed
      U0 = 2.0_wp * PI * r0 * Ucoef(1) * MAX(Cs(1) - C_min, 0.0_wp)
      DO i = 1, NR
        U0 = U0 + Hh(i) * Ucoef(i) * MAX(Cs(i) - C_min, 0.0_wp)
      END DO
      U_day = U_day + U0 * dt * Lv * dz * MP     ! g P m-2
      Cr = Cs
    END DO

    after = soil_P_domain(th, Ka, Lv, dz)
    soil_resid_max = MAX(soil_resid_max, ABS((after - before) - (add - U_day)))
    U_cum = U_cum + U_day
    Lv_prev = Lv ; depth_prev = dz

    ! uptake enters the cytosol, which is refilled to its set point first;
    ! the rest goes to the vacuole (tonoplast transport treated as fast
    ! relative to a day, design 6.1)
    Pvac = Pvac + U_day
    CALL rebalance_cytosol(W_live)
  END SUBROUTINE p_soil_day

  ! New roots and a deeper root zone (design 4.4): the representative profile
  ! is the length-weighted mix of the old depleted profile and an undepleted
  ! one at the domain mean; then rescaled so the soil P total is conserved,
  ! with any new soil volume entering at the domain mean.
  SUBROUTINE remap_profile(th, Ka, Lv, dz)
    REAL(wp), INTENT(IN) :: th, Ka, Lv, dz
    REAL(wp) :: Qold(NR), Qmean, Smean, fnew, xo(0:NR), tot_old, tot_new, s
    INTEGER :: i
    tot_old = soil_P_domain(th, Ka, Lv_prev, depth_prev)
    Qmean = 0.0_wp ; Smean = 0.0_wp
    DO i = 1, NR
      Qold(i) = Q_of(Cr(i), th, Ka)
      Qmean = Qmean + cell_vol(i) * Qold(i)
      Smean = Smean + cell_vol(i) * Ss(i)
    END DO
    Qmean = Qmean / (PI * (fr(NR)**2 - fr(0)**2))
    Smean = Smean / (PI * (fr(NR)**2 - fr(0)**2))
    fnew = MAX(0.0_wp, (Lv*dz - Lv_prev*depth_prev) / (Lv*dz))
    xo = fr
    CALL build_grid(Lv)
    ! the cells keep their index (same relative log position); mix with fresh soil
    DO i = 1, NR
      Qold(i) = (1.0_wp - fnew) * Qold(i) + fnew * Qmean
      Ss(i)   = (1.0_wp - fnew) * Ss(i)   + fnew * Smean
      Cr(i)   = C_of_Q(Qold(i), th, Ka)
    END DO
    ! the new soil volume (deeper roots) enters at the domain mean density:
    ! total = old total + mean density x added volume
    tot_new = soil_P_domain(th, Ka, Lv, dz)
    s = (tot_old + MAX(0.0_wp, dz - depth_prev) * (Qmean + rho_b*Smean) * MP) / MAX(tot_new, 1.0e-30_wp)
    DO i = 1, NR
      Cr(i) = C_of_Q(s * Q_of(Cr(i), th, Ka), th, Ka)
      Ss(i) = s * Ss(i)
    END DO
    Psoil_start = Psoil_start + MAX(0.0_wp, dz - depth_prev) * (Qmean + rho_b*Smean) * MP
  END SUBROUTINE remap_profile

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

  SUBROUTINE rebalance_cytosol(W_live)
    REAL(wp), INTENT(IN) :: W_live
    REAL(wp) :: cap, mv
    cap = c_cyt_star * v_cyt * MAX(W_live, 1.0e-6_wp) * MP
    mv  = MIN(Pvac, MAX(0.0_wp, cap - Pcyt))
    Pcyt = Pcyt + mv ; Pvac = Pvac - mv
    IF (Pcyt > cap) THEN                          ! surplus back to the vacuole
      Pvac = Pvac + (Pcyt - cap) ; Pcyt = cap
    END IF
    c_cyt = Pcyt / (v_cyt * MAX(W_live, 1.0e-6_wp) * MP)
    ! TPU capacity (design 6.7): TPU = TPU_max c/(K+c); the v1.0 value 0.11 Vcmax
    ! is TPU at c*, so TPU_max = 0.11 Vcmax (K+c*)/c*, and the scale is exact 1 at c*.
    p_tpu_scale = (c_cyt / (K_TPT + c_cyt)) * ((K_TPT + c_cyt_star) / c_cyt_star)
  END SUBROUTINE rebalance_cytosol

  ! ======================================================================= !
  ! P payment gate (design 6.3). Called from n_payment_gate_apply on the
  ! carbon-side potentials, before the nitrogen is paid. Returns the paid
  ! fractions p_pL, p_pS, p_pR (module variables) and the scaled increments.
  ! ======================================================================= !
  SUBROUTINE p_gate(dw_L, dw_S, dw_R)
    REAL(wp), INTENT(INOUT) :: dw_L, dw_S, dw_R
    REAL(wp) :: W_live
    REAL(wp) :: cyt_req, req_g, avail, floor_cyt, payL, payS, payR, sS, sR
    p_pL = 1.0_wp ; p_pS = 1.0_wp ; p_pR = 1.0_wp
    paid_today = 0.0_wp ; req_today = 0.0_wp
    IF (.NOT. p_on .OR. .NOT. initialized) RETURN
    W_live = W_live_today
    ! P per g of new DM (design 6.2): the measured composition of growing
    ! tissue - phospholipid membranes, nucleic acids (mostly ribosomal RNA, the
    ! growth-rate term: fast growth means much young, ribosome-rich tissue),
    ! phosphate esters - plus the cytosolic Pi the new cells hold.
    cyt_req = c_cyt_star * v_cyt * MP
    req_g = pP_lip + pP_NA_new + pP_ester + cyt_req
    floor_cyt = c_cyt_floor_frac * c_cyt_star * v_cyt * MAX(W_live, 1.0e-6_wp) * MP
    avail = Pvac + MAX(0.0_wp, Pcyt - floor_cyt)
    ! The composition above is that of lamina-like (protoplasm-rich) tissue. P,
    ! like N, sits in protoplasm, not in cell walls, so each organ's cost is
    ! scaled by its N construction concentration relative to the lamina (the
    ! model's own N-gate targets): one protoplasm N:P for all organs.
    sS = n4_stem_nmax / (n4_sln_target / csl_LMA)
    sR = npg_build_R  / (n4_sln_target / csl_LMA)
    payL = dw_L * req_g ; payS = dw_S * req_g * sS ; payR = dw_R * req_g * sR
    req_today = payL + payS + payR
    IF (req_today <= 0.0_wp) RETURN
    ! phloem priority: expanding blades, stem and chaff, root (design 6.3)
    IF (payL > 0.0_wp) p_pL = MIN(1.0_wp, avail / payL)
    avail = MAX(0.0_wp, avail - p_pL * payL)
    IF (payS > 0.0_wp) p_pS = MIN(1.0_wp, avail / payS)
    avail = MAX(0.0_wp, avail - p_pS * payS)
    IF (payR > 0.0_wp) p_pR = MIN(1.0_wp, avail / payR)
    dw_L = p_pL * dw_L ; dw_S = p_pS * dw_S ; dw_R = p_pR * dw_R
    paid_today = p_pL * payL + p_pS * payS + p_pR * payR
    p_pL_last = p_pL ; p_pS_last = p_pS
    ! commit: draw from the vacuole, then the cytosol above its floor; place the
    ! P in membranes + DNA + esters (structural), RNA, and the new cytosol
    CALL draw(paid_today)
    Pstr_L = Pstr_L + dw_L * (pP_lip + f_DNA*pP_NA_new + pP_ester)
    Pstr_S = Pstr_S + dw_S * sS * (pP_lip + f_DNA*pP_NA_new + pP_ester)
    Pstr_R = Pstr_R + dw_R * sR * (pP_lip + f_DNA*pP_NA_new + pP_ester)
    Prna   = Prna + (dw_L + dw_S*sS + dw_R*sR) * (1.0_wp - f_DNA) * pP_NA_new
    Pcyt   = Pcyt + (dw_L + dw_S*sS + dw_R*sR) * cyt_req
  END SUBROUTINE p_gate

  ! Draw P for construction: surplus ribosomes first, then vacuole, then cytosol.
  SUBROUTINE draw(amount)
    REAL(wp), INTENT(IN) :: amount
    REAL(wp) :: rest, t
    rest = MAX(0.0_wp, amount)
    t = MIN(rest, Pvac) ; Pvac = Pvac - t ; rest = rest - t
    t = MIN(rest, Pcyt) ; Pcyt = Pcyt - t ; rest = rest - t
  END SUBROUTINE draw

  ! ======================================================================= !
  ! End of day: senescence resorption, ribosome turnover, grain P, cytosol
  ! ======================================================================= !
  SUBROUTINE p_plant_end_day(W_green, W_stem, W_root, W_dead, W_grain)
    REAL(wp), INTENT(IN) :: W_green, W_stem, W_root, W_dead, W_grain
    REAL(wp) :: d_dead, frac, rel, lost, W_live, dg, dem, t, ptot
    IF (.NOT. p_on .OR. .NOT. initialized) RETURN
    W_live = W_green + W_stem + W_root
    IF (W_dead_prev < 0.0_wp) W_dead_prev = W_dead
    d_dead = MAX(0.0_wp, W_dead - W_dead_prev) ; W_dead_prev = W_dead
    ! senescing lamina: its structural, ribosomal and cytosolic P is hydrolyzed;
    ! eta_P of it is exported, the rest stays in the dead leaf
    IF (d_dead > 0.0_wp .AND. W_green + d_dead > 0.0_wp) THEN
      frac = d_dead / (W_green + d_dead)
      rel  = frac * Pstr_L
      Pstr_L = Pstr_L - rel
      rel  = rel + frac * Prna * W_green / MAX(W_live, 1.0e-9_wp)
      Prna = Prna - frac * Prna * W_green / MAX(W_live, 1.0e-9_wp)
      lost = (1.0_wp - eta_P) * rel
      Pdead = Pdead + lost
      Pvac  = Pvac + (rel - lost)
    END IF
    ! maturing tissue degrades part of its ribosomes (RNase); the Pi returns
    ! to the vacuole and the mature tissue keeps f_RNA_mature of its RNA
    t = (1.0_wp - f_DNA) * pP_NA_new * f_RNA_mature * W_live
    IF (Prna > t) THEN
      rel = (Prna - t) * (1.0_wp - EXP(-k_RNA))
      Prna = Prna - rel ; Pvac = Pvac + rel
    END IF
    ! grain P (phytate) demand, met from the free pools; reported if unmet
    dg = MAX(0.0_wp, W_grain - W_grain_prev) ; W_grain_prev = W_grain
    dem = dg * grainP_conc
    t = MIN(dem, Pvac) ; Pvac = Pvac - t ; Pgrain = Pgrain + t
    ! grain phytate is loaded from the phloem-mobile Pi (vacuole, which also
    ! receives the resorbed P of senescing tissue); the working cytosol of the
    ! leaves is not drained (design 6.8). Unmet demand is reported.
    Pgrain_unmet = Pgrain_unmet + (dem - t)
    CALL rebalance_cytosol(W_live)
    ! plant mass balance: everything in = seed + uptake
    ptot = Pstr_L + Pstr_S + Pstr_R + Prna + Pcyt + Pvac + Pgrain + Pdead
    plant_resid_max = MAX(plant_resid_max, ABS(ptot - (Pplant_start + U_cum)))
  END SUBROUTINE p_plant_end_day

  SUBROUTINE p_write_daily(iD, theta, T_soil)
    INTEGER,  INTENT(IN) :: iD
    REAL(wp), INTENT(IN) :: theta, T_soil
    REAL(wp) :: Ka, Cb, Sfb, Ssb, vt, ptot
    INTEGER :: i
    IF (.NOT. p_on .OR. .NOT. initialized) RETURN
    Ka = Ka_T(T_soil + 273.15_wp)
    Cb = 0.0_wp ; Sfb = 0.0_wp ; Ssb = 0.0_wp ; vt = 0.0_wp
    DO i = 1, NR
      Cb  = Cb  + cell_vol(i) * Cr(i)
      Sfb = Sfb + cell_vol(i) * Sf_of(Cr(i), Ka)
      Ssb = Ssb + cell_vol(i) * Ss(i)
      vt  = vt  + cell_vol(i)
    END DO
    ptot = Pstr_L + Pstr_S + Pstr_R + Prna + Pcyt + Pvac + Pgrain + Pdead
    WRITE (unit_p,'(I0,23(",",ES13.5))') iD, Cb/vt*1.0e3_wp, Cr(1)*1.0e3_wp, Sfb/vt*1.0e3_wp, Ssb/vt*1.0e3_wp, &
      b_of(Cb/vt, MAX(theta,0.02_wp), Ka), Lv_prev, I_rel, U_day, Pvac, Pcyt, c_cyt, Prna, &
      Pstr_L + Pstr_S + Pstr_R, Pgrain, Pdead, p_pL, p_pS, p_pR, p_tpu_scale, req_today, paid_today, &
      soil_resid_max, ABS(ptot - (Pplant_start + U_cum))
  END SUBROUTINE p_write_daily

  SUBROUTINE p_report()
    IF (.NOT. p_on) RETURN
    WRITE (*,'(A)') ' ---- P module season summary ----'
    WRITE (*,'(A,F10.4)') '  P uptake, g P m-2           : ', U_cum
    WRITE (*,'(A,F10.4)') '  fertilizer P, g P m-2       : ', fert_cum
    WRITE (*,'(A,F10.4)') '  grain P, g P m-2            : ', Pgrain
    WRITE (*,'(A,F10.4)') '  grain P unmet, g P m-2      : ', Pgrain_unmet
    WRITE (*,'(A,ES10.2)') '  max soil mass residual      : ', soil_resid_max
    WRITE (*,'(A,ES10.2)') '  max plant mass residual     : ', plant_resid_max
  END SUBROUTINE p_report

  SUBROUTINE p_close()
    IF (unit_p > 0) CLOSE (unit_p)
  END SUBROUTINE p_close

END MODULE PHOSPHORUS_MODULE
