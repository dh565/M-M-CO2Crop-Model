! ========================================================================== !
! M&M-CO2Wheat v1.0  |  SHARED.f90
!
! Purpose: The model's single global store of parameters and state. It declares
!          the physical constants, the measured plant, soil, canopy and soil-
!          nitrogen parameters with their sources, the driver-record switches of
!          input/driver_main.txt, every prognostic pool and every diagnostic
!          ledger, and six small helper routines (state synchronization, hourly
!          mass binning, Gauss-weighted hour-to-day summation, driver-comment
!          stripping and an event-reachability check). No process law is solved
!          here. The state it holds is written and read on all of the model's
!          time steps - the Richards adaptive sub-step, the hourly canopy and
!          energy-balance step, the daily crop and soil step, and the season
!          ledgers - and on all of its spatial units: the leaf and the sunlit /
!          shaded leaf class, the leaf rank and tiller axis, the plant organ,
!          the soil layer, and the point as a whole.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: TYPES, ONLY: wp, ik - the working real kind and the integer kind; the
!       module has no other dependency, so it sits at the bottom of the build.
! Output: the module IS the shared state. It publishes the physical and
!       empirical constants, the per-record configuration switches read from
!       driver_main.txt, the soil hydraulic and thermal fields, the canopy,
!       organ, nitrogen and mineral pools, the energy-balance and Richards
!       ledgers, and the output unit numbers. It is consumed by
!       hydro, et and th (soil water, evaporation, daylength),
!       energy_balance_module, radiation_module and soil_heat_module (the
!       surface energy balance and soil heat), solar, meteo and CANOPY_MODEL
!       (radiation, weather disaggregation and the hourly canopy), assim (leaf
!       gas exchange), growth_module, soiln_module and
!       mineral_module (crop carbon and nitrogen, residue nitrogen and grain
!       minerals), count_lines (input sizing), and main_MMCO2Crop, which fills
!       it from the decks and writes it out.
! Called by: USE SHARED appears in 15 of the 20 compiled units; csv_file.f90,
!       csv_file_1d.f90, phenology_module.f90 and TYPES.f90 do not use it.
!       The six contained routines are called from
!       main_MMCO2Crop.f90 (sync_hydro_step_state, strip_driver_comment,
!       check_event_reachability, dHour_to_dDay) and from hydro.f90 and
!       energy_balance_module.f90 (bin_series_add, bin_series_get).
! ========================================================================== !
! ========================================================================== !
MODULE SHARED
!----------------------------------------------------------------------!
! Global parameters, state and diagnostic ledgers of M&M-CO2Wheat v1.0.
!----------------------------------------------------------------------!
  USE TYPES, ONLY: wp, ik
  IMPLICIT NONE
  SAVE

  TYPE :: HydroStepState
    LOGICAL :: failed = .FALSE.
    LOGICAL :: min_dt_reached = .FALSE.
    INTEGER :: dt_halvings = 0
    INTEGER :: cap_active = 0
    REAL(wp) :: last_dt = 0.0_wp
    REAL(wp) :: mb_resid_m = 0.0_wp
    REAL(wp) :: mb_correction_m = 0.0_wp
    REAL(wp) :: top_flux_m = 0.0_wp
    REAL(wp) :: top_infil_m = 0.0_wp
    REAL(wp) :: bottom_flux_m = 0.0_wp
    REAL(wp) :: root_sink_m = 0.0_wp
    REAL(wp) :: flux_top_actual_m_d = 0.0_wp
    REAL(wp) :: flux_top_req_m_d = 0.0_wp
    REAL(wp) :: flux_top_cap_m_d = 0.0_wp
    REAL(wp) :: flux_top_applied_m_d = 0.0_wp
    REAL(wp) :: infil_cap_m_d = 0.0_wp
    REAL(wp) :: evap_cap_m_d = 0.0_wp
    REAL(wp) :: evap_rate_pot_m_d = 0.0_wp
    REAL(wp) :: evap_rate_day_m_d = 0.0_wp
    REAL(wp) :: sum_avail_m = 0.0_wp
    REAL(wp) :: k1_m_d = 0.0_wp
    REAL(wp) :: h1_m = 0.0_wp
    REAL(wp) :: theta1 = 0.0_wp
    REAL(wp) :: h_top_m = 0.0_wp
    REAL(wp) :: h_crit_atm_m = 0.0_wp
    REAL(wp) :: q34_m_d = 0.0_wp
    REAL(wp) :: k34_m_d = 0.0_wp
    REAL(wp) :: supply_step_mm = 0.0_wp
    REAL(wp) :: infil_step_mm = 0.0_wp
    REAL(wp) :: infil_from_supply_step_mm = 0.0_wp
    REAL(wp) :: infil_from_pond_step_mm = 0.0_wp
    REAL(wp) :: soil_evap_step_mm = 0.0_wp
    REAL(wp) :: pond_evap_step_mm = 0.0_wp
    REAL(wp) :: runoff_step_mm = 0.0_wp
    REAL(wp) :: ponding_start_step_mm = 0.0_wp
    REAL(wp) :: ponding_end_step_mm = 0.0_wp
  END TYPE HydroStepState

  TYPE :: HydroCumLedger
    REAL(wp) :: top_flux_req_mm = 0.0_wp
    REAL(wp) :: top_flux_act_mm = 0.0_wp
    REAL(wp) :: infiltration_mm = 0.0_wp
    REAL(wp) :: evaporation_mm = 0.0_wp
    REAL(wp) :: runoff_mm = 0.0_wp
    REAL(wp) :: root_uptake_pot_mm = 0.0_wp
    REAL(wp) :: root_uptake_act_mm = 0.0_wp
    REAL(wp) :: bottom_flux_mm = 0.0_wp
    REAL(wp) :: pond_evap_mm = 0.0_wp
    REAL(wp) :: storage_mm = 0.0_wp
  END TYPE HydroCumLedger

  !--------------------------------------------------------------------!
  ! 1. Mathematical & Physical Constants
  !--------------------------------------------------------------------!
  REAL(wp), PARAMETER :: zero = 0.0_wp, one = 1.0_wp
  REAL(wp), PARAMETER :: Param_pi = 3.141592653589793_wp

  ! Thermodynamics & Air Properties
  REAL(wp), PARAMETER :: const_SB    = 5.67E-8_wp       ! Stefan-Boltzman constant
  REAL(wp), PARAMETER :: gamma_const = 66.0_wp          ! Phsycrometric constant (Pa K-1)
  REAL(wp), PARAMETER :: rho         = 1.209_wp         ! Density of moist air (kg m-3)
  REAL(wp), PARAMETER :: Cp          = 1010.0_wp        ! Specific heat capacity (J kg-1 K-1)
  REAL(wp), PARAMETER :: mair = 28.9655_wp
  REAL(wp), PARAMETER :: gasc = 8.314510_wp             ! Gas constant (J mol-1 K-1)
  REAL(wp), PARAMETER :: rgas = 1000.0_wp * gasc / mair

  !--------------------------------------------------------------------!
  ! 2. Model Parameters
  !--------------------------------------------------------------------!
  ! Climate
  REAL(wp), PARAMETER :: P_temp = 1.5_wp                ! Offset 1.5h after sunrise

  ! Crop and soil optical properties. The released configuration takes the
  ! soil and canopy shortwave properties from radiation_module; alpha_mix_s is
  ! the Beer-Lambert blend of the two, retained for the bulk-radiation path.
  !   alpha_soil_dry_sandy = 0.28, dry bare sandy soil. Karnieli et al. (1996)
  !     Remote Sens Environ 56:1-12, Table 2 (measured 0.25-0.32, Negev).
  !   alpha_canopy_wheat   = 0.20, green cereal canopy. Monteith & Unsworth
  !     (2008) Principles of Environmental Physics, 4th ed., p. 48 Table 3.1
  !     ("grass/cereal crops 0.18-0.25").
  !   k_alpha_canopy       = 0.5, the extinction coefficient of the blend; the
  !     same Beer-Lambert form as k_dr in Goudriaan (1977).
  REAL(wp), PARAMETER :: alpha_soil_dry_sandy = 0.28_wp
  REAL(wp), PARAMETER :: alpha_canopy_wheat   = 0.20_wp
  REAL(wp), PARAMETER :: k_alpha_canopy       = 0.5_wp
  REAL(wp) :: alpha_mix_s = 0.25_wp   ! computed each hourly step from LAI

  ! LAI-dependent roughness length and displacement height.
  !   raupach_cd1 = 7.5 and raupach_cd = 0.2: Raupach (1994) Boundary-Layer
  !     Meteorol 71:211-216, Table 1 (cd is the mean drag coefficient per unit
  !     frontal area).
  !   frontal_lai_factor = 0.5: frontal-area index lambda_f = 0.5 LAI for the
  !     erectophile wheat canopy; Raupach et al. (1991) Boundary-Layer Meteorol
  !     56:351-369, Eq. 2.
  REAL(wp), PARAMETER :: raupach_cd1         = 7.5_wp
  REAL(wp), PARAMETER :: raupach_cd          = 0.2_wp
  REAL(wp), PARAMETER :: frontal_lai_factor  = 0.5_wp
  REAL(wp) :: frontal_area_index_s = 0.0_wp
  REAL(wp) :: disp_height_s        = 0.0_wp
  REAL(wp) :: rough_len_s          = 0.0_wp
  REAL(wp), PARAMETER :: sur_emiss = 1.0_wp
  REAL(wp), PARAMETER :: const_P_sr = 0.11_wp, const_P_par = 0.04_wp
  REAL(wp), PARAMETER :: const_alpha_sr = 0.5_wp, const_alpha_par = 0.8_wp
  REAL(wp), PARAMETER :: k_lw = 0.7_wp ! LW extinction for canopy (Campbell & Norman, 1998)
  REAL(wp), PARAMETER :: emiss_canopy = 0.98_wp, emiss_soil = 0.95_wp ! Typical emissivities

  ! Saturated bare-soil albedo of the component shortwave path; the dry
  ! end-member is alpha_soil_dry_sandy = 0.28 above.
  REAL(wp) :: alpha_soil_wet = 0.12_wp  ! saturated bare-soil albedo (Karnieli 1996 range ~0.10-0.14)
  REAL(wp), PARAMETER :: soil_k_t = 0.7_wp ! Soil thermal conductivity (W m-1 K-1), sandy soil
  REAL(wp), PARAMETER :: soil_ref_depth_m = 0.10_wp ! Reference depth for G (m)

  ! Aerodynamic
  REAL(wp), PARAMETER :: z0_soil = 0.004_wp
  REAL(wp), PARAMETER :: k_cons = 0.4_wp
  REAL(wp), PARAMETER :: Dmv = 2.47E-5_wp               ! Vapor diffusion (m2 s-1)
  ! ------------------------------------------------------------------------- !
  ! Photosynthesis (Farquhar-von Caemmerer-Berry) and respiration.
  ! Rubisco kinetics: Bernacchi et al. (2001) Plant Cell Environ 24:253-259,
  ! Eq. (8) p.254, Parameter = exp(c - dHa/(R Tk)), with c and dHa from Table 1
  ! p.256. Values at 25 C: Kc = 404.9 umol mol-1, Ko = 278.4 mmol mol-1,
  ! Gamma* = 42.75 umol mol-1. The unit basis is the mole fraction (p.253 and
  ! Table 1), which is also the basis of the model Ci and Oa, so no pressure
  ! conversion is applied; the published Kc is the apparent, mesophyll-
  ! conductance-inclusive constant (p.258), consistent with the Ci-based
  ! formulation used here. The Vcmax and Rd responses are normalized to 1.0 at
  ! 25 C (p.255) and applied as exp(dHa/R (1/298.15 - 1/Tk)) about the model's own
  ! 25 C values. Validity 10-40 C; no de-activation term is needed at or below
  ! 40 C (p.256).
  ! ------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: b01_c_Kc = 38.05_wp,  b01_dHa_Kc = 79430.0_wp
  REAL(wp), PARAMETER :: b01_c_Ko = 20.30_wp,  b01_dHa_Ko = 36380.0_wp
  REAL(wp), PARAMETER :: b01_c_Gs = 19.02_wp,  b01_dHa_Gs = 37830.0_wp
  REAL(wp), PARAMETER :: b01_dHa_Vc = 65330.0_wp
  REAL(wp), PARAMETER :: b01_dHa_Rd = 46390.0_wp
  REAL(wp), PARAMETER :: b01_Rgas = 8.314_wp    ! J mol-1 K-1 (B01 Eq. 8 R)
  ! ------------------------------------------------------------------------- !
  ! Jmax temperature response: Bernacchi, Pimentel & Long (2003) Plant Cell
  ! Environ 26:1419-1430, Eqs. (9)-(10) p.1421,
  !   Parameter = Parameter25 exp(c - dHa/(R Tk)) / [1 + exp((dS Tk - dHd)/(R Tk))].
  ! Parameters from Table 1 p.1424, the gas-exchange-derived Jmax column at
  ! growth temperature Tg = 14 C: c = 23.2, dHa = 57.5 kJ mol-1, dHd = 439.8
  ! kJ mol-1, dS = 1.4 kJ mol-1 K-1. The Tg row is selected by the measured
  ! mean air temperature of the growing phase, not fitted; the paper central
  ! result (Figs. 2-3) is that the Jmax response depends on Tg. Applied
  ! normalized, f(T)/f(25), about the model's own SLN-based Jmax25.
  ! ------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: b03_c_J   = 23.2_wp
  REAL(wp), PARAMETER :: b03_dHa_J = 57500.0_wp    ! J/mol
  REAL(wp), PARAMETER :: b03_dHd_J = 439800.0_wp   ! J/mol
  REAL(wp), PARAMETER :: b03_dS_J  = 1400.0_wp     ! J/mol/K

  ! Oxygen mole fraction of air, read by assim.f90, and the leaf absorptance
  ! constant, which no statement reads (kept for reference):
  REAL(wp), PARAMETER :: Oa = 210000.0_wp               ! O2 concentration
  REAL(wp), PARAMETER :: alpha = 0.8_wp                 ! Leaf absorbtion

  ! Nitrogen and photorespiration. Busch, Sage & Sage (2018) Nature Plants
  ! 4:46-54, doi:10.1038/s41477-017-0065-x, Eqs. (16)-(21) p.49.
  REAL(wp), PARAMETER :: alpha_g_max = 0.09_wp, alpha_s_max = 0.38_wp
  REAL(wp), PARAMETER :: f_photons = 0.23_wp
  ! ------------------------------------------------------------------------- !
  ! Stomatal conductance, Leuning (1995) form
  !   gsw = g0 + g1 An h(Ds) / (Ca - Gamma*),  h(Ds) = 1/(1 + Ds/D0).
  ! g1_leu = 10.58: the wheat refit of the stomatal database of Lin et al.
  ! (2015) Nature Climate Change 5:459-464, doi:10.1038/nclimate2550, with the
  ! h(Ds) term present (95 % CI +/- 0.34, RMSE 0.072).
  ! D0 = 1500 Pa: Leuning (1995) Plant Cell Environ 18:339-355,
  ! doi:10.1111/j.1365-3040.1995.tb00370.x, D0 defined at Eq. 7 p.341 and used
  ! at 1500 Pa across the six species rows of Table 2 p.347. The paper contains
  ! no wheat, so D0 is a cross-species value, not a wheat measurement.
  ! ------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: g0_leu = 0.01_wp, g1_leu = 10.58_wp
  ! ------------------------------------------------------------------------- !
  ! g1_tuzet = 5.776: the slope of the soil-plant-atmosphere (Tuzet) branch,
  !   gsw = g0 + f_psi(psi_L) g1_tuzet An / (Ca - Gamma*),
  ! which carries no h(Ds) term. Same source and same wheat subset as g1_leu
  ! above - Lin et al. (2015) Nature Climate Change 5:459-464 - refitted
  ! WITHOUT the h(Ds) term (g1 = 5.7763, 95 % CI +/- 0.240, RMSE 0.0926).
  ! Dropping h(Ds) is required by the source: Duursma & Medlyn (2012) Geosci
  ! Model Dev 5:919-940, p.922, "Tuzet et al. (2003) use a flexible psi_L
  ! modifier (in place of the f(D) function in Eq. 2)". Their own Table 3
  ! p.930 gives g1 (Tuzet model) = 5.5, independent corroboration of the
  ! refit. The leaf-water-potential pathway now carries the humidity response.
  ! ------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: g1_tuzet = 5.776_wp
  REAL(wp), PARAMETER :: D0 = 1500.0_wp
  ! ------------------------------------------------------------------------- !
  ! Vertical canopy leaf-nitrogen profile. Per unit cumulative leaf area index
  ! from the top, n(L) = n0 exp(-kN_canopy L), with the sunlit and shaded class
  ! means obtained by integration and the total leaf nitrogen conserved
  ! (de Pury & Farquhar 1997, Plant Cell Environ 20:537-557).
  ! kN_canopy = 0.41 is the product of two factors measured in the SAME field
  ! wheat canopies by Moreau et al. (2012) Plant Physiol 160(3):1479-1490,
  ! doi:10.1104/pp.112.199935, Table II: the ratio K_N/K_L = 0.96 and the light
  ! extinction coefficient K_L = 0.425, the mean of that table two controlled-
  ! field experiment averages (0.42 and 0.43). 0.96 x 0.425 = 0.408, written
  ! 0.41 at the two-decimal precision of the source table. Table II spans
  ! 0.37-1.85 for the ratio across cultivars and experiments.
  ! ------------------------------------------------------------------------- !
  REAL(wp) :: kN_canopy = 0.41_wp
  ! ------------------------------------------------------------------------- !
  ! Canopy-nitrogen base (structural) term n_b - driver record 76.
  ! Record 76 = 0 distributes the whole leaf-nitrogen pool exponentially,
  !   n(L) = n0 exp(-kN L),  n0 = N_leaf kN / (1 - exp(-kN LAI)).
  ! Record 76 = 1 subtracts a non-photosynthetic base first, which is the form
  ! the measurement takes: Bertheloot, Martre & Andrieu (2008) Plant Physiol
  ! 148:1707-1720, doi:10.1104/pp.108.124156, Eq. (2) p.1718,
  !   SLN = (SLN0 - nb) exp(-K_N F) + nb,
  ! with F the cumulative green area index above the lamina and nb "SLN at
  ! which light-saturated photosynthesis is null". Moreau et al. (2012) Plant
  ! Physiol 160(3):1479-1490, Eq. 4, print the identical structure.
  ! nb_canopy = 0.40 g N m-2 lamina is the value the wheat source itself adopts
  ! (p.1718, "For wheat, the average of the values reported for nb is 0.4 g N
  ! m-2 (Araus and Tapia, 1987; Bindraban, 1999; Dreccer et al., 2000)"), and
  ! it sits at the center of the 0.25-0.55 g N m-2 wheat range reported by
  ! Moreau et al. (2012). It is not fitted here.
  ! nb_eff = MIN(nb_canopy, N_leaf/LAI) keeps the profile continuous when the
  ! canopy mean falls to the base, where Eq. (2) would otherwise invert.
  ! ------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: nb_canopy = 0.40_wp   ! g N m-2 leaf lamina; Bertheloot 2008 M&M (N-3, N-3a)
  ! ------------------------------------------------------------------------- !
  ! Runtime switch for the base term. The initializer is never the value in
  ! force: main_MMCO2Crop.f90 reads driver record 76 unconditionally and
  ! assigns 0 at end of file, so the deck always decides.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: canopy_n_base_profile = 0  ! runtime flag (driver_main.txt rec 76; v1.0 implements only 1 = structural base term added back uniformly; any other value is refused)
  ! ------------------------------------------------------------------------- !
  ! Canopy nitrogen profile: flat top zone, exponential below. Driver record 94
  ! selects it. The exponential of Eq. (2) was fitted BELOW the flag-leaf
  ! collar only, and the source measured no gradient above it.
  ! Measurement: Bertheloot, Martre & Andrieu (2008) Plant Physiol
  ! 148:1707-1720, doi:10.1104/pp.108.124156.
  !   p.1710, verbatim: "For both cultivars, no SLN gradient was observed in
  !     the upper part of the canopy (Fig. 4)."
  !   p.1710, the fit domain: the K_N to K_L ratio was calculated below the
  !     flag-leaf collar and averaged 0.92 +/- 0.05 (Fig. 5).
  !   p.1710: the light extinction coefficient K_L averaged 0.42 +/- 0.01
  !     (kL_bertheloot below).
  !   p.1718: Eq. (1) I = I0 exp(-K_L F), Eq. (2) above, and the nb value.
  ! Corroboration from the light side: Bertheloot et al. (2008) Funct Plant
  ! Biol 35:781-796, p.785 Eq. (17), PAR_i = k PAR_0 T exp(-k F_i), with T the
  ! "PAR transmission coefficient of the organs above the flag lamina"
  ! (Table 1 p.784).
  ! The implemented profile, with Fc the flag-rank green area index and
  ! Lb = LAI - Fc, E = (1 - exp(-kN Lb))/kN:
  !   n(L) = n_top                                            0  <  L <= Fc
  !   n(L) = nb_eff + (n_top - nb_eff) exp(-kN (L - Fc))      Fc <  L <= LAI
  !   n_top = (N_leaf - nb_eff (Lb - E)) / (Fc + E)
  ! and the sunlit-class mean follows from the beam-weighted integral
  !   Nsun = n_top (1-exp(-kb Fc))/kb
  !        + nb_eff (exp(-kb Fc) - exp(-kb LAI))/kb
  !        + (n_top - nb_eff) exp(-kb Fc) (1-exp(-(kN+kb) Lb))/(kN+kb).
  ! At nb_eff = 0 and Fc = 0 the form collapses to the single exponential.
  ! The effective PAR fraction above the collar follows from Eq. (1) in closed
  ! form, f_PAR_eff = ln((P + a)/(P exp(-kL LAI) + a)) / (kL LAI).
  ! The source measured the flat zone at anthesis and for about 350-450 degree-
  ! Cd after it; the model applies it for the whole season, which is an
  ! extension of the measurement beyond its observed window.
  ! ------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: kL_bertheloot = 0.42_wp    ! -, wheat canopy PAR extinction coefficient per unit cumulative green area index; Bertheloot, Martre & Andrieu 2008 Plant Physiol 148:1707-1720 p.1710, "averaging 0.42 +- 0.01" over two cultivars (not significantly different, P = 0.12). Used only by the record-94 f_PAR limb. The cultivar-specific k values of Bertheloot et al. 2008 Funct Plant Biol Table 1 p.784 (0.54/0.70/0.62) are excluded by the one-executable rule, and p.786 of that paper states why they differ from p.1710: they are based on the surface area index of laminae AND stems rather than on leaf area alone
  INTEGER  :: canopy_n_flag_zone = 0             ! runtime flag (driver_main.txt rec 94; EOF-default 0). 0 = uniform profile. 1 = flat-top profile, composing with the deck rec 76. 2 = flat-top with nb_eff forced to 0, which isolates the zone from the base term. 3 = mode 1 plus the layered f_PAR. 4 = the layered f_PAR only
  ! ------------------------------------------------------------------------- !
  ! Diagnostic sunlit / shaded leaf temperature - driver record 95. The energy
  ! balance solves ONE bulk canopy node; the released configuration reports that
  ! node. Record 95 = 1 adds a per-class diagnostic temperature obtained by
  ! re-evaluating the same leaf energy balance
  !   Rn_leaf = rho Cp (T - Tac)/rb + (rho Cp/gamma)(esat(T) - e_ac)/(rb + rst)
  ! with per-leaf resistances (rb = rca LAI, rst = rcs LAI) and the class
  ! radiation split of Teh, C. B. S. (2006) Introduction to Mathematical
  ! Modeling of Crop Growth, pp. 68-69:
  !   Eq. [3.38] Q_sl = alpha (k_dr Q_dr + Qbar_p,df + Q_p,dr,alpha)
  !   Eq. [3.39] Q_sh = alpha (Qbar_p,df + Q_p,dr,alpha)
  !   Eq. [3.41] L_sl = (1 - exp(-k_dr L))/k_dr,  Eq. [3.42] L_sh = L - L_sl,
  ! both fluxes being per unit LEAF area. Net longwave is divided per unit leaf
  ! area equally between the classes. The diagnostic is written only; it feeds
  ! no flux, and the bulk solve is unchanged.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: tleaf_class_diag = 0             ! runtime flag (driver_main.txt rec 95; EOF-default 0). 0 = OFF. 1 = ON: diagnostic sunlit/shaded leaf temperature + output/hourly_tleaf_classes.csv.
  ! ------------------------------------------------------------------------- !
  ! Shortwave canopy-soil multiple scattering - driver record 96. Record 96 = 1
  ! closes the canopy-soil exchange to infinite reflection order instead of the
  ! single scatter of the released path. With tau_sw the beam transmittance,
  ! alpha_can the canopy albedo and alpha_soil the soil albedo,
  !   rho_c = (1 - tau_sw) alpha_can,  a_c = (1 - tau_sw)(1 - alpha_can),
  !   rho_s = alpha_soil,
  !   D = tau_sw S + rho_c B,  B = rho_s D
  !   =>  D = tau_sw S / (1 - rho_c rho_s),  B = rho_s D
  !   sw_can_abs  = a_c (S + B),   sw_soil_abs = (1 - rho_s) D
  !   sw_ref_can  = rho_c S,       sw_ref_soil = tau_sw B.
  ! The closure is algebra on the two albedos, so it introduces no new
  ! empirical constant. Canopy and soil albedos are those of Monteith &
  ! Unsworth (2008) p. 48 Table 3.1 and the Karnieli et al. (1996) soil range.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: sw_radiosity_mode = 0 ! runtime flag (driver_main.txt rec 96; EOF-default 0). 0 = OFF, the transitional single-scatter shortwave. 1 = ON: exact two-layer shortwave radiosity to infinite reflection order.
  ! ------------------------------------------------------------------------- !
  ! Grain-nitrogen closure under the leaf holding capacity - driver record 97.
  ! Source: Bertheloot, Martre & Andrieu (2008) Funct Plant Biol 35:781-796,
  ! doi:10.1071/FP08064, p.785:
  !   Eq. (11)  dN_g,tot = min(S_g,tot_pot(t), N_c,mob(t)), with the verbatim
  !             statement that "the only source of N for the grains is the"
  !             culm-common mobile pool,
  !   Eq. (12)  S_g,tot_pot(t) = gamma N_g,tot_pot(t) for t <= tau, constant
  !             thereafter,
  !   Eq. (13)  dN_c,mob/dt = -SUM_i dN_i,ph/dt.
  ! Mode 1 makes the lamina holding capacity bind on the SAME pass as the
  ! synthesis term, s_lam = MIN(s_lam, MAX(0, sln_max_N lai_basis - N_leaf)), so
  ! that end-of-day N_leaf <= sln_max_N lai_basis on every path. sln_max_N =
  ! 3.8 g N m-2 is a holding capacity, the maximum measured field-wheat Narea of
  ! Silva-Perez et al. (2020) J Exp Bot 71(7):2299-2311,
  ! doi:10.1093/jxb/erz439, Fig. 4B p.2305 (genotype-mean maximum 3.85).
  ! Mode 2 closes the grain deficit against the same Eq. (12) sink,
  !   deficit = MAX(0, gnf_sink_day - MAX(0, N_mobile_veg)),
  ! the max(0, demand - supply) idiom already used by both neighbouring terms.
  ! No new constant is introduced by either mode.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: itemF_Nclose_mode = 0              ! runtime flag (driver_main.txt rec 97; EOF-default 0). 0 = OFF; 1 = the capacity respected by the same-day synthesis; 2 = the grain-nitrogen demand taken net of the mobile pool, with a terminal close; 3 = both
  REAL(wp) :: fnc_syn_clip_day = 0.0_wp        ! g N m-2 d-1, Eqn-6 lamina synthesis REFUSED by the capacity headroom under mode 1/3 (diagnostic; 0 on every other mode)
  REAL(wp) :: fnc_syn_clip_cum = 0.0_wp        ! g N m-2, season total of the above (diagnostic)
  REAL(wp) :: fnc_mob_close_day = 0.0_wp       ! g N m-2 d-1, residual mobile N deposited into N_stem by (diagnostic; 0 on every other mode)
  REAL(wp) :: fnc_mob_close_cum = 0.0_wp       ! g N m-2, season total of the above (diagnostic)
  ! ------------------------------------------------------------------------- !
  ! Leaf-nitrogen holding capacity, where it binds - driver record 132.
  ! Mode 1 retires the post-anthesis path that returns refused synthesis N
  ! straight back into the lamina pool, so the capacity is respected rather
  ! than exceeded. Sources for the direction:
  !   Bertheloot, Martre & Andrieu (2008) Funct Plant Biol 35:781-796,
  !     doi:10.1071/FP08064, p.789 and p.791 (post-anthesis lamina N declines;
  !     the mobile pool is culm-common and the lamina is not its sink),
  !   Bertheloot, Martre & Andrieu (2008) Plant Physiol 148:1707-1720,
  !     doi:10.1104/pp.108.124156, p.1711 (lamina N concentration falls to a
  !     structural minimum after flowering).
  ! Mode 2 is mode 1 with a printed ceiling in place of a figure reading:
  ! sln_max_N_lnh = 2.5 g N m-2, Pask, A. J. D. (2009) PhD thesis, University
  ! of Nottingham, p.29 (printed; "wheat leaves contain about 4 % N as a
  ! fraction of dry matter"). At the model's own lamina mass per area
  ! csl_LMA = 48.5 g DM m-2 the same 2.5 g N m-2 reads as 5.16 %DM rather than
  ! the 4 %DM of the source, and both numbers are reported. The diagnostic
  ! file output/daily_leafn_hold.csv is written only when the record is on.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: leafn_hold_law = 0                ! runtime flag (driver_main.txt rec 132; EOF-default 0). 0 = OFF · 1 = retire the post-anthesis FPB Eqn-6 lamina re-synthesis (stem/chaff draw stays in the mobile pool, available to the grain) · 2 = mode 1 + re-sourced holding ceiling sln_max_N_lnh (2.5 g N m-2, Pask 2009 p.29) in place of the pixel-read sln_max_N (3.8) in the general post-anthesis capacity ladder.
  REAL(wp), PARAMETER :: sln_max_N_lnh = 2.5_wp ! g N m-2 leaf, PRINTED (not pixel-read) ceiling: Pask (2009) PhD thesis p.29, "Wheat leaves commonly accumulate contents in excess of 2.5 g m-2 (about 4% N in dry matter)"; used only under leafn_hold_law = 2
  REAL(wp) :: lnh_cap_in_force   = 0.0_wp       ! g N m-2 leaf, the holding-capacity ceiling ACTUALLY applied today (sln_max_N or sln_max_N_lnh; diagnostic, mirrors the active_cap_N local used in growth_module.f90)
  REAL(wp) :: lnh_syn_assim_off_day = 0.0_wp        ! g N m-2 d-1, the synthesis-from-assimilate term the un-gated law would have taken (offered; diagnostic)
  REAL(wp) :: lnh_syn_mob_off_day   = 0.0_wp        ! g N m-2 d-1, the synthesis-from-mobile-pool term the un-gated law would have taken (offered; diagnostic)
  REAL(wp) :: lnh_redirect_day = 0.0_wp                  ! g N m-2 d-1, offered synthesis withheld under mode >= 1: it stays in the mobile pool instead of entering N_leaf (diagnostic)
  REAL(wp) :: lnh_redirect_cum = 0.0_wp         ! g N m-2, season total of the above (diagnostic)
  LOGICAL  :: lnh_binding_day = .FALSE.         ! TRUE on a day the retirement actually removed a positive synthesis claim (diagnostic)
  ! ------------------------------------------------------------------------- !
  ! State-driven lamina nitrogen turnover - driver record 133. The order of
  ! claims on the day new nitrogen becomes two staged passes:
  !   Stage A, the grain is served first, from assimilated N and the mobile
  !     pool, which is the source order of Bertheloot et al. (2008) Funct Plant
  !     Biol 35:781-796 p.785 ("a culm-common mobile pool" with no stem-to-
  !     lamina term), of Bertheloot et al. (2008) Plant Physiol 148:1707-1720
  !     p.1711, and of Pask (2009) pp.209-210 (stem reserve N remobilizes to
  !     the grain);
  !   Stage B, the lamina synthesis term takes only what Stage A left,
  !     syn_realized = MIN(s_lam, leftover), and the unserved remainder
  !     refused = MAX(0, s_lam - syn_realized) is refused rather than drawn,
  ! with the lamina degradation term of Eq. (5) (first order, 0.008 degree-Cd-1)
  ! unchanged. The ceiling is sln_max_N_lnh = 2.5 g N m-2 (Pask 2009 p.29).
  ! output/daily_leafn_turnover.csv is written only when the record is on.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: leafn_turnover_law = 0            ! runtime flag (driver_main.txt rec 133; EOF-default 0). 0 = OFF · 1 = STATE-DRIVEN two-stage order-of-claims rule (Stage A: grain claim resolved first against current assimilation + mobile pool, drawing stem/chaff for the grain only; Stage B: lamina re-synthesis draws ONLY the leftover of assimilation+mobile, no stem/chaff of its own, unmet remainder refused and reported) · 2 = mode 1 + re-sourced holding ceiling sln_max_N_lnh (2.5 g N m-2, Pask 2009 p.29, same PARAMETER as record 132) in the general post-anthesis capacity ladder. REQUIRES leafn_hold_law (rec 132) = 0 (ERROR STOP otherwise: the two are alternatives, never stacked).
  REAL(wp) :: lto_cap_in_force          = 0.0_wp  ! g N m-2 leaf, the holding-capacity ceiling ACTUALLY applied today (mirrors lnh_cap_in_force; diagnostic)
  REAL(wp) :: lto_sen_N_out_day         = 0.0_wp  ! g N m-2 d-1, = d_lam: the Eqn-5 lamina degradation flux OUT of N_leaf today (senescence export this record keeps UNCHANGED; diagnostic)
  REAL(wp) :: lto_sen_N_resorb_day      = 0.0_wp  ! g N m-2 d-1, = d_lam: Eqn-5 degradation is 100% resorbed into N_mobile_veg by construction in this limb (no dead-tissue split here; diagnostic)
  REAL(wp) :: lto_s_lam_pot_day         = 0.0_wp  ! g N m-2 d-1, the lamina's EXISTING Eqn-6 synthesis POTENTIAL today (= s_lam entering Stage B, already record-97-clipped if that record is ON; diagnostic)
  REAL(wp) :: lto_grain_claim_day       = 0.0_wp  ! g N m-2 d-1, = sink: the grain's EXISTING FPB Eqn-11/12 claim today, UNCHANGED by this record (diagnostic)
  REAL(wp) :: lto_grain_served_assim_day  = 0.0_wp ! g N m-2 d-1, Stage A: grain's claim served from TODAY'S assimilation (diagnostic)
  REAL(wp) :: lto_grain_served_mobile_day = 0.0_wp ! g N m-2 d-1, Stage A: grain's claim served from the MOBILE pool (diagnostic)
  REAL(wp) :: lto_stem_chaff_draw_day   = 0.0_wp  ! g N m-2 d-1, Stage A ONLY: stem+chaff remobilization drawn FOR THE GRAIN (never for the lamina under this record; diagnostic)
  REAL(wp) :: lto_mob_before_day        = 0.0_wp  ! g N m-2, N_mobile_veg on entry to this day's claim resolution (diagnostic)
  REAL(wp) :: lto_mob_after_day         = 0.0_wp  ! g N m-2, N_mobile_veg after Stage A + Stage B have both drawn on it (diagnostic)
  REAL(wp) :: lto_leftover_day          = 0.0_wp  ! g N m-2 d-1, Stage B: the leftover of (assim + mobile) OFFERED to the lamina after the grain's claim is met (diagnostic)
  REAL(wp) :: lto_syn_realized_day      = 0.0_wp  ! g N m-2 d-1, Stage B: the lamina re-synthesis ACTUALLY credited into N_leaf today, = MIN(s_lam, leftover) (diagnostic)
  REAL(wp) :: lto_syn_realized_cum      = 0.0_wp  ! g N m-2, season total of the above (diagnostic)
  REAL(wp) :: lto_syn_from_assim_day    = 0.0_wp  ! g N m-2 d-1, the assim-funded share of lto_syn_realized_day (diagnostic split only; one resolved number, not two draws)
  REAL(wp) :: lto_syn_from_mobile_day   = 0.0_wp  ! g N m-2 d-1, the mobile-funded share of lto_syn_realized_day (diagnostic split only)
  REAL(wp) :: lto_refused_day           = 0.0_wp  ! g N m-2 d-1, = MAX(0, s_lam - syn_realized): unmet lamina re-synthesis, REFUSED, never drawn from anywhere (diagnostic)
  REAL(wp) :: lto_refused_cum           = 0.0_wp  ! g N m-2, season total of the above (diagnostic)
  LOGICAL  :: lto_binding_day = .FALSE.          ! TRUE on a day the lamina's potential synthesis exceeded the leftover (refused > 0; diagnostic)
  ! ------------------------------------------------------------------------- !
  ! Spike structural sink and grain number - driver record 98. The spike
  ! demand over the spike-growth window, per axis and per day, is
  !   D_spike(t) = SUM_c ct_dens(c) bgn_r_spike MAX(0, T_air - 4.5)
  !   f_storage  = MIN(D_spike / dm_gain, 1 - f_leaf - f_root),
  ! and the floret number an axis can carry follows from its spike mass above
  ! the structural requirement.
  ! Measured basis, Fischer, R. A. (1985) J Agric Sci Camb 105:447-461,
  ! doi:10.1017/S0021859600056495: spike dry weight at anthesis 0.323 g per
  ! spike (Table 4 p.451, cv. Egret), 0.010 g spike DW per competent floret
  ! (p.449), and a 300 degree-Cd spike-growth period (p.449), giving the
  ! per-axis rate bgn_r_spike = 0.323/300 = 1.0767e-3 g axis-1 degree-Cd-1.
  ! The spike composition shares are measured by Rivera-Amado et al. (2019)
  ! Field Crops Research 240:154-167, doi:10.1016/j.fcr.2019.04.016, p.160
  ! Table 5 and Table 7 (lemmas 0.237, glumes 0.204), and the causal direction
  ! - greater partitioning to the florets raises fruiting efficiency - is
  ! stated by Slafer et al. (2015) Food Energy Secur 4(2):92-109, p.96 and
  ! Fig. 3. Rivera-Amado (2019) p.160 also measures the spike partitioning
  ! index as the conserved trait (0.212-0.264, CV 5.9 %).
  ! Derived masses, with no new constant:
  !   bgn_m_struct_spklt = f_struct spike_dw / spikelets
  !   bgn_m_floret_floral = f_floral spike_dw / florets
  !   bgn_w_struct_axis  = spikelets bgn_m_struct_spklt
  !   florets(axis) = MIN(developmental limit, available floral mass / mass per
  !                       floret).
  ! ------------------------------------------------------------------------- !
  INTEGER  :: spike_sink_mech = 0              ! runtime flag (driver_main.txt rec 98; EOF-default 0). 0 = OFF · 1 = BOTH pieces · 2 = P1 only (spike sink priority) · 3 = P2 only (spike-composition co-limitation).
  REAL(wp) :: bgn_demand_day   = 0.0_wp        ! g m-2 d-1, today's potential spike growth = axes * bgn_r_spike * degCd (diagnostic + the P1 claim)
  REAL(wp) :: bgn_axes_day     = 0.0_wp        ! axes m-2 alive when the demand was formed (diagnostic)
  REAL(wp) :: bgn_dem_cum      = 0.0_wp        ! g m-2, window total of bgn_demand_day (diagnostic)
  REAL(wp) :: bgn_got_cum      = 0.0_wp        ! g m-2, window total the spike actually received under P1 (diagnostic)
  INTEGER  :: bgn_days_dem     = 0             ! days the DEMAND was the binding limit (diagnostic)
  INTEGER  :: bgn_days_room    = 0             ! days the lamina/root ROOM cap was the binding limit (diagnostic)
  INTEGER  :: bgn_days_win     = 0             ! days the P1 split was applied at all (diagnostic)
  INTEGER  :: bgn_floral_bound = 0             ! axis classes whose florets were set by the FLORAL-MASS limit, not the developmental limit (diagnostic)
  REAL(wp) :: bgn_ear_spike_dm = 0.0_wp        ! g m-2, spike DM held by EAR-BEARING classes only -- the like-for-like basis for Rivera-Amado FE (diagnostic)
  REAL(wp) :: bgn_fe_ear_basis = 0.0_wp        ! grains g-1 of EAR-BEARING spike DM; the like-for-like comparator for Rivera-Amado Fig. 2(c) (diagnostic)
  REAL(wp) :: bgn_dm_close_max = 0.0_wp        ! g DM m-2 d-1, max |organ-state sum AFTER - organ-state sum BEFORE - dm_gain| across the partition, read from STATE on BOTH sides. GATE D's check for the rec-98 carbon path: no shipped closure covered the partition, and deleting the f_stem renormalization moves this by ~13 orders of magnitude
  REAL(wp) :: bgn_agb_day      = 0.0_wp        ! g m-2, above-ground organ state at the END of today's allocation (green+dead+stem+storage); the denominator of the at-anthesis SPI
  REAL(wp) :: bgn_spi_anth     = 0.0_wp                ! -, spike partitioning index at anthesis = spike DM / above-ground DM. Falsification target against the measured band 0.212-0.264 (Rivera-Amado et al. 2019 Field Crops Res 240:154, p.160 Table 5, CV 5.9 %), which is the conserved trait of that study

  REAL(wp) :: cls_gai_flag = 0.0_wp            ! m2 m-2, ground-area GREEN blade area index carried by the TOP (flag) rank of every axis class -- the model's own zone depth Fc (revision-record item E-5; internal, not distributed). Written in the S4 area loop beside cls_gai, from the same areas, so Fc <= cls_gai identically. Identically 0 when the S4 organ chain is OFF, which degenerates revision-record item E-4 EXACTLY to the incumbent profile
  REAL(wp) :: canE_n_top   = 0.0_wp            ! g N m-2, the flat top-zone Narea actually used (diagnostic; 0 when rec 94 = 0 or the profile did not run)
  REAL(wp) :: canE_fPAR_top = 0.0_wp           ! -, the grain-nitrogen sink's canopy-top f_PAR, recorded every day the Eqn-6 limb runs (diagnostic, both modes)
  REAL(wp) :: canE_fPAR_eff = 0.0_wp           ! -, the canopy-mean (revision-record item E-6; internal, not distributed) f_PAR, computed every day the Eqn-6 limb runs EVEN WHEN NOT USED, so the re-measured note-(iii) error is available on every arm (diagnostic)

  ! Ci solution - driver record 77. Mode 0 keeps the literal clamp of Ci to
  ! [0.20 Ca, 0.90 Ca]; mode 1 removes it and lets the bisection find Ci with
  ! the physical floor only (A = 0 at Ci = Gamma*_alphaG).
  INTEGER  :: ci_solve_mode = 0

  ! Supply/demand factor basis - driver record 78. f_t = Tr_act/Tr_pot. Mode 0
  ! takes Tr_pot from the current, already water-limited demand pass; mode 1
  ! takes it from the first, unstressed pass of the day, so that f_t measures
  ! the stress rather than latching on its own previous value.
  INTEGER  :: ft_pot_unstressed = 0

  ! Where soil-water stress enters the stomatal law - driver record 79. Mode 0
  ! multiplies the WHOLE Leuning conductance by f_t, which scales g0 as well,
  !   gsw = f_t (g0 + g1 An h(Ds)/(Ca - Gamma*));
  ! mode 1 applies it to the slope only,
  !   gsw = g0 + f_t g1 An h(Ds)/(Ca - Gamma*),
  ! the form Leuning (1995) Plant Cell Environ 18:339-355, p.350, describes.
  INTEGER  :: ws_two_loci = 0
  ! ------------------------------------------------------------------------- !
  ! Soil-plant-atmosphere continuum (SPAC) stomatal branch - driver record 81.
  ! Leaf water potential replaces the empirical soil-water factor on the slope:
  !   (L) gsw   = g0 + g1_tuzet f_psi(psi_L) An / (Ca - Gamma*)
  !   (F) f_psi = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f - psi_L)))
  !   (E) E     = gsw Ds / P                     leaf transpiration, mol m-2 s-1
  !   (H) psi_L = psi_s,rw - E / K_plant
  !   (S) psi_s,rw = rho g SUM_j w_j psi_m,j     root-weighted soil potential
  ! (F) and (H) are Duursma & Medlyn (2012) Geosci Model Dev 5:919-940, Eq. 3
  ! p.922 and Eq. 4 p.922 (E_L = k_L (psi_S - psi_L)); Ca replaces the Ci of
  ! Tuzet et al. (2003) as those authors recommend (p.923). There is no h(Ds)
  ! term: p.922 states that the psi_L modifier stands in place of f(D).
  ! (S) is Couvreur, Vanderborght & Javaux (2012) Hydrol Earth Syst Sci
  ! 16:2957-2971, Eq. 12 p.2961, with matric head only (the gravity term is at
  ! most 0.01 MPa over the profile). The weighting w_j is selected by the
  ! record: 1 = the root fraction itself (the literal Eq. 12 form), 2 = root
  ! fraction times the Feddes alpha of the committed pass (deliverability),
  ! 3 = root fraction times K_j/Ks_j, the rhizosphere-conductance form of
  ! Gardner (1960) Soil Sci 89:63-73.
  ! psi_L is the unique root of g(psi_L) = psi_s - E(psi_L)/K_plant, solved per
  ! leaf class and per hour by bisection inside the Ci solve and required to
  ! converge to 1e-6 MPa.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: spac_mode = 0                   ! rec 81: 0 OFF; 1/2/3 psi_s weighting (see above)
  ! ------------------------------------------------------------------------- !
  ! K_plant, the soil-to-leaf hydraulic conductance of Eq. 4 above.
  ! Duursma & Medlyn (2012) Geosci Model Dev 5:919-940 define k_L as the
  ! leaf-specific hydraulic conductance (symbol table p.935) and tabulate the
  ! plant component k_P at 2 mmol m-2 s-1 MPa-1 (Table 2 p.928) and 1.35
  ! (Table 3 p.930) - both for trees.
  ! The coded value is a WHEAT measurement on exactly the estimator Eq. 4 is:
  ! Koehler, Ding, Ossola, Anderson, Carminati & Sadok (2026) Plant Physiol
  ! 201(3):kiag500, doi:10.1093/plphys/kiag500, Eq. 1, K_sp = TR / dPsi_leaf,
  ! measured as a two-flux-state Darcy slope between predawn and midday,
  ! 0.40 mg m-2 s-1 kPa-1 = 0.40/18.015 x 1000 = 22.2 mmol m-2 s-1 MPa-1
  ! = 0.0222 mol m-2 s-1 MPa-1.
  ! A published wheat alternative, retained as a sensitivity value and not the
  ! default, is K_leaf = 8.70 mmol m-2 s-1 MPa-1 (Wang, Zhang, Luo, Hu, Feng,
  ! Mo, Li, Gong et al. 2024, Plants 13(16):2263, doi:10.3390/plants13162263),
  ! which is measured against midday leaf potential rather than the predawn-to-
  ! midday span, so it is not the same quantity.
  ! Standing limitation, stated because it is real: Eq. 4 wants the series
  ! soil-to-leaf conductance, 1/k_L = 1/k_P + 1/k_soil, while K_sp is measured
  ! at a soil potential above -0.0016 MPa, where the rhizosphere resistance is
  ! negligible. A constant K_plant therefore omits the rhizosphere term that
  ! grows as the soil dries.
  ! ------------------------------------------------------------------------- !
  REAL(wp) :: spac_K_plant = 0.0222_wp                   ! mol m-2 s-1 MPa-1 (rec 81 third token; default 0.0222)

  ! sf, the steepness of f_psi: 3.2 MPa-1, Tuzet et al. (2003) as tabulated by
  ! Duursma & Medlyn (2012) Geosci Model Dev 5:919-940, p.928 (their fitted
  ! cherry value is 2.2 MPa-1, Table 3 p.930). Not determined for wheat.
  REAL(wp), PARAMETER :: spac_sf = 3.2_wp                ! MPa-1
  ! ------------------------------------------------------------------------- !
  ! psi_f, the leaf water potential at which f_psi = 0.5 (definition: Duursma &
  ! Medlyn 2012, symbol table p.935). Default -1.5 MPa, assumed bracket
  ! -1.3 to -1.9 MPa. The bracket is an inference from measured wheat water
  ! relations rather than a fitted gs-psi_L curve: wheat turgor-loss point
  ! -1.85 MPa (Guizani et al. 2023, Front Plant Sci 14:1232583, Table 2, six
  ! genotypes), midday leaf potential -1.3 to -1.65 MPa under moderate and
  ! -1.8 to -2.25 MPa under severe stress (Abid et al. 2018, Sci Rep 8:4615,
  ! Table 2), and the generic Tuzet value -1.9 MPa (Duursma & Medlyn p.928).
  ! Wang et al. (2024) Plants 13(16):2263 report a conductance breakpoint at
  ! midday potential -1.48 MPa; that is a different quantity from a
  ! half-closure potential and is not adopted as if it were the same.
  ! ------------------------------------------------------------------------- !
  REAL(wp) :: spac_psi_f = -1.5_wp                       ! MPa (rec 81 second token; default -1.5)
  REAL(wp), PARAMETER :: spac_P_Pa = 101325.0_wp         ! Pa, the EB's own gsw -> rst pressure
  REAL(wp), PARAMETER :: spac_rho_g = 0.00981_wp         ! MPa m-1, matric head (m) -> MPa
  INTEGER,  PARAMETER :: spac_max_iter = 60              ! psi_L solve iteration cap (ERROR STOP beyond)
  REAL(wp), PARAMETER :: spac_tol_MPa = 1.0e-6_wp        ! |psi_L - (psi_s - E/K)| convergence / assertion
  ! Committed per-layer state (twins of h1_committed_m; captured at main_MMCO2Crop.f90's two commit sites;
  ! allocated in main once NSoil_layers is known). Filled in BOTH modes (assignment only; no flux reads
  ! them when OFF). k_committed_m_d is r_k of the same accepted pass (weighting 3 only).
  REAL(wp), ALLOCATABLE :: h_committed_m(:)      ! m, matric head per layer, accepted Richards state
  REAL(wp), ALLOCATABLE :: k_committed_m_d(:)    ! m d-1, K per layer, same vintage
  LOGICAL  :: h_committed_vec_valid = .FALSE.
  LOGICAL  :: spac_pass_on = .TRUE.              ! .FALSE. only on the outer-iteration-1 (unstressed) demand pass
  REAL(wp) :: spac_psi_s_MPa = 0.0_wp            ! psi_s,rw of the day under the selected weighting
  REAL(wp) :: spac_psi_s_conv(3) = 0.0_wp        ! the three conventions, every day (diagnostic)
  INTEGER  :: spac_psi_s_fallback = 0            ! 1 = selected weighting had zero weight sum -> Couvreur literal used
  INTEGER  :: spac_psi_s_valid = 0               ! 0 = no committed profile yet (day 1) -> psi_s = 0, f_psi = 1
  ! Hourly census (output/spac_census.csv; one row per assim call on the writer pass when rec 81 >= 1,
  ! classes 1 sunlit / 2 shaded / 3 flag leaf; header in both modes) and the daily solve counters.
  INTEGER  :: unit_spac_census = -1
  LOGICAL  :: spac_census_open = .FALSE.
  INTEGER  :: unit_spac_daily = -1
  INTEGER  :: spac_iters_max_day = 0
  REAL(wp) :: spac_resid_max_day = 0.0_wp
  INTEGER  :: spac_calls_day = 0
  REAL(wp) :: trans_pot_unstressed_mm = 0.0_wp   ! daily sum of the unstressed series (mm/d)
  REAL(wp) :: trans_pot_stressed_mm   = 0.0_wp   ! DIAGNOSTIC: pass-2 (stressed) potential (mm/d)
  ! ------------------------------------------------------------------------- !
  ! Feddes p2 node, units and interpolation - driver record 80. Mode 0 compares
  ! the demand handed to Richards in m d-1 with the threshold 0.5, which the
  ! source defines in cm d-1 on the potential transpiration, so the high node
  ! p2h can never be selected. Mode 1, the released value, converts the demand
  ! to cm d-1 and interpolates in the unit of the manual:
  !   h3 = p2h                                     Tp >= r2H
  !   h3 = p2h + (p2l - p2h)/(r2H - r2L) (r2H - Tp) r2L < Tp < r2H
  !   h3 = p2l                                     Tp <= r2L
  ! Source: Simunek, Sejna, Saito, Sakai & van Genuchten (2013) The HYDRUS-1D
  ! Software Package, Version 4.17, Table 12.7 "Block G - Root water uptake
  ! information", record 5a, p. 207, with r2H = 0.5 and r2L = 0.1 cm d-1.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: feddes_p2_hydrus = 0                        ! rec 80: EOF-default 0 = the un-interpolated node
  REAL(wp), PARAMETER :: feddes_r2h_cm_d = 0.5_wp   ! r2H, HYDRUS-1D 4.17 Table 12.7 p.207 (cm/day)
  REAL(wp), PARAMETER :: feddes_r2l_cm_d = 0.1_wp   ! r2L, HYDRUS-1D 4.17 Table 12.7 p.207 (cm/day)
  ! Diagnostics written in BOTH modes to output/daily_feddes_node_diag.csv (read by nothing):
  ! snapshot at the entry of the last Richards call of the day (same instant as
  ! hydro_alpha_diag). "alt" = what the OTHER rule would select at the same pressure heads.
  REAL(wp) :: hydro_tp_sent_cm_d_diag = 0.0_wp    ! demand handed to Richards at that call (cm/d)
  REAL(wp) :: hydro_tp_pot_cm_d_diag  = 0.0_wp    ! daily unstressed potential (cm/d)
  REAL(wp) :: hydro_p2_used_m_diag    = 0.0_wp    ! p2 node the active rule selected (m)
  REAL(wp) :: hydro_p2_alt_m_diag     = 0.0_wp    ! p2 node the other rule would select (m)
  INTEGER  :: hydro_p2_class_diag     = 0         ! active rule: 0 = p2l, 1 = interpolated, 2 = p2h
  INTEGER  :: hydro_p2_class_alt_diag = 0         ! other rule, same coding
  REAL(wp), ALLOCATABLE :: hydro_alpha_alt_diag(:)  ! Feddes alpha under the other rule, same h
  ! Daylight-only mean Ci (umol/mol) of the last canopy pass, accumulated under the SAME
  ! is_day mask ((Rnc+Rns) > 0) that gsw_day / Tr_daylight already use in CANOPY_MODEL.
  ! Written to daily_gas_exchange.csv as an APPENDED column only when rec 78 = 1.
  REAL(wp) :: ci_daylight_mean_diag = 0.0_wp
  ! Diagnostic-only census of where the old literals bind / would bind
  ! (output/ci_clamp_census.csv). Written by assim on the hourly-CSV writer pass only.
  ! NO gate and NO model quantity reads any of these; they are write-only bookkeeping.
  INTEGER  :: unit_ci_census = -1
  LOGICAL  :: ci_census_open = .FALSE.        ! newunit= yields NEGATIVE units, so test this, not the sign
  INTEGER  :: ci_census_tag  = 0               ! 0 = do not write; 1 sunlit, 2 shaded, 3 flag-leaf diag
  INTEGER  :: ci_census_day  = 0, ci_census_j = 0
  REAL(wp) :: ci_census_hour = 0.0_wp
  ! ------------------------------------------------------------------------- !
  ! Limiting-limb census - driver record 116. Diagnostic and write-only: it
  ! records which of the three photosynthetic limitations binds, and feeds no
  ! gate, flux or model quantity. The limb algebra it reports is that of
  ! Busch, Sage & Sage (2018) "Plants increase CO2 uptake by assimilating
  ! nitrogen via the photorespiratory pathway", Nature Plants 4:46-54,
  ! doi:10.1038/s41477-017-0065-x:
  !   Eq. (6)  p.48  Wc = Vcmax C/(C + Kc(1 + O/Ko))
  !   Eq. (11) p.48  A  = min{Wc,Wj,Wp}(1 - Gstar_aG/C) - Rd
  !   Eq. (12) p.48  Gstar_aG = 0.5(1 - aG) O / S_c/o
  !   Eq. (13) p.49  Wj = J/(4 + (4 + 8aG + 4aS) Phi)
  !   Eq. (14) p.49  Wp = 3Tp/(1 - 0.5(1 + 3aG + 4aS) Phi)
  !   Eq. (15) p.49  beta = 3aG_max/(3aG_max + 2aS_max)
  !   Eqs. (16)-(21) p.49 and SI pp.4-5, aG and aS under a Rubisco, a RuBP and
  !     a triose-phosphate-use limitation.
  ! The upstream primaries of that formulation are Farquhar, von Caemmerer &
  ! Berry (1980) Planta 149:78-90, Harley & Sharkey (1991) Photosynth Res
  ! 27:169-178 (the triose-phosphate limitation), Sharkey (1985) Plant Physiol
  ! 78:71-75 and von Caemmerer (2000).
  ! alpha_g_max = 0.09 and alpha_s_max = 0.38 are the two values printed in
  ! Busch et al. (2018) Fig. 3a p.52; the Tp/Vcmax ratio 0.11 used with them is
  ! the ratio of that same panel fitted set (11.5/103.9 = 0.1107). No wheat
  ! measurement of that ratio was located, and it is reported as such.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: limb_census_mode = 0   ! driver rec 116; 0 = OFF (EOF-default), 1 = ON
  INTEGER :: unit_flag_diag = -1   ! top/flag-leaf validation diagnostic (mode 2; never feeds totals)
  LOGICAL :: flag_diag_open = .FALSE.

  ! Integration
  REAL(wp) :: tol_h = 1.0E-5_wp ! Richards convergence tolerance for pressure head (m)

  !--------------------------------------------------------------------!
  ! 3. Global State Variables (Allocatable Arrays)
  !--------------------------------------------------------------------!
  ! Input Data Arrays
  REAL(wp), ALLOCATABLE :: array_input_data(:,:)
  INTEGER, ALLOCATABLE :: Year(:), Month(:), Day(:)
  REAL(wp), ALLOCATABLE :: s_hr(:), Tmax_d(:), Tmin_d(:), RHmax_d(:), RHmin_d(:)
  REAL(wp), ALLOCATABLE :: w_d(:), WS_d(:), Irr_d(:), Rad_d(:)

  ! Soil Parameters & State
  REAL(wp), ALLOCATABLE :: zb(:), tk(:), dz(:)
  REAL(wp), ALLOCATABLE :: h_m(:), hm(:), hg(:), k_x(:), k_s(:)
  ! VGM parameters for soil hydraulic functions (van Genuchten 1980; Mualem 1976).
  ! theta_r: residual water content; vg_alpha: 1/m; vg_n, vg_m: shape parameters; vg_l: pore-connectivity.
  REAL(wp), ALLOCATABLE :: theta_r(:), vg_alpha(:), vg_n(:), vg_m(:), vg_l(:)
  REAL(wp), ALLOCATABLE :: theta_1500(:), theta_33(:), theta_s33(:), theta_0(:)
  REAL(wp), ALLOCATABLE :: vwc_sum(:), wc(:), root_frac(:)
  ! Compensatory-root-uptake diagnostics (per-layer day snapshot; allocated in main after NSoil_layers is known).
  REAL(wp), ALLOCATABLE :: hydro_root_dist_diag(:), hydro_alpha_diag(:)
  REAL(wp), ALLOCATABLE :: inflx(:), outflx(:), net_flx(:)
  REAL(wp), ALLOCATABLE :: theta(:)  ! Volumetric water content (state)

  ! Large Scale Grid (if enabled)
  REAL(wp), ALLOCATABLE :: lat(:)

  !--------------------------------------------------------------------!
  ! 4. Simulation Control & Scalar State
  !--------------------------------------------------------------------!
  INTEGER :: NDays, DOY, NTsteps, NSoil_layers, NinVar, iD, I, J
INTEGER :: k_avg_mode = 1 ! 0=arithmetic, 1=geometric, 2=upstream
  INTEGER :: hyd_model = 0  ! 0=VGM (van Genuchten-Mualem), 1=Brooks-Corey, 2=Campbell (Weihermüller et al. 2021, p.3)
  CHARACTER(LEN=200) :: file_name, file_name_swtr

  ! Single Point Location
  REAL(wp) :: lon_s, lat_s, lat_s_rad, lon_s_rad

  ! Plant / Canopy State
  REAL(wp) :: lai_s, gParam_lai, gParam_lai_max
  REAL(wp) :: gParam_sla, plant_density, gParam_leafN_frac
  REAL(wp) :: gParam_hgt, gParam_hgt_s, gParam_refhgt
  REAL(wp) :: root_dep, depth_lim
  REAL(wp) :: root_growth_rate, N_max_param
  ! Dynamic N_max parameters (Busch 2017, Nature Plants 4:46-54, p. 48-49, Eqs. 16-21).
  ! N_max_ref = reference value for non-limiting NO3- supply; read from driver_N.csv.
  ! N_max_param is updated daily: N_max_param = N_max_ref * N_soil_min/(Km_N_busch + N_soil_min).
  REAL(wp) :: N_max_ref     = 1.21_wp  ! umol m-2 s-1 (Busch 2017, Fig. 3 caption)
  REAL(wp) :: Km_N_busch   = 0.10_wp  ! g N m-2 — CALIBRATION PARAMETER; no empirical value in literature
  ! CO2 effect on root NO3- uptake.
  ! Mechanism: Gojon et al. (2022, Trends Plant Sci, p. 4-5, Fig. 1C):
  ! eCO2 → higher Ci → lower Phi_o,c = Oa/(Ci*tau) → less 2-OG shuttled to
  ! cytosol → less NADH for cytosolic nitrate reductase (NR) → NRT transporter
  ! downregulation → "eCO2 significantly reduced root NO3- uptake by ~50%".
  ! Scaling uses the model's own iteratively computed daily mean Ci (Ci_day):
  ! f_CO2_N = Ci_ref_N / Ci_day
  ! At Ci_day = Ci_ref_N (ambient CO2): f = 1.0 (no reduction).
  ! At Ci_day = 2*Ci_ref_N (doubled CO2): f = 0.5 (50% reduction, Gojon 2022).
  ! Ci_ref_N = intercellular CO2 the crop achieves at your reference ambient Ca.
  ! Determine by running the model at ambient Ca and reading Ci_day, OR from
  ! measured leaf gas exchange (Ci at ambient conditions for this cultivar).
  ! No Ci/Ca ratio is assumed; Ci_day is solved by FvCB+Leuning at every timestep.
  ! ------------------------------------------------------------------------- !
  REAL(wp) :: Ci_ref_N     = 300.0_wp ! umol mol-1 — CALIBRATION PARAMETER (see above)
  REAL(wp) :: N_mineral_rate = 0.015_wp! g N m-2 d-1 (soil organic N mineralization)
  REAL(wp) :: feddes_p0_m = -0.10_wp
  REAL(wp) :: feddes_p1_m = -0.25_wp
  REAL(wp) :: feddes_p2h_m = -4.0_wp
  REAL(wp) :: feddes_p2l_m = -6.0_wp
  REAL(wp) :: feddes_p3_m = -160.0_wp
  INTEGER :: DVS, DVS_max

  ! Daily/Current Step Met Variables
  REAL(wp) :: s_hr_s, Tmax_d_s, Tmin_d_s, RH_max_s, RH_min_s, temp_s, RH_s
  REAL(wp) :: w_d_s, WS_d_s, Irr_d_s, pr, pr_eff, f_pr
  INTEGER :: irrig_mode = 2 ! 1=above canopy (intercepted), 2=below canopy (no interception)
  INTEGER :: current_day = 0 ! Current simulation day index (for diagnostics).
  LOGICAL :: hydro_failed = .FALSE. ! Set true if Richards fails on a given day.
  INTEGER :: hydro_dt_halvings = 0 ! Count of dt reductions within a day.

  ! Richards failure bookkeeping. Events count every solver-failure occurrence
  ! in a call (any outer iteration); the evaporation flag marks a failure that
  ! occurred under evaporative (flux_top <= 0) demand.
  INTEGER :: hydro_fail_events_call = 0
  LOGICAL :: hydro_fail_evap = .FALSE.
  INTEGER :: hydro_evap_fallback_events = 0  ! P2-b head-BC retries taken
  INTEGER :: hydro_fail_evap_days = 0        ! accepted-pass evap-failure days (gate class)
  LOGICAL :: hydro_min_dt_reached = .FALSE. ! True if min_dt is hit.
  REAL(wp) :: hydro_last_dt = 0.0_wp ! Last dt used when exiting hydro.
  REAL(wp) :: hydro_q34_m_d = 0.0_wp ! Flux between layers 3->4 (m/day, positive downward).
  REAL(wp) :: hydro_k34_m_d = 0.0_wp ! Interface K between layers 3 and 4 (m/day).

  ! ---- Wet-front Richards diagnostics (written only, read by no logic) ---- !
  LOGICAL  :: wetfront_trace_enabled = .FALSE. ! True only while the trace CSV is open.
  INTEGER  :: wetfront_trace_unit    = -1      ! Unit for output/hydro_wetfront_trace.csv (-1 = closed).

  ! Wet-front descriptors written to the trace; no solver logic reads them.
  REAL(wp) :: wetfront_Se_trigger      = 0.995_wp    ! stay engaged up to near-saturation; the time-step collapse occurs near Se = 0.99
  REAL(wp) :: wetfront_max_dtheta      = 0.010_wp
  REAL(wp) :: wetfront_flux_eps_m_d    = 1.0e-8_wp
  REAL(wp) :: wetfront_dt_limit_d      = 0.0_wp
  LOGICAL  :: wetfront_risk_flag       = .FALSE.
  LOGICAL  :: wetfront_limiter_applied = .FALSE.
  REAL(wp) :: dTotrad, dTotrad_s, Ca_s, Ci_last
  REAL(wp) :: T_leaf_mean_total, T_air_mean_total
  REAL(wp) :: T_leaf_day, T_air_day, VPD_day
  REAL(wp) :: T_soil_ref, T_soil_day, Tsoil_surf_s

  ! Fluxes & Energy Balance Variables
  REAL(wp) :: Rn_day, Rn_s, G_day, G_s, LE_day, LE_s, H_day
  REAL(wp) :: LEc_s, LEs_s, Hc_s, Hs_s
  REAL(wp) :: Ac_day, As_day
  REAL(wp) :: trans_act, evapo_act, ET_act
  REAL(wp) :: outer_relax = 0.5_wp ! Outer EB↔Richards relaxation (0-1)
  REAL(wp) :: evap_rate_act = 0.0_wp ! Actual Es rate used in Richards (m/d)
  REAL(wp) :: t_i, e_i
  REAL(wp) :: L_Temp_h, sky_emiss

  ! Sky emissivity of the component longwave path. Brutsaert (1975) clear-sky
  ! form as republished with units by Crawford & Duchon (1999) J Appl Meteorol
  ! 38:474-480, Eq. 16 p.476 (T in kelvin, vapor pressure in millibars), with
  ! their cloud modification (Eq. 11, cloud fraction from Eq. 3).
  REAL(wp) :: sky_emiss_comp = 0.85_wp

  ! Three-node energy balance (energy_balance_module.f90): the canopy, soil and
  ! canopy-airspace temperature nodes and the explicit airspace vapor node,
  ! populated once per hourly substep.
  REAL(wp) :: Tcan_s = 0.0_wp        ! canopy/leaf source temperature (C) [= T_leaf_s at solution]
  REAL(wp) :: Tac_s  = 0.0_wp        ! canopy-airspace aerodynamic temperature (C)
  REAL(wp) :: e_ac_s = 0.0_wp        ! canopy-airspace vapor pressure (Pa)
  REAL(wp) :: VPD_can_s = 0.0_wp     ! canopy-airspace VPD = esat(Tcan)-e_ac (Pa)
  REAL(wp) :: Ds_leaf_s = 0.0_wp     ! leaf-surface VPD passed to assim at convergence (Pa)
  REAL(wp) :: h_soil_s = 1.0_wp      ! soil-surface relative humidity (0-1)
  ! Solved (pre-hydrology-scaling) component fluxes for closure diagnostics (W m-2).
  REAL(wp) :: LEc_solved_s = 0.0_wp, LEs_solved_s = 0.0_wp
  REAL(wp) :: Hc_solved_s  = 0.0_wp, Hs_solved_s  = 0.0_wp
  REAL(wp) :: Rnc_s = 0.0_wp, Rns_s = 0.0_wp
  ! Residuals (W m-2) and convergence telemetry for the nested solve.
  REAL(wp) :: F_can_s = 0.0_wp, F_soil_s = 0.0_wp, F_air_s = 0.0_wp, F_vap_s = 0.0_wp
  INTEGER  :: eb3_inner_iters = 0    ! inner (Tsoil,Tac,e_ac) Newton iterations
  INTEGER  :: eb3_outer_iters = 0    ! nested assim<->EB fixed-point iterations
  INTEGER  :: eb3_conv_flag = 0      ! 1 if both inner Newton and fixed point converged
  INTEGER  :: unit_eb_diag = -1      ! file unit for hourly_energy_balance.csv (newunit-> negative)
  LOGICAL  :: eb_diag_open = .FALSE. ! .TRUE. once hourly_energy_balance.csv is opened
  INTEGER  :: unit_rad = -1                    ! file unit for hourly_radiation.csv
  LOGICAL  :: rad_open = .FALSE.     ! .TRUE. once hourly_radiation.csv is opened

  ! Soil-heat diagnostics; the units live here so the energy-balance module can
  ! write rows, and are opened and closed in main_MMCO2Crop.f90.
  INTEGER  :: unit_soilheat_h = -1   ! file unit for hourly_soil_heat.csv
  LOGICAL  :: soilheat_h_open = .FALSE.
  INTEGER  :: unit_soilheat_d = -1   ! file unit for daily_soil_heat.csv
  LOGICAL  :: soilheat_d_open = .FALSE.
  INTEGER  :: eb3_implausible_count = 0 ! daytime substeps with Tcan outside [Tair-15,Tair+25]

  ! Monin-Obukhov stability correction of the above-canopy heat resistance raa
  ! (Paulson 1970, Businger-Dyer functions). Diagnostics for inspection:
  REAL(wp) :: raa_neutral_s = 0.0_wp ! neutral (uncorrected) raa, s/m
  REAL(wp) :: L_obukhov_s   = 0.0_wp ! Obukhov length, m (large +/- = near neutral)
  REAL(wp) :: zeta_s        = 0.0_wp ! stability parameter (z_ref - d)/L
  REAL(wp) :: psih_s        = 0.0_wp ! integrated heat stability function psi_h
  REAL(wp) :: psim_s        = 0.0_wp ! integrated moment. stability function psi_m
  REAL(wp) :: ustar_s       = 0.0_wp ! friction velocity used for raa, m/s
  INTEGER  :: stab_regime_s = 0      ! 1=unstable (MO-corrected raa); 0=neutral/stable (raa=neutral)

  ! ---- Hydrology demand/supply handshake ---- !
  ! Per-substep POTENTIAL profiles stored during the demand canopy pass, so
  ! that the daily Richards actuals can be redistributed over the hours by
  ! potential weight and the three node temperatures repartitioned with the
  ! actual latent fluxes.
  INTEGER, PARAMETER :: EB_MAXSTEP = 24
  INTEGER  :: eb_nsteps_s = 0
  REAL(wp) :: ebh_weight(EB_MAXSTEP)        = 0.0_wp ! integration weight of each substep
  REAL(wp) :: ebh_LEc_trans_pot(EB_MAXSTEP) = 0.0_wp ! potential canopy transpiration LE (W m-2)
  REAL(wp) :: ebh_LEc_int_pot(EB_MAXSTEP)   = 0.0_wp ! potential wet-canopy (interception) LE (W m-2)
  REAL(wp) :: ebh_LEs_pot(EB_MAXSTEP)       = 0.0_wp ! potential soil evaporation LE (W m-2)
  REAL(wp) :: ebh_Rnc(EB_MAXSTEP)           = 0.0_wp ! canopy net radiation (W m-2)
  REAL(wp) :: ebh_Rns(EB_MAXSTEP)           = 0.0_wp ! soil net radiation (W m-2)
  REAL(wp) :: ebh_G(EB_MAXSTEP)             = 0.0_wp ! soil heat flux (W m-2)
  REAL(wp) :: ebh_rca(EB_MAXSTEP)           = 0.0_wp ! leaf boundary-layer resistance (s/m)
  REAL(wp) :: ebh_rsa(EB_MAXSTEP)           = 0.0_wp ! soil-airspace resistance (s/m)
  REAL(wp) :: ebh_raa_neutral(EB_MAXSTEP)   = 0.0_wp ! neutral above-canopy resistance (s/m)
  REAL(wp) :: ebh_Tair(EB_MAXSTEP)          = 0.0_wp ! air temperature (C)
  REAL(wp) :: ebh_wind(EB_MAXSTEP)          = 0.0_wp ! wind speed (m/s)
  REAL(wp) :: ebh_hour(EB_MAXSTEP)          = 0.0_wp ! solar hour
  REAL(wp) :: ebh_Tcan_pot(EB_MAXSTEP)      = 0.0_wp ! potential (demand) canopy temperature (C)
  REAL(wp) :: ebh_lai(EB_MAXSTEP)           = 0.0_wp ! LAI at the substep (for bare-soil guard)
  REAL(wp) :: ebh_hgt(EB_MAXSTEP)           = 0.0_wp ! canopy height at the substep
  REAL(wp) :: ebh_SW_in(EB_MAXSTEP)         = 0.0_wp ! incoming shortwave (W m-2)
  REAL(wp) :: ebh_LW_down(EB_MAXSTEP)       = 0.0_wp ! sky downwelling longwave (W m-2)
  REAL(wp) :: ebh_theta_top(EB_MAXSTEP)     = 0.0_wp ! top-layer water content (for soil albedo)

  ! Per-substep solar geometry, so that the actual-flux repartition replays
  ! each hour's own sun instead of the last substep values of Kdr_s and
  ! SolarElev_s.
  REAL(wp) :: ebh_Kdr(EB_MAXSTEP)           = 0.0_wp ! beam extinction Kdr at the substep
  REAL(wp) :: ebh_SolarElev(EB_MAXSTEP)     = 0.0_wp ! solar elevation (rad) at the substep

  ! Per-substep diagnostic sunlit / shaded leaf-temperature row (driver record
  ! 95). Filled in the demand pass at the converged state and written from the
  ! committed writer, so the committed bulk canopy temperature sits on the same
  ! row. Nothing here is read by any flux, state or solver.
  ! slots: 1 T_sun  2 T_shade  3 f_sun  4 Rn_leaf_sun  5 Rn_leaf_shade  6 rb_leaf
  !        7 rs_sun 8 rs_shade 9 T_bulk_check 10 Tac 11 e_ac 12 Tcan_demand
  !        13 LAI   14 valid(1/0)  15 canopy latent cap gap (W m-2)
  !        16 fcap = LEc_solved/LEc_potential, the latent basis the bulk node
  !           was solved on and therefore the one the class solves carry
  REAL(wp) :: ebh_tlc(16, EB_MAXSTEP) = 0.0_wp
  INTEGER  :: unit_tlc = -1                    ! output/hourly_tleaf_classes.csv (written only when record 95 is on)
  LOGICAL  :: tlc_open = .FALSE.     ! .TRUE. only when rec 95 >= 1 -- Gate A FILE-SET clause
  REAL(wp) :: T_leaf_sun_s = 0.0_wp, T_leaf_shd_s = 0.0_wp

  ! Committed-state reporting. The canopy pass buffers its hourly energy-balance
  ! row (the demand state) instead of writing it; the repartition captures the
  ! per-substep ACTUAL state; the hourly file and the daily ledgers are then
  ! written from the actual state, with the demand values retained as
  ! explicitly suffixed _demand columns.
  REAL(wp) :: ebh_row(43, EB_MAXSTEP) = 0.0_wp     ! buffered demand-state hourly row
  REAL(wp) :: eba_Tac(EB_MAXSTEP) = 0.0_wp, eba_Tcan(EB_MAXSTEP) = 0.0_wp
  REAL(wp) :: eba_Tsoil(EB_MAXSTEP) = 0.0_wp, eba_Rnc(EB_MAXSTEP) = 0.0_wp
  REAL(wp) :: eba_Rns(EB_MAXSTEP) = 0.0_wp, eba_LEc(EB_MAXSTEP) = 0.0_wp
  REAL(wp) :: eba_LEs(EB_MAXSTEP) = 0.0_wp, eba_Hc(EB_MAXSTEP) = 0.0_wp
  REAL(wp) :: eba_Hs(EB_MAXSTEP) = 0.0_wp, eba_G(EB_MAXSTEP) = 0.0_wp
  ! Daily ACTUAL energy sums (MJ m-2 d-1), accumulated in repartition_day_actual.
  REAL(wp) :: Rn_day_act = 0.0_wp, G_day_act = 0.0_wp
  REAL(wp) :: LEc_day_act = 0.0_wp, LEs_day_act = 0.0_wp
  REAL(wp) :: Hc_day_act = 0.0_wp, Hs_day_act = 0.0_wp

  ! Committed per-substep soil-node trajectory, produced by advancing a WORKING
  ! COPY of the heat profile inside every repartition (the real profile still
  ! advances only on the accepted iteration). The demand solve of the next outer
  ! iteration conducts against this trajectory instead of a day-frozen node.
  REAL(wp) :: eba_node1(EB_MAXSTEP) = 0.0_wp
  LOGICAL  :: eba_node1_valid = .FALSE.
  ! Per-substep override consumed by soil_heat_surface_coeffs (armed by the
  ! demand pass when a trajectory from the previous outer iteration exists;
  ! DISARMED inside the repartition, which manages its own sequential state).
  LOGICAL  :: shf_node1_use = .FALSE.
  REAL(wp) :: shf_node1_val = 0.0_wp

  ! Per-substep supply cap on soil latent heat in the DEMAND energy state, armed
  ! from outer iteration 2 with the previous iteration committed LEs(J). The
  ! network uncapped value is kept separately as the POTENTIAL that drives the
  ! hydrology demand and sets the repartition scaling basis.
  LOGICAL  :: eb_les_cap_use = .FALSE.
  REAL(wp) :: eb_les_cap_val = 0.0_wp
  REAL(wp) :: LEs_unc_last = 0.0_wp     ! network LEs before the cap, last eb_fluxes call
  REAL(wp) :: LEs_pot_unc_s = 0.0_wp    ! uncapped LEs at the solved state (POTENTIAL)

  ! Transpiration-side mirror of the soil cap: the DEMAND energy state canopy
  ! latent heat is capped at the previous outer iteration committed canopy
  ! latent (delivered root-sink energy plus interception), while the uncapped
  ! network value remains the POTENTIAL for the hydrology demand and the
  ! repartition scaling basis.
  LOGICAL  :: eb_lec_cap_use = .FALSE.
  REAL(wp) :: eb_lec_cap_val = 0.0_wp
  REAL(wp) :: LEc_unc_last = 0.0_wp
  REAL(wp) :: LEc_pot_unc_s = 0.0_wp

  ! Surface matric head of the ACCEPTED (committed) Richards state, captured
  ! once per day after the accepted hydrology pass. This is the vintage the
  ! demand-pass surface humidity reads.
  REAL(wp) :: h1_committed_m = 0.0_wp
  LOGICAL  :: h1_committed_valid = .FALSE.

  ! Per-substep soil-evaporation DEMAND series (mm per window), filled by the
  ! demand pass with the same mass conversion as the daily term, so the series
  ! sums to the daily demand exactly. The hydrology substep loops extract each
  ! call window demand from it by overlap weighting instead of forcing Richards
  ! with the daily mean.
  REAL(wp) :: ebh_esoil_dem_mm(EB_MAXSTEP) = 0.0_wp

  ! Per-substep TRANSPIRATION demand series (mm per window), built and converted
  ! like the soil series above, so it sums to the daily transpiration demand
  ! exactly. The hydrology substep loops extract each call window demand by
  ! overlap weighting, so the Feddes root sink carries the diurnal shape rather
  ! than the daily mean.
  REAL(wp) :: ebh_etr_dem_mm(EB_MAXSTEP) = 0.0_wp

  ! The UNSTRESSED transpiration-demand series (driver record 78): a copy of
  ! ebh_etr_dem_mm taken after the first outer iteration canopy pass of the day,
  ! when f_t = 1. At record 78 = 1 it replaces the stressed series before every
  ! later hydrology forcing.
  REAL(wp) :: ebh_etr_dem_unstr_mm(EB_MAXSTEP) = 0.0_wp

  ! Soil-evaporation MASS that Richards delivered, binned into the substep
  ! windows by overlap weighting, which conserves the daily Richards total to
  ! round-off. The committed soil latent heat at substep J is
  !   esoil_bin_mm(J) L(Tair_J) / dt_J,
  ! a flux the soil produced rather than a scaled potential shape. The ponded
  ! evaporation bin is reporting only: pond latent heat is not an energy-balance
  ! stream in this version.
  REAL(wp) :: esoil_bin_mm(EB_MAXSTEP) = 0.0_wp
  REAL(wp) :: pondevap_bin_mm(EB_MAXSTEP) = 0.0_wp
  LOGICAL  :: esoil_bin_valid = .FALSE.

  ! Sent-vintage snapshot of the demand series (mm per window, residue factor
  ! applied), captured when hydrology is forced; the second pass refills the
  ! demand series afterwards, so the pairing diagnostic needs this copy.
  REAL(wp) :: esoil_sent_mm(EB_MAXSTEP) = 0.0_wp

  ! Transpiration delivered/sent pairing, the mirror of the soil-evaporation
  ! pair: the committed canopy transpiration latent heat at substep J is the
  ! ROOT-SINK MASS Richards delivered in that window (binned, conserving the
  ! daily uptake to round-off) converted with L(Tair_J). The supply/demand ratio
  ! at the transpiration node is then a diagnostic.
  REAL(wp) :: etr_bin_mm(EB_MAXSTEP) = 0.0_wp
  REAL(wp) :: etr_sent_mm(EB_MAXSTEP) = 0.0_wp

  ! Mobile vegetative nitrogen pool (g N m-2): phloem-exported N resorbed from
  ! senescing laminae before the tissue is written off. Bertheloot, Martre &
  ! Andrieu (2008) Funct Plant Biol 35:781-796, doi:10.1071/FP08064, pp.782-783,
  ! where the pool is culm-common. It is not in the leaves (it is excluded from
  ! Narea, from Vcmax and from green-leaf N) and not yet in the grain; the grain
  ! draws on it through the existing source ordering, and it is available to
  ! pre-anthesis construction.
  REAL(wp) :: N_mobile_veg = 0.0_wp
  REAL(wp) :: N_from_mobile_day = 0.0_wp

  ! The soil node-1 value the demand solve was armed with per substep
  ! (-999 when unarmed, i.e. on the first outer iteration).
  REAL(wp) :: node_armed_row(EB_MAXSTEP) = -999.0_wp
  INTEGER  :: unit_node_pair = -1
  LOGICAL  :: node_pair_open = .FALSE.
  INTEGER  :: unit_esoil_pair = -1
  LOGICAL  :: esoil_pair_open = .FALSE.
  INTEGER  :: les_bin_exceed_count = 0      ! substeps with binned LEs > substep potential
  REAL(wp) :: les_bin_exceed_max = 0.0_wp   ! max ratio binned/potential at exceedance
  ! Committed-state radiation ledger per substep (captured after each
  ! repartition_eb_actual return) + buffered demand radiation rows.
  REAL(wp) :: rad_row(22, EB_MAXSTEP) = 0.0_wp     ! demand-state hourly_radiation row
  REAL(wp) :: eba_rad(10, EB_MAXSTEP) = 0.0_wp     ! committed LW/Rn set (see writer)
  REAL(wp) :: eb3_outer_resid_day = 0.0_wp          ! final |dTcan_act| outer residual

  ! Daily POTENTIAL totals from the demand pass (mm/d water-equivalent).
  REAL(wp) :: LEc_int_pot_day = 0.0_wp
  ! Interception: same-day intercepted water (mm), evaporated + drip-to-soil (mm).
  REAL(wp) :: Ir_day_mm = 0.0_wp, Eintercept_act_mm = 0.0_wp, interception_drip_mm = 0.0_wp
  ! Actual (post-Richards repartition) daily energy state for diagnostics (W m-2 / C).
  REAL(wp) :: Tcan_act_day = 0.0_wp, Tac_act_day = 0.0_wp, Tsoil_act_day = 0.0_wp
  REAL(wp) :: Hc_act_day = 0.0_wp, Hs_act_day = 0.0_wp
  REAL(wp) :: LEc_trans_act_day = 0.0_wp, LEc_int_act_day = 0.0_wp, LEs_act_day = 0.0_wp
  REAL(wp) :: eb_canopy_resid_day = 0.0_wp   ! mean |Rnc - Hc - LEc_trans_act - LEc_int_act|
  INTEGER  :: eb3_outer_iters_day = 0
  LOGICAL  :: eb3_outer_converged = .FALSE.

  ! ---- Component radiation ledger, set by radiation_module per call ---- !
  REAL(wp) :: rad_SW_down = 0.0_wp, rad_SW_abs_can = 0.0_wp, rad_SW_abs_soil = 0.0_wp
  REAL(wp) :: rad_SW_ref_can = 0.0_wp, rad_SW_ref_soil = 0.0_wp
  REAL(wp) :: rad_LW_down = 0.0_wp, rad_LW_emit_can = 0.0_wp, rad_LW_emit_soil = 0.0_wp
  REAL(wp) :: rad_LW_can_to_soil = 0.0_wp, rad_LW_soil_to_can = 0.0_wp
  REAL(wp) :: rad_LW_can_net = 0.0_wp, rad_LW_soil_net = 0.0_wp
  REAL(wp) :: rad_LW_up_sky = 0.0_wp   ! A2: model LWout to sky (radiosity diagnostic)

  ! Component radiation is unconditional in this version: the reported net
  ! radiation, the daily closure, the flux-feasibility and plausibility guards
  ! and the daylight gate all use the component fluxes Rnc + Rns.
  REAL(wp) :: rad_alpha_soil = 0.0_wp, rad_alpha_mix = 0.0_wp
  REAL(wp) :: rad_tau_sw = 0.0_wp, rad_tau_lw = 0.0_wp
  INTEGER  :: eb_repartition_iters = 0   ! Stage-4 fixed-point re-close iterations (repartition)
  REAL(wp) :: eb_repartition_resid = 0.0_wp ! max |dT| at re-close convergence (C)

  ! Optional free-convection fallback (Beljaars 1995) for very-low-wind or
  ! non-convergent Monin-Obukhov cases. Driver record 40 accepts only 0, so the
  ! fallback is not reachable in this version; eb_zi_m and eb_gust_beta would be
  ! explicit inputs, never tuning knobs, if it were enabled.
  INTEGER  :: eb_freeconv_fallback = 0      ! 0=off (MO only), 1=Beljaars fallback allowed
  REAL(wp) :: eb_zi_m       = 1000.0_wp     ! boundary-layer depth for w* (m) [fallback only]
  REAL(wp) :: eb_gust_beta  = 1.2_wp        ! free-convection gustiness coeff [fallback only]

  ! ---- Mechanistic 1-D soil heat (energy_balance_mode == 2) ---- !
  ! All driver-set scalars below are end-of-file defaulted from driver_main.txt.
  !
  ! Volumetric heat capacities of the soil constituents (J m-3 K-1) = specific
  ! heat x density, from de Vries, D. A. (1963) "Thermal properties of soils",
  ! in van Wijk, W. R. (ed.), Physics of Plant Environment, North-Holland,
  ! Amsterdam, pp. 210-235, Table 7.1:
  !   mineral 0.75 J g-1 K-1 x 2.66 g cm-3 = 2.0e6;  organic 1.9 x 1.3 = 2.5e6;
  !   water   4.18 x 1.00 = 4.18e6;        air 1.0 x 0.00125 = 1.25e3.
  ! Mixing law: Cv = SUM_j x_j (rho c)_j over the constituent volume fractions.
  REAL(wp) :: Cv_mineral = 2.00e6_wp   ! mineral solids (de Vries 1963, Table 7.1)
  REAL(wp) :: Cv_org     = 2.50e6_wp   ! organic solids (de Vries 1963, Table 7.1)
  REAL(wp) :: Cv_water   = 4.18e6_wp   ! liquid water (de Vries 1963, Table 7.1)
  REAL(wp) :: Cv_air     = 1.25e3_wp   ! soil air (de Vries 1963, Table 7.1)
  ! ------------------------------------------------------------------------- !
  ! Thermal-conductivity end-members (W m-1 K-1) of the geometric-mean solids
  ! model lambda_solid = lambda_quartz^q lambda_other^(1-q), with
  ! lambda_sat = lambda_solid^(1-phi) lambda_water^phi. Johansen, O. (1975)
  ! "Thermal conductivity of soils", PhD thesis, Univ. of Trondheim (US Army
  ! CRREL Draft Translation 637, 1977); the end-member values as used in the
  ! generalized form of Cote, J. & Konrad, J.-M. (2005) Can Geotech J
  ! 42(2):443-458, doi:10.1139/t04-106, whose normalized conductivity is
  ! lambda(theta) = lambda_dry + (lambda_sat - lambda_dry) Ke.
  ! ------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: lambda_quartz = 7.7_wp  ! quartz mineral (Johansen 1975)
  REAL(wp), PARAMETER :: lambda_other  = 2.0_wp  ! non-quartz minerals (Johansen 1975)
  REAL(wp), PARAMETER :: lambda_water  = 0.57_wp ! liquid water (Johansen 1975)
  ! Driver-set thermal INPUTS (EOF-defaulted; see plan Section 8 / INPUT_GUIDE.md):
  REAL(wp) :: quartz_fraction = 0.60_wp ! quartz fraction of solids -> lambda_solid; for sandy soils
                                        ! often approximated by sand fraction (Peters-Lidard, C.D.
                                        ! et al. 1998, J. Atmos. Sci. 55:1209-1224, p. 1212). Sandy default.
  INTEGER  :: bottom_bc_type  = 0       ! 0 = zero-flux bottom (first pass); 1 = prescribed deep T
  REAL(wp) :: Tdeep_C         = 15.0_wp ! prescribed deep-soil T (C); used only if bottom_bc_type==1
  INTEGER  :: soil_T_init_method = 1    ! 0 = uniform Tdeep_C; 1 = damped-exponential air->deep
  ! ------------------------------------------------------------------------- !
  ! Sensor-anchored initial soil-temperature profile, read from the optional
  ! deck file input/tsoil_init.csv (z_m,T_C anchors; "!" comments). Piecewise
  ! linear between anchors and held flat beyond the ends. When the file is
  ! absent the profile follows soil_T_init_method; when it is PRESENT but
  ! malformed the run stops (see the reader in main_MMCO2Crop.f90), so these
  ! arrays are either unset or hold tsoil_init_n >= 2 strictly increasing
  ! anchors -- never a partly-read file.
  ! ------------------------------------------------------------------------- !
  LOGICAL  :: tsoil_init_from_file = .FALSE.
  INTEGER  :: tsoil_init_n = 0
  REAL(wp) :: tsoil_init_z(8) = 0.0_wp
  REAL(wp) :: tsoil_init_T(8) = 0.0_wp
  ! ------------------------------------------------------------------------- !
  ! Deep heat domain: the heat grid extends below the 1.20 m hydrology column so
  ! that the seasonal heat reservoir exists; with a zero-flux boundary at 1.20 m
  ! the whole column equilibrates to the surface forcing within two to three
  ! weeks, which the measured deep winter offset contradicts. The depth is a
  ! numerics choice, shown converged by a 2 / 4 / 8 m test. Water content of the
  ! deep extension is held at the bottom hydrology-layer value, which is a
  ! stated assumption rather than a measurement. The bottom boundary is zero
  ! flux at the extended depth; heat_domain_target_m = 0 keeps the 1.20 m grid.
  ! ------------------------------------------------------------------------- !
  REAL(wp) :: heat_domain_target_m = 4.0_wp
  INTEGER  :: n_heat_tot = 0                   ! total heat layers (set at init)
  REAL(wp), ALLOCATABLE :: tk_heat(:)          ! heat-grid layer thicknesses (m)

  ! Prognostic profile and the thermal properties derived from it:
  REAL(wp), ALLOCATABLE :: Tsoil_prof(:)   ! per-layer soil temperature (C)
  REAL(wp), ALLOCATABLE :: Cv_soil(:)      ! per-layer volumetric heat capacity (J m-3 K-1; de Vries 1963)
  REAL(wp), ALLOCATABLE :: lambda_soil(:)  ! per-layer thermal conductivity (W m-1 K-1; Cote-Konrad 2005)

  ! Assimilation & Stomatal Conductance
  REAL(wp) :: A_day, An_day, A_s, An_s, N_day, N_s
  REAL(wp) :: A_can_sun, A_can_shade, An_can_sun, An_can_shade
  REAL(wp) :: Assim_net_total, Assim_total
  REAL(wp) :: Rd_total, Rd_s_can_sun, Rd_s_can_shade
  REAL(wp) :: N_total, N_can_sun, N_can_shade
  REAL(wp) :: gsw_An, gsw_mm_day_3
  REAL(wp) :: rsa_s, raa_s, rca_s, rss_s, rst_p, rst_a, rcs_p
  REAL(wp) :: f_assim_reduction, rwt, rt
  REAL(wp) :: Total_trans, Total_assim, Biomass, Total_N_ass, Total_gsw_1

  ! Hourly Integration Arrays (5-point Gauss)
  REAL(wp) :: const_ABS(5), const_WGT(5), dHour(5)

  ! Solar & Radiation Geometry
  REAL(wp) :: Solar_Decl, dTsr, Daylength, SolarElev_s, th_s
  REAL(wp) :: dfrad_s, drrad_s
  REAL(wp) :: I_drrad_h_s, I_dfrad_h_s, Absor_IPAR_sunlit_h_s, Absor_IPAR_shaded_h_s
  REAL(wp) :: sunlit_lai_s, shaded_lai_s, Kdr_s, par, par_day
  REAL(wp) :: HourSW_net_s, HourLW_net_s

  ! Soil Water / Hydro Temporaries
  REAL(wp) :: vwc_cr, vwc_pwp, vwc_sat, vwc_root, vwc_fc
  REAL(wp) :: wc_root, wc_sat, wc_fc, wc_pwp
REAL(wp) :: f_t
  LOGICAL :: use_soil_water_stress = .TRUE.
  LOGICAL :: allow_drain = .TRUE.
  INTEGER :: drainage_type = 0     ! 0 = Silled/No-Flux, 1 = Free Drainage (Gravity)
  LOGICAL :: use_explicit_infiltration = .FALSE. ! Richards handles atmospheric BC directly
  LOGICAL :: meteo_debug = .FALSE.
  INTEGER :: top_bc_type = 0   ! 0=flux, 1=head (ponding) at soil surface
  REAL(wp) :: h_top = 0.0_wp   ! surface pressure head (m) when top_bc_type=1
  REAL(wp) :: flux_top_actual = 0.0_wp ! actual top flux (m/day) from Richards
  REAL(wp) :: ponding_mm = 0.0_wp
  REAL(wp) :: ponding_max_mm = 5.0_wp ! max surface ponding depth (mm) for sandy soil
  REAL(wp) :: hydro_mb_resid_m = 0.0_wp
  REAL(wp) :: hydro_mb_correction_m = 0.0_wp
  REAL(wp) :: hydro_top_flux_m = 0.0_wp
  REAL(wp) :: hydro_top_infil_m = 0.0_wp
  REAL(wp) :: hydro_bottom_flux_m = 0.0_wp
  REAL(wp) :: hydro_root_sink_m = 0.0_wp
  ! Diagnostics for evaporation cap inside Richards
  REAL(wp) :: hydro_sum_avail_m = 0.0_wp
  REAL(wp) :: hydro_evap_cap_m_d = 0.0_wp
  REAL(wp) :: hydro_evap_rate_pot_m_d = 0.0_wp
  REAL(wp) :: hydro_evap_rate_day_m_d = 0.0_wp
  REAL(wp) :: hydro_flux_top_m_d = 0.0_wp
  REAL(wp) :: hydro_flux_top_req_m_d = 0.0_wp
  REAL(wp) :: hydro_flux_top_cap_m_d = 0.0_wp
  REAL(wp) :: hydro_flux_top_applied_m_d = 0.0_wp
  REAL(wp) :: hydro_infil_cap_m_d = 0.0_wp
  REAL(wp) :: hydro_k1_m_d = 0.0_wp
  REAL(wp) :: hydro_h1_m = 0.0_wp
  REAL(wp) :: hydro_theta1 = 0.0_wp
  REAL(wp) :: hydro_h_top_m = 0.0_wp
  REAL(wp) :: hydro_h_crit_atm_m = 0.0_wp

  ! Darcy-deliverable surface evaporation toward the Kelvin equilibrium head
  ! (m d-1, non-negative), evaluated at the start of each evaporative-BC substep.
  REAL(wp) :: hydro_qmax_evap_m_d = 0.0_wp
  INTEGER :: hydro_cap_active = 0
  REAL(wp) :: hydro_supply_step_mm = 0.0_wp
  REAL(wp) :: hydro_infil_step_mm = 0.0_wp
  REAL(wp) :: hydro_infil_from_supply_step_mm = 0.0_wp
  REAL(wp) :: hydro_infil_from_pond_step_mm = 0.0_wp
  REAL(wp) :: hydro_soil_evap_step_mm = 0.0_wp
  REAL(wp) :: hydro_pond_evap_step_mm = 0.0_wp
  REAL(wp) :: hydro_runoff_step_mm = 0.0_wp
  REAL(wp) :: hydro_ponding_start_step_mm = 0.0_wp
  REAL(wp) :: hydro_ponding_end_step_mm = 0.0_wp
  TYPE(HydroStepState) :: hydro_step
  TYPE(HydroCumLedger) :: hydro_ledger
  REAL(wp)    :: drain_day = 0.0_wp
  INTEGER :: water_stress_mode = 2
  INTEGER :: integration_mode = 1 ! 1=Gauss-daylight, 2=full 24-hr
  INTEGER :: physics_mode = 0  ! Internal selector only. Gates the single-source

  ! physics_mode is an internal selector of the component-longwave path in
  ! meteo.f90 at >= 1. It is frozen at 0 and is never set from
  ! energy_balance_mode.
  !
  ! energy_balance_mode selects the surface energy-balance architecture
  ! (driver_main.txt record 39; see INPUT_GUIDE.md). This version accepts the
  ! value 2 only - the three-node canopy / soil / airspace balance coupled to
  ! the 1-D soil-heat solution - and main_MMCO2Crop.f90 rejects every other
  ! value at the driver read.
  INTEGER :: energy_balance_mode = 0
  LOGICAL :: wt_pres

  ! Richards Solver Arrays (Allocatable to prevent re-allocation)
  REAL(wp), ALLOCATABLE :: r_a(:), r_b(:), r_c(:), r_r(:) ! Matrix diagonals and RHS
  REAL(wp), ALLOCATABLE :: r_h(:)      ! Matric potential head (m)
  REAL(wp), ALLOCATABLE :: r_k(:)      ! Hydraulic conductivity (m/d)
  REAL(wp), ALLOCATABLE :: r_cap(:)    ! Specific moisture capacity dTheta/dh (1/m)
  REAL(wp), ALLOCATABLE :: r_theta(:)  ! Volumetric water content at previous time step
  REAL(wp), ALLOCATABLE :: sink_source(:) ! Sink/Source term (1/d)

  ! Root dynamics. Jarvis-style compensatory root water uptake, with omega the
  ! critical uptake index in (0,1]:
  !   omega = 1      -> no compensation, sink_i = Tr_demand r_i alpha_i
  !   0 < omega < 1  -> sink_i = Tr_demand r_i alpha_i / MAX(omega, F),
  !                     F = SUM_i r_i alpha_i.
  ! Set from driver_main.txt, defaulting to 1.0 at end of file.
  REAL(wp) :: root_comp_omega = 1.0_wp
  REAL(wp), PARAMETER :: root_comp_eps = 1.0e-9_wp   ! F<=eps -> zero uptake, no divide-by-zero
  REAL(wp) :: hydro_Fraw_diag = 0.0_wp, hydro_comp_denom_diag = 1.0_wp  ! day snapshot for diagnostics

  ! Residue-cover suppression of soil evaporation, a static cover fraction;
  ! driver records 52-54. At mode 0 the suppression factor is exactly 1.0.
  INTEGER  :: residue_evap_mode = 0
  REAL(wp) :: residue_Am_ha_t = 0.45_wp
  REAL(wp) :: residue_fsupp = 0.75_wp
  REAL(wp) :: f_res_evap = 1.0_wp
  REAL(wp) :: M_res_surf_t_ha = 0.0_wp
  ! ------------------------------------------------------------------------- !
  ! Carbon-reserve remobilization during grain filling; driver records 55-60.
  ! Pool size basis: Ruuska, Rebetzke, van Herwaarden, Richards, Fettell,
  ! Tabe & Jenkins (2006) Funct Plant Biol 33:799-809, doi:10.1071/FP06062.
  ! Their anthesis values of 108-203 mg g-1 (environments, Table 1 p.802) and
  ! 112-213 mg g-1 (genotypes, Table 2 p.802) are ABOVE-GROUND GRAB SAMPLES,
  ! not stem samples; the paper separates the two itself, p.803, where grab
  ! samples average 183 +/- 8 against 297 +/- 13 mg g-1 DW for stem samples.
  ! The stem water-soluble-carbohydrate fraction of stem dry mass is therefore
  ! read from the main-stem series of their Fig. 1 p.804. The same paper
  ! records the direction of the two interactions this model needs: the genetic
  ! correlation of reserve concentration with nitrogen concentration is
  ! rg = -0.80 +/- 0.11 (p.804), and reserve accumulation "is known to be
  ! increased by high CO2 (Dubois et al. 1990; Winzeler et al. 1990) and
  ! decreased by shading" (p.807).
  ! ------------------------------------------------------------------------- !
  INTEGER  :: carbon_reserve_mode = 0    ! rec 55: 0=OFF, 1=ON
  REAL(wp) :: stem_reserve_frac  = 0.20_wp ! rec 56: fraction of stem DM increment banked pre-anthesis (-)
  ! ------------------------------------------------------------------------- !
  ! reserve_remob_rate: the first-order draw rate of the reserve pool (d-1).
  ! No published wheat value for this rate was located. The two possible
  ! sources were read and do not supply one: Gebbing & Schnyder (1999) Plant
  ! Physiol 121:871-878, doi:10.1104/pp.121.3.871, measures the FRACTION of
  ! pre-anthesis carbon that reaches the grain (p.874, "on average, 83 % of the
  ! pre-anthesis C in grain carbohydrates was present as starch"; p.877, the
  ! lower bound >= 0.83) but prints no rate; Gebbing, Schnyder & Kuhbauch
  ! (1998) Plant Cell Environ 21:301-313 sampled only at anthesis and at
  ! maturity (p.303), so no rate can be read from it either. The value below is
  ! therefore an assumed first-order rate, not a measured constant.
  ! ------------------------------------------------------------------------- !
  REAL(wp) :: reserve_remob_rate = 0.05_wp ! first-order draw rate (d-1)
  REAL(wp) :: reserve_start_dvs2 = 1.00_wp ! rec 58: window start (dvs2; 1=anthesis, model-internal scale)
  REAL(wp) :: reserve_end_dvs2   = 2.00_wp ! rec 59: window end (dvs2; 2=maturity, model-internal scale)
  ! ------------------------------------------------------------------------- !
  ! reserve_conv_eff: grain dry matter per unit reserve dry matter. Gebbing &
  ! Schnyder (1999) Plant Physiol 121:871-878, doi:10.1104/pp.121.3.871, whose
  ! measured pre-anthesis carbon recovery in grain carbohydrate places 0.95
  ! within one standard error.
  ! ------------------------------------------------------------------------- !
  REAL(wp) :: reserve_conv_eff   = 0.95_wp ! g grain DM per g reserve DM
  ! ------------------------------------------------------------------------- !
  ! Shared stem-to-spike reallocation fraction applied over the spike-growth
  ! window (driver record 61); 0 leaves the partition unchanged.
  ! ------------------------------------------------------------------------- !
  REAL(wp) :: stem_spike_realloc_delta = 0.0_wp ! rec 61: fraction of f_stem moved to f_storage in the spike window
  REAL(wp) :: W_reserve = 0.0_wp         ! stem WSC sub-pool (g DM m-2), a labeled part of W_stem
  ! ------------------------------------------------------------------------- !
  ! Grain-nitrogen import termination state. Martre, Porter, Jamieson & Triboi
  ! (2003) Plant Physiol 133:1959-1967, doi:10.1104/pp.103.030585: grain
  ! nitrogen import stops before the end of dry-matter accumulation.
  ! ------------------------------------------------------------------------- !
  LOGICAL  :: grain_fill_ended = .FALSE.
  LOGICAL  :: canopy_senesced  = .FALSE.
  INTEGER  :: day_fill_ended = 0, day_canopy_dead = 0
  REAL(wp) :: assimN_surplus_postfill = 0.0_wp   ! N kept vegetative after termination
  ! B-1: post-maturity / post-fill carbon bookkeeping (mass-conserved, labeled).
  REAL(wp) :: unalloc_C_total = 0.0_wp           ! g CH2O m-2, dvs2>=2 unallocated
  REAL(wp) :: overflow_to_reserve_total = 0.0_wp ! g DM m-2, sink-full surplus -> stem WSC

  ! Spike-nitrogen split. N_storage remains the TOTAL spike nitrogen pool and
  ! every existing flow is unchanged; N_chaff tracks its chaff-structural share
  ! so that grain and chaff nitrogen can be reported separately.
  REAL(wp) :: N_chaff = 0.0_wp
  ! daily diagnostics for output/daily_carbon_reserve.csv (mode 1 only)
  REAL(wp) :: cres_to_grain_assim = 0.0_wp, cres_to_grain_res = 0.0_wp
  REAL(wp) :: cres_sink_demand = 0.0_wp, cres_sink_unfilled = 0.0_wp
  REAL(wp) :: cres_respired = 0.0_wp, cres_ledger_resid = 0.0_wp

  ! Soil-heat robustness guard: on a non-finite soil-heat residual the profile
  ! is HELD for that substep and the event is counted and reported.
  INTEGER :: soil_heat_fb_events = 0, soil_heat_fb_days = 0, soil_heat_fb_lastday = -1

  ! Canopy energy-balance solver guard: when the hourly canopy-temperature
  ! root-find has no sign change in its bracket, or a residual is not finite,
  ! the hour falls back to a flagged evaluation at air temperature and the event
  ! is counted and reported.
  INTEGER :: eb_noroot_events = 0, eb_noroot_days = 0, eb_noroot_lastday = -1
  INTEGER :: eb_seed_resets = 0

  ! Dry-surface head boundary-condition release. The atmospheric dry-head
  ! boundary condition (h_top = h_crit_atm; the atmospheric boundary condition
  ! of HYDRUS-1D, Simunek et al. 2013) is an evaporation limiter, so it is
  ! released as soon as the surface can again take the flux. Events are counted
  ! and reported.
  INTEGER :: hydro_bcswitch_events = 0, hydro_bcswitch_days = 0, hydro_bcswitch_lastday = -1

  ! Flux-side well-posedness fallback: on an hour whose accepted energy state
  ! carries a physically impossible latent-heat demand the hour reverts to
  ! equilibrium evaporation, LE_eq = delta/(delta+gamma) MAX(0, Rn) (Slatyer &
  ! McIlroy 1961), and the event is counted and reported rather than passing
  ! silently. The 600 W m-2 detection threshold is an assumed numerical
  ! allowance for advective enhancement; no published value was located for it.
  INTEGER :: eb_fluxfb_events = 0, eb_fluxfb_days = 0, eb_fluxfb_lastday = -1
  REAL(wp) :: root_beta = 0.96_wp      ! Jackson et al. (1996) parameter (0.9-0.98)
  REAL(wp) :: root_max_depth = 1.5_wp  ! Max root depth (m); runtime-set from driver_swtr.txt

  ! Miscellaneous working scalars.
  REAL(wp) :: n, a, b, ratio
  REAL(wp) :: dA, dB
  REAL(wp) :: temp_offset = 0.0_wp
  REAL(wp) :: T_leaf_s
  REAL(wp) :: HourTemp_s, HourNetRad_s, HourVP_s, Svp_s, VpdRef_s
  REAL(wp) :: ac_s, as_s, dLe_s
  REAL(wp) :: phi, J_e
  REAL(wp) :: gParam_rstA1, gParam_rstA2, gParam_leafwidth
  REAL(wp) :: a_wind, n_eddy, dry_layer_thick
  INTEGER :: use_dynamic_n_eddy
  REAL(wp) :: gParam_porosity, gParam_poredist, gParam_vwcsat, gParam_vwc01

  ! Plant biomass pools (g dry matter m-2 ground area)
  REAL(wp) :: W_green, W_dead, W_stem, W_root, W_storage
  ! Yield-component sub-pools: W_storage (reproductive) = W_grain + W_chaff
  REAL(wp) :: W_grain = 0.0_wp        ! grain dry mass (g m-2)
  REAL(wp) :: W_chaff = 0.0_wp        ! ear/chaff (non-grain spike) dry mass (g m-2)
  REAL(wp) :: kernel_number_m2 = 0.0_wp  ! grains per m2 (set at ear_end_dvs2)
  LOGICAL  :: grain_number_set = .FALSE.
  ! Cultivar yield-component coefficients (read from input/cultivar_yield.csv) — cv. Ruta defaults
  REAL(wp) :: fruiting_efficiency = 90.0_wp  ! grains per g ear DW (Slafer 2015; Mulero SI Ruta)
  REAL(wp) :: TKW_potential       = 31.0_wp  ! g per 1000 kernels (potential; Mulero SI Ruta ~31)

  ! Spike-window partition trait nodes, read per cultivar from
  ! input/cultivar_yield.csv so that a cultivar never needs a rebuild. The
  ! defaults are cv. Ruta: the onset 0.693 is the development stage 30 days
  ! before anthesis (Fischer 1985 J Agric Sci Camb 105:447, p.447), and the
  ! peak 0.474 is twice the measured chaff mass divided by the gross dry-matter
  ! gain over the window (Slafer et al. 2015 Food Energy Secur 4(2):92-109,
  ! p.96).
  REAL(wp) :: storage_onset_dvs2  = 0.693_wp
  REAL(wp) :: storage_peak_frac   = 0.474_wp
  REAL(wp) :: chaff_share         = 0.55_wp  ! fraction of ear-window reproductive DM to chaff (calibratable)
  REAL(wp) :: ear_end_dvs2        = 1.35_wp  ! dvs2 end of ear construction / grain-number set
  ! T2 leaf-partition window (cultivar/site; input/cultivar_partition.csv; defaults = WOFOST baseline)
  REAL(wp) :: fleaf_elong   = 0.2741_wp ! f_leaf AFGEN node at dvs2=0.646 (WOFOST Winter_wheat_106 conversion)
  REAL(wp) :: fleaf_heading = 0.0_wp    ! f_leaf AFGEN node at dvs2=0.95 (WOFOST: window closed at heading)
  REAL(wp) :: fleaf_end_dvs = 0.95_wp   ! dvs2 where f_leaf reaches 0 (>0.95 adds a closure node)
  REAL(wp) :: N_leaf
  REAL(wp) :: N_stem, N_root, N_storage, N_dead_leaf
  REAL(wp) :: N_assim_biochem_day

  ! Soil mineral N pool — state and parameters read from input/driver_N.csv
  ! and fertilizer schedule from input/N_fertilizer.csv.
  ! All values below are defaults; actual values are overridden by file read in main.
  REAL(wp) :: N_soil_min = 0.0_wp          ! current soil mineral NO3-N (g N m-2); == N_NO3_soil (leachable, taken up)
  ! Fertilizer availability pools (Task: 3-pool urea/NH4/NO3 model). N_soil_min IS the NO3 pool.
  REAL(wp) :: N_NH4_soil  = 0.0_wp         ! ammonium-N (g N m-2): retained (not leached), nitrifies to NO3
  REAL(wp) :: N_urea_soil = 0.0_wp         ! urea-N (g N m-2): hydrolyzes to NH4

  ! ------------------------------------------------------------------ !
  ! Mechanistic ammonia volatilization and denitrification. The surface
  ! ammoniacal pool, its pH, its partition coefficient and its bulk density are
  ! deck inputs read from input/driver_N.csv.
  REAL(wp) :: TAN_surf    = 0.0_wp  ! surface (0-z_surf_N) ammoniacal N pool, g N m-2:
                                    ! fed by surface urea hydrolysis + NH4-form fertilizer;
                                    ! drained by NH3 volatilization, nitrification, and
                                    ! infiltration advection to the bulk NH4 pool.
  REAL(wp) :: pH_surf     = 6.5_wp  ! measured surface-layer soil pH (AgMIP workbook
                                    ! "Soil characteristics" sheet, per site; driver_N pos 42)
  REAL(wp) :: Kd_NH4      = 3.0_wp  ! NH4+ solid/solution partition coefficient, L kg-1
                                    ! (measured range 1-4 for loams; Venterea et al. 2015,
                                    ! Methods of Soil Analysis; sensitivity-bounded, pos 43)
  REAL(wp) :: BD_surf     = 1.4_wp  ! surface-layer bulk density, Mg m-3 (workbook, pos 44)
  REAL(wp) :: z_surf_N    = 0.05_wp ! emission-layer thickness, m (measured emission depth
                                    ! 0-5 cm: Sommer et al. 2004, Adv Agron 82:557, sec III.B; pos 45)
  REAL(wp) :: Km_den      = 0.035_wp ! Michaelis constant for denitrification in NO3, g N L-1
                                    ! (measured 22-53 mg N L-1: Murray et al. 1989; pos 46)

  ! Nitrogen-uptake TIMING as one global process: component-bounded
  ! pre-anthesis demand plus grain-demand-driven post-anthesis uptake. All the
  ! constants below are global and taken from the literature; no site factors.
  REAL(wp), PARAMETER :: n4_sln_target = 1.5_wp   ! gN m-2 leaf area
      ! Jamieson & Semenov 2000, FCR 68:21, p.23: "Leaf is assumed to
      ! contain 1.5 g N m-2" (constant [N]a).
  REAL(wp), PARAMETER :: n4_stem_nmax = 0.015_wp  ! gN/gDM stem ceiling
      ! Jamieson & Semenov 2000 p.23: stem "may have a concentration of
      ! up to 1.5% N to store excess N" (min 0.5% while growing).
  REAL(wp), PARAMETER :: n4_nhi_target = 0.80_wp  ! target N harvest index
      ! Jamieson & Semenov 2000 p.24, hypothesis H1 (source-determined
      ! thermal-time trajectory to NHI 80%; beat sink-determined H2).

  ! Prospective increments (the previous day realized growth, an explicit
  ! one-day lag) and the stem structural minimum used by the priority
  ! allocation.
  REAL(wp), PARAMETER :: n5_stem_nstruct = 0.005_wp ! gN/gDM stem minimum
      ! Jamieson & Semenov 2000 p.23: stem "minimum requirement of 0.5%
      ! N while the plant is actively growing".
  REAL(wp) :: n5_lai_prev     = 0.0_wp
  REAL(wp) :: n5_dlai_prev    = 0.0_wp
  REAL(wp) :: n5_dwstem_prev  = 0.0_wp
  REAL(wp) :: n5_dwchaff_prev = 0.0_wp
  REAL(wp) :: n5_dwroot_prev  = 0.0_wp
  LOGICAL  :: n4_anth_snap = .FALSE.  ! anthesis crop-N snapshot taken?
  REAL(wp) :: n4_ag_anth   = 0.0_wp   ! above-ground crop N at anthesis
  REAL(wp) :: n4_grain_deficit_prev = 0.0_wp ! yesterday's unmet grain-N
      ! demand -> today's root uptake demand (explicit 1-day lag)
  REAL(wp) :: veg_N_demand_day = 0.0_wp        ! diagnostic (§5.4)
  REAL(wp) :: grain_N_demand_postanth_day = 0.0_wp
  REAL(wp) :: N_uptake_pre_cum  = 0.0_wp
  REAL(wp) :: N_uptake_post_cum = 0.0_wp

  ! Demand-driven grain-nitrogen remobilization state. The snapshots are frozen
  ! at the start of the source-limited filling phase (dvs2 = dvs_grain_stru_end).
  LOGICAL  :: n3_fill_snap = .FALSE. ! snapshot taken this season?
  REAL(wp) :: n3_dvs0      = 0.0_wp  ! dvs2 at snapshot (fill start)
  REAL(wp) :: n3_leaf0     = 0.0_wp  ! leaf N at fill start (g N m-2)
  REAL(wp) :: n3_stem0     = 0.0_wp  ! stem N at fill start (g N m-2)
  REAL(wp) :: n3_cum_leaf  = 0.0_wp  ! cumulative leaf N exported to grain
  REAL(wp) :: n3_cum_stem  = 0.0_wp  ! cumulative stem N exported to grain

  ! Measured per-organ remobilization efficiency ceilings (the fraction of the
  ! organ fill-start nitrogen exportable by maturity). Pask, Sylvester-Bradley,
  ! Jamieson & Foulkes (2012) Field Crops Res 126:104-118,
  ! doi:10.1016/j.fcr.2011.09.021: leaf laminae 0.76, true stem 0.48, leaf
  ! sheath 0.61, chaff 0.56. The model leaf pool is the laminae, so it takes
  ! 0.76; the model stem pool lumps true stem and sheath, and takes the
  ! conservative true-stem value 0.48.
  REAL(wp), PARAMETER :: n3_eff_leaf = 0.76_wp
  REAL(wp), PARAMETER :: n3_eff_stem = 0.48_wp

  ! ------------------------------------------------------------------ !
  ! Living-lamina nitrogen drawdown during grain filling and the leaf-nitrogen
  ! inflow capacity bound - driver record 62, which this version accepts at the
  ! value 0 only, so the incumbent lamina nitrogen path runs. The constants it
  ! declares are used by the records that are active, and their sources are:
  ! - k_leafN_TT = 0.0035 (degree-Cd)-1, the relative rate of green-lamina
  !   nitrogen depletion. Bertheloot, Martre & Andrieu (2008) Plant Physiol
  !   148:1707-1720, doi:10.1104/pp.108.124156, Eq. (5) with the fitted k
  !   reported p.1712 ("A k value of 0.0035 [degree-Cd]-1 was found to minimize
  !   the root mean square error"), a single k across laminae, on thermal time
  !   with base 0 C as their Materials and Methods state (p.1718).
  ! - sln_max_N = 3.8 g N m-2 leaf, the leaf nitrogen holding capacity: the
  !   maximum measured field-wheat Narea of Silva-Perez, Furbank, Condon &
  !   Evans (2020) "Genetic variation for photosynthetic capacity and efficiency
  !   in spring wheat", J Exp Bot 71(7):2299-2311, doi:10.1093/jxb/erz439,
  !   Fig. 4B p.2305 (genotype-mean maximum 3.85; the Narea axis spans 0-5 and
  !   the +N field cluster tops near 3.5). At the lamina mass per area used
  !   here, csl_LMA = 48.5 g DM m-2 (Sales et al. 2022), 3.8 g N m-2 is 7.8 %DM,
  !   above any measured wheat lamina nitrogen concentration, which is why the
  !   capacity is applied per unit AREA and not as a concentration.
  ! - SLN_min_N = 0.50 g N m-2, the lower bound: Bertheloot, Martre & Andrieu
  !   (2008) Plant Physiol 148:1707-1720, Results p.1711, where lamina nitrogen
  !   concentration "decreased to a minimum value" at the end of grain filling;
  !   their Fig. 2C p.1711 gives the derived span 0.41-0.51 g N m-2, so 0.50
  !   sits at the top of the measured range. The same paper's own structural
  !   value is nb = 0.4 g N m-2 (Methods p.1718).
  ! The capacity window the reservoir obeys is SLN_min_N LAI <= N_leaf <=
  ! sln_max_N LAI, and the excess above it is
  !   MAX(0, N_leaf - SLN_min_N LAI) + MAX(0, N_stem - n_struct_stem_N W_stem).
  ! ------------------------------------------------------------------------- !
  INTEGER  :: leafN_drawdown = 0      ! runtime flag (driver_main.txt rec 62; EOF-default 0)
  REAL(wp), PARAMETER :: k_leafN_TT = 0.0035_wp ! (degree-Cd)-1, base 0 C (Bertheloot 2008 PP 148: value p.1712, Eq.5 p.1718)
  REAL(wp), PARAMETER :: sln_max_N  = 3.8_wp    ! g N m-2 leaf (Silva-Perez 2020 JXB 71(7):2299 Fig.4B p.2305, measured max 3.85)
  REAL(wp), PARAMETER :: n_stem_resid_conc = 0.0035_wp ! g N g-1 DM stem residual (Bertheloot 2008 p.1711, quote above)
  REAL(wp) :: dTT_day_N = 0.0_wp      ! day's thermal time (degree-Cd, base 0 C) from the phenology clock
  REAL(wp) :: nJ_R0  = 0.0_wp         ! remobilizable shoot N at fill start (flag ON; formula above)
  REAL(wp) :: nJ_cum = 0.0_wp         ! cumulative drawn against nJ_R0 (flag ON)

  ! ------------------------------------------------------------------------- !
  ! Post-anthesis grain-nitrogen partition package - driver record 63, which
  ! this version accepts at the value 0 only; the constants it declares are read
  ! by the active organ-nitrogen records, and their sources are recorded here.
  ! The law, per organ g and per day, with the day thermal time dTT (degree-Cd,
  ! base 0 C):
  !   N_ph_g = MAX(0, N_g - c_struct_g M_g)            photosynthetic N
  !   D_g    = delta N_ph_g dTT                        degradation -> mobile
  !   S_g    = sigma_g M_g f_N f_PAR dTT               xylem resupply <- mobile
  !   f_N    = [N]mob / ([N]mob + k1)                  [N]mob = N_mobile_veg / DM
  !   f_PAR  = PAR / (PAR + k2)
  ! and the living-mass update M_alive_new = N_ph_g / (d_tp (N_ph_g/M_g)|anth),
  ! with N_to_dead = c_struct_g (M_g - M_alive_new).
  ! Constants:
  ! - n2_delta_TT = 0.008 (degree-Cd)-1, the relative degradation rate of
  !   photosynthetic nitrogen: Bertheloot, Andrieu, Fournier & Martre (2008)
  !   Funct Plant Biol 35:781-796, doi:10.1071/FP08064, Table 1 p.784, and
  !   Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot 108:1097-1109,
  !   doi:10.1093/aob/mcr125 (Bertheloot 2011b), Table 2, row delta (0.008 degC d-1).
  !
  ! REFERENCE KEY for the author-year citations used throughout this file and in
  ! main_MMCO2Crop.f90. Titles are given exactly as published:
  !   Bertheloot 2011a = Bertheloot, J., Cournede, P.-H. and Andrieu, B.:
  !     "NEMA, a functional-structural model of nitrogen economy within wheat
  !     culms after flowering. I. Model description", Ann. Bot. 108, 1085-1096,
  !     doi:10.1093/aob/mcr119, 2011.
  !   Bertheloot 2011b = Bertheloot, J., Wu, Q., Cournede, P.-H. and Andrieu, B.:
  !     "NEMA, a functional-structural model of nitrogen economy within wheat
  !     culms after flowering. II. Evaluation and sensitivity analysis",
  !     Ann. Bot. 108, 1097-1109, doi:10.1093/aob/mcr125, 2011.
  ! PROVENANCE: the post-anthesis organ-nitrogen turnover of records 154 and 167
  ! is adopted from the NEMA model of those two papers. It is a model
  ! formulation taken from that source, not a measurement of this experiment.
  ! Part 1 Table B2 states what v1.0 adopts from it and where it departs.
  !   The same paper p.791 gives the net rate as about half the gross one,
  !   which is the bridge to the 0.0035 degree-Cd-1 lamina fit above.
  ! - n2_sigma_* (g g-1 (degree-Cd)-1), the relative synthesis rate from the
  !   xylem: Bertheloot 2011b Table 2, row sigma, and Bertheloot et al. (2008) Funct
  !   Plant Biol Table 1 p.784, in the same equation and the same units
  !   (kg kg-1 [degC d]-1 == g g-1 degC d-1).
  ! - n2_k1_mob = 0.0018 g g-1, the Michaelis constant on the culm-common
  !   mobile-nitrogen concentration: Funct Plant Biol Table 1 p.784 and Bertheloot 2011b
  !   Table 2, which print the same number.
  ! - n2_k2_PAR = 8.64e5 J m-2 d-1, the Michaelis constant paired with incident
  !   PAR. Funct Plant Biol Table 1 prints 10 J m-2 s-1 and Bertheloot 2011b Table 2
  !   prints 864 000 J m-2 d-1; 10 x 86400 = 864000 exactly, so the two are the
  !   same constant in different units and the daily form is used.
  ! - n2_struct_lam = 0.0079 g N g-1 DM, the lamina structural nitrogen
  !   concentration: Bertheloot, Martre & Andrieu (2008) Plant Physiol
  !   148:1707-1720, doi:10.1104/pp.108.124156, Results p.1711.
  ! - the chaff nitrogen mass at anthesis, 2.3 g N m-2 ground for cv. Apache,
  !   is read from Fig. 1D p.1709 of that paper, with the 1.43-fold higher
  !   Isengrain value given verbatim on p.1711 and the chaff dry masses in its
  !   Table I p.1711 (132 +/- 7 and 202 +/- 14 g DM m-2).
  ! ------------------------------------------------------------------------- !
  INTEGER  :: grainN_partition = 0    ! runtime flag (driver_main.txt rec 63; EOF-default 0)
  REAL(wp), PARAMETER :: n2_delta_TT    = 0.008_wp    ! (degree-Cd)-1, uniform over entity types; Bertheloot et al. 2008 FPB Table 1 p.784 (DOI 10.1071/FP08064) == Bertheloot 2011b Table 2 (DOI 10.1093/aob/mcr125)
  ! ------------------------------------------------------------------------- !
  ! Lamina and sheath synthesis rate sigma. Bertheloot, Wu, Cournede & Andrieu
  ! (2011) Ann Bot 108:1097-1109, doi:10.1093/aob/mcr125, Table 2, row
  ! "for entities tp, g g-1 degCd-1, 0.00015". The earlier wheat estimate of
  ! the SAME quantity, 0.0015, is printed in Bertheloot et al. (2008) Funct
  ! Plant Biol 35:781-796, Table 1 p.784, ten times higher. Both are coded, so
  ! that the choice between two published estimates of one quantity is explicit
  ! and nothing is refitted.
  REAL(wp), PARAMETER :: n2_sigma_lam   = 1.5e-4_wp      ! g g-1 (degree-Cd)-1, laminae + sheaths; Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109, Table 2, row and footnote "*" (DOI 10.1093/aob/mcr125)
  REAL(wp), PARAMETER :: n2_sigma_stem  = 3.0e-5_wp      ! g g-1 (degree-Cd)-1, internodes + peduncle; same table, row and footnote "dagger" (DOI 10.1093/aob/mcr125). The lumped stem pool of this model is assigned to this group, which is the conservative choice; the sheath value 1.5e-4 is the alternative
  REAL(wp), PARAMETER :: n2_sigma_chaff = 3.0e-5_wp      ! g g-1 (degree-Cd)-1, chaff; the footnote "dagger Internodes, peduncle and chaff" of the same table names the chaff explicitly in this group (DOI 10.1093/aob/mcr125)
  REAL(wp), PARAMETER :: n2_k1_Nmob     = 0.0018_wp   ! g g-1, MM constant on [N]mob; Bertheloot et al. 2008 FPB Table 1 p.784 (DOI 10.1071/FP08064) == Bertheloot 2011b Table 2 (DOI 10.1093/aob/mcr125)
  REAL(wp), PARAMETER :: n2_k2_PAR      = 8.64e5_wp   ! J m-2 d-1, MM constant on PAR; Bertheloot 2011b Table 2 (DOI 10.1093/aob/mcr125) 864000 == Bertheloot et al. 2008 FPB Table 1 p.784 (DOI 10.1071/FP08064) 10 J m-2 s-1 x 86400 s d-1, exact unit conversion
  REAL(wp), PARAMETER :: n2_struct_lam  = 0.0079_wp      ! g N g-1 DM lamina structural N; Bertheloot, Martre & Andrieu 2008 Plant Physiol 148:1707-1720, p.1711 (DOI 10.1104/pp.108.124156). The rate 0.0035 (degree-Cd)-1, which has the same digits, is a different quantity and is on p.1712
  REAL(wp), PARAMETER :: n2_struct_stem = 0.0035_wp      ! g N g-1 DM sheath/internode/chaff/peduncle; Bertheloot, Martre & Andrieu 2008 Plant Physiol 148:1707-1720, p.1711 (DOI 10.1104/pp.108.124156)

  ! ========================================================================== !
  ! Chaff structural nitrogen - driver record 178. Mode 0, the released value,
  ! gives the chaff the STEM structural floor n2_struct_stem. Modes 1 and 2 give
  ! it its own measured floor; they are an inactive alternative, see Part 1
  ! Table 6.
  ! n2_struct_chaff = 5.10e-3 g N g-1 DM is the chaff's own measured structural
  ! nitrogen: Bertheloot, Wu, Cournede & Andrieu (2011) "NEMA, a functional-
  ! structural model of nitrogen economy within wheat culms after flowering",
  ! Ann Bot 108:1097-1109, doi:10.1093/aob/mcr125, Table 1 section (c)
  ! "Photosynthetic and structural N", from 1.07 mg N per 210 mg DM =
  ! 5.0952e-3, written 5.10e-3.
  REAL(wp), PARAMETER :: n2_struct_chaff = 5.10e-3_wp    ! g N g-1 DM chaff structural N = 1.07 mg N / 210 mg DM; Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109, Table 1(c), column "Chaff" (DOI 10.1093/aob/mcr125)

  ! The live floor the chaff limb reads. It is initialized FROM the stem floor,
  ! so at record 178 = 0 every chaff expression is the stem expression and
  ! nothing moves.
  REAL(wp) :: n2_struct_chaff_use = n2_struct_stem  ! g N g-1 DM, the floor the CHAFF limb applies; = n2_struct_stem at rec 178 = 0, = n2_struct_chaff at rec 178 = 1

  ! ========================================================================== !
  ! Chaff nitrogen remobilization floor - driver record 178 modes 2 and 3, an
  ! inactive alternative, see Part 1 Table 6. The released configuration uses
  ! neither. Two measured floors are coded so that the alternative is explicit:
  ! (a) n2_struct_chaff_pask = 3.80e-3 g N g-1 DM = 0.38 %DM, the lowest harvest
  !     chaff nitrogen concentration observed at the zero-nitrogen treatment in
  !     Pask, Sylvester-Bradley, Jamieson & Foulkes (2012) "Quantifying
  !     quantitative traits for enhancing radiation use efficiency and yield
  !     potential in wheat", Field Crops Research 126:104-118,
  !     doi:10.1016/j.fcr.2011.09.021, Table 4 p.110, column CH of the
  !     "N remobilisation efficiency (%)" block; their Section 2.9(i) pp.107-108
  !     states verbatim that "structural N was assumed to remain constant" and
  !     Results p.111 that the lowest value observed for each plant component
  !     was used.
  ! (b) the chaff remobilization efficiency of Barraclough, Lopez-Bellido &
  !     Hawkesford (2014) "Genotypic variation in the uptake, partitioning and
  !     remobilisation of nitrogen during grain-filling in wheat", Front Plant
  !     Sci 5:artic. 227, pp.244-247, recovered from their printed anthesis and
  !     harvest chaff pools as 1 - 19/(0.23 x 189) = 0.563 at the high nitrogen
  !     rate and 1 - 6/(0.30 x 46) = 0.565 at the nil rate.
  ! Pask Table 4 spans 0.510-0.633 for the chaff across site-seasons, and
  ! reports the nitrogen effect on chaff remobilization efficiency as
  ! significant (P < 0.001).
  ! ------------------------------------------------------------------------- !
  REAL(wp), PARAMETER :: n2_struct_chaff_pask = 3.80e-3_wp  ! g N g-1 DM CHAFF structural minimum = 0.38 % DM, the lowest harvest chaff N% observed at N0 in any plot of three site-seasons; Pask, Sylvester-Bradley, Jamieson & Foulkes 2012 Field Crops Res 126:104-118 Methods 2.9(i) p.107-108 and Results p.111 (DOI 10.1016/j.fcr.2011.09.021)
  REAL(wp), PARAMETER :: cnr_nre_chaff        = 0.56_wp     ! dimensionless post-anthesis chaff N-remobilization efficiency; Pask et al. 2012 Table 4 p.110 column CH (three-site Nrec mean 55.7 %, abstract "chaff (56%)") and Barraclough, Lopez-Bellido & Hawkesford 2014 Field Crops Res 156:242 p.247 (0.563 at N200, 0.565 at N0)
  REAL(wp) :: cnr_N_chaff_anth = 0.0_wp   ! g N m-2, the chaff's N on the FIRST live fill day, before any post-anthesis release; snapshot, never re-taken
  REAL(wp) :: cnr_floor_gN     = 0.0_wp   ! g N m-2, (1 - cnr_nre_chaff) x cnr_N_chaff_anth; IDENTICALLY 0.0 at rec 178 = 0 and 1, which is what makes those modes
  LOGICAL  :: cnr_snap_done    = .FALSE.  ! has the anthesis snapshot been taken this season?
  INTEGER  :: cnr_days_bound   = 0        ! number of fill days on which the NRE floor, not the concentration floor, was the binding one
  REAL(wp) :: cnr_held_gN     = 0.0_wp   ! g N m-2 STATE (not a sum): the current excess of the NRE floor over the concentration floor, i.e. the chaff N the ratio law is holding that the concentration floor alone would have released
  REAL(wp), PARAMETER :: n2_d_tp        = 0.4_wp      ! dimensionless tissue-death threshold; Bertheloot 2011b Table 2 (DOI 10.1093/aob/mcr125); mass-basis substitution stated in the header
  REAL(wp), PARAMETER :: n2_chaff_nmax  = 0.017_wp        ! g N g-1 DM chaff luxury ceiling; Bertheloot, Martre & Andrieu 2008 Plant Physiol 148:1707-1720, Fig. 1D p.1709 and Table I p.1711 (DOI 10.1104/pp.108.124156); a figure-level read
  REAL(wp), PARAMETER :: n2_par_frac    = 0.45_wp     ! PAR fraction of incoming global SW; the model's OWN existing constant, solar.f90:88 and :92 (no new number); STRUCTURAL re-use
  ! --- item-2 state (all identically zero when the flag is OFF --------- !
  REAL(wp) :: W_dead_stem = 0.0_wp    ! g DM m-2: dead stem/sheath + dead chaff dry mass (piece a)
  REAL(wp) :: N_dead_stem = 0.0_wp    ! g N m-2: structural N buried in that pool (piece a)
  LOGICAL  :: n2_anth_snap = .FALSE.  ! flowering snapshot of photosynthetic [N] taken?
  REAL(wp) :: n2_cph_stem_anth  = 0.0_wp ! g N g-1 DM, stem photosynthetic [N] at flowering (d_tp reference)
  REAL(wp) :: n2_cph_chaff_anth = 0.0_wp ! g N g-1 DM, chaff photosynthetic [N] at flowering (d_tp reference)
  ! --- item-2 daily diagnostics (evaluation ledgers; never fluxes ------ !
  REAL(wp) :: n2_deg_lam_day = 0.0_wp, n2_deg_stem_day = 0.0_wp, n2_deg_chaff_day = 0.0_wp
  REAL(wp) :: n2_syn_lam_day = 0.0_wp, n2_syn_stem_day = 0.0_wp, n2_syn_chaff_day = 0.0_wp
  REAL(wp) :: n2_dW_dead_stem_day = 0.0_wp, n2_dW_dead_chaff_day = 0.0_wp
  REAL(wp) :: n2_cum_deg_lam = 0.0_wp, n2_cum_deg_stem = 0.0_wp, n2_cum_deg_chaff = 0.0_wp

  ! ------------------------------------------------------------------------- !
  ! Per-culm leaf-rank chronology - driver record 64, accepted at the value 1,
  ! so the chronology below runs.
  ! THE CLOCK. tau is cumulative degree-Cd with base 0 C from emergence; the
  ! increment is the definition of McMaster, Wilhelm, Palic, Porter & Jamieson
  ! (2003) Ann Bot 91:697-705, doi:10.1093/aob/mcg074, Eqs. 1-2 p.699,
  !   GDD = (Tmax + Tmin)/2 - Tbase,  Tbase = 0 C,
  ! so calc_dtt returns MAX(0, Tavg).
  ! PRIMARY CHRONOLOGY (clr_HS), the Haun stage carrying the appearance-rate
  ! decline with leaf number:
  !   HS <  2 :  dHS/dtau = 1 / (P clr_seed_ratio)
  !   HS >= 2 :  dHS/dtau = (1/P) (HS/2)**(-0.3)
  ! with the phyllochron P from McMaster et al. (2003) Table 2 p.701 (default
  ! 73.9 degree-Cd) and the (HS/2)**(-0.3) decline and the seed-reserve ratio
  ! from Streck, Weiss & Baenziger (2003) Ann Bot 92:181-190, Eq. 5 and
  ! pp.183-184, doi:10.1093/aob/mcg121.
  ! LIGULATION. tau_lig(n) = tau_app(n) + LED, with the leaf elongation
  ! duration LED = 135 degree-Cd (base 0 C) of Bos, H. J. (1999) PhD thesis,
  ! Wageningen, Ch.3, doi:10.18174/198200.
  ! The alternative appearance intervals of Evers, J. B. (2006) PhD thesis,
  ! Wageningen, doi:10.18174/24302, Ch.3 Appendix B (52.17 degree-Cd for
  ! main-stem leaves 1-4 and 91.98 for later leaves) are coded as a comparison
  ! chronology and run without the Streck decline.
  ! Flag-leaf number is read on the Haun-stage convention of Abichou, M. (2016)
  ! PhD thesis (NNT tel-02316154), Eq. (1): HS(tau) = n + l/L with n the number
  ! of LIGULATED leaves, so the flag-leaf ligule corresponds to the Zadoks
  ! flag-leaf stage (Zadoks, Chang & Konzak 1974, Weed Res 14:415-421).
  ! ------------------------------------------------------------------------- !
  INTEGER, PARAMETER :: CLR_MAXRANK = 24   ! array bound only; far above any wheat FLN (Abichou 2016 doc p.66 panel spans FLN 8-14)
  INTEGER  :: canopy_leaf_ranks = 0        ! runtime flag (driver_main.txt rec 64; EOF-default 0)
  REAL(wp) :: clr_phyllochron_Cd = 73.9_wp    ! degree-Cd, air basis; McMaster, Wilhelm, Palic, Porter & Jamieson 2003 Ann Bot 91:697-705, Table 2 p.701 (DOI 10.1093/aob/mcg074), the mean of the two season means 73.0 and 74.8 for one spring cultivar; the per-planting-date bracket in the same table is 68.4-77.4. Cultivar-settable: set from phenology_cultivar.txt row Phyllochron_degCd, which defaults to this value. The quantity is cultivar-scoped, never site-specific
  REAL(wp), PARAMETER :: clr_seed_ratio     = 0.75_wp  ! -, leaves 1-2 run at P*clr_seed_ratio; Streck et al. 2003 Ann Bot 92:181-190 pp.183-184 (DOI 10.1093/aob/mcg121): R_LA,max 0.24 leaves d-1 for leaves > 2 against ~0.32 for the first two (seed reserves), 0.24/0.32
  REAL(wp), PARAMETER :: clr_fC_exp         = -0.3_wp  ! -, chronology-decline exponent; Streck et al. 2003 Eqn 5 p.183 (DOI 10.1093/aob/mcg121): f(C) = (HS/2)**(-0.3)
  REAL(wp), PARAMETER :: clr_fC_HS_ref      = 2.0_wp   ! -, Haun stage at which f(C) starts and equals 1; Streck et al. 2003 Eqn 5 p.183 (DOI 10.1093/aob/mcg121)
  REAL(wp), PARAMETER :: clr_LED_Cd         = 135.0_wp ! degree-Cd (Tbase 0), tip appearance to ligulation; Bos 1999 thesis Ch.3 (DOI 10.18174/198200) measured absolute range 123-147, midpoint; bracket swept in the S1 progress document
  REAL(wp), PARAMETER :: clr_P_evers_low    = 52.17_wp ! degree-Cd, main-stem leaves 1-4; Evers 2006 thesis Ch.3 Appendix B doc p.72 (DOI 10.18174/24302), 52.17 +/- 1.24
  REAL(wp), PARAMETER :: clr_P_evers_high   = 91.98_wp ! degree-Cd, all later main-stem leaves and tiller leaves; Evers 2006 thesis Ch.3 Appendix B doc p.72 (DOI 10.18174/24302), 91.98 +/- 2.42
  REAL(wp), PARAMETER :: clr_evers_break_HS = 4.0_wp   ! -, Haun stage at which the Evers split changes phyllochron; Evers 2006 thesis Ch.3 Appendix B doc p.72 (DOI 10.18174/24302), "leaves 1-4" against "all other"
  ! --- S1 state (all identically zero / sentinel when the flag is OFF --- !
  LOGICAL  :: clr_active   = .FALSE.  ! emergence reached (tau clock started)
  LOGICAL  :: clr_frozen   = .FALSE.  ! FLN read off at the anchored Booting stage
  REAL(wp) :: clr_tau_Cd   = 0.0_wp   ! degree-Cd base 0 since emergence
  REAL(wp) :: clr_HS       = 0.0_wp   ! Haun stage, primary chronology (McMaster P + Streck f(C))
  REAL(wp) :: clr_HS_ev    = 0.0_wp   ! Haun stage, Evers split-phyllochron variant (no f(C))

  ! Ligulation times per rank on the primary arm (-1 = not ligulated). Under the
  ! Haun-stage convention of Abichou (2016) Eq. (1), HS(tau) = n + l/L with n
  ! the number of ligulated leaves, the Haun-stage crossing IS the ligulation
  ! time of that rank.
  REAL(wp) :: clr_tau_lig(CLR_MAXRANK)    = -1.0_wp  ! degree-Cd of LIGULATION per rank, primary arm (-1 = not ligulated); the HS crossing itself
  REAL(wp) :: clr_tau_app(CLR_MAXRANK)    = -1.0_wp  ! degree-Cd of TIP APPEARANCE per rank, primary arm = tau_lig - LED (-1 = not set)
  REAL(wp) :: clr_tau_lig_ev(CLR_MAXRANK) = -1.0_wp  ! degree-Cd of ligulation per rank, Evers arm (-1 = not ligulated)
  REAL(wp) :: clr_tau_app_ev(CLR_MAXRANK) = -1.0_wp  ! degree-Cd of tip appearance per rank, Evers arm (-1 = not set)
  INTEGER  :: clr_day_app(CLR_MAXRANK)    = 0        ! model day index of tip appearance per rank, primary arm
  INTEGER  :: clr_day_lig(CLR_MAXRANK)    = 0        ! model day index of ligulation per rank, primary arm
  INTEGER  :: clr_n_app    = 0        ! ranks with a tip appeared (primary arm)
  INTEGER  :: clr_n_lig    = 0        ! ranks ligulated (primary arm)
  INTEGER  :: clr_n_app_ev = 0        ! ranks with a tip appeared (Evers arm)
  INTEGER  :: clr_n_lig_ev = 0        ! ranks ligulated (Evers arm)
  INTEGER  :: clr_FLN      = 0        ! EMERGENT final leaf number, ligule basis at the anchored Booting stage
  INTEGER  :: clr_FLN_tip  = 0        ! tip-basis leaf count at the anchored Booting stage (upper bound)
  INTEGER  :: clr_FLN_ev   = 0        ! same, Evers side-by-side arm
  INTEGER  :: clr_FLN_ev_tip = 0      ! same, Evers side-by-side arm, tip basis
  REAL(wp) :: clr_tau_boot = 0.0_wp   ! degree-Cd since emergence at the anchored Booting stage
  REAL(wp) :: clr_tau_app_flag = 0.0_wp  ! degree-Cd of flag-leaf tip appearance (primary arm)
  REAL(wp) :: clr_tau_lig_flag = 0.0_wp  ! degree-Cd of flag-leaf ligulation (primary arm)

  ! ------------------------------------------------------------------------- !
  ! Potential and realized blade area per rank - driver record 65, accepted at
  ! the value 1. Area = L W f_t per rank, with:
  ! - final blade length and maximum blade width by inversed rank read from the
  !   measured profiles of Abichou, M. (2016) PhD thesis (NNT tel-02316154),
  !   Ch.3 Fig. 8a (length, the NI treatment, median of the five cultivars) and
  !   Ch.2 Fig. 3e (width, the central series of the Soissons panel). These are
  !   figure-level reads: the thesis prints the polynomial forms (Table 2
  !   Eqns T2.1 and T2.4) but not all of their coefficients, so the measured
  !   profiles are used directly rather than reconstructed.
  ! - the form factor f_t, rank-resolved rather than flat: Dornbusch, Watt,
  !   Baccar, Fournier & Andrieu (2011) Ann Bot 107(5):865-873,
  !   doi:10.1093/aob/mcq181, Fig. 3A, wheat cv. Soissons - the same cultivar as
  !   the Abichou profiles.
  ! - the tiller dimension shifts s_Ti and the reduction coefficients from
  !   Abichou (2016) Table 4, together with the laws of Table 2.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: canopy_blade_area = 0     ! runtime flag (driver_main.txt rec 65; EOF-default 0)

  ! Final blade length by inversed rank (cm). Abichou (2016) Ch.3 Fig. 8a, the
  ! NI series, median of the five cultivars; a figure-level read, because the
  ! thesis does not print all coefficients of its Eqn T2.1 polynomial.
  REAL(wp), PARAMETER :: cba_l_irank_ref(12) = [16.8_wp, 20.8_wp, 20.3_wp, 17.3_wp, 14.5_wp, 11.3_wp, &
                                            11.0_wp, 11.2_wp,  9.0_wp,  9.0_wp,  9.0_wp,  9.0_wp]
  REAL(wp) :: cba_l_irank(12) = cba_l_irank_ref  ! WORKING profile (see note above); default = cba_l_irank_ref verbatim
  REAL(wp) :: Blade_L_max_cm = 20.8_wp  ! cm, cultivar scalar = cba_l_irank_ref(2) (current profile maximum); Phase D; default UNCHANGED; settable from phenology_cultivar.txt row Blade_L_max_cm (EOF-default keeps this literal)
  INTEGER,  PARAMETER :: cba_Nl_fit = 11 ! typical Nl_MS of the Fig. 8 panel culms; used ONLY as the preliminary top before the emergent FLN is frozen at booting

  ! Maximum blade width by inversed rank (cm). Abichou (2016) Ch.2 Fig. 3e, the
  ! central series of the 15-treatment Soissons panel; a figure-level read, the
  ! coefficients of Eqn T2.4 not being printed in the thesis.
  REAL(wp), PARAMETER :: cba_w_irank_ref(12) = [1.45_wp, 1.45_wp, 1.25_wp, 1.05_wp, 0.85_wp, 0.65_wp, &
                                            0.55_wp, 0.50_wp, 0.45_wp, 0.42_wp, 0.38_wp, 0.38_wp]
  REAL(wp) :: cba_w_irank(12) = cba_w_irank_ref  ! WORKING profile (see note above); default = cba_w_irank_ref verbatim
  REAL(wp) :: Blade_W_max_cm = 1.45_wp  ! cm, cultivar scalar = cba_w_irank_ref(1:2) (current profile maximum); Phase D; default UNCHANGED; settable from phenology_cultivar.txt row Blade_W_max_cm (EOF-default keeps this literal)
  ! Form factor anchors, wheat cv Soissons; Dornbusch et al. 2011 Ann Bot 107:865-873 (DOI 10.1093/aob/mcq181) Fig. 3A — see header for the quoted increments this reconstruction rests on
  REAL(wp), PARAMETER :: cba_ft_r1 = 0.90_wp, cba_ft_r4 = 0.75_wp, cba_ft_adult = 0.80_wp, cba_ft_flag = 0.75_wp
  INTEGER,  PARAMETER :: cba_ft_n4 = 4, cba_ft_n8 = 8   ! the rank breakpoints Dornbusch prints (1<=N<=4, 4<=N<=8, 8<=N<=Nmax-1)
  ! Tiller constants; Abichou 2016 Table 4 doc p.71 (verified verbatim)
  REAL(wp), PARAMETER :: cba_sTi(6) = [2.755_wp, 3.427_wp, 4.091_wp, 4.588_wp, 4.899_wp, 5.322_wp]  ! decimal phytomer shift of tillers T1-T6
  REAL(wp), PARAMETER :: cba_lb1    = 7.023_wp   ! cm, final length of the first blade of a tiller (constant for all tiller ranks)
  REAL(wp), PARAMETER :: cba_kw     = 0.877_wp   ! -, fertile-tiller blade-width reduction vs the main stem (Table 2 Eqn T2.10)
  REAL(wp), PARAMETER :: cba_kl_st  = 0.885_wp   ! -, sterile-tiller length reduction vs a fertile tiller (Eqn T2.13)
  REAL(wp), PARAMETER :: cba_kw_st  = 0.906_wp   ! -, sterile-tiller width reduction vs a fertile tiller (Eqn T2.16)
  ! --- S2 state (all identically zero when the flag is OFF ------------- !
  REAL(wp) :: cba_L_pot(CLR_MAXRANK)   = 0.0_wp  ! cm, final blade length per main-stem rank
  REAL(wp) :: cba_W_pot(CLR_MAXRANK)   = 0.0_wp  ! cm, maximum blade width per main-stem rank
  REAL(wp) :: cba_ft_rank(CLR_MAXRANK) = 0.0_wp  ! -, form factor per main-stem rank
  REAL(wp) :: cba_A_pot(CLR_MAXRANK)   = 0.0_wp  ! cm2, potential one-sided blade area per main-stem rank
  REAL(wp) :: cba_A(CLR_MAXRANK)       = 0.0_wp  ! cm2, realized (expanding) blade area per main-stem rank
  REAL(wp) :: cba_A_shoot_pot = 0.0_wp  ! cm2 per main stem, sum of A_pot over ranks that have appeared
  REAL(wp) :: cba_A_shoot     = 0.0_wp  ! cm2 per main stem, sum of realized A over ranks that have appeared
  INTEGER  :: cba_n_expanding = 0       ! ranks currently between tip appearance and ligulation

  ! ------------------------------------------------------------------------- !
  ! Axis (tiller) population and the ground-area basis - driver record 66,
  ! accepted at the value 1.
  ! - COHORT BIRTH TIMES: Abichou (2016) PhD thesis (NNT tel-02316154) Ch.1
  !   Table 6, the measured main-stem Haun stage at the emergence of tillers
  !   T1-T5. No measured emergence Haun stage is printed for T6, so the T6
  !   cohort is modeled rather than measured. These birth times are distinct
  !   from the Ch.2 Table 4 dimension shifts s_Ti and must not be interchanged.
  ! - BIRTH PROBABILITY: Evers, Vos, Andrieu & Struik (2006) Ann Bot
  !   97:649-658, doi:10.1093/aob/mcl020, Results, Fig. 5A, where the
  !   probability that a tiller bud grows out declines with the red:far-red
  !   signal of the canopy; the light interception that drives it is taken from
  !   the model's own green area index, P_int = 1 - exp(-k LAI).
  ! - SITE-LEVEL RATES: Bos (1999) PhD thesis, doi:10.18174/198200, Ch.2
  !   pp.14-17, whose 0.69 is a theoretical maximum of appearance rate and 0.67
  !   the measured winter-wheat value, both whole-plant quantities.
  ! - CESSATION BACKSTOP: Abichou (2016) Ch.3 Fig. 3, where axis density peaks
  !   and stops rising; the measurement basis is axes with at least one green
  !   leaf.
  ! - MORTALITY: Abichou (2016) Ch.1 Table 4 with the fits of its Section 3.6 -
  !   Eqn T4.2 for the start of the mortality phase at a main-stem Haun stage
  !   (r2 = 0.78, Fig. 10) and Eqn T4.4 for the normalized pattern
  !     Na(tau) = exp(-theta taun/(delta - taun)) (Ne - Namax) + Ne,
  !     taun = (tau - tau_SM)/(tau_EM - tau_SM).
  ! - PER-AXIS GEOMETRY: Abichou (2016) Ch.1 Table 3 Eqn T3.5, Nl_Ti = alpha
  !   Nl_MS - beta, and Ch.2 Table 2 Eqns T2.7a/b for the tiller dimensions.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: canopy_tillering = 0      ! runtime flag (driver_main.txt rec 66; EOF-default 0)
  INTEGER,  PARAMETER :: CT_MAXTIL = 5  ! primary tillers T1-T5; T6 has no measured emergence Haun stage (Abichou Ch.1 Table 6 stops at T5); omitted
  REAL(wp), PARAMETER :: ct_hs_emerg(CT_MAXTIL) = [2.5_wp, 3.2_wp, 4.3_wp, 5.5_wp, 6.9_wp]  ! main-stem Haun stage at emergence of T1-T5; Abichou 2016 Ch.1 Table 6 doc p.31
  REAL(wp), PARAMETER :: ct_pint_p1 = 0.40_wp   ! P_int below which the occurrence probability of a primary tiller is 1; Evers et al. 2006 Ann Bot 97:649-658 Fig. 5A
  REAL(wp), PARAMETER :: ct_pint_p0 = 0.55_wp   ! P_int at which it has fallen to 0; same figure, same sentence

  ! Canopy PAR interception coefficient used by the tillering light signal:
  ! 0.43, Moreau et al. (2012) Plant Physiol 160(3):1479-1490, Table II
  ! (0.42 and 0.43 in the two controlled-field experiments).
  REAL(wp), PARAMETER :: ct_kL      = 0.43_wp   ! measured wheat canopy PAR interception coefficient; Moreau et al. 2012 Plant Physiol 160(3):1479 Table II (0.42/0.43/0.46) — see item E-3 of the revision record (internal; not distributed)
  REAL(wp) :: ct_axmax    = 700.0_wp            ! axes m-2, empirical cessation cap; Abichou Ch.3 Fig. 3 doc p.100 (600-800; 600 excluded as central)
  REAL(wp) :: ct_die_frac = 0.632_wp                        ! -, tiller die fraction at 200 plants m-2; interpolated between the measured 0.57 at 100 and 0.67 at 262 plants m-2 of Evers (2006)
  REAL(wp), PARAMETER :: ct_hs_sm_lo = 5.5_wp, ct_hs_sm_hi = 8.0_wp   ! printed range of HS_MS at the start of tiller mortality; Abichou Ch.1 section 3.5
  REAL(wp), PARAMETER :: ct_nl_sm_lo = 10.0_wp, ct_nl_sm_hi = 14.0_wp ! the final-leaf-number range those endpoints were measured over; Abichou Ch.1 section 3.3
  REAL(wp), PARAMETER :: ct_theta = 3.32_wp, ct_delta = 1.77_wp       ! Eqn T4.4 mortality-pattern constants; Abichou Ch.1 section 3.6 (RMSE 0.10)
  REAL(wp), PARAMETER :: ct_nl_alpha = 0.942_wp, ct_nl_beta = 0.623_wp ! Eqn T3.5 tiller final-leaf-number law; Abichou Ch.1 section 3.3 (RMSE 0.42 leaves)
  INTEGER,  PARAMETER :: ct_nc_offset = 3       ! cohort number Nc_Ti = i + 3 (Ch.1 symbol table: Nc_T0 = 3, Nc_T1 = 4)
  INTEGER,  PARAMETER :: ct_internodes_elong = 4 ! elongated internodes, so r = Nl_MS - 4; Abichou Ch.2 doc p.66 ("generally 4 to 5"), swept as 5
  ! --- S3 state (all identically zero/false when the flag is OFF -------- !
  LOGICAL  :: ct_active   = .FALSE.
  LOGICAL  :: ct_ceased   = .FALSE.
  INTEGER  :: ct_cess_why = 0           ! 0 = tillering still open, 1 = Evers light law, 2 = empirical axis cap
  REAL(wp) :: ct_fpar     = 0.0_wp      ! P_int = 1 - exp(-k*LAI), the model's own like-for-like construction of Evers's measurement
  REAL(wp) :: ct_dens(0:CT_MAXTIL)      = 0.0_wp  ! axes m-2 now alive, index 0 = main stem
  REAL(wp) :: ct_dens_born(0:CT_MAXTIL) = 0.0_wp  ! axes m-2 at birth (before any mortality)
  REAL(wp) :: ct_p_occ(CT_MAXTIL)       = 0.0_wp  ! occurrence probability of each cohort (Evers Fig. 5A)
  REAL(wp) :: ct_pint_win(CT_MAXTIL)    = 0.0_wp  ! mean P_int over that cohort's window of opportunity
  REAL(wp) :: ct_win_sum(CT_MAXTIL)     = 0.0_wp  ! running sum for that window mean
  REAL(wp) :: ct_win_n(CT_MAXTIL)       = 0.0_wp  ! running day count for that window mean
  LOGICAL  :: ct_decided(CT_MAXTIL)     = .FALSE. ! window closed and p_occ fixed
  LOGICAL  :: ct_appeared(0:CT_MAXTIL)  = .FALSE.
  REAL(wp) :: ct_tau_birth(0:CT_MAXTIL) = -1.0_wp ! degCd from emergence at which each axis class appeared
  REAL(wp) :: ct_nl_axis(0:CT_MAXTIL)   = 0.0_wp  ! final leaf number of each axis class (Eqn T3.5)
  REAL(wp) :: ct_area_axis(0:CT_MAXTIL) = 0.0_wp  ! cm2 per axis, potential blade area of that class
  REAL(wp) :: ct_area_axis_real(0:CT_MAXTIL) = 0.0_wp ! cm2 per axis, realized (expanding) blade area
  REAL(wp) :: ct_area_axis_st(0:CT_MAXTIL)   = 0.0_wp ! cm2 per axis, potential area under the STERILE reductions (labeled sensitivity)
  REAL(wp) :: ct_axes_live = 0.0_wp     ! axes m-2 with at least one green leaf — the quantity Abichou's bound is measured on
  REAL(wp) :: ct_axes_born = 0.0_wp     ! axes m-2 ever born
  REAL(wp) :: ct_axes_max  = 0.0_wp     ! Namax, the plateau of tillering
  REAL(wp) :: ct_ears      = 0.0_wp     ! Ne, predicted final ear density (axes m-2 surviving the mortality phase)
  REAL(wp) :: ct_tau_sm = -1.0_wp, ct_tau_em = -1.0_wp  ! degCd at the start / end of the mortality phase
  REAL(wp) :: ct_hs_sm  = 0.0_wp        ! main-stem Haun stage at the start of mortality (Eqn T4.2 reconstruction)
  REAL(wp) :: ct_bai_pot = 0.0_wp       ! m2 m-2, GROUND-AREA potential blade area index = sum over classes of density * per-axis area
  REAL(wp) :: ct_bai     = 0.0_wp       ! m2 m-2, ground-area realized (expanding, no senescence until S4) blade area index
  REAL(wp) :: ct_fs_implied = 0.0_wp    ! -, Bos 1999 site filling implied by the model population, ln(axes per plant)/HS — REPORTED cross-check only

  ! ------------------------------------------------------------------------- !
  ! Blade senescence and the standing-dead pool - driver record 67, accepted at
  ! the value 1.
  ! THE GREEN-LEAF-NUMBER LAW: Abichou (2016) PhD thesis (NNT tel-02316154)
  ! Ch.2 Table 3. The green leaf number of an axis has four phases, the main
  ! stem being
  !   phi1  Gl = HS                                              tau  <= tau1
  !   phi2  Gl = (g2-g1)(tau-tau1)/(tau2-tau1) + g1              tau1 <= tau <= tau2
  !   phi3  Gl = (g3-g2)(tau-tau2)/(tau3-tau2) + g2              tau2 <= tau <= tau3
  !   phi4  Gl = alpha (tau-tau3)^3 + beta (tau-tau3) + g3       tau3 <= tau <= tau4
  ! with tau3 the ligulation of the FLAG leaf, and the senescent leaf number of
  ! an axis given by Ch.2 Eqn 3, Sl_axis(tau) = HS(tau) - Gl(tau).
  ! FITTED VALUES: Abichou (2016) Ch.2 Table S1. g1 is transferred as a Haun
  ! stage rather than as a degree-Cd time, because that is what the table
  ! prints.
  ! TILLER REDUCTION delta: Abichou (2016) Ch.2 Table 6.
  ! RANK-PROGRESSIVE LIFESPAN: Robert, Garin, Abichou, Houles, Pradal &
  ! Fournier (2018) Ann Bot 121:975-989, doi:10.1093/aob/mcx192, Results.
  ! Two measured sources disagree about the onset of senescence; the
  ! disagreement is reported by the gate rather than averaged away, and the
  ! source used is the one that measures this quantity directly
  ! directly (Abichou Ch.2 Fig. 2d).
  ! ------------------------------------------------------------------------- !
  INTEGER  :: canopy_leaf_senesc = 0    ! runtime flag (driver_main.txt rec 67; EOF-default 0)
  REAL(wp) :: cls_g1 = 4.9_wp           ! leaves, green leaf number (= Haun stage) at senescence onset; Abichou Ch.2 Table S1 doc p.86, median of the four D2+N main stems (4.7-5.2); panel 2.9-5.4
  REAL(wp) :: cls_g2 = 3.3_wp           ! leaves, minimum green leaf number, at ligulation of leaf Nl-4; same rows (3.1-3.6); panel 2.6-4.2
  REAL(wp) :: cls_g3 = 5.45_wp          ! leaves, maximum green leaf number, at flag ligulation; same rows (4.4-5.8); panel 4.0-6.6
  REAL(wp) :: cls_alpha = -4.1e-9_wp    ! (degCd)^-3, phi4 cubic coefficient; same rows (-2.9e-9..-8.8e-9); panel -2.6e-9..-9.0e-9
  REAL(wp) :: cls_beta  = -2.7e-3_wp    ! (degCd)^-1, phi4 linear coefficient; all four D2+N rows are -2.7e-3; panel -1.0e-3..-2.7e-3
  REAL(wp) :: cls_delta = 0.85_wp       ! -, g3_Ti/g3_MS; Abichou Ch.2 Table 6 doc p.71, median of the three D2 N+ values (0.79/0.85/0.96); bracket 0.79-1.00
  INTEGER,  PARAMETER :: cls_top_offset = 4  ! phi2 ends at the ligulation of leaf Nl-4, "the fourth leaf from the top" (Ch.2 doc p.64) — structural, from the equations themselves
  ! --- S4 state (all identically zero when the flag is OFF ------------- !
  LOGICAL  :: cls_active = .FALSE.
  REAL(wp) :: cls_gl(0:CT_MAXTIL)    = 0.0_wp  ! green leaf number per axis class (0 = main stem)
  REAL(wp) :: cls_sl(0:CT_MAXTIL)    = 0.0_wp  ! senescent leaf number per axis class, Eqn 3 = HS - Gl
  INTEGER  :: cls_phase(0:CT_MAXTIL) = 0       ! which of the four phases that axis is in
  REAL(wp) :: cls_tau1(0:CT_MAXTIL)  = -1.0_wp ! degCd at senescence onset of that axis
  REAL(wp) :: cls_tau2(0:CT_MAXTIL)  = -1.0_wp ! degCd at ligulation of its leaf Nl-4
  REAL(wp) :: cls_tau3(0:CT_MAXTIL)  = -1.0_wp ! degCd at ligulation of its flag leaf
  REAL(wp) :: cls_area_green(0:CT_MAXTIL) = 0.0_wp ! cm2 per axis, GREEN blade area
  REAL(wp) :: cls_area_dead(0:CT_MAXTIL)  = 0.0_wp ! cm2 per axis, STANDING DEAD blade area (never shed)
  REAL(wp) :: cls_gai      = 0.0_wp     ! m2 m-2, ground-area GREEN blade area index
  REAL(wp) :: cls_dai      = 0.0_wp     ! m2 m-2, ground-area STANDING DEAD blade area index
  REAL(wp) :: cls_dai_max  = 0.0_wp     ! m2 m-2, running maximum of the standing-dead index (monotonicity check)
  INTEGER  :: cls_mono_viol = 0         ! days on which the standing-dead index DECREASED (must stay 0)
  REAL(wp) :: cls_dai_pool = 0.0_wp          ! m2 m-2, the standing-dead pool proper: an accumulator of the area lost each day, monotone non-decreasing by construction. cls_dai is the instantaneous re-computation kept beside it for diagnosis
  REAL(wp) :: cls_dai_dip_max = 0.0_wp  ! m2 m-2, largest single-day DECREASE of the instantaneous recomputation — reported, see the S4 progress doc (it is the booting FLN rebase, not senescence running backwards)
  REAL(wp) :: cls_dai_prev = -1.0_wp    ! previous day's instantaneous dead index
  REAL(wp) :: cls_dai_culled = 0.0_wp   ! m2 m-2, blade area carried by axes that DIED in the S3 mortality phase, retained standing (never shed)
  REAL(wp) :: cls_dens_prev(0:CT_MAXTIL) = 0.0_wp ! previous day's axis densities, to detect culling and retire the culled axes' area
  REAL(wp) :: cls_gl_at_flag = -1.0_wp  ! green leaf number of the main stem AT flag ligulation (the g3 gate)
  REAL(wp) :: cls_hs_onset   = -1.0_wp  ! main-stem Haun stage at which senescence actually started (the E11 gate)

  ! ------------------------------------------------------------------------- !
  ! Non-lamina standing area: sheaths, exposed stem and ear - driver record 68,
  ! accepted at the value 1. The measurement these organs are compared against
  ! is the green plant area index of Abichou (2016) Ch.3, which counts all green
  ! phytoelements, not the laminae alone.
  ! - SHEATH LENGTH by inversed rank: Abichou (2016) PhD thesis (NNT
  !   tel-02316154) Ch.2 Fig. 3b, the medians of the 15-treatment Soissons
  !   panel, with the profile shape verified against the thesis's own description
  !   (Ch.2 p.66). A figure-level read: Table 2 Eqn T2.2 is a degree-5
  !   polynomial whose coefficients Table S2 does not print in full.
  ! - SHEATH DIAMETER: Abichou (2016) Ch.2 Table 2 Eqns T2.5a and T2.5b,
  !     Ds_MS(x) = (ds2 - ds1)/(r - 1) (x - 1) + ds1   for 1 <= x <= r,
  !   with ds1 = 0.124 cm a measured constant of the sub-model (Table 4), ds2
  !   the mean diameter of the upper sheaths (Table S2, bracket 0.21-0.49 cm)
  !   and r - 1 the rank of the lowest elongated internode, generally 4 to 5
  !   (Ch.2 p.66). Sub-model accuracy, Table 5: sheath length RMSE 3.675 cm
  !   (MAPE 11 %, n = 493), sheath diameter RMSE 0.067 cm (MAPE 18 %).
  ! - TILLER SHEATHS: Abichou (2016) Ch.2 Table 2 Eqns T2.8 and T2.11, with the
  !   same Table 4 shifts s_Ti the blade sub-model already uses, and the tiller
  !   diameter coefficient ko = MAX(rho1 pTi + rho2, rho3) with rho1 = -0.265,
  !   rho2 = 1.007 and rho3 = 0.702 (Table 4).
  ! - EXPOSED STEM: exposed_length(rank j) = MAX(0, Li(j) - Ls(j+1)) with the
  !   internode lengths Li read from Abichou (2016) Ch.2 Fig. 3c, the same
  !   figure family and panel as the sheath lengths.
  ! - EAR AREA and the ear window: Abichou (2016) Annex 2 ear table and Ch.3
  !   p.96, which measures the window from flag-leaf ligulation.
  ! - INDEPENDENT CROSS-CHECK: Zhang, Chen, Jing, Gao & Wang (2023) "Canopy
  !   structure, light interception and photosynthetic characteristics of wheat",
  !   Plants 12:3482, doi:10.3390/plants12193482, Table 1, which measures the
  !   area and the share of the green photosynthetic organs per stem, including
  !   the stem surface area L1 pi D1 + L2 pi D2.
  ! No field-measured source was located for the TIMING of sheath senescence;
  ! the model applies the blade-and-sheath extension description of Abichou
  ! (2016) Annex 2, and the gap is stated rather than filled with a number.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: canopy_sheath_ear = 0     ! runtime flag (driver_main.txt rec 68; EOF-default 0)

  ! Final SHEATH length by inversed rank (cm). Abichou (2016) Ch.2 Fig. 3b,
  ! medians of the 15-treatment Soissons panel, with rank 1 set at the midpoint
  ! of the 14-20 cm largest-sheath range the thesis prints; a figure-level read,
  ! because Eqn T2.2 and Table S2 do not supply all coefficients.
  REAL(wp), PARAMETER :: cse_ls_irank_ref(12) = [17.0_wp, 15.5_wp, 13.0_wp, 10.5_wp, 8.0_wp, 6.0_wp, &
                                              4.5_wp,  3.5_wp,  3.0_wp,  2.8_wp, 2.7_wp, 2.7_wp]
  REAL(wp) :: cse_ls_irank(12) = cse_ls_irank_ref  ! WORKING profile (see note above); default = cse_ls_irank_ref verbatim
  REAL(wp) :: Sheath_L_max_cm = 17.0_wp  ! cm, cultivar scalar = cse_ls_irank_ref(1) (panel maximum, monotone-to-top profile); default UNCHANGED; settable from phenology_cultivar.txt row Sheath_L_max_cm (EOF-default keeps this literal); range 14-20 cm, Abichou 2016 Ch.2 Fig. 3b doc p.67 (same verbatim sentence as the compiled default above)

  ! Final INTERNODE length by inversed rank (cm). Abichou (2016) Ch.2 Fig. 3c,
  ! medians of the same panel, zero below the fifth inversed rank ("generally 4
  ! to 5 elongated internodes", Ch.2 p.66). The PEDUNCLE is not included in
  ! this profile.
  REAL(wp), PARAMETER :: cse_li_irank(12) = [16.0_wp, 11.0_wp, 7.0_wp, 4.0_wp, 1.5_wp, 0.0_wp, &
                                              0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
  REAL(wp), PARAMETER :: cse_ds1  = 0.124_wp  ! cm, sheath diameter of main-stem leaf 1; Abichou Ch.2 Table 4 doc p.71, a fitted CONSTANT of the sub-model

  ! Mean diameter of the upper sheaths (cm): the median of the eight values
  ! printed in Abichou (2016) Table S2, measured bracket 0.21-0.49 cm. The
  ! compiled literal is kept as the reference value and the working copy is
  ! settable per cultivar.
  REAL(wp), PARAMETER :: cse_ds2_ref = 0.45_wp  ! cm, mean diameter of the upper sheaths; median of the eight values printed in Table S2 doc p.87 (measured bracket 0.36-0.64, swept)
  REAL(wp) :: cse_ds2 = cse_ds2_ref           ! WORKING scalar (see note above); default = cse_ds2_ref verbatim

  ! Tillering capacity as a genotype trait (axes per plant; 0 = no cap). The
  ! measured ceiling is Abichou (2016) Ch.3 Fig. 3, where tillering ceases at a
  ! maximum axis density of 600-800 axes m-2 at 150-250 plants m-2 for his
  ! winter cultivars.
  REAL(wp) :: Axes_max_per_plant = 0.0_wp     ! axes plant-1, cultivar scalar; 0 = no cap (identity)
  REAL(wp) :: Sheath_D_cm = 0.45_wp           ! cm, cultivar scalar = cse_ds2_ref; default UNCHANGED; settable from phenology_cultivar.txt row Sheath_D_cm (EOF-default keeps this literal); range 0.36-0.64 cm, already printed below as cse_ds2_lo/hi (Table S2 doc p.87)
  REAL(wp), PARAMETER :: cse_di2  = 0.385_wp  ! cm, mean diameter of the elongated internodes; median of the ten values printed in Table S2 doc p.87 (bracket 0.21-0.49). Flat with rank — Eqn T2.6's printed shape returns a negative diameter at the flag
  REAL(wp), PARAMETER :: cse_rho1 = -0.265_wp, cse_rho2 = 1.007_wp, cse_rho3 = 0.702_wp  ! -, tiller sheath/internode diameter reduction ko; Abichou Ch.2 Table 4 doc p.71 with the glossary definition at doc p.81
  INTEGER,  PARAMETER :: cse_persist_n = 2    ! the ADEL rule: sheath n senesces with blade n+2 (Annex 2 doc p.170). MODEL-ASSUMPTION; the n+1 and n+3 arms are emitted beside it every day
  REAL(wp), PARAMETER :: cse_adel_start = 0.4_wp, cse_adel_dur = 2.0_wp  ! phyllochrons; ADEL leaf-extension coordination, Annex 2 doc p.170 (blade+sheath starts 0.4 phyllochron after tip appearance and lasts 2 phyllochrons)
  REAL(wp) :: cse_a_ep = 15.0_wp                            ! cm2 per fertile culm, ear plus exposed peduncle; no numeric ear area is printed in the source thesis, whose ear table gives the schema only, so the 10-20 cm2 bracket is swept and reported
  REAL(wp), PARAMETER :: cse_a_ep_lo = 10.0_wp, cse_a_ep_hi = 20.0_wp
  REAL(wp), PARAMETER :: cse_a_ear_meas = 6.7_wp  ! cm2 PROJECTED per fertile culm, the EAR ALONE (no peduncle); MEASURED - Zhang et al. 2023 Plants 12(19):3482 DOI 10.3390/plants12193482 Table 1, 9.2-11.8 cm2 of "spike area" per stem on their halved-cylinder convention, inverted through their own Eqn 1; bracket 5.9-7.5. REPORTED beside A_ep, not substituted for it (see E-10)
  REAL(wp), PARAMETER :: cse_ep_win = 1.8_wp  ! phyllochrons from that axis's OWN flag ligulation to a fully emerged ear+peduncle; Abichou Ch.3 doc p.96 verbatim ("at 1333degCd (1.8 phyllochron after flag leaf ligulation)"), doc p.100 places 1333 degCd close to flowering
  REAL(wp), PARAMETER :: cse_half_cyl = 1.5707963267948966_wp  ! pi/2, silhouette L*D -> half TOTAL surface of a cylinder, the quantity a gap-fraction inversion returns (LI-2200C basis); a flat blade needs no such conversion
  ! --- S5 state (all identically zero when the flag is OFF ------------- !
  LOGICAL  :: cse_active = .FALSE.
  REAL(wp) :: cse_ls_pot(CLR_MAXRANK) = 0.0_wp  ! cm, final sheath length per main-stem rank
  REAL(wp) :: cse_ds_pot(CLR_MAXRANK) = 0.0_wp  ! cm, final sheath diameter per main-stem rank
  REAL(wp) :: cse_sh_green(0:CT_MAXTIL) = 0.0_wp ! cm2 per axis, GREEN sheath silhouette area
  REAL(wp) :: cse_sh_dead(0:CT_MAXTIL)  = 0.0_wp ! cm2 per axis, senesced sheath area, retained STANDING
  REAL(wp) :: cse_ep_axis(0:CT_MAXTIL)  = 0.0_wp ! cm2 per axis, ear + exposed peduncle
  REAL(wp) :: cse_st_axis(0:CT_MAXTIL)  = 0.0_wp ! cm2 per axis, EXPOSED internode (stem) silhouette
  REAL(wp) :: cse_sai      = 0.0_wp     ! m2 m-2, ground-area GREEN sheath index (n+2 central)
  REAL(wp) :: cse_sai_n1   = 0.0_wp     ! m2 m-2, same on the n+1 persistence arm (sensitivity)
  REAL(wp) :: cse_sai_n3   = 0.0_wp     ! m2 m-2, same on the n+3 persistence arm (sensitivity)
  REAL(wp) :: cse_sai_full = 0.0_wp     ! m2 m-2, upper-bound arm: sheath full at its OWN ligulation (no ADEL extension lag)
  REAL(wp) :: cse_sdi      = 0.0_wp     ! m2 m-2, instantaneous standing-DEAD sheath index
  REAL(wp) :: cse_sdi_pool = 0.0_wp     ! m2 m-2, the standing-dead sheath POOL (accumulator, monotone by construction — the S4 lesson)
  REAL(wp) :: cse_sdi_prev = -1.0_wp    ! previous day's instantaneous dead-sheath index
  REAL(wp) :: cse_sdi_dip_max = 0.0_wp  ! m2 m-2, largest single-day DECREASE of the instantaneous recomputation
  INTEGER  :: cse_mono_viol = 0         ! days on which the instantaneous dead-sheath index decreased
  REAL(wp) :: cse_eai      = 0.0_wp     ! m2 m-2, ear + exposed peduncle index at A_ep central
  REAL(wp) :: cse_eai_lo   = 0.0_wp     ! m2 m-2, same at the 10 cm2 bracket end
  REAL(wp) :: cse_eai_hi   = 0.0_wp     ! m2 m-2, same at the 20 cm2 bracket end
  REAL(wp) :: cse_eai_meas = 0.0_wp     ! m2 m-2, EAR ALONE on Zhang's measured projected area (no peduncle) — the independent cross-check arm
  REAL(wp) :: cse_iai      = 0.0_wp     ! m2 m-2, EXPOSED internode (stem) index — a RESULT of Li vs Ls, not an input
  REAL(wp) :: cse_nonlam   = 0.0_wp     ! m2 m-2, green non-lamina total = sheath + stem + ear
  REAL(wp) :: cse_ratio    = 0.0_wp     ! -, green non-lamina / green lamina (the 10-80 % gate)
  REAL(wp) :: cse_paig     = 0.0_wp     ! m2 m-2, GREEN plant area index on the Abichou dissection basis (silhouettes, no dead)
  REAL(wp) :: cse_pai_li   = 0.0_wp     ! m2 m-2, PAI on the LI-2200C basis (adds standing dead; cylinders on the pi/2 half-surface convention)
  REAL(wp) :: cse_paig_axis = 0.0_wp    ! cm2 per LIVE axis, the per-axis form of PAIg used for the density-free gate (ii)
  REAL(wp) :: cse_paig_flower = 0.0_wp  ! m2 m-2, PAIg recorded once at the anchored flowering date
  REAL(wp) :: cse_axis_flower = 0.0_wp  ! axes m-2 alive at that date
  REAL(wp) :: cse_ratio_flower = 0.0_wp ! -, non-lamina / lamina at that date
  REAL(wp) :: cse_pai_li_max = 0.0_wp   ! m2 m-2, seasonal maximum of the LI-2200C-basis PAI

  ! ------------------------------------------------------------------------- !
  ! The canopy becomes a carbon SINK - driver record 69, accepted at the value
  ! 1. The organ-level law replaces lai = W_green sla_from_dvs(dvs2), which is
  ! carbon-blind, by a demand for carbon,
  !   demand = dA LMA,
  ! where dA is the day increment of cumulative potential blade area and LMA the
  ! lamina mass per unit area, so that leaf area can only be built if the carbon
  ! to build it exists.
  ! LMA = 48.5 g DM m-2 is the one new constant. It is taken from the properly
  ! matched FIELD measurement: Sales, Molero, Reynolds & Carmo-Silva (2022)
  ! "Phenotypic variation in photosynthetic traits in wheat grown under field
  ! versus glasshouse conditions", J Exp Bot 73(10):3221-3237,
  ! doi:10.1093/jxb/erac096, Table 1 - field-grown wheat, two lines, two years,
  ! fully irrigated (range 47.0-53.0). Two alternatives were read and rejected
  ! as unmatched: Bos (1999) prints no lamina mass per area at all (its
  ! Appendix 2 Table 1.1 lists SLWmin 20 and SLWmax 50 g DM m-2 with no
  ! source), and Evers (2006) Ch.4 p.87 measures a glasshouse value.
  ! Stated limitation: Sales et al. measure the FLAG leaf and lower ranks are
  ! thinner; Robles-Zazueta, Molero, Pinto, Foulkes, Reynolds & Murchie (2022)
  ! Front Plant Sci 13:828451, doi:10.3389/fpls.2022.828451, Table 2, resolves
  ! the rank profile.
  ! The shade response of leaf-area expansion is demographic and measured:
  ! Evers (2006) PhD thesis, doi:10.18174/24302, Ch.3 Fig. 4 measures a 75 %
  ! shade treatment cutting tiller appearance (P < 0.001), and Ch.4 Table 2
  ! measures the density effect.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: canopy_sink_lai = 0       ! runtime flag (driver_main.txt rec 69; EOF-default 0)
  REAL(wp) :: csl_LMA = 48.5_wp         ! g DM m-2 of lamina; MEASURED, field-grown wheat flag-leaf lamina at ANTHESIS, mean of the two lines of Sales et al. 2022 J Exp Bot 73:3221 Table 1 (47.0 and 50.0 g m-2), DOI 10.1093/jxb/erac096, read at source. See item E-4 of the revision record (internal; not distributed) for the full history and for the flag-leaf rank caveat
  INTEGER,  PARAMETER :: CSL_FLN_CUT_MAX = 2   ! leaves; Evers 2006 thesis Ch.3 Fig. 4 doc p.57, DOI 10.18174/24302, read at source: measured maximum 75 %-shade response of the main-stem final leaf number (9 -> 8 in the mean, 9 -> 7 in 38 % of the D3 plants)
  ! --- S6 state (all identically zero / inert when the flag is OFF ------ !
  LOGICAL  :: csl_active = .FALSE.
  LOGICAL  :: csl_limit_active = .FALSE. ! the lamina backlog is open -> the DEMOGRAPHIC limiter is acting
  REAL(wp) :: csl_lai_sink   = 0.0_wp   ! m2 m-2, the SINK green lamina area index that becomes LAI
  REAL(wp) :: csl_pai_diag   = 0.0_wp   ! m2 m-2, PAI_diag_LI2200C_basis, emitted daily (cross-step invariant 3)
  REAL(wp) :: csl_area_built = 0.0_wp   ! m2 m-2, CUMULATIVE lamina area ever built (green + standing dead + culled); monotone
  REAL(wp) :: csl_area_prev  = 0.0_wp   ! m2 m-2, yesterday's value of the same
  REAL(wp) :: csl_dA_day     = 0.0_wp   ! m2 m-2 d-1, today's new lamina area
  REAL(wp) :: csl_demand_day = 0.0_wp   ! g DM m-2 d-1, today's lamina carbon demand INCLUDING the carried backlog
  REAL(wp) :: csl_alloc_day  = 0.0_wp   ! g DM m-2 d-1, lamina carbon actually debited today
  REAL(wp) :: csl_backlog    = 0.0_wp   ! g DM m-2, unmet lamina demand carried forward; HARD non-negative
  REAL(wp) :: csl_backlog_max = 0.0_wp  ! g DM m-2, largest backlog reached in the run
  INTEGER  :: csl_limited_days = 0      ! days on which the backlog was open
  INTEGER  :: csl_fln_cut    = 0        ! leaves removed from the emergent FLN by the demographic limiter (capped at CSL_FLN_CUT_MAX)
  INTEGER  :: csl_cohorts_aborted = 0   ! tiller cohorts aborted before their first ligulated leaf
  INTEGER  :: csl_births_suspended = 0  ! days on which a cohort birth was suppressed by the limiter
  REAL(wp) :: csl_sla_real   = 0.0_wp   ! m2 kg-1, realized SLA = LAI / W_green — REPORTED against Bos 2000 and Araus 1993 (see item E-4 of the revision record, internal; not distributed: it tests the constant, not the model)
  REAL(wp) :: csl_lai_Ncap   = 0.0_wp      ! m2 m-2, the nitrogen-capacity leaf-area basis, kept explicitly separate from the sink leaf area and printed beside it so the two cannot diverge silently
  REAL(wp) :: csl_fsum_err   = 0.0_wp   ! largest |1 - (f_stem+f_root+f_storage)| after renormalization (INVARIANT 6d; must stay at machine zero)

  ! ======================================================================== !
  ! Organ-level carbon gate, water-limited leaf expansion and organ-based LAI -
  ! driver records 82, 83 and 84, each accepted at the value 1 and each
  ! requiring records 64, 65, 66 and 69 = 1.
  ! Record 82, the carbon gate on tiller birth: the day surplus carbon
  !   C_sur(t) = MAX(0, dm_gain(t-1) - LMA dA_built(t-1))
  ! is compared with the cohort construction demand
  !   D_cohort(c) = Np LMA 1e-4 SUM(blade areas tipped within one window),
  ! and the born density is dens_born(c) = Np MIN(1, SUM_window C_sur /
  ! D_cohort(c)), with the empirical axis ceiling of Abichou (2016) Ch.3 Fig. 3
  ! (600-800 axes m-2) reported rather than imposed.
  ! Record 83, water-limited leaf elongation: the measured wheat response of
  ! leaf elongation rate to soil relative water content, Acker, Zhang, Wang,
  ! Yu & Zhou (2024) Ann Bot 134:787, doi:10.1093/aob/mcae114, Fig. 4B wheat
  ! panel, whose printed regressions LER = a + b SRWC are linear with no
  ! breakpoint. The model uses the normalized form f_w = (a + b SRWC)/(a + b),
  ! which is 1 at SRWC = 1, with a and b the mean of the two printed
  ! regressions, and SRWC_j = MIN(1, theta_j / theta_fc,j) on the field-capacity
  ! basis the source defines.
  ! Record 84, LAI from the organs: LAI = SUM_c dens(c) x green blade area(c),
  ! so that SLA = LAI / W_green becomes an output rather than an input.
  ! ------------------------------------------------------------------------- !
  INTEGER  :: tiller_carbon_gate    = 0   ! rec 82 (EOF-default 0); requires rec 64, 65, 66, 69 = 1
  INTEGER  :: leaf_expansion_water  = 0   ! rec 83 (EOF-default 0); requires rec 64, 65 = 1
  INTEGER  :: lai_from_organs       = 0   ! rec 84 (EOF-default 0); requires rec 64-69 = 1
  ! --- rec 82 state (all identically zero when OFF ------------------------- !
  REAL(wp) :: csl_dm_gain_day = 0.0_wp            ! g DM m-2 d-1, the day's dm_gain after the S6 split (set only when rec 69 is ON)
  REAL(wp) :: tcg_csur_day    = 0.0_wp            ! g DM m-2 d-1, max(0, dm_gain - LMA*dA_built) of the PREVIOUS day, used today
  REAL(wp) :: tcg_sup(CT_MAXTIL)   = 0.0_wp       ! g DM m-2, SUM of C_sur over the cohort's window of opportunity
  REAL(wp) :: tcg_dem(CT_MAXTIL)   = 0.0_wp       ! g DM m-2, D_cohort: organ-law carbon cost of one window of the new cohort's leaves
  REAL(wp) :: tcg_ratio(CT_MAXTIL) = 0.0_wp       ! -, min(1, sup/dem) at the birth decision
  REAL(wp) :: tcg_nwin(CT_MAXTIL)  = 0.0_wp       ! days in the window
  INTEGER  :: tcg_nrank(CT_MAXTIL) = 0            ! ranks counted in D_cohort
  INTEGER  :: tcg_env600_days = 0, tcg_env800_days = 0   ! days with live axes above the printed Abichou envelope bounds — REPORTED, nothing acts on them
  REAL(wp), PARAMETER :: tcg_axes_env_lo = 600.0_wp, tcg_axes_env_hi = 800.0_wp    ! axes m-2; Abichou 2016 Ch.3 Fig. 3 printed range, used as a reported envelope rather than a cap when record 82 is on
  REAL(wp), PARAMETER :: tcg_die_frac_lo = 0.57_wp, tcg_die_frac_hi = 0.67_wp   ! -, Evers 2006 thesis Ch.4 measured die fraction at 100 / 262 plants m-2 (via the tiller-environment design note (internal; not distributed)) — ENVELOPE, reported, NOT applied
  REAL(wp) :: tcg_ne_env_lo = 0.0_wp, tcg_ne_env_hi = 0.0_wp   ! ears m-2 the Evers band would imply from Namax — diagnostic only; the productive fraction is NOT DETERMINED
  ! --- rec 83 state (all identically zero / unity when OFF ------------------ !
  REAL(wp), PARAMETER :: lew_a200 = 1.6_wp, lew_b200 = 3.4_wp   ! Acker 2024 Fig. 4B wheat, 200 ppm: LER = 1.6 + 3.4*SRWC (cm tiller-1 d-1)
  REAL(wp), PARAMETER :: lew_a800 = 0.8_wp, lew_b800 = 4.7_wp   ! Acker 2024 Fig. 4B wheat, 800 ppm: LER = 0.8 + 4.7*SRWC
  REAL(wp) :: lew_srwc = 1.0_wp        ! -, root-fraction-weighted SRWC = theta/theta_fc (Acker Eq. 1 basis)
  REAL(wp) :: lew_se_rw = 1.0_wp       ! -, root-fraction-weighted effective saturation (theta-theta_r)/(theta_s-theta_r) — ALTERNATIVE basis, reported only
  REAL(wp) :: lew_theta_rw = 0.0_wp, lew_thetafc_rw = 0.0_wp   ! m3 m-3, root-weighted theta and theta_fc
  REAL(wp) :: lew_fw = 1.0_wp, lew_fw_200 = 1.0_wp, lew_fw_800 = 1.0_wp  ! -, expansion factor: mean line (APPLIED), and the two printed lines (bracket, reported)
  REAL(wp) :: lew_f(0:CT_MAXTIL, CLR_MAXRANK)    = 0.0_wp   ! -, accumulated realized fraction of A_pot under the arrest, per axis class and rank
  REAL(wp) :: lew_fclk(0:CT_MAXTIL, CLR_MAXRANK) = 0.0_wp   ! -, previous-day CLOCK fraction of the same rank (the S2/S3 linear LED ramp)
  REAL(wp) :: lew_area_axis(0:CT_MAXTIL) = 0.0_wp  ! cm2 per axis, arrested realized blade area (all ranks, before senescence)
  REAL(wp) :: lew_bai = 0.0_wp, lew_bai_clock = 0.0_wp   ! m2 m-2, arrested vs clock realized blade index over live axes (diagnostic pair)

  ! ======================================================================== !
  ! Tiller survival and fertile ear number - driver record 99, accepted at the
  ! value 4.
  ! The measured parts of the law:
  ! - the mortality WINDOW and the death ORDER are Abichou (2016) PhD thesis
  !   (NNT tel-02316154): the mortality phase and its normalized clock taun are
  !   defined on doc p.45, and p.29 states verbatim that the tillers that die
  !   are the last emerged, least developed ones, so the cull is applied
  !   youngest first; p.20-21 state that a proportion of tillers dies and that
  !   the rate of mortality is set by the competition for light.
  ! - the SIGNAL is the shortfall of assimilate supply against the axes' own
  !   construction demand. Light capture per cohort is Beer-Lambert on the
  !   cohorts above it,
  !     I_c = exp(-ct_kL SUM_{HS(c') > HS(c)} dens(c') A_real(c') 1e-4),
  !     w_c = dens(c) A_real(c) 1e-4 I_c,  S_c = dm_gain w_c / SUM w,
  !   so the shares sum to one and carbon is conserved; the demand is the
  !   Abichou blade geometry times the Sales et al. (2022) lamina mass per area.
  ! - the kill fraction is scaled by the normalized clock, kill(c) = taun
  !   (1 - MIN(1, SUM S_axis(c)/SUM D_axis(c))), so that nothing dies at
  !   taun = 0 (Abichou doc p.45), and it is monotone and irreversible.
  ! Falsification target for the resulting ear number, never a tuning target:
  ! Fischer, R. A. (1985) J Agric Sci Camb 105:447-461,
  ! doi:10.1017/S0021859600056495, whose measured spike densities are 436
  ! spikes m-2 (p.450, cv. Yecora 70), 598 control and 553 shaded at anthesis
  ! (Table 4 p.451, cv. Egret), and 319-470 across the eight Table 5 treatments
  ! (p.454).
  ! ------------------------------------------------------------------------- !
  INTEGER  :: tiller_survival_mech = 0   ! rec 99 (EOF-default 0); requires rec 64, 65, 66, 69 = 1
  REAL(wp) :: tsv_sup(0:CT_MAXTIL)  = 0.0_wp   ! g DM axis-1, cumulative own light-shared assimilate over the mortality window
  REAL(wp) :: tsv_dem(0:CT_MAXTIL)  = 0.0_wp   ! g DM axis-1, cumulative own potential lamina carbon over the same window
  REAL(wp) :: tsv_kill(0:CT_MAXTIL) = 0.0_wp   ! -, realized death fraction of the class (monotone, irreversible)
  REAL(wp) :: tsv_ilight(0:CT_MAXTIL) = 1.0_wp ! -, Beer light fraction reaching the class's own leaves
  REAL(wp) :: tsv_apot_prev(0:CT_MAXTIL) = 0.0_wp ! cm2 axis-1, yesterday's REALIZED blade area (today's increment is the axis's own lamina demand)
  REAL(wp) :: tsv_dens_pre(0:CT_MAXTIL)  = 0.0_wp ! axes m-2 entering the mortality window
  REAL(wp) :: tsv_taun = 0.0_wp                ! -, Abichou's normalized mortality clock, 0 at tau_sm and 1 at tau_em (doc p.45)
  INTEGER  :: tsv_days = 0                     ! days the mortality window has been open
  LOGICAL  :: tsv_frozen = .FALSE.             ! .TRUE. once tau >= tau_em: the fertile population is final

  ! Closure accumulators for the record-99 cull: the same green blade area
  ! removed by the cull, counted on both sides so the two must agree.
  REAL(wp) :: tsv_cullgreen_seen    = 0.0_wp   ! m2 m-2, green blade area removed by the rec-99 cull
  REAL(wp) :: tsv_cullgreen_charged = 0.0_wp   ! m2 m-2, the same area as charged to the leaf-death coupling

  ! ======================================================================== !
  ! Tiller carbon affordability - driver record 99 mode 2. The three sources
  ! divide cleanly: the source-sink PATTERN is the MIN(supply, demand) form
  ! printed by Yin & Schapendonk (2004) NJAS 51-4, Table 1 ("min(n_sup,
  ! n_dem)"); the mortality WINDOW and the death ORDER are Abichou (2016); and
  ! the ear-number envelope is Fischer (1985).
  ! Light shares telescope exactly:
  !   I_top(c) = exp(-k Lab_c);   f_c = I_top(c) (1 - exp(-k L_c));
  !   F_tot = SUM f_c = 1 - exp(-k L_tot);   s_c = f_c / F_tot,
  ! so SUM s_c = 1 identically and the un-intercepted fraction exp(-k L_tot) is
  ! not allocated to any axis.
  ! The supply currency is the growth supply of allocate_daily_biomass,
  !   MAX(0, net_assim_ch2o - rm_day)  (g CH2O m-2 d-1),
  ! i.e. canopy net photosynthesis converted at 30/44 minus whole-crop
  ! maintenance respiration (Teh 2006 Eqs. 7.3-7.4) - net of maintenance and
  ! gross of growth respiration, which is the currency Teh Eq. 7.7 converts to
  ! dry matter with the growth-conversion efficiency. The demand is
  ! dA_c = MAX(0, A_c(t) - A_c(t-1)) on the carbon-unlimited potential area
  ! times the lamina mass per area of Sales et al. (2022) J Exp Bot 73:3221,
  ! Table 1. Jiang, Mulero, Bonfil & Helman (2022) Plant Cell Environ,
  ! doi:10.1111/pce.14430, Table 7 is a falsification band only.
  REAL(wp) :: tca_gsup_day  = 0.0_wp           ! g CH2O m-2 d-1, TODAY's growth_supply, exported by allocate_daily_biomass (mode 2 only)
  REAL(wp) :: tca_gsup_prev = 0.0_wp           ! g CH2O m-2 d-1, G(t-1): the previous day's growth_supply, the supply shared out
  REAL(wp) :: tca_S(0:CT_MAXTIL)   = 0.0_wp    ! g CH2O axis-1, cumulative own light-shared supply over the window
  REAL(wp) :: tca_D(0:CT_MAXTIL)   = 0.0_wp    ! g CH2O axis-1, cumulative own lamina demand (dA x LMA x g_leaf) over the window
  REAL(wp) :: tca_phi(0:CT_MAXTIL) = 1.0_wp    ! -, affordability MIN(1, tca_S/tca_D); 1 = fully affordable
  REAL(wp) :: tca_share(0:CT_MAXTIL) = 0.0_wp  ! -, today's light share s_c (SUM = 1)
  REAL(wp) :: tca_ftot = 0.0_wp                ! -, F_tot = SUM f_c, the absorbed PAR fraction
  REAL(wp) :: tca_dens_pre(0:CT_MAXTIL)  = 0.0_wp  ! axes m-2 entering the window, captured per class at its first day with density
  REAL(wp) :: tca_area_prev(0:CT_MAXTIL) = 0.0_wp  ! cm2 axis-1, yesterday's organ-clock blade area (today's increment is the demand)
  REAL(wp) :: tca_unaff = 0.0_wp               ! axes m-2, the monotone population target MAX over days of n_unaff
  REAL(wp) :: tca_dead_c(0:CT_MAXTIL) = 0.0_wp ! axes m-2, PER-COHORT dead density, monotone and irreversible (Dev-10: a cohort's
                                               ! density never rises back once a younger cohort has been culled; P2 monotone)
  REAL(wp) :: tca_dead = 0.0_wp                ! axes m-2, SUM_c tca_dead_c: the dead density actually taken (>= tca_unaff)
  INTEGER  :: tca_days = 0                     ! window days stepped
  INTEGER  :: tca_days_skip = 0                ! window days skipped because L_tot = 0 (expected 0)
  LOGICAL  :: tca_frozen = .FALSE.             ! .TRUE. once tau >= tau_em: the fertile population is final
  ! reporting state for output/daily_tiller_afford.csv (written by the block, read by main only)
  REAL(wp) :: tca_hs(0:CT_MAXTIL)  = 0.0_wp    ! Haun stage per cohort (the vertical-order key)
  REAL(wp) :: tca_L(0:CT_MAXTIL)   = 0.0_wp    ! m2 m-2, leaf area index of the cohort
  REAL(wp) :: tca_sax(0:CT_MAXTIL) = 0.0_wp    ! g CH2O axis-1 d-1, today's S_ax
  REAL(wp) :: tca_dax(0:CT_MAXTIL) = 0.0_wp    ! g CH2O axis-1 d-1, today's D_ax
  REAL(wp) :: tca_ltot = 0.0_wp                ! m2 m-2, L_tot
  REAL(wp) :: tca_sum_s = 0.0_wp               ! -, SUM s_c as computed (identity check, 1 to machine precision)
  REAL(wp) :: tca_ident = 0.0_wp               ! g CH2O m-2 d-1, SUM_c S_c - G(t-1) as computed (identity check)
  REAL(wp) :: tca_taun = 0.0_wp                ! -, Abichou's normalized clock, REPORTED only (no ramp is applied)

  ! ======================================================================== !
  ! Tiller carbon affordability with source-timed survival - driver record 99
  ! mode 4, the released value. The normalized mortality clock and its shape
  ! parameters are Abichou (2016) Ch.1 Eqn T4.4, whose fitted theta = 3.32 and
  ! delta = 1.77 are printed with an RMSE of 0.10 (doc p.36).
  REAL(wp), PARAMETER :: tca4_theta = 3.32_wp    ! -, Abichou 2016 Ch.1 Eqn T4.4 theta; value doc p.36 / PDF p.53 ("3.32 and 1.77 with an RMSE of 0.10")
  REAL(wp), PARAMETER :: tca4_delta = 1.77_wp    ! -, Abichou 2016 Ch.1 Eqn T4.4 delta; value doc p.36 / PDF p.53 (same sentence); Eqn doc p.30 / PDF p.47
  REAL(wp) :: tca4_U = 0.0_wp                  ! axes m-2, U(t) = SUM_{c>=1} dens_pre (1 - phi_c): running estimate of the eventual deaths
  REAL(wp) :: tca4_F = 0.0_wp                  ! -, F(tau_n) = 1 - EXP(-theta tau_n/(delta - tau_n)); := 1 on the endpoint day
  REAL(wp) :: tca4_target = 0.0_wp             ! axes m-2, F(tau_n) U(t): today's capped death target (tca_unaff = MAX over days)
  LOGICAL  :: tca4_endpoint_day = .FALSE.      ! .TRUE. from the endpoint day on (the endpoint row of daily_tiller_afford4.csv)

  ! ======================================================================== !
  ! Stem structural demand - driver record 100. This version accepts the value
  ! 0, so the incumbent stem partition runs and the demand formulation is not
  ! reachable. Its sources are kept because the constants are shared with the
  ! height and internode records: the peduncle to penultimate-internode LENGTH
  ! ratio is measured by Rivera-Amado, Trujillo-Negrellos, Molero, Reynolds,
  ! Sylvester-Bradley & Foulkes (2019) Field Crops Research 240:154-167,
  ! doi:10.1016/j.fcr.2019.04.016, and the stem reserve concentration bracket
  ! by Cereal Research Communications 40(3):413-422,
  ! doi:10.1556/CRC.40.2012.3.10.
  INTEGER  :: stem_demand_mech = 0                ! rec 100 (EOF-default 0); 0 = OFF, 1 = the demand formulation, 2 = an identity diagnostic on the incumbent physics, 3 = mode 1 with the Eq. (6) height basis corrected, the peduncle simulated and height carbon-limited, 4 = a per-axis, per-degree-day, bounded structural stem-mass demand on printed rates and durations, with height as a diagnostic only (records 103 and 104)
  REAL(wp) :: sdm_dens      = 0.0_wp      ! plants m-2, set once from the driver's plant density
  REAL(wp) :: sdm_h_pot     = 0.0_wp      ! cm, today's potential telescoped collar height
  REAL(wp) :: sdm_h_pot_prev= 0.0_wp      ! cm, yesterday's

  ! Diagnostic state read off the model's own trajectory; none of it is a
  ! parameter.
  REAL(wp) :: sdm_dvs2_flag  = -1.0_wp    ! dvs2 on the flag-leaf-ligule day (mode 3); < 0 = peduncle window not open
  REAL(wp) :: sdm_f_ped      = 0.0_wp     ! -, fraction of the potential peduncle extension elapsed today
  REAL(wp) :: sdm_h_collar   = 0.0_wp     ! cm, the collar limb of the potential height (diagnostic)
  REAL(wp) :: sdm_Lped_pot   = 0.0_wp     ! cm, ARCHITECTURAL potential peduncle length = ratio x the model's own top internode
  REAL(wp) :: sdm_h_real     = 0.0_wp     ! cm, REALIZED whole-plant structural height -- CARBON-LIMITED, EMERGENT
  REAL(wp) :: sdm_Lped_emerg = 0.0_wp     ! cm, EMERGENT peduncle length = h_real - h_collar (>= 0)
  REAL(wp) :: sdm_dh_day     = 0.0_wp     ! cm d-1, today's height debt fed to the Eq.(6) inverse (mode 3)
  REAL(wp) :: sdm_frac_real  = 0.0_wp     ! -, realized / potential stem increment today (mode 3)
  REAL(wp) :: sdm_rm_stem_day = 0.0_wp    ! g CH2O m-2 d-1, the stem's OWN maintenance respiration
  REAL(wp) :: sdm_pot_dw_day  = 0.0_wp    ! g DM m-2 d-1, structural growth demand
  REAL(wp) :: sdm_dem_day     = 0.0_wp    ! g CH2O m-2 d-1, D_stem
  REAL(wp) :: sdm_avail_day   = 0.0_wp    ! g CH2O m-2 d-1, avail
  REAL(wp) :: sdm_take_day    = 0.0_wp    ! g CH2O m-2 d-1, take
  REAL(wp) :: sdm_dw_stem_day = 0.0_wp    ! g DM m-2 d-1, realized stem increment
  REAL(wp) :: sdm_legacy_dw_stem_day = 0.0_wp  ! g DM m-2 d-1, what the fraction table would have written
  REAL(wp) :: sdm_refused_day = 0.0_wp    ! g CH2O m-2 d-1, refused surplus (enters no organ)
  REAL(wp) :: sdm_refused_cum = 0.0_wp    ! g CH2O m-2, season total
  REAL(wp) :: sdm_dem_cum     = 0.0_wp, sdm_got_cum = 0.0_wp
  REAL(wp) :: sdm_ident_day   = 0.0_wp    ! g CH2O m-2 d-1, growth_supply - SUM g_i dW_i == THE IDENTITY GAP
  REAL(wp) :: sdm_ident_legacy_day = 0.0_wp  ! the same quantity on the incumbent increments of the same day
  REAL(wp) :: sdm_ident_max   = 0.0_wp, sdm_ident_legacy_max = 0.0_wp
  REAL(wp) :: sdm_close_max   = 0.0_wp    ! max |growth_supply - SUM g_i dW_i - refused| : the Gate-D closure of this path
  INTEGER  :: sdm_days_dem = 0, sdm_days_sup = 0, sdm_days_ref = 0  ! days demand-bound / supply-bound / with refusal

  ! ======================================================================== !
  ! Stem demand brackets - driver records 103 and 104, both accepted at the
  ! value 0, so neither bracket is selected in this version. The published
  ! brackets they would select are the internode elongation onset of Fournier,
  ! Andrieu, Ljutovac & Saint-Jean (2003), Proc. PMA03, and the internode
  ! duration and lineic mass of Gauthier, Barillot, Andrieu et al. (2020)
  ! J Exp Bot, Supplementary Data S1.
  INTEGER  :: sdm4_dur_bracket  = 0        ! rec 103 (EOF-default 0 = not selected); 1/2/3 = onset F with Delta_I lo/mid/hi, 4/5/6 = onset M with lo/mid/hi
  INTEGER  :: sdm4_rate_bracket = 0        ! rec 104 (EOF-default 0 = not selected); 1 = lo, 2 = mid, 3 = hi
  REAL(wp), PARAMETER :: sdm4_dur_lo   = 277.0_wp          ! degCd, Kirby 1988 FCR 18:127 folio 135, internodes 6 and 7 "about 277 degCd"
  REAL(wp), PARAMETER :: sdm4_dur_mid  = 4.0_wp * 87.7_wp  ! degCd, Kang 2008 folio 1104 (4 GC) x folio 1103 (87.7 degCd GC-1)
  REAL(wp), PARAMETER :: sdm4_dur_hi   = 476.0_wp          ! degCd, Kirby 1988 folio 135, peduncle
  REAL(wp), PARAMETER :: sdm4_wsc_struct = 0.7_wp          ! -, 1 - 0.30: Kirby 1988 folio 137 "up to 30% of its dry mass as soluble carbohydrate"; a model construction
  REAL(wp), PARAMETER :: sdm4_rate_I9   = 1.0_wp           ! mg DM degCd-1 internode-1, internode 9; Kirby 1988 folio 139 Fig. 8 caption
  REAL(wp), PARAMETER :: sdm4_rate_I8   = 1.7_wp           ! mg DM degCd-1 internode-1, internode 8 "at its most rapid"; Kirby 1988 folio 139 Fig. 8 caption / folio 140
  REAL(wp), PARAMETER :: sdm4_rate_stem = 3.8_wp           ! mg DM degCd-1 axis-1, "maximum growth rate of the whole stem"; Kirby 1988 folio 139; applied as a cap
  REAL(wp) :: sdm4_dur_Cd  = 0.0_wp        ! degCd, Delta_I of the selected rec-103 member
  INTEGER  :: sdm4_onset_k = 0             ! 0 = onset member F (Fournier collar rule), 2 = onset member M (Malvoisin (i-2).0 rule)
  REAL(wp) :: sdm4_r_mg    = 0.0_wp        ! mg DM degCd-1 internode-1, r_I of the selected rec-104 member
  REAL(wp) :: sdm4_cap_mg  = 0.0_wp        ! mg DM degCd-1 axis-1, r_cap of the selected rec-104 member
  REAL(wp) :: sdm4_tau_on(0:CT_MAXTIL, CLR_MAXRANK+1) = -1.0_wp  ! degCd, latched onset per axis class and rank (slot CLR_MAXRANK+1 = the peduncle; -1 = not yet)
  REAL(wp) :: sdm4_tau_se  = -1.0_wp       ! degCd, latched tau_SE (= ct_tau_sm on the day it turns >= 0)
  INTEGER  :: sdm4_x_first(0:CT_MAXTIL) = 0 ! first elongating rank per class (0 = none yet; diagnostic)
  REAL(wp) :: sdm4_pot_dw_today = 0.0_wp   ! g DM m-2 d-1, today's ceiling demand (the ONE value the mode-4 arm reads)
  INTEGER  :: sdm4_nl_ms       = 0         ! main-stem final leaf number used today (diagnostic)
  INTEGER  :: sdm4_n_open_ms   = 0         ! main-stem ranks (incl. peduncle) with an OPEN window at end of day
  REAL(wp) :: sdm4_n_open_wtd  = 0.0_wp    ! density-weighted open ranks per live axis, all classes
  REAL(wp) :: sdm4_dem_ms_mg   = 0.0_wp    ! mg DM axis-1 d-1, main-stem demand BEFORE the cap
  INTEGER  :: sdm4_cap_bound_today = 0     ! 1 if the whole-stem cap bound on any live class today
  INTEGER  :: sdm4_days_cap    = 0         ! season count of cap-bound days
  REAL(wp) :: sdm4_post_anth_dem_cum = 0.0_wp  ! g DM m-2, demand accumulated on days with dvs2 >= 1 (the top-two-rank residue)
  REAL(wp) :: sdm4_t_close_last = -1.0_wp  ! degCd, latest window close among the main stem's elongating ranks
  REAL(wp) :: sdm4_t_close_subtop = -1.0_wp ! degCd, latest window close among main-stem ranks BELOW the top internode (R5 rank clause)
  REAL(wp) :: sdm4_t_close_all = -1.0_wp   ! degCd, latest window close over ALL live axis classes (tiller windows close after the main stem's; R5 crop-wide clause)
  REAL(wp) :: sdm4_tau_anth    = -1.0_wp   ! degCd, tau on the anthesis day (first dvs2 >= 1; diagnostic)

  ! ======================================================================== !
  ! Fate of the surplus (refused) carbon - driver record 101, accepted at the
  ! value 0, so the incumbent fate runs. The law the record would select rests
  ! on Schnyder, H. (1993) "The role of carbohydrate storage and redistribution
  ! in the source-sink relations of wheat and barley during grain filling - a
  ! review", New Phytol 123:233-245,
  ! doi:10.1111/j.1469-8137.1993.tb03731.x, whose p.238 states the measured
  ! ceiling on stem water-soluble carbohydrate as a fraction of internode dry
  ! weight in wheat and barley, and reports that the store is increased by
  ! elevated CO2 (Dubois et al. 1990; Winzeler et al. 1990) and that ear removal
  ! after anthesis does not reduce assimilation (Apel et al. 1973).
  INTEGER  :: surplus_fate_mech = 0        ! rec 101 (EOF-default 0); 0 = OFF, 1 = the law, 2 = deposition closed at anthesis
  REAL(wp), PARAMETER :: sfl_g_wsc = 1.242_wp  ! g CH2O g-1 DM, CARBOHYDRATES; Teh 2006 Table 7.2 p.147 (Penning de Vries et al. 1983)
  REAL(wp) :: sfl_dw_wsc_day    = 0.0_wp   ! g DM m-2 d-1, refused carbon STORED as stem non-structural DM today
  REAL(wp) :: sfl_resid_day     = 0.0_wp   ! g CH2O m-2 d-1, REPORTED residual: refused carbon the capacity would not take
  REAL(wp) :: sfl_want_day      = 0.0_wp   ! g DM m-2 d-1, what today's refusal would have stored with no ceiling
  REAL(wp) :: sfl_room_day      = 0.0_wp   ! g DM m-2 d-1, today's closed-form capacity room
  REAL(wp) :: sfl_stored_dm_cum = 0.0_wp   ! g DM m-2, season total stored
  REAL(wp) :: sfl_stored_ch2o_cum = 0.0_wp ! g CH2O m-2, season total ACCOMMODATED (= sfl_g_wsc * stored DM)
  REAL(wp) :: sfl_resid_cum     = 0.0_wp   ! g CH2O m-2, season total REPORTED RESIDUAL
  REAL(wp) :: sfl_resid_pre     = 0.0_wp   ! g CH2O m-2, residual accrued at dvs2 < 1
  REAL(wp) :: sfl_resid_post    = 0.0_wp   ! g CH2O m-2, residual accrued at dvs2 >= 1
  REAL(wp) :: sfl_pool_ratio_max = 0.0_wp  ! max W_reserve / W_stem over the run (falsification metric, decl. (i))
  REAL(wp) :: sfl_ledger_day    = 0.0_wp   ! g CH2O m-2 d-1, supply - allocated - stored - respired - residual (EXACT BY CONSTRUCTION; never cited as closure)
  REAL(wp) :: sfl_ledger_max    = 0.0_wp   ! max |sfl_ledger_day|
  INTEGER  :: sfl_days_cap = 0, sfl_days_win = 0, sfl_days_store = 0  ! days capacity-bound / window-closed (mode 2) / with storage

  ! ======================================================================== !
  ! True-stem / leaf-sheath separation - driver record 102, accepted at the
  ! value 0, so the lumped stem pool runs. The separation the record would
  ! apply uses the sheath partitioning index measured by Rivera-Amado et al.
  ! (2019) Field Crops Research 240:154-167, doi:10.1016/j.fcr.2019.04.016
  ! (printed range 0.16 to 0.20).
  INTEGER  :: stem_sheath_sep_mode = 0     ! rec 102 (EOF-default 0); 0 = OFF, 1 = ShPI 0.16 (printed lower end), 2 = ShPI 0.20 (printed upper end)
  REAL(wp), PARAMETER :: ssp_shpi_lo = 0.16_wp   ! leaf-sheath DM partitioning index, printed LOWER end of the 26-cultivar range; Rivera-Amado 2019 folio 158
  REAL(wp), PARAMETER :: ssp_shpi_hi = 0.20_wp   ! the same range's printed UPPER end; folio 158
  INTEGER,  PARAMETER :: ssp_a7_offset_d = 7     ! d after anthesis: the source's ONE dissection stage is GS65 + 7 days (folio 156)
  REAL(wp) :: ssp_shpi        = 0.0_wp     ! the ShPI in force this run (set once from the mode)
  REAL(wp) :: W_ts            = 0.0_wp     ! g DM m-2, STRUCTURAL TRUE STEM (internodes + peduncle): W_stem - W_reserve - W_sheath_est. Read-only, DERIVED.
  REAL(wp) :: W_sheath_est    = 0.0_wp     ! g DM m-2, leaf-sheath ESTIMATE = MIN(ssp_shpi * AGDM, W_stem - W_reserve). Read-only, DERIVED.
  REAL(wp) :: ssp_ident_day   = 0.0_wp     ! g DM m-2, (W_ts + W_sheath_est + W_reserve) - W_stem today
  REAL(wp) :: ssp_ident_max   = 0.0_wp     ! max |ssp_ident_day| over the run
  REAL(wp) :: ssp_wstem_max   = 0.0_wp     ! max W_stem over the run (the G2 tolerance scale)
  INTEGER  :: ssp_cap_today   = 0          ! 1 when ssp_shpi * AGDM exceeded the pool and the cap bound
  INTEGER  :: ssp_cap_days    = 0          ! days on which the cap bound
  INTEGER  :: ssp_anth_day    = -1         ! model day of the first dvs2 >= 1.0 (-1 = not yet)
  INTEGER  :: ssp_a7_day      = -1         ! model day of the like-for-like snapshot (ssp_anth_day + ssp_a7_offset_d)
  LOGICAL  :: ssp_a7_done     = .FALSE.
  REAL(wp) :: ssp_ts_mg_axis  = 0.0_wp     ! mg per live axis, W_ts
  REAL(wp) :: ssp_sw_struct   = 0.0_wp     ! mg cm-1 per live axis, W_ts basis
  REAL(wp) :: ssp_sw_total    = 0.0_wp     ! mg cm-1 per live axis, (W_ts + W_reserve) basis = the source's total-DM basis
  REAL(wp) :: ssp_stpi_struct = 0.0_wp     ! W_ts / AGDM
  REAL(wp) :: ssp_stpi_total  = 0.0_wp     ! (W_ts + W_reserve) / AGDM -- compare printed StPI 0.32-0.41 (folio 158), FALSIFICATION only
  REAL(wp) :: ssp_shpi_eff    = 0.0_wp     ! W_sheath_est / AGDM (= ssp_shpi unless the cap bound)
  REAL(wp) :: ssp_share       = 0.0_wp     ! W_sheath_est / (W_sheath_est + W_ts + W_reserve) -- compare the DERIVED band 0.281-0.372, FALSIFICATION only
  REAL(wp) :: ssp_a7_wts = 0.0_wp, ssp_a7_wsh = 0.0_wp, ssp_a7_wres = 0.0_wp, ssp_a7_wstem = 0.0_wp   ! the A7 snapshot, g DM m-2
  REAL(wp) :: ssp_a7_agb = 0.0_wp, ssp_a7_h_cm = 0.0_wp, ssp_a7_axes = 0.0_wp                          ! g DM m-2, cm, axes m-2
  REAL(wp) :: ssp_a7_sw_struct = 0.0_wp, ssp_a7_sw_total = 0.0_wp                                     ! mg cm-1 per live axis
  ! --- rec 84 state ------------------------------------------------------- !
  REAL(wp) :: lfo_sla_out = 0.0_wp     ! m2 kg-1, SLA = LAI / W_green as an OUTPUT (identical to csl_sla_real; emitted to daily_growth.csv when rec 84 is ON)

  ! ---- Record 84 leaf-death coupling ---- !
  ! The nitrogen that leaves with a dead lamina is the structural remainder
  ! only: Bertheloot, Martre & Andrieu (2008) Plant Physiol 148:1707-1720,
  ! doi:10.1104/pp.108.124156, Results p.1711, and Bertheloot, Martre & Andrieu
  ! (2008) Funct Plant Biol 35(10):781-796, doi:10.1071/FP08064, p.782, where
  ! "it was assumed that all organs share a" culm-common mobile pool, so that
  ! the exportable fraction has already left before the tissue is written off.
  REAL(wp) :: cls_sl_prev(0:CT_MAXTIL) = 0.0_wp  ! yesterday's senesced-leaf number Sl per axis class (rec 84 only; the death increment is dSl by rank)
  REAL(wp) :: lfo_gai_died_day = 0.0_wp      ! m2 m-2 d-1, green blade area index lost today under the organ law (dead-fraction increments on live axes plus the green area of culled axes)
  REAL(wp) :: lfo_fdied_day    = 0.0_wp   ! -, died / (green_after + died) = fraction of today's pre-retirement green area that died; the fraction applied to W_green and N_leaf
  REAL(wp) :: lfo_sen_dm_day = 0.0_wp, lfo_sen_N_out_day = 0.0_wp, lfo_sen_N_dead_day = 0.0_wp, lfo_sen_N_resorb_day = 0.0_wp  ! today's shed lamina DM (g m-2), N out of N_leaf, N to the dead leaf, N resorbed to N_mobile_veg (g N m-2)
  REAL(wp) :: lfo_sen_dm_cum = 0.0_wp, lfo_sen_N_out_cum = 0.0_wp, lfo_sen_N_dead_cum = 0.0_wp, lfo_sen_N_resorb_cum = 0.0_wp  ! season cumulatives of the four

  ! ======================================================================== !
  ! Post-anthesis leaf-nitrogen capacity route - driver record 87, accepted at
  ! the value 1. Nitrogen in excess of the lamina holding capacity (the maximum
  ! measured field-wheat Narea of Silva-Perez et al. 2020, J Exp Bot
  ! 71(7):2299-2311, doi:10.1093/jxb/erz439, Fig. 4B p.2305) is routed to the
  ! culm-common mobile pool rather than retained in the lamina, which is the
  ! meta-process of Bertheloot, Martre & Andrieu (2008) Funct Plant Biol
  ! 35:781-796, doi:10.1071/FP08064.
  INTEGER  :: leafN_cap_route = 0     ! runtime flag (driver_main.txt rec 87; EOF-default 0)
  REAL(wp) :: lnc_lai_basis = 0.0_wp          ! m2 m-2, the capacity basis used today (= the LAI assim read)
  REAL(wp) :: lnc_excess_day = 0.0_wp         ! g N m-2 d-1, standing N above sln_max_N*LAI exported to N_mobile_veg
  REAL(wp) :: lnc_to_mobile_day = 0.0_wp      ! g N m-2 d-1, surplus inflow the lamina could not hold, routed to N_mobile_veg

  ! ======================================================================== !
  ! Sowing-date anchor - driver record 88. The deck's own sowing date is used as
  ! the anchor of the pre-emergence clock; at 0 the anchor is the compiled
  ! default of day 14.
  INTEGER  :: sowing_model_day = 0        ! runtime flag (driver_main.txt rec 88; EOF-default 0 = the compiled default day 14)
  INTEGER  :: iD_sowing        = 14   ! model day of SOWING = DAS 0; set from rec 88 when that is > 0

  ! ======================================================================== !
  ! Grain-nitrogen sink - driver record 86, accepted at the value 1.
  ! The grain nitrogen sink and its supply are those of Bertheloot, Martre &
  ! Andrieu (2008) Funct Plant Biol 35:781-796, doi:10.1071/FP08064, p.785:
  !   Eq. (11)  dN_g,tot = MIN(S_g,tot_pot, N_c,mob),  the mobile pool being
  !             the only source of nitrogen for the grains,
  !   Eq. (12)  S_g,tot_pot(t) = gamma N_g,tot_pot(t) up to the plateau time
  !             tau, constant thereafter,
  !   Eq. (6)   the lamina synthesis term, restored here so the sink and the
  !             synthesis compete for the same mobile pool,
  ! with the structural floors of Bertheloot, Martre & Andrieu (2008) Plant
  ! Physiol 148:1707-1720, doi:10.1104/pp.108.124156, p.1711.
  ! The per-kernel grain nitrogen envelope is a falsification target, never a
  ! tuning target: Martre, Porter, Jamieson & Triboi (2003) Plant Physiol
  ! 133:1959-1967, doi:10.1104/pp.103.030585, Table I p.1961, whose measured
  ! control values span 0.75 to 1.38 mg N per kernel. The duration constant,
  ! 250 degree-Cd of grain filling at which nitrogen import ceases, is the value
  ! Martre et al. report from earlier wheat measurements.
  INTEGER  :: grainN_sink_mech = 0                    ! runtime flag (driver_main.txt rec 86; EOF-default 0)
  REAL(wp), PARAMETER :: gnf_ng0_per_ear = 2.4e-3_wp  ! g N ear-1 at anthesis; Bertheloot 2008 FPB Table 1 p.784 (2.4e-6 kg) + p.786 "set at 2.4 mg (Girard 1997)"
  REAL(wp), PARAMETER :: gnf_g_rate      = 0.005_wp   ! [degCd]-1 relative rate of POTENTIAL grain N accumulation; FPB Table 1 p.784. DERIVED, NOT MEASURED
  REAL(wp), PARAMETER :: gnf_tau_degCd   = 250.0_wp   ! degCd base 0 C, duration of the endosperm cell-division phase; FPB Table 1 p.784, from Gleadow 1982 + Singh & Jenner 1982 (p.786)

  ! ---- The lamina synthesis term of Eq. (6) ---- !
  ! gnf_sigma_lam = 1.5e-3 (degree-Cd)-1 is the relative rate of lamina
  ! photosynthetic-nitrogen SYNTHESIS printed by Bertheloot et al. (2008) Funct
  ! Plant Biol 35:781-796, Table 1 p.784. The same quantity is re-estimated ten
  ! times lower, 1.5e-4, by Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot
  ! 108:1097-1109, doi:10.1093/aob/mcr125, Table 2; both published estimates are
  ! coded, and the record selects between them.
  REAL(wp), PARAMETER :: gnf_sigma_lam   = 1.5e-3_wp    ! (degree-Cd)-1 relative rate of lamina photosynthetic-N synthesis; Bertheloot, Martre & Andrieu 2008 Funct Plant Biol 35:781-796, Table 1 p.784 ("sigma ... 0.0015"). Used by the record-86 path; at record 154 >= 1 the Bertheloot 2011b Table 2 estimate n2_sigma_lam = 1.5e-4 is used instead
  REAL(wp), PARAMETER :: gnf_r_mob_ph    = 0.10_wp    ! dimensionless ratio of mobile to photosynthetic N, so N^ph = (N^tot - N^struct)/(1+r); Bertheloot 2008 FPB Eqn 1 p.783 + p.786 "with r set to 0.10 ... reached after the first time steps ... irrespective of the initial value"
  ! REPORTED envelopes (falsification targets, never tuning targets)
  REAL(wp), PARAMETER :: gnf_mgN_kernel_hi = 1.38_wp  ! mg N kernel-1, same table, widest value incl. ablation treatments
  ! ---- runtime state (all reset per run; no deck reads any of them ----
  LOGICAL  :: gnf_snap        = .FALSE.   ! first post-anthesis call taken?
  REAL(wp) :: gnf_ng0         = 0.0_wp    ! g N m-2, = gnf_ng0_per_ear * gnx_ears_m2
  REAL(wp) :: gnf_tt          = 0.0_wp    ! degCd base 0 C accumulated after anthesis
  REAL(wp) :: gnf_pot_cum     = 0.0_wp    ! g N m-2, ^pot N_g(t) -- FPB Eqn 12 integrated
  REAL(wp) :: gnf_sink_day    = 0.0_wp    ! g N m-2 d-1, ^pot S_g * dTT
  REAL(wp) :: gnf_supply_day  = 0.0_wp    ! g N m-2 d-1, assim + mobile pool + delta-limited organ potential
  REAL(wp) :: gnf_draw_lam    = 0.0_wp    ! g N m-2 d-1 taken from the lamina today
  REAL(wp) :: gnf_draw_stem   = 0.0_wp    ! g N m-2 d-1 taken from the stem today
  REAL(wp) :: gnf_draw_chaff  = 0.0_wp    ! g N m-2 d-1 taken from the chaff today
  REAL(wp) :: gnf_cum_lam     = 0.0_wp    ! g N m-2 season total from the lamina
  REAL(wp) :: gnf_cum_stem    = 0.0_wp    ! g N m-2 season total from the stem
  REAL(wp) :: gnf_cum_chaff   = 0.0_wp    ! g N m-2 season total from the chaff
  REAL(wp) :: gnf_cum_assim   = 0.0_wp    ! g N m-2 season total straight from assimilation
  REAL(wp) :: gnf_cum_mob     = 0.0_wp    ! g N m-2 season total from the common mobile pool
  INTEGER  :: gnf_days_fill   = 0         ! days with a live sink (post-anthesis, pre-termination)
  INTEGER  :: gnf_days_bind   = 0         ! days on which the SINK bound (sink <= supply)
  INTEGER  :: gnf_first_bind_dpa = -999   ! d after anthesis of the first sink-bound day; -999 = never
  INTEGER  :: gnf_last_bind_dpa  = -999   ! d after anthesis of the last sink-bound day; -999 = never
  INTEGER  :: gnf_bind_after_15d = 0      ! sink-bound days beyond 15 d post-anthesis -- Martre 2003 p.1963 says this should be RARE
  REAL(wp) :: gnf_nres_max    = 0.0_wp    ! g N m-2 d-1, max |organ-state N residual| -- GATE D, NOT a tautology
  REAL(wp) :: gnf_nres_day    = 0.0_wp    ! g N m-2 d-1, today's organ-state N residual
  ! ---- Eqn-6 repair state / diagnostics (-------------- !
  LOGICAL  :: gnf_eqn6        = .FALSE.   ! .TRUE. under rec 86 = 1, 2 or 4; .FALSE. under 5
  REAL(wp) :: gnf_deg_lam_day = 0.0_wp    ! g N m-2 d-1, FPB Eqn 5 lamina degradation credited to the pool today
  REAL(wp) :: gnf_syn_lam_day = 0.0_wp    ! g N m-2 d-1, FPB Eqn 6 lamina synthesis actually drawn from the pool today
  REAL(wp) :: gnf_syn_pot_day = 0.0_wp    ! g N m-2 d-1, FPB Eqn 6 POTENTIAL before the shared-pool scaling
  REAL(wp) :: gnf_fN_day      = 0.0_wp    ! dimensionless, [N]mob/([N]mob + k1) (FPB Eqn 6 / Eqn 7)
  REAL(wp) :: gnf_fPAR_day    = 0.0_wp    ! dimensionless, PAR/(PAR + k2) (FPB Eqn 6)
  REAL(wp) :: gnf_pool_scale  = 1.0_wp    ! dimensionless, the SHARED scaling applied to both claims (1 = pool met both)
  REAL(wp) :: gnf_cum_syn_lam = 0.0_wp    ! g N m-2 season total returned to the lamina by Eqn 6
  REAL(wp) :: gnf_cum_deg_lam = 0.0_wp    ! g N m-2 season total degraded out of the lamina by Eqn 5

  ! ======================================================================== !
  ! Basis of the grain-nitrogen initial condition - driver record 115, accepted
  ! at the value 0, so the per-ear basis runs. The per-kernel alternative is not
  ! reachable in this version. The weighed grain nitrogen at flowering behind
  ! both bases is Bertheloot, Martre & Andrieu (2008) Funct Plant Biol
  ! 35:781-796, doi:10.1071/FP08064.
  INTEGER  :: gnf_ng0_basis   = 0         ! runtime flag (driver_main.txt rec 115; EOF-default 0). 0 = PER EAR (incumbent) / 1 = PER KERNEL
  INTEGER  :: gnk_tier        = -1        ! which tier fixed the divisor: -1 = never evaluated, 0 = main stem, 1 = FALLBACK A, 2 = FALLBACK B
  INTEGER  :: gnk_cohort      = -1        ! c*, the divisor cohort index (-1 = none found)
  INTEGER  :: gnk_fb_events   = 0         ! FALLBACK counter (tiers 1 and 2). gnf_snap is one-off, so this is 0 or 1 BY CONSTRUCTION -- it is still counted, so a future multi-snap law cannot make it silent
  INTEGER  :: gnk_fb_day      = -1        ! model day the fallback fired (-1 = never); the deck is the run.log's own folder
  REAL(wp) :: gnk_divisor     = 0.0_wp    ! kernels per spike of cohort c*, the per-kernel divisor actually used
  REAL(wp) :: gnk_numer       = 0.0_wp    ! kernels m-2 at the snap = kernel_number_m2
  REAL(wp) :: gnk_ears_snap   = 0.0_wp    ! ears m-2 at the snap = gnx_ears_m2 (the incumbent basis, for the reported factor)
  REAL(wp) :: gnk_ng0_perear  = 0.0_wp    ! g N m-2 the INCUMBENT would have set on the same state -- REPORTED, never used
  REAL(wp) :: gnk_factor      = 0.0_wp    ! dimensionless (kernel_number_m2 / gnx_kn(c*)) / gnx_ears_m2 -- the whole effect of rec 115, REPORTED

  ! ======================================================================== !
  ! Grain number - driver record 89, accepted at the value 1.
  ! Grain number follows the spike dry mass at anthesis through a measured
  ! fruiting efficiency. Source: Fischer, R. A. (1985) "Number of kernels in
  ! wheat crops and the influence of solar radiation and temperature",
  ! J Agric Sci Camb 105:447-461, doi:10.1017/S0021859600056495:
  ! - gnx_spike_part = 0.44, the spike share of the above-ground dry-matter
  !   increment over the spike-growth period (p.452);
  ! - the spike-growth period itself, about 300 degree-Cd ending at anthesis
  !   (p.449);
  ! - the per-floret spike dry weight 0.010 g (p.449) and the spike dry weight
  !   at anthesis 0.323 g (Table 4 p.451, cv. Egret);
  ! - the measured ear densities and kernel numbers used as falsification
  !   envelopes and applied to nothing: 436 spikes m-2 (p.450, cv. Yecora 70),
  !   the Table 4 p.451 anthesis and maturity densities, and the Table 5 p.454
  !   treatment series.
  ! Fruiting efficiency as a cultivar trait is Slafer, Elia, Savin, Garcia,
  ! Terrile, Ferrante, Miralles & Gonzalez (2015) Food Energy Secur
  ! 4(2):92-109, doi:10.1002/fes3.59.
  INTEGER  :: grain_number_mech = 0    ! runtime flag (driver_main.txt rec 89; EOF-default 0)

  ! --- Fischer (1985) spike constants --- !
  REAL(wp), PARAMETER :: gnx_spike_part   = 0.44_wp    ! -, spike share of the ABOVE-GROUND DM increment over the spike growth period; p.452
  REAL(wp), PARAMETER :: gnx_m_floret     = 0.010_wp   ! g spike DM per competent floret at anthesis ("1 floret per 10 mg"); p.449, corroborated Table 4 p.451 (102.3/103.3 florets g-1)
  REAL(wp), PARAMETER :: gnx_kern_per_fl  = 0.78_wp    ! kernels per competent floret = (78 +- 2 kernels g-1 spike DW, p.447) / (100 florets g-1, p.449); envelope 0.78-0.92 (Table 4 p.451) REPORTED
  REAL(wp), PARAMETER :: gnx_tt_base_C    = 4.5_wp          ! degC, the base of the degree-day statement of Fischer 1985 p.449 ("300 day degrees above 4-5 degC"). This is the window-opening rule applied in growth_module.f90 update_spike_window
  REAL(wp), PARAMETER :: gnx_tt_ref_Cd    = 300.0_wp      ! degCd above 4.5 degC, the spike-growth-period duration of Fischer 1985 p.449, cv. Yecora 70. This is the threshold tested in update_spike_window; it is applied, not merely reported

  ! --- Fischer (1985) measured envelopes; applied to nothing --- !
  REAL(wp), PARAMETER :: gnx_env_spikes_t3 = 436.0_wp  ! spikes m-2, p.450 PROSE, unshaded control (Yecora 70); Table 3 is p.451 and prints slopes only
  REAL(wp), PARAMETER :: gnx_env_spklt_t3  = 19.4_wp   ! spikelets per spike, p.450 PROSE, unshaded control (Table 3 = p.451, slopes only)
  REAL(wp), PARAMETER :: gnx_env_kspk_t3   = 1.98_wp   ! kernels per spikelet, p.450 PROSE, unshaded control (Table 3 = p.451, slopes only)
  REAL(wp), PARAMETER :: gnx_env_spikes_t4 = 598.0_wp  ! spikes m-2 at anthesis, Table 4 p.451 (Egret)
  REAL(wp), PARAMETER :: gnx_env_flspk_t4  = 32.3_wp   ! florets per spike at anthesis, Table 4 p.451
  REAL(wp), PARAMETER :: gnx_env_kspk_t4   = 29.6_wp   ! kernels per spike at maturity, Table 4 p.451

  ! --- Derived spike quantities. No new constant: each is a product of the
  ! measured shares above with the measured spike dry weight --- !
  REAL(wp), PARAMETER :: bgn_spike_dw_ref    = gnx_env_flspk_t4 * gnx_m_floret        ! 0.323 g spike-1 at anthesis = 32.3 florets (Fischer Tab 4 p.451) x 10 mg per floret (p.449)
  REAL(wp), PARAMETER :: bgn_r_spike         = bgn_spike_dw_ref / gnx_tt_ref_Cd       ! 1.0766667e-3 g axis-1 (degCd above 4.5)-1 = 0.323 g / 300 degCd (Fischer p.449)
  REAL(wp), PARAMETER :: bgn_p_rachis        = 0.171_wp    ! Rivera-Amado 2019 p.160 Tab 7: rachis share of non-grain spike DM
  REAL(wp), PARAMETER :: bgn_p_awns          = 0.260_wp    ! same, awns (the highest single component)
  REAL(wp), PARAMETER :: bgn_p_glumes        = 0.204_wp    ! same, glumes
  REAL(wp), PARAMETER :: bgn_p_lemmas        = 0.237_wp    ! same, lemmas (FE r = +0.64 / +0.61)
  REAL(wp), PARAMETER :: bgn_p_paleas        = 0.127_wp    ! same, paleas (FE r = +0.43 trend)
  REAL(wp), PARAMETER :: bgn_f_struct        = bgn_p_rachis + bgn_p_awns + bgn_p_glumes   ! 0.635 STRUCTURAL

  ! The five printed component shares of Rivera-Amado et al. (2019) sum to
  ! 0.999 rather than 1.000; the printed values are kept and the residual is
  ! carried explicitly rather than renormalized silently.
  REAL(wp), PARAMETER :: bgn_f_floral_print  = bgn_p_lemmas + bgn_p_paleas               ! 0.364 as PRINTED by Rivera-Amado
  REAL(wp), PARAMETER :: bgn_f_floral        = 1.0_wp - bgn_f_struct                     ! 0.365 FLORAL, exhaustive with STRUCTURAL

  ! --- state: the per-axis spike pool and its endpoints ---------------- !
  LOGICAL  :: gnx_win_open        = .FALSE.   ! crop-level: the main stem's penultimate leaf has appeared
  INTEGER  :: gnx_win_open_day    = -1        ! model day it opened
  INTEGER  :: gnx_win_days        = 0         ! days the window has been open (to anthesis)
  REAL(wp) :: gnx_win_tt45        = 0.0_wp    ! degCd above 4.5 degC accumulated in the window — REPORTED against gnx_tt_ref_Cd
  REAL(wp) :: gnx_win_agb0        = 0.0_wp    ! above-ground DM at window open, g m-2
  REAL(wp) :: gnx_win_spike0      = 0.0_wp    ! spike (storage) DM at window open, g m-2
  REAL(wp) :: gnx_win_dagb        = 0.0_wp    ! above-ground DM increment over the window, g m-2
  REAL(wp) :: gnx_win_dspike      = 0.0_wp    ! spike DM increment over the window, g m-2
  REAL(wp) :: gnx_share_real      = 0.0_wp    ! realized window share = gnx_win_dspike / gnx_win_dagb — the item-D gate
  INTEGER  :: gnx_clamp_days      = 0         ! days the spike share had to be clamped because the lamina sink left no room
  REAL(wp) :: gnx_clamp_max       = 0.0_wp    ! largest clamp, in share units
  REAL(wp) :: gnx_pool(0:CT_MAXTIL)   = 0.0_wp ! g m-2 of spike DM accrued by each axis class
  REAL(wp) :: gnx_dens(0:CT_MAXTIL)   = 0.0_wp ! axes m-2 of that class at anthesis (snapshot)
  REAL(wp) :: gnx_wsp(0:CT_MAXTIL)    = 0.0_wp ! g per axis, spike DM at anthesis
  REAL(wp) :: gnx_nfl(0:CT_MAXTIL)    = 0.0_wp ! competent florets per spike
  REAL(wp) :: gnx_kn(0:CT_MAXTIL)     = 0.0_wp ! kernels per spike
  LOGICAL  :: gnx_is_ear(0:CT_MAXTIL) = .FALSE.! that class bears a grain-bearing ear (>= 1 kernel per spike)
  REAL(wp) :: gnx_unshared        = 0.0_wp    ! g m-2 of spike DM that fell on a day with no open axis (all credited to the main stem) — REPORTED
  REAL(wp) :: gnx_ears_m2         = 0.0_wp    ! ears m-2 — PREDICTED
  REAL(wp) :: gnx_florets_m2      = 0.0_wp    ! competent florets m-2 — PREDICTED
  REAL(wp) :: gnx_grains_per_ear  = 0.0_wp    ! kernels per ear — PREDICTED
  REAL(wp) :: gnx_spike_dw_anth   = 0.0_wp    ! g m-2, spike DW at anthesis (= W_chaff when the number is set)
  REAL(wp) :: gnx_sterile_C       = 0.0_wp    ! g m-2 of spike DM held by classes that bear no ear (stays in W_chaff)
  LOGICAL  :: gnx_set             = .FALSE.   ! the chain has fired
  REAL(wp) :: gnx_pool_scale      = 0.0_wp    ! rescale applied so SUM_c gnx_pool(c) = W_chaff exactly — REPORTED
  REAL(wp) :: gnx_dagb_day        = 0.0_wp    ! g m-2 d-1, today's ABOVE-GROUND DM increment = dm_gain*(1-f_root) — Fischer's basis
  REAL(wp) :: gnx_dspike_day      = 0.0_wp    ! g m-2 d-1, today's spike DM increment = dm_gain*f_storage
  REAL(wp) :: gnx_sum_tdu         = 0.0_wp    ! season-to-date SUM of the phenology module's own daily thermal-development units
  REAL(wp) :: gnx_sum_cd45        = 0.0_wp    ! season-to-date SUM of max(0, Tair_day - 4.5), degCd -- Fischer's own base, p.449
  REAL(wp) :: gnx_tt_ratio        = 0.0_wp    ! - , gnx_sum_cd45 / gnx_sum_tdu; the assumed unit conversion used to project the window; REPORTED
  REAL(wp) :: gnx_tdu_rem_open    = 0.0_wp    ! tdu still to anthesis on the day the window opened -- REPORTED
  INTEGER  :: gnx_day_now         = 0         ! model day, published by update_spike_window for the reported cross-check
  INTEGER  :: gnx_penult_day      = -1        ! model day the MAIN STEM reached Haun FLN-1 (Fischer's developmental form) -- REPORTED CROSS-CHECK, read by nothing

  ! ======================================================================== !
  ! Floret-site chain - driver record 90, accepted at the value 4.
  ! Floret primordia per spike, their survival to a fertile floret and the
  ! resulting grain number are carried per axis. Measured bases:
  ! - maximum floret site number per spikelet, about 9-11 in the field:
  !   Guo, Slafer & Schnurbusch (2016) J Exp Bot 67:4241,
  !   doi:10.1093/jxb/erw200 (MFS varies by about 20 %, from 8-10 primordia in
  !   the greenhouse and 9-11 in the field), corroborated by Backhaus,
  !   Boeven, Schnurbusch et al. (2023) J Exp Bot 74:5088,
  !   doi:10.1093/jxb/erad233 (interquartile 9-10, maximum 11);
  ! - the detillering evidence that floret survival, not primordium number, is
  !   the variable component: Guo & Schnurbusch (2015) J Exp Bot,
  !   doi:10.1093/jxb/erv303, section "Influence of detillering on maximum
  !   floret number";
  ! - the elevated-CO2 direction: Prieto, Ochagavia, Savin, Griffiths & Slafer
  !   (2020) J Exp Bot 71:2703, doi:10.1093/jxb/erz568, Discussion.
  INTEGER  :: floret_survival_mech = 0          ! runtime flag (driver_main.txt rec 90; EOF-default 0)
  REAL(wp), PARAMETER :: gnx_spikelets = 19.4_wp  ! spikelets per spike; Fischer 1985 p.450, cv. Yecora 70
  REAL(wp), PARAMETER :: gnx_mfs       = 10.0_wp  ! maximum floret primordia per spikelet; Guo & Schnurbusch 2015 JXB 66:5945 (envelope 9-12)

  ! --- Spike composition masses, derived from the measured shares above --- !
  REAL(wp), PARAMETER :: bgn_m_struct_spklt  = bgn_f_struct * bgn_spike_dw_ref / gnx_spikelets     ! 1.0572417e-2 g spikelet-1 (rachis+awns+glumes)
  REAL(wp), PARAMETER :: bgn_m_floret_floral = bgn_f_floral * bgn_spike_dw_ref / gnx_env_flspk_t4  ! 3.6400000e-3 g floret-1 (lemma+palea)
  REAL(wp), PARAMETER :: bgn_w_struct_axis   = gnx_spikelets * bgn_m_struct_spklt                  ! 0.20510492 g axis-1, the spike's own structural claim

  ! Maximum floret sites per spikelet; compared against the per-spikelet
  ! measurements of Guo et al. (2016) and Backhaus et al. (2023), not against a
  ! crop mean.
  REAL(wp) :: gnx_ninit          = 0.0_wp        ! floret primordia per spike = gnx_spikelets * gnx_mfs
  REAL(wp) :: gnx_surv(0:CT_MAXTIL)  = 0.0_wp    ! floret survival fraction of that class — PREDICTED
  REAL(wp) :: gnx_pool_pre(0:CT_MAXTIL) = 0.0_wp ! g m-2, the class's spike pool BEFORE any mortality redistribution — REPORTED
  INTEGER  :: gnx_axes_culled    = 0             ! classes culled by the fertile-tiller test
  REAL(wp) :: gnx_pool_moved     = 0.0_wp        ! g m-2 of spike DM returned from culled classes to survivors — REPORTED, carbon-neutral
  INTEGER  :: gnx_cull_passes    = 0             ! passes the monotone cull loop took
  REAL(wp) :: gnx_fe_realised    = 0.0_wp        ! kernels per g spike DW at anthesis — THE HEADLINE TEST; must vary across decks
  ! rec 90 = 4 state (no PARAMETER; set only inside the rec 90 = 4 branch, printed only when rec 90 = 4)
  REAL(wp) :: gnx_phi_sp(0:CT_MAXTIL) = 1.0_wp   ! ear-bearing FRACTION of the class = MIN(1, W_axis / W_pot) — PREDICTED; main stem 1
  REAL(wp) :: gnx_wpot_axis      = 0.0_wp        ! g per axis, W_pot = bgn_r_spike x gnx_win_tt45 (Fischer's unshaded spike over the realized window)
  REAL(wp) :: gnx_earless_m2     = 0.0_wp        ! axes m-2 standing WITHOUT a spike at anthesis = SUM_c ct_dens(c) (1 - phi_sp) — the Fischer p.452 class, REPORTED

  ! ======================================================================== !
  ! Grain filling - driver record 91, accepted at the value 1.
  ! The kernel growth rate follows the determinate sigmoid of Yin, Goudriaan,
  ! Lantinga, Vos & Spiertz (2003) "A flexible sigmoid function of determinate
  ! growth", Annals of Botany 91(3):361-371, doi:10.1093/aob/mcg029, Eq. (9)
  ! with the rate form of its Eq. (A1), integrated on thermal time above a base
  ! temperature.
  ! The base temperature and the filling duration are derived from the
  ! reciprocal-duration method of Fischer, R. A. (1985) J Agric Sci Camb
  ! 105:447-461, Table 6 footnote, applied to the mean filling duration of Yin
  ! et al. (2003); base 0 C for wheat is the value Li, Cao, Dai & Jing (2000)
  ! Crop Sci 40:1263, p.1264, attribute to Bauer et al. (1985) and that Martre
  ! et al. (2003) Plant Physiol 133:1959 also use.
  ! Reported envelopes, falsification targets and never tuning targets: the
  ! kernel-mass and duration spreads of Yin et al. (2003) Table 2 and of Yang &
  ! Zhang (2006) Table 3 p.228, and the single-kernel mass of Wang et al. (2017)
  ! Nature Plants 3:17102, Fig. 4b.
  INTEGER  :: grain_fill_mech = 0                 ! runtime flag (driver_main.txt rec 91; EOF-default 0)

  ! ======================================================================== !
  ! Cultivar final-kernel-mass ceiling - driver record 118. Mode 1 caps the
  ! integrated kernel mass at the cultivar potential, because the Yin et al.
  ! (2003) Ann Bot 91:361 rate form, doi:10.1093/aob/mcg029, integrates to a
  ! mass that is not itself bounded; the measured ceiling is a genotype
  ! property (Cook & Wardlaw 1977, Aust J Plant Physiol 4:785,
  ! doi:10.1071/PP9770785).
  INTEGER  :: kernel_mass_ceiling = 0             ! rec 118; 0 = OFF (EOF-default); 1 = cultivar w_max ceiling
  REAL(wp), PARAMETER :: gfb_cm_ref   = 1.928333_wp    ! mg kernel-1 d-1, mean c_m of Yin 2003 Table 2 regime T1
  REAL(wp), PARAMETER :: gfb_cm_slope = 0.1106667_wp   ! mg kernel-1 d-1 K-1 = (2.481667 - 1.928333)/5
  REAL(wp), PARAMETER :: gfb_t_ref    = 17.5_wp        ! degC, daily mean of Yin's 20/15 regime at a 12 h photoperiod — model assumption

  ! gfb_tbase and gfb_dfill are cultivar-file traits, read from
  ! input/cultivar_yield.csv, with the compiled values as defaults. The
  ! defaults are derived by the Fischer (1985) Table 6 reciprocal-duration
  ! method from four measured grain-filling durations (Yin et al. 2003;
  ! Baillot et al. 2018; Girousse et al. 2005; Bertheloot, Martre & Andrieu
  ! 2008 Plant Physiol 148:1707).
  REAL(wp) :: gfb_tbase    = 6.618733_wp    ! degC, DERIVED by Fischer 1985 Tab 6 reciprocal-duration method from Yin's mean t_e; DEFAULT, overridable by cultivar_yield.csv row 9 gfb_tbase_C
  REAL(wp) :: gfb_dfill    = 405.58109_wp   ! degCd above gfb_tbase, same derivation, identical from both regimes; DEFAULT, overridable by cultivar_yield.csv row 10 gfb_dfill_Cd
  ! --- provenance of the two values actually in force on this deck (reported in run.log ---
  LOGICAL  :: gfb_tbase_from_file = .FALSE. ! .TRUE. only when cultivar_yield.csv row 9 supplied it
  LOGICAL  :: gfb_dfill_from_file = .FALSE. ! .TRUE. only when cultivar_yield.csv row 10 supplied it
  REAL(wp), PARAMETER :: gfb_u        = 0.433143_wp    ! t_m/t_e, mean of ALL TWELVE rows of Yin 2003 Table 2

  ! Analytic and quadrature evaluations of the shape integral agree to the
  ! printed precision, so the compiled constant is the analytic value.
  REAL(wp), PARAMETER :: gfb_shape_int= 0.685644_wp    ! integral of g(x) over [0,1] at gfb_u; = mean rate / c_m. Analytic = quadrature = 0.685644546
  ! REPORTED envelopes (falsification, never tuning). Yin 2003 Table 2 spreads and Yang & Zhang 2006 Table 3 p.228.
  ! ---- state, all REPORTED in output/grain_filling_law.csv ----
  REAL(wp) :: gfb_tt            = 0.0_wp   ! degCd above gfb_tbase accumulated since the grain number was set
  INTEGER  :: gfb_day0          = 0        ! model day the grain number was set (= filling start)
  INTEGER  :: gfb_day_end       = 0        ! model day filling ended
  INTEGER  :: gfb_days          = 0        ! days on which the law was live
  INTEGER  :: gfb_end_bound     = 0        ! 0 = still filling, 1 = thermal duration, 2 = dvs2 >= 2 (phenology maturity), 3 = cultivar kernel-mass ceiling (rec 118)
  REAL(wp) :: gfb_tsum          = 0.0_wp   ! SUM Tair_day over the filling window (for the reported mean)
  REAL(wp) :: gfb_demand_tot    = 0.0_wp   ! g m-2, cumulative potential sink demand
  REAL(wp) :: gfb_supplied_tot  = 0.0_wp   ! g m-2, cumulative carbon actually delivered to grain (assimilate + reserve)
  REAL(wp) :: gfb_unmet_tot     = 0.0_wp   ! g m-2, cumulative UNMET demand = the source-limitation measure
  REAL(wp) :: gfb_reserve_end   = 0.0_wp   ! g m-2, W_reserve after supply on the LAST day the law was live
  REAL(wp) :: gfb_rate_max      = 0.0_wp   ! mg kernel-1 d-1, largest realized daily rate
  REAL(wp) :: gfb_cm_min        = 1.0e30_wp! mg kernel-1 d-1, smallest c_m(T) over the window
  REAL(wp) :: gfb_cm_max        = 0.0_wp   ! mg kernel-1 d-1, largest c_m(T) over the window
  LOGICAL  :: gfb_cm_clamped    = .FALSE.  ! .TRUE. if c_m(T) ever hit the zero clamp (would mean T below 0.08 degC)

  ! ======================================================================== !
  ! Pre-anthesis stem water-soluble-carbohydrate reserve - driver record 92,
  ! accepted at the value 1.
  ! The reserve is a bounded store deposited before anthesis and drawn during
  ! grain filling. Source: Schnyder, H. (1993) "The role of carbohydrate
  ! storage and redistribution in the source-sink relations of wheat and barley
  ! during grain filling - a review", New Phytol 123:233-245,
  ! doi:10.1111/j.1469-8137.1993.tb03731.x, p.238, which gives the measured
  ! ceiling as a fraction of internode dry weight in wheat and barley (Austin
  ! et al. 1977a,b, 1980b) - the quantity wsc_cap_frac holds. The genotypic and
  ! environmental span of the store itself is Ruuska et al. (2006) Funct Plant
  ! Biol 33:799, measured on stem samples (p.803).
  INTEGER  :: stem_wsc_mech = 0                       ! runtime flag (driver_main.txt rec 92; EOF-default 0)
  REAL(wp), PARAMETER :: wsc_cap_frac      = 0.50_wp  ! g WSC per g stem DM ceiling; Schnyder 1993 New Phytol 123:233 p.238(a)
  INTEGER,  PARAMETER :: wsc_accum_end_dpa = 15       ! d after anthesis at which NET WSC deposition stops; Schnyder 1993 p.238 "10-20 d" midpoint(b)
  REAL(wp), PARAMETER :: wsc_elong_dvs2    = 0.646_wp ! dvs2 of the WOFOST elongation node, USED ONLY BY THE rec 92 = 4 SENSITIVITY ARM(ii)

  ! Reported envelopes: falsification targets, never tuning targets.
  INTEGER  :: wsc_day_anth    = 0        ! model day the grain number was set (= anthesis on every shipped deck)
  REAL(wp) :: wsc_pool_anth   = -1.0_wp  ! g DM m-2, W_reserve on the anthesis day (-1 = not yet reached)
  REAL(wp) :: wsc_stem_anth   = -1.0_wp  ! g DM m-2, W_stem on the anthesis day
  REAL(wp) :: wsc_pool_peak   = 0.0_wp   ! g DM m-2, largest W_reserve of the run
  INTEGER  :: wsc_day_peak    = 0        ! model day of that peak
  REAL(wp) :: wsc_dep_pre     = 0.0_wp   ! g DM m-2, total LABELED pre-anthesis deposition
  REAL(wp) :: wsc_dep_post    = 0.0_wp   ! g DM m-2, total LABELED post-anthesis (sink-full surplus) deposition
  REAL(wp) :: wsc_refused_cap = 0.0_wp   ! g DM m-2, deposition refused by the CAPACITY (stays in W_stem, unlabeled)
  REAL(wp) :: wsc_refused_win = 0.0_wp   ! g DM m-2, deposition refused by the DEPOSITION CLOSE (stays in W_stem, unlabeled)
  REAL(wp) :: wsc_draw_tot    = 0.0_wp   ! g DM m-2, total drawn from the pool
  REAL(wp) :: wsc_draw_after16= 0.0_wp      ! g DM m-2, of that, drawn more than 16 d after anthesis. No published source was located for a 16-day threshold: it is in neither Gebbing & Schnyder 1999 Plant Physiol 121:871-878 nor Gebbing, Schnyder & Kuhbauch 1998 Plant Cell Environ 21:301-313, and the 1998 study sampled only at anthesis and at 35 d after anthesis (p.303). Diagnostic column only: this accumulator is read by two WRITE statements and by nothing else - no flux, state, ledger or gate reads it
  REAL(wp) :: wsc_gain_tot    = 0.0_wp   ! g DM m-2, total delivered to grain (= draw * reserve_conv_eff)
  REAL(wp) :: wsc_resp_tot    = 0.0_wp   ! g DM m-2, total conversion loss
  INTEGER  :: wsc_day_mob0    = 0        ! model day of the FIRST draw (0 = never drew)
  INTEGER  :: wsc_days_capbnd = 0        ! days on which the capacity bound a deposition
  INTEGER  :: wsc_days_winbnd = 0        ! days on which the deposition close bound a deposition
  REAL(wp) :: wsc_closure_max = 0.0_wp   ! g DM m-2 d-1, largest |organ-state closure residual|
  REAL(wp) :: wsc_dep_day     = 0.0_wp   ! g DM m-2 d-1, today's labeled deposition (daily file)

  ! ======================================================================== !
  ! Post-anthesis reserve law - driver records 105 and 106. Record 105 is
  ! accepted at the value 2 and record 106 at the value 0, so the bounded,
  ! deficit-independent release runs on the grain-filling clock.
  INTEGER  :: wsc_reserve_law = 0        ! rec 105 (EOF-default 0); 0 = OFF, 1 = E4 bounded draw, 2 = the REPAIR, deficit-independent release
  INTEGER  :: wsc_law_member  = 0                ! rec 106 (EOF-default 0 = not selected); 1 = the grain-filling clock, 2 = the record-59 clock, 3 = no ceiling
  REAL(wp) :: wsl_ddvs2       = 0.0_wp   ! today's dvs2 increment (main: dvs2_rate_growth), handed in ON-only for member 2
  ! ---- today's E4 terms (output/daily_wsc_law.csv ----
  REAL(wp) :: wsl_dtau_day    = 0.0_wp   ! dtau of the selected clock (degCd, member 1; dvs2 d-1, member 2)
  REAL(wp) :: wsl_tau_rem_day = 0.0_wp   ! tau_rem of the selected clock, remaining AFTER today
  REAL(wp) :: wsl_room_day    = 0.0_wp   ! g DM m-2 d-1, the unmet demand the draw was measured against
  REAL(wp) :: wsl_draw_def_day= 0.0_wp   ! g DM m-2 d-1, the rec 92 deficit draw BEFORE the ceiling
  REAL(wp) :: wsl_cap_day     = 0.0_wp   ! g DM m-2 d-1, cap_draw
  REAL(wp) :: wsl_draw_day    = 0.0_wp   ! g DM m-2 d-1, the realized draw
  INTEGER  :: wsl_bound_day   = 0        ! 0 none / 1 deficit / 2 pool / 3 stem / 4 ceiling (ceiling ONLY when cap_draw < draw_def)
  ! ---- run totals (run.log summary line ----
  INTEGER  :: wsl_days_ceil   = 0        ! days the ceiling strictly reduced the draw
  INTEGER  :: wsl_days_def    = 0        ! days the deficit was the binding term
  INTEGER  :: wsl_days_pool   = 0        ! days the pool was the binding term (pool emptied)
  INTEGER  :: wsl_days_stem   = 0        ! days W_stem was the binding term
  REAL(wp) :: wsl_frac_max    = 0.0_wp   ! largest one-day draw / W_reserve(before the draw)
  REAL(wp) :: wsl_wres_prev   = 0.0_wp   ! W_reserve at the previous daily row (ledger)
  REAL(wp) :: wsl_ledger_max  = 0.0_wp   ! max |W_reserve - (previous + dep - draw)| over the run

  ! ======================================================================== !
  ! Post-anthesis stem-reserve release - driver record 105 mode 2 with the
  ! window member of driver record 107, accepted at the value 2.
  ! The release is bounded by the pool and is independent of the grain deficit,
  ! which is what the measurements show: Schnyder (1993) New Phytol 123:233-245,
  ! doi:10.1111/j.1469-8137.1993.tb03731.x, p.239 reports that the store is
  ! mobilized when the ears are removed as well, so the draw is not driven by
  ! the grain demand, and Gebbing & Schnyder (1999) Plant Physiol 121:871-878,
  ! doi:10.1104/pp.121.3.871, p.877 gives the fraction of the pre-anthesis store
  ! that reaches the grain.
  INTEGER  :: wsc_rel_window_member = 0      ! rec 107 (EOF-default 0 = not selected); 1 = the grain-filling clock, 2 = the record-59 clock
  ! ---- today's mode-2 terms (output/daily_wsc_mode2.csv ----
  REAL(wp) :: wsl2_tau_day    = 0.0_wp   ! position on the SELECTED clock (gfb_tt member 1; dvs2 member 2)
  REAL(wp) :: wsl2_taurem_day = 0.0_wp   ! tau_rem of the selected clock, remaining AFTER today
  REAL(wp) :: wsl2_wres_pre   = 0.0_wp   ! g DM m-2, W_reserve immediately BEFORE the release
  REAL(wp) :: wsl2_rel_day    = 0.0_wp   ! g DM m-2 d-1, rel -- what LEFT the pool today (E1)
  REAL(wp) :: wsl2_relg_day   = 0.0_wp   ! g DM m-2 d-1, rel_g -- CASE 1 stream (E2)
  REAL(wp) :: wsl2_relf_day   = 0.0_wp   ! g DM m-2 d-1, rel_f -- CASE 2 stream (E2), NO efficiency
  REAL(wp) :: wsl2_gain_day   = 0.0_wp   ! g DM m-2 d-1, gain = rel_g * reserve_conv_eff
  REAL(wp) :: wsl2_lconv_day  = 0.0_wp   ! g DM m-2 d-1, loss_conv = rel_g * (1 - eff)
  REAL(wp) :: wsl2_cres_day   = 0.0_wp   ! g DM m-2 d-1, cres_respired of the day (the NAMED ledger slot)
  INTEGER  :: wsl2_bound_day  = 0        ! 0 none / 2 pool / 3 stem / 4 window ceiling (no deficit term on this path)
  ! ---- run totals: the TWO streams stay SEPARATE ----
  REAL(wp) :: wsl2_rel_tot    = 0.0_wp   ! g DM m-2, SUM rel
  REAL(wp) :: wsl2_relg_tot   = 0.0_wp   ! g DM m-2, SUM rel_g (CASE 1)
  REAL(wp) :: wsl2_relf_tot   = 0.0_wp   ! g DM m-2, SUM rel_f (CASE 2)
  REAL(wp) :: wsl2_gain_tot   = 0.0_wp   ! g DM m-2, SUM gain (delivered to grain)
  REAL(wp) :: wsl2_loss_conv_tot = 0.0_wp ! g DM m-2, SUM loss_conv (CASE 1 conversion loss, Gebbing eff)
  REAL(wp) :: wsl2_loss_fate_tot = 0.0_wp ! g DM m-2, SUM loss_fate (CASE 2 named respiration fate, UNITY)
  INTEGER  :: wsl2_days_zerodef  = 0     ! days the pool released with the grain deficit already at zero
  ! ---- conservation residual maxima (run.log, double precision ----
  REAL(wp) :: wsl2_wres_prev  = 0.0_wp   ! W_reserve at the previous daily row (I4)
  REAL(wp) :: wsl2_i1_max     = 0.0_wp   ! max |rel - rel_g - rel_f|
  REAL(wp) :: wsl2_i2_max     = 0.0_wp   ! max |gain + loss_conv + loss_fate - rel|
  REAL(wp) :: wsl2_i3_max     = 0.0_wp   ! max |cres_respired - (loss_conv + loss_fate)|
  REAL(wp) :: wsl2_i4_max     = 0.0_wp   ! max |W_reserve - (previous + dep - rel)|

  ! ======================================================================== !
  ! Post-anthesis reserve basis of canopy height - driver record 108, accepted
  ! at the value 0, so height is not rebased on the reserve draw. The reserve
  ! concentration bracket the alternative would use is printed in Cereal
  ! Research Communications 40(3):413-422, doi:10.1556/CRC.40.2012.3.10.
  INTEGER  :: height_reserve_basis = 0   ! rec 108 (EOF-default 0 = OFF); 1 = x_eff = W_stem - (W_reserve - hgt_wres_anth)
  REAL(wp) :: hgt_wres_anth = 0.0_wp     ! g DM m-2, W_reserve latched at the anthesis date (dvs2 < reserve_start_dvs2)

  ! ======================================================================== !
  ! Canopy height from internode elongation - driver records 109 to 113.
  ! Record 109 selects the law (accepted values 0 and 1) and records 110 to 113
  ! select its published members: the elongation window (110, values 0 to 6),
  ! the final internode lengths (111, values 0 to 5), the lineic mass (112,
  ! values 0 to 4) and the ear or peduncle-end length (113, values 0 to 3).
  ! The released configuration sets all five records to 0, so canopy height
  ! keeps the incumbent law. Every non-zero member below is an inactive
  ! alternative, see Part 1 Table 6.
  ! Height is built per inversed internode rank from a length profile, a lineic
  ! mass and an elongation window, with the ear length added at the top.
  INTEGER, PARAMETER :: INL_MAXIR = 11   ! array bound only: inversed ranks 0 (peduncle) .. 11
  INTEGER  :: height_internode_law  = 0    ! rec 109 (EOF-default 0 = OFF); 1 = canopy height from internode elongation
  INTEGER  :: hgt_int_window_member = 0  ! rec 110 (EOF-default 0 = not selected); 1-6 = onset arm x Delta_I
  INTEGER  :: hgt_int_length_member = 0  ! rec 111 (EOF-default 0 = not selected); 1-5 = the length-potential bracket
  INTEGER  :: hgt_int_lineic_member = 0  ! rec 112 (EOF-default 0 = not selected); 1-4 = the lineic-mass bracket
  INTEGER  :: hgt_int_ear_member    = 0  ! rec 113 (EOF-default 0 = not selected); 1-3 = the ear / peduncle-end bracket

  ! Record 111 = 1: Kirby, E. J. M. (1988) Field Crops Research 18:127-140,
  ! Table 3 folio 133, the upper measured internode lengths (glasshouse).
  ! Inactive alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_Lf_m1(0:INL_MAXIR) = [ 34.70_wp, 18.10_wp, 11.70_wp,  8.69_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp ]

  ! Record 111 = 2: Rivera-Amado, Trujillo-Negrellos, Molero, Reynolds,
  ! Sylvester-Bradley & Foulkes (2019) Field Crops Research 240:154-167,
  ! doi:10.1016/j.fcr.2019.04.016, Table 6 folio 161, the table's own printed
  ! means. Inactive alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_Lf_m2(0:INL_MAXIR) = [ 36.70_wp, 20.50_wp, 14.10_wp, 14.80_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp ]

  ! Record 111 = 3: Rivera-Amado et al. (2019) Table 6 folio 161, the table's own
  ! printed minima. Inactive alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_Lf_m3(0:INL_MAXIR) = [ 22.20_wp, 15.80_wp,  9.50_wp,  8.10_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp ]

  ! Record 111 = 4: Rivera-Amado et al. (2019) Table 6 folio 161, the table's own
  ! printed maxima. Inactive alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_Lf_m4(0:INL_MAXIR) = [ 44.20_wp, 23.30_wp, 16.50_wp, 19.50_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp ]

  ! Record 111 = 5: the continuity arm - the same internode-length profile
  ! cse_li_irank the sheath and stem geometry already uses (Abichou 2016 Ch.2
  ! Fig. 3c). Inactive alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_Lf_m5(0:INL_MAXIR) = [  0.00_wp, 16.00_wp, 11.00_wp,  7.00_wp, &
                                                     4.00_wp,  1.50_wp,  0.00_wp,  0.00_wp, &
                                                     0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp ]

  ! Record 110, the elongation window: onset arm times duration. The printed
  ! durations are Kirby (1988) Field Crops Research 18:127-140, folio 135
  ! (glasshouse). Inactive alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_dI_cd(3) = [ 277.0_wp, 360.0_wp, 476.0_wp ]  ! degCd, Kirby 1988 folio 135 PRINTED (glasshouse; see E3)
  INTEGER,  PARAMETER :: INL_ONSET_SHIFT = 2   ! phytomer ranks, the Malvoisin i-2 alpha-trap (folio 589-590); a RANK OFFSET, not a fitted constant

  ! Record 112, the lineic mass lambda_x (mg cm-1 per rank). Members 1 to 3 are
  ! the mean, minimum and maximum rows printed by Rivera-Amado et al. (2019)
  ! Field Crops Research 240:154-167, Table 6 folio 161. Inactive alternative,
  ! see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_lam_m1(0:INL_MAXIR) = [ 11.0_wp, 15.4_wp, 17.6_wp, 22.3_wp, &
                                                     22.3_wp, 22.3_wp, 22.3_wp, 22.3_wp, &
                                                     22.3_wp, 22.3_wp, 22.3_wp, 22.3_wp ]

  ! Record 112 = 2: the same table's own printed minimum row. Inactive
  ! alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_lam_m2(0:INL_MAXIR) = [  9.0_wp, 12.3_wp, 13.3_wp, 16.7_wp, &
                                                     16.7_wp, 16.7_wp, 16.7_wp, 16.7_wp, &
                                                     16.7_wp, 16.7_wp, 16.7_wp, 16.7_wp ]

  ! Record 112 = 3: the same table's own printed maximum row. Inactive
  ! alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_lam_m3(0:INL_MAXIR) = [ 13.2_wp, 18.0_wp, 21.0_wp, 27.6_wp, &
                                                     27.6_wp, 27.6_wp, 27.6_wp, 27.6_wp, &
                                                     27.6_wp, 27.6_wp, 27.6_wp, 27.6_wp ]

  ! Record 112 = 4: Gauthier, Barillot, Andrieu et al. (2020) J Exp Bot,
  ! Supplementary Data S1, Table of internode lineic masses. Inactive
  ! alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_lam_m4(0:INL_MAXIR) = [  7.0_wp, 14.0_wp, 16.0_wp, 17.0_wp, &
                                                     23.0_wp, 28.0_wp, 28.0_wp, 28.0_wp, &
                                                     28.0_wp, 28.0_wp, 28.0_wp, 28.0_wp ]

  ! Record 113, the ear or peduncle-end length. Member 1 is Kirby (1988) Field
  ! Crops Research 18:127-140, Table 3 folio 133 (ear length 101 mm,
  ! glasshouse); the remaining members span the measured 8-12 cm bracket.
  ! Inactive alternative, see Part 1 Table 6.
  REAL(wp), PARAMETER :: inl_Lear_kirby_cm = 10.1_wp  ! cm, Kirby 1988 Table 3 folio 133 C_ear = 101 mm (glasshouse); bracket 8-12 cm lives in rec 113, not here

  ! --- Internode state. Identically zero or sentinel when record 109 = 0 --- !
  LOGICAL  :: inl_ready = .FALSE.                     ! the previous day's structural stem mass has been latched
  LOGICAL  :: inl_ped_open = .FALSE.                  ! the peduncle window has opened at least once -- a monotone EVENT latch on an EXTRAPOLATED onset ( timing reading), NOT a height freeze and NOT a running maximum of height
  REAL(wp) :: inl_L(0:INL_MAXIR)      = 0.0_wp        ! cm, REALIZED internode length per inversed rank (0 = peduncle)
  REAL(wp) :: inl_m(0:INL_MAXIR)      = 0.0_wp        ! mg axis-1, structural mass bought into each rank, m_struct(x)
  REAL(wp) :: inl_dLpot(0:INL_MAXIR)  = 0.0_wp        ! cm d-1, E4 potential increment
  REAL(wp) :: inl_dLact(0:INL_MAXIR)  = 0.0_wp        ! cm d-1, E5' realized increment
  REAL(wp) :: inl_Lf_use(0:INL_MAXIR) = 0.0_wp        ! cm, the selected rec-111 final-length profile, as used today
  REAL(wp) :: inl_lam_use(0:INL_MAXIR)= 0.0_wp        ! mg cm-1, the selected rec-112 lineic profile, as used today
  REAL(wp) :: inl_tau_on(0:INL_MAXIR) = -1.0_wp       ! degCd, E1/E2 onset per inversed rank (-1 = no chronology yet)
  REAL(wp) :: inl_tau_end(0:INL_MAXIR)= -1.0_wp       ! degCd, E3 end per inversed rank (-1 = no chronology yet)
  INTEGER  :: inl_absrank(0:INL_MAXIR)= 0             ! the absolute phytomer rank each inversed rank maps to (0 = none; -1 = the peduncle)
  INTEGER  :: inl_bind(0:INL_MAXIR)   = 0             ! today's binding limb: 0 CLOSED, 1 POT, 2 MASS, 3 CAP
  REAL(wp) :: inl_h_b1 = 0.0_wp                       ! cm, BASIS 1: telescoped collar on realized L, peduncle EXCLUDED
  REAL(wp) :: inl_h_b2 = 0.0_wp                       ! cm, BASIS 2: B1 + realized peduncle length
  REAL(wp) :: inl_h_b3 = 0.0_wp                       ! cm, BASIS 3: B2 + L_ear -- the quantity handed to plant_height
  REAL(wp) :: inl_dM_struct = 0.0_wp                  ! g m-2 d-1, today's OWN structural stem increment (raw, sign kept)
  REAL(wp) :: inl_dM_offer  = 0.0_wp                  ! g m-2 d-1, MAX(0, inl_dM_struct) -- what the length pathway was offered
  REAL(wp) :: inl_dM_len    = 0.0_wp                  ! g m-2 d-1, SUM_x lambda_x dL(x) converted back to ground area -- what length actually cost
  REAL(wp) :: inl_e8_resid  = 0.0_wp                  ! g m-2 d-1, inl_dM_len - inl_dM_offer (<= 0 always; 0 when the MASS limb binds)
  REAL(wp) :: inl_e8_resid_max = 0.0_wp               ! g m-2 d-1, max |inl_e8_resid| over the MASS-limb days ONLY -- the G2'-(iii) number
  REAL(wp) :: inl_e8_over_max  = 0.0_wp               ! g m-2 d-1, max(0, inl_dM_len - inl_dM_offer) over ALL days -- must be 0: no mass may be created
  INTEGER  :: inl_days_mass    = 0                    ! days the MASS limb bound (every open rank mass-limited)
  INTEGER  :: inl_days_pot     = 0                    ! days at least one rank realized its full E4 potential (note sec.14 P3's signature)
  INTEGER  :: inl_days_thicken = 0                    ! days with inl_dM_offer > 0 and every window closed (E7 thickening)
  INTEGER  :: inl_open_days    = 0                    ! days with at least one open window
  REAL(wp) :: inl_Mstruct_prev = 0.0_wp               ! g m-2, yesterday's W_stem - W_reserve
  REAL(wp) :: inl_Mstruct_max  = 0.0_wp               ! g m-2, season-max structural stem mass -- the DECK'S OWN denominator of M_req_frac
  REAL(wp) :: inl_axes         = 0.0_wp               ! axes m-2, ct_axes_live as used today (rec 66 is an interlock)
  REAL(wp) :: inl_Lear_cm      = 0.0_wp               ! cm, L_ear in force today (0 before the peduncle window opens, or under rec 113 = 1)
  REAL(wp) :: inl_dI_use       = 0.0_wp               ! degCd, Delta_I in force (the rec-110 member)
  INTEGER  :: inl_onset_arm    = 0                    ! 1 = Fournier collar (E1), 2 = Malvoisin i-2

  ! ======================================================================== !
  ! Mineral share basis - driver record 114, accepted at the value 0, so the
  ! grain mineral share follows the dry-matter increment. The phloem-delivery
  ! basis the alternative would use is licensed as a ROUTE, not as a share, by
  ! Grusak & Graham (2005) Ann Bot 95:817, doi:10.1093/aob/mci085, and Waters,
  ! Uauy, Dubcovsky & Grusak (2009) J Exp Bot 60:4263, doi:10.1093/jxb/erp257.
  INTEGER  :: mineral_share_basis = 0    ! rec 114 (EOF-default 0 = OFF); 1 = the phloem-delivery share basis

  ! --- The delivery ledger. Each term is exported from growth_module, not
  ! invented here --- !
  REAL(wp) :: mnp_dw_leaf_day  = 0.0_wp  ! g DM m-2 d-1, P_leaf = dw_leaf (growth_module.f90:609)
  REAL(wp) :: mnp_dw_root_day  = 0.0_wp  ! g DM m-2 d-1, P_root = dw_root (growth_module.f90:610)
  REAL(wp) :: mnp_dw_stem_day  = 0.0_wp    ! g DM m-2 d-1, P_stem = the realized stem dry-matter increment, water-soluble carbohydrate excluded
  REAL(wp) :: mnp_dw_chaff_day = 0.0_wp  ! g DM m-2 d-1, P_chaff = to_chaff (growth_module.f90:2529)
  REAL(wp) :: mnp_dw_grain_day = 0.0_wp  ! g DM m-2 d-1, P_grain = to_grain (growth_module.f90:2530/:2836)
  ! --- ON-ONLY DIAGNOSTIC STATE. Identically zero when rec 114 = 0.
  REAL(wp) :: mnp_Ptot_day       = 0.0_wp ! g DM m-2 d-1, P_tot as used today
  REAL(wp) :: mnp_f_root_day     = 0.0_wp ! -, the realized shares as applied at mineral_module.f90:684-688
  REAL(wp) :: mnp_f_leaf_day     = 0.0_wp
  REAL(wp) :: mnp_f_stem_day     = 0.0_wp
  REAL(wp) :: mnp_f_chaff_day    = 0.0_wp
  REAL(wp) :: mnp_f_grain_day    = 0.0_wp
  REAL(wp) :: mnp_sumf_day       = 0.0_wp ! -, SUM_i f_i -- gate C1 says |this - 1| < 1e-12 EVERY day
  REAL(wp) :: mnp_ftot_resid_max = 0.0_wp ! -, season max |SUM_i f_i - 1|
  REAL(wp) :: mnp_legacy_dWt_day = 0.0_wp   ! g DM m-2 d-1, the dry-matter increment the incumbent law would have used; written only when the record is on
  INTEGER  :: mnp_days_share     = 0      ! days P_tot > 1e-9 and the delivery basis formed the shares
  INTEGER  :: mnp_days_fallback  = 0      ! days P_tot <= 1e-9 -> the UNCHANGED fl = 1 fallback
  INTEGER  :: mnp_days_vi        = 0            ! days on which P_tot <= 1e-9 while the incumbent increment exceeded 1e-9

  ! ======================================================================== !
  ! Canopy boundary-layer resistance and the wind input - driver records 70 to
  ! 74. Record 70 is accepted at the value 1, so the within-canopy wind
  ! attenuation of Teh, C. B. S. (2006) Introduction to Mathematical Modeling
  ! of Crop Growth, Eqs. 4.75-4.78 (p.92 and pp.95-96) is applied; record 71 is
  ! accepted at the value 0, so the daily mean wind is used, as Part 1 states.
  ! The attenuation coefficient domain is that of Teh Table 4.2 p.92, whose
  ! smallest tabulated value is 0.4 for citrus, and the factor
  ! alpha/(1 - exp(-alpha/2)) that Eq. 4.78 carries is applied explicitly.
  ! The wind measurement height is a required deck input: a negative value is
  ! refused by the reader rather than defaulted.
  INTEGER  :: rca_wind_attenuation = 0  ! runtime flag (driver_main.txt rec 70; EOF-default 0)
  INTEGER  :: wind_hourly_mode     = 0  ! runtime flag (driver_main.txt rec 71; EOF-default 0)
  ! -- Part A ---------------------------------------------------------- !
  REAL(wp) :: rwa_lai_blend = 1.0_wp    ! m2 m-2; driver rec 72. model assumption: a numerical
                                        ! domain boundary for Teh Eq. 4.75 uniform-leaf-area-density
                                        ! premise; it is not fitted and not taken from any source. It is
                                        ! exposed so that its inertness over the measured LAI range is provable.
  REAL(wp), PARAMETER :: RWA_ALPHA_MIN = 0.1_wp   ! guard only: alpha below this is outside Teh Table 4.2 p.92
                                                  ! (the smallest tabulated is citrus 0.4) and the alpha/f ratio
                                                  ! is then within 1.3 % of its alpha -> 0 limit anyway
  REAL(wp) :: rwa_att_factor  = 0.0_wp  ! alpha/(1-exp(-alpha/2)), the factor Eq. 4.78 carries and the code dropped
  REAL(wp) :: rwa_rca_ratio   = 0.0_wp  ! rca(ON)/rca(OFF) actually realized this substep
  REAL(wp) :: rwa_lim_worst   = 0.0_wp  ! worst second-order residual as a FRACTION of its analytic bound (W-3); < 1 = pass
  LOGICAL  :: rwa_lim_failed  = .FALSE. ! set by rwa_alpha0_limit_check if the alpha->0 identity is violated

  ! ======================================================================== !
  ! Full Teh canopy aerodynamic coupling r_a^c - driver record 85, accepted at
  ! the value 0, so the released configuration keeps the resistance form of the
  ! energy-balance module. The alternative is the full form of Teh (2006)
  ! pp.92-96, whose underlying source is Choudhury & Monteith (1988) Q J R
  ! Meteorol Soc 114:373-398, doi:10.1007/BF02033923.
  INTEGER  :: rca_teh_full = 0          ! runtime flag (driver_main.txt rec 85; EOF-default 0)
  REAL(wp) :: rtf_att_factor = 0.0_wp   ! DIAGNOSTIC: alpha/(1-exp(-alpha/2)) as coded this substep (rec 85 = 1)
  REAL(wp) :: rtf_rca_last   = 0.0_wp   ! DIAGNOSTIC: last-substep rca, s m-1 (rec 85 = 1)
  REAL(wp) :: rtf_rca_free_last = 0.0_wp ! DIAGNOSTIC: last-substep free-air rca, s m-1 (rec 85 = 1)
  ! -- Part B ---------------------------------------------------------- !
  REAL(wp) :: w_eb_s = 0.0_wp           ! m s-1, THE wind the energy balance uses this substep.
                                        ! = w_d_s exactly when wind_hourly_mode = 0.
  REAL(wp) :: wind_meas_hgt_m = -1.0_wp  ! m; driver rec 74. A negative value means the height was not supplied and the reader refuses it
  CHARACTER(LEN=256) :: wind_hourly_file = ''   ! driver rec 73; ignored when the flag is 0
  INTEGER, PARAMETER :: WIND_MAXREC = 60000     ! capacity guard; the WDF25 2 m package is 22752 ten-minute rows
  INTEGER  :: whr_n = 0                          ! intervals actually loaded
  INTEGER,  ALLOCATABLE :: whr_year(:), whr_doy(:)
  REAL(wp), ALLOCATABLE :: whr_hour(:), whr_u(:), whr_meas(:)   ! whr_meas: 1.0 measured, 0.0 estimated
  REAL(wp) :: whs_u_step  = 0.0_wp      ! m s-1, series mean over the current substep window
  REAL(wp) :: whs_meas_frac = 0.0_wp    ! measured fraction of the intervals in that window (W-7)
  INTEGER  :: whs_n_step  = 0           ! intervals found in that window; 0 -> daily fallback, counted
  INTEGER  :: whr_fallback_steps = 0    ! substeps that fell back to the daily mean, whole run
  INTEGER  :: whr_used_steps     = 0    ! substeps served from the series
  REAL(wp) :: whr_meas_frac_run  = 0.0_wp ! measured fraction over all served substeps

  ! ======================================================================== !
  ! Soil-surface resistance consistency - driver record 75, accepted at the
  ! value 1, so the soil-surface resistance reads the same saturated water
  ! content as the Richards solution rather than a separate value.
  INTEGER  :: rss_theta_s_consistency = 0  ! runtime flag (driver_main.txt rec 75; EOF-default 0)
  REAL(wp) :: rss_theta_s_used = 0.0_wp    ! the layer-1 theta_s actually used when the flag is on (diagnostic)
  REAL(wp) :: vgm_theta_s_L1  = 0.0_wp    ! layer-1 saturated water content AS READ FROM input/SoilLayer_vgm.csv and

  ! Published here so that the soil-surface resistance formula can read the same
  ! layer-1 quantities Richards uses.
  REAL(wp) :: N_volat_day = 0.0_wp  ! daily NH3-N volatilized (g N m-2 d-1)
  REAL(wp) :: N_denit_day = 0.0_wp  ! daily denitrified N (g N m-2 d-1)
  REAL(wp) :: N_volat_total = 0.0_wp, N_denit_total = 0.0_wp

  ! Urea hydrolysis and nitrification kinetics, first order and temperature
  ! dependent (d-1), read from input/driver_N.csv. Reference rates: Lasisi &
  ! Akinremi (2020) Can J Soil Sci, doi:10.1139/cjss-2020-0044, Table 3.
  REAL(wp) :: k_hyd_ref   = 0.12_wp        ! urea hydrolysis rate @15C (Lasisi & Akinremi 2020, CJSS, DOI 10.1139/cjss-2020-0044, Table 3); overwritten by driver
  REAL(wp) :: Q10_hyd     = 1.7_wp         ! urea hydrolysis Q10 (Lasisi & Akinremi 2020 Table 3: 0.07/0.12/0.20 d-1 at 5/15/25 C)
  REAL(wp) :: k_nit_ref   = 0.15_wp        ! nitrification rate at 25 C (PMC5364361; k1~0.186 d-1)
  REAL(wp) :: Q10_nit     = 2.5_wp         ! nitrification Q10 (Ea 61-202 kJ/mol)
  REAL(wp) :: T_ref_hyd   = 15.0_wp        ! reference temperature for k_hyd (C)
  REAL(wp) :: T_ref_nit   = 25.0_wp        ! reference temperature for k_nit (C)
  ! -------- General initial-N + organic-residue treatment module (read from driver_N.csv + --------
  ! input/N_organic_treatments.csv). Mass-conserving residue (fast/slow) + microbial pools; CENTURY-style
  ! first-order decomposition with C:N-controlled mineralization/immobilization.
  REAL(wp) :: N_soil_init_NH4 = 0.0_wp     ! g N m-2: initial soil NH4-N at sowing (NO3 init = N_soil_init)
  REAL(wp) :: residue_C_frac  = 0.42_wp    ! C fraction of residue dry mass (NRCS ~40%)
  REAL(wp) :: k_fast_ref = 0.04_wp         ! d-1 @T_ref_dec: metabolic residue decomposition (Parton 1988 ~14.8/yr)
  REAL(wp) :: k_slow_ref = 0.010_wp        ! d-1 @T_ref_dec: structural residue decomposition (Parton 1988 ~3.9/yr)
  REAL(wp) :: Q10_dec    = 2.0_wp          ! decomposition Q10
  REAL(wp) :: T_ref_dec  = 20.0_wp         ! reference soil temperature for decomposition (C)
  REAL(wp) :: theta_ref_dec = 0.30_wp      ! moisture (m3/m3) at which decomposition f_water=1
  REAL(wp) :: e_mic      = 0.40_wp         ! microbial carbon use efficiency (stoich ~0.34; range 0.25-0.55)
  REAL(wp) :: CN_mic     = 8.0_wp          ! microbial biomass C:N (CENTURY ~9; bacteria 7.4-fungi 13.4)
  REAL(wp) :: k_mic      = 0.011_wp        ! d-1 @T_ref_dec: microbial turnover (~95 d biomass turnover)
  REAL(wp) :: presow_Tsoil  = 15.0_wp      ! C: configurable pre-season mean soil temperature
  REAL(wp) :: presow_fwater = 0.7_wp       ! -: configurable pre-season moisture response factor
  ! Organic-N pools (g m-2): residue fast/slow (C and N) + microbial/immobilized (C and N).
  REAL(wp) :: res_fast_C = 0.0_wp, res_fast_N = 0.0_wp
  REAL(wp) :: res_slow_C = 0.0_wp, res_slow_N = 0.0_wp
  REAL(wp) :: soil_mic_C = 0.0_wp, soil_mic_N = 0.0_wp
  LOGICAL  :: org_treatments_active = .FALSE.   ! true if any residue treatment supplied
  ! Organic treatment schedule (read from input/N_organic_treatments.csv; ALLOCATABLE, generic)
  INTEGER  :: N_org_n = 0
  INTEGER,  ALLOCATABLE :: N_org_das(:), N_org_age(:)
  REAL(wp), ALLOCATABLE :: N_org_mass(:), N_org_Npct(:), N_org_Cpct(:), N_org_fast(:), N_org_CN(:)
  REAL(wp) :: N_soil_uptake_total = 0.0_wp ! cumulative season N uptake from soil
  REAL(wp) :: N_soil_leach_total = 0.0_wp  ! cumulative season N leaching loss
  INTEGER  :: N_fert_applied = 0           ! count of applications made so far
  ! N model parameters (read from input/driver_N.csv)
  REAL(wp) :: N_soil_init   = 2.0_wp   ! g N m-2 — initial soil mineral N
  REAL(wp) :: V_Nmax        = 0.30_wp    ! g N m-2 d-1 - not used by any statement; kept for driver_N.csv positional stability. It held a whole-crop maximum daily uptake (Barraclough 1986), which is not a root kinetic constant
  REAL(wp) :: Km_N          = 0.10_wp    ! g N m-2 - not used by any statement; kept for positional stability. The concentration-basis half-saturation Km_c_no3 is the quantity in force
  REAL(wp) :: trans_ref_N   = 1.0_wp   ! mm d-1 — reference transpiration for N scaling
  REAL(wp) :: SLN_min_N     = 0.50_wp  ! g N m-2 leaf — structural minimum SLN
  REAL(wp) :: k_remob_leaf_N = 0.04_wp ! d-1 — leaf N remobilization rate
  REAL(wp) :: k_remob_stem_N = 0.025_wp! d-1 — stem N remobilization rate
  REAL(wp) :: n_struct_stem_N = 0.005_wp! g N g-1 DM — structural stem N fraction
  REAL(wp) :: k_N_can       = 0.35_wp  ! m2 m-2 — leaf N canopy extinction coefficient
  REAL(wp) :: n_b_can       = 0.30_wp  ! g N m-2 leaf — minimum leaf N for L_N
  REAL(wp) :: f_postanth_stor = 0.70_wp! fraction — post-anthesis uptake to grain
  REAL(wp) :: f_postanth_lf  = 0.20_wp ! fraction — post-anthesis uptake to leaf
  REAL(wp) :: f_postanth_st  = 0.10_wp ! fraction — post-anthesis uptake to stem
  ! Root/shoot nitrate partition (Andrews 1986/1992; Gojon 2022 routing) — Process 2
  REAL(wp) :: f_shoot_min     = 0.20_wp ! fraction — min shoot share at low NO3 (Andrews 1986 p.511; 1992)
  REAL(wp) :: f_shoot_max     = 0.80_wp ! fraction — max shoot share at high NO3 (Andrews 1986 Table 1b; 1992)
  REAL(wp) :: K_partition_NO3 = 3.0_wp  ! mol m-3 — NO3 at half-max shoot share (Andrews 1992 transition 1-5)
  REAL(wp) :: r_root_NR       = 0.006_wp! g N g-1 DM d-1 — root NR capacity per root DM (derived, Andrews 1986 + V_Nmax whole-crop figure, Barraclough 1986 - corrected citation, see driver_N.csv)

  ! Root-scaled active nitrate uptake. The high-affinity transport maximum is
  ! the wheat root influx Vmax of Goyal & Huffaker (1986) Plant Physiol
  ! 82:1051, p.1055.
  REAL(wp) :: u_max_fw      = 2.6_wp   ! umol NO3- gFW-1 h-1 — wheat root HATS influx Vmax (Goyal & Huffaker 1986, Plant Physiol 82:1051, p.1055 from Fig.1c/d; DOI 10.1104/pp.82.4.1051)
  REAL(wp) :: fw_dw_root    = 18.6_wp  ! g FW g-1 DW — root fresh:dry ratio (root DM 5.4%; Khan et al. 2024, Heliyon 10:e29042, Table 3 control RFW 0.80 g / RDW 43 mg; DOI 10.1016/j.heliyon.2024.e29042; 35-d HYDROPONIC — closest readable same-species basis, FLAGGED)
  REAL(wp) :: Km_c_no3      = 0.027_wp ! mol m-3 — HATS half-saturation, solution-CONCENTRATION basis (Goyal & Huffaker 1986 p.1055, WHEAT); retires pool-basis Km_N
  REAL(wp) :: aroot_floor   = 0.28_wp  ! - — living-root capacity floor held beyond 20 DAA (Zheng et al. 2022, Front Plant Sci 12:814658, Tables 2-3: lower bound of RLD_frac x TTC-activity_frac at 20 DAA; DOI 10.3389/fpls.2021.814658; conservative extrapolation, stated)
  REAL(wp) :: aroot_decline_d = 0.031_wp ! d-1 - post-anthesis decline of a_root per DAY = (1 - 0.38)/20 d (Zheng 2022 Tables 2-3, composite midpoint 0.38 at 20 days after anthesis). The basis is days after anthesis because Zheng reports no site temperatures, so no thermal-time conversion can be made without inventing one

  ! Nitrogen-status factor f_Nstatus(NNI), piecewise linear on the nitrogen
  ! nutrition index. Nodes from Taulemesse, Le Gouis, Gouache, Gibon & Allard
  ! (2015) PLoS ONE 10(3):e0120291, doi:10.1371/journal.pone.0120291, Table 1
  ! p.7 (the nitrogen nutrition index at flowering of their N1, N4, N7 and N10
  ! treatments).
  REAL(wp) :: fNstat_nni(4) = (/ 0.50_wp, 1.12_wp, 1.53_wp, 1.95_wp /)   ! - — NNI nodes (Taulemesse 2015 Table 1 p.7, NNI at flowering, N1/N4/N7/N10)
  REAL(wp) :: fNstat_f(4)   = (/ 1.000_wp, 0.877_wp, 0.421_wp, 0.175_wp /) ! - — f_Nstatus at the nodes (Taulemesse 2015 p.8 uptake-flux sentence normalized to N1); f=1.000 for NNI<=0.50, f=0.175 for NNI>=1.95

  ! Layered soil mineral nitrogen and ammonium uptake: per-layer nitrate and
  ! ammonium states, upwind nitrate advection with the Richards water fluxes,
  ! and per-layer active uptake.
  INTEGER  :: soilN_layered = 1                ! driver_N.csv (positional): 1 = layered nitrate and ammonium states with upwind nitrate advection and per-layer uptake; 0 = the single-pool path
  REAL(wp) :: u_max_nh4_fw  = 3.5_wp   ! umol NH4+ gFW-1 h-1 — wheat root NH4+ transporter Vmax (Goyal & Huffaker 1986, Plant Physiol 82:1051-1056, p.1053 printed sentence: "A Km value of 0.05 mM and a Vmax of 3.5 umol (g fresh weight h)-1 was calculated for the NH4+ transporter"; DOI 10.1104/pp.82.4.1051; same intact-wheat-seedling system as the NO3 numbers)
  REAL(wp) :: Km_nh4        = 0.05_wp  ! mol m-3 — NH4+ transporter half-saturation (Goyal & Huffaker 1986 abstract p.1051 Km triple NO3/NO2/NH4 0.027/0.054/0.05 mM + p.1053)
  REAL(wp) :: n_min_depth_m = 0.20_wp  ! m — depth span receiving background mineralization, dz-weighted (FieldInfo: OM 0.4-0.55% is topsoil; a site assumption)
  REAL(wp) :: cfl_frac      = 0.5_wp   ! - — CFL guard: max fraction of a layer's NO3 advected out per sub-step (pure numerics, no physics)
  ! Layered states (g N m-2 per layer; allocated in main once NSoil_layers is known).
  ! Urea remains the scalar N_urea_soil (surface-applied, layer 1 only — design A);
  ! TAN_surf scalar semantics unchanged (Sommer 2004 volatilization).
  REAL(wp), ALLOCATABLE :: NO3_lay(:), NH4_lay(:)
  ! Daily read-only hydrology aggregates used by the nitrogen transport and
  ! uptake, filled by main_MMCO2Crop.f90 from the per-call accumulators below;
  ! the nitrogen code writes no water state:
  REAL(wp), ALLOCATABLE :: q_iface_day_mm(:)      ! signed water flux through the interface BELOW layer i (mm/d, + down); index NSoil_layers = bottom boundary (leaching interface)
  REAL(wp), ALLOCATABLE :: rootup_layer_day_mm(:) ! per-layer root water extraction (mm/d), for the mass-flow uptake term
  ! Per-hydrology-call accumulators, reset at entry and accumulated over the
  ! accepted solver substeps from the accepted state; diagnostics only, never
  ! read back by the water solver:
  REAL(wp), ALLOCATABLE :: hydro_q_iface_m(:), hydro_sink_layer_m(:)
  ! Layered-N daily diagnostics
  REAL(wp) :: U_act_NO3_day = 0.0_wp   ! realized NO3-N uptake (g N m-2 d-1)
  REAL(wp) :: U_act_NH4_day = 0.0_wp   ! realized NH4-N uptake (g N m-2 d-1)
  REAL(wp) :: N_adv_substeps_day = 0.0_wp ! advection CFL sub-steps taken today
  INTEGER  :: unit_soilN_lay = -1      ! daily_soilN_layers.csv unit (layered mode only)
  ! Grain N source-sink (Martre 2003) — Process 5
  REAL(wp) :: n_stru_grain      = 0.004_wp ! g N g-1 grain DM — structural/metabolic grain N (Martre 2003 Fig.3/Table I)
  REAL(wp) :: dvs_grain_stru_end = 1.40_wp ! dvs2 — end of cell-division/start of source-limited storage (Martre 2003 p.1963)
  ! Development-stage-dependent specific leaf area: it declines as the leaves
  ! thicken through development.
  REAL(wp) :: sla_dvs0   = 0.026_wp ! m2 g-1 — SLA at emergence (Bos 2000 Ann Bot 85:243 Table 1 upper; young thin leaves)
  REAL(wp) :: sla_dvs2   = 0.018_wp ! m2 g-1 — SLA at maturity (Araus 1993 Plant Physiol 101:185 lower; thick mature leaves)
  REAL(wp) :: sla_anchor = 0.022_wp ! m2 g-1 — profile value at anthesis (dvs2=1) = prior fixed input (consistency anchor)
  ! Fertilizer schedule (read from input/N_fertilizer.csv — ALLOCATABLE, generic)
  INTEGER  :: N_fert_n_apps = 0
  INTEGER,  ALLOCATABLE :: N_fert_DAS_arr(:)    ! DAS of each application
  REAL(wp), ALLOCATABLE :: N_fert_N_arr(:)      ! g N m-2 per application
  ! Per-application fertilizer form fractions (sum=1): NO3, NH4, urea. Default = urea.
  REAL(wp), ALLOCATABLE :: N_fert_fNO3_arr(:), N_fert_fNH4_arr(:), N_fert_furea_arr(:)

  ! Soil evaporation stage-1/2 (FAO-56 style)
  REAL(wp) :: TEW_top_mm, REW_top_mm, De_top_mm

  ! ======================================================================== !
  ! Soil humus (organic-matter) mineralization - driver record 93, accepted at
  ! the value 0, so the humus pool is not mineralized and the soil nitrogen
  ! supply is the deck initial pool plus the fertilizer and residue terms.
  ! The law the record would select is that of Brisson, Mary, Ripoche et al.
  ! (1998) "STICS: a generic model for the simulation of crops and their water
  ! and nitrogen balances. I. Theory and parameterization applied to wheat and
  ! corn", Agronomie 18(5-6):311-346, doi:10.1051/agro:19980501, Eq. (33)
  ! p.342, whose clay-content and depth terms are measured per deck and read
  ! from input/driver_N.csv rather than transplanted.
  ! The temperature response of that law has an optimum near 40 C in the
  ! published ranges (Burke et al. 2003; Davidson & Janssens 2006; Liu et al.
  ! 2018, optima 38.5-46.0 C, mean 42.4 C).
  INTEGER  :: soil_n_mech = 0                            ! runtime flag (driver_main.txt rec 93; EOF-default 0)

  ! Clay content, depth and initial humus nitrogen are MEASURED per deck
  ! (input/driver_N.csv); -1 marks a value the deck did not supply.
  REAL(wp) :: stics_argi_pct    = -1.0_wp    ! % clay, ARGI in Brisson 1998 eq. (33) p.342; MEASURED per deck (driver_N pos 65); -1 = not supplied
  REAL(wp) :: stics_calc_pct    = -1.0_wp    ! % CaCO3, CALC in eq. (33) p.342; INFERRED 0.0 on every deck, NOT MEASURED anywhere in this repo (driver_N pos 66)
  REAL(wp) :: stics_nhum0_gN_m2 = -1.0_wp    ! g N m-2, NHUM(0) over PROFHUM; MEASURED organic C / measured humus C:N 10.889 (driver_N pos 67); -1 = not supplied
  REAL(wp), PARAMETER :: stics_k2pot_num   = 0.21_wp     ! numerator of eq. (33), Brisson 1998 p.342
  REAL(wp), PARAMETER :: stics_k2pot_argi0 = 12.5_wp     ! clay offset of eq. (33), p.342
  REAL(wp), PARAMETER :: stics_k2pot_calc0 = 54.5_wp     ! limestone offset of eq. (33), p.342
  REAL(wp), PARAMETER :: stics_fth_slope   = 0.115_wp    ! eq. (42) exponent slope, p.342
  REAL(wp), PARAMETER :: stics_fth_tref    = 15.0_wp     ! degC, eq. (42) reference temperature, p.342 ("fixed at 15 degC, by convention", p.328)
  REAL(wp), PARAMETER :: stics_fh_min      = 0.2_wp      ! eq. (44) intercept, p.342
  REAL(wp), PARAMETER :: stics_fh_span     = 0.8_wp      ! eq. (44) slope, p.342
  REAL(wp), PARAMETER :: soiln_Tdec_cap_c  = 45.0_wp     ! degC, the ALREADY-APPROVED decomposition T-cap reused verbatim from soiln_module.f90 -- no new constant
  ! ---- state and diagnostics, all REPORTED in output/daily_humus_N.csv ----
  REAL(wp) :: NHUM_soil       = 0.0_wp   ! g N m-2, the humus organic-N pool (STICS NHUM). 0 until rec 93 >= 1 initializes it.
  REAL(wp) :: NHUM_soil_init  = 0.0_wp   ! g N m-2, the value it was initialized to (for the drawdown diagnostic)
  REAL(wp) :: humus_min_day   = 0.0_wp   ! g N m-2 d-1, today's SUM_L VMINH(L) (STICS VMINH)
  REAL(wp) :: humus_min_total = 0.0_wp   ! g N m-2, season sum (STICS QMINH)
  REAL(wp) :: humus_k2pot     = 0.0_wp   ! d-1, eq. (33) as evaluated in this run
  REAL(wp) :: humus_fth_day   = 0.0_wp   ! -, eq. (42) today
  REAL(wp) :: humus_fh_day    = 0.0_wp   ! -, eq. (44) today, depth-weighted mean over PROFHUM

  ! ======================================================================== !
  ! Post-anthesis spike-non-grain (chaff) sink - driver record 117, accepted at
  ! the value 0, so the chaff takes no post-anthesis sink. The published
  ! bracket the alternative would select is assembled in Slafer, Elia, Savin,
  ! Garcia, Terrile, Ferrante, Miralles & Gonzalez (2015) Food Energy Secur
  ! 4:92-109, doi:10.1002/fes3.59, folios 96-97, which relays the chaff shares
  ! of Gonzalez et al. (2011b) and Abbate et al. (2013), which are cited
  ! through Slafer et al. rather than directly. The front-loaded,
  ! first-order shape is that of Bertheloot, Martre & Andrieu (2008) Plant
  ! Physiol 148:1707-1720, doi:10.1104/pp.108.124156.
  INTEGER  :: chaff_growth_mode = 0      ! rec 117; 0 = OFF (EOF-default); 1-5 = bracket member
  ! ---- ON-ONLY STATE. Identically the initializers below when rec 117 = 0. ----
  REAL(wp) :: cgs_frac       = 0.0_wp   ! -, the selected member, set once at the driver read
  REAL(wp) :: cgs_W_anth     = 0.0_wp   ! g DM m-2, MODEL chaff at its OWN anthesis node
  REAL(wp) :: cgs_pot_tot    = 0.0_wp   ! g DM m-2, cgs_frac * cgs_W_anth
  REAL(wp) :: cgs_win        = 0.0_wp   ! dvs2, 2.0 - ear_end_dvs2
  REAL(wp) :: cgs_ddvs2      = 0.0_wp   ! dvs2 d-1, today's increment, handed in ON-only
  REAL(wp) :: cgs_demand_day = 0.0_wp   ! g DM m-2 d-1
  REAL(wp) :: cgs_taken_day  = 0.0_wp   ! g DM m-2 d-1, what the supply actually paid
  REAL(wp) :: cgs_unmet_day  = 0.0_wp   ! g DM m-2 d-1, demand - taken
  REAL(wp) :: cgs_demand_tot = 0.0_wp   ! g DM m-2
  REAL(wp) :: cgs_taken_tot  = 0.0_wp   ! g DM m-2, == the realized post-anthesis increment
  REAL(wp) :: cgs_unmet_tot  = 0.0_wp   ! g DM m-2
  REAL(wp) :: cgs_from_repro_tot   = 0.0_wp  ! g DM m-2, carbon taken from the day's reproductive stream (the grain's stream)
  REAL(wp) :: cgs_from_reserve_tot = 0.0_wp    ! g DM m-2, identically 0; printed so that the closure can be checked
  REAL(wp) :: cgs_repro_seen_tot   = 0.0_wp  ! g DM m-2, dW_repro summed over the window (the pool it competed in)
  INTEGER  :: cgs_days_win     = 0      ! days the sink was inside its window
  INTEGER  :: cgs_days_denied  = 0      ! days demand > 0 and taken < demand
  INTEGER  :: cgs_day_anth     = 0      ! model day of the anthesis node
  LOGICAL  :: cgs_anch_set     = .FALSE.

  ! ---- Double-precision closure of the chaff sink. The daily file prints the
  ! per-day residual and the run keeps its maximum ---- !
  REAL(wp) :: cgs_c6_resid   = 0.0_wp   ! g DM m-2 d-1, today's routine-level residual
  REAL(wp) :: cgs_c6_max     = 0.0_wp   ! g DM m-2 d-1, season max |cgs_c6_resid|
  REAL(wp) :: cgs_state_max  = 0.0_wp   ! g DM m-2 d-1, season max |dW_chaff - taken| on POST-ANTHESIS days only
  LOGICAL  :: cgs_post_today  = .FALSE. ! .TRUE. only on days the sink's own (post-anthesis) branch ran
  INTEGER  :: unit_cgs         = -1
  LOGICAL  :: cgs_open         = .FALSE.

  ! ======================================================================== !
  ! Organ nitrogen payment for structural growth - driver record 121, accepted
  ! at the value 2, so new structural dry matter is built at measured
  ! construction concentrations rather than at a floor.
  ! Sources: the organ structural concentrations of Bertheloot, Martre &
  ! Andrieu (2008) Plant Physiol 148:1707-1720, doi:10.1104/pp.108.124156;
  ! the shoot nitrogen-dry-matter relation of Ratjen & Kage (2013) Field Crops
  ! Res 149:338, doi:10.1016/j.fcr.2013.05.015; the dependence of
  ! photosynthetic capacity and leaf duration on leaf nitrogen in Lawlor (2002)
  ! J Exp Bot 53:773, doi:10.1093/jxb/53.370.773; and the root nitrogen
  ! requirement of Cooper & Clarkson (1989) J Exp Bot 40:753,
  ! doi:10.1093/jxb/40.7.753. Mode 2 introduces no new constant.
  INTEGER  :: n_payment_gate = 0                 ! rec 121; 0 OFF (EOF-default), 1 floors, 2 construction concentrations
  REAL(wp), PARAMETER :: npg_floor_L = 0.0079_wp ! g N g-1 lamina DM, Bertheloot et al. 2008 Plant Physiol 148:1707 structural lamina N 0.79 % DM
  REAL(wp), PARAMETER :: npg_floor_R = 0.005_wp  ! g N g-1 root DM, assumed (stem floor borrowed)
  REAL(wp), PARAMETER :: npg_build_R = 0.013_wp  ! g N g-1 root DM, Hicks 1928 New Phytol 27:108 root axis 1.29 %
  REAL(wp), PARAMETER :: npg_k_lux   = 0.10_wp   ! d-1, assumed release fraction of stem N above the structural floor
  REAL(wp) :: npg_n_assim_root_in = 0.0_wp       ! g N m-2 d-1, today's N_assim_root_day (root-reduced), set by main before allocation (ON only)
  REAL(wp) :: npg_n_assim_shoot_in = 0.0_wp      ! g N m-2 d-1, today's N_assim_shoot_day (shoot-reduced), set by main before allocation (ON only)
  REAL(wp) :: npg_pool_root_in = 0.0_wp          ! g N m-2 d-1, root-reduced N offered to the root tier first step (1)
  REAL(wp) :: npg_root_export = 0.0_wp           ! g N m-2 d-1, root-reduced N not needed by the root, exported to the shoot pool step (2)
  REAL(wp) :: npg_paid_R_root = 0.0_wp, npg_paid_R_shoot = 0.0_wp   ! g N m-2 d-1, root tier paid from root-reduced N / from the shoot pool (phloem return)
  REAL(wp) :: npg_paid_cum_R_root = 0.0_wp, npg_paid_cum_R_shoot = 0.0_wp
  LOGICAL  :: npg_live_today = .FALSE.           ! the gate ran today (vegetative window, allocation not skipped)
  REAL(wp) :: npg_c_L = 0.0_wp, npg_c_S = 0.0_wp, npg_c_R = 0.0_wp          ! tier concentrations in use, g N g-1 DM
  REAL(wp) :: npg_pool_mobile_in = 0.0_wp, npg_pool_assim_in = 0.0_wp      ! g N m-2 d-1
  REAL(wp) :: npg_lux_avail = 0.0_wp, npg_lux_used = 0.0_wp                ! g N m-2 d-1
  REAL(wp) :: npg_D_L = 0.0_wp, npg_D_S = 0.0_wp, npg_D_R = 0.0_wp          ! demands, g N m-2 d-1
  REAL(wp) :: npg_p_L = 1.0_wp, npg_p_S = 1.0_wp, npg_p_R = 1.0_wp          ! paid fractions
  REAL(wp) :: npg_paid_L = 0.0_wp, npg_paid_S = 0.0_wp, npg_paid_chaff = 0.0_wp, npg_paid_R = 0.0_wp
  REAL(wp) :: npg_pool_out = 0.0_wp, npg_mobile_left = 0.0_wp              ! handed to the ladder
  REAL(wp) :: npg_dw_leaf_pot = 0.0_wp, npg_dw_stem_pot = 0.0_wp, npg_dw_stor_pot = 0.0_wp, npg_dw_root_pot = 0.0_wp  ! g DM m-2 d-1
  REAL(wp) :: npg_dw_leaf = 0.0_wp, npg_dw_stem = 0.0_wp, npg_dw_stor = 0.0_wp, npg_dw_root = 0.0_wp                  ! realized
  REAL(wp) :: npg_unpaid_ch2o_day = 0.0_wp       ! g CH2O m-2 d-1, carbon of the unpaid increments
  REAL(wp) :: npg_dw_wsc_day = 0.0_wp            ! g DM m-2 d-1, deposited in the stem WSC reserve(a)
  REAL(wp) :: npg_resid_day = 0.0_wp             ! g CH2O m-2 d-1, explicit residual(b)
  REAL(wp) :: npg_lamina_unpaid_day = 0.0_wp     ! g DM m-2 d-1, lamina DM built as area today but not paid in N (REPORTED ledger)
  REAL(wp) :: npg_lamina_unpaid_cum = 0.0_wp     ! g DM m-2, cumulative 'area built ahead of its nitrogen' (Batch 2 closes it)
  REAL(wp) :: npg_unpaid_ch2o_cum = 0.0_wp, npg_wsc_cum = 0.0_wp, npg_resid_cum = 0.0_wp
  REAL(wp) :: npg_paid_cum_L = 0.0_wp, npg_paid_cum_S = 0.0_wp, npg_paid_cum_R = 0.0_wp
  INTEGER  :: npg_days_live = 0, npg_days_bind_L = 0, npg_days_bind_S = 0, npg_days_bind_R = 0
  REAL(wp) :: npg_ident_N_max = 0.0_wp           ! max |gate identity residual|, g N m-2

  ! ======================================================================== !
  ! Two-pool status-dependent construction - driver record 128, accepted at the
  ! value 2, so the construction pass runs with a floor and a target.
  ! The lamina floor npg_floor_L = 0.0079 g N g-1 DM is the structural lamina
  ! nitrogen concentration measured by Bertheloot, Martre & Andrieu (2008)
  ! Plant Physiol 148:1707-1720, doi:10.1104/pp.108.124156, p.1711.
  ! The shoot demand ceiling takes the marginal dilution-curve form of Brisson,
  ! Mary, Ripoche et al. (1998) Agronomie 18:311-346, doi:10.1051/agro:19980501,
  ! Eqs. (47), (50) and (51) - the formulation implemented in STICS - with the
  ! dilution coefficients of Justes, Mary, Meynard, Machet
  ! & Thelier-Huche (1994) "Determination of a critical nitrogen dilution curve
  ! for winter wheat crops", Ann Bot 74:397-407, doi:10.1006/anbo.1994.1133.
  ! The lamina mass per area used to convert between the area and the mass basis
  ! is csl_LMA = 48.5 g DM m-2 (Sales et al. 2022 J Exp Bot 73:3221, Table 1).
  INTEGER  :: n_construction_status = 0                    ! rec 128; 0 = OFF (EOF-default), 1 = the two-pass floor and target, 2 and 3 = diagnostic lamina-floor variants
  REAL(wp), PARAMETER :: ncs_floor3_SLN_martre = 0.8_wp  ! g N m-2, Martre et al. 2006 EJA 25:138, UNVERIFIED (not on disk) -- mode 3 only
  REAL(wp) :: ncs_c_floor_L_eff = 0.0_wp          ! g N g-1 lamina DM, the ACTIVE mode's lamina floor (mode 1 = npg_floor_L; the SAME value feeds pass 1, the target_eff MAX bound and the pass-2 subtraction, so every mode stays internally consistent)
  REAL(wp) :: ncs_W_t_ha = 0.0_wp                 ! above-ground DM, t ha-1, set by main before the gate call
  REAL(wp) :: ncs_W_eff = 0.0_wp                  ! = MAX(ncs_W_t_ha, 1.55) (Justes/Pask validity floor)
  REAL(wp) :: ncs_c_crit = 0.0_wp                  ! g N g-1 DM, shoot marginal critical-N rate
  REAL(wp) :: ncs_D_crit = 0.0_wp                  ! g N m-2 d-1, shoot critical-N allowance
  REAL(wp) :: ncs_c_target_L = 0.0_wp, ncs_c_target_S = 0.0_wp, ncs_c_target_R = 0.0_wp   ! g N g-1 DM, fixed targets
  REAL(wp) :: ncs_D_target_L = 0.0_wp, ncs_D_target_S = 0.0_wp                            ! g N m-2 d-1, fixed-target demands (shoot only)
  REAL(wp) :: ncs_f_crit = 1.0_wp                  ! shoot target scaling factor,
  REAL(wp) :: ncs_c_target_eff_L = 0.0_wp, ncs_c_target_eff_S = 0.0_wp, ncs_c_target_eff_R = 0.0_wp  ! g N g-1 DM, status-capped targets
  REAL(wp) :: ncs_D2_L = 0.0_wp, ncs_D2_S = 0.0_wp, ncs_D2_R = 0.0_wp                     ! g N m-2 d-1, pass-2 demand (target_eff - floor)
  REAL(wp) :: ncs_paid1_L = 0.0_wp, ncs_paid1_S = 0.0_wp, ncs_paid1_chaff = 0.0_wp, ncs_paid1_R = 0.0_wp  ! pass-1 totals, snapshot before pass 2
  REAL(wp) :: ncs_paid2_L = 0.0_wp, ncs_paid2_S = 0.0_wp, ncs_paid2_chaff = 0.0_wp, ncs_paid2_R = 0.0_wp  ! pass-2 increments
  REAL(wp) :: ncs_c_actual_L = 0.0_wp, ncs_c_actual_S = 0.0_wp, ncs_c_actual_R = 0.0_wp    ! g N g-1 DM, realized concentration
  REAL(wp) :: ncs_lamina_below_target_cum = 0.0_wp ! g DM m-2, cumulative lamina DM built with c_actual_L < c_target_eff_L
  INTEGER  :: ncs_days_fcrit = 0                   ! days f_crit < 1

  ! ======================================================================== !
  ! Spike and chaff as their own structural-nitrogen tier - driver record 138,
  ! accepted at the value 0, so the spike is served inside the existing tiers.
  ! The chaff structural concentration the alternative would use, 3.5e-3 g N
  ! g-1 DM, is measured by Bertheloot, Martre & Andrieu (2008) Plant Physiol
  ! 148:1707-1720, doi:10.1104/pp.108.124156.
  INTEGER  :: spike_n_priority = 0             ! rec 138; 0 OFF (EOF-default), 1 mode 1 (priority only, no new constant), 2 mode 2 (+ chaff floor 3.5e-3)
  REAL(wp), PARAMETER :: npg_c_chaff_mode2 = 3.5e-3_wp  ! g N g-1 DM, Bertheloot et al. 2008 Plant Physiol 148:1707 p.1711, chaff/internode/peduncle structural-N minimum (mode 2 only)
  REAL(wp) :: npg_c_chaff = 0.0_wp                            ! g N g-1 DM, the active chaff-tier floor (mode 1 = n5_stem_nstruct, mode 2 = npg_c_chaff_mode2)
  REAL(wp) :: npg_D_chaff = 0.0_wp, npg_D_stem = 0.0_wp          ! g N m-2 d-1, the two SPLIT tier demands (pass 1)
  REAL(wp) :: npg_p_chaff = 1.0_wp, npg_p_stem = 1.0_wp          ! -, the two SPLIT tier paid fractions (pass 1, against the floor only)
  REAL(wp) :: ncs_D2_chaff = 0.0_wp, ncs_D2_stem = 0.0_wp        ! g N m-2 d-1, the two SPLIT pass-2 (target) demands
  INTEGER  :: npg_days_bind_chaff = 0, npg_days_bind_stem = 0    ! days each SPLIT tier's pass-1 paid fraction < 1
  INTEGER  :: snp_viol_priority = 0             ! count of window days p_stem < 1 AND p_chaff < p_stem (priority-proof gate, MUST be 0)
  REAL(wp) :: snp_win_dw_stor_pot_sum = 0.0_wp   ! g DM m-2, SUM of the carbon-side chaff (dw_storage) potential over the Fischer pre-anthesis window
  REAL(wp) :: snp_win_p_chaff_sum = 0.0_wp, snp_win_p_stem_sum = 0.0_wp  ! -, SUM of daily p_chaff/p_stem over the window (for the window mean)
  INTEGER  :: snp_win_days = 0                   ! window days counted (dvs2 < 1, gnx_win_open, rec 138 >= 1)
  INTEGER  :: snp_win_days_pchaff_lt1 = 0, snp_win_days_pstem_lt1 = 0    ! window days p_chaff<1 / p_stem<1

  ! ======================================================================== !
  ! Diffusion to the root surface, the third nitrogen transport term - driver
  ! record 139, accepted at the value 1.
  ! Mass flow and active uptake are joined by diffusion to the root surface in
  ! the form of Roose & Fowler (2004) J Theor Biol 228:173-184,
  ! doi:10.1016/j.jtbi.2003.12.012, whose Eqs. (2.8)-(2.10) give the
  ! diffusion-plus-mass-flow transport to a root of given radius, with the root
  ! length density built from the measured specific root length
  ! nud_SRL = 120 m g-1 DM of Chen, Li, Zhang et al. (2020) BMC Plant Biol
  ! 20:185, doi:10.1186/s12870-020-02390-8, field wheat. The same STICS family
  ! of models uses the same construction (Brisson et al. 1998 Agronomie
  ! 18:311-346).
  INTEGER  :: n_uptake_diffusion = 0            ! rec 139; 0 OFF (EOF-default), 1 ON
  REAL(wp), PARAMETER :: nud_Dl_no3 = 1.9e-9_wp ! m2 s-1, free-solution diffusion coeff. of NO3- (and, reused, NH4+) in water -- a handbook-class value, no first-hand table obtained
  REAL(wp), PARAMETER :: nud_SRL    = 120.0_wp    ! m g-1 DM, specific root length; Chen et al. 2020, BMC Plant Biol 20, article 185, DOI 10.1186/s12870-020-02390-8, field wheat
  REAL(wp), PARAMETER :: nud_b_no3  = 1.0_wp    ! -, soil buffer power for NO3-/NH4+ -- an assumed default (b~theta is the non-adsorbed-ion convention, not applied; no WDF/Maricopa/Braunschweig soil-chemistry data to calibrate a site value, row-47 gap)
  REAL(wp) :: nud_s_mf_cum   = 0.0_wp           ! g N m-2, season cumulative REALIZED mass-flow-term uptake (post-cap, post-ratio)
  REAL(wp) :: nud_s_act_cum  = 0.0_wp           ! g N m-2, season cumulative REALIZED active-transporter-term uptake
  REAL(wp) :: nud_s_diff_cum = 0.0_wp           ! g N m-2, season cumulative REALIZED diffusion-term uptake
  INTEGER  :: nud_days_cap_bound = 0            ! days ANY layer's standing pool (not the rate terms mf+act+diff) was the binding cap
  REAL(wp), ALLOCATABLE :: nud_mf_lay(:), nud_diff_lay(:), nud_act_lay(:), nud_cap_lay(:), nud_real_lay(:)  ! g N m-2 d-1 per layer, today's potential mf/diff/act, the layer-pool-capped supply, and the ratio-realized uptake

  ! Record 139 mode 2 bounds the ACTIVE uptake term by the transport the soil
  ! can deliver, using the same Roose & Fowler (2004) J Theor Biol 228:173-184
  ! transport expression, so that the active term cannot exceed delivery.
  INTEGER, ALLOCATABLE :: nud_bound_lay(:)   ! 0/1 per layer, today: 1 if mode 2 bound the active term to delivery (either ion)
  INTEGER  :: nud_mode2_delivery_bound_days = 0   ! season count of days ANY layer's active term was delivery-bound, mode>=2 only

  ! Record 139 mode 3 supplies the sandy-soil impedance and buffer terms of the
  ! diffusion coefficient from measurements:
  ! - the impedance factor f_l(theta) = 1.1 (theta - theta_th), Olesen,
  !   Moldrup, Yamaguchi & Rolston (2001) "Constant slope impedance factor
  !   model for predicting the solute diffusion coefficient in unsaturated
  !   soil", Soil Sci 166:89-96, doi:10.1097/00010694-200102000-00002, Eq. (7a)
  !   p.94, fitted on 22 soils of 0-53 % clay;
  ! - the buffer power theta + rho_b Kd, with the worked value of Baldwin,
  !   Nye & Tinker (1973) Plant Soil 38:621-635, p.626.
  REAL(wp), PARAMETER :: nud_csif_slope = 1.1_wp    ! -, CSIF impedance slope, Olesen et al. 2001 Soil Sci 166:89-96 p.94 Eq.(7a), fitted 22 soils 0-53% clay
  REAL(wp), PARAMETER :: nud_th_c1 =  0.81_wp, nud_th_c2 = -0.90_wp, nud_th_c3 = -0.07_wp   ! theta_th(CF,SF,rho_b) Eq.(5) coefficients, Olesen et al. 2001 p.93
  REAL(wp), PARAMETER :: nud_th_c4 = -0.60_wp, nud_th_c5 =  0.22_wp, nud_th_c6 =  0.42_wp   ! (continued)
  REAL(wp), PARAMETER :: nud_r0 = 5.0e-4_wp         ! m (0.05 cm), root radius -- Baldwin, Nye & Tinker 1973 p.626 / Cushman 1979 Part II Fig.1, BORROWED worked-example value, non-wheat-specific
  REAL(wp) :: nud_clay_frac   = -1.0_wp   ! kg/kg, clay fraction, MEASURED/class-centroid per deck (driver_N pos 68); -1 = not supplied
  REAL(wp) :: nud_silt_frac   = -1.0_wp   ! kg/kg, silt fraction, MEASURED/class-centroid per deck (driver_N pos 69); -1 = not supplied
  REAL(wp) :: nud_bulk_density = -1.0_wp  ! Mg/m3, topsoil bulk density, MEASURED per deck (driver_N pos 70); -1 = not supplied
  REAL(wp) :: nud_theta_th    = 0.0_wp    ! m3/m3, threshold water content, computed once at init from Eq.(5)
  REAL(wp) :: nud_mode3_S_diff_cum = 0.0_wp   ! g N m-2, season cumulative REALIZED mode-3 diffusion-term share
  INTEGER  :: nud_mode3_delivery_bound_days = 0   ! season count of days ANY layer's active term was delivery-bound, mode>=3 only
  REAL(wp), ALLOCATABLE :: nud_fl3_lay(:), nud_rld3_lay(:), nud_r13_lay(:), nud_geom3_lay(:)  ! per-layer today's f_l, RLD, r1 (m), geometry factor, mode>=3 only
  REAL(wp), ALLOCATABLE :: nud_cno33_lay(:), nud_cnh43_lay(:)  ! per-layer today's solution concentration (mol m-3), mode>=3 only

  ! ======================================================================== !
  ! Mineral supply by the same three transport terms - driver record 140,
  ! accepted at the value 1, so grain minerals are supplied by mass flow
  ! (S_mf,e = Tr C_soln,e), by diffusion to the root surface with the same
  ! specific root length nud_SRL = 120 m g-1 DM (Chen et al. 2020 BMC Plant
  ! Biol 20:185), and by the element-specific active terms. The saturating
  ! zinc transport constant is a SEEDLING kinetic constant, never a field rate:
  ! Hacisalihoglu, Hart & Kochian (2001) Plant Physiol 125:456,
  ! doi:10.1104/pp.125.1.456.
  INTEGER :: mineral_supply_mech = 0   ! rec 140; 0 OFF (EOF-default), 1 ON

  ! ======================================================================== !
  ! Pre-anthesis lamina nitrogen remobilization - driver record 141, accepted
  ! at the value 1. The green lamina can export nitrogen down to its functional
  ! minimum SLN_min_N = 0.50 g N m-2 (Bertheloot, Martre & Andrieu 2008 Plant
  ! Physiol 148:1707-1720, doi:10.1104/pp.108.124156, p.1711) at the measured
  ! degradation rate 0.008 degree-Cd-1 (Bertheloot, Martre & Andrieu 2008 Funct
  ! Plant Biol 35:781-796, Table 1 p.784, doi:10.1071/FP08064). The root is a
  ! source of this nitrogen, not a sink for it (Cooper & Clarkson 1989 J Exp
  ! Bot 40:753, doi:10.1093/jxb/40.7.753).
  INTEGER  :: lamina_n_remob_pre = 0            ! rec 141; 0 OFF (EOF-default), 1 ON
  REAL(wp) :: lnr_excess_day = 0.0_wp           ! g N m-2, lamina N above SLN_min_N*cls_gai, today
  REAL(wp) :: lnr_rem_avail_day = 0.0_wp        ! g N m-2 d-1, first-order release at n2_delta_TT*dTT_day_N, today
  REAL(wp) :: lnr_rem_used_day = 0.0_wp         ! g N m-2 d-1, the part of rem_avail actually paid to stem/chaff gaps, today
  REAL(wp) :: lnr_rem_used_cum = 0.0_wp         ! g N m-2, season cumulative rem_used
  INTEGER  :: lnr_days_rem_used = 0             ! days rem_used_day > 0
  REAL(wp) :: lnr_win_rem_used_sum = 0.0_wp     ! g N m-2, SUM of rem_used over the Fischer pre-anthesis window

  ! ======================================================================== !
  ! Nitrogen-limited lamina expansion - driver record 142. The record is
  ! RESERVED at the value 0, so lamina expansion is not limited at a functional
  ! specific leaf nitrogen in this version, and no code for any other value
  ! ships. The switch itself is kept so that driver_main.txt stays positional.
  INTEGER  :: lamina_n_expansion = 0            ! rec 142; v1.0 implements only 0 (the mechanism is not part of the release); any other value is refused

  ! ======================================================================== !
  ! Whole-crop nitrogen-status-dependent construction price - driver record 143.
  ! The released configuration sets record 143 = 0, so the construction price is
  ! not modulated by the crop nitrogen status. Mode 1 is an inactive
  ! alternative, see Part 1 Table 6: it makes the price of new structural dry
  ! matter rise as the nitrogen nutrition index falls, using the same Justes et
  ! al. (1994) Ann Bot 74:397-407 dilution curve the construction ceiling
  ! already uses, and introduces no new constant.
  INTEGER  :: n_status_construction_price = 0   ! rec 143; 0 OFF (EOF-default), 1 ON
  REAL(wp) :: nsc_status_day = 0.0_wp            ! -, CLIP(nap_nni_veg,0,1) actually used today
  REAL(wp) :: nsc_price_L = 0.0_wp, nsc_price_S = 0.0_wp, nsc_price_chaff = 0.0_wp, nsc_price_R = 0.0_wp  ! g N g-1 DM, today's status-priced organ concentrations
  REAL(wp) :: nsc_cut_L_day = 0.0_wp, nsc_cut_S_day = 0.0_wp, nsc_cut_R_day = 0.0_wp   ! g DM m-2 d-1, DIAGNOSTIC-ONLY: DM refused today's relative to a floor-only counterfactual, MAX(0,p_floor,o - p_o)*dW_o,pot
  REAL(wp) :: nsc_cut_L_cum = 0.0_wp, nsc_cut_S_cum = 0.0_wp, nsc_cut_R_cum = 0.0_wp   ! g DM m-2, season cumulative
  REAL(wp) :: nsc_nni_veg_sum = 0.0_wp           ! -, SUM of nap_nni_veg over this record's own live (dvs2 < 1) days
  INTEGER  :: nsc_nni_veg_days = 0               ! live days counted in the mean above

  ! ======================================================================== !
  ! Carbon-coupled leaf-area expansion - driver record 144. The released
  ! configuration sets record 144 = 0. Mode 1 is an inactive alternative, see
  ! Part 1 Table 6: it makes the day leaf-area increment follow the carbon
  ! actually available for it rather than the potential area alone.
  ! Measured directions the alternative is judged against: the elevated-CO2
  ! leaf-area response of Wheeler, Batts, Ellis, Hadley & Morison (1996)
  ! J Agric Sci 127:37 (+30 to +45 % from 350 to 700 ppm); the measurement that
  ! elevated CO2 did NOT change tiller density in the same cultivars used here,
  ! Jiang, Mulero, Bonfil & Helman (2022) Plant Cell Environ,
  ! doi:10.1111/pce.14430; and the cross-species C3 leaf-area response of
  ! Poorter, Knopf, Wright et al. (2022) New Phytol 233:1560,
  ! doi:10.1111/nph.17802, Table 2a and Fig. 2.
  INTEGER  :: carbon_coupled_leaf_expansion = 0   ! rec 144; 0 OFF (EOF-default), 1 ON
  REAL(wp), PARAMETER :: cbc_fmax = 1.20_wp       ! -, ceiling on cbc_f; Masle 2000 Table II +15-20% cell flux, a model decision
  REAL(wp) :: cbc_f(0:CT_MAXTIL, CLR_MAXRANK)     = 1.0_wp  ! -, committed per-axis per-rank carbon bonus on cba_A_pot/a_pot, floored at 1, capped at cbc_fmax
  REAL(wp) :: cbc_fclk(0:CT_MAXTIL, CLR_MAXRANK)  = 0.0_wp  ! -, previous-day thermal-clock fraction of the rank (same convention as lew_fclk/lna_fclk)
  REAL(wp) :: cbc_surplus_day = 0.0_wp            ! -, (dm_gain-demand)/demand of the PREVIOUS day (set at the end of sink_lamina_split, rec-82 tcg_csur_day lag convention), read by today's accumulator
  REAL(wp) :: cbc_area_bonus_day = 0.0_wp         ! m2 m-2 ground, today's standing REALIZED area attributable to cbc_f > 1 (a LEVEL, not a flux)
  REAL(wp) :: cbc_area_bonus_cum = 0.0_wp         ! m2 m-2 ground, running MAX of the line above over the season -- the peak standing bonus area reached
  REAL(wp) :: cbc_f_mean_expanding = 1.0_wp       ! -, today's density-weighted mean cbc_f over ranks still expanding (0 < fclk < 1)
  REAL(wp) :: cbc_f_max_reached = 1.0_wp          ! -, running MAX of cbc_f(c,x) over every axis/rank cell, the whole season
  REAL(wp) :: cbc_laimax = 0.0_wp                 ! m2 m-2, running MAX of cls_gai (the model's own green area index) over the season
  INTEGER  :: cbc_days_surplus = 0                ! days cbc_surplus_day > 0 AND lew_fw >= 1 (the accumulator was actually incremented)
  LOGICAL  :: cbc_live_today = .FALSE.            ! diagnostic flag, this day's update_leaf_carbon_bonus call was live

  ! ======================================================================== !
  ! Bounded stem-luxury draw for the construction pass - driver record 135,
  ! accepted at the value 0, so the construction pass takes no second,
  ! luxury-nitrogen pass over the stem.
  INTEGER  :: n_payment_luxury_pass2 = 0   ! rec 135; 0 OFF (EOF-default), 1 ON
  REAL(wp) :: ncs_lux2_cap_day  = 0.0_wp   ! g N m-2 d-1, the day's bound = npg_k_lux x lux_remaining_after_pass1
  REAL(wp) :: ncs_lux2_used_day = 0.0_wp   ! g N m-2 d-1, lamina pass-2 luxury actually drawn today
  REAL(wp) :: ncs_lux2_used_cum = 0.0_wp   ! g N m-2, season cumulative
  INTEGER  :: ncs_lux2_days = 0            ! days the second tranche was drawn > 0

  ! ======================================================================== !
  ! Lamina area cannot stand ahead of its nitrogen - driver record 122,
  ! accepted at the value 2, so lamina area is both arrested at the nitrogen
  ! the lamina can carry and shed below the minimum functional specific leaf
  ! nitrogen. Direction: Lawlor (2002) J Exp Bot 53:773,
  ! doi:10.1093/jexbot/53.370.773, where a nitrogen shortage lowers both the
  ! photosynthetic capacity and the duration of the leaf.
  INTEGER  :: lai_n_arrest = 0                   ! rec 122; 0 OFF (EOF-default), 1 arrest by p_L, 2 arrest + shed below SLN_min_N
  REAL(wp) :: lna_f(0:CT_MAXTIL, CLR_MAXRANK)     = 0.0_wp  ! -, committed realized fraction of A_pot (water x N arrest)
  REAL(wp) :: lna_fprov(0:CT_MAXTIL, CLR_MAXRANK) = 0.0_wp  ! -, pre-commit fraction (today's increment at f_N = 1)
  REAL(wp) :: lna_fclk(0:CT_MAXTIL, CLR_MAXRANK)  = 0.0_wp  ! -, previous-day clock fraction of the rank
  REAL(wp) :: lna_dinc(0:CT_MAXTIL, CLR_MAXRANK)  = 0.0_wp  ! -, today's clock increment x f_w
  REAL(wp) :: lna_apot(0:CT_MAXTIL, CLR_MAXRANK)  = 0.0_wp  ! cm2, today's potential rank area
  LOGICAL  :: lna_live_today = .FALSE.
  REAL(wp) :: lna_pL_day = 1.0_wp                ! -, the p_L committed today
  REAL(wp) :: lna_inc_area_day = 0.0_wp          ! m2 m-2, today's pre-commit lamina area increment over live axes
  REAL(wp) :: lna_removed_day = 0.0_wp, lna_removed_cum = 0.0_wp        ! m2 m-2, unpaid area removed at commit
  REAL(wp) :: lna_removed_dm_day = 0.0_wp, lna_removed_dm_cum = 0.0_wp  ! g DM m-2 at csl_LMA (REPORT, not a flux)
  REAL(wp) :: lna_gai_before = 0.0_wp, lna_gai_after = 0.0_wp          ! m2 m-2, green area before / after the commit
  REAL(wp) :: lna_sl_bump(0:CT_MAXTIL) = 0.0_wp  ! leaves, persistent Sl bump per axis class (mode 2)
  REAL(wp) :: lna_sln_green_in = 0.0_wp          ! g N m-2, N_leaf / LAI_green read by the shed rule
  REAL(wp) :: lna_shed_need_day = 0.0_wp, lna_shed_day = 0.0_wp, lna_shed_cum = 0.0_wp  ! m2 m-2
  INTEGER  :: lna_days_arrest = 0, lna_days_shed = 0, lna_shed_unmet_days = 0, lna_expdead_rank_days = 0

  ! ======================================================================== !
  ! Nitrogen term of the tiller affordability test - driver record 123,
  ! accepted at the value 1, so the nitrogen cost of a new axis enters the
  ! record-99 ledger. There is no per-axis nitrogen state: the test is made on
  ! the canopy pools, and the measured direction (nitrogen shortage raises
  ! tiller mortality) is relayed by Alzueta, Abeledo, Mignone & Miralles
  ! (2012) Field Crops Res 124:248.
  INTEGER  :: tiller_n_afford = 0                 ! rec 123; 0 OFF (EOF-default), 1 N term in the rec 99 = 4 ledger
  REAL(wp) :: tna_DN(0:CT_MAXTIL) = 0.0_wp        ! g N axis-1, cumulative lamina construction N demand of the cohort
  REAL(wp) :: tna_SN(0:CT_MAXTIL) = 0.0_wp        ! g N axis-1, cumulative paid share (p_L x demand)
  REAL(wp) :: tna_phiN(0:CT_MAXTIL) = 1.0_wp      ! -, MIN(1, S_N / D_N)
  REAL(wp) :: tna_phiC(0:CT_MAXTIL) = 1.0_wp      ! -, the carbon-only phi before the MIN
  REAL(wp) :: tna_pL_used = 1.0_wp                ! -, the p_L read today (last gate)
  REAL(wp) :: tna_dN_day(0:CT_MAXTIL) = 0.0_wp    ! g N axis-1, today's demand
  INTEGER  :: tna_days_bind = 0                   ! window days on which phi_N < phi_C for at least one live cohort
  INTEGER  :: tna_days_step = 0                   ! window days stepped with the term active
  LOGICAL  :: tna_step_today = .FALSE.

  ! ======================================================================== !
  ! Nitrogen affordability and mortality priority - driver record 134, accepted
  ! at the value 0, so tiller mortality carries no separate nitrogen priority
  ! law. The nitrogen nutrition index the alternative would use is the Justes et
  ! al. (1994) Ann Bot 74:397-407 index computed in main_MMCO2Crop.f90; no
  ! source was located for a nitrogen-concentration threshold on tiller death
  ! (Backhaus et al. 2023, PMC10498016, relaying Ferrante et al. 2020, report a
  ! nitrogen DOSE, not a concentration), so no such threshold is coded.
  INTEGER  :: n_afford_priority_law = 0           ! rec 134; 0 OFF (EOF-default), 1 ON
  REAL(wp), PARAMETER :: nap_eps_gate    = 0.05_wp  ! -, no published source located, tiers 1-2 eligibility tolerance
  REAL(wp), PARAMETER :: nap_eps_gate_t3 = 0.15_wp  ! -, no published source located, tier-3 (reproductive) harder eligibility tolerance
  REAL(wp), PARAMETER :: nap_m_min       = 0.25_wp  ! -, no published source located, rate-damping floor
  REAL(wp), PARAMETER :: nap_nni_ctx_lo  = 0.60_wp  ! -, no published source located, context band lower node (max damping)
  REAL(wp), PARAMETER :: nap_nni_ctx_hi  = 0.90_wp  ! -, no published source located, context band upper node (no damping)
  REAL(wp) :: nap_num = 0.0_wp, nap_den = 0.0_wp  ! g N m-2, cumulative since window open, Component A
  REAL(wp) :: nap_phi_realized = 1.0_wp           ! -, MIN(1, nap_num/nap_den), or 1 if nap_den <= 0
  REAL(wp) :: nap_nni_veg = 1.0_wp                ! -, Component B context (vegetative-basis NNI)
  REAL(wp) :: nap_m = 1.0_wp                      ! -, rate multiplier in [nap_m_min, 1]
  LOGICAL  :: nap_gate12_open = .FALSE.           ! tiers 1-2 eligibility
  LOGICAL  :: nap_gate3_open  = .FALSE.           ! tier 3 eligibility (harder)
  REAL(wp) :: nap_tier3_cap = 1.0_wp              ! -, ceiling on gnx_phi_sp when nap_gate3_open
  REAL(wp) :: nap_f_Nstatus_day = 1.0_wp                    ! mirror of the f_Nstatus of main_MMCO2Crop.f90, set where f_Nstatus is computed
  REAL(wp) :: nap_f_Nstatus_at_gate = 1.0_wp      ! snapshot, taken before the tier-3 penalty
  REAL(wp) :: nap_N_mobile_veg_at_gate = 0.0_wp   ! snapshot, taken before the tier-3 penalty
  REAL(wp) :: nap_identity_resid = 0.0_wp         ! identically 0 by construction
  INTEGER  :: nap_days_live = 0                   ! window days this record has been active

  ! Test-day (pointwise) nitrogen-status factors, so that the affordability
  ! test reads the standing state on the day it is made rather than a cap
  ! carried over from an earlier day.
  REAL(wp) :: nap_phi_standing_testday = 1.0_wp   ! -, test-day (pointwise) MIN(1, npg_N_leaf/(n4_sln_target*cls_gai))
  REAL(wp) :: nap_m_testday            = 1.0_wp   ! -, same nap_m formula, evaluated with the test-day nap_nni_veg
  LOGICAL  :: nap_gate3_open_testday   = .FALSE.  ! tier-3 eligibility, test-day recompute
  REAL(wp) :: nap_tier3_cap_testday    = 1.0_wp      ! -, the test-day cap that gnx_set_grain_number reads
  LOGICAL  :: nap_tier3_cap_declared_inert = .FALSE.  ! TRUE only if cls_gai <= 0 on the test day (genuinely impossible; cap forced to 1, never faked)

  ! ======================================================================== !
  ! Standing basis of the affordability signal - driver record 136, accepted at
  ! the value 0, so the record-134 signal is unchanged.
  INTEGER  :: n_afford_standing_basis = 0   ! rec 136; 0 OFF, 1 = Part B (tiers 1-2 standing), 2 = + Part C (tier 3 standing)
  REAL(wp) :: nap_num_st = 0.0_wp, nap_den_st = 0.0_wp   ! g N m-2 (area-weighted), cumulative since window open, standing basis
  REAL(wp) :: nap_phi_standing = 1.0_wp     ! -, MIN(1, nap_num_st/nap_den_st), or 1 if nap_den_st <= 0
  REAL(wp) :: nap_phi_elig12 = 1.0_wp       ! -, the basis actually gating tiers 1-2 today (diagnostic)
  REAL(wp) :: nap_phi_elig3  = 1.0_wp       ! -, the basis actually gating tier 3 today (diagnostic)

  ! ======================================================================== !
  ! Nitrogen term on floret survival per spikelet - driver record 137, accepted
  ! at the value 0, so floret survival carries no nitrogen term. The measured
  ! nitrogen effect on floret survival the alternative would use is Ferrante,
  ! Savin & Slafer (2020) Field Crops Res 256:107908,
  ! doi:10.1016/j.fcr.2020.107908.
  INTEGER  :: floret_n_law = 0                      ! rec 137; 0 OFF, 1 ON
  REAL(wp) :: fnl_num = 0.0_wp, fnl_den = 0.0_wp     ! g N m-2 (area-weighted), cumulative since gnx_win_open, standing basis, THIS record's OWN sums
  REAL(wp) :: fnl_f_N = 1.0_wp                       ! -, MIN(1, fnl_num/fnl_den), or 1 if fnl_den <= 0; multiplies the hazard
  REAL(wp) :: fnl_r(0:CT_MAXTIL) = 0.0_wp            ! -, spike assimilate per primordium at the test day, per cohort, diagnostic
  REAL(wp) :: fnl_surv_c_only(0:CT_MAXTIL) = 0.0_wp  ! -, 1-EXP(-r), the carbon-only comparator, per cohort, diagnostic
  LOGICAL  :: fnl_step_today = .FALSE.               ! TRUE on days the window is open and the test has not yet fired (diagnostic gate)

  ! ======================================================================== !
  ! Root-front descent - driver record 124, accepted at the value 2, so the
  ! front descends at a measured constant rate of 9.6 mm d-1. The measured
  ! wheat penetration rates are Borg & Grimes (2007) Aust J Exp Agric 47:590,
  ! doi:10.1071/ea06071, with Gregory, McGowan, Biscoe & Hunter (1978) J Agric
  ! Sci 91:91, Barraclough & Leigh (1984) J Agric Sci 103:59 and
  ! Thorup-Kristensen, Cortasa & Loges (2009) Plant Soil 322:101 giving the
  ! 7.9-11.3 mm d-1 bracket that the three linear members span.
  INTEGER  :: root_front_law = 0                  ! rec 124; 0 Borg & Grimes (EOF-default), 1/2/3 linear 7.9/9.6/11.3 mm d-1
  REAL(wp), PARAMETER :: rfl_z0 = 0.10_wp                  ! m, emergence root-zone depth; an assumed floor
  REAL(wp), PARAMETER :: rfl_rz_mode(3) = (/ 7.9_wp, 9.6_wp, 11.3_wp /)  ! mm d-1, Kirkegaard & Lilley 2007 range: low, mid, high
  REAL(wp) :: rfl_rz = 0.0_wp                     ! mm d-1, the rate in use
  INTEGER  :: rfl_days_emerged = 0                ! days since emergence (dvs2 > 0), counted at the root-depth call
  REAL(wp) :: rfl_z_legacy = 0.0_wp, rfl_z_new = 0.0_wp, rfl_z_max = 0.0_wp    ! m, the incumbent and the new root-front depth, and the season maximum
  INTEGER  :: rfl_day_zmax = 0                    ! first model day on which z_root = z_max (0 = never)

  ! ======================================================================== !
  ! Per-kernel grain-nitrogen ceiling and potential - driver record 125. The
  ! released configuration sets record 125 = 1, so the per-kernel nitrogen
  ! ceiling applies and the kernel nitrogen potential is unchanged. Modes 2 and
  ! 3 are an inactive alternative, see Part 1 Table 6.
  ! The ceiling is the measured per-kernel nitrogen maximum of Martre, Porter,
  ! Jamieson & Triboi (2003) Plant Physiol 133:1959-1967,
  ! doi:10.1104/pp.103.030585, Table I p.1961 (1.38 mg N per kernel).
  INTEGER  :: kernel_potential_law = 0            ! rec 125; 0 OFF, 1 N ceiling, 2 ceiling + potential, 3 potential only
  REAL(wp) :: kpt_env_min = -1.0_wp, kpt_env_max = -1.0_wp   ! mg kernel-1 (= g per 1000), cultivar_yield.csv rows 7-8; -1 = absent
  REAL(wp) :: kpt_tt = 0.0_wp                     ! degCd base 0 since the grain number was set (window clock)
  REAL(wp) :: kpt_sup_cum = 0.0_wp, kpt_tsum = 0.0_wp   ! g m-2 reproductive supply over the window; SUM T over the window
  INTEGER  :: kpt_days = 0, kpt_day_set = 0, kpt_bound = 0   ! window days; day W_pot was set; 0 interior / 1 env_min / 2 env_max
  LOGICAL  :: kpt_set = .FALSE.                   ! W_pot set
  REAL(wp) :: kpt_sk_mean = 0.0_wp, kpt_tmean = 0.0_wp, kpt_dfill_days = 0.0_wp   ! mg kernel-1 d-1; degC; d
  REAL(wp) :: kpt_w_yin = 0.0_wp, kpt_w_raw = 0.0_wp, kpt_w_pot = 0.0_wp          ! mg kernel-1
  REAL(wp) :: gnc_cap = 0.0_wp, gnc_room_day = 0.0_wp, gnc_cut_day = 0.0_wp, gnc_cut_cum = 0.0_wp   ! g N m-2
  INTEGER  :: gnc_days_bound = 0                  ! days on which the ceiling cut the FPB sink

  ! ======================================================================== !
  ! Grain mineral sink bound - driver record 126, accepted at the value 1, so
  ! the grain mineral load is bounded by a measured maximum concentration. The
  ! measured grain mineral concentrations across nitrogen and genotype trials
  ! are those of Fan, Zhao, Fairweather-Tait, Poulton, Dunham & McGrath (2008)
  ! and Zhao, Su, Dunham, Rakszegi, Bedo, McGrath & Shewry (2009).
  INTEGER  :: mineral_grain_sink_cap = 0          ! rec 126; 0 OFF, 1 ON

  ! ======================================================================== !
  ! Fate of the grain-sink-full reproductive surplus - driver record 127,
  ! accepted at the value 1, so reproductive carbon the grain cannot take is
  ! carried as an explicit, labeled residual instead of disappearing.
  INTEGER  :: repro_surplus_fate = 0              ! rec 127; 0 OFF ( B-1(ii)), 1 refused surplus -> explicit residual
  REAL(wp) :: rsf_resid_day   = 0.0_wp            ! g DM m-2 d-1 refused (unlabeled) surplus booked to the residual today
  REAL(wp) :: rsf_resid_cum   = 0.0_wp            ! g DM m-2 cumulative
  REAL(wp) :: rsf_label_cum   = 0.0_wp            ! g DM m-2 cumulative surplus labeled stem WSC (stays in W_stem)
  REAL(wp) :: rsf_surplus_cum = 0.0_wp            ! g DM m-2 cumulative sink-full surplus
  INTEGER  :: rsf_days        = 0                 ! days with a refused (unlabeled) surplus > 0

  ! ======================================================================== !
  ! Fate of the sink-full surplus when the store is also full - driver record
  ! 129, accepted at the value 1 and requiring record 127 = 1: the surplus goes
  ! to the store first and the remainder is an explicit, labeled overflow loss,
  ! so that the carbon balance closes.
  INTEGER  :: sink_full_carbon_fate = 0      ! rec 129; 0 OFF (EOF-default), 1 storage-first/explicit-overflow-loss ON; REQUIRES rec 127 = 1
  REAL(wp) :: cfo_surplus_day = 0.0_wp       ! g DM m-2 d-1 today's sink-full surplus (diagnostic copy, rec 129 = 1 only)
  REAL(wp) :: cfo_wsc_day    = 0.0_wp        ! g DM m-2 d-1 today's surplus stored as WSC under the lifted window
  REAL(wp) :: cfo_wsc_cum    = 0.0_wp        ! g DM m-2 cumulative surplus stored as WSC under the lifted window (rec 129 = 1 only)
  REAL(wp) :: cfo_loss_day   = 0.0_wp        ! g DM m-2 d-1 unallocated-assimilate carbon loss today (NOT respiration, NOT stem structure)
  REAL(wp) :: cfo_loss_cum   = 0.0_wp        ! g DM m-2 cumulative
  INTEGER  :: cfo_days_loss  = 0             ! days with cfo_loss_day > 0
  INTEGER  :: cfo_days_store = 0             ! days with wsc stored under the lifted window > 0

  ! ======================================================================== !
  ! Grain nitrogen minimum concentration - driver record 130. The released
  ! configuration sets record 130 = 3, the soft structural floor below. Modes 1
  ! (cumulative) and 2 (guarded daily increment) are an inactive alternative,
  ! see Part 1 Table 6.
  INTEGER  :: grain_n_min_conc = 0                ! rec 130; 0 OFF (EOF-default), 1 cumulative, 2 guarded daily increment, 3 Martre (2003) soft structural floor; REQUIRES rec 127 = 1
  REAL(wp), PARAMETER :: gnm_min_conc = 0.015_wp  ! g N g-1 grain DM, Sinclair & Amir (1992) via Jamieson & Semenov (2000) FCR 68:21 p.23-24
  REAL(wp) :: gnm_N_grain_in      = 0.0_wp        ! g N m-2, N_grain entering today's call (= N_storage - N_chaff as of yesterday)
  REAL(wp) :: gnm_dN_grain_prev   = 0.0_wp        ! g N m-2 d-1, yesterday's grain N import (mode 2 one-day lag)
  REAL(wp) :: gnm_guard_dm        = 0.0_wp        ! g DM m-2, mode 2 guard threshold = 0.10 x kernel_number x TKW_potential/1000
  REAL(wp) :: gnm_room_dm_118     = 0.0_wp        ! g DM m-2 d-1, room diagnostic after rec 118/125 ceilings, before the floor
  REAL(wp) :: gnm_room_dm_125     = 0.0_wp        ! g DM m-2 d-1, same value (118 and 125-modes-2/3 are mutually exclusive)
  REAL(wp) :: gnm_room_dm_floor   = 0.0_wp        ! g DM m-2 d-1, the floor's own room (HUGE when rec 130 = 0 or the mode-2 guard has not been passed)
  REAL(wp) :: gnm_to_grain_before = 0.0_wp        ! g DM m-2 d-1, to_grain the existing ceilings alone would have allowed
  REAL(wp) :: gnm_to_grain_after  = 0.0_wp        ! g DM m-2 d-1, to_grain after the floor
  REAL(wp) :: gnm_refused_dm_day  = 0.0_wp, gnm_refused_dm_cum = 0.0_wp   ! g DM m-2 (d-1), the floor's OWN refused share
  LOGICAL  :: gnm_floor_bind_day  = .FALSE.
  INTEGER  :: gnm_days_bind       = 0

  ! Record 130 = 3, the soft structural grain-nitrogen floor, the released mode.
  ! Source: Martre, Porter, Jamieson & Triboi (2003) Plant Physiol
  ! 133:1959-1967, doi:10.1104/pp.103.030585, pp.1964-1965, whose regression on
  ! cv. Thesee gives the rate constant gnm_kcd = 8.44e-3 (degree-Cd)-1 from ONE
  ! treatment. Martre et al. state themselves that the structural-nitrogen
  ! concept they use comes from maize measurements (Engelen-Eigles et al. 2000),
  ! so the basis is assumed rather than wheat-specific, and the same citation
  ! pair (Martre et al. 2006; Sinclair & Amir 1992) underlies the grain-nitrogen
  ! sink of record 86.
  REAL(wp), PARAMETER :: gnm_kcd    = 8.44e-3_wp  ! (degCd)-1, Martre et al. 2003 p.1964-1965, regression fit cv. Thesee, ONE treatment
  REAL(wp), PARAMETER :: gnm_Nstru0 = 27.27e-6_wp ! g N grain-1 (= Nstru(Tt=0)), SAME single-cultivar fit
  REAL(wp), PARAMETER :: gnm_Der    = 450.0_wp    ! degCd, Martre et al. 2003 p.1964-1965, MAIZE-SOURCED (Engelen-Eigles et al. 2000)
  REAL(wp) :: gnm_Tt_degCd         = 0.0_wp       ! degCd after anthesis, the lagged gnf_tt value read this call
  REAL(wp) :: gnm_Nstru_pergrain   = 0.0_wp       ! g N grain-1, today's Eq. 3 structural/metabolic N demand per grain
  INTEGER  :: gnm_phase_day        = 0            ! 0 exponential/cell-division, 1 linear/endoreduplication, 2 frozen past Der
  REAL(wp) :: gnm_Nstru_total      = 0.0_wp       ! g N m-2, kernel_number * gnm_Nstru_pergrain
  REAL(wp) :: gnm_Nstru_total_prev = 0.0_wp       ! g N m-2, YESTERDAY's gnm_Nstru_total (one-day lag, diagnostic increment only)
  REAL(wp) :: gnm_dNstru_day       = 0.0_wp       ! g N m-2 d-1, today's structural N demand increment
  REAL(wp) :: gnm_c_stru_day       = 0.0_wp       ! dimensionless (g N g-1 DM), the Option C time-varying reference concentration
  REAL(wp) :: gnm_Nstru_paid       = 0.0_wp       ! g N m-2, cumulative MIN(structural demand increment, that day's grain N import); DIAGNOSTIC ONLY
  REAL(wp) :: gnm_Nsto             = 0.0_wp       ! g N m-2, gnm_N_grain_in - gnm_Nstru_paid; DIAGNOSTIC ONLY

  ! ======================================================================== !
  ! Site-level mineral availability - driver record 131, accepted at the value 1
  ! for iron, zinc, copper, manganese and calcium. The quantities it calibrates
  ! are soil availability parameters, and the enforced windows come from
  ! measurements:
  ! - the soil impedance factor, Barraclough & Tinker (1981) J Soil Sci
  !   32:225-236, doi:10.1111/j.1365-2389.1981.tb01702.x, p.234 Fig. 1;
  ! - the grain water-import ratio 0.7625 g H2O per g final grain and the
  !   calcium concentrations, Lizana, Riegel, Gomez, Herrera, Avila,
  !   McCormick & Calderini (2010) J Exp Bot 61:1147, doi:10.1093/jxb/erp380,
  !   Tables 2 and 3;
  ! - the xylem calcium concentration range, White & Broadley (2003) "Calcium
  !   in plants", Ann Bot 92:487, doi:10.1093/aob/mcg164.
  INTEGER :: mineral_supply_law = 0   ! rec 131; 0 OFF (EOF-default), 1 ON (Fe/Zn/Cu/Mn + Ca calibrated availability)

  ! ======================================================================== !
  ! Leaf maintenance-respiration basis - driver record 145, accepted at the
  ! value 1, so the leaf maintenance term is removed from the whole-crop
  ! maintenance sum and the leaf is charged through its own dark respiration Rd
  ! only, which is the quantity the gas-exchange formulation already resolves.
  ! No daytime light inhibition of Rd (the Kok effect; Atkin, Evans, Ball,
  ! Lambers & Pons 2000) is applied.
  INTEGER :: leaf_maintenance_basis = 0   ! rec 145; 0 OFF (EOF-default, the earlier Teh leaf term + Rd), 1 ON (leaf term removed from maint_raw, Rd only)

  ! ======================================================================== !
  ! PAR reflection basis - driver record 146, accepted at the value 1, so the
  ! PAR stream carries the single per-waveband reflection coefficient of Teh
  ! (2006) rather than the total-solar coefficient in addition to it.
  INTEGER :: par_reflection_basis = 0   ! rec 146; 0 OFF (EOF-default, the earlier double reflection), 1 ON (solar.f90:88,92 single Teh per-waveband coefficient)

  ! ======================================================================== !
  ! Nitrogen limit on retained green lamina area - driver record 147, accepted
  ! at the value 0, so the incumbent minimum-specific-leaf-nitrogen criterion
  ! runs. The criterion itself rests on SLN_min_N = 0.50 g N m-2, the
  ! end-of-grain-filling structural value of Bertheloot, Martre & Andrieu
  ! (2008) Plant Physiol 148:1707-1720, doi:10.1104/pp.108.124156, p.1711, and
  ! on the prescribed canopy value 1.5 g N m-2 of Jamieson & Semenov (2000)
  ! Field Crops Res 68:21-29, p.23, doi:10.1016/S0378-4290(00)00103-9. The
  ! measured green-area series it is judged against are those of Manderscheid,
  ! Dier, Erbs, Sickora & Weigel (2020) Open Data J Agric Res 6:28-33,
  ! doi:10.18174/odjar.v6i0.16397.
  INTEGER :: lamina_n_area_retention = 0   ! rec 147; 0 OFF (EOF-default: the incumbent SLN_min_N criterion runs verbatim, pre-change growth_module.f90:1811 / post-change :1872), 1 PRIMARY (n4_sln_target 1.5), 2 attribution arm (exponential profile), 3 refactor-identity (SLN_min_N), 4 sensitivity (0.8)
  REAL(wp), PARAMETER :: lar_sln_low   = 0.8_wp    ! g N m-2 leaf; mode 4 ONLY. Pask 2009 PhD Nottingham p.49, "a minimum SLN of 0.8 g m-2 (based on a specific leaf weight of 40 g m-2; Sinclair and Amir, 1992)"; the Sinclair & Amir 1992 primary is NOT ON DISK -- sensitivity bound only, never the primary
  REAL(wp), PARAMETER :: lar_nb_leaf   = 0.40_wp   ! g N m-2 leaf; mode 2 ONLY. The SAME NUMERAL as nb_leaf_v (assim.f90:188, a local PARAMETER not visible to GROWTH_MODULE); the wheat Vcmax25-SLN x-intercept, inside Bertheloot's 0.371-0.550 structural band. NOT a new constant
  REAL(wp) :: lar_sln_func      = 0.0_wp   ! g N m-2 leaf -- the functional SLN actually in use this run (mode 1: n4_sln_target; 3: SLN_min_N; 4: lar_sln_low; 2: unset, the profile form carries no single SLN)
  REAL(wp) :: lar_N_lam         = 0.0_wp   ! g N m-2 ground -- npg_N_leaf as read by the criterion today
  REAL(wp) :: lar_a_sup         = 0.0_wp   ! m2 m-2 -- the green area N_lam can hold (mode 1/3/4: N_lam/SLN_func; mode 2: L_N)
  REAL(wp) :: lar_need_day      = 0.0_wp   ! m2 m-2 -- MAX(0, cls_gai - lar_a_sup) before delivery
  REAL(wp) :: lar_shed_day      = 0.0_wp   ! m2 m-2 -- area the rank-ordered shed actually delivered today
  REAL(wp) :: lar_shed_cum      = 0.0_wp   ! m2 m-2 -- season sum of lar_shed_day
  REAL(wp) :: lar_gai_before    = 0.0_wp   ! m2 m-2 -- cls_gai as the criterion read it (before the rank-retirement loop)
  REAL(wp) :: lar_gai_after     = 0.0_wp   ! m2 m-2 -- cls_gai after the rank-retirement loop of the same day
  REAL(wp) :: lar_sln_in        = 0.0_wp   ! g N m-2 leaf -- N_lam / cls_gai before the shed
  REAL(wp) :: lar_sln_out       = 0.0_wp   ! g N m-2 leaf -- N_lam / cls_gai after the rank-retirement loop (same N_lam: the payment gate has not run yet)
  REAL(wp) :: lar_N_released_cum = 0.0_wp  ! g N m-2 -- resorbed N (N_carried_out - senesced_N) accumulated over the days this record's criterion fired. NOT separable from the organ law's own same-day senescence: the shed enters through Sl and the organ law returns ONE lfo_gai_died_day for the day. Labeled accordingly wherever it is printed
  INTEGER  :: lar_days_shed     = 0        ! count of days lar_shed_day > 0
  INTEGER  :: lar_days_unmet    = 0        ! count of days the delivery could not meet `need` (R6 of the plan risk register: the 60-iteration cap / Sl bounded by nt)

  ! ======================================================================== !
  ! Fate of the stem reserve - driver record 148, accepted at the value 1, so
  ! both limbs of the law run: the measured deposition limb and the measured
  ! release limb.
  ! Sources: Schnyder (1993) New Phytol 123:233-245,
  ! doi:10.1111/j.1469-8137.1993.tb03731.x, p.238, for the store and its
  ! ceiling; Gebbing & Schnyder (1999) Plant Physiol 121:871-878,
  ! doi:10.1104/pp.121.3.871, p.875, where 83 % of the total pre-anthesis carbon
  ! that reaches the grain is accounted for; and the measured stem
  ! water-soluble-carbohydrate concentrations 11.6-20.3 mg g-1 (1.2-2.0 %) of
  ! Manderscheid, Dier, Erbs, Sickora & Weigel (2020) Open Data J Agric Res
  ! 6:28-33, doi:10.18174/odjar.v6i0.16397, as the falsification target.
  INTEGER :: wsc_fate_law = 0              ! rec 148; 0 OFF (EOF-default), 1 THE LAW (term A + term B), 2 term A only, 3 term B only
  REAL(wp) :: wfl_a_refused_day  = 0.0_wp  ! g DM m-2 d-1 -- term A: sink-full surplus refused by the deposition close today (goes to cfo_loss_day)
  REAL(wp) :: wfl_a_refused_cum  = 0.0_wp  ! g DM m-2 -- season sum of wfl_a_refused_day
  INTEGER  :: wfl_a_days         = 0       ! days on which term A refused a deposition
  REAL(wp) :: wfl_b_released_day = 0.0_wp  ! g DM m-2 d-1 -- term B: pool released to cres_respired today
  REAL(wp) :: wfl_b_released_cum = 0.0_wp  ! g DM m-2 -- season sum of wfl_b_released_day
  INTEGER  :: wfl_b_days         = 0       ! days on which term B released
  INTEGER  :: wfl_gate_day       = 0       ! model day gfb_end_bound first became non-zero as seen by term B (0 = never)
  INTEGER  :: wfl_gate_bound     = 0       ! the value of gfb_end_bound on that day (1 thermal / 2 dvs2 / 3 kernel-mass / 4 rec 125)
  REAL(wp) :: wfl_pool_at_gate   = -1.0_wp ! g DM m-2 -- W_reserve on the gate day BEFORE term B's first release (-1 = gate never reached)
  INTEGER  :: wfl_days_grain_fed = 0       ! FALSIFIER READ-OUT: days term B delivered any carbon to the grain. MUST be 0 by the existing room = 0 arithmetic
  REAL(wp) :: wfl_grain_fed_cum  = 0.0_wp  ! g DM m-2 -- the mass of that, if it ever happens. MUST be 0.0

  ! The deposition-limb draw of the day, published for the record-105 ledger.
  REAL(wp) :: wfl_draw_def_day   = 0.0_wp  ! g DM m-2 d-1 -- today's deficit-limb draw, for the rec-105 = 2 I4 ledger
  REAL(wp) :: wfl_draw_def_cum   = 0.0_wp  ! g DM m-2 -- season sum (reported beside term B's released mass)
  INTEGER  :: wfl_draw_def_days  = 0       ! days the deficit limb drew under this record

  ! ======================================================================== !
  ! Surface water-balance conservation - driver record 149, accepted at the
  ! value 3. The top boundary flux is limited to what the surface can accept
  ! and the ledger residual is no longer zeroed at the minimum time step, so the
  ! true residual is carried and reported. The mixed-form mass-conservative
  ! discretization and its top-flux treatment are those of Celia, Bouloutas &
  ! Zarba (1990) Water Resour Res 26:1483, Eqs. 15-17.
  INTEGER :: wb_surface_limit = 0          ! rec 149; 0 OFF (EOF-default), 1 = no ledger zeroing at min_dt (tolerances unchanged), 2 = 1 + post-step conservative closure, 3 = 1 + consistent top control volume + no post-solve top-flux modification (conservative scheme)
  REAL(wp) :: mb_tol_step        = 1.0e-6_wp  ! m, per-substep accepted mass-balance ledger tolerance; unchanged in all modes of record 149
  REAL(wp) :: mb_tol_head_min_dt = 1.5e-5_wp ! m, Picard/Newton min-dt CONVERGENCE-FLOOR escape valve (hydro.f90:815,889) AND relaxed ledger-zeroing test (hydro.f90:1029, mode 0 only); UNCHANGED in ALL modes

  ! ======================================================================== !
  ! Phenology factor window - driver record 150, accepted at the value 1, which
  ! is the form the engine cites: McMaster, Wilhelm & Frank (2008) Ann Bot
  ! 102:561-569, in which the vernalization and photoperiod factors act as
  ! MIN(VF, DF) over the window from emergence to terminal spikelet, rather than
  ! as the product VF x DF on every day.
  INTEGER :: pheno_dev_factor_window = 0      ! rec 150; 0 = the product VF*DF on every day (EOF-default), 1 = McMaster 2008 as cited, MIN(VF,DF) from emergence to terminal spikelet only

  ! ======================================================================== !
  ! Fruiting efficiency as a cultivar trait - driver record 151, accepted at the
  ! value 1, so the cultivar fruiting-efficiency row sets the per-floret spike
  ! mass. Range: Slafer, Elia, Savin, Garcia, Terrile, Ferrante, Miralles &
  ! Gonzalez (2015) Food Energy Secur 4(2):92-109, doi:10.1002/fes3.59, where
  ! modern cultivars span about 50 to 200 grains per gram of spike, relaying
  ! Abbate et al. (2013), Gonzalez et al. (2011b) and Garcia et al. (2014).
  INTEGER  :: fe_trait_mode    = 0        ! rec 151; 0 OFF (Fischer constants for all cultivars, EOF-default), 1 = cultivar FE row sets m_floret,cv

  ! Census files - driver record 152, accepted at the value 1, which suppresses
  ! the two per-hour census files; the model state is unaffected either way.
  INTEGER  :: census_files     = 0        ! rec 152; 0 write (EOF-default), 1 suppress the two census files
  REAL(wp) :: gnx_floret_scale = 1.0_wp      ! rec 151: m_floret,cv / gnx_m_floret; exactly 1.0 when the record is 0

  ! ======================================================================== !
  ! Placement of the deck initial mineral nitrogen - driver record 153,
  ! accepted at the value 1, so the deck initial nitrogen is placed over the
  ! reported sampling depth of the deck rather than spread over the whole
  ! profile.
  INTEGER  :: n_init_placement = 0        ! rec 153; 0 = whole profile (EOF-default), 1 = reported sampling depth
  REAL(wp) :: n_init_depth_m   = -1.0_wp  ! driver_N.csv OPTIONAL row 71 N_soil_init_depth_m (m); < 0 = row absent -> whole profile

  ! ======================================================================== !
  ! Post-anthesis organ nitrogen turnover - driver record 154, accepted at the
  ! value 1, so the Bertheloot 2011a organ turnover runs, roots included.
  ! Sources: Bertheloot, Cournede & Andrieu (2011) "NEMA, a functional-
  ! structural model of nitrogen economy within wheat culms after flowering.
  ! I. Model description", Ann Bot 108:1085-1096, doi:10.1093/aob/mcr119, and
  ! Bertheloot, Wu, Cournede & Andrieu (2011) "II. Evaluation", Ann Bot
  ! 108:1097-1109, doi:10.1093/aob/mcr125, in which organ nitrogen degradation
  ! is first order as in Bertheloot et al. (2008a) and the mobile fraction of
  ! the remobilizable pool is 0.10. The construction replaces the earlier
  ! single-organ treatment, and the Eq. (6) synthesis term of Bertheloot,
  ! Martre & Andrieu (2008) Funct Plant Biol 35:781-796, doi:10.1071/FP08064,
  ! is resolved independently of it.
  INTEGER  :: organ_n_turnover_nema = 0   ! rec 154; 0 = OFF (EOF-default), 1 = Bertheloot 2011a organ turnover incl. roots, 2 = mode 1 without the root limb
  REAL(wp), PARAMETER :: n2_p_root = 0.1_wp   ! dimensionless, N influx into roots following root DM influx; Bertheloot 2011b Table 2 (DOI 10.1093/aob/mcr125) row "p r , p tp, i ... Dimensionless 0.1 Expt 1"
  ! --- rec 154 inputs handed in from main_MMCO2Crop.f90 (set every day; read only when the record is ON)
  REAL(wp) :: ont_W_root_in   = 0.0_wp    ! g DM m-2, today's root dry mass (W_root at the grain-N call)
  REAL(wp) :: ont_dw_root_day = 0.0_wp    ! g DM m-2 d-1, today's root DM increment max(0, W_root - W_root_prev_day)
  ! --- rec 154 ON-only daily diagnostics (evaluation ledgers; never fluxes ---
  REAL(wp) :: ont_pool_before = 0.0_wp, ont_uptake_credit = 0.0_wp
  REAL(wp) :: ont_D_lam = 0.0_wp, ont_D_stem = 0.0_wp, ont_D_chaff = 0.0_wp, ont_D_root = 0.0_wp
  REAL(wp) :: ont_sink_pot = 0.0_wp, ont_grain_import = 0.0_wp
  REAL(wp) :: ont_S_lam = 0.0_wp, ont_S_stem = 0.0_wp, ont_S_chaff = 0.0_wp, ont_S_root = 0.0_wp
  REAL(wp) :: ont_S_lam_pot = 0.0_wp, ont_s_scale = 1.0_wp, ont_f_N = 0.0_wp, ont_f_PAR = 0.0_wp
  REAL(wp) :: ont_dTT = 0.0_wp, ont_mob_after = 0.0_wp, ont_resid = 0.0_wp
  INTEGER  :: ont_bound_by = 0            ! 0 = sink-bound (grain got the full potential), 1 = pool-bound
  LOGICAL  :: ont_live_today = .FALSE.    ! this routine ran the turnover block today
  ! --- rec 154 season sums (run.log summary lines ---
  REAL(wp) :: ont_cum_D_lam = 0.0_wp, ont_cum_D_stem = 0.0_wp, ont_cum_D_chaff = 0.0_wp, ont_cum_D_root = 0.0_wp
  REAL(wp) :: ont_cum_S_lam = 0.0_wp, ont_cum_S_stem = 0.0_wp, ont_cum_S_chaff = 0.0_wp, ont_cum_S_root = 0.0_wp
  REAL(wp) :: ont_cum_grain = 0.0_wp, ont_cum_uptake = 0.0_wp, ont_resid_max = 0.0_wp
  REAL(wp) :: ont_cum_syn_clip = 0.0_wp
  INTEGER  :: ont_days = 0, ont_days_pool_bound = 0, ont_days_scaled = 0, ont_days_cap_clip = 0

  ! ======================================================================== !
  ! Regulation of root nitrogen-uptake demand - driver record 155, accepted at
  ! the value 1, so all three limbs run: (a) the construction-demand limb,
  ! (b) the mobile-nitrogen feedback limb and (c) the dilution-curve ceiling.
  ! Sources: the mobile-pool feedback of Bertheloot, Cournede & Andrieu (2011)
  ! Ann Bot 108:1085-1096, doi:10.1093/aob/mcr119, and Bertheloot, Wu, Cournede
  ! & Andrieu (2011) Ann Bot 108:1097-1109, doi:10.1093/aob/mcr125; the critical
  ! nitrogen dilution curve of Justes, Mary, Meynard, Machet & Thelier-Huche
  ! (1994) Ann Bot 74:397-407, doi:10.1006/anbo.1994.1133; and the demand form
  ! of Brisson, Mary, Ripoche et al. (1998) "STICS: a generic model for the
  ! simulation of crops and their water and nitrogen balances. I.", Agronomie
  ! 18:311-346, doi:10.1051/agro:19980501, which itself cites Justes et al.
  ! (1994) for the dilution curve.
  INTEGER  :: n_uptake_regulation = 0     ! rec 155; 0 = OFF (EOF-default), 1 = (a)+(b)+(c), 2 = (a) only, 3 = (b)+(c) only
  REAL(wp), PARAMETER :: n1_beta_N = 130.0_wp    ! dimensionless, E_N = exp(-beta_N [N]mob); Bertheloot, Cournede & Andrieu 2011 Ann Bot 108:1085-1096, Eq. (4) p.1090, with the value in Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109, Table 2 p.1101, row "beta C, beta N ... 2300, 130" (DOI 10.1093/aob/mcr125)
  ! --- Root-uptake constants of Bertheloot 2011b Table 2, recorded here and used by no
  ! statement of this record: they belong to the Bertheloot 2011a uptake equations, a
  ! separate construction --- !
  REAL(wp) :: nur_dwgreen_prev = 0.0_wp   ! g DM m-2 d-1, MAX(0, W_green - W_green_prev_day) of the PREVIOUS day
  ! --- rec 155 ON-only daily diagnostics (evaluation ledgers; never fluxes ---
  REAL(wp) :: nur_EN_day       = 1.0_wp   ! E_N of limb (b), dimensionless (1.0 when the limb is off)
  REAL(wp) :: nur_Nmob_day     = 0.0_wp   ! [N]mob = N_mobile_veg/(W_green+W_stem), g g-1
  REAL(wp) :: nur_cap_nmax_day = 0.0_wp   ! g N m-2, the Justes MAXIMUM-curve stock gap actually used as the outer cap (g_nmax_shoot)
  REAL(wp) :: nur_cap_nc_day   = 0.0_wp   ! g N m-2, the same stock gap on the CRITICAL curve 5.35 W^-0.442 / plateau 4.4 % -- DIAGNOSTIC ONLY, no statement reads it
  REAL(wp) :: nur_dem_before   = 0.0_wp   ! g N m-2 d-1, N_demand_day BEFORE limbs (a)/(b)
  REAL(wp) :: nur_dem_after    = 0.0_wp   ! g N m-2 d-1, N_demand_day AFTER limbs (a)/(b)
  REAL(wp) :: nur_D_max_day    = 0.0_wp   ! g N m-2 d-1, limb (a)'s maximum-curve marginal demand D_max_day
  ! --- rec 155 season sums (run.log summary line ---
  ! (season pre/post-anthesis realized uptake are the EXISTING N_uptake_pre_cum / N_uptake_post_cum; no new accumulator)
  REAL(wp) :: nur_EN_sum  = 0.0_wp, nur_EN_min = 1.0_wp     ! mean/min E_N over the live days
  REAL(wp) :: nur_cum_a_gain = 0.0_wp                       ! g N m-2, season sum of the limb-(a) demand increase actually granted
  INTEGER  :: nur_days_live = 0, nur_days_a_bind = 0, nur_days_c = 0

  ! ======================================================================== !
  ! Which nitrogen dilution curve caps construction - driver record 156,
  ! accepted at the value 1, the MAXIMUM pair. The curves are those of Justes,
  ! Mary, Meynard, Machet & Thelier-Huche (1994) Ann Bot 74:397-407,
  ! doi:10.1006/anbo.1994.1133, as used in the STICS formulation of Brisson,
  ! Mary, Ripoche et al. (1998) Agronomie 18:311-346,
  ! doi:10.1051/agro:19980501, whose demand a(1-b) W^-b IS the maximum-curve
  ! demand. The critical pair (5.35, 0.442) and the maximum pair (8.3, 0.44) are
  ! both printed by Justes et al. (1994).
  INTEGER  :: n_construction_ceiling_curve = 0   ! rec 156; 0 = OFF (EOF-default: critical pair 5.35/0.442), 1 = maximum pair 8.3/0.44
  ! --- rec 156 ON-only diagnostics: the coefficient pair actually used, written per day ---
  REAL(wp) :: ncc_a_used = 0.0_wp, ncc_b_used = 0.0_wp
  INTEGER  :: ncc_days_live = 0                  ! days the ceiling was evaluated under mode 1

  ! ======================================================================== !
  ! Luxury (storage) nitrogen in the green lamina - driver record 157, accepted
  ! at the value 3. The per-area lamina ceiling of the allocation ladder becomes
  ! the measured holding capacity rather than the prescribed canopy mean.
  ! Constants: csl_LMA = 48.5 g DM m-2 (Sales et al. 2022 J Exp Bot 73:3221,
  ! Table 1); the structural lamina floor 0.0079 g N g-1 DM (Bertheloot, Martre
  ! & Andrieu 2008 Plant Physiol 148:1707-1720, p.1711); and the measured
  ! maximum field-wheat Narea 3.8 g N m-2 (Silva-Perez et al. 2020 J Exp Bot
  ! 71(7):2299, Fig. 4B p.2305).
  INTEGER  :: lamina_n_storage = 0   ! rec 157; 0 = OFF (EOF-default: lamina ladder ceiling = n4_sln_target 1.5 g N m-2), 1 = ceiling = sln_max_N_lnh 2.5 g N m-2 (Pask 2009 p.29) with the incumbent order (lamina storage before stem luxury), 2 = same ceiling with stem/chaff luxury paid BEFORE lamina storage (documented alternative order), 3 = mode 1 PLUS the stem/chaff luxury N above the structural floor and the shoot nitrate store counted in the record-155 limb-(b) [N]mob numerator (the lamina store is EXCLUDED, on Bertheloot 2011a's own three-pool sentence)
  ! --- rec 157 ON-only diagnostics (ledgers, never fluxes ---
  REAL(wp) :: lns_sln_ceil_used = 0.0_wp   ! g N m-2 leaf, the per-area lamina ladder ceiling actually used today
  REAL(wp) :: lns_cap_day       = 0.0_wp   ! g N m-2, today's lamina storage headroom offered to the tier
  REAL(wp) :: lns_stored_day    = 0.0_wp   ! g N m-2, N actually stored in the lamina by the tier today
  REAL(wp) :: lns_stored_cum    = 0.0_wp   ! g N m-2, season sum of the above
  LOGICAL  :: lns_bind_day      = .FALSE.  ! .TRUE. on a day the tier took the whole headroom (capacity binding)
  INTEGER  :: lns_days_live = 0, lns_days_bind = 0, lns_days_store = 0

  ! Record 157 mode 3 also fixes what enters the numerator of the record-155
  ! mobile-nitrogen feedback. The estimate used is those papers' own: Bertheloot,
  ! Cournede & Andrieu (2011) Ann Bot 108(6):1085-1096,
  ! doi:10.1093/aob/mcr119, and Bertheloot, Wu, Cournede & Andrieu (2011) Ann
  ! Bot 108(6):1097-1109, doi:10.1093/aob/mcr125, which state verbatim that
  ! "mobile N mass was estimated as in Bertheloot et al. (2008a)" as a fixed
  ! fraction of the remobilizable nitrogen mass.
  REAL(wp) :: lns_nmob_num_day   = 0.0_wp  ! g N m-2, the [N]mob NUMERATOR actually used by rec 155 limb (b) today
  REAL(wp) :: lns_nmob_t_mobile  = 0.0_wp  ! g N m-2, term (i) N_mobile_veg
  REAL(wp) :: lns_nmob_t_no3     = 0.0_wp  ! g N m-2, term (ii) N_NO3_shoot (mode 3 only)
  REAL(wp) :: lns_nmob_t_stemlux = 0.0_wp  ! g N m-2, term (iii) stem N above the floor (mode 3 only)
  REAL(wp) :: lns_nmob_t_chaflux = 0.0_wp  ! g N m-2, term (iv) chaff N above the floor (mode 3 only)

  ! ======================================================================== !
  ! Basis of the nitrate drain rate - driver record 158, accepted at the value
  ! 1, so the nitrogen-supply term of the photorespiratory nitrogen assimilation
  ! reads the shoot de novo nitrogen-reduction RATE rather than the standing
  ! shoot nitrate store. Source of the supply term: Busch, Sage & Sage (2018)
  ! "Plants increase CO2 uptake by assimilating nitrogen via the
  ! photorespiratory pathway", Nature Plants 4:46-54,
  ! doi:10.1038/s41477-017-0065-x, which defines N_max as a rate.
  INTEGER  :: nitrate_drain_rate_basis = 0   ! rec 158; 0 = OFF (EOF-default: Busch N_max supply term = the STANDING shoot nitrate store N_NO3_shoot divided by the photoperiod and LAI, main_MMCO2Crop.f90:7622), 1 = the supply term is the model's own daily de novo nitrate-reduction DEMAND in the shoot (yesterday's N paid into lamina + stem + chaff protein plus the grain's draw on newly assimilated N), same unit chain, same Busch enzymatic ceiling N_max_ref above it (Busch et al. 2018 Nature Plants 4:46-54 p.49; SHARED.f90)

  ! --- Yesterday's shoot de novo nitrogen-reduction demand (g N m-2 d-1) --- !
  REAL(wp) :: ndb_rate_prev = 0.0_wp
  ! --- rec 158 ON-only diagnostics (ledgers, never fluxes ---
  REAL(wp) :: ndb_nmax_store_day = 0.0_wp  ! umol m-2leaf s-1, the store-basis supply term mode 0 would have used today
  REAL(wp) :: ndb_nmax_rate_day  = 0.0_wp  ! umol m-2leaf s-1, the rate-basis supply term mode 1 used today
  REAL(wp) :: ndb_ag_last        = 0.0_wp  ! -, last alpha_G computed in assim.f90 this day (Busch Eq. 18)
  REAL(wp) :: ndb_as_last        = 0.0_wp  ! -, last alpha_S computed in assim.f90 this day (Busch Eq. 19)
  REAL(wp) :: ndb_nmax_used_cum  = 0.0_wp  ! umol m-2leaf s-1, running sum of N_max_param for the season mean
  REAL(wp) :: ndb_nmax_store_cum = 0.0_wp  ! umol m-2leaf s-1, running sum of the store-basis term, same days
  REAL(wp) :: ndb_rate_cum       = 0.0_wp  ! g N m-2, season sum of the numerator flux
  INTEGER  :: ndb_days_live = 0            ! days the N_max block ran
  INTEGER  :: ndb_days_cap  = 0            ! days the Busch enzymatic ceiling N_max_ref bound instead of the supply term

  ! ======================================================================== !
  ! Mobile-nitrogen estimate - driver record 159, accepted at the value 0, so
  ! the record-155 mobile-nitrogen numerator stays whatever record 157 leaves
  ! it. The estimate the alternative would install is Bertheloot, Wu,
  ! Cournede & Andrieu (2011) Ann Bot 108(6):1097-1109,
  ! doi:10.1093/aob/mcr125, and Bertheloot, Cournede & Andrieu (2011) Ann Bot
  ! 108(6):1085-1096, doi:10.1093/aob/mcr119, section "Plant structure", where
  ! the mobile nitrogen mass is a fixed fraction of the remobilizable mass.
  INTEGER  :: mobile_n_nema_estimate = 0   ! rec 159; 0 = OFF (EOF-default: the rec-155 limb-(b) [N]mob numerator is whatever rec 157 formed), 1 = numerator = n2_mobile_frac * (the rec-154 remobilizable/N_ph pools of lamina + stem + chaff), Bertheloot 2011b's printed estimate "mobile N ... represented 10 % of the remobilizable N mass" (Bertheloot et al. 2011 Ann Bot 108:1097, DOI 10.1093/aob/mcr125; SHARED.f90), 2 = mode 1 PLUS the explicit transport pools N_mobile_veg + N_NO3_shoot (the DOCUMENTED alternative, which double-counts that model's single mobile pool)
  ! --- rec 159: the mobile fraction. ALIAS, not a new number ( ---
  REAL(wp), PARAMETER :: n2_mobile_frac = gnf_r_mob_ph   ! dimensionless, mobile N / remobilizable N = 0.10; Bertheloot 2011b (Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109, DOI 10.1093/aob/mcr125) "Mobile N mass was estimated as in Bertheloot et al. (2008 a ), assuming that it represented 10 % of the remobilizable N mass"; numerically and by provenance the SAME constant as gnf_r_mob_ph (Bertheloot et al. 2008 Funct Plant Biol 35:781 Eqn 1 p.783, value p.786)
  ! --- rec 159 ON-only diagnostics (ledgers, never fluxes ---
  REAL(wp) :: mne_nrem_lam   = 0.0_wp  ! g N m-2, remobilizable (N_ph) N of the lamina
  REAL(wp) :: mne_nrem_stem  = 0.0_wp  ! g N m-2, remobilizable (N_ph) N of the stem
  REAL(wp) :: mne_nrem_chaff = 0.0_wp  ! g N m-2, remobilizable (N_ph) N of the chaff
  REAL(wp) :: mne_nrem_day   = 0.0_wp  ! g N m-2, their sum
  REAL(wp) :: mne_num_day    = 0.0_wp  ! g N m-2, the [N]mob numerator this record installed today
  REAL(wp) :: mne_num_inc_day = 0.0_wp ! g N m-2, the numerator record 157 would have used today (the quantity replaced)
  REAL(wp) :: mne_num_cum    = 0.0_wp  ! g N m-2, running sum of mne_num_day over live days (for the season mean)
  REAL(wp) :: mne_num_inc_cum = 0.0_wp ! g N m-2, running sum of mne_num_inc_day over the same days
  INTEGER  :: mne_days_live  = 0       ! days this record formed a numerator

  ! ======================================================================== !
  ! Outer nitrogen cap below the dilution floor - driver record 160, accepted at
  ! the value 1, so the cap is the Justes et al. (1994) maximum curve evaluated
  ! at the floor rather than an undocumented plateau. Source: Justes, Mary,
  ! Meynard, Machet & Thelier-Huche (1994) Ann Bot 74:397-407,
  ! doi:10.1006/anbo.1994.1133.
  INTEGER  :: justes_nmax_floor_basis = 0   ! rec 160; 0 = OFF (EOF-default: the Nmax outer cap is the undocumented 6.6 % plateau below 1.55 t ha-1), 1 = the Nmax outer cap below the floor is held at the printed curve's own value at the floor, 8.3 * 1.55**(-0.44) = 6.844 %, the model's stated continuity construction (Justes et al. 1994 Ann Bot 74:397-407 DOI 10.1006/anbo.1994.1133 prints Eqn (17) Nmax = 8.3 DM^-0.44 on p.405, valid 1.5-14 t ha-1, and NO Nmax statement below 1.5 t ha-1; SHARED.f90)
  ! --- rec 160 ON-only diagnostics (counters, never fluxes ---
  REAL(wp) :: jnf_cap_pct    = 0.0_wp  ! %DM, the cap value mode 1 installed below the floor (8.3*1.55^-0.44)
  INTEGER  :: jnf_days_floor = 0       ! days the sub-floor branch (shoot DM <= 1.55 t ha-1) was taken
  INTEGER  :: jnf_days_bind  = 0       ! of those, days the Nmax stock gap actually bound the day's N demand

  ! ======================================================================== !
  ! Post-anthesis uptake regulation - driver record 163. The released
  ! configuration sets record 163 = 0, so post-anthesis potential uptake stays
  ! MIN(construction demand, supply). Mode 1 is an inactive alternative, see
  ! Part 1 Table 6: it removes that gate after anthesis, following the
  ! post-flowering uptake description of Bertheloot, Cournede & Andrieu (2011)
  ! Ann Bot 108(6):1085-1096, doi:10.1093/aob/mcr119.
  INTEGER  :: postanth_uptake_regulation = 0   ! rec 163; 0 = OFF (EOF-default: post-anthesis U_pot = MIN(construction demand, supply)), 1 = after anthesis the demand gate is removed and uptake potential = the existing soil-N supply law x the record-155 E_N, as Bertheloot 2011a Eqn 3-5 prints it with no demand term (Bertheloot, Cournede & Andrieu 2011 Ann Bot 108:1085, DOI 10.1093/aob/mcr119; SHARED.f90)
  ! --- rec 163 ON-only diagnostics (ledgers and counters, never fluxes ---
  REAL(wp) :: pur_gain_day   = 0.0_wp  ! g N m-2 d-1, today's increase in U_pot the removed gate granted
  REAL(wp) :: pur_cum_gain   = 0.0_wp  ! g N m-2, season sum of pur_gain_day
  REAL(wp) :: pur_EN_sum     = 0.0_wp  ! sum of E_N over the post-anthesis days the statement ran
  REAL(wp) :: pur_EN_min     = 1.0_wp  ! minimum E_N over those days
  INTEGER  :: pur_days_post  = 0       ! post-anthesis days on which the gate was removed
  INTEGER  :: pur_days_bind  = 0       ! of those, days on which the removal actually raised U_pot

  ! ======================================================================== !
  ! Lamina capacity on a profile basis - driver record 164. The released
  ! configuration sets record 164 = 0, so the record-157 ceiling is applied to
  ! the canopy MEAN specific leaf nitrogen. Mode 1 is an inactive alternative,
  ! see Part 1 Table 6: it applies the ceiling at the TOP of the profile, using
  ! the n_top/n_mean factor of the same canopy profile function. The canopy-mean
  ! basis of the published ceiling is Bertheloot, Martre & Andrieu (2008) Plant
  ! Physiol 148:1707-1720, doi:10.1104/pp.108.124156, Fig. 2A p.1711, whose
  ! axis is the mean canopy leaf nitrogen per unit area.
  INTEGER  :: lamina_capacity_profile_basis = 0   ! rec 164; 0 = OFF (EOF-default: the record-157 ceiling sln_max_N_lnh = 2.5 g N m-2 is applied as the CANOPY-MEAN lamina capacity), 1 = the same printed 2.5 g N m-2 (Pask 2009 Nottingham PhD p.29, "especially in the upper leaves of the canopy") applied at the TOP of the model's own existing exponential canopy N profile, so the canopy-mean capacity = 2.5 / (n_top/n_mean) from that profile's own daily state
  ! --- rec 164 ON-only diagnostics (ledgers and counters, never fluxes ---
  REAL(wp) :: lcp_factor_day    = 0.0_wp  ! n_top/n_mean actually used today (0 on a fallback day)
  REAL(wp) :: lcp_cap_day       = 0.0_wp  ! g N m-2 leaf, the canopy-mean lamina capacity in force today
  REAL(wp) :: lcp_factor_sum    = 0.0_wp  ! sum of lcp_factor_day over the days the profile basis was used
  INTEGER  :: lcp_days_profile  = 0       ! days the profile factor was used
  INTEGER  :: lcp_days_fallback = 0       ! days the mode-0 value was used because the profile was not defined

  ! ======================================================================== !
  ! Delay of the reserve release - driver record 165, accepted at the value 1,
  ! so the deficit limb of the reserve draw opens only after the measured
  ! accumulation window has closed, instead of at anthesis.
  ! Sources: Schnyder (1993) New Phytol 123:233-245,
  ! doi:10.1111/j.1469-8137.1993.tb03731.x, p.238, which states that the store
  ! accumulates in the stem and is drawn later; Xue, Zhu, Yang et al. (2006)
  ! "Molecular dissection of variation in carbohydrate metabolism related to
  ! water-soluble carbohydrate accumulation in stems of wheat", Crop Science
  ! 46:735-746, doi:10.2135/cropsci2005.04-0033; and the measured stem-reserve
  ! series of Manderscheid, Dier, Erbs, Sickora & Weigel (2020) Open Data
  ! J Agric Res 6, doi:10.18174/odjar.v6i0.16397, which the delay is judged
  ! against.
  INTEGER  :: wsc_release_delay_basis = 0   ! rec 165; 0 = OFF (EOF-default: the rec-92 deficit draw limb opens at anthesis, dvs2 >= reserve_start_dvs2, with no post-anthesis delay), 1 = the limb additionally requires the EXISTING conjunct wsc_day_anth > 0 .AND. (iday_in - wsc_day_anth) > wsc_accum_end_dpa -- the identical conjunct the rec-105 = 2 release limb already carries -- so mobilization cannot start before the printed 10-20 d accumulation window ends (Schnyder 1993 New Phytol 123:233 p.238; SHARED.f90)
  ! --- rec 165 ON-only diagnostics (counters and ledgers, never fluxes ---
  INTEGER  :: wrd_days_delayed   = 0       ! days the deficit limb WOULD have fired and was held by the accumulation window
  REAL(wp) :: wrd_room_held_cum  = 0.0_wp  ! g DM m-2, cumulative unmet grain demand the pool was not allowed to serve on those days
  REAL(wp) :: wrd_pool_held_max  = 0.0_wp  ! g DM m-2, largest W_reserve seen on a held day

  ! ======================================================================== !
  ! Grain calcium by xylem throughput - driver record 166, accepted at the value
  ! 1, so grain calcium is the calcium arriving in the ear xylem stream.
  ! Calcium is phloem-immobile and is not redistributed from the vegetative
  ! organs to the grain: White & Broadley (2003) "Calcium in plants", Ann Bot
  ! 92:487, doi:10.1093/aob/mcg164; Hocking (1994) J Plant Nutr 17:1289,
  ! doi:10.1080/01904169409364807, Table 4 p.1296, where calcium is below the
  ! paper's own 10 % redistribution threshold. The xylem calcium concentration
  ! and the ear water flux are the two measured factors (Lizana et al. 2010
  ! J Exp Bot 61:1147; Martin & McGowan 1979 J Agric Sci 93:485,
  ! doi:10.1017/S0021859600038268; Wolterbeek, van Luipen & de Bruin 1992
  ! J Exp Bot 43:1315, doi:10.1093/jxb/43.10.1315, Table 3 p.1319).
  INTEGER  :: grain_ca_xylem_throughput = 0  ! rec 166; 0 = OFF (EOF-default: grain Ca = grain_h2o_g_gDM x xylem_Ca_mgL_CAL, a product of two constants), 1 = the ear xylem-throughput route of.. with the campaign-calibrated xylem_Ca_mgL_CAL NOT read
  ! --- rec 166 ON-only diagnostics (ledgers, never fluxes ---
  REAL(wp) :: gcx_v_ear_cum    = 0.0_wp   ! L m-2, cumulative ear xylem water V_ear
  REAL(wp) :: gcx_v_net_cum    = 0.0_wp   ! L m-2, cumulative Lizana net grain water import
  REAL(wp) :: gcx_ca_grain_cum = 0.0_wp   ! mg m-2, cumulative Ca entering the grain by this route
  REAL(wp) :: gcx_ca_chaff_cum = 0.0_wp   ! mg m-2, cumulative Ca entering the chaff by this route
  REAL(wp) :: gcx_f_ear_max    = 0.0_wp   ! -, largest daily ear share of green area
  INTEGER  :: gcx_days_route   = 0        ! days the route moved Ca into the ear

  ! ======================================================================== !
  ! Basis of the lamina synthesis term - driver record 167, accepted at the
  ! value 2. The sheath mass is moved from the stem basis to the lamina basis,
  ! because the published synthesis rate was measured on the lamina-plus-sheath
  ! entity: Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot 108:1097-1109,
  ! doi:10.1093/aob/mcr125, Table 2 p.1101, and Bertheloot, Martre & Andrieu
  ! (2008) Funct Plant Biol 35:781-796, doi:10.1071/FP08064, Eqs. (5) and (6)
  ! p.784, with the structural lamina concentration 0.0079 g N g-1 DM measured
  ! by Bertheloot, Martre & Andrieu (2008) Plant Physiol 148:1707-1720,
  ! doi:10.1104/pp.108.124156, p.1711.
  INTEGER  :: lamina_synthesis_sheath_basis = 0   ! rec 167; 0 = OFF (EOF-default: the Bertheloot 2011a sigma_lam = 1.5e-4 keeps acting on the green LAMINA mass alone and sigma_stem = 3.0e-5 on the whole lumped W_stem, which contains the sheaths), 1 = M_basis = W_green + cse_sai*(W_green/lai) (area-derived sheath mass, model's own state, rec 68 = 1 required), 2 = M_basis = W_green*(1 + 272/420) from the PRINTED MEASURED masses of Bertheloot 2011b Table 1c p.1100; in both ON modes the same M_sheath is removed from the stem's sigma basis so no mass synthesizes twice
  REAL(wp), PARAMETER :: lsb_a_sheath = 0.6476190476190476_wp    ! -, sheath/lamina dry-mass ratio at flowering = 272/420; Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109, Table 1(c) p.1100 (DOI 10.1093/aob/mcr125): sheath structural 103+69+43+9 = 224 plus remobilizable 22+15+9+2 = 48 = 272 mg culm-1, over lamina structural 140+90+50+40 = 320 plus remobilizable 40+40+20+0 = 100 = 420 mg culm-1. A printed measured allometry, held fixed and never fitted
  ! --- rec 167 ON-only diagnostics (counters and ledgers, never fluxes ---
  INTEGER  :: lsb_days_applied     = 0       ! days the rebasing actually multiplied a non-zero lamina synthesis term
  INTEGER  :: lsb_days_no_area     = 0       ! mode 1 only: days on which lai <= 0 or cse_sai <= 0, so the mode-0 basis stood
  REAL(wp) :: lsb_msheath_cum      = 0.0_wp  ! g DM m-2, cumulative M_sheath actually used (divide by lsb_days_applied for the mean)
  REAL(wp) :: lsb_ratio_last       = 0.0_wp  ! -, M_sheath/W_green on the last day the rebasing applied
  REAL(wp) :: lsb_syn_lam_gain_cum = 0.0_wp  ! g N m-2, cumulative EXTRA lamina synthesis the rebasing asked for (before s_scale and before the capacity clip)
  REAL(wp) :: lsb_syn_stem_loss_cum = 0.0_wp ! g N m-2, cumulative stem synthesis the rebasing removed (before s_scale)

  ! --- Daily diagnostics appended to the lamina-nitrogen file when the record
  ! is on --- !
  REAL(wp) :: lsb_m_sheath_day     = 0.0_wp  ! g DM m-2, M_sheath actually moved from the stem basis into the lamina basis TODAY (0 on a day the rebasing did not apply)
  REAL(wp) :: lsb_ratio_day        = 1.0_wp  ! -, REALIZED M_basis/W_green TODAY (1 = the mode-0 lamina-only basis stood; mode 2 = 1 + 272/420 = 1.6476 by construction)

  ! ======================================================================== !
  ! Mineral diffusion impedance - driver record 168, accepted at the value 1, so
  ! the impedance factor of the mineral diffusion term is the measured one:
  ! Olesen, Moldrup, Yamaguchi & Rolston (2001) Soil Sci 166:89-96,
  ! doi:10.1097/00010694-200102000-00002, p.94, with the independent wheat-soil
  ! measurement of Barraclough & Tinker (1981) J Soil Sci 32:225-236,
  ! doi:10.1111/j.1365-2389.1981.tb01702.x, p.234 Fig. 1, as the corroboration.
  INTEGER  :: mineral_impedance_basis = 0   ! rec 168; 0 = OFF (EOF-default: f_l(theta) = MAX(0, 1.60*vwc_root - 0.17), the liquid impedance factor of Barraclough & Tinker 1981 J Soil Sci 32:225-236 DOI 10.1111/j.1365-2389.1981.tb01702.x Fig. 1, printed 1.58*theta - 0.17 (sieved Woburn soil; 0.99*theta - 0.17 at Rothamsted), site-blind), 1 = Olesen et al. 2001 Soil Sci 166:89-96 p.94 Eq.(7a)/(7b) CSIF, f_l = nud_csif_slope * MAX(0, vwc_root - nud_theta_th) with nud_theta_th from the deck's own clay/silt/bulk density by p.93 Eq.(5) -- MINERAL SUPPLY ONLY, record 139's nitrogen path untouched
  ! --- rec 168 ON-only diagnostics (ledgers, never fluxes ---
  REAL(wp) :: mib_fl_sum       = 0.0_wp   ! -, sum of the CSIF f_l over the element-days the live diffusion locus ran (mean = /mib_fl_days)
  REAL(wp) :: mib_fl_off_sum   = 0.0_wp   ! -, sum of the mode-0 pair over the SAME days, so the ratio is reported and not asserted
  INTEGER  :: mib_fl_days      = 0        ! ELEMENT-days the live diffusion locus computed an impedance at mode 1 (the element loop runs 5 non-Ca elements per day)
  INTEGER  :: mib_days_at_zero = 0        ! element-days theta <= theta_th, i.e. Olesen Eq.(7b) returned exactly 0
  REAL(wp) :: mib_vwc_min      = 1.0e30_wp ! m3/m3, smallest vwc_root seen on those days
  REAL(wp) :: mib_vwc_max      = -1.0e30_wp ! m3/m3, largest vwc_root seen on those days

  ! ======================================================================== !
  ! Zinc uptake kinetics - driver record 169, accepted at the value 1, so the
  ! saturating zinc transport term is bounded by the transport the soil
  ! delivers. The apparent Michaelis constant is the seedling value of
  ! Hacisalihoglu, Hart & Kochian (2001) Plant Physiol 125:456,
  ! doi:10.1104/pp.125.1.456, and the whole-season FIELD uptake capacity that
  ! scales it is Liu, Liu, Zhang, Chen & Zou (2019) Front Plant Sci 10:426,
  ! doi:10.3389/fpls.2019.00426, Tables 3, 4 and 6; the two numbers come from
  ! two papers and are cited separately.
  ! The grain-zinc regression of Carducci et al. (2025) is a falsification
  ! target for the resulting concentration, never a tuning target. The
  ! transporter route itself is described by Gojon, Cassan, Bach, Lejay & Martin
  ! (2022) Trends Plant Sci 27:882, doi:10.1016/j.tplants.2022.04.007.
  INTEGER  :: zn_uptake_kinetic_basis = 0   ! rec 169; 0 = OFF (EOF-default: S_hats = MIN(S_hats_raw, S_mf + S_diff), a MIN that binds on 7499/7499 active deck-days so zinc supply is exactly 2 x (mass flow + diffusion)), 1 = the surface-concentration kinetic of, Liu 2019 field capacity x Hacisalihoglu 2001 Table II p.459 high-affinity Km, a assumed CONSTRUCTION, with no MIN against S_mf + S_diff
  ! --- rec 169 ON-only diagnostics (counters and ledgers, never fluxes ---
  INTEGER  :: znk_days_active   = 0        ! days the mode-1 zinc limb was on the path
  INTEGER  :: znk_days_qdiff0   = 0        ! days q_diff <= 0, so c_surf fell back to the bulk solution concentration
  INTEGER  :: znk_days_csurf0   = 0        ! days c_surf <= 0, so the active term was 0.0 exactly
  REAL(wp) :: znk_sat_frac_sum  = 0.0_wp   ! sum of c_surf/(km+c_surf) over znk_days_active, for the season-mean saturation
  REAL(wp) :: znk_s_hats_cum    = 0.0_wp   ! mg m-2, cumulative mode-1 zinc active supply before the pool cap

  ! ======================================================================== !
  ! Records 170 and 173 are read and never used; they hold their positions so
  ! that the positional driver file keeps its numbering.
  INTEGER  :: reserved_rec_170 = 0   ! RESERVED (diffusion temperature, not built); read, never used

  ! Record 182 is read and never used; it holds its position in the positional
  ! driver file.
  INTEGER  :: reserved_rec_182 = 0   ! rec 182 RESERVED ( vegetative-N demand basis, NOT BUILT); read, never used

  ! ======================================================================== !
  ! Nitrogen basis of stem-reserve deposition - driver record 172, accepted at
  ! the value 0, so the pre-anthesis reserve bank keeps its constant deposition
  ! fraction. The nitrogen-status predictor the alternative would use is the
  ! nitrogen nutrition index on the Justes et al. (1994) Ann Bot 74:397-407
  ! curve, and the measured reserve-versus-nitrogen relation it is judged
  ! against is the shoot nitrogen and shoot dry matter of Manderscheid, Dier,
  ! Erbs, Sickora & Weigel (2020) Open Data J Agric Res 6,
  ! doi:10.18174/odjar.v6i0.16397.
  INTEGER  :: stem_reserve_deposition_n_basis = 0   ! rec 172; 0 = OFF (EOF-default: the pre-anthesis WSC bank keeps the constant record-56 stem_reserve_frac), 1 = the van Herwaarden 1997 PhD Ch.4 NNI law of
  REAL(wp), PARAMETER :: vhw_a_mgg = 435.26_wp   ! mg WSC g-1 stem DM at NNI = 0; mean of the three site intercepts of van Herwaarden 1997 PhD Ch.4 Tables 4.7/4.9/4.11 (447.79, 389.27, 468.72)
  REAL(wp), PARAMETER :: vhw_b_mgg = -311.26_wp  ! mg WSC g-1 stem DM per NNI unit; mean of the three site slopes of the same tables (-241.35, -281.47, -410.95)
  ! --- rec 172 ON-only diagnostics (ledgers, never fluxes ---
  REAL(wp) :: vhw_nni_last   = 0.0_wp    ! -, NNI used on the last deposition day
  REAL(wp) :: vhw_frac_last  = 0.0_wp    ! g g-1, f_eff on the last deposition day
  REAL(wp) :: vhw_frac_min   = 9.9_wp    ! g g-1, smallest f_eff used
  REAL(wp) :: vhw_frac_max   = -9.9_wp   ! g g-1, largest f_eff used
  INTEGER  :: vhw_days_live  = 0         ! deposition days on which the law was evaluated
  INTEGER  :: vhw_days_clamp_hi = 0      ! days f_eff hit the existing 0.35 record-56 ceiling
  INTEGER  :: vhw_days_clamp_lo = 0      ! days f_eff hit the existing 0.0 record-56 floor

  ! --- Inputs published by main_MMCO2Crop.f90 immediately before the call --- !
  REAL(wp) :: vhw_w_tha_in = 0.0_wp      ! t DM ha-1, today's above-ground DM on the nni_up basis
  REAL(wp) :: vhw_n_sh_in  = 0.0_wp      ! g N m-2, today's shoot N on the nni_up basis

  ! ======================================================================== !
  ! Mineral-nitrogen coupling - driver record 171, accepted at the value 1, the
  ! flux (daily) basis: the grain mineral load is bounded by the day grain
  ! nitrogen import, because the loading of iron and zinc into the grain is
  ! co-regulated with protein. Sources: Waters, Uauy, Dubcovsky & Grusak (2009)
  ! J Exp Bot 60:4263, doi:10.1093/jxb/erp257, pp.4267-4268, where the NAM genes
  ! co-regulate nitrogen, iron and zinc remobilization; Uauy, Brevis &
  ! Dubcovsky (2007) Plant Soil, doi:10.1007/s11104-007-9466-3, p.6; and the
  ! measured grain mineral-protein associations of Peleg et al. (2008), Kutman
  ! et al. (2010), Fan et al. (2008) and Zhao et al. (2009).
  INTEGER  :: mineral_n_coupling_basis = 0  ! rec 171; 0 = OFF (EOF-default: the module reads no N state on this path), 1 = FLUX/DAILY basis, 2 = POOL/CUMULATIVE basis

  ! ======================================================================== !
  ! Coupling of the leaf-water-potential solve to the root sink - driver record
  ! 173, accepted at the value 1: Richards is forced with the converged hourly
  ! transpiration series of the current outer pass, so the stomatal and the soil
  ! solutions see the same water.
  INTEGER  :: spac_sink_coupling = 0   ! rec 173; 0 = OFF (EOF-default), 1 = Richards is forced with the current outer pass's own converged SPAC hourly canopy transpiration
  ! --- rec 173 ON-only diagnostics (ledgers, never fluxes ---
  REAL(wp) :: ssc_relax_scale   = 1.0_wp  ! -, the relaxation factor applied to the hourly series this pass
  REAL(wp) :: ssc_resid_abs_max_mm = 0.0_wp ! mm, largest |E_pass2 - trans_act| over the season
  INTEGER  :: unit_ssc_daily    = -1
  LOGICAL  :: ssc_daily_open    = .FALSE.

  ! Leaf-water-potential solver robustness: extra bisection steps are allowed
  ! after the ordinary iteration count is exhausted, rather than stopping the
  ! run, because a leaf water potential in the measured wheat range (Scholander,
  ! Bradstreet, Hemmingsen & Hammel 1965, Science 148:339-346) is physically
  ! admissible even when the tolerance is slow to close.
  INTEGER,  PARAMETER :: spac_extra_iter = 60          ! extra psi_L steps allowed only after spac_max_iter is exhausted
  REAL(wp), PARAMETER :: spac_tol_fallback_MPa = 1.0e-5_wp  ! bounded acceptance threshold, MPa
  INTEGER  :: spac_fallback_accepts_run = 0            ! times (b) accepted a state; 0 on every deck of the v1.0 configuration

  ! ======================================================================== !
  ! Flag-leaf-number prediction - driver record 174, accepted at the value 0, so
  ! the top rank stays the measured panel median until the anchored booting
  ! stage. The Haun-stage chronology the prediction would use is the phyllochron
  ! of McMaster, Wilhelm, Palic, Porter & Jamieson (2003) Ann Bot 91:697-705,
  ! doi:10.1093/aob/mcg074, Table 2 p.701 (air basis, cultivar-settable), with
  ! the appearance-rate decline and the seed-reserve ratio of Streck, Weiss &
  ! Baenziger (2003) Ann Bot 92:181-190, Eq. 5 p.183, doi:10.1093/aob/mcg121.
  INTEGER  :: canopy_fln_prediction = 0   ! rec 174; 0 = OFF (EOF-default: n_top stays the panel median cba_Nl_fit = 11 until the Booting freeze and then jumps to the emergent clr_FLN, retro-dimensioning every finished blade), 1 = FLN_pred = FLOOR(clr_HS_at_tau(stage_thresh(4) - tdu_at_emerge)) is computed on the emergence day and used as n_top from emergence, the Booting freeze only CONFIRMING it
  ! --- rec 174 state: the prediction, its input and its confirmation ---
  REAL(wp) :: clr_tau_boot_target = -1.0_wp    ! degree-Cd since emergence at which the anchored Booting stage will be set = stage_thresh(4) - tdu_at_emerge; set by main_MMCO2Crop.f90 from the phenology state every day before update_canopy_leaf_ranks (-1 = not yet emerged or not available)
  INTEGER  :: clr_FLN_pred     = 0        ! rec 174 >= 1 only: the final leaf number predicted on the EMERGENCE day from clr_tau_boot_target (0 = not predicted; mode 0 leaves it 0)
  REAL(wp) :: cfp_HS_pred      = 0.0_wp   ! -, the closed-form Haun stage at clr_tau_boot_target that FLN_pred is the floor of (reported, never used as a rate)
  INTEGER  :: cfp_FLN_realized = 0        ! -, the emergent clr_FLN actually read off at the anchored Booting stage, recorded for the confirmation test
  INTEGER  :: cfp_mismatch     = 0                ! -, cfp_FLN_realized - clr_FLN_pred (0 = the prediction was confirmed); a non-zero value is reported and never corrected, because no finished blade is re-dimensioned

  ! ======================================================================== !
  ! Grain-nitrogen partition basis - driver record 175, accepted at the value 4.
  ! The handover of the grain nitrogen state at flowering is made on the
  ! measured grain nitrogen concentration at flowering rather than on a derived
  ! quantity. Sources: Bertheloot, Martre & Andrieu (2008) Funct Plant Biol
  ! 35:781-796, doi:10.1071/FP08064, p.786, verbatim, for the measured grain
  ! nitrogen at flowering; Bertheloot, Cournede & Andrieu (2011) Ann Bot
  ! 108:1085-1096, doi:10.1093/aob/mcr119, for the partition within the culm;
  ! and Busch, Sage & Sage (2018) Nature Plants 4:46-54,
  ! doi:10.1038/s41477-017-0065-x, for the definition of the assimilation-side
  ! nitrogen term the continuity is required of.
  INTEGER :: grain_n_partition_basis = 0   ! rec 175; 0 = OFF (EOF-default), 1 = restore the post-anthesis continuity record 158 declares (the record-154 limb's own grain import + organ synthesis enters ndb_rate_prev), 2 = 1 plus FPB's measured grain-N initial condition Ngtot(0) charged to the chaff pool on the first fill day, 3 = 2 plus the grain's unmet sink as one more claimant in the synthesis pro-rata solve (predicted INERT), 4 = 3 but the Ngtot(0) charge is paid against the grain's OWN dry mass at the measured flowering grain [N] of 2.0 % DM (Bertheloot 2011b Table 1(a): 2.4 mg N on 120 mg DM) instead of as a lump on a grain of zero mass
  REAL(wp) :: gnp_ndb_add_day   = 0.0_wp   ! g N m-2 d-1, mode >= 1: TODAY's added continuity term in ndb_rate_prev
  REAL(wp) :: gnp_ndb_add_cum   = 0.0_wp   ! g N m-2, season sum of the added continuity term
  INTEGER  :: gnp_days_ndb_add  = 0        ! days on which the added continuity term was > 0
  REAL(wp) :: gnp_g0_target     = 0.0_wp   ! g N m-2, mode >= 2: the law's own Ngtot(0) on the first fill day
  REAL(wp) :: gnp_g0_deficit    = 0.0_wp   ! g N m-2, mode >= 2: Ngtot(0) minus the grain's own N on that day
  REAL(wp) :: gnp_g0_take       = 0.0_wp   ! g N m-2, mode >= 2: what the chaff pool actually paid
  REAL(wp) :: gnp_g0_refused    = 0.0_wp   ! g N m-2, mode >= 2: the part the chaff structural floor refused (DOCUMENTED, never silent)
  LOGICAL  :: gnp_g0_done       = .FALSE.  ! mode >= 2: one-off guard, the initial condition is applied once per season

  ! Record 175 mode 4: the grain nitrogen mass at flowering,
  ! n2_grain_anth_conc = 0.020 g N g-1 DM = 2.4 mg N per 120 mg DM, measured by
  ! Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot 108:1097-1109,
  ! doi:10.1093/aob/mcr125, Table 1 section (a) "Grains".
  REAL(wp), PARAMETER :: n2_grain_anth_conc = 0.020_wp  ! g N g-1 DM, the MEASURED grain N concentration at flowering = 2.4 mg N / 120 mg DM; Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109 Table 1(a) "Grains" rows "Total N mass N_grain^tot mg 2.4" and "Total dry mass M_grain^tot mg 120" (DOI 10.1093/aob/mcr125); read ONLY on the rec 175 = 4 path
  REAL(wp) :: gnp_g0_paid_cum  = 0.0_wp   ! g N m-2, mode 4 ONLY: the cumulative Ngtot(0) charge paid by the chaff so far this season (mode <= 3 never writes it)
  REAL(wp) :: gnp_g0_frac      = 0.0_wp   ! dimensionless, mode 4 ONLY, REPORTED: gnp_g0_paid_cum / gnf_ng0, how much of the law's own initial condition the grain's dry mass has earned
  REAL(wp) :: gnp_m3_release_cum = 0.0_wp  ! g N m-2, mode 3: cumulative N the pro-rata solve released to the grain ( predicts 0)
  INTEGER  :: gnp_m3_days        = 0       ! mode 3: days on which that release was > 0 ( predicts 0)

  ! ======================================================================== !
  ! Vegetative nitrogen-demand basis - driver record 176, accepted at the value
  ! 0, so the record-155 feedback numerator is unchanged. The alternative reads
  ! the demand from the critical-curve construction of Justes, Mary, Meynard,
  ! Machet & Thelier-Huche (1994) Ann Bot 74:397-407,
  ! doi:10.1006/anbo.1994.1133, as used by Bertheloot, Cournede & Andrieu (2011)
  ! Ann Bot 108:1085-1096, doi:10.1093/aob/mcr119.
  INTEGER :: vegetative_n_demand_basis = 0  ! rec 176; 0 = OFF (EOF-default), 1 = the record-155 limb-(b) E_N numerator reads the CIRCULATING pool only (N_mobile_veg + N_NO3_shoot), dropping the stem and chaff luxury organ N that record 157 = 3 adds, per Bertheloot 2011a p.1090's own definition of N_mob^c, 2 = 1 plus the Justes 1994 Eqn 11 p.401 critical-curve STOCK deficit added to the pre-anthesis demand (still capped by the unchanged Eqn 17 p.405 maximum-curve stock gap)
  REAL(wp) :: vnd_num_drop_day   = 0.0_wp   ! g N m-2, mode >= 1: the organ luxury N dropped from the E_N numerator TODAY
  REAL(wp) :: vnd_num_drop_cum   = 0.0_wp   ! g N m-2, season sum of the dropped luxury N
  INTEGER  :: vnd_days_m1        = 0        ! days on which the dropped luxury N was > 0
  REAL(wp) :: vnd_crit_def_day   = 0.0_wp   ! g N m-2 d-1, mode >= 2: TODAY's critical-curve stock deficit added to the demand
  REAL(wp) :: vnd_crit_def_cum   = 0.0_wp   ! g N m-2, season sum of that stock deficit (as ASKED, not as served)
  INTEGER  :: vnd_days_crit_def  = 0        ! pre-anthesis days on which the stock deficit was > 0
  REAL(wp) :: vnd_nc_crit_pct_day = 0.0_wp  ! %DM, the Justes critical concentration the deficit was computed at today
  INTEGER  :: vnd_days_cap_bind  = 0        ! mode >= 2: pre-anthesis days on which the Justes MAXIMUM-curve stock cap actually bound the enlarged demand (0 before this record, by the record-155 block's own measurement)

  ! ======================================================================== !
  ! Driving temperature of the canopy clock - driver record 161. The released
  ! configuration sets record 161 = 0, so the canopy clock runs on screen-air
  ! degree-days. Mode 1 is an inactive alternative, see Part 1 Table 6: it runs
  ! the clock on soil temperature while the apex is below ground and on canopy
  ! temperature afterwards.
  ! The apex-temperature practice is stated by Jamieson, Brooking, Porter &
  ! Wilson (1995) Field Crops Res 41(1):35-44,
  ! doi:10.1016/0378-4290(94)00102-I, the primary source relayed by McMaster et
  ! al. (2003); Evers, Vos, Andrieu & Struik (2006) Ann Bot 97:649,
  ! doi:10.1093/aob/mcl020, p.650, report that their thermocouples measured the
  ! apex temperature directly. The terminal-spikelet window the clock drives is
  ! 400 degree-Cd after emergence (McMaster, Wilhelm & Frank 2008 Ann Bot
  ! 102:561-569).
  INTEGER  :: canopy_clock_driver = 0   ! rec 161; 0 = OFF (EOF-default: the canopy clock runs on screen-air degree-days), 1 = the canopy clock runs on Tsoil_prof(1) degree-days from emergence until the model's own terminal-spikelet marker, then on screen air (Abichou 2016 thesis HAL tel-02316154 printed p.23 and p.38; Evers et al. 2006 Ann Bot 97:650; Gauthier et al. 2020 JXB SI Table S1.6 p.28; SHARED.f90)

  ! ======================================================================== !
  ! Basis of the canopy clock increment - driver record 162, accepted at the
  ! value 1, so the canopy clock is advanced by the gated increment
  ! MIN(VF, DF) dTT over the window the source defines, rather than by the raw
  ! degree-day increment. Source: McMaster, Wilhelm & Frank (2008) Ann Bot
  ! 102(4):561-569, doi:10.1093/aob/mcn115, whose printed scope is stage
  ! progression from emergence, with the phyllochron of McMaster, Wilhelm,
  ! Palic, Porter & Jamieson (2003) Ann Bot 91(6):697-705, Table 2 p.701,
  ! doi:10.1093/aob/mcg074.
  INTEGER  :: canopy_clock_basis = 0   ! rec 162; 0 = OFF (EOF-default: the canopy clock is advanced by the raw degree-day increment current_dtt), 1 = the canopy clock is advanced by the VF/DF-gated development increment current_dtdu that the stage engine already uses (McMaster et al. 2008 Ann Bot 102:561-569 Eq. (1) p.563, DOI 10.1093/aob/mcn115; LABELED AN INFERENCE, the paper does not print that leaf appearance is gated; SHARED.f90)

  ! ======================================================================== !
  ! Chaff structural-nitrogen basis (driver record 178) and the per-kernel
  ! grain-nitrogen trait (driver record 177).
  ! Record 178: the released configuration sets it to 0, so the chaff keeps the
  ! STEM structural-nitrogen floor n2_struct_stem. Modes 1 and 2, which give the
  ! chaff its own measured floor, are an inactive alternative, see Part 1
  ! Table 6.
  ! Record 177 is accepted at the value 0, so the per-kernel grain-nitrogen
  ! ceiling of record 125 is the compiled constant rather than a cultivar trait.
  ! The measured basis of both floors is Bertheloot, Martre & Andrieu (2008)
  ! Funct Plant Biol 35:781-796, doi:10.1071/FP08064, Table 1 p.784 and p.786.
  INTEGER  :: chaff_structural_n_basis = 0   ! rec 178; 0 = OFF (EOF-default: the chaff keeps the STEM structural-N floor n2_struct_stem = 3.5e-3 g N g-1 DM), 1 = the chaff uses its OWN measured floor n2_struct_chaff = 5.10e-3 g N g-1 DM (1.07 mg N / 210 mg DM, Bertheloot, Wu, Cournede & Andrieu 2011 Ann Bot 108:1097-1109 Table 1(c) column "Chaff", DOI 10.1093/aob/mcr125) at the ten chaff loci !: 2 = the NRE LAW, N_chaff >= (1 - 0.56) x N_chaff(anthesis), with Pask 2012's measured chaff structural minimum 3.80e-3 g N g-1 DM as the hard floor; 3 = the same law on Bertheloot 2011b's 5.10e-3 floor

  INTEGER  :: grain_n_potential_trait = 0   ! rec 177; 0 = OFF (EOF-default: the record-125 per-kernel grain-N ceiling is the compiled gnf_mgN_kernel_hi = 1.38 mg N kernel-1 for every cultivar), 1 = the OPTIONAL name-keyed row `grainN_pot_mgN_kernel` of input/cultivar_yield.csv overrides it per cultivar, enforced inside the Martre et al. 2003 Table I control range
  ! --- rec 177 trait state (defined so that the OFF path is the SAME bit pattern ---
  REAL(wp) :: gnp_mgN_kernel = gnf_mgN_kernel_hi  ! mg N kernel-1, THE LIVE per-kernel grain-N ceiling the rec-125 path applies. Initialized FROM the compiled constant itself, so with no trait row gnc_cap is computed from the identical value as before this record
  LOGICAL  :: gnp_from_file  = .FALSE.            ! .TRUE. only when a cultivar file supplied the name-keyed row AND rec 177 = 1
  REAL(wp) :: gnp_row_value  = -1.0_wp            ! mg N kernel-1 as read from the file (-1 = no row seen); reported, so the run.log carries what the deck asked for as well as what it got

  ! Enforced calibration window of the per-kernel grain-nitrogen trait: the
  ! control minimum and maximum measured by Martre, Porter, Jamieson & Triboi
  ! (2003) Plant Physiol 133:1959-1967, doi:10.1104/pp.103.030585, Table I
  ! p.1961 (0.75 mg N per kernel for cv. Recital to 1.38 mg N per kernel).
  REAL(wp), PARAMETER :: gnp_range_lo  = 0.75_wp  ! mg N kernel-1, Martre 2003 Table I p.1961 CONTROL minimum (cv. Recital)
  REAL(wp), PARAMETER :: gnp_range_hi  = 1.11_wp  ! mg N kernel-1, Martre 2003 Table I p.1961 CONTROL maximum (cv. Renan)
  REAL(wp), PARAMETER :: gnp_margin    = 0.10_wp    ! -, the enforced tolerance: outside the range by more than this fraction is an ERROR STOP, inside it a run.log warning

  ! ======================================================================== !
  ! Nitrogen limit on lamina expansion - driver record 179, accepted at the
  ! value 1, the primary form: lamina expansion is limited by the running
  ! maximum of the model's own photosynthetic specific leaf nitrogen,
  ! max N^ph,a, which is Bertheloot, Martre & Andrieu (2008) Funct Plant Biol
  ! 35:781-796, doi:10.1071/FP08064, Eq. (8) and Table 1, both p.784. The
  ! structural floor it works against is n2_struct_lam = 0.0079 g N g-1 DM
  ! (Bertheloot, Martre & Andrieu 2008 Plant Physiol 148:1707, p.1711,
  ! doi:10.1104/pp.108.124156). No new basis and no new constant is introduced.
  ! The measured green-area and nitrogen series it is judged against are those
  ! of Manderscheid, Dier, Erbs, Sickora & Weigel (2020) Open Data J Agric Res
  ! 6:28-33, doi:10.18174/odjar.v6i0.16397.
  INTEGER :: lamina_expansion_n_limit = 0   ! rec 179; 0 OFF (EOF-default), 1 PRIMARY (Bertheloot 2008 FPB Eqn 8 minN^ph,a, running max), 2 arm (Eqn 8 plateau 1.114), 3 arm (Jamieson & Semenov 1.5 on total lamina N), 4 attribution arm (Yin l_N log form)
  REAL(wp), PARAMETER :: lnl_a_slope = 0.46_wp   ! dimensionless; slope a of Bertheloot, Martre & Andrieu 2008 Funct Plant Biol 35:781-796 Eqn 8, Table 1 p.784, doi 10.1071/FP08064; one value for cvs Apache/Isengrain/Thesee
  REAL(wp), PARAMETER :: lnl_b_max   = 1.9_wp    ! g N m-2 (Table 1 p.784 prints 1.9e-3 kg m-2); maxN^ph,a above which minN^ph,a stops rising
  REAL(wp), PARAMETER :: lnl_c_icpt  = 0.24_wp   ! g N m-2 (Table 1 p.784 prints 0.24e-3 kg m-2); intercept c of Eqn 8

  ! Green-area guard on the running maximum: below a very small green lamina
  ! area the per-area quantity is not defined, so the guard is an assumed
  ! domain limit and is not a literature constant. The equation it guards is
  ! Eq. (8) p.784 of Bertheloot, Martre & Andrieu (2008) Funct Plant Biol
  ! 35:781-796.
  REAL(wp), PARAMETER :: lnl_gai_max_min = 0.10_wp    ! m2 m-2 green lamina area; an assumed domain guard, not a literature constant: below this area the ratio N^ph/cls_gai is not the specific leaf nitrogen of any lamina and must not enter the running maximum of Eq. (8)
  INTEGER  :: lnl_days_maxupd  = 0        ! days the running maximum maxN^ph,a was actually updated under the guard (mode 1 only)
  REAL(wp) :: lnl_gai_at_maxupd = 0.0_wp  ! m2 m-2; cls_gai on the last day the running maximum was updated
  ! --- rec 179 state and ON-only diagnostics (always defined, so the OFF bit pattern is the SAME ---
  REAL(wp) :: lnl_maxNpha      = 0.0_wp   ! g N m-2 leaf; running maximum of the model's own photosynthetic SLN (Bertheloot Eqn 8 maxN^ph,a), monotone non-decreasing, mode 1 only
  REAL(wp) :: lnl_Nph_day      = 0.0_wp   ! g N m-2 ground; lamina photosynthetic N = npg_N_leaf - n2_struct_lam*W_green
  REAL(wp) :: lnl_slnph_day    = 0.0_wp   ! g N m-2 leaf; today's photosynthetic SLN = lnl_Nph_day / cls_gai
  REAL(wp) :: lnl_minNpha_day  = 0.0_wp   ! g N m-2 leaf; the threshold actually used today (mode-dependent)
  REAL(wp) :: lnl_a_sup_day    = 0.0_wp   ! m2 m-2; the green area the lamina N can hold at that threshold
  REAL(wp) :: lnl_allow_day    = 0.0_wp   ! m2 m-2; MAX(0, A_sup - cls_gai), the room above the standing green area
  REAL(wp) :: lnl_fN_day       = 1.0_wp   ! -; the area-basis paid fraction MIN(1, allow/increment)
  REAL(wp) :: lnl_refused_day  = 0.0_wp   ! m2 m-2; area this record refused TODAY BEYOND what records 121/122 already refused
  REAL(wp) :: lnl_refused_cum  = 0.0_wp   ! m2 m-2; season cumulative of the above
  INTEGER  :: lnl_days_refuse  = 0        ! days on which lnl_refused_day > 0
  LOGICAL  :: lnl_live_today   = .FALSE.  ! .TRUE. only on a day the record evaluated (rec >= 1 and the rec-122 increment ledger was live)

  ! ======================================================================== !
  ! Magnesium solution concentration - driver record 181, accepted at the value
  ! 0, so magnesium takes its solution concentration from the linear sorption
  ! partition. The measured soil-solution window the alternative would use is
  ! Barakat, M. A. et al. (2025) Scientific Reports 15:44670,
  ! doi:10.1038/s41598-025-29041-y, Table 2 p.3.
  INTEGER  :: mg_solution_basis = 0   ! rec 181; 0 = OFF (EOF-default: magnesium takes its solution concentration from the linear Kd partition (pool_mg_kg)/Kd, mineral_module.f90), 1 = magnesium joins the record-131 solution branch exactly as Fe/Zn/Cu/Mn do -- driver_minerals.csv column 17 avail_A_mgL_CAL on the Mg row IS the soil-solution Mg concentration (mg L-1), feeding the unchanged mass-flow and diffusion terms
  ! --- rec 181 enforced calibration window, Barakat 2025 Table 2 p.3, +- 1 SE ---
  REAL(wp), PARAMETER :: mgsb_range_lo = 7.7776_wp      ! mg L-1 = (0.47 - 0.15) mmol kg-1 x 24.305 g mol-1; Barakat et al. 2025 Sci Rep 15:44670, Table 2 p.3 ("Soluble Mg mmol kg-1 0.47 +- 0.15"), DOI 10.1038/s41598-025-29041-y
  REAL(wp), PARAMETER :: mgsb_range_hi = 15.0691_wp  ! mg L-1 = (0.47 + 0.15) mmol kg-1 x 24.305 g mol-1; same table row, same first-hand read
  REAL(wp), PARAMETER :: mgsb_margin   = 0.10_wp          ! -, the enforced tolerance, the same rule and the same 0.10 that record 177 uses: outside the range by more than this fraction is an ERROR STOP, inside it a run.log warning

  ! --- Ledgers only; never a flux and never read back into the physics --- !
  REAL(wp) :: mgsb_soln_used   = 0.0_wp   ! mg L-1, the concentration mode 1 actually applied (0.0 when OFF)
  LOGICAL  :: mgsb_warn_done   = .FALSE.  ! has the margin WARNING already been printed in this run?

  ! ======================================================================== !
  ! Buffer power of the mineral diffusion coefficient - driver record 183,
  ! accepted at the value 1, so the diffusion coefficient carries the measured
  ! buffer power theta + rho_b Kd rather than the water content alone.
  ! Sources: Baldwin, Nye & Tinker (1973) "Uptake of solutes by multiple root
  ! systems from soil. III. A model for calculating the solute uptake by a
  ! randomly dispersed root system developing in a finite volume of soil",
  ! Plant and Soil 38:621-635, doi:10.1007/BF00010701; Roose & Fowler (2004)
  ! J Theor Biol, doi:10.1016/j.jtbi.2003.10.008; and the wheat-soil impedance
  ! measurements of Barraclough & Tinker (1981) J Soil Sci 32:225-236,
  ! doi:10.1111/j.1365-2389.1981.tb01702.x.
  INTEGER  :: mineral_buffer_power_basis = 0   ! rec 183; 0 = OFF (EOF-default: the mineral diffusion coefficient De_e = Dl_e*theta*f_l/b_e divides by the driver_minerals.csv column-19 CONSTANT buffer_power_b, 1.0 on every shipped row and no published source located by record 140), 1 = ON: b_e = theta + rho_b*Kd_e is recomputed daily per element from the model's own soil water (vwc_root), the deck's own MEASURED topsoil bulk density (nud_bulk_density, driver_N.csv position 70) and the element's own driver_minerals.csv column-3 Kd -- the linear-isotherm buffer power of Baldwin, Nye & Tinker 1973 Plant Soil 38:621 p.622 ("b = buffer power of soil = dC/dCl") and p.625 ("Db = Dl.W.fl"), decomposed per Roose & Fowler 2004 J Theor Biol 228:173 Eq. (2.8) p.175
  REAL(wp), PARAMETER :: bpb_warn_b_lo = 1.0_wp  ! -, a computed b_e below this would mean record 183 INCREASES diffusion (the opposite of the mechanism); run.log WARNING,. The value 1.0 is the NON-ADSORBED limit b = theta + 0 ~ 1 of So & Nye 1989 p.744 Eq.(3), not a tuned threshold

  ! --- Ledgers, dimensioned by the fixed element order of mineral_module --- !
  REAL(wp) :: bpb_rho_b_used  = -1.0_wp        ! kg L-1, the bulk density mode 1 actually used (-1 = OFF / never set)
  LOGICAL  :: bpb_rho_fallback = .FALSE.       ! .TRUE. if the deck supplied no driver_N.csv position 70 and the module's own 1.63 kg L-1 basis was used
  REAL(wp) :: bpb_kd_used(6)  = 0.0_wp         ! L kg-1, the Kd each element's b_e was built from (Fe's is the derived operational 6.7/0.119)
  REAL(wp) :: bpb_b_sum(6)    = 0.0_wp         ! -, sum of b_e over the days the diffusion locus ran (season mean = /bpb_b_days)
  REAL(wp) :: bpb_b_min(6)    =  1.0e30_wp     ! -, season minimum of b_e
  REAL(wp) :: bpb_b_max(6)    = -1.0e30_wp     ! -, season maximum of b_e
  INTEGER  :: bpb_b_days(6)   = 0              ! days the diffusion locus computed a b_e for this element
  INTEGER  :: bpb_warn_lo(6)  = 0              ! days b_e < bpb_warn_b_lo (WARNING counter)

  ! ======================================================================== !
  ! Grain calcium loading - driver record 180, accepted at the value 1, the
  ! xylem-stream law: calcium reaches the grain in the xylem water and is
  ! divided at the ear by the measured partition, instead of being the product
  ! of a grain water import and a calibrated xylem concentration.
  ! Sources: Hocking (1994) J Plant Nutr 17:1289-1308,
  ! doi:10.1080/01904169409364807, Tables 1 and 4, for the head and grain
  ! shares of post-anthesis calcium arrival; White & Broadley (2003) Ann Bot
  ! 92:487-511, doi:10.1093/aob/mcg164, p.489 ("The delivery of Ca to" the
  ! shoot is in the xylem) and p.490 for the reported xylem calcium range;
  ! Atkinson, Else, Taylor & Dover (2011) J Exp Bot 62:2233-2250,
  ! doi:10.1093/jxb/err111; the ear water flux of Lizana et al. (2010) J Exp Bot
  ! 61:1147, doi:10.1093/jxb/erp380, Tables 2 and 3; and Barlow, Lee, Munns &
  ! Smart (1970) Aust J Biol Sci 23:783-791, doi:10.1071/BI9700783, Summary
  ! p.783.
  ! The elevated-CO2 mineral-decline band of Loladze (2014) eLife 3:e02245, p.6
  ! Fig. 2, is reported as a falsification target and applied to nothing.
  INTEGER  :: grain_ca_loading_basis = 0   ! rec 180; 0 = OFF (EOF-default: the record-166 calcium route stands), 1 = the xylem-stream loading law of (Ca = c_xyl_Ca x Tr_act, partitioned from anthesis by the Hocking 1994 measured fractions)

  ! --- The two measured ear partition constants (Hocking 1994) --- !
  REAL(wp), PARAMETER :: gcl_f_head_post  = 0.5827_wp    ! -, the HEAD's share of post-anthesis Ca arrival; Hocking 1994 J Plant Nutr 17:1289 Tables 1(b)+1(c) p.1294 with Table 2 p.1295
  REAL(wp), PARAMETER :: gcl_f_grain_post = 0.30549_wp   ! -, the GRAIN's share of post-anthesis Ca arrival = gcl_f_head_post x 0.52427; Hocking 1994 Table 1(c) p.1294 (grain 10.8 / (10.8 + 9.8))

  ! --- The enforced xylem calcium envelope (White & Broadley 2003) --- !
  REAL(wp), PARAMETER :: gcl_cxyl_lo     = 12.0_wp   ! mg L-1 (= 0.3 mmol L-1), White & Broadley 2003 Ann Bot 92:487 p.490 lower reported xylem [Ca]; Atkinson 1992 Table 3 p.1319 low entry 0.4 mol m-3
  REAL(wp), PARAMETER :: gcl_cxyl_hi     = 661.0_wp  ! mg L-1 (= 16.5 mmol L-1), White & Broadley 2003 p.490 upper reported xylem [Ca]
  REAL(wp), PARAMETER :: gcl_cxyl_ref    = 15.0_wp   ! mg L-1, the SOURCED DEFAULT written into every deck's driver_minerals.csv Ca row; reported, never substituted for a user value
  REAL(wp), PARAMETER :: gcl_cxyl_margin = 0.10_wp      ! -, the enforced tolerance, the same rule record 177 applies: outside the envelope by more than this fraction is an ERROR STOP, inside it a run.log warning

  ! --- The pre-registered elevated-CO2 band; reported, never applied --- !
  REAL(wp), PARAMETER :: gcl_eco2_lo_pct = -10.0_wp  ! %, Loladze 2014 eLife 3:e02245 p.6 Fig. 2 lower end of the eCO2 mineral decline
  REAL(wp), PARAMETER :: gcl_eco2_hi_pct = -6.5_wp   ! %, Loladze 2014 eLife 3:e02245 p.6 Fig. 2 upper end
  ! --- rec 180 ON-only diagnostics (ledgers, never fluxes --- !
  REAL(wp) :: gcl_ca_grain_cum = 0.0_wp   ! mg m-2, cumulative Ca delivered to the GRAIN terminus by this law
  REAL(wp) :: gcl_ca_chaff_cum = 0.0_wp   ! mg m-2, cumulative Ca delivered to the CHAFF terminus
  REAL(wp) :: gcl_ca_veg_cum   = 0.0_wp   ! mg m-2, cumulative Ca delivered to the VEGETATIVE terminus
  REAL(wp) :: gcl_tr_season_cum= 0.0_wp   ! L m-2, cumulative Tr_act seen by this law over the whole season
  REAL(wp) :: gcl_tr_post_cum  = 0.0_wp   ! L m-2, cumulative Tr_act from anthesis onwards
  REAL(wp) :: gcl_cxyl_used    = 0.0_wp   ! mg L-1, the c_xyl_Ca actually in force (echoed, so the run.log carries what the deck supplied)
  INTEGER  :: gcl_days_post    = 0        ! days the post-anthesis partition ran
  INTEGER  :: gcl_days_pool    = 0        ! days the labile Ca pool bounded the law (raw > pool)

  ! ======================================================================== !
CONTAINS

  SUBROUTINE sync_hydro_step_state()
    hydro_step%failed = hydro_failed
    hydro_step%min_dt_reached = hydro_min_dt_reached
    hydro_step%dt_halvings = hydro_dt_halvings
    hydro_step%cap_active = hydro_cap_active
    hydro_step%last_dt = hydro_last_dt
    hydro_step%mb_resid_m = hydro_mb_resid_m
    hydro_step%mb_correction_m = hydro_mb_correction_m
    hydro_step%top_flux_m = hydro_top_flux_m
    hydro_step%top_infil_m = hydro_top_infil_m
    hydro_step%bottom_flux_m = hydro_bottom_flux_m
    hydro_step%root_sink_m = hydro_root_sink_m
    hydro_step%flux_top_actual_m_d = flux_top_actual
    hydro_step%flux_top_req_m_d = hydro_flux_top_req_m_d
    hydro_step%flux_top_cap_m_d = hydro_flux_top_cap_m_d
    hydro_step%flux_top_applied_m_d = hydro_flux_top_applied_m_d
    hydro_step%infil_cap_m_d = hydro_infil_cap_m_d
    hydro_step%evap_cap_m_d = hydro_evap_cap_m_d
    hydro_step%evap_rate_pot_m_d = hydro_evap_rate_pot_m_d
    hydro_step%evap_rate_day_m_d = hydro_evap_rate_day_m_d
    hydro_step%sum_avail_m = hydro_sum_avail_m
    hydro_step%k1_m_d = hydro_k1_m_d
    hydro_step%h1_m = hydro_h1_m
    hydro_step%theta1 = hydro_theta1
    hydro_step%h_top_m = hydro_h_top_m
    hydro_step%h_crit_atm_m = hydro_h_crit_atm_m
    hydro_step%q34_m_d = hydro_q34_m_d
    hydro_step%k34_m_d = hydro_k34_m_d
    hydro_step%supply_step_mm = hydro_supply_step_mm
    hydro_step%infil_step_mm = hydro_infil_step_mm
    hydro_step%infil_from_supply_step_mm = hydro_infil_from_supply_step_mm
    hydro_step%infil_from_pond_step_mm = hydro_infil_from_pond_step_mm
    hydro_step%soil_evap_step_mm = hydro_soil_evap_step_mm
    hydro_step%pond_evap_step_mm = hydro_pond_evap_step_mm
    hydro_step%runoff_step_mm = hydro_runoff_step_mm
    hydro_step%ponding_start_step_mm = hydro_ponding_start_step_mm
    hydro_step%ponding_end_step_mm = hydro_ponding_end_step_mm
  END SUBROUTINE sync_hydro_step_state

  ! Overlap-weighted binning of one hydrology-call mass into the hourly substep
  ! windows, so that the daily total is conserved to round-off.
  SUBROUTINE bin_series_add(bins, t0, t1, mass_mm)
    REAL(wp), INTENT(INOUT) :: bins(:)
    REAL(wp), INTENT(IN)    :: t0, t1, mass_mm
    INTEGER  :: h, nb
    REAL(wp) :: bl, bu, ov, wlen
    nb = MIN(SIZE(bins), MAX(1, eb_nsteps_s))
    wlen = MAX(1.0e-12_wp, t1 - t0)
    DO h = 1, nb
      bl = REAL(h-1, wp) / REAL(nb, wp)
      bu = REAL(h,   wp) / REAL(nb, wp)
      ov = MAX(0.0_wp, MIN(t1, bu) - MAX(t0, bl))
      IF (ov > 0.0_wp) bins(h) = bins(h) + mass_mm * (ov / wlen)
    END DO
  END SUBROUTINE bin_series_add

  ! Overlap-weighted extraction of the demand mass lying in one hydrology-call
  ! window, the inverse of the binning above.
  FUNCTION bin_series_get(series, t0, t1) RESULT(mass_mm)
    REAL(wp), INTENT(IN) :: series(:)
    REAL(wp), INTENT(IN) :: t0, t1
    REAL(wp) :: mass_mm
    INTEGER  :: h, nb
    REAL(wp) :: bl, bu, ov
    nb = MIN(SIZE(series), MAX(1, eb_nsteps_s))
    mass_mm = 0.0_wp
    DO h = 1, nb
      bl = REAL(h-1, wp) / REAL(nb, wp)
      bu = REAL(h,   wp) / REAL(nb, wp)
      ov = MAX(0.0_wp, MIN(t1, bu) - MAX(t0, bl))
      IF (ov > 0.0_wp) mass_mm = mass_mm + series(h) * ov / (bu - bl)
    END DO
  END FUNCTION bin_series_get

  FUNCTION dHour_to_dDay(var_hour, constant_Gauss)
    REAL(wp), INTENT(IN) :: var_hour, constant_Gauss
    REAL(wp) :: dHour_to_dDay
    dHour_to_dDay = var_hour * 3600.0_wp * constant_Gauss * 24.0_wp * 1.0E-6_wp
  END FUNCTION dHour_to_dDay

  ! driver_main.txt records carry a trailing "!" comment column. A path is read
  ! as the whole remaining field, so the comment must be stripped first.
  SUBROUTINE strip_driver_comment(str)
    CHARACTER(LEN=*), INTENT(INOUT) :: str
    INTEGER :: ip
    ip = INDEX(str, '!')
    IF (ip > 0) str = str(1:ip-1)
    str = ADJUSTL(str)
    ip = LEN_TRIM(str)
    IF (ip > 0) str = str(1:ip)
  END SUBROUTINE strip_driver_comment

  ! Event-reachability refusal: a management event whose date falls outside the
  ! simulated period would be read and silently never applied, so it is refused
  ! instead. Safety and reproducibility only; no model law depends on it.
  SUBROUTINE check_event_reachability()
    INTEGER :: i_chk, d_chk
    IF (iD_sowing < 1 .OR. iD_sowing > NDays) THEN
      WRITE(*,'(A,I0,A,I0,A)') ' ERROR [R88-BOUND]: rec 88 sowing_model_day = ', iD_sowing, &
        ' is outside 1..', NDays, ' (days of weather in this deck). REFUSING TO RUN.'
      STOP 3
    END IF
    IF (ALLOCATED(N_fert_DAS_arr)) THEN
      DO i_chk = 1, N_fert_n_apps
        d_chk = N_fert_DAS_arr(i_chk) + iD_sowing
        IF (d_chk < 1 .OR. d_chk > NDays) THEN
          WRITE(*,'(A,I0,A,I0,A,I0,A)') ' ERROR [R88-BOUND]: fertiliser DAS ', &
            N_fert_DAS_arr(i_chk), ' would fire on model day ', d_chk, ', outside 1..', NDays, &
            '. The event would be SILENTLY SKIPPED. REFUSING TO RUN.'
          STOP 3
        END IF
      END DO
    END IF
    IF (ALLOCATED(N_org_das) .AND. ALLOCATED(N_org_age)) THEN
      DO i_chk = 1, N_org_n
        IF (N_org_age(i_chk) > 0) CYCLE
        d_chk = N_org_das(i_chk) + iD_sowing
        IF (d_chk < 1 .OR. d_chk > NDays) THEN
          WRITE(*,'(A,I0,A,I0,A,I0,A)') ' ERROR [R88-BOUND]: residue DAS ', N_org_das(i_chk), &
            ' would fire on model day ', d_chk, ', outside 1..', NDays, &
            '. The event would be SILENTLY SKIPPED. REFUSING TO RUN.'
          STOP 3
        END IF
      END DO
    END IF
  END SUBROUTINE check_event_reachability

END MODULE SHARED
  ! ======================================================================== !