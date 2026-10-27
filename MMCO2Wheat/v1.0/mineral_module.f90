!===============================================================================!
! M&M-CO2Wheat v1.0  |  mineral_module.f90
!
! Purpose: Grain minerals, six elements (Fe, Zn, Ca, Cu, Mg, Mn). One DAILY step
!          per element: soil-solution concentration, uptake by mass flow,
!          diffusion to the root surface and element-specific active terms,
!          allocation to the growing organs, remobilization from senescing tissue
!          to the grain, and a per-element mass balance. Spatial unit: one labile
!          soil pool over the whole 0-0.6 m profile plus the organ pools of the
!          crop, per element; every pool in mg element m-2 of ground.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: TYPES (working precision); SHARED, read-only, for the senescence dry-matter
!       fluxes (W_dead, W_green, W_stem, W_chaff and the stem/chaff death
!       ledgers), the living-root constants of the nitrogen-uptake law
!       (aroot_floor, aroot_decline_d), the specific root length nud_SRL, the
!       soil-wetness impedance constants (nud_csif_slope, nud_theta_th), the
!       measured topsoil bulk density nud_bulk_density, the grain nitrogen pools
!       (N_storage, N_chaff), the standing green areas of the organ-level canopy,
!       the measured calcium partition constants and the process switches of the
!       mineral records.
! Output: the labile soil pool M_soil_labile and the organ pools M_root, M_leaf,
!       M_stem, M_chaff, M_grain and M_dead; the daily diagnostics M_uptake_mf,
!       M_uptake_diff, M_uptake_ac, M_uptake_hats, M_remob_grain, M_soln,
!       M_cap_binding and M_balance_resid; and the season ledgers of the mineral
!       records. All are consumed by main_MMCO2Crop.f90, which converts the grain
!       pool to a concentration (mg kg-1) for output and writes the daily mineral
!       files. No carbon, nitrogen, water or energy state is written here.
! Called by: main_MMCO2Crop.f90 - read_mineral_params and init_mineral_pools once
!       per run, update_minerals_daily once per day, after grain filling and
!       before the senescence update.
!===============================================================================!
! ---------------------------------------------------------------------------- !
! THE LAW, PER ELEMENT AND PER DAY (all pools mg element m-2 of ground)
!  1. Labile soil pool over 0-0.6 m: M_soil = pool_mg_kg x 978 kg m-2
!     (bulk density 1.63 g cm-3 x 0.6 m).
!  2. Soil-solution concentration c: a site-level availability parameter for Fe,
!     Zn, Cu and Mn; a linear sorption partition c = (pool / soil mass) / Kd for
!     Mg; a solubility ceiling for Fe where no availability parameter is given;
!     and the xylem-sap concentration for Ca. Zero once the pool is exhausted.
!  3. Uptake = mass flow (transpiration x c) + diffusion to the root surface
!     (Barber-Cushman, with the soil-wetness impedance and the linear-isotherm
!     buffer power) + the saturating zinc transport term + the iron
!     phytosiderophore term, bounded by the labile pool, which is the only cap.
!  4. The day's uptake is allocated to the growing organs by their share of the
!     day's dry-matter increment. Calcium instead follows the transpiration
!     stream and is divided at the ear by the measured grain/chaff partition.
!  5. Grain loading from senescing tissue by the element's phloem export
!     fraction, bounded by the grain protein sink.
!  6. Mass balance: (organs + dead)_after - (organs + dead)_before - uptake = 0,
!     because only uptake adds atoms; senescence and grain loading are internal
!     transfers.
! Elevated CO2 has no term of its own here: it acts through transpiration, root
! mass and activity, the organ dry-matter increments and the senescence flux.
! Mineral shortage does not gate carbon or nitrogen growth.
! Atomic masses are the IUPAC/CIAAW standard atomic weights (Prohaska et al.
! 2022, Pure Appl Chem 94:573, Table 1, DOI 10.1515/pac-2019-0603).
! The senescence increment seen here is the PREVIOUS day's, because this routine
! runs before the senescence update in the daily order: a one-day lag, exact in
! the season total.
! Every element constant is read from input/driver_minerals.csv, whose header
! carries the source, page and DOI of each value.
! ---------------------------------------------------------------------------- !
MODULE MINERAL_MODULE
  USE TYPES, ONLY: wp
  ! SHARED is read-only here: this module writes no carbon, nitrogen, water or
  ! energy state. The senescence fluxes, the living-root constants of the
  ! nitrogen-uptake law, the soil transport constants and the grain nitrogen
  ! pools are imported as runtime values so that two input files cannot diverge.
  USE SHARED, ONLY: W_dead_sh => W_dead, W_green_sh => W_green, &
                    fw_dw_root, aroot_floor, aroot_decline_d, &
                    grainN_partition, W_stem_sh => W_stem, &
                    W_chaff_sh => W_chaff, &
                    n2_dW_dead_stem_day, n2_dW_dead_chaff_day, &
                    mineral_share_basis, &   ! record 114: allocation share basis (re-based by the caller; SHARED.f90)
                    mineral_grain_sink_cap, &   ! record 126: maximum grain concentration bound
                    mineral_supply_law, &   ! record 131: site-level solution availability
                    mineral_supply_mech, &   ! record 140: mass flow + diffusion + saturating zinc transport
                    grain_ca_xylem_throughput, &   ! record 166: ear xylem-throughput calcium route
                    grain_ca_loading_basis, &   ! record 180: post-anthesis calcium partition at the ear
                    gcl_f_head_post, gcl_f_grain_post, &   ! head and grain shares of post-anthesis Ca arrival; Hocking 1994 Tables 1(b)/1(c) p.1294 and Tables 2/3 p.1295
                    gcl_ca_grain_cum, gcl_ca_chaff_cum, gcl_ca_veg_cum, &   ! record 180 season ledgers
                    gcl_tr_season_cum, gcl_tr_post_cum, &
                    gcl_days_post, gcl_days_pool, &
                    mineral_n_coupling_basis, &    ! record 171: grain protein sink on mineral loading
                    ! The only nitrogen state this module reads, and it is read
                    ! read-only: the grain nitrogen pool is MAX(0, N_storage -
                    ! N_chaff), the same expression main_MMCO2Crop.f90 uses for
                    ! the grain nitrogen concentration and the grain protein, so
                    ! no new definition of grain nitrogen is introduced.
                    N_storage_sh => N_storage, N_chaff_sh => N_chaff, &       ! read-only
                    canopy_sheath_ear, &   ! record 166: the ear-area state's own switch
                    cls_gai, cse_sai, cse_iai, cse_eai_meas, &   ! standing green areas of the organ-level canopy, record 166
                    gcx_v_ear_cum, gcx_v_net_cum, gcx_ca_grain_cum, &        ! record 166 season ledgers
                    gcx_ca_chaff_cum, gcx_f_ear_max, gcx_days_route, &
                    nud_SRL, &   ! specific root length of the nitrogen-uptake law, reused unchanged: Chen et al. 2020, BMC Plant Biol 20:185, DOI 10.1186/s12870-020-02390-8
                    mineral_impedance_basis, &   ! record 168: soil-wetness impedance of the diffusion coefficient
                    nud_csif_slope, nud_theta_th, &   ! Olesen et al. 2001, Soil Sci 166:89-96, p.94 Eq. (7a) slope and p.93 Eq. (5) threshold, reused from the nitrogen-uptake law
                    mib_fl_sum, mib_fl_off_sum, mib_fl_days, &               ! record 168 season ledgers
                    mib_days_at_zero, mib_vwc_min, mib_vwc_max, &
                    zn_uptake_kinetic_basis, &   ! record 169: saturating root-surface zinc transport
                    mg_solution_basis, &   ! record 181: magnesium solution concentration; off in v1.0, magnesium keeping the sorption partition below
                    mineral_buffer_power_basis, &   ! record 183: linear-isotherm buffer power of the diffusion coefficient
                    nud_bulk_density, &   ! the deck's own measured topsoil bulk density, reused from the nitrogen-uptake law
                    bpb_rho_b_used, bpb_rho_fallback, bpb_kd_used, &   ! record 183 season ledgers
                    bpb_b_sum, bpb_b_min, bpb_b_max, bpb_b_days, &   ! record 183 season ledgers
                    bpb_warn_lo, bpb_warn_b_lo, &   ! record 183 warning counter and its non-adsorbed-limit threshold
                    znk_days_active, znk_days_qdiff0, znk_days_csurf0, &   ! record 169 season counters
                    znk_sat_frac_sum, znk_s_hats_cum, &   ! record 169 season ledgers
                    iD, &   ! model day index
                    mnp_Ptot_day, mnp_f_root_day, mnp_f_leaf_day, mnp_f_stem_day, &
                    mnp_f_chaff_day, mnp_f_grain_day, mnp_sumf_day, &
                    mnp_ftot_resid_max, mnp_legacy_dWt_day, &
                    mnp_days_share, mnp_days_fallback, mnp_days_vi
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: NMIN = 6
  CHARACTER(LEN=2), PARAMETER :: mineral_names(NMIN) = ['Fe','Zn','Ca','Cu','Mg','Mn']
  ! Fixed element order, named for the branch points below: calcium and
  ! magnesium are excluded from the Fe/Zn/Cu/Mn availability law, calcium taking
  ! its own xylem concentration instead.
  INTEGER, PARAMETER :: idx_Fe = 1, idx_Zn = 2, idx_Ca = 3, idx_Cu = 4, idx_Mg = 5, idx_Mn = 6

  ! Grain and chaff shares of the whole-plant calcium at maturity in an irrigated
  ! spring-wheat field crop: Hocking PJ (1994) J Plant Nutr 17:1289-1308,
  ! DOI 10.1080/01904169409364807, Table 1(c) p.1294 (roots 4.2, stems 27.1,
  ! leaves 48.1, chaff 9.8, grain 10.8 %). Used only as the ratio
  ! grain / (grain + chaff) = 0.524, the grain's share of the calcium the ear
  ! receives on the record-166 route. The same paper's Table 6 p.1297 (grain Ca
  ! 0.52, chaff Ca 1.51 mg g-1) with its Table 1(c) dry-matter shares gives 0.520
  ! independently, the two tables agreeing to 0.4 %. Nothing is fitted here.
  REAL(wp), PARAMETER :: gcx_hk_grain = 10.8_wp   ! % of whole-plant Ca in grain at maturity, Hocking 1994 Table 1(c) p.1294
  REAL(wp), PARAMETER :: gcx_hk_chaff =  9.8_wp   ! % of whole-plant Ca in chaff at maturity, same table

  ! Soil mass of the 0-0.6 m labile-pool layer (kg m-2): measured bulk density
  ! 1.63 g cm-3 x 0.6 m.
  REAL(wp), PARAMETER :: soil_mass_kg_m2 = 978.0_wp

  ! Standard atomic weights, order Fe,Zn,Ca,Cu,Mg,Mn (g mol-1 = ng nmol-1):
  ! IUPAC/CIAAW, Prohaska et al. 2022, Pure Appl Chem 94:573, Table 1,
  ! DOI 10.1515/pac-2019-0603.
  REAL(wp), PARAMETER :: atomic_mass(NMIN) = &
    [55.845_wp, 65.38_wp, 40.078_wp, 63.546_wp, 24.305_wp, 54.938_wp]

  ! Element parameters, read from input/driver_minerals.csv, whose header carries
  ! the source, page and DOI of every value.
  INTEGER  :: mm_modelled(NMIN)  = 0        ! 1 = modelled, 0 = sentinel output
  REAL(wp) :: mm_pool_mgkg(NMIN) = 0.0_wp   ! DTPA/exchangeable pool (mg kg-1 soil)
  REAL(wp) :: mm_kd(NMIN)        = 0.0_wp   ! L kg-1 (0 => ceiling mode, Fe)
  REAL(wp) :: mm_ceiling(NMIN)   = 0.0_wp   ! mg L-1 solubility ceiling on the solution concentration (Fe, Ca)
  REAL(wp) :: mm_vmax_hats(NMIN) = 0.0_wp   ! nmol gFW-1 h-1 (Zn)
  REAL(wp) :: mm_km_hats(NMIN)   = 1.0_wp   ! mg L-1
  REAL(wp) :: mm_vmax_lats(NMIN) = 0.0_wp   ! nmol gFW-1 h-1 (Zn)
  REAL(wp) :: mm_km_lats(NMIN)   = 1.0_wp   ! mg L-1
  REAL(wp) :: mm_ps_rate(NMIN)   = 0.0_wp   ! nmol gDW-1 d-1 (Fe, sufficient, carbonate-corrected)
  REAL(wp) :: mm_ps_rate_def(NMIN) = 0.0_wp ! nmol gDW-1 d-1 (Fe, DEFICIENT sensitivity rate)
  REAL(wp) :: mm_fe_def_switch(NMIN) = 0.0_wp  ! 1 = use the Fe-deficient rate instead of the Fe-sufficient one
  REAL(wp) :: mm_export(NMIN)    = 0.0_wp   ! senescence export fraction to grain (0-1)
  REAL(wp) :: mm_xylem(NMIN)     = 0.0_wp   ! 1 = xylem route instead of phloem allocation (Ca)
  REAL(wp) :: mm_gwr(NMIN)       = 0.0_wp   ! grain water import, g H2O per g final grain DM (Lizana 2010; Ca)
  ! -------------------------------------------------------------------- !
  ! Field-anchored active uptake per g of LIVING root dry matter, applied in
  ! update_minerals_daily as
  !   S_ac = S_ac + mm_field_rate(i) * MAX(0, W_root) * a_root   [mg m-2 d-1]
  ! Non-zero for zinc only; the initial 0.0 is a structural "off" default. The
  ! zinc value, 2.36e-3 mg Zn gDW-1 d-1, is the whole-season field uptake per
  ! unit of root dry matter of Liu, Liu, Zhang, Chen & Zou 2019, Front Plant Sci
  ! 10:426, DOI 10.3389/fpls.2019.00426, Tables 3, 4 and 6, and is pinned in the
  ! driver zinc row. Under record 169 it is the capacity of the saturating
  ! root-surface zinc term below.
  ! -------------------------------------------------------------------- !
  REAL(wp) :: mm_field_rate(NMIN) = 0.0_wp  ! mg gDW-1 d-1 field-anchored active uptake on living root
  REAL(wp) :: mm_grain_max(NMIN) = 0.0_wp   ! record 126: maximum grain concentration, mg kg-1 (0 = no bound); driver column 16
  ! Record 171: driver column 20, the per-element grain element:nitrogen
  ! stoichiometry r_e in mg element per g grain nitrogen (Fe 1.26, Zn 1.84;
  ! Peleg et al. 2008), 0 for the other four. 0.0 = NO BOUND, the same "declared
  ! absent" convention as record 126's column above: no first-hand wheat-grain
  ! source on disk prints this stoichiometry for them, and the campaign's own
  ! measured grain is a falsification target, not an input.
  REAL(wp) :: mm_grain_eN(NMIN) = 0.0_wp    ! record 171: mg element per g grain N (0 = no bound); driver column 20
  REAL(wp) :: mnc_refused(NMIN) = 0.0_wp    ! record 171 ledger: cumulative grain loading refused by the protein sink, mg m-2
  INTEGER  :: mnc_days(NMIN)    = 0         ! record 171 ledger: days the protein sink cut the loading
  REAL(wp) :: mnc_to_grain(NMIN) = 0.0_wp   ! record 171 ledger: cumulative grain loading admitted on days the sink was evaluated, mg m-2
  REAL(wp) :: M_cap_refused(NMIN) = 0.0_wp  ! record 126 ledger: cumulative grain loading refused by the bound, mg m-2
  INTEGER  :: M_cap_days(NMIN) = 0          ! record 126 ledger: days on which the bound cut the loading
  ! Record 131: driver columns 17 and 18, the site-level availability parameters,
  ! calibrated on one control treatment per site; 0.0 elsewhere and inert.
  REAL(wp) :: mm_avail_A(NMIN)       = 0.0_wp   ! mg L-1 solution availability A_e, Fe/Zn/Cu/Mn; driver column 17
  REAL(wp) :: mm_xylem_ca_cal(NMIN)  = 0.0_wp   ! mg L-1 xylem-sap Ca concentration; driver column 18
  ! Record 140: driver column 19, the soil buffer power b_e, a site-level
  ! parameter. Default 1.0; a non-positive or absent value falls back to it,
  ! enforced where the column is read in read_mineral_params below.
  REAL(wp) :: mm_b(NMIN)             = 1.0_wp   ! -, buffer power; column 19 buffer_power_b

  ! Free-solution diffusion coefficients of the dissolved ions at about 25 C and
  ! infinite dilution, order Fe,Zn,Ca,Cu,Mg,Mn (m2 s-1). These are DECLARED
  ! values: no primary measurement was found for them, and divalent
  ! transition-metal and alkaline-earth cations in dilute aqueous solution are
  ! tabulated in the 0.69-0.72e-9 m2 s-1 band. They are a physical property of
  ! the ion rather than a site parameter, so they stay PARAMETER here, unlike the
  ! buffer power above. Calcium's slot is unused (0.0): the diffusion term below
  ! excludes calcium, whose route is the transpiration stream.
  REAL(wp), PARAMETER :: mm_Dl_e(NMIN) = &
    [0.72e-9_wp, 0.70e-9_wp, 0.0_wp, 0.71e-9_wp, 0.71e-9_wp, 0.69e-9_wp]
  ! Depth of the labile pool (m): the same 0-0.6 m basis as soil_mass_kg_m2
  ! above. It stands in for the per-layer thickness of the nitrogen-uptake law,
  ! this module carrying no soil layers of its own, and is not a new soil
  ! measurement.
  REAL(wp), PARAMETER :: mm_depth_m = 0.6_wp
  ! Zinc high-affinity transport kinetics: Hacisalihoglu, Hart & Kochian 2001,
  ! Plant Physiol 125:456, DOI 10.1104/pp.125.1.456, Table II p.459, the -Zn
  ! BDME-10 row, which gives the highest high-affinity capacity in the table.
  ! The measurement is on SEEDLINGS and is never used as a field rate: only the
  ! apparent Km enters, as the saturation constant of the zinc term below, whose
  ! capacity comes from the field rate declared above. The seedling Vmax is
  ! therefore declared here and not used. Km is converted nM -> mg L-1 at the
  ! point of use (Km_nM x atomic_mass(Zn) x 1e-6), computed from the atomic-mass
  ! table above so the two cannot drift apart.
  REAL(wp), PARAMETER :: msm_zn_vmax_hats_nmol_gFW_h = 30.9_wp
  REAL(wp), PARAMETER :: msm_zn_km_hats_nM           = 2.3_wp

  ! Record 140: realized daily and season-cumulative supply terms (mg element
  ! m-2 [d-1]), 0.0 when the record is off. Realized = raw term x the
  ! pool-limitation factor, the same convention as M_uptake_mf and M_uptake_ac.
  REAL(wp) :: M_uptake_diff(NMIN) = 0.0_wp   ! today's realized diffusion-term uptake
  REAL(wp) :: M_uptake_hats(NMIN) = 0.0_wp   ! today's realized zinc-transport-term uptake (0 for the other elements)
  REAL(wp) :: msm_s_mf_cum(NMIN)   = 0.0_wp  ! season cumulative realized mass-flow uptake
  REAL(wp) :: msm_s_diff_cum(NMIN) = 0.0_wp  ! season cumulative realized diffusion-term uptake
  REAL(wp) :: msm_s_hats_cum(NMIN) = 0.0_wp  ! season cumulative realized zinc-transport-term uptake

  ! State (mg element m-2)
  REAL(wp) :: M_soil_labile(NMIN) = 0.0_wp
  REAL(wp) :: M_root(NMIN) = 0.0_wp, M_leaf(NMIN) = 0.0_wp, M_stem(NMIN) = 0.0_wp
  REAL(wp) :: M_chaff(NMIN) = 0.0_wp, M_grain(NMIN) = 0.0_wp
  ! Dead-tissue mineral pool. It is not a column of the daily mineral output
  ! file, but it is inside M_balance_resid, so the printed residual closes.
  REAL(wp) :: M_dead(NMIN) = 0.0_wp
  ! Daily diagnostics (mg m-2 d-1)
  REAL(wp) :: M_uptake_mf(NMIN) = 0.0_wp, M_uptake_ac(NMIN) = 0.0_wp
  REAL(wp) :: M_remob_grain(NMIN) = 0.0_wp, M_balance_resid(NMIN) = 0.0_wp
  REAL(wp) :: M_soln(NMIN) = 0.0_wp         ! record 131: today's soil-solution or xylem concentration, mg L-1
  ! Written 1 on a day when a grain bound cut the loading and 0 otherwise; the
  ! sentinel -999 for an element that is not modelled, like every other output.
  REAL(wp) :: M_cap_binding(NMIN) = 0.0_wp

  ! Internal trackers (module-own state; not C/N/water/energy state)
  REAL(wp) :: w_dead_prev = 0.0_wp     ! last-seen SHARED W_dead, for the senescence delta
  ! Record 171: yesterday's grain nitrogen pool, held as this module's own state
  ! so the daily increment is formed here and nothing is written to the nitrogen
  ! chain.
  REAL(wp) :: mnc_ngrain_prev = 0.0_wp ! g N m-2, last-seen MAX(0, N_storage - N_chaff)
  LOGICAL  :: first_call  = .TRUE.
  INTEGER  :: daa_min     = 0          ! completed days after anthesis (a_root law)
  LOGICAL  :: anth_seen   = .FALSE.    ! anthesis-day latch (crossing day = DAA 0)

  PUBLIC :: NMIN, mineral_names
  PUBLIC :: read_mineral_params, init_mineral_pools, update_minerals_daily
  PUBLIC :: M_soil_labile, M_root, M_leaf, M_stem, M_chaff, M_grain, M_dead
  PUBLIC :: M_uptake_mf, M_uptake_ac, M_remob_grain, M_balance_resid, M_cap_binding
  PUBLIC :: mm_grain_max, M_cap_refused, M_cap_days                      ! record 126 reporting
  PUBLIC :: mm_grain_eN, mnc_refused, mnc_days, mnc_to_grain              ! record 171 reporting (run log only; no output column added)
  PUBLIC :: M_soln, mm_avail_A, mm_xylem_ca_cal                          ! record 131 reporting
  PUBLIC :: M_uptake_diff, M_uptake_hats, msm_s_mf_cum, msm_s_diff_cum, msm_s_hats_cum, mm_b  ! record 140 reporting
  PUBLIC :: mm_kd, mm_ceiling, mm_pool_mgkg, mm_modelled                  ! record 183 reporting: the column-3 Kd, the column-4 ceiling and the column-2 pool from which the buffer power b_e = theta + rho_b*Kd_e is built (read-only in main_MMCO2Crop.f90)

CONTAINS

  ! Positional reader: comment lines and the 'name'-led header row are skipped;
  ! one row per element in the order Fe,Zn,Ca,Cu,Mg,Mn, the name field checked
  ! against mineral_names.
  SUBROUTINE read_mineral_params(fname, ios, ncols_min)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(OUT) :: ios
    ! Worst case (the minimum over the element rows) of the number of trailing
    ! columns this file actually carried: 20, 19, 18, 17, 16 or 15. It is a
    ! diagnostic of the file itself, never a default or fallback value, and no
    ! numeric path depends on it. OPTIONAL, so the existing caller need not ask
    ! for it.
    INTEGER, INTENT(OUT), OPTIONAL :: ncols_min
    INTEGER :: u, i, ncols_row, ncols_worst
    CHARACTER(LEN=1024) :: line
    CHARACTER(LEN=8) :: nm
    REAL(wp) :: v(20)
    ncols_worst = 20
    OPEN(newunit=u, file=TRIM(fname), status='old', iostat=ios)
    IF (ios /= 0) THEN
      IF (PRESENT(ncols_min)) ncols_min = 0
      RETURN
    END IF
    i = 0
    DO
      READ(u,'(A)',iostat=ios) line
      IF (ios /= 0) EXIT
      line = ADJUSTL(line)
      IF (LEN_TRIM(line) == 0 .OR. line(1:1) == '!') CYCLE
      IF (line(1:4) == 'name') CYCLE                 ! header row
      i = i + 1
      IF (i > NMIN) EXIT
      ! The trailing columns are OPTIONAL: the read is tried at the full width
      ! and falls back one column at a time to the 15-column form that predates
      ! them. A file written before a column existed therefore reads that column
      ! as 0.0, which is its "no bound" / "not supplied" value and the
      ! byte-identical off behavior of the record that owns it.
      v(16) = 0.0_wp; v(17) = 0.0_wp; v(18) = 0.0_wp; v(19) = 0.0_wp
      v(20) = 0.0_wp                                 ! record 171 column
      ncols_row = 20
      READ(line, *, iostat=ios) nm, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
                                v(9), v(10), v(11), v(12), v(13), v(14), v(15), v(16), v(17), v(18), v(19), v(20)
      IF (ios /= 0) THEN
        v(20) = 0.0_wp
        ncols_row = 19
        READ(line, *, iostat=ios) nm, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
                                  v(9), v(10), v(11), v(12), v(13), v(14), v(15), v(16), v(17), v(18), v(19)
      END IF
      IF (ios /= 0) THEN
        v(19) = 0.0_wp; v(20) = 0.0_wp
        ncols_row = 18
        READ(line, *, iostat=ios) nm, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
                                  v(9), v(10), v(11), v(12), v(13), v(14), v(15), v(16), v(17), v(18)
      END IF
      IF (ios /= 0) THEN
        v(17) = 0.0_wp; v(18) = 0.0_wp; v(19) = 0.0_wp; v(20) = 0.0_wp
        ncols_row = 17
        READ(line, *, iostat=ios) nm, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
                                  v(9), v(10), v(11), v(12), v(13), v(14), v(15), v(16), v(17)
      END IF
      IF (ios /= 0) THEN
        v(16) = 0.0_wp; v(17) = 0.0_wp; v(18) = 0.0_wp; v(19) = 0.0_wp; v(20) = 0.0_wp
        ncols_row = 16
        READ(line, *, iostat=ios) nm, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
                                  v(9), v(10), v(11), v(12), v(13), v(14), v(15), v(16)
      END IF
      IF (ios /= 0) THEN
        v(16) = 0.0_wp; v(17) = 0.0_wp; v(18) = 0.0_wp; v(19) = 0.0_wp; v(20) = 0.0_wp
        ncols_row = 15
        READ(line, *, iostat=ios) nm, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), &
                                  v(9), v(10), v(11), v(12), v(13), v(14), v(15)
      END IF
      IF (ios /= 0) EXIT
      ncols_worst = MIN(ncols_worst, ncols_row)
      IF (TRIM(nm) /= TRIM(mineral_names(i))) THEN   ! row order guard
        ios = 2
        EXIT
      END IF
      mm_modelled(i)   = NINT(v(1))
      mm_pool_mgkg(i)  = v(2);  mm_kd(i)        = v(3);  mm_ceiling(i)   = v(4)
      mm_vmax_hats(i)  = v(5);  mm_km_hats(i)   = v(6)
      mm_vmax_lats(i)  = v(7);  mm_km_lats(i)   = v(8)
      mm_ps_rate(i)    = v(9);  mm_ps_rate_def(i) = v(10)
      mm_fe_def_switch(i) = v(11)
      mm_export(i)     = MIN(1.0_wp, MAX(0.0_wp, v(12)))
      mm_xylem(i)      = v(13)
      mm_gwr(i)        = MAX(0.0_wp, v(14))
      ! Driver column 15, the field-anchored active-uptake rate; the MAX is a
      ! structural non-negativity floor. Value provenance: the driver zinc row.
      mm_field_rate(i) = MAX(0.0_wp, v(15))
      mm_grain_max(i)  = MAX(0.0_wp, v(16))      ! record 126: maximum grain concentration, mg kg-1
      mm_avail_A(i)      = MAX(0.0_wp, v(17))    ! record 131: site-level availability, mg L-1
      mm_xylem_ca_cal(i) = MAX(0.0_wp, v(18))    ! record 131: xylem calcium concentration, mg L-1
      ! Driver column 19, the buffer power b_e. A non-positive or absent value
      ! falls back to the declared default 1.0, never silently to 0.0, which
      ! would divide by zero in the effective diffusion coefficient below.
      mm_b(i) = v(19)
      IF (mm_b(i) <= 0.0_wp) mm_b(i) = 1.0_wp
      ! Driver column 20, the grain element:nitrogen stoichiometry r_e. A
      ! non-positive or absent value means NO BOUND - the record is then a
      ! provable no-op for that element - never a zero ceiling that would refuse
      ! all grain loading.
      mm_grain_eN(i) = MAX(0.0_wp, v(20))
    END DO
    CLOSE(u)
    IF (PRESENT(ncols_min)) ncols_min = ncols_worst
    IF (ios == 2) RETURN                             ! name-order mismatch
    ios = MERGE(0, 1, i == NMIN)
  END SUBROUTINE read_mineral_params

  SUBROUTINE init_mineral_pools()
    INTEGER :: i
    DO i = 1, NMIN
      IF (mm_modelled(i) == 1) THEN
        ! Labile pool: the measured extractable content times the soil mass.
        M_soil_labile(i) = mm_pool_mgkg(i) * soil_mass_kg_m2
        M_root(i) = 0.0_wp; M_leaf(i) = 0.0_wp; M_stem(i) = 0.0_wp
        M_chaff(i) = 0.0_wp; M_grain(i) = 0.0_wp; M_dead(i) = 0.0_wp
        M_uptake_mf(i) = 0.0_wp; M_uptake_ac(i) = 0.0_wp
        M_uptake_diff(i) = 0.0_wp; M_uptake_hats(i) = 0.0_wp   ! record 140
        M_remob_grain(i) = 0.0_wp; M_balance_resid(i) = 0.0_wp
        M_cap_binding(i) = 0.0_wp
      END IF
      msm_s_mf_cum(i) = 0.0_wp; msm_s_diff_cum(i) = 0.0_wp; msm_s_hats_cum(i) = 0.0_wp   ! record 140 season reset
      mnc_refused(i) = 0.0_wp; mnc_days(i) = 0; mnc_to_grain(i) = 0.0_wp                 ! record 171 season reset
    END DO
    w_dead_prev = 0.0_wp
    mnc_ngrain_prev = 0.0_wp             ! record 171 season reset
    ! Record 180: season reset of the calcium ledgers. Reset unconditionally;
    ! they are only ever written inside the record-180 branch and only ever read
    ! by the end-of-run summary, so an off run is unaffected.
    gcl_ca_grain_cum  = 0.0_wp; gcl_ca_chaff_cum = 0.0_wp; gcl_ca_veg_cum = 0.0_wp
    gcl_tr_season_cum = 0.0_wp; gcl_tr_post_cum  = 0.0_wp
    gcl_days_post     = 0;      gcl_days_pool    = 0
    first_call  = .TRUE.
    daa_min     = 0
    anth_seen   = .FALSE.
  END SUBROUTINE init_mineral_pools

  ! vwc_root enters the effective diffusion coefficient and its wetness
  ! impedance below; root_max_depth is part of the interface and is not used.
  ! dW_grain is the PHLOEM-SHARE grain term, the quantity record 114 re-bases in
  ! the caller. dW_grain_dm is the GRAIN DRY-MATTER INCREMENT and drives only the
  ! calcium xylem route, because the grain water import of Lizana et al. 2010
  ! (J Exp Bot 61:1147, DOI 10.1093/jxb/erp380, Tables 2+3) is defined per unit
  ! of FINAL GRAIN DRY MATTER, whichever carbon pool paid for it; the phloem
  ! delivery term is a strict subset of that increment and is the wrong driver
  ! for it.
  SUBROUTINE update_minerals_daily(trans_act, vwc_root, root_max_depth, W_root, W_grain, dvs2, &
                                   dW_root, dW_leaf, dW_stem, dW_chaff, dW_grain, dW_grain_dm)
    REAL(wp), INTENT(IN) :: trans_act, vwc_root, root_max_depth, W_root, W_grain, dvs2
    REAL(wp), INTENT(IN) :: dW_root, dW_leaf, dW_stem, dW_chaff, dW_grain
    REAL(wp), INTENT(IN) :: dW_grain_dm      ! the grain dry-matter increment, calcium xylem route only
    INTEGER :: i
    REAL(wp) :: soln, S_mf, S_ac, raw, uptake, bf, a_root, d_dead
    ! Record 140: diffusion and zinc-transport working scalars.
    REAL(wp) :: S_diff, S_hats, S_hats_raw, fl_theta_m, De_e, RLD_e, C_min_e, km_zn_mgL
    REAL(wp) :: fl_theta_msl   ! impedance for the retired record-131 active term
    REAL(wp) :: b_e_use, kd_e_use   ! record 183: the buffer power actually divided by, and the Kd it was built from
    REAL(wp) :: znk_q_diff, znk_c_surf   ! record 169: water-equivalent diffusive flux density and root-surface concentration
    REAL(wp) :: dWt, fr, fl, fs, fc, fg, before, after
    REAL(wp) :: dr, dl, ds, dc, dg, dg_dm, sfrac, moved, export, dead_part
    REAL(wp) :: d_dead_stem, d_dead_chaff   ! stem and chaff death fluxes (g DM m-2 d-1)
    REAL(wp) :: gwi, g_in
    REAL(wp) :: gcx_area_tot, gcx_f_ear, gcx_v_ear, gcx_v_thr, gcx_head_in, gcx_c_in, gcx_fg   ! record 166
    REAL(wp) :: gcl_g_in, gcl_c_in, gcl_v_in   ! record 180: today's calcium to grain / chaff / vegetative terminus, mg m-2 d-1
    REAL(wp) :: g_before, g_load, g_allow, g_excess, g_from_upt, g_from_exp, f_upt   ! record 126
    ! Record 171: the grain nitrogen pool, its daily increment and the
    ! protein-sink working set, all formed here from read-only SHARED state.
    REAL(wp) :: mnc_ngrain, mnc_dngrain, mnc_allow, mnc_excess                      ! record 171
    LOGICAL  :: growth_day

    ! ---- Senescence signal: the dead-leaf dry-matter transfer, read from
    !      SHARED; one day lagged, see the module header ------------------ !
    IF (first_call) THEN
      w_dead_prev = W_dead_sh          ! W_dead may be nonzero pre-season
      first_call  = .FALSE.
    END IF
    d_dead = MAX(0.0_wp, W_dead_sh - w_dead_prev)
    w_dead_prev = W_dead_sh
    ! Record 63: the stem and chaff death fluxes, read from the per-organ daily
    ! ledgers the organ-nitrogen turnover routine writes (g DM m-2 d-1). That
    ! routine runs later in the day than this one, so these hold the PREVIOUS
    ! day's values - the same explicit one-day senescence lag as d_dead above -
    ! and each day's flux is consumed exactly once, this routine being called
    ! once per day. Both are identically 0.0 when the record is off, which makes
    ! the stem and chaff export block below a no-op.
    d_dead_stem  = 0.0_wp
    d_dead_chaff = 0.0_wp
    IF (grainN_partition == 1) THEN
      d_dead_stem  = MAX(0.0_wp, n2_dW_dead_stem_day)
      d_dead_chaff = MAX(0.0_wp, n2_dW_dead_chaff_day)
    END IF

    ! ---- Living-root activity fraction after anthesis: the same law and the
    !      same anthesis latch as the nitrogen-uptake block of
    !      main_MMCO2Crop.f90, on the same driver constants (Zheng et al. 2022,
    !      Front Plant Sci 12:814658, DOI 10.3389/fpls.2021.814658): a_root = 1
    !      to anthesis, then MAX(aroot_floor, 1 - aroot_decline_d*DAA), the
    !      crossing day itself DAA = 0 ------------------------------------ !
    IF (dvs2 >= 1.0_wp) THEN
      IF (anth_seen) daa_min = daa_min + 1
      anth_seen = .TRUE.
    END IF
    a_root = 1.0_wp
    IF (daa_min > 0) a_root = MAX(aroot_floor, 1.0_wp - aroot_decline_d * REAL(daa_min, wp))

    ! Shares of the day's dry-matter increment, each clipped at zero, with the
    ! zero-delivery fallback below. dg carries the PHLOEM share of the grain and
    ! is the term record 114 re-bases in the caller; dg_dm is the grain
    ! dry-matter increment and is used only by the calcium xylem route, on either
    ! setting of record 114. Full provenance of the share basis, and the proof
    ! that every share is non-negative, at SHARED.f90.
      dr = MAX(0.0_wp, dW_root); dl = MAX(0.0_wp, dW_leaf); ds = MAX(0.0_wp, dW_stem)
      dc = MAX(0.0_wp, dW_chaff); dg = MAX(0.0_wp, dW_grain)
      dg_dm = MAX(0.0_wp, dW_grain_dm)
    dWt = dr + dl + ds + dc + dg
    growth_day = dWt > 1.0e-9_wp
    IF (growth_day) THEN
      fr = dr/dWt; fl = dl/dWt; fs = ds/dWt; fc = dc/dWt; fg = dg/dWt
    ELSE
      ! Zero-delivery fallback: the day's uptake accrues to the LEAF, the
      ! transpiration-stream terminus; the shares still sum to exactly 1 and
      ! nothing is divided by zero.
      fr = 0.0_wp; fl = 1.0_wp; fs = 0.0_wp; fc = 0.0_wp; fg = 0.0_wp
    END IF

    ! ---- Record 171: today's grain nitrogen, read-only ----------------- !
    ! N_grain = MAX(0, N_storage - N_chaff) is the same expression
    ! main_MMCO2Crop.f90 uses for the grain nitrogen concentration and the grain
    ! protein, so this record introduces no new definition of grain nitrogen and
    ! writes neither name. The increment is formed against this module's own
    ! saved previous value, so no nitrogen state is touched.
    mnc_ngrain  = MAX(0.0_wp, N_storage_sh - N_chaff_sh)
    mnc_dngrain = MAX(0.0_wp, mnc_ngrain - mnc_ngrain_prev)

    DO i = 1, NMIN
      IF (mm_modelled(i) /= 1) CYCLE

      before = M_root(i) + M_leaf(i) + M_stem(i) + M_chaff(i) + M_grain(i) + M_dead(i)
      g_before = M_grain(i)                       ! record 126: grain pool before today's loading

      ! ---- 2. Soil-solution concentration (mg/L) ---------------------- !
      ! Four laws, in the order the branches test them.
      !  Calcium under record 180: `soln` is the XYLEM-SAP concentration (driver
      !    column 18), not a soil concentration. Calcium is delivered to the ear
      !    by the transpiration stream, and only the extreme root tip and the
      !    regions where lateral roots are initiated load it, the endodermal
      !    Casparian band being absent or disrupted there: White PJ & Broadley MR
      !    (2003) Ann Bot 92:487-511, DOI 10.1093/aob/mcg164, p.489, "The
      !    delivery of Ca to the xylem is restricted to the extreme root tip and
      !    to regions in which lateral roots are being initiated". The soil
      !    calcite ceiling is therefore not a xylem value. The shipped
      !    15 mg L-1 lies at the lower end of the measured 12-660 mg L-1 range;
      !    the window is stated and enforced in main_MMCO2Crop.f90 before the
      !    first day runs.
      !  Fe, Zn, Cu and Mn under record 131: `soln` is the site-level
      !    availability parameter A_e (driver column 17), calibrated on one
      !    control treatment per site. It carries the solution chemistry of a
      !    site with no measured soil chemistry together with the root-diffusion
      !    geometry no specific-root-length state in this module could supply, so
      !    it is an availability index rather than a literal concentration.
      !  Magnesium: the linear sorption partition [M] = (pool / soil mass) / Kd,
      !    Kd in L kg-1, so the solution concentration falls in proportion as the
      !    pool depletes. The partition is a SOLID-PHASE equilibrium and carries
      !    no soil-moisture term.
      !  Iron's own column-4 ceiling: a fixed solution concentration while the
      !    pool lasts, the oxide-solubility and complexation control of a
      !    calcareous soil solution rather than a sorption partition.
      ! For calcium the active and diffusive terms are identically zero, so the
      ! day's raw supply is the xylem concentration times the day's
      ! transpiration, with no second term and no multiplier.
      IF (grain_ca_loading_basis == 1 .AND. i == idx_Ca) THEN
        soln = mm_xylem_ca_cal(i)
      ELSE IF (mineral_supply_law == 1 .AND. (i == idx_Fe .OR. i == idx_Zn .OR. i == idx_Cu .OR. i == idx_Mn)) THEN
        soln = mm_avail_A(i)
      ELSE IF (mm_kd(i) > 0.0_wp) THEN
        ! Linear operational partition: [M] = pool(mg/kg)/Kd(L/kg);
        ! declines proportionally as the pool depletes.
        soln = (M_soil_labile(i) / soil_mass_kg_m2) / mm_kd(i)
      ELSE
        ! Fe: oxide-solubility ceiling (fixed conc while the pool lasts)
        soln = mm_ceiling(i)
      END IF
      IF (M_soil_labile(i) <= 0.0_wp) soln = 0.0_wp
      M_soln(i) = soln                                     ! record 131: diagnostic, every mode

      ! ---- 3. Uptake --------------------------------------------------- !
      S_mf = MAX(0.0_wp, trans_act) * soln          ! mm/d x mg/L = mg m-2 d-1
      S_ac = 0.0_wp
      IF (mm_ps_rate(i) > 0.0_wp .AND. M_soil_labile(i) > 0.0_wp) THEN
        ! Iron: phytosiderophore-mediated Fe(III) uptake per g of root dry
        ! matter per day - Roemheld & Marschner 1986, Plant Physiol 80:175,
        ! Tables VI-VII, the Fe-sufficient rate with the calcareous (pH 7.5 plus
        ! CaCO3) correction, the Fe-deficient rate selectable in the driver and
        ! never blended. A solid-phase Fe(III) source, so the rate carries no
        ! Michaelis term in the solution iron. The paper's own experiments are in
        ! barley, wheat entering through its Table VII multi-species comparison,
        ! so this is a graminaceous transfer and not a wheat measurement.
        S_ac = S_ac + mm_ps_rate(i) * MAX(0.0_wp, W_root) * a_root      &
                      * atomic_mass(i) * 1.0e-6_wp
      END IF
      IF (mm_field_rate(i) > 0.0_wp .AND. M_soil_labile(i) > 0.0_wp .AND. mineral_supply_law /= 1) THEN
        ! Record 131: zinc's field-anchored active term is RETIRED when the
        ! record is on - the gate above adds only `.AND. mineral_supply_law /= 1`
        ! - because its capacity is taken up by the saturating root-surface term
        ! of record 169 below; it stays live and unchanged when the record is off.
        ! As coded, S_ac = S_ac + mm_field_rate(i) * MAX(0, W_root) * a_root,
        ! i.e. mg gDW-1 d-1 x g DW m-2 = mg m-2 d-1, the rate being in mg
        ! directly with no atomic-mass conversion, on the same living-root basis
        ! W_root x a_root as the iron term above and as the nitrogen-uptake law.
        ! Rate, source and derivation are declared at mm_field_rate above; supply
        ! feedback is the labile-pool cap and the separate mass-flow term only.
        S_ac = S_ac + mm_field_rate(i) * MAX(0.0_wp, W_root) * a_root
      END IF
      IF (mineral_supply_law == 1 .AND. (i == idx_Fe .OR. i == idx_Zn .OR. i == idx_Cu .OR. i == idx_Mn) &
          .AND. M_soil_labile(i) > 0.0_wp .AND. mineral_supply_mech < 1) THEN
        ! -------------------------------------------------------------- !
        ! Record 131: the per-root-mass active term, generalized from the form
        ! of iron's phytosiderophore line and zinc's field-rate line above -
        ! rate x MAX(0, W_root) x a_root - to all four record-131 elements, with
        ! the rate being the SAME availability parameter A_e the mass-flow term
        ! uses, scaled by the soil-wetness impedance f_l(theta). Each element
        ! therefore stays one-dimensional for the availability calibration. This
        ! is NOT a Michaelis-Menten term: there is no Km and no saturation, and
        ! A_e is a linear availability index, not a literal concentration.
        ! Retired in turn when record 140 is on - the only change to the gate
        ! above being the trailing `.AND. mineral_supply_mech < 1`, the same
        ! retirement technique this record used on zinc's field rate - and
        ! superseded by the diffusion and zinc-transport terms below. The
        ! impedance it uses is the one of record 168, so the same factor appears
        ! at both places. SHARED.f90.
        ! -------------------------------------------------------------- !
        IF (mineral_impedance_basis == 1) THEN
          fl_theta_msl = nud_csif_slope * MAX(0.0_wp, vwc_root - nud_theta_th)
        END IF
        S_ac = S_ac + mm_avail_A(i) * fl_theta_msl &
                      * MAX(0.0_wp, W_root) * a_root
      END IF
      ! ================================================================ !
      ! Record 140: DIFFUSION TO THE ROOT SURFACE for every modelled element
      ! except calcium, plus the SATURATING ZINC TRANSPORT term - the same
      ! three-component supply law as the nitrogen uptake. The full law, every
      ! constant's source and the unit-chain closure are at SHARED.f90
      ! Off: this block does not execute and S_diff and S_hats
      ! stay at the reset below, so `raw` is a bit-exact S_mf + S_ac + 0.0.
      ! ================================================================ !
      S_diff = 0.0_wp; S_hats = 0.0_wp
      IF (mineral_supply_mech >= 1 .AND. i /= idx_Ca .AND. M_soil_labile(i) > 0.0_wp) THEN
        ! Soil-wetness impedance of the diffusion path: f_l(theta) =
        ! 1.1 (theta - theta_th) above the threshold water content and 0 below
        ! it - Olesen, Moldrup, Yamaguchi & Rolston 2001, Soil Sci 166:89-96,
        ! DOI 10.1097/00010694-200102000-00002, printed p.94 Eq. (7a)/(7b) - with
        ! theta_th from the same paper's printed p.93 Eq. (5) and this deck's own
        ! measured clay fraction, silt fraction and topsoil bulk density. Both
        ! constants are those of the nitrogen-uptake law, reused verbatim. The
        ! line below, De_e = Dl theta f_l / b, is Olesen Eq. (8a) divided by the
        ! buffer power, so the functional form is the paper's own.
        IF (mineral_impedance_basis == 1) THEN
          fl_theta_m = nud_csif_slope * MAX(0.0_wp, vwc_root - nud_theta_th)
          ! Season ledgers only, never a flux: the alternative factor is
          ! accumulated on the same days so the run log can report the ratio.
          mib_fl_sum     = mib_fl_sum + fl_theta_m
          mib_fl_off_sum = mib_fl_off_sum + MAX(0.0_wp, 1.60_wp * vwc_root - 0.17_wp)
          mib_fl_days    = mib_fl_days + 1
          IF (vwc_root <= nud_theta_th) mib_days_at_zero = mib_days_at_zero + 1
          mib_vwc_min    = MIN(mib_vwc_min, vwc_root)
          mib_vwc_max    = MAX(mib_vwc_max, vwc_root)
        END IF
        ! Effective diffusion coefficient De_e (m2 s-1), the Barber-Cushman
        ! form: the free-solution coefficient times the water content times the
        ! impedance, divided by the buffer power.
        ! ---------------------------------------------------------------- !
        ! THE BUFFER POWER THAT De_e IS DIVIDED BY.
        ! Record 183 off: b_e_use is the driver column-19 constant this line has
        ! always divided by, 1.0 on every shipped row and declared unsourced by
        ! record 140 itself, so the divisor is bit-for-bit what it was.
        ! On: the LINEAR-ISOTHERM buffer power, recomputed every day for every
        ! element on this module's single declared 0-0.6 m layer,
        !     b_e = theta + rho_b * Kd_e
        ! Baldwin JP, Nye PH & Tinker PB (1973) Plant and Soil 38:621-635,
        ! DOI 10.1007/BF00010701, printed p.622 list of symbols: "b = buffer
        ! power of soil = dC/dCl", with C the total concentration of diffusible
        ! solute per unit soil volume and Cl its concentration in solution; and
        ! printed p.625, "Note that if only solution phase ions are mobile then
        ! Db = Dl.W.fl", which is D = Dl*theta*f_l/b and exactly the form of the
        ! line above - the FUNCTIONAL FORM DOES NOT CHANGE, only the value of b.
        ! The split of b into its solution and adsorbed parts is printed by
        ! Roose T & Fowler AC (2004) J Theor Biol 228:173-184,
        ! DOI 10.1016/j.jtbi.2003.10.008, Eq. (2.8) p.175, "d/dt[(b + phi*S)c]",
        ! "where b is the buffer power of the soil", with phi*S = theta: the
        ! solution part is theta and the adsorbed part is the remainder, which
        ! for the LINEAR isotherm that driver column 3 already asserts
        ! ([M]_soln = pool_mg_kg / Kd) is rho_b*Kd per unit soil VOLUME. That b
        ! is far above 1 for an adsorbed ion is printed in their Table 1 p.176,
        ! "b  Nutrient buffer power in soil (for phosphate)  239".
        ! NO NEW INPUT AND NO NEW CONSTANT: theta is the vwc_root this line
        ! already uses twice, rho_b is the deck's own measured bulk density
        ! already read for the nitrogen-uptake law, and Kd_e is the element's own
        ! column-3 value.
        b_e_use = mm_b(i)
        IF (mineral_buffer_power_basis == 1) THEN
          ! rho_b in kg L-1 (numerically Mg m-3 = g cm-3). The fallback is this
          ! module's own soil-mass basis, soil_mass_kg_m2 / (0.6 m x 1000) =
          ! 1.63 kg L-1, the same conversion the labile pool is built on; it is
          ! printed, never silent.
          IF (nud_bulk_density > 0.0_wp) THEN
            bpb_rho_b_used = nud_bulk_density
          ELSE
            bpb_rho_b_used = soil_mass_kg_m2 / (mm_depth_m * 1000.0_wp)
            bpb_rho_fallback = .TRUE.
          END IF
          ! Kd_e in L kg-1: driver column 3 where positive (Zn 645 after
          ! Duffner 2012 and Cu 1200 after Bravin et al. 2010, both operational;
          ! Mg 60 derived; Mn from the pH regression of Sheppard, Long &
          ! Sanipelli 2009, SKB report R-09-27, Table 5-29). Iron's column 3 is 0.0 by this file's own ceiling
          ! convention, so its Kd is recovered from the definition
          ! Kd = pool / [M]_soln on the iron row's two sourced numbers, giving an
          ! operational equivalent of the same class as zinc's and copper's.
          ! Calcium never reaches here.
          kd_e_use = mm_kd(i)
          IF (kd_e_use <= 0.0_wp .AND. mm_ceiling(i) > 0.0_wp) &
            kd_e_use = mm_pool_mgkg(i) / mm_ceiling(i)
          b_e_use = vwc_root + bpb_rho_b_used * kd_e_use
          ! Season ledgers only, never a flux and never read back into physics.
          bpb_kd_used(i) = kd_e_use
          bpb_b_sum(i)   = bpb_b_sum(i) + b_e_use
          bpb_b_min(i)   = MIN(bpb_b_min(i), b_e_use)
          bpb_b_max(i)   = MAX(bpb_b_max(i), b_e_use)
          bpb_b_days(i)  = bpb_b_days(i) + 1
          IF (b_e_use < bpb_warn_b_lo) bpb_warn_lo(i) = bpb_warn_lo(i) + 1
        END IF
        De_e  = mm_Dl_e(i) * vwc_root * fl_theta_m / b_e_use
        ! Root length density (m root per m3 of soil) from the specific root
        ! length of the nitrogen-uptake law over the 0-0.6 m pool depth, which
        ! stands in for that law's per-layer thickness.
        RLD_e = W_root * nud_SRL / mm_depth_m
        ! Minimum solution concentration C_min,e: the element's own apparent Km
        ! where wheat kinetics exist, which is zinc only, and 0.0 elsewhere, no
        ! wheat Km for Fe, Cu, Mg or Mn being available.
        C_min_e = 0.0_wp
        IF (i == idx_Zn) C_min_e = msm_zn_km_hats_nM * atomic_mass(idx_Zn) * 1.0e-6_wp
        ! Unit-chain closure (SHARED.f90, the lumped Barber-Cushman
        ! volumetric-sink-times-depth convention the nitrogen uptake also uses):
        ! RLD_e [m m-3] x De_e [m2 s-1] x mm_depth_m [m] x 86400 [s d-1] is a
        ! water-equivalent flux density in m d-1, and x 1000 converts m to mm,
        ! i.e. to the L m-2 d-1 of the mass-flow line above, so multiplying by
        ! the solution concentration in mg L-1 gives mg m-2 d-1. A structural
        ! unit conversion, no new physics.
        S_diff = MAX(0.0_wp, soln - C_min_e) * RLD_e * De_e * mm_depth_m * 86400.0_wp * 1.0e3_wp
        IF (i == idx_Zn) THEN
          ! Zinc's apparent Km, converted nM -> mg L-1: Km_nM x atomic mass
          ! (g mol-1) x 1e-6, since 1 nM is 1e-9 mol L-1.
          km_zn_mgL = msm_zn_km_hats_nM * atomic_mass(idx_Zn) * 1.0e-6_wp
          IF (zn_uptake_kinetic_basis >= 1) THEN
            ! ============================================================ !
            ! Record 169: THE ROOT-SURFACE-CONCENTRATION ZINC KINETIC.
            ! q_diff is the water-equivalent diffusive flux density that the
            ! diffusion line immediately above already multiplies its
            ! concentration gradient by, reused token for token (L m-2 d-1).
            ! c_surf is the bulk solution concentration MINUS the Barber-Cushman
            ! depletion that same line has already charged (S_diff / q_diff
            ! returns MAX(0, soln - C_min_e) in mg L-1), so the diffusion limb
            ! and the transporter limb speak about ONE root-surface
            ! concentration instead of two. Where q_diff is not positive no
            ! diffusive depletion has been charged and c_surf is the bulk
            ! concentration; that is also the guard against division by zero
            ! below the impedance threshold or at zero root mass.
            ! The term's CAPACITY is the field rate of Liu et al. 2019,
            ! Front Plant Sci 10:426, DOI 10.3389/fpls.2019.00426 - driver
            ! column 15 - on the same living-root basis W_root x a_root as
            ! iron's phytosiderophore term, SATURATED by the apparent Km of
            ! Hacisalihoglu et al. 2001, Plant Physiol 125:456,
            ! DOI 10.1104/pp.125.1.456, Table II p.459, converted just above.
            ! This is a DECLARED PAIRING of a measured shape with a measured
            ! scale: a seedling Km supplies the transporter's identity and a
            ! field rate its capacity, no field-wheat Michaelis-Menten zinc
            ! influx kinetic per unit root mass or surface area being available.
            ! Neither number is fitted to an observation and no new numeral is
            ! introduced: 86400 and 1.0e3 are the structural unit conversions of
            ! the diffusion line and 0.0 is a structural zero.
            ! Because the diffusion term is computed down to a minimum
            ! concentration equal to Km, c_surf equals Km on every day with
            ! diffusive supply and the term then runs at half its capacity.
            ! SHARED.f90.
            ! ============================================================ !
            znk_q_diff = RLD_e * De_e * mm_depth_m * 86400.0_wp * 1.0e3_wp
            IF (znk_q_diff > 0.0_wp) THEN
              znk_c_surf = soln - S_diff / znk_q_diff
            ELSE
              znk_c_surf = soln
              znk_days_qdiff0 = znk_days_qdiff0 + 1
            END IF
            IF (znk_c_surf > 0.0_wp) THEN
              S_hats = mm_field_rate(idx_Zn) * MAX(0.0_wp, W_root) * a_root &
                       * znk_c_surf / (km_zn_mgL + znk_c_surf)
              znk_sat_frac_sum = znk_sat_frac_sum + znk_c_surf / (km_zn_mgL + znk_c_surf)
            ELSE
              S_hats = 0.0_wp                     ! guard, stated above
              znk_days_csurf0 = znk_days_csurf0 + 1
            END IF
            znk_days_active = znk_days_active + 1
            znk_s_hats_cum  = znk_s_hats_cum + S_hats
          END IF
        END IF
      END IF
      S_ac = S_ac + S_hats
      raw    = S_mf + S_ac + S_diff
      uptake = MAX(0.0_wp, MIN(raw, M_soil_labile(i)))   ! pool is the only cap
      bf = 0.0_wp
      IF (raw > 0.0_wp) bf = uptake / raw
      M_uptake_mf(i) = S_mf * bf
      M_uptake_ac(i) = S_ac * bf
      M_uptake_diff(i) = S_diff * bf   ! record 140: 0.0 when off
      M_uptake_hats(i) = S_hats * bf   ! record 140: 0.0 when off or not zinc
      IF (mineral_supply_mech >= 1) THEN
        msm_s_mf_cum(i)   = msm_s_mf_cum(i)   + M_uptake_mf(i)
        msm_s_diff_cum(i) = msm_s_diff_cum(i) + M_uptake_diff(i)
        msm_s_hats_cum(i) = msm_s_hats_cum(i) + M_uptake_hats(i)
      END IF
      M_soil_labile(i) = M_soil_labile(i) - uptake

      ! ---- 4(i). Allocation of current uptake --------------------------- !
      ! The 0.5 test is the 0/1 driver-flag threshold on the xylem_route column,
      ! a structural comparison and not a physical constant.
      ! ---- Record 180: THE CALCIUM ROUTE, AS CODED -------------------- !
      ! Calcium is deposited irreversibly where the transpiration stream leaves
      ! the xylem. From anthesis the head's measured share of the day's arrival
      ! is divided at the grain-base xylem discontinuity between grain and chaff,
      ! and the remainder goes to the vegetative terminus:
      !   grain += f_grain_post                  x uptake   (dvs2 >= 1)
      !   chaff += (f_head_post - f_grain_post)  x uptake   (dvs2 >= 1)
      !   leaf  += (1 - f_head_post)             x uptake   (dvs2 >= 1)
      !   leaf  += uptake                                   (before anthesis)
      ! `uptake` is MIN(c_xyl,Ca x Tr_act, M_soil_labile), the day's arrival
      ! after the same labile-pool bound every element carries, so the three
      ! termini sum to it exactly and the element balance is unchanged in form.
      ! THE CONSTANTS, nothing fitted: Hocking PJ (1994) J Plant Nutr
      ! 17:1289-1308, DOI 10.1080/01904169409364807, Tables 1(b) and 1(c) p.1294
      ! with Table 2 p.1295 give f_head_post = 0.583, the head's share of
      ! post-anthesis calcium arrival, and Table 1(c)'s grain 10.8 / (10.8 + 9.8)
      ! = 0.524 gives f_grain_post = 0.305; the same paper's Table 6 p.1297
      ! reproduces the ratio independently. SHARED.f90.
      ! THE ROUTE: Zee SY & O'Brien TP (1970) Aust J Biol Sci 23:783-791,
      ! DOI 10.1071/BI9700783, pp.783-784 - the xylem is discontinuous at the
      ! grain base while the phloem is continuous, so the grain's OWN water
      ! import cannot carry a phloem-immobile element, which is why the grain
      ! water-import term and the ear-area share of record 166 below are
      ! superseded here.
      ! NO REMOBILIZATION: the calcium export fraction is 0.0 in every deck's
      ! driver, so the senescence-export limb below is already a no-op for
      ! calcium and this branch adds no clause to it (White & Broadley 2003
      ! p.488; Hocking 1994 Table 4 p.1296).
      ! Anthesis is dvs2 >= 1.0, this routine's own existing latch; the routine
      ! is not called after maturity, so the window closes there.
      IF (grain_ca_loading_basis == 1 .AND. i == idx_Ca) THEN
        gcl_tr_season_cum = gcl_tr_season_cum + MAX(0.0_wp, trans_act)
        IF (uptake < raw - 1.0e-12_wp) gcl_days_pool = gcl_days_pool + 1
        IF (dvs2 >= 1.0_wp) THEN
          gcl_tr_post_cum = gcl_tr_post_cum + MAX(0.0_wp, trans_act)
          gcl_days_post   = gcl_days_post + 1
          gcl_g_in = gcl_f_grain_post * uptake
          gcl_c_in = (gcl_f_head_post - gcl_f_grain_post) * uptake
          gcl_v_in = uptake - gcl_g_in - gcl_c_in
          M_grain(i) = M_grain(i) + gcl_g_in
          M_chaff(i) = M_chaff(i) + gcl_c_in
          M_leaf(i)  = M_leaf(i)  + gcl_v_in
          gcl_ca_grain_cum = gcl_ca_grain_cum + gcl_g_in
          gcl_ca_chaff_cum = gcl_ca_chaff_cum + gcl_c_in
          gcl_ca_veg_cum   = gcl_ca_veg_cum   + gcl_v_in
        ELSE
          ! Before anthesis there is no head: the whole day's arrival goes to
          ! the vegetative terminus, the leaf, this module carrying no
          ! organ-level transpiration state.
          M_leaf(i)      = M_leaf(i)      + uptake
          gcl_ca_veg_cum = gcl_ca_veg_cum + uptake
        END IF
      ELSE IF (mm_xylem(i) >= 0.5_wp) THEN
        ! Record 166: the earlier calcium route, superseded by record 180 above
        ! wherever that record is on. No phloem entry: the grain receives its
        ! xylem WATER IMPORT times the solution concentration, the import being
        ! 0.7625 g H2O per g of final grain dry matter (Lizana et al. 2010,
        ! J Exp Bot 61:1147, DOI 10.1093/jxb/erp380, Tables 2+3) times the GRAIN
        ! DRY-MATTER INCREMENT, since the ratio is defined per unit of final
        ! grain dry matter. A LOWER BOUND: transpirational grain-water throughput
        ! is excluded. bf scales for the pool cap - never binding for calcium -
        ! consistently with `uptake`.
        gwi  = mm_gwr(i) * dg_dm * 1.0e-3_wp         ! g/gDM x g/m2/d /1000 = L m-2 d-1
        IF (grain_ca_xylem_throughput == 1) THEN
          ! The ear xylem-throughput route. The grain's net import is wholly the
          ! grain's; the rest of the xylem water the ear receives is divided at
          ! the grain-base xylem discontinuity (Zee & O'Brien 1970) by the
          ! Hocking 1994 Table 1(c) p.1294 measured partition. The ear's share of
          ! the transpiration stream is its share of the model's own standing
          ! green area, on one basis, with the measured ear area - no new area
          ! constant and no multiplier.
          gcx_area_tot = MAX(0.0_wp, cls_gai) + MAX(0.0_wp, cse_sai) + &
                         MAX(0.0_wp, cse_iai) + MAX(0.0_wp, cse_eai_meas)
          gcx_f_ear = 0.0_wp
          IF (gcx_area_tot > 0.0_wp) gcx_f_ear = MAX(0.0_wp, cse_eai_meas) / gcx_area_tot
          gcx_v_ear = gcx_f_ear * MAX(0.0_wp, trans_act)            ! L m-2 d-1
          gcx_v_thr = MAX(0.0_wp, gcx_v_ear - gwi)                  ! L m-2 d-1
          gcx_fg    = gcx_hk_grain / (gcx_hk_grain + gcx_hk_chaff)  ! 0.524272
          g_in        = MIN(uptake, (gwi + gcx_fg * gcx_v_thr) * soln * bf)
          gcx_c_in    = MIN(MAX(0.0_wp, uptake - g_in), &
                            (1.0_wp - gcx_fg) * gcx_v_thr * soln * bf)
          gcx_head_in = g_in + gcx_c_in
          M_grain(i)  = M_grain(i) + g_in
          M_chaff(i)  = M_chaff(i) + gcx_c_in
          M_leaf(i)   = M_leaf(i)  + MAX(0.0_wp, uptake - gcx_head_in)
          ! Season ledgers (calcium is the only xylem-route element)
          gcx_v_ear_cum    = gcx_v_ear_cum    + gcx_v_ear
          gcx_v_net_cum    = gcx_v_net_cum    + gwi
          gcx_ca_grain_cum = gcx_ca_grain_cum + g_in
          gcx_ca_chaff_cum = gcx_ca_chaff_cum + gcx_c_in
          gcx_f_ear_max    = MAX(gcx_f_ear_max, gcx_f_ear)
          IF (gcx_head_in > 0.0_wp) gcx_days_route = gcx_days_route + 1
        END IF
      ELSE
        ! All-phloem entry: the day's uptake follows the organ growth shares.
        M_root(i)  = M_root(i)  + uptake * fr
        M_leaf(i)  = M_leaf(i)  + uptake * fl
        M_stem(i)  = M_stem(i)  + uptake * fs
        M_chaff(i) = M_chaff(i) + uptake * fc
        M_grain(i) = M_grain(i) + uptake * fg
      END IF

      ! ---- 4(ii). Senescence-coupled remobilization -------------------- !
      ! All grain entry other than the calcium xylem route is by PHLOEM, the
      ! xylem being discontinuous at the grain stalk (Garnett & Graham 2005,
      ! Ann Bot 95:817, DOI 10.1093/aob/mci085). This limb rides the model's own
      ! leaf-senescence machinery, the dead-leaf dry-matter transfer read from
      ! SHARED: when lamina dry matter senesces, the senescing tissue's minerals
      ! leave the leaf pool, the element-specific PHLOEM EXPORT FRACTION of
      ! driver column 12 goes to the grain and the remainder stays in dead
      ! tissue. Biological basis: silencing NAM abolishes mineral remobilization
      ! together with that of nitrogen (Waters et al. 2009, J Exp Bot 60:4263,
      ! DOI 10.1093/jxb/erp257, Table 1). The export fractions themselves are
      ! pinned in the driver with their sources - Waters et al. 2009 for Zn and
      ! Fe, Garnett & Graham 2005 for Cu and for the zero of Mn (their Table 3
      ! and Fig. 8 show no net manganese remobilization), and zero for calcium,
      ! which is phloem-immobile (White & Broadley 2003).
      export = 0.0_wp
      IF (d_dead > 0.0_wp .AND. M_leaf(i) > 0.0_wp) THEN
        ! Fraction of the PRE-senescence leaf pool that senesced today; W_green
        ! read from SHARED is post-senescence, so the pre-senescence mass is
        ! green + dead.
        sfrac = d_dead / MAX(d_dead + MAX(W_green_sh, 0.0_wp), 1.0e-9_wp)
        sfrac = MIN(1.0_wp, MAX(0.0_wp, sfrac))
        moved     = M_leaf(i) * sfrac              ! senescing tissue's minerals
        export    = mm_export(i) * moved           ! phloem export fraction
        dead_part = moved - export                 ! stays in dead tissue
        IF (W_grain > 1.0e-9_wp) THEN
          M_leaf(i)  = M_leaf(i)  - moved
          M_grain(i) = M_grain(i) + export
          M_dead(i)  = M_dead(i)  + dead_part
        ELSE
          ! No grain sink yet: the phloem-export share stays in the LIVING leaf
          ! and only the non-export share follows the dead dry matter.
          M_leaf(i) = M_leaf(i) - dead_part
          M_dead(i) = M_dead(i) + dead_part
          export = 0.0_wp
        END IF
      END IF
      ! ---- 4(iii). Record 63: the SAME senescence-coupled remobilization law,
      !      extended from the lamina to the stem and the chaff on the same
      !      per-element phloem export fractions (driver column 12), with each
      !      organ's senesced fraction formed against its OWN pre-senescence mass
      !      (living plus dying) and the same no-grain-sink fallback. Without it
      !      the stem's and chaff's minerals cannot move at all, mineral
      !      remobilization being wired to leaf death alone. No new constant.
      !      Off: the two death fluxes are 0 and the block is a no-op. ------ !
      IF (grainN_partition == 1) THEN
        IF (d_dead_stem > 0.0_wp .AND. M_stem(i) > 0.0_wp) THEN
          sfrac = d_dead_stem / MAX(d_dead_stem + MAX(W_stem_sh, 0.0_wp), 1.0e-9_wp)
          sfrac = MIN(1.0_wp, MAX(0.0_wp, sfrac))
          moved     = M_stem(i) * sfrac
          dead_part = moved - mm_export(i) * moved
          IF (W_grain > 1.0e-9_wp) THEN
            M_stem(i)  = M_stem(i)  - moved
            M_grain(i) = M_grain(i) + (moved - dead_part)
            M_dead(i)  = M_dead(i)  + dead_part
            export = export + (moved - dead_part)
          ELSE
            M_stem(i) = M_stem(i) - dead_part
            M_dead(i) = M_dead(i) + dead_part
          END IF
        END IF
        IF (d_dead_chaff > 0.0_wp .AND. M_chaff(i) > 0.0_wp) THEN
          sfrac = d_dead_chaff / MAX(d_dead_chaff + MAX(W_chaff_sh, 0.0_wp), 1.0e-9_wp)
          sfrac = MIN(1.0_wp, MAX(0.0_wp, sfrac))
          moved     = M_chaff(i) * sfrac
          dead_part = moved - mm_export(i) * moved
          IF (W_grain > 1.0e-9_wp) THEN
            M_chaff(i) = M_chaff(i) - moved
            M_grain(i) = M_grain(i) + (moved - dead_part)
            M_dead(i)  = M_dead(i)  + dead_part
            export = export + (moved - dead_part)
          ELSE
            M_chaff(i) = M_chaff(i) - dead_part
            M_dead(i)  = M_dead(i)  + dead_part
          END IF
        END IF
      END IF
      M_remob_grain(i) = export
      M_cap_binding(i) = 0.0_wp                    ! 0 unless a grain bound cuts the loading below
      ! ---- Record 126: THE GRAIN MINERAL CONCENTRATION BOUND (off: block not
      !      entered) ------------------------------------------------------- !
      ! The day's grain loading is limited to m_max,e x the day's grain
      ! dry-matter increment. Excess returns to its source - the uptake share to
      ! the leaf terminus, the export share to dead tissue - so the total is
      ! unchanged and the closure below still holds. The driver ships
      ! m_max = 0.0 for every element, which is no bound.
      IF (mineral_grain_sink_cap == 1 .AND. mm_grain_max(i) > 0.0_wp) THEN
        g_load  = M_grain(i) - g_before                        ! mg m-2 loaded today (uptake share + export [+ xylem])
        g_allow = mm_grain_max(i) * 1.0e-3_wp * MAX(0.0_wp, dg_dm)   ! mg kg-1 x g m-2 / 1000 = mg m-2
        IF (g_load > g_allow + 1.0e-15_wp) THEN
          g_excess   = g_load - g_allow
          g_from_exp = MAX(0.0_wp, export)
          g_from_upt = MAX(0.0_wp, g_load - g_from_exp)
          f_upt = 0.0_wp
          IF (g_from_upt + g_from_exp > 0.0_wp) f_upt = g_from_upt / (g_from_upt + g_from_exp)
          M_grain(i) = M_grain(i) - g_excess
          M_leaf(i)  = M_leaf(i)  + g_excess * f_upt
          M_dead(i)  = M_dead(i)  + g_excess * (1.0_wp - f_upt)
          M_remob_grain(i) = MAX(0.0_wp, export - g_excess * (1.0_wp - f_upt))
          M_cap_refused(i) = M_cap_refused(i) + g_excess
          M_cap_days(i)    = M_cap_days(i) + 1
          M_cap_binding(i) = 1.0_wp
        END IF
      END IF


      ! ---- Record 171: THE GRAIN PROTEIN SINK FOR PHLOEM-MOBILE MINERALS
      !      (off: block not entered) --------------------------------------- !
      ! THE LAW, as printed. Xia et al. 2020, Food Energy Secur 9:e243,
      ! DOI 10.1002/fes3.243, printed p.7: "Analysis of the speciation and
      ! localization of Zn in wheat grains indicates that Zn interacts with
      ! proteins, and therefore, grain protein content (as a physiological sink
      ! for Zn) may determine the sink capacity"; the same paper's printed
      ! pp.6-7 give the SHAPE, "zero or relatively low N supply results in low
      ! wheat grain Zn concentration, whereas optimal N supply achieves high
      ! grain Zn concentration, but excess N supply does not further increase
      ! grain Zn, grain N, or protein concentrations". Corroborated at first hand
      ! by Cakmak 2008, Plant Soil 302:1, printed p.6, and by Waters et al. 2009,
      ! J Exp Bot 60:4263, pp.4267-4268, which is already the cited basis of the
      ! senescence-export law in blocks 4(ii) and 4(iii) above.
      ! AS CODED. The day's grain loading of element e - the phloem share of
      ! current uptake plus the senescence export share - is bounded by
      !   allow = r_e x dN_grain(today)
      ! at the per-element stoichiometry r_e of driver column 20 (mg element per
      ! g grain nitrogen). Refused loading RETURNS TO ITS SOURCE in the
      ! proportions in which it arrived, the current-uptake share to the leaf
      ! terminus and the senescence share to dead tissue, which is exactly the
      ! routing record 126 uses, so the per-element daily balance is untouched:
      ! the block moves a flux between pools and creates or destroys nothing.
      ! DECLARED CONSTRUCTIONS, named as such because none is printed anywhere:
      !   (a) PHLOEM IMPORT IS IRREVERSIBLE. The bound limits the day's IMPORT;
      !       mineral already deposited in the grain is never withdrawn, so on a
      !       day when the ceiling has fallen below the standing pool the
      !       admissible loading is zero, not negative.
      !   (b) THE REFUSED SHARE'S DESTINATION is the pair of termini record 126
      !       declares, the model carrying no phloem pool in which refused
      !       mineral could wait.
      !   (c) r_e IS CONSTANT OVER THE FILLING PERIOD. Xia's sentence is about
      !       the mature grain's protein and no printed time course exists.
      ! The allowance is formed on the FLUX convention, against the day's grain
      ! nitrogen increment. An element whose r_e is 0.0 is skipped entirely, so
      ! the record is a provable no-op for it.
      IF (mineral_n_coupling_basis >= 1 .AND. mm_grain_eN(i) > 0.0_wp &
          .AND. W_grain > 1.0e-9_wp) THEN
        g_load = M_grain(i) - g_before                       ! mg m-2 loaded today
        IF (mineral_n_coupling_basis == 1) THEN
          mnc_allow = mm_grain_eN(i) * mnc_dngrain           ! mg (gN)-1 x g N m-2 = mg m-2
        ELSE
          ! RETIRED-OPTION REJECTION (v1.0): record 171 mode 2 had its own allowance
          ! formula removed by the v1.0 clean-up, which left mnc_allow undefined here.
          ! The startup validation in main_MMCO2Crop.f90 refuses record 171 /= 1, so
          ! this arm is unreachable in a run that started; it exists so that no path
          ! can read mnc_allow undefined.
          mnc_allow = 0.0_wp
          WRITE (*,'(A,I0)') ' REFUSING TO RUN: mineral loading reached with record 171 = ', &
            mineral_n_coupling_basis
          WRITE (*,'(A)') '   Only record 171 = 1 is implemented in M&M-CO2Wheat v1.0.'
          ERROR STOP 'rec 171: retired mineral_n_coupling_basis'
        END IF
        IF (g_load > mnc_allow + 1.0e-15_wp) THEN
          mnc_excess = g_load - mnc_allow
          g_from_exp = MAX(0.0_wp, M_remob_grain(i))         ! senescence share loaded today
          g_from_upt = MAX(0.0_wp, g_load - g_from_exp)      ! uptake (+ xylem) share
          f_upt = 0.0_wp
          IF (g_from_upt + g_from_exp > 0.0_wp) f_upt = g_from_upt / (g_from_upt + g_from_exp)
          M_grain(i) = M_grain(i) - mnc_excess
          M_leaf(i)  = M_leaf(i)  + mnc_excess * f_upt
          M_dead(i)  = M_dead(i)  + mnc_excess * (1.0_wp - f_upt)
          M_remob_grain(i) = MAX(0.0_wp, M_remob_grain(i) - mnc_excess * (1.0_wp - f_upt))
          mnc_refused(i) = mnc_refused(i) + mnc_excess
          mnc_days(i)    = mnc_days(i) + 1
          M_cap_binding(i) = 1.0_wp
        END IF
        mnc_to_grain(i) = mnc_to_grain(i) + MAX(0.0_wp, M_grain(i) - g_before)
      END IF

      ! ---- 5. Element mass balance: only uptake adds atoms, senescence and
      !      grain loading being internal transfers ---------------------- !
      after = M_root(i) + M_leaf(i) + M_stem(i) + M_chaff(i) + M_grain(i) + M_dead(i)
      M_balance_resid(i) = after - before - uptake
    END DO
    ! Record 171: carry today's grain nitrogen forward as this module's own
    ! saved previous value; nothing in SHARED is written.
    mnc_ngrain_prev = mnc_ngrain
  END SUBROUTINE update_minerals_daily

END MODULE MINERAL_MODULE