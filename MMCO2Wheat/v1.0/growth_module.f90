!===============================================================================!
! M&M-CO2Wheat v1.0  |  growth_module.f90
!
! Purpose: Crop carbon and nitrogen economy on a DAILY step. It partitions the
!          day's assimilate between lamina, true stem, root and spike, charges
!          maintenance and growth respiration, builds and senesces leaf area,
!          deposits and draws the stem water-soluble-carbohydrate reserve, pays
!          every structural increment in reduced nitrogen, turns organ nitrogen
!          over after flowering and loads the grain, sets grain number and
!          kernel mass, and descends the root front. Spatial units: the plant
!          organ and, where the organ-level canopy is active, the tiller axis,
!          the leaf rank and the individual blade or sheath; all pools are per
!          square meter of ground.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
!
! Uses: TYPES, ONLY: wp (working precision); SHARED, ONLY: <explicit list> -
!       the measured crop parameters (specific leaf area and its development
!       profile, fruiting efficiency, potential kernel mass, chaff share, the
!       partition trait nodes, the reserve and remobilization constants, the
!       structural and capacity nitrogen concentrations, the Abichou organ
!       geometry and the Fischer spike constants), every driver-record switch
!       this module obeys, and the ledgers it writes.
! Output: the organ dry-matter pools W_green, W_dead, W_stem, W_root,
!       W_storage, W_grain, W_chaff and the reserve W_reserve; leaf area index
!       and, under the organ law, the per-rank and per-axis green and dead
!       areas; canopy height; root depth; the organ nitrogen pools N_leaf,
!       N_stem, N_storage, N_chaff and N_dead_leaf and the mobile vegetative
!       pool N_mobile_veg; grain number and kernel mass; and the daily carbon
!       and nitrogen ledgers and closure residuals. Leaf area and leaf nitrogen
!       are consumed by CANOPY_MODEL and assim (the canopy nitrogen profile and
!       gas exchange), the root depth and the transpiration demand by hydro,
!       the organ dry-matter increments by mineral_module (the mineral
!       allocation shares), and every pool and ledger by main_MMCO2Crop, which
!       writes them out.
! Called by: main_MMCO2Crop.f90 only (USE GROWTH_MODULE, ONLY: ...), once per
!       day in the daily loop, in the order initialize_growth_state and the two
!       bootstrap routines at the start of the run, then per day
!       allocate_daily_biomass, update_structure, the organ-canopy routines
!       (update_canopy_leaf_ranks, update_canopy_blade_area,
!       update_canopy_tillering, update_leaf_expansion_water,
!       update_canopy_leaf_senesc, update_canopy_sheath_ear,
!       update_canopy_sink_state), apply_leaf_senescence,
!       allocate_yield_components_daily, apply_organN_turnover_daily,
!       allocate_grain_N_source_sink_daily and update_root_depth.

!===============================================================================!
MODULE GROWTH_MODULE
  USE TYPES, ONLY: wp
  USE SHARED, ONLY: SLN_min_N, k_remob_leaf_N, k_remob_stem_N, n_struct_stem_N, &
                    k_N_can, n_b_can, n_stru_grain, dvs_grain_stru_end, &
                    gParam_sla, sla_dvs0, sla_dvs2, sla_anchor, &
                    fruiting_efficiency, TKW_potential, chaff_share, ear_end_dvs2, &
                    storage_onset_dvs2, storage_peak_frac, stem_spike_realloc_delta, &
                    fleaf_elong, fleaf_heading, fleaf_end_dvs, &
                    carbon_reserve_mode, stem_reserve_frac, reserve_remob_rate, &
                    stem_reserve_deposition_n_basis, vhw_a_mgg, vhw_b_mgg, &
                    vhw_nni_last, vhw_frac_last, vhw_frac_min, vhw_frac_max, &
                    vhw_days_live, vhw_days_clamp_hi, vhw_days_clamp_lo, &
                    vhw_w_tha_in, vhw_n_sh_in, &
                    reserve_start_dvs2, reserve_end_dvs2, reserve_conv_eff, W_reserve, &
                    stem_wsc_mech, wsc_cap_frac, wsc_accum_end_dpa, wsc_elong_dvs2, &
                    wsc_reserve_law, wsc_law_member, wsl_ddvs2, wsl_dtau_day, wsl_tau_rem_day, &   ! rec 105/106
                    wsl_room_day, wsl_draw_def_day, wsl_cap_day, wsl_draw_day, wsl_bound_day, &
                    wsl_days_ceil, wsl_days_def, wsl_days_pool, wsl_days_stem, wsl_frac_max, &
                    wsc_rel_window_member, wsl2_tau_day, wsl2_taurem_day, wsl2_wres_pre, &   ! rec 105 = 2 / rec 107
                    wsl2_rel_day, wsl2_relg_day, wsl2_relf_day, wsl2_gain_day, wsl2_lconv_day, &
                    wsl2_cres_day, wsl2_bound_day, wsl2_rel_tot, wsl2_relg_tot, wsl2_relf_tot, &
                    wsl2_gain_tot, wsl2_loss_conv_tot, wsl2_loss_fate_tot, wsl2_days_zerodef, &
                    wsl2_i1_max, wsl2_i2_max, wsl2_i3_max, &
                    chaff_growth_mode, cgs_frac, cgs_W_anth, cgs_pot_tot, cgs_win, &   ! rec 117 POST-ANTHESIS CHAFF SINK
                    cgs_ddvs2, cgs_demand_day, cgs_taken_day, cgs_unmet_day, &
                    cgs_demand_tot, cgs_taken_tot, cgs_unmet_tot, cgs_from_repro_tot, &
                    cgs_repro_seen_tot, cgs_days_win, cgs_days_denied, cgs_day_anth, cgs_anch_set, &
                    cgs_c6_resid, cgs_c6_max, cgs_state_max, cgs_post_today, &
                    wsc_day_anth, wsc_pool_anth, wsc_stem_anth, wsc_pool_peak, wsc_day_peak, &
                    wsc_release_delay_basis, wrd_days_delayed, wrd_room_held_cum, &   ! rec 165
                    wrd_pool_held_max, &
                    wsc_dep_pre, wsc_dep_post, wsc_refused_cap, wsc_refused_win, &
                    wsc_draw_tot, wsc_draw_after16, wsc_gain_tot, wsc_resp_tot, &
                    wsc_day_mob0, wsc_days_capbnd, wsc_days_winbnd, wsc_closure_max, &
                    wsc_dep_day, &
                    n3_fill_snap, n3_dvs0, n3_leaf0, n3_stem0, n3_cum_leaf, n3_cum_stem, &
                    n3_eff_leaf, n3_eff_stem, &
                    leafN_drawdown, k_leafN_TT, sln_max_N, dTT_day_N, &
                    leafN_cap_route, lnc_lai_basis, lnc_excess_day, lnc_to_mobile_day, &   ! rec 87
                    itemF_Nclose_mode, fnc_syn_clip_day, fnc_syn_clip_cum, &            ! rec 97
                    leafn_hold_law, sln_max_N_lnh, lnh_cap_in_force, &                  ! rec 132
                    lnh_syn_assim_off_day, lnh_syn_mob_off_day, &
                    lnh_redirect_day, lnh_redirect_cum, lnh_binding_day, &
                    leafn_turnover_law, lto_cap_in_force, &                             ! rec 133
                    lto_sen_N_out_day, lto_sen_N_resorb_day, lto_s_lam_pot_day, &
                    lto_grain_claim_day, lto_grain_served_assim_day, lto_grain_served_mobile_day, &
                    lto_stem_chaff_draw_day, lto_mob_before_day, lto_mob_after_day, &
                    lto_leftover_day, lto_syn_realized_day, lto_syn_realized_cum, &
                    lto_syn_from_assim_day, lto_syn_from_mobile_day, &
                    lto_refused_day, lto_refused_cum, lto_binding_day, &
                    spike_sink_mech, bgn_r_spike, bgn_m_floret_floral, &                ! rec 98
                    bgn_w_struct_axis, bgn_demand_day, bgn_axes_day, bgn_dem_cum, &
                    bgn_got_cum, bgn_days_dem, bgn_days_room, bgn_days_win, &
                    bgn_floral_bound, bgn_ear_spike_dm, bgn_fe_ear_basis, &
                    bgn_dm_close_max, bgn_agb_day, bgn_spi_anth, &
                    bgn_spike_dw_ref, bgn_f_struct, bgn_f_floral, bgn_m_struct_spklt, &
                    fnc_mob_close_day, fnc_mob_close_cum, &                             ! rec 97
                    n_stem_resid_conc, nJ_R0, nJ_cum, &
                    N_mobile_veg, N_from_mobile_day, &
                    grain_fill_ended, canopy_senesced, assimN_surplus_postfill, &
                    N_chaff, chaff_share, &
                    grainN_partition, n4_nhi_target, n4_ag_anth, &
                    n2_delta_TT, n2_sigma_lam, n2_sigma_stem, &
                    n2_sigma_chaff, n2_k1_Nmob, n2_k2_PAR, n2_struct_lam, &
                    n2_struct_stem, n2_d_tp, n2_par_frac, dTotrad_s, &
                    W_dead_stem, N_dead_stem, n2_anth_snap, &
                    n2_cph_stem_anth, n2_cph_chaff_anth, &
                    n2_deg_lam_day, n2_deg_stem_day, n2_deg_chaff_day, &
                    n2_syn_lam_day, n2_syn_stem_day, n2_syn_chaff_day, &
                    n2_dW_dead_stem_day, n2_dW_dead_chaff_day, &
                    n2_cum_deg_lam, n2_cum_deg_stem, n2_cum_deg_chaff, &
                    unalloc_C_total, overflow_to_reserve_total, &
                    cres_to_grain_assim, cres_to_grain_res, cres_sink_demand, &
                    cres_sink_unfilled, cres_respired, cres_ledger_resid, &
                    CLR_MAXRANK, canopy_leaf_ranks, clr_phyllochron_Cd, &
                    clr_seed_ratio, clr_fC_exp, clr_fC_HS_ref, clr_LED_Cd, &
                    clr_P_evers_low, clr_P_evers_high, clr_evers_break_HS, &
                    clr_active, clr_frozen, clr_tau_Cd, clr_HS, clr_HS_ev, &
                    clr_tau_app, clr_tau_app_ev, clr_day_app, clr_day_lig, &
                    clr_n_app, clr_n_lig, clr_n_app_ev, clr_n_lig_ev, &
                    clr_FLN, clr_FLN_tip, clr_FLN_ev, clr_FLN_ev_tip, &
                    clr_tau_boot, clr_tau_app_flag, clr_tau_lig_flag, &
                    clr_tau_lig, clr_tau_lig_ev, &
                    canopy_blade_area, cba_l_irank, &
                    cba_Nl_fit, cba_w_irank, cba_ft_r1, cba_ft_r4, &
                    cba_ft_adult, cba_ft_flag, cba_ft_n4, cba_ft_n8, &
                    cba_sTi, cba_lb1, cba_kw, cba_kl_st, cba_kw_st, &
                    cba_L_pot, cba_W_pot, cba_ft_rank, cba_A_pot, cba_A, &
                    cba_A_shoot_pot, cba_A_shoot, cba_n_expanding, &
                    CT_MAXTIL, canopy_tillering, ct_hs_emerg, ct_pint_p1, &
                    ct_pint_p0, ct_kL, ct_axmax, ct_die_frac, ct_hs_sm_lo, &
                    ct_hs_sm_hi, ct_nl_sm_lo, ct_nl_sm_hi, ct_theta, ct_delta, &
                    ct_nl_alpha, ct_nl_beta, ct_nc_offset, ct_internodes_elong, &
                    ct_active, ct_ceased, ct_cess_why, ct_fpar, ct_dens, Axes_max_per_plant, &
                    ct_dens_born, ct_p_occ, ct_pint_win, ct_win_sum, ct_win_n, &
                    ct_decided, ct_appeared, ct_tau_birth, ct_nl_axis, &
                    ct_area_axis, ct_area_axis_real, ct_area_axis_st, &
                    ct_axes_live, ct_axes_born, ct_axes_max, ct_ears, &
                    tiller_survival_mech, tsv_sup, tsv_dem, tsv_kill, tsv_ilight, &
                    tsv_apot_prev, tsv_dens_pre, tsv_taun, tsv_days, tsv_frozen, &
                    tsv_cullgreen_seen, tsv_cullgreen_charged, &
                    tca_gsup_day, tca_gsup_prev, tca_S, tca_D, tca_phi, tca_share, &
                    tca_ftot, tca_dens_pre, tca_area_prev, tca_dead, tca_unaff, tca_dead_c, tca_days, &
                    tca_days_skip, tca_frozen, tca_hs, tca_L, tca_sax, tca_dax, &
                    tca_ltot, tca_sum_s, tca_ident, tca_taun, &
                    tca4_theta, tca4_delta, tca4_U, tca4_F, tca4_target, tca4_endpoint_day, &
                    stem_demand_mech, sdm_dens, sdm_h_pot, sdm_h_pot_prev, &
                    sdm_dvs2_flag, sdm_f_ped, sdm_h_collar, sdm_h_real, &
                    sdm_Lped_pot, sdm_Lped_emerg, sdm_dh_day, sdm_frac_real, &
                    sdm_rm_stem_day, sdm_pot_dw_day, sdm_dem_day, sdm_avail_day, &
                    sdm_take_day, sdm_dw_stem_day, sdm_legacy_dw_stem_day, &
                    sdm_refused_day, sdm_refused_cum, sdm_dem_cum, sdm_got_cum, &
                    sdm_ident_day, sdm_ident_legacy_day, sdm_ident_max, &
                    sdm_ident_legacy_max, sdm_close_max, &
                    sdm_days_dem, sdm_days_sup, sdm_days_ref, &
                    sdm4_dur_Cd, sdm4_onset_k, sdm4_r_mg, sdm4_cap_mg, &   ! rec 100 MODE 4
                    sdm4_tau_on, sdm4_tau_se, sdm4_x_first, sdm4_pot_dw_today, &
                    sdm4_nl_ms, sdm4_n_open_ms, sdm4_n_open_wtd, sdm4_dem_ms_mg, &
                    sdm4_cap_bound_today, sdm4_days_cap, sdm4_t_close_last, sdm4_t_close_subtop, &
                    sdm4_t_close_all, &
                    surplus_fate_mech, sfl_g_wsc, sfl_dw_wsc_day, sfl_resid_day, &   ! rec 101
                    mnp_dw_leaf_day, mnp_dw_root_day, mnp_dw_stem_day, &                 ! rec 114
                    mnp_dw_chaff_day, mnp_dw_grain_day, &
                    sfl_want_day, sfl_room_day, sfl_stored_dm_cum, sfl_stored_ch2o_cum, &
                    sfl_resid_cum, sfl_resid_pre, sfl_resid_post, sfl_pool_ratio_max, &
                    sfl_ledger_day, sfl_ledger_max, sfl_days_cap, sfl_days_win, sfl_days_store, &
                    ct_tau_sm, ct_tau_em, ct_hs_sm, ct_bai_pot, ct_bai, &
                    ct_fs_implied, &
                    canopy_leaf_senesc, cls_g1, cls_g2, cls_g3, cls_alpha, &
                    cls_beta, cls_delta, cls_top_offset, cls_active, cls_gl, &
                    cls_sl, cls_phase, cls_tau1, cls_tau2, cls_tau3, &
                    cls_area_green, cls_area_dead, cls_gai, cls_dai, &
                    cls_dai_max, cls_mono_viol, cls_dai_culled, cls_dens_prev, &
                    cls_dai_pool, cls_dai_dip_max, cls_dai_prev, &
                    cls_gl_at_flag, cls_hs_onset, cls_gai_flag, &
                    canopy_n_flag_zone, kL_bertheloot, canE_fPAR_top, canE_fPAR_eff, &
                    canopy_sheath_ear, cse_ls_irank, cse_li_irank, cse_ds1, cse_ds2, &
                    cse_di2, cse_rho1, cse_rho2, cse_rho3, cse_persist_n, &
                    cse_a_ear_meas, cse_eai_meas, &
                    cse_adel_start, cse_adel_dur, cse_a_ep, cse_a_ep_lo, cse_a_ep_hi, &
                    cse_ep_win, cse_half_cyl, cse_active, cse_ls_pot, cse_ds_pot, &
                    cse_sh_green, cse_sh_dead, cse_ep_axis, cse_st_axis, &
                    cse_sai, cse_sai_n1, cse_sai_n3, cse_sai_full, cse_sdi, &
                    cse_sdi_pool, cse_sdi_prev, cse_sdi_dip_max, cse_mono_viol, &
                    cse_eai, cse_eai_lo, cse_eai_hi, cse_iai, cse_nonlam, cse_ratio, &
                    cse_paig, cse_pai_li, cse_paig_axis, cse_paig_flower, &
                    cse_axis_flower, cse_ratio_flower, cse_pai_li_max, &
                    canopy_sink_lai, csl_LMA, CSL_FLN_CUT_MAX, csl_active, &
                    csl_limit_active, csl_lai_sink, csl_pai_diag, csl_area_built, &
                    csl_area_prev, csl_dA_day, csl_demand_day, csl_alloc_day, &
                    csl_backlog, csl_backlog_max, csl_limited_days, csl_fln_cut, &
                    n_payment_gate, npg_floor_L, npg_floor_R, npg_build_R, npg_k_lux, &   ! rec 121
                    npg_n_assim_root_in, npg_n_assim_shoot_in, npg_pool_root_in, npg_root_export, &
                    npg_paid_R_root, npg_paid_R_shoot, npg_paid_cum_R_root, npg_paid_cum_R_shoot, &
                    npg_live_today, npg_c_L, npg_c_S, npg_c_R, npg_pool_mobile_in, npg_pool_assim_in, &
                    npg_lux_avail, npg_lux_used, npg_D_L, npg_D_S, npg_D_R, npg_p_L, npg_p_S, npg_p_R, &
                    npg_paid_L, npg_paid_S, npg_paid_chaff, npg_paid_R, npg_pool_out, npg_mobile_left, &
                    npg_dw_leaf_pot, npg_dw_stem_pot, npg_dw_stor_pot, npg_dw_root_pot, &
                    npg_dw_leaf, npg_dw_stem, npg_dw_stor, npg_dw_root, &
                    npg_unpaid_ch2o_day, npg_dw_wsc_day, npg_resid_day, npg_lamina_unpaid_day, npg_lamina_unpaid_cum, &
                    npg_unpaid_ch2o_cum, npg_wsc_cum, npg_resid_cum, npg_paid_cum_L, npg_paid_cum_S, npg_paid_cum_R, &
                    npg_days_live, npg_days_bind_L, npg_days_bind_S, npg_days_bind_R, npg_ident_N_max, &
                    npg_N_leaf => N_leaf, npg_N_stem => N_stem, npg_N_root => N_root, npg_N_storage => N_storage, &
                    n4_sln_target, n4_stem_nmax, n5_stem_nstruct, &
                    lai_n_arrest, lna_f, lna_fprov, lna_fclk, lna_dinc, lna_apot, lna_live_today, &   ! rec 122
                    lna_pL_day, lna_inc_area_day, lna_removed_day, lna_removed_cum, lna_removed_dm_day, &
                    lna_removed_dm_cum, lna_gai_before, lna_gai_after, lna_sl_bump, lna_sln_green_in, &
                    lna_shed_need_day, lna_shed_day, lna_shed_cum, lna_days_arrest, lna_days_shed, &
                    lna_shed_unmet_days, lna_expdead_rank_days, SLN_min_N, &
                    tiller_n_afford, tna_DN, tna_SN, tna_phiN, tna_phiC, tna_pL_used, tna_dN_day, &   ! rec 123
                    tna_days_bind, tna_days_step, tna_step_today, &
                    n_afford_priority_law, nap_eps_gate, nap_eps_gate_t3, nap_m_min, &   ! rec 134
                    nap_nni_ctx_lo, nap_nni_ctx_hi, nap_num, nap_den, nap_phi_realized, &
                    nap_nni_veg, nap_m, nap_gate12_open, nap_gate3_open, nap_tier3_cap, &
                    nap_f_Nstatus_day, nap_f_Nstatus_at_gate, nap_N_mobile_veg_at_gate, &
                    nap_identity_resid, nap_days_live, &
                    nap_phi_standing_testday, nap_m_testday, nap_gate3_open_testday, &   ! rec 134 stale-cap fix
                    nap_tier3_cap_testday, nap_tier3_cap_declared_inert, &
                    root_front_law, rfl_z0, rfl_rz_mode, rfl_rz, rfl_days_emerged, &   ! rec 124
                    rfl_z_legacy, rfl_z_new, rfl_z_max, &
                    kernel_potential_law, kpt_env_min, kpt_env_max, kpt_tt, kpt_sup_cum, kpt_tsum, &   ! rec 125
                    repro_surplus_fate, rsf_resid_day, rsf_resid_cum, rsf_label_cum, rsf_surplus_cum, rsf_days, &   ! rec 127
                    sink_full_carbon_fate, cfo_surplus_day, cfo_wsc_day, cfo_wsc_cum, &   ! rec 129
                    cfo_loss_day, cfo_loss_cum, cfo_days_loss, cfo_days_store, &
                    grain_n_min_conc, gnm_min_conc, gnm_N_grain_in, gnm_dN_grain_prev, gnm_guard_dm, &   ! rec 130
                    gnm_room_dm_118, gnm_room_dm_125, gnm_room_dm_floor, gnm_to_grain_before, &
                    gnm_to_grain_after, gnm_refused_dm_day, gnm_refused_dm_cum, gnm_floor_bind_day, gnm_days_bind, &
                    gnm_kcd, gnm_Nstru0, gnm_Der, gnm_Tt_degCd, gnm_Nstru_pergrain, gnm_phase_day, &   ! rec 130 = 3
                    gnm_Nstru_total, gnm_Nstru_total_prev, gnm_dNstru_day, gnm_c_stru_day, &
                    gnm_Nstru_paid, gnm_Nsto, &
                    kpt_days, kpt_day_set, kpt_bound, kpt_set, kpt_sk_mean, kpt_tmean, kpt_dfill_days, &
                    kpt_w_yin, kpt_w_raw, kpt_w_pot, gnc_cap, gnc_room_day, gnc_cut_day, gnc_cut_cum, &
                    gnc_days_bound, gnf_mgN_kernel_hi, gnp_mgN_kernel, &   ! rec 177
                    csl_cohorts_aborted, csl_births_suspended, csl_sla_real, &
                    csl_lai_Ncap, csl_fsum_err, &
                    tiller_carbon_gate, leaf_expansion_water, lai_from_organs, &
                    csl_dm_gain_day, tcg_csur_day, tcg_sup, tcg_dem, tcg_ratio, tcg_nwin, &
                    tcg_nrank, tcg_env600_days, tcg_env800_days, tcg_axes_env_lo, &
                    tcg_axes_env_hi, tcg_die_frac_lo, tcg_die_frac_hi, tcg_ne_env_lo, &
                    tcg_ne_env_hi, lew_a200, lew_b200, lew_a800, lew_b800, lew_srwc, &
                    lew_se_rw, lew_theta_rw, lew_thetafc_rw, lew_fw, lew_fw_200, &
                    lew_fw_800, lew_f, lew_fclk, lew_area_axis, lew_bai, lew_bai_clock, &
                    lfo_sla_out, hydro_root_dist_diag, theta, theta_33, theta_r, theta_0, &
                    cls_sl_prev, lfo_gai_died_day, lfo_fdied_day, lfo_sen_dm_day, lfo_sen_N_out_day, &
                    lfo_sen_N_dead_day, lfo_sen_N_resorb_day, lfo_sen_dm_cum, lfo_sen_N_out_cum, &
                    lfo_sen_N_dead_cum, lfo_sen_N_resorb_cum, &
                    grain_number_mech, gnx_spike_part, gnx_m_floret, gnx_kern_per_fl, &
                    gnx_floret_scale, &                                               ! rec 151
                    gnx_tt_base_C, gnx_win_open, gnx_win_open_day, gnx_win_days, &
                    gnx_win_tt45, gnx_win_agb0, gnx_win_spike0, gnx_win_dagb, &
                    gnx_win_dspike, gnx_share_real, gnx_clamp_days, gnx_clamp_max, &
                    gnx_pool, gnx_dens, gnx_wsp, gnx_nfl, gnx_kn, gnx_is_ear, &
                    gnx_unshared, gnx_ears_m2, gnx_florets_m2, gnx_grains_per_ear, &
                    gnx_spike_dw_anth, gnx_sterile_C, gnx_set, gnx_pool_scale, &
                    gnx_dagb_day, gnx_dspike_day, gnx_sum_tdu, gnx_sum_cd45, &
                    gnx_tt_ratio, gnx_tdu_rem_open, gnx_penult_day, gnx_tt_ref_Cd, &
                    gnx_day_now, &
                    floret_survival_mech, gnx_spikelets, gnx_mfs, gnx_ninit, &
                    gnx_surv, gnx_pool_pre, gnx_axes_culled, gnx_pool_moved, &
                    gnx_cull_passes, gnx_fe_realised, &
                    gnx_phi_sp, gnx_wpot_axis, gnx_earless_m2, &                   ! rec 90 = 4
                    kernel_number_m2, gnf_ng0_basis, gnk_tier, gnk_cohort, &      ! rec 115 PER-KERNEL GRAIN-N SINK BASIS
                    gnk_fb_events, gnk_fb_day, gnk_divisor, gnk_numer, &
                    gnk_ears_snap, gnk_ng0_perear, gnk_factor, &
                    grain_fill_mech, gfb_cm_ref, gfb_cm_slope, gfb_t_ref, &   ! (rec 91), SHARED.f90-
                    kernel_mass_ceiling, &                                      ! rec 118
                    gfb_tbase, gfb_dfill, gfb_u, gfb_shape_int, &
                    gfb_tt, gfb_day0, gfb_day_end, gfb_days, gfb_end_bound, &
                    gfb_tsum, gfb_demand_tot, gfb_supplied_tot, gfb_unmet_tot, gfb_reserve_end, &
                    gfb_rate_max, gfb_cm_min, gfb_cm_max, gfb_cm_clamped, &
                    grainN_sink_mech, gnf_ng0_per_ear, gnf_g_rate, gnf_tau_degCd, &   ! (rec 86), SHARED.f90-
                    gnf_snap, gnf_ng0, gnf_tt, gnf_pot_cum, gnf_sink_day, &
                    gnf_supply_day, gnf_draw_lam, gnf_draw_stem, gnf_draw_chaff, &
                    gnf_cum_lam, gnf_cum_stem, gnf_cum_chaff, gnf_cum_assim, &
                    gnf_cum_mob, gnf_days_fill, gnf_days_bind, gnf_first_bind_dpa, &
                    gnf_last_bind_dpa, gnf_bind_after_15d, gnf_nres_max, gnf_nres_day, &
                    gnf_sigma_lam, gnf_r_mob_ph, gnf_eqn6, gnf_deg_lam_day, &      ! Eqn-6 repair, SHARED.f90
                    gnf_syn_lam_day, gnf_syn_pot_day, gnf_fN_day, gnf_fPAR_day, &
                    gnf_pool_scale, gnf_cum_syn_lam, gnf_cum_deg_lam, &
                    height_internode_law, hgt_int_window_member, hgt_int_length_member, &   ! rec 109-113
                    hgt_int_lineic_member, hgt_int_ear_member, INL_MAXIR, INL_ONSET_SHIFT, &
                    inl_Lf_m1, inl_Lf_m2, inl_Lf_m3, inl_Lf_m4, inl_Lf_m5, &
                    inl_lam_m1, inl_lam_m2, inl_lam_m3, inl_lam_m4, inl_dI_cd, inl_Lear_kirby_cm, &
                    inl_ready, inl_ped_open, inl_L, inl_m, inl_dLpot, inl_dLact, &
                    inl_Lf_use, inl_lam_use, inl_tau_on, inl_tau_end, inl_absrank, inl_bind, &
                    inl_h_b1, inl_h_b2, inl_h_b3, inl_dM_struct, inl_dM_offer, inl_dM_len, &
                    inl_e8_resid, inl_e8_resid_max, inl_e8_over_max, inl_days_mass, inl_days_pot, &
                    inl_days_thicken, inl_open_days, inl_Mstruct_prev, inl_Mstruct_max, &
                    inl_axes, inl_Lear_cm, inl_dI_use, inl_onset_arm, &
                    n_construction_status, ncs_c_floor_L_eff, ncs_floor3_SLN_martre, &        ! rec 128 modes 2/3
                   n_construction_ceiling_curve, ncc_a_used, ncc_b_used, ncc_days_live, &   ! rec 156
                    ncs_W_t_ha, ncs_W_eff, ncs_c_crit, ncs_D_crit, &   ! rec 128
                    ncs_c_target_L, ncs_c_target_S, ncs_c_target_R, ncs_D_target_L, ncs_D_target_S, &
                    ncs_f_crit, ncs_c_target_eff_L, ncs_c_target_eff_S, ncs_c_target_eff_R, &
                    ncs_D2_L, ncs_D2_S, ncs_D2_R, &
                    ncs_paid1_L, ncs_paid1_S, ncs_paid1_chaff, ncs_paid1_R, &
                    ncs_paid2_L, ncs_paid2_S, ncs_paid2_chaff, ncs_paid2_R, &
                    ncs_c_actual_L, ncs_c_actual_S, ncs_c_actual_R, &
                    ncs_lamina_below_target_cum, ncs_days_fcrit, &
                    n_payment_luxury_pass2, ncs_lux2_cap_day, ncs_lux2_used_day, &   ! rec 135
                    ncs_lux2_used_cum, ncs_lux2_days, &
                    n_afford_standing_basis, nap_num_st, nap_den_st, nap_phi_standing, &   ! rec 136
                    nap_phi_elig12, nap_phi_elig3, &
                    floret_n_law, fnl_num, fnl_den, fnl_f_N, fnl_r, fnl_surv_c_only, &   ! rec 137
                    fnl_step_today, &
                    spike_n_priority, npg_c_chaff_mode2, npg_c_chaff, npg_D_chaff, npg_D_stem, &   ! rec 138
                    npg_p_chaff, npg_p_stem, ncs_D2_chaff, ncs_D2_stem, &
                    npg_days_bind_chaff, npg_days_bind_stem, snp_viol_priority, &
                    snp_win_dw_stor_pot_sum, snp_win_p_chaff_sum, snp_win_p_stem_sum, &
                    snp_win_days, snp_win_days_pchaff_lt1, snp_win_days_pstem_lt1, &
                    lamina_n_remob_pre, lnr_excess_day, lnr_rem_avail_day, lnr_rem_used_day, &   ! rec 141
                    lnr_rem_used_cum, lnr_days_rem_used, lnr_win_rem_used_sum, &
                    n_status_construction_price, nsc_status_day, &   ! rec 143
                    nsc_price_L, nsc_price_S, nsc_price_chaff, nsc_price_R, &
                    nsc_cut_L_day, nsc_cut_S_day, nsc_cut_R_day, &
                    nsc_cut_L_cum, nsc_cut_S_cum, nsc_cut_R_cum, &
                    nsc_nni_veg_sum, nsc_nni_veg_days, &
                    carbon_coupled_leaf_expansion, cbc_fmax, cbc_f, cbc_fclk, &   ! rec 144
                    cbc_surplus_day, cbc_area_bonus_day, cbc_area_bonus_cum, &
                    cbc_f_mean_expanding, cbc_f_max_reached, cbc_laimax, &
                    cbc_days_surplus, cbc_live_today, &
                    leaf_maintenance_basis, &   ! rec 145
                    lamina_n_area_retention, lar_sln_low, lar_nb_leaf, kN_canopy, &   ! rec 147
                    lar_sln_func, lar_N_lam, lar_a_sup, lar_need_day, &
                    lar_shed_day, lar_shed_cum, lar_gai_before, lar_gai_after, &
                    lar_sln_in, lar_sln_out, lar_N_released_cum, &
                    lar_days_shed, lar_days_unmet, &
                    lamina_expansion_n_limit, lnl_a_slope, lnl_b_max, lnl_c_icpt, &   ! rec 179
                    lnl_maxNpha, lnl_Nph_day, lnl_slnph_day, lnl_minNpha_day, &
                    lnl_a_sup_day, lnl_allow_day, lnl_fN_day, lnl_refused_day, &
                    lnl_refused_cum, lnl_days_refuse, lnl_live_today, &
                    lnl_gai_max_min, lnl_days_maxupd, lnl_gai_at_maxupd, &   ! rec 179 green-area guard on maxN^ph,a
                    wsc_fate_law, wfl_a_refused_day, wfl_a_refused_cum, wfl_a_days, &   ! rec 148
                    wfl_b_released_day, wfl_b_released_cum, wfl_b_days, &
                    wfl_gate_day, wfl_gate_bound, wfl_pool_at_gate, &
                    wfl_days_grain_fed, wfl_grain_fed_cum, &
                    wfl_draw_def_day, wfl_draw_def_cum, wfl_draw_def_days, &
                    organ_n_turnover_nema, n2_p_root, &                        ! rec 154
                    ont_W_root_in, ont_dw_root_day, ont_pool_before, ont_uptake_credit, &
                    ont_D_lam, ont_D_stem, ont_D_chaff, ont_D_root, &
                    ont_sink_pot, ont_grain_import, &
                    ont_S_lam, ont_S_stem, ont_S_chaff, ont_S_root, ont_S_lam_pot, &
                    ont_s_scale, ont_f_N, ont_f_PAR, ont_dTT, ont_mob_after, ont_resid, &
                    ont_bound_by, ont_live_today, &
                    ont_cum_D_lam, ont_cum_D_stem, ont_cum_D_chaff, ont_cum_D_root, &
                    ont_cum_S_lam, ont_cum_S_stem, ont_cum_S_chaff, ont_cum_S_root, &
                    ont_cum_grain, ont_cum_uptake, ont_resid_max, ont_cum_syn_clip, &
                    ont_days, ont_days_pool_bound, ont_days_scaled, ont_days_cap_clip, &
                    lamina_synthesis_sheath_basis, lsb_a_sheath, &                      ! rec 167
                    lsb_days_applied, lsb_days_no_area, lsb_msheath_cum, &
                    lsb_ratio_last, lsb_syn_lam_gain_cum, lsb_syn_stem_loss_cum, &
                    lsb_m_sheath_day, lsb_ratio_day, &                                   ! rec 167 daily, wiring
                    canopy_fln_prediction, clr_tau_boot_target, clr_FLN_pred, &          ! rec 174
                    cfp_HS_pred, cfp_FLN_realized, cfp_mismatch, &
                    grain_n_partition_basis, &                                           ! rec 175
                    gnp_g0_target, gnp_g0_deficit, gnp_g0_take, gnp_g0_refused, &
                    gnp_g0_done, gnp_m3_release_cum, gnp_m3_days, &
                    gnp_g0_paid_cum, gnp_g0_frac, n2_grain_anth_conc, &                    ! rec 175 mode 4,
                    n2_struct_chaff, n2_struct_chaff_use, &                                 ! rec 178/,
                    n2_struct_chaff_pask, cnr_nre_chaff, cnr_N_chaff_anth, &          ! rec 178 modes 2/3,
                    cnr_floor_gN, cnr_snap_done, cnr_days_bound, cnr_held_gN, &            ! rec 178 modes 2/3,
                    chaff_structural_n_basis, &                                             ! rec 178,/BX
                    W_grain_state => W_grain                                                ! rec 175 mode 4: the grain's OWN dry mass, renamed on USE so it cannot shadow the W_grain DUMMY of allocate_yield_components_daily
  IMPLICIT NONE
  PRIVATE
  ! Crop growth, partition and organ economy.
  ! Primary construction, Teh, C. B. S. (2006) Introduction to Mathematical
  ! Modeling of Crop Growth:
  !   maintenance respiration summed over organs, Eq. [7.3] p.145
  !   age correction by green/total leaf mass,            Eq. [7.4] p.145
  !   maintenance coefficients,                           Table 7.1 p.146
  !   organ-weighted glucose requirement,                 Eq. [7.6] p.148
  !   daily dry-matter increment,                         Eq. [7.7] p.148
  !   glucose requirement coefficients,                   Table 7.4 p.149
  !   partition tables,                                   Listings 8.17-8.20 p.214
  ! Canopy height, Bakhshandeh, E., Soltani, A., Zeinali, E. &
  ! Kallate-Arabi, M. (2012) "Prediction of plant height by allometric
  ! relationships in field-grown wheat", Cereal Research Communications
  ! 40(3):413-422, doi:10.1556/CRC.40.2012.3.10: the segmented model Eq. (4)
  ! p.416 and the pooled plant-height / stem-dry-weight fit Eq. (6) p.419,
  ! Fig. 4.
  ! Root-front descent, Borg & Grimes (1986) as selected by driver record 124.
  ! Leaf senescence, Teh (2006) Listing 8.7 pp.193-194 with the nitrogen-driven
  ! term of Yin et al. (2004) NJAS 51:407-427, Appendix 2.
  PUBLIC :: bgn_set_spike_demand          ! rec 98 spike sink, SHARED.f90
  PUBLIC :: sdm4_set_demand               ! rec 100 MODE 4, SHARED..: the bounded per-axis stem mass demand
  PUBLIC :: inl_advance                   ! rec 109: one day of internode elongation + the E6 three-basis height assembly
  PUBLIC :: initialize_growth_state
  PUBLIC :: bootstrap_emergence_state
  PUBLIC :: bootstrap_nitrogen_state
  PUBLIC :: allocate_daily_biomass
  PUBLIC :: update_root_depth
  PUBLIC :: update_structure
  PUBLIC :: apply_leaf_senescence
  PUBLIC :: allocate_grain_N_source_sink_daily
  PUBLIC :: apply_organN_turnover_daily
  PUBLIC :: update_canopy_leaf_ranks
  PUBLIC :: update_canopy_blade_area
  PUBLIC :: cba_L_at
  PUBLIC :: update_canopy_tillering
  PUBLIC :: ct_L_tiller, ct_W_tiller, ct_nl_tiller, ct_hs_axis
  PUBLIC :: update_leaf_expansion_water   ! rec 83
  PUBLIC :: update_leaf_carbon_bonus      ! rec 144
  PUBLIC :: update_lamina_n_arrest_prov   ! rec 122 pre-commit per-rank increment (commit is internal)
  PUBLIC :: update_canopy_leaf_senesc
  PUBLIC :: update_canopy_sheath_ear
  PUBLIC :: cse_Ls_at, cse_Li_at, cse_Ds_ms, cse_ko_tiller
  PUBLIC :: update_canopy_sink_state
  PUBLIC :: allocate_yield_components_daily
  PUBLIC :: update_spike_axes             ! rec 89
  PUBLIC :: update_spike_window           ! rec 89
  PUBLIC :: sla_from_dvs
  PUBLIC :: partition_fractions
  PUBLIC :: sdm_h_pot_now, sdm_dw_from_dh   ! rec 100
  PUBLIC :: sdm_h_plant_now, sdm_Lped_pot_now, sdm_collar_done   ! rec 100 mode 3

  ! Wheat plant-height allometry, Bakhshandeh et al. (2012) Eq. (6) p.419 and
  ! Fig. 4. The pooled segmented fit up to anthesis is
  !   y = 286.3 x                        for x < 0.16
  !   y = 286.3 x 0.16 + 18.79 (x - 0.16) for x > 0.16
  ! with y the plant height (cm) and x the stem dry weight (g plant-1). Stem
  ! biomass is held as g m-2, so update_structure converts to g plant-1 with the
  ! plant density before the relation is applied.
  REAL(wp), PARAMETER :: stem_height_break_g_plant = 0.16_wp
  REAL(wp), PARAMETER :: stem_height_slope1_cm_g = 286.3_wp
  REAL(wp), PARAMETER :: stem_height_slope2_cm_g = 18.79_wp

  ! ------------------------------------------------------------------ !
  ! Height basis of driver record 100 mode 3. The published height is measured
  ! from the soil surface to the uppermost part of the plant (Bakhshandeh et al.
  ! 2012 p.415), while the quantity the model feeds the inverse relation is the
  ! telescoped COLLAR height, which omits the peduncle: the internode-length
  ! profile stops at the last leaf rank (Abichou 2016 Ch.2 Table 2 Eqn T2.3) and
  ! the thesis lists internode and peduncle as separate organ classes. Mode 3
  ! simulates the missing organ; it adds no observed length to the model.
  ! The two measured lengths below enter exactly one expression, the
  ! dimensionless peduncle to penultimate-internode ratio, and nowhere else:
  ! what the model uses is that ratio applied to its OWN internode length, and
  ! what it reports is a carbon-limited, emergent length.
  ! Source: Rivera-Amado, C., Trujillo-Negrellos, E., Molero, G., Reynolds,
  ! M. P., Sylvester-Bradley, R. & Foulkes, M. J. (2019) Field Crops Research
  ! 240:154-167, doi:10.1016/j.fcr.2019.04.016, Table 6 folio 161, block
  ! "Length (cm)", row "Overall mean": peduncle 36.7, internode 2 20.5,
  ! internode 3 14.1, remainder 14.8 cm; 26 elite spring wheat cultivars, two
  ! irrigated seasons, sampled seven days after anthesis (peduncle h2 = 0.91,
  ! CV = 7.6, measured bracket 22.2-44.2 cm).
  ! Rank correspondence, from the table title: "Int2: internode 2 top-down
  ! (peduncle -1)", so internode 2 is the one immediately below the peduncle,
  ! which is the model top internode. The ratio is applied at the rank at which
  ! it was measured.
  ! ------------------------------------------------------------------ !
  REAL(wp), PARAMETER :: sdm_tab_ped_len_cm  = 36.7_wp   ! PROVENANCE ONLY -- used ONLY in the ratio below
  REAL(wp), PARAMETER :: sdm_tab_int2_len_cm = 20.5_wp   ! PROVENANCE ONLY -- used ONLY in the ratio below
  ! THE ONLY NEW NUMBER THE MODEL USES: a DIMENSIONLESS measured length ratio.
  REAL(wp), PARAMETER :: sdm_ped_int_ratio = sdm_tab_ped_len_cm / sdm_tab_int2_len_cm   ! = 1.79024, -
  ! ------------------------------------------------------------------ !
  INTEGER,  PARAMETER :: stage_emergence = 1
  INTEGER,  PARAMETER :: stage_anthesis = 6
  ! Teh (2006) Table 7.1 p.146: maintenance respiration coefficients at 25 C,
  ! g CH2O g-1 DM d-1.
  REAL(wp), PARAMETER :: km_leaf_25 = 0.030_wp
  REAL(wp), PARAMETER :: km_stem_25 = 0.015_wp
  REAL(wp), PARAMETER :: km_root_25 = 0.015_wp
  REAL(wp), PARAMETER :: km_storage_25 = 0.010_wp
  ! Teh (2006) Table 7.4 p.149: glucose requirement for synthesis,
  ! g CH2O g-1 DM.
  REAL(wp), PARAMETER :: g_leaf = 1.463_wp
  REAL(wp), PARAMETER :: g_stem = 1.513_wp
  REAL(wp), PARAMETER :: g_root = 1.444_wp
  REAL(wp), PARAMETER :: g_storage = 1.415_wp
  ! Teh (2006), around Eq. [7.3], p. 145 and Listing 8.7, p. 194.
  REAL(wp), PARAMETER :: q10_maint = 2.0_wp
  ! 44 g CO2 -> 30 g CH2O, consistent with the textbook CH2O growth currency.
  REAL(wp), PARAMETER :: co2_to_ch2o = 30.0_wp / 44.0_wp

CONTAINS

  SUBROUTINE initialize_growth_state(green_biomass, dead_biomass, stem_biomass, root_biomass, storage_biomass, lai)
    REAL(wp), INTENT(OUT) :: green_biomass, dead_biomass, stem_biomass, root_biomass, storage_biomass
    REAL(wp), INTENT(OUT) :: lai

    ! State initialization only; no published equation is applied here.
    lai = 0.0_wp
    green_biomass = 0.0_wp
    dead_biomass = 0.0_wp
    stem_biomass = 0.0_wp
    root_biomass = 0.0_wp
    storage_biomass = 0.0_wp
  END SUBROUTINE initialize_growth_state

  SUBROUTINE bootstrap_emergence_state(phen_stage, plant_density, init_green_g_plant, init_stem_g_plant, &
                                       init_root_g_plant, init_storage_g_plant, &
                                       green_biomass, dead_biomass, stem_biomass, root_biomass, storage_biomass)
    INTEGER, INTENT(IN) :: phen_stage
    REAL(wp), INTENT(IN) :: plant_density, init_green_g_plant, init_stem_g_plant, init_root_g_plant, init_storage_g_plant
    REAL(wp), INTENT(INOUT) :: green_biomass, dead_biomass, stem_biomass, root_biomass, storage_biomass

    ! Crop-state bootstrap at emergence. This is deliberately not a seed-reserve or
    ! seedling growth model: the cultivar file supplies per-plant initial organ dry
    ! weights and this routine multiplies them by the plant density to obtain the
    ! g m-2 states, so that the initialization is consistent with the stand
    ! density. The values are not taken from Teh (2006) or Bakhshandeh et al.
    ! (2012).
    IF (phen_stage < stage_emergence) RETURN
    IF (green_biomass > 0.0_wp .OR. dead_biomass > 0.0_wp .OR. stem_biomass > 0.0_wp .OR. &
        root_biomass > 0.0_wp .OR. storage_biomass > 0.0_wp) RETURN
    IF (plant_density <= 0.0_wp) RETURN

    green_biomass = MAX(0.0_wp, init_green_g_plant) * plant_density
    dead_biomass = 0.0_wp
    stem_biomass = MAX(0.0_wp, init_stem_g_plant) * plant_density
    root_biomass = MAX(0.0_wp, init_root_g_plant) * plant_density
    storage_biomass = MAX(0.0_wp, init_storage_g_plant) * plant_density
  END SUBROUTINE bootstrap_emergence_state

  SUBROUTINE bootstrap_nitrogen_state(phen_stage, leafN_init_frac, stemN_init_frac, rootN_init_frac, storageN_init_frac, &
                                      green_biomass, dead_biomass, stem_biomass, root_biomass, storage_biomass, &
                                      N_leaf, N_dead_leaf, N_stem, N_root, N_storage, sla_bootstrap)
    INTEGER, INTENT(IN) :: phen_stage
    REAL(wp), INTENT(IN) :: leafN_init_frac, stemN_init_frac, rootN_init_frac, storageN_init_frac
    REAL(wp), INTENT(IN) :: green_biomass, dead_biomass, stem_biomass, root_biomass, storage_biomass
    REAL(wp), INTENT(INOUT) :: N_leaf, N_dead_leaf, N_stem, N_root, N_storage
    REAL(wp), INTENT(IN) :: sla_bootstrap  ! m2 g-1: SLA on the bootstrap day

    ! Nitrogen bootstrap: the organ nitrogen pools are initialized from the
    ! user-supplied organ nitrogen fractions, on the same density basis as the
    ! biomass bootstrap above.
    IF (phen_stage < stage_emergence) RETURN
    IF (N_leaf > 0.0_wp .OR. N_dead_leaf > 0.0_wp .OR. N_stem > 0.0_wp .OR. N_root > 0.0_wp .OR. N_storage > 0.0_wp) RETURN

    IF (leafN_drawdown == 1) THEN

      ! The bootstrap is a leaf-nitrogen INFLOW, so it is bounded by the same leaf
      ! holding capacity as every other inflow:
      !   N_leaf = MIN(leafN_init_frac W_green, sln_max_N LAI_init),
      !   LAI_init = sla_bootstrap W_green,
      ! which is the model's own leaf-area identity (update_structure below).
      ! sln_max_N = 3.8 g N m-2 is the maximum measured field-wheat Narea of
      ! Silva-Perez et al. (2020) J Exp Bot 71(7):2299-2311,
      ! doi:10.1093/jxb/erz439, Fig. 4B p.2305. No concentration output is capped;
      ! this bounds the initialization flux only.
      N_leaf = MIN(MAX(0.0_wp, leafN_init_frac) * green_biomass, &
                   sln_max_N * MAX(0.0_wp, sla_bootstrap) * green_biomass)
    ELSE
      N_leaf = MAX(0.0_wp, leafN_init_frac) * green_biomass
    END IF
    N_dead_leaf = MAX(0.0_wp, leafN_init_frac) * dead_biomass
    N_stem = MAX(0.0_wp, stemN_init_frac) * stem_biomass
    N_root = MAX(0.0_wp, rootN_init_frac) * root_biomass
    N_storage = MAX(0.0_wp, storageN_init_frac) * storage_biomass
  END SUBROUTINE bootstrap_nitrogen_state

  SUBROUTINE allocate_daily_biomass(dvs2, temp_day, net_assim_co2, green_biomass, dead_biomass, &
                                    stem_biomass, root_biomass, storage_biomass, rm_day, rg_day)
    REAL(wp), INTENT(IN) :: dvs2, temp_day, net_assim_co2
    REAL(wp), INTENT(INOUT) :: green_biomass, dead_biomass, stem_biomass, root_biomass, storage_biomass
    REAL(wp), INTENT(OUT) :: rm_day, rg_day

    REAL(wp) :: f_leaf, f_stem, f_root, f_storage
    REAL(wp) :: maint_raw, maint_factor, leaves_total
    REAL(wp) :: net_assim_ch2o, growth_supply, g_total, dm_gain
    REAL(wp) :: wsc_bank, wsc_room          ! (rec 92): pre-anthesis WSC labeling
    REAL(wp) :: dw_leaf, dw_stem, dw_root, dw_storage   ! rec 100: ORGAN INCREMENTS, g DM m-2 d-1

    ! Zero the day mineral delivery ledger at entry, before the two early returns
    ! below, so that a skipped branch can never carry yesterday's value.
    mnp_dw_leaf_day = 0.0_wp ; mnp_dw_root_day = 0.0_wp ; mnp_dw_stem_day = 0.0_wp
    mnp_dw_chaff_day = 0.0_wp ; mnp_dw_grain_day = 0.0_wp
    ! Driver record 121: the day flags of the construction-nitrogen gate, reset
    ! BEFORE the early RETURNs so that a skipped allocation day cannot hand
    ! yesterday's pool to the ladder.
    npg_live_today = .FALSE.
    npg_dw_wsc_day = 0.0_wp ; npg_resid_day = 0.0_wp ; npg_unpaid_ch2o_day = 0.0_wp ; npg_lamina_unpaid_day = 0.0_wp

    ! Dry-matter allocation, Teh (2006):
    !   the partition fractions F_i come from the tables of Listings 8.17-8.20
    !     p.214 (partition_fractions below);
    !   maintenance respiration before the age correction,
    !     RM0 = SUM(kM,i W_i),                              Eq. [7.3] p.145;
    !   age correction,
    !     RM = RM0 W_greenleaves / W_leaves,                Eq. [7.4] p.145;
    !   organ-weighted glucose requirement,
    !     G = SUM(F_i G_i),                                 Eq. [7.6] p.148;
    !   daily dry-matter increment,
    !     W_i,t+1 = W_i,t + F_i (Lambda_canopy - RM)/G dt,  Eq. [7.7] p.148,
    ! with kM,i from Table 7.1 p.146, G_i from Table 7.4 p.149 and dt one day.
    ! Canopy net assimilation enters as g CO2 m-2 d-1 and is converted to
    ! g CH2O m-2 d-1 by the factor 30/44 before Eq. [7.7].
    CALL partition_fractions(dvs2, f_leaf, f_stem, f_root, f_storage)

    ! Closure of the partition to a sum of one is taken from the STEM only: the
    ! lamina stays as tabled, inside the measured lamina partitioning band
    ! 0.18-0.23 of Rivera-Amado et al. (2019), and the root is untouched, so the
    ! spike ramp displaces stem allocation - the organ whose tabled share (stem
    ! plus sheath 0.76) lies above its measured band (0.48-0.61).
    f_stem = MAX(0.0_wp, 1.0_wp - f_leaf - f_root - f_storage)

    leaves_total = green_biomass + dead_biomass

    ! Eq. [7.3] p.145: the organ-summed maintenance demand before the age
    ! correction. Driver record 145 removes the leaf term, because leaf dark
    ! respiration is already charged once inside the Farquhar net assimilation
    ! through Rd; the stem, root and storage terms are the same either way.
    IF (leaf_maintenance_basis >= 1) THEN
      maint_raw = km_stem_25 * stem_biomass + &
                  km_root_25 * root_biomass + km_storage_25 * storage_biomass
    ELSE
      maint_raw = km_leaf_25 * green_biomass + km_stem_25 * stem_biomass + &
                  km_root_25 * root_biomass + km_storage_25 * storage_biomass
    END IF
    ! Q10 temperature correction of the maintenance term, Teh (2006) Listing 8.7
    ! p.194 with Eq. [7.3] p.145.
    maint_factor = q10_maint ** ((temp_day - 25.0_wp) / 10.0_wp)
    rm_day = maint_raw * maint_factor
    ! Eq. [7.4], p. 145: maintenance decreases with green/total leaf fraction.
    IF (leaves_total > 0.0_wp) rm_day = rm_day * (green_biomass / leaves_total)

    ! Convert canopy net CO2 assimilation to the CH2O growth currency before
    ! applying Eq. [7.7] p.148 of Teh (2006).
    net_assim_ch2o = MAX(0.0_wp, net_assim_co2) * co2_to_ch2o
    growth_supply = MAX(0.0_wp, net_assim_ch2o - rm_day)
    ! Driver record 99 modes 2 and 4: export today's supply in the growth currency.
    ! It is read as G(t-1) by tomorrow's update_canopy_tillering and moves no carbon.
    IF (tiller_survival_mech == 4) tca_gsup_day = growth_supply

    ! No structural allocation at or beyond physiological maturity (dvs2 >= 2, the
    ! phenology's own maturity). Post-maturity net assimilate of the residual green
    ! tissue is counted as UNALLOCATED carbon, a labeled diagnostic that is
    ! mass-conserved by entering no organ.
    IF (dvs2 >= 2.0_wp) THEN
      unalloc_C_total = unalloc_C_total + growth_supply
      rg_day = 0.0_wp
      RETURN
    END IF

    ! Eq. [7.6], p. 148: organ-fraction-weighted glucose requirement G.
    g_total = f_leaf * g_leaf + f_stem * g_stem + f_root * g_root + f_storage * g_storage
    IF (g_total <= 0.0_wp) THEN
      rg_day = 0.0_wp
      RETURN
    END IF

    ! Eq. [7.7], p. 148 with dt = 1 day.
    dm_gain = growth_supply / g_total

    ! ---------------------------------------------------------------- !
    ! The canopy becomes a sink (driver record 69). Lamina carbon stops being a
    ! partition FRACTION and becomes an AREA DEMAND: the organ law says how much
    ! new lamina area exists today and this step debits exactly the carbon that
    ! area costs. The remainder is partitioned to stem, root and storage on the
    ! retained fractions renormalized over those three organs, so they still sum
    ! to one and the carbon closure is unchanged. If the demand exceeds the supply
    ! the shortfall is carried as a non-negative BACKLOG and the demographic
    ! limiter acts on FUTURE demand; nothing already built is ever scaled back.
    ! The glucose requirement depends on the partition and the sink changes the
    ! partition, so it is re-solved ONCE with the realized fractions and the
    ! dry-matter gain is re-derived from it. That single corrective iteration is
    ! explicit; the carbon balance closes for any split because the growth
    ! respiration is supply minus gain.
    ! ---------------------------------------------------------------- !
    IF (canopy_sink_lai == 1 .AND. csl_active) THEN
      CALL sink_lamina_split(dm_gain, growth_supply, g_total, &
                             f_leaf, f_stem, f_root, f_storage)
      ! Driver record 82: export the day dry-matter gain for the next day
      ! supply-ratio birth gate (SHARED csl_dm_gain_day). Record 82 requires the
      ! organ law of record 69, so the export lives inside that block.
      csl_dm_gain_day = dm_gain
    END IF

    ! ---------------------------------------------------------------- !
    ! The spike is a sink with priority over the stem reserve, not a share of
    ! supply (driver record 98). Fischer (1985) J Agric Sci Camb 105:447-461,
    ! p.452, states that the stem water-soluble-carbohydrate pool "behaves as a
    ! low priority reserve ... to buffer other sinks such as growing spikes", that
    ! the measured partitioning coefficient "increased as crop growth was reduced
    ! by shading", and that 0.44 is the unshaded end-member only. The spike
    ! therefore carries a DEMAND of its own - potential spike growth per axis per
    ! degree-day - and takes it before the true stem, while the reserve bank,
    ! being a fraction of the stem's own allocation, is served last by construction.
    ! Priority order: lamina (the record-69 area sink above) -> spike -> true stem
    ! -> stem reserve bank. The glucose requirement and the dry-matter gain are
    ! re-solved once with the realized fractions, the same single corrective
    ! iteration the lamina sink performs.
    ! ---------------------------------------------------------------- !
    ! Driver record 100: the true stem stops being the algebraic residual. Up to
    ! this point the four organ increments are f_i dm_gain and therefore satisfy
    ! SUM_i g_i dW_i = growth_supply identically. Below, the stem states a DEMAND -
    ! its own maintenance respiration (Teh Table 7.1 p.146) plus a thermal-time
    ! structural growth demand (the Abichou Ch.2 Fig. 3b and 3c profiles through
    ! the inverted Bakhshandeh Eq. (6) p.419) - and that demand is resolved against
    ! the supply. What the stem does not demand is REFUSED: it enters no organ and
    ! is not renormalized onto leaf, root or storage.
    ! ---------------------------------------------------------------- !
    dw_leaf    = f_leaf    * dm_gain
    dw_stem    = f_stem    * dm_gain
    dw_root    = f_root    * dm_gain
    dw_storage = f_storage * dm_gain
    mnp_dw_leaf_day = dw_leaf             ! P_leaf
    mnp_dw_root_day = dw_root             ! P_root
    sdm_refused_day = 0.0_wp
    sdm_legacy_dw_stem_day = dw_stem

    ! ---------------------------------------------------------------- !
    ! Fate of the refused carbon (driver record 101). Schnyder, H. (1993) New
    ! Phytol 123:233-245, doi:10.1111/j.1469-8137.1993.tb03731.x, p.238, states
    ! that stem reserve storage "accommodates excess photosynthate when the demand
    ! by the ears and other sinks, e.g. respiration and extension growth of the
    ! peduncle, is less than the current net rate of photosynthesis", which is the
    ! record-100 refusal in the source's own words. The capacity is on the same
    ! page, "up to 50 % of the dry weight of internodes", which is wsc_cap_frac.
    ! The conversion 1.242 g CH2O g-1 DM for carbohydrates is Teh (2006) Table 7.2
    ! p.147. room is the closed-form solution of W_reserve + x <= wsc_cap_frac
    ! (S + x), because the stored fructan raises the stem dry weight the ceiling is
    ! expressed on. What the ceiling refuses is REPORTED, not deleted: the same
    ! page measures that wheat shows no negative photosynthetic feedback under sink
    ! removal (Wardlaw 1965; Austin & Edrich 1975; Rawson et al. 1976), so no
    ! assimilation-feedback coefficient is introduced.
    ! ---------------------------------------------------------------- !
    sfl_dw_wsc_day = 0.0_wp
    sfl_want_day   = 0.0_wp
    sfl_room_day   = 0.0_wp
    sfl_resid_day  = sdm_refused_day      ! rec 101 = 0 -> IDENTICAL to the shipped subtrahend
    sfl_stored_dm_cum   = sfl_stored_dm_cum + sfl_dw_wsc_day
    sfl_stored_ch2o_cum = sfl_stored_ch2o_cum + sfl_g_wsc * sfl_dw_wsc_day
    sfl_resid_cum       = sfl_resid_cum + sfl_resid_day
    IF (dvs2 < 1.0_wp) THEN
      sfl_resid_pre  = sfl_resid_pre  + sfl_resid_day
    ELSE
      sfl_resid_post = sfl_resid_post + sfl_resid_day
    END IF

    ! ---------------------------------------------------------------- !
    ! The organ nitrogen payment gate (driver record 121). It runs after the
    ! carbon side has produced the day potential increments and before the organ
    ! states are written. The realized increments are what the nitrogen could pay;
    ! the unpaid carbon goes to the stem reserve up to its capacity and then to an
    ! explicit residual.
    ! ---------------------------------------------------------------- !
    IF (n_payment_gate >= 1) THEN
      CALL n_payment_gate_apply(dvs2, stem_biomass, dw_leaf, dw_stem, dw_root, dw_storage)
      mnp_dw_leaf_day = dw_leaf             ! shares read the REALIZED increments
      mnp_dw_root_day = dw_root

      ! Driver record 122: commit today's lamina-area increment at the paid fraction
      ! p_L and remove the unpaid area from the standing green area this evening.
      ! Driver record 179 reads green_biomass read-only, to form the lamina
      ! PHOTOSYNTHETIC nitrogen on the Bertheloot basis.
      IF (lai_n_arrest >= 1) CALL lamina_n_arrest_commit(npg_p_L, green_biomass)
    END IF

    ! Record 100 mode 1 is the only mode that changes the dry-matter gain; the
    ! other modes leave the next statement an exact no-op, because the increments
    ! are still f_i dm_gain. The stored reserve is real dry matter and must be
    ! inside dm_gain, or the organ-state check below would fail by exactly that
    ! mass. Under record 121 the realized increments plus both reserve deposits ARE
    ! the day dry-matter gain, whatever records 100 and 101 did.
    IF (n_payment_gate >= 1) dm_gain = dw_leaf + dw_stem + dw_root + dw_storage + sfl_dw_wsc_day + npg_dw_wsc_day

    ! Growth respiration is the part of the carbohydrate not retained as dry
    ! matter (Teh 2006, Listing 8.7 p.194).
      rg_day = growth_supply - dm_gain
    ! Driver record 121: the explicit residual is NOT respiration; it is subtracted
    ! so that growth_supply = dm_gain + rg_day + sfl_resid_day + npg_resid_day.
    IF (n_payment_gate >= 1) rg_day = rg_day - npg_resid_day

    ! The daily carbon-fate ledger: supply = allocated + stored + respired +
    ! reported residual. It is exact by construction and is therefore never cited
    ! as closure; it exists so that the four destinations can be read side by side
    ! and so that a change which breaks the arithmetic is caught. The
    ! non-tautological closure test is the ORGAN-STATE check below, which fails by
    ! exactly the stored reserve if the stem update is deleted.
    ! Growth respiration is already inside the g_i factors (g_i > 1 is the
    ! carbohydrate consumed per g dry matter retained), so the carbohydrate-basis
    ! ledger has exactly three destinations and rg_day is its dry-matter-basis
    ! restatement:
    !   growth_supply = SUM_i g_i dW_i + g_wsc dW_wsc + residual
    !   growth_supply = dm_gain        + rg_day       + residual
    sfl_ledger_day = growth_supply - (g_leaf * dw_leaf + g_stem * dw_stem + &
                                      g_root * dw_root + g_storage * dw_storage) - &
                     sfl_g_wsc * sfl_dw_wsc_day - sfl_resid_day
    IF (n_payment_gate >= 1) sfl_ledger_day = sfl_ledger_day - sfl_g_wsc * npg_dw_wsc_day - npg_resid_day   ! rec 121(iii)
    sfl_ledger_max = MAX(sfl_ledger_max, ABS(sfl_ledger_day))

    ! The stem-demand identity, computed every day in every mode. Without the
    ! demand formulation growth_supply - SUM g_i dW_i is zero; with it, the
    ! difference is exactly the refused carbon.
    sdm_ident_day = growth_supply - (g_leaf * dw_leaf + g_stem * dw_stem + &
                                     g_root * dw_root + g_storage * dw_storage)
    sdm_ident_max = MAX(sdm_ident_max, ABS(sdm_ident_day))
    IF (n_payment_gate >= 1) THEN
      ! rec 121: the identity now carries the gate's unpaid carbon as well
      sdm_close_max = MAX(sdm_close_max, ABS(sdm_ident_day - sdm_refused_day - npg_unpaid_ch2o_day))
    ELSE
    sdm_close_max = MAX(sdm_close_max, ABS(sdm_ident_day - sdm_refused_day))
    END IF

    ! ---------------------------------------------------------------- !
    ! The partition's own carbon closure (driver record 98), read from ORGAN STATE
    ! on both sides of the four updates below and compared with the FLUX dm_gain.
    ! It is not algebraically zero for arbitrary inputs: it is zero only because
    ! the stem fraction is renormalized, and deleting that one line takes it to
    ! about 0.43 dm_gain. The two shipped closure columns of the reserve and sink
    ! files stay at machine precision when the renormalization is removed, which is
    ! why this check exists separately.
    ! ---------------------------------------------------------------- !
    BLOCK
      REAL(wp) :: bgn_sum0, bgn_sum1
      bgn_sum0 = green_biomass + dead_biomass + stem_biomass + root_biomass + storage_biomass
    ! P_stem is the REALIZED structural stem increment, read AFTER the record-100
    ! resolution and BEFORE the record-101 water-soluble-carbohydrate deposit, so
    ! that stored carbohydrate is excluded from the structural basis.
    mnp_dw_stem_day = dw_stem
    ! Eq. [7.7], p. 148: organ-level update W_i,t+1 = W_i,t + F_i * dm_gain.
    green_biomass = green_biomass + dw_leaf
    stem_biomass = stem_biomass + dw_stem

    ! The destination of the refused carbon (driver record 101). The stored fructan
    ! IS stem dry weight - Schnyder (1993) p.238 states the capacity as a
    ! percentage of the dry weight of internodes - and W_reserve is a labeled part
    ! of W_stem, the remobilization draw being W_stem = W_stem - draw, so a pool
    ! outside W_stem would be drawn out of mass that does not exist. Both updates
    ! sit inside the organ-state check that brackets this section: deleting the
    ! stem update alone makes that check fail by exactly the deposited mass.
    ! The same holds for the unpaid carbon stored under record 121.
    IF (n_payment_gate >= 1 .AND. npg_dw_wsc_day > 0.0_wp) THEN
      stem_biomass = stem_biomass + npg_dw_wsc_day
      W_reserve    = W_reserve    + npg_dw_wsc_day

      ! The record-121 reserve payment is labeled deposition for the reserve ledger,
      ! exactly as the record-101 deposit and the record-92 route already are.
      wsc_dep_day  = wsc_dep_day  + npg_dw_wsc_day
    END IF

    ! Carbon-reserve banking (driver record 55 mode 1, the pre-anthesis window
    ! only): a fraction of today's stem dry-matter gain is LABELED as reserve.
    ! Bookkeeping only - no mass moves. The fraction basis is the measured stem
    ! reserve concentration at anthesis of Ruuska et al. (2006) Funct Plant Biol
    ! 33:799-809; what is tested is the outcome, against the at-anthesis envelope
    ! of Schnyder (1993) New Phytol 123:233-245 p.238.
    ! Under driver record 92 the same arithmetic is bounded by the capacity
    ! wsc_cap_frac W_stem (Schnyder 1993 p.238, "up to 50 % of the dry weight of
    ! internodes"). Carbon above the capacity is simply NOT LABELED: it is already
    ! in the stem pool and stays there, so no mass moves and the carbon ledger is
    ! untouched.
    IF (carbon_reserve_mode == 1 .AND. dvs2 < reserve_start_dvs2) THEN

      ! This statement is outside the record-100 guard and is live on every deck, so
      ! it must keep the same expression whenever record 100 is not 1: the
      ! increments are the same quantity there, but (a b) c and a (b c) are not
      ! bit-equal in IEEE-754. When the stem demand is on, the increment is no
      ! longer f_stem dm_gain and the reserve must be banked on the realized stem
      ! increment.
      ! ================================================================ !
      ! Driver record 172 makes the deposition SHARE a function of the nitrogen
      ! nutrition index instead of a constant (van Herwaarden 1997 PhD thesis Ch.4).
      ! ================================================================ !
      IF (stem_demand_mech == 1 .OR. stem_demand_mech == 3 .OR. stem_demand_mech == 4 .OR. n_payment_gate >= 1) THEN
        wsc_bank = stem_reserve_frac * dw_stem      ! rec 121: the REALIZED stem increment (same rule as rec 100 ON)
      END IF
      IF (stem_wsc_mech == 1 .OR. stem_wsc_mech == 2 .OR. stem_wsc_mech == 4) THEN
        wsc_room = MAX(0.0_wp, wsc_cap_frac * stem_biomass - W_reserve)
        IF (wsc_bank > wsc_room) THEN
          wsc_refused_cap = wsc_refused_cap + (wsc_bank - wsc_room)
          wsc_days_capbnd = wsc_days_capbnd + 1
          wsc_bank = wsc_room
        END IF
      END IF
      W_reserve = W_reserve + wsc_bank
      IF (stem_wsc_mech >= 1) THEN
        wsc_dep_pre = wsc_dep_pre + wsc_bank
        wsc_dep_day = wsc_dep_day + wsc_bank
      END IF
    END IF
    root_biomass = root_biomass + dw_root
    storage_biomass = storage_biomass + dw_storage

    ! The reserve/stem ratio is read AFTER the record-55 bank, because that bank
    ! deposits without a capacity of its own whenever record 92 is 0. It is
    ! reported as a falsification metric.
      bgn_sum1 = green_biomass + dead_biomass + stem_biomass + root_biomass + storage_biomass
      IF (ABS(bgn_sum1 - bgn_sum0 - dm_gain) > bgn_dm_close_max) &
        bgn_dm_close_max = ABS(bgn_sum1 - bgn_sum0 - dm_gain)
      bgn_agb_day = green_biomass + dead_biomass + stem_biomass + storage_biomass
    END BLOCK

    ! Export the day increments on Fischer's own basis (driver record 89). The
    ! above-ground increment is taken as dm_gain (1 - f_root) from the REALIZED
    ! fractions rather than from the reported above-ground column, so that the
    ! realized window share is independent of the reporting basis.
    IF (grain_number_mech >= 1) THEN
      ! This statement is OUTSIDE the record-100 guard, so modes 0 and 2 evaluate the
      ! same expression as before.
      IF (stem_demand_mech == 1 .OR. stem_demand_mech == 3 .OR. stem_demand_mech == 4 .OR. n_payment_gate >= 1) THEN
        gnx_dagb_day   = dw_leaf + dw_stem + dw_storage   ! rec 121: realized increments (dm_gain now carries WSC)
      ELSE
        gnx_dagb_day   = dm_gain * (1.0_wp - f_root)
      END IF
      gnx_dspike_day = dw_storage
    END IF
  END SUBROUTINE allocate_daily_biomass

  ! ==================================================================== !
  ! Organ nitrogen payment for structural growth (driver record 121). Called from
  ! allocate_daily_biomass after the carbon side has produced the day potential
  ! organ increments and before the organ states are written. Each increment is
  ! paid in reduced nitrogen at its tier concentration c_k; the paid fraction p_k
  ! scales the realized increment; the carbon of the unpaid part goes to the stem
  ! reserve up to the Schnyder (1993) capacity and the remainder to an EXPLICIT
  ! residual (g CH2O m-2 d-1), never to an organ.
  ! ==================================================================== !
  SUBROUTINE n_payment_gate_apply(dvs2, stem_biomass, dw_leaf, dw_stem, dw_root, dw_storage)
    REAL(wp), INTENT(IN)    :: dvs2, stem_biomass
    REAL(wp), INTENT(INOUT) :: dw_leaf, dw_stem, dw_root, dw_storage
    REAL(wp) :: pool, pool_root, lux, take, lux_take, dws_pot, f_true_stem, from_pool, mobile_used
    REAL(wp) :: want, want_R, room, s_next, ident
    REAL(wp) :: lnr_excess, lnr_rem_avail, lnr_rem_left, lnr_gap, lnr_take   ! rec 141
    REAL(wp) :: nsc_pfloor_L, nsc_pfloor_S, nsc_pfloor_R                    ! rec 143: floor-only counterfactual paid fractions

    npg_dw_leaf_pot = dw_leaf ; npg_dw_stem_pot = dw_stem
    npg_dw_stor_pot = dw_storage ; npg_dw_root_pot = dw_root
    npg_dw_leaf = dw_leaf ; npg_dw_stem = dw_stem ; npg_dw_stor = dw_storage ; npg_dw_root = dw_root
    f_true_stem = 0.0_wp   ! rec 142: defensive init, compiler flow-analysis only (always set before use)
    npg_p_L = 1.0_wp ; npg_p_S = 1.0_wp ; npg_p_R = 1.0_wp
    npg_p_chaff = 1.0_wp ; npg_p_stem = 1.0_wp   ! rec 138
    npg_D_L = 0.0_wp ; npg_D_S = 0.0_wp ; npg_D_R = 0.0_wp
    npg_D_chaff = 0.0_wp ; npg_D_stem = 0.0_wp ; ncs_D2_chaff = 0.0_wp ; ncs_D2_stem = 0.0_wp   ! rec 138
    npg_paid_L = 0.0_wp ; npg_paid_S = 0.0_wp ; npg_paid_chaff = 0.0_wp ; npg_paid_R = 0.0_wp
    npg_paid_R_root = 0.0_wp ; npg_paid_R_shoot = 0.0_wp ; npg_root_export = 0.0_wp
    npg_lux_avail = 0.0_wp ; npg_lux_used = 0.0_wp
    npg_pool_mobile_in = 0.0_wp ; npg_pool_root_in = 0.0_wp ; npg_pool_assim_in = 0.0_wp
    npg_pool_out = 0.0_wp ; npg_mobile_left = 0.0_wp
    lnr_excess_day = 0.0_wp ; lnr_rem_avail_day = 0.0_wp ; lnr_rem_used_day = 0.0_wp   ! rec 141
    nsc_status_day = 0.0_wp ; nsc_price_L = 0.0_wp ; nsc_price_S = 0.0_wp               ! rec 143
    nsc_price_chaff = 0.0_wp ; nsc_price_R = 0.0_wp
    nsc_cut_L_day = 0.0_wp ; nsc_cut_S_day = 0.0_wp ; nsc_cut_R_day = 0.0_wp
    nsc_pfloor_L = 1.0_wp ; nsc_pfloor_S = 1.0_wp ; nsc_pfloor_R = 1.0_wp
    ! The vegetative window of the nitrogen ladder (dvs2 < 1). After anthesis the
    ! day reduced nitrogen is the grain-nitrogen source (Martre et al. 2003) and no
    ! structural organ is built from it, so the gate is inactive there.
    IF (dvs2 >= 1.0_wp) RETURN
    npg_live_today = .TRUE.
    npg_days_live  = npg_days_live + 1

    ! Tier concentrations, g N g-1 DM of the increment. Under driver record 128
    ! the first pass always pays the FLOOR, independent of the record-121 mode:
    ! record 128 replaces the single-tier choice. Its modes 2 and 3 change only the
    ! lamina floor; the stem and root floors are the mode-1 values in every mode.
    IF (n_construction_status >= 1) THEN
      SELECT CASE (n_construction_status)
      CASE (2)
        ncs_c_floor_L_eff = SLN_min_N / csl_LMA  ! model's own rec-122 minimum functional SLN
      END SELECT
      npg_c_L = ncs_c_floor_L_eff ; npg_c_S = n5_stem_nstruct ; npg_c_R = npg_floor_R
      ncs_c_target_L = n4_sln_target / csl_LMA
      ncs_c_target_S = n4_stem_nmax
      ncs_c_target_R = npg_build_R
    END IF
    ! ----------------------------------------------------------------------------------- !
    ! Driver record 142 overrides the lamina tier price onto a functional specific-
    ! leaf-nitrogen basis. Mode 1 is n4_sln_target/csl_LMA, the 1.5 g N m-2 leaf
    ! nitrogen per unit area of Jamieson & Semenov (2000) Field Crops Res 68:21-29,
    ! pp.23-24; mode 2 is SLN_min_N/csl_LMA, the model's own 0.50 g N m-2 functional
    ! minimum. The stem, chaff and root prices are untouched.
    ! Driver record 138 sets the chaff tier's own floor concentration: mode 1 is the
    ! Jamieson & Semenov stem floor, mode 2 the Bertheloot chaff floor.
    ! The order in which the pools pay:
    !   (1) today's root-reduced nitrate pays the root tier first, before any of it
    !       is exported to the shoot (Andrews 1986: root nitrate reductase reduces
    !       the share retained in the root; what the root does not keep goes up the
    !       xylem);
    !   (2) its surplus is exported and joins the SHOOT pool (the mobile vegetative
    !       pool plus today's shoot-reduced nitrogen), which pays lamina, then stem
    !       and chaff;
    !   (3) what the shoot pool still holds then tops up any unpaid root remainder,
    !       the phloem return of amino nitrogen to the root (Cooper & Clarkson 1989
    !       J Exp Bot 40:753, doi:10.1093/jxb/40.7.753, at the level of the
    !       mechanism; the step itself is a model construction);
    !   (4) the stem nitrogen above its structural floor, released at npg_k_lux
    !       d-1, is the last resort for the shoot tiers and the root remainder
    !       alike.
    ! A lamina-first order starves the root (root dry mass falls 22 to 59 %, the
    ! sign opposite to the measured shoot/root response), which is why the root is
    ! served first.
    ! ----------------------------------------------------------------------------------- !
    npg_pool_mobile_in = N_mobile_veg
    npg_pool_root_in   = npg_n_assim_root_in
    npg_pool_assim_in  = npg_n_assim_root_in + npg_n_assim_shoot_in
    N_mobile_veg = 0.0_wp
    lux = npg_k_lux * MAX(0.0_wp, npg_N_stem - n5_stem_nstruct * stem_biomass)
    npg_lux_avail = lux

    ! Demands on the carbon-side potentials.
    dws_pot = dw_stem + dw_storage                 ! true stem + spike structure (chaff): one tier

    IF (n_construction_status >= 1) THEN
      ! --------------------------------------------------------------------------------- !
      ! Shoot demand ceiling as the marginal form of the nitrogen dilution curve of
      ! Justes, Mary, Meynard, Machet & Thelier-Huche (1994) Ann Bot 74:397-407. 
      ! The above-ground dry mass is set by main_MMCO2Crop.f90 just before this call, on 
      ! the SAME basis as the live critical-concentration and nitrogen-nutrition-index 
      ! terms.
      ! --------------------------------------------------------------------------------- !
      ncs_W_eff  = MAX(ncs_W_t_ha, 1.55_wp) ! Justes/Pask validity floor, Pask 2009 p.150-151

      ! ================================================================== !
      ! Driver record 156 selects WHICH dilution curve caps the organ construction
      ! nitrogen targets. Mode 0 is the CRITICAL pair, taken on the ELSE line
      ! below. Mode 1 is the same marginal form a(1-b) W^-b on the same validity
      ! floor with the MAXIMUM pair (8.3, 0.44) - the model's own existing
      ! maximum-curve coefficients, which are the maximum-curve demand of
      ! Brisson et al. (1998) Agronomie 18:311-346, Eq. (47) p.342 (the formulation
      ! implemented in STICS), with the
      ! coefficients attributed in its Table VIII p.330 to Justes et al. (1994).
      ! No new constant is introduced. Justes et al. (1994) print Eq. (17) p.405,
      ! Nmax = 8.3 DM^-0.44, so 8.3 is the Justes value. That paper prints no
      ! maximum below 1.5 t ha-1, which is why the dry mass is held at or above 1.55
      ! t ha-1 and the derivative is always evaluated rather than the plateau.
      ! ================================================================== !
      IF (n_construction_ceiling_curve >= 1) THEN
        ncs_c_crit = (8.3_wp * (1.0_wp - 0.44_wp) * ncs_W_eff**(-0.44_wp)) * 0.01_wp   ! %DM -> g N g-1 DM
        ncc_a_used = 8.3_wp ; ncc_b_used = 0.44_wp
        ncc_days_live = ncc_days_live + 1
      END IF
      ncs_D_crit = ncs_c_crit * (dw_leaf + dws_pot)
      ncs_D_target_L = ncs_c_target_L * dw_leaf
      ncs_D_target_S = ncs_c_target_S * dws_pot
      ncs_f_crit = 1.0_wp
      IF (ncs_D_target_L + ncs_D_target_S > 0.0_wp) &
        ncs_f_crit = MIN(1.0_wp, ncs_D_crit / (ncs_D_target_L + ncs_D_target_S))
      IF (ncs_f_crit < 1.0_wp) ncs_days_fcrit = ncs_days_fcrit + 1
      ! The lower bound uses the ACTIVE mode's own lamina floor, so that the effective
      ! target never falls below the floor pass 1 has just paid, in every mode.
      ncs_c_target_eff_L = MAX(ncs_c_floor_L_eff, ncs_f_crit * ncs_c_target_L)
      ncs_c_target_eff_S = MAX(n5_stem_nstruct,   ncs_f_crit * ncs_c_target_S)
      ncs_c_target_eff_R = ncs_c_target_R                 ! root target not capped by the shoot curve
      ncs_D2_L = MAX(0.0_wp, ncs_c_target_eff_L * dw_leaf - ncs_c_floor_L_eff * dw_leaf)
        ncs_D2_S = MAX(0.0_wp, ncs_c_target_eff_S * dws_pot - n5_stem_nstruct   * dws_pot)
      ncs_D2_R = MAX(0.0_wp, ncs_c_target_eff_R * dw_root - npg_floor_R       * dw_root)
    END IF

    ! ==================================================================== !
    ! Driver record 143: a whole-crop nitrogen-status-dependent construction price
    ! for every structural organ. It overrides the lamina, stem, chaff and root
    ! prices onto a price that interpolates LINEARLY between the organ existing
    ! floor and its existing record-128 status-capped target, using the whole-shoot
    ! vegetative nitrogen nutrition index computed in main_MMCO2Crop.f90. No new
    ! constant. Mode 1 is an inactive alternative, see Part 1 Table 6.
    ! ==================================================================== !
    IF (n_status_construction_price >= 1) THEN
      nsc_status_day = MIN(1.0_wp, MAX(0.0_wp, nap_nni_veg))
      nsc_nni_veg_sum  = nsc_nni_veg_sum + nap_nni_veg
      nsc_nni_veg_days = nsc_nni_veg_days + 1
      nsc_price_L = npg_floor_L     + nsc_status_day * (ncs_c_target_eff_L - npg_floor_L)
      nsc_price_S = n5_stem_nstruct + nsc_status_day * (ncs_c_target_eff_S - n5_stem_nstruct)
      nsc_price_R = npg_floor_R     + nsc_status_day * (ncs_c_target_eff_R - npg_floor_R)
      npg_c_L = nsc_price_L
      npg_c_S = nsc_price_S
      npg_c_R = nsc_price_R
        ! Combined mode: the chaff has no separate tier, it shares the lamina
        ! construction concentration, which is the same treatment the combined floor
        ! applies.
        nsc_price_chaff = nsc_price_S
      ! --------------------------------------------------------------------------------- !
      ! Diagnostic-only floor counterfactual: what the first pass would pay at the
      ! FLOOR price alone, from the SAME day pool inputs, replicating steps (1) to
      ! (4) above only, in a self-contained block that touches no shared state and is
      ! never booked into any payment.
      ! --------------------------------------------------------------------------------- !
      BLOCK
        REAL(wp) :: d0_L, d0_chaff, d0_stem, d0_S, d0_R
        REAL(wp) :: pool0, pool0_root0, lux0, take0, luxtake0, want0, paid_R_root0
        d0_L = npg_floor_L * npg_dw_leaf_pot
          d0_chaff = 0.0_wp
          d0_stem  = n5_stem_nstruct * (npg_dw_stem_pot + npg_dw_stor_pot)
        d0_S = d0_chaff + d0_stem
        d0_R = npg_floor_R * npg_dw_root_pot
        lux0 = npg_lux_avail
        paid_R_root0 = 0.0_wp
        IF (d0_R > 0.0_wp) paid_R_root0 = MIN(npg_pool_root_in, d0_R)
        pool0_root0 = npg_pool_root_in - paid_R_root0
        pool0 = npg_pool_mobile_in + npg_n_assim_shoot_in + pool0_root0
        nsc_pfloor_L = 1.0_wp ; nsc_pfloor_S = 1.0_wp ; nsc_pfloor_R = 1.0_wp
        IF (d0_L > 0.0_wp) THEN
          take0 = MIN(pool0, d0_L) ; pool0 = pool0 - take0
          IF (take0 < d0_L) THEN
            luxtake0 = MIN(lux0, d0_L - take0) ; lux0 = lux0 - luxtake0 ; take0 = take0 + luxtake0
          END IF
          nsc_pfloor_L = take0 / d0_L
        END IF
        IF (d0_S > 0.0_wp) THEN
          take0 = MIN(pool0, d0_S) ; pool0 = pool0 - take0
          IF (take0 < d0_S) THEN
            luxtake0 = MIN(lux0, d0_S - take0) ; lux0 = lux0 - luxtake0 ; take0 = take0 + luxtake0
          END IF
          nsc_pfloor_S = take0 / d0_S
        END IF
        IF (d0_R > 0.0_wp) THEN
          want0 = d0_R - paid_R_root0
          take0 = MIN(pool0, want0) ; pool0 = pool0 - take0
          IF (take0 < want0) THEN
            luxtake0 = MIN(lux0, want0 - take0) ; lux0 = lux0 - luxtake0 ; take0 = take0 + luxtake0
          END IF
          nsc_pfloor_R = (paid_R_root0 + take0) / d0_R
        END IF
      END BLOCK
    END IF

    npg_D_L = npg_c_L * dw_leaf
      npg_D_S = npg_c_S * dws_pot
    npg_D_R = npg_c_R * dw_root

    ! step (1): root tier from root-reduced N, before export
    pool_root = npg_pool_root_in
    IF (npg_D_R > 0.0_wp) THEN
      take = MIN(pool_root, npg_D_R) ; pool_root = pool_root - take
      npg_paid_R_root = take
    END IF
    npg_root_export = pool_root
    ! step (2): the shoot pool pays lamina -> stem/chaff (luxury stem N last resort)
    pool = npg_pool_mobile_in + npg_n_assim_shoot_in + npg_root_export
    IF (npg_D_L > 0.0_wp) THEN
      take = MIN(pool, npg_D_L) ; pool = pool - take
      IF (take < npg_D_L) THEN
        lux_take = MIN(lux, npg_D_L - take) ; lux = lux - lux_take
        npg_lux_used = npg_lux_used + lux_take ; take = take + lux_take
      END IF
      npg_p_L = take / npg_D_L ; npg_paid_L = take
      IF (npg_p_L < 1.0_wp) npg_days_bind_L = npg_days_bind_L + 1
    END IF

    ! ==================================================================== !
    ! Driver record 142, the carbon re-offer and the structural-demand recompute.
    ! The lamina dry matter the nitrogen tier just refused is converted to
    ! carbohydrate at g_leaf and re-offered THE SAME DAY to the true stem, the
    ! spike and the root in proportion to their current carbon-side potentials; the
    ! structural nitrogen demands below are then recomputed on the enlarged
    ! potentials, before those tiers draw on the shoot nitrogen pool.
    ! ==================================================================== !
    IF (npg_D_S > 0.0_wp) THEN
      take = MIN(pool, npg_D_S) ; pool = pool - take
      IF (take < npg_D_S) THEN
        lux_take = MIN(lux, npg_D_S - take) ; lux = lux - lux_take
        npg_lux_used = npg_lux_used + lux_take ; take = take + lux_take
      END IF
      npg_p_S = take / npg_D_S
      f_true_stem    = dw_stem / dws_pot
      npg_paid_S     = take * f_true_stem
      npg_paid_chaff = take - npg_paid_S
      IF (npg_p_S < 1.0_wp) npg_days_bind_S = npg_days_bind_S + 1
    END IF
    ! step (3)+(4): root remainder from the shoot pool, then luxury stem N
    IF (npg_D_R > 0.0_wp) THEN
      want_R = npg_D_R - npg_paid_R_root
      take = MIN(pool, want_R) ; pool = pool - take
      IF (take < want_R) THEN
        lux_take = MIN(lux, want_R - take) ; lux = lux - lux_take
        npg_lux_used = npg_lux_used + lux_take ; take = take + lux_take
      END IF
      npg_paid_R_shoot = take
      npg_paid_R = npg_paid_R_root + npg_paid_R_shoot
      npg_p_R = npg_paid_R / npg_D_R
      IF (npg_p_R < 1.0_wp) npg_days_bind_R = npg_days_bind_R + 1
    END IF
    npg_pool_out = pool

    ! ==================================================================== !
    ! Driver record 141: pre-anthesis lamina nitrogen remobilization, the last
    ! resort for the stem and chaff structural tiers only, after the root
    ! remainder and the stem luxury nitrogen above and before the second pass.
    ! ==================================================================== !
    IF (lamina_n_remob_pre >= 1) THEN
      ! ---------------------------------------------------------------------------- !
      ! The same basis driver record 122 mode 2 uses (need = green area index minus
      ! N_leaf / SLN_min_N), inverted to a mass instead of an area, with the
      ! Eq. (5) rate of Bertheloot, Martre & Andrieu (2008) Funct Plant Biol
      ! 35:781-796 p.784 already used post-anthesis.
      ! ---------------------------------------------------------------------------- !
      lnr_excess    = MAX(0.0_wp, npg_N_leaf - SLN_min_N * cls_gai)
      lnr_rem_avail = MIN(lnr_excess, n2_delta_TT * MAX(0.0_wp, dTT_day_N) * lnr_excess)
      lnr_rem_left  = lnr_rem_avail
      ! ---------------------------------------------------------------------------- !
      ! Chaff first, then the true stem, which is the priority order of Fischer
      ! (1985) J Agric Sci Camb 105:447-461, p.452; the combined tier is split by the
      ! true-stem fraction exactly as the combined branch of driver record 121.
      ! ---------------------------------------------------------------------------- !
      IF (lnr_rem_left > 0.0_wp .AND. npg_D_S > 0.0_wp) THEN
        lnr_gap  = MAX(0.0_wp, npg_D_S - (npg_paid_S + npg_paid_chaff))
        lnr_take = MIN(lnr_rem_left, lnr_gap)
        IF (lnr_take > 0.0_wp) THEN
          npg_paid_S     = npg_paid_S     + lnr_take * f_true_stem
          npg_paid_chaff = npg_paid_chaff + lnr_take * (1.0_wp - f_true_stem)
          lnr_rem_left   = lnr_rem_left   - lnr_take
        END IF
        npg_p_S = (npg_paid_S + npg_paid_chaff) / npg_D_S
      END IF
      ! ---------------------------------------------------------------------------- !
      ! Debit the leaf nitrogen THE SAME DAY, so that the record-122 mode-2 shed test
      ! tomorrow sees the lowered standing nitrogen.
      ! ---------------------------------------------------------------------------- !
      lnr_rem_used_day = lnr_rem_avail - lnr_rem_left
      IF (lnr_rem_used_day > 0.0_wp) THEN
        npg_N_leaf       = npg_N_leaf       - lnr_rem_used_day
        lnr_rem_used_cum = lnr_rem_used_cum + lnr_rem_used_day
        lnr_days_rem_used = lnr_days_rem_used + 1
      END IF
      lnr_excess_day    = lnr_excess
      lnr_rem_avail_day = lnr_rem_avail
      IF (gnx_win_open .AND. dvs2 < 1.0_wp) lnr_win_rem_used_sum = lnr_win_rem_used_sum + lnr_rem_used_day
    END IF
    ! ------------------------------------------------------------------------------ !
    ! Driver record 143 diagnostic: the dry-matter cut relative to the floor-only
    ! counterfactual computed above, against the SAME day pool. Read here, after
    ! every first-pass payment this record price can affect, and before the second
    ! pass, which never touches the paid fractions. Reported only; nothing is
    ! booked and no state is fed back.
    ! ------------------------------------------------------------------------------ !
    IF (n_status_construction_price >= 1) THEN
      nsc_cut_L_day = MAX(0.0_wp, nsc_pfloor_L - npg_p_L) * npg_dw_leaf_pot
      nsc_cut_S_day = MAX(0.0_wp, nsc_pfloor_S - npg_p_S) * (npg_dw_stem_pot + npg_dw_stor_pot)
      nsc_cut_R_day = MAX(0.0_wp, nsc_pfloor_R - npg_p_R) * npg_dw_root_pot
      nsc_cut_L_cum = nsc_cut_L_cum + nsc_cut_L_day
      nsc_cut_S_cum = nsc_cut_S_cum + nsc_cut_S_day
      nsc_cut_R_cum = nsc_cut_R_cum + nsc_cut_R_day
    END IF
    ! ------------------------------------------------------------------------------ !
    ! Second pass (driver record 128 only): pay the remaining target demand from
    ! what the first pass left in the shoot pool, with no luxury draw and in the
    ! same order (lamina, then stem and chaff, then root, whose second-pass share
    ! is a shoot-pool phloem return exactly as in step (3) above). The paid
    ! fractions were fixed in the first pass against the floor and are not
    ! recomputed. The first-pass totals are snapshotted so that both are reported.
    ! ------------------------------------------------------------------------------ !
    ncs_paid2_L = 0.0_wp ; ncs_paid2_S = 0.0_wp ; ncs_paid2_chaff = 0.0_wp ; ncs_paid2_R = 0.0_wp
    ncs_lux2_cap_day = 0.0_wp ; ncs_lux2_used_day = 0.0_wp
    IF (n_construction_status >= 1) THEN
      ncs_paid1_L = npg_paid_L ; ncs_paid1_S = npg_paid_S
      ncs_paid1_chaff = npg_paid_chaff ; ncs_paid1_R = npg_paid_R
      pool = npg_pool_out
      IF (ncs_D2_L > 0.0_wp) THEN
        take = MIN(pool, ncs_D2_L) ; pool = pool - take
        ! -------------------------------------------------------------------------- !
        ! Driver record 135: a second, bounded luxury tranche for the lamina
        ! second-pass gap only, capped at the SAME rate the first pass used
        ! and applied again to whatever of the same stem-luxury pool the first pass
        ! left unused. No new parameter.
        ! -------------------------------------------------------------------------- !
        ncs_paid2_L = take
      END IF
      IF (ncs_D2_S > 0.0_wp) THEN
        take = MIN(pool, ncs_D2_S) ; pool = pool - take
        ncs_paid2_S     = take * f_true_stem
        ncs_paid2_chaff = take - ncs_paid2_S
      END IF
      IF (ncs_D2_R > 0.0_wp) THEN
        take = MIN(pool, ncs_D2_R) ; pool = pool - take
        ncs_paid2_R = take
      END IF
      npg_pool_out     = pool
      npg_paid_L       = npg_paid_L       + ncs_paid2_L
      npg_paid_S       = npg_paid_S       + ncs_paid2_S
      npg_paid_chaff   = npg_paid_chaff   + ncs_paid2_chaff
      npg_paid_R       = npg_paid_R       + ncs_paid2_R
      npg_paid_R_shoot = npg_paid_R_shoot + ncs_paid2_R
      ncs_lux2_used_cum = ncs_lux2_used_cum + ncs_lux2_used_day
      IF (ncs_lux2_used_day > 0.0_wp) ncs_lux2_days = ncs_lux2_days + 1
    END IF

    ! Book the nitrogen. N_chaff is the chaff-tagged subset of N_storage, so the
    ! chaff share is added to both.
    npg_N_leaf    = npg_N_leaf    + npg_paid_L
    npg_N_stem    = npg_N_stem    + npg_paid_S - npg_lux_used
    N_chaff       = N_chaff       + npg_paid_chaff
    npg_N_storage = npg_N_storage + npg_paid_chaff
    npg_N_root    = npg_N_root    + npg_paid_R
    ! ------------------------------------------------------------------------------ !
    ! Mobile-first accounting for the overflow rule: the shoot pool is spent mobile
    ! first, and the root's own share never touches the mobile pool. The lamina
    ! nitrogen released under record 141 is excluded here exactly as the stem
    ! luxury nitrogen already is, because neither came from the mobile or
    ! assimilated shoot pool.
    ! ------------------------------------------------------------------------------ !
    from_pool   = npg_paid_L + npg_paid_S + npg_paid_chaff + npg_paid_R_shoot - npg_lux_used - lnr_rem_used_day
    mobile_used = MIN(npg_pool_mobile_in, MAX(0.0_wp, from_pool))
    npg_mobile_left = npg_pool_mobile_in - mobile_used
    ! ------------------------------------------------------------------------------ !
    ! The gate nitrogen identity, machine zero by construction. The lamina release
    ! of record 141 is a new SOURCE; the sink side is unchanged in form because the
    ! stem and chaff payments already carry the top-up added above.
    ! ------------------------------------------------------------------------------ !
    ident = (npg_pool_mobile_in + npg_pool_assim_in + npg_lux_used + lnr_rem_used_day) - &
            (npg_paid_L + npg_paid_S + npg_paid_chaff + npg_paid_R + npg_pool_out)
    npg_ident_N_max = MAX(npg_ident_N_max, ABS(ident))
    IF (ABS(ident) > 1.0e-9_wp) THEN
      WRITE (*,'(A,G14.6)') ' ERROR rec 121: gate nitrogen identity broken, |resid| (g N m-2) = ', ABS(ident)
      ERROR STOP 'rec 121: nitrogen payment identity'
    END IF
    npg_paid_cum_L = npg_paid_cum_L + npg_paid_L
    npg_paid_cum_S = npg_paid_cum_S + npg_paid_S + npg_paid_chaff
    npg_paid_cum_R = npg_paid_cum_R + npg_paid_R
    npg_paid_cum_R_root  = npg_paid_cum_R_root  + npg_paid_R_root
    npg_paid_cum_R_shoot = npg_paid_cum_R_shoot + npg_paid_R_shoot

    ! Realized increments and the unpaid carbon.
    dw_leaf    = npg_p_L * npg_dw_leaf_pot
      dw_stem    = npg_p_S * npg_dw_stem_pot
      dw_storage = npg_p_S * npg_dw_stor_pot
    dw_root    = npg_p_R * npg_dw_root_pot
    npg_dw_leaf = dw_leaf ; npg_dw_stem = dw_stem ; npg_dw_stor = dw_storage ; npg_dw_root = dw_root
    ! The realized construction concentration, i.e. the nitrogen percentage of the
    ! tissue actually built today, and the cumulative lamina dry matter built below
    ! its status-capped target. By construction floor <= c_actual <= target_eff.
    IF (n_construction_status >= 1) THEN
      ncs_c_actual_L = 0.0_wp ; ncs_c_actual_S = 0.0_wp ; ncs_c_actual_R = 0.0_wp
      IF (dw_leaf > 0.0_wp) ncs_c_actual_L = npg_paid_L / dw_leaf
      IF (dw_stem + dw_storage > 0.0_wp) ncs_c_actual_S = (npg_paid_S + npg_paid_chaff) / (dw_stem + dw_storage)
      IF (dw_root > 0.0_wp) ncs_c_actual_R = npg_paid_R / dw_root
      IF (dw_leaf > 0.0_wp .AND. ncs_c_actual_L < ncs_c_target_eff_L - 1.0e-12_wp) &
        ncs_lamina_below_target_cum = ncs_lamina_below_target_cum + dw_leaf
    END IF
    npg_unpaid_ch2o_day = g_leaf * (npg_dw_leaf_pot - dw_leaf) + g_stem * (npg_dw_stem_pot - dw_stem) + &
                          g_storage * (npg_dw_stor_pot - dw_storage) + g_root * (npg_dw_root_pot - dw_root)

    ! Fate of the unpaid carbon: the stem water-soluble-carbohydrate reserve first
    ! (the record-101 closed form on wsc_cap_frac, the stored fructan raising the
    ! stem dry mass the ceiling is expressed on), then the explicit residual.
    IF (npg_unpaid_ch2o_day > 0.0_wp) THEN
      want   = npg_unpaid_ch2o_day / sfl_g_wsc
      s_next = stem_biomass + dw_stem + sfl_dw_wsc_day
      room   = MAX(0.0_wp, (wsc_cap_frac * s_next - W_reserve) / (1.0_wp - wsc_cap_frac))
      npg_dw_wsc_day = MIN(want, room)
      npg_resid_day  = MAX(0.0_wp, npg_unpaid_ch2o_day - sfl_g_wsc * npg_dw_wsc_day)
    END IF
    npg_unpaid_ch2o_cum = npg_unpaid_ch2o_cum + npg_unpaid_ch2o_day
    npg_wsc_cum   = npg_wsc_cum   + npg_dw_wsc_day
    npg_resid_cum = npg_resid_cum + npg_resid_day
    ! ------------------------------------------------------------------------------ !
    ! Lamina area already built today by the organ law is not scaled back, and its
    ! unpaid dry matter is NOT re-demanded through the carbon backlog: that carbon
    ! was refused for lack of NITROGEN, not of carbon. Re-demanding it lets the
    ! lamina sink take the whole assimilate on every binding day (stem and root
    ! increments zero, most of the supply leaking to the residual), a positive
    ! feedback that is an artifact of the coupling and not a measured process. The
    ! unpaid lamina dry matter is a reported ledger, "area built ahead of its
    ! nitrogen"; the nitrogen arrest of expansion that brings it to zero is driver
    ! record 122.
    ! ------------------------------------------------------------------------------ !
    npg_lamina_unpaid_day = npg_dw_leaf_pot - dw_leaf
    npg_lamina_unpaid_cum = npg_lamina_unpaid_cum + npg_lamina_unpaid_day
  END SUBROUTINE n_payment_gate_apply

  ! -------------------------------------------------------------------------------- !
  ! Driver record 122: the pre-commit per-rank increment of the lamina-area
  ! ramp. Called after update_leaf_expansion_water and before the senescence
  ! step. It mirrors the record-83 accumulator - the day clock increment of each
  ! rank, times the water factor when record 83 is on, at a nitrogen factor of
  ! one; the commit below applies the paid fraction.
  ! -------------------------------------------------------------------------------- !
  SUBROUTINE update_lamina_n_arrest_prov()
    INTEGER  :: c, x, nt, nl_ms, n_top
    REAL(wp) :: fclk, dclk, fw, a_p, a_r, a_s
    IF (lai_n_arrest == 0) RETURN
    lna_live_today = .FALSE.
    lna_inc_area_day = 0.0_wp ; lna_removed_day = 0.0_wp ; lna_removed_dm_day = 0.0_wp
    lna_shed_need_day = 0.0_wp ; lna_shed_day = 0.0_wp ; lna_pL_day = 1.0_wp
    lna_gai_before = cls_gai ; lna_gai_after = cls_gai
    lna_dinc = 0.0_wp ; lna_apot = 0.0_wp
    IF (.NOT. clr_active) RETURN
    lna_live_today = .TRUE.
    fw = 1.0_wp
    IF (leaf_expansion_water == 1) fw = lew_fw
    nl_ms = cfp_n_top_culm()
    n_top = MIN(MIN(MAX(clr_n_app, clr_n_lig), nl_ms), CLR_MAXRANK)
    DO x = 1, n_top
      a_p = cba_A_pot(x)
      fclk = 0.0_wp
      IF (a_p > 0.0_wp) fclk = MIN(1.0_wp, MAX(0.0_wp, cba_A(x) / a_p))
      dclk = MAX(0.0_wp, fclk - lna_fclk(0, x))
      lna_dinc(0, x)  = dclk * fw
      lna_apot(0, x)  = a_p
      lna_fprov(0, x) = MIN(1.0_wp, lna_f(0, x) + lna_dinc(0, x))
      lna_fclk(0, x)  = fclk
      lna_inc_area_day = lna_inc_area_day + ct_dens(0) * lna_dinc(0, x) * a_p * 1.0e-4_wp
    END DO
    DO c = 1, CT_MAXTIL
      IF (.NOT. ct_appeared(c)) CYCLE
      nt = MAX(1, NINT(ct_nl_tiller(c, nl_ms)))
      DO x = 1, MIN(nt, CLR_MAXRANK)
        CALL ct_rank_area(c, x, nt, nl_ms, a_p, a_r, a_s)
        fclk = 0.0_wp
        IF (a_p > 0.0_wp) fclk = MIN(1.0_wp, MAX(0.0_wp, a_r / a_p))
        dclk = MAX(0.0_wp, fclk - lna_fclk(c, x))
        lna_dinc(c, x)  = dclk * fw
        lna_apot(c, x)  = a_p
        lna_fprov(c, x) = MIN(1.0_wp, lna_f(c, x) + lna_dinc(c, x))
        lna_fclk(c, x)  = fclk
        lna_inc_area_day = lna_inc_area_day + ct_dens(c) * lna_dinc(c, x) * a_p * 1.0e-4_wp
      END DO
    END DO
  END SUBROUTINE update_lamina_n_arrest_prov

  ! -------------------------------------------------------------------------------- !
  ! Driver record 122, the commit: today's increment enters the realized fraction
  ! at the paid fraction p_L, and the unpaid area (1 - p_L) times the increment
  ! leaves the standing green area.
  ! -------------------------------------------------------------------------------- !
  SUBROUTINE lamina_n_arrest_commit(p_L, W_green_lam)
    REAL(wp), INTENT(IN) :: p_L
    REAL(wp), INTENT(IN) :: W_green_lam   ! rec 179: green-lamina DM, g m-2; READ ONLY, used by no statement when rec 179 = 0
    INTEGER  :: c, x
    REAL(wp) :: rem, pl
    IF (lai_n_arrest == 0 .OR. .NOT. lna_live_today) RETURN
    pl = MIN(1.0_wp, MAX(0.0_wp, p_L))
    ! ------------------------------------------------------------------------------ !
    ! Driver record 179: a nitrogen limit on the lamina-area expansion increment.
    ! The threshold is Eq. (8) of Bertheloot, Martre & Andrieu (2008) Funct Plant
    ! Biol 35:781-796, doi:10.1071/FP08064, Table 1 and Eq. (8) both p.784, on the
    ! running maximum of the lamina's own photosynthetic nitrogen per unit area.
    ! ------------------------------------------------------------------------------ !
    IF (lamina_expansion_n_limit >= 1) THEN
      lnl_live_today  = .TRUE.
      lnl_refused_day = 0.0_wp
      lnl_fN_day      = 1.0_wp
      ! ---------------------------------------------------------------------------- !
      ! The lamina PHOTOSYNTHETIC nitrogen, the basis on which the Eq. (8) threshold
      ! is defined, and the same expression driver record 154 uses, with
      ! n2_struct_lam = 0.0079 g N g-1 DM (Bertheloot, Martre & Andrieu 2008 Plant
      ! Physiol 148:1707-1720, p.1711, doi:10.1104/pp.108.124156).
      ! ---------------------------------------------------------------------------- !
      lnl_Nph_day   = MAX(0.0_wp, npg_N_leaf - n2_struct_lam * MAX(0.0_wp, W_green_lam))
      lnl_slnph_day = 0.0_wp
      IF (cls_gai > 1.0e-6_wp) lnl_slnph_day = lnl_Nph_day / cls_gai
      ! ---------------------------------------------------------------------------- !
      ! The running maximum of Eq. (8), monotone non-decreasing by construction, with
      ! its green-area guard: below a very small green lamina area the per-area ratio
      ! is not the specific leaf nitrogen of any lamina, so the maximum is not
      ! updated there. Source of the equation: Bertheloot, Martre & Andrieu (2008)
      ! Funct Plant Biol 35:781-796, p.784.
      ! ---------------------------------------------------------------------------- !
      IF (cls_gai >= lnl_gai_max_min .AND. lnl_slnph_day > lnl_maxNpha) THEN
        lnl_maxNpha       = lnl_slnph_day
        lnl_days_maxupd   = lnl_days_maxupd + 1
        lnl_gai_at_maxupd = cls_gai
      END IF
      SELECT CASE (lamina_expansion_n_limit)
      CASE (1)   ! PRIMARY: Bertheloot 2008 FPB Eqn 8, Table 1 p.784
        lnl_minNpha_day = lnl_a_slope * MIN(lnl_maxNpha, lnl_b_max) + lnl_c_icpt
        lnl_a_sup_day   = 0.0_wp
        IF (lnl_minNpha_day > 0.0_wp) lnl_a_sup_day = lnl_Nph_day / lnl_minNpha_day
      CASE DEFAULT   ! 4 = ATTRIBUTION ARM: Yin & Schapendonk 2004 App. 2 p.424 l_N
        lnl_minNpha_day = lar_nb_leaf
        lnl_a_sup_day   = 0.0_wp
        IF (kN_canopy > 0.0_wp .AND. lar_nb_leaf > 0.0_wp) &
          lnl_a_sup_day = (1.0_wp / kN_canopy) * &
            LOG(1.0_wp + kN_canopy * MAX(0.0_wp, npg_N_leaf) / lar_nb_leaf)
      END SELECT
      lnl_allow_day = MAX(0.0_wp, lnl_a_sup_day - cls_gai)
      IF (lna_inc_area_day > 1.0e-12_wp) THEN
        lnl_fN_day = MIN(1.0_wp, lnl_allow_day / lna_inc_area_day)
        lnl_refused_day = MAX(0.0_wp, (pl - MIN(pl, lnl_fN_day)) * lna_inc_area_day)
      END IF
      pl = MIN(pl, lnl_fN_day)
      lnl_refused_cum = lnl_refused_cum + lnl_refused_day
      IF (lnl_refused_day > 1.0e-12_wp) lnl_days_refuse = lnl_days_refuse + 1
    END IF
    lna_pL_day = pl
    rem = 0.0_wp
    DO c = 0, CT_MAXTIL
      IF (c > 0 .AND. .NOT. ct_appeared(c)) CYCLE
      DO x = 1, CLR_MAXRANK
        IF (lna_dinc(c, x) <= 0.0_wp) CYCLE
        lna_f(c, x) = MIN(1.0_wp, lna_f(c, x) + lna_dinc(c, x) * pl)
        rem = rem + ct_dens(c) * lna_dinc(c, x) * (1.0_wp - pl) * lna_apot(c, x) * 1.0e-4_wp
      END DO
    END DO
    lna_gai_before = cls_gai
    rem = MIN(rem, cls_gai)                      ! cannot remove more than the standing green area
    cls_gai        = cls_gai - rem
    csl_lai_sink   = cls_gai
    csl_area_built = MAX(0.0_wp, csl_area_built - rem)
    csl_area_prev  = MAX(0.0_wp, csl_area_prev  - rem)
    lna_gai_after  = cls_gai
    lna_removed_day    = rem
    lna_removed_cum    = lna_removed_cum + rem
    lna_removed_dm_day = rem * csl_LMA
    lna_removed_dm_cum = lna_removed_dm_cum + lna_removed_dm_day
    IF (rem > 0.0_wp) lna_days_arrest = lna_days_arrest + 1
  END SUBROUTINE lamina_n_arrest_commit

  ! Driver record 122 mode 2: shed the green area the nitrogen cannot support,
  ! oldest rank first, by raising the senesced-leaf number of every live axis.
  SUBROUTINE lna_shed_by_rank(nl_ms)
    INTEGER, INTENT(IN) :: nl_ms
    INTEGER  :: c, x, nt, it
    REAL(wp) :: need, avail, dens_area, fleft, a_p, a_r, a_s, dstep, shed, per_unit(0:CT_MAXTIL), frac_left(0:CT_MAXTIL)
    IF (lai_n_arrest /= 2) RETURN

    ! Per-day diagnostics, reset here so that a day on which the criterion does not
    ! fire prints zero rather than yesterday's value.
    DO c = 0, CT_MAXTIL
      IF (c > 0 .AND. .NOT. ct_appeared(c)) CYCLE
      IF (lna_sl_bump(c) <= 0.0_wp) CYCLE
      nt = MERGE(MIN(MIN(MAX(clr_n_app, clr_n_lig), nl_ms), CLR_MAXRANK), MAX(1, NINT(ct_nl_tiller(c, nl_ms))), c == 0)
      cls_sl(c) = MIN(REAL(nt, wp), cls_sl(c) + lna_sl_bump(c))
    END DO
    lna_sln_green_in = 0.0_wp
    need = 0.0_wp
    IF (cls_gai > 1.0e-6_wp) THEN
      lna_sln_green_in = npg_N_leaf / cls_gai
      ! ---------------------------------------------------------------------------- !
      ! Driver record 147: the area a lamina nitrogen pool can retain. Mode 0 takes
      ! the incumbent minimum-specific-leaf-nitrogen criterion. The two measured
      ! values the criterion rests on are SLN_min_N = 0.50 g N m-2, the
      ! end-of-grain-filling structural value of Bertheloot, Martre & Andrieu (2008)
      ! Plant Physiol 148:1707-1720, p.1711, and the prescribed canopy value
      ! 1.5 g N m-2 of Jamieson & Semenov (2000) Field Crops Res 68:21-29, p.23,
      ! doi:10.1016/S0378-4290(00)00103-9, whose Appendix supplies the form.
      ! ---------------------------------------------------------------------------- !
      IF (SLN_min_N > 0.0_wp) THEN
        need = MAX(0.0_wp, cls_gai - npg_N_leaf / SLN_min_N)
      END IF
    END IF
    ! ------------------------------------------------------------------------------ !
    ! Tier 1 of the shedding order (driver record 134): the LOCAL specific-leaf-
    ! nitrogen shortfall remains the mechanistic driver of which rank is
    ! unsupportable, and of the order in which ranks are shed.
    ! ------------------------------------------------------------------------------ !
    lna_shed_need_day = need
    shed = 0.0_wp
    IF (need <= 0.0_wp) RETURN
    DO it = 1, 60
      avail = 0.0_wp ; dens_area = 0.0_wp ; per_unit = 0.0_wp ; frac_left = 0.0_wp
      DO c = 0, CT_MAXTIL
        IF (c > 0 .AND. .NOT. ct_appeared(c)) CYCLE
        IF (ct_dens(c) <= 0.0_wp) CYCLE
        nt = MERGE(MIN(MIN(MAX(clr_n_app, clr_n_lig), nl_ms), CLR_MAXRANK), MAX(1, NINT(ct_nl_tiller(c, nl_ms))), c == 0)
        nt = MIN(nt, CLR_MAXRANK)
        x = INT(cls_sl(c)) + 1                   ! the oldest rank not yet fully dead
        IF (x > nt) CYCLE
        IF (c == 0) THEN
          a_p = cba_A_pot(x) ; a_r = cba_A(x)
        ELSE
          CALL ct_rank_area(c, x, nt, nl_ms, a_p, a_r, a_s)
        END IF
        IF (leaf_expansion_water == 1) a_r = lew_f(c, x) * a_p
        a_r = lna_fprov(c, x) * a_p
        fleft = REAL(x, wp) - cls_sl(c)          ! the green fraction of that rank, in (0, 1]
        IF (fleft <= 0.0_wp .OR. a_r <= 0.0_wp) CYCLE
        per_unit(c)  = ct_dens(c) * a_r * 1.0e-4_wp
        frac_left(c) = fleft
        avail     = avail     + per_unit(c) * fleft
        dens_area = dens_area + per_unit(c)
      END DO
      IF (dens_area <= 0.0_wp) EXIT
      IF (need >= avail) THEN
        DO c = 0, CT_MAXTIL
          IF (per_unit(c) <= 0.0_wp) CYCLE
          cls_sl(c)      = cls_sl(c)      + frac_left(c)
          lna_sl_bump(c) = lna_sl_bump(c) + frac_left(c)
        END DO
        shed = shed + avail ; need = need - avail
        IF (need <= 1.0e-12_wp) EXIT
      ELSE
        dstep = need / dens_area
        DO c = 0, CT_MAXTIL
          IF (per_unit(c) <= 0.0_wp) CYCLE
          cls_sl(c)      = cls_sl(c)      + dstep
          lna_sl_bump(c) = lna_sl_bump(c) + dstep
        END DO
        shed = shed + need ; need = 0.0_wp
        EXIT
      END IF
    END DO
    lna_shed_day = shed
    lna_shed_cum = lna_shed_cum + shed
    IF (shed > 0.0_wp) lna_days_shed = lna_days_shed + 1
    IF (need > 1.0e-9_wp) lna_shed_unmet_days = lna_shed_unmet_days + 1
    ! Driver record 147 delivery diagnostics. Here `need` is what the 60-iteration
    ! rank loop could NOT place, so a persistent gap means the delivery and not the
    ! retention criterion is the limit.
  END SUBROUTINE lna_shed_by_rank

  SUBROUTINE update_root_depth(dvs2, root_max_depth, depth_limit, root_depth)
    REAL(wp), INTENT(IN) :: dvs2, root_max_depth, depth_limit
    REAL(wp), INTENT(OUT) :: root_depth

    REAL(wp) :: dvs_frac, root_arg

    ! Root-front descent on the sinusoidal model of Borg & Grimes (1986), driven by
    ! the thermal-time development stage (dvs2: 0 at emergence, 1 at anthesis, 2 at
    ! maturity). Driver record 124 replaces it by a measured constant descent rate.
    dvs_frac = MAX(0.0_wp, MIN(1.0_wp, dvs2 / 2.0_wp))
    root_arg = 3.03_wp * dvs_frac - 1.47_wp
    root_depth = root_max_depth * (0.5_wp + 0.5_wp * SIN(root_arg))
    root_depth = MAX(root_depth, 0.1_wp)
    root_depth = MIN(root_depth, depth_limit)
    ! Driver record 124: linear descent from emergence at a measured rate,
    ! replacing the season envelope; the envelope value is kept as a diagnostic.
    IF (root_front_law >= 1) THEN
      rfl_z_legacy = root_depth
      rfl_rz       = rfl_rz_mode(MIN(3, MAX(1, root_front_law)))
      IF (dvs2 > 0.0_wp) rfl_days_emerged = rfl_days_emerged + 1
      rfl_z_max  = MIN(root_max_depth, depth_limit)
      root_depth = MIN(rfl_z_max, MAX(rfl_z0, rfl_z0 + rfl_rz * 1.0e-3_wp * REAL(rfl_days_emerged, wp)))
      rfl_z_new  = root_depth
    END IF
  END SUBROUTINE update_root_depth

  SUBROUTINE update_structure(sla_value, green_biomass, stem_biomass, plant_density, &
                              phen_stage, prev_height, lai, plant_height_state)
    REAL(wp), INTENT(IN) :: sla_value, green_biomass, stem_biomass, plant_density, prev_height
    INTEGER, INTENT(IN) :: phen_stage
    REAL(wp), INTENT(OUT) :: lai, plant_height_state

    REAL(wp) :: stem_biomass_g_plant, height_cm

    ! Leaf-area update. The incumbent identity is LAI = green biomass x specific
    ! leaf area; under the organ law (driver record 69) the leaf area is instead
    ! the organ-level green blade area and the specific leaf area becomes an
    ! output.
    IF (canopy_sink_lai == 1 .AND. csl_active) THEN
      lai = csl_lai_sink
      csl_sla_real = 0.0_wp
      IF (green_biomass > 1.0e-9_wp) csl_sla_real = lai / green_biomass * 1000.0_wp
      lfo_sla_out = csl_sla_real    ! (rec 84): SLA as an OUTPUT column (identical quantity)

      ! The nitrogen-capacity leaf-area basis is kept explicitly separate from the
      ! sink leaf area and printed beside it, so that the two cannot diverge
      ! silently.
      csl_lai_Ncap = green_biomass * sla_value
    ELSE
      lai = green_biomass * sla_value
    END IF
    plant_height_state = prev_height

    ! Height is held at zero before emergence. This stage gate is a model
    ! integration choice for the coupling with phenology, not a claim from
    ! Bakhshandeh et al. (2012).
    IF (phen_stage < stage_emergence) THEN
      plant_height_state = 0.0_wp
      RETURN
    END IF

    ! Height is frozen after the anthesis stage: the pooled plant-height /
    ! stem-dry-weight fit is reported valid up to the anthesis date (Bakhshandeh
    ! et al. 2012, Cereal Res Commun 40(3):413-422, folio 419,
    ! doi:10.1556/CRC.40.2012.3.10, printed on folio 413).
    IF (phen_stage > stage_anthesis) RETURN

    IF (plant_density <= 0.0_wp) THEN
      plant_height_state = 0.0_wp
      RETURN
    END IF

    ! Convert stem biomass from the ground-area basis (g m-2) to g plant-1 before
    ! applying the pooled segmented fit of Bakhshandeh et al. (2012), Eq. (6)
    ! p.419, Fig. 4.
    stem_biomass_g_plant = MAX(0.0_wp, stem_biomass) / plant_density
    IF (stem_biomass_g_plant < stem_height_break_g_plant) THEN
      height_cm = stem_height_slope1_cm_g * stem_biomass_g_plant
    ELSE
      height_cm = stem_height_slope1_cm_g * stem_height_break_g_plant + &
                  stem_height_slope2_cm_g * (stem_biomass_g_plant - stem_height_break_g_plant)
    END IF

    ! Height is constrained non-decreasing here, which is a model integration
    ! choice and not part of the fitted equation. The predictor x is the TOTAL
    ! oven-dry stem dry weight (Bakhshandeh et al. 2012, folio 415).
    plant_height_state = MAX(prev_height, height_cm / 100.0_wp)
  END SUBROUTINE update_structure

  SUBROUTINE apply_leaf_senescence(dvs2, dvs_rate, sla_value, lai, green_biomass, dead_biomass, N_leaf, N_dead_leaf)
    REAL(wp), INTENT(IN)    :: dvs2, dvs_rate, sla_value, lai
    REAL(wp), INTENT(INOUT) :: green_biomass, dead_biomass, N_leaf, N_dead_leaf

    REAL(wp) :: N_CV, L_N, d_N_leaves, d_age, d_shade, d_leaves_frac
    REAL(wp) :: leafN_conc, senesced_N
    LOGICAL  :: organ_clock   ! rec 84: the organ law is the one leaf-death clock
    ! ------------------------------------------------------------------------------ !
    ! Nitrogen-driven leaf-area terms, with the constants read from
    ! input/driver_N.csv into SHARED (k_N_can, n_b_can). Yin, Schapendonk, Kropff,
    ! van Oijen & Bindraban (2004) NJAS 51:407-427, Appendix 2 (symbols L_N, k_N
    ! and n_b) supply the form; Hirose & Werger (1987) Oecologia 72:520-526,
    ! Eq. (2), is the original derivation.
    ! 1. Nitrogen-driven senescence, Yin et al. (2004) Appendix 2, the r_S,LV form.
    ! ------------------------------------------------------------------------------ !
    d_N_leaves = 0.0_wp
    organ_clock = (lai_from_organs == 1 .AND. csl_active)
    lfo_sen_dm_day = 0.0_wp; lfo_sen_N_out_day = 0.0_wp; lfo_sen_N_dead_day = 0.0_wp; lfo_sen_N_resorb_day = 0.0_wp
    ! ------------------------------------------------------------------------------ !
    ! Under the organ law (driver record 84) the leaf area is the organ-level green
    ! blade area, so the leaf MASS and the leaf NITROGEN must die on the same clock
    ! as the AREA.
    ! ------------------------------------------------------------------------------ !
    IF (organ_clock) THEN
      IF (green_biomass > 1.0e-9_wp .AND. lfo_fdied_day > 0.0_wp) &
        d_N_leaves = MIN(green_biomass, green_biomass * MIN(1.0_wp, lfo_fdied_day))
    ELSE IF (N_leaf > 0.0_wp .AND. lai > 1.0e-6_wp .AND. green_biomass > 1.0e-6_wp) THEN
      N_CV = MAX(0.0_wp, N_leaf - n_b_can * lai)
      IF (N_CV > 0.0_wp) THEN
        L_N = (1.0_wp / k_N_can) * LOG(1.0_wp + k_N_can * N_CV / n_b_can)
      ELSE
        L_N = 0.0_wp
      END IF
      L_N = MAX(0.0_wp, MIN(lai, L_N))
      ! Shed biomass proportional to excess LAI above N-supported level.
      IF (lai > L_N) THEN
        d_N_leaves = (lai - L_N) / MAX(1.0e-9_wp, sla_value)
        d_N_leaves = MIN(d_N_leaves, green_biomass)
      END IF
    END IF
    ! ------------------------------------------------------------------------------ !
    ! 2. Age-driven and self-shading senescence, Teh (2006) Introduction to
    !    Mathematical Modeling of Crop Growth, Listing 8.7 pp.193-194.
    ! ------------------------------------------------------------------------------ !
    d_shade = 0.0_wp
    IF (lai > 4.0_wp) d_shade = MIN(0.03_wp, 0.03_wp * (lai - 4.0_wp) / 4.0_wp)

    d_age = 0.0_wp
    IF (dvs2 > 1.0_wp) d_age = dvs_rate / MAX(0.1_wp, 2.0_wp - dvs2)

    d_leaves_frac = MAX(d_age, d_shade)
    ! Teh (2006) fractional-rate leaves (proportional to green biomass).
    IF (.NOT. organ_clock) d_N_leaves = d_N_leaves + green_biomass * d_leaves_frac
    ! ------------------------------------------------------------------------------ !
    ! 3. Apply the senescence: transfer green to dead with PRE-ABSCISSION nitrogen
    !    resorption. The exported fraction is measured by 15N labeling in Kichey,
    !    Hirel, Heumez, Dubois & Le Gouis (2007) New Phytol 176:625-638, Table 2,
    !    and the residual that stays with the dead tissue is the structural lamina
    !    nitrogen of Bertheloot, Martre & Andrieu (2008) Plant Physiol
    !    148:1707-1720, doi:10.1104/pp.108.124156, Results p.1711, where the lamina
    !    nitrogen concentration falls to a minimum that is "assumed to be
    !    structural N"; with their Fig. 2C p.1711 (canopy lamina dry mass per unit
    !    area) this gives the derived span 0.41-0.51 g N m-2, of which the model
    !    uses 0.50.
    ! ------------------------------------------------------------------------------ !
    IF (d_N_leaves > 0.0_wp .AND. green_biomass > 1.0e-9_wp) THEN
      d_N_leaves = MIN(d_N_leaves, green_biomass)
      leafN_conc = N_leaf / MAX(1.0e-9_wp, green_biomass)
      BLOCK
        REAL(wp) :: N_carried, dA_sen, N_to_dead
        N_carried = d_N_leaves * leafN_conc
        dA_sen    = d_N_leaves * MAX(0.0_wp, sla_value)
        IF (organ_clock) dA_sen = MAX(0.0_wp, lfo_gai_died_day)   ! rec 84: the organ law's own dead area
        N_to_dead = MIN(N_carried, SLN_min_N * dA_sen)
        senesced_N = N_to_dead
      END BLOCK
      dead_biomass = dead_biomass + d_N_leaves
      green_biomass = MAX(0.0_wp, green_biomass - d_N_leaves)
      N_dead_leaf = N_dead_leaf + senesced_N
      ! ---------------------------------------------------------------------------- !
      ! The FULL proportional nitrogen leaves the leaf pool with the dying tissue;
      ! the resorbed, above-structural share goes to the MOBILE VEGETATIVE pool,
      ! which is the culm-common pool of Bertheloot, Martre & Andrieu (2008) Funct
      ! Plant Biol 35:781-796, pp.782-783, and not to the grain directly.
      ! ---------------------------------------------------------------------------- !
      BLOCK
        REAL(wp) :: N_carried_out
        N_carried_out = MIN(N_leaf, d_N_leaves * leafN_conc)
        N_leaf = MAX(0.0_wp, N_leaf - N_carried_out)
        N_mobile_veg = N_mobile_veg + MAX(0.0_wp, N_carried_out - senesced_N)
        IF (organ_clock) THEN   ! rec 84 diagnostics (daily_d2_organs.csv), nothing acts on them
          lfo_sen_dm_day = d_N_leaves; lfo_sen_N_out_day = N_carried_out; lfo_sen_N_dead_day = senesced_N
          lfo_sen_N_resorb_day = MAX(0.0_wp, N_carried_out - senesced_N)
          lfo_sen_dm_cum       = lfo_sen_dm_cum       + lfo_sen_dm_day
          lfo_sen_N_out_cum    = lfo_sen_N_out_cum    + lfo_sen_N_out_day
          lfo_sen_N_dead_cum   = lfo_sen_N_dead_cum   + lfo_sen_N_dead_day
          lfo_sen_N_resorb_cum = lfo_sen_N_resorb_cum + lfo_sen_N_resorb_day
        END IF

        ! Driver record 147 ledger: the nitrogen this day senescence released to the
        ! mobile vegetative pool, accumulated over the days the retention criterion
        ! was active.
      END BLOCK
    END IF
  END SUBROUTINE apply_leaf_senescence
    ! ------------------------------------------------------------------------------ !
  ! Post-flowering organ-nitrogen turnover and tissue death for the stem and
  ! chaff. Sources: Bertheloot, Martre & Andrieu (2008) Funct Plant Biol
  ! 35:781-796, doi:10.1071/FP08064, Eqs. (5), (6), (7) and (13), pp.784-785;
  ! Bertheloot, Martre & Andrieu (2008) Plant Physiol 148:1707-1720,
  ! doi:10.1104/pp.108.124156, p.1711 (the structural concentrations) and p.1712
  ! (the degradation rate); and Bertheloot, Wu, Cournede & Andrieu (2011) Ann Bot
  ! 108:1097-1109, doi:10.1093/aob/mcr125, Tables 1 and 2 (the organ-resolved
  ! rate constants). Per organ and per day, on the day thermal time dTT:
  !   N_ph   = MAX(0, N_organ - c_struct M_organ)      photosynthetic nitrogen
  !   D      = delta N_ph dTT                          degradation -> mobile
  !   S      = sigma M_organ f_N f_PAR dTT             synthesis <- mobile
  !   f_N    = [N]mob / ([N]mob + k1)                  Eq. (7) p.784
  !   f_PAR  = PAR / (PAR + k2)
  ! and the living mass is reduced so that the organ nitrogen per unit
  ! photosynthetic area stays at its flowering reference, the surplus structural
  ! nitrogen going to the dead pool.
    ! ------------------------------------------------------------------------------ !
  SUBROUTINE apply_organN_turnover_daily(dvs2, dTT, W_green, W_stem, W_chaff, W_storage, &
                                         N_leaf, N_stem, N_chaff, N_storage)
    REAL(wp), INTENT(IN)    :: dvs2, dTT, W_green
    REAL(wp), INTENT(INOUT) :: W_stem, W_chaff, W_storage
    REAL(wp), INTENT(INOUT) :: N_leaf, N_stem, N_chaff, N_storage
    REAL(wp) :: Mg_tot, conc_mob, f_N, PAR_day, f_PAR
    REAL(wp) :: nph_lam, nph_stem, nph_chaff
    REAL(wp) :: d_lam, d_stem, d_chaff, s_lam, s_stem, s_chaff, s_tot, s_scale

    n2_deg_lam_day = 0.0_wp; n2_deg_stem_day = 0.0_wp; n2_deg_chaff_day = 0.0_wp
    n2_syn_lam_day = 0.0_wp; n2_syn_stem_day = 0.0_wp; n2_syn_chaff_day = 0.0_wp
    n2_dW_dead_stem_day = 0.0_wp; n2_dW_dead_chaff_day = 0.0_wp
    IF (dvs2 < 1.0_wp) RETURN
    IF (dTT <= 0.0_wp) RETURN

    ! --- flowering snapshot of the reference nitrogen per unit photosynthetic
    !     area, taken on the first post-flowering call.
    IF (.NOT. n2_anth_snap) THEN
      n2_anth_snap = .TRUE.
      IF (W_stem  > 1.0e-9_wp) &
        n2_cph_stem_anth  = MAX(0.0_wp, N_stem  - n2_struct_stem * W_stem)  / W_stem
      IF (W_chaff > 1.0e-9_wp) &
        n2_cph_chaff_anth = MAX(0.0_wp, N_chaff - MAX(n2_struct_chaff_use * W_chaff, cnr_floor_gN)) / W_chaff
    END IF

    ! --- the two saturating factors of the xylem-resupply law (FPB Eqn 6).
    Mg_tot = MAX(0.0_wp, W_green) + MAX(0.0_wp, W_stem) + MAX(0.0_wp, W_chaff)
    conc_mob = 0.0_wp
    IF (Mg_tot > 1.0e-9_wp) conc_mob = MAX(0.0_wp, N_mobile_veg) / Mg_tot
    f_N   = conc_mob / (conc_mob + n2_k1_Nmob)
    PAR_day = n2_par_frac * MAX(0.0_wp, dTotrad_s) * 1.0e6_wp
    f_PAR = PAR_day / (PAR_day + n2_k2_PAR)

    ! --- photosynthetic-N state of each group (total N minus structural N).
    nph_lam   = MAX(0.0_wp, N_leaf  - n2_struct_lam  * MAX(0.0_wp, W_green))
    nph_stem  = MAX(0.0_wp, N_stem  - n2_struct_stem * MAX(0.0_wp, W_stem))
    nph_chaff = MAX(0.0_wp, N_chaff - MAX(n2_struct_chaff_use * MAX(0.0_wp, W_chaff), cnr_floor_gN))

    ! --- degradation (FPB Eqn 5), uniform delta, credited to the mobile pool.
    d_lam   = MIN(nph_lam,   n2_delta_TT * nph_lam   * dTT)
    d_stem  = MIN(nph_stem,  n2_delta_TT * nph_stem  * dTT)
    d_chaff = MIN(nph_chaff, n2_delta_TT * nph_chaff * dTT)
    N_leaf  = N_leaf  - d_lam
    N_stem  = N_stem  - d_stem
    N_chaff = N_chaff - d_chaff
    N_storage = N_storage - d_chaff          ! chaff N is a tagged subset of N_storage
    N_mobile_veg = N_mobile_veg + d_lam + d_stem + d_chaff
    n2_deg_lam_day = d_lam; n2_deg_stem_day = d_stem; n2_deg_chaff_day = d_chaff
    n2_cum_deg_lam   = n2_cum_deg_lam   + d_lam
    n2_cum_deg_stem  = n2_cum_deg_stem  + d_stem
    n2_cum_deg_chaff = n2_cum_deg_chaff + d_chaff

    ! --- xylem resupply (FPB Eqn 6), organ-type sigma, drawn from the pool.
    s_lam   = n2_sigma_lam   * MAX(0.0_wp, W_green) * f_N * f_PAR * dTT
    s_stem  = n2_sigma_stem  * MAX(0.0_wp, W_stem)  * f_N * f_PAR * dTT
    s_chaff = n2_sigma_chaff * MAX(0.0_wp, W_chaff) * f_N * f_PAR * dTT
    s_tot   = s_lam + s_stem + s_chaff
    IF (s_tot > 0.0_wp) THEN
      s_scale = MIN(1.0_wp, MAX(0.0_wp, N_mobile_veg) / s_tot)   ! supply limit, pro rata
      s_lam = s_lam * s_scale; s_stem = s_stem * s_scale; s_chaff = s_chaff * s_scale
      N_leaf  = N_leaf  + s_lam
      N_stem  = N_stem  + s_stem
      N_chaff = N_chaff + s_chaff
      N_storage = N_storage + s_chaff
      N_mobile_veg = MAX(0.0_wp, N_mobile_veg - (s_lam + s_stem + s_chaff))
      n2_syn_lam_day = s_lam; n2_syn_stem_day = s_stem; n2_syn_chaff_day = s_chaff
    END IF

    ! --- Bertheloot 2011a tissue death, stem group (piece a).
    IF (W_stem > 1.0e-9_wp .AND. n2_cph_stem_anth > 1.0e-12_wp) THEN
      BLOCK
        REAL(wp) :: nph, c_ph, thr, m_new, dM, n_to_dead
        nph  = MAX(0.0_wp, N_stem - n2_struct_stem * W_stem)
        c_ph = nph / W_stem
        thr  = n2_d_tp * n2_cph_stem_anth
        IF (c_ph < thr) THEN
          m_new = nph / thr
          dM    = MAX(0.0_wp, W_stem - m_new)
          n_to_dead = MIN(MAX(0.0_wp, N_stem), n2_struct_stem * dM)
          W_stem      = W_stem - dM
          N_stem      = N_stem - n_to_dead
          W_dead_stem = W_dead_stem + dM
          N_dead_stem = N_dead_stem + n_to_dead
          n2_dW_dead_stem_day = dM
        END IF
      END BLOCK
    END IF

    ! --- tissue death, chaff group. W_chaff is a subset of W_storage and N_chaff
    !     of N_storage, so both are decremented.
    IF (W_chaff > 1.0e-9_wp .AND. n2_cph_chaff_anth > 1.0e-12_wp) THEN
      BLOCK
        REAL(wp) :: nph, c_ph, thr, m_new, dM, n_to_dead
        nph  = MAX(0.0_wp, N_chaff - MAX(n2_struct_chaff_use * W_chaff, cnr_floor_gN))
        c_ph = nph / W_chaff
        thr  = n2_d_tp * n2_cph_chaff_anth
        IF (c_ph < thr) THEN
          m_new = nph / thr
          dM    = MAX(0.0_wp, W_chaff - m_new)
          n_to_dead = MIN(MAX(0.0_wp, N_chaff), n2_struct_chaff_use * dM)
          W_chaff   = W_chaff   - dM
          W_storage = W_storage - dM
          N_chaff   = N_chaff   - n_to_dead
          N_storage = N_storage - n_to_dead
          W_dead_stem = W_dead_stem + dM
          N_dead_stem = N_dead_stem + n_to_dead
          n2_dW_dead_chaff_day = dM
        END IF
      END BLOCK
    END IF
  END SUBROUTINE apply_organN_turnover_daily
  ! -------------------------------------------------------------------------------- !
  ! Grain nitrogen source-sink. Martre, Porter, Jamieson & Triboi (2003) Plant
  ! Physiol 133:1959-1967, doi:10.1104/pp.103.030585, measure grain nitrogen as
  ! SOURCE-limited during storage-protein accumulation (their Fig. 5: supply two
  ! to three times demand after about 15 days; p.1963, grain nitrogen is
  ! "regulated by the source of N and not by the activity of" the grain), so the
  ! grain draws on the vegetative pools rather than setting its own rate. The
  ! organ floors are the structural concentrations of Bertheloot, Martre &
  ! Andrieu (2008) Plant Physiol 148:1707-1720, p.1711.
  ! -------------------------------------------------------------------------------- !
  SUBROUTINE allocate_grain_N_source_sink_daily(dvs2, dW_storage_day, N_assim_avail, &
               lai, W_stem, W_green, W_chaff_gN, N_leaf, N_stem, N_storage, &
               grain_import_day, N_from_assim_day, N_remob_leaf_day, N_remob_stem_day, N_surplus_veg_day)
    REAL(wp), INTENT(IN)    :: dvs2, dW_storage_day, N_assim_avail, lai, W_stem
    REAL(wp), INTENT(IN)    :: W_green   ! g m-2 green leaf DM (same-day LAI identity basis)
    ! Chaff dry mass, needed by the chaff limb of the remobilization-efficiency
    ! floor.
    REAL(wp), INTENT(IN)    :: W_chaff_gN
    REAL(wp), INTENT(INOUT) :: N_leaf, N_stem, N_storage
    REAL(wp), INTENT(OUT)   :: grain_import_day, N_from_assim_day
    REAL(wp), INTENT(OUT)   :: N_remob_leaf_day, N_remob_stem_day, N_surplus_veg_day
    REAL(wp) :: source_tot, stru_demand
    REAL(wp) :: take_assim
    REAL(wp) :: lai_basis
    REAL(wp) :: active_cap_N   ! rec 132: sln_max_N, or sln_max_N_lnh under mode 2

    grain_import_day = 0.0_wp; N_from_assim_day = 0.0_wp
    N_remob_leaf_day = 0.0_wp; N_remob_stem_day = 0.0_wp; N_surplus_veg_day = 0.0_wp
    N_from_mobile_day = 0.0_wp
    IF (dvs2 < 1.0_wp) RETURN                                      ! grain fill only (post-anthesis)
    ! ------------------------------------------------------------------------------ !
    ! The leaf-nitrogen holding capacity applied in this routine. The capacity is
    ! sln_max_N = 3.8 g N m-2, the maximum measured field-wheat Narea of
    ! Silva-Perez et al. (2020) J Exp Bot 71(7):2299-2311,
    ! doi:10.1093/jxb/erz439, Fig. 4B p.2305. It is a holding capacity per unit
    ! AREA and never a concentration clip on an output.
    ! ------------------------------------------------------------------------------ !
    lai_basis = 0.0_wp
    ! ------------------------------------------------------------------------------ !
    ! Driver record 132: the ceiling applied throughout the general holding-
    ! capacity ladder of this routine and of route_surplus_bounded. Mode 2 uses the
    ! printed ceiling 2.5 g N m-2 of Pask (2009) PhD thesis, University of
    ! Nottingham, p.29.
    ! ------------------------------------------------------------------------------ !
    active_cap_N = sln_max_N
    ! ------------------------------------------------------------------------------ !
    ! Driver record 133 mode 2 re-selects the SAME parameter, so no new constant is
    ! introduced; records 132 and 133 are mutually exclusive and that is enforced
    ! at the driver read.
    ! ------------------------------------------------------------------------------ !
    lnh_cap_in_force = active_cap_N
    lto_cap_in_force = active_cap_N
    lnc_excess_day = 0.0_wp; lnc_to_mobile_day = 0.0_wp; lnc_lai_basis = 0.0_wp
    ! ------------------------------------------------------------------------------ !
    ! Driver record 87: the capacity basis is the leaf area index that the
    ! gas-exchange routine actually read, so that the capacity and the
    ! photosynthetic capacity are expressed on the same area.
    ! ------------------------------------------------------------------------------ !
    IF (leafN_drawdown == 1 .OR. leafN_cap_route == 1) THEN
      IF (leafN_cap_route == 1) THEN
        lai_basis = MAX(0.0_wp, lai)
      ELSE
        lai_basis = MAX(0.0_wp, W_green) * sla_from_dvs(dvs2)
      END IF
      lnc_lai_basis = lai_basis
      BLOCK
        REAL(wp) :: excess
        excess = MAX(0.0_wp, N_leaf - active_cap_N * lai_basis)
        N_leaf = N_leaf - excess
        N_mobile_veg = N_mobile_veg + excess
        lnc_excess_day = excess
      END BLOCK
    END IF
    ! ------------------------------------------------------------------------------ !
    ! The grain gets a sink (driver record 86). Source: Bertheloot, Martre &
    ! Andrieu (2008) Funct Plant Biol 35:781-796, doi:10.1071/FP08064, Eq. (11) and
    ! Eq. (12) p.785 and Table 1 p.784, on thermal time with base 0 C (p.782),
    ! which is the model's own daily nitrogen thermal time:
    !   dN_g,tot    = MIN(S_g,tot_pot, N_c,mob)                     Eq. (11)
    !   S_g,tot_pot = gamma N_g,tot_pot  up to the plateau time tau,
    !                 constant thereafter                           Eq. (12)
    ! and the potential grain nitrogen integrates to
    !   n_g = n_g0 exp(g tau) (1 + g (tt - tau))   for tt > tau.
    ! The lamina synthesis term of Eq. (6) p.784 is restored alongside it, so that
    ! the grain claim and the lamina claim compete for the same mobile pool.
    ! ------------------------------------------------------------------------------ !
    IF (grainN_sink_mech >= 1 .AND. grainN_sink_mech /= 3) THEN
      BLOCK
        REAL(wp) :: sink, take_assim_f, draw_mob, short
        REAL(wp) :: pot_lam, pot_stem, pot_chaff, dtt_f
        REAL(wp) :: nl0, ns0, nc0, nst0, nm0, resid
        REAL(wp) :: nph_lam, d_lam, s_lam, mg_tot, conc_mob, par_day, pool_av
        REAL(wp) :: claim_tot, need_f, act_syn, grain_pool, grain_assim_f
        ! Driver record 154: root nitrogen on ENTRY, captured unconditionally. It
        ! changes no flux and no output and is READ only by the mode-1 limb of the
        ! closure below.
        REAL(wp) :: nr0_154
        INTEGER  :: dpa_f

        ! --- one-off grain-nitrogen sink initialization on the first post-anthesis
        !     call. Driver record 115 selects the basis of the initial per-kernel or
        !     per-ear grain nitrogen.
        IF (.NOT. gnf_snap) THEN
          gnf_snap = .TRUE.
            gnf_ng0  = gnf_ng0_per_ear * MAX(0.0_wp, gnx_ears_m2)
          gnf_tt   = 0.0_wp
        END IF
        dtt_f   = MAX(0.0_wp, dTT_day_N)
        pot_lam = 0.0_wp ; pot_stem = 0.0_wp ; pot_chaff = 0.0_wp
        sink = 0.0_wp ; take_assim_f = 0.0_wp ; draw_mob = 0.0_wp ; short = 0.0_wp
        nph_lam = 0.0_wp ; d_lam = 0.0_wp ; s_lam = 0.0_wp ; mg_tot = 0.0_wp
        conc_mob = 0.0_wp ; par_day = 0.0_wp ; pool_av = 0.0_wp
        claim_tot = 0.0_wp ; need_f = 0.0_wp ; act_syn = 0.0_wp ; grain_pool = 0.0_wp
        grain_assim_f = 0.0_wp
        ! Organ state on ENTRY. The closure below is read from this and never from the
        ! fluxes themselves, which would be a tautology.
        nl0 = N_leaf ; ns0 = N_stem ; nc0 = N_chaff ; nst0 = N_storage ; nm0 = N_mobile_veg
        nr0_154 = npg_N_root                                  ! rec 154
        ont_live_today = .FALSE.
        gnf_draw_lam = 0.0_wp ; gnf_draw_stem = 0.0_wp ; gnf_draw_chaff = 0.0_wp
        gnf_sink_day = 0.0_wp ; gnf_supply_day = 0.0_wp
        gnf_deg_lam_day = 0.0_wp ; gnf_syn_lam_day = 0.0_wp ; gnf_syn_pot_day = 0.0_wp
        gnf_fN_day = 0.0_wp ; gnf_fPAR_day = 0.0_wp ; gnf_pool_scale = 1.0_wp
        fnc_syn_clip_day = 0.0_wp ; fnc_mob_close_day = 0.0_wp   ! rec 97 diagnostics
        ! ------------------------------------------------------------------------------ !
        IF (dvs2 >= 2.0_wp .OR. canopy_senesced) THEN
          ! ---------------------------------------------------------------------------- !
          ! Termination: assimilated nitrogen after the end of grain filling stays
          ! vegetative on the bounded ladder.
          ! ---------------------------------------------------------------------------- !
          N_surplus_veg_day = MAX(0.0_wp, N_assim_avail)
          IF (leafN_drawdown == 1 .OR. leafN_cap_route == 1) THEN
            CALL route_surplus_bounded(N_surplus_veg_day, lai_basis, N_leaf, active_cap_N)
          ELSE
            N_leaf = N_leaf + N_surplus_veg_day
          END IF
          assimN_surplus_postfill = assimN_surplus_postfill + N_surplus_veg_day
          ! ---------------------------------------------------------------------------- !
          ! Terminal close of the mobile vegetative pool (driver record 97 modes 2 and 3).
          ! The pool of Bertheloot, Martre & Andrieu (2008) Funct Plant Biol 35:781-796,
          ! p.782, is culm-common and in phloem and vacuolar transit, so what it still
          ! holds at the end of filling cannot be left in a vegetative organ.
          ! ---------------------------------------------------------------------------- !
        ELSE
          gnf_days_fill = gnf_days_fill + 1
          ! ---------------------------------------------------------------------------- !
          ! Driver record 178 modes 2 and 3: the anthesis snapshot of the chaff nitrogen
          ! and the remobilization-efficiency floor it defines. Taken here, on the first
          ! live filling day and BEFORE this day degradation, so that it is the chaff's own
          ! anthesis state. Inactive alternative, see Part 1 Table 6.
          ! ---------------------------------------------------------------------------- !
          IF (chaff_structural_n_basis >= 2 .AND. .NOT. cnr_snap_done) THEN
            cnr_snap_done     = .TRUE.
            cnr_N_chaff_anth  = MAX(0.0_wp, N_chaff)
            cnr_floor_gN      = (1.0_wp - cnr_nre_chaff) * cnr_N_chaff_anth
          END IF
          IF (chaff_structural_n_basis >= 2) THEN
            ! which of the two lower bounds is binding today, reported only
            cnr_held_gN = MAX(0.0_wp, cnr_floor_gN &
                        - n2_struct_chaff_use * MAX(0.0_wp, W_chaff_gN))
            IF (cnr_held_gN > 0.0_wp) cnr_days_bound = cnr_days_bound + 1
          END IF
          ! ---------------------------------------------------------------------------- !
          ! Driver record 175 modes 2 and above: the initial condition of the grain
          ! nitrogen. Bertheloot, Andrieu, Fournier & Martre (2008) Funct Plant Biol
          ! 35:781-796, p.786, states the measured grain nitrogen at flowering, and
          ! Table 1 p.784 prints the initial grain nitrogen mass of the equation set, so
          ! the handover is made on a measured state rather than on a grain that has no
          ! mass.
          ! ---------------------------------------------------------------------------- !
          IF (grain_n_partition_basis >= 4) THEN
            BLOCK
              REAL(wp) :: g0_want, g0_owed, g0_take
              gnp_g0_target = gnf_ng0
              g0_want = MIN(gnf_ng0, n2_grain_anth_conc * MAX(0.0_wp, W_grain_state))
              g0_owed = MAX(0.0_wp, g0_want - gnp_g0_paid_cum)
              g0_take = MIN(g0_owed, &
                        MAX(0.0_wp, N_chaff - MAX(n2_struct_chaff_use * MAX(0.0_wp, W_chaff_gN), cnr_floor_gN)))
              N_chaff         = N_chaff - g0_take
              gnp_g0_paid_cum = gnp_g0_paid_cum + g0_take
              ! ------------------------------------------------------------------------ !
              ! All three quantities are STATE and not sums: the take is what the chaff 
              ! has paid in total, the refused amount is the standing arrears, and the 
              ! third is the current claim.
              ! ------------------------------------------------------------------------ !
              gnp_g0_take    = gnp_g0_paid_cum
              gnp_g0_refused = MAX(0.0_wp, g0_want - gnp_g0_paid_cum)
              gnp_g0_deficit = MAX(0.0_wp, gnf_ng0 - gnp_g0_paid_cum)
              IF (gnf_ng0 > 0.0_wp) gnp_g0_frac = gnp_g0_paid_cum / gnf_ng0
              gnp_g0_done    = gnp_g0_paid_cum >= gnf_ng0
            END BLOCK
          END IF
          ! ---------------------------------------------------------------------------- !
          ! The thermal clock and the potential advance ONLY while the filling is live.
          ! The thermal time of Bertheloot et al. (2008) is time after anthesis DURING
          ! grain filling, and letting it run past the measured termination of Martre
          ! et al. (2003) Plant Physiol 133:1959-1967, p.1963, would extend the sink
          ! beyond the measurement.
          ! ---------------------------------------------------------------------------- !
          gnf_tt = gnf_tt + dtt_f
          ! ---------------------------------------------------------------------------- !
          ! The potential grain nitrogen, Eq. (12) of Bertheloot, Martre & Andrieu (2008)
          ! Funct Plant Biol 35:781-796 p.785 integrated. It is reported, and used as the
          ! post-anthesis uptake demand in main_MMCO2Crop.f90, so that the grain-nitrogen
          ! demand is defined once.
          ! ---------------------------------------------------------------------------- !
          IF (gnf_tt <= gnf_tau_degCd) THEN
            gnf_pot_cum = gnf_ng0 * EXP(gnf_g_rate * gnf_tt)
          ELSE
            gnf_pot_cum = gnf_ng0 * EXP(gnf_g_rate * gnf_tau_degCd) &
                        * (1.0_wp + gnf_g_rate * (gnf_tt - gnf_tau_degCd))
          END IF
          sink = gnf_g_rate * gnf_ng0 * EXP(gnf_g_rate * MIN(gnf_tt, gnf_tau_degCd)) * dtt_f
          ! ---------------------------------------------------------------------------- !
          ! Driver record 125 modes 1 and 2: the per-kernel grain-nitrogen ceiling,
          ! cap = kernel number x 1.38 mg N, the widest measured control value of Martre,
          ! Porter, Jamieson & Triboi (2003) Plant Physiol 133:1959-1967, Table I p.1961.
          ! ---------------------------------------------------------------------------- !
          gnc_cut_day = 0.0_wp
          IF ((kernel_potential_law == 1 .OR. kernel_potential_law == 2) .AND. kernel_number_m2 > 0.0_wp) THEN
            ! -------------------------------------------------------------------------- !
            ! The live per-kernel ceiling is that same 1.38 mg N per kernel unless a
            ! cultivar file supplies its own value under driver record 177, in which case
            ! the value is held inside the measured 0.75-1.38 mg N window of the same
            ! table.
            ! -------------------------------------------------------------------------- !
            gnc_cap      = kernel_number_m2 * gnp_mgN_kernel * 1.0e-3_wp
            gnc_room_day = MAX(0.0_wp, gnc_cap - MAX(0.0_wp, N_storage - N_chaff))
            IF (sink > gnc_room_day) THEN
              gnc_cut_day    = sink - gnc_room_day
              gnc_cut_cum    = gnc_cut_cum + gnc_cut_day
              gnc_days_bound = gnc_days_bound + 1
              sink           = gnc_room_day
            END IF
            gnf_pot_cum = MIN(gnf_pot_cum, gnc_cap)
          END IF
          pot_stem  = n2_delta_TT * MAX(0.0_wp, N_stem  - n2_struct_stem * MAX(0.0_wp, W_stem))     * dtt_f
          pot_chaff = n2_delta_TT * MAX(0.0_wp, N_chaff - MAX(n2_struct_chaff_use * MAX(0.0_wp, W_chaff_gN), cnr_floor_gN)) * dtt_f
          gnf_sink_day   = sink
          IF (gnf_eqn6) THEN
          IF (organ_n_turnover_nema >= 1) THEN
            ! -------------------------------------------------------------------------- !
            ! Post-anthesis organ nitrogen turnover (driver record 154). The organ-resolved
            ! degradation and synthesis constants are those of Bertheloot, Wu, Cournede &
            ! Andrieu (2011) Ann Bot 108:1097-1109, doi:10.1093/aob/mcr125, Table 2, and
            ! the equations are Eqs. (5), (6) and (12) of Bertheloot, Martre & Andrieu
            ! (2008) Funct Plant Biol 35:781-796, p.784-785.
            ! -------------------------------------------------------------------------- !
            BLOCK
              REAL(wp) :: nph_stem, nph_chaff, nrem_root
              REAL(wp) :: d_stem, d_chaff, d_root
              REAL(wp) :: s_stem, s_chaff, s_root, s_tot, s_scale
              REAL(wp) :: pool_entry, uptake_in, grain_in, s_lam_pot154
              REAL(wp) :: gnp_claim_resid, gnp_grain_extra   ! rec 175 mode 3
              REAL(wp) :: lsb_m_sheath, lsb_m_lam, lsb_m_stem   ! rec 167
              LOGICAL  :: root_limb
              root_limb  = (organ_n_turnover_nema == 1)
              pool_entry = MAX(0.0_wp, N_mobile_veg)
              uptake_in  = MAX(0.0_wp, N_assim_avail)
              d_stem = 0.0_wp ; d_chaff = 0.0_wp ; d_root = 0.0_wp
              s_stem = 0.0_wp ; s_chaff = 0.0_wp ; s_root = 0.0_wp
              s_lam_pot154 = 0.0_wp ; s_scale = 1.0_wp ; grain_in = 0.0_wp
              nph_stem = 0.0_wp ; nph_chaff = 0.0_wp ; nrem_root = 0.0_wp
              ! ------------------------------------------------------------------------ !
              ! ---- (a) DEGRADATION, unconditional. Eq. (5) p.784 of Bertheloot, Martre &
              !      Andrieu (2008) Funct Plant Biol 35:781-796, which is Eq. (8) of
              !      Bertheloot, Cournede & Andrieu (2011) Ann Bot 108:1085-1096. Each flux
              !      is bounded by the photosynthetic pool it draws on.
              ! ------------------------------------------------------------------------ !
              IF (dtt_f > 0.0_wp) THEN
                nph_lam   = MAX(0.0_wp, N_leaf  - n2_struct_lam  * MAX(0.0_wp, W_green))    &
                          / (1.0_wp + gnf_r_mob_ph)             ! FPB Eqn 1 p.783, r = 0.10 p.786
                nph_stem  = MAX(0.0_wp, N_stem  - n2_struct_stem * MAX(0.0_wp, W_stem))     &
                          / (1.0_wp + gnf_r_mob_ph)
                nph_chaff = MAX(0.0_wp, N_chaff - MAX(n2_struct_chaff_use * MAX(0.0_wp, W_chaff_gN), cnr_floor_gN)) &
                          / (1.0_wp + gnf_r_mob_ph)
                d_lam   = MIN(nph_lam,   n2_delta_TT * nph_lam   * dtt_f)
                d_stem  = MIN(nph_stem,  n2_delta_TT * nph_stem  * dtt_f)
                d_chaff = MIN(nph_chaff, n2_delta_TT * nph_chaff * dtt_f)
                N_leaf       = N_leaf    - d_lam
                N_stem       = N_stem    - d_stem
                N_chaff      = N_chaff   - d_chaff
                N_storage    = N_storage - d_chaff        ! chaff N is a tagged subset of N_storage
                N_mobile_veg = N_mobile_veg + d_lam + d_stem + d_chaff
                IF (root_limb) THEN
                  nrem_root = MAX(0.0_wp, npg_N_root - npg_floor_R * MAX(0.0_wp, ont_W_root_in))
                  d_root    = MIN(nrem_root, n2_delta_TT * nrem_root * dtt_f)
                  npg_N_root   = npg_N_root   - d_root
                  N_mobile_veg = N_mobile_veg + d_root
                END IF
              END IF
              gnf_deg_lam_day = d_lam
              gnf_cum_deg_lam = gnf_cum_deg_lam + d_lam
              ! ---- (b) UPTAKE CREDIT. With ONE pool, today's reduced nitrogen enters it and
              !      is spent from it, so there is no assimilation-first ladder and the
              !      vegetative surplus of the day is exactly zero on this path.
              N_mobile_veg   = N_mobile_veg + uptake_in
              gnf_supply_day = MAX(0.0_wp, N_mobile_veg)   ! the supply the grain sees

              ! ---- (c) GRAIN. The sink is the Eq. (12) potential rate including the
              !      per-kernel cut of driver record 125, both computed above.
              grain_in     = MIN(MAX(0.0_wp, sink), MAX(0.0_wp, N_mobile_veg))
              N_mobile_veg = N_mobile_veg - grain_in
              IF (sink > 1.0e-12_wp) THEN
                gnf_pool_scale = MIN(1.0_wp, grain_in / sink)
              ELSE
                gnf_pool_scale = 1.0_wp
              END IF
              ont_bound_by = MERGE(1, 0, grain_in < sink - 1.0e-12_wp)
              ! ------------------------------------------------------------------------ !
              ! ---- (d) SYNTHESIS from the REMAINING pool. Eq. (6) p.784 of Bertheloot,
              !      Martre & Andrieu (2008) Funct Plant Biol 35:781-796, which is Eq. (12)
              !      of Bertheloot, Cournede & Andrieu (2011) Ann Bot 108:1085-1096, with the
              !      mobile-nitrogen factor on the Eq. (7) p.784 basis and the light factor
              !      on the DAILY PAR integral, the basis on which Bertheloot, Wu, Cournede &
              !      Andrieu (2011) Ann Bot 108:1097-1109, Table 2, prints its constant.
              ! ------------------------------------------------------------------------ !
              mg_tot   = MAX(0.0_wp, W_green) + MAX(0.0_wp, W_stem)
              conc_mob = 0.0_wp
              IF (mg_tot > 1.0e-9_wp) conc_mob = MAX(0.0_wp, N_mobile_veg) / mg_tot
              gnf_fN_day = conc_mob / (conc_mob + n2_k1_Nmob)
              par_day    = n2_par_frac * MAX(0.0_wp, dTotrad_s) * 1.0e6_wp
              canE_fPAR_top = par_day / (par_day + n2_k2_PAR)
              canE_fPAR_eff = canE_fPAR_top
              IF (lai > 1.0e-6_wp .AND. par_day > 0.0_wp) &
                canE_fPAR_eff = LOG( (par_day + n2_k2_PAR) &
                                   / (par_day * EXP(-kL_bertheloot * lai) + n2_k2_PAR) ) &
                              / (kL_bertheloot * lai)
              IF (canopy_n_flag_zone == 3 .OR. canopy_n_flag_zone == 4) THEN
                gnf_fPAR_day = canE_fPAR_eff
              ELSE
                gnf_fPAR_day = canE_fPAR_top
              END IF
              IF (dtt_f > 0.0_wp) THEN

                ! ---- Driver record 167. Table 2 p.1101 of Bertheloot, Wu, Cournede & Andrieu
                !      (2011) Ann Bot 108:1097-1109 assigns the synthesis rate to the
                !      lamina-plus-sheath entity, so the sheath mass is moved from the stem
                !      basis to the lamina basis.
                lsb_m_sheath = 0.0_wp
                IF (lamina_synthesis_sheath_basis == 2) THEN
                  ! Mode 2: the printed measured masses of Bertheloot, Wu, Cournede & Andrieu
                  ! (2011) Ann Bot 108:1097-1109, Table 1c p.1100, sheath 272 mg and lamina
                  ! 420 mg per flowering culm.
                  lsb_m_sheath = lsb_a_sheath * MAX(0.0_wp, W_green)
                END IF
                ! The sheath cannot take more mass out of the stem than the stem has; the clamp
                ! keeps both bases non-negative and conserving.
                lsb_m_sheath = MIN(lsb_m_sheath, MAX(0.0_wp, W_stem))
                lsb_m_lam    = MAX(0.0_wp, W_green) + lsb_m_sheath
                lsb_m_stem   = MAX(0.0_wp, MAX(0.0_wp, W_stem) - lsb_m_sheath)
                s_lam   = n2_sigma_lam   * lsb_m_lam               * gnf_fN_day * gnf_fPAR_day * dtt_f
                s_stem  = n2_sigma_stem  * lsb_m_stem              * gnf_fN_day * gnf_fPAR_day * dtt_f
                s_chaff = n2_sigma_chaff * MAX(0.0_wp, W_chaff_gN) * gnf_fN_day * gnf_fPAR_day * dtt_f

                ! Driver record 167 ledgers: reported, never read back into the physics.
                IF (lamina_synthesis_sheath_basis >= 1) THEN
                  lsb_m_sheath_day = lsb_m_sheath
                  lsb_ratio_day    = 1.0_wp
                  IF (W_green > 0.0_wp) lsb_ratio_day = lsb_m_lam / W_green
                END IF
                IF (lamina_synthesis_sheath_basis >= 1 .AND. lsb_m_sheath > 0.0_wp) THEN
                  lsb_days_applied = lsb_days_applied + 1
                  lsb_msheath_cum  = lsb_msheath_cum + lsb_m_sheath
                  IF (W_green > 0.0_wp) lsb_ratio_last = lsb_m_sheath / W_green
                  lsb_syn_lam_gain_cum  = lsb_syn_lam_gain_cum  &
                    + n2_sigma_lam  * lsb_m_sheath * gnf_fN_day * gnf_fPAR_day * dtt_f
                  lsb_syn_stem_loss_cum = lsb_syn_stem_loss_cum &
                    + n2_sigma_stem * lsb_m_sheath * gnf_fN_day * gnf_fPAR_day * dtt_f
                END IF
                ! Eq. (10) of Bertheloot, Cournede & Andrieu (2011) Ann Bot 108:1085-1096: the
                ! phloem nitrogen influx into the roots follows the root DRY MASS influx and
                ! not thermal time. It is kept inside this guard so that the whole record is
                ! inert on a day with no thermal time.
                IF (root_limb) s_root = n2_p_root * conc_mob * MAX(0.0_wp, ont_dw_root_day)
              END IF
              gnf_syn_pot_day = s_lam
              s_lam_pot154    = s_lam

              ! The capacity ladder of driver records 87 and 97 stays and is still applied to
              ! the lamina synthesis, in the same statement as the incumbent path.
              IF ((itemF_Nclose_mode == 1 .OR. itemF_Nclose_mode == 3) .AND. &
                  (leafN_drawdown == 1 .OR. leafN_cap_route == 1)) THEN
                s_lam = MIN(s_lam, MAX(0.0_wp, sln_max_N * lai_basis - N_leaf))
                fnc_syn_clip_day = MAX(0.0_wp, s_lam_pot154 - s_lam)
                fnc_syn_clip_cum = fnc_syn_clip_cum + fnc_syn_clip_day
                IF (fnc_syn_clip_day > 0.0_wp) ont_days_cap_clip = ont_days_cap_clip + 1
                ont_cum_syn_clip = ont_cum_syn_clip + fnc_syn_clip_day
              END IF

              ! SUPPLY LIMIT, pro rata over every entity that asked. The grain influx is
              ! taken on the pool AS IT STANDS at that moment, which is the resolution of
              ! Bertheloot, Cournede & Andrieu (2011) Ann Bot 108:1085-1096, p.1090.
              gnp_claim_resid = 0.0_wp ; gnp_grain_extra = 0.0_wp
              IF (grain_n_partition_basis >= 3) gnp_claim_resid = MAX(0.0_wp, sink - grain_in)
              s_tot = s_lam + s_stem + s_chaff + s_root + gnp_claim_resid
              IF (s_tot > 0.0_wp) THEN
                s_scale = MIN(1.0_wp, MAX(0.0_wp, N_mobile_veg) / s_tot)
                IF (s_scale < 1.0_wp) ont_days_scaled = ont_days_scaled + 1
                s_lam = s_lam * s_scale ; s_stem = s_stem * s_scale
                s_chaff = s_chaff * s_scale ; s_root = s_root * s_scale
                gnp_grain_extra = gnp_claim_resid * s_scale
                N_leaf    = N_leaf    + s_lam
                N_stem    = N_stem    + s_stem
                N_chaff   = N_chaff   + s_chaff
                N_storage = N_storage + s_chaff
                IF (root_limb) npg_N_root = npg_N_root + s_root
                N_mobile_veg = N_mobile_veg &
                             - (s_lam + s_stem + s_chaff + s_root + gnp_grain_extra)
                IF (gnp_grain_extra > 0.0_wp) THEN
                  grain_in           = grain_in + gnp_grain_extra
                  gnp_m3_release_cum = gnp_m3_release_cum + gnp_grain_extra
                  gnp_m3_days        = gnp_m3_days + 1
                  IF (sink > 1.0e-12_wp) gnf_pool_scale = MIN(1.0_wp, grain_in / sink)
                  ont_bound_by = MERGE(1, 0, grain_in < sink - 1.0e-12_wp)
                END IF
              END IF
              act_syn         = s_lam
              gnf_syn_lam_day = s_lam
              gnf_cum_syn_lam = gnf_cum_syn_lam + s_lam

              ! ---- LEDGERS. The lamina ledger keeps one meaning throughout this routine:
              !      the lamina NET export, Eq. (5) minus the REALIZED Eq. (6). The stem and
              !      chaff ledgers carry their own.
              gnf_draw_lam   = d_lam - s_lam
              gnf_draw_stem  = d_stem
              gnf_draw_chaff = d_chaff
              grain_pool     = grain_in

              ! The ledger of today's assimilation already spent. Under a single pool that is
              ! all of it, because it entered the pool before any claim was served.
              take_assim_f  = uptake_in
              grain_assim_f = 0.0_wp          ! no separate assimilation claim exists on this path
              draw_mob      = grain_in        ! the pool is the grain's ONLY source here
              ! ---- ON-only diagnostics (mutate nothing further).
              ont_live_today    = .TRUE.
              ont_days          = ont_days + 1
              ont_dTT           = dtt_f
              ont_pool_before   = pool_entry
              ont_uptake_credit = uptake_in
              ont_D_lam = d_lam ; ont_D_stem = d_stem ; ont_D_chaff = d_chaff ; ont_D_root = d_root
              ont_sink_pot      = sink
              ont_grain_import  = grain_in
              ont_S_lam = s_lam ; ont_S_stem = s_stem ; ont_S_chaff = s_chaff ; ont_S_root = s_root
              ont_S_lam_pot     = s_lam_pot154
              ont_s_scale       = s_scale
              ont_f_N           = gnf_fN_day
              ont_f_PAR         = gnf_fPAR_day
              ont_mob_after     = N_mobile_veg
              ont_cum_D_lam = ont_cum_D_lam + d_lam ; ont_cum_D_stem = ont_cum_D_stem + d_stem
              ont_cum_D_chaff = ont_cum_D_chaff + d_chaff ; ont_cum_D_root = ont_cum_D_root + d_root
              ont_cum_S_lam = ont_cum_S_lam + s_lam ; ont_cum_S_stem = ont_cum_S_stem + s_stem
              ont_cum_S_chaff = ont_cum_S_chaff + s_chaff ; ont_cum_S_root = ont_cum_S_root + s_root
              ont_cum_grain  = ont_cum_grain  + grain_in
              ont_cum_uptake = ont_cum_uptake + uptake_in
              IF (ont_bound_by == 1) ont_days_pool_bound = ont_days_pool_bound + 1
            END BLOCK
            grain_import_day  = grain_pool
            N_storage         = N_storage + grain_import_day
            N_from_assim_day  = 0.0_wp               ! the pool is the only source on this path
            N_from_mobile_day = draw_mob
            N_remob_leaf_day  = MAX(0.0_wp, gnf_draw_lam)
            N_remob_stem_day  = gnf_draw_stem + gnf_draw_chaff
          ELSE
            ! ============================================================ !
            ! Eq. (4) = Eq. (6) - Eq. (5) for the LAMINA, both limbs present.
            ! ============================================================ !
            nph_lam = MAX(0.0_wp, N_leaf - n2_struct_lam * MAX(0.0_wp, W_green)) &
                    / (1.0_wp + gnf_r_mob_ph)                       ! Eqn 1 p.783, r = 0.10 p.786
            d_lam   = MIN(MAX(0.0_wp, N_leaf), n2_delta_TT * nph_lam * dtt_f)  ! Eqn 5 p.784
            N_leaf       = N_leaf       - d_lam
            N_mobile_veg = N_mobile_veg + d_lam            ! degradation feeds the common pool
            gnf_deg_lam_day = d_lam
            gnf_cum_deg_lam = gnf_cum_deg_lam + d_lam

            ! The two saturating factors of Eq. (6), evaluated on the state AFTER
            ! degradation: the same time step, as the conservation identity Eq. (13) p.785
            ! of Bertheloot, Martre & Andrieu (2008) Funct Plant Biol 35:781-796 requires.
            mg_tot   = MAX(0.0_wp, W_green) + MAX(0.0_wp, W_stem) + MAX(0.0_wp, W_chaff_gN)
            conc_mob = 0.0_wp
            IF (mg_tot > 1.0e-9_wp) conc_mob = MAX(0.0_wp, N_mobile_veg) / mg_tot   ! Eqn 7 p.784
            gnf_fN_day   = conc_mob / (conc_mob + n2_k1_Nmob)
            par_day      = n2_par_frac * MAX(0.0_wp, dTotrad_s) * 1.0e6_wp
            ! -------------------------------------------------------------------------- !
            ! Driver record 94: the light factor of Eq. (6) is evaluated at CANOPY-TOP PAR
            ! by the incumbent line, while Bertheloot, Martre & Andrieu (2008) Funct Plant
            ! Biol 35:781-796, Eq. (17) p.785, defines it on the per-lamina PAR after
            ! extinction, and their own flag lamina therefore sees less. The layered form
            ! uses the light profile of Bertheloot, Martre & Andrieu (2008) Plant Physiol
            ! 148:1707-1720, Eq. (1) p.1718, I(L) = I0 exp(-kL L).
            ! -------------------------------------------------------------------------- !
            canE_fPAR_top = par_day / (par_day + n2_k2_PAR)
            canE_fPAR_eff = canE_fPAR_top
            IF (lai > 1.0e-6_wp .AND. par_day > 0.0_wp) &
              canE_fPAR_eff = LOG( (par_day + n2_k2_PAR) &
                                 / (par_day * EXP(-kL_bertheloot * lai) + n2_k2_PAR) ) &
                            / (kL_bertheloot * lai)
            IF (canopy_n_flag_zone == 3 .OR. canopy_n_flag_zone == 4) THEN
              gnf_fPAR_day = canE_fPAR_eff
            ELSE
              gnf_fPAR_day = canE_fPAR_top
            END IF
            ! ---------------------------------------------------------------- !
            s_lam = gnf_sigma_lam * MAX(0.0_wp, W_green) * gnf_fN_day * gnf_fPAR_day * dtt_f
            gnf_syn_pot_day = s_lam
            ! -------------------------------------------------------------------------- !
            ! The leaf-nitrogen capacity and ordering (driver record 97 modes 1 and 3). The
            ! capacity is the maximum measured field-wheat Narea of Silva-Perez et al.
            ! (2020) J Exp Bot 71(7):2299-2311, doi:10.1093/jxb/erz439, Fig. 4B p.2305. It
            ! is applied on the SAME pass as the Eq. (6) synthesis flux, so that the
            ! capacity binds the inflow instead of being exceeded and corrected afterwards.
            ! -------------------------------------------------------------------------- !
            IF ((itemF_Nclose_mode == 1 .OR. itemF_Nclose_mode == 3) .AND. &
                (leafN_drawdown == 1 .OR. leafN_cap_route == 1)) THEN
              BLOCK
                REAL(wp) :: headroom, s_lam_pot
                s_lam_pot = s_lam
                headroom  = MAX(0.0_wp, sln_max_N * lai_basis - N_leaf)
                s_lam     = MIN(s_lam, headroom)
                fnc_syn_clip_day = MAX(0.0_wp, s_lam_pot - s_lam)
                fnc_syn_clip_cum = fnc_syn_clip_cum + fnc_syn_clip_day
              END BLOCK
            END IF
            ! -------------------------------------------------------------------------- !
            ! Driver record 132: retire the lamina re-synthesis after anthesis, so that
            ! refused synthesis nitrogen stays in the mobile pool instead of returning to
            ! the lamina.
            ! -------------------------------------------------------------------------- !
            lnh_syn_assim_off_day = 0.0_wp; lnh_syn_mob_off_day = 0.0_wp
            lnh_redirect_day = 0.0_wp; lnh_binding_day = .FALSE.
            ! -------------------------------------------------------------------------- !
            ! Driver record 133: state-driven lamina nitrogen re-synthesis. The grain
            ! Eq. (11) claim and the lamina Eq. (6) claim are served in order from the same
            ! pool, the grain first, and what the lamina cannot be served is refused rather
            ! than drawn from elsewhere.
            ! -------------------------------------------------------------------------- !
            pool_av   = MAX(0.0_wp, N_assim_avail) + MAX(0.0_wp, N_mobile_veg)
            claim_tot = sink + s_lam
            IF (claim_tot > pool_av) THEN
              need_f = claim_tot - pool_av
              gnf_draw_stem  = MIN(need_f, pot_stem)  ; need_f = need_f - gnf_draw_stem
              gnf_draw_chaff = MIN(need_f, pot_chaff)
              N_stem    = N_stem    - gnf_draw_stem
              N_chaff   = N_chaff   - gnf_draw_chaff
              N_storage = N_storage - gnf_draw_chaff   ! chaff N is tagged inside N_storage
              N_mobile_veg = N_mobile_veg + gnf_draw_stem + gnf_draw_chaff
              pool_av = pool_av + gnf_draw_stem + gnf_draw_chaff
            END IF
            IF (claim_tot > 1.0e-12_wp) THEN
              gnf_pool_scale = MIN(1.0_wp, pool_av / claim_tot)
            ELSE
              gnf_pool_scale = 1.0_wp
            END IF
            act_syn    = s_lam * gnf_pool_scale
            grain_pool = sink  * gnf_pool_scale
            gnf_syn_lam_day = act_syn
            gnf_cum_syn_lam = gnf_cum_syn_lam + act_syn
            gnf_supply_day  = pool_av
            ! -------------------------------------------------------------------------- !
            ! Debit the pool in the source order of Jamieson & Semenov (2000) Field Crops
            ! Res 68:21-29, p.24: today's assimilation first, then the standing mobile pool.
            ! The lamina's own draw is debited on the same footing as the grain.
            ! -------------------------------------------------------------------------- !
            take_assim_f  = MIN(MAX(0.0_wp, N_assim_avail), grain_pool)
            grain_assim_f = take_assim_f          ! the GRAIN's OWN share, kept separate: the
                                                  ! the ledger variable take_assim_f is raised
                                                  ! below by the lamina share, because the
                                                  ! surplus statement must exclude BOTH.
            draw_mob     = MIN(grain_pool - take_assim_f, MAX(0.0_wp, N_mobile_veg))
            N_mobile_veg = N_mobile_veg - draw_mob
            ! Whatever assimilation the grain did not take is available to the lamina
            ! before the standing pool is touched.
            BLOCK
              REAL(wp) :: syn_assim, syn_mob
              syn_assim = MIN(MAX(0.0_wp, N_assim_avail) - take_assim_f, act_syn)
              syn_mob   = MIN(act_syn - syn_assim, MAX(0.0_wp, N_mobile_veg))
              N_mobile_veg = N_mobile_veg - syn_mob
              N_leaf       = N_leaf       + syn_assim + syn_mob
              ! ------------------------------------------------------------------------ !
              ! The lamina ledger is the NET export, Eq. (5) minus Eq. (6), and is 
              ! therefore NEGATIVE on a day when synthesis exceeds degradation. Bertheloot, 
              ! Martre & Andrieu (2008) Funct Plant Biol 35:781-796, p.789, states that 
              ! during the phase of linear grain filling the lamina nitrogen declines, 
              ! which is the sign this ledger must carry.
              ! ------------------------------------------------------------------------ !
              gnf_draw_lam = d_lam - (syn_assim + syn_mob)
              ! The grain ledger records the assimilation the LAMINA took as spent, so the
              ! surplus statement below cannot hand it out twice.
              take_assim_f = take_assim_f + syn_assim
            END BLOCK
            grain_import_day  = grain_pool
            N_storage         = N_storage + grain_import_day
            N_from_assim_day  = grain_assim_f
            N_from_mobile_day = draw_mob
            N_remob_leaf_day  = MAX(0.0_wp, gnf_draw_lam)
            N_remob_stem_day  = gnf_draw_stem + gnf_draw_chaff
          END IF
          ! Driver record 154 mode 0 takes the ELSE arm above, which is the single-pool
          ! block.
          ELSE
            ! ---- driver record 86 mode 5 ----
            pot_lam   = n2_delta_TT * MAX(0.0_wp, N_leaf  - n2_struct_lam  * MAX(0.0_wp, W_green))    * dtt_f
            gnf_supply_day = MAX(0.0_wp, N_assim_avail) + MAX(0.0_wp, N_mobile_veg) &
                           + pot_lam + pot_stem + pot_chaff
            ! Eqn 11's min, drawn in the J&S source order.
            take_assim_f = MIN(MAX(0.0_wp, N_assim_avail), sink)
            short        = sink - take_assim_f
            draw_mob     = MIN(short, MAX(0.0_wp, N_mobile_veg))
            short        = short - draw_mob
            gnf_draw_lam   = MIN(short, pot_lam)   ; short = short - gnf_draw_lam
            gnf_draw_stem  = MIN(short, pot_stem)  ; short = short - gnf_draw_stem
            gnf_draw_chaff = MIN(short, pot_chaff)
            grain_assim_f  = take_assim_f          ! arm: the ladder has no competing claim
            ! --- move the N
            N_mobile_veg = N_mobile_veg - draw_mob
            N_leaf       = N_leaf  - gnf_draw_lam
            N_stem       = N_stem  - gnf_draw_stem
            N_chaff      = N_chaff - gnf_draw_chaff
            N_storage    = N_storage - gnf_draw_chaff     ! chaff N is tagged inside N_storage
            grain_import_day  = take_assim_f + draw_mob + gnf_draw_lam + gnf_draw_stem + gnf_draw_chaff
            N_storage         = N_storage + grain_import_day
            N_from_assim_day  = take_assim_f
            N_from_mobile_day = draw_mob
            N_remob_leaf_day  = gnf_draw_lam
            N_remob_stem_day  = gnf_draw_stem + gnf_draw_chaff
          END IF
          n3_cum_leaf = n3_cum_leaf + N_remob_leaf_day                   ! evaluation ledgers only
          n3_cum_stem = n3_cum_stem + gnf_draw_stem + gnf_draw_chaff
          gnf_cum_assim = gnf_cum_assim + grain_assim_f                 ! the GRAIN's share ONLY
          gnf_cum_mob   = gnf_cum_mob   + draw_mob
          gnf_cum_lam   = gnf_cum_lam   + gnf_draw_lam
          gnf_cum_stem  = gnf_cum_stem  + gnf_draw_stem
          gnf_cum_chaff = gnf_cum_chaff + gnf_draw_chaff
          ! ---------------------------------------------------------------------------- !
          ! Did the SINK bind, or the SOURCE? Martre, Porter, Jamieson & Triboi (2003)
          ! Plant Physiol 133:1959-1967, p.1963, measured grain nitrogen as "sink limited
          ! or colimited" for the first part of filling and source-limited thereafter, so
          ! both cases are distinguished and counted.
          ! ---------------------------------------------------------------------------- !
          dpa_f = gnf_days_fill

          ! Under the repaired law the grain claim competes with the lamina claim, so
          ! "the sink bound" means the shared pool met BOTH claims in full.
          IF (MERGE(gnf_pool_scale >= 1.0_wp, sink <= gnf_supply_day, gnf_eqn6)) THEN
            gnf_days_bind = gnf_days_bind + 1
            IF (gnf_first_bind_dpa < 0) gnf_first_bind_dpa = dpa_f
            gnf_last_bind_dpa = dpa_f
            IF (dpa_f > 15) gnf_bind_after_15d = gnf_bind_after_15d + 1
          END IF
          ! Assimilation the grain did not take stays vegetative, on the same bounded
          ! ladder as every other phase.
          N_surplus_veg_day = MAX(0.0_wp, N_assim_avail) - take_assim_f
          IF (leafN_drawdown == 1 .OR. leafN_cap_route == 1) THEN
            CALL route_surplus_bounded(N_surplus_veg_day, lai_basis, N_leaf, active_cap_N)
          ELSE
            N_leaf = N_leaf + N_surplus_veg_day
          END IF
        END IF

        ! ORGAN-STATE nitrogen closure: the five pools captured on entry against the
        ! same five now, minus the only external inflow of the day.
        resid = (N_leaf - nl0) + (N_stem - ns0) + (N_storage - nst0) &
              + (N_mobile_veg - nm0) - MAX(0.0_wp, N_assim_avail)

        ! Driver record 154 mode 1 is the only arm that moves ROOT nitrogen (the root
        ! limb of Bertheloot, Cournede & Andrieu 2011 Ann Bot 108:1085-1096, Eqs. (8)
        ! and (10)), so mode 1 alone adds the root pool to the closure.
        IF (organ_n_turnover_nema == 1) resid = resid + (npg_N_root - nr0_154)
        gnf_nres_day = resid
        ont_resid = resid
        IF (ABS(resid) > ont_resid_max) ont_resid_max = ABS(resid)
        IF (ABS(resid) > gnf_nres_max) gnf_nres_max = ABS(resid)
      END BLOCK
      RETURN
    END IF
    ! ------------------------------------------------------------------------------- !
    ! Demand-driven remobilization replaces a fixed first-order draw. Physiology:
    ! grain protein synthesis draws on the vegetative pools, and the measured 15N
    ! budgets of Kichey, Hirel, Heumez, Dubois & Le Gouis (2007) New Phytol
    ! 176:625-638 put 60-95 % of grain nitrogen at remobilized origin. The organ
    ! export ceilings are the measured remobilization efficiencies of Pask,
    ! Sylvester-Bradley, Jamieson & Foulkes (2012) Field Crops Res 126:104-118,
    ! doi:10.1016/j.fcr.2011.09.021, and the structural grain-nitrogen sink is that
    ! of Martre, Porter, Jamieson & Triboi (2003) Plant Physiol 133:1959-1967.
    ! ------------------------------------------------------------------------------- !
    stru_demand = n_stru_grain * MAX(0.0_wp, dW_storage_day) ! Martre structural N sink
    IF (dvs2 < dvs_grain_stru_end) THEN
      ! ----------------------------------------------------------------------------- !
      ! The sink-limited phase of Martre, Porter, Jamieson & Triboi (2003). The
      ! mobile pool is already-exported phloem nitrogen, available after today
      ! assimilation in the Jamieson & Semenov (2000) p.24 source order.
      ! ----------------------------------------------------------------------------- !
      source_tot = MAX(0.0_wp, N_assim_avail) + N_mobile_veg
      grain_import_day = MIN(source_tot, stru_demand)
      take_assim = MIN(grain_import_day, MAX(0.0_wp, N_assim_avail))
      N_from_assim_day = take_assim
      N_from_mobile_day = grain_import_day - take_assim
      N_mobile_veg = MAX(0.0_wp, N_mobile_veg - N_from_mobile_day)
      N_storage = N_storage + grain_import_day
      ! ----------------------------------------------------------------------------- !
      ! The chaff nitrogen share is set at the ear-construction window and not
      ! throughout filling: the chaff structure is built before anthesis (Slafer
      ! et al. 2015 Food Energy Secur 4(2):92-109, p.95; Fischer 1985 J Agric Sci
      ! Camb 105:447-461, p.447), so a share of the grain import taken after anthesis
      ! would be built from nitrogen the chaff never receives.
      ! ----------------------------------------------------------------------------- !
      N_surplus_veg_day = MAX(0.0_wp, N_assim_avail) - take_assim
      ! ----------------------------------------------------------------------------- !
      ! The vegetative surplus is routed through the bounded capacity ladder, so that
      ! nitrogen the grain did not take cannot raise the lamina above its holding
      ! capacity.
      ! ----------------------------------------------------------------------------- !
      IF (leafN_drawdown == 1 .OR. leafN_cap_route == 1) THEN   ! rec 87: capacity ladder independent of rec 62
        CALL route_surplus_bounded(N_surplus_veg_day, lai_basis, N_leaf, active_cap_N)
      ELSE
        N_leaf = N_leaf + N_surplus_veg_day
      END IF
      RETURN
    END IF
    ! ------------------------------------------------------------------------------- !
    ! Terminate the grain nitrogen import at the DEVELOPMENTAL end of grain filling
    ! (dvs2 >= 2, the phenology's own maturity), which is the measured termination of
    ! Martre, Porter, Jamieson & Triboi (2003) Plant Physiol 133:1959-1967.
    ! ------------------------------------------------------------------------------- !
    IF (dvs2 >= 2.0_wp .OR. canopy_senesced) THEN
      N_surplus_veg_day = MAX(0.0_wp, N_assim_avail)

      ! The same bounded ladder on the post-termination branch, so that the capacity
      ! holds on every path.
      IF (leafN_drawdown == 1 .OR. leafN_cap_route == 1) THEN   ! rec 87: capacity ladder independent of rec 62
        CALL route_surplus_bounded(N_surplus_veg_day, lai_basis, N_leaf, active_cap_N)
      ELSE
        N_leaf = N_leaf + N_surplus_veg_day
      END IF
      assimN_surplus_postfill = assimN_surplus_postfill + N_surplus_veg_day
      RETURN
    END IF
    ! --- source-limited phase: demand-driven linear draw ---
    IF (.NOT. n3_fill_snap) THEN
      n3_fill_snap = .TRUE.
      n3_dvs0  = dvs2
      n3_leaf0 = N_leaf
      n3_stem0 = N_stem
      n3_cum_leaf = 0.0_wp
      n3_cum_stem = 0.0_wp

      ! The lower bound of the lamina pool is SLN_min_N = 0.50 g N m-2: Bertheloot,
      ! Martre & Andrieu (2008) Plant Physiol 148:1707-1720,
      ! doi:10.1104/pp.108.124156, p.1711 and Fig. 2C p.1711, from which the derived
      ! span is 0.41-0.51 g N m-2.
      IF (leafN_drawdown == 1) THEN
        nJ_R0  = MAX(0.0_wp, N_leaf - SLN_min_N * lai_basis) &
               + MAX(0.0_wp, N_stem - n_stem_resid_conc * W_stem) &
               + N_mobile_veg
        nJ_cum = 0.0_wp
      END IF

      ! Under driver record 63 the organ structural floors are the measured
      ! concentrations of Bertheloot, Martre & Andrieu (2008) Plant Physiol
      ! 148:1707-1720, p.1711 (lamina 0.0079, sheath, internode, chaff and peduncle
      ! 0.0035 g N g-1 DM).
      IF (grainN_partition == 1) THEN
        nJ_R0  = MAX(0.0_wp, N_leaf  - n2_struct_lam  * MAX(0.0_wp, W_green)) &
               + MAX(0.0_wp, N_stem  - n2_struct_stem * MAX(0.0_wp, W_stem)) &
               + MAX(0.0_wp, N_chaff - MAX(n2_struct_chaff_use * MAX(0.0_wp, W_chaff_gN), cnr_floor_gN)) &
               + N_mobile_veg
        nJ_cum = 0.0_wp
      END IF
    END IF
    ! ------------------------------------------------------------------------------- !
    ! The grain draw under driver record 63. Once degradation and synthesis are
    ! both explicit, the net export of each organ is the difference of the two, and
    ! the grain is served from the culm-common pool rather than from the organs
    ! directly (Bertheloot, Martre & Andrieu 2008 Funct Plant Biol 35:781-796,
    ! p.785 and p.792).
    ! ------------------------------------------------------------------------------- !
    IF (grainN_partition == 1) THEN
      BLOCK
        REAL(wp) :: due, draw_mob

        ! The grain DAILY DEMAND is the source-determined thermal-time trajectory to a
        ! target nitrogen harvest index of Jamieson & Semenov (2000) Field Crops Res
        ! 68:21-29, p.24, hypothesis H1.
        due = MAX(0.0_wp, n4_nhi_target * n4_ag_anth * MIN(1.0_wp, dvs2 - 1.0_wp) &
                          - MAX(0.0_wp, N_storage - N_chaff))
        draw_mob = MIN(MAX(0.0_wp, due - MAX(0.0_wp, N_assim_avail)), MAX(0.0_wp, N_mobile_veg))
        N_mobile_veg = N_mobile_veg - draw_mob
        nJ_cum = nJ_cum + draw_mob
        N_from_mobile_day = draw_mob
        take_assim = MIN(MAX(0.0_wp, N_assim_avail), due)            ! J&S order: assim first
        N_from_assim_day = take_assim
        grain_import_day = take_assim + draw_mob
        N_storage = N_storage + grain_import_day
        ! evaluation ledgers: the organs' export today IS the B1 degradation
        N_remob_leaf_day = n2_deg_lam_day
        N_remob_stem_day = n2_deg_stem_day + n2_deg_chaff_day
        n3_cum_leaf = n3_cum_leaf + n2_deg_lam_day
        n3_cum_stem = n3_cum_stem + n2_deg_stem_day + n2_deg_chaff_day

        ! Assimilated nitrogen in excess of the grain demand stays vegetative, on the
        ! same bounded ladder the other two phases use.
        N_surplus_veg_day = MAX(0.0_wp, N_assim_avail) - take_assim
        IF (leafN_drawdown == 1 .OR. leafN_cap_route == 1) THEN   ! rec 87: capacity ladder independent of rec 62
          CALL route_surplus_bounded(N_surplus_veg_day, lai_basis, N_leaf, active_cap_N)
        ELSE
          N_leaf = N_leaf + N_surplus_veg_day
        END IF
      END BLOCK
      RETURN
    END IF
    IF (leafN_drawdown == 1) THEN
      ! ----------------------------------------------------------------------------- !
      ! The green-lamina nitrogen drawdown during grain filling. The rate is the
      ! first-order relative rate 0.0035 (degree-Cd)-1 fitted by Bertheloot, Martre &
      ! Andrieu (2008) Plant Physiol 148:1707-1720, doi:10.1104/pp.108.124156, p.1712
      ! and reported there as the value that minimizes the root mean square error;
      ! the same paper p.1715 states that the stem and sheath rate constant is lower
      ! but prints no value, so no stem rate is invented. Martre, Porter, Jamieson &
      ! Triboi (2003) Plant Physiol 133:1959-1967, doi:10.1104/pp.103.030585, p.1968,
      ! names green-tissue nitrogen as part of the remobilized source.
      ! ----------------------------------------------------------------------------- !
      BLOCK
        REAL(wp) :: prog, due, f_k, pot_leaf, pot_stem, draw_leaf, draw_stem, draw_mob
        prog = MIN(1.0_wp, MAX(0.0_wp, (dvs2 - n3_dvs0) / MAX(1.0e-6_wp, 2.0_wp - n3_dvs0)))
        due  = MAX(0.0_wp, nJ_R0 * prog - nJ_cum)
        f_k  = MIN(1.0_wp, k_leafN_TT * MAX(0.0_wp, dTT_day_N))
        pot_leaf = f_k * MAX(0.0_wp, N_leaf - SLN_min_N * lai_basis)
        pot_stem = f_k * MAX(0.0_wp, N_stem - n_stem_resid_conc * W_stem)
        draw_leaf = MIN(due, pot_leaf)
        draw_stem = MIN(due - draw_leaf, pot_stem)
        draw_mob  = MIN(due - draw_leaf - draw_stem, N_mobile_veg)
        N_remob_leaf_day = draw_leaf
        N_remob_stem_day = draw_stem
        N_from_mobile_day = draw_mob
        n3_cum_leaf = n3_cum_leaf + draw_leaf         ! evaluation ledger only
        n3_cum_stem = n3_cum_stem + draw_stem         ! evaluation ledger only
        nJ_cum = nJ_cum + draw_leaf + draw_stem + draw_mob
        N_mobile_veg = N_mobile_veg - draw_mob
      END BLOCK
      take_assim = MAX(0.0_wp, N_assim_avail)                        ! J&S order: assim first (unchanged)
      N_from_assim_day = take_assim
      grain_import_day = take_assim + N_from_mobile_day + N_remob_leaf_day + N_remob_stem_day
      N_leaf    = N_leaf    - N_remob_leaf_day
      N_stem    = N_stem    - N_remob_stem_day
      N_storage = N_storage + grain_import_day
      ! The surplus is zero at run time in this phase (the grain takes the full
      ! assimilation); it is kept on the bounded ladder so that every path is one.
      N_surplus_veg_day = MAX(0.0_wp, N_assim_avail) - take_assim
      CALL route_surplus_bounded(N_surplus_veg_day, lai_basis, N_leaf, active_cap_N)
      RETURN
    END IF
    BLOCK
      REAL(wp) :: prog, due_leaf, due_stem, cap_leaf, cap_stem
      prog = MIN(1.0_wp, MAX(0.0_wp, (dvs2 - n3_dvs0) / MAX(1.0e-6_wp, 2.0_wp - n3_dvs0)))
      ! Today's cumulative-due per organ (linear in dvs2), minus already drawn:
      due_leaf = MAX(0.0_wp, n3_eff_leaf * n3_leaf0 * prog - n3_cum_leaf)
      due_stem = MAX(0.0_wp, n3_eff_stem * n3_stem0 * prog - n3_cum_stem)
      ! Floors (SLN_min, Bertheloot 2008 PP 148 p.1711; Yin stem structural) bound today's draw:
      cap_leaf = MAX(0.0_wp, N_leaf - SLN_min_N * MAX(0.0_wp, lai))
      cap_stem = MAX(0.0_wp, N_stem - n_struct_stem_N * W_stem)
      N_remob_leaf_day = MIN(due_leaf, cap_leaf)
      N_remob_stem_day = MIN(due_stem, cap_stem)
      n3_cum_leaf = n3_cum_leaf + N_remob_leaf_day
      n3_cum_stem = n3_cum_stem + N_remob_stem_day
    END BLOCK
    take_assim = MAX(0.0_wp, N_assim_avail)                          ! J&S order: assim first
    N_from_assim_day = take_assim
    ! ------------------------------------------------------------------------------- !
    ! The mobile pool, leaf-origin and already exported, is drawn in full after
    ! assimilation, and counted toward the lamina cumulative so that the
    ! remobilization-efficiency cap of Pask et al. (2012) Field Crops Res
    ! 126:104-118 still bounds the TOTAL lamina export.
    ! ------------------------------------------------------------------------------- !
    N_from_mobile_day = N_mobile_veg
    N_mobile_veg = 0.0_wp
    n3_cum_leaf = n3_cum_leaf + N_from_mobile_day
    grain_import_day = take_assim + N_from_mobile_day + N_remob_leaf_day + N_remob_stem_day
    N_leaf    = N_leaf    - N_remob_leaf_day
    N_stem    = N_stem    - N_remob_stem_day
    N_storage = N_storage + grain_import_day

    ! Assimilation not imported by the grain (sink-limited phase) sustains leaf N.
    N_surplus_veg_day = MAX(0.0_wp, N_assim_avail) - take_assim
    ! Driver record 87: the fourth and last of the unbounded surplus statements.
    IF (leafN_cap_route == 1) THEN
      CALL route_surplus_bounded(N_surplus_veg_day, lai_basis, N_leaf, active_cap_N)
    ELSE
      N_leaf = N_leaf + N_surplus_veg_day
    END IF
  END SUBROUTINE allocate_grain_N_source_sink_daily
  ! --------------------------------------------------------------------------------- !
  ! Capacity-bounded routing of the vegetative surplus. The ceiling is the
  ! maximum measured field-wheat Narea of Silva-Perez, Furbank, Condon & Evans
  ! (2020) J Exp Bot 71(7):2299-2311, doi:10.1093/jxb/erz439, Fig. 4B p.2305,
  ! measured genotype-mean maximum 3.85 g N m-2.
  ! --------------------------------------------------------------------------------- !
  SUBROUTINE route_surplus_bounded(surplus, lai, N_leaf, cap_N)

    ! The optional ceiling argument lets a caller substitute the re-sourced ceiling
    ! of driver record 132 mode 2 for the default one.
    REAL(wp), INTENT(IN)    :: surplus, lai
    REAL(wp), INTENT(INOUT) :: N_leaf
    REAL(wp), INTENT(IN), OPTIONAL :: cap_N
    REAL(wp) :: to_leaf, cap_use
    cap_use = sln_max_N
    IF (PRESENT(cap_N)) cap_use = cap_N
    to_leaf = MIN(MAX(0.0_wp, surplus), &
                  MAX(0.0_wp, cap_use * MAX(0.0_wp, lai) - N_leaf))
    N_leaf = N_leaf + to_leaf
    N_mobile_veg = N_mobile_veg + MAX(0.0_wp, surplus) - to_leaf
    lnc_to_mobile_day = lnc_to_mobile_day + MAX(0.0_wp, surplus) - to_leaf   ! rec 87 diagnostic only
  END SUBROUTINE route_surplus_bounded
  ! --------------------------------------------------------------------------------- !
  ! Development-stage-dependent specific leaf area: it declines through
  ! development as the leaves thicken. Measured wheat ranges: Bos, Vos & Struik
  ! (2000) Ann Bot 85:243-249, Table 1 (winter wheat 20-28 m2 kg-1) and Araus,
  ! Reynolds & Acevedo (1993) Crop Sci 33:1273-1279.
  ! --------------------------------------------------------------------------------- !
  FUNCTION sla_from_dvs(dvs2) RESULT(sla)
    REAL(wp), INTENT(IN) :: dvs2
    REAL(wp) :: sla, d
    d = MIN(2.0_wp, MAX(0.0_wp, dvs2))
    sla = sla_dvs0 + (sla_dvs2 - sla_dvs0) * d / 2.0_wp   ! 0.026 -> 0.018, = 0.022 at dvs2=1
    sla = sla * (gParam_sla / MAX(1.0e-9_wp, sla_anchor)) ! preserve input calibratability
  END FUNCTION sla_from_dvs
  ! --------------------------------------------------------------------------------- !
  ! Yield-component structure: the day reproductive dry-matter increment (the
  ! storage-pool gain from allocate_daily_biomass) is split between grain and
  ! chaff.
  ! --------------------------------------------------------------------------------- !
  SUBROUTINE allocate_yield_components_daily(dvs2, dW_repro, W_stem, W_grain, W_chaff, &
                                             kernel_number, grain_number_set, dW_grain_out, &
                                             t_air_day_in, iday_in, W_storage)
    REAL(wp), INTENT(IN)    :: dvs2, dW_repro
    REAL(wp), INTENT(INOUT) :: W_stem, W_grain, W_chaff, kernel_number
    ! ------------------------------------------------------------------------------- !
    ! The reproductive pool W_storage, which allocate_daily_biomass credited with
    ! the whole storage share for this day, is debited here by whatever leaves it
    ! for the stem, so that W_storage = W_grain + W_chaff holds on every day.
    ! ------------------------------------------------------------------------------- !
    REAL(wp), INTENT(INOUT) :: W_storage
    LOGICAL,  INTENT(INOUT) :: grain_number_set
    REAL(wp), INTENT(OUT)   :: dW_grain_out

    ! The day screen-air daily mean (degC) and the model day, read only when the
    ! grain-filling law of driver record 91 is active.
    REAL(wp), INTENT(IN)    :: t_air_day_in
    INTEGER,  INTENT(IN)    :: iday_in
    REAL(wp) :: to_chaff, to_grain, grain_cap, room, draw, gain
    REAL(wp) :: gfb_x, gfb_cm, gfb_rate, gfb_before
    REAL(wp) :: gnm_Nstru_at_Dcd   ! rec 130 = 3, local: Nstru(Tt=Dcd), Eq. 3's exponential-phase endpoint
    ! Driver record 92. Organ states ON ENTRY, for the closure below; `wsc_c_in` is
    ! the assimilate this routine actually consumed, accumulated where it is taken.
    REAL(wp) :: wsc_wg0, wsc_ws0, wsc_wc0, wsc_c_in, wsc_label, wsc_room
    ! Driver record 117: the day reproductive stream AFTER the chaff sink has been
    ! paid. On mode 0 it is dW_repro and is never decremented.
    REAL(wp) :: dW_repro_use
    ! B1: the day's two cross-boundary fluxes, g DM m-2 d-1 (see the W_storage note above).
    REAL(wp) :: b1_repro_to_stem, b1_reserve_to_grain

    b1_repro_to_stem = 0.0_wp; b1_reserve_to_grain = 0.0_wp
    dW_grain_out = 0.0_wp
    cres_to_grain_assim = 0.0_wp; cres_to_grain_res = 0.0_wp
    cres_sink_demand = 0.0_wp; cres_sink_unfilled = 0.0_wp
    cres_respired = 0.0_wp; cres_ledger_resid = 0.0_wp
    wsc_wg0 = W_grain; wsc_ws0 = W_stem; wsc_wc0 = W_chaff
    wsc_c_in = 0.0_wp
    IF (repro_surplus_fate == 1) rsf_resid_day = 0.0_wp   ! rec 127: today's refused surplus
    IF (sink_full_carbon_fate == 1) THEN                   ! rec 129: today's overflow-loss / stored / surplus
      cfo_loss_day = 0.0_wp ; cfo_wsc_day = 0.0_wp ; cfo_surplus_day = 0.0_wp
    END IF
    ! ------------------------------------------------------------------------------- !
    ! Driver record 137: the window-cumulative standing nitrogen ratio that carries
    ! the nitrogen term of floret survival. The measured nitrogen effect on floret
    ! survival is Ferrante, Savin & Slafer (2020) Field Crops Res 256:107908,
    ! doi:10.1016/j.fcr.2020.107908.
    ! ------------------------------------------------------------------------------- !
    fnl_step_today = .FALSE.
    dW_repro_use = dW_repro
    cgs_post_today = .FALSE.
    ! ------------------------------------------------------------------------------- !
    ! Mode 0 keeps the original control flow (an early return on a zero increment).
    ! Mode 1 must continue, because the reserve route can move carbon on a day with
    ! no reproductive increment (Bidinger, Musgrave & Fischer 1977, Nature
    ! 270:431-433).
    ! ------------------------------------------------------------------------------- !
    IF (dW_repro <= 0.0_wp .AND. carbon_reserve_mode == 0) RETURN

    IF (dvs2 < ear_end_dvs2) THEN
      IF (dW_repro > 0.0_wp) THEN

        ! Ear-construction window: the chaff structure is built and a little grain
        ! begins. Before anthesis ALL reproductive dry matter is spike structure, and
        ! grain mass cannot begin before anthesis (Slafer et al. 2015 Food Energy Secur
        ! 4(2):92-109, p.95).
        IF (dvs2 < 1.0_wp) THEN
          to_chaff = dW_repro
        ELSE
          to_chaff = chaff_share * dW_repro
        END IF
        to_grain = dW_repro - to_chaff
        mnp_dw_chaff_day = to_chaff           ! P_chaff
        mnp_dw_grain_day = to_grain           ! P_grain
        W_chaff  = W_chaff + to_chaff
        W_grain  = W_grain + to_grain
        dW_grain_out = to_grain
        cres_to_grain_assim = to_grain
        wsc_c_in = wsc_c_in + dW_repro        ! assimilate this routine took
      END IF
    ELSE
      ! Grain-number set once at ear_end (fruiting efficiency x ear DW).
      IF (.NOT. grain_number_set) THEN

        ! Driver record 89 modes 1 and 3: grain number is PREDICTED from the spike
        ! carbon each axis actually accrued, through the two measured components of
        ! Fischer (1985) J Agric Sci Camb 105:447-461.
        IF (grain_number_mech == 1 .OR. grain_number_mech == 3) THEN
          CALL gnx_set_grain_number(MAX(0.0_wp, W_chaff), kernel_number)
        ELSE
          kernel_number = fruiting_efficiency * MAX(0.0_wp, W_chaff)
        END IF
        grain_number_set = .TRUE.

        ! Driver record 92: the anthesis anchor for the deposition close and for every
        ! days-after-anthesis diagnostic. It is the same node the grain number is set
        ! on.
        IF (stem_wsc_mech >= 1) THEN
          wsc_day_anth  = iday_in
          wsc_pool_anth = W_reserve
          wsc_stem_anth = W_stem
        END IF
        ! The thermal clock starts on the day the grain number is set, i.e. at anthesis
        ! (every shipped deck sets ear_end_dvs2 = 1.0).
        IF (grain_fill_mech >= 1) THEN
          gfb_tt = 0.0_wp ; gfb_day0 = iday_in ; gfb_days = 0 ; gfb_tsum = 0.0_wp
        END IF
      END IF
      ! ----------------------------------------------------------------------------- !
      ! Driver record 117: the post-anthesis spike-non-grain (chaff) growth sink. The
      ! published chaff share bracket is assembled in Slafer et al. (2015) Food
      ! Energy Secur 4:92-109, folios 96-97.
      ! ----------------------------------------------------------------------------- !
      cgs_demand_day = 0.0_wp ; cgs_taken_day = 0.0_wp ; cgs_unmet_day = 0.0_wp
      IF (grain_fill_mech >= 1 .AND. grain_fill_mech /= 3) THEN

        ! ------------------------- The rate x duration law ------------------------- !
        ! Yin, Goudriaan, Lantinga, Vos & Spiertz (2003) "A flexible sigmoid function
        ! of determinate growth", Annals of Botany 91(3):361-371,
        ! doi:10.1093/aob/mcg029, Eq. (9) with the rate form of Eq. (A1) and the
        ! measured durations of their Table 2; the base temperature and the filling
        ! duration are derived from that table by the reciprocal-duration method of
        ! Fischer (1985) J Agric Sci Camb 105:447-461, Table 6.
        ! --------------------------------------------------------------------------- !
        IF (gfb_end_bound == 0) gfb_tt = gfb_tt + MAX(0.0_wp, t_air_day_in - gfb_tbase)
        gfb_x  = gfb_tt / gfb_dfill
        gfb_cm = MAX(0.0_wp, gfb_cm_ref + gfb_cm_slope * (t_air_day_in - gfb_t_ref))
        IF (gfb_cm <= 0.0_wp) gfb_cm_clamped = .TRUE.
        IF (gfb_x > 1.0_wp) THEN
          gfb_rate = 0.0_wp
          IF (gfb_end_bound == 0) THEN
            gfb_day_end = iday_in ; gfb_end_bound = 1   ! also freezes the clock
          END IF
          IF (.NOT. grain_fill_ended .AND. grain_number_set) grain_fill_ended = .TRUE.
        ELSE IF (gfb_x <= 0.0_wp) THEN
          ! tau has not started to accrue: a filling day at or below the base
          ! temperature. The demand is zero TODAY and the window is NOT closed, since
          ! closing it here would end grain filling on the first cold day after anthesis.
          gfb_rate = 0.0_wp
        ELSE
          ! g(x) = ((1-x)/(1-u)) (x/u)**(u/(1-u)), Yin et al. (2003) Eq. (A1) with
          ! delta = 1, normalized so that g = 1 at x = u and g = 0 at x = 0 and x = 1.
          gfb_rate = gfb_cm * ((1.0_wp - gfb_x) / (1.0_wp - gfb_u)) &
                            * (gfb_x / gfb_u) ** (gfb_u / (1.0_wp - gfb_u))
        END IF
        ! Developmental maturity still ends filling if it arrives first; which of the
        ! two bounds fired is reported per deck and never assumed.
        IF (dvs2 >= 2.0_wp) THEN
          gfb_rate = 0.0_wp
          IF (gfb_end_bound == 0) THEN
            gfb_day_end = iday_in ; gfb_end_bound = 2
          END IF
          IF (.NOT. grain_fill_ended .AND. grain_number_set) grain_fill_ended = .TRUE.
        END IF
        room = kernel_number * gfb_rate / 1000.0_wp        ! g m-2 d-1 -- a FLUX, not a stock
        ! --------------------------------------------------------------------------- !
        ! Driver record 118 mode 1: the cultivar FINAL kernel mass ceiling. The
        ! measured cultivar maxima are Yin et al. (2003) Ann Bot 91:361, Table 2, and
        ! Sofield, Evans, Cook & Wardlaw (1977) Aust J Plant Physiol 4:785,
        ! doi:10.1071/PP9770785.
        ! --------------------------------------------------------------------------- !
        IF (kernel_mass_ceiling == 1) THEN
          grain_cap = kernel_number * TKW_potential / 1000.0_wp
          IF (gfb_end_bound == 3) THEN
            gfb_rate = 0.0_wp ; room = 0.0_wp
          ELSE IF (grain_number_set .AND. gfb_x > 0.0_wp) THEN
            room = MIN(room, MAX(0.0_wp, grain_cap - W_grain))
            IF (grain_cap - W_grain <= 1.0e-9_wp .AND. gfb_end_bound == 0) THEN
              gfb_rate = 0.0_wp ; room = 0.0_wp
              gfb_day_end = iday_in ; gfb_end_bound = 3   ! freezes the clock too
              IF (.NOT. grain_fill_ended) grain_fill_ended = .TRUE.
            END IF
          END IF
        END IF
        ! --------------------------------------------------------------------------- !
        ! Driver record 125 modes 2 and 3: the kernel-mass POTENTIAL, set at the end of
        ! the endosperm cell-division window. Inactive alternative, see Part 1 Table 6.
        ! --------------------------------------------------------------------------- !
        IF (kernel_potential_law >= 2 .AND. grain_number_set) THEN
          IF (.NOT. kpt_set) THEN
            kpt_tt      = kpt_tt + MAX(0.0_wp, t_air_day_in)
            kpt_sup_cum = kpt_sup_cum + MAX(0.0_wp, dW_repro_use)
            kpt_tsum    = kpt_tsum + t_air_day_in
            kpt_days    = kpt_days + 1
            IF (kpt_tt >= gnf_tau_degCd .AND. kernel_number > 0.0_wp) THEN
              kpt_tmean      = kpt_tsum / REAL(kpt_days, wp)
              kpt_sk_mean    = 1000.0_wp * kpt_sup_cum / REAL(kpt_days, wp) / kernel_number
              kpt_dfill_days = gfb_dfill / MAX(0.1_wp, kpt_tmean - gfb_tbase)
              kpt_w_yin      = gfb_shape_int * MAX(0.0_wp, gfb_cm_ref + gfb_cm_slope * (kpt_tmean - gfb_t_ref)) &
                             * kpt_dfill_days
              kpt_w_raw      = kpt_sk_mean * kpt_dfill_days
              kpt_w_pot      = MIN(kpt_env_max, MAX(kpt_env_min, kpt_w_raw))
              kpt_bound      = 0
              IF (kpt_w_raw < kpt_env_min) kpt_bound = 1
              IF (kpt_w_raw > kpt_env_max) kpt_bound = 2
              kpt_set = .TRUE. ; kpt_day_set = iday_in
            END IF
          END IF
          IF (kpt_set) THEN
            grain_cap = kernel_number * kpt_w_pot / 1000.0_wp
            IF (gfb_end_bound == 4) THEN
              gfb_rate = 0.0_wp ; room = 0.0_wp
            ELSE IF (gfb_x > 0.0_wp) THEN
              room = MIN(room, MAX(0.0_wp, grain_cap - W_grain))
              IF (grain_cap - W_grain <= 1.0e-9_wp .AND. gfb_end_bound == 0) THEN
                gfb_rate = 0.0_wp ; room = 0.0_wp
                gfb_day_end = iday_in ; gfb_end_bound = 4   ! freezes the clock too
                IF (.NOT. grain_fill_ended) grain_fill_ended = .TRUE.
              END IF
            END IF
          END IF
        END IF
        IF (gfb_rate > 0.0_wp) THEN
          gfb_days = gfb_days + 1 ; gfb_tsum = gfb_tsum + t_air_day_in
          gfb_demand_tot = gfb_demand_tot + room
          IF (gfb_rate > gfb_rate_max) gfb_rate_max = gfb_rate
          IF (gfb_cm < gfb_cm_min) gfb_cm_min = gfb_cm
          IF (gfb_cm > gfb_cm_max) gfb_cm_max = gfb_cm
        END IF
        ! --------------------------------------------------------------------------- !
        ! Driver record 125 modes 2 and above: the potential is set once, at the end of
        ! the cell-division window, and is not re-evaluated afterwards, so that a later
        ! change of grain number cannot move a potential that has already been set.
        ! Inactive alternative, see Part 1 Table 6.
        ! --------------------------------------------------------------------------- !
      ELSE
        ! Sink-limited grain fill: cannot exceed kernel_number x potential kernel weight.
        grain_cap = kernel_number * TKW_potential / 1000.0_wp
        room      = MAX(0.0_wp, grain_cap - W_grain)
        ! The unconstrained end of grain filling (Martre et al. 2003) is reached when
        ! the carbon sink is full. Latched once, with the day recorded for the log.
        IF (.NOT. grain_fill_ended .AND. room <= 1.0e-9_wp .AND. grain_number_set) THEN
          grain_fill_ended = .TRUE.
        END IF
      END IF
      ! ----------------------------------------------------------------------------- !
      ! Driver record 130: the grain nitrogen minimum concentration. The measured
      ! floor is 1.5 % nitrogen in the grain (Sinclair & Amir 1992, as relayed by
      ! Jamieson & Semenov 2000 Field Crops Res 68:21-29, pp.23-24).
      ! ----------------------------------------------------------------------------- !
      gnm_room_dm_118   = room
      gnm_room_dm_125   = room
      gnm_room_dm_floor = HUGE(1.0_wp)
      gnm_floor_bind_day = .FALSE.
      IF (grain_n_min_conc == 1) THEN
        ! --------------------------------------------------------------------------- !
        ! Mode 1, cumulative: the room is MAX(0, N_grain/min_conc - W_grain), with the
        ! grain nitrogen taken as the spike pool net of the chaff share. Inactive
        ! alternative, see Part 1 Table 6.
        ! --------------------------------------------------------------------------- !
        gnm_room_dm_floor = MAX(0.0_wp, gnm_N_grain_in / gnm_min_conc - W_grain)
        room = MIN(room, gnm_room_dm_floor)
      ELSE IF (grain_n_min_conc == 2) THEN
        ! --------------------------------------------------------------------------- !
        ! Mode 2, the guarded daily increment: dW_grain <= dN_grain/min_conc on
        ! yesterday's grain nitrogen import, an explicit one-day lag. Inactive
        ! alternative, see Part 1 Table 6.
        ! --------------------------------------------------------------------------- !
        gnm_guard_dm = 0.10_wp * kernel_number * TKW_potential / 1000.0_wp
        IF (W_grain > gnm_guard_dm) THEN
          gnm_room_dm_floor = MAX(0.0_wp, gnm_dN_grain_prev / gnm_min_conc)
          room = MIN(room, gnm_room_dm_floor)
        END IF
      ELSE IF (grain_n_min_conc == 3) THEN
        ! --------------------------------------------------------------------------- !
        ! Mode 3, the soft structural floor, the released mode. Martre, Porter,
        ! Jamieson & Triboi (2003) Plant Physiol 133:1959-1967,
        ! doi:10.1104/pp.103.030585, Eq. (3) p.1964, replaces the fixed ceiling by a
        ! structural nitrogen mass per grain that accumulates over thermal time.
        ! --------------------------------------------------------------------------- !
        gnm_Tt_degCd     = gnf_tt
        gnm_Nstru_at_Dcd = gnm_Nstru0 * EXP(gnm_kcd * gnf_tau_degCd)
        IF (gnm_Tt_degCd < gnf_tau_degCd) THEN
          ! exponential / cell-division phase: Nstru(Tt) = Nstru0*EXP(kcd*Tt)
          gnm_Nstru_pergrain = gnm_Nstru0 * EXP(gnm_kcd * gnm_Tt_degCd)
          gnm_phase_day      = 0
        ELSE IF (gnm_Tt_degCd <= gnm_Der) THEN
          ! ------------------------------------------------------------------------- !
          ! The linear, endoreduplication phase. The published linear-phase rate is the
          ! structural nitrogen at the end of cell division divided by that duration,
          ! times the mean temperature (Martre et al. 2003 p.1964); integrated over
          ! thermal time this is the expression below.
          ! ------------------------------------------------------------------------- !
          gnm_Nstru_pergrain = gnm_Nstru_at_Dcd * (gnm_Tt_degCd / gnf_tau_degCd)
          gnm_phase_day      = 1
        ELSE
          ! frozen past Der: Eq. 3's third branch, rate = 0.
          gnm_Nstru_pergrain = gnm_Nstru_at_Dcd * (gnm_Der / gnf_tau_degCd)
          gnm_phase_day      = 2
        END IF
        ! --------------------------------------------------------------------------- !
        ! Translate the absolute per-grain nitrogen mass into the SAME %DM-ceiling
        ! algebra that modes 1 and 2 use, through the grain capacity of driver records
        ! 118 and 125. That translation is this model's own and is not in Martre et al.
        ! (2003).
        ! --------------------------------------------------------------------------- !
        gnm_Nstru_total = kernel_number * gnm_Nstru_pergrain
        gnm_c_stru_day  = gnm_Nstru_total / MAX(1.0e-9_wp, grain_cap)
        ! Guard: divide only when the structural grain-nitrogen mass of the day is
        ! positive. The ratio starts near zero just after anthesis and only rises.
        IF (gnm_c_stru_day > 0.0_wp) THEN
          gnm_room_dm_floor = MAX(0.0_wp, gnm_N_grain_in / gnm_c_stru_day - W_grain)
          room = MIN(room, gnm_room_dm_floor)
        END IF
        ! Nitrogen-side diagnostic: no flux, no change to `room` beyond the statement
        ! above and no change to the import law of driver record 86.
        gnm_dNstru_day = MAX(0.0_wp, gnm_Nstru_total - gnm_Nstru_total_prev)
        gnm_Nstru_paid = gnm_Nstru_paid + MIN(gnm_dNstru_day, MAX(0.0_wp, gnm_dN_grain_prev))
        gnm_Nsto       = gnm_N_grain_in - gnm_Nstru_paid
        gnm_Nstru_total_prev = gnm_Nstru_total
      END IF
      IF (grain_n_min_conc >= 1) THEN
        gnm_to_grain_before = MIN(dW_repro_use, gnm_room_dm_118)
        gnm_to_grain_after  = MIN(dW_repro_use, room)
        gnm_refused_dm_day  = MAX(0.0_wp, gnm_to_grain_before - gnm_to_grain_after)
        gnm_refused_dm_cum  = gnm_refused_dm_cum + gnm_refused_dm_day
        gnm_floor_bind_day  = gnm_refused_dm_day > 1.0e-12_wp
        IF (gnm_floor_bind_day) gnm_days_bind = gnm_days_bind + 1
      END IF
      gfb_before = W_grain
      cres_sink_demand = room
      IF (dW_repro_use > 0.0_wp) THEN
        to_grain  = MIN(dW_repro_use, room)
        mnp_dw_grain_day = to_grain   ! P_grain, ASSIMILATE only (the reserve draw below is remob's, not this basis's)
        ! --------------------------------------------------------------------------- !
        ! Surplus when the grain sink is full goes to the STEM RESERVE pool, not to the
        ! chaff. Schnyder (1993) New Phytol 123:233-245, Summary p.233, measures the
        ! fructan pool as the destination of assimilate the grain cannot take.
        ! --------------------------------------------------------------------------- !
        BLOCK
          REAL(wp) :: surplus
          surplus = dW_repro_use - to_grain
          W_stem    = W_stem + surplus
          b1_repro_to_stem = b1_repro_to_stem + surplus   ! B1: leaves the reproductive pool
          ! ------------------------------------------------------------------------- !
          ! Driver record 92: the surplus still goes to the stem pool, so NO CARBON
          ! MOVES; what changes is how much of it may be LABELED as reserve. The window
          ! is that of Schnyder (1993) New Phytol 123:233-245, p.238, where the store
          ! "accumulate[s] in the stem until 10-20 d after anthesis", and Summary p.233.
          ! ------------------------------------------------------------------------- !
          wsc_label = surplus
          IF (stem_wsc_mech == 1 .OR. stem_wsc_mech == 2 .OR. stem_wsc_mech == 4) THEN
            ! ----------------------------------------------------------------------- !
            ! Driver record 129: while the grain sink is full - which is exactly what this
            ! block computes - the record-92 deposition WINDOW is lifted, because the
            ! measured accumulation continues while the sink is unable to take the
            ! assimilate (Schnyder 1993 New Phytol 123:233-245, p.238 and Summary p.233).
            ! ----------------------------------------------------------------------- !
            IF (sink_full_carbon_fate /= 1 .OR. wsc_fate_law == 1) THEN
              IF (wsc_day_anth > 0 .AND. (iday_in - wsc_day_anth) > wsc_accum_end_dpa) THEN
                IF (wsc_label > 0.0_wp) THEN
                  wsc_refused_win = wsc_refused_win + wsc_label
                  wsc_days_winbnd = wsc_days_winbnd + 1
                  ! rec 148 ON-only diagnostic of the term-A refusal
                  IF (wsc_fate_law == 1) THEN
                    wfl_a_refused_day = wfl_a_refused_day + wsc_label
                    wfl_a_refused_cum = wfl_a_refused_cum + wsc_label
                    wfl_a_days        = wfl_a_days + 1
                  END IF
                END IF
                wsc_label = 0.0_wp
              END IF
            END IF
            IF (repro_surplus_fate == 1) THEN
              ! rec 127: capacity read on the stem WITHOUT today's surplus (only the labeled part stays)
              wsc_room = MAX(0.0_wp, wsc_cap_frac * (W_stem - surplus) - W_reserve)
            ELSE
            wsc_room = MAX(0.0_wp, wsc_cap_frac * W_stem - W_reserve)
            END IF
            IF (wsc_label > wsc_room) THEN
              wsc_refused_cap = wsc_refused_cap + (wsc_label - wsc_room)
              wsc_days_capbnd = wsc_days_capbnd + 1
              wsc_label = wsc_room
            END IF
          END IF
          ! ------------------------------------------------------------------------- !
          ! Driver record 127: the UNLABELED remainder, refused by the record-92 window
          ! or capacity, does not become stem structure; it is removed from the stem pool
          ! again and booked to the explicit residual.
          ! ------------------------------------------------------------------------- !
          IF (repro_surplus_fate == 1) THEN
            IF (sink_full_carbon_fate == 1) THEN
              cfo_surplus_day = cfo_surplus_day + surplus
              cfo_wsc_day    = cfo_wsc_day   + wsc_label
              cfo_loss_day   = cfo_loss_day  + (surplus - wsc_label)
              cfo_loss_cum   = cfo_loss_cum  + (surplus - wsc_label)
              cfo_wsc_cum    = cfo_wsc_cum   + wsc_label
              IF (surplus - wsc_label > 0.0_wp) cfo_days_loss  = cfo_days_loss  + 1
              IF (wsc_label > 0.0_wp)           cfo_days_store = cfo_days_store + 1
            ELSE
              rsf_resid_day   = rsf_resid_day + (surplus - wsc_label)
              rsf_resid_cum   = rsf_resid_cum + (surplus - wsc_label)
            END IF
            rsf_label_cum   = rsf_label_cum + wsc_label
            rsf_surplus_cum = rsf_surplus_cum + surplus
            IF (surplus - wsc_label > 0.0_wp) rsf_days = rsf_days + 1
            W_stem = W_stem - (surplus - wsc_label)
          END IF
          W_reserve = W_reserve + wsc_label
          overflow_to_reserve_total = overflow_to_reserve_total + wsc_label
          IF (stem_wsc_mech >= 1) THEN
            wsc_dep_post = wsc_dep_post + wsc_label
            wsc_dep_day  = wsc_dep_day + wsc_label
          END IF
        END BLOCK
        W_grain   = W_grain + to_grain
        dW_grain_out = to_grain
        cres_to_grain_assim = to_grain
        room = room - to_grain
        wsc_c_in = wsc_c_in + dW_repro_use    ! assimilate this routine took (rec 117 books its own draw above)
      END IF
      ! ------------------------------------------------------------------------------- !
      ! Driver record 105 mode 2: deficit-independent window release with a
      ! deficit-decided split. Schnyder (1993) New Phytol 123:233-245, p.239, reports
      ! that the store is mobilized when the ears are trimmed or removed as well, so
      ! the release is not driven by the grain deficit; p.240 adds that the sink sets
      ! only what becomes of the released carbon. Case 1 keeps the conversion
      ! efficiency of Gebbing & Schnyder (1999) Plant Physiol 121:871-878, p.877, and
      ! case 2 books a named respiratory fate.
      ! ------------------------------------------------------------------------------- !
      IF (wsc_reserve_law == 2 .AND. carbon_reserve_mode == 1 .AND. &
          dvs2 >= reserve_start_dvs2 .AND. dvs2 <= reserve_end_dvs2 .AND. &
          ( wsc_fate_law == 1 .AND. gfb_end_bound /= 0 ) .AND. &
          W_reserve > 0.0_wp) THEN
        BLOCK
          REAL(wp) :: rel, rel_g, rel_f, gain2, loss_conv, loss_fate
          REAL(wp) :: dtau2, taurem2, cap2, room_pre, resid2
          room_pre = room
          ! Driver record 107: the release-window-end bracket. Both ends are clocks
          ! already read at this site, so no new state and no new constant enter.
            dtau2        = wsl_ddvs2
            taurem2      = reserve_end_dvs2 - dvs2
            wsl2_tau_day = dvs2
          ! RELEASE. taurem2 > dtau2 >= 0 on the division limb, so the quotient cannot
          ! divide by zero; a window already at or past its end releases the whole pool.
          IF (dtau2 >= taurem2) THEN
            cap2 = W_reserve
          ELSE
            cap2 = W_reserve * dtau2 / taurem2
          END IF
          rel = MIN(W_reserve, W_stem, cap2)
          ! which term bound today's release -- REPORTED, read by no flux
          IF (cap2 < MIN(W_reserve, W_stem)) THEN
            wsl2_bound_day = 4
          ELSE IF (W_reserve <= W_stem) THEN
            wsl2_bound_day = 2
          ELSE
            wsl2_bound_day = 3
          END IF
          ! SPLIT, case 1 reuses the room/efficiency bound and the conversion efficiency
          ! of the record-92 draw.
          rel_g     = MIN(rel, room / MAX(1.0e-9_wp, reserve_conv_eff))
          gain2     = rel_g * reserve_conv_eff
          loss_conv = rel_g * (1.0_wp - reserve_conv_eff)
          ! Case 2 applies NO conversion efficiency: the ear-removal fate and the
          ! sink-intact fate are different sink states.
          rel_f     = rel - rel_g
          loss_fate = rel_f
          resid2 = rel - rel_g - rel_f
          IF (ABS(resid2) > wsl2_i1_max) wsl2_i1_max = ABS(resid2)
          resid2 = gain2 + loss_conv + loss_fate - rel
          IF (ABS(resid2) > wsl2_i2_max) wsl2_i2_max = ABS(resid2)
          wsl2_wres_pre  = W_reserve
          wsl2_rel_day   = rel
          wsl2_relg_day  = rel_g
          wsl2_relf_day  = rel_f
          wsl2_gain_day  = gain2
          wsl2_lconv_day = loss_conv
          wsl2_taurem_day = taurem2
          ! STATE. The pool and the stem lose the FULL release; the grain takes only case
          ! 1 times the efficiency; everything else is the named ledger.
          W_reserve    = W_reserve - rel
          W_stem       = W_stem - rel
          W_grain      = W_grain + gain2
          b1_reserve_to_grain = b1_reserve_to_grain + gain2   ! B1: enters the reproductive pool
          dW_grain_out = dW_grain_out + gain2
          cres_to_grain_res = gain2
          cres_respired = rel - gain2
          resid2 = cres_respired - (loss_conv + loss_fate)
          IF (ABS(resid2) > wsl2_i3_max) wsl2_i3_max = ABS(resid2)
          wsl2_cres_day = cres_respired
          room = room - gain2
          ! the TWO streams stay SEPARATE in the totals as well as in the code
          wsl2_rel_tot  = wsl2_rel_tot + rel
          wsl2_relg_tot = wsl2_relg_tot + rel_g
          wsl2_relf_tot = wsl2_relf_tot + rel_f
          wsl2_gain_tot = wsl2_gain_tot + gain2
          wsl2_loss_conv_tot = wsl2_loss_conv_tot + loss_conv
          wsl2_loss_fate_tot = wsl2_loss_fate_tot + loss_fate
          IF (rel > 0.0_wp .AND. room_pre <= 0.0_wp) wsl2_days_zerodef = wsl2_days_zerodef + 1

          ! Driver record 148 term-B diagnostics, read by no flux. The gate day, the
          ! bound and the pool at the gate are latched on the first release.
          IF (wsc_fate_law == 1) THEN
            IF (wfl_gate_day == 0) THEN
              wfl_gate_day     = iday_in
              wfl_gate_bound   = gfb_end_bound
              wfl_pool_at_gate = wsl2_wres_pre
            END IF
            wfl_b_released_day = wfl_b_released_day + rel
            wfl_b_released_cum = wfl_b_released_cum + rel
            IF (rel > 0.0_wp) wfl_b_days = wfl_b_days + 1
            IF (gain2 > 0.0_wp) THEN
              wfl_days_grain_fed = wfl_days_grain_fed + 1
              wfl_grain_fed_cum  = wfl_grain_fed_cum + gain2
            END IF
          END IF

          ! The record-92 pool diagnostics are fed with the release and the gain, so that
          ! the reserve file keeps meaning what its header says on this path.
          IF (stem_wsc_mech >= 1) THEN
            wsc_draw_tot = wsc_draw_tot + rel
            wsc_gain_tot = wsc_gain_tot + gain2
            wsc_resp_tot = wsc_resp_tot + (rel - gain2)
            IF (wsc_day_mob0 == 0) wsc_day_mob0 = iday_in
          END IF
        END BLOCK
      ! ------------------------------------------------------------------------------- !
      ! Reserve remobilization: a demand-driven draw from the stem reserve pool when
      ! today's assimilate left sink room. The drawn carbon times the conversion
      ! efficiency reaches the grain and the remainder is respired. The window opens
      ! only after the accumulation window has closed (Schnyder 1993 New Phytol
      ! 123:233-245, p.238), and the measured internode series behind that timing is
      ! Xue, Waines, Zhu et al. (2006) Crop Sci 46:735-746,
      ! doi:10.2135/cropsci2005.04-0033, Table 2 p.738 with the text p.737, where the
      ! lower internodes peak about ten days after anthesis.
      ! ------------------------------------------------------------------------------- !
      ELSE IF ((wsc_reserve_law /= 2 .OR. wsc_fate_law == 1) .AND. &
          carbon_reserve_mode == 1 .AND. dvs2 >= reserve_start_dvs2 .AND. &
          dvs2 <= reserve_end_dvs2 .AND. room > 0.0_wp .AND. W_reserve > 0.0_wp .AND. &
          (wsc_day_anth > 0 .AND. (iday_in - wsc_day_anth) > wsc_accum_end_dpa)) THEN
        ! ----------------------------------------------------------------------------- !
        ! Driver record 92 modes 1, 3 and 4: the draw is set by the DEFICIT and not by
        ! a rate constant. Schnyder (1993) New Phytol 123:233-245, p.239, has the
        ! fructan pools supplementing the photosynthate supply to the grains, and
        ! Summary p.233 that the pools "provide" that supplement.
        ! ----------------------------------------------------------------------------- !
        IF (stem_wsc_mech == 1 .OR. stem_wsc_mech == 3 .OR. stem_wsc_mech == 4) THEN
          draw = MIN(W_reserve, W_stem, room / MAX(1.0e-9_wp, reserve_conv_eff))
        ELSE
          draw = MIN(reserve_remob_rate * W_reserve, W_reserve, W_stem, &
                     room / MAX(1.0e-9_wp, reserve_conv_eff))
        END IF
        ! ----------------------------------------------------------------------------- !
        ! Driver record 105 mode 1: the post-anthesis-only arm, a bounded draw on the
        ! grain-filling clock rather than on the deficit.
        ! ----------------------------------------------------------------------------- !
        IF (wsc_fate_law == 1) THEN
          wfl_draw_def_day = wfl_draw_def_day + draw
          wfl_draw_def_cum = wfl_draw_def_cum + draw
          IF (draw > 0.0_wp) wfl_draw_def_days = wfl_draw_def_days + 1
        END IF
        gain = draw * reserve_conv_eff
        W_reserve = W_reserve - draw
        W_stem    = W_stem - draw
        W_grain   = W_grain + gain
        b1_reserve_to_grain = b1_reserve_to_grain + gain   ! B1: enters the reproductive pool
        dW_grain_out = dW_grain_out + gain
        cres_to_grain_res = gain
        cres_respired = draw - gain
        ! ----------------------------------------------------------------------------- !
        ! The ledger residual is formed from the two independent terms of the day, so
        ! that it is not an identity in its own arguments.
        ! ----------------------------------------------------------------------------- !
        room = room - gain
        IF (stem_wsc_mech >= 1) THEN
          wsc_draw_tot = wsc_draw_tot + draw
          wsc_gain_tot = wsc_gain_tot + gain
          wsc_resp_tot = wsc_resp_tot + (draw - gain)
          IF (wsc_day_mob0 == 0) wsc_day_mob0 = iday_in
          ! --------------------------------------------------------------------------- !
          ! No published source was located for a 16-day threshold on the release. It is
          ! in neither Gebbing & Schnyder (1999) Plant Physiol 121:871-878, which was
          ! read page by page and whose only two "83 %" figures are on p.874 (the share
          ! of pre-anthesis carbon present as starch in grain carbohydrate) and p.877
          ! (the lower bound of the mobilized reserve), nor in Gebbing, Schnyder &
          ! Kuhbauch (1998) Plant Cell Environ 21:301-313, which sampled only at anthesis
          ! and at maturity. The accumulator below is therefore a DIAGNOSTIC column: it
          ! is read by two WRITE statements and by nothing else.
          ! --------------------------------------------------------------------------- !
          IF (wsc_day_anth > 0 .AND. (iday_in - wsc_day_anth) > 16) &
            wsc_draw_after16 = wsc_draw_after16 + draw
        END IF
      ! ------------------------------------------------------------------------------- !
      ! Driver record 165 diagnostic limb: reached only when every condition of the
      ! deficit limb above holds EXCEPT the accumulation-window delay, so that the
      ! days the delay held the release are counted.
      ! ------------------------------------------------------------------------------- !
      ELSE IF (wsc_release_delay_basis == 1 .AND. &
          (wsc_reserve_law /= 2 .OR. wsc_fate_law == 1) .AND. &
          carbon_reserve_mode == 1 .AND. dvs2 >= reserve_start_dvs2 .AND. &
          dvs2 <= reserve_end_dvs2 .AND. room > 0.0_wp .AND. W_reserve > 0.0_wp) THEN
        wrd_days_delayed  = wrd_days_delayed + 1
        wrd_room_held_cum = wrd_room_held_cum + room
        wrd_pool_held_max = MAX(wrd_pool_held_max, W_reserve)
      END IF
      cres_sink_unfilled = room
      ! ------------------------------------------------------------------------------- !
      ! Driver record 91 ledger. On this path the room is a daily FLUX, so the stock
      ! "sink room left" of the incumbent path has no meaning and is replaced by the
      ! flux form.
      ! ------------------------------------------------------------------------------- !
      IF (grain_fill_mech >= 1) THEN
        gfb_supplied_tot = gfb_supplied_tot + (W_grain - gfb_before)
        ! ----------------------------------------------------------------------------- !
        ! The demand and the unmet demand are FLUX quantities and exist only where the
        ! law is a flux. On driver record 91 mode 3 `room` is still a STOCK, so summing
        ! it would be meaningless and that arm reports them as not defined.
        ! ----------------------------------------------------------------------------- !
        IF (grain_fill_mech /= 3) gfb_unmet_tot = gfb_unmet_tot + MAX(0.0_wp, room)
        ! ----------------------------------------------------------------------------- !
        ! The reserve pool as it stands after THIS day supply, kept only while the law
        ! is live, so that grain_filling_law.csv reports the pool at the END OF FILLING
        ! and not the harvest-day value, which the record-105 release zeroes.
        ! ----------------------------------------------------------------------------- !
        IF (gfb_rate > 0.0_wp) gfb_reserve_end = W_reserve
      END IF
    END IF

    ! Book the day two cross-boundary fluxes on the reproductive pool, once. After
    ! this statement W_storage = W_grain + W_chaff holds on every day.
    IF (b1_repro_to_stem > 0.0_wp .OR. b1_reserve_to_grain > 0.0_wp) &
      W_storage = W_storage - b1_repro_to_stem + b1_reserve_to_grain

    ! ------------------------- The real closure ------------------------- !
    ! Not a tautology: every term is read from the ORGAN STATE captured on entry,
    ! so a deleted or duplicated update fails it.
    IF (stem_wsc_mech >= 1) THEN
      cres_ledger_resid = (W_grain - wsc_wg0) + (W_stem - wsc_ws0) &
                        + (W_chaff - wsc_wc0) + cres_respired - wsc_c_in
      IF (repro_surplus_fate == 1) cres_ledger_resid = cres_ledger_resid + rsf_resid_day   ! rec 127: named destination
      IF (sink_full_carbon_fate == 1) cres_ledger_resid = cres_ledger_resid + cfo_loss_day   ! rec 129: named destination
      IF (ABS(cres_ledger_resid) > wsc_closure_max) wsc_closure_max = ABS(cres_ledger_resid)
      IF (W_reserve > wsc_pool_peak) THEN
        wsc_pool_peak = W_reserve ; wsc_day_peak = iday_in
      END IF
    END IF

    ! Driver record 117: the SAME closure identity recomputed in double precision
    ! on the active path. It is not a tautology for the same reason.
  END SUBROUTINE allocate_yield_components_daily

  ! Driver record 89 modes 1 and 2: hold the measured spike share of Fischer
  ! (1985) J Agric Sci Camb 105:447-461, p.452, through the lamina
  ! renormalization of driver record 69, which would otherwise dilute it.
  SUBROUTINE gnx_hold_spike_share(f_leaf, f_stem, f_root, f_storage, fr, fo, count_it)
    REAL(wp), INTENT(IN)    :: f_leaf, fr, fo
    REAL(wp), INTENT(INOUT) :: f_stem, f_root, f_storage
    LOGICAL,  INTENT(IN)    :: count_it
    REAL(wp) :: room, want

    f_root = MIN(MAX(0.0_wp, fr), MAX(0.0_wp, 1.0_wp - f_leaf))
    room   = MAX(0.0_wp, 1.0_wp - f_leaf - f_root)
    want   = MAX(0.0_wp, fo)
    IF (want > room) THEN
      IF (count_it) THEN
        gnx_clamp_days = gnx_clamp_days + 1
        IF (want - room > gnx_clamp_max) gnx_clamp_max = want - room
      END IF
      f_storage = room
    ELSE
      f_storage = want
    END IF
    f_stem = 1.0_wp - f_leaf - f_root - f_storage
    IF (f_stem < 0.0_wp) f_stem = 0.0_wp
  END SUBROUTINE gnx_hold_spike_share
  ! --------------------------------------------------------------------------------- !
  ! Driver record 89 modes 1 and 3: distribute the day spike dry-matter increment
  ! among the live axis classes, so that ears per square meter and grains per ear
  ! stop being a single crop-level number. The window opens at penultimate-leaf
  ! appearance per axis (Fischer 1985 J Agric Sci Camb 105:447-461, p.449) and
  ! the resulting kernels per ear are compared with the measured 29.6 of that
  ! paper Table 4 p.451.
  ! --------------------------------------------------------------------------------- !
  SUBROUTINE update_spike_window(iday, tair_day, tdu_rem, dtdu_day)
    INTEGER,  INTENT(IN) :: iday
    REAL(wp), INTENT(IN) :: tair_day      ! screen air daily mean, degC (Fischer p.449 Materials and Methods)
    REAL(wp), INTENT(IN) :: tdu_rem       ! thermal-development units still to anthesis (>= 0)
    REAL(wp), INTENT(IN) :: dtdu_day      ! today's thermal-development increment

    gnx_day_now = iday
    gnx_sum_tdu  = gnx_sum_tdu  + MAX(0.0_wp, dtdu_day)
    gnx_sum_cd45 = gnx_sum_cd45 + MAX(0.0_wp, tair_day - gnx_tt_base_C)
    IF (gnx_sum_tdu > 0.0_wp) gnx_tt_ratio = gnx_sum_cd45 / gnx_sum_tdu
    IF (gnx_win_open) RETURN
    IF (gnx_sum_tdu <= 0.0_wp) RETURN
    IF (tdu_rem * gnx_tt_ratio <= gnx_tt_ref_Cd) THEN
      gnx_win_open     = .TRUE.
      gnx_win_open_day = iday
      gnx_tdu_rem_open = tdu_rem
    END IF
  END SUBROUTINE update_spike_window
  ! --------------------------------------------------------------------------------- !
  ! Driver record 98: today's potential spike growth, i.e. the spike DEMAND.
  ! Fischer (1985) J Agric Sci Camb 105:447-461: 32.3 florets per spike (Table 4
  ! p.451) times 10 mg spike dry weight per competent floret (p.449) over a
  ! 300 degree-Cd spike-growth period (p.449), giving the per-axis rate
  ! 1.0767e-3 g axis-1 degree-Cd-1. The base temperature 4.5 degC is that of the
  ! paper's own degree-day statement (p.449, "300 day degrees above 4-5 degC").
  ! --------------------------------------------------------------------------------- !
  SUBROUTINE bgn_set_spike_demand(tair_day)
    REAL(wp), INTENT(IN) :: tair_day

    bgn_demand_day = 0.0_wp
    bgn_axes_day   = 0.0_wp
    IF (spike_sink_mech /= 1 .AND. spike_sink_mech /= 2) RETURN
  END SUBROUTINE bgn_set_spike_demand

  SUBROUTINE update_spike_axes(iday, tair_day, dvs2)
    INTEGER,  INTENT(IN) :: iday
    REAL(wp), INTENT(IN) :: tair_day, dvs2
    INTEGER  :: c
    REAL(wp) :: w(0:CT_MAXTIL), wsum

    IF (.NOT. gnx_win_open) RETURN
    IF (dvs2 >= 1.0_wp) RETURN                 ! the spike growth period ENDS at anthesis

    gnx_win_days   = gnx_win_days + 1
    gnx_win_tt45   = gnx_win_tt45   + MAX(0.0_wp, tair_day - gnx_tt_base_C)   ! REPORTED vs Fischer's 300 degCd
    gnx_win_dagb   = gnx_win_dagb   + gnx_dagb_day
    gnx_win_dspike = gnx_win_dspike + gnx_dspike_day

    wsum = 0.0_wp
    DO c = 0, CT_MAXTIL
      w(c) = 0.0_wp
      IF (.NOT. ct_appeared(c)) CYCLE
      IF (ct_dens(c) <= 0.0_wp) CYCLE
      w(c) = ct_dens(c) * MAX(0.0_wp, ct_area_axis_real(c))   ! model construction
      wsum = wsum + w(c)
    END DO
    IF (wsum > 0.0_wp) THEN
      DO c = 0, CT_MAXTIL
        gnx_pool(c) = gnx_pool(c) + gnx_dspike_day * w(c) / wsum
      END DO
    ELSE
      ! No live axis carries any green lamina area on a day the crop window is open.
      ! Credited to the main stem and counted, never silent.
      gnx_pool(0)  = gnx_pool(0)  + gnx_dspike_day
      gnx_unshared = gnx_unshared + gnx_dspike_day
    END IF
  END SUBROUTINE update_spike_axes
  ! --------------------------------------------------------------------------------- !
  ! Driver record 89 modes 1 and 3: set the grain number at anthesis from the
  ! per-axis spike carbon through the two measured components of Fischer (1985)
  ! J Agric Sci Camb 105:447-461 -
  !   competent florets per spike = spike dry weight / 10 mg   (p.449; Table 4
  !     p.451 measures 102.3 and 103.3 florets g-1 on a second cultivar),
  !   kernels per unit spike dry weight at anthesis = 78 +/- 2 g-1 (p.447),
  ! so that grain number follows the spike mass each axis actually accrued.
  ! --------------------------------------------------------------------------------- !
  SUBROUTINE bgn_floral_limit(w_axis, nfl)
    REAL(wp), INTENT(IN)    :: w_axis
    REAL(wp), INTENT(INOUT) :: nfl

    IF (spike_sink_mech /= 1 .AND. spike_sink_mech /= 3) RETURN
  END SUBROUTINE bgn_floral_limit

  SUBROUTINE gnx_set_grain_number(W_chaff_anth, kernel_number)
    REAL(wp), INTENT(IN)  :: W_chaff_anth
    REAL(wp), INTENT(OUT) :: kernel_number
    INTEGER  :: c
    REAL(wp) :: psum, scale

    gnx_spike_dw_anth = W_chaff_anth

    ! The per-axis pools are a partition of the SAME carbon, so any drift - a spike
    ! increment that fell outside the per-axis window - is measured and reported
    ! rather than renormalized away.
    psum = 0.0_wp
    DO c = 0, CT_MAXTIL
      psum = psum + gnx_pool(c)
    END DO
    IF (psum > 1.0e-12_wp) THEN
      scale = W_chaff_anth / psum
    ELSE
      scale = 0.0_wp
    END IF
    gnx_pool_scale = scale

    gnx_ears_m2 = 0.0_wp ; gnx_florets_m2 = 0.0_wp ; gnx_sterile_C = 0.0_wp
    bgn_floral_bound = 0 ; bgn_ear_spike_dm = 0.0_wp ; bgn_fe_ear_basis = 0.0_wp
    gnx_axes_culled = 0 ; gnx_pool_moved = 0.0_wp ; gnx_cull_passes = 0
    gnx_ninit = gnx_spikelets * gnx_mfs           ! 19.4 x 10.0 = 194 primordia per spike
    kernel_number = 0.0_wp
    DO c = 0, CT_MAXTIL
      gnx_pool(c)     = gnx_pool(c) * scale
      gnx_pool_pre(c) = gnx_pool(c)
      gnx_dens(c)     = ct_dens(c)
      gnx_wsp(c)      = 0.0_wp ; gnx_nfl(c) = 0.0_wp ; gnx_kn(c) = 0.0_wp
      gnx_surv(c)     = 0.0_wp
      gnx_is_ear(c)   = .FALSE.
    END DO
      ! --------------------------------------------------------------------------------- !
      ! ----------- Driver record 90 mode 4: the spike-affordability ear test ----------- !
      ! --------------------------------------------------------------------------------- !
      ! The floret-site chain is run once per axis class and the ear-bearing fraction
      ! follows from the spike carbon the class could afford, inside the measured
      ! window of Fischer (1985) J Agric Sci Camb 105:447-461, p.449. The same paper
      ! p.452 measures that spikeless shoots are real in an unstressed field crop, so
      ! an earless class is an admissible outcome and not an error.
      ! --------------------------------------------------------------------------------- !
      BLOCK
        REAL(wp) :: base, ster, r, phi_new
        INTEGER  :: kill, pass
        gnx_wpot_axis = bgn_r_spike * gnx_win_tt45
        DO c = 0, CT_MAXTIL
          gnx_phi_sp(c) = MERGE(1.0_wp, 0.0_wp, ct_dens(c) > 0.0_wp)
        END DO
        DO pass = 0, CT_MAXTIL + 1
          gnx_cull_passes = pass
          ! (i) the spike carbon of the spikeless axes returns to the ear-bearers, in the
          !     proportion of their own pre-mortality pools.
          ster = 0.0_wp ; base = 0.0_wp
          DO c = 0, CT_MAXTIL
            base = base + gnx_pool_pre(c) * gnx_phi_sp(c)
            ster = ster + gnx_pool_pre(c) * (1.0_wp - gnx_phi_sp(c))
          END DO
          gnx_pool_moved = ster
          DO c = 0, CT_MAXTIL
            IF (gnx_phi_sp(c) <= 0.0_wp) THEN
              gnx_pool(c) = 0.0_wp
            ELSE IF (base > 0.0_wp) THEN
              gnx_pool(c) = gnx_pool_pre(c) * gnx_phi_sp(c) * (1.0_wp + ster / base)
            ELSE
              gnx_pool(c) = gnx_pool_pre(c) * gnx_phi_sp(c)
            END IF
          END DO
          ! (ii) floret survival on the redistributed pools, per EAR-BEARING axis
          DO c = 0, CT_MAXTIL
            gnx_wsp(c) = 0.0_wp ; gnx_surv(c) = 0.0_wp
            gnx_nfl(c) = 0.0_wp ; gnx_kn(c) = 0.0_wp
            IF (gnx_phi_sp(c) <= 0.0_wp) CYCLE
            gnx_wsp(c) = gnx_pool(c) / (ct_dens(c) * gnx_phi_sp(c))   ! g spike DM per ear-bearing axis = W'_c
            r = gnx_wsp(c) / (gnx_ninit * gnx_m_floret * gnx_floret_scale)   ! spike assimilate per primordium, in units of m_floret (rec 151 scale = 1 when OFF)
            ! rec 137: f_N,fl scales the hazard; f_N,fl = 1 ->.
            fnl_r(c) = r ; fnl_surv_c_only(c) = 1.0_wp - EXP(-r)
            gnx_surv(c) = fnl_surv_c_only(c)                          ! constant-hazard saturation
            gnx_nfl(c) = gnx_ninit * gnx_surv(c)                      ! fertile florets per spike
            CALL bgn_floral_limit(gnx_wsp(c), gnx_nfl(c))             ! B-GN P2 (rec 98 = 1 or 3)
            gnx_kn(c)  = gnx_kern_per_fl * gnx_nfl(c)                 ! Fischer p.447 / p.449
          END DO

          ! (iii) THE TEST: spike affordability on the PRE-mortality pool of the class,
          !       youngest first, with the main stem never earless.
          kill = -1
          DO c = CT_MAXTIL, 1, -1
            IF (gnx_phi_sp(c) <= 0.0_wp) CYCLE
            phi_new = 1.0_wp
            IF (gnx_wpot_axis > 0.0_wp) phi_new = MIN(1.0_wp, gnx_pool_pre(c) / ct_dens(c) / gnx_wpot_axis)

            ! Driver record 137: one nitrogen term on grain number and never two, so when
            ! the floret-survival nitrogen term is on, this ear-fraction tier does not also
            ! apply one.
            IF (phi_new < gnx_phi_sp(c)) THEN
              kill = c ; EXIT
            END IF
          END DO
          IF (kill < 0) EXIT
          gnx_phi_sp(kill)  = phi_new
          gnx_axes_culled   = gnx_axes_culled + 1
        END DO

        ! (iv) the fates at anthesis, exhaustive: dead, standing but spikeless (an
        !      ear-bearing fraction of zero included), or ear-bearing.
        gnx_earless_m2 = 0.0_wp
        DO c = 0, CT_MAXTIL
          gnx_sterile_C = gnx_sterile_C + gnx_pool_pre(c) * (1.0_wp - gnx_phi_sp(c))
          gnx_earless_m2 = gnx_earless_m2 + ct_dens(c) * (1.0_wp - gnx_phi_sp(c))
          IF (gnx_phi_sp(c) <= 0.0_wp) CYCLE
          gnx_is_ear(c)  = .TRUE.
          gnx_ears_m2    = gnx_ears_m2    + ct_dens(c) * gnx_phi_sp(c)
          gnx_florets_m2 = gnx_florets_m2 + ct_dens(c) * gnx_phi_sp(c) * gnx_nfl(c)
          kernel_number  = kernel_number  + ct_dens(c) * gnx_phi_sp(c) * gnx_kn(c)
        END DO
      END BLOCK
    gnx_grains_per_ear = 0.0_wp
    IF (gnx_ears_m2 > 0.0_wp) gnx_grains_per_ear = kernel_number / gnx_ears_m2
    gnx_fe_realised = 0.0_wp
    IF (W_chaff_anth > 0.0_wp) gnx_fe_realised = kernel_number / W_chaff_anth
    ! ----------------------------------------------------------------------------------- !
    ! The like-for-like fruiting-efficiency basis. Rivera-Amado et al. (2019) Field
    ! Crops Res 240:154-167, Fig. 2(c), and Slafer et al. (2015) Food Energy Secur
    ! 4(2):92-109, p.95, define fruiting efficiency on the dry weight of the
    ! SAMPLED spikes, so shoots that bear no ear are excluded; Fischer (1985)
    ! p.452 measures those at 5-10 % of the anthesis dry matter.
    ! ----------------------------------------------------------------------------------- !
    DO c = 0, CT_MAXTIL
      IF (gnx_is_ear(c)) bgn_ear_spike_dm = bgn_ear_spike_dm + gnx_pool(c)
    END DO
    IF (bgn_ear_spike_dm > 0.0_wp) bgn_fe_ear_basis = kernel_number / bgn_ear_spike_dm
    bgn_spi_anth = 0.0_wp
    IF (bgn_agb_day > 0.0_wp) bgn_spi_anth = W_chaff_anth / bgn_agb_day
    gnx_share_real = 0.0_wp
    IF (gnx_win_dagb > 0.0_wp) gnx_share_real = gnx_win_dspike / gnx_win_dagb
    gnx_set = .TRUE.
    IF (kernel_number /= kernel_number .OR. kernel_number > HUGE(1.0_wp)) THEN
      WRITE (*,'(A)') ' ERROR rec 89: predicted grain number is not finite'
      ERROR STOP 'rec 89: non-finite grain number'
    END IF
  END SUBROUTINE gnx_set_grain_number

  SUBROUTINE partition_fractions(dvs2, f_leaf, f_stem, f_root, f_storage)
    REAL(wp), INTENT(IN) :: dvs2
    REAL(wp), INTENT(OUT) :: f_leaf, f_stem, f_root, f_storage
    ! ----------------------------------------------------------------------------------- !
    ! Dry-matter partitioning against development stage.
    ! ----------------------------------------------------------------------------------- !
    ! PROVENANCE: these node values are another model's parameter tables, not a
    ! measurement of this experiment. They are the WOFOST winter-wheat
    ! partitioning tables (variety Winter_wheat_106: FRTB roots, FLTB leaves,
    ! FSTB stems, FOTB storage organs) expressed on a TOTAL dry matter basis,
    ! f_organ = (1 - FR) * F_organ(above-ground). Teh (2006) Listings 8.17-8.20,
    ! p.214, tabulates that conversion, and the leaf, stem and root nodes below
    ! are read off those listings:
    !   Listing 8.17 partlv.txt, 8.18 partst.txt, 8.19 partrt.txt, 8.20 partso.txt.
    ! The spike nodes alone are replaced further down by a measured window.
    !
    ! Root-table endpoint: afgen clamps at the last node rather than
    ! extrapolating, so the tabulated root-allocation fraction is held at 0.01
    ! for development stages above 1.10. No claim is made here about the source
    ! tables beyond that stage.
    ! ----------------------------------------------------------------------------------- !
    REAL(wp) :: x_leaf(7), y_leaf(7)
    INTEGER  :: n_leaf
    REAL(wp), PARAMETER :: x_stem(7) = [0.00_wp, 0.10_wp, 0.25_wp, 0.50_wp, 0.646_wp, 0.95_wp, 1.05_wp]
    REAL(wp), PARAMETER :: y_stem(7) = [0.175_wp, 0.175_wp, 0.198_wp, 0.435_wp, 0.6397_wp, 0.975_wp, 0.0_wp]
    REAL(wp), PARAMETER :: x_root(12) = [ &
      0.00_wp, 0.10_wp, 0.20_wp, 0.35_wp, 0.40_wp, 0.50_wp, &
      0.60_wp, 0.70_wp, 0.80_wp, 0.90_wp, 1.00_wp, 1.10_wp ]
    REAL(wp), PARAMETER :: y_root(12) = [ &
      0.5_wp, 0.5_wp, 0.4_wp, 0.22_wp, 0.17_wp, 0.13_wp, &
      0.10_wp, 0.07_wp, 0.05_wp, 0.03_wp, 0.02_wp, 0.01_wp ]
    ! ----------------------------------------------------------------------------------- !
    ! Spike allocation starts in the MEASURED spike-growth window and not at a
    ! nominal development stage. The window opens 30 days before anthesis (Fischer
    ! 1985 J Agric Sci Camb 105:447-461, p.447, where kernel number is set over
    ! that period) and the chaff mass at maturity against which the node height is
    ! set is the measured value relayed by Slafer et al. (2015) Food Energy Secur
    ! 4(2):92-109, p.96.
    ! ----------------------------------------------------------------------------------- !
    REAL(wp) :: x_storage(5), y_storage(5)

    x_storage = [0.00_wp, storage_onset_dvs2, 1.00_wp, 1.05_wp, 2.00_wp]
    y_storage = [0.0_wp, 0.0_wp, storage_peak_frac, 0.98_wp, 0.98_wp]

    x_leaf(1:6) = [0.00_wp, 0.10_wp, 0.25_wp, 0.50_wp, 0.646_wp, 0.95_wp]
    y_leaf(1:6) = [0.325_wp, 0.325_wp, 0.462_wp, 0.435_wp, fleaf_elong, fleaf_heading]
    IF (fleaf_end_dvs > 0.95_wp + 1.0e-9_wp) THEN
      n_leaf = 7; x_leaf(7) = fleaf_end_dvs; y_leaf(7) = 0.0_wp   ! closure node (taper)
    ELSE
      n_leaf = 6                                                   ! baseline shape, exact
    END IF

    f_leaf = afgen(dvs2, x_leaf(1:n_leaf), y_leaf(1:n_leaf))
    f_stem = afgen(dvs2, x_stem, y_stem)
    f_root = afgen(dvs2, x_root, y_root)
    f_storage = afgen(dvs2, x_storage, y_storage)
    ! ----------------------------------------------------------------------------------- !
    ! Driver record 89 modes 1 and 2: the spike partition comes from measurement.
    ! Both tabled nodes above are replaced - the onset by penultimate-leaf
    ! appearance (Fischer 1985 J Agric Sci Camb 105:447-461, p.449) and the height
    ! by the measured partitioning coefficient 0.44 (p.452).
    ! ----------------------------------------------------------------------------------- !
    IF (grain_number_mech == 1 .OR. grain_number_mech == 2) THEN
      IF (dvs2 < 1.0_wp) THEN
        IF (gnx_win_open) THEN
          f_storage = gnx_spike_part * (1.0_wp - f_root)
        ELSE
          f_storage = 0.0_wp
        END IF
      ELSE IF (dvs2 < 1.05_wp) THEN
        BLOCK
          REAL(wp) :: f_anth
          f_anth = gnx_spike_part * (1.0_wp - f_root)
          f_storage = f_anth + (0.98_wp - f_anth) * (dvs2 - 1.0_wp) / 0.05_wp
        END BLOCK
      END IF
    END IF

  END SUBROUTINE partition_fractions

  REAL(wp) FUNCTION afgen(x, xp, yp)
    REAL(wp), INTENT(IN) :: x
    REAL(wp), INTENT(IN) :: xp(:), yp(:)
    INTEGER :: i, n

    n = SIZE(xp)
    IF (x <= xp(1)) THEN
      afgen = yp(1)
      RETURN
    END IF
    IF (x >= xp(n)) THEN
      afgen = yp(n)
      RETURN
    END IF

    DO i = 1, n - 1
      IF (x <= xp(i + 1)) THEN
        afgen = yp(i) + (yp(i + 1) - yp(i)) * (x - xp(i)) / (xp(i + 1) - xp(i))
        RETURN
      END IF
    END DO
    afgen = yp(n)
  END FUNCTION afgen
  ! ------------------------------------------------------------------------------------- !
  ! Per-culm leaf-rank chronology (driver record 64). It writes only the rank
  ! chronology and the Haun stage; no carbon, nitrogen or water state is written
  ! here. The clock is cumulative degree-Cd with base 0 C from emergence, the
  ! definition of McMaster, Wilhelm, Palic, Porter & Jamieson (2003) Ann Bot
  ! 91:697-705, Eqs. 1-2 p.699, doi:10.1093/aob/mcg074.
  ! ------------------------------------------------------------------------------------- !
  SUBROUTINE update_canopy_leaf_ranks(istage_idx, dtt_day, iday)
    INTEGER,  INTENT(IN) :: istage_idx, iday
    REAL(wp), INTENT(IN) :: dtt_day

    REAL(wp) :: dtau, hs_prev, rate_P, fC
    INTEGER  :: n, n_prev, n_app_before
    ! ----------------------------------------------------------------------------------- !
    ! Pre-emergence: the tau clock is defined FROM emergence (Abichou 2016 p.58),
    ! so nothing accumulates and nothing is emitted.
    ! ----------------------------------------------------------------------------------- !
    IF (istage_idx < 1) RETURN
    ! The emergence day itself is tau = 0: start the clock, add no increment.
    IF (.NOT. clr_active) THEN
      clr_active = .TRUE.
      ! --------------------------------------------------------------------------------- !
      ! Driver record 174: the culm final leaf number is predicted here, on the
      ! emergence day, from the continuous solution of the same two-limb
      ! leaf-appearance law.
      ! --------------------------------------------------------------------------------- !
      RETURN
    END IF

    dtau = MAX(0.0_wp, dtt_day)
    clr_tau_Cd = clr_tau_Cd + dtau

    ! --- PRIMARY arm: McMaster degree-Cd phyllochron + Streck f(C -------- !
    hs_prev = clr_HS
    IF (clr_HS < clr_fC_HS_ref) THEN
      ! --------------------------------------------------------------------------------- !
      ! Leaves 1 and 2 run faster on the seed reserves (Streck, Weiss & Baenziger
      ! 2003 Ann Bot 92:181-190, pp.183-184, doi:10.1093/aob/mcg121: 0.32 against
      ! 0.24 leaves d-1), transferred as the RATIO on the phyllochron so that the
      ! faster limb stays outside the Eq. 5 domain (Haun stage >= 2).
      ! --------------------------------------------------------------------------------- !
      rate_P = clr_phyllochron_Cd * clr_seed_ratio
      clr_HS = clr_HS + dtau / rate_P
    ELSE
      fC = (clr_HS / clr_fC_HS_ref) ** clr_fC_exp
      clr_HS = clr_HS + dtau * fC / clr_phyllochron_Cd
    END IF
    ! ----------------------------------------------------------------------------------- !
    ! A Haun-stage crossing IS a LIGULATION under the convention of Abichou (2016)
    ! PhD thesis (NNT tel-02316154) Eq. (1), HS = n + l/L with n the number of
    ! ligulated leaves. The tip time precedes it by the leaf elongation duration.
    ! ----------------------------------------------------------------------------------- !
    n_prev = INT(hs_prev)
    DO n = n_prev + 1, MIN(INT(clr_HS), CLR_MAXRANK)
      IF (clr_tau_lig(n) < 0.0_wp) THEN
        clr_tau_lig(n) = clr_tau_Cd
        clr_day_lig(n) = iday
        clr_tau_app(n) = MAX(0.0_wp, clr_tau_Cd - clr_LED_Cd)
      END IF
    END DO
    clr_n_lig = MIN(INT(clr_HS), CLR_MAXRANK)
    ! ----------------------------------------------------------------------------------- !
    ! The alternative appearance intervals of Evers (2006) are cited through the
    ! secondary literature rather than from that thesis directly, and the arm that
    ! uses them is kept as a comparison chronology.
    ! ----------------------------------------------------------------------------------- !
    hs_prev = clr_HS_ev
    IF (clr_HS_ev < clr_evers_break_HS) THEN
      clr_HS_ev = clr_HS_ev + dtau / clr_P_evers_low
    ELSE
      clr_HS_ev = clr_HS_ev + dtau / clr_P_evers_high
    END IF
    n_prev = INT(hs_prev)
    DO n = n_prev + 1, MIN(INT(clr_HS_ev), CLR_MAXRANK)
      IF (clr_tau_lig_ev(n) < 0.0_wp) THEN
        clr_tau_lig_ev(n) = clr_tau_Cd
        clr_tau_app_ev(n) = MAX(0.0_wp, clr_tau_Cd - clr_LED_Cd)
      END IF
    END DO
    clr_n_lig_ev = MIN(INT(clr_HS_ev), CLR_MAXRANK)

    ! --------------- tipped-but-not-yet-ligulated ranks ------------------------------- !
    ! A rank tips one leaf elongation duration before it ligulates, so the tipped
    ! count runs ahead of the ligulated count by that interval.
    n_app_before = clr_n_app
    clr_n_app = clr_n_lig
    DO n = clr_n_lig + 1, CLR_MAXRANK
      IF (clr_tau_Cd >= clr_tau_lig_proj(n) - clr_LED_Cd) THEN
        clr_n_app = n
        clr_day_app(n) = iday
      ELSE
        EXIT
      END IF
    END DO
    ! ----------------------------------------------------------------------------------- !
    ! Demographic limiter (iii): while the lamina carbon backlog is open, one leaf
    ! is removed from the emergent final leaf number, to a maximum of two. Evers
    ! (2006) PhD thesis, doi:10.18174/24302, Ch.3 Fig. 4, measures a 75 % shade
    ! treatment cutting the main-stem final leaf number from 9 to 8.
    ! ----------------------------------------------------------------------------------- !
    clr_n_app_ev = clr_n_lig_ev
    DO n = clr_n_lig_ev + 1, CLR_MAXRANK
      IF (clr_HS_ev + clr_LED_Cd / clr_P_evers_high >= REAL(n, wp)) THEN
        clr_n_app_ev = n
      ELSE
        EXIT
      END IF
    END DO
    ! ----------------------------------------------------------------------------------- !
    ! --- EMERGENT final leaf number, read off at the anchored booting stage --- !
    ! Zadoks, Chang & Konzak (1974) Weed Res 14:415-421 define the flag-leaf-ligule
    ! stage, and Abichou (2016) anchors the final leaf number at flag ligulation in
    ! the same way.
    ! ----------------------------------------------------------------------------------- !
    IF (istage_idx >= 4 .AND. .NOT. clr_frozen) THEN
      clr_frozen      = .TRUE.
      clr_tau_boot    = clr_tau_Cd
      clr_FLN         = clr_n_lig
      clr_FLN_tip     = clr_n_app
      clr_FLN_ev      = clr_n_lig_ev
      clr_FLN_ev_tip  = clr_n_app_ev
      IF (clr_FLN >= 1) THEN
        clr_tau_lig_flag = clr_tau_lig(clr_FLN)
        clr_tau_app_flag = clr_tau_app(clr_FLN)
      END IF
      ! -------------------------------------------------------------------------------- !
      ! Driver record 174: booting CONFIRMS the prediction, it never rewrites it. The
      ! emergent final leaf number is still read off exactly as before.
      ! -------------------------------------------------------------------------------- !      
    END IF
  END SUBROUTINE update_canopy_leaf_ranks

  ! Projected ligulation thermal time of rank n on today's local phyllochron: exact
  ! for ranks that have already ligulated, a forward projection otherwise.
  REAL(wp) FUNCTION clr_tau_lig_proj(n)
    INTEGER, INTENT(IN) :: n
    REAL(wp) :: p_local
    IF (n >= 1 .AND. n <= CLR_MAXRANK) THEN
      IF (clr_tau_lig(n) >= 0.0_wp) THEN
        clr_tau_lig_proj = clr_tau_lig(n)
        RETURN
      END IF
    END IF
    IF (clr_HS < clr_fC_HS_ref) THEN
      p_local = clr_phyllochron_Cd * clr_seed_ratio
    ELSE
      p_local = clr_phyllochron_Cd * (clr_HS / clr_fC_HS_ref) ** (-clr_fC_exp)
    END IF
    clr_tau_lig_proj = clr_tau_Cd + (REAL(n, wp) - clr_HS) * p_local
  END FUNCTION clr_tau_lig_proj
  ! ----------------------------------------------------------------------------------- !
  ! Driver record 174: the CONTINUOUS solution of the same two-limb
  ! leaf-appearance law update_canopy_leaf_ranks integrates, with the phyllochron
  ! of McMaster, Wilhelm, Palic, Porter & Jamieson (2003) Ann Bot 91:697-705,
  ! Table 2 p.701 (air basis), doi:10.1093/aob/mcg074, and the
  ! appearance-rate limb and 0.75 seed-reserve ratio of Streck, Weiss &
  ! Baenziger (2003) Ann Bot 92:181-190, Eq. 5 p.183, doi:10.1093/aob/mcg121.
  ! ----------------------------------------------------------------------------------- !
  REAL(wp) FUNCTION clr_HS_at_tau(tau)
    REAL(wp), INTENT(IN) :: tau
    REAL(wp) :: a_exp, r_ref, p_ref, tau1

    clr_HS_at_tau = 0.0_wp
    IF (tau <= 0.0_wp) RETURN
    a_exp = -clr_fC_exp
    r_ref = clr_fC_HS_ref
    p_ref = clr_phyllochron_Cd
    IF (p_ref <= 0.0_wp) RETURN
    tau1 = r_ref * p_ref * clr_seed_ratio
    IF (tau <= tau1) THEN
      clr_HS_at_tau = tau / (p_ref * clr_seed_ratio)
    ELSE
      clr_HS_at_tau = ( r_ref ** (a_exp + 1.0_wp) &
                        + (a_exp + 1.0_wp) * r_ref ** a_exp * (tau - tau1) / p_ref ) &
                      ** (1.0_wp / (a_exp + 1.0_wp))
    END IF
  END FUNCTION clr_HS_at_tau

  ! Driver record 174: how many blades does this culm carry? One answer, used by
  ! every routine that loops over ranks, so that the chronology, the areas and
  ! the senescence cannot disagree about the top rank.
  INTEGER FUNCTION cfp_n_top_culm()
      cfp_n_top_culm = MERGE(clr_FLN, cba_Nl_fit, clr_frozen .AND. clr_FLN >= 1)
  END FUNCTION cfp_n_top_culm

  ! Potential and realized blade area per main-stem rank (driver record 65). It
  ! writes only the per-rank areas; no carbon, nitrogen or water state is written
  ! here.
  SUBROUTINE update_canopy_blade_area()
    INTEGER  :: n, irank, n_top
    REAL(wp) :: frac, tlig, tapp

    IF (.NOT. clr_active) RETURN

    cba_A_shoot_pot = 0.0_wp
    cba_A_shoot     = 0.0_wp
    cba_n_expanding = 0

    ! A culm carries as many blades as its final leaf number and no more. The
    ! chronology deliberately keeps emitting ranks past the flag leaf, because it
    ! has no terminal-spikelet signal of its own, so the top rank is imposed here.
    n_top = cfp_n_top_culm()
    DO n = 1, MIN(MIN(MAX(clr_n_app, clr_n_lig), n_top), CLR_MAXRANK)

      ! --- potential dimensions (final, rank-fixed) --------------------- !
      ! Both dimension profiles are indexed by INVERSED rank, counted down from the
      ! flag leaf, which is the convention of the measured profiles.
      IF (clr_frozen .AND. clr_FLN >= 1) THEN
        irank = clr_FLN - n + 1
      ELSE
        irank = cba_Nl_fit - n + 1
      END IF
      irank = MIN(MAX(irank, 1), SIZE(cba_w_irank))
      cba_L_pot(n) = cba_l_irank(irank)   ! Abichou Fig. 8a doc p.104, NI, 5-cultivar median
      cba_W_pot(n) = cba_w_irank(irank)   ! Abichou Fig. 3e doc p.67, Soissons panel central
      ! Form factor by ABSOLUTE rank (Dornbusch Fig. 3A breakpoints).
      cba_ft_rank(n) = cba_form_factor(n, n_top)
      ! Driver record 144: a carbon-coupled bonus on the potential ceiling only. See
      ! Part 1 Table 6; this is an inactive alternative.
      cba_A_pot(n) = cba_L_pot(n) * cba_W_pot(n) * cba_ft_rank(n) * cbc_f(0, n)

      ! --- realized area: linear expansion over the leaf elongation window --- !
      ! A blade expands from its tip appearance to its ligulation; the kinetic is
      ! taken as linear in thermal time over that window.
      tlig = clr_tau_lig_proj(n)
      tapp = MAX(0.0_wp, tlig - clr_LED_Cd)
      IF (clr_LED_Cd > 0.0_wp) THEN
        frac = (clr_tau_Cd - tapp) / clr_LED_Cd
      ELSE
        frac = 1.0_wp
      END IF
      IF (frac < 0.0_wp) frac = 0.0_wp
      IF (frac > 1.0_wp) frac = 1.0_wp
      cba_A(n) = cba_A_pot(n) * frac
      IF (frac > 0.0_wp .AND. frac < 1.0_wp) cba_n_expanding = cba_n_expanding + 1
      cba_A_shoot_pot = cba_A_shoot_pot + cba_A_pot(n)
      cba_A_shoot     = cba_A_shoot     + cba_A(n)
    END DO
  END SUBROUTINE update_canopy_blade_area

  ! Final blade length at a possibly NON-INTEGER absolute rank, by linear
  ! interpolation on the inversed-rank median profile of Abichou (2016) Ch.3
  ! Fig. 8a. Non-integer ranks are required because the tiller dimension shifts
  ! of Ch.2 Table 2 Eqn T2.7b and Table 4 are decimal.
  REAL(wp) FUNCTION cba_L_at(xabs, n_max)
    REAL(wp), INTENT(IN) :: xabs
    INTEGER,  INTENT(IN) :: n_max
    REAL(wp) :: ir, w
    INTEGER  :: i0, i1
    ir = REAL(n_max, wp) - xabs + 1.0_wp          ! inversed rank, decimal
    IF (ir <= 1.0_wp) THEN
      cba_L_at = cba_l_irank(1)
    ELSE IF (ir >= REAL(SIZE(cba_l_irank), wp)) THEN
      cba_L_at = cba_l_irank(SIZE(cba_l_irank))
    ELSE
      i0 = INT(ir) ; i1 = i0 + 1 ; w = ir - REAL(i0, wp)
      cba_L_at = cba_l_irank(i0) * (1.0_wp - w) + cba_l_irank(i1) * w
    END IF
  END FUNCTION cba_L_at

  ! Wheat form factor by absolute rank; Dornbusch, Watt, Baccar, Fournier &
  ! Andrieu (2011) Ann Bot 107:865-873, doi:10.1093/aob/mcq181, Fig. 3A, cv.
  ! Soissons, piecewise between the printed rank anchors.
  SUBROUTINE update_canopy_tillering(lai_now, dens_plant)
    REAL(wp), INTENT(IN) :: lai_now      ! the model's own green LAI, m2 m-2
    REAL(wp), INTENT(IN) :: dens_plant   ! plants m-2 (driver rec 8)
    INTEGER  :: c, nl_ms
    REAL(wp) :: t_open, t_close
    REAL(wp) :: nl_here

    IF (.NOT. clr_active) RETURN
    ct_active = .TRUE.

    ! Driver record 89: the DEVELOPMENTAL form of the spike-growth period of
    ! Fischer (1985) J Agric Sci Camb 105:447-461, p.449, "from penultimate leaf
    ! appearance until anthesis".
    IF (grain_number_mech >= 1 .AND. gnx_penult_day < 0) THEN
      BLOCK
        INTEGER :: n_top_gnx
        n_top_gnx = cfp_n_top_culm()
        IF (n_top_gnx >= 2 .AND. clr_HS >= REAL(n_top_gnx, wp) - 1.0_wp) gnx_penult_day = gnx_day_now
      END BLOCK
    END IF

    ! --- the like-for-like construction of the measured light interception: Beer
    ! on the model's own green leaf area index with the model's own extinction
    ! coefficient.
    ct_fpar = 1.0_wp - EXP(-ct_kL * MAX(0.0_wp, lai_now))

    ! Driver record 82: the assimilate surplus the birth gate reads is the PREVIOUS
    ! day dry-matter gain minus the lamina carbon that day built area cost.
    IF (tiller_carbon_gate == 1) THEN
      tcg_csur_day = MAX(0.0_wp, csl_dm_gain_day - csl_LMA * csl_dA_day)
    END IF

    ! Main stem: one axis per plant, present from emergence and never culled here.
    ! The mortality phase of Abichou (2016) regresses TILLERS, and the ear-bearing
    ! main stem is the floor of its Eqn T4.4.
    IF (.NOT. ct_appeared(0)) THEN
      ct_appeared(0)   = .TRUE.
      ct_dens_born(0)  = dens_plant
      ct_dens(0)       = dens_plant
      ct_tau_birth(0)  = clr_tau_Cd
    END IF

    nl_ms = cfp_n_top_culm()
    ct_nl_axis(0) = REAL(nl_ms, wp)

    DO c = 1, CT_MAXTIL
      ! --- window of opportunity for bud outgrowth, Evers et al. (2006) Fig. 5A:
      !     from the ligule appearance of the bud PARENT leaf (rank c) to the ligule
      !     appearance of the next leaf. The rank chronology supplies both.
      t_open  = clr_tau_lig_proj(c)
      t_close = clr_tau_lig_proj(c + 1)
      IF (.NOT. ct_decided(c)) THEN
        IF (clr_tau_Cd >= t_open .AND. clr_tau_Cd < t_close) THEN
          ct_win_sum(c) = ct_win_sum(c) + ct_fpar
          ct_win_n(c)   = ct_win_n(c) + 1.0_wp
          IF (tiller_carbon_gate == 1) THEN         ! D2-1: SUM_window C_sur
            tcg_sup(c)  = tcg_sup(c) + tcg_csur_day
            tcg_nwin(c) = tcg_nwin(c) + 1.0_wp
          END IF
        ELSE IF (clr_tau_Cd >= t_close) THEN
          IF (ct_win_n(c) > 0.0_wp) THEN
            ct_pint_win(c) = ct_win_sum(c) / ct_win_n(c)
          ELSE
            ct_pint_win(c) = ct_fpar
          END IF
          ! p = 1 below 0.40, linear to 0 at 0.55 and 0 beyond, the same interpolation
          ! the nitrogen status uses in SHARED.f90.
          IF (ct_pint_win(c) <= ct_pint_p1) THEN
            ct_p_occ(c) = 1.0_wp
          ELSE IF (ct_pint_win(c) >= ct_pint_p0) THEN
            ct_p_occ(c) = 0.0_wp
          ELSE
            ct_p_occ(c) = (ct_pint_p0 - ct_pint_win(c)) / (ct_pint_p0 - ct_pint_p1)
          END IF
          ct_decided(c) = .TRUE.
        END IF
      END IF

      ! Demographic limiter (ii): while the backlog is open, a cohort that has not
      ! yet produced a ligulated leaf is not continued. Cessation is checked here, on
      ! live axes, because that is where the measured cessation of Abichou (2016)
      ! Ch.1 Table 6 is defined.
      IF (.NOT. ct_appeared(c) .AND. ct_decided(c) .AND. clr_HS >= ct_hs_emerg(c)) THEN
        IF (ct_p_occ(c) <= 0.0_wp) THEN
          IF (.NOT. ct_ceased) THEN
            ct_ceased = .TRUE. ; ct_cess_why = 1     ! Evers light law
          END IF
        ! Cultivar TILLERING CAPACITY (axes per plant times plant density, axes m-2).
        ! Zero means no cap, so the statement is inert when the row is absent.
        ELSE IF (Axes_max_per_plant > 0.0_wp .AND. &
                 ct_axes_live >= Axes_max_per_plant * dens_plant) THEN
          IF (.NOT. ct_ceased) THEN
            ct_ceased = .TRUE. ; ct_cess_why = 3     ! cultivar tillering capacity
          END IF
        END IF

        ! Demographic limiter (i): while the lamina carbon backlog is open, no new
        ! cohort is opened. Evers (2006) PhD thesis, doi:10.18174/24302, Ch.4 Table 2,
        ! measures the same 75 % shade driving the probability of primary-tiller
        ! occurrence down to 0-0.3.
        IF (.NOT. ct_ceased) THEN
          ct_appeared(c)  = .TRUE.
          IF (tiller_carbon_gate == 1) THEN

            ! Driver record 82: supply-ratio birth. The red:far-red law above is TIMING
            ! only; the density ramp is Abichou (2016) Ch.1 Eqn T3.3, in which leaf x
            ! ligulates at x over the appearance rate after birth.
            BLOCK
              REAL(wp) :: dwin, per_leaf, nlt, apot_sum, lx, wx
              INTEGER  :: x, nt
              dwin = t_close - t_open
              nlt  = ct_nl_tiller(c, nl_ms)
              nt   = MAX(1, NINT(nlt))
              per_leaf = (clr_tau_lig_proj(nl_ms) - clr_tau_Cd) / nlt
              apot_sum = 0.0_wp ; tcg_nrank(c) = 0
              DO x = 1, nt
                IF (x > 1 .AND. REAL(x, wp) * per_leaf - clr_LED_Cd >= dwin) EXIT
                lx = ct_L_tiller(REAL(x, wp), c, nl_ms)
                wx = ct_W_tiller(REAL(x, wp), c, nl_ms)
                apot_sum = apot_sum + lx * wx * cba_form_factor(x, nt)
                tcg_nrank(c) = tcg_nrank(c) + 1
              END DO
              tcg_dem(c) = dens_plant * apot_sum * 1.0e-4_wp * csl_LMA
              IF (tcg_dem(c) > 0.0_wp) THEN
                tcg_ratio(c) = MIN(1.0_wp, tcg_sup(c) / tcg_dem(c))
              ELSE
                tcg_ratio(c) = 1.0_wp
              END IF
            END BLOCK
            ct_dens_born(c) = dens_plant * tcg_ratio(c)
          ELSE
            ct_dens_born(c) = dens_plant * ct_p_occ(c)
          END IF
          ! Truncate this cohort to the room left under the cultivar tillering capacity
          ! (zero means no cap: the statement is skipped).
          IF (Axes_max_per_plant > 0.0_wp) THEN
            ct_dens_born(c) = MIN(ct_dens_born(c), &
                                  MAX(0.0_wp, Axes_max_per_plant * dens_plant - ct_axes_live))
          END IF
          ct_dens(c)      = ct_dens_born(c)
          ct_tau_birth(c) = clr_tau_Cd
          ! final leaf number of this cohort, Eqn T3.5 (Abichou Ch.1 Table 3)
          ct_nl_axis(c) = ct_nl_alpha * REAL(nl_ms, wp) &
                        - ct_nl_beta  * REAL(c + ct_nc_offset, wp)
          IF (ct_nl_axis(c) < 1.0_wp) ct_nl_axis(c) = 1.0_wp
        END IF
      END IF
    END DO

    ! --- live axes and the plateau Namax -------------------------------- !
    ct_axes_live = 0.0_wp ; ct_axes_born = 0.0_wp
    DO c = 0, CT_MAXTIL
      ct_axes_live = ct_axes_live + ct_dens(c)
      ct_axes_born = ct_axes_born + ct_dens_born(c)
    END DO
    IF (ct_axes_live > ct_axes_max) ct_axes_max = ct_axes_live
    IF (tiller_carbon_gate == 1) THEN
      ! The axis density is reported against the envelope of Abichou (2016), 600-800
      ! axes m-2, and the die fraction against a reported ear envelope; the only hard
      ! failure is a non-finite density.
      IF (ct_axes_live > tcg_axes_env_lo) tcg_env600_days = tcg_env600_days + 1
      IF (ct_axes_live > tcg_axes_env_hi) tcg_env800_days = tcg_env800_days + 1
      tcg_ne_env_lo = ct_dens_born(0) + (ct_axes_max - ct_dens_born(0)) * (1.0_wp - tcg_die_frac_hi)
      tcg_ne_env_hi = ct_dens_born(0) + (ct_axes_max - ct_dens_born(0)) * (1.0_wp - tcg_die_frac_lo)
      IF (ct_axes_live /= ct_axes_live .OR. ct_axes_live > HUGE(1.0_wp)) THEN
        WRITE (*,'(A)') ' ERROR rec 82: tiller density is not finite'
        ERROR STOP 'rec 82: non-finite axis density'
      END IF
    END IF

    ! --- mortality phase (Abichou 2016 Ch.1 Table 4) --------------------- !
    ! Start: Eqn T4.2, a main-stem Haun stage linear in the main stem final leaf
    ! number.
    nl_here = REAL(nl_ms, wp)
    IF (nl_here <= ct_nl_sm_lo) THEN
      ct_hs_sm = ct_hs_sm_lo
    ELSE IF (nl_here >= ct_nl_sm_hi) THEN
      ct_hs_sm = ct_hs_sm_hi
    ELSE
      ct_hs_sm = ct_hs_sm_lo + (ct_hs_sm_hi - ct_hs_sm_lo) &
                 * (nl_here - ct_nl_sm_lo) / (ct_nl_sm_hi - ct_nl_sm_lo)
    END IF
    IF (ct_tau_sm < 0.0_wp .AND. clr_HS >= ct_hs_sm) ct_tau_sm = clr_tau_Cd
    ! End: Eqn T4.3, the ligulation of the flag leaf of the main stem.
    ct_tau_em = clr_tau_lig_proj(nl_ms)

    ! --- per-axis blade areas and the GROUND-AREA index ------------------ !
    ! The per-axis columns acquire a ground-area basis here, through the cohort
    ! densities.
    ct_bai_pot = 0.0_wp ; ct_bai = 0.0_wp
    IF (canopy_blade_area == 1) THEN
      ct_area_axis(0)      = cba_A_shoot_pot
      ct_area_axis_real(0) = cba_A_shoot
      ct_area_axis_st(0)   = 0.0_wp
      DO c = 1, CT_MAXTIL
        CALL ct_tiller_area(c, nl_ms, ct_area_axis(c), ct_area_axis_real(c), ct_area_axis_st(c))
      END DO
      DO c = 0, CT_MAXTIL
        ct_bai_pot = ct_bai_pot + ct_dens(c) * ct_area_axis(c)      * 1.0e-4_wp
        ct_bai     = ct_bai     + ct_dens(c) * ct_area_axis_real(c) * 1.0e-4_wp
      END DO
    END IF

    ! Driver record 99: tiller survival and the fertile-ear population, by
    ! per-axis light competition over the measured mortality phase of Abichou
    ! (2016) Ch.1 Table 4.
    IF (tiller_survival_mech == 4 .AND. canopy_blade_area == 1 .AND. &
             ct_tau_sm >= 0.0_wp .AND. ct_tau_em > ct_tau_sm .AND. .NOT. tca_frozen) THEN
      BLOCK
        INTEGER  :: cc, ii, jj, ord(0:CT_MAXTIL)
        REAL(wp) :: lab, rem, take, sup_ax, dem_ax, da_ax, f_abs(0:CT_MAXTIL)
        ! the supply shared out today is yesterday's growth_supply.
        tca_gsup_prev = tca_gsup_day
        IF (clr_tau_Cd >= ct_tau_em) THEN

          ! The window has closed: the shape factor becomes one and the matured estimate
          ! of the eventual deaths is removed in full, youngest first.
          DO cc = 0, CT_MAXTIL
            IF (tca_dens_pre(cc) <= 0.0_wp .AND. ct_dens(cc) > 0.0_wp) &
              tca_dens_pre(cc) = ct_dens(cc)
          END DO
          tca_taun    = 1.0_wp
          tca4_F      = 1.0_wp
          tca4_target = tca4_F * tca4_U
          IF (tca4_target > tca_unaff) tca_unaff = tca4_target
          rem = tca_unaff ; tca_dead = 0.0_wp
          DO cc = CT_MAXTIL, 1, -1
            take = MIN(tca_dens_pre(cc), rem)
            IF (take > tca_dead_c(cc)) tca_dead_c(cc) = take
            ct_dens(cc) = tca_dens_pre(cc) - tca_dead_c(cc)
            tca_dead    = tca_dead + tca_dead_c(cc)
            rem         = rem - take
          END DO
          ct_axes_live = 0.0_wp
          DO cc = 0, CT_MAXTIL
            ct_axes_live = ct_axes_live + ct_dens(cc)
          END DO
          ct_bai_pot = 0.0_wp ; ct_bai = 0.0_wp
          DO cc = 0, CT_MAXTIL
            ct_bai_pot = ct_bai_pot + ct_dens(cc) * ct_area_axis(cc)      * 1.0e-4_wp
            ct_bai     = ct_bai     + ct_dens(cc) * ct_area_axis_real(cc) * 1.0e-4_wp
          END DO
          ! the fertile population is final (pre-sterility survivors).
          ct_ears           = ct_axes_live
          tca4_endpoint_day = .TRUE.
          tca_frozen        = .TRUE.
        ELSE IF (clr_tau_Cd >= ct_tau_sm) THEN
          ! The population entering the window, captured PER CLASS on the first day the
          ! class has any density.
          DO cc = 0, CT_MAXTIL
            IF (tca_dens_pre(cc) <= 0.0_wp .AND. ct_dens(cc) > 0.0_wp) &
              tca_dens_pre(cc) = ct_dens(cc)
          END DO
          ! The normalized mortality clock of Abichou (2016) Ch.1 Table 4, read here
          ! (below 1 inside the window by the branch condition, 0 on window day 1).
          tca_taun = (clr_tau_Cd - ct_tau_sm) / (ct_tau_em - ct_tau_sm)
          ! Per-cohort development state and leaf area index.
          tca_ltot = 0.0_wp
          DO cc = 0, CT_MAXTIL
            tca_hs(cc) = ct_hs_axis(cc, nl_ms)
            tca_L(cc)  = ct_dens(cc) * ct_area_axis_real(cc) * 1.0e-4_wp
            tca_ltot   = tca_ltot + tca_L(cc)
          END DO
          IF (tca_ltot <= 0.0_wp) THEN
            tca_days_skip = tca_days_skip + 1
          ELSE
            tca_days = tca_days + 1
            ! Vertical order: Haun stage DESCENDING, ties by cohort index ASCENDING, a
            ! stable insertion sort from the index order.
            DO ii = 0, CT_MAXTIL
              ord(ii) = ii
            END DO
            DO ii = 1, CT_MAXTIL
              jj = ii
              DO WHILE (jj >= 1)
                IF (tca_hs(ord(jj)) <= tca_hs(ord(jj - 1))) EXIT
                cc = ord(jj) ; ord(jj) = ord(jj - 1) ; ord(jj - 1) = cc
                jj = jj - 1
              END DO
            END DO
            ! Exact layer absorption, top down. F_tot = SUM f_c telescopes to
            ! 1 - exp(-k L_tot), and the un-intercepted fraction is renormalized.
            lab = 0.0_wp ; tca_ftot = 0.0_wp
            DO ii = 0, CT_MAXTIL
              cc = ord(ii)
              f_abs(cc) = EXP(-ct_kL * lab) * (1.0_wp - EXP(-ct_kL * tca_L(cc)))
              lab       = lab + tca_L(cc)
              tca_ftot  = tca_ftot + f_abs(cc)
            END DO
            tca_sum_s = 0.0_wp
            DO cc = 0, CT_MAXTIL
              tca_share(cc) = f_abs(cc) / tca_ftot
              tca_sum_s     = tca_sum_s + tca_share(cc)
            END DO
            !supply, demand and the per-axis ledgers.
            tca_ident = 0.0_wp
            DO cc = 0, CT_MAXTIL
              tca_sax(cc) = 0.0_wp ; tca_dax(cc) = 0.0_wp
              IF (ct_dens(cc) <= 0.0_wp) CYCLE
              sup_ax = tca_share(cc) * tca_gsup_prev / ct_dens(cc)
              da_ax  = MAX(0.0_wp, ct_area_axis_real(cc) - tca_area_prev(cc))
              dem_ax = da_ax * 1.0e-4_wp * csl_LMA * g_leaf
              tca_S(cc)   = tca_S(cc) + sup_ax
              tca_D(cc)   = tca_D(cc) + dem_ax
              tca_sax(cc) = sup_ax
              tca_dax(cc) = dem_ax
              tca_ident   = tca_ident + tca_share(cc) * tca_gsup_prev
            END DO
            tca_ident = tca_ident - tca_gsup_prev
            DO cc = 0, CT_MAXTIL
              tca_phi(cc) = 1.0_wp
              IF (tca_D(cc) > 0.0_wp) tca_phi(cc) = MIN(1.0_wp, tca_S(cc) / tca_D(cc))
            END DO

            ! Driver record 123: the nitrogen term of the affordability test. The nitrogen
            ! demand is the cohort's own lamina increment, dA times the lamina mass per area,
            ! priced at the construction concentration of driver record 121.
            tna_step_today = .FALSE.
            IF (tiller_n_afford >= 1) THEN
              tna_step_today = .TRUE.
              tna_days_step  = tna_days_step + 1
              tna_pL_used    = MIN(1.0_wp, MAX(0.0_wp, npg_p_L))
              BLOCK
                LOGICAL :: bind
                bind = .FALSE.
                DO cc = 0, CT_MAXTIL
                  tna_dN_day(cc) = 0.0_wp
                  tna_phiC(cc)   = tca_phi(cc)
                  IF (ct_dens(cc) <= 0.0_wp) CYCLE
                    IF (g_leaf > 0.0_wp) tna_dN_day(cc) = tca_dax(cc) / g_leaf * npg_c_L
                    tna_DN(cc) = tna_DN(cc) + tna_dN_day(cc)
                    tna_SN(cc) = tna_SN(cc) + tna_pL_used * tna_dN_day(cc)
                  tna_phiN(cc) = 1.0_wp
                  IF (tna_DN(cc) > 0.0_wp) tna_phiN(cc) = MIN(1.0_wp, tna_SN(cc) / tna_DN(cc))
                    IF (tna_phiN(cc) < tca_phi(cc)) THEN
                      bind = .TRUE.
                      tca_phi(cc) = tna_phiN(cc)
                    END IF
                END DO
                IF (bind) tna_days_bind = tna_days_bind + 1
              END BLOCK
            END IF

            ! The running estimate of the eventual deaths, the fitted shape of Abichou
            ! (2016) Ch.1 Eqn T4.4, and the unchanged order: youngest first (Ch.1 p.29),
            ! the main stem never.
            tca4_U = 0.0_wp
            DO cc = 1, CT_MAXTIL
              tca4_U = tca4_U + tca_dens_pre(cc) * (1.0_wp - tca_phi(cc))
            END DO
            tca4_F      = 1.0_wp - EXP(-tca4_theta * tca_taun / (tca4_delta - tca_taun))
            tca4_target = tca4_F * tca4_U

            ! Driver record 134: by default the target is a monotone running maximum; when
            ! the record is on it tracks the current target directly.
            IF (tca4_target > tca_unaff) THEN
              tca_unaff = tca4_target
            END IF
            rem = tca_unaff ; tca_dead = 0.0_wp
            DO cc = CT_MAXTIL, 1, -1
              take = MIN(tca_dens_pre(cc), rem)
              IF (take > tca_dead_c(cc)) THEN
                tca_dead_c(cc) = take
              END IF
              ct_dens(cc) = tca_dens_pre(cc) - tca_dead_c(cc)
              tca_dead    = tca_dead + tca_dead_c(cc)
              rem         = rem - take
            END DO
            ! the canopy the crop carries is the canopy that survived.
            ct_axes_live = 0.0_wp
            DO cc = 0, CT_MAXTIL
              ct_axes_live = ct_axes_live + ct_dens(cc)
            END DO
            ct_bai_pot = 0.0_wp ; ct_bai = 0.0_wp
            DO cc = 0, CT_MAXTIL
              ct_bai_pot = ct_bai_pot + ct_dens(cc) * ct_area_axis(cc)      * 1.0e-4_wp
              ct_bai     = ct_bai     + ct_dens(cc) * ct_area_axis_real(cc) * 1.0e-4_wp
            END DO
          END IF
        END IF
      END BLOCK
    END IF
    ! Yesterday's potential per-axis blade area, so that tomorrow's demand is today
    ! increment. Written on EVERY call when driver record 99 is on, so the lag is
    ! the same on the first day of the window.
    IF (tiller_survival_mech >= 1) THEN
      DO c = 0, CT_MAXTIL
        tsv_apot_prev(c) = ct_area_axis_real(c)
        IF (tiller_survival_mech == 4) tca_area_prev(c) = ct_area_axis_real(c)
      END DO
    END IF

    ! Cross-check against Bos (1999): the site filling this population implies
    ! F_s = ln(axes per plant) / HS_MS. Diagnostic; nothing reads it.
    IF (clr_HS > 0.0_wp .AND. dens_plant > 0.0_wp .AND. ct_axes_born > 0.0_wp) THEN
      ct_fs_implied = LOG(ct_axes_born / dens_plant) / clr_HS
    END IF
  END SUBROUTINE update_canopy_tillering

  ! Final leaf number of primary tiller c on a main stem carrying nl_ms leaves.
  ! Abichou (2016) Ch.1 Table 3 Eqn T3.5, Nl_Ti = alpha Nl_MS - beta Nc_Ti.
  REAL(wp) FUNCTION ct_nl_tiller(c, nl_ms)
    INTEGER, INTENT(IN) :: c, nl_ms
    ct_nl_tiller = ct_nl_alpha * REAL(nl_ms, wp) - ct_nl_beta * REAL(c + ct_nc_offset, wp)
    IF (ct_nl_tiller < 1.0_wp) ct_nl_tiller = 1.0_wp
  END FUNCTION ct_nl_tiller

  ! Final blade length of tiller c at that tiller's own phytomer rank x. Abichou
  ! (2016) Ch.2 Table 2: below the shifted rank of the first elongated internode
  ! the profile follows the shifted main-stem profile.
  REAL(wp) FUNCTION ct_L_tiller(x, c, nl_ms)
    REAL(wp), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: c, nl_ms
    REAL(wp) :: r, xbrk, Lr
    r = REAL(nl_ms - ct_internodes_elong, wp)      ! Abichou doc p.66: "generally 4 to 5 elongated internodes"
    IF (r < 2.0_wp) r = 2.0_wp
    Lr   = cba_L_at(r, nl_ms)
    xbrk = r - cba_sTi(c)
    IF (xbrk > 1.0_wp .AND. x < xbrk) THEN
      ct_L_tiller = (Lr - cba_lb1) / (xbrk - 1.0_wp) * (x - 1.0_wp) + cba_lb1
    ELSE
      ct_L_tiller = cba_L_at(x - cba_sTi(c), nl_ms)
    END IF
  END FUNCTION ct_L_tiller

  ! Maximum blade width of tiller c at its own phytomer rank x. Abichou (2016)
  ! Ch.2 Table 2 Eqn T2.10, a LINEAR profile across the tiller's own ranks.
  REAL(wp) FUNCTION ct_W_tiller(x, c, nl_ms)
    REAL(wp), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: c, nl_ms
    REAL(wp) :: w1, wflg, nlt
    INTEGER  :: ir1
    ir1  = MIN(MAX(nl_ms, 1), SIZE(cba_w_irank))   ! inversed rank of main-stem leaf 1
    w1   = cba_w_irank(ir1)
    wflg = cba_w_irank(1)
    nlt  = ct_nl_tiller(c, nl_ms)
    IF (nlt > 1.0_wp) THEN
      ct_W_tiller = cba_kw * (wflg - w1) / (nlt - 1.0_wp) * (x - 1.0_wp) + cba_kw * w1
    ELSE
      ct_W_tiller = cba_kw * w1
    END IF
  END FUNCTION ct_W_tiller

  ! Potential, realized and sterile-variant blade area of ONE axis of tiller
  ! class c (cm2 per axis). The realized expansion uses that axis's OWN Haun-stage
  ! ramp, Abichou (2016) Ch.1 Eqn T3.3, linear from zero at its emergence.
  SUBROUTINE ct_tiller_area(c, nl_ms, a_pot, a_real, a_st)
    INTEGER,  INTENT(IN)  :: c, nl_ms
    REAL(wp), INTENT(OUT) :: a_pot, a_real, a_st
    INTEGER  :: x, nt
    REAL(wp) :: nlt, p1, r1, s1

    a_pot = 0.0_wp ; a_real = 0.0_wp ; a_st = 0.0_wp
    IF (.NOT. ct_appeared(c)) RETURN
    nlt = ct_nl_tiller(c, nl_ms)
    nt  = MAX(1, NINT(nlt))
    DO x = 1, nt
      CALL ct_rank_area(c, x, nt, nl_ms, p1, r1, s1)
      a_pot  = a_pot  + p1
      a_real = a_real + r1
      a_st   = a_st   + s1
    END DO
  END SUBROUTINE ct_tiller_area

  ! Potential, realized and sterile-variant area of ONE rank on ONE tiller class
  ! (cm2). The single implementation of tiller rank geometry: Abichou (2016) Ch.2
  ! Eqns T2.7a and T2.7b (length), T2.10 (width) and T2.13 and T2.16 (sterile),
  ! on that axis's OWN Haun ramp (Ch.1 Eqn T3.3).
  SUBROUTINE ct_rank_area(c, x, nt, nl_ms, a_pot, a_real, a_st)
    INTEGER,  INTENT(IN)  :: c, x, nt, nl_ms
    REAL(wp), INTENT(OUT) :: a_pot, a_real, a_st
    REAL(wp) :: nlt, lx, wx, ft, a_rate, t_lig, t_app, frac
    a_pot = 0.0_wp ; a_real = 0.0_wp ; a_st = 0.0_wp
    IF (.NOT. ct_appeared(c)) RETURN
    nlt = ct_nl_tiller(c, nl_ms)
    IF (ct_tau_em > ct_tau_birth(c)) THEN
      a_rate = nlt / (ct_tau_em - ct_tau_birth(c))     ! a_Ti of Eqn T3.3, leaves per degCd
    ELSE
      a_rate = 0.0_wp
    END IF
    lx = ct_L_tiller(REAL(x, wp), c, nl_ms)
    wx = ct_W_tiller(REAL(x, wp), c, nl_ms)
    ft = cba_form_factor(x, nt)
    ! Driver record 144: a carbon-coupled bonus on the potential ceiling only.
    ! Inactive alternative, see Part 1 Table 6.
    a_pot = lx * wx * ft * cbc_f(c, x)
    a_st  = (cba_kl_st * lx) * (cba_kw_st * wx) * ft
    IF (a_rate > 0.0_wp) THEN
      t_lig = ct_tau_birth(c) + REAL(x, wp) / a_rate
      t_app = MAX(0.0_wp, t_lig - clr_LED_Cd)
      frac  = (clr_tau_Cd - t_app) / clr_LED_Cd
      IF (frac < 0.0_wp) frac = 0.0_wp
      IF (frac > 1.0_wp) frac = 1.0_wp
      a_real = a_pot * frac
    END IF
  END SUBROUTINE ct_rank_area

  ! Haun stage of axis class c, capped at that axis's own final leaf number. Class
  ! 0 is the main stem; classes 1 and above follow Abichou (2016) Ch.1 Eqn T3.3.
  REAL(wp) FUNCTION ct_hs_axis(c, nl_ms)
    INTEGER, INTENT(IN) :: c, nl_ms
    REAL(wp) :: nlt, a_rate
    IF (c == 0) THEN
      ct_hs_axis = MIN(clr_HS, REAL(nl_ms, wp))
      RETURN
    END IF
    ct_hs_axis = 0.0_wp
    IF (.NOT. ct_appeared(c)) RETURN
    nlt = ct_nl_tiller(c, nl_ms)
    IF (ct_tau_em <= ct_tau_birth(c)) RETURN
    a_rate = nlt / (ct_tau_em - ct_tau_birth(c))
    ct_hs_axis = MIN(MAX(a_rate * (clr_tau_Cd - ct_tau_birth(c)), 0.0_wp), nlt)
  END FUNCTION ct_hs_axis

  ! Driver record 83: blade EXPANSION under root-zone soil water. Acker, Zhang,
  ! Wang, Yu & Zhou (2024) Ann Bot 134:787, doi:10.1093/aob/mcae114, Fig. 4B,
  ! wheat panel: leaf elongation rate is linear in soil relative water content
  ! over 0 to 1, with no significant CO2 effect.
  SUBROUTINE update_leaf_expansion_water()
    INTEGER  :: j, c, x, nt, nl_ms, n_top
    REAL(wp) :: wsum, srwc_sum, se_sum, th_sum, fc_sum, fr, fclk, dclk, a_p, a_r, a_s

    IF (.NOT. clr_active) RETURN

    ! --- root-fraction-weighted SRWC (Acker Eq. 1 basis: water / field capacity) --
    wsum = 0.0_wp ; srwc_sum = 0.0_wp ; se_sum = 0.0_wp ; th_sum = 0.0_wp ; fc_sum = 0.0_wp
    IF (ALLOCATED(hydro_root_dist_diag) .AND. ALLOCATED(theta) .AND. ALLOCATED(theta_33)) THEN
      DO j = 1, MIN(SIZE(hydro_root_dist_diag), SIZE(theta), SIZE(theta_33))
        fr = MAX(0.0_wp, hydro_root_dist_diag(j))
        IF (fr <= 0.0_wp) CYCLE
        wsum   = wsum + fr
        th_sum = th_sum + fr * theta(j)
        fc_sum = fc_sum + fr * theta_33(j)
        IF (theta_33(j) > 0.0_wp) srwc_sum = srwc_sum + fr * MIN(1.0_wp, theta(j) / theta_33(j))
        IF (theta_0(j) - theta_r(j) > 0.0_wp) &
          se_sum = se_sum + fr * MIN(1.0_wp, MAX(0.0_wp, (theta(j) - theta_r(j)) / (theta_0(j) - theta_r(j))))
      END DO
    END IF
    IF (wsum > 0.0_wp) THEN
      lew_srwc     = srwc_sum / wsum
      lew_se_rw    = se_sum / wsum
      lew_theta_rw = th_sum / wsum
      lew_thetafc_rw = fc_sum / wsum
    ELSE
      ! No root fraction yet, i.e. before the first hydrology call or before roots
      ! exist: the leaf is not arrested.
      lew_srwc = 1.0_wp ; lew_se_rw = 1.0_wp ; lew_theta_rw = 0.0_wp ; lew_thetafc_rw = 0.0_wp
    END IF
    ! --- expansion factor: the printed wheat regressions normalized at SRWC = 1 --
    lew_fw_200 = (lew_a200 + lew_b200 * lew_srwc) / (lew_a200 + lew_b200)
    lew_fw_800 = (lew_a800 + lew_b800 * lew_srwc) / (lew_a800 + lew_b800)
    lew_fw     = (0.5_wp * (lew_a200 + lew_a800) + 0.5_wp * (lew_b200 + lew_b800) * lew_srwc) &
                 / (0.5_wp * (lew_a200 + lew_a800) + 0.5_wp * (lew_b200 + lew_b800))

    nl_ms = cfp_n_top_culm()
    lew_area_axis = 0.0_wp ; lew_bai = 0.0_wp ; lew_bai_clock = 0.0_wp
    ! --- main stem: the clock fraction per rank, incremented under f_w ------
    n_top = MIN(MIN(MAX(clr_n_app, clr_n_lig), nl_ms), CLR_MAXRANK)
    DO x = 1, n_top
      IF (cba_A_pot(x) > 0.0_wp) THEN
        fclk = MIN(1.0_wp, MAX(0.0_wp, cba_A(x) / cba_A_pot(x)))
      ELSE
        fclk = 0.0_wp
      END IF
      dclk = MAX(0.0_wp, fclk - lew_fclk(0, x))
      lew_f(0, x)    = MIN(1.0_wp, lew_f(0, x) + dclk * lew_fw)
      lew_fclk(0, x) = fclk
      lew_area_axis(0) = lew_area_axis(0) + lew_f(0, x) * cba_A_pot(x)
    END DO
    ! --- tillers: the clock fraction per rank on each cohort's own ramp -----
    IF (canopy_tillering == 1) THEN
      DO c = 1, CT_MAXTIL
        IF (.NOT. ct_appeared(c)) CYCLE
        nt = MAX(1, NINT(ct_nl_tiller(c, nl_ms)))
        DO x = 1, MIN(nt, CLR_MAXRANK)
          CALL ct_rank_area(c, x, nt, nl_ms, a_p, a_r, a_s)
          IF (a_p > 0.0_wp) THEN
            fclk = MIN(1.0_wp, MAX(0.0_wp, a_r / a_p))
          ELSE
            fclk = 0.0_wp
          END IF
          dclk = MAX(0.0_wp, fclk - lew_fclk(c, x))
          lew_f(c, x)    = MIN(1.0_wp, lew_f(c, x) + dclk * lew_fw)
          lew_fclk(c, x) = fclk
          lew_area_axis(c) = lew_area_axis(c) + lew_f(c, x) * a_p
        END DO
      END DO
      DO c = 0, CT_MAXTIL
        lew_bai       = lew_bai       + ct_dens(c) * lew_area_axis(c)     * 1.0e-4_wp
        lew_bai_clock = lew_bai_clock + ct_dens(c) * ct_area_axis_real(c) * 1.0e-4_wp
      END DO
    END IF
  END SUBROUTINE update_leaf_expansion_water

  ! Driver record 144: carbon-coupled leaf-area expansion. Called after
  ! update_leaf_expansion_water, which supplies today's water factor, and before
  ! the senescence step. Inactive alternative, see Part 1 Table 6.
  SUBROUTINE update_leaf_carbon_bonus()
    INTEGER  :: c, x, nt, nl_ms, n_top
    REAL(wp) :: fclk, dclk, inc, a_p, a_r, a_s, dens_sum, f_wsum

    cbc_live_today = .FALSE.
    cbc_f_mean_expanding = 1.0_wp
    cbc_area_bonus_day = 0.0_wp
    IF (carbon_coupled_leaf_expansion < 1) RETURN
    IF (canopy_blade_area /= 1 .OR. canopy_tillering /= 1) RETURN
    IF (.NOT. clr_active) RETURN
    cbc_live_today = .TRUE.

    nl_ms = cfp_n_top_culm()
    n_top = MIN(MIN(MAX(clr_n_app, clr_n_lig), nl_ms), CLR_MAXRANK)
    dens_sum = 0.0_wp ; f_wsum = 0.0_wp

    ! --- main stem (axis 0): the clock fraction per rank ------------- !
    DO x = 1, n_top
      a_p = cba_A_pot(x)
      fclk = 0.0_wp
      IF (a_p > 0.0_wp) fclk = MIN(1.0_wp, MAX(0.0_wp, cba_A(x) / a_p))
      dclk = MAX(0.0_wp, fclk - cbc_fclk(0, x))
      cbc_fclk(0, x) = fclk
      IF (dclk > 0.0_wp .AND. lew_fw >= 1.0_wp .AND. cbc_surplus_day > 0.0_wp) THEN
        inc = dclk * cbc_surplus_day
        cbc_f(0, x) = MIN(cbc_fmax, MAX(1.0_wp, cbc_f(0, x) + inc))
      END IF
      IF (cbc_f(0, x) > cbc_f_max_reached) cbc_f_max_reached = cbc_f(0, x)
      IF (fclk > 0.0_wp .AND. fclk < 1.0_wp) THEN
        dens_sum = dens_sum + ct_dens(0)
        f_wsum   = f_wsum   + ct_dens(0) * cbc_f(0, x)
      END IF
      IF (cbc_f(0, x) > 1.0_wp) cbc_area_bonus_day = cbc_area_bonus_day + &
        ct_dens(0) * cba_A(x) * (cbc_f(0, x) - 1.0_wp) / cbc_f(0, x) * 1.0e-4_wp
    END DO

    ! --- tillers: the clock fraction per rank on each cohort ramp -- !
    IF (canopy_tillering == 1) THEN
      DO c = 1, CT_MAXTIL
        IF (.NOT. ct_appeared(c)) CYCLE
        nt = MAX(1, NINT(ct_nl_tiller(c, nl_ms)))
        DO x = 1, MIN(nt, CLR_MAXRANK)
          CALL ct_rank_area(c, x, nt, nl_ms, a_p, a_r, a_s)
          fclk = 0.0_wp
          IF (a_p > 0.0_wp) fclk = MIN(1.0_wp, MAX(0.0_wp, a_r / a_p))
          dclk = MAX(0.0_wp, fclk - cbc_fclk(c, x))
          cbc_fclk(c, x) = fclk
          IF (dclk > 0.0_wp .AND. lew_fw >= 1.0_wp .AND. cbc_surplus_day > 0.0_wp) THEN
            inc = dclk * cbc_surplus_day
            cbc_f(c, x) = MIN(cbc_fmax, MAX(1.0_wp, cbc_f(c, x) + inc))
          END IF
          IF (cbc_f(c, x) > cbc_f_max_reached) cbc_f_max_reached = cbc_f(c, x)
          IF (fclk > 0.0_wp .AND. fclk < 1.0_wp) THEN
            dens_sum = dens_sum + ct_dens(c)
            f_wsum   = f_wsum   + ct_dens(c) * cbc_f(c, x)
          END IF
          IF (cbc_f(c, x) > 1.0_wp) cbc_area_bonus_day = cbc_area_bonus_day + &
            ct_dens(c) * a_r * (cbc_f(c, x) - 1.0_wp) / cbc_f(c, x) * 1.0e-4_wp
        END DO
      END DO
    END IF

    IF (dens_sum > 0.0_wp) cbc_f_mean_expanding = f_wsum / dens_sum
    IF (cbc_area_bonus_day > cbc_area_bonus_cum) cbc_area_bonus_cum = cbc_area_bonus_day
    IF (cbc_surplus_day > 0.0_wp .AND. lew_fw >= 1.0_wp) cbc_days_surplus = cbc_days_surplus + 1
    IF (cls_gai > cbc_laimax) cbc_laimax = cls_gai
  END SUBROUTINE update_leaf_carbon_bonus

  ! Blade senescence and the standing-dead pool (driver record 67). It writes
  ! only the per-rank areas and the dead pool; no carbon, nitrogen or water state
  ! is written here. The law is Abichou (2016) Ch.2 Table 3 with its four green-
  ! leaf-number phases.
  SUBROUTINE update_canopy_leaf_senesc()
    INTEGER  :: c, x, nt, nl_ms
    REAL(wp) :: nl_axis, hs_c, gl, sl, t1, t2, t3, a_p, a_r, a_s
    REAL(wp) :: dead_frac, g1_here, gl_ms, dd, area_axis_tot, dai_new
    REAL(wp) :: df_prev, died   ! rec 84: yesterday's dead fraction of the rank; today's retired green area index

    IF (.NOT. clr_active) RETURN
    cls_active = .TRUE.

    nl_ms = cfp_n_top_culm()

    ! ---- main stem first: the tiller equations need its Gl and its g2 ----- !
    hs_c = ct_hs_axis(0, nl_ms)
    t2 = clr_tau_lig_proj(MAX(1, nl_ms - cls_top_offset))   ! ligulation of leaf Nl-4 (Eqn T3.2)
    t3 = clr_tau_lig_proj(nl_ms)                            ! ligulation of the flag leaf (Eqn T3.3)
    IF (cls_tau1(0) < 0.0_wp .AND. clr_HS >= cls_g1) cls_tau1(0) = clr_tau_Cd
    IF (cls_hs_onset < 0.0_wp .AND. cls_tau1(0) >= 0.0_wp) cls_hs_onset = clr_HS
    t1 = cls_tau1(0)
    cls_tau2(0) = t2 ; cls_tau3(0) = t3
    CALL cls_mlaw(clr_tau_Cd, hs_c, t1, t2, t3, cls_g1, cls_g2, cls_g3, &
                  1.0_wp, gl, cls_phase(0))
    cls_gl(0) = MIN(gl, hs_c)
    cls_sl(0) = MAX(0.0_wp, hs_c - cls_gl(0))
    gl_ms = cls_gl(0)
    IF (cls_gl_at_flag < 0.0_wp .AND. clr_tau_Cd >= t3 .AND. t3 > 0.0_wp) cls_gl_at_flag = cls_gl(0)

    ! ---- tillers (Eqns T3.5-T3.8), each on its own Haun ramp -------------- !
    DO c = 1, CT_MAXTIL
      IF (.NOT. ct_appeared(c)) CYCLE
      nl_axis = ct_nl_tiller(c, nl_ms)
      hs_c    = ct_hs_axis(c, nl_ms)
      ! tau2 = ligulation of THIS axis's leaf Nl_Ti - 4, on its own ramp
      IF (ct_tau_em > ct_tau_birth(c) .AND. nl_axis > 0.0_wp) THEN
        t2 = ct_tau_birth(c) + MAX(0.0_wp, nl_axis - REAL(cls_top_offset, wp)) &
             * (ct_tau_em - ct_tau_birth(c)) / nl_axis
      ELSE
        t2 = ct_tau_em
      END IF
      t3 = ct_tau_em                       ! sigma_i = 0, the synchronous case

      ! Onset: for EARLY tillers senescence starts when their green leaf number
      ! reaches that of the main stem (Abichou 2016 Ch.2 p.64); a late tiller can
      ! begin senescing inside the third phase, which is that thesis's own Fig. 2d
      ! case.
      IF (cls_tau1(c) < 0.0_wp .AND. hs_c >= gl_ms .AND. clr_tau_Cd < t2) cls_tau1(c) = clr_tau_Cd
      t1 = cls_tau1(c)
      cls_tau2(c) = t2 ; cls_tau3(c) = t3
      g1_here = MERGE(hs_c, cls_g1, t1 < 0.0_wp)
      CALL cls_mlaw(clr_tau_Cd, hs_c, t1, t2, t3, g1_here, cls_g2, cls_g3, &
                    cls_delta, gl, cls_phase(c))
      cls_gl(c) = MIN(gl, hs_c)
      cls_sl(c) = MAX(0.0_wp, hs_c - cls_gl(c))
    END DO
    ! Driver record 122 mode 2: retire the green area nitrogen cannot support,
    ! oldest rank first, by raising the senesced-leaf number; returns at once when
    ! the mode is not 2.
    CALL lna_shed_by_rank(nl_ms)

    ! ---- retire blade area by rank, from the BASE upward ----------------- !
    ! Abichou (2016) Ch.2 Eqn 3: the senesced leaf number is the Haun stage minus
    ! the green leaf number, a decimal count, so the partial rank is retired pro
    ! rata.
    cls_gai = 0.0_wp ; cls_dai = 0.0_wp ; died = 0.0_wp ; cls_gai_flag = 0.0_wp
    DO c = 0, CT_MAXTIL
      cls_area_green(c) = 0.0_wp ; cls_area_dead(c) = 0.0_wp
      IF (.NOT. ct_appeared(c)) CYCLE
      IF (c == 0) THEN
        nt = MIN(MIN(MAX(clr_n_app, clr_n_lig), nl_ms), CLR_MAXRANK)
      ELSE
        nt = MAX(1, NINT(ct_nl_tiller(c, nl_ms)))
      END IF
      sl = cls_sl(c)
      DO x = 1, nt
        IF (c == 0) THEN
          a_r = cba_A(x)                   ! S2's own realized main-stem rank area
          a_p = cba_A_pot(x)
        ELSE
          CALL ct_rank_area(c, x, nt, nl_ms, a_p, a_r, a_s)
        END IF
        ! Driver record 83: the realized area is the water-ARRESTED accumulation and
        ! not the clock ramp (update_leaf_expansion_water is called before this step).
        IF (leaf_expansion_water == 1) a_r = lew_f(c, x) * a_p
        ! Driver record 122: the realized area is the water times NITROGEN arrested
        ! accumulation, taken here before the commit.
        IF (lai_n_arrest >= 1) THEN
          a_r = lna_fprov(c, x) * a_p
          IF (lna_dinc(c, x) > 0.0_wp .AND. sl - REAL(x - 1, wp) > 0.0_wp) lna_expdead_rank_days = lna_expdead_rank_days + 1
        END IF
        dead_frac = sl - REAL(x - 1, wp)
        IF (dead_frac < 0.0_wp) dead_frac = 0.0_wp
        IF (dead_frac > 1.0_wp) dead_frac = 1.0_wp
        cls_area_green(c) = cls_area_green(c) + a_r * (1.0_wp - dead_frac)
        cls_area_dead(c)  = cls_area_dead(c)  + a_r * dead_frac
        ! The TOP rank of this axis class is its flag lamina.
        IF (x == nt) cls_gai_flag = cls_gai_flag &
                                  + ct_dens(c) * a_r * (1.0_wp - dead_frac) * 1.0e-4_wp
        ! Driver record 84: the green area this rank lost TODAY is the rank area times
        ! the change in dead fraction, on the axes that existed yesterday. The senesced
        ! leaf number of Abichou (2016) Ch.2 Eqn 3 is the only thing that moves it.
        IF (lai_from_organs == 1) THEN
          df_prev = cls_sl_prev(c) - REAL(x - 1, wp)
          IF (df_prev < 0.0_wp) df_prev = 0.0_wp
          IF (df_prev > 1.0_wp) df_prev = 1.0_wp
          IF (dead_frac > df_prev) died = died + cls_dens_prev(c) * a_r * (dead_frac - df_prev) * 1.0e-4_wp
        END IF
      END DO
      cls_gai = cls_gai + ct_dens(c) * cls_area_green(c) * 1.0e-4_wp
      cls_dai = cls_dai + ct_dens(c) * cls_area_dead(c)  * 1.0e-4_wp
    END DO

    ! Axes culled by the mortality phase do not vanish from the canopy: their blade
    ! area stays standing until harvest. Retiring it here keeps the standing-dead
    ! pool complete.
    DO c = 1, CT_MAXTIL
      dd = cls_dens_prev(c) - ct_dens(c)
      IF (dd > 0.0_wp) THEN
        area_axis_tot = cls_area_green(c) + cls_area_dead(c)
        cls_dai_culled = cls_dai_culled + dd * area_axis_tot * 1.0e-4_wp
        ! Driver record 99 closure, first of two limbs. Accumulated OUTSIDE the
        ! record-84 branch below, from the same increment but in a separate statement,
        ! so that it records the green area the cull ACTUALLY removed from the canopy.
        IF (tiller_survival_mech >= 1) &
          tsv_cullgreen_seen = tsv_cullgreen_seen + dd * cls_area_green(c) * 1.0e-4_wp

        ! Under driver record 84 the GREEN blades of culled axes die with the axis, and
        ! the same area is charged once to the organ law closure.
        IF (lai_from_organs == 1) THEN
          died = died + dd * cls_area_green(c) * 1.0e-4_wp
          IF (tiller_survival_mech >= 1) &
            tsv_cullgreen_charged = tsv_cullgreen_charged + dd * cls_area_green(c) * 1.0e-4_wp
        END IF
      END IF
      cls_dens_prev(c) = ct_dens(c)
    END DO
    cls_dens_prev(0) = ct_dens(0)
    ! Driver record 84: publish today's retired green area and the fraction of the
    ! pre-retirement green area it represents; apply_leaf_senescence sheds the same
    ! fraction of the green lamina mass and of the leaf nitrogen.
    IF (lai_from_organs == 1) THEN
      lfo_gai_died_day = died
      lfo_fdied_day = 0.0_wp
      IF (cls_gai + died > 0.0_wp) lfo_fdied_day = MIN(1.0_wp, died / (cls_gai + died))
      cls_sl_prev(0:CT_MAXTIL) = cls_sl(0:CT_MAXTIL)
    END IF

    ! The standing-dead POOL is an accumulator and not a re-computation: area that
    ! has been retired stays retired, so the pool rises by the day POSITIVE
    ! increment only.
    IF (cls_dai_prev >= 0.0_wp) THEN
      dai_new = cls_dai - cls_dai_prev
      IF (dai_new > 0.0_wp) THEN
        cls_dai_pool = cls_dai_pool + dai_new
      ELSE IF (dai_new < 0.0_wp) THEN
        cls_mono_viol = cls_mono_viol + 1
        IF (-dai_new > cls_dai_dip_max) cls_dai_dip_max = -dai_new
      END IF
    ELSE
      cls_dai_pool = cls_dai_pool + cls_dai
    END IF
    cls_dai_prev = cls_dai
    IF (cls_dai_pool + cls_dai_culled > cls_dai_max) cls_dai_max = cls_dai_pool + cls_dai_culled

    ! Driver record 147 diagnostics: the standing green area AFTER the rank
    ! retirement loop of this same day, and the mean green-lamina specific leaf
    ! nitrogen on it.
  END SUBROUTINE update_canopy_leaf_senesc

  ! The green-leaf-number law of Abichou (2016) Ch.2 Table 3, written once and
  ! used for both the main stem (scale 1) and the tillers (scale delta, Ch.2
  ! Table 6).
  SUBROUTINE cls_mlaw(tau, hs, tau1, tau2, tau3, g1, g2, g3, scale, gl, phase)
    REAL(wp), INTENT(IN)  :: tau, hs, tau1, tau2, tau3, g1, g2, g3, scale
    REAL(wp), INTENT(OUT) :: gl
    INTEGER,  INTENT(OUT) :: phase
    REAL(wp) :: dt
    IF (tau3 > 0.0_wp .AND. tau >= tau3) THEN
      dt = tau - tau3                                             ! phi4, Eqn T3.4 / T3.8
      gl = scale * (cls_alpha * dt**3 + cls_beta * dt + g3)
      phase = 4
    ELSE IF (tau2 > 0.0_wp .AND. tau >= tau2 .AND. tau3 > tau2) THEN
      gl = (scale * g3 - g2) * (tau - tau2) / (tau3 - tau2) + g2  ! phi3, Eqn T3.3 / T3.7
      phase = 3
    ELSE IF (tau1 >= 0.0_wp .AND. tau > tau1 .AND. tau2 > tau1) THEN
      gl = (g2 - g1) * (tau - tau1) / (tau2 - tau1) + g1          ! phi2, Eqn T3.2 / T3.6
      phase = 2
    ELSE
      gl = hs                                                     ! phi1, Eqn T3.1 / T3.5
      phase = 1
    END IF
    IF (gl < 0.0_wp) gl = 0.0_wp
  END SUBROUTINE cls_mlaw

  ! Non-lamina standing area: sheaths, exposed stem and ear with peduncle (driver
  ! record 68). It writes only the areas; no carbon, nitrogen or water state is
  ! written here. The geometry is Abichou (2016) Ch.2 Table 2 and Fig. 3b and 3c.
  SUBROUTINE update_canopy_sheath_ear(istage_idx)
    INTEGER, INTENT(IN) :: istage_idx
    INTEGER  :: c, x, nt, nl_ms
    REAL(wp) :: nl_axis, sl, ls, ds, li, a_sh, frac_ext
    REAL(wp) :: h_node, collar, collar_prev, base_vis, vis_len, bare_len
    REAL(wp) :: g_n2, g_n1, g_n3, t_lig, t_app, a_rate, dens
    REAL(wp) :: sh_n1, sh_n3, sh_full, ep_frac, t_flag, win, sdi_new
    REAL(wp) :: green_blade, dead_blade, sh_g, sh_d, st_g, ep_g
    REAL(wp) :: sai_n1_sum, sai_n3_sum, sai_full_sum

    IF (.NOT. clr_active) RETURN
    cse_active = .TRUE.

    nl_ms = cfp_n_top_culm()
    win   = cse_ep_win * clr_phyllochron_Cd        ! ear+peduncle emergence window, degCd

    ! ---- main-stem sheath profile, for the diagnostic dump ---------------- !
    DO x = 1, MIN(nl_ms, CLR_MAXRANK)
      cse_ls_pot(x) = cse_Ls_at(REAL(x, wp), nl_ms)
      cse_ds_pot(x) = cse_Ds_ms(REAL(x, wp), nl_ms)
    END DO

    ! ---- per axis class: sheath, exposed stem, ear ------------------------ !
    cse_sai = 0.0_wp ; cse_sdi = 0.0_wp ; cse_eai = 0.0_wp ; cse_iai = 0.0_wp
    cse_eai_lo = 0.0_wp ; cse_eai_hi = 0.0_wp ; cse_eai_meas = 0.0_wp
    sai_n1_sum = 0.0_wp ; sai_n3_sum = 0.0_wp ; sai_full_sum = 0.0_wp
    DO c = 0, CT_MAXTIL
      cse_sh_green(c) = 0.0_wp ; cse_sh_dead(c) = 0.0_wp
      cse_ep_axis(c)  = 0.0_wp ; cse_st_axis(c) = 0.0_wp
      IF (.NOT. ct_appeared(c)) CYCLE
      dens = ct_dens(c)
      IF (c == 0) THEN
        nt      = MIN(MIN(MAX(clr_n_app, clr_n_lig), nl_ms), CLR_MAXRANK)
        nl_axis = REAL(nl_ms, wp)
        a_rate  = 0.0_wp
      ELSE
        nl_axis = ct_nl_tiller(c, nl_ms)
        nt      = MAX(1, NINT(nl_axis))
        a_rate  = 0.0_wp
        IF (ct_tau_em > ct_tau_birth(c)) a_rate = nl_axis / (ct_tau_em - ct_tau_birth(c))
      END IF
      sl = cls_sl(c)
      sh_n1 = 0.0_wp ; sh_n3 = 0.0_wp ; sh_full = 0.0_wp
      h_node = 0.0_wp ; collar_prev = 0.0_wp
      DO x = 1, nt
        ! --- final dimensions of THIS phytomer (Eqns T2.2/T2.5ab, T2.8/T2.11)
        IF (c == 0) THEN
          ls = cse_Ls_at(REAL(x, wp), nl_ms)
          ds = cse_Ds_ms(REAL(x, wp), nl_ms)
          li = cse_Li_at(REAL(x, wp), nl_ms)
        ELSE
          ls = cse_Ls_at(REAL(x, wp) - cba_sTi(c), nl_ms)                  ! T2.8
          ds = cse_ko_tiller(c) * cse_Ds_ms(REAL(x, wp) - cba_sTi(c), nl_ms) ! T2.11
          li = cse_Li_at(REAL(x, wp) - cba_sTi(c), nl_ms)                  ! T2.9
        END IF

        ! --- NESTING: sheaths form a whorl, so only the OUTERMOST surface at a given
        !     height intercepts anything. Sheath x runs from its node upward, so the
        !     nested length is subtracted.
        h_node   = h_node + li
        base_vis = MAX(h_node, collar_prev)
        collar   = h_node + ls
        vis_len  = MAX(0.0_wp, collar - base_vis)
        bare_len = MAX(0.0_wp, h_node - collar_prev)
        collar_prev = collar
        a_sh = vis_len * ds                              ! visible sheath band
        ! --- ligulation of the BLADE of this phytomer, from the axis ramp
        IF (c == 0) THEN
          t_lig = clr_tau_lig_proj(x)
        ELSE IF (a_rate > 0.0_wp) THEN
          t_lig = ct_tau_birth(c) + REAL(x, wp) / a_rate
        ELSE
          t_lig = -1.0_wp
        END IF
        ! --- sequential extension: the sheath grows AFTER its own blade, from that
        !     blade ligulation to tau_app + (0.4 + 2) phyllochrons.
        frac_ext = 0.0_wp
        IF (t_lig > 0.0_wp) THEN
          t_app    = t_lig - clr_LED_Cd
          frac_ext = cse_ext_frac(clr_tau_Cd, t_lig, t_app)
        END IF
        ! --- ADEL persistence: sheath n senesces with blade n + cse_persist_n
        g_n2 = cse_clamp01(REAL(x + cse_persist_n, wp) - sl)
        g_n1 = cse_clamp01(REAL(x + 1, wp) - sl)
        g_n3 = cse_clamp01(REAL(x + 3, wp) - sl)
        cse_sh_green(c) = cse_sh_green(c) + a_sh * frac_ext * g_n2
        cse_sh_dead(c)  = cse_sh_dead(c)  + a_sh * frac_ext * (1.0_wp - g_n2)
        sh_n1   = sh_n1   + a_sh * frac_ext * g_n1
        sh_n3   = sh_n3   + a_sh * frac_ext * g_n3
        sh_full = sh_full + a_sh * g_n2                  ! upper-bound arm: no extension lag

        ! --- the exposed (bare) internode, the gap the whorl leaves uncovered. The
        !     coordination rule of Abichou (2016) Annex 2 extends the internode
        !     sequentially AFTER the sheath, and internodes do not senesce there, so
        !     the exposed length stays green.
        cse_st_axis(c) = cse_st_axis(c) + bare_len * cse_di2 * frac_ext
      END DO
      ! --- ear + exposed peduncle, ramped from THIS axis's flag ligulation --
      IF (c == 0) THEN
        t_flag = clr_tau_lig_proj(nl_ms)
      ELSE
        t_flag = ct_tau_em                               ! sigma_i = 0, the synchronous case
      END IF
      ep_frac = 0.0_wp
      IF (t_flag > 0.0_wp .AND. win > 0.0_wp) ep_frac = cse_clamp01((clr_tau_Cd - t_flag) / win)
      cse_ep_axis(c) = cse_a_ep * ep_frac

      cse_sai    = cse_sai    + dens * cse_sh_green(c) * 1.0e-4_wp
      cse_sdi    = cse_sdi    + dens * cse_sh_dead(c)  * 1.0e-4_wp
      cse_iai    = cse_iai    + dens * cse_st_axis(c)  * 1.0e-4_wp
      cse_eai    = cse_eai    + dens * cse_ep_axis(c)  * 1.0e-4_wp
      cse_eai_lo = cse_eai_lo + dens * cse_a_ep_lo * ep_frac * 1.0e-4_wp
      cse_eai_hi = cse_eai_hi + dens * cse_a_ep_hi * ep_frac * 1.0e-4_wp
      cse_eai_meas = cse_eai_meas + dens * cse_a_ear_meas * ep_frac * 1.0e-4_wp
      sai_n1_sum   = sai_n1_sum   + dens * sh_n1   * 1.0e-4_wp
      sai_n3_sum   = sai_n3_sum   + dens * sh_n3   * 1.0e-4_wp
      sai_full_sum = sai_full_sum + dens * sh_full * 1.0e-4_wp
    END DO
    cse_sai_n1 = sai_n1_sum ; cse_sai_n3 = sai_n3_sum ; cse_sai_full = sai_full_sum

    ! --- standing-dead sheath as a POOL, monotone by construction --------- !
    ! The same construction as the blade pool and for the same reason: the
    ! instantaneous re-computation is kept beside it for diagnosis only.
    IF (cse_sdi_prev >= 0.0_wp) THEN
      sdi_new = cse_sdi - cse_sdi_prev
      IF (sdi_new > 0.0_wp) THEN
        cse_sdi_pool = cse_sdi_pool + sdi_new
      ELSE IF (sdi_new < 0.0_wp) THEN
        cse_mono_viol = cse_mono_viol + 1
        IF (-sdi_new > cse_sdi_dip_max) cse_sdi_dip_max = -sdi_new
      END IF
    ELSE
      cse_sdi_pool = cse_sdi_pool + cse_sdi
    END IF
    cse_sdi_prev = cse_sdi

    ! --- the two BASES, kept apart ---------------------------------------- !
    green_blade = cls_gai
    dead_blade  = cls_dai_pool + cls_dai_culled
    sh_g = cse_sai ; sh_d = cse_sdi_pool ; st_g = cse_iai ; ep_g = cse_eai
    cse_nonlam = sh_g + st_g + ep_g
    cse_ratio  = 0.0_wp
    IF (green_blade > 1.0e-9_wp) cse_ratio = cse_nonlam / green_blade
    ! The destructive basis of Abichou (2016): GREEN organs only, projected areas,
    ! exactly the scans of leaves and stems of that thesis (Ch.3 p.95 and p.101).
    cse_paig = green_blade + cse_nonlam
    ! The gap-fraction basis of the LI-2200C: every STANDING element, with the
    ! cylindrical organs on the half-total-surface convention a gap-fraction
    ! inversion returns (pi/2 L D); flat blades enter unchanged.
    cse_pai_li = green_blade + dead_blade &
                 + cse_half_cyl * (sh_g + sh_d + st_g + ep_g)
    IF (cse_pai_li > cse_pai_li_max) cse_pai_li_max = cse_pai_li
    cse_paig_axis = 0.0_wp
    IF (ct_axes_live > 1.0e-9_wp) cse_paig_axis = cse_paig / ct_axes_live * 1.0e4_wp
    ! One snapshot at the anchored anthesis stage. The anchor is the model's own
    ! phenology and not a date derived from observations.
    IF (cse_paig_flower <= 0.0_wp .AND. istage_idx >= 6) THEN
      cse_paig_flower  = cse_paig
      cse_axis_flower  = ct_axes_live
      cse_ratio_flower = cse_ratio
    END IF
  END SUBROUTINE update_canopy_sheath_ear
  ! ----------------------------------------------------------------------------------- !
  ! Final SHEATH length at a possibly NON-INTEGER absolute rank, cm, by linear
  ! interpolation on the inversed-rank median profile of Abichou (2016) Ch.2
  ! Fig. 3b. Decimal ranks are required because Eqn T2.8 gives the tiller profile
  ! as a shift of the main-stem one.
  ! ----------------------------------------------------------------------------------- !
  REAL(wp) FUNCTION cse_Ls_at(xabs, n_max)
    REAL(wp), INTENT(IN) :: xabs
    INTEGER,  INTENT(IN) :: n_max
    cse_Ls_at = cse_prof_at(cse_ls_irank, xabs, n_max)
  END FUNCTION cse_Ls_at

  ! Final INTERNODE length at a possibly non-integer absolute rank, cm; Abichou
  ! (2016) Ch.2 Fig. 3c p.67, peduncle NOT included (glossary p.81).
  REAL(wp) FUNCTION cse_Li_at(xabs, n_max)
    REAL(wp), INTENT(IN) :: xabs
    INTEGER,  INTENT(IN) :: n_max
    cse_Li_at = cse_prof_at(cse_li_irank, xabs, n_max)
  END FUNCTION cse_Li_at

  ! Shared inversed-rank interpolator for the two sheath and internode profiles,
  ! written once.
  REAL(wp) FUNCTION cse_prof_at(prof, xabs, n_max)
    REAL(wp), INTENT(IN) :: prof(:), xabs
    INTEGER,  INTENT(IN) :: n_max
    REAL(wp) :: ir, w
    INTEGER  :: i0, i1
    ir = REAL(n_max, wp) - xabs + 1.0_wp          ! inversed rank, decimal
    IF (ir <= 1.0_wp) THEN
      cse_prof_at = prof(1)
    ELSE IF (ir >= REAL(SIZE(prof), wp)) THEN
      cse_prof_at = prof(SIZE(prof))
    ELSE
      i0 = INT(ir) ; i1 = i0 + 1 ; w = ir - REAL(i0, wp)
      cse_prof_at = prof(i0) * (1.0_wp - w) + prof(i1) * w
    END IF
  END FUNCTION cse_prof_at

  ! Final SHEATH diameter of the main stem at absolute rank x, cm; Abichou (2016)
  ! Ch.2 Table 2 Eqns T2.5a and T2.5b, a linear ramp to a plateau between the two
  ! measured end diameters.
  REAL(wp) FUNCTION cse_Ds_ms(xabs, n_max)
    REAL(wp), INTENT(IN) :: xabs
    INTEGER,  INTENT(IN) :: n_max
    REAL(wp) :: r, xx
    r  = REAL(n_max - ct_internodes_elong, wp)
    IF (r < 2.0_wp) r = 2.0_wp
    xx = xabs
    IF (xx < 1.0_wp) xx = 1.0_wp
    IF (xx <= r) THEN
      cse_Ds_ms = (cse_ds2 - cse_ds1) / (r - 1.0_wp) * (xx - 1.0_wp) + cse_ds1   ! T2.5a
    ELSE
      cse_Ds_ms = cse_ds2                                                        ! T2.5b
    END IF
  END FUNCTION cse_Ds_ms
  ! ----------------------------------------------------------------------------------- !
  ! The tiller sheath and internode diameter reduction factor; Abichou (2016)
  ! Ch.2 Eqn T2.11 with the glossary definition: the last emerged tiller takes
  ! the lower constant, every other tiller the lower-bounded linear form.
  ! ----------------------------------------------------------------------------------- !
  REAL(wp) FUNCTION cse_ko_tiller(c)
    INTEGER, INTENT(IN) :: c
    INTEGER  :: i, n_til, last
    REAL(wp) :: p
    n_til = 0 ; last = 0
    DO i = 1, CT_MAXTIL
      IF (ct_appeared(i)) THEN
        n_til = n_til + 1 ; last = i
      END IF
    END DO
    IF (n_til <= 0) THEN
      cse_ko_tiller = 1.0_wp ; RETURN
    END IF
    IF (c == last) THEN
      cse_ko_tiller = cse_rho3 ; RETURN
    END IF
    p = REAL(c, wp) / REAL(n_til, wp)
    cse_ko_tiller = MAX(cse_rho1 * p + cse_rho2, cse_rho3)
  END FUNCTION cse_ko_tiller

  ! Extension fraction of one sheath. The coordination rule of Abichou (2016)
  ! Annex 2 makes the sheath extend AFTER its own blade, so the window opens at
  ! the blade ligulation.
  REAL(wp) FUNCTION cse_ext_frac(tau, t_lig, t_app)
    REAL(wp), INTENT(IN) :: tau, t_lig, t_app
    REAL(wp) :: t_end, dur
    t_end = t_app + (cse_adel_start + cse_adel_dur) * clr_phyllochron_Cd
    dur   = t_end - t_lig
    IF (dur <= 0.0_wp) THEN
      cse_ext_frac = MERGE(1.0_wp, 0.0_wp, tau >= t_lig)
    ELSE
      cse_ext_frac = cse_clamp01((tau - t_lig) / dur)
    END IF
  END FUNCTION cse_ext_frac

  REAL(wp) FUNCTION cse_clamp01(v)
    REAL(wp), INTENT(IN) :: v
    cse_clamp01 = MIN(1.0_wp, MAX(0.0_wp, v))
  END FUNCTION cse_clamp01

  ! The sink split: today's NEW LAMINA AREA is converted into a carbon demand, what
  ! the supply allows is debited, and the rest is carried as a non-negative
  ! backlog.
  SUBROUTINE sink_lamina_split(dm_gain, growth_supply, g_total, &
                               f_leaf, f_stem, f_root, f_storage)
    REAL(wp), INTENT(INOUT) :: dm_gain, g_total
    REAL(wp), INTENT(IN)    :: growth_supply
    REAL(wp), INTENT(INOUT) :: f_leaf, f_stem, f_root, f_storage
    REAL(wp) :: dA, demand, alloc, rem, s_non, fs, fr, fo, err

    ! --- today's NEW lamina area, from the monotone CUMULATIVE built area ---
    dA = csl_area_built - csl_area_prev
    IF (dA < 0.0_wp) dA = 0.0_wp        ! cannot un-build; the pool is monotone
    csl_area_prev = csl_area_built
    csl_dA_day    = dA

    ! --- the carbon that area costs, plus whatever is still owed ----------- !
    demand = dA * csl_LMA + csl_backlog
    csl_demand_day = demand

    ! --- debit what the supply allows; the rest stays owed ----------------- !
    alloc = MIN(demand, MAX(0.0_wp, dm_gain))
    IF (alloc < 0.0_wp) alloc = 0.0_wp
    csl_alloc_day = alloc
    csl_backlog   = MAX(0.0_wp, demand - alloc)          ! HARD non-negative
    IF (csl_backlog > csl_backlog_max) csl_backlog_max = csl_backlog
    csl_limit_active = (csl_backlog > 0.0_wp)

    ! --- fractions that reproduce the split, renormalized over non-leaf ---- !
    IF (dm_gain > 0.0_wp) THEN
      f_leaf = alloc / dm_gain
    ELSE
      f_leaf = 0.0_wp
    END IF
    IF (f_leaf > 1.0_wp) f_leaf = 1.0_wp
    rem   = 1.0_wp - f_leaf
    fs = MAX(0.0_wp, f_stem) ; fr = MAX(0.0_wp, f_root) ; fo = MAX(0.0_wp, f_storage)
    s_non = fs + fr + fo
    IF (s_non > 0.0_wp) THEN
      f_stem = rem * fs / s_non ; f_root = rem * fr / s_non ; f_storage = rem * fo / s_non
    ELSE
      f_stem = rem ; f_root = 0.0_wp ; f_storage = 0.0_wp
    END IF
  ! ----------------------------------------------------------------------------------- !
    ! Driver record 89 modes 1 and 2: the renormalization of driver record 69 would
    ! dilute the measured partitioning coefficient of Fischer (1985), so the spike
    ! share is held.
  ! ----------------------------------------------------------------------------------- !
    IF (grain_number_mech == 1 .OR. grain_number_mech == 2) &
      CALL gnx_hold_spike_share(f_leaf, f_stem, f_root, f_storage, fr, fo, .FALSE.)
    err = ABS(1.0_wp - (f_leaf + f_stem + f_root + f_storage))
    IF (err > csl_fsum_err) csl_fsum_err = err           ! INVARIANT 6d

    ! --- ONE corrective iteration on the glucose requirement -------------- !
    ! Eq. [7.6] p.148 of Teh (2006) weights the glucose requirement by the
    ! partition, and the sink has changed the partition, so it is re-solved once.
    g_total = f_leaf * g_leaf + f_stem * g_stem + f_root * g_root + f_storage * g_storage
    IF (g_total > 0.0_wp) THEN
      dm_gain = growth_supply / g_total
      alloc = MIN(demand, MAX(0.0_wp, dm_gain))
      csl_alloc_day = alloc
      csl_backlog   = MAX(0.0_wp, demand - alloc)
      IF (csl_backlog > csl_backlog_max) csl_backlog_max = csl_backlog
      csl_limit_active = (csl_backlog > 0.0_wp)
      IF (dm_gain > 0.0_wp) THEN
        f_leaf = alloc / dm_gain
      ELSE
        f_leaf = 0.0_wp
      END IF
      IF (f_leaf > 1.0_wp) f_leaf = 1.0_wp
      rem = 1.0_wp - f_leaf
      IF (s_non > 0.0_wp) THEN
        f_stem = rem * fs / s_non ; f_root = rem * fr / s_non ; f_storage = rem * fo / s_non
      ELSE
        f_stem = rem ; f_root = 0.0_wp ; f_storage = 0.0_wp
      END IF
      IF (grain_number_mech == 1 .OR. grain_number_mech == 2) &
        CALL gnx_hold_spike_share(f_leaf, f_stem, f_root, f_storage, fr, fo, .TRUE.)
      err = ABS(1.0_wp - (f_leaf + f_stem + f_root + f_storage))
      IF (err > csl_fsum_err) csl_fsum_err = err
    END IF
    IF (csl_limit_active) csl_limited_days = csl_limited_days + 1
  ! ----------------------------------------------------------------------------------- !
    ! Driver record 144: today's lamina-sink carbon SURPLUS fraction, gated by this
    ! record's own flag. Inactive alternative, see Part 1 Table 6.
  ! ----------------------------------------------------------------------------------- !
    IF (carbon_coupled_leaf_expansion >= 1) THEN
      IF (csl_backlog <= 0.0_wp .AND. demand > 1.0e-9_wp .AND. dm_gain > demand) THEN
        cbc_surplus_day = (dm_gain - demand) / demand
      ELSE
        cbc_surplus_day = 0.0_wp
      END IF
    END IF
  END SUBROUTINE sink_lamina_split

  ! The sink state the rest of the day reads: the realized lamina area, the
  ! carbon backlog and the specific leaf area that follows from them.
  SUBROUTINE update_canopy_sink_state()
    IF (.NOT. cls_active) RETURN
    csl_active = .TRUE.
    ! The SINK leaf area index is the GREEN LAMINA area index, the like-for-like
    ! partner of a destructive leaf-area measurement and the quantity the light
    ! model reads.
    csl_lai_sink = cls_gai
    ! The instrument-basis diagnostic, emitted daily so that an observation can be
    ! compared on its own basis. It is never fed back into any flux.
    csl_pai_diag = cse_pai_li
    ! CUMULATIVE built lamina area: green plus standing dead plus the area carried
    ! by axes the mortality phase culled. Monotone non-decreasing by construction,
    ! so senescence never re-charges carbon and never refunds it.
    csl_area_built = cls_gai + cls_dai_pool + cls_dai_culled
  END SUBROUTINE update_canopy_sink_state

  REAL(wp) FUNCTION cba_form_factor(n, n_max)
    INTEGER, INTENT(IN) :: n, n_max
    REAL(wp) :: t
    IF (n_max >= 1 .AND. n >= n_max) THEN
      cba_form_factor = cba_ft_flag
    ELSE IF (n <= 1) THEN
      cba_form_factor = cba_ft_r1
    ELSE IF (n <= cba_ft_n4) THEN
      t = REAL(n - 1, wp) / REAL(cba_ft_n4 - 1, wp)
      cba_form_factor = cba_ft_r1 + (cba_ft_r4 - cba_ft_r1) * t
    ELSE IF (n <= cba_ft_n8) THEN
      t = REAL(n - cba_ft_n4, wp) / REAL(cba_ft_n8 - cba_ft_n4, wp)
      cba_form_factor = cba_ft_r4 + (cba_ft_adult - cba_ft_r4) * t
    ELSE
      cba_form_factor = cba_ft_adult
    END IF
  END FUNCTION cba_form_factor

  ! Driver record 100: the POTENTIAL telescoped collar (pseudostem) height of the
  ! main stem, cm, at today's Haun stage. The internode and sheath profiles are
  ! Abichou (2016) PhD thesis (NNT tel-02316154) Ch.2 Fig. 3c and Fig. 3b, and
  ! the telescoping is checked against that thesis separately measured Fig. 3d
  ! collar height. Carbon-independent.
  REAL(wp) FUNCTION sdm_h_pot_now()
    INTEGER  :: nl_ms, x, nt
    REAL(wp) :: hs, h, hp, fr
    sdm_h_pot_now = 0.0_wp
    IF (.NOT. clr_active) RETURN
    nl_ms = cfp_n_top_culm()
    IF (nl_ms < 1) RETURN
    hs = MIN(MAX(0.0_wp, clr_HS), REAL(nl_ms, wp))
    nt = MIN(INT(hs), nl_ms)
    fr = hs - REAL(nt, wp)
    h  = 0.0_wp ; hp = 0.0_wp
    DO x = 1, nt
      h  = h + cse_Li_at(REAL(x, wp), nl_ms)
      hp = MAX(hp, h + cse_Ls_at(REAL(x, wp), nl_ms))
    END DO
    IF (fr > 0.0_wp .AND. nt < nl_ms) THEN
      hp = MAX(hp, h + fr * (cse_Li_at(REAL(nt + 1, wp), nl_ms) + &
                             cse_Ls_at(REAL(nt + 1, wp), nl_ms)))
    END IF
    sdm_h_pot_now = hp
  END FUNCTION sdm_h_pot_now

  ! Driver record 100: height increment to stem dry-mass demand, g DM m-2 d-1,
  ! through the inverse of the plant-height relation of Bakhshandeh et al. (2012)
  ! Cereal Res Commun 40(3):413-422, Eq. (6) p.419 and Fig. 4.
  REAL(wp) FUNCTION sdm_dw_from_dh(dh, w_stem)
    REAL(wp), INTENT(IN) :: dh, w_stem
    REAL(wp) :: x
    sdm_dw_from_dh = 0.0_wp
    IF (dh <= 0.0_wp .OR. sdm_dens <= 0.0_wp) RETURN
    x = w_stem / sdm_dens
    IF (x < stem_height_break_g_plant) THEN
      sdm_dw_from_dh = dh / stem_height_slope1_cm_g * sdm_dens
    ELSE
      sdm_dw_from_dh = dh / stem_height_slope2_cm_g * sdm_dens
    END IF
  END FUNCTION sdm_dw_from_dh

  ! Driver record 100 mode 3: has the telescoped collar column run out of ranks?
  ! True from the FLAG-LEAF LIGULE onwards, where the Haun stage has reached the
  ! final leaf number.
  LOGICAL FUNCTION sdm_collar_done()
    INTEGER :: nl_ms
    sdm_collar_done = .FALSE.
    IF (.NOT. clr_active) RETURN
    nl_ms = cfp_n_top_culm()
    IF (nl_ms < 1) RETURN
    sdm_collar_done = (clr_HS >= REAL(nl_ms, wp))
  END FUNCTION sdm_collar_done

  ! Driver record 100 mode 4: today's BOUNDED per-axis structural stem-mass demand,
  ! called once per day from main_MMCO2Crop.f90.
  SUBROUTINE sdm4_set_demand(dtt_day)
    REAL(wp), INTENT(IN) :: dtt_day

    IF (stem_demand_mech /= 4) RETURN
  END SUBROUTINE sdm4_set_demand

  ! The PROJECTED onset of rank x on axis class c (degree-Cd; a negative value
  ! means not yet definable). Rank nl_c + 1 is the peduncle.
  REAL(wp) FUNCTION sdm4_onset_proj(x, c, nl_c, a_rate)
    INTEGER,  INTENT(IN) :: x, c, nl_c
    REAL(wp), INTENT(IN) :: a_rate
    sdm4_onset_proj = -1.0_wp
    IF (c == 0) THEN
      IF (x <= nl_c) THEN
        IF (x - sdm4_onset_k >= 1) THEN
          sdm4_onset_proj = clr_tau_lig_proj(x - sdm4_onset_k)
        ELSE
          sdm4_onset_proj = 0.0_wp
        END IF
      ELSE IF (sdm4_onset_k == 0) THEN
        sdm4_onset_proj = clr_tau_lig_proj(nl_c)                 ! D2: collar of the flag leaf (ADEL convention)
      ELSE
        sdm4_onset_proj = clr_tau_lig_proj(nl_c) - clr_LED_Cd    ! Malvoisin folio 590: flag-leaf TIP emergence
      END IF
    ELSE
      IF (a_rate <= 0.0_wp) RETURN
      IF (sdm4_onset_k == 0) THEN
        sdm4_onset_proj = ct_tau_birth(c) + REAL(MIN(x, nl_c), wp) / a_rate   ! D3: the synchronous ramp; peduncle with the flag
      ELSE IF (x <= nl_c) THEN
        IF (x - sdm4_onset_k >= 1) THEN
          sdm4_onset_proj = ct_tau_birth(c) + REAL(x - sdm4_onset_k, wp) / a_rate
        ELSE
          sdm4_onset_proj = 0.0_wp
        END IF
      ELSE
        sdm4_onset_proj = ct_tau_em - clr_LED_Cd
      END IF
    END IF
  END FUNCTION sdm4_onset_proj

  ! Driver record 100 mode 3: the ARCHITECTURAL POTENTIAL peduncle length, cm. It
  ! is not an observed length: it is the model's OWN top internode length scaled by
  ! the measured peduncle to penultimate-internode ratio of Rivera-Amado et al.
  ! (2019) Field Crops Res 240:154-167, Table 6 folio 161, applied at the rank at
  ! which it was measured (that table title defines internode 2 as "internode 2
  ! top-down").
  REAL(wp) FUNCTION sdm_Lped_pot_now()
    INTEGER :: nl_ms
    sdm_Lped_pot_now = 0.0_wp
    IF (.NOT. clr_active) RETURN
    nl_ms = cfp_n_top_culm()
    IF (nl_ms < 1) RETURN
    sdm_Lped_pot_now = sdm_ped_int_ratio * cse_Li_at(REAL(nl_ms, wp), nl_ms)
  END FUNCTION sdm_Lped_pot_now

  ! Driver record 100 mode 3: the POTENTIAL whole-plant structural height that
  ! the inverse of Eq. (6) is entitled to be fed, cm. The peduncle is the
  ! topmost internode (Rivera-Amado et al. 2019 Table 6 title, folio 161), and
  ! the published height is measured to the uppermost part of the plant
  ! (Bakhshandeh et al. 2012 p.415), which includes the free blade standing above
  ! the collar; the fit itself is reported valid up to the anthesis date.
  REAL(wp) FUNCTION sdm_h_plant_now(dvs2)
    REAL(wp), INTENT(IN) :: dvs2
    REAL(wp) :: f_ped, span
    sdm_h_collar = sdm_h_pot_now()
    sdm_Lped_pot = sdm_Lped_pot_now()
    sdm_f_ped    = 0.0_wp
    sdm_h_plant_now = sdm_h_collar
    IF (sdm_h_collar <= 0.0_wp) RETURN
    IF (sdm_dvs2_flag < 0.0_wp) RETURN          ! peduncle window has not opened yet
    span = 1.0_wp - sdm_dvs2_flag
    IF (span <= 0.0_wp) RETURN                  ! flag ligule at/after anthesis: no window
    f_ped = (dvs2 - sdm_dvs2_flag) / span
    f_ped = MIN(1.0_wp, MAX(0.0_wp, f_ped))
    sdm_f_ped = f_ped
    sdm_h_plant_now = sdm_h_collar + sdm_Lped_pot * f_ped
  END FUNCTION sdm_h_plant_now

  ! Driver records 109 to 113: SELECT the bracketed members. Pure selection -
  ! no arithmetic on the published values, no midpointing and no collapsing.
  ! Inactive alternative, see Part 1 Table 6.
  SUBROUTINE inl_select_members()
    SELECT CASE (hgt_int_length_member)
    CASE (1) ; inl_Lf_use = inl_Lf_m1     ! Kirby 1988 Table 3 folio 133 (GLASSHOUSE POTS)
    CASE (2) ; inl_Lf_use = inl_Lf_m2     ! Rivera-Amado 2019 Table 6 folio 161 mean (FIELD, 270-275 seeds m-2)
    CASE (3) ; inl_Lf_use = inl_Lf_m3     ! same table, PRINTED minimum row
    CASE (4) ; inl_Lf_use = inl_Lf_m4     ! same table, PRINTED maximum row
    CASE (5) ; inl_Lf_use = inl_Lf_m5     ! cse_li_irank (Abichou) -- the continuity arm, NO peduncle rank
    CASE DEFAULT ; inl_Lf_use = 0.0_wp    ! unreachable: rec 111 = 0 with rec 109 = 1 is an ERROR STOP
    END SELECT
    SELECT CASE (hgt_int_lineic_member)
    CASE (1) ; inl_lam_use = inl_lam_m1   ! Rivera Table 6 "Specific weight" mean (TOTAL DM incl. WSC, mg cm-1)
    CASE (2) ; inl_lam_use = inl_lam_m2   ! same, PRINTED minimum row
    CASE (3) ; inl_lam_use = inl_lam_m3   ! same, PRINTED maximum row
    CASE (4) ; inl_lam_use = inl_lam_m4   ! Gauthier 2020 Table S1.4 SLM_IN (STRUCTURAL, a DIFFERENT quantity)
    CASE DEFAULT ; inl_lam_use = 0.0_wp   ! unreachable: rec 112 = 0 with rec 109 = 1 is an ERROR STOP
    END SELECT
    ! rec 110 = 1-3 onset F (Fournier collar, E1); 4-6 onset M (Malvoisin i-2)
    inl_onset_arm = MERGE(1, 2, hgt_int_window_member <= 3)
    inl_dI_use    = inl_dI_cd(MOD(MAX(1, hgt_int_window_member) - 1, 3) + 1)
  END SUBROUTINE inl_select_members

  ! Driver record 109 mode 1: one day of internode elongation and the height
  ! assembly. Inactive alternative, see Part 1 Table 6.
  SUBROUTINE inl_advance(m_struct_now, axes_m2, dtt_day)
    REAL(wp), INTENT(IN) :: m_struct_now   ! g DM m-2, the model's OWN structural stem mass (W_stem - W_reserve)
    REAL(wp), INTENT(IN) :: axes_m2        ! axes m-2, ct_axes_live (rec 66 is an interlock of rec 109)
    REAL(wp), INTENT(IN) :: dtt_day        ! degCd d-1, the day's thermal time (phen_state%current_dtt)
    INTEGER  :: nl_ms, ir, xabs, x, nt, n_open, n_masslim
    REAL(wp) :: dtau, dM_ax, den, s_sh, cap, pot, ssum, hcol, tend, sum_len_mg

    ! The daily diagnostics are re-zeroed every day, so that a closed day prints
    ! zeros rather than yesterday's numbers.
    inl_dLpot = 0.0_wp ; inl_dLact = 0.0_wp ; inl_bind = 0
    inl_dM_struct = 0.0_wp ; inl_dM_offer = 0.0_wp
    inl_dM_len    = 0.0_wp ; inl_e8_resid = 0.0_wp
    IF (height_internode_law /= 1) RETURN
    inl_axes = axes_m2
    inl_Mstruct_max = MAX(inl_Mstruct_max, m_struct_now)
    IF (.NOT. inl_ready) THEN
      inl_Mstruct_prev = m_struct_now
      inl_ready        = .TRUE.
      RETURN
    END IF
    ! Today's OWN structural stem increment and nothing else. MAX(0,.) is a
    ! positivity guard, since the reserve and structural split can move the
    ! difference down on a day; the raw signed value is emitted beside it.
    inl_dM_struct    = m_struct_now - inl_Mstruct_prev
    inl_Mstruct_prev = m_struct_now
    inl_dM_offer     = MAX(0.0_wp, inl_dM_struct)
    IF (.NOT. clr_active) RETURN                 ! no chronology before emergence: nothing elongates
    nl_ms = cfp_n_top_culm()
    IF (nl_ms < 1) RETURN
    dtau = MAX(0.0_wp, dtt_day)
    CALL inl_select_members()

    ! ---- E1 / E2 / E3 / E4: the windows and the bounded daily potential ----
    DO ir = 0, INL_MAXIR
      inl_absrank(ir) = 0
      IF (inl_Lf_use(ir) <= 0.0_wp) CYCLE
      IF (ir == 0) THEN
        ! The peduncle opens at the FLAG COLLAR under BOTH onset arms.
        inl_absrank(ir) = -1
        inl_tau_on(ir)  = clr_tau_lig_proj(nl_ms)
        IF (hgt_int_ear_member == 3) THEN
          tend = inl_tau_on(ir) + cse_ep_win * clr_phyllochron_Cd    ! E3 peduncle member P2 (Abichou 1.8 phyllochrons)
        ELSE
          tend = inl_tau_on(ir) + inl_dI_use
        END IF
      ELSE

        ! Top-anchored mapping: inversed rank ir is the internode of phytomer
        ! nl_ms - ir + 1. The second onset arm reads the collar two phytomers below.
        xabs = nl_ms - ir + 1
        IF (inl_onset_arm == 2) xabs = xabs - INL_ONSET_SHIFT
        IF (xabs < 1) CYCLE
        inl_absrank(ir) = nl_ms - ir + 1
        inl_tau_on(ir)  = clr_tau_lig_proj(xabs)
        tend            = inl_tau_on(ir) + inl_dI_use
      END IF
      inl_tau_end(ir) = tend
      ! w(x) = 1 if and only if tau_on <= tau < tau_end. The window closes on its
      ! own; there is no freeze gate anywhere in this routine.
      IF (clr_tau_Cd >= inl_tau_on(ir) .AND. clr_tau_Cd < tend) THEN
        cap = MAX(0.0_wp, inl_Lf_use(ir) - inl_L(ir))
        pot = inl_Lf_use(ir) / inl_dI_use * dtau        ! rho_x = L_f(x) / Delta_I, the rectangular window
        IF (pot >= cap) THEN
          inl_dLpot(ir) = cap ; inl_bind(ir) = 3        ! CAP -- the limb the rate and duration bound makes load-bearing
        ELSE
          inl_dLpot(ir) = pot ; inl_bind(ir) = 1        ! POT
        END IF
      END IF
    END DO

    ! ---- E5': the mass-conserving share. NO new demand term, no allometry ----
    den = 0.0_wp
    DO ir = 0, INL_MAXIR
      den = den + inl_lam_use(ir) * inl_dLpot(ir)
    END DO
    IF (den > 0.0_wp) inl_open_days = inl_open_days + 1
    IF (den > 0.0_wp .AND. axes_m2 > 0.0_wp) THEN
      dM_ax      = inl_dM_offer * 1000.0_wp / axes_m2   ! g m-2 -> mg axis-1 (1000 mg per g)
      n_open     = 0
      n_masslim  = 0
      sum_len_mg = 0.0_wp
      DO ir = 0, INL_MAXIR
        IF (inl_dLpot(ir) <= 0.0_wp) CYCLE
        IF (inl_lam_use(ir) <= 0.0_wp) CYCLE
        n_open = n_open + 1
        s_sh   = inl_lam_use(ir) * inl_dLpot(ir) / den
        inl_dLact(ir) = s_sh * dM_ax / inl_lam_use(ir)
        IF (inl_dLact(ir) < inl_dLpot(ir)) THEN
          inl_bind(ir) = 2                              ! MASS -- the day's structural mass could not pay the potential
          n_masslim    = n_masslim + 1
        ELSE

          ! The length increment equals its potential when the day structural mass pays
          ! for the whole potential, and less otherwise; the excess mass is carried.
          inl_dLact(ir) = inl_dLpot(ir)
        END IF
        inl_L(ir)  = inl_L(ir) + inl_dLact(ir)
        inl_m(ir)  = inl_m(ir) + inl_lam_use(ir) * inl_dLact(ir)
        sum_len_mg = sum_len_mg + inl_lam_use(ir) * inl_dLact(ir)
      END DO
      inl_dM_len = sum_len_mg * axes_m2 / 1000.0_wp      ! mg axis-1 -> g m-2
      IF (n_open > 0 .AND. n_masslim == n_open) THEN
        inl_days_mass    = inl_days_mass + 1
        inl_e8_resid     = inl_dM_len - inl_dM_offer
        inl_e8_resid_max = MAX(inl_e8_resid_max, ABS(inl_e8_resid))   ! the G2'-(iii) number, MASS-limb days
      ELSE
        inl_days_pot = inl_days_pot + 1                  ! note sec.14 P3's own signature, counted not judged
        inl_e8_resid = inl_dM_len - inl_dM_offer
      END IF
      ! No mass may ever be CREATED: this must stay 0 on every day of every deck.
      inl_e8_over_max = MAX(inl_e8_over_max, inl_dM_len - inl_dM_offer)
    ELSE
      ! Every window closed, or no live axis: structural deposition adds mass without
      ! length, so the height holds with no freeze gate and no MAX.
      IF (inl_dM_offer > 0.0_wp) inl_days_thicken = inl_days_thicken + 1
      inl_e8_resid = -inl_dM_offer
    END IF

    ! ---- HEIGHT, THREE BASES, NONE PRIVILEGED ----
    ! The first basis is the telescoped collar evaluated on the same ranks, so that
    ! the three height definitions can be compared rather than silently mixed.
    nt   = MIN(INT(MIN(MAX(0.0_wp, clr_HS), REAL(nl_ms, wp))), nl_ms)
    ssum = 0.0_wp ; hcol = 0.0_wp
    DO x = 1, nt
      ir = nl_ms - x + 1
      IF (ir >= 1 .AND. ir <= INL_MAXIR) ssum = ssum + inl_L(ir)
      hcol = MAX(hcol, ssum + cse_Ls_at(REAL(x, wp), nl_ms))
    END DO
    inl_h_b1 = hcol                                      ! BASIS 1: collar / top-leaf, peduncle EXCLUDED
    inl_h_b2 = inl_h_b1 + inl_L(0)                       ! BASIS 2: + the realized peduncle (E2)
    ! The ear length enters at the flag-collar event. The latch is a monotone EVENT
    ! latch on an extrapolated onset, not a height freeze and not a running maximum
    ! of the height.
    IF (inl_tau_on(0) >= 0.0_wp .AND. clr_tau_Cd >= inl_tau_on(0)) inl_ped_open = .TRUE.
    inl_Lear_cm = 0.0_wp
    IF (hgt_int_ear_member >= 2 .AND. inl_ped_open) inl_Lear_cm = inl_Lear_kirby_cm
    inl_h_b3 = inl_h_b2 + inl_Lear_cm                    ! BASIS 3: + spike and awns -- the quantity handed to plant_height
  END SUBROUTINE inl_advance

END MODULE GROWTH_MODULE