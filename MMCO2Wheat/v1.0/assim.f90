!===============================================================================!
! M&M-CO2Wheat v1.0  |  assim.f90
!
! Purpose: Leaf gas exchange for ONE leaf class in ONE hour. Solves the coupled
!          photosynthesis-conductance system: the Farquhar et al. (1980) minimum
!          of the Rubisco-, electron-transport- and triose-phosphate-limited
!          carboxylation rates with the Busch et al. (2018) photorespiratory
!          nitrogen assimilation extension, the Leuning (1995) stomatal
!          conductance, and - where selected - the soil-plant-atmosphere leaf
!          water potential. Spatial unit: one leaf (sunlit, shaded or flag leaf).
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: SHARED (leaf temperature T_leaf_s, absorbed PAR, the Bernacchi and Busch
!       constants, the Leuning and soil-plant-atmosphere parameters, the leaf
!       nitrogen state and the configuration switches); TYPES (working precision).
! Output: Rd (day respiration), A_assim (gross assimilation), An_assim (net
!       assimilation), gsw_1 (stomatal conductance to water vapour), Ci
!       (intercellular CO2), N_total_assim (nitrate reduced in the chloroplast),
!       all per unit leaf area - consumed by CANOPY_MODEL, which integrates them
!       to the canopy and hands the canopy fluxes to the energy balance and to
!       growth_module; the soil-plant-atmosphere diagnostics and the leaf water
!       potential solve counters are written back into SHARED.
! Called by: CANOPY_MODEL.f90 - eval_canopy_assim (sunlit and shaded classes) and
!       tleaf_classes_diag (flag-leaf diagnostic call).
!
!===============================================================================!
SUBROUTINE assim (Rd, A_assim, An_assim, gsw_1, Ca, Ci, N_total_assim, Ds, Narea_in)
  USE SHARED
  USE TYPES, ONLY: wp
  IMPLICIT NONE

  ! Subroutine arguments.
  REAL(wp), INTENT(OUT)   :: Rd, A_assim, An_assim, gsw_1, Ci, N_total_assim
  REAL(wp), INTENT(IN)    :: Ca, Ds
  REAL(wp), INTENT(IN), OPTIONAL :: Narea_in   ! per-class Narea (g N m-2 leaf); mode-2 canopy-N gradient
!----------------------------------------------------------------------!
! Local control for Ci iteration.
!----------------------------------------------------------------------!
INTEGER :: iter_ci
INTEGER :: max_iter_ci
REAL(wp)    :: tol_ci
REAL(wp)    :: Ci_prev
REAL(wp)    :: gsw_leu
REAL(wp)    :: gsc
LOGICAL :: use_fixed_ci
REAL(wp)    :: D_vpd
REAL(wp), PARAMETER :: Ci_min_frac = 0.05_wp
REAL(wp), PARAMETER :: Ci_max_frac = 0.95_wp
! Bisection state for the flux-identity Ci solve. No bound is placed on Ci other
! than the bracket ends, which are themselves verified by evaluation.
REAL(wp)    :: ci_lo, ci_hi, f_ci
INTEGER     :: ci_stage, iters_used
LOGICAL     :: ci_converged
! Leaf water potential solve state. psi_L is the root of
! g(psi_L) = psi_s - E(psi_L)/K_plant - psi_L, found by bracketed regula falsi
! (Illinois) with a bisection safeguard wrapped around the Ci bisection; every
! evaluation is a full Ci solve at the trial f_psi.
LOGICAL  :: spac_on_call
INTEGER  :: ipsi, n_psi, spac_stage, spac_side
REAL(wp) :: fpsi, psiL, psi_s_loc, E_leaf, g_psi, psi_lo, psi_hi, g_lo, g_hi, psi_new
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Parameters.
!----------------------------------------------------------------------!
REAL(wp)    :: N_max        ! ceiling on chloroplast nitrate reduction
REAL(wp)    :: J_max_eff
!----------------------------------------------------------------------!
REAL(wp)    :: tau
REAL(wp)    :: K_C
REAL(wp)    :: K_O
REAL(wp)    :: Vc_max
REAL(wp)    :: beta
REAL(wp)    :: gamma_alph_g  , gamma
REAL(wp)    :: j_alpha_g_max , j_alpha_s_max
REAL(wp)    :: tp_alpha_g_max, tp_alpha_s_max
REAL(wp)    :: c_alpha_g_max , c_alpha_s_max
REAL(wp)    :: alpha_g
REAL(wp)    :: alpha_s
REAL(wp)    :: N_glycine_j  , N_serine_j , N_glycine_c, N_serine_c
REAL(wp)    :: N_glycine_p  , N_serine_p
REAL(wp)    :: N_total_light, N_total_Rub, N_total_TPU
REAL(wp)    :: Phi_o_c
REAL(wp)    :: TPU
REAL(wp)    :: dNum
REAL(wp)    :: dDen
REAL(wp)    :: A_light
REAL(wp)    :: A_Rubisco
REAL(wp)    :: A_sucrose
REAL(wp)    :: A_local
REAL(wp)    :: vcmax25_used   ! 25-C Vcmax, retained for Rd's own T-response
! The selected (minimum) carboxylation rate and its nitrogen-side counterpart,
! each hoisted into a name at the MIN below so the selection can be reported.
REAL(wp)    :: W_sel, W_nsel
!----------------------------------------------------------------------!
! Rubisco kinetics.
!   tau - CO2/O2 specificity factor                        (umol umol-1)
!   K_C - Michaelis-Menten constant for CO2                 (umol mol-1)
!   K_O - Michaelis-Menten constant for O2                  (umol mol-1)
! Bernacchi et al. (2001) Plant Cell Environ 24:253-259, in vivo Rubisco
! kinetics: Eq. (8) p.254 with the c and dHa of Table 1 p.256, on the
! mole-fraction basis of p.253 (Ko is published in mmol mol-1, hence the
! factor 1000 below).  tau is derived from Gamma* through Eq. (5) p.254
! (Gamma* = 0.5 O / tau), so gamma = Oa/(2 tau) below returns Bernacchi's
! Gamma* exactly.
!----------------------------------------------------------------------!
W_sel = 0.0_wp ; W_nsel = 0.0_wp
!----------------------------------------------------------------------!
BLOCK
  REAL(wp) :: Tk_b, gstar_b
  Tk_b    = T_leaf_s + 273.15_wp
  K_C     = EXP(b01_c_Kc - b01_dHa_Kc / (b01_Rgas * Tk_b))              ! umol/mol
  K_O     = 1000.0_wp * EXP(b01_c_Ko - b01_dHa_Ko / (b01_Rgas * Tk_b))  ! -> umol/mol
  gstar_b = EXP(b01_c_Gs - b01_dHa_Gs / (b01_Rgas * Tk_b))              ! umol/mol
  tau     = 0.5_wp * Oa / MAX(1.0e-9_wp, gstar_b)                       ! B01 Eq. (5)
END BLOCK
!----------------------------------------------------------------------!
! Rubisco capacity rate                                (umol m-2 sec-1).
!----------------------------------------------------------------------!
! Vcmax25 and Jmax25 from specific leaf nitrogen 
! (SLN = N_leaf/LAI, g N m-2 leaf):
!
! Vcmax25 = chi_vcmax * (SLN - nb_leaf_v)      Jmax25 = chi_jmax * Vcmax25
!
! The capacities rise linearly with leaf nitrogen per unit area above a
! structural offset that carries no carboxylation capacity.
!
! Source, wheat-specific and from measured gas-exchange A-Ci curves:
!   Silva-Perez et al. (2020) J Exp Bot 71(7):2299-2311,
!   doi:10.1093/jxb/erz439.  Triticum aestivum / turgidum; Vcmax25 normalized
!   to 25 C; Narea is total leaf N per unit area (SLN here). The measured
!   nitrogen use efficiency of carboxylation Peff = Vcmax25/Narea spans
!   26-107 umol CO2 m-2 s-1 (g N)-1 (p.2304 and Fig. 5 p.2306), and Fig. 4B
!   p.2305 shows Vcmax25 near 150 at Narea = 2.5 g N m-2.
!   chi_vcmax = 71 lies inside that measured span. The paper publishes no
!   fitted Vcmax25-Narea regression (Fig. 4B is a genotype-mean cloud with
!   phenotypic correlations of 0.11-0.61 by environment, Table 3), so 71 is a
!   value taken from inside the measured band and not a published slope.
!   chi_jmax = 1.8 is the low edge of the only wheat primary that publishes
!   both parameters, Wullschleger (1993) J Exp Bot 44:907, Table 1 pp.912-913
!   (seven Triticum rows, 1.81-3.22 at 20-25 C, un-normalized); the
!   cross-species 25-C ratio of Medlyn et al. (2002) Plant Cell Environ
!   25:1167 is 1.67 and contains no wheat.
!   nb_leaf_v = 0.40 g N m-2 is structural, non-Rubisco leaf nitrogen. It sits
!   just below the end-of-filling structural residual SLN_min = 0.50 g N m-2 of
!   Bertheloot et al. (2008) Plant Physiol 148:1707, p.1711 and Fig. 2C, and is
!   kept as a separate constant so that it cannot move the nitrogen-senescence
!   floor.
! The floors MAX(5,.) and MAX(10,.) prevent collapse at very low seedling SLN.
BLOCK
  REAL(wp), PARAMETER :: chi_vcmax = 71.0_wp   ! Peff (wheat), Silva-Perez et al. 2020 Fig. 5
  REAL(wp), PARAMETER :: chi_jmax  = 1.8_wp    ! Jmax25/Vcmax25, Wullschleger 1993 Table 1
  REAL(wp), PARAMETER :: nb_leaf_v = 0.40_wp   ! structural leaf N (g N m-2)
  REAL(wp) :: SLN_loc, Vcmax25_N, Jmax25_N
  SLN_loc = 0.0_wp
  IF (gParam_lai > 1.0e-6_wp) SLN_loc = N_leaf / gParam_lai
  ! Per-class Narea overrides canopy-mean SLN only if supplied and >0 (absent or <=0 -> unchanged).
  IF (PRESENT(Narea_in)) THEN
    IF (Narea_in > 0.0_wp) SLN_loc = Narea_in
  END IF
  Vcmax25_N = MAX(5.0_wp,  chi_vcmax            * MAX(0.0_wp, SLN_loc - nb_leaf_v))
  Jmax25_N  = MAX(10.0_wp, chi_jmax * chi_vcmax * MAX(0.0_wp, SLN_loc - nb_leaf_v))
  vcmax25_used = Vcmax25_N   ! Rd's 25-C anchor (the same Vcmax25 in both modes)
!----------------------------------------------------------------------!
  ! Water stress acts on the stomatal SLOPE only, never on the photosynthetic
  ! capacities: Leuning (1995) Plant Cell Environ 18:339-355 p.350.
  ! The slope locus is in the gsw block below, so Vcmax25_N and Jmax25_N are not 
  ! scaled here.
  ! Vc,max temperature response: Bernacchi et al. (2001) Table 1 p.256,
  ! normalized to 1.0 at 25 C (p.255), f(T) = exp(dHa_Vc/R (1/298.15 - 1/Tk))
  ! with dHa_Vc = 65.33 kJ mol-1.  No de-activation term is applied, because
  ! Bernacchi et al. found none in any parameter at T <= 40 C (p.256); the
  ! validity range of the fit is 10-40 C.
!----------------------------------------------------------------------!
  Vc_max = Vcmax25_N * EXP((b01_dHa_Vc / b01_Rgas) * (1.0_wp/298.15_wp - 1.0_wp/(T_leaf_s + 273.15_wp)))
!----------------------------------------------------------------------!
  ! Jmax temperature response: the peaked function of Bernacchi et al. (2003)
  ! Eq. (10) p.1421 with the Jmax,c coefficients of Table 1 p.1424, growth
  ! temperature Tg = 14 C row. It is normalized about the SLN-based Jmax25,
  ! i.e. applied as f(T)/f(25 C).
!----------------------------------------------------------------------!
  BLOCK
    REAL(wp) :: Tk_j, f_T, f_25
    Tk_j = T_leaf_s + 273.15_wp
    f_T  = EXP(b03_c_J - b03_dHa_J / (b01_Rgas * Tk_j)) / &
           (1.0_wp + EXP((b03_dS_J * Tk_j - b03_dHd_J) / (b01_Rgas * Tk_j)))
    f_25 = EXP(b03_c_J - b03_dHa_J / (b01_Rgas * 298.15_wp)) / &
           (1.0_wp + EXP((b03_dS_J * 298.15_wp - b03_dHd_J) / (b01_Rgas * 298.15_wp)))
    J_max_eff = Jmax25_N * f_T / f_25
  END BLOCK
END BLOCK
gamma = Oa / (2.0_wp * tau)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Outer solve on the leaf water potential psi_L (the soil-plant-atmosphere
! continuum of Tuzet et al. 2003 Plant Cell Environ 26:1097, in the form of
! Duursma & Medlyn 2012 Geosci Model Dev 5:919, Eq. 3 p.922).  On the
! unstressed demand pass (spac_pass_on = .FALSE.) the loop body runs ONCE with
! f_psi unused and exits after the Ci solve. Otherwise: stage 1 evaluates at
! psi_L = psi_s (the upper bracket end, where g <= 0 because E >= 0), stage 2
! at psi_L = psi_s - E(psi_s)/K_plant (the lower end, where g >= 0 by
! monotonicity), and stage 3 iterates regula falsi (Illinois) inside that
! bracket, bisecting whenever the secant point leaves it and every fifth step.
! It is converged when |g| < 1e-6 MPa at the evaluated state, so that gsw, An,
! Ci, E and psi_L are one consistent state; non-convergence stops the run. No
! bound is placed on psi_L other than the two verified bracket ends.
!----------------------------------------------------------------------!
spac_on_call = (spac_mode >= 1 .AND. spac_pass_on .AND. water_stress_mode == 1 .AND. use_soil_water_stress)
fpsi = 1.0_wp; psiL = 0.0_wp; psi_s_loc = 0.0_wp; E_leaf = 0.0_wp; g_psi = 0.0_wp
psi_lo = 0.0_wp; psi_hi = 0.0_wp; g_lo = 0.0_wp; g_hi = 0.0_wp
spac_stage = 0; spac_side = 0; n_psi = 0
IF (spac_on_call) THEN
  psi_s_loc = spac_psi_s_MPa
  psiL = psi_s_loc; psi_hi = psi_s_loc; spac_stage = 1
  fpsi = (1.0_wp + EXP(spac_sf * spac_psi_f)) / (1.0_wp + EXP(spac_sf * (spac_psi_f - psiL)))
END IF
!----------------------------------------------------------------------!
! The iteration budget is spac_max_iter + spac_extra_iter. A solve that 
! meets spac_tol_MPa exits at its own iteration; the extra steps are 
! reachable only by a solve that would otherwise have exhausted 
! spac_max_iter and stopped the run.
!----------------------------------------------------------------------!
psi_loop: DO ipsi = 1, spac_max_iter + spac_extra_iter
!----------------------------------------------------------------------!
! Solve for Ci using the Leuning stomatal conductance model.
!----------------------------------------------------------------------!
use_fixed_ci = .FALSE.
max_iter_ci = 50
tol_ci      = 1.0e-4_wp
Ci          = 0.70_wp * Ca
IF (use_fixed_ci) THEN
  Ci        = 0.76_wp * Ca
  max_iter_ci = 1
END IF
Ci          = MIN (Ca * Ci_max_frac, MAX (Ca * Ci_min_frac, Ci))
ci_converged = .FALSE.; ci_stage = 0; iters_used = 0
IF (ci_solve_mode == 1) THEN
  !--------------------------------------------------------------------!
  ! Ci is the root of the flux identity f(Ci) = An(Ci) - gsc(Ci)*(Ca - Ci).
  ! An rises with Ci and gsc*(Ca - Ci) falls with Ci, so f is increasing 
  ! and the root is unique. Lower end: below every admissible Gamma*_alphaG, 
  ! where A_local <= 0 and gsw = g0, so f < 0 there (checked by evaluation, 
  ! stage 2). Upper end: Ca, where f = An(Ca) (stage 1). If An(Ca) <= 0 
  ! (respiring leaf, stomata at g0) the identity's root lies ABOVE Ca and 
  ! the bracket is widened upward until f > 0 -- the exact root, not a 
  ! bound. Convergence is |f| < 1e-6 umol m-2 s-1; failure STOPs the run 
  ! and reports the row.
  !--------------------------------------------------------------------!
  max_iter_ci = 200
  ci_lo    = gamma * (1.0_wp - alpha_g_max)
  ci_hi    = Ca
  ci_stage = 1
  Ci       = ci_hi
END IF
!----------------------------------------------------------------------!

!======================================================================!
!======================================================================!
! 1) Gross photosynthesis limited by light (A_light = Wj) (umol m-2 s-1):
!    RuBP regeneration limits assimilation through the electron transport
!    rate J (J_e here).
!======================================================================!
!======================================================================!
DO iter_ci = 1, max_iter_ci
  Ci_prev = Ci
!----------------------------------------------------------------------!
! Phi_o_c is the ratio of oxygenation to carboxylation, Phi = O/(tau Ci)
! (Harley & Sharkey 1991 Photosynth Res 27:169; Busch et al. 2018 Nature
! Plants 4:46-54, doi:10.1038/s41477-017-0065-x).
!----------------------------------------------------------------------!
  Phi_o_c = (1.0_wp / tau) * (Oa / Ci)
!======================================================================!
! phi is the quantum yield of electron transport, 0.5 (1 - f) of Farquhar
! et al. (1980) Planta 149:78-90, Eq. (8): f = f_photons = 0.23 is the
! fraction of absorbed photons not captured by the light-harvesting
! complexes, and the factor 0.5 is the two photons needed to move one
! electron through linear transport.
!----------------------------------------------------------------------!
phi     = 0.5_wp * (1.0_wp - f_photons)
!----------------------------------------------------------------------!
! The non-rectangular hyperbola for J is that of Harley et al. (1992)
! Plant Cell Environ 15:271-282.  `par` arrives here ALREADY ABSORBED:
! solar.f90 applies the leaf PAR absorptance const_alpha_par when it builds
! Absor_IPAR_*, so absorptance is applied exactly once, there and not here.
!----------------------------------------------------------------------!
J_e     = (phi * par) &
                  / ( ( 1.0_wp + ((phi * par) / J_max_eff)**2 )**0.5_wp )
!----------------------------------------------------------------------!
! beta is the proportion of nitrogen assimilated as glycine, taking the
! ratio of glycine to serine export as constant.
!----------------------------------------------------------------------!
beta = (3.0_wp * alpha_g_max) &
                          / ( (3.0_wp * alpha_g_max) + (2.0_wp * alpha_s_max) )
!----------------------------------------------------------------------!
! N_max is the maximum rate of de novo nitrogen supply to the chloroplast
! (Busch et al. 2018 Eq. 13), read from driver_main.txt.
!----------------------------------------------------------------------!
N_max   = N_max_param ! umol m-2 s-1
!----------------------------------------------------------------------!
! alpha_g and alpha_s are the fractions of glycolate carbon leaving the
! photorespiratory pathway as glycine and as serine, which is how de novo
! nitrogen assimilation enters the light-limited rate Wj.  The two
! expressions below are Busch et al. (2018) Eqs. 18 and 19.
!----------------------------------------------------------------------!
alpha_g = 0.0_wp
alpha_s = 0.0_wp
!----------------------------------------------------------------------!
alpha_g = alpha_g_max
alpha_s = alpha_s_max
!----------------------------------------------------------------------!
IF (J_e.GT.(N_max * ( (2.0_wp * beta) + 6.0_wp)) ) THEN
  alpha_g       = 0.0_wp
  j_alpha_g_max = 0.0_wp
  j_alpha_g_max = ( 4.0_wp * N_max * beta * ( (1.0_wp / Phi_o_c) + 1.0_wp ) ) &
                / ( J_e - ( N_max * ( (2.0_wp * beta) + 6.0_wp ) )      )
  !-----------------------------------------------------!
  alpha_g       = MIN ( alpha_g_max, j_alpha_g_max )
  !-----------------------------------------------------!
  alpha_s       = 0.0_wp
  j_alpha_s_max = 0.0_wp
  j_alpha_s_max = ( 6.0_wp * N_max * (1.0_wp - beta) * ((1.0_wp / Phi_o_c) + 1.0_wp) ) &
                / ( J_e - ( N_max * ( (2.0_wp * beta) + 6.0_wp ) )           )
  !-----------------------------------------------------!
  alpha_s       = MIN ( alpha_s_max, j_alpha_s_max )
  !-----------------------------------------------------!
END IF
!----------------------------------------------------------------------!
! The CO2 compensation point gamma (umol umol-1), adjusted for de novo
! nitrogen assimilation with alpha_g after Busch et al. (2018) Eq. 12.
! Oa is the oxygen concentration of the air (210 000 umol mol-1, SHARED).
!----------------------------------------------------------------------!
! Record 158 diagnostic capture: the alpha_G and alpha_S that the Busch
! N_max produced on the RuBP-limited branch, so the gas-exchange and
! nitrogen diagnostic file can report the size of the electron drain.  The
! value stored is the last canopy element evaluated in the day, which is
! what the daily row prints; nothing reads it back.
!----------------------------------------------------------------------!
IF (nitrate_drain_rate_basis >= 1) THEN ; ndb_ag_last = alpha_g ; ndb_as_last = alpha_s ; END IF
gamma_alph_g  = (0.5_wp * (1.0_wp - alpha_g) * Oa ) &
                                            / tau
!----------------------------------------------------------------------!
dNum    = J_e
dDen    = 4.0_wp + ( ( 4.0_wp + (alpha_g * 8.0_wp) + (alpha_s * 4.0_wp) ) * Phi_o_c )
! // --------------------------------------------------------------- //
! // Light-limited carbon assimilation including de Novo N assimilation
! // using Busch et al. 2017 equations for FvCB
! // (A_light here is Wj in Busch etal. 2017).
! // --------------------------------------------------------------- //
A_light = dNum / dDen
! // --------------------------------------------------------------- //

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
! Calculate N-assimilation as glycine and serine (in umol m-2 s-1):
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
N_glycine_j   =        alpha_g * A_light * Phi_o_c
N_serine_j    = 0.75_wp * alpha_s * A_light * Phi_o_c
!----------------------------------------------------------------------!
! Total Nitrogen assimilated through the photorespiratory pathway
! (in umol m-2 s-1):
!-------------------------------------------------!
N_total_light = N_glycine_j + N_serine_j
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Free alpha_g and alpha_s variables to calculate for TPU limit...
!----------------------------------------------------------------------!
alpha_g = 0.0_wp
alpha_s = 0.0_wp
!----------------------------------------------------------------------!

!======================================================================!
!======================================================================!
! 2) Gross photosynthesis limited by Rubisco capacity (A_Rubisco = Wc):
!    RuBP supply is unlimited and carboxylation is the limit.
!======================================================================!
!======================================================================!
dNum      = 0.0_wp
dDen      = 0.0_wp
dNum      = Vc_max * Ci
dDen      = Ci + ( K_C * ( 1.0_wp + (Oa / K_O) ) )
A_Rubisco = dNum / dDen
!----------------------------------------------------------------------!
! alpha_g and alpha_s enter here for the nitrogen fluxes only; on this
! limb they do not change the carbon assimilation rate. Busch et al.
! (2018) Eqs. 16 and 17 and the Supplementary Information of that paper.
!----------------------------------------------------------------------!
c_alpha_g_max =     ( N_max * beta ) / ( A_Rubisco * Phi_o_c )
c_alpha_s_max =     ( 3.0_wp * N_max * (1.0_wp - beta) ) &
                                      / ( 2.0_wp * ( A_Rubisco * Phi_o_c ) )
!----------------------------------------------------------------------!
alpha_g       = MIN ( alpha_g_max, c_alpha_g_max )
alpha_s       = MIN ( alpha_s_max, c_alpha_s_max )
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
! Calculate N-assimilation as glycine and serine (in umol m-2 s-1):
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
N_glycine_c   =        alpha_g * A_Rubisco * Phi_o_c
N_serine_c    = 0.75_wp * alpha_s * A_Rubisco * Phi_o_c
!----------------------------------------------------------------------!
! Total Nitrogen assimilated through the photorespiratory pathway
! (in umol m-2 s-1):
!-------------------------------------------------!
N_total_Rub     = N_glycine_c + N_serine_c
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Reset alpha_g and alpha_s before the triose-phosphate limb.
!----------------------------------------------------------------------!
alpha_g = 0.0_wp
alpha_s = 0.0_wp
!----------------------------------------------------------------------!

!======================================================================!
!======================================================================!
! 3) Gross photosynthesis limited by triose-phosphate utilization
!    (A_sucrose = Wp), i.e. by the sucrose sink.
!======================================================================!
!======================================================================!
! Tp = 0.11 Vc_max. Busch et al. (2018) Nature Plants 4:46-54, the paper
! this routine already takes Eqs. 12, 13, 16-19 from, prints paired Tp and
! Vc_max values: Fig. 3a p.52 gives Tp = 11.5 with Vcmax = 103.9 (ratio
! 0.1107), Fig. 3c gives 0.1104 and Fig. 4 gives 0.1150, and the Harley
! cotton data span 0.0958-0.1247 over the specific leaf nitrogen this model
! realizes. All of them bracket 0.11; none of them is wheat, so 0.11 is
! taken as the readable evidence for the ratio and no wheat triose-phosphate
! measurement exists to replace it. Because Tp is a fixed multiple of
! Vc_max(T), the ratio Tp/Vc_max is temperature-invariant by construction,
! so the cool-temperature triose-phosphate limitation described by Busch et
! al. (2018) p.48 is not expressed by this formulation.
!----------------------------------------------------------------------!
TPU  = 0.11_wp * Vc_max
!----------------------------------------------------------------------!
! Under triose-phosphate limitation alpha_g and alpha_s follow Busch et al.
! (2018) Eq. 14, which carries the production of amino acids.
!----------------------------------------------------------------------!
tp_alpha_g_max = ( N_max * beta * ((2.0_wp / Phi_o_c) - 1.0_wp) ) &
               / ( (6.0_wp * TPU) + (3.0_wp * N_max * (2.0_wp - beta)) )

tp_alpha_s_max = ( ((3.0_wp/2.0_wp) * N_max) * (1.0_wp - beta) &
                             * ((2.0_wp / Phi_o_c) - 1.0_wp) ) &
               / ( (6.0_wp * TPU) + (3.0_wp * N_max * (2.0_wp - beta)) )
!----------------------------------------------------------------------!
alpha_g        = MIN ( alpha_g_max, tp_alpha_g_max )
alpha_s        = MIN ( alpha_s_max, tp_alpha_s_max )
!----------------------------------------------------------------------!
! Triose-phosphate-limited carbon assimilation including de novo nitrogen
! assimilation; A_sucrose is Wp of Busch et al. (2018).
!----------------------------------------------------------------------!
A_sucrose = (3.0_wp * TPU) &
      / (1.0_wp - (0.5_wp * (1.0_wp + 3.0_wp * alpha_g &
                                     + 4.0_wp * alpha_s) * Phi_o_c))

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
! Calculate N-assimilation as glycine and serine (in umol m-2 s-1):
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
N_glycine_p   =        alpha_g * A_sucrose * Phi_o_c
N_serine_p    = 0.75_wp * alpha_s * A_sucrose * Phi_o_c
!----------------------------------------------------------------------!
! Total Nitrogen assimilated through the photorespiratory pathway
! (in umol m-2 s-1):
!-------------------------------------------------!
N_total_TPU   = N_glycine_p + N_serine_p
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Reset alpha_g and alpha_s after the triose-phosphate limb.
!----------------------------------------------------------------------!
alpha_g = 0.0_wp
alpha_s = 0.0_wp
! ******************************************************************** !
! ********** Selection of the limiting carbon and N rates ************ !
! ******************************************************************** !
! Gross leaf photosynthesis (umol m-2 s-1) is the most limiting of the
! three rates (Farquhar et al. 1980; Busch et al. 2018).
!======================================================================!
! Carbon assimilation:
!======================================================================!
! The MIN is hoisted into a name: W_sel is the value multiplied below, so
! the selected limb is identified from the actual selection and never from a
! recomputation that could disagree with it.
!----------------------------------------------------------------------!
W_sel          = MIN ( ABS(A_light), ABS(A_Rubisco), ABS(A_sucrose) )
A_local        = W_sel &
                                     * ( 1.0_wp - (gamma_alph_g / Ci) )
!======================================================================!
! Nitrogen assimilation:
!======================================================================!
! The same construction on the nitrogen side.
!----------------------------------------------------------------------!
W_nsel         = MIN ( ABS(N_total_light), ABS(N_total_Rub), ABS(N_total_TPU) )
N_total_assim  = W_nsel
! ******************************************************************** !

! Day respiration. Rd carries its own temperature response from Bernacchi
! et al. (2001) Table 1 (dHa_Rd = 46.39 kJ mol-1, normalized at 25 C, p.255),
! not Vc_max's steeper 65.33 kJ mol-1. The 25-C ratio Rd25 = 0.0089 Vcmax25
! is de Pury & Farquhar (1997) Plant Cell Environ 20:537, as tabulated in the
! appendix of Buckley & Adams (2013) Plant Cell Environ 36:1547; it lies
! inside the measured wheat band of 0.006-0.013.
!----------------------------------------------------------------------!
Rd       = 0.0089_wp * vcmax25_used * &
           EXP((b01_dHa_Rd / b01_Rgas) * (1.0_wp/298.15_wp - 1.0_wp/(T_leaf_s + 273.15_wp)))
A_assim  = A_local
An_assim = A_local - Rd 
!----------------------------------------------------------------------!
! Leuning stomatal conductance with Leaf-Surface VPD limitation.
!----------------------------------------------------------------------!
D_vpd   = MAX (10.0_wp, Ds) 
IF (spac_on_call) THEN
  !--------------------------------------------------------------------!
  ! The leaf-water-potential form of the stomatal model: the modifier
  ! f_psi(psi_L) acts on the SLOPE (Leuning 1995 p.350), as in Tuzet et al.
  ! (2003) Plant Cell Environ 26:1097 written in the form of Duursma & Medlyn
  ! (2012) Geosci Model Dev 5:919-940 Eq. 3 p.922; g0 is untouched, and
  ! MAX(g0, .) is the floor at An < 0.
  !
  ! There is deliberately NO vapour-pressure-deficit function h(Ds) on this
  ! branch. In both sources the psi modifier REPLACES the D function rather
  ! than multiplying it: Duursma & Medlyn (2012) p.922, and Tuzet's own 
  ! equation, transcribed by Damour et al. (2010). Plant Cell Environ 
  ! 33:1419-1438 Table 2 Eq. 30 p.1425, is gs = 1.6 [g0 + f(Psi_l) a 
  ! Anet/(Ci - Gamma)], with no D term at all. Damour et al. p.1430 describe 
  ! Tuzet as "close to a modified version of Leuning's model where the slope 
  ! of the gs-Anet relationship is a function of Psi_l" -- 
  ! the slope, with psi instead of D.
  !
  ! The slope therefore matches the form: this branch uses g1_tuzet = 5.7763,
  ! least-squares-fitted to the wheat subset of Lin et al. (2015) Nature Clim
  ! Change 5:459 (n = 245, 95 % CI +/- 0.240, RMSE 0.0926) ON THIS form, i.e.
  ! without h(Ds), and within 4.8 % of the independently fitted Tuzet slope
  ! g1 = 5.5 of Duursma & Medlyn (2012) Table 3 p.930.  The Leuning branches
  ! below keep g1_leu = 10.58 because they DO evaluate h(Ds): a slope is a
  ! property of the functional form it was fitted to, so there is no global g1.
  !--------------------------------------------------------------------!
  gsw_leu = g0_leu + fpsi * g1_tuzet * (An_assim / MAX (1.0e-6_wp, (Ca - gamma_alph_g)))
  gsw_leu = MAX (g0_leu, gsw_leu)
  !--------------------------------------------------------------------!
ELSE IF (ws_two_loci == 1 .AND. water_stress_mode == 1 .AND. use_soil_water_stress) THEN
  !--------------------------------------------------------------------!
  ! Soil-water stress on the SLOPE, the single locus of Leuning (1995) p.350
  ! ("either the slope, a1, or the photosynthetic parameters"): g1 -> f_t g1
  ! with g0 untouched, gs = g0 + f_t g1 An h(Ds)/(Ca - Gamma*).
  !--------------------------------------------------------------------!
  gsw_leu = g0_leu + f_t * g1_leu * (An_assim / MAX (1.0e-6_wp, (Ca - gamma_alph_g))) &
            * (1.0_wp / (1.0_wp + (D_vpd / D0)))
  gsw_leu = MAX (g0_leu, gsw_leu)
END IF
gsc     = gsw_leu / 1.6_wp
!----------------------------------------------------------------------!
! Flux-identity residual at the Ci just evaluated; 
! An and gsc are consistent with it.
!----------------------------------------------------------------------!
f_ci       = An_assim - gsc * (Ca - Ci)
iters_used = iter_ci
IF (ABS (f_ci) < 1.0e-6_wp) THEN
  ci_converged = .TRUE.
ELSE
  SELECT CASE (ci_stage)
  CASE (1)                       ! evaluated at the upper end ci_hi
    IF (f_ci > 0.0_wp) THEN
      ci_stage = 2; Ci = ci_lo
    ELSE                         ! An(ci_hi) < 0: root above; widen upward
      ci_lo = ci_hi; ci_hi = 2.0_wp * ci_hi; Ci = ci_hi
    END IF
  CASE (2)                       ! evaluated at the lower end ci_lo
    IF (f_ci < 0.0_wp) THEN
      ci_stage = 3; Ci = 0.5_wp * (ci_lo + ci_hi)
    ELSE
      WRITE (*,'(A)') ' D1 CI SOLVE: f(ci_lo) >= 0, bracket invalid -- STOP'
      WRITE (*,'(A,3I6,I3,5ES14.6)') ' day,j,hour,class,Ca,ci_lo,ci_hi,f,par = ', &
        ci_census_day, ci_census_j, NINT(ci_census_hour), ci_census_tag, Ca, ci_lo, ci_hi, f_ci, par
      ERROR STOP 77
    END IF
  CASE (3)                       ! bisecting
    IF (f_ci > 0.0_wp) THEN
      ci_hi = Ci
    ELSE
      ci_lo = Ci
    END IF
    IF ((ci_hi - ci_lo) < 1.0e-12_wp * Ca) THEN
      WRITE (*,'(A)') ' D1 CI SOLVE: bracket collapsed before |f| < 1e-6 -- STOP'
      WRITE (*,'(A,3I6,I3,5ES14.6)') ' day,j,hour,class,Ca,ci_lo,ci_hi,f,par = ', &
        ci_census_day, ci_census_j, NINT(ci_census_hour), ci_census_tag, Ca, ci_lo, ci_hi, f_ci, par
      ERROR STOP 77
    END IF
    Ci = 0.5_wp * (ci_lo + ci_hi)
  END SELECT
END IF
IF (ci_converged) EXIT
END DO
IF (ci_solve_mode == 1 .AND. .NOT. ci_converged) THEN
  WRITE (*,'(A)') ' D1 CI SOLVE: no convergence within max_iter_ci -- STOP'
  WRITE (*,'(A,3I6,I3,5ES14.6)') ' day,j,hour,class,Ca,Ci,f,An,gsc = ', &
    ci_census_day, ci_census_j, NINT(ci_census_hour), ci_census_tag, Ca, Ci, f_ci, An_assim, gsc
  ERROR STOP 77
END IF
!----------------------------------------------------------------------!
! Leaf-water-potential residual, evaluated at the Ci-consistent state just solved.
!----------------------------------------------------------------------!
IF (.NOT. spac_on_call) EXIT psi_loop
n_psi  = ipsi
E_leaf = gsw_leu * MAX (0.0_wp, Ds) / spac_P_Pa                 ! (E) mol m-2 s-1
g_psi  = psi_s_loc - E_leaf / spac_K_plant - psiL               ! (H) residual, MPa
IF (ABS (g_psi) < spac_tol_MPa) EXIT psi_loop
SELECT CASE (spac_stage)
CASE (1)                       ! evaluated at psi_hi = psi_s: g = -E/K <= 0
  IF (g_psi > 0.0_wp) THEN
    WRITE (*,'(A)') ' SPAC PSI_L SOLVE: g(psi_s) > 0 (E < 0?) -- STOP'
    WRITE (*,'(A,3I6,I3,4ES14.6)') ' day,j,hour,class,psi_s,E,gsw,Ds = ', &
      ci_census_day, ci_census_j, NINT(ci_census_hour), ci_census_tag, psi_s_loc, E_leaf, gsw_leu, Ds
    ERROR STOP 81
  END IF
  g_hi = g_psi
  psi_lo = psi_s_loc - E_leaf / spac_K_plant                    ! lower bracket end
  spac_stage = 2; psiL = psi_lo
CASE (2)                       ! evaluated at psi_lo: g = (E(psi_s) - E(psi_lo))/K >= 0
  IF (g_psi < 0.0_wp) THEN
    WRITE (*,'(A)') ' SPAC PSI_L SOLVE: g(psi_lo) < 0, bracket invalid (E not monotone?) -- STOP'
    WRITE (*,'(A,3I6,I3,5ES14.6)') ' day,j,hour,class,psi_s,psi_lo,g_lo,g_hi,E = ', &
      ci_census_day, ci_census_j, NINT(ci_census_hour), ci_census_tag, psi_s_loc, psi_lo, g_psi, g_hi, E_leaf
    ERROR STOP 81
  END IF
  g_lo = g_psi; spac_stage = 3; spac_side = 0
  psi_new = (psi_lo * g_hi - psi_hi * g_lo) / (g_hi - g_lo)
  IF (.NOT. (psi_new > psi_lo .AND. psi_new < psi_hi)) psi_new = 0.5_wp * (psi_lo + psi_hi)
  psiL = psi_new
CASE (3)                       ! regula falsi (Illinois) with bisection safeguard
  IF (g_psi > 0.0_wp) THEN
    psi_lo = psiL; g_lo = g_psi
    IF (spac_side == -1) g_hi = 0.5_wp * g_hi
    spac_side = -1
  ELSE
    psi_hi = psiL; g_hi = g_psi
    IF (spac_side == 1) g_lo = 0.5_wp * g_lo
    spac_side = 1
  END IF
  IF ((psi_hi - psi_lo) < 1.0e-12_wp) THEN
    !--------------------------------------------------------------------!
    ! A collapsed bracket with |g| still above spac_tol_MPa means g is 
    ! discontinuous at the Ci solve's own tolerance, so no further iteration 
    ! can help. The state is accepted when the residual is below 
    ! spac_tol_fallback_MPa (10 Pa, three orders of magnitude below the
    ! pressure-chamber resolution of psi_L) and the acceptance is counted; 
    ! above that the run still stops.
    !--------------------------------------------------------------------!
    IF (ABS (g_psi) < spac_tol_fallback_MPa) THEN
      spac_fallback_accepts_run = spac_fallback_accepts_run + 1
      EXIT psi_loop
    END IF
    WRITE (*,'(A)') ' SPAC PSI_L SOLVE: bracket collapsed before |g| < 1e-6 MPa -- STOP'
    WRITE (*,'(A,3I6,I3,5ES14.6)') ' day,j,hour,class,psi_s,psi_lo,psi_hi,g,E = ', &
      ci_census_day, ci_census_j, NINT(ci_census_hour), ci_census_tag, psi_s_loc, psi_lo, psi_hi, g_psi, E_leaf
    ERROR STOP 81
  END IF
  psi_new = (psi_lo * g_hi - psi_hi * g_lo) / (g_hi - g_lo)
  IF (MOD (ipsi, 5) == 0 .OR. .NOT. (psi_new > psi_lo .AND. psi_new < psi_hi)) psi_new = 0.5_wp * (psi_lo + psi_hi)
  psiL = psi_new
END SELECT
fpsi = (1.0_wp + EXP(spac_sf * spac_psi_f)) / (1.0_wp + EXP(spac_sf * (spac_psi_f - psiL)))   ! (F)
END DO psi_loop
IF (spac_on_call) THEN
  ! The hydraulic identity psi_s - E/K_plant - psi_L = 0 is asserted at the final state.
  IF (ABS (psi_s_loc - E_leaf / spac_K_plant - psiL) >= spac_tol_MPa) THEN
    ! Bounded fallback acceptance below spac_tol_fallback_MPa, counted; otherwise the run stops.
    IF (ABS (psi_s_loc - E_leaf / spac_K_plant - psiL) < spac_tol_fallback_MPa) THEN
      spac_fallback_accepts_run = spac_fallback_accepts_run + 1
    ELSE
    WRITE (*,'(A)') ' SPAC PSI_L SOLVE: no convergence within spac_max_iter -- STOP'
    WRITE (*,'(A,3I6,I3,I4,5ES14.6)') ' day,j,hour,class,iters,psi_s,psi_L,g,E,gsw = ', &
      ci_census_day, ci_census_j, NINT(ci_census_hour), ci_census_tag, n_psi, psi_s_loc, psiL, g_psi, E_leaf, gsw_leu
    ERROR STOP 81
    END IF
  END IF
  IF (ci_census_tag > 0) THEN
    ! Daily leaf-water-potential solve counters, on the writer pass only.
    spac_calls_day = spac_calls_day + 1
    spac_iters_max_day = MAX (spac_iters_max_day, n_psi)
    spac_resid_max_day = MAX (spac_resid_max_day, ABS (g_psi))
  END IF
END IF
gsw_1 = gsw_leu

!----------------------------------------------------------------------!
END SUBROUTINE assim
!======================================================================!