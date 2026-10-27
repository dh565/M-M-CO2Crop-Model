!===============================================================================!
! M&M-CO2Wheat v1.0  |  soiln_module.f90
!
! Purpose: Crop-residue decomposition and the microbial carbon-nitrogen turnover
!          that feeds the soil nitrogen chain. Decays a fast (metabolic) and a
!          slow (structural) residue pool by first-order kinetics with
!          temperature and water factors, and resolves, from the carbon-use
!          efficiency and the microbial C:N ratio, whether each pool's decay
!          mineralizes nitrogen to ammonium or immobilizes it from the mineral
!          pools; microbial biomass then turns over and re-mineralizes its
!          nitrogen. Time step: daily. Spatial unit: the residue pools and the
!          mineral nitrogen of the soil as bulk stocks (g m-2).
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: TYPES (working precision wp); SHARED, ONLY the residue carbon and
!       nitrogen pools res_fast_C/N and res_slow_C/N, the microbial pools
!       soil_mic_C/N, the mineral pools N_NH4_soil (ammonium) and N_soil_min
!       (nitrate), and the decomposition parameters k_fast_ref, k_slow_ref,
!       Q10_dec, T_ref_dec, e_mic (carbon-use efficiency), CN_mic (microbial C:N)
!       and k_mic (microbial turnover rate).
! Output: the residue, microbial and mineral pools are updated in place, and the
!       day's fluxes are returned: decomposed carbon and nitrogen, respired
!       carbon, net mineralization to ammonium, immobilization from ammonium and
!       from nitrate, and microbial turnover nitrogen (all g m-2 d-1).
!       main_MMCO2Crop.f90 consumes them: the ammonium gained enters the
!       NH4 -> NO3 nitrification pathway and, in layered mode, is distributed
!       over the layered NH4/NO3 stocks; the respired carbon drives
!       denitrification; the fluxes are written to the soil-nitrogen balance and
!       its closure check.
! Called by: main_MMCO2Crop.f90 - once per day in the soil-nitrogen sequence, and
!       once per pre-sowing day of the residue spin-up.
!===============================================================================!
! Mass-conserving first-order decomposition on the CENTURY pool structure
! (Parton et al., 1988). The microbial carbon-use efficiency and the microbial
! C:N ratio follow doi:10.5194/gmd-11-4779-2018 and the microbial carbon-use
! efficiency review PMC12024887.
! Decomposition rate k(T,water) = k_ref * Q10_dec^((Tsoil-T_ref_dec)/10) * f_water.
!===============================================================================!
MODULE SOILN_MODULE
  USE TYPES, ONLY: wp
  USE SHARED, ONLY: res_fast_C, res_fast_N, res_slow_C, res_slow_N, &
    soil_mic_C, soil_mic_N, N_NH4_soil, N_soil_min, &
    k_fast_ref, k_slow_ref, Q10_dec, T_ref_dec, e_mic, CN_mic, k_mic
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: decompose_residue_day

CONTAINS

  ! One day of residue decomposition + microbial turnover. Operates on the SHARED
  ! pools; returns daily diagnostic fluxes. T_soil_c in C, f_water in [0,1].
  SUBROUTINE decompose_residue_day(T_soil_c, f_water, decomp_C, decomp_N, respired_C, &
                                   net_min_NH4, immob_NH4, immob_NO3, turnover_N)
    REAL(wp), INTENT(IN)  :: T_soil_c, f_water
    REAL(wp), INTENT(OUT) :: decomp_C, decomp_N, respired_C, net_min_NH4, immob_NH4, immob_NO3, turnover_N
    REAL(wp) :: fT, fw, kf, ks, kmt, dMicC, dMicN, T_dec

    decomp_C = 0.0_wp; decomp_N = 0.0_wp; respired_C = 0.0_wp
    net_min_NH4 = 0.0_wp; immob_NH4 = 0.0_wp; immob_NO3 = 0.0_wp; turnover_N = 0.0_wp

    ! Temperature response. Microbial decomposition is enzyme-mediated and has a
    ! temperature optimum rather than an indefinitely exponential response, so the
    ! temperature entering the Q10 is capped at 45 C: CENTURY uses Topt = 35 C and
    ! Tmax = 45 C (Burke et al., 2003); see also Davidson & Janssens (2006),
    ! Schipper et al. (2014) and Liu et al. (2018), whose measured optima span
    ! 38.5-46.0 C with a mean of 42.4 C. The cap applies to the temperature only,
    ! giving fT <= Q10_dec**((45-T_ref_dec)/10) = 4 at Q10 = 2 and T_ref = 25 C.
    ! Without it, soil-temperature spikes can decay the whole residue and
    ! microbial pools in a single day.
    T_dec = MIN(T_soil_c, 45.0_wp)
    fT = Q10_dec ** ((T_dec - T_ref_dec) / 10.0_wp)
    fw = MIN(1.0_wp, MAX(0.0_wp, f_water))
    kf = MAX(0.0_wp, k_fast_ref * fT * fw)
    ks = MAX(0.0_wp, k_slow_ref * fT * fw)

    CALL decay_pool(res_fast_C, res_fast_N, kf, decomp_C, decomp_N, respired_C, &
                    net_min_NH4, immob_NH4, immob_NO3)
    CALL decay_pool(res_slow_C, res_slow_N, ks, decomp_C, decomp_N, respired_C, &
                    net_min_NH4, immob_NH4, immob_NO3)

    ! Microbial turnover: biomass dies -> C respired (CO2), N re-mineralized to NH4.
    kmt = MAX(0.0_wp, k_mic * fT * fw)
    dMicC = MIN(soil_mic_C, kmt * soil_mic_C)
    dMicN = 0.0_wp
    IF (soil_mic_C > 1.0e-12_wp) dMicN = dMicC * soil_mic_N / soil_mic_C
    soil_mic_C = soil_mic_C - dMicC
    soil_mic_N = soil_mic_N - dMicN
    N_NH4_soil = N_NH4_soil + dMicN
    turnover_N = dMicN
    respired_C = respired_C + dMicC
  END SUBROUTINE decompose_residue_day

  ! One first-order decay step of a single residue pool, with the nitrogen
  ! bookkeeping that follows from it. The decayed carbon is split between
  ! microbial biomass (the fraction e_mic, the carbon-use efficiency) and
  ! respiration; the nitrogen that the new biomass requires is micCg/CN_mic, so
  ! the sign of the difference between the nitrogen released by the pool and that
  ! requirement decides between net mineralization and immobilization. Microbial
  ! C:N is preserved in every branch, so no nitrogen stress is hidden.
  SUBROUTINE decay_pool(rC, rN, k, dC_acc, dN_acc, resp_acc, min_acc, immNH4_acc, immNO3_acc)
    REAL(wp), INTENT(INOUT) :: rC, rN
    REAL(wp), INTENT(IN)    :: k
    REAL(wp), INTENT(INOUT) :: dC_acc, dN_acc, resp_acc, min_acc, immNH4_acc, immNO3_acc
    REAL(wp) :: dC, dN, micCg, micNneed, need, t1, t2, micNact, micCact

    IF (rC <= 1.0e-12_wp) RETURN
    dC = MIN(rC, k * rC)
    dN = dC * rN / rC
    rC = rC - dC; rN = rN - dN
    dC_acc = dC_acc + dC; dN_acc = dN_acc + dN

    micCg    = e_mic * dC                 ! C the microbes would assimilate
    micNneed = micCg / CN_mic             ! N required for that biomass
    IF (dN >= micNneed) THEN
      ! Net mineralization: residue supplies enough N; excess -> NH4.
      soil_mic_C = soil_mic_C + micCg
      soil_mic_N = soil_mic_N + micNneed
      N_NH4_soil = N_NH4_soil + (dN - micNneed)
      min_acc    = min_acc    + (dN - micNneed)
      resp_acc   = resp_acc   + (dC - micCg)
    ELSE
      ! Immobilization: draw the deficit from NH4 then NO3.
      need = micNneed - dN
      t1 = MIN(N_NH4_soil, need); N_NH4_soil = N_NH4_soil - t1; need = need - t1
      t2 = MIN(N_soil_min, need); N_soil_min = N_soil_min - t2; need = need - t2
      immNH4_acc = immNH4_acc + t1; immNO3_acc = immNO3_acc + t2
      micNact = dN + t1 + t2
      IF (need <= 0.0_wp) THEN
        ! Enough mineral N: full microbial growth.
        soil_mic_C = soil_mic_C + micCg
        soil_mic_N = soil_mic_N + micNneed
        resp_acc   = resp_acc   + (dC - micCg)
      ELSE
        ! N-limited: cap microbial C to the N actually available (preserving
        ! CN_mic) and respire the rest of the decayed carbon.
        micCact    = micNact * CN_mic
        soil_mic_C = soil_mic_C + micCact
        soil_mic_N = soil_mic_N + micNact
        resp_acc   = resp_acc   + (dC - micCact)
      END IF
    END IF
  END SUBROUTINE decay_pool

END MODULE SOILN_MODULE