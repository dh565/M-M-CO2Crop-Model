!===============================================================================!
! M&M-CO2Wheat v1.0  |  meteo.f90
!
! Purpose: Builds the hourly weather forcing from one daily weather record:
!          solar geometry, hourly global, direct and diffuse irradiance, hourly
!          air temperature and relative humidity (including the sunset anchor of
!          the night branch and the offset that holds the 24-hour air-temperature
!          mean at the daily mean), vapor pressure and vapor pressure deficit,
!          sky emissivity and bulk net radiation. Time step: the daily record
!          expanded to hourly (one call per hour of the day). Spatial unit:
!          point, at the above-canopy reference height.
!
! Author:      David Helman, The Hebrew University of Jerusalem
! Last update: 2026-10-20
!
!  Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
!  Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): 
!  non-commercial research use with citation; redistribution of modified versions 
!  and commercial use require the copyright holders' consent.
!
! Uses: SHARED (daily weather record, solar-time state computed by th, canopy
!       LAI, physical constants, and the hourly state this routine writes back),
!       TYPES (working precision). Calls no other routine.
! Output: SolarElev, dfrad_h, drrad_h -> solar (two-leaf PAR partition) and
!         CANOPY_MODEL; HourTemp, Svp_s, HourVP_s, VpdRef_met -> CANOPY_MODEL
!         gas exchange and energy balance; sky_emiss_comp -> CANOPY_MODEL
!         component longwave (ebh_LW_down / LW_down_s); dTotrad (daily global
!         irradiance, J m-2 d-1) -> growth_module and main_MMCO2Crop;
!         alpha_mix_s, HourSW_net_s, HourLW_net_s and HourNetRad -> the hourly
!         radiation and aerodynamics diagnostics.
! Called by: CANOPY_MODEL (canopy_model): once per hour over 24 hours to form
!            the air-temperature offset, then once per committed time step.

!===============================================================================!
SUBROUTINE meteo (SolarElev, dfrad_h, drrad_h, VpdRef_met, HourTemp, HourNetRad)
!------------------------------------------------------------------------------!
USE SHARED   ! Shared variables.
USE TYPES, ONLY: wp
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
IMPLICIT NONE
!------------------------------------------------------------------------------!
! Dummy arguments and local variables.
!------------------------------------------------------------------------------!
REAL(wp)    :: SolarElev
REAL(wp)    :: HourETRad
REAL(wp)    :: DayETRad
REAL(wp)    :: Sunset
REAL(wp)    :: Sunrise

REAL(wp)    :: dfrad
REAL(wp)    :: drrad
REAL(wp)    :: dfrad_h
REAL(wp)    :: drrad_h

REAL(wp)    :: Svp
REAL(wp)    :: HourVP
REAL(wp)    :: VpdRef_met
REAL(wp)    :: HourTemp, HourRH
REAL(wp)    :: HourNetRad

REAL(wp)    :: dTau
REAL(wp)    :: dTset, dTset_RH
REAL(wp)    :: dSlope, dSlope_RH

REAL(wp)    :: dTa
REAL(wp)    :: dTs
REAL(wp)    :: dSW_net
REAL(wp)    :: dLW_out
REAL(wp)    :: dLW_in
REAL(wp)    :: dLW_net

REAL(wp)    :: dAoB
REAL(wp)    :: dHa
REAL(wp)    :: d_parm
REAL(wp)    :: dE0
REAL(wp)    :: dPhi
REAL(wp)    :: dCoefA
REAL(wp)    :: dCoefB
REAL(wp)    :: dSinb
REAL(wp)    :: dR
REAL(wp)    :: dK

REAL(wp)    :: dETRad
REAL(wp)    :: dSinbe
REAL(wp)    :: dScc

REAL(wp)    :: dRelSunhr
REAL(wp)    :: dTotrad_h
REAL(wp)    :: dTrans
REAL(wp)    :: dDfFrac
REAL(wp)    :: dDayTotRad
REAL(wp)    :: dHourET

! *****************************************************************************!
! 1. Solar radiation.
! *****************************************************************************!
! The equation numbers (2.x) used through this routine are those of Teh (2006),
! Introduction to Mathematical Modeling of Crop Growth, Ch. 2, from which the
! daily-to-hourly disaggregation is transcribed.
! ---------------------------------------------------------------------------- !
! Solar elevation beta (radians) from the solar declination (computed in th),
! the site latitude and the hour angle                                   (2.7).
! Negative when the sun is below the horizon.
! ---------------------------------------------------------------------------- !
dA = zero
dB = zero
dA        = sin (Solar_Decl) * sin (lat_s_rad)
dB        = cos (Solar_Decl) * cos (lat_s_rad)
dHa       = Param_pi * (th_s - 12.0_wp) / 12.0_wp
SolarElev = asin (dA + ( dB * cos(dHa) ) )         ! values may be -ve
! ---------------------------------------------------------------------------- !
! Sunset solar time (hour)                                               (2.11).
! ---------------------------------------------------------------------------- !
Sunrise   = dTsr ! this is from th.f90
Sunset    = (24.0_wp - Sunrise)
! ---------------------------------------------------------------------------- !
! Hourly extraterrestrial solar irradiance (W m-2).
! ---------------------------------------------------------------------------- !
dETRad  = zero
! from Spencer (1971).                                                   (2.20).
d_parm  = 2.0_wp * Param_pi * ( (REAL(DOY, wp) - 1.0_wp)/365.0_wp )
dE0     = 1.00011_wp + 0.034221_wp * cos (d_parm) + 0.00128_wp * sin (d_parm) &
                  + 0.000719_wp * cos (2.0_wp * d_parm) + 0.000077_wp * sin (2.0_wp * d_parm)
dETRad  = 1370.0_wp * dE0 * sin (SolarElev)
IF (dETRad < zero) THEN
     dETRad = zero         !           don't allow -ve radiation values.
END IF
HourETRad   = dETRad
! ---------------------------------------------------------------------------- !
! Daily extraterrestrial solar irradiance (J m-2 day-1).
! ---------------------------------------------------------------------------- !
dA   = zero
dB   = zero
dE0  = zero
dA       = sin (Solar_Decl) * sin (lat_s_rad)
dB       = cos (Solar_Decl) * cos (lat_s_rad)
dAoB     = dA / dB
dSinbe   = (24.0_wp / Param_pi) * (dA * acos (-dAoB) + dB &
                * sqrt (1.0_wp - (dAoB * dAoB)))
dE0      = 1.0_wp + 0.033_wp * cos (2.0_wp * Param_pi * (REAL(DOY, wp) - 10.0_wp) / 365.0_wp)
dScc     = 1370.0_wp * dE0
DayETRad = 3600.0_wp * dScc * dSinbe
! ---------------------------------------------------------------------------- !
! Daily incoming shortwave irradiance (J m-2 day-1) from input (2.14).
! dfrad (diffuse irradiance), drrad (direct irradiance).
! ---------------------------------------------------------------------------- !
dRelSunhr = s_hr_s / Daylength
dTotrad   = dTotrad_s * 1.E6_wp ! measured daily global irradiance, MJ m-2 d-1 -> J m-2 d-1
IF (meteo_debug) THEN
  print *, '++++++++++++++++++++++'
  print *, 'dTotRad = ', dTotrad
  print *, '++++++++++++++++++++++'
END IF
dTrans    = dTotrad / DayETRad
dDfFrac   = one
!------------------------------------------------------------------------------!
! Daily diffuse fraction as a function of the atmospheric transmissivity
! dTrans = dTotrad / DayETRad (de Jong, 1980)                            (2.25).
! ---------------------------------------------------------------------------- !
IF (dTrans >= 0.75_wp) THEN
  dDfFrac = 0.23_wp
END IF

IF (dTrans < 0.75_wp .AND. dTrans >= 0.35_wp) THEN
  dDfFrac = 1.33_wp - 1.46_wp * dTrans
END IF

IF (dTrans < 0.35_wp .AND. dTrans >= 0.07_wp) THEN
  dA      = dTrans - 0.07_wp
  dDfFrac = 1.0_wp - 2.3_wp * dA * dA
END IF

IF (dTrans < 0.07_wp) THEN
  dDfFrac = 1.0_wp
END IF

dfrad     = dDfFrac * dTotrad
drrad     = dTotrad - dfrad
! ---------------------------------------------------------------------------- !
! Hourly solar irradiance (W m-2)                                   (2.27-2.30).
! dfrad (diffuse irradiance), drrad (direct irradiance).
! ---------------------------------------------------------------------------- !
dA   = zero
dB   = zero
dAoB = zero
dDayTotRad = dfrad + drrad
dA         = sin (Solar_Decl) * sin (lat_s_rad)
dB         = cos (Solar_Decl) * cos (lat_s_rad)
dAoB       = dA / dB
dPhi       = (Param_pi * dDayTotRad / 86400.0_wp) &
              / (dA * acos (-dAoB) + dB * sqrt (1.0_wp - dAoB * dAoB))
dCoefA     = -dB * dPhi
dCoefB     =  dA * dPhi
dTotrad_h  = dCoefA * cos (Param_pi * th_s / 12.0_wp) + dCoefB
dSinb      = sin (SolarElev)
! ---------------------------------------------------------------------------- !
! Hourly diffuse fraction from the hourly transmissivity and the sine of the
! solar elevation (Spitters et al., 1986)                                (2.31).
! ---------------------------------------------------------------------------- !
dDfFrac     = one                                        ! fraction diffuse.
! ---------------------------------------------------------------------------- !
IF (dTotrad_h > 0.0_wp .AND. dSinb >= 0.0_wp) THEN
  dR        = 0.847_wp - 1.61_wp * dSinb + 1.04_wp * dSinb * dSinb!               (2.22).
  dK        = (1.47_wp - dR) / 1.66_wp
  dHourET   = HourETRad                                  ! can't be zero.
  dTrans    = dTotrad_h / dHourET
  IF (dTrans > dK) THEN
    dDfFrac = dR
  END IF

  IF (dTrans <= dK .AND. dTrans > 0.35_wp) THEN
    dDfFrac = 1.47_wp - 1.66_wp * dTrans
  END IF

  IF (dTrans <= 0.35_wp .AND. dTrans > 0.22_wp) THEN
    dA      = dTrans - 0.22_wp
    dDfFrac = 1.0_wp - 6.4_wp * dA * dA
  END IF

END IF

IF (dTotrad_h < 0.0_wp) THEN
  dTotrad_h = 0.0_wp
END IF

dfrad_h  = dDfFrac * dTotrad_h
drrad_h  = dTotrad_h - dfrad_h
IF (meteo_debug) THEN
  print *, '$$$$$$$$$$$$$$$$$$$$$$$$       dTotrad = ', dTotrad_h, 'W m-2'
END IF

! ---------------------------------------------------------------------------- !

! *****************************************************************************!
! 2. Hourly wind speed (m s-1).
! *****************************************************************************!
! Wind speed is held at its daily mean: no hourly wind data are read.

! *****************************************************************************!
! 3. Hourly air temperature (deg C) and relative humidity (%).
! *****************************************************************************!
! Air temperature follows a sine curve from sunrise to sunset and a linear
! decline through the night, from the daily maximum and minimum         (2.47).
! Relative humidity follows the same curve in opposite phase, since it is
! negatively correlated with temperature.
! -----------------------------------------------------------------------------!
IF (th_s >= (Sunrise + P_temp) .AND. th_s <= Sunset) THEN
  ! Daytime branch.
  dTau     = Param_pi * ( (th_s - (Sunrise + P_temp)) / (Daylength) )
  HourTemp = Tmin_d_s + (Tmax_d_s - Tmin_d_s) * sin (dTau)
  HourRH   = RH_max_s - (RH_max_s - RH_min_s) * sin (dTau)
ELSE
  ! ---------------------------------------------------------------------------!
  ! Night branch. Hours before sunrise are advanced by 24 h so that the decline
  ! from the previous sunset is monotonic in th_s; the caller resets th_s at
  ! every hourly step.
  ! ---------------------------------------------------------------------------!
  IF (th_s < (Sunrise + P_temp)) THEN
    th_s = th_s + 24.
  END IF
  dTau      = Param_pi  * ( (Daylength - P_temp)  / (Daylength) )
  ! ---------------------------------------------------------------------------!
  ! Anchors of the night decline: the day curves evaluated at sunset. Daylength
  ! = Sunset - Sunrise, so this dTau is the daytime dTau at sunset and the two
  ! lines below reproduce the daytime temperature and humidity lines there.
  ! ---------------------------------------------------------------------------!
  dTset     = Tmin_d_s  + ( (Tmax_d_s  - Tmin_d_s) * sin (dTau) )
  ! ---------------------------------------------------------------------------!
  ! The humidity anchor is RH_max_s, the daytime amplitude term being anti-phased
  ! with temperature: relative humidity peaks near sunrise and reaches its
  ! minimum in the early afternoon (FAO-56; Allen et al., 1998, Ch. 3,
  ! ISBN 92-5-104219-5). Anchoring on RH_max_s makes hourly relative humidity
  ! continuous across sunset. This scheme disaggregates relative humidity rather
  ! than holding the vapor pressure e_a constant over the day.
  ! ---------------------------------------------------------------------------!
  dTset_RH  = RH_max_s  - ( (RH_max_s  - RH_min_s) * sin (dTau) )

  ! Linear interpolation from the sunset anchors to the next dawn.
  dSlope    = (Tmin_d_s - dTset)    / (Sunrise + P_temp + 24. - Sunset)
  dSlope_RH = (RH_max_s - dTset_RH) / (Sunrise + P_temp + 24. - Sunset)

  HourTemp  = dTset    + dSlope    * (th_s - Sunset)
  HourRH    = dTset_RH + dSlope_RH * (th_s - Sunset)
END IF
! Apply daily temperature offset so 24-hr mean matches input mean temp_s.
HourTemp = HourTemp + temp_offset
IF (HourRH > RH_max_s) THEN
  HourRH = RH_max_s
END IF
IF (meteo_debug) THEN
  print *, 'hourly temp ********************************************** = ', HourTemp, th_s, Tmin_d_s, Tmax_d_s
  print *, 'hourly RH   ********************************************** = ', HourRH  , th_s, RH_min_s, RH_max_s
END IF
! *****************************************************************************!
! 4. Saturated vapor pressure & VPD (hourly from temperature and RH).
! *****************************************************************************!
! Saturation vapor pressure e_s (Pa) from the hourly air temperature, Tetens
! (1930) equation                                                       (2.40).
! Both branches evaluate the same Tetens expression: the separate sub-zero form
! of Murray (1967) is not applied.
IF (HourTemp < 0.) THEN
  Svp  = 610.78 * exp (17.269 * (HourTemp / (HourTemp + 237.3)))
ELSE
  Svp  = 610.78 * exp (17.269 * (HourTemp / (HourTemp + 237.3)))
END IF
! Hourly vapor pressure e_a (Pa) from the hourly relative humidity, which is
! first clipped to [0, 100] %.
HourRH = MIN (100.0, MAX (0.0, HourRH))
HourVP = (HourRH / 100.) * Svp
! Export e_a and e_s to SHARED: they are the reference humidity of the
! three-node energy balance (e_ref = HourVP_s).
Svp_s    = Svp
HourVP_s = HourVP
! ---------------------------------------------------------------------------- !
! Effective sky emissivity for the component longwave path.
! Clear sky: Brutsaert (1975, Water Resour. Res. 11:742-744), Eq. (11); the
! formula and its units are also given by Crawford & Duchon (1999), Eq. 16,
! p. 476:
!   eps_clear = 1.24 * (e/T)^(1/7)   [T in kelvin, e in millibars]
! Cloud modification (Crawford & Duchon, 1999, Eq. 11 with Eq. 3):
!   eps_sky = clf + (1 - clf)*eps_clear,  clf = 1 - s,  clf clipped to [0,1],
!   s = Rs/Rso. Two deviations follow from the available site data:
!   (i) s is formed from the DAILY measured global radiation and carried
!       through the night hours, where Crawford & Duchon use daylight only;
!   (ii) Rso = 0.75*Ra (FAO-56, Allen et al. 1998, Eq. 37, dropping the
!       elevation term 2e-5*z, which is below 0.2 % at these sites).
! No coefficient here is fitted. sky_emiss_comp is consumed only by the
! component radiation path (CANOPY_MODEL ebh_LW_down / LW_down_s).
! ---------------------------------------------------------------------------- !
BLOCK
  REAL(wp) :: T_K_a1, e_mb_a1, eps_clear_a1, s_a1, clf_a1
  T_K_a1 = HourTemp + 273.15_wp
  e_mb_a1 = HourVP / 100.0_wp
  eps_clear_a1 = 1.24_wp * (MAX(1.0e-3_wp, e_mb_a1) / T_K_a1) ** (1.0_wp / 7.0_wp)
  s_a1 = 0.0_wp
  IF (DayETRad > 0.0_wp) s_a1 = (dTotrad / DayETRad) / 0.75_wp
  s_a1  = MIN(1.0_wp, MAX(0.0_wp, s_a1))
  clf_a1 = 1.0_wp - s_a1
  sky_emiss_comp = MIN(1.0_wp, clf_a1 + (1.0_wp - clf_a1) * eps_clear_a1)
END BLOCK
IF (meteo_debug) THEN
  WRITE (*,*) 'Svp hourly', Svp
  WRITE (*,*) 'VPD hourly', HourVP
END IF

! Vapor pressure deficit at the reference height (Pa).
VpdRef_met = Svp - HourVP
! *****************************************************************************!
! 5. Hourly bulk net radiation. th_s is local solar time.
! *****************************************************************************!
! LAI-dependent mixed surface albedo, Beer-Lambert blending of the soil and
! canopy albedos:
!   alpha_mix = alpha_soil*exp(-k*LAI) + alpha_can*(1 - exp(-k*LAI))
! At LAI = 0 this is 0.28 (bare sandy soil); at full canopy it tends to 0.20
! (wheat). Sources: Karnieli et al. (1996, RSE 56:1-12, Table 2) for the soil;
!          Monteith & Unsworth (2008, p. 48, Table 3.1) for the crop;
!          Goudriaan (1977) for the Beer-Lambert blending form.
! -----------------------------------------------------------------------------!
alpha_mix_s = alpha_soil_dry_sandy * EXP(-k_alpha_canopy * MAX(0.0_wp, gParam_lai)) + &
              alpha_canopy_wheat   * (1.0_wp - EXP(-k_alpha_canopy * MAX(0.0_wp, gParam_lai)))
dSW_net = (1.0_wp - alpha_mix_s) * (dfrad_h + drrad_h)
IF (meteo_debug) THEN
  print *, 'SWrad = ', dSW_net
END IF
dTa     = HourTemp + 273.15                 ! convert to kelvin
dTs     = dTa                               ! assume surface T = air T
dLW_out = sur_emiss * const_SB * (dTs** 4)  ! sur_emiss = 1, a black body

! -----------------------------------------------------------------------------!
! Bulk incoming longwave. The physics_mode >= 1 branch is INACTIVE in v1.0:
! physics_mode is held at 0 (main_MMCO2Crop.f90), so the branch taken is always
! the ELSE below. The active component longwave path uses sky_emiss_comp above.
! Sources of the inactive branch:
! - Brutsaert, W. (1975), On a derivable formula for long-wave radiation from
!   clear skies, Water Resour. Res. 11:742-744:
!   eps_clear = 1.24 * (e / Ta)^(1/7), Ta in K. Note that the 1.24 coefficient
!   is for e in millibars, whereas this branch supplies e in kPa; the
!   component path above uses millibars.
! - Allen et al. (1998), FAO-56, Eq. 39: cloudiness correction
!   (1 + 0.22 * (1 - Rs/Rso)^2).
! ---------------------------------------------------------------------------- !
IF (physics_mode >= 1) THEN
  BLOCK
    REAL(wp) :: ea_kpa, emiss_clear, emiss_raw
    REAL(wp) :: rs_rso, cloud_fac
    ea_kpa = HourVP / 1000.0_wp
    emiss_clear = 1.24_wp * (ea_kpa / dTa) ** (1.0_wp / 7.0_wp)
    cloud_fac = 1.0_wp
    ! Rs/Rso uses hourly shortwave over hourly extraterrestrial shortwave.
    IF (HourETRad > 0.0_wp .AND. dTotrad_h > 0.0_wp) THEN
      rs_rso = dTotrad_h / HourETRad
      rs_rso = MIN(1.0_wp, MAX(0.0_wp, rs_rso))
      cloud_fac = 1.0_wp + 0.22_wp * (1.0_wp - rs_rso)**2
    END IF
    emiss_raw = emiss_clear * cloud_fac
    sky_emiss = MIN(1.0_wp, MAX(0.0_wp, emiss_raw))
  END BLOCK
  dLW_in    = sky_emiss  * const_SB * (dTs** 4)
  dLW_net   = dLW_in - dLW_out
ELSE
  ! The branch that runs: empirical clear-sky emissivity 9.35E-6 * Ta^2 (Ta in
  ! K), with the net longwave then scaled by the day's relative sunshine
  ! duration to account for cloud.
  sky_emiss = 0.00000935 * (dTa ** 2)
  dLW_in    = sky_emiss  * const_SB * (dTs** 4)
  dLW_net   = dLW_in - dLW_out
  dLW_net   = dLW_net * ( 0.2 + (0.8 * (s_hr_s / Daylength)) )
END IF
IF (meteo_debug) THEN
  print *, 'LWrad = ', dLW_net
END IF
HourSW_net_s = dSW_net
HourLW_net_s = dLW_net
HourNetRad= dSW_net + dLW_net
!------------------------------------------------------------------------------!
END SUBROUTINE meteo
!==============================================================================!