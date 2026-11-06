#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Observed-simulated pairs for the energy-balance components, soil water, soil temperature,
leaf gas exchange and canopy temperature at every FACE deck that measured the quantity.

Model bases
  energy balance   hourly_radiation.csv / hourly_soil_heat.csv.  The day block is the n_day
                   highest-SW_down hours, n_day recovered per date from the deck's own forcing
                   (daily shortwave total / observed day-block mean), so both sides use one
                   window; the night block is the remaining 24 - n_day hours.  Composite upward
                   SW = SW_ref_can + SW_ref_soil, upward LW = LW_up_sky, net = Rn_total.
  albedo           mean of the hourly (SW_ref_can + SW_ref_soil) / SW_down ratios over the same window.
  soil heat flux   24-h mean (G_24h): the plate is 0.10 m below the surface and the model flux is
                   at the surface; the two are one quantity over a full diurnal cycle.
  Tsurf_rad        brightness temperature from each side's own upward/downward long-wave pair,
                   one emissivity (0.98), the same inversion on both sides.
  soil water       the time-weighted daily-mean layer file at every site, point-interpolated at a
                   point sensor's depth and thickness-weighted over an integrating sensor's span.
  soil temperature interpolation between layer mid-points at the sensor depths.
  gas exchange     flag_leaf_diag.csv, mean and sd across the hourly rows inside the campaign's
                   clock window (10:00-14:00 local time, converted to solar time per date); leaf
                   transpiration = gsw x Ds / Patm (Patm = 101325 Pa).
  canopy temperature  the hourly Tcan_C of hourly_energy_balance.csv over the same window.
"""
import csv
import datetime as _dt
import math
import os
import re
import sys
from collections import OrderedDict, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import face_results as P2
import tleaf_eval as TE

fl, isnan, rows_of = P2.fl, P2.isnan, P2.rows_of
NAN = float("nan")
DECKS = P2.DECKS
MAR_SWC_CSV = os.path.join(config.MARICOPA_OBS_DIR, "MARICOPA_SOIL_MOISTURE_AVG_OVER_REPS.csv")
PATM = 101325.0
TSOIL_MAX_VALID = 45.0         # a soil temperature above this at 20/55 cm under wheat is a logging fault
SB_CONST = 5.670374419e-8      # Stefan-Boltzmann, CODATA 2018
EPS_SURF = 0.98                # broadband surface emissivity used identically on both sides

# Fortran G14.6 exponent collapse: with a three-digit exponent the 'E' is dropped and the field
# reads '0.757898-254'; flg() repairs the exponent instead of returning NaN.
_GFIX = re.compile(r"^([+-]?[0-9.]+)([+-][0-9]{2,3})$")


def flg(x, default=0.0):
    v = fl(x)
    if not isnan(v):
        return v
    m = _GFIX.match((x or "").strip())
    if m:
        try:
            return float(m.group(1) + "E" + m.group(2))
        except Exception:
            return default
    return default


# --------------------------------------------------------------------- quantity registry
# key -> (label, unit, sites, loader-kind)
Q = OrderedDict([
    ("Rn_day",      ("Net radiation, day mean",         "W m-2",   ["WDF25", "WDF23"], "wdf")),
    ("Rn_night",    ("Net radiation, night mean",       "W m-2",   ["WDF25", "WDF23"], "wdf")),
    ("LWin_day",    ("Longwave in, day mean",           "W m-2",   ["WDF25", "WDF23"], "wdf")),
    ("LWin_night",  ("Longwave in, night mean",         "W m-2",   ["WDF25", "WDF23"], "wdf")),
    ("LWout_day",   ("Longwave out, day mean",          "W m-2",   ["WDF25", "WDF23"], "wdf")),
    ("LWout_night", ("Longwave out, night mean",        "W m-2",   ["WDF25", "WDF23"], "wdf")),
    ("SWin_day",    ("Shortwave in, day mean",          "W m-2",   ["WDF25", "WDF23"], "wdf")),
    ("SWout_day",   ("Shortwave out (reflected), day",  "W m-2",   ["WDF25", "WDF23"], "wdf")),
    ("albedo_day",  ("Albedo, day",                     "-",       ["WDF25", "WDF23"], "wdf")),
    ("G_day",       ("Soil heat flux G, day mean",      "W m-2",   ["WDF25"],          "wdf")),
    ("G_night",     ("Soil heat flux G, night mean",    "W m-2",   ["WDF25"],          "wdf")),
    ("SWC20",       ("Soil water content 20 cm",        "m3 m-3",  ["WDF25", "WDF23"], "wdf")),
    ("SWC55",       ("Soil water content 55 cm",        "m3 m-3",  ["WDF25", "WDF23"], "wdf")),
    ("Tsoil20",     ("Soil temperature 20 cm",          "deg C",   ["WDF25", "WDF23"], "wdf")),
    ("Tsoil55",     ("Soil temperature 55 cm",          "deg C",   ["WDF25", "WDF23"], "wdf")),
    ("d0_20",       ("Soil water content 0-20 cm",      "m3 m-3",  ["BRAUNSCHWEIG"],   "wdf")),
    ("d20_40",      ("Soil water content 20-40 cm",     "m3 m-3",  ["BRAUNSCHWEIG"],   "wdf")),
    ("MAR_SWC_0_30",     ("Soil water content 0-30 cm",      "m3 m-3",  ["MARICOPA"],       "mar")),
    ("MAR_SWC_30_50",     ("Soil water content 30-50 cm",     "m3 m-3",  ["MARICOPA"],       "mar")),
    ("Tleaf_FLIR",  ("Canopy T: thermal images vs T_can", "deg C", ["WDF25"],          "wdf")),
    ("Tsurf_rad",   ("Radiometric surface temperature (from outgoing LW), day", "deg C", ["WDF25", "WDF23"], "trad")),
    ("G_24h",       ("Soil heat flux G, 24-h mean",     "W m-2",   ["WDF25"],          "g24")),
    ("An",          ("Flag-leaf net assimilation An",   "umol m-2 s-1", ["WDF25", "WDF23"], "gx")),
    ("gs",          ("Stomatal conductance gsw",        "mol m-2 s-1",  ["WDF25", "WDF23"], "gx")),
    ("Tr",          ("Leaf transpiration Tr",           "mmol m-2 s-1", ["WDF25", "WDF23"], "gx")),
])

MAR_SPAN = {"MAR_SWC_0_30": (0.0, 0.30, "D30"), "MAR_SWC_30_50": (0.30, 0.50, "D50")}
GX_OBS = {"An": "An", "gs": "gs", "Tr": "Tr"}
GX_MODEL = {"An": "An_umol_m2_s", "gs": "gsw_mol_m2_s"}


def PKG():
    return P2.PKG


def model(deck):
    return TE.M(PKG(), deck)


_CACHE = {}

# Local clock window of the WheatDryFACE campaigns and its conversion to the deck's solar time.
TZ = {"WDF25": (2, "IL"), "WDF23": (2, "IL")}   # UTC offset of standard time, daylight-saving rule
CLOCK_WIN = (10.0, 14.0)                        # local clock hours of the campaign window
HALF_STEP = 0.5                                 # flag_leaf_diag is one point sample per whole solar hour


def _il_dst(d):
    """True if the Israel summer hour (UTC+3) is in force on date d: from the Friday before the
    last Sunday of March to the last Sunday of October."""
    last_sun_mar = max(day for day in range(25, 32)
                       if _dt.date(d.year, 3, day).weekday() == 6)
    start = _dt.date(d.year, 3, last_sun_mar) - _dt.timedelta(days=2)
    last_sun_oct = max(day for day in range(25, 32)
                       if _dt.date(d.year, 10, day).weekday() == 6)
    return start <= d <= _dt.date(d.year, 10, last_sun_oct)


def _deck_longitude(deck):
    """Longitude (deg east) from the deck's own driver_main.txt record 4."""
    key = ("lon", deck)
    if key in _CACHE:
        return _CACHE[key]
    v = NAN
    p = os.path.join(TE.deck_path(PKG(), deck), "input", "driver_main.txt")
    if os.path.isfile(p):
        lines = open(p).read().splitlines()
        if len(lines) >= 4:
            v = fl(lines[3].split("!")[0].strip())
    _CACHE[key] = v
    return v


def _solar_window(deck, date_iso, clock=CLOCK_WIN):
    """(solar_lo, solar_hi, source) for a local clock window on one date at one deck."""
    site = DECKS[deck]["site"]
    if site not in TZ:
        return clock[0], clock[1], "solar window (no time zone declared for this site)"
    utc_off, rule = TZ[site]
    lon = _deck_longitude(deck)
    if isnan(lon):
        return clock[0], clock[1], "solar window (deck longitude unreadable)"
    try:
        d = _dt.date(*[int(x) for x in date_iso.split("-")])
    except Exception:
        return clock[0], clock[1], "solar window (unparsable date)"
    n = d.timetuple().tm_yday
    B = math.radians(360.0 * (n - 81) / 364.0)
    eot = 9.87 * math.sin(2 * B) - 7.53 * math.cos(B) - 1.5 * math.sin(B)      # minutes
    dst = 1.0 if (rule == "IL" and _il_dst(d)) else 0.0
    shift = (4.0 * (lon - 15.0 * utc_off) + eot) / 60.0 - dst                  # hours
    return (clock[0] + shift, clock[1] + shift,
            "clock %04.1f-%04.1f %s converted to solar (lon %.4f E, L_std %d E, EoT %+.1f min, "
            "DST %s)" % (clock[0], clock[1], "IDT" if dst else "local standard", lon,
                         int(15 * utc_off), eot, "on" if dst else "off"))


def _flag_leaf_window(deck):
    """day -> {An, gs, Tleaf, Tr: mean over the window hours of flag_leaf_diag.csv; <key>_sd:
    sample sd across those hours; n_hours; window}.  A row is taken when its solar hour lies
    within the window widened by half a sampling step either side."""
    key = ("flw", deck)
    if key in _CACHE:
        return _CACHE[key]
    rows = rows_of(os.path.join(TE.deck_dir(PKG(), deck), "flag_leaf_diag.csv"))
    mm = model(deck)
    acc = defaultdict(lambda: defaultdict(list))
    wsrc = {}
    for r in rows:
        h = fl(r["hour"])
        d = int(fl(r["day"]))
        ds = mm.day2date.get(d)
        if ds is None:
            continue
        if d not in wsrc:
            wsrc[d] = _solar_window(deck, ds)
        s0, s1, src = wsrc[d]
        if s0 - HALF_STEP < h < s1 + HALF_STEP:
            an, gs, tl, dsp = (fl(r["An_umol_m2_s"]), fl(r["gsw_mol_m2_s"]),
                               fl(r["Tleaf_C"]), fl(r["Ds_Pa"]))
            acc[d]["An"].append(an)
            acc[d]["gs"].append(gs)
            acc[d]["Tleaf"].append(tl)
            acc[d]["Tr"].append(gs * max(dsp, 0.0) / PATM * 1000.0)
    out = {}
    for d, dd in acc.items():
        row = {}
        for k, v in dd.items():
            nn = len(v)
            mu = (sum(v) / nn) if nn else NAN
            row[k] = mu
            row[k + "_sd"] = (math.sqrt(sum((x - mu) ** 2 for x in v) / (nn - 1))
                              if nn > 1 else NAN)
        row["n_hours"] = len(dd.get("An", []))
        row["window"] = wsrc.get(d, (NAN, NAN, ""))
        out[d] = row
    _CACHE[key] = out
    return out


def _obs_day_hours(deck):
    """{date: number of whole hours in the observation file's own Day block}, recovered from
    the deck's own forcing: n_day = Rad_MJ * 1e6 / 3600 / SWin_day (the daily shortwave total
    of input_daily_meteo.csv is the integral of the same pyranometer whose Day-block mean the
    observation file prints)."""
    key = ("odh", deck)
    if key in _CACHE:
        return _CACHE[key]
    out = {}
    met = os.path.join(TE.deck_path(PKG(), deck), "input", "input_daily_meteo.csv")
    if os.path.isfile(met):
        rad = {}
        for r in rows_of(met):
            try:
                ds = _dt.date(int(fl(r["year"])), int(fl(r["month"])),
                              int(fl(r["day"]))).isoformat()
            except Exception:
                continue
            for k in r:
                if k.strip().startswith("Rad"):
                    v = fl(r[k])
                    if not isnan(v):
                        rad[ds] = v
        for d, o, s in TE.obs_series(deck, "SWin", "Day"):
            if d in rad and o and o > 0:
                nh = int(round(rad[d] * 1.0e6 / 3600.0 / o))
                if 1 <= nh <= 24:
                    out[d] = nh
    _CACHE[key] = out
    return out


def model_eb_obswindow(deck):
    """day -> dict of day/night means (same keys as tleaf_eval.model_eb_daynight), the day block
    being the n_day highest-SW_down hours of the model day (n_day from _obs_day_hours) and the
    night block the remaining hours; falls back to the SW_down > 0 split where no window is recovered."""
    key = ("ebow", deck)
    if key in _CACHE:
        return _CACHE[key]
    m = model(deck)
    nday = _obs_day_hours(deck)
    hrs = defaultdict(list)
    for r in rows_of(os.path.join(TE.deck_dir(PKG(), deck), "hourly_radiation.csv")):
        hrs[int(fl(r["day"]))].append((
            flg(r["SW_down"]), fl(r["solar_hour"]), flg(r["SW_down"]),
            flg(r["SW_ref_can"]) + flg(r["SW_ref_soil"]), flg(r["LW_down"]),
            flg(r["LW_up_sky"]), flg(r["Rn_total"])))
    ghr = defaultdict(dict)
    gp = os.path.join(TE.deck_dir(PKG(), deck), "hourly_soil_heat.csv")
    if os.path.isfile(gp):
        for r in rows_of(gp):
            ghr[int(fl(r["day"]))][round(fl(r["hour"]))] = flg(r["G_surface_W_m2"])
    out = {}
    for day, rr in hrs.items():
        ds = m.day2date.get(day)
        nh = nday.get(ds) if ds else None
        rs = sorted(rr, key=lambda x: -x[0])
        if nh is None:
            dayrows = [x for x in rr if x[0] > 0.0]
            nightrows = [x for x in rr if x[0] <= 0.0]
        else:
            dayrows, nightrows = rs[:nh], rs[nh:]

        def mn(rows, i):
            v = [x[i] for x in rows]
            return (sum(v) / len(v)) if v else NAN
        d = dict(SWin_day=mn(dayrows, 2), SWin_night=mn(nightrows, 2),
                 SWout_day=mn(dayrows, 3), SWout_night=mn(nightrows, 3),
                 LWin_day=mn(dayrows, 4), LWin_night=mn(nightrows, 4),
                 LWout_day=mn(dayrows, 5), LWout_night=mn(nightrows, 5),
                 Rn_day=mn(dayrows, 6), Rn_night=mn(nightrows, 6),
                 n_day_hours=(nh if nh is not None else len(dayrows)),
                 window_source=("obs-recovered" if nh is not None else "model SW_down>0"))
        gd = ghr.get(day, {})
        if gd:
            dh = [round(x[1]) for x in dayrows]
            nh2 = [round(x[1]) for x in nightrows]
            gday = [gd[h] for h in dh if h in gd]
            gnight = [gd[h] for h in nh2 if h in gd]
            gall = list(gd.values())
            d["G_day"] = (sum(gday) / len(gday)) if gday else NAN
            d["G_night"] = (sum(gnight) / len(gnight)) if gnight else NAN
            d["G_24h"] = (sum(gall) / len(gall)) if gall else NAN
        out[day] = d
    _CACHE[key] = out
    return out


def model_albedo_composite(deck):
    """day -> the model's composite daytime albedo: the mean of the hourly
    (SW_ref_can + SW_ref_soil) / SW_down ratios over the same day window as the energy-balance
    means, restricted to hours with SW_down > 0."""
    key = ("albcomp", deck)
    if key in _CACHE:
        return _CACHE[key]
    mm = model(deck)
    nday = _obs_day_hours(deck)
    acc = defaultdict(list)
    for r in rows_of(os.path.join(TE.deck_dir(PKG(), deck), "hourly_radiation.csv")):
        sw = flg(r["SW_down"])
        if sw > 0.0:
            acc[int(fl(r["day"]))].append(
                (sw, (flg(r["SW_ref_can"]) + flg(r["SW_ref_soil"])) / sw))
    out = {}
    for d, v in acc.items():
        if not v:
            continue
        ds = mm.day2date.get(d)
        nh = nday.get(ds) if ds else None
        vv = sorted(v, key=lambda x: -x[0])
        if nh is not None:
            vv = vv[:nh]
        out[d] = sum(x[1] for x in vv) / len(vv)
    _CACHE[key] = out
    return out


def _brightness_T(lw_up, lw_down):
    """Radiometric (brightness) surface temperature in C from an upward and a downward long-wave
    flux with one emissivity: T = ((LWup - (1-eps) LWdown)/(eps sigma))^0.25."""
    if isnan(lw_up) or isnan(lw_down):
        return NAN
    v = (lw_up - (1.0 - EPS_SURF) * lw_down) / (EPS_SURF * SB_CONST)
    if v <= 0.0:
        return NAN
    return v ** 0.25 - 273.15


def model_Tsurf_rad(deck):
    """day -> model composite radiometric surface temperature (C) over the day block of the
    energy-balance means, from LW_up_sky and LW_down."""
    key = ("tsr", deck)
    if key in _CACHE:
        return _CACHE[key]
    dn = model_eb_obswindow(deck)
    out = {d: _brightness_T(v.get("LWout_day", NAN), v.get("LWin_day", NAN))
           for d, v in dn.items()}
    _CACHE[key] = {d: v for d, v in out.items() if not isnan(v)}
    return _CACHE[key]


def model_Tcan_window(deck, col="Tcan_C"):
    """day -> (mean, sd, n) of the model's hourly canopy temperature over the campaign's clock
    window (converted to solar time per date), from hourly_energy_balance.csv."""
    key = ("tcanw", deck, col)
    if key in _CACHE:
        return _CACHE[key]
    mm = model(deck)
    acc = defaultdict(list)
    wsrc = {}
    p = os.path.join(TE.deck_dir(PKG(), deck), "hourly_energy_balance.csv")
    if os.path.isfile(p):
        for r in rows_of(p):
            if col not in r:
                continue
            d = int(fl(r["day"]))
            ds = mm.day2date.get(d)
            if ds is None:
                continue
            if d not in wsrc:
                wsrc[d] = _solar_window(deck, ds)
            s0, s1, _ = wsrc[d]
            h = fl(r["solar_hour"])
            if s0 - HALF_STEP < h < s1 + HALF_STEP:
                acc[d].append(flg(r[col], NAN))
    out = {}
    for d, v in acc.items():
        v = [x for x in v if not isnan(x)]
        if not v:
            continue
        mu = sum(v) / len(v)
        sd = (math.sqrt(sum((x - mu) ** 2 for x in v) / (len(v) - 1))
              if len(v) > 1 else NAN)
        out[d] = (mu, sd, len(v))
    _CACHE[key] = out
    return out


def _mar_swc_obs():
    """{TRNO: [(date, D30, D50)]} volumetric soil moisture of the Maricopa archive, treatment means."""
    key = ("marobs",)
    if key in _CACHE:
        return _CACHE[key]
    out = defaultdict(list)
    if os.path.isfile(MAR_SWC_CSV):
        lines = [l for l in open(MAR_SWC_CSV) if not l.startswith("#")]
        for r in csv.DictReader(lines):
            try:
                tr = int(float(r["TRNO"]))
            except Exception:
                continue
            out[tr].append((r["Date"], fl(r["D30"]), fl(r["D50"])))
    _CACHE[key] = out
    return out


def pairs(deck, qkey):
    """[(date, obs, obs_sd, model, model_sd)] on the dates the quantity was observed at this deck."""
    key = ("pairs", deck, qkey)
    if key in _CACHE:
        return _CACHE[key]
    site = DECKS[deck]["site"]
    lab, unit, sites, kind = Q[qkey]
    out = []
    if site in sites:
        if kind == "wdf":
            if qkey == "albedo_day":
                mdl = model_albedo_composite(deck)
                mm = model(deck)
                out = [(d, o, s, mdl[mm.date2day[d]])
                       for d, o, s in TE.obs_series(deck, "albedo", "Day")
                       if d in mm.date2day and mm.date2day[d] in mdl]
            elif qkey == "Tleaf_FLIR":
                mdl = model_Tcan_window(deck, "Tcan_C")
                mm = model(deck)
                out = [(d, o, s, mdl[mm.date2day[d]][0], mdl[mm.date2day[d]][1])
                       for d, o, s in TE.flir_obs_deck(deck)
                       if d in mm.date2day and mm.date2day[d] in mdl]
            elif TE.QTY.get(qkey, {}).get("kind") == "ebdn":
                mdl = model_eb_obswindow(deck)
                mm = model(deck)
                mk = TE.QTY[qkey]["mkey"]
                out = []
                for d, o, s in TE.obs_series(deck, TE.QTY[qkey]["obs_var"],
                                              TE.QTY[qkey]["period"]):
                    day = mm.date2day.get(d)
                    if day is None or day not in mdl:
                        continue
                    v = mdl[day].get(mk)
                    if v is not None and not isnan(v):
                        out.append((d, o, s, v))
            else:
                out = list(TE.obs_pairs_for(PKG(), deck, qkey))
        elif kind == "g24":
            # observed 24-h mean rebuilt from its own Day/Night blocks and the recovered window
            # length: (Day * n_day + Night * (24 - n_day)) / 24; its sd is the hour-weighted
            # arithmetic mean of the two block sds (the two blocks come from the same rings)
            mdl = model_eb_obswindow(deck)
            mm = model(deck)
            nday = _obs_day_hours(deck)
            od = {d: (o, s) for d, o, s in TE.obs_series(deck, "Soil heat flux", "Day")}
            on = {d: (o, s) for d, o, s in TE.obs_series(deck, "Soil heat flux", "Night")}
            for d in sorted(set(od) & set(on) & set(nday)):
                day = mm.date2day.get(d)
                if day is None or day not in mdl:
                    continue
                nh = nday[d]
                o = (od[d][0] * nh + on[d][0] * (24 - nh)) / 24.0
                s = ((abs(od[d][1]) * nh + abs(on[d][1]) * (24 - nh)) / 24.0
                     if not (isnan(od[d][1]) or isnan(on[d][1])) else NAN)
                v = mdl[day].get("G_24h")
                if v is not None and not isnan(v):
                    out.append((d, o, s, v))
        elif kind == "trad":
            mdl = model_Tsurf_rad(deck)
            mm = model(deck)
            lwi = {d: o for d, o, s in TE.obs_series(deck, "LWin", "Day")}
            for d, o, s in TE.obs_series(deck, "LWout", "Day"):
                day = mm.date2day.get(d)
                if day is None or day not in mdl or d not in lwi:
                    continue
                ot = _brightness_T(o, lwi[d])
                if not isnan(ot):
                    out.append((d, ot, NAN, mdl[day]))
        elif kind == "gx":
            m = model(deck)
            win = _flag_leaf_window(deck)
            for d, o, s in TE.obs_series(deck, GX_OBS[qkey]):
                day = m.date2day.get(d)
                if day is not None and day in win and not isnan(win[day][qkey]):
                    out.append((d, o, s, win[day][qkey],
                                win[day].get(qkey + "_sd", NAN)))
        elif kind == "mar":
            m = model(deck)
            z0, z1, col = MAR_SPAN[qkey]
            bnds = P2.layer_bounds(PKG(), deck)
            mrows = {int(fl(r["day"])): r for r in
                     rows_of(os.path.join(TE.deck_dir(PKG(), deck),
                                          "daily_swc_layers_mean.csv"))}
            tr = int(deck.split("_")[1])
            for d, v30, v50 in _mar_swc_obs().get(tr, []):
                o = v30 if col == "D30" else v50
                day = m.date2day.get(d)
                if day is None or isnan(o):
                    continue
                mv = P2.swc_span_mean(mrows, bnds, day, z0, z1)
                if mv is not None:
                    out.append((d, o, NAN, mv))
    out = [p for p in out if not isnan(p[1]) and p[3] is not None and not isnan(p[3])]
    # uniform 5-tuple: (date, obs, obs_sd, model, model_sd)
    out = [p if len(p) == 5 else (p[0], p[1], p[2], p[3], NAN) for p in out]
    if qkey.startswith("Tsoil"):
        out = [p for p in out if p[1] <= TSOIL_MAX_VALID]
    out.sort()
    _CACHE[key] = out
    return out


def model_course(deck, qkey):
    """[(date, model)] over the model's whole season."""
    lab, unit, sites, kind = Q[qkey]
    m = model(deck)
    if kind == "gx":
        win = _flag_leaf_window(deck)
        return [(m.day2date[d], win[d][qkey]) for d in sorted(win) if d in m.day2date]
    if qkey == "albedo_day":
        dd = model_albedo_composite(deck)
        return [(m.day2date[d], v) for d, v in sorted(dd.items()) if d in m.day2date]
    if qkey == "Tsurf_rad":
        dd = model_Tsurf_rad(deck)
        return [(m.day2date[d], v) for d, v in sorted(dd.items()) if d in m.day2date]
    if qkey == "Tleaf_FLIR":
        dd = model_Tcan_window(deck, "Tcan_C")
        return [(m.day2date[d], v[0]) for d, v in sorted(dd.items()) if d in m.day2date]
    if kind == "g24":
        dd = model_eb_obswindow(deck)
        return [(m.day2date[d], v["G_24h"]) for d, v in sorted(dd.items())
                if d in m.day2date and v.get("G_24h") is not None and not isnan(v["G_24h"])]
    if kind == "wdf" and TE.QTY[qkey]["kind"] == "ebdn":
        mk = TE.QTY[qkey]["mkey"]
        dd = model_eb_obswindow(deck)
        return [(m.day2date[d], v[mk]) for d, v in sorted(dd.items())
                if d in m.day2date and v.get(mk) is not None and not isnan(v[mk])]
    return []
