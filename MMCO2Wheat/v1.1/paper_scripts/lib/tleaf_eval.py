#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Observation series and model readers for the energy-balance, soil-water, soil-temperature and
canopy-temperature comparisons (WheatDryFACE and Braunschweig).
"""
import csv
import datetime as dt
import functools
import math
import os
import sys
from collections import OrderedDict, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import face_results as P2

fl, isnan, rows_of = P2.fl, P2.isnan, P2.rows_of
NAN = float("nan")
DECKS = P2.DECKS
deck_path, deck_dir = P2.deck_path, P2.deck_dir
layer_bounds, swc_span_mean, SwcBasis = P2.layer_bounds, P2.swc_span_mean, P2.SwcBasis


def mean(xs):
    xs = [x for x in xs if not isnan(x)]
    return (sum(xs) / len(xs)) if xs else NAN


# ===================================================================== obs readers
@functools.lru_cache(maxsize=None)
def wdf25_obs(deck):
    cv, co2, wat = DECKS[deck]["cultivar"], DECKS[deck]["co2"], DECKS[deck]["water"]
    return tuple(P2.W25.obs_for(cv, co2, wat))


@functools.lru_cache(maxsize=None)
def wdf23_obs(deck):
    cv, co2 = DECKS[deck]["cultivar"], DECKS[deck]["co2"]
    return tuple(P2.W23.obs_for(cv, co2))


SCALE = {"SWC 20cm": 0.01, "SWC 55cm": 0.01,   # per cent -> v/v
         "Soil heat flux": -1.0}              # the file's sign convention (negative by day) flipped to the
                                              # model's (positive into the soil by day)


@functools.lru_cache(maxsize=None)
def obs_series(deck, variable, period=None):
    """[(date, value, sd)] for one deck/variable (WheatDryFACE only), scaled."""
    site = DECKS[deck]["site"]
    rows = wdf25_obs(deck) if site == "WDF25" else (wdf23_obs(deck) if site == "WDF23" else [])
    sc = SCALE.get(variable, 1.0)
    out = []
    for r in rows:
        if r["variable"] != variable:
            continue
        if period is not None and r.get("period") != period:
            continue
        if period is None and r.get("period") is not None:
            continue
        v, sd = r["value"], r["std"]
        out.append((r["date"], (v * sc if not isnan(v) else NAN),
                   (sd * sc if sd is not None and not isnan(sd) else NAN)))
    return sorted(out)


@functools.lru_cache(maxsize=None)
def bra_soil_moisture(deck, variable):
    """[(date, value_v/v, sd_v/v, n_rings)] from the deck's observation file (d0_20 / d20_40)."""
    p = P2.bra_obs_path(deck)
    out = []
    if not os.path.isfile(p):
        return out
    for r in rows_of(p):
        if r["variable"] != variable:
            continue
        v, sd = fl(r["mean"]), fl(r["sd"])
        out.append((r["date"], v / 100.0, (sd / 100.0 if not isnan(sd) else NAN), fl(r["n_rings"])))
    return sorted(out)


@functools.lru_cache(maxsize=None)
def flir_obs_deck(deck):
    """[(date, obs_tleaf_C, obs_sd)] for one WDF25 deck: the thermal-image canopy temperature
    (10:00-14:00 window) from the Tleaf_flir / std columns of the treatment file."""
    cv, co2, wat = DECKS[deck]["cultivar"], DECKS[deck]["co2"], DECKS[deck]["water"]
    p = os.path.join(config.WDF25_V3_DIR, "S3_%s_%s_%s.csv" % (cv, co2, wat))
    out = []
    if not os.path.isfile(p):
        return out
    rows = list(csv.reader(open(p, encoding="utf-8-sig")))
    hdr = rows[1]
    if "Tleaf_flir" not in hdr:
        return out
    j = hdr.index("Tleaf_flir")
    for r in rows[2:]:
        if len(r) > j + 1 and r[j].strip():
            d = dt.datetime.strptime(r[0], "%m/%d/%Y").date().isoformat()
            out.append((d, fl(r[j]), fl(r[j + 1])))
    return sorted(out)


# ===================================================================== model readers
MODELS = {}


def M(pkg, deck):
    if deck not in MODELS:
        MODELS[deck] = P2.DeckModel(pkg, deck)
    return MODELS[deck]


@functools.lru_cache(maxsize=None)
def _swc_basis(pkg, deck):
    return SwcBasis(pkg, deck)


def model_swc(pkg, deck, label):
    """[(day, v/v)] for label in ('SWC20','SWC55','d0_20','d20_40') on the deck's SwcBasis."""
    sb = _swc_basis(pkg, deck)
    out = []
    for lab, z0, z1 in sb.labels:
        if lab != label:
            continue
        for day in sorted(sb.rows):
            v = sb.model_at(day, z0, z1)
            if v is not None:
                out.append((day, v))
    return out


@functools.lru_cache(maxsize=None)
def model_tsoil(pkg, deck, z):
    """Model soil temperature (C) at depth z (m), interpolated between layer mid-points of
    daily_tsoil_profile.csv."""
    bnds = layer_bounds(pkg, deck)
    if not bnds:
        return {}
    mid = [(bnds[i] + bnds[i + 1]) / 2.0 for i in range(10)]
    rows = rows_of(os.path.join(deck_dir(pkg, deck), "daily_tsoil_profile.csv"))
    out = {}
    for r in rows:
        day = int(fl(r["day"]))
        v = [fl(r["Tsoil_L%d" % (i + 1)]) for i in range(10)]
        if z <= mid[0]:
            out[day] = v[0]
        elif z >= mid[-1]:
            out[day] = v[-1]
        else:
            for i in range(9):
                if mid[i] <= z <= mid[i + 1]:
                    f = (z - mid[i]) / (mid[i + 1] - mid[i])
                    out[day] = v[i] + f * (v[i + 1] - v[i])
                    break
    return out


def model_radiation_aero(pkg, deck):
    return {int(fl(r["day"])): r for r in rows_of(os.path.join(deck_dir(pkg, deck), "daily_radiation_aero.csv"))}


@functools.lru_cache(maxsize=None)
def model_eb_daynight(pkg, deck):
    """day -> day-mean and night-mean of Rn_total, LW_down, LW_up_sky, SW_down and
    SW_ref_can + SW_ref_soil from hourly_radiation.csv, and of G_surface_W_m2 from
    hourly_soil_heat.csv; day/night split = SW_down > 0 / == 0."""
    out = defaultdict(lambda: dict(Rn_day=[], Rn_night=[], LWin_day=[], LWin_night=[],
                                   LWout_day=[], LWout_night=[], SWin_day=[], SWin_night=[],
                                   SWout_day=[], SWout_night=[]))
    for r in rows_of(os.path.join(deck_dir(pkg, deck), "hourly_radiation.csv")):
        day = int(fl(r["day"]))
        sw = fl(r["SW_down"])
        tag = "day" if sw > 0 else "night"
        out[day]["Rn_%s" % tag].append(fl(r["Rn_total"]))
        out[day]["LWin_%s" % tag].append(fl(r["LW_down"]))
        out[day]["LWout_%s" % tag].append(fl(r["LW_up_sky"]))
        out[day]["SWin_%s" % tag].append(sw)
        out[day]["SWout_%s" % tag].append(fl(r["SW_ref_can"]) + fl(r["SW_ref_soil"]))
    g = defaultdict(lambda: dict(G_day=[], G_night=[]))
    swp = os.path.join(deck_dir(pkg, deck), "hourly_radiation.csv")
    swrows = rows_of(swp)
    swflag = {}
    for r in swrows:
        swflag[(int(fl(r["day"])), fl(r["solar_hour"]))] = fl(r["SW_down"]) > 0
    for r in rows_of(os.path.join(deck_dir(pkg, deck), "hourly_soil_heat.csv")):
        day = int(fl(r["day"]))
        hh = fl(r["hour"])
        tag = "day" if swflag.get((day, hh), False) else "night"
        g[day]["G_%s" % tag].append(fl(r["G_surface_W_m2"]))
    means = {}
    for day in sorted(set(list(out.keys()) + list(g.keys()))):
        d = {}
        for k, v in out.get(day, {}).items():
            d[k] = mean(v) if v else NAN
        for k, v in g.get(day, {}).items():
            d[k] = mean(v) if v else NAN
        means[day] = d
    return means


@functools.lru_cache(maxsize=None)
def model_albedo_day(pkg, deck):
    """alpha_mix of daily_radiation_aero.csv (one value per day)."""
    rows = model_radiation_aero(pkg, deck)
    return {d: fl(r["alpha_mix"]) for d, r in rows.items()}


# ===================================================================== quantity registry
TOL_SWC = 0.02   # |bias| <= 0.02 v/v

QTY = OrderedDict([
    ("SWC20", dict(sites=["WDF25", "WDF23"], obs_var="SWC 20cm", period=None, mlab="SWC20",
                   unit="v/v", kind="swc", tol=TOL_SWC, basis="contrast=ratio")),
    ("SWC55", dict(sites=["WDF25", "WDF23"], obs_var="SWC 55cm", period=None, mlab="SWC55",
                   unit="v/v", kind="swc", tol=TOL_SWC, basis="contrast=ratio")),
    ("d0_20", dict(sites=["BRAUNSCHWEIG"], obs_var="d0_20", period=None, mlab="d0_20",
                   unit="v/v", kind="bra_swc", tol=TOL_SWC, basis="contrast=ratio")),
    ("d20_40", dict(sites=["BRAUNSCHWEIG"], obs_var="d20_40", period=None, mlab="d20_40",
                    unit="v/v", kind="bra_swc", tol=TOL_SWC, basis="contrast=ratio")),
    ("Tsoil20", dict(sites=["WDF25", "WDF23"], obs_var="Tsoil 20cm", period=None, depth=0.20,
                     unit="C", kind="tsoil", tol=None, basis="contrast=difference")),
    ("Tsoil55", dict(sites=["WDF25", "WDF23"], obs_var="Tsoil55 cm", period=None, depth=0.55,
                     unit="C", kind="tsoil", tol=None, basis="contrast=difference")),
    ("Rn_day", dict(sites=["WDF25", "WDF23"], obs_var="Rn", period="Day", mkey="Rn_day",
                    unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference")),
    ("Rn_night", dict(sites=["WDF25", "WDF23"], obs_var="Rn", period="Night", mkey="Rn_night",
                      unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference")),
    ("LWin_day", dict(sites=["WDF25", "WDF23"], obs_var="LWin", period="Day", mkey="LWin_day",
                      unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference")),
    ("LWin_night", dict(sites=["WDF25", "WDF23"], obs_var="LWin", period="Night", mkey="LWin_night",
                        unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference")),
    ("LWout_day", dict(sites=["WDF25", "WDF23"], obs_var="LWout", period="Day", mkey="LWout_day",
                       unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference")),
    ("LWout_night", dict(sites=["WDF25", "WDF23"], obs_var="LWout", period="Night", mkey="LWout_night",
                         unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference")),
    ("SWin_day", dict(sites=["WDF25", "WDF23"], obs_var="SWin", period="Day", mkey="SWin_day",
                      unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference")),
    ("SWout_day", dict(sites=["WDF25", "WDF23"], obs_var="SWout", period="Day", mkey="SWout_day",
                       unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference")),
    ("albedo_day", dict(sites=["WDF25", "WDF23"], obs_var="albedo", period="Day", kind="albedo",
                        unit="-", tol=None, basis="contrast=difference")),
    ("G_day", dict(sites=["WDF25"], obs_var="Soil heat flux", period="Day", mkey="G_day",
                   unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference",
                   note="Ruta rings 4/6 only; Zahir has no sensor")),
    ("G_night", dict(sites=["WDF25"], obs_var="Soil heat flux", period="Night", mkey="G_night",
                     unit="W/m2", kind="ebdn", tol=None, basis="contrast=difference",
                     note="Ruta rings 4/6 only; Zahir has no sensor")),
    ("Tleaf_FLIR", dict(sites=["WDF25"], kind="flir", unit="C", tol=None,
                        basis="contrast=difference")),
])


def obs_pairs_for(pkg, deck, qkey):
    """[(date, obs, obs_sd, model)] for one deck/quantity, using the registry above."""
    q = QTY[qkey]
    m = M(pkg, deck)
    if q["kind"] == "swc":
        obs = obs_series(deck, q["obs_var"], q["period"])
        mdict = dict(model_swc(pkg, deck, q["mlab"]))
        return [(d, o, s, mdict[m.date2day[d]]) for d, o, s in obs
                if d in m.date2day and m.date2day[d] in mdict]
    if q["kind"] == "bra_swc":
        obs = [(d, v, s) for d, v, s, n in bra_soil_moisture(deck, q["obs_var"])]
        mdict = dict(model_swc(pkg, deck, q["mlab"]))
        return [(d, o, s, mdict[m.date2day[d]]) for d, o, s in obs
                if d in m.date2day and m.date2day[d] in mdict]
    if q["kind"] == "tsoil":
        obs = obs_series(deck, q["obs_var"], q["period"])
        mdict = model_tsoil(pkg, deck, q["depth"])
        return [(d, o, s, mdict[m.date2day[d]]) for d, o, s in obs
                if d in m.date2day and m.date2day[d] in mdict]
    if q["kind"] == "ebdn":
        obs = obs_series(deck, q["obs_var"], q["period"])
        mdict = {day: v.get(q["mkey"]) for day, v in model_eb_daynight(pkg, deck).items()}
        return [(d, o, s, mdict[m.date2day[d]]) for d, o, s in obs
                if d in m.date2day and m.date2day[d] in mdict
                and mdict[m.date2day[d]] is not None and not isnan(mdict[m.date2day[d]])]
    if q["kind"] == "albedo":
        obs = obs_series(deck, "albedo", q["period"])
        mdict = model_albedo_day(pkg, deck)
        return [(d, o, s, mdict[m.date2day[d]]) for d, o, s in obs
                if d in m.date2day and m.date2day[d] in mdict]
    return []
