#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Maricopa hourly energy-balance / canopy-temperature pairs: observations (maricopa_eb_hourly.csv, extracted from the
Kimball et al. archive workbooks by extract_maricopa_eb.py) against the model's hourly_energy_balance.csv of each
Maricopa deck, matched on calendar date and hour (hour ending)."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

OBS = config.MARICOPA_EB_HOURLY
TRNO2DECK = {906: "MARICOPA_906", 908: "MARICOPA_908", 910: "MARICOPA_910", 912: "MARICOPA_912",
             913: "MARICOPA_913", 914: "MARICOPA_914", 915: "MARICOPA_915", 916: "MARICOPA_916"}
_cache = {}


def obs_hourly():
    if "obs" not in _cache:
        h = pd.read_csv(OBS, parse_dates=["date"])
        _cache["obs"] = h
    return _cache["obs"]


def model_hourly(pkg, deck):
    key = ("m", pkg, deck)
    if key in _cache:
        return _cache[key]
    ph = pd.read_csv(os.path.join(pkg, deck, "output", "daily_phenology.csv"))
    ph["date"] = pd.to_datetime(dict(year=ph.year.astype(int), month=ph.month.astype(int), day=ph.dom.astype(int)))
    d2d = dict(zip(ph.day.astype(int), ph.date))
    m = pd.read_csv(os.path.join(pkg, deck, "output", "hourly_energy_balance.csv"), skipinitialspace=True)
    m.columns = [c.strip() for c in m.columns]
    m["day"] = m["day"].astype(int)
    m["date"] = m["day"].map(d2d)
    # model rows are solar hours 0..23 (hour beginning); the archive labels the hour ending -> model hour + 1
    m["hour"] = (m["solar_hour"].round().astype(int) + 1)
    m = m[["date", "hour", "Tair_C", "Tcan_C", "Rn_total_W_m2", "H_total_W_m2", "LEc_W_m2", "LEs_W_m2", "G_W_m2"]].copy()
    m["LE"] = m["LEc_W_m2"] + m["LEs_W_m2"]
    _cache[key] = m
    return m


def pairs(pkg, trno):
    """Merged hourly frame for one treatment: obs tcan/tair/rn/h/le/g and model columns, same date and hour."""
    key = ("p", pkg, trno)
    if key in _cache:
        return _cache[key]
    o = obs_hourly(); o = o[o.trno == trno]
    m = model_hourly(pkg, TRNO2DECK[trno])
    p = o.merge(m, on=["date", "hour"], how="inner", suffixes=("", "_mod"))
    # plausibility screen on the archive's Bowen-ratio fluxes (single-hour failures of the ratio method): drop, do not patch
    for col, lo, hi in (("rn", -200, 1000), ("le", -150, 900), ("h", -300, 600), ("g", -200, 300), ("tcan", -10, 60)):
        p.loc[(p[col] < lo) | (p[col] > hi), col] = np.nan
    _cache[key] = p
    return p


def daily_means(pkg, trno, var_obs, var_mod, hours=(7, 18), night=False):
    """[(date, obs_mean, obs_sd_within, model_mean)] over the day (hours inclusive, hour-ending) or night (complement)."""
    p = pairs(pkg, trno).dropna(subset=[var_obs, var_mod])
    sel = p.hour.between(hours[0], hours[1])
    p = p[~sel] if night else p[sel]
    g = p.groupby("date")
    out = []
    for d, q in g:
        if len(q) < (6 if not night else 6):
            continue
        out.append((d, float(q[var_obs].mean()), float(q[var_obs].std(ddof=1)), float(q[var_mod].mean())))
    return out


def midday_dT(pkg, trno, hours=(11, 14)):
    """[(date, obs dT, obs sd across hours, model dT, obs Tcan, model Tcan, obs Tair, model Tair)], canopy minus air,
    midday window (hour ending 11..14)."""
    p = pairs(pkg, trno).dropna(subset=["tcan", "tair", "Tcan_C", "Tair_C"])
    p = p[p.hour.between(*hours)]
    out = []
    for d, q in p.groupby("date"):
        if len(q) < 3:
            continue
        out.append((d, float((q.tcan - q.tair).mean()), float((q.tcan - q.tair).std(ddof=1)), float((q.Tcan_C - q.Tair_C).mean()),
                    float(q.tcan.mean()), float(q.Tcan_C.mean()), float(q.tair.mean()), float(q.Tair_C.mean())))
    return out
