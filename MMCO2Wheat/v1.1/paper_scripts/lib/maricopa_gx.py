#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Maricopa FACE leaf gas exchange (Wall et al., 2000, 2006; archive workbook 'Physiology Gas Exchange Water Relations',
G. W. Wall) against the model's upper (flag) leaf over the same clock window, 10:00-14:00 MST, natural light on both sides."""
import datetime as dt
import math
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import ebwb_gx as EB
import maricopa_archive as MA

CLOCK = (10.0, 14.0)
LON, LSTD = -111.98, -105.0          # Maricopa Agricultural Center; MST (UTC-7, no daylight saving in Arizona)
_C = {}


def _obs():
    if "obs" in _C:
        return _C["obs"]
    d = MA.gas_exchange()
    d.columns = [c.strip() for c in d.columns]
    for c in ("A", "gs", "DH"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[(d.DH >= CLOCK[0]) & (d.DH <= CLOCK[1])].copy()
    d["date"] = pd.to_datetime(d.Date).dt.date.astype(str)
    _C["obs"] = d
    return d


def _solar(date_iso):
    n = dt.date(*map(int, date_iso.split("-"))).timetuple().tm_yday
    B = math.radians(360.0 * (n - 81) / 364.0)
    eot = 9.87 * math.sin(2 * B) - 7.53 * math.cos(B) - 1.5 * math.sin(B)
    sh = (4.0 * (LON - LSTD) + eot) / 60.0
    return CLOCK[0] + sh, CLOCK[1] + sh


def _model(deck):
    key = ("m", deck)
    if key in _C:
        return _C[key]
    mm = EB.model(deck)
    f = os.path.join(EB.TE.deck_dir(EB.PKG(), deck), "flag_leaf_diag.csv")
    m = pd.read_csv(f, skipinitialspace=True); m.columns = [c.strip() for c in m.columns]
    m["date"] = m.day.astype(int).map(mm.day2date)
    _C[key] = m
    return m


def pairs(deck, q):
    """[(date, obs mean over replicate plots, obs sd, model window mean)] for q in An, gs."""
    col = {"An": "A", "gs": "gs"}[q]; mcol = {"An": "An_umol_m2_s", "gs": "gsw_mol_m2_s"}[q]
    trno = int(deck.split("_")[1]); o = _obs(); o = o[o.TRNO == trno]
    if not len(o):
        return []
    m = _model(deck); out = []
    for d, g in o.groupby("date"):
        reps = g.groupby("Rep")[col].mean().dropna()
        if not len(reps):
            continue
        s0, s1 = _solar(d)
        q_ = m[(m.date == d) & (m.hour > s0 - 0.5) & (m.hour < s1 + 0.5)]
        if not len(q_):
            continue
        out.append((d, float(reps.mean()), float(reps.std(ddof=1)) if len(reps) > 1 else float("nan"), float(q_[mcol].mean())))
    return out
