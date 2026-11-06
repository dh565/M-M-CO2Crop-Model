#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Supplement S2, Fig. S8: WheatDryFACE 2024/25 LI-600 porometer record (stomatal conductance, transpiration, leaf
temperature) against the LI-6800 campaigns and the model, same dates, same 10:00-14:00 window.
Source: the WheatDryFACE 2024/25 treatment files (config.WDF25_PORO_DIR), columns gs_poro (mol m-2 s-1), Tr_poro (mmol m-2 s-1),
Tleaf_poro (C), DTleaf_air_poro (K), each with its std over the ring means."""
import datetime as dt
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import paperfigs as PF
from paperfigs import C_A, C_E, C_OBS, C_SIM, GRAY, W2, P2, EB, isnan, NAN, letter, save
import paperfigs2 as PF2

PORO = config.WDF25_PORO_DIR
NAMES = {"gs": "gs_poro (mol/m2/s)", "Tr": "Tr_poro (mmol/m2/s)", "Tleaf": "Tleaf_poro", "dT": "DTleaf_air_poro"}
DECKS8 = ["WDF25_%s_%s_%s" % (cv, co2, w) for cv in ("Ruta", "Zahir") for co2 in ("aCO2", "eCO2") for w in ("WW", "WL")]


def read_poro(deck):
    """{date: {q: (mean, sd)}} for one treatment file."""
    cv, co2, w = deck.split("_")[1:4]
    f = os.path.join(PORO, "S3_%s_%s_%s.csv" % (cv, co2, w))
    raw = pd.read_csv(f, header=None, dtype=str)
    hdr = raw.iloc[1].tolist()
    COLS = {q: (name, hdr.index(name)) for q, name in NAMES.items()}   # each value column is followed by its std column
    out = {}
    for _, r in raw.iloc[2:].iterrows():
        v = r[COLS["gs"][1]]
        if not isinstance(v, str) or not v.strip(): continue
        d = dt.datetime.strptime(r[0], "%m/%d/%Y").date().isoformat()
        out[d] = {q: (float(r[col]), float(r[col + 1]) if isinstance(r[col + 1], str) and r[col + 1].strip() else NAN) for q, (name, col) in COLS.items()}
    return out


def collect():
    """rows: deck, date, q, obs_poro, sd_poro, obs_li6800 (gs/Tr), model (same clock window)."""
    rows = []
    for dk in DECKS8:
        po = read_poro(dk)
        li = {q: {p[0]: (p[1], p[2]) for p in EB.pairs(dk, q)} for q in ("gs", "Tr")}
        win = EB._flag_leaf_window(dk); mm = EB.model(dk)
        fl = {p[0]: (p[1], p[3]) for p in EB.pairs(dk, "Tleaf_FLIR")}     # FLIR canopy T (obs, model committed Tcan)
        for d, rec in po.items():
            day = mm.date2day.get(d)
            mgs = win[day]["gs"] if day in win else NAN
            mtr = win[day]["Tr"] if (day in win and "Tr" in win[day]) else NAN
            dTm, Tcm, Tam = PF2.model_window(dk, d)
            rows.append(dict(deck=dk, date=d, cv=dk.split("_")[1], co2=dk.split("_")[2], water=dk.split("_")[3],
                             gs_poro=rec["gs"][0], gs_poro_sd=rec["gs"][1], gs_li=li["gs"].get(d, (NAN, NAN))[0], gs_mod=mgs,
                             Tr_poro=rec["Tr"][0], Tr_poro_sd=rec["Tr"][1], Tr_li=li["Tr"].get(d, (NAN, NAN))[0], Tr_mod=mtr,
                             dT_poro=rec["dT"][0], dT_poro_sd=rec["dT"][1], dT_mod=dTm, Tleaf_poro=rec["Tleaf"][0], Tleaf_flir=fl.get(d, (NAN, NAN))[0], Tcan_mod=Tcm))
    return pd.DataFrame(rows)


def figS8(out):
    df = collect()
    df.to_csv(os.path.join(out, "FigS8_porometer_pairs.csv"), index=False)
    fig, axes = plt.subplots(2, 2, figsize=(W2, 0.78 * W2))
    axes = axes.flat
    # (a-c) by date: porometer (filled circles), LI-6800 (open circles), model (open squares, dashed), ambient / elevated
    spec = [("gs", "$g_{sw}$ (mol m$^{-2}$ s$^{-1}$)"), ("Tr", "$T_r$ (mmol m$^{-2}$ s$^{-1}$)"), ("dT", "Leaf − air temperature (K)")]
    for k, (q, lab) in enumerate(spec):
        ax = axes[k]
        for co2, c in (("aCO2", C_A), ("eCO2", C_E)):
            g = df[df.co2 == co2].groupby("date")
            x = [dt.date.fromisoformat(d) for d in g.groups]
            mp = g[q + "_poro"].mean().values; sp = g[q + "_poro"].std(ddof=1).values
            ax.errorbar(x, mp, yerr=sp, fmt="o", color=c, mfc=c, ms=3.4, elinewidth=0.6, capsize=0, zorder=4)
            if q != "dT":
                ml = g[q + "_li"].mean().values
                ax.plot(x, ml, "o", color=c, mfc="white", ms=3.4, mew=0.8, zorder=3)
            mo = g[q + "_mod"].mean().values
            ax.plot(x, mo, "s--", color=c, mfc="white", ms=3.2, lw=0.9, zorder=2)
        ax.set_ylabel(lab); ax.xaxis.set_major_locator(mdates.DayLocator(interval=7)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax.tick_params(axis="x", labelsize=6.5)
        if q == "dT": ax.axhline(0, color=GRAY, lw=0.6)
        letter(ax, k)
    # (d) CO2 effect on gsw by date: porometer, LI-6800, model (mean over the 4 genotype x water pairs, bar = range)
    ax = axes[3]
    dates = sorted(df.date.unique()); x = np.arange(len(dates)); w = 0.26
    def eff(col):
        m, lo, hi = [], [], []
        for d in dates:
            s = df[df.date == d]; r = []
            for cv in ("Ruta", "Zahir"):
                for wt in ("WW", "WL"):
                    e = s[(s.cv == cv) & (s.water == wt) & (s.co2 == "eCO2")][col]; a = s[(s.cv == cv) & (s.water == wt) & (s.co2 == "aCO2")][col]
                    if len(e) and len(a) and not isnan(float(a.iloc[0])) and float(a.iloc[0]) > 0 and not isnan(float(e.iloc[0])):
                        r.append((float(e.iloc[0]) / float(a.iloc[0]) - 1) * 100)
            m.append(np.mean(r) if r else NAN); lo.append(min(r) if r else NAN); hi.append(max(r) if r else NAN)
        return np.array(m), np.array(lo), np.array(hi)
    for i, (col, c, lab) in enumerate((("gs_poro", C_OBS, "porometer (LI-600)"), ("gs_li", "#808080", "LI-6800"), ("gs_mod", C_SIM, "model"))):
        m, lo, hi = eff(col)
        ax.bar(x + (i - 1) * w, m, w, color=c, label=lab)
        ax.errorbar(x + (i - 1) * w, m, yerr=[m - lo, hi - m], fmt="none", ecolor="#404040", elinewidth=0.6, capsize=0)
    ax.axhline(0, color=GRAY, lw=0.6); ax.set_xticks(x); ax.set_xticklabels([dt.date.fromisoformat(d).strftime("%d %b") for d in dates], fontsize=6.5)
    ax.set_ylabel("Effect of elevated CO$_2$ on $g_{sw}$ (%)"); letter(ax, 3)
    ax.legend(fontsize=5.8, frameon=False, loc="lower left", ncol=3, bbox_to_anchor=(0.0, 1.0))
    h = [Line2D([], [], marker="o", ls="", color=C_A, mfc=C_A, label="porometer, ambient CO$_2$ (mean ± sd over the four treatments)"), Line2D([], [], marker="o", ls="", color=C_E, mfc=C_E, label="porometer, elevated CO$_2$"),
         Line2D([], [], marker="o", ls="", color=C_A, mfc="white", label="LI-6800, ambient CO$_2$"), Line2D([], [], marker="o", ls="", color=C_E, mfc="white", label="LI-6800, elevated CO$_2$"),
         Line2D([], [], marker="s", ls="--", color=C_A, mfc="white", label="model, ambient CO$_2$ (10:00–14:00)"), Line2D([], [], marker="s", ls="--", color=C_E, mfc="white", label="model, elevated CO$_2$")]
    fig.legend(handles=h, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.03), fontsize=5.8, frameon=False)
    fig.tight_layout(h_pad=1.4, w_pad=1.6, rect=(0, 0.04, 1, 1))
    save(fig, out, "FigS8_porometer_WDF25")
    # numbers for the text
    txt = []
    for q in ("gs", "Tr", "dT"):
        a = df[q + "_poro"].values; b = df[q + "_mod"].values; m = np.isfinite(a) & np.isfinite(b)
        txt.append("%s: n=%d, poro mean %.3f, model %.3f, r=%.2f, RMSD=%.3f, bias=%+.3f" % (q, m.sum(), a[m].mean(), b[m].mean(), np.corrcoef(a[m], b[m])[0, 1], math.sqrt(np.mean((b[m] - a[m]) ** 2)), np.mean(b[m] - a[m])))
        if q != "dT":
            c = df[q + "_li"].values; m2 = np.isfinite(a) & np.isfinite(c)
            txt.append("   poro vs LI-6800: n=%d, r=%.2f, LI mean %.3f, ratio poro/LI %.2f" % (m2.sum(), np.corrcoef(a[m2], c[m2])[0, 1], c[m2].mean(), a[m2].mean() / c[m2].mean()))
    for col in ("gs_poro", "gs_li", "gs_mod", "Tr_poro", "Tr_li", "Tr_mod", "dT_poro", "dT_mod"):
        m, lo, hi = eff(col); txt.append("CO2 effect %s by date: %s ; pooled mean %+.0f %%" % (col, ", ".join("%+.0f" % v for v in m), np.nanmean(m)))
    sys.stderr.write("\n".join(txt) + "\n")
    return df


