#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Supplement S2, Fig. S7: Maricopa soil temperature profiles (one instrumented plot per treatment, 1-40 cm, hourly;
Kimball et al. 1999/2017 energy-balance workbooks) against the model's daily soil temperature interpolated to the
sensor depths.  Panels (j, k) of Fig. 2 of Part 2 use obs_daily() and model_daily() of this module."""
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paperfigs as PF
import maricopa_archive as MA
from paperfigs import C_A, C_E, C_OBS, C_SIM, GRAY, W2, P2, isnan, NAN, letter, save, onetoone, place_stats, datefmt, stats_txt

SPEC = {"eb94": ("1993/94", 8, 10, {"CW": (31, 906), "FW": (37, 908)}),
        "eb96": ("1995/96", 8, 10, {"BH": (47, 910), "FH": (53, 912), "BL": (65, None), "FL": (71, None)}),
        "eb97": ("1996/97", 8, 10, {"BH": (41, 914), "FH": (47, 916), "BL": (59, 913), "FL": (65, 915)})}
DEPTHS = [0.01, 0.02, 0.05, 0.10, 0.20, 0.40]

def num(s):
    v = pd.to_numeric(s, errors="coerce"); return v.where(~v.isin([-99, -999]))

def obs_daily():
    """{trno: DataFrame(date, z, Tobs)} daily means of the hourly profile."""
    out = {}
    for tag, (season, hdr, d0, grp) in SPEC.items():
        df = MA.eb_workbook(tag).iloc[d0:].reset_index(drop=True)
        year = num(df[0]); doy = num(df[1])
        date = pd.to_datetime(year.astype("Int64").astype(str) + "-" + doy.astype("Int64").astype(str), format="%Y-%j", errors="coerce")
        for g, (c0, trno) in grp.items():
            if trno is None: continue
            rows = []
            for k, z in enumerate(DEPTHS):
                t = num(df[c0 + k])
                d = pd.DataFrame({"date": date, "T": t}).dropna()
                dd = d.groupby("date").agg(T=("T", "mean"), n=("T", "size")).reset_index()
                dd = dd[dd.n >= 20]; dd["z"] = z; rows.append(dd[["date", "z", "T"]])
            out[trno] = pd.concat(rows)
    return out

def model_daily(trno):
    dk = "MARICOPA_%d" % trno; m = P2.MODELS[dk]
    rr = P2.rows_of(os.path.join(m.dir, "daily_tsoil_profile.csv"))
    bnds = P2.layer_bounds(P2.PKG, dk) or P2.BNDS
    mid = [(bnds[i] + bnds[i + 1]) / 2.0 for i in range(10)]
    mm = PF.EB.model(dk)
    out = []
    for r in rr:
        v = [P2.fl(r["Tsoil_L%d" % (i + 1)]) for i in range(10)]; day = int(P2.fl(r["day"])); ds = mm.day2date.get(day)
        if ds is None: continue
        for z in DEPTHS:
            if z <= mid[0]: t = v[0]
            else:
                t = v[-1]
                for i in range(9):
                    if mid[i] <= z <= mid[i + 1]: t = v[i] + (z - mid[i]) / (mid[i + 1] - mid[i]) * (v[i + 1] - v[i]); break
            out.append((pd.Timestamp(ds), z, t))
    return pd.DataFrame(out, columns=["date", "z", "T"])

def figS7(out):
    O = obs_daily()
    panels = [(906, "1993/94 well-watered, ambient CO$_2$"), (908, "1993/94 well-watered, elevated CO$_2$"),
              (910, "1995/96 high N, ambient CO$_2$"), (912, "1995/96 high N, elevated CO$_2$"),
              (914, "1996/97 high N, ambient CO$_2$"), (916, "1996/97 high N, elevated CO$_2$"),
              (913, "1996/97 low N, ambient CO$_2$"), (915, "1996/97 low N, elevated CO$_2$")]
    fig, axes = plt.subplots(3, 4, figsize=(W2, 0.72 * W2)); axes = list(axes.flat)
    cols = {0.10: C_OBS, 0.40: C_OBS}; ls = {0.10: "-", 0.40: "--"}
    XA = {0.10: [], 0.40: []}; YA = {0.10: [], 0.40: []}; stats = {}
    for k, (trno, ttl) in enumerate(panels):
        ax = axes[k]; o = O[trno]; m = model_daily(trno)
        for z in (0.10, 0.40):
            oo = o[o.z == z].sort_values("date"); mmz = m[m.z == z].sort_values("date")
            ax.plot(mmz.date, mmz["T"], ls[z], color=C_SIM, lw=0.9, zorder=2)
            ax.plot(oo.date, oo["T"], ls[z], color=C_OBS, lw=0.9, zorder=3)
            j = oo.merge(mmz, on="date", suffixes=("_o", "_m"))
            XA[z] += list(j.T_o); YA[z] += list(j.T_m)
        datefmt(ax); ax.set_title(ttl, pad=3, fontsize=6.6); ax.set_ylabel("$T_{soil}$ (°C)"); letter(ax, k)
    for z, ax in ((0.10, axes[8]), (0.40, axes[9])):
        X = np.array(XA[z]); Y = np.array(YA[z])
        ax.plot(X, Y, "D", mfc="none", mec=C_SIM, mew=0.45, ms=2.4, alpha=0.6)
        onetoone(ax, list(X), list(Y))
        lo, hi = np.nanmin(np.r_[X, Y]) - 1, np.nanmax(np.r_[X, Y]) + 1; ax.set_xlim(lo, hi); ax.set_ylim(lo, hi + 0.22 * (hi - lo))
        ax.text(0.04, 0.86, stats_txt(list(X), list(Y)), transform=ax.transAxes, ha="left", va="top", fontsize=5.6, zorder=7, bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.9))
        ax.set_xlabel("Observed $T_{soil}$ (°C)"); ax.set_ylabel("Simulated $T_{soil}$ (°C)"); ax.set_title("All treatments, %d cm" % int(z * 100), pad=3, fontsize=6.6); letter(ax, 8 + (z > 0.1))
        stats[z] = (len(X), np.corrcoef(X, Y)[0, 1], math.sqrt(np.mean((Y - X) ** 2)), np.mean(Y - X))
    # (k) depth profile of the mean bias per treatment
    ax = axes[10]; BTOP = []
    for trno, ttl in panels:
        o = O[trno]; m = model_daily(trno); b = []
        for z in DEPTHS:
            j = o[o.z == z].merge(m[m.z == z], on="date", suffixes=("_o", "_m")); b.append(np.mean(j.T_m - j.T_o) if len(j) else NAN)
        BTOP.append([round(v, 2) for v in b[:2]])
        ax.plot(b, [z * 100 for z in DEPTHS], "o-", ms=2.5, lw=0.8, color=C_E if trno in (908, 912, 916, 915) else C_A, alpha=0.8)
    ax.set_ylim(-2, 43); ax.invert_yaxis(); ax.axvline(0, color=GRAY, lw=0.6); ax.set_xlabel("Mean bias, model − observed (K)"); ax.set_ylabel("Depth (cm)"); ax.set_title("Bias by depth, every treatment", pad=3, fontsize=6.6); letter(ax, 10)
    axes[11].axis("off")
    h = [Line2D([], [], color=C_OBS, ls="-", label="observed, 10 cm"), Line2D([], [], color=C_SIM, ls="-", label="model, 10 cm"),
         Line2D([], [], color=C_OBS, ls="--", label="observed, 40 cm"), Line2D([], [], color=C_SIM, ls="--", label="model, 40 cm"),
         Line2D([], [], color=C_A, marker="o", ls="-", label="(k) ambient CO$_2$"), Line2D([], [], color=C_E, marker="o", ls="-", label="(k) elevated CO$_2$")]
    fig.legend(handles=h, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02), fontsize=6.3)
    fig.tight_layout(h_pad=1.2, w_pad=1.1, rect=(0, 0.04, 1, 1))
    save(fig, out, "FigS7_soil_temperature_MARICOPA")
    sys.stderr.write("S7 stats: " + "; ".join("%d cm: n=%d r=%.2f RMSD=%.2f bias=%+.2f" % (int(z * 100), *stats[z]) for z in stats) + "\n")
    return stats
