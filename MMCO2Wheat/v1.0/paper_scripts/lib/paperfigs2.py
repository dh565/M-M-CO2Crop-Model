#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Part 2 figures that pool WheatDryFACE and Maricopa: Fig. 1 (energy balance), Fig. 3 (leaf gas exchange) and
Fig. 10 (canopy temperature).  Shares style, loaders and helpers with paperfigs.py."""
import csv
import datetime as dt
import math
import os
import sys
from collections import defaultdict, OrderedDict

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
import maricopa_eb as ME
from paperfigs import (C_A, C_E, C_OBS, C_SIM, SITE_C, SITE_M, LEG_N, GRAY, W1, W2, CM, DECKS, FACE, EB, P2, isnan, NAN,
                       letter, save, onetoone, place_stats, datefmt, d2d)

V3 = os.path.join(config.WDF25_V3_DIR, "S3_%s.csv")
MAR_SEASONS = [("1993/94 wet", 906, 908), ("1995/96 high N", 910, 912), ("1996/97 high N", 914, 916), ("1996/97 low N", 913, 915)]


# ------------------------------------------------------------------ WDF25 FLIR: canopy minus air per treatment and date
def flir_dT_by_treatment():
    """{date: {treatment: (dT_obs, sd, Tleaf_obs, sd)}} from the WheatDryFACE 2024/25 treatment files (DTleaf_air_flir, Tleaf_flir)."""
    out = defaultdict(dict)
    for cv in ("Ruta", "Zahir"):
        for co2 in ("aCO2", "eCO2"):
            for w in ("WW", "WL"):
                tr_ = "%s_%s_%s" % (cv, co2, w)
                rows = list(csv.reader(open(V3 % tr_, encoding="utf-8-sig"))); hdr = rows[1]
                i = hdr.index("DTleaf_air_flir"); j = hdr.index("Tleaf_flir")
                for r in rows[2:]:
                    if r[i]:
                        d = dt.datetime.strptime(r[0], "%m/%d/%Y").date().isoformat()
                        out[d][tr_] = (float(r[i]), float(r[i + 1]), float(r[j]), float(r[j + 1]))
    return out


def model_window(deck, date_iso, hours=(10, 14)):
    """(Tcan - Tair, Tcan, Tair) model means over the solar-hour window on that date, from hourly_energy_balance.csv."""
    pkg = P2.PKG
    key = ("mw", deck)
    if key not in ME._cache:
        ph = pd.read_csv(os.path.join(pkg, deck, "output", "daily_phenology.csv"))
        ph["date"] = pd.to_datetime(dict(year=ph.year.astype(int), month=ph.month.astype(int), day=ph.dom.astype(int)))
        d2 = dict(zip(ph.day.astype(int), ph.date.dt.date.astype(str)))
        m = pd.read_csv(os.path.join(pkg, deck, "output", "hourly_energy_balance.csv"), skipinitialspace=True)
        m.columns = [c.strip() for c in m.columns]
        m["date"] = m["day"].astype(int).map(d2)
        ME._cache[key] = m
    m = ME._cache[key]
    q = m[(m.date == date_iso) & (m.solar_hour >= hours[0]) & (m.solar_hour <= hours[1])]
    if not len(q): return NAN, NAN, NAN
    return float((q.Tcan_C - q.Tair_C).mean()), float(q.Tcan_C.mean()), float(q.Tair_C.mean())


N_REP = 3  # WheatDryFACE rings per treatment (Part 2, Table 2)


def ci95(vals, within_var):
    """95 % interval of a mean over treatment pairs: between-pair spread + mean within-pair replicate variance."""
    k = len(vals)
    if k == 0: return 0.0
    between = (np.var(vals, ddof=1) / k) if k > 1 else 0.0
    return 1.96 * math.sqrt(between + within_var / k)


# ------------------------------------------------------------------ Fig 10: canopy temperature
def fig_tcan(out):
    fig = plt.figure(figsize=(W2, 1.12 * W2))
    gs = fig.add_gridspec(3, 6, height_ratios=[1.55, 1, 1.05], hspace=0.55, wspace=0.9)
    trt = flir_dT_by_treatment(); dates = sorted(trt)
    # (a) WDF25 FLIR: canopy - air by date, aCO2 / eCO2, obs (filled, mean +- sd of 4 treatments) and model (open squares)
    ax = fig.add_subplot(gs[0, :4])
    wdf_eff = []; XO, YM, SITE = [], [], []
    for co2, c in (("aCO2", C_A), ("eCO2", C_E)):
        xo, yo, eo, ym = [], [], [], []
        for d in dates:
            cc = [k for k in trt[d] if co2 in k]
            vo = [trt[d][k][0] for k in cc]; vm = [model_window("WDF25_" + k, d)[0] for k in cc]
            xo.append(dt.date.fromisoformat(d)); yo.append(np.mean(vo)); eo.append(np.std(vo, ddof=1)); ym.append(np.mean(vm))
            XO += vo; YM += vm; SITE += ["WDF25"] * len(vo)
        ax.errorbar(xo, yo, yerr=eo, fmt="o", color=c, mfc=c, ms=3.6, elinewidth=0.7, capsize=0, label="thermal images, %s CO$_2$" % ("ambient" if co2 == "aCO2" else "elevated"))
        ax.plot(xo, ym, "s--", color=c, mfc="white", ms=3.4, lw=0.9, label="model, %s CO$_2$" % ("ambient" if co2 == "aCO2" else "elevated"))
    for d in dates:
        pe = [(k, k.replace("eCO2", "aCO2")) for k in trt[d] if "eCO2" in k]
        do = [trt[d][e][0] - trt[d][a][0] for e, a in pe]
        wv = np.mean([(trt[d][e][1] ** 2 + trt[d][a][1] ** 2) / N_REP for e, a in pe])
        dm = [model_window("WDF25_" + e, d)[0] - model_window("WDF25_" + a, d)[0] for e, a in pe]
        wdf_eff.append((d, np.mean(do), ci95(do, wv), np.mean(dm), 1.96 * np.std(dm, ddof=1) / math.sqrt(len(dm))))
    ax.axhline(0, color=GRAY, lw=0.6)
    ax.set_ylabel("Canopy − air temperature (K)"); ax.xaxis.set_major_locator(mdates.DayLocator(interval=7)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.set_title("WheatDryFACE 2024/25: thermal images, 10:00–14:00", pad=3, fontsize=7.5)
    y0_, y1_ = ax.get_ylim(); ax.set_ylim(y0_, y1_ + 0.45 * (y1_ - y0_)); ax.legend(fontsize=5.5, ncol=2, loc="upper center", bbox_to_anchor=(0.58, 1.0), frameon=False); letter(ax, 0)
    # (b) 1:1 of midday canopy - air, all sites (WDF25 treatments x dates; Maricopa days, Feb-Apr)
    axb = fig.add_subplot(gs[0, 4:])
    Xm, Ym = [], []
    for lab, ta, te in MAR_SEASONS:
        for tr in (ta, te):
            for d, o, _, m, *_ in ME.midday_dT(P2.PKG, tr):
                if d.month in (2, 3, 4): Xm.append(o); Ym.append(m)
    axb.plot(Xm, Ym, SITE_M["MARICOPA"], mfc="none", mec=SITE_C["MARICOPA"], mew=0.45, ms=2.4, alpha=0.5, zorder=2)
    axb.plot(XO, YM, SITE_M["WDF25"], mfc=SITE_C["WDF25"], mec=SITE_C["WDF25"], mew=0.6, ms=3.2, zorder=4)
    def _st(x, y):
        x = np.array(x, float); y = np.array(y, float); m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
        return "n = %d, r = %.2f, RMSD = %.2f K, bias %+.2f K" % (len(x), np.corrcoef(x, y)[0, 1], np.sqrt(np.mean((y - x) ** 2)), np.mean(y - x))
    sys.stderr.write("FIG10B WDF %s | MAR %s\n" % (_st(XO, YM), _st(Xm, Ym)))
    onetoone(axb, XO + Xm, YM + Ym)
    def _sts(x, y):
        x = np.array(x, float); y = np.array(y, float); m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
        return "r %.2f, RMSD %.1f K, n %d" % (np.corrcoef(x, y)[0, 1], np.sqrt(np.mean((y - x) ** 2)), len(x))
    axb.text(0.04, 0.865, "WDF: %s\nMaricopa: %s" % (_sts(XO, YM), _sts(Xm, Ym)), transform=axb.transAxes, ha="left", va="top", fontsize=4.6, zorder=7, bbox=dict(facecolor="white", edgecolor="none", pad=0.4, alpha=0.9))
    hb = [Line2D([], [], marker=SITE_M["WDF25"], ls="", mfc=SITE_C["WDF25"], mec=SITE_C["WDF25"], label="WheatDryFACE"),
          Line2D([], [], marker=SITE_M["MARICOPA"], ls="", mfc="none", mec=SITE_C["MARICOPA"], label="Maricopa")]
    axb.legend(handles=hb, loc="lower right", fontsize=5.8, handletextpad=0.3, borderpad=0.3, frameon=False, ncol=1)
    axb.set_xlabel("Observed canopy − air (K)"); axb.set_ylabel("Simulated canopy − air (K)"); axb.set_title("Midday, both sites", pad=3, fontsize=7); letter(axb, 1)
    # (c-e) Maricopa: midday canopy - air, weekly means from 1 Feb, obs vs model, aCO2/eCO2, three high-N seasons
    mar_eff = []; axs = []
    for k, (lab, ta, te) in enumerate(MAR_SEASONS[:3]):
        ax = fig.add_subplot(gs[1, 2 * k:2 * k + 2]); axs.append(ax)
        for trno, c in ((ta, C_A), (te, C_E)):
            df = pd.DataFrame(ME.midday_dT(P2.PKG, trno), columns=["date", "o", "osd", "m", "tco", "tcm", "tao", "tam"])
            df = df[df.date.dt.month.isin([2, 3, 4, 5])]
            wk = df.set_index("date").resample("7D").agg(o=("o", "mean"), m=("m", "mean"), n=("o", "size"), osd=("o", "std")).dropna()
            wk = wk[wk.n >= 3]
            ax.errorbar(wk.index, wk.o, yerr=wk.osd, fmt="o", color=c, mfc=c, ms=2.8, elinewidth=0.5, capsize=0)
            ax.plot(wk.index, wk.m, "s--", color=c, mfc="white", ms=2.6, lw=0.8)
        ax.axhline(0, color=GRAY, lw=0.6); datefmt(ax)
        ax.set_title("Maricopa %s" % lab.replace("wet", "well-watered"), pad=3, fontsize=7); letter(ax, 2 + k)
        if k == 0: ax.set_ylabel("Canopy − air temperature (K)")
    lo = min(a.get_ylim()[0] for a in axs); hi = max(a.get_ylim()[1] for a in axs)
    for a in axs: a.set_ylim(lo, hi)
    for lab, ta, te in MAR_SEASONS:
        a = {d: (o, m) for d, o, _, m, *_ in ME.midday_dT(P2.PKG, ta)}
        e = {d: (o, m) for d, o, _, m, *_ in ME.midday_dT(P2.PKG, te)}
        ds = [d for d in a if d in e and d.month in (2, 3, 4)]
        do = [e[d][0] - a[d][0] for d in ds]; dm = [e[d][1] - a[d][1] for d in ds]
        mar_eff.append((lab, np.mean(do), 1.96 * np.std(do, ddof=1) / math.sqrt(len(do)), np.mean(dm), 1.96 * np.std(dm, ddof=1) / math.sqrt(len(dm)), len(ds)))
    # (f) effect of eCO2 on midday canopy temperature: WDF25 dates + Maricopa seasons
    ax = fig.add_subplot(gs[2, :])
    short = {"1993/94 wet": "1993/94\nwell-watered", "1995/96 high N": "1995/96\nhigh N", "1996/97 high N": "1996/97\nhigh N", "1996/97 low N": "1996/97\nlow N"}
    labels = [dt.date.fromisoformat(d).strftime("%d %b") for d, *_ in wdf_eff] + [short.get(l, l) for l, *_ in mar_eff]
    O = [w[1] for w in wdf_eff] + [m[1] for m in mar_eff]; OE = [w[2] for w in wdf_eff] + [m[2] for m in mar_eff]
    M = [w[3] for w in wdf_eff] + [m[3] for m in mar_eff]; MEe = [w[4] for w in wdf_eff] + [m[4] for m in mar_eff]
    x = np.arange(len(labels)); w = 0.38
    ax.bar(x - w / 2, O, w, yerr=OE, color=C_OBS, error_kw=dict(elinewidth=0.7, capsize=0), label="observed (mean ± 95 % CI)")
    ax.bar(x + w / 2, M, w, yerr=MEe, color=C_SIM, error_kw=dict(elinewidth=0.7, capsize=0), label="model (mean ± 95 % CI over the same pairs or days)")
    ax.axhline(0, color=GRAY, lw=0.6); ax.axvline(len(wdf_eff) - 0.5, color=GRAY, lw=0.6, ls=":")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6.5)
    yl0_, yl1_ = ax.get_ylim(); ax.set_ylim(yl0_, yl1_ + 0.18 * (yl1_ - yl0_))
    ax.text((len(wdf_eff) - 1) / 2.0, 0.97, "WheatDryFACE 2024/25 (thermal images, by date)", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.5)
    ax.text(len(wdf_eff) + (len(mar_eff) - 1) / 2.0, 0.97, "Maricopa (infrared thermometers, Feb–Apr)", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.5)
    ax.set_ylabel("Elevated − ambient CO$_2$ (K)"); ax.legend(fontsize=6.5, loc="upper right", bbox_to_anchor=(1.0, 0.9))
    ax.set_title("Effect of elevated CO$_2$ on midday canopy temperature", pad=3, fontsize=7.5)
    letter(ax, 5)
    save(fig, out, "FigN10_canopy_temperature_all_sites")
    return wdf_eff, mar_eff


# ------------------------------------------------------------------ Fig 1: energy balance, WDF + Maricopa
def fig_eb(out):
    keys = [("Rn_day", "Net radiation, day (W m$^{-2}$)", True), ("Rn_night", "Net radiation, night (W m$^{-2}$)", True),
            ("LE_day", "Latent heat flux, day (W m$^{-2}$)", False), ("H_day", "Sensible heat flux, day (W m$^{-2}$)", False),
            ("G_24h", "Soil heat flux, 24-h (W m$^{-2}$)", True), ("LWout_day", "Outgoing longwave, day (W m$^{-2}$)", True),
            ("albedo_day", "Albedo, day (–)", True), ("Tsurf_rad", "Surface temperature, day (°C)", True)]
    fig, axes = plt.subplots(2, 4, figsize=(W2, 0.55 * W2))
    marv = {"Rn_day": ("rn", "Rn_total_W_m2", False), "Rn_night": ("rn", "Rn_total_W_m2", True), "LE_day": ("le", "LE", False),
            "H_day": ("h", "H_total_W_m2", False), "G_24h": ("g", "G_W_m2", None)}
    for k, (ax, (q, lab, wdf)) in enumerate(zip(axes.flat, keys)):
        X, Y = [], []
        if q in marv:
            vo, vm, night = marv[q]; xs, ys = [], []
            for lab_, ta, te in MAR_SEASONS:
                for tr in (ta, te):
                    if night is None:
                        dd = ME.daily_means(P2.PKG, tr, vo, vm, hours=(1, 24))
                    else:
                        dd = ME.daily_means(P2.PKG, tr, vo, vm, night=night)
                    xs += [v[1] for v in dd]; ys += [v[3] for v in dd]
            if xs:
                ax.plot(xs, ys, SITE_M["MARICOPA"], mfc="none", mec=SITE_C["MARICOPA"], mew=0.45, ms=2.4, alpha=0.55, zorder=2); X += xs; Y += ys
        if wdf:
            xs, ys = [], []
            for d in [x for x in FACE if DECKS[x]["site"] in ("WDF25", "WDF23")]:
                pp = EB.pairs(d, q); xs += [p[1] for p in pp]; ys += [p[3] for p in pp]
            if xs:
                ax.plot(xs, ys, SITE_M["WDF25"], mfc="none", mec=SITE_C["WDF25"], mew=0.5, ms=2.2, alpha=0.85, zorder=3); X += xs; Y += ys
        if not wdf:
            ax.text(0.97, 0.04, "Maricopa only", transform=ax.transAxes, ha="right", va="bottom", fontsize=6, color=GRAY)
        onetoone(ax, X, Y); place_stats(ax, X, Y)
        ax.set_title(lab, pad=3, fontsize=7); ax.set_xlabel("Observed"); ax.set_ylabel("Simulated"); letter(ax, k)
    h = [Line2D([], [], marker=SITE_M[s], ls="", mfc="none", mec=SITE_C[s], label=LEG_N[s]) for s in ("WDF25", "MARICOPA")]
    fig.legend(handles=h, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.03), fontsize=5.8, frameon=False)
    fig.tight_layout(h_pad=1.4, w_pad=1.2, rect=(0, 0.03, 1, 1))
    save(fig, out, "FigN01_energy_balance_all_sites")


# ------------------------------------------------------------------ Fig 3: gas exchange refigure (campaign dates, both sides with sd)
def fig_gx(out):
    keys = [("An", "$A_n$ (µmol m$^{-2}$ s$^{-1}$)"), ("gs", "$g_{sw}$ (mol m$^{-2}$ s$^{-1}$)"), ("Tr", "$T_r$ (mmol m$^{-2}$ s$^{-1}$)")]
    seasons = [("WDF23", "WDF 2022/23 (low N)"), ("WDF25", "WDF 2024/25 (high N)")]
    fig, axes = plt.subplots(5, 3, figsize=(W2, 1.50 * W2))
    k = 0
    # rows 1-2: per campaign date, treatment means: obs (filled, +-sd over treatments) and model (open square, window mean)
    for i, (site, ttl) in enumerate(seasons):
        for j, (q, lab) in enumerate(keys):
            ax = axes[i, j]
            for co2, c in (("aCO2", C_A), ("eCO2", C_E)):
                O = defaultdict(list); M = defaultdict(list)
                for d in [x for x in FACE if DECKS[x]["site"] == site and DECKS[x]["co2"] == co2]:
                    for p in EB.pairs(d, q):
                        O[p[0]].append(p[1]); M[p[0]].append(p[3])
                ds = sorted(O); xo = [d2d(x) for x in ds]
                if not ds: continue
                ax.errorbar(xo, [np.mean(O[x]) for x in ds], yerr=[np.std(O[x], ddof=1) if len(O[x]) > 1 else 0 for x in ds],
                            fmt="o", color=c, mfc=c, ms=3.4, elinewidth=0.6, capsize=0, zorder=4)
                ax.errorbar(xo, [np.mean(M[x]) for x in ds], yerr=[np.std(M[x], ddof=1) if len(M[x]) > 1 else 0 for x in ds],
                            fmt="s--", color=c, mfc="white", ms=3.2, lw=0.8, elinewidth=0.6, capsize=0, zorder=3)
            datefmt(ax); ax.set_ylim(bottom=0)
            ax.set_title("%s: %s" % (["$A_n$", "$g_{sw}$", "$T_r$"][j], ttl), pad=3, fontsize=7); ax.set_ylabel(lab); letter(ax, k); k += 1
    for j in range(3):
        top = max(axes[0, j].get_ylim()[1], axes[1, j].get_ylim()[1]); axes[0, j].set_ylim(0, top); axes[1, j].set_ylim(0, top)
    # rows 3-4: eCO2 / aCO2 effect per campaign date, one season per row
    for i, (site, ttl) in enumerate(seasons):
        for j, (q, lab) in enumerate(keys):
            ax = axes[2 + i, j]
            RO, RM, SEo, SEm, xs, LO, HI = [], [], [], [], [], [], []
            A = defaultdict(dict)
            for d in [x for x in FACE if DECKS[x]["site"] == site]:
                for p in EB.pairs(d, q): A[p[0]][d] = (p[1], p[3], p[2])
            for x in sorted(A):
                pe = [(e, e.replace("eCO2", "aCO2")) for e in A[x] if "eCO2" in e and e.replace("eCO2", "aCO2") in A[x]]
                if not pe: continue
                ro = [A[x][e][0] / A[x][a][0] for e, a in pe if A[x][a][0] > 0]; rm = [A[x][e][1] / A[x][a][1] for e, a in pe if A[x][a][1] > 0]
                if not ro or not rm: continue
                wv = np.mean([(A[x][e][0] / A[x][a][0]) ** 2 * ((A[x][e][2] / A[x][e][0]) ** 2 + (A[x][a][2] / A[x][a][0]) ** 2) / N_REP
                              for e, a in pe if A[x][a][0] > 0 and A[x][e][0] > 0 and not (isnan(A[x][e][2]) or isnan(A[x][a][2]))] or [0.0])
                xs.append(d2d(x)); RO.append((np.mean(ro) - 1) * 100); SEo.append(ci95(ro, wv) * 100); RM.append((np.mean(rm) - 1) * 100)
                SEm.append((np.std(rm, ddof=1) if len(rm) > 1 else 0.0) * 100)
                LO.append((min(rm) - 1) * 100); HI.append((max(rm) - 1) * 100)
            # model spread drawn as in Figs 7-9: thick translucent red bar over the min-max range of the eCO2/aCO2 ratios
            # across treatment pairs, plus a thin red bar of +- sd across pairs
            for x0, lo, hi in zip(xs, LO, HI):
                ax.plot([x0, x0], [lo, hi], color=C_SIM, lw=2.6, alpha=0.45, solid_capstyle="butt", zorder=2)
            ax.errorbar(xs, RO, yerr=SEo, fmt="o", color=C_OBS, ms=3.4, elinewidth=0.6, capsize=0, zorder=4)
            ax.errorbar(xs, RM, yerr=SEm, fmt="s--", color=C_SIM, mfc="white", ms=3.2, lw=0.8, elinewidth=0.8, capsize=0, zorder=3)
            ax.axhline(0, color=GRAY, lw=0.6); datefmt(ax)
            ax.set_ylabel("Effect of elevated CO$_2$ on %s (%%)" % ["$A_n$", "$g_{sw}$", "$T_r$"][j]); ax.set_title("CO$_2$ effect, %s" % ttl.split(" (")[0], pad=3, fontsize=7); letter(ax, k); k += 1
    for j in range(3):
        lo = min(axes[2, j].get_ylim()[0], axes[3, j].get_ylim()[0]); hi = max(axes[2, j].get_ylim()[1], axes[3, j].get_ylim()[1])
        axes[2, j].set_ylim(lo, hi); axes[3, j].set_ylim(lo, hi)
    # row 5: Maricopa leaf gas exchange at natural light (Wall et al., 2000, 2006), 10:00-14:00 MST
    MAR = [("1993/94 dry", 905, 907), ("1993/94 wet", 906, 908), ("1995/96 high N", 910, 912), ("1996/97 high N", 914, 916), ("1996/97 low N", 913, 915)]
    for j, (q, lab) in enumerate(keys[:2]):
        # (m, n) pooled comparison for both sites that measured leaf gas exchange: WheatDryFACE campaign
        # means (both seasons, purple) and Maricopa dates (gold); site colours as in Figs 1, 2 and 6
        ax = axes[4, j]; Xw, Yw, Xm, Ym = [], [], [], []
        for _, ta, te in MAR:
            for tr in (ta, te):
                pp = EB.pairs("MARICOPA_%d" % tr, q)
                ax.errorbar([p[1] for p in pp], [p[3] for p in pp], xerr=[0 if isnan(p[2]) else p[2] for p in pp], fmt=SITE_M["MARICOPA"], color=SITE_C["MARICOPA"], mfc="none", ms=2.6, mew=0.6, elinewidth=0.4, capsize=0, alpha=0.75, zorder=2)
                Xm += [p[1] for p in pp]; Ym += [p[3] for p in pp]
        for dk in [x for x in FACE if DECKS[x]["site"] in ("WDF25", "WDF23")]:
            pp = EB.pairs(dk, q)
            ax.errorbar([p[1] for p in pp], [p[3] for p in pp], xerr=[0 if isnan(p[2]) else p[2] for p in pp], fmt=SITE_M["WDF25"], color=SITE_C["WDF25"], mfc="none", ms=2.6, mew=0.6, elinewidth=0.4, capsize=0, alpha=0.85, zorder=3)
            Xw += [p[1] for p in pp]; Yw += [p[3] for p in pp]
        X = Xw + Xm; Y = Yw + Ym
        onetoone(ax, X, Y)
        def _st3(x, y):
            x = np.array(x, float); y = np.array(y, float); m_ = np.isfinite(x) & np.isfinite(y); x, y = x[m_], y[m_]
            return "n = %d, r = %.2f, RMSD = %.2f" % (len(x), np.corrcoef(x, y)[0, 1], np.sqrt(np.mean((y - x) ** 2)))
        place_stats(ax, X, Y, txt="WDF: %s\nMaricopa: %s" % (_st3(Xw, Yw), _st3(Xm, Ym)), fs=4.4)
        sys.stderr.write("FIG3%s %s WDF %s | Maricopa %s\n" % ("mn"[j], q, _st3(Xw, Yw), _st3(Xm, Ym)))
        ax.set_xlabel("Observed " + lab.split(" (")[0]); ax.set_ylabel("Simulated " + lab.split(" (")[0])
        ax.set_title("%s: both sites, 10:00–14:00" % ["$A_n$", "$g_{sw}$"][j], pad=3, fontsize=7); letter(ax, k); k += 1
    ax = axes[4, 2]; xl = []; xi = 0
    for q, nm in (("An", "$A_n$"), ("gs", "$g_{sw}$")):
        for lab_, ta, te in MAR:
            pa = {p[0]: p for p in EB.pairs("MARICOPA_%d" % ta, q)}; pe = {p[0]: p for p in EB.pairs("MARICOPA_%d" % te, q)}
            cm = sorted(set(pa) & set(pe))
            if not cm: continue
            oa = np.mean([pa[x][1] for x in cm]); oe = np.mean([pe[x][1] for x in cm]); ma = np.mean([pa[x][3] for x in cm]); me = np.mean([pe[x][3] for x in cm])
            rr = [pe[x][1] / pa[x][1] for x in cm if pa[x][1] > 0]
            ci = 1.96 * np.std(rr, ddof=1) / math.sqrt(len(rr)) * 100 if len(rr) > 1 else 0
            ax.bar(xi - 0.2, (oe / oa - 1) * 100, 0.38, color=C_OBS, yerr=ci, error_kw=dict(lw=0.6, capsize=0))
            ax.bar(xi + 0.2, (me / ma - 1) * 100, 0.38, color=C_SIM)
            xl.append({"1993/94 dry": "1993/94 WL", "1993/94 wet": "1993/94 WW", "1995/96 high N": "1995/96 high N", "1996/97 high N": "1996/97 high N", "1996/97 low N": "1996/97 low N"}[lab_]); xi += 1
        if q == "An": ax.axvline(xi - 0.5, color=GRAY, lw=0.5, ls=":")
    ax.text(2, 0.97, "$A_n$", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.5)
    ax.text(7, 0.03, "$g_{sw}$", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.5)
    ax.axhline(0, color=GRAY, lw=0.6); ax.set_xticks(range(len(xl))); ax.set_xticklabels(xl, fontsize=5.4, rotation=60, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Effect of elevated CO$_2$ (%)"); ax.set_title("CO$_2$ effect, Maricopa (mean over dates)", pad=3, fontsize=7); letter(ax, k); k += 1
    # legend in three columns, one per panel block (a–f: seasonal course; g–l and o: CO2 effect; m, n: both sites pooled)
    h = [Line2D([], [], marker="o", ls="", color=C_A, label="observed, ambient CO$_2$ (± sd over treatments)"), Line2D([], [], marker="o", ls="", color=C_E, label="observed, elevated CO$_2$"),
         Line2D([], [], marker="s", ls="--", color=C_A, mfc="white", label="model, ambient CO$_2$ (10:00–14:00, ± sd)"), Line2D([], [], marker="s", ls="--", color=C_E, mfc="white", label="model, elevated CO$_2$"),
         Line2D([], [], marker="o", ls="", color=C_OBS, label="observed effect ± 95 % CI"), Line2D([], [], marker="s", ls="--", color=C_SIM, mfc="white", label="model effect (thin: ± sd over pairs; thick: range)"),
         Line2D([], [], marker=SITE_M["WDF25"], ls="", mfc="none", mec=SITE_C["WDF25"], label="WheatDryFACE, campaign means ± sd"), Line2D([], [], marker=SITE_M["MARICOPA"], ls="", mfc="none", mec=SITE_C["MARICOPA"], label="Maricopa, dates ± sd")]
    for hh, x0, ttl in ((h[:4], 0.17, "(a–f)"), (h[4:6], 0.52, "(g–l, o)"), (h[6:], 0.84, "(m, n)")):
        lg = fig.legend(handles=hh, loc="upper center", ncol=1, bbox_to_anchor=(x0, 0.044), fontsize=5.8, frameon=False, handletextpad=0.4, title=ttl, title_fontsize=6)
        lg._legend_box.align = "left"
    fig.tight_layout(h_pad=1.0, w_pad=1.3, rect=(0, 0.045, 1, 1))
    save(fig, out, "Fig03_flag_leaf_gas_exchange_all_sites")


