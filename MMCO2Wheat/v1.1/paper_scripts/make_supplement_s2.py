#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Supplement S2 of Part 2: Figs. S1-S8 and Tables S1-S5, from the reference-outputs archive.

Writes to OUT_DIR/supplement_s2/:
  FigS1_soil_water_WDF23 .. FigS4_soil_water_MARICOPA   soil water course of every treatment, by site
  FigS5_harvest_by_treatment                             grain yield, above-ground DM and kernel number, all 38 treatments
  FigS6_genotype_contrast_WDF                            Zahir / Ruta at WheatDryFACE (lib/paperfigs.figN_genotype)
  FigS7_soil_temperature_MARICOPA                        Maricopa soil temperature profiles (lib/maricopa_tsoil.figS7)
  FigS8_porometer_WDF25                                  WheatDryFACE 2024/25 porometer record (lib/wdf_poro.figS8)
  TableS1_WDF23.csv .. TableS4_MARICOPA.csv              simulated / observed harvest and canopy variables by treatment
  TableS5_closures_by_simulation.csv                     numerical behaviour and closures per simulation (from run.log)
  s2_tables.json                                         Tables S1-S5 as one JSON (T1 by site, T2)

    python3 make_supplement_s2.py [--only S1,S5,T,...]   (T = the tables)
"""
import argparse
import csv
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
import config
import paperfigs as PF
from paperfigs import DECKS, FACE, SITE_N, SITE_C, W2, GRAY, C_OBS, C_SIM, letter, save, isnan, NAN

SITES = ["WDF23", "WDF25", "BRAUNSCHWEIG", "MARICOPA"]
TABLE_OF_SITE = {"WDF23": "TableS1_WDF23", "WDF25": "TableS2_WDF25", "BRAUNSCHWEIG": "TableS3_BRAUNSCHWEIG", "MARICOPA": "TableS4_MARICOPA"}


def decks_of(site):
    return [d for d in FACE if DECKS[d]["site"] == site]


def mt(x):
    return x.replace(" kg N ha−1", " kg N ha$^{-1}$")


def tlabel(d):
    m = DECKS[d]
    if m["site"].startswith("WDF"):
        return "%s, %s%s" % (m["cultivar"], "ambient" if m["co2"] == "aCO2" else "elevated", ", water-limited" if m.get("water") == "WL" else (", well-watered" if m["site"] == "WDF25" else ""))
    if m["site"] == "BRAUNSCHWEIG":
        nl = m["nlevel"].split("(")[1].split(")")[0].replace(" kgN/ha", "")
        nl = nl.replace(" split", " kg N ha−1 (warming control)") if "split" in nl else nl + " kg N ha−1"
        return "%s, %s, %s" % (m["season"], nl, "ambient" if m["co2"] == "aCO2" else "elevated")
    t = m["treatment"]
    w = "water-limited" if t.startswith("Dry") else ("low N" if "Low_N" in t else ("well-watered" if m["season"] == "1993/94" else "high N"))
    return "%s, %s, %s" % (m["season"], w, "ambient" if "aCO2" in t else "elevated")


# ------------------------------------------------------------------ Tables S1-S4: harvest metrics by treatment
QS = [("grain_yield", "Grain yield (g m−2)", "%.0f"), ("AGB_harvest", "Above-ground DM (g m−2)", "%.0f"), ("kernel_number", "Kernels (m−2)", "%.0f"),
      ("TKW", "Kernel mass (g)", "%.1f"), ("LAI_peak", "Peak LAI", "%.2f"), ("grain_N_pct", "Grain N (%)", "%.2f"), ("crop_N_harvest", "Above-ground N (g N m−2)", "%.1f")]


def fmt(v, f):
    return "–" if v is None or isnan(v) else f % v


def tables_s1_s4(out):
    P2 = PF.P2
    T1 = {}
    for site in SITES:
        rows = []
        for d in decks_of(site):
            row = [tlabel(d)]
            for q, lab, f in QS:
                if q == "crop_N_harvest":
                    o, s = PF._cropN_obs(d); mv = P2.CROPN[d][0] if isinstance(P2.CROPN.get(d), (tuple, list)) else P2.CROPN.get(d, NAN)
                else:
                    o, s = PF.obs(d, q); mv = PF.mod(d, q)
                row.append("%s / %s" % (fmt(mv, f), fmt(o, f)))
            rows.append(row)
        T1[site] = rows
        with open(os.path.join(out, TABLE_OF_SITE[site] + ".csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["treatment"] + [l + " model / observed" for _, l, _ in QS])
            for r in rows:
                w.writerow(r)
    return T1


# ------------------------------------------------------------------ Table S5: closures per simulation from run.log
def grab(txt, pat, cast=float):
    m = re.search(pat, txt)
    return cast(m.group(1)) if m else None


def table_s5(out):
    T2 = []
    for d in FACE:
        txt = open(os.path.join(config.REF_OUTPUTS, d, "run.log"), errors="ignore").read()
        tot = grab(txt, r"unconverged days = \d+ of (\d+)", int)
        nun = grab(txt, r"unconverged days = (\d+) of", int)
        T2.append([tlabel(d), SITE_N[DECKS[d]["site"]],
                   "%d of %d" % (tot - nun, tot) if tot is not None else "–",
                   fmt(grab(txt, r"max resid = +([0-9.E+-]+) K"), "%.2f"),
                   grab(txt, r"Richards non-convergence events over ALL outer passes: (\d+) event", str) or "0",
                   fmt(abs(grab(txt, r"Total Hydro MB: +([0-9.E+-]+)") or 0), "%.3f"),
                   "%.1e" % grab(txt, r"heat-storage closure residual\| = +([0-9.E+-]+)"),
                   "%.1e" % grab(txt, r"CARBON CLOSURE max\|organ-state - dm_gain\| = +([0-9.E+-]+)"),
                   "%.1e" % grab(txt, r"gate N identity max\|resid\| = +([0-9.E+-]+)")])
    with open(os.path.join(out, "TableS5_closures_by_simulation.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["simulation", "site", "days converged to 0.05 K", "max residual (K)", "Richards events", "seasonal water residual (mm)", "heat residual (J m-2)", "carbon closure (g m-2)", "N payment identity (g N m-2)"])
        for r in T2:
            w.writerow(r)
    return T2


# ------------------------------------------------------------------ Figures S1-S4: soil water course of every treatment
SWQ = {"WDF23": "SWC20", "WDF25": "SWC20", "BRAUNSCHWEIG": "d0_20", "MARICOPA": "MAR_SWC_0_30"}


def fig_soil_water(out, k):
    site = SITES[k]
    ds = decks_of(site); n = len(ds); nc = 4; nr = int(np.ceil(n / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(W2, 0.22 * W2 * nr), sharey=True)
    axes = np.atleast_2d(axes)
    for j, d in enumerate(ds):
        ax = axes[j // nc, j % nc]
        try:
            PF._swc_course(ax, d, SWQ[site], legend=False)
        except Exception:
            ax.text(0.5, 0.5, "no record", transform=ax.transAxes, ha="center", color=GRAY, fontsize=6)
        ttl = mt(tlabel(d))
        if site == "BRAUNSCHWEIG": ttl = ttl.replace(" (warming control)", "\n(warming control)")
        ax.set_title(ttl, fontsize=6.0, pad=2)
        if site == "MARICOPA": ax.set_ylim(0, 0.5)
        letter(ax, j, x=-0.22)
    if site == "BRAUNSCHWEIG":   # common x-range per season
        for row0, row1 in ((0, 1), (2, 3)):
            lo = min(axes[r, c].get_xlim()[0] for r in (row0, row1) for c in range(nc)); hi = max(axes[r, c].get_xlim()[1] for r in (row0, row1) for c in range(nc))
            for r in (row0, row1):
                for c in range(nc): axes[r, c].set_xlim(lo, hi)
    for j in range(n, nr * nc): axes[j // nc, j % nc].axis("off")
    fig.supylabel("Soil water content (m$^3$ m$^{-3}$)", fontsize=7, x=0.005)
    hs = [Line2D([], [], color=C_OBS, lw=0.9, label="observed")] + ([] if site == "MARICOPA" else [Patch(color=C_OBS, alpha=0.15, label="observed ± sd (replicate plots)")]) + [Line2D([], [], color=C_SIM, lw=0.9, label="model")]
    fig.legend(handles=hs, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.01 if nr > 1 else -0.05), fontsize=6.5)
    fig.tight_layout(h_pad=1.2, w_pad=1.0, rect=(0.02, 0.03 if nr > 2 else 0.06, 1, 1))
    save(fig, out, "FigS%d_soil_water_%s" % (k + 1, site))


# ------------------------------------------------------------------ Figure S5: harvest variables by treatment, all 38
def fig_s5(out):
    ORD = [d for site in SITES for d in decks_of(site)]      # 2022/23, 2024/25, Braunschweig, Maricopa (main-paper order)
    fig, axes = plt.subplots(3, 1, figsize=(W2, 0.9 * W2), sharex=True)
    x = np.arange(len(ORD)); cols = [SITE_C[DECKS[d]["site"]] for d in ORD]
    for ax, (q, lab, f) in zip(axes, QS[:3]):
        O = [PF.obs(d, q)[0] for d in ORD]; S = [PF.obs(d, q)[1] for d in ORD]; Mv = [PF.mod(d, q) for d in ORD]
        ax.bar(x - 0.2, [0 if isnan(v) else v for v in O], 0.4, color=cols, alpha=0.45, label="observed")
        ax.errorbar(x - 0.2, [np.nan if isnan(v) else v for v in O], yerr=[0 if isnan(s) else s for s in S], fmt="none", ecolor=GRAY, elinewidth=0.5, capsize=0)
        ax.bar(x + 0.2, [0 if isnan(v) else v for v in Mv], 0.4, color=cols, label="model")
        ax.set_ylabel(lab.replace("DM", "dry matter").replace(" (", "\n(").replace("m−2", "m$^{-2}$")); ax.set_xlim(-0.6, len(ORD) - 0.4)
    for k_, ax in enumerate(axes): letter(ax, k_, x=-0.06)
    axes[-1].set_xticks(x); axes[-1].set_xticklabels([mt(tlabel(d)) for d in ORD], rotation=90, fontsize=5)
    i0 = 0
    for site in SITES:
        n_ = len(decks_of(site)); xm = i0 + (n_ - 1) / 2.0
        axes[0].text(xm, 1.02, SITE_N[site].replace("WheatDryFACE ", "WheatDryFACE\n"), transform=axes[0].get_xaxis_transform(), ha="center", va="bottom", fontsize=6.3, color=SITE_C[site])
        if i0:
            for ax in axes: ax.axvline(i0 - 0.5, color=GRAY, lw=0.5, ls=":")
        i0 += n_
    fig.text(0.5, -0.005, "light bars: observed (± sd over replicate plots); solid bars: model", fontsize=6.5, va="top", ha="center", color="#404040")
    fig.tight_layout(h_pad=0.8, rect=(0, 0, 1, 0.97))
    save(fig, out, "FigS5_harvest_by_treatment")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="comma list among S1,S2,S3,S4,S5,S6,S7,S8,T (all tables) or T5 (Table S5 alone, archive only)")
    a = ap.parse_args()
    config.require(config.REF_OUTPUTS, "reference-outputs archive")
    out = config.out_subdir("supplement_s2")
    want = set(a.only.split(",")) if a.only else {"S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "T"}
    if want == {"T5"}:
        T2 = table_s5(out)
        print("rows T2", len(T2)); return
    PF.init(config.REF_OUTPUTS)
    if "T" in want:
        T1 = tables_s1_s4(out)
        T2 = table_s5(out)
        json.dump({"T1": T1, "T2": T2, "QS": [(q, l) for q, l, _ in QS]}, open(os.path.join(out, "s2_tables.json"), "w"))
        print("rows T1", sum(len(v) for v in T1.values()), "T2", len(T2))
    for k in range(4):
        if "S%d" % (k + 1) in want:
            fig_soil_water(out, k)
    if "S5" in want:
        fig_s5(out)
    if "S6" in want:
        PF.figN_genotype(out)
    if "S7" in want:
        import maricopa_tsoil
        maricopa_tsoil.figS7(out)
    if "S8" in want:
        import wdf_poro
        wdf_poro.figS8(out)
    print("written to", out)


if __name__ == "__main__":
    main()
