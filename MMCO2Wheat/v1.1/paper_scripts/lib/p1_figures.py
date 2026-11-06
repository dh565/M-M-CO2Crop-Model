#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Part 1 figures: Fig. 1 (model structure, no data) and Figs. 2-5 (Sect. 14, from the one-at-a-time sensitivity runs
summarised in all_metrics.json by run_sensitivity.py).

  fig1(out)        P1_Fig01_model_structure
  fig2(out, M)     P1_Fig02_driver_responses          (drivers one at a time)
  fig3(out, M)     P1_Fig03_alt_cascade_elasticity    (CO2 response curve and elasticities)
  fig4(out, M)     P1_Fig04_CO2_interactions          (550 vs 420 across water, N, temperature)
  fig5(out, M)     P1_Fig05_alt_matrix_tornado        (parameter sensitivity indices and yield tornado)
M = the dict of all_metrics.json: run name -> metrics (see run_sensitivity.metrics).
"""
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paperfigs import letter, save, W2, W1, CM, GRAY, C_OBS, C_SIM, C_A, C_E

INK = "#333333"; BAR = "#4d4d4d"
TRUNC = {"T-2", "co2_550_T-2", "P_thermal_time_thresholds_x1.2", "P550_thermal_time_thresholds_x1.2"}  # runs that do not reach maturity before the end of the weather record
M = {}
ctrl = {}


def set_metrics(metrics):
    """Bind the all_metrics dict used by figs. 2-5 and by the Table 5 script."""
    global M, ctrl
    M = metrics
    ctrl = M["ctrl"]


def rel(name, q): return 100.0 * (M[name][q] / ctrl[q] - 1.0) if name in M else np.nan
def dif(name, q): return M[name][q] - ctrl[q] if name in M else np.nan


def fig1(out):
    """Part 1, Fig. 1: structure of M&M-CO2Wheat v1.0 (hourly and daily process bands)."""
    import matplotlib as _mpl
    with _mpl.rc_context(_mpl.rcParamsDefault):
        plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans"],
                             "pdf.fonttype": 42, "savefig.dpi": 400})
        _fig1_draw(out)


def _fig1_draw(OUT):
    INK = "#1a1a1a"; MUTE = "#6b6b6b"
    # two-hue palette: hourly physics in a single blue, daily crop processes in a single green; connectors in ink
    BLUE = "#1f4e79"      # hourly band and boxes
    ORANGE = "#2f6b3a"    # daily band and boxes (green)
    C_FORCE = C_RAD = C_GAS = C_EB = C_SOIL = BLUE
    C_PHEN = C_CAN = C_N = C_C = C_GRAIN = ORANGE
    def tint(c, f=0.90):
        r, g, b = int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)
        return "#%02x%02x%02x" % tuple(int(v + (255 - v) * f) for v in (r, g, b))

    fig = plt.figure(figsize=(17 * CM, 10.4 * CM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 170); ax.set_ylim(2, 98); ax.axis("off")

    BW, BH = 29.0, 17.5
    GAP = 5.0
    X0 = 3.0
    xs = [X0 + i * (BW + GAP) for i in range(5)]
    Y_H = 68.0
    Y_D = 22.0


    def box(x, y, title, lines, ec, w=BW, h=BH):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=tint(ec), ec=ec, lw=1.4))
        ax.add_patch(Rectangle((x + 0.9, y + h - 4.9), w - 1.8, 3.4, fc=ec, ec="none", zorder=2))
        ax.text(x + w / 2, y + h - 3.2, title, ha="center", va="center", fontsize=6.3, fontweight="bold", color="white", zorder=3)
        for i, l in enumerate(lines):
            ax.text(x + w / 2, y + h - 7.9 - 3.3 * i, l, ha="center", va="center", fontsize=5.6, color=INK)


    def arr(p, q, color=INK, lw=0.9):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=7, lw=lw, color=color, shrinkA=0, shrinkB=0))


    def lab(x, y, s, color=INK, ha="center", va="center", fs=5.8, bg=True):
        ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va=va,
                bbox=dict(fc="white", ec="none", pad=0.8) if bg else None, zorder=5)


    # ---------------- band frames and titles
    ax.add_patch(Rectangle((1, Y_H - 7.5), 168, BH + 16.5, fc="none", ec=BLUE, lw=1.0, ls=(0, (4, 3))))
    ax.text(3, Y_H + BH + 6.6, "HOURLY  ·  inside the daily coupling loop (up to five passes, until the committed canopy temperature converges)",
            fontsize=6.4, color=BLUE, fontweight="bold", va="center")
    ax.add_patch(Rectangle((1, Y_D - 5), 168, BH + 11, fc="none", ec=ORANGE, lw=1.0, ls=(0, (4, 3))))
    ax.text(3, Y_D + BH + 3.6, "DAILY  ·  crop processes advanced once per day after the loop has converged (call order left to right)",
            fontsize=6.4, color=ORANGE, fontweight="bold", va="center")

    # ---------------- hourly row
    box(xs[0], Y_H, "Forcing  (3.1)", ["daily weather → hourly", "T, RH, S↓, L↓, wind, rain", "solar geometry"], C_FORCE)
    box(xs[1], Y_H, "Radiation  (3.2–3.4)", ["sunlit / shaded leaf area", "absorbed PAR", "leaf-N profile"], C_RAD)
    box(xs[2], Y_H, "Gas exchange  (4–5)", ["FvCB + de novo N assimilation", "stomatal conductance with", "leaf water potential"], C_GAS)
    box(xs[3], Y_H, "Energy balance  (6)", ["canopy – air – soil network", "R$_\\mathrm{n}$, H, λE, G", "canopy and soil temperature"], C_EB)
    box(xs[4], Y_H, "Soil water, heat  (7)", ["Richards equation, 10 layers", "(48–240 sub-steps per day)", "heat conduction, hourly"], C_SOIL)
    ym = Y_H + BH / 2
    for i in range(4):
        arr((xs[i] + BW, ym), (xs[i + 1], ym))
    lab(xs[1] + BW + GAP / 2, ym + 2.2, "PAR, T", bg=True, fs=4.8)
    lab(xs[2] + BW + GAP / 2, ym + 2.2, "g$_\\mathrm{s}$, A$_\\mathrm{n}$", bg=True, fs=4.8)
    lab(xs[3] + BW + GAP / 2, ym + 2.2, "λE", bg=True, fs=4.8)
    # feedbacks (blue) routed below the hourly boxes
    yb1 = Y_H - 2.6
    x_from = xs[4] + BW * 0.40; x_to = xs[3] + BW * 0.60
    ax.plot([x_from, x_from, x_to, x_to], [Y_H, yb1, yb1, yb1], color=C_SOIL, lw=1.0)
    arr((x_to, yb1), (x_to, Y_H), color=C_SOIL, lw=1.0)
    lab((x_from + x_to) / 2, yb1, "committed T$_\\mathrm{r}$, E", color=C_SOIL)
    yb2 = Y_H - 5.6
    x_from = xs[4] + BW * 0.15; x_to = xs[2] + BW * 0.55
    ax.plot([x_from, x_from, x_to, x_to], [Y_H, yb2, yb2, yb2], color=C_SOIL, lw=1.0)
    arr((x_to, yb2), (x_to, Y_H), color=C_SOIL, lw=1.0)
    lab((x_from + x_to) / 2, yb2, "root-zone water potential ψ$_\\mathrm{s}$", color=C_SOIL)
    # canopy temperature back to the leaf, above the boxes
    yt = Y_H + BH + 2.6
    x_from = xs[3] + BW * 0.50; x_to = xs[2] + BW * 0.75
    ax.plot([x_from, x_from, x_to, x_to], [Y_H + BH, yt, yt, yt], color=C_EB, lw=1.0)
    arr((x_to, yt), (x_to, Y_H + BH), color=C_EB, lw=1.0)
    lab((x_from + x_to) / 2, yt, "T$_\\mathrm{can}$, leaf-to-air vapor deficit", color=C_EB)

    # ---------------- daily row
    box(xs[0], Y_D, "Phenology  (8)", ["thermal clock with", "vernalization and photoperiod", "eight stages"], C_PHEN)
    box(xs[1], Y_D, "Canopy  (9.3–9.4)", ["leaf appearance, expansion,", "senescence per rank and axis", "tillering, axis survival"], C_CAN)
    box(xs[2], Y_D, "Nitrogen  (10)", ["soil mineral N; uptake demand", "and regulation; nitrate reduction", "N payment; organ-N turnover"], C_N)
    box(xs[3], Y_D, "Carbon, growth  (9.1–9.5)", ["respiration, partitioning", "stem reserves and surplus", "structure only when N is paid"], C_C)
    box(xs[4], Y_D, "Spike, grain  (9.6, 11–12)", ["kernel number and filling", "grain nitrogen; minerals", "Fe, Zn, Cu, Mn, Mg, Ca"], C_GRAIN)
    yd = Y_D + BH / 2
    for i in range(4):
        arr((xs[i] + BW, yd), (xs[i + 1], yd), color=ORANGE)

    # ---------------- exchanges between the bands: straight vertical connectors in the gutter
    Ytop = Y_D + BH + 6.0          # top of the daily frame
    Ybot = Y_H - 7.5               # bottom of the hourly frame
    yl = (Ytop + Ybot) / 2
    # every label centered on its own arrow, over the arrow line (white box); two neighbours staggered in height
    x = xs[1] + BW * 0.5
    arr((x, Ytop), (x, Ybot), color=C_CAN, lw=1.0); lab(x, yl + 2.6, "LAI, leaf N (day t − 1)", color=C_CAN)
    x = xs[2] + BW * 0.25
    arr((x, Ybot), (x, Ytop), color=C_GAS, lw=1.0); lab(x, yl - 2.6, "daily A$_\\mathrm{n}$ and N-assimilation capacity", color=C_GAS)
    x = xs[4] + BW * 0.5
    arr((x, Ybot), (x, Ytop), color=C_SOIL, lw=1.0); lab(x, yl, "soil water, root uptake, T$_\\mathrm{soil}$", color=C_SOIL)
    x = xs[0] + BW * 0.5
    arr((x, Ybot), (x, Ytop), color=C_EB, lw=1.0); lab(x, yl, "daily-mean T$_\\mathrm{can}$", color=C_EB)

    ax.text(85, 8.0, "Outputs: hourly energy-balance and gas-exchange ledgers · daily organ, nitrogen, water, soil and grain states · yield, yield components and grain quality at harvest",
            ha="center", fontsize=5.9, color=MUTE)

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, "P1_Fig01_model_structure.%s" % ext), bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print("  wrote P1_Fig01_model_structure")


# ------------------------------------------------------------ Fig 2: drivers one at a time
def fig2(OUT, metrics=None):
    if metrics is not None: set_metrics(metrics)
    groups = [("CO$_2$", [("350", "co2_350"), ("550", "co2_550"), ("700", "co2_700"), ("850", "co2_850")]),
              ("Temp.", [("−2 K", "T-2"), ("+2 K", "T+2"), ("+4 K", "T+4")]),
              ("RH", [("−15 pts", "RH-15"), ("+15 pts", "RH+15")]),
              ("Rad.", [("×0.8", "RADx0.8"), ("×1.2", "RADx1.2")]),
              ("Irrig.", [("×0.5", "IRRx0.50"), ("×0.75", "IRRx0.75"), ("×1.25", "IRRx1.25"), ("×1.5", "IRRx1.50")]),
              ("N fert.", [("×0.25", "FERTx0.25"), ("×0.5", "FERTx0.50"), ("×1.5", "FERTx1.50"), ("×2", "FERTx2.00")])]
    outs = [("yield_g_m2", "Grain yield (%)", rel), ("AGB_g_m2", "Above-ground DM (%)", rel), ("ET_mm", "Evapotranspiration (%)", rel),
            ("Tr_mm", "Transpiration (%)", rel), ("grain_protein_pct", "Grain protein (%)", rel), ("dT_midday_K", "Midday canopy − air (K)", dif)]
    fig, axes = plt.subplots(3, 2, figsize=(W2, 0.95 * W2), sharex=True)
    labels, xs, seps, gx = [], [], [], []
    x = 0
    for gi, (gname, items) in enumerate(groups):
        for lab, _ in items: labels.append(lab); xs.append(x); x += 1
        seps.append(x - 0.5); gx.append((x - len(items) / 2 - 0.5, gname)); x += 0.6
    for k, (ax, (q, lab, fn)) in enumerate(zip(axes.flat, outs)):
        names = [name for _, items in groups for _, name in items]
        vals = [fn(name, q) for name in names]
        cols = ["white" if n in TRUNC else BAR for n in names]; hat = ["////" if n in TRUNC else "" for n in names]
        bars = ax.bar(xs, vals, 0.72, color=cols, edgecolor=BAR, linewidth=0.6)
        for b, h in zip(bars, hat): b.set_hatch(h)
        ax.axhline(0, color=INK, lw=0.6)
        for sx in seps[:-1]: ax.axvline(sx + 0.3, color="#cccccc", lw=0.5)
        ax.set_ylabel(lab); letter(ax, k)
        ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_xlim(-0.7, x - 0.6 + 0.1)
    for ax in axes[2]:
        for gxx, gname in gx: ax.text(gxx, -0.30, gname, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6, color=INK)
    fig.tight_layout(h_pad=0.8, w_pad=1.2, rect=(0, 0.03, 1, 1))
    save(fig, OUT, "P1_Fig02_driver_responses")

# ------------------------------------------------------------ Fig 3 alt: cascade (normalized) + elasticity chain
def fig3(OUT, metrics=None):
    if metrics is not None: set_metrics(metrics)
    levels = [(280, "co2_280"), (350, "co2_350"), (420, "ctrl"), (550, "co2_550"), (700, "co2_700"), (850, "co2_850")]
    # three groups: gas exchange (black), canopy and yield (grey), grain quality (red); symbols and fill distinguish
    # members within a group
    GX = [("An_ground_daylight_umol_m2_s", "Canopy assimilation", "#1a1a1a", "o", True, 3.8),
          ("gsw_day_mol_m2_s", "Stomatal conductance", "#1a1a1a", "s", True, 3.6),
          ("Tr_mm", "Transpiration", "#1a1a1a", "^", False, 4.0)]
    # grain yield is not drawn (it coincides with kernel number)
    GY = [("LAI_peak", "Peak LAI", "#7a7a7a", "D", True, 3.4),
          ("kernels_m2", "Kernel number", "#7a7a7a", "o", True, 3.6),
          ("WUE_yield_per_ET_g_mm", "Yield per unit ET", "#7a7a7a", "P", False, 4.2)]
    GQ = [("grain_N_pct", "Grain N concentration", C_SIM, "s", True, 3.6),
          ("grain_Zn_mg_kg", "Grain Zn concentration", C_SIM, "D", True, 3.4),
          ("grain_Fe_mg_kg", "Grain Fe concentration", C_SIM, "^", False, 4.0)]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(W2, 0.54 * W2), gridspec_kw=dict(width_ratios=[1.3, 1]))
    xv = [c for c, _ in levels]
    handles = {}
    for grp, items in (("Leaf gas exchange", GX), ("Canopy and yield", GY), ("Grain quality", GQ)):
        hs = []
        for q, lab, col, mk, filled, ms in items:
            yv = [100 * (M[n][q] / ctrl[q] - 1) for _, n in levels]
            ln, = ax.plot(xv, yv, marker=mk, color=col, ms=ms, lw=1.0, mfc=col if filled else "white", mew=0.9, label=lab,
                          zorder=3 if filled else 4)
            hs.append(ln)
        handles[grp] = hs
    ax.axhline(0, color=GRAY, lw=0.6); ax.axvline(420, color=GRAY, lw=0.6, ls=":")
    ax.set_xlabel("Atmospheric CO$_2$ (µmol mol$^{-1}$)"); ax.set_ylabel("Change from 420 µmol mol$^{-1}$ (%)")
    ax.set_xlim(260, 870); ax.margins(y=0.12)
    ax.text(0.03, 0.965, "(a)", transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    # three grouped legends below panel (a)
    x0 = 0.0
    for k, (grp, hs) in enumerate(handles.items()):
        lg = ax.legend(handles=hs, title=grp, fontsize=5.6, title_fontsize=6.0, loc="upper left", bbox_to_anchor=(x0, -0.2), frameon=False, handlelength=1.5, handletextpad=0.4, borderpad=0.2, labelspacing=0.25)
        lg._legend_box.align = "left"
        ax.add_artist(lg); x0 += 0.40 if k == 0 else 0.37
    ax.text(0.0, -0.60, "Grain yield follows kernel number exactly (the kernel-mass bound binds; Fig. 3b).", transform=ax.transAxes, ha="left", fontsize=5.6, color=GRAY)
    # (b) elasticities d ln Y / d ln CO2 between 420 and 550, ordered along the causal chain
    order = [("An_ground_daylight_umol_m2_s", "Assimilation"), ("gsw_day_mol_m2_s", "Conductance"), ("Tr_mm", "Transpiration"),
             ("dT_midday_K", "Canopy − air (K per ln CO$_2$)"), ("LAI_peak", "Peak LAI"), ("AGB_g_m2", "Above-ground DM"),
             ("kernels_m2", "Kernel number"), ("TKW_g", "Kernel mass"), ("yield_g_m2", "Grain yield"), ("N_uptake_gN_m2", "N uptake"),
             ("grain_N_pct", "Grain N concentration"), ("grain_Zn_mg_kg", "Grain Zn concentration"), ("grain_Fe_mg_kg", "Grain Fe concentration"), ("WUE_yield_per_ET_g_mm", "Yield per unit ET")]
    dl = math.log(550 / 420)
    vals = []
    for q, lab in order:
        if q == "dT_midday_K": vals.append((M["co2_550"][q] - ctrl[q]) / dl)
        else: vals.append(math.log(M["co2_550"][q] / ctrl[q]) / dl if ctrl[q] else 0.0)
    y = np.arange(len(order))[::-1]
    cols = [C_SIM if v > 0 else "#4d4d4d" for v in vals]
    bx.barh(y, vals, color=cols, height=0.62)
    for yi, v, (q, lab) in zip(y, vals, order):
        bx.text(v + (0.02 if v >= 0 else -0.02), yi, "%+.2f" % v, va="center", ha="left" if v >= 0 else "right", fontsize=6)
    bx.axvline(0, color=INK, lw=0.6); bx.set_yticks(y); bx.set_yticklabels([l for _, l in order], fontsize=6.5)
    bx.set_xlabel("Elasticity to CO$_2$, d ln Y / d ln CO$_2$ (420 → 550)"); bx.set_xlim(-1.3, 1.4)
    bx.text(0.03, 0.965, "(b)", transform=bx.transAxes, fontsize=8, fontweight="bold", va="top")
    fig.tight_layout(w_pad=2.0)
    save(fig, OUT, "P1_Fig03_alt_cascade_elasticity")

# ------------------------------------------------------------ Fig 4: 550 vs ambient across water, N, temperature
def fig4(OUT, metrics=None):
    if metrics is not None: set_metrics(metrics)
    panels = [("Water supply (× irrigation)", [(0.5, "IRRx0.50", "co2_550_IRRx0.50"), (0.75, "IRRx0.75", "co2_550_IRRx0.75"), (1.0, "ctrl", "co2_550"), (1.25, "IRRx1.25", "co2_550_IRRx1.25")]),
              ("N fertilizer (× 212 kg N ha$^{-1}$)", [(0.25, "FERTx0.25", "co2_550_FERTx0.25"), (0.5, "FERTx0.50", "co2_550_FERTx0.50"), (1.0, "ctrl", "co2_550"), (1.5, "FERTx1.50", "co2_550_FERTx1.50")]),
              ("Air temperature change (K)", [(-2, "T-2", "co2_550_T-2"), (0, "ctrl", "co2_550"), (2, "T+2", "co2_550_T+2"), (4, "T+4", "co2_550_T+4")])]
    outs = [("yield_g_m2", "Grain yield", "o"), ("AGB_g_m2", "Above-ground DM", "s"), ("ET_mm", "Evapotranspiration", "^"), ("grain_N_pct", "Grain N", "D")]
    fig, axes = plt.subplots(1, 3, figsize=(W2, 0.36 * W2), sharey=True)
    cols = {"yield_g_m2": "#1a1a1a", "AGB_g_m2": "#7f7f7f", "ET_mm": "#b0b0b0", "grain_N_pct": "#4d4d4d"}
    for k, (ax, (xl, items)) in enumerate(zip(axes, panels)):
        ax.axhspan(10, 20, color="#e8e8e8", lw=0, zorder=0)
        if k == 0: ax.text(0.5, 0.955, "shaded: range of mean wheat FACE yield\nresponses at ~550 µmol mol$^{-1}$ (+10 to +20 %)", transform=ax.transAxes, fontsize=5.4, color=GRAY, ha="center", va="top")
        for q, lab, mk in outs:
            ok = [(x, a, e) for x, a, e in items if a in M and e in M and M[a][q]]; xv = [x for x, a, e in ok]; yv = [100 * (M[e][q] / M[a][q] - 1) for x, a, e in ok]
            ax.plot(xv, yv, mk + "-", color=cols[q], ms=3.6, lw=0.9, mfc="white" if q in ("ET_mm",) else cols[q], label=lab)
            for x, a, e in ok:
                if a in TRUNC or e in TRUNC: ax.plot([x], [100 * (M[e][q] / M[a][q] - 1)], mk, color=cols[q], ms=3.6, mfc="white", mew=0.9)
        ax.axhline(0, color=INK, lw=0.6); ax.set_xlabel(xl); letter(ax, k)
        if k == 0: ax.set_ylabel("Effect of 550 vs 420 µmol mol$^{-1}$ CO$_2$ (%)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=6.5, frameon=True, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(w_pad=1.5, rect=(0, 0.08, 1, 1))
    save(fig, OUT, "P1_Fig04_CO2_interactions")


def _vik():
    """The diverging 'vik' colour map of Crameri (cmcrameri package)."""
    from cmcrameri import cm as cmc
    return cmc.vik


# ------------------------------------------------------------ Fig 5 alt: sensitivity matrix + yield tornado
PARAMS = [("fruiting_efficiency", "Fruiting efficiency"), ("TKW_potential", "Potential kernel mass"), ("plant_density", "Plant density"),
          ("N_mineral_rate", "Soil N mineralization rate"), ("leaf_width", "Leaf width"), ("SLN_min", "Minimum SLN"),
          ("u_max_fw", "Root uptake u$_{max}$"), ("N_max_ref", "Nitrate assimilation N$_{max}$"),
          ("N_soil_init", "Initial soil mineral N")]
SOUT = [("yield_g_m2", "Yield"), ("kernels_m2", "Kernels"), ("AGB_g_m2", "AGB"), ("LAI_peak", "LAI"), ("ET_mm", "ET"), ("grain_N_pct", "Grain N"), ("N_uptake_gN_m2", "N uptake"), ("grain_Zn_mg_kg", "Grain Zn"), ("grain_Fe_mg_kg", "Grain Fe"), ("dT_midday_K", "Canopy − air")]
def S(pn, q, f):
    n = "P_%s_x%.1f" % (pn, f)
    if n not in M: return np.nan
    if q == "dT_midday_K": return (M[n][q] - ctrl[q]) / (f - 1)          # K per unit relative change
    return ((M[n][q] - ctrl[q]) / ctrl[q]) / (f - 1)
def fig5(OUT, metrics=None):
    if metrics is not None: set_metrics(metrics)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(W2, 0.5 * W2), gridspec_kw=dict(width_ratios=[1.55, 1]))
    Z = np.array([[np.nanmean([S(pn, q, 0.8), S(pn, q, 1.2)]) for q, _ in SOUT] for pn, _ in PARAMS])
    disp = Z.copy(); disp[:, SOUT.index(("dT_midday_K", "Canopy − air"))] *= 1.0
    norm = TwoSlopeNorm(vmin=-1.3, vcenter=0, vmax=1.3)
    im = ax.imshow(np.clip(disp, -1.3, 1.3), cmap=_vik(), norm=norm, aspect="auto")
    for i in range(len(PARAMS)):
        for j in range(len(SOUT)):
            v = Z[i, j]
            if abs(v) >= 0.05: ax.text(j, i, "%+.2f" % v, ha="center", va="center", fontsize=5.8, color="white" if abs(v) > 0.8 else INK)
    ax.set_xticks(range(len(SOUT))); ax.set_xticklabels([l for _, l in SOUT], rotation=40, ha="right", rotation_mode="anchor", fontsize=6.5)
    ax.set_yticks(range(len(PARAMS))); ax.set_yticklabels([l for _, l in PARAMS], fontsize=6.5); ax.tick_params(length=0)
    for k in range(len(SOUT) - 1): ax.axvline(k + 0.5, color="white", lw=0.8)
    for k in range(len(PARAMS) - 1): ax.axhline(k + 0.5, color="white", lw=0.8)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02); cb.set_label("S$_{Y,P}$ (mean of ±20 %)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    ax.set_title("Normalized sensitivity index", fontsize=7.5, pad=4); ax.text(-0.42, 1.03, "(a)", transform=ax.transAxes, fontsize=8, fontweight="bold")
    # (b) tornado for yield
    rows = []
    for pn, lab in PARAMS:
        lo = rel("P_%s_x0.8" % pn, "yield_g_m2"); hi = rel("P_%s_x1.2" % pn, "yield_g_m2")
        rows.append((lab, lo, hi))
    rows.sort(key=lambda r: abs(r[2] - r[1]))
    y = np.arange(len(rows))
    for yi, (lab, lo, hi) in zip(y, rows):
        # same sign: draw the longer bar first so the shorter one stays visible on top
        pairs = [(hi, C_SIM, "+20 %"), (lo, "#4d4d4d", "−20 %")]
        pairs.sort(key=lambda t: -abs(t[0]))
        for v, c, l in pairs:
            bx.barh(yi, v, color=c, height=0.6, label=l if yi == 0 else None, zorder=3 if abs(v) == min(abs(hi), abs(lo)) else 2)
        bx.text(max(hi, 0) + 0.6, yi, "%+.1f / %+.1f" % (lo, hi), ha="left", va="center", fontsize=5.6, color=INK)
    bx.set_yticks(y); bx.set_yticklabels([r[0] for r in rows], fontsize=6.5); bx.axvline(0, color=INK, lw=0.6)
    bx.set_xlabel("Change in grain yield (%)"); bx.set_xlim(-24, 44); bx.legend(fontsize=6, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.16), frameon=True)
    bx.set_title("Yield, ±20 % per parameter", fontsize=7.5, pad=4); bx.text(-0.55, 1.03, "(b)", transform=bx.transAxes, fontsize=8, fontweight="bold")
    fig.tight_layout(w_pad=2.5)
    save(fig, OUT, "P1_Fig05_alt_matrix_tornado")

