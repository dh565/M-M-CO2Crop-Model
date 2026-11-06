#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Part 2 figures of M&M-CO2Wheat v1.0 at the three FACE experiments: shared style, observation
access and the figure functions (Figs. 2, 4-9 and 11 of Part 2; Fig. S6 of Supplement S2).

Every model number is read from the reference-outputs archive through lib/face_results.py and
lib/ebwb_gx.py; init(archive) loads the archive once.  Output: PDF + PNG (300 dpi), GMD widths.
"""
import csv
import datetime as dt
import math
import os
import sys
from collections import OrderedDict, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import face_results as P2
import ebwb_gx as EB
import tleaf_eval as TE
import maricopa_gx as MGX
import maricopa_archive as MA

_EB_pairs = EB.pairs


def _pairs_all_sites(deck, qkey):
    """Maricopa leaf gas exchange (Wall et al., 2000, 2006) joins the WheatDryFACE campaigns."""
    if deck.startswith("MARICOPA_") and qkey in ("An", "gs"):
        return MGX.pairs(deck, qkey)
    return _EB_pairs(deck, qkey)


EB.pairs = _pairs_all_sites

NAN = float("nan")
isnan = lambda x: x is None or (isinstance(x, float) and math.isnan(x))

# ----------------------------------------------------------------- style
plt.rcParams.update({
    # style after GMD evaluation papers: Helvetica/Arial, full box frame, black observations and red simulations,
    # bold centered panel titles, framed legends
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold", "mathtext.default": "regular",
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "axes.titleweight": "bold", "legend.fontsize": 6.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "figure.titlesize": 13, "axes.formatter.use_mathtext": True,
    "axes.linewidth": 0.7, "axes.edgecolor": "black", "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.8, "ytick.major.size": 2.8, "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": True, "axes.spines.right": True,
    "lines.linewidth": 1.0, "lines.markersize": 4,
    "legend.frameon": True, "legend.fancybox": False, "legend.edgecolor": "black", "legend.framealpha": 1.0,
    "legend.handlelength": 1.8, "legend.borderpad": 0.4,
    "patch.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
})
CM = 1 / 2.54
W1, W15, W2 = 8.3 * CM, 12.0 * CM, 17.0 * CM          # GMD column widths

C_A, C_E = "#0072B2", "#D55E00"                         # aCO2 blue / eCO2 orange (Okabe-Ito)
C_OBS = "black"
C_MOD = C_SIM = "#C8102E"                               # simulated where CO2 is not the contrast (obs black, sim red)
SITE_C = {"WDF25": "#7B2D8E", "WDF23": "#7B2D8E", "BRAUNSCHWEIG": "#2CA02C", "MARICOPA": "#D4A017"}   # one color per site
SITE_M = {"WDF25": "o", "WDF23": "o", "BRAUNSCHWEIG": "^", "MARICOPA": "D"}   # one symbol per site (seasons pooled)
LEG_SITES = ["WDF25", "BRAUNSCHWEIG", "MARICOPA"]
LEG_N = {"WDF25": "WheatDryFACE 2022/23, 2024/25", "BRAUNSCHWEIG": "Braunschweig 2013/14, 2014/15", "MARICOPA": "Maricopa 1993–97"}
SITE_N = {"WDF25": "WheatDryFACE 2024/25", "WDF23": "WheatDryFACE 2022/23",
          "BRAUNSCHWEIG": "Braunschweig 2013/14, 2014/15", "MARICOPA": "Maricopa 1993–97"}
LS_W = {"WW": "-", "WL": "--"}
GRAY = "#8a8a8a"

DECKS = P2.DECKS
FACE = [d for d in P2.FACE_DECKS]


def d2d(s):
    return dt.date.fromisoformat(s)


def datefmt(ax, interval=None):
    lo, hi = ax.get_xlim()
    if interval is None and (hi - lo) < 75:                      # short span: fortnightly day ticks
        ax.xaxis.set_major_locator(mdates.DayLocator(bymonthday=(1, 15)))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax.tick_params(axis="x", labelsize=6.5); return
    if interval is None:
        interval = 1 if (hi - lo) < 200 else (2 if (hi - lo) < 300 else 3)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=interval))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.tick_params(axis="x", labelsize=6.5)


def letter(ax, k, x=None, y=None, inside=True):
    """Panel letter: bold "(a)" inside the frame at top-left (x, y kept for call compatibility)."""
    if inside:
        ax.text(0.03, 0.965, "(%s)" % "abcdefghijklmnopqrstuvwxyz"[k], transform=ax.transAxes, fontsize=8,
                fontweight="bold", va="top", ha="left", zorder=6,
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))
    else:
        ax.text(x if x is not None else -0.14, y if y is not None else 1.04, "(%s)" % "abcdefghijklmnopqrstuvwxyz"[k],
                transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def save(fig, out, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, "%s.%s" % (stem, ext)), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("  wrote", stem)


def onetoone(ax, x, y, unit=""):
    xs = [v for v in list(x) + list(y) if not isnan(v)]
    if not xs:
        return
    lo, hi = min(xs), max(xs)
    pad = 0.06 * (hi - lo if hi > lo else abs(hi) or 1)
    lo, hi = lo - pad, hi + pad
    if lo > 0 and lo < 0.25 * hi:
        lo = 0
    ax.plot([lo, hi], [lo, hi], color="black", lw=0.5, zorder=0)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    # identical tick positions on both axes of a 1:1 panel
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))


def stats_txt(x, y):
    p = [(a, b) for a, b in zip(x, y) if not isnan(a) and not isnan(b)]
    if len(p) < 2:
        return "n = %d" % len(p)
    o = np.array([a for a, b in p]); m = np.array([b for a, b in p])
    bias = (m - o).mean(); rmse = math.sqrt(((m - o) ** 2).mean())
    r = np.corrcoef(o, m)[0, 1] if o.std() > 0 and m.std() > 0 else NAN
    om = o.mean()
    rel = " (%.0f %%)" % (rmse / om * 100) if om > 0 and (o.min() >= 0) and rmse / om < 1.5 else ""
    return "n = %d, r = %.2f\nRMSD = %s%s" % (len(p), r, fmt(rmse), rel)


def place_stats(ax, X, Y, txt=None, XE=None, fs=5.6):
    """Put the statistics text in the emptier off-diagonal corner (top-left under the panel letter, or bottom-right),
    counting plotted points and their x error bars in the two candidate boxes."""
    txt = txt if txt is not None else stats_txt(X, Y)
    def nrm(v, a, b): return (v - a) / (b - a) if b > a else 0.5
    def count():
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        tl = br = 0
        for i, (o, m) in enumerate(zip(X, Y)):
            if isnan(o) or isnan(m): continue
            e = XE[i] if (XE is not None and i < len(XE) and not isnan(XE[i])) else 0.0
            ny = nrm(m, y0, y1)
            for ox in (o - e, o, o + e):
                nx = nrm(ox, x0, x1)
                if nx < 0.76 and ny > 0.64: tl += 1
                if nx > 0.28 and ny < 0.30: br += 1
        return tl, br
    tl, br = count()
    for _ in range(3):                      # both corners occupied: widen the square range so the cloud shrinks to the lower left
        if min(tl, br) == 0: break
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        lo, hi = min(x0, y0), max(x1, y1); hi = lo + (hi - lo) * 1.25
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        for l in ax.get_lines():
            if len(l.get_xdata()) == 2 and l.get_zorder() == 0:
                l.set_data([lo, hi], [lo, hi])
        tl, br = count()
    if tl <= br:
        ax.text(0.04, 0.865, txt, transform=ax.transAxes, va="top", ha="left", fontsize=fs, zorder=6, bbox=dict(facecolor="white", edgecolor="none", pad=0.4, alpha=0.85))
    else:
        ax.text(0.97, 0.03, txt, transform=ax.transAxes, va="bottom", ha="right", fontsize=fs, zorder=6, bbox=dict(facecolor="white", edgecolor="none", pad=0.4, alpha=0.85))


def fmt(v):
    if isnan(v):
        return "–"
    a = abs(v)
    if a >= 100: return "%.0f" % v
    if a >= 10: return "%.1f" % v
    if a >= 1: return "%.2f" % v
    return "%.3f" % v



# ----------------------------------------------------------------- init
def init(pkg):
    """Load the archive: model scalars and series of every deck, the observation registry and the contrast rows."""
    P2.PKG = os.path.abspath(pkg)
    sys.stderr.write("loading %s ...\n" % pkg)
    P2.MODELS = P2.load_arm(P2.PKG, FACE)
    P2.OBS_REG = P2.Obs()
    P2.PHEN = {d: P2.phenology_rows(P2.PKG, d) for d in FACE}
    P2.SWC = {d: P2.swc_pairs(P2.PKG, P2.OBS_REG, d) for d in FACE}
    P2.CAPS = {d: P2.grain_fill_law(P2.PKG, d) for d in FACE}
    P2.EARS = {d: P2.ear_number(P2.PKG, d) for d in FACE}
    P2.GRAIN_MIN = {d: P2.grain_minerals(P2.PKG, d) for d in FACE}
    P2.CROPN = {d: P2.crop_N_derived(P2.PKG, d) for d in FACE}
    P2.AGB_ANTH = {d: P2.agb_at_dated_sample(P2.PKG, d, P2.OBS_REG) for d in FACE}
    P2.MAR_EARS = {}
    lines = [l for l in open(os.path.join(config.MARICOPA_OBS_DIR, "MARICOPA_YECORA_ROJO_MEASURED.csv"))
             if not l.startswith("#")]
    for r in csv.DictReader(lines):
        dk = "MARICOPA_%s" % r["TRNO"]
        if dk in DECKS:
            P2.MAR_EARS[dk] = P2.fl(r["ears_m2_PnAM"])
    P2.ORGAN_N_OBS = {}
    for d in FACE:
        s = DECKS[d]["site"]
        if s == "WDF25":
            cv, co2, wat = d.split("_")[1:4]
            for r in P2.W25.obs_for(cv, co2, wat):
                if r["variable"] == "[N]":
                    P2.ORGAN_N_OBS[(d, (r["organ"] or "").strip())] = (P2.fl(r["value"]), P2.fl(r["std"]))
        elif s == "WDF23":
            cv, co2 = d.split("_")[1:3]
            for r in P2.W23.obs_for(cv, co2):
                if r["variable"] == "[N]":
                    P2.ORGAN_N_OBS[(d, (r["organ"] or "").strip())] = (P2.fl(r["value"]), P2.fl(r["std"]))
        elif s == "BRAUNSCHWEIG":
            rr = P2.rows_of(P2.bra_obs_path(d))
            dates = sorted({x["date"] for x in rr if x["variable"] == "SUMDW"})
            if dates:
                hd = dates[-1]
                for x in rr:
                    if x["date"] == hd and x["variable"].startswith("N% "):
                        P2.ORGAN_N_OBS[(d, x["variable"])] = (P2.fl(x["mean"]), P2.fl(x["sd"]))
    # Maricopa grain protein from the archive quality workbook (Kimball et al., 2001)
    _mar_quality_into_registry()
    # treatment contrasts: the Braunschweig N1/N3 pairs of Fig. 9b and the CO2 pairs of both seasons in the order of Fig. 7c
    P2.CONTRASTS = [c for c in P2.CONTRASTS if not (c[1] == "N" and c[2] == "BRAUNSCHWEIG")] + [
        ("Braunschweig N1/N3 %s" % lab, "N", "BRAUNSCHWEIG", "BRAUNSCHWEIG_" + n, "BRAUNSCHWEIG_" + h)
        for lab, n, h in (("aCO2 2014", "101", "103"), ("eCO2 2014", "104", "106"), ("aCO2 2015", "113", "115"), ("eCO2 2015", "116", "118"))]
    if not any(c[3] == "BRAUNSCHWEIG_116" and c[1] == "eCO2" for c in P2.CONTRASTS):
        P2.CONTRASTS = [c for c in P2.CONTRASTS if not (c[1] == "eCO2" and c[2] == "BRAUNSCHWEIG")] + [
            ("Braunschweig eCO2/aCO2 %s" % nl, "eCO2", "BRAUNSCHWEIG", "BRAUNSCHWEIG_" + ea, "BRAUNSCHWEIG_" + aa)
            for nl, ea, aa in (("N1 2014", "104", "101"), ("N2 2014", "105", "102"), ("N3 2014", "106", "103"), ("N4 2014", "110", "107"),
                               ("N1 2015", "116", "113"), ("N2 2015", "117", "114"), ("N3 2015", "118", "115"), ("N4 2015", "122", "119"))]
    # one observation source for Figs. 6-11 and the tables: the registry used by the contrast tally is synchronised
    # with obs() (archive harvested yield and HI at Maricopa, consolidated-table peak LAI, kernel number with sd)
    for d in FACE:
        for q_ in P2.CQ:
            if q_.startswith("grain_") and q_[6:] in P2.MINERALS: continue
            v_, s_ = obs(d, q_)
            if not isnan(v_):
                o_ = P2.OBS_REG.get(d, q_)
                if (not o_) or isnan(o_["value"]) or abs(o_["value"] - v_) > 1e-9 * max(1.0, abs(v_)):
                    P2.OBS_REG.put(d, q_, v_, s_, "harvest", "synchronised with the figure observation source (obs())")
    P2.CROWS = P2.contrast_rows(P2.OBS_REG)
    _sign_by_ttest()
    sys.stderr.write("  %d decks, %d contrast rows\n" % (len(P2.MODELS), len(P2.CROWS)))
    # grain Ca is computed by the mineral module and carries a response sign; the Braunschweig grain Fe/Zn
    # observations live in the consolidated observation table (as drawn in Fig. 8) and are scored too
    for r in P2.CROWS:
        q = r["quantity"]
        if q == "grain_Ca" and r.get("prescribed", "").startswith("model grain Ca is a CONSTANT"):
            r["cls"], r["resolved"] = P2.sign_class(r["R_model"], r["R_obs"], r["R_obs_sd"]); r["prescribed"] = ""
        if r["site"] == "BRAUNSCHWEIG" and q in ("grain_Fe", "grain_Zn") and isnan(r["R_obs"]):
            ro, sd, rm = _ratio_rows2(r["num"], r["den"], q)
            if not isnan(ro) and not isnan(rm):
                r["R_obs"], r["R_obs_sd"], r["R_model"] = ro, sd, rm
                r["cls"], r["resolved"] = P2.sign_class(rm, ro, sd); r["prescribed"] = ""


MAR_QMAP = {'MARICOPA_905': (1994, 1), 'MARICOPA_906': (1994, 2), 'MARICOPA_907': (1994, 3), 'MARICOPA_908': (1994, 4),
            'MARICOPA_910': (1996, 6), 'MARICOPA_912': (1996, 8),
            'MARICOPA_913': (1997, 5), 'MARICOPA_914': (1997, 6), 'MARICOPA_915': (1997, 7), 'MARICOPA_916': (1997, 8)}


def _mar_quality_into_registry():
    """Maricopa grain protein from the archive quality workbook (Kimball et al., 2001; whole-grain protein,
    converted to dry basis); grain N = protein / 5.7."""
    try:
        x = MA.quality_sheet().copy()
    except Exception as e:
        sys.stderr.write('Maricopa quality workbook not read: %s\n' % e); return
    x.columns = ['LOC', 'YR', 'T', 'C', 'I', 'N', 'R', 'TWT', 'TWTkg', 'HRD', 'WPROT', 'WMIST']
    x = x.dropna(subset=['YR'])
    mst = x.WMIST.astype(float).where(x.WMIST.astype(float) > 0, 9.55)   # 1995/96-1996/97 moisture not recorded: 1993-94 mean
    x['PDRY'] = x.WPROT.astype(float) / (1 - mst / 100.0)
    for dk, (yr, t) in MAR_QMAP.items():
        g = x[(x.YR == yr) & (x['T'] == t)].PDRY
        if len(g) == 0: continue
        m, sd = float(g.mean()), (float(g.std(ddof=1)) if len(g) > 1 else NAN)
        src = 'Maricopa FACE quality data (Kimball et al., 2001), whole-grain protein on dry basis, n=%d' % len(g)
        P2.OBS_REG.put(dk, 'grain_protein_pct', m, sd, 'harvest', src)
        P2.OBS_REG.put(dk, 'grain_N_pct', m / 5.7, sd / 5.7 if sd == sd else NAN, 'harvest', src + '; N = protein / 5.7')


NREP = {'WDF25': 3, 'WDF23': 3, 'BRAUNSCHWEIG': 3, 'MARICOPA': 4}


def _sign_by_ttest():
    """An observed contrast counts as resolved when it is statistically significant (two-sided t test on the
    treatment means, P < 0.05, n replicate plots per treatment), the criterion the 95 % intervals of Figs. 7-9 express."""
    from scipy import stats
    for r in P2.CROWS:
        if r.get('cls') == 'C' and (r.get('prescribed') or isnan(r.get('R_obs', NAN))): continue
        rm, ro, sd = r['R_model'], r['R_obs'], r['R_obs_sd']
        if isnan(ro) or isnan(rm): continue
        n = NREP.get(r['site'], 3); t = stats.t.ppf(0.975, 2 * n - 2)
        sig = (abs(ro - 1) > t * sd / math.sqrt(n)) if (not isnan(sd) and sd > 0) else (abs(ro - 1) >= P2.SIGN_RESOLVE_PCT)
        same = (rm - 1) * (ro - 1) >= 0
        r['resolved'] = sig
        if sig and abs(rm - 1) < P2.SIGN_FLOOR and abs(rm - 1) < 0.5 * abs(ro - 1): r['cls'] = 'N'   # significant observed effect, simulated < 5 % and < half of it
        else: r['cls'] = 'agree' if same else ('A' if sig else 'B')


def obs(deck, q):
    o = P2.OBS_REG.get(deck, q)
    if q == "TKW" and DECKS[deck]["site"] == "WDF25":
        # a TKW row of the treatment file (with its replicate sd) takes precedence over the registry value
        cv, co2, wat = deck.split("_")[1:4]
        for r in P2.W25.obs_for(cv, co2, wat):
            if r["variable"] == "TKW" and r["value"] is not None:
                return r["value"], (r["std"] if r["std"] is not None else NAN)
        if o and not isnan(o["value"]):
            return o["value"], o.get("sd", NAN)
        return NAN, NAN
    if q == "kernel_number" and DECKS[deck]["site"] == "WDF25":
        # derived from the TKW returned above; sd in quadrature
        gy, gs = obs(deck, "grain_yield"); tk, ts = obs(deck, "TKW")
        if not isnan(gy) and not isnan(tk) and tk > 0:
            kn = gy / tk * 1000.0
            rel2 = sum((sd_ / v) ** 2 for v, sd_ in ((gy, gs), (tk, ts)) if not isnan(sd_) and v > 0)
            return kn, (kn * math.sqrt(rel2) if rel2 > 0 else NAN)
    if DECKS[deck]["site"] == "MARICOPA" and q in ("grain_yield", "HI"):
        # archive values: harvested yield HWAH (kg/ha -> g/m2) and the archive harvest index HIAM
        tr = deck.split("_")[1]
        for r in csv.DictReader(open(os.path.join(config.MARICOPA_OBS_DIR, "maricopa_archive_summary.csv"))):
            if r["TRNO"] == tr:
                if q == "grain_yield":
                    return float(r["HWAH"]) / 10.0, (float(r["HWAH_sd"]) / 10.0 if r["HWAH_sd"] else NAN)
                return (float(r["HIAM"]) if r["HIAM"] else NAN), NAN
    if not o:
        return _cv_obs(deck, q)
    return o["value"], o["sd"]


_CV = {}


def _cv_obs(deck, q):
    """Fallback: observation from the consolidated observation table (e.g. Braunschweig grain Fe/Zn,
    Maricopa peak LAI), which the loaders above do not carry."""
    if not _CV:
        try:
            for r in P2.table_rows(_consolidated_csv()):
                try: v = float(r["obs"])
                except (ValueError, TypeError): continue
                sd = NAN
                try: sd = float(r["obs_sd"]) if r["obs_sd"] not in ("", None) else NAN
                except ValueError: pass
                _CV.setdefault((r["deck"], r["quantity"]), (v, sd))
        except Exception as e:
            print("consolidated view unavailable:", e); _CV["__none__"] = None
    return _CV.get((deck, q), (NAN, NAN))


def mod(deck, q):
    m = P2.MODELS.get(deck)
    if not m:
        return NAN
    if q.startswith("grain_") and q[6:] in P2.MINERALS:
        return P2.GRAIN_MIN[deck].get(q[6:], NAN)
    return m.v.get(q, NAN)


def fig02(out):
    """Soil water content and soil temperature, time series at representative decks + pooled 1:1 all sites."""
    fig = plt.figure(figsize=(W2, 0.84 * W2))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.05], hspace=0.62, wspace=0.5, bottom=0.12)
    # rows 1-2 (site order: WDF 2022/23, WDF 2024/25, Braunschweig, Maricopa)
    panels = [("WDF23_Ruta_aCO2", "SWC20"), ("WDF23_Ruta_aCO2", "SWC55"), ("WDF25_Ruta_aCO2_WW", "SWC20"), ("WDF25_Ruta_aCO2_WL", "SWC20"),
              ("WDF25_Ruta_aCO2_WW", "SWC55"), ("WDF25_Ruta_aCO2_WL", "SWC55"), ("BRAUNSCHWEIG_103", "d0_20"), ("MARICOPA_906", "MAR_SWC_0_30")]
    titles = ["WDF 2022/23, 20 cm", "WDF 2022/23, 55 cm", "WDF 2024/25 WW, 20 cm", "WDF 2024/25 WL, 20 cm",
              "WDF 2024/25 WW, 55 cm", "WDF 2024/25 WL, 55 cm", "Braunschweig 2013/14\n320 kg N ha$^{-1}$, 0–20 cm", "Maricopa 1993/94\nwell-watered, 0–30 cm"]
    k = 0
    for j, (d, q) in enumerate(panels):
        ax = fig.add_subplot(gs[j // 4, j % 4])
        _swc_course(ax, d, q, legend=False)
        if d.startswith("WDF23") or d.startswith("MARICOPA"): datefmt(ax, interval=2)
        ax.set_title(titles[j], pad=3, fontsize=7)
        if j % 4 == 0: ax.set_ylabel("SWC (m$^3$ m$^{-3}$)")
        letter(ax, k, x=-0.25); k += 1
    # row 3: pooled 1:1 SWC all sites, Tsoil 20, Tsoil 55, and WDF Tsoil course
    ax = fig.add_subplot(gs[2, 0])
    XY = {}
    for site, qs in (("MARICOPA", ["MAR_SWC_0_30", "MAR_SWC_30_50"]), ("BRAUNSCHWEIG", ["d0_20", "d20_40"]),
                     ("WDF25", ["SWC20", "SWC55"]), ("WDF23", ["SWC20", "SWC55"])):
        X, Y = [], []
        for d in [x for x in FACE if DECKS[x]["site"] == site]:
            for q in qs:
                pp = EB.pairs(d, q)
                X += [p[1] for p in pp]; Y += [p[3] for p in pp]
        z = {"MARICOPA": 2, "BRAUNSCHWEIG": 5, "WDF25": 4, "WDF23": 4}[site]   # Braunschweig drawn on top
        ax.plot(X, Y, SITE_M[site], ms=2.0, mfc="none", mec=SITE_C[site], mew=0.45, alpha=0.55 if site == "MARICOPA" else 0.8, zorder=z, label=LEG_N.get(site, "_nolegend_"))
        key = "WDF25" if site.startswith("WDF") else site; XY.setdefault(key, ([], []))[0].extend(X); XY[key][1].extend(Y)
    onetoone(ax, *_pool_xy(ax))
    STATS_I = {}
    for key, nm in (("WDF25", "WheatDryFACE"), ("BRAUNSCHWEIG", "Braunschweig"), ("MARICOPA", "Maricopa")):
        x, y = np.array(XY[key][0]), np.array(XY[key][1])
        if len(x) > 2: STATS_I[nm] = (len(x), float(np.corrcoef(x, y)[0, 1]), float(np.sqrt(np.mean((y - x) ** 2))))
    print("  Fig2i per-site SWC stats (n, r, RMSE):", STATS_I)
    ax.set_xlabel("Observed SWC (m$^3$ m$^{-3}$)"); ax.set_ylabel("Simulated SWC (m$^3$ m$^{-3}$)")
    ax.set_title("All sites, all depths", pad=3); letter(ax, k, x=-0.25); k += 1
    hh, ll = ax.get_legend_handles_labels(); fig.legend(hh[::-1], ll[::-1], loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.015), fontsize=5.8, frameon=False)
    # (j, k) soil temperature, both sites that measured it: WheatDryFACE 20 / 55 cm and Maricopa 10 / 40 cm
    try:
        import maricopa_tsoil as _MT
        _MO = _MT.obs_daily(); _MM = {t: _MT.model_daily(t) for t in _MO}
    except Exception as _e:
        print("Maricopa Tsoil not available for Fig 2:", _e); _MO = {}; _MM = {}
    for j, (q, zmar, ttl) in enumerate([("Tsoil20", 0.10, "$T_{soil}$, upper sensor\n(WDF 20 cm, Maricopa 10 cm)"), ("Tsoil55", 0.40, "$T_{soil}$, deeper sensor\n(WDF 55 cm, Maricopa 40 cm)")]):
        ax = fig.add_subplot(gs[2, 1 + j])
        X, Y = [], []; Xw, Yw = [], []; Xm, Ym = [], []
        for d in [x for x in FACE if DECKS[x]["site"] in ("WDF25", "WDF23")]:
            pp = EB.pairs(d, q)
            o = [p[1] for p in pp]; m = [p[3] for p in pp]; Xw += o; Yw += m
        for t in _MO:
            oo = _MO[t][_MO[t].z == zmar]; mm_ = _MM[t][_MM[t].z == zmar]
            jn = oo.merge(mm_, on="date", suffixes=("_o", "_m")); Xm += list(jn.T_o); Ym += list(jn.T_m)
        ax.plot(Xm, Ym, SITE_M["MARICOPA"], ms=2.2, mfc="none", mec=SITE_C["MARICOPA"], mew=0.45, alpha=0.6, zorder=2)
        ax.plot(Xw, Yw, SITE_M["WDF25"], ms=2.2, mfc="none", mec=SITE_C["WDF25"], mew=0.5, alpha=0.8, zorder=3)
        X = Xw + Xm; Y = Yw + Ym
        onetoone(ax, X, Y)
        def _st2(x, y):
            x = np.array(x, float); y = np.array(y, float)
            return "r %.2f, RMSD %.1f K, n %d" % (np.corrcoef(x, y)[0, 1], np.sqrt(np.mean((y - x) ** 2)), len(x))
        ax.text(0.97, 0.03, "WDF: %s\nMaricopa: %s" % (_st2(Xw, Yw), _st2(Xm, Ym) if Xm else "no record"), transform=ax.transAxes, ha="right", va="bottom", fontsize=4.0, zorder=7, bbox=dict(facecolor="white", edgecolor="none", pad=0.5, alpha=0.9))
        print("  Fig2%s Tsoil stats WDF %s | Maricopa %s" % ("jk"[j], _st2(Xw, Yw), _st2(Xm, Ym) if Xm else "-"))
        ax.set_xlabel("Observed $T_{soil}$ (°C)"); ax.set_ylabel("Simulated $T_{soil}$ (°C)")
        ax.set_title(ttl, pad=3, fontsize=7); letter(ax, k, x=-0.25); k += 1
    ax = fig.add_subplot(gs[2, 3])
    d = "WDF25_Ruta_aCO2_WW"
    for q, ls in (("Tsoil20", "-"), ("Tsoil55", "--")):
        pp = EB.pairs(d, q)
        if pp:
            dates = [d2d(p[0]) for p in pp]
            ax.fill_between(dates, [p[1] - (p[2] if not isnan(p[2]) else 0) for p in pp], [p[1] + (p[2] if not isnan(p[2]) else 0) for p in pp],
                            color=C_OBS, alpha=0.15, lw=0)
            ax.plot(dates, [p[1] for p in pp], ls, color=C_OBS, lw=0.8, label=("observed ± sd (a–h; l: 20 cm)" if q == "Tsoil20" else "observed ± sd, 55 cm (l)"))
            ax.plot(dates, [p[3] for p in pp], ls, lw=0.9, color=C_SIM, label=("model (a–h; l: 20 cm)" if q == "Tsoil20" else "model, 55 cm (l)"))
    ax.set_title("$T_{soil}$ course, WDF 2024/25\nRuta, ambient CO$_2$, WW", pad=3, fontsize=7); ax.set_ylabel("$T_{soil}$ (°C)")
    datefmt(ax)
    hh, ll = ax.get_legend_handles_labels(); fig.legend(hh, ll, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.005), fontsize=5.8, frameon=False, handletextpad=0.3, columnspacing=1.2)
    letter(ax, k, x=-0.25)
    save(fig, out, "Fig02_soil_water_temperature")


def _pool_xy(ax):
    X, Y = [], []
    for l in ax.get_lines():
        X += list(l.get_xdata()); Y += list(l.get_ydata())
    return X, Y


def _swc_course(ax, d, q, legend=True):
    pp = EB.pairs(d, q)
    m = P2.MODELS[d]
    if q in ("SWC20", "SWC55", "d0_20", "d20_40", "MAR_SWC_0_30", "MAR_SWC_30_50") and pp:
        dates = [d2d(p[0]) for p in pp]
        ax.fill_between(dates, [p[1] - (p[2] if not isnan(p[2]) else 0) for p in pp],
                        [p[1] + (p[2] if not isnan(p[2]) else 0) for p in pp], color=C_OBS, alpha=0.15, lw=0)
        ax.plot(dates, [p[1] for p in pp], color=C_OBS, lw=0.8, label="observed (band: ± sd)")
        ax.plot(dates, [p[3] for p in pp], color=C_MOD, lw=0.9, label="model")
    datefmt(ax)
    ax.set_ylim(0, 0.42)                                  # same range in every soil-water panel
    if legend and d == "WDF25_Ruta_aCO2_WW" and q == "SWC20":
        ax.legend(loc="upper right", fontsize=5.5, borderpad=0.3)


def _scatter_sites(ax, q, decks=None):
    X, Y, XE = [], [], []
    for site in SITE_C:
        xs, ys, es = [], [], []
        for d in (decks or FACE):
            if DECKS[d]["site"] != site:
                continue
            o, s = obs(d, q)
            if q == "crop_N_harvest":
                o, s = _cropN_obs(d)
                mv = P2.CROPN[d][0] if isinstance(P2.CROPN.get(d), (tuple, list)) else P2.CROPN.get(d, NAN)
            else:
                mv = mod(d, q)
            if isnan(o) or isnan(mv):
                continue
            pres = (q == "TKW" and str(P2.CAPS[d][0]).strip() == "KERNEL_MASS_CAP")
            xs.append(o); ys.append(mv); es.append(s)
            if pres:
                ax.plot(o, mv, "o", ms=2.8, mfc=GRAY, mec=GRAY, zorder=4)
        if xs:
            ax.errorbar(xs, ys, xerr=[e if not isnan(e) else 0 for e in es], fmt=SITE_M[site], ms=2.8, mfc="none",
                        mec=SITE_C[site], ecolor=SITE_C[site], elinewidth=0.4, capsize=0, mew=0.6)
            X += xs; Y += ys; XE += es
    onetoone(ax, X, Y)
    place_stats(ax, X, Y, XE=XE)


def _cropN_obs(d):
    o = P2.OBS_REG.get(d, "crop_N_harvest")
    if o:
        return o["value"], o["sd"]
    return NAN, NAN


def _wdf_height_obs(d):
    site = DECKS[d]["site"]
    if site == "WDF25":
        cv, co2, wat = d.split("_")[1:4]
        rows = P2.W25.obs_for(cv, co2, wat)
    else:
        cv, co2 = d.split("_")[1:3]
        rows = P2.W23.obs_for(cv, co2)
    out = [(r["date"], r["value"], r["std"] if r["std"] is not None else NAN) for r in rows
           if r["variable"].startswith("Height") and r["value"] is not None]
    return sorted(out)


def _bra_series(d, var):
    rr = P2.rows_of(P2.bra_obs_path(d))
    out = []
    for x in rr:
        if x["variable"] == var:
            out.append((x["date"], P2.fl(x["mean"]), P2.fl(x["sd"])))
    return sorted(out)


RESP_AGRO = [("An", "Leaf $A_n$"), ("gs", "Leaf $g_{sw}$"), ("Tr", "Leaf $T_r$"),
             ("AGB_harvest", "Above-ground DM"), ("spike", "Spike DM"), ("W_stem", "Stem DM"), ("W_leaf_tot", "Leaf DM"),
             ("grain_yield", "Grain yield"), ("kernel_number", "Kernel number"), ("HI", "Harvest index"),
             ("LAI_peak", "Peak LAI"), ("height_max", "Height")]


RESP_QUAL = [("grain_protein_pct", "Grain protein"), ("grain_N_pct", "Grain N"), ("leaf_N_pct", "Leaf N"), ("stem_N_pct", "Stem N"),
             ("grain_Fe", "Grain Fe"), ("grain_Zn", "Grain Zn"), ("grain_Mn", "Grain Mn"), ("grain_Cu", "Grain Cu"),
             ("grain_Mg", "Grain Mg"), ("grain_Ca", "Grain Ca")]


def _ratio_rows(num, den, q):
    """(R_obs, R_obs_sd, R_model) for a deck pair; gas exchange from the EB window pairs."""
    if q in ("An", "gs", "Tr"):
        pn = {p[0]: p for p in EB.pairs(num, q)}; pd = {p[0]: p for p in EB.pairs(den, q)}
        common = sorted(set(pn) & set(pd))
        if not common:
            return NAN, NAN, NAN
        # ratio of the season means over the common campaign dates (a mean of per-date ratios is dominated by
        # the near-zero winter conductances); SE from the per-date ratios' spread, scaled to the mean ratio
        oa = np.mean([pd[x][1] for x in common]); oe = np.mean([pn[x][1] for x in common])
        ma = np.mean([pd[x][3] for x in common]); me = np.mean([pn[x][3] for x in common])
        ro = oe / oa if oa else NAN; rm = me / ma if ma else NAN
        rr = [pn[x][1] / pd[x][1] for x in common if pd[x][1]]
        se = np.std(rr, ddof=1) / math.sqrt(len(rr)) * (ro / np.mean(rr)) if len(rr) > 1 and np.mean(rr) else NAN
        return float(ro), float(se), float(rm)
    for r in P2.CROWS:
        if r["num"] == num and r["den"] == den and r["quantity"] == q:
            return r["R_obs"], r["R_obs_sd"], r["R_model"]
    # fall back to raw values
    on, osd = obs(num, q); od, dsd = obs(den, q); mn, md = mod(num, q), mod(den, q)
    ro = on / od if not isnan(on) and not isnan(od) and od else NAN
    rm = mn / md if not isnan(mn) and not isnan(md) and md else NAN
    sd = abs(ro) * math.sqrt((osd / on) ** 2 + (dsd / od) ** 2) if not isnan(ro) and not isnan(osd) and not isnan(dsd) and on and od else NAN
    return ro, sd, rm


def _wdf_pairs():
    return OrderedDict([
        ("2024/25", [("WDF25_%s_eCO2_%s" % (cv, w), "WDF25_%s_aCO2_%s" % (cv, w), "%s %s" % (cv, w)) for cv in ("Ruta", "Zahir") for w in ("WW", "WL")]),
        ("2022/23", [("WDF23_%s_eCO2" % cv, "WDF23_%s_aCO2" % cv, cv) for cv in ("Ruta", "Zahir")])])


MAR_DAILY = os.path.join(config.MARICOPA_OBS_DIR, "maricopa_archive_daily.csv")


LS_CO2 = {"aCO2": C_A, "eCO2": C_E}


def _mar_daily(trno, var):
    out = []
    for r in csv.DictReader(open(MAR_DAILY)):
        if int(r["TRNO"]) == trno and r[var + "_mean"]:
            out.append((r["DATE"], float(r[var + "_mean"]), float(r[var + "_std"]) if r[var + "_std"] else NAN))
    return sorted(out)


def _model_daily(d, col, scale=1.0):
    m = P2.MODELS[d]
    return [(m.day2date[int(P2.fl(r["day"]))], P2.fl(r[col]) * scale) for r in m.g if int(P2.fl(r["day"])) in m.day2date]


def _course_panel(ax, decks, obs_fn, col, scale=1.0, ls_by_water=True):
    ymax = 0.0
    for d in decks:
        c = LS_CO2[DECKS[d]["co2"]]
        ls = "--" if (ls_by_water and DECKS[d]["water"] == "WL") else "-"
        ser = _model_daily(d, col, scale)
        ax.plot([d2d(x) for x, _ in ser], [v for _, v in ser], color=c, lw=0.9, ls=ls)
        ob = obs_fn(d)
        if ob:
            mk = "o" if DECKS[d]["water"] == "WW" else "^"
            ax.errorbar([d2d(x) for x, _, _ in ob], [v for _, v, _ in ob], yerr=[(s if not isnan(s) else 0) for _, _, s in ob],
                        fmt=mk, ms=2.4, color=c, mfc="white", mew=0.7, elinewidth=0.5, capsize=0)
            ymax = max([ymax] + [v + (s if not isnan(s) else 0) for _, v, s in ob if not isnan(v)])
        ymax = max([ymax] + [v for _, v in ser if not isnan(v)])
    datefmt(ax); ax.set_ylim(bottom=0)
    return ymax


def _lai_obs(d):
    site = DECKS[d]["site"]
    if site in ("WDF25", "WDF23"):
        lt = P2.OBS_REG.lai_ts.get(d)
        return lt[0] if lt else []
    if site == "BRAUNSCHWEIG":
        return _bra_series(d, "LGrAI")
    return _mar_daily(int(d.split("_")[1]), "LAID")


def _agb_obs(d):
    site = DECKS[d]["site"]
    if site == "BRAUNSCHWEIG":
        return _bra_series(d, "SUMDW")
    if site == "MARICOPA":
        return [(x, v / 10.0, s / 10.0 if not isnan(s) else NAN) for x, v, s in _mar_daily(int(d.split("_")[1]), "CWAD")]
    # WDF: harvest value only
    o, s = obs(d, "AGB_harvest")
    if isnan(o):
        return []
    m = P2.MODELS[d]
    last = m.day2date[max(m.day2date)]
    return [(last, o, s)]


PANELS = [  # (title, decks) — every site × cultivar/N level/season; aCO2 and eCO2 in each panel, WL dashed
    ("WDF 2022/23 Ruta", ["WDF23_Ruta_aCO2", "WDF23_Ruta_eCO2"]),
    ("WDF 2022/23 Zahir", ["WDF23_Zahir_aCO2", "WDF23_Zahir_eCO2"]),
    ("WDF 2024/25 Ruta", ["WDF25_Ruta_aCO2_WW", "WDF25_Ruta_eCO2_WW", "WDF25_Ruta_aCO2_WL", "WDF25_Ruta_eCO2_WL"]),
    ("WDF 2024/25 Zahir", ["WDF25_Zahir_aCO2_WW", "WDF25_Zahir_eCO2_WW", "WDF25_Zahir_aCO2_WL", "WDF25_Zahir_eCO2_WL"]),
    ("Braunschweig 2013/14\n40 kg N ha$^{-1}$", ["BRAUNSCHWEIG_101", "BRAUNSCHWEIG_104"]),
    ("Braunschweig 2013/14\n180 kg N ha$^{-1}$", ["BRAUNSCHWEIG_102", "BRAUNSCHWEIG_105"]),
    ("Braunschweig 2013/14\n320 kg N ha$^{-1}$", ["BRAUNSCHWEIG_103", "BRAUNSCHWEIG_106"]),
    ("Braunschweig 2013/14\n180 kg N ha$^{-1}$ (warming control)", ["BRAUNSCHWEIG_107", "BRAUNSCHWEIG_110"]),
    ("Braunschweig 2014/15\n35 kg N ha$^{-1}$", ["BRAUNSCHWEIG_113", "BRAUNSCHWEIG_116"]),
    ("Braunschweig 2014/15\n200 kg N ha$^{-1}$", ["BRAUNSCHWEIG_114", "BRAUNSCHWEIG_117"]),
    ("Braunschweig 2014/15\n320 kg N ha$^{-1}$", ["BRAUNSCHWEIG_115", "BRAUNSCHWEIG_118"]),
    ("Braunschweig 2014/15\n202 kg N ha$^{-1}$ (warming control)", ["BRAUNSCHWEIG_119", "BRAUNSCHWEIG_122"]),
    ("Maricopa 1993/94\nwell-watered, water-limited", ["MARICOPA_906", "MARICOPA_908", "MARICOPA_905", "MARICOPA_907"]),
    ("Maricopa 1995/96 well-watered", ["MARICOPA_910", "MARICOPA_912"]),
    ("Maricopa 1996/97 high N", ["MARICOPA_914", "MARICOPA_916"]),
    ("Maricopa 1996/97 low N", ["MARICOPA_913", "MARICOPA_915"]),
]


def _course_legend(fig, y=-0.01):
    h = [Line2D([], [], color=C_A, label="model, ambient CO$_2$"), Line2D([], [], color=C_E, label="model, elevated CO$_2$"),
         Line2D([], [], color=C_A, ls="--", label="model, water-limited (dashed, CO$_2$ color)"),
         Line2D([], [], marker="o", ls="", color=C_OBS, mfc="white", label="observed ± sd (well-watered)"),
         Line2D([], [], marker="^", ls="", color=C_OBS, mfc="white", label="observed ± sd (water-limited)")]
    fig.legend(handles=h, loc="lower center", ncol=5, bbox_to_anchor=(0.5, y), fontsize=5.8, handletextpad=0.4, columnspacing=1.2, frameon=False)


def figN_canopy(out):
    fig, axes = plt.subplots(4, 4, figsize=(W2, 1.02 * W2), sharey=True); ymax = 0.0
    for k, (ax, (ttl, decks)) in enumerate(zip(axes.flat, PANELS)):
        ymax = max(ymax, _course_panel(ax, decks, _lai_obs, "LAI_m2_m2"))
        ax.set_title(ttl, pad=3, fontsize=7); letter(ax, k)
        if k % 4 == 0: ax.set_ylabel("Green LAI (m$^2$ m$^{-2}$)")
    axes[0, 0].set_ylim(0, 1.06 * ymax)
    _course_legend(fig)
    fig.tight_layout(h_pad=1.2, w_pad=1.0, rect=(0, 0.03, 1, 1))
    save(fig, out, "FigN04_leaf_area_all_sites")


def figN_biomass(out):
    fig, axes = plt.subplots(4, 4, figsize=(W2, 1.02 * W2), sharey=True); ymax = 0.0
    for k, (ax, (ttl, decks)) in enumerate(zip(axes.flat, PANELS)):
        ymax = max(ymax, _course_panel(ax, decks, _agb_obs, "AGB_g_m2"))
        ax.set_title(ttl, pad=3, fontsize=7); letter(ax, k)
        if k % 4 == 0: ax.set_ylabel("Above-ground DM (g m$^{-2}$)")
    axes[0, 0].set_ylim(0, 1.06 * ymax)
    _course_legend(fig)
    fig.tight_layout(h_pad=1.2, w_pad=1.0, rect=(0, 0.03, 1, 1))
    save(fig, out, "FigN05_biomass_all_sites")


def figN_harvest(out):
    qq = [("LAI_peak", "Peak LAI (m$^2$ m$^{-2}$)"), ("AGB_harvest", "Above-ground DM (g m$^{-2}$)"),
          ("grain_yield", "Grain yield (g m$^{-2}$)"), ("kernel_number", "Kernels (m$^{-2}$)"),
          ("TKW", "Kernel weight, TKW (g)"), ("HI", "Harvest index (–)"),
          ("grain_N_pct", "Grain N (% DM)"), ("grain_protein_pct", "Grain protein (% DM)"),
          ("height_max", "Height (m)"), ("crop_N_harvest", "Above-ground N (g N m$^{-2}$)"),
          ("grain_Fe", "Grain Fe (mg kg$^{-1}$)"), ("grain_Zn", "Grain Zn (mg kg$^{-1}$)")]
    fig, axes = plt.subplots(3, 4, figsize=(W2, 0.8 * W2))
    for k, (ax, (q, lab)) in enumerate(zip(axes.flat, qq)):
        _scatter_sites(ax, q)
        ax.set_title(lab, pad=3); ax.set_xlabel("Observed"); ax.set_ylabel("Simulated"); letter(ax, k, x=-0.25)
    h = [Line2D([], [], marker=SITE_M[s], ls="", mfc="none", mec=SITE_C[s], label=LEG_N[s]) for s in LEG_SITES]
    h.append(Line2D([], [], marker="o", ls="", mfc=GRAY, mec=GRAY, label="(e) kernel mass at the prescribed potential (WDF 2024/25)"))
    fig.legend(handles=h, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.01), fontsize=5.8, frameon=False)
    fig.tight_layout(h_pad=1.2, w_pad=1.0, rect=(0, 0.03, 1, 1))
    save(fig, out, "FigN06_harvest_all_sites")


BRA_ECO2 = [("BRAUNSCHWEIG_104", "BRAUNSCHWEIG_101", "N1 2014"), ("BRAUNSCHWEIG_105", "BRAUNSCHWEIG_102", "N2 2014"),
            ("BRAUNSCHWEIG_106", "BRAUNSCHWEIG_103", "N3 2014"), ("BRAUNSCHWEIG_110", "BRAUNSCHWEIG_107", "N4 2014"),
            ("BRAUNSCHWEIG_116", "BRAUNSCHWEIG_113", "N1 2015"), ("BRAUNSCHWEIG_117", "BRAUNSCHWEIG_114", "N2 2015"),
            ("BRAUNSCHWEIG_118", "BRAUNSCHWEIG_115", "N3 2015"), ("BRAUNSCHWEIG_122", "BRAUNSCHWEIG_119", "N4 2015")]


MAR_ECO2 = [("MARICOPA_907", "MARICOPA_905", "Dry 93/94"), ("MARICOPA_908", "MARICOPA_906", "Wet 93/94"),
            ("MARICOPA_912", "MARICOPA_910", "Wet 95/96"), ("MARICOPA_916", "MARICOPA_914", "High-N 96/97"), ("MARICOPA_915", "MARICOPA_913", "Low-N 96/97")]


AGRO_ALL = [("An", "Leaf $A_n$"), ("gs", "Leaf $g_{sw}$"), ("Tr", "Leaf $T_r$"), ("LAI_peak", "Peak LAI"), ("height_max", "Height"),
            ("AGB_harvest", "Above-ground DM"), ("W_stem", "Stem DM"), ("W_leaf_tot", "Leaf DM"), ("spike", "Spike DM"),
            ("kernel_number", "Kernel number"), ("TKW", "Kernel mass"), ("grain_yield", "Grain yield"), ("HI", "Harvest index")]


QUAL_ALL = [("crop_N_harvest", "Above-ground N"), ("grain_N_pct", "Grain N"), ("grain_protein_pct", "Grain protein"), ("leaf_N_pct", "Leaf N"), ("stem_N_pct", "Stem N"),
            ("grain_Fe", "Grain Fe"), ("grain_Zn", "Grain Zn"), ("grain_Mn", "Grain Mn"), ("grain_Cu", "Grain Cu"), ("grain_Mg", "Grain Mg"), ("grain_Ca", "Grain Ca")]


def _consolidated_csv():
    return config.CONSOLIDATED_OBS


def _ratio_rows2(num, den, q):
    if q == "crop_N_harvest":
        on, _ = _cropN_obs(num); od, _ = _cropN_obs(den)
        mn = P2.CROPN[num][0] if isinstance(P2.CROPN.get(num), (tuple, list)) else P2.CROPN.get(num, NAN)
        md = P2.CROPN[den][0] if isinstance(P2.CROPN.get(den), (tuple, list)) else P2.CROPN.get(den, NAN)
        ro = on / od if not isnan(on) and not isnan(od) and od else NAN
        rm = mn / md if not isnan(mn) and not isnan(md) and md else NAN
        return ro, NAN, rm
    if DECKS.get(num, {}).get("site", num.split("_")[0]) == "BRAUNSCHWEIG" and q in ("grain_Fe", "grain_Zn"):
        # observations from the consolidated observation table (grain quality sheet), model from the mineral module;
        # a treatment without a row in that table (no grain mineral observation) contributes no ratio
        def g(dk):
            o_, s_ = _cv_obs(dk, q)
            if isnan(o_):
                return NAN, NAN, NAN
            return o_, s_, P2.GRAIN_MIN[dk].get(q[6:], NAN)
        on, osd, mn = g(num); od, dsd, md = g(den)
        ro = on / od if not isnan(on) and od else NAN; rm = mn / md if not isnan(mn) and md else NAN
        sd = abs(ro) * math.sqrt((osd / on) ** 2 + (dsd / od) ** 2) if not isnan(ro) and not isnan(osd) and not isnan(dsd) and on and od else NAN
        return ro, sd, rm
    if DECKS[num]["site"] == "MARICOPA" and q in ("grain_yield", "HI", "LAI_peak"):
        on, osd = obs(num, q); od, dsd = obs(den, q); mn, md = mod(num, q), mod(den, q)
        ro = on / od if not isnan(on) and od else NAN; rm = mn / md if not isnan(mn) and md else NAN
        sd = abs(ro) * math.sqrt((osd / on) ** 2 + (dsd / od) ** 2) if not isnan(ro) and not isnan(osd) and not isnan(dsd) and on and od else NAN
        return ro, sd, rm
    return _ratio_rows(num, den, q)


def _forest2(ax, items, pairs, title, xlabel="Effect of elevated CO$_2$ (%)", show_pairs=True):
    """Forest panel: observed mean effect over the treatment pairs with its 95 % interval, simulated mean and range."""
    rows = []
    for q, lab in items:
        RO, RM, SE = [], [], []
        for num, den, tag in pairs:
            ro, se, rm = _ratio_rows2(num, den, q)
            if not isnan(ro): RO.append(ro); SE.append(se)
            if not isnan(rm): RM.append(rm)
        pe = None
        rows.append((lab, RO, SE, RM, pe))          # same rows in every panel; empty where not measured
    n = len(rows); y = np.arange(n)[::-1]
    for i, (lab, RO, SE, RM, pe) in enumerate(rows):
        if RM:
            lo, hi = (min(RM) - 1) * 100, (max(RM) - 1) * 100
            ax.plot([lo, hi], [y[i], y[i]], color=C_SIM, lw=3.0, alpha=0.35, solid_capstyle="butt", zorder=2)
            ax.plot((np.mean(RM) - 1) * 100, y[i], "s", ms=3.4, color=C_SIM, zorder=4)
        if pe is not None:
            eff, lo, hi, sig = pe
            if lo is not None:
                ax.errorbar(eff, y[i], xerr=[[eff - lo], [hi - eff]], fmt="none", ecolor=C_OBS, elinewidth=0.7, capsize=0, zorder=5)
            ax.plot(eff, y[i], "o", ms=3.6, color=C_OBS, mfc=C_OBS if sig else "white", mew=0.8, zorder=6)
        elif RO:
            m = np.mean(RO); k = len(RO)
            within = sum((x if not isnan(x) else 0.0) ** 2 for x in SE) / k ** 2   # replicate dispersion inside the treatment pairs
            between = (np.var(RO, ddof=1) / k) if k > 1 else 0.0                      # spread between treatment pairs
            e = 1.96 * math.sqrt(within + between)                                    # 95 % interval of the mean effect
            ax.errorbar((m - 1) * 100, y[i], xerr=e * 100, fmt="o", ms=3.4, color=C_OBS, mfc=C_OBS, mew=0.8, elinewidth=0.7, capsize=0, zorder=5)
            if show_pairs and len(RO) > 1:
                ax.plot([(v - 1) * 100 for v in RO], [y[i] + 0.28] * len(RO), "|", ms=3, color=C_OBS, alpha=0.5, zorder=3)
    ax.axvline(0, color=GRAY, lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.set_ylim(-0.6, n - 0.4)
    for i, (lab, RO, SE, RM, pe) in enumerate(rows):
        if not RO and not RM:
            ax.text(0.985, y[i], "not measured", fontsize=5.2, color=GRAY, ha="right", va="center", style="italic",
                    transform=ax.get_yaxis_transform())
        elif not RO:
            ax.text(0.985, y[i], "no observation", fontsize=5.2, color=GRAY, ha="right", va="center", style="italic",
                    transform=ax.get_yaxis_transform(), zorder=7, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    if any((not RO) for _, RO, _, _, _ in rows):          # reserve a right margin for the "no obs." / "not measured" labels
        lo, hi = ax.get_xlim(); ax.set_xlim(lo, hi + 0.24 * (hi - lo))
    ax.set_xlabel(xlabel); ax.set_title(title, pad=4, fontsize=7.5, loc="left")
    ax.grid(axis="x", color="#e6e6e6", lw=0.5); ax.set_axisbelow(True)


def _same_x(axes):
    lo = min(a.get_xlim()[0] for a in axes); hi = max(a.get_xlim()[1] for a in axes)
    for a in axes: a.set_xlim(lo, hi)


def _forest_legend(fig, y=-0.03, pub=True):
    h = [Line2D([], [], marker="o", ls="", color=C_OBS, label="observed: mean over treatment pairs ± 95 % CI (ticks: single pairs)"),
         Line2D([], [], marker="s", ls="-", lw=3, alpha=0.5, color=C_SIM, label="simulated: mean over treatment pairs (bar: range)")]
    fig.legend(handles=h, loc="lower center", ncol=2, bbox_to_anchor=(0.5, y), fontsize=5.8, frameon=False)


def figN_eco2_agro(out):
    pb = _wdf_pairs()
    fig, axes = plt.subplots(2, 2, figsize=(W2, 0.95 * W2)); axes = list(axes.flat)
    _forest2(axes[0], AGRO_ALL, pb["2022/23"], "WheatDryFACE 2022/23, low N (2 genotypes)")
    _forest2(axes[1], AGRO_ALL, pb["2024/25"], "WheatDryFACE 2024/25, high N (2 genotypes × 2 water levels)")
    _forest2(axes[2], AGRO_ALL, BRA_ECO2, "Braunschweig 2013/14 + 2014/15 (4 N levels × 2 seasons)")
    _forest2(axes[3], AGRO_ALL, MAR_ECO2, "Maricopa 1993–97 (5 water/N treatments × seasons)")
    _same_x(axes)
    for k, ax in enumerate(axes): letter(ax, k, x=-0.42, y=1.02)
    _forest_legend(fig, y=0.0)
    fig.tight_layout(w_pad=4.0, h_pad=2.0, rect=(0, 0.02, 1, 1))
    save(fig, out, "FigN07_eCO2_response_agronomy_all_sites")


def figN_eco2_quality(out):
    pb = _wdf_pairs()
    fig, axes = plt.subplots(2, 2, figsize=(W2, 0.85 * W2)); axes = list(axes.flat)
    _forest2(axes[0], QUAL_ALL, pb["2022/23"], "WheatDryFACE 2022/23, low N (2 genotypes)")
    _forest2(axes[1], QUAL_ALL, pb["2024/25"], "WheatDryFACE 2024/25, high N (2 genotypes × 2 water levels)")
    _forest2(axes[2], QUAL_ALL, BRA_ECO2, "Braunschweig 2013/14 + 2014/15 (4 N levels × 2 seasons)")
    _forest2(axes[3], QUAL_ALL, MAR_ECO2, "Maricopa 1993–97 (5 water/N treatments × seasons)")
    _same_x(axes)
    for k, ax in enumerate(axes): letter(ax, k, x=-0.42, y=1.02)
    fig.text(0.5, 0.035, "Maricopa: grain N and protein from the harvest samples; grain minerals and organ N were not measured there.",
             ha="center", fontsize=5.8, color="#404040")
    _forest_legend(fig, y=0.0)
    fig.tight_layout(w_pad=4.0, h_pad=2.0, rect=(0, 0.02, 1, 1))
    save(fig, out, "FigN08_eCO2_response_quality_all_sites")


def figN_water_N(out):
    wl = [("WDF25_%s_%s_WL" % (cv, co2), "WDF25_%s_%s_WW" % (cv, co2), "%s %s" % (cv, co2)) for cv in ("Ruta", "Zahir") for co2 in ("aCO2", "eCO2")]
    mar_w = [("MARICOPA_905", "MARICOPA_906", "aCO2"), ("MARICOPA_907", "MARICOPA_908", "eCO2")]
    bra_n = [("BRAUNSCHWEIG_101", "BRAUNSCHWEIG_103", "N1/N3 aCO2"), ("BRAUNSCHWEIG_104", "BRAUNSCHWEIG_106", "N1/N3 eCO2")]
    bra_n2 = [("BRAUNSCHWEIG_102", "BRAUNSCHWEIG_103", "N2/N3 aCO2"), ("BRAUNSCHWEIG_105", "BRAUNSCHWEIG_106", "N2/N3 eCO2")]
    mar_n = [("MARICOPA_913", "MARICOPA_914", "aCO2"), ("MARICOPA_915", "MARICOPA_916", "eCO2")]
    items = AGRO_ALL + QUAL_ALL
    fig, axes = plt.subplots(2, 2, figsize=(W2, 1.45 * W2)); axes = list(axes.flat)
    bra_n = bra_n + [("BRAUNSCHWEIG_113", "BRAUNSCHWEIG_115", "N1/N3 aCO2 2015"), ("BRAUNSCHWEIG_116", "BRAUNSCHWEIG_118", "N1/N3 eCO2 2015")]
    _forest2(axes[0], items, wl, "WheatDryFACE 2024/25: WL / WW", xlabel="Effect of water limitation (%)")
    _forest2(axes[1], items, bra_n, "Braunschweig: 40 (35) / 320 kg N ha$^{-1}$", xlabel="Effect of low N supply (%)")
    _forest2(axes[2], items, mar_w, "Maricopa 1993/94: WL / WW", xlabel="Effect of water limitation (%)")
    _forest2(axes[3], items, mar_n, "Maricopa 1996/97: low / high N", xlabel="Effect of low N supply (%)")
    _same_x(axes[:2]); _same_x(axes[2:])
    for k, ax in enumerate(axes): letter(ax, k, x=-0.42, y=1.02)
    _forest_legend(fig, y=0.0)
    fig.tight_layout(w_pad=4.0, h_pad=2.0, rect=(0, 0.02, 1, 1))
    save(fig, out, "FigN09_water_and_N_responses_all_sites")


MINERALS = {"grain_Fe", "grain_Zn", "grain_Mn", "grain_Cu", "grain_Mg", "grain_Ca"}


def figN_signs_site2(out):
    from collections import Counter
    kinds = [("eCO2", "Elevated CO$_2$"), ("water", "Water limitation"), ("N", "Nitrogen supply"), ("cultivar", "Genotype")]
    sites = ["WDF23", "WDF25", "BRAUNSCHWEIG", "MARICOPA"]
    cnt = {}
    for r in P2.CROWS:
        if r["cls"] == "C": continue
        grp = "min" if r["quantity"] in MINERALS else "agr"
        cnt.setdefault((r["kind"], r["site"], grp), Counter())[r["cls"]] += 1
    fig, axes = plt.subplots(1, 4, figsize=(W2, 0.40 * W2), sharey=True)
    cols = {"agree": "#4d4d4d", "B": "#b0b0b0", "N": "#f0a8ab", "A": C_SIM}
    off = {"agr": 0.19, "min": -0.19}
    for k, (ax, (kind, ttl)) in enumerate(zip(axes, kinds)):
        ys = np.arange(len(sites))[::-1]
        for yi, site in zip(ys, sites):
            any_ = False
            for grp in ("agr", "min"):
                c = cnt.get((kind, site, grp), Counter()); tot = c["agree"] + c["A"] + c["B"] + c["N"]
                if tot == 0: continue
                any_ = True; left = 0; y = yi + off[grp]
                for cls in ("agree", "B", "N", "A"):
                    w = c[cls] / tot
                    ax.barh(y, w, left=left, color=cols[cls], height=0.34, edgecolor="white", lw=0.4, hatch=None if grp == "agr" else "////", alpha=1.0)
                    if c[cls]: ax.text(left + w / 2, y, str(c[cls]), ha="center", va="center", fontsize=5.4, color="white" if cls != "B" else "#1a1a1a")
                    left += w
                ax.text(1.03, y, "%d" % tot, va="center", fontsize=5.2, color=GRAY)
            if not any_:
                ax.text(0.5, yi, "not tested", fontsize=5.6, color=GRAY, style="italic", va="center", ha="center")
        ax.set_yticks(ys); ax.set_yticklabels([SITE_N[s].replace("WheatDryFACE", "WDF") for s in sites])
        ax.set_xlim(0, 1.3); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "50", "100"]); ax.set_xlabel("Share of contrasts (%)")
        ax.set_title(ttl, pad=3, loc="right"); letter(ax, k, x=0.0, y=1.02, inside=False)
        ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    h = [Patch(color=cols["agree"], label="same sign"), Patch(color=cols["B"], label="opposite sign, observed effect not significant"), Patch(color=cols["N"], label="observed effect significant, simulated < 5 %"), Patch(color=cols["A"], label="opposite sign, observed effect significant")]
    fig.legend(handles=h, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.07), fontsize=5.8, frameon=False)
    fig.text(0.5, -0.13, "Solid bars: the twelve harvest and canopy variables listed in Sect. 5.  Hatched bars: the six grain minerals.  Numbers at right: contrasts with an observation.", ha="center", fontsize=6.0, color="#404040")
    fig.tight_layout(w_pad=1.5)
    save(fig, out, "FigN11_sign_agreement_site_two_bars")


def figN_genotype(out):
    """Supplement S2, Fig. S6: Zahir over Ruta at WheatDryFACE, both seasons."""
    g23 = [("WDF23_Zahir_%s" % c, "WDF23_Ruta_%s" % c, c) for c in ("aCO2", "eCO2")]
    g25 = [("WDF25_Zahir_%s_%s" % (c, w), "WDF25_Ruta_%s_%s" % (c, w), c + w) for c in ("aCO2", "eCO2") for w in ("WW", "WL")]
    items = AGRO_ALL + QUAL_ALL
    fig, axes = plt.subplots(1, 2, figsize=(W2, 0.8 * W2)); axes = list(axes.flat)
    _forest2(axes[0], items, g23, "WheatDryFACE 2022/23: Zahir / Ruta (2 CO$_2$ levels)", xlabel="Genotype difference (%)")
    _forest2(axes[1], items, g25, "WheatDryFACE 2024/25: Zahir / Ruta (2 CO$_2$ × 2 water)", xlabel="Genotype difference (%)")
    _same_x(axes)
    for k, ax in enumerate(axes): letter(ax, k, x=-0.42, y=1.02)
    _forest_legend(fig, y=0.0)
    fig.tight_layout(w_pad=4.0, rect=(0, 0.02, 1, 1))
    save(fig, out, "FigS6_genotype_contrast_WDF")


