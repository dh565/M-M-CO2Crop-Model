#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Part 2, Tables 4 and 5 and the numbers quoted in the text, from the reference-outputs archive.

Writes to OUT_DIR/part2_tables/:
  part2_numbers.json                        every number of this script (scalars by deck, agreement statistics by
                                            site and variable, time-series statistics, treatment contrasts,
                                            sign classes, canopy-temperature and phenology summaries)
  table4_agreement_by_site_variable.csv     Table 4: n, observed and simulated means, bias, RMSD, fraction within
                                            1 sd (or 15 %) by site and variable (harvest and canopy scalars)
  table4_series_statistics.csv              Table 4, time-series rows: n, bias, RMSD, r by site and quantity
  table5_numerical_checks.csv               Table 5: coupling-loop passes, convergence and conservation residuals
  table5_numerical_checks.json              the same, unformatted
"""
import csv
import json
import math
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
import config
import paperfigs as PF

NAN = float("nan")
isnan = PF.isnan


# ============================================================================== Table 4 and text numbers
def st(X, Y):
    p = [(a, b) for a, b in zip(X, Y) if not isnan(a) and not isnan(b)]
    if len(p) < 2:
        return dict(n=len(p))
    o = np.array([a for a, b in p]); m = np.array([b for a, b in p])
    return dict(n=len(p), bias=float((m - o).mean()), rmse=float(math.sqrt(((m - o) ** 2).mean())),
                r=float(np.corrcoef(o, m)[0, 1]) if o.std() > 0 and m.std() > 0 else NAN, obs_mean=float(o.mean()))


Q = ["LAI_peak", "AGB_harvest", "W_leaf_tot", "W_stem", "spike", "chaff", "height_max", "grain_yield", "kernel_number", "TKW", "HI",
     "grain_protein_pct", "grain_N_pct", "leaf_N_pct", "stem_N_pct", "crop_N_harvest", "ears_m2", "grain_Fe", "grain_Zn", "grain_Cu", "grain_Mg", "grain_Mn", "grain_Ca"]


def numbers():
    P2, EB, DECKS = PF.P2, PF.EB, PF.DECKS
    out = {}
    # ---- 1. all decks: scalar table (obs, obs_sd, model)
    scal = {}
    for d in PF.FACE:
        scal[d] = {}
        for q in Q:
            if q == "crop_N_harvest":
                o, s = PF._cropN_obs(d); mv = P2.CROPN[d][0] if isinstance(P2.CROPN.get(d), (tuple, list)) else P2.CROPN.get(d, NAN)
            elif q == "ears_m2":
                o, s, mv = P2.MAR_EARS.get(d, NAN), NAN, P2.EARS.get(d, NAN)
                if isinstance(mv, (tuple, list)): mv = mv[0]
            else:
                o, s = PF.obs(d, q); mv = PF.mod(d, q)
                if q == "grain_yield" and isnan(o) and DECKS[d]["site"] == "MARICOPA":
                    o, s = PF.obs(d, "grain_yield_HWAH")
            scal[d][q] = (o, s, mv)
    out["scalars"] = {d: {q: [None if isnan(x) else round(x, 4) for x in v] for q, v in dd.items()} for d, dd in scal.items()}

    # ---- 2. Table 4: by variable and site: n, bias, RMSE (%), fraction within 1 sd (or 15 % where no sd)
    t3 = []
    for site in ["WDF25", "WDF23", "BRAUNSCHWEIG", "MARICOPA"]:
        decks = [d for d in PF.FACE if DECKS[d]["site"] == site]
        for q in Q:
            O, M, W = [], [], []
            for d in decks:
                o, s, mv = scal[d][q]
                if isnan(o) or isnan(mv): continue
                O.append(o); M.append(mv)
                tol = s if (not isnan(s) and s > 0) else 0.15 * abs(o)
                W.append(abs(mv - o) <= tol)
            if O:
                s_ = st(O, M)
                t3.append(dict(site=site, quantity=q, n=len(O), bias=s_.get("bias", NAN), rmse=s_.get("rmse", NAN),
                               nrmse=100 * s_.get("rmse", NAN) / s_["obs_mean"] if s_.get("obs_mean") else NAN,
                               frac=sum(W) / len(W), obs_mean=float(np.mean(O)), mod_mean=float(np.mean(M)), ratio=float(np.mean(M) / np.mean(O))))
    out["table3"] = t3

    # ---- 3. time-series families via the EB pairs
    fam = {}
    for site in ["WDF25", "WDF23", "BRAUNSCHWEIG", "MARICOPA"]:
        decks = [d for d in PF.FACE if DECKS[d]["site"] == site]
        for q in EB.Q:
            X, Y = [], []
            for d in decks:
                pp = EB.pairs(d, q)
                X += [p[1] for p in pp]; Y += [p[3] for p in pp]
            if X:
                fam[(site, q)] = st(X, Y)
    out["series"] = {"%s|%s" % k: v for k, v in fam.items()}

    # ---- 4. contrasts (eCO2) WDF25/WDF23 means per quantity, incl. gas exchange
    def contrast_summary(pairs, items):
        res = {}
        for q, lab in items:
            RO, RM, SE = [], [], []
            for num, den, tag in pairs:
                ro, se, rm = PF._ratio_rows(num, den, q)
                if not isnan(ro): RO.append(ro); SE.append(se)
                if not isnan(rm): RM.append(rm)
            res[q] = dict(obs=float(np.mean(RO)) if RO else NAN, obs_se=float(np.nanmean(SE)) if RO and any(not isnan(x) for x in SE) else NAN,
                          model=float(np.mean(RM)) if RM else NAN, model_min=min(RM) if RM else NAN, model_max=max(RM) if RM else NAN,
                          obs_pairs=[round(x, 3) for x in RO], model_pairs=[round(x, 3) for x in RM])
        return res
    pb = PF._wdf_pairs()
    out["eco2_wdf25"] = contrast_summary(pb["2024/25"], PF.RESP_AGRO + PF.RESP_QUAL)
    out["eco2_wdf23"] = contrast_summary(pb["2022/23"], PF.RESP_AGRO + PF.RESP_QUAL)
    wl = [("WDF25_%s_%s_WL" % (cv, co2), "WDF25_%s_%s_WW" % (cv, co2), "%s %s" % (cv, co2)) for cv in ("Ruta", "Zahir") for co2 in ("aCO2", "eCO2")]
    out["water_wdf25"] = contrast_summary(wl, PF.RESP_AGRO)
    bra = [("N1", "104", "101"), ("N2", "105", "102"), ("N3", "106", "103"), ("N4", "110", "107"), ("N1 2015", "116", "113"), ("N2 2015", "117", "114"), ("N3 2015", "118", "115"), ("N4 2015", "122", "119")]
    out["eco2_bra"] = {nl: contrast_summary([("BRAUNSCHWEIG_%s" % e, "BRAUNSCHWEIG_%s" % a, nl)], PF.RESP_AGRO + PF.RESP_QUAL) for nl, e, a in bra}
    mar = [("Dry 93/94", "907", "905"), ("Wet 93/94", "908", "906"), ("Wet 95/96", "912", "910"), ("Wet 96/97", "916", "914"), ("Low-N 96/97", "915", "913")]
    out["eco2_mar"] = {lab: contrast_summary([("MARICOPA_%s" % e, "MARICOPA_%s" % a, lab)], PF.RESP_AGRO) for lab, e, a in mar}
    out["water_mar"] = {lab: contrast_summary([("MARICOPA_%s" % n, "MARICOPA_%s" % d_, lab)], PF.RESP_AGRO) for lab, n, d_ in (("aCO2 93/94", "905", "906"), ("eCO2 93/94", "907", "908"))}

    # ---- 5. sign classes from the contrast rows by site and kind
    cls = {}
    for r in P2.CROWS:
        key = (r["site"], r["kind"])
        cls.setdefault(key, Counter())[r["cls"]] += 1
    out["sign_classes"] = {"%s|%s" % k: dict(v) for k, v in cls.items()}
    out["classA_rows"] = [dict(contrast=r["contrast"], q=r["quantity"], ro=r["R_obs"], rm=r["R_model"]) for r in P2.CROWS if r["cls"] == "A"]

    # ---- 6. canopy temperature, thermal images
    O = []; M = []; D = {}
    for d in [x for x in PF.FACE if DECKS[x]["site"] == "WDF25"]:
        for p in EB.pairs(d, "Tleaf_FLIR"):
            O.append(p[1]); M.append(p[3]); D.setdefault(p[0], {})[d] = (p[1], p[3])
    out["flir"] = st(O, M)
    dob, dmo = [], []
    for x, dd in D.items():
        for e in dd:
            a = e.replace("eCO2", "aCO2")
            if "eCO2" in e and a in dd:
                dob.append(dd[e][0] - dd[a][0]); dmo.append(dd[e][1] - dd[a][1])
    out["flir_eco2_minus_aco2"] = dict(obs=float(np.mean(dob)), obs_se=float(np.std(dob, ddof=1) / math.sqrt(len(dob))), model=float(np.mean(dmo)), n=len(dob))

    # ---- 7. phenology
    ph = {}
    for d in PF.FACE:
        for st_, md, od, diff, scored, src in P2.PHEN[d]:
            if od and not isnan(diff):
                ph.setdefault(DECKS[d]["site"], []).append(abs(diff))
    out["phenology_abs_dev_days"] = {s: dict(n=len(v), mean=float(np.mean(v)), max=float(max(v)), within2=float(np.mean([x <= 2 for x in v]))) for s, v in ph.items()}
    return out


# ============================================================================== Table 5
SITE_GROUPS = [("WDF25", "WDF25_"), ("WDF23", "WDF23_"), ("BRA", "BRAUNSCHWEIG_"), ("MAR", "MARICOPA_")]
SITE_TITLES = {"WDF25": "WheatDryFACE 2024/25", "WDF23": "WheatDryFACE 2022/23", "BRA": "Braunschweig 2014, 2015", "MAR": "Maricopa 1993-97"}


def _grab(txt, pat):
    m = re.search(pat, txt)
    return float(m.group(1)) if m else None


def table5(ref):
    """Numerical behaviour and conservation checks over the FACE simulations, by site group.

    Per simulation, from daily_eb_actual.csv: outer_iters (coupling-loop passes per day), outer_converged (day converged
    to 0.05 K within five passes), outer_resid_K (remaining change of the daily-mean committed canopy temperature on the
    last pass), canopy_resid_W_m2; from run.log: the seasonal soil-water residual (WB Residual, daily sum), the soil-heat
    storage closure residual per sub-step, the organ-state carbon closure, the nitrogen-payment identity residual, the
    organ-state N closure and the Richards failed-day statement; from daily_nitrogen.csv: N_balance_resid_gN_m2; from
    daily_minerals.csv: balance_resid_mg_m2.  Maxima are of absolute values over all days and simulations of the group.
    """
    decks = sorted(d for d in os.listdir(ref) if os.path.isdir(os.path.join(ref, d, "output")))
    out = {}
    for g, prefix in SITE_GROUPS:
        dk = [d for d in decks if d.startswith(prefix)]
        eb, wb, heat, carbon, npay, norg, nbal, mineral = [], [], [], [], [], [], [], []
        rich = 0
        for d in dk:
            e = pd.read_csv(os.path.join(ref, d, "output", "daily_eb_actual.csv"), skipinitialspace=True)
            e.columns = [c.strip() for c in e.columns]; eb.append(e)
            txt = open(os.path.join(ref, d, "run.log"), errors="ignore").read()
            wb.append(abs(_grab(txt, r"WB Residual \(daily sum, mm\):\s*([-0-9.E+]+)")))
            heat.append(_grab(txt, r"heat-storage closure residual\| =\s*([-0-9.E+]+)"))
            carbon.append(_grab(txt, r"CARBON CLOSURE max\|organ-state - dm_gain\| =\s*([-0-9.E+]+)"))
            npay.append(_grab(txt, r"gate N identity max\|resid\| =\s*([-0-9.E+]+)"))
            norg.append(_grab(txt, r"max \|organ-state N closure\| =\s*([-0-9.E+]+)"))
            rich += 0 if re.search(r"Richards had no failed days", txt) else 1
            n = pd.read_csv(os.path.join(ref, d, "output", "daily_nitrogen.csv"), skipinitialspace=True)
            n.columns = [c.strip() for c in n.columns]; nbal.append(n["N_balance_resid_gN_m2"].abs().max())
            mi = pd.read_csv(os.path.join(ref, d, "output", "daily_minerals.csv"), skipinitialspace=True)
            mi.columns = [c.strip() for c in mi.columns]; mineral.append(mi["balance_resid_mg_m2"].abs().max())
        E = pd.concat(eb); un = E[E.outer_converged < 0.5]
        out[g] = dict(n=len(dk), canopy_resid=float(E.canopy_resid_W_m2.abs().max()), wb=max(wb), heat=max(heat), carbon=max(carbon),
                      npay=max(npay), norg=max(norg), richards=rich, nbal=float(max(nbal)), mineral=float(max(mineral)),
                      passes_mean=float(E.outer_iters.mean()), passes_max=float(E.outer_iters.max()), conv_frac=float(E.outer_converged.mean()),
                      unconv_days=int(len(un)), ndays=int(len(E)), resid_unconv_median=float(un.outer_resid_K.median()),
                      resid_unconv_max=float(un.outer_resid_K.max()), resid_unconv_p90=float(un.outer_resid_K.quantile(0.9)),
                      cap_days_frac=float((E.outer_iters >= 5).mean()))
    return out


def sci(v, digits=1):
    s = "%.*e" % (digits, v); m, e = s.split("e")
    return "%s x 10^%d" % (m, int(e))


def table5_rows(T5):
    G = ["WDF25", "WDF23", "BRA", "MAR"]
    return [["Simulations"] + ["%d" % T5[s]["n"] for s in G],
            ["Coupling-loop passes per day, mean (max 5)"] + ["%.1f" % T5[s]["passes_mean"] for s in G],
            ["Days converged to 0.05 K within five passes"] + ["%.2f" % T5[s]["conv_frac"] for s in G],
            ["Remaining change on unconverged days, median (max) (K)"] + ["%.2f (%.1f)" % (T5[s]["resid_unconv_median"], T5[s]["resid_unconv_max"]) for s in G],
            ["Canopy energy-balance residual, max (W m-2)"] + [sci(T5[s]["canopy_resid"]) for s in G],
            ["Seasonal soil-water residual, max (mm)"] + ["%.3f" % T5[s]["wb"] for s in G],
            ["Soil-heat storage residual per sub-step, max (J m-2)"] + [sci(T5[s]["heat"]) for s in G],
            ["Organ-state carbon closure, max (g DM m-2)"] + [sci(T5[s]["carbon"]) for s in G],
            ["Nitrogen-payment identity, max (g N m-2)"] + [sci(T5[s]["npay"]) for s in G],
            ["Plant nitrogen balance per day, max (g N m-2)"] + [sci(T5[s]["nbal"]) for s in G],
            ["Grain mineral balance, max (mg m-2)"] + [sci(T5[s]["mineral"]) for s in G],
            ["Richards-solver failed days"] + ["%d" % T5[s]["richards"] for s in G]]


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="'T4' (Table 4 and the text numbers; needs the observation data) or 'T5' (Table 5; archive only)")
    a = ap.parse_args()
    config.require(config.REF_OUTPUTS, "reference-outputs archive")
    out = config.out_subdir("part2_tables")
    if a.only != "T4":
        T5 = table5(config.REF_OUTPUTS)
        json.dump(T5, open(os.path.join(out, "table5_numerical_checks.json"), "w"), indent=1)
        with open(os.path.join(out, "table5_numerical_checks.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["Check"] + [SITE_TITLES[g] for g in ("WDF25", "WDF23", "BRA", "MAR")])
            for r in table5_rows(T5):
                w.writerow(r)
    if a.only == "T5":
        print("tables written to", out)
        return
    PF.init(config.REF_OUTPUTS)
    N = numbers()
    json.dump(N, open(os.path.join(out, "part2_numbers.json"), "w"), indent=1, default=lambda x: None if isnan(x) else x)
    with open(os.path.join(out, "table4_agreement_by_site_variable.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["site", "quantity", "n", "obs_mean", "mod_mean", "bias", "rmse", "nrmse_pct", "frac_within_1sd_or_15pct", "ratio_mod_obs"])
        for r in N["table3"]:
            w.writerow([r["site"], r["quantity"], r["n"], r["obs_mean"], r["mod_mean"], r["bias"], r["rmse"], r["nrmse"], r["frac"], r["ratio"]])
    with open(os.path.join(out, "table4_series_statistics.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["site", "quantity", "n", "obs_mean", "bias", "rmse", "r"])
        for k, v in N["series"].items():
            site, q = k.split("|")
            w.writerow([site, q, v.get("n"), v.get("obs_mean"), v.get("bias"), v.get("rmse"), v.get("r")])
    print("tables written to", out)


if __name__ == "__main__":
    main()
