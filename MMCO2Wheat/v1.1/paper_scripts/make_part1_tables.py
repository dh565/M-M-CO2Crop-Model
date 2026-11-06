#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Part 1, Table 5 (normalized sensitivity indices of the perturbed input parameters) and the numbers quoted in
Sects. 13.3 and 14, from SENS_DIR/all_metrics.json (run_sensitivity.py).

Writes to OUT_DIR/part1_tables/:
  table5_parameter_sensitivity.csv   S_Y,P = (dY/Y) / (dP/P), mean of the -20 % and +20 % perturbations, per parameter
                                     and output (canopy - air temperature: dY in K per unit dP/P); the yield change of
                                     the two perturbations and the change of the relative CO2 effect on yield
  part1_numbers.txt                  the printed report (CO2 chain, elasticities, drivers, CO2 x water / N / T,
                                     parameter indices, inert-parameter check)
"""
import csv
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import config

PARAMS = ["SLA", "plant_density", "leaf_width", "root_growth_rate", "N_max_ref", "fruiting_efficiency", "TKW_potential", "chaff_share", "N_mineral_rate", "u_max_fw", "k_N_canopy", "N_soil_init", "SLN_min", "thermal_time_thresholds"]
SOUT = ["yield_g_m2", "kernels_m2", "AGB_g_m2", "LAI_peak", "ET_mm", "grain_N_pct", "N_uptake_gN_m2", "grain_Zn_mg_kg", "grain_Fe_mg_kg", "dT_midday_K"]


def report(M, fh):
    ctrl = M["ctrl"]
    P = lambda *a: print(*a, file=fh)

    def rel(n, q): return 100.0 * (M[n][q] / ctrl[q] - 1.0) if n in M else float('nan')
    def dif(n, q): return M[n][q] - ctrl[q] if n in M else float('nan')
    P("runs", len(M), "missing:", [n for n in ("co2_280",) if n not in M])
    P("ctrl:", {k: round(v, 3) if isinstance(v, float) else v for k, v in ctrl.items()})
    P("\n== CO2 chain (vs 420)")
    for q in ("An_ground_daylight_umol_m2_s", "gsw_day_mol_m2_s", "Tr_mm", "ET_mm", "LAI_peak", "AGB_g_m2", "kernels_m2", "TKW_g", "yield_g_m2", "N_uptake_gN_m2", "grain_N_pct", "grain_protein_pct", "grain_Zn_mg_kg", "grain_Fe_mg_kg", "WUE_yield_per_ET_g_mm"):
        P("  %-30s ctrl %9.3f | " % (q, ctrl[q]) + " ".join("%s %+.1f%% (%.3f)" % (c, rel("co2_%d" % c, q), M.get("co2_%d" % c, {}).get(q, float('nan'))) for c in (280, 350, 550, 700, 850)))
    P("  dT_midday_K ctrl %.3f | " % ctrl["dT_midday_K"] + " ".join("%s %+.2f K" % (c, dif("co2_%d" % c, "dT_midday_K")) for c in (280, 350, 550, 700, 850)))
    dl = math.log(550 / 420)
    P("\n== elasticities 420->550")
    for q in ("An_ground_daylight_umol_m2_s", "gsw_day_mol_m2_s", "Tr_mm", "LAI_peak", "AGB_g_m2", "kernels_m2", "TKW_g", "yield_g_m2", "N_uptake_gN_m2", "grain_N_pct", "grain_Zn_mg_kg", "grain_Fe_mg_kg", "WUE_yield_per_ET_g_mm"):
        P("  %-30s %+.2f" % (q, math.log(M["co2_550"][q] / ctrl[q]) / dl if ctrl[q] else 0))
    P("  dT K per ln CO2 %+.2f  (per doubling %+.2f)" % ((M["co2_550"]["dT_midday_K"] - ctrl["dT_midday_K"]) / dl, (M["co2_550"]["dT_midday_K"] - ctrl["dT_midday_K"]) / dl * math.log(2)))
    P("\n== drivers")
    for n in ("T-2", "T+2", "T+4", "RH-15", "RH+15", "RADx0.8", "RADx1.2", "IRRx0.50", "IRRx0.75", "IRRx1.25", "IRRx1.50", "FERTx0.00", "FERTx0.25", "FERTx0.50", "FERTx1.50", "FERTx2.00"):
        if n not in M: P("  ", n, "MISSING"); continue
        P("  %-10s yield %+5.1f AGB %+5.1f ET %+5.1f Tr %+5.1f gs %+5.1f An %+5.1f LAI %+5.1f prot %+5.1f Nupt %+5.1f TKW %.1f dT %+.2f K anth %+d d mat %+d d" % (
            n, rel(n, "yield_g_m2"), rel(n, "AGB_g_m2"), rel(n, "ET_mm"), rel(n, "Tr_mm"), rel(n, "gsw_day_mol_m2_s"), rel(n, "An_ground_daylight_umol_m2_s"), rel(n, "LAI_peak"), rel(n, "grain_protein_pct"), rel(n, "N_uptake_gN_m2"), M[n]["TKW_g"], dif(n, "dT_midday_K"), M[n]["anthesis_day"] - ctrl["anthesis_day"], M[n]["maturity_day"] - ctrl["maturity_day"]))
    P("\n== CO2 x water / N / T: 550 effect (%) on yield, AGB, ET, grain N")

    def eff(e, a, q): return 100 * (M[e][q] / M[a][q] - 1) if e in M and a in M and M[a][q] else float("nan")
    for lab, pairs in (("water", [("co2_550_IRRx%.2f" % f, "IRRx%.2f" % f if f != 1 else "ctrl") for f in (0.5, 0.75, 1.0, 1.25)]),
                       ("water700", [("co2_700_IRRx%.2f" % f, "IRRx%.2f" % f if f != 1 else "ctrl") for f in (0.5, 1.0)]),
                       ("N", [("co2_550_FERTx%.2f" % f, "FERTx%.2f" % f if f != 1 else "ctrl") for f in (0.0, 0.25, 0.5, 1.0, 1.5)]),
                       ("T", [("co2_550_T%+d" % k, "T%+d" % k if k != 0 else "ctrl") for k in (-2, 0, 2, 4)])):
        for e, a in pairs:
            P("  %-8s %-22s yield %+5.1f AGB %+5.1f ET %+5.1f grainN %+5.1f prot %+5.1f" % (lab, e, eff(e, a, "yield_g_m2"), eff(e, a, "AGB_g_m2"), eff(e, a, "ET_mm"), eff(e, a, "grain_N_pct"), eff(e, a, "grain_protein_pct")))
    P("\n== parameter indices (mean of -20/+20), and yield -20/+20 %")

    def S(pn, q, f):
        n = "P_%s_x%.1f" % (pn, f)
        if n not in M: return float('nan')
        if q == "dT_midday_K": return (M[n][q] - ctrl[q]) / (f - 1)
        return ((M[n][q] - ctrl[q]) / ctrl[q]) / (f - 1)

    def Seff(pn, f):
        a = "P_%s_x%.1f" % (pn, f); e = "P550_%s_x%.1f" % (pn, f)
        if a not in M or e not in M: return float('nan')
        return ((M[e]["yield_g_m2"] / M[a]["yield_g_m2"] - 1) - (M["co2_550"]["yield_g_m2"] / ctrl["yield_g_m2"] - 1)) / (f - 1.0)
    P("  %-24s " % "param" + " ".join("%8s" % q[:8] for q in SOUT) + "   yld-20 yld+20  Seff-20 Seff+20")
    rows = []
    for pn in PARAMS:
        idx = [np.nanmean([S(pn, q, 0.8), S(pn, q, 1.2)]) for q in SOUT]
        rows.append([pn] + idx + [rel("P_%s_x0.8" % pn, "yield_g_m2"), rel("P_%s_x1.2" % pn, "yield_g_m2"), Seff(pn, 0.8), Seff(pn, 1.2)])
        P("  %-24s " % pn + " ".join("%+8.2f" % v for v in idx) + "   %+5.1f %+5.1f   %+.2f %+.2f" % tuple(rows[-1][-4:]))
    P("\n== inert check: params whose +-20 % runs change no output by > 1e-9")
    inert = {}
    for pn in PARAMS:
        mx = 0
        for f in (0.8, 1.2):
            n = "P_%s_x%.1f" % (pn, f)
            if n in M: mx = max(mx, max(abs(M[n][q] - ctrl[q]) for q in ctrl if isinstance(ctrl[q], (int, float))))
        inert[pn] = mx
        P("  %-24s max |dY| %.3g" % (pn, mx))
    return rows, inert


def main():
    M = json.load(open(config.require(os.path.join(config.SENS_DIR, "all_metrics.json"), "all_metrics.json of the sensitivity runs")))
    out = config.out_subdir("part1_tables")
    with open(os.path.join(out, "part1_numbers.txt"), "w") as fh:
        rows, inert = report(M, fh)
    with open(os.path.join(out, "table5_parameter_sensitivity.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parameter"] + ["S_%s" % q for q in SOUT] + ["yield_change_pct_x0.8", "yield_change_pct_x1.2", "dCO2effect_per_dP_x0.8", "dCO2effect_per_dP_x1.2", "max_abs_output_change_inert_check"])
        for r in rows:
            w.writerow(r + [inert[r[0]]])
    print("tables written to", out)


if __name__ == "__main__":
    main()
