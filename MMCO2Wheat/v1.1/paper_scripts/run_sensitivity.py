#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Part 1, Sect. 14: the one-at-a-time sensitivity runs of M&M-CO2Wheat v1.0 on the WheatDryFACE 2024/25 cv. Ruta
ambient-CO2 well-watered deck (the reference simulation).  The control input is the archive deck
REF_OUTPUTS/WDF25_Ruta_aCO2_WW/input; each run is a copy of it with ONE edit, run with the model executable
(MODEL_EXE, built from the code package), and reduced to metrics.json.  All runs are collected in SENS_DIR/all_metrics.json,
the input of make_part1_figures.py (Figs. 2-5) and make_part1_tables.py (Table 5).

Experiments (one-at-a-time):
  E1 drivers  : CO2 (280..850), air temperature (+-2, +4 K on Tmax and Tmin), relative humidity (-15/+15 points, x0.9, x1.1,
                capped at 100), radiation (x0.8..x1.2), water (irrigation x0.5, x0.75, x1.25, x1.5), N fertilizer (x0..x2)
  E2 CO2 curve: 280, 350, 420, 550, 700, 850
  E3 550 vs 420 across water (irrigation x0.5, 0.75, 1, 1.25), N (fertilizer x0, 0.25, 0.5, 1, 1.5), temperature (-2, 0, +2, +4 K)
  E4 parameters x0.8, 0.9, 1.1, 1.2 (and x0.8, 1.2 at 550): SLA, plant density, leaf width, root growth rate, N_max_ref,
                fruiting efficiency, TKW potential, chaff share, N_mineral_rate, u_max_fw, k_N_canopy, N_soil_init, SLN_min,
                thermal-time thresholds (all stages)

    MODEL_EXE=/path/to/mmco2crop python3 run_sensitivity.py [run names...]      (125 runs, about 1-2 min each)
"""
import csv
import json
import math
import os
import re
import shutil
import statistics as st
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import config

EXE = config.MODEL_EXE
CTRL = os.path.join(config.REF_OUTPUTS, "WDF25_Ruta_aCO2_WW", "input")
ROOT = os.path.join(config.SENS_DIR, "runs")

# driver_main.txt records edited (record N is line N of the file; see INPUT_GUIDE.md of the code package)
REC_SLA, REC_DENSITY, REC_LEAF_WIDTH, REC_CO2, REC_ROOT_RATE = 7, 8, 16, 22, 28


# ---------------------------------------------------------------- editors
def edit_driver_main(path, line_idx0, new_value):
    L = open(path).read().split("\n")
    s = L[line_idx0]; m = re.match(r"^(\s*)(\S+)(.*)$", s)
    L[line_idx0] = "%s%s%s" % (m.group(1), new_value, m.group(3))
    open(path, "w").write("\n".join(L))

def driver_main_index(path, key):
    """0-based line index of a driver_main.txt record: record number (int, record N is line N) or a comment key (str)."""
    if isinstance(key, int):
        return key - 1
    L = open(path).read().split("\n")
    for i, s in enumerate(L):
        if key in s: return i
    raise KeyError(key)

def edit_csv_param(path, pname, factor=None, value=None):
    L = open(path).read().split("\n"); done = False
    for i, s in enumerate(L):
        if s.startswith("!") or not s.strip(): continue
        parts = s.split(",")
        if len(parts) > 1 and parts[1].strip() == pname:
            v = float(parts[0]); nv = value if value is not None else v * factor
            parts[0] = ("%.6g" % nv); L[i] = ",".join(parts); done = True; break
    assert done, pname
    open(path, "w").write("\n".join(L))

def edit_meteo(path, fn):
    rows = list(csv.DictReader(open(path))); hdr = list(rows[0].keys())
    for r in rows: fn(r)
    w = csv.DictWriter(open(path, "w", newline=""), fieldnames=hdr); w.writeheader(); w.writerows(rows)

def edit_fert(path, factor):
    L = open(path).read().split("\n"); out = []; n = None
    for s in L:
        if s.startswith("!") or not s.strip(): out.append(s); continue
        if n is None: n = int(s.strip()); out.append(s); continue
        p = s.split(","); p[1] = "%.4g" % (float(p[1]) * factor); out.append(",".join(p))
    open(path, "w").write("\n".join(out))

def edit_phen_scale(path, factor):
    L = open(path).read().split("\n")
    for i, s in enumerate(L):
        p = s.split()
        if len(p) >= 2 and p[1] in ("Emergence", "3_leaves", "Tillering", "Booting", "Heading", "Anthesis", "Grain_filling", "Maturity"):
            L[i] = "%.1f   \t%s" % (float(p[0]) * factor, p[1])
    open(path, "w").write("\n".join(L))

def temp_shift(dK):
    def f(r):
        r["Tair_max (C) "] = "%.3f" % (float(r["Tair_max (C) "]) + dK); r["Tair_min (C)"] = "%.3f" % (float(r["Tair_min (C)"]) + dK)
    return f
def rh_scale(fac):
    def f(r):
        for k in ("RHmax (%)", "RHmin (%)"): r[k] = "%.3f" % min(100.0, float(r[k]) * fac)
    return f
def rad_scale(fac):
    def f(r): r["Rad (MJ m-2 d-1)"] = "%.4f" % (float(r["Rad (MJ m-2 d-1)"]) * fac)
    return f
def irr_scale(fac):
    def f(r): r["irrigation"] = "%.3f" % (float(r["irrigation"]) * fac)
    return f

# ---------------------------------------------------------------- experiment registry: name -> list of (kind, args)
RUNS = {}
def add(name, *edits): RUNS[name] = edits
add("ctrl")
for c in (280, 350, 550, 700, 850): add("co2_%d" % c, ("dm", REC_CO2, c))
for dK in (-2, 2, 4): add("T%+d" % dK, ("met", temp_shift(dK)))
for f in (0.9, 1.1): add("RHx%.1f" % f, ("met", rh_scale(f)))
def rh_shift(dp):
    def f(r):
        for k in ("RHmax (%)", "RHmin (%)"): r[k] = "%.3f" % min(100.0, max(1.0, float(r[k]) + dp))
    return f
for dp in (-15, 15): add("RH%+d" % dp, ("met", rh_shift(dp)))
for f in (0.8, 0.9, 1.1, 1.2): add("RADx%.1f" % f, ("met", rad_scale(f)))
for f in (0.5, 0.75, 1.25, 1.5): add("IRRx%.2f" % f, ("met", irr_scale(f)))
for f in (0.0, 0.25, 0.5, 1.5, 2.0): add("FERTx%.2f" % f, ("fert", f))
# E3: 550 x water / N / T
for f in (0.5, 0.75, 1.0, 1.25): add("co2_550_IRRx%.2f" % f, ("dm", REC_CO2, 550), ("met", irr_scale(f)))
for f in (0.0, 0.25, 0.5, 1.0, 1.5): add("co2_550_FERTx%.2f" % f, ("dm", REC_CO2, 550), ("fert", f))
for f in (0.5, 1.0): add("co2_700_IRRx%.2f" % f, ("dm", REC_CO2, 700), ("met", irr_scale(f)))
for dK in (-2, 0, 2, 4): add("co2_550_T%+d" % dK, ("dm", REC_CO2, 550), ("met", temp_shift(dK)))
# E4: parameters +-20 % (and +-10 %)
PARAMS = [("SLA", ("dm", REC_SLA)), ("plant_density", ("dm", REC_DENSITY)), ("leaf_width", ("dm", REC_LEAF_WIDTH)),
          ("root_growth_rate", ("dm", REC_ROOT_RATE)), ("N_max_ref", ("dn", "N_max_ref")), ("fruiting_efficiency", ("cy", "fruiting_efficiency")),
          ("TKW_potential", ("cy", "TKW_potential")), ("chaff_share", ("cy", "chaff_share")), ("N_mineral_rate", ("dn", "N_mineral_rate")),
          ("u_max_fw", ("dn", "u_max_fw")), ("k_N_canopy", ("dn", "k_N_canopy")), ("N_soil_init", ("dn", "N_soil_init")),
          ("SLN_min", ("dn", "SLN_min")), ("thermal_time_thresholds", ("phen", None)), ("leaf_width_dummy", None)]
PARAMS = [p for p in PARAMS if p[1] is not None]
for pn, spec in PARAMS:
    for f in (0.8, 0.9, 1.1, 1.2):
        add("P_%s_x%.1f" % (pn, f), (spec[0] + "_scale", spec[1], f))
    for f in (0.8, 1.2):
        add("P550_%s_x%.1f" % (pn, f), ("dm", REC_CO2, 550), (spec[0] + "_scale", spec[1], f))

def dm_value(key):
    L = open(os.path.join(CTRL, "driver_main.txt")).read().split("\n")
    return float(L[driver_main_index(os.path.join(CTRL, "driver_main.txt"), key)].split()[0])

def prepare(name):
    d = os.path.join(ROOT, name)
    if os.path.exists(os.path.join(d, "metrics.json")): return d
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d); shutil.copytree(CTRL, os.path.join(d, "input")); os.makedirs(os.path.join(d, "output"))
    inp = os.path.join(d, "input")
    for e in RUNS[name]:
        kind = e[0]
        if kind == "dm": edit_driver_main(os.path.join(inp, "driver_main.txt"), driver_main_index(os.path.join(inp, "driver_main.txt"), e[1]), e[2])
        elif kind == "dm_scale":
            p = os.path.join(inp, "driver_main.txt"); i = driver_main_index(p, e[1]); v = float(open(p).read().split("\n")[i].split()[0]); edit_driver_main(p, i, "%.6g" % (v * e[2]))
        elif kind == "met": edit_meteo(os.path.join(inp, "input_daily_meteo.csv"), e[1])
        elif kind == "fert": edit_fert(os.path.join(inp, "N_fertilizer.csv"), e[1])
        elif kind == "dn_scale": edit_csv_param(os.path.join(inp, "driver_N.csv"), e[1], factor=e[2])
        elif kind == "cy_scale": edit_csv_param(os.path.join(inp, "cultivar_yield.csv"), e[1], factor=e[2])
        elif kind == "phen_scale": edit_phen_scale(os.path.join(inp, "phenology_cultivar.txt"), e[2])
        else: raise ValueError(kind)
    return d

def metrics(d):
    o = os.path.join(d, "output")
    def rows(f): return list(csv.DictReader(open(os.path.join(o, f))))
    g = rows("daily_growth.csv"); yc = rows("daily_yield_components.csv"); fl = rows("daily_fluxes.csv"); ge = rows("daily_gas_exchange.csv")
    nn = rows("daily_nitrogen.csv"); eb = rows("daily_eb_actual.csv"); ph = rows("daily_phenology.csv")
    F = lambda r, k: float(r[k])
    m = dict(yield_g_m2=F(yc[-1], "grain_yield_g_m2"), AGB_g_m2=F(g[-1], "AGB_g_m2"), kernels_m2=F(yc[-1], "kernel_number_m2"), TKW_g=F(yc[-1], "TKW_g"),
             HI=F(yc[-1], "harvest_index"), LAI_peak=max(F(r, "LAI_m2_m2") for r in g),
             ET_mm=sum(F(r, "ET_total_act_mm") for r in eb), Tr_mm=sum(F(r, "Tr_act_mm") for r in eb), Esoil_mm=sum(F(r, "Esoil_act_mm") for r in eb),
             An_season_gC=sum(F(r, "An_gCO2_m2_d") for r in ge) * 12.0 / 44.0,
             grain_protein_pct=F(nn[-1], "grain_protein_pct"), grain_N_pct=100 * F(nn[-1], "grainN_conc_gN_gDM"), N_grain_gN_m2=F(nn[-1], "N_grain_gN_m2"),
             N_uptake_gN_m2=sum(F(r, "N_soil_uptake_day_gN_m2") for r in nn))
    # midday canopy - air and gsw over the canopy-closed window (anthesis-45 d .. anthesis), from the hourly EB file
    anth = [int(float(r["day"])) for r in ph if r["stage_name"].strip() == "Anthesis"]
    a = anth[0] if anth else int(0.6 * len(g))
    H = csv.DictReader(open(os.path.join(o, "hourly_energy_balance.csv")), skipinitialspace=True)
    dT, Tc = [], []
    for r in H:
        day = int(float(r["day"])); hr = float(r["solar_hour"])
        if a - 45 <= day <= a and 10 <= hr <= 14:
            dT.append(float(r["Tcan_C"]) - float(r["Tair_C"])); Tc.append(float(r["Tcan_C"]))
    m["dT_midday_K"] = st.mean(dT) if dT else float("nan"); m["Tcan_midday_C"] = st.mean(Tc) if Tc else float("nan")
    gs = [F(r, "gsw_day_mol_m2_s") for r in ge if a - 45 <= int(float(r["day"])) <= a]; m["gsw_day_mol_m2_s"] = st.mean(gs)
    anp = [F(r, "An_ground_umol_m2_s_daylight") for r in ge if a - 45 <= int(float(r["day"])) <= a]; m["An_ground_daylight_umol_m2_s"] = st.mean(anp)
    m["WUE_yield_per_ET_g_mm"] = m["yield_g_m2"] / m["ET_mm"]; m["anthesis_day"] = a; m["maturity_day"] = len(g)
    # grain mineral concentrations at harvest (mg kg-1) and grain mineral content (mg m-2)
    try:
        mi = rows("daily_minerals.csv"); last = max(int(r["day"]) for r in mi)
        for r in mi:
            if int(r["day"]) == last:
                m["grain_%s_mg_kg" % r["mineral"].strip()] = F(r, "grain_mg_kg"); m["grain_%s_mg_m2" % r["mineral"].strip()] = F(r, "grain_mg_m2")
    except FileNotFoundError:
        pass
    return m

def run(name):
    d = prepare(name)
    if os.path.exists(os.path.join(d, "metrics.json")): return name, "cached"
    with open(os.path.join(d, "run.log"), "w") as lf:
        rc = subprocess.call([EXE], cwd=d, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0: return name, "exit %d" % rc
    m = metrics(d); json.dump(m, open(os.path.join(d, "metrics.json"), "w"), indent=1)
    # trim the big hourly files to save disk
    for f in ("substep_diag_all.csv", "hydro_wetfront_trace.csv", "hourly_esoil_pairing.csv", "hourly_node_pairing.csv", "hourly_radiation.csv", "hourly_soil_heat.csv", "diurnal_energy_balance.csv"):
        p = os.path.join(d, "output", f)
        if os.path.exists(p): os.remove(p)
    return name, "ok"


if __name__ == "__main__":
    config.require(CTRL, "control deck input (REF_OUTPUTS/WDF25_Ruta_aCO2_WW/input)")
    if not EXE or not os.path.isfile(EXE):
        raise SystemExit("set MODEL_EXE to the model executable built from the code package")
    os.makedirs(ROOT, exist_ok=True)
    names = sys.argv[1:] or list(RUNS)
    print(len(names), "runs")
    with ThreadPoolExecutor(max_workers=int(os.environ.get("SENS_WORKERS", "2"))) as ex:
        for name, status in ex.map(run, names):
            print(name, status, flush=True)
    M = {n: json.load(open(os.path.join(ROOT, n, "metrics.json"))) for n in RUNS if os.path.exists(os.path.join(ROOT, n, "metrics.json"))}
    dest = os.path.join(config.SENS_DIR, "all_metrics.json")
    json.dump(M, open(dest, "w"), indent=1)
    print("saved", len(M), "runs to", dest)
