#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Deck registry, model-output readers and observation registry for the FACE evaluation (Part 2).

No model is run here.  Every model number is read from the reference-outputs archive
(config.REF_OUTPUTS, one folder per deck with input/ and output/); every observation is read
from the treatment files of the three experiments or from the derived observation tables
described in observations/README_observations.md.
"""
import csv
import datetime as dt
import math
import os
import sys
from collections import OrderedDict, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
from wdf_ingest import W25, W23

NAN = float("nan")

# ----------------------------------------------------------------- fixed basis constants
# Soil layer thicknesses / boundaries of the 10-layer profile used when a deck's own
# driver_swtr.txt cannot be read.
THICK = [0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.18, 0.18, 0.18, 0.18]
BNDS = [0, 0.02, 0.06, 0.12, 0.20, 0.32, 0.48, 0.66, 0.84, 1.02, 1.20]
MIDP = [(BNDS[i] + BNDS[i + 1]) / 2 for i in range(10)]

SIGN_FLOOR = 0.05         # smallest simulated effect counted as a resolved sign (|R_model - 1|)
SIGN_RESOLVE_PCT = 0.15   # observed effect counted as resolved where no sd exists (|R_obs - 1|)

# ----------------------------------------------------------------- deck registry (ORDERED)
SITES = ["WDF25", "WDF23", "BRAUNSCHWEIG", "MARICOPA"]
DECKS = OrderedDict()   # deck -> dict(site, treatment, cultivar, co2, water, nlevel, season)


def _reg(deck, site, treatment, cultivar, co2, water, nlevel, season):
    DECKS[deck] = dict(deck=deck, site=site, treatment=treatment, cultivar=cultivar,
                       co2=co2, water=water, nlevel=nlevel, season=season)


for cv in ("Ruta", "Zahir"):
    for co2 in ("aCO2", "eCO2"):
        for wat in ("WW", "WL"):
            _reg("WDF25_%s_%s_%s" % (cv, co2, wat), "WDF25", "%s %s %s" % (cv, co2, wat),
                 cv, co2, wat, "ample (212 kgN/ha)", "2024/25")
for cv in ("Ruta", "Zahir"):
    for co2 in ("aCO2", "eCO2"):
        _reg("WDF23_%s_%s" % (cv, co2), "WDF23", "%s %s" % (cv, co2), cv, co2, "WW",
             "low (55 kgN/ha)", "2022/23")
# Braunschweig 2013/14: treatment numbers of the data set (Manderscheid et al., 2020), N level, CO2, N fertilizer
_BRA = [("101", "N1", "aCO2", "40 kgN/ha"), ("102", "N2", "aCO2", "180"),
        ("103", "N3", "aCO2", "320"), ("107", "N4", "aCO2", "180 split"),
        ("104", "N1", "eCO2", "40"), ("105", "N2", "eCO2", "180"),
        ("106", "N3", "eCO2", "320"), ("110", "N4", "eCO2", "180 split")]
for trno, nl, co2, fert in _BRA:
    _reg("BRAUNSCHWEIG_%s" % trno, "BRAUNSCHWEIG", "%s %s %s" % (nl, co2, fert),
         "Batis", co2, "WW", "%s (%s)" % (nl, fert), "2013/14")
# Braunschweig 2014/15: season N totals from each deck's own N_fertilizer.csv
# (N1 = 35, N2 = 200, N3 = 320, N4 = 202 as a 50/50 urea + ammonium-nitrate split)
_BRA15 = [("113", "N1", "aCO2", "35 kgN/ha"), ("114", "N2", "aCO2", "200"),
          ("115", "N3", "aCO2", "320"), ("119", "N4", "aCO2", "202 split"),
          ("116", "N1", "eCO2", "35"), ("117", "N2", "eCO2", "200"),
          ("118", "N3", "eCO2", "320"), ("122", "N4", "eCO2", "202 split")]
for trno, nl, co2, fert in _BRA15:
    _reg("BRAUNSCHWEIG_%s" % trno, "BRAUNSCHWEIG", "%s %s %s" % (nl, co2, fert),
         "Batis", co2, "WW", "%s (%s)" % (nl, fert), "2014/15")
# Maricopa: treatment numbers (TRNO) of the Kimball et al. (2017) archive
_MAR = [("905", "Dry", "aCO2", "High_N", "1993/94"), ("906", "Wet", "aCO2", "High_N", "1993/94"),
        ("907", "Dry", "eCO2", "High_N", "1993/94"), ("908", "Wet", "eCO2", "High_N", "1993/94"),
        ("910", "Wet", "aCO2", "High_N", "1995/96"), ("912", "Wet", "eCO2", "High_N", "1995/96"),
        ("913", "Wet", "aCO2", "Low_N", "1996/97"), ("914", "Wet", "aCO2", "High_N", "1996/97"),
        ("915", "Wet", "eCO2", "Low_N", "1996/97"), ("916", "Wet", "eCO2", "High_N", "1996/97")]
for trno, wat, co2, nl, seas in _MAR:
    _reg("MARICOPA_%s" % trno, "MARICOPA", "%s %s %s %s" % (wat, nl, co2, seas),
         "Yecora Rojo", co2, "WW" if wat == "Wet" else "WL", nl, seas)

FACE_DECKS = list(DECKS)


# ----------------------------------------------------------------- generic helpers
def fl(x, default=NAN):
    try:
        v = float(str(x).strip())
        return v if v == v else default
    except Exception:
        return default


def isnan(x):
    return x is None or (isinstance(x, float) and x != x)


def rows_of(path):
    if not os.path.isfile(path):
        return []
    with open(path, encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def table_rows(path):
    """Rows of a shipped observation table, skipping its '#' provenance header."""
    if not os.path.isfile(path):
        return []
    with open(path, encoding="utf-8-sig") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def bra_obs_path(deck):
    """The Braunschweig observation file of one deck (Manderscheid et al., 2020)."""
    return os.path.join(config.BRAUNSCHWEIG_OBS_DIR, "%s_obs.csv" % deck)


# ================================================================= model-side extraction
# peak LAI = max over daily_growth LAI_m2_m2; AGB / yield / organ masses = the last row of
# daily_growth / daily_yield_components; leaf total = W_green + W_dead (dead lamina is retained
# in the simulated above-ground biomass, as it is in the harvested biomass).
MQ = ["LAI_peak", "LAI_peak_date", "AGB_harvest", "grain_yield", "TKW", "kernel_number", "HI",
      "grain_protein_pct", "grain_N_pct", "stem_N_pct", "leaf_N_pct", "height_max", "chaff",
      "W_stem", "W_leaf_tot", "spike", "stem_WSC_mg_g"]


def _kv(path):
    out = {}
    if not os.path.isfile(path):
        return out
    for ln in open(path):
        if ln.startswith("#") or "," not in ln:
            continue
        t = ln.split(",")
        out[t[0].strip()] = t[1].strip()
    return out


def deck_path(pkg, deck):
    """Deck root inside the archive: <pkg>/<deck>/ (a site-rooted layout <pkg>/<SITE>/<deck>/ is also accepted)."""
    a = os.path.join(pkg, DECKS[deck]["site"], deck)
    if os.path.isdir(os.path.join(a, "output")):
        return a
    return os.path.join(pkg, deck)


def deck_dir(pkg, deck):
    return os.path.join(deck_path(pkg, deck), "output")


def isod(r):
    return "%04d-%02d-%02d" % (fl(r["year"]), fl(r["month"]), fl(r["dom"]))


class DeckModel(object):
    """The model series and scalars one deck provides, read from its output files."""

    def __init__(self, pkg, deck):
        self.deck, self.pkg = deck, pkg
        d = deck_dir(pkg, deck)
        self.dir = d
        self.g = rows_of(os.path.join(d, "daily_growth.csv"))
        self.y = rows_of(os.path.join(d, "daily_yield_components.csv"))
        self.n = rows_of(os.path.join(d, "daily_nitrogen.csv"))
        self.p = rows_of(os.path.join(d, "daily_phenology.csv"))
        self.sw = rows_of(os.path.join(d, "daily_swc_layers.csv"))
        self.wsc = rows_of(os.path.join(d, "daily_stem_wsc.csv"))
        self.res = rows_of(os.path.join(d, "daily_carbon_reserve.csv"))
        self.fld = rows_of(os.path.join(d, "flag_leaf_diag.csv"))
        self.min = rows_of(os.path.join(d, "daily_minerals.csv"))
        self.law = _kv(os.path.join(d, "stem_wsc_law.csv"))
        self.ok = bool(self.g and self.p)
        self.date2day = {isod(r): int(fl(r["day"])) for r in self.p}
        self.day2date = {int(fl(r["day"])): isod(r) for r in self.p}
        self._swc = {int(fl(r["day"])): [fl(r["layer%d_mm" % (i + 1)]) / (THICK[i] * 1000.0)
                                         for i in range(10)] for r in self.sw}
        self.v = self._scalars() if self.ok else {}
        self.stages = self._stages() if self.ok else {}

    # --- series
    def series(self, rows, col):
        return [(int(fl(r["day"])), fl(r[col])) for r in rows if col in r]

    def lai(self):
        return self.series(self.g, "LAI_m2_m2")

    def agb(self):
        return self.series(self.g, "AGB_g_m2")

    def swc_at(self, day, z):
        if day not in self._swc:
            return None
        v = self._swc[day]
        if z <= MIDP[0]:
            return v[0]
        if z >= MIDP[-1]:
            return v[-1]
        for i in range(9):
            if MIDP[i] <= z <= MIDP[i + 1]:
                f = (z - MIDP[i]) / (MIDP[i + 1] - MIDP[i])
                return v[i] + f * (v[i + 1] - v[i])
        return v[-1]

    def swc_series(self, z):
        return [(d, self.swc_at(d, z)) for d in sorted(self._swc)]

    def organ_at_day(self, day):
        """(stem, leaf_tot, spike, chaff, grain, AGB) at one model day."""
        gi = {int(fl(r["day"])): r for r in self.g}
        yi = {int(fl(r["day"])): r for r in self.y}
        if day not in gi:
            return None
        gr, yr = gi[day], yi.get(day, {})
        return dict(W_stem=fl(gr["W_stem_g_m2"]),
                    W_leaf_tot=fl(gr["W_green_g_m2"]) + fl(gr["W_dead_g_m2"]),
                    spike=fl(yr.get("W_spike_g_m2", "nan")),
                    chaff=fl(yr.get("W_chaff_g_m2", "nan")),
                    grain=fl(yr.get("W_grain_g_m2", "nan")),
                    AGB_harvest=fl(gr["AGB_g_m2"]))

    # --- scalars
    def _scalars(self):
        gl, yl = self.g[-1], (self.y[-1] if self.y else {})
        nl = self.n[-1] if self.n else {}
        laiv = [(fl(r["LAI_m2_m2"]), int(fl(r["day"]))) for r in self.g]
        laipk = max(laiv)
        wstem = fl(gl["W_stem_g_m2"])
        rend = fl(self.law.get("W_reserve_end_g_m2", "nan"))
        agb = fl(gl["AGB_g_m2"])
        gy = fl(gl["GrainYield_g_m2"])
        return dict(
            LAI_peak=laipk[0],
            LAI_peak_date=self.day2date.get(laipk[1], ""),
            AGB_harvest=agb,
            grain_yield=gy,
            TKW=fl(yl.get("TKW_g", "nan")),
            kernel_number=fl(yl.get("kernel_number_m2", "nan")),
            HI=(gy / agb if agb > 0 else NAN),
            grain_protein_pct=fl(nl.get("grain_protein_pct", "nan")),
            grain_N_pct=fl(nl.get("grainN_conc_gN_gDM", "nan")) * 100.0,
            stem_N_pct=fl(nl.get("stemN_conc_gN_gDM", "nan")) * 100.0,
            leaf_N_pct=fl(nl.get("wholeleafN_conc_gN_gDM", "nan")) * 100.0,
            height_max=max(fl(r["plant_height_m"]) for r in self.g),
            chaff=fl(yl.get("W_chaff_g_m2", "nan")),
            W_stem=wstem,
            W_leaf_tot=fl(gl["W_green_g_m2"]) + fl(gl["W_dead_g_m2"]),
            spike=fl(yl.get("W_spike_g_m2", "nan")),
            stem_WSC_mg_g=(rend / wstem * 1000.0) if wstem > 0 else NAN,
        )

    def _stages(self):
        out = OrderedDict()
        for r in self.p:
            s = r["stage_name"].strip()
            if s and s not in out:
                out[s] = isod(r)
        return out


def load_arm(pkg, decks):
    out = OrderedDict()
    for d in decks:
        m = DeckModel(pkg, d)
        if m.ok:
            out[d] = m
    return out


# ================================================================= observation side
# Model stage name -> the observation file's own stage label, per site.  Phenology is compared
# by calendar date: the dated 'Phenological stage' column at WheatDryFACE, the dated stage rows
# at Braunschweig, and the Obs_*_DAP entries + sowing date at Maricopa.
WDF25_STAGE_MAP = [("Emergence", "Emergence"), ("3_leaves", "3 leaves"),
                   ("Tillering", "Tillering"), ("Booting", "Booting"), ("Heading", "Heading"),
                   ("Anthesis", "Anthesis"), ("Grain_filling", "Grain filling start"),
                   ("Maturity", "Dry kernel")]
WDF23_STAGE_MAP = [("Emergence", "Emergence"), ("3_leaves", "2 to 3 leaves"),
                   ("Tillering", "Tillering"), ("Booting", "Booting"), ("Heading", "Heading"),
                   ("Anthesis", "Flowering"), ("Grain_filling", "Early dough"),
                   ("Maturity", "Ripening")]
DAP_STAGE_MAP = [("Emergence", "Obs_Emergence_DAP"), ("3_leaves", "Obs_3_leaves_DAP"),
                 ("Tillering", "Obs_Tillering_DAP"), ("Booting", "Obs_Booting_DAP"),
                 ("Heading", "Obs_Heading_DAP"), ("Anthesis", "Obs_Anthesis_DAP"),
                 ("Grain_filling", "Obs_Grain_filling_DAP"), ("Maturity", "Obs_Maturity_DAP")]
BRA_STAGE_MAP = [("Emergence", "emergence"), ("3_leaves", "three_leaf_M2020"),
                 ("Anthesis", "anthesis"), ("Maturity", "grain_maturity")]
# Stages with an observation whose basis is established at each site.  At Braunschweig the
# `three_leaf_M2020` row (2014-03-31) is not a 3-leaf date on the model's convention and is not
# paired; at Maricopa only anthesis and maturity are direct archive stage records.
SCORED_STAGES = {
    "WDF25": set(m for m, _ in DAP_STAGE_MAP),
    "WDF23": set(m for m, _ in DAP_STAGE_MAP),
    "BRAUNSCHWEIG": {"Emergence", "Anthesis", "Maturity"},
    "MARICOPA": {"Anthesis", "Maturity"},
}


class Obs(object):
    """obs[deck][quantity] -> dict(value, sd, basis, source, note, date)."""

    def __init__(self):
        self.d = defaultdict(dict)
        self._swc_obs = {}     # (deck, label) -> list of (date, value_vv, sd_vv)
        self.lai_ts = {}       # deck -> ([(date, value, sd)], basis label)
        self.phen = {}
        self._from_registry()
        self._wdf25_organs()
        self._wdf23_derived()
        self._bra_organs_and_N()
        self._mar_scalars()
        self._wdf_tkw_and_derived_kernels()
        self._mar_replicate_sd()
        self._collect_swc()

    def put(self, deck, q, value, sd, basis, source, note="", date=""):
        self.d[deck][q] = dict(value=value, sd=sd, basis=basis, source=source,
                               note=note, date=date)

    def get(self, deck, q):
        return self.d.get(deck, {}).get(q)

    # -- the per-treatment observation registry (harvest and canopy scalars with their sd),
    #    derived from the treatment files of the three experiments; two tables: the open sites
    #    (shipped) and WheatDryFACE (under OBS_DIR)
    def _from_registry(self):
        for p in (config.OBS_REGISTRY_OPEN, config.OBS_REGISTRY_WDF):
            if not os.path.isfile(p):
                sys.stderr.write("observation registry not found: %s\n" % p)
            for r in table_rows(p):
                dk, q = r["deck"], r["quantity"]
                if dk not in DECKS:
                    continue
                sd = fl(r["obs_sd"])
                self.put(dk, q, fl(r["obs"]), (sd if not isnan(sd) else NAN), r["basis"],
                         os.path.basename(p),
                         note=("harvest index obs = obs grain / obs AGB" if q == "HI" else ""))

    # -- WDF25 harvest organ DM, all 8 decks (t/ha x 100 -> g m-2); grain minerals; LAI series
    def _wdf25_organs(self):
        mapping = [("Stem", "W_stem"), ("Leaf", "W_leaf_tot"), ("Spike (head)", "spike"),
                   ("Chaff", "chaff"), ("Grain", "grain_yield")]
        for dk in [d for d in DECKS if DECKS[d]["site"] == "WDF25"]:
            cv, co2, wat = dk.split("_")[1:4]
            rows = W25.obs_for(cv, co2, wat)
            bio = {(r["organ"] or "").strip(): r for r in rows if r["variable"] == "Biomass"}
            tot, totv = 0.0, 0.0
            for org, q in mapping:
                r = bio.get(org)
                if not r:
                    continue
                v = fl(r["value"]) * 100.0
                s = fl(r["std"]) * 100.0
                if q in ("W_stem", "W_leaf_tot", "spike"):
                    tot += v
                    totv += (0.0 if isnan(s) else s * s)
                if self.get(dk, q) is None:
                    self.put(dk, q, v, s, "sd",
                             "WDF25/S3_%s_%s_%s.csv harvest %s Biomass block "
                             "(t ha-1 x 100)" % (cv, co2, wat, org), date=r["date"])
            if self.get(dk, "AGB_harvest") is None and tot > 0:
                self.put(dk, "AGB_harvest", tot, math.sqrt(totv), "sd",
                         "WDF25 Stem + Spike + Leaf harvest Biomass blocks (quadrature sd)")
            for mn in ("Fe", "Mg", "Ca", "Cu", "Mn", "Zn"):
                mr = [r for r in rows if r["variable"] == mn]
                if mr:
                    self.put(dk, "grain_%s" % mn, fl(mr[0]["value"]), fl(mr[0]["std"]), "sd",
                             "WDF25/S3_%s_%s_%s.csv Grain %s" % (cv, co2, wat, mn),
                             date=mr[0]["date"])
            lai = [(r["date"], fl(r["value"]), fl(r["std"])) for r in rows
                   if r["variable"] == "LAI"]
            self.lai_ts[dk] = (sorted(lai), "WDF25 campaign LAI")

    # -- WDF23: AGB and straw are derived from the harvest block (Straw + Spike)
    def _wdf23_derived(self):
        for dk in [d for d in DECKS if DECKS[d]["site"] == "WDF23"]:
            cv, co2 = dk.split("_")[1:3]
            rows = W23.obs_for(cv, co2)
            bio = {(r["organ"] or "").strip(): r for r in rows if r["variable"] == "Biomass"}
            st, sp = bio.get("Straw"), bio.get("Spike")
            if st and sp:
                v = (fl(st["value"]) + fl(sp["value"])) * 100.0
                s = math.sqrt((fl(st["std"], 0) * 100.0) ** 2 + (fl(sp["std"], 0) * 100.0) ** 2)
                self.put(dk, "AGB_harvest", v, s, "sd",
                         "DERIVED: WDF23/WDF23_%s_%s.csv harvest Straw + Spike Biomass "
                         "(t ha-1 x 100), sd in quadrature" % (cv, co2),
                         note="DERIVED obs", date=st["date"])
                self.put(dk, "straw", fl(st["value"]) * 100.0, fl(st["std"], 0) * 100.0, "sd",
                         "WDF23/WDF23_%s_%s.csv harvest Straw Biomass "
                         "(leaf + stem, t ha-1 x 100)" % (cv, co2), date=st["date"])
            gr = bio.get("Grain")
            if gr and self.get(dk, "HI") is None and self.get(dk, "AGB_harvest"):
                a = self.get(dk, "AGB_harvest")
                self.put(dk, "HI", fl(gr["value"]) * 100.0 / a["value"], NAN, "pct",
                         "DERIVED obs grain / DERIVED obs AGB", note="DERIVED on both limbs")
            for mn in ("Fe", "Mg", "Ca", "Cu", "Mn", "Zn"):
                mr = [r for r in rows if r["variable"] == mn]
                if mr:
                    self.put(dk, "grain_%s" % mn, fl(mr[0]["value"]), fl(mr[0]["std"]), "sd",
                             "WDF23/WDF23_%s_%s.csv Grain %s" % (cv, co2, mn),
                             date=mr[0]["date"])
            lai = [(r["date"], fl(r["value"]), fl(r["std"])) for r in rows
                   if r["variable"] == "LAI"]
            self.lai_ts[dk] = (sorted(lai), "WDF23 campaign LAI")

    # -- WDF25 kernel mass: least-square means by cultivar x CO2 x water of the 2024/25 season
    #    (Mulero et al., submitted; Supplementary Table S4B), no sd.  Kernel number of both
    #    WheatDryFACE seasons = harvest grain DM / TKW x 1000, sd propagated from the limbs that carry one.
    WDF25_TKW_S4B = {("Ruta", "aCO2", "WW"): 30.0, ("Ruta", "eCO2", "WW"): 29.3,
                     ("Zahir", "aCO2", "WW"): 43.0, ("Zahir", "eCO2", "WW"): 41.8,
                     ("Ruta", "aCO2", "WL"): 32.6, ("Ruta", "eCO2", "WL"): 32.9,
                     ("Zahir", "aCO2", "WL"): 44.2, ("Zahir", "eCO2", "WL"): 44.5}

    def _wdf_tkw_and_derived_kernels(self):
        for dk in [d for d in DECKS if DECKS[d]["site"] in ("WDF25", "WDF23")]:
            site = DECKS[dk]["site"]
            if site == "WDF25" and self.get(dk, "TKW") is None:
                cv, co2, wat = dk.split("_")[1:4]
                tk = self.WDF25_TKW_S4B.get((cv, co2, wat))
                if tk is not None:
                    self.put(dk, "TKW", tk, NAN, "pct",
                             "Mulero et al. (submitted), Supplementary Table S4B, least-square mean, no sd",
                             note="the simulated kernel mass of this treatment is bounded by the potential kernel mass input")
            tkr, gyr = self.get(dk, "TKW"), self.get(dk, "grain_yield")
            if tkr and gyr and self.get(dk, "kernel_number") is None and tkr["value"] > 0:
                kn = gyr["value"] / tkr["value"] * 1000.0
                rel2 = 0.0
                for o in (gyr, tkr):
                    if not isnan(o["sd"]) and o["value"] > 0:
                        rel2 += (o["sd"] / o["value"]) ** 2
                sd = kn * math.sqrt(rel2) if rel2 > 0 else NAN
                self.put(dk, "kernel_number", kn, sd, "sd" if not isnan(sd) else "pct",
                         "DERIVED: observed harvest grain DM / observed TKW x 1000 (%s; %s); "
                         "sd propagated in quadrature from the limbs that carry one"
                         % (gyr["source"][:60], tkr["source"][:60]), note="DERIVED obs")

    # -- Maricopa replicate dispersion: the by-replicate sheet of the archive gives the sd of
    #    HWAH (yield), CWAM (AGB) and LAIX; HWUM / H#UM / P#AM carry no replicate values.
    def _mar_replicate_sd(self):
        p = os.path.join(config.MARICOPA_OBS_DIR, "MARICOPA_YECORA_ROJO_REPS_SD.csv")
        if not os.path.exists(p):
            return
        lines = [l for l in open(p) if not l.startswith("#")]
        src = (" + replicate sd from MARICOPA_YECORA_ROJO_REPS_SD.csv "
               "(sheet Obs_summary_by_Reps, sample sd over replicate plots)")
        for r in csv.DictReader(lines):
            dk = "MARICOPA_%s" % r["TRNO"]
            if dk not in DECKS:
                continue
            for q, mk, sk in (("AGB_harvest", "CWAM_mean", "CWAM_sd"),
                              ("grain_yield", "HWAH_mean", "HWAH_sd"),
                              ("LAI_peak", "LAIX_mean", "LAIX_sd")):
                o = self.get(dk, q)
                m, sd = fl(r.get(mk, "")), fl(r.get(sk, ""))
                if o is None or isnan(sd) or isnan(m) or m <= 0 or not isnan(o["sd"]):
                    continue
                # relative sd of the replicate set applied to the registered value
                o["sd"] = o["value"] * (sd / m)
                o["basis"] = "sd"
                o["source"] = o["source"] + src
                o["note"] = (o["note"] + "; " if o["note"] else "") + \
                            "sd = value x (replicate sd / replicate mean), n = %s" % r["n_reps"]

    # -- Braunschweig: organ DM, crop N and the green-lamina area index series, all 16 decks
    def _bra_organs_and_N(self):
        for dk in [d for d in DECKS if DECKS[d]["site"] == "BRAUNSCHWEIG"]:
            p = bra_obs_path(dk)
            rr = rows_of(p)
            dates = sorted({r["date"] for r in rr if r["variable"] == "SUMDW"})
            if not dates:
                continue
            hd = dates[-1]
            at = {r["variable"]: r for r in rr if r["date"] == hd}

            def q(v, nm, qq, scale=1.0, note=""):
                r = at.get(v)
                if r and self.get(dk, qq) is None:
                    self.put(dk, qq, fl(r["mean"]) * scale, fl(r["sd"]) * scale, "sd",
                             "braunschweig/%s_obs.csv %s at %s "
                             "(Manderscheid et al. 2020 doi 10.18174/odjar.v6i0.16397)"
                             % (dk, v, hd), note=note, date=hd)
            q("STDW", "stem", "W_stem")
            q("LDeDW", "leaf", "W_leaf_tot", note="at maturity the observed lamina is dead "
                                                  "lamina only (no LGrDW row at this date)")
            q("EADW", "spike", "spike")
            q("CHDW", "chaff", "chaff")
            q("SUMDW", "AGB", "AGB_harvest")
            q("N plant", "cropN", "crop_N_harvest",
              note="model side is a DERIVED organ-N sum")
            q("GRDW", "grain", "grain_yield")
            q("TGW", "TKW", "TKW")
            q("GRNO", "kernels", "kernel_number")
            q("CP conc", "protein", "grain_protein_pct",
              note="obs = the file's own CP conc column")
            q("N% GR", "grainN", "grain_N_pct")
            q("N% ST", "stemN", "stem_N_pct")
            q("N% LDe", "leafN", "leaf_N_pct",
              note="at maturity the observed lamina is dead lamina only (N% LDe)")
            q("WSC", "wsc", "stem_WSC_mg_g")
            # the dated pre-anthesis SUMDW sample (the one immediately before harvest)
            dsum = sorted({r["date"] for r in rr if r["variable"] == "SUMDW"})
            if len(dsum) >= 2 and self.get(dk, "AGB_anthesis") is None:
                prev = [r for r in rr if r["date"] == dsum[-2] and r["variable"] == "SUMDW"]
                if prev:
                    self.put(dk, "AGB_anthesis", fl(prev[0]["mean"]), fl(prev[0]["sd"]), "sd",
                             "braunschweig/%s_obs.csv SUMDW at %s "
                             "(Manderscheid et al. 2020 doi 10.18174/odjar.v6i0.16397)"
                             % (dk, dsum[-2]), date=dsum[-2])
            # LGrAI = green lamina area index, the like-for-like basis for the model's LAI
            lai = sorted((r["date"], fl(r["mean"]), fl(r["sd"]))
                         for r in rr if r["variable"] == "LGrAI")
            if lai and self.get(dk, "LAI_peak") is None:
                top = max(lai, key=lambda t: t[1])
                self.put(dk, "LAI_peak", top[1], top[2], "sd",
                         "braunschweig/%s_obs.csv peak LGrAI, at %s "
                         "(Manderscheid et al. 2020 doi 10.18174/odjar.v6i0.16397)"
                         % (dk, top[0]),
                         note="peak of the dated green-lamina area index series",
                         date=top[0])
            self.lai_ts[dk] = (lai, "Braunschweig LGrAI (green lamina area index)")

    # -- Maricopa: scalars from the measured evidence table of the archive, no dispersion
    def _mar_scalars(self):
        p = os.path.join(config.MARICOPA_OBS_DIR, "MARICOPA_YECORA_ROJO_MEASURED.csv")
        lines = [l for l in open(p) if not l.startswith("#")]
        rr = list(csv.DictReader(lines))
        srcs = ("MARICOPA_YECORA_ROJO_MEASURED.csv "
                "(Kimball et al. 2017 ODJAR 3:29-38 doi 10.18174/odjar.v3i1.15826); no "
                "dispersion in any column")
        for r in rr:
            dk = "MARICOPA_%s" % r["TRNO"]
            if dk not in DECKS:
                continue
            cw = fl(r["CWAM_kg_ha"]) / 10.0
            if not isnan(cw) and self.get(dk, "AGB_harvest") is None:
                self.put(dk, "AGB_harvest", cw, NAN, "pct", srcs + " [CWAM_kg_ha / 10]")
            tk = fl(r["TKW_measured_g"])
            kn = fl(r["grain_number_m2"])
            if not isnan(tk) and self.get(dk, "TKW") is None:
                self.put(dk, "TKW", tk, NAN, "pct", srcs + " [TKW_measured_g = HWUM x 1000]")
            if not isnan(kn) and self.get(dk, "kernel_number") is None:
                self.put(dk, "kernel_number", kn, NAN, "pct",
                         srcs + " [grain_number_m2 = H#UM x P#AM]")
            ea = fl(r["ears_m2_PnAM"])
            if not isnan(ea) and self.get(dk, "ears_m2") is None:
                self.put(dk, "ears_m2", ea, NAN, "pct",
                         srcs + " [ears_m2_PnAM = P#AM, panicle (ear) number at maturity]")
            if not isnan(tk) and not isnan(kn) and self.get(dk, "grain_yield") is None:
                self.put(dk, "grain_yield", tk * kn / 1000.0, NAN, "pct",
                         srcs + " [DERIVED = TKW_measured_g x grain_number_m2 / 1000]",
                         note="DERIVED obs")

    # -- SWC observation series, by site
    def _collect_swc(self):
        for dk in [d for d in DECKS if DECKS[d]["site"] == "WDF25"]:
            cv, co2, wat = dk.split("_")[1:4]
            rows = W25.obs_for(cv, co2, wat)
            for var, lab in (("SWC 20cm", "SWC20"), ("SWC 55cm", "SWC55")):
                self._swc_obs[(dk, lab)] = sorted(
                    (r["date"], fl(r["value"]) * 0.01, fl(r["std"], 0) * 0.01)
                    for r in rows if r["variable"] == var)
        for dk in [d for d in DECKS if DECKS[d]["site"] == "WDF23"]:
            cv, co2 = dk.split("_")[1:3]
            rows = W23.obs_for(cv, co2)
            for var, lab in (("SWC 20cm", "SWC20"), ("SWC 55cm", "SWC55")):
                self._swc_obs[(dk, lab)] = sorted(
                    (r["date"], fl(r["value"]) * 0.01, fl(r["std"], 0) * 0.01)
                    for r in rows if r["variable"] == var)
        for dk in [d for d in DECKS if DECKS[d]["site"] == "BRAUNSCHWEIG"]:
            rr = rows_of(bra_obs_path(dk))
            for var, lab in (("d0_20", "d0_20"), ("d20_40", "d20_40")):
                self._swc_obs[(dk, lab)] = sorted(
                    (r["date"], fl(r["mean"]) * 0.01, fl(r["sd"], 0) * 0.01)
                    for r in rr if r["variable"] == var)

    def swc_obs(self, deck, label):
        return self._swc_obs.get((deck, label), [])


# ================================================================= soil water and phenology
def layer_bounds(pkg, deck):
    """The deck's own soil-layer boundary depths, read from its driver_swtr.txt."""
    p = os.path.join(deck_path(pkg, deck), "input", "driver_swtr.txt")
    if not os.path.isfile(p):
        return None
    lines = open(p).read().splitlines()
    vals = []
    for ln in lines[15:26]:
        tok = ln.split("!")[0].strip()
        v = fl(tok)
        if not isnan(v):
            vals.append(abs(v))
    return vals if len(vals) == 11 else None


def swc_span_mean(mrows, bnds, day, z0, z1):
    """Thickness-weighted mean v/v over [z0, z1] for one model day (used where the observation
    is itself a depth-span integral)."""
    r = mrows.get(day)
    if r is None or not bnds:
        return None
    num = den = 0.0
    for i in range(10):
        a, b = bnds[i], bnds[i + 1]
        lo, hi = max(a, z0), min(b, z1)
        if hi <= lo:
            continue
        th = b - a
        theta = fl(r["layer%d_mm" % (i + 1)]) / (th * 1000.0)
        num += theta * (hi - lo)
        den += (hi - lo)
    return num / den if den > 0 else None


class SwcBasis(object):
    """Per-site soil-water comparison basis: point sensors at 0.20 / 0.55 m at WheatDryFACE
    (interpolated between layer mid-points), depth-span integrals 0-20 / 20-40 cm at
    Braunschweig; the model side is the time-weighted daily-mean layer file in both cases."""

    def __init__(self, pkg, deck):
        site = DECKS[deck]["site"]
        d = deck_dir(pkg, deck)
        self.bnds = layer_bounds(pkg, deck)
        self.site = site
        if site in ("WDF25", "WDF23"):
            self.file = "daily_swc_layers_mean.csv"
            self.labels = [("SWC20", 0.20, 0.20), ("SWC55", 0.55, 0.55)]
            self.kind = "point"
        elif site == "BRAUNSCHWEIG":
            self.file = "daily_swc_layers_mean.csv"
            self.labels = [("d0_20", 0.0, 0.20), ("d20_40", 0.20, 0.40)]
            self.kind = "span"
        else:
            self.file = None
            self.labels = []
            self.kind = "none"
        self.rows = ({int(fl(r["day"])): r for r in rows_of(os.path.join(d, self.file))}
                     if self.file else {})

    def model_at(self, day, z0, z1):
        if self.kind == "point":
            r = self.rows.get(day)
            if r is None or not self.bnds:
                return None
            mid = [(self.bnds[i] + self.bnds[i + 1]) / 2.0 for i in range(10)]
            v = [fl(r["layer%d_mm" % (i + 1)]) / ((self.bnds[i + 1] - self.bnds[i]) * 1000.0)
                 for i in range(10)]
            z = z0
            if z <= mid[0]:
                return v[0]
            if z >= mid[-1]:
                return v[-1]
            for i in range(9):
                if mid[i] <= z <= mid[i + 1]:
                    f = (z - mid[i]) / (mid[i + 1] - mid[i])
                    return v[i] + f * (v[i + 1] - v[i])
            return v[-1]
        if self.kind == "span":
            return swc_span_mean(self.rows, self.bnds, day, z0, z1)
        return None


def swc_pairs(pkg, obs, deck):
    """[(label, [(date, obs, obs_sd, model)], bias, n)] for one deck."""
    m = MODELS[deck]
    sb = SwcBasis(pkg, deck)
    out = []
    for lab, z0, z1 in sb.labels:
        series = []
        for date, ov, osd in obs.swc_obs(deck, lab):
            day = m.date2day.get(date)
            if day is None:
                continue
            mv = sb.model_at(day, z0, z1)
            if mv is None:
                continue
            series.append((date, ov, osd, mv))
        bias = (sum(s[3] - s[1] for s in series) / len(series)) if series else NAN
        out.append((lab, series, bias, len(series)))
    return out, sb


def phenology_rows(pkg, deck):
    """[(stage, model_date, obs_date, diff_d, scored, source)] by calendar date."""
    m = MODELS[deck]
    site = DECKS[deck]["site"]
    out = []
    obs_map, src = {}, ""
    if site == "WDF25":
        cv, co2, wat = deck.split("_")[1:4]
        sd = {s["stage"]: s["date"] for s in W25.stage_dates(cv, co2, wat)}
        obs_map = {k: sd.get(v) for k, v in WDF25_STAGE_MAP}
        src = "WDF25/S3_%s_%s_%s.csv 'Phenological stage' dated column" % (cv, co2, wat)
    elif site == "WDF23":
        cv, co2 = deck.split("_")[1:3]
        sd = {s["stage"]: s["date"] for s in W23.stage_dates(cv, co2)}
        obs_map = {k: sd.get(v) for k, v in WDF23_STAGE_MAP}
        src = "WDF23/WDF23_%s_%s.csv 'Phenological stage' dated column" % (cv, co2)
    elif site == "BRAUNSCHWEIG":
        rr = rows_of(bra_obs_path(deck))
        dd = {r["variable"]: r["date"] for r in rr if r["unit"] == "date"}
        obs_map = {k: dd.get(v) for k, v in BRA_STAGE_MAP}
        obs_map["3_leaves"] = None      # see SCORED_STAGES
        src = "braunschweig/%s_obs.csv dated stage rows" % deck
    else:
        p = os.path.join(deck_path(pkg, deck), "input", "phenology_cultivar.txt")
        kv, sow = {}, None
        if os.path.isfile(p):
            toks = [l.split() for l in open(p) if l.strip()]
            d = {}
            for t in toks:
                if len(t) >= 2:
                    d[t[1].strip()] = t[0].strip()
            try:
                sow = dt.date(int(fl(d["Sowing_Year"])), int(fl(d["Sowing_Month"])),
                              int(fl(d["Sowing_Day"])))
            except Exception:
                sow = None
            kv = d
        for mk, ok in DAP_STAGE_MAP:
            v = fl(kv.get(ok, "nan"))
            # DAP -> calendar date = sowing + DAP days; a DAP <= 0 is a placeholder, not a measurement
            obs_map[mk] = ((sow + dt.timedelta(days=int(v))).isoformat()
                           if sow and not isnan(v) and v > 0 else None)
        src = ("<archive>/%s/input/phenology_cultivar.txt Sowing_Y/M/D + Obs_*_DAP "
               "(DAP from sowing, converted to a calendar date)" % deck)
    for mk, _ in DAP_STAGE_MAP:
        md = m.stages.get(mk)
        od = obs_map.get(mk)
        scored = mk in SCORED_STAGES.get(site, set())
        diff = NAN
        if md and od:
            diff = (dt.date(*[int(x) for x in md.split("-")])
                    - dt.date(*[int(x) for x in od.split("-")])).days
        out.append((mk, md or "", od or "", diff, scored, src))
    return out


# ================================================================= treatment contrasts
# Quantities carried through every contrast.
CQ = ["AGB_harvest", "grain_yield", "TKW", "kernel_number", "HI", "LAI_peak", "height_max",
      "grain_protein_pct", "grain_N_pct", "leaf_N_pct", "stem_N_pct", "chaff",
      "grain_Fe", "grain_Zn", "grain_Mn", "grain_Cu", "grain_Mg", "grain_Ca"]

# (label, kind, site, numerator deck, denominator deck).  WDF25 carries the full factorial:
# eCO2/aCO2 inside each water x cultivar combination, WL/WW inside each CO2 x cultivar combination,
# Zahir/Ruta inside each CO2 x water combination.
CONTRASTS = []


def _c(lab, kind, site, num, den):
    CONTRASTS.append((lab, kind, site, num, den))


for cv in ("Ruta", "Zahir"):
    for wat in ("WW", "WL"):
        _c("WDF25 eCO2/aCO2 %s %s" % (cv, wat), "eCO2", "WDF25",
           "WDF25_%s_eCO2_%s" % (cv, wat), "WDF25_%s_aCO2_%s" % (cv, wat))
for cv in ("Ruta", "Zahir"):
    for co2 in ("aCO2", "eCO2"):
        _c("WDF25 WL/WW %s %s" % (cv, co2), "water", "WDF25",
           "WDF25_%s_%s_WL" % (cv, co2), "WDF25_%s_%s_WW" % (cv, co2))
for co2 in ("aCO2", "eCO2"):
    for wat in ("WW", "WL"):
        _c("WDF25 Zahir/Ruta %s %s" % (co2, wat), "cultivar", "WDF25",
           "WDF25_Zahir_%s_%s" % (co2, wat), "WDF25_Ruta_%s_%s" % (co2, wat))
for cv in ("Ruta", "Zahir"):
    _c("WDF23 eCO2/aCO2 %s" % cv, "eCO2", "WDF23",
       "WDF23_%s_eCO2" % cv, "WDF23_%s_aCO2" % cv)
for co2 in ("aCO2", "eCO2"):
    _c("WDF23 Zahir/Ruta %s" % co2, "cultivar", "WDF23",
       "WDF23_Zahir_%s" % co2, "WDF23_Ruta_%s" % co2)
for nl, ea, aa in (("N1", "104", "101"), ("N2", "105", "102"),
                   ("N3", "106", "103"), ("N4", "110", "107")):
    _c("Braunschweig eCO2/aCO2 %s" % nl, "eCO2", "BRAUNSCHWEIG",
       "BRAUNSCHWEIG_%s" % ea, "BRAUNSCHWEIG_%s" % aa)
for lab, num, den in (("N1/N3", "101", "103"), ("N2/N3", "102", "103"),
                      ("N4/N2", "107", "102")):
    _c("Braunschweig %s aCO2" % lab, "N", "BRAUNSCHWEIG",
       "BRAUNSCHWEIG_%s" % num, "BRAUNSCHWEIG_%s" % den)
for lab, num, den in (("Dry 93/94", "907", "905"), ("Wet 93/94", "908", "906"),
                      ("Wet 95/96", "912", "910"), ("Wet 96/97", "916", "914"),
                      ("Low-N 96/97", "915", "913")):
    _c("Maricopa eCO2/aCO2 %s" % lab, "eCO2", "MARICOPA",
       "MARICOPA_%s" % num, "MARICOPA_%s" % den)
for lab, num, den in (("aCO2 93/94", "905", "906"), ("eCO2 93/94", "907", "908")):
    _c("Maricopa Dry/Wet %s" % lab, "water", "MARICOPA",
       "MARICOPA_%s" % num, "MARICOPA_%s" % den)
for lab, num, den in (("aCO2 96/97", "913", "914"), ("eCO2 96/97", "915", "916")):
    _c("Maricopa Low-N/High-N %s" % lab, "N", "MARICOPA",
       "MARICOPA_%s" % num, "MARICOPA_%s" % den)


def sign_class(rm, ro, ro_sd):
    """Sign agreement of a simulated (rm) and observed (ro) treatment ratio.  Returns (class, resolved_flag):
    'agree' = same sign; 'A' = opposite sign with both effects resolved; 'B' = opposite sign, unresolved;
    'C' = no observation."""
    if isnan(ro):
        return "C", False
    resolved = ((not isnan(ro_sd) and ro_sd > 0 and abs(ro - 1.0) > ro_sd)
                or (isnan(ro_sd) and abs(ro - 1.0) >= SIGN_RESOLVE_PCT))
    same = (rm - 1.0) * (ro - 1.0) >= 0.0
    if same:
        return "agree", resolved
    if abs(rm - 1.0) >= SIGN_FLOOR and resolved:
        return "A", resolved
    return "B", resolved


def contrast_rows(obs):
    out = []
    for lab, kind, site, num, den in CONTRASTS:
        for q in CQ:
            mn = MODELS[num].v.get(q) if num in MODELS else NAN
            md = MODELS[den].v.get(q) if den in MODELS else NAN
            if q.startswith("grain_") and q[6:] in ("Fe", "Zn", "Mn", "Cu", "Mg", "Ca"):
                mn = GRAIN_MIN[num].get(q[6:], NAN)
                md = GRAIN_MIN[den].get(q[6:], NAN)
            rm = (mn / md) if (not isnan(mn) and not isnan(md) and md != 0) else NAN
            on, od = obs.get(num, q), obs.get(den, q)
            ro = ro_sd = NAN
            if on and od and not isnan(on["value"]) and not isnan(od["value"]) and od["value"]:
                ro = on["value"] / od["value"]
                if not isnan(on["sd"]) and not isnan(od["sd"]) and on["value"] and od["value"]:
                    ro_sd = abs(ro) * math.sqrt((on["sd"] / on["value"]) ** 2
                                                + (od["sd"] / od["value"]) ** 2)
            cls, resolved = sign_class(rm, ro, ro_sd)
            if isnan(rm):
                continue
            # A quantity whose simulated value is prescribed cannot carry a response sign (its
            # ratio is 1 by construction); such rows are marked class C here.
            pres_reason = ""
            if q == "grain_Ca":
                cls, resolved = "C", False
                pres_reason = ("model grain Ca is a CONSTANT (xylem-sap concentration input): "
                               "the ratio is 1.000 by construction")
            elif q == "TKW" and site == "WDF25":
                cls, resolved = "C", False
                pres_reason = ("WDF25 TKW is bounded by the potential kernel mass input: "
                               "the ratio is 1.000 by construction")
            out.append(dict(contrast=lab, kind=kind, site=site, num=num, den=den, quantity=q,
                            R_model=rm, R_obs=ro, R_obs_sd=ro_sd, cls=cls, resolved=resolved,
                            prescribed=pres_reason))
    return out


# ================================================================= per-deck extra reads
MINERALS = ["Fe", "Zn", "Mn", "Cu", "Mg", "Ca"]


def grain_minerals(pkg, deck):
    """Model grain mineral concentration (mg kg-1) on the last modeled day.
    A negative grain_mg_kg is the module's own not-modelled sentinel (-999)."""
    rr = rows_of(os.path.join(deck_dir(pkg, deck), "daily_minerals.csv"))
    if not rr:
        return {}
    last = max(int(fl(r["day"])) for r in rr)
    out = {}
    for r in rr:
        if int(fl(r["day"])) == last:
            v = fl(r["grain_mg_kg"])
            out[r["mineral"].strip()] = (NAN if (not isnan(v) and v < 0) else v)
    return out


def grain_fill_law(pkg, deck):
    """(filling_end_bound, TKW_potential_input_g, kernel_mass_ceiling_flag) from grain_filling_law.csv."""
    p = os.path.join(deck_dir(pkg, deck), "grain_filling_law.csv")
    bound, tkwp, r118 = "", NAN, ""
    for ln in open(p) if os.path.isfile(p) else []:
        t = ln.split(",")
        if len(t) < 2:
            continue
        k = t[0].strip()
        if k == "filling_end_bound":
            bound = t[1].strip()
        elif k == "TKW_potential_input_g":
            tkwp = fl(t[1])
        elif k == "rec118_kernel_mass_ceiling":
            r118 = t[1].strip()
    return bound, tkwp, r118


def gnx_summary(pkg, deck, key):
    """One value from the key/value summary block of grain_number_chain.csv."""
    p = os.path.join(deck_dir(pkg, deck), "grain_number_chain.csv")
    if not os.path.isfile(p):
        return NAN
    for ln in open(p):
        t = [x.strip() for x in ln.split(",")]
        if len(t) >= 2 and t[0] == key:
            return fl(t[1])
    return NAN


def ear_number(pkg, deck):
    """Simulated ear count (m-2) = the `ears_m2` summary row of grain_number_chain.csv (the
    spikelet-weighted ear count; the observed counterpart, P#AM, exists at Maricopa only)."""
    return gnx_summary(pkg, deck, "ears_m2")


def agb_at_dated_sample(pkg, deck, obs):
    """Simulated AGB on the calendar date of the observed dated pre-anthesis sample
    (Braunschweig only: the SUMDW sample preceding the final harvest); NaN elsewhere."""
    o = obs.get(deck, "AGB_anthesis")
    if o is None or not o.get("date"):
        if DECKS[deck]["site"] != "BRAUNSCHWEIG":
            return NAN
        rr = rows_of(bra_obs_path(deck))
        ds = sorted({r["date"] for r in rr if r["variable"] == "SUMDW"})
        if len(ds) < 2:
            return NAN
        date = ds[-2]
    else:
        date = o["date"]
    m = MODELS[deck]
    day = m.date2day.get(date)
    if day is None:
        return NAN
    org = m.organ_at_day(day)
    return org["AGB_harvest"] if org else NAN


def crop_N_derived(pkg, deck):
    """Derived above-ground N = N_leaf + N_dead_leaf + N_stem + N_grain + N_chaff on the last day."""
    rr = rows_of(os.path.join(deck_dir(pkg, deck), "daily_nitrogen.csv"))
    if not rr:
        return NAN
    r = rr[-1]
    ks = ["N_leaf_gN_m2", "N_dead_leaf_gN_m2", "N_stem_gN_m2", "N_grain_gN_m2", "N_chaff_gN_m2"]
    if not all(k in r for k in ks):
        return NAN
    return sum(fl(r[k], 0.0) for k in ks)


# ================================================================= run-time state (filled by paperfigs.init)
PKG = None
MODELS = {}
OBS_REG = None
PHEN = {}
SWC = {}
CAPS = {}
EARS = {}
GRAIN_MIN = {}
CROPN = {}
AGB_ANTH = {}
MAR_EARS = {}
ORGAN_N_OBS = {}
CROWS = []
