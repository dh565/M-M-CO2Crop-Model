#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Reader for the WheatDryFACE treatment files (seasons 2024/25 and 2022/23).

Each treatment file has a two-row header: row 1 carries group labels (Night/Day for the
radiation blocks, organ names Grain/Stem/Spike/Leaf/Chaff (and Straw in 2022/23) for the
harvest block, blank for the daily and campaign columns), row 2 the variable names, each
value column followed by its "Std" column.  Genotype, CO2 level and water treatment are
taken from the file name.

The reader returns tidy rows: date, dap, genotype, co2, water, variable, value, std, period,
organ, units, category, stage, source_file.  One derived row is appended for 2024/25, the
harvest index HI = Biomass[Grain] / (Biomass[Stem] + Biomass[Leaf] + Biomass[Spike]), with
the standard deviation propagated in quadrature (Spike = Grain + Chaff in these files, so
Stem + Leaf + Spike counts every above-ground organ once).

  W25.obs_for(genotype, co2, water) / W25.stage_dates(...)   season 2024/25
  W23.obs_for(genotype, co2)        / W23.stage_dates(...)   season 2022/23 (no water treatment)
"""
import csv
import datetime
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

# variable -> (units, category)
UNITS = {
    "An": ("umol/m2/s", "gasexchange_campaign"), "gs": ("mol/m2/s", "gasexchange_campaign"),
    "Tr": ("mmol/m2/s", "gasexchange_campaign"),
    "Tleaf_flir": ("C", "campaign"), "DTleaf_air_flir": ("C", "campaign"),
    "LAI": ("m2/m2", "campaign"), "Height": ("cm", "campaign"),
    "SWC 20cm": ("v/v", "daily_sensor"), "SWC 55cm": ("v/v", "daily_sensor"),
    "Tsoil 20cm": ("C", "daily_sensor"), "Tsoil55 cm": ("C", "daily_sensor"),
    "SWC gravimetric_5cm": ("g/g", "field_sample"), "SWC gravimetric_20cm": ("g/g", "field_sample"),
    "SWC gravimetric_55cm": ("g/g", "field_sample"), "SWC volumetric_5cm": ("v/v", "field_sample"),
    "SWC volumetric_20cm": ("v/v", "field_sample"), "SWC volumetric_55cm": ("v/v", "field_sample"),
    "albedo": ("unitless", "radiation_daynight"), "Rn": ("W/m2", "radiation_daynight"),
    "LWin": ("W/m2", "radiation_daynight"), "LWout": ("W/m2", "radiation_daynight"),
    "SWin": ("W/m2", "radiation_daynight"), "SWout": ("W/m2", "radiation_daynight"),
    "Soil heat flux": ("W/m2", "radiation_daynight"),
    "Biomass": ("t/ha", "harvest"), "[C]": ("%DM", "harvest"), "[N]": ("%DM", "harvest"),
    "Fe": ("mg/kg", "harvest"), "Mg": ("mg/kg", "harvest"), "Ca": ("mg/kg", "harvest"),
    "Cu": ("mg/kg", "harvest"), "Mn": ("mg/kg", "harvest"), "Zn": ("mg/kg", "harvest"),
    "Protein%": ("%", "harvest"),
    "P+I": ("mm", "daily_forcing"),
    "Tmax": ("C", "daily_forcing"), "Tmin": ("C", "daily_forcing"), "Tmean": ("C", "daily_forcing"),
    "TKW": ("g/1000 grains", "harvest"),
    "HI": ("-", "harvest_derived"),
}
ALIASES = {"P+I_(mm)": ("P+I", "mm"), "P+I_(mm)_WL": ("P+I", "mm")}
ORGANS_25 = {"Grain", "Stem", "Spike (head)", "Leaf", "Chaff"}
ORGANS_23 = ORGANS_25 | {"Straw"}

HI_NUM_ORGAN = "Grain"
HI_AGB_ORGANS = ("Stem", "Leaf", "Spike")
HI_N_REPLICATES = 3   # ring replicates per CO2 level


def derived_hi(recs):
    """(HI, HI_std, GY, GY_std, AGB, AGB_std, date, n) from the tidy rows of one treatment, or None."""
    b = {}
    for r in recs:
        if r["variable"] == "Biomass" and r["organ"]:
            b[r["organ"]] = (r["value"], r["std"], r["date"])
    need = (HI_NUM_ORGAN,) + HI_AGB_ORGANS
    if any(o not in b for o in need):
        return None
    gy, gys, gdate = b[HI_NUM_ORGAN]
    agb = sum(b[o][0] for o in HI_AGB_ORGANS)
    agbs = math.sqrt(sum((b[o][1] or 0.0) ** 2 for o in HI_AGB_ORGANS))
    if not gy or not agb:
        return None
    hi = gy / agb
    his = hi * math.sqrt(((gys or 0.0) / gy) ** 2 + (agbs / agb) ** 2)
    return (hi, his, gy, gys, agb, agbs, gdate, HI_N_REPLICATES)


def _hi_record(recs, tr, fname, sowing):
    d = derived_hi(recs)
    if d is None:
        return None
    hi, his, _gy, _gys, _agb, _agbs, gdate, _n = d
    units, cat = UNITS["HI"]
    return dict(date=gdate, dap=(_parse_date(gdate) - sowing).days + 1,
                genotype=tr["genotype"], co2=tr["co2"], water=tr["water"],
                variable="HI", value=hi, std=his, period=None, organ=None,
                units=units, category=cat, stage="Derived (harvest)",
                source_file=os.path.basename(fname))


def _parse_date(s):
    s = s.strip()
    if not s:
        return None
    for fmt in ("%m/%d/%Y", "%Y-%m-%d", "%d/%m/%Y"):
        try:
            return datetime.datetime.strptime(s, fmt).date()
        except ValueError:
            pass
    return None


def _fnum(s):
    s = (s or "").strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _column_map(row_groups, row_vars):
    """[(idx, group, variable, std_idx)] pairing each value column with its Std column."""
    n = len(row_vars)
    groups, cur = [], ""
    for x in row_groups:
        x = (x or "").strip()
        if x:
            cur = x
        groups.append(cur)
    cols = []
    i = 2   # skip Date(0), Phenological stage(1)
    while i < n:
        var = (row_vars[i] or "").strip()
        if var == "" or var.lower() == "std":
            i += 1
            continue
        std_idx = i + 1 if (i + 1 < n and (row_vars[i + 1] or "").strip().lower() == "std") else None
        cols.append((i, groups[i] if i < len(groups) else "", var, std_idx))
        i = (std_idx + 1) if std_idx else (i + 1)
    return cols


def load_file(fname, tr, sowing, organs):
    """Tidy rows of one treatment file."""
    rows = list(csv.reader(open(fname, encoding="utf-8-sig")))
    cols = _column_map(rows[0], rows[1])
    out = []
    for r in rows[2:]:
        if not r or not r[0].strip():
            continue
        d = _parse_date(r[0])
        if d is None:
            continue
        dap = (d - sowing).days + 1
        stage = r[1].strip() if len(r) > 1 else ""
        for idx, group, var, std_idx in cols:
            if idx >= len(r):
                continue
            val = _fnum(r[idx])
            if val is None:
                continue
            std = _fnum(r[std_idx]) if (std_idx is not None and std_idx < len(r)) else None
            if var in ALIASES:
                var, _al_units = ALIASES[var]
                units_dict, cat = UNITS.get(var, ("", "unknown"))
                out.append(dict(date=d.isoformat(), dap=dap, genotype=tr["genotype"], co2=tr["co2"],
                                water=tr["water"], variable=var, value=val, std=std, period=None,
                                organ=None, units=_al_units or units_dict, category=cat, stage=stage,
                                source_file=os.path.basename(fname)))
                continue
            if var.endswith(")") and "(" in var:
                base = var[:var.index("(")].strip()
                inline_units = var[var.index("(") + 1:-1].strip()
            else:
                base = var
                inline_units = None
            var = base
            units_dict, cat = UNITS.get(base, ("", "unknown"))
            units = inline_units or units_dict
            period = group if group in ("Night", "Day") else None
            organ = group if group in organs else None
            out.append(dict(date=d.isoformat(), dap=dap, genotype=tr["genotype"], co2=tr["co2"],
                            water=tr["water"], variable=var, value=val, std=std, period=period,
                            organ=("Spike" if organ == "Spike (head)" else organ),
                            units=units, category=cat, stage=stage, source_file=os.path.basename(fname)))
    _hi = _hi_record(out, tr, fname, sowing)
    if _hi is not None:
        out.append(_hi)
    return out


def _stage_dates(fname, sowing):
    """First calendar date of every 'Phenological stage' label, including rows without data."""
    rows = list(csv.reader(open(fname, encoding="utf-8-sig")))
    seen = {}
    for r in rows[2:]:
        if not r or not r[0].strip():
            continue
        d = _parse_date(r[0])
        st = r[1].strip() if len(r) > 1 else ""
        if d is None or not st:
            continue
        if st not in seen:
            has_data = any((c or "").strip() for c in r[2:])
            seen[st] = dict(stage=st, date=d.isoformat(), dap=(d - sowing).days + 1, has_data=has_data)
    return sorted(seen.values(), key=lambda x: x["date"])


class _W25(object):
    """WheatDryFACE 2024/25: files S3_<genotype>_<CO2>_<water>.csv; day 1 = 2024-11-18 (field preparation)."""
    SOWING = datetime.date(2024, 11, 18)

    def __init__(self, datadir=None):
        self.datadir = datadir or config.WDF25_DIR

    def path(self, genotype, co2, water):
        return os.path.join(self.datadir, "S3_%s_%s_%s.csv" % (genotype, co2, water))

    def obs_for(self, genotype, co2, water):
        return load_file(self.path(genotype, co2, water), dict(genotype=genotype, co2=co2, water=water),
                         self.SOWING, ORGANS_25)

    def stage_dates(self, genotype, co2, water):
        return _stage_dates(self.path(genotype, co2, water), self.SOWING)


class _W23(object):
    """WheatDryFACE 2022/23: files WDF23_<genotype>_<CO2>.csv; day 1 = 2022-12-07; no water treatment
    (the water label "WW" is carried so that the rows have the same keys as 2024/25)."""
    SOWING = datetime.date(2022, 12, 7)

    def __init__(self, datadir=None):
        self.datadir = datadir or config.WDF23_DIR

    def path(self, genotype, co2):
        return os.path.join(self.datadir, "WDF23_%s_%s.csv" % (genotype, co2))

    def obs_for(self, genotype, co2, water=None):
        return load_file(self.path(genotype, co2), dict(genotype=genotype, co2=co2, water="WW"),
                         self.SOWING, ORGANS_23)

    def stage_dates(self, genotype, co2, water=None):
        return _stage_dates(self.path(genotype, co2), self.SOWING)


W25, W23 = _W25(), _W23()
