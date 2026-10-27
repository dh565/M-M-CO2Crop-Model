#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Maricopa FACE hourly canopy temperature and energy-balance observations -> one tidy CSV per modelled treatment.

Sources: 'Weather Canopy & Soil Temperatures Energy Balance 1994 hourly.ods' (1993/94), '... CO2 1996 15-min.ods'
(1995/96), '... CO2 O3 1997 15-min.ods' (1996/97); Kimball et al. (1999) Water Resour. Res. 35:1179 and Kimball et
al. (2017).  Plot -> TRNO mapping: 1993/94 CW -> 906 (wet aCO2), FW -> 908 (wet eCO2); 1995/96 BH -> 910, FH -> 912;
1996/97 BH -> 914, FH -> 916, BL -> 913, FL -> 915.  -99/-999 = missing.  15-min rows are averaged to the hour
(hour ending); replicate plots are averaged and the sd across replicate plots is kept.

The shipped observations/maricopa/maricopa_eb_hourly.csv is the output of this script.  Run it to rebuild that
file from the archive:  python3 lib/extract_maricopa_eb.py [output.csv]
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import maricopa_archive as MA

SPEC = {
 "eb94": dict(season="1993/94", hdr=8, data0=10, tair=4, has_min=False,
              tcan={"CW": [15, 16, 17, 18], "FW": [19, 20, 21, 22]},
              rn={"CW": [56, 57, 58, 59], "FW": [60, 61, 62, 63]},
              h={"CW": [67, 68, 69, 70], "FW": [71, 72, 73, 74]},
              g={"CW": [78, 79], "FW": [80, 81]},
              le={"CW": [84, 85, 86, 87], "FW": [88, 89, 90, 91]},
              trno={"CW": 906, "FW": 908}),
 "eb96": dict(season="1995/96", hdr=8, data0=10, tair=5, has_min=True,
              tcan={"BH": [31, 32], "FH": [33, 34], "BL": [35, 36], "FL": [37, 38]},
              rn={"BH": [81, 82], "FH": [83, 84], "BL": [85, 86], "FL": [87, 88]},
              h={"BH": [92, 93], "FH": [94, 95], "BL": [96, 97], "FL": [98, 99]},
              g={"BH": [102], "FH": [103], "BL": [104], "FL": [105]},
              le={"BH": [109, 110], "FH": [111, 112], "BL": [113, 114], "FL": [115, 116]},
              trno={"BH": 910, "FH": 912}),
 "eb97": dict(season="1996/97", hdr=8, data0=10, tair=5, has_min=True,
              tcan={"BH": [31, 32], "FH": [33, 34], "BL": [35, 36], "FL": [37, 38]},
              rn={"BH": [73, 74], "FH": [75, 76], "BL": [77, 78], "FL": [79, 80]},
              h={"BH": [82, 83], "FH": [84, 85], "BL": [86, 87], "FL": [88, 89]},
              g={"BH": [91], "FH": [92], "BL": [93], "FL": [94]},
              le={"BH": [96, 97], "FH": [98, 99], "BL": [100, 101], "FL": [102, 103]},
              trno={"BH": 914, "FH": 916, "BL": 913, "FL": 915}),
}


def num(s):
    v = pd.to_numeric(s, errors="coerce")
    return v.where(~v.isin([-99, -999, -99.0, -999.0]))


def build():
    out = []
    for tag, sp in SPEC.items():
        df = MA.eb_workbook(tag).iloc[sp["data0"]:].reset_index(drop=True)
        year = num(df[0]); doy = num(df[1]); hour = num(df[2])
        if sp["has_min"]:
            minute = num(df[3]); hour = hour + (minute >= 59).astype(float) * 0
            # HOUR h with MINUTE 15, 30, 45 and HOUR h+1 MINUTE 0 belong to the hour ending at h+1
            hour_end = np.where(minute == 0, hour, hour + 1)
        else:
            hour_end = hour
        base = pd.DataFrame({"year": year, "doy": doy, "hour": hour_end, "tair": num(df[sp["tair"]])})
        for grp, trno in sp["trno"].items():
            d = base.copy(); d["trno"] = trno; d["season"] = sp["season"]
            for var in ("tcan", "rn", "h", "g", "le"):
                cols = sp[var][grp]
                block = pd.concat([num(df[c]) for c in cols], axis=1)
                d[var] = block.mean(axis=1, skipna=True)
                d[var + "_sd"] = block.std(axis=1, ddof=1) if len(cols) > 1 else np.nan
                d[var + "_n"] = block.notna().sum(axis=1)
            out.append(d)
    tidy = pd.concat(out, ignore_index=True).dropna(subset=["year", "doy", "hour"])
    tidy = tidy[(tidy.hour >= 1) & (tidy.hour <= 24)]
    agg = {c: "mean" for c in tidy.columns if c not in ("year", "doy", "hour", "trno", "season")}
    hourly = tidy.groupby(["season", "trno", "year", "doy", "hour"], as_index=False).agg(agg)
    hourly["date"] = pd.to_datetime(hourly.year.astype(int).astype(str) + "-" + hourly.doy.astype(int).astype(str), format="%Y-%j")
    return hourly


if __name__ == "__main__":
    dest = sys.argv[1] if len(sys.argv) > 1 else os.path.join(config.out_subdir("observations_rebuilt"), "maricopa_eb_hourly.csv")
    hourly = build()
    hourly.to_csv(dest, index=False)
    print("wrote", dest)
    print(hourly.groupby(["season", "trno"]).agg(n=("tcan", "size"), tcan=("tcan", lambda x: x.notna().sum()), rn=("rn", lambda x: x.notna().sum()),
                                                    h=("h", lambda x: x.notna().sum()), le=("le", lambda x: x.notna().sum()), g=("g", lambda x: x.notna().sum()),
                                                    d0=("date", "min"), d1=("date", "max")))
