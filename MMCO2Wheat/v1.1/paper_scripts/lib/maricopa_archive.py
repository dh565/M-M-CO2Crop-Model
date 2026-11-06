#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Readers for the Maricopa FACE archive workbooks (Kimball et al., 2017, Open Data Journal for
Agricultural Research 3:29-38; data doi 10.7910/DVN/YRTQ3Z), read from config.MARICOPA_ARCHIVE_DIR.

The .ods workbooks take one to two minutes each to parse; the parsed frames are cached as
pickles under config.CACHE_DIR (inside OUT_DIR, never inside the archive).
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


def _cached(name, loader):
    os.makedirs(config.CACHE_DIR, exist_ok=True)
    cache = os.path.join(config.CACHE_DIR, name + ".pkl")
    if os.path.isfile(cache):
        return pd.read_pickle(cache)
    df = loader()
    df.to_pickle(cache)
    return df


def eb_workbook(tag):
    """Raw sheet (header=None) of the hourly / 15-min energy-balance workbook 'eb94', 'eb96' or 'eb97'."""
    path = config.require(config.MARICOPA_EB_ODS[tag], "Maricopa energy-balance workbook %s" % tag)
    return _cached(tag, lambda: pd.read_excel(path, engine="odf", header=None))


def gas_exchange():
    """Sheet FACE_93-97_data of 'Physiology Gas Exchange Water Relations.ods' (Wall et al., 2000, 2006)."""
    path = config.require(config.MARICOPA_GX_ODS, "Maricopa gas-exchange workbook")
    return _cached("maricopa_gx", lambda: pd.read_excel(path, "FACE_93-97_data", engine="odf"))


def quality_sheet():
    """Rows 11 onward, first 12 columns, of 'Quality Wheat Grain Baking 1993-1997.ods' (Kimball et al., 2001)."""
    path = config.require(config.MARICOPA_QUALITY_ODS, "Maricopa grain-quality workbook")
    return _cached("maricopa_quality", lambda: pd.read_excel(path, engine="odf", header=None, skiprows=10).iloc[:, :12])
