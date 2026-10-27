#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""All paths used by the figure and table scripts, in one place.

Every path can be overridden by an environment variable of the same name.  Nothing in
this folder is written to; results go to OUT_DIR, which defaults to a sibling folder
``generated_results/`` next to ``paper_scripts/``.

  REF_OUTPUTS  the reference-outputs archive MMCO2Wheat_v1.0_FACE_reference_outputs
               (38 deck folders, each with input/, output/ and run.log)
  OBS_DIR      the observation data that are not redistributed here (see README.md,
               section "Observation data")
  OUT_DIR      where figures, tables and caches are written
  SENS_DIR     where the Part 1 sensitivity runs and their all_metrics.json live
  MODEL_EXE    the model executable built from the code package (Part 1 Figs. 2-5 only)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)

REF_OUTPUTS = os.environ.get("REF_OUTPUTS", os.path.join(PARENT, "MMCO2Wheat_v1.0_FACE_reference_outputs"))
OBS_DIR = os.environ.get("OBS_DIR", os.path.join(PARENT, "observation_data"))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(PARENT, "generated_results"))
SENS_DIR = os.environ.get("SENS_DIR", os.path.join(OUT_DIR, "sensitivity_runs"))
MODEL_EXE = os.environ.get("MODEL_EXE", "")

# Observation tables redistributed with the scripts (derived from the open Braunschweig and
# Maricopa data sets; see observations/README_observations.md)
SHIPPED_OBS = os.path.join(HERE, "observations")
BRAUNSCHWEIG_OBS_DIR = os.path.join(SHIPPED_OBS, "braunschweig")           # BRAUNSCHWEIG_<trno>_obs.csv
MARICOPA_OBS_DIR = os.path.join(SHIPPED_OBS, "maricopa")                   # MARICOPA_*.csv, maricopa_archive_*.csv
OBS_REGISTRY_OPEN = os.path.join(SHIPPED_OBS, "obs_registry_braunschweig_maricopa.csv")
CONSOLIDATED_OBS = os.path.join(SHIPPED_OBS, "consolidated_obs_braunschweig_maricopa.csv")
MARICOPA_EB_HOURLY = os.path.join(MARICOPA_OBS_DIR, "maricopa_eb_hourly.csv")

# Observation files the reader supplies under OBS_DIR (relative layout documented in README.md)
WDF_DIR = os.path.join(OBS_DIR, "WheatDryFACE")
WDF25_DIR = os.path.join(WDF_DIR, "WDF25")       # S3_<genotype>_<CO2>_<water>.csv, season 2024/25, one file per treatment
WDF25_V3_DIR = WDF25_DIR                         # the same files: thermal-image columns (Tleaf_flir, DTleaf_air_flir)
WDF25_PORO_DIR = WDF25_DIR                       # the same files: porometer columns (gs_poro, Tr_poro, Tleaf_poro, DTleaf_air_poro)
WDF23_DIR = os.path.join(WDF_DIR, "WDF23")       # WDF23_<genotype>_<CO2>.csv, season 2022/23
OBS_REGISTRY_WDF = os.path.join(WDF_DIR, "obs_registry_wheatdryface.csv")
MARICOPA_ARCHIVE_DIR = os.path.join(OBS_DIR, "Maricopa_FACE_archive")   # the .ods workbooks of Kimball et al. (2017)
MARICOPA_QUALITY_ODS = os.path.join(MARICOPA_ARCHIVE_DIR, "Quality Wheat Grain Baking 1993-1997.ods")
MARICOPA_GX_ODS = os.path.join(MARICOPA_ARCHIVE_DIR, "Physiology Gas Exchange Water Relations.ods")
MARICOPA_EB_ODS = {
    "eb94": os.path.join(MARICOPA_ARCHIVE_DIR, "Weather Canopy & Soil Temperatures Energy Balance 1994 hourly.ods"),
    "eb96": os.path.join(MARICOPA_ARCHIVE_DIR, "Weather Canopy & Soil Temperatures Energy Balance CO2 1996 15-min.ods"),
    "eb97": os.path.join(MARICOPA_ARCHIVE_DIR, "Weather Canopy & Soil Temperatures Energy Balance CO2 O3 1997 15-min.ods"),
}

# Caches of parsed .ods workbooks (slow to read); written under OUT_DIR, never under the archive
CACHE_DIR = os.path.join(OUT_DIR, "_cache")


def out_subdir(name):
    """OUT_DIR/<name>, created on demand."""
    p = os.path.join(OUT_DIR, name)
    os.makedirs(p, exist_ok=True)
    return p


def require(path, what):
    """Raise a clear error when an input the scripts need is absent."""
    if not os.path.exists(path):
        raise FileNotFoundError("%s not found at %s (set the environment variable described in config.py / README.md)" % (what, path))
    return path
