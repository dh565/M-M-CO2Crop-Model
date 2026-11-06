#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Part 2, Figs. 1-11: draws every figure from the reference-outputs archive and writes it twice,
under its descriptive name (OUT_DIR/part2_figures/<stem>.pdf|png) and as f01..f11 (OUT_DIR/part2_figures/gmd/).

    python3 make_part2_figures.py [--only 1,2,...]

Figure -> function (file stem):
   1  paperfigs2.fig_eb         FigN01_energy_balance_all_sites
   2  paperfigs.fig02           Fig02_soil_water_temperature
   3  paperfigs2.fig_gx         Fig03_flag_leaf_gas_exchange_all_sites
   4  paperfigs.figN_canopy     FigN04_leaf_area_all_sites
   5  paperfigs.figN_biomass    FigN05_biomass_all_sites
   6  paperfigs.figN_harvest    FigN06_harvest_all_sites
   7  paperfigs.figN_eco2_agro  FigN07_eCO2_response_agronomy_all_sites
   8  paperfigs.figN_eco2_quality FigN08_eCO2_response_quality_all_sites
   9  paperfigs.figN_water_N    FigN09_water_and_N_responses_all_sites
  10  paperfigs2.fig_tcan       FigN10_canopy_temperature_all_sites
  11  paperfigs.figN_signs_site2 FigN11_sign_agreement_site_two_bars
"""
import argparse
import os
import shutil
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
import config
import paperfigs as PF
import paperfigs2 as PF2

FIGURES = OrderedDict([
    (1, (PF2.fig_eb, "FigN01_energy_balance_all_sites")),
    (2, (PF.fig02, "Fig02_soil_water_temperature")),
    (3, (PF2.fig_gx, "Fig03_flag_leaf_gas_exchange_all_sites")),
    (4, (PF.figN_canopy, "FigN04_leaf_area_all_sites")),
    (5, (PF.figN_biomass, "FigN05_biomass_all_sites")),
    (6, (PF.figN_harvest, "FigN06_harvest_all_sites")),
    (7, (PF.figN_eco2_agro, "FigN07_eCO2_response_agronomy_all_sites")),
    (8, (PF.figN_eco2_quality, "FigN08_eCO2_response_quality_all_sites")),
    (9, (PF.figN_water_N, "FigN09_water_and_N_responses_all_sites")),
    (10, (PF2.fig_tcan, "FigN10_canopy_temperature_all_sites")),
    (11, (PF.figN_signs_site2, "FigN11_sign_agreement_site_two_bars")),
])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="comma list of figure numbers")
    a = ap.parse_args()
    config.require(config.REF_OUTPUTS, "reference-outputs archive")
    out = config.out_subdir("part2_figures")
    gmd = config.out_subdir("part2_figures/gmd")
    PF.init(config.REF_OUTPUTS)
    nums = [int(x) for x in a.only.split(",")] if a.only else list(FIGURES)
    for n in nums:
        fn, stem = FIGURES[n]
        print("Fig. %d  %s" % (n, stem))
        fn(out)
        for ext in ("pdf", "png"):
            shutil.copy(os.path.join(out, "%s.%s" % (stem, ext)), os.path.join(gmd, "f%02d.%s" % (n, ext)))
    print("figures written to", out)


if __name__ == "__main__":
    main()
