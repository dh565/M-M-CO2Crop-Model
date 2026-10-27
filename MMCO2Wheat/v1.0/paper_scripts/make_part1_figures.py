#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent.
"""Part 1, Figs. 1-5.  Fig. 1 is a schematic (no data).  Figs. 2-5 read SENS_DIR/all_metrics.json, written by
run_sensitivity.py (125 model runs with MODEL_EXE).  Output: OUT_DIR/part1_figures/<stem>.pdf|png and gmd/f01..f05.

    python3 make_part1_figures.py [--only 1,2,...]

Figure -> function (file stem):
   1  p1_figures.fig1   P1_Fig01_model_structure
   2  p1_figures.fig2   P1_Fig02_driver_responses
   3  p1_figures.fig3   P1_Fig03_alt_cascade_elasticity
   4  p1_figures.fig4   P1_Fig04_CO2_interactions
   5  p1_figures.fig5   P1_Fig05_alt_matrix_tornado
"""
import argparse
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "lib"))
import config
import p1_figures as P1

FIGURES = {1: (P1.fig1, "P1_Fig01_model_structure"), 2: (P1.fig2, "P1_Fig02_driver_responses"),
           3: (P1.fig3, "P1_Fig03_alt_cascade_elasticity"), 4: (P1.fig4, "P1_Fig04_CO2_interactions"),
           5: (P1.fig5, "P1_Fig05_alt_matrix_tornado")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="comma list of figure numbers")
    a = ap.parse_args()
    out = config.out_subdir("part1_figures")
    gmd = config.out_subdir("part1_figures/gmd")
    nums = [int(x) for x in a.only.split(",")] if a.only else [1, 2, 3, 4, 5]
    if any(n > 1 for n in nums):
        metrics = json.load(open(config.require(os.path.join(config.SENS_DIR, "all_metrics.json"), "all_metrics.json of the sensitivity runs")))
        P1.set_metrics(metrics)
    for n in nums:
        fn, stem = FIGURES[n]
        print("Fig. %d  %s" % (n, stem))
        fn(out)
        for ext in ("pdf", "png"):
            shutil.copy(os.path.join(out, "%s.%s" % (stem, ext)), os.path.join(gmd, "f%02d.%s" % (n, ext)))
    print("figures written to", out)


if __name__ == "__main__":
    main()
