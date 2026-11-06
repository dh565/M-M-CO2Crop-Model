"""Timing of the low-P effect against Rodriguez et al. (field, Bray-1 5.5 mg/kg).
Usage: compare_rodriguez.py <run_dir from run_co2xP.sh>
Targets (evaluation only, never model inputs):
  Rodriguez, Andrade & Goudriaan 2000 FCR 67:227, abstract: at 61 DAE, -70 % above-ground
  biomass and -73 % LAI (non-shaded), low-P vs P-fertilized.
  Rodriguez, Pomar & Goudriaan 1998 Plant Soil 202:149, abstract: -25 % shoot growth at 44 DAE
  (sand culture), longer phyllochron, delayed or suppressed tillers.
"""
import csv, sys
R = sys.argv[1]
def ser(tag, f, c):
    return {int(r["day"]): float(r[c]) for r in csv.DictReader(open(f"{R}/{tag}_420/output/{f}"))}
e = min(ser("P40", "daily_P.csv", "Pvac_gP_m2"))
print(f"emergence day {e}")
print("DAE  | biomass loss P8 / P5 / P3 | LAI loss P5 | axes m-2 P40 / P5 | Haun stage P40 / P5")
for dae in (20, 30, 44, 61, 80):
    d = e + dae
    B = {t: ser(t, "daily_fluxes_full.csv", "Biomass_g_m2")[d] for t in ("P40", "P8", "P5", "P3")}
    L = {t: ser(t, "daily_d2_organs.csv", "LAI_used").get(d, float("nan")) for t in ("P40", "P5")}
    A = {t: ser(t, "daily_canopy_tillering.csv", "axes_live_m2").get(d, float("nan")) for t in ("P40", "P5")}
    H = {t: ser(t, "daily_canopy_ranks.csv", "HS").get(d, float("nan")) for t in ("P40", "P5")}
    f = lambda a, b: 100 * (1 - a / b) if b > 0 else float("nan")
    print(f"{dae:4d} | {f(B['P8'],B['P40']):5.0f} {f(B['P5'],B['P40']):5.0f} {f(B['P3'],B['P40']):5.0f} % "
          f"| {f(L['P5'],L['P40']):5.0f} % | {A['P40']:6.0f} {A['P5']:6.0f} | {H['P40']:5.2f} {H['P5']:5.2f}")
