"""Soil mineral N at 0-40 cm at maturity at Braunschweig, without and with the crop-driven immobilization (Part 2,
Sect. 5: "lowers the root-mean-square difference of these eight values from ... g N m-2 without the immobilization to
... g N m-2").
The eight decks are the deficient (N1) and adequate (N2) nitrogen treatments at both CO2 levels in both seasons.
Without: each deck is copied from decks/ and run with driver record 197 (rhizo_immob) set to 0 at run time; no other
record or file is changed. With: the same decks of REF_OUTPUTS (record 197 = "1 0.5", as shipped).
Model: NO3 + NH4 of soil-N layers 1-4 (0-0.40 m, 0.10 m each at Braunschweig, driver_swtr.txt) of
daily_soilN_layers.csv on the maturity day of daily_phenology.csv (the last simulated day if no maturity is reached).
Observed: Dier et al. (2019), Field Crops Res. 234, 107-118, Fig. 5b, values read from the figure (g N m-2).
RMSD over the eight decks.
usage: python soil_nmin_immobilization.py <paper_scripts dir> <mmco2crop> <work dir> <out txt> [workers]
       (REF_OUTPUTS from env; the runs without immobilization are written to <work dir>/<deck>)"""
import csv, math, os, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
SCR, EXE, WORK, OUT = sys.argv[1:5]; NW = int(sys.argv[5]) if len(sys.argv) > 5 else 8
sys.path.insert(0, SCR)
import config
DECKS_DIR = os.path.join(os.path.dirname(os.path.abspath(SCR)), "decks")
REC_RHIZO = 197
NMIN_OBS = {"BRAUNSCHWEIG_101": 2.0, "BRAUNSCHWEIG_104": 2.0, "BRAUNSCHWEIG_113": 2.3, "BRAUNSCHWEIG_116": 1.3,   # N1
            "BRAUNSCHWEIG_102": 1.7, "BRAUNSCHWEIG_105": 1.4, "BRAUNSCHWEIG_114": 2.4, "BRAUNSCHWEIG_117": 1.8}   # N2
LAB = {"BRAUNSCHWEIG_101": "2014 aCO2 N1", "BRAUNSCHWEIG_104": "2014 eCO2 N1", "BRAUNSCHWEIG_113": "2015 aCO2 N1",
       "BRAUNSCHWEIG_116": "2015 eCO2 N1", "BRAUNSCHWEIG_102": "2014 aCO2 N2", "BRAUNSCHWEIG_105": "2014 eCO2 N2",
       "BRAUNSCHWEIG_114": "2015 aCO2 N2", "BRAUNSCHWEIG_117": "2015 eCO2 N2"}


def run_off(deck):
    w = os.path.join(WORK, deck); shutil.rmtree(w, ignore_errors=True); os.makedirs(os.path.join(w, "output"))
    shutil.copytree(os.path.join(DECKS_DIR, deck, "input"), os.path.join(w, "input"))
    p = os.path.join(w, "input", "driver_main.txt"); L = open(p).read().split("\n")
    L[REC_RHIZO - 1] = "0" + L[REC_RHIZO - 1][len(L[REC_RHIZO - 1].split("!")[0].rstrip()):]
    open(p, "w").write("\n".join(L))
    with open(os.path.join(w, "run.log"), "w") as fh:
        rc = subprocess.call([os.path.abspath(EXE)], cwd=w, stdout=fh, stderr=subprocess.STDOUT)
        fh.write("exit=%d\n" % rc)
    if rc: raise RuntimeError("%s: exit %d" % (deck, rc))


def nmin040(d):
    o = os.path.join(d, "output")
    ph = list(csv.DictReader(open(os.path.join(o, "daily_phenology.csv")), skipinitialspace=True))
    s = {int(float(r["day"])): r for r in csv.DictReader(open(os.path.join(o, "daily_soilN_layers.csv")), skipinitialspace=True)}
    mat = next((int(float(r["day"])) for r in ph if r["stage_name"].strip() == "Maturity"), max(s))
    r = s[min(mat, max(s))]
    return sum(float(r["NO3_l%d" % i]) + float(r["NH4_l%d" % i]) for i in range(1, 5))


os.makedirs(WORK, exist_ok=True)
with ThreadPoolExecutor(NW) as ex: list(ex.map(run_off, NMIN_OBS))
out = ["%-18s %-14s %6s %8s %8s" % ("deck", "treatment", "obs", "without", "with")]
dv = {"without": [], "with": []}
for deck, ob in NMIN_OBS.items():
    off = nmin040(os.path.join(WORK, deck)); on = nmin040(os.path.join(config.REF_OUTPUTS, deck))
    dv["without"].append(off - ob); dv["with"].append(on - ob)
    out.append("%-18s %-14s %6.2f %8.2f %8.2f" % (deck, LAB[deck], ob, off, on))
rmsd = {k: math.sqrt(sum(x * x for x in v) / len(v)) for k, v in dv.items()}
out.append("RMSD (g N m-2): without immobilization %.2f, with %.2f" % (rmsd["without"], rmsd["with"]))
print("\n".join(out))
open(OUT, "w").write("\n".join(out) + "\n")
