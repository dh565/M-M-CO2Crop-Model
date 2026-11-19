# recentre.py <runs dir> <decks dir>: grain-mineral calibration step. Column 17 avail_A_mgL_CAL of driver_minerals.csv is
# multiplied by obs / model of the control-deck grain concentration (daily_minerals.csv grain_mg_kg, last day) and written
# to every deck of the group (minerals are MG0-isolated and grain concentration is linear in avail_A). Controls and obs:
#   WDF25_* <- WDF25_Ruta_aCO2_WW Cu 4.5267, Mn 30.5467;  WDF23_* <- WDF23_Ruta_aCO2 Zn 15.08, Cu 2.25, Mn 16.3133;
#   BRAUNSCHWEIG_* <- Zn, mean factor of the 320 N controls 103 (obs 23.1069) and 115 (obs 25.0145).  Maricopa: none.
import sys, os, csv, glob
R, D = sys.argv[1:3]
def model(deck, mineral):
    rows = [r for r in csv.DictReader(open(os.path.join(R, deck, "output", "daily_minerals.csv")), skipinitialspace=True)]
    rows = [r for r in rows if r["mineral"].strip() == mineral]
    last = max(int(r["day"]) for r in rows); v = [r for r in rows if int(r["day"]) == last][0]
    return float(v["grain_mg_kg"]), sum(float(r["cap_binding"]) > 0 for r in rows)
def col17(path, mineral):
    for l in open(path):
        if l.startswith(mineral + ","): return float(l.split(",")[17])
GROUPS = [("WDF25_", "WDF25_Ruta_aCO2_WW", {"Cu": [("WDF25_Ruta_aCO2_WW", 4.5267)], "Mn": [("WDF25_Ruta_aCO2_WW", 30.5467)]}),
          ("WDF23_", "WDF23_Ruta_aCO2", {"Zn": [("WDF23_Ruta_aCO2", 15.08)], "Cu": [("WDF23_Ruta_aCO2", 2.25)], "Mn": [("WDF23_Ruta_aCO2", 16.3133)]}),
          ("BRAUNSCHWEIG_", "BRAUNSCHWEIG_103+115", {"Zn": [("BRAUNSCHWEIG_103", 23.1069), ("BRAUNSCHWEIG_115", 25.0145)]})]
for pre, ctrl, spec in GROUPS:
    decks = sorted(os.path.basename(p) for p in glob.glob(os.path.join(D, pre + "*")))
    ref = os.path.join(D, decks[0], "input", "driver_minerals.csv")
    notes, new = [], {}
    OLD = {mn: col17(ref, mn) for mn in spec}
    for mn, ctl in spec.items():
        old = OLD[mn]; fs, parts = [], []
        for dk, obs in ctl:
            m, cb = model(dk, mn); fs.append(obs / m); parts.append((dk, m, obs, obs / m, cb))
        f = sum(fs) / len(fs); new[mn] = float("%.7g" % (old * f))
        if len(parts) == 1:
            dk, m, obs, ff, cb = parts[0]
            notes.append("!   %s model %.4f obs %.4f factor %.5f cap_binding_days %d old %.7g new %.7g" % (mn, m, obs, ff, cb, old, new[mn]))
        else:
            notes.append("!   %s mean factor of decks %s old %.7g new %.7g" % (mn, " and ".join("%s (%.5f; model %.4f obs %.4f)" % (p[0][-3:], p[3], p[1], p[2]) for p in parts), old, new[mn]))
    hdr = "! Grain-mineral calibration step: column 17 avail_A_mgL_CAL rescaled by obs/model on control deck %s, closed-form rule (minerals are MG0-isolated, grain conc linear in avail_A):" % ctrl
    for dk in decks:
        p = os.path.join(D, dk, "input", "driver_minerals.csv"); L = open(p).read().split("\n"); out = []
        for l in L:
            mn = l.split(",")[0]
            if not l.startswith("!") and mn in new:
                t = l.split(","); assert abs(float(t[17]) - OLD[mn]) < 1e-12, (dk, mn); t[17] = "%.7g" % new[mn]; l = ",".join(t)
            out.append(l)
        open(p, "w").write("\n".join([hdr] + notes + out))
    print(pre, len(decks), "decks"); print(hdr); print("\n".join(notes))
