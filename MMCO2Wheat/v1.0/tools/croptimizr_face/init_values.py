# init_values.py <decks> <initial_values.csv>: write the initial values of the calibration to every deck listed
# (phenology thresholds and Blade_L_max_cm in phenology_cultivar.txt, fruiting_efficiency in cultivar_yield.csv,
# column 17 avail_A_mgL_CAL of the re-centred minerals in driver_minerals.csv). Values are written as given.
# The anthesis and grain-fill values equal A0 and G0 of targets.csv, Blade_L_max_cm L0 and fruiting_efficiency FE0.
import sys, csv, os
D, F = sys.argv[1:3]
edits = {}
for r in csv.DictReader(open(F)):
    edits.setdefault((r["deck"], r["file"]), {})[r["name"]] = r["value"]
for (deck, fn), new in sorted(edits.items()):
    p = os.path.join(D, deck, "input", fn); out = []; done = set()
    for l in open(p).read().split("\n"):
        if l and not l.startswith("!"):
            if fn == "phenology_cultivar.txt":
                t = l.split()
                if len(t) >= 2 and t[1] in new: l = new[t[1]] + l[len(t[0]):]; done.add(t[1])
            elif fn == "cultivar_yield.csv" and ",fruiting_efficiency" in l:
                l = new["fruiting_efficiency"] + l[l.index(","):]; done.add("fruiting_efficiency")
            elif fn == "driver_minerals.csv":
                k = l.split(",")[0] + ":avail_A_mgL_CAL"
                if k in new: t = l.split(","); t[17] = new[k]; l = ",".join(t); done.add(k)
        out.append(l)
    assert done == set(new), (deck, fn, set(new) - done)
    open(p, "w").write("\n".join(out))
print(len(edits), "deck files set to the initial values")
