# CO2, water and N effects (treatment ratios): observations vs the model, one or more simulation sets.
# Usage: OBS_DIR=<obs_latest> python3 effects.py <paper_scripts dir> <out.md> label=path [label=path ...]
import os, sys, math, csv
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, os.path.join(sys.argv[1], "lib"))
import face_results as F
OUT = sys.argv[2]
ARMS = [a.split("=", 1) for a in sys.argv[3:]]
obs = F.Obs()
arms = {lab: F.load_arm(p, list(F.DECKS)) for lab, p in ARMS}
Q = [("grain_yield", "grain yield"), ("AGB_harvest", "biomass at harvest"), ("kernel_number", "kernels m-2"),
     ("crop_N_harvest", "crop N"), ("grain_N_pct", "grain N %"), ("LAI_peak", "peak LAI")]
nan = float("nan")


def mv(arm, d, q):
    m = arm.get(d)
    if not m:
        return nan
    if q == "crop_N_harvest":
        nl = m.n[-1]
        return sum(F.fl(nl[k]) for k in ("N_leaf_gN_m2", "N_dead_leaf_gN_m2", "N_stem_gN_m2", "N_storage_gN_m2"))
    return m.v.get(q, nan)


def ov(d, q):
    o = obs.get(d, q)
    return o["value"] if o and o.get("value") is not None else nan


D = F.DECKS
key = lambda d, drop: tuple(D[d][k] for k in ("site", "season", "cultivar", "co2", "water", "nlevel") if k != drop)
pairs = {"CO2 (eCO2 / aCO2)": [], "Water (limited / well-watered)": [], "N (low / high)": []}
for t in D:
    for b in D:
        if D[t]["co2"] == "eCO2" and D[b]["co2"] == "aCO2" and key(t, "co2") == key(b, "co2"):
            pairs["CO2 (eCO2 / aCO2)"].append((t, b))
        if D[t]["water"] == "WL" and D[b]["water"] == "WW" and key(t, "water") == key(b, "water"):
            pairs["Water (limited / well-watered)"].append((t, b))
for y, sets in (("13/14", (("101", "102"), ("104", "105"), ("103", "102"), ("106", "105"))),
                ("14/15", (("113", "114"), ("116", "117"), ("115", "114"), ("118", "117")))):
    for a, b in sets:
        pairs["N (low / high)"].append(("BRAUNSCHWEIG_" + a, "BRAUNSCHWEIG_" + b))
pairs["N (low / high)"] += [("MARICOPA_913", "MARICOPA_914"), ("MARICOPA_915", "MARICOPA_916")]


def lab(t, b):
    s = D[t]["site"].replace("BRAUNSCHWEIG", "BRA").replace("MARICOPA", "MAR")
    if D[t]["co2"] != D[b]["co2"]:
        x = "%s %s %s %s %s" % (s, D[t]["season"], D[t]["cultivar"], D[t]["water"], D[t]["nlevel"].split(" (")[0])
    elif D[t]["water"] != D[b]["water"]:
        x = "%s %s %s %s" % (s, D[t]["season"], D[t]["cultivar"], D[t]["co2"])
    else:
        x = "%s %s %s %s/%s" % (s, D[t]["season"], D[t]["co2"], D[t]["nlevel"].split(" (")[0], D[b]["nlevel"].split(" (")[0])
    return x


names = [l for l, _ in ARMS]
o, summ, rows = [], [], []
P = o.append
P("# CO2, water and N effects: observed vs model (M&M-CO2Wheat v1.0)\n")
P("Effect = treatment / control ratio. Obs = latest observation set. Models: " + ", ".join(names) + ".")
P("'Error' = mean |model effect - observed effect| in % points over the pairs with an observation. '-' = not measured.\n")
for eff, pl in pairs.items():
    P("## " + eff + "\n")
    for q, qn in Q:
        rr = []
        for t, b in sorted(pl):
            r = [ov(t, q) / ov(b, q)] + [mv(arms[n], t, q) / mv(arms[n], b, q) for n in names]
            rr.append((lab(t, b), r))
            rows.append([eff, qn, lab(t, b)] + ["%.4f" % x for x in r])
        if all(math.isnan(x) for _, r in rr for x in r[1:]):
            continue
        P("**%s**\n" % qn)
        P("| pair | obs | " + " | ".join(names) + " |")
        P("|---|---|" + "---|" * len(names))
        for l, r in rr:
            f = lambda x: "-" if math.isnan(x) else "%+.0f%%" % ((x - 1) * 100)
            P("| %s | %s | %s |" % (l, f(r[0]), " | ".join(f(x) for x in r[1:])))
        err = []
        for i in range(len(names)):
            e = [abs(r[i + 1] - r[0]) * 100 for _, r in rr if not math.isnan(r[0]) and not math.isnan(r[i + 1])]
            err.append((sum(e) / len(e), len(e)) if e else (nan, 0))
        if err[0][1]:
            P("| **Error (n=%d)** | | %s |" % (err[0][1], " | ".join("%.1f" % e for e, _ in err)))
            summ.append((eff, qn, err[0][1], [e for e, _ in err]))
        P("")
P("## Summary: mean error of the effect (% points; lower = closer to observed)\n")
P("| effect | variable | n | " + " | ".join(names) + " |")
P("|---|---|---|" + "---|" * len(names))
for eff, qn, n, e in summ:
    P("| %s | %s | %d | %s |" % (eff, qn, n, " | ".join("%.1f" % x for x in e)))
open(OUT, "w").write("\n".join(o) + "\n")
with open(OUT.replace(".md", ".csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["effect", "variable", "pair", "obs"] + names)
    w.writerows(rows)
print("\n".join(o))
