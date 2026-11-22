"""Share of grain-filling days on which the protein-coupled bound sets grain Fe loading (Part 2, Sect. 6.5 and 7:
"58 % of grain-filling days (80 % at Braunschweig)"). Per deck: the "days the sink bound" count of the Fe line of
"MINERAL SINK (mode 1)" in run.log, over the days with grain growth (W_grain_g_m2 rising, daily_yield_components.csv).
Sums over the decks of each site and over all decks.
usage: python fe_bound_days.py <paper_scripts dir> <out txt>   (REF_OUTPUTS from env)"""
import csv, glob, os, re, sys
SCR = sys.argv[1]; sys.path.insert(0, SCR)
import config
tot = {}
for d in sorted(glob.glob(os.path.join(config.REF_OUTPUTS, "*", ""))):
    deck = os.path.basename(os.path.dirname(d)); lg = os.path.join(d, "run.log")
    if not os.path.exists(lg): continue
    m = re.search(r"MINERAL SINK \(mode 1\) Fe .*?days the sink bound =\s*(\d+)", open(lg).read())
    if not m: continue
    W = [float(r["W_grain_g_m2"]) for r in csv.DictReader(open(os.path.join(d, "output", "daily_yield_components.csv")))]
    gf = sum(1 for a, b in zip(W, W[1:]) if b > a)
    t = tot.setdefault(deck.split("_")[0][:5], [0, 0]); t[0] += int(m.group(1)); t[1] += gf
a = sum(v[0] for v in tot.values()); b = sum(v[1] for v in tot.values())
out = ["%-6s %5d of %5d days  %.0f %%" % (s, v[0], v[1], 100 * v[0] / v[1]) for s, v in sorted(tot.items())]
out.append("%-6s %5d of %5d days  %.0f %%" % ("all", a, b, 100 * a / b))
print("\n".join(out))
open(sys.argv[2], "w").write("\n".join(out) + "\n")
