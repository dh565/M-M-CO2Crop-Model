"""CO2 x water and CO2 x nitrogen interactions of the observations (Part 2, Sect. 6.7 and the CO2 x supply sentences):
two-way ANOVA of the interaction term from the treatment means, SDs and replicate numbers of each design (balanced:
SS_AB from the cell means x n; SS_error = sum over cells of (n - 1) sd^2). Observed values only (paperfigs.obs).
Designs: WheatDryFACE 2024/25 CO2 x water per genotype (n = 3 rings), Braunschweig 2014 and 2015 CO2 x nitrogen (3 levels,
n = 3 rings), Maricopa 1993/94 CO2 x water and 1996/97 CO2 x nitrogen (n = 4 replicates).
usage: python anova_summary.py <paper_scripts dir> <out txt>   (REF_OUTPUTS, OBS_DIR, OUT_DIR from env)"""
import sys, os
SCR = sys.argv[1]; sys.path.insert(0, SCR); sys.path.insert(0, os.path.join(SCR, 'lib'))
import config
import paperfigs as PF
from scipy import stats
PF.init(config.REF_OUTPUTS)


def anova(cells, n):
    """cells[(i, j)] = (mean, sd); i = CO2 level, j = supply level. F and P of the i x j interaction."""
    I = sorted({k[0] for k in cells}); J = sorted({k[1] for k in cells})
    gm = sum(m for m, s in cells.values()) / len(cells)
    mi = {i: sum(cells[i, j][0] for j in J) / len(J) for i in I}; mj = {j: sum(cells[i, j][0] for i in I) / len(I) for j in J}
    ss_ab = n * sum((cells[i, j][0] - mi[i] - mj[j] + gm) ** 2 for i in I for j in J)
    df_ab = (len(I) - 1) * (len(J) - 1); ss_e = sum((n - 1) * s ** 2 for m, s in cells.values()); df_e = len(cells) * (n - 1)
    F = (ss_ab / df_ab) / (ss_e / df_e)
    return F, stats.f.sf(F, df_ab, df_e)


D = {"WDF25 Ruta CO2 x water": ({("a", w): "WDF25_Ruta_aCO2_%s" % w for w in ("WW", "WL")} | {("e", w): "WDF25_Ruta_eCO2_%s" % w for w in ("WW", "WL")}, 3),
     "WDF25 Zahir CO2 x water": ({("a", w): "WDF25_Zahir_aCO2_%s" % w for w in ("WW", "WL")} | {("e", w): "WDF25_Zahir_eCO2_%s" % w for w in ("WW", "WL")}, 3),
     "BRA2014 CO2 x N": ({("a", 1): "BRAUNSCHWEIG_101", ("a", 2): "BRAUNSCHWEIG_102", ("a", 3): "BRAUNSCHWEIG_103",
                          ("e", 1): "BRAUNSCHWEIG_104", ("e", 2): "BRAUNSCHWEIG_105", ("e", 3): "BRAUNSCHWEIG_106"}, 3),
     "BRA2015 CO2 x N": ({("a", 1): "BRAUNSCHWEIG_113", ("a", 2): "BRAUNSCHWEIG_114", ("a", 3): "BRAUNSCHWEIG_115",
                          ("e", 1): "BRAUNSCHWEIG_116", ("e", 2): "BRAUNSCHWEIG_117", ("e", 3): "BRAUNSCHWEIG_118"}, 3),
     "MAR9394 CO2 x water": ({("a", "ww"): "MARICOPA_906", ("a", "wl"): "MARICOPA_905", ("e", "ww"): "MARICOPA_908", ("e", "wl"): "MARICOPA_907"}, 4),
     "MAR9697 CO2 x N": ({("a", "hi"): "MARICOPA_914", ("a", "lo"): "MARICOPA_913", ("e", "hi"): "MARICOPA_916", ("e", "lo"): "MARICOPA_915"}, 4)}
out = []
for name, (dd, n) in D.items():
    for q in ("LAI_peak", "grain_N_pct", "kernel_number", "AGB_harvest", "grain_yield", "TKW", "HI"):
        cells = {k: PF.obs(v, q) for k, v in dd.items()}
        if any(m != m or s != s for m, s in cells.values()):
            out.append("%-24s %-14s not observed in every cell" % (name, q)); continue
        F, p = anova(cells, n)
        out.append("%-24s %-14s F=%.2f P=%.3f" % (name, q, F, p))
print("\n".join(out))
open(sys.argv[2], "w").write("\n".join(out) + "\n")
