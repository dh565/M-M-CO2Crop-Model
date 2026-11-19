"""Treatment-level Table 4 series rows (method of table4_agg.py, extended to all series quantities).
Per deck: obs and model averaged over the paired dates of ebwb_gx.pairs; n = number of treatments (decks with >=1 pair);
Obs/Sim = mean of treatment means; Bias = mean(Sim-Obs); RMSD over treatment means;
1sigma = fraction with |Sim-Obs| <= mean daily obs sd (15 % of |Obs| where no sd);
soil-water rows: |Sim-Obs| <= 0.02 m3 m-3 (the rule of the published Table 4 soil-water rows).
Maricopa soil temperature 10/40 cm: lib/maricopa_tsoil obs_daily/model_daily, paired by date, 15 % rule.
usage: python table4_series.py <paper_scripts dir> <out json>   (REF_OUTPUTS, OBS_DIR, OUT_DIR from env)"""
import sys, os, json, math
SCR = sys.argv[1]; sys.path.insert(0, SCR); sys.path.insert(0, os.path.join(SCR, 'lib'))
import config
import paperfigs as PF
PF.init(config.REF_OUTPUTS)
import ebwb_gx as E
def mean(x):
    x = list(x); return sum(x)/len(x)
KEYS = ["SWC20","SWC55","Tsoil20","Tsoil55","Rn_day","Rn_night","LWout_day","LWout_night","SWout_day","albedo_day",
        "G_24h","Tsurf_rad","d0_20","d20_40","MAR_SWC_0_30","MAR_SWC_30_50","G_day","G_night","LWin_day","LWin_night","SWin_day"]
SWC = {'SWC20','SWC55','d0_20','d20_40','MAR_SWC_0_30','MAR_SWC_30_50'}; SWC_TOL = 0.02
sites = {}
for d, info in E.DECKS.items(): sites.setdefault(info["site"], []).append(d)
out = {}
for site in sorted(sites):
    for q in KEYS:
        rows = []
        for deck in sorted(sites[site]):
            pr = E.pairs(deck, q)
            if not pr: continue
            sd = [p[2] for p in pr if not (isinstance(p[2], float) and math.isnan(p[2]))]
            rows.append((deck, len(pr), mean(p[1] for p in pr), mean(p[3] for p in pr), mean(sd) if sd else float('nan')))
        if not rows: continue
        n = len(rows); ob = [r[2] for r in rows]; md = [r[3] for r in rows]
        bias = mean(m-o for o, m in zip(ob, md)); rmsd = math.sqrt(mean((m-o)**2 for o, m in zip(ob, md)))
        within = [abs(m-o) <= (SWC_TOL if q in SWC else (s if not math.isnan(s) else 0.15*abs(o))) for o, m, s in zip(ob, md, [r[4] for r in rows])]
        within_sd = [abs(m-o) <= (s if not math.isnan(s) else 0.15*abs(o)) for o, m, s in zip(ob, md, [r[4] for r in rows])]
        out[f"{site}/{q}"] = dict(n=n, obs=mean(ob), sim=mean(md), bias=bias, rmsd=rmsd, frac1sd=sum(within)/n,
                                  npairs=sum(r[1] for r in rows), frac1sd_obs_sd_rule=sum(within_sd)/n, per_deck={r[0]: r[1:] for r in rows})
        print(f"{site:12s} {q:14s} n={n} Obs={mean(ob):.4g} Sim={mean(md):.4g} Bias={bias:+.4g} RMSD={rmsd:.4g} 1sd={sum(within)/n:.2f} pairs={out[f'{site}/{q}']['npairs']}")
try:
    import maricopa_tsoil as MT
    O = MT.obs_daily()
    for z in (0.10, 0.40):
        rows = []
        for trno in sorted(O):
            o = O[trno]; m = MT.model_daily(trno)
            j = o[o.z == z].merge(m[m.z == z], on="date", suffixes=("_o", "_m"))
            if len(j): rows.append(("MARICOPA_%d" % trno, len(j), float(j.T_o.mean()), float(j.T_m.mean()), float('nan')))
        n = len(rows); ob = [r[2] for r in rows]; md = [r[3] for r in rows]
        bias = mean(m-o for o, m in zip(ob, md)); rmsd = math.sqrt(mean((m-o)**2 for o, m in zip(ob, md)))
        within = [abs(m-o) <= 0.15*abs(o) for o, m in zip(ob, md)]
        q = "MAR_Tsoil%d" % int(z*100)
        out["MARICOPA/" + q] = dict(n=n, obs=mean(ob), sim=mean(md), bias=bias, rmsd=rmsd, frac1sd=sum(within)/n, frac1sd_obs_sd_rule=sum(within)/n,
                                    npairs=sum(r[1] for r in rows), per_deck={r[0]: r[1:] for r in rows})
        print(f"MARICOPA     {q:14s} n={n} Obs={mean(ob):.4g} Sim={mean(md):.4g} Bias={bias:+.4g} RMSD={rmsd:.4g} 1sd={sum(within)/n:.2f}")
except Exception as e:
    print("Maricopa Tsoil failed:", repr(e))
json.dump(out, open(sys.argv[2], "w"), indent=1)
