"""Table 4 gas-exchange and thermal-image rows that make_part2_tables.py does not write: per deck, obs and model averaged over the
observation dates of ebwb_gx.pairs; n = number of treatments; Bias = mean(Sim-Obs); RMSD over treatment means;
1-sigma = |Sim-Obs| <= mean daily sd (15 % where no sd). Prints the leaf gas-exchange rows of every site and the WDF25 Tleaf_FLIR row.
Also prints the pooled per-date statistics (n, r) and the per-date CO2 effects used in Fig. 3g-l / Sect. 6.6.1.
usage: python table4_agg.py <paper_scripts dir> <out json>   (REF_OUTPUTS, OBS_DIR, OUT_DIR from env)"""
import sys, os, json, math
SCR = sys.argv[1]; sys.path.insert(0, SCR); sys.path.insert(0, os.path.join(SCR, 'lib'))
import config
import paperfigs as PF
PF.init(config.REF_OUTPUTS)
import ebwb_gx as E
def mean(x):
    x = list(x); return sum(x)/len(x)
def corr(x, y):
    mx, my = mean(x), mean(y); sx = math.sqrt(sum((a-mx)**2 for a in x)); sy = math.sqrt(sum((b-my)**2 for b in y))
    return sum((a-mx)*(b-my) for a, b in zip(x, y))/(sx*sy) if sx and sy else float('nan')
sites = {}
for d, info in E.DECKS.items(): sites.setdefault(info["site"], []).append(d)
keys = [k for k, v in E.Q.items() if v[3] == "gx"] + ["Tleaf_FLIR"]
out = {}
for site in sorted(sites):
    for q in keys:
        rows = []; pooled = []
        for deck in sorted(sites[site]):
            pr = E.pairs(deck, q)
            if not pr: continue
            ob = [p[1] for p in pr]; md = [p[3] for p in pr]; sd = [p[2] for p in pr if not (isinstance(p[2], float) and math.isnan(p[2]))]
            rows.append((deck, len(pr), mean(ob), mean(md), mean(sd) if sd else float('nan'), [p[0] for p in pr]))
            pooled += [(deck, p[0], p[1], p[3]) for p in pr]
        if not rows: continue
        n = len(rows); ob = [r[2] for r in rows]; md = [r[3] for r in rows]
        bias = mean(m - o for o, m in zip(ob, md)); rmsd = math.sqrt(mean((m - o)**2 for o, m in zip(ob, md)))
        within = [abs(m - o) <= (s if not math.isnan(s) else 0.15*abs(o)) for (o, m, s) in zip(ob, md, [r[4] for r in rows])]
        px = [p[2] for p in pooled]; py = [p[3] for p in pooled]
        rec = dict(n_treat=n, obs=mean(ob), sim=mean(md), bias=bias, rmsd=rmsd, frac1sd=sum(within)/n,
                   pooled_n=len(pooled), pooled_r=corr(px, py), pooled_bias=mean(b-a for a, b in zip(px, py)),
                   pooled_rmsd=math.sqrt(mean((b-a)**2 for a, b in zip(px, py))),
                   dates=sorted({p[1] for p in pooled}), per_deck={r[0]: dict(n=r[1], obs=r[2], sim=r[3], sd=r[4]) for r in rows})
        out[f"{site}/{q}"] = rec
        print(f"{site:6s} {q:12s} n={n} Obs={mean(ob):.4g} Sim={mean(md):.4g} Bias={bias:+.4g} RMSD={rmsd:.4g} 1sd={sum(within)/n:.3f} | pooled n={len(pooled)} r={corr(px,py):.3f} bias={rec['pooled_bias']:+.3g} rmsd={rec['pooled_rmsd']:.3g} dates={len(rec['dates'])}")
# per-date CO2 effects at WDF (Fig. 3g-l method: eCO2/aCO2 - 1 of the mean over the decks sharing genotype/water on each date)
for site in ("WDF25", "WDF23"):
    if site not in sites: continue
    for q in [k for k in keys if k != "Tleaf_FLIR"]:
        bydate = {}
        for deck in sites[site]:
            for p in E.pairs(deck, q):
                bydate.setdefault(p[0], {}).setdefault("e" if "eCO2" in deck else "a", []).append((p[1], p[3]))
        eff = []
        for d in sorted(bydate):
            g = bydate[d]
            if "e" in g and "a" in g:
                oa, ma = mean([x[0] for x in g["a"]]), mean([x[1] for x in g["a"]]); oe, me = mean([x[0] for x in g["e"]]), mean([x[1] for x in g["e"]])
                eff.append((d, 100*(oe/oa-1), 100*(me/ma-1)))
        out[f"{site}/{q}/co2_effect_by_date"] = eff
        print(f"{site} {q} CO2 effect by date (obs%, sim%):", ["%s %+.1f/%+.1f" % e for e in eff])
        if eff: print(f"   range obs {min(e[1] for e in eff):+.1f}..{max(e[1] for e in eff):+.1f}  sim {min(e[2] for e in eff):+.1f}..{max(e[2] for e in eff):+.1f}")
json.dump(out, open(sys.argv[2], "w"), indent=1, default=str)
