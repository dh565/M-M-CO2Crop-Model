"""Part 2 Table 4 rows that no number file covers (stage dates, WDF23 straw DM, chaff N, Braunschweig DM at anthesis and
stem WSC, Maricopa harvest index and Maricopa energy-balance rows), computed with the paper scripts' own loaders.
usage: python t4extra.py <paper_scripts dir> <out json>   (REF_OUTPUTS, OBS_DIR, OUT_DIR from env, as config.py)

Definitions (each one reproduces the submitted v1.0 cells from the published Z2 outputs + published obs; see README.md):
- Stage dates: dates of emergence, anthesis, maturity (Maricopa: anthesis, maturity) from face_results.phenology_rows
  (calendar-date pairing); n = number of paired dates; 1sigma = share of dates with |diff| <= 2 d.
  |dDelta| reported two ways: abs_mean_listed (mean over the n listed dates = caption text) and abs_mean_scored
  (mean over all SCORED_STAGES = part2_numbers.json phenology_abs_dev_days 'mean').
- Treatment rows: per treatment (obs, obs_sd, model); n, mean obs, mean model, bias, RMSD over treatments;
  1sigma = share with |model - obs| <= obs sd (15 % of |obs| where no sd) - make_part2_tables.py rule.
  straw (WDF23)         obs registry 'straw' (harvest Straw biomass), model W_stem + W_leaf_tot (green + dead) at the last day
  chaff N % DM          obs ORGAN_N_OBS chaff ('Chaff' / Braunschweig 'N% CH'), model chaffN_conc_gN_gDM x 100 at the last day
  chaff N g N m-2       obs = chaff N % x observed chaff DM / 100 (sd in quadrature), model N_chaff_gN_m2 at the last day
  AGB at anthesis (BRA) obs 'AGB_anthesis' (SUMDW sample before the harvest), model AGB on that date (P2.AGB_ANTH)
  stem WSC (BRA)        obs 'stem_WSC_mg_g', model W_reserve_end / W_stem x 1000 (paperfigs.mod)
  HI (Maricopa)         obs = HI column of observations/obs_registry_braunschweig_maricopa.csv (grain / AGB of the measured
                        table, i.e. the registry value BEFORE paperfigs.init synchronises it with obs()=archive HIAM);
                        model grain / AGB; 15 % rule.  The HIAM variant (what the agreement csv uses) is also given.
- Maricopa energy balance: pooled plot-days of paperfigs2.MAR_SEASONS treatments (906, 908, 910, 912, 914, 916, 913, 915),
  maricopa_eb.daily_means: day = hour-ending 07-18, night = complement, 24-h = hours 1-24, >= 6 valid hours, Bowen-ratio
  plausibility screen; n = plot-days, Obs/Sim = mean over plot-days, Bias, RMSD over plot-days; no 1sigma.
"""
import sys, os, json, math, csv
SCR = sys.argv[1]; sys.path.insert(0, SCR); sys.path.insert(0, os.path.join(SCR, 'lib'))
import config
import paperfigs as PF
PF.init(config.REF_OUTPUTS)
import paperfigs2 as PF2
import maricopa_eb as ME
P2 = PF.P2; NAN = float('nan'); isnan = PF.isnan; fl = P2.fl
site = lambda s: [d for d in PF.FACE if PF.DECKS[d]['site'] == s]
out = {}


def trow(key, triples):
    t = [(d, o, s, m) for d, o, s, m in triples if not isnan(o) and not isnan(m)]
    n = len(t)
    if not n: out[key] = dict(n=0); return
    O = [x[1] for x in t]; M = [x[3] for x in t]
    tol = lambda o, s: s if (not isnan(s) and s > 0) else 0.15 * abs(o)
    out[key] = dict(n=n, obs=sum(O) / n, sim=sum(M) / n, bias=sum(m - o for o, m in zip(O, M)) / n,
                    rmsd=math.sqrt(sum((m - o) ** 2 for o, m in zip(O, M)) / n),
                    frac1sd=sum(abs(m - o) <= tol(o, s) for _, o, s, m in t) / n,
                    per_deck={d: [o, s, m] for d, o, s, m in t})


# ---- stage dates
for S_, stages in (('WDF23', ('Emergence', 'Anthesis', 'Maturity')), ('WDF25', ('Emergence', 'Anthesis', 'Maturity')),
                   ('BRAUNSCHWEIG', ('Emergence', 'Anthesis', 'Maturity')), ('MARICOPA', ('Anthesis', 'Maturity'))):
    lst, scored, per = [], [], {}
    for d in site(S_):
        for st, md, od, diff, sc, src in P2.PHEN[d]:
            if od and not isnan(diff):
                if sc: scored.append(abs(diff))
                if st in stages: lst.append(diff); per.setdefault(d, {})[st] = [md, od, diff]
    a = [abs(x) for x in lst]
    out[S_ + '/stage_dates'] = dict(n=len(lst), abs_mean_listed=sum(a) / len(a), abs_mean_scored=sum(scored) / len(scored),
                                    n_scored=len(scored), bias=sum(lst) / len(lst), frac1sd=sum(x <= 2 for x in a) / len(a),
                                    n_within2=sum(x <= 2 for x in a), per_deck=per)

# ---- WDF23 straw
trow('WDF23/straw', [(d, *PF.obs(d, 'straw'), P2.MODELS[d].v['W_stem'] + P2.MODELS[d].v['W_leaf_tot']) for d in site('WDF23')])

# ---- chaff N
for S_, org in (('WDF23', 'Chaff'), ('WDF25', 'Chaff'), ('BRAUNSCHWEIG', 'N% CH')):
    pct, gm = [], []
    for d in site(S_):
        nl = P2.MODELS[d].n[-1]
        o, s = P2.ORGAN_N_OBS.get((d, org), (NAN, NAN))
        pct.append((d, o, s, fl(nl['chaffN_conc_gN_gDM']) * 100.0))
        oc, sc = PF.obs(d, 'chaff')
        og = o * oc / 100.0 if not isnan(o) and not isnan(oc) else NAN
        sg = og * math.sqrt((s / o) ** 2 + (sc / oc) ** 2) if not isnan(og) and not isnan(s) and not isnan(sc) else NAN
        gm.append((d, og, sg, fl(nl['N_chaff_gN_m2'])))
    trow(S_ + '/chaffN_pct', pct); trow(S_ + '/chaffN_g', gm)

# ---- Braunschweig DM at anthesis, stem WSC
trow('BRAUNSCHWEIG/AGB_anthesis', [(d, *PF.obs(d, 'AGB_anthesis'), P2.AGB_ANTH[d]) for d in site('BRAUNSCHWEIG')])
trow('BRAUNSCHWEIG/stem_WSC', [(d, *PF.obs(d, 'stem_WSC_mg_g'), PF.mod(d, 'stem_WSC_mg_g')) for d in site('BRAUNSCHWEIG')])

# ---- Maricopa HI: registry csv value (submitted) and archive HIAM via obs() (agreement csv)
reg = {}
for r in csv.DictReader(l for l in open(config.OBS_REGISTRY_OPEN) if not l.startswith('#')):
    if r['quantity'] == 'HI' and r['site'] == 'MARICOPA':
        reg[r['deck']] = (fl(r['obs']), fl(r['obs_sd']) if r['obs_sd'] else NAN)
trow('MARICOPA/HI_registry', [(d, *reg.get(d, (NAN, NAN)), PF.mod(d, 'HI')) for d in site('MARICOPA')])
trow('MARICOPA/HI_archive_HIAM', [(d, *PF.obs(d, 'HI'), PF.mod(d, 'HI')) for d in site('MARICOPA')])

# ---- Maricopa energy balance, pooled plot-days
marv = {"Rn_day": ("rn", "Rn_total_W_m2", False), "Rn_night": ("rn", "Rn_total_W_m2", True), "LE_day": ("le", "LE", False),
        "H_day": ("h", "H_total_W_m2", False), "G_24h": ("g", "G_W_m2", None)}
for q, (vo, vm, night) in marv.items():
    X, Y = [], []
    for lab_, ta, te in PF2.MAR_SEASONS:
        for tr in (ta, te):
            dd = ME.daily_means(P2.PKG, tr, vo, vm, hours=(1, 24)) if night is None else ME.daily_means(P2.PKG, tr, vo, vm, night=night)
            X += [v[1] for v in dd]; Y += [v[3] for v in dd]
    n = len(X)
    out['MARICOPA/' + q] = dict(n=n, obs=sum(X) / n, sim=sum(Y) / n, bias=sum(y - x for x, y in zip(X, Y)) / n,
                                rmsd=math.sqrt(sum((y - x) ** 2 for x, y in zip(X, Y)) / n))

for k, v in out.items():
    print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a != 'per_deck'})
json.dump(out, open(sys.argv[2], 'w'), indent=1, default=str)
