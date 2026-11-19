#!/usr/bin/env python3
# maricopa_soil_ll.py <archive .ods> <decks dir>
# Rewrites input/SoilLayer_vgm.csv of the MARICOPA_* decks in place (water audit section 7).
# Per layer, lower limit = min(archive SLLL, driest rep-averaged measured water content of any plot
# sharing the same rings in the same season, i.e. same Exp_year and CO2 treatment ME).
# Where the lower limit falls, the van Genuchten curve is rebuilt the way the v1.0 curve was built from
# the archive table: theta_s = SLSAT, theta_r = LL/2, curve through SLDUL at -33 kPa and LL at -1.5 MPa.
# Layers whose lower limit does not fall keep their v1.0 parameters.
import sys, csv, os, zipfile
import xml.etree.ElementTree as ET

T = '{urn:oasis:names:tc:opendocument:xmlns:table:1.0}'
O = '{urn:oasis:names:tc:opendocument:xmlns:office:1.0}'


def load(path):
    root = ET.fromstring(zipfile.ZipFile(path).read('content.xml'))
    out = {}
    for t in root.iter(T + 'table'):
        R = []
        for r in t.iter(T + 'table-row'):
            vals = []
            for c in r:
                n = int(c.get(T + 'number-columns-repeated', '1'))
                v = c.get(O + 'value') or ''.join(c.itertext())
                vals += [v] * min(n, 40)
            while vals and vals[-1] == '':
                vals.pop()
            if vals:
                R.append(vals)
        out[t.get(T + 'name')] = R
    return out


def theta(h, tr, ts, a, n):
    m = 1 - 1 / n
    return tr + (ts - tr) * (1 + (a * h) ** n) ** (-m)


def alpha_for(n, LL, tr, ts, hLL):
    m = 1 - 1 / n
    Se = (LL - tr) / (ts - tr)
    return (Se ** (-1 / m) - 1) ** (1 / n) / hLL


def solve(LL, DUL, SAT, hD, hL):
    tr = LL / 2
    f = lambda n: theta(hD, tr, SAT, alpha_for(n, LL, tr, SAT, hL), n) - DUL
    lo, hi = 1.05, 4.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    n = (lo + hi) / 2
    return tr, SAT, alpha_for(n, LL, tr, SAT, hL), n


sh = load(sys.argv[1]); decks = sys.argv[2]
arch = {int(r[1]): (float(r[2]), float(r[3]), float(r[4])) for r in sh['Soil_layers'][1:] if r[0] == 'AZMC920001'}
meta = {r[0]: (r[16], r[4]) for r in sh['Metadata'][1:]}       # Exp_year, ME (CO2)
bot = [5, 15, 30, 50, 60, 70, 90, 130, 170, 210]                # deck layer bottoms, cm
# archive depth classes: D30 = 0-30 cm, D50 = 30-50 cm, then 20 cm steps
probe = {0: ['D30'], 1: ['D30'], 2: ['D30'], 3: ['D50'], 4: ['D70'], 5: ['D70'], 6: ['D90'],
         7: ['D110', 'D130'], 8: ['D150', 'D170'], 9: ['D190', 'D210']}
A = sh['Soil_moisture_Avg_over_Reps']; h = A[0]


def lay_series(trno):
    out = {k: [] for k in probe}
    for r in A[1:]:
        if r[0] != trno:
            continue
        for k, cols in probe.items():
            v = [r[h.index(c)] for c in cols]
            if all(x not in ('', 'NA') for x in v):
                out[k].append(sum(map(float, v)) / len(v))
    return out


mins = {t: {k: min(s) for k, s in lay_series(t).items() if s} for t in meta}
hD, hL = 33 / 9.80665, 1500 / 9.80665
for d in sorted(os.listdir(decks)):
    if not d.startswith('MARICOPA_'):
        continue
    t = d.split('_')[1]
    grp = [g for g in mins if meta[g] == meta[t]]
    f = os.path.join(decks, d, 'input', 'SoilLayer_vgm.csv')
    rows = list(csv.DictReader(open(f)))
    for k, b in enumerate(bot):
        LL, DUL, SAT = arch[b]
        meas = min(mins[g][k] for g in grp if k in mins[g])
        if meas < LL - 1e-9:
            tr, ts, a, n = solve(meas, DUL, SAT, hD, hL)
            rows[k].update({'theta_r': '%.4f' % tr, 'alpha_m-1': '%.4f' % a, 'n': '%.4f' % n})
            print('%s L%d LL %.3f -> %.3f (group %s)' % (t, k + 1, LL, meas, '+'.join(grp)))
    with open(f, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator='\n')
        w.writeheader(); w.writerows(rows)
