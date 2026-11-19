#!/usr/bin/env python3
"""FACE-deck model runner for CroptimizR (port of Yaron's AgMIP croptimizr_run.py, 29 Oct 2026).

  --deck D --params @p.json --out sim.csv --work DIR [--log evals.csv]
Parameters (missing ones keep the deck value):
  Anthesis (deg Cd), GrainFill_TDU (= Maturity - Anthesis): stages 3_leaves..Heading follow Anthesis and Grain_filling
  keeps its share of the anthesis-to-maturity interval, deck-relative (calibrate_agmip.pheno_params);
  fruiting_efficiency (cultivar_yield.csv); Blade_L_max_cm (phenology_cultivar.txt).
Outputs (constant over the run except Date/LAI): anthesis_dayf, maturity_dayf = fractional model day the threshold is
crossed (calibrate_agmip.frac_day); yield_final (g/m2, last day); lai_peak (max LAI).
Decks are read from $FCAL_DECKS, the executable from $FCAL_EXE.
"""
import argparse, csv, json, os, shutil, subprocess, tempfile, time
from pathlib import Path

DECKS = Path(os.environ["FCAL_DECKS"]); EXE = os.environ["FCAL_EXE"]
SCALED = ("3_leaves", "Tillering", "Booting", "Heading")


def F(x):
    try:
        return float(x)
    except ValueError:  # Fortran drops the 'E' in tiny exponents
        return 0.0


def read_pheno(path):
    v = {}
    for l in path.read_text().splitlines():
        t = l.split()
        if len(t) >= 2 and not l.startswith("!"):
            v[t[1]] = float(t[0])
    return v


def set_rows(path, new):
    out = []
    for l in path.read_text().splitlines():
        t = l.split()
        if len(t) >= 2 and not l.startswith("!") and t[1] in new:
            l = f"{new[t[1]]:.3f}" + l[len(t[0]):]
        out.append(l)
    path.write_text("\n".join(out) + "\n")


def set_fe(path, fe):
    out = []
    for l in path.read_text().splitlines():
        if not l.startswith("!") and ",fruiting_efficiency" in l:
            l = f"{fe:.4f}" + l[l.index(","):]
        out.append(l)
    path.write_text("\n".join(out) + "\n")


def pheno_params(r0, A, G):
    A0, M0 = r0["Anthesis"], r0["Maturity"]
    p = {s: r0[s] * A / A0 for s in SCALED}
    p["Anthesis"] = A
    p["Grain_filling"] = A + (r0["Grain_filling"] - A0) / (M0 - A0) * G
    p["Maturity"] = A + G
    return p


def rd(p):
    with open(p) as f:
        r = csv.reader(f); h = [x.strip() for x in next(r)]
        return [dict(zip(h, [x.strip() for x in row])) for row in r]


def frac_day(ph, stage, thr):
    days = [int(r["day"]) for r in ph if r["stage_name"] == stage]
    c = {int(r["day"]): F(r["cum_tdu"]) for r in ph}
    if not days:  # not reached: extrapolate at the run's mean thermal-time rate
        last = max(c); return last + max(thr - c[last], 0) / max(c[last] / len(c), 1.0)
    d = min(days); c0, c1 = c.get(d - 1), c.get(d)
    if c0 is None or c1 is None or c1 <= c0:
        return float(d)
    return d - 1 + min(max((thr - c0) / (c1 - c0), 0.0), 1.0)


def run(deck, params, out_csv, work, log):
    src = DECKS / deck / "input"; r0 = read_pheno(src / "phenology_cultivar.txt")
    A = params.get("Anthesis", r0["Anthesis"]); G = params.get("GrainFill_TDU", r0["Maturity"] - r0["Anthesis"])
    Path(work).mkdir(parents=True, exist_ok=True)
    rdir = Path(tempfile.mkdtemp(prefix=f"{deck}_", dir=work)); t0 = time.time()
    try:
        shutil.copytree(src, rdir / "input"); (rdir / "output").mkdir()
        new = pheno_params(r0, A, G)
        if "Blade_L_max_cm" in params:
            new["Blade_L_max_cm"] = params["Blade_L_max_cm"]
        set_rows(rdir / "input" / "phenology_cultivar.txt", new)
        if "fruiting_efficiency" in params:
            set_fe(rdir / "input" / "cultivar_yield.csv", params["fruiting_efficiency"])
        subprocess.run([EXE], cwd=rdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1800)
        ph = rd(rdir / "output" / "daily_phenology.csv"); g = rd(rdir / "output" / "daily_growth.csv")
        kn = F(rd(rdir / "output" / "daily_yield_components.csv")[-1]["kernel_number_m2"])
    finally:
        shutil.rmtree(rdir, ignore_errors=True)
    a, m = frac_day(ph, "Anthesis", A), frac_day(ph, "Maturity", A + G)
    y = F(g[-1]["GrainYield_g_m2"]); lai = max(F(r["LAI_m2_m2"]) for r in g)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["Date", "LAI", "anthesis_dayf", "maturity_dayf", "yield_final", "lai_peak", "kernel_number"])
        for r in ph:
            w.writerow([f"{int(r['year']):04d}-{int(r['month']):02d}-{int(r['dom']):02d}",
                        next((x["LAI_m2_m2"] for x in g if x["day"] == r["day"]), "0"), a, m, y, lai, kn])
    if log:
        new_file = not Path(log).exists()
        with open(log, "a", newline="") as f:
            w = csv.writer(f)
            if new_file:
                w.writerow(["time", "deck", "Anthesis", "GrainFill_TDU", "fruiting_efficiency", "Blade_L_max_cm",
                            "anthesis_dayf", "maturity_dayf", "yield_final", "lai_peak", "seconds"])
            w.writerow([time.strftime("%H:%M:%S"), deck, A, G, params.get("fruiting_efficiency", ""),
                        params.get("Blade_L_max_cm", ""), round(a, 3), round(m, 3), round(y, 2), round(lai, 3),
                        round(time.time() - t0)])


def target_days(deck, ref_run, dates):
    """Model day of each observed date minus 0.5 (the whole simulated day equals the observed one)."""
    ph = rd(Path(ref_run) / deck / "output" / "daily_phenology.csv")
    dayof = {f"{int(r['year']):04d}-{int(r['month']):02d}-{int(r['dom']):02d}": int(r["day"]) for r in ph}
    return {k: dayof[v] - 0.5 for k, v in dates.items() if v}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--deck"); ap.add_argument("--params", default="{}"); ap.add_argument("--out")
    ap.add_argument("--work", default=tempfile.gettempdir()); ap.add_argument("--log")
    ap.add_argument("--targets", nargs=3, metavar=("TARGETS_JSON", "REF_RUN", "OUT_CSV"))
    a = ap.parse_args()
    if a.targets:
        T = json.load(open(a.targets[0])); rows = []
        for d, dates in T["phen"].items():
            r0 = read_pheno(DECKS / d / "input" / "phenology_cultivar.txt"); td = target_days(d, a.targets[1], dates)
            ph = rd(Path(a.targets[1]) / d / "output" / "daily_phenology.csv")
            rows.append(dict(deck=d, kind="phen", start_date=f"{int(ph[0]['year']):04d}-{int(ph[0]['month']):02d}-{int(ph[0]['dom']):02d}",
                             anthesis_dayf=td["Anthesis"], maturity_dayf=td["Maturity"], yield_final="", lai_peak="",
                             A0=r0["Anthesis"], G0=r0["Maturity"] - r0["Anthesis"], FE0="", L0=r0.get("Blade_L_max_cm")))
        for d, v in T["cult"].items():
            r0 = read_pheno(DECKS / d / "input" / "phenology_cultivar.txt")
            fe = [l.split(",")[0] for l in (DECKS / d / "input" / "cultivar_yield.csv").read_text().splitlines()
                  if not l.startswith("!") and ",fruiting_efficiency" in l][0]
            ph = rd(Path(a.targets[1]) / d / "output" / "daily_phenology.csv")
            rows.append(dict(deck=d, kind="cult", start_date=f"{int(ph[0]['year']):04d}-{int(ph[0]['month']):02d}-{int(ph[0]['dom']):02d}",
                             anthesis_dayf="", maturity_dayf="", yield_final=v["grain_yield"], lai_peak=v["LAI_peak"],
                             A0=r0["Anthesis"], G0=r0["Maturity"] - r0["Anthesis"], FE0=float(fe), L0=r0.get("Blade_L_max_cm")))
        with open(a.targets[2], "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print(open(a.targets[2]).read())
    else:
        p = json.loads(Path(a.params[1:]).read_text() if a.params.startswith("@") else a.params)
        run(a.deck, p, a.out, a.work, a.log)
